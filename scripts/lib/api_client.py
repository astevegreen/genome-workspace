"""Cached external API clients for gnomAD v4 and Ensembl REST (gene lookup & VEP)."""

from __future__ import annotations

import datetime
import json
import os
import time
from pathlib import Path
from typing import Any

import requests

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
ENSEMBL_API = "https://rest.ensembl.org"
ENSEMBL_VEP_URL = f"{ENSEMBL_API}/vep/human/region"

GNOMAD_VARIANT_QUERY = """
query VariantInfo($variantId: String!, $datasetId: DatasetId!) {
  variant(variantId: $variantId, dataset: $datasetId) {
    variant_id
    exome {
      ac
      an
      af
      homozygote_count
      populations {
        id
        ac
        an
      }
    }
    genome {
      ac
      an
      af
      homozygote_count
      populations {
        id
        ac
        an
      }
    }
  }
}
"""

GNOMAD_GENE_QUERY = """
query GeneConstraint($geneSymbol: String!, $referenceGenome: ReferenceGenomeId!) {
  gene(gene_symbol: $geneSymbol, reference_genome: $referenceGenome) {
    symbol
    gene_id
    name
    gnomad_constraint {
      pLI
      oe_lof_upper
      oe_mis_upper
    }
  }
}
"""


class CachedAnnotationClient:
    """Disk-backed JSON cache for gnomAD v4 & Ensembl API lookups with TTL & retries."""

    def __init__(
        self,
        cache_path: Path,
        sleep_sec: float = 0.35,
        max_age_days: int = 365,
        refresh: bool | None = None,
    ) -> None:
        self.cache_path = cache_path
        self.sleep_sec = sleep_sec
        self._dirty = False
        now_iso = datetime.datetime.now(datetime.timezone.utc).isoformat()
        self._cache: dict[str, Any] = {
            "_meta": {"created_at": now_iso, "updated_at": now_iso},
            "gnomad_variants": {},
            "gnomad_genes": {},
            "ensembl_vep": {},
        }
        force_refresh = (
            refresh
            if refresh is not None
            else os.environ.get("REFRESH_API_CACHE", "").strip() in ("1", "true", "yes")
        )
        if not force_refresh and cache_path.exists():
            try:
                loaded = json.loads(cache_path.read_text(encoding="utf-8"))
                meta = loaded.get("_meta") if isinstance(loaded.get("_meta"), dict) else {}
                created_str = meta.get("created_at")
                expired = False
                if created_str and max_age_days > 0:
                    try:
                        created_dt = datetime.datetime.fromisoformat(created_str)
                        if created_dt.tzinfo is None:
                            created_dt = created_dt.replace(tzinfo=datetime.timezone.utc)
                        age_days = (datetime.datetime.now(datetime.timezone.utc) - created_dt).days
                        if age_days >= max_age_days:
                            expired = True
                    except Exception:
                        pass
                if not expired:
                    if meta:
                        self._cache["_meta"] = meta
                    for k in ("gnomad_variants", "gnomad_genes", "ensembl_vep"):
                        if isinstance(loaded.get(k), dict):
                            self._cache[k] = loaded[k]
                else:
                    self._dirty = True
            except Exception:
                pass
        elif force_refresh:
            self._dirty = True

    def save(self) -> None:
        if not self._dirty:
            return
        self._cache["_meta"]["updated_at"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.cache_path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(self._cache, indent=2), encoding="utf-8")
        tmp.replace(self.cache_path)
        self._dirty = False

    def query_gnomad_variant(self, variant_id: str) -> dict[str, Any]:
        """Fetch gnomAD v4 AF, popmax, and homozygote count with explicit status."""
        if variant_id.startswith("MT-"):
            return {
                "gnomad_af": None,
                "gnomad_popmax_af": None,
                "gnomad_popmax_pop": "",
                "gnomad_hom": None,
                "gnomad_source": "",
                "gnomad_status": "not_applicable",
            }
        cached = self._cache["gnomad_variants"].get(variant_id)
        if cached and cached.get("gnomad_status") in ("ok", "absent"):
            return dict(cached)

        out: dict[str, Any] = {
            "gnomad_af": None,
            "gnomad_popmax_af": None,
            "gnomad_popmax_pop": "",
            "gnomad_hom": None,
            "gnomad_source": "",
            "gnomad_status": "api_error",
        }
        resp = None
        for attempt in range(4):
            try:
                resp = requests.post(
                    GNOMAD_API,
                    json={
                        "query": GNOMAD_VARIANT_QUERY,
                        "variables": {"variantId": variant_id, "datasetId": "gnomad_r4"},
                    },
                    timeout=20,
                )
                time.sleep(self.sleep_sec)
                if resp.status_code == 200:
                    break
                if resp.status_code in (429, 500, 502, 503, 504):
                    time.sleep(1.5 * (2 ** attempt))
                    continue
                return out
            except Exception:
                time.sleep(1.5 * (2 ** attempt))
        if resp is None or resp.status_code != 200:
            return out
        try:
            payload = resp.json()
            if "errors" in payload and not payload.get("data", {}).get("variant"):
                # gnomAD returns an error message "Variant not found" when truly absent
                errs = " ".join(str(e.get("message", "")) for e in payload["errors"]).lower()
                if "not found" in errs:
                    out["gnomad_status"] = "absent"
                    self._cache["gnomad_variants"][variant_id] = dict(out)
                    self._dirty = True
                return out

            data = (payload.get("data") or {}).get("variant")
            if not data:
                out["gnomad_status"] = "absent"
                self._cache["gnomad_variants"][variant_id] = dict(out)
                self._dirty = True
                return out

            ex = data.get("exome")
            ge = data.get("genome")
            ac_tot = (ex["ac"] if ex else 0) + (ge["ac"] if ge else 0)
            an_tot = (ex["an"] if ex else 0) + (ge["an"] if ge else 0)
            hom_tot = (ex.get("homozygote_count", 0) if ex else 0) + (
                ge.get("homozygote_count", 0) if ge else 0
            )
            if an_tot > 0:
                out["gnomad_af"] = round(ac_tot / an_tot, 6)
                out["gnomad_hom"] = hom_tot
                out["gnomad_source"] = "+".join(
                    s for s, obj in (("exome", ex), ("genome", ge)) if obj
                )
                out["gnomad_status"] = "ok"
            else:
                out["gnomad_status"] = "absent"

            pop_ac: dict[str, int] = {}
            pop_an: dict[str, int] = {}
            for block in (ex, ge):
                if not block:
                    continue
                for p in block.get("populations") or []:
                    pid = p["id"]
                    if pid in ("XX", "XY") or "_" in pid:
                        continue
                    pop_ac[pid] = pop_ac.get(pid, 0) + (p["ac"] or 0)
                    pop_an[pid] = pop_an.get(pid, 0) + (p["an"] or 0)
            best_af = -1.0
            best_pop = ""
            for pid, an in pop_an.items():
                if an >= 100:
                    af = pop_ac[pid] / an
                    if af > best_af:
                        best_af = af
                        best_pop = pid
            if best_af >= 0:
                out["gnomad_popmax_af"] = round(best_af, 6)
                out["gnomad_popmax_pop"] = best_pop

            self._cache["gnomad_variants"][variant_id] = dict(out)
            self._dirty = True
        except Exception:
            pass
        return out

    def query_gnomad_gene(self, gene_symbol: str) -> dict[str, Any]:
        """Fetch gnomAD gene constraint (pLI/LOEUF) with retry/backoff & Ensembl REST fallback."""
        out: dict[str, Any] = {"gene_name": "", "ensembl_id": "", "pLI": None, "LOEUF": None}
        if not gene_symbol:
            return out
        cached = self._cache["gnomad_genes"].get(gene_symbol)
        if cached:
            return dict(cached)

        got_response = False
        for attempt in range(3):
            try:
                resp = requests.post(
                    GNOMAD_API,
                    json={
                        "query": GNOMAD_GENE_QUERY,
                        "variables": {"geneSymbol": gene_symbol, "referenceGenome": "GRCh38"},
                    },
                    timeout=15,
                )
                time.sleep(self.sleep_sec)
                if resp.status_code == 200:
                    got_response = True
                    g = (resp.json().get("data") or {}).get("gene")
                    if g:
                        out["gene_name"] = g.get("name") or ""
                        out["ensembl_id"] = g.get("gene_id") or ""
                        c = g.get("gnomad_constraint") or {}
                        if c.get("pLI") is not None:
                            out["pLI"] = round(float(c["pLI"]), 3)
                        if c.get("oe_lof_upper") is not None:
                            out["LOEUF"] = round(float(c["oe_lof_upper"]), 3)
                    break
                if resp.status_code in (429, 500, 502, 503, 504):
                    time.sleep(1.5 * (2 ** attempt))
                    continue
                break
            except Exception:
                time.sleep(1.5 * (2 ** attempt))

        if not out["gene_name"]:
            for attempt in range(3):
                try:
                    r = requests.get(
                        f"{ENSEMBL_API}/lookup/symbol/homo_sapiens/{gene_symbol}",
                        headers={"Accept": "application/json"},
                        timeout=10,
                    )
                    time.sleep(self.sleep_sec)
                    if r.status_code == 200:
                        got_response = True
                        j = r.json()
                        out["gene_name"] = (j.get("description") or "").split(" [Source:")[0]
                        out["ensembl_id"] = j.get("id") or out["ensembl_id"]
                        break
                    if r.status_code in (429, 500, 502, 503, 504):
                        time.sleep(1.5 * (2 ** attempt))
                        continue
                    break
                except Exception:
                    time.sleep(1.5 * (2 ** attempt))

        if got_response:
            self._cache["gnomad_genes"][gene_symbol] = dict(out)
            self._dirty = True
        return out

    def query_ensembl_vep_batch(
        self,
        variants: list[tuple[str, str, str, int, str, str]],
    ) -> dict[str, dict[str, Any]]:
        """
        Batch-annotate variants via Ensembl VEP POST /vep/human/region with retry/backoff.
        Input tuples: (variant_id, target_gene, chrom, pos, ref, alt).
        Returns {variant_id: {impact, consequence, hgvsc, hgvsp, vep_status}}.
        """
        results: dict[str, dict[str, Any]] = {}
        to_fetch: list[tuple[str, str, str]] = []  # (vid, target_gene, vep_line)

        for vid, gene, chrom, pos, ref, alt in variants:
            cache_key = f"{vid}|{gene}"
            cached = self._cache["ensembl_vep"].get(cache_key)
            if cached is not None and cached.get("vep_status", "ok") == "ok":
                entry = dict(cached)
                entry.setdefault("vep_status", "ok")
                results[vid] = entry
            else:
                to_fetch.append((vid, gene, f"{chrom} {pos} {vid} {ref} {alt} . . ."))

        gene_by_vid = {vid: gene for vid, gene, _ in to_fetch}
        for i in range(0, len(to_fetch), 100):
            batch = to_fetch[i : i + 100]
            vep_lines = [x[2] for x in batch]
            resp = None
            for attempt in range(4):
                try:
                    resp = requests.post(
                        ENSEMBL_VEP_URL,
                        headers={"Content-Type": "application/json", "Accept": "application/json"},
                        json={"variants": vep_lines, "canonical": 1, "hgvs": 1},
                        timeout=35,
                    )
                    time.sleep(0.25)
                    if resp.status_code == 200:
                        break
                    if resp.status_code in (429, 500, 502, 503, 504):
                        time.sleep(1.5 * (2 ** attempt))
                        continue
                    break
                except Exception:
                    time.sleep(1.5 * (2 ** attempt))

            if resp is None or resp.status_code != 200:
                for vid, _, _ in batch:
                    results[vid] = {
                        "impact": "",
                        "consequence": "",
                        "hgvsc": "",
                        "hgvsp": "",
                        "vep_status": "api_error",
                    }
                continue

            try:
                returned_vids: set[str] = set()
                for item in resp.json():
                    vid = item.get("id", "")
                    target_gene = gene_by_vid.get(vid, "")
                    best_tc = None
                    for tc in item.get("transcript_consequences") or []:
                        if tc.get("gene_symbol") == target_gene and tc.get("canonical") == 1:
                            best_tc = tc
                            break
                    if best_tc is None:
                        for tc in item.get("transcript_consequences") or []:
                            if tc.get("gene_symbol") == target_gene:
                                best_tc = tc
                                break
                    info = (
                        {
                            "impact": best_tc.get("impact", ""),
                            "consequence": ",".join(best_tc.get("consequence_terms") or []),
                            "hgvsc": best_tc.get("hgvsc", ""),
                            "hgvsp": best_tc.get("hgvsp", ""),
                            "vep_status": "ok",
                        }
                        if best_tc
                        else {
                            "impact": "",
                            "consequence": "",
                            "hgvsc": "",
                            "hgvsp": "",
                            "vep_status": "ok",
                        }
                    )
                    results[vid] = info
                    self._cache["ensembl_vep"][f"{vid}|{target_gene}"] = dict(info)
                    returned_vids.add(vid)
                for vid, gene, _ in batch:
                    if vid not in returned_vids:
                        empty = {
                            "impact": "",
                            "consequence": "",
                            "hgvsc": "",
                            "hgvsp": "",
                            "vep_status": "ok",
                        }
                        results[vid] = empty
                        self._cache["ensembl_vep"][f"{vid}|{gene}"] = empty
                self._dirty = True
            except Exception:
                for vid, _, _ in batch:
                    results[vid] = {
                        "impact": "",
                        "consequence": "",
                        "hgvsc": "",
                        "hgvsp": "",
                        "vep_status": "api_error",
                    }
        return results


def format_af(af: float | None, status: str | None = None) -> str:
    """Format allele frequency while never confusing API errors / unqueried with 'absent'."""
    if status == "api_error":
        return "API error"
    if status == "not_applicable":
        return "N/A (MT)"
    if af is None:
        if status == "absent":
            return "absent"
        return "not queried"
    if af == 0:
        return "0"
    if af < 0.0001:
        return f"{af:.2e}"
    return f"{af:.4%}"
