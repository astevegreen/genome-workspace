# Storage and Backup: BAM to CRAM, and a Second Copy

Follow-up item 6. The owner runs every step here. Nothing in the pipeline deletes or moves a BAM, and nothing changes
the `data/alignments/<label>/` symlinks.

## What the conversion does

`scripts/bam_to_cram.sh` reads a BAM and writes a verified CRAM next to the BAM's real file. The symlink is followed,
so the CRAM lands in the same folder on the external drive. Alongside it go `<name>.cram.crai`, `<name>.cram.md5` and
`<name>.cram.verify.txt`. The CRAM is written as `<name>.part.cram` and renamed only after every check passes:

| Check | What it proves |
|---|---|
| `samtools quickcheck` | the CRAM is complete (EOF block present) |
| `flagstat` BAM = CRAM | same read totals by flag class |
| `idxstats` BAM = CRAM | same mapped / unmapped reads on every contig |
| `samtools checksum -a` BAM = CRAM | every record identical, in order: flags, sequence, qualities, positions, MAPQ/CIGAR, mates, tags |
| region read counts BAM = CRAM | the CRAM index works (UGT1A1, SMN, HLA, CYP2D6, FMR1, chrY, chrM, a decoy contig) |

### Why no reference (`no_ref`), and why CRAM 3.0

The BAMs were aligned to GRCh38 "no-alt + hs38d1 decoy" (2,581 contigs: 2,385 decoys, EBV, and UCSC-named unplaced
contigs). `data/reference/GRCh38.primary_assembly.fa.bgz` has only 195 contigs, and about 1.8% of the reads lie on
contigs it lacks.

Converting against that FASTA makes samtools 1.22 silently switch the whole file to `embed_ref=2`. The resulting
`.crai` index gives every slice an alignment span of 0. Quickcheck, flagstat, idxstats and even the record checksums
all pass, but region queries miss reads: in a test, 9,235 of 19,234 reads in a 100 kb window, and `bcftools mpileup`
called no sites at all. The script therefore refuses a `--reference` that lacks any of the BAM's contigs.

The default is a **no_ref CRAM**:
- it stores the bases themselves, so it decodes without any reference FASTA;
- its index is correct;
- it is about 66% of the BAM's size. A reference-based CRAM against the exact alignment reference would be about 59%,
  but that FASTA would then have to be kept as carefully as the CRAMs.

Estimated sizes are 117 GB → about 77 GB and 89 GB → about 59 GB. The CRAM is written as **version 3.0**, because the
ExpansionHunter v5 binary cannot decode samtools' default 3.1 codecs.

## Step 3 and Step 10 on a CRAM

These were checked on a 3-million-read slice (HLA, UGT1A1, CYP2D6, SMN, FMR1, HTT, plus the ~3,000 depth-normalisation
regions Cyrius and SMNCopyNumberCaller use), comparing the BAM with a CRAM made by this script:

| Tool | CRAM 3.0 (this script) | Notes |
|---|---|---|
| Step 3 identity genotyping (`bcftools mpileup`) | identical | |
| Step 3 genetic sex (index read counts) | identical | `get_index_statistics()` returns zeros for any CRAM, so Step 3 now uses `samtools idxstats` for CRAMs |
| UGT1A1 TA-repeat counter | identical | |
| ExpansionHunter | identical | fails on CRAM 3.1; Step 10 now reports `UNSUPPORTED_CRAM` instead of crashing |
| SMNCopyNumberCaller | identical | Step 10 passes `--reference` for CRAMs (needed only for reference-based ones) |
| Cyrius | identical | same |
| T1K | runs via a temporary BAM; HLA-A, HLA-B, HLA-DQA1 and HLA-DQB1 identical | see below |
| SV verification, phasing checks (pysam) | open the CRAM with the reference when given | |

T1K's `bam-extractor` is built on the samtools 0.1.19 API: it cannot read CRAM, and it opens its input twice, so a
pipe does not work either. For a CRAM, Step 10 writes `t1k/t1k_input.bam` holding the reads T1K's whole-file scan
would pick: the HLA gene regions, every decoy / unplaced / EBV contig, and the unmapped reads.

On the full BAM this reproduced 96% of T1K's candidate reads. HLA-A, HLA-B, HLA-DQA1 and HLA-DQB1 (the genes the card,
PharmCAT and the HLA associations use) came out identical, qualities included. Two genes the pipeline does not use
differed slightly: an HLA-DRB1 quality, and a low-quality (q = 3) second HLA-DRB3 allele.

Step 10 re-runs every caller when the alignment file changes. So once the pipeline is pointed at a CRAM (optional,
below), compare the new `bam_blind_spots.json` with a copy kept from the BAM run.

## Run checklist (owner)

Run in WSL from the repo root, one sample at a time, in a terminal that stays open: each sample reads the BAM about
three times and takes a few hours on a USB drive.

1. [ ] D: is mounted and both links resolve: `ls -lL data/alignments/*/`.
2. [ ] `samtools --version` shows 1.20 or newer (`checksum` is needed; 1.22 is installed).
3. [ ] Convert and verify each sample. Use the link path; the CRAM goes next to the real BAM on D: (about 140 GB in
       total; the script checks for free space):

   ```bash
   cd /mnt/c/Antigravity/Genome
   bash scripts/bam_to_cram.sh --threads 8 data/alignments/<label>/<name>.bam
   ```

4. [ ] Each run ends with `VERIFIED -> ...cram`. Read `<name>.cram.verify.txt`: every line says identical or PASS,
       and the region counts match. On `VERIFICATION FAILED`, the `.part.cram` and the check outputs are kept for
       inspection. Nothing is renamed and the BAM is untouched.
5. [ ] Make the second copy off D: (an internal disk, another external drive, or encrypted cloud storage) of:
   - [ ] each `<name>.cram`, `.cram.crai`, `.cram.md5` and `.cram.verify.txt`;
   - [ ] the original VCFs from the raw source folders listed in `config/samples.tsv`;
   - [ ] `data/snapshots/` (the annual "what changed" snapshots: the only record of earlier years);
   - [ ] `config/samples.tsv` and `config/paths.env`.

   ```bash
   rsync -av --progress "/mnt/d/<folder>/<name>.cram"* "/mnt/<backup-drive>/genome-backup/<label>/"
   ```

6. [ ] Verify the copy on the destination, not the source:

   ```bash
   cd "/mnt/<backup-drive>/genome-backup/<label>/" && md5sum -c <name>.cram.md5
   ```

7. [ ] If the copy is in the cloud, it is encrypted before upload (e.g. `age`, `gpg`, or an encrypted archive):
       genomes identify the person and their relatives.
8. [ ] Write down where the second copy is and when it was last verified; repeat step 6 at each annual checkup.
9. [ ] Only after steps 3-8: decide whether to keep the BAMs. Deleting them is your call and is never done by a script.
       With the BAMs gone, the CRAM plus its second copy is the only alignment, and Step 10 needs the CRAM.

### Optional: point the pipeline at the CRAM

The pipeline uses the first `*.bam` in `data/alignments/<label>/`, then `*.cram`. To switch a sample, replace the BAM
links with CRAM links. This removes only the links, not the files.

```bash
cd /mnt/c/Antigravity/Genome
cp results/reports/<label>/bam_blind_spots.json /tmp/<label>.bam_blind_spots.from_bam.json
rm data/alignments/<label>/<name>.bam data/alignments/<label>/<name>.bam.bai          # symlinks only
ln -s "/mnt/d/<folder>/<name>.cram" "data/alignments/<label>/<name>.cram"
ln -s "/mnt/d/<folder>/<name>.cram.crai" "data/alignments/<label>/<name>.cram.crai"
.venv-wsl/bin/python scripts/03_sample_identity_qc.py --sample <label>
.venv-wsl/bin/python scripts/10_bam_blind_spots.py --sample <label>
```

Then compare the calls in `results/reports/<label>/bam_blind_spots.json` with the saved copy. Expect the same SMN,
CYP2D6, repeat and UGT1A1 results and the same HLA-A/B/DQ alleles.
