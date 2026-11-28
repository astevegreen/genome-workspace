# Genome workspace

The owner's personal 30x WGS (GRCh38) pipeline: carrier screening, ACMG actionable variants, pharmacogenomics,
HLA, and an annual "what changed" report. Python 3.12, run in WSL2 Ubuntu with `uv` and `.venv-wsl/`.
Runbook: `docs/pipeline.md`. Build items and owner decisions: `docs/follow-ups.md`.

## Privacy (hard rules; the repository is PUBLIC)
- Only code goes into git. `data/`, `results/`, `config/samples.tsv` and `config/paths.env` hold personal genomic
  data and stay gitignored; never read `config/paths.env`.
- Tests, commits, PR text, issue comments and handoff notes carry no genotypes, personal variants, findings, report
  contents, or paths that carry a personal name. Tests use synthetic data.
- Before every commit, read `git diff --cached` for personal data. Widen `.gitignore` as needed; never narrow it.

## Commands
- Tests, from Windows: `wsl -e bash -lc "cd /mnt/c/Antigravity/Genome && .venv-wsl/bin/python -m pytest tests/ -q"`.
  A worktree has no `.venv-wsl/`: cd into the worktree and run the main checkout's
  `/mnt/c/Antigravity/Genome/.venv-wsl/bin/python`. The suite takes about a minute.
- CI (`.github/workflows/ci.yml`) runs the same suite with `uv sync --frozen` on Ubuntu.
- Pipeline entry point: `bash scripts/run_annual_checkup.sh [--refresh-refs]` (every sample, then the couple);
  one sample: `bash scripts/run_sample.sh <label>`. Single steps: `docs/pipeline.md` section 2.
- Handoff tools: `node tools/handoff-notes.mjs` (Windows `node`); their tests: `node --test tools/` (a few minutes;
  not in CI).

## Before debugging a failed run as a code bug
- Raw VCFs live on the owner's Google Drive (paths in `config/samples.tsv`): check that they are reachable.
- BAMs and CRAMs in `data/alignments/<label>/` are symlinks to the external drive `D:`. With D: unplugged, Step 10
  skips; that is expected, not a bug.

## Git
- One branch per change from `origin/main`, `<type>/<slug>` (`feat` `fix` `chore` `docs` `test` `refactor`), in its
  own worktree; the main checkout may hold someone else's work. Conventional Commits, as in `git log`.
- Land every change as a pull request into `main`, squash-merged once CI is green.
- Force-pushes are blocked. After a squash-merge, update a branch by merging `main` into it, never by rebasing.

## What waits on the owner
- Anything medical or interpretive: what a variant, call or report means, and how a report words it. Ask, with
  your recommendation; never settle it yourself.
- Scope: a new analysis or step not already on `docs/follow-ups.md`.
- Tool installs: the owner runs `scripts/setup_bam_tools.sh` and `scripts/setup_pharmcat.sh` (they download
  binaries). Say which script to run and why, then wait.
- Money, identity and accounts; deleting data, reports, snapshots or history.
Decide technical questions yourself and give the reason in the PR.

## Finishing
End every task with:
- **What changed** - the behaviour, in plain English, one bullet per change the owner would notice
- **Tests** - what proves it (pytest count, CI on the PR); what could not be verified and why
- **Follow-ups** - what you added to `docs/follow-ups.md` (leave the line out when none)
- **Unsure about** - what you tried to settle and could not; or "nothing"
