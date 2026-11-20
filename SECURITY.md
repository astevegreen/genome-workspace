# Security & Privacy Policy

## Privacy Warning for Contributors & Users

This repository contains **reproducible bioinformatics scripts only**. Never commit, upload, or include in issues/pull requests:
- Raw or processed genomic files (`*.vcf*`, `*.gvcf*`, `*.bam*`, `*.cram*`, `*.fastq*`, `*.bed*`)
- Personal clinical, carrier screening, pharmacogenomic, or polygenic risk reports (`results/`)
- Local sample registries or path configurations (`config/samples.tsv`, `config/paths.env`, `.env`)

All of these paths and extensions are excluded by `.gitignore`, and raw source directories are treated as strictly read-only.

## Reporting a Vulnerability

If you discover a security or privacy vulnerability in this codebase, please report it privately via **GitHub Private Vulnerability Reporting** (`Security` -> `Report a vulnerability` on this repository) rather than opening a public issue.

## Clinical Disclaimer

This software is intended strictly for educational and research-grade local exploration. It is **not** a medical device or a substitute for CLIA-certified clinical genetic testing or genetic counseling.
