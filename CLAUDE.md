# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

Edi-Set-Flow is a Nextflow pipeline for detecting RNA editing sites (A-to-I via ADAR, C-to-T via APOBEC) in bulk RNA-seq and running differential analysis across experimental groups.

## Running the Pipeline

**Standard run (from a published revision):**
```bash
nextflow run bahlolab/Edi-Set-Flow \
    -revision 26.03-beta.2 \
    -profile hg38,singularity \
    -resume \
    --input sample_manifest.csv \
    --outdir esf_results
```

**Local development run:**
```bash
nextflow run main.nf \
    -profile hg38,singularity \
    -resume \
    --input misc/test_manifest.csv \
    --outdir esf_test_results
```

**Testing local `edisetr` changes** (forces reinstall from `edisetr/` source rather than using the bundled container version):
```bash
nextflow run main.nf -profile hg38,singularity -resume --input misc/test_manifest.csv --install_edisetr true
```

Always pass `-resume` to reuse cached Nextflow task outputs. Reference resources download to `esf_resources/` and are reused across runs.

## Architecture

### Top-level entry points

- **`main.nf`** — Defines all pipeline parameters with defaults, calls `WfEdiSetFlow.preflight_checks()`, then invokes the `ESF` workflow. All parameter documentation lives here.
- **`lib/WfEdiSetFlow.groovy`** — Groovy helper class: CSV/TSV parsing (`read_csv`, `read_tsv`) and preflight parameter validation.
- **`nextflow.config`** — Profile definitions for genomes (`hg38`, `mm10`, `mm39`), container engines (`singularity`, `apptainer`, `docker`), and per-process container image + resource labels.

### Workflow structure

`workflows/esf.nf` is the main workflow. It reads the sample manifest, resolves FASTQs (local paths or via `FASTERQDUMP` for SRA accessions), then calls five subworkflows in sequence:

```
SETUP → ALIGN → DISCO (COUNT pass 1) → COUNT (COUNT pass 2) → ANNOT → STATS
```

**`subworkflows/local/setup.nf`** — Downloads and caches all reference data: reference genome (bgzipped FASTA + index), UCSC repeat masker, GTF (Gencode), REDIportal catalog, dbSNP and gnomAD exclusion lists, and builds STAR or BWA-MEM2 indices. Also computes `target_regions` (union of GTF/REDIportal/custom BED minus exclusions) and splits them into `n_intervals` genomic intervals for parallelism.

**`subworkflows/local/align.nf`** — Optional fastp trimming → STAR or BWA-MEM2 alignment → samtools filtering (dedup, mapq) → mosdepth coverage → automatic strand inference via RSeQC `infer_experiment.py`. Emits per-sample tuples of `(sample, bam, bai, strand, coverage_bed, coverage_idx)`. Errors if stranded and unstranded samples are mixed. Samples supplying a pre-aligned `bam` in the manifest bypass fastp/alignment: they are name-collated (`COLLATE`, so `samtools fixmate -m` is valid) and mixed into the aligned-BAM channel between the aligner and SAMTOOLS, then follow the identical downstream path.

**`subworkflows/local/count.nf`** — Used **twice** with different options, aliased as `DISCO` (discovery) and `COUNT` (counting):
1. `WHERE`: intersects each sample's mosdepth-callable regions with target regions to get per-sample input BEDs for JACUSA2.
2. `JACUSA2`: runs JACUSA2 `call-1` per sample. JACUSA2's internal multithreading has a known stochastic bug, so the script manually splits the BED into per-CPU shards and runs parallel bash processes.
3. `TO_VCF`: converts JACUSA2 BED output to VCF.
4. `MERGE`: merges per-sample VCFs within each genomic interval using `bcftools merge`, computes allele counts/frequencies/PASS rates, and applies site-level filters.
5. `GATHER`: concatenates the interval-sharded VCFs into a single output.

The DISCO pass runs loose filters and emits a `sites_bed` (candidate sites); COUNT uses those sites as targets (`by_id=true`), runs stricter filters, and retains per-sample `NALT`/`NREF` genotype fields.

**`subworkflows/local/annot.nf`** — Runs Ensembl VEP per-interval (parallelised), merges, then runs `vcfanno` to add REDIportal membership and repeat masker context.

**`subworkflows/local/stats.nf`** — MultiQC aggregates QC files from fastp/STAR/mosdepth/RSeQC. `SCATTER` shards the annotated VCF. `REPORT` renders `bin/esf_report.qmd` via Quarto using the `edisetr` R package to fit per-site GLMs and produce the HTML report + CSV outputs.

### `edisetr` R package (`edisetr/`)

Companion R package bundled in the Docker image (`ghcr.io/bahlolab/edi-set-flow:<tag>`). Key functions:

- `read_edisites()` — reads sharded VCF/CSV data into a tidy data frame
- `fit_edisites()` — fits per-site GLMs using `group` + optional fixed effects; supports `quasibinomial`, `binomial`, `linear`, `arcsine` families; returns summary, contrasts, margins, and ANOVA tables with BH FDR correction

The Quarto report (`bin/esf_report.qmd`) calls these functions and generates the interactive HTML output. All `--report_*` params are forwarded as Quarto execute-params.

### `bin/` helper scripts

- `gtf_to_bed.awk` — converts GTF to BED for target region computation
- `jacusa_to_vcf.awk` — converts JACUSA2 BED output to VCF format
- `infer_strand.py` — parses RSeQC output to emit `FIRSTstrand`, `SECONDSTRAND`, or `UNSTRANDED`
- `splinter.awk` — splits target regions BED into N balanced genomic intervals

### Docker image

`docker/Dockerfile` builds from `rocker/verse:4.4.3`, installs CRAN dependencies, then installs `edisetr` from `edisetr/`. Published as `ghcr.io/bahlolab/edi-set-flow:<tag>` via `.github/workflows/docker-build.yml`. The `edisetr` label in `nextflow.config` points to this image; all other tools use BioContainers images.

## Key Configuration

Process resources and container images are set via labels in `nextflow.config`. To override resources for a specific tool:
```groovy
process {
    withLabel: 'jacusa2' { cpus = 16; memory = '80 GB' }
}
```