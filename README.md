# Edi-Set-Flow - Detecting RNA editing at scale
A robust pipeline for RNA editing detection and differential analysis in bulk RNA-seq

## Overview
<p align="center"><img src="img/Edi-Set-Flow.png"/></p>

## Getting Started

### 1) Requirements
- [Nextflow](https://www.nextflow.io/) (≥ 22.03.0)
- Container engine, any of:
    - [Apptainer](https://apptainer.org/)
    - [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/index.html)
    - [Docker](https://www.docker.com/)

### 2) Sample Manifest (CSV format)
- **Required Columns:**
    - "sample_id" - uniuqe identified for sample
    - "fastq1" - path to first fastq file
    - "fastq2" - path to second fastq file (if paired data)
    - "group" - experimental condition of interest
- **Optional Columns:**
    - Arbitrary covariates to be include in GLM as fixed effects
    - Note: Columns with numbers will be interpreted as numeric by the GLM, otherwise they will be treated as factors
- **Example:**
     ```
    sample_id,fastq1,fastq2,group,sex,age
    SRR1311086,/PATH/TO/SRR1311086_1.fastq.gz,/PATH/TO/SRR1311086_2.fastq.gz,cortex,male,50
    SRR1477080,/PATH/TO/SRR1477080_1.fastq.gz,/PATH/TO/SRR1477080_2.fastq.gz,cerebellum,female,60
    SRR1085825,/PATH/TO/SRR1085825_1.fastq.gz,/PATH/TO/SRR1085825_2.fastq.gz,hippocampus,male,50
    ...
    ```

### 3) Run Edi-Set-Flow
- **Example:**
    ```{bash}
    nextflow run bahlolab/Edi-Set-Flow \
        -r 25.07-beta.2 \
        -profile singularity \
        -resume \
        --input sample_manifest.csv \
        --report_fixed_effects sex,age \
        --outdir esf_results
    ```
- **Notes**:
    - Use `-profile apptainer` or `-profile docker` to use Apptainer or Docker
    - Currently only the human GRCh38 reference and Gencode v48 gene annotations are supported
    - Resources (e.g. reference genome) are downloaded and stored in directory 'esf_resources', use `--resource_dir /path/to/resource_dir` to override and share between runs
    
## Output Files

| File | Description |
|------|-------------|
| **`EdiSetFlow.report.html`** | Interactive Edi-Set-Flow HTML report &mdash; [see example (GTEx brain)](https://bahlolab.github.io/Edi-Set-Flow/) |
| **`EdiSetFlow.sample_counts.csv.gz`** | Reference and alternate allele counts per site and sample |
| **`EdiSetFlow.sample_summary.csv.gz`** | Summary statistics (median depth & editing rate) per sample |
| **`EdiSetFlow.site_summary.csv.gz`** | Per-site statistics and annotation (e.g. gene, region, REDIPORTAL, consequence) |
| **`EdiSetFlow.glm_summary.csv.gz`** | GLM coefficients, standard errors, and significance per site |
| **`EdiSetFlow.glm_anova.csv.gz`** | ANOVA test results for each model term per site |
| **`EdiSetFlow.glm_contrasts.csv.gz`** | Pairwise group comparisons for differential editing |
| **`EdiSetFlow.glm_margins.csv.gz`** | Estimated marginal editing rates per group |
| **`multiqc_report.html`** | MultiQC summary report (fastp, STAR, mosdepth, etc.) |


## References
> **Piechotta, M., Naarmann-de Vries, I. S., Wang, Q., Altmüller, J. & Dieterich, C. (2022)**  
> _RNA modification mapping with JACUSA2._  
> _Genome Biology_, **23**(1), 115.

> **D’Addabbo, P., Cohen-Fultheim, R., Twersky, I., Fonzino, A., Silvestris, D. A., Prakash, A., … & Picardi, E. (2025)**  
> _REDIportal: toward an integrated view of the A-to-I editing._  
> _Nucleic Acids Research_, **53**(D1), D233–D242.

> **McLaren, W., Gil, L., Hunt, S. E., Riat, H. S., Ritchie, G. R., Thormann, A., … & Cunningham, F. (2016)**  
> _The Ensembl Variant Effect Predictor._  
> _Genome Biology_, **17**, 1–14.

> **da Veiga Leprevost, F., Grüning, B. A., Alves Aflitos, S., Röst, H. L., Uszkoreit, J., Barsnes, H., … & Perez-Riverol, Y. (2017)**  
> _BioContainers: an open-source and community-driven framework for software standardization._  
> _Bioinformatics_, **33**(16), 2580–2582.

> **Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E. & Notredame, C. (2017)**  
> _Nextflow enables reproducible computational workflows._  
> _Nature Biotechnology_, **35**(4), 316–319.
