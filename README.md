# Human hereditary cancer workflow

Variant calling in human hereditary cancer samples.



## Introduction

<!---This section of documentation typically contains a list of things the workflow can perform also any other intro.--->

This workflow is currently released to support a beta product and is not for general use. If you are interested in the workflow please contact your sales representative.

This repository contains a Nextflow workflow for analysing variation in human hereditary cancer samples. Specifically this workflow can perform the following:

* diploid variant calling
* structural variant calling
* analysis of modified base calls



## Compute requirements

Recommended requirements:

+ CPUs = 32
+ Memory = 128GB

Minimum requirements:

+ CPUs = 16
+ Memory = 32GB

Approximate run time: 40 minutes per sample

ARM processor support: False




## Install and run

<!---Nextflow text remains the same across workflows, update example cmd and demo data sections.--->
These are instructions to install and run the workflow on command line. You can also access the workflow via the [EPI2ME application](https://labs.epi2me.io/downloads/).

The workflow uses [Nextflow](https://www.nextflow.io/) to manage compute and software resources, therefore nextflow will need to be installed before attempting to run the workflow.

The workflow can currently be run using either [Docker](https://www.docker.com/products/docker-desktop) or
[Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker or singularity is installed. This is controlled by the [`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles) parameter as exemplified below.

It is not required to clone or download the git repository in order to run the workflow.
More information on running EPI2ME workflows can be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow. This will pull the repository in to the assets folder of nextflow and provide a list of all parameters available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-hereditary-cancer –-help
```
A demo dataset is provided for testing of the workflow. Please note, you will also need the Hereditary Cancer BED file to run this workflow. Please contact your sales to obtain this. The demo data can be downloaded using:
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-hereditary-cancer/wf-hereditary-cancer-demo.tar.gz
tar -xzvf wf-hereditary-cancer-demo.tar.gz
```
The workflow can be run with the demo data using:
```
nextflow run epi2me-labs/wf-hereditary-cancer \
--bam wf-hereditary-cancer-demo/demo.bam \
--bed <hereditary_cancer_bed>
-profile standard
```
For further information about running a workflow on the cmd line see https://labs.epi2me.io/wfquickstart/




## Related protocols

<!---Hyperlinks to any related protocols that are directly related to this workflow, check the community for any such protocols.--->

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.

Find related protocols in the [Nanopore community](https://community.nanoporetech.com/docs/).



## Input example

<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
The `--bam` input parameter for this workflow accepts a path to a single BAM file, or a folder containing multiple BAM files for the sample. A sample name can be supplied with `--sample`.

```
(i)                     (ii)    
input_reads.bam     ─── input_directory
                        ├── reads0.bam
                        └── reads1.bam
```



## Input parameters

### Workflow Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sv | boolean | Call for structural variants. | If this option is selected, structural variant calling will be carried out using Sniffles2. | True |
| snp | boolean | Call for small variants | If this option is selected, small variant calling will be carried out using Clair3. | True |
| mod | boolean | Enable output of modified calls to a bedMethyl file [requires input BAM with Ml and Mm tags] | This option is automatically selected and aggregation of modified calls with be carried out using modkit if Ml and Mm tags are found. Disable this option to prevent output of a bedMethyl file. | True |


### Main options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_name | string | Sample name to be displayed in workflow outputs. |  | SAMPLE |
| bam | string | BAM or unaligned BAM (uBAM) files for the sample to use in the analysis. | This accepts one of two cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files. A sample name can be supplied with `--sample`. |  |
| bed | string | A BED file enumerating regions to process for variant calling. | Please contact your sales representative to obtain the hereditary cancer BED file |  |
| basecaller_cfg | string | Name of the model to use for selecting a small variant calling model. | The workflow will attempt to find the basecaller model from the headers of your input data. If the model cannot be found in the header, it must be provided with this option as the basecaller model is required for small variant calling. The basecaller model is used to automatically select the appropriate small variant calling model. The model list shows all models that are compatible for small variant calling with this workflow. You should select 'custom' to override the basecaller_cfg with clair3_model_path. | dna_r10.4.1_e8.2_400bps_hac@v4.3.0 |
| out_dir | string | Directory for output of all workflow results. |  | output |
| store_dir | string | Where to store initial download of reference genome. | Reference genome will be downloaded as part of the workflow and saved in this location, on subsequent runs it will use this as the reference. |  |


### Multiprocessing Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| threads | integer | Set max number of threads to use for more intense processes (limited by config executor cpus) |  | 4 |
| ubam_map_threads | integer | Set max number of threads to use for aligning reads from uBAM (limited by config executor cpus) |  | 8 |
| ubam_sort_threads | integer | Set max number of threads to use for sorting and indexing aligned reads from uBAM (limited by config executor cpus) |  | 3 |
| ubam_bam2fq_threads | integer | Set max number of threads to use for uncompressing uBAM and generating FASTQ for alignment (limited by config executor cpus) |  | 1 |
| merge_threads | integer | Set max number of threads to use for merging alignment files (limited by config executor cpus) |  | 4 |
| modkit_threads | integer | Total number of threads to use in modkit modified base calling (limited by config executor cpus) |  | 4 |






## Outputs

Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Report of the alignment statistics | {{ alias }}.wf-human-alignment-report.html | Report summarising the results of the alignment statistics for the sample. | per-sample |
| JSON file of some base statistics | {{ alias }}.stats.json | This JSON file contains base statistics on the reads, mappings, SNPs and SVs for the sample. | per-sample |
| Report of the SNP workflow | {{ alias }}.wf-human-snp-report.html | Report summarising the results of the SNP subworkflow for the sample. | per-sample |
| Report of the SV workflow | {{ alias }}.wf-human-sv-report.html | Report summarising the results of the SV subworkflow for the sample. | per-sample |
| Short variant VCF | {{ alias }}.wf_snp.vcf.gz | VCF file with the SNPs for the sample. | per-sample |
| Structural variant VCF | {{ alias }}.wf_sv.vcf.gz | VCF file with the SVs for the sample. | per-sample |
| Modified bases BEDMethyl | {{ alias }}.wf_mods.bedmethyl.gz | BED file with the aggregated modification counts for the sample. | per-sample |
| Modified bases BEDMethyl (haplotype 1) | {{ alias }}.wf_mods.1.bedmethyl.gz | BED file with the aggregated modification counts for haplotype 1 of the sample. | per-sample |
| Modified bases BEDMethyl (haplotype 2) | {{ alias }}.wf_mods.2.bedmethyl.gz | BED file with the aggregated modification counts for haplotype 2 of the sample. | per-sample |
| Modified bases BEDMethyl (ungrouped) | {{ alias }}.wf_mods.ungrouped.bedmethyl.gz | BED file with the aggregated modification counts of non-haplotagged reads for the sample. | per-sample |
| Haplotagged alignment file | {{ alias }}.haplotagged.bam | BAM file with the haplotagged reads for the sample. | per-sample |
| Haplotagged alignment file index | {{ alias }}.haplotagged.bam.bai | The index of the resulting BAM file with the haplotagged reads for the sample. | per-sample |
| Mean coverage for each region | {{ alias }}.regions.bed.gz | The mean coverage in the individual regions of the genome in BED format. | per-sample |
| Coverage per region above the given thresholds | {{ alias }}.thresholds.bed.gz | The BED reporting the number of bases in each region that are covered at or above each threshold values (1x, 10x, 20x and 30x). | per-sample |
| Distribution of the proportion of total bases covered by a given coverage value | {{ alias }}.mosdepth.global.dist.txt | The cumulative distribution indicating the proportion of total bases covered by a given coverage value, both genome-wide and by sequence. | per-sample |
| Mean coverage per sequence and target region | {{ alias }}.mosdepth.summary.txt | The summary of mean depths per chromosome and within specified regions per chromosome. | per-sample |
| BEDgraph of the single-base coverage | {{ alias }}.per-base.bedgraph.gz | The single-base coverage of the genome in BED graph format. | per-sample |
| Gene level coverage summary | SAMPLE.gene_summary.tsv | A table where each gene of the input BED file has columns describing the percentage of positions along the gene region that are covered to a given threshold, and a column with the average coverage. | per-sample |




## Pipeline overview

<!---High level numbered list of main steps of the workflow and hyperlink to any tools used. If multiple workflows/different modes perhaps have subheadings and numbered steps. Use nested numbering or bullets where required.--->
The workflow is composed of 3 distinct subworkflows, each enabled by a command line option:

* [SNP calling](#3-small-variant-calling-with-clair3): `--snp`
* [SV calling](#4-structural-variant-sv-calling-with-sniffles2): `--sv`
* [Analysis of modified bases](#5-modified-base-calling-with-modkit): `--mod`

All of these options are enabled in this workflow by default.

### 1. Input and data preparation

The workflow relies on one primary input file:
1. Sequencing data for the sample in the form of a single [BAM file](https://samtools.github.io/hts-specs/SAMv1.pdf), either aligned or unaligned.

On first run, the workflow will download the following reference FASTA file: [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz). This is stored locally and will be used in subsequent runs.

The input BAM file can be generated using the [wf-basecalling](https://github.com/epi2me-labs/wf-basecalling/) workflow, which is up to date with the current dorado releases and models.

### 2. Data QC and pre-processing
The workflow starts by performing multiple checks of the input BAM file, as well as computing:
1. Depth of sequencing with [mosdepth](https://github.com/brentp/mosdepth);
2. Read alignment statistics with [fastcat](https://github.com/epi2me-labs/fastcat).

### 3. Small variant calling with Clair3

The workflow implements a deconstructed version of [Clair3](https://github.com/HKU-BAL/Clair3) (v1.0.4) to call germline variants. The appropriate model can be provided with the `--basecaller_cfg` option. The two models compatible with this workflow are `dna_r10.4.1_e8.2_400bps_hac@v4.3.0` and `dna_r10.4.1_e8.2_400bps_sup@v4.3.0`.
This workflow takes advantage of the parallel nature of Nextflow, providing optimal efficiency in high-performance, distributed systems. The workflow will automatically call small variants (SNPs and indels), collect statistics, annotate them with [SnpEff](https://pcingola.github.io/SnpEff/) (and additionally for SNPs, ClinVar details), and create a report summarising the findings.

The workflow will perform phasing of variants by using the `--phased` option. This will lead the workflow to use [whatshap](https://whatshap.readthedocs.io/) to perform phasing of the variants, with the option to use [longphase](https://github.com/twolinin/longphase) instead by setting `--use_longphase true`. The phasing will also generate a GFF file with the annotation of the phase blocks, enabling the visualisation of these blocks in genome browsers.

### 4. Structural variant (SV) calling with Sniffles2

The workflow allows for calling of SVs using long-read sequencing data with [Sniffles2](https://github.com/fritzsedlazeck/Sniffles).
The workflow will perform SV calling, filtering and generation of a report.
The SV workflow takes an optional `--tr_bed` option to specify tandem repeats in the reference sequence --- see the [sniffles](https://github.com/fritzsedlazeck/Sniffles) documentation for more information.
SVs will be phased using `--phased`.

### 5. Modified base calling with modkit

Modified base calling can be performed by specifying `--mod`. The workflow will call modified bases using [modkit](https://github.com/nanoporetech/modkit). 
The workflow will automatically check whether the files contain the appropriate `MM`/`ML` tags, required for running [modkit pileup](https://nanoporetech.github.io/modkit/intro_bedmethyl.html). If the tags are not found, the workflow will not run the individual analysis, but will still run the other subworkflows requested by the user.
The default behaviour of the workflow is to run modkit with the `--cpg --combine-strands` options set. It is possible to report strand-aware modifications by providing `--force_strand`, which will trigger modkit to run in default mode. The resulting bedMethyl will include modifications for each site on each strand separately.
The modkit run can be fully customized by providing `--modkit_args`. This will override any preset, and allow full control over the run of modkit.
Haplotype-resolved aggregated counts of modified bases can be obtained with the `--phased` option. This will generate three distinct BEDMethyl files with the naming pattern `{{ alias }}.wf_mods.{{ haplotype }}.bedmethyl.gz`, where `haplotype` can be `1`, `2` or `ungrouped`.

### 6. Phasing variants

Variant phasing is switched on simply using the `--phased` option.
By default, the workflow uses [whatshap](https://whatshap.readthedocs.io/) to perform phasing of the variants, with the option to use [longphase](https://github.com/twolinin/longphase) instead by setting `--use_longphase true`.
The workflow will automatically turn on the necessary phasing processes based on the selected subworkflows.
The behaviour of the phasing is summarised in the below table:

|         |        |         |            | Phased SNP VCF | Phased SV VCF | Joint SV+SNP phased VCF | Phased bedMethyl |
|---------|--------|---------|------------|----------------|---------------|-------------------------|------------------|
| `--snp` | `--sv` | `--mod` | `--phased` |                |               |          &check;        |       &check;    |
| `--snp` | `--sv` |         | `--phased` |                |               |          &check;        |                  |
| `--snp` |        |         | `--phased` |     &check;    |               |                         |                  |
|         | `--sv` |         | `--phased` |                |     &check;   |                         |                  |
|         |        | `--mod` | `--phased` |                |               |                         |       &check;    |

The joint physical phasing of SNP and SVs can only be performed with [longphase](https://github.com/twolinin/longphase) by selecting the options: `--phased --snp --sv`. Setting `--use_longphase false` will not disable the final joint phasing with longphase.

In some circumstances, users may wish to keep the separate VCF files before joint phasing. This can be done with `--output_separate_phased`.

### 7. Variant annotation
Annotation will be performed automatically by the SNP and SV subworkflows, and can be disabled by the user with `--annotation false`. The workflow will annotate the variants using [SnpEff](https://pcingola.github.io/SnpEff/), and currently only support the human hg19 and hg38 genomes. Additionally, the workflow will add the [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) annotations for the SNP variants.




## Troubleshooting

<!---Any additional tips.--->
+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).



## FAQ's

<!---Frequently asked questions, pose any known limitations as FAQ's.--->

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-template/issues) page or start a discussion on the [community](https://community.nanoporetech.com/).



## Related blog posts

+ [Importing third-party workflows into EPI2ME Labs](https://labs.epi2me.io/nexflow-for-epi2melabs/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.



