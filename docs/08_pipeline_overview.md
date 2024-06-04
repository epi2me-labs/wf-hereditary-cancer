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
