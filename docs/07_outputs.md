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
