# eqtl
## Generate normalized expression in BED format

The expression data are normalized as follows:

1. Read counts are normalized between samples using TMM.

2. Genes are selected based on the following expression thresholds:
 - $\geq 0.1$ TPM in $\geq 20\%$ samples AND
 - $\geq 6$ reads (unnormalized) in $\geq 20\%$ samples

3. Each gene is inverse normal transformed across samples.

```{r,eval=FALSE}
eqtl_prepare_expression.py ${tpm_gct} ${counts_gct} ${annotation_gtf}
${sample_participant_lookup}$ ${vcf_chr_list}$ ${prefix}$
--tpm_threshold 0.1 
--count_threshold 6 
--sample_frac_threshold 0.2 
--normalization_method tmm
```

The file ${vcf_chr_list} lists the chromosomes in the VCF, and can be generated using

tabix --list-chroms {vcf} > {vcf_chr_list}
The file ${sample_participant_lookup} must contain two columns, sample_id and participant_id, mapping IDs in the expression files to IDs in the VCF (these can be the same).

This step generates the following BED file and index:

```{r,eval=FALSE}
${prefix}.expression.bed.gz
${prefix}.expression.bed.gz.tbi
```

## Calculate PEER factors
```{r,eval=FALSE}
Rscript run_PEER.R ${prefix}.expression.bed.gz ${prefix} ${num_peer}
```

The number of PEER factors was selected as function of sample size (N):

 - 15 factors for N $< 150$
- 30 factors for 150 $\leq$ N $<$ 250
- 45 factors for 250 $\leq$ N $<$ 350
- 60 factors for N $\geq$ 350

This step will generate 3 files:

```{r,eval=FALSE}
${prefix}.PEER_residuals.txt
${prefix}.PEER_alpha.txt
${prefix}.PEER_covariates.txt
```

## Combine covariates

This step generates a combined covariates file, containing genotype PCs, PEER factors, and additional explicit covariates (e.g., genotyping platform).

```{r,eval=FALSE}
combine_covariates.py ${prefix}.PEER_covariates.txt ${prefix}
    --genotype_pcs ${genotype_pcs}
    --add_covariates ${add_covariates}
```

The covariate files should have one covariate per row, with an identifier in the first column, and a header line with sample identifiers. This step will generate the file ${prefix}.combined_covariates.txt

## Run FastQTL
```{r,eval=FALSE}
# nominal pass
run_FastQTL_threaded.py ${vcf} ${prefix}.expression.bed.gz ${prefix} 
    --covariates ${prefix}.combined_covariates.txt 
    --window 1e6 --chunks 100 --threads 16

# permutation pass
run_FastQTL_threaded.py ${vcf} ${prefix}.expression.bed.gz ${prefix} 
    --covariates ${prefix}.combined_covariates.txt 
    --window 1e6 --chunks 100 --threads 16 
    --permute 1000 10000 
```

The following files will be generated:

```{r, eval=FALSE}
${prefix}.allpairs.txt.gz
${prefix}.egenes.txt.gz
```
