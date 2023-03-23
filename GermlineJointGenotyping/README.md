<img src="https://raw.githubusercontent.com/oicr-gsi/analysisProcedures/resources/oicr_logo.png" width=25% height=25%>

# Germline Joint Genotyping   

### Introduction
Analysis of Whole Genome or Exome sequence for Germline variation is handled by the Haplotype caller workflow. Haplotype caller identified potential variant sites in each sample and saves the results as a genomic vcf (gvcf) with non-variant sites represented as blocks. This allows the generation of mutisample vcf files with each site showing information from all samples. The Germline Joint Genotyping procedure takes in gvcf files from a set of samples and utilizes information from all samples to detect variants with greater sensitivity and accuracy.

### Summary
A predefined set of gvcf files are provided as input to the joint Genotyping workflow. This generates a single, multi-sample vcf files with information from all samples where a call is made in any one sample.  Variant score recalibration produces a more accurate assessment of the germline variation.

### Inputs  
- gVCF files: A set of gVCF files across which joint genotypes are to be called.
- genome resources: as specified in the WDL workflow.

### Procedure
 1. Identify the set of gvcf files produced by the haplotype caller workflow
 2. Construct an olive that will pull these files as inputs to the joint genotyping workflow. This workflow is available in vidarr stage. The olive can be set up as an .actnow.
 3. Collect and release data
 
### Deliverables
 - VCF file.
 - the raw joint genotyped vcf + tbi index
 - the recalibrated joint gentotyped vcf + tbi index
 - a calling summary metrics file, a single summary record for the joint calling
 - a calling detailed metrics file, one record per input vcf

### Resources
- GATK Haplotype Caller : https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
- workflow : https://github.com/oicr-gsi/gatk-haloptype-caller
- GVCF File Format : https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format
- GATK GenotypeGVCFs : https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs
- workflow : https://github.com/oicr-gsi/gatk-genotype-GVCFs
- Joint Genotyping : https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants
- VQSR : https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR-
- VCF specification : https://samtools.github.io/hts-specs/VCFv4.2.pdf


