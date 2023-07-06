<img src="https://oicr.on.ca/wp-content/themes/oicr/assets/img/logo.svg" width=25% height=25%>

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
 2. Construct an olive that will pull these files as inputs to the joint genotyping workflow. This workflow is available in vidarr stage. The olive can be set up as an .actnow.  The olive on this page can be modified to specify the project, with other modifications to limit the input data as needed, for example when we only process certain cases that have been released to the time point. An olive template "JointGenotyping.shesmu" can be found on this repo.
 3. The olive will invoke the running of genotypeGVCFs workflow (https://github.com/oicr-gsi/gatk-genotype-GVCFs/blob/master/genotypeGVCFs.wdl), which will run the joint genotyping tool, GenotypeGVCFs, and also run tools for consolidate GVCFs and VariantRecalibrator (VQSR) included in gatk best practise. Once done processing, the workflow will produce a joint genotyped, variant quality score recalibrated outputs. 
 4. Annotate the workflow with VEP. The vcf output of step 3 is a single recalibrated.vcf.gz file, use this as input for variantEffectPredictor. This should be done using the WDL file and can likely be set up as another actnow, since this will only need to run once for one joint analysis. 
 3. Collect and release data.
 
### Deliverables
 - VCF file.
 - the raw joint genotyped vcf + tbi index
 - the recalibrated joint gentotyped vcf + tbi index
 - a calling summary metrics file, a single summary record for the joint calling
 - a calling detailed metrics file, one record per input vcf
 - md5sum for the result files

### Static action running mode
Since joint genetyping usually performed for specific project at certain time points, it wouldn't be desirable to use olive for pipeline running mode, which will generate actions continuously. We can use .actnow file to generate static actions (https://github.com/oicr-gsi/shesmu/blob/master/actnow.md). 
For this purpose write a shesmu file, but don't push this .shesmu file as olive to shesmu stage, instead use it in simulator to generate .actnow file and push the .actnow file to shesmu stage ( where we usually put olive). shesmu will find this file but only generate static actions specified in .actnow.
Before this .vidarrworkflow needs to be set up as usual (for example genotypeGVCFs.vidarrworkflow under vidarr/stage).

### Resources
- GATK Haplotype Caller : https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
- workflow : https://github.com/oicr-gsi/gatk-haloptype-caller
- GVCF File Format : https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format
- GATK GenotypeGVCFs : https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs
- workflow : https://github.com/oicr-gsi/gatk-genotype-GVCFs
- Joint Genotyping : https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants
- VQSR : https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR-
- VCF specification : https://samtools.github.io/hts-specs/VCFv4.2.pdf


