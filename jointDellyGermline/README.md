<img src="https://oicr.on.ca/wp-content/themes/oicr/assets/img/logo.svg" width=25% height=25%>

# DELLY Germline SV and CNV Calling  

### Introduction
DELLY performs germline structural variant (SV) and copy number variation (CNV) calling using a comprehensive tool designed for the detection of germline SVs and CNVs from next-generation sequencing data. Leveraging sophisticated algorithms, including discordant read pairs, split reads, and read depth analysis, DELLY identifies a wide range of genomic variations, including deletions, duplications, inversions, amplifications, and deletions. The Joint Delly Germline procedure takes in whole genome bam files and vcf files from DELLY, for a set of samples (ideally 20+ unrelated samples) and utilizes information from all samples to detect variants with greater sensitivity and accuracy.

### Summary
A predefined set of vcf, bam, and bai files are provided as input to the Joint Delly Germline workflow. This generates a joint SV vcf and a joint CNV vcf with information from all samples where a call is made in any one sample. 

### Inputs  
- VCF files: A set of VCF files that have been individually processed through DELLY for SV calling
- BAM files: A set of whole genome bam files for all samples to be genotyped
- BAI files: A set of bam index files to accompany all bam file inputs
- genome resources: as specified in the WDL workflow

### Procedure
 1. Identify the set of samples you want to perform germline SV and CNV calling on
 2. Construct an olive that will pull these files as inputs to the joint genotyping workflow. This workflow is available in vidarr stage. The olive can be set up as an .actnow. The olive on this page can be modified to specify the project, with other modifications to limit the input data as needed; for example we may only process cases that have been released up to a certain time point. An olive template "jointDellyGermline.shesmu" can be found in this repo.
 3. The olive will invoke the dellyGermline workflow (https://github.com/oicr-gsi/dellyGermline), which will run the joint genotyping tool, and produce two vcf output files - one with SV calls across all samples being analyzed, and one with CNV calls across all samples being analyzed.
 
### Deliverables
 - Germline SV vcf
 - Germline CNV vcf

### Static action running mode
Joint Delly Germline genotyping will be performed for specific projects at certain timepoints, and thus it will not be initiated from a shesmu olive that runs actions continuously. Instead, a .actnow file can be used to generate static actions (https://github.com/oicr-gsi/shesmu/blob/master/actnow.md). The .shesmu file in this repository can be updated to target a specific project or subset of samples, and the edited script can be loaded into shesmu's simulator to generate a .actnow file. To do this, simply simulate the actions and then click "download" on the actions page. This .actnow olive can then be added to shesmu stage, where shesmu will identify the file and generate only the static actions specified in the .actnow olive.

### Resources
- DELLY : https://github.com/dellytools/delly
- OICR delly Workflow : https://github.com/oicr-gsi/delly
- OICR bamMergePreprocessing Workflow: https://github.com/oicr-gsi/bamMergePreprocessing
- OICR dellyGermline Workflow : https://github.com/oicr-gsi/dellyGermline
