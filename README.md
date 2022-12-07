File description:  

 - `kallisto_pipeline.nf` - the main quantification pipeline  
 - `kallisto_pipeline_SRR.nf` - quantification pipeline for files downloaded from SRA  
 - `kallisto_velocity_pipeline.nf` - the quantification pipeline for RNA velocity  
  
  
Notes on pipeline:  

 - while `kallisto_pipeline.nf` detects read1 and read2 as all files in a folder that have "R1."/"R2.", `kallisto_pipeline_SRR.nf` detects the same reads as files with "_2."/"_3.". This is the only different between these pipelines  
 - reads are passed as a list of all fastq files in a directory (e.g. "*.fastq.gz"), so ideally all reads to be used in the quantification of a given sample should be in the same directory  
 - currently, the %MT detection assumes a gene format like "geneID-genename", and tries to detect MT genes based on predefined gene names  
 - the transcriptome index is created in the same folder where the fasta file is located. if it has been created in advance, you still need to pass a fasta file path, as well as the index file name  
 - running the pipeline for Visium data also requires installing spaceranger. this is necessary to do a mock run, where the image outputs are retrieved from  
 - at the moment, the mock spaceranger run relies on a quantification against a human transcriptome index THAT IS HARDCODED IN THE CODE AND NEEDS TO BE ADAPTED   
   
  
Example commands:  

 - 10xv3:  
 ```
 nextflow run kallisto_pipeline.nf --transcriptome /links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/full_transcript/Amex_Tv47_Gv6.full_transcript_artificial.fa --transindex Amex_Tv47_Gv6.full_transcript_artificial.kalid --t2g /links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/full_transcript/AmexT_v47_artificial_genenames_t2g.txt --white /links/groups/treutlein/USERS/tomasgomes/gene_refs/other/737K-arc-v1.txt --samplename "a1_1_GEX" --outdir /links/groups/treutlein/USERS/tomasgomes/data/axolotl/ --protocol 10xv3 --reads "/local1/DATA/sequencing/20210518_P1567_ASHLEY_10X_axolotl_whole_pallium_multiome/raw/a1_1_GEX/*.fastq.gz"
 ```
 
 - Smart-seq2 data (requires preparing data table):  
 ```
 nextflow run kallisto_pipeline.nf --transcriptome /links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/cDNA_transcripts/AmexT_v47_artificial.fa --transindex AmexT_v47_artificial.kalid --t2g /links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/cDNA_transcripts/AmexT_v47_artificial_genenames_t2g.txt --samplename "Gerber_plate" --outdir /links/groups/treutlein/USERS/tomasgomes/data/axolotl/ --protocol plate --reads /links/groups/treutlein/USERS/tomasgomes/projects/axolotl/data/raw/Gerber_allcells/kallisto_batch.txt
 ```
 
 - Visium:  
 ```
 nextflow run kallisto_pipeline.nf --transcriptome /links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/cDNA_transcripts/AmexT_v47_artificial.fa --transindex AmexT_v47_artificial.kalid --t2g /links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/cDNA_transcripts/AmexT_v47_artificial_genenames_t2g.txt --white /links/groups/treutlein/USERS/tomasgomes/gene_refs/other/visium-v1_whitelist_kallisto.txt --samplename "A1_limb" --outdir /links/groups/treutlein/USERS/tomasgomes/data/axolotl/ --protocol visiumv1 --reads "/links/groups/treutlein/DATA/sequencing/20200821_P1288_ASHLEY_VISIUM_axolotl_visium_control_11dpa/raw/A1_Animal1_Control/*.fastq.gz" --images "V19S23-109" --imagear "A1" --imagef "/links/groups/treutlein/DATA/sequencing/20200821_P1288_ASHLEY_VISIUM_axolotl_visium_control_11dpa/image/A1_large_image1.jpg" --imageal "/links/groups/treutlein/DATA/sequencing/20200821_P1288_ASHLEY_VISIUM_axolotl_visium_control_11dpa/alignment_files/V19S23-109-A1.json"
 ```
 
 - RNA velocity:  
 ```
 nextflow run kallisto_velocity_pipeline.nf --transcriptome /links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/introns/AmexT_v47_artificial_introns.fa --transindex AmexT_v47_artificial_introns.kalid --t2g /links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/introns/AmexT_v47_artificial_introns_genenames_t2g.txt --white /links/groups/treutlein/USERS/tomasgomes/gene_refs/other/10xv3_whitelist.txt --samplename "a1_2_GEX" --outdir /local1/USERS/tomasgomes/kallisto_velocity/ --protocol 10xv3 --reads "/links/groups/treutlein/DATA/sequencing/20210518_P1567_ASHLEY_10X_axolotl_whole_pallium_multiome/raw/a1_2_GEX/*.fastq.gz"
 ```
  
  
Files in `other_files/`:  

 - `instructions_Amex_transcriptome.txt` - (older) instructions on how to format the axolotl transcriptome to fit with kallisto, as well as overall use of the pipeline  
 - `nextflowenv.yaml` - conda environment used to run the pipeline  
 - `kallisto_batch.txt` - example batch file for Smart-seq2  
 - `barcode_whitelists/` - several files with barcode whitelists vor various 10x protocols (compressed) 
 - `axolotl_genome_cut.config` - an example config file to index the axolotl (cut) genome for cellranger-ARC
 - `chrNameLength.txt` - example file with chromosome name and length (a STAR indexing output)
 

Files in `other_scripts/`:  

 - `convertCoords.py` - convert chromosome coordinates in a BED file from a chopped down genome back into the original
  - example: `python3 convertCoords.py chrNameLength.txt input.bed output.bed`
 - `convertCoords_indiv.py` - convert specific chromosome coordinates from a chopped down genome back into the original
  - example: `python3 convertCoords_indiv.py chrNameLength.txt chr1ps3:213318726-213320220`
 - `cutChroms.py` - script to cut chromosomes into a specific maximum length. also prevents cutting in the middle of genes, including giving additional space around its ends
  - requires preparation:
  ```
  # split genome into individual chromosomes
  faidx --split-files ../AmexG_v6.DD.corrected.round2.chr.fa
  # return chromosome sizes
  faidx ../AmexG_v6.DD.corrected.round2.chr.fa -i chromsizes > AmexG_v6.DD.corrected.round2.chr.sizes
  ```
  - example: `python cutChroms.py AmexG_v6_chr/AmexG_v6.DD.corrected.round2.chr.sizes ../AmexT_v47.FULL_corr_chr.gtf AmexG_v6_chr/AmexT_v47.FULL_corr_chr_cut.gtf AmexG_v6_chr/cutcoord.txt AmexG_v6_chr 500000000 2000`
 - `get_t2g_gtf.py` - script to extract gene ID and gene names from gtf file, and compile a `t2g.txt` file for kallisto. Returns gene ID and gene name versions, as well as a full annotation version.
  - this is taylored to the axolotl GTF file
  - example: `python get_t2g_gtf.py AmexT_v47.FULL_corr_chr.gtf t2g.txt t2g_genenames.txt t2g_note.txt`
 - `get_t2g_gtf_introns.py` - script to extract gene ID and gene names from gtf file, for use in a gtf file with introns (for kallisto RNA velocity)
  - example: `python get_t2g_gtf.py AmexT_v47.FULL_corr_chr.gtf AmexT_v47.introns_t2g.txt AmexT_v47.introns_genenames_t2g.txt AmexT_v47.introns_t2g_note.txt AmexT_v47.introns_list.txt`
 - `getTranscriptFastaHead.py` - simplifies fasta header to only have gene ID
  - example: `python getTranscriptFastaHead.py AmexT_v47_CDS_reformatted.fa AmexT_v47.FULL.fa`
 - `getExonFastaHead.py` - simplifies fasta header to only have gene ID, for a fasta with exon sequences. each gene name will get an added number (.1, .2, ...) depending on the exon number
  - I couldn't find the original use case for this...
  - example: `python getExonFastaHead.py AmexT_v47_CDS_reformatted.fa AmexT_v47.FULL.fa exon_list.txt`
 - `getExont2g.py` - make a `t2g.txt` file from a GTF, for exons
  - example: `python getExont2g.py AmexT_v47.FULL_corr_chr.gtf exons_t2g.txt`
 - `list_annotations.py` - create a new file with gene annotations
  - I couldn't find the original use case for this, and not entirely sure what the input and output is...
 - `simplify_gtf_ids.py` - simplify the gene annotations in a gtf file
  - I couldn't find the original use case for this
 - `simplifyExonsDesc.py` - simplify exon description
  - I couldn't find the original use case for this, and not entirely sure what the input and output is...
 
 
 
