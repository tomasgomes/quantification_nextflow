File description:
 - `kallisto_pipeline.nf` - the main quantification pipeline
 - `kallisto_pipeline_SRR.nf` - quantification pipeline for files downloaded from SRA
 - `kallisto_velocity_pipeline.nf` - the quantification pipeline for RNA velocity


Notes on pipeline:
 - while `kallisto_pipeline.nf` detects read1 and read2 as all files in a folder that have "R1."/"R2.", `kallisto_pipeline_SRR.nf` detects the same reads as files with "_2."/"_3.". This is the only different between these pipelines
 - currently, the %MT detection assumes a gene format like "geneID-genename", and tries to detect MT genes based on predefined gene names
 - the transcriptome index is created in the same folder where the fasta file is located. if it has been created in advance, you still need to pass a fasta file path, as well as the index file name
 - running the pipeline for Visium data also requires installing spaceranger. this is necessary to do a mock run, where the image outputs are retrieved from
 - at the momentl, the mock spaceranger run relies on a quantification against a human transcriptome index THAT IS HARDCODED IN THE CODE AND NEEDS TO BE ADAPTED 
 

Example commands:
 - 10xv3:
 
 - Smart-seq2 data (requires preparing data table)
 
 - Visium
 

Files in `other_files/`:
 - `instructions_Amex_transcriptome.txt` - (older) instructions on how to format the axolotl transcriptome to fit with kallisto, as well as overall use of the pipeline
 - `nextflowenv.yaml` - conda environment used to run the pipeline
 - `kallisto_batch.txt` - example batch file for Smart-seq2
 - `barcode_whitelists/` - several files with barcode whitelists vor various 10x protocols (compressed)
 
 
 
 
