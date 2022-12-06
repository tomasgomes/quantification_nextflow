#!~/bin/miniconda3/envs/nextflowenv/bin nextflow

/*
 * Example command:
 *
 * nextflow run kallisto_velocity_pipeline.nf --transcriptome /links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/introns/AmexT_v47_artificial_introns.fa --transindex AmexT_v47_artificial_introns.kalid --t2g /links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/introns/AmexT_v47_artificial_introns_genenames_t2g.txt --white /links/groups/treutlein/USERS/tomasgomes/gene_refs/other/10xv3_whitelist.txt --samplename "a1_2_GEX" --outdir /local1/USERS/tomasgomes/kallisto_velocity/ --protocol 10xv3 --reads "/links/groups/treutlein/DATA/sequencing/20210518_P1567_ASHLEY_10X_axolotl_whole_pallium_multiome/raw/a1_2_GEX/*.fastq.gz"
 *
 */



/*
 * Defines necessary parameters
 * reads is used for path for fastq files, but also for the batch file for plate-based data
 */
params.samplename = false
//params.samplename = "GER006_10x"
params.outdir = "./"
params.protocol = false
//params.protocol = "10xv2"
params.transcriptome = false
//params.transcriptome = "/links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/cDNA_transcripts/AmexT_v47_artificial.fa"
params.transindex = false
//params.transindex = "AmexT_v47_artificial.kalid"
params.reads = false
//params.reads = "/links/groups/treutlein/USERS/tomasgomes/projects/axolotl/data/raw/Gerber_all10x/GER006_10x/*.fastq.gz" 
params.white = false
//params.white = "/links/groups/treutlein/USERS/tomasgomes/gene_refs/other/10xv2_whitelist.txt"
params.t2g = false
//params.t2g = "/links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/cDNA_transcripts/AmexT_v47_artificial_genenames_t2g.txt"
params.intronsfile = "/links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/introns/AmexT_v47.introns_list.txt"
params.exonsfile = "/links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/introns/AmexT_v47_artificial_list.txt"



/*
 * Input parameters validation and preparation
 */
//----- TRANSCRIPTOME -----//
//if(!params.transindex) exit 1, "Transcriptome index file path is mandatory (will be created if it does not exist)."
transindex_file = file(params.transindex)
if(!transindex_file.exists() and params.transcriptome!="") transcriptome_file = file(params.transcriptome)
else exit 1, "Missing transcriptome file or transcriptome index file."

//----- READS -----//
R1files = Channel
    .from(params.reads.tokenize())
    .flatMap{ files(it) }
    .filter{it =~ /.*R1.*/}
    .toSortedList()
    .flatten()
    //.view()
R2files = Channel
    .from(params.reads.tokenize())
    .flatMap{ files(it) }
    .filter{it =~ /.*R2.*/}
    .toSortedList()
    .flatten()
    //.view()
R1files
    .merge(R2files){ a, b -> tuple(a,b) }
    .flatten()
    //.view()
    .set{read_files_kallisto}
if(params.protocol=='plate'){
    Channel.fromPath(params.reads)
        .set{batch_kal}
} else{ batch_kal = "" }

//----- OTHER FILES -----//
if(!params.samplename) exit 1, "Please provide a name for this sample"

if(!params.protocol) exit 1, "Please provide an adequate protocol"

if(!params.white && params.protocol!='plate'){
    exit 1, "Barcode whitelist is mandatory for Chromium and Visium runs."
} else if (params.protocol!='plate'){
    Channel.fromPath(params.white)
        .set{bc_wl_kal}
} else{ bc_wl_kal = "" }

if(!params.t2g && params.protocol!='plate'){
    exit 1, "Transcriptome-to-gene reference is required for quantification."
} else if(params.protocol!='plate'){
    Channel.fromPath(params.t2g)
        .set{t2g_kal}
    t2g_plate = ""
} else{
    Channel.fromPath(params.t2g)
        .set{t2g_plate}
    t2g_kal = ""
}

if(file("${params.outdir}${params.samplename}").exists()) println "Output folder already exists. No files will be overwritten, but execution may fail."

Channel.fromPath(params.intronsfile)
        .set{intronsfile_l}
Channel.fromPath(params.exonsfile)
        .set{exonsfile_l}


/*
 * Step 0. For some reason pseudoalignment doesn't work if it is the first step
 */
process starting{
 
    script:
    """
    echo Pseudoalignment doesn't work as a first step, so we'll have this here.
    """
}
 
/*
 * Step 1. Builds the transcriptome index, if it doesn't exist
 */
process index {

    // when index does not exist, it is stored in the transcriptome dir
    //publishDir transcriptome_file.getParent(), mode: 'copy'
    storeDir transcriptome_file.getParent()
    
    input:
    file transcriptome_file
    
    output:
    file params.transindex into transcriptome_index
  
    script:
    """
    kallisto index -i ${params.transindex} -k 31 --make-unique ${transcriptome_file}
    """
}

/*
 * Step 2. Do pseudoalignment with kallisto
 */
process pseudoal {

    //publishDir params.outdir, mode: 'copy'
    storeDir params.outdir

    input:
    file index from transcriptome_index
    file reads from read_files_kallisto.collect()

    output:
    file params.samplename into kallisto_bus_to_sort
    file "${params.samplename}/output.bus"
    file "${params.samplename}/matrix.ec"
    file "${params.samplename}/transcripts.txt"
    
    when: params.protocol!='plate'

    script:
    if(params.protocol=='visiumv1')
        """
        kallisto bus \\
            -i $index \\
            -o ${params.samplename} \\
            -x 0,0,16:0,16,28:1,0,0 \\
            -t 24 \\
            $reads
        """
    else
        """
        kallisto bus \\
            -i $index \\
            -o ${params.samplename} \\
            -x ${params.protocol} \\
            -t 24 \\
            $reads
        """
}

/*
 * Step 3. Correct barcodes and sort
 * A future version of kallisto may include options to correct the UMI
 */
process corrsort {

    //publishDir params.outdir, mode: 'copy'
    storeDir params.outdir

    input:
    file outbus from kallisto_bus_to_sort
    file white from bc_wl_kal.collect()

    output:
    file "${outbus}/output.cor.sort.bus"
    file outbus into kal_sort_to_capt

    when: params.protocol!='plate'

    script:
    """
    bustools correct -w $white -o ${outbus}/output.cor.bus ${outbus}/output.bus
    bustools sort -o ${outbus}/output.cor.sort.bus -t 20 ${outbus}/output.cor.bus
    """
}

/*
 * Step 4. Subset introns and exons
 */
process captureInEx {

    //publishDir params.outdir, mode: 'copy'
    storeDir params.outdir

    input:
    file outbus from kal_sort_to_capt
    file intronsfile from intronsfile_l.collect()
    file exonsfile from exonsfile_l.collect()

    output:
    file "${outbus}/spliced.bus"
    file "${outbus}/unspliced.bus"
    file outbus into kal_capt_to_count
    
    when: params.protocol!='plate'

    script:
    """
    bustools capture -s -x -o ${outbus}/spliced.bus -c $exonsfile -e ${outbus}/matrix.ec -t ${outbus}/transcripts.txt ${outbus}/output.cor.sort.bus
bustools capture -s -x -o ${outbus}/unspliced.bus -c $intronsfile -e ${outbus}/matrix.ec -t ${outbus}/transcripts.txt ${outbus}/output.cor.sort.bus
    """
}

/*
 * Step 5. Obtain the counts
 */
process countbus {
    //publishDir params.outdir, mode: 'copy'
    storeDir params.outdir

    input:
    file outbus from kal_capt_to_count
    file t2g from t2g_kal.collect()

    output:
    file "${outbus}/unspliced.mtx"
    file "${outbus}/spliced.mtx"

    when: params.protocol!='plate'

    script:
    """
    bustools count --em -g $t2g -e ${outbus}/matrix.ec -t ${outbus}/transcripts.txt --genecounts -o ${outbus}/unspliced ${outbus}/unspliced.bus
    bustools count --em -g $t2g -e ${outbus}/matrix.ec -t ${outbus}/transcripts.txt --genecounts -o ${outbus}/spliced ${outbus}/spliced.bus
    """
}




