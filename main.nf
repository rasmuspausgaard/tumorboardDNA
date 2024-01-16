#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

date=new Date().format( 'yyMMdd' )
user="$USER"
runID="${date}.${user}"


params.rundir                           ="${launchDir.baseName}" 
params.gatkTEMP                         ="${launchDir.baseName}/gatkTEMP"
params.server                           ="lnx01"
params.genome                           ="hg38" 
params.outdir                           ="TN_WES_results"
params.panel                            ="WES_2"    // set ROI to full WES
params.fastq                            =null
params.cram                             =null
params.fastqInput                       =null
params.help                             =false
params.pcgr_tumor                       =null
params.qualimap                         =null
params.samplesheet                      =null
params.jointgeno                        =null
params.hg38v1                           =null
params.hg38v2                           =null
params.skipQC                           =null
params.archiveStorage                   =null
params.gatk                             ="new"
//outdir_full_path= "${launchDir}/${params.outdir}/"

runtype = "TN_WES"



switch (params.server) {
    case 'lnx01':
        modules_dir="/home/mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules";
        subworkflow_dir="/home/mmaj/scripts_lnx01/nextflow_lnx01/dsl2/subworkflows";
        dataArchive="/lnx01_data2/shared/dataArchive";        
    break;
    case 'kga01':
        modules_dir="/home/mmaj/LNX01_mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules";
        subworkflow_dir="/home/mmaj/LNX01_mmaj/scripts_lnx01/nextflow_lnx01/dsl2/subworkflows";
        dataArchive="/data/shared/dataArchive";
    break;
}



def helpMessage() {
    log.info"""

    Generel info:
    Requires a samplesheet containing 5 columns in specific order (tab separated), without headerline:
    1) caseID, 2) NPN normal WES, 3) NPN tumor WES, 4) NPN tumor RNA, 5) PCGR tumor value

    Example samplesheet:

    johnDoe 112217976652	111184925465	111184925473    23

    The script will automatically look for fastq or ubam files in subfolders at /lnx01_data2/shared/dataArchive/. This location contains read-only access to the data archive, containing all FastQ files. Theres no need to copy or move any raw data (FastQ or uBAM).

    The user can point to a specific folder containing raw data (FastQ) using the --fastq option 
    This is only needed if raw data only exists outside the data archive (e.g. if data are in personal folders).
    
    NEW ANALYSIS SETUP: UMI pipeline using the --umibam option.
    using the --umibam option, the user can point to a folder containing unmapped bam (uBAM) files. This analysis uses the UMI information to generate consensus reads based on the UMIs. The UMI pipeline does not use fastq files at all.

    Usage:

    Main options:
      --help                print this help message
      
      --genome              hg19 or hg38
                                Default: hg38

      --outdir              Select which folder to write output to.
                                Default: TN_WES_results

      --samplesheet         path to case samplesheet. Can contain multiple patients/cases (each case in a separate line). 

      --server              Select which server the analysis is performed on (kga01 or lnx01)
                                Default: lnx01

      --fastq               Path to dir with fastq files
                                Default: data storage dirs at lnx01 server or kga01 server

    """.stripIndent()
}
if (params.help) exit 0, helpMessage()


def errorMessage1() {

    log.info"""

    USER INPUT ERROR: If no samplesheet is selected, the user needs to point to a folder containing relevant fastq or CRAM files... 
    Run the script with the --help parameter to see available options
    
    """.stripIndent()
}

if (!params.samplesheet && !params.fastq && !params.cram) exit 0, errorMessage1()

def errorMessage2() {

    log.info"""

    USER INPUT ERROR: Choose either fastq or CRAM as input... Not both. 
    Run the script with the --help parameter to see available options
    
    """.stripIndent()
}

if (params.fastq && params.cram) exit 0, errorMessage2()


///////////////////////////// SAMPLESHEET CHANNELS /////////////////////////////

// Samplesheet cols (fixed order)
// 0: CaseID, 1: WES.blood, 2: WES.tumor, 3: RNA tumor, 4: pcgr_tumor_code

////////////////////////////////////////////////////////////////////////////////

channel.fromPath(params.samplesheet)
    .splitCsv(sep:'\t')
    .map { row -> tuple(row[1], row[0])}
    .set { normalID_caseID }
//above: Normal sampleID (NPN), caseID

channel.fromPath(params.samplesheet)
    .splitCsv(sep:'\t')
    .map { row -> tuple(row[2], row[0])}
    .set { tumorID_caseID }

//above: tumor sampleID (NPN), caseID

channel.fromPath(params.samplesheet)
    .splitCsv(sep:'\t')
    .map { row -> tuple(row[0], row[4])} 
    .set { caseID_pcgrID }



channel.fromPath(params.samplesheet)
    .splitCsv(sep:'\t')
    .map { row -> tuple(row[0], row[1]+"_EV8")}
    .set { caseID_normalID } // use for Mutect2 --normal

///////////////// END: SAMPLESHEET CHANNELS ////////////////////////


////////////////// INPUT DATA (FASTQ) CHANNELS ///////////////////

if (params.fastq) {
    params.reads = "${params.fastq}/**{.,-}{EV8}{.,-}*R{1,2}*{fq,fastq}.gz"
}


if (!params.cram && !params.fastq && params.fastqInput) {
    params.reads="${dataArchive}/{lnx01,kga01_novaRuns,tank_kga_external_archive}/**/*{.,-}{EV8}{.,-}*R{1,2}*{fq,fastq}.gz"
}

if (!params.cram && params.fastqInput) {
    channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .map { it -> [it[0], file(it[1][0]),file(it[1][1])] }
    .set { read_pairs_ch }
    // above sampleID, r1, r2
    normalID_caseID
    .join(read_pairs_ch)
    .map {tuple(it[1],it[0]+"_EV8",it[2],it[3],"NORMAL")}
    .set { NN1 }
    //above: caseid, NPN_sampletype(NORMAL), normal R1, normal R2

    tumorID_caseID
    .join(read_pairs_ch)
    .map {tuple(it[1],it[0]+"_EV8",it[2],it[3],"TUMOR")}
    .set { CF1 }
    //above: caseid, sampleID, ,R1, R2, type

    NN1.concat(CF1)
    .set { case_fastq_input_ch }
    //above: NN2 and CF2 in the same channel (same structure as NN2 and CF2)

    //case_fastq_input_ch.view()

}

////////////////// INPUT DATA (CRAM) CHANNELS ///////////////////

if (!params.cram && !params.fastqInput && !params.fastq) {
    cramfiles="${dataArchive}/{lnx01,kga01_novaRuns,tank_kga_external_archive}/**/*{_,-}{EV8}*.cram"
    craifiles="${dataArchive}/{lnx01,kga01_novaRuns,tank_kga_external_archive}/**/*{_,-}{EV8}*.crai"
}

if (params.cram ) {
    cramfiles="${params.cram}/*{_,-}{EV8}*.cram"
    craifiles="${params.cram}/*{_,-}{EV8}*.crai"
}

if (!params.fastqInput) {
    Channel
    .fromPath(cramfiles)
    .map { tuple(it.baseName.tokenize('_').get(0),it) }
    .set { sampleID_cram }
    // above: sampleID, sampleCRAM
    Channel
    .fromPath(craifiles)
    .map { tuple(it.baseName.tokenize('_').get(0),it) }
    .set { sampleID_crai }
    // above: sampleID, sampleCRAI

    // Join with samplesheet:
    normalID_caseID // sampleID normal, caseID
    .join(sampleID_cram).join(sampleID_crai)
    .map {tuple(it[1],it[0]+"_EV8", it[2],it[3],"NORMAL")}
    .set { cram_normal }
    //above structure: caseID, NPN_EV8, CRAM, CRAI, NORMAL
    
    tumorID_caseID
    .join(sampleID_cram).join(sampleID_crai)
    .map {tuple(it[1],it[0]+"_EV8",it[2],it[3],"TUMOR")}
    .set { cram_tumor }
    //above structure: caseID, NPN_EV8, CRAM, CRAI, TUMOR
    
    cram_normal.concat(cram_tumor)
    .set { case_npn_cram_crai_ch }
    // caseID, NPN, CRAM, CRAI

     case_npn_cram_crai_ch
    .filter{it =~ /NORMAL/}
    .set { normals_ch }

    case_npn_cram_crai_ch 
    .filter{it =~ /TUMOR/}
    .set { tumor_ch }
    
    normals_ch
    .join(tumor_ch)
    .set { tumorNormal_cram_ch } 

}


log.info """\

========================================================
KGA Vejle paired tumor-normalWES nextflow pipeline v3
========================================================
results     : $params.outdir
user        : $user
rundir      : $params.rundir
runtype     : $runtype
runID       : $date.$user
"""




 //   {msisensor_input; mutect2_input; sequenza_input; accucopy_input;tumor_normal_bams4;facets_input}


include { 

         inputFiles_symlinks_cram;
         tb_cram_bam;
         tb_haplotypecaller;
         SUB_DNA_PREPROCESS;
         SUB_DNA_QC;
         SUB_DNA_TUMOR_NORMAL
        } from "/data/shared/analyseScripts/modules/tumorBoard.modules.v1.nf" 

//from "./modules/tumorBoard.modules.v1.nf"


workflow DNA_TUMOR_NORMAL {
    take:
    tumorNormal_bam_ch
    main:
    mutect2(tumorNormal_bam_ch)
    msisensor(tumorNormal_bam_ch)
    sequenza(tumorNormal_bam_ch)
    sequenza_R_output(sequenza.out)
    pcgr_v103(mutect2.out.mutect2_tumorPASS.join(caseID_pcgrID))
    pcgr_v141(mutect2.out.mutect2_tumorPASS.join(caseID_pcgrID))
    emit:    
    mutect2_out=mutect2.out.mutect2_vcf

}



workflow {
    if (params.fastqInput) {

        SUB_DNA_PREPROCESS(case_fastq_input_ch)
        
        if (!params.skipQC) {
            SUB_DNA_QC(SUB_DNA_PREPROCESS.out.finalBam)
        }
        
        tb_haplotypecaller(SUB_DNA_PREPROCESS.out.finalBam)

        SUB_DNA_PREPROCESS.out.finalBam
        .filter{it =~ /NORMAL/}  
        .set {normal_ch }

        SUB_DNA_PREPROCESS.out.finalBam
        .filter{it =~ /TUMOR/}  
        .set {tumor_ch }

        normal_ch.join(tumor_ch)
        .set { tumorNormal_bam_ch }

        SUB_DNA_TUMOR_NORMAL(tumorNormal_bam_ch, caseID_pcgrID)
    }

    if (!params.fastqInput) {
        inputFiles_symlinks_cram(case_npn_cram_crai_ch)
        tb_haplotypecaller(case_npn_cram_crai_ch)  // caseid, npn, cram, crai, type
        tb_cram_bam(case_npn_cram_crai_ch)
        
        if (!params.skipQC) {
            SUB_DNA_QC(tb_cram_bam.out.bam)
        }
        
        tb_cram_bam.out.bam
        .filter{it =~ /NORMAL/}
        .set { bam_normals_ch }

        tb_cram_bam.out.bam
        .filter{it =~ /TUMOR/}
        .set { bam_tumor_ch }

        bam_normals_ch
        .join(bam_tumor_ch)
        .set { tumorNormal_bam_ch }   
        // above structure: tuple val(caseID), val(sampleID_normal), path(bamN), path(baiN),val(typeN), val(sampleID_tumor),path(bamT), path(baiT),val(typeT)

        SUB_DNA_TUMOR_NORMAL(tumorNormal_bam_ch, caseID_pcgrID)
    }
}

/*


//------------------------------------------------------------------------//
//---------------------- Tumor-normal based analysis ---------------------//
//-----------------------------------------------------------------------//


process cfDNA_facets_snp_pileup {
    errorStrategy 'ignore'
    tag "$caseID"
    publishDir "${caseID}/${params.outdir}/facets", mode: 'copy'

    input:
    tuple val(caseID), val(sampleID_normal), path(bamN), path(baiN), val(sampleID_tumor),path(bamT), path(baiT) from facets_input
   
    output:
    tuple val(caseID), path("*.snp_pileup.gz") into facets_snppileup
   
    when:
    !params.germline_only
   
    script:
    """
    singularity run -B ${s_bind} ${simgpath}/facetssuite.sif snp-pileup-wrapper.R \
    -vcf ${dbsnp} \
    -n ${bamN} \
    -t ${bamT} \
    -o ${caseID}.facets
    """
}


process facets_output {
    errorStrategy 'ignore'
    tag "$caseID"
    publishDir "${caseID}/${params.outdir}/facets", mode: 'copy'

    input:
    tuple val(caseID), path(facets_pileup) from facets_snppileup

    output:
    path("${caseID}.facets/*")
    
    when:
    !params.germline_only
    
    script:
    """
    singularity run -B ${s_bind} ${simgpath}/facetssuite.sif run-facets-wrapper.R \
    -f ${facets_pileup} \
    -s ${caseID}.facets \
    -g ${params.genome} \
    -e -pc 1000 -c 500
    """
}



*/



