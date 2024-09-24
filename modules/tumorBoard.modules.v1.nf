#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


date=new Date().format( 'yyMMdd' )
user="$USER"
runID="${date}.${user}"

multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"
//////////////////////////// SWITCHES ///////////////////////////////// 

switch (params.gatk) {

    case 'danak':
    gatk_image="gatk419.sif";
    break;
    case 'new':
    gatk_image="gatk4400.sif";
    break;
    case 'latest':
    gatk_image="gatk4500.sif";
    default:
    gatk_image="gatk4400.sif";
    break;
}


switch (params.server) {

    case 'lnx02':
        s_bind="/data/:/data/,/lnx01_data2/:/lnx01_data2/,/fast/:/fast/,/lnx01_data3/:/lnx01_data3/";
        simgpath="/data/shared/programmer/simg";
        params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/wgs_splitinterval_BWI_subdivision3/*.interval_list";
        tmpDIR="/fast/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/${gatk_image} gatk";
        multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"
        tank_storage="/home/mmaj/tank.kga/data/data.storage.archive/";
        //modules_dir="/home/mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules/";
    break;
    case 'lnx01':
        s_bind="/data/:/data/,/lnx01_data2/:/lnx01_data2/";
        simgpath="/data/shared/programmer/simg";
        params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/wgs_splitinterval_BWI_subdivision3/*.interval_list";
        tmpDIR="/data/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/${gatk_image} gatk";
        multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"
        tank_storage="/home/mmaj/tank.kga2/data/data.storage.archive/";
        modules_dir="/home/mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules/";
    break;
    case 'kga01':
        simgpath="/data/shared/programmer/simg";
        s_bind="/data/:/data/";
        tmpDIR="/data/TMP/TMP.${user}/";
        params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/wgs_splitinterval_BWI_subdivision3/*.interval_list";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/${gatk_image} gatk";
        tank_storage="/home/mmaj/tank.kga/data/data.storage.archive/";
        modules_dir="/home/mmaj/LNX01_mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules/";
    break;
}




switch (params.genome) {
    case 'hg19':
        assembly="hg19"
        // Genome assembly files:
        genome_fasta = "/data/shared/genomes/hg19/human_g1k_v37.fasta"
        genome_fasta_fai = "/data/shared/genomes/hg19/human_g1k_v37.fasta.fai"
        genome_fasta_dict = "/data/shared/genomes/hg19/human_g1k_v37.dict"
        genome_version="V1"
        break;


    case 'hg38':
        assembly="hg38"
        smncaller_assembly="38"
        // Genome assembly files:
        if (params.hg38v1) {
        genome_fasta = "/data/shared/genomes/hg38/GRCh38.primary.fa"
        genome_fasta_fai = "/data/shared/genomes/hg38/GRCh38.primary.fa.fai"
        genome_fasta_dict = "/data/shared/genomes/hg38/GRCh38.primary.dict"
        genome_version="V1"
        cnvkit_germline_reference_PON="/data/shared/genomes/hg38/inhouse_DBs/hg38v1_primary/cnvkit/wgs_germline_PON/jgmr_45samples.reference.cnn"
        cnvkit_inhouse_cnn_dir="/data/shared/genomes/hg38/inhouse_DBs/hg38v1_primary/cnvkit/wgs_persample_cnn/"
        inhouse_SV="/data/shared/genomes/hg38/inhouse_DBs/hg38v1_primary/"
        }
        
        if (params.hg38v2){
        genome_fasta = "/data/shared/genomes/hg38/ucsc.hg38.NGS.analysisSet.fa"
        genome_fasta_fai = "/data/shared/genomes/hg38/ucsc.hg38.NGS.analysisSet.fa.fai"
        genome_fasta_dict = "/data/shared/genomes/hg38/ucsc.hg38.NGS.analysisSet.dict"
        genome_version="V2"
        }

        // Current hg38 version (v3): NGC with masks and decoys.
        if (!params.hg38v2 && !params.hg38v1){
        genome_fasta = "/data/shared/genomes/hg38/GRCh38_masked_v2_decoy_exclude.fa"
        genome_fasta_fai = "/data/shared/genomes/hg38/GRCh38_masked_v2_decoy_exclude.fa.fai"
        genome_fasta_dict = "/data/shared/genomes/hg38/GRCh38_masked_v2_decoy_exclude.dict"
        genome_version="V3"
        cnvkit_germline_reference_PON="/data/shared/genomes/hg38/inhouse_DBs/hg38v3_primary/cnvkit/hg38v3_109samples.cnvkit.reference.cnn"
        cnvkit_inhouse_cnn_dir="/data/shared/genomes/hg38/inhouse_DBs/hg38v3_primary/cnvkit/wgs_persample_cnn/"
        inhouse_SV="/data/shared//genomes/hg38/inhouse_DBs/hg38v3_primary/"
        }

        // Gene and transcript annotation files:

        gencode_gtf = "/data/shared/genomes/hg38/gene.annotations/gencode.v36.annotation.gtf"
        gencode_gff3 = "/data/shared/genomes/hg38/gene.annotations/gencode.v36.annotation.gff3"
     
        //Program  files:
        msisensor_list="/data/shared/genomes/hg38/program_DBs/msisensor/hg38_msisensor_scan.txt"
        
        accucopy_config="/data/shared/genomes/hg38/accucopy/accucopy.docker.nextflow.conf"
        cnvradar_anno="/data/shared/genomes/hg38/program_DBs/cnvradar/All_20180418.vcf.gz"
        cnvradar_anno_idx="/data/shared/genomes/hg38/program_DBs/cnvradar/All_20180418.vcf.gz.tbi"
        cnvradar_ROI="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.SM.bed" 

        cnvradar_roisum_dir="/data/shared/genomes/hg38/program_DBs/cnvradar/inhouse_roi_summaries/"
       


        //Structural variants
        delly_exclude="/data/shared/genomes/hg38/program_DBs/delly/human.hg38.excl.tsv"
        
        smoove_exclude="/data/shared/genomes/hg38/interval.files/smoove/smoove.hg38.excluderegions.bed"
        smoove_gff="/data/shared/genomes/hg38/gene.annotations/GRCh38_latest_genomic.gff.gz"



        //Repeat Expansions:
        expansionhunter_catalog="/data/shared/genomes/hg38/program_DBs/expansionHunter/expansionHunter_hg38_stripy.variant_catalog.json"
        hipSTR_bed="/data/shared/genomes/hg38/interval.files/STRs/GRCh38.hipstr_reference.bed"

        // Somatic calling files:
        gatk_wgs_pon="/data/shared/genomes/hg38/program_DBs/GATK/somatic/somatic-hg38_1000g_pon.hg38.vcf.gz"
        mutect_gnomad="/data/shared/genomes/hg38/program_DBs/GATK/somatic/somatic-hg38_af-only-gnomad.hg38.vcf.gz"
        gatk_contamination_ref="/data/shared/genomes/hg38/program_DBs/GATK/somatic/somatic-hg38_small_exac_common_3.hg38.vcf.gz"


        pcgr_data_dir="/data/shared/genomes/hg38/program_DBs/PCGR/"


        // Program indexes:
        pcgr_assembly="grch38"
        sequenza_cg50_wig="/data/shared/genomes/hg38/program_DBs/sequenza/GRCh38.primary.cg50.sequenza.wig.gz"


        // Regions & variants:
        qualimap_ROI="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.20bp.SM.6col.bed"
        gencode_exons_ROI="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.SM.bed"

        ROI="/data/shared/genomes/hg38/interval.files/exome.ROIs/211130.hg38.refseq.gencode.fullexons.50bp.SM.bed"
        
        //ROI="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.20bp.SM.bed"

        callable_regions="/data/shared/genomes/hg38/interval.files/GATK.hg38.callable.regions.bed"
        manta_callable_regions="/data/shared/genomes/hg38/interval.files/manta/GATK.hg38.callable.regions.bed.gz"

        dbsnp="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
        KGindels="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz"
        KGindels_idx="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"

        KGmills="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
        KGmills_idx="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
        KG_p1_High_snps="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz"

        hapmap="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz"
        omni="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz"
        WES_ROI="/data/shared/genomes/hg38/interval.files/exome.ROIs/211130.hg38.refseq.gencode.fullexons.50bp.SM.bed"

        break;
}

switch (params.panel) {
    case "WES_2":
        ROI="${WES_ROI}";
        panelID="WES"
    break;

    case "WES":
        ROI="${WES_ROI}";
        panelID="WES_subpanel"
    break;

    default: 
        ROI="${WES_ROI}";
        panelID="WES"
    break;
}


if (!params.archiveStorage) {
outputDir="${params.outdir}/"
}

if (params.archiveStorage) {
outputDir="${tank_storage}/alignedData/${params.genome}/${params.outdir}/"
}


Channel
    .fromPath(params.intervals_list)
    .map { it -> tuple(it.baseName,it)}
    .set { HTC_interval_list }


// SYMLINK PROCESSES

process inputFiles_symlinks_fq{
    errorStrategy 'ignore'
    publishDir "${caseID}/${outputDir}/fq_symlinks/", mode: 'link', pattern:'*.{fastq,fq}.gz'
    input:
    tuple val(caseID), val(sampleID), path(r1),path(r2), val(type)// from read_input2
    
    output:
    tuple path(r1),path(r2)
    script:
    """
    """
}



process inputFiles_symlinks_cram{
    errorStrategy 'ignore'
    publishDir "${caseID}/${outputDir}/cram_TN_symlinks/", mode: 'link', pattern: '*.{ba,cr}*'
    publishDir "${caseID}/${outputDir}/variantcalls/Alignment_symlinks/", mode: 'link', pattern: "*.{ba,cr}*"

    input:

    tuple val(caseID), val(sampleID), path(cram), path(crai),val(type)    
    output:
    path("${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.symedit.cram")
    path("${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.symedit.cram.crai")
    script:
    """
    mv ${cram} ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.symedit.cram
    mv ${crai} ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.symedit.cram.crai
    """
}


///////////////////////// PREPROCESS MODULES ////////////////////////////////////

process tb_fastq_to_ubam {
    errorStrategy 'ignore'
    tag "$sampleID"
    cpus 20
    maxForks 10

    input:
    tuple val(caseID),val(sampleID), path(r1), path(r2), val(type)

    output:
    tuple val(caseID),val(sampleID), path("${sampleID}.unmapped.from.fq.bam"),val(type)
    //tuple path(r1),path(r2)
    
    script:
    """
    ${gatk_exec} FastqToSam \
    -F1 ${r1} \
    -F2 ${r2} \
    -SM ${sampleID} \
    -PL illumina \
    -PU KGA_PU \
    -RG KGA_RG \
    -O ${sampleID}.unmapped.from.fq.bam
    """
}

process tb_markAdapters {

    input:
    tuple val(caseID),val(sampleID), path(uBAM),val(type)
    
    output:
    tuple val(caseID),val(sampleID), path("${sampleID}.ubamXT.bam"), path("${sampleID}.markAdapterMetrics.txt"), val(type)
    
    script:

    """
    ${gatk_exec} MarkIlluminaAdapters \
    -I ${uBAM} \
    -O ${sampleID}.ubamXT.bam \
    --TMP_DIR ${tmpDIR} \
    -M ${sampleID}.markAdapterMetrics.txt
    """


}

process tb_align {
    tag "$sampleID"

    maxForks 5
    errorStrategy 'ignore'
    cpus 20

    input:
    tuple val(caseID),val(sampleID), path(uBAM), path(metrics),val(type)

    output:
    tuple val(caseID),val(sampleID), path("${sampleID}.${type}.${params.genome}.${genome_version}.QNsort.BWA.clean.bam"),val(type)
    
    script:
    """
    ${gatk_exec} SamToFastq \
    -I ${uBAM} \
    -INTER \
    -CLIP_ATTR XT \
    -CLIP_ACT 2 \
    -NON_PF \
    -F /dev/stdout \
    |  singularity run -B ${s_bind} ${simgpath}/bwa0717.sif bwa mem \
    -t ${task.cpus} \
    -p \
    ${genome_fasta} \
    /dev/stdin \
    | ${gatk_exec} MergeBamAlignment \
    -R ${genome_fasta} \
    -UNMAPPED ${uBAM} \
    -ALIGNED /dev/stdin \
    -MAX_GAPS -1 \
    -ORIENTATIONS FR \
    -SO queryname \
    -O ${sampleID}.${type}.${params.genome}.${genome_version}.QNsort.BWA.clean.bam
    """
}


process tb_markDup_v2_bam_cram {
    errorStrategy 'ignore'
    maxForks 6
    tag "$sampleID"
    publishDir "${caseID}/${outputDir}/CRAM/", mode: 'copy', pattern: "*.BWA.MD.cr*"
    
    input:
    tuple val(caseID),val(sampleID), path(bam), val(type)
    
    output:
    tuple val(caseID), val(sampleID), path("${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.bam"), path("${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD*bai"), val(type), emit: markDup_bam

    tuple val(caseID),val(sampleID),  path("${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.cram"), path("${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD*crai"),val(type), emit: markDup_cram


    script:
    """
    samtools view -h ${bam} \
    | samblaster | sambamba view -t 8 -S -f bam /dev/stdin | sambamba sort -t 8 --tmpdir=/data/TMP/TMP.${user}/ -o ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.bam /dev/stdin
    sambamba index ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.bam

    samtools view \
    -T ${genome_fasta} \
    -C \
    -o ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.cram ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.bam

    samtools index ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.cram

    """
}

process tb_markDup_v3_cram {
    errorStrategy 'ignore'
    maxForks 6
    tag "$sampleID"
    publishDir "${caseID}/${outputDir}/CRAM/", mode: 'copy', pattern: "*.BWA.MD.cr*"
    publishDir "${caseID}/${outputDir}/variantcalls/Alignment_symlinks/", mode: 'link', pattern: "*.BWA.MD.cr*"

    input:
    tuple val(caseID),val(sampleID), path(bam), val(type)
    
    output:
    tuple val(caseID),val(sampleID),  path("${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.cram"), path("${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD*crai"),val(type), emit: markDup_output
    
    script:
    """
    samtools view -h ${bam} \
    | samblaster | sambamba view -t 8 -S -f bam /dev/stdin | sambamba sort -t 8 --tmpdir=/data/TMP/TMP.${user}/ -o /dev/stdout /dev/stdin \
    |  samtools view \
    -T ${genome_fasta} \
    -C \
    -o ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.cram -

    samtools index ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.cram
    """
}

process tb_cram_bam {

   // publishDir "${caseID}/${params.outdir}/cram_TN_symlinks/", mode: 'link', pattern: '*.symedit*'
    input:
    tuple val(caseID), val(sampleID), path(cram), path(crai),val(type)

    output:
    tuple val(caseID), val(sampleID), path("*.bam"), path("*.bai"), val(type), emit:bam
    path("${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.symedit.cram")
    path("${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.symedit.cram.crai")
    script:
    """

    samtools view \
    -b \
    -T ${genome_fasta} \
    -o ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.bam ${cram}

    samtools index ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.bam
    mv ${cram} ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.symedit.cram
    mv ${crai} ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.symedit.cram.crai
    """
}


///////////////////////////////// QC MODULES ////////////////////////////
// TB QC input channel structure:
// tuple val(caseID), val(sampleID),  path(aln), path(aln_index),val(type)

process tb_samtools {
    errorStrategy 'ignore'
    tag "$sampleID"
    publishDir "${caseID}/${outputDir}/QC/", mode: 'copy'

    input:
    tuple val(caseID), val(sampleID),  path(bam), path(bai),val(type) 
    output:
    path("${sampleID}.samtools.sample.stats.txt")

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/samtools.sif samtools stats \
    ${bam} > ${sampleID}.samtools.sample.stats.txt
    """
}

process tb_qualimap {
    errorStrategy 'ignore'
    tag "$sampleID"
    publishDir "${caseID}/${outputDir}/QC/bamQC", mode: 'copy'

    cpus 10
    input:
    tuple val(caseID), val(sampleID),  path(bam), path(bai),val(type)
    //path(targetBED) from ch_qualimap_ROI
    output:
    path ("${sampleID}/")

    script:
    use_bed = qualimap_ROI ? "-gff ${qualimap_ROI}" : ''

    """
    unset DISPLAY
    singularity run -B ${s_bind} ${simgpath}/qualimap.sif \
    qualimap --java-mem-size=5G bamqc \
    -nt ${task.cpus} \
    -outdir ${sampleID} \
    -bam ${bam} $use_bed -sd -sdmode 0
    """
}

process tb_fastqc_bam {
    errorStrategy 'ignore'
    tag "$sampleID"
    cpus 1
    publishDir "${caseID}/${outputDir}/QC/", mode: 'copy',saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename" }
    input:
    tuple val(caseID), val(sampleID),  path(bam), path(bai),val(type)
    
    output:
    path "*_fastqc.{zip,html}"
    
    script:
    """
    singularity run -B ${s_bind} ${simgpath}/fastqc.sif --quiet --threads ${task.cpus} ${bam}
    """
}

process multiQC {
    errorStrategy 'ignore'
    publishDir "${launchDir}", mode: 'copy'
    //publishDir "${launchDir}/*/${params.outdir}", mode: 'copy'
    input:
    path("*_fastqc.*") 
    path("${sn}.${sampleID_type}.samtools.sample.stats.txt")
    path("bamQC/*")
    output:
    path ("*.multiQC.report.html")
    script:
    """
    singularity run -B ${s_bind} ${simgpath}/multiqc.sif \
    -c ${multiqc_config} \
    -n ${date}.TN_WES.multiQC.report.html \
    -f -q  ${launchDir}/*/${outputDir}/QC/
    """
}



/////////////////////////////// VARIANT CALLING MODULES ///////////////////////

process tb_haplotypecaller {
    errorStrategy 'ignore'
    cpus 4
    tag "$sampleID"
    publishDir "${caseID}/${outputDir}/variantcalls/", mode: 'copy', pattern: "*.haplotypecaller.*"
    publishDir "${caseID}/${outputDir}/variantcalls/gvcf/", mode: 'copy', pattern: "*.g.*"
    input:
    tuple val(caseID), val(sampleID), path(cram), path(crai),val(type) 
    
    output:
    tuple val(caseID), val(sampleID),  path("${caseID}.${sampleID}.${type}.haplotypecaller.vcf.gz"), path("${caseID}.${sampleID}.${type}.haplotypecaller.vcf.gz.tbi"),emit: sample_gvcf
    path("${caseID}.${sampleID}.${type}.g.*")
    path("${caseID}.${sampleID}.${type}.HCbamout.*")
    path("${crai}")
    path("${cram}")
    script:
    """
    ${gatk_exec} --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=30" HaplotypeCaller \
    -I ${cram} \
    -R ${genome_fasta} \
    -ERC GVCF \
    -L ${ROI} \
    --smith-waterman FASTEST_AVAILABLE \
    --native-pair-hmm-threads 30 \
    -pairHMM FASTEST_AVAILABLE \
    --dont-use-soft-clipped-bases \
    -O ${caseID}.${sampleID}.${type}.g.vcf.gz \
    -bamout ${caseID}.${sampleID}.${type}.HCbamout.bam
    
    ${gatk_exec} GenotypeGVCFs \
    -R ${genome_fasta} \
    -V ${caseID}.${sampleID}.${type}.g.vcf.gz \
    -O ${caseID}.${sampleID}.${type}.haplotypecaller.vcf.gz \
    -G StandardAnnotation \
    -G AS_StandardAnnotation
    """
}

process mutect2 {
    tag "$caseID"
    publishDir "${caseID}/${outputDir}/variantcalls/", mode: 'copy', pattern: "*.{VarSeq.*,PASSonly.*}"
    publishDir "${caseID}/${outputDir}/tumorBoard_files/", mode: 'copy', pattern: "*.for.VarSeq.{vcf,idx}"

    publishDir "${caseID}/${outputDir}/variantcalls/mutect2_bamout", mode: 'copy', pattern: "*.bamout.*"

    publishDir "${caseID}/${outputDir}/QC/mutect2_filtering/", mode: 'copy', pattern: "*.{table,stats,tsv}"

    input:
    tuple val(caseID), val(sampleID_normal), path(bamN), path(baiN),val(typeN), val(sampleID_tumor),path(bamT), path(baiT),val(typeT)
    output:
    tuple val(caseID), val(sampleID_tumor), path("${caseID}.mutect2.PASSonly.vcf.gz"), path("${caseID}.mutect2.PASSonly.vcf.gz.tbi"),emit: mutect2_allPASS 
  
    tuple val(caseID), path("${caseID}.mutect2.tumor.PASSonly.vcf.gz"), path("${caseID}.mutect2.tumor.PASSonly.vcf.gz.tbi"),emit: mutect2_tumorPASS 
    tuple val(caseID), path("${caseID}.mutect2.for.VarSeq.vcf.gz"), path("${caseID}.mutect2.for.VarSeq.vcf.gz.tbi"), emit: mutect2_vcf
    
    path("*.{tsv,table,stats}")
    path("*tumor.PASSonly.*")
    path("*.bamout.*")

    script:
    """
    ${gatk_exec} --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=30" Mutect2 \
    -R ${genome_fasta} \
    -I ${bamT} \
    -I ${bamN} \
    -normal ${sampleID_normal} \
    --germline-resource ${mutect_gnomad} \
    --panel-of-normals ${gatk_wgs_pon} \
    -L ${ROI} \
    --dont-use-soft-clipped-bases \
    --native-pair-hmm-threads 30 \
    -pairHMM FASTEST_AVAILABLE \
    --smith-waterman FASTEST_AVAILABLE \
    -O ${caseID}.mutect2.raw.vcf.gz \
    -bamout ${caseID}.mutect2.bamout.bam \
    --f1r2-tar-gz ${caseID}.within.f1r2.tar.gz
    
    ${gatk_exec} LearnReadOrientationModel \
    -I ${caseID}.within.f1r2.tar.gz \
    -O ${caseID}.within.ROmodel.tar.gz
    
    ${gatk_exec} GetPileupSummaries -I ${bamT} \
    -V ${gatk_contamination_ref} \
    -L ${gatk_contamination_ref} \
    -O ${caseID}.within.getpileupsummaries.table
    
    ${gatk_exec} CalculateContamination \
    -I ${caseID}.within.getpileupsummaries.table \
    -tumor-segmentation ${caseID}.segments.table \
    -O ${caseID}.contamination.table
    
    ${gatk_exec} FilterMutectCalls \
    -V ${caseID}.mutect2.raw.vcf.gz \
    -R ${genome_fasta} \
    --tumor-segmentation ${caseID}.segments.table \
    --contamination-table ${caseID}.contamination.table \
    --min-allele-fraction 0.001 \
    -O ${caseID}.mutect2.for.VarSeq.vcf.gz
    
    ${gatk_exec} SelectVariants -R ${genome_fasta} \
    -V ${caseID}.mutect2.for.VarSeq.vcf.gz \
    --exclude-filtered \
    -O ${caseID}.mutect2.PASSonly.vcf.gz
    
    ${gatk_exec} SelectVariants -R ${genome_fasta} \
    -V ${caseID}.mutect2.for.VarSeq.vcf.gz \
    --exclude-filtered -xl-sn ${sampleID_normal} --exclude-non-variants \
    -O ${caseID}.mutect2.tumor.PASSonly.vcf.gz
    """
}

process strelka2 {
    errorStrategy 'ignore'
    tag "$caseID"
    publishDir "${caseID}/${outputDir}/variantcalls/strelka2", mode: 'copy'
    cpus 10

    input: 
    tuple val(caseID), val(sampleID_normal), path(bamN), path(baiN),val(typeN), val(sampleID_tumor),path(bamT), path(baiT),val(typeT)

    output:
    path("*.strelka2.*")
    
    script:
    """
    singularity run -B ${s_bind} ${simgpath}/strelka2_2.9.10.sif /tools/strelka2/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam ${bamN} \
    --tumorBam ${bamT} \
    --referenceFasta  ${genome_fasta} \
    --exome \
    --runDir strelka

    singularity run -B ${s_bind} ${simgpath}/strelka2_2.9.10.sif python2 strelka/runWorkflow.py \
    -j ${task.cpus} \
    -m local

    mv strelka/results/variants/somatic.indels.vcf.gz ${caseID}.strelka2.somatic.indels.vcf.gz
    mv strelka/results/variants/somatic.indels.vcf.gz.tbi ${caseID}.strelka2.somatic.indels.vcf.gz.tbi
    mv strelka/results/variants/somatic.snvs.vcf.gz ${caseID}.strelka2.somatic.snvs.vcf.gz
    mv strelka/results/variants/somatic.snvs.vcf.gz.tbi ${caseID}.strelka2.somatic.snvs.vcf.gz.tbi

    ${gatk_exec} MergeVcfs \
    -I ${caseID}.strelka2.somatic.indels.vcf.gz \
    -I ${caseID}.strelka2.somatic.snvs.vcf.gz \
    -O ${caseID}.strelka2.merged.vcf 
    """
}

process msisensor {
    errorStrategy 'ignore'
    tag "$caseID"
    publishDir "${caseID}/${outputDir}/MSIsensor/", mode: 'copy'
    publishDir "${caseID}/${outputDir}/tumorBoard_files/", mode: 'copy', pattern: "*_msi"

    input: 
    tuple val(caseID), val(sampleID_normal), path(bamN), path(baiN),val(typeN), val(sampleID_tumor),path(bamT), path(baiT),val(typeT)

    output:
    path("*_msi*")
 
    script:
    """
    msisensor msi \
    -d ${msisensor_list} \
    -n ${bamN} -t ${bamT} \
    -e ${ROI} \
    -o ${caseID}_msi
    """
}
/*
process sequenza {
    errorStrategy 'ignore'
    tag "$caseID"
    
    input:
    tuple val(caseID), val(sampleID_normal), path(bamN), path(baiN),val(typeN), val(sampleID_tumor),path(bamT), path(baiT),val(typeT)
    output:
    tuple val(caseID), path("${caseID}.seqz.final.gz") 

    script:
    """
    sequenza-utils bam2seqz -n ${bamN} -t ${bamT} \
    --fasta ${genome_fasta} \
    -gc ${sequenza_cg50_wig} \
    -o ${caseID}.seqz.phase1.gz
    sequenza-utils seqz_binning --seqz ${caseID}.seqz.phase1.gz \
    -w 50 -o ${caseID}.seqz.final.gz 
    """
}

process sequenza_R_output {
    errorStrategy 'ignore'
    tag "$caseID"
    publishDir "${caseID}/${params.outdir}/", mode: 'copy'
    publishDir "${caseID}/${params.outdir}/tumorBoard_files/", mode: 'copy', pattern: "*_{segments,alternative_fit,genome_view}.{txt,pdf}"
    input:
    tuple val(caseID),  path(seqz)

    output:
    path("sequenza/*")


    script:
    """
    #!/usr/bin/env Rscript
    library(readr)
    readr::local_edition(1)
    library(sequenza)
    t1 = sequenza.extract("${seqz}",verbose=F)
    cp = sequenza.fit(t1)
    sequenza.results(sequenza.extract = t1, cp.table = cp, sample.id = "${caseID}", out.dir = "sequenza" )
    """
}
*/
process sequenza_conda {
    errorStrategy 'ignore'
    tag "$caseID"
    
    input:
    tuple val(caseID), val(sampleID_normal), path(bamN), path(baiN),val(typeN), val(sampleID_tumor),path(bamT), path(baiT),val(typeT)
    output:
    tuple val(caseID), path("${caseID}.seqz.final.gz") 
    
    script:
    """
    sequenza-utils bam2seqz \
    -n ${bamN} -t ${bamT} \
    --fasta ${genome_fasta} \
    -gc ${sequenza_cg50_wig} \
    -o ${caseID}.seqz.phase1.gz
    sequenza-utils seqz_binning --seqz ${caseID}.seqz.phase1.gz \
    -w 50 -o ${caseID}.seqz.final.gz 
    """
}

process sequenza_R_output_conda {
    errorStrategy 'ignore'
    tag "$caseID"
    publishDir "${caseID}/${outputDir}/", mode: 'copy'
    publishDir "${caseID}/${outputDir}/tumorBoard_files/", mode: 'copy', pattern: "*_{segments,alternative_fit,genome_view}.{txt,pdf}"

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sequenzaEnv'
    input:
    tuple val(caseID),  path(seqz)

    output:
    path("sequenza_conda/*")

    script:
    """
    #!/usr/bin/env Rscript
    library(sequenza)
    t1 = sequenza.extract("${seqz}",verbose=F)
    cp = sequenza.fit(t1)
    sequenza.results(sequenza.extract = t1, cp.table = cp, sample.id = "${caseID}", out.dir = "sequenza_conda" )
    """
}


process pcgr_v103 {
    tag "$caseID"
    errorStrategy 'ignore'
    publishDir "${caseID}/${outputDir}/PCGR/", mode: 'copy', pattern: "*.pcgr_acmg.*"
    publishDir "${caseID}/${outputDir}/tumorBoard_files/", mode: 'copy', pattern: "*.flexdb.html"
    input:
    tuple val(caseID),  path(vcf), path(idx), val(pcgr_tumor)

    output:
    path("*.pcgr_acmg.*")

    script:
    //tumorsite=${pcgr_tumor} ? "--tumor_site ${pcgr_tumor}" : ""
    """
    singularity run -B ${s_bind} ${simgpath}/pcgr103.sif pcgr \
    --input_vcf ${vcf} \
    --pcgr_dir ${pcgr_data_dir} --output_dir . \
    --genome_assembly ${pcgr_assembly} \
    --sample_id ${caseID}_TMB_all \
    --min_mutations_signatures 100 \
    --all_reference_signatures \
    --estimate_tmb --estimate_msi_status \
    --no_docker \
    --exclude_dbsnp_nonsomatic \
    --assay WES \
    --tumor_site ${pcgr_tumor} \
    --estimate_signatures \
    --include_trials


    singularity run -B ${s_bind} ${simgpath}/pcgr103.sif pcgr \
    --input_vcf ${vcf} \
    --pcgr_dir ${pcgr_data_dir} --output_dir . \
    --genome_assembly ${pcgr_assembly} \
    --sample_id ${caseID}_TMB_NonSyn \
    --min_mutations_signatures 100 \
    --all_reference_signatures \
    --estimate_tmb --estimate_msi_status \
    --no_docker \
    --exclude_dbsnp_nonsomatic \
    --assay WES \
    --tumor_site ${pcgr_tumor} \
    --estimate_signatures \
    --tmb_algorithm nonsyn \
    --include_trials
    """
}


process pcgr_v141 {
    tag "$caseID"
    errorStrategy 'ignore'
    publishDir "${caseID}/${outputDir}/PCGR141/", mode: 'copy', pattern: "*.pcgr_acmg.*"
    //publishDir "${caseID}/${params.outdir}/tumorBoard_files/", mode: 'copy', pattern: "*.flexdb.html"
    input:
    tuple val(caseID),  path(vcf), path(idx), val(pcgr_tumor)

    output:
    path("*.pcgr_acmg.*")
    
    script:
    //tumorsite=${pcgr_tumor} ? "--tumor_site ${pcgr_tumor}" : ""
    """
    singularity run -B ${s_bind} ${simgpath}/pcgr141.sif pcgr \
    --input_vcf ${vcf} \
    --pcgr_dir ${pcgr_data_dir} --output_dir . \
    --genome_assembly ${pcgr_assembly} \
    --sample_id ${caseID}_TMB_all \
    --min_mutations_signatures 100 \
    --all_reference_signatures \
    --estimate_tmb --estimate_msi_status \
    --exclude_dbsnp_nonsomatic \
    --assay WES \
    --tumor_site ${pcgr_tumor} \
    --estimate_signatures \
    --include_trials


    singularity run -B ${s_bind} ${simgpath}/pcgr141.sif pcgr \
    --input_vcf ${vcf} \
    --pcgr_dir ${pcgr_data_dir} --output_dir . \
    --genome_assembly ${pcgr_assembly} \
    --sample_id ${caseID}_TMB_NonSyn \
    --min_mutations_signatures 100 \
    --all_reference_signatures \
    --estimate_tmb --estimate_msi_status \
    --exclude_dbsnp_nonsomatic \
    --assay WES \
    --tumor_site ${pcgr_tumor} \
    --estimate_signatures \
    --tmb_algorithm nonsyn \
    --include_trials
    """
}



/////////// SUBWORKFLOWS

workflow SUB_DNA_PREPROCESS {

    take:
    case_fastq_input_ch     // caseid, sampleID, ,R1, R2, type
   
    main:
    inputFiles_symlinks_fq(case_fastq_input_ch)
    tb_fastq_to_ubam(case_fastq_input_ch)
    tb_markAdapters(tb_fastq_to_ubam.out)
    tb_align(tb_markAdapters.out)
    tb_markDup_v2_bam_cram(tb_align.out)

    emit:
    finalBam=tb_markDup_v2_bam_cram.out.markDup_bam //caseID, sampleID, BAM, BAI,type
    finalCram=tb_markDup_v2_bam_cram.out.markDup_cram //caseID, sampleID, CRAM, CRAI,type

}

workflow SUB_DNA_QC {
    take:
    case_sid_cram_crai
    main:
    tb_samtools(case_sid_cram_crai)
    tb_qualimap(case_sid_cram_crai)
    tb_fastqc_bam(case_sid_cram_crai)
    multiQC(tb_samtools.out.collect(),tb_qualimap.out.collect(),tb_fastqc_bam.out.collect())

}

workflow SUB_DNA_TUMOR_NORMAL {
    take:
    tumorNormal_bam_ch
    caseID_pcgrID
    main:
    mutect2(tumorNormal_bam_ch)
    strelka2(tumorNormal_bam_ch)
    msisensor(tumorNormal_bam_ch)
  //  sequenza(tumorNormal_bam_ch)
   // sequenza_R_output(sequenza.out)
    sequenza_conda(tumorNormal_bam_ch)
    sequenza_R_output_conda(sequenza_conda.out)
    pcgr_v103(mutect2.out.mutect2_tumorPASS.join(caseID_pcgrID))
    pcgr_v141(mutect2.out.mutect2_tumorPASS.join(caseID_pcgrID))
    emit:    
    mutect2_out=mutect2.out.mutect2_vcf

}
