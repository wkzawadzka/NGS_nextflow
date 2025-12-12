#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.samples     = ['DRR656897','DRR656905']
params.out_dir     = 'results'
params.ref_dir     = 'data/reference'

include { fastqc as fastqc1 } from './modules/QC'
include { fastqc as fastqc2 } from './modules/QC'
include { multiqc_before } from './modules/QC'
include { multiqc_after } from './modules/QC'

// structure:
// data/
//     reference/
//         here fasta + vcf + indexes
//     DRR656897/
//         fastq1 fastq2
//     DRR656905/
//         fastq1 fastq2
//     bam/
// results/
//     qc/
//       before/    
//       after/
//     ...

// ---------------------------------------------------------
//  DOWNLOAD DATA
// ---------------------------------------------------------
process downloadReference {
    // note: Files are copied into the specified directory in an asynchronous manner, so they may not be immediately available in the publish directory at the end of the process execution.
    publishDir "${params.ref_dir}", mode: 'copy',  overwrite: false

    output:  
    path "Gallus_gallus_all.fa", emit: ref_fa
    path "gallus_gallus.vcf.gz",  emit: ref_vcf
    path "gallus_gallus.vcf.gz.csi", emit: ref_vcfi
    path "Gallus_gallus_all.fa.fai", emit: ref_fai
    path "Gallus_gallus_all.dict", emit: ref_dict

    script:
    """
    # if [ ! -f "${params.ref_dir}/Gallus_gallus_all.fa.gz" ]; then
    # download chromosomes
    for i in {1..39}; do
        wget -P . \
        "https://ftp.ensembl.org/pub/release-115/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.\${i}.fa.gz"
    done

    # concat into one genome
    zcat *.fa.gz > Gallus_gallus_all.fa
    # zcat *.fa.gz | gzip > Gallus_gallus_all.fa.gz

    # download VCF
    wget https://ftp.ensembl.org/pub/release-115/variation/vcf/gallus_gallus/gallus_gallus.vcf.gz
    wget https://ftp.ensembl.org/pub/release-115/variation/vcf/gallus_gallus/gallus_gallus.vcf.gz.csi

    # create index file
    samtools faidx Gallus_gallus_all.fa

    # create dict file .dict
    gatk CreateSequenceDictionary -R Gallus_gallus_all.fa
    """
}

process downloadSamples {
    tag "${sample_id}"

    input:
    val sample_id

    output:
    tuple val(sample_id), path("*_1.fastq.gz"), path("*_2.fastq.gz")

    publishDir "data/${sample_id}", mode: 'copy', overwrite: false

    script:
    """
    # if [ ! -f "data/${sample_id}/${sample_id}_1.fastq.gz" ]; then
    fastq-dump -X 1000000 --gzip --split-files ${sample_id} 
    """
}


// ---------------------------------------------------------
// INDEX REFERENCE
// ---------------------------------------------------------
process bwa_index {
    tag "BWA index"

    input:
    path ref_fa

    output:
    path "Gallus_gallus_all.*", emit: bwa_index

    publishDir "${params.ref_dir}", mode: 'copy', overwrite: false

    script:
    """
    # if [ ! -f "Gallus_gallus_all.fa.bwt" ]; then
    echo "[INFO] Indexing reference genome with BWA..."
    bwa index ${ref_fa}
    """
}


// ---------------------------------------------------------
//  TRIMMOMATIC
// ---------------------------------------------------------
process trimmomatic {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id),
          path("${sample_id}_trimmed_1P.fastq.gz"),
          path("${sample_id}_trimmed_2P.fastq.gz")

    script:
    """
    trimmomatic PE -threads 4 \
        ${fastq1} ${fastq2} \
        ${sample_id}_trimmed_1P.fastq.gz ${sample_id}_trimmed_1U.fastq.gz \
        ${sample_id}_trimmed_2P.fastq.gz ${sample_id}_trimmed_2U.fastq.gz \
        TRAILING:20 LEADING:20 SLIDINGWINDOW:4:20 MINLEN:60 
    """
}

// ---------------------------------------------------------
//  ALLIGNMENT
// ---------------------------------------------------------
process bwa_mem {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)
    path ref_fa
    path bwa_index // just so it is dependedent on the index being created

    output:
    tuple val(sample_id), path("${sample_id}_aligned.bam"), emit: bam

    publishDir "data/bam", mode: 'copy', overwrite: true

    script:
    """
    echo "[INFO] Running allignment with bwa mem..."
    bwa mem -R '@RG\\tID:HTVHTCCXY_8\\tPL:ILLUMINA\\tPU:HTVHTCCXY.8\\tSM:${sample_id}\\tLB:unknown' \
    ${ref_fa} ${fastq1} ${fastq2} | samtools view -b -o ${sample_id}_aligned.bam
    """
}

// ---------------------------------------------------------
//  POST-PROCESSING BAM
// ---------------------------------------------------------
process samtools_postprocess {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_aligned_final.bam"), path("${sample_id}_aligned_final.bam.bai"), emit: final_bam

    publishDir "${params.out_dir}/bam", mode: 'copy', overwrite: true

    script:
    """
    # 6.1 sort by query name
    samtools sort -n ${bam} -o ${sample_id}_temp_sort.bam

    # 6.2 fixmate (add CIGAR + mate info)
    samtools fixmate -mc ${sample_id}_temp_sort.bam ${sample_id}_temp_fixmate.bam

    # 6.3 sort again by coordinates
    samtools sort ${sample_id}_temp_fixmate.bam -o ${sample_id}_temp_coord.bam

    # 6.4 markdup and remove duplicates
    samtools markdup -r ${sample_id}_temp_coord.bam ${sample_id}_aligned_final.bam

    # remove intermediate temp files
    rm -f ${sample_id}_temp_*.bam

    # 6.5 index BAM
    samtools index ${sample_id}_aligned_final.bam

    """
}

// ---------------------------------------------------------
//  STATS
// ---------------------------------------------------------
process samtools_flagstat {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(bam), path(bam_bai)

    output:
    path "${sample_id}_flagstat.txt", emit: flagstat

    publishDir "${params.out_dir}/flagstat", mode: 'copy', overwrite: true

    script:
    """
    samtools flagstat ${bam} > ${sample_id}_flagstat.txt
    """
}

process samtools_coverage {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(bam), path(bam_bai)

    output:
    tuple path("${sample_id}_avg_depth.txt"),
          path("${sample_id}_coverage.txt"),
          emit: coverage_stats

    publishDir "${params.out_dir}/coverage", mode: 'copy', overwrite: true

    script:
    """
    echo "[INFO] Running depth1 calulation for ${sample_id}"
    # genome wide:
    # samtools depth -a ${bam} > ${sample_id}_depth.txt 
    # chromosome, position, depth
    # 3 = liczba readów pokrywających ta pozycje
    # liczba wierszy = liczba pozycji (count)
    # also where depth = 0!

    echo "[INFO] Running depth2 calulation for ${sample_id}"
    samtools depth -a ${bam}  | awk '{sum+=\$3; count++} END {print sum/count}' > ${sample_id}_avg_depth.txt

    # awk '{sum+=\$3; count++} END {print sum/count}' ${sample_id}_depth.txt > ${sample_id}_avg_depth.txt
    echo "[INFO] Running coverage calulation for ${sample_id}"
    samtools coverage ${bam} -o ${sample_id}_coverage.txt
    """
}

process qualimap_bamqc {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(bam), path(bam_bai)

    output:
    path "qualimap_${sample_id}", emit: qualimap_dir

    publishDir "${params.out_dir}/qualimap", mode: 'copy', overwrite: true

    script:
    """
    qualimap bamqc -bam ${bam} -outdir qualimap_${sample_id} --java-mem-size=4G
    """
}

process picard_markdup {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(bam), path(bam_bai)

    output:
    path "${sample_id}_marked.bam", emit: marked_bam
    path "${sample_id}_marked_metrics.txt", emit: dup_metrics

    publishDir "${params.out_dir}/picard", mode: 'copy', overwrite: true

    script:
    """
    picard MarkDuplicates I=${bam} O=${sample_id}_marked.bam M=${sample_id}_marked_metrics.txt REMOVE_DUPLICATES=false
    samtools index ${sample_id}_marked.bam
    """
}

// ---------------------------------------------------------
//  BQSR
// ---------------------------------------------------------
process base_recalibration {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(bam), path(bam_bai)
    path ref_fa
    path ref_vcf
    path ref_fai
    path ref_dict


    output:
    tuple val(sample_id), path("${sample_id}_recal.table"), emit: recal_table
    tuple val(sample_id), path("${sample_id}_bqsr.bam"), emit: bqsr_bam
    tuple val(sample_id), path("${sample_id}_recal_after.table"), emit: recal_table_after

    publishDir "${params.out_dir}/bqsr", mode: 'copy', overwrite: true

    script:
    """ 
    # index VCF
    gatk IndexFeatureFile -I ${ref_vcf}

    echo "[INFO] running BaseRecalibrator for ${sample_id} (before bqsr)"
    gatk BaseRecalibrator \
        -I ${bam} \
        -R ${ref_fa} \
        --known-sites ${ref_vcf} \
        -O ${sample_id}_recal.table

    echo "[INFO] running ApplyBQSR for ${sample_id}"
    gatk ApplyBQSR \
        -I ${bam} \
        -R ${ref_fa} \
        --bqsr-recal-file ${sample_id}_recal.table \
        -O ${sample_id}_bqsr.bam

    echo "[INFO] running BaseRecalibrator for ${sample_id} (after bqsr)"
    gatk BaseRecalibrator \
        -I ${sample_id}_bqsr.bam \
        -R ${ref_fa} \
        --known-sites ${ref_vcf} \
        -O ${sample_id}_recal_after.table

    """
}

process analyze_covariates {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(recal_table)
    tuple val(sample_id_), path(recal_table_after)

    output:
    path "${sample_id}_AnalyzeCovariates.pdf", emit: analyze_pdf

    publishDir "${params.out_dir}", mode: 'copy', overwrite: true

    script:
    """
    gatk AnalyzeCovariates \
        -before ${recal_table} \
        -after ${recal_table_after} \
        -plots ${sample_id}_AnalyzeCovariates.pdf
    """
}

// ---------------------------------------------------------
//  VARIANT CALLING
// ---------------------------------------------------------

process haplotype_caller {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(bam)
    path ref_fa
    path ref_fai
    path ref_dict

    output:
    path "gatk_${sample_id}.vcf.gz", emit: vcf

    publishDir "${params.out_dir}/variants", mode: 'copy', overwrite: true

    script:
    """
    # zwykły VCF daje tylko miejsca gdzie coś się zmieniło
    # GVCF daje także informacje o miejscach bez wariantu
    # czyli GVCF to taki “półprodukt”, który pozwala później połączyć dane z wielu ptaków/osobników
    gatk HaplotypeCaller \
        -R ${ref_fa} \
        -I ${bam} \
        -O gatk_${sample_id}.vcf.gz \
        -ERC GVCF
    """
}

process genomics_db_import {
    tag "GenomicsDBImport"

    input:
    path vcfs

    output:
    path "gatk_database", emit: gatk_database

    publishDir "${params.out_dir}/variants", mode: 'copy', overwrite: true

    script:
    """
    echo "[INFO] Creating chromosome list"
    seq 1 39 > chromosome.list

    echo "[INFO] Running GenomicsDBImport"

    for v in ${vcfs}; do
        if [ ! -f "\${v}.idx" ]; then
            gatk IndexFeatureFile -I \$v
        fi
    done

    gatk GenomicsDBImport \
        --genomicsdb-workspace-path gatk_database \
        -L chromosome.list \
        ${vcfs.collect{ "-V " + it }.join(' ')}
    """
}

process joint_genotyping {
    tag "GenotypeGVCFs"

    input:
    path gatk_database
    path ref_fa
    path ref_fai
    path ref_dict

    output:
    path "gatk_joint.vcf.gz", emit: joint_vcf

    publishDir "${params.out_dir}/variants", mode: 'copy', overwrite: true

    script:
    """
    gatk GenotypeGVCFs \
        -R ${ref_fa} \
        -V gendb://gatk_database \
        -O gatk_joint.vcf.gz
    """
}

process variant_call_bcftools {
    tag "bcftools joint calling"

    input:
    path(bams)        // list of BAM files
    path ref_fa

    output:
    path "bcftools_joint.vcf.gz", emit: bcftools_vcf

    publishDir "${params.out_dir}/variants", mode: 'copy', overwrite: true

    script:
    """
    echo "[INFO] Running joint variant calling with bcftools on ${bams.join(', ')}"
    
    bcftools mpileup -Ou -f ${ref_fa} ${bams.join(' ')} |
    bcftools call -vmO z -o bcftools_joint.vcf.gz
    """
}

// ---------------------------------------------------------
//  FILTERING VARIANTS
// ---------------------------------------------------------

process filter_gatk_variants {
    tag "Filter GATK variants"

    input:
    path gatk_vcf
    path ref_fa
    path ref_fai
    path ref_dict

    output:
    path "gatk_snps.vcf.gz", emit: gatk_snps
    path "gatk_variant_counts.txt", emit: variant_counts

    publishDir "${params.out_dir}/variants", mode: 'copy', overwrite: true

    script:
    """
    gatk IndexFeatureFile -I ${gatk_vcf}

    num_before=\$(bcftools view -H ${gatk_vcf} | wc -l)
    echo "Before filtration: $num_before" > gatk_variant_counts.txt

    # https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
    gatk VariantFiltration \
        -R ${ref_fa} \
        -V ${gatk_vcf} \
        -O gatk_filtered.vcf.gz \
        --filter-name "QualByDepth_le_2" \
        --filter-expression "QD < 2.0" \
        --filter-name "FisherStrandBias_over_60" \
        --filter-expression "FS > 60.0" \
        --filter-name "SOR_bias_gt_3" \
        --filter-expression "SOR > 3.0" \
        --filter-name "RootMeansSquared_MeanQuality_lt_40" \
        --filter-expression "MQ < 40.0" \

    # exclude filtered
    gatk SelectVariants \
        -V gatk_filtered.vcf.gz \
        -R ${ref_fa} \
        --exclude-filtered \
        -O gatk_filtered2.vcf.gz

    num_v=\$(bcftools view -H gatk_filtered2.vcf.gz | wc -l)
    echo "After GATK hard filtration: $num_v" >> gatk_variant_counts.txt

    # remove indels
    # bcftools view -v snps gatk_filtered2.vcf.gz -Oz -o gatk_snps.vcf.gz
    # lub można grep po ref/... czy więcej niż 1 character = nie snp

    gatk SelectVariants \
        -V gatk_filtered2.vcf.gz \
        -select-type SNP \
        -O gatk_snps.vcf.gz

    num_snps=\$(bcftools view -H gatk_snps.vcf.gz | wc -l)
    echo "[INFO] Number of SNPs for GATK: \$num_snps"
    echo "After SNP filtration: $num_snps" >> gatk_variant_counts.txt

    """
}

process filter_bcftools_variants {
    tag "Filter Bcftools variants"

    input:
    path bcftools_vcf
    path ref_fa

    output:
    path "bcftools_snps.vcf.gz", emit: bcftools_snps
    path "bcftools_variant_counts.txt", emit: variant_counts

    publishDir "${params.out_dir}/variants", mode: 'copy', overwrite: true

    script:
    """
    echo "[INFO] starting bcftools filtering on ${bcftools_vcf}"
    before=\$(bcftools view -H ${bcftools_vcf}  | wc -l)
    echo "Before any filtration: $before" > bcftools_variant_counts.txt

    # Remove 2 or more alleles
    bcftools view -i 'N_ALT=1' ${bcftools_vcf} -Oz -o step1.vcf.gz
    bcftools index step1.vcf.gz
    after_all=\$(bcftools view -H step1.vcf.gz  | wc -l)
    echo "After removal of 2 or more alleles: $after_all" >> bcftools_variant_counts.txt

    # Filter by QUAL
    bcftools view -i 'QUAL>=30' step1.vcf.gz -Oz -o step2.vcf.gz
    bcftools index step2.vcf.gz

    # saving qality scores, nice to do mean, meadian etc stats plots e.g. in R later on
    bcftools query -f '%QUAL\\n' step2.vcf.gz > bcftools_quals.txt
    after_qual=\$(bcftools view -H step2.vcf.gz  | wc -l)
    echo "After quality >=30 filtration: $after_qual" >> bcftools_variant_counts.txt

    # min & max coverage depth
    # the DP value in the INFO is the sum of the DP value over all samples in your vcf at this position.
    # The DP value in the FORMAT column in the read depth for the given sample at this position. 
    # If you just have one sample this value should be equal.
    bcftools view -i 'INFO/DP>=2 && INFO/DP<=10' step2.vcf.gz -Oz -o step3.vcf.gz
    bcftools index step3.vcf.gz
    after_dp=\$(bcftools view -H step2.vcf.gz  | wc -l)
    echo "After INFO/DP filtration: $after_dp" >> bcftools_variant_counts.txt

    # remove variants with the same basepair position (if two variants have the same bp position both are remove)
    # - Resolves issue with SNP and INDEL calls at same position
    bcftools query -f '%CHROM\\t%POS\\n' step3.vcf.gz | sort | uniq -c > counts.txt
    awk '\$1>1 {print \$2"\\t"(\$3-1)"\\t"\$3}' counts.txt > dup_positions.bed
    rm counts.txt

    # remove duplicated positions if needed
    if [ -s dup_positions.bed ]; then
        bcftools view -T ^dup_positions.bed step3.vcf.gz -Oz -o bcftools_filtered.vcf.gz
    else
        cp step3.vcf.gz bcftools_filtered.vcf.gz
    fi
    bcftools index bcftools_filtered.vcf.gz

    after_dup=\$(bcftools view -H bcftools_filtered.vcf.gz  | wc -l)
    echo "After removal of variants with same basepair position: $after_dup" >> bcftools_variant_counts.txt

    # remove indels
    bcftools view -v snps bcftools_filtered.vcf.gz -Oz -o bcftools_snps.vcf.gz
    bcftools index bcftools_snps.vcf.gz

    num_snps=\$(bcftools view -H bcftools_snps.vcf.gz | wc -l)
    echo "[INFO] After removal of indels: number of SNPs for BCFtools: \$num_snps"
    echo "After SNP filtration: $num_snps" >> bcftools_variant_counts.txt
    """
}



// ------------------- WORKFLOW -------------------
workflow {
    samples_id_ch = Channel.from(params.samples)

    // download & index reference
    ref = downloadReference()
    bwa_index_ch = bwa_index(ref.ref_fa)

    // download samples
    samples = downloadSamples(samples_id_ch)

    // QC before trimming
    fastqc_before = fastqc1(samples, "before")
    multiqc_before_report = multiqc_before(fastqc_before.qc_reports.collect())

    // trimming
    samples_trimmed = trimmomatic(samples)

    // FastQC after trimming
    fastqc_after = fastqc2(samples_trimmed, "after")

    // allignment
    aligned_bams = bwa_mem(samples_trimmed, ref.ref_fa, bwa_index_ch.bwa_index)

    // postprocess BAMs
    final_bams = samtools_postprocess(aligned_bams.bam)

    // stats, QC after
    flagstat = samtools_flagstat(final_bams)
    coverage = samtools_coverage(final_bams)
    qualimaps = qualimap_bamqc(final_bams)
    picards   = picard_markdup(final_bams)
    multiqc_after(fastqc_after.qc_reports.collect(), flagstat.flagstat.collect(), coverage.coverage_stats.collect(), qualimaps.qualimap_dir.collect(), picards.dup_metrics.collect())

    // BQSR
    bqsr_tables = base_recalibration(final_bams, ref.ref_fa, ref.ref_vcf, ref.ref_fai, ref.ref_dict)
    analyze_covariates(bqsr_tables.recal_table, bqsr_tables.recal_table_after)

    // ---------------------------------------
    // VARIANT CALLING -> GATK
    // ---------------------------------------
    haplotype_caller_bqsr = haplotype_caller(bqsr_tables.bqsr_bam, ref.ref_fa, ref.ref_fai, ref.ref_dict)
    genomics_db = genomics_db_import(haplotype_caller_bqsr.vcf.collect())
    gatk_vcf = joint_genotyping(genomics_db.gatk_database, ref.ref_fa, ref.ref_fai, ref.ref_dict)

    // ---------------------------------------
    // VARIANT CALLING -> Bcftools
    // ---------------------------------------
    bams_list = bqsr_tables.bqsr_bam.map { sample_id, bam -> bam }.collect()
    bcftools_vcf = variant_call_bcftools(bams_list, ref.ref_fa)

    // ---------------------------------------
    // FILTER
    // ---------------------------------------
    gatk_snps = filter_gatk_variants(gatk_vcf.joint_vcf, ref.ref_fa, ref.ref_fai, ref.ref_dict)
    bcftools_snps = filter_bcftools_variants(bcftools_vcf.bcftools_vcf, ref.ref_fa)
}
