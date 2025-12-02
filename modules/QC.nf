// ---------------------------------------------------------
//  QUALITY CONTROL
// ---------------------------------------------------------
process fastqc {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)
    val before_or_after

    output:
    path "*_fastqc.zip", emit: qc_reports

    publishDir "${params.out_dir}/qc/${before_or_after}", mode: 'copy', overwrite: true

    script:
    """
    echo "[INFO] Running FastQC on sample ${sample_id} (${before_or_after} trimming)"
    fastqc -o . ${fastq1} ${fastq2} --noextract
    """
}

// ---------------------------------------------------------
// MULTIQC BEFORE TRIMMING
// ---------------------------------------------------------
process multiqc_before {
    tag "MultiQC before"

    input:
    path fastqc_zips

    output:
    path "multiqc_before.html"

    publishDir "${params.out_dir}", mode: 'copy', overwrite: true

    script:
    """
    echo "[INFO] FastQC files: ${fastqc_zips}"
    multiqc ${fastqc_zips} -o . -n multiqc_before.html
    """
}

// ---------------------------------------------------------
// MULTIQC AFTER TRIMMING 
// ---------------------------------------------------------
process multiqc_after {
    tag "MultiQC after"

    input:
    path fastqc_zips
    path flagstat  
    path coverage  
    path qualimap   
    path dups      

    output:
    path "multiqc_after.html"

    publishDir "${params.out_dir}", mode: 'copy', overwrite: true

    script:
    """
    echo "[INFO] FastQC files: ${fastqc_zips}"
    echo "[INFO] Flagstat files: ${flagstat}"
    echo "[INFO] Coverage files: ${coverage}"
    echo "[INFO] Qualimap directories: ${qualimap}"
    echo "[INFO] Picard duplicates files: ${dups}"

    multiqc ${fastqc_zips} \
        ${flagstat} \
        ${coverage} \
        ${qualimap} \
        ${dups} \
        -o . -n multiqc_after.html
    """
}