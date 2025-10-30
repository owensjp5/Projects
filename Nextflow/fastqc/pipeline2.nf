params.fastq="/Users/jackowens/Desktop/Projects/Nextflow/fastqc/fastq/*.fastq.gz"

params.qc_report="/Users/jackowens/Desktop/Projects/Nextflow/fastqc/fastqc_report"

process QualityControl {

publishDir("${params.qc_report}", mode:'copy')

input:
 path fastq

output:
 path "*"

script:
"""
fastqc $fastq
"""

}

workflow {

fastq_ch=Channel.fromPath(params.fastq)
QualityControl(fastq_ch)
QualityControl.out.view()

}