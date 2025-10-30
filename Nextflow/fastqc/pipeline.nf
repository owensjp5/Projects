params.fastq="/Users/jackowens/Desktop/Projects/Nextflow/fastqc/fastq/*.fastq.gz"

process QualityControl {

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