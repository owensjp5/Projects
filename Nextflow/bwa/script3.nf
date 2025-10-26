params.ref='/Users/jackowens/Desktop/Projects/Nextflow/ref/sequence.fasta'
params.index_dir='/Users/jackowens/Desktop/Projects/Nextflow/bwa/index_dir'

process index {

publishDir("${params.index_dir}", mode: 'copy')

input:
    path genome

output:
    path "*"

script:
"""
bwa index $genome

"""

}

workflow {

    ref_ch=Channel.fromPath(params.ref)

    index(ref_ch)
    index.out.view()

}
