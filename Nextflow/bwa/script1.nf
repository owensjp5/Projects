params.ref='/Users/jackowens/Desktop/Projects/Nextflow/ref/sequence.fasta'

process index {

input:
    path genome

script:
"""
bwa index $genome

"""

}

workflow {

    ref_ch=Channel.fromPath(params.ref)

    index(ref_ch)

}