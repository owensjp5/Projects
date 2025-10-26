params.linkfile="/Users/jackowens/Desktop/Projects/Nextflow/ref/links.txt"
params.fastqdir="/Users/jackowens/Desktop/Projects/Nextflow/spades/fastq"

params.read1="/Users/jackowens/Desktop/Projects/Nextflow/spades/fastq/ERR3335404_1.fastq.gz"
params.read2="/Users/jackowens/Desktop/Projects/Nextflow/spades/fastq/ERR3335404_2.fastq.gz"

params.SPADES_OUTPUT="/Users/jackowens/Desktop/Projects/Nextflow/spades/SPADES_OUTPUT"



process download {

publishDir("${params.fastqdir}", mode: 'copy')

input:
 path linkfile

output:
 path "*", emit: outputfile

script:
"""
cat $linkfile | xargs -P 2 -I{} wget '{}'
"""

}

process assemble {

publishDir("${params.SPADES_OUTPUT}", mode: 'copy')

input:
 path read1
 path read2

output:
 path "*", emit: spades_output

script:

"""
echo "${read1.simpleName}" | cut -d'_' -f1 | xargs -I{} spades.py --careful -1 "$read1" -2 "$read2" -o "{}"
"""

}


workflow {

link_ch=Channel.fromPath(params.linkfile)

download(link_ch)
download.out.outputfile.view()


read1_ch=Channel.fromPath(params.read1)
read2_ch=Channel.fromPath(params.read2)

assemble(read1_ch,read2_ch)
assemble.out.spades_output.view()

}
