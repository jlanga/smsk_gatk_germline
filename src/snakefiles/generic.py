rule faidx_fasta:
    input: "{filename}.fasta"
    output: "{filename}.fasta.fai"
    shell: "samtools faidx {input}"

rule dict_fasta:
    input: "{filename}.fasta"
    output: "{filename}.dict"
    shell: "samtools dict {input} > {output}"

rule bai:
    input: "{filename}.bam"
    output: "{filename}.bam.bai"
    shell: "samtools index {input}"
