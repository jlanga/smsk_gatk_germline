shell.prefix("set -euo pipefail;")
configfile: "src/config.yaml"

max_threads = 96

samples = [x for x in config["samples"]]
libraries = {sample: library for sample in samples for library in config["samples"][sample]}

snakefiles = "src/snakefiles/"

include: snakefiles + "folders.py"
include: snakefiles + "clean.py"
include: snakefiles + "generic.py"
include: snakefiles + "raw.py"
include: snakefiles + "map.py"
include: snakefiles + "variants.py"

rule all:
    input:
        ## raw
        # raw + "reference.fasta",
        # raw + "reference.dict",
        # maps + "reference.fasta.alt",
        # ## map
        # expand(
        #     maps + "{sample}/lib1.raw.bam",
        #     sample=samples,
        # ),
        # expand(
        #     maps + "{sample}/lib1.dedup.bam",
        #     sample = samples,
        # ),
        # expand(
        #     maps + "{sample}/lib1.recalibrated.bam",
        #     sample = samples,
        # ),
        expand(
            maps + "{sample}.bam",
            sample = samples
        ),
        ## Variants
        # expand(
        #     variants + "{sample}.haplotypecaller.g.vcf",
        #     sample = samples
        # ),
        # expand(
        #     variants + "{sample}.snp.recal",
        #     sample = samples
        # ),
        # expand(
        #     variants + "{sample}.snp.tranches",
        #     sample = samples
        # ),
        # expand(
        #     variants + "{sample}.snps_recalibrated.g.vcf",
        #     sample = samples
        # ),
        # expand(
        #     variants + "{sample}.indel.tranches",
        #     sample = samples
        # ),
        # expand(
        #     variants + "{sample}.indels_recalibrated.g.vcf",
        #     sample = samples
        # ),
