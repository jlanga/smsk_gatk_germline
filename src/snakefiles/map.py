rule map_bwa_index:
    input:
        raw + "reference.fasta"
    output:
        mock = touch(maps + "reference"),
        aux = expand(
            maps + "reference.{extension}",
            extension = "amb ann bwt pac sa".split(" ")
        )
    params:
        prefix = maps + "reference"
    log: maps + "bwa_index.log"
    shell: "bwa index -p {params.prefix} {input} 2> {log} 1>&2"


rule map_link_alt:
    input: config["alternative"]
    output: maps + "reference.alt"
    shell: "ln --symbolic $(readlink --canonicalize {input}) {output}"


rule map_bwa:
    """

    Note: The -M flag causes BWA to mark shorter split hits as secondary
    (essential for Picard compatibility). https://gatkforums.broadinstitute.org/gatk/discussion/2799/howto-map-and-mark-duplicates
    """
    input:
        mock = maps + "reference",
        reference = raw + "reference.fasta",
        forward = raw + "{sample}/{library}_1.fq.gz",
        reverse = raw + "{sample}/{library}_2.fq.gz",
        ref_files = expand(
            maps + "reference.{extension}",
            extension = "amb ann bwt pac sa alt".split(" ")
        )
    output:
        bam = maps + "{sample}/{library}.raw.bam",
    params:
        group = "group1",
        sample = "{sample}",
        platform = "illumina",
        library = "{library}",
        unit = "unit1"
    log:
        maps + "{sample}/{library}.bwa.log"
    benchmark:
        maps + "{sample}/{library}.bwa.json"
    threads:
        max_threads
    shell:
        "(bwa mem "
            "-M "
            '-R \'@RG\\tID:{params.group}\\tSM:{params.sample}\\tPL:{params.platform}\\tLB:{params.library}\\tPU:{params.unit}\' '
            "-t {threads} "
            "{input.mock} "
            "{input.forward} "
            "{input.reverse} "
        "| picard SortSam "
             "INPUT=/dev/stdin "
             "OUTPUT={output.bam} "
             "SORT_ORDER=coordinate "
        ") 2> {log}"



rule map_markduplicates:
    input:
        bam=maps + "{sample}/{library}.raw.bam",
    output:
        bam = maps + "{sample}/{library}.dedup.bam",
        dupstats = maps + "{sample}/{library}_library.dupstats"
    params:
        memory="1g"
    log:
        maps + "{sample}/{library}.markduplicates.log"
    shell:
        "picard MarkDuplicates "
            "INPUT={input.bam} "
            "OUTPUT={output.bam} "
            "METRICS_FILE={output.dupstats} "
            "ASSUME_SORT_ORDER=coordinate "
            "VALIDATION_STRINGENCY=SILENT "
            "COMPRESSION_LEVEL=1 "
            "REMOVE_DUPLICATES=true "
            "QUIET=false "
        "2> {log}"



rule map_bqsr_pre:
    input:
        bam = maps + "{sample}/{library}.dedup.bam",
        reference=raw + "reference.fasta",
        vcf = "common_all_20170710.vcf.gz"  # TODO
    output:
        pre_table = maps + "{sample}/{library}.pre.table"
    shell:
        "gatk -T BaseRecalibrator "
            "--input_file {input.bam} "
            "-R {input.reference} "
            "--knownSites {input.vcf} "
            "--out {output.pre_table} "
        #"2> {log}"


rule map_bqsr_post:
    input:
        bam = maps + "{sample}/{library}.dedup.bam",
        reference = raw + "reference.fasta",
        vcfs="common_all_20170710.vcf.gz",  # TODO
        pre_table = maps + "{sample}/{library}.pre.table"
    output:
        post_table = maps + "{sample}/{library}.post.table"
    shell:
        "gatk -T BaseRecalibrator "
            "--input_file {input.bam} "
            "-R {input.reference} "
            "--knownSites {input.vcf} "
            "--BQSR {input.table} "
            "--out {output.post_table} "
        #"2> {log}"

rule map_bqsr_analyzecovariates:
    input:
        reference = raw + "reference.fasta",
        pre_table = maps + "{sample}/{library}.pre.table",
        post_table = maps + "{sample}/{library}.post.table"
    output:
        pdf = maps + "{sample}/{library}.bqsr.pdf"
    shell:
        "gatk -T AnalyzeCovariates "
            "-R {input.reference} "
            "-before {input.pre_table} "
            "-after {input.post_table} "
            "-plots {output.pdf} "
        "2> {log}"

rule map_bqsr_printreads:
    input:
        bam = maps + "{sample}/{library}.dedup.bam",
        bai = maps + "{sample}/{library}.dedup.bam.bai",
        reference = raw + "reference.fasta",
        pre_table = maps + "{sample}/{library}.pre.table",
    output:
        bam = maps + "{sample}/{library}.recalibrated.bam",
    shell:
        "gatk -T PrintReads "
            "-R {input.reference} "
            "--input_file {input.bam} "
            "-BQSR {input.pre_table} "
            "-o {output.bam}"



def get_library_files_from_sample(wildcards):
    """ TODO: needs improvement/simplification
    Return the list of libraries corresponding to a population and chromosome.
    """
    files = [
        maps + wildcards.sample + "/" + library + ".dedup.bam" \
        for library in config["samples"][wildcards.sample]
    ]
    return files


rule map_merge_libraries:
    input:
        bams = get_library_files_from_sample
    output:
        bam = maps + "{sample}.bam"
    shell:
        "samtools merge {output.bam} {input.bams}"
