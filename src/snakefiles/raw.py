rule raw_link_genome:
    input: config["reference"]
    output: raw + "reference.fasta"
    shell: "ln --symbolic $(readlink --canonicalize {input}) {output}"


rule raw_link_pe_sample:
    """
    Make a link to the original file, with a prettier name than default.
    """
    input:
        forward= lambda wildcards: config["samples"][wildcards.sample][wildcards.library]["forward"],
        reverse= lambda wildcards: config["samples"][wildcards.sample][wildcards.library]["reverse"]
    output:
        forward= raw + "{sample}/{library}_1.fq.gz",
        reverse= raw + "{sample}/{library}_2.fq.gz"
    shell:
        "ln --symbolic $(readlink --canonicalize {input.forward}) {output.forward}; "
        "ln --symbolic $(readlink --canonicalize {input.reverse}) {output.reverse}"
