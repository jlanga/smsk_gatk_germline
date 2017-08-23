rule variants_haplotypecaller:
    input:
        bam = maps + "{sample}.bam",
        reference = raw + "reference.fasta",
        fai = raw + "reference.fasta.fai",
        dictionary = raw + "reference.dict",
        bai = maps + "{sample}.bam.bai"
    output:
        gvcf = variants + "{sample}.haplotypecaller.g.vcf",
        idx =  variants + "{sample}.haplotypecaller.g.vcf.idx",
    params:
        stand_call_conf = 30,
        genotyping_mode = "DISCOVERY",
    log:
        variants + "{sample}.haplotypecaller.log"
    shell:
        "gatk -T HaplotypeCaller "
            "-R {input.reference} "
            "-I {input.bam} "
            "--genotyping_mode {params.genotyping_mode} "
            "-stand_call_conf {params.stand_call_conf} "
            "--emitRefConfidence GVCF "
            "-o {output.vcf} "
        "2> {log} 1>&2"


rule variants_snp_variantrecalibrator:
    input:
        reference = raw + "reference.fasta",
        gvcf = variants + "{sample}.haplotypecaller.g.vcf",
        hapmap = "data/known_variants/hapmap_3.3.hg38.vcf.gz",
        omni = "data/known_variants/1000G_omni2.5.hg38.vcf.gz",
        thougen = "data/known_variants/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        dbsnp = "data/known_variants/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
    output:
        recal = variants + "{sample}.snp.recal",
        tranches = variants + "{sample}.snp.tranches",
        rplot = variants + "{sample}.snp.R"
    params:
        tranches = None
    log:
        variants + "{sample}.snp.variantrecalibrator.log"
    shell:
        "gatk -T VariantRecalibrator "
            "-R {input.reference} "
            "-input {input.gvcf} "
            "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} "
            "-resource:omni,known=false,training=true,truth=true,prior=12.0 {input.omni} "
            "-resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.thougen} "
            "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} "
            "-an DP "
            "-an QD "
            "-an FS "
            "-an SOR "
            "-an MQ "
            "-an MQRankSum "
            "-an ReadPosRankSum "
            "-an InbreedingCoeff "
            "-mode SNP "
            "-tranche 100.0 "
            "-tranche 99.9 "
            "-tranche 99.0 "
            "-tranche 90.0 "
            "-recalFile {output.recal} "
            "-tranchesFile {output.tranches} "
            "-rscriptFile {output.rplot} "
        "2> {log}"




rule variants_snp_applyrecalibration:
    input:
        reference = raw + "reference.fasta",
        gvcf = variants + "{sample}.haplotypecaller.g.vcf",
        recal = variants + "{sample}.snp.recal",
        tranches = variants + "{sample}.snp.tranches",
    output:
        gvcf = variants + "{sample}.snps_recalibrated.g.vcf"
    params:
        ts_filter_level = 99.0
    log:
        variants + "{sample}.snp.applyrecalibration.log"
    shell:
        "gatk -T ApplyRecalibration "
            "-R {input.reference} "
            "-input {input.gvcf} "
            "-mode SNP "
            "--ts_filter_level {params.ts_filter_level} "
            "-recalFile {input.recal} "
            "-tranchesFile {input.tranches} "
            "-o {output.gvcf} "


rule variants_indel_variantrecalibrator:
    input:
        reference = raw + "reference.fasta",
        gvcf = variants + "{sample}.snps_recalibrated.g.vcf",
        mills = "data/known_variants/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        dbsnp_indels = "data/known_variants/Homo_sapiens_assembly38.known_indels.vcf.gz",
    output:
        recal = variants + "{sample}.indel.recal",
        tranches = variants + "{sample}.indel.tranches",
        rplot = variants + "{sample}.indel.R"
    params:
        max_gaussians = 4
    log:
        variants + "{sample}.snp.variantrecalibrator.log"
    shell:
        "gatk -T VariantRecalibrator "
            "-R {input.reference} "
            "-input {input.gvcf} "
            "-resource:mills,known=false,training=true,truth=true,prior=12.0 {input.mills} "
            "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp_indels} "
            "-an DP "
            "-an QD "
            "-an FS "
            "-an SOR "
            "-an MQ "
            "-an MQRankSum "
            "-an ReadPosRankSum "
            "-an InbreedingCoeff "
            "-mode INDEL "
            "-tranche 100.0 "
            "-tranche 99.9 "
            "-tranche 99.0 "
            "-tranche 90.0 "
            "--maxGaussians {params.max_gaussians} "
            "-recalFile {output.recal} "
            "-tranchesFile {output.tranches} "
            "-rscriptFile {output.rplot} "
        "2> {log}"


rule variants_indel_applyrecalibration:
    input:
        reference = raw + "reference.fasta",
        gvcf = variants + "{sample}.snps_recalibrated.g.vcf",
        recal = variants + "{sample}.indel.recal",
        tranches = variants + "{sample}.indel.tranches",
    output:
        gvcf = variants + "{sample}.indels_recalibrated.g.vcf"
    params:
        ts_filter_level = 99.0
    log:
        variants + "{sample}.indel.applyrecalibration.log"
    shell:
        "gatk -T ApplyRecalibration "
            "-R {input.reference} "
            "-input {input.gvcf} "
            "-mode INDEL "
            "--ts_filter_level {params.ts_filter_level} "
            "-recalFile {input.recal} "
            "-tranchesFile {input.tranches} "
            "-o {output.gvcf} "
