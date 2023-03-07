# CHR = list(range(1, 23))
CHR = 21

HAPS = 'exome_haps'

def fun_vcf_exome(wildcards):
    vcf = None
    if wildcards.haps == 'exome_haps':
        vcf = 's3://<s3-bucket-path>'
    return(vcf)

rule all:
    input:
        # VCF of haps
        expand('output/rfmix/input/{haps}/haps.chr{chr}.vcf.fixed.gz', haps = HAPS, chr = CHR)

#-------------------
# VCFs
#-------------------

rule cp_vcfgz:
    threads: 1000
    wildcard_constraints: haps = '(?!wgs_nodup).+' # not starts with wgs_nodup
    output: 'output/rfmix/input/{haps}/haps.chr{chr}.vcf.gz'
    params: lambda wildcards: fun_vcf_exome(wildcards)
    shell: 'echo {params} && aws s3 cp {params} {output}'


rule vcfgz_dups:
    output: 'output/rfmix/input/wgs_nodup_haps/dups.chr{chr}.txt'
    input: 'output/rfmix/input/wgs_haps/haps.chr{chr}.vcf.fixed.gz',
    shell:
        """
        zcat {input} | cut -f3 > {output}
        Rscript -e '
            library(tidyverse)
            v = read_lines("{output}")
            stopifnot(any(duplicated(v)))
            dups = v[duplicated(v)]
            write_lines(dups, "{output}")
        '
        """


rule cp_vcfgz_nodup:
    threads: 16
    wildcard_constraints: haps = 'wgs_nodup.+' # starts with wgs_nodup
    output: 'output/rfmix/input/{haps}/haps.chr{chr}.vcf.gz'
    input: vcfgz = 'output/rfmix/input/wgs_haps/haps.chr{chr}.vcf.gz',
        dups = rules.vcfgz_dups.output,
    shell: 'bcftools view --threads {threads} -e "ID=@{input.dups}" -Oz {input.vcfgz} > {output}'

#-------------------
# Index & fixed
#-------------------

rule index_vcfgz:
    threads: 16
    output: 'output/rfmix/input/{haps}/haps.chr{chr}.vcf.gz.tbi'
    input: 'output/rfmix/input/{haps}/haps.chr{chr}.vcf.gz'
    shell: 'bcftools index --threads {threads} -f -t {input}'

rule fixed_vcfgz:
    threads: 16
    output: 'output/rfmix/input/{haps}/haps.chr{chr}.vcf.fixed.gz'
    input: vcfgz = 'output/rfmix/input/{haps}/haps.chr{chr}.vcf.gz', 
        index = rules.index_vcfgz.output
    shell: 'bcftools view --threads {threads} -H -G {input.vcfgz} | gzip > {output}'
