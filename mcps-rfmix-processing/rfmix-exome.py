import tempfile

# CHR = list(range(1, 23))
# CHR = list(range(2, 23))
CHR = 15
# CHR = 'X'
# batch size = 500
DICT_BATCHES_EXOME = {"1": 2412, "2": 1784, "3": 1423, "4": 992, "5": 1089, "6": 1197,
        "7": 1179, "8": 891, "9": 1059, "10": 1015, "11": 1420, "12": 1293,
        "13": 453, "14": 754, "15": 834, "16": 1164, "17": 1399, "18": 398,
        "19": 1595, "20": 627, "21": 270, "22": 564,
        "X": 1200}
# batch size = 500
DICT_BATCHES_WGS = {"1": 20649, "2": 22543, "3": 18833, "4": 18239, "5": 16953, "6": 16079,
        "7": 15106, "8": 14592, "9": 11330, "10": 12705, "11": 12788, "12": 12457,
        "13": 9164, "14": 8331, "15": 7572, "16": 8444, "17": 7492, "18": 7250,
        "19": 5806, "20": 5933, "21": 3316, "22": 3468,
        "X": 16000}

PANEL = '1kg_hgdp'; PARAM = 'rbase'
VCF = 'clean_wgs_nodup'

def fun_tmpdir(wildcards):
    suffix = '-' + wildcards.chr + '-' + wildcards.vcf
    with tempfile.TemporaryDirectory(suffix = suffix) as t:
        tempdir = t
    tempdir_sub = tempdir.split('/')[-1]
    tempdir_local = 'tmp/' + tempdir_sub
    return(tempdir_local)

def fun_nbatches(wildcards, vcf = None, seq = False):
    # query vcf
    if 'vcf' in wildcards.keys():
        qvcf = wildcards.vcf 
    elif vcf is not None:
        qvcf = vcf
    else:
        print('vcf is underfined')
        exit(1)

    if qvcf in ('raw', 'raw_rerun', 'clean', 'clean_rerun'):
        nb = DICT_BATCHES_EXOME.get(str(wildcards.chr))
    elif qvcf in ('raw_wgs_nodup', 'clean_wgs_nodup'):
        nb = DICT_BATCHES_WGS.get(str(wildcards.chr))
    else:
        print('vcf is unknown')
        exit(1)
    seqb = list(range(1, nb + 1))
    if(seq):
        return(seqb)
    else:
        return(nb)

def fun_vcf_haps(wildcards, ext = '.gz'):
    vcf = None
    # disable WGS for a while (indexing is going on)
    if wildcards.vcf in ('raw', 'clean'):
        vcf = 'output/rfmix/input/exome_haps/haps.chr' + wildcards.chr + '.vcf'
    elif wildcards.vcf in ('raw_rerun', 'clean_rerun'):
        vcf = 'output/rfmix/input/exome_rerun_haps/haps.chr' + wildcards.chr + '.vcf'
    elif wildcards.vcf in ('raw_wgs_nodup', 'clean_wgs_nodup'):
        vcf = 'output/rfmix/input/wgs_nodup_haps/haps.chr' + wildcards.chr + '.vcf'
    vcf = vcf + ext

    return(vcf)

# array probs (raw or clean)
def fun_vcf_probs(wildcards, ext = '.gz'):
    if wildcards.vcf in ('raw', 'raw_rerun', 'raw_wgs_nodup'):
        vcf = 'raw'
    elif wildcards.vcf in ('clean', 'clean_rerun', 'clean_wgs_nodup'):
        vcf = 'clean'
    else:
        print('vcf is unknown')
        exit(1)
    vcf = 'output/rfmix/processed/panel-' + wildcards.panel + '/param-' + wildcards.param + '/probs-vcf-' + vcf + '/probs.chr' + wildcards.chr + '.vcf'
    vcf = vcf + ext

    return(vcf)

def fun_rfmix_k(wildcards):
    k = 3
    return(k)

rule all:
    input:
        # VCF
        # expand('output/rfmix/exome/panel-{panel}/param-{param}/probs-{vcf}/probs.chr{chr}.vcf.gz', panel = PANEL, param = PARAM, vcf = VCF, chr = CHR),
        expand('output/rfmix/exome/panel-{panel}/param-{param}/probs-{vcf}/probs.chr{chr}.vcf.fixed.gz', panel = PANEL, param = PARAM, vcf = VCF, chr = CHR)
        # sorted VCF
        # expand('output/rfmix/exome/panel-{panel}/param-{param}/probs-{vcf}/probs.sorted.chr{chr}.vcf.fixed.gz', panel = PANEL, param = PARAM, vcf = VCF, chr = CHR),

#-------------------------------------
# 1. Copy large files from S3 
# - done in rfmix-input.py
#-------------------------------------

#--------------------------------------------------
# Split into batches, run every batch & merge
#--------------------------------------------------

# extract variants from VCF & check for duplicates
rule variants:
    output: 'output/rfmix/exome/panel-{panel}/param-{param}/probs-{vcf}/chr-{chr}/haps_exome.chr{chr}.variants.txt'
    input: lambda wildcards: fun_vcf_haps(wildcards, '.fixed.gz')
    shell: 
        """
        echo 'input: {input}'
        zcat {input} | cut -f3 > {output}
        Rscript -e '
            library(tidyverse)
            v = read_lines("{output}")
            stopifnot(!any(duplicated(v)))
        '
        """

rule batch_splits:
    output: 'output/rfmix/exome/panel-{panel}/param-{param}/probs-{vcf}/chr-{chr}/splits.tsv.gz'
    input: rules.variants.output
    params: nbatches = lambda wildcards: fun_nbatches(wildcards),
    shell:
        """
        Rscript -e '
            library(tidyverse); library(glue); library(bedpca)
            n = read_lines("{input}") %>% length 
            stopifnot(!is.na(n)); stopifnot(length(n) == 1)
            cat(" - n", n, "\n")
            nb = {params.nbatches} %>% as.numeric
            stopifnot(!is.na(nb)); stopifnot(length(nb) == 1)
            # sp = bedpca:::split_seq(n, nb = nb)
            beg1 = seq(1, n, length = nb + 1) %>% as.integer
            stopifnot(beg1[1] == 1); stopifnot(tail(beg1, 1) == n); stopifnot(length(beg1) == nb + 1)
            beg = head(beg1, -1) # remove the last element
            end1 = tail(beg, -1) # remove the first element
            end = c(end1 - 1, n)
            sp = tibble(beg = beg, end = end) %>% mutate(size = end - beg + 1)
            stopifnot(nrow(sp) == nb); stopifnot(all(sp$size > 0)); stopifnot(sum(sp$size) == n)
            cat(" - #batches", nrow(sp), "\n")
            write_tsv(sp, "{output}")
        '
        """

rule batch_vcf:
    # threads: 2 # for more RAM
    output: temp('output/rfmix/exome/panel-{panel}/param-{param}/probs-{vcf}/chr-{chr}/{b}.batch.vcf'),
    input: variants = rules.variants.output,
        splits = rules.batch_splits.output,
        haps = lambda wildcards: fun_vcf_haps(wildcards, '.gz'),
        haps_index = lambda wildcards: fun_vcf_haps(wildcards, '.gz.tbi'),
        haps_fixed = lambda wildcards: fun_vcf_haps(wildcards, '.fixed.gz'),
        probs = lambda wildcards: fun_vcf_probs(wildcards, '.gz'),
        probs_index = lambda wildcards: fun_vcf_probs(wildcards, '.gz.tbi'),
        probs_fixed = lambda wildcards: fun_vcf_probs(wildcards, '.fixed.gz'),
    params: k = lambda wildcards: fun_rfmix_k(wildcards)
    script: 'scripts/rfmix/exome/02c-interpolate-lai-vcf.R'

rule batch_vcf_gz:
    threads: 4
    output: 'output/rfmix/exome/panel-{panel}/param-{param}/probs-{vcf}/chr-{chr}/{b}.batch.vcf.gz'
    input: rules.batch_vcf.output
    shell: 'bcftools view --threads {threads} -Oz {input} > {output}'

rule batch_vcf_index:
    threads: 4
    output: 'output/rfmix/exome/panel-{panel}/param-{param}/probs-{vcf}/chr-{chr}/{b}.batch.vcf.gz.tbi'
    input: rules.batch_vcf_gz.output
    shell: 'bcftools index --threads {threads} -f -t {input}'

## Merge VCF batches
rule batch_merge_raw:
    threads: 1000
    output: vcfgz = 'output/rfmix/exome/panel-{panel}/param-{param}/probs-raw/probs.chr{chr}.vcf.gz',
        index = 'output/rfmix/exome/panel-{panel}/param-{param}/probs-raw/probs.chr{chr}.vcf.gz.tbi',
    input: vcfgz = lambda wildcards: expand('output/rfmix/exome/panel-' + wildcards.panel + '/param-' + wildcards.param + '/probs-raw/chr-' + wildcards.chr + '/{b}.batch.vcf.gz', b = fun_nbatches(wildcards, vcf = 'raw', seq = True)),
        index = lambda wildcards: expand('output/rfmix/exome/panel-' + wildcards.panel + '/param-' + wildcards.param + '/probs-raw/chr-' + wildcards.chr + '/{b}.batch.vcf.gz.tbi', b = fun_nbatches(wildcards, vcf = 'raw', seq = True)),
    params: 
        nbatches = lambda wildcards: fun_nbatches(wildcards, vcf = 'raw'),
        vcf_dir = 'output/rfmix/exome/panel-{panel}/param-{param}/probs-raw/chr-{chr}/',
        flist = 'output/rfmix/exome/panel-{panel}/param-{param}/probs-raw/chr-{chr}/vcf.list',
    shell: 
        """
            echo ' - write {params.nbatches} batches / vcfs to a file'
            Rscript -e '
                library(tidyverse); library(glue)
                files = sapply(seq({params.nbatches}), function(x) glue("{params.vcf_dir}/{{x}}.batch.vcf.gz"))
                write_lines(files, "{params.flist}")
            '
            echo ' - concatinating batches / vcfs by bcftools'
            bcftools concat --threads {threads} -f {params.flist} -Oz > {output.vcfgz}
            rm {params.flist}
            echo ' - indexing merged vcf'
            bcftools index --threads {threads} -f -t {output.vcfgz}
        """

rule batch_merge_clean:
    threads: 16
    output: 'output/rfmix/exome/panel-{panel}/param-{param}/probs-clean/probs.chr{chr}.vcf.gz'
    input: lambda wildcards: expand('output/rfmix/exome/panel-' + wildcards.panel + '/param-' + wildcards.param + '/probs-clean/chr-' + wildcards.chr + '/{b}.batch.vcf.gz.tbi', b = fun_nbatches(wildcards, vcf = 'clean', seq = True)),
    params: 
        nbatches = lambda wildcards: fun_nbatches(wildcards, vcf = 'clean'),
        vcf_dir = 'output/rfmix/exome/panel-{panel}/param-{param}/probs-clean/chr-{chr}/',
        flist = 'output/rfmix/exome/panel-{panel}/param-{param}/probs-clean/chr-{chr}/vcf.list',
    shell: 
        """
            echo ' - write {params.nbatches} batches / vcfs to a file'
            Rscript -e '
                library(tidyverse); library(glue)
                files = sapply(seq({params.nbatches}), function(x) glue("{params.vcf_dir}/{{x}}.batch.vcf.gz"))
                write_lines(files, "{params.flist}")
            '
            echo ' - merging batches / vcfs by vcftools'
            bcftools concat --threads {threads} -f {params.flist} -Oz > {output}
            rm {params.flist}
        """
    # shell: "bcftools concat --threads {threads} -Oz {input} > {output}"

rule batch_merge_raw_wgs_nodup:
    threads: 1000
    output: vcfgz = 'output/rfmix/exome/panel-{panel}/param-{param}/probs-raw_wgs_nodup/probs.chr{chr}.vcf.gz',
        index = 'output/rfmix/exome/panel-{panel}/param-{param}/probs-raw_wgs_nodup/probs.chr{chr}.vcf.gz.tbi'
    input: lambda wildcards: expand('output/rfmix/exome/panel-' + wildcards.panel + '/param-' + wildcards.param + '/probs-raw_wgs_nodup/chr-' + wildcards.chr + '/{b}.batch.vcf.gz.tbi', b = fun_nbatches(wildcards, vcf = 'raw_wgs_nodup', seq = True)),
    params: 
        nbatches = lambda wildcards: fun_nbatches(wildcards, vcf = 'raw_wgs_nodup'),
        vcf_dir = 'output/rfmix/exome/panel-{panel}/param-{param}/probs-raw_wgs_nodup/chr-{chr}/',
        flist = 'output/rfmix/exome/panel-{panel}/param-{param}/probs-raw_wgs_nodup/chr-{chr}/vcf.list',
    shell: 
        """
            echo ' - write {params.nbatches} batches / vcfs to a file'
            Rscript -e '
                library(tidyverse); library(glue)
                files = sapply(seq({params.nbatches}), function(x) glue("{params.vcf_dir}/{{x}}.batch.vcf.gz"))
                write_lines(files, "{params.flist}")
            '
            echo ' - merging batches / vcfs by vcftools'
            bcftools concat --threads {threads} -f {params.flist} -Oz > {output.vcfgz}
            rm {params.flist}
            echo ' - indexing merged vcf'
            bcftools index --threads {threads} -f -t {output.vcfgz}
        """


rule batch_merge_clean_wgs_nodup:
    threads: 16
    output: 'output/rfmix/exome/panel-{panel}/param-{param}/probs-clean_wgs_nodup/probs.chr{chr}.vcf.gz'
    input: lambda wildcards: expand('output/rfmix/exome/panel-' + wildcards.panel + '/param-' + wildcards.param + '/probs-clean_wgs_nodup/chr-' + wildcards.chr + '/{b}.batch.vcf.gz.tbi', b = fun_nbatches(wildcards, vcf = 'clean_wgs_nodup', seq = True)),
    params: 
        nbatches = lambda wildcards: fun_nbatches(wildcards, vcf = 'clean_wgs_nodup'),
        vcf_dir = 'output/rfmix/exome/panel-{panel}/param-{param}/probs-clean_wgs_nodup/chr-{chr}/',
        flist = 'output/rfmix/exome/panel-{panel}/param-{param}/probs-clean_wgs_nodup/chr-{chr}/vcf.list',
    shell: 
        """
            echo ' - write {params.nbatches} batches / vcfs to a file'
            Rscript -e '
                library(tidyverse); library(glue)
                files = sapply(seq({params.nbatches}), function(x) glue("{params.vcf_dir}/{{x}}.batch.vcf.gz"))
                write_lines(files, "{params.flist}")
            '
            echo ' - merging batches / vcfs by vcftools'
            bcftools concat --threads {threads} -f {params.flist} -Oz > {output}
            rm {params.flist}
        """

rule batch_merge_raw_rerun:
    threads: 1000
    output: vcfgz = 'output/rfmix/exome/panel-{panel}/param-{param}/probs-raw_rerun/probs.chr{chr}.vcf.gz',
        index = 'output/rfmix/exome/panel-{panel}/param-{param}/probs-raw_rerun/probs.chr{chr}.vcf.gz.tbi',
    input: vcfgz = lambda wildcards: expand('output/rfmix/exome/panel-' + wildcards.panel + '/param-' + wildcards.param + '/probs-raw_rerun/chr-' + wildcards.chr + '/{b}.batch.vcf.gz', b = fun_nbatches(wildcards, vcf = 'raw_rerun', seq = True)),
        index = lambda wildcards: expand('output/rfmix/exome/panel-' + wildcards.panel + '/param-' + wildcards.param + '/probs-raw_rerun/chr-' + wildcards.chr + '/{b}.batch.vcf.gz.tbi', b = fun_nbatches(wildcards, vcf = 'raw_rerun', seq = True)),
    params: 
        nbatches = lambda wildcards: fun_nbatches(wildcards, vcf = 'raw_rerun'),
        vcf_dir = 'output/rfmix/exome/panel-{panel}/param-{param}/probs-raw_rerun/chr-{chr}/',
        flist = 'output/rfmix/exome/panel-{panel}/param-{param}/probs-raw_rerun/chr-{chr}/vcf.list',
    shell: 
        """
            echo ' - write {params.nbatches} batches / vcfs to a file'
            Rscript -e '
                library(tidyverse); library(glue)
                files = sapply(seq({params.nbatches}), function(x) glue("{params.vcf_dir}/{{x}}.batch.vcf.gz"))
                write_lines(files, "{params.flist}")
            '
            echo ' - concatinating batches / vcfs by bcftools'
            bcftools concat --threads {threads} -f {params.flist} -Oz > {output.vcfgz}
            rm {params.flist}
            echo ' - indexing merged vcf'
            bcftools index --threads {threads} -f -t {output.vcfgz}
        """

rule batch_merge_clean_rerun:
    threads: 1000
    output: vcfgz = 'output/rfmix/exome/panel-{panel}/param-{param}/probs-clean_rerun/probs.chr{chr}.vcf.gz',
        index = 'output/rfmix/exome/panel-{panel}/param-{param}/probs-clean_rerun/probs.chr{chr}.vcf.gz.tbi',
    input: vcfgz = lambda wildcards: expand('output/rfmix/exome/panel-' + wildcards.panel + '/param-' + wildcards.param + '/probs-clean_rerun/chr-' + wildcards.chr + '/{b}.batch.vcf.gz', b = fun_nbatches(wildcards, vcf = 'clean_rerun', seq = True)),
        index = lambda wildcards: expand('output/rfmix/exome/panel-' + wildcards.panel + '/param-' + wildcards.param + '/probs-clean_rerun/chr-' + wildcards.chr + '/{b}.batch.vcf.gz.tbi', b = fun_nbatches(wildcards, vcf = 'clean_rerun', seq = True)),
    params: 
        nbatches = lambda wildcards: fun_nbatches(wildcards, vcf = 'clean_rerun'),
        vcf_dir = 'output/rfmix/exome/panel-{panel}/param-{param}/probs-clean_rerun/chr-{chr}/',
        flist = 'output/rfmix/exome/panel-{panel}/param-{param}/probs-clean_rerun/chr-{chr}/vcf.list',
    shell: 
        """
            echo ' - write {params.nbatches} batches / vcfs to a file'
            Rscript -e '
                library(tidyverse); library(glue)
                files = sapply(seq({params.nbatches}), function(x) glue("{params.vcf_dir}/{{x}}.batch.vcf.gz"))
                write_lines(files, "{params.flist}")
            '
            echo ' - concatinating batches / vcfs by bcftools'
            bcftools concat --threads {threads} -f {params.flist} -Oz > {output.vcfgz}
            rm {params.flist}
            echo ' - indexing merged vcf'
            bcftools index --threads {threads} -f -t {output.vcfgz}
        """

#---------------------------
# Post-process final VCFs
#---------------------------

rule merged_vcf_fixed_raw:
    threads: 4
    output: 'output/rfmix/exome/panel-{panel}/param-{param}/probs-raw/probs.chr{chr}.vcf.fixed.gz'
    input: vcfgz = rules.batch_merge_raw.output.vcfgz
    shell: 'bcftools view --threads {threads} -H -G {input.vcfgz} | gzip > {output}'

rule merged_vcf_index_clean:
    threads: 16
    output: 'output/rfmix/exome/panel-{panel}/param-{param}/probs-clean/probs.chr{chr}.vcf.gz.tbi'
    input: rules.batch_merge_clean.output
    shell: 'bcftools index --threads {threads} -f -t {input}'

rule merged_vcf_fixed_clean:
    threads: 16
    output: 'output/rfmix/exome/panel-{panel}/param-{param}/probs-clean/probs.chr{chr}.vcf.fixed.gz'
    input: vcfgz = rules.batch_merge_clean.output, index = rules.merged_vcf_index_clean.output
    shell: 'bcftools view --threads {threads} -H -G {input.vcfgz} | gzip > {output}'

rule merged_vcf_fixed_raw_wgs_nodup:
    threads: 16
    output: 'output/rfmix/exome/panel-{panel}/param-{param}/probs-raw_wgs_nodup/probs.chr{chr}.vcf.fixed.gz'
    input: vcfgz = rules.batch_merge_raw_wgs_nodup.output.vcfgz
    shell: 'bcftools view --threads {threads} -H -G {input.vcfgz} | gzip > {output}'

rule merged_vcf_index_clean_wgs_nodup:
    threads: 16
    output: 'output/rfmix/exome/panel-{panel}/param-{param}/probs-clean_wgs_nodup/probs.chr{chr}.vcf.gz.tbi'
    input: rules.batch_merge_clean_wgs_nodup.output
    shell: 'bcftools index --threads {threads} -f -t {input}'

rule merged_vcf_fixed_clean_wgs_nodup:
    threads: 16
    output: 'output/rfmix/exome/panel-{panel}/param-{param}/probs-clean_wgs_nodup/probs.chr{chr}.vcf.fixed.gz'
    input: vcfgz = rules.batch_merge_clean_wgs_nodup.output, index = rules.merged_vcf_index_clean_wgs_nodup.output
    shell: 'bcftools view --threads {threads} -H -G {input.vcfgz} | gzip > {output}'

rule merged_vcf_fixed_raw_rerun:
    threads: 4
    output: 'output/rfmix/exome/panel-{panel}/param-{param}/probs-raw_rerun/probs.chr{chr}.vcf.fixed.gz'
    input: vcfgz = rules.batch_merge_raw_rerun.output.vcfgz
    shell: 'bcftools view --threads {threads} -H -G {input.vcfgz} | gzip > {output}'

rule merged_vcf_fixed_clean_rerun:
    threads: 4
    output: 'output/rfmix/exome/panel-{panel}/param-{param}/probs-clean_rerun/probs.chr{chr}.vcf.fixed.gz'
    input: vcfgz = rules.batch_merge_clean_rerun.output.vcfgz
    shell: 'bcftools view --threads {threads} -H -G {input.vcfgz} | gzip > {output}'

#---------------------------
# Sorted VCF after merging
#---------------------------

rule merged_sorted_vcf:
    output: 'output/rfmix/exome/panel-{panel}/param-{param}/probs-{vcf}/probs.sorted.chr{chr}.vcf.gz'
    params: vcfgz = 'output/rfmix/exome/panel-{panel}/param-{param}/probs-{vcf}/probs.chr{chr}.vcf.gz',
        tmpdir = lambda wildcards: fun_tmpdir(wildcards),
        mem = ' 27G ',
    shell: "mkdir -p {params.tmpdir} && bcftools sort --max-mem {params.mem} -T {params.tmpdir} {params.vcfgz} -Oz -o {output} && rm -rf {params.tmpdir}"
    # shell: 'mkdir -p {params.tmpdir} && echo {params.tmpdir}'

rule merged_sorted_vcf_index:
    threads: 16
    output: 'output/rfmix/exome/panel-{panel}/param-{param}/probs-{vcf}/probs.sorted.chr{chr}.vcf.gz.tbi'
    input: rules.merged_sorted_vcf.output
    shell: 'bcftools index --threads {threads} -f -t {input}'

rule merged_sorted_vcf_fixed:
    threads: 16
    output: 'output/rfmix/exome/panel-{panel}/param-{param}/probs-{vcf}/probs.sorted.chr{chr}.vcf.fixed.gz'
    input: vcfgz = rules.merged_sorted_vcf.output, index = rules.merged_sorted_vcf_index.output
    shell: 'bcftools view --threads {threads} -H -G {input.vcfgz} | gzip > {output}'

