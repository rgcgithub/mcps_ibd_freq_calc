# Processing of RFMix v2 output (rfmix re-run fixref)
# 1. convert RFMix output from tsv too vcf format
#    - ref. panel/run: mais_k3
#    - rfmix param: rbase (base + reanalyze)
#    - processing flags: raw, trimmed
# 2.
#    - Copy from S3, index, fixed
#    - Check concordance with BIM
# 2. Prepare input VCF of Ref. samples

NBATCHES = 40

# CHR = list(range(1, 23)) + ['X']
# CHR = list(range(1, 23))
CHR = 15

PANEL = '1kg_hgdp'; PARAM = 'rbase'

# VCF = 'raw'
VCF = 'clean'

def fun_rfmix_fb_s3(wildcards):
    fb = None
    if wildcards.panel == '1kg_hgdp':
        if wildcards.param == 'rbase':
            fb = 's3:<s3-bucket-path'
    return(fb)

def fun_rfmix_k(wildcards):
    k = None
    if wildcards.panel == '1kg_hgdp':
        k = 3
    return(k)

rule all:
    input:
        # Probs. VCFs
        expand('output/rfmix/processed/panel-{panel}/param-{param}/probs-vcf-{vcf}/probs.chr{chr}.vcf.fixed.gz', panel = PANEL, param = PARAM, vcf = VCF, chr = CHR),
        # Local anc. (lanc)
        # expand('output/rfmix/processed/panel-{panel}/param-{param}/lanc-vcf-{vcf}/lanc.unrel.chr{chr}.tsv.gz', panel = PANEL, param = PARAM, vcf = VCF, chr = CHR),
        # expand('output/rfmix/processed/panel-{panel}/param-{param}/ganc.tsv.gz', panel = PANEL, param = PARAM, vcf = VCF),

#--------------------
# Global anc
#--------------------

rule ganc_avr:
    output: 'output/rfmix/processed/ganc.avr.tsv.gz'
    input: unrel = 'output/MC150K/rel2/v2-king-ibdseg/degree-3rd/ivs.unrel.samples.txt',
    script: 'scripts/rfmix/process/02-ganc.R'

rule ganc_samples:
    output: 'output/rfmix/processed/panel-{panel}/param-{param}/ganc.tsv.gz'
    script: 'scripts/rfmix/process/02b-ganc-samples.R'

#--------------------
# Local anc. (lanc)
#--------------------

rule lanc_unrel:
    threads: 1000
    output: 'output/rfmix/processed/panel-{panel}/param-{param}/lanc-vcf-{vcf}/lanc.unrel.chr{chr}.tsv.gz'
    input: vcfgz = 'output/rfmix/processed/panel-{panel}/param-{param}/probs-vcf-{vcf}/probs.chr{chr}.vcf.gz',
        unrel = "output/MC150K/rel2/v2-king-ibdseg/degree-3rd/ivs.unrel.samples.txt",
        index = 'output/rfmix/processed/panel-{panel}/param-{param}/probs-vcf-{vcf}/probs.chr{chr}.vcf.gz.tbi',
        # fixed = 'output/rfmix/processed/panel-{panel}/param-{param}/probs-vcf-{vcf}/probs.chr{chr}.vcf.fixed.gz',
    params: k = lambda wildcards: fun_rfmix_k(wildcards)
    script: 'scripts/rfmix/process/01-derive-lanc.R'

#---------------------------
# Probs VCFs by batches
#---------------------------

## 1. Copy *.fb.tsv, extract meta/nvariants

rule vcf_cp_fb:
    output: temp('output/rfmix/processed/panel-{panel}/param-{param}/probs-vcf-{vcf}/tmp.chr{chr}.fb.tsv')
    params: fb_s3 = lambda wildcards: fun_rfmix_fb_s3(wildcards),
        fbgz_clean = 'output/rfmix/processed/panel-{panel}/param-{param}/fb-tsv-clean/fb.chr{chr}.tsv.gz',
    run:
        if wildcards.vcf == 'raw':
            shell('echo raw && aws s3 cp {params.fb_s3} {output}')
        else:
            shell('echo clean && zcat {params.fbgz_clean} > {output}')

rule vcf_meta_fb:
    output: 'output/rfmix/processed/panel-{panel}/param-{param}/probs-vcf-{vcf}/fb.meta.chr{chr}.tsv'
    input: rules.vcf_cp_fb.output
    shell: "cut {input} -f1-4 | tail -n+2 > {output}"

rule vcf_nvariants_fb:
    output: 'output/rfmix/processed/panel-{panel}/param-{param}/probs-vcf-{vcf}/fb.nvariants.chr{chr}.txt'
    input: rules.vcf_meta_fb.output
    shell:
        """
            Rscript -e '
                library(tidyverse)
                n = read_tsv("{input}") %>% nrow
                write_lines(n, "{output}")
            '
        """

## 2. Compress & zndex *.fb.tsv

rule vcf_gz_fb:
    output: temp('output/rfmix/processed/panel-{panel}/param-{param}/probs-vcf-{vcf}/tmp.chr{chr}.fb.tsv.gz')
    input: rules.vcf_cp_fb.output
    shell: "gzip -c {input} > {output}"

rule vcf_zindex_fbgz:
    output: temp('output/rfmix/processed/panel-{panel}/param-{param}/probs-vcf-{vcf}/tmp.chr{chr}.fb.tsv.gz.zindex')
    input: rules.vcf_gz_fb.output
    shell: "zindex {input} --regex '^.{{1}}' -v"

# 3. Split *.fb.tsv by batches

rule vcf_splits:
    input: rules.vcf_nvariants_fb.output
    output: temp('output/rfmix/processed/panel-{panel}/param-{param}/probs-vcf-{vcf}/batches/chr-{chr}/splits.tsv.gz')
    shell:
        """
        Rscript -e '
            library(tidyverse); library(glue); library(bedpca)
            n = read_lines("{input}") %>% as.numeric
            stopifnot(!is.na(n)); stopifnot(length(n) == 1)
            cat(" - n", n, "\n")
            sp = bedpca:::split_seq(n, nb = {NBATCHES}, correct = TRUE)
            cat(" - #batches", nrow(sp), "\n")
            write_tsv(sp, "{output}")
        '
        """

rule vcf_nbatches:
    output: temp('output/rfmix/processed/panel-{panel}/param-{param}/probs-vcf-{vcf}/batches/chr-{chr}/nbatches.txt')
    input: rules.vcf_splits.output
    shell:
        """
        Rscript -e '
            library(tidyverse); library(glue); library(bedpca)
            sp = fread("{input}")
            write_lines(nrow(sp), "{output}")
            '
        """

rule vcf_fb_batch:
    output: temp('output/rfmix/processed/panel-{panel}/param-{param}/probs-vcf-{vcf}/batches/chr-{chr}/batch{b}.fb.tsv')
    input: splits = rules.vcf_splits.output,
        fb = rules.vcf_gz_fb.output, fb_zindex = rules.vcf_zindex_fbgz.output,
    shell: 
        """
            Rscript -e '
                library(tidyverse); library(glue)
                b = as.numeric("{wildcards.b}")
                sp = read_tsv("{input.splits}")
                stopifnot(b <= nrow(sp))
                cols = seq(sp$beg[b], sp$end[b])
                cols = c(1, 2, 2 + cols) # 2 header lines
                cols_str = paste(cols, collapse = " ")
                cmd = glue("zq {input.fb} --line {{cols_str}} > {output}")
                cat(" - cmd:", cmd, "\n")
                ret = system(cmd)
                stopifnot(ret == 0)
            '
        """

# 4. Process every *.fb.tsv batch, i.e. convert tsv to vcf

rule vcf_covert_batch:
    threads: 4 # for more RAM
    output: temp('output/rfmix/processed/panel-{panel}/param-{param}/probs-vcf-{vcf}/batches/chr-{chr}/batch{b}.vcf')
    input: fb = rules.vcf_fb_batch.output, splits = rules.vcf_splits.output,
    shell: 
        """
            Rscript -e '
                library(tidyverse); library(bedpca)
                ret = rfmix_convert_vcf("{input.fb}", "{output}", bs = 50, verbose = 2)
                sp = read_tsv("{input.splits}")

                b = as.numeric("{wildcards.b}")
                batch_size  = sp$size[b]
                n_variants = vcftools_n_variants("{output}")
                stopifnot(batch_size == n_variants)
            '
        """

rule vcf_gz_batch:
    output: temp('output/rfmix/processed/panel-{panel}/param-{param}/probs-vcf-{vcf}/batches/chr-{chr}/batch{b}.vcf.gz')
    input: rules.vcf_covert_batch.output
    shell: "bcftools view -Oz {input} > {output}"

## 5. Merge batches

rule vcf_merge_batches:
    threads: 16
    output: 'output/rfmix/processed/panel-{panel}/param-{param}/probs-vcf-{vcf}/probs.chr{chr}.vcf.gz'
    input: expand('output/rfmix/processed/panel-{{panel}}/param-{{param}}/probs-vcf-{{vcf}}/batches/chr-{{chr}}/batch{b}.vcf.gz', b = list(range(1, NBATCHES + 1)))
    shell: "bcftools concat --threads {threads} -Oz {input} > {output}"

rule vcf_index:
    threads: 16
    output: 'output/rfmix/processed/panel-{panel}/param-{param}/probs-vcf-{vcf}/probs.chr{chr}.vcf.gz.tbi'
    input: rules.vcf_merge_batches.output
    shell: 'bcftools index --threads {threads} -f -t {input}'

rule vcf_fixed:
    output: 'output/rfmix/processed/panel-{panel}/param-{param}/probs-vcf-{vcf}/probs.chr{chr}.vcf.fixed.gz'
    input: vcfgz = rules.vcf_merge_batches.output, index = rules.vcf_index.output
    shell: "Rscript -e 'library(bedpca); vcftools_fixed(\"{input.vcfgz}\", write_file = TRUE)'"

