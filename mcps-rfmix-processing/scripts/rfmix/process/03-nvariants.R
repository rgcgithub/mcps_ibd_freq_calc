library(tidyverse)
library(data.table)
library(glue)
library(unglue)
library(bedpca)

## args
no_snakemake = !exists("snakemake")
testing = no_snakemake

# testing
if(no_snakemake) {
  panel = 'rgc_1kg_hgdp'
  param = 'base2'
} else {
}

## extract #sites from RFMix log (if available)
if(panel == 'rgc_1kg_hgdp') {
  tab_log = lapply(1:22, function(x) {
    log = glue('rfmix/ref-1kg_hgdp-K3/run-{param}/chr-{x}/log.txt') %>% read_lines
    n = grep('common', log, value = TRUE) %>% 
      unglue_vec('{}...{x} SNPs') %>% as.integer
    tibble(chr = x, src = 'log', n = n)
  }) %>% bind_rows
} else {
  tab_log = tibble()
}

tab_vcf = lapply(1:22, function(x) {
  vcf = glue('output/rfmix/processed/panel-{panel}/param-{param}/probs-vcf-raw/probs.chr{x}.vcf.gz')
  print(vcf)
  n = vcftools_n_variants(vcf)
  tibble(chr = x, src = 'vcf', n = n)
}) %>% bind_rows

tab = tibble(chr = as.character(tab_log$chr), n_log = tab_log$n, n_vcf = tab_vcf$n) 
tab = bind_rows(tab, tibble(chr = 'all', n_log = sum(tab$n_log), n_vcf = sum(tab$n_vcf)))
