library(tidyverse)
library(data.table)
library(glue)

library(bedpca)

## par.
vals_chr = c(1:22, 'X')
# vals_chr = 21:22

## data
if(!exists('len')) {
  len = lapply(vals_chr, function(chr) {
    cat(chr)
    lapply(c('exome', 'wgs'), function(src) {
      f_vcf = switch(src,
        'exome' = glue('output/rfmix/input/exome_haps/haps.chr{chr}.vcf.gz'),
        'wgs' = glue('output/rfmix/input/wgs_haps/haps.chr{chr}.vcf.gz'),
        stop())
      size = switch(src,
        'exome' = 500, 'wgs' = 500, stop())
      if(file.exists(f_vcf)) {
        n = vcftools_n_variants(f_vcf)
        sp = split_seq(n, size = size)
        tibble(chr = chr, n = n, size = size, nb = nrow(sp), src = src)
      } else {
        tibble()
      }
    }) %>% bind_rows
  }) %>% bind_rows
}

