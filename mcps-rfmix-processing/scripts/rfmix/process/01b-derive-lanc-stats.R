library(tidyverse)
library(data.table)
library(glue)
library(rex)

library(bedpca)

library(parallel)
cores <- detectCores() - 1

## par
verbose = 1
bsize = 1

## args
no_snakemake = !exists("snakemake")
testing = no_snakemake

# testing
if(no_snakemake) {
  f_probs = "output/rfmix/processed/panel-oxford_mais_k3/param-rbase/probs-vcf-raw/probs.chr21.vcf.gz"
  f_unrel = "output/MC150K/rel2/v2-king-ibdseg/degree-3rd/ivs.unrel.samples.txt"
  f_out = "tmp.tsv.gz"
} else {
  f_probs = snakemake@input[[1]]
  f_unrel = snakemake@input[[2]]
  k_str = snakemake@params[[1]]
  f_out = snakemake@output[[1]]
}

K = as.numeric(k_str)
stopifnot(!is.na(K))

# read probs
if(verbose) { cat(" - reading probs (VCF)\n") }
if(verbose) { cat("  -- f_probs:", f_probs, "\n") }
vcf_probs = bigdat_bcf(f_probs)

fixed = vcf_probs$ftab
vals_anc = with(fixed[1, ], paste(REF, ALT, sep = ",")) %>% str_split(",") %>% unlist
stopifnot(length(vals_anc) == K)

# merge samples: exome & array (RFMix)
samples_probs = vcf_probs$rownames()
# samples = samples_probs 
samples_unrel = read_lines(f_unrel)
# stopifnot(all(samples_unrel %in% samples_probs))
samples = intersect(samples_unrel, samples_probs)
if(testing) {
  # samples = head(samples)
}
n = length(samples)

# vcf_probs$subset(keep = samples)

# split by batches
sp = vcf_probs$split(bsize)

# function to extract dosages of local ancestry
fun_str2num = function(x) strsplit(x, ",|\\|") %>% lapply(as.numeric) %>% do.call(rbind, .)

nb = nrow(sp)
if(testing) {
  nb = 10
}
out <- mclapply(seq(nb), function(b) {
  cat(" -", b, "/", nb, "\n")
  P = vcf_probs$batch(b, format = "AP") 
  P = P[samples, , drop = FALSE]
  stopifnot(nrow(P) == n)
  stopifnot(all(rownames(P) %in% samples))

  # extract dosages
  stopifnot(ncol(P) == 1)
  # mat_probs = fun_str2num(P[, 1])
  # mat_sum = colSums(mat_probs)
  mat_sum = fun_str2num(P[, 1]) %>% colSums 
  dos = sapply(seq(K), function(k) {
    cols = k + c(0, K)
    mat_sum[cols[1]] + mat_sum[cols[2]]
  })
  mat_dos = matrix(dos, nrow = 1)
  tab_dos = as_tibble(mat_dos)
  colnames(tab_dos) = vals_anc

  bind_cols(
    tibble(batch = b),
    tab_dos)
}, mc.cores = cores)
stopifnot(length(out) == nb)

# join individual site table together
tab = bind_rows(out)
# convert sum of prob. to proportions of anc.
tab = mutate_at(tab, vals_anc, function(x) x / (2*n))

## site annotations
# > fixed$ID %>% head(2)
# [1] "site_index_0_cm_1.78755" "site_index_5_cm_1.90107"
# extract position
pat <- rex(anything, "_", anything, "_", anything, "_", anything, "_", capture(anything))
annot = select(fixed, CHROM, POS, ID) %>%
  mutate(CM = re_matches(ID, pat) %>% .[[1]] %>% as.numeric) %>%
  select(-ID) %>%
  mutate(batch = seq(n())) %>% select(batch,  everything())

tab = left_join(annot, tab, by = "batch") %>% select(-batch)

## write resilts
write_tsv(tab, f_out)
