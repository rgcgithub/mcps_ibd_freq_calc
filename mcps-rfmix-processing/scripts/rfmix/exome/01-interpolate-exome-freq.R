# 13:19464376:T:A
# singletons in pgen (unphased)
# doubletons in vcf (phased)

library(tidyverse)
library(data.table)
library(glue)
library(unglue)

library(bedpca)

## par
vals_thr = c(0.5, 0.9, 0.99)
digits = 3 # signif. digits

## args
args = commandArgs(trailingOnly = TRUE)
no_args = is.na(args[1])

cat(" - args:", args, "\n")

f_variants = args[1]
f_unrel = args[2]
f_haps = args[3]
f_probs = args[4]
f_out = args[5]

# testing
if(no_args) {
  chr = 21
  b = 1
  f_variants = glue("output/lai2/phase-run_10222021_singletons/batches/{chr}/variants/variants.{b}.txt")
  f_haps = glue("output/lai/exome_phased/run_10222021_singletons/{chr}.vcf.gz")
  f_probs = glue("output/lai2/phase-run_10222021_singletons/rfmix-mais_1kg_hgdp_k7/vcf_probs/{chr}/{b}.batch.vcf.gz")
  f_unrel = "output/lai/samples/qc2.array.samples.txt"
  f_out = "tmp.tsv.gz"
}

# define k, vals_anc based on dirname_oxford
dirname_oxford = f_probs %>% dirname %>% dirname %>% dirname %>% basename
K = switch(dirname_oxford,
  "rfmix-mais_1kg_hgdp_k7" = 7,
  "rfmix-mais_1kg_hgdp_k3" = 3,
  stop("define K from dirname_oxford"))

vals_anc = switch(dirname_oxford,
  #reference_panel_population:    AFRICA  EUROPE  MEXICO_C        MEXICO_N        MEXICO_NW       MEXICO_S        MEXICO_SE
  "rfmix-mais_1kg_hgdp_k7" = c("AFRICA", "EUROPE", "MEXICO_C", "MEXICO_N", "MEXICO_NW", "MEXICO_S", "MEXICO_SE"),
  #reference_panel_population:    AFRICA  AMERICA EUROPE
  "rfmix-mais_1kg_hgdp_k3" = c("AFR", "MEX", "EUR"),
  stop("define K from dirname_oxford"))

## process parameters
variants_batch = read_lines(f_variants)
unrel <- read_lines(f_unrel)

## read probs
vcf_probs = bigdat_bcf(f_probs)
stopifnot(ncol(vcf_probs) == length(variants_batch))
fixed_probs = vcf_probs$fixed(variants_batch, elem = 1) # throws an error if any `variants_batch` is missing

## read haplotypes
vcf_haps = bigdat_bcf(f_haps)
fixed = vcf_haps$fixed(variants_batch, elem = 1) # throws an error if any `variants_batch` is missing

# fixed (exome)
samples_haps = vcf_haps$rownames()
samples_probs = vcf_probs$rownames()
samples = intersect(samples_haps, samples_probs)
samples <- samples[samples %in% unrel]
variants = fixed$ID
stopifnot(length(variants) == length(variants_batch))
stopifnot(all(variants %in% variants_batch))

## testing
if(no_args) {
  variants = head(variants, 5)
  # variants = variants[5] # singleton
}

n = length(samples)
p = length(variants)
vals_chr = fixed$CHR
stopifnot(all(vals_chr == vals_chr[1]))
chr = vals_chr[1]
anc = paste0("anc", 1:3)

## data matrices
H = vcf_haps[samples, variants]
P = vcf_probs[samples, variants, format = "AP"]
# check H and P are sync
stopifnot(all(rownames(H) == rownames(P)))
stopifnot(all(colnames(H) == colnames(P)))
# print object sizes
object.size(H) %>% print(u = "a")
object.size(P) %>% print(u = "a")

# functions to compute freq.
fun_freq_thr = function(cols, thr)
{
  lapply(thr, function(thr) {
    h = with(probs, c(h1, h2))
    p = c(probs[[cols[1]]], probs[[cols[2]]])
    ind = (p > thr)
    tibble(thr = thr, na = sum(ind), counts = sum(h[ind]), freq = mean(h[ind]))
  }) %>% bind_rows
}
fun_freq_avr = function(cols) 
{
  counts = sum(probs$h1 * probs[[cols[1]]] + probs$h2 * probs[[cols[2]]])
  len = sum(probs[[cols[1]]] + probs[[cols[2]]])
  tibble(na = len, counts = counts, freq = counts / len)
}

out = list()
for(i in seq_along(variants)) {
# freq = lapply(seq(2), function(i) {
  cat(" - freq. variant", i, "/", length(variants), "\n")
  v = variants[i]
  h = H[, v] # vector of "0|1", ...
  p = P[, v] # vector of "0,0.031,0.969|0,0.062,0.938", ...

  # haps
  h1 = substr(h, 1, 1) %>% as.numeric
  h2 = substr(h, 3, 3) %>% as.numeric
  haps = tibble(h1 = h1, h2 = h2)
  is_carrier = (h1 | h2)
  is_singleton = (sum(is_carrier) == 1)
  if(is_singleton) {
    cat("  -- singleton\n")
  }

  ## case 1: not a singleton
  if(!is_singleton) {
    # probs
    probs = str_split(p, ",|\\|") %>%
      lapply(as.numeric) %>% do.call(rbind, .) %>% as_tibble %>%
      bind_cols(haps)

    n = nrow(probs)
    counts_raw = with(probs, h1 + h2) %>% sum 
    freq_raw = counts_raw / (2*n)

    header = tibble(variant = v, singleton = 0)
    tab_header = header
    # freq. by thresholding probs
    tab_thr = lapply(seq(K), function(k) {
      offset = c(0, K)
      fun_freq_thr(k + offset, vals_thr) %>% mutate(anc = vals_anc[k], method = "thr")
    }) %>% bind_rows
    # freq. by avr.
    tab_avr = lapply(seq(K), function(k) {
      offset = c(0, K)
      fun_freq_avr(k + offset) %>% mutate(anc = vals_anc[k], method = "avr")
    }) %>% bind_rows
    # freq. raw
    tab_raw = tibble(anc = "ALL", method = "raw", counts = counts_raw, freq = freq_raw) %>% mutate(na = 2*n)

    # join all freq. results
    tab = bind_cols(tab_header, bind_rows(tab_thr, tab_avr, tab_raw))
    # add dummy columns
    tab = mutate(tab, p_ahom = NA, ghom = NA, ghom = NA)
  } else {
    # case 2: singleton
    pos = which(is_carrier)
    stopifnot(length(pos) == 1)
    gen = h1[pos] + h2[pos] # 1 or 2
    stopifnot(gen %in% c(1, 2))
    ghom = ifelse(gen == 2, 1, 0)

    # probs
    probs = str_split(p, ",|\\|") %>%
      lapply(as.numeric) %>% do.call(rbind, .) %>% as_tibble 
    probs1 = probs[pos, seq(K)] %>% as.numeric
    probs2 = probs[pos, K + seq(K)] %>% as.numeric

    counts = (0.5*gen)*(probs1 + probs2)
    len = sapply(seq(K), function(k) {
      sum(probs[[k]] + probs[[K + k]])
    })
    tab_anc = tibble(variant = v, singleton = 1L) %>%
      bind_cols(tibble(na = len, counts = counts, freq = counts / len)) %>%
      bind_cols(tibble(anc = vals_anc, method = "avr")) %>%
      bind_cols(tibble(ghom = ghom, p_ahom = probs1*probs2))

    n = length(h1)
    tab_raw = tibble(variant = v, singleton = 1L,
      na = 2*n, counts = gen, freq = 1 / (2*n), anc = "ALL", method = "raw",
        ghom = ghom, p_ahom = sum(probs1*probs2))

    tab = bind_rows(tab_anc, tab_raw)
    # add dummy columns
    tab = mutate(tab, thr = NA)
  }

  # check columns
  cols = c("variant", "singleton", "na", "counts", "freq",
    "anc", "method", "thr", "ghom", "p_ahom")
  stopifnot(all(cols %in% names(tab)))
  stopifnot(all(names(tab) %in% cols))
  # apply signif
  cols_counts = c("na", "counts")
  cols_probs = c("freq", "p_ahom")
  tab = mutate_at(tab, c(cols_counts, cols_probs), as.numeric)
  tab = mutate_at(tab, cols_counts, round, digits = digits)
  tab = mutate_at(tab, cols_probs, signif, digits = digits)
  # order columns
  tab = select(tab, all_of(cols))
  stopifnot(all(colnames(tab) == cols))

  # return
  out[[i]] = tab
}
# }) %>% bind_rows
freq = bind_rows(out)

## save
write_tsv(freq, f_out, na = "")
