library(tidyverse)
library(data.table)
library(glue)
library(unglue)

library(bedpca)



## par
digits = 3
verbose = 1

## args
no_snakemake = !exists("snakemake")
testing = no_snakemake

if(!no_snakemake) {
  ## 1 thread for use by data.table
  setDTthreads(1)
  ## Sleep for 1-5 min
  min_sleep = sample(seq(5), 1)
  Sys.sleep(60*min_sleep) # in sec
}

# input variables: variants, haps, probs, k, out
if(no_snakemake) {
  # K = 7
  # f_variants = 'output/rfmix/exome/panel-oxford_mais_k7/param-rbase/probs-raw_wgs_nodup/chr-1/haps_exome.chr1.variants.txt'
  # f_splits = 'output/rfmix/exome/panel-oxford_mais_k7/param-rbase/probs-raw_wgs_nodup/chr-1/splits.tsv.gz'
  # f_haps = 'output/rfmix/input/wgs_nodup_haps/haps.chr1.vcf.gz'
  # f_probs = 'output/rfmix/processed/panel-oxford_mais_k7/param-rbase/probs-vcf-raw/probs.chr1.vcf.gz'
  # b_str = '7985'
  # k_str = '7'

  # Chr X, K = 3
  f_variants = 'output/rfmix/exome/panel-rgc_1kg_hgdp/param-rbase/probs-raw/chr-X/haps_exome.chrX.variants.txt'
  f_splits = 'output/rfmix/exome/panel-rgc_1kg_hgdp/param-rbase/probs-raw/chr-X/splits.tsv.gz'
  f_haps = 'output/rfmix/input/exome_haps/haps.chrX.vcf.gz'
  f_probs = 'output/rfmix/processed/panel-rgc_1kg_hgdp/param-rbase/probs-vcf-raw/probs.chrX.vcf.gz'
  b_str = '1'
  k_str = '3'

  f_out = 'tmp.vcf'
} else {
  f_variants = snakemake@input[['variants']]
  f_splits = snakemake@input[['splits']]
  f_haps = snakemake@input[['haps']]
  f_probs = snakemake@input[['probs']]
  b_str = snakemake@wildcards[['b']]
  k_str = snakemake@params[['k']]
  f_out = snakemake@output[[1]]
}

K = as.numeric(k_str)
stopifnot(!is.na(K))

## variants
b = as.integer(b_str)
sp = read_tsv(f_splits)
stopifnot(b <= nrow(sp))
ind = seq(sp$beg[b], sp$end[b])
v = read_lines(f_variants)
stopifnot(max(ind) <= length(v))
variants_batch = v[ind]
rm(v); gc()

## read haplotypes
vcf_haps = bigdat_bcf(f_haps)
gc()
fixed = vcf_haps$fixed(variants_batch, elem = 1) # throws an error if any `variants_batch` is missing
gc()

## samples (exome/wgs)
samples_exome = vcf_haps$rownames()
rm(vcf_haps); gc()
# case of WGS sample IDs
if(substr(samples_exome[1], 1, 4) == "MCPS") {
  f_mapid = "input/MCPS-Freeze150K-SampleInfo-Genotypying-QCPass.csv"
  s_patid = samples_exome
  # manual adjustment
  s_patid[s_patid == "MCPS_MCPSRGN057303_MEXB059604"] = "MCPS_MCPSRGN057303_MEXB059604-DUP"
  map = fread(f_mapid, colClasses = "character") %>% as_tibble
  names(map) = c("patid", "ID")
  tab = tibble(patid = s_patid) %>% left_join(map, by = "patid")
  stopifnot(!any(is.na(tab$ID)))
  samples_exome = tab$ID
}

# fixed (exome)
variants = fixed$ID
stopifnot(length(variants) == length(variants_batch))
stopifnot(all(variants %in% variants_batch))

p = length(variants)
vals_chr = fixed$CHR
stopifnot(all(vals_chr == vals_chr[1]))
chr = vals_chr[1]

# read probs
if(verbose) { cat(" - reading probs (VCF)\n") }
if(verbose) { cat("  -- f_probs:", f_probs, "\n") }
vcf_probs = bigdat_bcf(f_probs)
gc()
fixed_array = vcf_probs$ftab
# fix of CHROM name for Chr. X: X -> 23 
# phased exomes encode Chr. X as 23 
fixed_array = mutate(fixed_array, CHROM = ifelse(CHROM == 'X', '23', CHROM))

anc_str1 = fixed_array$REF[1]
anc_str2 = fixed_array$ALT[1]

# merge samples: exome & array (RFMix)
samples_probs = vcf_probs$rownames()
samples = intersect(samples_exome, samples_probs)
n = length(samples)
stopifnot(n > 0)

# merge variants: fixed (exome) & fixed_array
if(verbose) { cat(" - merging fixed of exome and array\n") }
t = bind_rows(mutate(fixed, elem = 0), fixed_array) %>%
  select(CHROM, POS, ID, elem, index) 
setkey(t, POS)
setorder(t, POS)
t = mutate(t, row = seq(n()))

t0 = filter(t, elem == 0)
rmin = min(t0$row)
rmax = max(t0$row)
rmin = max(0, rmin - 1)
rmax = min(nrow(t), rmax + 1)

t = t[row >= rmin & row <= rmax]

t1 = filter(t, elem == 1) 
setorder(t1, POS)
variants_array = t1$ID
P = vcf_probs[samples, variants_array, format = "AP"]
rm(vcf_probs); gc()
stopifnot(all(rownames(P) == samples))
stopifnot(all(colnames(P) == variants_array))

## convert matrix of probs. (P) to hapltoypes P1/P2
fun_hap1 = function(x) unglue_vec(x, "{x}|{}")
fun_hap2 = function(x) unglue_vec(x, "{}|{x}")
# P1 = apply(P, 2, fun_hap1)
# P2 = apply(P, 2, fun_hap2)

fun_str2num = function(x) strsplit(x, ",|\\|") %>% lapply(as.numeric) %>% do.call(rbind, .)
fun_num2str = function(X) 
{
  h1 = X[, seq(K)] %>% round(digits) %>% apply(1, paste, collapse = ",")
  h2 = X[, K + seq(K)] %>% round(digits) %>% apply(1, paste, collapse = ",")
  paste(h1, h2, sep = "|")
}

# write header
if(verbose) { cat(" - writing VCF header\n") }
lines_header = c(
  "##fileformat=VCFv4.2",
  glue("##contig=<ID={chr}>"),
  "##FORMAT=<ID=AP,Number=G,Type=String,Description=\"Allele probabilities\">")
line_colnames = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
line_samples = paste(samples, collapse = "\t")
line_names = paste0(line_colnames, line_samples)

write_lines(lines_header, f_out)
write_lines(line_names, f_out, append = TRUE)

# write data records to VCF (probs)
if(verbose) { cat(" - writing VCF records\n") }
for(i in seq_along(variants)) {
# for(i in seq(2)) {
  if(verbose) cat(" - writing variant ", i, "/", p, "\n")
  pos = as.integer(fixed$POS[i])
  t1_left = filter(t1, POS <= pos)
  t1_right = filter(t1, POS >= pos)
  t1_left = tail(t1_left, 1)
  t1_right = head(t1_right, 1)
  left = t1_left$ID 
  right = t1_right$ID
  w_left = (pos - t1_left$POS)
  w_right = (t1_right$POS - pos) 
  if(nrow(t1_left) > 0 & nrow(t1_right) > 0) {
    # case 0: left == right
    if(left == right) {
      P_target = P[, left] %>% fun_str2num
    } else {
      # case 1: left & right boundaries exist
      if(verbose > 1) { cat("  -- case 1\n") }
      stopifnot(w_left >= 0)
      stopifnot(w_right >= 0)
      P_left = P[, left] %>% fun_str2num
      P_right = P[, right] %>% fun_str2num
      P_target = (w_left * P_left + w_right * P_right) / (w_left + w_right)
    }
  } else if(nrow(t1_left) > 0 & nrow(t1_right) == 0) {
    # case 2: left only
    if(verbose > 1) { cat("  -- case 2\n") }
    P_target = P[, left] %>% fun_str2num
  } else if(nrow(t1_left) == 0 & nrow(t1_right) > 0) {
    # case 3: right only
    if(verbose > 1) { cat("  -- case 3\n") }
    P_target = P[, right] %>% fun_str2num
  } else {
    stop("both left and right are missing")
    # str_dat = paste(rep(".,.,.|.,.,.", n), collapse = "\t")
  }
  stopifnot(!any(is.na(P_target)))
  probs_target = fun_num2str(P_target)
  str_dat = paste(probs_target, collapse = "\t")

  str_info = glue("{chr}\t{as.integer(pos)}\t{fixed$ID[i]}\t{anc_str1}\t{anc_str2}\t.\t.\t.\tAP\t")
  str_record = paste0(str_info, str_dat)
  # if(verbose > 1) { cat("'", str_record, "'", "\n") }
  if(verbose > 1) { cat("'", str_info, "'", "\n") }
  if(verbose > 1) { cat("'", substr(str_dat, 1, 100), "'", "\n") }
  write_lines(str_record, f_out, append = TRUE)
}

## check #variants
stopifnot(length(variants_batch) == vcftools_n_variants(f_out))
