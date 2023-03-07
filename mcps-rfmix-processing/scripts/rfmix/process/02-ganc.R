library(tidyverse)
library(data.table)
library(glue)

## args
no_snakemake = !exists("snakemake")
testing = no_snakemake

# testing
if(no_snakemake) {
  f_unrel = "output/MC150K/rel2/v2-king-ibdseg/degree-3rd/ivs.unrel.samples.txt"
  f_out = "tmp.tsv.gz"
} else {
  f_unrel = snakemake@input[[1]]
  f_out = snakemake@output[[1]]
}

K = 3
s_unrel = read_lines(f_unrel)
# Sex code (‘1’ = male, ‘2’ = female, ‘0’ = unknown)
fam = "output/MC150K/qc-array/final/qc2.fam" %>% 
  fread %>% select(V1, V5) %>% rename(id = V1, sex = V5)

## Chr X
q = 'runs/chrX-fixref-rfmix/output/rfmix/run-base0/out.rfmix.Q' %>% fread
names(q)[1] = 'id'
q = mutate(q, id = as.character(id))
table(s_unrel %in% q$id) # 51 IDs are out
q = filter(q, id %in% s_unrel)

qm = left_join(q, fam, by = 'id') %>% filter(sex == 1) %>% select(-sex)
qf = left_join(q, fam, by = 'id') %>% filter(sex == 2) %>% select(-sex)

# q_avr = (colMeans(qm[, -1]) + 2*colMeans(qf[, -1]))/1.5
q_avr = colMeans(qm[, -1])
ganc_chrx = tibble(chr = "X", anc = names(q_avr), ganc_mean = q_avr)

## autosomes
fun_file_q = function(chr) glue('rfmix_fixref_mais_k3/mcps.rfmix.chr{chr}.rfmix.Q')
ganc_chr_auto = lapply(1:22, function(x) {
  f = fun_file_q(x)
  q = fread(f)
  names(q)[1] = 'id'
  stopifnot(all(s_unrel %in% q$id))
  q = filter(q, id %in% s_unrel)
  q_avr = colMeans(q[, -1])

  f = glue('~/data/maps/chr{x}.b38.gmap.gz') 
  chr_len = fread(f) %$% cM %>% range %>% diff

  tibble(chr = x, anc = names(q_avr), ganc_mean = q_avr, chr_len = chr_len)
}) %>% bind_rows

## all auto chr
ganc_auto = group_split(ganc_chr_auto, anc) %>% 
  lapply(function(x) tibble(ganc_mean = sum(x$ganc_mean * x$chr_len) / sum(x$chr_len), anc = unique(x$anc))) %>% bind_rows %>%
  mutate(chr = "auto")

ganc = bind_rows(
  mutate(ganc_chr_auto, chr = as.character(chr)),
  ganc_auto,
  ganc_chrx
)

## write resilts
# write_tsv(tab, f_out)
