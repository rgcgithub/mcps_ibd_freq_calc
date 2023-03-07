library(tidyverse)
library(data.table)
library(glue)

## par
vals_anc = c('AFRICA', 'AMERICA', 'EUROPE')
vals_chr = 1:22

## args
no_snakemake = !exists("snakemake")
testing = no_snakemake

if(testing) {
  vals_chr = 21:22
}

# testing
# output: 'output/rfmix/processed/panel-{panel}/param-{param}/ganc.tsv.gz'
if(no_snakemake) {
  panel = 'rgc_1kg_hgdp'; param = 'rbase'
  out = "tmp.tsv.gz"
} else {
  panel = snakemake@wildcards[['panel']]
  param = snakemake@wildcards[['param']]
  f_out = snakemake@output[[1]]
}


path_rfmix = switch(panel,
  'rgc_1kg_hgdp' = glue('rfmix/ref-1kg_hgdp-K3/run-{param}/'),
    stop('panel'))
fun_rfmix_q = function(chr) {
  glue('{path_rfmix}/chr-{chr}/out.rfmix.Q')
}

# Sex code (‘1’ = male, ‘2’ = female, ‘0’ = unknown)
fam = "output/MC150K/qc-array/final/qc2.fam" %>% 
  fread %>% select(V2, V5) %>% rename(id = V2, sex = V5)
s_fam = fam$id %>% as.character

f_unrel = "output/MC150K/rel2/v2-king-ibdseg/degree-3rd/ivs.unrel.samples.txt"
s_unrel = read_lines(f_unrel)

## autosomes
out = list()
chr_len_total = 0

for(x in vals_chr) {
  print(x)
  file_q = fun_rfmix_q(x)
  stopifnot(file.exists(file_q))
  q = fread(file_q)
  names(q)[1] = 'id'

  q_mcps = mutate(q, ref = grepl('-', id)) 
  stop()
  q_mcps = mutate(q, ref = grepl('-', id)) %>% filter(!ref) %>% select(-ref)
  stopifnot(nrow(q_mcps) == length(s_fam))
  stopifnot(all(s_fam %in% q_mcps$id))
  stopifnot(colnames(q_mcps) == c('id', vals_anc))

  # re-order ids in Q based on s_fam
  dat = tibble(id = s_fam) %>% left_join(q_mcps) %>% select(-id)
  stopifnot(ncol(dat) == 3)
  stopifnot(nrow(dat) == length(s_fam))

  f = glue('~/data/maps/chr{x}.b38.gmap.gz') 
  chr_len = fread(f) %$% cM %>% range %>% diff
  dat_chr = dat %>% mutate_all(function(x) chr_len*x)

  out[[x]] = dat_chr
  chr_len_total = chr_len_total + chr_len
}

## average
dat = out[[vals_chr[1]]]
for(i in seq(2, length(vals_chr))) {
  dat = dat + out[[vals_chr[i]]]
}
dat = dat / chr_len_total
dat = as_tibble(dat)
dat = mutate(dat, ID = s_fam) %>%
  mutate(UNREL = (ID %in% s_unrel)) %>%
  select(ID, UNREL, everything())


## write results
write_tsv(dat, f_out)
