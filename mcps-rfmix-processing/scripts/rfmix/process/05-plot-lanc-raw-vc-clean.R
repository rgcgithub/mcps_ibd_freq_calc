library(tidyverse)
library(data.table)
library(glue)

library(cowplot)
theme_set(theme_minimal_grid(8))

## par.
K = 7
vals_chr = c(15, 19)
levs_anc = c("AFRICA", "AMERICA", "EUROPE")
cols = c(AFRICA = "#FDE725FF", EUROPE = "#238A8DFF", AMERICA = "#94D840FF",
  "MEXICO_N" = "#E31A1C", "MEXICO_NW" = "#FF7F00",
  "MEXICO_C" = "#FDBF6F", "MEXICO_S" = "#FB9A99", "MEXICO_SE" = "#B3367AFF")

if(K == 7) { 
  levs_anc = c("AFRICA", "MEXICO_C", "MEXICO_N", "MEXICO_NW", "MEXICO_S", "MEXICO_SE", "EUROPE", "AMERICA") 
}

## raw
lanc = lapply(vals_chr, function(x) {
  glue('output/rfmix/processed/panel-oxford_mais_k{K}/param-rbase/lanc-vcf-raw/lanc.unrel.chr{x}.tsv.gz') %>%
    fread %>% mutate(chr = x)
}) %>% bind_rows
if(K == 7) { lanc = mutate(lanc, AMERICA = 1 - AFRICA - EUROPE) }
tab = pivot_longer(lanc, cols = levs_anc, names_to = "ANC", values_to = "PROP")

# fixed
lanc = lapply(vals_chr, function(x) {
  glue('output/rfmix/processed/panel-oxford_mais_k{K}/param-rbase/lanc-vcf-clean/lanc.unrel.chr{x}.tsv.gz') %>%
    fread %>% mutate(chr = x)
}) %>% bind_rows
if(K == 7) { lanc = mutate(lanc, AMERICA = 1 - AFRICA - EUROPE) }
ftab = pivot_longer(lanc, cols = levs_anc, names_to = "ANC", values_to = "PROP")

fun_plot = function(tab)
  tab %>% mutate(chr = paste("Chr.", chr)) %>%
    # filter(POS > 130e6) %>%
    ggplot(aes(POS, PROP, group = ANC)) +
      geom_line(aes(color = ANC), size = 0.5) +
      geom_point(aes(color = ANC), size = 1) +
      scale_color_manual(values = cols[levs_anc], drop = TRUE) +
      scale_x_continuous(labels = scales::comma) +
      scale_y_continuous(labels = scales::percent) +
      labs(x = NULL, y = NULL) +
      facet_wrap(~ chr, scales = 'free_x')

g = plot_grid(
  fun_plot(tab) + labs(title = 'Raw') + theme(legend.position = 'none'),
  fun_plot(ftab) + labs(title = 'Clean') + theme(legend.position = 'none'),
  nrow = 2, rel_widths = c(1, 1.2))

## save plot
ggsave('tmp.png', g, dpi = 100)

