library(tidyverse)
library(data.table)
library(glue)

library(cowplot)
theme_set(theme_minimal_grid(8))

## par.
levs_anc = c("AFRICA", "AMERICA", "EUROPE")
cols = c(AFRICA = "#FDE725FF", EUROPE = "#238A8DFF", AMERICA = "#94D840FF")


## chr 15
lanc = glue('output/rfmix/processed/panel-rgc_1kg_hgdp/param-rbase/lanc-vcf-raw/lanc.unrel.chr15.tsv.gz') %>% fread
tab15 = pivot_longer(lanc, cols = levs_anc, names_to = "ANC", values_to = "PROP")

tab15a = filter(tab15, ANC == 'AFRICA')

# fixed version
flanc = lanc[-seq(3), ]
ftab15 = pivot_longer(flanc, cols = levs_anc, names_to = "ANC", values_to = "PROP")

ftab15a = filter(ftab15, ANC == 'AFRICA')

fun_plot15 = function(tab)
  ggplot(tab, aes(POS, PROP, group = ANC)) +
    geom_line(aes(color = ANC), size = 0.5) +
    geom_point(aes(color = ANC), size = 1) +
    scale_color_manual(values = cols[levs_anc], drop = TRUE) +
    scale_x_continuous(limits = c(19e6, 30e6), labels = scales::comma) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = NULL)

g15 = plot_grid(
  fun_plot15(tab15) + labs(title = 'Raw', subtitle = 'Oxford, 1KG/HGDP/MAIS Ref., K = 3, rbase parameter set') + theme(legend.position = 'none'),
  fun_plot15(tab15a) + labs(title = 'Raw', subtitle = 'AFR ancestry only'),
  fun_plot15(ftab15) + labs(title = 'Clean', subtitle = 'first n = 1 CRF points removed') + theme(legend.position = 'none'),
  fun_plot15(ftab15a) + labs(title = 'Clean', subtitle = 'AFR ancestry only'),
  ncol = 2, rel_widths = c(1, 1.2))

## save plot
ggsave('tmp.png', g15, dpi = 100)

