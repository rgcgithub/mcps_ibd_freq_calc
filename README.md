This repository contains scripts code to compute ancestry-specific allele
frequencies from Ziyatdinov, A., Torres, J., Alegre-Diaz, J. et a. (2022).
Genotyping, sequencing and analysis of 140,000 adults from the Mexico City
Prospective Study. bioRxiv, 2022-06.
<https://doi.org/10.1101/2022.06.26.495014>

The scripts in `mcps-rfmix-processing` extrapolate local ancestry estimates
from RFMix assignments using array variants to WES/WGS sites.

The program in `calc-ibd-freq` uses these estimates to compute raw and
relatedness corrected allele frequencies.
