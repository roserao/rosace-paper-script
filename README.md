# rosace-paper-script

Scripts and public datasets used to perform data analysis and generate plots for the paper "Rosace: a robust deep mutational scanning analysis framework employing position and mean-variance shrinkage".

For the full paper, see [link]().

<p align="left">
  <img src="/image/rosace_logo.png" width="150">
</p>

## Folder: realdata-main

#### Description

It includes raw and processed OCT1 data, processed MET data, and scripts to generate Figure 3, 4, S2, S3, and S4.

#### References

1. OCT1: Yee SW, Macdonald C, Mitrovic D, et al. The full spectrum of OCT1 (SLC22A1) mutations bridges transporter biophysics to drug pharmacogenomics. Preprint. bioRxiv. 2023;2023.06.06.543963. Published 2023 Jun 7. doi:10.1101/2023.06.06.543963
2. MET: Estevam GO, Linossi EM, Macdonald CB, et al. Conserved regulatory motifs in the juxtamembrane domain and kinase N-lobe revealed through deep mutational scanning of the MET receptor tyrosine kinase domain. eLife. 2023. doi: 10.7554/eLife.91619.1

#### Script

* run_fig3.R: format inputs and run Rosace, Enrich2, mustcan, and simple linear regression on all subsets of replicates in OCT1 drug cytotoxictiy screen from raw count (data/OCT1/count/1SM73)

* plot_fig3.R: summarise OCT1 functional score result from all methods and generate two plots in Figure 3 (rank-fdr and sensitivity).

* plot_fig4.R: generate Rosette object for OCT1 data and run a round of simulation. Plot Figures 4C (distribution of real vs. simulated functional score), 4D (estimated dispersion), and 4E (hierarchical clusteringo of mutant and histograms of score with clustering). Additionally, plot the Figures S2, S3, and S4.

* plot_fig4_met.R: repeat the procedure in plot_fig4.R on MET dataset.



## Folder: realdata-supplement

#### Description

It includes downloaded public data (BRCA1, CARD11, and MSH2) and scripts to generate Figure S1.

#### References
1. BRCA1: Findlay GM, Daza RM, Martin B, et al. Accurate classification of BRCA1 variants with saturation genome editing. Nature. 2018;562(7726):217-222. doi:10.1038/s41586-018-0461-z
2. CARD11: Meitlis I, Allenspach EJ, Bauman BM, et al. Multiplexed Functional Assessment of Genetic Variants in CARD11. Am J Hum Genet. 2020;107(6):1029-1043. doi:10.1016/j.ajhg.2020.10.015
3. MSH2: Jia X, Burugula BB, Chen V, et al. Massively parallel functional testing of MSH2 missense variants conferring Lynch syndrome risk. Am J Hum Genet. 2021;108(1):163-175. doi:10.1016/j.ajhg.2020.12.003

#### Script

* preprocess_{gene}.R: create cleaned Rosace object from the downloaded public excel sheets.

* workflow/scripts/prepare_method.R: source workflow/scripts/mutscan_utils.R for mutscan "calculateRelativeFC" function. Format inputs for Enrich2 and DiMSum and run mutscan.

* run.sh: run Enrich2 and DiMSum.

* workflow/scripts/compare_method.R: source workflow/scripts/compare_method_utils.R to extract results from Enrich2 output to a Rosace Score object. 

* analyze_{gene}.R: using the public validation set from the papers, generate plots in Figure S1 (rank-fdr, sensitibity, and specificity).



## Folder: simulation-snakemake

#### Description

A Snakemake pipeline to run Rosette simulation and benchmark methods (DiMSum, Enrich2, mutscan, and simple lienar regression).

Warning: this version of pipeline used an old version of Rosace, whose interface is slightly different from the most updated version (https://github.com/pimentellab/rosace). The key change that would screw up things is the interface of score data frame in Score object related to "workflow/scripts/compare_method.R": [variant, estimate, and test statistics] compared to [variant**s**, estimate, **standard error**, and test statistics].

#### Note

* flag index: pos (OCT1, position-favored Rosette), neg (MET, position-favored Rosette), posxxx (OCT1, default Rosette), and negxxx (MET, default Rosette).

#### Input

* Rosette object: 1SM73 (OCT1 drug screen) and MET. Manually generated beforehand with tested configs from plot_fig4.R in folder realdata-main.

#### Output 

* Folder results/sim_{flag}: files too large. Not included in the github. Each subfolder is a condition for the number of replicate and selection round. Each sub-subfolder is a round of simulation, which contains both the raw and summarized results for all methods.
* Folder results/summary_{flag}: summarize to the sub-folder level mentioned above (condition with number of replicate and selection round).
* Folder results/output_{flag}: summarize to the flag index level.



## Folder: simulation-analysis

#### Description

Take data input from folder simulation-snakemake/output_{neg, negxxx, pos, posxxx}. 

#### Script

* plot_fig5.R: plot Figure 5 and Figure S6 (rank-fdr).

* plot_fig6.R

* plot_supfig.R: plot Figures S5 (correlation) and S7 (run time and resources).

