# rosace-paper-script

Scripts and public datasets used to perform data analysis and generate plots for the paper "Rosace: a robust deep mutational scanning analysis framework employing position and mean-variance shrinkage".

For the full paper, see [link]().

<p align="left">
  <img src="/image/rosace_logo.png" width="150">
</p>

## Folder: realdata-snakemake

#### Description

This folder is a snakemake pipeline to perform real-data analaysis in the paper (OCT1, MET, MSH2, CARD11, BRCA1, BRCA1-RING, and Cohesin). It include raw and processed data, scripts to generate Figure 3, 4, S2-17, and summarized data frames and plots. 

#### References

1. OCT1: Yee SW, Macdonald C, Mitrovic D, et al. The full spectrum of OCT1 (SLC22A1) mutations bridges transporter biophysics to drug pharmacogenomics. Preprint. bioRxiv. 2023;2023.06.06.543963. Published 2023 Jun 7. doi:10.1101/2023.06.06.543963
2. MET: Estevam GO, Linossi EM, Macdonald CB, et al. Conserved regulatory motifs in the juxtamembrane domain and kinase N-lobe revealed through deep mutational scanning of the MET receptor tyrosine kinase domain. eLife. 2023. doi: 10.7554/eLife.91619.1
3. BRCA1: Findlay GM, Daza RM, Martin B, et al. Accurate classification of BRCA1 variants with saturation genome editing. Nature. 2018;562(7726):217-222. doi:10.1038/s41586-018-0461-z
4. CARD11: Meitlis I, Allenspach EJ, Bauman BM, et al. Multiplexed Functional Assessment of Genetic Variants in CARD11. Am J Hum Genet. 2020;107(6):1029-1043. doi:10.1016/j.ajhg.2020.10.015
5. MSH2: Jia X, Burugula BB, Chen V, et al. Massively parallel functional testing of MSH2 missense variants conferring Lynch syndrome risk. Am J Hum Genet. 2021;108(1):163-175. doi:10.1016/j.ajhg.2020.12.003
6. BRCA1-RING: Starita LM, Young DL, Islam M, Kitzman JO, Gullingsrud J, Hause RJ, Fowler DM, Parvin JD, Shendure J, Fields S. Massively Parallel Functional Analysis of BRCA1 RING Domain Variants. Genetics. 2015 Jun;200(2):413-22. doi: 10.1534/genetics.115.175802. Epub 2015 Mar 30. Erratum in: Genetics. 2017 Dec;207 (4):1713. PMID: 25823446; PMCID: PMC4492368.
7. Cohesin: Kowalsky CA, Whitehead TA. Determination of binding affinity upon mutation for type I dockerin-cohesin complexes from Clostridium thermocellum and Clostridium cellulolyticum using deep sequencing. Proteins. 2016 Dec;84(12):1914-1928. doi: 10.1002/prot.25175. Epub 2016 Oct 26. PMID: 27699856.

#### Results

* Folder results/{OCT1, MET, BRCA1, CARD11, MSH2, BRCA1-RING, Cohesin}: analysis on the complete data.
* Folder replicates/{OCT1, MET}/{rep1_1, rep1_2, rep1_3, rep2_1, rep2_2, rep2_3, rep3_1}: replicate-level analysis on OCT1 and MET only. 

Each folder contains method subfolders {rosace, rosace_nopos, dimsum, enrich2, mutscan} that are not included in the github due to their large sizes. The summarised output is in "analysis/rosace_full.rda", which contains the outputs of all methods.

Each folder contains subfolders {analysis, plot}. The "analysis" folder stores method comparsion with several metrics available (synonymous mutations, EVE, clinvar, and AlphaMissense) that are used to make the plots in the "plot" folder.

Moreover, the replicate folder has a subfolder "raw" with partial-replicate data that are being analyzed in the folder.


#### Script

All the following files are in the folder workflow/scripts.

* prepare_method.R: format inputs for Enrich2, DiMSum and mutscan; run mustcan.

* run_rosace.R: run rosace and rosace (no position mode). 

* compare_method.R: run the naive method, and summarise all method outputs into one rosace object in the subfolder "analyis".

* analysis.R: compute metrics and generate plots.

* parse_replicate.R: prepare input for replicate-level analysis. Specifically, generate partial-replicate rosace object in the subfolder "raw".

* summarize_replicate.R: summarize replicate-level metrics from analaysis and generate plots in the replicates/{OCT1, MET}/summary.

* plot_mainfig.R: plot main Figure 3 and 4.


## Folder: simulation-snakemake

#### Description

A Snakemake pipeline to run Rosette simulation and benchmark methods (DiMSum, Enrich2, mutscan, Rosace, Rosace(nopos), and Naive). A full run of simdms simulation (from the Enrich2 paper) is alos included.

#### simdms Simulation (Enrich2)

The archived version of simdms on zenodo (10.5281/zenodo.546311) is downloaded in workflow/simdms. The simulation is run on the default parameters. 

- results/simdms/{growth, binding}_simulation_depth_200: raw output of simdms simualtion.

- workflow/scripts/create_rosace_simdms.R: create a rosace object from the ouput of simdms.
  - results/simdms/{growth, binding}/{clean, jackpot, reperror}: each simdms simulation have three modes named clean, jackpot, and reperror. For each mode, a rosace object is created.

#### Rosette Simulation

The input of rosette simulation is the rosace object in the data folder (OCT1, MET, and CARD11). 

- workflow/scripts/rosette.R: genereate rosette object for {OCT1, MET, CARD11} and compare the difference between real data and simluation. The result is infolders "results/rosette/{OCT1, MET, CARD11}". The plots in the folders are arranged in Figure S18-20.

- workflow/scripts/run_rosette.R: for each dataset, run rosette simulation in 8 different modes. These raw simulation folders are not uploaded on github due to size.
  - results/{OCT1, MET, CARD11}/{pos, nopos}/growth_rep{1, 3}_rd{1, 3}_clean: "pos" and "nopos" indicates the position-favored rosette or the original rosette simulation. Furthermore, simulation is run with 1 or 3 number of replicates and 1 or 3 number of selection rounds.
  - Each folder mentioned above contains 10 simulation rounds. 

#### Analysis Pipeline

For each folder mentioned above, similar analysis pipeline is performed. The scripts are in "workflow/scripts".

- prepare_method_simdms.R and prepare_method.R: format inputs for Enrich2, DiMSum and mutscan; run mustcan.

- run_rosace_simdms.R and run_rosace.R: run rosace and rosace (no position mode). 

- compare_method_simdms.R and compare_method.R: run the naive method, and summarise all method outputs into one rosace object in the subfolder "analyis".

- analysis_simdms.R and analysis.R: compute metrics and generate plots.

#### Plots

- plot_gene.R: gene-level (OCT1, MET, CARD11) summary plot to generate Figure S21-S26.

- analysis_simdms.R: contains scripts used to generate Figure S27.

- plot_mianfig.R: scripts used to generate Figure 5 and 6.

