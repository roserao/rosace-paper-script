# Rscript workflow/scripts/prepare_method.R --data CARD11
# Rscript workflow/scripts/prepare_method.R --data MSH2
# Rscript workflow/scripts/prepare_method.R --data BRCA1-RING
# Rscript workflow/scripts/prepare_method.R --data Cohesin

# Rscript workflow/scripts/run_rosace.R --data BRCA1

# conda activate enrich2
# enrich_cmd results/CARD11/enrich2/config.json ratios wt --no-plots
# enrich_cmd results/MSH2/enrich2/config.json ratios wt --no-plots
# enrich_cmd results/BRCA1/enrich2/config.json WLS wt --no-plots
# enrich_cmd results/OCT1/enrich2/config.json WLS wt --no-plots
# enrich_cmd results/BRCA1-RING/enrich2/config.json WLS wt --no-plots

# Rscript workflow/scripts/compare_method.R --data MSH2
# Rscript workflow/scripts/compare_method.R --data BRCA1

# Rscript workflow/scripts/analysis.R --d CARD11 
# Rscript workflow/scripts/analysis.R --d OCT1
# Rscript workflow/scripts/analysis.R --d MSH2
# Rscript workflow/scripts/analysis.R --d BRCA1
# Rscript workflow/scripts/analysis.R --d MET


rule prepare_method:
    input: 
        "data/{data}/rosace.rda"
    output:
        "results/{data}/mutscan/score_edgeR.tsv",
        "results/{data}/mutscan/score_limma.tsv",
        "results/{data}/enrich2/config.json",
        "results/{data}/dimsum/dimsum_count.tsv"
    shell:
        """
        Rscript workflow/scripts/prepare_method.R \
            --data {wildcards.data} 
        """

rule run_dimsum:
    input: 
        "results/{data}/dimsum/dimsum_count.tsv"
    output:
        "results/{data}/dimsum/DiMSum_Project/DiMSum_Project_fitness_replicates.RData"
    conda:  
        'dimsum'
    shell:
        """
        module unload R
        REF="$(cat results/{wildcards.data}/dimsum/dimsum_ref.txt)"
        DiMSum --startStage 4 \
            --countPath results/{wildcards.data}/dimsum/dimsum_count.tsv \
            --experimentDesignPath results/{wildcards.data}/dimsum/exp_design.tsv \
            --wildtypeSequence $REF \
            --outputPath results/{wildcards.data}/dimsum/ \
            --mutagenesisType random \
            --maxSubstitutions 6 \
            --sequenceType noncoding
        """

rule run_rosace: 
    input:
        "data/{data}/rosace.rda"
    output:
        "results/{data}/rosace/rosace_eval.rda"
    log: 
        stdout = "logs/{data}/rosace.stdout",
        stderr =  "logs/{data}/rosace.stderr"
    threads: 4
    benchmark:
        "results/{data}/rosace/benchmark.txt"
    shell:
        """
        Rscript workflow/scripts/run_rosace.R \
            --data {wildcards.data} \
            > {log.stdout} 2> {log.stderr}
        """

rule run_rosace_nopos: 
    input:
        "data/{data}/rosace.rda"
    output:
        "results/{data}/rosace_nopos/rosace_eval.rda"
    log: 
        stdout = "logs/{data}/rosace_nopos.stdout",
        stderr =  "logs/{data}/rosace_nopos.stderr"
    threads: 4
    benchmark:
        "results/{data}/rosace_nopos/benchmark.txt"
    shell:
        """
        Rscript workflow/scripts/run_rosace.R \
            --data {wildcards.data} --mode N \
            > {log.stdout} 2> {log.stderr}
        """

rule run_enrich2:
    input:
        expand("results/{data}/enrich2/config.json", data = v_data)
    output:
        expand("results/{data}/enrich2/results/tsv/simulation_exp/main_identifiers_scores.tsv", data = v_data),
        expand("results/{data}/enrich2/results/tsv/simulation_exp/main_identifiers_scores_pvalues_wt.tsv", data = v_data)
    conda:  
        'enrich2'
    shell:
        """
        enrich_cmd results/MET/enrich2/config.json WLS wt --no-plots
        enrich_cmd results/OCT1/enrich2/config.json WLS wt --no-plots
        enrich_cmd results/CARD11/enrich2/config.json ratios wt --no-plots
        enrich_cmd results/MSH2/enrich2/config.json ratios wt --no-plots
        enrich_cmd results/BRCA1/enrich2/config.json WLS wt --no-plots
        enrich_cmd results/BRCA1-RING/enrich2/config.json WLS wt --no-plots
        enrich_cmd results/Cohesin/enrich2/config.json ratios complete --no-plots
        """

rule compare_method:
    input: 
        "results/{data}/mutscan/score_edgeR.tsv",
        "results/{data}/mutscan/score_limma.tsv",
        "results/{data}/rosace/rosace_eval.rda",
        "results/{data}/rosace_nopos/rosace_eval.rda",
        "results/{data}/dimsum/DiMSum_Project/DiMSum_Project_fitness_replicates.RData",
        "results/{data}/enrich2/results/tsv/simulation_exp/main_identifiers_scores.tsv",
        "results/{data}/enrich2/results/tsv/simulation_exp/main_identifiers_scores_pvalues_wt.tsv"
    output:
        "results/{data}/analysis/rosace_full.rda",
        "results/{data}/analysis/effects.tsv"
    shell:
        """
        Rscript workflow/scripts/compare_method.R --data {wildcards.data} 
        """

rule analysis:
    input: 
        "results/{data}/analysis/effects.tsv",
        "data/sysdata.rda"
    output:
        expand("results/{{data}}/plot/{plot}", plot = ["rank_fdr_syn.png", "gnomad_roc.png"])
    shell:
        """
        Rscript workflow/scripts/analysis.R --data {wildcards.data} 
        """
# "data/{data}/gnomAD_cleaned.tsv"

# Rscript workflow/scripts/analysis.R --data Cohesin
# Rscript workflow/scripts/analysis.R --data Cohesin
# Rscript workflow/scripts/analysis.R --data Cohesin
# Rscript workflow/scripts/analysis.R --data Cohesin
# Rscript workflow/scripts/analysis.R --data Cohesin
# Rscript workflow/scripts/analysis.R --data Cohesin



