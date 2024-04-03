rule parse_replicate:
    input: 
        "data/{data}/rosace.rda"
    output:
        expand("replicates/{{data}}/{mode}/raw/rosace.rda", mode = v_mode)
    shell:
        """
        Rscript workflow/scripts/parse_repliacte.R --data {wildcards.data} 
        """

# Rscript workflow/scripts/parse_replicate.R --data MET
# Rscript workflow/scripts/parse_replicate.R --data OCT1

rule prepare_method_replicate:
    input: 
        "replicates/{data}/{mode}/raw/rosace.rda"
    output:
        "replicates/{data}/{mode}/mutscan/score_edgeR.tsv",
        "replicates/{data}/{mode}/mutscan/score_limma.tsv",
        "replicates/{data}/{mode}/enrich2/config.json",
        "replicates/{data}/{mode}/dimsum/dimsum_count.tsv"
    shell:
        """
        Rscript workflow/scripts/prepare_method.R \
            --data {wildcards.data} --mode {wildcards.mode}
        """

# Rscript workflow/scripts/prepare_method.R --data OCT1 --mode rep2_2

rule run_dimsum_replicate:
    input: 
        "replicates/{data}/{mode}/dimsum/dimsum_count.tsv"
    output:
        "replicates/{data}/{mode}/dimsum/DiMSum_Project/DiMSum_Project_fitness_replicates.RData"
    conda:  
        'dimsum'
    shell:
        """
        module unload R
        REF="$(cat replicates/{wildcards.data}/{wildcards.mode}/dimsum/dimsum_ref.txt)"
        DiMSum --startStage 4 \
            --countPath replicates/{wildcards.data}/{wildcards.mode}/dimsum/dimsum_count.tsv \
            --experimentDesignPath replicates/{wildcards.data}/{wildcards.mode}/dimsum/exp_design.tsv \
            --wildtypeSequence $REF \
            --outputPath replicates/{wildcards.data}/{wildcards.mode}/dimsum/ \
            --mutagenesisType random \
            --maxSubstitutions 6 \
            --sequenceType noncoding
        """

rule run_enrich2_replicate:
    input:
        "replicates/{data}/{mode}/enrich2/config.json"
    output:
        "replicates/{data}/{mode}/enrich2/results/tsv/simulation_exp/main_identifiers_scores.tsv",
        "replicates/{data}/{mode}/enrich2/results/tsv/simulation_exp/main_identifiers_scores_pvalues_wt.tsv"
    conda:  
        'enrich2'
    shell:
        """
        enrich_cmd replicates/{wildcards.data}/{wildcards.mode}/enrich2/config.json WLS wt --no-plots
        """

rule run_rosace_replicate: 
    input:
        "replicates/{data}/{mode}/raw/rosace.rda"
    output:
        "replicates/{data}/{mode}/rosace/rosace_eval.rda"
    log: 
        stdout = "logs/{data}/{mode}/rosace.stdout",
        stderr =  "logs/{data}/{mode}/rosace.stderr"
    threads: 4
    benchmark:
        "replicates/{data}/{mode}/rosace/benchmark.txt"
    shell:
        """
        Rscript workflow/scripts/run_rosace.R \
            --data {wildcards.data} --rep {wildcards.mode} \
            > {log.stdout} 2> {log.stderr}
        """

# Rscript workflow/scripts/run_rosace.R --data OCT1 --rep rep1_1

rule run_rosace_nopos_replicate: 
    input:
        "replicates/{data}/{mode}/raw/rosace.rda"
    output:
        "replicates/{data}/{mode}/rosace_nopos/rosace_eval.rda"
    log: 
        stdout = "logs/{data}/{mode}/rosace_nopos.stdout",
        stderr =  "logs/{data}/{mode}/rosace_nopos.stderr"
    threads: 4
    benchmark:
        "replicates/{data}/{mode}/rosace_nopos/benchmark.txt"
    shell:
        """
        Rscript workflow/scripts/run_rosace.R \
            --data {wildcards.data} --mode N --rep {wildcards.mode} \
            > {log.stdout} 2> {log.stderr}
        """

rule compare_method_replicate:
    input: 
        "replicates/{data}/{mode}/mutscan/score_edgeR.tsv",
        "replicates/{data}/{mode}/mutscan/score_limma.tsv",
        "replicates/{data}/{mode}/rosace/rosace_eval.rda",
        "replicates/{data}/{mode}/dimsum/DiMSum_Project/DiMSum_Project_fitness_replicates.RData",
        "replicates/{data}/{mode}/enrich2/results/tsv/simulation_exp/main_identifiers_scores.tsv",
        "replicates/{data}/{mode}/enrich2/results/tsv/simulation_exp/main_identifiers_scores_pvalues_wt.tsv"
    output:
        "replicates/{data}/{mode}/analysis/rosace_full.rda",
        "replicates/{data}/{mode}/analysis/effects.tsv"
    shell:
        """
        Rscript workflow/scripts/compare_method.R --data {wildcards.data} --mode {wildcards.mode}
        """

# Rscript workflow/scripts/compare_method.R --data MET --mode rep1_1
# Rscript workflow/scripts/compare_method.R --data OCT1 --mode rep1_1

rule analysis_replicate:
    input: 
        "replicates/{data}/{mode}/analysis/effects.tsv",
        "data/sysdata.rda",
        "data/{data}/gnomAD_cleaned.tsv"
    output:
        expand("replicates/{{data}}/{{mode}}/plot/{file}", file = ["stop_corr.png"])
    shell:
        """
        Rscript workflow/scripts/analysis.R --data {wildcards.data} --mode {wildcards.mode}
        """



# Rscript workflow/scripts/analysis.R --data MET --mode rep1_1

# Rscript workflow/scripts/analysis.R --data OCT1 --mode rep1_1
# Rscript workflow/scripts/analysis.R --data OCT1 --mode rep1_2
# Rscript workflow/scripts/analysis.R --data OCT1 --mode rep1_3
# Rscript workflow/scripts/analysis.R --data OCT1 --mode rep2_1
# Rscript workflow/scripts/analysis.R --data OCT1 --mode rep2_2
# Rscript workflow/scripts/analysis.R --data OCT1 --mode rep2_3
# Rscript workflow/scripts/analysis.R --data OCT1 --mode rep3_1

# Rscript workflow/scripts/summarize_replicate.R --data OCT1 
# Rscript workflow/scripts/summarize_replicate.R --data MET