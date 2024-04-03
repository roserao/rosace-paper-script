rule run_simdms:
    input: 
        "workflow/simdms/config.json"
    output: 
        expand("results/simdms/{data}/{mode}/rosace/rosace.rds", data = simdms_data, mode = simdms_mode),
        expand("results/simdms/{data}/{mode}/tsv/counts.tsv", data = simdms_data, mode = simdms_mode)
    conda:
        "enrich2"
    shell: 
        """
        python workflow/simdms/simdms.py workflow/simdms/config.json
        python workflow/scripts/simdms_h5.py 
        enrich_cmd results/simdms/binding_simulation_depth_200/config.json WLS wt --no-plots
        enrich_cmd results/simdms/growth_simulation_depth_200/config.json WLS wt --no-plots
        Rscript workflow/scripts/create_rosace_simdms.R 
        """

rule prepare_simdms:
    input: 
        "results/simdms/{data}/{mode}/rosace/rosace.rds",
        "results/simdms/{data}/{mode}/tsv/counts.tsv"
    output:
        "results/simdms/{data}/{mode}/mutscan/score_edgeR.tsv",
        "results/simdms/{data}/{mode}/mutscan/score_limma.tsv",
        "results/simdms/{data}/{mode}/dimsum/dimsum_count.tsv"
    shell:
        """
        Rscript workflow/scripts/prepare_method_simdms.R --data {wildcards.data} --mode {wildcards.mode}
        """

# Rscript workflow/scripts/prepare_method_simdms.R --data binding --mode clean 

rule run_dimsum_simdms:
    input: 
        "results/simdms/{data}/{mode}/dimsum/dimsum_count.tsv"
    output:
        "results/simdms/{data}/{mode}/dimsum/DiMSum_Project/DiMSum_Project_fitness_replicates.RData"
    conda:  
        'dimsum'
    shell:
        """
        module unload R
        REF="$(cat results/simdms/{wildcards.data}/{wildcards.mode}/dimsum/dimsum_ref.txt)"
        DiMSum --startStage 4 \
            --countPath results/simdms/{wildcards.data}/{wildcards.mode}/dimsum/dimsum_count.tsv \
            --experimentDesignPath results/simdms/{wildcards.data}/{wildcards.mode}/dimsum/exp_design.tsv \
            --wildtypeSequence $REF \
            --outputPath results/simdms/{wildcards.data}/{wildcards.mode}/dimsum/ \
            --mutagenesisType random \
            --maxSubstitutions 6 \
            --sequenceType noncoding
        """

rule run_rosace_simdms:
    input:
        "results/simdms/{data}/{mode}/rosace/rosace.rds"
    output:
        "results/simdms/{data}/{mode}/rosace/rosace_eval.rds"
    log: 
        stdout = "logs/{data}/{mode}/rosace.stdout",
        stderr = "logs/{data}/{mode}/rosace.stderr"
    threads: 4
    benchmark:
        "results/simdms/{data}/{mode}/rosace/benchmark.txt"
    shell:
        """
        Rscript workflow/scripts/run_rosace_simdms.R --data {wildcards.data} --mode {wildcards.mode} > {log.stdout} 2> {log.stderr}
        """

rule compare_method_simdms:
    input: 
        "results/simdms/{data}/{mode}/rosace/rosace_eval.rds",
        "results/simdms/{data}/{mode}/dimsum/DiMSum_Project/DiMSum_Project_fitness_replicates.RData",
        "results/simdms/{data}/{mode}/mutscan/score_edgeR.tsv",
        "results/simdms/{data}/{mode}/mutscan/score_limma.tsv"
    output:
        "results/simdms/{data}/{mode}/analysis/effects.tsv"
    shell:
        """
        Rscript workflow/scripts/compare_method_simdms.R --data {wildcards.data} --mode {wildcards.mode} 
        """

# Rscript workflow/scripts/compare_method_simdms.R --data growth --mode clean
# Rscript workflow/scripts/compare_method_simdms.R --data growth --mode jackpot
# Rscript workflow/scripts/compare_method_simdms.R --data growth --mode reperror
# Rscript workflow/scripts/compare_method_simdms.R --data binding --mode clean
# Rscript workflow/scripts/compare_method_simdms.R --data binding --mode jackpot
# Rscript workflow/scripts/compare_method_simdms.R --data binding --mode reperror


rule analysis_simdms:
    input: 
        "results/simdms/{data}/{mode}/analysis/effects.tsv"
    output:
        "results/simdms/{data}/{mode}/plot/corr.png"
    shell:
        """
        Rscript workflow/scripts/analysis_simdms.R --data {wildcards.data} --mode {wildcards.mode} 
        """

# Rscript workflow/scripts/analysis_simdms.R --data growth --mode clean
# Rscript workflow/scripts/analysis_simdms.R --data growth --mode jackpot
# Rscript workflow/scripts/analysis_simdms.R --data growth --mode reperror
# Rscript workflow/scripts/analysis_simdms.R --data binding --mode clean
# Rscript workflow/scripts/analysis_simdms.R --data binding --mode jackpot
# Rscript workflow/scripts/analysis_simdms.R --data binding --mode reperror

