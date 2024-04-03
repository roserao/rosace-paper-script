rule generate_rosette:
    input: 
        "data/{data}/rosace.rda"
    output: 
        "results/rosette/{data}/rosette.rda"
    shell: 
        """
        Rscript workflow/scripts/rosette.R --data {wildcards.data}
        """
# Rscript workflow/scripts/rosette.R --data OCT1
# Rscript workflow/scripts/rosette.R --data MET
# Rscript workflow/scripts/rosette.R --data CARD11

rule generate_simulation:
    input: 
         "results/rosette/{data}/rosette.rda"
    output: 
        expand("results/{{data}}/{{pos}}/growth_rep{{rep}}_rd{{rd}}_clean/sim{sim}/rosace/rosace.rds",
            sim = range(1, 1 + nsim)),
        expand("results/{{data}}/{{pos}}/growth_rep{{rep}}_rd{{rd}}_clean/sim{sim}/enrich2/config.json",
            sim = range(1, 1 + nsim))
    shell: 
        """
        Rscript workflow/scripts/run_rosette.R --data {wildcards.data} --sim {nsim} --rep {wildcards.rep} --round {wildcards.rd} --pos {wildcards.pos}
        """
# Rscript workflow/scripts/run_rosette.R --data CARD11 --sim 1 --rep 3 --round 3 --pos pos

# rule all:
#     input:
#         expand("results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/rosace/rosace.rds",
#             data = v_data, pos = v_pos, rep = v_nrep, rd = v_nround, sim = range(1, 1 + nsim))


rule prepare_method:
    input: 
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/rosace/rosace.rds"
    output:
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/mutscan/score_edgeR.tsv",
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/mutscan/score_limma.tsv",
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/dimsum/dimsum_count.tsv"
    shell:
        """
        Rscript workflow/scripts/prepare_method.R --data {wildcards.data} --sim {wildcards.sim} --rep {wildcards.rep} --round {wildcards.rd} --pos {wildcards.pos}
        """

# Rscript workflow/scripts/prepare_method.R --data MET --sim 1 --rep 3 --round 1 --pos nopos
# Rscript workflow/scripts/prepare_method.R --data CARD11 --sim 1 --rep 3 --round 3 --pos pos

rule run_enrich_wls:
    input:
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/enrich2/config.json"
    output:
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/enrich2/results/tsv/simulation_exp/main_identifiers_scores.tsv",
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/enrich2/results/tsv/simulation_exp/main_identifiers_scores_pvalues_wt.tsv"
    wildcard_constraints:
        rd = "[2-3]"
    conda:  
        'enrich2'
    shell:
        """
        enrich_cmd results/{wildcards.data}/{wildcards.pos}/growth_rep{wildcards.rep}_rd{wildcards.rd}_clean/sim{wildcards.sim}/enrich2/config.json \
            WLS wt --no-plots
        """

# Example: enrich_cmd results/sim_pos/growth_rep2_rd2_clean/sim1/enrich2/config.json WLS wt --no-plots

rule run_enrich_ratios:
    input:
        "results/{data}/{pos}/growth_rep{rep}_rd1_clean/sim{sim}/enrich2/config.json"
    output:
        "results/{data}/{pos}/growth_rep{rep}_rd1_clean/sim{sim}/enrich2/results/tsv/simulation_exp/main_identifiers_scores.tsv",
        "results/{data}/{pos}/growth_rep{rep}_rd1_clean/sim{sim}/enrich2/results/tsv/simulation_exp/main_identifiers_scores_pvalues_wt.tsv"
    conda:  
        'enrich2'
    shell:
        """
        enrich_cmd results/{wildcards.data}/{wildcards.pos}/growth_rep{wildcards.rep}_rd1_clean/sim{wildcards.sim}/enrich2/config.json \
            ratios wt --no-plots
        """

rule run_dimsum:
    input: 
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/dimsum/dimsum_count.tsv"
    output:
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/dimsum/DiMSum_Project/DiMSum_Project_fitness_replicates.RData"
    conda:  
        'dimsum'
    shell:
        """
        module unload R
        REF="$(cat results/{wildcards.data}/{wildcards.pos}/growth_rep{wildcards.rep}_rd{wildcards.rd}_clean/sim{wildcards.sim}/dimsum/dimsum_ref.txt)"
        DiMSum --startStage 4 \
            --countPath results/{wildcards.data}/{wildcards.pos}/growth_rep{wildcards.rep}_rd{wildcards.rd}_clean/sim{wildcards.sim}/dimsum/dimsum_count.tsv \
            --experimentDesignPath results/{wildcards.data}/{wildcards.pos}/growth_rep{wildcards.rep}_rd{wildcards.rd}_clean/sim{wildcards.sim}/dimsum/exp_design.tsv \
            --wildtypeSequence $REF \
            --outputPath results/{wildcards.data}/{wildcards.pos}/growth_rep{wildcards.rep}_rd{wildcards.rd}_clean/sim{wildcards.sim}/dimsum/ \
            --mutagenesisType random \
            --maxSubstitutions 6 \
            --sequenceType noncoding
        """

rule run_rosace:
    input:
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/rosace/rosace.rds"
    output:
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/rosace/rosace_eval.rds",
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/rosace/benchmark.txt"
    log: 
        stdout = "logs/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/rosace.stdout",
        stderr = "logs/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/rosace.stderr"
    threads: 4
    benchmark:
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/rosace/benchmark.txt"
    shell:
        """
        Rscript workflow/scripts/run_rosace.R --data {wildcards.data} --sim {wildcards.sim} --rep {wildcards.rep} --round {wildcards.rd} --pos {wildcards.pos} --mode P > {log.stdout} 2> {log.stderr}
        """

rule run_rosace_nopos:
    input:
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/rosace/rosace.rds"
    output:
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/rosace_nopos/rosace_eval.rds"
    log: 
        stdout = "logs/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/rosace_nopos.stdout",
        stderr = "logs/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/rosace_nopos.stderr"
    threads: 4
    benchmark:
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/rosace_nopos/benchmark.txt"
    shell:
        """
        Rscript workflow/scripts/run_rosace.R --data {wildcards.data} --sim {wildcards.sim} --rep {wildcards.rep} --round {wildcards.rd} --pos {wildcards.pos} --mode N > {log.stdout} 2> {log.stderr}
        """

# Rscript workflow/scripts/run_rosace.R --data MET --sim 1 --rep 1 --round 3 --pos nopos --mode N
# Rscript workflow/scripts/run_rosace.R --data MET --sim 1 --rep 1 --round 3 --pos nopos --mode P

rule compare_method:
    input:
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/rosace_nopos/rosace_eval.rds",
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/rosace/rosace_eval.rds",
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/dimsum/DiMSum_Project/DiMSum_Project_fitness_replicates.RData",
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/enrich2/results/tsv/simulation_exp/main_identifiers_scores.tsv",
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/enrich2/results/tsv/simulation_exp/main_identifiers_scores_pvalues_wt.tsv",
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/mutscan/score_edgeR.tsv",
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/mutscan/score_limma.tsv"
    output:
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/analysis/effects.tsv"
    shell:
        """
        Rscript workflow/scripts/compare_method.R --data {wildcards.data} --sim {wildcards.sim} --rep {wildcards.rep} --round {wildcards.rd} --pos {wildcards.pos}
        """

# Rscript workflow/scripts/compare_method.R --data MET --sim 1 --rep 1 --round 3 --pos nopos 


rule analysis_sim:
    input:
        "results/{data}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/analysis/effects.tsv"
    output:
        expand("results/{{data}}/{{pos}}/growth_rep{{rep}}_rd{{rd}}_clean/sim{{sim}}/analysis/{file}",
            file = ["corr.tsv", "fdr.tsv", "rankfdr.tsv", "rankfdrsense.tsv"])
    shell:
        """
        Rscript workflow/scripts/analysis_sim.R --data {wildcards.data} --sim {wildcards.sim} --rep {wildcards.rep} --round {wildcards.rd} --pos {wildcards.pos}
        """

# Rscript workflow/scripts/analysis_sim.R --data MET --sim 1 --rep 1 --round 3 --pos nopos 

rule analysis_gene:
    input:
        expand("results/{{data}}/{pos}/growth_rep{rep}_rd{rd}_clean/sim{sim}/analysis/{file}",
            pos = v_pos, rep = v_nrep, rd = v_nround, sim = range(1, 1 + nsim),
            file = ["corr.tsv", "fdr.tsv", "rankfdr.tsv", "rankfdrsense.tsv"])
    output:
        expand("results/{{data}}/summary/{file}",
            file = ["corr.tsv", "fdr.tsv", "rankfdr.tsv", "rankfdrsense.tsv", "benchmark.tsv"])
    shell:
        """
        Rscript workflow/scripts/analysis_gene.R --data {wildcards.data}
        """

rule plot_gene:
    input:
        expand("results/{{data}}/summary/{file}",
            file = ["corr.tsv", "fdr.tsv", "rankfdr.tsv", "rankfdrsense.tsv", "benchmark.tsv"])
    output:
        expand("results/{{data}}/plot/{file}",
            file = ["corr.png", "fdr.png", "rankfdr.png", "time.png"])
    shell:
        """
        Rscript workflow/scripts/plot_gene.R --data {wildcards.data}
        """

# Rscript workflow/scripts/analysis_gene.R --data MET 
# Rscript workflow/scripts/plot_gene.R --data MET 
# Rscript workflow/scripts/plot_gene.R --data OCT1
# Rscript workflow/scripts/plot_gene.R --data CARD11

