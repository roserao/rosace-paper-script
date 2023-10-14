rule run_enrich_wls:
    """
    run Enrich2 on the case of round 2 or 3 with weighted least square method
    """
    input:
        "results/sim_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/enrich2/config.json"
    output:
        "results/sim_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/enrich2/results/tsv/simulation_exp/main_identifiers_scores.tsv",
        "results/sim_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/enrich2/results/tsv/simulation_exp/main_identifiers_scores_pvalues_wt.tsv"
    wildcard_constraints:
        rd = "[2-3]"
    conda:  
        'enrich2'
    shell:
        """
        enrich_cmd results/sim_{wildcards.sel}/{wildcards.mode}_rep{wildcards.rep}_rd{wildcards.rd}_{wildcards.cond}/sim{wildcards.sim}/enrich2/config.json \
            WLS wt --no-plots
        """

# Example: enrich_cmd results/sim_pos/growth_rep2_rd2_clean/sim1/enrich2/config.json WLS wt --no-plots

rule run_enrich_ratios:
    """
    run Enrich2 on the case of round 1 with ratios
    """
    input:
        "results/sim_{sel}/{mode}_rep{rep}_rd1_{cond}/sim{sim}/enrich2/config.json"
    output:
        "results/sim_{sel}/{mode}_rep{rep}_rd1_{cond}/sim{sim}/enrich2/results/tsv/simulation_exp/main_identifiers_scores.tsv",
        "results/sim_{sel}/{mode}_rep{rep}_rd1_{cond}/sim{sim}/enrich2/results/tsv/simulation_exp/main_identifiers_scores_pvalues_wt.tsv"
    conda:  
        'enrich2'
    shell:
        """
        enrich_cmd results/sim_{wildcards.sel}/{wildcards.mode}_rep{wildcards.rep}_rd1_{wildcards.cond}/sim{wildcards.sim}/enrich2/config.json \
            ratios wt --no-plots
        """

# Example: enrich_cmd results/sim_pos/growth_rep2_rd1_clean/sim1/enrich2/config.json ratios wt --no-plots

rule run_rosace: 
    input:
        "results/sim_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/rosace/rosace.rds"
    output:
        "results/sim_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/rosace/rosace_eval.rds"
    log: 
        stdout = "results/logs_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/rosace.stdout",
        stderr =  "results/logs_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/rosace.stderr"
    threads: 4
    benchmark:
        "results/sim_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/rosace/benchmark.txt"
    shell:
        """
        Rscript workflow/scripts/run_rosace.R \
            --mode {wildcards.mode} \
            --rep {wildcards.rep} \
            --round {wildcards.rd} \
            --cond {wildcards.cond} \
            --sim {wildcards.sim} \
            --flag {wildcards.sel} \
            > {log.stdout} 2> {log.stderr}
        """

# Example: Rscript workflow/scripts/run_rosace.R --mode growth --rep 1 --round 1 --cond clean --sim 1 --flag pos

# rule all:
#     input:
#         expand("results/sim_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/rosace/rosace_eval.rds",
#             sel = v_sel, mode = v_mode, rep = v_nrep, rd = v_nround,
#             cond = v_cond, sim = range(1, 1 + nsim)),
#         expand("results/sim_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/enrich2/results/tsv/simulation_exp/{output}",
#             sel = v_sel, mode = v_mode, rep = v_nrep, rd = v_nround, 
#             cond = v_cond, sim = range(1, 1 + nsim),
#             output = ["main_identifiers_scores.tsv", "main_identifiers_scores_pvalues_wt.tsv"])

rule run_mutscan: 
    input:
        "results/sim_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/rosace/rosace.rds"
    output:
        "results/sim_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/mutscan/score_edgeR.tsv",
        "results/sim_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/mutscan/score_limma.tsv"
    shell:
        """
        Rscript workflow/scripts/run_mutscan.R \
            --mode {wildcards.mode} \
            --rep {wildcards.rep} \
            --round {wildcards.rd} \
            --cond {wildcards.cond} \
            --sim {wildcards.sim} \
            --flag {wildcards.sel}
        """

# Example: Rscript workflow/scripts/run_mutscan.R --mode growth --rep 1 --round 1 --cond clean --sim 1 --flag pos # fail
# Example: Rscript workflow/scripts/run_mutscan.R --mode growth --rep 3 --round 3 --cond clean --sim 1 --flag pos # success

# rule all:
#     input:
#         expand("results/sim_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/mutscan/{output}",
#             sel = v_sel, mode = v_mode, rep = v_nrep, rd = v_nround, 
#             cond = v_cond, sim = range(1, 1 + nsim),
#             output = ['score_edgeR.tsv', 'score_limma.tsv'])


rule prepare_dimsum:
    """
    prepare DiMSum input on the case of round 1 
    """
    input:
        "results/sim_{sel}/{mode}_rep{rep}_rd1_{cond}/sim{sim}/rosace/rosace.rds"
    output:
        expand("results/sim_{{sel}}/{{mode}}_rep{{rep}}_rd1_{{cond}}/sim{{sim}}/dimsum/{output}",
            output = ["dimsum_count.tsv", "dimsum_ref.txt", "exp_design.tsv"])
    shell:
        """
        Rscript workflow/scripts/prepare_dimsum.R \
            --mode {wildcards.mode} \
            --rep {wildcards.rep} \
            --cond {wildcards.cond} \
            --sim {wildcards.sim} \
            --flag {wildcards.sel}
        """

# Example: Rscript workflow/scripts/prepare_dimsum.R --mode growth --rep 1 --cond clean --sim 1 --flag neg

rule run_dimsum:
    """
    run DiMSum on the case of round 1 
    """
    input: 
        expand("results/sim_{{sel}}/{{mode}}_rep{{rep}}_rd1_{{cond}}/sim{{sim}}/dimsum/{output}",
            output = ["dimsum_count.tsv", "dimsum_ref.txt", "exp_design.tsv"])
    output:
        "results/sim_{sel}/{mode}_rep{rep}_rd1_{cond}/sim{sim}/dimsum/DiMSum_Project/DiMSum_Project_fitness_replicates.RData"
    conda:  
        'dimsum'
    shell:
        """
        module unload R
        REF="$(cat results/sim_{wildcards.sel}/{wildcards.mode}_rep{wildcards.rep}_rd1_{wildcards.cond}/sim{wildcards.sim}/dimsum/dimsum_ref.txt)"
        DiMSum --startStage 4 \
            --countPath results/sim_{wildcards.sel}/{wildcards.mode}_rep{wildcards.rep}_rd1_{wildcards.cond}/sim{wildcards.sim}/dimsum/dimsum_count.tsv \
            --experimentDesignPath results/sim_{wildcards.sel}/{wildcards.mode}_rep{wildcards.rep}_rd1_{wildcards.cond}/sim{wildcards.sim}/dimsum/exp_design.tsv \
            --wildtypeSequence $REF \
            --outputPath results/sim_{wildcards.sel}/{wildcards.mode}_rep{wildcards.rep}_rd1_{wildcards.cond}/sim{wildcards.sim}/dimsum/ \
            --mutagenesisType random \
            --maxSubstitutions 6 \
            --sequenceType noncoding
        """

# Example: 
# module unload R
# REF="$(cat results/sim_neg/growth_rep1_rd1_clean/sim1/dimsum/dimsum_ref.txt)"
# DiMSum --startStage 4 \
#     --countPath results/sim_neg/growth_rep1_rd1_clean/sim1/dimsum/dimsum_count.tsv \
#     --experimentDesignPath results/sim_neg/growth_rep1_rd1_clean/sim1/dimsum/exp_design.tsv \
#     --wildtypeSequence $REF \
#     --outputPath results/sim_neg/growth_rep1_rd1_clean/sim1/dimsum/ \
#     --mutagenesisType random \
#     --maxSubstitutions 6 \
#     --sequenceType noncoding

# rule all:
#     input:
#         expand("results/sim_{sel}/{mode}_rep{rep}_rd1_{cond}/sim{sim}/dimsum/DiMSum_Project/DiMSum_Project_fitness_replicates.RData",
#            sel = v_sel, mode = v_mode, rep = v_nrep, 
#            cond = v_cond, sim = range(1, 1 + nsim))        
