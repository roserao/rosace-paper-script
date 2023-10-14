rule compare_method:
    """
    compute correlation, fdr, and rankfdr for SLR, ENRICH2, ROSACE, mutscan, and DIMSUM
    """
    input:
        "results/sim_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/rosace/rosace_eval.rds",
        "results/sim_{sel}/{mode}_rep{rep}_rd1_{cond}/sim{sim}/dimsum/DiMSum_Project/DiMSum_Project_fitness_replicates.RData",
        "results/sim_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/enrich2/results/tsv/simulation_exp/main_identifiers_scores.tsv",
        "results/sim_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/enrich2/results/tsv/simulation_exp/main_identifiers_scores_pvalues_wt.tsv",
        "results/sim_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/mutscan/score_edgeR.tsv",
        "results/sim_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/mutscan/score_limma.tsv"
    output:
        expand("results/sim_{{sel}}/{{mode}}_rep{{rep}}_rd{{rd}}_{{cond}}/sim{{sim}}/analysis/{output}",
            output = ['corr.tsv', 'fdr.tsv', 'rankfdr.tsv', 'rankfdrsense.tsv'])
    shell:
        """
        Rscript workflow/scripts/compare_method.R \
            --mode {wildcards.mode} \
            --rep {wildcards.rep} \
            --round {wildcards.rd} \
            --cond {wildcards.cond} \
            --sim {wildcards.sim} \
            --flag {wildcards.sel}
        """

# Example: Rscript workflow/scripts/compare_method.R --mode growth --rep 1 --round 1 --cond clean --sim 1 --flag pos
# Example: Rscript workflow/scripts/compare_method.R --mode growth --rep 3 --round 1 --cond clean --sim 3 --flag negxxx
# Example: Rscript workflow/scripts/compare_method.R --mode growth --rep 3 --round 3 --cond clean --sim 1 --flag neg

# rule all:
#     input:
#         expand("results/sim_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/analysis/{output}",
#             sel = v_sel, mode = v_mode, rep = v_nrep, rd = v_nround, 
#             cond = v_cond, sim = range(1, 1 + nsim),
#             output = ['corr.tsv', 'fdr.tsv', 'rankfdr.tsv', 'rankfdrsense.tsv'])

rule summarise_condition:
    """
    summarise all rounds of simulation in one data frame for each condition
    """
    input: 
        expand("results/sim_{{sel}}/{{mode}}_rep{{rep}}_rd{{rd}}_{{cond}}/sim{sim}/analysis/{output}",
            sim = range(1, 1 + nsim),
            output = ['corr.tsv', 'fdr.tsv', 'rankfdr.tsv', 'rankfdrsense.tsv'])
    output:
        expand("results/summary_{{sel}}/{{mode}}_rep{{rep}}_rd{{rd}}_{{cond}}/{output}",
            output = ['corr.tsv', 'fdr.tsv', 'rankfdr.tsv', 'rankfdrsense.tsv', "benchmark.tsv"])
    shell:
        """
        Rscript workflow/scripts/summarise_condition.R \
            --mode {wildcards.mode} \
            --rep {wildcards.rep} \
            --round {wildcards.rd} \
            --cond {wildcards.cond} \
            --flag {wildcards.sel}
        """

# rule all:
#     input:
#         expand("results/summary_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/{output}",
#             sel = v_sel, mode = v_mode, rep = v_nrep, rd = v_nround, cond = v_cond, 
#             output = ['corr.tsv', 'fdr.tsv', 'rankfdr.tsv', 'rankfdrsense.tsv', "benchmark.tsv"])

# Example: Rscript workflow/scripts/summarise_condition.R --mode growth --rep 1 --round 1 --cond clean --flag negxxx

rule summarise_run:
    """
    summarise all simulation needed for one run (pos, neg, posxxx, negxxx)
    """
    input: 
        expand("results/summary_{{sel}}/{{mode}}_rep{rep}_rd{rd}_{{cond}}/{output}",
            rep = v_nrep, rd = v_nround,
            output = ['corr.tsv', 'fdr.tsv', 'rankfdr.tsv', 'rankfdrsense.tsv', "benchmark.tsv"])
    output:
        expand("results/output_{{sel}}/{{mode}}_{{cond}}/{output}",
            output = ['corr.tsv', 'fdr.tsv', 'rankfdr.tsv', 'rankfdrsense.tsv', "benchmark.tsv"])
    shell:
        """
        Rscript workflow/scripts/summarise_run.R \
            --mode {wildcards.mode} \
            --cond {wildcards.cond} \
            --flag {wildcards.sel}
        """

# Example: Rscript workflow/scripts/summarise_run.R --mode growth --cond clean --flag neg

# rule all:
#     input:
#         expand("results/output_{sel}/{mode}_{cond}/{output}",
#             sel = v_sel, mode = v_mode, cond = v_cond,
#             output = ['corr.tsv', 'fdr.tsv', 'rankfdr.tsv', 'rankfdrsense.tsv', "benchmark.tsv"])