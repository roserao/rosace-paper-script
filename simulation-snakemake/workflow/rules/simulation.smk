rule generate_simulation_pos:
    """
    'mode', 'rep', and 'rd' are wildcards.
    'sim' should be range(n_sim).
    'cond' should be 'clean' and 'reperror'.
    """
    input: 
        "data/rosette_1SM73.RData"
    output: 
        expand("results/sim_pos/{{mode}}_rep{{rep}}_rd{{rd}}_{cond}/sim{sim}/rosace/rosace.rds",
            cond = v_cond, sim = range(1, 1 + nsim)),
        expand("results/sim_pos/{{mode}}_rep{{rep}}_rd{{rd}}_{cond}/sim{sim}/enrich2/config.json",
            cond = v_cond, sim = range(1, 1 + nsim))
    shell: 
        """
        Rscript workflow/scripts/generate_simulation.R \
            --mode {wildcards.mode} \
            --sim {nsim} \
            --rep {wildcards.rep} \
            --round {wildcards.rd} \
            --flag pos \
            --pos T
        """

rule generate_simulation_posxxx:
    input: 
        "data/rosette_1SM73.RData"
    output: 
        expand("results/sim_posxxx/{{mode}}_rep{{rep}}_rd{{rd}}_{cond}/sim{sim}/rosace/rosace.rds",
            cond = v_cond, sim = range(1, 1 + nsim)),
        expand("results/sim_posxxx/{{mode}}_rep{{rep}}_rd{{rd}}_{cond}/sim{sim}/enrich2/config.json",
            cond = v_cond, sim = range(1, 1 + nsim))
    shell: 
        """
        Rscript workflow/scripts/generate_simulation.R \
            --mode {wildcards.mode} \
            --sim {nsim} \
            --rep {wildcards.rep} \
            --round {wildcards.rd} \
            --flag pos \
            --pos F
        """

rule generate_simulation_neg:
    input: 
        "data/rosette_MET.RData"
    output: 
        expand("results/sim_neg/{{mode}}_rep{{rep}}_rd{{rd}}_{cond}/sim{sim}/rosace/rosace.rds",
            cond = v_cond, sim = range(1, 1 + nsim)),
        expand("results/sim_neg/{{mode}}_rep{{rep}}_rd{{rd}}_{cond}/sim{sim}/enrich2/config.json",
            cond = v_cond, sim = range(1, 1 + nsim))
    shell: 
        """
        Rscript workflow/scripts/generate_simulation.R \
            --mode {wildcards.mode} \
            --sim {nsim} \
            --rep {wildcards.rep} \
            --round {wildcards.rd} \
            --flag neg \
            --pos T
        """

rule generate_simulation_negxxx:
    input: 
        "data/rosette_MET.RData"
    output: 
        expand("results/sim_negxxx/{{mode}}_rep{{rep}}_rd{{rd}}_{cond}/sim{sim}/rosace/rosace.rds",
            cond = v_cond, sim = range(1, 1 + nsim)),
        expand("results/sim_negxxx/{{mode}}_rep{{rep}}_rd{{rd}}_{cond}/sim{sim}/enrich2/config.json",
            cond = v_cond, sim = range(1, 1 + nsim))
    shell: 
        """
        Rscript workflow/scripts/generate_simulation.R \
            --mode {wildcards.mode} \
            --sim {nsim} \
            --rep {wildcards.rep} \
            --round {wildcards.rd} \
            --flag neg \
            --pos F
        """



# rule all:
#     input:
#         expand("results/sim_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/rosace/rosace.rds",
#             sel = v_sel, mode = v_mode, rep = v_nrep, rd = v_nround, 
#             cond = v_cond, sim = range(1, 1 + nsim)),
#         expand("results/sim_{sel}/{mode}_rep{rep}_rd{rd}_{cond}/sim{sim}/enrich2/config.json",
#             sel = v_sel, mode = v_mode, rep = v_nrep, rd = v_nround, 
#             cond = v_cond, sim = range(1, 1 + nsim))

