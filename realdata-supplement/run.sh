
Rscript workflow/scripts/prepare_method.R --key CARD11
Rscript workflow/scripts/prepare_method.R --key MSH2
Rscript workflow/scripts/prepare_method.R --key BRCA1

conda activate enrich2
enrich_cmd results/CARD11/enrich2/config.json ratios wt --no-plots
enrich_cmd results/MSH2/enrich2/config.json ratios wt --no-plots
enrich_cmd results/BRCA1/enrich2/config.json WLS wt --no-plots
conda deactivate

conda activate dimsum
module unload R
REF="$(cat results/CARD11/dimsum/dimsum_ref.txt)"
DiMSum --startStage 4 \
    --countPath results/CARD11/dimsum/dimsum_count.tsv \
    --experimentDesignPath results/CARD11/dimsum/exp_design.tsv \
    --wildtypeSequence $REF \
    --outputPath results/CARD11/dimsum/ \
    --mutagenesisType random \
    --maxSubstitutions 6 \
    --sequenceType noncoding

REF="$(cat results/MSH2/dimsum/dimsum_ref.txt)"
DiMSum --startStage 4 \
    --countPath results/MSH2/dimsum/dimsum_count.tsv \
    --experimentDesignPath results/MSH2/dimsum/exp_design.tsv \
    --wildtypeSequence $REF \
    --outputPath results/MSH2/dimsum/ \
    --mutagenesisType random \
    --maxSubstitutions 6 \
    --sequenceType noncoding
conda deactivate

Rscript workflow/scripts/compare_method.R --key CARD11
Rscript workflow/scripts/compare_method.R --key MSH2
Rscript workflow/scripts/compare_method.R --key BRCA1