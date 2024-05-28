# This is the main file for how data was treated with regards to metabolites

# data analysis (from peak intensity tables), were done using MetaboanalystR on a laptop. (metaboanalystR was a headache getting on the cluster)

# I used the script: "analyzing_metabolome.r".
# overall, I normalized to sum (was done in excel before reading it into R for LCMS data)
# and used pareto scaling
# This was done for both NMR and LCMS data.
# then I did Anova with fdr and exported the results.
