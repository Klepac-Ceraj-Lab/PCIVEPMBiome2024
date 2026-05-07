sens_entropy = perform_fsea_prevalence_sensitivity(
    joinpath(Base.pwd(), "manuscript", "gene_glm", "lms_Model3_entropy.csv"),
    joinpath(Base.pwd(), "manuscript", "gene_glm", "fsea_sensitivity_consolidated_Model3_entropy.csv");
    prevalence_col = :Prevalence,
    cutoffs = 0.01:0.01:0.20,
    should_consolidate = true
)

sens_visual = perform_fsea_prevalence_sensitivity(
    joinpath(Base.pwd(), "manuscript", "gene_glm", "lms_Model3_visual.csv"),
    joinpath(Base.pwd(), "manuscript", "gene_glm", "fsea_sensitivity_consolidated_Model3_visual.csv");
    prevalence_col = :Prevalence,
    cutoffs = 0.01:0.01:0.20,
    should_consolidate = true
)

sens_product = perform_fsea_prevalence_sensitivity(
    joinpath(Base.pwd(), "manuscript", "gene_glm", "lms_Model3_product.csv"),
    joinpath(Base.pwd(), "manuscript", "gene_glm", "fsea_sensitivity_consolidated_Model3_product.csv");
    prevalence_col = :Prevalence,
    cutoffs = 0.01:0.01:0.20,
    should_consolidate = true
)

using CairoMakie

plot_fsea_sensitivity(sens_visual; genesets = [
    "Acetate synthesis",
    "Butyrate synthesis",
    "ClpB",
    "GABA synthesis",
    "Glutamate degradation",
    "Glutamate synthesis",
    "Menaquinone synthesis",
    "Propionate degradation",
    "Propionate synthesis",
    "Quinolinic acid degradation",
    "Quinolinic acid synthesis",
    "S-Adenosylmethionine synthesis",
    "Tryptophan synthesis",
    "p-Cresol synthesis"
], type = :es)