#####
# FSEA codeblock
#####

## Model 1
perform_fsea(
    joinpath(Base.pwd(), "manuscript", "gene_glm_previousmodel", "lms_Model1_child_sex.csv"),
    joinpath(Base.pwd(), "manuscript", "FSEA", "fsea_consolidated_Model1_child_sex.csv");
    should_consolidate=true
)

perform_fsea(
    joinpath(Base.pwd(), "manuscript", "gene_glm_previousmodel", "lms_Model1_mbiome_age.csv"),
    joinpath(Base.pwd(), "manuscript", "FSEA", "fsea_consolidated_Model1_mbiome_age.csv");
    should_consolidate=true
)

## Model 2
perform_fsea(
    joinpath(Base.pwd(), "manuscript", "gene_glm_previousmodel", "lms_Model2_visual.csv"),
    joinpath(Base.pwd(), "manuscript", "FSEA", "fsea_consolidated_Model2_visual.csv");
    should_consolidate=true
)
perform_fsea(
    joinpath(Base.pwd(), "manuscript", "gene_glm_previousmodel", "lms_Model2_entropy.csv"),
    joinpath(Base.pwd(), "manuscript", "FSEA", "fsea_consolidated_Model2_entropy.csv");
    should_consolidate=true
)

## Model 3
perform_fsea(
    joinpath(Base.pwd(), "manuscript", "gene_glm_previousmodel", "lms_Model3_visual.csv"),
    joinpath(Base.pwd(), "manuscript", "FSEA", "fsea_consolidated_Model3_visual.csv");
    should_consolidate=true
)
perform_fsea(
    joinpath(Base.pwd(), "manuscript", "gene_glm_previousmodel", "lms_Model3_entropy.csv"),
    joinpath(Base.pwd(), "manuscript", "FSEA", "fsea_consolidated_Model3_entropy.csv");
    should_consolidate=true
)
perform_fsea(
    joinpath(Base.pwd(), "manuscript", "gene_glm_previousmodel", "lms_Model3_product.csv"),
    joinpath(Base.pwd(), "manuscript", "FSEA", "fsea_consolidated_Model3_product.csv");
    should_consolidate=true
)