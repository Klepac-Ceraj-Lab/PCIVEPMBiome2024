#####
# FSEA codeblock
#####

## Model 1
perform_fsea( "manuscript/FSEA/lms_Model1_child_sex.csv", "manuscript/FSEA/fsea_consolidated_Model1_child_sex.csv"; should_consolidate=true)
perform_fsea( "manuscript/FSEA/lms_Model1_mbiome_age.csv", "manuscript/FSEA/fsea_consolidated_Model1_mbiome_age.csv"; should_consolidate=true)

## Model 2
perform_fsea( "manuscript/FSEA/lms_Model2_child_sex.csv", "manuscript/FSEA/fsea_consolidated_Model2_child_sex.csv"; should_consolidate=true)
perform_fsea( "manuscript/FSEA/lms_Model2_mbiome_age.csv", "manuscript/FSEA/fsea_consolidated_Model2_mbiome_age.csv"; should_consolidate=true)
perform_fsea( "manuscript/FSEA/lms_Model2_visual.csv", "manuscript/FSEA/fsea_consolidated_Model2_visual.csv"; should_consolidate=true)
perform_fsea( "manuscript/FSEA/lms_Model2_entropy.csv", "manuscript/FSEA/fsea_consolidated_Model2_entropy.csv"; should_consolidate=true)

## Model 2.1 (Entropy)
perform_fsea( "manuscript/FSEA/lms_Model21_child_sex.csv", "manuscript/FSEA/fsea_consolidated_Model21_child_sex.csv"; should_consolidate=true)
perform_fsea( "manuscript/FSEA/lms_Model21_mbiome_age.csv", "manuscript/FSEA/fsea_consolidated_Model21_mbiome_age.csv"; should_consolidate=true)
perform_fsea( "manuscript/FSEA/lms_Model21_entropy.csv", "manuscript/FSEA/fsea_consolidated_Model21_entropy.csv"; should_consolidate=true)

## Model 2.2 (Visual)
perform_fsea( "manuscript/FSEA/lms_Model22_child_sex.csv", "manuscript/FSEA/fsea_consolidated_Model22_child_sex.csv"; should_consolidate=true)
perform_fsea( "manuscript/FSEA/lms_Model22_mbiome_age.csv", "manuscript/FSEA/fsea_consolidated_Model22_mbiome_age.csv"; should_consolidate=true)
perform_fsea( "manuscript/FSEA/lms_Model22_visual.csv", "manuscript/FSEA/fsea_consolidated_Model22_visual.csv"; should_consolidate=true)

## Model 3
perform_fsea( "manuscript/FSEA/lms_Model3_child_sex.csv", "manuscript/FSEA/fsea_consolidated_Model3_child_sex.csv"; should_consolidate=true)
perform_fsea( "manuscript/FSEA/lms_Model3_visual.csv", "manuscript/FSEA/fsea_consolidated_Model3_visual.csv"; should_consolidate=true)
perform_fsea( "manuscript/FSEA/lms_Model3_entropy.csv", "manuscript/FSEA/fsea_consolidated_Model3_entropy.csv"; should_consolidate=true)
perform_fsea( "manuscript/FSEA/lms_Model3_product.csv", "manuscript/FSEA/fsea_consolidated_Model3_product.csv"; should_consolidate=true)