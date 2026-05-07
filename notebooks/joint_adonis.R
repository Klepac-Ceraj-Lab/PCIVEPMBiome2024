#####
# Code for Figure 3 - Community-level descriptive analysis
# Marginal effects in `adonis2` PERMANOVAs
#####

library(magrittr)
library(dplyr)
library(vegan)
library(here)
library(data.table)
library(tibble)

extract_term <- function(fit, term, analysis) {

  tab <- as.data.frame(fit) %>%
    rownames_to_column("term")

  tab %>%
    filter(term == !!term) %>%
    transmute(
      analysis = analysis,
      term = term,
      R2 = R2,
      pval = `Pr(>F)`
    )
}

functional_mdata_cols = c(
    "subject_id",
    "sample",
    "datasource",
    "child_sex",
    "pci_assess_age",
    "mbiome_sample_age",
    "maternal_education",
    "feeding_state",
    "delivery_mode",
    "MaternalEntropy",
    "InfantVisAtt"
)

taxonomic_mdata_cols = c(
    functional_mdata_cols,
    "alpha_shannon"
)

taxonomic_data <- read.csv(
    fs::path(
        here::here(),
        "data",
        "adonis",
        "taxonmic_inputs.csv"
    )
)

functional_data <- fread(
    fs::path(
        here::here(),
        "data",
        "adonis",
        "functional_inputs.csv"
    ),
    sep = ","
)

taxonomic_metadata <- select(taxonomic_data, all_of(taxonomic_mdata_cols))
taxonomic_abundances <- select(taxonomic_data, !all_of(taxonomic_mdata_cols))

functional_metadata <- select(functional_data, all_of(functional_mdata_cols))
functional_abundances <- select(functional_data, !all_of(functional_mdata_cols))

set.seed(42)
D_taxa <- vegdist(as.matrix(taxonomic_abundances), method = "bray")
D_unirefs <- vegdist(as.matrix(functional_abundances), method = "bray")

## Fits for Taxa
# Separate models
fit_visual_only_taxa <- adonis2(D_taxa ~ InfantVisAtt,data = taxonomic_metadata,permutations = 999)
fit_feed_only_taxa <- adonis2(D_taxa ~ feeding_state,data = taxonomic_metadata,permutations = 999)
# Joint model
fit_joint_taxa <- adonis2(D_taxa ~ feeding_state + InfantVisAtt,data = taxonomic_metadata,permutations = 999)
# Marginal (conditional) effects
fit_margin_taxa <- adonis2(D_taxa ~ feeding_state + InfantVisAtt,data = taxonomic_metadata,permutations = 999,by = "margin")

adonis_summary_taxa <- bind_rows(
  extract_term(fit_visual_only_taxa,"Model","visual_only"),     # Unique effect for visual
  extract_term(fit_feed_only_taxa,"Model","feeding_only"),     # Unique effect for feed
  extract_term(fit_joint_taxa,"Model","joint_model"),                  # Joint effect
  extract_term(fit_margin_taxa,"InfantVisAtt","visual_after_feeding"), # Marginal/conditional effect for visual after feeding
  extract_term(fit_margin_taxa,"feeding_state","feeding_after_visual") # Marginal/conditional effect for feeding after visual
) %>% mutate(profile = "taxa", .before = 1)
adonis_summary_taxa

## Fits for Unirefs
# Separate models
fit_visual_only_unirefs <- adonis2(D_unirefs ~ InfantVisAtt,data = functional_metadata,permutations = 999)
fit_feed_only_unirefs <- adonis2(D_unirefs ~ feeding_state,data = functional_metadata,permutations = 999)
# Joint model
fit_joint_unirefs <- adonis2(D_unirefs ~ feeding_state + InfantVisAtt,data = functional_metadata,permutations = 999)
# Marginal (conditional) effects
fit_margin_unirefs <- adonis2(D_unirefs ~ feeding_state + InfantVisAtt,data = functional_metadata,permutations = 999,by = "margin")

adonis_summary_unirefs <- bind_rows(
  extract_term(fit_visual_only_unirefs,"Model","visual_only"),     # Unique effect for visual
  extract_term(fit_feed_only_unirefs,"Model","feeding_only"),     # Unique effect for feed
  extract_term(fit_joint_unirefs,"Model","joint_model"),                  # Joint effect
  extract_term(fit_margin_unirefs,"InfantVisAtt","visual_after_feeding"), # Marginal/conditional effect for visual after feeding
  extract_term(fit_margin_unirefs,"feeding_state","feeding_after_visual") # Marginal/conditional effect for feeding after visual
) %>% mutate(profile = "unirefs", .before = 1)
adonis_summary_unirefs

## Variance partition
varpart_taxa <- varpart(
  D_taxa,
  ~ InfantVisAtt,
  ~ feeding_state,
  data = taxonomic_metadata
)

dbrda_feed_taxa <- dbrda(
  D_taxa ~ feeding_state + Condition(InfantVisAtt),
  data = taxonomic_metadata
)
anova(dbrda_feed_taxa, permutations = 999)

dbrda_vob_taxa <- dbrda(
  D_taxa ~ InfantVisAtt + Condition(feeding_state),
  data = taxonomic_metadata
)
anova(dbrda_vob_taxa, permutations = 999)

## Variance partition
varpart_unirefs <- varpart(
  D_unirefs,
  ~ InfantVisAtt,
  ~ feeding_state,
  data = functional_metadata
)

dbrda_feed_unirefs <- dbrda(
  D_unirefs ~ feeding_state + Condition(InfantVisAtt),
  data = functional_metadata
)
anova(dbrda_feed_unirefs, permutations = 999)

dbrda_vob_unirefs <- dbrda(
  D_unirefs ~ InfantVisAtt + Condition(feeding_state),
  data = functional_metadata
)
anova(dbrda_vob_unirefs, permutations = 999)

sig_code <- function(p) {
  case_when(
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    p <= 0.1   ~ ".",
    TRUE       ~ ""
  )
}

write.csv(
    x = mutate(
        bind_rows(adonis_summary_taxa, adonis_summary_unirefs),
        sig_code = sig_code(pval)
    ),
    file = fs::path(
        here::here(),
        "manuscript",
        "tables",
        "adonis_results.csv"
    ),
    row.names = FALSE
)

## Printing results

varpart_taxa
anova(dbrda_feed_taxa, permutations = 999)
anova(dbrda_vob_taxa, permutations = 999)

varpart_unirefs
anova(dbrda_feed_unirefs, permutations = 999)
anova(dbrda_vob_unirefs, permutations = 999)