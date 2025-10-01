## 1. Load data:

input_df <- read.csv("manuscript/tables/mdata_table_in_bb3.csv")

input_df$MaternalEntropy <- scale(input_df$MaternalEntropy)
input_df$InfantVisAtt <- scale(input_df$InfantVisAtt)

entropy_cutoff <- 0.0
visatt_cutoff <- 0.0

input_df$Quadrant <- with(input_df, ifelse(
  MaternalEntropy > entropy_cutoff & InfantVisAtt > visatt_cutoff, "Q1",
  ifelse(MaternalEntropy <= entropy_cutoff & InfantVisAtt > visatt_cutoff, "Q2",
  ifelse(MaternalEntropy <= entropy_cutoff & InfantVisAtt <= visatt_cutoff, "Q3", "Q4")
)))

input_df$Quadrant <- factor(input_df$Quadrant, levels = c("Q1", "Q2", "Q3", "Q4"))

input_df$Q1_vs_others <- ifelse(input_df$Quadrant == "Q1", "Q1", "Other")
input_df$Q2_vs_others <- ifelse(input_df$Quadrant == "Q2", "Q2", "Other")
input_df$Q3_vs_others <- ifelse(input_df$Quadrant == "Q3", "Q3", "Other")
input_df$Q4_vs_others <- ifelse(input_df$Quadrant == "Q4", "Q4", "Other")

t.test(Bifidobacterium_breve[Bifidobacterium_breve > 0] ~ Q1_vs_others[Bifidobacterium_breve > 0], data = input_df)
t.test(Bifidobacterium_breve[Bifidobacterium_breve > 0] ~ Q2_vs_others[Bifidobacterium_breve > 0], data = input_df)
t.test(Bifidobacterium_breve[Bifidobacterium_breve > 0] ~ Q3_vs_others[Bifidobacterium_breve > 0], data = input_df)
t.test(Bifidobacterium_breve[Bifidobacterium_breve > 0] ~ Q4_vs_others[Bifidobacterium_breve > 0], data = input_df)

t.test(Bifidobacterium_longum[Bifidobacterium_longum > 0] ~ Q1_vs_others[Bifidobacterium_longum > 0], data = input_df)
t.test(Bifidobacterium_longum[Bifidobacterium_longum > 0] ~ Q2_vs_others[Bifidobacterium_longum > 0], data = input_df)
t.test(Bifidobacterium_longum[Bifidobacterium_longum > 0] ~ Q3_vs_others[Bifidobacterium_longum > 0], data = input_df)
t.test(Bifidobacterium_longum[Bifidobacterium_longum > 0] ~ Q4_vs_others[Bifidobacterium_longum > 0], data = input_df)

mean(input_df$Bifidobacterium_breve[input_df$Quadrant == "Q1"])
mean(input_df$Bifidobacterium_breve[input_df$Quadrant != "Q1"])

t.test(Bifidobacterium_breve[Bifidobacterium_breve > 0] ~ Q1_vs_others[Bifidobacterium_breve > 0], data = )

#####
# Data-driven approach
#####

optimize_quadrant_split <- function(data, A, B, lower = 0.3, upper = 0.7, step = 0.01) {
  # Generate grid of candidate quantile thresholds
  quantile_grid <- expand.grid(
    A_q = seq(lower, upper, by = step),
    B_q = seq(lower, upper, by = step)
  )

  imbalance_scores <- numeric(nrow(quantile_grid))

  for (i in seq_len(nrow(quantile_grid))) {
    A_cut <- quantile(data[[A]], quantile_grid$A_q[i], na.rm = TRUE)
    B_cut <- quantile(data[[B]], quantile_grid$B_q[i], na.rm = TRUE)

    quadrant <- with(data, ifelse(
      data[[A]] > A_cut & data[[B]] > B_cut, "Q1",
      ifelse(data[[A]] <= A_cut & data[[B]] > B_cut, "Q2",
      ifelse(data[[A]] <= A_cut & data[[B]] <= B_cut, "Q3", "Q4")))
    )

    sizes <- table(quadrant)[c("Q1", "Q3", "Q4")]
    imbalance_scores[i] <- var(as.numeric(sizes))  # variance of group sizes
  }

  best_idx <- which.min(imbalance_scores)
  best_split <- quantile_grid[best_idx, ]

  return(list(
    A_quantile = best_split$A_q,
    B_quantile = best_split$B_q,
    A_cut = quantile(data[[A]], best_split$A_q, na.rm = TRUE),
    B_cut = quantile(data[[B]], best_split$B_q, na.rm = TRUE),
    imbalance_score = imbalance_scores[best_idx]
  ))
}

best_split <- optimize_quadrant_split(input_df, "MaternalEntropy", "InfantVisAtt")
A_cut <- best_split$A_cut
B_cut <- best_split$B_cut


input_df$Quadrant_opt <- with(input_df, ifelse(
  MaternalEntropy > A_cut & InfantVisAtt > B_cut, "Q1",
  ifelse(MaternalEntropy <= A_cut & InfantVisAtt > B_cut, "Q2",
  ifelse(MaternalEntropy <= A_cut & InfantVisAtt <= B_cut, "Q3", "Q4")))
)

ggplot(input_df, aes(x = MaternalEntropy, y = InfantVisAtt)) +
  geom_point(aes(color = Quadrant_opt)) +
  geom_vline(xintercept = A_cut, linetype = "dashed") +
  geom_hline(yintercept = B_cut, linetype = "dashed") +
  theme_minimal()

input_df$Q1_vs_others <- ifelse(input_df$Quadrant_opt == "Q1", "Q1", "Other")
input_df$Q2_vs_others <- ifelse(input_df$Quadrant_opt == "Q2", "Q2", "Other")
input_df$Q3_vs_others <- ifelse(input_df$Quadrant_opt == "Q3", "Q3", "Other")
input_df$Q4_vs_others <- ifelse(input_df$Quadrant_opt == "Q4", "Q4", "Other")

anova(Bifidobacterium_breve ~ Q1_vs_others, data = input_df)
wilcox.test(Bifidobacterium_breve ~ Q2_vs_others, data = input_df)
wilcox.test(Bifidobacterium_breve ~ Q3_vs_others, data = input_df)
wilcox.test(Bifidobacterium_breve ~ Q4_vs_others, data = input_df)

kruskal.test(Bifidobacterium_breve ~ Quadrant_opt, data = input_df)
kruskal.test(Bifidobacterium_longum ~ Quadrant_opt, data = input_df)

manova(cbind(Bifidobacterium_breve, Bifidobacterium_longum) ~ Quadrant_opt, data = input_df)

# Example data
X <- input_df[, c("Bifidobacterium_breve", "Bifidobacterium_longum")]
dist_matrix <- dist(X, method = "euclidean")  # or other metric
adonis2(dist_matrix ~ Quadrant_opt, data = input_df, permutations = 999)

assign_bifido_class <- function(bif_duo) {
  if (bif_duo[1] > 0.75) {
    return("dom_breve")
  } else if (bif_duo[2] > 0.75) {
    return("dom_longum")
  } else if ((bif_duo[1] > 0.05) & (((bif_duo[1] + 0.01) / (bif_duo[2] + 0.01)) > 5)) {
    return("more_breve")
  } else if ((bif_duo[2] > 0.05) & (((bif_duo[2] + 0.01) / (bif_duo[1] + 0.01)) > 5)) {
    return("more_longum")
  } else {
    return("null")
  }
}

assign_bifido_class_df <- function(input_df) {
  apply(
    cbind(
      input_df$Bifidobacterium_breve,
      input_df$Bifidobacterium_longum
    ),
    1,
    assign_bifido_class
  )

}

input_df$bifido_class <- assign_bifido_class_df(input_df)

fit <- manova(cbind(InfantVisAtt, MaternalEntropy) ~ bifido_class, data = input_df)
summary(fit)
summary.aov(fit)