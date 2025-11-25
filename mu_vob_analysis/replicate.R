# Reproduce analyses for MU_VOB.csv in R (SPSS-mode)
# - Type III SS for IC model (car::Anova)
# - Non-robust SEs (no sandwich/HC3) to mirror SPSS
# - sd() in R already uses n-1 (sample SD), matching SPSS for z and ±1SD

library(tidyverse)
library(car)        # for Anova(type=3)

DATA <- "MU_VOB.csv"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) DATA <- args[1]
df <- read.csv(DATA, check.names = FALSE)

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# Column aliases
age <- if ("age_pci_days_NoMissingValues" %in% names(df)) "age_pci_days_NoMissingValues" else "age_pci_days"
epds <- if ("epds_score_3m" %in% names(df)) "epds_score_3m" else if ("epds_score_3m_1" %in% names(df)) "epds_score_3m_1" else NA

df$MU <- to_num(df$MU)
if ("VOB" %in% names(df)) df$VOB <- to_num(df$VOB)
if (!"VOB" %in% names(df) && "VOB_raw" %in% names(df)) df$VOB <- to_num(df$VOB_raw)
df[[age]] <- to_num(df[[age]])
if (!is.na(epds)) df[[epds]] <- to_num(df[[epds]])
df$child_sex <- to_num(df$child_sex)
df$IC <- to_num(df$glitterwands_score_t5)
if ("VOB_For_Age_Residual" %in% names(df)) {
  df$VOBres <- to_num(df$VOB_For_Age_Residual)
} else if ("VOB_resid" %in% names(df)) {
  df$VOBres <- to_num(df$VOB_resid)
} else {
  df$VOBres <- NA_real_
}

dir.create("results_spss_R", showWarnings = FALSE)

# ---- Correlations with series-mean imputation ----
corr_cols <- c("MU","VOB", age)
if (!is.na(epds)) corr_cols <- c(corr_cols, epds)
if ("Maternal_Education" %in% names(df)) {
  df$Maternal_Education <- to_num(df$Maternal_Education)
  corr_cols <- c(corr_cols, "Maternal_Education")
}
corr_df <- df[,corr_cols]
for (nm in names(corr_df)) {
  m <- mean(corr_df[[nm]], na.rm=TRUE)
  corr_df[[nm]][is.na(corr_df[[nm]])] <- m
}
R <- cor(corr_df, use="everything", method="pearson")
P <- matrix(NA_real_, nrow=ncol(corr_df), ncol=ncol(corr_df))
for (i in 1:ncol(corr_df)) for (j in 1:ncol(corr_df)) {
  P[i,j] <- cor.test(corr_df[,i], corr_df[,j])$p.value
}
write.csv(R, file="results_spss_R/correlations_r.csv")
write.csv(P, file="results_spss_R/correlations_p.csv")

# ---- ANCOVA: VOB ~ age + MU (Type II) ----
anc <- na.omit(df[, c("VOB", age, "MU")])
fit <- lm(as.formula(paste("VOB ~", age, "+ MU")), data=anc)
a2 <- Anova(fit, type=2)
ss_res <- tail(a2$`Sum Sq`, 1)
a2$eta_p2 <- a2$`Sum Sq`/(a2$`Sum Sq` + ss_res)
write.csv(as.data.frame(a2), file="results_spss_R/ANCOVA_typeII.csv")

# ---- Age-only + residual correlation ----
age_fit <- lm(as.formula(paste("VOB ~", age)), data=anc)
a_age <- Anova(age_fit, type=2)
write.csv(as.data.frame(a_age), file="results_spss_R/AgeOnly_ANOVA.csv")
resid_vob <- resid(age_fit)
mu_aligned <- anc$MU[as.numeric(names(resid_vob))]
r <- cor(resid_vob, mu_aligned)
p <- cor.test(resid_vob, mu_aligned)$p.value
writeLines(sprintf("r(resid VOB, MU) = %.3f, p = %.3g, df = %d", r, p, length(resid_vob)-2), con="results_spss_R/ResidualCorr.txt")

# ---- 2x2 mean split chi-square ----
if (all(is.na(df$VOBres))) {
  age_fit2 <- lm(as.formula(paste("VOB ~", age)), data=df)
  df$VOBres <- resid(age_fit2)
}
sub <- na.omit(df[, c("MU","VOBres")])
sub$MU_bin <- sub$MU > mean(sub$MU, na.rm=TRUE)
sub$VOB_bin <- sub$VOBres > mean(sub$VOBres, na.rm=TRUE)
ct <- table(sub$MU_bin, sub$VOB_bin)
chis <- chisq.test(ct, correct=FALSE)
phi <- sqrt(chis$statistic / sum(ct))
write.csv(as.data.frame.matrix(ct), file="results_spss_R/ChiSquare_counts.csv")
write.csv(data.frame(Chi2=chis$statistic, df=chis$parameter, p=chis$p.value, phi=phi), file="results_spss_R/ChiSquare_results.csv", row.names=FALSE)

# ---- Inhibitory control model (Type III, non-robust) ----
ic <- df[, c("IC","child_sex", epds, "VOBres", age, "MU")]
ic <- ic[!is.na(ic$IC) & !is.na(ic$VOBres) & !is.na(ic[[age]]) & !is.na(ic$MU), ]
ic$child_sex[is.na(ic$child_sex)] <- mean(ic$child_sex, na.rm=TRUE)
if (!is.na(epds)) ic[[epds]][is.na(ic[[epds]])] <- mean(ic[[epds]], na.rm=TRUE)

form <- as.formula(paste("IC ~ child_sex +", epds, "+ VOBres +", age, "+ MU +", age, ":MU + VOBres:", age, "+ VOBres:MU + VOBres:", age, ":MU"))
fit_ic <- lm(form, data=ic)
a3 <- Anova(fit_ic, type=3)
ss_err3 <- tail(a3$`Sum Sq`, 1)
a3$eta_p2 <- a3$`Sum Sq`/(a3$`Sum Sq` + ss_err3)
write.csv(as.data.frame(a3), file="results_spss_R/IC_TypeIII_ANOVA.csv", row.names=TRUE)

# ---- Simple slopes of MU at ±1 SD of age and VOBres (non-robust) ----
coefs <- coef(fit_ic); nm <- names(coefs)
age_mean <- mean(ic[[age]], na.rm=TRUE); age_sd <- sd(ic[[age]], na.rm=TRUE)       # sample SD
v_mean <- mean(ic$VOBres, na.rm=TRUE); v_sd <- sd(ic$VOBres, na.rm=TRUE)          # sample SD
V <- vcov(fit_ic)  # non-robust covariance

slope <- function(a, v) {
  coef <- coefs["MU"] + coefs[paste(age,":MU",sep="")] * a + coefs["VOBres:MU"] * v + coefs[paste("VOBres:",age,":MU",sep="")] * a * v
  g <- rep(0, length(nm)); names(g) <- nm
  g["MU"] <- 1
  g[paste(age,":MU",sep="")] <- a
  g["VOBres:MU"] <- v
  g[paste("VOBres:",age,":MU",sep="")] <- a*v
  se <- sqrt(as.numeric(t(g) %*% V %*% g))
  t <- coef/se
  p <- 2*pt(abs(t), df=fit_ic$df.residual, lower.tail=FALSE)
  c(coef=coef, se=se, p=p)
}

combos <- expand.grid(age_k=c(-1,1), v_k=c(-1,1))
rows <- list()
for (i in 1:nrow(combos)) {
  a_val <- age_mean + combos$age_k[i]*age_sd
  v_val <- v_mean + combos$v_k[i]*v_sd
  out <- slope(a_val, v_val)
  rows[[i]] <- data.frame(Age.level=ifelse(combos$age_k[i]<0,"-1 SD age","+1 SD age"),
                          VOB.level=ifelse(combos$v_k[i]<0,"-1 SD VOB-for-age","+1 SD VOB-for-age"),
                          b.MU.to.IC=out["coef"], SE=out["se"], p=out["p"])
}
write.csv(bind_rows(rows), file="results_spss_R/IC_simple_slopes_raw_SPSS.csv", row.names=FALSE)

# ---- Standardized DV (sample SD) ----
ic$ICz <- as.numeric(scale(ic$IC, center=TRUE, scale=TRUE))  # R scale uses sample SD
fitz <- lm(update(form, ICz ~ .), data=ic)
Vz <- vcov(fitz)
coefz <- coef(fitz); nmz <- names(coefz)

slopez <- function(a, v) {
  coef <- coefz["MU"] + coefz[paste(age,":MU",sep="")] * a + coefz["VOBres:MU"] * v + coefz[paste("VOBres:",age,":MU",sep="")] * a * v
  g <- rep(0, length(nmz)); names(g) <- nmz
  g["MU"] <- 1; g[paste(age,":MU",sep="")] <- a; g["VOBres:MU"] <- v; g[paste("VOBres:",age,":MU",sep="")] <- a*v
  se <- sqrt(as.numeric(t(g) %*% Vz %*% g))
  t <- coef/se; p <- 2*pt(abs(t), df=fitz$df.residual, lower.tail=FALSE)
  c(coef=coef, se=se, p=p)
}
rows <- list()
for (i in 1:nrow(combos)) {
  a_val <- age_mean + combos$age_k[i]*age_sd
  v_val <- v_mean + combos$v_k[i]*v_sd
  out <- slopez(a_val, v_val)
  rows[[i]] <- data.frame(Age.level=ifelse(combos$age_k[i]<0,"-1 SD age","+1 SD age"),
                          VOB.level=ifelse(combos$v_k[i]<0,"-1 SD VOB-for-age","+1 SD VOB-for-age"),
                          b_z.MU.to.ICz=out["coef"], SEz=out["se"], p=out["p"])
}
write.csv(bind_rows(rows), file="results_spss_R/IC_simple_slopes_standardized_SPSS.csv", row.names=FALSE)

cat("Done. Outputs in results_spss_R/\n")
