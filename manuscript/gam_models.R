## install.packages("mgcv")
library(mgcv)

## 1. Load data:

input_df <- read.csv("manuscript/tables/mdata_table_in_bb3.csv")

input_df$MaternalEntropy <- scale(input_df$MaternalEntropy)
input_df$InfantVisAtt <- scale(input_df$InfantVisAtt)

gam_Bbreve_joint <- gam(Bifidobacterium_breve ~ s(MaternalEntropy, InfantVisAtt), data = input_df)
gam_Blongum_joint <- gam(Bifidobacterium_longum ~ s(MaternalEntropy, InfantVisAtt), data = input_df)
gam_Bbreve_sep <- gam(Bifidobacterium_breve ~ s(MaternalEntropy) + s(InfantVisAtt), data = input_df)
gam_Blongum_sep <- gam(Bifidobacterium_longum ~ s(MaternalEntropy) + s(InfantVisAtt), data = input_df)
gam_Bbreve_int <- gam(Bifidobacterium_breve ~ s(MaternalEntropy) + s(InfantVisAtt) + s(MaternalEntropy, InfantVisAtt), data = input_df)
gam_Blongum_int <- gam(Bifidobacterium_longum ~ s(MaternalEntropy) + s(InfantVisAtt) + s(MaternalEntropy, InfantVisAtt), data = input_df)

summary(gam_Bbreve_joint)
summary(gam_Blongum_joint)
summary(gam_Bbreve_sep)
summary(gam_Blongum_sep)
summary(gam_Bbreve_int)
summary(gam_Blongum_int)

input_df$Bdiff <- input_df$Bifidobacterium_breve - input_df$Bifidobacterium_longum
input_df$Blogratio <- log2( (input_df$Bifidobacterium_breve + 2*.Machine$double.eps ) / (input_df$Bifidobacterium_longum + 2*.Machine$double.eps ) )
input_df$Bbrevedom <- input_df$Bifidobacterium_breve / (input_df$Bifidobacterium_longum + input_df$Bifidobacterium_breve + 2*.Machine$double.eps)
input_df$Blongumdom <- input_df$Bifidobacterium_longum / (input_df$Bifidobacterium_longum + input_df$Bifidobacterium_breve + 2*.Machine$double.eps)

gam_Bdiff_joint <- gam(Bdiff ~ s(MaternalEntropy, InfantVisAtt), data = input_df)
gam_Bdiff_sep <- gam(Bdiff ~ s(MaternalEntropy) + s(InfantVisAtt), data = input_df)
gam_Bdiff_int <- gam(Bdiff ~ s(MaternalEntropy) + s(InfantVisAtt) + s(MaternalEntropy, InfantVisAtt), data = input_df)

summary(gam_Bdiff_joint)
summary(gam_Bdiff_sep)
summary(gam_Bdiff_int)

gam_Blogratio_joint <- gam(Blogratio ~ s(MaternalEntropy, InfantVisAtt), data = input_df)
gam_Blogratio_sep <- gam(Blogratio ~ s(MaternalEntropy) + s(InfantVisAtt), data = input_df)
gam_Blogratio_int <- gam(Blogratio ~ s(MaternalEntropy) + s(InfantVisAtt) + s(MaternalEntropy, InfantVisAtt), data = input_df)

summary(gam_Blogratio_joint)
summary(gam_Blogratio_sep)
summary(gam_Blogratio_int)

gam_Bbrevedom_joint <- gam(Bbrevedom ~ s(MaternalEntropy, InfantVisAtt), data = input_df)
gam_Bbrevedom_sep <- gam(Bbrevedom ~ s(MaternalEntropy) + s(InfantVisAtt), data = input_df)
gam_Bbrevedom_int <- gam(Bbrevedom ~ s(MaternalEntropy) + s(InfantVisAtt) + s(MaternalEntropy, InfantVisAtt), data = input_df)

summary(gam_Bbrevedom_joint)
summary(gam_Bbrevedom_sep)
summary(gam_Bbrevedom_int)

gam_Blongumdom_joint <- gam(Blongumdom ~ s(MaternalEntropy, InfantVisAtt), data = input_df)
gam_Blongumdom_sep <- gam(Blongumdom ~ s(MaternalEntropy) + s(InfantVisAtt), data = input_df)
gam_Blongumdom_int <- gam(Blongumdom ~ s(MaternalEntropy) + s(InfantVisAtt) + s(MaternalEntropy, InfantVisAtt), data = input_df)

summary(gam_Blongumdom_joint)
summary(gam_Blongumdom_sep)
summary(gam_Blongumdom_int)

#####
# Re-scaling the space
#####

input_df$radius <- scale(sqrt(input_df$MaternalEntropy^2 + input_df$InfantVisAtt^2))
input_df$angle <- atan2(input_df$MaternalEntropy, input_df$InfantVisAtt)

gam_Bbreve_polar_joint <- gam(Bifidobacterium_breve ~ s(radius, angle), data = input_df)
gam_Blongum_polar_joint <- gam(Bifidobacterium_longum ~ s(radius, angle), data = input_df)
gam_Bbreve_polar_sep <- gam(Bifidobacterium_breve ~ s(radius) + s(angle), data = input_df)
gam_Blongum_polar_sep <- gam(Bifidobacterium_longum ~ s(radius) + s(angle), data = input_df)
gam_Bbreve_polar_int <- gam(Bifidobacterium_breve ~ s(radius) + s(angle) + s(radius, angle), data = input_df)
gam_Blongum_polar_int <- gam(Bifidobacterium_longum ~ s(radius) + s(angle) + s(radius, angle), data = input_df)

summary(gam_Bbreve_polar_joint)
summary(gam_Blongum_polar_joint)
summary(gam_Bbreve_polar_sep)
summary(gam_Blongum_polar_sep)
summary(gam_Bbreve_polar_joint)
summary(gam_Blongum_polar_joint)

#####
# Inverting the problem!
#####

gam_MatEnt_joint <- gam(MaternalEntropy ~ s(Bifidobacterium_breve, Bifidobacterium_longum), data = input_df)
gam_InfVisAtt_joint <- gam(InfantVisAtt ~ s(Bifidobacterium_breve, Bifidobacterium_longum), data = input_df)
gam_MatEnt_sep <- gam(MaternalEntropy ~ s(Bifidobacterium_breve) + s(Bifidobacterium_longum), data = input_df)
gam_InfVisAtt_sep <- gam(InfantVisAtt ~ s(Bifidobacterium_breve) + s(Bifidobacterium_longum), data = input_df)
gam_MatEnt_int <- gam(MaternalEntropy ~ s(Bifidobacterium_breve) + s(Bifidobacterium_longum)  s(Bifidobacterium_breve, Bifidobacterium_longum), data = input_df)
gam_InfVisAtt_int <- gam(InfantVisAtt ~ s(Bifidobacterium_breve) + s(Bifidobacterium_longum) + s(Bifidobacterium_breve, Bifidobacterium_longum), data = input_df)

summary(gam_MatEnt_joint)
summary(gam_InfVisAtt_joint)
summary(gam_MatEnt_sep)
summary(gam_InfVisAtt_sep)
summary(gam_MatEnt_int)
summary(gam_InfVisAtt_int)