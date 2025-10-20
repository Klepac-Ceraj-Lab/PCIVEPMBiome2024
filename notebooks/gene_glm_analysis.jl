using DataFrames
using GLM
using StatsModels
using Tables
using CSV
using CategoricalArrays
using StatsBase
using ProgressMeter

# Input data for all GLM models

indf = DataFrame(get(unstratified_unirefs_filtered))
indf.InfantVisAtt = Leap.rangenormalize(indf.InfantVisAtt)
indf.MaternalEntropy = Leap.rangenormalize(indf.MaternalEntropy)
indf.mbiome_sample_age = Leap.rangenormalize(indf.mbiome_sample_age)

## convenience
featlist = collect(featurenames(unstratified_unirefs_filtered))

# MODEL 1

outfiles1 = (
    "manuscript/FSEA/lms_Model1_mbiome_age.csv",
    "manuscript/FSEA/lms_Model1_child_sex.csv"
)

## Make a progress bar so I'm not driven mad by the expectation
p = Progress(length(featlist))

## Actually fit the models using `ThreadsX.map`
res1 = ThreadsX.map(featlist) do feat
    next!(p)
    ab = vec(abundances(unstratified_unirefs_filtered[feat, :]))
    fit_one(indf, ab,
        @formula(bug ~ mbiome_sample_age + child_sex),
        ProbitLink(),
        ["mbiome_sample_age", "child_sex"];
        feature = string(feat))
end

## res1 is Vector{Vector{NamedTuple}} (length 2 per feature). Split & write:
write_results(outfiles1[1], getindex.(res1, 1))
write_results(outfiles1[2], getindex.(res1, 2))

# MODEL 2

outfiles2 = (
    "manuscript/FSEA/lms_Model2_mbiome_age.csv",
    "manuscript/FSEA/lms_Model2_child_sex.csv",
    "manuscript/FSEA/lms_Model2_visual.csv",
    "manuscript/FSEA/lms_Model2_entropy.csv"
)

p = Progress(length(featlist))
res2 = ThreadsX.map(featlist) do feat
    next!(p)
    ab = vec(abundances(unstratified_unirefs_filtered[feat, :]))
    fit_one(indf, ab,
        @formula(bug ~ mbiome_sample_age + child_sex + InfantVisAtt + MaternalEntropy),
        ProbitLink(),
        ["mbiome_sample_age", "child_sex", "InfantVisAtt", "MaternalEntropy"];
        feature = string(feat))
end

write_results(outfiles2[1], getindex.(res2, 1))
write_results(outfiles2[2], getindex.(res2, 2))
write_results(outfiles2[3], getindex.(res2, 3))
write_results(outfiles2[4], getindex.(res2, 4))

# MODEL 3

outfiles3 = (
    "manuscript/FSEA/lms_Model3_mbiome_age.csv",
    "manuscript/FSEA/lms_Model3_child_sex.csv",
    "manuscript/FSEA/lms_Model3_visual.csv",
    "manuscript/FSEA/lms_Model3_entropy.csv",
    "manuscript/FSEA/lms_Model3_product.csv"
)

p = Progress(length(featlist))
res3 = ThreadsX.map(featlist) do feat
    next!(p)
    ab = vec(abundances(unstratified_unirefs_filtered[feat, :]))
    fit_one(indf, ab,
        @formula(bug ~ mbiome_sample_age + child_sex + InfantVisAtt + MaternalEntropy + MaternalEntropy & InfantVisAtt),
        ProbitLink(),
        ["mbiome_sample_age", "child_sex", "InfantVisAtt", "MaternalEntropy", "MaternalEntropy & InfantVisAtt"];
        feature = string(feat))
end

write_results(outfiles3[1], getindex.(res3, 1))
write_results(outfiles3[2], getindex.(res3, 2))
write_results(outfiles3[3], getindex.(res3, 3))
write_results(outfiles3[4], getindex.(res3, 4))
write_results(outfiles3[5], getindex.(res3, 5))

# BONUS: MODEL 2 SUBTYPE 1 (ENTROPY ONLY)

outfiles21 = (
    "manuscript/FSEA/lms_Model21_mbiome_age.csv",
    "manuscript/FSEA/lms_Model21_child_sex.csv",
    "manuscript/FSEA/lms_Model21_entropy.csv"
)

p = Progress(length(featlist))
res21 = ThreadsX.map(featlist) do feat
    next!(p)
    ab = vec(abundances(unstratified_unirefs_filtered[feat, :]))
    fit_one(indf, ab,
        @formula(bug ~ mbiome_sample_age + child_sex + MaternalEntropy),
        ProbitLink(),
        ["mbiome_sample_age", "child_sex", "MaternalEntropy"];
        feature = string(feat))
end

write_results(outfiles21[1], getindex.(res21, 1))
write_results(outfiles21[2], getindex.(res21, 2))
write_results(outfiles21[3], getindex.(res21, 3))

# BONUS: MODEL 2 SUBTYPE 2 (VISUAL ONLY)

outfiles22 = (
    "manuscript/FSEA/lms_Model22_mbiome_age.csv",
    "manuscript/FSEA/lms_Model22_child_sex.csv",
    "manuscript/FSEA/lms_Model22_visual.csv"
)

p = Progress(length(featlist))
res22 = ThreadsX.map(featlist) do feat
    next!(p)
    ab = vec(abundances(unstratified_unirefs_filtered[feat, :]))
    fit_one(indf, ab,
        @formula(bug ~ mbiome_sample_age + child_sex + InfantVisAtt),
        ProbitLink(),
        ["mbiome_sample_age", "child_sex", "InfantVisAtt"];
        feature = string(feat))
end

write_results(outfiles22[1], getindex.(res22, 1))
write_results(outfiles22[2], getindex.(res22, 2))
write_results(outfiles22[3], getindex.(res22, 3))