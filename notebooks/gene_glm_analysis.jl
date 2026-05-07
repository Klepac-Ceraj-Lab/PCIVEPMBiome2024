using DataFrames
using GLM
using StatsModels
using Tables
using CSV
using CategoricalArrays
using StatsBase
using ProgressMeter

# Input data for all GLM models

indf = deepcopy(unimdata)

indf.mbiome_sample_age = Leap.rangenormalize(indf.mbiome_sample_age)
indf.maternal_education = Leap.rangenormalize(indf.maternal_education)

indf.feeding_state = Float64.(indf.feeding_state .== "ExclBreastFed") 
indf.delivery_mode = Float64.(indf.delivery_mode .== "Caesarean")

indf.MaternalEntropy = Leap.rangenormalize(indf.MaternalEntropy)
indf.InfantVisAtt = Leap.rangenormalize(indf.InfantVisAtt)
## convenience
featlist = collect(featurenames(unstratified_unirefs_filtered))

# MODEL 1
outfiles1 = (
    joinpath(Base.pwd(),"manuscript", "gene_glm", "lms_Model1_mbiome_age.csv"),
    joinpath(Base.pwd(),"manuscript", "gene_glm", "lms_Model1_child_sex.csv"),
    joinpath(Base.pwd(),"manuscript", "gene_glm", "lms_Model1_maternal_education.csv"),
    joinpath(Base.pwd(),"manuscript", "gene_glm", "lms_Model1_feeding_state.csv"),
    joinpath(Base.pwd(),"manuscript", "gene_glm", "lms_Model1_delivery_mode.csv")
)

## Make a progress bar so I'm not driven mad by the expectation
p = Progress(length(featlist))

## Actually fit the models using `ThreadsX.map`
res1 = ThreadsX.map(featlist) do feat
    next!(p)
    ab = vec(abundances(unstratified_unirefs_filtered[feat, :]))
    fit_one(indf, ab,
        @formula(bug ~ mbiome_sample_age + child_sex + maternal_education + feeding_state + delivery_mode),
        ProbitLink(),
        ["mbiome_sample_age", "child_sex", "maternal_education", "feeding_state", "delivery_mode"];
        feature = string(feat))
end

## res1 is Vector{Vector{NamedTuple}} (length 2 per feature). Split & write:
for i in 1:length(outfiles1)
    write_results(outfiles1[i], getindex.(res1, i))
end

# MODEL 2
outfiles2 = (
    joinpath(Base.pwd(),"manuscript", "gene_glm", "lms_Model2_mbiome_age.csv"),
    joinpath(Base.pwd(),"manuscript", "gene_glm", "lms_Model2_child_sex.csv"),
    joinpath(Base.pwd(),"manuscript", "gene_glm", "lms_Model2_maternal_education.csv"),
    joinpath(Base.pwd(),"manuscript", "gene_glm", "lms_Model2_feeding_state.csv"),
    joinpath(Base.pwd(),"manuscript", "gene_glm", "lms_Model2_delivery_mode.csv"),
    joinpath(Base.pwd(),"manuscript", "gene_glm", "lms_Model2_visual.csv"),
    joinpath(Base.pwd(),"manuscript", "gene_glm", "lms_Model2_entropy.csv")
)

p = Progress(length(featlist))
res2 = ThreadsX.map(featlist) do feat
    next!(p)
    ab = vec(abundances(unstratified_unirefs_filtered[feat, :]))
    fit_one(indf, ab,
        @formula(bug ~ mbiome_sample_age + child_sex + maternal_education + feeding_state + delivery_mode + InfantVisAtt + MaternalEntropy),
        ProbitLink(),
        ["mbiome_sample_age", "child_sex", "maternal_education", "feeding_state", "delivery_mode", "InfantVisAtt", "MaternalEntropy"];
        feature = string(feat))
end

for i in 1:length(outfiles2)
    write_results(outfiles2[i], getindex.(res2, i))
end

# MODEL 3
outfiles3 = (
    joinpath(Base.pwd(),"manuscript", "gene_glm", "lms_Model3_mbiome_age.csv"),
    joinpath(Base.pwd(),"manuscript", "gene_glm", "lms_Model3_child_sex.csv"),
    joinpath(Base.pwd(),"manuscript", "gene_glm", "lms_Model3_maternal_education.csv"),
    joinpath(Base.pwd(),"manuscript", "gene_glm", "lms_Model3_feeding_state.csv"),
    joinpath(Base.pwd(),"manuscript", "gene_glm", "lms_Model3_delivery_mode.csv"),
    joinpath(Base.pwd(),"manuscript", "gene_glm", "lms_Model3_visual.csv"),
    joinpath(Base.pwd(),"manuscript", "gene_glm", "lms_Model3_entropy.csv"),
    joinpath(Base.pwd(),"manuscript", "gene_glm", "lms_Model3_product.csv")
)

p = Progress(length(featlist))
res3 = ThreadsX.map(featlist) do feat
    next!(p)
    ab = vec(abundances(unstratified_unirefs_filtered[feat, :]))
    fit_one(indf, ab,
        @formula(bug ~ mbiome_sample_age + child_sex + maternal_education + feeding_state + delivery_mode + InfantVisAtt + MaternalEntropy + MaternalEntropy & InfantVisAtt),
        ProbitLink(),
        ["mbiome_sample_age", "child_sex", "maternal_education", "feeding_state", "delivery_mode", "InfantVisAtt", "MaternalEntropy", "MaternalEntropy & InfantVisAtt"];
        feature = string(feat))
end

for i in 1:length(outfiles3)
    write_results(outfiles3[i], getindex.(res3, i))
end