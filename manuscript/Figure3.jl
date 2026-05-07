#####
# This script generates Figure 1 for PCI-VEP-Mbiome Manuscript
#####
using CSV
using DataFrames
using Chain
using MultivariateStats
using HypothesisTests
using MultipleTesting
using Distributions
using CategoricalArrays
using Microbiome
using Distances
using Dictionaries
using Clustering
using Diversity
using DecisionTree
using Random
using JLD2
using Statistics
using CairoMakie
using ThreadsX
using GLM
using CodecZlib
using ColorSchemes
using Leap
using PCIVEPMBiome2024
using Arrow
using SparseArrays
using StatsBase
using ColorSchemes

#####
# 1. Loading data
#####
paper_sample_set = readlines(joinpath(Base.pwd(), "data", "kept_complete_samples.txt"))

mdata_cols, khula_pci_mbiome_data,
  unstratified_unirefs, stratified_unirefs,
  unstratified_unirefs_filtered, stratified_unirefs_filtered = load_pcimbiome_data(paper_sample_set)

unimdata = DataFrame(get(unstratified_unirefs_filtered))

## Optionally export inputs for `R::adonis2`

CSV.write(
    joinpath(Base.pwd(), "data", "adonis", "taxonmic_inputs.csv"),
    khula_pci_mbiome_data
)

CSV.write(
    joinpath(Base.pwd(), "data", "adonis", "functional_inputs.csv"),
    hcat(
        DataFrame(get(unstratified_unirefs_filtered)),
        DataFrame(collect(abundances(unstratified_unirefs)'), featurenames(unstratified_unirefs))
    )
)

#####
# Distance analysis
#####

unidm = Microbiome.braycurtis(unstratified_unirefs)
uni_MDS_results = fit(MDS, unidm; maxoutdim = 20, distances=true)
uni_MDS_columns = DataFrame(:MDS1 => uni_MDS_results.U[:,1], :MDS2 => uni_MDS_results.U[:,2], :MDS3 => uni_MDS_results.U[:,3], :MDS4 => uni_MDS_results.U[:,4], :MDS5 => uni_MDS_results.U[:,5])
uni_MDS_variances = uni_MDS_results.λ ./ sum(uni_MDS_results.λ)

spedm = Distances.pairwise(BrayCurtis(), Matrix(select(khula_pci_mbiome_data, Not(mdata_cols))), dims = 1)
spe_MDS_results = fit(MDS, spedm; maxoutdim = 20, distances=true)
spe_MDS_columns = DataFrame(:MDS1 => spe_MDS_results.U[:,1], :MDS2 => spe_MDS_results.U[:,2], :MDS3 => spe_MDS_results.U[:,3], :MDS4 => spe_MDS_results.U[:,4], :MDS5 => spe_MDS_results.U[:,5])
spe_MDS_variances = spe_MDS_results.λ ./ sum(spe_MDS_results.λ)

## Printing some useful stats for text
spects = sum(Matrix(select(khula_pci_mbiome_data, Not(mdata_cols))) .> 0.0, dims = 2)
gfcts = sum(unstratified_unirefs.abundances .> 0.0, dims = 1)
println("---------------- PROFILE STATS -----------------")
println("Species count quantiles: $(quantile(spects))")
println("Species count MEAN/SD: $(mean(spects)), $(std(spects))")
println("Gene Functions count quantiles: $(quantile(gfcts))")
println("Gene Functions count MEAN/SD: $(mean(gfcts)), $(std(gfcts))")
println("--------------- MAIN MATTER STATS ---------------")
println("GENE/FUNCTIONAL PCoA Variance explained by 3 PCs: $(sum(uni_MDS_variances[1:3]))")
println("GENE/FUNCTIONAL PCoA Variance explained by PC1: $(uni_MDS_variances[1])")
println("GENE/FUNCTIONAL PCoA Variance explained by PC2: $(uni_MDS_variances[2])")
println("GENE/FUNCTIONAL PCoA Variance explained by PC3: $(uni_MDS_variances[3])")
println("-------------------------------------------------")
println("SPECIES/TAXONOMIC PCoA Variance explained by 3 PCs: $(sum(spe_MDS_variances[1:3]))")
println("SPECIES/TAXONOMIC PCoA Variance explained by PC1: $(spe_MDS_variances[1])")
println("SPECIES/TAXONOMIC PCoA Variance explained by PC2: $(spe_MDS_variances[2])")
println("SPECIES/TAXONOMIC PCoA Variance explained by PC3: $(spe_MDS_variances[3])")
println("-------------------------------------------------")
println("SPECIES/TAXONOMIC PCoA rsquare between PC1 and B_breve: $(cor(spe_MDS_results.U[:,1], khula_pci_mbiome_data.Bifidobacterium_breve)^2)")
println("SPECIES/TAXONOMIC PCoA rsquare between PC1 and B_longum: $(cor(spe_MDS_results.U[:,1], khula_pci_mbiome_data.Bifidobacterium_longum)^2)")
println("SPECIES/TAXONOMIC PCoA rsquare between PC1 and VOB: $(cor(spe_MDS_results.U[:,1], khula_pci_mbiome_data.InfantVisAtt)^2)")
println("-------------------------------------------------")
println("SPECIES/TAXONOMIC PCoA rsquare between PC2 and B_breve: $(cor(spe_MDS_results.U[:,2], khula_pci_mbiome_data.Bifidobacterium_breve)^2)")
println("SPECIES/TAXONOMIC PCoA rsquare between PC2 and B_longum: $(cor(spe_MDS_results.U[:,2], khula_pci_mbiome_data.Bifidobacterium_longum)^2)")
println("SPECIES/TAXONOMIC PCoA rsquare between PC2 and B_bifidum: $(cor(spe_MDS_results.U[:,2], khula_pci_mbiome_data.Bifidobacterium_bifidum)^2)")
println("SPECIES/TAXONOMIC PCoA rsquare between PC2 and R_gnavus: $(cor(spe_MDS_results.U[:,2], khula_pci_mbiome_data.Ruminococcus_gnavus)^2)")
println("SPECIES/TAXONOMIC PCoA rsquare between PC2 and E_ramosum: $(cor(spe_MDS_results.U[:,2], khula_pci_mbiome_data.Erysipelatoclostridium_ramosum)^2)")
println("SPECIES/TAXONOMIC PCoA rsquare between PC2 and E_coli: $(cor(spe_MDS_results.U[:,2], khula_pci_mbiome_data.Escherichia_coli)^2)")
println("SPECIES/TAXONOMIC PCoA rsquare between PC2 and VOB: $(cor(spe_MDS_results.U[:,2], khula_pci_mbiome_data.InfantVisAtt)^2)")

# PERMANOVAs
mdlabels = [ "Age at collection", "Child sex", "Maternal education", "Feeding status", "Delivery mode", "Maternal\nunpredictability", "Visual orienting\nbehavior (VOB)"]

spe_pmn = permanovas(
    [ spedm ], 
    [
        khula_pci_mbiome_data.mbiome_sample_age,
        khula_pci_mbiome_data.child_sex,
        khula_pci_mbiome_data.maternal_education,
        khula_pci_mbiome_data.feeding_state,
        khula_pci_mbiome_data.delivery_mode,
        khula_pci_mbiome_data.MaternalEntropy,
        khula_pci_mbiome_data.InfantVisAtt
    ]; commlabels = ["Taxa"], mdlabels
)

uni_pmn = permanovas(
    [ unidm ], 
    [
        unimdata.mbiome_sample_age,
        unimdata.child_sex,
        unimdata.maternal_education,
        unimdata.feeding_state,
        unimdata.delivery_mode,
        unimdata.MaternalEntropy,
        unimdata.InfantVisAtt
    ]; commlabels = ["UniRef90s"], mdlabels
)

## Create Figure 3

fig = Figure(size = (1000,1000))
Top_Row = GridLayout(fig[1,:]; alignmode = Inside())
Bottom_Row = GridLayout(fig[2,:]; alignmode = Inside())
A_Subfig = GridLayout(Top_Row[:,1]; alignmode = Mixed(left = 10, right = 60))
BCD_Subfig = GridLayout(Top_Row[:,2]; alignmode = Inside())
EF_Subfig = GridLayout(Bottom_Row[1,1]; alignmode = Outside())

## Panel A: PERMANOVAS
axA = Axis(
    A_Subfig[1, 1];
    xticklabelsize = 16,
    yticklabelsize = 16,
    title = "PERMANOVAs",
)

pmnA = plot_permanovas!(axA, vcat(spe_pmn, uni_pmn))

## Panels BCD: ADONIS and Shannon and Visual vs Feed

## B: PERMANOVA rsquares
adonis_outputs = CSV.read(
    joinpath(Base.pwd(), "manuscript", "tables", "adonis_results.csv"),
    DataFrame;
    stringtype = String
)

plotdf = subset(adonis_outputs, :profile => ByRow(==("taxa")))

# Desired order
analysis_order = [
    "visual_only",
    "visual_after_feeding",
    "feeding_only",
    "feeding_after_visual",
    "joint_model"
]

# Pretty labels
label_map = Dict(
    "visual_only" => "VOB only",
    "feeding_only" => "Feeding only",
    "joint_model" => "Joint model",
    "visual_after_feeding" => "VOB | Feeding\n(VOB marginal)",
    "feeding_after_visual" => "Feeding | VOB\n(Feeding marginal)"
)

plotdf.analysis = categorical(
    plotdf.analysis;
    ordered=true,
    levels=analysis_order
)

sort!(plotdf, :analysis)

adonis_x = 1:nrow(plotdf)
adonis_y = plotdf.R2

axB = Axis(
    BCD_Subfig[1, 1:2],
    ylabel = "Explained variance (R²)",
    xlabel = "",
    xticks = (adonis_x, [label_map[a] for a in plotdf.analysis]),
    title = "Combined PERMANOVAs for feeding and VOB",
    xticklabelrotation = pi/6
)
hidedecorations!(
    axB;
    label = false, ticklabels = false, ticks = false,
    grid = true, minorgrid = false, minorticks = false
)

bpB = barplot!(axB, adonis_x, adonis_y)

# significance labels
tx = text!(axB,
    collect(zip(collect(adonis_x), adonis_y .+ 0.003));
    text = string.(plotdf.sig_code),
    align = (:center, :bottom),
    fontsize = 18
)

ylims!(axB, 0.0, maximum(adonis_y) + 0.03)

## C,D: Boxplots
xs_mapper = Dict(
    "ExclBreastFed" => 1,
    "Mixed" => 2,
    "ExclFormulaFed" => 3
)

feeding_cmap = Dict(
    "ExclBreastFed" => (:blue, 0.7),
    "Mixed" => (:purple, 0.7),
    "ExclFormulaFed" => (:red, 0.7),
)

axC = Axis(
    BCD_Subfig[2, 1];
    xticklabelsize = 16,
    yticklabelsize = 16,
    xticks = (1:3, ["Excl.\nBF", "Mixed", "Excl.\nFF"]),
    title = "Shannon index by\nfeeding practice"
)
hidedecorations!(
    axC;
    label = false, ticklabels = false, ticks = false,
    grid = true, minorgrid = false, minorticks = false
)
ylims!(axC, [0.0, 3.3])

axD = Axis(
    BCD_Subfig[2, 2];
    xticklabelsize = 16,
    yticklabelsize = 16,
    xticks = (1:3, ["Excl.\nBF", "Mixed", "Excl.\nFF"]),
    title = "Infant VOB by\nfeeding practice"
)
hidedecorations!(
    axD;
    label = false, ticklabels = false, ticks = false,
    grid = true, minorgrid = false, minorticks = false
)
ylims!(axD, [0.0, 1.5e-4])

bpB1 = boxplot!(
    axC,
    [ xs_mapper[el] for el in khula_pci_mbiome_data.feeding_state ],
    khula_pci_mbiome_data.alpha_shannon;
    color = [ feeding_cmap[el] for el in khula_pci_mbiome_data.feeding_state ]
)

shannon_bf_mx_pval = HypothesisTests.pvalue(
    UnequalVarianceTTest(
        khula_pci_mbiome_data.alpha_shannon[khula_pci_mbiome_data.feeding_state .== "ExclBreastFed"],
        khula_pci_mbiome_data.alpha_shannon[khula_pci_mbiome_data.feeding_state .== "Mixed"]
    )
)

shannon_mx_ff_pval = HypothesisTests.pvalue(
    UnequalVarianceTTest(
        khula_pci_mbiome_data.alpha_shannon[khula_pci_mbiome_data.feeding_state .== "Mixed"],
        khula_pci_mbiome_data.alpha_shannon[khula_pci_mbiome_data.feeding_state .== "ExclFormulaFed"]
    )
)

shannon_bf_ff_pval = HypothesisTests.pvalue(
    UnequalVarianceTTest(
        khula_pci_mbiome_data.alpha_shannon[khula_pci_mbiome_data.feeding_state .== "ExclBreastFed"],
        khula_pci_mbiome_data.alpha_shannon[khula_pci_mbiome_data.feeding_state .== "ExclFormulaFed"]
    )
)

bpB2 = boxplot!(
    axD,
    [ xs_mapper[el] for el in khula_pci_mbiome_data.feeding_state ],
    khula_pci_mbiome_data.InfantVisAtt;
    color = [ feeding_cmap[el] for el in khula_pci_mbiome_data.feeding_state ]
)

visual_bf_mx_pval = HypothesisTests.pvalue(
    UnequalVarianceTTest(
        khula_pci_mbiome_data.InfantVisAtt[khula_pci_mbiome_data.feeding_state .== "ExclBreastFed"],
        khula_pci_mbiome_data.InfantVisAtt[khula_pci_mbiome_data.feeding_state .== "Mixed"]
    )
)

visual_mx_ff_pval = HypothesisTests.pvalue(
    UnequalVarianceTTest(
        khula_pci_mbiome_data.InfantVisAtt[khula_pci_mbiome_data.feeding_state .== "Mixed"],
        khula_pci_mbiome_data.InfantVisAtt[khula_pci_mbiome_data.feeding_state .== "ExclFormulaFed"]
    )
)

visual_bf_ff_pval = HypothesisTests.pvalue(
    UnequalVarianceTTest(
        khula_pci_mbiome_data.InfantVisAtt[khula_pci_mbiome_data.feeding_state .== "ExclBreastFed"],
        khula_pci_mbiome_data.InfantVisAtt[khula_pci_mbiome_data.feeding_state .== "ExclFormulaFed"]
    )
)

## Panel C: Gene PCoA

axE = Axis(
    EF_Subfig[1,1],
    xlabel = "MDS1 ("*string(round(100*uni_MDS_variances[1]; digits = 2))*"%)",
    ylabel = "MDS3 ("*string(round(100*uni_MDS_variances[3]; digits = 2))*"%)",
    title = "Functional Profiles (Genes/UniRef90s)",
    aspect = 1.0
)
hidedecorations!(
    axE;
    label = false,
    ticklabels = false,
    ticks = false,
    grid = true,
    minorgrid = false,
    minorticks = false
)
scE = scatter!(
    axE,
    uni_MDS_columns[:,1],
    uni_MDS_columns[:,3];
    markersize = 12,
    color = unimdata.InfantVisAtt,
    colormap = :viridis
)

## Panel C: Species PCoA

axF = Axis(
    EF_Subfig[1,2],
    xlabel = "MDS1 ("*string(round(100*spe_MDS_variances[1]; digits = 2))*"%)",
    ylabel = "MDS2 ("*string(round(100*spe_MDS_variances[2]; digits = 2))*"%)",
    title = "Taxonomic Profiles (species)",
    aspect = 1.0
)
hidedecorations!(
    axF;
    label = false,
    ticklabels = false,
    ticks = false,
    grid = true,
    minorgrid = false,
    minorticks = false
)
scF = scatter!(
    axF,
    spe_MDS_columns[:,1],
    spe_MDS_columns[:,2];
    markersize = 12,
    color = khula_pci_mbiome_data.InfantVisAtt,
    colormap = :viridis
)

cbD = Colorbar(
    EF_Subfig[1,3], scF;
    label = "Visual Orienting Behavior",
    vertical = true,
    tellheight = true, tellwidth = true,
    alignmode = Outside())

## Layout

Label(A_Subfig[1, 1, TopLeft()], "a", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(BCD_Subfig[1, :, TopLeft()], "b", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(BCD_Subfig[2, 1, TopLeft()], "c", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(BCD_Subfig[2, 2, TopLeft()], "d", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(EF_Subfig[1, 1, TopLeft()], "e", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(EF_Subfig[1, 2, TopLeft()], "f", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())

colgap!(fig.layout, 0)
rowgap!(fig.layout, 0)

rowsize!(fig.layout, 1, Relative(0.50))
rowsize!(fig.layout, 2, Relative(0.50))

rowsize!(BCD_Subfig, 1, Relative(0.40))
rowsize!(BCD_Subfig, 2, Relative(0.60))

fig

save(joinpath(Base.pwd(),"manuscript", "figures", "Figure3.png"), fig)
save(joinpath(Base.pwd(),"manuscript", "figures", "Figure3.eps"), fig)
save(joinpath(Base.pwd(),"manuscript", "figures", "Figure3.pdf"), fig)

#####
# Exporting Supplementary Figures
#####

include(joinpath(Base.pwd(), "manuscript", "FigureS4.jl"))
include(joinpath(Base.pwd(), "manuscript", "FigureS5.jl"))
include(joinpath(Base.pwd(), "manuscript", "FigureS6.jl"))