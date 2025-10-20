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

#####
# 1. Loading data
#####
include(joinpath(Base.pwd(), "notebooks", "load_data.jl"))

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

## Create Figure 1

fig = Figure(size = (1000,1100))
Left_Col = GridLayout(fig[:,1]; alignmode = Inside())
Right_Col = GridLayout(fig[:,2]; alignmode = Inside())
A_Subfig = GridLayout(Left_Col[1,:]; alignmode = Mixed(left = 10, right = 60))
B_Subfig = GridLayout(Right_Col[1,:]; alignmode = Inside())
C_Subfig = GridLayout(Left_Col[2,:]; alignmode = Inside())
D_Subfig = GridLayout(Right_Col[2,:]; alignmode = Inside())

## Panel A: PERMANOVAS

mdlabels = [ "Age at collection", "Child sex", "Maternal education", "Maternal\nunpredictability", "Visual orienting\nbehavior (VOB)"]

spe_pmn = permanovas(
    [ spedm ], 
    [
        khula_pci_mbiome_data.mbiome_sample_age,
        khula_pci_mbiome_data.child_sex,
        khula_pci_mbiome_data.MaternalEducation,
        khula_pci_mbiome_data.MaternalEntropy,
        khula_pci_mbiome_data.InfantVisAtt
    ]; commlabels = ["Taxa"], mdlabels
)

uni_pmn = permanovas(
    [ unidm ], 
    [
        unimdata.mbiome_sample_age,
        unimdata.child_sex,
        unimdata.MaternalEducation,
        unimdata.MaternalEntropy,
        unimdata.InfantVisAtt,
    ]; commlabels = ["UniRef90s"], mdlabels
)

axA = Axis(
    A_Subfig[1, 1];
    xticklabelsize = 16,
    yticklabelsize = 16,
    title = "PERMANOVAs",
)

plot_permanovas!(axA, vcat(spe_pmn, uni_pmn))

## Panel B: Gene PCoA

axB = Axis(
    B_Subfig[1,1],
    xlabel = "MDS1 ("*string(round(100*uni_MDS_variances[1]; digits = 2))*"%)",
    ylabel = "MDS3 ("*string(round(100*uni_MDS_variances[3]; digits = 2))*"%)",
    title = "Functional Profiles (Genes/UniRef90s)",
    aspect = 1.0
)
hidedecorations!(
    axB;
    label = false,
    ticklabels = false,
    ticks = false,
    grid = true,
    minorgrid = false,
    minorticks = false
)
scB = scatter!(
    axB,
    uni_MDS_columns[:,1],
    uni_MDS_columns[:,3];
    markersize = 12,
    color = unimdata.InfantVisAtt,
    colormap = :viridis
)

## Panel C: PCI Topology and Bifidobacteria

bifido_cdict = Dict(
    "dom_breve" => (:purple4, 1.0),
    "more_breve" => (:darkorchid2, 0.8),
    "null" => (:gray, 0.6),
    "more_longum" => (:orange, 0.8),
    "dom_longum" => (:darkorange2, 1.0)
)

function assign_bifido_class(tt::T where T<: Tuple)

    if (tt[1] > 0.75)
        return "dom_breve"
    elseif (tt[2] > 0.75)
        return "dom_longum"
    elseif ( (tt[1] > 0.05) & ( ( (tt[1] + 0.01) / (tt[2] + 0.01) ) > 5 ) )
        return "more_breve"
    elseif ( (tt[2] > 0.05) & ( ( (tt[2] + 0.01) / (tt[1] + 0.01) ) > 5 ) )
        return "more_longum"
    else
        return "null"
    end

end

function assign_bifido_class(df::DataFrame)

    zipobj = zip( df[:, "Bifidobacterium_breve"], df[:, "Bifidobacterium_longum"] )

    map(assign_bifido_class, zipobj)

end

axC = Axis(
    C_Subfig[1,1],
    xlabel = "Maternal unpredictability",
    ylabel = "Visual Orienting Behavior (VOB) x10^-5",
    yticks = ([0.0, 0.00005, 0.00010], ["0.0", "5.0", "10.0"]),
    title = "Bifidobacteria across VOB x unpredictability levels",
    aspect = 1.0
)
hidedecorations!(
    axC;
    label = false,
    ticklabels = false,
    ticks = false,
    grid = true,
    minorgrid = false,
    minorticks = false
)
scD = scatter!(
    axC,
    khula_pci_mbiome_data.MaternalEntropy,
    khula_pci_mbiome_data.InfantVisAtt;
    markersize = 12,
    color = [ bifido_cdict[el] for el in assign_bifido_class(khula_pci_mbiome_data) ],
    # color = khula_pci_mbiome_data.Bifidobacterium_breve,
    colormap = :viridis
)

lgC = Legend(
    C_Subfig[2,1],
    [
        MarkerElement(; marker=:circle, color=bifido_cdict["dom_breve"], strokewidth=1),
        MarkerElement(; marker=:circle, color=bifido_cdict["more_breve"], strokewidth=1),
        MarkerElement(; marker=:circle, color=bifido_cdict["more_longum"], strokewidth=1),
        MarkerElement(; marker=:circle, color=bifido_cdict["dom_longum"], strokewidth=1),
        MarkerElement(; marker=:circle, color=bifido_cdict["null"], strokewidth=1)
    ],
    [
        rich("Dominated (>70%) by ", rich("Bifidobacterium breve", font = :italic)),
        rich(rich("B. breve", font = :italic), " above 5% and ", rich("B. breve", font = :italic), "-to-", rich("B. longum", font = :italic), " ratio > 5.0" ),
        rich(rich("B. longum", font = :italic), " above 5% and ", rich("B. longum", font = :italic), "-to-", rich("B. breve", font = :italic), " ratio > 5.0" ),
        rich("Dominated (>70%) by ", rich("Bifidobacterium longum", font = :italic)),
        rich("Not characterized by ", rich("B. breve", font = :italic), " or ", rich("B. longum", font = :italic))
    ];
    tellheight = false,
    tellwidth = false,
    nbanks = 1
)

## Panel C: Species PCoA

axD = Axis(
    D_Subfig[1,1],
    xlabel = "MDS1 ("*string(round(100*spe_MDS_variances[1]; digits = 2))*"%)",
    ylabel = "MDS2 ("*string(round(100*spe_MDS_variances[2]; digits = 2))*"%)",
    title = "Taxonomic Profiles (species)",
    aspect = 1.0
)
hidedecorations!(
    axD;
    label = false,
    ticklabels = false,
    ticks = false,
    grid = true,
    minorgrid = false,
    minorticks = false
)
scD = scatter!(
    axD,
    spe_MDS_columns[:,1],
    spe_MDS_columns[:,2];
    markersize = 12,
    color = khula_pci_mbiome_data.InfantVisAtt,
    colormap = :viridis
)

cbD = Colorbar(Right_Col[3,:], scD; label = "Visual Orienting Behavior", vertical = false)

## Layout

Label(Left_Col[1, 1, TopLeft()], "a", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(B_Subfig[1, 1, TopLeft()], "b", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(C_Subfig[1, 1, TopLeft()], "c", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(D_Subfig[1, 1, TopLeft()], "d", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())

rowsize!(Left_Col, 1, Relative(0.35))
rowsize!(Left_Col, 2, Relative(0.65))
rowsize!(Right_Col, 1, Relative(0.48))
rowsize!(Right_Col, 2, Relative(0.48))
rowsize!(Right_Col, 3, Relative(0.04))
rowgap!(Right_Col, 10)

rowsize!(C_Subfig, 1, Relative(0.70))
rowsize!(C_Subfig, 2, Relative(0.30))
rowgap!(C_Subfig, 10)

colgap!(fig.layout, 2)

fig

save(joinpath(Base.pwd(),"manuscript", "figures", "Figure1.png"), fig)
save(joinpath(Base.pwd(),"manuscript", "figures", "Figure1.eps"), fig)
save(joinpath(Base.pwd(),"manuscript", "figures", "Figure1.pdf"), fig)

#####
# Exporting Supplementary Figures
#####

include(joinpath(Base.pwd(), "manuscript", "FigureS1.jl"))
include(joinpath(Base.pwd(), "manuscript", "FigureS2.jl"))
include(joinpath(Base.pwd(), "manuscript", "FigureS3.jl"))