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

Random.seed!(0)

#####
# 1. Loading data
#####
include(joinpath(Base.pwd(), "notebooks", "load_data.jl"))

#####
# Sourcing LMs and FSEA analysis from separate file
#####
include(joinpath(Base.pwd(), "notebooks", "rf_analysis.jl"))
include(joinpath(Base.pwd(), "notebooks", "maaslin3_analysis.jl"))
include(joinpath(Base.pwd(), "notebooks", "gene_glm_analysis.jl"))
include(joinpath(Base.pwd(), "notebooks", "fsea_analysis.jl"))
include(joinpath(Base.pwd(), "notebooks", "fetch_neuroactive_abundance.jl"))

#####
# Figure plotting block
#####

fig = Figure(; size = (1200, 1300))

ABC_Subfig = GridLayout(fig[1,:]; alignmode = Outside())
DE_Subfig = GridLayout(fig[2,:]; alignmode = Outside())
FG_Subfig = GridLayout(fig[3,:]; alignmode = Outside())

AB_Subfig = GridLayout(ABC_Subfig[:,1]; alignmode = Inside())
C_Subfig = GridLayout(ABC_Subfig[:,2]; alignmode = Inside())

AB_Subfig_plots = GridLayout(AB_Subfig[1,:]; alignmode = Inside())
AB_Subfig_legends = GridLayout(AB_Subfig[2,:]; alignmode = Inside())

C_Subfig_plots = GridLayout(C_Subfig[1,:]; alignmode = Inside())
C_Subfig_legends = GridLayout(C_Subfig[2,:]; alignmode = Inside())

D_Subfig = GridLayout(DE_Subfig[1,:]; alignmode = Inside())
E_Subfig = GridLayout(DE_Subfig[2,:]; alignmode = Inside())

F_Subfig = GridLayout(FG_Subfig[1,1]; alignmode = Inside())
G_Subfig = GridLayout(FG_Subfig[1,2:3]; alignmode = Inside())


## A+B - Abundance - Outcome scatterplots

axA = Axis(
    AB_Subfig_plots[1,1],
    xlabel = "arcsin-transformed abundance",
    ylabel = rich("normalized VOB x10", superscript("-5") ),
    yticks = ([0.0, 0.00005, 0.00010], ["0", "5", "10"]),
    title = rich("Bifidobacterium breve", font = :italic),
    aspect = 1.0
)
hidexdecorations!(axA; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = true, minorticks = true )
hideydecorations!(axA; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = true, minorticks = true )

scA = scatter!(
    axA,
    asin.(sqrt.(khula_pci_mbiome_data[:, :Bifidobacterium_breve])),
    khula_pci_mbiome_data.InfantVisAtt, ## Remember this is already autoscaled!
    color = map( x -> (x == 0 ? :red : :blue), khula_pci_mbiome_data[:, :Bifidobacterium_breve])
)

ax_B = Axis(
    AB_Subfig_plots[1,2],
    xlabel = "arcsin-transformed abundance",
    ylabel = rich("normalized VOB x10", superscript("-5") ),
    yticks = ([0.0, 0.00005, 0.00010], ["0", "5", "10"]),
    title = rich("Bifidobacterium longum", font = :italic),
    aspect = 1.0
)
hidexdecorations!(ax_B; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = true, minorticks = true )
hideydecorations!(ax_B; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = true, minorticks = true )

sc_B = scatter!(
    ax_B,
    asin.(sqrt.(khula_pci_mbiome_data[:, :Bifidobacterium_longum])),
    khula_pci_mbiome_data.InfantVisAtt, ## Remember this is already autoscaled!
    color = map( x -> (x == 0 ? :red : :blue), khula_pci_mbiome_data[:, :Bifidobacterium_longum])
)

Legend(
    AB_Subfig_legends[1, :],
    [
        MarkerElement(; marker=:circle, color=:red, strokewidth=1),
        MarkerElement(; marker=:circle, color=:blue, strokewidth=1),
    ],
    [
        "Absent from sample",
        "Present in sample"
    ];
    tellheight = false,
    tellwidth = false,
    nbanks = 2
)

## C - LM vs RF

function attribute_colors(lm_qvalue, rf_cumImp, plot_colorset = [:dodgerblue, :orange, :green, :purple])
    if (lm_qvalue & rf_cumImp)
        return plot_colorset[3]
    elseif (!lm_qvalue & rf_cumImp)
        return plot_colorset[1]
    elseif (lm_qvalue & !rf_cumImp)
        return plot_colorset[2]
    else
        return plot_colorset[4]
    end
end

### Comparison btw RF and Linear
ax_comparescatter = Axis(
    C_Subfig_plots[:,2],
    xlabel = "-log(p) for LM coefficients"
)

ax_comparescatter_aux = Axis(
    C_Subfig_plots[:,1],
    ylabel = "Relative Importance for RF"
)

linkyaxes!(ax_comparescatter, ax_comparescatter_aux)
hidexdecorations!(ax_comparescatter; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = true, minorticks = true )
hideydecorations!(ax_comparescatter)
hidexdecorations!(ax_comparescatter_aux)
hideydecorations!(ax_comparescatter_aux; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = true, minorticks = true )

t10 = ColorSchemes.tableau_10
cumulative_importance_threshold = 0.6
plot_colorset = [t10[1], t10[3], t10[7], (:white, 0.)]

## 1. calculate the importances and the cumulative sum, excluding the relevant variables
rf_model_importances = hpimportances(regression_VOB_NoEnt_FullCV, sort(report_regression_merits(regression_VOB_NoEnt_FullCV), :Val_RMSE_mean)[1, :Hyperpar_Idx]; change_hashnames=false)
rf_model_importances.cumulativeWeightedImportance = cumsum(rf_model_importances.weightedImportance)

## 2. get the data for linear models, from `Figure2-calculations.jl`
lm_coefs = select(sort(linear_lm_pvals[2], [:pvalue, :qvalue]), [:feature, :coef, :pvalue, :qvalue])

## 3. Join the data
plot_comparative_df = innerjoin(rf_model_importances, lm_coefs, on = [ :variable => :feature ])

## 4. Attribute color to each point
is_significative_lm = map(q-> q < 0.2 ? true : false, plot_comparative_df.qvalue)
is_over_threshold = map(ci-> ci <= cumulative_importance_threshold ? true : false, plot_comparative_df.cumulativeWeightedImportance)
point_colors = map((a, b) -> attribute_colors(a, b, plot_colorset), is_significative_lm, is_over_threshold)
# 5. Plot scatterplot
scatter!(
    ax_comparescatter,
    log.(plot_comparative_df.pvalue).*(-1), plot_comparative_df.weightedImportance,
    color = point_colors,
    marker = [ ( (el > 0) ? :utriangle : :dtriangle) for el in plot_comparative_df.coef],
    strokewidth=1, markersize = 12
)

not_in_lms = subset(outerjoin(rf_model_importances, lm_coefs, on = [ :variable => :feature ]), :coef => x -> ismissing.(x))
not_in_lms = subset(not_in_lms, :variable => x ->  x .!= "mbiome_sample_age")
not_in_lms = subset(not_in_lms, :variable => x ->  x .!= "alpha_shannon")

scatter!(
    ax_comparescatter_aux,
    ones(nrow(not_in_lms)),
    not_in_lms.weightedImportance,
    color = [ ( (el < cumulative_importance_threshold) ? plot_colorset[1] : plot_colorset[4]) for el in not_in_lms.cumulativeWeightedImportance],
    strokewidth=1, markersize = 12
)

Legend(
    C_Subfig_legends[1, 1],
    [
        MarkerElement(; marker=:dtriangle, color=:white, strokewidth=1),
        MarkerElement(; marker=:utriangle, color=:white, strokewidth=1),
        MarkerElement(; marker=:circle, color=:white, strokewidth=1),
    ],
    [
        "Negative ass. with VOB",
        "Positive ass. with VOB",
        "Removed from LMs\nby abundance filter"
    ];
    tellheight = false,
    tellwidth = false,
    margin = (0.0f0, 32.0f0, 0.0f0, 0.0f0),
    nbanks = 1
)

Legend(
    C_Subfig_legends[1, 2],
    [
        MarkerElement(; marker=:circle, color=plot_colorset[1], strokewidth=1),
        MarkerElement(; marker=:circle, color=plot_colorset[3], strokewidth=1),
        MarkerElement(; marker=:circle, color=plot_colorset[4], strokewidth=1),
    ],
    [
        "> 60% ranked importance in RF and\nq < 0.2 in LM",
        "> 60% ranked importance in RF",
        "Not important or significant"
    ];
    tellheight = false,
    tellwidth = false,
    margin = (0.0f0, 75.0f0, 0.0f0, 0.0f0),
    nbanks = 1
)

aux_annot = annotations(
    [ replace(el, "_" => "\n") for el in plot_comparative_df.variable ],
    Point.(collect(zip(log.(plot_comparative_df.pvalue).*(-1), plot_comparative_df.weightedImportance))),
    color = point_colors, fontsize = 8, rotation = pi/8, align = (:center, :center)
)
save("manuscript/figures/aux_annot.png", aux_annot)
# Middle part - FSEA plots

datadir = "manuscript/FSEA/"

## Middle-top (Unique effects)

res = NamedTuple(map(("entropy", "visual")) do mod
    lms = CSV.read(joinpath(datadir, "lms_Model2_$mod.csv"), DataFrame)
    Symbol(mod) => (;
        lms = subset(lms, :z => x -> .!isnan.(x)),  # remove missing z (these are the intercepts)
        na = Leap.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), lms.feature)),
        fsea = CSV.read(joinpath(datadir, "fsea_consolidated_Model2_$mod.csv"), DataFrame)
    )
end)

#-
sigs = unique(subset(vcat([res[k].fsea for k in 1:2]...), "qvalue"=> ByRow(<(0.2))).geneset)
sigs_y = Dict(gs=> i for (i, gs) in enumerate(sigs))
#-

prettytitles = Dict(
    :entropy => "Maternal Unpredictability",
    :visual => "Infant VOB"
)

for (i, k) in enumerate(keys(res))

    @show mdz = median(res[k].lms.z)
    lms = res[k].lms
    na = res[k].na
    fsea = res[k].fsea
    
    ax = Axis(
        DE_Subfig[1,i];
        xlabel = "z statistic",
        ylabel = "geneset",
        title="$(prettytitles[k])")

    xlims!(ax, (-3.5, 3.5))
    ylims!(ax, (0.5, 4.5))

    hidedecorations!(ax, label = false, ticklabels = false, ticks = false, minorticks = false)
    i == 1 && hideydecorations!(ax, ticklabels = false, ticks = false, minorticks = false)
    i != 1 && hideydecorations!(ax, ticks = false, minorticks = false)

    # if i == 1
        ax.yticks = (1:length(sigs), sigs)
    # end

    vlines!(ax, [mdz]; ymin = 0.0, ymax = 8.0, linestyle = :dash, color = :gray)

    for sig in sigs
        xs = na[sig]
        i = findfirst(==(sig), fsea.geneset)
        y = sigs_y[sig]
        med = median(lms.z[xs])
        color = fsea.qvalue[i] > 0.2 ? :gray :
                    # med < 0. ? :dodgerblue : :darkorange3
                    med < mdz ? :dodgerblue : :darkorange3
        violin!(ax, fill(y, length(xs)), lms.z[xs] ; orientation=:horizontal, color=(color, 0.3))
        scatter!(ax, lms.z[xs], y .+ rand(Normal(0., 0.1), length(xs)); color)
        lines!(ax, [med, med], y .+ [- 0.4, 0.4]; color)
    end
end

Legend(## Vertical version, to the side
    DE_Subfig[1,3],
    vcat(
        [
            MarkerElement(; color = c, marker = :circle) for c in (:dodgerblue, :darkorange3, :gray)
        ],
        [
            LineElement(; color = :gray, linestyle = :dash)
        ]
    ),
    [
        "q < 0.2, set median < effect median", 
        "q < 0.2, set median > effect median", 
        "not significant", 
        "effect median z"
    ],
    orientation = :vertical
)

## Middle-bottom (Interaction effects)

res = NamedTuple(map(("entropy", "visual", "product")) do mod
    lms = CSV.read(joinpath(datadir, "lms_Model3_$mod.csv"), DataFrame)
    Symbol(mod) => (;
        lms = subset(lms, :z => x -> .!isnan.(x)),  # remove missing z (these are the intercepts)
        na = Leap.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), lms.feature)),
        fsea = CSV.read(joinpath(datadir, "fsea_consolidated_Model3_$mod.csv"), DataFrame)
    )
end)

#-
sigs = unique(subset(vcat([res[k].fsea for k in 1:3]...), "qvalue"=> ByRow(<(0.2))).geneset)
sigs_y = Dict(gs=> i for (i, gs) in enumerate(sigs))
#-

prettytitles = Dict(
    :entropy => "Maternal Unpredictability",
    :visual => "Infant VOB",
    :product => "Interaction term"
)

for (i, k) in enumerate(keys(res))

    mdz = median(res[k].lms.z)
    lms = res[k].lms
    na = res[k].na
    fsea = res[k].fsea
    
    ax = Axis(
        DE_Subfig[2,i];
        xlabel = "z statistic",
        ylabel = "geneset",
        title="$(prettytitles[k])")

    xlims!(ax, (-3.5, 3.5))
    ylims!(ax, (0.5, 6.5))

    hidedecorations!(ax, label = false, ticklabels = false, ticks = false, minorticks = false)
    i == 1 && hideydecorations!(ax, ticklabels = false, ticks = false, minorticks = false)
    i != 1 && hideydecorations!(ax, ticks = false, minorticks = false)

    # if i == 1
        ax.yticks = (1:length(sigs), sigs)
    # end

    vlines!(ax, [mdz]; ymin = 0.0, ymax = 8.0, linestyle = :dash, color = :gray)

    for sig in sigs
        xs = na[sig]
        i = findfirst(==(sig), fsea.geneset)
        y = sigs_y[sig]
        med = median(lms.z[xs])
        color = fsea.qvalue[i] > 0.2 ? :gray :
                    # med < 0. ? :dodgerblue : :darkorange3
                    med < mdz ? :dodgerblue : :darkorange3
        violin!(ax, fill(y, length(xs)), lms.z[xs] ; orientation=:horizontal, color=(color, 0.3))
        scatter!(ax, lms.z[xs], y .+ rand(Normal(0., 0.1), length(xs)); color)
        lines!(ax, [med, med], y .+ [- 0.4, 0.4]; color)
    end
end

Label(DE_Subfig[1, 0, Left()], "Unique effects models", rotation = pi/2)
Label(DE_Subfig[2, 0, Left()], "Interaction models", rotation = pi/2)

## Organism barplot
ab_cmap = :PuBu
ab_cgrad = cgrad(ab_cmap, collect(0.0:0.2:0.8), categorical=false, rev=false)
## Plotting the GENUS version as a heatmap, only aggregated genesets
plot_geneset_genusab_sub = deepcopy(plot_geneset_genusab_table)

# Define key genera you want to show explicitly
focus_genera = [
    "g__Bifidobacterium",
    "g__Klebsiella",
    "g__Escherichia",
    "g__Bacteroides"
]

# Normalize the abundance *within each geneset* so that total abundance per geneset = 1
geneset_grouped = groupby(plot_geneset_genusab_sub, :geneset)
for g in geneset_grouped
    total = sum(g.ab)
    if total > 0
        g.ab .= g.ab ./ total
    end
end
plot_geneset_genusab_sub = vcat(geneset_grouped...)  # reassemble

# Collapse all genera not in your focus list into "Other"
plot_geneset_genusab_sub.gen = ifelse.(
    in.(plot_geneset_genusab_sub.gen, Ref(focus_genera)),
    plot_geneset_genusab_sub.gen,
    "Other"
)

# Aggregate abundance for duplicates (e.g. same gene/genus after pooling)
plot_geneset_genusab_sub = DataFrames.combine(
    groupby(plot_geneset_genusab_sub, [:geneset, :gen]),
    # groupby(plot_geneset_genusab_sub, [:geneset, :gene, :gen]),
    :ab => sum => :ab
)

# Transform and pivot to wide format
plot_geneset_genusab_sub.gen = [ replace(el, "g__" => "") for el in plot_geneset_genusab_sub.gen ]
plot_geneset_genusab_sub.ab_log = log10.((plot_geneset_genusab_sub.ab .+ 1e-8) * 1e8)
plot_df_wide = unstack(plot_geneset_genusab_sub, [ :geneset ], :gen, :ab)
# plot_df_wide = unstack(plot_geneset_genusab_sub, [:geneset, :gene], :gen, :ab)

# Fill missing values with zero
for i in 1:nrow(plot_df_wide)
    for j in 1:ncol(plot_df_wide)
        if ismissing(plot_df_wide[i, j])
            plot_df_wide[i, j] = 0.0
        end
    end
end

# Sort rows by geneset/gene and columns alphabetically (so "Other" will come last)
sort!(plot_df_wide, :geneset)
genus_cols = sort(names(plot_df_wide, Not(:geneset, :Other)))
plot_df_wide = plot_df_wide[:, vcat(["geneset"], genus_cols, ["Other"])]

gen_mat = Matrix(plot_df_wide[:, Not([:geneset])])
genera = names(plot_df_wide, Not([:geneset]))
genera = vcat( [ rich(el, font = :italic) for el in genera[1:(end-1)] ], [ rich(genera[end]) ] )
ylabs = plot_df_wide.geneset

# Plot
ax_genus_aggregate = Axis(
    F_Subfig[1, 1],
    xticklabelrotation = pi/5,
    yticks = (1:length(ylabs), ylabs),
    xticks = (1:length(genera), genera),
    xlabel = "Species",
    ylabel = "Gene set",
    title = "Gene set abundance\ncarried by select genera",
    yreversed = true,
    xreversed = false
)

hm_genus_aggregate = heatmap!(ax_genus_aggregate, gen_mat'; colormap = ab_cgrad, colorrange = [0.0, 0.75])

## Plotting the SPECIES version as a heatmap, only aggregated genesets
plot_geneset_speciesab_sub = deepcopy(plot_geneset_speciesab_table)

# Define key species you want to show explicitly
focus_species = [
    "s__Bifidobacterium_longum",
    "s__Bifidobacterium_breve",
    "s__Bifidobacterium_bifidum",
    "s__Klebsiella_pneumoniae",
    "s__Escherichia_coli",
    "s__Bacteroides_vulgatus",
]

# Normalize the abundance *within each geneset* so that total abundance per geneset = 1
geneset_grouped = groupby(plot_geneset_speciesab_sub, :geneset)
for g in geneset_grouped
    total = sum(g.ab)
    if total > 0
        g.ab .= g.ab ./ total
    end
end
plot_geneset_speciesab_sub = vcat(geneset_grouped...)  # reassemble

# Collapse all species not in your focus list into "Other"
plot_geneset_speciesab_sub.spec = ifelse.(
    in.(plot_geneset_speciesab_sub.spec, Ref(focus_species)),
    plot_geneset_speciesab_sub.spec,
    "Other"
)

# Aggregate abundance for duplicates (e.g. same gene/species after pooling)
plot_geneset_speciesab_sub = DataFrames.combine(
    groupby(plot_geneset_speciesab_sub, [:geneset, :spec]),
    # groupby(plot_geneset_speciesab_sub, [:geneset, :gene, :spec]),
    :ab => sum => :ab
)

# Transform and pivot to wide format
plot_geneset_speciesab_sub.spec = [ replace(el, "s__" => "") for el in plot_geneset_speciesab_sub.spec ]
plot_geneset_speciesab_sub.ab_log = log10.((plot_geneset_speciesab_sub.ab .+ 1e-8) * 1e8)
plot_df_wide = unstack(plot_geneset_speciesab_sub, [ :geneset ], :spec, :ab)
# plot_df_wide = unstack(plot_geneset_speciesab_sub, [:geneset, :gene], :spec, :ab)

# Fill missing values with zero
for i in 1:nrow(plot_df_wide)
    for j in 1:ncol(plot_df_wide)
        if ismissing(plot_df_wide[i, j])
            plot_df_wide[i, j] = 0.0
        end
    end
end

sort!(plot_df_wide, :geneset)
species_cols = sort(names(plot_df_wide, Not(:geneset, :Other)))
plot_df_wide = plot_df_wide[:, vcat(["geneset"], species_cols, ["Other"])]
spe_mat = Matrix(plot_df_wide[:, Not([:geneset])])
species = names(plot_df_wide, Not([:geneset]))
species = [ replace(el, "_" => " ") for el in species ]
species = vcat( [ rich(el, font = :italic) for el in species[1:(end-1)] ], [ rich(species[end]) ] )
ylabs = plot_df_wide.geneset

ax_species_aggregate = Axis(
    G_Subfig[1, 1],
    xticklabelrotation = pi/5,
    yticks = (1:length(ylabs), ylabs),
    xticks = (1:length(species), species),
    xlabel = "Species",
    title = "Gene set abundance\ncarried by select species",
    yreversed = true,
    xreversed = false
)

hm_species_aggregate = heatmap!(ax_species_aggregate, spe_mat'; colormap = ab_cgrad, colorrange = [0.0, 0.75])

Colorbar(G_Subfig[1,2], hm_genus_aggregate; label = "Normalized carried abundance")
# Colorbar(G_Subfig[1,2], hm_species_aggregate; label = "Normalized carried abundance")

linkyaxes!(ax_genus_aggregate, ax_species_aggregate)
hideydecorations!(ax_species_aggregate; label = true, ticklabels = true, ticks = false, grid = true, minorgrid = true, minorticks = false )
hideydecorations!(ax_genus_aggregate; label = true, ticklabels = false, ticks = false, grid = true, minorgrid = true, minorticks = false )

Label(FG_Subfig[1, 0, Left()], "Gene set", rotation = pi/2)
## Fixing global layout and labeling panels

colsize!(ABC_Subfig, 1, Relative(0.50))
colsize!(ABC_Subfig, 2, Relative(0.50))

rowsize!(AB_Subfig, 1, Relative(0.7))
rowsize!(AB_Subfig, 2, Relative(0.3))

colgap!(AB_Subfig, 3)

rowsize!(C_Subfig, 1, Relative(0.7))
rowsize!(C_Subfig, 2, Relative(0.3))

colsize!(C_Subfig_plots, 1, Relative(0.1))
colsize!(C_Subfig_plots, 2, Relative(0.9))

colsize!(DE_Subfig, 0, Relative(0.01))
colsize!(DE_Subfig, 1, Relative(0.33))
colsize!(DE_Subfig, 2, Relative(0.33))
colsize!(DE_Subfig, 3, Relative(0.33))
rowsize!(DE_Subfig, 1, Relative(0.4))
rowsize!(DE_Subfig, 2, Relative(0.6))

colsize!(FG_Subfig, 0, Relative(0.01))
colsize!(FG_Subfig, 1, Relative(0.40))
colsize!(FG_Subfig, 2, Relative(0.59))

rowsize!(fig.layout, 1, Relative(0.33))
rowsize!(fig.layout, 2, Relative(0.37))
rowsize!(fig.layout, 3, Relative(0.30))

Label(AB_Subfig_plots[1, 1, TopLeft()], "a", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(AB_Subfig_plots[1, 2, TopLeft()], "b", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(C_Subfig[1, 1, TopLeft()], "c", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(DE_Subfig[1, 1, TopLeft()], "d", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(DE_Subfig[2, 1, TopLeft()], "e", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(F_Subfig[1, 1, TopLeft()], "f", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(G_Subfig[1, 1, TopLeft()], "g", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())

fig

## Export Figure
save("manuscript/figures/Figure2.png", fig)
save("manuscript/figures/Figure2.eps", fig)
save("manuscript/figures/Figure2.pdf", fig)
