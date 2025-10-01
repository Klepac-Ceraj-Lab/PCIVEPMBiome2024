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

khula_rf_data = deepcopy(khula_pci_mbiome_data)
khula_rf_data = filter_prevalence(khula_rf_data, 0.05)
khula_rf_data.InfantVisAtt = Leap.rangenormalize(khula_rf_data.InfantVisAtt)
select!(khula_rf_data, mdata_cols, :)
select!(khula_rf_data, :subject_id, :sample, :datasource, :pci_assess_age, :InfantVisAtt, :)

# regression_VOB_NoEnt_FullCV = Leap.probe_regression_randomforest( ## Quick version that performs regression with the final hyperparameter set. 
#     "regression_VOB_NoEnt_FullCV",
#     khula_rf_data,
#     identity,
#     collect(findfirst(names(khula_rf_data) .== "mbiome_sample_age"):ncol(khula_rf_data)),
#     :InfantVisAtt;
#     split_strat = "subject",
#     custom_input_group = nothing,
#     unique_col = :sample,
#     n_folds = 3,
#     n_replicas = 5,
#     n_rngs = 3,
#     tuning_space = (; #PRODUCTION
#         maxnodes_range = [ -1, 2, 4, 6, 8, 10 ],
#         nodesize_range = [ 3, 5, 7, 9, 11 ],
#         min_samples_split = [ 2, 3, 4 ],
#         sampsize_range = [ 0.7, 0.8, 0.9 ],
#         mtry_range = [ -1, 0, 10, 15, 20, 25, 30 ],
#         ntrees_range = [ 64, 128, 256 ]
#     )
# )
# sort(report_regression_merits(regression_VOB_NoEnt_FullCV), :Val_RMSE_mean)
# JLD2.@save joinpath(Base.pwd(), "manuscript", "models", "regression_VOB_NoEnt_FullCV_Results.jld") regression_VOB_NoEnt_FullCV

# regression_VOB_WtEnt_FullCV = Leap.probe_regression_randomforest( ## Quick version that performs regression with the final hyperparameter set. 
#     "regression_VOB_WtEnt_FullCV",
#     khula_rf_data,
#     identity,
#     collect(findfirst(names(khula_rf_data) .== "MaternalEntropy"):ncol(khula_rf_data)),
#     :InfantVisAtt;
#     split_strat = "subject",
#     custom_input_group = nothing,
#     unique_col = :sample,
#     n_folds = 3,
#     n_replicas = 5,
#     n_rngs = 3,
#     tuning_space = (; #PRODUCTION
#         maxnodes_range = [ -1, 2, 4, 6, 8, 10 ],
#         nodesize_range = [ 3, 5, 7, 9, 11 ],
#         min_samples_split = [ 2, 3, 4 ],
#         sampsize_range = [ 0.7, 0.8, 0.9 ],
#         mtry_range = [ -1, 0, 10, 15, 20, 25, 30 ],
#         ntrees_range = [ 64, 128, 256 ]
#     )
# )
# sort(report_regression_merits(regression_VOB_WtEnt_FullCV), :Val_RMSE_mean)
# JLD2.@save joinpath(Base.pwd(), "manuscript", "models", "regression_VOB_WtEnt_FullCV.jld") regression_VOB_WtEnt_FullCV

# hpimportances(regression_VOB_WtEnt_FullCV)
# khula_selection_rf = select(khula_rf_data, ["subject_id", "sample", "datasource", "InfantVisAtt", "MaternalEntropy", "mbiome_sample_age", "Bifidobacterium_breve", "Bifidobacterium_longum", "Streptococcus_salivarius", "alpha_shannon", "Enterococcus_faecalis", "Escherichia_coli", "Streptococcus_mitis", "Enterococcus_gallinarum", "Veillonella_atypica", "Clostridium_neonatale", "Bacteroides_vulgatus", "Bifidobacterium_bifidum", "Ruminococcus_gnavus", "Lactobacillus_reuteri", "Streptococcus_parasanguinis", "Bifidobacterium_kashiwanohense", "Enterococcus_avium", "Veillonella_dispar", "Enterococcus_faecium", "Klebsiella_variicola", "Collinsella_aerofaciens", "Klebsiella_pneumoniae", "Clostridioides_difficile", "Parabacteroides_distasonis", "Veillonella_parvula", "Actinomyces_sp_oral_taxon_181", "Flavonifractor_plautii", "Megasphaera_micronuciformis"])

# regression_VOB_SelectionNoEnt_FullCV = Leap.probe_regression_randomforest( ## Quick version that performs regression with the final hyperparameter set. 
#     "regression_VOB_SelectionNoEnt_FullCV",
#     khula_selection_rf,
#     identity,
#     collect(findfirst(names(khula_selection_rf) .== "mbiome_sample_age"):ncol(khula_selection_rf)),
#     :InfantVisAtt;
#     split_strat = "subject",
#     custom_input_group = nothing,
#     unique_col = :sample,
#     n_folds = 3,
#     n_replicas = 5,
#     n_rngs = 3,
#     tuning_space = (; #PRODUCTION
#         maxnodes_range = [ -1, 2, 4, 6, 8, 10 ],
#         nodesize_range = [ 3, 5, 7, 9, 11 ],
#         min_samples_split = [ 2, 3, 4 ],
#         sampsize_range = [ 0.7, 0.8, 0.9 ],
#         mtry_range = [ -1, 0, 10, 15, 20, 25 ],
#         ntrees_range = [ 64, 128, 256 ]
#     )
# )
# sort(report_regression_merits(regression_VOB_SelectionNoEnt_FullCV), :Val_RMSE_mean)
# JLD2.@save joinpath(Base.pwd(), "manuscript", "models", "regression_VOB_SelectionNoEnt_FullCV.jld") regression_VOB_SelectionNoEnt_FullCV

# regression_VOB_SelectionWtEnt_FullCV = Leap.probe_regression_randomforest( ## Quick version that performs regression with the final hyperparameter set. 
#     "regression_VOB_SelectionWtEnt_FullCV",
#     khula_selection_rf,
#     identity,
#     collect(findfirst(names(khula_selection_rf) .== "MaternalEntropy"):ncol(khula_selection_rf)),
#     :InfantVisAtt;
#     split_strat = "subject",
#     custom_input_group = nothing,
#     unique_col = :sample,
#     n_folds = 3,
#     n_replicas = 5,
#     n_rngs = 3,
#     tuning_space = (; #PRODUCTION
#         maxnodes_range = [ -1, 2, 4, 6, 8, 10 ],
#         nodesize_range = [ 3, 5, 7, 9, 11 ],
#         min_samples_split = [ 2, 3, 4 ],
#         sampsize_range = [ 0.7, 0.8, 0.9 ],
#         mtry_range = [ -1, 0, 10, 15, 20, 25, 30 ],
#         ntrees_range = [ 64, 128, 256 ]
#     )
# )
# sort(report_regression_merits(regression_VOB_SelectionWtEnt_FullCV), :Val_RMSE_mean)
# JLD2.@save joinpath(Base.pwd(), "manuscript", "models", "regression_VOB_SelectionWtEnt_FullCV.jld") regression_VOB_SelectionWtEnt_FullCV

JLD2.@load joinpath(Base.pwd(), "manuscript", "models", "regression_VOB_NoEnt_FullCV_Results.jld") regression_VOB_NoEnt_FullCV
JLD2.@load joinpath(Base.pwd(), "manuscript", "models", "regression_VOB_WtEnt_FullCV_Results.jld") regression_VOB_WtEnt_FullCV

#####
# Sourcing LMs and FSEA analysis from separate file
#####
# include(joinpath(Base.pwd(), "notebooks", "maaslin3_analysis.jl"))
# include(joinpath(Base.pwd(), "notebooks", "gene_glm_analysys.jl"))
# include(joinpath(Base.pwd(), "notebooks", "fsea_analysis.jl"))

#####
# Figure plotting block
#####

fig = Figure(; size = (1200, 1200))

ABC_Subfig = GridLayout(fig[1,:]; alignmode = Outside())
DE_Subfig = GridLayout(fig[2,:]; alignmode = Outside())
F_Subfig = GridLayout(fig[3,:]; alignmode = Outside())

AB_Subfig = GridLayout(ABC_Subfig[:,1]; alignmode = Inside())
C_Subfig = GridLayout(ABC_Subfig[:,2]; alignmode = Inside())

AB_Subfig_plots = GridLayout(AB_Subfig[1,:]; alignmode = Inside())
AB_Subfig_legends = GridLayout(AB_Subfig[2,:]; alignmode = Inside())

C_Subfig_plots = GridLayout(C_Subfig[1,:]; alignmode = Inside())
C_Subfig_legends = GridLayout(C_Subfig[2,:]; alignmode = Inside())

D_Subfig = GridLayout(DE_Subfig[1,:]; alignmode = Outside())
E_Subfig = GridLayout(DE_Subfig[2,:]; alignmode = Outside())

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
    nbanks = 1
)

aux_annot = annotations(
    plot_comparative_df.variable,
    Point.(collect(zip(log.(plot_comparative_df.pvalue).*(-1), plot_comparative_df.weightedImportance))),
    color = point_colors, fontsize = 6
)

# Middle part - FSEA plots

datadir = "manuscript/FSEA/"

## Middle-top (Unique effects)

res = NamedTuple(map(("entropy", "visual")) do mod
    lms = CSV.read(joinpath(datadir, "lms_Model2_$mod.csv"), DataFrame)
    Symbol(mod) => (;
        lms = lms,
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

    mdz = median(res[k].lms.z)
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

    if i == 1
        ax.yticks = (1:length(sigs), sigs)
    end

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
        lms = lms,
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

    if i == 1
        ax.yticks = (1:length(sigs), sigs)
    end

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

organisms_df = CSV.read("manuscript/FSEA/exports_neuroactives/unirefannotations.csv", DataFrame; stringtype = String)
# genus_level_data = select(khula_pre_data_genus, Not(mdata_cols[mdata_cols ∈ Ref(names(khula_pre_data_genus))]))

# genus_level_abundances = sort(DataFrame(:gen => names(genus_level_data), :ab => mean.(eachcol(genus_level_data))), :ab)
# genus_level_abundances.cumab = cumsum(genus_level_abundances.ab)

counts_df = DataFrames.combine(
    groupby(organisms_df, [:geneset, :organism]),
    nrow => :count,
    :organism => (x -> map(y -> replace(y, "g__" => ""), x)) => :clean_organism
)

counts_df = DataFrames.combine(
    groupby(counts_df, [:geneset]),
    :clean_organism => identity => :organism,
    :count => ( x -> x ./ sum(x) ) => :proportion
)

choice_genera = Dict(
    "Bifidobacterium" => (stack_pos = 1, bug_color = (:blue, 0.9)),
    "Bacteroides" => (stack_pos = 2, bug_color = (:orange, 0.9)),
    "Escherichia" => (stack_pos = 3, bug_color = (:red, 0.9)),
    "Klebsiella" => (stack_pos = 4, bug_color = (:green, 0.9)),
    # "Veillonella" => (stack_pos = 5, bug_color = (:purple, 0.9)),
    "Enterococcus" => (stack_pos = 6, bug_color = (:yellow, 0.9))
)

axis_names = Dict(
    "Acetate synthesis" => "Acetate\nsynthesis",
    "Butyrate synthesis" => "Butyrate\nsynthesis",
    "Menaquinone synthesis" => "Menaquinone\nsynthesis",
    "Quinolinic acid degradation" => "Quinolinic acid\ndegradation",
    "GABA synthesis" => "GABA\nsynthesis",
    "Glutamate synthesis" => "Glutamate\nsynthesis",
    "Tryptophan synthesis" => "Tryptophan\nsynthesis"
)

counts_df[:, :grouped_organism] .= ifelse.(counts_df.organism .∈ Ref(keys(choice_genera)), counts_df.organism, "Other")
counts_df = DataFrames.combine(groupby(counts_df, [:geneset, :grouped_organism]), :proportion => sum => :proportion)
counts_df = DataFrames.combine(groupby(counts_df, [:geneset]), groupindices => :xs, :proportion => identity => :proportion, :grouped_organism => identity => :grouped_organism)

push!(choice_genera, "Other" => (stack_pos = 12, bug_color = (:gray, 0.4)))

axE = Axis(
    F_Subfig[1,1],
    xticks = (
        collect(1:length(unique(counts_df.xs))),
        [ axis_names[subset(counts_df, :xs => (x -> x .== el)).geneset[1]] for el in 1:length(keys(axis_names)) ]
    ),
    ylabel = "Genus proportion/contribution"
)
hidexdecorations!(axE; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = true, minorticks = true )
hideydecorations!(axE; label = false, ticklabels = false, ticks = false, grid = true, minorgrid = true, minorticks = true )

bpEg = barplot!(
    axE,
    counts_df.xs,# .+ 1.2,
    counts_df.proportion,
    stack = [ choice_genera[el][:stack_pos] for el in counts_df.grouped_organism ],
    color = [ choice_genera[el][:bug_color] for el in counts_df.grouped_organism ],
    gap = 0.2
)

Legend(
    F_Subfig[1, 2],
    [
        PolyElement(; color=choice_genera["Bifidobacterium"].bug_color, strokewidth=1),
        PolyElement(; color=choice_genera["Bacteroides"].bug_color, strokewidth=1),
        PolyElement(; color=choice_genera["Escherichia"].bug_color, strokewidth=1),
        PolyElement(; color=choice_genera["Klebsiella"].bug_color, strokewidth=1),
        # PolyElement(; color=choice_genera["Veillonella"].bug_color, strokewidth=1),
        PolyElement(; color=choice_genera["Enterococcus"].bug_color, strokewidth=1),
        PolyElement(; color=choice_genera["Other"].bug_color, strokewidth=1)
    ],
    [
        "Bifidobacterium",
        "Bacteroides",
        "Escherichia",
        "Klebsiella",
        # "Veillonella",
        "Enterococcus",
        "Other"
    ];
    orientation = :vertical,
    tellheight = false,
    tellwidth = false,
    nbanks = 1
)

## Fixing global layout and labeling panels

colsize!(ABC_Subfig, 1, Relative(0.5))
colsize!(ABC_Subfig, 2, Relative(0.5))

rowsize!(AB_Subfig, 1, Relative(0.8))
rowsize!(AB_Subfig, 2, Relative(0.2))

rowsize!(C_Subfig, 1, Relative(0.8))
rowsize!(C_Subfig, 2, Relative(0.2))

colsize!(C_Subfig_plots, 1, Relative(0.1))
colsize!(C_Subfig_plots, 2, Relative(0.9))

colsize!(DE_Subfig, 0, Relative(0.01))
colsize!(DE_Subfig, 1, Relative(0.33))
colsize!(DE_Subfig, 2, Relative(0.33))
colsize!(DE_Subfig, 3, Relative(0.33))
rowsize!(DE_Subfig, 1, Relative(0.4))
rowsize!(DE_Subfig, 2, Relative(0.6))

colsize!(F_Subfig, 1, Relative(0.80))
colsize!(F_Subfig, 2, Relative(0.20))

rowsize!(fig.layout, 1, Relative(0.35))
rowsize!(fig.layout, 2, Relative(0.40))
rowsize!(fig.layout, 3, Relative(0.25))

Label(ABC_Subfig[1, 1, TopLeft()], "a", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(ABC_Subfig[1, 2, TopLeft()], "b", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(ABC_Subfig[1, 3, TopLeft()], "c", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(DE_Subfig[1, 1, TopLeft()], "d", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(DE_Subfig[2, 1, TopLeft()], "e", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())
Label(F_Subfig[1, 1, TopLeft()], "f", fontsize = 22, font = :bold, halign = :right, alignmode = Outside())

## Export Figure
save("manuscript/figures/Figure2.png", fig)
save("manuscript/figures/Figure2.eps", fig)
save("manuscript/figures/Figure2.pdf", fig)
