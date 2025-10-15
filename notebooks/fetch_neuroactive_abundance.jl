#####
# Getting the neuraoactive genes
#####

using DataFrames
using CSV
using ProgressMeter
using Statistics
using CairoMakie

#####
# Part 1: Getting the gene-species table
#####
gene_spec_table = CSV.read(joinpath(Base.pwd(), "data", "long_stratified_gfs.csv"), DataFrame)

# Transform the "spec_raw" column by splitting each string at the dot,
# extracting the first two parts into new columns.
gene_spec_table.gen = getindex.(split.(gene_spec_table.spec_raw, "."), 1)
gene_spec_table.spec = getindex.(split.(gene_spec_table.spec_raw, "."), 2)
# Optionally drop the raw column because it is no longer needed:
select!(gene_spec_table, Not(:spec_raw))

# Normalizing abundance values to sum to 1 per sample
persample_grouped = groupby(gene_spec_table, :sample)
for g in persample_grouped
    total_ab = sum(g.ab)
    g.ab .= g.ab ./ total_ab
end
gene_spec_table_normalized = vcat(persample_grouped...)
# Now, gene_spec_table_normalized has normalized abundance values per sample.
# We can now safely compute the mean abundance of each gene per genus
# HOWEVER, If we simply perform
# DataFrames.combine(groupby(gene_spec_table_normalized, [:gene, :gen]), :ab => x -> mean(x))
# We will get only the gene abundances _when present_.
# i.e. if a single sample has `UniRef90_A0A376TR46 g__Escherichia``,
# the mean gene abundance would be the abundance in that sample only.
# What we should do instead is divide the sample-wise relabs
# by the total number of samples and then sum them.
gene_spec_table_normalized.ab = gene_spec_table_normalized.ab ./ length(unique(gene_spec_table.sample))

abs_grouped_gene_genus = groupby(gene_spec_table_normalized, [:gene, :gen])

abs_grouped_gene_species = groupby(gene_spec_table_normalized, [:gene, :spec])

## ATTENTION-TOGGLE!
## OPTION 1: Counting ABUNDANCE of the genes by carrier
gene_genus_ab_table = DataFrames.combine(abs_grouped_gene_genus, :ab => ( x -> sum(x) ) => :ab)
gene_species_ab_table = DataFrames.combine(abs_grouped_gene_species, :ab => ( x -> sum(x) ) => :ab)
## OPTION 2: Counting PRESENCE of the genes
# gene_genus_prev_table = DataFrames.combine(abs_grouped_gene_genus, :ab => ( x -> length(x) ) => :ab)
# gene_species_prev_table = DataFrames.combine(abs_grouped_gene_species, :ab => ( x -> length(x) ) => :ab)

# Now that all wrangling was performed with the whole profile,
# We can finally remove the functions filtered out by prevalence

passed_filter = Set(featurenames(unstratified_unirefs_filtered))
retain_genus_rows = [ el ∈ passed_filter for el in gene_genus_ab_table.gene ]
retain_species_rows = [ el ∈ passed_filter for el in gene_species_ab_table.gene ]

gene_genus_ab_table = gene_genus_ab_table[retain_genus_rows, :]
gene_species_ab_table = gene_species_ab_table[retain_species_rows, :]

#####
# Part 2: gathering the neuroactive genes and respective genesets
#####

# Get the list of neuroactive KOs
neuroactivepath = datafiles("gbm.txt")
map_ko_uniref_path = datafiles("map_ko_uniref90.txt.gz")
neuroactivekos = Leap.get_neuroactive_kos(neuroactivepath; consolidate = true)

# Map all the KOs to UNIREFS
kos2uniref = Dictionary{String, Vector{String}}()
for line in eachline(GzipDecompressorStream(open(map_ko_uniref_path)))
    line = split(line, '\t')
    insert!(kos2uniref, line[1], map(x-> String(match(r"UniRef90_(\w+)", x).captures[1]), line[2:end]))
end

# Convert the KOs into UNIREFS for the first list. Product is a Dict of pathways and unirefs.
neuroactive_unirefs = Dictionary{String, Vector{String}}()
for na in keys(neuroactivekos)
    searchfor = Iterators.flatten([kos2uniref[ko] for ko in neuroactivekos[na] if ko in keys(kos2uniref)]) |> Set |> (x -> string.(x))
    insert!(neuroactive_unirefs, na, searchfor)
end

# Collecting the UniRefs of interest in a single dataframe
gene_geneset_lines = Set{Tuple{String, String}}()
for interesting_geneset in
    [
        "Acetate synthesis",
        "Butyrate synthesis",
        "Menaquinone synthesis",
        "Quinolinic acid degradation",
        "GABA synthesis",
        "Glutamate synthesis",
        "Tryptophan synthesis",
    ]

    @showprogress for corresponding_uniref in neuroactive_unirefs[interesting_geneset]
        if ("UniRef90_"*corresponding_uniref ∈ gene_genus_ab_table.gene)
            push!(gene_geneset_lines, (interesting_geneset, "UniRef90_"*corresponding_uniref))
        end
    end

end

gene_geneset_table = DataFrame(collect(gene_geneset_lines), ["geneset", "gene"])

geneset_genes_grouped = groupby(gene_geneset_table, [:geneset, :gene])

## GENUS version:
joined_genus_abs = DataFrame[]
@showprogress for sg in geneset_genes_grouped
    push!(joined_genus_abs, innerjoin(sg, gene_genus_ab_table, on = :gene))
end
geneset_genusab_table = vcat(joined_genus_abs...)

# Normalize geneset abs to 1.0
geneset_genusabnormalized_table = groupby(geneset_genusab_table, :geneset)
for g in geneset_genusabnormalized_table
    sum_ab = sum(g.ab)
    g.proprotional_ab .= g.ab ./ sum_ab
end

geneset_genusabnormalized_table = vcat(geneset_genusabnormalized_table...)
plot_geneset_genusab_table = deepcopy(geneset_genusabnormalized_table)
final_geneset_genusabnormalized_table = DataFrames.combine(groupby(geneset_genusabnormalized_table, [:geneset, :gen]), :ab => ( x -> sum(x) ) => :total_ab)

sort!(final_geneset_genusabnormalized_table, [:geneset, :gen])
sort!(final_geneset_genusabnormalized_table, :total_ab; rev = true)

CSV.write(joinpath(Base.pwd(), "manuscript", "FSEA", "carrier_abundance_genesets_genus.csv"), final_geneset_genusabnormalized_table)

## SPECIES version:
joined_species_abs = DataFrame[]
@showprogress for sg in geneset_genes_grouped
    push!(joined_species_abs, innerjoin(sg, gene_species_ab_table, on = :gene))
end
geneset_speciesab_table = vcat(joined_species_abs...)

# Normalize geneset abs to 1.0
geneset_speciesabnormalized_table = groupby(geneset_speciesab_table, :geneset)
for g in geneset_speciesabnormalized_table
    sum_ab = sum(g.ab)
    g.proprotional_ab .= g.ab ./ sum_ab
end

geneset_speciesabnormalized_table = vcat(geneset_speciesabnormalized_table...)
plot_geneset_speciesab_table = deepcopy(geneset_speciesabnormalized_table)
final_geneset_speciesabnormalized_table = DataFrames.combine(groupby(geneset_speciesabnormalized_table, [:geneset, :spec]), :proprotional_ab => ( x -> sum(x) ) => :proportional_ab)

sort!(final_geneset_speciesabnormalized_table, [:geneset, :spec])
sort!(final_geneset_speciesabnormalized_table, :proportional_ab; rev = true)

CSV.write(joinpath(Base.pwd(), "manuscript", "FSEA", "carrier_abundance_genesets_species.csv"), final_geneset_speciesabnormalized_table)

#####
# Plotting the GENUS version as a heatmap, only aggregated genesets
#####
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

mat = Matrix(plot_df_wide[:, Not([:geneset])])
species = names(plot_df_wide, Not([:geneset]))
species = vcat( [ rich(el, font = :italic) for el in species[1:(end-1)] ], [ rich(species[end]) ] )
ylabs = plot_df_wide.geneset

# Plot
combined_gen_spec_fig =  Figure(size = (1200, 600))
Left_Subfig = GridLayout(combined_gen_spec_fig[1,1]; alignmode = Inside())
Right_Subfig = GridLayout(combined_gen_spec_fig[1,2]; alignmode = Inside())

ax_genus_aggregate = Axis(
    Left_Subfig[1, 1],
    xticklabelrotation = pi/5,
    yticks = (1:length(ylabs), ylabs),
    xticks = (1:length(species), species),
    xlabel = "Species",
    ylabel = "Gene set",
    title = "Gene set carried abundance across select genera",
    yreversed = true,
    xreversed = false
)

hm_genus_aggregate = heatmap!(ax_genus_aggregate, mat'; colormap = :viridis)

Colorbar(Left_Subfig[1,2], hm_genus_aggregate; label = "Normalized carried abundance")

#####
# Plotting the SPECIES version as a heatmap, only aggregated genesets
#####
plot_geneset_speciesab_sub = deepcopy(plot_geneset_speciesab_table)

# Define key species you want to show explicitly
focus_species = [
    "s__Bifidobacterium_longum",
    "s__Bifidobacterium_breve",
    "s__Bifidobacterium_bifidum",
    "s__Bifidobacterium_kashiwanohense",
    "s__Bifidobacterium_pseudocatenulatum",
    "s__Klebsiella_pneumoniae",
    "s__Klebsiella_michiganensis",
    "s__Escherichia_coli",
    "s__Bacteroides_vulgatus",
    "s__Bacteroides_fragilis"
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

# Sort rows by geneset/gene and columns alphabetically (so "Other" will come last)
# sort!(plot_df_wide, [:geneset, :gene])
# species_cols = sort(names(plot_df_wide, Not([:geneset, :gene])))
# plot_df_wide = plot_df_wide[:, vcat(["geneset", "gene"], species_cols)]
sort!(plot_df_wide, :geneset)
species_cols = sort(names(plot_df_wide, Not(:geneset, :Other)))
plot_df_wide = plot_df_wide[:, vcat(["geneset"], species_cols, ["Other"])]

# Matrix + labels
# mat = Matrix(plot_df_wide[:, Not([:geneset, :gene])])
# species = names(plot_df_wide, Not([:geneset, :gene]))
# genes = plot_df_wide.gene
# genesets = plot_df_wide.geneset
# ylabs = ["$(genesets[i]) | $(genes[i])" for i in eachindex(genes)]

mat = Matrix(plot_df_wide[:, Not([:geneset])])
species = names(plot_df_wide, Not([:geneset]))
species = [ replace(el, "_" => " ") for el in species ]
species = vcat( [ rich(el, font = :italic) for el in species[1:(end-1)] ], [ rich(species[end]) ] )
ylabs = plot_df_wide.geneset

# Plot
ax_species_aggregate = Axis(
    Right_Subfig[1, 1],
    xticklabelrotation = pi/5,
    yticks = (1:length(ylabs), ylabs),
    xticks = (1:length(species), species),
    xlabel = "Species",
    ylabel = "Gene set",
    title = "Gene set carried abundance across select species",
    yreversed = true,
    xreversed = false
)

hm_species_aggregate = heatmap!(ax_species_aggregate, mat'; colormap = :viridis)

Colorbar(Right_Subfig[1,2], hm_species_aggregate; label = "Normalized carried abundance")

colsize!(combined_gen_spec_fig.layout, 1, Relative(0.33))

combined_gen_spec_fig

### Attic/Graveyard
# #####
# # What is B3DS10 about?
# #####

# unihead = DataFrame(get(unstratified_unirefs_filtered))
# unihead.mygene = vec(abundances(unstratified_unirefs_filtered["UniRef90_B3DS10", :]))
# cor(unihead.mygene, khula_pci_mbiome_data.Bifidobacterium_longum)

# #####
# # Are Bbreve and Blongum carrying differental amounts of genesets?
# #####
# df_bif = filter(:spec => s -> occursin(r"s__Bifidobacterium_(longum|breve)", s), plot_geneset_speciesab_sub)

# # Keep only neuroactive genesets
# neuroactive_sets = unique(df_bif.geneset[occursin.("neuroactive", df_bif.geneset)])  # or define your own list
# df_bif = filter(:geneset => ∈(neuroactive_sets), df_bif)

# # Aggregate total abundance per geneset × species
# df_sum = combine(groupby(df_bif, [:geneset, :spec]), :ab => sum => :total_abundance)

# # Pivot species to columns (longum and breve)
# df_wide = unstack(df_sum, :spec, :total_abundance)

# # Replace missings with 0 (some genesets might not be carried by one of them)
# foreach(c -> coalesce!(df_wide[!, c], 0.0), names(df_wide, Not(:geneset)))

# # Optional: compute relative contribution of each species
# df_wide.total = coalesce.(df_wide.s__Bifidobacterium_longum, 0.0) .+
#                 coalesce.(df_wide.s__Bifidobacterium_breve, 0.0)
# df_wide.longum_frac = df_wide.s__Bifidobacterium_longum ./ df_wide.total
# df_wide.breve_frac = df_wide.s__Bifidobacterium_breve ./ df_wide.total

# # Clean up any NaNs (genesets with 0 total)
# replace!(df_wide.longum_frac, NaN => 0.0)
# replace!(df_wide.breve_frac, NaN => 0.0)

# # Sort by whichever species dominates
# sort!(df_wide, :longum_frac, rev=true)

# # Preview
# first(df_wide, 10)

# #####
# # Bonus: proportional carrier count
# #####

# geneset_geneprev_table = vcat(joined_abs...)
# rename!(geneset_geneprev_table, :ab => :count)
# geneset_geneprev_table = DataFrames.combine(groupby(geneset_geneprev_table, [:geneset, :gen]), :count => ( x -> length(x) ) => :total_count)
# # Normalize geneset prevs to 1.0
# geneset_geneprevnormalized_table = groupby(geneset_geneprev_table, :geneset)
# for g in geneset_geneprevnormalized_table
#     sum_count = sum(g.total_count)
#     g.total_count .= g.total_count ./ sum_count
# end
# final_geneset_geneprevnormalized_table = vcat(geneset_geneprevnormalized_table...)
# sort!(final_geneset_geneprevnormalized_table, :total_count; rev = true)
