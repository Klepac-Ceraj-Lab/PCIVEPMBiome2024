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
subset!(gene_spec_table, :sample => x -> x .∈ Ref(paper_sample_set))

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
neuroactivepath = joinpath(Base.pwd(), "data", "gbm.txt")
map_ko_uniref_path = joinpath(Base.pwd(), "data", "map_ko_uniref90.txt.gz")
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
        "Menaquinone synthesis",
        "Quinolinic acid degradation",
        "GABA synthesis",
        "Glutamate synthesis",
        "ClpB",
        "Tryptophan synthesis"
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

CSV.write(joinpath(Base.pwd(), "manuscript", "gene_glm", "carrier_abundance_genesets_species.csv"), final_geneset_speciesabnormalized_table)