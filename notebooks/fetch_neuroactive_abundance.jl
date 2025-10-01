
# #####
# # Getting the neuraoactive genes
# #####

# using HTTP
# using JSON
using ProgressMeter

#####
# Part 1: Getting the gene-species table
#####

genefamilies_samples = samplenames(unirefs_filtered_nonzeroentropy)
# Use a Set to store unique pairs for O(1) lookup
ab_lines = Set{Tuple{String, String, String, Float64}}()

@showprogress for this_genefamilies_sample in genefamilies_samples

    this_genefamilies_filepath = @chain readdir("/grace/sequencing/processed/mgx/humann/main") begin
        filter(x -> contains(x, this_genefamilies_sample), _)
        filter(x -> contains(x, r"genefamilies"), _)
    end

    # Open the file and iterate over each line.
    open("/grace/sequencing/processed/mgx/humann/main/" * first(this_genefamilies_filepath)) do file
        for line in eachline(file)
            if occursin("s__", line)  # Slightly faster alternative to contains for literal substrings.
                # Only split the line as much as needed.
                first_field, ab = split(line, "\t", limit=2)
                parts = split(first_field, "|", limit=3)
                if length(parts) >= 2
                    ll = (this_genefamilies_sample, parts[1], parts[2], parse(Float64, ab))
                    push!(ab_lines, ll)
                end
            end
        end
    end
end

gene_spec_table = DataFrame(collect(ab_lines), ["sample", "gene", "spec_raw", "ab"])
# Transform the "spec_raw" column by splitting each string at the dot,
# extracting the first two parts into new columns.
gene_spec_table.gen = getindex.(split.(gene_spec_table.spec_raw, "."), 1)
gene_spec_table.spec = getindex.(split.(gene_spec_table.spec_raw, "."), 2)
# Optionally drop the raw column because it is no longer needed:
select!(gene_spec_table, Not(:spec_raw))
# Normalizing abundance values to sum to 1 per sample
grouped = groupby(gene_spec_table, :sample)
for g in grouped
    total_ab = sum(g.ab)
    g.ab .= g.ab ./ total_ab
end
gene_spec_table = vcat(grouped...)
# Now, gene_spec_table has normalized abundance values per sample.
# We can now safely compute the mean abundance of each gene per genus
grouped = groupby(gene_spec_table, [:gene, :gen])
gene_genab_table = DataFrames.combine(grouped, :ab => x -> mean(x))

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
        if ("UniRef90_"*corresponding_uniref ∈ gene_genab_table.gene)
            push!(gene_geneset_lines, (interesting_geneset, "UniRef90_"*corresponding_uniref))
        end
    end

end

gene_geneset_table = DataFrame(collect(gene_geneset_lines), ["geneset", "gene"])

gg = groupby(gene_geneset_table, [:geneset, :gene])
joined_abs = DataFrame[]
@showprogress for sg in gg
    push!(joined_abs, innerjoin(sg, gene_spec_table, on = :gene))
end
geneset_geneab_table = vcat(joined_abs...)
geneset_geneab_table = DataFrames.combine(groupby(geneset_geneab_table, [:geneset, :gen]), :ab_function => ( x -> sum(x) ) => :total_ab)

# Normalize geneset abs to 1.0
geneset_geneabnormalized_table = groupby(geneset_geneab_table, :geneset)
for g in geneset_geneabnormalized_table
    sum_ab = sum(g.total_ab)
    g.total_ab .= g.total_ab ./ sum_ab
end
final_geneset_geneabnormalized_table = vcat(geneset_geneabnormalized_table...)
sort!(final_geneset_geneabnormalized_table, [:geneset, :gen])


organisms_df = @chain DataFrame(organism_matches) begin
    # unique( [:uniref_id, :organism] )
    sort([:geneset, :organism])
    # subset(:geneset => x -> x .== "Tryptophan synthesis")
end

CSV.write("manuscript/FSEA/exports_neuroactives/unirefabundances.csv", organisms_df)