# #####
# # Getting the neuraoactive genes
# #####

# using HTTP
# using JSON
using ProgressMeter

genefamilies_samples = samplenames(unirefs_filtered_nonzeroentropy)
# Use a Set to store unique pairs for O(1) lookup
unique_pairs = Set{Tuple{String, String}}()

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
                first_field = split(line, "\t", limit=2)[1]
                parts = split(first_field, "|", limit=3)
                if length(parts) >= 2
                    pair = (parts[1], parts[2])
                    # The set automatically ignores duplicates.
                    if !(pair in unique_pairs)
                        #println("Adding $pair to the list!")
                        push!(unique_pairs, pair)
                    else
                        #println("$pair already in the list! Skipping...")
                    end
                end
            end
        end
    end
end

gene_spec_table = DataFrame(collect(unique_pairs), ["gene", "spec_raw"])
# Transform the "spec_raw" column by splitting each string at the dot,
# extracting the first two parts into new columns.
gene_spec_table.gen = getindex.(split.(gene_spec_table.spec_raw, "."), 1)
gene_spec_table.spec = getindex.(split.(gene_spec_table.spec_raw, "."), 2)
# Optionally drop the raw column if you no longer need it:
select!(gene_spec_table, Not(:spec_raw))

# Now build the dictionary by grouping the table by gene.
# For each gene, collect a vector of Dicts with keys "spec" and "gen".
gene_spec_dict = Dict{String, Vector{Dict{String, String}}}()

for group in groupby(gene_spec_table, :gene)
    gene = group.gene[1]  # every row in the group shares the same gene
    gene_spec_dict[gene] = [
        Dict("spec" => spec, "gen" => gen) for (spec, gen) in zip(group.spec, group.gen)
    ]
end

organism_matches = Vector{NamedTuple}()

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
        

    comm_features = featurenames(unirefs_filtered_nonzeroentropy)
    interest_features = comm_features[findall("UniRef90_".*neuroactive_unirefs[interesting_geneset] .∈ Ref(comm_features))]
    export_table = @chain comm2wide(unirefs_filtered_nonzeroentropy[interest_features, :]) begin
        select!([names(DataFrame(get(unirefs_filtered_nonzeroentropy)))..., interest_features...])
        sort(:subject_id)
    end

    open("manuscript/FSEA/exports_neuroactives/"*replace(interesting_geneset, " " => "_")*"_features.dat", "w") do fio
        for el in sort(interest_features)
            println(fio, el)
        end
    end

    for this_uniref_id in sort(interest_features)
        if this_uniref_id ∈ keys(gene_spec_dict)
            for el in gene_spec_dict[this_uniref_id]
                push!(organism_matches, (; geneset = interesting_geneset, organism = el["gen"]))
            end
        end
    end
end

organisms_df = @chain DataFrame(organism_matches) begin
    # unique( [:uniref_id, :organism] )
    sort([:geneset, :organism])
    # subset(:geneset => x -> x .== "Tryptophan synthesis")
end

CSV.write("manuscript/FSEA/exports_neuroactives/unirefannotations.csv", organisms_df)