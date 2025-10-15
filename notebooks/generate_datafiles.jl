#####
# The purpose of this script is to read the source of truth locally on `vkclab-hopper`
# and generate the data files that will be shared to Dryad and
# loaded by `DataToolkit.jl` .toml data collection configuration files.
# It should not be sourced by individuals intending to reproduce the analysis, which should
# rather use the exporta herein generated.
#####

using CSV
using DataFrames
using Chain
using CategoricalArrays
using Microbiome
using Distances
using Dictionaries
using CodecZlib
using Leap
using PCIVEPMBiome2024
using Arrow
using SparseArrays
using StatsBase
using ProgressMeter

## Step 1: Derive a core set of samples that contain complete metadata

khula_visits = @chain CSV.read(joinpath(Base.pwd(), "ext", "khula_ages_visit_reltable.csv"), DataFrame; stringtype = String) begin
    select!(Not([:datacolor]))
    subset!(:visit => x -> x .== "3mo")
    subset!(:site => x -> x .== "ZAF")
end

## Samples excluded from the study
samps_to_delete = [ ] ## Left empty on purpose - moving filtering down to `load_data.jl` for transparency

khula_pci_mdata = @chain CSV.read("ext/PCI_HIV_MEDS_EPDS.csv", DataFrame; stringtype = String) begin
    select(:subject_id, :MaternalEntropy, :InfantVisAtt, :age_pci_days_1ReplSeriesMean, :MomEducationContinuous, :child_sex)
    dropmissing()
    transform!(:subject_id => (x -> "khula-" .* x) => :subject_id)
    transform!(:age_pci_days_1ReplSeriesMean => (x -> ( x ./ 30.5) ) => :age_pci_days_1ReplSeriesMean)
    innerjoin(khula_visits, _, on = :individualID => :subject_id)
    rename!([
        :specimenID => :sample,
        :individualID => :subject_id,
        :samplingAge => :mbiome_sample_age,
        :age_pci_days_1ReplSeriesMean => :pci_assess_age,
        :MomEducationContinuous => :MaternalEducation,
        :site => :datasource,
        # :epds_score_en_prenataldepression => :EPDS_dep_prenatal
    ])
    transform!([:pci_assess_age, :mbiome_sample_age] => ((x,y) -> abs.(x .- y)) => :age_diff; renamecols = false)
    subset!([:pci_assess_age, :mbiome_sample_age] => (x,y) -> abs.(x .- y) .< 1.0)
    # subset!([:pci_assess_age, :mbiome_sample_age] => (x,y) -> abs.(x .- y) .< 3.0)
    subset!(:MaternalEntropy => x -> x .>= 0.0)
    subset!(:sample => (x -> (x .∉ Ref(samps_to_delete))))
    dropmissing()
end

### On SEP/2024, the data for InfVisAttn had a small correction. The following code block deals with such alteration.
rename!(khula_pci_mdata, :InfantVisAtt => :OldInfantVisAtt)

new_infvisattn = @chain CSV.read("ext/BVisualAdjusted.csv", DataFrame; stringtype = String) begin
    dropmissing()
    transform!(:subject_id => (x -> replace.(x, "_" => "-")) => :subject_id)
    transform!(:subject_id => (x -> "khula-" .* x) => :subject_id)
    rename!(["NT_b_visual_active_counts_ordinal_U" => "VisualActiveCounts", "bVisual_Adjusted_CorrectedLength" => "InfantVisAtt"])
    subset(:VisualActiveCounts => x -> x .!= " " )
    transform!(:VisualActiveCounts => (x -> tryparse.(Float64, x) ) => :VisualActiveCounts)
    transform!(:InfantVisAtt => (x -> tryparse.(Float64, x) ) => :InfantVisAtt)
    dropmissing()
end

khula_pci_mdata = innerjoin(khula_pci_mdata, new_infvisattn, on = :subject_id)
transform!(khula_pci_mdata, [:MaternalEntropy, :InfantVisAtt] => ((x, y) -> ( Leap.rangenormalize(x) .* Leap.rangenormalize(y) ) ) => :Product_Single)
### End of handling InfVisAttn data alteration

sort(setdiff(khula_visits.individualID, khula_pci_mdata.subject_id))
sort(setdiff(khula_pci_mdata.subject_id, khula_visits.individualID))

mdata_cols =  ["subject_id", "sample", "datasource", "pci_assess_age", "MaternalEntropy", "InfantVisAtt", "mbiome_sample_age", "child_sex", "MaternalEducation"]

focus_mdata = select(
    khula_pci_mdata,
    mdata_cols
)

## Step 2: Using this core set of samples, collect the functional profiles matching the sample IDs from the full set.

function my_load(::UnirefProfiles, timepoint_metadata::DataFrame)
    comm = Leap.read_gfs_arrow(; arrow_path = "/brewster/guilherme/scratch/genefunctions/genefamilies_HM36.arrow")
    insert!(comm, timepoint_metadata; namecol=:sample)

    valid_samples = intersect(name.(comm.samples), timepoint_metadata.sample)
    comm = comm[:, valid_samples]

    if isempty(comm.samples)
        @warn "No overlapping samples between profile and metadata!"
    end

    return comm
end

unirefs = my_load(Leap.UnirefProfiles(), focus_mdata) # this can take a bit

## Step 3: lock in the included samples by the intersection of
# (1) full metadata and (2) functional profile available.

paper_sample_set = sort(samplenames(unirefs))

## Write the sample names on the "data" folder
# Note: the sample filter for the outlier S. mitis sample was moved to load_data so it could be included with the metadata.
# Hence, this line should not be run anymore - it ran once to export the core set of 93 samples.
# open(joinpath(Base.pwd(), "data", "samples.txt"), "w") do fio
#     for this_sample in paper_sample_set
#         println(fio, this_sample)
#     end
# end

## Intermission:
# Now that we have the full sample set, we can once again export the uniref.arrow file
# containing only the samples that we intend to make available with the manuscript!
## Run this only once every update to Humann!
#
# for (this_name, should_stratify) in [
#     ("genefamilies_HM36_unstratified.arrow", false)
#     ("genefamilies_HM36_stratified.arrow", true)
# ]
#     Leap.write_gfs_arrow(;
#         kind = "genefamilies",
#         file_glob_pattern = r"SEQ\d+_S\d+_genefamilies",
#         file_root = "/home/guilherme/thunderbay-mount/processed/internal/humann/main",
#         reads_root = "/home/guilherme/thunderbay-mount/processed/internal/kneaddata/read_counts.csv",
#         suffix_pattern = r"_S\d+.*$",
#         output_path = joinpath(Base.pwd(), "data"),
#         output_name = this_name,
#         names = false,
#         stratified = should_stratify,
#         sample_filter = paper_sample_set
#     )
# end

## Step 4: with the locked-in set, now we can load and save the master metadata table, and load nad save the taxonomic profiles...

locked_mdata = sort(subset(focus_mdata, :sample => x -> x .∈ Ref(paper_sample_set)), :sample)
CSV.write(joinpath(Base.pwd(), "data", "study_metadata.csv"), locked_mdata)

#####
# Taxonomic Profiles
#####
khula_tax_data_species = @chain Leap.load_raw_metaphlan(
    "/home/guilherme/Repos/PCIVEPMBiome2024/ext/profiles";
    replace_pattern = r"SEQ0\d+_S\d+_mpa_v31_CHOCOPhlAn_201901_profile",
    cleanup_pattern = "_mpa_v31_CHOCOPhlAn_201901_profile.tsv") begin
    filter(t-> !ismissing(taxrank(t)), _[:, samplenames(_)])
    filter(t-> taxrank(t) == :species, _[:, samplenames(_)])
    comm2wide()
    select(Not([:sample, :file]))
    subset!(:sample_base => x -> x .∈ Ref(paper_sample_set))
    dropmissing()
end

khula_tax_data_genera = @chain Leap.load_raw_metaphlan(
    "/home/guilherme/Repos/PCIVEPMBiome2024/ext/profiles";
    replace_pattern = r"SEQ0\d+_S\d+_mpa_v31_CHOCOPhlAn_201901_profile",
    cleanup_pattern = "_mpa_v31_CHOCOPhlAn_201901_profile.tsv") begin
    filter(t-> !ismissing(taxrank(t)), _[:, samplenames(_)])
    filter(t-> taxrank(t) == :genus, _[:, samplenames(_)])
    comm2wide()
    select(Not([:sample, :file]))
    subset!(:sample_base => x -> x .∈ Ref(paper_sample_set))
    dropmissing()
end

## Normalizing to 1.0 (after removing unknowns)
for i in 1:nrow(khula_tax_data_species)
    row_sum = sum(khula_tax_data_species[i,2:end])
    for j in 2:ncol(khula_tax_data_species)
        khula_tax_data_species[i,j] = khula_tax_data_species[i,j] / row_sum
    end
end

for i in 1:nrow(khula_tax_data_genera)
    row_sum = sum(khula_tax_data_genera[i,2:end])
    for j in 2:ncol(khula_tax_data_genera)
        khula_tax_data_genera[i,j] = khula_tax_data_genera[i,j] / row_sum
    end
end

CSV.write(joinpath(Base.pwd(), "data", "species_data_MP3.csv"), khula_tax_data_species)
CSV.write(joinpath(Base.pwd(), "data", "genus_data_MP3.csv"), khula_tax_data_genera)

#####
# Bonus: assembling and digesting the stratified abundances for Figure 2 barplots!
#####

unstratified_unirefs = load_custom_unirefs(UnirefProfiles(), study_mdata; custom_arrow_path = joinpath(Base.pwd(), "data", "genefamilies_HM36_unstratified.arrow"))
relativeabundance!(unstratified_unirefs)

genefamilies_samples = samplenames(unstratified_unirefs)
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
CSV.write(joinpath(Base.pwd(), "data", "long_stratified_gfs.csv"), gene_spec_table)
