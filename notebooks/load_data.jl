#####
# This notebook loads the microbial data for PCI-VEP-MBiome
#####

khula_visits = @chain CSV.read("ext/khula_ages_visit_reltable.csv", DataFrame; stringtype = String) begin
    select!(Not([:datacolor]))
    subset!(:visit => x -> x .== "3mo")
    subset!(:site => x -> x .== "ZAF")
end

samps_to_delete = [ "SEQ00447", "SEQ00449" ] ## Exclude due to outlier S.mitis high-leverage effect
# samps_to_delete = [ "SEQ00447" ] ## Exclude due to outlier S.mitis high-leverage effect

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
# select!(khula_pci_mdata, Not(:OldInfantVisAtt))
transform!(khula_pci_mdata, [:MaternalEntropy, :InfantVisAtt] => ((x, y) -> ( Leap.rangenormalize(x) .* Leap.rangenormalize(y) ) ) => :Product_Single)
### End of handling InfVisAttn data alteration

sort(setdiff(khula_visits.individualID, khula_pci_mdata.subject_id))
sort(setdiff(khula_pci_mdata.subject_id, khula_visits.individualID))

## 2025-01: Color by class

khula_pci_mdata = @chain CSV.read("/home/guilherme/Repos/Leap-notebooks/notebooks/2023-09-13-PCICollab/Catx2.csv", DataFrame; stringtype = String) begin
    transform!(:ID => (x -> "khula-" .* x) => :ID)
    transform!(:Cat => (x -> map(y -> (ismissing(y) ? "NA" : string(y)), x)) => :Cat)
    innerjoin(khula_pci_mdata, _ , on = :subject_id => :ID)
end

## 2025-02: Behavioural results
khula_pci_mdata = @chain CSV.read("GlitterWand.csv", DataFrame; stringtype = String) begin
    select!(Not([:NT_overall_conditional_entropy_U, :bVisual_Adjusted_CorrectedLength_U]))
    transform!(:subject_id => (x -> "khula-" .* x) => :subject_id)
    leftjoin(khula_pci_mdata, _ , on = :subject_id => :subject_id)
end

mdata_cols =  ["subject_id", "sample", "datasource", "pci_assess_age", "MaternalEntropy", "InfantVisAtt", "mbiome_sample_age", "child_sex", "MaternalEducation"]

focus_mdata = select(
    khula_pci_mdata,
    mdata_cols
)

# ## Run this only once every update to Humann!
# ## Remember this line depends on "scratch" location!
# Leap.write_gfs_arrow(; 
#     kind = "genefamilies",
#     file_glob_pattern = r"SEQ\d+_S\d+_genefamilies", 
#     file_root = "/home/guilherme/thunderbay-mount/processed/internal/humann/main",
#     reads_root = "/home/guilherme/thunderbay-mount/processed/internal/kneaddata/read_counts.csv",
#     suffix_pattern = r"_S\d+.*$",
#     output_path = scratchfiles("genefunctions"),
#     output_name = "genefamilies_HM40.arrow",
#     names = false,
#     stratified = false,
#     sample_filter = nothing
# )
# # ## Took 1h30 last time it ran on the entire folder

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

## Prevalence filters
unirefs_filtered_nonzeroentropy = let filt = get(unirefs, :MaternalEntropy) .>= 0.0
    # keepuni = vec(prevalence(unirefs[:, filt]) .> 0.0001) # Take note of prevalence filter used for the calculations!
    keepuni = vec(prevalence(unirefs[:, filt]) .> 0.025) # Take note of prevalence filter used for the calculations!
    # keepuni = vec(prevalence(unirefs[:, filt]) .> 0.05) # Take note of prevalence filter used for the calculations!
    # keepuni = vec(prevalence(unirefs[:, filt]) .> 0.10) # Take note of prevalence filter used for the calculations!
    unirefs[keepuni, filt]
end
relativeabundance!(unirefs_filtered_nonzeroentropy)

unimdata = DataFrame(get(unirefs))

#####
# Taxonomic Profiles
#####
khula_tax_data_species = @chain Leap.load_raw_metaphlan(
    "/home/guilherme/Repos/PCIVEPMBiome2024/ext/profiles";
    replace_pattern = r"SEQ0\d+_S\d+_mpa_v31_CHOCOPhlAn_201901_profile",
    cleanup_pattern = "_mpa_v31_CHOCOPhlAn_201901_profile.tsv") begin
    filter(t-> !ismissing(taxrank(t)), _[:, samplenames(_)])
    # filter(t-> taxrank(t) == :subspecies, _[:, samplenames(_)])
    filter(t-> taxrank(t) == :species, _[:, samplenames(_)])
    comm2wide()
    select(Not([:sample, :file]))
    # leftjoin(unimdata, _, on = :sample => :sample_base )
    filter_prevalence(0.001)
    # dropmissing()
end

khula_tax_data_genera = @chain Leap.load_raw_metaphlan(
    "/home/guilherme/Repos/PCIVEPMBiome2024/ext/profiles";
    replace_pattern = r"SEQ0\d+_S\d+_mpa_v31_CHOCOPhlAn_201901_profile",
    cleanup_pattern = "_mpa_v31_CHOCOPhlAn_201901_profile.tsv") begin
    filter(t-> !ismissing(taxrank(t)), _[:, samplenames(_)])
    filter(t-> taxrank(t) == :genus, _[:, samplenames(_)])
    comm2wide()
    select(Not([:sample, :file]))
    # leftjoin(unimdata, _, on = :sample => :sample_base )
    filter_prevalence(0.001)
    # dropmissing()
end
#####
# Alpha Diversity
#####

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

focus_mdata_sp = innerjoin(
    focus_mdata,
    DataFrame(
        :sample => khula_tax_data_species.sample_base,
        :alpha_shannon => map(x -> shannon(collect(x)), eachrow(khula_tax_data_species[:, 2:end]))
    ),
    on = :sample
)

focus_mdata_ge = innerjoin(
    focus_mdata,
    DataFrame(
        :sample => khula_tax_data_genera.sample_base,
        :alpha_shannon => map(x -> shannon(collect(x)), eachrow(khula_tax_data_genera[:, 2:end]))
    ),
    on = :sample
)

append!(mdata_cols, ["alpha_shannon"])

khula_pci_mbiome_data = innerjoin(
    focus_mdata_sp, 
  khula_tax_data_species, on = :sample => :sample_base)
CSV.write("manuscript/tables/mdata_table_in_bb3.csv", khula_pci_mbiome_data)

khula_pci_mbiome_data_genera = innerjoin(
    focus_mdata_ge, 
  khula_tax_data_genera, on = :sample => :sample_base)
CSV.write("manuscript/tables/mdata_table_genera_in_bb3.csv", khula_pci_mbiome_data_genera)

## Printing stats
println("Total number of participants: $(nrow(khula_pci_mbiome_data))")
println("Total number of males: $(sum(khula_pci_mbiome_data.child_sex .== 1)) ($(round(sum(khula_pci_mbiome_data.child_sex .== 1)/nrow(khula_pci_mbiome_data)*100; digits = 1))%)")##0 is female, 1 is male
println("Total number of females: $(sum(khula_pci_mbiome_data.child_sex .== 0)) ($(round(sum(khula_pci_mbiome_data.child_sex .== 0)/nrow(khula_pci_mbiome_data)*100; digits = 1))%)")##0 is female, 1 is male
println("Range of ages: $(ceil(minimum(khula_pci_mbiome_data.mbiome_sample_age)*30.5))-$(ceil(maximum(khula_pci_mbiome_data.mbiome_sample_age)*30.5))")##0 is female, 1 is male
