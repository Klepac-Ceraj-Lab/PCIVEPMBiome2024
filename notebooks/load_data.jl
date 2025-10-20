#####
# Common dataloading script
#####

using Chain, DataFrames, CSV, PCIVEPMBiome2024

# 1. Set of samples considered for the manuscript
paper_sample_set = readlines(joinpath(Base.pwd(), "data", "samples.txt"))

# 2. Core metadata

mdata_cols =  ["subject_id", "sample", "datasource", "pci_assess_age", "MaternalEntropy", "InfantVisAtt", "mbiome_sample_age", "child_sex", "MaternalEducation"]

study_mdata = @chain CSV.read(joinpath(Base.pwd(), "data", "study_metadata.csv"), DataFrame) begin
    select!(mdata_cols)
    subset!(:sample => x -> x .∈ Ref(paper_sample_set))
    unique!(:sample)
end

append!(mdata_cols, ["alpha_shannon"]) ## That will be added to the taxonomic profiles below...

# 3. Taxonomic profiles, species-level

khula_tax_data_species = @chain CSV.read(joinpath(Base.pwd(), "data", "species_data_MP3.csv"), DataFrame) begin
    subset!(:sample_base => x -> x .∈ Ref(paper_sample_set))
    insertcols!(_, 2, :alpha_shannon => map(x -> shannon(collect(x)), eachrow(_[:, 2:end])))
    unique!(:sample_base)
    filter_prevalence(0.0001)  ## Effectively removes all-zero taxa
end
khula_pci_mbiome_data = innerjoin(study_mdata, khula_tax_data_species, on = :sample => :sample_base)

# 4. Taxonomic profiles, genus-level

khula_tax_data_genus = @chain CSV.read(joinpath(Base.pwd(), "data", "genus_data_MP3.csv"), DataFrame) begin
    subset!(:sample_base => x -> x .∈ Ref(paper_sample_set))
    insertcols!(_, 2, :alpha_shannon => map(x -> shannon(collect(x)), eachrow(_[:, 2:end])))
    unique!(:sample_base)
    filter_prevalence(0.0001)  ## Effectively removes all-zero taxa
end
khula_pci_mbiome_data_genus = innerjoin(study_mdata, khula_tax_data_genus, on = :sample => :sample_base)

# 5. Functional UniRef90 profiles, unstratified.

unstratified_unirefs = load_custom_unirefs(UnirefProfiles(), study_mdata; custom_arrow_path = joinpath(Base.pwd(), "data", "genefamilies_HM36_unstratified.arrow"))
relativeabundance!(unstratified_unirefs)
## Prevalence filters
unstratified_unirefs_filtered = let filt = get(unstratified_unirefs, :MaternalEntropy) .>= 0.0
    keepuni = vec(prevalence(unstratified_unirefs[:, filt]) .> 0.025) # Take note of prevalence filter used for the calculations!
    unstratified_unirefs[keepuni, filt]
end

unimdata = DataFrame(get(unstratified_unirefs_filtered))

# 6. Functional UniRef90 profiles, stratified.

stratified_unirefs = load_custom_unirefs(UnirefProfiles(), study_mdata; custom_arrow_path = joinpath(Base.pwd(), "data", "genefamilies_HM36_stratified.arrow"))
relativeabundance!(stratified_unirefs)
## Prevalence filters
stratified_unirefs_filtered = let filt = get(stratified_unirefs, :MaternalEntropy) .>= 0.0
    keepuni = vec(prevalence(stratified_unirefs[:, filt]) .> 0.025) # Take note of prevalence filter used for the calculations!
    stratified_unirefs[keepuni, filt]
end

# 7. Printing Stats

println("Total number of participants: $(nrow(khula_pci_mbiome_data))")
println("Total number of males: $(sum(khula_pci_mbiome_data.child_sex .== 1)) ($(round(sum(khula_pci_mbiome_data.child_sex .== 1)/nrow(khula_pci_mbiome_data)*100; digits = 1))%)")##0 is female, 1 is male
println("Total number of females: $(sum(khula_pci_mbiome_data.child_sex .== 0)) ($(round(sum(khula_pci_mbiome_data.child_sex .== 0)/nrow(khula_pci_mbiome_data)*100; digits = 1))%)")##0 is female, 1 is male
println("Range of ages: $(ceil(minimum(khula_pci_mbiome_data.mbiome_sample_age)*30.5))-$(ceil(maximum(khula_pci_mbiome_data.mbiome_sample_age)*30.5))")##0 is female, 1 is male

let
    taxprofile_richness = [ sum(el .> 0.0) for el in eachrow(Matrix(khula_tax_data_species[:, 3:end])) ]
    println("Species richness: $(floor(Int64, mean(taxprofile_richness))) (SD = $(floor(Int64, Statistics.std(taxprofile_richness))))")

    geneprofile_richness = [ sum(el .> 0.0) for el in eachcol(abundances(stratified_unirefs)) ]
    println("Gene richness: $(floor(Int64, mean(geneprofile_richness))) (SD = $(floor(Int64, Statistics.std(geneprofile_richness))))")

    geneprofile_richness_filtered = [ sum(el .> 0.0) for el in eachcol(abundances(stratified_unirefs_filtered)) ]
    println("Gene richness: $(floor(Int64, mean(geneprofile_richness_filtered))) (SD = $(floor(Int64, Statistics.std(geneprofile_richness_filtered))))")

end