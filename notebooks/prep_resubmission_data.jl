#####
# This script generates Figure 1 for PCI-VEP-Mbiome Manuscript
#####
using CSV
using DataFrames
using Chain
using CategoricalArrays
using Microbiome
using PCIVEPMBiome2024

paper_sample_set = readlines(joinpath(Base.pwd(), "data", "initial_samples.txt"))
khula_pci_mbiome_data, unstratified_unirefs_filtered, stratified_unirefs_filtered = load_pcimbiome_data(paper_sample_set)

feeding_data = CSV.read(joinpath(Base.pwd(), "ext", "2026-02-20-KHULA_Feeding.csv"), DataFrame; stringtype = String)

delivery_data = @chain CSV.read(joinpath(Base.pwd(), "ext", "2026-04-10-KhulaSpace_Guo_SubjectMetadata.csv"), DataFrame; stringtype = String) begin
    transform!(:subject_id => ( x -> "khula-" .* x) => :subject_id; renamecols = false)
end

added_mdata_table = @chain khula_pci_mbiome_data begin
    innerjoin(feeding_data, on = :sample)
    innerjoin(delivery_data, on = :subject_id)
    select(
        [
            :subject_id,
            :sample,
            :datasource,
            :pci_assess_age,
            :mbiome_sample_age,
            :child_sex,
            :MaternalEducation,
            :feeding_state,
            :delivery_mode,
            :MaternalEntropy,
            :InfantVisAtt,
            :alpha_shannon,
        ], :
    )
    dropmissing()
end

println(setdiff(khula_pci_mbiome_data.subject_id, added_mdata_table.subject_id))

CSV.write(joinpath(Base.pwd(), "data", "study_metadata.csv"), added_mdata_table)