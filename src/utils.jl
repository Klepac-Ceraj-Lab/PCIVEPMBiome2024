function load_custom_unirefs(::UnirefProfiles, timepoint_metadata::DataFrame; custom_arrow_path = "data/genefamilies_HM36_stratified.arrow")
    comm = Leap.read_gfs_arrow(; arrow_path = custom_arrow_path)
    insert!(comm, timepoint_metadata; namecol=:sample)

    valid_samples = intersect(name.(comm.samples), timepoint_metadata.sample)
    comm = comm[:, valid_samples]

    if isempty(comm.samples)
        @warn "No overlapping samples between profile and metadata!"
    end

    return comm
end