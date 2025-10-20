khula_rf_data = deepcopy(khula_pci_mbiome_data)
khula_rf_data = filter_prevalence(khula_rf_data, 0.05)
khula_rf_data.InfantVisAtt = Leap.rangenormalize(khula_rf_data.InfantVisAtt)

select!(khula_rf_data, mdata_cols, :)
select!(khula_rf_data, :subject_id, :sample, :datasource, :pci_assess_age, :InfantVisAtt, :)

if isfile(joinpath(Base.pwd(), "manuscript", "models", "regression_VOB_NoEnt_FullCV_Results.jld"))
    JLD2.@load joinpath(Base.pwd(), "manuscript", "models", "regression_VOB_NoEnt_FullCV_Results.jld") regression_VOB_NoEnt_FullCV
else
    regression_VOB_NoEnt_FullCV = Leap.probe_regression_randomforest(
        "regression_VOB_NoEnt_FullCV",
        khula_rf_data,
        identity,
        collect(findfirst(names(khula_rf_data) .== "mbiome_sample_age"):ncol(khula_rf_data)),
        :InfantVisAtt;
        split_strat = "subject",
        custom_input_group = nothing,
        unique_col = :sample,
        n_folds = 3,
        n_replicas = 5,
        n_rngs = 3,
        tuning_space = (; #PRODUCTION
            maxnodes_range = [ -1, 2, 4, 6, 8, 10 ],
            nodesize_range = [ 3, 5, 7, 9, 11 ],
            min_samples_split = [ 2, 3, 4 ],
            sampsize_range = [ 0.7, 0.8, 0.9 ],
            mtry_range = [ -1, 0, 10, 15, 20, 25, 30 ],
            ntrees_range = [ 64, 128, 256 ]
        )
    )
    sort(report_regression_merits(regression_VOB_NoEnt_FullCV), :Val_RMSE_mean)
    JLD2.@save joinpath(Base.pwd(), "manuscript", "models", "regression_VOB_NoEnt_FullCV_Results.jld") regression_VOB_NoEnt_FullCV
end

if isfile(joinpath(Base.pwd(), "manuscript", "models", "regression_VOB_WtEnt_FullCV_Results.jld"))
    JLD2.@load joinpath(Base.pwd(), "manuscript", "models", "regression_VOB_WtEnt_FullCV_Results.jld") regression_VOB_WtEnt_FullCV
else
    regression_VOB_WtEnt_FullCV = Leap.probe_regression_randomforest(
        "regression_VOB_WtEnt_FullCV",
        khula_rf_data,
        identity,
        collect(findfirst(names(khula_rf_data) .== "MaternalEntropy"):ncol(khula_rf_data)),
        :InfantVisAtt;
        split_strat = "subject",
        custom_input_group = nothing,
        unique_col = :sample,
        n_folds = 3,
        n_replicas = 5,
        n_rngs = 3,
        tuning_space = (; #PRODUCTION
            maxnodes_range = [ -1, 2, 4, 6, 8, 10 ],
            nodesize_range = [ 3, 5, 7, 9, 11 ],
            min_samples_split = [ 2, 3, 4 ],
            sampsize_range = [ 0.7, 0.8, 0.9 ],
            mtry_range = [ -1, 0, 10, 15, 20, 25, 30 ],
            ntrees_range = [ 64, 128, 256 ]
        )
    )
    sort(report_regression_merits(regression_VOB_WtEnt_FullCV), :Val_RMSE_mean)
    JLD2.@save joinpath(Base.pwd(), "manuscript", "models", "regression_VOB_WtEnt_FullCV_Results.jld") regression_VOB_WtEnt_FullCV
end

#####
# Now with MaternalEntropy
#####
select!(khula_rf_data, :subject_id, :sample, :datasource, :pci_assess_age, :MaternalEntropy, :)

if isfile(joinpath(Base.pwd(), "manuscript", "models", "regression_ENT_NoVobFullCV_Results.jld"))
    JLD2.@load joinpath(Base.pwd(), "manuscript", "models", "regression_ENT_NoVobFullCV_Results.jld") regression_ENT_NoVobFullCV
else
    regression_ENT_NoVobFullCV = Leap.probe_regression_randomforest(
        "regression_ENT_NoVobFullCV",
        khula_rf_data,
        identity,
        collect(findfirst(names(khula_rf_data) .== "mbiome_sample_age"):ncol(khula_rf_data)),
        :MaternalEntropy;
        split_strat = "subject",
        custom_input_group = nothing,
        unique_col = :sample,
        n_folds = 3,
        n_replicas = 5,
        n_rngs = 3,
        tuning_space = (; #PRODUCTION
            maxnodes_range = [ -1, 2, 4, 6, 8, 10 ],
            nodesize_range = [ 3, 5, 7, 9, 11 ],
            min_samples_split = [ 2, 3, 4 ],
            sampsize_range = [ 0.7, 0.8, 0.9 ],
            mtry_range = [ -1, 0, 10, 15, 20, 25, 30 ],
            ntrees_range = [ 64, 128, 256 ]
        )
    )
    sort(report_regression_merits(regression_ENT_NoVobFullCV), :Val_RMSE_mean)
    JLD2.@save joinpath(Base.pwd(), "manuscript", "models", "regression_ENT_NoVobFullCV_Results.jld") regression_ENT_NoVobFullCV
end

if isfile(joinpath(Base.pwd(), "manuscript", "models", "regression_ENT_WtVobFullCV_Results.jld"))
    JLD2.@load joinpath(Base.pwd(), "manuscript", "models", "regression_ENT_WtVobFullCV_Results.jld") regression_ENT_WtVobFullCV
else
    regression_ENT_WtVobFullCV = Leap.probe_regression_randomforest(
        "regression_ENT_WtVobFullCV",
        khula_rf_data,
        identity,
        collect(findfirst(names(khula_rf_data) .== "InfantVisAtt"):ncol(khula_rf_data)),
        :MaternalEntropy;
        split_strat = "subject",
        custom_input_group = nothing,
        unique_col = :sample,
        n_folds = 3,
        n_replicas = 5,
        n_rngs = 3,
        tuning_space = (; #PRODUCTION
            maxnodes_range = [ -1, 2, 4, 6, 8, 10 ],
            nodesize_range = [ 3, 5, 7, 9, 11 ],
            min_samples_split = [ 2, 3, 4 ],
            sampsize_range = [ 0.7, 0.8, 0.9 ],
            mtry_range = [ -1, 0, 10, 15, 20, 25, 30 ],
            ntrees_range = [ 64, 128, 256 ]
        )
    )
    sort(report_regression_merits(regression_ENT_WtVobFullCV), :Val_RMSE_mean)
    JLD2.@save joinpath(Base.pwd(), "manuscript", "models", "regression_ENT_WtVobFullCV_Results.jld") regression_ENT_WtVobFullCV
end