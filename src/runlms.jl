function runlms(
    indf::DataFrame, feature_vector::Vector{String},# outfile::String,
    modcols=[ "mbiome_sample_age", "child_sex", "maternal_education", "feeding_state", "delivery_mode", "InfantVisAtt", "MaternalEntropy" ];
    prevalence_threshold = 0.1,
    abundance_threshold = 1.0,
    model_kind = :linear, # OR :logistic
    bug_preproc = :none, # [:none, :autoscale, :log2, :arcsin]
    bugrank = false,
    sample_id_col::Symbol = :sample,
    outdir::String = joinpath("manuscript", "figures", "glm"),
    export_cooks::Bool = true
    )

    isdir(outdir) && throw(error("Output directory already exists!"))
    mkpath(outdir)

    ## Accessory function for cooks D
    function cooks_distance_lm(this_model)
        X = modelmatrix(this_model)
        r = residuals(this_model)

        n, p = size(X)
        σ2 = sum(r .^ 2) / (n - p)

        H = X * inv(X' * X) * X'
        h = LinearAlgebra.diag(H)

        D = (r .^ 2 ./ (p * σ2)) .* (h ./ (1 .- h).^2)

        return D
    end

    f1=@formula(bug ~ mbiome_sample_age + child_sex + maternal_education + feeding_state + delivery_mode + MaternalEntropy)
    f2=@formula(bug ~ mbiome_sample_age + child_sex + maternal_education + feeding_state + delivery_mode + InfantVisAtt)
    f3=@formula(bug ~ mbiome_sample_age + child_sex + maternal_education + feeding_state + delivery_mode + MaternalEntropy + InfantVisAtt + InfantVisAtt*MaternalEntropy)

    # lmresults = DataFrame(ThreadsX.map(feature_vector) do this_feature
    lmresults = DataFrame(map(feature_vector) do this_feature

    @debug this_feature

        default_nt_f1 = (; Name = "MaternalEntropy", feature = this_feature, coef = NaN, std_err = NaN, stat = NaN, pvalue = NaN, lower_95=NaN, upper_95=NaN, rsquare=NaN, kind="f1_"*string(model_kind)*"_taxa")
        default_nt_f2 = (; Name = "InfantVisAtt", feature = this_feature, coef = NaN, std_err = NaN, stat = NaN, pvalue = NaN, lower_95=NaN, upper_95=NaN, rsquare=NaN, kind="f2_"*string(model_kind)*"_taxa")
        default_nt_f3_MatEnt = (; Name = "MaternalEntropy", feature = this_feature, coef = NaN, std_err = NaN, stat = NaN, pvalue = NaN, lower_95=NaN, upper_95=NaN, rsquare=NaN, kind="f3_"*string(model_kind)*"_taxa")
        default_nt_f3_InfVis = (; Name = "InfantVisAtt", feature = this_feature, coef = NaN, std_err = NaN, stat = NaN, pvalue = NaN, lower_95=NaN, upper_95=NaN, rsquare=NaN, kind="f3_"*string(model_kind)*"_taxa")
        default_nt_f3_Interact = (; Name = "InfantVisAtt & MaternalEntropy", feature = this_feature, coef = NaN, std_err = NaN, stat = NaN, pvalue = NaN, lower_95=NaN, rsquare=NaN, upper_95=NaN, kind="f3_"*string(model_kind)*"_taxa")

        ab = vec(indf[:, this_feature])
        over0 = ab .> 0
        if (sum(over0) / size(indf, 1) < prevalence_threshold) 
            return(
                default_nt_f1, default_nt_f2, default_nt_f3_MatEnt, default_nt_f3_InfVis, default_nt_f3_Interact
            )
        elseif ( mean(ab[over0]) < abundance_threshold )
            return(
                default_nt_f1, default_nt_f2, default_nt_f3_MatEnt, default_nt_f3_InfVis, default_nt_f3_Interact
            )       
        end

        df = select(indf, vcat([sample_id_col], Symbol.(modcols)))

        if model_kind == :linear

            println("Linear regression for $this_feature")

            df.bug = ab
            df = subset(df, :bug => (x -> (x .> 0.0)))

            if ( ( !all(map( x -> length(unique(x)), eachcol(select(df, modcols))) .> 1) ) | ( nrow(df) <= length(modcols) ) )
                return(
                    default_nt_f1, default_nt_f2, default_nt_f3_MatEnt, default_nt_f3_InfVis, default_nt_f3_Interact
                )
            end

            if bugrank == true

                df.bug = invperm(sortperm(df.bug))

            else

                if bug_preproc == :autoscale
                    df.bug = Leap.autonormalize(df.bug)
                elseif bug_preproc == :log2
                    df.bug = log2.(df.bug .+ 1e-5)
                elseif bug_preproc == :log10
                    df.bug = log10.((df.bug ./ 100.0) .+ 1e-5)
                elseif bug_preproc == :rclr
                    robust_geomean = exp(mean(log.(df.bug[over0])))
                    df.bug = map(x -> (log(x .+ 1e-5) .- log(robust_geomean)), df.bug)
                elseif bug_preproc == :arcsin
                    # df.bug = asin.(df.bug ./ 100)
                    df.bug = asin.(sqrt.(df.bug))
                    # df.bug = asinh.(df.bug)
                end
            
            end

            # lmodf1 = glm(f1, df, GLM.Normal())
            # lmodf2 = glm(f2, df, GLM.Normal())
            # lmodf3 = glm(f3, df, GLM.Normal())
            # lmodf1 = glm(f1, df, GLM.Gamma(); maxiter = 1000)
            # lmodf2 = glm(f2, df, GLM.Gamma(); maxiter = 1000)
            # lmodf3 = glm(f3, df, GLM.Gamma(); maxiter = 1000)
            # lmodf1 = glm(f1, df, GLM.Normal(), LogLink())
            # lmodf2 = glm(f2, df, GLM.Normal(), LogLink())
            # lmodf3 = glm(f3, df, GLM.Normal(), LogLink())
            # lmodf1 = glm(f1, df, GLM.Gamma(), LogLink(); maxiter = 1000)
            # lmodf2 = glm(f2, df, GLM.Gamma(), LogLink(); maxiter = 1000)
            # lmodf3 = glm(f3, df, GLM.Gamma(), LogLink(); maxiter = 1000)
            lmodf1 = lm(f1, df; dropcollinear=true)
            lmodf2 = lm(f2, df; dropcollinear=true)
            lmodf3 = lm(f3, df; dropcollinear=true)

            if export_cooks
            
                try
                    cooks_df = DataFrame(
                        sample_id = df[!, sample_id_col],
                        feature = this_feature,
                        cooksD_f1 = cooks_distance_lm(lmodf1),
                        cooksD_f2 = cooks_distance_lm(lmodf2),
                        cooksD_f3 = cooks_distance_lm(lmodf3)
                    )

                    CSV.write(
                        joinpath(outdir, "CooksD_$(this_feature).csv"),
                        cooks_df
                    )
                catch e
                end
            
            end

            ## Figure generation step
            fig = Figure()
            ax = Axis(
                fig[1, 1],
                xlabel="InfantVisAtt",
                ylabel="LOG2(abundance)", 
                title="Bug vs InfantVisAtt ($(this_feature))"
            )
            # Plot the data points
            scatter!(ax, df.InfantVisAtt, df.bug, color = (:blue, 0.6), markersize=8)
            save(
                joinpath(outdir, "Linear_$(this_feature).png"),
                fig
            )    

        elseif model_kind == :logistic

            println("Logistic regression for $this_feature")

            df.bug = over0
            pseudo_ones = @chain deepcopy(df) begin
                transform!(:bug => (x -> 1.0) => :bug)
            end
            pseudo_zeros = @chain deepcopy(df) begin
                transform!(:bug => (x -> 0.0) => :bug)
            end

            combined_wts = vcat(
                repeat([ 1 ], nrow(df)),
                repeat([ 0.1 ], nrow(pseudo_ones)),
                repeat([ 0.1 ], nrow(pseudo_zeros))
            )

            df = vcat(df, pseudo_ones, pseudo_zeros)

            lmodf1 = glm(f1, df, Binomial(), ProbitLink(); wts = combined_wts, dropcollinear=true)
            lmodf2 = glm(f2, df, Binomial(), ProbitLink(); wts = combined_wts, dropcollinear=true)
            lmodf3 = glm(f3, df, Binomial(), ProbitLink(); wts = combined_wts, dropcollinear=true)

            ## PLOTTING BLOCK
            # Create a range for InfantVisAtt
            InfantVisAtt_vals = range(minimum(df.InfantVisAtt), stop=maximum(df.InfantVisAtt), length=200)
            # Create a DataFrame for predictions
            df_pred = DataFrame(
                :mbiome_sample_age => fill(mean(df.mbiome_sample_age), length(InfantVisAtt_vals)),
                :child_sex => fill(mean(df.child_sex), length(InfantVisAtt_vals)),
                :InfantVisAtt => InfantVisAtt_vals
            )
            # Get predicted probabilities from the model
            pred_probs = GLM.predict(lmodf2, df_pred)
            ## Figure generation step
            fig = Figure()
            ax = Axis(fig[1, 1], xlabel="InfantVisAtt", ylabel="Predicted probability", 
                      title="Bug vs InfantVisAtt (Other predictors fixed)")           
            # Plot the predicted logistic curve
            lines!(ax, InfantVisAtt_vals, pred_probs, color=:red, linewidth=2)
            # Optionally, overlay the original data points
            scatter!(ax, df.InfantVisAtt, df.bug, color=:blue, markersize=8)
            save(
                joinpath(outdir, "Logistic_$(this_feature).png"),
                fig
            )         

        else
            throw(ArgumentError("No model_kind: $model_kind"))
        end

        ct1 = DataFrame(coeftable(lmodf1))[7:7, :]
        @assert ct1.Name[1] == "MaternalEntropy"

        ct2 = DataFrame(coeftable(lmodf2))[7:7, :]
        @assert ct2.Name[1] == "InfantVisAtt"

        ct3 = DataFrame(coeftable(lmodf3))[7:7, :]
        @assert ct3.Name[1] == "MaternalEntropy"

        ct4 = DataFrame(coeftable(lmodf3))[8:8, :]
        @assert ct4.Name[1] == "InfantVisAtt"

        ct5 = DataFrame(coeftable(lmodf3))[9:9, :]
        @assert ct5.Name[1] == "InfantVisAtt & MaternalEntropy"

        function postformat(ct, featnmm, this_model)
            ct.feature = [ featnmm ]
            rename!(ct, "Lower 95%"=> "lower_95", "Upper 95%"=> "upper_95", "Coef."=> "coef", "Std. Error"=>"std_err")
            model_kind == :linear && rename!(ct, "Pr(>|t|)"=>"pvalue", "t"=>"stat")
            model_kind == :logistic && rename!(ct, "Pr(>|z|)"=>"pvalue", "z"=>"stat")
            # rename!(ct, "Pr(>|z|)"=>"pvalue", "z"=>"stat")
            ct.rsquare = ( (model_kind == :linear) ? [r2(this_model)] : [NaN] )
            select!(ct, Cols(:Name, :feature, :))
            return(ct[1, :])
        end

        ct1_fmt = postformat(ct1, this_feature, lmodf1)
        ct2_fmt = postformat(ct2, this_feature, lmodf2)
        ct3_fmt = postformat(ct3, this_feature, lmodf3)
        ct4_fmt = postformat(ct4, this_feature, lmodf3)
        ct5_fmt = postformat(ct5, this_feature, lmodf3)

        return (
            (; ct1_fmt... , kind="f1_"*string(model_kind)*"_taxa"),
            (; ct2_fmt... , kind="f2_"*string(model_kind)*"_taxa"),
            (; ct3_fmt... , kind="f3_"*string(model_kind)*"_taxa"),
            (; ct4_fmt... , kind="f3_"*string(model_kind)*"_taxa"),
            (; ct5_fmt... , kind="f3_"*string(model_kind)*"_taxa")
        )

        end)

        # subset!(lmresults, "pvalue"=> ByRow(!isnan))
        # DataFrames.transform!(lmresults, :pvalue => (col-> MultipleTesting.adjust(collect(col), BenjaminiHochberg())) => :qvalue)
        # sort!(lmresults, :qvalue)

        # CSV.write(outfile, lmresults)
    lmresults
end