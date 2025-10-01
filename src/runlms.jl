function runlms(
    indf::DataFrame, feature_vector::Vector{String},# outfile::String,
    modcols=[ "mbiome_sample_age", "child_sex", "InfantVisAtt", "MaternalEntropy" ];
    prevalence_threshold = 0.1,
    abundance_threshold = 1.0,
    model_kind = :linear, # OR :logistic
    bug_preproc = :none, # [:none, :autoscale, :log2, :arcsin]
    bugrank = false
    )

    f1=@formula(bug ~ mbiome_sample_age + child_sex + MaternalEntropy)
    f2=@formula(bug ~ mbiome_sample_age + child_sex + InfantVisAtt)
    f3=@formula(bug ~ mbiome_sample_age + child_sex + MaternalEntropy + InfantVisAtt + InfantVisAtt*MaternalEntropy)

    # lmresults = DataFrame(ThreadsX.map(feature_vector) do this_feature
    lmresults = DataFrame(map(feature_vector) do this_feature

    @debug this_feature

        default_nt_f1 = (; Name = "MaternalEntropy", feature = this_feature, coef = NaN, std_err = NaN, stat = NaN, pvalue = NaN, lower_95=NaN, upper_95=NaN, kind="f1_"*string(model_kind)*"_unirefs")
        default_nt_f2 = (; Name = "InfantVisAtt", feature = this_feature, coef = NaN, std_err = NaN, stat = NaN, pvalue = NaN, lower_95=NaN, upper_95=NaN, kind="f2_"*string(model_kind)*"_unirefs")
        default_nt_f3_MatEnt = (; Name = "MaternalEntropy", feature = this_feature, coef = NaN, std_err = NaN, stat = NaN, pvalue = NaN, lower_95=NaN, upper_95=NaN, kind="f3_"*string(model_kind)*"_unirefs")
        default_nt_f3_InfVis = (; Name = "InfantVisAtt", feature = this_feature, coef = NaN, std_err = NaN, stat = NaN, pvalue = NaN, lower_95=NaN, upper_95=NaN, kind="f3_"*string(model_kind)*"_unirefs")
        default_nt_f3_Interact = (; Name = "InfantVisAtt & MaternalEntropy", feature = this_feature, coef = NaN, std_err = NaN, stat = NaN, pvalue = NaN, lower_95=NaN, upper_95=NaN, kind="f3_"*string(model_kind)*"_unirefs")

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

        df = indf[!, modcols]

        if model_kind == :linear

            println("Linear regression for $this_feature")

            df.bug = ab
            df = subset(df, :bug => (x -> (x .> 0.0)))

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
                joinpath("manuscript", "figures", "glm", "Linear_$(this_feature).png"),
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
                joinpath("manuscript", "figures", "glm", "Logistic_$(this_feature).png"),
                fig
            )         

        else

            throw(ArgumentError("No model_kind: $model_kind"))    

        end

        # # ————————————————
        # # helper: median p-value by subsampling
        # function bootstrap_median_p(formula, data; row_idx::Int, B=100, frac=0.8)
        #     n = nrow(data)
        #     ps = Float64[]
        #     for _ in 1:B
        #         idx = sample(1:n, round(Int,frac*n), replace=false)
        #         sub = data[idx, :]
        #         m = glm(formula, sub, GLM.Gamma(), LogLink(); maxiter=1000)
        #         ct = DataFrame(coeftable(m))
        #         push!(ps, ct[row_idx, :"Pr(>|z|)"])
        #     end
        #     return median(ps)
        # end
        # # ————————————————
        # # bootstrap and overwrite p-values in each fit
        # # f1’s coeff “MaternalEntropy” is row 4
        # p1_med = bootstrap_median_p(f1, df; row_idx=4, B=200, frac=0.5)
        # df1 = DataFrame(coeftable(lmodf1))
        # df1[4, :"Pr(>|z|)"] = p1_med

        # # f2’s coeff “InfantVisAtt” is row 4
        # p2_med = bootstrap_median_p(f2, df; row_idx=4, B=200, frac=0.5)
        # df2 = DataFrame(coeftable(lmodf2))
        # df2[4, :"Pr(>|z|)"] = p2_med

        # # f3: MaternalEntropy row 4, InfantVisAtt row 5, Interaction row 6
        # p3_m = bootstrap_median_p(f3, df; row_idx=4, B=200, frac=0.5)
        # p3_i = bootstrap_median_p(f3, df; row_idx=5, B=200, frac=0.5)
        # p3_x = bootstrap_median_p(f3, df; row_idx=6, B=200, frac=0.5)
        # df3 = DataFrame(coeftable(lmodf3))
        # df3[4, :"Pr(>|z|)"] = p3_m
        # df3[5, :"Pr(>|z|)"] = p3_i
        # df3[6, :"Pr(>|z|)"] = p3_x
        # # ————————————————

        # # the rest of your postformat logic, but now using df1, df2, df3
        # ct1 = df1[4:4, :]
        # ct2 = df2[4:4, :]
        # ct3 = df3[4:4, :]
        # ct4 = df3[5:5, :]
        # ct5 = df3[6:6, :]

        ct1 = DataFrame(coeftable(lmodf1))[4:4, :]
        @assert ct1.Name[1] == "MaternalEntropy"

        ct2 = DataFrame(coeftable(lmodf2))[4:4, :]
        @assert ct2.Name[1] == "InfantVisAtt"

        ct3 = DataFrame(coeftable(lmodf3))[4:4, :]
        @assert ct3.Name[1] == "MaternalEntropy"

        ct4 = DataFrame(coeftable(lmodf3))[5:5, :]
        @assert ct4.Name[1] == "InfantVisAtt"

        ct5 = DataFrame(coeftable(lmodf3))[6:6, :]
        @assert ct5.Name[1] == "InfantVisAtt & MaternalEntropy"

        function postformat(ct, featnmm)
            ct.feature = [ featnmm ]
            rename!(ct, "Lower 95%"=> "lower_95", "Upper 95%"=> "upper_95", "Coef."=> "coef", "Std. Error"=>"std_err")
            model_kind == :linear && rename!(ct, "Pr(>|t|)"=>"pvalue", "t"=>"stat")
            model_kind == :logistic && rename!(ct, "Pr(>|z|)"=>"pvalue", "z"=>"stat")
            # rename!(ct, "Pr(>|z|)"=>"pvalue", "z"=>"stat")
            select!(ct, Cols(:Name, :feature, :))
            return(ct[1, :])
        end

        ct1_fmt = postformat(ct1, this_feature)
        ct2_fmt = postformat(ct2, this_feature)
        ct3_fmt = postformat(ct3, this_feature)
        ct4_fmt = postformat(ct4, this_feature)
        ct5_fmt = postformat(ct5, this_feature)

        return (
            (; ct1_fmt... , kind="f1_"*string(model_kind)*"_unirefs"),
            (; ct2_fmt... , kind="f2_"*string(model_kind)*"_unirefs"),
            (; ct3_fmt... , kind="f3_"*string(model_kind)*"_unirefs"),
            (; ct4_fmt... , kind="f3_"*string(model_kind)*"_unirefs"),
            (; ct5_fmt... , kind="f3_"*string(model_kind)*"_unirefs")
        )

        end)

        # subset!(lmresults, "pvalue"=> ByRow(!isnan))
        # DataFrames.transform!(lmresults, :pvalue => (col-> MultipleTesting.adjust(collect(col), BenjaminiHochberg())) => :qvalue)
        # sort!(lmresults, :qvalue)

        # CSV.write(outfile, lmresults)
    lmresults
end