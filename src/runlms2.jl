using DataFrames, GLM, StatsModels, CairoMakie, ThreadsX

function make_default(feature, name, kind)
    return (; Name = name, feature, coef = NaN, std_err = NaN, stat = NaN,
              pvalue = NaN, lower_95 = NaN, upper_95 = NaN, kind)
end

function postformat(ct::DataFrame, feature, model_kind)
    ct.feature = [feature]
    rename!(ct, "Lower 95%"=>"lower_95", "Upper 95%"=>"upper_95",
             "Coef."=>"coef", "Std. Error"=>"std_err")
    # model_kind == :linear && rename!(ct, "Pr(>|t|)"=>"pvalue", "t"=>"stat")
    # model_kind == :logistic && rename!(ct, "Pr(>|z|)"=>"pvalue", "z"=>"stat")
    rename!(ct, "Pr(>|z|)"=>"pvalue", "z"=>"stat") ## for GLM on both
    select!(ct, Cols(:Name, :feature, :))
    return ct[1, :]
end

function runlms2(indf::DataFrame, feature_vector::Vector{String};
    modcols=["mbiome_sample_age", "child_sex", "InfantVisAtt", "MaternalEntropy"],
    prevalence_threshold=0.1, abundance_threshold=1.0,
    model_kind=:linear, bug_preproc=:none, include_sex=true,
    save_figures=true, fig_path="manuscript/figures/glm")

    f1=@formula(bug ~ mbiome_sample_age + child_sex + MaternalEntropy)
    f2=@formula(bug ~ mbiome_sample_age + child_sex + InfantVisAtt)
    f3=@formula(bug ~ mbiome_sample_age + child_sex + MaternalEntropy + InfantVisAtt + InfantVisAtt*MaternalEntropy)

    return DataFrame(ThreadsX.map(feature_vector) do feature
        ab = vec(indf[:, feature])
        over0 = ab .> 0
        if mean(over0) < prevalence_threshold || mean(ab[over0]) < abundance_threshold
            return (
                make_default(feature, "MaternalEntropy", "f1_$(model_kind)_unirefs"),
                make_default(feature, "InfantVisAtt", "f2_$(model_kind)_unirefs"),
                make_default(feature, "MaternalEntropy", "f3_$(model_kind)_unirefs"),
                make_default(feature, "InfantVisAtt", "f3_$(model_kind)_unirefs"),
                make_default(feature, "InfantVisAtt & MaternalEntropy", "f3_$(model_kind)_unirefs")
            )
        end

        df = indf[!, modcols]
        df.bug = model_kind == :logistic ? over0 : ab

        if model_kind == :linear
            df = subset(df, :bug => x -> x .> 0.0)
            if bug_preproc == :autoscale
                df.bug = (df.bug .- mean(df.bug)) ./ std(df.bug)
            elseif bug_preproc == :log2
                df.bug = log2.((df.bug ./ 100.0) .+ 1e-5)
            elseif bug_preproc == :log10
                df.bug = log10.((df.bug ./ 100.0) .+ 1e-5)
            elseif bug_preproc == :rclr
                robust_geomean = exp(mean(log.(df.bug[over0])))
                df.bug = log.(df.bug .+ 1e-5) .- log(robust_geomean)
            elseif bug_preproc == :arcsin
                df.bug = asin.(sqrt.(df.bug ./ 100))
            end

            lmodf1 = glm(f1, df, GLM.Gamma(), LogLink())
            lmodf2 = glm(f2, df, GLM.Gamma(), LogLink())
            lmodf3 = glm(f3, df, GLM.Gamma(), LogLink())

            # lmodf1, lmodf2, lmodf3 = lm.(Ref.((f1, f2, f3)), Ref(df); dropcollinear=true)

            if save_figures
                fig = Figure(); ax = Axis(fig[1, 1], xlabel="InfantVisAtt", ylabel="LOG2(abundance)",
                                          title="Bug vs InfantVisAtt ($feature)")
                scatter!(ax, df.InfantVisAtt, df.bug, color=(:blue, 0.6), markersize=8)
                save(joinpath(fig_path, "Linear_$(feature).png"), fig)
            end

        elseif model_kind == :logistic
            pseudo_ones = DataFrames.transform(deepcopy(df), :bug => x -> 1.0)
            pseudo_zeros = DataFrames.transform(deepcopy(df), :bug => x -> 0.0)
            df = vcat(df, pseudo_ones, pseudo_zeros)
            wts = vcat(fill(1.0, size(indf, 1)), fill(0.1, size(pseudo_ones, 1) + size(pseudo_zeros, 1)))

            lmodf1, lmodf2, lmodf3 = map(f -> glm(f, df, Binomial(), ProbitLink(); wts=wts, dropcollinear=true), (f1, f2, f3))

            if save_figures
                visvals = range(minimum(df.InfantVisAtt), stop=maximum(df.InfantVisAtt), length=200)
                df_pred = DataFrame(
                    mbiome_sample_age=fill(mean(df.mbiome_sample_age), 200),
                    child_sex=fill(mean(df.child_sex), 200),
                    InfantVisAtt=visvals
                )
                pred = GLM.predict(lmodf2, df_pred)
                fig = Figure(); ax = Axis(fig[1, 1], xlabel="InfantVisAtt", ylabel="Predicted probability",
                                          title="Bug vs InfantVisAtt")
                lines!(ax, visvals, pred, color=:red, linewidth=2)
                scatter!(ax, df.InfantVisAtt, df.bug, color=:blue, markersize=8)
                save(joinpath(fig_path, "Logistic_$(feature).png"), fig)
            end
        else
            throw(ArgumentError("Unsupported model_kind: $model_kind"))
        end

        ct = coeftable.(Ref.((lmodf1, lmodf2, lmodf3))) .|> DataFrame
        names = ["MaternalEntropy", "InfantVisAtt", "MaternalEntropy", "InfantVisAtt", "InfantVisAtt & MaternalEntropy"]
        idxs = include_sex ? [4, 4, 4, 5, 6] : [3, 3, 3, 4, 5]

        formatted = [postformat(ct[i][idxs[i]:idxs[i], :], feature, model_kind) for i in 1:5]
        formatted .= [merge(f, (; kind="f$(div(i-1, 2)+1)_$(model_kind)_unirefs")) for (i, f) in enumerate(formatted)]
        return tuple(formatted...)
    end)
end
