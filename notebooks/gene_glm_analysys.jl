#####
# Running LMs
#####

let # MODEL 1
    indf = DataFrame(get(unirefs_filtered_nonzeroentropy))
    indf.InfantVisAtt = Leap.rangenormalize(indf.InfantVisAtt)
    indf.MaternalEntropy = Leap.rangenormalize(indf.MaternalEntropy)
    indf.mbiome_sample_age = Leap.rangenormalize(indf.mbiome_sample_age)
    outfiles = (
        "manuscript/FSEA/lms_Model1_mbiome_age.csv",
        "manuscript/FSEA/lms_Model1_child_sex.csv"
        )
    lmresults = DataFrame(ThreadsX.map(features(unirefs_filtered_nonzeroentropy)) do feat

    try
        ab = vec(abundances(unirefs_filtered_nonzeroentropy[feat, :]))
        df = indf[!, :]
        df.bug = Int8.(ab .> 0)
    
        ## MODEL 1
        mod = GLM.glm(
            @formula(
                bug ~
                mbiome_sample_age +
                child_sex
            ), df, Binomial(), ProbitLink(); dropcollinear=true
        )

        ct1 = DataFrames.Tables.rowtable(coeftable(mod))[2]
        @assert ct1.Name == "mbiome_sample_age"
            
        ct2 = DataFrames.Tables.rowtable(coeftable(mod))[3]
        @assert ct2.Name == "child_sex"
            
        return (
            (; ct1..., feature = string(feat), Cor = cor(df.bug, df.mbiome_sample_age), kind="unirefs"),
            (; ct2..., feature = string(feat), Cor = cor(df.bug, df.child_sex), kind="unirefs")
        )

    catch
        ct1 = (;Name = "mbiome_sample_age", var"Coef." = NaN, var"Std. Error" = NaN, z = NaN, var"Pr(>|z|)" = NaN, var"Lower 95%" = NaN, var"Upper 95%" = NaN)
        ct2 = (;Name = "child_sex", var"Coef." = NaN, var"Std. Error" = NaN, z = NaN, var"Pr(>|z|)" = NaN, var"Lower 95%" = NaN, var"Upper 95%" = NaN)
        
            return (
                (; ct1..., feature = string(feat), Cor = NaN, kind="unirefs"),
                (; ct2..., feature = string(feat), Cor = NaN, kind="unirefs")
            )

        end

    end)

    for (idx, outfile) in enumerate(outfiles)

        this_df = copy(DataFrame(lmresults[!, idx]))
        select!(this_df, Cols(:feature, :Name, :))
        rename!(this_df, "Pr(>|z|)"=>"pvalue");
        sort!(this_df, :pvalue); 

        CSV.write(outfile, this_df)

    end

end

let # MODEL 2
    indf = DataFrame(get(unirefs_filtered_nonzeroentropy))
    indf.InfantVisAtt = Leap.rangenormalize(indf.InfantVisAtt)
    indf.MaternalEntropy = Leap.rangenormalize(indf.MaternalEntropy)
    indf.mbiome_sample_age = Leap.rangenormalize(indf.mbiome_sample_age)
    outfiles = (
        "manuscript/FSEA/lms_Model2_child_sex.csv",
        "manuscript/FSEA/lms_Model2_visual.csv",
        "manuscript/FSEA/lms_Model2_entropy.csv"
        )
    lmresults = DataFrame(ThreadsX.map(features(unirefs_filtered_nonzeroentropy)) do feat

    try
        ab = vec(abundances(unirefs_filtered_nonzeroentropy[feat, :]))
        df = indf[!, :]
        df.bug = Int8.(ab .> 0)
    
        ## MODEL 2
        mod = GLM.glm(
            @formula(
                bug ~
                mbiome_sample_age +
                child_sex +
                InfantVisAtt +
                MaternalEntropy
            ), df, Binomial(), ProbitLink(); dropcollinear=true
        )

        ct2 = DataFrames.Tables.rowtable(coeftable(mod))[3]
        @assert ct2.Name == "child_sex"

        ct3 = DataFrames.Tables.rowtable(coeftable(mod))[4]
        @assert ct3.Name == "InfantVisAtt"

        ct4 = DataFrames.Tables.rowtable(coeftable(mod))[5]
        @assert ct4.Name == "MaternalEntropy"
            
        return (
            (; ct2..., feature = string(feat), Cor = cor(df.bug, df.child_sex), kind="unirefs"),
            (; ct3..., feature = string(feat), Cor = cor(df.bug, df.InfantVisAtt), kind="unirefs"),
            (; ct4..., feature = string(feat), Cor = cor(df.bug, df.MaternalEntropy), kind="unirefs")
        )

    catch
        ct2 = (;Name = "child_sex", var"Coef." = NaN, var"Std. Error" = NaN, z = NaN, var"Pr(>|z|)" = NaN, var"Lower 95%" = NaN, var"Upper 95%" = NaN)
        ct3 = (;Name = "InfantVisAtt", var"Coef." = NaN, var"Std. Error" = NaN, z = NaN, var"Pr(>|z|)" = NaN, var"Lower 95%" = NaN, var"Upper 95%" = NaN)
        ct4 = (;Name = "MaternalEntropy", var"Coef." = NaN, var"Std. Error" = NaN, z = NaN, var"Pr(>|z|)" = NaN, var"Lower 95%" = NaN, var"Upper 95%" = NaN)

            return (
                (; ct2..., feature = string(feat), Cor = NaN, kind="unirefs"),
                (; ct3..., feature = string(feat), Cor = NaN, kind="unirefs"),
                (; ct4..., feature = string(feat), Cor = NaN, kind="unirefs")
            )

        end

    end)

    for (idx, outfile) in enumerate(outfiles)

        this_df = copy(DataFrame(lmresults[!, idx]))
        select!(this_df, Cols(:feature, :Name, :))
        rename!(this_df, "Pr(>|z|)"=>"pvalue");
        sort!(this_df, :pvalue); 
    
        CSV.write(outfile, this_df)

    end

end

let # MODEL 3
    indf = DataFrame(get(unirefs_filtered_nonzeroentropy))
    indf.InfantVisAtt = Leap.rangenormalize(indf.InfantVisAtt)
    indf.MaternalEntropy = Leap.rangenormalize(indf.MaternalEntropy)
    indf.mbiome_sample_age = Leap.rangenormalize(indf.mbiome_sample_age)
    outfiles = (
        "manuscript/FSEA/lms_Model3_child_sex.csv",
        "manuscript/FSEA/lms_Model3_visual.csv",
        "manuscript/FSEA/lms_Model3_entropy.csv",
        "manuscript/FSEA/lms_Model3_product.csv",
        )
    lmresults = DataFrame(ThreadsX.map(features(unirefs_filtered_nonzeroentropy)) do feat

    try
        ab = vec(abundances(unirefs_filtered_nonzeroentropy[feat, :]))
        df = indf[!, :]
        df.bug = Int8.(ab .> 0)
    
        ## MODEL 2
        mod = GLM.glm(
            @formula(
                bug ~
                mbiome_sample_age +
                child_sex +
                InfantVisAtt +
                MaternalEntropy +
                MaternalEntropy&InfantVisAtt
            ), df, Binomial(), ProbitLink(); dropcollinear=true
        )

        ct2 = DataFrames.Tables.rowtable(coeftable(mod))[3]
        @assert ct2.Name == "child_sex"

        ct3 = DataFrames.Tables.rowtable(coeftable(mod))[4]
        @assert ct3.Name == "InfantVisAtt"

        ct4 = DataFrames.Tables.rowtable(coeftable(mod))[5]
        @assert ct4.Name == "MaternalEntropy"
            
        ct5 = DataFrames.Tables.rowtable(coeftable(mod))[6]
        @assert ct5.Name == "MaternalEntropy & InfantVisAtt"
            
        return (
            (; ct2..., feature = string(feat), Cor = cor(df.bug, df.child_sex), kind="unirefs"),
            (; ct3..., feature = string(feat), Cor = cor(df.bug, df.InfantVisAtt), kind="unirefs"),
            (; ct4..., feature = string(feat), Cor = cor(df.bug, df.MaternalEntropy), kind="unirefs"),
            (; ct5..., feature = string(feat), Cor = cor(df.bug, (df.MaternalEntropy .* df.InfantVisAtt)), kind="unirefs")
        )

    catch
        ct2 = (;Name = "child_sex", var"Coef." = NaN, var"Std. Error" = NaN, z = NaN, var"Pr(>|z|)" = NaN, var"Lower 95%" = NaN, var"Upper 95%" = NaN)
        ct3 = (;Name = "InfantVisAtt", var"Coef." = NaN, var"Std. Error" = NaN, z = NaN, var"Pr(>|z|)" = NaN, var"Lower 95%" = NaN, var"Upper 95%" = NaN)
        ct4 = (;Name = "MaternalEntropy", var"Coef." = NaN, var"Std. Error" = NaN, z = NaN, var"Pr(>|z|)" = NaN, var"Lower 95%" = NaN, var"Upper 95%" = NaN)
        ct5 = (;Name = "MaternalEntropy & InfantVisAtt", var"Coef." = NaN, var"Std. Error" = NaN, z = NaN, var"Pr(>|z|)" = NaN, var"Lower 95%" = NaN, var"Upper 95%" = NaN)

            return (
                (; ct2..., feature = string(feat), Cor = NaN, kind="unirefs"),
                (; ct3..., feature = string(feat), Cor = NaN, kind="unirefs"),
                (; ct4..., feature = string(feat), Cor = NaN, kind="unirefs"),
                (; ct5..., feature = string(feat), Cor = NaN, kind="unirefs")
            )

        end

    end)

    for (idx, outfile) in enumerate(outfiles)

        this_df = copy(DataFrame(lmresults[!, idx]))
        select!(this_df, Cols(:feature, :Name, :))
        rename!(this_df, "Pr(>|z|)"=>"pvalue");
        sort!(this_df, :pvalue); 
    
        CSV.write(outfile, this_df)

    end

end
