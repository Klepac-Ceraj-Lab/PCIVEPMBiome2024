## Get a NamedTuple with coef table row by term name; returns nothing if absent
function coefrow_by_term(mod, term::AbstractString)
    ctdf = DataFrame(coeftable(mod))
    # term column is usually :Name
    hasproperty(ctdf, :Name) || return nothing
    idx = findfirst(==(term), ctdf.Name)
    isnothing(idx) && return nothing
    nt = DataFrames.Tables.rowtable(ctdf)[idx]
    return nt
end

## Fit one instance of GLM (to be run by a single thread)
function fit_one(df, ab, formula::FormulaTerm, link::Link, wanted_terms::Vector{String}; feature::AbstractString)
    df = copy(df)
    df.bug = ab .> 0

    default_ret = [(; Name=t, var"Coef."=NaN, var"Std. Error"=NaN, z=NaN,
                   var"Pr(>|z|)"=NaN, var"Lower 95%"=NaN, var"Upper 95%"=NaN,
                   var"Prevalence" = mean(df.bug),
                   feature=feature, kind="unirefs") for t in wanted_terms]

    # skip degenerate responses
    if all(df.bug .== 0) || all(df.bug .== 1) || !all(map( x -> length(unique(x)), eachcol(select(df, wanted_terms[wanted_terms .∈ Ref(names(df))]))) .> 1)
        return default_ret
    end

    # fit
    mod = try 
        GLM.glm(formula, df, Binomial(), link; dropcollinear=false, maxiter=100)
    catch e
        return default_ret
    end

    outs = NamedTuple[]
    for t in wanted_terms
        r = coefrow_by_term(mod, t)
        if r === nothing
            r = (Name=t, var"Coef."=NaN, var"Std. Error"=NaN, z=NaN,
                 var"Pr(>|z|)"=NaN, var"Lower 95%"=NaN, var"Upper 95%"=NaN)
        end
        push!(outs, (; r..., var"Prevalence" = mean(df.bug), feature=feature, kind="unirefs"))
    end
    return outs
end

# Write a vector of NamedTuples to CSV with consistent ordering/renaming
function write_results(path, vec_nt::Vector{<:NamedTuple})
    df = DataFrame(vec_nt)
    # Ensure :Term exists (older GLM -> :Term), if not, try to derive from row names
    rename!(df, Symbol("Pr(>|z|)") => :pvalue)
    # Put :feature, :Name up front
    if hasproperty(df, :feature) && hasproperty(df, :Name)
        select!(df, :feature, :Name, All())
    end
    sort!(df, :pvalue, by=x->(ismissing(x) || isnan(x)) ? Inf : x)
    CSV.write(path, df)
end
