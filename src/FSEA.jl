#####
# Structs
#####

struct FSEAResult
    pvalue::Float64
    nfeatures::Int
    setranks::Vector{Int}
    es::Float64
end

abstract type FSEATest end

struct MWU <: FSEATest end

struct Permutation <: FSEATest
    nperm::Int
    seed::Int
end

Permutation(; nperm=1000, seed=1234) = Permutation(nperm, seed)

pvalue(r::FSEAResult) = r.pvalue
enrichment_score(r::FSEAResult) = r.es

function Base.show(io::IO, ::MIME"text/plain", fr::FSEAResult)
    println(io,
        "FSEA Result\n",
        "  n features: ", fr.nfeatures, "\n",
        "  n in-set:   ", length(fr.setranks), "\n",
        "  p-value:    ", fr.pvalue, "\n",
        "  E.S.:       ", fr.es
    )
end

#####
# Utilities
#####

_isvalid(x) = !ismissing(x) && isfinite(x)

function _prepare(vals, fset_idx)
    keep = findall(_isvalid, vals)

    vals2 = Float64.(vals[keep])

    mask = falses(length(vals))
    mask[fset_idx] .= true

    fset2 = findall(mask[keep])

    return vals2, fset2
end

function _validate_fset!(vals2, fset2)
    n = length(vals2)
    k = length(fset2)

    if !(0 < k < n)
        throw(ArgumentError(
            "Gene set must contain at least one valid feature and fewer than all valid features. Got $k of $n."
        ))
    end

    return nothing
end

"""
    feature_ranks(vals; high_is_enriched=true)

Return ranks of each feature after sorting `vals`.

If `high_is_enriched=true`, larger values are placed toward the top of the ranked list,
so positive ES means enrichment among high positive statistics.
"""
function feature_ranks(vals; high_is_enriched=true)
    ord = sortperm(vals; rev=high_is_enriched, alg=MergeSort)
    return invperm(ord)
end

#####
# Enrichment score
#####

function enrichment_score(ranks::AbstractVector{<:Integer}, nfeatures::Integer)
    nr = length(ranks)

    if nr == 0 || nr == nfeatures
        throw(ArgumentError("Gene set size must be between 1 and nfeatures - 1."))
    end

    setscore = 1 / nr
    notscore = -1 / (nfeatures - nr)

    hitset = falses(nfeatures)
    hitset[ranks] .= true

    score = 0.0
    maxscore = -Inf
    minscore = Inf

    @inbounds for i in 1:nfeatures
        score += hitset[i] ? setscore : notscore
        maxscore = max(maxscore, score)
        minscore = min(minscore, score)
    end

    return abs(maxscore) > abs(minscore) ? maxscore : minscore
end

#####
# Backends
#####

function _run_fsea(::MWU, vals, fset_idx; high_is_enriched=true)
    vals2, fset2 = _prepare(vals, fset_idx)
    _validate_fset!(vals2, fset2)

    inset_mask = falses(length(vals2))
    inset_mask[fset2] .= true

    in_set = vals2[inset_mask]
    out_set = vals2[.!inset_mask]

    mwu = MannWhitneyUTest(in_set, out_set)

    ranks = feature_ranks(vals2; high_is_enriched=high_is_enriched)
    fset_ranks = ranks[fset2]

    es = enrichment_score(fset_ranks, length(vals2))

    return FSEAResult(
        HypothesisTests.pvalue(mwu),
        length(vals2),
        fset_ranks,
        es
    )
end

function _run_fsea(perm::Permutation, vals, fset_idx; high_is_enriched=true)
    vals2, fset2 = _prepare(vals, fset_idx)
    _validate_fset!(vals2, fset2)

    nfeatures = length(vals2)

    ranks = feature_ranks(vals2; high_is_enriched=high_is_enriched)
    fset_ranks = ranks[fset2]

    es_true = enrichment_score(fset_ranks, nfeatures)
    k = length(fset_ranks)

    hits = ThreadsX.count(1:perm.nperm) do b
        rng = Xoshiro(perm.seed + b)
        rranks = sample(rng, 1:nfeatures, k; replace=false)
        es_perm = enrichment_score(rranks, nfeatures)
        abs(es_perm) ≥ abs(es_true)
    end

    p = (hits + 1) / (perm.nperm + 1)

    return FSEAResult(p, nfeatures, fset_ranks, es_true)
end

#####
# Front door
#####

fsea(test::FSEATest, vals, fset_idx; kwargs...) =
    _run_fsea(test, vals, fset_idx; kwargs...)

fsea(vals, fset_idx; kwargs...) =
    fsea(MWU(), vals, fset_idx; kwargs...)

#####
# High-level API
#####

function perform_fsea(
    infile::String,
    outfile::String;
    should_consolidate=false,
    test::FSEATest=MWU(),
    min_geneset_size::Int=10,
    zcol::Symbol=:z,
    featurecol::Symbol=:feature,
    pcol::Symbol=:pvalue,
    clamp_z::Union{Nothing,Tuple{Real,Real}}=(-7, 7),
    high_is_enriched::Bool=true
)
    df = CSV.read(infile, DataFrame)

    # Keep rows with valid model p-values and valid ranking statistic
    df = subset(
        df,
        pcol => ByRow(_isvalid),
        zcol => ByRow(_isvalid)
    )

    Ts = Float64.(df[!, zcol])

    if clamp_z !== nothing
        lo, hi = clamp_z
        Ts = clamp.(Ts, lo, hi)
    end

    genes = replace.(String.(df[!, featurecol]), "UniRef90_" => "")

    neuroactive = Leap.getneuroactive(
        genes;
        consolidate=should_consolidate
    )

    valid = findall(_isvalid, Ts)

    results = ThreadsX.map(collect(keys(neuroactive))) do gs
        raw_ixs = neuroactive[gs]
        ixs = intersect(raw_ixs, valid)

        if length(ixs) < min_geneset_size
            return (;
                geneset=gs,
                enrichment=NaN,
                pvalue=NaN,
                n_in_set=length(ixs),
                n_features=length(valid)
            )
        end

        try
            res = fsea(
                test,
                Ts,
                ixs;
                high_is_enriched=high_is_enriched
            )

            return (;
                geneset=gs,
                enrichment=res.es,
                pvalue=res.pvalue,
                n_in_set=length(res.setranks),
                n_features=res.nfeatures
            )
        catch err
            @warn "FSEA failed for geneset $gs" exception=(err, catch_backtrace())

            return (;
                geneset=gs,
                enrichment=NaN,
                pvalue=NaN,
                n_in_set=length(ixs),
                n_features=length(valid)
            )
        end
    end

    out = DataFrame(results)
    subset!(out, :pvalue => ByRow(_isvalid))

    out.qvalue = adjust(out.pvalue, BenjaminiHochberg())

    sort!(out, [:qvalue, :pvalue])

    CSV.write(outfile, out)

    return out
end

function perform_fsea_prevalence_sensitivity(
    infile::String,
    outfile::String;
    prevalence_col::Symbol = :prevalence,
    cutoffs = 0.0:0.05:0.50,
    should_consolidate = false,
    min_geneset_size::Int = 10,
    test::FSEATest = MWU()
)
    df = subset(CSV.read(infile, DataFrame), :pvalue => ByRow(x -> !ismissing(x) && !isnan(x)))

    @assert prevalence_col in propertynames(df) "Column $(prevalence_col) not found in input file."

    Ts_raw = map(x -> isfinite(x) ? clamp(x, -7, 7) : NaN, df.z)

    genes = map(f -> replace(f, "UniRef90_" => ""), df.feature)
    neuroactive = Leap.getneuroactive(genes; consolidate = should_consolidate)

    all_results = DataFrame()

    for cutoff in cutoffs
        pass_prev = map(x -> !ismissing(x) && x >= cutoff, df[!, prevalence_col])
        
        # Preserve original indexing so neuroactive indices remain valid.
        Ts = copy(Ts_raw)
        Ts[.!pass_prev] .= NaN

        valid = findall(!isnan, Ts)

        results = ThreadsX.map(collect(keys(neuroactive))) do gs
            ixs = neuroactive[gs]

            if isempty(ixs)
                return (;
                    prevalence_cutoff = cutoff,
                    geneset = gs,
                    enrichment = NaN,
                    pvalue = NaN,
                    nfeatures_tested = length(valid),
                    n_in_set = 0
                )
            end

            ixs_valid = intersect(ixs, valid)

            if length(ixs_valid) < min_geneset_size ||
               length(ixs_valid) >= length(valid)
                return (;
                    prevalence_cutoff = cutoff,
                    geneset = gs,
                    enrichment = NaN,
                    pvalue = NaN,
                    nfeatures_tested = length(valid),
                    n_in_set = length(ixs_valid)
                )
            end

            fes = fsea(test, Ts, ixs_valid)

            (;
                prevalence_cutoff = cutoff,
                geneset = gs,
                enrichment = fes.es,
                pvalue = fes.pvalue,
                nfeatures_tested = fes.nfeatures,
                n_in_set = length(fes.setranks)
            )
        end

        tmp = DataFrame(results)
        subset!(tmp, :pvalue => ByRow(!isnan))

        if nrow(tmp) > 0
            tmp.qvalue = adjust(tmp.pvalue, BenjaminiHochberg())
        else
            tmp.qvalue = Float64[]
        end

        append!(all_results, tmp)
    end

    sort!(all_results, [:prevalence_cutoff, :pvalue])
    CSV.write(outfile, all_results)

    return all_results
end


function plot_fsea_sensitivity(df; genesets = unique(df.geneset), type = :es)
    sub = subset(df, :geneset => ByRow(in(genesets)))

    fig = Figure(size = (900, 500))
    ax = Axis(
        fig[1, 1],
        xlabel = "Prevalence cutoff",
        ylabel = "Enrichment score",
        title = "FSEA sensitivity across prevalence filters"
    )

    # Generate distinct colors from a continuous colormap
    cmap = cgrad(:tab20, length(genesets), categorical = true)

    legend_elements = []
    legend_labels = String[]

    for (i, gs) in enumerate(genesets)
        sdf = subset(sub, :geneset => ByRow(==(gs)))
        if type == :es
            lines!(ax, sdf.prevalence_cutoff, sdf.enrichment, label = string(gs), color = cmap[i])
            scatter!(ax, sdf.prevalence_cutoff, sdf.enrichment, color = cmap[i])
        else
            lines!(ax, sdf.prevalence_cutoff, sdf.pvalue, label = string(gs), color = cmap[i])
            scatter!(ax, sdf.prevalence_cutoff, sdf.pvalue, color = cmap[i])
        end

        push!(
            legend_elements,
            LineElement(color = cmap[i], linewidth = 3)
        )
        push!(legend_labels, string(gs))

    end
    
    Legend(
        fig[1, 2],
        legend_elements,
        legend_labels;
        framevisible = false,
        tellheight = false,
        tellwidth = true
    )

    fig
end