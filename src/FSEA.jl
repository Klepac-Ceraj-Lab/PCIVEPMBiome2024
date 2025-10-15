#####
# FSEA-related functions
#####

# Structs
struct FSEAResult
    pvalue::Float64
    nfeatures::Int
    setranks::Vector{Int}
    es::Float64
end

abstract type FSEATest end
struct Permutation <: FSEATest; nperm::Int; end
Permutation() = Permutation(1000)
struct MWU <: FSEATest end

# Helpers and aliases
pvalue(r::FSEAResult) = r.pvalue
enrichment_score(r::FSEAResult) = r.es

function Base.show(io::IO, ::MIME"text/plain", fr::FSEAResult)
    println(io, "FSEA Result of dataset with\n",
        "  n features: ", fr.nfeatures, "\n",
        "  n in-set:  ", length(fr.setranks), "\n",
        "  p-value:   ", fr.pvalue, "\n",
        "  E.S.:      ", fr.es)
end

function _prepare(vals, fset_idx)
    keep = findall(x -> isfinite(x) && !ismissing(x), vals)
    vals2 = vals[keep]
    mask = falses(length(vals))
    mask[fset_idx] .= true
    fset2 = findall(mask[keep])
    return vals2, fset2
end

# Enrichment score
function enrichment_score(ranks::Vector{Int}, nfeatures::Int)
    nr = length(ranks)
    setscore  = 1 / nr
    notscore  = -1 / (nfeatures - nr)

    hitset = falses(nfeatures)
    hitset[ranks] .= true

    score = 0.0
    maxscore = -Inf
    minscore = Inf
    for i in 1:nfeatures
        score += hitset[i] ? setscore : notscore
        maxscore = max(maxscore, score)
        minscore = min(minscore, score)
    end
    # ES is the max deviation from 0
    return abs(maxscore) > abs(minscore) ? maxscore : minscore
end

# Backends
function _run_fsea(::MWU, vals, fset_idx)
    vals2, fset2 = _prepare(vals, fset_idx)
    @assert 0 < length(fset2) < length(vals2)

    in_set  = vals2[fset2]
    out_set = vals2[Not(fset2)]
    mwu = MannWhitneyUTest(in_set, out_set)

    ranks = invperm(sortperm(vals2))
    es = enrichment_score(ranks[fset2], length(vals2))
    return FSEAResult(HypothesisTests.pvalue(mwu), length(vals2), ranks[fset2], es)
end

function _run_fsea(perm::Permutation, vals, fset_idx)
    vals2, fset2 = _prepare(vals, fset_idx)
    @assert 0 < length(fset2) < length(vals2)

    nfeatures = length(vals2)
    ranks = invperm(sortperm(vals2))
    fset_ranks = ranks[fset2]

    es_true = enrichment_score(fset_ranks, nfeatures)

    hits = ThreadsX.count(1:perm.nperm) do _
        rranks = sample(1:nfeatures, length(fset_ranks); replace=false)
        es_perm = enrichment_score(rranks, nfeatures)
        abs(es_perm) ≥ abs(es_true)
    end
    p = (hits + 1) / (perm.nperm + 1)
    return FSEAResult(p, nfeatures, fset_ranks, es_true)
end

# Front door
fsea(test::FSEATest, vals, fset_idx) = _run_fsea(test, vals, fset_idx)
fsea(vals, fset_idx) = fsea(MWU(), vals, fset_idx)

# High-level API
function perform_fsea(infile::String, outfile::String; should_consolidate=false)
    df = subset(CSV.read(infile, DataFrame), "pvalue" => ByRow(!isnan))

    Ts = df.z
    Ts = map(x -> isfinite(x) ? clamp(x, -7, 7) : NaN, Ts)

    genes = map(f -> replace(f, "UniRef90_" => ""), df.feature)
    neuroactive = Leap.getneuroactive(genes; consolidate=should_consolidate)

    results = ThreadsX.map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]

        if isempty(ixs) || length(ixs) < 10
            return (; geneset=gs, enrichment=NaN, pvalue=NaN)
        end

        valid = findall(!isnan, Ts)
        ixs = intersect(ixs, valid)
        if isempty(ixs)
            return (; geneset=gs, enrichment=NaN, pvalue=NaN)
        end

        fes = fsea(MWU(), Ts, ixs)
        (; geneset=gs, enrichment=fes.es, pvalue=fes.pvalue)
    end

    tmp = DataFrame(results)
    subset!(tmp, :pvalue => ByRow(!isnan))
    tmp.qvalue = adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, [:qvalue, :pvalue])
    CSV.write(outfile, tmp)
    tmp
end