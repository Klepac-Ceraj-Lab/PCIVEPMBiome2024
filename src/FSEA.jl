abstract type FSEATest end
struct Permutation <: FSEATest; nperm::Int; end
Permutation() = Permutation(1000)
struct MWU <: FSEATest end

struct FSEAResult
    pvalue::Float64
    nfeatures::Int
    setranks::Vector{Int}
end

pvalue(r::FSEAResult) = r.pvalue

function Base.show(io::IO, ::MIME"text/plain", fr::FSEAResult)
    println(io, "FSEA Result of dataset with\n",
                "  n features: ", fr.nfeatures, "\n",
                "  n in-set:  ", length(fr.setranks), "\n",
                "  p-value:   ", fr.pvalue)
end

# Helpers
function _prepare(vals, fset_idx)
    keep = findall(isfinite, vals)
    vals2 = vals[keep]
    mask = falses(length(vals)); mask[fset_idx] .= true
    fset2 = findall(mask[keep])
    return vals2, fset2
end

# Backends
function _run_fsea(::MWU, vals, fset_idx)
    vals2, fset2 = _prepare(vals, fset_idx)
    @assert 0 < length(fset2) < length(vals2)
    in_set  = collect(skipmissing(vals2[fset2]))
    out_set = collect(skipmissing(vals2[Not(fset2)]))
    mwu = MannWhitneyUTest(in_set, out_set)  # two-sided
    ranks = invperm(sortperm(vals2))
    return FSEAResult(HypothesisTests.pvalue(mwu), length(vals2), ranks[fset2])
end

function _run_fsea(perm::Permutation, vals, fset_idx)
    vals2, fset2 = _prepare(vals, fset_idx)
    @assert 0 < length(fset2) < length(vals2)
    ranks = invperm(sortperm(vals2))
    fset_ranks = ranks[fset2]
    nfeatures = length(vals2)
    es_true = enrichment_score(fset_ranks, nfeatures)
    hits = ThreadsX.count(1:perm.nperm) do _
        rranks = sort(sample(1:nfeatures, length(fset_ranks); replace=false))
        es_perm = enrichment_score(rranks, nfeatures)
        abs(es_perm) ≥ abs(es_true)
    end
    p = (hits + 1) / (perm.nperm + 1)
    return FSEAResult(p, nfeatures, fset_ranks)
end

# Front door
fsea(test::FSEATest, vals, fset_idx) = _run_fsea(test, vals, fset_idx)
fsea(vals, fset_idx) = fsea(MWU(), vals, fset_idx)

_es_at_pos(rank, idx, setscore, notscore) = idx * setscore + (rank - idx) * notscore 

function enrichment_score(ranks, nfeatures)
    nr = length(ranks)
    setscore =  -1 / nr
    notscore = 1 / (nfeatures - nr)

    ranks = sort(ranks)
    scores = ThreadsX.map(i-> _es_at_pos(ranks[i], i, setscore, notscore), eachindex(ranks))
    candidates = [scores; scores .- setscore]
    (_, i) = findmax(abs, candidates)
    return candidates[i]
end

"""
    enrichment_score(result::FSEAResult)

Calculates the enrichment score (E.S.) for a feature set
in an [`FSEAResult`](@ref).
See [Mootha et. al. (2013)](https://doi.org/10.1038/ng1180)
for details about how this is calculated.
"""
enrichment_score(result::FSEAResult) = enrichment_score(result.setranks, result.nfeatures)