module PCIVEPMBiome2024

    using GLM
    using CSV
    using DataFrames
    using Chain
    using MultivariateStats
    using HypothesisTests
    using MultipleTesting
    using Distributions
    using CategoricalArrays
    using Microbiome
    using Distances
    using Dictionaries
    using Clustering
    using Diversity
    using DecisionTree
    using Random
    using JLD2
    using Statistics
    using CairoMakie
    using ThreadsX
    using GLM
    using CodecZlib
    using Leap

    include("utils.jl")

    export load_custom_unirefs

    include("gene_glm_src.jl")

    export coefrow_by_term
    export fit_one
    export write_results

    include("runlms.jl")
    
    export runlms
    export runlms_nosex

    include("FSEA.jl")

    export FSEATest
    export FSEAResult
    export Permutation
    export MWU
    export pvalue
    export fsea
    export enrichment_score
    export perform_fsea

end
