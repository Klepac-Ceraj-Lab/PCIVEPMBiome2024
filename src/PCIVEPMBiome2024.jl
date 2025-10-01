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

end
