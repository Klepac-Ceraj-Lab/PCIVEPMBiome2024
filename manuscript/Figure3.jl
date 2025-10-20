#####
# This script generates Figure 1 for PCI-VEP-Mbiome Manuscript
#####
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
using ColorSchemes
using Leap
using PCIVEPMBiome2024

Random.seed!(0)

#####
# 1. Loading data
#####
include(joinpath(Base.pwd(), "notebooks", "load_data.jl"))
