module NCutSegmentation
export comp_weights, comp_segmentation, ncut

using DataStructures
using LinearAlgebra
using Statistics
using Distances
using KrylovKit
using PyCall

"""
    comp_weights(coords, pref_dir, r_max, σ_r, κ)

lorem ipsum `x = 1`
"""
function comp_weights(coords::Matrix{Int}, pref_dir::Vector{Float64}, r_max::Float64, σ_r::Float64, κ::Float64)
    spatial_dist = pairwise(SqEuclidean(), coords, dims=1)
    w_spatial = exp.(-0.5 * spatial_dist / σ_r^2)
    tuning_dist = pairwise(Minkowski(1), pref_dir[:,[CartesianIndex()]], dims=1)
    w_tuning = exp.(κ * (cos.(tuning_dist) .- 1))
    weights = w_spatial .* w_tuning
    weights[spatial_dist .> r_max^2] .= 0
    return weights
end


function comp_weights(coords::Matrix{Int}, pref_dir::Vector{Float64}, intensity::Vector{Float64}, r_max::Float64, σ_r::Float64, κ::Float64, σ_i::Float64)
    spatial_dist = pairwise(SqEuclidean(), coords, dims=1)
    w_spatial = exp.(-0.5 * spatial_dist / σ_r^2)
    tuning_dist = pairwise(Minkowski(1), pref_dir[:,[CartesianIndex()]], dims=1)
    w_tuning = exp.(κ * (cos.(tuning_dist) .- 1))
    intensity_dist = pairwise(Minkowski(1), intensity[:,[CartesianIndex()]], dims=1)
    w_intensity = exp.(-0.5 * intensity_dist / σ_i^2)
    weights = w_spatial .* w_tuning
    weights[spatial_dist .> r_max^2] .= 0
    return weights
end


function comp_segmentation(coords::Matrix{Int}, pref_dir::Vector{Float64}, r_max::Float64, σ_r::Float64, κ::Float64, min_size::Int, max_size::Int, cut_threshold::Float64)
    W = comp_weights(coords, pref_dir, r_max, σ_r, κ)
    return ncut(W, min_size, max_size, cut_threshold)
end


function comp_segmentation(coords::Matrix{Int}, pref_dir::Vector{Float64}, intensity::Vector{Float64}, r_max::Float64, σ_r::Float64, κ::Float64, σ_i::Float64, min_size::Int, max_size::Int, cut_threshold::Float64)
    W = comp_weights(coords, pref_dir, intensity, r_max, σ_r, κ, σ_i)
    return ncut(W, min_size, max_size, cut_threshold)
end


function ncut(weights::Matrix{Float64}, min_size::Int, max_size::Int, cut_threshold::Float64)
    n = size(weights)[1]
    stack = Stack{Vector{Int32}}()
    push!(stack, 1:n |> collect)
    d_component = sum(weights, dims=2)[:,1]

    m = 200
    cutvals = zeros(Float64, m)
    labels = ones(Int32, m, n)
    partition_sizes = zeros(Int, m, 2)
    label = 1
    cut_count = 1

    while ~isempty(stack)

        vertices = pop!(stack)
        n_vertices = length(vertices)
        W = weights[vertices, vertices]

        # compute Symmetric normalized Laplacian
        d = d_component[vertices]
        D = Diagonal(d)
        c = 1 ./ d
        replace!(c, Inf => 0)
        B = Diagonal(sqrt.(c))
        A = I - B * W * B

        if size(W)[1] > 1500
            # use Lanczos algorithm for large graphs
            vals, vecs, info = eigsolve(A, 2, :SR, Float64; issymmetric=true, tol=1e-6, krylovdim=50, maxiter=100, verbosity=0)
            v = vecs[2]
        else
            v = real.(eigvecs(A)[:,2])
        end

        order = sortperm(v)
        d_sorted = d[order]
        b = cumsum(d_sorted)[1:end-1] ./ cumsum(d_sorted[end:-1:1])[end-1:-1:1]

        Y = ones(Float64, n_vertices, n_vertices-1)
        @inbounds for i in 1:n_vertices-1
            Y[order[i+1:end],i] .= -b[i]
        end

        p = sum(Y .* (D * Y), dims=1)[:]
        q = sum(Y .* (W * Y), dims=1)[:]
        cuts = (p .- q) ./ p
        mincut_index = argmin(cuts)
        cut = cuts[mincut_index]

        if ((cut > cut_threshold) && (n_vertices <= max_size)) || (mincut_index < min_size && n_vertices - mincut_index < min_size)
            continue
        end

        cutvals[cut_count] = cut

        partition1 = vertices[order[1:mincut_index]]   # @ view
        partition2 = vertices[order[mincut_index+1:n_vertices]]
        #partitions = Array[partition1, partition2]
        partition_sizes[cut_count, 1] = mincut_index
        partition_sizes[cut_count, 2] = n_vertices - mincut_index

        cut_count += 1
        labels[cut_count, :] .= labels[cut_count-1, :]

        if length(partition1) < min_size
            labels[cut_count, partition1] .= 0
            push!(stack, partition2)
        elseif length(partition2) < min_size
            labels[cut_count, partition2] .= 0
            push!(stack, partition1)
        else
            labels[cut_count, partition1] .= maximum(labels[cut_count, :]) + 1
            push!(stack, partition1)
            push!(stack, partition2)
        end

        if cut_count - 1 > m
            cutvals = vcat(cutvals, zeros(Float64, m))
            labels = vcat(labels, ones(Int, m, n))
            partition_sizes = vcat(partition_sizes, zeros(Int, m, 2))
        end
    end

    return labels[1:cut_count, :], cutvals[1:cut_count-1], partition_sizes[1:cut_count-1,:]
end

end
