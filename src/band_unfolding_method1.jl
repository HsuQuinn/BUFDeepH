# band_unfolding.jl version 0.1.0
# only for honeycomb lattice now, but can be easily modified to other lattices, in _find function.
# purpose for moire band unfolding, maybe exist some bugs for other situations.
# Hsu Quinn
using DelimitedFiles, LinearAlgebra, JSON
using HDF5
using ArgParse
using SparseArrays
using Pardiso, Arpack, LinearMaps
using JLD
using PyPlot
using Statistics

# Const
default_dtype = Complex{Float64}


M_trans = [1 -2 0; 2 3 0; 0 0 1]        # To do : Read M
N = 2                                   # number of g-vectors in |p+g> , equal to 2(N+1)^2
max_iter = 400
num_band = 400                         # bands num in LC
fermi_level_shift = 0.0                 # fermi level shift for deep band
num_points = 17
mid_z_coord = 19.0                      # z coordinate of the middle plane, a lazy input for developer.
y_down = -7                             # y axis down limit
y_up   = 0                              # y axis up limit


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--input_dir", "-i"
            help = "path of rlat.dat, orbital_types.dat, site_positions.dat, hamiltonians_pred.h5, and overlaps.h5"
            arg_type = String
            default = "../test/deephdata"
        "--output_dir", "-o"
            help = "path of output .png file"
            arg_type = String
            default = "./"
    end
    return parse_args(s)
end
parsed_args = parse_commandline()

function preprocess()
    if isfile(joinpath(parsed_args["input_dir"],"info.json"))
        spinful = JSON.parsefile(joinpath(parsed_args["input_dir"],"info.json"))["isspinful"]
        fermi_level = JSON.parsefile(joinpath(parsed_args["input_dir"],"info.json"))["fermi_level"] + fermi_level_shift
    else
        spinful = false
        fermi_level = 0.0
    end
    site_positions = readdlm(joinpath(parsed_args["input_dir"], "site_positions.dat"))
    nsites = size(site_positions, 2)
    orbital_types_f = open(joinpath(parsed_args["input_dir"], "orbital_types.dat"), "r")
    site_norbits = zeros(nsites)
    orbital_types = Vector{Vector{Int64}}()
    for index_site = 1:nsites
        orbital_type = parse.(Int64, split(readline(orbital_types_f)))
        push!(orbital_types, orbital_type)
    end
    site_norbits = (x->sum(x .* 2 .+ 1)).(orbital_types) * (1 + spinful)
    norbits = sum(site_norbits)
    site_norbits_cumsum = cumsum(site_norbits)
    orbital_position = []
    for i in 1:nsites
        site_pos = site_positions[:, i]
        num_orbitals = sum(orbital_types[i] .* 2 .+ 1)
        for _ in 1:num_orbitals
            push!(orbital_position, site_pos)
        end
    end
    if isfile(joinpath(parsed_args["input_dir"], "sparse_matrix.jld"))
        @info string("read sparse matrix from ", parsed_args["input_dir"], "/sparse_matrix.jld")
        H_R = load(joinpath(parsed_args["input_dir"], "sparse_matrix.jld"), "H_R")
        S_R = load(joinpath(parsed_args["input_dir"], "sparse_matrix.jld"), "S_R")
    else
        @info "read h5"
        begin_time = time()
        hamiltonians_pred = _create_dict_h5(joinpath(parsed_args["input_dir"], "hamiltonians_pred.h5"))
        overlaps = _create_dict_h5(joinpath(parsed_args["input_dir"], "overlaps.h5"))
        println("Time for reading h5: ", time() - begin_time, "s")

        I_R = Dict{Vector{Int64}, Vector{Int64}}()
        J_R = Dict{Vector{Int64}, Vector{Int64}}()
        H_V_R = Dict{Vector{Int64}, Vector{default_dtype}}()
        S_V_R = Dict{Vector{Int64}, Vector{default_dtype}}()

        @info "construct sparse matrix in the format of COO"
        begin_time = time()
        for key in collect(keys(hamiltonians_pred))
            hamiltonian_pred = hamiltonians_pred[key]
            if (key ∈ keys(overlaps))
                overlap = overlaps[key]
            else
                # continue
                overlap = zero(hamiltonian_pred)
            end
            if spinful
                overlap = vcat(hcat(overlap,zeros(size(overlap))),hcat(zeros(size(overlap)),overlap)) # the readout overlap matrix only contains the upper-left block # TODO maybe drop the zeros?
            end
            R = key[1:3]; atom_i=key[4]; atom_j=key[5]

            @assert (site_norbits[atom_i], site_norbits[atom_j]) == size(hamiltonian_pred)
            @assert (site_norbits[atom_i], site_norbits[atom_j]) == size(overlap)
            if !(R ∈ keys(I_R))
                I_R[R] = Vector{Int64}()
                J_R[R] = Vector{Int64}()
                H_V_R[R] = Vector{default_dtype}()
                S_V_R[R] = Vector{default_dtype}()
            end
            for block_matrix_i in 1:site_norbits[atom_i]
                for block_matrix_j in 1:site_norbits[atom_j]
                    coo_i = site_norbits_cumsum[atom_i] - site_norbits[atom_i] + block_matrix_i
                    coo_j = site_norbits_cumsum[atom_j] - site_norbits[atom_j] + block_matrix_j
                    push!(I_R[R], coo_i)
                    push!(J_R[R], coo_j)
                    push!(H_V_R[R], hamiltonian_pred[block_matrix_i, block_matrix_j])
                    push!(S_V_R[R], overlap[block_matrix_i, block_matrix_j])
                end
            end
        end
        println("Time for constructing sparse matrix in the format of COO: ", time() - begin_time, "s")

        @info "convert sparse matrix to the format of CSC"
        begin_time = time()
        H_R = Dict{Vector{Int64}, SparseMatrixCSC{default_dtype, Int64}}()
        S_R = Dict{Vector{Int64}, SparseMatrixCSC{default_dtype, Int64}}()

        for R in keys(I_R)
            H_R[R] = sparse(I_R[R], J_R[R], H_V_R[R], norbits, norbits)
            S_R[R] = sparse(I_R[R], J_R[R], S_V_R[R], norbits, norbits)
        end
        println("Time for converting to the format of CSC: ", time() - begin_time, "s")

        save(joinpath(parsed_args["input_dir"], "sparse_matrix.jld"), "H_R", H_R, "S_R", S_R)
    end
    return H_R, S_R, norbits, fermi_level, orbital_position
end

function _create_dict_h5(filename::String)
    fid = h5open(filename, "r")
    T = eltype(fid[keys(fid)[1]])
    d_out = Dict{Array{Int64,1}, Array{T, 2}}()
    for key in keys(fid)
        data = read(fid[key])
        nk = map(x -> parse(Int64, convert(String, x)), split(key[2 : length(key) - 1], ','))
        d_out[nk] = permutedims(data)
    end
    close(fid)
    return d_out
end

function construct_linear_map(H, S)
    ps = MKLPardisoSolver()
    set_matrixtype!(ps, Pardiso.COMPLEX_HERM_INDEF)
    pardisoinit(ps)
    fix_iparm!(ps, :N)
    H_pardiso = get_matrix(ps, H, :N)
    b = rand(ComplexF64, size(H, 1))
    set_phase!(ps, Pardiso.ANALYSIS)
    pardiso(ps, H_pardiso, b)
    set_phase!(ps, Pardiso.NUM_FACT)
    pardiso(ps, H_pardiso, b)
    return (
        LinearMap{ComplexF64}(
            (y, x) -> begin
                set_phase!(ps, Pardiso.SOLVE_ITERATIVE_REFINE)
                pardiso(ps, y, H_pardiso, S * x)
            end,
            size(H, 1);
            ismutating=true
        ),
        ps
    )
end

function SolveHk(k)
    kx = k[1]
    ky = k[2]
    kz = k[3]
    H_R, S_R, norbits, fermi_level, orbital_position = preprocess()
    begin_time = time()
    H_k = spzeros(default_dtype, norbits, norbits)
    S_k = spzeros(default_dtype, norbits, norbits)
    for R in keys(H_R)
        H_k += H_R[R] * exp(im*2π*([kx, ky, kz]⋅R))
        S_k += S_R[R] * exp(im*2π*([kx, ky, kz]⋅R))
    end
    lm, ps = construct_linear_map(H_k - (fermi_level) * S_k, S_k)
    println("Time for matrix factorization: ", time() - begin_time, "s")
    egval_inv, X = eigs(lm, nev=num_band, which=:LM, ritzvec=true, maxiter=max_iter)
    set_phase!(ps, Pardiso.RELEASE_ALL)
    eigvectors = S_k * X
    pardiso(ps)
    egval = real(1 ./ egval_inv)
    println("Time for diagonal: ", time() - begin_time, "s")
    return egval, eigvectors, orbital_position
end



Glat = readdlm(joinpath(parsed_args["input_dir"], "rlat.dat"))
glat = M_trans' * Glat
g_vecs = []
for n in -N:N
    for m in -N:N
        g = M_trans' * [n, m, 0]
        push!(g_vecs, g)
    end
end

point_K = [1/3,   2/3,    0.0]
point_G = [0.0,   0.0,    0.0]
point_M = [0.5,   0.5,    0.0]
KK = M_trans' * point_K
GG = M_trans' * point_G
MM = M_trans' * point_M
path_KG = [KK + t * (GG - KK) for t in range(0, 1, length=num_points)]
path_GM = [GG + t * (MM - GG) for t in range(0, 1, length=num_points)]
p_points = vcat(path_KG, path_GM[2:end]) 
#p_points = [[0,0,0]]

all_energy_levels = []
all_weights = []

for (idx_p, p) in enumerate(p_points)
    println("Calculating No$idx_p p-points: ")
    begin_time = time()
    energy_levels = []
    total_weights = []
    val, vec, pos = SolveHk(p)
    for n in 1:size(val)[1]
        push!(energy_levels, val[n])  
    end
    for n in 1:size(val)[1]
        weight = 0
        for g in g_vecs
            k = p + g 
            pnk = 0
            for a in 1:size(pos)[1]
                k_car = (k' * Glat')'
                if pos[a][3] > mid_z_coord
                    pnk += exp(-1.0im * dot(k_car, pos[a])) * vec[a,n]
                end
            end
            weight += norm(pnk)^2 
        end
        push!(total_weights, weight)
    end
    total_weights = (total_weights .- minimum(total_weights)) ./ (maximum(total_weights) - minimum(total_weights))
    println("Time for unfolding: ", time() - begin_time, "s")
    push!(all_energy_levels, energy_levels)
    push!(all_weights, total_weights)
end


figure(figsize=(6, 10))
for (idx_p, (energy_level_avg, total_weight_sum)) in enumerate(zip(all_energy_levels, all_weights))
    for (energy, weight) in zip(energy_level_avg, total_weight_sum)
        scatter(idx_p, energy, color="red", s=weight * 10, alpha=0.6)  
        println("idx_p:$idx_p ")
        println(energy)
        println(weight)
    end
end

xlabel("K->G->M")
ylabel("Energy Level")
title("band_unfolding")
xlim(0, length(p_points) + 1)
ylim(y_down, y_up)
savefig("BandUnFolDing_method1.png")
show()
close()