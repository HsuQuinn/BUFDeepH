# band_unfolding.jl
#
# BUFDeepH v3.1.0
# Formalized from the 2025-11-19 moire unfolding experiment.
using DelimitedFiles, LinearAlgebra, JSON
using HDF5
using ArgParse
using SparseArrays
using Pardiso, Arpack, LinearMaps
using JLD
using Statistics

const BUFDeepH_VERSION = "3.1.0"
const default_dtype = Complex{Float64}
const M_TRANS_389 = [8 -9 0; 9 17 0; 0 0 1]
const M_TRANS_943 = [3 -4 0; 4 7 0; 0 0 1]

Base.@kwdef struct UnfoldingConfig
    input_dir::String
    output_dir::String
    output_file::String
    m_trans::Matrix{Int64}
    n_gvec::Int64
    max_iter::Int64
    num_band::Int64
    fermi_shift::Float64
    num_points::Int64
    mid_z::Float64
    point_start::Int64
    point_end::Int64
    part_index::Int64
    num_parts::Int64
end

function parse_m_trans(value::String)
    normalized = lowercase(strip(value))
    if normalized in ["3.89", "389", "4deg", "4degree"]
        return copy(M_TRANS_389)
    elseif normalized in ["9.43", "943", "10deg", "10degree"]
        return copy(M_TRANS_943)
    end

    matrix_text = replace(value, "[" => " ", "]" => " ")
    tokens = split(matrix_text, r"[,;\s]+", keepempty=false)
    if length(tokens) != 9
        error("--m-trans must be a preset (3.89 or 9.43) or 9 integers, e.g. '8 -9 0; 9 17 0; 0 0 1'")
    end
    return reshape(parse.(Int64, tokens), 3, 3)'
end

function parse_commandline()
    s = ArgParseSettings(
        description = "Band unfolding for moire systems. Version $(BUFDeepH_VERSION)."
    )
    @add_arg_table! s begin
        "--input_dir", "-i"
            help = "Path containing rlat.dat, orbital_types.dat, site_positions.dat, hamiltonians_pred.h5, and overlaps.h5"
            arg_type = String
            default = "../test/deephdata"
        "--output_dir", "-o"
            help = "Output directory"
            arg_type = String
            default = "./"
        "--output-file"
            help = "Output triplet data file name"
            arg_type = String
            default = "data.txt"
            dest_name = "output_file"
        "--m-trans"
            help = "Moire transformation matrix preset (3.89, 9.43) or 9 integers"
            arg_type = String
            default = "3.89"
            dest_name = "m_trans"
        "--n-gvec", "-N"
            help = "G-vector range N; loops over n,m in -N:N"
            arg_type = Int
            default = 9
            dest_name = "n_gvec"
        "--num-band"
            help = "Number of eigen bands to solve"
            arg_type = Int
            default = 5200
            dest_name = "num_band"
        "--max-iter"
            help = "Maximum Arpack iterations"
            arg_type = Int
            default = 1500
            dest_name = "max_iter"
        "--num-points"
            help = "Number of points on K-G and G-M segments"
            arg_type = Int
            default = 25
            dest_name = "num_points"
        "--fermi-shift"
            help = "Fermi-level shift in eV"
            arg_type = Float64
            default = -2.25
            dest_name = "fermi_shift"
        "--mid-z"
            help = "Only orbitals above this z coordinate contribute to the unfolded weight"
            arg_type = Float64
            default = 18.0
            dest_name = "mid_z"
        "--point-start"
            help = "1-based first k-path point to run; 0 means all/partition default"
            arg_type = Int
            default = 0
            dest_name = "point_start"
        "--point-end"
            help = "1-based last k-path point to run; 0 means all/partition default"
            arg_type = Int
            default = 0
            dest_name = "point_end"
        "--part-index"
            help = "1-based partition index for batch or Slurm array execution; 0 disables partitioning"
            arg_type = Int
            default = 0
            dest_name = "part_index"
        "--num-parts"
            help = "Total number of k-path partitions"
            arg_type = Int
            default = 1
            dest_name = "num_parts"
    end
    args = parse_args(s)
    return UnfoldingConfig(
        input_dir = args["input_dir"],
        output_dir = args["output_dir"],
        output_file = args["output_file"],
        m_trans = parse_m_trans(args["m_trans"]),
        n_gvec = args["n_gvec"],
        max_iter = args["max_iter"],
        num_band = args["num_band"],
        fermi_shift = args["fermi_shift"],
        num_points = args["num_points"],
        mid_z = args["mid_z"],
        point_start = args["point_start"],
        point_end = args["point_end"],
        part_index = args["part_index"],
        num_parts = args["num_parts"],
    )
end

function validate_config(config::UnfoldingConfig)
    config.n_gvec >= 0 || error("--n-gvec must be non-negative")
    config.max_iter > 0 || error("--max-iter must be positive")
    config.num_band > 0 || error("--num-band must be positive")
    config.num_points > 1 || error("--num-points must be greater than 1")
    config.num_parts > 0 || error("--num-parts must be positive")
    config.part_index >= 0 || error("--part-index must be non-negative")
    config.part_index <= config.num_parts || error("--part-index cannot exceed --num-parts")
    config.point_start >= 0 || error("--point-start must be non-negative")
    config.point_end >= 0 || error("--point-end must be non-negative")
    return config
end

function _create_dict_h5(filename::String)
    fid = h5open(filename, "r")
    first_key = collect(keys(fid))[1]
    T = eltype(fid[first_key])
    d_out = Dict{Array{Int64,1}, Array{T, 2}}()
    for key in keys(fid)
        data = read(fid[key])
        nk = map(x -> parse(Int64, convert(String, x)), split(key[2:length(key)-1], ','))
        d_out[nk] = permutedims(data)
    end
    close(fid)
    return d_out
end

function preprocess(input_dir::String, fermi_shift::Float64)
    if isfile(joinpath(input_dir, "info.json"))
        info = JSON.parsefile(joinpath(input_dir, "info.json"))
        spinful = info["isspinful"]
        fermi_level = info["fermi_level"] + fermi_shift
    else
        spinful = false
        fermi_level = 0.0
    end

    site_positions = readdlm(joinpath(input_dir, "site_positions.dat"))
    nsites = size(site_positions, 2)
    orbital_types_f = open(joinpath(input_dir, "orbital_types.dat"), "r")
    orbital_types = Vector{Vector{Int64}}()
    for _ in 1:nsites
        orbital_type = parse.(Int64, split(readline(orbital_types_f)))
        push!(orbital_types, orbital_type)
    end
    close(orbital_types_f)

    site_norbits = (x -> sum(x .* 2 .+ 1)).(orbital_types) * (1 + spinful)
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

    sparse_matrix_file = joinpath(input_dir, "sparse_matrix.jld")
    if isfile(sparse_matrix_file)
        @info string("read sparse matrix from ", sparse_matrix_file)
        H_R = load(sparse_matrix_file, "H_R")
        S_R = load(sparse_matrix_file, "S_R")
    else
        @info "read h5"
        begin_time = time()
        hamiltonians_pred = _create_dict_h5(joinpath(input_dir, "hamiltonians_pred.h5"))
        overlaps = _create_dict_h5(joinpath(input_dir, "overlaps.h5"))
        println("Time for reading h5: ", time() - begin_time, "s")

        I_R = Dict{Vector{Int64}, Vector{Int64}}()
        J_R = Dict{Vector{Int64}, Vector{Int64}}()
        H_V_R = Dict{Vector{Int64}, Vector{default_dtype}}()
        S_V_R = Dict{Vector{Int64}, Vector{default_dtype}}()

        @info "construct sparse matrix in the format of COO"
        begin_time = time()
        for key in collect(keys(hamiltonians_pred))
            hamiltonian_pred = hamiltonians_pred[key]
            overlap = (key in keys(overlaps)) ? overlaps[key] : zero(hamiltonian_pred)
            if spinful
                overlap = vcat(
                    hcat(overlap, zeros(size(overlap))),
                    hcat(zeros(size(overlap)), overlap)
                )
            end
            R = key[1:3]
            atom_i = key[4]
            atom_j = key[5]

            @assert (site_norbits[atom_i], site_norbits[atom_j]) == size(hamiltonian_pred)
            @assert (site_norbits[atom_i], site_norbits[atom_j]) == size(overlap)
            if !(R in keys(I_R))
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

        save(sparse_matrix_file, "H_R", H_R, "S_R", S_R)
    end
    return H_R, S_R, norbits, fermi_level, orbital_position
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
            ismutating = true
        ),
        ps
    )
end

function SolveHk(k, H_R, S_R, norbits, fermi_level, config::UnfoldingConfig)
    kx = k[1]
    ky = k[2]
    kz = k[3]
    begin_time = time()
    H_k = spzeros(default_dtype, norbits, norbits)
    S_k = spzeros(default_dtype, norbits, norbits)
    for R in keys(H_R)
        H_k += H_R[R] * exp(im * 2π * ([kx, ky, kz] ⋅ R))
        S_k += S_R[R] * exp(im * 2π * ([kx, ky, kz] ⋅ R))
    end
    lm, ps = construct_linear_map(H_k - fermi_level * S_k, S_k)
    println("Time for matrix factorization: ", time() - begin_time, "s")
    egval_inv, X = eigs(lm, nev = config.num_band, which = :LM, ritzvec = true, maxiter = config.max_iter)
    set_phase!(ps, Pardiso.RELEASE_ALL)
    eigvectors = S_k * X
    pardiso(ps)
    egval = real(1 ./ egval_inv)
    println("Time for diagonal: ", time() - begin_time, "s")
    return egval, eigvectors
end

function build_kpath(m_trans::Matrix{Int64}, num_points::Int64)
    point_K = [1 / 3, 2 / 3, 0.0]
    point_G = [0.0, 0.0, 0.0]
    point_M = [0.5, 0.5, 0.0]
    KK = m_trans' * point_K
    GG = m_trans' * point_G
    MM = m_trans' * point_M
    path_KG = [KK + t * (GG - KK) for t in range(0, 1, length = num_points)]
    path_GM = [GG + t * (MM - GG) for t in range(0, 1, length = num_points)]
    return vcat(path_KG, path_GM[2:end])
end

function partition_indices(total_points::Int64, part_index::Int64, num_parts::Int64)
    part_index > 0 || return 1, total_points
    base = div(total_points, num_parts)
    extra = rem(total_points, num_parts)
    start_idx = 1 + (part_index - 1) * base + min(part_index - 1, extra)
    stop_idx = start_idx + base - 1 + (part_index <= extra ? 1 : 0)
    return start_idx, min(stop_idx, total_points)
end

function select_point_range(config::UnfoldingConfig, total_points::Int64)
    start_idx, stop_idx = partition_indices(total_points, config.part_index, config.num_parts)
    if config.point_start > 0
        start_idx = config.point_start
    end
    if config.point_end > 0
        stop_idx = config.point_end
    end
    1 <= start_idx <= total_points || error("selected point start is outside 1:$total_points")
    start_idx <= stop_idx <= total_points || error("selected point end is outside $(start_idx):$total_points")
    return start_idx, stop_idx
end

function build_g_vectors(m_trans::Matrix{Int64}, n_gvec::Int64)
    g_vecs = []
    for n in -n_gvec:n_gvec
        for m in -n_gvec:n_gvec
            push!(g_vecs, m_trans' * [n, m, 0])
        end
    end
    return g_vecs
end

function run_unfolding(config::UnfoldingConfig)
    validate_config(config)
    mkpath(config.output_dir)

    Glat = readdlm(joinpath(config.input_dir, "rlat.dat"))
    g_vecs = build_g_vectors(config.m_trans, config.n_gvec)
    pp_points = build_kpath(config.m_trans, config.num_points)
    start_idx, stop_idx = select_point_range(config, length(pp_points))
    p_points = pp_points[start_idx:stop_idx]

    println("BUFDeepH band unfolding v$(BUFDeepH_VERSION)")
    println("Running k-path points $(start_idx):$(stop_idx) / $(length(pp_points))")
    H_R, S_R, norbits, fermi_level, orbital_position = preprocess(config.input_dir, config.fermi_shift)
    pos = orbital_position

    all_energy_levels = []
    all_weights = []

    for (local_idx, p) in enumerate(p_points)
        global_idx = start_idx + local_idx - 1
        println("Calculating local point $(local_idx), global point $(global_idx): ")
        begin_time = time()
        energy_levels = []
        total_weights = []
        val, vec = SolveHk(p, H_R, S_R, norbits, fermi_level, config)
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
                    if pos[a][3] > config.mid_z
                        pnk += exp(-1.0im * dot(k_car, pos[a])) * vec[a, n]
                    end
                end
                weight += norm(pnk)^2
            end
            push!(total_weights, weight)
        end
        println("Time for unfolding: ", time() - begin_time, "s")
        push!(all_energy_levels, energy_levels)
        push!(all_weights, total_weights)
    end

    outfile = joinpath(config.output_dir, config.output_file)
    open(outfile, "w") do io
        for (idx_p, (energy_level_avg, total_weight_sum)) in enumerate(zip(all_energy_levels, all_weights))
            for (energy, weight) in zip(energy_level_avg, total_weight_sum)
                println(io, "idx_p:$idx_p ")
                println(io, energy)
                println(io, weight)
            end
        end
    end
    println("Wrote ", outfile)
    println("over")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_unfolding(parse_commandline())
end
