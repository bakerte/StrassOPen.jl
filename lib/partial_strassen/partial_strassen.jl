#########################################################################
#
#           Strassen Open Implementation (StrassOPen)
#                          v1.0
#
#########################################################################
# Made by Thomas E. Baker, Kiana Gallagher, and Aaron Dayton and « les qubits volants » (2024)
# See accompanying license with this program
# This code is native to the julia programming language (v1.10.4+)
#

#include("partial_strassOPen.jl")
using LinearAlgebra
using BenchmarkTools

# SINGLE THREADING
#LinearAlgebra.BLAS.set_num_threads(1)


function init_subblock!(
    Z::AbstractMatrix, A::AbstractMatrix,
    A_row_inds::UnitRange, A_col_inds::UnitRange,
	sign::Bool, region_end_row::Integer, region_end_col::Integer
)
    nrows = length(A_row_inds)
    ncols = length(A_col_inds)


	if sign # +
		@inbounds for j = 1:ncols
			jA = A_col_inds[j]

			@simd for i = 1:nrows
				iA = A_row_inds[i]

				Z[i, j] = A[iA, jA]
			end
		end
	else # -
		@inbounds for j = 1:ncols
			jA = A_col_inds[j]

			@simd for i = 1:nrows
				iA = A_row_inds[i]

				Z[i, j] = -A[iA, jA]
			end
		end
	end

	# Padding:

	# Bottom padding
	for j = 1 : region_end_col
		for i = nrows+1 : region_end_row
			Z[i, j] = 0
		end
	end

	# Side padding
	for j = ncols+1 : region_end_col
		for i = 1 : region_end_row 
			Z[i, j] = 0
		end
	end

	# Bottom corner padding
	for j = ncols+1 : region_end_col
		for i = nrows+1 : region_end_row
			Z[i, j] = 0
		end
	end
end


function add_subblock!(
    Z::AbstractMatrix, A::AbstractMatrix,
    A_row_inds::UnitRange, A_col_inds::UnitRange,
	sign::Bool
)
    nrows = length(A_row_inds)
    ncols = length(A_col_inds)

	if sign # +
		@inbounds for j = 1:ncols
			jA = A_col_inds[j]

			@simd for i = 1:nrows
				iA = A_row_inds[i]

				Z[i, j] += A[iA, jA]
			end
		end
	else # -
		@inbounds for j = 1:ncols
			jA = A_col_inds[j]

			@simd for i = 1:nrows
				iA = A_row_inds[i]

				Z[i, j] -= A[iA, jA]
			end
		end
	end
end


function mul_subblocks!(
    C::AbstractMatrix, C_row_inds::UnitRange, C_col_inds::UnitRange,
    A::AbstractMatrix, A_row_inds::UnitRange, A_col_inds::UnitRange,
    B::AbstractMatrix, B_row_inds::UnitRange, B_col_inds::UnitRange
)
	@views mul!(C[C_row_inds, C_col_inds], A[A_row_inds, A_col_inds], B[B_row_inds, B_col_inds])
end


# Handle the C buffer shape:
# There are 3 regions to the C buffer
#  _____ _____             _____ _____ _____
# | 1   | 2   |           | 1   | 2   | 3   |
# |_____|_____|   ---->   |_____|_____|_____|
# | 3   |
# |_____|
#
# We must index it as depicted on the right above,
# so elements in region 3 of the buffer must be shifted to the side accordingly.
# The second boolean of the if statements corresponds to checking this condition
function find_buffer_indices(eq, level, outerA_dict, outerB_dict, recombination_index=0)
    #########
    # NOTE: #
    #########
    # We use dictionary starts below to avoid any issues that may arise from the dictionary ranges ending prematurely
    # as they were initialized for matrices A and B, such that they cannot go out of bounds,
    # but that may be a location required when indexing C or Cbuffer
    if recombination_index == 0
        if eq.buffer_bool && ( outerA_dict[(eq.C_block.row_block, level)].start > length(outerA_dict[(0, 1)]) )

            C_row_start = ( outerA_dict[(eq.C_block.row_block, level)].start - length(outerA_dict[(0, 1)]))
            C_rows = C_row_start : (C_row_start + length(outerA_dict[(0, level)]) - 1)

            C_col_start = 2 * length(outerB_dict[(0, 1)]) + outerB_dict[(eq.C_block.col_block, level)].start
            C_cols = C_col_start : (C_col_start  + length(outerB_dict[(0, level)]) - 1)

            return C_rows, C_cols
        else

            C_rows = (outerA_dict[(eq.C_block.row_block, level)].start) : (outerA_dict[(eq.C_block.row_block, level)].start + length(outerA_dict[(0, level)]) - 1)
            C_cols = (outerB_dict[(eq.C_block.col_block, level)].start) : (outerB_dict[(eq.C_block.col_block, level)].start + length(outerB_dict[(0, level)]) - 1)

            return C_rows, C_cols
        end

    else # recombination_index != 0
        # This is the recombination case for the recombination step where 7 matrices get combined
        # This only gets called explicitly for buffer matrixes!!!
        if ( outerA_dict[(eq.recombination_blocks[recombination_index].row_block, level)].start > length(outerA_dict[(0, 1)]) )

            C_row_start = ( (outerA_dict[(eq.recombination_blocks[recombination_index].row_block, level)].start - length(outerA_dict[(0, 1)])))
            C_rows = C_row_start : (C_row_start + length(outerA_dict[(0, level)]) - 1)

            C_col_start = 2 * length(outerB_dict[(0, 1)]) + outerB_dict[(eq.recombination_blocks[recombination_index].col_block, level)].start
            C_cols = C_col_start : (C_col_start  + length(outerB_dict[(0, level)]) - 1)

            return C_rows, C_cols
        else

            C_rows = (outerA_dict[(eq.recombination_blocks[recombination_index].row_block, level)].start) : (outerA_dict[(eq.recombination_blocks[recombination_index].row_block, level)].start + length(outerA_dict[(0, level)]) - 1)
            C_cols = (outerB_dict[(eq.recombination_blocks[recombination_index].col_block, level)].start) : (outerB_dict[(eq.recombination_blocks[recombination_index].col_block, level)].start + length(outerB_dict[(0, level)]) - 1)

            return C_rows, C_cols
        end
    end
end


function add_and_multiply_A_and_B_submatrices!(
    A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix,
    Abuffer::AbstractMatrix, Bbuffer::AbstractMatrix, Cbuffer::AbstractMatrix,
    eq::Equation,
	outerA_dict::Dict{Tuple{Int64, Int64}, UnitRange{Int64}},
	inner_dict::Dict{Tuple{Int64, Int64}, UnitRange{Int64}},
	outerB_dict::Dict{Tuple{Int64, Int64}, UnitRange{Int64}}
)
    # The level must be the same for all submatrices upon the function call
	level = eq.A_ops[1].arg.level

	s = 1
	for op in eq.A_ops

		if s == 1
			# Replace elements in buffer matrix
			init_subblock!(
				Abuffer,
				A,
                outerA_dict[(op.arg.row_block, level)],
                inner_dict[(op.arg.col_block, level)],
                op.func,
                length(outerA_dict[(0, level)]),
                length(inner_dict[(0, level)])
			)
		else
			# Add to elements in buffer matrix
			add_subblock!(
				Abuffer,
				A,
                outerA_dict[(op.arg.row_block, level)],
                inner_dict[(op.arg.col_block, level)],
                op.func
			)
		end
	s += 1
	end

	s = 1
	for op in eq.B_ops

		if s == 1
			# Replace elements in buffer matrix
			init_subblock!(
				Bbuffer,
				B,
                inner_dict[(op.arg.row_block, level)],
                outerB_dict[(op.arg.col_block, level)],
                op.func,
                length(inner_dict[(0, level)]),
                length(outerB_dict[(0, level)])
			)
		else
			# Add to elements in buffer matrix
			add_subblock!(
				Bbuffer,
				B,
                inner_dict[(op.arg.row_block, level)],
                outerB_dict[(op.arg.col_block, level)],
                op.func
			)
		end
		s += 1
	end




    C_rows, C_cols = find_buffer_indices(eq, level, outerA_dict, outerB_dict)

	# Matrix multiply
	mul_subblocks!(
		eq.buffer_bool ? Cbuffer : C,
        C_rows,
        C_cols,
		Abuffer,
        outerA_dict[(0, level)],
        inner_dict[(0, level)],
		Bbuffer,
        inner_dict[(0, level)],
        outerB_dict[(0, level)]
	)
	return
end



function add_M_submatrices!(
    C::AbstractMatrix,
    Cbuffer::AbstractMatrix,
    eq::Equation,
	n_dict::Dict{Tuple{Int64, Int64}, UnitRange{Int64}},
	m_dict::Dict{Tuple{Int64, Int64}, UnitRange{Int64}}
)
    # All levels must be the same
    level = eq.recombination_blocks[1].level

    #########
    # NOTE: #
    #########
    # We use dictionary starts below to avoid any issues that may arise from the dictionary ranges ending prematurely
    # as they were initialized for matrices A and B, such that they cannot go out of bounds,
    # but that may be a location required when indexing C or Cbuffer
    rl = length(n_dict[(0, level)])-1
    cl = length(m_dict[(0, level)])-1

    # M1
    M1_rows, M1_cols = find_buffer_indices(eq, level, n_dict, m_dict, 1)

    # M2
    M2_rows = (n_dict[(eq.recombination_blocks[2].row_block, level)].start) : (n_dict[(eq.recombination_blocks[2].row_block, level)].start + rl)
    M2_cols = (m_dict[(eq.recombination_blocks[2].col_block, level)].start) : (m_dict[(eq.recombination_blocks[2].col_block, level)].start + cl)
    
    # M3
    M3_rows, M3_cols = find_buffer_indices(eq, level, n_dict, m_dict, 3)

    # M4
    M4_rows = (n_dict[(eq.recombination_blocks[4].row_block, level)].start) : (n_dict[(eq.recombination_blocks[4].row_block, level)].start + rl)
    M4_cols = (m_dict[(eq.recombination_blocks[4].col_block, level)].start) : (m_dict[(eq.recombination_blocks[4].col_block, level)].start + cl)

    # M5
    M5_rows = (n_dict[(eq.recombination_blocks[5].row_block, level)].start) : (n_dict[(eq.recombination_blocks[5].row_block, level)].start + rl)
    M5_cols = (m_dict[(eq.recombination_blocks[5].col_block, level)].start) : (m_dict[(eq.recombination_blocks[5].col_block, level)].start + cl)

    # M6
    M6_rows = (n_dict[(eq.recombination_blocks[6].row_block, level)].start) : (n_dict[(eq.recombination_blocks[6].row_block, level)].start + rl)
    M6_cols = (m_dict[(eq.recombination_blocks[6].col_block, level)].start) : (m_dict[(eq.recombination_blocks[6].col_block, level)].start + cl)

    # M7
    M7_rows, M7_cols = find_buffer_indices(eq, level, n_dict, m_dict, 7)

    row_inds = n_dict[(0, level)]
    col_inds = m_dict[(0, level)]
	@inbounds for j in col_inds
		@simd for i in row_inds
			# M6 += M1 - M2 + M3
			C[M6_rows[i], M6_cols[j]] += Cbuffer[M1_rows[i], M1_cols[j]] - C[M2_rows[i], M2_cols[j]] +  Cbuffer[M3_rows[i], M3_cols[j]]
			# M2 += M4
			C[M2_rows[i], M2_cols[j]] += C[M4_rows[i], M4_cols[j]]
			# M4 += M1 - M5 + M7
			C[M4_rows[i], M4_cols[j]] += Cbuffer[M1_rows[i], M1_cols[j]] - C[M5_rows[i], M5_cols[j]] + Cbuffer[M7_rows[i], M7_cols[j]]
			# M5 += M3
			C[M5_rows[i], M5_cols[j]] += Cbuffer[M3_rows[i], M3_cols[j]]
		end
	end
end



function initialize_block_map_dictionaries(n, k, m, depth)
	# Have a dictionary for the outer dimension of A (n),
	# inner dimension (k), and outer dimension of B (m)
	#
	#        Dict{Tuple{row/col, level}, corresponding unit range}
	n_dict = Dict{Tuple{Int64, Int64}, UnitRange{Int64}}()
	k_dict = Dict{Tuple{Int64, Int64}, UnitRange{Int64}}()
	m_dict = Dict{Tuple{Int64, Int64}, UnitRange{Int64}}()

	for level in 1:depth
		curr_n = cld(n, (2^depth)) * 2^(depth-level)
		curr_k = cld(k, (2^depth)) * 2^(depth-level)
		curr_m = cld(m, (2^depth)) * 2^(depth-level)

        # Put base block size into dict
        n_dict[(0, level)] = 1:curr_n
        k_dict[(0, level)] = 1:curr_k
        m_dict[(0, level)] = 1:curr_m


		for c in 1:2^level
            # if-else statements make sure that we do not go out of bounds
            if ( curr_n * c ) < n
			    n_dict[(c, level)] = ( curr_n * (c-1) + 1 ) : ( curr_n * c )
            else
			    n_dict[(c, level)] = ( curr_n * (c-1) + 1 ) : ( n )
            end

            if ( curr_k * c ) < k
			    k_dict[(c, level)] = ( curr_k * (c-1) + 1 ) : ( curr_k * c )
            else
			    k_dict[(c, level)] = ( curr_k * (c-1) + 1 ) : ( k )
            end

            if ( curr_m * c ) < m
			    m_dict[(c, level)] = ( curr_m * (c-1) + 1 ) : ( curr_m * c )
            else
			    m_dict[(c, level)] = ( curr_m * (c-1) + 1 ) : ( m )
            end
		end
	end

	return n_dict, k_dict, m_dict
end





function partial_strassen(A::AbstractMatrix, B::AbstractMatrix, eqs::Vector{Equation}, n::Int64)
    # Use a dictionary for indexing the matrices
    outerA_dict::Dict{Tuple{Int64, Int64}, UnitRange{Int64}},
        inner_dict::Dict{Tuple{Int64, Int64}, UnitRange{Int64}},
        outerB_dict::Dict{Tuple{Int64, Int64}, UnitRange{Int64}} =
            initialize_block_map_dictionaries(size(A,1), size(A,2), size(B,2), n)

	# If mod(size(A,1), 2^n) == 0, we do not need padding,
	# otherwise pad so that the matrix dimensions are divisible by 2^n
	padding_rows = mod(size(A,1), 2^n) == 0 ? 0 : 2^n - mod(size(A,1), 2^n)
	padding_cols = mod(size(B,2), 2^n) == 0 ? 0 : 2^n - mod(size(B,2), 2^n)

	# Intialize C and use it as a buffer throughout 
	C = Array{Float64}(undef, size(A,1)+padding_rows, size(B,2)+padding_cols)

	# In single threaded environment this is all we need
	Abuffer = Array{Float64}(undef, cld(size(A,1), 2^n)*2^(n-1), cld(size(A,2), 2^n)*2^(n-1))
	Bbuffer = Array{Float64}(undef, cld(size(B,1), 2^n)*2^(n-1), cld(size(B,2), 2^n)*2^(n-1))
	Cbuffer = Array{Float64}(undef, cld(size(A,1), 2^n)*2^(n-1), 3*cld(size(B,2), 2^n)*2^(n-1))

	for eq in eqs
		add_and_multiply_A_and_B_submatrices!(A, B, C, Abuffer, Bbuffer, Cbuffer, eq, outerA_dict, inner_dict, outerB_dict)

		if eq.recombine
			add_M_submatrices!(C, Cbuffer, eq, outerA_dict, outerB_dict)
		end
	end
	return view(C, 1:size(A,1), 1:size(B,2))
end


# Benchmark times:
#	square_size = 2^14
#		BLAS: 288.103591 seconds (3 allocations: 2.000 GiB, 0.04% gc time)
#		PS, n=1: 256.748429 seconds (12 allocations: 4.500 GiB, 0.01% gc time)
#		PS, n=2: 253.987697 seconds (12 allocations: 4.500 GiB, 0.03% gc time)
#		PS, n=3: 248.663672 seconds (12 allocations: 4.500 GiB, 0.02% gc time)
function test(n::Int64, outer_A::Int64, inner::Int64, outer_B::Int64)

	A::AbstractMatrix = rand(outer_A, inner)
	B::AbstractMatrix = rand(inner, outer_B)

	# A::AbstractMatrix = zeros(outer_A, inner)
	# B::AbstractMatrix = zeros(inner, outer_B)
    # for i in 1:length(A)
    #     A[i] = i
    # end
    # for i in 1:length(B)
    #     B[i] = i
    # end


    eqs::Vector{Equation} = []
	create_partial_strassen_equations(eqs, n)
    result = partial_strassen(A, B, eqs, n)
    @btime partial_strassen($A, $B, $eqs, $n)
    @btime $A*$B

	println(norm(result - A*B));
end
# To run test:
# 	test(2, 2^10, 2^10, 2^10)




function testmany(ns, range, time=false)

	# For precompilation
    temp_eqs::Vector{Equation} = []
	create_partial_strassen_equations(temp_eqs, 1)
	result = partial_strassen(rand(10,10), rand(10,10), temp_eqs, 1);

	for n in ns
		for outerA in range
			for inner in range
				for outerB in range
					A = rand(outerA, inner)
					B = rand(inner, outerB)

                    eqs::Vector{Equation} = []
					create_partial_strassen_equations(eqs, n)
					result = partial_strassen(A, B, eqs, n)
					diff = norm(result - A*B)
					if !isapprox(diff, 0.0, atol=1e-8)
						println("ERROR:")
						println("  diff   = ", diff)
						println("  n      = ", n)
						println("  outerA = ", outerA)
						println("  inner  = ", inner)
						println("  outerB = ", outerB)
						println()
					end

					if time
						println("BLAS:")
						@time A*B
						println("PS:")
						@time partial_strassen(A, B, eqs, n);
					end
				end
			end
		end
	end
end
