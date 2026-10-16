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

mutable struct PSBlock
	row_block :: Int64
	col_block :: Int64
	level     :: Int64
end



# The Operation struct has:
# 	- an operation to be performed (true = +) ; (false = -)
#	- a location in memory to be used as the argument in the operation
#
# This struct is used for matrix storage, addition, and subtraction
mutable struct Operation
	func :: Bool
	arg  :: PSBlock
end



struct Equation
  A_ops                ::Vector{Operation}
  B_ops                ::Vector{Operation}
  C_block              ::PSBlock
  buffer_bool          ::Bool
  recombine            ::Bool
  recombination_blocks :: Vector{PSBlock}
end






#####################
# Helper Functions  #
#####################

# Rotate the L-shaped buffer block regions
# This code is verbose, but it works.
function rotate_blocks!(M1::PSBlock, M3::PSBlock, M7::PSBlock, curr_rotation::Int64, prev_rotation::Int64)
	if prev_rotation == 0
		if curr_rotation == 90
			## _____      _____
			# |M1 M3| -> |M7 M1|
			# |M7|          |M3|
			##
			M3.row_block = M7.row_block
			M1.col_block = M3.col_block
			M7.row_block = M1.row_block
		elseif curr_rotation == -90
			## _____      __
			# |M1 M3| -> |M3|__
			# |M7|       |M1 M7|
			##
			M7.col_block = M3.col_block
			M1.row_block = M7.row_block
			M3.col_block = M1.col_block
		end
	elseif prev_rotation == 90
		if curr_rotation == 90
			## _____         __
			# |M7 M1| ->  __|M7|
			#    |M3|    |M3 M1|
			##
			M3.col_block = M7.col_block
			M1.row_block = M3.row_block
			M7.col_block = M1.col_block
		elseif curr_rotation == -90
			## _____      _____ 
			# |M7 M1| -> |M1 M3|
			#    |M3|	    |M7|
			##
			M7.row_block = M3.row_block
			M1.col_block = M7.col_block
			M3.row_block = M1.row_block
		end
	elseif prev_rotation == 180
		if curr_rotation == 90
			##    __      __
			#  __|M7| -> |M3|__
			# |M3 M1|    |M1 M7|
			##
			M3.row_block = M7.row_block
			M1.col_block = M3.col_block
			M7.row_block = M1.row_block
		elseif curr_rotation == -90
			##    __      _____
			#  __|M7| -> |M7 M1|
			# |M3 M1|       |M3|
			##
			M7.col_block = M3.col_block
			M1.row_block = M7.row_block
			M3.col_block = M1.col_block
		end

	elseif prev_rotation == 270
		if curr_rotation == 90
			## __         _____
			# |M3|__  -> |M1 M3|
			# |M1 M7|    |M7|
			##
			M3.col_block = M7.col_block
			M1.row_block = M3.row_block
			M7.col_block = M1.col_block
		elseif curr_rotation == -90
			## __            __
			# |M3|__  ->  __|M7|
			# |M1 M7|    |M3 M1|
			##
			M7.row_block = M3.row_block
			M1.col_block = M7.col_block
			M3.row_block = M1.row_block
		end
	end
end


# Apply a general map to M1, M3, and M7 in the ***STANDARD ORIENTATION***
# 
#  ___________
# |M1   |M3   |
# |_____|_____|
# |M7   |
# |_____|
# If M1, M3, and M7 are rotated, undo the rotation such that they
# are in the standard orientation, apply this map,
# then rotate the result back to get your answer.
function master_map(M1, M3, M7, coord)
	new_level = M1.level+1
	if coord == 1
		#  ___________
		# |11 12|__   |
		# |21|   __|__|
		# |  |__|
		# |_____|
		block_11 = PSBlock(2*M1.row_block-1, 2*M1.col_block-1, new_level) # M1
		block_12 = PSBlock(2*M1.row_block-1, 2*M1.col_block,   new_level) # M3
		block_21 = PSBlock(2*M1.row_block,   2*M1.col_block-1, new_level) # M7

		return block_11, block_12, block_21
	elseif coord == 2
		#  ___________
		# |   __|13 14|
		# |__|   __|24|
		# |  |__|
		# |_____|
		block_14 = PSBlock(2*M3.row_block-1, 2*M3.col_block,   new_level) # M1
		block_24 = PSBlock(2*M3.row_block,   2*M3.col_block,   new_level) # M3
		block_13 = PSBlock(2*M3.row_block-1, 2*M3.col_block-1, new_level) # M7

		return block_14, block_24, block_13
	elseif coord == 3
		#  ___________
		# |   __|__   |
		# |__|   __|__|
		# |31|__|
		# |41 42|
		block_41 = PSBlock(2*M7.row_block,   2*M7.col_block-1, new_level) # M1
		block_31 = PSBlock(2*M7.row_block-1, 2*M7.col_block-1, new_level) # M3
		block_42 = PSBlock(2*M7.row_block,   2*M7.col_block,   new_level) # M7

		return block_41, block_31, block_42
	elseif coord == 4
		#  ___________
		# |   __|__   |
		# |__|22 23|__|
		# |  |32|
		# |_____|
		# Note: The numbering scheme is a little different for this one, but it is correct.
		block_22 = PSBlock(2*M1.row_block,   2*M1.col_block,   new_level) # M1
		block_23 = PSBlock(2*M3.row_block,   2*M3.col_block-1, new_level) # M3
		block_32 = PSBlock(2*M7.row_block-1, 2*M7.col_block,   new_level) # M7

		return block_22, block_23, block_32
	end
end



# Rotate blocks and then place them down into their new location.
function rotate_shrink_and_place_blocks(M1::PSBlock, M3::PSBlock, M7::PSBlock, coord::Int64, block_rotation::Int64)
	# Temporarily set M1, M3, and M7 into the standard rotation of 0 degrees 
	temp_rotation = block_rotation
	while temp_rotation != 0
		rotate_blocks!(M1, M3, M7, 90, temp_rotation)
		temp_rotation = (temp_rotation + 90) % 360
	end

	# Find the new locations of M1, M3, and M7 in the non-rotated frame
	new_M1, new_M3, new_M7 = master_map(M1, M3, M7, coord)

	# Get the new block rotation
	if coord == 2
		block_rotation = (block_rotation + 90) % 360
		temp_rotation = 90
	elseif coord == 3
		block_rotation = (block_rotation + 270) % 360
		temp_rotation = 270
	else
		temp_rotation = 0
	end

	# Rotate the blocks back into the correct orientation
	while temp_rotation != block_rotation

		rotate_blocks!(new_M1, new_M3, new_M7, 90, temp_rotation)
		temp_rotation = (temp_rotation + 90) % 360

		if coord == 4
			# In this case, the rotation simply puts the block where you want it.
			# No shifting needed
			# Do nothing
		else
			# In this case, the blocks must be shifted over upon each rotation
			# to ensure they hold the correct location in their parent block.
			if temp_rotation == 90
				new_M1.col_block += 2
				new_M3.col_block += 2
				new_M7.col_block += 2
			elseif temp_rotation == 180
				new_M1.row_block += 2
				new_M3.row_block += 2
				new_M7.row_block += 2
			elseif temp_rotation == 270
				new_M1.col_block -= 2
				new_M3.col_block -= 2
				new_M7.col_block -= 2
			elseif temp_rotation == 0
				new_M1.row_block -= 2
				new_M3.row_block -= 2
				new_M7.row_block -= 2
			end
		end
	end

	return new_M1, new_M3, new_M7, block_rotation
end


# Buffer coordinates are as follows:
#
#		 ___________
#		|1  __|__ 2 |
#		|__| 4 __|__|
# 		|3 |__|
#		|_____|
#
# Just imagine recursively splitting each of these L-shaped blocks
# while applying the same numbering scheme.
function find_buffer_coordinates(row, col, depth)
	buffer_coordinates = []
	for d in depth-1 : -1 : 0
		if row <= 2^d && col <= 2^d
			# Quadrant (1,1)
			push!(buffer_coordinates, 1)

			# No shift required
		elseif row <= 2^d && col > 2^d
			# Quadrant (1,2)
			push!(buffer_coordinates, 2)

			# Shift col location for next quadrant check
			col = col - 2^d
		elseif row > 2^d && col <= 2^d
			# Quadrant (2,1)
			push!(buffer_coordinates, 3)

			# Shift row location for next quadrant check
			row = row - 2^d
		elseif row > 2^d && col > 2^d
			# Quadrant (2,2)
			push!(buffer_coordinates, 4)

			# Shift both locations for next quadrant check
			row = row - 2^d
			col = col - 2^d
		end
	end

	return buffer_coordinates
end


# Find the quadrant which the block location is in,
# from this the buffer location may be found.
function find_buffer_blocks(block_location::PSBlock)
	row = block_location.row_block
	col = block_location.col_block
	depth = block_location.level

	buffer_coordinates = find_buffer_coordinates(row, col, depth)

	##
	# Now that the buffer coordinate has been obtained
	# map it to a unique L-shaped block in the buffer.
	#
	# M1, M3, and M7 will always be contained within the 
	# set of blocks they were assigned at the previous step.
	#
	# PSBlock rotation is retained over the iterations of the loop
	##
	M1 = PSBlock(1,1,1)
	M3 = PSBlock(1,2,1)
	M7 = PSBlock(2,1,1)
	block_rotation = 0
	for coord in buffer_coordinates
		# Find the new L-shaped region within the parent region
		M1, M3, M7, block_rotation = rotate_shrink_and_place_blocks(M1, M3, M7, coord, block_rotation)
	end

	return M1, M3, M7
end

# C blocks are arranged as
#  _________
# | M4 | M5 |
# |____|____|
# | M2 | M6 |
# |____|____|
function find_C_blocks(block_location::PSBlock)
	row = block_location.row_block
	col = block_location.col_block
	depth = block_location.level

	new_base_row = 2*row - 1
	new_base_col = 2*col - 1
	new_depth = depth + 1

	M2 = PSBlock(new_base_row+1, new_base_col,   new_depth)
	M4 = PSBlock(new_base_row,  new_base_col,    new_depth)
	M5 = PSBlock(new_base_row,  new_base_col+1,  new_depth)
	M6 = PSBlock(new_base_row+1, new_base_col+1, new_depth)

	return M2, M4, M5, M6
end

# This is the single block version of the function for when just a single
# sub-block with no new addition or subtraction operations is involved
# (as is the case for parts of M2, M4, M5, and M6)
function find_new_ops(curr_ops::Vector{Operation}, quadrant::Tuple{Int64, Int64})
	if isempty(curr_ops)
		return [Operation(true, PSBlock(quadrant[1], quadrant[2], 1))]
	else
		# Copy the operations list to prevent overwriting data
		new_ops = deepcopy(curr_ops)
		for op in new_ops
			curr_block = op.arg

			curr_block.row_block = (2*curr_block.row_block - 1) + (quadrant[1] - 1)
			curr_block.col_block = (2*curr_block.col_block - 1) + (quadrant[2] - 1)
			curr_block.level += 1

			op.arg = curr_block

		end

		return new_ops
	end
end



# This is the two block version of the function for when two sub-blocks
# must be added or subtracted
function find_new_ops(curr_ops::Vector{Operation}, quadrant_1::Tuple{Int64, Int64}, new_func::Bool, quadrant_2::Tuple{Int64, Int64})
	
	if isempty(curr_ops)
		return [
			Operation(true,     PSBlock(quadrant_1[1], quadrant_1[2], 1)),
			Operation(new_func, PSBlock(quadrant_2[1], quadrant_2[2], 1)),
		]
	else
		# Create a new operation list as it will double in size
		new_ops::Vector{Operation} = []

		# Copy the operations list to prevent overwriting data
		copy_ops = deepcopy(curr_ops)

		for op in copy_ops
			curr_block = op.arg 

			# Make the second operation first before the curr_block data is overwritten
			other_block = PSBlock(
				(2*curr_block.row_block - 1) + (quadrant_2[1] - 1),
				(2*curr_block.col_block - 1) + (quadrant_2[2] - 1),
				curr_block.level + 1
			)

			# Same sign     -> true +
			# opposite sign -> false -
			this_func = (new_func == op.func)
			other_op = Operation(this_func, other_block)
	
			# Make the first operation
			curr_block.row_block = (2*curr_block.row_block - 1) + (quadrant_1[1] - 1)
			curr_block.col_block = (2*curr_block.col_block - 1) + (quadrant_1[2] - 1)
			curr_block.level += 1

			# Push the first operation onto the list
			push!(new_ops, op)

			# Push the second operation onto the list
			push!(new_ops, other_op)
		end

		return new_ops
	end
end


function create_partial_strassen_equations(
    all_equations        :: Vector{Equation},
	remaining_depth      :: W,
	curr_level           :: W                  = 1,
	block_location       :: PSBlock 		   = PSBlock(1,1,0),
	curr_A_operations    :: Vector{Operation}  = Array{Operation}(undef, 0),
	curr_B_operations    :: Vector{Operation}  = Array{Operation}(undef, 0),
	buffer_bool          :: Bool               = true,
    recombine            :: Bool               = false,
    recombination_blocks :: Vector{PSBlock}    = Array{PSBlock}(undef, 0)
) where W <: Integer

	if remaining_depth > 0
		# In this case, we have a submatrix which must be formed by recombination

		# Get child node locations
		M1_block, M3_block, M7_block = find_buffer_blocks(block_location)
		M2_block, M4_block, M5_block, M6_block = find_C_blocks(block_location)

        #=M4=# 
		create_partial_strassen_equations(
            all_equations,
			remaining_depth-1,
			curr_level+1,
			M4_block,
			find_new_ops(curr_A_operations, (2,2)),               # A22
			find_new_ops(curr_B_operations, (2,1), false, (1,1)), # B21 - B11
			false,
            false
		)
        #=M5=# 
		create_partial_strassen_equations(
            all_equations,
			remaining_depth-1,
			curr_level+1,
			M5_block,
			find_new_ops(curr_A_operations, (1,1), true, (1,2)), # A11 + A12
			find_new_ops(curr_B_operations, (2,2)),              # B22
			false,
            false
		)
        #=M2=# 
		create_partial_strassen_equations(
            all_equations,
			remaining_depth-1,
			curr_level+1,
			M2_block,
			find_new_ops(curr_A_operations, (2,1), true, (2,2)), # A21 + A22
			find_new_ops(curr_B_operations, (1,1)),              # B11
			false,
            false
		)
        #=M6=# 
		create_partial_strassen_equations(
            all_equations,
			remaining_depth-1,
			curr_level+1,
			M6_block,
			find_new_ops(curr_A_operations, (2,1), false, (1,1)), # A21 - A11
			find_new_ops(curr_B_operations, (1,1), true, (1,2)),  # B11 + B12
			false,
            false
		)
        #=M1=# 
		create_partial_strassen_equations(
            all_equations,
			0,
			curr_level+1,
			M1_block,
			find_new_ops(curr_A_operations, (1,1), true, (2,2)), # A11 + A22
			find_new_ops(curr_B_operations, (1,1), true, (2,2)), # B11 + B22
			true,
            false
		)
        #=M3=# 
		create_partial_strassen_equations(
            all_equations,
			0,
			curr_level+1,
			M3_block,
			find_new_ops(curr_A_operations, (1,1)),               # A11
			find_new_ops(curr_B_operations, (1,2), false, (2,2)), # B12 - B22
			true,
            false
		)
        #=M7=# 
		create_partial_strassen_equations(
            all_equations,
			0,
			curr_level+1,
			M7_block,
			find_new_ops(curr_A_operations, (1,2), false, (2,2)), # A12 - A22
			find_new_ops(curr_B_operations, (2,1), true, (2,2)),  # B21 + B22
			true,
            true,
            [M1_block, M2_block, M3_block, M4_block, M5_block, M6_block, M7_block]
		)

	else # remaining_depth == 0
		# In this case, the bottom level has been reached.
		# Create an equation for multiplication of submatrices
        new_equation = Equation(curr_A_operations, curr_B_operations, block_location, buffer_bool, recombine, recombination_blocks)

        push!(all_equations, new_equation)
	end
end


function display_operation_list(ops)
    for op in ops
        print(op.func, " ", op.arg, " ")
    end
    println()
end



function display_partial_strassen_multiplication_equations(depth)
    eqs::Vector{Equation} = []
	create_partial_strassen_equations(eqs, depth)

    for eq in eqs
        println("A:")
        display_operation_list(eq.A_ops)
        println("B:")
        display_operation_list(eq.B_ops)
        println("C:")
        println(eq.C_block)
        println("buffer bool:", eq.buffer_bool)
        println("recombine:", eq.recombine)
        println()
    end
end


