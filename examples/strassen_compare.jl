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

using StrassOPen
import LinearAlgebra


#using BenchmarkTools


m = 12

b = 2^m

A = rand(b,b)
B = rand(b,b)

C = A*B
@time A*B



#you may choose to increase the number of threads for a fair comparison
#num_threads = 4
#LinearAlgebra.BLAS.set_num_threads(num_threads)



C = A*B
@time A*B

#=
if num_threads > 1
   #However, the number of threads provided to BLAS is often automatically twice the input number on most machines (hyper-threading); for the Strassen it will always be the input number
   LinearAlgebra.BLAS.set_num_threads(2*LinearAlgebra.BLAS.get_num_threads())
end
=#

checkC = StrassOPen.strassen(A,B,n=1)
@time StrassOPen.strassen(A,B,n=1)

nlevel = 2

checkC = StrassOPen.strassen(A,B,n=nlevel)
@time StrassOPen.strassen(A,B,n=nlevel)




LinearAlgebra.norm(C-checkC)
