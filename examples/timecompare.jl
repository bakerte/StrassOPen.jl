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
#using BenchmarkTools

function regmatmul(A,B)
  C = Array{Float64,2}(undef,size(A,1),size(B,2))
  for x = 1:size(A,1)
    for y = 1:size(B,2)
      C[x,y] = 0
      for z = 1:size(A,2)
            C[x,y] += A[x,z]*B[z,y]
      end
    end
  end
  return C
end

m = 1028

A = rand(m,m)
B = rand(m,m)

@time regmatmul(A,B);
@time A*B;
