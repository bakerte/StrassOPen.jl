using BenchmarkTools

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

@btime regmatmul(A,B);
@btime A*B;
