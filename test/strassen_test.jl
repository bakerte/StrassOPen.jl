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

import Serialization
file = testpath * "dict" * file_extension
if isfile(file)
  performancevals = Serialization.deserialize(file)
else
  performancevals = Dict()
end

println("#            +-------------------+")
println("#>-----------|     strassen.jl   |-----------<")
println("#            +-------------------+")

fulltest = true

p = 12
m = 2^p

A = rand(m,m)
B = rand(m,m)

checkC = A*B
@time A*B

C = strassen(A,B)
@time strassen(A,B)


import LinearAlgebra


testval = "LinearAlgebra.norm(strassen(A,B)-checkC) < 1E-8"
fulltest &= testfct(testval,"strassen(Array,Array)",performancevals)


tcheck = A'*B
#tCcheck = strassen('T','N',A,B)

testval = "LinearAlgebra.norm(strassen('T','N',A,B)-tcheck) < 1E-8"
fulltest &= testfct(testval,"strassen('T','N',Array,Array)",performancevals)


tcheck = A*B'
#tCcheck = strassen('N','T',A,B)

testval = "LinearAlgebra.norm(strassen('N','T',A,B)-tcheck) < 1E-8"
fulltest &= testfct(testval,"strassen('N','T',Array,Array)",performancevals)



A = rand(ComplexF64,m,m)
B = rand(ComplexF64,m,m)

tcheck = A'*B
#tCcheck = 

testval = "LinearAlgebra.norm(strassen('C','N',A,B)-tcheck) < 1E-8"
fulltest &= testfct(testval,"strassen('C','N',Array,Array)",performancevals)

tcheck = A*B'
#tCcheck = 

testval = "LinearAlgebra.norm(strassen('N','C',A,B)-tcheck) < 1E-8"
fulltest &= testfct(testval,"strassen('N','C',Array,Array)",performancevals)

#A = rand(m,m)
#B = rand(ComplexF64,m,m)


Serialization.serialize(file,performancevals)

