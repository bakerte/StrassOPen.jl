#########################################################################
#
#           Strassen Open Implementation (StrassOPen)
#                          v1.0
#
#########################################################################
# Made by Thomas E. Baker, Kiana Gallagher, and Aaron Dayton and « les qubits volants » (2025)
# See accompanying license with this program
# This code is native to the julia programming language (v1.11.0+)
#

"""
StrassOPen  (version 1.0)\n
(made for julia v1.11.0+ (October 14, 2024), see included license)

Code: https://github.com/bakerte/StrassOPen.jl

Documentation: A. Dayton, K. Gallagher, and T.E. Baker, "forthcoming"\n

Funding for this program is graciously provided by:
   + Institut quantique (Université de Sherbrooke)
   + Département de physique, Université de Sherbrooke
   + Canada First Research Excellence Fund (CFREF)
   + Institut Transdisciplinaire d'Information Quantique (INTRIQ)
   + US-UK Fulbright Commission (Bureau of Education and Cultural Affairs from the United States Department of State)
   + Department of Physics, University of York
   + Canada Research Chair in Quantum Computing for Modelling of Molecules and Materials
   + Department of Physics & Astronomy, University of Victoria
   + Department of Chemistry, University of Victoria
   + Faculty of Science, University of Victoria
   + National Science and Engineering Research Council (NSERC)

# Warning:

We recommend not defining `using LinearAlgebra` to avoid conflicts.  Instead, define
```julia
import LinearAlgebra
```
and define functions as `LinearAlgebra.svd` to use functions from that package.

"""
module StrassOPen

   const libdir = @__DIR__
   const libpath = libdir*"/lib/"

   files = ["types.jl","imports.jl","exports.jl"]
   for w = 1:length(files)
      include(libpath*files[w])
   end

   println("initialization properly loaded...more functionality coming soon")

#=
   files = ["libalg.jl","blas.jl","matmul.jl"]
   const mathpath = libpath*"libalg/"
   for w = 1:length(files)
      include(mathpath*files[w])
   end

   files = ["strassOPen.jl","Strassen.jl"]
   const strassenpath = libpath*"strassen/"
   for w = 1:length(files)
      include(strassenpath*files[w])
   end

   const testpath = libdir*"/test/"
   include(testpath*"alltest.jl")
=#

end

#using .StrassOPen
