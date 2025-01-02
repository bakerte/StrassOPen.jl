#########################################################################
#
#           Tensor Linear Algebra Package (TENPACK)
#                          v1.0
#
#########################################################################
# Made by Thomas E. Baker and « les qubits volants » (2024)
# See accompanying license with this program
# This code is native to the julia programming language (v1.10.4+)
#

#         +---------------+
#>--------|     mult!     |---------<
#         +---------------+

for (gemm, elty) in
  ((:dgemm_,:Float64),
   (:sgemm_,:Float32),
   (:zgemm_,:ComplexF64),
   (:cgemm_,:ComplexF32))
  @eval begin
        # SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        # *     .. Scalar Arguments ..
        #       DOUBLE PRECISION ALPHA,BETA
        #       INTEGER K,LDA,LDB,LDC,M,N
        #       CHARACTER TRANSA,TRANSB
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
    function matmul!(transA::AbstractChar, transB::AbstractChar,
                  alpha::Union{($elty), Bool},
                  A::AbstractArray{$elty,N},
                  B::AbstractArray{$elty,M},
                  beta::Union{($elty), Bool},
                  C::AbstractArray{$elty,G},
                  m::Integer,ka::Integer,kb::Integer,n::Integer) where {N,M,G}
        lda = max(1,transA == 'N' ? m : ka)
        ldb = max(1,transB == 'N' ? ka : n)
        ldc = max(1,m)
        ccall((@blasfunc($gemm), libblastrampoline), Cvoid,
            (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{$elty}, Ptr{$elty}, Ref{BlasInt},
            Ptr{$elty}, Ref{BlasInt}, Ref{$elty}, Ptr{$elty},
            Ref{BlasInt}, Clong, Clong),
            transA, transB, m, n,
            ka, alpha, A,  lda,
            B, ldb, beta, C,
            ldc, 1, 1)
        C
    end

    function matmul(transA::AbstractChar, transB::AbstractChar, alpha::($elty), A::AbstractArray{$elty,N},B::AbstractArray{$elty,M},m::Integer,ka::Integer,kb::Integer,n::Integer) where {N,M}
        C = Array{($elty),2}(undef,m,n)
        matmul!(transA, transB, alpha, A, B, zero($elty), C,m,ka,kb,n)
        return C
    end
    
    function matmul(transA::AbstractChar, transB::AbstractChar, alpha::($elty), A::AbstractArray{$elty,N},B::AbstractArray{$elty,M},beta::($elty),C::AbstractArray{$elty,G},m::Integer,ka::Integer,kb::Integer,n::Integer) where {N,M,G}
        newC = copy(C)
        matmul!(transA, transB, alpha, A, B, beta, newC,m,ka,kb,n)
        return newC
    end

  end
end

