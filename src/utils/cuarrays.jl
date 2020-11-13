# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

import .CuArrays

# Aovid `mat + a * I` with CuArrays which involes scalar operations and is slow
safe_diagm(mat::CuArrays.CuArray, a::T) where {T} = CuArrays.CuArray{T}(diagm(0 => fill(a, size(mat, 1))))
