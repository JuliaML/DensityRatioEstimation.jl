import .CuArrays

# Aovid `Kdede + λ * I` with CuArrays which involes scalar operations and is slow
safe_diagm(Kdede::T, λ) where {T<:CuArrays.CuArray} = T(diagm(0 => fill(λ, size(Kdede, 1)))) 
