module SymMaps
# symmetry reduce map for meshes

using ..TreeMeshes

export SymMap, MappedData

function _find_in(x, arr::AbstractArray; atol=1e-6, rtol=1e-6)
    # return index if in, return 0 otherwise
    for yi in 1:length(arr)
        y = arr[yi]
        if isapprox(x, y, atol=atol, rtol=rtol)
            return yi
        end
    end

    return 0
end

struct SymMap{T,N}
    map::Vector{Int}
    reduced_length::Int
    _vals::Vector{T}
    inv_map::Vector{Vector{Int}}

    function SymMap{T}(tg::TreeGrid, density; atol=1e-6, rtol=1e-6) where {T}
        map = zeros(Int, length(tg))
        reduced_vals = []
        inv_map = []
        for pi in 1:length(tg)
            # println(pi, " ", p)
            p = tg[pi]
            val = density(p)
            # println(val)
            pos = _find_in(val, reduced_vals; atol=atol, rtol=rtol)
            if pos == 0
                push!(reduced_vals, val)
                push!(inv_map, [pi,])
                map[pi] = length(reduced_vals)
            else
                push!(inv_map[pos], pi)
                map[pi] = pos
            end
        end

        return new{T,length(tg)}(map, length(reduced_vals), reduced_vals, inv_map)
    end
end

struct MappedData{T,N} <: AbstractArray{T,N}
    smap::SymMap{T,N}
    data::Vector{T}

    function MappedData(smap::SymMap{T,N}) where {T,N}
        data = zeros(T, smap.reduced_length)
        return new{T,N}(smap, data)
    end
end

Base.length(md::MappedData) = length(md.smap.map)
Base.size(md::MappedData) = (length(md),)
# index and iterator
Base.getindex(md::MappedData, i) = md.data[md.smap.map[i]]
function Base.setindex!(md::MappedData, x, i)
    md.data[md.smap.map[i]] = x
end
Base.firstindex(md::MappedData) = 1
Base.lastindex(md::MappedData) = length(tg)

Base.iterate(md::MappedData) = (md[1], 1)
Base.iterate(md::MappedData, state) = (state >= length(md)) ? nothing : (md[state+1], state + 1)

end