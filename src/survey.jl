"""
    Sample{M}

A point `(z, p...)` of a locus.
"""
struct Sample{M}
    z::ComplexF64
    p::NTuple{M,Float64}
end

"""
    Branch{M} <: AbstractVector{Sample{M}}

The loci over `M` parameters that following connects, as their [`Sample`](@ref)s sorted
by parameters. `b.winding` is their multiplicity, negative for poles, and 0 for a seed
whose kind could not be read (a branch of that one sample); `b.degree` how many of them
the branch holds at one parameter point: 1 for a single-valued `z*(p)`.
`b.uncertified` lists the stretches `(p1, p2)` over which the branch could not be
certified as one locus ([What is certified](@ref)).
"""
struct Branch{M} <: AbstractVector{Sample{M}}
    samples::Vector{Sample{M}}
    winding::Int
    uncertified::Vector{NTuple{2,NTuple{M,Float64}}}
    degree::Int
end

Base.size(b::Branch) = size(b.samples)
Base.getindex(b::Branch, i::Int) = b.samples[i]

_kind(b::Branch) = b.winding > 0 ? "Zero" : b.winding < 0 ? "Pole" : "Indeterminate"

function Base.summary(io::IO, b::Branch{M}) where {M}
    print(io, "Branch{", M, "}(", b.winding, " ", _kind(b), ", ", length(b), " samples")
    b.degree == 1 || print(io, ", degree ", b.degree)
    isempty(b.uncertified) || print(io, ", uncertified on ", _count(length(b.uncertified), "stretch", "es"))
    print(io, ")")
end
Base.show(io::IO, b::Branch) = summary(io, b)

"""
    Survey{M} <: AbstractVector{Branch{M}}

The [`Branch`](@ref)es [`findloci`](@ref) found over `M` parameters. `s.deficit` lists the
knots where the zeros minus poles in the box, ``N(p)``, differ from the sum of the windings
``w_b`` of the branches there ([The count](@ref)): `p => n` with ``n = N(p) - \\sum_b w_b``, or
`p => nothing` where ``N(p)`` cannot be read.
"""
struct Survey{M} <: AbstractVector{Branch{M}}
    branches::Vector{Branch{M}}
    nevals::Int
    deficit::Vector{Pair{NTuple{M,Float64},Union{Nothing,Int}}}
end

nevals(s::Survey) = s.nevals
Base.size(s::Survey) = size(s.branches)
Base.getindex(s::Survey, i::Int) = s.branches[i]
Base.IndexStyle(::Type{<:Survey}) = IndexLinear()

function Base.summary(io::IO, s::Survey{M}) where {M}
    nz = count(b -> b.winding > 0, s.branches)
    np = count(b -> b.winding < 0, s.branches)
    print(io, "Survey{", M, "} with ", _count(length(s.branches), "branch", "es"), " (",
        nz, " zero, ", np, " pole")
    (i = count(b -> b.winding == 0, s.branches)) > 0 && print(io, ", ", i, " indeterminate")
    print(io, "), ", s.nevals, " evals")
    (k = count(d -> !isnothing(last(d)), s.deficit)) > 0 && print(io, ", count off at ", _count(k, "knot"))
    (k = count(d -> isnothing(last(d)), s.deficit)) > 0 && print(io, ", unchecked at ", _count(k, "knot"))
    (k = count(b -> !isempty(b.uncertified), s.branches)) > 0 && print(io, ", ", _count(k, "branch", "es"), " uncertified")
    return nothing
end

_count(n, word, plural="s") = string(n, " ", word, n == 1 ? "" : plural)
