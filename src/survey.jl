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

One locus `z*(p)` followed over `M` parameters, as [`Sample`](@ref)s sorted by parameter.

- `b.winding`: the multiplicity; 0 for a lone sample whose kind could not be read.
- `b.degree`: how many loci the branch joins at one parameter point
  ([Several parameters](@ref)); 1 for a single-valued `z*(p)`.
- `b.uncertified`: stretches `(p1, p2)` where another locus came within `zfloor`, so the
  branch may have switched to it, and `(p, p)` where following stopped at `p`
  ([What is certified](@ref)).
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

The [`Branch`](@ref)es over `M` parameters. `s.deficit` lists
`p => n` at the knots where ``n = N(p) - \\sum_b w_b``, the zeros minus poles in the box
less the windings of the branches there, is nonzero ([The count](@ref)), and `p => nothing`
where ``N(p)`` cannot be read. ``n > 0`` means the seed missed loci at `p`. `s.nevals`
counts the evaluations of `f(p...)(z)`.
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
