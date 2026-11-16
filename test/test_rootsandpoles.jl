# Every locus of RootsAndPoles.jl's test problems, of the right winding, polished to the
# reference's digits: multiple ones included (`_moment_mean`).
include("rootsandpoles_problems.jl")

function matches(r, loci)
    length(r.points) == length(loci) && all(r.converged) && isempty(r.unresolved) || return false
    return all(loci) do (z, w)
        i = argmin(abs.(r.points .- z))
        r.windings[i] == w && abs(r.points[i] - z) < 2e-9 * max(1, abs(z))
    end
end

@testset "$(p.name)" for p in RP_PROBLEMS
    @test matches(findzeropoles(p.f, p.box, AAA()), p.loci)
    alg, box = p.ap
    box = something(box, p.box)
    r = findzeropoles(p.f, box, alg)
    @test matches(r, filter(l -> inbox(box, l[1]), p.loci)) && r.deficit == 0
end
