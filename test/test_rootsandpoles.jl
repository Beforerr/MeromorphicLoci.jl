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

@testset "a box edge whose phase turns near 2π per vertex spacing is counted" begin
    # along Re z = 2.5 the waveguide's phase turns ≈ 23 rad per unit, past 2π per half edge
    p = only(filter(p -> p.f === wvgd, RP_PROBLEMS))
    @test findzeropoles(p.f, p.box, ArgumentPrinciple(zres = 0.02, minsep = 1.0)).deficit == 0
end
