module TestExamplePassiveHouse

using Test
using Trixi
using TermiteMoundInducedAirflowTrixi

include("test_trixi.jl")

EXAMPLES_DIR = pkgdir(TermiteMoundInducedAirflowTrixi, "examples", "passive_house")

# Start with a clean environment: remove Trixi.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

@testset "PassiveHouse1D" begin
#! format: noindent

@trixi_testset "passive_house.jl" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "passive_house.jl"))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(Trixi.rhs!, semi, sol, 1000000)
end

end # PasiiveHouse1D

end # module