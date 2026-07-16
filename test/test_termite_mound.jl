module TestExampleTermiteMound

using Test
using Trixi
using TermiteMoundInducedAirflowTrixi

include("test_trixi.jl")

EXAMPLES_DIR = pkgdir(TermiteMoundInducedAirflowTrixi, "examples", "termite_mound")

# Start with a clean environment: remove Trixi.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

@testset "TermiteMound1D" begin
#! format: noindent

@trixi_testset "termite_mound.jl" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "termite_mound.jl"))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(Trixi.rhs!, semi, sol, 1000)
end

end # TermiteMound1D

end # module