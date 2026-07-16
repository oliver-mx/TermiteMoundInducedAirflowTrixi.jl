module TestExamplesUpstream

using Test
using TermiteMoundInducedAirflowTrixi

include("test_trixi.jl")

EXAMPLES_DIR = pkgdir(TermiteMoundInducedAirflowTrixi, "examples")

# Start with a clean environment: remove output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

# Run upstream tests for each mesh and dimension to test compatibility with Trixi.jl
@testset "Upstream tests" begin
#! format: noindent

# Run tests for TreeMesh1D
@trixi_testset "TreeMesh1D: termite_mound.jl" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "termite_mound",
                                 "termite_mound.jl"))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(Trixi.rhs!, semi, sol, 10000)
end

@trixi_testset "TreeMesh1D: passive_house.jl" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "passive_house",
                                 "passive_house.jl"))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(Trixi.rhs!, semi, sol, 1000000)
end

end # Upstream tests

end # module