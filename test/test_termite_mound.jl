module TestExampleTermiteMound

using Test
using Trixi
using TermiteMoundInducedAirflowTrixi

include("test_trixi.jl")

EXAMPLES_DIR = pkgdir(TermiteMoundInducedAirflowTrixi, "test", "examples")

# Start with a clean environment: remove Trixi.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

@testset "TermiteMound1D" begin
#! format: noindent

@trixi_testset "test_termite.jl" begin
    @test_trixi_include(
        joinpath(EXAMPLES_DIR, "test_termite.jl"),
        l2=[
            2.28880456e-03,
            3.64482034e-01,
            1.98390017e-04,
            2.66448807e-05,
            3.99059200e-17,
        ],
        linf=[
            4.65607362e-03,
            3.70896581e-01,
            1.98390017e-04,
            5.50018656e-05,
            2.22044605e-16,
        ],
        atol = 1e-4
    )
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(Trixi.rhs!, semi, sol, 500_000)
end
@trixi_testset "test_reference.jl" begin
    @test_trixi_include(
        joinpath(EXAMPLES_DIR, "test_reference.jl"),
        l2=[
            2.24311704e-03,
            6.98593919e-02,
            1.47349259e-04,
            1.56925914e-04,
            3.99059200e-17,
        ],
        linf=[
            4.86510680e-03,
            7.50852786e-02,
            1.47349259e-04,
            2.39882029e-04,
            2.22044605e-16,
        ],
        atol = 1e-4
    )
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(Trixi.rhs!, semi, sol, 500_000)
end

end # TermiteMound1D

end # module
