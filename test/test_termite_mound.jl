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
            2.13984688e-03,
            3.41899299e-01,
            7.03547907e-04,
            1.47956067e-05,
            3.99059200e-17,
        ],
        linf=[
            4.28105399e-03,
            3.47697851e-01,
            7.03547907e-04,
            3.70136309e-05,
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
            1.95925763e-03,
            1.76925676e-02,
            7.17164354e-04,
            5.31868089e-05,
            3.99059200e-17,
        ],
        linf=[
            4.28060276e-03,
            2.19855918e-02,
            7.17164354e-04,
            8.80297913e-05,
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
