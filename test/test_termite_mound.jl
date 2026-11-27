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

# "splitting" method: first-order Godunov (Lie-Trotter) operator splitting
# (`GodunovSplittingCallback`).
# Both elixirs simulate one full day. The "errors" l2/linf are the difference between the
# final state after one day and the initial condition, which is flat (constant rho, v1 and p,
# only Ti has a spatial profile). Hence the l2/linf of v1 (~0.26-0.27) measure how far the
# velocity has developed away from the constant initial velocity, not a numerical error.
# Reference values: one day simulated with `cfl = 0.05` (time-converged to ~1e-6, see the run
# record at the end of the elixirs). The tests run the elixirs with `cfl = 3.0`, whose O(dt)
# deviation (~6e-5 to 9e-5 in v1) is below `atol`.
# The values were updated on 2026-11-27 after the fix of `x_c = (L - 0.3)/L` (last soil section
# of 0.3 m, before: `L - 0.09`), which shortens the mound region (`I_w`, `T_u`) and changes `k_i`.
@trixi_testset "elixir_termite_splitting.jl" begin
    @test_trixi_include(
        joinpath(EXAMPLES_DIR, "elixir_termite_splitting.jl"),
        l2=[2.85065479e-03, 2.62533517e-01, 2.46944669e-03, 9.18868740e-05],
        linf=[4.81726039e-03, 2.66092963e-01, 2.46944669e-03, 1.76871225e-04],
        atol = 2e-4
    )
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(Trixi.rhs!, semi, sol, 500_000)
end
@trixi_testset "elixir_reference_splitting.jl" begin
    @test_trixi_include(
        joinpath(EXAMPLES_DIR, "elixir_reference_splitting.jl"),
        l2=[2.03338087e-03, 2.64640909e-01, 1.95860867e-03, 2.73873479e-04],
        linf=[3.97670479e-03, 2.67574474e-01, 1.95860867e-03, 4.67148540e-04],
        atol = 2e-4
    )
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(Trixi.rhs!, semi, sol, 500_000)
end

# Manufactured solution with variable cross section A(x) = 1 + 0.3 sin(2πx + π/2), so that all
# A and A_x terms of the equations are active: l2/linf are true errors w.r.t. the exact solution.
# They are dominated by the O(dt) splitting error (first order in dt, and hence in h at fixed cfl).
@trixi_testset "elixir_manufactured_splitting.jl" begin
    @test_trixi_include(
        joinpath(EXAMPLES_DIR, "elixir_manufactured_splitting.jl"),
        l2=[5.70286753e-05, 1.37150107e-02, 4.88006589e-06, 6.28086294e-09],
        linf=[9.17286217e-05, 1.82033360e-02, 4.88006589e-06, 9.43274192e-09]
    )
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(Trixi.rhs!, semi, sol, 500_000)
end
# Levels 3, 4, 5 at fixed cfl: the splitting error O(dt) = O(h) limits the order to one.
# Only rho and v1 are checked: the errors of p and Ti are orders of magnitude smaller and their
# EOCs are not monotone.
@trixi_testset "elixir_manufactured_splitting.jl (convergence)" begin
    using Trixi: convergence_test, calc_mean_convergence
    eocs, _ = convergence_test(
        @__MODULE__,
        joinpath(EXAMPLES_DIR, "elixir_manufactured_splitting.jl"),
        3,
        initial_refinement_level = 3,
    )
    eoc_mean = calc_mean_convergence(eocs)
    @test all(eoc_mean[:l2][1:2] .> 0.9)
    @test all(eoc_mean[:linf][1:2] .> 0.9)
end


# "coupled" method: coupled method of lines (`CoupledSourceTerms`, no operator splitting).
# Same cases as above. Without the O(dt) splitting error the elixirs with `cfl = 3.0` already
# reproduce the time-converged (cfl = 0.05) splitting reference values above up to ~1.4e-6, which
# is the remaining O(dt) splitting error of the reference itself (coupled: cfl 3 vs 0.375 differ
# by < 1e-9; splitting minus coupled halves with cfl: v1 1.3e-5 at cfl 0.375 -> ~1.7e-6 at 0.05).
@trixi_testset "elixir_termite_coupled.jl" begin
    @test_trixi_include(
        joinpath(EXAMPLES_DIR, "elixir_termite_coupled.jl"),
        l2=[2.85065846e-03, 2.62532095e-01, 2.46945124e-03, 9.18867874e-05],
        linf=[4.81726397e-03, 2.66091545e-01, 2.46945124e-03, 1.76870908e-04]
    )
    @test_allocations(Trixi.rhs!, semi, sol, 500_000)
end
@trixi_testset "elixir_reference_coupled.jl" begin
    @test_trixi_include(
        joinpath(EXAMPLES_DIR, "elixir_reference_coupled.jl"),
        l2=[2.03338307e-03, 2.64639869e-01, 1.95861169e-03, 2.73873619e-04],
        linf=[3.97670731e-03, 2.67573436e-01, 1.95861169e-03, 4.67148477e-04]
    )
    @test_allocations(Trixi.rhs!, semi, sol, 500_000)
end
# Manufactured solution (variable cross section): without the splitting error the errors are
# 5-6 orders of magnitude smaller than with the splitting.
@trixi_testset "elixir_manufactured_coupled.jl" begin
    @test_trixi_include(
        joinpath(EXAMPLES_DIR, "elixir_manufactured_coupled.jl"),
        l2=[1.39504267e-08, 7.38850513e-08, 5.22777609e-11, 2.08686057e-11],
        linf=[1.31886604e-07, 1.33669501e-07, 5.22779597e-11, 1.90103711e-10]
    )
    @test_allocations(Trixi.rhs!, semi, sol, 500_000)
end
# Levels 3, 4, 5 at fixed cfl: without splitting the full order p + 1 = 5 of the DG
# discretization is obtained for rho, v1 and Ti. The (tiny, ~5e-11) error of p0 is limited by
# the 4th-order RK time error (dt ∝ h at fixed cfl), hence the separate threshold for p0.
@trixi_testset "elixir_manufactured_coupled.jl (convergence)" begin
    using Trixi: convergence_test, calc_mean_convergence
    eocs, _ = convergence_test(
        @__MODULE__,
        joinpath(EXAMPLES_DIR, "elixir_manufactured_coupled.jl"),
        3,
        initial_refinement_level = 3,
    )
    eoc_mean = calc_mean_convergence(eocs)
    @test all(eoc_mean[:l2][[1, 2, 4]] .> 4)
    @test all(eoc_mean[:linf][[1, 2, 4]] .> 4)
    @test eoc_mean[:l2][3] > 3.5
    @test eoc_mean[:linf][3] > 3.5
end

end # TermiteMound1D

end # module
