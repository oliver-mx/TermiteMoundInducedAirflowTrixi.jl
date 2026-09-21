module TestUnit

using TermiteMoundInducedAirflowTrixi
using Test
using Trixi, OrdinaryDiffEqLowStorageRK, Interpolations, FastGaussQuadrature
using Trixi: AbstractEquations, @muladd
import Interpolations: Line
import Trixi:
    flux_ranocha,
    ln_mean,
    inv_ln_mean,
    flux,
    varnames,
    cons2cons,
    cons2prim,
    prim2cons,
    cons2entropy,
    max_abs_speeds
import TermiteMoundInducedAirflowTrixi:
    temp2unscaled, temp2scaled, vel2unscaled, vel2scaled, Get_initial_Ti, T_u

include("test_trixi.jl")

# Start with a clean environment: remove TrixiShallowWater.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

# Run various unit (= non-elixir-triggered) tests
@testset "Unit tests" begin
#! format: noindent

@testset "Unit conversions" begin
    @timed_testset "Temperature" begin
        let equations = TermiteMoundEquations1D()
            t = 1.01
            @test t ≈ temp2scaled(temp2unscaled(t, equations), equations)
            @test t ≈ temp2unscaled(temp2scaled(t, equations), equations)
            z = Get_initial_Ti(0.55, 2.65)
            t = z(0.25)
            @test t ≈ temp2scaled(temp2unscaled(t, equations), equations)
            @test t ≈ temp2unscaled(temp2scaled(t, equations), equations)
        end
    end

    @timed_testset "Veloctiy" begin
        let equations = TermiteMoundEquations1D()
            v = 1.01
            @test v ≈ vel2scaled(vel2unscaled(v, equations), equations)
            @test v ≈ vel2unscaled(vel2scaled(v, equations), equations)
        end
    end
end # @testset "Unit conversions"

@testset "Consistency checks" begin
    @timed_testset "TME flux" begin
        let equations = TermiteMoundEquations1D()
            u_ll = SVector(0.1, 1.0, 0.0, 0.0, 0.0)
            u_rr = SVector(0.1, 1.0, 0.0, 0.0, 0.0)
            @test flux_ranocha(u_ll, u_rr, 1, equations) ≈
                  Trixi.flux(u_ll, 1, equations)
            u_ll = SVector(0.1, -1.0, 0.0, 0.0, 0.0)
            u_rr = SVector(0.1, -1.0, 0.0, 0.0, 0.0)
            @test flux_ranocha(u_ll, u_rr, 1, equations) ≈
                  Trixi.flux(u_ll, 1, equations)
        end
    end

    @timed_testset "Boundary temperatures" begin
        let equations = TermiteMoundEquations1D()
            x_soil = 0.5 * (0.0 + equations.xa)
            x_flute = 0.5 * (equations.xa + equations.xb)
            x_chimney = 0.5 * (equations.xb + equations.xc)
            x_soil_2 = 0.5 * (equations.xc + 1.0)
            t1 = 1.5 .* 86400 ./ equations.tᵣ
            t2 = 2.5 .* 86400 ./ equations.tᵣ
            z = Get_initial_Ti(0.55, 2.65)
            Ti = z(0.25)
            @test T_u(t1, x_soil, Ti, equations) ≈ T_u(t2, x_soil_2, Ti, equations)
            @test T_u(t1, x_flute, Ti, equations) ≈ T_u(t2, x_flute, Ti, equations)
            @test T_u(t1, x_chimney, Ti, equations) ≈ T_u(t2, x_chimney, Ti, equations)
        end
    end
end # @testset "Consistency checks"

end # @testset "Unit tests"
end # module
