module TestUnit

using TermiteMoundInducedAirflowTrixi
using Test
using Trixi, OrdinaryDiffEqLowStorageRK
import TermiteMoundInducedAirflowTrixi:
    temp2unscaled, temp2scaled, vel2unscaled, vel2scaled, Get_initial_Ti, T_u

include("test_trixi.jl")

# constant alias (closures must not capture a non-constant module binding: that would make
# the calls dynamic and allocate)
const TMI_CONST = TermiteMoundInducedAirflowTrixi

# Equations with the default scalar parameters (mound r = 0.6 m, h = 2 m) and the termite mound
# model data; `cross_section` can be replaced
function default_equations(; cross_section = CrossSectionConstant())
    (_, _, _, tᵣ, _, xa, xb, xc, _, h, L, _, _, _, _, T_ref, _, _, _, _) =
        termite_parameters(0.6, 2.0)
    return TermiteMoundEquations1D(;
        cross_section,
        height_profile = HeightProfileTermiteMound(; xa, xb, L, h),
        boundary_temperature = BoundaryTemperatureTermiteMound(; xa, xb, xc, tᵣ, T_ref),
    )
end

# Start with a clean environment: remove Trixi.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

# Run various unit (= non-elixir-triggered) tests
@testset "Unit tests" begin
#! format: noindent

@testset "Unit conversions" begin
    @timed_testset "Temperature" begin
        let equations = default_equations()
            t = 1.01
            @test t ≈ temp2scaled(temp2unscaled(t, equations), equations)
            @test t ≈ temp2unscaled(temp2scaled(t, equations), equations)
            z = Get_initial_Ti(0.55, 2.65)
            t = z(0.25)
            @test t ≈ temp2scaled(temp2unscaled(t, equations), equations)
            @test t ≈ temp2unscaled(temp2scaled(t, equations), equations)
        end
    end

    @timed_testset "Velocity" begin
        let equations = default_equations()
            v = 1.01
            @test v ≈ vel2scaled(vel2unscaled(v, equations), equations)
            @test v ≈ vel2unscaled(vel2scaled(v, equations), equations)
        end
    end
end # @testset "Unit conversions"

@testset "Consistency checks" begin
    @timed_testset "TME flux" begin
        let equations = default_equations()
            u_ll = SVector(0.1, 1.0, 0.0, 0.0)
            u_rr = SVector(0.1, 1.0, 0.0, 0.0)
            @test flux_ranocha(u_ll, u_rr, 1, equations) ≈
                  Trixi.flux(u_ll, 1, equations)
            @test flux_hllc(u_ll, u_rr, 1, equations) ≈ Trixi.flux(u_ll, 1, equations)
            u_ll = SVector(0.1, -1.0, 0.0, 0.0)
            u_rr = SVector(0.1, -1.0, 0.0, 0.0)
            @test flux_ranocha(u_ll, u_rr, 1, equations) ≈
                  Trixi.flux(u_ll, 1, equations)
            @test flux_hllc(u_ll, u_rr, 1, equations) ≈ Trixi.flux(u_ll, 1, equations)
        end
    end

    @timed_testset "Boundary temperatures" begin
        let equations = default_equations()
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
    @timed_testset "Ti decoupled from boundary conditions" begin
        # Ti carries no flux: the periodic coupling of the mesh, which only
        # acts through interface fluxes, can never reach it
        let equations = default_equations()
            u_ll = SVector(1.02, -0.8, 1.0, 1.01)
            u_rr = SVector(0.98, -0.7, 1.0, 1.03)
            for f in (
                flux_ranocha(u_ll, u_rr, 1, equations),
                flux_hllc(u_ll, u_rr, 1, equations),
                flux_hllc(u_ll, u_rr, SVector(-1.0), equations),
                Trixi.flux(u_ll, 1, equations),
            )
                @test f[4] == 0
            end
        end
        # Ti only enters the model on (xa, xc] via T_u: perturbing Ti outside of it must
        # leave rho, v1, p0 and the inner Ti bitwise unchanged
        (
            γ,
            k_i,
            k_w,
            tᵣ,
            uᵣ,
            xa,
            xb,
            xc,
            r,
            h,
            L,
            β,
            η,
            Fr²,
            ρₕ₀,
            T_ref,
            t_ref,
            T0,
            v0,
            Ti_LI,
        ) = termite_parameters(1.5, 2.5)
        equations = TermiteMoundEquations1D(;
            γ,
            k_i,
            k_w,
            tᵣ,
            uᵣ,
            xa,
            xb,
            xc,
            r,
            h,
            L,
            β,
            η,
            Fr²,
            ρₕ₀,
            T_ref,
            t_ref,
            T0,
            v0,
            Ti_LI,
            cross_section = CrossSectionConstant(),
            height_profile = HeightProfileTermiteMound(; xa, xb, L, h),
            boundary_temperature = BoundaryTemperatureTermiteMound(;
                xa,
                xb,
                xc,
                tᵣ,
                T_ref,
            ),
        )
        dg = DGSEM(
            polydeg = 4,
            surface_flux = flux_ranocha,
            volume_integral = VolumeIntegralFluxDifferencing(flux_ranocha),
        )
        mesh = TreeMesh(
            (0.0,),
            (1.0,),
            initial_refinement_level = 3,
            n_cells_max = 1000,
            periodicity = true,
        )
        semi = SemidiscretizationHyperbolic(
            mesh,
            equations,
            TermiteMoundInitialCondition,
            dg,
            source_terms = source_terms,
            boundary_conditions = boundary_condition_periodic,
        )
        ode = semidiscretize(semi, (0.0, 0.5))
        x = vec(semi.cache.elements.node_coordinates[1, :, :])
        outside = (x .<= equations.xa) .| (x .>= equations.xc)
        Ti_idx = 4:4:length(ode.u0)
        u0_perturbed = copy(ode.u0)
        u0_perturbed[Ti_idx[outside]] .+= 0.05
        run(u0) = solve(
            ODEProblem(ode.f, u0, ode.tspan, semi),
            CarpenterKennedy2N54(williamson_condition = false);
            dt = 1.0,
            save_everystep = false,
            callback = CallbackSet(
                StepsizeCallback(cfl = 0.4),
                GodunovSplittingCallback(
                    CarpenterKennedy2N54(williamson_condition = false),
                ),
            ),
        ).u[end]
        u_ref = run(ode.u0)
        u_pert = run(u0_perturbed)
        not_outer_Ti = trues(length(u_ref))
        not_outer_Ti[Ti_idx[outside]] .= false
        @test u_ref[not_outer_Ti] == u_pert[not_outer_Ti]
        @test u_ref[Ti_idx[outside]] != u_pert[Ti_idx[outside]]

        # same for the "coupled" method (coupled method of lines)
        semi_coupled = SemidiscretizationHyperbolic(
            mesh,
            equations,
            TermiteMoundInitialCondition,
            dg,
            source_terms = CoupledSourceTerms(source_terms),
            boundary_conditions = boundary_condition_periodic,
        )
        ode_coupled = semidiscretize(semi_coupled, (0.0, 0.5))
        run_coupled(u0) = solve(
            ODEProblem(ode_coupled.f, u0, ode_coupled.tspan, semi_coupled),
            CarpenterKennedy2N54(williamson_condition = false);
            dt = 1.0,
            save_everystep = false,
            callback = CallbackSet(
                AlgebraicVelocityCallback(),
                StepsizeCallback(cfl = 0.4),
            ),
        ).u[end]
        u_coupled_ref = run_coupled(ode_coupled.u0)
        u_coupled_pert = run_coupled(u0_perturbed)
        @test u_coupled_ref[not_outer_Ti] == u_coupled_pert[not_outer_Ti]
        @test u_coupled_ref[Ti_idx[outside]] != u_coupled_pert[Ti_idx[outside]]

        # plots show Ti only inside the mound (xa, xc), all other variables everywhere
        for pd in (
            PlotData1D(u_ref, semi),
            PlotData1D(
                (x, equations) -> TermiteMoundInitialCondition(x, 0.0, equations),
                semi,
            ),
        )
            @test collect(pd.variable_names) == ["rho", "v1", "p", "Ti"]
            inside = equations.xa .< pd.x .< equations.xc
            @test all(isnan, pd.data[.!inside, 4])
            @test !any(isnan, pd.data[inside, 4])
            @test !any(isnan, pd.data[:, 1:3])
        end
    end

    @timed_testset "Cross section" begin
        A = TermiteMoundInducedAirflowTrixi.A
        A_x = TermiteMoundInducedAirflowTrixi.A_x
        let equations = default_equations()
            @test A(0.3, equations) == 1.0
            @test A_x(0.3, equations) == 0.0
        end
        let equations = default_equations(cross_section = CrossSectionSine(0.3))
            xs = range(0.0, 1.0, length = 101)
            @test A(0.0, equations) ≈ 1.3
            @test maximum(A.(xs, Ref(equations))) ≈ A(0.0, equations)
            @test A(0.0, equations) ≈ A(1.0, equations)
            @test A_x(0.0, equations) ≈ A_x(1.0, equations) atol = 1e-12
            dx = 1e-6
            for x in (0.1, 0.37, 0.8)
                @test A_x(x, equations) ≈
                      (A(x + dx, equations) - A(x - dx, equations)) / (2 * dx) rtol =
                    1e-6
            end
        end
        @test_throws ArgumentError CrossSectionSine(1.0)
        @test_throws ArgumentError CrossSectionSine(-1.2)
    end

    @timed_testset "Model data: user-defined vs. predefined" begin
        # the generic wrappers of user-defined functions must reproduce the predefined
        # termite mound model data bitwise
        TMI = TermiteMoundInducedAirflowTrixi
        equations = default_equations(cross_section = CrossSectionSine(0.3))
        cs = equations.cross_section
        hp = equations.height_profile
        bt = equations.boundary_temperature
        cs_user = CrossSection(x -> TMI_CONST.A(x, cs), x -> TMI_CONST.A_x(x, cs))
        hp_user = HeightProfile(x -> TMI_CONST.h_x(x, hp))
        bt_user = BoundaryTemperature((t, x, Ti) -> TMI_CONST.T_u(t, x, Ti, bt))
        xs = range(0.0, 1.0, length = 1001)
        Ti = 1.01
        for t in (0.0, 3.7, 500.0)
            @test TMI.A.(xs, Ref(cs_user)) == TMI.A.(xs, Ref(cs))
            @test TMI.A_x.(xs, Ref(cs_user)) == TMI.A_x.(xs, Ref(cs))
            @test TMI.h_x.(xs, Ref(hp_user)) == TMI.h_x.(xs, Ref(hp))
            @test TMI.T_u.(t, xs, Ti, Ref(bt_user)) == TMI.T_u.(t, xs, Ti, Ref(bt))
            # time-only part evaluated once, as in the solvers
            for b in (bt, bt_user)
                time_data = TMI.T_u_time(t, b)
                @test TMI.T_u_node.(Ref(time_data), xs, Ti, Ref(b)) ==
                      TMI.T_u.(t, xs, Ti, Ref(b))
            end
        end
        @test CrossSectionConstant() isa CrossSectionConstant
        @test TMI.A(0.4, CrossSectionConstant()) == 1.0
        @test TMI.A_x(0.4, CrossSectionConstant()) == 0.0

        # a coupled simulation with user-defined closures: same result, allocation-free RHS
        equations_user = TermiteMoundEquations1D(;
            cross_section = cs_user,
            height_profile = hp_user,
            boundary_temperature = bt_user,
        )
        dg = DGSEM(
            polydeg = 3,
            surface_flux = flux_hllc,
            volume_integral = VolumeIntegralFluxDifferencing(flux_ranocha),
        )
        mesh = TreeMesh(
            (0.0,),
            (1.0,),
            initial_refinement_level = 3,
            n_cells_max = 1000,
            periodicity = true,
        )
        rhs_of(eqs) = begin
            semi = SemidiscretizationHyperbolic(
                mesh,
                eqs,
                TermiteMoundInitialCondition,
                dg,
                source_terms = CoupledSourceTerms(source_terms),
                boundary_conditions = boundary_condition_periodic,
            )
            u_ode = compute_coefficients(0.0, semi)
            du_ode = similar(u_ode)
            Trixi.rhs!(du_ode, u_ode, semi, 2.0)
            @test (@allocated Trixi.rhs!(du_ode, u_ode, semi, 2.0)) == 0
            du_ode
        end
        @test rhs_of(equations_user) == rhs_of(equations)
    end
    @timed_testset "Coupled RHS" begin
        # The coupled RHS must use the same spatial operators as the ODE step of the "splitting"
        # method (`v_rhs`, `p0_rhs`, `dg_cumulative_integral`), evaluated at the current state
        TMI = TermiteMoundInducedAirflowTrixi
        for (cross_section, ode_source_terms) in (
            (CrossSectionConstant(), nothing),
            (CrossSectionSine(0.3), manufactured_source_terms_ode),
        )
            equations = default_equations(; cross_section)
            dg = DGSEM(
                polydeg = 4,
                surface_flux = flux_hllc,
                volume_integral = VolumeIntegralFluxDifferencing(flux_ranocha),
            )
            mesh = TreeMesh(
                (0.0,),
                (1.0,),
                initial_refinement_level = 3,
                n_cells_max = 1000,
                periodicity = true,
            )
            semi = SemidiscretizationHyperbolic(
                mesh,
                equations,
                ManufacturedInitialCondition,
                dg,
                source_terms = CoupledSourceTerms(
                    manufactured_source_terms;
                    ode_source_terms,
                ),
                boundary_conditions = boundary_condition_periodic,
            )
            t = 1.3
            u_ode = compute_coefficients(t, semi)
            du_ode = similar(u_ode)
            Trixi.rhs!(du_ode, u_ode, semi, t)
            @test (@allocated Trixi.rhs!(du_ode, u_ode, semi, t)) == 0
            _, _, _, cache = Trixi.mesh_equations_solver_cache(semi)
            u = Trixi.wrap_array(u_ode, semi)
            du = Trixi.wrap_array(du_ode, semi)
            ops = TMI.dg_operators(dg, cache)
            x = cache.elements.node_coordinates[1, :, :]
            A_nodes = TMI.A.(x, Ref(equations))
            first_node = CartesianIndex(1, ops.order[1])

            # algebraic velocity
            rho = Array(u[1, :, :])
            p0 = Array(u[3, :, :])
            Ti = Array(u[4, :, :])
            T = p0 ./ rho
            Tu = T_u.(t, x, Ti, Ref(equations))
            s_p, s_v = TMI.ode_source_terms(ode_source_terms, x, t, equations)
            Q =
                TMI.add_source(
                    -TMI.I_w.(x, Ref(equations)) .* equations.k_w ./ sqrt.(A_nodes) .*
                    (T .- Tu),
                    s_p,
                ) ./ (equations.γ .* p0[first_node])
            c = sum(Q .* A_nodes .* ops.W) / sum(A_nodes .* ops.W)
            v = A_nodes[first_node] * u[2, first_node]
            v1 = (v .+ TMI.dg_cumulative_integral(A_nodes .* (Q .- c), ops)) ./ A_nodes
            u_alg = copy(u)
            TMI.compute_algebraic_velocity!(
                u_alg,
                t,
                equations,
                ode_source_terms,
                cache.termite,
            )
            @test u_alg[2, :, :] ≈ v1 rtol = 1e-12
            @test u_alg[[1, 3, 4], :, :] == u[[1, 3, 4], :, :]
            # the manufactured initial condition satisfies the constraint
            if cross_section isa CrossSectionSine
                @test maximum(abs, u_alg[2, :, :] .- u[2, :, :]) < 1e-3
            end

            # ODE right-hand sides
            v1_dx = TMI.dg_derivative(v1, ops)
            A_x_nodes = TMI.A_x.(x, Ref(equations))
            dp0 = TMI.p0_rhs(
                p0,
                v1,
                v1_dx,
                T,
                Tu,
                A_nodes,
                A_x_nodes,
                inv.(sqrt.(A_nodes)),
                TMI.I_w.(x, Ref(equations)),
                ops.W,
                equations,
                s_p,
            )
            dv = TMI.v_rhs(rho, v1, v1_dx, x, A_nodes, ops.W, equations, s_v)
            @test all(du[3, :, :] .== du[3, first_node])
            @test du[3, first_node] ≈ dp0 rtol = 1e-10 atol = 1e-14
            @test du[2, :, :] ≈ dv ./ A_nodes rtol = 1e-10
        end
    end
end # @testset "Consistency checks"

end # @testset "Unit tests"
end # module
