using TermiteMoundInducedAirflowTrixi
using Trixi, OrdinaryDiffEqLowStorageRK, Plots
using Printf

###############################################################################
# Semidiscretization
#
# Manufactured solution of the TME with variable cross section A(x) = 1 + 0.3 sin(2πx + π/2).
# The numerical solution of the "coupled" method is compared with the exact solution.

(γ, k_i, k_w, tᵣ, uᵣ, xa, xb, xc, r, h, L, β, η, Fr², ρₕ₀, T_ref, t_ref, T0, v0, Ti_LI) =
    termite_parameters(1.5, 2.5)

cross_section = CrossSectionSine(0.3)
height_profile = HeightProfileTermiteMound(; xa, xb, L, h)
boundary_temperature = BoundaryTemperatureTermiteMound(; xa, xb, xc, tᵣ, T_ref)

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
    cross_section,
    height_profile,
    boundary_temperature,
)

initial_condition = ManufacturedInitialCondition
volume_flux = flux_ranocha
surface_flux = flux_hllc

dg = DGSEM(
    polydeg = 4,
    surface_flux = surface_flux,
    volume_integral = VolumeIntegralFluxDifferencing(volume_flux),
)

mesh = TreeMesh(
    (0.0,),
    (1.0,),
    initial_refinement_level = 5,
    n_cells_max = 1000,
    periodicity = true,
)

semi = SemidiscretizationHyperbolic(
    mesh,
    equations,
    initial_condition,
    dg,
    source_terms = CoupledSourceTerms(
        manufactured_source_terms;
        ode_source_terms = manufactured_source_terms_ode,
    ),
    boundary_conditions = boundary_condition_periodic,
);

###############################################################################
# ODE solvers, callbacks etc.

# one time period of the manufactured solution
tspan = (0.0, 10.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()
algebraic_velocity_callback = AlgebraicVelocityCallback()
stepsize_callback = StepsizeCallback(cfl = 0.5)
analysis_callback = AnalysisCallback(semi, interval = 100_000, uEltype = real(dg))

callbacks = CallbackSet(
    summary_callback,
    algebraic_velocity_callback,
    stepsize_callback,
    analysis_callback,
)

###############################################################################
# run the simulation

sol = solve(
    ode,
    CarpenterKennedy2N54(williamson_condition = false);
    dt = 1.0,
    ode_default_options()...,
    saveat = (tspan[end] - tspan[1]) / 100,
    callback = callbacks,
);

###############################################################################
# compare with the exact solution and create plots
#
# Nodal values are plotted directly (not via `PlotData1D`, which hides `Ti` outside the mound;
# the manufactured `Ti` is defined on the whole domain).

x_nodes = vec(semi.cache.elements.node_coordinates[1, :, :])
perm = sortperm(x_nodes)
x = x_nodes[perm]
x_fine = range(0.0, 1.0, length = 401)
var_names = ("rho", "v1", "p", "Ti")

numerical(i, var) = vec(Trixi.wrap_array(sol.u[i], semi)[var, :, :])[perm]
exact(t, var, xs) =
    [ManufacturedInitialCondition(SVector(xi), t, equations)[var] for xi in xs]

# errors over time (max norm at the nodes, absolute and relative to the range of the exact field)
err = [
    maximum(abs.(numerical(i, var) .- exact(sol.t[i], var, x))) for
    i in eachindex(sol.t), var = 1:4
]
println("\nMax nodal error over the simulation (absolute | relative to |exact|):")
for var = 1:4
    scale = maximum(abs.(exact(last(sol.t), var, x_fine)))
    @printf(
        "  %-4s  final %.3e   max over t %.3e   (%.2e relative)\n",
        var_names[var],
        err[end, var],
        maximum(err[:, var]),
        maximum(err[:, var]) / scale
    )
end

# time history of the spatially constant / boundary quantities: p0 and v = A(0) u(x = 0)
A_0 = TermiteMoundInducedAirflowTrixi.A(0.0, equations)
p0_num = [Trixi.wrap_array(sol.u[i], semi)[3, 1, 1] for i in eachindex(sol.t)]
v_num = [A_0 * numerical(i, 2)[1] for i in eachindex(sol.t)]
p0_ex = [exact(t, 3, (0.0,))[1] for t in sol.t]
v_ex = [A_0 * exact(t, 2, (0.0,))[1] for t in sol.t]

# solution: numerical (markers) vs exact (lines); same layout as the error plots
# (rho, v1 / p, Ti), p0 is shown as a time history since it is constant in space
anim = @animate for i in eachindex(sol.t)
    t = sol.t[i]
    plot_rho, plot_v1, plot_Ti = map((1, 2, 4)) do var
        plot(
            x_fine,
            exact(t, var, x_fine),
            label = "exact",
            lw = 2,
            xlabel = "x",
            title = @sprintf("%s, t = %.2f", var_names[var], t)
        )
        scatter!(x, numerical(i, var), label = "numerical", ms = 2)
    end
    plot_p = plot(
        sol.t[1:i],
        p0_ex[1:i],
        label = "exact",
        lw = 2,
        xlims = tspan,
        xlabel = "t",
        title = @sprintf("p, t = %.2f", t)
    )
    scatter!(plot_p, sol.t[1:i], p0_num[1:i], label = "numerical", ms = 2)
    plot(plot_rho, plot_v1, plot_p, plot_Ti, layout = (2, 2), size = (1000, 700))
end
gif(anim, "sol_coupled.gif", fps = 10);

# error history and v(t)
plot(
    plot(
        sol.t,
        max.(err, eps()),
        label = permutedims(collect(var_names)),
        yscale = :log10,
        xlabel = "t",
        title = "max nodal error",
    ),
    plot(sol.t, [v_ex v_num], label = ["v exact" "v numerical"], xlabel = "t", title = "v"),
    layout = (1, 2),
    size = (1000, 400),
)
savefig("err_coupled.png");

# pointwise error numerical - exact
anim_err = @animate for i in eachindex(sol.t)
    t = sol.t[i]
    plots = map(1:4) do var
        plot(
            x,
            numerical(i, var) .- exact(t, var, x),
            label = "",
            xlabel = "x",
            title = @sprintf("error %s, t = %.2f", var_names[var], t)
        )
    end
    plot(plots..., layout = (2, 2), size = (1000, 700))
end
gif(anim_err, "err2_coupled.gif", fps = 10);
