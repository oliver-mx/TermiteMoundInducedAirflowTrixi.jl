using TermiteMoundInducedAirflowTrixi
using Trixi, OrdinaryDiffEqLowStorageRK, Plots

###############################################################################
# Semidiscretization
#
# Full 24h simulation of the reference termite mound using the "coupled" method.

(γ, k_i, k_w, tᵣ, uᵣ, xa, xb, xc, r, h, L, β, η, Fr², ρₕ₀, T_ref, t_ref, T0, v0, Ti_LI) =
    termite_parameters(0.6, 2.0)

cross_section = CrossSectionConstant()
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

initial_condition = TermiteMoundInitialCondition
volume_flux = flux_ranocha
surface_flux = flux_hllc

dg = DGSEM(
    polydeg = 4,
    surface_flux = flux_hllc,
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
    source_terms = CoupledSourceTerms(source_terms),
    boundary_conditions = boundary_condition_periodic,
);

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 1.0) .* 86400 ./ equations.tᵣ
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()
algebraic_velocity_callback = AlgebraicVelocityCallback()
stepsize_callback = StepsizeCallback(cfl = 3.0)

callbacks = CallbackSet(summary_callback, algebraic_velocity_callback, stepsize_callback)

###############################################################################
# run the simulation

sol = solve(
    ode,
    CarpenterKennedy2N54(williamson_condition = false);
    dt = 1.0,
    ode_default_options()...,
    saveat = (tspan[end]-tspan[1])/24/2,
    callback = callbacks,
);

###############################################################################
# create plots

pd_init = PlotData1D((x, equations) -> initial_condition(x, last(tspan), equations), semi)
plot(pd_init);
savefig("initial_coupled.png");

# density and velocity together, time in hours
anim = @animate for i ∈ 1:length(sol.t)
    pd_i = PlotData1D(sol.u[i], semi)
    plot(
        plot(pd_i["rho"]),
        plot(pd_i["v1"]),
        layout = (1, 2),
        size = (1400, 500),
        plot_title = string(
            "t = ",
            round(sol.t[i] * equations.tᵣ / 3600, digits = 1),
            " h",
        ),
    );
end
gif(anim, "evolution_coupled.gif", fps = 5);

pd = PlotData1D(sol);
plot(pd);
savefig("final_coupled.png");
