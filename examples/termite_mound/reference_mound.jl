using TermiteMoundInducedAirflowTrixi
using Trixi,
    Trixi2Vtk,
    OrdinaryDiffEqLowStorageRK,
    Interpolations,
    QuadGK,
    FastGaussQuadrature,
    Plots
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
    TermiteMoundInitialCondition, source_terms, termite_parameters

###############################################################################
# Semidiscretization

(γ, k_i, k_w, tᵣ, uᵣ, xa, xb, xc, r, h, L, β, η, Fr², ρₕ₀, T_ref, t_ref, T0, v0, Ti_LI) =
    termite_parameters(0.6, 2.0)
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
)

initial_condition = TermiteMoundInitialCondition
volume_flux = flux_ranocha
surface_flux = flux_ranocha

dg = DGSEM(
    polydeg = 4,
    surface_flux = flux_ranocha,
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
    source_terms = source_terms,
    boundary_conditions = boundary_condition_periodic,
);

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 1.0) .* 86400 ./ equations.tᵣ
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()
stepsize_callback = StepsizeCallback(cfl = 0.3)
update_velocity_callback =
    UpdateVelocityCallback(CarpenterKennedy2N54(williamson_condition = false))
amr_controller = ControllerThreeLevel(
    semi,
    IndicatorMax(semi, variable = (u, equations) -> u[2]),
    base_level = 5,
    max_level = 6,
    max_threshold = -0.35,
)
amr_callback = AMRCallback(
    semi,
    amr_controller,
    interval = 1,
    adapt_initial_condition = true,
    adapt_initial_condition_only_refine = true,
)

callbacks = CallbackSet(summary_callback, stepsize_callback, update_velocity_callback)

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
#plot(pd_init)

anim_rho = @animate for i ∈ 1:length(sol.t)
    pd_rho = PlotData1D(sol.u[i], semi)
    plot(pd_rho["rho"])
end
#gif(anim_rho, "density.gif", fps = 5)

anim_v = @animate for i ∈ 1:length(sol.t)
    pd_v = PlotData1D(sol.u[i], semi)
    plot(pd_v["v1"])
end
#gif(anim_v, "velocity.gif", fps = 5)

pd = PlotData1D(sol)
plot(pd)
