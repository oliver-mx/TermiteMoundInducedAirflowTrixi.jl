using TermiteMoundInducedAirflowTrixi
using Trixi, OrdinaryDiffEqLowStorageRK, Interpolations, FastGaussQuadrature
using Trixi: AbstractEquations, @muladd
import Interpolations: Line
import Trixi:
    flux_ranocha,
    flux_hllc,
    ln_mean,
    inv_ln_mean,
    flux,
    varnames,
    cons2prim,
    prim2cons,
    cons2entropy,
    cons2cons,
    max_abs_speeds,
    max_abs_speed_naive

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
surface_flux = flux_hllc

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

tspan = (0.0, 75.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()
stepsize_callback = StepsizeCallback(cfl = 0.5)
update_velocity_callback =
    UpdateVelocityCallback(CarpenterKennedy2N54(williamson_condition = false))
analysis_callback = AnalysisCallback(semi, interval = 100_000, uEltype = real(dg))

callbacks = CallbackSet(
    summary_callback,
    stepsize_callback,
    update_velocity_callback,
    analysis_callback,
)

###############################################################################
# run the simulation

sol = solve(
    ode,
    CarpenterKennedy2N54(williamson_condition = false);
    dt = 1.0,
    ode_default_options()...,
    callback = callbacks,
);

#────────────────────────────────────────────────────────────────────────────────────────────────────
# Simulation running 'TermiteMoundEquations1D' with DGSEM(polydeg=4) and cfl=0.25
#────────────────────────────────────────────────────────────────────────────────────────────────────
#timesteps:             106415                run time:       3.07041996e+01 s
# Δt:             3.95991094e-04                └── GC time:    3.44538290e+00 s (11.221%)
# sim. time:      7.50000000e+01 (100.000%)     time/DOF/rhs!:  2.27529539e-08 s
#                                               PID:            3.55300156e-07 s
#DOFs per field:           160                alloc'd memory:        100.502 MiB
#elements:                  32                device memory:           0.000 MiB
#
# Variable:       rho              v1               p                Ti               x_var
# L2 error:       2.24311704e-03   6.98593919e-02   1.47349259e-04   1.56925914e-04   3.99059200e-17
# Linf error:     4.86510680e-03   7.50852786e-02   1.47349259e-04   2.39882029e-04   2.22044605e-16
# ∑∂S/∂U ⋅ Uₜ :  -2.36750442e-06
#────────────────────────────────────────────────────────────────────────────────────────────────────
