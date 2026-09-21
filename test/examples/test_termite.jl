using TermiteMoundInducedAirflowTrixi
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

###############################################################################
# Semidiscretization

(γ, k_i, k_w, tᵣ, uᵣ, xa, xb, xc, r, h, L, β, η, Fr², ρₕ₀, T_ref, t_ref, T0, v0, Ti_LI) =
    termite_parameters(1.5, 2.5)
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

tspan = (0.0, 30.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()
stepsize_callback = StepsizeCallback(cfl = 0.2)
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
# Simulation running 'TermiteMoundEquations1D' with DGSEM(polydeg=4)
#────────────────────────────────────────────────────────────────────────────────────────────────────
#timesteps:              79973                run time:       2.24382625e+01 s
# Δt:             1.22227279e-04                └── GC time:    2.52116590e+00 s (11.236%)
# sim. time:      3.00000000e+01 (100.000%)     time/DOF/rhs!:  2.27753431e-08 s
#                                               PID:            3.50668024e-07 s
#DOFs per field:           160                alloc'd memory:         83.859 MiB
#elements:                  32                device memory:           0.000 MiB
#
# Variable:       rho              v1               p0               Ti               x_var
# L2 error:       2.13984688e-03   3.41899299e-01   7.03547907e-04   1.47956067e-05   3.99059200e-17
# Linf error:     4.28105399e-03   3.47697851e-01   7.03547907e-04   3.70136309e-05   2.22044605e-16
# ∑∂S/∂U ⋅ Uₜ :  -4.69099982e-07
#────────────────────────────────────────────────────────────────────────────────────────────────────
