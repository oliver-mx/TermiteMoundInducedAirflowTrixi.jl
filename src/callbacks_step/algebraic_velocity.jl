# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

@doc raw"""
    AlgebraicVelocityCallback()

Callback of the "coupled" method (coupled method of lines, see [`CoupledSourceTerms`](@ref)).
At initialization and after every time step, it writes the algebraic velocity
```math
u = \frac{v + I_0}{A}, \qquad I_0(x) = \int_{x_{min}}^x A \, (Q - c) \, dy,
```
determined by the current state ``(\rho, v, p_0, T_i)``, into the `v1` field of the solution.
Here ``Q = -I_w k_w (T - T_u) / (\sqrt{A} \gamma p_0)`` (plus an optional source ``s_p``),
``c = \int A Q \, dx / \int A \, dx`` and ``T = p_0 / \rho``. The auxiliary velocity
``v = A(x_{min}) \, v_1(x_{min})`` is stored at the first node, where ``I_0 = 0``.

The callback does not change the solution of the coupled system: the right-hand side reads
only ``v`` (first node), ``\rho``, ``p_0`` and ``T_i``, and recomputes ``u`` from the
velocity constraint in every Runge-Kutta stage. The callback is needed because the
`StepsizeCallback` (CFL condition), the `AnalysisCallback` and the output read the `v1` field
directly; without it they would see outdated velocities. At initialization it replaces the
initial `v1` (e.g. a constant ``v_0``) by the velocity that satisfies the constraint.

### Usage
Use it together with `source_terms = CoupledSourceTerms(...)` and place it **before** the
`StepsizeCallback` and the `AnalysisCallback`:
```julia
algebraic_velocity_callback = AlgebraicVelocityCallback()
callbacks = CallbackSet(summary_callback, algebraic_velocity_callback, stepsize_callback,
                        analysis_callback)
```

!!! note
    In contrast to the [`GodunovSplittingCallback`](@ref) of the "splitting" method, this
    callback does not advance ``v`` and ``p_0`` in time; this is done by the time integrator.
    An error is thrown if the semidiscretization was not created with
    [`CoupledSourceTerms`](@ref).
"""
function AlgebraicVelocityCallback()
    algebraic_velocity = AlgebraicVelocity()
    return DiscreteCallback(
        algebraic_velocity,
        algebraic_velocity;
        initialize = initialize_algebraic_velocity!,
        save_positions = (false, false),
    )
end

struct AlgebraicVelocity end

function Base.show(io::IO, ::DiscreteCallback{<:Any,<:AlgebraicVelocity})
    print(io, "AlgebraicVelocityCallback()")
end

function Base.show(
    io::IO,
    ::MIME"text/plain",
    cb::DiscreteCallback{<:Any,<:AlgebraicVelocity},
)
    if get(io, :compact, false)
        show(io, cb)
    else
        Trixi.summary_box(io, "AlgebraicVelocityCallback")
    end
end

function initialize_algebraic_velocity!(
    cb::DiscreteCallback{Condition,Affect!},
    u,
    t,
    integrator,
) where {Condition,Affect!<:AlgebraicVelocity}
    cb.affect!(integrator)
    return nothing
end

# condition: called after every time step
(::AlgebraicVelocity)(u, t, integrator) = true

# affect!
function (::AlgebraicVelocity)(integrator)
    semi = integrator.p
    mesh, equations, dg, cache = Trixi.mesh_equations_solver_cache(semi)
    source_terms = semi.source_terms
    if !(source_terms isa CoupledSourceTerms)
        error(
            "AlgebraicVelocityCallback requires `source_terms = CoupledSourceTerms(...)`",
        )
    end
    u = Trixi.wrap_array(integrator.u, semi)
    compute_algebraic_velocity!(
        u,
        integrator.t,
        equations,
        source_terms.ode_source_terms,
        cache.termite,
    )
    return nothing
end
end # @muladd
