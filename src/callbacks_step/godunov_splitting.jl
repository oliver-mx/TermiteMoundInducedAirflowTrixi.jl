# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent
@doc raw"""
    GodunovSplittingCallback(solver; source_terms = nothing)

Callback of the "splitting" method: a first-order Godunov (Lie-Trotter) operator splitting of
the TME. Every time step ``[t_{prev}, t_{now}]`` with ``\Delta t = t_{now} - t_{prev}``
consists of two sub-steps:
1. **PDE step** (Trixi): ``\rho`` (transport PDE) and ``T_i`` (ODE at every node, a source
   term in Trixi) are advanced from ``t_{prev}`` to ``t_{now}`` with ``u`` and ``p_0`` frozen (the `v1` and `p` components of the Trixi right-hand side
   are zero).
2. **ODE step** (this callback, called after every time step):
   - ``v`` is advanced by a predictor-corrector (Heun) step: explicit Euler predictor
     ``v^* = v_{prev} + \Delta t \, v'_{prev}``, then the trapezoidal rule
     ```math
     v = v_{prev} + \frac{\Delta t}{2} \left( v'_{prev} + v'(\rho_{now}, u^*) \right),
     ```
     where ``u^* = u_{prev} + (v^* - v_{prev})/A`` is the shifted velocity profile,
   - the velocity is computed from the velocity constraint,
     ```math
     u = \frac{v + I_0}{A}, \qquad I_0(x) = \int_{x_{min}}^x A \, (Q - c) \, dy,
     ```
     with ``Q = -I_w k_w (T - T_u) / (\sqrt{A} \gamma p_0)`` (plus an optional source
     ``s_p``) and ``c = \int A Q \, dx / \int A \, dx``, evaluated at the state at
     ``t_{prev}``,
   - ``p_0`` is advanced by the trapezoidal rule with the rates [`p0_rhs`](@ref) at
     ``t_{prev}`` and at ``t_{now}`` (new ``\rho``, ``T_i``, ``u``),
   - ``u`` is written into the `v1` field and ``p_0`` into the `p` field at all nodes.

The auxiliary velocity ``v = A(x_{min}) \, v_1(x_{min})`` is stored at the first node, where
``I_0 = 0``. Integrals and derivatives use the DG operators of the solver (LGL quadrature,
elementwise derivative and cumulative integral, see [`dg_operators`](@ref)).
The splitting error is of first order in the time step.

### Parameters
- `solver`: ODE solver used for the integration of the (linearly interpolated) rates of
  ``v`` and ``p_0``. Since these rates are linear in time, any Runge-Kutta method gives the
  trapezoidal rule. We recommend `CarpenterKennedy2N54` from `OrdinaryDiffEqLowStorageRK`.
- `source_terms`: optional additional sources `(x, t, equations) -> SVector(s_p, s_v)` of the
  ``p_0`` and ``v`` equations (e.g. [`manufactured_source_terms_ode`](@ref)); `nothing`
  (default) for no sources.

### Usage
Use it with the plain Trixi `source_terms` of ``\rho`` and ``T_i`` (e.g. [`source_terms`](@ref)):
```julia
godunov_splitting_callback =
    GodunovSplittingCallback(CarpenterKennedy2N54(williamson_condition = false))
callbacks = CallbackSet(summary_callback, stepsize_callback, godunov_splitting_callback)
```

!!! note
    If the callback is not used, the quantities ``u``, ``v`` and ``p_0`` remain constant in
    time. For the coupled method of lines without splitting error, see
    [`CoupledSourceTerms`](@ref) and [`AlgebraicVelocityCallback`](@ref).
"""
mutable struct GodunovSplittingCallback{ODE_solver,SourceTerms}
    solver::ODE_solver
    source_terms::SourceTerms
end

function Base.show(io::IO, cb::DiscreteCallback{<:Any,<:GodunovSplittingCallback})
    @nospecialize cb # reduce precompilation time

    godunov_splitting_callback = cb.affect!
    print(
        io,
        "GodunovSplittingCallback(solver = ",
        godunov_splitting_callback.solver,
        ")",
    )
end

function Base.show(
    io::IO,
    ::MIME"text/plain",
    cb::DiscreteCallback{<:Any,<:GodunovSplittingCallback},
)
    @nospecialize cb # reduce precompilation time

    if get(io, :compact, false)
        show(io, cb)
    else
        godunov_splitting_callback = cb.affect!

        setup = ["Solver" => godunov_splitting_callback.solver]
        Trixi.summary_box(io, "GodunovSplittingCallback", setup)
    end
end

function GodunovSplittingCallback(solver; source_terms = nothing)
    godunov_splitting_callback =
        GodunovSplittingCallback{typeof(solver),typeof(source_terms)}(
            solver,
            source_terms,
        )

    DiscreteCallback(
        godunov_splitting_callback,
        godunov_splitting_callback;
        save_positions = (false, false),
    )
end

# condition: called after every time step
(::GodunovSplittingCallback)(u, t, integrator) = true

@inline function (godunov_splitting_callback::GodunovSplittingCallback)(integrator)
    semi = integrator.p
    mesh, equations, dg, cache = Trixi.mesh_equations_solver_cache(semi)
    ops = dg_operators(dg, cache)
    @unpack W, order = ops
    # nodal fields as (node, element) matrices
    x = cache.elements.node_coordinates[1, :, :]
    u_prev = Trixi.wrap_array(integrator.uprev, semi)
    u_now = Trixi.wrap_array(integrator.u, semi)
    #-----------------------
    A_nodes = A.(x, Ref(equations))
    A_x_nodes = A_x.(x, Ref(equations))
    A_sqrt_A_inv = inv.(sqrt.(A_nodes))
    I_w_nodes = I_w.(x, Ref(equations))
    # first node of the leftmost element, x = x_min: the cumulative integral I_0 vanishes
    # here by construction, so A * v1 at this node is exactly v
    first_node = CartesianIndex(1, order[1])
    #-----------------------
    t_prev = integrator.tprev
    rho_prev = Array(u_prev[1, :, :])
    v1_prev = Array(u_prev[2, :, :])
    p0_prev = Array(u_prev[3, :, :])
    Ti_prev = Array(u_prev[4, :, :])
    v1_dx_prev = dg_derivative(v1_prev, ops)
    v_prev = A_nodes[first_node] * v1_prev[first_node]
    T_prev = p0_prev ./ rho_prev
    Tu_prev = T_u.(t_prev, x, Ti_prev, Ref(equations)) .* equations.t_ref
    s_p_prev, s_v_prev =
        ode_source_terms(godunov_splitting_callback.source_terms, x, t_prev, equations)
    #-----------------------
    p0_dt_prev = p0_rhs(
        p0_prev,
        v1_prev,
        v1_dx_prev,
        T_prev,
        Tu_prev,
        A_nodes,
        A_x_nodes,
        A_sqrt_A_inv,
        I_w_nodes,
        W,
        equations,
        s_p_prev,
    )
    v_dt_prev = v_rhs(rho_prev, v1_prev, v1_dx_prev, x, A_nodes, W, equations, s_v_prev)
    #-----------------------
    t_now = integrator.t
    #-----------------------
    rho = Array(u_now[1, :, :])
    if any(isnan, rho)
        error("NaN detected in rho!")
    end
    p0 = Array(u_now[3, :, :])
    Ti = Array(u_now[4, :, :])
    T = p0 ./ rho
    Tu = T_u.(t_now, x, Ti, Ref(equations)) .* equations.t_ref
    s_p, s_v =
        ode_source_terms(godunov_splitting_callback.source_terms, x, t_now, equations)
    #-----------------------
    # predictor for v: shift the velocity profile, u = (v + I_0) / A
    v_exp = v_prev + integrator.dt * v_dt_prev
    v1_exp = v1_prev .+ (v_exp - v_prev) ./ A_nodes
    v1_dx_exp = dg_derivative(v1_exp, ops)
    v_dt = v_rhs(rho, v1_exp, v1_dx_exp, x, A_nodes, W, equations, s_v)
    LI = linear_interpolation(
        [t_prev, t_now],
        [v_dt_prev, v_dt],
        extrapolation_bc = Line(),
    )
    prob = ODEProblem((u, p, t) -> LI(t), v_prev, (t_prev, t_now))
    sol = solve(prob, godunov_splitting_callback.solver, dt = integrator.dt)
    v = sol.u[end]
    #-----------------------
    # u = (v + I_0) / A with I_0(x) = ∫_{x_min}^x A (Q - c) dy and c = ∫ A Q / ∫ A, which follows
    # from (A u)_x = -A (p0_t + I_w k_w (T - T_u) / √A) / (γ p0) with p0_t = γ p0 c
    # (optional source term s_p of the p0 equation enters the constraint as well)
    Q =
        add_source(
            - I_w_nodes .* equations.k_w .* A_sqrt_A_inv .* (T_prev .- Tu_prev),
            s_p_prev,
        ) ./ (equations.γ .* p0_prev[first_node])
    c = sum(Q .* A_nodes .* W) / sum(A_nodes .* W)
    I_0 = dg_cumulative_integral(A_nodes .* (Q .- c), ops)
    v1 = (v .+ I_0) ./ A_nodes
    v1_dx = dg_derivative(v1, ops)
    #------------------------
    p0_dt = p0_rhs(
        p0,
        v1,
        v1_dx,
        T,
        Tu,
        A_nodes,
        A_x_nodes,
        A_sqrt_A_inv,
        I_w_nodes,
        W,
        equations,
        s_p,
    )
    LI = linear_interpolation(
        [t_prev, t_now],
        [p0_dt_prev, p0_dt],
        extrapolation_bc = Line(),
    )
    prob = ODEProblem((u, p, t) -> LI(t), p0_prev[first_node], (t_prev, t_now))
    sol = solve(prob, godunov_splitting_callback.solver, dt = integrator.dt)
    p0_new = sol.u[end]
    #-----------------------
    u_now[2, :, :] .= v1
    u_now[3, :, :] .= p0_new
    #-----------------------
    return integrator
end

@doc raw"""
    v_rhs(rho, v1, v1_dx, x, A_nodes, W, equations, s_v = nothing)

Right-hand side of the ODE for the auxiliary velocity ``v``,
```math
v' = \frac{\int_0^1 F_v \, dx}{\int_0^1 \rho / A \, dx}, \qquad
F_v = -\rho u u_x - \frac{\beta}{A} \rho u \left(\eta + (1-\eta)|u|\right)
      - \frac{h_x}{Fr^2} (\rho - \rho_{h_0}) \; (+ s_v).
```
The nodal fields `rho`, `v1` (``u``), `v1_dx` (``u_x``), `x` and `A_nodes` are
`(node, element)` matrices; both integrals are evaluated with the DG quadrature weights `W`.
`s_v` is an optional nodal source (`nothing`: no source).
"""
function v_rhs(rho, v1, v1_dx, x, A_nodes, W, equations, s_v = nothing)
    @unpack β, η, ρₕ₀, Fr² = equations
    f = @. - rho * v1 * v1_dx - β / A_nodes * rho * v1 * (η + (1 - η) * abs(v1)) -
       h_x(x, (equations,)) * (rho - ρₕ₀) / Fr²
    return sum(add_source(f, s_v) .* W) / sum(rho ./ A_nodes .* W)
end

@doc raw"""
    p0_rhs(p0, v1, v1_dx, T, Tu, A_nodes, A_x_nodes, A_sqrt_A_inv, I_w_nodes, W, equations, s_p = nothing)

Right-hand side of the ODE for the spatially constant pressure ``p_0``, the ``A``-weighted
domain average
```math
p_0' = \frac{\int_0^1 A f_p \, dx}{\int_0^1 A \, dx}, \qquad
f_p = -\gamma p_0 u_x - \gamma \frac{A_x}{A} u p_0 - \frac{I_w k_w}{\sqrt{A}} (T - T_u)
      \; (+ s_p).
```
This follows from integrating the pointwise ``p_0`` equation over the domain: the term
``-\gamma p_0 (A u)_x / A`` vanishes after multiplication by ``A`` and integration
(periodicity). The nodal fields are `(node, element)` matrices; the integrals are evaluated
with the DG quadrature weights `W`. `s_p` is an optional nodal source (`nothing`: no source).
"""
function p0_rhs(
    p0,
    v1,
    v1_dx,
    T,
    Tu,
    A_nodes,
    A_x_nodes,
    A_sqrt_A_inv,
    I_w_nodes,
    W,
    equations,
    s_p = nothing,
)
    @unpack γ, k_w = equations
    f = @. - γ * p0 * v1_dx - γ * A_x_nodes / A_nodes * v1 * p0 -
       I_w_nodes * k_w * A_sqrt_A_inv * (T - Tu)
    return sum(add_source(f, s_p) .* A_nodes .* W) / sum(A_nodes .* W)
end

@doc raw"""
    ode_source_terms(source_terms, x, t, equations)

Evaluates the optional sources of the ODE step at the nodes `x` (`(node, element)` matrix) and
time `t`. `source_terms(x, t, equations)` returns `SVector(s_p, s_v)`:
- `s_p` is added pointwise to the right-hand side of the ``p_0`` equation ([`p0_rhs`](@ref))
  and hence also to ``Q`` in the velocity constraint,
- `s_v` is added to the integrand of the ``v`` equation ([`v_rhs`](@ref)).

Returns the nodal matrices `(s_p, s_v)`, or `(nothing, nothing)` if
`source_terms === nothing` (no sources, no extra work).
"""
ode_source_terms(::Nothing, x, t, equations) = (nothing, nothing)
function ode_source_terms(source_terms, x, t, equations)
    s = source_terms.(x, t, Ref(equations))
    return (getindex.(s, 1), getindex.(s, 2))
end

# add an optional nodal source `s` to the nodal field `f`
add_source(f, ::Nothing) = f
add_source(f, s) = f .+ s
end # @muladd
