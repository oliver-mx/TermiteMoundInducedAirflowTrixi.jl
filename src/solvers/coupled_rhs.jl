# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# "coupled" method: coupled method of lines for the TME. The transport PDE for rho and the ODE
# for Ti at every node (Trixi) and the ODEs for v and p0 are one ODE system, integrated together by OrdinaryDiffEq (no operator
# splitting). The velocity u is the algebraic variable of the system: it is eliminated by the
# velocity constraint in every right-hand side evaluation.
#
# The ODE state is the usual Trixi state (rho, v1, p0, Ti):
# - p0 is spatially constant. Every node gets the same dp0/dt, so it stays exactly uniform.
# - v is stored in the v1 field: v = A(x_min) v1(x_min) (the cumulative integral I_0 vanishes at
#   the first node). Every node gets dv1/dt = (dv/dt) / A, i.e. only the first node carries
#   information, the other v1 entries are recomputed from the state.
# - In every RHS evaluation (every Runge-Kutta stage) the velocity u = (v + I_0) / A is
#   computed from the current state and inserted into a copy of the state, on which the
#   Trixi RHS (rho, Ti) and the ODE right-hand sides (v, p0) are evaluated.
# - `AlgebraicVelocityCallback` writes the algebraic velocity u into the v1 field after every
#   step (a pure function of the state, no splitting), so that the step size control, analysis
#   and output see the true velocity.

@doc raw"""
    CoupledSourceTerms(source_terms; ode_source_terms = nothing)

Select the "coupled" method, a coupled method of lines for the TME, when passed as
`source_terms` to `SemidiscretizationHyperbolic`. The ODEs for ``v`` and ``p_0`` are part of
the right-hand side `Trixi.rhs!` and are integrated together with the transport PDE for
``\rho`` and the ODE for ``T_i`` at every node by OrdinaryDiffEq (no operator splitting). The
velocity ``u`` is the algebraic variable of the system: it is computed from the velocity
constraint in every right-hand side evaluation (every Runge-Kutta stage).

In every right-hand side evaluation:
1. the state is copied into a buffer and the velocity ``u = (v + I_0)/A`` is computed from it
   ([`compute_algebraic_velocity!`](@ref)),
2. the Trixi right-hand side of ``\rho`` and ``T_i`` is evaluated on this state with the inner
   `source_terms`,
3. the right-hand sides of ``v`` and ``p_0`` are evaluated ([`ode_rhs`](@ref)). Every node
   gets ``\mathrm{d}p_0/\mathrm{d}t`` in the `p` field and ``(\mathrm{d}v/\mathrm{d}t)/A`` in the
   `v1` field (``v = A(x_{min}) \, v_1(x_{min})`` is stored at the first node).

Operators and work arrays are precomputed when the semidiscretization is created, so the
right-hand side does not allocate.

# Arguments
- `source_terms`: source terms of the ``\rho`` and ``T_i`` equations (e.g.
  [`source_terms`](@ref)).
- `ode_source_terms`: optional additional sources `(x, t, equations) -> SVector(s_p, s_v)` of
  the ``p_0`` and ``v`` equations (e.g. [`manufactured_source_terms_ode`](@ref)), or `nothing`
  (default) for no sources. These sources are only used for the manufactured solution.

# Usage
Use it together with the [`AlgebraicVelocityCallback`](@ref), placed before the
`StepsizeCallback`:
```julia
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, dg,
                                    source_terms = CoupledSourceTerms(source_terms),
                                    boundary_conditions = boundary_condition_periodic)
algebraic_velocity_callback = AlgebraicVelocityCallback()
callbacks = CallbackSet(summary_callback, algebraic_velocity_callback, stepsize_callback)
```

!!! note
    For the "splitting" method, pass the plain `source_terms` instead and use the
    [`GodunovSplittingCallback`](@ref).
"""
struct CoupledSourceTerms{SourceTerms,ODESourceTerms}
    source_terms::SourceTerms
    ode_source_terms::ODESourceTerms
end

function CoupledSourceTerms(source_terms; ode_source_terms = nothing)
    return CoupledSourceTerms(source_terms, ode_source_terms)
end

# Precomputed operators and work buffers of the coupled RHS, stored as `cache.termite`.
# They are only used if the semidiscretization is created with `CoupledSourceTerms`.
function Trixi.create_cache(
    mesh::TreeMesh{1},
    equations::TermiteMoundEquations1D,
    dg::DGSEM,
    RealT,
    uEltype,
)
    cache = invoke(
        Trixi.create_cache,
        Tuple{TreeMesh{1},Any,Trixi.DG,Any,Any},
        mesh,
        equations,
        dg,
        RealT,
        uEltype,
    )
    return (;
        cache...,
        termite = create_termite_cache(mesh, equations, dg, cache, uEltype),
    )
end

function create_termite_cache(mesh, equations, dg, cache, uEltype)
    ops = dg_operators(dg, cache)
    x = cache.elements.node_coordinates[1, :, :]
    A_nodes = A.(x, Ref(equations))
    AW = A_nodes .* ops.W
    n_nodes, n_elements = size(x)
    u_buffer =
        Vector{uEltype}(undef, Trixi.nvariables(equations) * n_nodes * n_elements)
    return (;
        ops.W,
        D = Matrix(ops.D),
        ops.S,
        ops.inverse_jacobian,
        ops.order,
        first_node = CartesianIndex(1, ops.order[1]),
        x,
        A_nodes,
        A_inv = inv.(A_nodes),
        A_x_A = A_x.(x, Ref(equations)) ./ A_nodes,
        A_sqrt_A_inv = inv.(sqrt.(A_nodes)),
        I_w_nodes = I_w.(x, Ref(equations)),
        h_x_nodes = h_x.(x, Ref(equations)),
        AW,
        sum_AW = sum(AW),
        u_buffer,
        v1 = zeros(uEltype, n_nodes, n_elements),
        v1_dx = zeros(uEltype, n_nodes, n_elements),
        Q = zeros(uEltype, n_nodes, n_elements),
        T = zeros(uEltype, n_nodes, n_elements),
        Tu = zeros(uEltype, n_nodes, n_elements),
        s_p = zeros(uEltype, n_nodes, n_elements),
        s_v = zeros(uEltype, n_nodes, n_elements),
    )
end

# Coupled right-hand side: algebraic velocity u, Trixi RHS for (rho, Ti), ODE RHS for (v, p0).
function Trixi.rhs!(
    du,
    u,
    t,
    mesh::TreeMesh{1},
    equations::TermiteMoundEquations1D,
    boundary_conditions,
    source_terms::CoupledSourceTerms,
    dg::DGSEM,
    cache,
)
    tc = cache.termite
    u_alg = Trixi.wrap_array(tc.u_buffer, mesh, equations, dg, cache)
    copyto!(u_alg, u)
    compute_algebraic_velocity!(u_alg, t, equations, source_terms.ode_source_terms, tc)

    invoke(
        Trixi.rhs!,
        Tuple{Any,Any,Any,TreeMesh{1},Any,Any,Any,Trixi.DG,Any},
        du,
        u_alg,
        t,
        mesh,
        equations,
        boundary_conditions,
        source_terms.source_terms,
        dg,
        cache,
    )

    dv, dp0 = ode_rhs(u_alg, equations, tc)
    @inbounds for e in axes(tc.v1, 2), i in axes(tc.v1, 1)
        du[2, i, e] = dv * tc.A_inv[i, e]
        du[3, i, e] = dp0
    end
    return nothing
end

@doc raw"""
    compute_algebraic_velocity!(u, t, equations, ode_source_terms, tc)

Overwrite the `v1` field of the state vector `u` (the node values of the velocity ``u``,
except at the first node) with the algebraic velocity ``u = (v + I_0(x)) / A`` determined by
the current state, where
```math
v = A(x_{min}) v_1(x_{min}), \quad I_0(x) = \int_{x_{min}}^x A (Q - c) \, dy, \quad
Q = \frac{-I_w k_w (T - T_u)/\sqrt{A} + s_p}{\gamma p_0}, \quad
c = \frac{\int A Q \, dx}{\int A \, dx}.
```
- The source ``s_p`` of the ``p_0`` equation is only used for the manufactured solution
  (``s_p = 0`` otherwise).
- ``T_u`` is evaluated with the time-only part computed once (`T_u_time`, `T_u_node`, see
  [`T_u`](@ref)).
- ``I_0`` is the exact integral of the elementwise interpolation polynomial of ``A (Q - c)``
  (elements from left to right), so ``I_0 = 0`` at the first node, and `v1` there stays
  exactly the stored value.
- The work arrays `Q`, `T`, `Tu`, `s_p`, `s_v`, `v1` and `v1_dx` (DG derivative of the
  velocity) of the cache `tc` are filled for [`ode_rhs`](@ref).

Used by the right-hand side of [`CoupledSourceTerms`](@ref) and by the
[`AlgebraicVelocityCallback`](@ref).
"""
function compute_algebraic_velocity!(u, t, equations, ode_source_terms::F, tc) where {F}
    @unpack γ, k_w = equations
    @unpack x, A_nodes, A_inv, A_sqrt_A_inv, I_w_nodes, AW, sum_AW, first_node = tc
    @unpack Q, T, Tu, s_p, s_v, v1 = tc
    p0 = u[3, first_node]
    v = A_nodes[first_node] * u[2, first_node]
    inv_γp0 = inv(γ * p0)

    T_u_time_t = T_u_time(t, equations)
    c = zero(eltype(Q))
    @inbounds for e in axes(Q, 2), i in axes(Q, 1)
        T[i, e] = p0 / u[1, i, e]
        Tu[i, e] = T_u_node(T_u_time_t, x[i, e], u[4, i, e], equations)
        s_p[i, e], s_v[i, e] = ode_sources(ode_source_terms, x[i, e], t, equations)
        Q[i, e] =
            (
                -I_w_nodes[i, e] * k_w * A_sqrt_A_inv[i, e] * (T[i, e] - Tu[i, e]) +
                s_p[i, e]
            ) * inv_γp0
        c += Q[i, e] * AW[i, e]
    end
    c /= sum_AW

    # v1 = (v + I_0) / A, with the cumulative integral I_0 of g = A (Q - c), integrating the
    # elementwise interpolation polynomial of g exactly (elements from left to right)
    @unpack S, inverse_jacobian, order = tc
    offset = zero(eltype(Q))
    @inbounds for e in order
        jacobian = inv(inverse_jacobian[e])
        for i in axes(Q, 1)
            integral = zero(eltype(Q))
            for j in axes(Q, 1)
                integral += S[i, j] * A_nodes[j, e] * (Q[j, e] - c)
            end
            v1[i, e] = (v + offset + jacobian * integral) * A_inv[i, e]
        end
        # I_0 at the right end of this element = offset of the next element
        integral = zero(eltype(Q))
        for j in axes(Q, 1)
            integral += S[end, j] * A_nodes[j, e] * (Q[j, e] - c)
        end
        offset += jacobian * integral
    end
    # keep the first node exactly (it stores v)
    v1[first_node] = u[2, first_node]
    @inbounds for e in axes(v1, 2), i in axes(v1, 1)
        u[2, i, e] = v1[i, e]
    end

    # elementwise derivative of the velocity
    @unpack D, v1_dx = tc
    @inbounds for e in axes(v1, 2)
        for i in axes(v1, 1)
            derivative = zero(eltype(v1))
            for j in axes(v1, 1)
                derivative += D[i, j] * v1[j, e]
            end
            v1_dx[i, e] = inverse_jacobian[e] * derivative
        end
    end
    return u
end

@doc raw"""
    ode_rhs(u, equations, tc)

Right-hand sides of the ODEs for ``v`` and ``p_0``,
```math
\frac{\mathrm{d}v}{\mathrm{d}t} = \frac{\int_0^1 (F_v + s_v) \, dx}{\int_0^1 \rho / A \, dx},
\qquad
\frac{\mathrm{d}p_0}{\mathrm{d}t} = \frac{\int_0^1 A (f_p + s_p) \, dx}{\int_0^1 A \, dx},
```
with the integrands ``F_v`` and ``f_p`` as in [`v_rhs`](@ref) and [`p0_rhs`](@ref) (same
formulas, computed here in one loop over the nodes with the DG quadrature weights). The state
`u` must have the velocity computed by [`compute_algebraic_velocity!`](@ref), which also fills
the work arrays of `tc`.

# Returns
`(dv, dp0)`.
"""
function ode_rhs(u, equations, tc)
    @unpack γ, k_w, β, η, ρₕ₀, Fr² = equations
    @unpack W, AW, sum_AW, A_inv, A_x_A, A_sqrt_A_inv, I_w_nodes, h_x_nodes = tc
    @unpack T, Tu, s_p, s_v, v1, v1_dx, first_node = tc
    p0 = u[3, first_node]
    int_v = zero(eltype(v1))
    int_rho_A = zero(eltype(v1))
    int_p0 = zero(eltype(v1))
    @inbounds for e in axes(v1, 2), i in axes(v1, 1)
        rho = u[1, i, e]
        vel = v1[i, e]
        vel_x = v1_dx[i, e]
        f_v =
            -rho * vel * vel_x -
            β * A_inv[i, e] * rho * vel * (η + (1 - η) * abs(vel)) -
            h_x_nodes[i, e] * (rho - ρₕ₀) / Fr² + s_v[i, e]
        int_v += f_v * W[i, e]
        int_rho_A += rho * A_inv[i, e] * W[i, e]
        f_p =
            -γ * p0 * vel_x - γ * A_x_A[i, e] * vel * p0 -
            I_w_nodes[i, e] * k_w * A_sqrt_A_inv[i, e] * (T[i, e] - Tu[i, e]) +
            s_p[i, e]
        int_p0 += f_p * AW[i, e]
    end
    return int_v / int_rho_A, int_p0 / sum_AW
end

# pointwise optional ODE sources (s_p, s_v)
@inline ode_sources(::Nothing, x, t, equations) = (zero(x), zero(x))
@inline function ode_sources(ode_source_terms::F, x, t, equations) where {F}
    s = ode_source_terms(x, t, equations)
    return (s[1], s[2])
end
end # @muladd
