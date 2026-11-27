# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

@doc raw"""
    manufactured_fields(x, t, equations::TermiteMoundEquations1D)

Exact fields and derivatives of the manufactured solution of the TME at the point `x`
(a scalar) and time `t`, in scaled variables.

# Forced problem
On the periodic domain ``x \in [0, 1]``, the manufactured solution solves the TME with
additional forcing terms ``s_\rho, s_T, s_p, s_v``:
```math
\begin{aligned}
\rho_t + (\rho u)_x &= -\frac{A_x}{A} \rho u + s_\rho, \\
T_{i,t} &= k_i (T - T_i) + s_T, \\
p_0'(t) &= -\gamma p_0 \left(u_x + \frac{A_x}{A} u\right)
           - \frac{I_w k_w}{\sqrt{A}} \left(T - T_u(t, x, T_i)\right) + s_p
           \qquad \text{for all } x, \\
\left(\int_0^1 \frac{\rho}{A} \, dx\right) v'(t) &= \int_0^1 \left[ -\rho u u_x
           - \frac{\beta}{A} \rho u \left(\eta + (1-\eta)|u|\right)
           - \frac{h_x}{Fr^2} (\rho - \rho_{h_0}) + s_v \right] dx, \\
T &= \frac{p_0}{\rho}, \qquad A(0) \, u(t, 0) = v(t).
\end{aligned}
```
The ``p_0`` equation holds pointwise. Since ``p_0`` does not depend on ``x``, it is
equivalent to the velocity constraint ``A u = v + \int_0^x A (Q - c) \, dy`` (with the
additional term ``s_p / (\gamma p_0)`` in ``Q``) together with the ODE ``p_0' = \gamma p_0 c``. The sources are given in
[`manufactured_source_terms`](@ref) (``s_\rho, s_T``) and
[`manufactured_source_terms_ode`](@ref) (``s_p, s_v``).

# Data
The model data of `equations` are used unchanged, so the manufactured solution works with any
cross section, height profile and boundary temperature:
- cross section ``A(x)``, ``A_x(x)`` (`equations.cross_section`, see [`A`](@ref)), where
  the examples use [`CrossSectionSine`](@ref)`(0.3)`, i.e.
  ``A(x) = 1 + 0.3 \sin(2\pi x + \pi/2)``, so that all ``A`` and ``A_x`` terms are active,
- the slope ``h_x`` (`equations.height_profile`, see [`h_x`](@ref)) and the heat transfer
  factor ``I_w`` ([`I_w`](@ref)),
- the boundary temperature ``T_u(t, x, T_i)`` (`equations.boundary_temperature`, see
  [`T_u`](@ref)).
For the manufactured solution, ``A`` must be ``1``-periodic and positive.

# Manufactured fields
With ``\omega = 2\pi/10`` (one period on ``t \in [0, 10]``) and ``k = 2\pi``:
```math
\begin{aligned}
\rho &= \rho_0 \left(1 + \varepsilon_\rho \sin(kx + \omega t)\right), &
p_0 &= 1 + \varepsilon_p \sin(\omega t), \\
v &= v_0 \left(1 + \varepsilon_v \sin(\omega t)\right), &
T_i &= T_{i,0} \left(1 + \varepsilon_T \cos(kx - \omega t)\right), \\
u &= \frac{v + w}{A}, &
w &= \varepsilon_u \left(1 + \tfrac{1}{2} \sin(\omega t)\right) \sin(kx),
\end{aligned}
```
with ``\varepsilon_\rho = 10^{-3}``, ``\varepsilon_p = 10^{-3}``, ``\varepsilon_v = 0.2``,
``\varepsilon_u = 0.1``, ``\varepsilon_T = 2 \cdot 10^{-3}``. The base values are those of
the physical initial state: ``\rho_0 = 1/T_0`` and ``T_{i,0} = T_0`` with the scaled `T0`,
and ``v_0`` the scaled `v0`.

These choices have the following properties:
- all fields are smooth and ``1``-periodic in ``x``, and ``p_0`` only depends on ``t``,
- ``A u = v + w`` is periodic and ``w(t, 0) = 0``, so ``A(0) u(t, 0) = v(t)`` holds exactly,
- ``u < 0`` everywhere (no sign change in ``|u|`` or in `flux_ranocha`),
- ``u_x = (w_x - A_x u)/A`` follows from ``(A u)_x = w_x``.

!!! note
    ``\varepsilon_\rho`` must stay small. A velocity perturbation ``\delta v`` advects the
    density gradient (``\delta\rho_t \approx -\delta v \, \rho_x / A``), which feeds back into
    ``v`` through the buoyancy term ``h_x \delta\rho / Fr^2`` (``Fr^2 \sim 5 \cdot 10^{-5}``).
    For ``\varepsilon_\rho = 10^{-2}``, this makes the manufactured trajectory unstable during
    part of each period (top-heavy density stratification), and the error grows to ``O(1)``
    independently of the mesh size and the time step.

# Returns
A `NamedTuple` with the fields and the derivatives needed for the source terms:
`rho, rho_t, rho_x, p0, p0_t, v, v_t, u, u_x, Ti, Ti_t`, and `A_`, `A_x_` (``A`` and
``A_x`` at `x`).
"""
@inline function manufactured_fields(x, t, equations::TermiteMoundEquations1D)
    ε_ρ = 0.001
    ε_p = 0.001
    ε_v = 0.2
    ε_u = 0.1
    ε_T = 0.002
    ω = 2 * pi / 10
    k = 2 * pi

    rho_0 = inv(temp2scaled(equations.T0, equations))
    Ti_0 = temp2scaled(equations.T0, equations)
    v_0 = vel2scaled(equations.v0, equations)

    s_t, c_t = sincos(ω * t)
    s_x, c_x = sincos(k * x)
    s_θ, c_θ = sincos(k * x + ω * t)
    s_φ, c_φ = sincos(k * x - ω * t)

    rho = rho_0 * (1 + ε_ρ * s_θ)
    rho_t = rho_0 * ε_ρ * ω * c_θ
    rho_x = rho_0 * ε_ρ * k * c_θ

    p0 = 1 + ε_p * s_t
    p0_t = ε_p * ω * c_t

    v = v_0 * (1 + ε_v * s_t)
    v_t = v_0 * ε_v * ω * c_t

    w = ε_u * (1 + 0.5f0 * s_t) * s_x
    w_x = ε_u * (1 + 0.5f0 * s_t) * k * c_x

    A_ = A(x, equations)
    A_x_ = A_x(x, equations)
    u = (v + w) / A_
    # (A u)_x = w_x
    u_x = (w_x - A_x_ * u) / A_

    Ti = Ti_0 * (1 + ε_T * c_φ)
    Ti_t = Ti_0 * ε_T * ω * s_φ

    return (; rho, rho_t, rho_x, p0, p0_t, v, v_t, u, u_x, Ti, Ti_t, A_, A_x_)
end

@doc raw"""
    ManufacturedInitialCondition(x, t, equations::TermiteMoundEquations1D)

Exact (manufactured) solution of the forced TME (see [`manufactured_fields`](@ref)) as the
state vector `(rho, v1, p, Ti)` = ``(\rho, u, p_0, T_i)`` at the position `x` (vector with one
entry) and time `t`.
Used as the initial condition (``t = 0``) and as the reference solution for the errors.

At ``t = 0`` the initial condition is
```math
\begin{aligned}
\rho(0, x) &= \rho_0 \left(1 + \varepsilon_\rho \sin(2\pi x)\right), &
u(0, x) &= \frac{v_0 + \varepsilon_u \sin(2\pi x)}{A(x)}, \\
p_0(0) &= 1, &
T_i(0, x) &= T_{i,0} \left(1 + \varepsilon_T \cos(2\pi x)\right),
\end{aligned}
```
so that ``v(0) = A(0) \, u(0, 0) = v_0``. The constants are defined in
[`manufactured_fields`](@ref).

Use it together with [`manufactured_source_terms`](@ref) (sources of ``\rho`` and ``T_i``) and
[`manufactured_source_terms_ode`](@ref) (sources of ``p_0`` and ``v``), for the "splitting"
([`GodunovSplittingCallback`](@ref)) and the "coupled" ([`CoupledSourceTerms`](@ref)) methods.

# Returns
`SVector(rho, v1, p, Ti)`.
"""
@inline function ManufacturedInitialCondition(x, t, equations::TermiteMoundEquations1D)
    m = manufactured_fields(x[1], t, equations)
    return SVector(m.rho, m.u, m.p0, m.Ti)
end

@doc raw"""
    manufactured_source_terms(u, x, t, equations::TermiteMoundEquations1D)

Source terms of the ``\rho`` and ``T_i`` equations of the forced TME (handled by Trixi, see
[`manufactured_fields`](@ref)): the model source terms
[`source_terms`](@ref), evaluated on the numerical state `u` at the position `x` (vector with
one entry), plus the manufactured sources
``s_\rho`` and ``s_T``, evaluated on the exact fields.

The forced problem reads as follows (the ``p_0`` and ``v`` sources are given by
[`manufactured_source_terms_ode`](@ref)):
```math
\begin{aligned}
\rho_t + (\rho u)_x &= -\frac{A_x}{A} \rho u + s_\rho, \\
T_{i,t} &= k_i (T - T_i) + s_T, \\
p_0' &= -\gamma p_0 \left(u_x + \frac{A_x}{A} u\right)
        - \frac{I_w k_w}{\sqrt{A}} (T - T_u) + s_p \qquad \text{for all } x, \\
\left(\int_0^1 \frac{\rho}{A} \, dx\right) v' &= \int_0^1 \left[ F_v + s_v \right] dx,
\end{aligned}
```
with ``T = p_0/\rho`` and the momentum integrand
``F_v = -\rho u u_x - \frac{\beta}{A} \rho u (\eta + (1-\eta)|u|) - \frac{h_x}{Fr^2}(\rho - \rho_{h_0})``.
Inserting the exact fields ``\hat\rho, \hat u, \hat p_0, \hat T_i`` (with
``\hat T = \hat p_0 / \hat\rho``) gives the sources
```math
\begin{aligned}
s_\rho &= \hat\rho_t + \hat\rho_x \hat u + \hat\rho \hat u_x + \frac{A_x}{A} \hat\rho \hat u, \\
s_T &= \hat T_{i,t} - k_i (\hat T - \hat T_i).
\end{aligned}
```

# Returns
`source_terms(u, x, t, equations) + SVector(s_rho, 0, 0, s_Ti)`. The `v1` and `p` components
stay zero: they are updated by the ODE part of the system ([`GodunovSplittingCallback`](@ref)
or [`CoupledSourceTerms`](@ref)).
"""
@inline function manufactured_source_terms(u, x, t, equations::TermiteMoundEquations1D)
    m = manufactured_fields(x[1], t, equations)
    T = m.p0 / m.rho
    s_rho = m.rho_t + m.rho_x * m.u + m.rho * m.u_x + m.A_x_ / m.A_ * m.rho * m.u
    s_Ti = m.Ti_t - equations.k_i * (T - m.Ti)
    return source_terms(u, x, t, equations) + SVector(s_rho, 0.0, 0.0, s_Ti)
end

@doc raw"""
    manufactured_source_terms_ode(x, t, equations::TermiteMoundEquations1D)

Manufactured sources ``(s_p, s_v)`` of the ``p_0`` and ``v`` equations of the forced TME (see
[`manufactured_fields`](@ref)) at the position `x` (a scalar) and time `t`. With the exact
fields ``\hat\rho, \hat u, \hat p_0, \hat v, \hat T_i``, ``\hat T = \hat p_0 / \hat\rho`` and
``\hat T_u = T_u(t, x, \hat T_i)``:
```math
\begin{aligned}
s_p &= \hat p_0' - \left[ -\gamma \hat p_0 \hat u_x - \gamma \frac{A_x}{A} \hat u \hat p_0
       - \frac{I_w k_w}{\sqrt{A}} \left(\hat T - \hat T_u\right) \right], \\
s_v &= \frac{\hat\rho}{A} \hat v' - \hat F_v, \qquad
\hat F_v = -\hat\rho \hat u \hat u_x - \frac{\beta}{A} \hat\rho \hat u
           \left(\eta + (1-\eta)|\hat u|\right) - \frac{h_x}{Fr^2} (\hat\rho - \rho_{h_0}).
\end{aligned}
```
- ``s_p`` is added pointwise to the right-hand side of the ``p_0`` equation, and hence, as
  ``s_p / (\gamma p_0)``, also to ``Q`` in the velocity constraint. It must depend on ``x``: only then can the manufactured
  velocity ``\hat u`` be chosen freely and smoothly. ``s_p`` is uniquely determined by the
  fields.
- ``s_v`` is added to the integrand of the ``v`` equation. Only its integral matters. The
  pointwise choice above makes the discrete integrand equal to ``\hat\rho \hat v' / A`` at
  every node, so the ``v`` equation is reproduced without quadrature error.

# Usage
Pass it as `source_terms` to [`GodunovSplittingCallback`](@ref) ("splitting") or as
`ode_source_terms` to [`CoupledSourceTerms`](@ref) ("coupled").

# Returns
`SVector(s_p, s_v)`.
"""
@inline function manufactured_source_terms_ode(x, t, equations::TermiteMoundEquations1D)
    @unpack γ, k_w, β, η, ρₕ₀, Fr² = equations
    m = manufactured_fields(x, t, equations)
    @unpack rho, p0, u, u_x, A_, A_x_ = m
    T = p0 / rho
    Tu = T_u(t, x, m.Ti, equations)
    # p0 equation: s_p = p0_t - f, f = right-hand side as in `p0_rhs`
    f_p =
        - γ * p0 * u_x - γ * A_x_ / A_ * u * p0 -
        I_w(x, equations) * k_w * inv(sqrt(A_)) * (T - Tu)
    s_p = m.p0_t - f_p
    # v equation: s_v = rho/A v_t - f, f = integrand as in `v_rhs`
    f_v =
        - rho * u * u_x - β / A_ * rho * u * (η + (1 - η) * abs(u)) -
        h_x(x, equations) * (rho - ρₕ₀) / Fr²
    s_v = rho / A_ * m.v_t - f_v
    return SVector(s_p, s_v)
end
end # @muladd
