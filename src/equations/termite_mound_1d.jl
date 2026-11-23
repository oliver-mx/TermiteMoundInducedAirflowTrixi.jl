# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent
@doc raw"""
    TermiteMoundEquations1D(; γ, k_i, k_w, tᵣ, uᵣ, xa, xb, xc, r, h, L, β, η, Fr², ρₕ₀,
                            T_ref, t_ref, T0, v0, Ti_LI,
                            cross_section, height_profile, boundary_temperature)

Termite mound equations (TME) in one space dimension: a low Mach number model of the airflow
through the closed flow channel of a termite mound (soil → flute → chimney → soil),
parametrized by ``x \in [0, 1]`` (periodic). All quantities are scaled.

### Equations
```math
\begin{aligned}
\rho_t + (\rho u)_x &= -\frac{A_x}{A} \rho u, \\
T_{i,t} &= k_i (T - T_i), \\
v'(t) &= \frac{1}{\int_0^1 \frac{\rho}{A} \, dx} \int_0^1 \left[ -\rho u u_x
         - \frac{\beta}{A} \rho u \left(\eta + (1-\eta)|u|\right)
         - \frac{h_x}{Fr^2} (\rho - \rho_{h_0}) \right] dx, \\
p_0'(t) &= -\gamma p_0 \left(u_x + \frac{A_x}{A} u\right)
           - \frac{I_w k_w}{\sqrt{A}} \left(T - T_u(t, x, T_i)\right)
           \qquad \text{for all } x, \\
p_0 &= \rho T.
\end{aligned}
```
Since ``p_0`` does not depend on ``x``, the ``p_0`` equation determines the velocity. With
``Q = -I_w k_w (T - T_u)/(\sqrt{A} \gamma p_0)`` and ``c = \int_0^1 A Q \, dx / \int_0^1 A \, dx``:
```math
p_0' = \gamma p_0 c, \qquad
u(t, x) = \frac{v(t)}{A(x)} + \frac{1}{A(x)} \int_0^x A (Q - c) \, dy,
```
so that ``v = A(0) u(t, 0)`` is the volume flux at ``x = 0``.

### Unknowns and state vector
The unknowns are the air density ``\rho``, the air velocity ``u``, the leading-order pressure
``p_0`` (spatially constant), the mound material temperature ``T_i`` and the auxiliary
velocity ``v``. The state vector at every node is `(rho, v1, p, Ti)` =
``(\rho, u, p_0, T_i)``; ``v`` is not stored separately (see [`GodunovSplittingCallback`](@ref)
and [`CoupledSourceTerms`](@ref)).

### Model data
The cross section, the height profile and the boundary temperature are passed as objects
(required keyword arguments, see below); the heat transfer factor is fixed:
- cross section ``A(x)`` and ``A_x(x)`` (`cross_section`, see [`A`](@ref)),
- slope of the height profile ``h_x(x)`` (`height_profile`, see [`h_x`](@ref)),
- boundary temperature ``T_u(t, x, T_i)`` (`boundary_temperature`, see [`T_u`](@ref)),
- spatially variable heat transfer coefficient ``\alpha_w I_w(x)`` ([`I_w`](@ref)).

### Usage
```julia
(γ, k_i, k_w, tᵣ, uᵣ, xa, xb, xc, r, h, L, β, η, Fr², ρₕ₀, T_ref, t_ref, T0, v0, Ti_LI) =
    termite_parameters(1.5, 2.5)
cross_section = CrossSectionConstant()
height_profile = HeightProfileTermiteMound(; xa, xb, L, h)
boundary_temperature = BoundaryTemperatureTermiteMound(; xa, xb, xc, tᵣ, T_ref)
equations = TermiteMoundEquations1D(; γ, k_i, k_w, tᵣ, uᵣ, xa, xb, xc, r, h, L, β, η, Fr²,
                                    ρₕ₀, T_ref, t_ref, T0, v0, Ti_LI,
                                    cross_section, height_profile, boundary_temperature)
```

### Keyword arguments
The scalar arguments are usually created by [`termite_parameters`](@ref). Their defaults are
those of `termite_parameters(0.6, 2.0)` (mound radius 0.6 m, height 2 m), except
`v0 = -2.219` and a simple three-point `Ti_LI`.
- `γ`: adiabatic exponent.
- `k_i`: heat exchange rate between air and mound material (``T_i`` equation).
- `k_w`: wall heat transfer coefficient (``p_0`` equation, contains ``\alpha_w``).
- `tᵣ`: reference time ``L/u_r`` (s). `uᵣ`: reference velocity (m/s).
- `xa`, `xb`, `xc`: scaled positions of the transitions soil → flute, flute → chimney (top
  of the mound) and chimney → soil (see [`mound_geometry`](@ref)).
- `r`, `h`: radius and height of the mound (m). `L`: length of the flow channel (m).
- `β`: friction coefficient. `η`: viscosity.
- `Fr²`: squared Froude number ``u_r^2/(gL)``. `ρₕ₀`: scaled reference density of the
  buoyancy term.
- `T_ref`: reference temperature (K). `t_ref`: its inverse ``1/T_{ref}`` (1/K), used to scale
  ``T_u`` (which is computed in K).
- `T0`: initial air temperature (°C). `v0`: initial velocity (cm/s).
- `Ti_LI`: initial mound material temperature (°C) as an interpolation object on ``[0, 1]``.
- `cross_section` (required): cross section of the flow channel, e.g.
  [`CrossSectionConstant`](@ref) (``A \equiv 1``), [`CrossSectionSine`](@ref) or a
  user-defined [`CrossSection`](@ref).
- `height_profile` (required): slope of the height profile, e.g.
  [`HeightProfileTermiteMound`](@ref) or a user-defined [`HeightProfile`](@ref).
- `boundary_temperature` (required): boundary (wall) temperature, e.g.
  [`BoundaryTemperatureTermiteMound`](@ref) or a user-defined [`BoundaryTemperature`](@ref).

!!! note
    The positions `xa`, `xb`, `xc` and the length `L` are used both by the equations
    (``I_w``, plotting of ``T_i``) and by the termite mound height profile and boundary
    temperature. Create all of them from the same [`termite_parameters`](@ref) call.
"""
struct TermiteMoundEquations1D{
    RealT<:Real,
    CrossSectionType,
    HeightProfileType,
    BoundaryTemperatureType,
} <: AbstractEquations{1,4}
    γ::RealT
    inv_gamma_minus_one::RealT
    k_i::RealT
    k_w::RealT
    tᵣ::RealT
    uᵣ::RealT
    xa::RealT
    xb::RealT
    xc::RealT
    r::RealT
    h::RealT
    L::RealT
    β::RealT
    η::RealT
    Fr²::RealT
    ρₕ₀::RealT
    T_ref::RealT
    t_ref::RealT
    T0::RealT
    v0::RealT
    Ti_LI::Interpolations.Extrapolation
    cross_section::CrossSectionType
    height_profile::HeightProfileType
    boundary_temperature::BoundaryTemperatureType
    function TermiteMoundEquations1D(;
        γ = 1.3987529976019184,
        k_i = 0.0006167609978013738,
        k_w = 5.285934848600675,
        tᵣ = 101.17763390064094,
        uᵣ = 0.05,
        xa = 0.1326025065794086,
        xb = 0.5453540646624214,
        xc = 0.9822095068780948,
        r = 0.6,
        h = 2.0,
        L = 5.058881695032047,
        β = 36.397915028380886,
        η = 0.93,
        Fr² = 5.037516457669399e-5,
        ρₕ₀ = 0.982,
        T_ref = 297.76980029566033,
        t_ref = 0.0033582989242263127,
        T0 = 27.0,
        v0 = -2.219,
        Ti_LI = linear_interpolation(
            [0.0, 0.5344827593241117, 1.0],
            [27.53453, 28.98715, 27.55033],
        ),
        cross_section,
        height_profile,
        boundary_temperature,
    )
        γ, inv_gamma_minus_one = promote(γ, inv(γ - 1.0))
        new{
            typeof(β),
            typeof(cross_section),
            typeof(height_profile),
            typeof(boundary_temperature),
        }(
            γ,
            inv_gamma_minus_one,
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
    end
end
# Reference for the TME:
#- Marx, Oliver P. and Gasser, Ingenuin and Schmidgall, Annika (2026)
#  ...
#  [DOI: ...](https://doi.org/...)

@doc raw"""
    termite_parameters(r, height; g=9.81, R=8.314, M=0.0289652, c_v=20.85, α_w = 1.1, D0 = 0.05, ρₕ₀ = 0.982, η = 0.93, uᵣ = 0.05, ρᵣ = 1.17, pᵣ = 1, c_v_ = 790, κ = 0.184, T_ref = 297.76980029566033, t_ref = 0.0033582989242263127, T0 = 27.0, v0 = -5.0)

Computes the (scaled) parameters of the termite mound model for a mound with radius `r` and
height `height` (both in m).

The main steps are:
- geometry: channel length ``L`` and transitions ``x_a, x_b, x_c`` ([`mound_geometry`](@ref)),
- reference scales: ``x_r = L``, ``t_r = L/u_r``, ``Fr^2 = u_r^2/(gL)``,
- adiabatic exponent ``\gamma = (c_v + R)/c_v`` with the specific gas constant ``R/M`` and
  ``c_v/M``,
- friction coefficient ``\beta = \lambda_w x_r \sqrt{\pi}/(4\sqrt{A_r})`` with the pipe
  friction coefficient ``\lambda_w = 64/Rey``, the Reynolds number ``Rey`` and the cross section
  ``A_r = \pi (D_0/2)^2``,
- wall heat transfer coefficient
  ``k_w = (\gamma - 1) x_r \alpha_w T_r \sqrt{\pi}/(\sqrt{A_r} u_r p_r)``,
- heat exchange rate ``k_i = t_r \kappa L_i/(c_i m)`` with the thermal conductivity
  ``\kappa`` (`κ`), the specific heat capacity ``c_i`` (`c_v_`) and the mass ``m`` of the
  mound material, and the length ``L_i = x_r (x_c - x_a)`` of the channel inside the mound,
- initial mound material temperature profile ``T_i`` ([`Get_initial_Ti`](@ref)).

### Parameters
- `r`: radius of the mound (m).
- `height`: height of the mound (m).

### Optional parameters
- `g`: gravitational acceleration (default: 9.81 m/s²).
- `R`: universal gas constant (default: 8.314 J/(mol·K)).
- `M`: molar mass of air (default: 0.0289652 kg/mol).
- `c_v`: molar heat capacity of air at constant volume (default: 20.85 J/(mol·K)).
- `α_w`: wall heat transfer coefficient (default: 1.1). It is made variable in space by the
  factor [`I_w`](@ref) in the equations.
- `D0`: diameter of the flow channel (default: 0.05 m).
- `ρₕ₀`: scaled reference density of the buoyancy term (default: 0.982).
- `η`: viscosity (default: 0.93).
- `uᵣ`: reference velocity (default: 0.05 m/s).
- `ρᵣ`: reference density (default: 1.17 kg/m³).
- `pᵣ`: reference pressure (default: 1 Pa).
- `c_v_`: specific heat capacity of the mound material (default: 790 J/(kg·K)).
- `κ`: thermal conductivity of the mound material (default: 0.184 W/(m·K)).
- `T_ref`: reference temperature (default: 297.76980029566033 K).
- `t_ref`: inverse reference temperature ``1/T_{ref}`` (default: 0.0033582989242263127 1/K).
- `T0`: initial air temperature (default: 27.0 °C).
- `v0`: initial velocity (default: -5.0 cm/s).

### Returns
A 20-tuple `(γ, k_i, k_w, tᵣ, uᵣ, xa, xb, xc, r, h, L, β, η, Fr², ρₕ₀, T_ref, t_ref, T0, v0,
Ti_LI)` of the scalar keyword arguments of [`TermiteMoundEquations1D`](@ref). The model data
`cross_section`, `height_profile` and `boundary_temperature` are not included and are passed
separately (e.g. [`HeightProfileTermiteMound`](@ref) and
[`BoundaryTemperatureTermiteMound`](@ref), built from the returned `xa`, `xb`, `xc`, `L`, `h`,
`tᵣ` and `T_ref`).
"""
function termite_parameters(
    r,
    height;
    g = 9.81,
    R = 8.314,
    M = 0.0289652,
    c_v = 20.85,
    α_w = 1.1,
    D0 = 0.05,
    ρₕ₀ = 0.982,
    η = 0.93,
    uᵣ = 0.05,
    ρᵣ = 1.17,
    pᵣ = 1,
    c_v_ = 790,
    κ = 0.184,
    T_ref = 297.76980029566033,
    t_ref = 0.0033582989242263127,
    T0 = 27.0,
    v0 = -5.0,
)
    h = height
    L, xa, xb, xc = mound_geometry(r, h)
    Aᵣ = (D0 / 2)^2 * pi
    xᵣ = L
    tᵣ = L / uᵣ
    Re = 9000 * 2 * sqrt(Aᵣ) / (xᵣ * sqrt(pi))
    R = R / M
    cᵥ = c_v / M
    γ = (cᵥ + R) / cᵥ
    Fr² = (uᵣ^2) / (g * L)
    𝑇ᵣ = pᵣ / (ρᵣ * R)
    λ_w = 64 / Re
    β = (λ_w * xᵣ * sqrt(pi)) / (2 * 2 * sqrt(Aᵣ))
    k_w = (γ-1) * (xᵣ * α_w * 𝑇ᵣ * sqrt(pi)) / (sqrt(Aᵣ) * uᵣ * pᵣ)
    mass = (1742.448 * (1/3 * pi * r * r * h)) ./ 8
    Li = (xᵣ*(xc-xa))
    k_i = tᵣ * κ / (c_v_ * mass/Li)
    Ti_LI = Get_initial_Ti(r, h)
    return (
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
end

@doc raw"""
    mound_geometry(r, h)

Geometry of the flow channel of a mound with radius `r` and height `h` (both in m):
```math
\begin{aligned}
L &= h + 0.3 + \sqrt{0.09 + r^2} + \sqrt{r^2 + h^2}, \\
x_a &= \frac{\sqrt{0.09 + r^2}}{L}, \qquad
x_b = \frac{\sqrt{0.09 + r^2} + \sqrt{r^2 + h^2}}{L}, \qquad
x_c = \frac{L - 0.09}{L}.
\end{aligned}
```

### Returns
`(L, xa, xb, xc)`: the channel length `L` (m) and the scaled positions `xa` (soil → flute),
`xb` (flute → chimney, top of the mound) and `xc` (chimney → soil).
"""
function mound_geometry(r, h)
    L = h + 0.3 + sqrt(0.09 + r * r) + sqrt(r * r + h * h)
    xa = sqrt(0.09 + r * r) / L
    xb = (sqrt(0.09 + r * r) + sqrt(r * r + h * h)) / L
    xc = (L - 0.09) / L
    return L, xa, xb, xc
end

@doc raw"""
    Trixi.varnames(u, ::TermiteMoundEquations1D)

Names of the four variables of the state vector, used for output and plotting.

### Returns
`("rho", "v1", "p", "Ti")`: density ``\rho``, velocity ``u``, leading-order pressure ``p_0``
and mound material temperature ``T_i``.
"""
function varnames(u, ::TermiteMoundEquations1D)
    return ("rho", "v1", "p", "Ti")
end

@doc raw"""
    Trixi.cons2prim(u, ::TermiteMoundEquations1D)

Primitive variables of the state `u`. The state `(rho, v1, p, Ti)` is already stored in
primitive variables, so this is the identity.

### Returns
`u` unchanged.
"""
function Trixi.cons2prim(u, equations::TermiteMoundEquations1D)
    return u
end

@doc raw"""
    Trixi.PlotData1D(u, mesh::TreeMesh, equations::TermiteMoundEquations1D, solver::DGSEM, cache; kwargs...)

Same as Trixi's `PlotData1D`, except for ``T_i``. The mound material temperature is only
physically meaningful inside the mound, ``x_a < x < x_c`` (see [`T_u`](@ref)). Outside, the
plot data of `Ti` are set to `NaN`, so that plots only show the relevant part. The solution
itself is not changed.

This method is used by `PlotData1D(u_ode, semi)`, `PlotData1D(func, semi)` and
`PlotData1D(sol)`.
"""
function Trixi.PlotData1D(
    u,
    mesh::TreeMesh,
    equations::TermiteMoundEquations1D,
    solver::DGSEM,
    cache;
    kwargs...,
)
    pd = invoke(
        Trixi.PlotData1D,
        Tuple{Any,TreeMesh,Any,Any,Any},
        u,
        mesh,
        equations,
        solver,
        cache;
        kwargs...,
    )
    iTi = findfirst(==("Ti"), pd.variable_names)
    if !isnothing(iTi)
        outside = .!(equations.xa .< pd.x .< equations.xc)
        pd.data[outside, iTi] .= NaN
    end
    return pd
end

@doc raw"""
    Trixi.prim2cons(u, ::TermiteMoundEquations1D)

Converts primitive variables `(rho, v1, p, Ti)` to the variables of the state vector. Both are
identical, so this is the identity.

### Returns
`u` unchanged.
"""
function Trixi.prim2cons(u, equations::TermiteMoundEquations1D)
    return u
end

@doc raw"""
    Trixi.cons2entropy(u, ::TermiteMoundEquations1D)

Entropy variables of the state `u = (rho, v1, p, Ti)`, used by the `AnalysisCallback` (e.g.
for the output ``\sum \partial S/\partial U \cdot U_t``). The first three components are the
entropy variables of the compressible Euler equations, evaluated with ``(\rho, u, p_0)``:
```math
w_1 = \frac{\gamma - s}{\gamma - 1} - \frac{\rho u^2}{2 p_0}, \qquad
w_2 = \frac{\rho u}{p_0}, \qquad
w_3 = -\frac{\rho}{p_0}, \qquad
s = \log p_0 - \gamma \log \rho.
```
The fourth component is ``T_i``. This is a formal adaptation for the analysis output; it is
not an entropy of the TME.

### Returns
`SVector(w1, w2, w3, Ti)`.
"""
function Trixi.cons2entropy(u, equations::TermiteMoundEquations1D)
    rho, v1, p, Ti = u
    rho_v1 = rho * v1
    rho_e_total = p * equations.inv_gamma_minus_one + 0.5f0 * (rho_v1 * v1)
    v_square = v1^2
    s = log(p) - equations.γ * log(rho)
    rho_p = rho / p
    w1 = (equations.γ - s) * equations.inv_gamma_minus_one - 0.5f0 * rho_p * v_square
    w2 = rho_p * v1
    w3 = -rho_p
    return SVector(w1, w2, w3, Ti)
end

function Trixi.cons2cons(u, equations::TermiteMoundEquations1D)
    return u
end

@doc raw"""
    TermiteMoundInitialCondition(x, t, equations::TermiteMoundEquations1D)

Initial state of the termite mound simulations (independent of `t`):
```math
\rho = \frac{1}{T_0}, \qquad u = v_0, \qquad p_0 = 1, \qquad T_i = T_{i,LI}(x),
```
with the scaled initial air temperature ``T_0`` (from `T0` in °C), the scaled initial velocity
``v_0`` (from `v0` in cm/s) and the scaled initial mound material temperature profile
`Ti_LI` (in °C). Hence ``T = p_0/\rho = T_0`` everywhere.

The constant velocity ``u = v_0`` does not satisfy the velocity constraint in general. The
"splitting" method replaces it after the first time step, the "coupled" method at
initialization ([`AlgebraicVelocityCallback`](@ref)).

### Parameters
- `x`: position (vector with one entry).
- `t`: time (not used).
- `equations`: the [`TermiteMoundEquations1D`](@ref).

### Returns
`SVector(rho, v1, p, Ti)`.
"""
@inline function TermiteMoundInitialCondition(x, t, equations::TermiteMoundEquations1D)
    RealT = eltype(x)
    rho = inv(temp2scaled(equations.T0, equations))
    v1 = vel2scaled(equations.v0, equations)
    p0 = 1.0
    Ti_LI = equations.Ti_LI
    Ti = temp2scaled(Ti_LI(x[1]), equations)
    return SVector(rho, v1, p0, Ti)
end

@doc raw"""
    source_terms(u, x, t, equations::TermiteMoundEquations1D)

Source terms of the ``\rho`` and ``T_i`` equations of the TME (handled by Trixi):
```math
s = \left( -\frac{A_x}{A} \rho u, \; 0, \; 0, \; k_i (T - T_i) \right), \qquad T = \frac{p_0}{\rho}.
```
The first component is the cross section term of the ``\rho`` equation, the last one the heat
exchange between air and mound material. The `v1` and `p` components are zero: ``u`` and
``p_0`` are updated by the ODE part of the system ([`GodunovSplittingCallback`](@ref) or
[`CoupledSourceTerms`](@ref)).

### Parameters
- `u`: state vector `(rho, v1, p, Ti)`.
- `x`: position (vector with one entry).
- `t`: time (not used).
- `equations`: the [`TermiteMoundEquations1D`](@ref).

### Returns
`SVector(du1, 0, 0, du4)`.
"""
@inline function source_terms(u, x, t, equations::TermiteMoundEquations1D)
    rho, v1, p0, Ti = u
    x1 = x[1]
    T = p0 / rho
    du1 = - A_x(x1, equations) / A(x1, equations) * rho * v1
    du4 = equations.k_i * (T - Ti)
    return SVector(du1, 0.0, 0.0, du4)
end

@doc raw"""
    temp2scaled(y, equations::TermiteMoundEquations1D)

Converts a temperature `y` in °C into a scaled temperature, ``(y + 273.15)/T_{ref}``.
Works for scalars and arrays. Inverse of [`temp2unscaled`](@ref).

### Returns
The scaled temperature.
"""
@inline function temp2scaled(y, equations::TermiteMoundEquations1D)
    z = 273.15
    y = y .+ z
    x = y ./ equations.T_ref
    return x
end

@doc raw"""
    temp2unscaled(x, equations::TermiteMoundEquations1D)

Converts a scaled temperature `x` into a temperature in °C, ``x \, T_{ref} - 273.15``.
Works for scalars and arrays. Inverse of [`temp2scaled`](@ref).

### Returns
The temperature in °C.
"""
@inline function temp2unscaled(x, equations::TermiteMoundEquations1D)
    z = 273.15
    return x .* equations.T_ref .- z
end

@doc raw"""
    vel2scaled(y, equations::TermiteMoundEquations1D)

Converts a velocity `y` in cm/s into a scaled velocity, ``y/(100 \, u_r)`` with the reference
velocity ``u_r`` in m/s. Works for scalars and arrays. Inverse of [`vel2unscaled`](@ref).

### Returns
The scaled velocity.
"""
@inline function vel2scaled(y, equations::TermiteMoundEquations1D)
    return inv(equations.uᵣ * 100) .* y
end

@doc raw"""
    vel2unscaled(x, equations::TermiteMoundEquations1D)

Converts a scaled velocity `x` into a velocity in cm/s, ``100 \, u_r \, x``. Works for
scalars and arrays. Inverse of [`vel2scaled`](@ref).

### Returns
The velocity in cm/s.
"""
@inline function vel2unscaled(x, equations::TermiteMoundEquations1D)
    return x .* equations.uᵣ * 100
end

@doc raw"""
    A(x, equations::TermiteMoundEquations1D)
    A_x(x, equations::TermiteMoundEquations1D)

Cross section ``A(x)`` and its derivative ``A_x(x)`` of the flow channel at the scalar position
`x`, given by `equations.cross_section` (see [`CrossSectionConstant`](@ref),
[`CrossSectionSine`](@ref) and [`CrossSection`](@ref)).
"""
@inline A(x, equations::TermiteMoundEquations1D) = A(x, equations.cross_section)
@inline A_x(x, equations::TermiteMoundEquations1D) = A_x(x, equations.cross_section)

@doc raw"""
    h_x(x, equations::TermiteMoundEquations1D)

Slope ``h_x(x)`` of the height profile of the flow channel at the scalar position `x` (used in
the buoyancy term of the ``v`` equation), given by `equations.height_profile` (see
[`HeightProfileTermiteMound`](@ref) and [`HeightProfile`](@ref)).
"""
@inline h_x(x, equations::TermiteMoundEquations1D) = h_x(x, equations.height_profile)

@doc raw"""
    T_u(t, x, Ti, equations::TermiteMoundEquations1D)

Boundary (wall) temperature ``T_u`` in K at the scaled time `t`, the scalar position `x` and
the scaled mound material temperature `Ti` at `x`, given by `equations.boundary_temperature`
(see [`BoundaryTemperatureTermiteMound`](@ref) and [`BoundaryTemperature`](@ref)). Multiply by
`t_ref` ``= 1/T_{ref}`` to obtain the scaled temperature.

For efficiency, the solvers evaluate the time-only part once for all nodes,
`time_data = T_u_time(t, equations)`, and then `T_u_node(time_data, x, Ti, equations)` at
every node; both give the same result as `T_u(t, x, Ti, equations)`.
"""
@inline T_u(t, x, Ti, equations::TermiteMoundEquations1D) =
    T_u(t, x, Ti, equations.boundary_temperature)
@inline T_u_time(t, equations::TermiteMoundEquations1D) =
    T_u_time(t, equations.boundary_temperature)
@inline function T_u_node(time_data, x, Ti, equations::TermiteMoundEquations1D)
    return T_u_node(time_data, x, Ti, equations.boundary_temperature)
end

@doc raw"""
    I_w(x, equations::TermiteMoundEquations1D)

Spatial factor ``I_w(x)`` of the wall heat transfer coefficient. **It makes the heat transfer
coefficient ``\alpha_w`` variable in space**: the effective coefficient at ``x`` is
``\alpha_w I_w(x)``. Since ``k_w`` is proportional to ``\alpha_w`` (see
[`termite_parameters`](@ref)), the heat transfer term of the ``p_0`` equation reads
``-I_w(x) \, k_w (T - T_u)/\sqrt{A}``.
```math
I_w(x) = \begin{cases}
  0.5 & x \le x_a \text{ or } x > x_c \quad \text{(soil)}, \\
  1 & x_a < x \le x_c \quad \text{(mound: flute and chimney)}.
\end{cases}
```
So the heat transfer in the soil sections is half of that inside the mound.

### Parameters
- `x`: position (scalar).
- `equations`: the [`TermiteMoundEquations1D`](@ref).
"""
@inline function I_w(x, equations::TermiteMoundEquations1D)
    if x ≤ equations.xa
        return 0.5
    elseif x ≤ equations.xc
        return 1.0
    else
        return 0.5
    end
end

@doc raw"""
    Trixi.flux(u, orientation::Integer, ::TermiteMoundEquations1D)

Physical flux of the TME at a single point, ``f(u) = (\rho u, 0, 0, 0)``. Only the ``\rho``
equation is in conservation form. ``u`` and ``p_0`` are updated by the ODE part of the system,
and the ``T_i`` equation has no flux.

### Parameters
- `u`: state vector `(rho, v1, p, Ti)`.
- `orientation`: orientation (always `1` in 1D).
- `equations`: the [`TermiteMoundEquations1D`](@ref).
"""
@inline function Trixi.flux(u, orientation::Integer, equations::TermiteMoundEquations1D)
    rho, v1, _, _ = u
    f1 = rho * v1
    return SVector(f1, 0.0, 0.0, 0.0)
end

@doc raw"""
    flux_ranocha(u_ll, u_rr, orientation::Integer, ::TermiteMoundEquations1D)

Two-point volume flux for the flux differencing volume integral (`VolumeIntegralFluxDifferencing`),
built from logarithmic means as in Ranocha's flux for the Euler equations. Only the mass flux
is nonzero:
```math
f_1 = \rho^{\ln} \, \bar u, \qquad
\bar u = \begin{cases}
  u^{\ln} & u_{ll}, u_{rr} > 0, \\
  -(-u)^{\ln} & u_{ll}, u_{rr} < 0, \\
  \frac{1}{2}(u_{ll} + u_{rr}) & \text{otherwise},
\end{cases}
```
where ``a^{\ln}`` is the logarithmic mean of the left and right values (Trixi's `ln_mean`).
The logarithmic mean of the velocity is only defined for values of the same sign, so the
arithmetic mean is used if the sign changes. The flux is consistent: ``f_1 = \rho u`` for
``u_{ll} = u_{rr}``.

### Parameters
- `u_ll`, `u_rr`: left and right state vectors `(rho, v1, p, Ti)`.
- `orientation`: orientation (always `1` in 1D).
- `equations`: the [`TermiteMoundEquations1D`](@ref).

### Returns
`SVector(f1, 0, 0, 0)`.
"""
@inline function flux_ranocha(
    u_ll,
    u_rr,
    orientation::Integer,
    equations::TermiteMoundEquations1D,
)
    rho_ll, v1_ll, _, _ = u_ll
    rho_rr, v1_rr, _, _ = u_rr
    rho_mean = ln_mean(rho_ll, rho_rr)
    v1_mean = if v1_ll > 0.0 && v1_rr > 0.0
        ln_mean(v1_ll, v1_rr)
    elseif v1_ll < 0.0 && v1_rr < 0.0
        - ln_mean(-v1_ll, -v1_rr)
    else
        0.5f0 * (v1_ll + v1_rr)
    end
    f1 = rho_mean * v1_mean
    return SVector(f1, 0.0, 0.0, 0.0)
end

@doc raw"""
    Trixi.flux_hllc(u_ll, u_rr, orientation, ::TermiteMoundEquations1D)

Surface flux: the mass component of the HLLC flux (HLL with contact) of the compressible
Euler equations by E. F. Toro, evaluated with ``(\rho, u, p_0)`` (``p_0`` as pressure, total
energy ``p_0/(\gamma - 1) + \rho u^2/2``, sound speed ``\sqrt{\gamma p_0/\rho}``).
- [Lecture slides](http://www.prague-sum.com/download/2012/Toro_2-HLLC-RiemannSolver.pdf)
- Signal speeds (Einfeldt estimates with Roe averages):
  [DOI: 10.1137/S1064827593260140](https://doi.org/10.1137/S1064827593260140)

Only the mass component is computed, since only the ``\rho`` equation has a flux. It needs:
- the signal speeds ``S_L = \min(u_{ll} - c_{ll}, \tilde u - \tilde c)`` and
  ``S_R = \max(u_{rr} + c_{rr}, \tilde u + \tilde c)`` with the Roe averages ``\tilde u`` and
  ``\tilde c`` (the latter from the Roe-averaged enthalpy, hence the total energy),
- the contact speed ``S^*`` and the density ``\rho^*`` of the star region on the side of the
  interface.

The flux is then
```math
f_1 = \begin{cases}
  \rho_{ll} u_{ll} & S_L \ge 0, \\
  \rho_{rr} u_{rr} & S_R \le 0, \\
  \rho_{ll} u_{ll} + S_L (\rho^*_{ll} - \rho_{ll}) & S_L \le 0 \le S^*, \\
  \rho_{rr} u_{rr} + S_R (\rho^*_{rr} - \rho_{rr}) & \text{otherwise},
\end{cases}
\qquad \rho^*_K = \rho_K \frac{S_K - u_K}{S_K - S^*}.
```
The momentum and energy components of the Euler HLLC flux are not computed. The upwinding
through the acoustic signal speeds provides the numerical dissipation at element interfaces.

A second method with a `normal_direction::AbstractVector` rotates the velocity to the normal
direction, calls this method and scales the result by the length of the normal vector.

### Parameters
- `u_ll`, `u_rr`: left and right state vectors `(rho, v1, p, Ti)`.
- `orientation`: orientation (always `1` in 1D).
- `equations`: the [`TermiteMoundEquations1D`](@ref).

### Returns
`SVector(f1, 0, 0, 0)`.
"""
function Trixi.flux_hllc(
    u_ll,
    u_rr,
    orientation::Integer,
    equations::TermiteMoundEquations1D,
)
    rho_ll, v1_ll, p0_ll, _ = u_ll
    rho_rr, v1_rr, p0_rr, _ = u_rr

    # Sound speeds
    c_ll = sqrt(equations.γ * p0_ll / rho_ll)
    c_rr = sqrt(equations.γ * p0_rr / rho_rr)

    # Left and right mass fluxes rho * u
    f1_ll = rho_ll * v1_ll
    f1_rr = rho_rr * v1_rr

    # Roe averages of the velocity and the enthalpy (only needed for the signal speeds)
    rho_e_total_ll =
        p0_ll * equations.inv_gamma_minus_one + 0.5f0 * (rho_ll * v1_ll * v1_ll)
    rho_e_total_rr =
        p0_rr * equations.inv_gamma_minus_one + 0.5f0 * (rho_rr * v1_rr * v1_rr)
    sqrt_rho_ll = sqrt(rho_ll)
    sqrt_rho_rr = sqrt(rho_rr)
    sum_sqrt_rho = sqrt_rho_ll + sqrt_rho_rr
    vel_roe = (sqrt_rho_ll * v1_ll + sqrt_rho_rr * v1_rr) / sum_sqrt_rho
    ekin_roe = 0.5f0 * vel_roe^2
    H_ll = (rho_e_total_ll + p0_ll) / rho_ll
    H_rr = (rho_e_total_rr + p0_rr) / rho_rr
    H_roe = (sqrt_rho_ll * H_ll + sqrt_rho_rr * H_rr) / sum_sqrt_rho
    c_roe = sqrt((equations.γ - 1) * (H_roe - ekin_roe))

    # Signal speeds (Einfeldt)
    Ssl = min(v1_ll - c_ll, vel_roe - c_roe)
    Ssr = max(v1_rr + c_rr, vel_roe + c_roe)

    # Mass component of the HLLC flux
    if Ssl >= 0
        f1 = f1_ll
    elseif Ssr <= 0
        f1 = f1_rr
    else
        sMu_L = Ssl - v1_ll
        sMu_R = Ssr - v1_rr
        # speed of the contact wave
        SStar =
            (p0_rr - p0_ll + rho_ll * v1_ll * sMu_L - rho_rr * v1_rr * sMu_R) /
            (rho_ll * sMu_L - rho_rr * sMu_R)
        if Ssl <= 0 <= SStar
            # density in the left star region
            densStar = rho_ll * sMu_L / (Ssl - SStar)
            f1 = f1_ll + Ssl * (densStar - rho_ll)
        else
            # density in the right star region
            densStar = rho_rr * sMu_R / (Ssr - SStar)
            f1 = f1_rr + Ssr * (densStar - rho_rr)
        end
    end
    return SVector(f1, 0.0, 0.0, 0.0)
end

function Trixi.flux_hllc(
    u_ll,
    u_rr,
    normal_direction::AbstractVector,
    equations::TermiteMoundEquations1D,
)
    norm_ = abs(normal_direction[1])
    normal_direction_unit = normal_direction[1] * inv(norm_)
    f = Trixi.flux_hllc(
        SVector(u_ll[1], normal_direction_unit * u_ll[2], u_ll[3], u_ll[4]),
        SVector(u_rr[1], normal_direction_unit * u_rr[2], u_rr[3], u_rr[4]),
        1,
        equations,
    )
    return SVector(f[1], 0.0, 0.0, 0.0) * norm_
end

@doc raw"""
    Trixi.max_abs_speeds(u, equations::TermiteMoundEquations1D)

Maximum absolute wave speed at a single point, ``|u| + c`` with the sound speed
``c = \sqrt{\gamma p_0/\rho}``. Used by the `StepsizeCallback` for the CFL condition. The
`v1` field must therefore contain the current velocity ``u`` (see
[`AlgebraicVelocityCallback`](@ref) for the "coupled" method).

### Parameters
- `u`: state vector `(rho, v1, p, Ti)`.
- `equations`: the [`TermiteMoundEquations1D`](@ref).

### Returns
`(abs(v1) + c,)` (a 1-tuple in 1D).
"""
@inline function Trixi.max_abs_speeds(u, equations::TermiteMoundEquations1D)
    rho, v1, p0, _ = u
    c = sqrt(equations.γ * p0 / rho)
    return (abs(v1) + c,)
end

@doc raw"""
    Trixi.max_abs_speed_naive(u_ll, u_rr, orientation::Integer, equations::TermiteMoundEquations1D)

Estimate of the maximum wave speed between two states for local Lax-Friedrichs-type
dissipation:
```math
\lambda_{max} = \max(|u_{ll}|, |u_{rr}|) + \max(c_{ll}, c_{rr}), \qquad
c = \sqrt{\gamma p_0/\rho}.
```

### Parameters
- `u_ll`, `u_rr`: left and right state vectors `(rho, v1, p, Ti)`.
- `orientation`: orientation (always `1` in 1D).
- `equations`: the [`TermiteMoundEquations1D`](@ref).

### Returns
The wave speed estimate ``\lambda_{max}``.
"""
@inline function Trixi.max_abs_speed_naive(
    u_ll,
    u_rr,
    orientation::Integer,
    equations::TermiteMoundEquations1D,
)
    rho_ll, v1_ll, p0_ll, _ = u_ll
    rho_rr, v1_rr, p0_rr, _ = u_rr
    v_mag_ll = abs(v1_ll)
    v_mag_rr = abs(v1_rr)
    c_ll = sqrt(equations.γ * p0_ll / rho_ll)
    c_rr = sqrt(equations.γ * p0_rr / rho_rr)
    return max(v_mag_ll, v_mag_rr) + max(c_ll, c_rr)
end

# Parameters of the initial Ti profiles on the grid r ∈ TI_RADII, h ∈ TI_HEIGHTS (r-major):
# Ti(xa), Ti(xb) (maximum), Ti(xc) in °C and the curvature k of the chimney branch.
# Least-squares fits to quasi-periodic states of multi-day simulations.
const TI_RADII = range(0.3, 1.5, length = 5)
const TI_HEIGHTS = range(0.5, 2.5, length = 5)
const TI_PARAMS = [
    28.5629 29.1405 28.2099 2.70  # r = 0.3, h = 0.5
    27.9362 29.2659 27.7738 2.12  # r = 0.3, h = 1.0
    27.7044 29.3670 27.6229 2.00  # r = 0.3, h = 1.5
    27.5723 29.4020 27.5393 1.79  # r = 0.3, h = 2.0
    27.4933 29.4153 27.5071 1.72  # r = 0.3, h = 2.5
    28.3231 28.6427 28.1880 0.21  # r = 0.6, h = 0.5
    28.1104 28.8462 27.7681 1.56  # r = 0.6, h = 1.0
    27.8786 28.9473 27.6172 1.44  # r = 0.6, h = 1.5
    27.7464 28.9823 27.5336 1.23  # r = 0.6, h = 2.0
    27.6675 28.9956 27.5014 1.16  # r = 0.6, h = 2.5
    28.2534 28.7113 27.9918 0.09  # r = 0.9, h = 0.5
    28.1411 28.6864 27.8398 0.20  # r = 0.9, h = 1.0
    28.0170 28.8476 27.6526 0.85  # r = 0.9, h = 1.5
    27.9209 28.8660 27.5965 0.77  # r = 0.9, h = 2.0
    27.8420 28.8793 27.5643 0.70  # r = 0.9, h = 2.5
    28.2178 28.7448 27.9002 0.05  # r = 1.2, h = 0.5
    28.1478 28.7332 27.7618 0.13  # r = 1.2, h = 1.0
    28.0893 28.8205 27.6681 0.49  # r = 1.2, h = 1.5
    28.0065 28.8344 27.6180 0.46  # r = 1.2, h = 2.0
    27.9332 28.8399 27.5842 0.46  # r = 1.2, h = 2.5
    28.1954 28.7674 27.8552 0.05  # r = 1.5, h = 0.5
    28.1575 28.7612 27.7262 0.09  # r = 1.5, h = 1.0
    27.8388 28.5208 27.3742 0.28  # r = 1.5, h = 1.5
    28.0540 28.8197 27.6209 0.30  # r = 1.5, h = 2.0
    27.9961 28.8307 27.5987 0.32  # r = 1.5, h = 2.5
]

@doc raw"""
    Ti_profile(x, Ta, Tb, Tc, k, xa, xb, xc)

Initial mound material temperature (°C) at position `x`:
```math
T_i(x) = \begin{cases}
  T_a & x \le x_a, \\
  T_a + (T_b - T_a) \dfrac{x - x_a}{x_b - x_a} & x_a < x \le x_b \quad \text{(flute)}, \\[1ex]
  T_c + (T_b - T_c) \dfrac{e^{-k\eta} - e^{-k}}{1 - e^{-k}},
    \quad \eta = \dfrac{x - x_b}{x_c - x_b} & x_b < x < x_c \quad \text{(chimney)}, \\[1ex]
  T_c & x \ge x_c.
\end{cases}
```
The profile rises linearly from `Ta` at `xa` to the maximum `Tb` at the top of the mound `xb`,
then decays to `Tc` at `xc` with curvature `k`. The decay is linear for `k < 1e-8`. Outside of
``[x_a, x_c]``, the profile is constant; there ``T_i`` does not enter the model.
"""
function Ti_profile(x, Ta, Tb, Tc, k, xa, xb, xc)
    if x <= xa
        return Ta
    elseif x <= xb
        return Ta + (Tb - Ta) * (x - xa) / (xb - xa)
    elseif x < xc
        η = (x - xb) / (xc - xb)
        g = k < 1e-8 ? 1 - η : (exp(-k * η) - exp(-k)) / (1 - exp(-k))
        return Tc + (Tb - Tc) * g
    else
        return Tc
    end
end

@doc raw"""
    Get_initial_Ti(r, h; n_samples = 513)

Initial profile of the mound material temperature (°C) of a mound with radius `r` and height
`h` (both in m). The steps are:
1. The four profile parameters `(Ta, Tb, Tc, k)` of [`Ti_profile`](@ref) are bilinearly
   interpolated in `(r, h)` from the table `TI_PARAMS` on the grid
   ``r \in \{0.3, 0.6, \dots, 1.5\}``, ``h \in \{0.5, 1.0, \dots, 2.5\}``. Outside of the grid,
   `r` and `h` are clamped to it.
2. The positions `xa`, `xb`, `xc` follow from [`mound_geometry`](@ref).
3. The profile is sampled at `n_samples` equidistant points in ``[0, 1]`` plus `xa`, `xb`,
   `xc` (so that the kinks are represented exactly).

The table contains least-squares fits to quasi-periodic states of multi-day simulations.

### Parameters
- `r`: radius of the mound (m).
- `h`: height of the mound (m).
- `n_samples`: number of equidistant sample points (default: 513).

### Returns
A linear interpolation object on ``[0, 1]`` (constant extrapolation), stored as `Ti_LI` in
[`TermiteMoundEquations1D`](@ref).
"""
function Get_initial_Ti(r, h; n_samples = 513)
    rc = clamp(r, first(TI_RADII), last(TI_RADII))
    hc = clamp(h, first(TI_HEIGHTS), last(TI_HEIGHTS))
    i = min(searchsortedlast(TI_RADII, rc), length(TI_RADII) - 1)
    j = min(searchsortedlast(TI_HEIGHTS, hc), length(TI_HEIGHTS) - 1)
    s = (rc - TI_RADII[i]) / (TI_RADII[i+1] - TI_RADII[i])
    t = (hc - TI_HEIGHTS[j]) / (TI_HEIGHTS[j+1] - TI_HEIGHTS[j])
    row(i, j) = TI_PARAMS[(i-1)*length(TI_HEIGHTS)+j, :]
    Ta, Tb, Tc, k =
        (1 - s) * (1 - t) * row(i, j) +
        s * (1 - t) * row(i + 1, j) +
        (1 - s) * t * row(i, j + 1) +
        s * t * row(i + 1, j + 1)
    _, xa, xb, xc = mound_geometry(r, h)
    x = sort!(unique!(vcat(collect(range(0.0, 1.0, length = n_samples)), xa, xb, xc)))
    return linear_interpolation(
        x,
        Ti_profile.(x, Ta, Tb, Tc, k, xa, xb, xc),
        extrapolation_bc = Interpolations.Flat(),
    )
end

end # @muladd
