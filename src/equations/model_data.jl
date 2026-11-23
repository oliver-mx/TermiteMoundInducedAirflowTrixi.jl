# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Model data of the TME that are passed to `TermiteMoundEquations1D` in the elixir:
# - cross section A(x) and its derivative A_x(x)           (`cross_section = ...`)
# - slope of the height profile h_x(x)                     (`height_profile = ...`)
# - boundary (wall) temperature T_u(t, x, Ti)              (`boundary_temperature = ...`)
# For each, there is a predefined termite mound version and a generic wrapper of user-defined
# functions. New types only need to implement the functions listed in the docstrings.

####################################################################################################
# Cross section

@doc raw"""
    CrossSectionConstant()

Constant cross section ``A \equiv 1``, ``A_x \equiv 0`` (the default of the termite mound
simulations). Pass it as `cross_section` to [`TermiteMoundEquations1D`](@ref).

See also [`CrossSectionSine`](@ref) and [`CrossSection`](@ref).
"""
struct CrossSectionConstant end

@inline A(x, ::CrossSectionConstant) = 1.0
@inline A_x(x, ::CrossSectionConstant) = 0.0

@doc raw"""
    CrossSectionSine(amplitude)

Cross section
```math
A(x) = 1 + a \sin(2\pi x + \pi/2), \qquad A_x(x) = 2\pi a \cos(2\pi x + \pi/2),
```
with ``a`` = `amplitude`. ``A`` is ``1``-periodic and largest at ``x = 0``; ``|a| < 1`` is
required for ``A > 0``. For ``a = 0``, it returns exactly ``A = 1`` and ``A_x = 0``.
Pass it as `cross_section` to [`TermiteMoundEquations1D`](@ref).
"""
struct CrossSectionSine{RealT<:Real}
    amplitude::RealT
    function CrossSectionSine(amplitude::RealT) where {RealT<:Real}
        abs(amplitude) < 1 ||
            throw(ArgumentError("|amplitude| < 1 is required for A > 0"))
        return new{RealT}(amplitude)
    end
end

@inline function A(x, cross_section::CrossSectionSine)
    a = cross_section.amplitude
    iszero(a) && return one(a)
    return 1 + a * sin(2 * pi * x + 0.5f0 * pi)
end

@inline function A_x(x, cross_section::CrossSectionSine)
    a = cross_section.amplitude
    iszero(a) && return zero(a)
    return 2 * pi * a * cos(2 * pi * x + 0.5f0 * pi)
end

@doc raw"""
    CrossSection(A, A_x)

User-defined cross section from two functions `A(x)` and `A_x(x)` of the scalar position
``x \in [0, 1]``. Pass it as `cross_section` to [`TermiteMoundEquations1D`](@ref), e.g.
```julia
cross_section = CrossSection(x -> 1 + 0.2 * cos(2pi * x), x -> -0.4pi * sin(2pi * x))
```
Requirements (not checked):
- ``A > 0`` on ``[0, 1]``,
- ``A`` and ``A_x`` are ``1``-periodic (the domain is periodic) and smooth,
- `A_x` is the exact derivative of `A`.
"""
struct CrossSection{AFunction,AxFunction}
    A::AFunction
    A_x::AxFunction
end

@inline A(x, cross_section::CrossSection) = cross_section.A(x)
@inline A_x(x, cross_section::CrossSection) = cross_section.A_x(x)

####################################################################################################
# Height profile

@doc raw"""
    HeightProfileTermiteMound(; xa, xb, L, h)

Slope ``h_x`` of the height profile of the flow channel of a termite mound (scaled). It is
piecewise constant:
```math
h_x(x) = \begin{cases}
  \dfrac{0.3}{x_a L} & x < x_a \quad \text{(soil, rising)}, \\[1ex]
  \dfrac{h}{(x_b - x_a) L} & x_a \le x < x_b \quad \text{(flute, rising to the top)}, \\[1ex]
  -1 & x \ge x_b \quad \text{(chimney and soil, descending)}.
\end{cases}
```
The arguments are the scaled transitions `xa`, `xb`, the channel length `L` (m) and the
height `h` (m) of the mound, as returned by [`termite_parameters`](@ref).
Pass it as `height_profile` to [`TermiteMoundEquations1D`](@ref).
"""
struct HeightProfileTermiteMound{RealT<:Real}
    xa::RealT
    xb::RealT
    L::RealT
    h::RealT
end

function HeightProfileTermiteMound(; xa, xb, L, h)
    return HeightProfileTermiteMound(promote(xa, xb, L, h)...)
end

@inline function h_x(x, height_profile::HeightProfileTermiteMound)
    @unpack xa, xb, L, h = height_profile
    dy = if x < xa
        0.3 / (xa * L)
    elseif x < xb
        h / (xb - xa) / L
    else
        -1.0
    end
    return dy
end

@doc raw"""
    HeightProfile(h_x)

User-defined slope of the height profile from a function `h_x(x)` of the scalar position
``x \in [0, 1]`` (scaled). Pass it as `height_profile` to [`TermiteMoundEquations1D`](@ref),
e.g.
```julia
height_profile = HeightProfile(x -> x < 0.5 ? 1.0 : -1.0)
```
"""
struct HeightProfile{HxFunction}
    h_x::HxFunction
end

@inline h_x(x, height_profile::HeightProfile) = height_profile.h_x(x)

####################################################################################################
# Boundary temperature
#
# Interface: `T_u(t, x, Ti, boundary_temperature)` returns the boundary temperature in K at the
# scaled time `t`, position `x` and scaled mound temperature `Ti`. For efficiency, the solvers
# split the evaluation into a time-only part, computed once for all nodes,
#   time_data = T_u_time(t, boundary_temperature)
#   T_u_node(time_data, x, Ti, boundary_temperature)
# The fallbacks below (time_data = t) work for every type that implements `T_u`.

@inline T_u_time(t, boundary_temperature) = t
@inline T_u_node(t, x, Ti, boundary_temperature) = T_u(t, x, Ti, boundary_temperature)

@doc raw"""
    BoundaryTemperatureTermiteMound(; xa, xb, xc, tᵣ, T_ref)

Boundary (wall) temperature ``T_u`` of the flow channel of a termite mound, defined piecewise:
```math
T_u(t, x) = \begin{cases}
  T_{soil}(t) & x \le x_a \text{ or } x \ge x_c \quad \text{(soil)}, \\
  \frac{1}{2}\left(T_{air}(t) + T_i(t, x)\right) & x_a < x \le x_b \quad \text{(flute)}, \\
  T_i(t, x) & x_b < x < x_c \quad \text{(chimney)}.
\end{cases}
```
It depends on the mound material temperature ``T_i``, the soil temperature ``T_{soil}`` (at
30 cm depth, [`T_soil`](@ref)) and the ambient air temperature ``T_{air}`` ([`T_air`](@ref)),
which are daily cycles. The arguments are the scaled transitions `xa`, `xb`, `xc`, the reference
time `tᵣ` (s) and the reference temperature `T_ref` (K), as returned by
[`termite_parameters`](@ref). Pass it as `boundary_temperature` to
[`TermiteMoundEquations1D`](@ref).

``T_{air}`` and ``T_{soil}`` only depend on time; the solvers compute them once per time point
([`T_air_soil`](@ref)) instead of once per node.
"""
struct BoundaryTemperatureTermiteMound{RealT<:Real}
    xa::RealT
    xb::RealT
    xc::RealT
    tᵣ::RealT
    T_ref::RealT
end

function BoundaryTemperatureTermiteMound(; xa, xb, xc, tᵣ, T_ref)
    return BoundaryTemperatureTermiteMound(promote(xa, xb, xc, tᵣ, T_ref)...)
end

@doc raw"""
    T_air_soil(t_var, boundary_temperature::BoundaryTemperatureTermiteMound)

Ambient air and soil temperature `(T_air, T_soil)` in K at the scaled time `t_var`. The time is
converted into the hour of the day: ``t \mapsto (t \bmod t_{day}) / (t_{day}/24)`` with the
scaled length of a day ``t_{day} = 86400/t_r``.
"""
@inline function T_air_soil(
    t_var,
    boundary_temperature::BoundaryTemperatureTermiteMound,
)
    @unpack tᵣ = boundary_temperature
    t_var = mod(t_var, 86400 / tᵣ)
    t_var = t_var/(86400 / tᵣ/24)
    return (T_air(t_var), T_soil(t_var))
end

@inline function T_u_time(t, boundary_temperature::BoundaryTemperatureTermiteMound)
    return T_air_soil(t, boundary_temperature)
end

# T_u from the precomputed (T_air, T_soil)
@inline function T_u_node(
    T_air_soil::Tuple,
    x,
    Ti,
    boundary_temperature::BoundaryTemperatureTermiteMound,
)
    @unpack xa, xb, xc, T_ref = boundary_temperature
    c0, c1 = T_air_soil
    c2 = 0.5f0 * (c0 + Ti * T_ref)
    c3 = Ti * T_ref
    return ifelse((x ≤ xa) || (x ≥ xc), c1, ifelse(x ≤ xb, c2, c3))
end

@inline function T_u(t, x, Ti, boundary_temperature::BoundaryTemperatureTermiteMound)
    return T_u_node(T_u_time(t, boundary_temperature), x, Ti, boundary_temperature)
end

@doc raw"""
    BoundaryTemperature(T_u)

User-defined boundary (wall) temperature from a function `T_u(t, x, Ti)` of the scaled time
`t`, the scalar position `x` and the scaled mound material temperature `Ti` at `x`. It must
return the temperature **in K** (it is scaled by `t_ref` ``= 1/T_{ref}`` in the solvers).
Pass it as `boundary_temperature` to [`TermiteMoundEquations1D`](@ref), e.g. a constant wall
temperature of 300 K:
```julia
boundary_temperature = BoundaryTemperature((t, x, Ti) -> 300.0)
```
"""
struct BoundaryTemperature{TuFunction}
    T_u::TuFunction
end

@inline T_u(t, x, Ti, boundary_temperature::BoundaryTemperature) =
    boundary_temperature.T_u(t, x, Ti)

####################################################################################################
# Daily cycles of the air and soil temperature

@doc raw"""
    Fourier(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, t_var)

Evaluates the truncated Fourier series with a period of 24 hours
```math
F(t) = a_0 + \sum_{i=1}^{5} \left[ a_{2i-1} \cos\left(\frac{i \pi t}{12}\right)
       + a_{2i} \sin\left(\frac{i \pi t}{12}\right) \right],
```
with the time ``t`` in hours. Used for the daily cycles of [`T_air`](@ref) and
[`T_soil`](@ref).

### Parameters
- `a0`, …, `a10`: Fourier coefficients ``a_0, \dots, a_{10}``.
- `t_var`: time of day in hours.

### Returns
The value ``F(t)``.
"""
@inline function Fourier(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, t_var)
    return a0 +
           a1 * cos(pi * t_var / 12.0) +
           a2 * sin(pi * t_var / 12.0) +
           a3 * cos(pi * t_var / 6.0) +
           a4 * sin(pi * t_var / 6.0) +
           a5 * cos(pi * t_var / 4.0) +
           a6 * sin(pi * t_var / 4.0) +
           a7 * cos(pi * t_var / 3.0) +
           a8 * sin(pi * t_var / 3.0) +
           a9 * cos(pi * t_var / 2.4) +
           a10 * sin(pi * t_var / 2.4)
end

@doc raw"""
    T_air(t_var)

Ambient air temperature ``T_{air}`` as a daily cycle, given by the Fourier series
[`Fourier`](@ref). The coefficients are a Fourier series fit of the mean air temperature in
Figure 2 of:
- Singh, R. K. and Sharma, R. V. (2017)
  Numerical analysis for ground temperature variation
  [DOI: 10.1186/s40517-017-0082-z](https://doi.org/10.1186/s40517-017-0082-z)

### Parameters
- `t_var`: time of day in hours (see [`T_air_soil`](@ref) for the conversion from the scaled
  time).

### Returns
The air temperature in K.
"""
@inline function T_air(t_var)
    a0 = 300.1;
    a1 = -2.627;
    a10 = -0.03607;
    a2 = -6.092;
    a3 = 0.4451;
    a4 = 0.2248;
    a5 = -0.5704;
    a6 = 0.2365;
    a7 = 0.02813;
    a8 = -0.4064;
    a9 = -0.07545
    return Fourier(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, t_var)
end

@doc raw"""
    T_soil(t_var)

Soil temperature ``T_{soil}`` at 30 cm depth as a daily cycle, given by the Fourier series
[`Fourier`](@ref). The coefficients are a Fourier series fit of the mean soil temperature at
30 cm depth in Figures 8 and 9 of:
- Singh, R. K. and Sharma, R. V. (2017)
  Numerical analysis for ground temperature variation
  [DOI: 10.1186/s40517-017-0082-z](https://doi.org/10.1186/s40517-017-0082-z)

### Parameters
- `t_var`: time of day in hours (see [`T_air_soil`](@ref) for the conversion from the scaled
  time).

### Returns
The soil temperature in K.
"""
@inline function T_soil(t_var)
    a0 = 302.9;
    a1 = 1.413;
    a10 = -0.005461;
    a2 = -0.09037;
    a3 = -0.08303;
    a4 = -0.2159;
    a5 = -0.01125;
    a6 = 0.01467;
    a7 = -0.0008037;
    a8 = -0.01079;
    a9 = -0.002245
    return Fourier(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, t_var)
end
end # @muladd
