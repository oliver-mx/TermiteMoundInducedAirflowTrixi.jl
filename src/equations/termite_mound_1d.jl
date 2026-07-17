# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent
@doc raw"""
    TermiteMoundEquations1D()

TermiteMoundEquations1D (TME) in one space dimension. The equations are given by
```math
	\begin{equation}
		\frac{\partial \rho}{\partial t} + \frac{\partial \left( \rho u\right)}{\partial x} = - \frac{\textrm{A}_x}{\textrm{A}}\rho u 
	\end{equation}
    \begin{equation}
		\frac{\partial v}{\partial t} = \frac{1}{\int_0^1 \frac{\rho}{A} \, dy} \left[ 
		\int_0^1 -\rho u u_x  - \rho u \left( \beta \eta - \beta\left(1-\eta\right) \vert u \vert \right) - \frac{\textrm{h}_x}{Fr^2}(\rho - \rho_{h_0}) \, dy \right] 
	\end{equation}
    \begin{equation}
		\frac{\partial p_0}{\partial t} + \gamma p_0 \frac{\partial u}{\partial x}= - \gamma \frac{\textrm{A}_x}{\textrm{A}} u p_0 - \frac{k_w}{\textrm{A} \sqrt{\textrm{A}}} \left(T - \textrm{T}_\textrm{u} \right) 
	\end{equation}
    \begin{equation}
		\frac{\partial T_i}{\partial t} = k_i \left( T - T_i \right) 
	\end{equation}
    \begin{equation}
		p_0 = \rho T 
	\end{equation}
    \begin{equation}
		u(t,x) = \frac{v(t)}{A(x)} + \frac{1}{A(x)} \int_0^x A(y)\left[ -\frac{\partial  p_0}{\partial t} - \frac{k_w}{\textrm{A} \sqrt{\textrm{A}}} \left(T - \textrm{T}_\textrm{u} \right)\right] \, dy 
	\end{equation}
```
The unknown quantities of the TME are the air density ``ρ``, the air velocity ``u``, the leading order pressure ``p_0``, the internal temperature of the mound ``T_i``, and the auxiliary velocity variable ``v``.
Due to the single spatial dimensional, the velocity ``u`` is referred to as ``v1`` in the implementation.
We denote the (possibly) variable cross section function ``A(x)``, the height profile function ``h(x)`` and the boundary temperature function ``T_u(t,x)``.
The boundary temperature function is defined by
```math
    \begin{align}
		\textrm{T}_\textrm{u}(t,x)  = \begin{cases}
			\textrm{T}_\textrm{soil}(t) & ,\text{if } x \in [0,x_a] \, \cup \, (x_c, 1] \\
			\frac{\textrm{T}_\textrm{air}(t) + T_i(t,x) }{2} & ,\text{if } x \in (x_a,x_b] \\
			T_i(t,x) & ,\text{if } x \in (x_b,x_c] \\
		\end{cases}
	\end{align}
```
It depends on the internal mound temperature ``T_i``, the soil temperature (at 30cm depth) ``T_{soil}`` and the ambient air temperature ``T_{air}``. 
The parameters ``x_a``, ``x_b`` and ``x_c`` are fixed spatial locations along the flow channel.

The viscosity parameter is denoted by ``η`` and the dimensionless parameters are ``k_w`` and ``k_i``.

Reference for the TME:
- Marx, Oliver P. and Gasser, Ingenuin and Annika Schmidgall (2026)
  ...
  [DOI: ...](https://doi.org/...)
"""

struct TermiteMoundEquations1D{RealT<:Real} <: AbstractEquations{1,5}
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
    function TermiteMoundEquations1D(;
        γ=1.3987529976019184,
        k_i=0.0004184703194089407,
        k_w=4.550553359854458,
        tᵣ=87.10175873820496,
        uᵣ=0.05,
        xa=0.17221939177753745,
        xb=0.5680821161575328,
        xc=0.9793345160180965,
        r=0.6874197366854203,
        h=1.5810403656581307,
        L=4.355087936910248,
        β=26.9749915565812,
        η=0.93,
        Fr²=5.851592474205325e-5,
        ρₕ₀=0.982,
        T_ref=297.76980029566033,
        t_ref=0.0033582989242263127,
        T0=27.0,
        v0=-2.219,
        Ti_LI=LinearInterpolation([0.0, 0.5, 1.0], [27.6, 28.8, 27.5]),
    )
        γ, inv_gamma_minus_one = promote(γ, inv(γ - 1.0))
        new{typeof(β)}(
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
        )
    end
end

function Trixi.varnames(u, ::TermiteMoundEquations1D)
    ("rho", "v1", "p0", "Ti", "x_var")
end

function Trixi.cons2prim(u, ::TermiteMoundEquations1D)
    return u
end

function Trixi.prim2cons(prim, ::TermiteMoundEquations1D)
    return prim
end

function Trixi.cons2entropy(u, ::TermiteMoundEquations1D)
    return u
end

function Trixi.cons2cons(u, ::TermiteMoundEquations1D)
    return u
end

@inline function TermiteMoundInitialCondition(x, t, equations::TermiteMoundEquations1D)
    RealT = eltype(x)
    rho = inv(temp2scaled(equations.T0, equations))
    v1 = vel2scaled(equations.v0, equations)
    p0 = 1.0
    Ti_LI = equations.Ti_LI
    Ti = temp2scaled(Ti_LI(x[1]), equations)
    return SVector(rho, v1, p0, Ti, x[1])
end

@inline function source_terms(u, x, t, equations::TermiteMoundEquations1D)
    rho, v1, p0, Ti, _ = u
    x_var = x[1]
    T = p0 / rho
    du1 = - A_x(x_var, equations) / A(x_var, equations) * rho * v1
    du4 = equations.k_i * (T - Ti)
    return SVector(du1, 0.0, 0.0, du4, 0.0)
end

@inline function temp2scaled(y, equations::TermiteMoundEquations1D)
    z = 273.15
    y = y .+ z
    x = y ./ equations.T_ref
    return x
end

@inline function vel2scaled(y, equations::TermiteMoundEquations1D)
    return inv(equations.uᵣ * 100) .* y
end

@inline function time2scaled(y, equations::TermiteMoundEquations1D)
    return y * (3600/equations.tᵣ)
end

@inline function space2unscaled(x, equations::TermiteMoundEquations1D)
    return equations.L .* x
end

@inline function temp2unscaled(x, equations::TermiteMoundEquations1D)
    z = 273.15
    return x .* equations.T_ref .- z
end

@inline function vel2unscaled(x, equations::TermiteMoundEquations1D)
    return x .* equations.uᵣ * 100
end

@inline function time2unscaled(x, equations::TermiteMoundEquations1D)
    return x ./ (3600/equations.tᵣ)
end

@inline function h(x, equations::TermiteMoundEquations1D)
    #xb = (sqrt(0.09 + equations.r*r) + sqrt(r*r + h*h)) / L
    LI = LinearInterpolation(
        [0.0, equations.xa, equations.xb, 1.0],
        [0.0, 0.3, equations.h + 0.3, 0.0],
    )
    return LI(x)
end

@inline function h_x(x, equations::TermiteMoundEquations1D)
    dy = if x < equations.xa
        0.3 / (equations.xa * equations.L)
    elseif x < equations.xb
        equations.h / (equations.xb - equations.xa) / equations.L
    else
        -1.0
    end
    return dy
end

@inline function A(x, ::TermiteMoundEquations1D)
    return 1.0
end

@inline function A_x(x, ::TermiteMoundEquations1D)
    return 0.0
end

@inline function T_u(t_var, x_var, Ti, equations::TermiteMoundEquations1D)
    t_var = mod(t_var, 86400 / equations.tᵣ)
    t_var = t_var/(86400 / equations.tᵣ/24)
    c0 = T_air(t_var, equations)
    c1 = T_soil(t_var, equations)
    c2 = 0.5f0 * (c0 + Ti * equations.T_ref)
    c3 = Ti * equations.T_ref
    return ifelse(
        (x_var ≤ equations.xa) || (x_var ≥ equations.xc),
        c1,
        ifelse(x_var ≤ equations.xb, c2, c3),
    )
end

@inline function T_u_dt(t_var, x_var, Ti_dt, equations::TermiteMoundEquations1D)
    t_var = mod(t_var, 86400 / equations.tᵣ)
    t_var = t_var/(86400 / equations.tᵣ / 24)
    c0_dt = dt_T_air(t_var, equations)
    c1_dt = dt_T_soil(t_var, equations)
    c2_dt = 0.5f0 * (c0_dt + Ti_dt)
    c3_dt = Ti_dt
    return ifelse(
        (x_var ≤ equations.xa) || (x_var ≥ equations.xc),
        c1_dt,
        ifelse(x_var ≤ equations.xb, c2_dt, c3_dt),
    )
end

@inline function LinearInterpolation2(x, y, ::TermiteMoundEquations1D)
    iunique_indices = unique!(collect(zip(x, y)))[2]
    x_unique = []
    y_unique = []
    seen = Dict{typeof(x[1]),Bool}()
    for (xi, yi) in zip(x, y)
        if haskey(seen, xi)
            continue
        else
            push!(x_unique, xi)
            push!(y_unique, yi)
            seen[xi] = true
        end
    end
    LI = LinearInterpolation(x_unique, y_unique)
    return LI
end

@inline function Fourir(
    a0,
    a1,
    a2,
    a3,
    a4,
    a5,
    a6,
    a7,
    a8,
    a9,
    a10,
    t_var,
    ::TermiteMoundEquations1D,
)
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

@inline function dt_Fourir(
    a1,
    a2,
    a3,
    a4,
    a5,
    a6,
    a7,
    a8,
    a9,
    a10,
    t_var,
    ::TermiteMoundEquations1D,
)
    return - a1 .* sin(pi * t_var / 12.0) .* (pi / 12.0) +
           a2 * cos(pi * t_var / 12.0) * (pi / 12.0) -
           a3 * sin(pi * t_var / 6.0) * (pi / 6.0) +
           a4 * cos(pi * t_var / 6.0) * (pi / 6.0) -
           a5 * sin(pi * t_var / 4.0) * (pi / 4.0) +
           a6 * cos(pi * t_var / 4.0) * (pi / 4.0) -
           a7 * sin(pi * t_var / 3.0) * (pi / 3.0) +
           a8 * cos(pi * t_var / 3.0) * (pi / 3.0) -
           a9 * sin(pi * t_var / 2.4) * (pi / 2.4) +
           a10 * cos(pi * t_var / 2.4) .* (pi / 2.4)
end

@inline function bspline2linear(nodes, vals, t, ti, ::TermiteMoundEquations1D)
    itp = Interpolations.scale(
        interpolate(hcat(nodes, vals), (BSpline(Cubic(Natural(OnGrid()))), NoInterp())),
        t,
        1:2,
    )
    nodesitp, valsitp = [itp(t, 1) for t in ti], [itp(t, 2) for t in ti]
    return LinearInterpolation(nodesitp, valsitp)
end

@inline function T_air(t_var, equations::TermiteMoundEquations1D)
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
    return Fourir(
        a0,
        a1,
        a2,
        a3,
        a4,
        a5,
        a6,
        a7,
        a8,
        a9,
        a10,
        t_var,
        equations::TermiteMoundEquations1D,
    )
end

@inline function dt_T_air(t_var, equations::TermiteMoundEquations1D)
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
    return dt_Fourir(
        a1,
        a2,
        a3,
        a4,
        a5,
        a6,
        a7,
        a8,
        a9,
        a10,
        t_var,
        equations::TermiteMoundEquations1D,
    )
end

@inline function T_soil(t_var, equations::TermiteMoundEquations1D)
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
    return Fourir(
        a0,
        a1,
        a2,
        a3,
        a4,
        a5,
        a6,
        a7,
        a8,
        a9,
        a10,
        t_var,
        equations::TermiteMoundEquations1D,
    )
end

@inline function dt_T_soil(t_var, equations::TermiteMoundEquations1D)
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
    return dt_Fourir(
        a1,
        a2,
        a3,
        a4,
        a5,
        a6,
        a7,
        a8,
        a9,
        a10,
        t_var,
        equations::TermiteMoundEquations1D,
    )
end

@inline function I_w(x_var, equations::TermiteMoundEquations1D)
    if x_var ≤ equations.xa
        return 0.5
    elseif x_var ≤ equations.xc
        return 1.0
    else
        return 0.5
    end
end

@inline function unique_idx(x, ::TermiteMoundEquations1D)
    y=1:1:length(x)
    iunique_indices = unique!(collect(zip(x, y)))[2]
    x_unique = Vector{typeof(x[1])}()
    y_unique = Vector{Int}()
    seen = Dict{typeof(x[1]),Bool}()
    for (xi, yi) in zip(x, y)
        if haskey(seen, xi)
            continue
        else
            push!(x_unique, xi)
            push!(y_unique, yi)
            seen[xi] = true
        end
    end
    return x_unique, y_unique
end

# Calculate 1D flux for a single point
@inline function Trixi.flux(u, orientation::Integer, ::TermiteMoundEquations1D)
    rho, v1, _, _, _ = u
    f1 = rho * v1
    return SVector(f1, 0.0, 0.0, 0.0, 0.0)
end

@inline function flux_ranocha(
    u_ll,
    u_rr,
    orientation::Integer,
    ::TermiteMoundEquations1D,
)
    rho_ll, v1_ll, _, _, _ = u_ll
    rho_rr, v1_rr, _, _, _ = u_rr
    rho_mean = ln_mean(rho_ll, rho_rr)
    v1_mean = if v1_ll > 0.0 && v1_rr > 0.0
        ln_mean(v1_ll, v1_rr)
    elseif v1_ll < 0.0 && v1_rr < 0.0
        - ln_mean(-v1_ll, -v1_rr)
    else
        0.5f0 * (v1_ll + v1_rr)
    end
    f1 = rho_mean * v1_mean
    return SVector(f1, 0.0, 0.0, 0.0, 0.0)
end

@inline function Trixi.max_abs_speeds(u, equations::TermiteMoundEquations1D)
    rho, v1, p0, _, _ = u
    c = sqrt(equations.γ * p0 / rho)
    return (abs(v1) + c,)
end
end # @muladd
