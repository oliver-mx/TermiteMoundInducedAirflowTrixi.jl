# Mathematics

The **1D termite mound model** originates from the Euler equations of gas dynamics in a low
Mach number regime. This page summarizes the derivation. More details can be found in the
corresponding publication.

## Scaled Euler equations

A scaled version of the Euler equations in one spatial dimension with the ideal gas law reads
```math
\begin{aligned}
\frac{\partial \rho}{\partial t} + \frac{\partial (\rho u)}{\partial x}
  &= -\frac{A_x}{A} \rho u, \\
\frac{\partial (\rho u)}{\partial t}
  + \frac{\partial}{\partial x} \left[ \rho u^2 + \frac{1}{\epsilon} p \right]
  &= -\frac{A_x}{A} \rho u^2
     - \frac{\beta}{A} \rho u \big(\eta + (1-\eta)|u|\big)
     - \frac{1}{Fr^2} \rho h_x, \\
\frac{\partial}{\partial t} \left[ \rho T
  + \epsilon (\gamma - 1) \left( \frac{\rho u^2}{2} + \frac{\rho h}{Fr^2} \right) \right] \\
  + \frac{\partial}{\partial x} \left[ \rho u T + (\gamma - 1) \left( \epsilon
    \left( \frac{\rho u^2}{2} + \frac{\rho u h}{Fr^2} \right) + u p \right) \right]
  &= -\frac{A_x}{A} \rho u T - \frac{k_w}{\sqrt{A}} (T - T_u) \\
  &\hphantom{=} - (\gamma - 1) \frac{A_x}{A}
     \left( \epsilon \left( \frac{\rho u^2}{2} + \frac{\rho u h}{Fr^2} \right) + u p \right), \\
p &= \rho T.
\end{aligned}
```

**Unknowns:**
- density ``\rho``,
- velocity ``u``,
- pressure ``p``,
- temperature ``T``.

**Variables:**
- space ``x``,
- time ``t``.

**Functions (model data):**
- cross section of the flow channel ``A``,
- height profile ``h``,
- boundary temperature ``T_u``.

**Parameters:**
- adiabatic exponent ``\gamma = \frac{c_v + R}{c_v}``,
- Mach number ``M = \sqrt{\frac{\rho_r u_r^2}{\gamma p_r}}``,
- Froude number ``Fr = \frac{u_r}{\sqrt{g L}}``,
- ``\epsilon = \gamma M^2 = \frac{\rho_r u_r^2}{p_r}``,
- pipe friction coefficient ``\lambda_w = \frac{64}{Rey}``,
- Reynolds number ``Rey = \frac{u D_0 \rho}{\eta}``,
- viscosity ``\eta``,
- friction coefficient ``\beta = \frac{\lambda_w x_r \sqrt{\pi}}{4 \sqrt{A_r}}``,
- wall heat transfer coefficient
  ``k_w = (\gamma - 1) \frac{x_r \alpha_w T_r \sqrt{\pi}}{u_r p_r \sqrt{A_r}}``.

The boundary temperature is defined piecewise by
```math
T_u(t, x) = \begin{cases}
  T_{soil}(t) & x \in [0, x_a] \cup (x_c, 1] \quad \text{(soil)}, \\
  \frac{1}{2} \left( T_{air}(t) + T_i(t, x) \right) & x \in (x_a, x_b] \quad \text{(flute)}, \\
  T_i(t, x) & x \in (x_b, x_c] \quad \text{(chimney)}.
\end{cases}
```
It depends on the mound material temperature ``T_i``, the soil temperature ``T_{soil}`` (at
30 cm depth) and the ambient air temperature ``T_{air}``. The functions ``T_{air}`` and
``T_{soil}`` (daily cycles) and the positions ``x_a``, ``x_b`` and ``x_c`` are known.
The mound material temperature ``T_i`` is an additional unknown. Its evolution is described by
```math
\frac{\partial T_i}{\partial t} = k_i (T - T_i).
```

**New unknown:**
- mound material temperature ``T_i``.

**New parameter:**
- heat exchange rate ``k_i = \frac{t_r \kappa L_i}{c_i m}`` with the thermal conductivity
  ``\kappa``, the specific heat capacity ``c_i`` and the mass ``m`` of the mound material,
  and the length ``L_i`` of the channel inside the mound.

## Asymptotic model

In the scaled equations, the parameter ``\epsilon`` is problematic. For a termite mound,
``\epsilon`` is small and the equations become stiff and expensive to solve. Therefore, we
consider the asymptotic expansion with leading order ``O(\epsilon^0)``
```math
\begin{aligned}
\rho &= \rho_0 + O(\epsilon), \\
u &= u_0 + O(\epsilon), \\
T &= T_0 + O(\epsilon), \\
p &= p_0 + p_h + \epsilon p_1 + O(\epsilon^2).
\end{aligned}
```
Inserting this expansion gives
```math
\begin{aligned}
\frac{\partial \rho}{\partial t} + \frac{\partial (\rho u)}{\partial x}
  &= -\frac{A_x}{A} \rho u, \\
\frac{\partial (\rho u)}{\partial t}
  + \frac{\partial}{\partial x} \left[ \rho u^2 + p_1 \right]
  &= -\frac{A_x}{A} \rho u^2
     - \frac{\beta}{A} \rho u \big(\eta + (1-\eta)|u|\big) \\
  &\hphantom{=} - \frac{h_x}{Fr^2} (\rho - \rho_{h_0}), \\
\frac{\partial p_0}{\partial t} + \gamma p_0 \frac{\partial u}{\partial x}
  &= -\gamma \frac{A_x}{A} u p_0 - \frac{k_w}{\sqrt{A}} (T - T_u), \\
\frac{\partial p_0}{\partial x} &= 0, \\
p_0 &= \rho T, \\
\frac{\partial T_i}{\partial t} &= k_i (T - T_i).
\end{aligned}
```

**Unknowns:**
- density ``\rho``,
- velocity ``u``,
- leading-order pressure ``p_0``,
- first-order pressure ``p_1``,
- air temperature ``T``,
- mound material temperature ``T_i``.

Due to the closed flow channel, we impose periodic boundary conditions. As a result,
``p_0`` is constant in space. We assume appropriate initial states ``\rho^0(x)``,
``u^0(x)``, ``p_0^0`` and ``T_i^0(x)``. The boundary conditions are ``p_1(t, 0) = 0``,
``p_1(t, 1) = 0`` and ``A(0) \rho(t, 0) u(t, 0) = A(1) \rho(t, 1) u(t, 1)``.

## Stationary case

In the stationary case, we make two assumptions for the mound material temperature. First,
``T_i`` is constant in ``x``. Second, the ordinary differential equation (ODE) for ``T_i`` is
replaced by an algebraic equation, assuming that the heat flux of the stationary system is in
equilibrium. The stationary equations are
```math
\begin{aligned}
\frac{\partial (\rho u)}{\partial x} &= -\frac{A_x}{A} \rho u, \\
\frac{\partial}{\partial x} \left[ \rho u^2 + p_1 \right]
  &= -\frac{A_x}{A} \rho u^2
     - \frac{\beta}{A} \rho u \big(\eta + (1-\eta)|u|\big) \\
  &\hphantom{=} - \frac{h_x}{Fr^2} (\rho - \rho_{h_0}), \\
\gamma p_0 \frac{\partial u}{\partial x}
  &= -\gamma \frac{A_x}{A} u p_0 - \frac{k_w}{\sqrt{A}} (T - T_u), \\
\frac{\partial p_0}{\partial x} &= 0, \\
p_0 &= \rho T, \\
0 &= \int_0^1 \frac{k_w}{\gamma \sqrt{A}} (T - T_u) \, dx.
\end{aligned}
```
With appropriate initial and boundary conditions, ``p_0 = 1``.

**Unknowns:**
- density ``\rho``,
- velocity ``u``,
- first-order pressure ``p_1``,
- air temperature ``T``,
- mound material temperature ``T_i`` (a scalar).

The stationary system can be solved as a differential-algebraic equation (DAE).

## Reformulation for the numerical solution

The ``p_0`` equation of the asymptotic model can be written as an equation for the spatial
variation of the velocity,
```math
\frac{\partial u}{\partial x} = -\frac{A_x}{A} u + Q, \qquad
Q = \frac{1}{\gamma p_0} \left[ -\frac{\mathrm{d} p_0}{\mathrm{d} t}
    - \frac{k_w}{\sqrt{A}} (T - T_u) \right],
```
i.e. ``(A u)_x = A Q``. Integrating from ``0`` to ``x`` gives the **velocity constraint**
```math
u(t, x) = \frac{v(t)}{A(x)} + \frac{1}{A(x)} \int_0^x \frac{A(y)}{\gamma p_0}
  \left[ -\frac{\mathrm{d} p_0}{\mathrm{d} t} - \frac{k_w}{\sqrt{A}} (T - T_u) \right] dy,
```
with the new variable ``v(t) = A(0) u(t, 0)``, the volume flux at ``x = 0``. Since ``A u`` is
periodic, the integral over the whole domain vanishes. This determines the time derivative of
the pressure:
```math
\frac{\mathrm{d} p_0}{\mathrm{d} t} = -\gamma p_0 \,
  \frac{\int_0^1 \frac{k_w \sqrt{A}}{\gamma p_0} (T - T_u) \, dx}{\int_0^1 A \, dx}.
```
The variable ``v`` satisfies the ODE
```math
\frac{\mathrm{d} v}{\mathrm{d} t} = \frac{1}{\int_0^1 \frac{\rho}{A} \, dy}
  \int_0^1 \left[ -\rho u u_x - \frac{\beta}{A} \rho u \big(\eta + (1-\eta)|u|\big)
  - \frac{h_x}{Fr^2} (\rho - \rho_{h_0}) \right] dy,
```
which is obtained by integrating the momentum balance of the asymptotic model over the
domain. The periodic boundary conditions of ``p_1`` eliminate the variable ``p_1``.

Finally, the reformulated asymptotic model, i.e. the **1D termite mound model**, reads
```math
\begin{aligned}
\frac{\partial \rho}{\partial t} + \frac{\partial (\rho u)}{\partial x}
  &= -\frac{A_x}{A} \rho u, \\
\frac{\partial T_i}{\partial t} &= k_i (T - T_i), \\
\frac{\partial p_0}{\partial t} + \gamma p_0 \frac{\partial u}{\partial x}
  &= -\gamma \frac{A_x}{A} u p_0 - \frac{k_w}{\sqrt{A}} (T - T_u), \\
\frac{\mathrm{d} v}{\mathrm{d} t}
  &= \frac{1}{\int_0^1 \frac{\rho}{A} \, dy} \int_0^1 \Big[ -\rho u u_x
     - \frac{\beta}{A} \rho u \big(\eta + (1-\eta)|u|\big) \\
  &\hphantom{= \frac{1}{\int_0^1 \frac{\rho}{A} \, dy} \int_0^1 \Big[}
     - \frac{h_x}{Fr^2} (\rho - \rho_{h_0}) \Big] \, dy, \\
p_0 &= \rho T, \\
u(t, x) &= \frac{v(t)}{A(x)} + \frac{1}{A(x)} \int_0^x \frac{A(y)}{\gamma p_0}
  \left[ -\frac{\mathrm{d} p_0}{\mathrm{d} t} - \frac{k_w}{\sqrt{A}} (T - T_u) \right] dy.
\end{aligned}
```

**Unknowns:**
- density ``\rho``,
- velocity ``u``,
- auxiliary velocity ``v``,
- leading-order pressure ``p_0``,
- air temperature ``T``,
- mound material temperature ``T_i``.
