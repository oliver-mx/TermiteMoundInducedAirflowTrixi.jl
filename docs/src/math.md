The **1D Termite Mound Model** originates from the Euler equations of gas dynamics in a low Mach number regime.

## Scaled Euler equations

A scaled version of the Euler equations with one spatial dimension and ideal gas law is given by 
```math
	\begin{aligned}
		\frac{\partial \rho}{\partial t} + \frac{\partial \left( \rho u\right)}{\partial x} &= - \frac{\textrm{A}_x}{\textrm{A}}\rho u \\
		\frac{\partial \left(\rho u\right)}{\partial t} + \frac{\partial}{\partial x} \left[\rho u^2 + \frac{1}{\epsilon}p\right] &= - \frac{\textrm{A}_x}{\textrm{A}}\rho u^2- \rho u \left( \beta \eta + \beta\left(1-\eta\right) \vert u \vert \right) - \frac{1}{Fr^2} \rho \textrm{h}_x \\
		\frac{\partial}{\partial t}\left[ \rho T + \epsilon \left(\gamma -1\right) \left(\frac{\rho u^2}{2} + \frac{1}{Fr^2} \rho \textrm{h} \right) \right] + \frac{\partial}{\partial x}\left[\rho u T + \left(\gamma -1 \right) \left( \epsilon \left( \frac{\rho u^2}{2} + \frac{1}{Fr^2} \rho u \textrm{h} \right) + up\right) \right] &= - \frac{\textrm{A}_x}{\textrm{A}}\rho u T - \left(\gamma -1\right) \frac{\textrm{A}_x}{\textrm{A}} \left( \epsilon\left(\frac{\rho u^2}{2} + \frac{1}{Fr^2} \rho u \textrm{h}\right) + up \right) - \frac{k_w}{\textrm{A} \sqrt{\textrm{A}}} \left(T - \textrm{T}_\textrm{u} \right) \\
		p &= \rho T.
    \end{aligned}
```

Unknowns:

- Density ``\rho``

- Velocity ``u``

- Pressure ``p``

- Temperature ``p``

Variables:

- Space ``x``

- Time ``t``

Functions:

- Flow channel cross section ``A``

- Height profile ``h``

- Boundary temperature function ``T_u``

Parameters:

-  Diabatic constant ``\gamma = \frac{c_v + R}{c_v}``

-  Mach number ``M= \sqrt{\frac{\rho_r u_{r}^2}{\gamma p_r}}``

-  Froude number ``Fr = \sqrt{\frac{u_r}{g L}}``

- ``\epsilon = \gamma M^2 = \frac{\rho_r u_{r}^2}{p_r}``

- Pipe friction coefficient ``\lambda_w = \frac{64}{Re}``

- Reynolds-number ``Re = \frac{u D_0 \rho}{\eta}``

- Viscosity ``eta``

- ``\beta = \frac{\lambda_w x_r \sqrt{\pi}}{u \sqrt{A_r A}}``

- ``k_w = \left(\gamma-1\right) \frac{x_r \alpha_w T_r \sqrt{\pi}}{u_r p_r \sqrt{A_r}}``

The boundary temperature function is defined by
```math
    \begin{align*}
		\textrm{T}_\textrm{u}(t,x)  &= \begin{cases}
			\textrm{T}_\textrm{soil}(t) & ,\text{if } x \in [0,x_a] \, \cup \, (x_c, 1]  \quad \text{(soil)} \\
			\frac{\textrm{T}_\textrm{air}(t) + T_i(t,x) }{2} & ,\text{if } x \in (x_a,x_b]  \quad \text{(flute)} \\
			T_i(t,x) & ,\text{if } x \in (x_b,x_c]  \quad \text{(chimney).} \\
		\end{cases}
	\end{align*}
```
It depends on the internal mound temperature ``T_i``, the soil temperature (at 30cm depth) ``T_{soil}`` and the ambient air temperature ``T_{air}``. 
The functions ``T_{air}`` and ``T_{soil}`` and the parameters ``x_a``, ``x_b`` and ``x_c`` are known.
The function ``T_i`` becomes an additional unknown to our system.
Its evolution is described by 
```math
    \begin{align*}
		\frac{\partial T_i}{\partial t} &= k_i \left( T - T_i \right).
	\end{align*}
```

New unknown:

- Termite mound material temperature ``T_i``

New parameters:

-  ``k_i = \frac{t_r \kappa T_r}{c_i M_i}``

## Asymptotic model

In the scaled equations the parameter ``\epsilon`` is problematic.
When modelling a termite mound ``\epsilon`` will be small and the model becomes expensive to solve.
Therefore, we consider the following asmyptotic expansion with leading order ``O(\epsilon^0)``

```math
    \begin{align*}
		\rho &= \rho_0 + O(\epsilon) \\
		u &= u_0 + O(\epsilon) \\
		T &= T_0 + O(\epsilon) \\
		p &= {\color{blue}p_0} + p_h + \epsilon {\color{red}p_1} + O(\epsilon^2)
	\end{align*}
```

Inserting this expansion gives 

```math
	\begin{aligned}
		\frac{\partial \rho}{\partial t} + \frac{\partial \left( \rho u\right)}{\partial x} &= - \frac{\textrm{A}_x}{\textrm{A}}\rho u \\
		\frac{\partial \left( \rho u \right)}{\partial t} + \frac{\partial}{\partial x} \left[\rho u^2 + {\color{red}p_1}\right] &= - \frac{\textrm{A}_x}{\textrm{A}}\rho u^2- \rho u \left( \beta \eta + \beta\left(1-\eta\right) \vert u \vert \right) - \frac{\textrm{h}_x}{Fr^2}(\rho - \rho_{h_0}) \\
		\frac{\partial {\color{blue}p_0}}{\partial t} + \gamma {\color{blue}p_0} \frac{\partial u}{\partial x} &= - \gamma \frac{\textrm{A}_x}{\textrm{A}} u {\color{blue}p_0} - \frac{k_w}{\textrm{A} \sqrt{\textrm{A}}} \left(T - \textrm{T}_\textrm{u} \right) \\
		\frac{\partial {\color{blue}p_0}}{\partial x} &= 0 \\
		{\color{blue}p_0} &= \rho T \\
		\frac{\partial T_i}{\partial t} &= k_i \left( T - T_i \right).		
    \end{aligned}
```

Unknowns:

- Density ``\rho``

- Velocity ``u``

- Leading order pressure ``{\color{blue}p_0}``

- First order pressure ``{\color{red}p_1}``

- Air temperature ``T``

- Material temperature ``T_i``

Due to the flow channel design we impose periodic boundary conditions.
As a result ``{\color{blue}p_0}`` is constant in space. 
We assume to have appropriate initial states given by ``\rho^0(x)``, ``u^0(x)``, ``p_0^0`` and ``T_i^0(x)``, respectively.
The boundary conditions are ``p_1(t, 0) = 0``,  ``p_1(t, 1) = 0`` and ``\textrm{A}(0)\rho(t,0)u(t,0) = \textrm{A}(1)\rho(t,1)u(t,1)``.

## Stationary case

In the stationary case we make two assumptions for the internal mound temperature.
Firstly, we assume ``T_i`` is constant in ``x``.
Secondly, the Ordinary Differential equation (ODE) for ``T_i`` is replaced by an algebraic equation, assuming the heat flux of the stationary system is in an equilibrium state.
The stationary equations are

```math
	\begin{aligned}
		\frac{\partial \left( \rho u\right)}{\partial x} &= - \frac{\textrm{A}_x}{\textrm{A}}\rho u \\
		\frac{\partial}{\partial x} \left[\rho u^2 + {\color{red}p_1}\right] &= - \frac{\textrm{A}_x}{\textrm{A}}\rho u^2- \rho u \left( \beta \eta + \beta\left(1-\eta\right) \vert u \vert \right) - \frac{\textrm{h}_x}{Fr^2}(\rho - \rho_{h_0}) \\
		\gamma {\color{blue}p_0} \frac{\partial u}{\partial x} &= - \gamma \frac{\textrm{A}_x}{\textrm{A}} u {\color{blue}p_0} - \frac{k_w}{\textrm{A} \sqrt{\textrm{A}}} \left(T - \textrm{T}_\textrm{u} \right) \\
		\frac{\partial {\color{blue}p_0}}{\partial x} &= 0 \\
		{\color{blue}p_0} &= \rho T \\
		0 &= \int_0^1  \frac{k_w}{\gamma A \sqrt{A}}  \left( T - T_u \right).	
    \end{aligned}
```

With appropriate initial and boundary conditions we have ``{\color{blue}p_0}=1``.

Unknowns:

- Density ``\rho``

- Velocity ``u``

- Pressure ``{\color{red}p_1}``

- Temperature ``T``

- Material temperature ``T_i`` (scalar)

The system can be implemented and solved as a differential algebraic equation (DAE).

## Numerical reformulation

We make the assumption that the spatial variation of the velocity satisfies 

```math
	\begin{aligned}
		\frac{\partial\left( u\right)}{\partial x} &= - \frac{\textrm{A}_x}{\textrm{A}} u + Q,
    \end{aligned}
```
 where ``Q = \frac{1}{\gamma {\color{blue}p_0}} \left[ -\frac{\partial {\color{blue}p_0}}{\partial t} - \frac{k_w}{\textrm{A} \sqrt{\textrm{A}}} \left(T - \textrm{T}_\textrm{u} \right)\right]``.

 If ``Q`` is time independent one can write 

```math
	\begin{aligned}
		u(t,x) &= \frac{v(t)}{A(x)} + \frac{1}{A(x)} \int_0^x A(y)\left[ -\frac{\partial  {\color{blue}p_0}}{\partial t} - \frac{k_w}{\textrm{A} \sqrt{\textrm{A}}} \left(T - \textrm{T}_\textrm{u} \right)\right] \, dy.
    \end{aligned}
```

The new variable ``v`` is time dependent and it must satisfy

```math
	\begin{aligned}
		\frac{\partial v}{\partial t} &= \frac{1}{\int_0^1 \frac{\rho}{A} \, dy} \left[ \int_0^1 -\rho u u_x  - \rho u \left( \beta \eta - \beta\left(1-\eta\right) \vert u \vert \right) - \frac{\textrm{h}_x}{Fr^2}(\rho - \rho_{h_0}) \, dy \right].
    \end{aligned}
```

The ODE for ``v`` is obtained from integrating the momentum balance od the asymptotic model.
The periodic boundary condition of ``{\color{red}p_1}`` is used in this set to eliminate the variable ``{\color{red}p_1}``.

Finally the reformulated asymptotic model, i.e. **1D Termite Mound Model**, reads 


```math
	\begin{aligned}
		\frac{\partial \rho}{\partial t} + \frac{\partial \left( \rho u\right)}{\partial x} &= - \frac{\textrm{A}_x}{\textrm{A}}\rho u \\
	    \frac{\partial v}{\partial t} &= \frac{1}{\int_0^1 \frac{\rho}{A} \, dy} \left[ \int_0^1 -\rho u u_x  - \rho u \left( \beta \eta - \beta\left(1-\eta\right) \vert u \vert \right) - \frac{\textrm{h}_x}{Fr^2}(\rho - \rho_{h_0}) \, dy \right] \\
        \frac{\partial {\color{blue}p_0}}{\partial t} + \gamma {\color{blue}p_0} \frac{\partial u}{\partial x} &= - \gamma \frac{\textrm{A}_x}{\textrm{A}} u {\color{blue}p_0} - \frac{k_w}{\textrm{A} \sqrt{\textrm{A}}} \left(T - \textrm{T}_\textrm{u} \right) \\
        \frac{\partial T_i}{\partial t} &= k_i \left( T - T_i \right) \\
        {\color{blue}p_0} &= \rho T \\
        u(t,x) &= \frac{v(t)}{A(x)} + \frac{1}{A(x)} \int_0^x A(y)\left[ -\frac{\partial  {\color{blue}p_0}}{\partial t} - \frac{k_w}{\textrm{A} \sqrt{\textrm{A}}} \left(T - \textrm{T}_\textrm{u} \right)\right] \, dy.
    \end{aligned}
```

Unknowns:

- Density ``\rho``

- Velocity ``u``

- Velocity ``v``

- Leading order pressure ``{\color{blue}p_0}``

- Air temperature ``T``

- Material temperature ``T_i``

More details about the models can be found in the corresponding publication: