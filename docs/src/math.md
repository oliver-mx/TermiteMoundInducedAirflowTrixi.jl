This thermo fluid dynamic model describes the airflow dynamics inside of termite mounds. The one-dimensional mathematicalmodel originates from the Euler equations of gas dynamics in a low Mach number regime. 

```math
    \begin{equation}
		\frac{\partial \rho}{\partial t} + \frac{\partial \left( \rho u\right)}{\partial x} = - \frac{\textrm{A}_x}{\textrm{A}}\rho u \\
		\frac{\partial v}{\partial t} = \frac{1}{\int_0^1 \frac{\rho}{A} \, dy} \left[ 
		\int_0^1 -\rho u u_x  - \rho u \left( \beta \eta - \beta\left(1-\eta\right) \vert u \vert \right) - \frac{\textrm{h}_x}{Fr^2}(\rho - \rho_{h_0}) \, dy \right] \\
		\frac{\partial p_0}{\partial t} + \gamma p_0 \frac{\partial u}{\partial x}= - \gamma \frac{\textrm{A}_x}{\textrm{A}} u p_0 - \frac{k_w}{\textrm{A} \sqrt{\textrm{A}}} \left(T - \textrm{T}_\textrm{u} \right) \\
		\frac{\partial T_i}{\partial t} = k_i \left( T - T_i \right) \\
		p_0 = \rho T \\
		u(t,x) = \frac{v(t)}{A(x)} + \frac{1}{A(x)} \int_0^x A(y)\left[ -\frac{\partial  p_0}{\partial t} - \frac{k_w}{\textrm{A} \sqrt{\textrm{A}}} \left(T - \textrm{T}_\textrm{u} \right)\right] \, dy 
	\end{equation}
```