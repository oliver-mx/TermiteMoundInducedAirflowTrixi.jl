# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
    """
        UpdateVelocityCallback

        Callback performs explicit time steps for the velocity v and pressure p0.
        If the callback is not used, both quantities will remain constant in time!

        CarpenterKennedy2N54() is used to perform the time integration for v and p0.
        The air velocity u(t,x) (i.e. "v1") is computed using the integrated v.
          
    """

    mutable struct UpdateVelocityCallback{Vis_count}
        a::Vis_count
    end

    function Base.show(io::IO, cb::DiscreteCallback{<:Any, <:UpdateVelocityCallback})
        @nospecialize cb # reduce precompilation time

        update_velocity_callback = cb.affect!
        @unpack a, b = update_velocity_callback
        print(io, "UpdateVelocityCallback(a = ", a, ")")
    end

    function Base.show(io::IO, ::MIME"text/plain",
                       cb::DiscreteCallback{<:Any, <:UpdateVelocityCallback})
        @nospecialize cb # reduce precompilation time

        if get(io, :compact, false)
            show(io, cb)
        else
            update_velocity_callback = cb.affect!

            setup = [
                "Vis Count" => update_velocity_callback.a#,
            #"Error" => update_velocity_callback.b
            ]
            Trixi.summary_box(io, "UpdateVelocityCallback", setup)
        end
    end

    function UpdateVelocityCallback(; a = 1::Int)
        # Convert plain real numbers to functions for unified treatment
        a_conv = isa(a, Real) ? Returns(a) : a
        update_velocity_callback = UpdateVelocityCallback{typeof(a_conv)}(a_conv)

        DiscreteCallback(condition, update_velocity_callback;
                         save_positions = (false, false))
    end

    # callback always activated
    @inline function condition(u, t, integrator)
        return true
    end

    # This method is called as callback during the time integration.
    @inline function (update_velocity_callback::UpdateVelocityCallback)(integrator)
        """
        "update_velocity_callback" updates the velocity "v1" and pressure "p0" after each time step.
        This function is essential to solve the reformulated asymptotic model.

        Note: If the callback is not used, both "v1" and "p0" will remain constant in time!

        The CarpenterKennedy2N54() method is used to perform the time integration for v.
        Then, the air velocity u(t,x) (i.e. "v1") is computed using the integrated value of v.
        Finally, CarpenterKennedy2N54() is used again for "p0" and the state of the integrator is updated.
        """
        #-----------------------
        original_nodes = integrator.sol.prob.p.cache.elements.node_coordinates
        nodes = Float64.(vec(original_nodes[1, :, :]))
        L_nodes = length(nodes)
        t = LinRange(nodes[1], nodes[end], L_nodes)
        ti = t[1]:0.01:t[end]
        #-----------------------
        equations = integrator.sol.prob.p.equations
        tspan = integrator.sol.prob.tspan
        #-----------------------
        # TODO: use weights from integrator and remove FastGaussQuadrature dependence
        # try: write functions for quantities that require integration.
        # semi = integrator.sol.prob.p
        # integrate(f,semi)
        #
        # or: wrap_array to get elements for each variable
        # integrate_via_indices
        nodes_unique, idx_unique = unique_idx(nodes, equations)
        L_nodes_unique = length(nodes_unique)
        x, w = gausslegendre(L_nodes_unique)
        weights_unique = w ./ 2
        weights = zeros(L_nodes)
        weights[idx_unique] = weights_unique
        #-----------------------
        A_nodes = [A(x, equations) for x in nodes]
        A_x_nodes = [A_x(x, equations) for x in nodes]
        A_nodes_inv = [inv(A(x, equations)) for x in nodes]
        A_sqrt_A_inv = [inv(sqrt(A(x, equations))) for x in nodes]
        I_w_nodes = [I_w(x, equations) for x in nodes]
        #-----------------------
        t_prev = integrator.tprev
        rho_prev = integrator.uprev[1:5:(end - 4)]
        v1_prev = integrator.uprev[2:5:(end - 3)]
        v1_prev_LI = bspline2linear(nodes, v1_prev, t, ti, equations)
        v1_dx_prev = [Interpolations.gradient(v1_prev_LI, nodes[i])[1] for i in 1:L_nodes]
        v_prev = v1_prev[1] * A(0, equations)
        Ti_prev = integrator.uprev[4:5:(end - 1)]
        p0_prev = integrator.uprev[3:5:(end - 2)]
        T_prev = p0_prev ./ rho_prev
        Tu_prev = [T_u(t_prev, x, y, equations) * equations.t_ref
                   for (x, y) in zip(nodes, Ti_prev)]
        #-----------------------
        # \frac{\partial p_0}{\partial t} = - \gamma p_0 \frac{\partial u}{\partial x} - \gamma \frac{\textrm{A}_x}{\textrm{A}} u p_0 - \frac{k_w}{\textrm{A} \sqrt{\textrm{A}}} \left(T - \textrm{T}_\textrm{u} \right)
        p0_dt_prev = .- equations.γ .* p0_prev .* v1_dx_prev .-
                     equations.γ .* A_x_nodes ./ A_nodes .* v1_prev .* p0_prev .-
                     I_w_nodes .* equations.k_w .* A_sqrt_A_inv .* (T_prev .- Tu_prev)
        p0_dt_prev = sum(p0_dt_prev) / length(p0_dt_prev)
        #-----------------------
        # \frac{\partial v}{\partial t} = \frac{1}{\int_0^1 \frac{\rho}{A} \, dy} \left[ \int_0^1 -\rho u u_x  - \rho u \left( \beta \eta - \beta\left(1-\eta\right) \vert u \vert \right) - \frac{\textrm{h}_x}{Fr^2}(\rho - \rho_{h_0}) \, dy \right]
        I_inv = inv.(sum((rho_prev .* A_nodes_inv) .* weights))
        beta = [equations.β ./ A(x, equations) for x in nodes]
        f = - rho_prev .* v1_prev .* v1_dx_prev -
            beta .* rho_prev .* v1_prev .*
            (equations.η .- (1 - equations.η) .* abs.(v1_prev)) -
            h_x.(nodes, equations) .* (rho_prev .- equations.ρₕ₀) ./ equations.Fr²
        F = sum((f) .* weights)
        v_dt_prev = I_inv .* F
        #-----------------------
        t_now = integrator.t
        #-----------------------
        rho = integrator.u[1:5:(end - 4)]
        if any(isnan, rho)
            println(rho)
            error("NaN detected in rho!")
        end
        Ti = integrator.u[4:5:(end - 1)]
        p0 = integrator.u[3:5:(end - 2)]
        T = p0 ./ rho
        Tu = [T_u(t_now, x, y, equations) * equations.t_ref for (x, y) in zip(nodes, Ti)]
        #-----------------------
        ### Explixit velocity time step via CarpenterKennedy2N54() ###
        direction_prev = sign(v_prev)
        v_exp = v_prev + integrator.dt * v_dt_prev
        direction_exp = sign(v_exp)
        v_dt = if direction_prev == direction_exp
            v1_exp = v1_prev .+ (v_exp - v_prev) / A(0.0, equations)
            v1_dx_exp = v1_dx_prev
            I_inv = inv.(sum((rho .* A_nodes_inv) .* weights))
            f = - rho .* v1_exp .* v1_dx_exp -
                beta .* rho .* v1_exp .*
                (equations.η .- (1 - equations.η) .* abs.(v1_exp)) -
                h_x.(nodes, equations) .* (rho .- equations.ρₕ₀) ./ equations.Fr²
            F = sum((f) .* weights)
            I_inv .* F
        else
            v1_exp = reverse(- v1_prev .+ (v_exp + v_prev) / A(0.0, equations))
            v1_exp_LI = bspline2linear(nodes, v1_exp, t, ti, equations)
            v1_dx_exp = [Interpolations.gradient(v1_exp_LI, nodes[i])[1] for i in 1:L_nodes]
            I_inv = inv.(sum((rho .* A_nodes_inv) .* weights))
            f = - rho .* v1_exp .* v1_dx_exp -
                beta .* rho .* v1_exp .*
                (equations.η .- (1 - equations.η) .* abs.(v1_exp)) -
                h_x.(nodes, equations) .* (rho .- equations.ρₕ₀) ./ equations.Fr²
            F = sum((f) .* weights)
            I_inv .* F
        end
        LI = LinearInterpolation([t_prev, t_now], [v_dt_prev, v_dt],
                                 extrapolation_bc = Line())
        prob = ODEProblem((u, p, t) -> LI(t), v_prev, (t_prev, t_now))
        sol = solve(prob,
                    CarpenterKennedy2N54(williamson_condition = integrator.sol.alg.williamson_condition),
                    dt = integrator.dt)
        v = sol.u[end]
        #-----------------------
        ### update u(t,x) ###
        Q = (- I_w_nodes .* equations.k_w .* A_sqrt_A_inv .* (T_prev .- Tu_prev)) ./
            (equations.γ .* p0_prev[1])
        c = sum((Q .* A_nodes) .* weights)
        I_0 = cumsum((Q .* A_nodes .- c) .* weights)
        v1 = v .* A_nodes_inv .+ A_nodes_inv .* I_0
        v1_LI = LinearInterpolation2(nodes, v1, equations)
        v1_dx = [Interpolations.gradient(v1_LI, nodes[i])[1] for i in 1:L_nodes]
        #------------------------
        p0_dt = .- equations.γ .* p0 .* v1_dx .- equations.γ .* A_x_nodes ./ A_nodes .* v1 .* p0 .-
                I_w_nodes .* equations.k_w .* A_sqrt_A_inv .* (T .- Tu)
        p0_dt = sum(p0_dt) / length(p0_dt)
        LI = LinearInterpolation([t_prev, t_now], [p0_dt_prev, p0_dt],
                                 extrapolation_bc = Line())
        prob = ODEProblem((u, p, t) -> LI(t), p0_prev[1], (t_prev, t_now))
        sol = solve(prob, CarpenterKennedy2N54(williamson_condition = false),
                    dt = integrator.dt)
        p0 = sol.u[end]
        #-----------------------
        integrator.u[2:5:(end - 3)] = v1
        integrator.u[3:5:(end - 2)] .= p0
        #-----------------------
        #update_velocity_callback.a()
        #new_a = 0
        #update_velocity_callback.a = isa(new_a, Real) ? Returns(new_a) : new_a
        #-----------------------
        return integrator
    end
end # @muladd