# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin

    @doc raw"""
        UpdateVelocityCallback(solver)

    ### Parameters
        - `solver`: ODE solver.

    The `UpdateVelocityCallback` called after every time step.
    In this callback explicit time steps for the velocity ``v`` and pressure ``p_0`` are performed.
    For the solver we reccomend `CarpenterKennedy2N54` from `OrdinaryDiffEqLowStorageRK`.
    After the time integration, the velocity ``u`` is computed and the integrator is updated.

    !!! note
        If the callback is not used, the quantities ``u``, ``v`` and ``p_0`` will remain constant in time!
    """
    mutable struct UpdateVelocityCallback{ODE_solver,Vis_count}
        solver::ODE_solver
        a::Vis_count
    end

    function Base.show(io::IO, cb::DiscreteCallback{<:Any,<:UpdateVelocityCallback})
        @nospecialize cb # reduce precompilation time

        update_velocity_callback = cb.affect!
        @unpack solver, a = update_velocity_callback
        print(io, "UpdateVelocityCallback(solver = ", solver, ", a = ", a, ")")
    end

    function Base.show(
        io::IO,
        ::MIME"text/plain",
        cb::DiscreteCallback{<:Any,<:UpdateVelocityCallback},
    )
        @nospecialize cb # reduce precompilation time

        if get(io, :compact, false)
            show(io, cb)
        else
            update_velocity_callback = cb.affect!

            setup = ["Solver" => update_velocity_callback.solver]
            Trixi.summary_box(io, "UpdateVelocityCallback", setup)
        end
    end

    function UpdateVelocityCallback(solver; a = 1::Int)
        # Convert plain real numbers to functions for unified treatment
        solver_conv = isa(solver, Real) ? Returns(solver) : solver
        a_conv = isa(a, Real) ? Returns(a) : a
        update_velocity_callback =
            UpdateVelocityCallback{typeof(solver_conv),typeof(a_conv)}(solver_conv, a_conv)

        DiscreteCallback(
            condition,
            update_velocity_callback;
            save_positions = (false, false),
        )
    end

    @doc raw"""
        condition(u, t, integrator)

    ### Parameters
    - `u`: State vector.
    - `t`: Time scalar.
    - `integrator`: ODE integrator.

    Ensures that `UpdateVelocityCallback` is used after every time step.

    ### Returns
    True.
    """
    @inline function condition(u, t, integrator)
        return true
    end

    @inline function (update_velocity_callback::UpdateVelocityCallback)(integrator)
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
        rho_prev = integrator.uprev[1:5:(end-4)]
        v1_prev = integrator.uprev[2:5:(end-3)]
        v1_prev_LI = bspline2linear(nodes, v1_prev, t, ti, equations)
        v1_dx_prev =
            [Interpolations.gradient(v1_prev_LI, nodes[i])[1] for i = 2:(L_nodes-1)]
        v1_dx_prev = vcat(
            Interpolations.gradient(v1_prev_LI, 0.0001)[1],
            vcat(v1_dx_prev, Interpolations.gradient(v1_prev_LI, 0.9999)[1]),
        )
        v_prev = v1_prev[1] * A(0, equations)
        Ti_prev = integrator.uprev[4:5:(end-1)]
        p0_prev = integrator.uprev[3:5:(end-2)]
        T_prev = p0_prev ./ rho_prev
        Tu_prev = [
            T_u(t_prev, x, y, equations) * equations.t_ref for (x, y) in zip(nodes, Ti_prev)
        ]
        #-----------------------
        p0_dt_prev =
            .- equations.γ .* p0_prev .* v1_dx_prev .-
            equations.γ .* A_x_nodes ./ A_nodes .* v1_prev .* p0_prev .-
            I_w_nodes .* equations.k_w .* A_sqrt_A_inv .* (T_prev .- Tu_prev)
        p0_dt_prev = sum(p0_dt_prev) / length(p0_dt_prev)
        #-----------------------
        I_inv = inv.(sum((rho_prev .* A_nodes_inv) .* weights))
        beta = [equations.β ./ A(x, equations) for x in nodes]
        f =
            - rho_prev .* v1_prev .* v1_dx_prev -
            beta .* rho_prev .* v1_prev .*
            (equations.η .- (1 - equations.η) .* abs.(v1_prev)) -
            h_x.(nodes, equations) .* (rho_prev .- equations.ρₕ₀) ./ equations.Fr²
        F = sum((f) .* weights)
        v_dt_prev = I_inv .* F
        #-----------------------
        t_now = integrator.t
        #-----------------------
        rho = integrator.u[1:5:(end-4)]
        if any(isnan, rho)
            error("NaN detected in rho!")
        end
        Ti = integrator.u[4:5:(end-1)]
        p0 = integrator.u[3:5:(end-2)]
        T = p0 ./ rho
        Tu = [T_u(t_now, x, y, equations) * equations.t_ref for (x, y) in zip(nodes, Ti)]
        #-----------------------
        direction_prev = sign(v_prev)
        v_exp = v_prev + integrator.dt * v_dt_prev
        direction_exp = sign(v_exp)
        v_dt = if direction_prev == direction_exp
            v1_exp = v1_prev .+ (v_exp - v_prev) / A(0.0, equations)
            v1_dx_exp = v1_dx_prev
            I_inv = inv.(sum((rho .* A_nodes_inv) .* weights))
            f =
                - rho .* v1_exp .* v1_dx_exp -
                beta .* rho .* v1_exp .*
                (equations.η .- (1 - equations.η) .* abs.(v1_exp)) -
                h_x.(nodes, equations) .* (rho .- equations.ρₕ₀) ./ equations.Fr²
            F = sum((f) .* weights)
            I_inv .* F
        else
            v1_exp = reverse(- v1_prev .+ (v_exp + v_prev) / A(0.0, equations))
            v1_exp_LI = bspline2linear(nodes, v1_exp, t, ti, equations)
            v1_dx_exp =
                [Interpolations.gradient(v1_exp_LI, nodes[i])[1] for i = 1:L_nodes]
            I_inv = inv.(sum((rho .* A_nodes_inv) .* weights))
            f =
                - rho .* v1_exp .* v1_dx_exp -
                beta .* rho .* v1_exp .*
                (equations.η .- (1 - equations.η) .* abs.(v1_exp)) -
                h_x.(nodes, equations) .* (rho .- equations.ρₕ₀) ./ equations.Fr²
            F = sum((f) .* weights)
            I_inv .* F
        end
        LI = linear_interpolation(
            [t_prev, t_now],
            [v_dt_prev, v_dt],
            extrapolation_bc = Line(),
        )
        prob = ODEProblem((u, p, t) -> LI(t), v_prev, (t_prev, t_now))
        sol = solve(prob, update_velocity_callback.solver, dt = integrator.dt)
        v = sol.u[end]
        #-----------------------
        Q =
            (- I_w_nodes .* equations.k_w .* A_sqrt_A_inv .* (T_prev .- Tu_prev)) ./
            (equations.γ .* p0_prev[1])
        c = sum((Q .* A_nodes) .* weights)
        I_0 = cumsum((Q .* A_nodes .- c) .* weights)
        v1 = v .* A_nodes_inv .+ A_nodes_inv .* I_0
        v1_LI = linear_interpolation2(nodes, v1, equations)
        v1_dx = [Interpolations.gradient(v1_LI, nodes[i])[1] for i = 1:L_nodes]
        #------------------------
        p0_dt =
            .- equations.γ .* p0 .* v1_dx .- equations.γ .* A_x_nodes ./ A_nodes .* v1 .* p0 .-
            I_w_nodes .* equations.k_w .* A_sqrt_A_inv .* (T .- Tu)
        p0_dt = sum(p0_dt) / length(p0_dt)
        LI = linear_interpolation(
            [t_prev, t_now],
            [p0_dt_prev, p0_dt],
            extrapolation_bc = Line(),
        )
        prob = ODEProblem((u, p, t) -> LI(t), p0_prev[1], (t_prev, t_now))
        sol = solve(prob, update_velocity_callback.solver, dt = integrator.dt)
        p0 = sol.u[end]
        #-----------------------
        integrator.u[2:5:(end-3)] = v1
        integrator.u[3:5:(end-2)] .= p0
        #-----------------------
        #update_velocity_callback.a()
        #new_a = 0
        #update_velocity_callback.a = isa(new_a, Real) ? Returns(new_a) : new_a
        #-----------------------
        return integrator
    end
end # @muladd
