## Tutorial for the 1D model

In this section, we describe how to use **TermiteMoundInducedAirflowTrixi.jl** with **Trixi.jl**. This tutorial will guide you through setting up and running a 1D blood flow simulation, including mesh creation, boundary conditions, numerical fluxes, and visualization of results.

### Packages 

Before starting, ensure that the required packages are loaded:
```julia
julia> using Pkg

julia> Pkg.add(["TermiteMoundInducedAirflowTrixi", "Trixi", "Trixi2Vtk", "OrdinaryDiffEqLowStorageRK", "Plots"])

julia> using TermiteMoundInducedAirflowTrixi, Trix, Trixi2Vtk, OrdinaryDiffEqLowStorageRK, Plots          
```

#### Semidiscretization

First, we need to choose the equation that describes the dynamics. Is optional to set specific parameter values:
```julia
julia> equations = TermiteMoundEquations1D(;η=0.93)
```

Next, we prescribe the initial state, volume and surface flux of all variable in our equation system:
```julia
julia> initial_condition = TermiteMoundInitialCondition
julia> volume_flux = flux_ranocha
julia> equations = surface_flux = flux_ranocha
```

The polynomial degree for the `DGSEM` can be specified now:
```
julia> dg = DGSEM(polydeg = 3,surface_flux = flux_ranocha,volume_integral = VolumeIntegralFluxDifferencing(volume_flux),)
```

We define a one-dimensional Tree mesh, which discretizes the (scaled) spatial domain:
```
julia> mesh = TreeMesh((0.0,),(1.0,), initial_refinement_level = 5, n_cells_max = 1000, periodicity = true,)
```

The semidiscretization is now defined:
```
julia> semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, dg, source_terms = source_terms, boundary_conditions = boundary_condition_periodic,);
```

#### ODE solvers and callbacks

We want to simulate a full day. In this case we set the time span to:
```
julia> mesh = tspan = (0.0, 1.0).* 86400 ./ equations.tᵣ
julia> ode = semidiscretize(semi, tspan)
```

We want to use a `Summary`, `Stepsize` and `AMRCallback` in addition to the mandatory `UpdateVelocityCallback`:
```
julia> summary_callback = SummaryCallback()
julia> stepsize_callback = StepsizeCallback(cfl = 0.5)
julia> amr_controller = ControllerThreeLevel(semi, IndicatorMax(semi, variable = (u, equations) -> u[2]), base_level = 5, max_level = 6, max_threshold = -0.35)
julia> amr_callback = AMRCallback(semi, amr_controller, interval = 1, adapt_initial_condition = true, adapt_initial_condition_only_refine = true)
julia> callbacks = CallbackSet(summary_callback, stepsize_callback, amr_callback, update_velocity_callback)
```

#### Run the simulation
We use the `CarpenterKennedy2N54` time integration method from `OrdinaryDiffEq` to solve the system of equations:
 ```
julia> sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false); dt = 1.0, ode_default_options()...,callback = callbacks, );
```

#### Plot the results
The initial state of the flow quantities can be visualized by:
 ```
julia> pd = PlotData1D((x, equations) -> initial_condition(x, last(tspan), equations), semi); plot(pd)
```

The final state of the flow quantities can be visualized via:
 ```
julia> pd = PlotData1D(sol); plot(pd)
```
