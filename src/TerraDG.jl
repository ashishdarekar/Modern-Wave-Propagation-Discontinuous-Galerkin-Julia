module TerraDG
using WriteVTK
using Printf
using Logging
using LinearAlgebra
import YAML

include("configuration.jl")
include("basis.jl")
include("equations.jl")
include("grid.jl")
include("kernels/surface.jl")
include("kernels/volume.jl")
include("kernels/filtering.jl")
include("kernels/time.jl")
include("plotters.jl")
include("error_writer.jl")
include("global_matrices.jl")
include("kernels/limiter.jl")

"""
    evaluate_rhs(eq, scenario, filter, globals, du, dofs, grid)

Evalutes the right-hand-side of the equation `eq` for
scenario `scenario`, with filter `filter`,
collection of global matrices `globals`, update
`du`, degrees of freedom `dofs` and grid `grid`.

Updates `du` in place.
"""
function evaluate_rhs(eq, scenario, filter, globals, du, dofs, grid, config)
    buffers_face = BuffersFaceIntegral(grid.basis, get_ndofs(eq))
    buffers_volume = BuffersVolume(grid.basis, get_ndofs(eq))
    linearindex = vec(collect(Iterators.product(1:grid.basis.order, 1:grid.basis.order))) #linear indexing
    var = get_ndofs(eq)
    reference_massmatrix = massmatrix(grid.basis, grid.basis.dimensions)
    du .= 0.0
    maxeigenval = -Inf

    ∇ = globals.reference_derivative_matrix
    for i in eachindex(grid.cells)
        @views cell = grid.cells[i]
        @views data = dofs[:,:, cell.dataidx]
        @views flux = grid.flux[:,:,cell.dataidx]

        #evaluate_flux(eq, data, flux)

        #ATTEMPT 1
        """
        if config.source == "without_source"
            evaluate_flux(eq, data, flux)
        else
            evaluate_flux_with_source(eq, data, flux, cell, grid.basis)
        end
        """
       
        x, y = globalposition(cell, (0.5, 0.5))
        
        if config.scenario_name == "stiff_inclusion"
            evaluate_flux_with_stiffinclusion(eq, scenario, data, flux, (x,y))
        else
            evaluate_flux(eq, data, flux)
        end

        # Volume matrix is zero for FV/order=1
        if length(grid.basis.quadpoints) > 1
            @views evaluate_volume(globals, buffers_volume, flux, grid.basis, inverse_jacobian(cell), volume(cell), du[:,:,cell.dataidx])
            #@info "after VI"
            #@info du[:,:,cell.dataidx]
        end
    end
    for i in eachindex(grid.cells)
        @views cell = grid.cells[i]
        @views data = dofs[:,:, cell.dataidx]
        @views flux = grid.flux[:,:,cell.dataidx]
        elem_massmatrix = volume(cell) * reference_massmatrix
        inv_massmatrix = inv(elem_massmatrix)

        # Here we also need to compute the maximum eigenvalue of each cell
        # and store it for each cell (needed for timestep restriction later!)
        faces = [left, top, right, bottom]
        for (i, neigh) in enumerate(cell.neighbors)
            @views dofsneigh = dofs[:,:,neigh.dataidx]
            @views fluxneigh = grid.flux[:,:,neigh.dataidx]

            # @info data
            # Project dofs and flux of own cell to face
            project_to_faces(globals, data, flux, buffers_face.dofsface, buffers_face.fluxface, faces[i])
            facetypeneigh = cell.facetypes[i]

            #@info data

            if facetypeneigh == regular
                # Project neighbors to faces
                # Neighbor needs to project to opposite face
                faceneigh = globals.oppositefaces[faces[i]]
                project_to_faces(globals, dofsneigh, fluxneigh, buffers_face.dofsfaceneigh, buffers_face.fluxfaceneigh, faceneigh)
                #@info "from if"
                #@info data
            else
                #@info "from else"
                #@info data
                @assert(facetypeneigh == boundary)
                normalidx = globals.normalidxs[faces[i]]

                # For boundary cells, we operate directly on the dofsface
                evaluate_boundary(eq, scenario, faces[i], normalidx, buffers_face.dofsface, buffers_face.dofsfaceneigh)
                # Evaluate flux on face directly
                # Note: When extrapolating, this is not exact!
                # The error is given by the commutation error of face projection and flux!
                
                #evaluate_flux(eq, buffers_face.dofsfaceneigh, buffers_face.fluxfaceneigh)
             
                # ATTEMPT 1
                """
                if config.source == "without_source"
                    evaluate_flux(eq, buffers_face.dofsfaceneigh, buffers_face.fluxfaceneigh)
                else
                    evaluate_flux_with_source(eq, buffers_face.dofsfaceneigh, buffers_face.fluxfaceneigh, neigh, grid.basis)
                end
                """
                x, y = globalposition(cell, (0.5, 0.5))
                if config.scenario_name == "stiff_inclusion"
                    evaluate_flux_with_stiffinclusion(eq, scenario, buffers_face.dofsfaceneigh, buffers_face.fluxfaceneigh, (x,y))
                else
                    evaluate_flux(eq, buffers_face.dofsfaceneigh, buffers_face.fluxfaceneigh)
                end
                
                #@info data
            end
            @views cureigenval = evaluate_face_integral(eq, globals, buffers_face, cell, faces[i], du[:,:,cell.dataidx], config.riemann)
            maxeigenval = max(maxeigenval, cureigenval)
        end

        @views du[:,:,cell.dataidx] = inv_massmatrix * @views du[:,:,cell.dataidx] #
    end
    grid.maxeigenval = maxeigenval
end

"""
    main(configfile::String)

Runs a DG-simulation with configuration from `configfile`.
"""
function main(configfile::String)
    config = Configuration(configfile)
    filter = make_filter(config)
    eq = make_equation(config)
    scenario = make_scenario(config)
    grid = make_grid(config, eq, scenario)
    integrator = make_timeintegrator(config, grid)

    #@info size(grid.dofs)
    #@info grid.dofs[:,:,40]
    #@info grid.cells[5].size
    #@info grid.cells[5].facetypes

    @info "Initialising global matrices"
    globals = GlobalMatrices(grid.basis, filter, grid.basis.dimensions)
    @info "Initialised global matrices"

    filename = "output/plot"

    # Init everything
    for cell in grid.cells
        @views interpolate_initial_dofs(eq, scenario, grid.dofs[:,:,cell.dataidx],cell,grid.basis)
    end
    plotter = VTKPlotter(eq, scenario, grid, filename)

    grid.time = 0
    timestep = 0
    next_plotted = config.plot_start

    while grid.time < config.end_time
        if timestep > 0

            time_start = time()
            if config.equation_type == "non_linear"
                minmod_limiter(eq, scenario, grid.dofs, grid, grid.basis, get_ndofs(eq), 2)
            end
            dt = 1/(config.order^2+1) * config.cellsize[1] * config.courant * 1/grid.maxeigenval
            # Only step up to either end or next plotting
            dt = min(dt, next_plotted-grid.time, config.end_time - grid.time)
            @assert dt > 0

            @info "Running timestep" timestep dt grid.time
            step(integrator, grid, dt) do du, dofs, time
                # slope limiting before each timestep
                if config.equation_type == "non_linear"
                    minmod_limiter(eq, scenario, grid.dofs, grid, grid.basis, get_ndofs(eq), 2)
                end
                evaluate_rhs(eq, scenario, filter, globals, du, dofs, grid, config)
                #Here is the end of the PROBLEM A of Fractional Step method for a problem with source 
            
                #Now we need to solve the PROBLEM B of Fractional Step Method for a problem with source 
                #We will use Forward Euler for this (order=1)
                if config.source == "with_source"
                    
                    #@info "from outside" 
                    #@info sp[1](1,2,3)
                    for i in eachindex(grid.cells)
                        @views cell = grid.cells[i]
                        @views data = dofs[:,:, cell.dataidx]
                        sp = Array{Function,1}
                        globalx,globaly = globalposition(cell, (0.5, 0.5))
                        sp = evaluate_source(eq, grid.time, (globalx,globaly))
                        # Finding global Coordinates of Cell
#Calculation is specifically for order=1 and ndof=5
                        
                        #@info globalx
                        #@info globaly
                        #globalz = 0.0 # as problem is 2 Dimensional
                        #@info size(data)
                        #@info data

#Calculation is specifically for order=1 and ndof=5
                        data[1,1] = data[1,1] + ( dt* (sp[1](globalx,globaly,grid.time)[1][1]) )
                        data[1,2] = data[1,2] + ( dt* (sp[2](globalx,globaly,grid.time)[1][1]) )
                        data[1,3] = data[1,3] + ( dt* (sp[3](globalx,globaly,grid.time)[1][1]) )
                        data[1,4] = data[1,4] + ( dt* (sp[4](globalx,globaly,grid.time)[1][1]) )
                        data[1,5] = data[1,5] + ( dt* (sp[5](globalx,globaly,grid.time)[1][1]) )

                        #@info data
                    end

                end
            end
            grid.time += dt
            time_end = time()
            time_elapsed = time_end - time_start
            @info "Timestep took" time_elapsed
        else
            # Compute initial eigenvalue (needed for dt)
            grid.maxeigenval = -1
            for cell in grid.cells
                @views celldata = grid.dofs[:,:,cell.dataidx]
                for normalidx=1:2
                    cureigenval = max_eigenval(eq, celldata, normalidx)
                    grid.maxeigenval = max(grid.maxeigenval, cureigenval)
                end
            end
        end

        if abs(grid.time - next_plotted) < 1e-10
            # slope limiting before plotting
            @info "Writing output" grid.time
            if config.equation_type == "non_linear"
                minmod_limiter(eq, scenario, grid.dofs, grid, grid.basis, get_ndofs(eq), 2)
            end
            plot(plotter)
            next_plotted = grid.time + config.plot_step
        end
        #break # For sibson to run only for t=0
        timestep += 1
    end
    save(plotter)
    if is_analytical_solution(eq, scenario)
        evaluate_error(eq, scenario, grid, grid.time)
    end
end

end