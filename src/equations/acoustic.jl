struct Acoustic <: Equation end
@declare_dofs Acoustic [:u, :v, :pressure, :rho, :lambda]
struct GaussianWave <: Scenario
end

function is_periodic_boundary(equation::Acoustic, scenario::GaussianWave)
    # GaussianWave scenario doesn't require periodic boundary conditions
    false
end

function evaluate_boundary(eq::Acoustic, scenario::GaussianWave, face, normalidx, dofsface, dofsfaceneigh)
# dofsface and dofsfaceneigh have shape (num_2d_quadpoints, dofs)
# you need to set dofsfaceneigh
    dofsfaceneigh .= dofsface
    if normalidx == 1
        dofsfaceneigh[:, 1] = -dofsface[:, 1]
    else
        dofsfaceneigh[:, 2] = -dofsface[:, 2]
    end
end


function get_initial_values(eq::Acoustic, scenario::GaussianWave, global_position; t=0.0)
    x, y = global_position
    if is_analytical_solution(eq, scenario)
        println("No analytical solution implemented")
    else
        [ 0, 0, exp(-100*(x - 0.5)^2 - 100*(y - 0.5)^2), 1.0, if x<= 0.5 1/5 else 1 end]
    end
end

function is_analytical_solution(equation::Acoustic, scenario::GaussianWave)
    false
end


function evaluate_flux(eq::Acoustic, celldofs, cellflux)
    # size(cellflux)=(2*order*order, ndof)=(2,5)
    # size(celldofs)=(order*order, ndof)=(1,5)

    order_sq, ndof = size(celldofs)
    for j=1:order_sq
        # explicit names of dofs for clarity
        vx = celldofs[j, 1]
        vy = celldofs[j, 2]
        pressure = celldofs[j, 3]
        density = celldofs[j, 4]
        bulk = celldofs[j, 5]

        # here flux is computed elementwise
        cellflux[j,1]= pressure / density
        cellflux[j + order_sq, 1] = 0

        cellflux[j,2] = 0
        cellflux[j + order_sq, 2] = pressure / density

        cellflux[j,3] = bulk*vx
        cellflux[j + order_sq, 3] = bulk*vy

        cellflux[j,4] = 0
        cellflux[j + order_sq, 4] = 0

        cellflux[j,5] = 0
        cellflux[j + order_sq, 5] = 0

    end

end

function max_eigenval(eq::Acoustic, celldata, normalidx)
    # sqrt( K_0 / rho_0 )
    maxeigenval = 0.
    for i in 1:size(celldata)[1]
        vx, vy, pressure, density, bulk = celldata[i, :]
        maxeigenval_new = sqrt(bulk / density)
        if maxeigenval_new > maxeigenval
            maxeigenval = maxeigenval_new
        end
    end
    maxeigenval
end
