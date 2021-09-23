struct Euler <: Equation end
@declare_dofs Euler [:rhou, :rhov, :rho, :e]

struct GaussianWaveeuler <: Scenario end
struct Sod_shock_tube <: Scenario end


function is_periodic_boundary(equation::Euler, scenario::GaussianWaveeuler)
    # GaussianWave scenario with periodic boundary conditions
    true
end

function is_periodic_boundary(equation::Euler, scenario::Sod_shock_tube)
    # sod_shock_tube scenario doesn't require periodic boundary conditions
    false
end

# for Sod_shock_tube Scenario
function evaluate_boundary(eq::Euler, scenario::Sod_shock_tube, face, normalidx, dofsface, dofsfaceneigh)
# dofsface and dofsfaceneigh have shape (num_2d_quadpoints, dofs)
# you need to set dofsfaceneigh
    dofsfaceneigh .= dofsface
    if normalidx == 1
        dofsfaceneigh[:, 1] = -dofsface[:, 1]
    else
        dofsfaceneigh[:, 2] = -dofsface[:, 2]
    end
end

function evaluate_energy(eq::Euler, rhou, rhov, rho, p)
    # value of gamma is 1.4, hence gamma-1 = 0.4
    return ( (p/0.4)  +  (1/(2*rho)) * (rhou^2 + rhov^2)  )
end

function evaluate_pressure(eq::Euler, rhou, rhov, rho, e)
    # value of gamma is 1.4, hence gamma-1 = 0.4
    return ( (0.4) * ( e - (1/(2*rho)) * (rhou^2 + rhov^2) ))
end

function get_initial_values(eq::Euler, scenario::GaussianWaveeuler, global_position; t=0.0)
    x, y = global_position
    if is_analytical_solution(eq, scenario)
        println("No analytical solution implemented")
    else
        p = exp(-100 * (x - 0.5)^2 - 100 *(y - 0.5)^2) + 1
        [0, 0, 1.0, evaluate_energy(eq, 0, 0, 1.0, p)]
    end
end

function get_initial_values(eq::Euler, scenario::Sod_shock_tube, global_position; t=0.0)
    x, y = global_position
    if is_analytical_solution(eq, scenario)
        println("No analytical solution implemented")
    else
        [ 0, 0, if x<=0.5 0.125 else 1.0 end, if x<= 0.5 evaluate_energy(eq, 0, 0, 0.125, 0.1) else evaluate_energy(eq, 0, 0, 1.0, 1.0) end ]
    end
end

function is_analytical_solution(equation::Euler, scenario::GaussianWaveeuler)
    false
end

function is_analytical_solution(equation::Euler, scenario::Sod_shock_tube)
    false
end

function evaluate_flux(eq::Euler, celldofs, cellflux)
    # size(cellflux)=(2*order*order, ndof)=(2,5)
    # size(celldofs)=(order*order, ndof)=(1,5)

    order_sq, ndof = size(celldofs)
    for j=1:order_sq
        # explicit names of dofs for clarity
        rhou = celldofs[j, 1]
        rhov = celldofs[j, 2]
        density = celldofs[j, 3]
        e = celldofs[j, 4]

        p=evaluate_pressure(eq,rhou,rhov,density,e)

        # here flux is computed elementwise
        cellflux[j,1]= (1/density)*rhou^2 + p
        cellflux[j + order_sq, 1] = (1/density)*rhou*rhov

        cellflux[j,2] = (1/density)*rhou*rhov
        cellflux[j + order_sq, 2] = (1/density)*rhov^2 + p

        cellflux[j,3] = rhou
        cellflux[j + order_sq, 3] = rhov

        cellflux[j,4] = (rhou/density)*(e+p)
        cellflux[j + order_sq, 4] = (rhov/density)*(e+p)

    end

end

function max_eigenval(eq::Euler, celldata, normalidx)
    # v_n + c, where c = sqrt(gamma*p /rho)
    maxeigenval = 0.
    #@info size(celldata)
    for i in 1:size(celldata)[1]
        rhou, rhov, density, e = celldata[i, :]
        pressure = evaluate_pressure(eq,rhou,rhov,density,e)
        #@info i
        #@info pressure
        #@info density
        #@info e
        #@info rhou
        #@info rhov
        c = sqrt(1.4*pressure/density)
        if normalidx == 1
            vn = rhou / density
        else
            vn = rhov / density
        end
        maxeigenval_new = vn + c
        if maxeigenval_new > maxeigenval
            maxeigenval = maxeigenval_new
        end
    end
    maxeigenval
end
