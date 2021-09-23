struct Advection <: Equation end
@declare_dofs Advection [:rho1, :rho2, :rho3]

struct PlanarWaves <: Scenario end
struct Sibson <: Scenario end


function is_periodic_boundary(equation::Advection, scenario::PlanarWaves)
    true
end

function is_periodic_boundary(equation::Advection, scenario::Sibson)
    false
end

function get_initial_values(eq::Advection, scenario::PlanarWaves, global_position; t=0.0)
    x, y = global_position
    # supports t ≠ 0 if there is an analytical solution
    if is_analytical_solution(eq, scenario)
        [ sin( 2*pi* ( (x + y) -  (2*t) ) ), sin( 2*pi* (y - t) ), 1.0] #For all t , analytical solution
    else
        [ sin(2*pi* (x + y)), sin(2*pi*y), 1.0] #at t=0
    end
end

function get_initial_values(eq::Advection, scenario::Sibson, global_position; t=0.0)
    x, y = global_position
    # supports t ≠ 0 if there is an analytical solution
    if is_analytical_solution(eq, scenario)
        println("No analytical solution yet")
    else
        [ cos( 4*π*sqrt( (x-1.25)^2 + (y-1.25)^2 ) ), 1.0, 1.0]
    end
end

function is_analytical_solution(equation::Advection, scenario::PlanarWaves)
    true
end

function is_analytical_solution(equation::Advection, scenario::Sibson)
    false
end

function evaluate_flux(eq::Advection, celldofs, cellflux)
    s = AdvectionShortcuts()

    # you can use s.rho1 == 1 etc to simplify stuff, these symbols can useful for simplicity
    # size(cellflux)=(2*order*order, ndof)=(2,3)
    # size(celldofs)=(order*order, ndof)=(1,3)

    order_sq, ndof = size(celldofs)
    for j=1:order_sq
        for i=1:ndof
            cellflux[j,i] = celldofs[j,i]
            cellflux[j + order_sq, i] = celldofs[j,i]
        end
    end
    #@info celldofs
    end

function max_eigenval(eq::Advection, celldata, normalidx)
    # Is actually correct!
    1.0
end
