module velocity
use rossby_wave_attributes
use admin
implicit none
real :: eps = 0.01 !Placeholder value

contains
    function streamfunction(x,y,t,phase1) result(psi)
        real :: x, y, t, k, l, phase1, psi
        real, dimension(64) :: k_array, l_array
        integer :: counter_l, counter_k

        !Initialisation
        k_array = linspace(-2.0,2.0,64)
        l_array = linspace(-2.0,2.0,64)
        psi = 0.
        phase1 = 0.0

        do counter_k = 1, 64
            k = k_array(counter_k)
            do counter_l = 1, 64
                l = l_array(counter_l)
                psi = psi + amplitude(k,l)*cos(k*x + l*y - dispersion(k,l)*t + phase1)
            end do
        end do
                
    end function streamfunction
    
    function potential(x,y,t,phase2) result(phi)
        real :: x, y, t, k, l, phase2, phi
        real, dimension(64) :: k_array, l_array
        integer :: counter_l, counter_k

        !Initialisation
        k_array = linspace(-2.0,2.0,64)
        l_array = linspace(-2.0,2.0,64)
        phi = 0.
        phase2 = 0.0

        do counter_k = 1, 64
            k = k_array(counter_k)
            do counter_l = 1, 64
                l = l_array(counter_l)
                phi = phi + amplitude(k,l)*cos(k*x + l*y - dispersion(k,l)*t + phase2)
            end do
        end do
    end function potential
    
    function velocity_field(x,y,t,phase1,phase2) result(velocity)
        real :: x, y, t, k, l, phase1, phase2, velocity(2)
        real, dimension(64) :: k_array, l_array
        integer :: counter_l, counter_k

        !Initialisation
        k_array = linspace(-2.0,2.0,64)
        l_array = linspace(-2.0,2.0,64)
        velocity = (/0.,0./)
        phase1 = 0.0
        phase2 = 0.0

        do counter_k = 1, 64
            k = k_array(counter_k)
            do counter_l = 1, 64
                l = l_array(counter_l)
                velocity(1) = velocity(1) + l*(1-eps)*amplitude(k,l)*sin(k*x + l*y - dispersion(k,l)*t + phase1) 
                velocity(1) = velocity(1)- k*eps*amplitude(k,l)*sin(k*x + l*y - dispersion(k,l)*t + phase2)
                velocity(2) = velocity(2) - k*(1-eps)*amplitude(k,l)*sin(k*x + l*y - dispersion(k,l)*t + phase1) 
                velocity(2) = velocity(2) - l*eps*amplitude(k,l)*sin(k*x + l*y - dispersion(k,l)*t + phase2)
            end do
        end do
                
    end function velocity_field
    
    function velocity_divergence(x,y,t,phase2) result(div)
        real :: x, y, t, k, l, phase2, div
        real, dimension(64) :: k_array, l_array
        integer :: counter_l, counter_k

        !Initialisation
        k_array = linspace(-2.0,2.0,64)
        l_array = linspace(-2.0,2.0,64)
        div = 0.
        phase2 = 0.0

        do counter_k = 1, 64
            k = k_array(counter_k)
            do counter_l = 1, 64
                l = l_array(counter_l)
                div = div - (k**2 + l**2)*amplitude(k,l)*cos(k*x + l*y - dispersion(k,l)*t + phase2)
            end do
        end do
    end function velocity_divergence

end module velocity