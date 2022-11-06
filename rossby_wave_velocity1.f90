program rossby_wave_velocity
    use rossby_wave_attributes
    use admin
    implicit none

    print*, velocity_field(0.0,0.0,0.0)

    contains 
    function velocity_field(x,y,t) result(velocity)
        !Returns the velocity field at point (x,y) at time t
        !upper and lower k,l are -2 and 2
        real :: x,y,t
        real :: k,l,phase1
        real, dimension(64) :: k_array, l_array
        integer :: counter_l, counter_k
        real, dimension(2) :: velocity

        !Initialisation
        k_array = linspace(-2.0,2.0,64)
        l_array = linspace(-2.0,2.0,64)
        velocity = (/0.0,0.0/)
        phase1 = 0.0

        do counter_k = 1, 64
            k = k_array(counter_k)
            do counter_l = 1, 64
                l = l_array(counter_l)
                velocity(1) = velocity(1) + l*amplitude(k,l)*sin(k*x + l*y - dispersion(k,l)*t + phase1)
                velocity(2) = velocity(2) - k*amplitude(k,l)*sin(k*x + l*y - dispersion(k,l)*t + phase1)
            end do
        end do

    end function velocity_field
end program rossby_wave_velocity