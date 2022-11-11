program rossby_wave_velocity
    use rossby_wave_attributes
    use admin
    implicit none

    print*,velocity_field(1.,0.,0.)
    call average_speed(5.0,5.0,64,0.)

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
                velocity(1) = velocity(1) + 2*pi*l*amplitude(k,l)*sin(2*pi*k*x + 2*pi*l*y - dispersion(k,l)*t + phase1)
                velocity(2) = velocity(2) - 2*pi*k*amplitude(k,l)*sin(2*pi*k*x + 2*pi*l*y - dispersion(k,l)*t + phase1)
            end do
        end do

    end function velocity_field

    subroutine average_speed(x_endpoint, y_endpoint, grid_resolution, t)
    !x_endpoint = 5.0, y_endpoint = 5.0, grid_resolution = 64
        real :: x_endpoint, y_endpoint, t, spatial_speed_average
        integer :: grid_resolution, counter_x, counter_y
        real, dimension(64) :: x_array, y_array
        real :: x,y

        !Initialisation
        x_array = linspace(-x_endpoint,x_endpoint,grid_resolution)
        y_array = linspace(-y_endpoint,y_endpoint,grid_resolution)
        spatial_speed_average = 0

        do counter_x = 1,grid_resolution
            x = x_array(counter_x)
            do counter_y = 1,grid_resolution
                y = y_array(counter_y)
                spatial_speed_average = spatial_speed_average + dot_product(velocity_field(x,y,t), velocity_field(x,y,t))
            end do
        end do
        spatial_speed_average = spatial_speed_average/(grid_resolution**2)
    end subroutine average_speed

end program rossby_wave_velocity