program test
use velocity
implicit none

real :: x,y,t,phase1(5,5),phase2

!x = 1.
!y = 0.
!t = 0.
phase1(1,:) = (/.1,.4,.8,.3,0./)
phase1(2,:) = (/.2,.6,.34,.76,.92/)
phase1(3,:) = (/.04,.37,.08,.5,.8/)
phase1(4,:) = (/1.,.28,.31,.93,.64/)
phase1(5,:) = (/.39,.56,.7,.9,.4/)
phase2 = 0.

!print*,dispersion(1.,0.)
!call test_func(1.,1.,-1.,0.,0.5,-1.,21,phase1,phase2)
call velocity_csv(81,phase1,phase2)
!print*,velocity_field((300./63.) -5.,(320./63.) -5.,0.,phase1,phase2)

contains

    subroutine test_func(k1,l1,k2,l2,k3,l3, timesteps, phase1, phase2)
    !x_endpoint = 5.0, y_endpoint = 5.0, grid_resolution = 64
        real :: k1, l1, k2, l2, k3, l3, t, phase1, phase2, velocity(2), t_array(timesteps)
        integer :: grid_resolution, counter_x, counter_y, counter_t, timesteps
        real, dimension(64) :: x_array, y_array
        real :: x,y

        !Initialisation
        x_array = linspace(-5.0,5.0,64)
        y_array = linspace(-5.0,5.0,64)
        t_array = linspace(0.,20.,timesteps)
        print*,t_array

        open(1, file = 'data.csv')

        do counter_t = 1,timesteps
        print*,counter_t
        t = t_array(counter_t)
        do counter_x = 1,64
        x = x_array(counter_x)
        do counter_y = 1,64
            y = y_array(counter_y)
            velocity(1) = l1*(1-eps)*amplitude(k1,l1)*sin(2*pi*k1*x + 2*pi*l1*y - dispersion(k1,l1)*t + 2*pi*phase1) 
            velocity(1) = velocity(1) + l2*(1-eps)*amplitude(k2,l2)*sin(2*pi*k2*x + 2*pi*l2*y - dispersion(k2,l2)*t + 2*pi*phase1) 
            velocity(1) = velocity(1) + l3*(1-eps)*amplitude(k3,l3)*sin(2*pi*k3*x + 2*pi*l3*y - dispersion(k3,l3)*t + 2*pi*phase1) 
            !velocity(1) = velocity(1)- k*eps*amplitude(k,l)*sin(k*x + l*y - dispersion(k,l)*t + phase2)
            velocity(2) = - k1*(1-eps)*amplitude(k1,l1)*sin(2*pi*k1*x + 2*pi*l1*y - dispersion(k1,l1)*t + 2*pi*phase1) 
            velocity(2) = velocity(2) - k2*(1-eps)*amplitude(k2,l2)*sin(2*pi*k2*x + 2*pi*l2*y - dispersion(k2,l2)*t + 2*pi*phase1) 
            velocity(2) = velocity(2) - k3*(1-eps)*amplitude(k3,l3)*sin(2*pi*k3*x + 2*pi*l3*y - dispersion(k3,l3)*t + 2*pi*phase1) 
            !velocity(2) = velocity(2) - l*eps*amplitude(k,l)*sin(k*x + l*y - dispersion(k,l)*t + phase2)
            write(1,*) velocity(1), ',', velocity(2)
        end do
        end do
        end do
        close(1) 
    end subroutine test_func
    
    function velocity_test(x,y,t,phase1,phase2) result(velocity)
        real :: x, y, t, k, l, phase1(5,5), phase2, velocity(2), phi1
        real, dimension(5) :: k_array, l_array
        integer :: counter_l, counter_k

        !Initialisation
        k_array = linspace(-2.0,2.0,5)
        l_array = linspace(-2.0,2.0,5)
        velocity = (/0.,0./)
        phase1 = 0.0
        phase2 = 0.0

        do counter_k = 1, 5
            k = k_array(counter_k)
            do counter_l = 1, 5
                l = l_array(counter_l)
                phi1 = 2*pi*phase1(counter_k,counter_l)
                velocity(1) = velocity(1) + 2*pi*l*(1-eps)*amplitude(k,l)*sin(2*pi*k*x + 2*pi*l*y - dispersion(k,l)*t + phi1) 
                velocity(1) = velocity(1) - 2*pi*k*eps*amplitude(k,l)*sin(2*pi*k*x + 2*pi*l*y - dispersion(k,l)*t + phase2)
                velocity(2) = velocity(2) - 2*pi*k*(1-eps)*amplitude(k,l)*sin(2*pi*k*x + 2*pi*l*y - dispersion(k,l)*t + phi1) 
                velocity(2) = velocity(2) - 2*pi*l*eps*amplitude(k,l)*sin(2*pi*k*x + 2*pi*l*y - dispersion(k,l)*t + phase2)
            end do
        end do
                
    end function velocity_test

    subroutine velocity_csv(timesteps, phase1, phase2)
    !x_endpoint = 5.0, y_endpoint = 5.0, grid_resolution = 64
        real :: k1, l1, k2, l2, k3, l3, t, phase1(5,5), phase2, velocity(2), t_array(timesteps)
        integer :: grid_resolution, counter_x, counter_y, counter_t, timesteps
        real, dimension(64) :: x_array, y_array
        real :: x,y

        !Initialisation
        x_array = linspace(-5.0,5.0,64)
        y_array = linspace(-5.0,5.0,64)
        t_array = linspace(0.,timesteps-1.,timesteps)
        print*,t_array

        open(1, file = 'data.csv')

        do counter_t = 1,timesteps
        print*,counter_t
        t = t_array(counter_t)
        do counter_x = 1,64
        x = x_array(counter_x)
        do counter_y = 1,64
            y = y_array(counter_y)
            velocity = velocity_test(x,y,t,phase1,phase2)
            write(1,*) velocity(1), ',', velocity(2)
        end do
        end do
        end do
        close(1) 
    end subroutine velocity_csv

end program test