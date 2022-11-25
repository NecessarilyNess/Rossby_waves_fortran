program test
use velocity
use random
implicit none

real(dp) :: phase1(5,5), phase2(3,3), phase3(9,9)
integer :: num_waves
real(dp), allocatable :: k(:), l(:), phase(:), ring(:)

num_waves = 10

phase1(1,:) = (/.1_dp,.4_dp,.8_dp,.3_dp,0._dp/)
phase1(2,:) = (/.2_dp,.6_dp,.34_dp,.76_dp,.92_dp/)
phase1(3,:) = (/.04_dp,.37_dp,.08_dp,.5_dp,.8_dp/)
phase1(4,:) = (/1._dp,.28_dp,.31_dp,.93_dp,.64_dp/)
phase1(5,:) = (/.39_dp,.56_dp,.7_dp,.9_dp,.4_dp/)

phase2(1,:) = (/.1_dp,.4_dp,.8_dp/)
phase2(2,:) = (/.2_dp,.6_dp,.34_dp/)
phase2(3,:) = (/.04_dp,.37_dp,.08_dp/)

phase3(1,:) = (/.1_dp,.4_dp,.8_dp,.3_dp,0._dp,.31_dp,.41_dp,.59_dp,.26_dp/)
phase3(2,:) = (/.2_dp,.6_dp,.34_dp,.76_dp,.92_dp,.53_dp,.58_dp,.97_dp,.93_dp/)
phase3(3,:) = (/.04_dp,.37_dp,.08_dp,.5_dp,.8_dp,.23_dp,.84_dp,.62_dp,.64_dp/)
phase3(4,:) = (/1._dp,.28_dp,.31_dp,.93_dp,.64_dp,.33_dp,.83_dp,.27_dp,.95_dp/)
phase3(5,:) = (/.39_dp,.56_dp,.7_dp,.9_dp,.4_dp,.02_dp,.88_dp,.41_dp,.97_dp/)
phase3(6,:) = (/.16_dp,.93_dp,.99_dp,.37_dp,.51_dp,.05_dp,.82_dp,.09_dp,.74_dp/)
phase3(7,:) = (/.94_dp,.45_dp,.92_dp,.30_dp,.78_dp,.16_dp,.40_dp,.62_dp,.86_dp/)
phase3(8,:) = (/.20_dp,.89_dp,.98_dp,.62_dp,.80_dp,.34_dp,.82_dp,.53_dp,.42_dp/)
phase3(9,:) = (/.11_dp,.70_dp,.67_dp,.98_dp,.21_dp,.48_dp,.08_dp,.65_dp,.13_dp/)

allocate(k(num_waves))
allocate(l(num_waves))
allocate(phase(num_waves))
allocate(ring(num_waves))

!k = (/1._dp,-.2_dp,-1._dp/)
!l = (/0._dp,1._dp,-.7_dp/)
phase = random_numbers(num_waves)
ring = random_numbers(num_waves)
k = cos(2*pi*ring)
l = sin(2*pi*ring)
print*,phase
print*,ring
print*,k
print*,l

!print*,dispersion(1.,0.)
!call test_func(1.,1.,-1.,0.,0.5,-1.,21,phase1,phase2)
!call velocity_csv(21,0._dp,0._dp)
!print*,velocity_field((300./63.) -5.,(320./63.) -5.,0.,phase1,phase2)
!call waves_csv(1._dp,1._dp,64,0._dp)
!call streamfunction1_csv(21,-1._dp,-.7_dp,0._dp)
call streamfunctionn_csv(21,k,l,phase,num_waves)

contains

    subroutine test_func(k1,l1,k2,l2,k3,l3, timesteps, phase1, phase2)
        real(dp) :: x, y, k1, l1, k2, l2, k3, l3, t, phase1, phase2, velocity(2), t_array(timesteps)
        integer :: grid_resolution, counter_x, counter_y, counter_t, timesteps
        real(dp), dimension(64) :: x_array, y_array

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
    
    function sin_problem(k,l,x,y,t,phase) result(sin1)
        real(dp) :: x, y, t, k, l, phase, sin1
        !real :: x, y, t, k, l, phase, sin1

        !sin1 = sin(2*pi*k*x + 2*pi*l*y - dispersion(k,l)*t + phase)
        sin1 = sin(2*pi*k*x + 2*pi*l*y - dispersion(k,l)*t + phase)
                
    end function sin_problem
    
    function velocity_test(x,y,t,phase1,phase2) result(velocity)
        real(dp) :: x, y, t, k, l, phase1, phase2, velocity(2), k_array(4), l_array(4), phi1
        !real :: x, y, t, k, l, phase1, phase2, velocity(2)
        !real, dimension(64) :: k_array, l_array
        integer :: counter_l, counter_k

        !Initialisation
        k_array = linspace(-2.0,2.0,4)
        l_array = linspace(-2.0,2.0,4)
        velocity = (/0.,0./)

        do counter_k = 1, 4
            k = k_array(counter_k)
            do counter_l = 1, 4
                l = l_array(counter_l)
                !phi1 = 2*pi*phase1(counter_k,counter_l)
                velocity(1) = velocity(1) + 2*pi*l*(1-eps)*amplitude(k,l)*sin_problem(k,l,x,y,t,phase1)
                !velocity(1) = velocity(1) + 2*pi*l*amplitude(k,l)*sin_problem(k,l,x,y,t,phase1)
                velocity(1) = velocity(1) - 2*pi*k*eps*amplitude(k,l)*sin_problem(k,l,x,y,t,phase2)
                velocity(2) = velocity(2) - 2*pi*k*(1-eps)*amplitude(k,l)*sin_problem(k,l,x,y,t,phase1)
                !velocity(2) = velocity(2) - 2*pi*k*amplitude(k,l)*sin_problem(k,l,x,y,t,phase1)
                velocity(2) = velocity(2) - 2*pi*l*eps*amplitude(k,l)*sin_problem(k,l,x,y,t,phase2)
            end do
        end do
                
    end function velocity_test

    subroutine velocity_csv(timesteps, phase1, phase2)
        real(dp) :: x, y, t, phase1, phase2, velocity(2), t_array(timesteps), x_array(64), y_array(64)
        !real :: x, y, t, phase1, phase2, velocity(2), t_array(timesteps)
        !real, dimension(64) :: x_array, y_array
        integer :: grid_resolution, counter_x, counter_y, counter_t, timesteps

        !Initialisation
        x_array = linspace(-5.0,5.0,64)
        y_array = linspace(-5.0,5.0,64)
        t_array = linspace(0.,timesteps-1.,timesteps)

        open(1, file = 'data.csv')

        do counter_t = 1,timesteps
        print*,counter_t
        t = t_array(counter_t)
        do counter_y = 1,64
        y = y_array(counter_y)
        do counter_x = 1,64
        x = x_array(counter_x)
            velocity = velocity_test(x,y,t,phase1,phase2)
            write(1,*) velocity(1), ',', velocity(2)
        end do
        end do
        end do
        close(1) 
    end subroutine velocity_csv

    subroutine streamfunction1_csv(timesteps, k, l, phase1)
        real(dp) :: x, y, t, k, l, phase1, psi, t_array(timesteps), x_array(64), y_array(64)
        !real :: x, y, t, phase1, phase2, velocity(2), t_array(timesteps)
        !real, dimension(64) :: x_array, y_array
        integer :: grid_resolution, counter_x, counter_y, counter_t, timesteps

        !Initialisation
        x_array = linspace(-5.0,5.0,64)
        y_array = linspace(-5.0,5.0,64)
        t_array = linspace(0.,timesteps-1.,timesteps)

        open(1, file = 'streamfunction.csv')

        do counter_t = 1,timesteps
        print*,counter_t
        t = t_array(counter_t)
        do counter_x = 1,64
        x = x_array(counter_x)
        do counter_y = 1,64
        y = y_array(counter_y)
            psi = cos(2*pi*k*x + 2*pi*l*y - dispersion(k,l)*t + phase1)
            write(1,*) psi
        end do
        end do
        end do
        close(1) 
    end subroutine streamfunction1_csv

    subroutine streamfunctionn_csv(timesteps, k, l, phase1, n)
        real(dp) :: x, y, t, k(n), l(n), phase1(n), psi, t_array(timesteps), x_array(64), y_array(64)
        !real :: x, y, t, phase1, phase2, velocity(2), t_array(timesteps)
        !real, dimension(64) :: x_array, y_array
        integer :: grid_resolution, counter_x, counter_y, counter_t, counter_n, timesteps, n

        !Initialisation
        x_array = linspace(-5.0,5.0,64)
        y_array = linspace(-5.0,5.0,64)
        t_array = linspace(0.,timesteps-1.,timesteps)

        open(1, file = 'streamfunctionn.csv')

        do counter_t = 1,timesteps
        print*,counter_t
        t = t_array(counter_t)
        do counter_x = 1,64
        x = x_array(counter_x)
        do counter_y = 1,64
        y = y_array(counter_y)
        psi = 0.
        do counter_n = 1,n
            psi = psi + cos(2*pi*k(counter_n)*x + 2*pi*l(counter_n)*y - dispersion(k(counter_n),l(counter_n))*t + phase1(counter_n))
        end do
        write(1,*) psi
        end do
        end do
        end do
        close(1) 
    end subroutine streamfunctionn_csv

    subroutine waves_csv(x, y, num, phase)
        real(dp) :: x, y, phase, velocity(2), k_array(num), l_array(num), k, l, wave1, wave2
        !real :: x, y, t, phase1, phase2, velocity(2), t_array(timesteps)
        !real, dimension(64) :: x_array, y_array
        integer :: grid_resolution, counter_k, counter_l, num

        !Initialisation
        k_array = linspace(-2.,2.,num)
        l_array = linspace(-2.,2.,num)
        velocity = (/0.,0./)

        open(1, file = 'waves3.csv')
            write(1,*) x, ',', y, ',', phase, ',', num

        do counter_k = 1,num
        k = k_array(counter_k)
        do counter_l = 1,num
            l = l_array(counter_l)
            wave1 = 2*pi*l*amplitude(k,l)*sin(2*pi*k*x + 2*pi*l*y - dispersion(k,l)*0 + phase)
            wave2 = - 2*pi*k*amplitude(k,l)*sin(2*pi*k*x + 2*pi*l*y - dispersion(k,l)*0 + phase)
            velocity(1) = velocity(1) + 2*pi*l*amplitude(k,l)*sin(2*pi*k*x + 2*pi*l*y - dispersion(k,l)*0 + phase)
            velocity(2) = velocity(2) - 2*pi*k*amplitude(k,l)*sin(2*pi*k*x + 2*pi*l*y - dispersion(k,l)*0 + phase)
            write(1,*) k, ',', l, ',', wave1, ',', wave2
        end do
        end do
        write(1,*) velocity(1), ',', velocity(2)
        close(1) 
    end subroutine waves_csv

end program test