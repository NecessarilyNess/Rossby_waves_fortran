program average_velocity
 implicit none
 integer :: k, l, x, y, a, t
 real :: u, v, x1, y1, k1, l1, phase1, amplitude_func, dispersion_func, s, b, t1
 
 a = 63
 b = 63.

 u = 0.
 v = 0.
 s = 0.
 !t = 5E11
 phase1 = 0
 do t=0,50
 do x=0,a
 do y=0,a
 do k=0,a
 do l=0,a
    x1 = (10*x/b - 5.)
    y1 = (10*y/b - 5.)
    k1 = (4*k/b - 2.)
    l1 = (4*l/b - 2.)
    t1 = t*1E8
    !print*,'l1 is',l1,'k1 is',k1,'vel is',amplitude_func(k1,l1)
    u = u + l1*amplitude_func(k1,l1)*sin(k1*x1 + l1*y1 - dispersion_func(k1,l1)*t1 + phase1)
    v = v - k1*amplitude_func(k1,l1)*sin(k1*x1 + l1*y1 - dispersion_func(k1,l1)*t1 + phase1)
    
 enddo
 enddo
    s = s + (u**2 + v**2)**(0.5)
    u = 0.
    v = 0.

 enddo
 enddo
 !print*,'time =',t1
 print*,s/((b+1)**2),','
 s = 0.
 
 !print*,'Average speed is',s/((b+1)**2)
 enddo
end program average_velocity