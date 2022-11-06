module velocity
use rossby_wave_attributes
implicit none
real :: eps = 0.01 !Placeholder value

contains
    function streamfunction(x,y,t,phase1) result(psi)
        real :: x, y, t, phase1, psi
                
        psi = amplitude(k1,l1)*cos(k1*x1 + l1*y1 - dispersion(k1,l1)*t1 + phase1)
    end function streamfunction
    
    function potential(x,y,t,phase2) result(phi)
        real :: x, y, t, phase2, phi
                
        phi = amplitude(k1,l1)*cos(k1*x1 + l1*y1 - dispersion(k1,l1)*t1 + phase2)
    end function potential
    
    function velocity_field(x,y,t,phase1,phase2) result(v)
        real :: x, y, t, phase1, phase2, v(2)
        v = (0,0)
                
        v(1) = l1*(1-eps)*amplitude(k1,l1)*sin(k1*x1 + l1*y1 - dispersion(k1,l1)*t1 + phase1) - k1*eps*amplitude(k1,l1)*sin(k1*x1 + l1*y1 - dispersion(k1,l1)*t1 + phase2)
        v(2) = -k1*(1-eps)*amplitude(k1,l1)*sin(k1*x1 + l1*y1 - dispersion(k1,l1)*t1 + phase1) - l1*eps*amplitude(k1,l1)*sin(k1*x1 + l1*y1 - dispersion(k1,l1)*t1 + phase2)
    end function velocity_field
    
    function velocity_divergence(x,y,t,phase2) result(div)
        real :: x, y, t, phase2, div
                
        div = -(k1**2 + l1**2)*amplitude(k1,l1)*cos(k1*x1 + l1*y1 - dispersion(k1,l1)*t1 + phase2)
    end function velocity_divergence

end module velocity