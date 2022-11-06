module rossby_wave_attributes
use constants
implicit none

contains
    function amplitude(k, l) result(amp)
        real :: k, l, amp

        amp = exp(-k**2 - l**2) * (k**2 + l**2)
    
    end function amplitude

    function dispersion(k, l) result(omega)
        real :: k, l, omega

        omega = -beta * k/(k**2 + l**2 + Rd**(-2))
 
    end function dispersion
end module rossby_wave_attributes