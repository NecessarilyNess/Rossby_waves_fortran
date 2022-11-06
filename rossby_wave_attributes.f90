module rossby_wave_attributes
use constants
implicit none

contains
    function amplitude_func(k, l) result(amp)
        real :: k, l, amp

        amp = exp(-k**2 - l**2) * (k**2 + l**2)
    
    end function amplitude_func

    function dispersion_func(k, l) result(omega)
        real :: k, l, omega

        omega = -beta * k/(k**2 + l**2 + Rd**(-2))
 
    end function dispersion_func
end module rossby_wave_attributes