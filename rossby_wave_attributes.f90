module rossby_wave_attributes
use constants
implicit none
integer,parameter :: dp = selected_real_kind(p=15,r=200)

contains
    function amplitude(k, l) result(amp)
        real(dp) :: k, l, amp

        amp = exp(-k**2 - l**2) * (k**2 + l**2)
    
    end function amplitude

    function dispersion(k, l) result(omega)
        real(dp) :: k, l, omega

        omega = -beta * k/(k**2 + l**2 + Rd**(-2))
 
    end function dispersion
end module rossby_wave_attributes