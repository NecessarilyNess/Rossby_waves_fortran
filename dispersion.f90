function dispersion_func(k, l) result(omega)
    real :: k, l, omega, beta, Rd
    beta = 2E-11
    Rd = 1E5

    omega = -beta * k/(k**2 + l**2 + Rd**(-2))
 
end function dispersion_func