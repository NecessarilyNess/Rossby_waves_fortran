module random
    use rossby_wave_attributes
    implicit none

    contains
    function random_numbers(num) result(r)
        real(dp), allocatable :: r(:)
        integer :: num

        allocate(r(num))
        call random_number(r)
    end function random_numbers
end module random