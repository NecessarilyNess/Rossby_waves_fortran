function amplitude_func(k, l) result(amp)
    real :: k, l, amp
    !call get_command_argument(1, k)
    !if (len_trim(k) == 0) exit
    !read(k1,'(f9.0)')k
    !print*,k

    !call get_command_argument(2, l)
    !if (len_trim(l1) == 0) exit
    !read(l1,'(f9.0)')l
    !print*,l

    amp = exp(-k**2 - l**2) * (k**2 + l**2)
 
end function amplitude_func