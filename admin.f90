module admin
    !use velocity
    implicit none

    contains
    function linspace(start_point, end_point, num_points) result(points_array)
        real :: start_point, end_point, increment, new_point
        integer :: num_points, counter
        real, dimension(:), allocatable :: points_array
        
        allocate(points_array(num_points))
        increment = (end_point - start_point)/(num_points-1)
        new_point = start_point
        do counter = 1, num_points
            points_array(counter) = new_point
            new_point = new_point + increment
        end do

    end function linspace
end module admin