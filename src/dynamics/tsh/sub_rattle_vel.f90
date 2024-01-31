! Correct velocity for fixing the distance between atom pairs by using RATTLE
subroutine sub_rattle_vel(n_pair, atom_a, atom_b, dis, &
                            n_atom, coor_x, coor_y, coor_z, &
                                    vel_x, vel_y, vel_z, &
                                    mass, &
                                    tol, file_md_out, n_loop)
implicit none
integer, intent(in) :: n_atom, n_pair, n_loop
double precision, intent(in) :: tol

double precision, intent(in), dimension(n_pair) :: dis
double precision, intent(in), dimension(n_atom) :: mass
integer, intent(in), dimension(n_pair) :: atom_a, atom_b


double precision, intent(in), dimension(n_atom) :: coor_x, &
                                                   coor_y, &
                                                   coor_z
                        
double precision, intent(inout), dimension(n_atom) :: vel_x, &
                                                      vel_y, &
                                                      vel_z

integer :: i, j
integer :: a, b
logical :: loop 
double precision :: xx, yy, zz
double precision :: diff
double precision :: rma, rmb, constr
integer, dimension(n_pair) :: check
integer :: times, file_md_out, time1, time2
double precision :: dtime

loop = .true.
times = 0

write (file_md_out, *) "doing RATTLE for velocity..."
call system_clock(time1)
do while ( loop )
    times = times + 1
    do i=1, n_pair
        a = atom_a(i)
        b = atom_b(i)
        rma = 1.0 / mass(a)
        rmb = 1.0 / mass(b)
        xx = coor_x(a) - coor_x(b)
        yy = coor_y(a) - coor_y(b)
        zz = coor_z(a) - coor_z(b)
        diff = xx * (vel_x(a) - vel_x(b)) + yy * (vel_y(a) - vel_y(b)) + zz * (vel_z(a) - vel_z(b))

        if ( abs(diff) > tol ) then 
            constr = diff / ((rma + rmb) * dis(i))
            vel_x(a) = vel_x(a) - constr * xx * rma
            vel_y(a) = vel_y(a) - constr * yy * rma
            vel_z(a) = vel_z(a) - constr * zz * rma
            vel_x(b) = vel_x(b) + constr * xx * rmb
            vel_y(b) = vel_y(b) + constr * yy * rmb
            vel_z(b) = vel_z(b) + constr * zz * rmb
            check(i) = 0
        else
            check(i) = 1
        endif 
    enddo 
    
    if ( times .gt. n_loop ) then 
        exit
    endif 

    loop = .false.
    do j = 1, n_pair
        if ( check(j) .eq. 0 ) then 
            loop = .true.
            exit
        endif 
    enddo 
enddo

call system_clock(time2)
dtime = (time2 - time1) / 1000
write (file_md_out, *)  "correction of velocities"
write (file_md_out, *)  "times of loop : ", times, "Total time of correction of velocity : ", dtime, "seconds"
if ( times .gt. n_loop ) then 
    write (file_md_out, *) "Loop time have exceeded the limit. ", "(", n_loop, ")"
    close (file_md_out)
    stop
endif 
return 

end 

