! Correct coordinate for fixing the distance between atom pairs by using RATTLE
subroutine sub_rattle_coor(n_pair, atom_a, atom_b, dis, &
                            n_atom, coor_x, coor_y, coor_z, &
                            old_coor_x, old_coor_y, old_coor_z, &
                            vel_x, vel_y, vel_z, &
                            mass, dt, &
                            tol, file_md_out, n_loop)
implicit none
integer, intent(in) :: n_atom, n_pair, n_loop
double precision, intent(in) :: tol, dt

double precision, intent(in), dimension(n_pair) :: dis
double precision, intent(in), dimension(n_atom) :: mass
integer, intent(in), dimension(n_pair) :: atom_a, atom_b


double precision, intent(inout), dimension(n_atom) :: coor_x, &
                                                      coor_y, &
                                                      coor_z, &
                                                      vel_x, &
                                                      vel_y, &
                                                      vel_z

double precision, intent(in), dimension(n_atom) :: old_coor_x, &
                                                   old_coor_y, &
                                                   old_coor_z

integer :: i, j
integer :: a, b
logical :: loop 
double precision :: xx, yy, zz
double precision, dimension(n_pair) :: old_xx, old_yy, old_zz
double precision :: diff
double precision :: rma, rmb, constr
integer, dimension(n_pair) :: check

integer :: times, file_md_out, time1, time2
double precision :: dtime 
loop = .true.
times = 0

write (file_md_out, *) "doing RATTLE for coordinate..."
call system_clock(time1)

do i=1, n_pair
    a = atom_a(i)
    b = atom_b(i)
    old_xx(i) = old_coor_x(a) - old_coor_x(b)
    old_yy(i) = old_coor_y(a) - old_coor_y(b)
    old_zz(i) = old_coor_z(a) - old_coor_z(b)
enddo

do while ( loop )
    times = times + 1
    ! write(file_md_out, *) "current times of loop : ", times
    do i=1, n_pair
        a = atom_a(i)
        b = atom_b(i)
        rma = 1.0 / mass(a)
        rmb = 1.0 / mass(b)
        xx = coor_x(a) - coor_x(b)
        yy = coor_y(a) - coor_y(b)
        zz = coor_z(a) - coor_z(b)
        diff = xx * xx + yy * yy + zz * zz - dis(i)

        if ( abs(diff) > tol) then 
            ! write(file_md_out, *) times, a, b, diff, (diff + dis(i)) ** 0.5 
            constr = diff / (2.0 * dt * &
                (old_xx(i) * xx + old_yy(i) * yy + old_zz(i) * zz) * (rma + rmb))
            coor_x(a) = coor_x(a) - constr * dt * old_xx(i) * rma
            coor_y(a) = coor_y(a) - constr * dt * old_yy(i) * rma 
            coor_z(a) = coor_z(a) - constr * dt * old_zz(i) * rma 
            coor_x(b) = coor_x(b) + constr * dt * old_xx(i) * rmb
            coor_y(b) = coor_y(b) + constr * dt * old_yy(i) * rmb 
            coor_z(b) = coor_z(b) + constr * dt * old_zz(i) * rmb 
            vel_x(a) = vel_x(a) - constr * old_xx(i) * rma
            vel_y(a) = vel_y(a) - constr * old_yy(i) * rma
            vel_z(a) = vel_z(a) - constr * old_zz(i) * rma
            vel_x(b) = vel_x(b) + constr * old_xx(i) * rmb
            vel_y(b) = vel_y(b) + constr * old_yy(i) * rmb
            vel_z(b) = vel_z(b) + constr * old_zz(i) * rmb
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
write (file_md_out, *)  "correction of coordinates"
write (file_md_out, *)  "times of loop : ", times, "Total time : ", dtime, "seconds"
if ( times .gt. n_loop ) then 
    write (file_md_out, *) "Loop time have exceeded the limit. ", "(", n_loop, ")"
    close (file_md_out)
    stop
endif 
return 

end subroutine sub_rattle_coor


