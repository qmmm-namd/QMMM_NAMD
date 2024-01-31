subroutine sub_get_c_n_elec_2_traj_adjusted(x, p, gammas, c_n)
use mod_main, only : label_debug, n_state, ggamma
implicit none

! --- arguments ---
real(kind=8) x(n_state)
real(kind=8) p(n_state)
real(kind=8) gammas(n_state)
real(kind=8) c_n(n_state, n_state)

! --- local variables ---
integer i, j


do i = 1, n_state
  c_n(i,i) = 0.5 * (x(i)**2 + p(i)**2 - gammas(i) )
enddo

do i = 1, n_state-1
  do j = i+1, n_state
    c_n(i,j) = 0.5 * ( x(i) * x(j) + p(i) * p(j) )
    c_n(j,i) = c_n(i,j)
  enddo
enddo

if (label_debug >= 2)then
open(10, file='c_n.dat', access='append')
do i = 1, n_state
  write(10, *) c_n(i,:)
enddo
close(10)
endif

return
end
