subroutine sub_get_dia_sub_get_dia_V_elec_elec(x_mode, p_mode, &
                                   freq,k_tun, exci_e, &
                                   lambda0, lambda1, &
                                   V)
use mod_main
implicit none                            

! --- arguments ---
real(kind=8) exci_e(n_state)
real(kind=8), dimension(n_state, n_state) :: V
real(kind=8), dimension(n_mode_HO, n_state, n_state) :: lambda0, lambda1
real(kind=8), dimension(n_mode_HO) :: x_mode, &
                                   p_mode
real(kind=8) k_tun(n_state,n_mode_HO)
real(kind=8) freq(n_state,n_mode_HO)

! -- local variables ---
integer i, j, k
real(kind=8) tmp, x_tmp, p_tmp, freq_tmp

V = 0.d0

do i = 1, n_state
  tmp = exci_e(i)
!  write(*,*) i, tmp
  do j = 1, n_mode_HO
    x_tmp = x_mode(j)
    freq_tmp = freq(i,j)
    tmp = tmp + 0.5d0 * freq_tmp * x_tmp * x_tmp &
              + k_tun(i,j) * x_tmp
!    write(*,*) x_tmp, p_tmp, freq_tmp, k_tun(i,j)
  enddo
  V(i,i) =  tmp
enddo

do i = 1, n_state
  do j = 1, n_state
    if(i == j) cycle
    do k = 1, n_mode_HO
      V(i, j) = V(i, j) + lambda0(k, i, j) + lambda1(k, i, j) * x_mode(k)
    enddo
  enddo
enddo

!open(10, file='E1_E2.dat', access='append')
!write(10,*) V(1,1), V(2,2)
!close(10)

return
end
