subroutine sub_calc_dV_eff_div_dx_nuc_no_aver_PES_adia(dV_eff_div_dx_nuc, &
                           x_elec, p_elec, &
                           grad)
use mod_main
implicit none
! --- arguments ---
real(kind=8), dimension(n_mode) :: dV_eff_div_dx_nuc
real(kind=8), dimension(n_state, n_mode) :: grad
real(kind=8), dimension(n_state) :: x_elec, p_elec

! --- local variables ---
integer i_mode, i_state, j_state
real(kind=8) dtmp, dtmp2, dn_state_inv
real(kind=8) x_elec_tmp_i, p_elec_tmp_i
real(kind=8) x_elec_tmp_j, p_elec_tmp_j
real(kind=8) grad_i, grad_j
real(kind=8), dimension(n_state,n_state) :: c_n


dn_state_inv = 1.d0/dble(n_state)

call sub_get_c_n_elec_2(x_elec, p_elec, c_n)

! 1/n_state * sum(grad_i)
do i_mode = 1, n_mode
  dtmp = 0.d0
  do i_state = 1, n_state
    dtmp = dtmp + grad(i_state,i_mode) * c_n(i_state, i_state)
  enddo
  dV_eff_div_dx_nuc(i_mode) = dtmp
enddo


return
end