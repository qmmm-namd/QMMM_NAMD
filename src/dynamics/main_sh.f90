       Program main

   
       implicit none
       include 'param.def'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  - Define all input files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       integer :: file_ini_coor,  file_ini_vel
              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
!!!!!!!!-- Define all variables to read the geometry----
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      integer :: n_atom
      double precision, allocatable, dimension(:) ::  coor_x, &
                                                      coor_y, &
                                                      coor_z
      double precision, allocatable, dimension(:) ::  vel_x, &
                                                      vel_y, &
                                                      vel_z

      character*2, allocatable, dimension(:)      ::  atom_label




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  - Define all input files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       integer :: file_dynvar_in


       integer :: i, j, it, i_state
       integer :: ntime, n_sav_stat, n_sav_traj
       integer :: n_state, ntime_ele
       
       integer :: qm_method, dyn_method, qm_package, label_ZN
     
       integer, allocatable, dimension(:) ::  md_state
       character (len=256) md_state_list
       
      

       double precision :: time,  dtime, cor_dec
       integer :: seed_random
 
       integer :: label_read_velocity
       integer :: label_nac_phase      
       integer :: label_reject_hops
       integer :: label_restart

       !!!!!!!!!!!!!!!! qmmm, rattle, bomd
       integer :: label_qmmm=0, label_rattle=0, label_qmmm_bomd=0
       character (len=256) :: python
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       double precision ::   hop_e   
        
       complex (kind=8), allocatable, dimension(:, :) :: rho
       integer :: index_state

       namelist /control/ dyn_method, label_ZN, ntime, dtime, & 
            n_sav_stat, n_sav_traj, ntime_ele, &
            qm_method, n_state, md_state_list, i_state, &
            seed_random, cor_dec, label_nac_phase, label_reject_hops, &
            hop_e, label_read_velocity, label_restart, &
            label_qmmm, label_rattle, python


       double precision :: gamma0, temperature
       namelist /langevin/ gamma0, temperature
    
       !!!!!!!!!!!!!!!!!!!!!!!!!! qmmm and shake
       integer :: n_atom_qm, n_atom_fro, n_pair, &
              qm_method_qmmm, mm_method_qmmm, qmmm_nac, label_LA=0, CS_scheme, &
              qmmm_scheme, n_link_atom, n_tot_atom, n_qm_la, n_loop_shake=1000, &
              qmmm_label_qm
       double precision :: rattle_tol=10E-9, LA_dis
       character*2 :: LA_ele

       namelist /qmmm/ qm_method_qmmm, mm_method_qmmm, label_qmmm_bomd, qmmm_nac, &
              label_LA, LA_ele, LA_dis, CS_scheme, qmmm_scheme, n_link_atom, qmmm_label_qm
       namelist /shake/ rattle_tol, n_loop_shake
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       write (*,*) "------------------------"
       write (*,*) "JADE (version 1.0 alpha)"
       write (*,*) "========================"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!   ---  Initialize the parameters for dynamics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

       write (*,*) "Initialization"
       file_dynvar_in=21       
       open(unit=file_dynvar_in, file="dyn.inp")      
       write (*,*) "dyn.inp open to read.."
       read(file_dynvar_in, nml = control)

       if (dyn_method == 201) then
          write(*, *) "Langevin Method is used to Control Temperature"
          read(file_dynvar_in, nml = langevin)
       endif       

       close(file_dynvar_in)

       allocate (md_state(n_state))       
       read(md_state_list, *) (md_state(i),i=1,n_state)
       write(*,nml=control)

       write(*,nml=langevin)
       dtime = dtime/TOFS


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
       ! read shake and qmmm      
        
       if (label_rattle .eq. 1) then 
              open(unit=file_dynvar_in, file="dyn.inp")
              read(file_dynvar_in, nml = shake)
              close(file_dynvar_in)

              write(*, *) "Rattle method is used to fix bond distance"
              write(*, nml=shake)

            open(unit=1111, file='atom_pair')
            read(1111, *) n_pair
            close(1111)

            endif 

       
       if (label_qmmm .eq. 1) then 
              open(unit=file_dynvar_in, file="dyn.inp")
              read(file_dynvar_in, nml = qmmm)
              close(file_dynvar_in)

              write(*, *) "Electronic calculation is used by QM/MM"
              write(*, nml=qmmm)

            open(unit=1111, file='qmmm_index')
            read(1111, *) n_atom_qm, n_atom_fro
            close(1111)

       endif 
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       write (*,*) " Finish  the initialization"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  -   Find how many atoms !!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     if (dyn_method /= 11) then

       file_ini_coor=11
       open(unit=file_ini_coor, file="stru_xyz.in")
       read (file_ini_coor, *) n_atom
       close(11)
       


       allocate (coor_x(n_atom))
       allocate (coor_y(n_atom))
       allocate (coor_z(n_atom))
       allocate (atom_label(n_atom))

       allocate (vel_x(n_atom))
       allocate (vel_y(n_atom))
       allocate (vel_z(n_atom))
       allocate (rho(n_state, n_state))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! --   Initialize the coordinates  and momentum !!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      
      
        file_ini_coor=11
        open(unit=file_ini_coor, file="stru_xyz.in")
        call sub_read_coor ( n_atom, coor_x, coor_y, coor_z, &
                         atom_label, file_ini_coor)
        close(11)
       !write(*,*)"atom_label(1)  ", atom_label(1)
      

       vel_x =0.d0
       vel_y =0.d0
       vel_z =0.d0

!     Read  the initial velocity if necessary
      if (label_read_velocity .eq. 1) then
          file_ini_vel=12
          open(unit=file_ini_vel, file="vel_xyz.in")
          call sub_read_ini_vel (n_atom, vel_x, vel_y, vel_z, &
                             atom_label, file_ini_vel)
          close(file_ini_vel)
      endif      

   
      if (label_restart .eq. 0 ) then
         do i=1, n_state
            if (md_state(i)  .eq.  i_state) then
               index_state=i
            endif
         enddo

         do i=1, n_state
            do j=1, n_state
               rho(i,j)=0.d0
            enddo
         enddo

         rho(index_state, index_state) = 1.d0


         call system ("rm -rf *.dat")
         call system ("rm -rf *.out")
         call system ("rm -rf QC_TMP")
!         call system ("rm -rf ZN QC_TMP")

       ! remove old wfu, in molpro calculation of QM/MM
        call system ("rm -rf old_wfu")
       ! remove qmmm work directory and temporary directories
        call system ("rm -rf QMMM_TMP")
        call system ("rm -rf old_fort11")
        call system ("rm -rf old_chk")
        ! remove json files
        call system ("rm -rf *.json")
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      endif

   


      if (label_restart .eq. 1 ) then
 
          call sub_read_restart (n_atom, &
                             coor_x, coor_y, coor_z, &
                             vel_x, vel_y, vel_z, &
                             atom_label, &
                             n_state, &
                             rho, &
                             index_state, &
                             it, time )

      endif



       it = 0
       time = 0
       !write(*,*)"atom_label(1)  ", atom_label(1)
       call  sub_write_current_geom (n_atom, &
                            coor_x, coor_y, coor_z, &
                            vel_x, vel_y, vel_z, &
                            atom_label, &
                            n_state, &
                            rho, &
                            index_state, & 
                            it, time )

       deallocate (coor_x)
       deallocate (coor_y)
       deallocate (coor_z)
       deallocate (atom_label)

       deallocate (vel_x)
       deallocate (vel_y)
       deallocate (vel_z)
!     endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!   - Dynamics Methods !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       



if (label_qmmm .ne. 1) then 
       if ( (n_state .gt. 1)  .and. &
            ( dyn_method .eq. 1 ) )  then
              write (*,*) "Surface Hopping with analytical NAC"
          
              if (label_ZN .eq. 0) then 
                     write (*,*) "Surface Hopping with tully method"
                     call  sub_sh_ana_nac   ( n_atom, &
                                   i_state, & 
                                   n_state, md_state, &
                                   ntime, dtime, ntime_ele, &
                                   n_sav_stat, n_sav_traj, &
                                   label_ZN, &
                                   qm_method,&
                                   seed_random, &
                                   cor_dec, &
                                   label_nac_phase, &
                                   label_reject_hops, &
                                   hop_e )
              
              else
                     write (*,*) "Surface Hopping with ZN method"
                     write (*,*) "hop_e", hop_e
                     call  sub_sh_zn        ( n_atom, &
                                      i_state, & 
                                      n_state, md_state, &
                                      ntime, dtime, ntime_ele, &
                                      n_sav_stat, n_sav_traj, &
                                      label_ZN, &
                                      qm_method,&
                                      seed_random, &
                                      label_reject_hops, &
                                      hop_e )
              endif

       endif

       
       if ( (n_state .gt. 1)  .and. &
            (dyn_method .eq. 2 ) )  then
              write (*,*) "Surface Hopping with numerical NAC"
              write (*,*) "RANDOM", seed_random
              call  sub_sh_num_nac  ( n_atom, &
                                   i_state, & 
                                   n_state, md_state, &
                                   ntime, dtime, ntime_ele, &
                                   n_sav_stat, n_sav_traj, &
                                   qm_method,&
                                   seed_random, &
                                   cor_dec, &
                                   label_nac_phase, &
                                   label_reject_hops, &
                                   hop_e )
       endif

       if ( (n_state .gt. 1)  .and. &
            (dyn_method .eq. 201 ) )  then
              write (*,*) "Surface Hopping with numerical NAC and langevin thermostat"
              write (*,*) "RANDOM", seed_random
              call  sub_lang_sh_num_nac  ( n_atom, &
                                   i_state, &                           
                                   n_state, md_state, &
                                   ntime, dtime, ntime_ele, &
                                   n_sav_stat, n_sav_traj, &
                                   qm_method,&
                                   seed_random, &
                                   cor_dec, &
                                   label_nac_phase, &
                                   label_reject_hops, &
                                   hop_e, &
                                   gamma0, &
                                   temperature )
       endif


       
       if ( (n_state .gt. 1)  .and. &
            (dyn_method .eq. 3 ) )  then
              write (*,*) "ab initio dynamics at single surface (BOMD)."
!       call  takasuka          ( n_atom, &
!                                 coor_x, coor_y, coor_z, &
!                                 vel_x, vel_y, vel_z, &
!                                 atom_label, &
!                                 n_states, md_state )
       endif


       if (  dyn_method .eq. 5  )  then
              write (*,*) "Surface Hopping with SOC, only for first spin excited states."
              call  sub_sh_ana_nac ( n_atom, &
                                   i_state, &                                 
                                   n_state, md_state, &
                                   ntime, dtime, ntime_ele, &
                                   n_sav_stat, n_sav_traj, &
                                   qm_method,&
                                   seed_random, &
                                   cor_dec, &
                                   label_nac_phase, &
                                   label_reject_hops, &
                                   hop_e )
       endif

       if ( (n_state .gt. 1)  .and. &
            ( dyn_method .eq. 11 ) )  then
              write (*,*) "Mapping dynamics with Meyer-Miller Model."
              write (*,*) "J. Chem. Phys. 147, 064112 (2017)."
              call  sub_mapping_ana_nac()

       endif
      

! QMMM calculation
else
       ! specify total atom with link atom 
       if (label_LA == 1) then 
              n_tot_atom = n_atom + n_link_atom
              n_qm_la = n_atom_qm + n_link_atom
       else 
              n_tot_atom = n_atom 
              n_qm_la = n_atom_qm
       endif 
          
       if (label_ZN .eq. 0) then 
              write (*,*) "Surface Hopping with tully method for QM/MM"
              call sub_sh_ana_nac_qmmm ( n_atom, &
                                   i_state, & 
                                   n_state, md_state, &
                                   ntime, dtime, ntime_ele, &
                                   n_sav_stat, n_sav_traj, &
                                   label_ZN, &
                                   qm_method,&
                                   seed_random, &
                                   cor_dec, &
                                   label_nac_phase, &
                                   label_reject_hops, &
                                   hop_e, &
                                   label_rattle, rattle_tol, n_pair, &
                                   n_atom_qm, n_atom_fro, label_qmmm_bomd, n_loop_shake, &
                                   label_LA, n_link_atom, n_tot_atom, n_qm_la)

       else 
              write (*,*) "Surface Hopping with ZN method for QM/MM"
              write (*,*) "hop_e", hop_e
              call  sub_sh_zn_qmmm ( n_atom, &
                                   i_state, & 
                                   n_state, md_state, &
                                   ntime, dtime, ntime_ele, &
                                   n_sav_stat, n_sav_traj, &
                                   label_ZN, &
                                   qm_method,&
                                   seed_random, &
                                   label_reject_hops, &
                                   hop_e, & 
                                   label_rattle, rattle_tol, n_pair, &
                                   n_atom_qm, n_atom_fro, label_qmmm_bomd, n_loop_shake, &
                                   label_LA, n_link_atom, n_tot_atom, n_qm_la)
       endif 

endif


       deallocate (md_state)


       end 
