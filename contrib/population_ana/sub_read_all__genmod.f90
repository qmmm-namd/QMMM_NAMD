        !COMPILER-GENERATED INTERFACE MODULE: Thu Jul 16 14:19:50 2015
        MODULE SUB_READ_ALL__genmod
          INTERFACE 
            SUBROUTINE SUB_READ_ALL(N_AO,N_ATOM,N_STATE,TYPE_INPUT,CI_1,&
     &BASIS,S_AO_TO_MO,S_AO_OVERLAP)
              INTEGER(KIND=4), INTENT(IN) :: N_STATE
              INTEGER(KIND=4), INTENT(IN) :: N_ATOM
              INTEGER(KIND=4), INTENT(IN) :: N_AO
              INTEGER(KIND=4), INTENT(IN) :: TYPE_INPUT
              REAL(KIND=8), INTENT(INOUT) :: CI_1(N_STATE,N_AO,N_AO)
              INTEGER(KIND=4), INTENT(INOUT) :: BASIS(N_ATOM)
              REAL(KIND=8), INTENT(INOUT) :: S_AO_TO_MO(N_AO,N_AO)
              REAL(KIND=8), INTENT(INOUT) :: S_AO_OVERLAP(N_AO,N_AO)
            END SUBROUTINE SUB_READ_ALL
          END INTERFACE 
        END MODULE SUB_READ_ALL__genmod
