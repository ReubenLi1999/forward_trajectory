        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar 12 21:45:20 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE IAU_PN06__genmod
          INTERFACE 
            SUBROUTINE IAU_PN06(DATE1,DATE2,DPSI,DEPS,EPSA,RB,RP,RBP,RN,&
     &RBPN)
              REAL(KIND=8) :: DATE1
              REAL(KIND=8) :: DATE2
              REAL(KIND=8) :: DPSI
              REAL(KIND=8) :: DEPS
              REAL(KIND=8) :: EPSA
              REAL(KIND=8) :: RB(3,3)
              REAL(KIND=8) :: RP(3,3)
              REAL(KIND=8) :: RBP(3,3)
              REAL(KIND=8) :: RN(3,3)
              REAL(KIND=8) :: RBPN(3,3)
            END SUBROUTINE IAU_PN06
          END INTERFACE 
        END MODULE IAU_PN06__genmod
