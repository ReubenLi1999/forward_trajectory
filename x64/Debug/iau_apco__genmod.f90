        !COMPILER-GENERATED INTERFACE MODULE: Sun Mar 29 09:40:35 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE IAU_APCO__genmod
          INTERFACE 
            SUBROUTINE IAU_APCO(DATE1,DATE2,EBPV,EHP,X,Y,S,THETA,ELONG, &
     &PHI,HM,XP,YP,SP,REFA,REFB,ASTROM)
              REAL(KIND=8) :: DATE1
              REAL(KIND=8) :: DATE2
              REAL(KIND=8) :: EBPV(3,2)
              REAL(KIND=8) :: EHP(3)
              REAL(KIND=8) :: X
              REAL(KIND=8) :: Y
              REAL(KIND=8) :: S
              REAL(KIND=8) :: THETA
              REAL(KIND=8) :: ELONG
              REAL(KIND=8) :: PHI
              REAL(KIND=8) :: HM
              REAL(KIND=8) :: XP
              REAL(KIND=8) :: YP
              REAL(KIND=8) :: SP
              REAL(KIND=8) :: REFA
              REAL(KIND=8) :: REFB
              REAL(KIND=8) :: ASTROM(30)
            END SUBROUTINE IAU_APCO
          END INTERFACE 
        END MODULE IAU_APCO__genmod
