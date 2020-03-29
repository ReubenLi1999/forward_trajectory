        !COMPILER-GENERATED INTERFACE MODULE: Sun Mar 29 09:40:31 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE IAU_APCI__genmod
          INTERFACE 
            SUBROUTINE IAU_APCI(DATE1,DATE2,EBPV,EHP,X,Y,S,ASTROM)
              REAL(KIND=8) :: DATE1
              REAL(KIND=8) :: DATE2
              REAL(KIND=8) :: EBPV(3,2)
              REAL(KIND=8) :: EHP(3)
              REAL(KIND=8) :: X
              REAL(KIND=8) :: Y
              REAL(KIND=8) :: S
              REAL(KIND=8) :: ASTROM(30)
            END SUBROUTINE IAU_APCI
          END INTERFACE 
        END MODULE IAU_APCI__genmod
