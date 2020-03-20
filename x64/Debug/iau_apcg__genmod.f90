        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar 12 21:45:38 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE IAU_APCG__genmod
          INTERFACE 
            SUBROUTINE IAU_APCG(DATE1,DATE2,EBPV,EHP,ASTROM)
              REAL(KIND=8) :: DATE1
              REAL(KIND=8) :: DATE2
              REAL(KIND=8) :: EBPV(3,2)
              REAL(KIND=8) :: EHP(3)
              REAL(KIND=8) :: ASTROM(30)
            END SUBROUTINE IAU_APCG
          END INTERFACE 
        END MODULE IAU_APCG__genmod
