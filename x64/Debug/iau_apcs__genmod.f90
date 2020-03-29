        !COMPILER-GENERATED INTERFACE MODULE: Sun Mar 29 09:40:24 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE IAU_APCS__genmod
          INTERFACE 
            SUBROUTINE IAU_APCS(DATE1,DATE2,PV,EBPV,EHP,ASTROM)
              REAL(KIND=8) :: DATE1
              REAL(KIND=8) :: DATE2
              REAL(KIND=8) :: PV(3,2)
              REAL(KIND=8) :: EBPV(3,2)
              REAL(KIND=8) :: EHP(3)
              REAL(KIND=8) :: ASTROM(30)
            END SUBROUTINE IAU_APCS
          END INTERFACE 
        END MODULE IAU_APCS__genmod
