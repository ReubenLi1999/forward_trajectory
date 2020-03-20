        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar 12 21:45:27 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE IAU_ATCIQN__genmod
          INTERFACE 
            SUBROUTINE IAU_ATCIQN(RC,DC,PR,PD,PX,RV,ASTROM,N,B,RI,DI)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: RC
              REAL(KIND=8) :: DC
              REAL(KIND=8) :: PR
              REAL(KIND=8) :: PD
              REAL(KIND=8) :: PX
              REAL(KIND=8) :: RV
              REAL(KIND=8) :: ASTROM(30)
              REAL(KIND=8) :: B(8,N)
              REAL(KIND=8) :: RI
              REAL(KIND=8) :: DI
            END SUBROUTINE IAU_ATCIQN
          END INTERFACE 
        END MODULE IAU_ATCIQN__genmod
