        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar 12 21:45:11 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE IAU_ATICQN__genmod
          INTERFACE 
            SUBROUTINE IAU_ATICQN(RI,DI,ASTROM,N,B,RC,DC)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: RI
              REAL(KIND=8) :: DI
              REAL(KIND=8) :: ASTROM(30)
              REAL(KIND=8) :: B(8,N)
              REAL(KIND=8) :: RC
              REAL(KIND=8) :: DC
            END SUBROUTINE IAU_ATICQN
          END INTERFACE 
        END MODULE IAU_ATICQN__genmod
