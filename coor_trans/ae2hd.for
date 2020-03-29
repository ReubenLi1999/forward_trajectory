      SUBROUTINE iau_AE2HD (AZ, EL, PHI, HA, DEC)


      IMPLICIT NONE

      double precision AZ, EL, PHI, HA, DEC

      double precision SA, CA, SE, CE, SP, CP, X, Y, Z, R

      SA = SIN(AZ)
      CA = COS(AZ)
      SE = SIN(EL)
      CE = COS(EL)
      SP = SIN(PHI)
      CP = COS(PHI)

      X = - CA*CE*SP + SE*CP
      Y = - SA*CE
      Z = CA*CE*CP + SE*SP

      R = SQRT(X*X + Y*Y)
      IF ( R.EQ.0D0 ) THEN
         HA = 0D0
      ELSE
         HA = ATAN2(Y,X)
      END IF
      DEC = ATAN2(Z,R)



      END
