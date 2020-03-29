      SUBROUTINE iau_A2AF ( NDP, ANGLE, SIGN, IDMSF )


      IMPLICIT NONE

      INTEGER NDP
      double precision ANGLE
      CHARACTER SIGN*(*)
      INTEGER IDMSF(4)

      double precision D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

      double precision F
      PARAMETER ( F = 15D0/D2PI )

      CALL iau_D2TF ( NDP, ANGLE*F, SIGN, IDMSF )



      END
