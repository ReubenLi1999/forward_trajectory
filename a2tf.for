      SUBROUTINE iau_A2TF ( NDP, ANGLE, SIGN, IHMSF )


      IMPLICIT NONE

      INTEGER NDP
      double precision ANGLE
      CHARACTER SIGN*(*)
      INTEGER IHMSF(4)

      double precision D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

      CALL iau_D2TF ( NDP, ANGLE/D2PI, SIGN, IHMSF )



      END
