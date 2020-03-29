      SUBROUTINE iau_AB ( PNAT, V, S, BM1, PPR )

      IMPLICIT NONE
      double precision PNAT(3), V(3), S, BM1, PPR(3)

      double precision SRS
      PARAMETER ( SRS = 1.97412574336D-08 )

      INTEGER I
      double precision PDV, W1, W2, R2, W, P(3), R


      CALL iau_PDP ( PNAT, V, PDV )
      W1 = 1D0 + PDV/(1D0+BM1)
      W2 = SRS / S
      R2 = 0D0
      DO 1 I=1,3
         W = PNAT(I)*BM1 + W1*V(I) + W2*(V(I)-PDV*PNAT(I))
         P(I) = W
         R2 = R2 + W*W
 1    CONTINUE
      R = SQRT ( R2 )
      DO 2 I=1,3
         PPR(I) = P(I) / R
 2    CONTINUE


      END
