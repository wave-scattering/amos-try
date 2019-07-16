      SUBROUTINE ZDIV(AR, AI, BR, BI, CR, CI)
C***BEGIN PROLOGUE  ZDIV
C***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
C
C     DOUBLE PRECISION COMPLEX DIVIDE C=A/B.
C
C***ROUTINES CALLED  AZABS
C***END PROLOGUE  ZDIV
      real(8) AR, AI, BR, BI, CR, CI, BM, CA, CB, CC, CD
      real(8) AZABS
      BM = 1.0D0/AZABS(BR,BI)
      CC = BR*BM
      CD = BI*BM
      CA = (AR*CC+AI*CD)*BM
      CB = (AI*CC-AR*CD)*BM
      CR = CA
      CI = CB
      RETURN
      END
