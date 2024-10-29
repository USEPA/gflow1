       CHARACTER(1) BLCHRT,BLRES,BLSPD,BLINE
c  This is a fix by VAK -- F90 type-checking requires that
c  arguments to subroutines match.  The "implicit equivalencing"
c  of the I*2 array over the C*1 array is baaaad news in F90
       INTEGER(2) IBLCHRT(2)
       EQUIVALENCE (BLCHART,IBLCHRT(1))
       INTEGER(4) IBLU
       REAL(8) XA,YA,XB,YB
       LOGICAL LBON
       COMMON/BLPLOT1/ BLCHRT,BLRES,BLSPD,BLINE
       COMMON/BLPLOT2/ IBLU,LBON
       COMMON/MARGC/ XA,YA,XB,YB
     