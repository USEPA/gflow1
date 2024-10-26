C     Last change:  HMH   4 Apr 2007    3:17 pm
c       This file contains the following subroutines and functions:
c
c       gvmatsize updates the matrix array dimensions
c       GVCZI     generates reference control point
c       GVMAT     generates matrix coefficients for reference point
c       GVMATCORRECT provides an update for the matrix column for the reference point
c       GVKNO     generates known vector for reference point
C       GVSUB     substitutes solution vector (element) into integration constant
c       GVUPDATE  routine is a dummy routine for the update of the knownvector
c
c
c ------------------------------------------------------------------------------------------------------
c
      subroutine gvmatsize (M,N)
c
c ------------------------------------------------------------------------------------------------------
c
c     Routine updates the matrix array sizes M and N.
c
      implicit none
      INTEGER M,N
c
      M=M+1
      N=N+1
      return
c
      end subroutine
c
c -------------------------------------------------------------------------------------
c
      subroutine gvkeep(drb,j,m)
c
c --------------------------------------------------------------------------------------
c
c
c     Add an element to drb array and set to zero to represent the constant.
c     The drb array will contain any strength differences between the
c     solution and the actual stored strength parameters, see lskeepsigma and lkkeepstrength routines.
c
      INTEGER j,m
      REAL(8) drb
      DIMENSION drb(m)
         j=j+1
         drb(j)=0.0d0
      return
      end subroutine
C
C-------------------------------------------------------------------------------------------------------
C
	SUBROUTINE GVCZI (CZC,M,N,DRFAC,ITYPE)
C
C-------------------------------------------------------------------------------------------------------
C
C     REFERENCE POINT EQUATION
      IMPLICIT NONE
      INTEGER(4) M,N,ITYPE
      REAL(8) DRFAC
      COMPLEX(8) CZC
      INCLUDE 'gvcom.inc'
      INCLUDE 'tracom.inc'
      DIMENSION CZC(*),DRFAC(4,*),ITYPE(*)
      M=M+1  ! add one equation for the reference point
      N=N+1  ! add one unknown (constant potential) for the reference point
      CZC(M)=CREFZ
      ITYPE(M)=1
      DRFAC(1,M)=1.0D00
      DRFAC(4,M)=1.0D0
      RETURN
      END
C
C-------------------------------------------------------------------------------------------------------
C
      SUBROUTINE gvgenmat (DRA,M,N,J,DRFAC,ITYPE)
C
C-------------------------------------------------------------------------------------------------------
C
C     REFERENCE POINT EQUATION AND CONTINUITY EQUATION
      IMPLICIT NONE
      INTEGER(4) M,N,J,ITYPE,I
      LOGICAL lset1
      REAL(8) DRA,DRFAC
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      DIMENSION DRA(M,*),DRFAC(4,*),ITYPE(*)
C
C ITYPE=+1 potential specified at CZ (Note: at doublet nodes ITYPE=+1,
C          but potential is expressed in terms of local strenght parameters.)
C       -1 difference in potential specified: PHI(CZ)-PHI(CZA)
C       +2 stream function specified at CZ
C       -2 flow across line between CZ and CZA, positive to the left when at CZ
C          note: because PSI is continuous for line doublets, this is the same as
C          the difference bwteen CZ and CZA: PSI(CZ)-PSI(CZA)
C       +3 discharge component normal to a line 
C       +4 discharge component parallel to a line
C        5 continuity equation: provide total discharge
C        6 request for zero matrix coefficient
C
c
c      Note: for the initial matrix itype=+/-1 and itype=6 are interpreted as potential specified
c
      J=J+1
!      WRITE (ILUME,1000) N,J
      DO 10 I=1,M
      lset1=ABS(itype(i)).eq.1.or.itype(i).eq.6
      DRA(I,J)=0.0   ! for all type of equations, except for potential specified conditions
      IF (lset1) DRA(I,J)=1.0d0 ! provide 1.0 for potential specified equations
  10  CONTINUE
      RETURN
 1000 FORMAT ('+Generating',I4,' equations, doing equation #: ',I4)      
      END
C-------------------------------------------------------------------------------------------------------
C
      SUBROUTINE gvmatcorrect (DRA,M,N,J,DRFAC,ITYPE)
C
C-------------------------------------------------------------------------------------------------------
C
C     REFERENCE POINT EQUATION AND CONTINUITY EQUATION
      IMPLICIT NONE
      INTEGER(4) M,N,J,ITYPE,I
      REAL(8) DRA,DRFAC
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      DIMENSION DRA(M,*),DRFAC(4,*),ITYPE(*)
C
C ITYPE=+1 potential specified at CZ (Note: at doublet nodes ITYPE=+1,
C          but potential is expressed in terms of local strenght parameters.)
C       -1 difference in potential specified: PHI(CZ)-PHI(CZA)
C       +2 stream function specified at CZ
C       -2 flow across line between CZ and CZA, positive to the left when at CZ
C          note: because PSI is continuous for line doublets, this is the same as
C          the difference bwteen CZ and CZA: PSI(CZ)-PSI(CZA)
C       +3 discharge component normal to a line 
C       +4 discharge component parallel to a line
C        5 continuity equation: provide total discharge
C        6 request for zero matrix coefficient
C
      J=J+1
!      WRITE (ILUME,1000) N,J
      DO 10 I=1,M
      DRA(I,J)=0.0   ! for all type of equations, except potential specified conditions
      IF (ITYPE(I).EQ.1) DRA(I,J)=1.0d0 ! provide 1.0 for potential specified equations
  10  CONTINUE
      RETURN
 1000 FORMAT ('+Correcting',I4,' equations, doing equation #: ',I4)
      END
C

C
C-------------------------------------------------------------------------------------------------------
C
      SUBROUTINE GVKNO (DRB,J)
C
C-------------------------------------------------------------------------------------------------------
C
C     REFERENCE POINT EQUATION
      IMPLICIT NONE
      INTEGER(4) J
      REAL(8) DRB,RFTOP,RFPOT,RFPOTH
      COMPLEX(8) CZ
      INCLUDE 'gvcom.inc'
      INCLUDE 'com3d.inc'
      INCLUDE 'tracom.inc'
      DIMENSION DRB(*)
      J=J+1
      CZ=CREFZ
      R3DZ=RFTOP(CZ)
      DRB(J)=DRB(J)+RFPOTH(RHEAD0,CZ)-RFPOT(CZ)
      RETURN
      END
C
C-------------------------------------------------------------------------------------------------------
C
      SUBROUTINE GVSUB (DRB,J)
C
C-------------------------------------------------------------------------------------------------------
C
C     REFERENCE POINT EQUATION AND CONTINUITY EQUATION
      IMPLICIT NONE
      INTEGER(4) J
      REAL(8) DRB
      INCLUDE 'gvcom.inc'
      INCLUDE 'tracom.inc'
      DIMENSION DRB(*)
      J=J+1
      RPOTC=RPOTC+DRB(J)
      RETURN
      END
c
c ------------------------------------------------------------------------------------------------------
c
      subroutine gvupdate(drscr,j,m,n)
c
c ------------------------------------------------------------------------------------------------------
c
      implicit none
      INTEGER j,m,n
      REAL(8) drscr
      DIMENSION drscr(n)
c
c     no update required
c
      j=j+1
c
      end subroutine
c
c ------------------------------------------------------------------------------------------------------
c
      subroutine gvupdate_check (j,m,n)
c
c ------------------------------------------------------------------------------------------------------
c
      implicit none
      INTEGER j,m,n
      REAL(8) drscr
      DIMENSION drscr(n)
c
c     no update required, nothing to check
c
      j=j+1
c
      end subroutine
