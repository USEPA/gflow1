C     Last change:  HMH   4 Apr 2007    8:35 pm
c     This file contains the following routines and functions:
c
c     CFUFOM             calculates the complex potential due to uniform flow
c     GVQI               calculates the discharge vector due to uniform flow
c     RFCONPT            adds integration constant to the real discharge potential
c     SetInitialConstant called in SOLUT when nsol=0 to initialize integration constant for discharge potential
C
c -------------------------------------------------------------------------------------------------------
c
      COMPLEX(8) FUNCTION CFUFOM(CZ)
C
c -------------------------------------------------------------------------------------------------------
c
C
C     Calculates contribution to potential due to uniform flow.
C
      IMPLICIT NONE
      COMPLEX(8) CZ
      INCLUDE 'GVCOM.INC'
      INCLUDE 'TRACOM.INC'
      CFUFOM=(0.0,0.0)
      IF (RQ0.EQ.0.0) RETURN
      CFUFOM=-RQ0*CZ*EXP(-CUNALP) 
      RETURN
      END
c          
C-------------------------------------------------------------------------------------------------------
C     
	SUBROUTINE GVQI (CZ,RVI)
C
C-------------------------------------------------------------------------------------------------------
C    
C     Discharge due to uniform flow
C
      IMPLICIT NONE
      REAL(8) RVI
      COMPLEX(8) CZ,CDUM
      INCLUDE 'GVCOM.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'      
      DIMENSION RVI(3)
      IF (ABS(CZ).GT.1.0E30) THEN
      WRITE (ILUER,1000) CZ
 1000 FORMAT (' ***WARNING in GVQI, CZ out of range: ',2G11.4)
      ENDIF      
      IF (RQ0.EQ.0) RETURN
      CDUM=RQ0*EXP(-CUNALP)
      RVI(1)=RVI(1)+REAL(CDUM)
      RVI(2)=RVI(2)-AIMAG(CDUM)
      RETURN
      END
C
C -------------------------------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFCONPT(CZ)
C
c -------------------------------------------------------------------------------------------------------
c
C     Calculates constant potential.
C
      IMPLICIT NONE
      COMPLEX(8) CZ
      INCLUDE 'GVCOM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'      
C
      IF (ABS(CZ).GT.1.0E30) THEN
      WRITE (ILUER,1000) CZ
 1000 FORMAT (' ***WARNING in RFCONPT, CZ out of range: ',2G11.4)
      ENDIF      
      RFCONPT=RPOTC
      RETURN
      END
c
c ---------------------------------------------------------------------------------------------------------------
c
      subroutine SetInitialConstant ()
c
c ---------------------------------------------------------------------------------------------------------------
c
c     Routine is called in SOLUT when nsol=0 to set the initial integration constant equal to
c     the potential calculated at the reference point minus the initial result from the potential function.
c     This strategy forces the initial value of the calculated potential at the reference point (calculated
c     with all the unknown strength parameters equal to zero) to equal the specified potential at the
c     reference point.
c     This is arbitrary, but results in a better first GW solution.
c
      implicit none
      REAL(8) rfpoth,rfpot
      include 'gvcom.inc'
c
      rpotc=rfpoth (rhead0,crefz)-rfpot(crefz)
c
      return
      end subroutine


