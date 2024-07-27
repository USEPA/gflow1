C     Last change:  HMH  14 Mar 2003   12:06 pm
C     This file contains the following routines and functions
C
C     
C     Logical FUNCTION lw3info        routine returns information for the 2D well
C     LOGICAL function lw3discharge   TRUE if there are discharge specified partially penetrating wells
C     LOGICAL function lw3head        TRUE if there are discharge specified partially penetrating wells
c
c ----------------------------------------------------------------------------
c
       LOGICAL function lw3discharge()
c
c ----------------------------------------------------------------------------
c
C     TRUE if there are discharge specified partially penetrating wells
C
      IMPLICIT NONE
      INTEGER(4) IW
      INCLUDE 'W3COM.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'MATCH.INC'
      lw3discharge=.false.
      IF (NW3.EQ.0) RETURN
      DO 10 IW=1,NW3
      if (lw3q(iw)) lw3discharge=.true.
  10  continue
      return
      end
c
c ----------------------------------------------------------------------------
c
       LOGICAL function lw3head()
c
c ----------------------------------------------------------------------------
c
C     TRUE if there are discharge specified partially penetrating wells
C
      IMPLICIT NONE
      INTEGER(4) IW
      INCLUDE 'W3COM.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'MATCH.INC'
      lw3head=.false.
      IF (NW3.EQ.0) RETURN
      DO 10 IW=1,NW3
      if (.not.lw3q(iw)) lw3head=.true.
  10  continue
      return
      end
c
c -------------------------------------------------------------------------------------------------
c
      LOGICAL function lw3info (cz,czwell,rwell,rqwell,awell)
c
c ----------------------------------------------------------------------------
c
c     Function is TRUE when well at CZ.
c     If true data for well is returned.
C     Routine is called by "lwellinfo" in gftrace.for
c
      implicit none
      INTEGER(4) iw,ien,ist,i
      REAL(8) rwell,rqwell,rdis
      COMPLEX(8) cz,czwell,cz0
      CHARACTER(16) awell
      include 'w3com.inc'
      include 'lusys.inc'
c
      lw3info=.false.
      if (nw3.eq.0) return
      do iw=1,nw3
        IEN=IPNT(IW+1)-1
        IST=IPNT(IW)
      cz0=CMPLX(rw3st(1,ist),rw3st(2,ist))
      rdis=ABS(cz0-cz)
      if (rdis.lt.rw3rad(iw)) then   ! within the radius, close enough
        czwell=cz0
        rwell=rw3rad(iw)
        IEN=IPNT(IW+1)-1
        IST=IPNT(IW)
        rqwell=0.0          ! start of discharge calculations
        DO I=IEN,IST,-1
          rqwell=rqwell+RW3S(I)*RW3L(I)
        END do              ! end of discharge calculations
        awell=aw3lab(iw)
      lw3info=.true.
      end if
      end do
      return
      end

