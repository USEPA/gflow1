C     Last change:  HMH  18 Dec 2006    2:17 pm
c --------------------------------------------------------------------------------
c
c     This file contains the following routines or functions
c
c     RFNFDB       returns flow across a line due to all line doublets and line dipoles
c     RFNFDBGAM    coefficient function for flow across line due to added exfiltration rate for inhom.
c     RFNFDBGAMI   coefficient function for flow across line when inside inhomogeneity
c     RFNFDBGAMO   coefficient function for flow across line when outside inhomogeneity
c     DBINTERSECT  calculates intersection of line with line doublet and the flow in the branch cut
c     DFDBAREA     calculates area cut of inhomogeneity by a line
c
c ----------------------------------------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFNFDB(CZ1,CZ2)
c
c ----------------------------------------------------------------------------------------------------------
c
C
C     Returns flow across line between CZ1 and CZ2 due to all line doublets and line dipoles.
C     Includes the contribution of specified exfiltration inside inhomogeneity domains.
C     Includes flow in the branch cut outside the inhomogeneity domain with exfiltration.
C
      IMPLICIT NONE
      INTEGER I
      REAL(8) RFNFDBGAM
      COMPLEX(8) CZ1,CZ2,CDBOM
      INCLUDE 'DBCOM.INC'
      INCLUDE 'TRACOM.INC'
      include 'lusys.inc'
C
      RFNFDB=0.0
      IF (NDB.EQ.0) RETURN
      RFNFDB=RFNFDB+AIMAG(CDBOM(CZ1))-AIMAG(CDBOM(CZ2)) ! calculate flow due to line doublets/dipoles
      DO I=1,NDBSTR
      IF (RDBGAM(I).NE.0.0) THEN ! add flow due to exfiltration rate inside inhomogeneities
      RFNFDB=RFNFDB+RFNFDBGAM(I,CZ1,CZ2)
      END IF
      END DO
      RETURN
      end
C
C ----------------------------------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFNFDBGAM(I,CZ1,CZ2)
c
c ----------------------------------------------------------------------------------------------------------
c
C
C     Returns the flow across a line betweeen CZ1 and CZ2 due to the exfiltration rate (sink density) of 1
c     inside the inhomogeneity domain I  Will only be called when RDBGAM(I) not equal to zero. This condition
C     precludes encountering open strings.
C
      IMPLICIT NONE
      LOGICAL L1,L2,LFDBIN
      INTEGER I
      REAL(8) RFNFDBGAMI,RFNFDBGAMO,RBRANCH
      COMPLEX(8) CZ1,CZ2,CZ3
      INCLUDE 'DBCOM.INC'
      INCLUDE 'TRACOM.INC'
      include 'lusys.inc'
C
      RFNFDBGAM=0.0
      CALL DBPREP (CZ2)
      L2=LFDBIN (I)
      CALL DBPREP (CZ1)
      L1=LFDBIN (I)
      IF (L1.AND.L2) THEN ! both points inside, calculate POISSON contribution
      RFNFDBGAM = RFNFDBGAM + RFNFDBGAMI(I,CZ1,CZ2)
C     No correction for one or more sections of the segment being on the outside!!!
      RETURN
      END IF
      IF (.NOT.L1.AND..NOT.L2) THEN ! both points outside, add flow in branch cut outside
      RFNFDBGAM = RFNFDBGAM + RFNFDBGAMO(I,CZ1,CZ2)
c     No correction for one or more sections of the segment being on the inside!!!
      RETURN
      END IF
      CALL DBINTERSECT(I,CZ1,CZ2,CZ3,RBRANCH) ! intersection CZ3 and flow RBRANCH in branch cut
      RFNFDBGAM = RFNFDBGAM + RBRANCH ! add flow in branch cut along line dipole at CZ3
      IF (L1) THEN ! CZ1 inside and CZ2 outside
      RFNFDBGAM = RFNFDBGAM + RFNFDBGAMI(I,CZ1,CZ3) ! add POISSON contribution
      RFNFDBGAM = RFNFDBGAM + RFNFDBGAMO(I,CZ3,CZ2) ! add flow in branch cut outside
      ELSE         ! CZ2 inside and CZ1 outside
      RFNFDBGAM = RFNFDBGAM + RFNFDBGAMO(I,CZ1,CZ3) ! add flow in branch cut outside
      RFNFDBGAM = RFNFDBGAM + RFNFDBGAMI(I,CZ3,CZ2) ! add POISSON contribution
      END IF
c     No correction for one or more segments being inside or outside!!!
      RETURN
      END
c
c ----------------------------------------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFNFDBGAMI(ISTRING,CZ1,CZ2)
c
c ----------------------------------------------------------------------------------------------------------
c
C
C     Returns the flow across a line between CZ1 and CZ2, that is considered INSIDE domain I
c     Adds flow in branch cut that extends from first node, in case line segment crosses it.
C
      IMPLICIT NONE
      INTEGER ISTRING
      REAL(8) RFAREATRIANGLE,DFDBAREA
      COMPLEX(8) CZ1,CZ2
      INCLUDE 'DBCOM.INC'
      INCLUDE 'TRACOM.INC'
      include 'lusys.inc'
      RFNFDBGAMI=RDBGAM(ISTRING)*
     &(RFAREATRIANGLE(CZ1,CZ2,CDBZ0(ISTRING))-DFDBAREA(ISTRING,CZ1,CZ2))
      RETURN
      END
c
c ----------------------------------------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFNFDBGAMO(ISTRING,CZ1,CZ2)
c
c ----------------------------------------------------------------------------------------------------------
c
C
C     Returns the flow across a line between CZ1 and CZ2, that is considered OUTSIDE domain I
C
      IMPLICIT NONE
      INTEGER ISTRING,INOD1
      REAL(8) RFBRANCH,DFDBAREA
      COMPLEX(8) CZ1,CZ2
      INCLUDE 'DBCOM.INC'
      INCLUDE 'TRACOM.INC'
      include 'lusys.inc'
C     Flow due to line doublets/dipoles already incorporated
c     in Delta-PSI calculated in function RFNFDB.
C     Add flow in branch cut if it is being crossed.
      INOD1=IDBSTA(ISTRING) ! start of branch cut outside the domain
      RFNFDBGAMO=
     & RDBGAM(ISTRING)*(DBAREA(ISTRING)*RFBRANCH(CDBZ(INOD1),CZ1,CZ2)-
     & DFDBAREA(ISTRING,CZ1,CZ2))
      RETURN
      END
c
c ----------------------------------------------------------------------------------------------------------
c
      SUBROUTINE DBINTERSECT (ISTRING,CZ1,CZ2,CZ3,RBRANCH)
c
c ----------------------------------------------------------------------------------------------------------
c
c
c     Routine calculates the intersection CZ3 and the flow in the branch cut along the dipole at CZ3
C     Input:   I       domain number
c              CZ1     starting point of line segment
c              CZ2     end point of line segment
c     Output:  CZ3     intersection of line segment and line dipole in string
c              RBRANCH flow in branch cut of line dipole at CZ3
c
      IMPLICIT NONE
      INTEGER INOD,INOD1,INODL,INODM1,ISTRING,I0,I1
      LOGICAL LBRANCH
      REAL(8) RBRANCH,RDIS,RDISMIN,RS,RSIGN,DS, D_ONE
      COMPLEX(8) CZ1,CZ2,CZ3,CZC1,CZC2
      INCLUDE 'DBCOM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      D_ONE = 1.0
      RBRANCH=0.0
      CZC1=0.5*(CZ1+CZ2) ! center of line segment
      RDISMIN=1.0e21
      INOD1=IDBSTA(ISTRING)
      INODL=INOD1+NDBSTI(ISTRING)-1
      INODM1=INODL
      DO INOD=INOD1,INODL ! first line dipole in loop is last line dipole in string
      CZC2=0.5*(CDBZ(INODM1)+CDBZ(INOD)) ! center of line dipole
      RDIS=ABS(CZC1-CZC2)
      IF (RDIS.LT.RDISMIN) THEN  ! find closest line dipole
      I0=INODM1
      I1=INOD
      RDISMIN=RDIS
      END IF
      INODM1=INOD
      END DO
      CALL BRANCHCUT ! try if intersection on closest line dipole
     &     (CZ1,CZ2,CDBZ(I0),CDBZ(I1),CZ3,RBRANCH,LBRANCH)
      IF (LBRANCH) THEN
        IF (I1.EQ.INOD1) THEN ! on last line dipole in string (is first in loop)
        DS=-DBAREA(ISTRING)*RDBGAM(ISTRING) ! set dipole strength at total discharge
        ELSE
        DS=AIMAG(CDDBSS(I1))
        END IF
        RS=(DS-AIMAG(CDDBSS(I0)))
        RSIGN=SIGN(D_ONE,RBRANCH)  ! preserve sign
        RBRANCH=ABS(RBRANCH)/ABS(CDBZ(I1)-CDBZ(I0)) ! make rbranch vary along CZS-CZE between 0 and 1
        RBRANCH=(1.0-RBRANCH)*RS+(AIMAG(CDDBSS(I0))) ! calculate jump in PSI across line dipole at intersection
        RBRANCH=RSIGN*RBRANCH ! apply proper sign
        RETURN
      ENDIF
!     Intersection appears not to be on closest line dipole, look further
      INODM1=INODL
      DO INOD=INOD1,INODL
      IF (INOD.NE.I1) THEN ! don't try closest line dipole again
        CALL BRANCHCUT ! find intersection on one of the other line dipoles
     &       (CZ1,CZ2,CDBZ(INODM1),CDBZ(INOD),CZ3,RBRANCH,LBRANCH)
        IF (LBRANCH) THEN
        IF (INOD.EQ.INOD1) THEN ! on last line dipole in string (is first in loop)
        DS=-DBAREA(ISTRING)*RDBGAM(ISTRING) ! set dipole strength at total discharge
        ELSE
        DS=AIMAG(CDDBSS(INOD))
        END IF
          RS=(DS-AIMAG(CDDBSS(INODM1)))
          RSIGN=SIGN(D_ONE,RBRANCH) !preserve sign
          RBRANCH=ABS(RBRANCH)/ABS(CDBZ(INOD)-CDBZ(INODM1)) ! make rbranch vary along CZS-CZE between 0 and 1
          RBRANCH=(1.0-RBRANCH)*RS+(AIMAG(CDDBSS(INODM1))) !  calculate jump in PSI across line dipole at intersection
          RBRANCH=RSIGN*RBRANCH ! apply proper sign
          INODM1=INOD
          RETURN
        ENDIF
      END IF
      INODM1=INOD
      END DO
      END SUBROUTINE
c
c ----------------------------------------------------------------------------------------------------------
c
      REAL*8 FUNCTION DFDBAREA(ISTRING,CZ1,CZ2)
c
c ----------------------------------------------------------------------------------------------------------
c
C
C     Function returns the area that the line CZ1-CZ2 cuts from an inhomogeneity domain if CZ1 and CZ2
C     are both outside the domain. If CZ1 and CZ2 are both inside the domain the area of subdomains cut
C     of by the line CZ1 and CZ2 which do not belong to the domain is returned.
C
      IMPLICIT NONE
      INTEGER ISTRING,INOD,INOD1,INODM1,INODL
      LOGICAL LAREA, LINSECTLINE
      REAL(8) DA,DFA,DS, D_ONE
      COMPLEX(8) CZ1,CZ2,CZZ1,CZZ2,CZ0
      INCLUDE 'DBCOM.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
C
      DFDBAREA=0.0
      LAREA=.FALSE.
      DS=1.0 ! test eliminating DS, FIX LATER !!!!!!
      INOD1=IDBSTA(ISTRING)
      INODL=INOD1+NDBSTI(ISTRING)-1
      INODM1=INODL
      DO INOD=INOD1,INODL
      IF (LINSECTLINE(CDBZ(INODM1),CDBZ(INOD),CZ1,CZ2,CZ0)) THEN
        IF (LAREA) THEN ! already accumulating area, finish.
           CZZ1=CZZ2
           CZZ2=CZ0-CZ1
           DA=DA+DFA(CZZ1,CZZ2)
           DFDBAREA=DFDBAREA+DS*DA
           LAREA=.FALSE.
        ELSE            ! start accumulating new area
           CZZ1=CZ0-CZ1
           CZZ2=CDBZ(INOD)-CZ1
           DA=DFA(CZZ1,CZZ2)
           D_ONE=1.0
           DS=SIGN(D_ONE,DA)
           LAREA=.TRUE.
        END IF
      ELSE
        IF (LAREA) THEN ! no intersection, but continue to accumulate area
          CZZ1=CZZ2
          CZZ2=CDBZ(INOD)-CZ1
          DA=DA+DFA(CZZ1,CZZ2)
        END IF
      END IF
      IF ((INOD.EQ.INOD1).AND.LAREA) THEN
      DA=DA-DBAREA(ISTRING) ! correct for branch cut
      END IF
      INODM1=INOD
      END DO
      RETURN
      END

