C     Last change:  HMH  18 Dec 2006    2:17 pm
c --------------------------------------------------------------------------------
c
c     This file contains the following routines or functions
c
c     CDBOM  complex potential due to all line doublets
c     DBQI   discharge vector due to all line doublets
c     CFDBF  coefficient function F for discharge potential
c     CFDBG  coefficient function G for discharge potential
c     CFDBS  coefficient function S for discharge potential
c     CFDBFD coefficient function F for derivative of discharge potential
c     CFDBGD coefficient function G for derivative of discharge potential
c     CFDBSD coefficient function S for derivative of discharge potential
c     RFDBPT driver to add all contributions due to added exfiltration to discharge potential
c     RFDBRP returns contribution to the discharge potential due to added exfiltration for one inhomogeneity
c
c -----------------------------------------------------------------------------------------------------------
c
      COMPLEX(8) FUNCTION CDBOM (CZ)
c
c -----------------------------------------------------------------------------------------------------------
C
C     Routine returns the complex potential due to all line doublet strings.
C
      IMPLICIT NONE
      INTEGER(4) ISTR,INOD1,INODL,INOD,INODM1,INODP1
      REAL(8) RQ,DQS,roffset
      COMPLEX(8) CZ,CFDBF,CFDBS,CFDBG,CDDS,CDDBOM,CZ0,CBZ,
     &           cdds1,cdds2
      INCLUDE 'DBCOM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      CDBOM=(0.0,0.0)
      IF (NDB.EQ.0) RETURN
      CALL DBPREP (CZ)
      CZ0=CZ                   ! start logic to offset point (2/7/05 modification)
c      write (iluer,1004) ldbmatrix,dbeps
c 1004 format ('DBOM4: ldbmatrix= ',l2,' dbeps=',d14.7)
      if (.not.ldbmatrix) then
      DO ISTR=1,NDBSTR
      INOD1=IDBSTA(ISTR)
      INODL=INOD1+NDBSTI(ISTR)-1
      DO INOD=INOD1,INODL
      IF (INOD.EQ.INODL.AND.IDOMAINTYPE(ISTR).NE.2) THEN
      INODP1=INOD1
      ELSE
      INODP1=INOD+1
      END IF
c      IF (LDBNOD(inod).OR.LDBNOD(INODP1)) THEN   ! offset to avoid singular potential
      IF (LDBNOD(inod)) THEN   ! offset to avoid singular potential (use only first node)
c        write (iluer,1002)
c 1002 format (' Invoked offset logic in CDBOM.')
        roffset=2.1d+0*dbeps/ABS(cdbz21(inod)) ! roffset is 5% over the value of DBEPS
        CBZ=CMPLX(-1.0d+0,roffset) ! move inside domain
        CZ=0.5*(CDBZ21(INOD)*CBZ+CDBZ(INOD)+CDBZ(INODP1))  ! map back to physical plane
c        write (iluer,1003) cz0,cz
c 1003 format ('old point=',d14.7,d14.7,' new point=',d14.7,d14.7)
        CALL DBPREP (CZ)  ! regenerate line doublet data
        GOTO 7  ! once offset, stop looking for point to coincide with node
      ENDIF
      END DO
      END DO
  7   CONTINUE                 ! end logic to offset point
      endif
      CDDBOM=(0.0D0,0.0D0)
      DO 20 ISTR=1,NDBSTR
      INOD1=IDBSTA(ISTR)
      INODL=INOD1+NDBSTI(ISTR)-1
      do inod=inod1,inodl
       IF (RDBTFACNOD(INOD1).EQ.1.0E+21.AND.RDBGAM(ISTR).EQ.0.0) GOTO 20 ! skip
       inodp1=inod+1
       if (inod.eq.inodl) inodp1=inod1
       cdds1=cddbss(inod)+dbstr(inod)
       cdds2=cddbss(inodp1)+dbsen(inod) ! note: cddbss only defined at starting nodes
       cddbom=cddbom+cdds1*cfdbf(inod)+cdds2*cfdbg(inod)
       DQS=DBQSTS(INOD)+DBQSTR(INOD)
       CDDBOM=CDDBOM+DQS*CFDBS(INOD)
      end do
C      IF (IDOMAINTYPE(ISTR).EQ.2) THEN ! open string
C       CDDBOM=CDDBOM+DBQSTR(INOD1)*CFDBS(INOD1) ! note: F function not needed for first line doublet (strength zero)
C       DO INOD=INOD1+1,INODL
C       CDDBOM=CDDBOM+DBSTR(INOD)*(CFDBG(INOD-1)+CFDBF(INOD))
C       CDDBOM=CDDBOM+DBQSTR(INOD)*CFDBS(INOD)
C       END DO ! note: G function not needed for last line doublet (strength zero)
C      ELSE                             ! closed string
C       INODM1=INODL
C       IF (RDBTFACNOD(INOD1).EQ.1.0E+21.AND.RDBGAM(ISTR).EQ.0.0) GOTO 20 ! skip
C       DO 10 INOD=INOD1,INODL
C       CDDS=CDDBSS(INOD)+DBSTR(INOD)
C       CDDBOM=CDDBOM+CDDS*(CFDBG(INODM1)+CFDBF(INOD))
C       DQS=DBQSTS(INOD)+DBQSTR(INOD)
C       CDDBOM=CDDBOM+DQS*CFDBS(INOD)
C       INODM1=INOD
C  10   CONTINUE
      IF (RDBGAM(ISTR).NE.0.0) THEN
        RQ=RDBGAM(ISTR)*DBAREA(ISTR)
        IF (LDBNOD(INOD1)) THEN      ! remove singularity caused by dipole string   (should not be true after 2/7/05)
         IF (.not.ldbmatrix) write (iluer,1001)
 1001 format (' Error: found point to coincide with node in CDBOM.')
          CDDBOM=CDDBOM+RQ*(LOG(CDBZ(INOD1)-CDBZ(INODL))-1.0)/D2PI
        ELSE
          CDDBOM=CDDBOM-CDI*RQ*CFDBG(INODL)
          CDDBOM=CDDBOM+RQ*LOG(CZ-CDBZ(INOD1))/D2PI ! branch cut from first node (for added exfiltration)
        ENDIF
      ENDIF
C      ENDIF
  20  CONTINUE
      CZ=CZ0          ! restore original point in case of an offset
      CDBOM=CDDBOM
      RETURN
      END
c
C --------------------------------------------------------------------------------------------------------
c
      SUBROUTINE DBQI (CZ,RQI)
c
c --------------------------------------------------------------------------------------------------------
C
C     Routine generates the discharge vector components RQI(1) and RQI(2).
C     Vertical flow due to recharge is being added.
C
      IMPLICIT NONE
      LOGICAL lfinterface
      INTEGER(4) ISTR,INOD1,INODL,INOD,INODP1,INODM1
      REAL(8) R3DUM,RBLOC,RHGHT,RQI,RFBASE,RFHGHT,DQS,RQ,
     &        rfinterface
      COMPLEX(8) CZ,CDDUM,CFDBFD,CFDBSD,CFDBGD,CDDS,CBZ,CZ0,
     &           cdds1,cdds2
      INCLUDE 'DBCOM.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
      DIMENSION RQI(3)
C
      IF (NDB.EQ.0.0) RETURN
      CALL DBPREP (CZ)
      CZ0=CZ
      DO ISTR=1,NDBSTR
      INOD1=IDBSTA(ISTR)
      INODL=INOD1+NDBSTI(ISTR)-1
      DO INOD=INOD1,INODL
      IF (INOD.EQ.INODL.AND.IDOMAINTYPE(ISTR).NE.2) THEN
      INODP1=INOD1
      ELSE
      INODP1=INOD+1
      END IF
      IF (LDBNOD(INOD)) THEN   ! offset to avoid singular discharge
        cz=cdbz(inod)+CMPLX(1.1d+0*dbeps,0.0d+0) ! move point horizontally
        CALL DBPREP (CZ)  ! regenerate line doublet data
        GOTO 7  ! once offset, stop looking for point to coincide with node
      ENDIF
      END DO
      END DO
  7   CONTINUE
      CDDUM=(0.0D0,0.0D0)
      R3DUM=0.0
      RBLOC=RFBASE(CZ)
      RHGHT=RFHGHT(CZ)
      if (lfinterface()) then
      rbloc=MAX(rbloc,rfinterface(cz))
      end if
      IF (RHGHT.LT.0.001) THEN
        WRITE (ILUER,1000)
        RETURN
      ENDIF
      DO 20 ISTR=1,NDBSTR
      INOD1=IDBSTA(ISTR)
      INODL=INOD1+NDBSTI(ISTR)-1
      IF (RDBTFACNOD(INOD1).EQ.1.0E+21.AND.RDBGAM(ISTR).EQ.0.0) GOTO 20 ! skip
      do inod=inod1,inodl
       inodp1=inod+1
       if (inod.eq.inodl) inodp1=inod1
       cdds1=cddbss(inod)+dbstr(inod)
       cdds2=cddbss(inodp1)+dbsen(inod) ! note: cddbss only defined at starting nodes
       cddum=cddum+cdds1*cfdbfd(inod)+cdds2*cfdbgd(inod)
       DQS=DBQSTS(INOD)+DBQSTR(INOD)
       CDDUM=CDDUM+DQS*CFDBSD(INOD)
      end do
c      IF (IDOMAINTYPE(ISTR).EQ.2) THEN ! open string
c      CDDUM=CDDUM+DBQSTR(INOD1)*CFDBSD(INOD1) ! note: F function not needed for first line doublet (strength zero)
c      DO INOD=INOD1+1,INODL
c      CDDUM=CDDUM+DBSTR(INOD)*(CFDBGD(INOD-1)+CFDBFD(INOD))
c      CDDUM=CDDUM+DBQSTR(INOD)*CFDBSD(INOD)
c      END DO ! note: G function not needed for last line doublet (strength zero)
c      ELSE                             ! closed string
c      INODM1=INODL
c      IF (RDBTFACNOD(INOD1).EQ.1.0E+21.AND.RDBGAM(ISTR).EQ.0.0) GOTO 20 ! skip
c      DO 10 INOD=INOD1,INODL
c      CDDS=CDDBSS(INOD)+DBSTR(INOD)
c      CDDUM=CDDUM+CDDS*(CFDBGD(INODM1)+CFDBFD(INOD))
c      DQS=DBQSTS(INOD)+DBQSTR(INOD)
c      CDDUM=CDDUM+DQS*CFDBSD(INOD)
c      INODM1=INOD
c  10  CONTINUE
      IF (RDBGAM(ISTR).NE.0.0) THEN
      RQ=RDBGAM(ISTR)*DBAREA(ISTR)
c conditional statement disabled in view of offset
c      IF (.NOT.LDBNOD(INOD1)) THEN ! remove singularity caused by dipole string
      CDDUM=CDDUM-CDI*RQ*CFDBGD(INODL)
      CDDUM=CDDUM-RQ/(CZ-CDBZ(INOD1))/D2PI
c      ENDIF
      IF (LDBINS(ISTR)) THEN
        CDDUM=CDDUM-0.5*RDBGAM(ISTR)*CONJG(CZ-CDBZ0(ISTR))
        R3DUM=R3DUM+(R3DZ-RBLOC)*RDBGAM(ISTR)/RHGHT ! ----- add vertical flow
      ENDIF
      ENDIF
c      ENDIF
  20  CONTINUE
      CZ=CZ0  ! undo possible offset applied to avoid singular discharge
      RQI(1)=RQI(1)+CDDUM
      RQI(2)=RQI(2)-AIMAG(CDDUM)
      RQI(3)=RQI(3)+R3DUM
      RETURN
 1000 FORMAT (' ***ERROR in INHOMOGENEITY module: zero aquifer height.')      
      END
c
c ------------------------------------------------------------------------------------------------------------
c
      COMPLEX(8) FUNCTION CFDBF (IEL)
c
c ------------------------------------------------------------------------------------------------------------
c
C
C     POTENTIAL due to the F function for line doublet # IEL
C     Prior to this function the routine DBPREP(CZ) should have
C     been called to assign values to variables in common.
C
      IMPLICIT NONE
      INTEGER(4) IEL,I
      REAL(8) DCOF
      COMPLEX(8) CDDUM
      INCLUDE 'DBCOM.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
      DIMENSION DCOF(21)
      DATA DCOF /-1.0D0,
     &           0.3333333333333333D0,-0.3333333333333333D0,
     &           0.2000000000000000D0,-0.2000000000000000D0,
     &           0.1428571428571428D0,-0.1428571428571428D0,
     &           0.1111111111111111D0,-0.1111111111111111D0,
     &           0.0909090909090909D0,-0.0909090909090909D0,
     &           0.0769230769230769D0,-0.0769230769230769D0,
     &           0.0666666666666666D0,-0.0666666666666666D0,
     &           0.0588235294117647D0,-0.0588235294117647D0,
     &           0.0526315789473684D0,-0.0526315789473684D0,
     &           0.0476190476190476D0,-0.0476190476190476D0/
C
      IF (.NOT.LDBFAR(IEL)) THEN
            CDDUM=CMPLX(CDBZEL(IEL))-1.0D0
            CFDBF=(1.0D0+0.5D0*CDDUM*CDDBLN(IEL))*CDI/D2PI
            IF (LDBNOD(IEL)) CFDBF=CFDBF-CDDPAD(iel)*CDI/D2PI
      ELSE
            CFDBF=-CMPLX(CDBZIV(1,IEL))
            DO 10 I=2,NDBTRM,2
            CDDUM=CMPLX(CDBZIV(I,IEL)-CDBZIV(I+1,IEL))
            CFDBF=CFDBF+DCOF(I)*CDDUM
  10        CONTINUE
            CFDBF=-CFDBF*CDI/D2PI
      ENDIF
      RETURN
      END
c
c ------------------------------------------------------------------------------------------------------------
c
      COMPLEX(8) FUNCTION CFDBG (IEL)
c
c ------------------------------------------------------------------------------------------------------------
c
C
C     POTENTIAL due to the G function for line doublet # IEL
C     Prior to this function the routine DBPREP(CZ) should have
C     been called to assign values to variables in common.
C
      IMPLICIT NONE
      INTEGER(4) IEL,I
      REAL(8) DCOF
      COMPLEX(8) CDDUM
      INCLUDE 'DBCOM.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
      DIMENSION DCOF(21)
      DATA DCOF /1.0D0,
     &           0.3333333333333333D0,0.3333333333333333D0,
     &           0.2000000000000000D0,0.2000000000000000D0,
     &           0.1428571428571428D0,0.1428571428571428D0,
     &           0.1111111111111111D0,0.1111111111111111D0,
     &           0.0909090909090909D0,0.0909090909090909D0,
     &           0.0769230769230769D0,0.0769230769230769D0,
     &           0.0666666666666666D0,0.0666666666666666D0,
     &           0.0588235294117647D0,0.0588235294117647D0,
     &           0.0526315789473684D0,0.0526315789473684D0,
     &           0.0476190476190476D0,0.0476190476190476D0/
C
      IF (.NOT.LDBFAR(IEL)) THEN
            CDDUM=CMPLX(CDBZEL(IEL))+1.0D0
            CFDBG=-(1.0D0+0.5D0*CDDUM*CDDBLN(IEL))*CDI/D2PI
      ELSE
            CFDBG=CMPLX(CDBZIV(1,IEL))
            DO 10 I=2,NDBTRM,2
            CDDUM=CMPLX(CDBZIV(I,IEL)+CDBZIV(I+1,IEL))
            CFDBG=CFDBG+DCOF(I)*CDDUM
  10        CONTINUE
            CFDBG=CFDBG*CDI/D2PI
      ENDIF
      RETURN
      END
c
c ------------------------------------------------------------------------------------------------------------
c
      COMPLEX(8) FUNCTION CFDBS (IEL)
c
c ------------------------------------------------------------------------------------------------------------
c
C
C     POTENTIAL due to the S function for line doublet # IEL (quadratic term)
C     Prior to this function the routine DBPREP(CZ) should have
C     been called to assign values to variables in common.
C
      IMPLICIT NONE
      INTEGER(4) IEL,I
      REAL(8) DCOF
      COMPLEX(8) CDDUM,CBZ
      INCLUDE 'DBCOM.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
      DIMENSION DCOF(21)
      DATA DCOF /0.3333333333333333D0,0.0D0,
     &           0.0666666666666666D0,0.0D0,
     &           0.0285714285714285D0,0.0D0,
     &           0.0158730158730158D0,0.0D0,
     &           0.0101010101010101D0,0.0D0,
     &           0.0069930069930069D0,0.0D0,
     &           0.0051282051282051D0,0.0D0,
     &           0.0039215686274500D0,0.0D0,
     &           0.0030959752321981D0,0.0D0,
     &           0.0025062656641604D0,0.0D0,
     &           0.0020703933747412D0/
C
      IF (.NOT.LDBFAR(IEL)) THEN
            CBZ=CDBZEL(IEL)
            CDDUM=CMPLX(CBZ*CBZ)-1.0D0
            CFDBS=(CMPLX(CBZ)+0.5D0*CDDUM*CDDBLN(IEL))*CDI/DPI
      ELSE
            CDDUM=CMPLX(CDBZIV(1,IEL))
            CFDBS=DCOF(1)*CDDUM
            DO 10 I=3,NDBTRM,2
            CDDUM=CMPLX(CDBZIV(I,IEL))
            CFDBS=CFDBS+DCOF(I)*CDDUM
  10        CONTINUE
            CFDBS=CFDBS*2.0D0*CDI/DPI
      ENDIF
      RETURN
      END
c
c ------------------------------------------------------------------------------------------------------------
c
      COMPLEX(8) FUNCTION CFDBFD (IEL)
c
c ------------------------------------------------------------------------------------------------------------
c
C
C     DERIVATIVE of the F function for line doublet # IEL
C     (Function W=-dOmega/dz)
C     Prior to this function the routine DBPREP(CZ) should have
C     been called to assign values to variables in common.
C
      IMPLICIT NONE
      INTEGER(4) IEL,I
      REAL(8) DCOF
      INCLUDE 'DBCOM.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
      DIMENSION DCOF(21)
      DATA DCOF /
     &             0.0000000000000000D+00, 0.1000000000000000D+01,
     &            -0.6666666666666666D+00, 0.1000000000000000D+01,
     &            -0.8000000000000000D+00, 0.1000000000000000D+01,
     &            -0.8571428571428571D+00, 0.1000000000000000D+01,
     &            -0.8888888888888888D+00, 0.1000000000000000D+01,
     &            -0.9090909090909091D+00, 0.1000000000000000D+01,
     &            -0.9230769230769231D+00, 0.1000000000000000D+01,
     &            -0.9333333333333333D+00, 0.1000000000000000D+01,
     &            -0.9411764705882353D+00, 0.1000000000000000D+01,
     &            -0.9473684210526315D+00, 0.1000000000000000D+01,
     &            -0.9523809523809523D+00/
C
      IF (.NOT.LDBFAR(IEL)) THEN
            CFDBFD=CMPLX(CDBZEL(IEL))-1.0D0
            CFDBFD=-CFDBFD/(CMPLX(CDBZEL(IEL))+1.0D0)
            CFDBFD=(CFDBFD+CDDBLN(IEL)+1.0D0)
            CFDBFD=-CFDBFD*CDI/D2PI/CMPLX(CDBZ21(IEL))
      ELSE
            CFDBFD=0.0D0
            DO 10 I=2,NDBTRM
            CFDBFD=CFDBFD+DCOF(I)*CMPLX(CDBZIV(I,IEL))
  10        CONTINUE
            CFDBFD=CFDBFD*CDI/DPI/CMPLX(CDBZ21(IEL))
      ENDIF
      RETURN
      END
c
c ------------------------------------------------------------------------------------------------------------
c
      COMPLEX(8) FUNCTION CFDBGD (IEL)
c
c ------------------------------------------------------------------------------------------------------------
c
C
C     DERIVATIVE of the G function for line doublet # IEL
C     (Function W=-dOmega/dz)
C     Prior to this function the routine DBPREP(CZ) should have
C     been called to assign values to variables in common.
C
      IMPLICIT NONE
      INTEGER(4) IEL,I
      REAL(8) DCOF
      INCLUDE 'DBCOM.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
      DIMENSION DCOF(21)
      DATA DCOF /
     &             0.0000000000000000D+00,-0.1000000000000000D+01,
     &            -0.6666666666666666D+00,-0.1000000000000000D+01,
     &            -0.8000000000000000D+00,-0.1000000000000000D+01,
     &            -0.8571428571428571D+00,-0.1000000000000000D+01,
     &            -0.8888888888888888D+00,-0.1000000000000000D+01,
     &            -0.9090909090909091D+00,-0.1000000000000000D+01,
     &            -0.9230769230769231D+00,-0.1000000000000000D+01,
     &            -0.9333333333333333D+00,-0.1000000000000000D+01,
     &            -0.9411764705882353D+00,-0.1000000000000000D+01,
     &            -0.9473684210526315D+00,-0.1000000000000000D+01,
     &            -0.9523809523809523D+00/
C
      IF (.NOT.LDBFAR(IEL)) THEN
            CFDBGD=CMPLX(CDBZEL(IEL))+1.0D0
            CFDBGD=CFDBGD/(CMPLX(CDBZEL(IEL))-1.0D0)
            CFDBGD=CFDBGD+CDDBLN(IEL)-1.0D0
            CFDBGD=CDI*CFDBGD/D2PI/CMPLX(CDBZ21(IEL))
      ELSE
            CFDBGD=0.0D0
            DO 10 I=2,NDBTRM
            CFDBGD=CFDBGD+DCOF(I)*CMPLX(CDBZIV(I,IEL))
  10        CONTINUE
            CFDBGD=-CFDBGD*CDI/DPI/CMPLX(CDBZ21(IEL))
      ENDIF
      RETURN
      END
c
c ------------------------------------------------------------------------------------------------------------
c
      COMPLEX(8) FUNCTION CFDBSD (IEL)
c
c ------------------------------------------------------------------------------------------------------------
c
C
C     DERIVATIVE of the S function for line doublet # IEL (quadratic term)
C     (Function W=-dOmega/dz)
C     Prior to this function the routine DBPREP(CZ) should have
C     been called to assign values to variables in common.
C
      IMPLICIT NONE
      INTEGER(4) IEL,I
      REAL(8) DCOF
      COMPLEX(8) CBZ
      INCLUDE 'DBCOM.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
      DIMENSION DCOF(21)
      DATA DCOF /
     &             0.0000000000000000D+00,-0.3333333333333333D+00,
     &             0.0000000000000000D+00,-0.2000000000000000D+00,
     &             0.0000000000000000D+00,-0.1428571428571428D+00,
     &             0.0000000000000000D+00,-0.1111111111111111D+00,
     &             0.0000000000000000D+00,-0.9090909090909091D-01,
     &             0.0000000000000000D+00,-0.7692307692307693D-01,
     &             0.0000000000000000D+00,-0.6666666666666667D-01,
     &             0.0000000000000000D+00,-0.5882352941176471D-01,
     &             0.0000000000000000D+00,-0.5263157894736842D-01,
     &             0.0000000000000000D+00,-0.4761904761904762D-01,
     &             0.0000000000000000D+00/
C
      IF (.NOT.LDBFAR(IEL)) THEN
            CBZ=CDBZEL(IEL)
            CFDBSD=(2.0D0+CMPLX(CBZ)*CDDBLN(IEL))
            CFDBSD=-2.0D0*CDI*CFDBSD/DPI/CMPLX(CDBZ21(IEL))
      ELSE
            CFDBSD=0.0D0
            DO 10 I=2,NDBTRM,2
            CFDBSD=CFDBSD+DCOF(I)*CMPLX(CDBZIV(I,IEL))
  10        CONTINUE
            CFDBSD=-CFDBSD*4.0D0*CDI/DPI/CMPLX(CDBZ21(IEL))
      ENDIF
      RETURN
      END
c
c ------------------------------------------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFDBPT (CZ)
c
c ------------------------------------------------------------------------------------------------------------
c
C
C     Function returns the real potential inside one or more inhomogeneities
C     due to added recharge.
C
      IMPLICIT NONE
      INTEGER(4) ISTR
      REAL(8) RFDBRP
      COMPLEX(8) CZ
      INCLUDE 'DBCOM.INC'
      INCLUDE 'LUSYS.INC'
C
      RFDBPT=0.0
      IF (NDB.EQ.0) RETURN
      DO 10 ISTR=1,NDBSTR
      IF (RDBGAM(ISTR).EQ.0.0) GOTO 10
      IF (LDBINS(ISTR)) RFDBPT=RFDBPT+RFDBRP(CZ,ISTR)
  10  CONTINUE
      RETURN
      END
c
c ------------------------------------------------------------------------------------------------------------
c
      REAL FUNCTION RFDBRP (CZ,ISTR)
c
c ------------------------------------------------------------------------------------------------------------
c
C      
C     Function returns the inside potential due added exfiltration for the
C     inhomogeneity ISTR. Note: recharge is negative.
C
      IMPLICIT NONE
      INTEGER(4) ISTR
      REAL(8) RAD
      COMPLEX(8) CZ
      INCLUDE 'DBCOM.INC'
      INCLUDE 'LUSYS.INC'
C
      RFDBRP=0.0
      IF (ISTR.LE.0.OR.ISTR.GT.NDBSTR) THEN
      WRITE (ILUER,1000) ISTR,NDBSTR
      RETURN
      ENDIF
      RAD=ABS(CZ-CDBZ0(ISTR))
      RAD=RAD*RAD
      RFDBRP=0.25*RDBGAM(ISTR)*(RAD-RDBRAD(ISTR))
      RETURN
 1000 FORMAT (' *** ERROR in RFDBRP: ISTR out of range.',/,
     &' ISTR is ',I3,' Current number of strings is ',I3)
      END

