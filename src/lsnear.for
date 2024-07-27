c       This file contains the following routines and functions:
c
c	SUBROUTINE LSNEAR        handles particle tracing near line-sinks
c       LOGICAL FUNCTION LINSID  true if inside a box around the linesink
c
C
C----------------------------------------------------------------------------------------------------------
C
      SUBROUTINE LSNEAR(CZ,CZNEW,R3DZ,RZNEW,CZEULER,RZEULER,LEULER,
     &                  LREDO)
C
C----------------------------------------------------------------------------------------------------------
C
C     Last change:  HMH  18 Dec 2006    2:56 pm
C
C     Routine checks for point near line sink and adjusts stepsize,
C     prevents oscillating of pathline around the line sink, and
C     sets L3DEND when end of pathline.
C
C     INPUT:
C
C     CZ       current point of pathline
C     CZNEW    projected position of next point
C     R3DZ     elevation of current point
C     RZNEW    projected elevation of next point
C     LREDO    some module wants to recalculate CZNEW & RZNEW
C
C     OUTPUT:
C
C     CZNEW    position of new pathline point
C     RZNEW    elevation of new pathline point
C     RSTEP    distance between current and new point
C
      IMPLICIT NONE
      INTEGER(4) I,ILSCLOSE,ILS
      REAL(8) R3DZ,RBASE,RFBASE,RTOP,RHGHT,RFHGHT,RWIDTH,RZONEL,
     &     RZONEW,RSD0,RSTEP,REAL,AIMAG,RF3DSP,RVI2,RS,
     &     RY0,RY0NEW,RYY1,RYY2,RX0,RX0NEW,RX,RNI,RQNMIN,
     &     RQNPLUS,RZ,RZNEW,R1,R2,RDZ,RLENGTH,RZEULER,RZNEW1
      COMPLEX(8) CZ,CZ1,CZ2,CZNEW,CFBIGZ,CBIGZ,CZ0,CZEULER,
     &        CBIGZNEW,CUNIT,CBIG1,CBIG2,CI,CN,CS,CZZ0,CZNEW1
      LOGICAL L3DEND,L3DREV,LINSID,L1,L2,LJUMP,LEULER,lredo
      INCLUDE 'lscom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      EXTERNAL VELOC
      DIMENSION RVI2(3),RNI(3)
      DATA CI /(0.0,1.0)/
      DATA ILSCLOSE,RVI2 /0,3*0.0/
      DATA LJUMP /.FALSE./
      IF (NLS.EQ.0) RETURN
      LCLOSE=.FALSE.
      LEULER=.FALSE.
      ILS=1
      CALL GETSTEP (RSD0,RSTEP,L3DEND,L3DREV)   ! Get data
      IF (L3DEND) RETURN                        ! streamline ended elsewhere
      RHGHT=RFHGHT(CZ)
      RBASE=RFBASE(CZ)
      RTOP=RBASE+RHGHT
c      write (iluout,1001) cz,cznew,
c     &r3dz,rznew,rsd0,rstep,l3dend,ilsclose
c 1001 format (' lsnear1: cz,cznew ',4g14.7,/,
c     &' r3dz,rznew,rsd0,rstep,l3dend,ilsclose ',/,
c     &4G11.4,l3,i3)            
      IF (ILSCLOSE.NE.0) THEN ! check if still near previous line sink  ----(1)
c      write (iluout,1001) cz,cznew,
c     &r3dz,rznew,rsd0,rstep,l3dend,ilsclose      
        CZ1=CLSZS(ILSCLOSE)
        CZ2=CLSZE(ILSCLOSE)
        RWIDTH=RLSWID(ILSCLOSE)
        RZONEL=RSD0                 ! "alert box" (2*rsd0+rwidth) wide
        RZONEW=0.5*RWIDTH+RSD0      !             (2*rsd0+length) long
c        write (iluout,1003) ilsclose        
        IF (LINSID(CZ,RZONEW,RZONEL,CZ1,CZ2)) THEN ! -----------------(3)
c          write (iluout,1003) ilsclose        
          IF (LJUMP) THEN           ! jump already taken in previous call
            LJUMP=.FALSE.
            GOTO 50
          ENDIF
          GOTO 15
c 1003     format (' lsnear3: ilsclose ',i3)                
        ELSE
          ILSCLOSE=0
          LJUMP=.FALSE.
        ENDIF
      ENDIF
  5   CONTINUE    
c      write (iluout,1002) ils,ilsclose,cz,cznew
c 1002 format (' lsnear2: ils,ilsclose,cz,cznew ',i3,2x,i3,4g11.4)       
       DO 10 I=ILS,NLS                                            ! Loop over all line sinks
       IF (I.EQ.ILSCLOSE) GOTO 10
       CZ1=CLSZS(I)
       CZ2=CLSZE(I)
       CBIG1=CFBIGZ(CZ1,CZ,CZNEW)
       CBIG2=CFBIGZ(CZ2,CZ,CZNEW)
       RYY1=AIMAG(CBIG1)
       RYY2=AIMAG(CBIG2)
       IF (RYY1*RYY2.LT.0.0) THEN       ! crossing the line sink line ----(2)
          RWIDTH=RLSWID(I) ! test only for closeness if path crosses line sink
          RZONEL=RSD0
          RZONEW=0.5*RWIDTH+RSD0
          IF (LINSID(CZ,RZONEW,RZONEL,CZ1,CZ2)) THEN ! Close to line sink 
            ILSCLOSE=I
            ILS=I+1
c            write (iluout,1002) ils,ilsclose,cz,cznew
            GOTO 15
          ENDIF
        ENDIF
  10    CONTINUE
c        write (iluout,1111) cz,cznew,r3dz,rznew,rstep
        RETURN      ! not close to any line sink, return
  15  RS=RLSIG(ILSCLOSE)
c      write (iluout,1004) ilsclose,rs
c 1004 format (' lsnear4: ilsclose,rs ',i3,g14.7)      
      IF (RS.EQ.0.0) THEN   ! line sink not active, return ---------------(4)
c        write (iluout,1004) ilsclose,rs      
c        write (iluout,1111) cz,cznew,r3dz,rznew,rstep
        GOTO 5  ! also close to other line sink?
      ENDIF
C      
c      RZONEL=0.0 
c      RZONEW=0.5*RWIDTH
c      CALL PLOTWIDTH (RZONEW+RSD0,RZONEL+RSD0,CZ1,CZ2) ! debugging
c      CALL DRAW (CZ,0)                                 ! 
c      CALL PLOTWIDTH (RZONEW,RZONEL,CZ1,CZ2)           !
c      CALL DRAW (CZ,0)                                 ! debugging
C      
 40   CBIGZ=CFBIGZ(CZ,CZ1,CZ2)
      CBIGZNEW=CFBIGZ(CZNEW,CZ1,CZ2)
      RY0=AIMAG(CBIGZ)
      RY0NEW=AIMAG(CBIGZNEW)
      L1=RY0*RY0NEW.LE.0.0
      CBIG1=CFBIGZ(CZ1,CZ,CZNEW)
      CBIG2=CFBIGZ(CZ2,CZ,CZNEW)
      RYY1=AIMAG(CBIG1)
      RYY2=AIMAG(CBIG2)      
      L2=RYY1*RYY2.LE.0.0
c      write (iluout,1005) L1,L2,LEULER
c 1005 format (' lsnear5: l1,l2,leuler ',3l5)
      IF (L1.AND.L2) THEN  ! intersecting the line sink -------------------(5)
c        write (iluout,1005) L1,L2,LEULER      
C       ------------------------------------------------calculate CZ0
        RX0=REAL(CBIGZ)
        RX0NEW=REAL(CBIGZNEW)
        RX=(RX0NEW*RY0-RX0*RY0NEW)/(RY0-RY0NEW)
        CZ0=0.5*(RX*(CZ2-CZ1)+CZ1+CZ2)
C       --------------------------------------calculate RQNMIN and RQNPLUS
        RLENGTH=ABS(CZ2-CZ1)
        CUNIT=(CZ2-CZ1)/RLENGTH
        CN=CI*CUNIT
        RNI(1)=REAL(CN)
        RNI(2)=AIMAG(CN)
        RNI(3)=0.0
        CS=CZNEW-CZ
        CS=CS/ABS(CS)
        CZZ0=CZ0+0.001*RLENGTH*CS
        CALL DISCH (CZZ0,RVI2)
        RQNPLUS=RF3DSP(RNI,RVI2)
c        write (iluout,1051) czz0,rvi2
c 1051   format (' lsnear5a: czz0,rvi2 ',5g11.4)        
        CZZ0=CZ0-0.001*RLENGTH*CS
        CALL DISCH (CZZ0,RVI2)
        RQNMIN=RF3DSP(RNI,RVI2)
c        write (iluout,1052) czz0,rvi2
c 1052   format (' lsnear5b: czz0,rvi2 ',5g11.4)                
c        write (iluout,1006) cz0,rqnmin,rqnplus
        IF (RQNMIN*RQNPLUS.LE.0.0) THEN ! ---------------------------------(6)
c          write (iluout,1006) cz0,rqnmin,rqnplus
c 1006     format (' lsnear6: cz0,rqnmin,rqnplus ',4g14.7)
            L1=.NOT.L3DREV.AND.RS.GT.0.0
            L2=L3DREV.AND.RS.LT.0.0
c            write (iluout,1007) l3drev,l1,l2,rs
c 1007       format (' lsnear7: l3drev,l1,l2,rs ',3l5,g14.7)            
            IF (L1.OR.L2) THEN      ! end at line sink --------------(7)
c              write (iluout,1007) l3drev,l1,l2,rs
              RZNEW=RTOP
              RSTEP=ABS(CZ-CZ0)
              CZNEW=CZ0
              L3DEND=.TRUE.
              iFlag=2
              iElementType=2
              IF (LEN(ALSLAB(ILSCLOSE)).GT.0)
     &        aElementLabel=ALSLAB(ILSCLOSE)
              GOTO 50            
            ELSE
              RDIRI(1)=REAL(CUNIT)   ! store line sink orientation
              RDIRI(2)=AIMAG(CUNIT)  ! velocity component parallel to
              RDIRI(3)=0.0           ! line sink should be used
              LCLOSE=.TRUE.
              GOTO 50
            ENDIF
        ELSE
C         ---------------------------------------------calculate R3DZ        
          R1=ABS(CZNEW-CZ)
          R2=ABS(CZ0-CZ)
          R3DZ=R3DZ+R2/R1*(RZNEW-R3DZ)
          RZ=R3DZ-RBASE
          RDZ=RQNMIN/RQNPLUS*RZ
          RZNEW=RBASE+RDZ
c          write (iluout,1008) rz,rbase,rtop,rznew
          IF (RZNEW.GT.RTOP) THEN ! ---------------------------------------(8)
c            write (iluout,1008) rz,rbase,rtop,rznew
c 1008       format (' lsnear8: rz,rbase,rtop,rznew ',4g14.7)
            RZNEW=RTOP
            L3DEND=.TRUE.
            iFlag=2
            iElementType=2
            IF (LEN(ALSLAB(ILSCLOSE)).GT.0)
     &      aElementLabel=ALSLAB(ILSCLOSE)
            CZNEW=CZ0
            RSTEP=ABS(CZ0-CZ)
            GOTO 50
          ELSE
            LJUMP=.TRUE.
            GOTO 50
          ENDIF
        ENDIF
      ELSE
        IF (.NOT.LEULER) THEN   ! calculate new L1 and L2 based on Euler step
          CZNEW1=CZNEW
          RZNEW1=RZNEW
          CZNEW=CZEULER
          RZNEW=RZEULER
          LEULER=.TRUE.
          GOTO 40
        ELSE
          CZNEW=CZNEW1
          RZNEW=RZNEW1
          LEULER=.FALSE.
          GOTO 5    ! also close to other line sink?
        ENDIF
      ENDIF
  50  continue
      CALL SETSTEP (RSTEP,L3DEND)
c      write (iluout,1111) cz,cznew,r3dz,rznew,rstep
c 1111 format (' leaving lsnear with: cz,cznew ',4g11.4,/,
c     & ' r3dz,rznew,rstep ',3g11.4)
      RETURN
      END
C
C----------------------------------------------------------------------------------------------------
C
      LOGICAL FUNCTION LINSID(CZ,RZONEW,RZONEL,CZ1,CZ2)
C
C----------------------------------------------------------------------------------------------------
C
C
C     Checks for CZ inside domain of width 2*RZONEW and 
C     length of 2*RZONEL surrounding the linesink between
C     CZ1 and CZ2.
C      
      IMPLICIT COMPLEX (C), LOGICAL (L)
      INCLUDE 'lusys.inc'
C      CHARACTER(1) ACHAR
      DIMENSION CZI(5)
      DATA CI /(0.0,1.0)/
C
      CZDIR=CZ2-CZ1
      CZDIR=CZDIR/ABS(CZDIR)
      CZI(1)=CZ1-CZDIR*(RZONEL+CI*RZONEW)
      CZI(2)=CZ2+CZDIR*(RZONEL-CI*RZONEW)
      CZI(3)=CZI(2)+CZDIR*(2.0*CI*RZONEW)
      CZI(4)=CZI(1)+CZDIR*(2.0*CI*RZONEW)
      CZI(5)=CZI(1)
      RTETA=0.0
      DO 10 I=1,4
      CDUM=(CZI(I)-CZ)/(CZI(I+1)-CZ)
      RX=REAL(CDUM)
      RY=AIMAG(CDUM)
      RDUM=ATAN2(RY,RX)
C      WRITE (ILUOUT,1004) CDUM
C 1004 FORMAT (' LINSID4: CDUM ',2G11.4)      
      RTETA=RTETA+RDUM
C      WRITE (ILUOUT,1003) RDUM,RTETA
C 1003 FORMAT (' LINSID3: RDUM, RTETA ',2G11.4)      
  10  CONTINUE
      LINSID=ABS(RTETA).GT.1.0
C      WRITE (ILUOUT,1001) (I,CZI(I),I=1,5)
C 1001 FORMAT (' LINSID1: CZI(',I1,')=',2G11.4)
C      WRITE (ILUOUT,1002) CZ,CZDIR,RZONEW,RZONEL,RTETA,LINSID
C 1002 FORMAT (' LINSID2: CZ,CZDIR ',4G11.4,/,
C     &        ' LINSID2: RZONEW,RZONEL,RTETA,LINSID ',3G11.4,L5)
C      READ (ILUIN,1000) ACHAR
C      IF (ACHAR.EQ.'Q'.OR.ACHAR.EQ.'q') RETURN
C 1000 FORMAT (A1)            
      RETURN
      END


