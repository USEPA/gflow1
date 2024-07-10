C     Last change:  HMH  18 Dec 2006    2:21 pm
c     This file contains the following routines and functions
c
c    SUBROUTINE DINEAR    aborts streamline tracing when at 3D sink disc
c
c
c
c -----------------------------------------------------------
c
      SUBROUTINE DINEAR (CZ,CZNEW,R3DZ,RZNEW,LREDO)
c
c -----------------------------------------------------------
c
C
C     Routine checks for point near disc and adjusts stepsize
C     and sets L3DEND when end of pathline
C
C     INPUT:
C     CZ    current point of pathline
C     RVI1  unit vector in the direction of flow at CZ
C     LREDO some module wants to redo the calculation of CZNEW & RZNEW
C
C     Routine is called in PREDCOR
C     No action is taken in case stepsize has already been adjusted
C     in another L..NEAR routine, or in case L3DEND is already true.
C
      IMPLICIT NONE
      INTEGER(4) I,IDICLOSE,ID
      LOGICAL L3DEND,L1,L2,L3DREV,lredo
      REAL(8) RSD0,RSTEP,R3DZ,RS,RDIST,RZDIS,RZNEW,R1,RFHGHT,RAD
      COMPLEX(8) CZ,CZ0,CZNEW
      INCLUDE 'DICOM.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
C
      DATA IDICLOSE /0/
      IF (NDIS.EQ.0.or.lredo) RETURN
      CALL GETSTEP (RSD0,RSTEP,L3DEND,L3DREV) ! get data
      IF (L3DEND) RETURN     ! streamline ended elsewhere
      ID=1
c      write (iluout,1001) idiclose,cz,cznew,r3dz,rznew
      IF (IDICLOSE.NE.0) THEN       ! check if still near disc IDICLOSE --(1)
c      write (iluout,1001) idiclose,cz,cznew,r3dz,rznew
c 1001 format (' dinear1: idiclose,cz,cznew ',i5,4g11.4,/,
c     &' r3dz,rznew ',2g11.4)
        CZ0=CMPLX(RDIZ(1,IDICLOSE),RDIZ(2,IDICLOSE))
        RTAU=ABS(CZ0-CZNEW)
        RAD=RDIR(IDICLOSE)
c        write (iluout,1003) cz0,rtau,rad,rsd0
        IF (RTAU.GT.RAD+RSD0) GOTO 5         ! outside disc area -----------(3)
c        write (iluout,1003) cz0,rtau,rad,rsd0
c 1003   format (' dinear3: cz0,rtau,rad,rsd0 ',5g11.4)                
        RDIST=ABS(RZNEW-RDIZ(3,IDICLOSE))
        IF (RDIST.GT.RSD0) GOTO 5  ! not near disc
        GOTO 15
      ENDIF
  5   IF (ID.GT.NDIS) RETURN        ! no more discs to check
        DO 10 I=ID,NDIS
        CZ0=CMPLX(RDIZ(1,I),RDIZ(2,I))
        RTAU=ABS(CZ0-CZNEW)
        RAD=RDIR(I)
c        write (iluout,1002) i,cz0,rtau,rad,rsd0        
        IF (RTAU.GT.RAD+RSD0) GOTO 10         ! outside disc area ------(2)
c 1002   format (' dinear2: i,cz0,rtau,rad,rsd0 ',i3,5g11.4)                
        RDIST=ABS(RZNEW-RDIZ(3,I))
        IF (RDIST.GT.RSD0) GOTO 10  ! not near disc
c        write (iluout,1002) i,cz0,rtau,rad,rsd0        
        IDICLOSE=I
        ID=I+1
        GOTO 15
  10    CONTINUE      
      RETURN                        ! not close to any disc
  15  RS=RDIS(IDICLOSE)
c      write (iluout,1004) idiclose,rs
      L1=.NOT.L3DREV.AND.RS.LT.0.0
      L2=L3DREV.AND.RS.GT.0.0
      IF (RS.EQ.0.0.OR.L1.OR.L2) GOTO 5  ! cannot end at disc -------------(4)
c      write (iluout,1004) idiclose,rs
c 1004 format (' dinear4: idiclose,rs ',i5,g14.7)      
      R1=0.1*RFHGHT(CZNEW)          ! stepsize reduced to 0.1 * aquifer height
      RSTEP=MIN(R1,RSTEP)      ! make sure existing stepsize is not enlarged
      RAD=RDIR(IDICLOSE)
c      write (iluout,1005) idiclose,rtau,rad
      IF (RTAU.LE.RAD) THEN ! ----------------------------------------------(5)
c      write (iluout,1005) idiclose,rtau,rad
c 1005 format (' dinear5: idiclose,rtau,rad ',i5,2g11.4)      
        RZDIS=RDIZ(3,IDICLOSE)
        L1=(RZDIS-RZNEW)*(RZDIS-R3DZ).LT.0.0 ! passed through disc
        RDIST=ABS(RZDIS-R3DZ)
        L2=RDIST.LT.0.1*RSTEP   ! arrived at disc
c        write (iluout,1006) l1,l2,rdist
        IF (L1.OR.L2) THEN   ! crossed disc -------------------------------(6)
c        write (iluout,1006) l1,l2,rdist
c 1006   format (' dinear6: l1,l2,rdist ',2l5,g11.4)        
          L3DEND=.TRUE.
          iFlag=2
          iElementType=6
          IF (LEN(ADILAB(IDICLOSE)).GT.0)
     &        aElementLabel=ADILAB(IDICLOSE)
          RZNEW=RZDIS
          RSTEP=RDIST           
        ENDIF
      ELSE
c        write (iluout,1007) rtau,rad,rstep
        IF (RTAU.LE.RAD+RSTEP) THEN   ! ------------------------------------(7)
c          write (iluout,1007) rtau,rad,rstep
c 1007     format (' dinear7: rtau,rad,rstep ',3g14.7)          
          L3DEND=.TRUE.
          iFlag=2
          iElementType=6
          IF (LEN(ADILAB(IDICLOSE)).GT.0)
     &        aElementLabel=ADILAB(IDICLOSE)
          RZNEW=RZDIS
          RSTEP=ABS(CZ-CZ0)-RAD
        ENDIF
      ENDIF
      CALL SETSTEP (RSTEP,L3DEND)
      IF (.NOT.L3DEND) GOTO 5 ! see if close to another disc
      RETURN
      END

