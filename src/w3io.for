C     Last change:  HMH   8 Aug 2005    6:47 pm
C     This file cotains the following routines and functions
C
C     SUBROUTINE W3EXTRACT        print data for partially penetrating wells (NOT USED).
C     SUBROUTINE w3extracthead    writes data for head specified wells to .xtr file
C     SUBROUTINE w3extractdisch   writes data for discharge specified wells to .xtr file
C     SUBROUTINE W3DATIO          rites all well data to a ".dat" file
C     SUBROUTINE W3IO             reads or writes contents of WLCOM.INC from or to .sol file
C     SUBROUTINE W3OUT            writes contents of WLCOM formatted to ILUOUT for debugging purposes
C
c
c ---------------------------------------------------------------------------------------------
c
           SUBROUTINE W3EXTRACT (ILU)    ! CURRENTLY NOT USED
c
c ---------------------------------------------------------------------------------------------
c
C
C     Print data for partially penetrating wells.
C
      IMPLICIT NONE
      INTEGER(4) ILU,IEN,IST,IW,I,NPT
      REAL(8) RQ,RHDSUM,RBAS,RFBASE,RHDI,RFHEAD,RERR,RHDAV,RHEDS
      COMPLEX(8) CZ
      INCLUDE 'w3com.inc'
      INCLUDE 'com3d.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
      DIMENSION RHDI(100)
      IF (NW3.EQ.0) RETURN
      DO 50 IW=1,NW3
      IEN=IPNT(IW+1)-1
      IST=IPNT(IW)
      RQ=0.0
      RHDSUM=0.0
      NPT=0
      cz=CMPLX(rw3cp(1,ist),rw3cp(2,ist))
      rbas=rfbase(cz)
      DO 10 I=IEN,IST,-1      ! calculate heads and discharge
      CZ=CMPLX(RW3CP(1,I),RW3CP(2,I))
      R3DZ=RW3CP(3,I)
      RQ=RQ+RW3S(I)*RW3L(I)
      NPT=NPT+1
      RHDI(NPT)=RFHEAD(CZ)
      RHDSUM=RHDSUM+RHDI(NPT)
  10  CONTINUE      
      IF (LW3Q(IW)) THEN            ! discharge specified well
      RERR=100.0
      IF (RW3Q(IW).NE.0.0) RERR=ABS((RW3Q(IW)-RQ)/RW3Q(IW))*100.0
      RHDAV=RHDSUM/NPT
      WRITE (ILU,1000) IW,RHDAV,RW3Q(IW),RQ,RERR
      NPT=0.0
      WRITE (ILU,4000)
      DO 20 I=IEN,IST,-1
      NPT=NPT+1
      RERR=100.0
      IF (RHDAV.NE.0.0) RERR=ABS((RHDI(NPT)-RHDAV)/(RHDAV-rbas))*100.0
      WRITE (ILU,4050) RW3CP(1,I),RW3CP(2,I),RW3CP(3,I),RHDI(NPT)
  20  CONTINUE
      ELSE                    ! head specified well
      WRITE (ILU,1100) IW,RQ
      WRITE (ILU,3100)
      RHEDS=RW3HED(IW)
      NPT=0
      DO 30 I=IEN,IST,-1
      NPT=NPT+1
      RERR=100.0
      IF (RHEDS.NE.0.0) RERR=ABS((RHEDS-RHDI(NPT))/(RHEDS-rbas))*100.0
      WRITE (ILU,4100) RW3CP(1,I),RW3CP(2,I),RW3CP(3,I),RHEDS,
     &RHDI(NPT),RERR
  30  CONTINUE
      ENDIF
      WRITE (ILU,5000)
      IF (LW3DRT(IW)) THEN
      WRITE (ILU,6000) RW3S(IEN)
      IEN=IEN-1
      ENDIF
      write (ilu,7000)
      DO 40 I=IEN,IST,-1
      WRITE (ILU,7050) RW3ST(3,I),RW3S(I)
  40  CONTINUE
  50  CONTINUE
      RETURN
C
 1000 FORMAT ('! partially penetrating well  (discharge specified) ',/,
     &        '* well number average calc. head  spec. disch.  ',
     &        'calc. disch     %error',/,I5,10x,4(E14.7,2x))
 1100 FORMAT ('! partially penetrating well  (head specified)',/,
     &        '* well number    calc. discharge   ',/,i5,10x,e14.7)
 3100 FORMAT (
     &'*    x                y              z           specified head '
     &,' calculated head  error [%]')
 4000 FORMAT ('*     x              y             z             ',
     &'calc. head')
 4050 FORMAT (5(E14.7,2X))
 4100 FORMAT (6(E14.7,2X))
 5000 FORMAT ('* strength parameter of double root function ')
 6000 FORMAT (e14.7)
 7000 FORMAT ('* z-value     linear strength parameter')
 7050 FORMAT (2(E14.7,2X))
      END
c
c ---------------------------------------------------------------------------------------------
c
           SUBROUTINE W3EXTRACTHEAD (ILU)
c
c ---------------------------------------------------------------------------------------------
c
C
C     Print data for head specified partially penetrating wells.
C
      IMPLICIT NONE
      INTEGER(4) ILU,IEN,IST,IW,NPT,I
      LOGICAL LW3HEAD
      REAL(8) RQ,RHDSUM,RBAS,RFBASE,RHDI,RFHEAD,RHED,RHEDS,RERR
      COMPLEX(8) CZ
      INCLUDE 'w3com.inc'
      INCLUDE 'com3d.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
      DIMENSION RHDI(100)
      if (lw3head()) then
      write (ilu,2000)
      DO 50 IW=1,NW3
      IEN=IPNT(IW+1)-1
      IST=IPNT(IW)
      RQ=0.0
      RHDSUM=0.0
      NPT=0
      cz=CMPLX(rw3cp(1,ist),rw3cp(2,ist))
      rbas=rfbase(cz)
      DO 10 I=IEN,IST,-1      ! calculate heads and discharge
      CZ=CMPLX(RW3CP(1,I),RW3CP(2,I))
      R3DZ=RW3CP(3,I)
      RQ=RQ+RW3S(I)*RW3L(I)
      NPT=NPT+1
      RHDI(NPT)=RFHEAD(CZ)
      RHDSUM=RHDSUM+RHDI(NPT)
  10  CONTINUE      
      IF (.NOT.LW3Q(IW)) THEN ! HEAD specified well, print data
        RHED=RHDSUM/NPT
        RHEDS=RW3HED(IW)
        NPT=0
        DO 20 I=IEN,IST,-1
        NPT=NPT+1
        RERR=100.0
        IF (RHEDS.NE.0.0) RERR=ABS((RHEDS-RHDI(NPT))/(RHEDS-rbas))*100.0
   20   CONTINUE
        WRITE (ILU,3000) CZ,RW3RAD(IW),RQ,RHED,RERR,AW3LAB(IW)
      ENDIF
  50  CONTINUE
      endif
      RETURN
C
 2000 format ('! head specified partially penetrating wells',/,
     & '*     x1            y1          radius         Q',
     & '             head            % error in head   label')
 3000 format (6(e14.7,','),a16)
      END
c
c ---------------------------------------------------------------------------------------------
c
           SUBROUTINE W3EXTRACTDISCH (ILU)
c
c ---------------------------------------------------------------------------------------------
c
C
C     Print data for discharge specified partially penetrating wells. (called in WLEXTRACTDISCH)
C
      IMPLICIT NONE
      INTEGER(4) ILU,IW,IEN,IST,NPT,I
      LOGICAL LW3DISCHARGE
      REAL(8) RQ,RHDSUM,RBAS,RFBASE,RHDI,RFHEAD,RHED
      COMPLEX(8) CZ
      INCLUDE 'w3com.inc'
      INCLUDE 'com3d.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
      DIMENSION RHDI(100)
      if (lw3discharge()) then
      write (ilu,2000)
      DO 50 IW=1,NW3
      IEN=IPNT(IW+1)-1
      IST=IPNT(IW)
      RQ=0.0
      RHDSUM=0.0
      NPT=0
      cz=CMPLX(rw3cp(1,ist),rw3cp(2,ist))
      rbas=rfbase(cz)
      DO 10 I=IEN,IST,-1      ! calculate heads and discharge
      CZ=CMPLX(RW3CP(1,I),RW3CP(2,I))
      R3DZ=RW3CP(3,I)
      RQ=RQ+RW3S(I)*RW3L(I)
      NPT=NPT+1
      RHDI(NPT)=RFHEAD(CZ)
      RHDSUM=RHDSUM+RHDI(NPT)
  10  CONTINUE      
      IF (LW3Q(IW)) THEN ! DISCHARGE specified well, print data
        RHED=RHDSUM/NPT
        WRITE (ILU,3000) CZ,RW3RAD(IW),RQ,RHED,AW3LAB(IW)
      ENDIF
  50  CONTINUE
      endif
      RETURN
C
 2000 format ('! discharge specified partially penetrating wells',/,
     & '*     x1            y1          radius         Q',
     & '             head            label')

 3000 format (5(e14.7,','),a16)
      END
c
c ---------------------------------------------------------------------------------------------
c
      SUBROUTINE W3DATIO (ILU)
c
c ---------------------------------------------------------------------------------------------
c
C
C     Routine writes all partially penetrating well data to a ".dat" file.
C     Routine is called in DATIO.
C
      IMPLICIT NONE
      INTEGER(4) ILU,I,IBOT,ITOP
      LOGICAL LHEAD,LDISCHARGE
      REAL(8) RZB,RZT
      COMPLEX(8) CZ
      INCLUDE 'w3com.inc'
      INCLUDE 'lusys.inc'
C
      IF (NW3.EQ.0) RETURN
      WRITE (ILU,1000)
      LHEAD=.FALSE.
      LDISCHARGE=.FALSE.
      DO 10 I=1,NW3
      IF (LW3Q(I)) LDISCHARGE=.TRUE.
      IF (.NOT.LW3Q(I)) LHEAD=.TRUE.
  10  CONTINUE      
      IF (LHEAD) THEN
        WRITE (ILU,2000)
        DO 20 I=1,NW3
        IF (LW3Q(I)) GOTO 20
        IBOT=IPNT(I)
        ITOP=IPNT(I+1)-1
        IF (LW3DRT(I)) ITOP=ITOP-1
        CZ=CMPLX(RW3ST(1,IBOT),RW3ST(2,IBOT))
        RZB=RW3ST(3,IBOT)
        RZT=RW3ST(3,ITOP)
        WRITE (ILU,3000) CZ,RZB,RZT,RW3HED(I),RW3RAD(I),AW3LAB(I)
  20    CONTINUE
      ENDIF     
      IF (LDISCHARGE) THEN
        WRITE (ILU,4000)
        DO 30 I=1,NW3
        IF (.NOT.LW3Q(I)) GOTO 30
        IBOT=IPNT(I)
        ITOP=IPNT(I+1)-1
        IF (LW3DRT(I)) ITOP=ITOP-1
        CZ=CMPLX(RW3ST(1,IBOT),RW3ST(2,IBOT))
        RZB=RW3ST(3,IBOT)
        RZT=RW3ST(3,ITOP)
        WRITE (ILU,5000) CZ,RZB,RZT,RW3Q(I),RW3RAD(I),AW3LAB(I)
  30    CONTINUE      
      ENDIF
      WRITE (ILU,6000)
      RETURN
C      
 1000 FORMAT (' ppwell')
 2000 FORMAT (' head',/,'*    x        y      zbottom     ztop',
     &        '       head       radius    label')
 3000 FORMAT (4(F9.0,1X),G14.7,1X,F5.1,1X,A16)
 4000 FORMAT (' discharge',/,'*    x         y      zbottom    ztop',
     &        '     discharge    radius    label')
 5000 FORMAT (4(F9.0,1X),G14.7,1X,F5.1,1X,A16)     
 6000 FORMAT (' quit')
      END
c
c ---------------------------------------------------------------------------------------------
c
      SUBROUTINE W3IO (ICODE,ILU,RVERSION,ierr)
c
c ---------------------------------------------------------------------------------------------
c
C
C     Routine reads or writes contents of W3COM.INC common blocks to
C     an external file.
C     ICODE=41 write
C     ICODE=73 read
C
      IMPLICIT NONE
      INTEGER(4) ICODE,ILU,IERR,NPNT,NWORD
      REAL(8) RVERSION
      INCLUDE 'w3com.inc'
      INCLUDE 'lusys.inc'
C
      if (ierr.ne.0) return
      CALL BUFIO4 (NW3,1,ILU,ICODE,IERR)
      IF (NW3.EQ.0) RETURN
      CALL BUFIOR (RW3RAD,NW3,ILU,ICODE,IERR)
      CALL BUFIOR (RW3HED,NW3,ILU,ICODE,IERR)
      CALL BUFIOR (RW3Q,NW3,ILU,ICODE,IERR)
      CALL BUFIO4 (IPNT,NW3+1,ILU,ICODE,IERR)
      CALL BUFIOL (LW3Q,NW3,ILU,ICODE,IERR)
      CALL BUFIOL (LW3DRT,NW3,ILU,ICODE,IERR)
      CALL BUFIOL (LW3BOT,NW3,ILU,ICODE,IERR)
      CALL BUFIOL (LW3TOP,NW3,ILU,ICODE,IERR)
      CALL BUFIOA (AW3LAB,NW3,ILU,ICODE,IERR)                  
      NPNT=IPNT(NW3+1)
      NWORD=NPNT*3
      CALL BUFIOR (RW3ST,NWORD,ILU,ICODE,IERR)
      CALL BUFIOR (RW3CP,NWORD,ILU,ICODE,IERR)
      CALL BUFIOR (RW3S,NPNT,ILU,ICODE,IERR)
      CALL BUFIOR (RW3L,NPNT,ILU,ICODE,IERR)
      CALL BUFIO4 (NLIMAG,3,ILU,ICODE,IERR)
      CALL BUFIOR (RW3ALP,1,ILU,ICODE,IERR)
      IF (RVERSION.LE.10.0) RETURN
      CALL BUFIOR (DW3POT,NPNT,ILU,ICODE,IERR)
C
      RETURN
      END
c
c ---------------------------------------------------------------------------------------------
c
      SUBROUTINE W3OUT (ICALL)
c
c ---------------------------------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER(4) ICALL,IEN,I,J
      CHARACTER(1) ADUM
      INCLUDE 'w3com.inc'
      INCLUDE 'lusys.inc'
      SAVE
      IEN=IPNT(NW3+1)-1
      WRITE (ILUOUT,1000) NW3,NLIMAG,NW3RH,RW3ALP
      WRITE (ILUOUT,2000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000)
      READ (ILUIN,2500) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE (ILUOUT,1001) (I,IPNT(I),RW3RAD(I),RW3HED(I),RW3Q(I),
     &                     I=1,10)
      WRITE (ILUOUT,1002) (I,LW3Q(I),LW3DRT(I),LW3BOT(I),LW3TOP(I),
     &                     I=1,10)
      WRITE (ILUOUT,2000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000)
      READ (ILUIN,2500) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN    
      WRITE (ILUOUT,1003) (I,(RW3ST(J,I),J=1,3),I=1,IEN)
      WRITE (ILUOUT,2000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000)
      READ (ILUIN,2500) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN     
      WRITE (ILUOUT,1004) (I,(RW3CP(J,I),J=1,3),I=1,IEN)
      WRITE (ILUOUT,2000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000)
      READ (ILUIN,2500) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN           
      WRITE (ILUOUT,1005) (I,RW3S(I),I=1,IEN)
      WRITE (ILUOUT,2000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000)
      READ (ILUIN,2500) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN     
      WRITE (ILUOUT,1006) (I,RW3L(I),I=1,IEN)
      WRITE (ILUOUT,2000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000)
      READ (ILUIN,2500) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE (ILUOUT,1007) (I,AW3LAB(I),I=1,IEN)
      WRITE (ILUOUT,2000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000)
      READ (ILUIN,2500) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN      
 1000 FORMAT (' W3OUT:NW3,NLIMAG,NW3RH,RW3ALP ',3I4,G14.7)
 1001 FORMAT (' W3OUT:I,IPNT,RW3RAD,RW3HED,RW3Q ',2I3,3G11.4)
 1002 FORMAT (' W3OUT:I,LW3Q,LW3DRT,LW3BOT,LW3TOP ',I4,4L4)
 1003 FORMAT (' W3OUT:I,RW3ST(1-3,I) ',I3,3G14.7)
 1004 FORMAT (' W3OUT:I,RW3CP(1-3,I) ',I3,3G14.7) 
 1005 FORMAT (' W3OUT:I,RW3S(I) ',I3,G14.7)
 1006 FORMAT (' W3OUT:I,RW3L(I) ',I3,G14.7)
 1007 FORMAT (' W3OUT:I,AW3LAB(I) ',I3,1X,A16)
 2000 FORMAT (' W3OUT: ICALL=',I3,' press <Enter> to continue')
 2500 FORMAT (A1)
      RETURN
      END
    

