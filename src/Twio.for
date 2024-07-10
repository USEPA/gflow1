C     Last change:  HMH   8 Aug 2005    6:47 pm
C     This file contains the following routines and functions
C
C     SUBROUTINE twextract        writes data for transient wells to .xtr file
C     SUBROUTINE TWDATIO          writes all transient well data to a ".dat" file
C     SUBROUTINE TWIO             reads or writes contents of TWCOM.INC from or to .sol file
C     SUBROUTINE TWOUT            writes contents of TWCOM formatted to ILUOUT for debugging purposes
C
C ---------------------------------------------------------------------------------
C
      subroutine twextract (ilu)
C
C ---------------------------------------------------------------------------------
C
c
c     Routine writes data for Theis wells to .xtr file.
c     Routine is called in EXTRACT.FOR
c
      implicit none
      include 'twcom.inc'
      include 'lusys.inc'
      COMPLEX(8) cz
      REAL(8) rh,rfhead,rhead
      integer(4) ilu,i
      LOGICAL lnew
c
      if (ntw.gt.0) then       ! write out Theis wells
       write (ilu,1000)
       lnew=.true.
       do i=1,ntw
       cz=ctwz(i)
       rhead=rfhead(cz)
      if (lnew) write(ilu,2000) ctwz(i),rtwr0(i),rtwq(i),rhead,atwlab(i)
       lnew=atwlab(i+1).ne.atwlab(i)  ! test for multiple occurances due to cyclic pumping.
       END do
      endif
      return
 1000 format ('! Theis wells',/,
     & '*     x1            y1          radius         Q',
     & '             head            label')
 2000 format (5(e14.7,','),a16)
      end
C
C ---------------------------------------------------------------------------------
C
      SUBROUTINE TWDATIO (ILU)
C
C ---------------------------------------------------------------------------------
C
C
C     Routine writes all transient well data to an input (.dat) file
C     Routine is called by DATIO
C
      IMPLICIT NONE
      INTEGER(4) ILU,I
      INCLUDE 'TWCOM.INC'
      INCLUDE 'LUSYS.INC'
C
      IF (NTW.EQ.0) RETURN
      WRITE (ILU,1000)
      DO 10 I=1,NTW
      WRITE (ILU,2000) CTWZ(I),RTWQ(I),RTWST(I),RTWH0(I),RTWSTO(I),
     &                 RTWR0(I),ATWLAB(I)
  10  CONTINUE
      WRITE (ILU,3000)
      RETURN
C
 1000 FORMAT (' twell',/,' discharge',/,
     &        '*   x         y     discharge   t0    h0        Ss     '
     &       ,'radius  label')
 2000 FORMAT (3(F9.1,1X),F5.1,F9.1,G11.4,F5.1,1X,A16)
 3000 FORMAT (' quit')
      END
C
C ---------------------------------------------------------------------------------
C
      SUBROUTINE TWIO (ICODE,ILU,RVERSION,ierr)
C
C ---------------------------------------------------------------------------------
C
C    
C     Routine reads or writes contents of TWCOM.INC to an
C     external file.
C     ICODE=41 write
C     ICODE=73 read
C
      IMPLICIT NONE
      INTEGER(4) ICODE,ILU,IERR
      REAL(8) RVERSION
      INCLUDE 'TWCOM.INC'
c
      if (ierr.ne.0) return
      CALL BUFIO4 (NTW,1,ILU,ICODE,IERR)
      IF (NTW.EQ.0) RETURN
      CALL BUFIOR (RTWCT,3,ILU,ICODE,IERR)
      CALL BUFIOC (CTWZ,NTW,ILU,ICODE,IERR)  
      CALL BUFIOR (RTWR0,NTW,ILU,ICODE,IERR)          
      CALL BUFIOR (RTWST,NTW,ILU,ICODE,IERR)
      CALL BUFIOR (RTWQ,NTW,ILU,ICODE,IERR)
      CALL BUFIOR (RTWH0,NTW,ILU,ICODE,IERR)
      CALL BUFIOR (RTWSTO,NTW,ILU,ICODE,IERR)
      CALL BUFIOA (ATWLAB,NTW,ILU,ICODE,IERR)
      IF (RVERSION.EQ.1.0) RETURN
      IF (RVERSION.EQ.2.0) RETURN
      RETURN
      END
C
C ---------------------------------------------------------------------------------
C
      SUBROUTINE TWOUT (ICALL)
C
C ---------------------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER(4) ICALL,I
      CHARACTER(1) ADUM
      INCLUDE  'TWCOM.INC'
      INCLUDE  'LUSYS.INC'
      WRITE (ILUOUT,1002) NTW,RTWCT,RTWRD0,RTWST0,RPI4
      IF (NTW.EQ.0) RETURN
      DO 10 I=1,NTW
      WRITE (ILUOUT,1006) I,CTWZ(I),RTWQ(I),RTWST(I),RTWH0(I),
     &                     ATWLAB(I)
  10  CONTINUE
      WRITE (ILUOUT,2000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000) ICALL
      READ (ILUIN,3000) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE (ILUOUT,1007)
      DO 20 I=1,NTW
      WRITE (ILUOUT,1008) I,RTWR0(I),RTWSTO(I),ATWLAB(I)
  20  CONTINUE
      WRITE (ILUOUT,2000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000) ICALL
      READ (ILUIN,3000) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
 1002 FORMAT (' TW: NTW ',I3,/,
     & ' RTWCT,RTWRD0,RTWST0,RPI4 ',4(1X,G11.4),/
     & ' no.      x           y            Q           t          h0   '
     & ,'    label')
 1006 FORMAT (' ',I3,1X,2G11.4,3(1X,G11.4),1X,A16)
 1007 FORMAT (' no.      R0         S       label')
 1008 FORMAT (' ',I3,2(1X,G11.4),1X,A16)
 2000 FORMAT (' TW: ICALL= ',I3,' press <Enter> to continue.')
 3000 FORMAT (A1)
      RETURN
      END


