C     Last change:  HMH   9 Aug 2005    7:27 pm
c
c
c     SUBROUTINE BUFIN      initializes buffer for binary IO
c     SUBROUTINE BUFEX      flushes buffer at the end of binary IO
c     SUBROUTINE BUFIO2     adds integer*2 words to buffer for binary IO
c     SUBROUTINE BUFIO4     adds integer*4 words to buffer for binary IO
c     SUBROUTINE BUFIOR     adds real*4 words to buffer for binary IO (currently adds real*8)
c     SUBROUTINE BUFIOA     adds character words to buffer for binary IO
c     SUBROUTINE BUFIOC     adds complex*4 words to buffer for binary IO (currently complex*8)
c     SUBROUTINE BUFIOD     adds real*8 words to buffer for binary IO
c     SUBROUTINE BUFIOCD    adds complex*8 words to buffer for binary IO
c     SUBROUTINE BUFIOL     adds logical*4 words to buffer for binary IO
c     SUBROUTINE BUFIO      adds words to buffer for binary IO (called by other BUFIO routines)
c     SUBROUTINE XFERW2     routine swaps bytes of data between two arrays (called by BUFIO)
c     SUBROUTINE SYSIOW     performs binary read or write operation

C
C --------------------------------------------------------------------------
C
      SUBROUTINE BUFIN (ICODE,ILU,IERR)
C
C --------------------------------------------------------------------------
C
C
C     First call BUFIN
C     Call BUFIO2/4/R/C/L for IO of buffers of the type
C     integer(2)/(4)/real/complex/logical
C     Last call BUFEX
C
      IMPLICIT NONE
      INTEGER(4) ICODE,ILU,IERR,IBUF,IBYTES,NBYTES,NBTSUM,NBMAX
      COMMON /BUFCOM/ IBUF(256),IBYTES,NBYTES,NBTSUM
      INCLUDE 'LUSYS.INC'
C
      DATA NBMAX /1024/
C
      NBTSUM=0  ! total number of bytes read or written
      NBYTES=NBMAX  ! current number of bytes read or written
      IBYTES=0  ! current number of bytes in buffer IBUF
      IERR=0    ! error code
      IF(ICODE.EQ.73) CALL SYSIOW (ICODE,ILU,IBUF,NBYTES,NBMAX,IERR)
      RETURN
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE BUFEX (ICODE,ILU,IERR)
C
C --------------------------------------------------------------------------
C
C
c     last call BUFEX
c
      IMPLICIT NONE
      INTEGER(4) ICODE,ILU,IERR,IBUF,IBYTES,NBYTES,NBTSUM
      COMMON /BUFCOM/ IBUF(256),IBYTES,NBYTES,NBTSUM
      INCLUDE 'LUSYS.INC'
C
      IF (IERR.NE.0) RETURN
       NBTSUM=NBTSUM+IBYTES
      IF (ICODE.EQ.41) GOTO 10
       WRITE (ILUME,1000) NBTSUM
      RETURN
  10  CALL SYSIOW(ICODE,ILU,IBUF,IBYTES,NBYTES,IERR)
      ENDFILE(ILU,IOSTAT=IERR)
      WRITE (ILUME,2000) NBTSUM
      RETURN
 1000 FORMAT (' Total number of bytes read: ',I7/)
 2000 FORMAT (' Total number of bytes written: ',I7/)
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE BUFIO2 (IARR2,NWORD,ILU,ICODE,IERR)
C
C --------------------------------------------------------------------------
C
C
      IMPLICIT NONE
      INTEGER(2) IARR2,IBUF2
      INTEGER(4) NWORD,ILU,IERR,NBYTES,NBTSUM,ICODE,IBYTES,NBYT
      COMMON /BUFCOM/ IBUF2(512),IBYTES,NBYTES,NBTSUM
      INCLUDE 'LUSYS.INC'
C
      DIMENSION IARR2(*)
C
      IF (IERR.NE.0) RETURN
      NBYT=2*NWORD
      CALL BUFIO(IARR2,NBYT,ILU,ICODE,IERR)
      IF(IERR.EQ.100) WRITE(ILUER,1000) IBYTES
      RETURN
 1000 FORMAT (' ***ERROR: IBYTES out of range in BUFIO2: IBYTES=',I5/)
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE BUFIO4 (IARR4,NWORD,ILU,ICODE,IERR)
C
C --------------------------------------------------------------------------
C
C
      IMPLICIT NONE
      INTEGER(4) IARR4,NWORD,ILU,ICODE,IERR,
     &           IBUF4,IBYTES,NBYTES,NBTSUM,
     &           NBYT
      COMMON /BUFCOM/ IBUF4(256),IBYTES,NBYTES,NBTSUM
      INCLUDE 'LUSYS.INC'
C
      NBYT=4*NWORD
      CALL BUFIO (IARR4,NBYT,ILU,ICODE,IERR)
      IF (IERR.EQ.100) WRITE (ILUER,1000) IBYTES
      RETURN
 1000 FORMAT (' ***ERROR: IBYTES out of range in BUFIO4: IBYTES=',I5/)
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE BUFIOR (RARR,NWORD,ILU,ICODE,IERR)
C
C --------------------------------------------------------------------------
C
c
      IMPLICIT NONE
      INTEGER(4) NWORD,ILU,ICODE,IERR
      REAL(8) RARR
      call bufiod (RARR,NWORD,ILU,ICODE,IERR) !  fix to implement double precision
C
c      COMMON /BUFCOM/ RBUF(256),IBYTES,NBYTES,NBTSUM
c      INCLUDE 'LUSYS.INC'
C
c      NBYT=4*NWORD
c      CALL BUFIO (RARR,NBYT,ILU,ICODE,IERR)
c      IF (IERR.EQ.100) WRITE (ILUER,1000) IBYTES
c      RETURN
c 1000 FORMAT (' ***ERROR: IBYTES out of range in BUFIOR: IBYTES=',I5/)
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE BUFIOA (AARR,NWORD,ILU,ICODE,IERR)
C
C --------------------------------------------------------------------------
C
C
      IMPLICIT NONE
      INTEGER(2) I2BUF(16)
      INTEGER(4) NWORD,ILU,ICODE,IERR,I,J,
     &           IBUF4,IBYTES,NBYTES,NBTSUM,IDUM
      CHARACTER(16) AARR(*),ADUM
      COMMON /BUFCOM/ IBUF4(256),IBYTES,NBYTES,NBTSUM
      INCLUDE 'LUSYS.INC'
C
      IF (IERR.NE.0) RETURN
      DO 30 I=1,NWORD
      IF (ICODE.EQ.41) THEN    ! convert to INTEGER(2) and write
      ADUM=AARR(I)
      DO 10 J=1,16
      I2BUF(J)=ICHAR(ADUM(J:J))
  10  CONTINUE
      CALL BUFIO (I2BUF,32,ILU,ICODE,IERR)
      IF (IERR.EQ.100) THEN
      WRITE (ILUER,1000) IBYTES
      RETURN
      ENDIF
      ELSE                     ! read and convert to CHARACTER(16)
      CALL BUFIO (I2BUF,32,ILU,ICODE,IERR)
      IF (IERR.EQ.100) THEN
      WRITE (ILUER,1000) IBYTES
      RETURN
      ENDIF
      DO 20 J=1,16
      IDUM=I2BUF(J)
      ADUM(J:J)=CHAR(IDUM)
  20  CONTINUE
      AARR(I)=ADUM
      ENDIF
  30  CONTINUE      
      RETURN            
 1000 FORMAT (' ***ERROR: IBYTES out of range in BUFIOA: IBYTES=',I5/)
      END      
C
C --------------------------------------------------------------------------
C
      SUBROUTINE BUFIOC (CARR,NWORD,ILU,ICODE,IERR)
C
C --------------------------------------------------------------------------
C
c
      IMPLICIT NONE
      INTEGER(4) NWORD,ILU,ICODE,IERR
      COMPLEX(8) CARR
      call buficd (CARR,NWORD,ILU,ICODE,IERR) !  fix to implement double precision
C
c      IMPLICIT COMPLEX (C)
c      COMPLEX CARR
c      COMMON /BUFCOM/ CBUF(128),IBYTES,NBYTES,NBTSUM
c      INCLUDE 'LUSYS.INC'
C
c      IF (IERR.NE.0) RETURN
c      NBYT=8*NWORD
c      CALL BUFIO(CARR,NBYT,ILU,ICODE,IERR)
c      IF(IERR.EQ.100) WRITE(ILUER,1000) IBYTES
c      RETURN
c 1000 FORMAT (' ***ERROR: IBYTES out of range in BUFIOC: IBYTES=',I5/)
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE BUFIOD (DARR,NWORD,ILU,ICODE,IERR)
C
C --------------------------------------------------------------------------
C
C
      IMPLICIT NONE
      INTEGER(4) NWORD,ILU,ICODE,IERR,
     &           IBYTES,NBYTES,NBTSUM,NBYT
      REAL(8) DARR,DBUF
      COMMON /BUFCOM/ DBUF(128),IBYTES,NBYTES,NBTSUM
      INCLUDE 'LUSYS.INC'
C
      IF (IERR.NE.0) RETURN
      NBYT=8*NWORD
      CALL BUFIO(DARR,NBYT,ILU,ICODE,IERR)
      IF(IERR.EQ.100) WRITE(ILUER,1000) IBYTES
      RETURN
 1000 FORMAT (' ***ERROR: IBYTES out of range in BUFIOD: IBYTES=',I5/)
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE BUFICD (CDARR,NWORD,ILU,ICODE,IERR)
C
C --------------------------------------------------------------------------
C
C
      IMPLICIT NONE
      INTEGER(4) NWORD,ILU,ICODE,IERR,
     &           IBYTES,NBYTES,NBTSUM,NBYT
      COMPLEX(8) CDARR,CDBUF
      COMMON /BUFCOM/ CDBUF(64),IBYTES,NBYTES,NBTSUM
      INCLUDE 'LUSYS.INC'
C
      IF (IERR.NE.0) RETURN
      NBYT=16*NWORD
      CALL BUFIO(CDARR,NBYT,ILU,ICODE,IERR)
      IF(IERR.EQ.100) WRITE(ILUER,1000) IBYTES
      RETURN
 1000 FORMAT (' ***ERROR: IBYTES out of range in BUFICD: IBYTES=',I5/)
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE BUFIOL (LARR,NWORD,ILU,ICODE,IERR)
C
C --------------------------------------------------------------------------
C
C
      IMPLICIT NONE
      INTEGER(4) NWORD,ILU,ICODE,IERR,
     &           IBYTES,NBYTES,NBTSUM,NBYT
      LOGICAL LARR,LBUF
      COMMON /BUFCOM/ LBUF(256),IBYTES,NBYTES,NBTSUM
      INCLUDE 'LUSYS.INC'
C
      NBYT=4*NWORD
      CALL BUFIO (LARR,NBYT,ILU,ICODE,IERR)
      IF (IERR.EQ.100) WRITE (ILUER,1000) IBYTES
      RETURN
 1000 FORMAT (' ***ERROR: IBYTES out of range in BUFIOL: IBYTES=',I5/)
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE BUFIO (IARR,NBYT,ILU,ICODE,IERR)
C
C --------------------------------------------------------------------------
C
C
      IMPLICIT NONE
      INTEGER(4) IARR,NBYT,ILU,ICODE,IERR,
     &           IBUF4,IBYTES,NBYTES,NBTSUM,
     &           NBYTL,NSWAP,IARST,IST
      DIMENSION IARR(*)
      COMMON /BUFCOM/ IBUF4(256),IBYTES,NBYTES,NBTSUM
      INCLUDE 'LUSYS.INC'
C
      IF (IERR.NE.0) RETURN
      NBYTL=NBYT
      NSWAP=0
      IARST=1
  10  IF (NBYTL.EQ.0) RETURN
      IST=IBYTES+1
      IARST=IARST+NSWAP
      IF (IST.GT.NBYTES) GOTO 20
      NSWAP=NBYTES-IST+1
      IF (NSWAP.GE.NBYTL) NSWAP=NBYTL
      IST=IST-1
      IARST=IARST-1
      IF (ICODE.EQ.73) CALL XFERW2 (NSWAP,IST,IBUF4,IARST,IARR)
      IF (ICODE.EQ.41) CALL XFERW2 (NSWAP,IARST,IARR,IST,IBUF4)
      IST=IST+1
      IARST=IARST+1
      IBYTES=IBYTES+NSWAP
      IF (IBYTES.LT.NBYTES) RETURN
      CALL SYSIOW (ICODE,ILU,IBUF4,NBYTES,NBYTES,IERR)
      NBTSUM=NBTSUM+NBYTES
      IF (IERR.NE.0) RETURN
      IBYTES=0
      NBYTL=NBYTL-NSWAP
      GOTO 10
  20  IERR=100
      RETURN
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE XFERW2(NBYTES,ISTIN,I2BFIN,ISTOUT,I2BFOU)
C
C --------------------------------------------------------------------------
C
C
      IMPLICIT NONE
      INTEGER(2) I2BFIN(*),I2BFOU(*)
      INTEGER(4) NBYTES,ISTIN,ISTOUT,NWORDS,ISTWI,JOUT,JIN
C
      NWORDS=NBYTES/2
      ISTWI=ISTIN/2+1
      JOUT=ISTOUT/2
      DO 100 JIN=ISTWI,ISTWI+NWORDS-1
      JOUT=JOUT+1
      I2BFOU(JOUT)=I2BFIN(JIN)
 100  CONTINUE
      RETURN
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE SYSIOW(ICODE,ILU,I2BUF,NBYTES,NBMAX,IERR)
C
C --------------------------------------------------------------------------
C
C
      IMPLICIT NONE
      INTEGER(2) I2BUF(*)
      INTEGER(4) ICODE,ILU,NBYTES,NBMAX,IERR,
     &           NWORDS,NWMAX,I
      INCLUDE 'LUSYS.INC'
C
      NWORDS=NBYTES/2
      NWMAX=NBMAX/2
      IF(ICODE.EQ.73) GOTO 2000
C     write
      IF(NWORDS.LT.NWMAX) THEN
        DO 50 I=NWORDS+1,NWMAX
          I2BUF(I)=0
 50     CONTINUE
      ENDIF
      WRITE(ILU,IOSTAT=IERR) (I2BUF(I),I=1,NWMAX)
      IF (IERR.NE.0) then
      WRITE (ILUER,3000) IERR
        AMESS(1)='No write access to solution file.'
        AMESS(2)='Make sure the directory and file are write enabled.'
        CALL HALT(2)  ! stop program execution for batch version
      endif
      RETURN
C     read
 2000 READ(ILU,IOSTAT=IERR)  (I2BUF(I),I=1,NWMAX)
      IF (IERR.NE.0) WRITE (ILUER,3000) IERR
      RETURN
 3000 FORMAT (' ***IO ERROR in SYSIOW: IERR=',I5)      
      END
