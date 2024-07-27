C     Last change:  HMH   8 Aug 2005    6:47 pm
c     This file contains the following routines and functions
c
c     SUBROUTINE DIIO         reads or writes common block data to .SOL file
c     SUBROUTINE DIDATIO      writes 3D disc data to .DAT file
c     SUBROUTINE DIOUT        writes common block data (formatted) to disk for debugging purposes
c     subroutine diextracthead  writes head specified 3D disc data to .XTR file
c         entry diextractdisch  writes discharge specified 3D disc data to .XTR file
c
c
c
c ----------------------------------------------------------------------
c
      SUBROUTINE DIIO (ICODE,ILU,RVERSION,ierr)
c
c ----------------------------------------------------------------------
c
C
C     Routine reads or writes contents of DICOMN.INC common blocks
C     to an external file.
C     ICODE=41 write
C     ICODE=73 read
C
      IMPLICIT NONE
      INTEGER(4) ICODE,ILU,IERR,NWORD
      REAL(8) RVERSION
      INCLUDE 'dicom.inc'
C
      if (ierr.ne.0) return
      CALL BUFIO4 (NDIS,1,ILU,ICODE,IERR)
      IF (NDIS.EQ.0) RETURN
      NWORD=3*NDIS
      CALL BUFIOR (RDIZ,NWORD,ILU,ICODE,IERR)
      CALL BUFIOR (RDIR,NDIS,ILU,ICODE,IERR)
      CALL BUFIOR (RDIS,NDIS,ILU,ICODE,IERR)
      CALL BUFIOR (RDIH,NDIS,ILU,ICODE,IERR)
      CALL BUFIOR (RDICPT,NDIS,ILU,ICODE,IERR)
      CALL BUFIOR (RDIA,NDIS,ILU,ICODE,IERR)
      CALL BUFIOL (LDIH,NDIS,ILU,ICODE,IERR)
      CALL BUFIOA (ADILAB,NDIS,ILU,ICODE,IERR)
      CALL BUFIO4 (NDIMAG,2,ILU,ICODE,IERR)
      IF (RVERSION.EQ.1.0) RETURN
      IF (RVERSION.EQ.2.0) RETURN
      CALL BUFIOR (RDIERR,NDIS,ILU,ICODE,IERR)
C
      RETURN
      END
c
c ----------------------------------------------------------------------
c
      SUBROUTINE DIDATIO (ILU)
c
c ----------------------------------------------------------------------
c
C
C     Routine writes all 3D sinkdisc data to a ".dat" file.
C     Routine is called in DATIO.
C
      IMPLICIT NONE
      INTEGER(4) ILU,I
      LOGICAL LHEAD,LDISCHARGE
      REAL(8) RZ
      COMPLEX(8) CZ
      INCLUDE 'dicom.inc'
      INCLUDE 'lusys.inc'
C
      IF (NDIS.EQ.0) RETURN
      WRITE (ILU,1000)
      LHEAD=.FALSE.
      LDISCHARGE=.FALSE.
      DO 10 I=1,NDIS
      IF (.NOT.LDIH(I)) LDISCHARGE=.TRUE.
      IF (LDIH(I)) LHEAD=.TRUE.
  10  CONTINUE      
      IF (LHEAD) THEN
        WRITE (ILU,2000)
        DO 20 I=1,NDIS
        IF (.NOT.LDIH(I)) GOTO 20
        CZ=CMPLX(RDIZ(1,I),RDIZ(2,I))
        RZ=RDIZ(3,I)
        WRITE (ILU,3000) CZ,RZ,RDIR(I),RDIH(I),RDICPT(I),ADILAB(I)
  20    CONTINUE
      ENDIF     
      IF (LDISCHARGE) THEN
        WRITE (ILU,4000)
        DO 30 I=1,NDIS
        IF (LDIH(I)) GOTO 30
        CZ=CMPLX(RDIZ(1,I),RDIZ(2,I))
        RZ=RDIZ(3,I)
        WRITE (ILU,5000) CZ,RZ,RDIR(I),RDIS(I),ADILAB(I)
  30    CONTINUE      
      ENDIF
      WRITE (ILU,6000)
      RETURN
C      
 1000 FORMAT (' sd3d')
 2000 FORMAT (' head',/,'*    x       y          z         radius  ',
     &        '   head     off set  label')
 3000 FORMAT (4(F9.0,1X),G14.7,1X,F8.0,1X,A16)
 4000 FORMAT (' discharge',/,'*    x       y         z        radius',
     &        '  exfiltr.rate    label')
 5000 FORMAT (4(F9.0,1X),G14.7,1X,A16)     
 6000 FORMAT (' quit')
      END
c
c ----------------------------------------------------------------------
c
      SUBROUTINE DIOUT (ICALL)
c
c ----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER(4) ICALL,I,J
      CHARACTER(1) ADUM
      INCLUDE 'dicom.inc'
      INCLUDE 'lusys.inc'
      SAVE
      WRITE (ILUOUT,1002) (I,(RDIZ(J,I),J=1,3),I=1,NDIS)
      WRITE (ILUOUT,2000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000)
      READ (ILUIN,3000) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE (ILUOUT,1004) (I,RDIR(I),I=1,NDIS)
      WRITE (ILUOUT,1006) (I,RDIS(I),I=1,NDIS)
      WRITE (ILUOUT,2000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000)
      READ (ILUIN,3000) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN      
      WRITE (ILUOUT,1008) (I,RDIH(I),I=1,NDIS)
      WRITE (ILUOUT,1010) (I,RDICPT(I),I=1,NDIS)
      WRITE (ILUOUT,1015) (I,LDIH(I),I=1,NDIS)
      WRITE (ILUOUT,2000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000)
      READ (ILUIN,3000) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN      
      WRITE (ILUOUT,1012) (I,RDIA(I),I=1,NDIS)
      WRITE (ILUOUT,1013) (I,ADILAB(I),I=1,NDIS)
      WRITE (ILUOUT,1014) NDIS,NDIMAG,NDIRH
      WRITE (ILUOUT,1016) RDISR,RDISR2,RTMR,RTPR
      WRITE (ILUOUT,1018) RTDN,ROOT,RTARA,RK2
      WRITE (ILUOUT,1020) REPS,IEL
      WRITE (ILUOUT,2000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000)
      READ (ILUIN,3000) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN      
 1002 FORMAT (120(' I,RDIZ(1-3,I) ',(I3,E12.4,E12.4,E12.4)/))
 1004 FORMAT (10(' I,RDIR(I) ',4(I3,G11.4)/))
 1006 FORMAT (10(' I,RDIS(I) ',4(I3,G11.4)/))
 1008 FORMAT (10(' I,RDIH(I) ',4(I3,G11.4)/))
 1010 FORMAT (10(' I,RDICPT(I) ',4(I3,G11.4)/))
 1012 FORMAT (10(' I,RDIA(I) ',4(I3,G11.4)/))
 1013 FORMAT (10(' I,ADILAB(I) ',3(I3,1X,A16)/)) 
 1014 FORMAT (' NDIS,NDIMAG,NDIRH ',3I5)
 1015 FORMAT (10(' I,LDIH(I) ',4(I3,L5)/))
 1016 FORMAT (' RDISR,RDISR2,RTMR,RTPR ',4G11.4)
 1018 FORMAT (' RTDN,ROOT,RTARA,RK2 ',4G11.4)
 1020 FORMAT (' REPS,IEL ',G11.4,I5)
 2000 FORMAT (' DIOUT: ICALL= ',I3,' press <Enter> to continue.')
 3000 FORMAT (A1)
      RETURN
      END
C
C -------------------------------------------------------------------------
C
      subroutine diextracthead (ilu)
C
C -------------------------------------------------------------------------
C
c
c     Routine writes data for HEAD SPECIFIED sink discs to .xtr file.
c     Routine is called in EXTRACT.FOR
c
      implicit none
      integer(4) ilu,i
      logical ldisc
      REAL(8) rh,rfhead
      COMPLEX(8) cz
      INCLUDE 'dicom.inc'
      INCLUDE 'lusys.inc'
c
      ldisc=.false.
      if (ndis.gt.0) then       ! look if there are head specified sink discs
        do 5 i=1,ndis
        if (ldih(i)) ldisc=.true.
  5     continue 
      endif 
      if (.not.ldisc) return  ! there are no head specified sink discs
      write (ilu,2000)
      do 10 i=1,ndis
      if (ldih(i)) write(ilu,3000) rdiz(1,i),rdiz(2,i),rdiz(3,i),
     &             rdir(i),rdis(i),rdih(i),rdicpt(i),rdierr(i),adilab(i)
  10  continue
      return
C ------------------------------
      entry diextractdisch (ilu)
C ------------------------------
c
c     Routine writes data for DISCHARGE SPECIFIED sink discs to .xtr file.
c     Routine is called in EXTRACT.FOR
c
      ldisc=.false.
      if (ndis.gt.0) then       ! look if there are discharge specified sink discs
        do 15 i=1,ndis
        if (.not.ldih(i)) ldisc=.true.
  15    continue 
      endif 
      if (.not.ldisc) return ! there are no discharge specified sink discs
      write (ilu,5000)
      do 20 i=1,ndis
      if (.not.ldih(i)) then
      cz=CMPLX(rdiz(1,i),rdiz(2,i))  ! provide calculated head at discharge specified sink disc
      rh=rfhead(cz)            
      write(ilu,3050) rdiz(1,i),rdiz(2,i),rdiz(3,i),
     &             rdir(i),rdis(i),rh,adilab(i)
      endif
  20  continue
      return
 2000 format ('! head specified sink discs (3d)',/,
     & '*     x1            y1          z1             radius   ',
     & '     sink density    head       dist. to coll. pnt    ',
     & '% error in head     label')
 3000 format (8(e14.7,','),2x,a16)
 3050 format (6(e14.7,','),2x,a16)
 5000 format ('! discharge specified sink discs (3d)',/,
     & '*     x1            y1            z1            radius   ',
     & '     sink density  calc. head           label')
      end


      