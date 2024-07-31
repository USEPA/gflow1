C     Last change:  HMH   8 Aug 2005    6:47 pm
C     This file contains the following routines and functions
C
C     SUBROUTINE wlextracthead    writes data for head specified wells to .xtr file
C          ENTRY wlextractdisch   writes data for discharge specified wells to .xtr file
C     SUBROUTINE WLDATIO          rites all well data to a ".dat" file
C     SUBROUTINE WLIO             reads or writes contents of WLCOM.INC from or to .sol file
C     SUBROUTINE WLOUT            writes contents of WLCOM formatted to ILUOUT for debugging purposes
C
C
C -------------------------------------------------------------------------
C
      subroutine wlextracthead (ilu)
C
C -------------------------------------------------------------------------
C
c
c     Routine writes data for HEAD SPECIFIED wells to .xtr file.
c     Routine is called in EXTRACT.FOR
c
      implicit none
      INCLUDE 'wlcom.inc'
      INCLUDE 'lusys.inc'
      COMPLEX(8) cz
      REAL(8) rh,rfhead
      INTEGER(4) ilu,i
      logical lwell,lw3head,lw3discharge
c
      lwell=.false.
      if (nwl.gt.0) then       ! look if there are head specified 2D wells
        do 5 i=1,nwl
        if (lwlh(i)) lwell=.true.
  5     continue
      endif
      if (.not.lwell) return  ! there are no head specified wells
      write (ilu,2000)
      do 10 i=1,nwl
      if (lwlh(i)) write(ilu,3000) cwlz(i),rwlr(i),rwlq(i),rwlh(i),
     & rwlerrh(i),awlab(i)
  10  continue
      return
C ------------------------------
      entry wlextractdisch (ilu)
C ------------------------------
c
c     Routine writes data for DISCHARGE SPECIFIED wells to .xtr file.
c     Routine is called in EXTRACT.FOR
c
      lwell=.false.
      if (nwl.gt.0) then       ! look if there are discharge specified 2D wells
        do 15 i=1,nwl    
        if (.not.lwlh(i)) lwell=.true.
  15    continue
      endif
      if (.not.lwell) return  ! there are no discharge specified wells
      write (ilu,5000)
      do 20 i=1,nwl
      if (.not.lwlh(i)) then
      cz=cwlz(i)                ! provide calculated head at discharge specified well
      if (rwlh(i).EQ.-9999.0) then
      rh=-9999.0
      else
      rh=rfhead(cz)
      end if
      write(ilu,3050) cwlz(i),rwlr(i),rwlq(i),rh,awlab(i)
      endif
  20  continue
      return
 2000 format ('! head specified wells',/,
     & '*     x1            y1          radius         Q',
     & '             head            % error in head   label')
 3000 format (6(e14.7,','),a16)
 3050 format (5(e14.7,','),a16)
 5000 format ('! discharge specified wells',/,
     & '*     x1            y1          radius         Q',
     & '             head            label')
      end
C
C -------------------------------------------------------------------------
C
      SUBROUTINE WLDATIO (ILU)
C
C -------------------------------------------------------------------------
C
C
C     Routine writes all well data to a ".dat" file.
C     Routine is called in DATIO.
C
      IMPLICIT NONE
      INTEGER(4) ILU,I
      LOGICAL LWH,LWS
      INCLUDE 'wlcom.inc'
      INCLUDE 'lusys.inc'
C
      IF (NWL.EQ.0) RETURN
      WRITE (ILU,1000)
      DO 15 I=1,NWL     ! check for presence of head and disch. spec. wells
      IF (LWLH(I)) LWH=.TRUE.
      IF (.NOT.LWLH(I)) LWS=.TRUE.
  15  CONTINUE
      IF (LWH) THEN      ! write head specified wells
        WRITE (ILU,2000)
        DO 22 I=1,NWL
        IF (LWLH(I)) 
     &  WRITE (ILU,3000) CWLZ(I),RWLH(I),RWLR(I),AWLAB(I)
   22   CONTINUE
      ENDIF
      IF (LWS) THEN      ! write discharge specified wells
        WRITE (ILU,4000)
        DO 24 I=1,NWL
        IF (.NOT.LWLH(I))
     &  WRITE (ILU,3000) CWLZ(I),RWLQ(I),RWLR(I),AWLAB(I)
   24   CONTINUE
      ENDIF
      WRITE (ILU,5000)
      RETURN
C      
 1000 FORMAT (' well')
 2000 FORMAT (' head',/,'*      x             y             head      ',
     &        ' radius      label')     
 3000 FORMAT (3(G14.7,1X),G11.4,1X,A16)
 4000 FORMAT (' discharge',/,'*      x             y         discharge',
     &        '      radius      label')     
 5000 FORMAT (' quit')
      END
C
C -------------------------------------------------------------------------
C
      SUBROUTINE WLIO (ICODE,ILU,RVERSION,ierr)
C
C -------------------------------------------------------------------------
C
C
C     Routine reads or writes contents of WLCOM.INC common blocks to
C     an external file.
C     ICODE=41 write
C     ICODE=73 read
C
      IMPLICIT NONE
      INTEGER(4) ICODE,ILU,IERR
      REAL(8) RVERSION
      INCLUDE 'wlcom.inc'
c
      if (ierr.ne.0) return
      CALL BUFIO4 (NWL,1,ILU,ICODE,IERR)
      IF (NWL.EQ.0) RETURN
      CALL BUFIOR (RPI2,2,ILU,ICODE,IERR)
      CALL BUFIOC (CWLZ,NWL,ILU,ICODE,IERR)
      CALL BUFIOR (RWLQ,NWL,ILU,ICODE,IERR)
      CALL BUFIOR (RWLH,NWL,ILU,ICODE,IERR)
      CALL BUFIOR (RWLR,NWL,ILU,ICODE,IERR)
      CALL BUFIOA (AWLAB,NWL,ILU,ICODE,IERR)
      CALL BUFIOL (LWLH,NWL,ILU,ICODE,IERR)
      IF (RVERSION.LE.2.0) RETURN
      CALL BUFIOR (RWLERRH,NWL,ILU,ICODE,IERR)
      IF (RVERSION.LE.10.0) RETURN
      CALL BUFIOR (DWLPOT,NWL,ILU,ICODE,IERR)
      RETURN
      END
C
C -------------------------------------------------------------------------
C
      SUBROUTINE WLOUT (ICALL)
C
C -------------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER(4) ICALL,I,J
      INCLUDE 'wlcom.inc'
      INCLUDE 'lusys.inc'
      CHARACTER(1) ADUM
      SAVE
      WRITE (ILUOUT,1002) NWL,RPI2
      WRITE (ILUOUT,1004) (I,CWLZ(I),AWLAB(I),I=1,NWL)
      WRITE (ILUOUT,2000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000)
      READ (ILUIN,3000) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE (ILUOUT,1006) (I,RWLQ(I),I=1,NWL)
      WRITE (ILUOUT,2000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000)
      READ (ILUIN,3000) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE (ILUOUT,1008) (I,RWLH(I),I=1,NWL)
      WRITE (ILUOUT,2000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000)
      READ (ILUIN,3000) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE (ILUOUT,1010) (I,RWLR(I),I=1,NWL)
      WRITE (ILUOUT,2000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000)
      READ (ILUIN,3000) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE (ILUOUT,1012) (I,LWLH(I),I=1,NWL)
      WRITE (ILUOUT,2000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000)
      READ (ILUIN,3000) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE (ILUOUT,1014) (I,(RWLD0(J,I),J=1,3),I=1,NWL)
      WRITE (ILUOUT,2000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000)
      READ (ILUIN,3000) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
 1002 FORMAT (' WLOUT:NWL,RPI2 ',I3,G11.4)
 1004 FORMAT (20(' WLOUT:I,CWLZ(I) ',I3,G11.4,G11.4,2X,A16/))
 1006 FORMAT (10(' WLOUT:I,RWLQ(I) ',4(I3,G11.4)/))
 1008 FORMAT (10(' WLOUT:I,RWLH(I) ',4(I3,G11.4)/))
 1010 FORMAT (10(' WLOUT:I,RWLR(I) ',4(I3,G11.4)/))
 1012 FORMAT (10(' WLOUT:I,LWLH(I) ',4(I3,L3)/))
 1014 FORMAT (120(' WLOUT:I,RWLD0(1-3,I) ',(I3,G11.4,G11.4,G11.4)/))
 2000 FORMAT (' WLOUT: ICALL= ',I3,' press <Enter> to continue.')
 3000 FORMAT (A1)
      RETURN
      END
