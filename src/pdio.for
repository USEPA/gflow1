C     Last change:  HMH   8 Aug 2005    6:47 pm
c      This file contains the following routines and functions
c
c      SUBROUTINE PDIO           reads and writes contents of PDCOM to .SOL file
c      SUBROUTINE PDDATIO        writes (2D) sink disc data to .DAT file
c      subroutine pdextracthead  writes head specified disc data to .XTR file
c          entry pdextractdisch  writes discharge specified disc data to .XTR file
c      SUBROUTINE PDOUT          writes common blocks (formatted) for debugging purposes
c
c
c
c
c -------------------------------------------------------------------
c
        SUBROUTINE PDIO (ICODE,ILU,RVERSION,ierr)
c
c -------------------------------------------------------------------
c
C 
C       Routine reads or writes contents of COMMON /PDCOM/ to
C       an external file.
C       ICODE=41 write
C       ICODE=73 read
C
        IMPLICIT NONE
        INTEGER(4) ICODE,ILU,IERR
        LOGICAL LDUM(4)
        REAL(8) RVERSION
        INCLUDE 'PDCOM.INC'
C
      if (ierr.ne.0) return
        CALL BUFIO4 (NPD,1,ILU,ICODE,IERR)
        IF (NPD.EQ.0) RETURN
        CALL BUFIO4 (NPDHD,2,ILU,ICODE,IERR)
        IF (ICODE.EQ.41) CALL BUFIOL (LPDHD,4,ILU,ICODE,IERR) ! to maintain
        IF (ICODE.EQ.73) CALL BUFIOL (LDUM,4,ILU,ICODE,IERR) ! plot command
        CALL BUFIOR (RPDC0,2,ILU,ICODE,IERR)
        CALL BUFIO4 (IPDHD,NPDHD,ILU,ICODE,IERR)
        CALL BUFIO4 (IPDRC,NPDRC,ILU,ICODE,IERR)
        CALL BUFIOL (LPDGIV,NPD,ILU,ICODE,IERR)
        CALL BUFIOA (APDLAB,NPD,ILU,ICODE,IERR)
        CALL BUFIOC (CPDZ,NPD,ILU,ICODE,IERR)
        CALL BUFIOR (RPDS,NPD,ILU,ICODE,IERR)
        CALL BUFIOR (RPDH,NPD,ILU,ICODE,IERR)
        CALL BUFIOR (RPDR,NPD,ILU,ICODE,IERR)
        CALL BUFIOR (RPDRES,NPD,ILU,ICODE,IERR)
        CALL BUFIOR (RPDEP,NPD,ILU,ICODE,IERR)
        CALL BUFIOR (RPDSB,NPD,ILU,ICODE,IERR)
        IF (RVERSION.EQ.1.0) RETURN
        IF (RVERSION.EQ.2.0) RETURN
        RETURN
        END
c
c -------------------------------------------------------------------
c
      SUBROUTINE PDDATIO (ILU)
c
c -------------------------------------------------------------------
c
C
C     Routine writes all 2D sinkdisc data to a ".dat" file.
C     Routine is called in DATIO.
C
      IMPLICIT NONE
      INTEGER(4) ILU,I,II
      REAL(8) RESISTANCE,RDEPTH,RSTOP
      COMPLEX(8) CZRAD
      INCLUDE 'PDCOM.INC'
      INCLUDE 'LUSYS.INC'
C
      IF (NPD.EQ.0) RETURN
      RESISTANCE=-1.0
      RDEPTH=-1.0
      WRITE (ILU,1000)
      IF (NPDHD.NE.0) THEN
        WRITE (ILU,2000)
        DO 22 II=1,NPDHD
        I=IPDHD(II)
        CZRAD=CPDZ(I)+RPDR(I)
        IF (RPDRES(I).NE.RESISTANCE) THEN
          RESISTANCE=RPDRES(I)
          WRITE (ILU,2600) RESISTANCE
        ENDIF
        IF (RPDEP(I).NE.RDEPTH) THEN
          RDEPTH=RPDEP(I)
          WRITE (ILU,2700) RDEPTH
        ENDIF
        WRITE (ILU,4000) CPDZ(I),CZRAD,RPDH(I),APDLAB(I)
   22   CONTINUE
      ENDIF
      IF (NPDRC.NE.0) THEN
        WRITE (ILU,5000)
        DO 24 II=1,NPDRC
        I=IPDRC(II)
        CZRAD=CPDZ(I)+RPDR(I)
        RSTOP=RPDS(I)-RPDSB(I)
        WRITE (ILU,6000) CPDZ(I),CZRAD,RSTOP,RPDSB(I),APDLAB(I)
   24   CONTINUE
      ENDIF
      WRITE (ILU,7000)
      RETURN
C      
 1000 FORMAT (' sinkdisc')
 2000 FORMAT (' head',/,'*      x       y        xr        yr   ',
     &        '    head        label')
 2600 FORMAT (' resistance ',G14.7)
 2700 FORMAT (' depth ',G14.7)     
 4000 FORMAT (4(F9.0,1X),G14.7,1X,A16)
 5000 FORMAT (' discharge',/,'*    x       y        xr        yr   ',
     &        '   top disch.  bottom disch. label')
 6000 FORMAT (4(F9.0,1X),2(G11.4,1X),A16)
 7000 FORMAT (' quit')
      END
C
C -------------------------------------------------------------------------
C
      subroutine pdextracthead (ilu)
C
C -------------------------------------------------------------------------
C
c
c     Routine writes data for HEAD SPECIFIED sink discs to .xtr file.
c     Routine is called in EXTRACT.FOR
c
      implicit none
      integer(4) ilu,i,ipd
      REAL(8) rh,rfhead
      COMPLEX(8) cz
      include 'pdcom.inc'
      include 'lusys.inc'
c
      if (npdhd.eq.0) return  ! there are no head specified wells
      write (ilu,2000)
      do ipd=1,npdhd
      i=ipdhd(ipd)
      write(ilu,3000)cpdz(i),rpdr(i),rpds(i),rpdh(i),rpdres(i),rpdep(i),
     &               apdlab(i)
      END do
      return
C ------------------------------
      entry pdextractdisch (ilu)
C ------------------------------
c
c     Routine writes data for DISCHARGE SPECIFIED sink discs to .xtr file.
c     Routine is called in EXTRACT.FOR
c
      if (npdrc.eq.0) return  ! there are no discharge specified sink discs
      write (ilu,5000)
      do ipd=1,npdrc
      i=ipdrc(ipd)
      cz=cpdz(i)                ! provide calculated head at discharge specified sink disc
      rh=rfhead(cz)            
      write(ilu,3050) cpdz(i),rpdr(i),rh,rpds(i),rpdsb(i),apdlab(i)
      END do                         
      return
 2000 format ('! head specified sink discs (2d)',/,
     & '*     x1            y1          radius       sink density',
     & '       head          resistance      depth      ',
     & '       label')
 3000 format (7(e14.7,','),2x,a16)
 3050 format (6(e14.7,','),2x,a16)
 5000 format ('! discharge specified sink discs (2d)',/,
     & '*     x1            y1          radius         calc. head  ',
     & '  top sink dens. bottom sink dens.     label')
      end
C
C -------------------------------------------------------------------------
C
      SUBROUTINE PDOUT (ICALL)
C
C -------------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER(4) ICALL,I
      CHARACTER(1) ADUM
      INCLUDE 'PDCOM.INC'
      INCLUDE 'LUSYS.INC'
      SAVE
      WRITE (ILUOUT,1000) NPD,NPDHD,NPDRC
      WRITE (ILUOUT,1010) LPDHD,LPDRC,LPDRCH,LPDPER
      WRITE (ILUOUT,1020) RPI2,RPDC0,RPDD0
      WRITE (ILUOUT,2000) ICALL      
      IF (LUOUTFILE) WRITE (ILUME,2000) ICALL
      READ (ILUIN,3000) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN      
      WRITE (ILUOUT,1030) (I,IPDHD(I),I=1,NPDHD)
      WRITE (ILUOUT,1040) (I,IPDRC(I),I=1,NPDRC)
      WRITE (ILUOUT,1050) (I,LPDGIV(I),I=1,NPD)
      WRITE (ILUOUT,2000) ICALL      
      IF (LUOUTFILE) WRITE (ILUME,2000) ICALL
      READ (ILUIN,3000) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN      
      WRITE (ILUOUT,1060) (I,APDLAB(I),I=1,NPD)
      WRITE (ILUOUT,2000) ICALL      
      IF (LUOUTFILE) WRITE (ILUME,2000) ICALL
      READ (ILUIN,3000) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN      
      WRITE (ILUOUT,1070) (I,CPDZ(I),RPDS(I),RPDSB(I),I=1,NPD)
      WRITE (ILUOUT,2000) ICALL      
      IF (LUOUTFILE) WRITE (ILUME,2000) ICALL
      READ (ILUIN,3000) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN      
      WRITE (ILUOUT,1080) (I,RPDR(I),RPDH(I),RPDRES(I),RPDEP(I),I=1,NPD)
      WRITE (ILUOUT,2000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000) ICALL
      READ (ILUIN,3000) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      RETURN
 1000 FORMAT (' NPD,NPDHD,NPDRC ',3I5)
 1010 FORMAT (' LPDHD,LPDRC,LPDRCH,LPDPER ',5L5)
 1020 FORMAT (' RPI2,RPDC0,RPDD0 ',3G11.4)
 1030 FORMAT ( 5(' IPDHD(',I3,')=',I3))
 1040 FORMAT ( 5(' IPDRC(',I3,')=',I3))
 1050 FORMAT ( 5(' LPDGIV(',I3,')=',L2))
 1060 FORMAT ( 2(' APDLAB(',I3,')=',A16))
 1070 FORMAT (' I,CPDZ,RPDS,RPDSB ',I3,2G11.4,1X,G11.4,1X,G11.4)
 1080 FORMAT (' I,RPDR,RPDH,RPDRES,RPDEP ',I3,1X,G11.4,1X,G11.4,
     &1X,G11.4,1X,G11.4)
 2000 FORMAT (' PDOUT: CALL #',I3)
 3000 FORMAT (A1)
      END


