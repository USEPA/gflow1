C     Last change:  HMH   8 Aug 2005    6:41 pm
c       This file contains the following routines and functions:
c
C	subroutine lsextracthead    writes data for head specified line sinks and drains to .xtr file
C	SUBROUTINE LSDATIO          writes all line sink data to a ".dat" file
C	SUBROUTINE LSIO             reads or writes contents of LSCOM.INC from or to .sol file
C	SUBROUTINE LSNKOUT          writes contents of LSCOM formatted to ILUOUT for debugging purposes
c
C
C-----------------------------------------------------------------------------------------------------------
C
      subroutine lsextracthead (ilu)
C
C-----------------------------------------------------------------------------------------------------------
C
c     Routine writes data for head specified line sinks and drains to .xtr file
c     Routine is called in EXTRACT.FOR
c
      implicit none
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
      integer(4) ilu,i,ils,ils1,ils2,istr
      LOGICAL ldrains,lgalleries,lheadonly
      REAL(8) rhc,rfhedp
      COMPLEX(8) cz
c
      if (nlsh.eq.0) return
      lheadonly=.false.
      ldrains=.false.
      lgalleries=.false.
      do i=1,nlsh      ! see what is there
      ils=klspth(i)
      if (.not.lsdrain(ils).and..not.lsgallery(ils)) lheadonly=.true. ! there are line-sinks that are neither drains not galleries
      if (lsdrain(ils)) ldrains=.true. ! there are drains
      if (lsgallery(ils)) lgalleries=.true. ! there are galleries
      end do
      if (lheadonly) then ! list line sinks that are no drains or galleries
      write (ilu,2000)
      do i=1,nlsh
      ils=klspth(i)
      if (.not.lsdrain(ils).and..not.lsgallery(ils)) then    ! list head specified line-sinks only
      cz=(clszs(ils)+clsze(ils))/2.0
      rhc=rfhedp(rlspot(ils),cz)
      write (ilu,3000) clszs(ils),clsze(ils),rlsh(ils),rhc,rlsig(ils),
     &                 rlswid(ils),rlsres(ils),rlsdep(ils),
     &                 rlsbf(ils),rlsof(ils),rlserror(ils),alslab(ils)
      end if
      END do
      end if
      if (ldrains) then ! list drains
       write (ilu,2100)
       do i=1,nlsh
       ils=klspth(i)
       if (lsdrain(ils)) then         ! include only drains
        cz=(clszs(ils)+clsze(ils))/2.0
        rhc=rfhedp(rlspot(ils),cz)
        write (ilu,3100)clszs(ils),clsze(ils),rlsh(ils),rhc,rlsig(ils),
     &                rlswid(ils),rlsres(ils),rlsdep(ils),rlserror(ils),
     &                alslab(ils)
       end if
       end do
      end if
      if (lgalleries) then ! list galleries
       write (ilu,2200)
       do i=1,nlsh
       ils=klspth(i)
       if (lsgallery(ils)) then         ! include only galleries
        cz=(clszs(ils)+clsze(ils))/2.0
        rhc=rfhedp(rlspot(ils),cz)
        write (ilu,3100)clszs(ils),clsze(ils),rlsh(ils),rhc,rlsig(ils),
     &                rlswid(ils),rlsres(ils),rlsdep(ils),rlserror(ils),
     &                alslab(ils)
       end if
       end do
      end if
      if (nlstring.gt.0) then     ! write linkages for stream networks
      write (ilu,3201)
      write (ilu,3200) nlstring
      do 20 istr=1,nlstring
      ils1=klstrend(istr)
      ils2=klsdn(ils1)
      if (ils2.eq.0) then
      write (ilu,3500) alslab(ils1)
      else
      write (ilu,3600) alslab(ils1),alslab(ils2)
      endif
  20  continue      
      endif
      return
c     --------------------------           
      entry lsextractdisch (ilu)
c     --------------------------
c
      if (nlsig.eq.0) return
      write (ilu,5000)
      do 30 i=1,nlsig
      ils=klspts(i)
      write (ilu,6000) clszs(ils),clsze(ils),rlsig(ils),alslab(ils)
  30  continue
      return
 2000 format ('! head specified line sinks',/,
     & '*     x1             y1             x2             y2       ',
     & ' spec. head     calc. head      discharge         width     ',
     & ' resistance       depth          baseflow     overlandflow  ',
     & ' % error BC       label')
 2100 format ('! drains',/,
     & '*     x1             y1             x2             y2       ',
     & ' drain elev.     calc. head      discharge         width     ',
     & ' resistance       depth          % error BC        label')
 2200 format ('! galleries',/,
     & '*     x1             y1             x2             y2       ',
     & ' gallery elev.  calc. wlevel     discharge         width     ',
     & ' resistance       depth          % error BC        label')
 3000 format (13(e14.7,','),a16)
 3100 format (11(e14.7,','),a16)
 3200 format ('* there are ',I3,' stream sections linked as follows:')
 3201 format ('! stream linkages')
 3500 format (1x,a16,', END_LINESINK')
 3600 format (1x,a16,',',a16)
 5000 format ('! discharge specified line sinks',/,
     & '*     x1             y1             x2             y2       ',
     & '    discharge       label')
 6000 format (5(e14.7,','),a16)
      end
C
C-----------------------------------------------------------------------------------------------------------
C
      SUBROUTINE LSDATIO (ILU)
C
C-----------------------------------------------------------------------------------------------------------
C
C     Routine writes all line sink data to a ".dat" file.
C     Routine is called in DATIO.
C
      IMPLICIT NONE
      REAL(8) RESISTANCE,RDEPTH
      INTEGER(4) ILU,ISTR,ISTART,IEND,IDUM,I,II,ICOUNT
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
C
      IF (NLS.EQ.0) RETURN
c
c debugging, REMOVE
      do i=1,nls
      write (iluer,1001) i,alslab(i)
 1001 format ('lsdatio: i,alslab(i) ',i4,2x,a16)
      end do
c
      RESISTANCE=-1.0
      RDEPTH=-1.0
      WRITE (ILU,1000)
      IF (NLSTRING.GT.0) THEN
        CALL LSOPENEND         ! fix openends or an extra head spec. ls
c                                is written (see line 52)
        DO 23 ISTR=1,NLSTRING
        ISTART=KLSTRING(ISTR)
        IEND=KLSTREND(ISTR)
        IF (RLSOVHW(ISTR).EQ.0.0) THEN
          WRITE (ILU,1500)
        ELSE
          WRITE (ILU,2000) RLSOVHW(ISTR)
        ENDIF
        IF (KLSDN(IEND).EQ.0) WRITE (ILU,2500)
        IF (RLSOFST(ISTR).NE.0.0) WRITE (ILU,3000) RLSOFST(ISTR)
        IF (IEND.LT.ISTART) THEN  ! ensure proper order of line sinks
          IDUM=ISTART
          ISTART=IEND
          IEND=IDUM
        ENDIF
        WRITE (ILU,3050)
        DO 22 I=ISTART,IEND
        IF (RLSRES(I).NE.RESISTANCE) THEN
          RESISTANCE=RLSRES(I)
          WRITE (ILU,2600) RESISTANCE
        ENDIF
        IF (RLSDEP(I).NE.RDEPTH) THEN
          RDEPTH=RLSDEP(I)
          WRITE (ILU,2700) RDEPTH
        ENDIF
        WRITE (ILU,3500) CLSZS(I),RLSH(I),RLSWID(I),ALSLAB(I)
   22   CONTINUE
        WRITE (ILU,4000) CLSZE(IEND)
   23   CONTINUE
      ENDIF
      IF (NLSH.NE.0) THEN
        ICOUNT=0
        DO 24 II=1,NLSH
        IF (KLSDN(II).NE.II) GOTO 24 ! skip, part of network
        IF (ICOUNT.EQ.0) WRITE (ILU,4500)
        ICOUNT=ICOUNT+1
        I=KLSPTH(II)
        IF (RLSRES(I).NE.RESISTANCE) THEN
          RESISTANCE=RLSRES(I)
          WRITE (ILU,2600) RESISTANCE
        ENDIF
        IF (RLSDEP(I).NE.RDEPTH) THEN
          RDEPTH=RLSDEP(I)
          WRITE (ILU,2700) RDEPTH
        ENDIF
         WRITE (ILU,5000) CLSZS(I),CLSZE(I),RLSH(I),RLSWID(I),ALSLAB(I)
   24   CONTINUE
      ENDIF
      IF (NLSIG.NE.0) THEN
        WRITE (ILU,6000)
        DO 25 II=1,NLSIG
        I=KLSPTS(II)
        WRITE (ILU,6500) CLSZS(I),CLSZE(I),RLSIG(I),RLSWID(I),ALSLAB(I)
   25   CONTINUE
      ENDIF
      WRITE (ILU,7000)
      RETURN
C      
 1000 FORMAT (' linesink')
 1500 FORMAT (' stream')
 2000 FORMAT (' stream ',F9.0)
 2500 FORMAT (' end')
 2600 FORMAT (' resistance ',G14.7)
 2700 FORMAT (' depth ',G14.7)
 3000 FORMAT (' overlandflow ',F9.0)
 3050 FORMAT ('*    x         y      head    ',
     &        '    width        label')
 3500 FORMAT (2(F9.0,1X),2(G14.7,1X),A16)
 4000 FORMAT (2(F9.0,1X))     
 4500 FORMAT (' head',/,'*  x1        y1        x2        y2       ',
     &'head       width        label')
 5000 FORMAT (4(F9.0,1X),G14.7,1X,F8.3,1X,A16)
 6000 FORMAT (' discharge',/,'*  x1        y1        x2        y2    ',
     &        '  exf. rate    width       label')
 6500 FORMAT (4(F9.0,1X),G14.7,1X,F8.3,1X,A16)
 7000 FORMAT (' quit')
      END
C
C-----------------------------------------------------------------------------------------------------------
C
      SUBROUTINE LSIO (ICODE,ILU,RVERSION,ierr)
C
C-----------------------------------------------------------------------------------------------------------
C
C     Last change:  HH    6 Dec 1999    2:48 pm
C     Routine reads or writes contents of LSCOM.INC to an 
C     external file.
C     ICODE=41 write
C     ICODE=73 read
C
      IMPLICIT NONE
      INTEGER(4) ICODE,ILU,IERR,NWORD
      LOGICAL LDUM(2)
      REAL(8) RVERSION
      INCLUDE 'lscom.inc'
C
      if (ierr.ne.0) return
      CALL BUFIO4 (NLS,1,ILU,ICODE,IERR)
      IF (NLS.EQ.0) RETURN
      CALL BUFIOC (CLSZS,NLS,ILU,ICODE,IERR)
      CALL BUFIOC (CLSZE,NLS,ILU,ICODE,IERR)
      CALL BUFIOR (RLSIG,NLS,ILU,ICODE,IERR)
      CALL BUFIOR (RLSPOT,NLS,ILU,ICODE,IERR)
      CALL BUFIOR (RLSH,NLS,ILU,ICODE,IERR)
      CALL BUFIOR (RLSRES,NLS,ILU,ICODE,IERR)
      CALL BUFIO4 (NLSIG,2,ILU,ICODE,IERR)
      CALL BUFIO4 (KLSPTS,NLSIG,ILU,ICODE,IERR)
      CALL BUFIO4 (KLSPTH,NLSH,ILU,ICODE,IERR)
      CALL BUFIOR (RLSWID,NLS,ILU,ICODE,IERR)
      CALL BUFIOA (ALSLAB,NLS,ILU,ICODE,IERR)
      CALL BUFIOR (RLSW0,3,ILU,ICODE,IERR)
      CALL BUFIOR (RLSDEP,NLS,ILU,ICODE,IERR)
      CALL BUFIOL (LSGIV,NLS,ILU,ICODE,IERR)
      IF (ICODE.EQ.41) CALL BUFIOL (LSHIRE,2,ILU,ICODE,IERR) ! to maintain
      IF (ICODE.EQ.73) CALL BUFIOL (LDUM,2,ILU,ICODE,IERR) ! highlight comm.
      CALL BUFIO4 (KLSUP,NLS,ILU,ICODE,IERR)
      CALL BUFIO4 (KLSDN,NLS,ILU,ICODE,IERR)
      CALL BUFIOR (RLSBF,NLS,ILU,ICODE,IERR)
      CALL BUFIOR (RLSOF,NLS,ILU,ICODE,IERR)
      CALL BUFIOR (RLSOFSIG,NLS,ILU,ICODE,IERR)
      CALL BUFIOR (RLSBFMX,2,ILU,ICODE,IERR)
      CALL BUFIOL (LSBASE,2,ILU,ICODE,IERR)      
      CALL BUFIO4 (NLSTRING,1,ILU,ICODE,IERR)
      IF (NLSTRING.EQ.0) GOTO 10
      CALL BUFIO4 (KLSTRING,NLSTRING,ILU,ICODE,IERR)
      CALL BUFIO4 (KLSTREND,NLSTRING,ILU,ICODE,IERR)
      CALL BUFIOR (RLSOFST,NLSTRING,ILU,ICODE,IERR)
      CALL BUFIOR (RLSOVHW,NLSTRING,ILU,ICODE,IERR)
  10  IF (RVERSION.EQ.1.0) RETURN
      IF (RVERSION.EQ.2.0) RETURN
      CALL BUFIOL (LSDRAIN,NLS,ILU,ICODE,IERR)
      CALL BUFIOR (RLSERROR,NLS,ILU,ICODE,IERR)
      IF (RVERSION.EQ.3.0) RETURN
      CALL BUFIOC (CLSCONST,NLS,ILU,ICODE,IERR)
      IF (RVERSION.EQ.4.0) RETURN
      CALL BUFIOR (RLSHMIN,NLS,ILU,ICODE,IERR)
      CALL BUFIOR (RLSQ,NLS,ILU,ICODE,IERR)
      CALL BUFIOL (LSGALLERY,NLS,ILU,ICODE,IERR)
      CALL BUFIOL (LSINLET,NLS,ILU,ICODE,IERR)
      CALL BUFIOL (LSOUTLET,NLS,ILU,ICODE,IERR)
      CALL BUFIOL (LSLAKE,NLS,ILU,ICODE,IERR)
      IF (NLSTRING.GT.0) THEN
       CALL BUFIOA (ALSTBLFILENAME,NLSTRING,ILU,ICODE,IERR)
       CALL BUFIO4 (NLSTAB,NLSTRING,ILU,ICODE,IERR)
       CALL BUFIO4 (NLSTABLENGTH,NLSTRING,ILU,ICODE,IERR)
       NWORD=NLSTBLMAX*NLSTBLSIZEMAX
       CALL BUFIOR (RLSTAGE,NWORD,ILU,ICODE,IERR)
       CALL BUFIOR (RLSTABLE,NWORD,ILU,ICODE,IERR)
       CALL BUFIOR (RLSH1,NLSTBLMAX,ILU,ICODE,IERR) ! wrong range, keep for compatibility with rversion 5 - 7
       CALL BUFIOR (RLSH2,NLSTBLMAX,ILU,ICODE,IERR) ! wrong range, keep for compatibility with rversion 5 - 7
       CALL BUFIOR (RLSQ1,NLSTBLMAX,ILU,ICODE,IERR) ! wrong range, keep for compatibility with rversion 5 - 7
       CALL BUFIOR (RLSQ2,NLSTBLMAX,ILU,ICODE,IERR) ! wrong range, keep for compatibility with rversion 5 - 7
      ENDIF
      IF (RVERSION.EQ.5.0) RETURN
       CALL BUFIOR (RLSH0,NLSTMX,ILU,ICODE,IERR)
       CALL BUFIOR (RLSEVAP,NLSTMX,ILU,ICODE,IERR)
      IF (RVERSION.EQ.6.0) RETURN
      IF (RVERSION.EQ.7.0) RETURN
       CALL BUFIO4 (nlakeiterations,1,ILU,ICODE,IERR)
      if (nlstring.gt.0) then
       CALL BUFIOR (RLSH1,NLSTRING,ILU,ICODE,IERR)
       CALL BUFIOR (RLSH2,NLSTRING,ILU,ICODE,IERR)
       CALL BUFIOR (RLSQ1,NLSTRING,ILU,ICODE,IERR)
       CALL BUFIOR (RLSQ2,NLSTRING,ILU,ICODE,IERR)
      end if
      if (rversion.le.11) return
       call bufio4 (ilsbound,nls,ilu,icode,ierr)
      RETURN
      END
C
C-----------------------------------------------------------------------------------------------------------
C
      SUBROUTINE LSNKOUT (ICALL)
C
C-----------------------------------------------------------------------------------------------------------
C
c     Output for debugging
c
      IMPLICIT NONE
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
      INTEGER(4) ICALL,I
      COMPLEX(8) CZ
      CHARACTER(1) ADUM
      SAVE
      NLS=NLSIG+NLSH
      WRITE (ILUOUT,1002) NLS,NLSMAX,NLSIG,NLSH,ILSPLOT,LSBASE,LSEND,
     &                    LSFIRST
      WRITE (ILUOUT,1004) RLSW0,RLSC0,RLSD0,RLSBFMX,RLSOFMX
      WRITE (ILUOUT,2000) ICALL
c      IF (LUOUTFILE) WRITE (ILUME,2000)
c      READ (ILUIN,2500) ADUM
c      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE (ILUOUT,1006) (I,CLSZS(I),CLSZE(I),I=1,NLS)
      WRITE (ILUOUT,2000) ICALL
c      IF (LUOUTFILE) WRITE (ILUME,2000)
c      READ (ILUIN,2500) ADUM
c      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      DO 10 I=1,NLS
      CZ=0.5*(CLSZS(I)+CLSZE(I))
      WRITE (ILUOUT,1007) CZ,I,RLSPOT(I)
  10  CONTINUE      
      WRITE (ILUOUT,2000) ICALL
c      IF (LUOUTFILE) WRITE (ILUME,2000)
c      READ (ILUIN,2500) ADUM
c      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE (ILUOUT,1010) (I,RLSIG(I),RLSOFSIG(I),RLSH(I),I=1,NLS)
      WRITE (ILUOUT,2000) ICALL
c      IF (LUOUTFILE) WRITE (ILUME,2000)
c      READ (ILUIN,2500) ADUM
c      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE (ILUOUT,1015) (I,RLSRES(I),RLSWID(I),RLSDEP(I),
     &                    I=1,NLS)
      WRITE (ILUOUT,2000) ICALL
c      IF (LUOUTFILE) WRITE (ILUME,2000)
c      READ (ILUIN,2500) ADUM
c      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE (ILUOUT,1016) (I,KLSPTS(I),KLSPTH(I),LSGIV(I),
     &                    ALSLAB(I),I=1,NLS)
      WRITE (ILUOUT,2000) ICALL
c      IF (LUOUTFILE) WRITE (ILUME,2000)
c      READ (ILUIN,2500) ADUM
c      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE(ILUOUT,2900)(I,KLSUP(I),KLSDN(I),RLSBF(I),RLSOF(I),I=1,NLS)
      WRITE (ILUOUT,2000) ICALL
c      IF (LUOUTFILE) WRITE (ILUME,2000)
c      READ (ILUIN,2500) ADUM
c      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE (ILUOUT,3000) NLSTRING
      IF (NLSTRING.EQ.0) RETURN
      WRITE (ILUOUT,3100) (I,KLSTRING(I),KLSTREND(I),RLSOFST(I),
     &                    RLSOVHW(I),I=1,NLSTRING)
      WRITE (ILUOUT,2000) ICALL
c      IF (LUOUTFILE) WRITE (ILUME,2000)
c      READ (ILUIN,2500) ADUM
c      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE (ILUOUT,3000) NLSTRING
      IF (NLSTRING.EQ.0) RETURN
      RETURN
 1002 FORMAT (' LS: NLS,NLSMAX,NLSIG,NLSH,ILSPLOT ',5I5,/,
     &        ' LSBASE,LSEND,LSFIRST ',3L5)
 1004 FORMAT (' LS: RLSW0,RLSC0,RLSD0,RLSBFMX,RLSOFMX ',5(G11.4))
 1006 FORMAT (1(' LS:I,CLSZS,CLSZE ',I3,2(G11.4,G11.4,1X)))
 1007 FORMAT (' LS:CZ0=(',2G14.7,') RLSPOT(',I3,')=',G14.7)
 1010 FORMAT (1(' LS:I,RLSIG,RLSOFSIG,RLSH ',I3,1X,G14.7,G14.7,G11.4)) 
 1015 FORMAT (1(' LS:I,RLSRES,RLSWID,RLSDEP ',I3,3(G11.4)))
 1016 FORMAT (1(' LS:I,KLSPTS,KLSPTH,LSGIV,ALSLAB ',3(I3),L3,2X,A16))
 2000 FORMAT (' LS: ICALL=',I3,' press <Enter> to continue.')
 2500 FORMAT (A1)
 2900 FORMAT (1(' LS:I,KLSUP,KLSDN,RLSBF,RLSOF ',3I5,2G11.4))
 3000 FORMAT (' LS: NLSTRING=',I3)
 3100 FORMAT (1(' LS: I,KLSTRING,KLSTREND,RLSOFST,RLSOVHW '
     &      ,3I4,1X,G11.4,2X,G11.4))
      END
