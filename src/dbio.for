C     Last change:  HMH  18 Dec 2006    2:17 pm
c --------------------------------------------------------------------------------
c
c     This file contains the following routines or functions
c
c     DBIO          write/reads contents of common blocks to disk
c     INHOMEXTRACT  writes all line doublet data to file *.xtr
c     DBDATIO       writes all inhomogeneity and slurry wall (horizontal barrier) data to *.dat file
c     DBOUT         formatted write of all commomblock data to a specified logical unit (debugging) 
c       
c
c --------------------------------------------------------------------------------------------------------
c
      SUBROUTINE DBIO (ICODE,ILU,RVERSION,ierr)
c
c --------------------------------------------------------------------------------------------------------
c
C
C     Routine writes/reads contents of doublet common blocks to disk
C
      IMPLICIT NONE
      INTEGER(4) ICODE,ILU,IERR
      LOGICAL LDUM
      REAL(8) RVERSION
      INCLUDE 'dbcom.inc'
c
      if (ierr.ne.0) return
      CALL BUFIO4 (NDB,1,ILU,ICODE,IERR)
      IF (NDB.EQ.0) RETURN
      CALL BUFIO4 (NDBSTR,1,ILU,ICODE,IERR)
      CALL BUFIO4 (NDBSTI,NDBSTR,ILU,ICODE,IERR)
      CALL BUFIO4 (IDBSTA,NDBSTR,ILU,ICODE,IERR)
      CALL BUFIOC (CDBZ,NDB,ILU,ICODE,IERR)
      CALL BUFIOD (DBSTR,NDB,ILU,ICODE,IERR)
      CALL BUFIOD (DBQSTR,NDB,ILU,ICODE,IERR)
      CALL BUFIOR (RDBTFACNOD,NDB,ILU,ICODE,IERR)
      CALL BUFIOA (ADBLAB,NDB,ILU,ICODE,IERR)
      CALL BUFIOR (RDBK,NDBSTR,ILU,ICODE,IERR)
      CALL BUFIOR (RDBP,NDBSTR,ILU,ICODE,IERR)
      CALL BUFIOR (RDBB,NDBSTR,ILU,ICODE,IERR)
      CALL BUFIOR (RDBT,NDBSTR,ILU,ICODE,IERR)
      CALL BUFIOD (DBAVS,NDBSTR,ILU,ICODE,IERR)
      CALL BUFIOL (LDBINS,NDBSTR,ILU,ICODE,IERR)
      IF (ICODE.EQ.41) CALL BUFIOL (LDBPLT,1,ILU,ICODE,IERR) ! to maintain
      IF (ICODE.EQ.73) CALL BUFIOL (LDUM,1,ILU,ICODE,IERR) ! plot command
      CALL BUFIOR (RDBGAM,NDBSTR,ILU,ICODE,IERR)
      CALL BUFIOD (DBAREA,NDBSTR,ILU,ICODE,IERR)      
      CALL BUFIOR (RDBRAD,NDBSTR,ILU,ICODE,IERR)      
      CALL BUFIOC (CDBZ0,NDBSTR,ILU,ICODE,IERR)      
      CALL BUFICD (CDDBSS,NDB,ILU,ICODE,IERR)
      CALL BUFIOD (DBQSTS,NDB,ILU,ICODE,IERR)      
      IF (RVERSION.EQ.1.0) RETURN
      IF (RVERSION.EQ.2.0) RETURN
      CALL BUFIO4 (IDOMAINTYPE,NDBSTR,ILU,ICODE,IERR)
      CALL BUFIOR (RDBW,NDBSTR,ILU,ICODE,IERR)
      CALL BUFIOR (RDBTFACCTR,NDB,ILU,ICODE,IERR)
      CALL BUFIOR (RDBBFACNOD,NDB,ILU,ICODE,IERR)
      CALL BUFIOR (RDBBFACCTR,NDB,ILU,ICODE,IERR)
      CALL BUFIOR (RDBERRSTRT,NDB,ILU,ICODE,IERR)
      CALL BUFIOR (RDBERRCNTR,NDB,ILU,ICODE,IERR)
      CALL BUFIOR (RDBTFACCTR,NDB,ILU,ICODE,IERR)
      CALL BUFIOR (RDBBFACNOD,NDB,ILU,ICODE,IERR)
      CALL BUFIOR (RDBBFACCTR,NDB,ILU,ICODE,IERR)
      IF (RVERSION.LE.8) RETURN
      CALL BUFIOR (RDBHNOD,NDB,ILU,ICODE,IERR)
      CALL BUFIOR (RDBHCTR,NDB,ILU,ICODE,IERR)
      IF (RVERSION.LE.10) RETURN
      CALL BUFIOR (DBKNON,NDB,ILU,ICODE,IERR)
      CALL BUFIOR (DBKNOC,NDB,ILU,ICODE,IERR)
      if (rversion.le.11) return
      call bufior (dbeps,1,ilu,icode,ierr)
      CALL BUFIOR (rdboffsetinhom,1,ILU,ICODE,IERR)
      call bufiol (ldbrechargeonly,ndbstr,ilu,icode,ierr)
      call bufio4 (idbcode,ndb,ilu,icode,ierr)
      call bufior (dbsen,ndb,ilu,icode,ierr)
      CALL BUFIOR (rdberrend,NDB,ILU,ICODE,IERR)
      CALL BUFIOR (rdbhend,NDB,ILU,ICODE,IERR)
      CALL BUFIOR (rdbtfacend,NDB,ILU,ICODE,IERR)
      CALL BUFIOR (rdbbfacend,NDB,ILU,ICODE,IERR)
      CALL BUFIOR (dbknoe,NDB,ILU,ICODE,IERR)
      CALL BUFIOR (rdbki,NDB,ILU,ICODE,IERR)
      CALL BUFIOR (rdbko,NDB,ILU,ICODE,IERR)
      CALL BUFIOR (rdbbi,NDB,ILU,ICODE,IERR)
      CALL BUFIOR (rdbbo,NDB,ILU,ICODE,IERR)
      RETURN
      END
c
c --------------------------------------------------------------------------------------------------------
c
      subroutine inhomextract(ilu)
c
c --------------------------------------------------------------------------------------------------------
c
c
c     Routine allows data from inhomogeneity module to be extracted
c
      IMPLICIT NONE
      INTEGER(4) ILU,ISTR,INOD1,INODL,INODM1,INOD
      COMPLEX(8) CFDBF,CFDBS,CFDBG,CDDS,CDDBOM,CZCNTR
      INCLUDE 'dbcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      IF (NDB.EQ.0)  RETURN
      write (ilu,5000)
      DO ISTR=1,NDBSTR
      INOD1=IDBSTA(ISTR)
      INODL=INOD1+NDBSTI(ISTR)-1
      INODM1=INODL
c
c      hydraulic conductivity inhomogeneities
c
      if (idomaintype(istr).eq.1) then
      write (ilu,1000)
      write (ilu,1010) istr,rdbk(istr),rdbb(istr),rdbgam(istr),
     &                 rdbp(istr),dbarea(istr),dbavs(istr)
      write (ilu,1020)
      DO INOD=INOD1,INODL
      if (inod.eq.inodl) then
      czcntr=0.5*(cdbz(inod)+cdbz(inod1))
      else
      czcntr=0.5*(cdbz(inod)+cdbz(inod+1))
      endif
      write (ilu,1030) inod,cdbz(inod),czcntr,dbstr(inod),dbqstr(inod),
     &                 cddbss(inod),dbqsts(inod),
     &               rdbtfacnod(inod),rdberrstrt(inod),rdberrcntr(inod),  ! note: not writing the error at the last collocation point if it exists.
     &                 adblab(inod)
      END do
      endif
c
c      Open slurry walls
c
      if (idomaintype(istr).eq.2) THEN
      write (ilu,2000)
      write (ilu,2010) istr,rdbk(istr),rdbw(istr),rdbp(istr),rdbb(istr),
     &                 adblab(inod1)
      write (ilu,2020)
      DO INOD=INOD1,INODL
      czcntr=0.5*(cdbz(inod)+cdbz(inod+1))
      write (ilu,2030) inod,cdbz(inod),czcntr,dbstr(inod),dbqstr(inod),
     &                 rdberrcntr(inod),rdberrstrt(inod),adblab(inod)
      END DO
      end if
c
c      Closed slurry walls
c
      if (idomaintype(istr).eq.3) THEN
      write (ilu,3000)
      write (ilu,2010) istr,rdbk(istr),rdbw(istr),rdbp(istr),rdbb(istr),
     &                 adblab(inod1)
      write (ilu,2020)
      DO INOD=INOD1,INODL
      if (inod.eq.inodl) then
      czcntr=0.5*(cdbz(inodl)+cdbz(inod1))
      else
      czcntr=0.5*(cdbz(inod)+cdbz(inod+1))
      end if
      write (ilu,2030) inod,cdbz(inod),czcntr,dbstr(inod),dbqstr(inod),
     &                 rdberrcntr(inod),rdberrstrt(inod),adblab(inod)
      END do
      end if
      END DO
      return
c
 1000 format ('! transmissivity inhomogeneity domain.',/,
     &'* string# conductivity  bottom elev.   extract. rate   porosity',
     &        '      domain area (L^2)   average potential jump')
 1010 format (i7,4(',',E14.7),2(',',E21.14))
 1020 format ('* node #       x1             y1             xc        ',
     &       '     yc       str. pot. jump node  parab. str. pot. jump',
     &       '  real str. extr. node  imag. str. extr. node',
     &        ' real str. extr. ',
     &        'cntr   Ti/(Ti-T0)   delta P node %  delta P cntr % ',
     &        ' label')
 1030 format (i7,4(',',E14.7),5(',',E21.14),3(',',E14.7),',',a16)
 2000 format ('! open slurry wall.',/,
     &'* string# conductivity       width       porosity       ',
     &        'wall bottom elev.   label')
 2010 format (i7,4(',',E14.7),2x,a16)
 2020 format ('* node #       x1             y1             xc        ',
     &       '     yc       str. pot. jump node  parab. str. pot. jump',
     &      'delta flow cntr % delta flow node %  label')
 2030 format (i7,4(',',E14.7),2(',',E21.14),2(',',E14.7),',',3x,a16)
 3000 format ('! closed slurry wall.',/,
     &'* string# conductivity       width       porosity       ',
     &        'wall bottom elev.   label')
 5000 format ('! line doublet strings.')
c
C     Note: For "Open Slurry Walls" the last node is not printed, but it is used in calculating the center.
c     The doublet strength for the last node of an "Open Slurry Wall" is equal to zero.
c     The "Delta Flow" values are the flows across the center or node integration intervals as follows:
c     from the doublet strenght minus the actual flow found across that line.
c
      end subroutine
c
c
c --------------------------------------------------------------------------------------------------------
c
      SUBROUTINE DBDATIO (ILU)
c
c --------------------------------------------------------------------------------------------------------
c
C
C     Routine writes all inhomogeneity data to a ".dat" file.
C     Routine is called in DATIO.
C
      IMPLICIT NONE
      INTEGER(4) ILU,ISTR,INOD1,INODL,INOD
      INCLUDE 'dbcom.inc'
      INCLUDE 'lusys.inc'
C
      IF (NDB.EQ.0) RETURN
      WRITE (ILU,1000)
      DO ISTR=1,NDBSTR
      INOD1=IDBSTA(ISTR)
      INODL=INOD1+NDBSTI(ISTR)-1
      IF (IDOMAINTYPE(ISTR).EQ.1) THEN ! hydraulic conductivity inhomogeneities
      WRITE (ILU,2000) RDBK(ISTR),RDBGAM(ISTR),RDBP(ISTR)
      DO INOD=INOD1,INODL
      WRITE (ILU,3000) CDBZ(INOD),ADBLAB(INOD)
      END DO
      ENDIF
      IF (IDOMAINTYPE(ISTR).EQ.2) THEN ! open slurry walls
      WRITE (ILU,4000) RDBK(ISTR),RDBW(ISTR),RDBP(ISTR),RDBB(ISTR)
      DO INOD=INOD1,INODL+1
      WRITE (ILU,3000) CDBZ(INOD),ADBLAB(INOD)
      END DO
      END IF
      IF (IDOMAINTYPE(ISTR).EQ.3) THEN ! closed slurry walls
      WRITE (ILU,5000) RDBK(ISTR),RDBW(ISTR),RDBP(ISTR),RDBB(ISTR)
      DO INOD=INOD1,INODL
      WRITE (ILU,3000) CDBZ(INOD),ADBLAB(INOD)
      END DO
      END IF
      END DO
      WRITE (ILU,6000)
      RETURN
C      
 1000 FORMAT (' inhomogeneity')
 2000 FORMAT ('*        hydraul. cond.   added exf.rate   porosity',/,
     &        ' inhom ',3(2X,G14.7))
 3000 FORMAT (2(F14.4,2X),A16)
 4000 FORMAT ('*                 hydraul. cond.    width       ',
     &        'porosity       bottom elev.',/,
     &        ' slurry open   ',4(1X,G14.7))
 5000 FORMAT ('*                 hydraul. cond.    width       ',
     &        'porosity       bottom elev.',/,
     &        ' slurry closed ',4(1X,G14.7))
 6000 FORMAT (' quit')
      END
c
c --------------------------------------------------------------------------------------------------------
c
      SUBROUTINE DBOUT (ICALL)
c
c --------------------------------------------------------------------------------------------------------
c
C
C     Routine writes the contents of DBCOM common blocks
C      
      IMPLICIT NONE
      INTEGER(4) ICALL,I,J
      CHARACTER(1) ADUM
      INCLUDE 'dbcom.inc'
      INCLUDE 'lusys.inc'
C
      WRITE (ILUOUT,1000) NDB,NDBSTR,NDBZMX,NDBTRM,NDBSMX,LDBPLT,
     &                    DPI,D2PI,CDI,DBEPS,CDDPAD
      WRITE (ILUOUT,1010) (I,NDBSTI(I),I=1,NDBSTR)
      WRITE (ILUOUT,1020) (I,IDBSTA(I),I=1,NDBSTR)
      WRITE  (ILUOUT,2000)  ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000) ICALL
      READ (ILUIN,3000)  ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE (ILUOUT,1030) (I,RDBK(I),I=1,NDBSTR)
      WRITE (ILUOUT,1040) (I,RDBP(I),I=1,NDBSTR)
      WRITE (ILUOUT,1050) (I,RDBB(I),I=1,NDBSTR)
      WRITE (ILUOUT,1060) (I,RDBT(I),I=1,NDBSTR)
      WRITE (ILUOUT,1070) (I,LDBINS(I),I,DBAVS(I),I=1,NDBSTR)
      WRITE  (ILUOUT,2000)  ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000) ICALL
      READ (ILUIN,3000)  ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE (ILUOUT,1080) (I,CDBZ(I),I,DBSTR(I),I=1,NDB)      
      WRITE  (ILUOUT,2000)  ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000) ICALL
      READ (ILUIN,3000)  ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE (ILUOUT,1090) (I,CDBZEL(I),I,DBQSTR(I),I=1,NDB)
      WRITE  (ILUOUT,2000)  ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000) ICALL
      READ (ILUIN,3000)  ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE (ILUOUT,1100) (I,LDBFAR(I),I,LDBNOD(I),I,ADBLAB(I),
     &I,CDDBLN(I),I=1,NDB)
      WRITE  (ILUOUT,2000)  ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000) ICALL
      READ (ILUIN,3000)  ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE (ILUOUT,1110) (I,CDBZ21(I),I,RDBTFACNOD(I),I=1,NDB)
      WRITE  (ILUOUT,2000)  ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000) ICALL
      READ (ILUIN,3000)  ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      DO 10 I=1,NDB
      WRITE (ILUOUT,1120) (J,I,CDBZIV(J,I),J=1,NDBTRM)
  10  CONTINUE      
      WRITE  (ILUOUT,2000)  ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000) ICALL
      READ (ILUIN,3000)  ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      DO 20 I=1,NDBSTR
      WRITE (ILUOUT,1130) I,RDBGAM(I),DBAREA(I),CDBZ0(I),RDBRAD(I)
  20  CONTINUE
      WRITE  (ILUOUT,2000)  ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000) ICALL
      READ (ILUIN,3000)  ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      DO 30 I=1,NDB
      WRITE (ILUOUT,1140) I,CDDBSS(I),I,DBQSTS(I)
  30  CONTINUE
      WRITE  (ILUOUT,2000)  ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000) ICALL
      READ (ILUIN,3000)  ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN          
      RETURN
C
 1000 FORMAT (' DBOUT: NDB,NDBSTR,NDBZMX,NDBTRM,NDBSMX,LDBPLT ',5I5,L3,
     &/,' DPI,D2PI=',2(D24.16),/,' CDI=',2(D24.16),' DBEPS=',D24.16,/,
     &' CDDPAD=',2(D24.16))
 1010 FORMAT (' DBOUT:',3(' NDBSTI(',I3,')=',I5)/)
 1020 FORMAT (' DBOUT:',3(' IDBSTA(',I3,')=',I5)/)
 1030 FORMAT (' DBOUT:',3(' RDBK(',I3,')=',G14.7)/)
 1040 FORMAT (' DBOUT:',3(' RDBP(',I3,')=',G14.7)/)
 1050 FORMAT (' DBOUT:',3(' RDBB(',I3,')=',G14.7)/)
 1060 FORMAT (' DBOUT:',3(' RDBT(',I3,')=',G14.7)/)
 1070 FORMAT (' DBOUT:',3(' LDBINS(',I3,')=',L2,
     &                   ' DBAVS(',I3,')=',D24.16))
 1080 FORMAT (' DBOUT: CDBZ(',I3,')=',2G11.4,' DBSTR(',I3,')=',
     &          D24.16,/)
 1090 FORMAT (' DBOUT: CDBZEL(',I3,')=',2G11.4,' DBQSTR(',I3,')=',
     &          D24.16,/)
 1100 FORMAT (' DBOUT: LDBFAR(',I3,')=',L3,' LDBNOD(',I3,')=',L3,
     &' ADBLAB(',I3,')=',A16,/,' CDDBLN(',I3,')=',2(D24.16),/)
 1110 FORMAT (' DBOUT: CDBZ21(',I3,')=',2G11.4,' RDBTFACNOD(',I3,')=',
     &          G11.4,/)
 1120 FORMAT (' DBOUT:',2(' CDBZIV(',I2,',',I3,')=',2G11.4),/)
 1130 FORMAT (' DBOUT: RDBGAM(',I3,')=',G11.4,' DBAREA=',D24.16,/
     &        ' CDBZ0=',2G11.4,' RDBRAD=',G11.4)
 1140 FORMAT (' DBOUT: CDDBSS(',I3,')=',2(D24.16),/,
     &        '        DBQSTS(',I3,')=',D24.16)
 2000 FORMAT (' DBOUT: ICALL=',I3)
 3000 FORMAT (A1)
      END


