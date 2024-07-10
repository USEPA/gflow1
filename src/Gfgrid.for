C     Last change:  HMH   6 May 2008    3:38 pm
c     This file contains the following routines and functions
c
c     BLOCK DATA GRIDINIT   set commonblock variables for GRID.INC
c     SUBROUTINE GRINIT     grid initialization routine (for contour plots)
c     SUBROUTINE GRIDIN     driver for filling grid with function values
c     SUBROUTINE GFFILL     adds function values to grid
c          ENTRY PFILL      same as GFFILL, but plots a dot map
c     SUBROUTINE MINMAX     determines minimum and maximum values in the grid
c     SUBROUTINE GRIDIO     reads or writes contents of common /GRID/ to .SOL file
c     SUBROUTINE RAIO       reads or writes the (grid) array RA to a binary .GRI file
c     SUBROUTINE GRIO       reads or writes contents of common /GRID/ to .GRI file
c     SUBROUTINE GRTEMPIO   reads data from .GRI files into dummy storage for subtracting grids
c     SUBROUTINE RATEMPIO   reads or writes the array (grid) RA2 to a file
c     SUBROUTINE DIFFGRID   reads a grid and subtracts subtracts it from the current grid
c     SUBROUTINE GRIDOUT    writes the contents of common /GRID/ formatted for debugging purposes
c
c -------------------------------------------------------------------------------
c
      BLOCK DATA GRIDINIT
      IMPLICIT NONE
      INCLUDE 'GRID.INC'
      DATA NXMAX,NYMAX /2*200/ ! see array in /BASEG/, see GRIDIN.FOR
      DATA RX1,RY1,RX2,RY2 /4*0.0/
      DATA NX,NY /2*0/
      DATA IX1,IY1,IX2,IY2 /-1,-1,1,1/
      DATA LHEAD /.TRUE./
      DATA LPHI,LPSI,LDISCH,LFLOWNET,LBASE,LTOP /6*.FALSE./
      DATA LHEIGHT,LPOROSITY,LRECHARGE,LINTELV /4*.FALSE./
      DATA LEAKANCE,LEAKAGE /2*.FALSE./
      END
c
c
c
c -------------------------------------------------------------------------------
c
      SUBROUTINE GRINIT (RA,ISIZE)
c
c -------------------------------------------------------------------------------
c
C                                        
C     Routine initializes the matrix RA and generates the x, y, z value
C     arrays: RXG, RYG, R3DCPZ.
C      
      IMPLICIT NONE
      INTEGER(4) ISIZE,I,J,NINCX,NINCY,INUM
      LOGICAL LSET,LDBOUTSIDE
      REAL(8) RA,RDELX,RDELY,RDELZ
      INCLUDE 'GRID.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
      DIMENSION RA(ISIZE,*)
      RWX1=0.0
      RWX2=0.0
      RWY1=0.0
      RWY2=0.0
      LGRCOM=.FALSE.
      LSET=NX.GT.NXMAX
      IF(LSET) THEN
        LSET=.FALSE.
        WRITE (ILUER,1000) NXMAX,NX,NXMAX
        NX=NXMAX
        IF (LUCON) CALL INLINE
      ENDIF
 3    NINCX=NX-1
      RDELX=(RX2-RX1)/FLOAT(NINCX)
      NY=INT((RY2-RY1)/RDELX)
      IF(NY.LE.NYMAX) GOTO 5
      INUM=NX*NYMAX
      NX=INUM/NY
      LSET=.TRUE.
      GOTO 3
 5    NINCY=NY-1
      IF (LSET) THEN
        LSET=.FALSE.
        WRITE (ILUER,2000) NYMAX,NX,NY
        IF (LUCON) CALL INLINE
      ENDIF
      IF (L3DPL.AND.L3DVER) THEN
      RDELX=(R3DX2-R3DX1)/FLOAT(NINCX)
      RDELY=(R3DY2-R3DY1)/FLOAT(NINCX)
      RDELZ=(R3DZ2-R3DZ1)/FLOAT(NINCY)
      DO 20 I=1,NX
      RXG(I)=(I-1)*RDELX+R3DX1
      RYG(I)=(I-1)*RDELY+R3DY1
      DO 10 J=1,NY
      R3DCPZ(J)=(J-1)*RDELZ+R3DZ1
      RA(I,J)=0.0
  10  CONTINUE
  20  CONTINUE
      RETURN
      ENDIF
      RDELY=(RY2-RY1)/FLOAT(NINCY)
      DO 40 I=1,NX
      RXG(I)=(I-1)*RDELX+RX1
      DO 30 J=1,NY
      RYG(J)=(J-1)*RDELY+RY1
      IF (LDBOUTSIDE(RXG(I),RYG(J))) THEN ! flag point if outside model area
      RA(I,J)=-9999.0
      ELSE
      RA(I,J)=0.0
      END IF
 30   CONTINUE
 40   CONTINUE
      RETURN
 1000 FORMAT (' ***WARNING: in grid module:',/,
     &' Maximum number of HORIZONTALPOINTS = ',I4,/,
     &' Specified ',I4,' points reduced to ',I4,/,
     &' Press any key to continue.')
 2000 FORMAT (' ***WARNING: in grid module:',/,
     &' Too many vertical points generated due to aspect ratio of plot!'
     &,/,' Number of vertical points is limited to ',I4,/,
     &' New choice of HORIZONTALPOINT = ',I3,' which leads to ',I4,
     &' vertical points.',/,' Press any key to continue.')
      END
c
c -------------------------------------------------------------------------------
c
      SUBROUTINE GRIDIN (RA,IMXSZE)
c
c -------------------------------------------------------------------------------
c
C
C     Grid module for generating grids and grid IO
C      
      IMPLICIT NONE
      INTEGER(4) IMXSZE
      LOGICAL LESCAP
      REAL(8) RA,RFHEAD,RFPOT,RFPSI,RFDSCH,rflklowerhead,
     &        RFBASE,RFTOP,RFHGHT,RFPERM,RFPOR,RFDBGM,
     &        RA2
      COMPLEX(8) CZ
      character(1) adum
      DIMENSION ra(imxsze,*)
      EXTERNAL RFHEAD,RFPOT,RFPSI,RFDSCH,RFINTERFACE,
     &         RFBASE,RFTOP,RFHGHT,RFPERM,RFPOR,RFDBGM,
     &         rflklowerhead
      COMMON /BASEG/ RA2(200,200)       ! for flownet or drawdown plots
      INCLUDE 'GRID.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
C
  10  CZ=(0.0,0.0)
      CALL GGRPAR(LESCAP,RA,IMXSZE)
      IF (LESCAP) RETURN
      CALL GRINIT (RA,IMXSZE)
      NLEVELS=0
      NTEXTS=0
      IF (LHEAD) CALL PFILL (RA,IMXSZE,RFHEAD)
c      IF (LPHI) CALL PFILL (RA,IMXSZE,RFPOT)
       if (lphi) call pfill (ra,imxsze,rflklowerhead) ! temoporary use of "plot potential" option
      IF (LPSI) CALL PFILL (RA,IMXSZE,RFPSI)
      IF (LDISCH) CALL PFILL (RA,IMXSZE,RFDSCH)
      IF (LINTELV) CALL PFILL (RA,IMXSZE,RFINTERFACE)
      IF (LFLOWNET) THEN
      CALL PFILL (RA,IMXSZE,RFPOT)
      CALL GRINIT (RA2,200)
      CALL GFFILL (RA2,200,RFPSI)
      ENDIF
      IF (LBASE) CALL PFILL (RA,IMXSZE,RFBASE)
      IF (LTOP) CALL PFILL (RA,IMXSZE,RFTOP)
      IF (LHEIGHT) CALL PFILL (RA,IMXSZE,RFHGHT)
      IF (LCONDUCTIVITY) CALL PFILL (RA,IMXSZE,RFPERM)
      IF (LPOROSITY) CALL PFILL (RA,IMXSZE,RFPOR)
      IF (LRECHARGE) CALL PFILL (RA,IMXSZE,RFDBGM)
      RWX1=RX1
      RWX2=RX2
      RWY1=RY1
      RWY2=RY2
      LGRCOM=.TRUE.
      GOTO 10
      END
c
c -------------------------------------------------------------------------------
c
      SUBROUTINE GFFILL (RA,ISIZE,RFUNC)
c
c -------------------------------------------------------------------------------
c
C
C     Routine fills the matrix RA with values generated by RFUNC.
C     Entry PFILL invokes the plotting of dots at x,y,z in the domain
C     defined by the window parameters.
C      
      IMPLICIT NONE
      INTEGER(4) ISIZE,I,J
      LOGICAL LDISPL
      REAL(8) RA,RFUNC,RDUM
      COMPLEX(8) CZ
      INCLUDE 'GRID.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
      EXTERNAL RFUNC
      DIMENSION RA(ISIZE,*)
  1   LDISPL=.FALSE.
      GOTO 2
C     ----------------------------
      ENTRY PFILL (RA,ISIZE,RFUNC)
C     ----------------------------
      IF (.NOT.LSCRN.OR..NOT.LDOTMAP) GOTO 1
      LDISPL=.TRUE.
  2   CONTINUE
      if (lucon) WRITE (ilume,1000)
      IF (L3DPL.AND.L3DVER) THEN
        DO 20 I=1,NX
        CZ=CMPLX(RXG(I),RYG(I))
        DO 10 J=1,NY
        IF (RA(I,J).NE.-9999.0) THEN   ! skip point outside model area
          R3DZ=R3DCPZ(J)
          RDUM=0.0
          RDUM=RDUM+RFUNC(CZ)
          if (lucon) WRITE (ILUME,2000)
          RA(I,J)=RA(I,J)+RDUM
        END IF
  10    CONTINUE
  20    CONTINUE

        RETURN
      ENDIF
      R3DZ=R3DH
      IF (L3DPL.AND.L3DHOR) R3DZ=R3DZ1
      DO 40 I=1,NX
      DO 30 J=1,NY
      IF (RA(I,J).NE.-9999.0) THEN ! skip points outside model domain
        CZ=CMPLX(RXG(I),RYG(J))
        RDUM=0.0
        RDUM=RDUM+RFUNC(CZ)
        if (lucon) WRITE (ILUME,2000)
        RA(I,J)=RA(I,J)+RDUM
        IF (LDISPL) THEN
        ENDIF
      END IF
  30  CONTINUE
  40  CONTINUE
        RETURN
 1000 FORMAT (' Gridding in progress, press Ctrl-Break to abort.',/)      
 2000 FORMAT ('+') 
 3000 FORMAT ('+Press <F7> for hardcopy or <Esc> to return to menu.')
      END
c
c -------------------------------------------------------------------------------
c
      SUBROUTINE MINMAX (RA,ISIZE)
c
c -------------------------------------------------------------------------------
c
C      
C     Routine assigns the extreme values of RA to RMAX and RMIN.
C
      IMPLICIT NONE
      INTEGER(4) ISIZE,I,J
      REAL(8) RA,R
      INCLUDE 'GRID.INC'
      INCLUDE 'TRACOM.INC'
      DIMENSION RA(ISIZE,*)
      RMIN=1.E31
      RMAX=-1.E31
      DO I=1,NX
       DO J=1,NY
        R=RA(I,J)
        IF (R.NE.-9999) THEN ! skip points outside model area
          RMIN=MIN(R,RMIN)
          RMAX=MAX(R,RMAX)
        END IF
       END DO
      END DO
      RETURN
      END
c
c -------------------------------------------------------------------------------
c
      SUBROUTINE GRIDIO (ICODE,ILU,RVERSION,ierr)
c
c -------------------------------------------------------------------------------
c
C
C     Routine reads or writes contents of common /GRID/ to
C     an external file.
C     ICODE=41 write
C     ICODE=73 read
C     Routine used only when reading or writing .SOL files
C
      IMPLICIT NONE
      INTEGER(4) ICODE,ILU,IERR
      REAL(8) RVERSION
      INCLUDE 'GRID.INC'
      INCLUDE 'LUSYS.INC'
c
      if (ierr.ne.0) return
      CALL BUFIOA (ATITLE,1,ILU,ICODE,IERR)
      CALL BUFIOR (RX1,6,ILU,ICODE,IERR)
      CALL BUFIOR (RSTEP,1,ILU,ICODE,IERR)
      CALL BUFIOL (LPHI,3,ILU,ICODE,IERR)
      CALL BUFIOL (LRESP,25,ILU,ICODE,IERR)
      CALL BUFIO4 (ILAYER,1,ILU,ICODE,IERR)
      CALL BUFIOL (LSCRN,5,ILU,ICODE,IERR)
      CALL BUFIOR (RVFAC,3,ILU,ICODE,IERR)
      CALL BUFIOL (LCURS,1,ILU,ICODE,IERR)
      CALL BUFIOR (RDS0,1,ILU,ICODE,IERR)
      CALL BUFIOR (RXMN,4,ILU,ICODE,IERR)
      CALL BUFIOL (LDISCH,1,ILU,ICODE,IERR)
      IF (ICODE.EQ.73) LGRCOM=.FALSE.
      IF (RVERSION.EQ.1.0) RETURN
      IF (RVERSION.EQ.2.0) RETURN
      RETURN
      END
c
c -------------------------------------------------------------------------------
c            
      SUBROUTINE RAIO (ICODE,ILU,RA,ISIZE)
c
c -------------------------------------------------------------------------------
c
C
C     Routine reads or writes the array RA to a file
C
      IMPLICIT NONE
      INTEGER(4) ICODE,ILU,ISIZE,IERR,J
      REAL(8) RA
      INCLUDE 'GRID.INC'
      DIMENSION RA(ISIZE,*)
      IF (NX.LE.0) RETURN
      DO 10 J=1,NY
      CALL BUFIOR (RA(1,J),NX,ILU,ICODE,IERR)
  10  CONTINUE
      RETURN
      END
c
c -------------------------------------------------------------------------------
c      
      SUBROUTINE GRIO (ICODE,ILU)
c
c -------------------------------------------------------------------------------
c
C
C     Routine reads or writes contents of common /GRID/ to
C     an external file.
C     ICODE=41 write
C     ICODE=73 read
C     Routine used only when reading or writing .GRI files
C
      IMPLICIT NONE
      INTEGER(4) ICODE,ILU,IERR
      INCLUDE 'GRID.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'LUSYS.INC'
      CALL BUFIOR (RX1,6,ILU,ICODE,IERR)
      CALL BUFIO4 (NX,2,ILU,ICODE,IERR)
      CALL BUFIOR (RSCALE,201,ILU,ICODE,IERR) ! Note, incomplete and 
c      unnecessary write of the arrays RXG(1000) and RYG(1000) !!!
      CALL BUFIO4 (IX1,4,ILU,ICODE,IERR)
      CALL BUFIOL (LPHI,3,ILU,ICODE,IERR)
      CALL BUFIO4 (NDOTX,4,ILU,ICODE,IERR)
      IF (ICODE.EQ.41) CALL BUFIO4 (ILAYER,1,ILU,ICODE,IERR)      
      IF (ICODE.EQ.73) CALL BUFIO4 (IGRIDCODE,1,ILU,ICODE,IERR)
      CALL BUFIOL (LSCRN,5,ILU,ICODE,IERR)
      CALL BUFIOR (RVFAC,3,ILU,ICODE,IERR)
      CALL BUFIOL (LCURS,1,ILU,ICODE,IERR)
      CALL BUFIOR (RDS0,1,ILU,ICODE,IERR)
      CALL BUFIOR (RWX1,4,ILU,ICODE,IERR)
      CALL BUFIOR (RXMN,4,ILU,ICODE,IERR)
      CALL BUFIOL (LDISCH,1,ILU,ICODE,IERR)
      CALL BUFIOR (R3DX1,6,ILU,ICODE,IERR)
      CALL BUFIOL (L3DPL,3,ILU,ICODE,IERR)
      IF (ICODE.EQ.73) LGRCOM=.TRUE.
      RETURN
      END
c
c -------------------------------------------------------------------------------
c      
      SUBROUTINE GRTEMPIO (ICODE,ILU,LRET)
c
c -------------------------------------------------------------------------------
c
C
C     Routine reads data from .GRI files into dummy storage
C     locations for comparison with previously read data
C     Routine is called in PLOPAR (grid module) as part of the
C     SUBTRACT command.
C
C     ICODE=41 write
C     ICODE=73 read
C
      IMPLICIT NONE
      INTEGER(4) IDUM1(2),IDUM2(4),IDUM3(4),IDUM4,ICODE,ILU,IERR
      LOGICAL LDUM1(3),LDUM2(5),LDUM3,LDUM4,LRET,LDUM5(3)
      REAL(8) RDUM1(6),RDUM2(201),RDUM3(3),RDUM4,RDUM5(4),RDUM6(4),
     &     RDUM7(6)
      INCLUDE 'GRID.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'COM3D.INC'
      IF (ICODE.NE.73) THEN
       WRITE (ILUER,1000) ICODE
       LRET=.TRUE.
       RETURN
      ENDIF      
      CALL BUFIOR (RDUM1,6,ILU,ICODE,IERR)
      CALL BUFIO4 (IDUM1,2,ILU,ICODE,IERR)
      CALL BUFIOR (RDUM2,201,ILU,ICODE,IERR)
      CALL BUFIO4 (IDUM2,4,ILU,ICODE,IERR)
      CALL BUFIOL (LDUM1,3,ILU,ICODE,IERR)
      CALL BUFIO4 (IDUM3,4,ILU,ICODE,IERR)
      CALL BUFIO4 (IDUM4,1,ILU,ICODE,IERR)
      CALL BUFIOL (LDUM2,5,ILU,ICODE,IERR)
      CALL BUFIOR (RDUM3,3,ILU,ICODE,IERR)
      CALL BUFIOL (LDUM3,1,ILU,ICODE,IERR)
      CALL BUFIOR (RDUM4,1,ILU,ICODE,IERR)
      CALL BUFIOR (RDUM5,4,ILU,ICODE,IERR)
      CALL BUFIOR (RDUM6,4,ILU,ICODE,IERR)
      CALL BUFIOL (LDUM4,1,ILU,ICODE,IERR)
      CALL BUFIOR (RDUM7,6,ILU,ICODE,IERR)
      CALL BUFIOL (LDUM5,3,ILU,ICODE,IERR)
      LRET=.FALSE.
      LGRCOM=.FALSE.
      IF ((L3DPL.NEQV.LDUM5(1)).OR.(L3DHOR.NEQV.LDUM5(2)).OR.
     &   (L3DVER.NEQV.LDUM5(3))) LRET=.TRUE.
      IF (R3DX1-RDUM7(1)+R3DY1-RDUM7(2)+R3DZ1-RDUM7(3)+R3DX2-RDUM7(4)+
     &    R3DY2-RDUM7(5)+R3DZ2-RDUM7(6).NE.0.0) LRET=.TRUE.
      IF (RX1-RDUM1(1)+RY1-RDUM1(2)+RX2-RDUM1(3)+RY2-RDUM1(4).NE.0.0)
     &LRET=.TRUE.
      IF (NX.NE.IDUM1(1).AND.NY.NE.IDUM1(2)) LRET=.TRUE.
      IF(RWX1-RDUM5(1)+RWY1-RDUM5(2)+RWX2-RDUM5(3)+RWY2-RDUM5(4).NE.0.0)
     &LRET=.TRUE.
      IF (ICODE.EQ.73.AND..NOT.LRET) LGRCOM=.TRUE.
      RETURN
 1000 FORMAT (' ***ERROR in GRTEMPIO: icode=',i3,' should be 73.')
      END
c
c -------------------------------------------------------------------------------
c
      SUBROUTINE RATEMPIO (ICODE,ILU)
c
c -------------------------------------------------------------------------------
c
C
C     Routine reads or writes the array RA to a file
C
      IMPLICIT NONE
      INTEGER(4) ICODE,ILU,IERR,J
      REAL(8) RA2
      INCLUDE 'GRID.INC'
      COMMON /BASEG/ RA2(200,200)      
      IF (NX.LE.0) RETURN
      DO 10 J=1,NY
      CALL BUFIOR (RA2(1,J),NX,ILU,ICODE,IERR)
  10  CONTINUE
      RETURN
      END
c
c -------------------------------------------------------------------------------
c
            SUBROUTINE DIFFGRID (RA,ISIZE)
c
c -------------------------------------------------------------------------------
c
C
C     Routine subtracts the grid in RA from the grid in RA2
c
c     This has been reversed from the original DOS version of GFLOW to
c     work with the GUI (8/17/04)
C      
      IMPLICIT NONE
      INTEGER(4) I,J,ISIZE
      REAL(8) RA,RA2
      INCLUDE 'GRID.INC'
      INCLUDE 'LUSYS.INC'
      COMMON /BASEG/ RA2(200,200)       ! for flownet or drawdown plots      
      DIMENSION RA(ISIZE,*)
C
      IF (.NOT.LGRCOM) THEN
       WRITE (ILUER,1000)
       RETURN
      ENDIF
      IF (NX.LE.0.OR.NX.GT.200) THEN
       WRITE (ILUER,2000) NX
       RETURN
      ENDIF
      IF (NY.LE.0.OR.NY.GT.200) THEN
       WRITE (ILUER,3000) NY
       RETURN
      ENDIF      
      DO 20 J=1,NY
      DO 10 I=1,NX
      RA(I,J)=RA2(I,J)-RA(I,J)
  10  CONTINUE
  20  CONTINUE      
C      
      CALL MINMAX (RA,ISIZE)      
      RETURN
 1000 FORMAT (' ***ERROR in DIFFGRID: no grid or incompatible grids.')      
 2000 FORMAT (' ***ERROR in DIFFGRID: illegal NX:',I5)
 3000 FORMAT (' ***ERROR in DIFFGRID: illegal NY:',I5) 
      END
c
c -------------------------------------------------------------------------------
c
      SUBROUTINE GRIDOUT (ICALL)
c
c -------------------------------------------------------------------------------
c
C
C     Routine writes the contents of common /GRID/ for debugging purposes
C
      IMPLICIT NONE
      INTEGER(4) ICALL,IERR,I
      CHARACTER(1) ADUM
      INCLUDE 'GRID.INC'
      INCLUDE 'LUSYS.INC'
C
      WRITE (ILUOUT,1000) RX1,RX2,RY1,RY2
      WRITE (ILUOUT,1005) RMAX,RMIN
      WRITE (ILUOUT,1010) NX,NY,RSCALE,IX1,IY1,IX2,IY2
      WRITE (ILUOUT,1020) (RXG(I),I=1,NX)
      WRITE (ILUOUT,2000) ICALL
      WRITE (ILUOUT,1030) (RYG(I),I=1,NY)
      WRITE (ILUOUT,1040) RSTEP,LPHI,LPSI,LHEAD,LFLOWNET,
     &                    NDOTX,NDOTY,NXMAX,NYMAX
      WRITE (ILUOUT,1041) ILAYER
      WRITE (ILUOUT,1042) IGRIDCODE
      WRITE (ILUOUT,1050) (LRESP(I),I=1,25)
      WRITE (ILUOUT,1060) LSCRN,LPRINT,LPLOT,LHPLOT,LGRID
      WRITE (ILUOUT,1065) RVFAC,RXSCALE,RYSCALE
      WRITE (ILUOUT,1066) RXMN,RYMN,RXMX,RYMX
      WRITE (ILUOUT,1067) LCURS,RDS0
      WRITE (ILUOUT,1070) LSURF,ISFCNT,ISFREC,CZBL0
      WRITE (ILUOUT,1075) RWX1,RWY1,RWX2,RWY2,LGRCOM,LDISCH
      WRITE (ILUOUT,1080) ATITLE
      WRITE (ILUOUT,2000) ICALL
      DO 10 I=1,NGRWIN
      WRITE(ILUOUT,1090)I,CZWIN(1,I),CZWIN(2,I),I,CZWIN(3,I),CZWIN(4,I)
  10  CONTINUE
      WRITE (ILUOUT,2000) ICALL
      WRITE (ILUOUT,1100) LFLOWNET,LDOTMAP,LSUBTRACT,LSLIDES,NLEVELS,
     &NTEXTS
       WRITE (ILUOUT,2000) ICALL
      RETURN
 1000 FORMAT (' GROUT: RX1,RX2,RY1,RY2 ',4(G11.4,2X))
 1005 FORMAT (' GROUT: RMAX,RMIN ',2(G11.4,2X))
 1010 FORMAT (' GROUT: NX,NY,RSCALE',2I5,G11.4,/
     &        ' GROUT: IX1,IY1,IX2,IY2 ',4I7)
 1020 FORMAT (20(' GROUT: RXG ',5(G11.4,2X)/))
 1030 FORMAT (20(' GROUT: RYG ',5(G11.4,2X)/))
 1040 FORMAT (' GROUT: RSTEP,LPHI,LPSI,LHEAD,LFLOWNET ',G11.4,2X,4L4,/
     &        ' GROUT: NDOTX,NDOTY,NXMAX,NYMAX ',4(I5,2X))
 1041 FORMAT (' ILAYER   =',I10)     
 1042 FORMAT (' IGRIDCODE=',I10)      
 1050 FORMAT (5(' GROUT: LRESP ',5(L3,2X)/))
 1060 FORMAT (' GROUT: LSCRN,LPRINT,LPLOT,LHPLOT,LGRID ',5(L3,2X))
 1065 FORMAT (' GROUT: RVFAC,RXSCALE,RYSCALE ',3(G11.4))
 1066 FORMAT (' GROUT: RXMN,RYMN,RXMX,RYMX ',4G11.4)
 1067 FORMAT (' GROUT: LCURS,RDS0 ',L3,3X,G11.4)
 1070 FORMAT (' GROUT: LSURF,ISFCNT,ISFREC,CZBL0 ',L3,2I5,2G11.4)
 1075 FORMAT (' GROUT: RWX1,RWY1,RWX2,RWY2 ',4G11.4,/,
     &        '        LGRCOM,LDISCH ',2L3)
 1080 FORMAT (' GROUT: ATITLE ',A16)
 1090 FORMAT (' GROUT: CZWIN (1-2,',I3,') ',4G11.4,/,
     &        '        CZWIN (3-4,',I3,') ',4G11.4)
 1100 FORMAT (' GROUT: LFLOWNET,LDOTMAP,LSUBTRACT,LSLIDES ',4L3,/,
     &        '        NLEVELS,NTEXTS ',3I5)     
 2000 FORMAT (' GRIDOUT: ICALL=',I3,' press <Enter> to continue.')
 3000 FORMAT (A1)
      END

