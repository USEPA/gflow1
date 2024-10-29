C     Last change:  HMH  10 May 2010   10:04 am
c       This file contains the following routines and functions:
c
C	BLOCK DATA W3DAT    initializes ppwell common blocks
C	SUBROUTINE W3IN     handles input of ppwells as defined in the .dat file
C
c -------------------
c w3mat.for contains:
c -------------------
c
C    SUBROUTINE W3CZC         generates control point vector for ppwells
C    SUBROUTINE W3MAT         generates matrix coefficients for ppwells
C    SUBROUTINE W3KNO         generates known vector for ppwells
C    SUBROUTINE W3SUB         substitutes solution vector into strengths parameters for ppwells
C    REAL FUNCTION RFW3MATC   coefficient for ppwell strength parameters   
c
c -------------------
c w3fun.for contains:
c -------------------
c
C     REAL FUNCTION RFW3PT        Specific discharge potential contribution of all part. penetrating wells
C     COMPLEX FUNCTION CFW3OM     Calculates ppwell contribution to potential in 2D zone
C     COMPLEX FUNCTION CFW3OMCOF  Coefficient function for ppwell in 2D zone
C     REAL FUNCTION RFLSF         coefficient function F
C     REAL FUNCTION RFLSG         coefficient function G
C     REAL FUNCTION RFLSE         coefficient function E
C     SUBROUTINE W3GENC           generation of common constants of F,G,E pot. coef. fcns.
C     REAL FUNCTION RFW3CP        coefficient function  (potential) for node INODE of line sink string (well) IW
C     REAL FUNCTION RFW3IP        coefficient due to the F, G, or E functions
C     SUBROUTINE W3QCO            coefficient function for specific discharge of the well.
C     SUBROUTINE W3QIM            contributions to QTAU and QZ due to the F, G, and E functions
C     SUBROUTINE W3QI2D           specific discharge vector due to all ppwells
C          ENTRY W3QI3D
C     SUBROUTINE W3QTZF           returns the radial and vertical specific discharge component for F function
C     SUBROUTINE W3QTZG           returns the radial and vertical specific discharge component for G function
C     SUBROUTINE W3QTZE           returns the radial and vertical specific discharge component for E function
c
c ---------------------
c w3nflow.for contains:
c ---------------------
c
C     REAL FUNCTION RFNFW3    returns the flow across a line due to all wells
C     REAL FUNCTION RFNFW3CO  coefficient function for the flow across a line due to well i.
c
c ------------------
c w3io.for contains:
c ------------------
c
C     SUBROUTINE W3EXTRACT        print data for partially penetrating wells (NOT USED).
C     SUBROUTINE w3extracthead    writes data for head specified wells to .xtr file
C     SUBROUTINE w3extractdisch   writes data for discharge specified wells to .xtr file
C     SUBROUTINE W3DATIO          rites all well data to a ".dat" file
C     SUBROUTINE W3IO             reads or writes contents of WLCOM.INC from or to .sol file
C     SUBROUTINE W3OUT            writes contents of WLCOM formatted to ILUOUT for debugging purposes
c
c ---------------------
c w3check.for contains:
c ---------------------
c
C     SUBROUTINE W3ERROR   maximum error in the boundary conditions specified at ppwells
c
c --------------------
c w3near.for contains:
c --------------------
c
C     subroutine W3NEAR   Routine sets the logical L3DEND=.true. when trace within well radius.
c
c --------------------
c w3serv.for contains:
c --------------------
c
C     Logical FUNCTION lw3info        routine returns information for the 2D well
C     LOGICAL function lw3discharge   TRUE if there are discharge specified partially penetrating wells
C     LOGICAL function lw3head        TRUE if there are discharge specified partially penetrating wells
c
c ---------------------
c w3extra.for contains:
c ---------------------
c
C     SUBROUTINE W3PLOT        plots all partially penetrating wells
C     SUBROUTINE W3CUR         facilitates graphical display or modification of well attributes
C     SUBROUTINE W3CHEC        check of boundary conditions
C
C
C ------------------------------------------------------------------------------------------------
C
      BLOCK DATA W3DAT
C
C ------------------------------------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER(4) NIMAG
      LOGICAL L2D,L3D,L3DH,LW3MAT,LDRT,LEQPT,LEQPB
      REAL(8) RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &        RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU,
     &        RQZF,RQZG,RQZE,RA,RB,RC,RD,RK,GENOLD,
     &        RXOLD,RF,RG,RE
      INCLUDE 'w3com.inc'
      COMMON/W3PASS/RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &              RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU
      COMMON/FAR3D/ L2D,L3D,L3DH
      COMMON/W3MATT/ RQZF,RQZG,RQZE,LW3MAT
      COMMON/GENC/RA,RB,RC,RD,RK,GENOLD(4),LDRT
      COMMON/W3OLD/RXOLD(3,NW3WMX),RF,RG,RE      
      COMMON/IMAG/ NIMAG,LEQPT,LEQPB
      SAVE
      DATA RPI1,RPI2,RPI4,RPI8/3.141592654,6.283185307,
     & 12.56637061,25.13274123/
      DATA RW3RAD,RW3HED,RW3Q /NW3WMX*0.0,NW3WMX*0.0,NW3WMX*0.0/
      DATA AW3LAB /NW3WMX*'                '/
      DATA RW3ST /NW3SMX*0.0,NW3SMX*0.0,NW3SMX*0.0/
      DATA RW3CP /NW3SMX*0.0,NW3SMX*0.0,NW3SMX*0.0/
      DATA RW3S,RW3L /NW3SMX*0.0,NW3SMX*0.0/
      DATA IPNT,NW3,NW3RH,NLIMAG /NW3WMX*0,0,10,10/
      DATA LW3Q,LW3DRT,LW3BOT,LW3TOP /NW3WMX*.FALSE.,NW3WMX*.FALSE.,
     &                                NW3WMX*.FALSE.,NW3WMX*.FALSE./
      DATA RZT,RZB,RHO2,RT,RH,RH2,R2H,RTAU2,RTAU/9*0.0/
      DATA RW3ALP/0.04/
      DATA LW3MAT /.FALSE./
      DATA GENOLD,LDRT /4*1.0E30,.FALSE./
      DATA RXOLD/NW3WMX*1.E30,NW3WMX*1.E30,NW3WMX*1.E30/ 
C     NIMAG controls # of images (see W3QCO and RFW3CP)
C     if LEQPT aquifer top equipotential and aquifer base no-flow boundary
C     if LEQPB aquifer base equipotential and aquifer top no-flow boundary
      DATA NIMAG,LEQPT,LEQPB /3,.FALSE.,.FALSE./
      DATA rw3d0 /nw3wmx*0.0,nw3wmx*0.0,nw3wmx*0.0/
      DATA RZONEMLTP /3.0/
      END
C
C ------------------------------------------------------------------------------------------------
C
      SUBROUTINE W3IN (LSOL,RA,IRA,RSCR)
C
C ------------------------------------------------------------------------------------------------
C
C
C     Input routine for partially penetrating wells.
C
      IMPLICIT NONE
      INTEGER(4) IRA,INODE,JUMP,IZPT,I,IST,IEN,IW,NCP,IENCP,II
      LOGICAL LSOL,LBAD,LW3QLOC,LOUT1,LOUT2
      REAL(8) RA,RSCR,RZPT,RW3R0,RTOP,RFTOP,RBASE,RFBASE,RFPOR,RFPERM,
     &        RDUM1,RVAR,RDUM2,RDUM3,RDUM4,RDUM5,RDUM6,RXWELL,RYWELL,
     &        RDUM,RLSUM,RL,RSINIT,RH,RZST,RZEN,RDL
      COMPLEX(8) CZ,CDUM
      CHARACTER(1) AWORD(28)
      INCLUDE 'w3com.inc'
      INCLUDE 'match.inc'
      INCLUDE 'com3d.inc'
      INCLUDE 'lusys.inc'
      DIMENSION RZPT(NW3SMX),RA(IRA,*),RSCR(*)
      DATA INODE,RW3R0 /1,1.0/
      DATA AWORD /
     .            '?',' ',
     .            'H','E','A','D',' ',
     .            'D','I','S','C',' ',     
     .            'C','U','R','S',' ',
     .            'R','A','D','I',' ',
     .            'Q','U','I','T',' ',
     .            ATERM/
      LERROR=.FALSE.
      LMISS=.FALSE.     
C      CALL CLEARSCREEN    ! NOT AVAILABLE IN BATCH MODE
   10 IF (LERROR.OR.LMISS) WRITE (ILUER,3000) ALINE2
      if (lucon) WRITE (ILUME,1000) NW3WMX,NW3SMX,RW3R0
      LERROR=.FALSE.
      LMISS=.FALSE.
      CALL INLINE
   11 CALL MATCH (AWORD,1,JUMP,LBAD)
C      CALL CLEARSCREEN    ! NOT AVAILABLE IN BATCH MODE
      IF (.NOT.LBAD) GOTO 15
      GOTO (10,13),JUMP
   13 WRITE(ILUER,4000) ALINE2
      LERROR=.FALSE.
      GOTO 10
   15 GOTO (100,200,300,400,500,600),JUMP
C
  100 AFILE='W3HLP.HLP      '
C      CALL HELP                                 ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
C
C   Head Specified Wells
C
  200 LW3QLOC=.FALSE.
      GOTO 350
C
C     Discharge Specified Wells
C
  300 LW3QLOC=.TRUE.
      GOTO 350
  340 WRITE (ILUER,6000)
  350 IF (.NOT.LW3QLOC.and.lucon) WRITE(ILUME,5000)       ! enter well data
      IF (LW3QLOC.and.lucon) WRITE (ILUME,7000)
      CALL INLINE
      RDUM1=RVAR(1)
      IF (LERROR) THEN
      LERROR=.FALSE.
      GOTO 11
      ENDIF
      RDUM2=RVAR(2)
      CZ=CMPLX(rdum1,rdum2)
      RTOP=RFTOP(CZ)
      RBASE=RFBASE(CZ)
      R3DH=RTOP-RBASE
      R3DPOR=RFPOR(CZ)
      R3DK=RFPERM(CZ)
c      write (iluer,1001) cz,rbase
 1001 format (' w3mod1: cz,rbase ',3(E14.7))
      RDUM3=RVAR(3)  ! z-value well screen
      RDUM4=RVAR(4)  ! z-value well screen
      RDUM5=RVAR(5)
      IF (LERROR) GOTO 340
      IF (RDUM3.EQ.RDUM4) THEN      ! zero length of well screen
        WRITE (ILUER,4700)
        GOTO 340
      ENDIF
      LOUT1=RDUM3.LT.RBASE.OR.RDUM3.GT.RTOP
      LOUT2=RDUM4.LT.RBASE.OR.RDUM4.GT.RTOP
      IF (LOUT1.OR.LOUT2) THEN      ! one or both points outside aquifer
        WRITE (ILUER,6500) RDUM3,RDUM4,RBASE,RTOP
        GOTO 340
      ENDIF
      RDUM6=RVAR(6)
      IF (LERROR) RDUM6=RW3R0   ! radius of well
      CALL GETFN (7)
      NW3=NW3+1
      IF (.NOT.LERROR.OR..NOT.LMISS) AW3LAB(NW3)=AFILE  ! label of well
      LERROR=.FALSE.
      LMISS=.FALSE.
      LSOL=.FALSE.
      RXWELL=RDUM1
      RYWELL=RDUM2
      IF (.NOT.LW3QLOC) RW3HED(NW3)=RDUM5   ! set head
      IF (LW3QLOC) RW3Q(NW3)=RDUM5          ! set discharge
      IF (RDUM6.LT.1.0E-10) THEN
        WRITE (ILUER,3500) RDUM6,RW3R0
        RW3RAD(NW3)=RW3R0    ! if well radius too small, set to default
      ELSE
        RW3RAD(NW3)=RDUM6
      ENDIF
      LW3Q(NW3)=LW3QLOC    ! set flag for discharge specified
      IPNT(NW3)=INODE      ! point to first node on well
  450 IZPT=2               ! hard wire a single line-sink at well center
      IF (RDUM4.LT.RDUM3) THEN
      RZPT(1)=RDUM4
      RZPT(2)=RDUM3
      ELSE                 ! make sure bottom of well screen is stored first
      RZPT(1)=RDUM3
      RZPT(2)=RDUM4
      ENDIF
c
c       This logic allowed for wells consisting of strings of line-sinks at the well axis.
c       Currently not in use. Use of 1 line-sink along well axis only.
c
C     WRITE(ILUME,8000)       ! enter z points along well
C      IZPT=0
C  460 WRITE (ILUME,8010)
C      CALL INLINE
C      IF (ALINE(1).EQ.'Q'.OR.ALINE(1).EQ.'q') GOTO 465
C      RDUM=RVAR(1)
C      IF(LERROR)THEN
C      WRITE(ILUER,6000)
C      GOTO 460
C      ENDIF
C      IF (RDUM.LE.RZPT(IZPT).AND.IZPT.GT.0) THEN ! check for proper order
C      WRITE (ILUER,8100)
C      GOTO 450
C      ENDIF
C      IF (RDUM.LT.RBASE.OR.RDUM.GT.RTOP) THEN
C      WRITE (ILUER,6500)
C      GOTO 460
C      ENDIF
C      IZPT=IZPT+1
C      IF (IZPT+INODE.GT.NW3SMX) THEN
C      WRITE (ILUER,2000)
C      GOTO 450
C      ENDIF
C      RZPT(IZPT)=RDUM      
C      GOTO 460
C  465 IF (IZPT.LT.2) THEN     ! need at least two points, skip
C      WRITE (ILUER,4600) IZPT
C      IF (NW3.EQ.0) GOTO 10
C      RW3HED(NW3)=0.0
C      RW3Q(NW3)=0.0
C      RW3RAD(NW3)=0.0
C      LW3Q(NW3)=.FALSE.
C      NW3=NW3-1
C      GOTO 10
C      ENDIF      
      DO 469 I=1,IZPT         ! create screen sections (line sinks)
      RW3ST(1,INODE)=RXWELL   ! currently izpt=2, thus only one screen section
      RW3ST(2,INODE)=RYWELL
      RW3ST(3,INODE)=RZPT(I)       
      INODE=INODE+1
  469 CONTINUE    
  470 IPNT(NW3+1)=INODE
      IST=IPNT(NW3)      ! not used, statement may be deleted?
      IEN=IPNT(NW3+1)-1  ! not used, statement may be deleted?
      CZ=CMPLX(RXWELL,RYWELL)
      CDUM=CZ+CMPLX(0.0,RW3RAD(NW3)) ! setting window extents in GFLOW1
      CALL PLWIND (CDUM)
      CDUM=CZ-CMPLX(0.0,RW3RAD(NW3))
      CALL PLWIND (CDUM)
      CDUM=CZ+RW3RAD(NW3)
      CALL PLWIND (CDUM)
      CDUM=CZ-RW3RAD(NW3)
      CALL PLWIND (CDUM)
      LW3BOT(NW3)=RW3ST(3,IST).EQ.RFBASE(CZ)
      LW3TOP(NW3)=RW3ST(3,IEN-1).EQ.RFTOP(CZ)
      IF (LW3BOT(NW3).AND.LW3TOP(NW3)) GOTO 10
c
c       Note: double root flag set true for all cases, except when well screen extends
c             from aquifer bottom to aquifer top. In this case the extra storage location
c             for a third collocation point is not needed.
c
        
      LW3DRT(NW3)=.TRUE.
      INODE=INODE+1
      IPNT(NW3+1)=INODE
      RW3ST(1,INODE-1)=RXWELL
      RW3ST(2,INODE-1)=RYWELL
      IF (LW3QLOC) GOTO 300
      GOTO 200
C
C     Cursor
C      
  400 IF (.NOT.LSOL) THEN
      WRITE (ILUER,4500)
C      CALL TONE                                ! NOT AVAILABLE IN BATCH MODE
      ENDIF
C      CALL W3CUR (RA,IRA,RSCR,LSOL)            ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
C
C     Radius
C
  500 RDUM=RVAR(2)
      IF (LERROR) GOTO 10
      IF (RDUM.LT.1.0E-10) THEN
        WRITE (ILUER,3500) RDUM,RW3R0
        GOTO 10
      ENDIF
      RW3R0=RDUM
      GOTO 10
C
C     Quit                (calculate length parameters belonging to strength parameters)
C
  600 IF (NW3.EQ.0) GOTO 825
      DO 820 IW=1,NW3
      IST=IPNT(IW)
      IEN=IPNT(IW+1)-1
      NCP=IEN-IST+1
      IENCP=IEN
      IF (LW3DRT(IW)) IEN=IEN-1
      RLSUM=0.0
      DO 810 II=IST,IEN       ! calc. length parameters for linear strength
      RL=0.0
      IF (II.GT.IST) RL=RL+RW3ST(3,II)-RW3ST(3,II-1)
      IF (II.LT.IEN) RL=RL+RW3ST(3,II+1)-RW3ST(3,II)
      RW3L(II)=RL/2.0
      RLSUM=RLSUM+RW3L(II)
  810 CONTINUE
      IF (LW3Q(IW)) THEN       ! set initial strength to equal total discharge
            RSINIT=RW3Q(IW)/RLSUM
            DO 812 II=IST,IEN
            RW3S(II)=RSINIT
  812       CONTINUE      
      ENDIF
      IF (LW3DRT(IW)) THEN    ! calc. length parameter for singular strength
          RL=RW3ST(3,IEN)-RW3ST(3,IST)
          RH=RL/2.0
          RW3L(IEN+1)=RH*3.141592654
      ENDIF
C
C   Generate control points
C
      RXWELL=RW3ST(1,IST)
      RYWELL=RW3ST(2,IST)
      RZST=RW3ST(3,IST)+RW3ALP*(RW3ST(3,IST+1)-RW3ST(3,IST))
      RZEN=RW3ST(3,IEN)-RW3ALP*(RW3ST(3,IEN)-RW3ST(3,IEN-1))
      RDL=(RZEN-RZST)/(NCP-1)
      RW3CP(1,IST)=RXWELL+RW3RAD(IW)
      RW3CP(2,IST)=RYWELL
      RW3CP(3,IST)=RZST
      IST=IST+1
      DO 815 I=IST,IENCP
      RW3CP(1,I)=RXWELL+RW3RAD(IW)
      RW3CP(2,I)=RYWELL
      RW3CP(3,I)=RW3CP(3,I-1)+RDL
  815 CONTINUE
  820 CONTINUE
C  
  825 RETURN 
C 
 1000 FORMAT (' -------- PARTIALLY PENETRATING WELL module ----------'/
     &' Maximum number of wells:   ',I4,/
     &' Maximum number of sections:',I4,' (total for all wells)'/
     &' <F1> = Help'/' HEAD'/' DISCHARGE'/' RADIUS ',G11.4,/
     &' <Esc> or QUIT '/' >')
 2000 FORMAT (' ***ERROR: too many well screen sections:'/
     &' redo, entering fewer z-values.')
 3000 FORMAT (' ***ILLEGAL OR MISSING PARAMETER(S) in ppwell module:',/,
     &        ' ',80A1)
 3500 FORMAT (' ***ERROR in well module: radius of ',G11.4,
     &' is too small.',/,' Current default value used: ',G14.7,/)     
 4000 FORMAT (' ***ILLEGAL COMMAND in ppwell module:',/,' ',80A1)
 4500 FORMAT (' ***WARNING: no valid solution, data likely in error.')
 4600 FORMAT (' ***ERROR: need at least two points to define well',
     &' screen, only ',I2,' found.'/' WELL NOT ENTERED !'/)
 4700 FORMAT (' ***ERROR: well screen length cannot be zero,',
     &' well not entered.'/)
 5000 FORMAT (' X   Y   Z1  Z2  HEAD   [RADIUS] [LABEL]'/' >')
 6000 FORMAT (' ***ERROR: repeat input',/)
 6500 FORMAT (' ***ERROR: well screen outside aquifer!',/,
     &' end points of well screen: ',G11.4,1X,G11.4,/,
     &' aquifer bottom=',G11.4,' aquifer top=',G11.4,/,
     &' Well has not been entered!',/)
 7000 FORMAT (' X   Y   Z1  Z2  DISCHARGE  [RADIUS] [LABEL]'/' >')
 8000 FORMAT (
     &' Enter z-values of well screen sections, starting at the bottom'/
     &' Type RETURN to terminate input.'/)
 8010 FORMAT ('+>') 
 8100 FORMAT (' *** ERROR: z-value out of order, give new value.')
      END
