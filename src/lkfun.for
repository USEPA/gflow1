C     Last change:  HMH  17 Nov 2008    4:52 am
c
c
c   lk_cell_construct_actual  return the four corner points for a grid cell at i,j
c   cflk_omega                complex function that calls ENTRY cflk_omega_sub in LKIN
c                             which calls:
c   cflk_omega_actual         to return the complex potential at cz
c   cflk_omega_par            returns the complex potential due to a line-doublet
c   lkqi                      call ENTRY cflk_W_sub in LKIN which calls:
c   cflk_W_sub_actual         to return Q_i and q_3 for the grid
c   cflk_W_par                returns the function W (negative derivative of complex potential) due to a line-doublet
c
c   cflk_subomega             complex function that call entry cflk_subomega_sub in LKIN
c                             which calls:
c   cflk_subomega_actual      to return complex potential at CZ for sub-cells
c   cflk_subomegaMODFLOWcell  called by cflk_subomega_actual returns all sub-cell contributions for a single MODFLOW cell
c   lksubqi                   calls entry cflk_subW_sub in LKIN which calls:
c   lk_subW_sub_actual        to return Q_i and q_3 for sub-cells
c   cflk_omega_coefficient    coefficient function for a MODFLOW cell or sub_cell
c
c --------------------------------------------------------------
c
      subroutine lk_cell_construct_actual(i,j,rdx,rdy,cz0,cz1,cz2,cz3,
     &                                    rlkdeltax,rlkdeltay)
c
c     Routine is called by "lk_construct" in LKIN.
c
c     Input:   i and j
c
c     Output;  cz0, cz1, cz2, and cz3
c
c     rx0,ryo is lower left of cell =cz0
c     rx1,ry1 is upper right of cell =cz2
c
c     cell structure:
c                               Origin-----> columns
c     cz3-----------cz2         ||
c      |             |          ||
c      |             |          \/
c      |             |
c      |             |         rows
c     cz0-----------cz1
c
c
      implicit none
      INTEGER i,j
      REAL(8) rx0,ry0,rx1,ry1,rdx,rdy,
     &         rlkdeltax,rlkdeltay
      COMPLEX(8) cz0,cz1,cz2,cz3
      DIMENSION rlkdeltax(*),rlkdeltay(*)
      INCLUDE 'lkcom.inc'
      include 'lusys.inc'
c
      if (lkleakage) then
      rx1=rlkx0+SUM(rlkdeltax(1:j))
      ry0=rlky0-SUM(rlkdeltay(1:i))
      rx0=rx1-rlkdeltax(j)
      ry1=ry0+rlkdeltay(i)
      cz0=CMPLX(rx0,ry0)
      cz1=CMPLX(rx1,ry0)
      cz2=CMPLX(rx1,ry1)
      cz3=CMPLX(rx0,ry1)
      rdx=rlkdeltax(j)
      rdy=rlkdeltay(i)
      endif
      end subroutine
c
c ---------------------------------------------------------------
c
      complex(8) function cflk_omega (cz)
c
c     Return potential due to all leakage elements in MODFLOW grid
c
      implicit none
      COMPLEX(8) cz,cfunc
      include 'lkcom.inc'
      cflk_omega=CMPLX(0.0,0.0)
      if (lkleakage) then
       call cflk_omega_sub (cz,cfunc)
       cflk_omega=cfunc
      end if
      return
      end
c
c ----------------------------------------------------------------
c
      subroutine cflk_omega_actual (cz,cfunc,nrow,ncol,
     &           rlkdeltax,rlkdeltay,rlkleakage,rlkrecharge,rlkconst)
c
c     This routine is called by entry "clk_omega" in LKIN.
c     Includes the contributions of the MODFLOW cells, but not yet of the sub-grids.
c     Note: The complex potential consists of the true "omega" for the line-sinks and
c           line-doublets plus the Poisson term inside the element added to the real part of "omega."
c
      implicit none
      INTEGER nrow,ncol,i,j,i0,j0
      REAL(8) rlkdeltax,rlkdeltay,rlkleakage,rlkrecharge,rlkconst,rpot,
     &        rx,ry,rdx,rdy,rx0,ry0,rx1,ry1,rx2,ry2,rs,rsm,rsp
      COMPLEX(8) cz,cfunc,cz1,cz2,comls,cflk_omega_par,cdum
      DIMENSION rlkdeltax(ncol),rlkdeltay(nrow),rlkleakage(nrow,ncol),
     &          rlkrecharge(nrow,ncol),rlkconst(nrow,ncol+1)
      include 'lkcom.inc'
      include 'lusys.inc'
c
c                   >> include all line-sinks (along vertical grid lines)
c
c      write (iluer,1100)
c 1100 format (//,' Entering routine "cflk_omega_actual."',/)
      cfunc=CMPLX(0.0,0.0)
      rx1=rlkx0
      ry1=rlky0
      do j=1,ncol+1
      if (j.ne.ncol+1) rdx=rlkdeltax(j)
      do i=1,nrow
       cz1=CMPLX(rx1,ry1)
       rdy=rlkdeltay(i)
       ry2=ry1-rdy
       cz2=CMPLX(rx1,ry2)
       if (j.eq.1) then
         rsm=0.0
         rsp=(rlkleakage(i,1)+rlkrecharge(i,1))*rdx*0.5
       else if (j.eq.ncol+1) then
         rsm=(rlkleakage(i,ncol)+rlkrecharge(i,ncol))*rdx*0.5 ! rdx still rlkdeltax(ncol)
         rsp=0.0
       else
         rsm=(rlkleakage(i,j-1)+rlkrecharge(i,j-1))*rlkdeltax(j-1)*0.5
         rsp=(rlkleakage(i,j)+rlkrecharge(i,j))*rdx*0.5
       end if
       rs=rsm+rsp
       cdum=comls(cz,cz1,cz2)+rlkconst(i,j)
       cfunc=cfunc+rs*cdum ! note: COMLS must be given a far-field expansion
c       write (iluer,1001) i,j,rs,cz,cz1,cz2,cdum,cfunc
c 1001 format (' cflk_omega_actual1: i,j,rs / cz,cz1,cz2 / cdum, cfunc ',
c     &/,2(i5),1x,d14.7,/,6(d14.7),/,4(d14.7))
       ry1=ry2
      end do
      ry1=rlky0
      rx1=rx1+rdx
      enddo
c                   >>  include all line-doublets (along horizontal grid lines)
      rx1=rlkx0
      ry1=rlky0
      do i=1,nrow+1
      if (i.ne.nrow+1) rdy=rlkdeltay(i)
      do j=1,ncol
       cz1=CMPLX(rx1,ry1)
       rdx=rlkdeltax(j)
       rx2=rx1+rdx
       cz2=CMPLX(rx2,ry1)
       IF (i.eq.1) then
         rsm=0.0
         rsp=-(rlkleakage(1,j)+rlkrecharge(1,j))*rdx*rdx*0.125
       else if (i.eq.nrow+1) then
         rsm=-(rlkleakage(nrow,j)+rlkrecharge(nrow,j))*rdx*rdx*0.125
         rsp=0.0
       else
          rsm=-(rlkleakage(i-1,j)+rlkrecharge(i-1,j))*rdx*rdx*0.125
          rsp=-(rlkleakage(i,j)+rlkrecharge(i,j))*rdx*rdx*0.125
       end if
       rs=rsp-rsm ! note, due to MODFLOW row numbering when i increases ry1 and ry2 decrease
       cdum=cflk_omega_par(cz,cz1,cz2)  ! this leads to different sign than in sub-cell routines
       cfunc=cfunc+rs*cdum
c       write (iluer,1002) i,j,rs,cz,cz1,cz2,cdum,cfunc
c 1002 format (' cflk_omega_actual2: i,j,rs / cz,cz1,cz2 / cdum, cfunc ',
c     &/,2(i5),1x,d14.7,/,6(d14.7),/,4(d14.7))
       rx1=rx2
      end do
      rx1=rlkx0
      ry1=ry1-rdy
      end do
c
c                  >> include the Poisson term
c
      call lkfindcell (cz,i0,j0,rx0,ry0)
      if (i0.gt.0.and.j0.gt.0) then  ! CZ is inside cell i0,j0
       rs=(rlkleakage(i0,j0)+rlkrecharge(i0,j0))  ! add Poisson term for cell
       rdx=0.5*rlkdeltax(j0) ! rdx is now half width of cell
       rx=REAL(cz)
       rpot=0.5*rs*(rx-rx0-rdx)*(rx-rx0+rdx)
       cfunc=cfunc+CMPLX(rpot,0.0)
c       write (iluer,1003) i0,j0,rs,rdx,rpot,cfunc
c 1003 format (' cflk_omega_actual3: i0,j0,rs,rdx,rpot,cfunc ',/,
c     & 2(i3,1x),1x,5(d14.7,1x))
      end if
c                        !!!! no sub-grid contribution yet
      return
      end subroutine
c
c ------------------------------------------------------------------------------------------------------------
c
      COMPLEX(8) FUNCTION cflk_omega_par (cz,cz1,cz2)
c
c ------------------------------------------------------------------------------------------------------------
c
C
C     POTENTIAL due to the quadratic term for line doublet.
C     Repeat of the function CFDBS(IEL), but not using DBPREP
C     Function is called in clk_omega_actual in file LKFUN.FOR
C
      IMPLICIT NONE
      INTEGER(4) I,ndbtrm
      LOGICAL lfar
      REAL(8) DCOF,dpi
      COMPLEX(8) cz,cz1,cz2,cdum1,cdum2,CBZ,cfbigz,cdi
      INCLUDE 'LUSYS.INC'
      DIMENSION DCOF(21)
      DATA ndbtrm,cdi,dpi /17,(0.0,1.0),.3141592653589793D+1/
      DATA DCOF /0.3333333333333333D0,0.0D0,
     &           0.0666666666666666D0,0.0D0,
     &           0.0285714285714285D0,0.0D0,
     &           0.0158730158730158D0,0.0D0,
     &           0.0101010101010101D0,0.0D0,
     &           0.0069930069930069D0,0.0D0,
     &           0.0051282051282051D0,0.0D0,
     &           0.0039215686274500D0,0.0D0,
     &           0.0030959752321981D0,0.0D0,
     &           0.0025062656641604D0,0.0D0,
     &           0.0020703933747412D0/
C
            cbz=cfbigz(cz,cz1,cz2)
            cflk_omega_par=CMPLX(0.0,0.0)
            lfar=ABS(cbz).gt.6.0 ! switch to farfield expansion at 3 times line-doublet length
      IF (.NOT.lfar) THEN
            if (abs(cbz-1.0).lt.0.001) cbz=cbz-0.001 ! shift onto line-doublet away from vertex
            if (abs(cbz+1.0).lt.0.001) cbz=cbz+0.001 ! shift onto line-doublet away from vertex
            cdum1=cbz*cbz-1.0D0
            cdum2=(cbz-1.0)/(cbz+1.0)
            cdum2=LOG((cbz-1.0)/(cbz+1.0))
            cflk_omega_par=(cbz+0.5D0*cdum1*cdum2)*CDI/DPI
      ELSE
            cdum1=1.0/cbz
            cflk_omega_par=DCOF(1)*cdum1
            cdum2=cdum1*cdum1
            DO 10 I=3,ndbtrm,2
            cdum1=cdum1*cdum2
            cflk_omega_par=cflk_omega_par+DCOF(I)*cdum2
  10        CONTINUE
            cflk_omega_par=cflk_omega_par*2.0D0*CDI/DPI
      ENDIF
      RETURN
      END
c
c ---------------------------------------------------------------
c
      subroutine lkqi (cz,rqi)
c
c     Return discharge vector due to all leakage elements in MODFLOW grid
c     Note: rqi(1) and rqi(2) are the discharge vector components, while
c           rqi(3) is the vertical component of the specific discharge.
c
      implicit none
      REAL(8) rqi,rq3
      COMPLEX(8) cz,cfunc
      DIMENSION rqi(3)
      include 'lkcom.inc'
      if (lkleakage) then
       rq3=0.0
       call cflk_W_sub (cz,cfunc,rq3)
       rqi(1)=rqi(1)+REAL(cfunc)
       rqi(2)=rqi(2)-AIMAG(cfunc)
       rqi(3)=rqi(3)+rq3
      endif
      return
      end
c
c ----------------------------------------------------------------
c
      subroutine cflk_W_actual (cz,cfunc,rq3,nrow,ncol,
     &           rlkdeltax,rlkdeltay,rlkleakage,rlkrecharge)
c
c     This routine is called by entry "clk_W_sub" in LKIN.
c     Includes the contributions of the MODFLOW cells, but not yet of the sub-grids.
c
      implicit none
      INTEGER nrow,ncol,i,j,i0,j0
      REAL(8) rlkdeltax,rlkdeltay,rlkleakage,rlkrecharge,rqx,
     &        rx,ry,rdx,rdy,rx0,ry0,rx1,ry1,rx2,ry2,rs,
     &        rq3,rfbase,rfhght,rbloc,rhght,rsm,rsp
      COMPLEX(8) cz,cfunc,cz1,cz2,cdls,cflk_W_par,cdum
      DIMENSION rlkdeltax(ncol),rlkdeltay(nrow),rlkleakage(nrow,ncol),
     &          rlkrecharge(nrow,ncol)
      include 'lkcom.inc'
      include 'com3d.inc'
      include 'lusys.inc'
c
c                   >> include all line-sinks (along vertical grid lines)
c
c
c      write (iluer,1100)
c 1100 format (//,' Entering routine "cflk_W_actual."',/)
      cfunc=CMPLX(0.0,0.0)
      rx1=rlkx0
      ry1=rlky0
      do j=1,ncol+1
      if (j.ne.ncol+1) rdx=rlkdeltax(j)
      do i=1,nrow
       cz1=CMPLX(rx1,ry1)
       rdy=rlkdeltay(i)
       ry2=ry1-rdy
       cz2=CMPLX(rx1,ry2)
       if (j.eq.1) then
         rsm=0.0
         rsp=(rlkleakage(i,1)+rlkrecharge(i,1))*rdx*0.5
       else if (j.eq.ncol+1) then
         rsm=(rlkleakage(i,ncol)+rlkrecharge(i,ncol))*rdx*0.5
         rsp=0.0
       else
         rsm=(rlkleakage(i,j-1)+rlkrecharge(i,j-1))*rlkdeltax(j-1)*0.5
         rsp=(rlkleakage(i,j)+rlkrecharge(i,j))*rdx*0.5
       end if
       rs=rsm+rsp
       cdum=-cdls(cz,cz1,cz2) ! we need the negative derivative of the complex potential
       cfunc=cfunc+rs*cdum ! note: COMLS must be given a far-field expansion
c       write (iluer,1001) i,j,rs,cz,cz1,cz2,cdum,cfunc
c 1001 format (' cflk_W_actual1: i,j,rs / cz,cz1,cz2 / cdum, cfunc ',
c     &/,2(i5),1x,d14.7,/,6(d14.7),/,4(d14.7))
       ry1=ry2
      end do
      ry1=rlky0
      rx1=rx1+rdx
      enddo
c                   >>  include all line-doublets (along horizontal grid lines)
      rx1=rlkx0
      ry1=rlky0
      do i=1,nrow+1
      if (i.ne.nrow+1) rdy=rlkdeltay(i)
      do j=1,ncol
       cz1=CMPLX(rx1,ry1)
       rdx=rlkdeltax(j)
       rx2=rx1+rdx
       cz2=CMPLX(rx2,ry1)
       IF (i.eq.1) then
         rsm=0.0
         rsp=-(rlkleakage(1,j)+rlkrecharge(1,j))*rdx*rdx*0.125
       else if (i.eq.nrow+1) then
         rsm=-(rlkleakage(nrow,j)+rlkrecharge(nrow,j))*rdx*rdx*0.125
         rsp=0.0
       else
          rsm=-(rlkleakage(i-1,j)+rlkrecharge(i-1,j))*rdx*rdx*0.125
          rsp=-(rlkleakage(i,j)+rlkrecharge(i,j))*rdx*rdx*0.125
       end if
       rs=rsp-rsm  ! note, due to MODFLOW row numbering when i increase ry1 and ry2 decrease
       cdum=cflk_W_par(cz,cz1,cz2)    ! This leads to different signs than in subcell routines!
       cfunc=cfunc+rs*cdum
c       write (iluer,1002) i,j,rs,cz,cz1,cz2,cdum,cfunc
c 1002 format (' cflk_W_actual2: i,j,rs / cz,cz1,cz2 / cdum, cfunc ',
c     &/,2(i5),1x,d14.7,/,6(d14.7),/,4(d14.7))
       rx1=rx2
      end do
      rx1=rlkx0
      ry1=ry1-rdy
      end do
c
c                    >>  add Poisson term
c
      call lkfindcell (cz,i0,j0,rx0,ry0)
      rq3=0.0
      if (i0.gt.0.and.j0.gt.0) then  ! CZ is inside cell i0,j0
       rs=(rlkleakage(i0,j0)+rlkrecharge(i0,j0))  ! add Poisson term for cell
       rx=REAL(cz)
       rqx=-rs*(rx-rx0)
       cfunc=cfunc+CMPLX(rqx,0.0)
       rbloc=rfbase(cz)
       rhght=rfhght(cz)
c       rq3=rs*((r3dz-rbloc)/rhght-1.0)  ! add vertical flux due to leakage OLD STATEMENT
       rq3=(r3dz-rbloc)/rhght*rs-rlkleakage(i0,j0) ! Compare Haitjema (1995), but....
c
c      Equation (3.312) contains Nt and Nb, which are both positive for inflow into the aquifer.
c      The variables rlkrecharge and rlkleakage are both positive for outflow from the aquifer.
c      Consequently, the signs on Nt (=rlkrecharge) and Nb (=rlkleakage) are reversed in the line
c      of code above.
c
c       write (iluer,1003) i0,j0,rs,r3dz,rbloc,rhght
c 1003 format (' cflk_W_actual3: i0,j0,rs,r3dz,rbloc,rhght ',/,
c     & 2(i3,1x),1x,4(d14.7,1x))
      end if
c                      !!!! no sub-grid contribution yet
      return
      end subroutine
c
c ------------------------------------------------------------------------------------------------------------
c
      COMPLEX(8) FUNCTION cflk_W_par (cz,cz1,cz2)
c
c ------------------------------------------------------------------------------------------------------------
c
C
C     Negative DERIVATIVE of the quadratic function for line doublet.
C     Repeat of the function CFDBSD (IEL), but not using DBPREP
C     (Function W=-dOmega/dz)
C
      INTEGER(4) I,ndbtrm
      LOGICAL lfar
      REAL(8) DCOF,dpi
      COMPLEX(8) cz,cz1,cz2,cdum1,cdum2,CBZ,cfbigz,cdi
      INCLUDE 'LUSYS.INC'
      DIMENSION DCOF(21)
      DATA ndbtrm,cdi,dpi /17,(0.0,1.0),.3141592653589793D+1/
      DATA DCOF /
     &             0.0000000000000000D+00,-0.3333333333333333D+00,
     &             0.0000000000000000D+00,-0.2000000000000000D+00,
     &             0.0000000000000000D+00,-0.1428571428571428D+00,
     &             0.0000000000000000D+00,-0.1111111111111111D+00,
     &             0.0000000000000000D+00,-0.9090909090909091D-01,
     &             0.0000000000000000D+00,-0.7692307692307693D-01,
     &             0.0000000000000000D+00,-0.6666666666666667D-01,
     &             0.0000000000000000D+00,-0.5882352941176471D-01,
     &             0.0000000000000000D+00,-0.5263157894736842D-01,
     &             0.0000000000000000D+00,-0.4761904761904762D-01,
     &             0.0000000000000000D+00/
C
            cbz=cfbigz(cz,cz1,cz2)
            cflk_W_par=CMPLX(0.0,0.0)
            lfar=ABS(cbz).gt.6.0 ! switch to farfield expansion at 3 times line-doublet length
      IF (.NOT.lfar) THEN
            if (abs(cbz-1.0).lt.0.001) cbz=cbz-0.001 ! shift onto line-doublet away from vertex
            if (abs(cbz+1.0).lt.0.001) cbz=cbz+0.001 ! shift onto line-doublet away from vertex
            cdum2=LOG((cbz-1.0)/(cbz+1.0))
            cdum1=(2.0D0+cbz*cdum2)
            cflk_W_par=-2.0D0*CDI*cdum1/DPI/(cz2-cz1)
      ELSE
c            write (iluer,1001)
c 1001 format (' cflk_W_par1: Farfield expansion in use.')
            cflk_W_par=CMPLX(0.0,0.0)
            cdum1=1.0/cbz
            cdum2=cdum1*cdum1
            DO 10 I=2,NDBTRM,2
            cflk_W_par=cflk_W_par+DCOF(I)*cdum2
            cdum2=cdum2*cdum2
  10        CONTINUE
            cflk_W_par=-cflk_W_par*4.0D0*CDI/DPI/(cz2-cz1)
      ENDIF
      RETURN
      END
c
c ---------------------------------------------------------------
c
      complex(8) function cflk_subomega (cz)
c
c     Return potential due to all subgrid leakage elements (sub-cells) in MODFLOW cells
c
      implicit none
      COMPLEX(8) cz,cfunc
      include 'lkcom.inc'
      include 'lusys.inc'
      cflk_subomega=CMPLX(0.0,0.0)
      if (lkleakage) then
       call cflk_subomega_sub (cz,cfunc)
       cflk_subomega=cfunc
c      write (iluer,1001) REAL(cfunc)
c 1001 format (' cflk_subomega1: real(cfunc) ',d14.7)
      end if
      return
      end

c
c --------------------------------------------------------------------------------
c
      subroutine cflk_subomega_actual (cz,cfunc,nrow,ncol,nbuf,
     &  rlkdeltax,rlkdeltay,ilkresolution,rlksubleakage,ilkpointer,
     &  rlkconst)
c
c    Routine returns the complex potential due to all sub-cells in all MODFLOW cells.
c    These calculations are organized per MODFLOW cell.
c
c    The Poisson term, added when CZ inside the sub-cell, is added to the real part
c    of the complex potention CFUNC
c
c
      implicit none
      INTEGER i,j,i0,j0,nrow,ncol,nbuf,ilkresolution,ires,ilkpointer,
     &        istart
      LOGICAL lx,ly,lseries
      REAL(8) rlkdeltax,rlkdeltay,rlksubleakage,rx,ry,rx0,ry0,rdx,rdy,
     &        rpotconst,rlkconst,rx1,rx2,ry1,ry2,rxx,ryy
      COMPLEX(8) cz,cz0,cflk_subomegaMODFLOWcell,cfunc
      DIMENSION rlkdeltax(ncol),rlkdeltay(nrow),ilkpointer(nrow,ncol),
     &          ilkresolution(nrow,ncol),rlkconst(nrow,ncol+1),
     &          rlksubleakage(nbuf)
      include 'lusys.inc'
      include 'lkcom.inc'
c
      cfunc=CMPLX(0.0d0,0.0d0)
      lseries=.false.
      ry=rlky0
      rx=rlkx0
      do i=2,nrow-1
      ry=ry-rlkdeltay(i-1) ! top side of cell
       do j=2,ncol-1
       rx=rx+rlkdeltax(j-1) ! left side of cell
       ires=ilkresolution(i,j)
       if (ires.gt.1) THEN ! MODFLOW cell has sub-cells, include these
         rdx=rlkdeltax(j)
         rdy=rlkdeltay(i)
         rx0=rx+0.5*rdx  ! center of cell
         ry0=ry-0.5*rdy  ! center of cell
         istart=ilkpointer(i,j)
         rpotconst=rlkconst(i,j)/ires
c              decide on use of series; do so if CZ more than one MODFLOW cell away
         rx1=rx-rlkdeltax(j-1)
         rx2=rx+rdx+rlkdeltax(j+1)
         ry1=ry-rdy-rlkdeltay(i+1)
         ry2=ry+rlkdeltay(i-1)
         rxx=REAL(cz)
         ryy=AIMAG(cz)
         lx=rxx.gt.rx2.or.rxx.lt.rx1
         ly=ryy.gt.ry2.or.ryy.lt.ry1
         lseries=lx.or.ly
c         lseries=.FALSE.                                            ! ********** forced series OFF --- debugging
         cfunc=cfunc+cflk_subomegaMODFLOWcell(cz,rx0,ry0,ires,
     &             rdx,rdy,nbuf,rlksubleakage,istart,rpotconst,lseries)
       end if
       end do
       rx=rlkx0
      end do
      RETURN
      end
c
c ----------------------------------------------------------------------------------
c
      COMPLEX(8) function cflk_subomegaMODFLOWcell(cz,rxc,ryc,ires,
     &                rdx,rdy,nbuf,rlksubleakage,istart,rconst,lseries)
c
c     Function returns the complex potential contribution of all sub-cells inside
c     the currently selected MODFLOW cell. The Poisson term is added to the real part of "omega."
c
c     Input:
c              cz         location where potential is to be evaluated
c              rxc, ryc   center of MODFLOW cell
c              rdx,rdy    MODFLOW cell dimensions
c              ires       sub-cell resolution (ires*ires sub-cells)
c              istart     starting address of sub-cell leakages in linear array: rlksubleakage
c              lseries    use series expansion; logical is set in cflk_subomega_actual
c
c     Note: We must organize sub-cells per MODFLOW cell with disregard of neighboring MODFLOW cells,
c           since those may or may not have sub-cells or have sub-cells at a different resolution.
c
c
c    Sub-cell order in a MODFLOW cell, thus order of leakages, upper heads, lower heads, etc.
c
c               _____________________
c              |       |      |      |
c              |  7    |  8   |  9   |
c              |_______|______|______|
c              |       |      |      |
c              |  4    |  5   |  6   |
c              |       |      |      |
c              |_______|______|______|
c              |       |      |      |
c              |  1    |  2   |  3   |
c              |       |      |      |
c              |_______|______|______|
c

      implicit none
      INTEGER i,j,i0,j0,ires,nbuf,istart,jstart,ncells,ierr,ncolpts,
     &        nterms
      LOGICAL l1,l2,lseries,l3,l4,l5,l6
      REAL(8) rx,ry,rxc,ryc,rx0,ry0,rx1,ry1,rx2,ry2,rdx,rdy,rdxx,rdyy,
     &        rlksubleakage,rs,rsm,rsp,rpot,rpotconst,rl,rconst,rsall,
     &        rphim,Rf
      COMPLEX(8) cz,cz1,cz2,cfunc,comls,cflk_omega_par,ca,cas,cdum,cz0,
     &           czm,cbz
      DIMENSION rlksubleakage(nbuf)
      ALLOCATABLE ca(:,:),cas(:),rsall(:),czm(:),rphim(:)
      include 'lusys.inc'
      include 'lkcom.inc'
      save    ! to keep the allocated arrays
c
      if (istart+ires*ires-1.gt.nbuf) THEN ! bail out, buffer overrun
       cflk_subomegaMODFLOWcell=CMPLX(0.0d0,0.0d0)
       write (iluer,1000) istart,ires,nbuf
       return
      endif
      if (.not.lseries) THEN ! return the complete function
       rdxx=rdx/ires
       rdyy=rdy/ires  ! sub-cell dimensions
       rpotconst=rconst
       rx0=rxc-0.5*rdx
       ry0=ryc-0.5*rdy  ! lower-left of MODFLOW cell
       cfunc=CMPLX(0.0d0,0.0d0)
c
c       >> include all line-sinks along vertical grid lines
c
       do j=1,ires+1
         rx1=rx0+(j-1)*rdxx
         rx2=rx1
         do i=1,ires
           ry1=ry0+(i-1)*rdyy
           ry2=ry0+i*rdyy
           cz1=CMPLX(rx1,ry1)
           cz2=CMPLX(rx2,ry2)
           jstart=istart+(i-1)*ires
           if (j.eq.1) then
             rsm=0.0d0
             rsp=rlksubleakage(jstart)*0.5*rdxx
           else if (j.eq.ires+1) then
             rsm=rlksubleakage(jstart+ires-1)*0.5*rdxx
             rsp=0.0d0
           else
             rsm=rlksubleakage(jstart+j-2)*0.5*rdxx
             rsp=rlksubleakage(jstart+j-1)*0.5*rdxx
           end if
           rs=rsm+rsp
           cfunc=cfunc+rs*(comls(cz,cz2,cz1)+rpotconst) ! rpotconst is fraction of rlkconst(irow,jcol)
c                                                       to make effect vanish "near" reference point
         end do                                         ! cz1 and cz2 interchanged to force line-sink to point upward
       end do                                           ! this is important to get the branch cuts to point upward
c
c       >> include all line-doublets along horizontal lines
c
       do i=1,ires+1
         ry1=ry0+(i-1)*rdyy
         ry2=ry1
         do j=1,ires
           rx1=rx0+(j-1)*rdxx
           rx2=rx0+j*rdxx
           cz1=CMPLX(rx1,ry1)
           cz2=CMPLX(rx2,ry2)
           jstart=istart+(i-1)*ires
           if (i.eq.1) then ! signs of rsm & rsp differs from MODFLOW grid routines due to difference of row numbering.
             rsm=0.0d0      ! here when i increases so does ry1 and ry2
             rsp=rlksubleakage(jstart+j-1)*rdxx*rdxx*0.125
           else if (i.eq.ires+1) then
             rsm=rlksubleakage(jstart-ires+j-1)*rdxx*rdxx*0.125
             rsp=0.0d0
           else
             rsm=rlksubleakage(jstart-ires+j-1)*rdxx*rdxx*0.125
             rsp=rlksubleakage(jstart+j-1)*rdxx*rdxx*0.125
           end if
            rs=rsm-rsp
           cfunc=cfunc-rs*cflk_omega_par(cz,cz1,cz2) ! here when i increases so does ry1 and ry2
         end do
       end do
c
c       >> include Poisson term
c
        rx=REAL(cz)
        ry=aimag(cz)
        l1=rx.lt.rxc-0.5*rdx.or.rx.gt.rxc+0.5*rdx
        l2=ry.lt.ryc-0.5*rdy.or.ry.gt.ryc+0.5*rdy
        if (.not.l1.and..not.l2) THEN ! we are inside the MODFLOW cell, now find the sub-cell.
          i0=1
          do i=1,ires
           ry2=ry0+(i-1)*rdyy
           if (ry.gt.ry2) then
             i0=i
           end if
          end do
          rx1=rx0
          j0=1
          do j=1,ires
            rx2=rx0+(j-1)*rdxx
            if (rx.gt.rx2) then
              rx1=rx2+0.5*rdxx ! need to be at center of sub-cell
              j0=j
            end if
          end do
          jstart=istart+(i0-1)*ires+j0-1
          rs=rlksubleakage(jstart)
          rdxx=0.5*rdxx ! rdxx is now half the sub-cell width
          rx=rx-rx1 ! translate origin to sub-cell center
          rpot=0.5*rs*(rx-rdxx)*(rx+rdxx)
          cfunc=cfunc+CMPLX(rpot,0.0d0)
        end if
c
       cflk_subomegaMODFLOWcell=cfunc ! done, now return
c
      else                                  ! return the series expansion !!!!
c
       ncells=ires*ires
       nterms=nlkterms  ! set in lkdata in lkmod.for (currently: 10)
       ncolpts=2*nterms ! set much larger than nterms
       l1=ALLOCATED(ca)
       l2=ALLOCATED(cas)
       l3=ALLOCATED(rsall)
       l4=ALLOCATED(czm)
       l5=ALLOCATED(rphim)
       l6=.false.
       if (l2) l6=size(rsall).eq.ncells ! note: we assume all MODFLOW cells are of the same size, then
c      if the resolution of the sub-grids are the same we do not have to regenerate!!
       if (l1.and.l2.and.l3.and.l4.and.l5.and.l6) then
c        OK, we have the correct arrays and of the correct size
       else ! must allocate and generate CAS
        if (l1) deallocate (ca)
        if (l2) deallocate (cas)
        if (l3) deallocate (rsall)
        if (l4) deallocate (czm)
        if (l5) deallocate (rphim)
        ALLOCATE (ca(nterms,ncells),cas(nterms),rsall(ncells),
     &            czm(ncolpts),rphim(ncolpts),stat=ierr)
        if (ierr.ne.0) then
          write (iluer,2000) nterms,ncells
          AMESS(1)='Failed to allocate series arrays.'
          AMESS(2)='Verify resources and rerun.'
          CALL HALT(2) ! stop program execution for batch version
        end if
        call lkgensubcellcoefmatrix (ca,czm,rphim,ncolpts,nterms,ncells,
     &               rxc,ryc,rdx,rdy,Rf,ires,rconst)
       end if
       rsall(1:ncells)=rlksubleakage(istart:istart+ncells-1) ! get the sub-cell leakages for this MODFLOW cell
       cas=MATMUL(ca,rsall) ! calculate the composite series coefficients
       cfunc=(0.0d0,0.0d0)
       cdum=1.0d0
       cz0=CMPLX(rxc,ryc)
       cbz=(cz-cz0)/Rf
       do i=1,nterms
        cfunc=cfunc+cas(i)*cdum
        cdum=cdum/cbz
       end do
       cflk_subomegaMODFLOWcell=cfunc ! done, now return
      end if
c
c
c     Note: If the sum of the sub-cell leakages does not add up to zero then we must add a logarithmic term.
c           This has not yet been implemented.
c
c
      return
 1000 format (' ***ERROR in cflk_subomegaMODFLOWcell:',
     &' impending buffer overrun. ',/,'istart, ires, nbuf= '3i5)
 2000 format (' ***ERROR in cflk_subomegaMODFLOWcell:',/,
     &' cannot allocate CA, CAS, RSALL, CZM, and RPHIM for ',i4,
     &' terms and ',i4, 'sub-cells',/,' Execution aborted.')
      end
c
c ---------------------------------------------------------------
c
      subroutine lksubqi (cz,rqi)
c
c     Return discharge vector due to all sub-cells in the MODFLOW grid
c     Note: rqi(1) and rqi(2) are the discharge vector components, while
c           rqi(3) is the vertical component of the specific discharge.c
c
c     Note: Only the sub-cell contributions of 9 MODFLOW cells are included;
c           from the MODFLOW cell inwhich CZ is located and 8 surrounding cells.
c
c
      implicit none
      REAL(8) rqi,rq3
      COMPLEX(8) cz,cfunc
      DIMENSION rqi(3)
      include 'lkcom.inc'
      include 'lusys.inc'
      if (lkleakage) then
      rq3=0.0
      call lk_subW_sub (cz,cfunc,rq3)
      rqi(1)=rqi(1)+REAL(cfunc)
      rqi(2)=rqi(2)-AIMAG(cfunc)
      rqi(3)=rqi(3)+rq3
      end if
c      write (iluer,1001) rq3,rqi(3)
c 1001 format (' lksubqi1: rq3,rqi(3) ',2(d14.7))
      return
      end
c
c --------------------------------------------------------------------------------
c
      subroutine lk_subW_actual (cz,cfunc,rq3,nrow,ncol,nbuf,
     &  rlkdeltax,rlkdeltay,ilkresolution,rlksubleakage,ilkpointer)
c
c     Routine adds all sub-cell contributions to the W function and RQ3
c     The calculations are organized by MODFLOW cell.
c
c     A contribution to RQ3 is added when CZ is inside the sub-cell.
c
      implicit none
      INTEGER i,j,i0,j0,nrow,ncol,nbuf,ilkresolution,ires,ilkpointer,
     &        istart
      LOGICAL lx,ly,lseries
      REAL(8) rlkdeltax,rlkdeltay,rlksubleakage,rx,ry,rx0,ry0,rdx,rdy,
     &        rq3,rq3MF,rx1,rx2,ry1,ry2,rxx,ryy
      COMPLEX(8) cz,cz0,cfunc,cfuncMF
      DIMENSION rlkdeltax(ncol),rlkdeltay(nrow),ilkpointer(nrow,ncol),
     &          ilkresolution(nrow,ncol),rlksubleakage(nbuf)
      include 'lusys.inc'
      include 'lkcom.inc'
c
      cfunc=CMPLX(0.0d0,0.0d0)
      rq3=0.0d0
      ry=rlky0
      rx=rlkx0
      do i=2,nrow-1
      ry=ry-rlkdeltay(i-1) ! top side of cell
       do j=2,ncol-1
       rx=rx+rlkdeltax(j-1) ! left side of cell
       ires=ilkresolution(i,j)
       if (ires.gt.1) THEN ! MODFLOW cell has sub-cells, include these
         rdx=rlkdeltax(j)
         rdy=rlkdeltay(i)
         rx0=rx+0.5*rdx  ! center of cell
         ry0=ry-0.5*rdy  ! center of cell
         istart=ilkpointer(i,j)
c              decide on use of series; do so if CZ more than one MODFLOW cell away
         rx1=rx-rlkdeltax(j-1)
         rx2=rx+rdx+rlkdeltax(j+1)
         ry1=ry-rdy-rlkdeltay(i+1)
         ry2=ry+rlkdeltay(i-1)
         rxx=REAL(cz)
         ryy=AIMAG(cz)
         lx=rxx.gt.rx2.or.rxx.lt.rx1
         ly=ryy.gt.ry2.or.ryy.lt.ry1
         lseries=lx.or.ly
         call lk_subWMODFLOWcell(cz,cfuncMF,rq3MF,rx0,ry0,ires,
     &             rdx,rdy,nbuf,rlksubleakage,istart,lseries)
         cfunc=cfunc+cfuncMF
         rq3=rq3+rq3MF
       end if
       end do
       rx=rlkx0
      end do
      RETURN
      end
c
c ----------------------------------------------------------------------------------
c
      subroutine lk_subWMODFLOWcell(cz,cfunc,rq3,rxc,ryc,ires,
     &                 rdx,rdy,nbuf,rlksubleakage,istart,lseries)
c
c     Function returns the W-function for all sub-cells inside
c     the currently selected MODFLOW cell. The Poisson term contribution Q1 and Q2 is added
c     as Q1-iQ2 to the complex W-function. The Poisson term also produces the vertical
c     specific discharge vector component rq3.
c
c     Input:
c              cz         location where potential is to be evaluated
c              rxc, ryc   center of MODFLOW cell
c              rdx,rdy    MODFLOW cell dimensions
c              ires       sub-cell resolution (ires*ires sub-cells)
c              istart     starting address of sub-cell leakages in linear array: rlksubleakage
c              lseries    use farfield series expansion instead of function
c
c     Output:
c
c              cfunc      complex W-function
c              rq3        vertical specific discharge component
c
c
c     Note: We must organize sub-cells per MODFLOW cell with disregard of neighboring MODFLOW cells,
c           since those may or may not have sub-cells or have sub-cells at a different resolution.
c
c
c    Sub-cell order in a MODFLOW cell, thus order of leakages, upper heads, lower heads, etc.
c
c               _____________________
c              |       |      |      |
c              |  7    |  8   |  9   |
c              |_______|______|_______
c              |       |      |      |
c              |  4    |  -5- |  6   |
c              |       |      |      |
c              |_______|______|_______
c              |       |      |      |
c              |  1    |  2   |  3   |
c              |       |      |      |
c              |_______|______|______|
c

      implicit none
      INTEGER i,j,i0,j0,ires,nbuf,istart,jstart,ncells,ierr,ncolpts,
     &        nterms
      LOGICAL l1,l2,lseries,l3,l4,l5,l6
      REAL(8) rx,ry,rxc,ryc,rx0,ry0,rx1,ry1,rx2,ry2,rdx,rdy,rdxx,rdyy,
     &        rlksubleakage,rs,rsm,rsp,rqx,rq3,rfhght,rhght,rbloc,
     &        rfbase,rconst,rsall,rphim,Rf
      COMPLEX(8) cz,cz1,cz2,cfunc,cdls,cflk_W_par,ca,cas,cz0,czm,cbz,
     &           cdum
      DIMENSION rlksubleakage(nbuf)
      ALLOCATABLE ca(:,:),cas(:),rsall(:),czm(:),rphim(:)
      include 'lusys.inc'
      include 'com3d.inc'
      include 'lkcom.inc'
      cfunc=CMPLX(0.0d0,0.0d0)
      rq3=0.0d0
      if (istart+ires*ires-1.gt.nbuf) THEN ! bail out, buffer overrun
       write (iluer,1000) istart,ires,nbuf
       return
      endif
      if (.not.lseries) THEN ! return the complete function
       rdxx=rdx/ires
       rdyy=rdy/ires  ! sub-cell dimensions
       rx0=rxc-0.5*rdx
       ry0=ryc-0.5*rdy  ! lower-left of MODFLOW cell
c
c       >> include all line-sinks along vertical grid lines
c
       do j=1,ires+1
         rx1=rx0+(j-1)*rdxx
         rx2=rx1
         do i=1,ires
           ry1=ry0+(i-1)*rdyy
           ry2=ry0+i*rdyy
           cz1=CMPLX(rx1,ry1)
           cz2=CMPLX(rx2,ry2)
           jstart=istart+(i-1)*ires
           if (j.eq.1) then
             rsm=0.0d0
             rsp=rlksubleakage(jstart)*0.5*rdxx
           else if (j.eq.ires+1) then
             rsm=rlksubleakage(jstart+ires-1)*0.5*rdxx
             rsp=0.0d0
           else
             rsm=rlksubleakage(jstart+j-2)*0.5*rdxx
             rsp=rlksubleakage(jstart+j-1)*0.5*rdxx
           end if
           rs=rsm+rsp
           cfunc=cfunc-rs*cdls(cz,cz1,cz2)  ! Note the minus sign: W=-dO/dz
         end do
       end do
c
c       >> include all line-doublets along horizontal lines
c
       do i=1,ires+1
         ry1=ry0+(i-1)*rdyy
         ry2=ry1
         do j=1,ires
           rx1=rx0+(j-1)*rdxx
           rx2=rx0+j*rdxx
           cz1=CMPLX(rx1,ry1)
           cz2=CMPLX(rx2,ry2)
           jstart=istart+(i-1)*ires
           if (i.eq.1) then ! sign differs from MODFLOW grid routines due to different row numbering.
             rsm=0.0d0      ! here when i increases so does ry1 and ry2
             rsp=rlksubleakage(jstart+j-1)*rdxx*rdxx*0.125
           else if (i.eq.ires+1) then
             rsm=rlksubleakage(jstart-ires+j-1)*rdxx*rdxx*0.125
             rsp=0.0d0
           else
             rsm=rlksubleakage(jstart-ires+j-1)*rdxx*rdxx*0.125
             rsp=rlksubleakage(jstart+j-1)*rdxx*rdxx*0.125
           end if
           rs=rsm-rsp
           cfunc=cfunc+rs*cflk_W_par(cz,cz1,cz2)
         end do
       end do
c
c       >> include Poisson term
c
        rx=REAL(cz)
        ry=aimag(cz)
        l1=rx.lt.rxc-0.5*rdx.or.rx.gt.rxc+0.5*rdx
        l2=ry.lt.ryc-0.5*rdy.or.ry.gt.ryc+0.5*rdy
        if (.not.l1.and..not.l2) THEN ! we are inside the MODFLOW cell, now find the sub-cell.
          i0=1
          do i=1,ires
           ry2=ry0+(i-1)*rdyy
           if (ry.gt.ry2) then
             i0=i
           end if
          end do
          rx1=rx0
          j0=1
          do j=1,ires
            rx2=rx0+(j-1)*rdxx
            if (rx.gt.rx2) then
              rx1=rx2+0.5*rdxx ! need to be at center of sub-cell
              j0=j
            end if
          end do
          jstart=istart+(i0-1)*ires+j0-1   ! we found the sub-cell i0,j0 - set address for linear buffers
          rs=rlksubleakage(jstart)
          rdxx=0.5*rdxx ! rdxx is now half the sub-cell width
          rx=rx-rx1 ! translate origin to sub-cell center
          rqx=-rs*rx
          cfunc=cfunc+CMPLX(rqx,0.0d0)
c
c     now calculate rq3, the vertical component of the specific discharge vector
c
          rbloc=rfbase(cz)
          rhght=rfhght(cz)
          rq3=rs*((r3dz-rbloc)/rhght-1.0d0)  ! there is no recharge term in sub-cells (part of MODFLOW cells)
        end if          ! done
c
      else              ! use farfield series
c
       rconst=0.0d0
       ncells=ires*ires
       nterms=nlkterms  ! set in lkdata in lkmod.for (currently: 10)
       ncolpts=2*nterms ! set much larger than nterms
       l1=ALLOCATED(ca)
       l2=ALLOCATED(cas)
       l3=ALLOCATED(rsall)
       l4=ALLOCATED(czm)
       l5=ALLOCATED(rphim)
       l6=.false.
       if (l2) l6=size(rsall).eq.ncells ! note: we assume all MODFLOW cells are of the same size, then
c      if the resolution of the sub-grids are the same we do not have to regenerate!!
       if (l1.and.l2.and.l3.and.l4.and.l5.and.l6) then
c        OK, we have the correct arrays and of the correct size
       else ! must allocate and generate CAS
        if (l1) deallocate (ca)
        if (l2) deallocate (cas)
        if (l3) deallocate (rsall)
        if (l4) deallocate (czm)
        if (l5) deallocate (rphim)
        ALLOCATE (ca(nterms,ncells),cas(nterms),rsall(ncells),
     &            czm(ncolpts),rphim(ncolpts),stat=ierr)
        if (ierr.ne.0) then
          write (iluer,2000) nterms,ncells
          AMESS(1)='Failed to allocate series arrays.'
          AMESS(2)='Verify resources and rerun.'
          CALL HALT(2) ! stop program execution for batch version
        end if
        call lkgensubcellcoefmatrix (ca,czm,rphim,ncolpts,nterms,ncells,
     &               rxc,ryc,rdx,rdy,Rf,ires,rconst)
       end if
       rsall(1:ncells)=rlksubleakage(istart:istart+ncells) ! get the sub-cell leakages for this MODFLOW cell
       cas=MATMUL(ca,rsall) ! calculate the composite series coefficients
       cz0=CMPLX(rxc,ryc)
       cbz=(cz-cz0)/Rf
       cdum=1.0/cbz
       do i=2,nterms
        cdum=cdum/cbz
        cfunc=cfunc+(i-1)*cas(i)*cdum
       end do
       cfunc=cfunc/Rf ! W=-dO/dZ*dZ/dz  --> dZ/dz=1/Rf
                  ! done, now return
      end if
c
      return
 1000 format (' ***ERROR in lk_subWMODFLOWcell:',
     &' impending buffer overrun. ',/,'istart, ires, nbuf= '3i5)
 2000 format (' ***ERROR in lk_subWMODFLOWcell:',/,
     &' cannot allocate CA, CAS, RSALL, CZM, and RPHIM for ',i4,
     &' terms and ',i4, 'sub-cells',/,' Execution aborted.')
      end
c
c -------------------------------------------------------------------------------
c
      COMPLEX(8) function cflk_omega_coefficient(cz,
     &                            cz0,rdxx,rdyy,rconst,linside)
c
c     Function rerturn the complex coefficient for a MODFLOW cell or sub-cell.
c
c     cz        point at which the potential is being calculated
c     cz0       is the lower left corner point of the cell or sub-cell
c     rdxx      is the width in x=direction of the cell
c     rdyy      is the width in y-direction of the cell
c     rconst    is a real constant added to the potential (for instance to make it vanish as a remote point)
c     linside   TRUE in case cz is located inside the cell.
c
      implicit none
      LOGICAL linside
      REAL(8) rconst,rdxx,rdyy,rx,ry
      COMPLEX(8) cz,cz0,cz1,cz2,cz3,cdum,cdum1,cdum2,comls,
     &           cflk_omega_par
      include 'lusys.inc'
c
c     cell structure:
c
c      ^ y
c      |
c      |
c     cz3-----------cz2
c      |             |
c      |             |
c      |             |
c      |             |
c     cz0-----------cz1  ----> x
c
      cz1=cz0+rdxx
      cz2=cz1+CMPLX(0.0,rdyy)
      cz3=cz2-rdxx
      cdum1=comls(cz,cz3,cz0)+rconst
      cdum1=cdum1+comls(cz,cz2,cz1)+rconst
      cdum2=cflk_omega_par(cz,cz0,cz1)-cflk_omega_par(cz,cz3,cz2)
      cflk_omega_coefficient=0.5*rdxx*cdum1+0.125*rdxx*rdxx*cdum2
      if (linside) then
      rx=REAL(cz-0.5*(cz0+cz2))
      ry=AIMAG(cz-0.5*(cz0+cz2))
      cdum=0.5*(rx-0.5*rdxx)*(rx+0.5*rdxx)
c      write (iluer,1002) rx,ry,cdum
c 1002 format (' cflk_omega_coefficient2: rx,ry,cdum=',4(d14.7))
      cflk_omega_coefficient=cflk_omega_coefficient+cdum
      end if
c
      return
      end
