C     Last change:  HMH  24 Mar 2008   12:59 pm
c
c   lkfindcell_actual  finds cell location and cell center in which CZ occurs
c   rflklowerhead      calls rfklowerhead_sub (entry in LKIN) to return lower head at CZ
c   rflklowerhead_sub_actual called by entry rflklowerhead_sub in LKIN
c   rflkleakage      calls rflkleakage_sub (entry in LKIN) to return leakage at CZ
c   rflkleakage_sub_actual called by entry rflkleakage_sub in LKIN
c   rflksubcellleakage returns leakage in sub-cells by calling rflksubcellleakage_sub (entry in LKIN)
c   rflksubcellleakage_sub_actual  called by rflksubcellleakage_sub (entry in LKIN)
c   rflkrecharge        returns the recharge in the MODFLOW cells by calling rflkrecharge_sub (entry in LKIN)
c   rflkrecharge_sub_actual called by rflkrecharge_sub (entry in LKIN)
c   rflkresistance      calls rflkresistance_sub (entry in LKIN) to return resistance at CZ
c   rflkresistance_sub_actual called by entry rflkresistance_sub in LKIN
c   lkinside            logical function that is true if CZ inside cell
c   lkfincludesubgrid   logical function to communicate lkincludesubgrid to outside routines
c   lkgensubcellcoefmatrix  routine generates the coefficients for the far-field series of subcells
c
c -------------------------------------------------------------------------
c
      subroutine lkfindcell_actual (cz,i0,j0,rx0,ry0,
     &                              rlkdeltax,rlkdeltay,nrow,ncol)
c
c    Routine returns the cell location and cell center in which CZ occurs
c    Input: CZ and the arrays rlkdeltax(ncol) and rlkdeltay(nrow)
c    Output: i0, j0 (row and column) and rx0, ry0 (center of cell)
c    i0 and/or j0 are set to zero when CZ is outside the grid.
c
c    Routine is called by entry lkfindcell(cz,i0,j0,rx0,ry0)
c
c
      implicit none
      INTEGER i,j,i0,j0,nrow,ncol
      REAL(8) rx,ry,rx0,ry0,rx1,ry1,rx2,ry2,rdx,rdy,rlkdeltax,rlkdeltay
      COMPLEX(8) cz
      DIMENSION rlkdeltax(ncol),rlkdeltay(nrow)
      INCLUDE 'lkcom.inc'
      include 'lusys.inc'
c
      i0=0
      ry=aimag(cz)
      ry1=rlky0
      do i=1,nrow       ! look below which horizontal grid line
       rdy=rlkdeltay(i)
       ry2=ry1-rdy
       if (ry.lt.ry1) then ! below upper gridline of cell, CZ may be in cell
        i0=i
        ry0=0.5*(ry1+ry2) ! y-value at center of cell
       endif
       ry1=ry2
      end do
      if (ry.lt.ry2) then ! below lower grid boundary, CZ is outside grid
        i0=0
      endif
      j0=0
      rx=REAL(cz)
      rx1=rlkx0
      do j=1,ncol       ! look to the right of which vertical grid line
       rdx=rlkdeltax(j)
       rx2=rx1+rdx
       if (rx.gt.rx1) then ! to the right of left-hand gridline of cell, CZ may be in cell
        j0=j
        rx0=0.5*(rx1+rx2)  ! x-value at center of cell
       endif
       rx1=rx2
      end do
      if (rx.gt.rx2) then  ! to the right of the right-hand boundary of the grid, CZ is outside grid.
        j0=0
      endif

      end subroutine
c
c -------------------------------------------------------------------------
c
      REAL(8) function rflklowerhead (cz)
c
c     Function returns lower head at CZ
c
      implicit none
      REAL(8) rfunc
      COMPLEX(8) cz
c
      call rflklowerhead_sub (cz,rfunc)
      rflklowerhead=rfunc
      return
c
      end
c
c --------------------------------------------------------------------------
c
      subroutine rflklowerhead_sub_actual (cz,rfunc,nrow,
     &                  ncol,nbuf,ilkpointer,rlksubheadlower,
     &                  rlkdeltax,rlkdeltay,rlkheadlower,ilkresolution)
c
c     Routine is called by rflklowerhead_sub in LKIN
c     Lower head in cell or sub-cell in which CZ occurs is returned.
c
c
      implicit none
      INTEGER nrow,ncol,i0,j0,ilkresolution,ilkpointer,ires,i,j,
     &        icell,jcell,nbuf,iadd
      REAL(8) rfunc,rlkdeltax,rlkdeltay,rlkheadlower,rx0,ry0,
     &        rlksubheadlower,rdxx,rdyy,rxc,ryc,rx,ry,rx2,ry2
      COMPLEX(8) cz
      DIMENSION rlkdeltax(ncol),rlkdeltay(nrow),rlkheadlower(nrow,ncol),
     &          ilkresolution(nrow,ncol),ilkpointer(nrow,ncol),
     &          rlksubheadlower(nbuf)
      INCLUDE 'lusys.inc'
      include 'lkcom.inc'
c
      rfunc=-9999.0  ! flag for being outside the grid.
      call lkfindcell (cz,icell,jcell,rxc,ryc)
      if (icell.gt.0.and.jcell.gt.0) then ! OK we are inside a MODFLOW cell
        iadd=ilkpointer(icell,jcell)
        if (iadd.gt.0) THEN ! MODFLOW cell contains subcells; assume interpolated lower heads
          ires=ilkresolution(icell,jcell)
          rdxx=rlkdeltax(jcell)/ires
          rdyy=rlkdeltay(icell)/ires
          rx0=rxc-rlkdeltax(jcell)/2.0
          ry0=ryc-rlkdeltay(icell)/2.0
          rx=REAL(cz)
          ry=AIMAG(cz)
          i0=1
          do i=1,ires
           ry2=ry0+(i-1)*rdyy
           if (ry.gt.ry2) then
             i0=i
           end if
          end do
          j0=1
          do j=1,ires
            rx2=rx0+(j-1)*rdxx
            if (rx.gt.rx2) then
              j0=j
            end if
          end do
          iadd=iadd+(i0-1)*ires+j0-1
          rfunc=rlksubheadlower(iadd)
        else
          rfunc=rlkheadlower(icell,jcell)  ! no interpolation
        endif
      end if
c
      return
      END subroutine
c
c -------------------------------------------------------------------------
c
      REAL(8) function rflkleakage (cz)
c
c     Function returns leakage at CZ
c
      implicit none
      REAL(8) rfunc
      COMPLEX(8) cz
      include 'lkcom.inc'
c
      rfunc=0.0d0
      if (lkleakage) then
      call rflkleakage_sub (cz,rfunc)
      endif
      rflkleakage=rfunc
      return
c
      end

c
c --------------------------------------------------------------------------
c
      subroutine rflkleakage_sub_actual (cz,rfunc,nrow,
     &                    ncol,rlkdeltax,rlkdeltay,rlkleakage)
c
c     Routine is called by rflkleakage_sub in LKIN
c     Leakage in cell in which CZ occurs is returned.
c     NOTE: no sub-grids in this version
c
c
      implicit none
      INTEGER nrow,ncol,i0,j0
      REAL(8) rfunc,rlkdeltax,rlkdeltay,rlkleakage,rx0,ry0
      COMPLEX(8) cz
      DIMENSION rlkdeltax(ncol),rlkdeltay(nrow),rlkleakage(nrow,ncol)
      INCLUDE 'lusys.inc'
      include 'lkcom.inc'
c
      rfunc=0.0   ! default zero leakage outside the grid.
      call lkfindcell (cz,i0,j0,rx0,ry0)
      if (i0.gt.0.and.j0.gt.0) then
        rfunc=rlkleakage(i0,j0)
      end if
c
      return
      END subroutine
c
c -------------------------------------------------------------------------
c
      REAL(8) function rflksubcellleakage (cz)
c
c     Function returns leakage at CZ
c
      implicit none
      REAL(8) rfunc
      COMPLEX(8) cz
      include 'lkcom.inc'
c
      rfunc=0.0d0
      if (lkleakage) then
      call rflksubcellleakage_sub (cz,rfunc)
      endif
      rflksubcellleakage=rfunc
      return
c
      end

c
c --------------------------------------------------------------------------
c
      subroutine rflksubcellleakage_sub_actual (cz,rfunc,nrow,
     &                    ncol,nbuf,ilkpointer,rlksubleakage,
     &                    rlkdeltax,rlkdeltay,ilkresolution)
c
c     Routine is called by rflkleakage_sub in LKIN
c     Leakage in cell in which CZ occurs is returned.
c     NOTE: no sub-grids in this version
c
c
      implicit none
      INTEGER nrow,ncol,icell,jcell,i0,j0,ilkpointer,nbuf,
     &        i,j,iadd,ilkresolution,ires
      REAL(8) rfunc,rlksubleakage,rx0,ry0,rx2,ry2,rx,ry,
     &        rlkdeltax,rlkdeltay,rdxx,rdyy,rxc,ryc
      COMPLEX(8) cz
      DIMENSION ilkpointer(nrow,ncol),rlksubleakage(nbuf)
      DIMENSION rlkdeltax(ncol),rlkdeltay(nrow)
      DIMENSION ilkresolution(nrow,ncol)
      INCLUDE 'lusys.inc'
      include 'lkcom.inc'
c
      rfunc=0.0   ! default zero leakage outside the grid.
      call lkfindcell (cz,icell,jcell,rxc,ryc)
      if (icell.gt.0.and.jcell.gt.0) then ! inside a MODFLOW cell
        iadd=ilkpointer(icell,jcell)
c      write (iluer,1001) icell,jcell,rxc,ryc,iadd
c 1001 format (' rflksubcellleakage_sub_actual1: ',
c     & 'icell,jcell,rxc,ryc,iadd ',/,
c     &          2(I5),2(d14.7),i5)
        if (iadd.gt.0) THEN ! MODFLOW cell contains subcells
          ires=ilkresolution(icell,jcell)
          rdxx=rlkdeltax(jcell)/ires
          rdyy=rlkdeltay(icell)/ires
          rx0=rxc-rlkdeltax(jcell)/2.0
          ry0=ryc-rlkdeltay(icell)/2.0
          rx=REAL(cz)
          ry=AIMAG(cz)
c          write (iluer,1002) rx,ry,rx0,ry0,rdxx,rdyy
c 1002 format (' rflksubcellleakage_sub_actual2: rx,ry,rx0,ry0,rdxx,rdyy'
c     &,/,6(d14.7))
          i0=1
          do i=1,ires
           ry2=ry0+(i-1)*rdyy
           if (ry.gt.ry2) then
             i0=i
           end if
          end do
          j0=1
          do j=1,ires
            rx2=rx0+(j-1)*rdxx
            if (rx.gt.rx2) then
              j0=j
            end if
          end do
          iadd=iadd+(i0-1)*ires+j0-1
          rfunc=rlksubleakage(iadd)
c          write (iluer,1003) iadd,i0,j0,icell,jcell,rfunc
c 1003 format (' rflksubcellleakage_sub_actual3: ',
c     &'iadd,i0,j0,icell,jcell,rfunc ',5(I5),d14.7)
        end if
      end if
c
      return
      END subroutine

c
c -------------------------------------------------------------------------
c
      REAL(8) function rflkrecharge (cz)
c
c     Function returns recharge in "leakage" grid at CZ
c
      implicit none
      REAL(8) rfunc
      COMPLEX(8) cz
      include 'lkcom.inc'
c
      rfunc=0.0d0
      if (lkleakage) then
      call rflkrecharge_sub (cz,rfunc)
      endif
      rflkrecharge=rfunc
      return
c
      end

c
c --------------------------------------------------------------------------
c
      subroutine rflkrecharge_sub_actual (cz,rfunc,nrow,
     &                    ncol,rlkdeltax,rlkdeltay,rlkrecharge)
c
c     Routine is called by rflkrecharge_sub in LKIN
c     Recharge in cell in which CZ occurs is returned.
c     NOTE: This is the recharge as obtained from MODFLOW (no subdivision -> no subgrid recharge)
c
      implicit none
      INTEGER nrow,ncol,i0,j0
      REAL(8) rfunc,rlkdeltax,rlkdeltay,rlkleakage,rlkrecharge,rx0,ry0
      COMPLEX(8) cz
      DIMENSION rlkdeltax(ncol),rlkdeltay(nrow),rlkrecharge(nrow,ncol)
      INCLUDE 'lusys.inc'
      include 'lkcom.inc'
c
      rfunc=0.0   ! default zero recharge outside the grid.
      call lkfindcell (cz,i0,j0,rx0,ry0)
      if (i0.gt.0.and.j0.gt.0) then
        rfunc=rlkrecharge(i0,j0)
      end if
c
      return
      END subroutine

c
c -------------------------------------------------------------------------
c
      REAL(8) function rflkresistance (cz)
c
c     Function returns resistance at CZ.
c
      implicit none
      REAL(8) rfunc
      COMPLEX(8) cz
c
      call rflkresistance_sub (cz,rfunc)
      rflkresistance=rfunc
      return
c
      end
c
c --------------------------------------------------------------------------
c
      subroutine rflkresistance_sub_actual (cz,rfunc,nrow,
     &                    ncol,rlkdeltax,rlkdeltay,rlkresist)
c
c     Routine is called by rflkresistance_sub in LKIN
c     Resistance in cell in which CZ occurs is returned.
c     NOTE: no sub-grid used in this routine.
c
c
      implicit none
      INTEGER nrow,ncol,i0,j0
      REAL(8) rfunc,rlkdeltax,rlkdeltay,rlkresist,rx0,ry0
      COMPLEX(8) cz
      DIMENSION rlkdeltax(ncol),rlkdeltay(nrow),rlkresist(nrow,ncol)
      INCLUDE 'lusys.inc'
      include 'lkcom.inc'
c
      rfunc=1.0d+30 ! default infinite resistance outside the grid.
      call lkfindcell (cz,i0,j0,rx0,ry0)
      if (i0.gt.0.and.j0.gt.0) then
        rfunc=rlkresist(i0,j0)
      end if
c
      return
      END subroutine
c
c --------------------------------------------------------------------------
c
      LOGICAL function lkinside(cz,cz0,cz1,cz2,cz3)       ! designed for lkgenmat, but not used there
c
c     True if cz inside cell, which includes the upper and right-hand boundary.
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
      LOGICAL l1,l2
      REAL(8) rx,rx1,rx2,ry,ry1,ry2
      COMPLEX(8) cz,cz0,cz1,cz2,cz3
      include 'lkcom.inc'
c
      lkinside=.false.
      if (lkleakage) then
      rx=REAL(cz)
      ry=AIMAG(cz)
      rx1=REAL(cz0)
      ry1=AIMAG(cz0)
      rx2=REAL(cz2)
      ry2=AIMAG(cz2)
      l1=rx.gt.rx1.and.rx.le.rx2
      l2=ry.gt.ry1.and.ry.le.ry2
      lkinside=l1.and.l2
      end if
      return
      end
c
c ------------------------------------------------------------------------------------
c
      LOGICAL function lflkincludesubgrid()
c
c     True if sub-cell strengths are to be calculated.
c
      implicit none
      INTEGER nsolOut
      LOGICAL lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut
      CHARACTER(8) aBasenameOut
      CHARACTER(16)aDateTimeOut
      include 'lkcom.inc'
c
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
c
c
      lflkincludesubgrid=lkincludesubgrid.and.nsolOut.ge.ilkiniter
c
      return
      end
c
c ---------------------------------------------------------------------------------------
c
      subroutine lkgensubcellcoefmatrix (ca,czm,rphim,ncolpts,
     &                 nterms,ncells,rxc,ryc,rdx,rdy,Rf,ires,rconst)
c
c     called in cflk_subomegaMODFLOWcell in lkfun.for
c     Routine controls the generation of a coefficient matrix for sub-cell farfield expansions.
c     All sub-cells are expanded about the same point; the center of the MODFLOW cell to which they belong.
c     This routine generates a matrix of coefficients for all sub-cells using a strength = 1
c     The calling routine will multiply this matrix with a vector of sub-cell strenghts to
c     obtain the final coefficients for the far-field expansion of the entire sub-grid in the MODFLOW cell.
c
c     ca(nterms,ncells) matrix of complex constants for the series  (Note: ca is allocated in
c                       cflk_subomegaMODFLOWcell in lkfun.for)
c     czm(ncolpts)      Zm as defined below (work array allocated in
c                       cflk_subomegaMODFLOWcell in lkfun.for)
c     rphim(ncolpts)    Phi at m-th collocation point (work array allocated in
c                       cflk_subomegaMODFLOWcell in lkfun.for)
c
c     Omega(Z)=ca(1,n) + ca(2,n)/(Z) + ca(3,n)/(Z)^2 + ca(4,n)/(Z)^3 + .........
c
c     Z=(z-z0)/Rf
c
c     cz0   =cmplx(rxc,ryc)
c     czm   m-th control point (later the point where Omega is to be calculated).
c     Rf    radius outside of which the series convergences.
c     ires  sub-grid resolution (there are ires^2 sub-cells).
c     rconst average constant for the line-sink fuctions in the cells.
c
c     coefficient generation based on equations 12 - 18 in "The superblock
c     approach for the analytic element method" O.D.L. Strack, I. Jankovic, and R. Barnes
c     Journal of Hydrology 226 (1999) 179-187.
c
c     Note the logarithmic term in eq. 12, and 16-18 has been left off as the sub-cell strengths are
c     expected to add up to zero; QT=0. Moreover, the logarithmic terms in 16-18 are zero in view of 14.
c
      implicit none
      integer ires,nterms,ncells,ncolpts,
     &        isubcell,jsubcell,ncell,m,n
      LOGICAL linside,ldeb1
      REAL(8) rxc,ryc,rdx,rdy,Rf,rteta,rphim,rconst,
     &        rdxx,rdyy,rxlowleft,rylowleft,rar,rai,
     &        rpsi_func,rpsi_series,rcosine,rsine,rconstant
      COMPLEX(8) ca,czm,cz0,cz,cteta,czlowleft,cflk_omega_coefficient,
     &           czterm,comega,cdum,cbz,clog
      DIMENSION rphim(ncolpts),czm(ncolpts),ca(nterms,ncells)
      DATA ldeb1 /.false./
      include 'lusys.inc'
c
c      write (iluer,1110) rxc,ryc
c 1110 format (' lkgensubcellcoefmatrix0: entered with cz0=',2(d14.6))
      Rf=1.9*rdx ! assumes that all MODFLOW cells are squares of the same size !! **limitation**
      cz0=cmplx(rxc,ryc)
      linside=.false.
      rconstant=rconst
c
      rdxx=rdx/ires
      rdyy=rdy/ires
      ncell=0
      do isubcell=1,ires
       rylowleft=AIMAG(cz0)-rdy/2+(isubcell-1)*rdyy
       do jsubcell=1,ires
        rxlowleft=REAL(cz0)-rdx/2+(jsubcell-1)*rdxx
        czlowleft=CMPLX(rxlowleft,rylowleft)
        do m=1,ncolpts
         rteta=6.2831853*(m-1)/ncolpts
         cteta=CMPLX(0.0,rteta)
         czm(m)=cz0+Rf*exp(cteta)
         rphim(m)=REAL(cflk_omega_coefficient(czm(m),
     &                     czlowleft,rdxx,rdyy,rconstant,linside))
        end do
        ncell=ncell+1
        rar=0.0
        do m=1,ncolpts
         rar=rar+rphim(m)
        end do
        rar=rar/ncolpts ! real part of first coefficient
        rai=0.0         ! imaginary part of first coefficient to be determined later in this routine
        ca(1,ncell)=CMPLX(rar,rai)
        do n=1,nterms-1  ! nterms=Ns+1 with the series from 0 - Ns
         rar=0.0
         rai=0.0
         do m=1,ncolpts
          rteta=6.2831853*(m-1)/ncolpts
          rcosine=COS(n*rteta)
          rsine=SIN(n*rteta)
          rar=rar+rphim(m)*rcosine
          rai=rai+rphim(m)*rsine
         end do
         rar=2*rar/ncolpts
         rai=2*rai/ncolpts
         ca(n+1,ncell)=CMPLX(rar,rai)
        end do
c
c       now find the imaginary part of ca(1,ncell) using continuity of the stream function across the
c       circle with radius Rf
c
        comega=(0.0,0.0)
        cbz=(czm(1)-cz0)/Rf ! evaluate series at the first control point
        czterm=1.0
        do n=1,nterms
         comega=comega+ca(n,ncell)*czterm
         czterm=czterm/cbz
        end do
        rpsi_series=AIMAG(comega)
        rpsi_func=AIMAG(cflk_omega_coefficient(czm(1),
     &                     czlowleft,rdxx,rdyy,rconstant,linside))
        ca(1,ncell)=ca(1,ncell)+CMPLX(0.0,rpsi_func-rpsi_series)
c
        if (ldeb1.and.ncell.eq.18) then
c
c       test series at some point  when ldeb1 is set true
c
        cteta=CMPLX(0.0,1.0)
        cz=cz0+2.0*Rf*EXP(cteta) ! test at some point on circle
c        cz=czm(1)
        cbz=(cz-cz0)/Rf
        write (iluer,1001) rdx,Rf,cz,cz0,czlowleft,cbz
 1001 format (' lkgensubcellcoefmatrix1: rdx,Rf,cz,cz0,czlowleft,cbz ',/
     &         ,10(d14.7))
        czterm=1.0
        comega=0.0
        do n=1,nterms
         comega=comega+ca(n,ncell)*czterm
         czterm=czterm/cbz
        end do
        clog=rdxx*rdyy/6.2831853*LOG(cbz)
        comega=comega+clog ! must add this term since strength = 1 (not zero)
        cdum=cflk_omega_coefficient(cz,
     &                     czlowleft,rdxx,rdyy,rconstant,linside)
        write (iluer,1002) clog,cdum,comega
 1002 format ('  lkgensubcellcoefmatrix2: clog,cdum,comega ',
     &         3(2(d14.7),2x))
        ldeb1=.false.   ! write only once
        endif
       end do
      end do
c
      return
      end subroutine

