C     Last change:  HMH  17 Nov 2008    8:32 am
c
c --------------------------------------------------------------------------
c
      REAL(8) function rfnflk (cz1,cz2)
c
c     Routine returns the flow across a line between cz1 and cz2 due to leakage
c     and recharge grid, both MODFLOW grid and subgrids.
c
c     NOTE: Subgrid contributions are not yet included!!!!!!!
c
      implicit none
      REAL(8) rfunc
      COMPLEX(8) cz1,cz2
      include 'lusys.inc'
c      write (iluer,1001) cz1,cz2
c 1001 format (' rfnflk1: Entering with cz1 and cz2= ',4(d14.7))
      rfnflk=0.0d0
      call lkfn_sub(cz1,cz2,rfunc)
      rfnflk=rfunc
c      write (iluer,1002) rfunc
c 1002 format (' rfnflk2: rfunc=',d14.7)
c
      return
      end
c
c -----------------------------------------------------------------------------
c
      subroutine lkfnf_sub_actual (cz1,cz2,rfunc,nrow,ncol,rlkdeltax,
     &                             rlkdeltay,rlkleakage,rlkrecharge)
c
c     Routine assigns the flow that crosses the line cz1 - cz2 due to all MODFLOW grid cells to rfunc
c     This is done in three steps:
c     Step 1 calculate all Laplace contributions as PSI_1 - PSI_2
c     Step 2 correct for all branch cuts along line-sinks (vertical grid lines) in and above the grid.
c     Step 3 add all flows across the line due to the Poisson contributions (recharge and leakage)
c

      implicit none
      INTEGER nrow,ncol,i,j,icase
      LOGICAL lvinsect,lhinsect
      REAL(8) rlkdeltax,rlkdeltay,rlkleakage,rlkrecharge,rfunc, d_one,
     &        rx,ry,rdx,rdy,rx0,ry0,rx1,ry1,rx2,ry2,rxg1,ryg1,rxg2,ryg2,
     &       rs,rsm,rsp,rsign,rxx1,ryy1,rxx2,ryy2,rdeltapsi,ra,rm,rmsign
      COMPLEX(8) cz,cz1,cz2,cz0,czz1,czz2,czseg1,czseg2,cflk_omega,
     &           czj,czi
      DIMENSION rlkdeltax(ncol),rlkdeltay(nrow),rlkleakage(nrow,ncol),
     &          rlkrecharge(nrow,ncol)
      include 'lkcom.inc'
      include 'lusys.inc'
c
c      write (iluer,1001)
c 1001 format (' lkfnf_sub_actual1: Entering.')
      rxg1=rlkx0    ! set grid dimensions
      ryg2=rlky0
      rxg2=rxg1+SUM(rlkdeltax(1:ncol))
      ryg1=ryg2-SUM(rlkdeltay(1:nrow))
      rx1=REAL(cz1) ! set line segment end points
      ry1=aimag(cz1)
      rx2=REAL(cz2)
      ry2=aimag(cz2)
      rsign=1.0       ! for now assume line entered from left to right
      if (rx1.gt.rx2) THEN ! reverse points, but remember to change sign on flow across line
      rx1=REAL(cz2)
      ry1=aimag(cz2)
      rx2=REAL(cz1)
      ry2=aimag(cz1)
      rsign=-1.0      ! found that line is entered from right to left, multiplyer for nowmal flow
      end if
      czz1=CMPLX(rx1,ry1) ! redefine line segment end points
      czz2=CMPLX(rx2,ry2) ! use local variables so as to not change the input parameters cz1 and cz2
      if (rx1.ne.rx2) then ! calculate slope "rm"
       rm=(ry2-ry1)/(rx2-rx1)
      ELSEIF (ry1.lt.ry2) then ! line segment is vertical
       rm=1.0d30 ! vertically upward
      else
       rm=-1.0d30 ! vertically downward
      end if
      d_one = 1.0
      rmsign=SIGN(d_one,rm)
      if (rx2.le.rxg1) then
       icase=1 ! line segment entirely to the left of the grid (no branch cuts)
      elseif (rx1.ge.rxg2) then
       icase=1 ! line segment entirely to the right of the grid (no branch cuts)
      elseif (ry1.le.ryg1.and.ry2.le.ryg1) then
       icase=1 ! line segment entirely below the grid  (no branch cuts)
      elseif (ry1.ge.ryg2.and.ry2.ge.ryg2) then
       icase=2 ! line segment entirely above the grid (may cross branch cuts!!)
      else
       icase=3 ! at least part of line segment is inside the grid (may cross branch cuts
               ! and add Poisson terms)
      endif
c      write (iluer,1002) rm,rsign,icase
c 1002 format (' lkfnf_sub_actual2: rm, rsign, icase ',2(d14.7),i3)
c
      rfunc=aimag(cflk_omega(czz1)-cflk_omega(czz2)) ! psi1-psi2 for all cases *** STEP 1 ***
c
c      write (iluer,1003) rfunc
c 1003 format(' lkfnf_sub_actual3: rfunc after PSI_1 - PSI_2: ',d14.7)
      if (icase.ne.1) then  ! cases 2 and 3: add branch cut due to line-sinks along vertical grid lines
      rx=rxg1
      do j=1,ncol+1
      if (rx1.lt.rx.and.rx2.gt.rx) THEN ! we are intersecting a branch cut on gridline j
       ry=ry1+rm*(rx-rx1) ! get y-value of intersect
       call lkadd_deltapsi (j,ry,rdeltapsi,nrow,ncol,rlkdeltax,
     &                      rlkdeltay,rlkleakage,rlkrecharge)
c
       rfunc=rfunc+rdeltapsi  ! adding branch cuts                              *** STEP 2 ***
c
c      write (iluer,1004) icase,j,rx,ry,rdeltapsi
c 1004 format (' lkfnf_sub_actual4: icase,j,rx,ry,rdeltapsi ',
c     &         2(i4),3(D14.7))
      end if
      rx=rx+rlkdeltax(j)
      end do
      end if
c
      if (icase.eq.3) then  ! must add poisson term(s)
c
c      First deal with possibility that parts of the line element are outside of the grid.
c      For this case we will redefine the end points on the grid boundary and ignore the outside
c      sections, since they do not have a Poisson contribution.
c
       if (rx1.lt.rxg1) then ! czz1 to left of the grid
        ry=ry1+rm*(rxg1-rx1)
        if (ryg1.le.ry.and.ryg2.ge.ry) THEN ! intersection point on left-hand grid boundary,
         czz1=CMPLX(rxg1,ry)                ! ignore outside segment
        end if
       end if
       if (rx2.gt.rxg2) THEN ! czz2 to the right of the grid
        ry=ry1+rm*(rxg2-rx1)
        if (ryg1.lt.ry.and.ryg2.gt.ry) THEN ! intersection point on right-hand grid boundary,
         czz2=CMPLX(rxg2,ry)                ! ignore outside segment
        end if
       end if
       if (rm.gt.0) then
        if (ry1.lt.ryg1) then ! czz1 below the grid
         rx=rx1+(ryg1-ry1)/rm
         if (rxg1.lt.rx.and.rxg2.gt.rx) THEN ! intersection point on bottom of grid,
          czz1=CMPLX(rx,ryg1)                ! ignore outside segment
         end if
        end if
        if (ry2.gt.ryg2) THEN ! czz2 above the grid
         rx=rx1+(ryg2-ry1)/rm
         if (rxg1.lt.rx.and.rxg2.gt.rx) THEN ! intersection point on top of grid,
          czz2=CMPLX(rx,ryg2)                ! ignore outside segment
         end if
        end if
       else  ! (rm.lt.0)
        if (ry1.gt.ryg2) THEN ! czz1 above the grid
         rx=rx1+(ryg2-ry1)/rm
         if (rxg1.lt.rx.and.rxg2.gt.rx) THEN ! intersection point on top of grid,
          czz1=CMPLX(rx,ryg2)                ! ignore outside secgment
         end if
        end if
        if (ry2.lt.ryg1) THEN ! czz2 below the grid
         rx=rx1+(ryg1-ry1)/rm
         if (rxg1.lt.rx.and.rxg2.gt.rx) THEN ! intersection point on bottom of grid,
          czz2=CMPLX(rx,ryg1)                ! ignore outside segment
         end if
        end if
       end if
       rx1=REAL(czz1)   ! for the purpose of adding the Poisson contribution
       ry1=AIMAG(czz1)  ! the line element is chopped to fit entirely inside the grid
       rx2=REAL(czz2)   ! (possibly with points on the grid boundary)
       ry2=AIMAG(czz2)
c      write (iluer,1005) czz1,czz2
c 1005 format (' lkfnf_sub_actual5: czz1,czz2 after clipping to grid ',
c     &         4(d14.7))
       czseg1=czz1
c
c     Now break up the line element into sections that fit inside a single cell
c
  10  continue
      lvinsect=.false.
      lhinsect=.false.
c      Find the nearest intersection with a vertical grid line
       rx=rxg2-rlkdeltax(ncol)
       do j=ncol,2,-1 ! loop left to right, so we end with the nearest intersection
       if (rx1.lt.rx.and.rx2.gt.rx) THEN ! yes, we are intersecting vertical gridline j
        ry=ry1+rm*(rx-rx1) ! get y-value of intersect
        czj=CMPLX(rx,ry) ! make new tentative end point
        lvinsect=.true.
       endif
       rx=rx-rlkdeltax(j)
       enddo
c      Find nearest intersection with a horizontal grid line
       if (rm.gt.0) THEN ! line slopes upward
       ry=ryg2-rlkdeltay(1)
       do i=2,nrow  ! loop top to bottom, so we end with nearest intersection
        if (ry1.lt.ry.and.ry2.gt.ry) THEN ! we are intersecting horizontal line i
         rx=rx1+(ry-ry1)/rm ! get x-value of intersect
         czi=CMPLX(rx,ry) ! make new tentative end point
         lhinsect=.true.
        end if
        ry=ry-rlkdeltay(i)
       enddo
       else ! (rm.lt.0)  ! line slopes downward
       ry=ryg1+rlkdeltay(nrow)
       do i=nrow,2,-1 ! loop from bottom to top, so we end up with nearest intersection
        if (ry1.gt.ry.and.ry2.lt.ry) THEN ! we are intersecting horizontal line i
         rx=rx1+(ry-ry1)/rm ! get x-value of intersect
         czi=CMPLX(rx,ry) ! make new tentative end point
         lhinsect=.true.
        end if
        ry=ry+rlkdeltay(i)
       enddo
       end if
      if (lvinsect.and.lhinsect) then ! two new potential end points, select closest to czz1
c      now which is closest?
       if (cdabs(czj-czz1).lt.cdabs(czi-czz1)) then
        czseg2=czj
       else
        czseg2=czi
       end if
      else if (lvinsect) then ! select point on vertical grid line
        czseg2=czj
      else if (lhinsect) then ! select point on horizontal grid line
        czseg2=czi
      else                    ! keep original end point (no intersect)
        czseg2=czz2
      end if
c      write (iluer,1006) lvinsect,lhinsect,czseg1,czseg2
c 1006 format (' lkfnf_sub_actual6: lvinsect,lhinsect,czseg1,czseg2 ',
c     &           2(l3),4(d14.7))
c
c      now add Poisson contribution for the line segment czseg1 - czseg2
c
      cz0=0.5*(czseg1+czseg2)
      call lkfindcell_actual (cz0,i,j,rx0,ry0,
     &                        rlkdeltax,rlkdeltay,nrow,ncol)
      rxx1=REAL(czseg1)
      ryy1=AIMAG(czseg1)
      rxx2=REAL(czseg2)
      ryy2=AIMAG(czseg2)
      ra=0.5*((rxx1-rx0)+(rxx2-rx0))*ABS(ryy2-ryy1)
      rfunc=rfunc-rmsign*(rlkleakage(i,j)+rlkrecharge(i,j))*ra ! add Poisson contribution *** STEP 3 ***
c      note: recharge and leakage are negative when entering the aquifer
c      note: when we slope downward,the area must be given a different sign to keep positive flow convention
c      write (iluer,1007) i,j,ra
c 1007 format (' lkfnf_sub_actual7: i,j,ra ',2(i3),d14.7)
      rx1=REAL(czseg2) ! drop segment just considered
      ry1=AIMAG(czseg2)
      czseg1=czseg2
      if (cdabs(czseg2-czz2).gt.0.0001) GOTO 10 ! loop will terminate when last segment has zero length
      endif
c
c     now correct for sign
c
      rfunc=rsign*rfunc
c
      return
c
      END subroutine
c
c ---------------------------------------------------------------------------------
c
      subroutine lkadd_deltapsi (j,ry,rdeltapsi,nrow,ncol,
     &                           rlkdeltax,rlkdeltay,
     &                           rlkleakage,rlkrecharge)
c
c
c
      implicit none
      INTEGER nrow,ncol,i,j,i0,j0
      REAL(8) rlkdeltax,rlkdeltay,rlkleakage,rlkrecharge,rdeltapsi,
     &        rx,ry,rs,rsm,rsp,rytop,rybot,rdx
      DIMENSION rlkdeltax(ncol),rlkdeltay(nrow),rlkleakage(nrow,ncol),
     &          rlkrecharge(nrow,ncol)
      include 'lkcom.inc'
      include 'lusys.inc'
c
c      write (iluer,1001)
c 1001 format (' lkadd_deltapsi1: Entering.')
      rdeltapsi=0.0
      rytop=rlky0-SUM(rlkdeltay(1:nrow)) ! will be used as first bottom y
      if (j.ne.ncol+1) rdx=rlkdeltax(j)
      do i=nrow,1,-1 ! start at bottom of grid and sum line-sink discharges
      rybot=rytop
      rytop=rybot+rlkdeltay(i)
      if (ry.gt.rybot) then ! more of branch cut to be included
       if (j.eq.1) then
        rsm=0.0
        rsp=(rlkleakage(i,1)+rlkrecharge(i,1))*rdx*0.5
       else if (j.eq.ncol+1) then
        rsm=(rlkleakage(i,ncol)+rlkrecharge(i,ncol))*rlkdeltax(ncol)*0.5
        rsp=0.0
       else
        rsm=(rlkleakage(i,j-1)+rlkrecharge(i,j-1))*rlkdeltax(j-1)*0.5
        rsp=(rlkleakage(i,j)+rlkrecharge(i,j))*rdx*0.5
       end if
       rs=rsm+rsp
       if (ry.lt.rytop) then ! delta psi due to that portion of the line-sink below y
        rdeltapsi=rdeltapsi+(ry-rybot)*rs
       else                  ! delta psi due to entire line-sink
        rdeltapsi=rdeltapsi+rlkdeltay(i)*rs ! same as (rytop-rybot) of course
       end if
c       write (iluer,1002) j,i,rdeltapsi
c 1002 format (' lkadd_deltapsi2: j,i,rdeltapsi '2(i3),D14.7)
      end if
      END do
c
      return
      end subroutine
c
c -----------------------------------------------------------------------------------
c
      REAL(8) function rfnflksub (cz1,cz2)
c
c     This is a separate function which must be called separately from "rfnflk"
c     This is so, because it is called separately in DBMAT
c     Note:  ilkresolution may be > 1 only to average upper heads
c            ensure that subgrids are actually being included: lkincludesubgrid=.true.
c
      implicit none
      REAL(8) rfunc
      COMPLEX(8) cz1,cz2
      include 'lkcom.inc'
      include 'lusys.inc'
c
      rfunc=0.0d0
c
      if (lkincludesubgrid) then ! Yes, subgrids are included
       call lknfsub_sub (cz1,cz2,rfunc)
      end if
c
      rfnflksub=rfunc
c
      return
      end
c
c ------------------------------------------------------------------------------
c
      subroutine lkfnfsub_sub_actual (cz1,cz2,rfunc,nrow,ncol,rlkdeltax,
     &                             rlkdeltay,nbuf,rlksubleakage,
     &                             ilkresolution,ilkpointer)
c
c     Routine assigns the flow that crosses the line cz1 - cz2 due to all SUBGRID cells to rfunc
c     This is done in two steps:
c     Step 1 calculate all Laplace contributions of SUBGRID functions as PSI_1 - PSI_2
c     Step 2 add all flows across the line due to the Poisson contributions (SUBGRID leakage)
c     NOTE: Jumps in PSI are avoided by offsetting the intersection points on grid lines.
c

      implicit none
      INTEGER nrow,ncol,i,j,icase,nbuf,ilkresolution,ilkpointer
      LOGICAL lvinsect,lhinsect,linsubgrid
      REAL(8) rlkdeltax,rlkdeltay,rfunc,rtol,
     &        rx,ry,rdx,rdy,rx0,ry0,rx1,ry1,rx2,ry2,rxg1,ryg1,rxg2,ryg2,
     &        rs,rsm,rsp,rsign,rxx1,ryy1,rxx2,ryy2,rdeltapsi,ra,rm,
     &        rmsign,rlksubleakage,rPoisson, d_one
      COMPLEX(8) cz,cz1,cz2,cz0,czz0,czz1,czz2,czseg1,czseg2,
     &           cflk_subomega,czj,czi,coffset
      DIMENSION rlkdeltax(ncol),rlkdeltay(nrow),
     &          ilkresolution(nrow,ncol),
     &          ilkpointer(nrow,ncol),rlksubleakage(nbuf)
      include 'lkcom.inc'
      include 'lusys.inc'
c
c      write (iluer,1001)
c 1001 format (' lkfnfsub_sub_actual1: Entering.')
      rxg1=rlkx0    ! set MODFLOW grid dimensions
      ryg2=rlky0
      rxg2=rxg1+SUM(rlkdeltax(1:ncol))
      ryg1=ryg2-SUM(rlkdeltay(1:nrow))
      rx1=REAL(cz1) ! set line segment end points
      ry1=aimag(cz1)
      rx2=REAL(cz2)
      ry2=aimag(cz2)
      rsign=1.0       ! for now assume line entered from left to right
      if (rx1.gt.rx2) THEN ! reverse points, but remember to change sign on flow across line
       rx1=REAL(cz2)
       ry1=aimag(cz2)
       rx2=REAL(cz1)
       ry2=aimag(cz1)
       rsign=-1.0      ! found that line is entered from right to left, multiplyer for normal flow
      end if
c
      rdx=rxg2-rxg1
c      write (iluer,1011) rdx
c 1011 format (' lkfnfsub_sub_actual11: rdx ',d14.7)
      rdx=rdx/ncol
c      write (iluer,1012) rdx
c 1012 format (' lkfnfsub_sub_actual12: rdx ',d14.7)
      rtol=0.00001*rdx  ! used to offset points
c      write (iluer,1013) rtol
c 1013 format (' lkfnfsub_sub_actual13: rtol ',d14.7)
c
      d_one = 1.0
      czz1=CMPLX(rx1,ry1) ! redefine line segment end points
      czz2=CMPLX(rx2,ry2) ! use local variables so as to not change the input parameters cz1 and cz2
      if (rx1.ne.rx2) then ! calculate slope "rm"
       rm=(ry2-ry1)/(rx2-rx1)
      ELSEIF (ry1.lt.ry2) then ! line segment is vertical
       rm=1.0d30 ! vertically upward
      else
       rm=-1.0d30 ! vertically downward
      end if
      rmsign=SIGN(d_one,rm) ! save sign of rm to be used at end of routine to correct sign of normal flow
      if (rx2.le.rxg1) then
       icase=1 ! line segment entirely to the left of the MODFLOW grid (no subgrid branch cuts)
      elseif (rx1.ge.rxg2) then
       icase=1 ! line segment entirely to the right of the MODFLOW grid (no subgrid branch cuts)
      elseif (ry1.le.ryg1.and.ry2.le.ryg1) then
       icase=1 ! line segment entirely below the grid  (no subgrid branch cuts)
      elseif (ry1.ge.ryg2.and.ry2.ge.ryg2) then
       icase=2 ! line segment entirely above the grid (no subgrid branch cuts, provided no subgrids in
      else     !                                       perimeter MODFLOW cells!)
       icase=3 ! at least part of line segment is inside the grid (may cross subgrid branch cuts
               ! and may need to add Poisson terms for subgrid cells)
      endif
c      write (iluer,1002) rxg1,rxg2,rm,rsign,rtol,icase
c 1002 format (' lkfnfsub_sub_actual2: rxg1,rxg2,rm,rsign,rtol,icase ',
c     &                        5(d14.7),i3)
c
      if (icase.ne.3) THEN ! just add Psi1 - Psi2

c
      rfunc=aimag(cflk_subomega(czz1)-cflk_subomega(czz2)) !                        *** STEP 1 ***
c
c      write (iluer,1003) rfunc
c 1003 format(' lkfnfsub_sub_actual3: rfunc after PSI_1 - PSI_2: ',d14.7)
c
      else  ! (icase=3) may need to add Poisson terms and Psi1 - Psi2 due to subgrids
c
c      First deal with possibility that parts of the line element are outside of the MODFLOW grid.
c      For this case we will redefine the end points on the MODFLOW grid boundary and ignore the outside
c      sections (no branch cuts or Poisson contributions).
c
       if (rx1.lt.rxg1) then ! czz1 to left of the MODFLOW grid
        ry=ry1+rm*(rxg1-rx1)
        if (ryg1.le.ry.and.ryg2.ge.ry) THEN ! intersection point on left-hand MODFLOW grid boundary
         czz0=CMPLX(rxg1-rtol,ry)
         rfunc=rfunc+aimag(cflk_subomega(czz1)-cflk_subomega(czz0))
         czz1=CMPLX(rxg1+rtol,ry)
        end if
       end if
       if (rx2.gt.rxg2) THEN ! czz2 to the right of the grid
        ry=ry1+rm*(rxg2-rx1)
        if (ryg1.lt.ry.and.ryg2.gt.ry) THEN ! intersection point on right-hand MODFLOW grid boundary
         czz0=CMPLX(rxg2+rtol,ry)
         rfunc=rfunc+aimag(cflk_subomega(czz0)-cflk_subomega(czz2))
         czz2=CMPLX(rxg2-rtol,ry)
        end if
       end if
       if (rm.gt.0) then  ! line element slopes upward
        if (ry1.lt.ryg1) then ! czz1 below the grid
         rx=rx1+(ryg1-ry1)/rm
         if (rxg1.lt.rx.and.rxg2.gt.rx) THEN ! intersection point on bottom of MODFLOW grid
          czz0=CMPLX(rx,ryg1-rtol)
         rfunc=rfunc+aimag(cflk_subomega(czz1)-cflk_subomega(czz0))
          czz1=CMPLX(rx,ryg1+rtol)
         end if
        end if
        if (ry2.gt.ryg2) THEN ! czz2 above the grid
         rx=rx1+(ryg2-ry1)/rm
         if (rxg1.lt.rx.and.rxg2.gt.rx) THEN ! intersection point on top of MODFLOW grid
          czz0=CMPLX(rx,ryg2+rtol)
         rfunc=rfunc+aimag(cflk_subomega(czz0)-cflk_subomega(czz2))
          czz2=CMPLX(rx,ryg2-rtol)
         end if
        end if
       else  ! (rm.lt.0)   line element slopes downward
        if (ry1.gt.ryg2) THEN ! czz1 above the grid
         rx=rx1+(ryg2-ry1)/rm
         if (rxg1.lt.rx.and.rxg2.gt.rx) THEN ! intersection point on top of MODFLOW grid
          czz0=CMPLX(rx,ryg2+rtol)
         rfunc=rfunc+aimag(cflk_subomega(czz1)-cflk_subomega(czz0))
          czz1=CMPLX(rx,ryg2-rtol)
         end if
        end if
        if (ry2.lt.ryg1) THEN ! czz2 below the grid
         rx=rx1+(ryg1-ry1)/rm
         if (rxg1.lt.rx.and.rxg2.gt.rx) THEN ! intersection point on bottom of MODFLOW grid
          czz0=CMPLX(rx,ryg1-rtol)
         rfunc=rfunc+aimag(cflk_subomega(czz0)-cflk_subomega(czz2))
          czz2=CMPLX(rx,ryg1+rtol)
         end if
        end if
       end if
       rx1=REAL(czz1)   ! for the purpose of adding the Poisson contributions and branch cuts
       ry1=AIMAG(czz1)  ! the line element is chopped to fit entirely inside the MODFLOW grid
       rx2=REAL(czz2)   ! (no points on the MODFLOW grid boundary; used offset)
       ry2=AIMAG(czz2)
c      write (iluer,1005) czz1,czz2
c 1005 format(' lkfnfsub_sub_actual5: czz1,czz2 after clipping to grid ',
c     &         4(d14.7))
       czseg1=czz1
c
c     Now break up the line element into sections czseg1-czseg2 that fit inside a single MODFLOW cell
c
  10  continue
      lvinsect=.false.
      lhinsect=.false.
c      Find the nearest intersection with a vertical MODFLOW grid line
       rx=rxg2-rlkdeltax(ncol)
       do j=ncol,2,-1 ! loop left to right, so we end with the nearest intersection
       if (rx1.lt.rx.and.rx2.gt.rx) THEN ! yes, we are intersecting vertical MODFLOW gridline j
        ry=ry1+rm*(rx-rx1) ! get y-value of intersect
        czj=CMPLX(rx,ry) ! make new tentative end point
        lvinsect=.true.
       endif
       rx=rx-rlkdeltax(j)
       enddo
c      Find nearest intersection with a horizontal MODFLOW grid line
       if (rm.gt.0) THEN ! line slopes upward
       ry=ryg2-rlkdeltay(1)
       do i=2,nrow  ! loop top to bottom, so we end with nearest intersection
        if (ry1.lt.ry.and.ry2.gt.ry) THEN ! we are intersecting horizontal MODFLOW gridline i
         rx=rx1+(ry-ry1)/rm ! get x-value of intersect
         czi=CMPLX(rx,ry) ! make new tentative end point
         lhinsect=.true.
        end if
        ry=ry-rlkdeltay(i)
       enddo
       else ! (rm.lt.0)  ! line slopes downward
       ry=ryg1+rlkdeltay(nrow)
       do i=nrow,2,-1 ! loop from bottom to top, so we end up with nearest intersection
        if (ry1.gt.ry.and.ry2.lt.ry) THEN ! we are intersecting horizontal MODFLOW gridline i
         rx=rx1+(ry-ry1)/rm ! get x-value of intersect
         czi=CMPLX(rx,ry) ! make new tentative end point
         lhinsect=.true.
        end if
        ry=ry+rlkdeltay(i)
       enddo
       end if
      if (lvinsect.and.lhinsect) then ! two new potential end points, select closest to czz1
c      now which is closest?
       if (cdabs(czj-czz1).lt.cdabs(czi-czz1)) then
        czseg2=czj
        coffset=CMPLX(rtol,0.0)
       else
        czseg2=czi
        coffset=rmsign*CMPLX(0.0,rtol)
       end if
      else if (lvinsect) then ! select point on vertical MODFLOW grid line
        czseg2=czj
        coffset=CMPLX(rtol,0.0)
      else if (lhinsect) then ! select point on horizontal MODFLOW grid line
        czseg2=czi
        coffset=rmsign*CMPLX(0.0,rtol)
      else                    ! keep original end point (no intersect)
        czseg2=czz2
        coffset=CMPLX(0.0,0.0)
      end if
c      write (iluer,1006) lvinsect,lhinsect,czseg1,czseg2
c 1006 format (' lkfnfsub_sub_actual6: lvinsect,lhinsect,czseg1,czseg2 ',
c     &           2(l3),4(d14.7))
c
c      now see if me must add subgrid contributions for the
c      line segment czseg1 - czseg2 inside MODFLOW cell i,j
c
      cz0=0.5*(czseg1+czseg2)
      call lkfindcell_actual (cz0,i,j,rx0,ry0,   ! returns MODFLOW cell i,j
     &                        rlkdeltax,rlkdeltay,nrow,ncol)
      linsubgrid=ilkresolution(i,j).gt.1   ! contains a subgrid (Poisson contributions and branch cut)
c
c      write (iluer,1007) linsubgrid
c 1007 format (' lkfnfsub_sub_actual7: linsubgrid ',l3)
c
      if (linsubgrid) then
        czz0=czseg2-coffset
        call lkadd_subPsiPoisson (czseg1,czz0,rm,i,j, ! include Poisson terms        *** STEP 3 ***
     &                          nrow,ncol,rlkdeltax,rlkdeltay,
     &                          rfunc,nbuf,rtol,
     &                          ilkresolution,ilkpointer,rlksubleakage)
      end if
      rx1=REAL(czseg2) ! drop segment just considered
      ry1=AIMAG(czseg2)
      czseg1=czseg2+coffset
      if (cdabs(czseg2-czz2).gt.0.0001) GOTO 10 ! loop will terminate when last segment has zero length
      endif
c
c     now correct for sign
c
      rfunc=rsign*rfunc
c
c      write (iluer,1008) rsign,rfunc
c 1008 format (' lkfnfsub_sub_actual8: leaving with rsign and rfunc ',
c     &             2(d14.7))
      return
c
      END subroutine
c
c -----------------------------------------------------------------------------------------
c
      subroutine lkadd_subPsiPoisson  (cz1,cz2,rm,irow,jcol,
     &                          nrow,ncol,rlkdeltax,rlkdeltay,
     &                          rfunc,nbuf,rtol,
     &                          ilkresolution,ilkpointer,rlksubleakage)
c
c     Add Psi1 - Psi2 and Poisson terms for subgrid in MODFLOW cell i,j
c     which the segment czseg1 - czseg2 intersects.
c
c
      implicit none
      INTEGER nrow,ncol,i,j,i0,j0,istart,jstart,irow,jcol,nbuf,
     &        ilkresolution,ilkpointer,ires
      LOGICAL lvinsect,lhinsect
      REAL(8) rlkdeltax,rlkdeltay,rlksubleakage,rfunc,
     &        rx,ry,rx0,ry0,rx1,ry1,rx2,ry2,rxx1,ryy1,rxx2,ryy2,
     &        rxcenter,rdx,rdy,rdxx,rdyy,ra,rtol,rPoisson,rpsi,
     &        rs,rytop,rybot,rm,rmsign, d_one
      COMPLEX(8) cz1,cz2,czseg1,czseg2,czj,czi,coffset,czz0,
     &           cflk_subomega
      DIMENSION rlkdeltax(ncol),rlkdeltay(nrow),rlksubleakage(nbuf),
     &          ilkresolution(nrow,ncol),ilkpointer(nrow,ncol)
      include 'lkcom.inc'
      include 'lusys.inc'
c
c      write (iluer,1001) irow,jcol,rtol
c 1001 format (' lkadd_subPoisson1: entering with irow, jcol, rtol = ',
c     &          2(i3),d14.7)
      istart=ilkpointer(irow,jcol)
      czseg1=cz1
      czseg2=cz2
      rx1=REAL(czseg1)
      ry1=AIMAG(czseg1)
      rx2=REAL(czseg2)
      ry2=AIMAG(czseg2)
      ires=ilkresolution(irow,jcol)
      rdx=rlkdeltax(jcol)
      rdy=rlkdeltay(irow)
      rx0=rlkx0+SUM(rlkdeltax(1:jcol-1)) ! x-value of left grid line of MODFLOW cell
      ry0=rlky0-SUM(rlkdeltay(1:irow))  ! y-value of lower grid line of MODFLOW cell
      rdxx=rdx/ires
      rdyy=rdy/ires
      d_one = 1.0
      rmsign=SIGN(d_one,rm)
c      write (iluer,1002) rx0,ry0
c 1002 format (' lkadd_subPoisson2: rx0,ry0 ',2(d14.7))
c      write (iluer,1003) istart,
c     &                 rlksubleakage(istart),rlksubleakage(istart+1),
c     &                 rlksubleakage(istart+2),rlksubleakage(istart+3)
c 1003 format(' lkadd_subPoisson3:istart,rlksubleakage(istart-istart+3)',
c     &         /,i3,4(d14.7)                    )
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
c
  10  continue
      lvinsect=.false.
      lhinsect=.false.
c     find the nearest intersection with a vertical subgrid line
      rx=rx0+rdx
      do j=ires,1,-1 ! loop left to right so we end at the nearest intersection
       if (rx1.lt.rx) THEN ! we may have found subcell inwhich czseg1 occurs
        rxcenter=rx-0.5*rdxx ! needed for Poisson contribution calculations
        j0=j   ! needed to calculate address for subcell leakage
       end if
       if (rx1.lt.rx.and.rx2.gt.rx) THEN ! yes, we intersect vertical subgrid line j
        ry=ry1+rm*(rx-rx1) ! get y-value of intersect
        czj=CMPLX(rx,ry) ! make new tentative end point
        lvinsect=.true.
       end if
       rx=rx-rdxx
      end do
c     find the nearest intersection with a horizontal subgrid line
      if (rm.gt.0) THEN ! line slopes upward
       ry=ry0+rdy ! start at top of upper most subgrid cell
       do i=ires,1,-1 ! loop from top to bottom so we end with the nearest intersection
       if (ry1.lt.ry) then ! we may have found the subcell in which czseg1 occurs
        i0=i ! needed to calculate address for subcell leakage
       end if
       if (ry1.lt.ry.and.ry2.gt.ry) then ! yes, we are intersecting horizontal grid line i
        rx=rx1+(ry-ry1)/rm ! get x-value of end point
        czi=CMPLX(rx,ry)
        lhinsect=.true.
       end if
       ry=ry-rdyy
       end do
      else ! (rm.lt.0) line slopes downward
       ry=ry0+rdyy ! start at top of lower most subgrid cell
       do i=1,ires ! loop from bottom to top so we end with the nearest intersection
       if (ry1.gt.ry) then ! we may have found the subcell inwhich czseg1 occurs
        i0=i   ! needed to calculate address for subcell leakage
       end if
       if (ry1.gt.ry.and.ry2.lt.ry) THEN ! we are intersecting horizontal line i
        rx=rx1+(ry-ry1)/rm ! get x-value of intersect
        czi=CMPLX(rx,ry)
        lhinsect=.true.
       end if
       ry=ry+rdyy
       end do
      end if
      if (lvinsect.and.lhinsect) then ! two new potential end points, select closest to czseg1
c      now which one is closest ?
       if (cdabs(czj-czseg1).lt.cdabs(czi-czseg1)) THEN ! czj wins
         czseg2=czj
         coffset=CMPLX(rtol,0.0)
       else                                             ! czi wins
         czseg2=czi
         coffset=rmsign*CMPLX(0.0,rtol)
       end if
      elseif (lvinsect) THEN ! select point on vertical grid
       czseg2=czj
       coffset=CMPLX(rtol,0.0)
      elseif (lhinsect) THEN ! select point on horizontal grid
       czseg2=czi
       coffset=rmsign*CMPLX(0.0,rtol)
      else
       czseg2=cz2
       coffset=CMPLX(0.0,0.0)
      end if
c
         czz0=czseg2-coffset
         rpsi=aimag(cflk_subomega(czseg1)-cflk_subomega(czz0))
         rfunc=rfunc+rpsi
c      write (iluer,1004) czseg1,czz0,rpsi
c 1004 format (' lkadd_subPoisson4: czseg1,czz0,rpsi ',/,5(d14.7))
c
      jstart=istart+(i0-1)*ires+j0-1 ! find the subleakage address
      rs=rlksubleakage(jstart)
      rxx1=rx1
      ryy1=ry1
      rxx2=REAL(czseg2)
      ryy2=AIMAG(czseg2)
      ra=0.5*((rxx1-rxcenter)+(rxx2-rxcenter))*ABS(ryy2-ryy1)
      rPoisson=-rmsign*rs*ra
      rfunc=rfunc+rPoisson
c      write (iluer,1005) i0,j0,irow,jcol,rxcenter,rxx1,rxx2,ryy1,ryy2,
c     &                   rs,ra
c 1005 format (' lkadd_subPoisson5: i0,j0,irow,jrow,rxcenter ',4(i3),
c     &d14.7,/,' rxx1,rxx2,ryy1,ryy2,rs,ra ',6(d14.7))
c
      rx1=REAL(czseg2)
      ry1=AIMAG(czseg2)
      czseg1=CMPLX(rx1,ry1)+coffset
      if (cdabs(czseg2-cz2).gt.0.0001) GOTO 10 ! loop will terminate when last segment has zero length
c
c      write (iluer,1006) rPoisson
c 1006 format (' lkadd_subPoisson6: leaving with rPoisson=',d14.7)
      return
c
      end subroutine

