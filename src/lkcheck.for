C     Last change:  HMH  15 Nov 2010    8:09 pm
c
c  lkset_check      calculate the leakage in a grid of points inside a MODFLOW cell and sum over all cells
c  lkerror_actual          reports error at collocation points
c
c --------------------------------------------------------------------
c
      subroutine lkset_check (i1,j1,i2,j2,ires,
     & nrow,ncol,
     & rlkdeltax,rlkdeltay,rlkleakage,rlkheadlower,
     & rlkresist,ilkresolution,lwrite)
c
c     calculate total leakage in a subset of the MODFLOW cells and calculate the % error.
c     The total leakage is calculated in two ways:
c      1) add all leakage parameters times their cell sizes.
c      2) calculate the leakage from the upper and lower heads and the resistance in a grid
c         of points inside each cell and multiply each by its representative area.
c     The % error is defined as the error between the sink density and the sum of the leakage
c     calculated at a grid of points in the MODFLOW cell.
c
c
c     Input:
c             i1,j1 starting row and column (upper left corner of grid)
c             i2,j2 ending row and column (lower right corner of grid)
c             ires  resolution of leakage computations
c             lwrite   true if the %error is to be written to the file basename.ler
c             Note: when lwrite=.true. the calculations are performed for the entire MODFLOW grid.
c
c
c
      implicit none
      INTEGER i1,j1,i2,j2,ires,nrow,ncol,i,j,ii,jj,
     &        nsolOut,ilu,ierr,ilkresolution,icode
      LOGICAL lwrite,
     &        lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut
      REAL(8) rlkdeltax,rlkdeltay,rlkleakage,rlkheadlower,rlkresist,
     &        rx1,ry1,rxx1,ryy1,rdx,rdy,rdxx,rdyy,rx0,rq,rdq,resist,
     &        rheadlower,rheadupper,rfhead,res,scratch,rleak,rsink,
     &        rbuf
      COMPLEX(8) cz,cz0,cz1,cz2,cz3
      CHARACTER*8 aBasenameOut
      CHARACTER*16 aDateTimeOut
      ALLOCATABLE scratch(:),rbuf(:,:)
      DIMENSION rlkdeltax(ncol),rlkdeltay(nrow),
     & rlkleakage(nrow,ncol),rlkheadlower(nrow,ncol),
     & rlkresist(nrow,ncol),ilkresolution(nrow,ncol)
      include 'lkcom.inc'
      include 'lusys.inc'
      include 'match.inc'
c
c ------------ check data integrity
c
      lkcheck=.true.
      if (lwrite) then ! use entire grid if basename.ler is to be written.
       i1=1
       i2=nrow
       j1=1
       j2=ncol
      else
       if (i1.lt.1.or.i1.gt.nrow) lkcheck=.false.
       if (j1.lt.1.or.j1.gt.ncol) lkcheck=.false.
       if (i2.lt.i1.or.j2.lt.j1) lkcheck=.false.
       if (.not.lkcheck) then
        write (iluer,1000) i1,j1,i2,j2
       end if
      end if
      if (ires.lt.1.or.ires.gt.1000) lkcheck=.false.
      if (.not.lkcheck) then
        write (iluer,2000) ires
      end if
c
      if (lkcheck.and.lwrite) then  ! open file for % errors
        ilu=2
       call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
       afile=TRIM(aBasenameOut)//".ler"
       OPEN (UNIT=ILU,FILE=AFILE,STATUS='UNKNOWN',IOSTAT=IERR)
       if (ierr.ne.0) then  ! failed to open file
        write (iluer,3000) afile
        return
       endif
      end if
      if (lwrite) then  ! allocate a scratch array to write the percent errors
       ALLOCATE (scratch(ncol))
      end if
c
c      write (iluer,1001) i1,i2,j1,j2,ires,lkcheck,lwrite
c 1001 format (' lkset_check1: i1,i2,j1,j2,ires,lcheck,lwrite ',/,
c     & 5(I4),2x,2(l4))
      if (lkcheck) then
      if (i1.eq.1) then
       ry1=rlky0
      else
       ry1=rlky0-SUM(rlkdeltay(1:i1-1))
      end if
      if (j1.eq.1) then
       rx1=rlkx0
      else
       rx1=rlkx0+SUM(rlkdeltax(1:j1-1))
      end if
      rx0=rx1
      rQ_leakage=0.0d0  ! calculate the total leakage due to the MODFLOW grid cells.
      rQ_target_leakage=0.0d0
      do i=i1,i2
       rdy=rlkdeltay(i)
       do j=j1,j2
        rdx=rlkdeltax(j)
        rsink=rlkleakage(i,j)
        rQ_leakage=rQ_leakage+rsink*rdx*rdy
c       ! now do the grid inside this cell
        res=float(ires)
        rdxx=rdx/res
        rdyy=rdy/res
        rxx1=rx1+0.5*rdxx
        ryy1=ry1-0.5*rdyy
c        rheadlower=rlkheadlower(i,j) ! NOTE: currently only a constant lower head will work!!!
c
      if (ALLOCATED(rbuf)) deallocate (rbuf)
      ALLOCATE (rbuf(ires,ires),stat=ierr)
      if (ierr.ne.0) then
       write (iluer,5000) ires,ires
        AMESS(1)='Error in routine lkset_check.'
        AMESS(2)='Failed to allocate rbuf(ires,ires).'
        AMESS(3)='Verify resources and rerun.'
        CALL HALT(3) ! stop program execution for batch version
      endif
c
c
        ! Note: we are forcing a uniform grid in x- and y-direction
        call lk_cell_construct (1,1,rdx,rdy,cz0,cz1,cz2,cz3)
        call lkinterpolateheads(nrow,ncol,i,j,rdx,rdy,rlkheadlower,ires,
     &                               rbuf,icode,iluer)
c
        resist=rlkresist(i,j)
        rq=0.0d0
c        write (iluer,1002) i,j,rsink,res,rdx,rdy,rdxx,rdyy,rxx1,ryy1
c 1002 format (' lkset_check2: i,j,rsink / ',
c     &  'res,rdx,rdy,rdxx,rdyy,rxx1,ryy1 ',2(I5),2x,d14.7,/,6(d14.7))
        do jj=1,ires
         do ii=1,ires
          cz=CMPLX(rxx1,ryy1)
          rheadupper=rfhead(cz)
          rheadlower=rbuf(ii,jj)
          rdq=(rheadupper-rheadlower)/resist
          rq=rq+rdq
          write (iluer,1003)i,j,cz,rheadupper,rheadlower,resist,rdq,rq
 1003 format (' lkset_check3: '
     &' i,j,cz,rheadupper,rheadlower,resist,rdq,rq',/,
     &  2(i4,1x),7(d14.7))
          ryy1=ryy1-rdyy
         end do
         ryy1=ry1-0.5*rdyy
         rxx1=rxx1+rdxx
        end do
        rQ_target_leakage=rQ_target_leakage+rq*rdxx*rdyy
        if (lwrite) then
         rleak=rq/res/res
         scratch(j)=100*ABS(rleak-rsink)/(0.5*(ABS(rleak)+ABS(rsink))) ! % error in leakage
c         if (i.eq.1.and.j.eq.1) then
c         write (iluer,1004) ires,res,rq,rleak,rsink,scratch(j)
c 1004    format (' lkset_check4: ires,res,rq,rleak,rsink,scratch ',/,
c     &         i5,2x,5(d14.7))
c         END if
        end if
        if (ilkresolution(i,j).gt.0.and.rleak*rsink.lt.0.0d0) then ! opposite sign, report
          write (iluer,4000) i,j,rsink,rleak
        end if
        rx1=rx1+rdx
       end do
       rx1=rx0
       ry1=ry1-rdy
       if (lwrite) then
         write (ilu,*) scratch
       end if
      end do
      end if
c
      if (lwrite) then
       close (ilu)
       deallocate (scratch)
      end if
      return
 1000 format (' ***ERROR in lkset_check: Illegal array parameters.',/,
     &' i1, j1, i2, j2 = ',4(i5,1x))
 2000 format (' ***ERROR in lkset_check: resolution out of range.',/,
     &' ires = ',i5,', but must be between 1 and 1000.')
 3000 format (' ERROR in lkset_check:',/,
     &        ' The file ',a16,' is not found.',/,
     &        ' No % errors are written.')
 4000 format (' Warning in LKSET_CHECK: leakage rate for cell(',i3,
     &',',i3,') has opposite sign of calculated leakage.',/,
     &' cell=',d14.7,' calculated=',d14.7)
 5000 format ('***ERROR in LKSET_CHECK: failed to allocate rbuf('
     &,i5,',',i5,').')
      end subroutine
C
C---------------------------------------------------------------------------------------------------------
C
      SUBROUTINE lkerror_actual (rermax,
     &     nrow,ncol,rlkpot,rlkheadlower,rlkheadupper,rlkleakage,
     &     rlkresist,rlkdeltax,rlkdeltay,ilkresolution,rlkpercenterror,
     &     rlksubpotupper,rlksubheadlower,ilkpointer,rlksubleakage,nbuf)
C
C---------------------------------------------------------------------------------------------------------
c
c     update the potential at collocation points of leakage elements
c
c
      IMPLICIT NONE
      INTEGER(4) i,j,ii,jj,nrow,ncol,n,nbuf,iad,ires0,ires2,itr,
     &           irow,jcol,ncell,ilkresolution,iresolution,iad0,ierr,
     &           ilkpointer,istart
      LOGICAL lwrite,l1,lflkincludesubgrid
      CHARACTER(256) amessage
      REAL(8) rlkresist,rlkpot,rdx,rdy,rx,ry,rx1,ry1,rfhedp,
     &        rlkleakage,rlkheadlower,rpotlower,rlkerror,rermax,
     &        rhupper,rhlower,rc,rs,rcalculatedleak,rpercenterror,
     &        rlkdeltax,rlkdeltay,rconverge_leakage,rlkheadupper,
     &        rlksubpotupper,rlksubheadlower,rlkpercenterror,
     &        rdxx,rdyy,resolution,rcalculatedleak_sum,rtrans,rt,rb,rk,
     &        rtransaverage,rh,rh0,rftop,rfbase,rfperm,rcloc,
     &        rlksubleakage,rs0
      REAL(8) rhupperc,rhlowerc,rhupper_temp,rdeltahupper
      COMPLEX(8) CZ,cz0,cz1,cz2,cz3
      DIMENSION rlkresist(nrow,ncol),rlkheadupper(nrow,ncol),
     &          rlkpot(nrow,ncol),rlkleakage(nrow,ncol),
     &          rlkheadlower(nrow,ncol),rlkdeltax(ncol),
     &          rlkdeltay(nrow),ilkresolution(nrow,ncol),
     &          rlkpercenterror(nrow,ncol),ilkpointer(nrow,ncol),
     &          rlksubpotupper(nbuf),rlksubheadlower(nbuf),
     &          rlksubleakage(nbuf)
      ALLOCATABLE rtrans(:)
      INCLUDE 'LKCOM.INC'
      INCLUDE 'LUSYS.INC'
      include 'TRACOM.INC'
      include 'match.inc'
      save
c
      DATA rconverge_leakage /1.0d0/  ! hardwire convergence criteria to 1% (NOTE: must be added to converge file)
c
      ires0=0
      rlkerror=0.0
      ncell=0
      iad=0
      if (lkleakage.and.lksolving) then
c      rx1=rlkx0
c      ry1=rlky0
      do i=1,nrow
c      rdy=rlkdeltay(i)
c      ry=ry1-0.5*rdy
       do j=1,ncol
c       rdx=rlkdeltax(j)
c       rx=rx1+0.5*rdx
c      lwrite=(i.eq.20.and.j.eq.20).OR.(i.eq.3.and.j.eq.3)   !debugging
       lwrite=.false.
c
        rhupperc=rlkheadupper(i,j)  !  upper head at cell center
        rhlowerc=rlkheadlower(i,j)  !  lower head at cell center
        rc=rlkresist(i,j)           !  resistance between GFLOW and MODFLOW layer
        rs0=rlkleakage(i,j)         !  average leakage rate for MODFLOW cell
        rhupper_temp=rhlowerc+rc*rs0
        rdeltahupper=rhupper_temp-rhupperc
        rlkheadupper(i,j)=rhupperc+rlkrelax*rdeltahupper ! Warning: rlkheadupper(i,j) now contains the NEW upperhead for MODFLOW
c                                                                  Do not access again below!!
        iresolution=ilkresolution(i,j)
        select case (iresolution)
c
            case (0)                                ! subgrid resolution = 0   (no subcells)
c             No further action needed.
c             rc=rlkresist(i,j)
c             rs=rlkleakage(i,j)
c             rlkheadupper(i,j)=rlkheadlower(i,j)+rs*rc ! calculate a representative upper head
c
            case (1)
             ncell=ncell+1                          ! subgrid resolution = 1   (no subcells)
c             cz=CMPLX(rx,ry)
c             rhupper=RFHEDP(rlkpot(i,j),CZ) !old statement, see next
c             rhupper=rlkheadupper(i,j)
c             rdum=rfhead(cz) ! debugging, delete
c             rhlower=rlkheadlower(i,j)
c             rc=rlkresist(i,j)
c             rs=rlkleakage(i,j)
             rcalculatedleak=(rhupperc-rhlowerc)/rc
             rpercenterror=ABS(rcalculatedleak-rs0)
             rpercenterror=rpercenterror/(0.5*(ABS(rcalculatedleak)+
     &                                                      ABS(rs0)))
             rpercenterror=rpercenterror*100.0
            if (rs0*rcalculatedleak.gt.0.0d0) then ! avoid 200% errors from opposite signs
               rlkerror=MAX(rlkerror,rpercenterror)
            else
              write (iluer,2000) i,j,rs0,rcalculatedleak
            end if
      if (lwrite)
     &write (iluer,1001) i,j,rhupperc,rhlowerc,rs0,rc,rcalculatedleak
 1001 format (' lkerror1: rhupperc,rhlowerc,rs0,rc,rcalculatedleak'
     &         ,/,2i5,5(d14.7))
c
c                                                           ---------- case 2:
c
            case (2:)                                       ! subgrid resolution => 2    (subcells present)
              ncell=ncell+1
              if (lflkincludesubgrid()) THEN ! calculate conditions at sub-cell centers
c               rc=rlkresist(i,j)
c               rs0=rlkleakage(i,j)
               iresolution=ilkresolution(i,j)
               istart=ilkpointer(i,j)
               call lk_cell_construct(i,j,rdx,rdy,cz0,cz1,cz2,cz3)
               resolution=real(iresolution)
               rdxx=rdx/resolution
               rdyy=rdy/resolution
               rx=0.5*rdxx
               ry=0.5*rdyy
               cz=cz0+CMPLX(rx,ry)
               rcalculatedleak=0.0d0
               rs=0.0d0
               do ii=1,iresolution
                do jj=1,iresolution
                iad=istart+(ii-1)*iresolution+jj-1
                rhupper=rfhedp(rlksubpotupper(iad),cz)  ! local upper head at subcell center
                rhlower=rlksubheadlower(iad)            ! local lower head at subcell center (from interpolator)
c                rcalculatedleak=rcalculatedleak+(rhupper-rhlower)/rc
                rcalculatedleak=(rhupper-rhlower)/rc     ! ??? strange, left over from editing??
                rs=rs0+rlksubleakage(iad)   ! rs0 contains average leakage for the MODFLOW cell and
c                                             rlksubleakage(iad) contains the deviation of the average leakage
c                                             for the subcell (ii,jj).
                rpercenterror=ABS(rcalculatedleak-rs)
                rpercenterror=rpercenterror/(0.5*(ABS(rcalculatedleak)+
     &                                                      ABS(rs)))
                rpercenterror=rpercenterror*100.0
                if (rs*rcalculatedleak.gt.0.0d0) then ! avoid 200% errors from opposite signs
                 rlkerror=MAX(rlkerror,rpercenterror)  ! For unconfined flow this will never be accurate!!
                else
                 write (iluer,2000) i,j,rs,rcalculatedleak
                end if
c
                cz=cz+CMPLX(rdxx,0.0d0)
                end do
                rx=0.5*rdxx
                ry=ry+rdyy
                cz=cz0+CMPLX(rx,ry)
               end do
C               rlkheadupper(i,j)=rlkheadlower(i,j)+rs0*rc ! calculate a representative upper head     >>>>>>>>> replaced "rs" by "rs0" !!!!!
c
              else                       ! subgrid not included, consider average conditions or condition at MODFLOW cell center
              if (lkaveragehead) then ! average heads over sub-grid centers of MODFLOW cell
c              rc=rlkresist(i,j)
c              rs=rlkleakage(i,j)
              call lk_cell_construct(i,j,rdx,rdy,cz0,cz1,cz2,cz3)
              resolution=real(iresolution)
              rdxx=rdx/resolution
              rdyy=rdy/resolution
              if (iresolution.gt.ires0) then ! allocate or reallocate array "rtrans"
                ires0=iresolution            ! make sure that array "rtrans" is big enough
                l1=ALLOCATED(rtrans)
                if (l1) DEALLOCATE(rtrans)
                ires2=ires0*ires0
                ALLOCATE (rtrans(ires2),stat=ierr)
                if (ierr.ne.0) then
                  call iostat_msg (ierr,amessage)
                  write (ilume,1112) amessage
 1112             format (a132)
                  deallocate (rtrans)
                  write (ilume,8001) ires2
 8001     format ('***ERROR in LKERROR: cannot allocate rtrans(',i5,')')
                  AMESS(1)='Error in LKERROR routine.'
                  AMESS(2)='Failed to allocate a scratch array.'
                  AMESS(3)='Execution has been aborted.'
                  CALL HALT(3)   ! stop program execution for batch version
              end if
              end if
c
c            generate T_i and T
c
              rx=0.5*rdxx
              ry=0.5*rdyy
              cz=cz0+CMPLX(rx,ry)
              itr=0
              iad0=iad
              rtransaverage=0.0d0
              do ii=1,iresolution
               do jj=1,iresolution
                itr=itr+1
                iad=iad+1 ! next sub-cell address
                rt=rftop(cz)
                rb=rfbase(cz)
                rk=rfperm(cz)
                rh0=rt-rb
                rhupper=rfhedp(rlksubpotupper(iad),cz)  ! upper head at subcell center
                rhlower=rlksubheadlower(iad)            ! lower head at subcell center (from interpolator)
                if (rhupper.lt.rt) then
                 rh=0.5*(rhupper+rhlower)-rb
                else
                 rh=rh0
                end if
                rtrans(itr)=rk*rh
                rtransaverage=rtransaverage+rtrans(itr)
      if (lwrite) write (iluer,1002)
     &cz,rhupper,rhlower,rh,rk,rtrans(itr),rtransaverage
 1002 format (
     &' lkerror2: cz,rhupper,rhlower,rh,rk,rtrans(itr),rtransaverage'
     &         ,/,2(d14.7),2x,6(d14.7))
               cz=cz+CMPLX(rdxx,0.0d0)
               end do
               rx=0.5*rdxx
               ry=ry+rdyy
               cz=cz0+CMPLX(rx,ry)
              end do
              rtransaverage=rtransaverage/ires2
c
c            generate c_i, L_i, and L
c
              rx=0.5*rdxx
              ry=0.5*rdyy
              cz=cz0+CMPLX(rx,ry)
              itr=0
              iad=iad0
              rcalculatedleak_sum=0.0d0
              do ii=1,iresolution
               do jj=1,iresolution
                iad=iad+1 ! next sub-cell address
                itr=itr+1
                rhupper=rfhedp(rlksubpotupper(iad),cz)   ! local upper head at subcell center
                rhlower=rlksubheadlower(iad)             ! local lower head at subcell center (from interpolator)
c                rcloc=rc*rtransaverage/rtrans(itr)      ! blocked, we used correct condition at sub-cell centers
                rcloc=rc ! may be varied in the future
                rcalculatedleak=(rhupper-rhlower)/rcloc
                rcalculatedleak_sum=rcalculatedleak_sum+rcalculatedleak
      if (lwrite)
     &write (iluer,1003) cz,rhupper,rhlower,rs,rcloc,rcalculatedleak
 1003 format (' lkerror3: cz,rhupper,rhlower,rs,rc,rcalculatedleak'
     &         ,/,7(d14.7))
               cz=cz+CMPLX(rdxx,0.0d0)
               end do
               rx=0.5*rdxx
               ry=ry+rdyy
               cz=cz0+CMPLX(rx,ry)
              end do
              rcalculatedleak=rcalculatedleak_sum/ires2
              rpercenterror=ABS(rcalculatedleak-rs)
              rpercenterror=rpercenterror/(0.5*(ABS(rcalculatedleak)+
     &                                                      ABS(rs)))
              rpercenterror=rpercenterror*100.0
      if (lwrite)
     &write (iluer,1004) i,j,rcalculatedleak,rpercenterror
 1004 format (' lkerror4: i,j,rcalculatedleak (average), rpercenterror '
     &         ,2i5,2(d14.7))
              if (rs*rcalculatedleak.gt.0.0d0) then ! avoid 200% errors from opposite signs
                rlkerror=MAX(rlkerror,rpercenterror)  ! For unconfined flow this will never be accurate!!
              else
                write (iluer,2000) i,j,rs,rcalculatedleak
              end if
c
c              rlkheadupper(i,j)=rlkheadlower(i,j)+rs*rc ! calculate a representative upper head  >> here "rs" is MODFLOW cell leakage
c      write (iluer,1005) i,j,rlkheadupper(i,j),rlkheadlower(i,j),rs,rc
c 1005 format ('lkerror5: i,j,upperhead,lowerhead,rs,rc ',2i5,4(d14.7))
c              Note: When using upper heads as input for MODFLOW, we want the upperhead, lower head and
c                    leakage to be consistent with each other.
c
          else                                     ! sub-cells, but no average head specified
c             rhupper=rlkheadupper(i,j)             ! treat as case 1
c             rdum=rfhead(cz) ! debugging, delete
c             rhlower=rlkheadlower(i,j)
c             rc=rlkresist(i,j)
c             rs=rlkleakage(i,j)
             rcalculatedleak=(rhupperc-rhlowerc)/rc
             rpercenterror=ABS(rcalculatedleak-rs0)
             rpercenterror=rpercenterror/(0.5*(ABS(rcalculatedleak)+
     &                                                      ABS(rs0)))
             rpercenterror=rpercenterror*100.0
             if (rs0*rcalculatedleak.gt.0.0d0) then ! avoid 200% errors from opposite signs
               rlkerror=MAX(rlkerror,rpercenterror)
            else
              write (iluer,2000) i,j,rs0,rcalculatedleak
            end if
c
          endif
          endif
c
            case default
               write (iluer,999) i,j,iresolution
 999  format (' ***ERROR in LKERROR_ACTUAL:',/,
     & 'iresolution(',i3,','i3,')=',i6,' which is invalid.',/,
     & 'Program execution has been aborted.')
               AMESS(1)='Subgrid resolution out of range.'
               AMESS(2)='Stored in file: "basename.grs"'
               AMESS(3)='Correct and rerun.'
               CALL HALT(3) ! stop program execution for batch version
        end select
c
c       store the percent error
c
        rlkpercenterror(i,j)=rpercenterror
c
        rx1=rx1+rdx
       end do
       rx1=rlkx0
       ry1=ry1-rdy
      end do
      if (.not.lucon) WRITE (ILUME,1000) ncell,rlkerror
      write (*,1000) ncell,rlkerror
      RERMAX=MAX(RERMAX,rlkerror)
      endif
      if (rlkerror.ge.rconverge_leakage) lquit=.FALSE. ! do not yet abort iterations
      RETURN
 1000 FORMAT (' ',I4,' leakage cells in matrix:           max. error=',
     &        E11.4,' %')
 2000 format (' Warning in LKERROR: leakage rate for cell(',i3,',',i3,
     &') has opposite sign of calculated leakage.',/,
     &' cell=',d14.7,' calculated=',d14.7)
      END

