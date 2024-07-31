C     Last change:  HMH  16 Apr 2015    3:11 am
c
c     lkmatsize_actual       determines the number of control points and strength parameters for the leakage grid.
c     lkczc_actual           collocation point generator
c     lkgenmat               generate matrix coefficients
c     lkmatcorrect_actual    correct matrix coefficients for leakage condition (same as line-sinks)
c     lkkno_actual           generation of the known vector
c     lkupdate_actual        update conditions at collocation points
c     lksub_actual           substition of the leakage parameters
c     lksolve_actual         Gauss-Seidel solver for leakage
c     lkcombineequations_actual  combine equations for sub-cells into a single average equation for the MODFLOW cell
c     lksubgridsolve         subroutine is called in SOLUT to generate a new solution for all sub-cell leakage rates
c     lksubgridsolve_actual  called from entry lksubgridsolve_sub in LKIN
c
c
c ---------------------------------------------------------------------------------------------
c
c     Solution procedure:
c     Add all MODFLOW cells with a ilkresolution(i,j)>0 to the matrix.
c     Set as a boundary condition the potential at the center of the cell if ilkresolution(i,j)=1,
c     else use the average potential at all sub-grid cell centers (ilkresolution(i,j)>1)
c     This formulation will be designed in such a manner that the resulting leakage is equal to
c     or very close to the average leakage over the cell.
c     To calculate sub-grid leakages, solve for these in groups of nine cells
c     and store the difference with the average leakage for the MODFLOW cell.
c     Calculation of the field variables (heads, discharges, velocities) is done by evaluating all
c     analytic elements (including all MODFLOW cells) and then adding the sub-grid contributions of
c     the cell in which the point occurs and of the 8 surrounding cells.
c
c     Note: formulate the procedure above in terms of leakage increments from the existing leakage.
c
c -------------------------------------------------------------------------------------------------------
c
      subroutine lkmatsize_actual (M,N,nrow,ncol,ilkresolution)
c
c     Routine add the number of equations and number of unkowns to the counters M and N
c     Each cell may have a number of sub-cells defined by ilkresolution(irow,jcol).
c     If lkaveragehead=.true. overspecification is used by introducing 1 equation per sub-cell.
c     If lkaveragehead=.false. no overspecification is used; there will be only one equation per cell.
c
c     This routine is called by ENTRY lkmatsize(M,N) in LKIN in the file LKMOD.FOR
c
c    Note on the meaning of ilkresolution
c
c    ilkresolution(i,j)=0  cell is to be excluded from the matrix (initial leakage value to be kept)
c    ilkresolution(i,j)=1  cell has NO subgrid, only one collocation point at the center.
c    ilkresolution(i,j)=2  cell has a 2*2 subgrid, hence 4 sub-cells.
c    ilkresolution(i,j)=3  cell has a 3*3 subgrid, hence 9 sub-cells.
c            etc.
c
c    For ilkresolution>1 either "lkaveragehead" or "lkincludesubgrid" must be true!
c
c
      implicit none
      INTEGER M,N,nrow,ncol,ilkresolution,i,j,iresolution
      DIMENSION ilkresolution(nrow,ncol)
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
c
      if (lkleakage.and.lksolving) then
      do i=1,nrow
        do j=1,ncol
          iresolution=ilkresolution(i,j)
          select case (iresolution)
              case (0)
c              no equation, skip
              case (1)
               M=M+1
               N=N+1
              case (2:)
               if (lkaveragehead) then ! use average head condition, requires 1 equation per sub-cell
                M=M+iresolution*iresolution ! later, in lkmatcorrect_actual these equations will be combined into 1
                N=N+1
               else
                M=M+1
                N=N+1
                write (iluer,1000)
 1000 format (' WARNING: ilkresolution > 1, but lkaveragehead=false!')
               endif
              case default
               write (iluer,2000) i,j,iresolution
 2000 format (' ***ERROR in LKMATSIZE_ACTUAL:',/,
     & 'iresolution(',i4,','i4,')=',i6,' which is invalid.',/,
     & 'Program execution has been aborted.')
               AMESS(1)='Subgrid resolution out of range.'
               AMESS(2)='Stored in file: "basename.grs"'
               AMESS(3)='Correct and rerun.'
               CALL HALT(3) ! stop program execution for batch version
          end select
        end do
      end do
      end if
      return
c
      end subroutine
c
C---------------------------------------------------------------------------------------------------------
C
      SUBROUTINE lkczc_actual (CZI,M,N,mtot,DRFAC,CALPH,ITYPE,
     &             lDirectFromDisk,ltimer,
     &           nrow,ncol,rlkdeltax,rlkdeltay,
     &           ilkresolution,rlkleakage)
C
C---------------------------------------------------------------------------------------------------------
C
c     Routine generates the solution arrays: czi, itype, and drfac.
c
c
c
c
c
      IMPLICIT NONE
      INTEGER(4) M,N,ITYPE,i,j,nsolOut,nrow,ncol,
     &           irow,jcol,ii,jj,mtot,
     &           ilkresolution,iresolution
      LOGICAL lDirectFromDisk,
     &        lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut,ltimer
      REAL(8) DRFAC,rlkdeltax,rlkdeltay,rx,ry,rx1,ry1,rdx,rdy,
     &        resolution,rdxx,rdyy,rxx,ryy,rlkleakage
      COMPLEX(8) CZI,CALPH,CZ0,CZOLD,CZ0FIRST,
     &           cz1,cz2,cz3,cz
      CHARACTER(8) aBasenameOut
      CHARACTER(16)aDateTimeOut
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      DIMENSION CZI(mtot),DRFAC(4,mtot),CALPH(mtot),ITYPE(mtot)
      DIMENSION rlkdeltax(ncol),rlkdeltay(nrow),
     &          ilkresolution(nrow,ncol),rlkleakage(nrow,ncol)
c
      if (lkleakage.and.lksolving) then
c
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
c
      ilkmin=m+1 ! first equation for the leakage grid
      rx1=rlkx0
      ry1=rlky0
      do i=1,nrow
       rdy=rlkdeltay(i)
       ry=ry1-0.5*rdy
       do j=1,ncol
        iresolution=ilkresolution(i,j)
        rdx=rlkdeltax(j)
        rx=rx1+0.5*rdx
       select case (iresolution)
           case (0)
c           no equation, skip
c
           case (1)             ! add one equation to matrix
            m=m+1 ! one equation per cell
            n=n+1 ! one unknown leakage rate per cell
            czi(m)=CMPLX(rx,ry)
            calph(m)=(0.0,0.0)
            drfac(1,m)=1.0d0
            itype(m)=1
c            if (nsolOut.eq.0) then ! zero out the leakage in case first solve
c              rlkleakage(i,j)=0.0
c            end if
c
           case (2:)
           if (lkaveragehead) then
             n=n+1  ! add only one unkown, the leakage for the MODFLOW cell
             resolution=REAL(iresolution)
             call lk_cell_construct(i,j,rdx,rdy,cz0,cz1,cz2,cz3)
c
c     MODFLOW cell structure:
c                               Origin-----> columns
c     cz3-----------cz2         ||
c      |             |          ||
c      |             |          \/
c      |             |
c      |             |         rows
c     cz0-----------cz1
c
             rdxx=rdx/resolution
             rdyy=rdy/resolution
             rxx=0.5*rdxx
             ryy=0.5*rdyy
             cz=cz0+CMPLX(rxx,ryy)
c
c     example of subgrid control point generation:
c
c      -----------------------
c     |                       |
c     |   7       8       9   |
c     |                       |
c     |   4       5       6   |
c     |                       |
c     |   1       2       3   |
c     |                       |
c      -----------------------
c
             do ii=1,iresolution
             do jj=1,iresolution
               m=m+1  ! add an equation for every sub-cell
               czi(m)=cz
c      write (iluer,1001) ii,jj,m,rx,ry,rdxx,rdyy,cz
c 1001 format ('lkczc1: ii,jj,m,rx,ry,rdxx,rdyy,cz',/,
c     &         3i5,6(d14.7))
               calph(m)=(0.0,0.0)
               drfac(1,m)=1.0d0
               itype(m)=1
               cz=cz+CMPLX(rdxx,0.0d0)
             end do
               rxx=0.5*rdxx
               ryy=ryy+rdyy
               cz=cz0+CMPLX(rxx,ryy)
             end do
           else    ! resolution >1, but no average head requested, treat as case 1
             m=m+1 ! one equation per cell
             n=n+1 ! one unknown leakage rate per cell
             czi(m)=CMPLX(rx,ry)
             calph(m)=(0.0,0.0)
             drfac(1,m)=1.0d0
             itype(m)=1
c            if (nsolOut.eq.0) then ! zero out the leakage in case first solve
c              rlkleakage(i,j)=0.0
c            end if
           endif
c
           case default
            write (iluer,1000) i,j,iresolution
 1000 format (' ***ERROR in LKCZC_ACTUAL:',/,
     & 'iresolution(',i3,','i3,')=',i6,' which is invalid.',/,
     & 'Program execution has been aborted.')
            AMESS(1)='Subgrid resolution out of range.'
            AMESS(2)='Stored in file: "basename.grs"'
            AMESS(3)='Correct and rerun.'
            CALL HALT(3) ! stop program execution for batch version
       end select
c
        rx1=rx1+rdx
       end do
       rx1=rlkx0
       ry1=ry1-rdy
      end do
      end if
      ilkmax=m ! last equation for the leakage grid
      return
      end
C
C---------------------------------------------------------------------------------------------------------
C
      SUBROUTINE lkgenmat_actual (DRA,CZI,M,N,J,CALPH,ITYPE,
     &                  nrow,ncol,rlkconst,ilkresolution)
C
C---------------------------------------------------------------------------------------------------------
c
c     Generate initial matrix coefficients (without corrections)
c
c     ra(i,j) i=row number (equation number) j=element number (column number)
c
c     IMPORTANT
c     The idea is to let every cell generate a potential coefficient at its own collocation point
C     so the matrix equations can be used in "lkupdate_head" for fast recalculation of all potentials
C     at cell collocation points. The actual matrix coefficients necessary for the matrix solution
c     will be created by modifying these initial matrix coefficients in "lkmatcorrect."
c
      IMPLICIT NONE
      INTEGER(4) M,N,J,ITYPE,I,IEQS,IEQ,nrow,ncol,
     &           irow,jcol,ilkresolution,iresolution
      LOGICAL LNEG,linside
      REAL(8) DRA,rlkconst,rdx,rdy,rx,ry,rconst
      COMPLEX(8) CZI,CALPH,CZ1,CZ2,CZ3,CZ0,CZ,CZA,
     &           comega,cdum1,cdum2,comls,cflk_omega_par,cdum,
     &           cflk_omega_coefficient
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      DIMENSION DRA(M,N),CZI(M),CALPH(M),ITYPE(M)
      DIMENSION rlkconst(nrow,ncol+1),ilkresolution(nrow,ncol)
c
      if (lkleakage.and.lksolving) then
      do irow=1,nrow
       do jcol=1,ncol
       iresolution=ilkresolution(irow,jcol)
       select case (iresolution)
           case (0)
c           no equation, skip
           case (1:)               ! add a column to the matrix
       call lk_cell_construct (irow,jcol,rdx,rdy,cz0,cz1,cz2,cz3)
c       write (iluer,1001) irow,jcol,cz0,cz1,cz2,cz3,rdx,rdy
c 1001 format (' lkgenmat1: irow,jcol='2i5,/,
c     & 'cz0,cz1,cz2,cz3,rdx,rdy -->',/,4(2(d14.7),2x),d14.7,2x,d14.7)
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
       J=J+1    ! next column in coefficient matrix
       DO 10 I=1,M    ! add as many equations as there are control points generated in LKCZC_actual
       CZ=CZI(I)
c                      generate comega at CZ due to the cell.
c       cdum1=comls(cz,cz3,cz0)+comls(cz,cz2,cz1)+
c     &       rlkconst(irow,jcol)+rlkconst(irow,jcol+1)
c       cdum2=cflk_omega_par(cz,cz0,cz1)-cflk_omega_par(cz,cz3,cz2)
c       comega=0.5*rdx*cdum1+0.125*rdx*rdx*cdum2
       rconst=(rlkconst(irow,jcol)+rlkconst(irow,jcol+1))*0.5
       rx=REAL(cz-0.5*(cz0+cz2))
       ry=aimag(cz-0.5*(cz0+cz2))
       linside=ABS(rx).lt.0.5*rdx.and.ABS(ry).lt.0.5*rdx
       comega=cflk_omega_coefficient(cz,
     &                            cz0,rdx,rdy,rconst,linside)
c       if (ABS(rx).lt.0.5*rdx.and.ABS(ry).lt.0.5*rdx) THEN ! inside cell, add Poisson term
c        comega=comega+0.5*(rx-0.5*rdx)*(rx+0.5*rdx)
c      cdum=0.5*(rx-0.5*rdx)*(rx+0.5*rdx)
c      write (iluer,1002) rx,ry,cdum
c 1002 format (' lkgenmat2: rx,ry,cdum=',4(d14.7))
c       endif                                                    !!! REPLACE
c        write (iluer,1004) i,j,cz,cdum1,cdum2,comega
c 1004 format (' lkgenmat4: i,j,cz,cdum1,cdum2,comega-->'
c     &        ,i3,1x,i3,4(2(d14.7),2x))
c
       CZA=CALPH(I)
       IEQS=ITYPE(I)
C
C ITYPE=+1 potential specified at CZ
C       -1 difference in potential specified: PHI(CZ)-PHI(CZA)
C       +2 stream function specified at CZ
C       -2 flow across the line between CZ & CZA, positive from left to right when at CZ
C       +3 discharge component normal to the unit vector CZA (rotated to the left)
C       +4 discharge component parallel to the unit vector CZA
C        5 continuity equation: provide total discharge
C        6 request for zero matrix element
C
       LNEG=IEQS.LT.0
       IEQ=IABS(IEQS)
       GOTO (1,2,3,4,5,1),IEQ  ! note itype 6 is interpreted here as itype=1, to get desired potential
c                               coefficient at the collocation points.
   1   DRA(I,J)=DRA(I,J)+REAL(comega) ! provide potential at CZ
c                       note: handle difference in PHI in klmatcorrect
       GOTO 9
   2   IF (LNEG) THEN  !!! NOT YET PROVIDED, NECESSARY FOR HORIZONTAL BARRIERS!!!!
c       DRA(I,J)=DRA(I,J)+RFNFLSCO(IAD,CZ,CZA) ! provide flow across CZ & CZA (sink density 1)
       ELSE
        DRA(I,J)=DRA(I,J)+AIMAG(comega)  ! provide PSI at CZ
       END IF
       GOTO 9
   3   CONTINUE ! not provided
       GOTO 9
   4   CONTINUE ! not provided
       GOTO 9
   5   DRA(I,J)=DRA(I,J)+rdx*rdy    ! provide total discharge for continuity equation
   9   continue
   10  CONTINUE
           case default
            write (iluer,1000) i,j,iresolution
 1000 format (' ***ERROR in LKGENMAT_ACTUAL:',/,
     & 'iresolution(',i3,','i3,')=',i6,' which is invalid.',/,
     & 'Program execution has been aborted.')
            AMESS(1)='Subgrid resolution out of range.'
            AMESS(2)='Stored in file: "basename.grs"'
            AMESS(3)='Correct and rerun.'
            CALL HALT(3) ! stop program execution for batch version
       end select
       end do
      end do
      end if
      RETURN
      END
C
C---------------------------------------------------------------------------------------------------------
C
      SUBROUTINE lkmatcorrect_actual (DRA,CZI,M,N,J,CALPH,ITYPE,
     &              nrow,ncol,rlkconst,rlkresist,rlkpot,rlkheadlower,
     &              rlksubheadlower,ilkresolution,rlksubpotupper,
     &              rlksublambda2init,nbuf)
C
C---------------------------------------------------------------------------------------------------------
c
c     Apply corrections to matrix coefficients for Cauchy boundary condition
c
c     ra(i,j) i=row number (equation number) j=element number (column number)
c
      IMPLICIT NONE
      INTEGER(4) M,N,J,ITYPE,I,IEQS,IEQ,nrow,ncol,ires2,nbuf,
     &           irow,jcol,nsolOut,ilkresolution,iresolution,iad
      LOGICAL LNEG,lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           leakage_equation
      REAL(8) DRA,rlkconst,rlkresist,rdx,rdy,rh,rx,ry,
     &        rfperm,rfbase,rfhght,rfhead,RK0,RH0,RHED0,RB0,RP0,
     &        lambdasquared,rlkpot,rb,rfhedp,rlkheadlower,
     &        rlksubpotupper,rlksublambda2init,rlksubheadlower
      COMPLEX(8) CZI,CALPH,CZ1,CZ2,CZ3,CZ0,CZ,CZA,
     &           comega,cdum1,cdum2,comls,cflk_omega_par,
     &           czref
      CHARACTER(8) aBasenameOut
      CHARACTER(16)aDateTimeOut
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      DIMENSION DRA(M,N),CZI(M),CALPH(M),ITYPE(M)
      DIMENSION rlkconst(nrow,ncol+1),rlkresist(nrow,ncol),
     &          rlkpot(nrow,ncol),rlkheadlower(nrow,ncol),
     &          ilkresolution(nrow,ncol),rlksubpotupper(nbuf),
     &          rlksublambda2init(nbuf),rlksubheadlower(nbuf)
c
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      call GVPAR (RK0,RH0,RHED0,RB0,RP0,CZref)
c
      if (lkleakage.and.lksolving) then
      iad=0  ! address counter for 1D sub-cell arrays
c      write (iluer,1001) ilkmin,ilkmax
c 1001 format (' lkmatcorrect_actual1: ilkmin,ilkmax ',2i6)
      do irow=1,nrow
       do jcol=1,ncol
       iresolution=ilkresolution(irow,jcol)
c       write (iluer,1003) irow,jcol,iresolution
c 1003 format ('lkmatcorrect3: irow,jcol,iresolution ',3i5)
       select case (iresolution)
c
c
           case (0)
c           skip, no equation
c
c
           case (1)               ! correct equation                          ---- case (1) -----
       call lk_cell_construct (irow,jcol,rdx,rdy,cz0,cz1,cz2,cz3)
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
       J=J+1    ! next column in coefficient matrix
       DO I=1,M
       leakage_equation=i.ge.ilkmin.and.i.le.ilkmax ! check that we are in the range of the leakage equations
       CZ=CZI(I)
       CZA=CALPH(I)
       IEQS=ITYPE(I)
       if (ieqs.EQ.-1) then
c                      generate comega at CZA due to the cell.
       cdum1=comls(cza,cz3,cz0)+comls(cza,cz2,cz1)+
     &       rlkconst(irow,jcol)+rlkconst(irow,jcol+1)
       cdum2=cflk_omega_par(cza,cz3,cz2)-cflk_omega_par(cza,cz0,cz1)
       comega=0.5*rdx*cdum1+0.125*rdx*rdx*cdum2
       rx=REAL(cza-0.5*(cz0+cz2))
       if (ABS(rx).lt.0.5*rdx) THEN ! inside cell, add Poisson term
        rx=0.5*rdx-rx
        comega=comega+0.5*rx*(rdx-rx)
       end if
       END if
c
C
C ITYPE=+1 potential specified at CZ
C       -1 difference in potential specified: PHI(CZ)-PHI(CZA)
C       +2 stream function specified at CZ
C       -2 flow across the line between CZ & CZA, positive from left to right when at CZ
C       +3 discharge component normal to the unit vector CZA (rotated to the left)
C       +4 discharge component parallel to the unit vector CZA
C        5 continuity equation: provide total discharge
C        6 request for zero matrix element
C
       LNEG=IEQS.LT.0
       IEQ=ABS(IEQS)
       if (ieq.eq.1.and.leakage_equation) then   ! this assures that we are not hitting a coinciding control point of anoter element type
        if (ABS(cz-0.5*(cz0+cz2)).lt.0.01) then  ! at MODFLOW cell center, apply lambda squared term
         rb=rfbase(cz)
         if (nsolOut.eq.1) then ! nsol incremented after generation of czi and calls of update routine in SOLUT!
           rh=rhed0-rb  ! use head at reference point as the initial head (no solution yet)
         else  ! at least one solution completed
           rh=rfhedp(rlkpot(irow,jcol),cz)-rb
         endif
         if (rh.lt.rh0) then ! unconfined --> rh= (upperhead+lowerhead)/2
           rh=0.5*(rh+rlkheadlower(irow,jcol)-rb)
         else
           rh=rh0  ! confined --> rh= aquifer thickness
         end if
c      write (iluer,1001) nsolOut,rh
c 1001 format (' lkmatcorrect_actual1: nsolOut=',i5,' rh=',d14.7)
         lambdasquared=rfperm(cz)*rh*rlkresist(irow,jcol)
         DRA(I,J)=DRA(I,J)-lambdasquared
        end if
        IF (LNEG) THEN
          DRA(I,J)=DRA(I,J)-REAL(comega) ! subtract potential at CZA
        ENDIF
       endif
       if (ieq.eq.6) then ! zero matrix coefficient requested
        DRA(I,J)=0.0D0
       end if
       enddo
c
c
c
       case (2:)      ! apply correction to all subcells                      ---- case (2:) ----
c
       call lk_cell_construct (irow,jcol,rdx,rdy,cz0,cz1,cz2,cz3)
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
       J=J+1    ! next column in coefficient matrix (only one column for this MODFLOW cell, regardless of number of subcells)
      if (lkaveragehead) then
       DO I=1,M
       leakage_equation=i.ge.ilkmin.and.i.le.ilkmax ! check that we are in the range of the leakage equations
       CZ=CZI(I)
       CZA=CALPH(I)
       IEQS=ITYPE(I)
       if (ieqs.EQ.-1) then
c                      generate comega at CZA due to the cell.
       cdum1=comls(cza,cz3,cz0)+comls(cza,cz2,cz1)+
     &       rlkconst(irow,jcol)+rlkconst(irow,jcol+1)
       cdum2=cflk_omega_par(cza,cz3,cz2)-cflk_omega_par(cza,cz0,cz1)
       comega=0.5*rdx*cdum1+0.125*rdx*rdx*cdum2
c
       rx=REAL(cza-0.5*(cz0+cz2))
       if (ABS(rx).lt.0.5*rdx) THEN ! inside cell, add Poisson term
        rx=0.5*rdx-rx
        comega=comega+0.5*rx*(rdx-rx)
       end if
       END if
c
C
C ITYPE=+1 potential specified at CZ
C       -1 difference in potential specified: PHI(CZ)-PHI(CZA)
C       +2 stream function specified at CZ
C       -2 flow across the line between CZ & CZA, positive from left to right when at CZ
C       +3 discharge component normal to the unit vector CZA (rotated to the left)
C       +4 discharge component parallel to the unit vector CZA
C        5 continuity equation: provide total discharge
C        6 request for zero matrix element
C
       LNEG=IEQS.LT.0
       IEQ=ABS(IEQS)
       if (ieq.eq.1.and.leakage_equation) then  ! only when in range of leakage equations check for lambda squared correction
       rx=REAL(cz-0.5*(cz0+cz2))
       ry=AIMAG(cz-0.5*(cz0+cz2))
       rb=rfbase(cz)
       if (ABS(rx).lt.0.5*rdx.and.ABS(ry).lt.0.5*rdy) THEN ! at a sub-cell center, apply lambda squared term
c       write (iluer,1004)
c 1004  format ('lkmatcorrect_actual4: we came to here.')
           iad=iad+1             ! find address for sub-cell arrays
           if (iad.gt.nbuf) then
            WRITE(iluer,2000) iad,nbuf
           end if
         if (nsolOut.eq.1) then ! nsol incremented after generation of czi and calls of update routine in SOLUT!
           rh=rhed0-rb  ! use head at reference point as the initial head (no solution yet)
         else  ! at least one solution completed
           rh=rfhedp(rlksubpotupper(iad),cz)-rb
         endif
         if (rh.lt.rh0) then ! unconfined --> rh= (upperhead+lowerhead)/2
           rh=0.5*(rh+rlksubheadlower(iad)-rb)
         else
           rh=rh0  ! confined --> rh= aquifer thickness
         end if
c      write (iluer,1004) nsolOut,iad,rh
c 1004 format('lkmatcorrect_actual4: nsolOut=',i5,'iad=',i5,' rh=',d14.7)
         lambdasquared=rfperm(cz)*rh*rlkresist(irow,jcol)
         rlksublambda2init(iad)=lambdasquared ! keep for use in lkkno() !!!! **** Hmmm.. we will always have an inaccurate lambda
c                                             Perhaps do matrix correction and decomposition twice. Also for other non-linear
c                                             equations that would be affected by an inaccurate first head and potential.
c      write (iluer,1002) iad,lambdasquared
c 1002 format ('lkmatcorrect_actual2: rlksublambda2init(',i4,')=',d14.7)
         DRA(I,J)=DRA(I,J)-lambdasquared     ! note that this is an initial value, we will use updated values in lkkno()
        end if
        IF (LNEG) THEN
          DRA(I,J)=DRA(I,J)-REAL(comega) ! subtract potential at CZA
        ENDIF
       endif
       if (ieq.eq.6) then ! zero matrix coefficient requested
        DRA(I,J)=0.0D0
       end if
       enddo
c
      else ! sub-cells, but no average heads specified, treat as case 1
c
       DO I=1,M
       leakage_equation=i.ge.ilkmin.and.i.le.ilkmax ! check that we are in the range of the leakage equations
       CZ=CZI(I)
       CZA=CALPH(I)
       IEQS=ITYPE(I)
       if (ieqs.EQ.-1) then
c                      generate comega at CZA due to the cell.
       cdum1=comls(cza,cz3,cz0)+comls(cza,cz2,cz1)+
     &       rlkconst(irow,jcol)+rlkconst(irow,jcol+1)
       cdum2=cflk_omega_par(cza,cz3,cz2)-cflk_omega_par(cza,cz0,cz1)
       comega=0.5*rdx*cdum1+0.125*rdx*rdx*cdum2
       rx=REAL(cza-0.5*(cz0+cz2))
       if (ABS(rx).lt.0.5*rdx) THEN ! inside cell, add Poisson term
        rx=0.5*rdx-rx
        comega=comega+0.5*rx*(rdx-rx)
       end if
       END if
c
C
C ITYPE=+1 potential specified at CZ
C       -1 difference in potential specified: PHI(CZ)-PHI(CZA)
C       +2 stream function specified at CZ
C       -2 flow across the line between CZ & CZA, positive from left to right when at CZ
C       +3 discharge component normal to the unit vector CZA (rotated to the left)
C       +4 discharge component parallel to the unit vector CZA
C        5 continuity equation: provide total discharge
C        6 request for zero matrix element
C
       LNEG=IEQS.LT.0
       IEQ=ABS(IEQS)
       if (ieq.eq.1.and.leakage_equation) then   ! this assures that we are not hitting a coinciding control point of anoter element type
        if (ABS(cz-0.5*(cz0+cz2)).lt.0.01) then  ! at MODFLOW cell center, apply lambda squared term
         rb=rfbase(cz)
         if (nsolOut.eq.1) then ! nsol incremented after generation of czi and calls of update routine in SOLUT!
           rh=rhed0-rb  ! use head at reference point as the initial head (no solution yet)
         else  ! at least one solution completed
           rh=rfhedp(rlkpot(irow,jcol),cz)-rb
         endif
         if (rh.lt.rh0) then ! unconfined --> rh= (upperhead+lowerhead)/2
           rh=0.5*(rh+rlkheadlower(irow,jcol)-rb)
         else
           rh=rh0  ! confined --> rh= aquifer thickness
         end if
c      write (iluer,1001) nsolOut,rh
c 1001 format (' lkmatcorrect_actual1: nsolOut=',i5,' rh=',d14.7)
         lambdasquared=rfperm(cz)*rh*rlkresist(irow,jcol)
         DRA(I,J)=DRA(I,J)-lambdasquared
        end if
        IF (LNEG) THEN
          DRA(I,J)=DRA(I,J)-REAL(comega) ! subtract potential at CZA
        ENDIF
       endif
       if (ieq.eq.6) then ! zero matrix coefficient requested
        DRA(I,J)=0.0D0
       end if
       enddo
      endif

c
c         invalid resolution
c
           case default
            write (iluer,1000) irow,jcol,iresolution
            AMESS(1)='Subgrid resolution out of range.'
            AMESS(2)='Stored in file: "basename.grs"'
            AMESS(3)='Correct and rerun.'
            CALL HALT(3) ! stop program execution for batch version
c
c
       end select
       end do
      end do
      end if
      RETURN
 1000 format (' ***ERROR in LKMATCORRECT_ACTUAL:',/,
     & 'iresolution(',i3,','i3,')=',i6,' which is invalid.',/,
     & 'Program execution has been aborted.')
 2000 format (' ***ERROR in LKMATCORRECT_ACTUAL:',/,
     & 'subarray address',i7,' exceeds subarray buffer',i7)
      END
C
C---------------------------------------------------------------------------------------------------------
C
      SUBROUTINE lkkno_actual (DRB,J,CZI,n,
     &      nrow,ncol,rlkconst,rlkresist,rlkpot,rlkleakage,rlkheadlower,
     &     ilkresolution,rlksubpotupper,rlksubheadlower,
     &     rlksublambda2init,rlksubleakage,nbuf)
C
C---------------------------------------------------------------------------------------------------------
c
c     generate known vector for leakage elements
c
      IMPLICIT NONE
      INTEGER(4) J,nrow,ncol,iad,nbuf,n,
     &           irow,jcol,nsolOut,ii,jj,
     &           ilkresolution,iresolution
      LOGICAL LNEG,lDirectFromDiskOut,lErrorReportOut,
     &        lsolOut,loadsolOut,linalreadyOut,lwrite,lflkincludesubgrid
      REAL(8) DRB,rlkconst,rlkresist,rlkpot,rdx,rdy,rdxx,rdyy,rx,ry,
     &        rfperm,rfbase,RK0,RH0,RHED0,RB0,RP0,resolutionsquared,
     &        lambdasquared,rlkleakage,rlkheadlower,resolution,resist,
     &        rhead,rfhead,rk,rh,rb,rpotlower,rfpoth,rfhght,rfhedp,
     &        rheadlower,rlksubpotupper,rlksubheadlower,
     &        rlksublambda2init,rlksubleakage,rsubleak,rpotupper,
     &        lambda2ratio
      COMPLEX(8) CZI,CALPH,CZ1,CZ2,CZ3,CZ0,CZ,CZA,
     &           comega,cdum1,cdum2,comls,cflk_omega_par,
     &           czref
      CHARACTER(8) aBasenameOut
      CHARACTER(16)aDateTimeOut
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      DIMENSION DRB(n),CZI(n)
      DIMENSION rlkconst(nrow,ncol+1),rlkresist(nrow,ncol),
     &          rlkpot(nrow,ncol),rlkleakage(nrow,ncol),
     &          rlkheadlower(nrow,ncol),ilkresolution(nrow,ncol),
     &          rlksubpotupper(nbuf),rlksubheadlower(nbuf),
     &          rlksublambda2init(nbuf),rlksubleakage(nbuf)
c
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      call GVPAR (RK0,RH0,RHED0,RB0,RP0,CZref)
c
      if (lkleakage.and.lksolving) then
      iad=0
      do irow=1,nrow
       do jcol=1,ncol
                                                       lwrite=.false.
       iresolution=ilkresolution(irow,jcol)
       select case (iresolution)
c
       case (0)  ! no equation, skip
c
       case (1)  ! no sub-cells, only one unknown
c
       J=J+1
       CZ=CZI(j)
       rk=rfperm(cz)
       rb=rfbase(cz)
       rheadlower=rlkheadlower(irow,jcol)-rb   ! head: \phi
       if (nsolOut.eq.1) then ! nsol incremented after generation of czi and calls of update routine in SOLUT!
         rh=rhed0-rb  ! use head at reference point as the initial head (no solution yet)
       else
         rh=rfhedp(rlkpot(irow,jcol),cz)-rb  ! use calculated head
       endif
       if (rh.lt.rh0) then ! unconfined --> rh= (upperhead+lowerhead)/2
         rh=0.5*(rh+rheadlower)
         rpotlower=0.5*rk*rheadlower*rheadlower
       else
         rh=rh0  ! confined --> rh= aquifer thickness
         rpotlower=rk*rheadlower*rh0-0.5*rk*rh0*rh0
       end if
       lambdasquared=rk*rh*rlkresist(irow,jcol)
       drb(j)=drb(j)+rpotlower-rlkpot(irow,jcol)+
     &                lambdasquared*rlkleakage(irow,jcol)
c
       case (2:)              ! sub-cells present, calculate weighted average Known value (thus one value per MODFLOW cell)
c
       j=j+1
       if (lflkincludesubgrid().and.lkexclude) then
           drb(j)=0.0d0
       else
           rsubleak=0.0d0         ! S_i^k*
       if (lkaveragehead) then ! average heads over sub-grid centers of MODFLOW cell
           resolution=REAL(iresolution)
           resolutionsquared=resolution*resolution
           call lk_cell_construct(irow,jcol,rdx,rdy,cz0,cz1,cz2,cz3)
           rdxx=rdx/resolution
           rdyy=rdy/resolution
           rpotlower=0.0d0        ! Phi_i^l*
           rpotupper=0.0d0        ! Phi_i^u*
           lambdasquared=0.0d0    ! lambda_i*
           resist=rlkresist(irow,jcol)
           rx=0.5*rdxx
           ry=0.5*rdyy
           cz=cz0+CMPLX(rx,ry)
           do ii=1,iresolution
           do jj=1,iresolution
             rk=rfperm(cz)
             rb=rfbase(cz)
             iad=iad+1 ! next sub-cell address
             rheadlower=rlksubheadlower(iad)-rb   ! phi^lower
             if (nsolOut.eq.1) then
               rh=rhed0-rb  ! use head at reference point as the initial head (no solution yet)
5             else
               rh=rfhedp(rlksubpotupper(iad),cz)-rb ! use calculated head = phi^upper
             endif
             if (rh.lt.rh0) then ! unconfined --> rh= (upperhead+lowerhead)/2
               rh=0.5*(rh+rheadlower)
               rpotlower=rpotlower+
     &             0.5*rk*rheadlower*rheadlower/rlksublambda2init(iad)
             else
               rh=rh0  ! confined --> rh= aquifer thickness
               rpotlower=rpotlower+
     &       (rk*rheadlower*rh0-0.5*rk*rh0*rh0)/rlksublambda2init(iad)
             end if
             rpotupper=rpotupper+
     &         rlksubpotupper(iad)/rlksublambda2init(iad)
             lambda2ratio=rk*rh*resist/rlksublambda2init(iad) ! lambda2-current / lambda2-initial
             rsubleak=rsubleak+
     &        rlksubleakage(iad)*lambda2ratio
             lambdasquared=lambdasquared+lambda2ratio
     &
      if (lwrite)
     & write (iluer,1001) iad,rlksublambda2init(iad),lambda2ratio,
     &                    rk,rh,resist
 1001 format (' lkkno_actual1: iad,rlksublambda2init(iad),lambda2ratio,'
     & 'rk,rh,resist',/,i3,5(d14.7))
             cz=cz+CMPLX(rdxx,0.0d0)
           end do
             rx=0.5*rdxx
             ry=ry+rdyy
             cz=cz0+CMPLX(rx,ry)
           end do
           rpotupper=rpotupper/resolutionsquared         ! weighted average upper potential
           rpotlower=rpotlower/resolutionsquared         ! weighted average lower potential
           lambdasquared=lambdasquared/resolutionsquared ! weighted average lambda
           rsubleak=rsubleak/resolutionsquared           ! weighted average sub-cell leakage
      if (lwrite) write(iluer,1002) irow,jcol,iad,
     &rpotlower,rpotupper,rsubleak,lambdasquared
 1002 format (' lkkno_actual2: irow,jcol,iad,rpotlower,rpotupper,',
     &        'rsubleak,lambdasquared',/,3i3,4(d14.7))
       else
           CZ=CZI(j)                    ! treat as case (1) (repeat same logic)
           rk=rfperm(cz)
           rb=rfbase(cz)
           rheadlower=rlkheadlower(irow,jcol)-rb   ! head: \phi
           if (nsolOut.eq.1) then ! nsol incremented after generation of czi and calls of update routine in SOLUT!
             rh=rhed0-rb  ! use head at reference point as the initial head (no solution yet)
           else
             rh=rfhedp(rlkpot(irow,jcol),cz)-rb  ! use calculated head
           endif
             if (rh.lt.rh0) then ! unconfined --> rh= (upperhead+lowerhead)/2
               rh=0.5*(rh+rheadlower)
               rpotlower=0.5*rk*rheadlower*rheadlower
             else
               rh=rh0  ! confined --> rh= aquifer thickness
               rpotlower=rk*rheadlower*rh0-0.5*rk*rh0*rh0
            end if
            lambdasquared=rk*rh*rlkresist(irow,jcol)
            rpotupper=rlkpot(irow,jcol)
       end if
c
       drb(j)=drb(j)+rpotlower-rpotupper+
     &                lambdasquared*rlkleakage(irow,jcol)+rsubleak
c
       end if
c
       case default
        write (iluer,1000) irow,jcol,iresolution
 1000 format (' ***ERROR in LKKNO_ACTUAL:',/,
     & 'iresolution(',i3,','i3,')=',i6,' which is invalid.',/,
     & 'Program execution has been aborted.')
        AMESS(1)='Subgrid resolution out of range.'
        AMESS(2)='Stored in file: "basename.grs"'
        AMESS(3)='Correct and rerun.'
        CALL HALT(3) ! stop program execution for batch version
       END select
c
       end do
      end do
      end if
      RETURN
      END
C
C---------------------------------------------------------------------------------------------------------
C
      SUBROUTINE lkupdate_actual (DRSCR,j,m,n,czi,lsubinclude,isubsign,
     &      nrow,ncol,rlkpot,ilkresolution,rlksubpotupper,
     &      rlkheadlower,rlkheadupper,nbuf)
C
C---------------------------------------------------------------------------------------------------------
c
c     update the potential at collocation points of leakage elements
c
c     drscr(j)       increment of potential resulting from last matrix solution (j is number of matrix equation) (input)
c     m,n            number of equations and number of unlnowns, respectively (rows and columns of matrix) (input)
c     czi(j)         collocation point for j-th equation (input)
c     lsubinclude    ...
c     isubsign       ...
c     nrow,ncol      size of MODFLOW grid (input)
c     rlkpot         array of potentials at the centers of MODFLOW cells --> to be updated here (input & output)
c     ilkresolution  array of MODFLOW cell subgrid resolutions, see below (input)
c     rlksubpotupper array of uper potentials at sub-cell centers --> to be updated here (input & output)
c     rlkheadlower   array of lower heads at MODFLOW cell centers (input)
c     rlkheadupper   array of uper heads at MODFLOW cell centers (output)
c     nbuf           array size for rlksubpotupper
c
c    Note on the meaning of ilkresolution
c
c    ilkresolution(i,j)=0  cell is to be excluded from the matrix (initial leakage value to be kept)
c    ilkresolution(i,j)=1  cell has NO subgrid, only one collocation point at the center.
c    ilkresolution(i,j)=2  cell has a 2*2 subgrid, hence 4 sub-cells.
c    ilkresolution(i,j)=3  cell has a 3*3 subgrid, hence 9 sub-cells.
c            etc.
c
c   NOTE: rlkheadlower has been used in a test, it is not a part of this procedure.
c
      IMPLICIT NONE
      INTEGER(4) j,nrow,ncol,m,n,iad,nbuf,i1,isubsign,
     &           irow,jcol,nsolOut,ii,jj,iresolution,ilkresolution
      LOGICAL lDirectFromDiskOut,lErrorReportOut,
     &        lsolOut,loadsolOut,linalreadyOut,lsubinclude,
     &        laddsubcells
      REAL(8) DRSCR,rlkpot,rdx,rdy,rfpot,rdum,rlkheadupper,
     &        rpotaverage,resolution,rx,ry,rdxx,rdyy,rpot,rfhedp,
     &        rlksubpotupper,rhcenter,rlkheadlower,rfpoth,rpotother,
     &        rsubsign,RK0,RH0,RHED0,RB0,RP0
      COMPLEX(8) CZI,CZ,cz0,cz1,cz2,cz3,cflk_subomega,CZref
      CHARACTER(8) aBasenameOut
      CHARACTER(16)aDateTimeOut
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      DIMENSION DRSCR(m),CZI(m)
      DIMENSION rlkpot(nrow,ncol),ilkresolution(nrow,ncol),
     &          rlksubpotupper(nbuf),rlkheadlower(nrow,ncol),
     &          rlkheadupper(nrow,ncol)
c
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      call GVPAR (RK0,RH0,RHED0,RB0,RP0,CZref)
c
      if (lkleakage.and.lksolving) then
      laddsubcells=lsubinclude
      rsubsign=real(isubsign)
      iad=0
      do irow=1,nrow
       do jcol=1,ncol
        iresolution=ilkresolution(irow,jcol)
        resolution=REAL(iresolution)
c      write (iluer,990) iresolution,nsolOut
c 990  format (' lkupdate990: iresolution, nsolOut=',2i5)
        select case (iresolution)
        case (0)
c        no equation, skip
c
        case (1)
c
         J=J+1
         CZ=CZI(j)
         if (nsolOut.eq.0) then
          rlkpot(irow,jcol)=rfpot(cz)
c      rpotother=rfpot(cz)   !!!!! debugging only !!!!!!!!!!!!!!!!!!!!!!!!!
         else
          rlkpot(irow,jcol)=rlkpot(irow,jcol)+drscr(j) ! drscr calculated in "initialupdate" in GFMAT.FOR
          if (laddsubcells) then
            rlkpot(irow,jcol)=rlkpot(irow,jcol)+
     &                        rsubsign*REAL(cflk_subomega(cz))     ! sub-cell contributions not in DRSCR(j)
          end if
c      rpotother=rfpot(cz)   !!!!! debugging only !!!!!!!!!!!!!!!!!!!!!!!!!
         end if
c      write (iluer,991) irow,jcol,cz,rpotother,rlkpot(irow,jcol)
c 991  format ('lkupdate991: irow,jcol,cz,rpotother,rlkpot(irow,jcol) '
c     &        ,/,2(i5),2(d14.7),2x,2(d14.7))
         rlkheadupper(irow,jcol)=rfhedp(rlkpot(irow,jcol),cz) ! store upper head
c      rpot=rfpot(cz)   !!!!!!! debugging only !!!!!!!!!!!!!!!!!!!!!!!!!!
c      write (iluer,1001) nsolOut,irow,jcol,rpot,rlkpot(irow,jcol),
c     &                   cz,rlkheadupper(irow,jcol)
c 1001 format (' lkupdate_actual1: nsolOut,irow,jcol,rpot,'
c     &        ' rlkpot(irow,jcol)'
c     &,/,3i4,2(1x,d14.7),/,' cz,rlkheadupper ',2(d14.7),2x,d14.7)
c
c
        case (2:)
c
         if (lkaveragehead) then ! average heads over sub-grid centers of MODFLOW cell
           call lk_cell_construct(irow,jcol,rdx,rdy,cz0,cz1,cz2,cz3)
           rdxx=rdx/resolution
           rdyy=rdy/resolution
           rpotaverage=0.0d0
c       rpotother=0.0d0         !!!!! need this for debugging
           rx=0.5*rdxx
           ry=0.5*rdyy
           cz=cz0+CMPLX(rx,ry)
           do ii=1,iresolution
           do jj=1,iresolution
             J=J+1
             iad=iad+1 ! next sub-cell address
             if (nsolOut.eq.0) then
               rpot=rfpot(cz)
               rpotaverage=rpotaverage+rpot
c      rpotother=rpotother+rfpot(cz) !debugging, remove!!!!!!!!!!!!!!!
             else
               rpot=rlksubpotupper(iad)+drscr(j)
               if (laddsubcells) then
                  rpot=rpot+rsubsign*REAL(cflk_subomega(cz))     ! sub-cell contributions not in DRSCR(j)
c                                                             In terms of the analysis, here is where we add the Gi function.
               end if
c      write (iluer,992) cz,j,drscr(j)
c 992  format ('lkupdate992: cz,j,drscr ',2(d14.7),2x,i4,2x,d14.7)
               rpotaverage=rpotaverage+rpot
c      rpotother=rpotother+rfpot(cz) !debugging, remove!!!!!!!!!!!!!!!
             endif
             rlksubpotupper(iad)=rpot ! store new upper potentials at sub-cell centers (use this in lkkno() )
c      write(iluer,993)nsolOut,irow,jcol,ii,jj,cz,iad,rlksubpotupper(iad)
c 993  format ('lkupdate993: nsolOut,irow,jcol,ii,jj,cz,iad,'
c     &        ' rlksubpotupper(iad) '
c     &        ,/,5(i4),1x,2(d14.7),i4,1x,d14.7)
             cz=cz+CMPLX(rdxx,0.0d0)
           end do
             rx=0.5*rdxx
             ry=ry+rdyy
             cz=cz0+CMPLX(rx,ry)
           end do
c NOTE: rlkpot(irow,jcol) is not used when subgrid and average head are specified, but it is here the basis for generating
c       the average upper head to be handed to MODFLOW. This must be revisited!!!
           rlkpot(irow,jcol)=rpotaverage/resolution/resolution
           rlkheadupper(irow,jcol)=rfhedp(rlkpot(irow,jcol),cz) ! store average upper head
c      rpotother=rpotother/resolution/resolution !!! debugging only !!!
c
          else
           j=j+1                        ! sub-grid specified, but no average heads requested
           CZ=CZI(j)                    ! treat as case (1)
           if (nsolOut.eq.0) then
            rlkpot(irow,jcol)=rfpot(cz)
c      rpotother=rfpot(cz)   !!!!! debugging only !!!!!!!!!!!!!!!!!!!!!!!!!
           else
            rlkpot(irow,jcol)=rlkpot(irow,jcol)+drscr(j)
            if (laddsubcells) then
             rlkpot(irow,jcol)=rlkpot(irow,jcol)+
     &                         rsubsign*REAL(cflk_subomega(cz))     ! sub-cell contributions not in DRSCR(j)
            end if
c      rpotother=rfpot(cz)   !!!!! debugging only !!!!!!!!!!!!!!!!!!!!!!!!!
           end if
           rlkheadupper(irow,jcol)=rfhedp(rlkpot(irow,jcol),cz) ! store upper head
          endif
c      write(iluer,994)irow,jcol,cz,rpotother,rlkpot(irow,jcol)
c 994  format ('lkupdate994: irow,jcol,cz,rpotother,rlkpot(irow,jcol) '
c     &        ,/,2(i5),2(d14.7),2x,2(d14.7))
c
        case default
        write (iluer,1000) irow,jcol,iresolution
 1000 format (' ***ERROR in LKUPDATE_ACTUAL:',/,
     & 'iresolution(',i3,','i3,')=',i6,' which is invalid.',/,
     & 'Program execution has been aborted.')
        AMESS(1)='Subgrid resolution out of range.'
        AMESS(2)='Stored in file: "basename.grs"'
        AMESS(3)='Correct and rerun.'
        CALL HALT(3) ! stop program execution for batch version
        END select
c
c      debugging
c
c        rdum=rfpot(cz)
c        if (irow.eq.16.and.jcol.eq.16) THEN  ! pick cell 16,16
c      write (iluer,1001) nsolOut,irow,jcol,cz,rlkpot(irow,jcol),rdum
c 1001 format (' lkupdate_actual1: ',
c     &'nsol,irow,jcol,cz,rlkpot(irow,jcol),rfpot(cz) -->',/,
c     & 3i5,2x,2(d14.7),2x,2(d14.7,2x))
c       end if
c
        end do
      end do
      end if
      RETURN
      END
C
C---------------------------------------------------------------------------------------------------------
C
      SUBROUTINE lkupdate_check_actual (j,m,n,czi,
     &      nrow,ncol,rlkpot,ilkresolution,rlksubpotupper,nbuf)
C
C---------------------------------------------------------------------------------------------------------
c
c     routine is for debugging purposes only.
C     A comparison is made between rlkpot and rfpot
c
c     m,n            number of equations and number of unlnowns, respectively (rows and columns of matrix) (input)
c     czi(j)         collocation point for j-th equation (input)
c     nrow,ncol      size of MODFLOW grid (input)
c     rlkpot         array of potentials at the centers of MODFLOW cells --> to be updated here (input & output)
c     ilkresolution  array of MODFLOW cell subgrid resolutions, see below (input)
c     rlksubpotupper array of uper potentials at sub-cell centers --> to be updated here (input & output)
c     nbuf           array size for rlksubpotupper
c
c    Note on the meaning of ilkresolution
c
c    ilkresolution(i,j)=0  cell is to be excluded from the matrix (initial leakage value to be kept)
c    ilkresolution(i,j)=1  cell has NO subgrid, only one collocation point at the center.
c    ilkresolution(i,j)=2  cell has a 2*2 subgrid, hence 4 sub-cells.
c    ilkresolution(i,j)=3  cell has a 3*3 subgrid, hence 9 sub-cells.
c            etc.
c
c   NOTE: rlkheadlower has been used in a test, it is not a part of this procedure.
c
      IMPLICIT NONE
      INTEGER(4) j,nrow,ncol,m,n,iad,nbuf,i1,isubsign,
     &           irow,jcol,nsolOut,ii,jj,iresolution,ilkresolution
      LOGICAL lDirectFromDiskOut,lErrorReportOut,
     &        lsolOut,loadsolOut,linalreadyOut
      REAL(8) rlkpot,rdx,rdy,rfpot,rdum,
     &        resolution,rx,ry,rdxx,rdyy,
     &        rlksubpotupper,rpotother,
     &        RK0,RH0,RHED0,RB0,RP0
      COMPLEX(8) CZI,CZ,cz0,cz1,cz2,cz3,CZref
      CHARACTER(8) aBasenameOut
      CHARACTER(16)aDateTimeOut
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      DIMENSION CZI(m)
      DIMENSION rlkpot(nrow,ncol),ilkresolution(nrow,ncol),
     &          rlksubpotupper(nbuf)
c
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      call GVPAR (RK0,RH0,RHED0,RB0,RP0,CZref)
c
      if (lkleakage.and.lksolving) then
      iad=0
      do irow=1,nrow
       do jcol=1,ncol
        iresolution=ilkresolution(irow,jcol)
        resolution=REAL(iresolution)
        select case (iresolution)
        case (0)
c        no equation, skip
c
        case (1)
c
         J=J+1
         CZ=CZI(j)
          rpotother=rfpot(cz)
      write (ilume,991) j,rpotother,rlkpot(irow,jcol)
 991  format('lkupdate_check1: j,rpotother,rlkpot(irow,jcol)'
     &        ,/,i4,2x,2(d14.7))
c
        case (2:)
c
         if (lkaveragehead) then ! average heads over sub-grid centers of MODFLOW cell
           call lk_cell_construct(irow,jcol,rdx,rdy,cz0,cz1,cz2,cz3)
           rdxx=rdx/resolution
           rdyy=rdy/resolution
           rx=0.5*rdxx
           ry=0.5*rdyy
           cz=cz0+CMPLX(rx,ry)
           do ii=1,iresolution
           do jj=1,iresolution
             J=J+1
             iad=iad+1 ! next sub-cell address
             rpotother=rfpot(cz)
      write (ilume,992) j,rpotother,rlksubpotupper(iad)
 992  format ('lkupdate_check2: j,rpotother,rlksubpotupper(iad) ',/,
     & i4,2x,2(d14.7))
             cz=cz+CMPLX(rdxx,0.0d0)
           end do
             rx=0.5*rdxx
             ry=ry+rdyy
             cz=cz0+CMPLX(rx,ry)
           end do
          else
           j=j+1                        ! sub-grid specified, but no average heads requested
           CZ=CZI(j)                    ! treat as case (1)
           rpotother=rfpot(cz)
      write(ilume,993) j,rpotother,rlkpot(irow,jcol)
 993  format('lkupdate_check3: j,rpotother,rlkpot(irow,jcol)'
     &        ,/,i4,2x,2(d14.7))
          endif
c
        case default
        write (iluer,1000) irow,jcol,iresolution
 1000 format (' ***ERROR in LKUPDATE_ACTUAL:',/,
     & 'iresolution(',i3,','i3,')=',i6,' which is invalid.',/,
     & 'Program execution has been aborted.')
        AMESS(1)='Subgrid resolution out of range.'
        AMESS(2)='Stored in file: "basename.grs"'
        AMESS(3)='Correct and rerun.'
        CALL HALT(3) ! stop program execution for batch version
        END select
c
        end do
      end do
      end if
      RETURN
      END

C
C---------------------------------------------------------------------------------------------------------
c
      SUBROUTINE lksub_actual (DRB,j,nrow,ncol,rlkleakage,ilkresolution)
C
C---------------------------------------------------------------------------------------------------------
c
c     substitute leakages in rlkleakage array
c
c
      IMPLICIT NONE
      INTEGER(4) J,nrow,ncol,irow,jcol,
     &           ilkresolution, iresolution
      REAL(8) DRB,rdx,rdy,rlkleakage
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      DIMENSION DRB(*)
      DIMENSION rlkleakage(nrow,ncol),ilkresolution(nrow,ncol)
c
      if (lkleakage.and.lksolving) then
      do irow=1,nrow
       do jcol=1,ncol
       iresolution=ilkresolution(irow,jcol)
       select case (iresolution)
           case (0)
c           no equation, skip
           case (1)
            J=J+1
            rlkleakage(irow,jcol)=rlkleakage(irow,jcol)+drb(j)
           case (2:)
            j=j+1
            rlkleakage(irow,jcol)=rlkleakage(irow,jcol)+drb(j) ! like for case (1)
           case default
        write (iluer,1000) irow,jcol,iresolution
 1000 format (' ***ERROR in LKSUB_ACTUAL:',/,
     & 'iresolution(',i3,','i3,')=',i6,' which is invalid.',/,
     & 'Program execution has been aborted.')
        AMESS(1)='Subgrid resolution out of range.'
        AMESS(2)='Stored in file: "basename.grs"'
        AMESS(3)='Correct and rerun.'
        CALL HALT(3) ! stop program execution for batch version
       end select
       end do
      end do
      end if
      RETURN
      END
c
c ----------------------------------------------------------------------------------
c
      subroutine lksolve_actual (nrow,ncol,nbuf,               ! NOT IN USE !!
     &                 rlkdeltax,rlkdeltay,rlkheadlower,
     &                 rlkleakage,
     &                 rlkresist,rlksubheadlower,
     &                 rlksubleakage,ilkresolution)
c
c
c     THIS ROUTINE IS NOT UP-TO-DATE, e.g. does not use ilkresolution as a flag for leakage calculations
c     per cell.
c
c
c
c     Gauss-Seidel solution to leakages only (assuming all other analytic elements have known strengths)
c     Routine is called by ENTRY lksolve() in the recursive routine LKIN in file LKMOD.FOR
c
      implicit none
      INTEGER nrow,ncol,nbuf,ilkresolution,i,j,iter
      REAL(8)  rlkdeltax,rlkdeltay,rlkheadlower,
     &         rlkleakage,
     &         rlkresist,rlksubheadlower,
     &         rlksubleakage,rx,ry,rx1,ry1,rfhead,rds,rheadupper
      COMPLEX(8) cz
      DIMENSION rlkdeltax(ncol),rlkdeltay(nrow),rlkheadlower(nrow,ncol),
     &                 rlkleakage(nrow,ncol),
     &                 rlkresist(nrow,ncol),
     &                 ilkresolution(nrow,ncol),
     &                 rlksubheadlower(nbuf),rlksubleakage(nbuf)
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      write (iluer,1000)
 1000 format (' Warning: routine lksolve-actual is called,',
     &         '  but is not up-to-date!')
c
c      write (iluer,1001) lkleakage,nlkiter, rlkrelax
c 1001 format (//,' lksolve_actual1: entering lksolve ',/,
c     & ' lkleakage,nlkiter, rlkrelax ',l4,2x,i5,2x,d14.7)
      if (lkleakage.and.nlkiter.gt.0) then
       do iter=1,nlkiter
        rx1=rlkx0
        do j=1,ncol
         rx=rx1+0.5*rlkdeltax(j)
         ry1=rlky0
         do i=1,nrow
          ry=ry1-0.5*rlkdeltay(i)
c
          cz=CMPLX(rx,ry)
          rheadupper=rfhead(cz)
      rds=(rheadupper-rlkheadlower(i,j))/rlkresist(i,j)-rlkleakage(i,j)
          rlkleakage(i,j)=rlkleakage(i,j)+rlkrelax*rds              ! Gauss Seidel method if enabled
c
c      write (iluer,1002) i,j,cz,rheadupper,rlkheadlower(i,j),
c     &                   rlkresist(i,j),rds,rlkleakage(i,j)
c 1002 format (' lksolve_actual2:'
c     &        ' i,j,cz,rheadupper,rheadlower,resist,rds,rlkleakage '
c     &        ,/,2(i4),7(d14.7))
c
          ry1=ry1-rlkdeltay(i)
         end do
         rx1=rx1+rlkdeltax(j)
        end do
       end do
      end if
      end subroutine
c
c --------------------------------------------------------------------------------------------------------------------------
c
      subroutine lkcombineequations_actual (dra,czi,drfac,itype,lskip,
     &                                 m,mout,n,ilkresolution,nrow,ncol,
     &                                 rlksublambda2init,nbuf)
c
c     Routine combines the equations for all sub-cells in a MODFLOW cell
c     into one equation. This will make the matrix square.
c
c     The number of incoming equations is "m" and the number of outgoing equations
c     is "n" (unless other routines also produced more equations than unknowns).
c     The value of "m" is unaltered, but the new number of equations is given to "mout"
c
      implicit none
      LOGICAL lskip,lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           lflkincludesubgrid
      INTEGER m,n,nrow,ncol,irow,jcol,i,isteps,ist,nsolOut,
     &        ilkresolution,iresolution,itype,mloc,mout,iad,
     &        nbuf
      CHARACTER(8) aBasenameOut
      CHARACTER(16)aDateTimeOut
      REAL(8) dra,drfac,rlksublambda2init
      COMPLEX(8) czi
      DIMENSION dra(m,n),drfac(4,m),czi(m),itype(m),lskip(m),
     &          ilkresolution(nrow,ncol),rlksublambda2init(nbuf)
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
c
      mloc=m  ! do not change "m" !!
      if (lkleakage.and.lksolving.and.lkaveragehead) then
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
       i=ilkmin    ! first equation for the leakage grid
       iad=0       ! address of subgrid arrays
       do irow=1,nrow
        do jcol=1,ncol
         iresolution=ilkresolution(irow,jcol)
         select case (iresolution)
             case (0)
               ! skip, inactive cell
c
             case (1)
               i=i+1   ! move to next equation, no sub-cells
c
             case (2:)
             isteps=iresolution*iresolution  ! number of equations to be combined into one
             if (lflkincludesubgrid().and.lkexclude) then ! force leakage zero.  Note: ONLY AFTER .. ITERATIONS
              dra(i,1:n)=0.0d0
              dra(i,i)=1.0d01
             else                       ! use average of sub-cell equations
c

c
              do ist=1,isteps               ! divide equations by lambda^2_i
              iad=iad+1
              if (iad.gt.nbuf) then
                write (iluer,2000) nbuf
                AMESS(1)='Address for subgrid buffer out of range.'
                AMESS(2)='Program aborted. Contact Tech Support.'
                call halt(2)
              end if
              dra(i+ist-1,1:n)=dra(i+ist-1,1:n)/rlksublambda2init(iad)
              end do
              do ist=1,isteps-1             ! sum the modified equations
               dra(i,1:n)=dra(i,1:n)+dra(i+ist,1:n)
              end do
              dra(i,1:n)=dra(i,1:n)/isteps  ! generate weighted average equation
             end if
             i=i+1
             mloc=mloc-isteps+1
             do ist=i,mloc                  ! shift relevant arrays
              dra(ist,1:n)=dra(ist+isteps-1,1:n)
              czi(ist)=czi(ist+isteps-1)
              drfac(1:4,ist)=drfac(1:4,ist+isteps-1)
              itype(ist)=itype(ist+isteps-1)
              lskip(ist)=lskip(ist+isteps-1)
             end do
c
             case default
        write (iluer,1000) irow,jcol,iresolution
 1000 format (' ***ERROR in LKCOMBINEEQUATIONS_ACTUAL:',/,
     & 'iresolution(',i3,','i3,')=',i6,' which is invalid.',/,
     & 'Program execution has been aborted.')
        AMESS(1)='Subgrid resolution out of range.'
        AMESS(2)='Stored in file: "basename.grs"'
        AMESS(3)='Correct and rerun.'
        CALL HALT(3) ! stop program execution for batch version
         end select
        end do
       end do
      end if
      mout=mloc
c      write (iluer,1002) nrow,ncol,m,n,mloc
c 1002 format (' lcombineequations_actual2: nrow,ncol,m,n,mloc ',5i6)
c
 2000 format ('***Error in LKCOMBINEEQUATIONS_ACTUAL:',/,
     &' Address of array "rlksublambda2init" exceeds buffer size',i5,/,
     &' Program execution has been aborted. Contact Tech Support.')
      end subroutine
c
c ----------------------------------------------------------------------------------------------------------
c
      subroutine lksubgridsolve (nouter,dra,drb,drscr,n,m,imxsze,nstr)
c
c     Routine uses a Gauss-Seidel iteration scheme with successive over-relaxation
c     to solve for all sub-cell leakages
c     This routine is called in SOLUT in GFMAT.FOR
c
c     ninner =  ninner*ires = number of loops over the sub-cells in each MODFLOW cell
c     nouter =  number of loops over all MODFLOW cells (with sub-cells)
c     drb    =  solution vector with elements for modflow cells with subgrids set to zero in lksub
c     n      =  size of drb
c
      implicit none
      INTEGER nouter,nouter_local,n,m,imxsze,nstr
      LOGICAL lflkincludesubgrid
      REAL(8) dra,drb,drscr
      DIMENSION dra(imxsze,nstr),drb(imxsze),drscr(imxsze)
      INCLUDE 'lkcom.inc'
c
      if (lkleakage.and.lflkincludesubgrid()) then
      nouter_local=nouter    ! this way ninner and nouter cannot accidentally be reset
        call lksubgridsolve_sub (nouter_local,dra,drb,drscr,n,m,
     &                           imxsze,nstr) ! entry in LKIN in LKMOD.FOR
      end if
c
      end subroutine
c
c ------------------------------------------------------------------------------------------------------------
c
      subroutine lksubgridsolve_actual1 (ninner_base,nouter,drb,n,     ! CURRENTLY NOT IN USE
     &         nrow,ncol,nbuf,rlkleakage,rlkresist,ilkresolution,
     &         ilkpointer,rlksubpotupper,rlksubheadlower,rlksubleakage,
     &         rlkconst)
c
c     Routine is called by entry lksubgridsolve_sub in LKIN , which is called by
c     lksubgridsolve (above) to execute the actual solution procedure.
c     Note: this is the first version that directly applies Gauss-Seidel, without a matrix.
c
      implicit none
      integer i,j,irow,jcol,inner,iouter,istart,jstart,i0,j0,nrow,ncol,
     &        nbuf,ilkresolution,ilkpointer,ierr,itemp,ninner,nouter,
     &        ires,iadd,n,ninner_base
      LOGICAL lastloop,lflkincludesubgrid
      REAL(8) rx,ry,rxx,ryy,rdx,rdxx,rdy,rdyy,rtemp,rc,rs,rsMODFLOW,
     &        deltars,rlkleakage,rlkresist,rpot,rheadupper,rheadlower,
     &        rlksubpotupper,rlksubheadlower,rlksubleakage,rfhedp,relax,
     &        rs_sum,deltasMF,drb,rpot_sub,rtems,rlkconst
      COMPLEX(8) cz,cz0,cz1,cz2,cz3,cztemp,cflk_subomega
      DIMENSION rlkleakage(nrow,ncol),
     &          rlkresist(nrow,ncol),ilkresolution(nrow,ncol),
     &          ilkpointer(nrow,ncol),
     &          rlksubpotupper(nbuf),rlksubheadlower(nbuf),
     &          rlksubleakage(nbuf),
     &          drb(n)
      ALLOCATABLE rtems(:,:),rtemp(:,:),cztemp(:,:)
      parameter (itemp=25,relax=0.02)   !!!!!!!!!!!!! NOTE: select under- or over-relaxation.
      save
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
c
      if (lkleakage.and.lflkincludesubgrid()) then
       if (.not.ALLOCATED(rtemp)) then ! not yet allocated, do so
       ALLOCATE (rtems(itemp,itemp),rtemp(itemp,itemp),
     &           cztemp(itemp,itemp),stat=ierr)
       if (ierr.ne.0) then
        write (iluer,1000)
         AMESS(1)='Failed to allocate scratch arrays in '
         AMESS(2)='lksubgridsolve_actual. Verify resources and rerun.'
         CALL HALT(2) ! stop program execution for batch version
       end if
       end if
        do iouter=1,nouter
        iadd=ilkmin-1 ! ilkmin is first leakage in drb
        lastloop=iouter.eq.nouter
        do irow=1,nrow
         do jcol=1,ncol
         ires=ilkresolution(irow,jcol)
         if (ires.gt.0) iadd=iadd+1
         if (ires.gt.1) THEN ! MODFLOW cell contains sub-cells, solve for leakages
           ires=ilkresolution(irow,jcol)
           istart=ilkpointer(irow,jcol)
           rc=rlkresist(irow,jcol)
           if (rc.le.0.0d0) THEN ! quit, invalid resistance
             write (iluer,2000)
             AMESS(1)='Illegal MODFLOW cell resistance in '
             AMESS(2)='lksubgridsolve_actual. Program aborted.'
             CALL HALT(2) ! stop program execution for batch version
           end if
           rsMODFLOW=rlkleakage(irow,jcol)
           if (ires.gt.itemp) THEN ! increase temporary arrays
            deallocate (rtems,rtemp,cztemp)
            ALLOCATE (rtems(itemp,itemp),rtemp(itemp,itemp),
     &                cztemp(itemp,itemp),stat=ierr)
            if (ierr.ne.0) then
             write (iluer,1000)
             AMESS(1)='Failed to allocate scratch arrays in '
           AMESS(2)='lksubgridsolve_actual. Verify resources and rerun.'
             CALL HALT(2) ! stop program execution for batch version
            end if
           end if
c
c         fill temporary arrays
c
          call lk_cell_construct(irow,jcol,rdx,rdy,cz0,cz1,cz2,cz3)
c      write (iluer,9999) irow,jcol,rdx,rdy,cz0
c 9999 format (' lksubgridsolve9999: i,j,rdx,rdy,cz0 ',2i4,4(d14.7))
c
c     MODFLOW cell structure:
c                               Origin-----> columns
c     cz3-----------cz2         ||
c      |             |          ||
c      |             |          \/
c      |             |
c      |             |         rows
c     cz0-----------cz1
c
          rdxx=rdx/ires
          rdyy=rdy/ires
          rxx=0.5*rdxx
          ryy=0.5*rdyy
          cz=cz0+CMPLX(rxx,ryy)
c
c     example of sub-cell arrangement in a MODFLOW cell
c     numbers are the sequence of data in the subgrid buffers
c
c      -----------------------
c     |                       |
c     |   7       8       9   |
c     |                       |
c     |   4       5       6   |
c     |                       |
c     |   1       2       3   |
c     |                       |
c      -----------------------
c
c      write (iluer,1001) irow,jcol,ires,istart
c 1001 format (' lksubgridsolve_actual1: irow,jcol,ires,istart=',
c     &          4i5)
          do i=1,ires
           do j=1,ires
             cztemp(i,j)=cz
             jstart=istart+(i-1)*ires+j-1
             rtemp(i,j)=rlksubpotupper(jstart)
             rtems(i,j)=rlksubleakage(jstart)
c      write (iluer,1002) i,j,jstart,rdxx,rdyy,cz,cztemp(i,j),rtemp(i,j)
c 1002 format (' lksubgridsolve_actual2: i,j,jstart,rdxx,rdyy,cz,'
c     &        'cztemp(i,j),rtemp(i,j)=',/,
c     &                                  2i4,2x,i4,2x,7(d14.7))
             cz=cz+CMPLX(rdxx,0.0d0)
           end do
           rxx=0.5*rdxx
           ryy=ryy+rdyy
           cz=cz0+CMPLX(rxx,ryy)
          end do
c
c         end filling temporary arrays, start inner loop for Gauss-Seidel
c
      write (iluer,2001) iouter,irow,jcol,relax
 2001 format ('lksubgridsolve_actual<2001>: iouter=',i4,
     &' 9 subcell leakages for cell ',2i4,' and relax=',d14.7)
         ninner=ninner_base*ires
         do inner=1,ninner     ! inner Gauss-Seidel loop
          rs_sum=0.0d0
          do i=1,ires
           do j=1,ires
            jstart=istart+(i-1)*ires+j-1
            cz=cztemp(i,j)
            rpot_sub=REAL(cflk_subomega(cz))
            rpot=rtemp(i,j)+rpot_sub !  add current sub-cell contributions
c            rpot=rtemp(i,j) ! do not add current sub-cell contributions
            rheadupper=rfhedp(rpot,cz)
            rheadlower=rlksubheadlower(jstart)
            rs=rtems(i,j)
            deltars=(rheadupper-rheadlower)/rc-rs-rsMODFLOW
c            rlksubleakage(jstart)=rs+relax*deltars  ! update sub-cell leakage, to be used in next calculation of potential (rpot).
            rtems(i,j)=rs+relax*deltars ! rtems used to temporary store newly computed sub-cell leakage
      write (iluer,1003) i,j,jstart,rpot_sub,rheadupper,rheadlower,rs,
     &                   rsMODFLOW,deltars,rtems(i,j)
 1003 format(' lksubgridsolve_actual3: i,j,jstart,rpot_sub,rheadupper,'
     &        'rheadlower,rs,rsMODFLOW,deltars,rtems(i,j)',/,
     &        3i4,7(d14.7))
             if (lastloop.and.inner.eq.ninner) then
               rs_sum=rs_sum+rlksubleakage(jstart)
             end if
           end do
          end do
          do i=1,ires
           do j=1,ires
            jstart=istart+(i-1)*ires+j-1
            rlksubleakage(jstart)=rtems(i,j)
           end do
          end do
      write (iluer,2002) inner,istart,rlksubleakage(istart:istart+8)
 2002 format ('lksubgridsolve_actual2002:',2i5,2x,9(d14.7))
         end do
c
c       end of inner loop for Gauss-Seidel
c
c       if last iteration, update MODFLOW cell leakage and store difference in DRB
c
         if (lastloop) then
           deltasMF=rs_sum/ires/ires
           drb(iadd)=deltasMF   ! used in initialupdate followed by --update calls to update BC values
           rlkleakage(irow,jcol)=rsMODFLOW+deltasMF
      write (iluer,1004) irow,jcol,iadd,rsMODFLOW,deltasMF,
     &                   rlkleakage(irow,jcol)
 1004 format (' lksubgridsolve_actual4: irow,jcol,iadd,rsMODFLOW,'
     &        ' deltasMF, rlkleakage ',
     &          3i4,1x,3(d14.7))
           do i=1,ires
            do j=1,ires
             jstart=istart+(i-1)*ires+j-1
             rlksubleakage(jstart)=rlksubleakage(jstart)-deltasMF
      write (iluer,1005) i,j,jstart,rlksubleakage(jstart)
 1005 format (' lksubgridsolve_actual5: i,j,jstart,rlksubleakage',
     &        3i4,1x,d14.7)
            end do
           end do
         end if
c
c        end of last iteration updates
c
         end if
         end do
        end do
       end do
      end if
      return
 1000 format (' ***ERROR in lksubgridsolution_actual:',/,
     &    'Failed to allocate arrays rtemp and cztemp, program halted.')
 2000 format (' ***ERROR in lksubgridsolution_actual:',/,
     &    'Illegal resistance value in lksubgridsolution_actual, '
     &    'program halted.')
      end subroutine
c
c ------------------------------------------------------------------------------------------------------------
c
      subroutine lksubgridsolve_actual2 (nouterin,
     &         dra,drb,drscr,n,m,imxsze,nstr,
     &         nrow,ncol,nbuf,rlkleakage,rlkresist,ilkresolution,
     &         ilkpointer,rlksubpotupper,rlksubheadlower,rlksubleakage,
     &         rlkconst)
c
c     Routine is called by entry lksubgridsolve_sub in LKIN , which is called by
c     lksubgridsolve (above) to execute the actual solution procedure.
c     Note: this is the second version that uses a coefficient matrix.
c
      implicit none
      integer i,j,irow,jcol,inner,iouter,istart,jstart,i0,j0,nrow,ncol,
     &        nbuf,ilkresolution,ilkpointer,ierr,itemp,nouter,
     &        nouterin,ires,iadd,n,m,imat,ieq,jeq,ii,jj,neq,
     &        ipiv,nscrsize,ist,ien,nres,imxsze,nstr
      LOGICAL lflkincludesubgrid,lmat,linside
      REAL(8) rx,ry,rxx,ryy,rdx,rdxx,rdy,rdyy,rtemp,rc,rs,rsMODFLOW,
     &        deltars,rlkleakage,rlkresist,rpot,rheadupper,rheadlower,
     &        rlksubpotupper,rlksubheadlower,rlksubleakage,rfhedp,relax,
     &        rs_sum,deltasMF,dra,drb,drscr,rpot_sub,rtems,rmat,rkno,
     &        rsol,rfbase,rlkconst,rconst,rdum1,rk,rh,rlambdasquared,
     &        rfperm,rfhght,rcond,scr1,scr2,rdum2,rb,rlab2,rh0,rftop,
     &        rhlower,rpotlower,rhupper,rscal,rsmat,rfhead,rdum
      COMPLEX(8) cz,cz0,cz1,cz2,cz3,cztemp,cflk_subomega,cdum1,cdum2,
     &           comega,comls,cflk_omega_par,cflk_omega_coefficient
      DIMENSION rlkleakage(nrow,ncol),
     &          rlkresist(nrow,ncol),ilkresolution(nrow,ncol),
     &          ilkpointer(nrow,ncol),
     &          rlksubpotupper(nbuf),rlksubheadlower(nbuf),
     &          rlksubleakage(nbuf),rlkconst(nrow,ncol+1),
     &          dra(imxsze,nstr),drb(imxsze),drscr(imxsze)
      DATA itemp /3/
      ALLOCATABLE rtems(:,:),rtemp(:,:),cztemp(:,:),ipiv(:)
      ALLOCATABLE rmat(:,:),rkno(:),rlab2(:,:),scr1(:),scr2(:)
      parameter (relax=0.02)   !!!!!!!!!!!!! NOTE: select under- or over-relaxation.   Not used at present!
      save
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
c
      nouter=nouterin
      if (lkleakage.and.lflkincludesubgrid()) then
       if (.not.ALLOCATED(rtemp)) then ! arrays not yet allocated, do so
       imat=itemp*itemp
       nscrsize=max(imat,1+(imat*imat)/31)
       ALLOCATE (rtems(itemp,itemp),rtemp(itemp,itemp),
     &           cztemp(itemp,itemp),rlab2(itemp,itemp),rmat(imat,imat),
     &           rkno(imat),scr1(nscrsize),scr2(imat),ipiv(imat),
     &           stat=ierr)
c
c     rtems(i,j)      is the strength parameter of the sub-cell i,j
c     rtemp(i,j)      is the upper potential of the sub-cell i,j
c     cztemp(i,j)     is the collocaton point (center) of the sub-cell i,j
c     rlab2(i,j)      is the lambda squared value for sub-cell i,j
c     rmat(imat,imat) are the matrix coefficients (initial, corrected and decomposed) for the subcells.
c     rkno(imat)      is the known vector
c     scr1(nscrsize)  is a scratch array for decomp routine
c     scr2(imat)      is a scratch array for decomp routine
c     ipiv(imat)      is the pivit vector generated by decomp and used by solve routine
c
c     just to test
       if (ierr.ne.0) then
        write (iluer,1000)
         AMESS(1)='Failed to allocate scratch arrays in '
         AMESS(2)='lksubgridsolve_actual. Verify resources and rerun.'
         CALL HALT(2) ! stop program execution for batch version
       end if
       end if
        do iouter=1,nouter
        iadd=ilkmin-1 ! ilkmin is first leakage in drb
c        write (iluer,1001) rlksubleakage(1:9),rlksubleakage(10:18)
c 1001 format (' rlksubgridsolve_actual1: rlksubleakage(1:9) ',9(d14.7),/
c     &        '                         rlksubleakage(10:18) ',9(d14.7))
        do irow=1,nrow    ! loop over all MODFLOW cells
         do jcol=1,ncol
         ires=ilkresolution(irow,jcol)
         if (ires.gt.0) iadd=iadd+1
         if (ires.gt.1) THEN ! MODFLOW cell contains sub-cells, solve for leakages
           nres=ires*ires
           istart=ilkpointer(irow,jcol)
           rc=rlkresist(irow,jcol)
           if (rc.le.0.0d0) THEN ! quit, invalid resistance
             write (iluer,2000)
             AMESS(1)='Illegal MODFLOW cell resistance in '
             AMESS(2)='lksubgridsolve_actual. Program aborted.'
             CALL HALT(2) ! stop program execution for batch version
           end if
           rsMODFLOW=rlkleakage(irow,jcol)
           if (ires.ne.itemp) THEN ! make arrays fit the number of equations.
            deallocate (rmat)
            deallocate (rtems)
            deallocate (rtemp)
            deallocate (cztemp)
            deallocate (rkno)
            deallocate (rlab2)
            deallocate (scr1)
            deallocate (scr2)
            deallocate (ipiv)
            itemp=ires
            imat=itemp*itemp
            nscrsize=max(imat,1+(imat*imat)/31)
       ALLOCATE (rtems(itemp,itemp),rtemp(itemp,itemp),
     &           cztemp(itemp,itemp),rlab2(itemp,itemp),rmat(imat,imat),
     &           rkno(imat),scr1(nscrsize),scr2(imat),ipiv(imat),
     &           stat=ierr)
            if (ierr.ne.0) then
             write (iluer,1000)
             AMESS(1)='Failed to allocate scratch arrays in '
           AMESS(2)='lksubgridsolve_actual. Verify resources and rerun.'
             CALL HALT(2) ! stop program execution for batch version
            end if
           end if
c
c         fill temporary arrays
c
          call lk_cell_construct(irow,jcol,rdx,rdy,cz0,cz1,cz2,cz3)
c      write (iluer,9999) irow,jcol,rdx,rdy,cz0
c 9999 format (' lksubgridsolve9999: i,j,rdx,rdy,cz0 ',2i4,4(d14.7))
c
c     MODFLOW cell structure:
c                               Origin-----> columns
c     cz3-----------cz2         ||
c      |             |          ||
c      |             |          \/
c      |             |
c      |             |         rows
c     cz0-----------cz1
c
          cz=0.5*(cz0+cz2) ! *****Strange, seems inconsequential
          rdxx=rdx/ires
          rdyy=rdy/ires
          rxx=0.5*rdxx
          ryy=0.5*rdyy
          cz=cz0+CMPLX(rxx,ryy)
c
c     example of sub-cell arrangement in a MODFLOW cell
c     numbers are the sequence of data in the subgrid buffers
c
c      -----------------------
c     |                       |
c     |   7       8       9   |
c     |                       |
c     |   4       5       6   |
c     |                       |
c     |   1       2       3   |
c     |                       |
c      -----------------------
c
          do i=1,ires
           do j=1,ires
             cztemp(i,j)=cz
             jstart=istart+(i-1)*ires+j-1
             rtemp(i,j)=rlksubpotupper(jstart)+REAL(cflk_subomega(cz)) ! rlksubpotupper() does not yet contain
c                                                           the sub-cell contributions, which are constantly being updated.
             rtems(i,j)=rlksubleakage(jstart)
             rh=rfhedp(rtemp(i,j),cz)
             if (rh.lt.rftop(cz)) THEN ! unconfined
              rh=0.5*(rh+rlksubheadlower(jstart))-rfbase(cz)
             else                ! confined
              rh=rftop(cz)-rfbase(cz)
             end if
             rk=rfperm(cz)
             rlab2(i,j)=rk*rh*rc
c      write (iluer,1002) i,j,jstart,rdxx,rdyy,cz,cztemp(i,j),rtemp(i,j),
c     &                   rh,rk,rc,rlab2(i,j)
c 1002 format (' lksubgridsolve_actual2: i,j,jstart,rdxx,rdyy,cz,'
c     &        'cztemp(i,j),rtemp(i,j)=',/,
c     &                                  2i4,2x,i4,2x,7(d14.7),/,
c     &        'rh,rk,rc,rlab2(i,j)',4(d14.7))
             cz=cz+CMPLX(rdxx,0.0d0)
           end do
           rxx=0.5*rdxx
           ryy=ryy+rdyy
           cz=cz0+CMPLX(rxx,ryy)
          end do
c
c     construct initial coeficient matrix (Aij*Sj=PHIi)  ---- use of DIRECT SOLUTION (Decomp and Solve)
c
          rconst=rlkconst(irow,jcol)/ires ! same as in cflk_subomega
          ieq=0
          do i=1,ires
           do j=1,ires
            ieq=ieq+1
            cz=cztemp(i,j) ! center of first sub-cell
            jeq=0
            do ii=1,ires
             do jj=1,ires
             jeq=jeq+1
             cz0=cztemp(ii,jj)-CMPLX(0.5*rdxx,0.5*rdyy) ! lower left corner of subcell
             linside=ieq.eq.jeq
             comega=cflk_omega_coefficient(cz,
     &                            cz0,rdxx,rdyy,rconst,linside)
             rmat(ieq,jeq)=REAL(comega)
             end do
            end do
           end do
          end do
c
c        correct the matrix diagonal coefficient
c
         ieq=0
         do i=1,ires
         do j=1,ires
         ieq=ieq+1
         rmat(ieq,ieq)=rmat(ieq,ieq)-rlab2(i,j)
         end do
         end do
c
c         generate known vector  (iterative refinement)
c
          ieq=0
          do ii=1,ires
           do jj=1,ires
            ieq=ieq+1
            cz=cztemp(ii,jj)
            jstart=istart+(ii-1)*ires+jj-1
            rk=rfperm(cz)
            rh0=rftop(cz)-rfbase(cz)
            rdum1=0.5*rk*rh0*rh0
            rhlower=rlksubheadlower(jstart)-rfbase(cz)
            rpot=rtemp(ii,jj) ! upper potential
            if (rpot.lt.rdum1) THEN ! upper is unconfined,
              rpotlower=0.5*rk*rhlower*rhlower
            else
              rpotlower=rk*rh0*rhlower-rdum1
            end if
            rs=rsMODFLOW+rtems(ii,jj)
            rkno(ieq)=rpotlower-rpot+rs*rlab2(ii,jj)
           end do
          end do
c
c         solve using Decomp and Solve          NOTE: debugging hardwired to ires=3, imat=9
c
          call decomp (imat,rmat,rcond,ipiv,scr1,scr2)
c          write (iluer,5111) rcond
c 5111 format ('lksubgridsolve_actual5111: rcond',d14.7)
          call solve (imat,rmat,rkno,ipiv)
c
c         check last sub-cell condition prior to solving
c          cz=cztemp(ires,ires)
c          rhupper=rfhedp(rtemp(ires,ires),cz)-rfbase(cz)
c          WRITE(iluer,3454) rhupper,rhlower
c 3454 format ('lksubgridsolve_actual2 3454 no update: ',
c     &    'rhupper,rhlower ',2(d14.7))
c          rc=rlkresist(irow,jcol)
c          rscal=(rhupper-rhlower)/rc
c          rsmat=rtems(ires,ires)+rsMODFLOW
c          write (iluer,3455) irow,jcol,ires,rscal,rsmat     ! are only equal if the MODFLOW leakages are accurate, but they are probably not....
c 3455 format('lksubgridsolve_actual2 3455 no update: ',
c     &      'irow,jcol,ires,rscal,rsmat:',
c     &       3i5,(2(D14.7)))
c
c
c         Now update all boundary conditions due to the change in strenght of this MODFLOW cell
c
          drscr(1:m)=0.0d0 ! initialize the array with strength supplements
          rs_sum=SUM(rkno(1:nres))
          deltasMF=rs_sum/nres
          rkno(1:nres)=rkno(1:nres)-deltasMF ! this will make sum of sub-cell strenghts equal to zero.
          drscr(1:m)=dra(1:m,iadd)*deltasMF    ! in lieu of an initialupdate() call
c          write (iluer,3456) irow,jcol,iadd
c 3456 format ('lksubgridsolve_actual2 3456: irow,jcol,iadd',3I8)
c          write (iluer,3457)drscr(1:m)
c 3457 format ('lksubgridsolve_actual2 3457: drscr',10(d14.7))
c
          call update (drscr)  !------------------ entry in SOLUT in file GFMAT
c
          ieq=0      ! store solution in sub-cell strength arrays
          do ii=1,ires
           do jj=1,ires
            ieq=ieq+1
            jstart=istart+(ii-1)*ires+jj-1
            rlksubleakage(jstart)=rlksubleakage(jstart)+rkno(ieq) ! iterative refinement
           end do
          end do
          rlkleakage(irow,jcol)=rsMODFLOW+deltasMF ! update MODFLOW cell strength array
c
c         end direct solution procedure and B.C. update for sub-cells in MODFLOW cell (irow,jcol)
c
c         check last sub-cell solution
c          cz=cztemp(ires,ires)
c          rhupper=rfhead(cz)-rfbase(cz)
c          rdum=rlksubpotupper(jstart)+REAL(cflk_subomega(cz))
c          rdum=rfhedp(rdum,cz)-rfbase(cz)
c          write (iluer,3459) rhupper,rdum   ! if equal, the update process seems to work!
c 3459 format ('lksubgridsolve_actual2 3459 after update: ',
c     &      'rhupper,rdum:',2(D14.7))
c          rc=rlkresist(irow,jcol)
c          rscal=(rhupper-rhlower)/rc
c          rsmat=rlksubleakage(jstart)+rlkleakage(irow,jcol)
c          write (iluer,3458) irow,jcol,ires,rscal,rsmat     ! are only equal if the MODFLOW leakages are accurate, but they probably are not....
c 3458 format('lksubgridsolve_actual2 3458 after update: ',
c     &      'irow,jcol,ires,rscal,rsmat:',
c     &       3i5,(2(D14.7)))
c

         end if
         end do
        end do
c
c        end of loop over all MODFLOW cells
c
       end do
c
c        end of iterations for sub-cell strength calculations
c
      end if
      return
 1000 format (' ***ERROR in lksubgridsolution_actual:',/,
     &    'Failed to allocate arrays rtemp and cztemp, program halted.')
 2000 format (' ***ERROR in lksubgridsolution_actual:',/,
     &    'Illegal resistance value in lksubgridsolution_actual, '
     &    'program halted.')
      end subroutine

