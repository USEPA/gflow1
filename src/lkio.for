C     Last change:  HMH  13 Nov 2017    5:40 pm
c
c lkextract         routine writes sink density and target leakage to basename.xtr file.
c lkio              routine writes leakage common blocks and allocatable arrays to basename.sol file
c lkio_array_actual routine writes the allocatable arrays (called in LKIN by entry LKIO_array)
c lk_write_filename_vlk routine writes the current sink densities (leakage) to basename.vlk
c lk_write_filename_rs  routine writes the current resistances in the MODFLOW grid to the file basename.rst
c lk_write_filename_uph routine writes upper heads at cell centers to basename.uph
c lk_scan_filename_del  scan the file basename.del to get the grid dimensions to allocate arrays
c lk_read_filename_del  now read the file to get the x and y spacings
c lk_scan_filename_grs  scan the file basename.grs to get the total number of subcells to allocate arrays
c lk_read_filename_grs  now read and store the grid resolutions per cell.
c lk_read_filename_fmt  read the heads in the lower aquifer (upper MODFLOW aquifer) at the centers of the MODFLOW cells
c lk_read_filename_rla  read the resistances for the lower aquifer (upper MODFLOW aquifer) for the MODFLOW cells
c lk_read_filename_rst  read the total resistances (for the lower aquifer (upper MODFLOW aquifer) and the GFLOW aquifer)
c lk_read_filename_vlk  read the leakages into the lower aquifer (upper MODFLOW aquifer) for the MODFLOW cells
c lkreadnewleakages_actual     reads only those leakages that could have been changed by MODFLOW
c lk_read_filename_rta  read recharge array to add recharges from MODFLOW
c lk_read_filename_top  read aquifer top elevations and use to replace -900 flags in lower heads.
c lk_write_filename_ler write errors at collocation points to file
c lkgfmfmonitor_actual  write (append) data to file during successive GFLOW calls in GF-MF conj. sol. process


c
c -----------------------------------------------------------------------------------
c
      subroutine lkextract (ilu)
c
c     Currentky only used to write total leakage from a sub-set of MODFLOW cells
c     and the total leakage calculated from upper, lower head and resistance.
c
      implicit none
      INTEGER ilu
      INCLUDE 'lkcom.inc'
      if (lkcheck) then
       write (ilu,1000)
       write (ilu,2000) rQ_leakage,rQ_target_leakage
      end if
      return
 1000 format ('* leakage   sum cells     integrated calculated leakage')
 2000 format (' leakage ',2(d14.7,2x))
      end subroutine
C
c -----------------------------------------------------------------------------------
c
      subroutine lkio (ICODE,ILU,RVERSION,ierr)
c
c     Routine is called by BIO in GFIO.FOR
c
C       Routine reads or writes contents of COMMON /LKCOM/ to
C       an external file.
C       ICODE=41 write
C       ICODE=73 read
C
C
      IMPLICIT NONE
      INTEGER ICODE,ILU,IERR,
     &        nrow,ncol,nbuf
      LOGICAL LDUM      
      REAL(8) RVERSION
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      save
C
      if (ierr.ne.0) return
      IF (RVERSION.le.12.0) RETURN ! execute only if in version 12.0 or higher
      call bufio4 (nlkrowsize,1,ilu,icode,ierr)
      if (nlkrowsize.eq.0) RETURN ! no leakage grid available
      call bufio4 (nlkcolsize,1,ilu,icode,ierr)
      call bufio4 (nlkbufsize,1,ilu,icode,ierr)
      call bufio4 (nlkiter,1,ilu,icode,ierr)
      call bufiol (lkleakage,1,ilu,icode,ierr)
      call bufiol (lkincludesubgrid,1,ilu,icode,ierr)
      call bufiol (lkcheck,1,ilu,icode,ierr)
      call bufiol (lksolving,1,ilu,icode,ierr)
      call bufiol (lkaveragehead,1,ilu,icode,ierr)
      call bufiol (lkupdateresistances,1,ilu,icode,ierr)
      call bufiod (rlkx0,1,ilu,icode,ierr)
      call bufiod (rlky0,1,ilu,icode,ierr)
      call bufiod (rlkrelax,1,ilu,icode,ierr)
      if (icode.eq.73) then
        nrow=nlkrowsize
        ncol=nlkcolsize
        nbuf=nlkbufsize
        call lkallocate(nrow,ncol,nbuf) ! make detour through LKIN
      endif
      call lkio_array (icode,ilu,rversion)
      RETURN
      end
c
c ------------------------------------------------------------------------------------
c
      subroutine lkio_array_actual (ICODE,ILU,RVERSION,nrow,ncol,nbuf,
     &           rlkdeltax,rlkdeltay,rlkheadlower,rlkleakage,rlkresist,
     &           rlkconst,rlksubheadlower,rlksubleakage,ilkresolution,
     &           rlkpot,rlkrecharge,rlksubpotupper,ilkpointer,
     &           rlksublambda2init)
c
c     Routine is called by lkio_array in LKIN (file LKMOD.FOR)
c
C       Routine reads or writes leakage arrays to
C       an external file, the basename.sol file.
C       ICODE=41 write
C       ICODE=73 read
C
C
      IMPLICIT NONE
      INTEGER ICODE,ILU,IERR,i,
     &        nrow,ncol,nbuf,ilkresolution,ntot,ilkpointer
      LOGICAL LDUM      
      REAL RVERSION,
     &     rlkdeltax,rlkdeltay,rlkheadlower,rlkleakage,rlkresist,
     &     rlkconst,rlksubheadlower,rlksubleakage,rlkpot,rlkrecharge,
     &     rlksubpotupper,rlksublambda2init
      DIMENSION rlkdeltax(ncol),rlkdeltay(nrow),rlkheadlower(nrow,ncol),
     &          rlkleakage(nrow,ncol),rlkresist(nrow,ncol),
     &          rlkconst(nrow,ncol+1),rlkpot(nrow,ncol),
     &          rlksubheadlower(nbuf),rlksubleakage(nbuf),
     &          ilkresolution(nrow,ncol),rlkrecharge(nrow,ncol),
     &          rlksubpotupper(nbuf),ilkpointer(nrow,ncol),
     &          rlksublambda2init(nbuf)
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      save
C
      ntot=nrow*ncol
      call bufiod (rlkheadlower,ntot,ilu,icode,ierr)
      call bufiod (rlkleakage,ntot,ilu,icode,ierr)
      call bufiod (rlkresist,ntot,ilu,icode,ierr)
      call bufio4 (ilkresolution,ntot,ilu,icode,ierr)
      call bufiod (rlkpot,ntot,ilu,icode,ierr)
      call bufiod (rlkdeltax,ncol,ilu,icode,ierr)
      call bufiod (rlkdeltay,nrow,ilu,icode,ierr)
      call bufiod (rlksubheadlower,nbuf,ilu,icode,ierr)
      call bufiod (rlksubleakage,nbuf,ilu,icode,ierr)
      ntot=nrow*(ncol+1)
      call bufiod(rlkconst,ntot,ilu,icode,ierr)
      if (rversion.ge.14) then
      ntot=nrow*ncol
       call bufiod (rlkrecharge,ntot,ilu,icode,ierr)
      end if
      if (rversion.ge.15) then
       call bufiod (rlksubpotupper,nbuf,ilu,icode,ierr)
      END if
      if (rversion.ge.16) then
       call bufio4 (ilkpointer,ntot,ilu,icode,ierr)
      end if
      if (rversion.ge.17) then
       call bufiod (rlksublambda2init,nbuf,ilu,icode,ierr)
      end if
c
c
c                  Start debugging
c
c
      RETURN ! bypass debugging
c
c
      write (iluer,1122) icode,nrow,ncol,nlkrowsize,nlkcolsize,nbuf,ntot
     &                   ,rlkx0,rlky0
      write (iluer,1123) (i,rlkdeltay(i),i=1,nlkrowsize)
      write (iluer,1124) (i,rlkdeltax(i),i=1,nlkcolsize)
 1122 format (//,' lkio: icode=',i5,/,
     &' nrow,ncol,nlkrowsize,nlkcolsize,nbuf, ntot ',6(i5,2x),/,
     &' rlkx0,rlky0 ',2(d14.7),/)
 1123 format (' rlkdeltay(',i4,')',d14.7)
 1124 format (' rlkdeltax(',i4,')',d14.7)
      write (iluer,1125) nlkrowsize,nlkcolsize
 1125 format (//,' ilkresolution: nlkrowsize,nlkcolsize ',2(i5,2x))
      do i=1,nlkrowsize
      write (iluer,*) ilkresolution(i,1:nlkcolsize)
      end do
      write (iluer,1126) nlkrowsize,nlkcolsize
 1126 format (//,' rlkheadlower: nlkrowsize,nlkcolsize ',2(i5,2x))
      do i=1,nlkrowsize
      write (iluer,*) rlkheadlower(i,1:nlkcolsize)
      end do
      write (iluer,1127) nlkrowsize,nlkcolsize
 1127 format (//,' rlkresist: nlkrowsize,nlkcolsize ',2(i5,2x))
      do i=1,nlkrowsize
      write (iluer,*) rlkresist(i,1:nlkcolsize)
      end do
      write (iluer,1128) nlkrowsize,nlkcolsize
 1128 format (//,' rlkleakage: nlkrowsize,nlkcolsize ',2(i5,2x))
      do i=1,nlkrowsize
      write (iluer,*) rlkleakage(i,1:nlkcolsize)
      end do
c
c                 End debugging
c
      RETURN
      end
c
c ---------------------------------------------------------------------------------
c
      subroutine lk_write_filename_vlk (nrow,ncol,rlkleakage,
     &                                  ilkresolution)
c
c     This routine writes the current leakage into the lower aquifer
c     Routine is called by the entry lkwrite_leakage in LKIN in file LKMOD.FOR
c
      implicit none
      INTEGER nsolOut,i,j,ilu,itemp,ivar,ierr,ilen,icharacter,
     &        nrow,ncol,ilkresolution
      LOGICAL lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut
      REAL(8) rdum,rvar,rlkleakage,scratch,rincrement,rfac,rfac1,
     &        raverageleakage
      CHARACTER*8 aBasenameOut
      CHARACTER*16 aDateTimeOut
      DIMENSION rlkleakage(nrow,ncol),ilkresolution(nrow,ncol)
      ALLOCATABLE scratch(:,:)
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
c
       if (.not.lksolving) return
c
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      ilu=2
      afile=TRIM(aBasenameOut)//".vlk"
      OPEN (UNIT=ILU,FILE=AFILE,STATUS='OLD',IOSTAT=IERR)
      if (ierr.ne.0) then  ! failed to open file
        write (iluer,1000) afile
        return
      endif
      if (rlkrelax.le.1.0) THEN ! apply relaxation factor to leakages before writing.
        ALLOCATE (scratch(nrow,ncol))
        do i=1,nrow
          read (ilu,*) scratch(i,1:ncol) ! read old leakages again
        end do
        do i=1,nrow
          do j=1,ncol
            if (ilkresolution(i,j).gt.0) then  ! keep leakage zero if resolution=0
             rincrement=rlkleakage(i,j)-scratch(i,j) ! calculate leakage increment
c             if (i.eq.19.and.j.eq.17) then
c              write (iluer,1001) rlkleakage(i,j),scratch(i,j),rincrement
c 1001 format ('lk_write_filename_vlk1: rlkleakage(i,j),scratch(i,j),',
c     &'rincrement,afile',3(d14.7),a16)
c             end if
c
c Option 1: The statements below are for a straight forward relaxation , but keeping 0 leakages 0 (see condition above)
c
                scratch(i,j)=scratch(i,j)+rlkrelax*rincrement
c                if (i.eq.19.and.j.eq.17) then
c                  write (iluer,1002) scratch(i,j),rlkrelax,rincrement
c 1002 format ('lk_write_filename_vlk2: scratch(i,j),rlkrelax,',
c     &'rincrement ',3(d14.7))
c                end if
c
c Option 2: The statements below are for a leakage change limitation based on the average leakage
c
c             rfac=MIN(ABS(rincrement),rlkrelax*raverageleakage)
c             rfac=SIGN(rfac,rincrement)
c             scratch(i,j)=scratch(i,j)+rfac
c             if (ABS(rfac).eq.rlkrelax*raverageleakage) then
c               write (iluer,2000) i,j,rincrement,rfac
c             end if
c
c Option 3: The statements below are for a leakage limitation as a percentage of the old leakage on a cell by cell basis
c
c             rfac1=rincrement/scratch(i,j)           ! calculate leakage increment as fraction of old leakage
c             rfac=MIN(ABS(rfac1),rlkrelax)     ! limit change in leakage to rlkrelax  (% change is rlkrelax*100)
c             rfac=SIGN(rfac,rfac1)                  ! restore sign
c             scratch(i,j)=(1.0+rfac)*scratch(i,j) ! this implies a variable relaxation factor of rlkrelax/rfac1
c for the case that abs(rfac1) > rlkrelax, else the relaxation factor is 1.
c             if (ABS(rfac).eq.rlkrelax) then ! do this only for options 2 & 3
c               write (iluer,2000) i,j,rfac1,rfac
c             end if
            endif
          end do
        end do
        REWIND ilu
c
c     now write the relaxed leakage values row by row, but keep the actual leakages in rlkleakage
c
        do i=1,nrow
          write (ilu,*) scratch(i,1:ncol) ! list directed write
        end do
c
        deallocate (scratch)
c
      else
c
c     write the actual leakage values row by row
c
c
        do i=1,nrow
          write (ilu,*) rlkleakage(i,1:ncol) ! list directed write
        end do
      endif
      close (ilu)
      return
c
 1000 format (' ERROR in lk_write_filename_vlk:',/,
     &        ' The file ',a16,' is not found.',/,
     &        ' No leakage values are written.')
 2000 format (' WARNING: leakage change in cell(',i3,',',i3,') ',
     &'has been limited from ',d14.7,' to ',d14.7)
c 3000 format (200(d14.7))
      end subroutine
c
c ---------------------------------------------------------------------------------
c
      subroutine lk_write_filename_rst (nrow,ncol,rlkresist)
c
c     This routine writes the current resistances to a file
c     Routine is called by the entry lkwrite_resistance in LKIN in file LKMOD.FOR
c
      implicit none
      INTEGER nsolOut,i,j,ilu,itemp,ivar,ierr,ilen,icharacter,
     &        nrow,ncol
      LOGICAL lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut
      REAL(8) rdum,rvar,rlkresist
      CHARACTER*8 aBasenameOut
      CHARACTER*16 aDateTimeOut
      DIMENSION rlkresist(nrow,ncol)
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
c
       if (.not.lksolving) return     ! do not write resistances unless we are solving
c
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      ilu=2
      afile=TRIM(aBasenameOut)//".rst"
      OPEN (UNIT=ILU,FILE=AFILE,STATUS='UNKNOWN',IOSTAT=IERR)
      if (ierr.ne.0) then  ! failed to open file
        write (iluer,1000) afile
        return
      endif
c
c     write the leakage values row by row
c
c
      do i=1,nlkrowsize
       write (ilu,*) rlkresist(i,1:nlkcolsize) ! list directed write
      end do
c
      close (ilu)
      return
c
 1000 format (' ERROR in lk_write_filename_resist:',/,
     &        ' The file ',a16,' is not found.',/,
     &        ' No leakage values are written.')
      end subroutine
c

c
c ---------------------------------------------------------------------------------
c
      subroutine lk_write_filename_uph (nrow,ncol,
     &                   rlkdeltax,rlkdeltay,rlkheadupper,ilkresolution)
c
c     This routine writes the current upper heads to disk whenever a leakage grid is present
c     Routine is called by the entry lkwrite_upperheads in LKIN in file LKMOD.FOR
c
c     When a leakage grid is present and the leakage is being solved for there are two options:
c     1) The upper heads in
c     the rlkheadupper array are set in LKCHECK in such a manner that the current lower head,
c     upperhead, and newly calculated leakage are consistent. These upper heads are read by MODFLOW for recalculating
c     the lower heads (this is done instead of passing the new leakage rates themselves).
c     2) The upper heads are calculated as the average of mres*mres heads in each MODFLOW cell (except those on the grid boundary)
c     regardless of the subgrid resolution. These average upper heads are then passed to MODFLOW.
c     This option is activated in this routine by simply overwriting the values calculated in LKCHECK, see option 1) above.
c
c     When the leakage grid is present, but the leakage is not being solved for, the upper heads in
c     the rlkheadupper array are calculated in this routine as follows
c     res=0 or res=1  calculate head at MODFLOW cell centers and store in rlkheadupper(i,j)
c     res>1           calculate heads at sub-cell centers and store average in rlkheadupper(i,j)
c
      implicit none
      INTEGER nsolOut,i,j,ii,jj,ilu,ierr,nrow,ncol,iresolution,
     &        ilkresolution,mres
      LOGICAL lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut,loption2
      REAL(8) rdum,rvar,rlkdeltax,rlkdeltay,resolution,
     &        rx,ry,rdx,rdy,rdxx,rdyy,scratch,rfhead,rlkheadupper
      complex(8) cz,cz0,cz1,cz2,cz3
      CHARACTER*8 aBasenameOut
      CHARACTER*16 aDateTimeOut
      DIMENSION rlkdeltax(ncol),rlkdeltay(nrow),rlkheadupper(nrow,ncol)
      DIMENSION ilkresolution(nrow,ncol)
      ALLOCATABLE scratch(:)
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'

c                                       loption2=.TRUE.  ! implement option 2 (default is option 1), see comments above.
                                       loption2=.false.  ! implement option 1 (default is option 1), see comments above.
c
      if (lkleakage.and.lksolving) then ! write upperheads only when we solve for leakages  (used to be always - changed 4/16/2015)
        call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
        ilu=2
        afile=TRIM(aBasenameOut)//".uph"
        OPEN (UNIT=ILU,FILE=AFILE,STATUS='UNKNOWN',IOSTAT=IERR)
        if (ierr.ne.0) then  ! failed to open file
          write (iluer,1000) afile
          return
        endif
c
c     Fill rlkheadupper array for the case that no leakages are being calculated
c
      if (.not.lksolving) then
c
       RETURN ! this is new! We are not writing upper heads unless we solve for leakages (change on 4/16/2015)
c
        do i=1,nrow
        do j=1,ncol
         call lk_cell_construct(i,j,rdx,rdy,cz0,cz1,cz2,cz3)
         iresolution=ilkresolution(i,j)
         resolution=REAL(iresolution)
         select case (iresolution)
c
             case (:1) ! resolution is 0 or 1
c
              cz=0.5*(cz0+cz2)  ! MODFLOW cell center
              rlkheadupper(i,j)=rfhead(cz)
c
             case (2:) ! resolution is 2 or larger
c
               rdxx=rdx/resolution
               rdyy=rdy/resolution
               rx=0.5*rdxx
               ry=0.5*rdyy
               cz=cz0+CMPLX(rx,ry)
               rdum=0.0d0
               do ii=1,iresolution
               do jj=1,iresolution
                rdum=rdum+rfhead(cz)
                cz=cz+CMPLX(rdxx,0.0d0)
               end do
               rx=0.5*rdxx
               ry=ry+rdyy
               cz=cz0+CMPLX(rx,ry)
               end do
                rlkheadupper(i,j)=rdum/resolution/resolution
c
             case default
c
             write (iluer,2000) iresolution,i,j    ! incorrect resolution number
c
         end select
c
       end do
       end do
c
      elseif (loption2) then ! we are solving, now overrule the calculation of upperheads in LKCHECK (option 2 in description of this routine).
        mres=0
        do i=1,nrow     ! find largest subgrid resolution
        do j=1,ncol
         iresolution=ilkresolution(i,j)
         if (iresolution.gt.mres) then
          mres=iresolution
         end if
        end do
        end do
c
       if (mres.gt.1) then ! we have subgrids
        do i=2,nrow-1    ! update all cells, except the boundary cells.
        do j=2,ncol-1
         call lk_cell_construct(i,j,rdx,rdy,cz0,cz1,cz2,cz3)
               rdxx=rdx/mres  ! use larges subgrid resolution
               rdyy=rdy/mres
               rx=0.5*rdxx
               ry=0.5*rdyy
               cz=cz0+CMPLX(rx,ry)
               rdum=0.0d0
               do ii=1,mres
               do jj=1,mres
                rdum=rdum+rfhead(cz)
                cz=cz+CMPLX(rdxx,0.0d0)
               end do
               rx=0.5*rdxx
               ry=ry+rdyy
               cz=cz0+CMPLX(rx,ry)
               end do
                rlkheadupper(i,j)=rdum/mres/mres
        end do
        end do
       else   ! no subgrids in any of the MODFLOW cells, do not overwrite rlkheadupper
             write (iluer,3000)
       endif
      end if
c
c     write the upper heads values row by row
c
c
       do i=1,nrow
         write (ilu,4000) rlkheadupper (i,1:ncol)
       end do
       close (ilu)
      endif
      return
c
 1000 format (' ERROR in lk_write_filename_uph:',/,
     &        ' The file ',a16,' could not be opened.',/,
     &        ' No upper heads are written.')
 2000 format (' ERROR in lk_write_filename_uph:',/,
     &        ' invalid resolution ',i5,'in cell (',i5,';',i5')')
 3000 format (' Warning: No subgrids found in lk_write_filename_uph.',/,
     &        ' Upperheads are not averaged in MODFLOW cells.')
 4000 format (200(d14.7))
      end subroutine
c
c --------------------------------------------------------------------
c
      subroutine lk_scan_filename_del(nlkrowsize,nlkcolsize)
c
c     Routine scans the MODFLOW grid dimensions to allocate the arrays for
c     storing data for the MODFLOW grid.
c
      implicit none
      INTEGER nsolOut,i,ilu,itemp,ivar,ierr,ilen,icharacter,
     &        nlkcolsize,nlkrowsize
      LOGICAL lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut
      REAL(8) rdum,rvar
      CHARACTER*8 aBasenameOut
      CHARACTER*16 aDateTimeOut
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
c
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      ilu=2
      afile=TRIM(aBasenameOut)//".del"
      itemp=iluin ! temporarily store iluin and replace for this read
      iluin=ilu
      OPEN (UNIT=ILU,FILE=AFILE,STATUS='OLD',IOSTAT=IERR)
      if (ierr.ne.0) then  ! failed to open file
        write (iluer,1000) afile
        AMESS(1)='The file basename.del is not found.'
        AMESS(2)='Correct and rerun.'
        CALL HALT(2) ! stop program execution for batch version
      endif
      call inline
      nlkrowsize=ivar(1) ! number of rows
      if (lmiss.or.lerror) then
        write (iluer,2000)
        AMESS(1)='Incorrect syntax in basename.del'
        AMESS(2)='Correct and rerun.'
        CALL HALT(2) ! stop program execution for batch version
      endif
      do i=1,nlkrowsize
       call inline
       rdum=rvar(1) ! dummy read to check syntax and forward to # rows
       if (lmiss.or.lerror) then
        write (iluer,2000)
        AMESS(1)='Incorrect syntax in basename.del'
        AMESS(2)='Correct and rerun.'
        CALL HALT(2) ! stop program execution for batch version
       endif
      end do
c
      call inline
      nlkcolsize=ivar(1) ! number of columns
      if (lmiss.or.lerror) then
        write (iluer,2000)
        AMESS(1)='Incorrect syntax in basename.del'
        AMESS(2)='Correct and rerun.'
        CALL HALT(2) ! stop program execution for batch version
      endif
      do i=1,nlkcolsize
       call inline
       rdum=rvar(1) ! dummy read to check syntax
       if (lmiss.or.lerror) then
        write (iluer,2000)
        AMESS(1)='Incorrect syntax in basename.del'
        AMESS(2)='Correct and rerun.'
        CALL HALT(2) ! stop program execution for batch version
       endif
      end do
c
      close (iluin)
      iluin=itemp ! restore iluin
      return
c
 1000 format (' ERROR in lk_scan_filename_del:',/,
     &        ' The file ',a16,' is not found.',/,
     &        ' program execution is aborted.')
 2000 format (' ERROR in lk_scan_filename_del:',/,
     &        ' Incorrect syntax in file ',a16,/,
     &        ' program execution is aborted.')
      end subroutine
c
c --------------------------------------------------------------------
c
      subroutine lk_read_filename_del(nrow,ncol,rlkdeltax,rlkdeltay)
c
c     Routine reads the MODFLOW grid dimensions.
c
      implicit none
      INTEGER nsolOut,i,ilu,itemp,ivar,ierr,ilen,icharacter,
     &        nrow,ncol
      LOGICAL lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut
      CHARACTER*8 aBasenameOut
      CHARACTER*16 aDateTimeOut
      REAL(8) rlkdeltax,rlkdeltay,rvar
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
      Dimension rlkdeltax(ncol),rlkdeltay(nrow)

      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      ilu=2
      afile=TRIM(aBasenameOut)//".del"
      itemp=iluin ! temporarily store iluin and replace for this read
      iluin=ilu   ! necessary because we are using the routine "inline"
      OPEN (UNIT=ILU,FILE=AFILE,STATUS='OLD',IOSTAT=IERR)
      if (ierr.ne.0) then  ! failed to open file
        write (iluer,1000) afile
        AMESS(1)='The file basename.del is mot found.'
        AMESS(2)='Correct and rerun.'
        CALL HALT(2) ! stop program execution for batch version
      endif
      call inline
      nlkrowsize=ivar(1) ! number of rows
      if (lmiss.or.lerror) then
        write (iluer,2000)
        AMESS(1)='Incorrect syntax in basename.del'
        AMESS(2)='Correct and rerun.'
        CALL HALT(2) ! stop program execution for batch version
      endif
      do i=1,nlkrowsize
       call inline
       rlkdeltay(i)=rvar(1)
       if (lmiss.or.lerror) then
        write (iluer,2000)
        AMESS(1)='Incorrect syntax in basename.del'
        AMESS(2)='Correct and rerun.'
        CALL HALT(2) ! stop program execution for batch version
       endif
      end do
c
      call inline
      nlkcolsize=ivar(1) ! number of columns
      if (lmiss.or.lerror) then
        write (iluer,2000)
        AMESS(1)='Incorrect syntax in basename.del'
        AMESS(2)='Correct and rerun.'
        CALL HALT(2) ! stop program execution for batch version
      endif
      do i=1,nlkcolsize
       call inline
       rlkdeltax(i)=rvar(1)
       if (lmiss.or.lerror) then
        write (iluer,2000)
        AMESS(1)='Incorrect syntax in basename.del'
        AMESS(2)='Correct and rerun.'
        CALL HALT(2) ! stop program execution for batch version
       endif
      end do
c
      close (iluin)
      iluin=itemp ! restore iluin
      return
c
 1000 format (' ERROR in lk_read_filename_del:',/,
     &        ' The file ',a16,' is not found.',/,
     &        ' program execution is aborted.')
 2000 format (' ERROR in lk_read_filename_del:',/,
     &        ' Incorrect syntax in file ',a16,/,
     &        ' program execution is aborted.')
c
      end subroutine
c
c ---------------------------------------------------------------------------------
c
      subroutine  lk_scan_filename_grs(nlkrowsize,nlkcolsize,nlkbufsize,
     &                                 lnofile)
c
c     Routine reads the sub-grid resolution for the MODFLOW cells
c     to calculate the required array sizes for storing sub-grid data.
c
c     Note: the resolution is the dimension of a square grid of sub-cells, like:
c      ilkresolution(i,j)=0  cell is to be excluded from the matrix (initial leakage value to be kept)
c      ilkresolution(i,j)=1  cell has NO subgrid, only one collocation point at the center.
c      ilkresolution(i,j)=2  cell has a 2*2 subgrid, hence 4 sub-cells.
c      ilkresolution(i,j)=3  cell has a 3*3 subgrid, hence 9 sub-cells.
c                 etc.
c
c      In sizing the arrays for sub-grids, a resolution of 1 is interpreted as 1 sub-cell that equals
c      the MODFLOW cells.
c      NOTE: This storage space is probably not needed and may be eliminated at some point.
c
      implicit none
      INTEGER nsolOut,i,j,ilu,itemp,ivar,ierr,ilen,icharacter,
     &        nlkcolsize,nlkrowsize,nlkbufsize,
     &        iresolution
      LOGICAL lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut,
     &        lnofile
      REAL(8) rdum,rvar
      CHARACTER*8 aBasenameOut
      CHARACTER*16 aDateTimeOut
      ALLOCATABLE iresolution(:)
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
c
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      nlkbufsize=0
      ilu=2
      afile=TRIM(aBasenameOut)//".grs"
      itemp=iluin ! temporarily store iluin and replace for this read
      iluin=ilu
      OPEN (UNIT=ILU,FILE=AFILE,STATUS='OLD',IOSTAT=IERR)
      if (ierr.ne.0) then  ! failed to open file
        lnofile=.TRUE.
        close (iluin)
        iluin=itemp ! restore iluin
        return
      endif
      lnofile=.FALSE.
c      write (iluer,1001)
c 1001 format (' lk_scan_filename_grs has been successfully opened.')
c
c     read resolution numbers row by row and sum to assess size of subgrid buffers.
c
      ALLOCATE (iresolution(nlkcolsize),stat=ierr)
      if (ierr.ne.0) then
       write (iluer,2000) nlkcolsize
        AMESS(1)='Failed to allocate array: iresolution.'
        AMESS(2)='Verify resources and rerun.'
        CALL HALT(2) ! stop program execution for batch version
      end if
c
      do i=1,nlkrowsize
       read (ilu,*) iresolution ! list directed read expects a set of integers in a record that is equal
c                                 to the dimension of the array iresolution
c       write (*,1002) iresolution(1),iresolution(2),iresolution(19)
c 1002 format (' lk_scan_filename_grs2: iresolution(1,2,19) ',3i5)
       do j=1,nlkcolsize
        if (iresolution(j).gt.1) then
          nlkbufsize=nlkbufsize+iresolution(j)*iresolution(j)
        end if
       end do
      end do
c
      close (iluin)
      iluin=itemp ! restore iluin
      return
c
 1000 format (' ERROR in lk_scan_filename_grs:',/,
     &        ' The file ',a16,' is not found.',/,
     &        ' program execution is aborted.')
 2000 format (' ERROR in lk_scan_filename_grs:',/,
     &        ' Allocation of iresolution failed.',/,
     & ' array size = ',I5,' Program execution is aborted.')
 3000 format (i3)
      end subroutine
c ---------------------------------------------------------------------------------
c
      subroutine lk_read_filename_grs (nrow,ncol,ilkresolution,
     &                                 ilkpointer,lnofile)
c
c     Routine reads the sub-grid resolutions provided by the GFLOW GUI.
c     The routine also creates a pointer array "ilkpointer" to provide the starting address
c     for a MODFLOW cell of data in sub-cell arrays. If there are no sub-cells in MODFLOW cell i,j
c     then the value of ilkpointer(i,j) is set to zero.

c
      implicit none
      INTEGER nsolOut,i,j,ilu,itemp,ivar,ierr,ilen,icharacter,
     &        ilkresolution,ilkpointer,nrow,ncol,ibuf,ires
      LOGICAL lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut,lnofile
      REAL(8) rdum,rvar
      CHARACTER*8 aBasenameOut
      CHARACTER*16 aDateTimeOut
      DIMENSION ilkpointer(nlkrowsize,nlkcolsize)
      DIMENSION ilkresolution(nrow,ncol)
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
c
      lknosubgridfile=lnofile ! set this in lkcom.inc
      if (lknosubgridfile) then
        ilkresolution(1:nlkrowsize,1:nlkcolsize)=0
        ilkpointer(1:nlkrowsize,1:nlkcolsize)=0
      else   ! read subgrid resolutions and assign pointers
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      ilu=2
      afile=TRIM(aBasenameOut)//".grs"
      itemp=iluin ! temporarily store iluin and replace for this read
      iluin=ilu
      OPEN (UNIT=ILU,FILE=AFILE,STATUS='OLD',IOSTAT=IERR)
      if (ierr.ne.0) then  ! failed to open file
        write (iluer,1000) afile
        AMESS(1)='The file basename.grs is mot found.'
        AMESS(2)='Correct and rerun.'
        CALL HALT(2) ! stop program execution for batch version
      endif
c
c     read resolution numbers row by row and store in ilkresolution
c
c
      ibuf=1
      do i=1,nlkrowsize
       read (ilu,*) ilkresolution(i,1:nlkcolsize) ! list directed read expects a set of integers in a record that is equal
c                                 to the dimension of the array iresolution
c       write (*,1002) iresolution(1),iresolution(2),iresolution(19)
c 1002 format (' lk_scan_filename_grs2: iresolution(1,2,19) ',3i5)
       do j=1,nlkcolsize
        ires=ilkresolution(i,j)
        ilkpointer(i,j)=0
        if (ires.gt.1) then
          ilkpointer(i,j)=ibuf
          ibuf=ibuf+ires*ires
        end if
       end do
      end do
c
      close (iluin)
      iluin=itemp ! restore iluin
      ENDIF  ! end of lknosubgridfile loop
      return
c
 1000 format (' ERROR in lk_read_filename_grs:',/,
     &        ' The file ',a16,' is not found.',/,
     &        ' program execution is aborted.')
      end subroutine
c
c ---------------------------------------------------------------------------------
c
      subroutine lk_read_filename_fmt (nrow,ncol,rlkheadlower)
c
c     Routine reads the heads in the lower aquifer (MODFLOW layer underneath the GFLOW aquifer).
c
      implicit none
      INTEGER nsolOut,i,j,ilu,itemp,ivar,ierr,ilen,icharacter,
     &        nrow,ncol
      LOGICAL lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut
      REAL(8) rdum,rvar,rlkheadlower
      CHARACTER*8 aBasenameOut
      CHARACTER*16 aDateTimeOut
      DIMENSION rlkheadlower(nrow,ncol)
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
c
      if (.not.lksolving) then
       rlkheadlower(1:nrow,1:ncol)=0.0d0
       return         ! default lower heads and return. Do not read .fmt file.
      end if

      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      ilu=2
      afile=TRIM(aBasenameOut)//".fmt"
      itemp=iluin ! temporarily store iluin and replace for this read
      iluin=ilu
      OPEN (UNIT=ILU,FILE=AFILE,STATUS='OLD',IOSTAT=IERR)
      if (ierr.ne.0) then  ! failed to open file
        write (iluer,1000) afile
        AMESS(1)='The file basename.fmt is mot found.'
        AMESS(2)='Correct and rerun.'
        CALL HALT(2) ! stop program execution for batch version
      endif
c
c     read lower heads row by row and store
c
c
      do i=1,nlkrowsize
       read (ilu,*) rlkheadlower(i,1:nlkcolsize) ! list directed read expects a set of reals in a record that is equal
c                                 to the first dimension of the array rlkheadlower
      end do
c
      close (iluin)
      iluin=itemp ! restore iluin
c
c     ! check for "dry cells" and replace flag (-900) by aquifer top
c
      call lk_read_filename_top (nrow,ncol,rlkheadlower)
c
c
c     debugging
c
c      do i=1,nlkrowsize
c      do j=1,nlkcolsize
c      if (i.eq.20.and.j.eq.18) THEN ! cell with the well present
c         write (iluer,1001) i,j,rlkheadlower(i,j)
c      end if
c      enddo
c      enddo
c 1001 format(' lk_read_filename_fmt1: rlkheadlower(',i3,',',i3,')='
c     &       ,d14.7)
c

      return
c
 1000 format (' ERROR in lk_read_filename_fmt:',/,
     &        ' The file ',a16,' is not found.',/,
     &        ' program execution is aborted.')
      end subroutine
c
c ---------------------------------------------------------------------------------
c
      subroutine lk_read_filename_rla (nrow,ncol,rlkresist)
c
c     This routine reads the resistance for the (hypothetical) resistance layer based
c     on the aquifer properties of the lower (MODFLOW layer underneath the GFLOW aquifer).
c
      implicit none
      INTEGER nsolOut,i,j,ilu,itemp,ivar,ierr,ilen,icharacter,
     &        nrow,ncol
      LOGICAL lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut
      REAL(8) rdum,rvar,rlkresist
      CHARACTER*8 aBasenameOut
      CHARACTER*16 aDateTimeOut
      DIMENSION rlkresist(nrow,ncol)
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
c
      if (.not.lksolving) then
       do i=1,nlkrowsize
        do j=1,nlkcolsize
         rlkresist(i,j)=0.0d0
        end do
       end do
       return         ! default resistances and return. Do not read .rla file.
      end if

      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      ilu=2
      afile=TRIM(aBasenameOut)//".rla"
      itemp=iluin ! temporarily store iluin and replace for this read
      iluin=ilu
      OPEN (UNIT=ILU,FILE=AFILE,STATUS='OLD',IOSTAT=IERR)
      if (ierr.ne.0) then  ! failed to open file
        write (iluer,1000) afile
        AMESS(1)='The file basename.rla is mot found.'
        AMESS(2)='Correct and rerun.'
        CALL HALT(2) ! stop program execution for batch version
      endif
c
c     read resistances row by row and store
c
c
      do i=1,nlkrowsize
       read (ilu,*) rlkresist(i,1:nlkcolsize) ! list directed read expects a set of reals in a record that is equal
c                                 to the first dimension of the array rlkresist
      end do
c
      close (iluin)
      iluin=itemp ! restore iluin
      return
c
 1000 format (' ERROR in lk_read_filename_rla:',/,
     &        ' The file ',a16,' is not found.',/,
     &        ' program execution is aborted.')
      end subroutine
c
c ---------------------------------------------------------------------------------
c
      subroutine lk_read_filename_rst (nrow,ncol,rlkresist)
c
c     This routine reads the resistance for the (hypothetical) resistance layer based
c     on the aquifer properties of both the lower (MODFLOW) layer and the GFLOW aquifer.
c
      implicit none
      INTEGER nsolOut,i,j,ilu,itemp,ivar,ierr,ilen,icharacter,
     &        nrow,ncol
      LOGICAL lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut
      REAL(8) rdum,rvar,rlkresist
      CHARACTER*8 aBasenameOut
      CHARACTER*16 aDateTimeOut
      DIMENSION rlkresist(nrow,ncol)
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
c
      if (.not.lksolving) then
       rlkresist(1:nrow,1:ncol)=0.0d0
       return         ! default lower heads and return. Do not read .fmt file.
      end if

      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      ilu=2
      afile=TRIM(aBasenameOut)//".rst"
      itemp=iluin ! temporarily store iluin and replace for this read
      iluin=ilu
      OPEN (UNIT=ILU,FILE=AFILE,STATUS='OLD',IOSTAT=IERR)
      if (ierr.ne.0) then  ! failed to open file
        write (iluer,1000) afile
        AMESS(1)='The file basename.rst is mot found.'
        AMESS(2)='Correct and rerun.'
        CALL HALT(2) ! stop program execution for batch version
      endif
c
c     read resistances row by row and store
c
c
      do i=1,nlkrowsize
       read (ilu,*) rlkresist(i,1:nlkcolsize) ! list directed read expects a set of reals in a record that is equal
c                                 to the first dimension of the array rlkresist
      end do
c
      close (iluin)
      iluin=itemp ! restore iluin
      return
c
 1000 format (' ERROR in lk_read_filename_rst:',/,
     &        ' The file ',a16,' is not found.',/,
     &        ' program execution is aborted.')
      end subroutine
c
c ---------------------------------------------------------------------------------
c
      subroutine lk_read_filename_vlk (nrow,ncol,rlkleakage)
c
c     This routine reads the initial leakage into the lower aquifer
c     (MODFLOW layer underneath the GFLOW aquifer). This leakage is provided
c     by MODFLOW (based on effective heads and effective leakages) provided by
c     the GFLOW GUI
c
      implicit none
      INTEGER nsolOut,i,j,ilu,itemp,ivar,ierr,ilen,icharacter,
     &        nrow,ncol
      LOGICAL lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut
      REAL(8) rdum,rvar,rlkleakage
      CHARACTER*8 aBasenameOut
      CHARACTER*16 aDateTimeOut
      DIMENSION rlkleakage(nrow,ncol)
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
c
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      ilu=2
      afile=TRIM(aBasenameOut)//".vlk"
      itemp=iluin ! temporarily store iluin and replace for this read
      iluin=ilu
      OPEN (UNIT=ILU,FILE=AFILE,STATUS='OLD',IOSTAT=IERR)
      if (ierr.ne.0) then  ! failed to open file
        write (iluer,1000) afile
        AMESS(1)='The file basename.vlk is mot found.'
        AMESS(2)='Correct and rerun.'
        CALL HALT(2) ! stop program execution for batch version
      endif
c
c     read leakages row by row and store
c
c
      do i=1,nlkrowsize
       read (ilu,*) rlkleakage(i,1:nlkcolsize) ! list directed read
      end do
c
      close (iluin)
      iluin=itemp ! restore iluin
c
      return
c
 1000 format (' ERROR in lk_read_filename_vlk:',/,
     &        ' The file ',a16,' is not found.',/,
     &        ' program execution is aborted.')
      end subroutine
c
c --------------------------------------------------------------------------------
c
      subroutine lkreadnewleakages_actual(nrow,ncol,
     &                         ilkresolution,rlkleakage,drb,jad,m)
c
c     This routine is called once in SOLUT after a solution has been read from disk and
c     new additional matrix solutions are to be created. Since we may have had a MODFLOW solution
c     preceed this, we must read the updated leakages in active cells and store the difference in
c     the array DRB, similarly as done for modificatons of the line-sink strengths.
c     Note: MODFLOW should not modify leakages in cells with a zero resolution, thus inactive for the GFLOW-
c           MODFLOW interactions. Either way, GFLOW will not read leakages in inactive cells upon a
c           load *.sol command.
c
      implicit none
      INTEGER nsolOut,i,j,ilu,itemp,ivar,ierr,ilen,icharacter,
     &        nrow,ncol,ilkresolution,jad,m
      LOGICAL lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut
      REAL(8) rdum,rvar,rlkleakage,rtemp,drb
      CHARACTER*8 aBasenameOut
      CHARACTER*16 aDateTimeOut
      DIMENSION rlkleakage(nrow,ncol),ilkresolution(nrow,ncol),
     &          rtemp(ncol),drb(m)
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
c
      if (lkleakage.and.lksolving) then   ! only do so if we are actually solving for leakages
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      ilu=2
      afile=TRIM(aBasenameOut)//".vlk"
      itemp=iluin ! temporarily store iluin and replace for this read
      iluin=ilu
      OPEN (UNIT=ILU,FILE=AFILE,STATUS='OLD',IOSTAT=IERR)
      if (ierr.ne.0) then  ! failed to open file
        write (iluer,1000) afile
        AMESS(1)='The file basename.vlk is mot found.'
        AMESS(2)='Correct and rerun.'
        CALL HALT(2) ! stop program execution for batch version
      endif
c
c     read leakages row by row and store
c
c
      do i=1,nlkrowsize
       read (ilu,*) rtemp(1:nlkcolsize) ! list directed read
       do j=1,nlkcolsize
         if (ilkresolution(i,j).gt.0) then ! only replace leakage value for active cells
           jad=jad+1
           drb(jad)=rtemp(j)-rlkleakage(i,j) ! store the difference in drb
c           write (iluer,1001) jad,drb(jad),rtemp(j),rlkleakage(i,j)
c 1001 format (' lkreadnewleakages1:',
c     &' jad,drb(jad),rtemp(j),rlkleakage(i,j)',/,i5,3(2x,d14.7))
           rlkleakage(i,j)=rtemp(j)          ! now store the new leakage value in rlkleakage(j,j)
         end if
       end do
      end do
c
      close (iluin)
      iluin=itemp ! restore iluin
      endif
c
      return
c
 1000 format (' ERROR in lkreadleakages:',/,
     &        ' The file ',a16,' is not found.',/,
     &        ' program execution is aborted.')
      end subroutine
c
c ---------------------------------------------------------------------------------
c
      subroutine lk_read_filename_rta (nrow,ncol,rlkrecharge)
c
c     This routine reads the recharge values used by MODFLOW for the upper layer.
c
c     *.rta   "Recharge at Top of Aquifer"
c
c     The input of recharge values is optional, if the file "basename.rta" is not
c     found a warning will be issued and the recharge rates will be set to 0.0
c
      implicit none
      INTEGER nsolOut,i,j,ilu,itemp,ivar,ierr,ilen,icharacter,
     &        nrow,ncol
      LOGICAL lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut
      REAL(8) rdum,rvar,rlkrecharge
      CHARACTER*8 aBasenameOut
      CHARACTER*16 aDateTimeOut
      DIMENSION rlkrecharge(nrow,ncol)
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
c
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      ilu=2
      afile=TRIM(aBasenameOut)//".rta"
      itemp=iluin ! temporarily store iluin and replace for this read
      iluin=ilu
      OPEN (UNIT=ILU,FILE=AFILE,STATUS='OLD',IOSTAT=IERR)
      if (ierr.ne.0) then  ! failed to open file
        write (iluer,1000) afile
        rlkrecharge(1:nlkrowsize,1:nlkcolsize)=0.0d0 ! set recharge to 0.0
      else
c
c     read leakages row by row and store
c
c
        do i=1,nlkrowsize
         read (ilu,*) rlkrecharge(i,1:nlkcolsize) ! list directed read expects a set of reals in a record that is equal
c                                 to the first dimension of the array rlkrecharge
        end do
c
        close (iluin)
      endif
      iluin=itemp ! restore iluin
c
c     Reverse the sign since recharges in MODFLOW are positive, while GFLOW will use them as sink densities.
c
      do i=1,nlkrowsize
        do j=1,nlkcolsize
          rlkrecharge(i,j)=-rlkrecharge(i,j)
        end do
      end do
      return
c
 1000 format (' Warning in lk_read_filename_rta:',/,
     &        ' The file ',a16,' is not found.',/,
     &        ' Recharge in grid is set to 0.0.')
      end subroutine
c
c ---------------------------------------------------------------------------------
c
      subroutine lk_read_filename_top (nrow,ncol,rlkheadlower)
c
c     Routine reads the aquifer top elevations of the lower aquifer (MODFLOW layer underneath the GFLOW aquifer).
c     All flags -900 in the lower head array (dry MODFLOW cells) are replaced by the aquifer top of that cell.
c     Routine is called in lk_read_filename_fmt
c
      implicit none
      INTEGER nsolOut,i,j,ilu,itemp,ivar,ierr,ilen,icharacter,
     &        nrow,ncol
      LOGICAL lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut
      REAL(8) rdum,rvar,rlkheadlower,scratch
      CHARACTER*8 aBasenameOut
      CHARACTER*16 aDateTimeOut
      DIMENSION rlkheadlower(nrow,ncol)
      ALLOCATABLE scratch(:)
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
c
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      ilu=2
      afile=TRIM(aBasenameOut)//".top"
      itemp=iluin ! temporarily store iluin and replace for this read
      iluin=ilu
      OPEN (UNIT=ILU,FILE=AFILE,STATUS='OLD',IOSTAT=IERR)
      if (ierr.ne.0) then  ! failed to open file
        write (iluer,1000) afile
        AMESS(1)='The file basename.top is mot found.'
        AMESS(2)='Correct and rerun.'
        CALL HALT(2) ! stop program execution for batch version
      endif
c
c     read lower heads row by row and store
c
c
      ALLOCATE (scratch(ncol))
c
      do i=1,nrow
       read (ilu,*) scratch(1:ncol) ! list directed read expects a set of reals in a record that is equal
c                                 to the dimension of the array scratch
       do j=1,ncol
         if (rlkheadlower(i,j).EQ.-900.0) THEN ! dry cell, replace flag by aquifer top
           rlkheadlower(i,j)=scratch(j)
           write (iluer,2000) i,j,scratch(j)
         end if
       end do
      end do
c
      close (iluin)
      iluin=itemp ! restore iluin
      deallocate (scratch)
      return
c
 1000 format (' ERROR in lk_read_filename_top:',/,
     &        ' The file ',a16,' is not found.',/,
     &        ' program execution is aborted.')
 2000 format (' NOTE: cell(',i3,',',i3,') is dry. Flag replaced by',/,
     & ' aquifer top elevation =',d14.7)
      end subroutine
c
c ---------------------------------------------------------------------------------
c
      subroutine lk_write_filename_ler (nrow,ncol,rlkpercenterror)
c
c     This routine writes the errors at collocation points to a file
c     Routine is called SOLUT in GFMOD.FOR
c
      implicit none
      INTEGER nsolOut,i,j,ilu,itemp,ivar,ierr,ilen,icharacter,
     &        nrow,ncol
      LOGICAL lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut
      REAL(8) rdum,rvar,rlkpercenterror
      CHARACTER*8 aBasenameOut
      CHARACTER*16 aDateTimeOut
      DIMENSION rlkpercenterror(nrow,ncol)
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
c
      if (.not.lksolving) RETURN ! do not write leakage errors if no leakage solution
c
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      ilu=2
      afile=TRIM(aBasenameOut)//".ler"
      OPEN (UNIT=ILU,FILE=AFILE,STATUS='UNKNOWN',IOSTAT=IERR)
      if (ierr.ne.0) then  ! failed to open file
        write (iluer,1000) afile
        return
      endif
c
c     write the percent errors row by row
c
c
      do i=1,nlkrowsize
       write (ilu,*) rlkpercenterror(i,1:nlkcolsize) ! list directed write
      end do
c
      close (ilu)
      return
c
 1000 format (' ERROR in lk_write_filename_ler:',/,
     &        ' Could not open or create file ',A16,/,
     &        ' No errors have been written.')
      end subroutine
c
c -----------------------------------------------------------------------------------
c
      subroutine lkgfmfmonitor_actual (nrow,ncol,
     &                   rlkheadlower,rlkheadupper,rlkleakage)
c
c     Monitor conjunctive GFLOW - MODFLOW solution performance
c     Original routine lkgfmfmonitor is called in Setup at end of solution process.
c
      implicit none
      INTEGER nrow,ncol,i,j,idt
      CHARACTER (LEN=10) :: atime,adate,azone
      REAL(8) rlkheadlower,rlkheadupper,rlkleakage
      DIMENSION rlkheadlower(nrow,ncol),rlkheadupper(nrow,ncol),
     &                    rlkleakage(nrow,ncol),idt(8)
      INCLUDE 'lusys.inc'
      INCLUDE 'lkcom.inc'
c
      if (.not.lksolving) RETURN ! Do not write monitoring file if no leakage solution, hence no iterations between GFLOW and MODFLOW
c

      i=8    ! to be adjusted as desired
      j=9    ! to be adjusted as desired
      call date_and_time (adate,atime,azone,idt)
      write (ilutmp,2000) adate,atime,i,j
      if (i.lt.nrow.and.j.lt.ncol) THEN ! make sure we are inside the current MODFLOW grid
      write (iluer,1001)
 1001 format ('lkgfmfmonitor_actual1: writing to gfmfmonitor.dat')
      write (ilutmp,1000) rlkheadlower(i,j),rlkheadupper(i,j),
     &                    rlkleakage(i,j)
      end if
c
      return
 1000 format (' hlower',d14.7,' hupper',d17.7,' leakage',d14.7)
 2000 format (' date: ',a10,' time: ',a10,' i,j=',2(i4))
      end subroutine


