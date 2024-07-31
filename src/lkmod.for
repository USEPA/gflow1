C     Last change:  HMH  16 Apr 2015    2:54 am
c  lkdat            block data routine to initialize common LKCOM
c  lkin             recursive subroutine to enter data and pass allocatable arrays through entry statements
c
c  ENTRY statements in LKIN:
c  lk_cell_construct   calls lk_cell_construct_actual to provide four corner points of the cell for given row and column
c  lkfindcell          calls lkfindcell_actual to provide the row and column of the cell in which a point occurs
c  cflk_omega_sub      calls cflk_omega_actual to return the complex potential for the entire grid (note, Poisson contribution added to real part)
c  cflk_W_sub          calls cflk_W_actual to return the W-function (negative derivative of complex potential)
c  lkio_array          calls lkio_array_actual to write the allocatable arrays to basefilename.sol
c  rflklowerhead_sub   call rflklowerhead_sub_actual to report lower heads at cell centers
c  rflkleakage_sub     calls rflkleakage_sub_actual to report leakages
c  rflkresistance_sub  calls rflkresistance_sub_actual to report resistances
c  lksolve             calls lksolve_actual to solve for leakages using Gauss-Seidel
c  lkwrite_leakage     call lk_write_filename_vlk writes leakage values back to "basefilename.vlk"
c  lkwrite_error       calls lkset_check to write errors in leakage to the "basefilename.ler"
c  lkczc               calls lkczc_actual to make control point array
c
c  additional subroutines:
c  lksuborigin      substitute origin parameters in common after a new grid is being defined.
c  lkset_solve      substitute solution parameters in common
c  lksetconst       calculate constants for the line-sink potential so it vanishes at the reference point
c  lflksolving           passing lksolving
c  lflkleakage            passing lkleakage
c  lksetsubheadlower     assign interpolated lower heads to rlksubheadlower(nbuf)
c  lkinterpolateheads    interpolation routine for lower heads
c
      block data lkdat
      INCLUDE 'lkcom.inc'
      DATA nlkcolsize,nlkrowsize,nlkbufsize,nlkiter/4*0/
      DATA nlkterms /10/
      DATA ilkiniter /5/ ! iteration at which to start subgrid calculations if specified.
      DATA lkleakage,lkincludesubgrid,lkcheck,lksolving,lkexclude,
     &     lknosubgridfile
     &        /.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE./
      DATA rlkrelax /1.0d0/
      DATA rQ_leakage,rQ_target_leakage /0.0,0.0/
      DATA ci /(0.0,1.0)/
c
      end
c
c -------------------------------------------------------------------
      recursive subroutine LKIN (LSOL,m,dra,drscr,drb,drfac,
     &                           czi,calph,itype,nstr)
c
c     Routine reads MODFLOW grid structure and GFLOW subgrid information
c     to set up a leakege grid for coupling GFLOW and MODFLOW
c
c     This recursive routine is also the basis for entries to pass
c     the adjustable arrays.
c
      implicit none
      INTEGER i,j,m,n,nlkcolsize,nlkrowsize,nlkbufsize,ilkresolution,
     &        imatsize,ierr,idum,jdum,ilu,icode,nrow,ncol,nbuf,mout,
     &        i0,j0,jump,idum1,idum2,idum3,idum4,idum5,ivar,itype,
     &        iloc,jloc,nsolOut,ilkpointer,ninner,nouter,isubsign,mtot,
     &        imxsze,nstr
      LOGICAL lsol,lreadsol,l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,
     &        l14,l15,l16,lksubinclude,lkinout,
     &        lflksolving,lfleakage,lDirectFromDisk,ltimer,lskip,lbad,
     &        lsolOut,loadsolOut,linalreadyOut,lErrorReportOut,
     &        lDirectfromDiskOut,lnofile
      CHARACTER(1) AWORD(21)
      CHARACTER(16) aBasenameOut,aDateTimeOut
      REAL(8) rlkdeltax,rlkdeltay,rlkheadlower,rlkleakage,rlkheadupper,
     &        rlkpot,rlkresist,rlkrecharge,rlksubpotupper,
     &        rlkconst,rlksubheadlower,rlksubleakage,rlkpercenterror,
     &        rdum1,rdum2,rvar,rx0,ry0,DRSCR,rermax,rlksublambda2init,
     &        rscr1,rfunc,rversion,rq3,rqi,dra,drfac,drb,rdy,rdx
      COMPLEX(8) cz,cfunc,cz0,cz1,cz2,cz3,czi,calph
      DIMENSION DRA(m,*),rqi(3),drfac(4,*),lskip(*),           ! forced DRA with fixed dimension 1
     &          czi(*),drb(*),calph(*),itype(*),
     &          DRSCR(*)
      ALLOCATABLE rlkdeltax(:),rlkdeltay(:),rlkheadlower(:,:),
     &            rlkleakage(:,:),rlkpot(:,:),rlkheadupper(:,:),
     &            rlkresist(:,:),rlkconst(:,:),rlkpercenterror(:,:),
     &            ilkresolution(:,:),rlkrecharge(:,:),
     &            rlksubheadlower(:),rlksubleakage(:),
     &            rlksubpotupper(:),rlksublambda2init(:),ilkpointer(:,:)
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
      DATA AWORD/ 'O','R','I','G',' ',
     .            'S','O','L','V',' ',
     .            'C','H','E','C',' ',
     .            'Q','U','I','T',' ',
     .            ATERM/
      save      ! IMPORTANT TO KEEP THE ALLOCATABLE ARRAYS!!!!
C
      LERROR=.FALSE.
      LMISS=.FALSE.
 10   IF (LERROR.OR.LMISS) WRITE(ILUER,9001) ALINE2
      LERROR=.FALSE.
      LMISS=.FALSE.
      CALL INLINE
 11   CALL MATCH(AWORD,1,JUMP,LBAD)
      IF(.NOT.LBAD) GOTO 15
      GOTO (10,13), JUMP
 13   WRITE(ILUER,9002) ALINE2
      LERROR=.FALSE.
      GOTO 10
 15   GOTO(100,200,300,400),JUMP
c
c  ------------------------------------------------ start of origin command --------------------
c
c
c     read origin of grid (read as coordinates of the LOWER left corner of the grid)
c
c     NOTE: THIS COMMAND MUST BE THE FIRST COMMAND GIVEN IN THE LEAKAGE MODULE!!
c
c
c
  100 continue
      rdum1=rvar(2)
      rdum2=rvar(3)
      if (lerror.or.lmiss) then
        write (iluer,5000)
        return
      end if ! values will be assigned after the grid data has been read in.
c
c     scan files for array dimensions
c
      call lk_scan_filename_del(nlkrowsize,nlkcolsize)
      call lk_scan_filename_grs(nlkrowsize,nlkcolsize,nlkbufsize,
     &                          lnofile)
      nlkrowsize=MAX(nlkrowsize,1)
      nlkcolsize=MAX(nlkcolsize,1)
      nlkbufsize=MAX(nlkbufsize,1)
c
      lreadsol=.false.
  50  continue
c
c     deallocate any existing arrays
c
      l1=ALLOCATED(rlkdeltax)
      l2=ALLOCATED(rlkdeltay)
      l3=ALLOCATED(rlkheadlower)
      l4=ALLOCATED(rlkleakage)
      l5=ALLOCATED(rlkresist)
      l6=ALLOCATED(rlkconst)
      l7=ALLOCATED(ilkresolution)
      l8=ALLOCATED(rlksubheadlower)
      l9=ALLOCATED(rlksubleakage)
      l10=ALLOCATED(rlkpot)
      l11=ALLOCATED(rlkrecharge)
      l12=ALLOCATED(rlksubpotupper)
      l13=ALLOCATED(rlkheadupper)
      l14=ALLOCATED(rlkpercenterror)
      l15=ALLOCATED(ilkpointer)
      l16=ALLOCATED(rlksublambda2init)
      if (l1) deallocate(rlkdeltax)
      if (l2) deallocate(rlkdeltay)
      if (l3) deallocate(rlkheadlower)
      if (l4) deallocate(rlkleakage)
      if (l5) deallocate(rlkresist)
      if (l6) deallocate(rlkconst)
      if (l7) deallocate(ilkresolution)
      if (l8) deallocate(rlksubheadlower)
      if (l9) deallocate(rlksubleakage)
      if (l10) deallocate(rlkpot)
      if (l11) deallocate(rlkrecharge)
      if (l12) deallocate(rlksubpotupper)
      if (l13) DEALLOCATE(rlkheadupper)
      if (l14) DEALLOCATE(rlkpercenterror)
      if (l15) deallocate(ilkpointer)
      if (l16) DEALLOCATE(rlksublambda2init)
c
c     allocate arrays
c
      ALLOCATE (
     & rlkdeltay(nlkrowsize),
     & rlkdeltax(nlkcolsize),
     & rlkheadlower(nlkrowsize,nlkcolsize),
     & rlkheadupper(nlkrowsize,nlkcolsize),
     & rlkpot(nlkrowsize,nlkcolsize),
     & rlkleakage(nlkrowsize,nlkcolsize),
     & rlkresist(nlkrowsize,nlkcolsize),
     & rlkconst(nlkrowsize,nlkcolsize+1),
     & ilkresolution(nlkrowsize,nlkcolsize),
     & rlkrecharge(nlkrowsize,nlkcolsize),
     & rlkpercenterror(nlkrowsize,nlkcolsize),
     & ilkpointer(nlkrowsize,nlkcolsize),
     & stat=ierr)
      if (ierr.ne.0) then
       write (iluer,1000) nlkrowsize,nlkcolsize
        AMESS(1)='Failed to allocate MODFLOW grid arrays.'
        AMESS(2)='Verify resources and rerun.'
        CALL HALT(2) ! stop program execution for batch version
      end if
       ALLOCATE (rlksubheadlower(nlkbufsize),rlksubleakage(nlkbufsize),
     & rlksubpotupper(nlkbufsize),rlksublambda2init(nlkbufsize),
     & stat=ierr)
       if (ierr.ne.0) then
        write (iluer,2000) nlkbufsize
        AMESS(1)='Failed to allocate subgrid buffers.'
        AMESS(2)='Verify resources and rerun.'
        CALL HALT(2) ! stop program execution for batch version
       end if
       if (nlkbufsize.gt.0) then
        rlksubleakage(1:nlkbufsize)=0.0d0 ! initialize sub-leakages to zero
       end if
c             default all arrays to zero
c
      rlkdeltay(1:nlkrowsize)=0.0d0
      rlkdeltax(1:nlkcolsize)=0.0d0
      rlkheadlower(1:nlkrowsize,1:nlkcolsize)=0.0d0
      rlkheadupper(1:nlkrowsize,1:nlkcolsize)=0.0d0
      rlkpot(1:nlkrowsize,1:nlkcolsize)=0.0d0
      rlkleakage(1:nlkrowsize,1:nlkcolsize)=0.0d0
      rlkresist(1:nlkrowsize,1:nlkcolsize)=0.0d0
      rlkconst(1:nlkrowsize,1:nlkcolsize+1)=0.0d0
      ilkresolution(1:nlkrowsize,1:nlkcolsize)=0
      rlkrecharge(1:nlkrowsize,1:nlkcolsize)=0.0d0
      rlkpercenterror(1:nlkrowsize,1:nlkcolsize)=0.0d0
      ilkpointer(1:nlkrowsize,1:nlkcolsize)=0
c
      if (lreadsol) RETURN ! exit point for entry "lkallocate"
c
c     read files for data
c
      call lk_read_filename_del
     &                   (nlkrowsize,nlkcolsize,rlkdeltax,rlkdeltay)
      call lk_read_filename_grs(nlkrowsize,nlkcolsize,ilkresolution,
     &                            ilkpointer,lnofile)
      call lk_read_filename_vlk(nlkrowsize,nlkcolsize,rlkleakage)
      call lk_read_filename_rta(nlkrowsize,nlkcolsize,rlkrecharge)
      call lksuborigin(rdum1,rdum2,
     &                 nlkrowsize,nlkcolsize,nlkbufsize,rlkdeltay)  ! now set grid origin
      call lksetconst (nlkrowsize,nlkcolsize,rlkdeltax,rlkdeltay,   ! next calculate constants
     &                 rlkconst)
      call lksetsubheadlower(nlkrowsize,nlkcolsize,rlkheadlower,
     &                       ilkresolution,
     &                       ilkpointer,rlksubheadlower,nlkbufsize)
c
c
c                  Start debugging
c
c
c
       GOTO 10 ! bypass debugging
c
c
      write (iluer,1122) nlkrowsize,nlkcolsize
      write (iluer,1123) (i,rlkdeltay(i),i=1,nlkrowsize)
      write (iluer,1124) (i,rlkdeltax(i),i=1,nlkcolsize)
 1122 format (//,' nlkrowsize,nlkcolsize ',2(i5,2x))
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
      write (iluer,1129) nlkrowsize,nlkcolsize
 1129 format (//,' rlkrecharge: nlkrowsize,nlkcolsize ',2(i5,2x))
      do i=1,nlkrowsize
      write (iluer,*) rlkrecharge(i,1:nlkcolsize)
      end do
      write (iluer,1130) nlkrowsize,nlkcolsize
 1130 format (//,' ilkpointer: nlkrowsize,nlkcolsize ',2(i5,2x))
      do i=1,nlkrowsize
      write (iluer,*) ilkpointer(i,1:nlkcolsize)
      end do
c
c                 End debugging
c
c
c
      GOTO 10 ! -------------------------------------- end of origin command ----------------
c
c ----------- Solve command ---------------------
c
 200  continue
      idum1=ivar(2) ! 1 for include in matrix solution
      idum2=ivar(3) ! 1 for calculate average upper head based on subgrid collocation points
      idum3=ivar(4) ! 1 for update resistances, starting with the resistance for the lower aquifer only
      idum4=ivar(5) ! 1 include iterations for the subgrids
      rdum1=rvar(6) ! relaxation parameter
      if (lerror.or.lmiss) then
       write (iluer,6000)
      end if
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      if (idum1.eq.1) then ! calculate leakages in matrix - load heads and resistances
      call lk_read_filename_fmt(nlkrowsize,nlkcolsize,rlkheadlower)
      call lk_read_filename_rst(nlkrowsize,nlkcolsize,rlkresist)
        call lksetsubheadlower(nlkrowsize,nlkcolsize,rlkheadlower,
     &                       ilkresolution,
     &                       ilkpointer,rlksubheadlower,nlkbufsize)
        IF (.not.loadsolOut)  call lk_read_filename_vlk(nlkrowsize,
     &                                         nlkcolsize,rlkleakage) ! else read selected leakages in SOLUT
      end if                                                          ! after lkkeepstrenghts has been called
      call lkset_solve (idum1,idum2,idum3,idum4,rdum1)
      if (idum3.eq.1) THEN ! read resistances due to lower aquifer parameters
       call lk_read_filename_rla (nlkrowsize,nlkcolsize,rlkresist)  ! thus overwrite the earlier read values from the .rst file
      end if
      if (idum4.eq.1.and.lnofile) then
        write (iluer,8000)
      end if
      GOTO 10
c
c ----------- Check command --------------------
c
 300  continue
      idum1=ivar(2)
      idum2=ivar(3)
      idum3=ivar(4)
      idum4=ivar(5)
      idum5=ivar(6)
      if (lerror.or.lmiss) then
       write (iluer,7000)
      end if
      call lk_read_filename_fmt(nlkrowsize,nlkcolsize,rlkheadlower) ! overwrite lower heads by the most recent values from MODFLOW (if conjunctive solution available)
      call lkset_check (idum1,idum2,idum3,idum4,idum5,
     & nlkrowsize,nlkcolsize,
     & rlkdeltax,rlkdeltay,rlkleakage,rlkheadlower,
     & rlkresist,ilkresolution,.false.)
      GOTO 10
c
c ----------- Quit command---------------------
c
 400  continue
      RETURN
c
c
c    ---------------------  ENTRY statements -------------------------
c
      ENTRY lkallocate (nrow,ncol,nbuf)
      nlkrowsize=nrow
      nlkcolsize=ncol
      nlkbufsize=nbuf
      lreadsol=.true.
      GOTO 50  ! will return at end of allocations.
c
      ENTRY lk_cell_construct(idum,jdum,rdx,rdy,cz0,cz1,cz2,cz3)   ! used in lkgenmat
      if (lfleakage()) then
      CALL  lk_cell_construct_actual(idum,jdum,rdx,rdy,cz0,cz1,cz2,cz3,
     &                                    rlkdeltax,rlkdeltay)
      end if
      return
c
      ENTRY lkfindcell (cz,i0,j0,rx0,ry0)  ! returns the cell in which CZ occurs
      call lkfindcell_actual (cz,i0,j0,rx0,ry0,
     &                        rlkdeltax,rlkdeltay,nlkrowsize,nlkcolsize)
      return
c
      ENTRY cflk_omega_sub (cz,cfunc)    ! calculate complex potential at CZ
      if (lfleakage()) then
      call cflk_omega_actual (cz,cfunc,nlkrowsize,nlkcolsize,
     &             rlkdeltax,rlkdeltay,rlkleakage,rlkrecharge,rlkconst)
      endif
      return
c
      ENTRY cflk_W_sub (cz,cfunc,rq3)    ! calculate W-function at CZ
      if (lfleakage()) then
      call cflk_W_actual (cz,cfunc,rq3,nlkrowsize,nlkcolsize,
     &                    rlkdeltax,rlkdeltay,rlkleakage,rlkrecharge)
      endif
      return
c
      ENTRY lkio_array (ICODE,ILU,RVERSION)  ! write leakage arrays to *.sol
      call lkio_array_actual (ICODE,ILU,RVERSION,nlkrowsize,nlkcolsize,
     &                 nlkbufsize,rlkdeltax,rlkdeltay,rlkheadlower,
     &                 rlkleakage,rlkresist,rlkconst,rlksubheadlower,
     &                 rlksubleakage,ilkresolution,rlkpot,rlkrecharge,
     &                 rlksubpotupper,ilkpointer,rlksublambda2init)
      return
c
      ENTRY rflklowerhead_sub (cz,rfunc) ! report the lower head
      rfunc=0.0
      if (lfleakage()) then
      call rflklowerhead_sub_actual (cz,rfunc,nlkrowsize,nlkcolsize,
     &                   nlkbufsize,ilkpointer,rlksubheadlower,
     &                   rlkdeltax,rlkdeltay,rlkheadlower,ilkresolution)
      endif
      return
c
c
      ENTRY rflkleakage_sub (cz,rfunc) ! report the MODFLOW cell leakage
      rfunc=0.0
      if (lfleakage()) then
      call rflkleakage_sub_actual (cz,rfunc,nlkrowsize,nlkcolsize,
     &                    rlkdeltax,rlkdeltay,rlkleakage)
      endif
      return
c
c
      ENTRY rflksubcellleakage_sub (cz,rfunc) ! report the sub-cell leakage
      rfunc=0.0
      if (lfleakage()) then
      call rflksubcellleakage_sub_actual (cz,rfunc,nlkrowsize,
     &                  nlkcolsize,nlkbufsize,ilkpointer,rlksubleakage,
     &                  rlkdeltax,rlkdeltay,ilkresolution)
      endif
      return
c
c
      ENTRY rflkrecharge_sub (cz,rfunc) ! report the recharge
      rfunc=0.0
      if (lfleakage()) then
      call rflkrecharge_sub_actual (cz,rfunc,nlkrowsize,nlkcolsize,
     &                    rlkdeltax,rlkdeltay,rlkrecharge)
      endif
      return
c
c

      ENTRY rflkresistance_sub (cz,rfunc) ! report the resistance
      rfunc=0.0
      if (lfleakage()) then
      call rflkresistance_sub_actual (cz,rfunc,nlkrowsize,nlkcolsize,
     &                    rlkdeltax,rlkdeltay,rlkresist)
      endif
      return
c
      ENTRY lksolve ()   ! solve for leakage using a Gauss-Seidel iterative scheme
      if (lfleakage()) then
      call lksolve_actual (nlkrowsize,nlkcolsize,nlkbufsize,
     &                 rlkdeltax,rlkdeltay,rlkheadlower,
     &                 rlkleakage,rlkresist,rlksubheadlower,
     &                 rlksubleakage,ilkresolution)
      endif
      return
c
      ENTRY lkwrite_leakage() ! write leakage values back to "basefilename.vlk"
      if (lfleakage()) then
      call lk_write_filename_vlk (nlkrowsize,nlkcolsize,rlkleakage,
     &                            ilkresolution)
      endif
      return
c
      ENTRY lkwrite_error()  ! write % errors in leakage at cell centers to "basefilename.ler"
      if (lfleakage().and.lflksolving()) then
       call lk_write_filename_ler(nlkrowsize,nlkcolsize,rlkpercenterror)
      endif
      return
c
      ENTRY lktotal_Q() ! calculate the total leakage for the grid
      if (lfleakage()) then
       call lkset_check (1,1,nlkrowsize,nlkcolsize,1,
     & nlkrowsize,nlkcolsize,
     & rlkdeltax,rlkdeltay,rlkleakage,rlkheadlower,
     & rlkresist,ilkresolution,.false.)
      endif
      return
c
      ENTRY lkczc(CZI,M,N,mtot,DRFAC,CALPH,ITYPE,
     &            lDirectFromDisk,ltimer)
      if (lfleakage()) then
       call lkczc_actual (CZI,M,N,mtot,DRFAC,CALPH,ITYPE,
     &             lDirectFromDisk,ltimer,
     &           nlkrowsize,nlkcolsize,rlkdeltax,rlkdeltay,
     &           ilkresolution,rlkleakage)
      endif
      return
c
      ENTRY lkgenmat (DRA,CZI,M,N,J,CALPH,ITYPE)
      if (lfleakage()) then
       call lkgenmat_actual (DRA,CZI,M,N,J,CALPH,ITYPE,
     &                  nlkrowsize,nlkcolsize,rlkconst,
     &                  ilkresolution)
      endif
      return
c
      ENTRY lkmatcorrect (DRA,CZI,M,N,J,CALPH,ITYPE)
      if (lfleakage()) then
       call lkmatcorrect_actual (DRA,CZI,M,N,J,CALPH,ITYPE,
     &    nlkrowsize,nlkcolsize,rlkconst,rlkresist,rlkpot,rlkheadlower,
     &    rlksubheadlower,ilkresolution,rlksubpotupper,
     &    rlksublambda2init,nlkbufsize)
      end if
      return
c
      ENTRY lkkno (DRB,J,CZI,n)
      if (lfleakage()) then
      call lkkno_actual (DRB,J,CZI,n,
     &      nlkrowsize,nlkcolsize,rlkconst,rlkresist,
     &      rlkpot,rlkleakage,rlkheadlower,
     &      ilkresolution,rlksubpotupper,rlksubheadlower,
     &      rlksublambda2init,rlksubleakage,nlkbufsize)
      end if
      return
c
      ENTRY lkupdate (DRSCR,j,M,N,czi,lksubinclude,isubsign)
      if (lfleakage()) then
      call lkupdate_actual (DRSCR,j,m,n,czi,lksubinclude,isubsign,
     &      nlkrowsize,nlkcolsize,rlkpot,ilkresolution,
     &      rlksubpotupper,rlkheadlower,rlkheadupper,nlkbufsize)
      end if
      return
c
      ENTRY lkupdate_check (j,m,n,czi)
      if (lfleakage()) then
         call lkupdate_check_actual (j,m,n,czi,
     &   nlkrowsize,nlkcolsize,rlkpot,ilkresolution,rlksubpotupper,nbuf)
      end if
      return
c
      ENTRY lksub (DRB,j)
      if (lfleakage()) then
       call lksub_actual (DRB,j,nlkrowsize,nlkcolsize,rlkleakage,
     &       ilkresolution)
      end if
      return
c
      ENTRY lkcombineequations (dra,czi,drfac,itype,lskip,m,mout,n)
      if (lfleakage()) then
       call lkcombineequations_actual (dra,czi,drfac,itype,lskip,m,mout,
     &      n,ilkresolution,nlkrowsize,nlkcolsize,rlksublambda2init,
     &      nlkbufsize)
      else
       mout=m
      end if
      return
c
      ENTRY lkerror (rermax)
      if (lfleakage()) then
      call lkerror_actual (rermax,
     &      nlkrowsize,nlkcolsize,rlkpot,rlkheadlower,rlkheadupper,
     &      rlkleakage,rlkresist,rlkdeltax,rlkdeltay,ilkresolution,
     &      rlkpercenterror,rlksubpotupper,rlksubheadlower,ilkpointer,
     &      rlksubleakage,nlkbufsize)
      end if
      return
c
      ENTRY lkwrite_upperheads()
      if (lfleakage()) then
       call lk_write_filename_uph (nlkrowsize,nlkcolsize,
     &                   rlkdeltax,rlkdeltay,rlkheadupper,ilkresolution)
      end if
      return
c
      ENTRY lkwrite_resistances()
      if (lfleakage()) then
        call lk_write_filename_rst (nlkrowsize,nlkcolsize,rlkresist)
      end if
      return
c
      ENTRY lkmatsize(M,N)
      if (lfleakage()) then
         call lkmatsize_actual(M,N,nlkrowsize,nlkcolsize,ilkresolution)
      end if
      return
c
      ENTRY lkreadnewleakages(drb,j,m)
      if (lfleakage()) then
         call lkreadnewleakages_actual(nlkrowsize,nlkcolsize,
     &                         ilkresolution,rlkleakage,drb,j,m)
      end if
      return
c
      ENTRY cflk_subomega_sub (cz,cfunc)
      if (lfleakage()) then
         call cflk_subomega_actual (cz,cfunc,nlkrowsize,nlkcolsize,
     &  nlkbufsize,rlkdeltax,rlkdeltay,ilkresolution,rlksubleakage,
     &  ilkpointer,rlkconst)
      end if
      return
c
      ENTRY lk_subW_sub (cz,cfunc,rq3)
      if (lfleakage()) then
         call lk_subW_actual (cz,cfunc,rq3,nlkrowsize,nlkcolsize,
     &  nlkbufsize,
     &  rlkdeltax,rlkdeltay,ilkresolution,rlksubleakage,ilkpointer)
      end if
      return
c
      ENTRY lkfn_sub (cz1,cz2,rfunc)
      if (lfleakage()) then
        call lkfnf_sub_actual (cz1,cz2,rfunc,nlkrowsize,nlkcolsize,
     &                       rlkdeltax,rlkdeltay,rlkleakage,rlkrecharge)
      end if
      return
c
      ENTRY lknfsub_sub (cz1,cz2,rfunc)
      if (lfleakage()) then
        call lkfnfsub_sub_actual (cz1,cz2,rfunc,nlkrowsize,nlkcolsize,
     &                           rlkdeltax,rlkdeltay,nlkbufsize,
     &                           rlksubleakage,ilkresolution,ilkpointer)
      end if
      return
c
c
c      ENTRY lksubgridsolve_sub (ninner,nouter,drb,n)
c      if (lfleakage()) then
c         call lksubgridsolve_actual1 (ninner,nouter,drb,n,nlkrowsize,
c     &         nlkcolsize,nlkbufsize,
c     &         rlkleakage,rlkresist,ilkresolution,
c     &         ilkpointer,rlksubpotupper,rlksubheadlower,rlksubleakage)
c      end if
c      return
c
      ENTRY lksubgridsolve_sub (nouter,
     &                          dra,drb,drscr,n,m,imxsze,nstr)
      if (lfleakage()) then
         call lksubgridsolve_actual2 (nouter,
     &         dra,drb,drscr,n,m,imxsze,nstr,
     &         nlkrowsize,nlkcolsize,nlkbufsize,
     &         rlkleakage,rlkresist,ilkresolution,
     &         ilkpointer,rlksubpotupper,rlksubheadlower,rlksubleakage,
     &         rlkconst)
      end if
      return
c
      ENTRY lkgfmfmonitor()
      if (lfleakage()) then
         call lkgfmfmonitor_actual (nlkrowsize,nlkcolsize,
     &          rlkheadlower,rlkheadupper,rlkleakage)
      end if
      return
c
 1000 format (' ERROR in LKIN at LKIN:',/,
     &        ' Allocation of MODFLOW grid arrays failed.',/,
     & ' nx, ny of MODFLOW grid = ',2(I5),/,
     & ' Program execution is aborted.')
 2000 format (' ERROR in LKIN at LKIN:',/,
     &        ' Allocation of sub-grid arrays failed.',/,
     & ' buffersize = ',I9,/,
     & ' Program execution is aborted.')
 5000 format (' ERROR in LKIN: no grid origin found.',/,
     &        ' No leakage grid applied.',/,
     &        ' Other Leakage commands are ignored.')
 6000 format (' ERROR in LKIN: missing or illegal Solve parameters.')
 7000 format (' ERROR in LKIN: missing or illegal Check parameters.')
 8000 format (' WARNING in LKIN: missing *.grs file - no subgrids.')
 9001 FORMAT(' ***ILLEGAL OR MISSING PARAMETERS in leakage module:',/,
     &       ' ',80A1,/)
 9002 FORMAT(' *** ILLEGAL COMMAND in leakage module:',/,
     &       ' ',80A1,/)
      end subroutine
c
c --------------------------------------------------------------------
c
      subroutine lksuborigin(rdum1,rdum2,
     &                       nrow,ncol,nbuf,rlkdeltay)
c
c     set grid origin and store array sizes in common
c
      implicit none
      INTEGER nrow,ncol,nbuf
      REAL(8) rdum1,rdum2,rlkdeltay
      DIMENSION rlkdeltay(nrow)
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
c
      lkleakage=.TRUE. ! authorize leakage calulations
      nlkrowsize=nrow
      nlkcolsize=ncol
      nlkbufsize=nbuf
      rlkx0=rdum1
      rlky0=rdum2+SUM(rlkdeltay(1:nrow)) ! redefine origin for UPPER left corner of grid
      return
      end subroutine
c
c --------------------------------------------------------------------
c
      subroutine lkset_solve (idum1,idum2,idum3,idum4,rdum1)
c
c     store solution parameters in common
c
      implicit none
      INTEGER idum1,idum2,idum3,idum4
      REAL(8) rdum1
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
c
      lksolving=idum1.eq.1
      lkaveragehead=idum2.eq.1
      lkupdateresistances=idum3.eq.1
      lkincludesubgrid=idum4.eq.1
      rlkrelax=rdum1
c
c      write (iluer,1001) lksolving,lkaveragehead,lkupdateresistances,
c     &                   lkincludesubgrid
c 1001 format (' lkset_solve1: lksolving,lkaveragehead,',
c     &        ' lkupdateresistances,lkincludesubgrid',4(l3))
      return
      end subroutine
c
c --------------------------------------------------------------------
c
      subroutine lksetconst (nrow,ncol,rlkdeltax,rlkdeltay,rlkconst)
c
c     Routine generates constants for line-sinks to ensure that the line-sink function
c     does not contribute to the head at the reference point.
c     This is necessary for the Gauss-Seidel solution procedure.
c
      implicit none
      INTEGER nrow,ncol,i,j
      REAL(8) rlkdeltax,rlkdeltay,rlkconst,
     &        rx1,ry1,rx2,ry2,rdx,rdy,
     &        RK0,RH0,RHED0,RB0,RP0
      COMPLEX cz0,cz1,cz2,comls
      DIMENSION rlkdeltax(ncol),rlkdeltay(nrow),rlkconst(nrow,ncol+1)
      INCLUDE 'lkcom.inc'
c
      call GVPAR (RK0,RH0,RHED0,RB0,RP0,CZ0)
c
      rx1=rlkx0
      ry1=rlky0
      do j=1,ncol+1
      if (j.ne.ncol+1) rdx=rlkdeltax(j)
      do i=1,nrow
       cz1=CMPLX(rx1,ry1)
       rdy=rlkdeltay(i)
       ry2=ry1-rdy
       cz2=CMPLX(rx1,ry2)
       rlkconst(i,j)=-comls(cz0,cz1,cz2)
       ry1=ry2
      end do
      ry1=rlky0
      rx1=rx1+rdx
      enddo
c
      end subroutine
c
c ---------------------------------------------------------------------------------------------------
c
      LOGICAL function lflksolving()
c
c     passing lksolving
c
      implicit none
      INCLUDE 'lkcom.inc'
c
      lflksolving=lksolving
c
      return
      end
c
c ---------------------------------------------------------------------------------------------------
c
      LOGICAL function lfleakage()
c
c     passing lkleakage
c
      implicit none
      INCLUDE 'lkcom.inc'
c
      lfleakage=lkleakage
c
      return
      end
c
c ----------------------------------------------------------------------------------------------------
c
      subroutine lksetsubheadlower (nrow,ncol,rlkheadlower,
     &                       ilkresolution,
     &                       ilkpointer,rlksubheadlower,nbuf)
c
c     Routine generates the lower heads at sub-cell centers.

c     For MODFLOW grids smaller than 5*5 the heads at sub-cell centers are set equal to the head at
c     the corresponding MODFLOW cell center.
c
c
      implicit none
      LOGICAL lout
      integer nrow,ncol,nbuf,ilkresolution,i,j,iad,ires,
     &        ilkpointer,istart,ii,jj,ierr,icode
      REAL(8) rlkheadlower,rlksubheadlower,rbuf,rdum,rdx,rdy
      COMPLEX(8) cz0,cz1,cz2,cz3
      DIMENSION rlkheadlower(nrow,ncol),rlksubheadlower(nbuf),
     &          ilkresolution(nrow,ncol),ilkpointer(nrow,ncol)
      ALLOCATABLE rbuf(:,:)
      INCLUDE 'lkcom.inc'
      INCLUDE 'lusys.inc'
      save
c
      if (lknosubgridfile) return
c
c      write (iluer,1011) nrow,ncol,nbuf
c 1011 format (' lksetsubheadlower11: nrow,ncol,nbuf=',3i10)
      if (nrow.ge.5.and.ncol.ge.5) then
        iad=0
        do i=2,nrow-1 ! assuming that the perimeter cells do not contain sub-cells.
        do j=2,ncol-1 ! assuming that the perimeter cells do not contain sub-cells.
c        lout=i.eq.5.and.j.eq.5    ! debugging
        lout=.FALSE.   ! turn debugging off
        ires=ilkresolution(i,j)
        if (ires.gt.1) then
      if (ALLOCATED(rbuf)) deallocate (rbuf)
      ALLOCATE (rbuf(ires,ires),stat=ierr)
      if (ierr.ne.0) then
       write (iluer,2000) ires,ires
        AMESS(1)='Error in routine lkinterpolateheads.'
        AMESS(2)='Failed to allocate rbuf(ires,ires).'
        AMESS(3)='Verify resources and rerun.'
        CALL HALT(3) ! stop program execution for batch version
      end if
c
        ! Note: we are forcing a uniform grid in x- and y-direction    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call lk_cell_construct (1,1,rdx,rdy,cz0,cz1,cz2,cz3)
c
        call lkinterpolateheads(nrow,ncol,i,j,rdx,rdy,rlkheadlower,ires,
     &                               rbuf,icode,iluer)
c
          select case (icode)
c
              case (0)
c                      successful interpolation, now assign values to rlksubheadlower()
              istart=ilkpointer(i,j)
      IF (lout) write (iluer,1002) i,j,ires,istart
 1002 format ('lksetsubheadlower2: i,j,ires,istart ',4I5)
              do ii=1,ires
                do jj=1,ires
                  iad=iad+1
                  rlksubheadlower(iad)=rbuf(ii,jj)
                end do
      if (lout) write (iluer,1003) ii,rbuf(ii,1:3)
 1003 format ('lksetsubheadlower3: ii,rbuf(ii,1:3) ',i3,1x,3(d14.7))
              end do
c
              case (3000)
c                         error report for DBICD3
c                         just assign head at center of MODFLOW cell center to all sub-cell centers
              write (iluer,1000)
              istart=ilkpointer(i,j)
              rdum=rlkheadlower(i,j)
              do ii=1,ires
                do jj=1,ires
                  iad=iad+1
                  rlksubheadlower(iad)=rdum
                end do
              end do
c
              case default
          end select

        end if
        end do
        end do
      else         ! less than 5*5 grid, just assign head at MODFLOW cell center to all sub-cell centers.
        do i=1,nrow
          do j=1,ncol
            ires=ilkresolution(i,j)
            if (ires.gt.1) then
              istart=ilkpointer(i,j)
              rdum=rlkheadlower(i,j)
              do ii=1,ires
                do jj=1,ires
                  iad=istart+(ii-1)*ires+jj-1
                  rlksubheadlower(iad)=rdum
c      write (iluer,1008) i,j,ii,jj,istart,iad,rlksubheadlower(iad)
c 1008 format ('lksetsubheadlower8: i,j,ii,jj,istart,iad,rlksubheadlower'
c     &' ',6i3,1x,d14.7)
                end do
              end do
            end if
          end do
        end do
      end if
c
      return
c
 1000 format (' ***ERROR in lksetsubheadlower: interpolation of lower',
     &' heads failed. Contact technical support.',/,
     &' Heads at sub-cell centers set to corresponding head at MODFLOW',
     &' cell center.')
 2000 format ('***ERROR in LKINTERPOLATEHEADS: failed to allocate rbuf('
     &,i5,',',i5,').')
      end subroutine
c
c ---------------------------------------------------------------------------------------------------
c

      subroutine lkinterpolateheads(nrow,ncol,i,j,rdx,rdy,rlkheadlower,
     &ires,rbuf,icode,iluer)
!---------------------------------
!  Bilinear interpolator for GFMF coupling
!  Limitations:
!  1)  This interpolator cannot be used for cells on the model edge
!      (every cell being interpolated must have 4 adjacent cells)
!  2)  The minimum subgrid resolution for a cell i,j is ires=1
!  3)  When the center points for all cells in the original model are connected,
!      a rectangular grid must be formed for this model to be valid
!  4)  Interpolated values may not fall below the lowest head
!      nor rise above the highest head in the original model
!
!      Coded and tested by Daniel Abrams, November 2010
!
!---------------------------------
!
      implicit none
      LOGICAL lout
      integer nrow,ncol,i,j,ii,jj,ires,icode,iluer
      REAL(8) rlkheadlower,rbuf,rdx,rdy,ii_real,jj_real,ires_real,
     &xdist,ydist
      Dimension rlkheadlower(nrow,ncol)
      Dimension rbuf(ires,ires)

!--------------------------------------------------
!  Check that i and j are not on the boundary of the model , otherwise interpolation cannot proceed
!--------------------------------------------------
c      lout=i.eq.5.and.j.eq.5  ! debugging active for cell i=5, j=5
      lout=.FALSE. ! turn off debugging
      icode=0
      if(i.gt.1.and.i.lt.nrow.and.j.gt.1.and.j.lt.ncol)then
	!------------------------
	!  Cycle through subgrid
	!------------------------
      do ii=1,ires
      	do jj = 1,ires
      		!-------------------------
      		!  Need to do math with integers, change to reals
      		!-------------------------
      	ii_real = REAL(ii)
      	jj_real = REAL(jj)
      	ires_real = REAL(ires)
      		!----------------------------------
      		!  how far is the center of subcell ii,jj from the center of the original cell
      		!  Note, if the subcell falls above the original cell center, y distance is negative
      		!  If the subcell falls to the left of the original cell center, x distance is negative
      		!----------------------------------
      	ydist = (ii_real-0.5)*rdy/ires_real - 0.5*rdy  !
      	xdist = (jj_real-0.5)*rdx/ires_real - 0.5*rdx
      if (lout) write (iluer,1001) ii,jj,ydist,xdist
 1001 format (' lkinterpolateheads1 ',2I5,2(d14.7))
       		!----------------------------------
      		!  if cell has odd number of subcells, center subcell will have same head as original cell
      		!----------------------------------
      	if(xdist.eq.(0.).and.ydist.eq.(0.))then
          rbuf(ii,jj) = rlkheadlower(i,j)
			!-------------------------------------
			!  when the point is the same in x but negative in ydistance, we look at the head above (decrease in i)
			!-------------------------------------
	else if(xdist.eq.(0.).and.ydist.lt.(0.))then
	  rbuf(ii,jj) = rlkheadlower(i,j)+(rlkheadlower(i+1,j)-
     &                  rlkheadlower(i,j))*abs(ydist)/rdy
			!---------------------------------------
			!  when the point is the same in x but positive in ydistance, we look at the head below (increase in i)
			!---------------------------------------
	else if(xdist.eq.(0.).and.ydist.gt.(0.))then
	  rbuf(ii,jj) = rlkheadlower(i,j)+(rlkheadlower(i-1,j)-
     &                  rlkheadlower(i,j))*abs(ydist)/rdy
			!------------------------------------------
			!  when the point is the same in y but positive in xdistance, we look at the head to the right (increase in j)
			!------------------------------------------
	else if(xdist.gt.(0.).and.ydist.eq.(0.))then
	  rbuf(ii,jj) = rlkheadlower(i,j)+(rlkheadlower(i,j+1)-
     &                  rlkheadlower(i,j))*xdist/rdx
			!-------------------------------------------
			!  when the point is the same in y but negative in xdistance, we look at the head to the left (decrease in j)
			!-------------------------------------------
	else if(xdist.lt.(0.).and.ydist.eq.(0.))then
	  rbuf(ii,jj) = rlkheadlower(i,j)+(rlkheadlower(i,j-1)-
     &                  rlkheadlower(i,j))*abs(xdist)/rdx
			!---------------------------------------
			!  Interpolate head in subcell based on 4 cells (if the subcell is in the upper right quadrant)
			!---------------------------------------
	else if(xdist.gt.(0.).and.ydist.lt.(0.))then
	  rbuf(ii,jj) = (abs((rdx-xdist)*(rdy+ydist)*rlkheadlower(i,j))+
     &                  abs((xdist)*(rdy+ydist)*rlkheadlower(i,j+1))+
     &                  abs((xdist)*(ydist)*rlkheadlower(i+1,j+1))+
     &  	        abs((rdx-xdist)*(ydist)*rlkheadlower(i+1,j)))/
     &                  (rdx*rdy)
			!------------------------
			!  verify that calculated head is between 4 heads used to calculate it, if not, assign error as icode = 3000
			!-------------------------
	 if (rbuf(ii,jj).lt.min(rlkheadlower(i,j),rlkheadlower(i+1,j+1),
     &                    rlkheadlower(i,j+1),rlkheadlower(i+1,j)))then
	  	icode = 3000
	  else if (rbuf(ii,jj).gt.max(rlkheadlower(i,j),
     &                    rlkheadlower(i+1,j+1),rlkheadlower(i,j+1),
     &                    rlkheadlower(i+1,j)))then
	  	icode = 3000
	  end if
			!-----------------------------------
			!  Repeat for remaining quadrants
			!-----------------------------------
	else if(xdist.gt.(0.).and.ydist.gt.(0.))then
	  rbuf(ii,jj) = (abs((rdx-xdist)*(rdy-ydist)*rlkheadlower(i,j))+
     &	  abs((xdist)*(rdy-ydist)*rlkheadlower(i,j+1))+
     &    abs((xdist)*(ydist)*rlkheadlower(i-1,j+1))+
     &	  abs((rdx-xdist)*(ydist)*rlkheadlower(i-1,j)))/(rdx*rdy)
			!------------------------
c
	 if (rbuf(ii,jj).lt.min(rlkheadlower(i,j),rlkheadlower(i-1,j+1),
     &                    rlkheadlower(i-1,j),rlkheadlower(i,j+1)))then
	  	icode = 3000
	  else if (rbuf(ii,jj).gt.max(rlkheadlower(i,j),
     &                  rlkheadlower(i-1,j+1),rlkheadlower(i-1,j),
     &                  rlkheadlower(i,j+1)))then
	  	icode = 3000
	  end if
	else if(xdist.lt.(0.).and.ydist.gt.(0.))then
	  rbuf(ii,jj) = (abs((rdx+xdist)*(rdy-ydist)*rlkheadlower(i,j))+
     &	  abs((xdist)*(rdy-ydist)*rlkheadlower(i,j-1))+
     &    abs((xdist)*(ydist)*rlkheadlower(i-1,j-1))+
     &	  abs((rdx+xdist)*(ydist)*rlkheadlower(i-1,j)))/(rdx*rdy)
	  if(rbuf(ii,jj).lt.min(rlkheadlower(i,j),rlkheadlower(i-1,j-1),
     &     rlkheadlower(i,j-1),rlkheadlower(i-1,j)))then
	  	icode = 3000
	  else if (rbuf(ii,jj).gt.max(rlkheadlower(i,j),
     &       rlkheadlower(i-1,j-1),rlkheadlower(i,j-1),
     &       rlkheadlower(i-1,j)))then
	  	icode = 3000
	  end if
	else if(xdist.lt.(0.).and.ydist.lt.(0.))then
	  rbuf(ii,jj) = (abs((rdx+xdist)*(rdy+ydist)*rlkheadlower(i,j))+
     &	  abs((xdist)*(rdy+ydist)*rlkheadlower(i,j-1))+
     &    abs((xdist)*(ydist)*rlkheadlower(i+1,j-1))+
     &	  abs((rdx+xdist)*(ydist)*rlkheadlower(i+1,j)))/(rdx*rdy)
	  if(rbuf(ii,jj).lt.min(rlkheadlower(i,j),rlkheadlower(i+1,j-1),
     &     rlkheadlower(i+1,j),rlkheadlower(i,j-1)))then
	  	icode = 3000
	  else if (rbuf(ii,jj).gt.max(rlkheadlower(i,j),
     &     rlkheadlower(i+1,j-1),rlkheadlower(i+1,j),
     &     rlkheadlower(i,j-1)))then
	  	icode = 3000
	  end if
	end if
	end do
      end do
!------------------------------------------
!  if initial conditions not met (cell not on boundary), assign an error code 3000
!------------------------------------------
      else
	  icode = 3000
      end if
!
      return
!
      end subroutine
c
c ----------------------------------------------------------------------------------------------
c
c    NOTE: THIS ROUTINE MAY CAUSE PROBLEMS IN A NON-SQUARE GRID - TO BE TESTED IF USED AGAIN.
c
c
!      subroutine lkinterpolateheads(nrow,ncol,i,j,rlkheadlower,ires,
!     &                               rbuf,icode)
c
c     For a MODFLOW grid of 5*5 or larger the heads at sub-cell centers will be iterpolated
c     from MODFLOW cell centers in 25 surrounding MODFLOW cells using a "B-spline two-dimensional
c     interpolation." This interpolation is performed by the routines DBICD3 and DBIFD3 of
c     the SSL2 library of Lahey Fortran.
c
c
c     nrow by ncol MODFLOW grid.
c     lower heads at centers stored in rlkheadlower(:,:)
c     seek heads at sub-cell centers for the i,j th MODFLOW cell
c     the subgrid has a resolution of ires by ires
c     interpolated head for sub-cell ii,jj is stored in rbuf(ii,jj)
c
c     MODFLOW grid must be at least 5 by 5 for interpolation to occur, else the head at
c     the MODFLOW cell center is given to all sub-cell centers.
c
c     icode (output) 0 if interpolation successful, else an error code (3000).
c
c
!      implicit none
!      integer nrow,ncol,i,j,iad,isub,ires,npts,mres,
!     &        ix,iy,ii,jj,morder,icode,k,iswx,iswy,
!     &        itemp,is,ie,js,je,jad,ierr
!      LOGICAL l1,l2,linterpolate
!      REAL(8) rlkheadlower,rbuf,rtest,
!     &        rdx,rdxx,rdy,rdyy,rx,ry,resolution,rxx,ryy,
!     &        rxc,ryc,rhed,rcoef,rknots,rscr,rfinterp
!      DIMENSION rlkheadlower(nrow,ncol)
!      parameter (morder=3,mres=10,k=5,iswx=0,iswy=0,npts=5)
!      DIMENSION rcoef(k,k),rknots(8), ! rknots over-dimensioned: calculated dimension = 6
!     &          rscr(30), ! rscr over-dimensioned: calculated = 27
!     &          rxc(npts),ryc(npts),rhed(k,k),rbuf(ires,ires)
!      INCLUDE 'lkcom.inc'
!      INCLUDE 'lusys.inc'
!      save
c
!      if (nrow.ge.5.and.ncol.ge.5) then   ! OK MODFLOW grid large enough to do interpolation
!      is=MAX(1,i-2)
!      js=MAX(1,j-2)
!      ie=MIN(nrow,i+2)
!      je=MIN(ncol,j+2)
!      if (is.eq.1) ie=is+4
!      if (ie.eq.nrow) is=nrow-4
!      if (js.eq.1) je=js+4
!      if (je.eq.ncol) js=ncol-4
!      iad=0
!      rtest=rlkheadlower(i,j)
!      linterpolate=.false.
!      do ii=ie,is,-1         ! setup matrix of heads and associated coordinates
!        iad=iad+1
!        call lk_cell_construct (ii,js,rdx,rdy,cz0,cz1,cz2,cz3)
!        ryc(iad)=0.5*aimag(cz0+cz2)
!        jad=0
!        do jj=js,je
!          jad=jad+1
!          if (ii.eq.ie) then
!            call lk_cell_construct (ie,jj,rdx,rdy,cz0,cz1,cz2,cz3)
!            rxc(jad)=0.5*REAL(cz0+cz2)
!          end if
!          rhed(jad,iad)=rlkheadlower(ii,jj)  ! note rhed(x,y)
!          if (ABS(rhed(jad,iad)-rtest).gt.0.0001) linterpolate=.TRUE. ! interpolate if there are variations in head
!        end do
!      end do
!          if (linterpolate) then  ! found variations in lower head, DO interpolate
c
!          call DBICD3(RXC,npts,RYC,npts,RHED,K,morder,     ! generate interpolation coefficients
!     &                RCOEF,RKNOTS,RSCR,ICODE)
c
!          select case (icode)
!              case (0)
c
!          call lk_cell_construct(i,j,rdx,rdy,cz0,cz1,cz2,cz3)
!          resolution=real(ires)
!          rdxx=rdx/resolution
!          rdyy=rdy/resolution
!          rx=0.5*rdxx
!          ry=0.5*rdyy
!          cz=cz0+CMPLX(rx,ry)
!          do ii=1,ires
!            do jj=1,ires
c
!        ! set up the interpolation process here; calling DBIFD3
!              rxx=REAL(cz)
!              ryy=AIMAG(cz)
!              ix=1
!              iy=1
!              do itemp=2,5
!                if (rxx.ge.rxc(itemp)) then
!                  ix=itemp
!                end if
!                if (ryy.ge.ryc(itemp)) then
!                  iy=itemp
!                end if
!              end do
!              call DBIFD3 (rxc,npts,ryc,npts,morder,rcoef,k,rknots,
!     &                     iswx,rxx,ix,iswy,ryy,iy,rfinterp,rscr,icode)
!              rbuf(ii,jj)=rfinterp
!              cz=cz+CMPLX(rdxx,0.0d0)
!            end do
!            rx=0.5*rdxx
!            ry=ry+rdyy
!            cz=cz0+CMPLX(rx,ry)
!          end do
!         end select
!         else   !            found no variations in lower head, DO NOT interpolate
!         write (iluer,1005)
! 1005 format ('lkinterpolateheads5: ',
!     &        'Found no variation lower heads, hence no interpolation.')
!              do ii=1,ires
!                do jj=1,ires
!                  rbuf(ii,jj)=rtest
!                end do
!              end do
!         endif
!      else             ! grid is less than 5*5 return head at center of cell
!              rfinterp=rlkheadlower(i,j)
!              do ii=1,ires
!                do jj=1,ires
!                  rbuf(ii,jj)=rfinterp
!                end do
!              end do
!      end if
c
c     DEBUGGING
c
c      write (iluer,1009) i,j
c 1009 format ('lkinterpolateheads9: MODFLOW cell i,j 'i3,2x,i3)
c      do ii=ires,1,-1
c      write (iluer,1010) ii,rbuf(ii,1:ires)
c 1010 format ('lkinterpolateheads10: rbuf(',i1,',1:ires=',10(d14.7))
c      end do
!
c
!      return
c
c
!       end subroutine
