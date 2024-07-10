C     Last change:  HMH  18 Jul 2010    6:40 pm
c     This file contains the following routines and functions
c
c     SUBROUTINE SOLUT           driver for groundwater flow solution process
c     initialupdate              multiplies the matrix times the solution vector to get increments in the BC
c     subroutine condition       normalizes the matrix to diagonal elements of 1
c     subroutine adjust          adjusts the known vector based on normalized matrix
c     SUBROUTINE DECOMP          decomposes the matrix for gaussian elimination
c     SUBROUTINE SOLVE           forward elimination and back substitution (use after DECOMP)
c     subroutine Gauss_Seidel    double sweep Gauss-Seidel matrix solver
c     CorrectKnowVector          <not in use>
c     SpecialSolve               <not in use>
c     Sherman_Morrison           modification of solution vector based on line-sinks added or dropped
c     gfread_convergence_file    reads file with GW solution convergence criteria
c     reshapematrix              reorganizes elements in matrix after it has been made square
c
c
c --------------------------------------------------------------------------------
c
      SUBROUTINE SOLUT
     &(DRA,IMXSZE,nstr,DRSCR,DRB,DRFAC,CALPH,CZI,ITYPE,IPIV,NITER)
c
c --------------------------------------------------------------------------------
c
C                                                
C   Groundwater flow solution driver
c
c   if niter=0 and lErrorReport=.false. only an error report will be issued.
c
c
C      
      IMPLICIT NONE
      INTEGER(4) IMXSZE,ITYPE,IPIV,NITER,I,ITER,M,N,N1,J,igs,neffective,
     &           iticks,iticks1,iticks2,iticks3,iticks4,ItickStart,
     &           ItickEnd,ierr,ntemp,mout,nstr,isubiter
      LOGICAL LSOLVEFIRST,LSURFWAT,lkeepsigmadone,lskip,ltimer,
     &        lnotlskip,lfinterface,lnotczi,lnotczi_short,
     &        lflkincludesubgrid,lupdatecheck
      REAL(8) DRA,DRSCR,DRB,DRFAC,DRCOND,RDUM,RFPOT,dlabda,dum,dv,dalpha
      COMPLEX(8) CZ,CALPH,CZI,cz1,cz2,cz3,czi_short,czi_long,cztest,
     &           czreset
      CHARACTER(256) amessage
      ALLOCATABLE lskip(:),czi_short(:),czi_long(:)
      INCLUDE 'MAIN.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'MATCH.INC'
      INCLUDE 'TRACOM.INC'
      DIMENSION DRA(IMXSZE,nstr),DRSCR(imxsze),DRB(imxsze),
     &          DRFAC(4,imxsze), ! drscr may be larger
     &          CALPH(imxsze),CZI(imxsze),ITYPE(imxsze),IPIV(imxsze)
      DATA lnotczi_short,lnotczi,lnotlskip,lkeepsigmadone
     &        /.true.,.TRUE.,.TRUE.,.false./
      DATA czreset /(1.0d21,1.0d21)/
      SAVE  ! important to keep data in arrays when reentering
c
c
      ltimer=.false.      ! control execution time reporting for development purposes
      lupdatecheck=.false.
c
c --- check input data and allocate some scratch arrays if this is the first execution
c
      if (imxsze.ne.Number_equations) then
          write (iluer,5001) imxsze,Number_equations
 5001 format (' ***Error in Solut: imxsze=',i6,'Number_equations=',
     &i6,/,' these must be equal. Program aborted.')
          AMESS(1)='Error in SOLUT routine.'
          AMESS(2)='Declared dimension of matrix array not equal to'
          AMESS(3)='the number of equations, see error.log file.'
          CALL HALT(3)   ! stop program execution for batch version
      end if
c
      if (nstr.ne.Number_strengths) then
          write (iluer,5002) nstr,Number_strengths
 5002 format (' ***Error in Solut: nstr=',i6,'Number_strengths=',
     &i6,/,' these must be equal. Program aborted.')
          AMESS(1)='Error in SOLUT routine.'
          AMESS(2)='Declared dimension of matrix array not equal to'
          AMESS(3)='the number of strengths, see error.log file.'
          CALL HALT(3)   ! stop program execution for batch version
      end if
c
      if (lnotczi) then
        ALLOCATE (czi_short(imxsze),czi_long(imxsze),stat=ierr)   ! allocate scratch arrays
        if (ierr.ne.0) then
          call iostat_msg (ierr,amessage)
          write (ilume,1111) amessage
 1111     format (a132)
          deallocate (czi_short,czi_long)
          write (ilume,8000) imxsze
          AMESS(1)='Error in SOLUT routine.'
          AMESS(2)='Failed to allocate a scratch array.'
          AMESS(3)='Execution has been aborted.'
          CALL HALT(3)   ! stop program execution for batch version
        end if
        lnotczi=.false.
      end if
c
      if (lfinterface().and.lDirectFromDisk) then
        lDirectFromDisk=.FALSE. ! force this in the presence of interface flow.
        write (iluer,1200)
      endif
      call setldbmatrix_true()     ! enable control points on vertices (no point shifts in DBCOM)
      if (lnotlskip) then
        ALLOCATE (lskip(imxsze),stat=ierr)   ! allocate an extra scratch array, when solving over disc
        if (ierr.ne.0) then
          call iostat_msg (ierr,amessage)
          write (ilume,1112) amessage
 1112     format (a132)
          deallocate (lskip)
          write (ilume,8001) imxsze
          AMESS(1)='Error in SOLUT routine.'
          AMESS(2)='Failed to allocate a scratch array. Try to solve'
          AMESS(3)='without using decomposed matrices stored to disk.'
          CALL HALT(3)   ! stop program execution for batch version
        end if
        lnotlskip=.false.
      end if
c
c ---- if called with NITER=0 just do a check on meeting boundary conditions
c
      if (niter.eq.0.and..not.lErrorReport) then
        if (.not.lucon) WRITE (ILUME,1500)
        write (*,1500)
        RERMAX=0.0
        write (ilume,9070)
        write (*,9070)
        CALL GVERROR (RERMAX)   ! Note: make sure that the known vector arrays (for all elements) are current
        CALL WLERROR (RERMAX)
        CALL LSERROR (RERMAX)
c        CALL PDERROR (RERMAX)
        CALL DBERROR (RERMAX)
        CALL W3ERROR (RERMAX)
        CALL LKERROR (RERMAX)
c        CALL DIERROR (RERMAX)
        write (ilume,7000)
        call setldbmatrix_false()   ! enable point shift logic in doublet routines (no calculations at vertices)
        return
      end if
c
c ---- NITER>0 get ready to generate a solution ---------------------
c
      call timer(iticks1)
      DO I=1,IMXSZE     ! INITIALIZE ARRAYS
      DRSCR(I)=0.0D00
      if (nsol.eq.0) DRB(I)=0.0D00  ! else it may contain line-sink or leakage element strenghts differences
      DRFAC(1,I)=0.0D00
      DRFAC(2,I)=0.0D00
      DRFAC(3,I)=0.0D00
      DRFAC(4,I)=0.0D00
      CALPH(I)=(0.0D00,0.0D00)
      CZI(I)=(0.0D00,0.0D00)
      ITYPE(I)=0
      if (.not.lDirectFromDisk) IPIV(I)=0
      END DO
      dra(1:Number_equations,1:Number_strengths)=0.0d0
      if (loadsol.and..not.lkeepsigmadone) then ! we entered SOLUT for the first time after a load *.sol
       drb(1:nstr)=0.0d0 ! initialize so only changes line-sink strengths in LSSTREAM and LSCZC and
                     ! changes in leakage rates due to import of *.vlk file will be in DRB
       j=0            ! array counter
       call lskeepsigma (drb,j,imxsze)  ! new line-sinks strengths in LSSTREAM and LSCZC will lead
c                                        to strength differences to be used to update all known vector
c     <normally called at end>           arrays, see statement at end of LSCZC and calls of update routines
c                                        after the known vectors have been updated.
       call dbkfac()
       call dbkeep(drb,j,nstr)        ! dummy routine to update array counter, no changes in doublet strengths
       call wlkeep(drb,j,nstr)        ! dummy routine to update array counter, no changes in well discharges
       call w3keep(drb,j,nstr)        ! dummy routine to update array counter, no changes in well discharges
       call lkreadnewleakages(drb,j,nstr)  ! get new leakages from *.vlk and store difference in drb
c Note: this is one step, while for line-sink strength we first store old values and in LSCZC store differences.
       call gvkeep(drb,j,nstr)        ! dummy routine to update array counter, no changes in given constant
       IF (J.NE.nstr) then
        write (ilume,2950) nstr,J ! j is counter of known vector elements that have been updated
        AMESS(1)='Keep strengths array update out of sync.'
        AMESS(2)='Contact technical support.'
        CALL HALT(2) ! stop program execution for batch version
       endif
       lkeepsigmadone=.true.
      end if
      call timer (iticks2)
      iticks=iticks2-iticks1
      if (ltimer) write (ilume,9002) iticks
 9002 format(' Array initialization execution time=               ',i10,
     &' E-2 seconds.')
c
      LSOLVEFIRST=NSOL.EQ.0
      IF (LSOLVEFIRST) THEN
        CALL DBCHECKBOTTOM ()  ! ensure no aquifer bottoms are set above the aquifer top
        CALL LSPREP()          ! generate line-sink specific constants for the potential function
        call SetInitialConstant() ! force potential at reference point as specified when all strenght parameters are still zero
      ENDIF
c
c     ------------------------------------------
c     start groundwater flow solution iterations
c     ------------------------------------------
c
  5   DO  30 ITER=1,NITER
      call timer(ItickStart)
      M=0
      N=0
      lskip(1:imxsze)=.false.
      call dbkfac ! update inhomogeneity domain correction factors (before nsol>0)
c
c     ----------------------------------------------------------------------
c     construction of collocation point vector (and update line-sink status)
c     ----------------------------------------------------------------------
c
      if (LerrorReport) write (ilume,9010)
      if (LerrorReport) write (*,9010)
 9010 format (' generating collocation point vector')
      call timer (iticks1)
      N1=N
c      if (lsolvefirst) then
c        ntemp=1
c      end if
c      write (ilume,1001) (i,drb(i),i=1,ntemp)
c 1001 format ('solut1: drb(',i3,')=',d14.7)
      CALL LSCZC (CZI,M,N,DRB,DRFAC,CALPH,ITYPE,  ! must be first in the list, see comment before "lskeepsigma" call
     &            lsurfwat,lDirectFromDisk,lskip,ltimer)
      !write (iluer,1001) M,N
 1001 format (' solut1: LS:  M,N ',2(i5))
      CALL DBCZC (CZI,M,N,DRFAC,CALPH,ITYPE)
      !write (iluer,1002) M,N
 1002 format (' solut2: DB: M,N ',2(i5))
      CALL WLCZC (CZI,M,N,DRFAC,CALPH,ITYPE)
      !write (iluer,1003) M,N
 1003 format (' solut3: WL: M,N ',2(i5))
      CALL W3CZC (CZI,M,N,DRFAC,CALPH,ITYPE)
      !write (iluer,1004) M,N
 1004 format (' solut4: W3: M,N ',2(i5))
      CALL LKCZC (CZI,M,N,Number_equations,
     &            DRFAC,CALPH,ITYPE,lDirectFromDisk,ltimer)
      !write (iluer,1005) M,N
 1005 format (' solut5: LK: M,N ',2(i5))
      CALL GVCZI (CZI,M,N,DRFAC,ITYPE)
      !write (iluer,1006) M,N
 1006 format (' solut6: GV: M,N ',2(i5))
c
      czi_long(1:m)=czi(1:m) ! keep for update routines
      if (loadsol.and.lnotczi_short) then
        call lkcombineequations (dra,czi,drfac,itype,lskip,m,mout,n) ! just to shrink czi
        czi_short(1:n)=czi(1:n) ! store czi for square matrix
        lnotczi_short=.false.
      end if
c
c     keep one control point to follow the potential calculations.
c
      if (lsolvefirst) then
      cztest=czi(1)
      end if
c
c      Provide a check of M=Number_equations and N=Number_strengths !!!!! TO BE DONE
c
      if (m.ne.Number_equations.or.n.ne.Number_strengths) then
          write (iluer,1900) M,Number_equations,N,Number_strengths
          AMESS(1)='Matrix is sized incorrectly.'
          AMESS(2)='See ERROR.LOG file for details.'
          AMESS(3)='Contact technical support.'
          CALL HALT(3) ! stop program execution for batch version
      end if
c
      call timer (iticks2)
      iticks=iticks2-iticks1
      if (ltimer) write (ilume,9011) iticks
 9011 format(' Collocation point vector execution time=           ',i10,
     &' E-2 seconds.')
c
c                        --------------------------------------
c                        update known vector arrays of elements
c                        --------------------------------------
c
c Note: this must be done because line-sink strengths may have been changed in LSCZC
c       such changes will be stored in drb
c
      call timer (iticks1)
      if (nsol.gt.0) then  ! load initial matrix and use to store known vector array updates in DRSCR
        call LoadMatrix (dra,m,n,lErrorReport,ltimer)
        call initialupdate (dra,drb,drscr,m,n,lErrorReport)
      endif
c           Note: drscr(M) will be used to update M known vector elements!!
c                 we are not adding sub-cell contributions at this point since
c                 they should already be included
      j=0
c      write (iluer,2345)
c 2345 format (' solut2345: we are now updating after collocation vector'
c     &,/,'has been (re)created and line-sink strengths may have been '
c     & 'altered.')
      call lsupdate (drscr,j,m,n,.FALSE.,0)
      call dbupdate (drscr,j,m,n,czi_long,calph,.FALSE.,0)
      call wlupdate (drscr,j,m,n,czi_long,.FALSE.,0)
      call w3update (DRSCR,J,M,N,czi_long,CALPH,ITYPE,.FALSE.,0)
      call lkupdate (drscr,j,m,n,czi_long,.FALSE.,0)
      call gvupdate (drscr,j,m,n)
        IF (J.NE.m) then
          write (ilume,2900) N,J  ! j is counter of known vector elements that are updated.
          AMESS(1)='Known vector array update out of sync.'
          AMESS(2)='Contact technical support.'
          CALL HALT(2) ! stop program execution for batch version
        endif
c
      if (lupdatecheck) then   ! if true, check that the BC are calculated correctly
c
      write (ilume,2905)
c
      j=0
      call lsupdate_check(j,m,n)
      call dbupdate_check(j,m,n,czi_long,calph,itype)
      call wlupdate_check(j,m,n,czi_long)
      call w3update_check(J,m,n,czi_long,ITYPE)
      call lkupdate_check(j,m,n,czi_long)
      call gvupdate_check(j,m,n)
        IF (J.NE.m) then
          write (ilume,2901) N,J  ! j is counter of known vector elements that are checked.
          AMESS(1)='Known vector array check out of sync.'
          AMESS(2)='Contact technical support.'
          CALL HALT(2) ! stop program execution for batch version
        endif
c
      end if ! End of debugging of known vector update routines
c
      call timer(iticks2)
      iticks=iticks2-iticks1
      if (ltimer) write (ilume,9044) iticks
 9044 format(' Known vector arrays update execution time=         ',i10,
     &' E-2 seconds.')
c
      NSOL=NSOL+1
c      if (lDirectFromDisk) then
c      calculate the number of actual unknown strength parameters, "neffective", since
c      some line-sinks may have become "known" and removed from the matrix
c      solution process. This decision is made in LSCZC and recorded in array "lskip"
         neffective=0
         do i=1,n                                     !!!!! LOOK AT THIS, HOW DO WE USE "M"
           if (.not.lskip(i)) neffective=neffective+1
         end do
         if (.not.lucon) WRITE (ILUME,1100) NSOL,neffective
         write (*,1000) nsol,neffective
c      else
c         if (.not.lucon) WRITE (ILUME,1100) NSOL,N
c         write (*,1000) nsol,n
c      end if
      IF (NSOL.GT.1) LSOLVEFIRST=.FALSE.
c       call matout (1,DRA,DRB,CZI,CALPH,DRFAC,M,N,ITYPE)
      IF (m.GT.IMXSZE) THEN
        WRITE (ILUER,4000) IMXSZE
        AMESS(1)='Too many equations for this version of GFLOW1.EXE'
        AMESS(2)='Reduce the size of your problem and try again.'
        CALL HALT(2) ! stop program execution for batch version
        NSOL=-1
        RETURN
      ENDIF
c
c     ----------------------
c     construction of matrix
c     ----------------------
c
      IF (LSOLVEFIRST.or.lWrongFile) THEN ! matrix construction will only be done once
        call timer (iticks1)
        DRA(1:M,1:N)=0.0D0 ! initialize matrix
        call timer (iticks2)
        iticks=iticks2-iticks1
        if (ltimer) write (ilume,9012) iticks
 9012 format(' Initializing of matrix coefficients time=          ',i10,
     &  ' E-2 seconds.')
        J=0
        if (LerrorReport) write (ilume,9020)
        if (LerrorReport) write (*,9020)
 9020   format (' generating matrix coefficients')
        call timer (iticks1)
        CALL LSGENMAT (DRA,CZI,M,N,J,DRFAC,CALPH,ITYPE,
     &                 lDirectFromDisk,ltimer)
        CALL DBGENMAT (DRA,CZI,M,N,J,DRFAC,CALPH,ITYPE)
c        call matrixrowout  (dra,m,500,'matrix row 500 db')     !!!!!!!!!! debugging
        CALL WLGENMAT (DRA,CZI,M,N,J,DRFAC,CALPH,ITYPE)
c        call matrixrowout  (dra,m,500,'matrix row 500 wl')     !!!!!!!!!! debugging
        CALL W3GENMAT (DRA,CZI,M,N,J,CALPH,ITYPE)
c        call matrixrowout  (dra,m,500,'matrix row 500 w3')     !!!!!!!!!! debugging
        CALL LKGENMAT (DRA,CZI,M,N,J,CALPH,ITYPE)
c        call matrixrowout  (dra,m,500,'matrix row 500 lk')     !!!!!!!!!! debugging
        CALL GVGENMAT (DRA,M,N,J,DRFAC,ITYPE)
C        CALL PDMAT (DRA,CZI,N,J,DRFAC,CALPH,ITYPE)
C        CALL DIMAT (DRA,CZI,N,J,DRFAC,CALPH,ITYPE)
c        CALL LKMAT (DRA,CZI,N,J,DRFAC,CALPH,ITYPE)
        call timer (iticks2)
        iticks=iticks2-iticks1
        if (ltimer) write (ilume,9013) iticks
 9013 format(' Matrix generation execution time=                  ',i10,
     &  ' E-2 seconds.')
c         call matout (2,DRA,DRB,CZI,CALPH,DRFAC,M,N,ITYPE)
C
c        call matrixrowout  (dra,m,500,'matrix row 500  ')
c        call matrixrowout (dra,m,1145,'matrix row 1145 ')
c
      if (m.ne.Number_equations.or.n.ne.Number_strengths) then
          write (ilume,3900) M,Number_equations,N,Number_strengths
          AMESS(1)='Matrix is sized incorrectly.'
          AMESS(2)='See ERROR.LOG file for details.'
          AMESS(3)='Contact technical support.'
          CALL HALT(3) ! stop program execution for batch version
      end if
c
        call WriteMatrix (DRA,M,N,lErrorReport,ltimer) ! store coefficient matrix  !!! not now, Debugging!!!!
      ENDIF
      IF (.NOT.lDirectFromDisk.or.(LSOLVEFIRST.or.lWrongFile)) THEN
c
c                         ---------------------------
c                         correct matrix coefficients
c                         ---------------------------
c
c     (resistance line-sinks, galleries, Q-spec. 3D well, etc.)
c
        if (nsol.gt.1) call LoadMatrix (dra,m,n,lErrorReport,ltimer)
        j=0
        call lsmatcorrect (dra,czi,drb,drscr,m,n,j,drfac,calph,itype,
     &                     lDirectFromDisk,ltimer)
        CALL DBMATcorrect (DRA,CZI,M,N,J,DRFAC,CALPH,ITYPE)
        CALL WLMATcorrect (DRA,CZI,M,N,J,DRFAC,CALPH,ITYPE)
        CALL W3MATcorrect (DRA,CZI,M,N,J,CALPH,ITYPE)
c       write (iluer,1234)
 1234 format (' solut1234: entering lkmatcorrect.')
        CALL LKMATcorrect (DRA,CZI,M,N,J,CALPH,ITYPE)
c       write (iluer,1235)
 1235 format (' solut1235: leaving lkmatcorrect.')
        CALL gvmatcorrect (DRA,M,N,J,DRFAC,ITYPE)
c
c                    remove equations when strength became prescribed (in lsczc or genstreamflow)
c
      if (.NOT.(lsolvefirst.and.lDirectFromDisk))
     &call ls_remove_equations (dra,m,n,ltimer) ! do not remove equations first time around when
c                                               the matrix is still to be decomposed.
c
c         call matout (3,DRA,DRB,CZI,CALPH,DRFAC,M,N,ITYPE)
c
c
c
        !write (iluer,1236)
 1236 format (' solut1236: entering lkcombineequations.')
                                                                   ! WARNING --> DEBUGGING
c
        call lkcombineequations (dra,czi,drfac,itype,lskip,m,mout,n) ! combine equations for subcells into one per MODFLOW cell
c
        czi_short(1:n)=czi(1:n) ! store czi for square matrix
        !write (iluer,1237)                                     ! from here we have n equations
 1237 format (' solut1237: leaving lkcombineequations.')
c
c
c
c         call matout (4,DRA,DRB,CZI,CALPH,DRFAC,M,N,ITYPE)  ! still stored as DRA(m,n)
        IF (mout.NE.N) then
          write (ilume,3910) mout,N  ! this should result in a square matrix
          AMESS(1)='Matrix is not square.'
          AMESS(2)='Contact technical support.'
          CALL HALT(2) ! stop program execution for batch version
        endif
c
c        call comparematrices (dra,m,n,'Before reshape  ')
c        call matrixrowout  (dra,m,500,'matrix row 500  ')
c        call matrixrowout (dra,m,1145,'matrix row 1145 ')
c
c        amess(1)='compared matrices, now stop.'
c        call halt (1)                                   !!!!!!!!stop here for now.
c
        !write (iluer,1238)
 1238 format (' solut1238: entering reshapematrix.')
        call reshapematrix (dra,m,n)                    ! from here the matrix is stored as DRA(n,n)
        !write (iluer,1239)
 1239 format (' solut1239: leaving reshapematrix.')
c
c       From this point on the correct first dimension of DRA is N!!
c
c         call matout (5,DRA,DRB,CZI,CALPH,DRFAC,n,N,ITYPE)
c
c        call WriteMatrix (DRA,n,n,lErrorReport,ltimer) ! store corrected matrix this time !!! debugging
c        call comparematrices (dra,n,n,'After reshape   ')
c        call matrixrowout  (dra,n,500,'matrix row 500  ')
c        call matrixrowout (dra,n,1145,'matrix row 1145 ')
c
c        amess(1)='compared matrices, now stop.'
c        call halt (1)                                   !!!!!!!!stop her for now.
c
c
      ENDIF ! end of matrix construction
c                      ----------------------------
c                      construction of known vector
c                      ----------------------------
  20  CONTINUE
      if (LerrorReport) write (ilume,9030)
      if (LerrorReport) write (*,9030)
 9030 format (' generating known vector')
      call timer (iticks1)
      drb(1:n)=0.0d0 ! initialize known vector
      J=0
      CALL LSKNO (DRB,J,czi_short,ITYPE,lDirectFromDisk)
      CALL DBKNO (DRB,N,J,czi_short,CALPH,ITYPE)
      CALL WLKNO (DRB,J,czi_short)                                     ! remember: DRA(n,n)
      CALL W3KNO (DRB,J,czi_short,CALPH,ITYPE)
      CALL LKKNO (DRB,J,czi_short,N)
      CALL GVKNO (DRB,J)
C      CALL PDKNO (DRB,J,CZI)
C      CALL DIKNO (DRB,J,CZI)
C      CALL LKKNO (DRB,J,CZI)
c
      call timer (iticks2)
      iticks=iticks2-iticks1
      if (ltimer) write (ilume,9031) iticks
 9031 format(' Knownvector execution time=                        ',i10,
     &       ' E-2 seconds.')
c        call matout (6,DRA,DRB,CZI,CALPH,DRFAC,n,N,ITYPE)
      if (n.gt.10) then
      end if
      IF (J.NE.N) then
        write (ilume,4100) N,J  ! known vector elements must equal unknowns
        AMESS(1)='Known vector does not have proper length.'
        AMESS(2)='Contact technical support.'
        CALL HALT(2)  ! stop program execution for batch version
      endif
c
      if ((lsolvefirst.or.lWrongFile).or..NOT.lDirectFromDisk) then
c                                    --------------------
c                                    matrix decomposition
c                                    --------------------
        if (LerrorReport) write (ilume,9040)
        if (LerrorReport) write (*,9040)
 9040   format (' decomposing the matrix')
        call timer (iticks1)
        CALL DECOMP (N,DRA,DRCOND,IPIV,DRSCR,DRFAC)
        call timer (iticks2)
        iticks=iticks2-iticks1
        if (ltimer) write (ilume,9041) iticks
 9041 format(' DECOMP execution time=                             ',i10,
     &         ' E-2 seconds.')
        IF (.NOT.LUCON.and.lErrorReport) WRITE (ILUME,6000) DRCOND
        if (lErrorReport) WRITE (*,6000) DRCOND
        IF (DRCOND.GE.1.D30) THEN
          WRITE (ILUER,4500)
          AMESS(1)='Matrix is singular.'
          AMESS(2)='Contact technical support.'
          CALL HALT(2)  ! stop program execution for batch version
        ENDIF
        if (lDirectFromDisk) then
           call WriteDecompMatrix(dra,ipiv,N,lErrorReport,ltimer) ! note: ipiv must be preserved
           lWrongFile=.FALSE. ! new files have been written
        endif
c                               --------------------------------
c                               solve after matrix decomposition
c                               --------------------------------
        if (LerrorReport) write (ilume,9050)
        if (LerrorReport) write (*,9050)
 9050  format (' forward elimination and back substitution (solve)')
        call timer (iticks1)
        if (.not.lDirectFromDisk) then   ! could be on first iteration
           CALL SOLVE (N,DRA,DRB,IPIV) ! forward elimination and back substitution
        else  ! this occurs when both "lsolvefirst" and "ldirectfromdisk" are true
c
c        Replaced Solve by Sherman-Morrison procedure in case there are already equations to be
c        removed on during the first iteration, for instance when "drains" are included.
c
          do i=1,n              ! set LSKIP and adjust known vector DRB
           if (lskip(i)) drb(i)=0.0d0
          end do
c
c     create solution based on stored matrix and skipped equations
c
          call Sherman_Morrison (n,dra,drb,ipiv,lskip,
     &                       lErrorReport,ltimer)
c
          call timer (iticks2)
          iticks=iticks2-iticks1
          if (ltimer) write (ilume,9051) iticks
 9051 format(' SOLVE execution time=                              ',i10,
     &          ' E-2 seconds.')
        endif
c         call vectorout (drb,n,'solution vector2')
      else     ! solve using decomposed matrix from disk
c
c                               -----------------------------------------
        if (lDirectFromDisk) then ! solve using decomposed matrix from disk
c                               -----------------------------------------
c
c     Note: SOLVE is called inside the Sherman-Morrison routine.
c           Both the decomposed and the original matrixes are also loaded in the Sherman_Morrison routine.
c     WARNING: m=m_original here and only "n" rows of the original coefficient matrix is loaded
c              in the Sherman-Morrison routine. OK as long as LK-routines are last (save for GV-routines).
c
          call timer(iticks1)
          do i=1,n              ! set LSKIP and adjust known vector DRB
           if (lskip(i))  drb(i)=0.0d0
c         write (ilume,1002) i,lskip(i)                                  !   DEBUGGING!!!!
c 1002 format ('SOLUT2: i,lskip(i) ',i3,1x,l3)
          end do
c
c     create solution based on stored matrix and skipped equations
c
          call Sherman_Morrison (n,dra,drb,ipiv,lskip,
     &                       lErrorReport,ltimer)
          call timer(iticks2)
          iticks=iticks2-iticks1
          if (ltimer) write (ilume,9045) iticks
 9045 format(' DirectFromDisk (Sherman_Morrison) execution time=  ',i10,
     &          ' E-2 seconds.')
        endif
      endif
c       call matout (7,DRA,DRB,CZI,CALPH,DRFAC,n,N,ITYPE)                ! still DRA(n,n)
c
c                        -----------------------------------
c                        substitution of the solution vector
c                        -----------------------------------
      if (LerrorReport) write (ilume,9060)
      if (LerrorReport) write (*,9060)
 9060 format (' substituting the solution vector')
      call timer (iticks1)
      J=0
      CALL LSSUB (DRB,J)
      CALL DBSUB (DRB,J)
      CALL WLSUB (DRB,J)
      CALL W3SUB (DRB,J)
      CALL LKSUB (DRB,J)
      CALL GVSUB (DRB,J)
C      CALL PDSUB (DRB,J)
C      CALL DISUB (DRB,J)
C      CALL LKSUB (DRB,J)
      call timer (iticks2)
      iticks=iticks2-iticks1
      if (ltimer) write (ilume,9061) iticks
 9061 format(' Substitution of solution vector execution time=    ',i10,
     &' E-2 seconds.')
c       call matout (8,DRA,DRB,CZI,CALPH,DRFAC,n,N,ITYPE)                    ! still DRA(n,n)
      IF (J.NE.N) then
        write (ilume,4200) J,N ! lenght of solution vector must equal unknowns
        AMESS(1)='Error in substitution of solution vector.'
        AMESS(2)='Contact technical support.'
        CALL HALT(2)   ! stop program execution for batch version
      endif
c                --------------------------------------------------------------
c                update known vector arrays of elements for use in error report
c                 -------------------------------------------------------------
c
c               WARNING: Make sure DRB is not changed in SUB routines or before the calls below!!
c
      call timer (iticks1)
      call LoadMatrix (dra,m,n,lErrorReport,ltimer)  !!! NOW we have DRA(m,n) M X N matrix
      call initialupdate (dra,drb,drscr,m,n,lErrorReport) ! use initial matrix to generate modifications of known vector arrays
c
c           Note: drscr(M) will be used to update M known vector elements!!
c                 we are removing sub-cell contributions at this point (multiplyer = -1).
      j=0
c      write (iluer,2346)
c 2346 format (' solut2346: we are now updating after solution vector'
c     &,/,'has been (re)created and before error report ')
      call lsupdate (drscr,j,m,n,.TRUE.,-1)
      call dbupdate (drscr,j,m,n,czi_long,calph,.TRUE.,-1)
      call wlupdate (drscr,j,m,n,czi_long,.TRUE.,-1)
      call w3update (DRSCR,J,M,N,czi_long,CALPH,ITYPE,.TRUE.,-1)
      call lkupdate (drscr,j,m,n,czi_long,.TRUE.,-1)
      call gvupdate (drscr,j,m,n)
        IF (J.NE.m) then
          write (ilume,2900) N,J ! j is counter of known vector elements that have been updated
          AMESS(1)='Known vector array update out of sync.'
          AMESS(2)='Contact technical support.'
          CALL HALT(2) ! stop program execution for batch version
        endif
      call timer(iticks2)
      iticks=iticks2-iticks1
      if (ltimer) write (ilume,9044) iticks
c
      if (lflkincludesubgrid()) then
c
       call timer(iticks1)
c
c              -------------------------------------------
c              solve sub-cells separately
c              -------------------------------------------
c
c       First initialize DRB so that only updated MODFLOW cell leakages will affect the BC
c
        drb(1:n)=0.0d0
        isubiter=4
      write (ilume,7500) isubiter
      write (*,7500) isubiter
       call lksubgridsolve (isubiter,dra,drb,drscr,n,m,imxsze,nstr) ! to test number of iterations=1
c
c      write (iluer,4567) n,drb(1:n)
c 4567 format (' Solut4567: n, drb(1:n) ',i4,2x,10(d14.7))
      call timer(iticks2)
      iticks=iticks2-iticks1
      if (ltimer) write (ilume,9088) iticks
 9088 format(' Solution process of subgrids    execution time=    ',i10,
     &' E-2 seconds.')
c
c
c     At this point include the (new) sub-cell leakages, multiplier = +1  (only reason for these calls)
c
        drscr(1:m)=0.0d0  ! no other potential contributions.
        j=0
        call lsupdate (drscr,j,m,n,.TRUE.,1)
        call dbupdate (drscr,j,m,n,czi_long,calph,.TRUE.,1)
        call wlupdate (drscr,j,m,n,czi_long,.TRUE.,1)
        call w3update (drscr,J,m,n,czi_long,CALPH,ITYPE,.TRUE.,1)
        call lkupdate (drscr,j,m,n,czi_long,.TRUE.,1)
        call gvupdate (drscr,j,m,n)
        IF (J.NE.m) then
          write (ilume,2900) N,J ! j is counter of known vector elements that have been updated
          AMESS(1)='Known vector array update out of sync.'
          AMESS(2)='Contact technical support.'
          CALL HALT(2) ! stop program execution for batch version
        endif
      endif
c
c     Store the current line-sink strengths RLSSIG in DRB for use in LSCZC
c
         drb(1:m)=0.0d0 ! initialize so only changes line-sink strengths in LSSTREAM and LSCZC and
         j=0   ! this implies that LSCZC is the first routine in the list!!
         call lskeepsigma (drb,j,m)  ! new line-sinks strengths in LSSTREAM and LSCZC will lead
c                                   to strength differences to be used to update all known vector
c                                   arrays, see statement at end of LSCZC and calls of update routines
c                                   after the known vectors have been updated.
c Note: At this point there is no need for an lkkeepstrenghts call, we will read the strengths only once
c       for each solver call.
c
c             call matout (9,DRA,DRB,CZI,CALPH,DRFAC,M,N,ITYPE)    ! matrix is DRA(m,n)
c      write (ilume,1003) (i,drb(i),i=1,m)
c 1003 format ('solut3: drb(',i3,')=',d14.7)
c
      if (lupdatecheck) then ! if true check if boundary conditions are calculated correctly
c
      write (ilume,2910)
c
      j=0
      call lsupdate_check(j,m,n)
      call dbupdate_check(j,m,n,czi,calph,itype)
      call wlupdate_check(j,m,n,czi)
      call w3update_check(J,m,n,CZI,ITYPE)
      call lkupdate_check(j,m,n,czi)
      call gvupdate_check(j,m,n)
        IF (J.NE.m) then
          write (ilume,2901) N,J  ! j is counter of known vector elements that are checked.
          AMESS(1)='Known vector array check out of sync.'
          AMESS(2)='Contact technical support.'
          CALL HALT(2) ! stop program execution for batch version
        endif
      call timer(iticks2)
      iticks=iticks2-iticks1
      if (ltimer) write (ilume,9044) iticks
c
      end if ! End of debugging of known vector update routines
c
c
c                                ---------------------------
c                                cleanup and error reporting
c                                ---------------------------
      call timer(iticks1)
      CALL DBPREP (czreset) ! reset the "previous point" memory in DBPREP
      RDUM=RFPOT(czreset)   ! reset the "previous point" memory in RFPOT
      call timer(iticks2)
      iticks=iticks2-iticks1
      if (ltimer) write (ilume,9082) iticks
 9082 format(' cleanup                             execution time=',i10,
     &' E-2 seconds.')
      lquit=.false.
      if (lErrorReport) then ! do not generate report if fastsolve requested
       if (.not.lucon) WRITE (ILUME,1500)
       write (*,1500)
       RERMAX=0.0
       write (ilume,9070)
       write (*,9070)
 9070  format (' checking boundary conditions')
       call timer (iticks1)
       lquit=.TRUE. ! anticipate that GW solution may have converged.
c                      lquit may be set FALSE in any of the XXERROR routines.
c       write (ilume,9999) nsol,minimum_iterations,lquit
       CALL GVERROR (RERMAX)
c       write (ilume,9999) nsol,minimum_iterations,lquit
       CALL WLERROR (RERMAX)
c       write (ilume,9999) nsol,minimum_iterations,lquit
       CALL LSERROR (RERMAX)
c       write (ilume,9999) nsol,minimum_iterations,lquit
c       CALL PDERROR (RERMAX)
       CALL DBERROR (RERMAX)
c       write (ilume,9999) nsol,minimum_iterations,lquit
       CALL W3ERROR (RERMAX)
c       CALL DIERROR (RERMAX)
       CALL LKERROR (RERMAX)
c       write (ilume,9999) nsol,minimum_iterations,lquit
       if (nsol.lt.minimum_iterations) lquit=.false. ! don't abort yet
c       write (ilume,9999) nsol,minimum_iterations,lquit
c 9999  format (' solut: nsol,minimum_iterations,lquit ',i3,i3,l3)
       call timer (iticks2)
       iticks=iticks2-iticks1
       if (ltimer) write (ilume,9071) iticks
 9071  format(' Check BC execution time=                          ',i10,
     & ' E-2 seconds.')
      endif
      write (ilume,7000)
      IF (lquit) GOTO 31
      call timer(ItickEnd)
      iticks=ItickEnd-ItickStart
      if (ltimer) write (ilume,9081) iticks
 9081 format(' Iteration total execution time=                    ',i10,
     &' E-2 seconds.',/)
  30  CONTINUE
  31  continue
      IF (LSOLVEFIRST.AND.LSURFWAT) THEN  ! redo solution to check for percolating surface waters
        if (.not.lucon) WRITE (ILUME,2000)
        write (*,2000)
c       NSOL=NSOL-1   ! do not reset because it fools the logic in lsczc in not updating heads
        LSOLVEFIRST=.FALSE.
        GOTO 5
      ENDIF
      call setldbmatrix_false() ! enable point shift in line doublet routines
      RETURN
c
c     entry for update calls needed in lksubgridsolve
c
      ENTRY update (drscr)
c
c      call initialupdate (dra,drb,drscr,m,n,lErrorReport) ! use initial matrix to generate modifications of known vector arrays
c
c      write (iluer,3457) drscr(1:m)
c 3457 format (' solut3457: drscr',10(d14.7))
      j=0
c      write (ilume,3456)
c 3456 format ('solut3456: we are now updating after a sub-cell solution'
c     & ' for a single MODFLOW cell'
c     &,/,' We are in the ENTRY update in SOLUT.')
      call timer(iticks3)
      call lsupdate (drscr,j,m,n,.FALSE.,0)
      call dbupdate (drscr,j,m,n,czi_long,calph,.FALSE.,0)
      call wlupdate (drscr,j,m,n,czi_long,.FALSE.,0)
      call w3update (drscr,J,m,n,czi_long,CALPH,ITYPE,.FALSE.,0)
      call lkupdate (drscr,j,m,n,czi_long,.FALSE.,0)
      call gvupdate (drscr,j,m,n)
        IF (J.NE.m) then
          write (ilume,2900) N,J  ! j is counter of known vector elements that are updated.
          AMESS(1)='Known vector array update out of sync.'
          AMESS(2)='Contact technical support.'
          CALL HALT(2) ! stop program execution for batch version
        endif
       call timer (iticks4)
       iticks=iticks4-iticks3
       if (ltimer) write (ilume,9077) iticks
 9077 format(' Known vector arrays update (no A*b) execution time=',i10,
     &' E-2 seconds.')
      return
c
 1000 format (/,' Solution # ',I3,'. Number of equations: ',I4)
 1100 FORMAT (/,' Solution # ',I3,'. Number of equations: ',I4)
 1200 format (/,' ERROR in SOLUT: In the presence of iterface flow ',
     & 'the "save decomposed matrix to disk" option should not be used.'
     &,/,' This option has been disabled.')
 1500 FORMAT (/,' Maximum errors in boundary conditions:
     &      ')
 1900 format (' ***ERROR in SOLUT during colloc. vector construction',/,
     &' Number of equations ',i6,' not equal to array dimension ',i6,/,
     &' Number of strengths ',i6,' not equal to array dimension ',i6,/,
     &' Execution aborted.')
 2000 FORMAT (' Redo solution to check for percolating surface waters.')
 2900 FORMAT (' ***ERROR in SOLUT: known vector update out of sync.',\,
     &' Number of equations=',I5,' Number updates=',I5)
 2901 FORMAT (' ***ERROR in SOLUT: known vector check out of sync.',\,
     &' Number of equations=',I5,' Number updates=',I5)
 2905 format (' Check on update routines prior to solution process.')
 2910 format (' Check on update routines prior to BC check process.')
 2950 FORMAT (' ***ERROR in SOLUT: keep strengths array out of sync.',\,
     &' Number of equations=',I5,' Number strengths=',I5)
 3900 format (' ***ERROR in SOLUT during matrix construction.',/,
     &' Number of equations ',i6,' not equal to array dimension ',i6,/,
     &' Number of strengths ',i6,' not equal to array dimension ',i6,/,
     &' Execution aborted.')
 3910 format (' ***ERROR in SOLUT during matrix correction.',\,
     &' Number of equations=',I5,' Number of unknowns=',I5)

 4000 FORMAT (' ***ERROR in SOLUT: too many equations. (max.=',I4,')')
 4100 format (' ***ERROR is SOLUT: known vector mismatch.',\,
     &' Number of equations=',I5,' Length of known vector=',I5)
 4200 format (' ***ERROR is SOLUT: solution vector mismatch.',\,
     &' Number of equations=',I5,' Length of solution vector=',I5)
 4500 FORMAT (' ***ERROR in SOLUT: singular matrix.',/,
     &' Likely cause coinciding control points, e.g. center of a',
     &' line sink',/,' on top of a line doublet center or end.',
     &' Correct input data.')
 6000 FORMAT (' Matrix condition number = ',E11.4)
 7000 format (' ')
 7500 format (' Starting subgrid iterations, number of iteration=',i4)
 8000 format (' Error in SOLUT: Failed to allocate CZKEEP(',i4,').')
 8001 format (' Error in SOLUT: Failed to allocate LSKIP(',i4,').')
      END
c
c --------------------------------------------------------------------------------
c
      subroutine initialupdate(dra,drb,drscr,m,n,lErrorReport)
c
c --------------------------------------------------------------------------------
c
c     Routine generates corrections on the known vector arrays of the analytic elements
c
c     IMPORTANT: DRA must contain the initial coefficient matrix and DRB must contain
c                the strength increments (resulting from the latest solution)
c                DRSCR will get the potential or normal flow increments to be used
c                in the varies "__update" routines called next.
c
      implicit none
      INTEGER m,n,i,
     &        nsolOut
      LOGICAL lErrorReport,
     &        lsolOut,loadsolOut,linalreadyOut,lErrorReportOut,
     &        lDirectfromDiskOut
      CHARACTER(8) aBasenameOut
      CHARACTER(16)aDateTimeOut
      REAL(8) dra,drb,drscr
      DIMENSION dra(m,n),drb(n),drscr(m)
      include 'lusys.inc'
c
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      if (nsolOut.eq.0) return
c
      drscr=MATMUL(dra,drb)
c
c      write (ilume,1001) (i,drscr(i),i=1,m)
c 1001 format ('initialupdate1: drscr(',i3,')=',d14.7)
      end subroutine
c
c --------------------------------------------------------------------------------
c
      subroutine condition(DRA,DRSCR,m,n)          ! not used
c
c --------------------------------------------------------------------------------
c
      implicit none
      INTEGER(4) i,j,n,m
      REAL(8) DRA,DRSCR
      DIMENSION DRA(m,*),DRSCR(*)
      do i=1,m
      DRSCR(i)=DRA(i,i)
      do j=1,n
      DRA(i,j)=DRA(i,j)/DRSCR(i)
      end do
      end do
      return
      end subroutine
c
c --------------------------------------------------------------------------------
c
      subroutine adjust(DRB,DRSCR,n)             ! not used
c
c --------------------------------------------------------------------------------
c
      implicit none
      INTEGER(4) i,n
      REAL(8) DRB,DRSCR
      DIMENSION DRB(*),DRSCR(*)
      do i=1,n
      DRB(i)=DRB(i)/DRSCR(i)
      end do
      return
      end subroutine
c
c --------------------------------------------------------------------------------
c
       SUBROUTINE DECOMP (N,A,C,P,W,RN)
c
c --------------------------------------------------------------------------------
c
C
C     DECOMPOSES A REAL MATRIX BY GAUSSIAN ELIMINATION
C     AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     USE SOLVE TO COMPUTE SOLUTIONS TO LINEAR SYSTEMS.
C
C     INPUT-
C
C       N=ORDER OF THE MATRIX.
C
C       A=MATRIX TO BE TRIANGULARIZED.
C
C     OUTPUT-
C
C        VERSION OF A LOWER TRIANGULAR MATRIX I-L SO THAT
C        (PERMUTED MATRIX)*A=L*U
C
C       C=AN ESTIMATE OF THE CONDITION OF A.
C        FOR THE LINEAR SYSTEM A*X=B, CHANGES IN A AND B
C        MAY CAUSE CHANGES C TIMES AS LARGE IN X.
C        IF C+1.0 .EQ. C, A IS SINGULAR TO WORKING PRECISION.
C        C IS SET TO 1.0E+32 IF EXACT SINGULARITY IS DETECTED.
C
C       P=THE PIVOT VECTOR.
C        P(K)=THE INDEX OF THE K-TH PIVOT ROW
C        P(N)=(-1)**(THE NUMBER OF INTERCHANGES)
C
C
C     WORK SPACE-
C      THE VECTOR W MUST BE DECLARED AND INCLUDED IN THE CALL.
C      ITS INPUT CONTENTS ARE IGNORED. ITS OUTPUT CONTENTS ARJ_
C       USUALLY UNIMPORTANT.
C
C     THE DETERMINANT OF A CAN BE OBTAINED ON OUTPUT BY
C       D(A)=P(N)*A(1,1)*A(2,2)* ... *A(N,N)
C
C***********************************************************************
C
C     PROGRAM HISTORY:
C     TRANSCRIBED FROM - "COMPUTER METHODS FOR MATHEMATICAL COMPUTATION"
C     AUTHORS - G.E.FORSYTHE, M.A.MALCOLM, AND C.B.MOLER
C     BY - T. GRAY CURTIS
C     DATE - 1 NOVEMBER, 1979
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER(4) I,J,K,KB,KM1,KP1,M,N,P(N),NM1
      REAL(8) A(N,N),AN,C,E,RN(N),T,W(N),YN,ZN
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
c
c        call matrixrowout  (a,n,938,'ra in decomp 938')
c        call matrixrowout (a,n,1145,'ra in decom 1145')
c
      P(N)=1
      IF(N.EQ.1)GOTO 2
      NM1=N-1
C
C     COMPUTE L1-NORMS OF A
      DO 8 I=1,N
      RN(I)=0.
  8   CONTINUE
      AN=0.
      DO 10 J=1,N
      T=0.
      DO 9 I=1,N
      T=T+ABS(A(I,J))
      RN(I)=RN(I)+ABS(A(I,J))
  9   CONTINUE
      IF(T.GT.AN)AN=T
 10   CONTINUE
      DO 11 I=1,N
      RN(I)=1./RN(I)
 11   CONTINUE
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      DO 24 K=1,NM1
      KP1=K+1
!      WRITE (ILUME,1000) KP1
C
C     FIND PIVOT
C
      M=K
      DO 20 I=KP1,N
      IF (ABS(A(I,K))*RN(I).GT.ABS(A(M,K))*RN(M)) M=I
 20   CONTINUE
      P(K)=M
      IF(M.NE.K)P(N)=-P(N)
      T=A(M,K)
      A(M,K)=A(K,K)
      A(K,K)=T
C
C     SKIP STEP IF PIVOT IS ZERO
C
      IF(T.EQ.0.)GOTO 24
C
C     COMPUTE MULTIPLIERS
C
      DO 21 I=KP1,N
      A(I,K)=-A(I,K)/T
 21   CONTINUE
C
C     INTERCHANGE AND ELIMINATE BY COLUMNS
C
      DO 23 J=KP1,N
      T=A(M,J)
      A(M,J)=A(K,J)
      A(K,J)=T
      IF(T.EQ.0.)GOTO 23
      DO 22 I=KP1,N
      A(I,J)=A(I,J)+A(I,K)*T
 22   CONTINUE
 23   CONTINUE
 24   CONTINUE
c                 The remainder of the logic is to calculate the
c                 matrix condition number. This adds almost very
c                 little time to the overall execution time of
c                 this routine. (usually less than 1%)
C
C     C=(1-NORM OF A)*(AN ESTIMATE OF 1-NORM OF A-INVERSE)
C     ESTIMATE OBTAINED BY ONE STEP OF INVERSE ITERATION
C     FOR THE SMALL SINGULAR VECTOR.
C     THIS INVOLVES SOLVING TWO SYSTEMS OF EQUATIONS,
C     (A-TRANSPOSE)*Y=E AND A*Z=Y WHERE E IS A VECTOR OF +1 OR -1
C     CHOSEN TO CAUSE GROWTH OF Y.
C     ESTIMATE=(1-NORM OF Z)/(1-NORM OF Y)
C
C     SOLVE (A-TRANSPOSE)*Y=E
C
      DO 31 K=1,N
      T=0.
      IF(K.EQ.1)GOTO 1
      KM1=K-1
      DO 30 I=1,KM1
      T=T+A(I,K)*W(I)
 30   CONTINUE
  1   E=1.
      IF(T.LT.0.)E=-1.
      IF(A(K,K).EQ.0.)GOTO 3
      W(K)=-(E+T)/A(K,K)
 31   CONTINUE
      DO 41 KB=1,NM1
      K=N-KB
      T=0.
      KP1=K+1
      DO 40 I=KP1,N
      T=T+A(I,K)*W(K)
 40   CONTINUE
      W(K)=T
      M=P(K)
      IF(M.EQ.K)GOTO 41
      T=W(M)
      W(M)=W(K)
      W(K)=T
 41   CONTINUE
C
      YN=0.
      DO 50 I=1,N
      YN=YN+ABS(W(I))
 50   CONTINUE
C
C     SOLVE A*Z=Y
C     FOR THE VECTOR Z,
C     AND COMPUTE ITS NORM, ZN
C
      CALL SOLVE(N,A,W,P)
C
      ZN=0.
      DO 60 I=1,N
      ZN=ZN+ABS(W(I))
 60   CONTINUE
C
C     ESTIMATE CONDITION
C
      C=AN*ZN/YN
      IF(C.LT.1.0)C=1.
      RETURN
C
C     1-BYE-1
C
  2   C=1.0
      IF(A(1,1).NE.0.)RETURN
C
C     EXACT SINGULARITY
C
  3   C=1.0E+32
      RETURN
 1000 FORMAT ('+Matrix solution, doing equation #: ',I3)      
      END
c
c --------------------------------------------------------------------------------
c
      SUBROUTINE SOLVE(N,A,B,P)
c
c --------------------------------------------------------------------------------
c
C***********************************************************************
C     SOLUTION OF LINEAR SYSTEM A*X=B
C     DO NOT USE IF DECOMP HAS DETECTED SINGULARITY
C
C     INPUT-
C       N=ORDER OF MATRIX
C       A=TRANGULARIZED MATRIX OBTAINED FROM DECOMP
C       B=RIGHT HAND SIDE VECTOR
C       P=PIVOT VECTOR OBTAINED FROM DECOMP
C
C    OUTPUT-
C       B=SOLUTION VECTOR, X.
C**********************************************************************
C
C     PROGRAM HISTORY:
C     TRANSCRIBED FROM - "COMPUTER METHODS FOR MATHEMATICAL COMPUTATION"
C     AUTHORS - G.E.FORSYTHE, M.A.MALCOLM, AND C.B.MOLER
C     BY - T.GRAY CURTIS
C     DATE - 1 NOVEMBER, 1979
C
C*********************************************************************
      IMPLICIT NONE
      INTEGER(4) I,K,KB,KP1,M,N,P(N),NM1,NMKB
      REAL(8) A(N,N),B(N),T
      INCLUDE 'TRACOM.INC'
C
C     FORWARD ELIMINATION
C
      IF(N.EQ.1)GOTO 1
      NM1=N-1
      DO 11 K=1,NM1
      KP1=K+1
      M=P(K)
      T=B(M)
      B(M)=B(K)
      B(K)=T
      DO 11 I=KP1,N
      B(I)=B(I)+A(I,K)*T
 11   CONTINUE
C
C     BACK SUBSTITUTION
C
      DO 21 KB=1,NM1
      NMKB=N-KB
      K=NMKB+1
      B(K)=B(K)/A(K,K)
      T=-B(K)
      DO 21 I=1,NMKB
      B(I)=B(I)+A(I,K)*T
 21   CONTINUE
  1   B(1)=B(1)/A(1,1)
      RETURN
      END
c
c --------------------------------------------------------------------------------
c
      subroutine Gauss_Seidel(n,dra,drb,drfac,niter,drs)        ! not used
c
c --------------------------------------------------------------------------------
c
c     !!!!! NOT YET TESTED !!!!!
c
c     Routine performs a Gauss Seidel iteration with selective relaxation.
c     rdfac(4,i) is an array with relaxation factors, which is specific to each equation
c
c     Input:
c            imxsize array dimension of dra(.,.)
c            n     number of equations
c            dra   coefficient matrix
c            drb   known vector
c            drfac drfac(4,i) serves as the relaxation factor for the ith equation
c            niter number of iterations (both forward and backward sweep)
c
c     Output:
c            drs   solution vector.
c
c
      implicit none
      INTEGER(4) n,niter,istart,iend,istep,itemp,kend,
     &           i,j,k
      REAL(8) dra,drb,db,drs,rdels,drfac
      DIMENSION dra(n,n),drb(*),drs(*),drfac(4,*)
      include 'lusys.inc'
c
c      write (iluer,999) (i,(dra(i,j),j=1,n),i=1,n)    ! for 5 equations only
c 999  format (' GaussSeidel ',/,5('dra(',i3,',j)=',5(d14.7,1x),/))
      kend=niter
      istart=1
      iend=n
      istep=1
c      write (iluer,1000) (i,drs(i),i=1,n)
 1000 format (' GaussSeidel0: i,drs ',i3,d14.7)
      do k=1,kend
c      WRITE(*,1001) k
c      write (iluer,1001) k
 1001 FORMAT(' GaussSeidel1: iteration #=',i3)
c      write (iluer,1003) kend,relax
 1003 format (' GaussSeidel3: kend,relax ',i4,d14.7)
        do i=istart,iend,istep
c          if (drfac(4,i).ne.0.0d0) then    ! may save time ..
            db=drb(i)
            do j=1,n
             db=db-dra(i,j)*drs(j)
             end do
             rdels=db/dra(i,i)
             drs(i)=drs(i)+drfac(4,i)*rdels
c             write (iluer,1002) i,rdels,drs(i)
 1002 format (' GaussSeidel2: i,rdels,drs(i) ',i5,2(d14.7))
c          end if
        end do
      end do
c
c     debugging
c
c      do i=1,n
c      WRITE(iluer,1005) i,drfac(4,i),drs(i)
c      end do
 1005 format (' GaussSeidel5: i,drfac,drs ',I4,1x,2(d14.7,1x))
      return
      end subroutine
c
c -----------------------------------------------------------------------------------
c
      subroutine CorrectKnownVector (n,drb,drfac,dra,drscr)     ! not used
c
c -----------------------------------------------------------------------------------
c
c    Routine forces the known vector element the same as the equation for skipped linesinks.
c    The strength parameters in DRSCR(i) are the strength corrections coming out of SOLVE
c    When DRFAC(4,i)=0 the strength correction should be zero, if not reset RDB(i), and resolve
c
c    INPUT:
c        n       size of (square) matrix
c        drb     known vector
c        drfac   if drfac(4,i)=0.0 equation number i should be skipped
c        dra     coefficient matrix
c        drscr   solution vector
c
c    OUTPUT:
c        drb     all elements for which drfac(4,i)=0.0 are corrected
c

      implicit none
      INTEGER(4) n,i,j
      REAL(8) drb,drfac,dra,drscr
      DIMENSION drb(*),drfac(4,*),dra(n,*),drscr(*)
      include 'lusys.inc'
c
      do i=1,n
      write (iluer,1000) i,drscr(i),drb(i)
 1000 format(' CorrectKnownVector0: i,drscr(i),drb(i) ',
     &        i4,1x,d14.7,1x,d14.7)
      end do
      do i=1,n
      if (drfac(4,i).eq.0.0d0) THEN ! line sink equation to be skiped
        write (iluer,1001) i,drscr(i) ! should be close to zero when converged
 1001   format (' CorrectKnownVector1: i,drscr(i) ',i4,1x,d14.7)
        drb(i)=0.0d0
        do j=1,n
        drb(i)=drb(i)+dra(i,j)*drscr(j)
        end do
        drb(i)=drb(i)-dra(i,i)*drscr(i)
c        drb(i)=0.5*drb(i) ! try relaxation
        write (iluer,1002) i,drb(i)
 1002   format (' CorrectKnownVector2: i,drb(i) ',i4,1x,d14.7)
      end if
      end do
      end subroutine
c
c --------------------------------------------------------------------------------
c
      SUBROUTINE SpecialSolve(N,A,B,P,drfac)       ! not used
c
c --------------------------------------------------------------------------------
c
c     Special experimental routine for DirectFromDisk solutions
c
C***********************************************************************
C     SOLUTION OF LINEAR SYSTEM A*X=B
C     DO NOT USE IF DECOMP HAS DETECTED SINGULARITY
C
C     INPUT-
C       N=ORDER OF MATRIX
C       A=TRANGULARIZED MATRIX OBTAINED FROM DECOMP
C       B=RIGHT HAND SIDE VECTOR
C       P=PIVOT VECTOR OBTAINED FROM DECOMP
C
C    OUTPUT-
C       B=SOLUTION VECTOR, X.
C**********************************************************************
C
C     PROGRAM HISTORY:
C     TRANSCRIBED FROM - "COMPUTER METHODS FOR MATHEMATICAL COMPUTATION"
C     AUTHORS - G.E.FORSYTHE, M.A.MALCOLM, AND C.B.MOLER
C     BY - T.GRAY CURTIS
C     DATE - 1 NOVEMBER, 1979
C
C*********************************************************************
      IMPLICIT NONE
      INTEGER(4) I,K,KB,KP1,M,N,P(N),NM1,NMKB
      REAL(8) A(N,N),B(N),T,drfac
      DIMENSION drfac(4,*)
      INCLUDE 'TRACOM.INC'
C
C     FORWARD ELIMINATION
C
      IF(N.EQ.1)GOTO 1
      NM1=N-1
      DO 11 K=1,NM1
      KP1=K+1
      M=P(K)
      T=B(M)
      B(M)=B(K)
      B(K)=T
      DO 11 I=KP1,N
      B(I)=B(I)+A(I,K)*T
 11   CONTINUE
C
C     BACK SUBSTITUTION
C
      DO 22 KB=1,NM1
      NMKB=N-KB
      K=NMKB+1
      if (drfac(4,k).eq.0.0d0) then
        B(K)=0.0d0 ! force s to zero when no equation
      else
        B(K)=B(K)/A(K,K)
        T=-B(K)
        DO 21 I=1,NMKB
        B(I)=B(I)+A(I,K)*T
 21     CONTINUE
      endif
 22   continue
  1   B(1)=B(1)/A(1,1)
      RETURN
      END
c
c -------------------------------------------------------------------------------------------------------
c
      SUBROUTINE Sherman_Morrison (n,dra,drb,ipiv,lskip,
     &                             lErrorReport,ltimer)
c
c -------------------------------------------------------------------------------------------------------
c
c
c     Routine generates a solution using an existing decomposed (LU) matrix and
c     the repeated application of the Sherman Morrison formula to remove several
c     equations from the original set.
c
c     Input:
c            n       size of the original set of equations
c            dra     n*n matrix used to store the original coefficient matrix and the LU matrix
c            drb     known vector for the original set of equations. Note: drb(i)=0.0 if equation i is removed.
c            ipiv    pivot vector for the original set of equations.
c            lskip   lskip(i)=.true. is equation i is to be removed
c            lErrorReport=.true. report activity
c            ltimer  ltimer=.true. write execution time reports
c
c     Output:
c            drb     solution vector
c
c
c
      implicit none
      LOGICAL lskip,ltimer,lErrorReport
      INTEGER(4) n,ipiv,i,j,nskip,ipoint,ierr
      REAL(8) dra,drb,drv,dru,druarray,dralpha,drlabda
      DIMENSION dra(n,n),drb(n),ipiv(n),lskip(n)
      ALLOCATABLE drv(:),dru(:),druarray(:,:),ipoint(:)
      include 'lusys.inc'
c
c     allocate scratch arrays
c
      nskip=0
      do i=1,n
      if (lskip(i)) then
        nskip=nskip+1
      end if
      end do
      write (*,2000) nskip
      if (lErrorReport)write (ilume,2000) nskip
      ALLOCATE (drv(n),dru(n),druarray(n,nskip),ipoint(nskip),STAT=ierr)
      if (ierr.ne.0) then
        deallocate (drv,dru,druarray,ipoint)
        write (ilume,1000) n,nskip
        AMESS(1)='Error in Sherman_Morrison routine.'
        AMESS(2)='Failed to allocate scratch arrays. Try to solve'
        AMESS(3)='without using decomposed matrices stored to disk.'
        CALL HALT(3)   ! stop program execution for batch version
      end if
c
c     set up pointer array to equations to be removed
c
      nskip=0
      do i=1,n
        if (lskip(i)) then
          nskip=nskip+1
          ipoint(nskip)=i
        end if
      end do
c
c     create initial solutions
c
      call LoadDecompMatrix (dra,ipiv,n,lErrorReport,ltimer)
c
      call Solve  (n,dra,drb,ipiv) ! dra now contains vector X superscript 1
c
      do i=1,nskip
        dru(1:n)=0.0d0
        dru(ipoint(i))=1.0d0        ! dru is now the vector U superscript i
        call Solve (n,dra,dru,ipiv) ! dru now contains vector U' superscript i subscript 1
        druarray(1:n,i)=dru(1:n)    ! store for later use
!         Note, this may be simplified by initializing druarray(1:n,i) and calling it directly with Solve.
!         In fact, there may be Solve routines that can be handed DRU(N,NSKIP) all at once and return all
!         NSKIP vectors DRU.
      end do   ! we are done calling Solve, now load the coefficient matrix
      call LoadMatrix (dra,n,n,lErrorReport,ltimer)   ! Note: this is tricky!!
c      We are loading only "n" rows from the matrix.
c      We are picking some rows from the matrix for which lskip=.true.
c      However, "lskip(i>n)" should be false (no line-sinks).
c      This must be tested.
      do i=1,nskip
        drv(1:n)=dra(ipoint(i),1:n)
        drv(ipoint(i))=1.0d100
        dru(1:n)=druarray(1:n,i)
        drlabda=DOT_PRODUCT(drv,dru)
        dralpha=DOT_PRODUCT(drv,drb)/(1.0d0 - drlabda)
        drb(1:n)=drb(1:n)+dralpha*dru(1:n) ! solution vector X after equation i has been removed
        if (i.lt.nskip) then ! generate new U' vectors
          do j=i+1,nskip
            dru(1:n)=druarray(1:n,j)
            dralpha=DOT_PRODUCT(drv,dru)/(1.0d0 - drlabda)
            druarray(1:n,j)=druarray(1:n,j) + dralpha*druarray(1:n,i)
          end do
        end if
      end do
      return
c
 1000 format (' ERROR in Sherman_Morrison routine: ',/,
     &        'Cannot allocate a ',i4,' by ',i4,' scratch array.')
 2000 format (' Number of equations removed are ',i4)
      END SUBROUTINE
c
C ---------------------------------------------------------------------------------
C
      subroutine gfread_convergence_file
c
c     Routine reads the file "converge.tab" if it exists. If the table does not exists the
c     convergence criteria are set to 0.0, which means they will never be met.
c     If the table does exists it passes a minimim number of iterations to the routine SOLUT
c     and passes the relevant convergence criteria to the routines XXERROR. If all convergence
c     criteria are met and the minimum number of iterations is met, the solution procedure
c     will be terminated and the current solution saved, etc.
c
c     Example and format of the file "converge.tab":
c
c     * comment statement (* in first column)
c     *analytic element group    convergence criterion
c        minimum iterations        3
c        reference point           0.00001
c        linesinks dirichlet       0.00001
c        linesinks resistance      0.01
c        inhomogeneity domains     0.001
c        barrier noflow            0.001
c        barrier resistance        0.001
c        well 3D                   0.001
c        well 2D                   0.001
c        lake water balance        0.1
c     quit
c
c     Note: The criteria for 2D wells apply to head specified wells, while for
c           3D wells it applies to both head and discharge specified wells.
c
c     If all values are set to 0.0 its effect is that no convergence criteria are
c     set and that the maximum number of iterations (as set by the GUI) will be applied.
c

      implicit none
      INTEGER ilu,itemp,i,ierr,nconvergearraylength,igftablelength
      LOGICAL lret,lnoconvergefile
      REAL(8) rdum,rvar,rconvergearray
      parameter (nconvergearraylength=10)
      DIMENSION rconvergearray(nconvergearraylength)
      include 'lusys.inc'
      include 'match.inc'
      include 'main.inc'
      include 'tracom.inc'
c
      lnoconvergefile=.false.
      igftablelength=0
      if (.not.lErrorReport) then ! no error reporting, cannot use convergence criteria
      write (ilume,1500)
      GOTO 25
      endif
      ilu=2
      afile='converge.tab'
      itemp=iluin ! temporarily store iluin and replace for this read
c      call opafil(ilu,-2,lret)
        iluin=ilu
      OPEN (UNIT=ILU,FILE=AFILE,STATUS='OLD',ERR=20,IOSTAT=IERR)
  10    call inline
        if (aline(1).eq.'Q'.or.aline(1).eq.'q') GOTO 20 ! end of data in file
        rdum=rvar(3) ! numeric data is always third item on the line
        if (.not.lerror) then
          igftablelength=igftablelength+1
          if (igftablelength.le.nconvergearraylength) then
            rconvergearray(igftablelength)=rdum
            GOTO 10
          else
            write (iluer,2000)
            igftablelength=0
          end if
        end if
        if (lerror) then
          write (iluer,3000)
          lerror=.false.
        end if
  20    close (iluin)
        iluin=itemp ! restore iluin
        if (ierr.ne.0) then  ! failed to open file, assume it is not present
          write (ilume,1000)
          lnoconvergefile=.true.
        endif
  25    if (igftablelength.eq.0.or.lnoconvergefile
     &   .or..not.lErrorReport) THEN ! set criteria to 0.0
          minimum_iterations=0
          rconverge_reference=0.0
          rconverge_linesinks=0.0
          rconverge_linesinks_resistance=0.0
          rconverge_inhomogeneities=0.0
          rconverge_barriers_noflow=0.0
          rconverge_barriers_resistance=0.0
          rconverge_well_3D=0.0
          rconverge_well_2D=0.0
          rconverge_lake_waterbalance=0.0
        end if
      if (.not.lnoconvergefile.and.
     &    igftablelength.eq.nconvergearraylength
     &     .and.lErrorReport) then ! store criteria
        minimum_iterations=INT(rconvergearray(1))
        rconverge_reference=rconvergearray(2)
        rconverge_linesinks=rconvergearray(3)
        rconverge_linesinks_resistance=rconvergearray(4)
        rconverge_inhomogeneities=rconvergearray(5)
        rconverge_barriers_noflow=rconvergearray(6)
        rconverge_barriers_resistance=rconvergearray(7)
        rconverge_well_3D=rconvergearray(8)
        rconverge_well_2D=rconvergearray(9)
        rconverge_lake_waterbalance=rconvergearray(10)
        write (ilume,4000)
        WRITE(ilume,5000)  minimum_iterations
        WRITE(ilume,6000)  rconverge_reference
        WRITE(ilume,7000)  rconverge_linesinks
        WRITE(ilume,8000)  rconverge_linesinks_resistance
        WRITE(ilume,9000)  rconverge_inhomogeneities
        WRITE(ilume,9200)  rconverge_barriers_noflow
        WRITE(ilume,9300)  rconverge_barriers_resistance
        WRITE(ilume,9500)  rconverge_well_3D
        WRITE(ilume,9700)  rconverge_well_2D
        write(ilume,9800)  rconverge_lake_waterbalance
      end if
      return
 1000 format (/,' No file "converge.tab" found.',/,
     &        ' No convergence criteria will be applied.',/)
 1500 format (/,' Error reporting suppressed for fast solve.',/,
     &        ' Cannot apply convergence criteria.',/)
 2000 format(' ***ERROR in the file "converge.tab"!',/,
     &       'Too many convergence factors in file',/,
     &       'No convergence criteria will be applied')
 3000 format (' ***ERROR in gfread_convergence_file:',/,
     &        'Illegal syntax in the file "converge.tab"')
 4000 format(/,' Found file "converge.tab"',
     &       ' The following convergence criteria will be used:',/)
 5000 format (' minimum iterations             ',i5)
 6000 format (' reference point                ',d14.7)
 7000 format (' linesinks dirichlet            ',d14.7)
 8000 format (' linesinks resistance           ',d14.7)
 9000 format (' inhomogeneity domains          ',d14.7)
 9200 format (' barriers noflow                ',d14.7)
 9300 format (' barriers resistance            ',d14.7)
 9500 format (' wells 3D                       ',d14.7)
 9700 format (' wells 2D                       ',d14.7)
 9800 format (' lake water balance             ',d14.7,/)
      end subroutine
c
c ----------------------------------------------------------------------------
c
      subroutine reshapematrix (dra,m,n)
c
c     The array DRA(M,N), with M>N contains the square matrix A(N,N)
c     This routine shifts the elements in DRA so that DRA(N,N) contains A(N,N)
c     Example:
c     DRA(6,4)
c                  |1 1 1 1|
c                  |2 2 2 2|
c                  |3 3 3 3|
c                  |4 4 4 4|
c                  |* * * *|
c                  |* * * *|
c
c     DRA(*) =  |1234**1234**1234**1234**|
c
c     after reshape:
c
c     DRA(*) =  |1234123412341234|
c
c     so that
c
c     DRA(4,4)
c
c                  |1 1 1 1|
c                  |2 2 2 2|
c                  |3 3 3 3|
c                  |4 4 4 4|
c
c
      implicit none
      INTEGER i,j,m,n,ioffset1,ioffset2
      REAL(8) dra
      DIMENSION dra(*)
      include 'lusys.inc'
c
      if (m.gt.n) then  ! else assume m=n and we are done
        write (iluer,1001) m,n
 1001 format (' reshapematrix1: m=',i6,' n=',i6)
        ioffset1=0
        ioffset2=0
        do i=1,n-1
         ioffset1=ioffset1+n
         ioffset2=ioffset2+m
         do j=1,n
          dra(j+ioffset1)=dra(j+ioffset2)
         end do
        end do
      end if
      return
      end subroutine
c
c ------------------------------------------------------------------------
c
      subroutine vectorout (rvec,nsize,aname)
c
c     Routine write a vector for debugging purposes
c
      implicit none
      INTEGER nsize,i
      CHARACTER*16 aname
      REAL(8) rvec
      DIMENSION rvec(*)
      include 'lusys.inc'
c
      if (nsize.gt.0) then
        write (iluer,1000) aname
 1000   format (' From "vectorout": writing ',a16)
        do i=1,nsize,10
          write (iluer,2000) rvec(i:i+9)
 2000     format (10(d14.7))
        end do
      else
        write (iluer,3000) nsize
 3000   format ('From "vectorout": nsize=',i5,' is out of range.')
      end if
      return
c
      end subroutine
c
c ------------------------------------------------------------------------
c
      subroutine complexvectorout (cvec,nsize,aname)
c
c     Routine write a vector for debugging purposes
c
      implicit none
      INTEGER nsize,i
      CHARACTER*16 aname
      complex(8) cvec
      DIMENSION cvec(*)
      include 'lusys.inc'
c
      if (nsize.gt.0) then
        write (iluer,1000) aname
 1000   format (' From "complexvectorout": writing ',a16)
        do i=1,nsize,5
          write (iluer,2000) cvec(i:i+4)
 2000     format (5(2(d14.7),3x))
        end do
      else
        write (iluer,3000) nsize
 3000   format('From "complexvectorout": nsize=',i5,' is out of range.')
      end if
      return
c
      end subroutine
c
c ------------------------------------------------------------------------
c
      subroutine diagonalelementsout (dra,nsize,aname)
c
c     Routine writes the diagonal elements of square matrix dra for debugging purposes
c
      implicit none
      INTEGER nsize,i,j,nend
      CHARACTER*16 aname
      REAL(8) dra
      DIMENSION dra(nsize,nsize)
      include 'lusys.inc'
      if (nsize.gt.0) then
        write (iluer,1000) aname
 1000   format (' From "diagonalelementsout": writing ',a16)
        nend=(nsize/10)*10
        do i=1,nend,10
          write (iluer,2000) (dra(j,j),j=i,i+9)
 2000     format (10(d14.7))
        end do
        write (iluer,2000) (dra(j,j),j=nend+1,nsize)
      else
        write (iluer,3000) nsize
 3000   format('From "diagonalelementsout": nsize=',i5,' out of range.')
      end if
      return
c
      end subroutine
c
c ------------------------------------------------------------------------
c
      subroutine matrixrowout (dra,nsize,irow,aname)
c
c     Routine writes the diagonal elements of square matrix dra for debugging purposes
c
      implicit none
      INTEGER nsize,i,j,irow
      CHARACTER*16 aname
      REAL(8) dra
      DIMENSION dra(nsize,nsize)
      include 'lusys.inc'
      if (nsize.gt.0) then
        write (iluer,1000) aname
 1000   format (' From "matrixrowout": writing ',a16)
          write (iluer,2000) dra(irow,1:nsize)
 2000     format (10(d14.7))
      else
        write (iluer,3000) nsize
 3000   format('From "matrixrowout": nsize=',i5,' out of range.')
      end if
      return
c
      end subroutine
c
c --------------------------------------------------------------------------
c
      subroutine comparematrices (dra,isize,jsize,aname)
c
c     read matrices from disk and compare to each other
c     hand current matrix, read old matrix.
c
      implicit none
      INTEGER isize,jsize,i,j,ierr
      LOGICAL lErrorReport,ltimer
      CHARACTER*16 aname
      REAL(8) dra,scratch,rdiff
      DIMENSION dra(isize,jsize)
      ALLOCATABLE scratch(:,:)
      include 'lusys.inc'
c
      lErrorReport=.true.
      ltimer=.false.
      write (iluer,1000) aname,isize
 1000 format (' From "comparematrices" ',A16,' size=',i5)
      ALLOCATE (scratch(isize,jsize),STAT=ierr)
      if (ierr.ne.0) then
        deallocate (scratch)
        write (ilume,1000) isize
        AMESS(1)='Error in comparematrices routine.'
        AMESS(2)='Failed to allocate scratch matrix. Abort.'
        CALL HALT(2)   ! stop program execution for batch version
      end if
      call LoadMatrix (scratch,jsize,jsize,lErrorReport,ltimer) ! only read a square matrix
      DO i=1,isize
      do j=1,jsize
       rdiff=ABS(dra(i,j)-scratch(i,j))
       if (rdiff.gt.0.00001d0) then
         WRITE (iluer,2000) i,j,rdiff
 2000 format (' matrix elements (',i4,',',i4') differ by ',d14.7)
       end if
      end do
      END DO
c
      return
c
      end subroutine


