C     Last change:  HMH   8 Aug 2005    2:13 pm
c       This file contains the following routines and functions:
c
C	SUBROUTINE NOCONJUNCTIVE   sets flag for conjuctive solutions false
c       SUBROUTINE LSPREP          calculate width, CLSCONST, and cancel recharge on lakes
c       real function rflsnearflow returns the streamflow in the nearest line-sink
c       subroutine ls_setwidth     routine sets the representative width parameter for resistance line-sinks
c       subroutine ls_cancellakerecharge  add recharge inhomogeneity to cancel lake recharge
C
C-----------------------------------------------------------------------------------------------------
C    
      SUBROUTINE NOCONJUNCTIVE
C
C-----------------------------------------------------------------------------------------------------
C
      IMPLICIT NONE
      INCLUDE 'lscom.inc'     
C
      LSBASE=.FALSE.
C
      RETURN
      END
c
c ----------------------------------------------------------------------------------------------------
c
      subroutine lsprep()
c
c ----------------------------------------------------------------------------------------------------
c
c     Routine prepares line-sink functions for solution process by:
c     1) calculating the correct width parameter
c        --> call ls_setwidth
c     2) calculating line-sink specific constants for the potential function
C     that makes the contribution of each line-sink vanish at the reference point.
c     3) adding a recharge inhomogeneity for each lake feature that cancels aquifer recharge
c        --> call ls_cancellakerecharge
c     Routine is called once in SOLUT prior to the first solution step.
c
      implicit none
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
      INTEGER ils
      REAL RK0,RH0,RHED0,RB0,RP0
      COMPLEX czs,cze,cz0,comls
c
      call ls_setwidth() ! calculate effective leake zones for line-sinks representing surface waters with bottom resistance
      call GVPAR (RK0,RH0,RHED0,RB0,RP0,CZ0)
      do ils=1,nls
      czs=CLSZS(ils)
      cze=CLSZE(ils)
      clsconst(ils)=-comls(cz0,czs,cze) ! calculate approriate constant
      end do
      call ls_cancellakerecharge()
      return
c
      end subroutine
c
c ------------------------------------------------------------------------------------
c
      REAL function rflsnearflow (cdum)
c
c ------------------------------------------------------------------------------------
c
c     Function looks for nearest line-sink end point and returns the
c     streamflow=overlandflow+baseflow
c     If no line-sink is found zero is returned.
c
      implicit none
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
      INTEGER ils,ilsclose
      REAL rdist,rdistmin
      COMPLEX cdum
      rflsnearflow=0.0d0
      ilsclose=0
      if (nls.gt.0) then
       rdistmin=1.0d+30
       do ils=1,nls
        rdist=ABS(cdum-clsze(ils))
        if (rdist.lt.rdistmin) then
         rdistmin=rdist
         ilsclose=ils
        end if
       end do
      end if
      if (ilsclose.gt.0) then
      rflsnearflow=rlsbf(ilsclose)+rlsof(ilsclose)
      endif
      return
c
      end
c
c ------------------------------------------------------------------------------------
c
      REAL function rflsnearlake (cdum)
c
c ------------------------------------------------------------------------------------
c
c     Function looks for nearest line-sink (of a lake) and returns the specified head,
c     which has been solved for using two estimated lake stages.
c     The user is responsible for placing the "lake stage" near a line-sink boundary of a lake feature.
c     If no line-sink is found zero is returned.
c     NOTE: the distance to the nearest line-sink is returned, whether a lake or not!
c
      implicit none
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
      INTEGER ils,ilsclose
      REAL rdist,rdistmin
      COMPLEX cdum
      rflsnearlake=0.0d0
      if (nls.gt.0) then
       rdistmin=1.0d+30
       do ils=1,nls
        rdist=ABS(cdum-clsze(ils))
        if (rdist.lt.rdistmin) then
         rdistmin=rdist
         ilsclose=ils
        end if
       end do
      end if
      rflsnearlake=rlsh(ilsclose)
      return
c
      end
c

c -----------------------------------------------------------------------------------------------
c
      subroutine ls_setwidth ()
c
c -----------------------------------------------------------------------------------------------
c
c     Routine scans all head specified line-sinks with a resistance set larger than zero and
c     calculates the appropriate line-sink width (see "Dealing with resistance to flow into surface waters.pdf"
c     The process is controlled by the code ilsbound(ils), which is an integer added to the "width" command in the
c     line-sink input module. The actions are as follows:
c
c     ilsbound=0    accept the current width as the effective leakage zone (no action)
c     ilsbound=1    calculate the effective leakage zone assuming the line-sink is along the center line of a stream
c     ilsbound=2    calculate the effective leakage zone assuming the line-sink is along the boundary of a stream
c
c     Note: in the case of ilsbound=0 it is assumed that the width parameter provided by the GUI does
c           already reflect the effective leakage zone. This is true for data sets produced before version 2.1.0
c           of the GUI (before April 2005). After calculating and storing the proper leakage zone into RLSWID(ils)
c           the code ilsbound is set to zero to ensure that the effective leakage zone is not being recalculated
c           in case ls_setwidth is called again.
c
c     This routine is called in LSPREP, which is called in SOLUT.
c
c
      implicit none
      INTEGER ils,i,iflag
      REAL rfbase,rfperm,raquifertop,raquiferbase,rhlocal,
     &     rklocal,rtransmissivity,rlambda,rwidth,rzone
      COMPLEX cz
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
      if (nlsh.gt.0) then
       do i=1,nlsh
        ils=klspth(i)
c        write (iluer,1001) ils,ilsbound(ils),rlsres(ils),rlswid(ils),
c     &                     rlsdep(ils)
c 1001 format (' ls_setwidth1: ils,ilsbound,rlsres,rlswid,',
c     &        ' rlsdep ',/,2(I4),3(d14.7))
        if (rlsres(ils).gt.0.0d0.and.ilsbound(ils).gt.0) then  ! line-sink has a resistance, calculate correct leakage zone
          cz=0.5d0*(clszs(ils)+clsze(ils))
c          write (iluer,1002) cz
c 1002 format (' ls_setwidth2: cz ',2(d14.7))
          raquifertop=rlsh(ils)-rlsdep(ils)
          raquiferbase=rfbase(cz)
          rhlocal=raquifertop-raquiferbase
          rklocal=rfperm(cz)
          rtransmissivity=rhlocal*rklocal
          rlambda=SQRT(rlsres(ils)*rtransmissivity)
          rwidth=0.5d0*rlswid(ils) ! this is half the width of the stream or lake
          rzone=rlambda*TANH(rwidth/rlambda)
          if (ilsbound(ils).eq.1) then
            rzone=2.0d0*rzone   ! combine leakage zones on both sides of stream or lake into one.
          end if
          rlswid(ils)=rzone ! store the effective leakage zone for use during the solution process
          ilsbound(ils)=0   ! flag the line-sink as effective leakage zone already stored in rlswidth(ils)
c        write (iluer,1003) raquifertop,raquiferbase,rtransmissivity,
c     &                     rklocal,rlambda,rwidth, rzone
c 1003 format (' ls_setwidth3: raquifertop,raquiferbase,rtransmissivity,'
c     &       ,'rklocal,rlambda,rwidth, rzone',/,6(d14.7))
        end if
       end do
      end if
      end subroutine
c
c ------------------------------------------------------------------------
c
      subroutine ls_cancellakerecharge()
c
c ------------------------------------------------------------------------
c
c     Routine adds a recharge only inhomogeneity for each lake feature it finds.
c     The recharge will be set equal to the negative recharge found inside the lake.
c     The inhomogeneity is created by the routine "db_cancellakerecharge"
c
      implicit none
      INTEGER i,istr,iad,inext,isize
      CHARACTER(16) alabtemp,alabfirst
      COMPLEX(8) cztemp
      parameter (isize=500)
      DIMENSION cztemp(isize),alabtemp(isize)
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
c
      do istr=1,nlstring
       iad=klstring(istr)
       i=0
       if (lslake(iad)) THEN ! found lake
        alabfirst=alslab(iad)
  5     i=i+1
        if (i.le.isize) then ! OK space in temporary arrays
         cztemp(i)=clszs(iad)
         alabtemp(i)=alslab(iad)
         inext=klsdn(iad)
         if (inext.gt.0) then ! next line-sink
          iad=inext
          GOTO 5
         else                 ! at end of string, create inhomogeneity
          call db_cancellakerecharge(i,cztemp,alabtemp)
         end if
        else                  ! exceeding temporary array size
         write (iluer,1000) isize,alabfirst
        endif
       endif
      enddo
      return
 1000 format (' ***ERROR in ls_cancellakerecharge: ',
     &'Maximum array size of ',i5,' exceeded.',/,
     &' Lake with first label ',a16,' does not have recharge canceled!')
      end subroutine

