C     Last change:  HMH  20 Oct 2017   11:22 am
c      This file contains the following routines and functions:
c
c      RFTOP      returns elevation of the confining layer at CZ
c      RFHEDP     returns the HEAD belonging to the potential RPOT
C      rfhfp      calculates the head from the potential
c      GVPAR      returns given commonblock data as arguments
c      RFBASE     returns the aquifer base at CZ
c      RFPOTH     returns the POTENTIAL belonging to the head RHEDIN
c      RFINTERFACE     returns the interface elevation
C      LFINTERFACE     .TRUE. is interface flow being solved for
C      INTERFACEDATA   returns interface parameters
c      RFPERM     returns the hydraulic conductivity at CZ
c      RFPOR      returns the porosity at CZ
c      RFHGHT     returns the SATURATED AQUIFER THICKNESS at CZ
c      RGVREFDIST returns the distance from CZ to the reference point
C
C ---------------------------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFTOP (CZ)
C
C ---------------------------------------------------------------------------------------------------
C
C     Function returns elevation of the confining layer measured
C     with respect to the datum. The current version of GFLOW
C     has a constant elevation of the aquifer top (CZ is not used).
C
      IMPLICIT NONE
      COMPLEX(8) CZ
      INCLUDE 'gvcom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc' 
      RFTOP=RBASE+RH
      RETURN
      END
C
C ---------------------------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFHEDP (RPOT,CZ)
C
C ---------------------------------------------------------------------------------------------------
C
C     Function returns the HEAD belonging to the potential RPOT.
c     function calls real function rfhfp to conduct the actual conversion
C
      IMPLICIT NONE
      REAL(8) RPOT,RFHDP0,RLPOT0,RKLOC,RFPERM,RHLOC,RFBASE,RFTOP,RFHFP,
     &        RBLOC
      COMPLEX(8) CZ,CZ0
      INCLUDE 'gvcom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      DATA RLPOT0,RFHDP0,CZ0 /2*1.E21,(1.E21,1.E21)/
      RFHEDP=RFHDP0
      IF (RPOT.EQ.RLPOT0.AND.CZ.EQ.CZ0) RETURN
      RKLOC=RFPERM(CZ)
      RBLOC=RFBASE(CZ)
      RHLOC=RFTOP(CZ)-RBLOC
      RFHEDP=rfhfp(rpot,rkloc,rhloc,rbloc)+RBLOC
      RFHDP0=RFHEDP
      RLPOT0=RPOT
      CZ0=CZ  
      RETURN
      END
c
c --------------------------------------------------------------------------------------------------------
c
      REAL(8) function rfhfp (rpot,rkloc,rhloc,rbloc)
c
c     Function performs the actual conversion of a potential to a head
c     the head is measured with respect to the local aquifer base.
c     This function is called by RFHEDP, DBKFAC and RFDBER.
C     INPUT:
c            RPOT = local discharge potential
c            RKLOC= local hydraulic conductivity
c            RHLOC= local total aquifer thickness (saturated or not)
c            RBLOC= local aquifer base elevation
c
      implicit none
      REAL(8) rpot,rkloc,rhloc,rcond,rbloc,
     &        RHS,RALPHA,RBETA,RHED,RCCONFINED,RCUNCONFINED,
     &        RTEST1,RTEST2
      INCLUDE 'gvcom.inc'
      INCLUDE 'lusys.inc'
        RCOND=0.5D0*RKLOC*RHLOC*RHLOC
c        write (iluer,1001) rpot,rkloc,rhloc,rbloc
c 1001 format (' rfhfp1: rpot,rkloc,rhloc,rbloc ',4(d14.7))
      if (linterface) THEN ! ---------------- interface flow ------------
        RHS=RSEALEVEL-RBLOC
        RTEST1=RKLOC*RGVFAC1*RHLOC*RHS-RCOND
        RTEST2=0.5D0*RKLOC*(RGVFAC1*RHS)*(RGVFAC1*RHS)
c      write (iluer,1002) rcond,rtest1,rtest2
c 1002 format (' rfhfp2: rcond,rtest1,rtest2 ',3(d14.7))
       if (rpot.ge.rcond) then ! ---------------- confined flow -------------------
          RCCONFINED=RKLOC*RGVFAC1*RHLOC*RHS-
     &               0.5D0*RKLOC*RGVFAC1*RHLOC*RHLOC
          RALPHA=RGVFAC2/RGVFAC1
          RBETA=-RGVFAC2*RHS+RHLOC
        IF (RPOT.GE.RTEST1) THEN ! confined no interface
          rfhfp=(RPOT+RCOND)/(RKLOC*RHLOC)
c      write (iluer,1003) rfhfp
c 1003 format (' rfhfp3 confined, no interface: rfhfp=',d14.7 )
        ELSE IF (RPOT.GT.RCCONFINED) THEN ! confined with interface
         rfhfp=SQRT(2.0D0*(RPOT-RCCONFINED)/(RKLOC*RALPHA))-RBETA/RALPHA
c      write (iluer,1004) ralpha,rbeta,rcconfined,rfhfp
c 1004 format(' rfhfp4 confined+interface: ralpha,rbeta,rcconfined,',
c     &       ' rfhfp=',/,5(d14.7))
        ELSE ! head lower than head at the coast
         rfhfp=-RBETA/RALPHA ! set head at the coast
        endif
       else    ! --------------- unconfined flow ----------------------
          RCUNCONFINED=0.5D0*RKLOC*RGVFAC1*RHS*RHS
        IF (RPOT.GE.RTEST2) THEN ! unconfined no interface
          rfhfp=SQRT(2.0*RPOT/RKLOC)
c      write (iluer,1005) rfhfp
c 1005 format (' rfhfp5 unconfined, no interface: rfhfp=',d14.7 )
        ELSEIF (RPOT.GT.RCUNCONFINED) THEN ! unconfined with interface
          RALPHA=RGVFAC2
          RBETA=-RGVFAC2*RHS
          rfhfp=SQRT(2.0D0*(RPOT-RCUNCONFINED)/(RKLOC*RALPHA))
     &         -RBETA/RALPHA
c      write (iluer,1006) ralpha,rbeta,rcunconfined,rfhfp
c 1006 format(' rfhfp6 unconfined+interface ralpha,rbeta,rcunconfined,',
c     &       ' rfhfp=',/,5(d14.7))
        ELSE ! head is below sea level, no fresh water present
           rfhfp=RHS ! set head at sea level
        ENDIF
       endif
      else  ! -------------------  original logic without interface --------
        if (rpot.le.0.0) then        ! head at aquifer base
          rfhfp=0.0
        else if (rpot.lt.rcond) then  ! head below aquifer top: unconfined conditions
          rfhfp=SQRT(2.0*RPOT/RKLOC)
        else                     ! head at or above aquifer top: confined conditions
          rfhfp=(RPOT+RCOND)/(RKLOC*RHLOC)
        endif
      endif
      return
      end
C
C ---------------------------------------------------------------------------------------------------
C
      SUBROUTINE GVPAR (RK0,RH0,RHED0,RB0,RP0,CZ0)
C
C ---------------------------------------------------------------------------------------------------
C
C     Routine returns the data from the given commonblock, see GVCOM.INC
C
      IMPLICIT NONE
      REAL(8) RK0,RH0,RHED0,RB0,RP0
      COMPLEX(8) CZ0
      INCLUDE 'gvcom.inc'
C
      RK0=RK
      RH0=RH
      RHED0=RHEAD0
      RB0=RBASE
      RP0=RPOR
      CZ0=CREFZ
      RETURN
      END SUBROUTINE
C
C ---------------------------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFBASE (CZ)
C
C ---------------------------------------------------------------------------------------------------
C
C     Function returns the elevation of the aquifer base.
C
      IMPLICIT NONE
      LOGICAL Linside_inhom
      REAL(8) RT,RBO,RBI,RP,RKO,RKI
      COMPLEX(8) CZ
      INCLUDE 'gvcom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      IF (Linside_inhom (CZ,RKI,RP,RBI,RT)) THEN
c      write (iluer,1001)
 1001 format (' rfbase1: inside')
      RFBASE=RBI
      ELSE
c      write (iluer,1002)
 1002 format (' rfbase2: outside')
      RFBASE=RBASE
      ENDIF
c      write (iluer,1003) cz,rfbase
 1003 format (' rfbase3: rfbase=',3(E14.7))
      RETURN
      END
C
C ---------------------------------------------------------------------------------------------------
C
      REAL FUNCTION RFPOTH (RHEDIN,CZ)
C
C ---------------------------------------------------------------------------------------------------
C
C     Function returns the POTENTIAL belonging to the head RHED.
C
      IMPLICIT NONE
      REAL(8) RHEDIN,RFTOP,RFPOTH_confined,RFPOTH_unconfined
      COMPLEX(8) CZ
      INCLUDE 'gvcom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
C
      IF (RHEDIN.GE.RFTOP(CZ)) THEN ! confined flow
      RFPOTH=RFPOTH_confined (RHEDIN,CZ)
      ELSE ! unconfined flow
      RFPOTH=RFPOTH_unconfined (RHEDIN,CZ)
      ENDIF
      RETURN
      END
C
C ---------------------------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFPOTH_confined (RHEDIN,CZ)
C
C ---------------------------------------------------------------------------------------------------
C
C     Function returns the POTENTIAL belonging to the head RHED under confined conditions.
C
      IMPLICIT NONE
      REAL(8) RFPERM,RFTOP,RFBASE,RHEDIN,RHEDTIP,RALPHA,RBETA,RCOND,
     &        RKLOC,RHLOC,RBLOC,RHED,RHS,RCCONFINED,RCUNCONFINED,RDUM
      COMPLEX(8) CZ
      INCLUDE 'gvcom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      RFPOTH_confined=0.0D0
      RKLOC=RFPERM(CZ)
      RBLOC=RFBASE(CZ)
      RHLOC=RFTOP(CZ)-RBLOC
      RHED=RHEDIN-RBLOC
      IF (LINTERFACE) THEN ! --------- interface flow -------------
        RHS=RSEALEVEL-RBLOC
        RHEDTIP=RGVFAC1*RHS
c      write (ilume,1001) rhed,rhs,rgvfac1,rgvfac2,rhedtip
c 1001 FORMAT(' rfpoth1: rhed,rhs,rgvfac1,rgvfac2,rhedtip ',5(d14.7))
          IF (RHED.GE.RHEDTIP) THEN ! no interface
            RCOND=0.5D0*RKLOC*RHLOC*RHLOC
            RFPOTH_confined=RKLOC*RHLOC*RHED-RCOND
c      write (ilume,1002) rcond,rfpoth
c 1002 format (' rfpoth2 confined, no interface: rcond,rfpoth=',2(d14.7))
          ELSE ! interface present
            RALPHA=RGVFAC2/RGVFAC1
            RBETA=-RGVFAC2*RHS+RHLOC
            RCCONFINED=RKLOC*RGVFAC1*RHLOC*RHS
     &                 -0.5D0*RKLOC*RGVFAC1*RHLOC*RHLOC
            RDUM=(RHED+RBETA/RALPHA)
            RFPOTH_confined=0.5D0*RKLOC*RALPHA*RDUM*RDUM+RCCONFINED
c      write (ilume,1003) ralpha,rbeta,rcconfined,rdum,rfpoth
c 1003 format (' rfpoth3 confined+interface: ralpha,rbeta,rcconfined,',
c     &        ' rdum,rfpoth=',/,5(d14.7))
          END IF
      ELSE  ! --------------- original logic without interface -------------
          RCOND=0.5D0*RKLOC*RHLOC*RHLOC
          RFPOTH_confined=RKLOC*RHLOC*RHED-RCOND
      ENDIF
      RETURN
      END
C
C ---------------------------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFPOTH_unconfined (RHEDIN,CZ)
C
C ---------------------------------------------------------------------------------------------------
C
C     Function returns the POTENTIAL belonging to the head RHED under unconfined conditions.
C
      IMPLICIT NONE
      REAL(8) RFPERM,RFTOP,RFBASE,RHEDIN,RHEDTIP,RALPHA,RBETA,RCOND,
     &        RKLOC,RHLOC,RBLOC,RHED,RHS,RCCONFINED,RCUNCONFINED,RDUM
      COMPLEX(8) CZ
      INCLUDE 'gvcom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      RFPOTH_unconfined=0.0D0
      RKLOC=RFPERM(CZ)
      RBLOC=RFBASE(CZ)
      RHLOC=RFTOP(CZ)-RBLOC
      RHED=RHEDIN-RBLOC
      IF (LINTERFACE) THEN ! --------- interface flow -------------
        RHS=RSEALEVEL-RBLOC
        RHEDTIP=RGVFAC1*RHS
c      write (ilume,1001) rhed,rhs,rgvfac1,rgvfac2,rhedtip
c 1001 FORMAT(' rfpoth1: rhed,rhs,rgvfac1,rgvfac2,rhedtip ',5(d14.7))
          IF (RHED.GE.RHEDTIP) THEN ! no interface
            RFPOTH_unconfined=0.5D0*RKLOC*RHED*RHED
c      write (ilume,1004) rfpoth
c 1004 format (' rfpoth4 unconfined, no interface: rfpoth=',d14.7)
          ELSE ! interface present
            RALPHA=RGVFAC2
            RBETA=-RGVFAC2*RHS
            RCUNCONFINED=0.5D0*RKLOC*RGVFAC1*RHS*RHS
            RDUM=(RHED+RBETA/RALPHA)
            RFPOTH_unconfined=0.5D0*RKLOC*RALPHA*RDUM*RDUM+RCUNCONFINED
c      write (ilume,1005) ralpha,rbeta,rcunconfined,rdum,rfpoth
c 1005 format(' rfpoth5 unconfined+interface: ralpha,rbeta,rcunconfined,'
c     &        ,' rdum,rfpoth=',/,5(d14.7))
          END IF
      ELSE  ! --------------- original logic without interface -------------
        IF (RHED.LE.0.0D0) RETURN
          RFPOTH_unconfined=0.5D0*RKLOC*RHED*RHED
      ENDIF
      RETURN
      END
C
C ----------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFINTERFACE(CZ)
C
C     Routine returns the elevation of the interface  
      IMPLICIT NONE
      INCLUDE 'gvcom.inc'
      INCLUDE 'lusys.inc'
      REAL(8) RHS,RFBASE,RFHEAD,RBLOC,RHED
      COMPLEX(8) CZ
      RBLOC=RFBASE(CZ)
      rfinterface=rbloc            ! Henk on 9/1/17
      if (.not.linterface) return  ! Henk on 9/1/17
      RHS=RSEALEVEL-RBLOC
      RHED=RFHEAD(CZ)-RBLOC
      RFINTERFACE=RGVFAC2*RHS-RGVFAC2/RGVFAC1*RHED
c     This is the interface elevation relative to the aquifer base, now
c     get the interface elevation relative to the reference plane (same as the heads).
      RFINTERFACE=RFINTERFACE+RBLOC
      END
c
c ------------------------------------------------------------------------------------
c
      LOGICAL function lfinterface()
c
c     True if interface flow is being solved.
c
      IMPLICIT NONE
      INCLUDE 'gvcom.inc'
      INCLUDE 'lusys.inc'
      lfinterface=linterface
      return
      end
c
c ------------------------------------------------------------------------------------
c
      subroutine interfacedata (rhss,rsgs,rsgf,rfac1,rfac2)
c
c     Returns data pertaing to interface flow
c
c       rhss = local average sea level (RSEALEVEL)
c       rsgs  = salt water specific gravity (RSPECIFICGRAVITYSALT)
c       rsgf  = fresh water specific gravity (RSPECIFICGRAVITYFRESH)
c       rfac1 = rsgs/rsgf (RGVFAC1)
c       rfac2 = rsgs/(rsgs-rsgf) (RGVFAC2)
      IMPLICIT NONE
      INCLUDE 'gvcom.inc'
      INCLUDE 'lusys.inc'
      REAL(8) rhss,rsgs,rsgf,rfac1,rfac2
      RHSS=RSEALEVEL
      RSGS=RSPECIFICGRAVITYSALT
      RSGF=RSPECIFICGRAVITYFRESH
      RFAC1=RGVFAC1
      RFAC2=RGVFAC2
      return
      end subroutine
C
C ---------------------------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFPERM (CZ)
C
C ---------------------------------------------------------------------------------------------------
C
C     Function returns the HYDRAULIC CONDUCTIVITY at CZ.
C
      IMPLICIT NONE
      INCLUDE 'gvcom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      LOGICAL Linside_inhom
      REAL(8) RKI,RKO,RP,RBI,RBO,RT
      COMPLEX(8) CZ
      IF (Linside_inhom (CZ,RKI,RP,RBI,RT)) THEN
      RFPERM=RKI
      ELSE
      RFPERM=RK
      ENDIF
c      if (abs(cz).lt.10.0) rfperm=1.0 ! temporary statement for circ. inh.
      RETURN
      END
C
C ---------------------------------------------------------------------------------------------------
C
      REAL FUNCTION RFPOR (CZ)
C
C ---------------------------------------------------------------------------------------------------
C
C     Function returns the POROSITY at CZ.
C
      IMPLICIT NONE
      INCLUDE 'gvcom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      LOGICAL Linside_inhom
      REAL(8) RKI,RKO,RP,RBI,RBO,RT
      COMPLEX(8) CZ
      IF (Linside_inhom (CZ,RKI,RP,RBI,RT)) THEN
      RFPOR=RP
      ELSE
      RFPOR=RPOR
      ENDIF
      RETURN
      END
C
C ---------------------------------------------------------------------------------------------------
C
      REAL FUNCTION RFHGHT (CZ)
C
C ---------------------------------------------------------------------------------------------------
C
C     Function returns the SATURATED AQUIFER THICKNESS at CZ.
C     In case of interface flow this is the distance between the
C     fresh-salt water interface and the saturated aquifer top.
C
      IMPLICIT NONE
      INCLUDE 'gvcom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc' 
      REAL(8) RHEAD,RTOP,RFHEAD,RFTOP,RFBASE,RBLOC,RFINTERFACE
      COMPLEX(8) CZ
      RHEAD=RFHEAD(CZ)
      RTOP=RFTOP(CZ)
      RBLOC=RFBASE(CZ)
      if (linterface) then
       rbloc=MAX(rbloc,rfinterface(cz))
      end if
      IF (RHEAD.LT.RTOP) THEN
       RFHGHT=RHEAD-RBLOC   ! unconfined
      ELSE
       RFHGHT=RTOP-RBLOC      ! confined
      ENDIF
      RETURN
      END
C
C ---------------------------------------------------------------------------------------------------
C
       REAL FUNCTION RGVREFDIST (CZ)
C
C ---------------------------------------------------------------------------------------------------
C
      IMPLICIT NONE
      COMPLEX CZ
      INCLUDE 'gvcom.inc'
      RGVREFDIST=ABS(CREFZ-CZ)
      rgvrefdist=MAX(1.0,rgvrefdist)
      RETURN
      END
