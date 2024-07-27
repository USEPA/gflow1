C     Last change:  HMH   1 Sep 2017    2:08 pm
c --------------------------------------------------------------------------------
c
c     This file contains the following routines or functions
c
c     DBKFAC generates the jump conditions for the transmissivity
c     DBMATSIZE updates the matrix array dimensions
c     DBCZC  generates control points and sets equation type
c     DBMAT  generates matrix equations
c     DBKNO  generates the known vector
c     DBSUB  substitutes the solution vector in the common blocks
c     RFINTLABLEFT  integral for linear strength along line doublet
c     RFINTLABRIGHT integral for linear strength along line doublet
c     RFINTLABMU    integral for quadratic strength along line doublet
c     RFDBNFCONDITION generates normal flow condition for matrix coefficients
c     RFDBCONDUCTANCE calculates the conductance for a slurry wall
c     ldbbottomOFF    converts bottom jumps to k jumps for stability reasons (temporarily)
c     ldbbottomON     restores base jumps
c     dbupdate        update the known vector for the line-doublets after an iteration using the matrix
c     dbupdate_check  debugging function to check the working of dbupdate
c     dbgentbfac      called by DBKFAC to set RDBTFAC... parameters
c       
c
C ------------------------------------------------------------------------------------------------------
c
      SUBROUTINE DBKFAC
c
c ------------------------------------------------------------------------------------------------------
C
C     Routine generates:
c     For inhomogeneity domains   the factor Ti/(Ti-To) used in DBMAT and DBKNO
C                                 and the factor TiTo/(Ti-To)*(bi-bo)/2 used in DBKNO
c     For slurry walls            the conductance at the first node and the center of a line doublet
C
      IMPLICIT NONE
      INTEGER(4) ISTR,INOD1,INODL,INOD,INODP1
      INTEGER(4) InterfaceCase
      LOGICAL lfinterface,lshift1,lshift2
      REAL(8)RBO,RBI,RKO,RKI,rhcol,rfcol,rgcol,
     &        RC1,RC2,RFDBCONDUCTANCE
      REAL(8) rzInterfaceElevation,rfinterface,rCInterfaceI,
     &        rCInterfaceO,rhss,rsgs,rsgf,rfac1,rfac2
      COMPLEX(8) CZ,CZA,CZB,cfsmallz
      INCLUDE 'dbcom.inc'
      INCLUDE 'main.inc'
      INCLUDE 'lusys.inc'
C
      IF (NDB.EQ.0) RETURN
      IF (NDBSTR.EQ.0) THEN
      WRITE (ILUER,1000) NDB
      RETURN
      ENDIF
      DO ISTR=1,NDBSTR
      INOD1=IDBSTA(ISTR)
      INODL=INOD1+NDBSTI(ISTR)-1
      IF (IDOMAINTYPE(ISTR).EQ.1.and..not.ldbrechargeonly(istr)) then !    ---------- inhomogeneity domains ----------------
          DO INOD=INOD1,INODL
          IF (IDBCODE(INOD).GT.0) THEN  ! only when we need to make equations in the matrix
C
           RKI=RDBKI(INOD)   ! get inside and outside properties for line doublet starting at CDBZ(INOD)
           RKO=RDBKO(INOD)
           RBI=RDBBI(INOD)
           RBO=RDBBO(INOD)
c
           lshift1=idbcode(inod).eq.2
           inodp1=inod+1
           if (inod.eq.inodl) inodp1=inod1
           lshift2=idbcode(inodp1).eq.2.or.idbcode(inodp1).eq.0
c
c           --------------- first collocation point on line-doublet INOD
c
            cz=cdbz(inod)
            if (lshift1) then
             cz=CMPLX(-1.0d0+rdboffsetinhom,0.0d0) ! make BIGZ
             cz=cfsmallz(cz,cdbz(inod),cdbz(inodp1)) ! convert to z in physical plane
            end if
            rhcol=rdbhnod(inod) ! set in DBUPDATE
            call dbgentbfac (cz,rhcol,rki,rko,rbi,rbo,rfcol,rgcol)
            rdbtfacnod(inod)=rfcol
            rdbbfacnod(inod)=rgcol
c
c ----------------------- center collocation point on line-doublet INOD
c
            cz=0.5d0*(cdbz(inod)+cdbz(inodp1))
            rhcol=rdbhctr(inod) ! set in DBUPDATE
            call dbgentbfac (cz,rhcol,rki,rko,rbi,rbo,rfcol,rgcol)
            rdbtfacctr(inod)=rfcol
            rdbbfacctr(inod)=rgcol
c
c ---------------------- last collocation point on line-doublet INOD (if exists)
c
          if (lshift2) then
            cz=CMPLX(1.0d0-rdboffsetinhom,0.0d0) ! make BIGZ
            cz=cfsmallz(cz,cdbz(inod),cdbz(inodp1)) ! convert to z in physical plane
            rhcol=rdbhend(inod) ! set in DBUPDATE
            call dbgentbfac (cz,rhcol,rki,rko,rbi,rbo,rfcol,rgcol)
            rdbtfacend(inod)=rfcol
            rdbbfacend(inod)=rgcol
          end if
c
          ENDIF
          END DO
      ENDIF          ! -------------------- end inhomogeneity domains ----------------------
      IF (IDOMAINTYPE(ISTR).EQ.2) THEN
c
c      Open slurry wall
c
        DO INOD=INOD1,INODL
! first interval on center of line doublet
          IF (INOD.EQ.INOD1) THEN
            CZA=CDBZ(INOD)+0.001D0*(CDBZ(INOD+1)-CDBZ(INOD)) ! first point near start of string
          ELSE
            CZA=CZB
          END IF
          RC1=RFDBCONDUCTANCE(ISTR,CZA)
          IF (INOD.EQ.INODL) THEN
            CZB=CDBZ(INOD)+0.999D0*(CDBZ(INOD+1)-CDBZ(INOD)) ! last point near end of string
            RC2=RFDBCONDUCTANCE(ISTR,CZB)
            RDBTFACCTR(INOD)=0.5D0*(RC1+RC2) ! average conductance for interval on center of line doublet
           ELSE
            CZB=CDBZ(INOD)+0.75D0*(CDBZ(INOD+1)-CDBZ(INOD))
            RC2=RFDBCONDUCTANCE(ISTR,CZB)
            RDBTFACCTR(INOD)=0.5D0*(RC1+RC2) ! average conductance for interval on center of line doublet
 ! second interval straddles the line doublet end point
            CZB=CDBZ(INOD+1)+0.25D0*(CDBZ(INOD+2)-CDBZ(INOD+1))
 ! average conductance for interval that straddles vertex INOD+1
            RDBTFACNOD(INOD+1)=0.5D0*(RC2+RFDBCONDUCTANCE(ISTR,CZB))  ! not done for INOD=INODL
           END IF
        END DO
      END IF
      IF (IDOMAINTYPE(ISTR).EQ.3) THEN
c
c      Closed slurry wall
c
        CZA=CDBZ(INODL)+(1.0D0-RDBOFFSETN)*(CDBZ(INOD1)-CDBZ(INODL)) ! first interval straddles first node
        CZB=CDBZ(INOD1)+RDBOFFSETN*(CDBZ(INOD1+1)-CDBZ(INOD1))
        RC1=RFDBCONDUCTANCE(ISTR,CZA)
        RC2=RFDBCONDUCTANCE(ISTR,CZB)
        RDBTFACNOD(INOD1)=0.5D0*(RC1+RC2)  ! average conductance
        DO INOD=INOD1,INODL-2
        CZA=CDBZ(INOD)+RDBOFFSETC*(CDBZ(INOD+1)-CDBZ(INOD))  ! next interval on center of line doublet
        CZB=CDBZ(INOD)+(1.0D0-RDBOFFSETC)*(CDBZ(INOD+1)-CDBZ(INOD))
        RC1=RFDBCONDUCTANCE(ISTR,CZA)
        RC2=RFDBCONDUCTANCE(ISTR,CZB)
        RDBTFACCTR(INOD)=0.5D0*(RC1+RC2)   ! average conductance
        CZA=CDBZ(INOD)+(1.0D0-RDBOFFSETN)*(CDBZ(INOD+1)-CDBZ(INOD))
        CZB=CDBZ(INOD+1)+RDBOFFSETN*(CDBZ(INOD+2)-CDBZ(INOD+1)) ! next interval straddles end node of line doublet
        RC1=RFDBCONDUCTANCE(ISTR,CZA)
        RC2=RFDBCONDUCTANCE(ISTR,CZB)
        RDBTFACNOD(INOD+1)=0.5D0*(RC1+RC2)  ! average conductance
        END DO
        CZA=CDBZ(INODL-1)+RDBOFFSETC*(CDBZ(INODL)-CDBZ(INODL-1))  ! next interval on center of one but last line doublet
        CZB=CDBZ(INODL-1)+(1.0D0-RDBOFFSETC)*
     &           (CDBZ(INODL)-CDBZ(INODL-1))
        RC1=RFDBCONDUCTANCE(ISTR,CZA)
        RC2=RFDBCONDUCTANCE(ISTR,CZB)
        RDBTFACCTR(INODL-1)=0.5D0*(RC1+RC2)   ! average conductance
        CZA=CDBZ(INODL-1)+(1.0D0-RDBOFFSETN)*
     &           (CDBZ(INODL)-CDBZ(INODL-1))  ! next interval straddles last node (INODL)
        CZB=CDBZ(INODL)+RDBOFFSETN*(CDBZ(INOD1)-CDBZ(INODL))
        RC1=RFDBCONDUCTANCE(ISTR,CZA)
        RC2=RFDBCONDUCTANCE(ISTR,CZB)
        RDBTFACNOD(INODL)=0.5D0*(RC1+RC2)     ! average conductance
        CZA=CDBZ(INODL)+RDBOFFSETC*(CDBZ(INOD1)-CDBZ(INODL))  ! last interval on center of last line doublet
        CZB=CDBZ(INODL)+(1.0D0-RDBOFFSETC)*(CDBZ(INOD1)-CDBZ(INODL))
        RC1=RFDBCONDUCTANCE(ISTR,CZA)
        RC2=RFDBCONDUCTANCE(ISTR,CZB)
        RDBTFACCTR(INODL)=0.5D0*(RC1+RC2)    ! average conductance
      END IF
      enddo
c      write (ilume,1003) rdbk(1),rdbb(1),rdbk(2),rdbb(2)
c 1003 format (' dbkfac1: rdbk(1),rdbb(1),rdbk(2),rdbb(2) ',4(d14.7))
      RETURN
 1000 FORMAT (' ***ERROR in DBKFAC: ',I4,' doublet nodes, but no ',
     &        'domains encountered!')
 2000 FORMAT (' ***ERROR in DBKFAC: the point ',2F11.1,' does not ',
     &'occur inside string ',I3,/,
     &' Likely cause is a sharp corner in the inhomogeneity boundary.')
 3000 format (' ***Error in DBKFAC: InterfaceCase=0 (not set)')
      END
c
c ---------------------------------------------------------------------------------------------------
c
      subroutine dbmatsize (M,N)
c
c ---------------------------------------------------------------------------------------------------
c
c     updates the matrix array dimensions M and N
c
c     M is the number of equations to be updated
c     N is the number of strength parameters to be updated
c
c     Currently equations are added for collocation points only, thus no overspecification.
c     The routine will contribute to the asymmetry of the matrix when overspecification is
c     implemented. This will require a change in the logic for this routine.
c
c
      implicit none
      INTEGER M,N,inod,inod1,inodl,istr,inodp1,N0
      LOGICAL ls,lc,le,lshift1,lshift2
      INCLUDE 'dbcom.inc'
      INCLUDE 'lusys.inc'
      ndbstrengths=0      ! number of doublet strength parameters used in "dbkeep" routine
      IF (NDB.EQ.0) RETURN
      N0=N
      DO 20 ISTR=1,NDBSTR
      INOD1=IDBSTA(ISTR)
      INODL=INOD1+NDBSTI(ISTR)-1  ! last node of closed string, one but last node of open string
      IF (IDOMAINTYPE(ISTR).EQ.1) THEN
      IF (ldbrechargeonly(istr)) GOTO 20 !  "recharge only", skip domain
c
c      Hydraulic conductivity inhomogeneity
c
        DO INOD=INOD1,INODL           ! we are on line-doublet between INOD and INODP1 of string ISTR
        if (idbcode(inod).gt.0) THEN ! line-doublet is in matrix
         INODP1=INOD+1
         IF (INOD.EQ.INODL) INODP1=INOD1
         lshift1=idbcode(inod).eq.2
         lshift2=idbcode(inodp1).eq.2.or.idbcode(inodp1).eq.0
         ls=rdbtfacnod(inod).ne.1.0d21
         lc=rdbtfacctr(inod).ne.1.0d21
         le=rdbtfacend(inod).ne.1.0d21
         if (ls) then ! add collocation point at or near first node
          M=M+1
          N=N+1
         endif
         if (lc) then ! add collocation point between nodes (center of line-doublet)
          M=M+1
          N=N+1
         endif
         if (lshift2.and.le) then ! add a collocation point at the end of the line-doublet
          M=M+1
          N=N+1
         end if
        endif
        END DO
      ENDIF
      IF (IDOMAINTYPE(ISTR).EQ.2) THEN
c
c      Open slurry wall
c
        DO INOD=INOD1,INODL
        M=M+1
        N=N+1            ! first interval on center of line doublet
        IF (INOD.EQ.INODL) THEN
c       skip
        ELSE
        M=M+1            ! second interval straddles the line doublet end point
        N=N+1
        END IF
        END DO
      END IF
      IF (IDOMAINTYPE(ISTR).EQ.3) THEN
c
c      Closed slurry wall
c
        M=M+1
        N=N+1
        DO INOD=INOD1,INODL-2
        M=M+2
        N=N+2
        END DO
        M=M+3
        N=N+3
      endif
  20  CONTINUE
      ndbstrengths=N-N0 ! number of doublet strength parameters used in "dbkeep" routine
      return
c
      end subroutine
c
c ---------------------------------------------------------------------------------------------------
c
      subroutine dbkeep(drb,j,m)
c
c     add number of doublet strengths to array counter and set corresponding
c     values of drb array equal to 0. The drb array will contain any strength differences between the
c     solution and the actual stored strength parameters, see lskeepsigma and lkkeepstrength routines.
c
      implicit none
      INTEGER i,j,m
      REAL(8) drb
      DIMENSION drb(m)
      INCLUDE 'dbcom.inc'
      if (ndbstrengths.gt.0) then  ! we have line-doublets in the matrix
       do i=1,ndbstrengths
         j=j+1
         drb(j)=0.0d0
       end do
      endif
      return
      end subroutine
c
c ---------------------------------------------------------------------------------------------------
c
c
      SUBROUTINE DBCZC (CZI,M,N,DRFAC,CALPH,ITYPE)
c
c ---------------------------------------------------------------------------------------------------
C
C     Routine generates the control points for all doublet strings
c
C
      IMPLICIT NONE
      LOGICAL lshift1,lshift2,ls,lc,le
      INTEGER(4) M,N,ITYPE,ISTR,INOD1,INODL,INOD,INODP1
      REAL(8) DRFAC
      COMPLEX(8) CZI,CALPH,CZZ,CZ,CFSMALLZ
      INCLUDE 'dbcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'

      DIMENSION CZI(*),DRFAC(4,*),CALPH(*),ITYPE(*)
      IF (NDB.EQ.0) RETURN
      DO 20 ISTR=1,NDBSTR
      INOD1=IDBSTA(ISTR)
      INODL=INOD1+NDBSTI(ISTR)-1  ! last node of closed string, one but last node of open string
      IF (IDOMAINTYPE(ISTR).EQ.1) THEN
      IF (ldbrechargeonly(istr)) GOTO 20 !  "recharge only", skip domain
c
c      Hydraulic conductivity inhomogeneity
c
        DO INOD=INOD1,INODL           ! we are on line-doublet between INOD and INODP1 of string ISTR
        if (idbcode(inod).gt.0) THEN ! line-doublet is in matrix
         INODP1=INOD+1
         IF (INOD.EQ.INODL) INODP1=INOD1
         lshift1=idbcode(inod).eq.2
         lshift2=idbcode(inodp1).eq.2.or.idbcode(inodp1).eq.0
         ls=rdbtfacnod(inod).ne.1.0d21
         lc=rdbtfacctr(inod).ne.1.0d21
         le=rdbtfacend(inod).ne.1.0d21
         if (ls) then ! add first collocation point
          M=M+1 ! add an equation for each collocation point
          N=N+1 ! add an strength parameter for each collocation point
c
c Note: The current implementation uses "collocation points" with one equation and one strength parameter for each.
c       This logic may be modified to use "control points" with one equation for each control point,
c       but the same number of strength parameters as currently in use.
c
          if (lshift1) then
            cz=CMPLX(-1.0d0+rdboffsetinhom,1.0d-8)
            CZI(M)=CFSMALLZ (CZ,CDBZ(INOD),CDBZ(INODP1))  ! first collocation point is offset from line-doublet starting node
          else
            CZI(M)=CDBZ(INOD)  ! first collocation point is line doublet starting point
          end if
          ITYPE(M)=1
          IF (ABS(RDBTFACNOD(INOD)).GT.1.0D4)  ITYPE(M)=6 ! Ti=To, but bi<>bo
          DRFAC(1,M)=1.0D0
          CALPH(M)=(0.0D0,0.0D0)
         endif
         if (lc) then ! add second collocation point
          M=M+1
          N=N+1
          CZZ=CDBZ(INOD)+CDBZ(INODP1)
          CZI(M)=0.5*CZZ    ! second collocation point is line doublet mid point
          ITYPE(M)=1
          IF (ABS(RDBTFACCTR(INOD)).GT.1.0D4) ITYPE(M)=6  ! Ti=To, but bi<>bo
          DRFAC(1,M)=1.0D0
          CALPH(M)=(0.0D0,0.0D0)
         endif
         if (lshift2.and.le) then ! add a third collocation point at the end of the line-doublet
          M=M+1
          N=N+1
          cz=CMPLX(+1.0d0-rdboffsetinhom,1.0d-8)
          CZI(M)=CFSMALLZ (CZ,CDBZ(INOD),CDBZ(INODP1))  ! extra collocation point is offset from line-doublet end node
          ITYPE(M)=1
          IF (ABS(RDBTFACNOD(INOD)).GT.1.0D4)  ITYPE(M)=6 ! Ti=To, but bi<>bo
          DRFAC(1,M)=1.0D0
          CALPH(M)=(0.0D0,0.0D0)
         end if
        endif
        END DO
      ENDIF
      IF (IDOMAINTYPE(ISTR).EQ.2) THEN
c
c      Open slurry wall
c
        DO INOD=INOD1,INODL
        M=M+1
        N=N+1            ! first interval on center of line doublet
        DRFAC(1,M)=1.0D00
        ITYPE(M)=-2
        IF (INOD.EQ.INOD1) THEN
        CZI(M)=CDBZ(INOD)+0.001*(CDBZ(INOD+1)-CDBZ(INOD)) ! first point near start of string
        ELSE
        CZI(M)=CALPH(M-1)
        END IF
        IF (INOD.EQ.INODL) THEN
        CALPH(M)=CDBZ(INOD)+0.999*(CDBZ(INOD+1)-CDBZ(INOD)) ! last point near end of string
        ELSE
        CALPH(M)=CDBZ(INOD)+0.75*(CDBZ(INOD+1)-CDBZ(INOD))
        M=M+1            ! second interval straddles the line doublet end point
        N=N+1
        DRFAC(1,M)=1.0D00
        ITYPE(M)=-2
        CZI(M)=CALPH(M-1)
        CALPH(M)=CDBZ(INOD+1)+0.25*(CDBZ(INOD+2)-CDBZ(INOD+1))
        END IF
        END DO
      END IF
      IF (IDOMAINTYPE(ISTR).EQ.3) THEN
c
c      Closed slurry wall     (rdboffsetn=0.5 and rdboffsetc=0.05)
c                             these are set in Block Data DBDAT
C                             current settings create overlapping integration segments
c
        M=M+1
        N=N+1
        DRFAC(1,M)=1.0D00
        ITYPE(M)=-2
        CZI(M)=CDBZ(INODL)+(1.0-RDBOFFSETN)*(CDBZ(INOD1)-CDBZ(INODL)) ! first interval straddles first node
        CALPH(M)=CDBZ(INOD1)+RDBOFFSETN*(CDBZ(INOD1+1)-CDBZ(INOD1))
        DO INOD=INOD1,INODL-2
        M=M+1
        N=N+1
        DRFAC(1,M)=1.0D00
        ITYPE(M)=-2
        CZI(M)=CDBZ(INOD)+RDBOFFSETC*(CDBZ(INOD+1)-CDBZ(INOD))  ! next interval on center of line doublet
        CALPH(M)=CDBZ(INOD)+(1.0-RDBOFFSETC)*(CDBZ(INOD+1)-CDBZ(INOD))
        M=M+1
        N=N+1
        DRFAC(1,M)=1.0D00
        ITYPE(M)=-2
        CZI(M)=CDBZ(INOD)+(1.0-RDBOFFSETN)*(CDBZ(INOD+1)-CDBZ(INOD))
        CALPH(M)=CDBZ(INOD+1)+RDBOFFSETN*(CDBZ(INOD+2)-CDBZ(INOD+1)) ! next interval straddles end node of line doublet
        END DO
        M=M+1
        N=N+1
        DRFAC(1,M)=1.0D00
        ITYPE(M)=-2
        CZI(M)=CDBZ(INODL-1)+RDBOFFSETC*(CDBZ(INODL)-CDBZ(INODL-1))  ! next interval on center of one but last line doublet
        CALPH(M)=CDBZ(INODL-1)+(1.0-RDBOFFSETC)*
     &           (CDBZ(INODL)-CDBZ(INODL-1))
        M=M+1
        N=N+1
        DRFAC(1,M)=1.0D00
        ITYPE(M)=-2
        CZI(M)=CDBZ(INODL-1)+(1.0-RDBOFFSETN)*
     &           (CDBZ(INODL)-CDBZ(INODL-1))  ! next interval straddles last node (INODL)
        CALPH(M)=CDBZ(INODL)+RDBOFFSETN*(CDBZ(INOD1)-CDBZ(INODL))
        M=M+1
        N=N+1
        DRFAC(1,M)=1.0D00
        ITYPE(M)=-2
        CZI(M)=CDBZ(INODL)+RDBOFFSETC*(CDBZ(INOD1)-CDBZ(INODL))  ! last interval on center of last line doublet
        CALPH(M)=CDBZ(INODL)+(1.0-RDBOFFSETC)*(CDBZ(INOD1)-CDBZ(INODL))
      END IF
  20  CONTINUE
      RETURN
      END
c
C ---------------------------------------------------------------------------------------------------------------
c
      SUBROUTINE DBGENMAT (DRA,CZI,M,N,J0,DRFAC,CALPH,ITYPE)        ! 7/7/03 for T inhomogeneities only
c
c ---------------------------------------------------------------------------------------------------------------
C
C     Routine generates the initial matrix coefficients for the line doublets
C     Two columns are generated at the same time for DBSTR and DBQSTR
C     of each line doublet.
c     Note: for inhomogeneities with common boundaries 3 columns may be generated per line-doublet,
c     for DBSTR, DBQSTR, and DBSTE
c
c     These initial coeeficients create equations that either generate the potential or tottal flow
c     across a line segment.
c
C
      IMPLICIT NONE
      INTEGER(4) M,N,J0,ITYPE,ISTR,INOD1,INODL,INODM1,INOD,I,J,
     &           IEQS,IEQ,INODP1,jlook
      LOGICAL LNEG,l1,l2,l3,lout,ldominsidedom,lshift1,lshift2,ls,lc,le
      REAL(8) DRA,DRFAC,RCZTEST1,RXA1,RCZTEST2,RXA2,RCZATEST2,RCZATEST1,
     &        RXB2,RXB1,RCORRECTION,RFINTMU,RFINTLABLEFT,RFINTLABRIGHT,
     &        RFPERM
      COMPLEX(8) CZI,CALPH,CFDBF,CFDBG,CFDBS,CFDBFD,CFDBGD,CFDBSD,CZ,
     &           CZA,CDUM
      INCLUDE 'dbcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      DIMENSION DRA(M,*),CZI(*),DRFAC(4,*),CALPH(*),ITYPE(*)
c
      IF (NDB.EQ.0) RETURN
      DO 30 ISTR=1,NDBSTR
      INOD1=IDBSTA(ISTR)
      IF (ldbrechargeonly(istr)) GOTO 30 ! k inside = k outside or "recharge only", skip domain
      INODL=INOD1+NDBSTI(ISTR)-1
      INODM1=INODL
      IF (IDOMAINTYPE(ISTR).EQ.1) THEN
C
C     Hydraulic conductivity inhomogeneity
C
      DO INOD=INOD1,INODL
      if (idbcode(inod).gt.0) THEN ! include line-doublet in the matrix
      inodm1=inod-1
      if (inod.eq.inod1) inodm1=inodl
      inodp1=inod+1
      if (inod.eq.inodl) inodp1=inod1
      lshift1=idbcode(inod).eq.2   ! shift first collocation point away from first node
      lshift2=idbcode(inodp1).eq.2.or.idbcode(inodp1).eq.0 ! add an extra collocation point shated away from the last node
      ls=rdbtfacnod(inod).ne.1.0d21
      lc=rdbtfacctr(inod).ne.1.0d21
      le=rdbtfacend(inod).ne.1.0d21
      DO I=1,M
      J=J0
      CZ=CZI(I)
      CZA=CALPH(I)
      IEQS=ITYPE(I)
C
C ITYPE=+1 potential specified at CZ (Note: at doublet nodes ITYPE=+1,
C          but potential is expressed in terms of local strenght parameters.)
C       -1 difference in potential specified: PHI(CZ)-PHI(CZA) !! to be done in DBMATCORRECT
C       +2 stream function specified at CZ
C       -2 flow across line between CZ and CZA, positive to the left when at CZ
C          note: because PSI is continuous for line doublets, this is the same as
C          the difference bwteen CZ and CZA: PSI(CZ)-PSI(CZA)
C       +3 discharge component normal to the unit vector CZA (rotated to the left)
C       +4 discharge component parallel to the unit vector CZA
C        5 continuity equation: provide total discharge
C        6 generate a zero matrix coefficient. Used to create an element specific equation. !! in DBMATCORRECT
C
      LNEG=IEQS.LT.0
      IEQ=IABS(IEQS)
      CALL DBFILL (INOD1,INODM1,INODL,CZ)  ! next three calls to calculate parameters for line-doublet
      CALL DBFILL (INOD1,INOD,INODL,CZ)
      call dbfill (inod1,inodp1,inodl,cz)
      GOTO (1,2,3,4,5,1),IEQ            ! itype=6 send to itype is 1
C
   1  if (ls) then
       CDUM=CFDBG(INODM1)+CFDBF(INOD) ! provide PHI influence of DBSTR(inod) at CZ
       if (lshift1) cdum=cfdbf(inod)  ! if shifted F function is the only influence of DBSTR(inod)
       J=J+1
       DRA(I,J)=REAL(CDUM)
      endif
      if (lc) then
       J=J+1
       CDUM=CFDBS(INOD)              ! provide PHI influence of DBQSTR(inod) at CZ
       DRA(I,J)=REAL(CDUM)
      endif
      if (lshift2.and.le) THEN ! also provide PHI influence of DBSEN(inod)
        cdum=cfdbg(inod)
        J=J+1
        DRA(I,J)=REAL(CDUM)
      end if
      GOTO 9
  2   if (ls) then
       CDUM=CFDBG(INODM1)+CFDBF(INOD) ! provide PSI influence of DBSTR(inod) at CZ
       if (lshift1) cdum=cfdbf(inod)
       J=J+1
       DRA(I,J)=AIMAG(CDUM)
      endif
      if (lc) then
       J=J+1
       CDUM=CFDBS(INOD)
       DRA(I,J)=AIMAG(CDUM)
      endif
      if (lshift2.and.le) then
        cdum=cfdbg(inod)
        J=J+1
        DRA(I,J)=AIMAG(CDUM)
      end if
      IF (LNEG) THEN                 ! subtract PSI at CZA
       j=j0
       CALL DBFILL (INOD1,INODM1,INODL,CZA)
       CALL DBFILL (INOD1,INOD,INODL,CZA)
       call dbfill (inod1,inodp1,inodl,cza)
       if (ls) then
        CDUM=CFDBG(INODM1)+CFDBF(INOD)
        if (lshift1) cdum=cfdbf(inod)  ! if shifted F function is the only influence of DBSTR(inod)
        j=j+1
        DRA(I,J)=DRA(I,J)-AIMAG(CDUM)
       endif
       if (lc) then
        J=J+1
        CDUM=CFDBS(INOD)
        DRA(I,J)=DRA(I,J)-AIMAG(CDUM)
       endif
       if (lshift2.and.le) then
         cdum=cfdbg(inod)
         J=J+1
         DRA(I,J)=DRA(i,j)-AIMAG(CDUM)
       end if
      ENDIF
      GOTO 9
  3   if (ls) then
       CDUM=CFDBGD(INODM1)+CFDBFD(INOD)  ! provide Q normal to unit vector CZA
       if (lshift1) cdum=cfdbfd(inod)
       J=J+1
       DRA(I,J)=AIMAG(CDUM*CONJG(CZA))
      endif
      if (lc) then
       CDUM=CFDBSD(INOD)
       J=J+1
       DRA(I,J)=AIMAG(CDUM*CONJG(CZA))
      endif
      if (lshift2.and.le) then
        cdum=cfdbgd(inod)
        J=J+1
        DRA(I,J)=AIMAG(CDUM*CONJG(cza))
      end if
      GOTO 9
  4   if (ls) then
       CDUM=CFDBGD(INODM1)+CFDBFD(INOD)  ! provide Q parallel to unit vector CZA
       if (lshift1) cdum=cfdbfd(inod)
       J=J+1
       DRA(I,J)=REAL(CDUM*CONJG(CZA))
      endif
      if (lc) then
       CDUM=CFDBSD(INOD)
       J=J+1
       DRA(I,J)=REAL(CDUM*CONJG(CZA))
      endif
      if (lshift2.and.le) then
        cdum=cfdbgd(inod)
        J=J+1
        DRA(I,J)=real(CDUM*CONJG(cza))
      end if
      GOTO 9
  5   if (ls) j=j+1 ! continuity equation, no contribution from DBSTR
      if (lc) j=j+1 ! continuity equation, no contribution from DBQSTR
      if (lshift2.and.le) j=j+1 ! continuity equation, no contribution from DBSEN
      GOTO 9
   9  J=J0
      if (ls) then
       j=j+1
c       DRA(I,J)=DRA(I,J)*DRFAC(1,I) do no longer, other multiplication used in lkcombine equations
      endif
      if (lc) then
       j=j+1
c       DRA(I,J)=DRA(I,J)*DRFAC(1,I)
      endif
      if (lshift2.and.le) then
        j=j+1
c        DRA(I,J)=DRA(I,J)*DRFAC(1,I)
      endif
      END DO
      if (ls) j0=j0+1
      if (lc) j0=j0+1
      if (lshift2.and.le) j0=j0+1
      end if
      END DO
      ENDIF
      IF (IDOMAINTYPE(ISTR).EQ.2) THEN
C
C     Open slurry wall
C
      DO INOD=INOD1,INODL
      DO I=1,M
      J=J0
      CZ=CZI(I)
      CZA=CALPH(I)
      IEQS=ITYPE(I)
C
C ITYPE=+1 potential specified at CZ (Note: at doublet nodes ITYPE=+1,
C          but potential is expressed in terms of local strenght parameters.)
C       -1 difference in potential specified: PHI(CZ)-PHI(CZA)
C       +2 stream function specified at CZ
C       -2 flow across line between CZ and CZA, positive to the left when at CZ
C          note: because PSI is continuous for line doublets, this is the same as
C          the difference bwteen CZ and CZA: PSI(CZ)-PSI(CZA)
C       +3 discharge component normal to the unit vector CZA (rotated to the left)
C       +4 discharge component parallel to the unit vector CZA
C        5 continuity equation: provide total discharge
C        6 generate a zero matrix coefficient. Used to create an element specific equation.
C
      LNEG=IEQS.LT.0
      IEQ=IABS(IEQS)
      IF (INOD.GT.INOD1) THEN
        CALL DBFILL (INOD1,INOD-1,0,CZ)
        RCZTEST1=ABS(AIMAG(CDDBLN(INOD-1)))
        RXA1=0.5*(REAL(CDBZEL(INOD-1))+1.0)
      ENDIF
      CALL DBFILL (INOD1,INOD,0,CZ)
        RCZTEST2=ABS(AIMAG(CDDBLN(INOD)))
      GOTO (11,12,13,14,15,11),IEQ
C
   11 IF (INOD.GT.INOD1) THEN
        CDUM=CFDBG(INOD-1)+CFDBF(INOD) ! provide potential at CZ
        J=J+1
        DRA(I,J)=REAL(CDUM)
      ENDIF
      CDUM=CFDBS(INOD)
      J=J+1
      DRA(I,J)=REAL(CDUM)
      GOTO 19
  12  IF (INOD.GT.INOD1) THEN
        CDUM=CFDBG(INOD-1)+CFDBF(INOD)
        J=J+1
        DRA(I,J)=AIMAG(CDUM)    ! provide PSI at CZ due to DBSTR(INOD)
      ENDIF
      J=J+1
      CDUM=CFDBS(INOD)
      DRA(I,J)=AIMAG(CDUM)      ! provide PSI at CZ due to DBQSTR(INOD)
      IF (LNEG) THEN      ! flow across line through CZ and CZA condition
        CALL DBFILL (INOD1,INOD,0,CZA)
        IF (INOD.GT.INOD1) THEN    ! no element before node INOD1
          J=J-1
          CALL DBFILL (INOD1,INOD-1,0,CZA)
          CDUM=CFDBG(INOD-1)+CFDBF(INOD)
          DRA(I,J)=DRA(I,J)-AIMAG(CDUM)   ! subtract PSI at CZA due to DBSTR(INOD)
        J=J+1
        ENDIF
        CALL DBFILL (INOD1,INOD,0,CZA)
        CDUM=CFDBS(INOD)
        DRA(I,J)=DRA(I,J)-AIMAG(CDUM)     ! subtract PSI at CZA due to DBQSTR(INOD)
      ENDIF
      GOTO 19
  13  IF (INOD.GT.INOD1) THEN
        CDUM=CFDBGD(INOD-1)+CFDBFD(INOD)  ! provide Q normal to unit vector CZA
        J=J+1
        DRA(I,J)=+AIMAG(CDUM*CONJG(CZA))
      ENDIF
      CDUM=CFDBSD(INOD)
      J=J+1
      DRA(I,J)=AIMAG(CDUM*CONJG(CZA))
      GOTO 19
  14  IF (INOD.GT.INOD1) THEN
        CDUM=CFDBGD(INOD-1)+CFDBFD(INOD)  ! provide Q parallel to unit vector CZA
        J=J+1
        DRA(I,J)=REAL(CDUM*CONJG(CZA))
      ENDIF
      CDUM=CFDBSD(INOD)
      J=J+1
      DRA(I,J)=REAL(CDUM*CONJG(CZA))
      GOTO 19
  15  J=J+2                             ! continuity equation, no contribution from DBSTR & DBQSTR
      GOTO 19
  19  continue
      END DO
      IF (INOD.EQ.INOD1) THEN
      J0=J0+1
      ELSE
      J0=J0+2
      END IF
      END DO
      END IF
      IF (IDOMAINTYPE(ISTR).EQ.3) THEN
C
C     Closed slurry wall
C
      DO INOD=INOD1,INODL
      if (inod.eq.inod1) then
      inodm1=inodl
      else
      inodm1=inod-1
      end if
      if (inod.eq.inodl) then
      inodp1=inod1
      else
      inodp1=inod+1
      end if
      DO I=1,M
      J=J0
      CZ=CZI(I)
      CZA=CALPH(I)
      IEQS=ITYPE(I)
C
C ITYPE=+1 potential specified at CZ (Note: at doublet nodes ITYPE=+1,
C          but potential is expressed in terms of local strenght parameters.)
C       -1 difference in potential specified: PHI(CZ)-PHI(CZA)
C       +2 stream function specified at CZ
C       -2 flow across line between CZ and CZA, positive to the left when at CZ
C          note: because PSI is continuous for line doublets, this is the same as
C          the difference bwteen CZ and CZA: PSI(CZ)-PSI(CZA)
C       +3 discharge component normal to the unit vector CZA (rotated to the left)
C       +4 discharge component parallel to the unit vector CZA
C        5 continuity equation: provide total discharge
C        6 zero matrix element requested
C
      LNEG=IEQS.LT.0
      IEQ=IABS(IEQS)
      CALL DBFILL (INOD1,INODM1,INODL,CZ)
      CALL DBFILL (INOD1,INOD,INODL,CZ)
      GOTO (21,22,23,24,25,21),IEQ
C
   21 CDUM=CFDBG(INODM1)+CFDBF(INOD) ! provide potential at CZ
      J=J+1
      DRA(I,J)=REAL(CDUM)
      CDUM=CFDBS(INOD)
      J=J+1
      DRA(I,J)=REAL(CDUM)
      GOTO 29
  22  CDUM=CFDBG(INODM1)+CFDBF(INOD) ! provide PSI at CZ
      J=J+1
      DRA(I,J)=AIMAG(CDUM)
      J=J+1
      CDUM=CFDBS(INOD)
      DRA(I,J)=AIMAG(CDUM)
      IF (LNEG) THEN                 ! subtract PSI at CZA
        J=J-1
        CALL DBFILL (INOD1,INODM1,INODL,CZA)
        CALL DBFILL (INOD1,INOD,INODL,CZA)
        CDUM=CFDBG(INODM1)+CFDBF(INOD)
        DRA(I,J)=DRA(I,J)-AIMAG(CDUM)
        J=J+1
        CALL DBFILL (INOD1,INOD,INODL,CZA) ! necessary because of RFPERM(CZ) calls above  MAY NOT, TEST THIS OUT
        CDUM=CFDBS(INOD)
        DRA(I,J)=DRA(I,J)-AIMAG(CDUM)
      ENDIF
      GOTO 29
  23  CDUM=CFDBGD(INODM1)+CFDBFD(INOD)  ! provide Q normal to unit vector CZA    ! currently CZA is not a unit vector
      J=J+1
      DRA(I,J)=AIMAG(CDUM*CONJG(CZA))
      CDUM=CFDBSD(INOD)
      J=J+1
      DRA(I,J)=AIMAG(CDUM*CONJG(CZA))
      GOTO 29
  24  CDUM=CFDBGD(INODM1)+CFDBFD(INOD)  ! provide Q parallel to unit vector CZA
      J=J+1
      DRA(I,J)=REAL(CDUM*CONJG(CZA))
      CDUM=CFDBSD(INOD)
      J=J+1
      DRA(I,J)=REAL(CDUM*CONJG(CZA))
      GOTO 29
  25  J=J+2                             ! continuity equation, no contribution from DBSTR & DBQSTR
      GOTO 29
  29  continue
      END DO
      J0=J0+2
!      WRITE (ILUME,1000) N,J0
      END DO
      END IF
  30  CONTINUE
C ----------------reset CZ0 in DBPREP since DBFILL calls have changed
C                 common block contents  
      CZ=(1.0E21,1.0E21)
      CALL DBPREP (CZ)
      RETURN
 1000 FORMAT ('+Generating',I4,' equations, doing equation #: ',I4)
      END
c
C ---------------------------------------------------------------------------------------------------------------
c
      SUBROUTINE DBMATCORRECT (DRA,CZI,M,N,J0,DRFAC,CALPH,ITYPE)  !  ! 7/7/03 for T inhomogeneities only
c
c ---------------------------------------------------------------------------------------------------------------
C
C     Routine corrects the matrix coefficients for the line doublets
C     Two columns are generated at the same time for DBSTR and DBQSTR
C     of each line doublet.
c     Note: for inhomogeneities with common boundaries 3 columns may be generated per line-doublet,
c     for DBSTR, DBQSTR, and DBSTE
c
C
      IMPLICIT NONE
      INTEGER(4) M,N,J0,ITYPE,ISTR,INOD1,INODL,INODM1,INOD,I,J,
     &           IEQS,IEQ,INODP1
      LOGICAL LNEG,lsamedomain,ls,lc,le,lsm,lcm,lem,
     &        lshift0,lshift1,lshift2
      REAL(8) DRA,DRFAC,RCZTEST1,RXA1,RCZTEST2,RXA2,RCZATEST2,RCZATEST1,
     &        RXB2,RXB1,RCORRECTION,RFINTMU,RFINTLABLEFT,RFINTLABRIGHT,
     &        RFPERM
      COMPLEX(8) CZI,CALPH,CFDBF,CFDBG,CFDBS,CFDBFD,CFDBGD,CFDBSD,CZ,
     &           CZA,CDUM
      INCLUDE 'dbcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      DIMENSION DRA(M,*),CZI(*),DRFAC(4,*),CALPH(*),ITYPE(*)
      IF (NDB.EQ.0) RETURN
      DO 30 ISTR=1,NDBSTR
      INOD1=IDBSTA(ISTR)
      IF (ldbrechargeonly(istr)) GOTO 30 !  "recharge only", skip domain
      INODL=INOD1+NDBSTI(ISTR)-1
      INODM1=INODL
      IF (IDOMAINTYPE(ISTR).EQ.1) THEN
C
C     Hydraulic conductivity inhomogeneity
C
      DO INOD=INOD1,INODL            ! loop over line-doublets
c
c corrections apply to collocation points in current line-doublet and, for case 1
c (collocation point at first node) also apply to collocation points on previous line-doublet.
c
      if (idbcode(inod).gt.0) then   ! line-doublet is indeed in the matrix
      inodm1=inod-1
      if (inod.eq.inod1) inodm1=inodl
      inodp1=inod+1
      if (inod.eq.inodl) inodp1=inod1
      lshift0=idbcode(inodm1).eq.2
      lshift1=idbcode(inod).eq.2   ! shift first collocation point away from first node
      lshift2=idbcode(inodp1).eq.2.or.idbcode(inodp1).eq.0 ! add an extra collocation point shated away from the last node
      ls=rdbtfacnod(inod).ne.1.0d21 ! Ti<>To at this collocation point
      lc=rdbtfacctr(inod).ne.1.0d21 ! Ti<>To at this collocation point
      le=rdbtfacend(inod).ne.1.0d21 ! Ti<>To at this collocation point
      lsm=rdbtfacnod(inodm1).ne.1.0d21.and.idbcode(inodm1).gt.0  ! collocation point is indeed in matrix
      lcm=rdbtfacctr(inodm1).ne.1.0d21.and.idbcode(inodm1).gt.0  ! collocation point is indeed in matrix
      lem=rdbtfacend(inodm1).ne.1.0d21.and.idbcode(inodm1).gt.0  ! collocation point is indeed in matrix
      DO I=1,M
      J=J0
      CZ=CZI(I)
      CZA=CALPH(I)
      IEQS=ITYPE(I)
C
C ITYPE=+1 potential specified at CZ (Note: at doublet nodes ITYPE=+1,
C          but potential is expressed in terms of local strenght parameters.)
C       -1 difference in potential specified: PHI(CZ)-PHI(CZA)
C       +2 stream function specified at CZ
C       -2 flow across line between CZ and CZA, positive to the left when at CZ
C          note: because PSI is continuous for line doublets, this is the same as
C          the difference bwteen CZ and CZA: PSI(CZ)-PSI(CZA)
C       +3 discharge component normal to the unit vector CZA (rotated to the left)
C       +4 discharge component parallel to the unit vector CZA
C        5 continuity equation: provide total discharge
C        6 generate a zero matrix coefficient. Used to create an element specific equation.
C
      LNEG=IEQS.LT.0
      IEQ=IABS(IEQS)
      CALL DBFILL (INOD1,INODM1,INODL,CZ)
      CALL DBFILL (INOD1,INOD,INODL,CZ)
      call dbfill (inod1,inodp1,inodl,cz)
      GOTO (1,9,9,9,9,6),IEQ              ! 2,3,4,5 replaced by 9, since no corrections neeeded
C
   1  if (ls) THEN ! there is indeed a DBSTR(inod) in the matrix
      select case (idbcode(inod))
                case (1)   ! DBSTR for line-doublet INOD is also DBSEN for line-doublet INODM1
       J=J+1                ! corrections for DBSTR
       IF (LDBNOD(INOD).and.ls) DRA(I,J)=DRA(I,J)-RDBTFACNOD(INOD) ! if at first node of line-doublet INOD
       IF (ABS(CDBZEL(INODM1)).LT.1.0E-4.and.lcm)                   ! if at center of line-doublet INODM1
     & DRA(I,J)=DRA(I,J)-0.5*RDBTFACCTR(INODM1)
       IF (ABS(CDBZEL(INOD)).LT.1.0E-4.and.lc)                     ! if at center of line-doublet INOD
     & DRA(I,J)=DRA(I,J)-0.5*RDBTFACCTR(INOD)
       if (lshift0) then ! if at (shifted) collocation point near first node of line-doublet INODM1
        if (ABS(cdbzel(inodm1)+1.0d0-rdboffsetinhom).lt.1.0d-4.and.lsm)
     &  dra(i,j)=dra(i,j)-0.5*rdboffsetinhom*rdbtfacnod(inodm1)
       end if
       if (lshift2) then              ! if at (shifted) collocation point near last node of line-doublet INOD
        if (ABS(cdbzel(inod)-1.0d0+rdboffsetinhom).lt.1.0d-4.and.le)
     &  dra(i,j)=dra(i,j)-0.5*rdboffsetinhom*rdbtfacend(inod)
       end if
               case (2)   ! DBSTR(INOD) applies to line-doublet INOD only
       j=j+1                      ! corrections for DBSTR
       if (ABS(cdbzel(inod)+1.0d0-rdboffsetinhom).lt.1.0d-4.and.ls) ! if at (shifted) collocation point near first node of line-doublet INOD
     & dra(i,j)=dra(i,j)-0.5*(2.0d0-rdboffsetinhom)*rdbtfacnod(inod)
       IF (ABS(CDBZEL(INOD)).LT.1.0E-4.and.lc)                     ! if at center of line-doublet INOD
     & DRA(I,J)=DRA(I,J)-0.5*RDBTFACCTR(INOD)
       if (lshift2) then
        if (ABS(cdbzel(inod)-1.0d0+rdboffsetinhom).lt.1.0d-4.and.le) ! if at (shifted) collocation point near last node of line-doublet INOD
     &  dra(i,j)=dra(i,j)-0.5*rdboffsetinhom*rdbtfacend(inod)
       endif
      end select
      endif
c
      if (lc) then ! there is indeed a DBQSTR(inod) in the matrix
       J=J+1                 ! corrections for DBQSTR
       IF (ABS(CDBZEL(INOD)).LT.1.0E-4.and.lc) ! if at center of line-doublet
     & DRA(I,J)=DRA(I,J)-RDBTFACCTR(INOD)
       if (lshift1) then              ! if at (shifted) collocation point near first node of line-doublet INOD
        if (ABS(cdbzel(inod)+1.0d0-rdboffsetinhom).lt.1.0d-4.and.ls)
     &  dra(i,j)=dra(i,j)-(2.0d0*rdboffsetinhom-
     &           rdboffsetinhom*rdboffsetinhom)*rdbtfacnod(inod)
       end if
       if (lshift2) then              ! if at (shifted) collocation point near last node of line-doublet INOD
        if (ABS(cdbzel(inod)-1.0d0+rdboffsetinhom).lt.1.0d-4.and.le)
     &  dra(i,j)=dra(i,j)-(2.0d0*rdboffsetinhom-
     &           rdboffsetinhom*rdboffsetinhom)*rdbtfacend(inod)
       end if
      endif
c
      if (lshift2.and.le) then   ! there is indeed a DBSEN(inod) in the matrix
       j=j+1                      ! corrections for DBSEN
       if (lshift1) then
        if (ABS(cdbzel(inod)+1.0d0-rdboffsetinhom).lt.1.0d-4.and.ls) ! if at (shifted) collocation point near first node of line-doublet INOD
     &  dra(i,j)=dra(i,j)-0.5*rdboffsetinhom*rdbtfacnod(inod)
       end if
       IF (ABS(CDBZEL(INOD)).LT.1.0E-4.and.lc)                     ! if at center of line-doublet INOD
     & DRA(I,J)=DRA(I,J)-0.5*RDBTFACCTR(INOD)
       if (ABS(cdbzel(inod)-1.0d0+rdboffsetinhom).lt.1.0d-4.and.le) ! if at (shifted) collocation point near last node of line-doublet INOD
     & dra(i,j)=dra(i,j)-0.5*(2.0d0-rdboffsetinhom)*rdbtfacend(inod)
      end if
c
      IF (LNEG) THEN                 ! subtract potential at CZA, not yet done in DBGENMAT
       j=j0
       CALL DBFILL (INOD1,INODM1,INODL,CZA)
       CALL DBFILL (INOD1,INOD,INODL,CZA)
       call dbfill (inod1,inodp1,inodl,cza)
       if (ls) then
        j=j+1
        CDUM=CFDBG(INODM1)+CFDBF(INOD)
        if (lshift1) cdum=cfdbf(inod)
        DRA(I,J)=DRA(I,J)-REAL(CDUM)
       end if
       if (lc) then
        J=J+1
        CDUM=CFDBS(INOD)
        DRA(I,J)=DRA(I,J)-REAL(CDUM)
       end if
       if (lshift2.and.le) then
        cdum=cfdbg(inod)
        j=j+1
        dra(i,j)=dra(i,j)-REAL(cdum)
       end if
      ENDIF
      GOTO 9
c
  6   if (ls) THEN ! there is indeed a DBSTR(inod) in the matric
      select case (idbcode(inod))   ! equation for the case that Ti=To, but Bi<>Bo
             case (1) ! DBSTR applies to both line-doublet INOD and line-doublet INODM1
       j=j+1          ! corrections for DBSTR
       DRA(I,J)=0.0
       IF (LDBNOD(INOD).and.ls) DRA(I,J)=1.0D0  ! if at first node of line-doublet INOD
       IF (ABS(CDBZEL(INODM1)).LT.1.0E-4.and.lcm) ! if at center of previous line-doublet (INODM1)
     & DRA(I,J)=0.5D0
       IF (ABS(CDBZEL(INOD)).LT.1.0E-4.and.lc) ! if at center of current line-doublet (INOD)
     & DRA(I,J)=0.5D0
       if (lshift0) then ! if at (shifted) collocation point near first node of line-doublet INODM1
       if (ABS(cdbzel(inodm1)+1.0d0-rdboffsetinhom).lt.1.0d-4.and.lsm) ! if at (shifted) collocation point near first node of previous line-doublet (INODM1)
     & dra(i,j)=0.5d0*rdboffsetinhom
       end if
       if (lshift2) then              ! if at (shifted) collocation point near last node of line-doublet INOD
       if (ABS(cdbzel(inod)-1.0d0+rdboffsetinhom).lt.1.0d-4.and.le)
     & dra(i,j)=0.5d0*rdboffsetinhom
       end if
            case (2)   ! DBSTR(INOD) applies to line-doublet INOD only
       j=j+1                      ! corrections for DBSTR
       if (ABS(cdbzel(inod)+1.0d0-rdboffsetinhom).lt.1.0d-4.and.ls) ! if at (shifted) collocation point near first node of line-doublet INOD
     & dra(i,j)=0.5d0*(2.0d0-rdboffsetinhom)
       IF (ABS(CDBZEL(INOD)).LT.1.0E-4.and.lc)                     ! if at center of line-doublet INOD
     & DRA(I,J)=0.5D0
       if (ABS(cdbzel(inod)-1.0d0+rdboffsetinhom).lt.1.0d-4.and.le) ! if at (shifted) collocation point near last node of line-doublet INOD
     & dra(i,j)=0.5d0*rdboffsetinhom
      end select
      endif
c
      if (lc) THEN ! there is indeed a DBQSTR(inod) in the matrix
       J=J+1                 ! corrections for DBQSTR
       IF (ABS(CDBZEL(INOD)).LT.1.0E-4.and.le)  ! if at center of line-doublet INOD
     & DRA(I,J)=1.0d0
       if (lshift1) then              ! if at (shifted) collocation point near first node of line-doublet INOD
        if (ABS(cdbzel(inod)+1.0d0-rdboffsetinhom).lt.1.0d-4.and.ls)
     &  dra(i,j)=(2.0d0*rdboffsetinhom-
     &           rdboffsetinhom*rdboffsetinhom)
       end if
       if (lshift2) then              ! if at (shifted) collocation point near last node of line-doublet INOD
        if (ABS(cdbzel(inod)-1.0d0+rdboffsetinhom).lt.1.0d-4.and.le)
     &  dra(i,j)=(2.0d0*rdboffsetinhom-
     &           rdboffsetinhom*rdboffsetinhom)
       end if
      endif
c
      if (lshift2.and.le) then   ! there is indeed a DBSEN(inod) in the matrix
       j=j+1                      ! corrections for DBSEN
       if (lshift1) then
        if (ABS(cdbzel(inod)+1.0d0-rdboffsetinhom).lt.1.0d-4.and.ls) ! if at (shifted) collocation point near first node of line-doublet INOD
     &  dra(i,j)=0.5*rdboffsetinhom
       end if
       IF (ABS(CDBZEL(INOD)).LT.1.0E-4.and.lc)                     ! if at center of line-doublet INOD
     & DRA(I,J)=0.5d0
       if (lshift2) then
        if (ABS(cdbzel(inod)-1.0d0+rdboffsetinhom).lt.1.0d-4.and.le) ! if at (shifted) collocation point near last node of line-doublet INOD
     &  dra(i,j)=0.5d0*(2.0d0-rdboffsetinhom)
       endif
      end if
   9  continue
      END DO
      if (ls) j0=j0+1
      if (lc) j0=j0+1
      if (lshift2.and.le) j0=j0+1
      endif
      END DO
      ENDIF
      IF (IDOMAINTYPE(ISTR).EQ.2) THEN
C
C     Open slurry wall
C
      DO INOD=INOD1,INODL
      DO I=1,M
      J=J0
      CZ=CZI(I)
      CZA=CALPH(I)
      IEQS=ITYPE(I)
C
C ITYPE=+1 potential specified at CZ (Note: at doublet nodes ITYPE=+1,
C          but potential is expressed in terms of local strenght parameters.)
C       -1 difference in potential specified: PHI(CZ)-PHI(CZA)
C       +2 stream function specified at CZ
C       -2 flow across line between CZ and CZA, positive to the left when at CZ
C          note: because PSI is continuous for line doublets, this is the same as
C          the difference bwteen CZ and CZA: PSI(CZ)-PSI(CZA)
C       +3 discharge component normal to the unit vector CZA (rotated to the left)
C       +4 discharge component parallel to the unit vector CZA
C        5 continuity equation: provide total discharge
C        6 generate a zero matrix coefficient. Used to create an element specific equation.
C
      LNEG=IEQS.LT.0
      IEQ=IABS(IEQS)
      IF (INOD.GT.INOD1) THEN
        CALL DBFILL (INOD1,INOD-1,0,CZ)
        RCZTEST1=ABS(AIMAG(CDDBLN(INOD-1)))
        RXA1=0.5*(REAL(CDBZEL(INOD-1))+1.0)
      ENDIF
      CALL DBFILL (INOD1,INOD,0,CZ)
        RCZTEST2=ABS(AIMAG(CDDBLN(INOD)))
      GOTO (11,12,13,14,15,16),IEQ
C
   11 IF (INOD.GT.INOD1) THEN
        J=J+1
      ENDIF
      J=J+1
      IF (LNEG) THEN                 ! subtract potential at CZA
       J=J-1
       IF (INOD.GT.INOD1) THEN
         CALL DBFILL (INOD1,INOD-1,0,CZA)
       ENDIF
       CALL DBFILL (INOD1,INOD,0,CZA)
       IF (INOD.GT.INOD1) THEN
         CDUM=CFDBG(INOD-1)+CFDBF(INOD)
         DRA(I,J)=DRA(I,J)-REAL(CDUM)
       ENDIF
       J=J+1
       CDUM=CFDBS(INOD)
       DRA(I,J)=DRA(I,J)-REAL(CDUM)
      ENDIF
      GOTO 19
  12  IF (INOD.GT.INOD1) THEN
        J=J+1
      ENDIF
      J=J+1
      IF (LNEG) THEN      ! flow across line through CZ and CZA condition
        RXA2=0.5*(REAL(CDBZEL(INOD))+1.0)
        CALL DBFILL (INOD1,INOD,0,CZA)
        RCZATEST2=ABS(AIMAG(CDDBLN(INOD)))
        RCZATEST1=0.0
        RXB2=0.5*(REAL(CDBZEL(INOD))+1.0)
        IF (INOD.GT.INOD1) THEN    ! no element before node INOD1
          J=J-1
          CALL DBFILL (INOD1,INOD-1,0,CZA)
          RCZATEST1=ABS(AIMAG(CDDBLN(INOD-1)))
          RXB1=0.5*(REAL(CDBZEL(INOD-1))+1.0)
c          CDUM=CFDBG(INOD-1)+CFDBF(INOD)
c          DRA(I,J)=DRA(I,J)-AIMAG(CDUM)   ! subtract PSI at CZA due to DBSTR(INOD)
          IF (RCZTEST1.LT.3.0.AND.RCZATEST1.GT.3.0) THEN ! only CZA on line doublet between node INOD-1 & INOD
            RCORRECTION=RDBTFACNOD(INOD-1)/RFPERM(CZA)*
     &                  ABS(CDBZ21(INOD-1))*RFINTLABLEFT(0.0,RXB1)
            DRA(I,J)=DRA(I,J)-RCORRECTION
           END IF
          IF (RCZTEST1.GT.3.0.AND.RCZATEST1.GT.3.0) THEN ! CZ & CZA on line doublet between node INOD-1 & INOD
            RCORRECTION=RDBTFACCTR(INOD-1)/RFPERM(CZ)*
     &                  ABS(CDBZ21(INOD-1))*RFINTLABLEFT(RXA1,RXB1)
            DRA(I,J)=DRA(I,J)-RCORRECTION
          END IF
          IF (RCZTEST1.GT.3.0.AND.RCZATEST2.GT.3.0) THEN ! CZ between INOD-1 & INOD, CZA between INOD & INOD+1
            RCORRECTION=RDBTFACNOD(INOD)/RFPERM(CZ)*
     &                  (ABS(CDBZ21(INOD-1))*RFINTLABLEFT(RXA1,1.0)+
     &                   ABS(CDBZ21(INOD))*RFINTLABRIGHT(0.0,RXB2))
            DRA(I,J)=DRA(I,J)-RCORRECTION
          END IF
          IF (RCZTEST2.GT.3.0.AND.RCZATEST2.GT.3.0) THEN ! CZ & CZA on line doublet between node INOD & INOD+1
            RCORRECTION=RDBTFACCTR(INOD)/RFPERM(CZ)*
     &                  ABS(CDBZ21(INOD))*RFINTLABRIGHT(RXA2,RXB2)
            DRA(I,J)=DRA(I,J)-RCORRECTION
          END IF
          IF (RCZTEST2.GT.3.0.AND.RCZATEST2.LT.3.0) THEN ! only CZ on line doublet between node INOD & INOD+1
            RCORRECTION=RDBTFACNOD(INOD+1)/RFPERM(CZ)*
     &                  ABS(CDBZ21(INOD))*RFINTLABRIGHT(RXA2,1.0)
            DRA(I,J)=DRA(I,J)-RCORRECTION
          END IF
        J=J+1
        ENDIF
        CALL DBFILL (INOD1,INOD,0,CZA) ! necessary because of RFPERM(CZ) call above
c        CDUM=CFDBS(INOD)
c        DRA(I,J)=DRA(I,J)-AIMAG(CDUM)     ! subtract PSI at CZA due to DBQSTR(INOD)
        IF (RCZTEST2.LT.3.0.AND.RCZATEST2.GT.3.0) THEN ! only CZA on line doublet between INOD & INOD+1
          RCORRECTION=RDBTFACNOD(INOD)/RFPERM(CZA)*
     &               ABS(CDBZ21(INOD))*RFINTMU(0.0,RXB2)
          DRA(I,J)=DRA(I,J)-RCORRECTION
        END IF
        IF (RCZTEST2.GT.3.0.AND.RCZATEST2.GT.3.0) THEN ! CZ & CZA on line doublet between node INOD & INOD+1
          RCORRECTION=RDBTFACCTR(INOD)/RFPERM(CZ)*
     &                ABS(CDBZ21(INOD))*RFINTMU(RXA2,RXB2)
          DRA(I,J)=DRA(I,J)-RCORRECTION
        END IF
        IF (RCZTEST2.GT.3.0.AND.RCZATEST2.LT.3.0) THEN ! only CZ on line doublet between node INOD & INOD+1
          RCORRECTION=RDBTFACNOD(INOD+1)/RFPERM(CZ)*
     &               ABS(CDBZ21(INOD))*RFINTMU(RXA2,1.0)
          DRA(I,J)=DRA(I,J)-RCORRECTION
        END IF
      ENDIF
      GOTO 19
  13  IF (INOD.GT.INOD1) THEN
        J=J+1
      ENDIF
      J=J+1
      GOTO 19
  14  IF (INOD.GT.INOD1) THEN
        J=J+1
      ENDIF
      J=J+1
      GOTO 19
  15  J=J+2
      GOTO 19
  16  j=j+1          ! zero matrix elements requested
      DRA(I,J)=0.0
      J=J+1
      DRA(I,J)=0.0
      GOTO 19
  19  continue
      END DO
      IF (INOD.EQ.INOD1) THEN
      J0=J0+1
      ELSE
      J0=J0+2
      END IF
      END DO
      END IF
      IF (IDOMAINTYPE(ISTR).EQ.3) THEN
C
C     Closed slurry wall
C
      DO INOD=INOD1,INODL
      if (inod.eq.inod1) then
      inodm1=inodl
      else
      inodm1=inod-1
      end if
      if (inod.eq.inodl) then
      inodp1=inod1
      else
      inodp1=inod+1
      end if
      DO I=1,M
      J=J0
      CZ=CZI(I)
      CZA=CALPH(I)
      IEQS=ITYPE(I)
C
C ITYPE=+1 potential specified at CZ (Note: at doublet nodes ITYPE=+1,
C          but potential is expressed in terms of local strenght parameters.)
C       -1 difference in potential specified: PHI(CZ)-PHI(CZA)
C       +2 stream function specified at CZ
C       -2 flow across line between CZ and CZA, positive to the left when at CZ
C          note: because PSI is continuous for line doublets, this is the same as
C          the difference bwteen CZ and CZA: PSI(CZ)-PSI(CZA)
C       +3 discharge component normal to the unit vector CZA (rotated to the left)
C       +4 discharge component parallel to the unit vector CZA
C        5 continuity equation: provide total discharge
C        6 zero matrix element requested
C
      LNEG=IEQS.LT.0
      IEQ=IABS(IEQS)
      CALL DBFILL (INOD1,INODM1,INODL,CZ)
      CALL DBFILL (INOD1,INOD,INODL,CZ)
      RCZTEST1=ABS(AIMAG(CDDBLN(INODM1)))
      RXA1=0.5*(REAL(CDBZEL(INODM1))+1.0)
      RCZTEST2=ABS(AIMAG(CDDBLN(INOD)))
      RXA2=0.5*(REAL(CDBZEL(INOD))+1.0)
      GOTO (21,22,23,24,25,26),IEQ
C
   21 continue
      J=J+2
      IF (LNEG) THEN                 ! subtract potential at CZA
      J=J-1
      CALL DBFILL (INOD1,INODM1,INODL,CZA)
      CALL DBFILL (INOD1,INOD,INODL,CZA)
      CDUM=CFDBG(INODM1)+CFDBF(INOD)
      DRA(I,J)=DRA(I,J)-REAL(CDUM)
      J=J+1
      CDUM=CFDBS(INOD)
      DRA(I,J)=DRA(I,J)-REAL(CDUM)
      ENDIF
      GOTO 29
  22  j=j+2
      IF (LNEG) THEN                 ! subtract PSI at CZA
        J=J-1
        CALL DBFILL (INOD1,INODM1,INODL,CZA)
        CALL DBFILL (INOD1,INOD,INODL,CZA)
        RCZATEST1=ABS(AIMAG(CDDBLN(INODM1)))
        RXB1=0.5*(REAL(CDBZEL(INODM1))+1.0)
        RCZATEST2=ABS(AIMAG(CDDBLN(INOD)))
        RXB2=0.5*(REAL(CDBZEL(INOD))+1.0)
c        CDUM=CFDBG(INODM1)+CFDBF(INOD)
c        DRA(I,J)=DRA(I,J)-AIMAG(CDUM)
        ! CORRECTION?
          IF (RCZTEST1.LT.3.0.AND.RCZATEST1.GT.3.0) THEN ! only CZA on line doublet between node INOD-1 & INOD
            RCORRECTION=RDBTFACNOD(INODM1)/RFPERM(CZA)*
     &                  ABS(CDBZ21(INODM1))*RFINTLABLEFT(0.0,RXB1)
            DRA(I,J)=DRA(I,J)-RCORRECTION
           END IF
          IF (RCZTEST1.GT.3.0.AND.RCZATEST1.GT.3.0) THEN ! CZ & CZA on line doublet between node INOD-1 & INOD
            RCORRECTION=RDBTFACCTR(INODM1)/RFPERM(CZ)*
     &                  ABS(CDBZ21(INODM1))*RFINTLABLEFT(RXA1,RXB1)
            DRA(I,J)=DRA(I,J)-RCORRECTION
          END IF
          IF (RCZTEST1.GT.3.0.AND.RCZATEST2.GT.3.0) THEN ! CZ between INOD-1 & INOD, CZA between INOD & INOD+1
            RCORRECTION=RDBTFACNOD(INOD)/RFPERM(CZ)*
     &                  (ABS(CDBZ21(INODM1))*RFINTLABLEFT(RXA1,1.0)+
     &                   ABS(CDBZ21(INOD))*RFINTLABRIGHT(0.0,RXB2))
            DRA(I,J)=DRA(I,J)-RCORRECTION
          END IF
          IF (RCZTEST2.GT.3.0.AND.RCZATEST2.GT.3.0) THEN ! CZ & CZA on line doublet between node INOD & INOD+1
            RCORRECTION=RDBTFACCTR(INOD)/RFPERM(CZ)*
     &                  ABS(CDBZ21(INOD))*RFINTLABRIGHT(RXA2,RXB2)
            DRA(I,J)=DRA(I,J)-RCORRECTION
          END IF
          IF (RCZTEST2.GT.3.0.AND.RCZATEST2.LT.3.0) THEN ! only CZ on line doublet between node INOD & INOD+1
            RCORRECTION=RDBTFACNOD(INODP1)/RFPERM(CZ)*
     &                  ABS(CDBZ21(INOD))*RFINTLABRIGHT(RXA2,1.0)
            DRA(I,J)=DRA(I,J)-RCORRECTION
          END IF
        J=J+1
        CALL DBFILL (INOD1,INOD,INODL,CZA) ! necessary because of RFPERM(CZ) calls above
c        CDUM=CFDBS(INOD)
c        DRA(I,J)=DRA(I,J)-AIMAG(CDUM)
        ! CORRECTION?
        IF (RCZTEST2.LT.3.0.AND.RCZATEST2.GT.3.0) THEN ! only CZA on line doublet between INOD & INOD+1
          RCORRECTION=RDBTFACNOD(INOD)/RFPERM(CZA)*
     &               ABS(CDBZ21(INOD))*RFINTMU(0.0,RXB2)
          DRA(I,J)=DRA(I,J)-RCORRECTION
        END IF
        IF (RCZTEST2.GT.3.0.AND.RCZATEST2.GT.3.0) THEN ! CZ & CZA on line doublet between node INOD & INOD+1
          RCORRECTION=RDBTFACCTR(INOD)/RFPERM(CZ)*
     &                ABS(CDBZ21(INOD))*RFINTMU(RXA2,RXB2)
          DRA(I,J)=DRA(I,J)-RCORRECTION
        END IF
        IF (RCZTEST2.GT.3.0.AND.RCZATEST2.LT.3.0) THEN ! only CZ on line doublet between node INOD & INOD+1
          RCORRECTION=RDBTFACNOD(INODP1)/RFPERM(CZ)*
     &               ABS(CDBZ21(INOD))*RFINTMU(RXA2,1.0)
          DRA(I,J)=DRA(I,J)-RCORRECTION
        END IF
      ENDIF
      GOTO 29
  23  j=j+2
      GOTO 29
  24  j=j+2
      GOTO 29
  25  J=J+2                             ! continuity equation, no contribution from DBSTR & DBQSTR
      GOTO 29
  26  j=j+1          ! zero matrix elements requested
      DRA(I,J)=0.0
      J=J+1
      DRA(I,J)=0.0
      GOTO 29
  29  continue
      END DO
      J0=J0+2
!      WRITE (ILUME,1000) N,J0
      END DO
      END IF
  30  CONTINUE
C ----------------reset CZ0 in DBPREP since DBFILL calls have changed
C                 common block contents  
      CZ=(1.0E21,1.0E21)
      CALL DBPREP (CZ)
      RETURN
 1000 FORMAT ('+Generating',I4,' equations, doing equation #: ',I4)
      END
c
c
C ----------------------------------------------------------------------------------------------------------
c
      SUBROUTINE DBKNO (DRB,N,J,CZI,CALPH,ITYPE)
c
c ----------------------------------------------------------------------------------------------------------
C
C     Routine generates the known vector for the doublets
C      
      IMPLICIT NONE
      INTEGER(4) N,J,ITYPE,ISTR,INOD1,INODL,INOD,INODP1,INODM1,I,
     &           NHEADCONDITIONS
      LOGICAL LFDBIN,lshift1,lshift2,ls,lc,le
      REAL(8) DRB,RFPOT,RXA,RXB,RFINTLABLEFT,RFINTLABRIGHT,RFINTMU,
     &        RFPERM,RFNORMALFLOW,rlocs,rfdbslocal,rx,
     &        rdum1,rdum2
      COMPLEX(8) CZI,CALPH,CZRESTART,CZ,CZA
      INCLUDE 'dbcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      DATA CZRESTART /(1.0E21,1.0E21)/
      DIMENSION DRB(*),CZI(*),CALPH(*),ITYPE(*)
      IF (NDB.EQ.0) RETURN
      DO 20 ISTR=1,NDBSTR
      INOD1=IDBSTA(ISTR)
      IF (ldbrechargeonly(istr)) GOTO 20 ! "recharge only", skip domain
      INODL=INOD1+NDBSTI(ISTR)-1
      IF (IDOMAINTYPE(ISTR).EQ.1) THEN
c
c     Hydraulic conductivity inhomogeneity
c
      DO INOD=INOD1,INODL
      if (idbcode(inod).gt.0) then ! line-doublet is indeed included in the matrix
      inodp1=inod+1
      if (inod.eq.inodl) inodp1=inod1
      lshift1=idbcode(inod).eq.2
      lshift2=idbcode(inodp1).eq.2.or.idbcode(inodp1).eq.0
      ls=rdbtfacnod(inod).ne.1.0d21
      lc=rdbtfacctr(inod).ne.1.0d21
      le=rdbtfacend(inod).ne.1.0d21
      if (ls) then ! first collocation point does indeed exist
       J=J+1
       CZ=CZI(J)
       call dbprep(cz) ! the next two tests are for debugging purposes only
       if (lshift1) then
         if (ABS(cdbzel(inod)+1.0d0-rdboffsetinhom).gt.1.0d-4) then ! error, not at correct collocation point
           write (iluer,1100) cz,inod,istr
         end if
         rx=-1.0d0+rdboffsetinhom
         rlocs=rfdbslocal(rx,dbstr(inod),dbqstr(inod),dbsen(inod))
       else
         if (.not.ldbnod(inod)) THEN ! error, not at correct collocation point
           write (iluer,1200) cz,inod,istr
         end if
         rlocs=dbstr(inod)
       end if
       IF (ITYPE(J).EQ.6) THEN ! Ti=To, while Bi=/Bo
        DRB(J)=RDBBFACNOD(INOD)-rlocs
       ELSE                    ! Ti=/To
c        rdum1=0.0   !    ******************************************************* DEBUGGING
c        rdum1=rfpot(cz)
        rdum2=dbknon(inod)
c        write (ilume,1001) inod,rdum1,rdum2
c 1001   format ('dbkno1: inod,rfpot,dbknon ',i3,1x,2(d14.7))
        DRB(J)=DRB(J)-rdum2
        DRB(J)=DRB(J)+rlocs*RDBTFACNOD(INOD) ! correct for previous strength times transmissivity factor
        DRB(J)=DRB(J)-RDBTFACNOD(INOD)*RDBBFACNOD(INOD) ! subtract base jump factor, see RDBKFAC
       END IF
      endif
      if (lc) then
       J=J+1
       CZ=CZI(J)
       call dbprep(cz) ! the next test is for debugging purposes only
       if (ABS(cdbzel(inod)).gt.1.0d-4) THEN !error, not at correct collocation point
         write (iluer,1300) cz,inod,istr
       end if
       rx=0.0d0
       rlocs=rfdbslocal(rx,dbstr(inod),dbqstr(inod),dbsen(inod))
       IF (ITYPE(J).EQ.6) THEN ! Ti=To, while Bi=/Bo
        DRB(J)=RDBBFACCTR(INOD)-rlocs
       ELSE                     ! Ti=/To
c        rdum1=rfpot(cz) !    ********************************************************* DEBUGGING
        rdum2=dbknoc(inod)
c        write (ilume,1002) inod,rdum1,rdum2
c 1002   format ('dbkno2: inod,rfpot,dbknoc ',i3,1x,2(d14.7))
        DRB(J)=DRB(J)-rdum2
        DRB(J)=DRB(J)+rlocs*RDBTFACCTR(INOD) ! correct for previous strength times transmissivity factor
        DRB(J)=DRB(J)-RDBTFACCTR(INOD)*RDBBFACCTR(INOD) ! subtract base jump factor, see RDBKFAC
       END IF
      endif
      if (lshift2.and.le) then  ! there is one more collocation point near the end of the line-doublet
       J=J+1
       CZ=CZI(J)
       call dbprep(cz) ! the next test is for debugging purposes only
       if (ABS(cdbzel(inod)-1.0d0+rdboffsetinhom).gt.1.0d-4) THEN ! error, not at correct collocation point
        write (iluer,1400) cz,inod,istr
       end if
       rx=1.0d0-rdboffsetinhom
       rlocs=rfdbslocal(rx,dbstr(inod),dbqstr(inod),dbsen(inod))
       IF (ITYPE(J).EQ.6) THEN ! Ti=To, while Bi=/Bo
        DRB(J)=RDBBFACEND(INOD)-rlocs
       ELSE                    ! Ti=/To
c        rdum1=0.0   !    ******************************************************* DEBUGGING
c        rdum1=rfpot(cz)
        rdum2=dbknoe(inod)
c        write (ilume,1001) inod,rdum1,rdum2
c 1001   format ('dbkno1: inod,rfpot,dbknon ',i3,1x,2(d14.7))
        DRB(J)=DRB(J)-rdum2
        DRB(J)=DRB(J)+rlocs*RDBTFACEND(INOD) ! correct for previous strength times transmissivity factor
        DRB(J)=DRB(J)-RDBTFACEND(INOD)*RDBBFACEND(INOD) ! subtract base jump factor, see RDBKFAC
       END IF
      end if
      endif
      END DO
      ENDIF
      IF (IDOMAINTYPE(ISTR).EQ.2) THEN
C
C     Open slurry wall
C
      DO INOD=INOD1,INODL
      IF (INOD.EQ.INOD1) THEN
      INODM1=0
      ELSE
      INODM1=INOD-1
      END IF
      J=J+1
      CZ=CZI(J)
      CZA=CALPH(J)                 ! integration interval on line doublet between nodes INOD and INODP1
      CALL DBFILL (INOD1,INOD,0,CZ)
      RXA=0.5*(REAL(CDBZEL(INOD))+1.0)
      CALL DBFILL (INOD1,INOD,0,CZA)   ! note: call DBFILL with a zero end node to signal an open string
      RXB=0.5*(REAL(CDBZEL(INOD))+1.0)
      CALL DBPREP(CZRESTART)  ! reset DBPREP because arrays DBCOM are inconsistent after a single DBFILL call
      IF (INOD.EQ.INOD1) THEN
      DRB(J)=DRB(J)+RDBTFACCTR(INOD)/RFPERM(CZ)*
     &ABS(CDBZ21(INOD1))*(DBQSTR(INOD1)*RFINTMU(RXA,RXB)+
     &DBSTR(INOD+1)*RFINTLABLEFT(RXA,RXB))    ! correct for old strength terms
      ELSE
      DRB(J)=DRB(J)+RDBTFACCTR(INOD)/RFPERM(CZ)*
     &ABS(CDBZ21(INOD))*(DBQSTR(INOD)*RFINTMU(RXA,RXB)+
     &DBSTR(INOD)*RFINTLABRIGHT(RXA,RXB)+
     &DBSTR(INOD+1)*RFINTLABLEFT(RXA,RXB))      ! correct for old strength terms
      END IF
c      rdum1=rfnormalflow(cz,cza) ! ********************************************** DEBUGGING
      rdum2=dbknoc(inod)
c      write (ilume,1003) inod,rdum1,rdum2
c 1003 format ('dbkno3: inod,rfnormalflow,dbknoc ',i3,1x,2(d14.7))
      DRB(J)=DRB(J)-rdum2         ! subtract normal flow across interval
      IF (INOD.NE.INODL) THEN
      J=J+1
      CZ=CZI(J)
      CZA=CALPH(J)                 ! integration interval across node INOD+1
      CALL DBFILL (INOD1,INOD,0,CZ)
      RXA=0.5*(REAL(CDBZEL(INOD))+1.0)
      CALL DBFILL (INOD1,INOD+1,0,CZA)
      RXB=0.5*(REAL(CDBZEL(INOD+1))+1.0)
      CALL DBPREP(CZRESTART)  ! reset CZ0 in DBPREP because DBFILL changed data in common
      DRB(J)=DRB(J)+RDBTFACNOD(INOD+1)/RFPERM(CZ)*(
     &ABS(CDBZ21(INOD))*(DBQSTR(INOD)*RFINTMU(RXA,1.0)+
     &DBSTR(INOD)*RFINTLABRIGHT(RXA,1.0)+
     &DBSTR(INOD+1)*RFINTLABLEFT(RXA,1.0))+
     &ABS(CDBZ21(INOD+1))*(DBQSTR(INOD+1)*RFINTMU(0.0,RXB)+
     &DBSTR(INOD+1)*RFINTLABRIGHT(0.0,RXB)+
     &DBSTR(INOD+2)*RFINTLABLEFT(0.0,RXB)))    ! correct for old strength terms
c      rdum1=RFNORMALFLOW(CZ,CDBZ(INOD+1))+RFNORMALFLOW(CDBZ(INOD+1),CZA)
      rdum2=dbknon(inod)
c      write (ilume,1004) inod,rdum1,rdum2
c 1004 format ('dbkno4: inod,rfnormalflow,dbknon ',i3,1x,2(d14.7))
      DRB(J)=DRB(J)-rdum2  ! subtract normal flow across interval
      END IF
      END DO
      ENDIF
      IF (IDOMAINTYPE(ISTR).EQ.3) THEN
C
C     Closed slurry wall
C
      IF (RDBK(ISTR).EQ.0.0.OR.RDBW(ISTR).EQ.0.0) THEN ! test for head specified boundaries inside wall
        NHEADCONDITIONS=0
        DO I=1,N
          IF (ITYPE(I).EQ.1) THEN
            CALL DBPREP(CZI(I))
            IF (LFDBIN(ISTR)) THEN
             NHEADCONDITIONS=NHEADCONDITIONS+1
            END IF
          END IF
        END DO
      IF (NHEADCONDITIONS.EQ.0) THEN
      AMESS(1)='Problem in solution procedure when processing slurry'
      AMESS(2)='wall with label:'
      AMESS(3)=ADBLAB(INOD1)
      AMESS(4)='At least one head specified boundary condition is '
      AMESS(5)=
     &  'needed inside a closed slurry wall with zero conductance.'
      CALL HALT(5)  ! end program execution and produce fatal error message
      END IF
      END IF
      DO INOD=INOD1,INODL
      IF (INOD.EQ.INOD1) THEN
      INODM1=INODL
      ELSE
      INODM1=INOD-1
      END IF
      IF (INOD.EQ.INODL) THEN
      INODP1=INOD1
      ELSE
      INODP1=INOD+1
      END IF
      J=J+1      ! integration section across node INOD
      CZ=CZI(J)
      CZA=CALPH(J)
      CALL DBFILL (INOD1,INODM1,INODL,CZ)
      RXA=0.5*(REAL(CDBZEL(INODM1))+1.0)
      CALL DBFILL (INOD1,INOD,INODL,CZA)
      RXB=0.5*(REAL(CDBZEL(INOD))+1.0)
      CALL DBPREP(CZRESTART)   ! reset CZ0 in DBPREP because DBFILL changed data in common
      DRB(J)=DRB(J)+RDBTFACNOD(INOD)/RFPERM(CZ)*(
     &ABS(CDBZ21(INODM1))*(DBSTR(INODM1)*RFINTLABRIGHT(RXA,1.0)+
     &DBSTR(INOD)*RFINTLABLEFT(RXA,1.0)+
     &DBQSTR(INODM1)*RFINTMU(RXA,1.0))+
     &ABS(CDBZ21(INOD))*(DBSTR(INOD)*RFINTLABRIGHT(0.0,RXB)+
     &DBSTR(INODP1)*RFINTLABLEFT(0.0,RXB)+
     &DBQSTR(INOD)*RFINTMU(0.0,RXB))) ! correct for effect of old strength terms
c      rdum1=RFNORMALFLOW(CZ,CDBZ(INOD))+RFNORMALFLOW(CDBZ(INOD),CZA)
      rdum2=dbknon(inod)
c      write (ilume,1005) inod,rdum1,rdum2
c 1005 format ('dbkno5: inod,rfnormalflow,dbknon ',i3,1x,2(d14.7))
      DRB(J)=DRB(J)-rdum2  ! subtract normal flow across interval
      J=J+1     ! integration section on line doublet between INOD and INOD+1
      CZ=CZI(J)
      CZA=CALPH(J)
      CALL DBFILL (INOD1,INOD,INODL,CZ)
      RXA=0.5*(REAL(CDBZEL(INOD))+1.0)
      CALL DBFILL (INOD1,INOD,INODL,CZA)
      RXB=0.5*(REAL(CDBZEL(INOD))+1.0)
      CALL DBPREP(CZRESTART)  ! reset CZ0 in DBPREP because DBFILL changed data in common
      DRB(J)=DRB(J)+RDBTFACCTR(INOD)/RFPERM(CZ)*(
     &ABS(CDBZ21(INOD))*(DBSTR(INOD)*RFINTLABRIGHT(RXA,RXB)+
     &DBSTR(INODP1)*RFINTLABLEFT(RXA,RXB)+
     &DBQSTR(INOD)*RFINTMU(RXA,RXB)))        ! correct for the old strength term
c      rdum1=rfnormalflow(cz,cza) ! ********************************************** DEBUGGING
      rdum2=dbknoc(inod)
c      write (ilume,1006) inod,rdum1,rdum2
c 1006 format ('dbkno6: inod,rfnormalflow,dbknoc ',i3,1x,2(d14.7))
      DRB(J)=DRB(J)-rdum2         ! subtract normal flow across interval
      END DO
      END IF
  20  CONTINUE
      RETURN
 1100 format (' ***ERROR in DBKNO: collocation point ',2(d14.7),/,
     &' is not on shifted first collocation point of line-doublet ',i3,/
     &' of inhomogeneity domain ',i3)
 1200 format (' ***ERROR in DBKNO: collocation point ',2(d14.7),/,
     &' is not on starting node of line-doublet ',i3,/
     &' of inhomogeneity domain ',i3)
 1300 format (' ***ERROR in DBKNO: collocation point ',2(d14.7),/,
     &' is not on center collocation point of line-doublet ',i3,/
     &' of inhomogeneity domain ',i3)
 1400 format (' ***ERROR in DBKNO: collocation point ',2(d14.7),/,
     &' is not on shifted last collocation point of line-doublet ',i3,/
     &' of inhomogeneity domain ',i3)
      END
c
C -------------------------------------------------------------------------------------------------------
c
      SUBROUTINE DBSUB (DRB,J)
C
c -------------------------------------------------------------------------------------------------------
c
C     Routine substitutes the solution vector into DBSTR and DBQSTR.
C     The average strength ,DBAVS, is also computed or updated.
C
      IMPLICIT NONE
      LOGICAL lshift1,lshift2,ls,lc,le
      INTEGER(4) J,ISTR,INOD1,INODL,J0,INOD,INODM1,inodp1
      REAL(8) DRB,DAVADD
      INCLUDE 'dbcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      DIMENSION DRB(*)
      IF (NDB.EQ.0) RETURN
      j0=0
      DO 30 ISTR=1,NDBSTR
      INOD1=IDBSTA(ISTR)
      IF (ldbrechargeonly(istr)) GOTO 30 !  "recharge only", skip domain
      INODL=INOD1+NDBSTI(ISTR)-1
      IF (IDOMAINTYPE(ISTR).EQ.1) THEN
C
C     Hydraulic conductivity inhomogeneity
C
C --------------add new strengths
      DO INOD=INOD1,INODL
      if (idbcode(inod).gt.0) THEN ! line-doublet is indeed in the matrix
       inodm1=inod-1
       if (inod.eq.inod1) inodm1=inodl
       inodp1=inod+1
       if (inod.eq.inodl) inodp1=inod1
       lshift1=idbcode(inod).eq.2
       lshift2=idbcode(inodp1).eq.2.or.idbcode(inodp1).eq.0
       ls=rdbtfacnod(inod).ne.1.0d21
       lc=rdbtfacctr(inod).ne.1.0d21
       le=rdbtfacend(inod).ne.1.0d21
       if (ls) then
        j=j+1
        DBSTR(INOD)=DBSTR(INOD)+DRB(J)
        if (.not.lshift1) then
         dbsen(inodm1)=dbstr(inod)
        end if
       end if
       if (lc) then
        j=j+1
        DBQSTR(INOD)=DBQSTR(INOD)+DRB(J)
       end if
       if (lshift2.and.le) then
        j=j+1
        dbsen(inod)=dbsen(inod)+drb(j)
       end if
      endif
      END DO
      ENDIF
      IF (IDOMAINTYPE(ISTR).EQ.2) THEN
C
C     Open slurry wall
C
      DBSTR(INOD1)=0.0D0   ! first node has zero strength
      DO INOD=INOD1,INODL
      inodm1=inod-1
      IF (INOD.NE.INOD1) THEN
      J=J+1
      DBSTR(INOD)=DBSTR(INOD)+DRB(J)
      dbsen(inodm1)=dbstr(inod) ! no common boundaries, this is always true
      ENDIF
      J=J+1
      DBQSTR(INOD)=DBQSTR(INOD)+DRB(J)
      END DO
      DBSTR(INODL+1)=0.0D0   ! last node has zero strength
      dbsen(inodl)=0.0d0
      END IF
      IF (IDOMAINTYPE(ISTR).EQ.3) THEN
C
C     Closed slurry wall
C
      DO INOD=INOD1,INODL
      inodm1=inod-1
      if (inod.eq.inod1) inodm1=inodl
      J=J+1
      DBSTR(INOD)=DBSTR(INOD)+DRB(J) ! first strength for first node
      dbsen(inodm1)=dbstr(inod) ! no common boundaries, this is always true
      J=J+1
      DBQSTR(INOD)=DBQSTR(INOD)+DRB(J)
      END DO
      END IF
  30  CONTINUE
      RETURN
      END
c
c -------------------------------------------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFINTLABLEFT(RX1,RX2)
c
c -------------------------------------------------------------------------------------------------------------
c
C
C     Function returns the integral of the linear strength component along the line doublet.
C     Integration is performed from RX1 to RX2, between 0 and 1, with the strength 0 at RX=0 and
C     with the strength 1 at RX=1.
C
      IMPLICIT NONE
      REAL(8) RX1,RX2
C
      RFINTLABLEFT=0.5*(RX2*RX2-RX1*RX1)
      RETURN
      END
c
c -------------------------------------------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFINTLABRIGHT(RX1,RX2)
c
c -------------------------------------------------------------------------------------------------------------
c
C
C     Function returns the integral of the linear strength component along the line doublet.
C     Integration is performed from RX1 to RX2, between 0 and 1, with the strength 1 at RX=0 and
C     with the strength 0 at RX=1.
C
      IMPLICIT NONE
      REAL(8) RX1,RX2
C
      RFINTLABRIGHT=RX2-RX1-0.5*(RX2*RX2-RX1*RX1)
      RETURN
      END
C
C -------------------------------------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFINTMU(RX1,RX2)
c
c -------------------------------------------------------------------------------------------------------------
c
C
C     Function returns the integral of the quadratic strength component along the line doublet.
C     Integration is performed from RX1 to RX2, between 0 and 1, with the strength 0 at RX=0 and
C     RX=1, and the strength is 1 at RX=0.5
C
      IMPLICIT NONE
      REAL(8) RX1,RX2
C
      RFINTMU=(RX2*RX2*RX2-RX1*RX1*RX1)/3.0-(RX2*RX2-RX1*RX1)/2.0
      RFINTMU=RFINTMU*4.0
      RETURN
      END
C
C --------------------------------------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFDBNFCONDITION (ISTR,CZA,CZB)
c
c -------------------------------------------------------------------------------------------------------------
c
C
C     Function returns the normal flow condition for the section between CZA and CZB of string ISTR.
C     CZA and CZB must be on the string and the string must be a slurry wall.
C     The flow condition is the doublet strength (delta PHI) times the conductance of the wall divided by
C     the hydraulic conductivity of the aquifer.
C
      IMPLICIT NONE
      INTEGER(4) ISTR,INOD,INOD1,INODL,INODP1
      LOGICAL LA,LB
      REAL(8) RFPERM,RXA,RXB,RCONDUCT,
     &     RFINTLABLEFT,RFINTLABRIGHT,RFINTMU,rkcz
      COMPLEX(8) CZ,CZA,CZB
      INCLUDE 'main.inc'
      INCLUDE 'dbcom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      RFDBNFCONDITION=0.0
      IF (IDOMAINTYPE(ISTR).EQ.1) THEN
      WRITE (ILUER,1000) ISTR
      RETURN
      END IF
      INOD1=IDBSTA(ISTR)
      INODL=INOD1+NDBSTI(ISTR)-1
      IF (IDOMAINTYPE(ISTR).EQ.2) THEN
c
c     Open slurry wall
c
      DO INOD=INOD1,INODL   ! loop over line doublets
      CALL DBFILL (INOD1,INOD,0,CZA)
      LA=ABS(AIMAG(CDDBLN(INOD))).GT.3.0
      IF (LA) THEN  ! CZA on line doublet
        CZ=CZA
        RXA=0.5*(REAL(CDBZEL(INOD))+1.0)
        RCONDUCT=RDBTFACNOD(INOD+1)
      ELSE
        RXA=0.0     ! CZA not on line doublet, but if CZB on line doublet start integral at starting node
      ENDIF
      CALL DBFILL (INOD1,INOD,0,CZB)
      LB=ABS(AIMAG(CDDBLN(INOD))).GT.3.0
      IF (LB) THEN  ! CZB on line doublet
        CZ=CZB
        RXB=0.5*(REAL(CDBZEL(INOD))+1.0)
        RCONDUCT=RDBTFACNOD(INOD)
      ELSE
        RXB=1.0     ! CZB not on line doublet, but if CZA on line doublet end integral at ending node
      END IF
      IF (LA.AND.LB) RCONDUCT=RDBTFACCTR(INOD)
      IF (LA.OR.LB) THEN ! all or part of section on line doublet
      RFDBNFCONDITION=RFDBNFCONDITION+
     & RCONDUCT/RFPERM(CZ)*ABS(CDBZ21(INOD))*
     & (DBSTR(INOD)*RFINTLABRIGHT(RXA,RXB)+
     & DBSTR(INOD+1)*RFINTLABLEFT(RXA,RXB)+
     & DBQSTR(INOD)*RFINTMU(RXA,RXB))
      END IF
      END DO
      END IF
      IF (IDOMAINTYPE(ISTR).EQ.3) THEN
c
c     Closed slurry wall
c
      DO INOD=INOD1,INODL
      IF (INOD.EQ.INODL) THEN
      INODP1=INOD1
      ELSE
      INODP1=INOD+1
      END IF
      CALL DBFILL (INOD1,INOD,INODL,CZA)
      LA=ABS(AIMAG(CDDBLN(INOD))).GT.3.0
      IF (LA) THEN    ! CZA on line doublet
        CZ=CZA
        RXA=0.5*(REAL(CDBZEL(INOD))+1.0)
        RCONDUCT=RDBTFACNOD(INODP1)
      ELSE
        RXA=0.0       ! CZA not on line doublet, but if CZB is then start integration interval at node.
      ENDIF
      CALL DBFILL (INOD1,INOD,INODL,CZB)
      LB=ABS(AIMAG(CDDBLN(INOD))).GT.3.0
      IF (LB) THEN   ! CZB on line doublet
        CZ=CZB
        RXB=0.5*(REAL(CDBZEL(INOD))+1.0)
        RCONDUCT=RDBTFACNOD(INOD)
      ELSE
        RXB=1.0      ! CZB not on line doublet, but if CZA is then end integration interval at node
      END IF
      IF (LA.AND.LB) RCONDUCT=RDBTFACCTR(INOD)
      IF (LA.OR.LB) THEN ! all or part of section on line doublet
c      rkcz=rfperm(cz)
c      write (iluer,1001) rkcz
c 1001 format (' RFDBNFCONDITION1: rkcz=',d14.7)
      RFDBNFCONDITION=RFDBNFCONDITION+
     & RCONDUCT/RFPERM(CZ)*ABS(CDBZ21(INOD))*
     & (DBSTR(INOD)*RFINTLABRIGHT(RXA,RXB)+
     & DBSTR(INODP1)*RFINTLABLEFT(RXA,RXB)+
     & DBQSTR(INOD)*RFINTMU(RXA,RXB))
      END IF
      END DO
      END IF
      RETURN
 1000 FORMAT (' ***ERROR in RFDBNFCONDITION: string ',i3,
     &        ' is not a slurry wall')
      END
c
c ------------------------------------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFDBCONDUCTANCE (ISTR,CZ)
c
c ------------------------------------------------------------------------------------------------------
c
C     Function calculates the conductance of the slurry wall ISTR at point CZ
C     It is assumed that the slurry wall has constant properties, but the aquifer
C     bottom may vary, resulting in a varying gap between the bottom of the slurry wall
C     and the aquifer bottom.
C
      IMPLICIT NONE
      INTEGER(4) ISTR
      REAL(8) RFPERM,RFBASE,RFHGHT,RGAP,RHEIGHT,RWETWALL,
     &     rht,rftop,rkcz
      COMPLEX(8) CZ
      INCLUDE 'main.inc'
      INCLUDE 'dbcom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      IF (ABS(RDBT(ISTR)).EQ.9999.0) THEN ! inside or outside wall, must be zero conductance.
      RFDBCONDUCTANCE=0.0
      RETURN
      ENDIF
      IF (RDBW(ISTR).LT.1.0E-10) THEN  ! zero wall thickness, makes for infinite conductance!?
      WRITE (ILUER,1000) ISTR,RDBW(ISTR)   ! issue warning
      RFDBCONDUCTANCE=1.0E+21
      RETURN
      END IF
      IF (RDBB(ISTR).EQ.-9999.0D+0) THEN
        RGAP=0.0D+0
      ELSE
        RGAP=RDBB(ISTR)-RFBASE(CZ)
        RGAP=MAX(0.0,RGAP)   ! gap underneath slurry wall is zero when slurry wall ends at or below aquifer base
      END IF
      rkcz=rfperm(cz)
c      write (iluer,1001) rkcz
c 1001 format (' RFDBCONDUCTANCE1: rkcz=',d14.7)
      IF (NSOL.GT.1) THEN ! subsequent iterations
      RHEIGHT=RFHGHT(CZ)  ! head will be used as aquifer top for unconfined flow !! THIS MAKES FOR NON-LINEARITY
        IF (RHEIGHT.le.1.0E-10) THEN
          RFDBCONDUCTANCE=0.0
        ELSE
           RWETWALL=RHEIGHT-RGAP
           RFDBCONDUCTANCE=(RWETWALL*RDBK(ISTR)+RGAP*rkcz)/
     &                (RDBW(ISTR)*RHEIGHT)
c Note: if rgap=0 the term rheight cancels and the conductance does not vary between iteration: linear equation
c       if rgap>0 the term rheight does not cancel. For unconfined flow this means that the conductance
c       may change between iterations (because the head, thus rheight changes). This results in a non-linear equation
        END IF
      ELSE               ! first iteration, no solution yet.
        if (rgap.gt.0.0) then
          rht=rftop(cz)-rfbase(cz)  ! initially use aquifer thickness for heigth of wall
          rfdbconductance=((rht-rgap)*rdbk(istr)+rgap*rkcz)/
     &    (rdbw(istr)*rht)
        else
          RFDBCONDUCTANCE=RDBK(ISTR)/RDBW(ISTR)
        endif
      END IF
      RETURN
 1000 FORMAT (' ***ERROR in RFDBCONDUCTANCE: slurry wall width in ',
     &'of string ',I5,/,' is too small: width= ',E14.7)
      END
c
c -------------------------------------------------------------------------------------------------
c
      LOGICAL function ldbbottomOFF ()
c
c     True if inhomogeneities are found with a bottom jump AND an estimated average head. However, false when
c     the interface position that belongs with the estimated head is below the aquifer base.
c     The jump in bottom elevation is replaced by a jump in hydraulic conductivity for solution stability.
c     The logical function ldbbottomON will restore the bottom jump after which some extra iterations are made.
c
      implicit none
      INTEGER i,istr,inod1,inodl,inod
      LOGICAL lindom,lfinterface,lInterfaceExempt
      REAL(8) rki,rko,rbi,rbo,rp,rt,rh1,rhss,rfac1,rfac2,rsgs,rsgf
      COMPLEX(8) c1,c2,cz
      INCLUDE 'dbcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'main.inc'
c
c      write (iluer,1002)
c 1002 format (' Entering LDBBOTTOMOFF')
      ldbbottomOFF=.false.
c
c
      RETURN ! testing the exclusion of this logic.
c
c     Use head estimate, but do not replaced bottom jumps by k jumps.
c
c

      DO istr=1,ndbstr
      IF (IDOMAINTYPE(ISTR).EQ.1) then !    ---------- inhomogeneity domains ----------------
        IF (LINDOM (ISTR,RKI,RKO,RP,RBI,RBO,RT)) THEN
          rdbktemp(istr)=-9999.0  ! flag all inhomogeneties that are not altered
          if (nsol.eq.0.and.rbi.ne.rbo.and.rdbh(istr).NE.-9999.0) THEN ! found bottom jump with estimated head
            if (lfinterface()) then
              call interfacedata (rhss,rsgs,rsgf,rfac1,rfac2)
              lInterfaceExempt=rdbh(istr).lt.rfac1*(rhss-rbi)+rbi ! interface is above aquifer base
              lInterfaceExempt=.true. ! do not use conditional logic, since head estimate may be too high
            else
              lInterfaceExempt=.false.
            end if
            if (.not.lInterfaceExempt) then ! do not use k-jump when interface above aquifer base
              rdbbtemp(istr)=rdbb(istr)  ! save original aquifer bottom setting
              rdbktemp(istr)=rdbk(istr)  ! save original aquifer conductivity setting
              rdbb(istr)=rbo ! remove bottom jump
              rh1=rt
              rh1=MIN(rh1,rdbh(istr))
              rdbk(istr)=rki*(rh1-rbi)/(rh1-rbo) ! replace by suitable k jump
              ldbbottomOFF=.true.
c              write (iluer,1001) istr,rdbbtemp(istr),rdbktemp(istr),
c     &                           rdbb(istr),rdbk(istr),rdbh(istr)
c 1001 format(' LDBBOTTOMOFF: found bottom inhomogeneity with estimated',
c     &        ' head and calculated equivalent k jump.',/,' istr, '
c     & 'rdbbtemp(istr),rdbktemp(istr),rdbb(istr),rdbk(istr),rdbh(istr)',
c     &  i5,5(d14.7))
            endif
          endif
        endif
      endif
      end do
C
c     set properties
c
      DO ISTR=1,NDBSTR
       IF (IDOMAINTYPE(ISTR).EQ.1) THEN
        INOD1=IDBSTA(ISTR)
        INODL=INOD1+NDBSTI(ISTR)-1
        DO INOD=INOD1,INODL
         CALL SET_DBPROPERTIES (ISTR,INOD) ! store inside and outside k and b values
        END DO
       END IF
      END DO
      return
      end
c
c -------------------------------------------------------------------------------------------------
c
      LOGICAL function ldbbottomON ()
c
c     True if inhomogeneities are found with a bottom jump AND an estimated average head.
c     The jump in bottom elevation was replaced by a jump in hydraulic conductivity in the
c     logical function ldbbottomOFF.
c     The inside bottom and inside conductivity are now restored to the original settings.
c
      implicit none
      INTEGER i,istr,inod1,inodl,inod
      REAL(8) rki,rko,rbi,rbo,rp,rt
      COMPLEX(8) c1,c2,cz
      INCLUDE 'dbcom.inc'
      INCLUDE 'lusys.inc'
c
c      write (iluer,1002)
c 1002 format (' Entering LDBBOTTOMON')
      ldbbottomON=.false.
c
c
      RETURN ! testing the exclusion of this logic.
c
c     Use head estimate, but bottom jumps not replaced by k jumps.
c
c
      DO istr=1,ndbstr
      IF (IDOMAINTYPE(ISTR).EQ.1) then !    ---------- inhomogeneity domains ----------------
        if (rdbktemp(istr).NE.-9999.0) then ! found modified bottom inhomogeneity, restore
          rdbb(istr)=rdbbtemp(istr)
          rdbk(istr)=rdbktemp(istr)
          ldbbottomON=.true.
c          write (iluer,1001) istr,rdbbtemp(istr),rdbktemp(istr),
c     &                       rdbb(istr),rdbk(istr),rdbh(istr)
c 1001 format(' LDBBOTTOMON: found bottom inhomogeneity with estimated',
c     &        ' head and restored bottom jump.',/,' istr, '
c     & 'rdbbtemp(istr),rdbktemp(istr),rdbb(istr),rdbk(istr),rdbh(istr)',
c     &  i5,5(d14.7))
        end if
      endif
      end do
C
c     set properties
c
      DO ISTR=1,NDBSTR
       IF (IDOMAINTYPE(ISTR).EQ.1) THEN
        INOD1=IDBSTA(ISTR)
        INODL=INOD1+NDBSTI(ISTR)-1
        DO INOD=INOD1,INODL
         CALL SET_DBPROPERTIES (ISTR,INOD) ! store inside and outside k and b values
        END DO
       END IF
      END DO
      return
      end
c
c ------------------------------------------------------------------------------------------
c
      subroutine dbupdate(drscr,j,m,n,czi,calph,lsubinclude,isubsign)
c
c ------------------------------------------------------------------------------------------
c
c     Routine substitutes potentials or integrated flows in known vector arrays for line-doublets
c
c
c     NOTE: make sure that DRA contains uncorrected matrix coefficients
c           make sure that DRB contains the latest solution vector
c
c
      implicit none
      INTEGER j,m,n,nsolOut,istr,inod,inod1,inodl,inodm1,
     &        icount,inods,inode,inodp1,isubsign
      LOGICAL lfirst,lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut,
     &        lshift1,lshift2,ls,lc,le,
     &        lsubinclude,laddsubcells
      REAL(8) drscr,rfpot,rfnormalflow,rdum1,rdum2,sloc,rbigx,
     &        rt,rbo,rbi,rp,rko,rki,rfhead,rfhedp,rhloc,rftop,rfhfp,
     &        rfdbslocal,rpot,rfnflksub,rsubsign
      COMPLEX(8) czi,calph,cz,czz,cza,cfsmallz,cflk_subomega
      CHARACTER(8) aBasenameOut
      CHARACTER(16)aDateTimeOut
      DIMENSION drscr(m),czi(m),calph(m)
      INCLUDE 'dbcom.inc'
      INCLUDE 'lusys.inc'
c
      if (ndb.EQ.0) return  ! no line-doublets present
      laddsubcells=lsubinclude
      rsubsign=REAL(isubsign)
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      if (nsolOut.eq.0.) then
c
c      -------------------------------------------
c      first solution, call RFPOT and RFNORMALFLOW
c      -------------------------------------------
c
c      Note: these values are not yet meaningful (no strength parameters calculated yet),
c            but will be build on by successive solutions.
c      Note: the heads rdbhnod, rdbhctr and rdbhend (for use in DBKAC are set in DBcollocation_PREP
c
c      write (ilume,1001) nsolOut
c 1001 format (' dbupdate1: nsolOut=',i3,
c     &        ' set arrays to rfpot or rfnormalflow.')
      DO 20 ISTR=1,NDBSTR
      INOD1=IDBSTA(ISTR)
      INODL=INOD1+NDBSTI(ISTR)-1  ! last node of closed string, one but last node of open string
      IF (IDOMAINTYPE(ISTR).EQ.1) THEN
      IF (LDBRECHARGEONLY(ISTR)) GOTO 20 !  "recharge only", skip domain
c
c      Hydraulic conductivity inhomogeneity
c
        DO INOD=INOD1,INODL    ! calculate inside potentials at collocation points
        if (idbcode(inod).gt.0) then ! line-doublet is indeed in the matrix
         inodp1=inod+1
         if (inod.eq.inodl) inodp1=inod1
         lshift1=idbcode(inod).eq.2
         lshift2=idbcode(inodp1).eq.2.or.idbcode(inodp1).eq.0
         ls=rdbtfacnod(inod).ne.1.0d21
         lc=rdbtfacctr(inod).ne.1.0d21
         le=rdbtfacend(inod).ne.1.0d21
         if (ls) then
          j=j+1 ! first collocation point
          cz=cdbz(inod)
          if (lshift1) then ! collocation point shifted away from first node
           cz=CMPLX(-1.0d0+rdboffsetinhom,0.0d0)
           cz=cfsmallz(cz,cdbz(inod),cdbz(inodp1))
          end if
          dbknon(inod)=rfpot(cz)
         endif
         if (lc) then
          j=j+1 ! second (center) collocation point
          CZ=0.5*(CDBZ(INOD)+CDBZ(INODP1))
          dbknoc(inod)=rfpot(cz)
         endif
         if (lshift2.and.le) then  ! there is also a collocation point near the end
           j=j+1
           cz=CMPLX(1.0d0-rdboffsetinhom,0.0d0)
           cz=cfsmallz(cz,cdbz(inod),cdbz(inodp1))
           dbknoe(inod)=rfpot(cz)
         end if
        endif
        END DO
      ENDIF
      IF (IDOMAINTYPE(ISTR).EQ.2) THEN
c
c      Open slurry wall
c
      DO INOD=INOD1,INODL
      J=J+1
      CZ=CZI(J)
      CZA=CALPH(J)                 ! integration interval on line doublet between nodes INOD and INODP1
      dbknoc(inod)=RFNORMALFLOW(CZ,CZA)         ! normal flow across interval
      IF (INOD.NE.INODL) THEN
      J=J+1
      CZ=CZI(J)
      CZA=CALPH(J)                 ! integration interval across node INOD+1
      dbknon(inod)=RFNORMALFLOW(CZ,CDBZ(INOD+1))+
     &            RFNORMALFLOW(CDBZ(INOD+1),CZA)  ! normal flow across interval
      END IF
      END DO
      END IF
      IF (IDOMAINTYPE(ISTR).EQ.3) THEN
c
c      Closed slurry wall     (rdboffsetn=0.5 and rdboffsetc=0.05)
c                             these are set in Block Data DBDAT
C                             current settings create overlapping integration segments
c
      DO INOD=INOD1,INODL
      J=J+1      ! integration section across node INOD
      CZ=CZI(J)
      CZA=CALPH(J)
      dbknon(inod)=RFNORMALFLOW(CZ,CDBZ(INOD))+
     &RFNORMALFLOW(CDBZ(INOD),CZA)
c
      J=J+1     ! integration section on line doublet between INOD and INOD+1
      CZ=CZI(J)
      CZA=CALPH(J)
      dbknoc(inod)=rfnormalflow(cz,cza)
      END DO
      END IF
  20  continue
c
      else
c
c     -----------------------------------------------------
c     successive iterations, update the known vector arrays
c     -----------------------------------------------------
c
c      write (ilume,1002) nsolOut
c 1002 format (' dbupdate2: nsolOut=',i3,
c     &        ' update arrays drscr from Aij*Sj.')
      DO 30 ISTR=1,NDBSTR
      INOD1=IDBSTA(ISTR)
      INODL=INOD1+NDBSTI(ISTR)-1  ! last node of closed string, one but last node of open string
      IF (IDOMAINTYPE(ISTR).EQ.1) THEN
      IF (ldbrechargeonly(istr)) GOTO 30 ! "recharge only", skip domain
c
c      Hydraulic conductivity inhomogeneity
c
        icount=0
        rdbpotaverage(istr)=0.0D0
        DO INOD=INOD1,INODL
        if (idbcode(inod).gt.0) THEN ! skip if not in matrix
         inodp1=inod+1
         if (inod.eq.inodl) inodp1=inod1
         lshift1=idbcode(inod).eq.2
         lshift2=idbcode(inodp1).eq.2.or.idbcode(inodp1).eq.0
         ls=rdbtfacnod(inod).ne.1.0d21
         lc=rdbtfacctr(inod).ne.1.0d21
         le=rdbtfacend(inod).ne.1.0d21
         rbo=rdbbo(inod)
         rbi=rdbbi(inod)
         rko=rdbko(inod)
         if (ls) THEN ! collocation point at start of line-doublet
          j=j+1           ! collocation point at or near first node
          cz=czi(j)
          dbknon(inod)=dbknon(inod)+drscr(j)
          if (laddsubcells) then
            dbknon(inod)=dbknon(inod)+rsubsign*REAL(cflk_subomega(cz))
          end if
c      rpot=rfpot(cz) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! for debugging only
c      write (iluer,1001) inod,rpot,dbknon(inod)
c 1001 format (' dbupdate1: inod,rpot,dbknon(inod)',i5,2(1x,d14.7))
           if (rbi.le.rbo) then
           rdberrstrt(inod)=dbknon(inod) ! store inside potential for use in RFDBER
           rdbhnod(inod)=rfhedp(rdberrstrt(inod),cz) ! store inside head for use in DBKFAC
          else
           if (lshift1) then
            rbigx=-1.0d0+rdboffsetinhom
            sloc=rfdbslocal (rbigx,dbstr(inod),dbqstr(inod),dbsen(inod))
           else
            sloc=dbstr(inod)
           end if
           rdberrstrt(inod)=dbknon(inod)-sloc ! store outside potential for use in RFDBER
           rhloc=rftop(cz)-rbo
           rdbhnod(inod)=rfhfp(rdberrstrt(inod),rko,rhloc,rbo)+rbo ! store outside head for use in DBKFAC
          end if
          icount=icount+1
         endif
         if (lc) then  ! collocation point at center of line-doublet
          j=j+1          ! collocation point at the center
          cz=czi(j)
          dbknoc(inod)=dbknoc(inod)+drscr(j)
          if (laddsubcells) then
            dbknoc(inod)=dbknoc(inod)+rsubsign*REAL(cflk_subomega(cz))
          end if
c      rpot=rfpot(cz) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! for debugging only
c      write (iluer,1002) inod,rpot,dbknoc(inod)
c 1002 format (' dbupdate2: inod,rpot,dbknoc(inod)',i5,2(1x,d14.7))
          if (rbi.le.rbo) then
           rdberrcntr(inod)=dbknoc(inod) ! store inside potential for use in RFDBER
           rdbhctr(inod)=rfhedp(rdberrcntr(inod),cz) ! store inside head for use in DBKFAC
          else
           rbigx=0.0d0
           sloc=rfdbslocal (rbigx,dbstr(inod),dbqstr(inod),dbsen(inod))
           rdberrcntr(inod)=dbknoc(inod)-sloc ! store outside potential for use in RFDBER
           rhloc=rftop(cz)-rbo
           rdbhctr(inod)=rfhfp(rdberrcntr(inod),rko,rhloc,rbo)+rbo ! store outside head for use in DBKFAC
          end if
          icount=icount+1
         endif
         rdberrend(inod)=0.0d0
         if (lshift2.and.le) then ! there is also a collocation point near the line-doublet end
          j=j+1
          cz=czi(j)
          dbknoe(inod)=dbknoe(inod)+drscr(j)
          if (laddsubcells) then
            dbknoe(inod)=dbknoe(inod)+rsubsign*REAL(cflk_subomega(cz))
          end if
c      rpot=rfpot(cz) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! for debugging only
c      write (iluer,1003) inod,rpot,dbknoe(inod)
c 1003 format (' dbupdate3: inod,rpot,dbknoe(inod)',i5,2(1x,d14.7))
          if (rbi.le.rbo) then
           rdberrend(inod)=dbknoe(inod) ! store inside potential for use in RFDBER
           rdbhend(inod)=rfhedp(rdberrend(inod),cz) ! store inside head for use in DBKFAC
          else
           rbigx=1.0d0-rdboffsetinhom
           sloc=rfdbslocal (rbigx,dbstr(inod),dbqstr(inod),dbsen(inod))
           rdberrend(inod)=dbknoe(inod)-sloc ! store outside potential for use in RFDBER
           rhloc=rftop(cz)-rbo
           rdbhend(inod)=rfhfp(rdberrend(inod),rko,rhloc,rbo)+rbo ! store outside head for use in DBKFAC
          end if
          icount=icount+1
         end if
         rdbpotaverage(istr)=rdbpotaverage(istr)+
     &                 rdberrstrt(inod)+rdberrcntr(inod)+rdberrend(inod)
        endif
        END DO
        if (icount.ne.0) then
         rdbpotaverage(istr)=rdbpotaverage(istr)/icount ! calculate average potential for use in RFDBER
        end if
      ENDIF
      IF (IDOMAINTYPE(ISTR).EQ.2) THEN
c
c      Open slurry wall
c
        DO INOD=INOD1,INODL
          j=j+1
          cz=czi(j)
          cza=calph(j)
          dbknoc(inod)=dbknoc(inod)+drscr(j)
          if (laddsubcells) then
            dbknoc(inod)=dbknoc(inod)+rsubsign*rfnflksub(cz,cza) ! sub-cell contributions not in drscr(j)
          end if
        if (inod.ne.inodl) then
          j=j+1
          cz=czi(j)
          cza=calph(j)
          dbknon(inod)=dbknon(inod)+drscr(j)
          if (laddsubcells) then
            dbknon(inod)=dbknon(inod)+rsubsign*rfnflksub(cz,cza) ! sub-cell contributions not in drscr(j)
          end if
        END IF
        END DO
      END IF
      IF (IDOMAINTYPE(ISTR).EQ.3) THEN
c
c      Closed slurry wall     (rdboffsetn=0.5 and rdboffsetc=0.05)
c                             these are set in Block Data DBDAT
C                             current settings create overlapping integration segments
c
        DO INOD=INOD1,INODL
        j=j+1
        cz=czi(j)
        cza=calph(j)
        dbknon(inod)=dbknon(inod)+drscr(j)
        if (laddsubcells) then
          dbknon(inod)=dbknon(inod)+rsubsign*rfnflksub(cz,cza) ! sub-cell contributions not in drscr(j)
        end if
        j=j+1
        cz=czi(j)
        cza=calph(j)
        dbknoc(inod)=dbknoc(inod)+drscr(j)
        if (laddsubcells) then
          dbknoc(inod)=dbknoc(inod)+rsubsign*rfnflksub(cz,cza) ! sub-cell contributions not in drscr(j)
        end if
        END DO
      END IF
  30  continue
      endif
      return
c
      end subroutine
c
C ----------------------------------------------------------------------------------------------------------
c
      SUBROUTINE dbupdate_check (j,m,n,czi,calph,itype)
c
c ----------------------------------------------------------------------------------------------------------
C
C     Routine compares the boundary value correction terms with rfpot or rfnormalflow
C      
      IMPLICIT NONE
      INTEGER(4) J,m,N,ISTR,INOD1,INODL,INOD,INODP1,INODM1,I,
     &           NHEADCONDITIONS,itype
      LOGICAL LFDBIN,ls,lc,le,lshift2
      REAL(8) RFPOT,RXA,RXB,RFINTLABLEFT,RFINTLABRIGHT,RFINTMU,
     &        RFPERM,RFNORMALFLOW,
     &        rdum1,rdum2
      COMPLEX(8) CZI,CALPH,CZ,CZA
      INCLUDE 'dbcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      DIMENSION CZI(m),CALPH(m),itype(m)
      IF (NDB.EQ.0) RETURN
      DO 20 ISTR=1,NDBSTR
      INOD1=IDBSTA(ISTR)
      IF (ldbrechargeonly(istr)) GOTO 20 !  "recharge only", skip domain
      INODL=INOD1+NDBSTI(ISTR)-1
      IF (IDOMAINTYPE(ISTR).EQ.1) THEN
c
c     Hydraulic conductivity inhomogeneity
c
      DO INOD=INOD1,INODL
      IF (idbcode(inod).gt.0) THEN ! skip if not in matrix
       inodp1=inod+1
       if (inod.eq.inodl) inodp1=inod1
       lshift2=idbcode(inodp1).eq.2.or.idbcode(inodp1).eq.0
       ls=rdbtfacnod(inod).ne.1.0d21
       lc=rdbtfacctr(inod).ne.1.0d21
       le=rdbtfacend(inod).ne.1.0d21
       if (ls) then
        J=J+1
        CZ=CZI(J)
        rdum1=rfpot(cz)
        rdum2=dbknon(inod)
        write (ilume,1001) j,rdum1,rdum2
 1001 format ('dbupdate_check1: j,rfpot,dbknon ',i4,2x,2(d14.7))
       end if
       if (ls) then
        J=J+1
        CZ=CZI(J)
        rdum1=rfpot(cz)
        rdum2=dbknoc(inod)
        write (ilume,1002) j,rdum1,rdum2
 1002 format ('dbupdate_check2: j,rfpot,dbknoc ',i4,2x,2(d14.7))
       end if
       if (le.and.lshift2) then
        J=J+1
        CZ=CZI(J)
        rdum1=rfpot(cz)
        rdum2=dbknoe(inod)
        write (ilume,1003) j,rdum1,rdum2
 1003 format ('dbupdate_check3: j,rfpot,dbknoe ',i4,2x,2(d14.7))
       end if
      END IF
      END DO
      ENDIF
      IF (IDOMAINTYPE(ISTR).EQ.2) THEN
C
C     Open slurry wall
C
      DO INOD=INOD1,INODL
      J=J+1
      CZ=CZI(J)
      CZA=CALPH(J)                 ! integration interval on line doublet between nodes INOD and INODP1
      rdum1=rfnormalflow(cz,cza)
      rdum2=dbknoc(inod)
      write (ilume,1004) j,rdum1,rdum2
 1004 format ('dbupdate_check4: j,rfnormalflow,dbknoc ',
     &          i4,2x,2(d14.7))
      IF (INOD.NE.INODL) THEN
      J=J+1
      CZ=CZI(J)
      CZA=CALPH(J)                 ! integration interval across node INOD+1
      rdum1=RFNORMALFLOW(CZ,CDBZ(INOD+1))+RFNORMALFLOW(CDBZ(INOD+1),CZA)
      rdum2=dbknon(inod)
      write (ilume,1005) j,rdum1,rdum2
 1005 format ('dbupdate_check5: j,rfnormalflow,dbknon ',
     &         i4,2x,2(d14.7))
      END IF
      END DO
      ENDIF
      IF (IDOMAINTYPE(ISTR).EQ.3) THEN
C
C     Closed slurry wall
C
      DO INOD=INOD1,INODL
      J=J+1      ! integration section across node INOD
      CZ=CZI(J)
      CZA=CALPH(J)
      rdum1=RFNORMALFLOW(CZ,CDBZ(INOD))+RFNORMALFLOW(CDBZ(INOD),CZA)
      rdum2=dbknon(inod)
      write (ilume,1006) j,rdum1,rdum2
 1006 format ('dbupdate_check6: j,rfnormalflow,dbknon ',
     &         i4,2x,2(d14.7))
c
      J=J+1     ! integration section on line doublet between INOD and INOD+1
      CZ=CZI(J)
      CZA=CALPH(J)
      rdum1=rfnormalflow(cz,cza)
      rdum2=dbknoc(inod)
      write (ilume,1007) j,rdum1,rdum2
 1007 format ('dbupdate_check7: j,rfnormalflow,dbknoc ',
     &          i4,2x,2(d14.7))
      END DO
      END IF
  20  CONTINUE
      RETURN
      END
c
c-----------------------------------------------------------------------------------------------
c
      subroutine dbgentbfac (czcol,rhcol,rki,rko,rbi,rbo,rfcol,rgcol)
c
c-----------------------------------------------------------------------------------------------
c
c    Routine generates the matrix and knownfactor correction factors at control point CZCOL
c    Routine is called in DBFAC
c
c    Input:
c
c    czcol    collocation point on line-doublet
c    rhcol    head at collocation point
c    rki,rko  hydraulic conductivity inside and outside of domain to which line-doublet belongs
c    rbi,rbo  aquifer base elevation inside and outside of domain to which line-doublet belongs
c
c    Output:
c
c    rfcol    transmissivity factor used in matrix correction and known vector routines
c    rgcol    base jump factor used in matrix correction and known vector routines
c
c    Note: For readability there is some repetition of logic!
c
      IMPLICIT NONE
      INTEGER(4) InterfaceCase,nsolOut
      LOGICAL lfinterface,
     &        lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut
      REAL(8) rhit,rftop,rhcol,rki,rko,rbi,rbo,rfcol,rgcol,
     &        rtop,rhloci,rhloco,rh1,rhi1,rho1,
     &        rti1,rto1,rg1,rf1
      REAL(8) rzInterfaceElevation,rfinterface,rCInterfaceI,
     &        rCInterfaceO,rhss,rsgs,rsgf,rfac1,rfac2
      CHARACTER*16 aBasenameOut,aDateTimeOut
      COMPLEX(8) cz,czcol
      INCLUDE 'dbcom.inc'
      INCLUDE 'lusys.inc'
c         write (iluer,1001) czcol,rhcol,rki,rko,rbi,rbo,rfcol,rgcol
c 1001 format (' dbgentbfac1: czcol,rhcol,rki,rko,rbi,rbo,rfcol,rgcol',/,
c     &        2(d14.7),2x,7(d14.7))

C
                cz=czcol
C
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      if (lfinterface().and.nsolout.gt.0) then
        IF ((RKI.NE.RKO.OR.RBI.NE.RBO)) THEN
           if (rbi.ne.rbo) then
                rtop=rftop(cz)
                call interfacedata (rhss,rsgs,rsgf,rfac1,rfac2)
c         write (iluer,1002) rhss,rsgs,rsgf,rfac1,rfac2
c 1002 format (' dbgentbfac2: rhss,rsgs,rsgf,rfac1,rfac2',/,
c     &        5(d14.7))
                rhloci=rtop-rbi
                rhloco=rtop-rbo
                rh1=rhcol
                InterfaceCase=0 ! must be reset below
c
c               InterfaceCase = 0    Error, no case selected
c               InterfaceCase = 1    Interface below lowest aquifer base (use regular logic)
c               InterfaceCase = 2    Interface above highest base (submerged base jump)
c               InterfaceCase = 3    Interface is below outside base, but above inside base
c               InterfaceCase = 4    Interface is below inside base, but above outside base
c
                if (rbi.lt.rbo) then  ! base jumps down
                  rzInterfaceElevation=
     &            rfac2*(rhss-rbi)-rfac2/rfac1*(rh1-rbi)+rbi
                  if (rzInterfaceElevation.le.rbi) InterfaceCase=1
                  if (rzInterfaceElevation.ge.rbo) InterfaceCase=2
                  if (rzInterfaceElevation.gt.rbi.and.
     &                rzInterfaceElevation.lt.rbo) InterfaceCase=3
                else                  ! base jumps up
                  rzInterfaceElevation=
     &            rfac2*(rhss-rbo)-rfac2/rfac1*(rh1-rbo)+rbo
                  if (rzInterfaceElevation.le.rbo) InterfaceCase=1
                  if (rzInterfaceElevation.ge.rbi) InterfaceCase=2
                  if (rzInterfaceElevation.gt.rbo.and.
     &                rzInterfaceElevation.lt.rbi) InterfaceCase=4
                end if
c        write (iluer,1003) InterfaceCase,rzInterfaceElevation
c 1003   format (' DBGENTBFAC3: InterfaceCase=',i3,/,
c     &          '              Interface elevation=',d14.7)
c
                select case (InterfaceCase)
c
                  case (0)   ! error, case not selected
                  write (iluer,3000)
c
                  case (1)   ! no interface, add no delta C
                  rhi1=rh1-rbi
                  rhi1=MIN(rhi1,rhloci)
                  RHI1=MAX(rhi1,0.0d0)  ! avoid negative aquifer thickness
                  rho1=rh1-rbo
                  rho1=MIN(rho1,rhloco)
                  RHO1=MAX(rho1,0.0d0)  ! avoid negative aquifer thickness
                  RG1=0.5d0*rko*RHO1*(rbo-rbi)
c
                  case (2)   ! submerged interface, use k factor and add delta C across base jump
                  RHI1=1.0d0
                  RHO1=1.0d0
                if (rh1.lt.rtop) then
                  rCInterfaceI=0.5d0*rki*rfac1*(rhss-rbi)*(rhss-rbi)
                  rCInterfaceO=0.5d0*rko*rfac1*(rhss-rbo)*(rhss-rbo)
                else
                  rCInterfaceI=rki*rfac1*(rtop-rbi)*(rhss-rbi)-
     &                         0.5d0*rki*rfac1*(rtop-rbi)*(rtop-rbi)
                  rCInterfaceO=rko*rfac1*(rtop-rbo)*(rhss-rbo)-
     &                         0.5d0*rko*rfac1*(rtop-rbo)*(rtop-rbo)
                end if
                  RG1=rko/rki*rCInterfaceI-rCInterfaceO
c
                  case (3)  ! interface inside, T factor for above interface part of jump, add delta C for inside
                if (rh1.lt.rtop) then
                 RHI1=rh1-rzInterfaceElevation
                 RHO1=rh1-rbo
                 RHI1=MAX(RHI1,0.0d0)
                 RHO1=MAX(RHO1,0.0d0)
                  rCInterfaceI=0.5d0*rki*rfac1*(rhss-rbi)*(rhss-rbi)
                rCInterfaceO=0.5d0*rki*rfac1*(rhss-rzInterfaceElevation)
     &                       *(rhss-rzInterfaceElevation)
                else
                 RHI1=rtop-rzInterfaceElevation
                 RHO1=rtop-rbo
                  rCInterfaceI=rki*rfac1*(rtop-rbi)*(rhss-rbi)-
     &                         0.5d0*rki*rfac1*(rtop-rbi)*(rtop-rbi)
                  rCInterfaceO=rki*rfac1*(rtop-rzInterfaceElevation)
     &                         *(rhss-rzInterfaceElevation)-
     &                      0.5d0*rki*rfac1*(rtop-rzInterfaceElevation)
     &                         *(rtop-rzInterfaceElevation)
                end if
                  RG1=0.5d0*rko*rho1*(rbo-rzInterfaceElevation)+
     &               rko*rho1/(rki*rhi1)*(rCInterfaceI-rCInterfaceO)
c
                  case (4)  ! interface outside, add delta C for outside
                if (rh1.lt.rtop) then
                 RHI1=rh1-rbi
                 RHO1=rh1-rzInterfaceElevation
                 RHI1=MAX(RHI1,0.0d0)
                 RHO1=MAX(RHO1,0.0d0)
                rCInterfaceI=0.5d0*rko*rfac1*(rhss-rzInterfaceElevation)
     &                         *(rhss-rzInterfaceElevation)
                  rCInterfaceO=0.5d0*rko*rfac1*(rhss-rbo)*(rhss-rbo)
                else
                 RHI1=rtop-rbi
                 RHO1=rtop-rzInterfaceElevation
                  rCInterfaceI=rko*rfac1*(rtop-rzInterfaceElevation)
     &                         *(rhss-rzInterfaceElevation)-
     &                      0.5d0*rko*rfac1*(rtop-rzInterfaceElevation)
     &                         *(rtop-rzInterfaceElevation)
                  rCInterfaceO=rko*rfac1*(rtop-rbo)*(rhss-rbo)-
     &                         0.5d0*rko*rfac1*(rtop-rbo)*(rtop-rbo)
                end if
                  RG1=0.5d0*rko*rho1*(rzInterfaceElevation-rbi)+
     &                rCInterfaceI-rCInterfaceO
                end select
          else
              RHI1=1.0D0  ! No aquifer base jump, generate a conductivity factor only
              RHO1=1.0D0
          endif
            RTI1=RKI*RHI1
            RTO1=RKO*RHO1
            if (rti1.ne.rto1) then
              RF1=RTI1/(RTI1-RTO1)
            else               ! handle singularity, make rf1 very big
              rf1=rti1/1.0d-10 ! do not flag as zero jump, since the RG1 and RG2 term may still be non-zero.
            end if
        ELSE
            RF1=1.0D21  ! flag for inside T = outside T, or inside k < 0 (flag for recharge inhom. only)
            RG1=0.0D0
        ENDIF
c
      else   ! --------------------------------  original logic without interface flow --------------
          IF ((RKI.NE.RKO.OR.RBI.NE.RBO)) THEN
            IF (RBI.NE.RBO) THEN  ! Jump in aquifer base, figure jump in aquifer thickness.
                rtop=rftop(cz)
                rhloci=rtop-rbi
                rhloco=rtop-rbo
                rh1=rhcol
                rhi1=rh1-rbi
                rhi1=MIN(rhi1,rhloci) ! determine saturated aquifer thickness on the inside
                RHI1=MAX(RHI1,0.0D0)   ! prevent negative saturated thickness
                rho1=rh1-rbo
                rho1=MIN(rho1,rhloco) ! determine saturated aquifer thicknes on the outside
                RHO1=MAX(RHO1,0.0D0)   ! prevent negative saturated thickness
            ELSE
              RHI1=1.0D0  ! No aquifer base jump, generate a conductivity factor only
              RHO1=1.0D0
            ENDIF
            RTI1=RKI*RHI1
            RTO1=RKO*RHO1
            if (rti1.ne.rto1) then
              RF1=RTI1/(RTI1-RTO1)
            else               ! handle singularity, make rf1 very big
              rf1=rti1/1.0d-10 ! do not flag as zero jump, since the RG1 and RG2 term may still be non-zero.
            end if
            RG1=RTO1*0.5*(RBO-RBI)
          ELSE
            RF1=1.0D21  ! flag for inside T = outside T, or inside k < 0 (flag for recharge inhom. only)
            RG1=0.0D0
          ENDIF
      endif
c
          RFCOL=RF1
          RGCOL=RG1
c
      RETURN
 3000 format (' ***ERROR in DBGENTBFAC: ',
     & 'case for interface flow not selected.')
      end subroutine

