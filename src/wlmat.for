C     Last change:  HMH   4 Apr 2007    4:40 pm
C    This file contains the following routines and functions:
C
c    wlmatsize         update the matrix array dimensions
C    WLCZC             generates control point vector for wells
C    WLMAT             generates matrix coefficients for wells
C    WLKNO             generates known vector for wells
c    wlupdate          updates the potential at collocation points
c    wlupdate_check    check on wlupdate using rfpot
C    WLSUB             substitutes solution vector into strengths parameters for wells
C
c
c --------------------------------------------------------------------------------------
c
      subroutine wlmatsize (M,N)
c
c --------------------------------------------------------------------------------------
c
c     Routine adds the number of equations and strength parameters
c     to the matrix array dimensions M and N, respectively.
c
      implicit none
      INTEGER i,M,N,N0
      INCLUDE 'wlcom.inc'
c
      nwldischarges=0
      if (nwl.gt.0) then
        N0=N
        do i=1,nwl
          if (lwlh(i)) then
            M=M+1
            N=N+1
          end if
        end do
        nwldischarges=N-N0
      end if
      return
c
      end subroutine
c
c -------------------------------------------------------------------------------------
c
      subroutine wlkeep(drb,j,m)
c
c --------------------------------------------------------------------------------------
c
c
c     add number of well discharges to array counter and set corresponding
c     values of drb array equal to 0. The drb array will contain any strength differences between the
c     solution and the actual stored strength parameters, see lskeepsigma and lkkeepstrength routines.
c
      INTEGER i,j,m
      REAL(8) drb
      DIMENSION drb(m)
      include 'wlcom.inc'
      if (nwldischarges.gt.0) then  ! we have wells in the matrix
       do i=1,nwldischarges
         j=j+1
         drb(j)=0.0d0
       end do
      endif
      return
      end subroutine
c
C --------------------------------------------------------------------------------------
C
      SUBROUTINE WLCZC (CZI,M,N,DRFAC,CALPH,ITYPE)
C
C --------------------------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER(4) M,N,ITYPE,IW,ICPT
      REAL(8) DRFAC,RTOL
      COMPLEX(8) CZI,CALPH
      INCLUDE 'WLCOM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      DIMENSION CZI(*),DRFAC(4,*),CALPH(*),ITYPE(*)
      IF (NWL.EQ.0) RETURN
      DO 10 IW=1,NWL
      IF (LWLH(IW)) THEN
            M=M+1 ! add one equation for each head specified well
            N=N+1 ! add one strength parameter (unknown well pumping rate) for each head specified well
            CZI(N)=CWLZ(IW)+CMPLX(0.0d0,RWLR(IW)) ! make sure is same as in WLERROR
            ITYPE(N)=1
            DRFAC(1,N)=1.0D00
            DRFAC(4,N)=1.0D0
            CALPH(N)=(0.0D0,0.0D0)
C --------------------------------  Check for nearby control points
      RTOL=10.0*RWLR(IW)   ! 10 times well radius
      DO 5 ICPT=1,M-1
      IF (ABS(CZI(M)-CZI(ICPT)).LT.RTOL)
     & WRITE (ILUER,1000) IW,AWLAB(IW),ICPT,CZI(ICPT)
   5  CONTINUE
      ENDIF
  10  CONTINUE
      RETURN
 1000 FORMAT (' ***WARNING: well ',I3,' with label ',A16,/
     &' may be too close to control point # ',I3,' =',2G11.4)      
      END
C
C --------------------------------------------------------------------------------------
C
      SUBROUTINE WLGENMAT (DRA,CZI,M,N,J,DRFAC,CALPH,ITYPE)
C
C --------------------------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER(4) N,M,J,ITYPE,IW,IEQS,IEQ,I
      LOGICAL LNEG
      REAL(8) DRA,DRFAC,RFNFWLCO
      COMPLEX(8) CZ,CZA,CZI,CALPH,CFWLOMC,CFWLQC
      INCLUDE 'WLCOM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      DIMENSION DRA(M,*),CZI(*),CALPH(*),ITYPE(*),DRFAC(4,*)
      IF (NWL.EQ.0) RETURN
      DO 20 IW=1,NWL
      IF (LWLH(IW)) THEN
            J=J+1
            DO 10 I=1,M
            CZ=CZI(I)
            CZA=CALPH(I)
            IEQS=ITYPE(I)
C
C ITYPE=+1 potential specified at CZ
C       -1 difference in potential specified: PHI(CZ)-PHI(CZA)
C       +2 stream function specified at CZ
C       -2 flow across line between CZ & CZA, positive from left to right when at CZ
C       +3 discharge component normal to the unit vector CZA (rotated to the left)
C       +4 discharge component parallel to the unit vector CZA
C        5 continuity equation: provide total discharge
C        6 request zero matrix coefficient
C
            LNEG=IEQS.LT.0
            IEQ=ABS(IEQS)
            GOTO (1,2,3,4,5,1),IEQ
  1         DRA(I,J)=REAL(CFWLOMC(CZ,IW)) ! provide potential at CZ
            GOTO 9
  2         IF (LNEG) THEN
            DRA(I,J)=RFNFWLCO(IW,CZ,CZA) ! provide flow across CZ & CZA
            ELSE
            DRA(I,J)=AIMAG(CFWLOMC(CZ,IW)) ! provide PSI at CZ
            END IF
            GOTO 9
  3         DRA(I,J)=AIMAG(CONJG(CZA)*CFWLQC(CZ,IW))  ! provide Q normal to unit vector CZA
            GOTO 9
  4         DRA(I,J)=REAL(CONJG(CZA)*CFWLQC(CZ,IW))   ! provide Q parallel to unit vector CZA
            GOTO 9
  5         DRA(I,J)=1.0    ! provide total discharge for continuity equation
            GOTO 9
  9         CONTINUE
  10        CONTINUE
      ENDIF
  20  CONTINUE
      RETURN
 1000 FORMAT ('+Generating',I4,' equations, doing equation #: ',I4)      
      END
C
C --------------------------------------------------------------------------------------
C
      SUBROUTINE WLMATCORRECT (DRA,CZI,M,N,J,DRFAC,CALPH,ITYPE)
C
C --------------------------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER(4) M,N,J,ITYPE,IW,IEQS,IEQ,I
      LOGICAL LNEG
      REAL(8) DRA,DRFAC,RFNFWLCO
      COMPLEX(8) CZ,CZA,CZI,CALPH,CFWLOMC,CFWLQC
      INCLUDE 'WLCOM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      DIMENSION DRA(M,*),CZI(*),CALPH(*),ITYPE(*),DRFAC(4,*)
      IF (NWL.EQ.0) RETURN
      DO 20 IW=1,NWL
      IF (LWLH(IW)) THEN
            J=J+1
!      WRITE (ILUME,1000) N,J
            DO 10 I=1,M
            CZ=CZI(I)
            CZA=CALPH(I)
            IEQS=ITYPE(I)
C
C ITYPE=+1 potential specified at CZ
C       -1 difference in potential specified: PHI(CZ)-PHI(CZA)
C       +2 stream function specified at CZ
C       -2 flow across line between CZ & CZA, positive from left to right when at CZ
C       +3 discharge component normal to the unit vector CZA (rotated to the left)
C       +4 discharge component parallel to the unit vector CZA
C        5 continuity equation: provide total discharge
C        6 request zero matrix coefficient
C
            LNEG=IEQS.LT.0
            IEQ=ABS(IEQS)
            GOTO (1,9,9,9,9,6),IEQ
  1         IF (LNEG) DRA(I,J)=DRA(I,J)-REAL(CFWLOMC(CZA,IW)) ! subtract potential at CZA !! to be handled in correction routine
            GOTO 9
  6         DRA(I,J)=0.0    ! provide total discharge for continuity equation
            GOTO 9
  9         CONTINUE
  10        CONTINUE
      ENDIF
  20  CONTINUE
      RETURN
 1000 FORMAT ('+Generating',I4,' equations, doing equation #: ',I4)      
      END

C
C --------------------------------------------------------------------------------------
C
      SUBROUTINE WLKNO (DRB,J,CZI)
C
C --------------------------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER(4) J,IW
      REAL(8) DRB,RPOT,RFPOT,RPOTWL,RFPOTH,rdum1,rdum2
      COMPLEX(8) CZI
      INCLUDE 'WLCOM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      DIMENSION DRB(*),CZI(*)
      IF (NWL.EQ.0) RETURN
      DO 10 IW=1,NWL
      IF (LWLH(IW)) THEN
            J=J+1
c            rdum1=rfpot(czi(j))
            rdum2=dwlpot(iw)
c            write (ilume,1001) iw,rdum1,rdum2,rwlh(iw)
c 1001 format ('wlkno1: iw,rfpot,dwlpot,rwlh ',
c     &        i3,2x,d14.7,2x,d14.7,2x,d14.7)
            rpotwl=rfpoth(rwlh(iw),czi(j))
            DRB(J)=DRB(J)+rpotwl-rdum2
      ENDIF
  10  CONTINUE
      RETURN
       END
C
C --------------------------------------------------------------------------------------
C
      SUBROUTINE WLUPDATE (DRSCR,J,M,N,czi,lincludesub,isubsign)
C
C --------------------------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER(4) J,IW,M,N,nsolOut,isubsign
      LOGICAL lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut,
     &        lincludesub,laddsubcells
      REAL(8) DRSCR,RPOT,RFPOT,rsubsign
      COMPLEX(8) cz,CZI,cflk_subomega
      CHARACTER(8) aBasenameOut
      CHARACTER(16)aDateTimeOut
      INCLUDE 'WLCOM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      DIMENSION DRSCR(*),CZI(*)
      IF (NWL.EQ.0) RETURN
      laddsubcells=lincludesub
      rsubsign=REAL(isubsign)
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      DO 10 IW=1,NWL
      IF (LWLH(IW)) THEN
            J=J+1
            cz=czi(j)
            IF (nsolOut.eq.0) THEN
               dwlpot(iw)=rfpot(cz)
            else
               dwlpot(iw)=dwlpot(iw)+drscr(j)
               if (laddsubcells) then
                 dwlpot(iw)=dwlpot(iw)+rsubsign*REAL(cflk_subomega(cz)) ! sub-cells not in drscr(j)
               end if
            endif
      ENDIF
  10  CONTINUE
      RETURN
       END
C
C --------------------------------------------------------------------------------------
C
      SUBROUTINE WLUPDATE_CHECK (J,M,N,czi)
C
C --------------------------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER(4) J,IW,N,M
      REAL(8) DRSCR,RFPOT,rdum1,rdum2
      COMPLEX(8) CZI
      INCLUDE 'WLCOM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      DIMENSION CZI(m)
      IF (NWL.EQ.0) RETURN
      DO 10 IW=1,NWL
      IF (LWLH(IW)) THEN
            J=J+1
            RDUM1=RFPOT(CZI(J))
            RDUM2=DWLPOT(IW)
            write (ilume,1000) j,rdum1,rdum2
      ENDIF
  10  CONTINUE
      RETURN
 1000 format ('wlupdate_check: j,rfpot,dwlpot ',
     &          i4,2x,2(d14.6,1x))
       END
C
C --------------------------------------------------------------------------------------
C
      SUBROUTINE WLSUB (DRB,J)
C
C --------------------------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER(4) J,IW
      REAL(8) DRB
      INCLUDE 'WLCOM.INC'
      INCLUDE 'TRACOM.INC'
      DIMENSION DRB(*)
      IF (NWL.EQ.0) RETURN
      DO 10 IW=1,NWL
      IF (LWLH(IW)) THEN
            J=J+1
            RWLQ(IW)=RWLQ(IW)+DRB(J)
      ENDIF
  10  CONTINUE
      RETURN
      END
