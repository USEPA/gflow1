C     Last change:  HMH   6 Jan 2001    3:21 pm
C     This file contains the following routines and functions
C
C     subroutine W3NEAR   Routine sets the logical L3DEND=.true. when trace within well radius.
c
c     ---------------------------------------------------------------------------------
c
      subroutine W3NEAR (CZ,CZNEW,RZ0,RZNEW,LREDO)
c
c     ---------------------------------------------------------------------------------
c
c     Routine sets the logical L3DEND=.true. when trace within well radius.
c
      implicit none
      INTEGER(4) i,iw,ist,ien
      LOGICAL l3dend,l3drev,lsrce,lredo,lscrn,
     &        lstartunder,lstartover,lendunder,lendover
      REAL(8) rz0,rznew,rsd0,rstep,rq,rdis,rzone,rstmax,rdisi(3),
     &     rtest(3),rscalp,rfsp3d,radw,rlength,rdisend
      COMPLEX(8) cz,cznew,cz0,cdis,cdisend
      include 'w3com.inc'
      include 'tracom.inc'
      include 'lusys.inc'
c
      if (nw3.eq.0) return   ! no partially penetrating wells, skip
      call getstep (rsd0,rstep,l3dend,l3drev)
      if (l3dend) RETURN ! streamline ended elsewhere, skip
      do iw=1,nw3
      ist=ipnt(iw)
      ien=ipnt(iw+1)-1
      radw=rw3rad(iw)
        rq=0.0 ! start calculating the discharge
        do i=ien,ist,-1
        rq=rq+rw3s(i)*rw3l(i)
        END do ! end calculating discharge
        lsrce=rq.lt.0.0
        if (rq.eq.0.0) RETURN ! well does not pump, skip
        if (l3drev.and.lsrce.or..not.l3drev.and..not.lsrce) then !may be approaching the well
          cz0=CMPLX(rw3st(1,ist),rw3st(2,ist))
          cdis=cz-cz0
          rdis=ABS(cdis)
          rzone=radw+rsd0
          if (rdis.lt.rzone) then ! starting point near the well, start checking end point
            rlength=abs(cznew-cz)
            rlength=SQRT(rlength*rlength+(rznew-rz0)*(rznew-rz0))
            if (rlength.gt.rdis+radw) THEN ! step too large, slow down
              rstep=rdis
              call setstep (rstep,l3dend)
              lredo=.true.
              return
            endif
          endif
          if (.not.lredo) then
            cdisend=cznew-cz0
            rdisend=abs(cdisend)
            lscrn=rw3st(3,ist).lt.rznew.and.rw3st(3,ien-1).gt.rznew
            if (rdisend.lt.radw.and.lscrn) THEN ! end point inside the well, stop the trace
              l3dend=.true.
            end if
          endif
          if (.not.lredo.and..not.l3dend) then
            lstartunder=rz0.lt.rw3st(3,ist)
            lstartover=rz0.gt.rw3st(3,ien-1)
            lendunder=rznew.lt.rw3st(3,ist)
            lendover=rznew.gt.rw3st(3,ien-1)
            if (lstartunder.and.lendover.or.
     &      lstartover.and.lendunder) then  !  crossed the well in vertical direction, stop trace
              l3dend=.true.
              if (lstartunder) rznew=rw3st(3,ist)
              if (lstartover)  rznew=rw3st(3,ien-1)
            end if
          endif
          if (.not.lredo.and..not.l3dend) then
            rdisi(1)=REAL(cdis)
            rdisi(2)=AIMAG(cdis)
            rdisi(3)=0.0
            rtest(1)=REAL(cdisend)
            rtest(2)=AIMAG(cdisend)
            rtest(3)=0.0
            rscalp=rfsp3d(rdisi,rtest)
            if (rscalp.lt.0.0) THEN ! crossed the well in a horizontal plane, stop the trace
              l3dend=.true.
              if (lendunder) rznew=rw3st(3,ist)
              if (lendover)  rznew=rw3st(3,ien-1)
            end if
          end if
        end if
        if (l3dend) then
          iFlag=2
          iElementType=4
          if (LEN(AW3LAB(IW)).GT.0)
     &    aElementLabel=AW3LAB(IW)
          rstep=0.1*radw
          call setstep (rstep,l3dend)
          return
        endif
      end do
      return
      end subroutine
