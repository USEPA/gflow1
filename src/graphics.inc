c     GRAPHICS.INC
      real(4)      VVRBIG,VVRMARKHT,VVDOWN,VVUP,VVCUML,VVCURX,VVCURY,
     &             VVDASHL,VVDASHSP,VVDOTL,VVDOTSP
      integer(4)   VVIGRERR,VVITYPE,VVIPEN,VVMARKHT,VVMARKTP,VVNBUTN,
     &             VVILPEN,VVIMPEN,VVIMTYP,VVILINWID,VVILINTYP,
     &             VVIDDPEN,VVCOLTBL(16),VVSYMTBL(6),VVMARKFAC
      logical(1)   VVLDEVOPN,VVLCURADDR,VVLHCOPY
      character(1) VVAMSCHR
      common/GR/   VVRBIG,VVRMARKHT,VVDOWN,VVUP,VVCUML,VVCURX,VVCURY,
     &             VVDASHL,VVDASHSP,VVDOTL,VVDOTSP,
     &             VVIGRERR,VVITYPE,VVIPEN,VVMARKHT,VVMARKTP,VVNBUTN,
     &             VVILPEN,VVIMPEN,VVIMTYP,VVILINWID,VVILINTYP,
     &             VVIDDPEN,
     &             VVCOLTBL,VVSYMTBL,VVMARKFAC,
     &             VVLDEVOPN,VVLCURADDR,VVLHCOPY,
     &             VVAMSCHR
c
c     Error Codes
      integer(4) VVMYGRAPH        ! General flag returned for device handles
      parameter (VVMYGRAPH=-9999)
      integer(4) VVFAIL           ! Failure return code
      parameter (VVFAIL=-1)
      integer(4) VVOK             ! OK return code
      parameter (VVOK=0)
      integer(4) VVNOTIMP         ! Not yet implemented
      parameter (VVNOTIMP=1)
      integer(4) VVNOTOPN         ! Device not open
      parameter (VVNOTOPN=2)            
      integer(4) VVNOMOUSE        ! No Mouse available
      parameter (VVNOMOUSE=3) 
      integer(4) VVNOTCUR         ! Not in cursor addressing mode
      parameter (VVNOTCUR=4) 
      integer(4) VVNODEV          ! No such device
      parameter (VVNODEV=5) 
      integer(4) VVNODVR          ! No such driver for Graphoria hardcopy
      parameter (VVNODVR=6) 
c
      integer(4) VVPENUP          ! For plot()
      parameter (VVPENUP=3)
      integer(4) VVPENDOWN
      parameter (VVPENDOWN=2)
