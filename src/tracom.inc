      INTEGER(4) iFlag,iElementType,
     &           iticks_near,iticks_near1,iticks_near2,
     &           iticks_trace,iticks_trace1,iticks_trace2,
     &           inear_call,minimum_iterations
      LOGICAL LCLOSE,LTRACEOUT,LGRAPHICS,LBOUNDARY,lquit,
     $           lakedone
      REAL(8) RDIRI,RCONCTR,rconverge_reference,
     &        rconverge_linesinks,rconverge_linesinks_resistance,
     &        rconverge_inhomogeneities,rconverge_barriers_noflow,
     &        rconverge_barriers_resistance,
     &        rconverge_well_3D,rconverge_well_2D,
     &        rconverge_lake_waterbalance
      complex(8) ctr0,ctr1,ctr2,ctrnew
      CHARACTER(16) aElementLabel
      COMMON /TRACOM/ LCLOSE,RDIRI(3),
     &           iticks_near,iticks_near1,iticks_near2,
     &           iticks_trace,iticks_trace1,iticks_trace2,
     &           inear_call
      COMMON /TRACO2/ LTRACEOUT,RCONCTR,iFlag,
     &                iElementType,aElementLabel
      COMMON /TRACO3/ LGRAPHICS,LBOUNDARY,ctr0,ctr1,ctr2,ctrnew
      common /traco4/ minimum_iterations,rconverge_reference,
     &           rconverge_linesinks,rconverge_linesinks_resistance,
     &           rconverge_inhomogeneities,rconverge_barriers_noflow,
     &           rconverge_barriers_resistance,
     &           rconverge_well_3D,rconverge_well_2D,
     &           rconverge_lake_waterbalance,lquit,lakedone
C
C     LCLOSE        .TRUE. when close to a linear feature (line sink)
C     LTRACEOUT     .TRUE. when writing a streamline file for GUI
C     LBOUNDARY     .TRUE. when streamlines should be stopped at an infinitely long line (Simple WHPA infinitely WhAEM)
C     RDIRI          unit vector along linear feature
C     RCONCTR        solute concentration calculated in TRACE
C     iticks_near    cummulative time used by "..near()" subroutine calls
C     iticks_near1   time before "..near()" call
C     iticks_near2   time after "..near()" call
C     iticks_trace   cummulative time used by pathline tracing
C     iticks_trace1  time before pathline trace
C     iticks_trace2  time after pathline trace
C     inear_call     number of "..near()" calls during a single trace
c     Ctr0             vector from start of streamline to point 1 on
c                    boundary, divided by vector from point 1 to 2
C     Ctr1             point 1 of the boundary
C     Ctr2             vector from point 1 to 2
C     CtrNEW           vector from nwe pathline point to point 1 on
c                    boundary divided by vector from point 1 to 2
C     iFlag          indicates reason for trace termination, see below
c     iElementType   indicates element type where trace stops, see below
c     aElementLabel  label of element where trace ends
C     LGRAPHICS     .TRUE. graphics screen activated.
C     
C     currently included in linesink discharge function only
C     and, of course, in PREDCOR
c
c     iFlag=-1 Error, iFlag not set.
c     iFlag=0  Requested time of travel is exceeded.
C     iFlag=1  Particle left window.
C     iFlag=2  Particle stopped at a boundary (see iElementType).
C     iFlag=3  Particle entered a region where the aquifer is dry.
C     iFlag=4  Particle is trapped in a stagnation point.
C     iFlag=5  Particle crossed aquifer top (entered layer above?).
C     iFlag=6  Particle crossed aquifer bottom (entered layer below?).
c     iFlag=7  Particle stopped by CTR-Break key.
c     iFlag=8  Particle stopped because of infinite velocity (singularity).
c
c     iElementType=-1 No element was reached.
c     iElementType=1  Particle terminated at a well.
c     iElementType=2  Particle terminated at a linesink.
c     iElementType=3  Particle terminated at a sink disc.
c     iElementType=4  Particle terminated at a partially penetrating well.
c     iElementType=5  Particle terminated at a Theis well.
c     iElementType=6  Particle terminated at a 3D sink disc.
c     iElementType=7  Particle terminated at a stream boundary (Simple WHPA in WhAEM)
c
c     aElementLabel="_NO_ELEMENT_" when iFlag /= 2
c     aElementLabel="_NO_LABEL_"   when iFlag=2, but element has no label.
c
c     Below are convergence parameters from the file "converge.tab", which are
c     used in the routines XXERROR and lslakewaterbalance called in SOLUT.
c
c     minimum_iterations
c     rconverge_reference
c     rconverge_linesinks
c     rconverge_linesinks_resistance
c     rconverge_inhomogeneities
c     rconverge_barriers_noflow
c     rconverge_barriers_resistance
c     rconverge_well_3D
c     rconverge_well_2D
c     rconverge_lake_waterbalance
c
c     lquit         is true when all analytic elements have met their convergence criterion
c     lakedone      is true when all lakedone features have met their water balance criteria
c
c
c-------------------------------------------------------------
c
c      Format of POINTS command lines in the trace module:
c
c      POINTS
c       x  y  z  iswitch
c       x  y  z  iswitch
c            etc.      
c
c       iswitch=1  forward trace    iswitch=-1 backward trace
c
c
c
C      