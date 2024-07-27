      INTEGER(4) IPNT,NW3,NLIMAG,NW3RH,NW3WMX,NW3SMX,nw3strengths
      LOGICAL LW3Q,LW3DRT,LW3BOT,LW3TOP
      REAL(8) RW3RAD,RW3HED,RW3Q,RW3ST,RW3CP,RW3S,RW3L,RW3ALP,
     &     RZONEMLTP,rw3d0,dw3pot
      CHARACTER(16) AW3LAB
c
      PARAMETER (NW3WMX=200)  ! maximum number of partially penetrating wells
c
      PARAMETER (NW3SMX=4*NW3WMX)  ! 4 control points per well
      COMMON /W3COMN/ RW3RAD(NW3WMX),RW3HED(NW3WMX),RW3Q(NW3WMX),
     &                IPNT(NW3WMX),LW3Q(NW3WMX),LW3DRT(NW3WMX),
     &                LW3BOT(NW3WMX),LW3TOP(NW3WMX),AW3LAB(NW3WMX),
     &                RW3ST(3,NW3SMX),RW3CP(3,NW3SMX),RW3S(NW3SMX),
     &                RW3L(NW3SMX),NLIMAG,NW3,NW3RH,RW3ALP,
     &                rw3d0(3,nw3wmx),DW3POT(NW3SMX),nw3strengths
      COMMON /W3ADD1/ RZONEMLTP
C
C     NW3WMX      maximum number of wells
C     NW3SMX      maximum total number of strength parameters for all wells
C     RW3RAD(I)   radius of well number I
C     RW3HED(I)   specified head at well number I
C     RW3Q(I)     total discharge of well number I (L^3/T)
C     IPNT(I)     pointer to starting line sink of well number I
C     LW3Q(I)     .TRUE. if well number I is discharge specified
C     LW3DRT(I)   .TRUE. if well number I is a double root element
C     LW3BOT(I)   .TRUE. if well number I touches the aquifer bottom
C     LW3TOP(I)   .TRUE. if well number I touches the aquifer top
C     AW3LAB(I)   label of well number I
C     RW3ST(3,J)  coordinates of line sinks for well 
C     RW3CP(3,J)  coordinates of control points of well
C     RW3S(J)     strength parameters of line sink
C     RW3L(J)     length of line sink
C     NLIMAG      number of image aquifer (above and below actual aquifer)
C     NW3         total number of wells
C     NW3RH       currently not in use
C     RW3ALP      factor to calculate off set of control point from well end
C     rw3d0(3,i)  vector (3d) from well center to point CZ (used in W3NEAR)
C     DW3POT      potential at control points
C     nw3strengths number of 3D well strength parameters in the matrix
C     RZONEMLTP   factor times aquifer thickness is multiplied for 3D zone
C

