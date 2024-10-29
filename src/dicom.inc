      REAL(8)RDIZ,RDIR,RDIS,RDIH,RDICPT,RDIA,RDIERR,RTAU,RTAU2,RDISR,
     &     RDISR2,RTMR,RTPR,RTDN,ROOT,RTARA,RK,RK2,REPS,RPI,RPI2,RPI4
      INTEGER(4) NDIS,NDIMAG,NDIRH,NDIMAX,IEL,N_DI_INIT
      LOGICAL LDIH
      CHARACTER(16) ADILAB
C
      PARAMETER (NDIMAX=30)             ! maximum number of 3D sink discs
c
      PARAMETER (N_DI_INIT=8*NDIMAX)
      COMMON /DISC/ RDIZ(3,NDIMAX),RDIR(NDIMAX),RDIS(NDIMAX),
     &              RDIH(NDIMAX),RDICPT(NDIMAX),RDIA(NDIMAX),
     &              RDIERR(NDIMAX),
     &              ADILAB(NDIMAX),LDIH(NDIMAX),NDIS,NDIMAG,NDIRH
      COMMON /DIEX/ RTAU,RTAU2,RDISR,RDISR2,RTMR,RTPR,RTDN,ROOT,
     &              RTARA,RK,RK2,REPS,RPI,RPI2,RPI4
      COMMON/IEL/ IEL     
C
C     NDIMAX      maximum number of discs
C     RDIZ(3,I)   coordinates of the center of disc I
C     RDIR(I)     radius of disc I
C     RDIS(I)     exfiltration rate of disc I (L/T)
C     RDIH(I)     specified head at the center of disc I
C     RDICPT(I)   distance from center of disc I where head is specified
C     RDIA(I)     area of disc I
C     RDIERR(I)   relative error in head condition for disc I
C     ADILAB(I)   label of disc I
C     LDIH(I)     .TRUE. if disc I is head specified
C     NDIS        number of discs
C     NDIMAG      number of times the aquifer and its lower or upper image
C                 are imaged upward and downward
C     NDIRH       number of aquifer thicknesses away from the rim of the
C                 disc where a switch is made to a two dimensional solution
C     RTAU        horizontal distance from the center of the current disc
C     RTAU2       RTAU squared
C     RDISR       radius of the current disc
C     RDISR2      RDISR squared
C     RTMR        RTAU-RDISR
C     RTRP        RTAU+RDISR
C     RTDN        RTMR/RTPR
C     ROOT        SQRT(RTAU^2+RETA^2+RDISR^2+RTARA)
C     RTARA       2.0*RTAU*RDISR
C     RK          modulus of the elliptic integrals
C     RK2         RK squared
C     REPS        0.01
C     RPI         3.1415.......
C     RPI2        2*3.1415.....
C     RPI4        4*3.1415.....

