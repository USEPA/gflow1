      INTEGER(4) NDB,NDBSTR,NDBSTI,IDBSTA,NDBZMX,NDBTRM,NDBSMX,
     &           IDOMAINTYPE,istr1old,istr2old,IDBCODE,ndbstrengths
      COMPLEX(8) CDBZ,CDBZEL,CDBZIV,CDBZ21,CDBZ0
      COMPLEX(8) CDDBLN,CDI,CDDPAD,CDDBSS
      REAL(8) RDBK,RDBP,RDBB,RDBT,RDBGAM,RDBRAD,RDBH,
     &     RDBERRSTRT,RDBERRCNTR,rdberrend,RDBW,RDBOFFSETN,
     &     RDBOFFSETC,
     &     RDBTFACNOD,RDBTFACCTR,RDBBFACNOD,RDBBFACCTR,
     &     rdbbtemp,rdbktemp,rdbhnod,rdbhctr,rdbhend,dbknon,dbknoc,
     &     rdbpotaverage,RDBKI,RDBKO,RDBBI,RDBBO,rdboffsetinhom,
     &     rdbtfacend,rdbbfacend,dbknoe
      CHARACTER(16) ADBLAB
      LOGICAL LDBPLT,LDBFAR,LDBNOD,LDBINS,ldbmatrix,ldbrechargeonly
      REAL(8) DBSTR,DBSEN,DBQSTR,DBAVS,DPI,D2PI,DBEPS,DBQSTS,DBAREA
   
c
c Note: NDBZMX and NDBSMX must be even or they will cause misalignment 
c       of double precision variables!!!
C
      PARAMETER (NDBZMX=3000,NDBTRM=17,NDBSMX=1000)
C
      COMMON /DOUBL/ NDB,NDBSTR,NDBSTI(NDBSMX),IDBSTA(NDBSMX),
     &               IDOMAINTYPE(NDBSMX),ldbrechargeonly(ndbsmx),
     &               CDBZ(NDBZMX),DBSTR(NDBZMX),DBSEN(NDBZMX),
     &               DBQSTR(NDBZMX),RDBERRSTRT(NDBZMX),
     &               RDBERRCNTR(NDBZMX),rdberrend(ndbzmx),
     &               CDBZEL(NDBZMX),CDBZIV(NDBTRM,NDBZMX),
     &               CDDBLN(NDBZMX),CDBZ21(NDBZMX),RDBTFACNOD(NDBZMX),
     &               RDBTFACCTR(NDBZMX),RDBBFACNOD(NDBZMX),
     &               RDBBFACCTR(NDBZMX),rdbhnod(ndbzmx),rdbhctr(ndbzmx),
     &               rdbhend(ndbzmx),rdbtfacend(ndbzmx),
     &               rdbbfacend(ndbzmx),
     &               LDBFAR(NDBZMX),LDBNOD(NDBZMX),ADBLAB(NDBZMX),
     &               RDBK(NDBSMX),RDBP(NDBSMX),RDBB(NDBSMX),
     &               RDBT(NDBSMX),DBAVS(NDBSMX),LDBINS(NDBSMX),
     &               DPI,D2PI,CDI,DBEPS,CDDPAD(ndbzmx),LDBPLT,
     &               RDBW(NDBSMX),RDBH(NDBSMX),RDBOFFSETN,RDBOFFSETC,
     &               rdbbtemp(ndbsmx),rdbktemp(ndbsmx),
     &               dbknon(ndbzmx),dbknoc(ndbzmx),dbknoe(ndbzmx),
     &               rdbpotaverage(ndbsmx),
     &               istr1old,istr2old,ldbmatrix,
     &               IDBCODE(NDBZMX),ndbstrengths,
     &               RDBKI(NDBZMX),RDBKO(NDBZMX),RDBBI(NDBZMX),
     &               RDBBO(NDBZMX),rdboffsetinhom
      COMMON /DBRECHC/RDBGAM(NDBSMX),DBAREA(NDBSMX),RDBRAD(NDBSMX),
     &                CDBZ0(NDBSMX),CDDBSS(NDBZMX),DBQSTS(NDBZMX)
C
C     NDBZMX  Maximum number of line doublets.
C     NDBSMX  Maximum number of strings.
C     NDBTRM  Current number of terms in farfield expansion.
C             Note: functions have a maximum of 21 coefficients defined!
C     ---------------------------
C     COMMON /DOUBL/
C     ---------------------------
C     NDB     Number of line doublets in use.
C     NDBSTR  Number of line doublet strings in use.
C     NDBSTI  Number of line doublets in i th string.
C     IDBSTA  Starting address of i th string.
C     IDOMAINTYPE = 1 hydraulic conductivity inhomogeneity
C                   2 open slurry wall
C                   3 closed slurry wall
C     ldbrechargeonly   = true when they domain with idomaintype=1 is a recharge only domain
C     CDBZ    First node of i th line doublet.
C     DBSTR   Strength parameter (double precision) for first node of i th line doublet.
C     DBSEN   Strength parameter (double precision) for last node of i th line doublet.
C             (for jump in T only!)
C     DBQSTR  Strength parameter (double precision) of the quadratic term for i th line doublet. (for jump in T only!)
C     RDBERRSTRT Percent error in the heads inside and outside the domain at first node (not used for slurry walls)
C     RDBERRCNTR Percent error in the heads inside and outside the domain at center of line doublet (not used for slurry walls)
C     rdberrend  percent error in the heads inside and outside the domain at the end of a line doublet (not used for slurry walls)
C     CBDZEL  Big Z for the i th line doublet.
C     CDBZIV  (j,i) 1/Z to the j th power for the i th line doublet.
C     CDDBLN  Ln(Z-1/Z+1) for i th line doublet. Note: CDDBLN=(0.0,0.0)  when
C             either Z-1=(0.0,0.0) or Z+1=(0.0,0.0)!
C     CDBZ21  z2-z1 for i th element.     
C     RDBTFACNOD  when idomaintype=1:  Ti/(Ti-To) at first collocation point, where To is transmissivity outside the
C                 i th string and Ti the transmissivity inside the string (not used for slurry walls).
c                 when idomaintype=2 or 3: average conductance of slurry wall for interval across node.
C     RDBTFACCTR  when idomaintype=1:  Ti/(Ti-To) at center, where To is transmissivity outside the
C                 i th string and Ti the transmissivity inside the string (not used for slurry walls).
c                 when idomaintype=2 or 3: average conductance of slurry wall for interval on line doublet.
C     RDBTFACEND  when idomaintype=1:  Ti/(Ti-To) at last collocation point, where To is transmissivity outside the
C                 i th string and Ti the transmissivity inside the string (not used for slurry walls).
c                 when idomaintype=2 or 3: average conductance of slurry wall for interval across node.
c     RDBBFACNOD  Ti*To*0.5*(b0-bi)/(Ti-To) at first collocation point used in DBKNO.
c     RDBBFACCTR  Ti*To*0.5*(b0-bi)/(Ti-To) at center used in DBKNO.
c     RDBBFACend  Ti*To*0.5*(b0-bi)/(Ti-To) at last collocation point used in DBKNO.
C     RDBHNOD  head at first collocation point used for base jump inhomogeneities (set in DBERROR)
C     RDBHCTR  head at center of line-doublet for base jump inhomogeneities (set in DBERROR)
C     RDBHEND  head at last collocation point of line-doublet for base jump in inhomogeneities (set in DBERROR)
C     LDBFAR  .TRUE. if CABS(CDBZEL(I)).GT.3.0 (farfield expansion invoked
C             for the i th line doublet).     
C     LDBNOD  .TRUE. if CZ coincides with the i th node.
C     ADBLAB  Label to identify individual line doublets.
C     RDBK    Hydraulic conductivity inside the i th string or inside i th slurry wall
C     RDBP    Porosity inside the i th string or inside the i th slurry wall.
C     RDBB    Elevation of the aquifer base inside the i th string or elevation of the bottom of the i th slurry wall
C     RDBT    Elevation of the aquifer top, presently not used.
C     RDBW    Width of i th slurry wall (only used for slurry walls).
C     RDBH    Estimated average head for the i th inhomogeneity domain, only used if base elevation jumps.
C     RDBOFFSETN distance from starting point of line doublet to first collocation point
C             that defines the integration interval for a node.
C     RDBOFFSETC distance from starting point of line doublet to first collocation point
C             that defines the integration interval for line doublet center.
C     DBAVS   Average strength (double precision) of the i th string (not used for slurry walls).
C             The average strength is subtracted from the nodal strength
C             parameters and directly added to the potential inside the string.
C     LBDINS  .TRUE. if CZ is inside string ISTR
C     DPI     3.1415..... (double precision).
C     D2PI    6.283...... (double precision).
C     CDI     (0.0,0.1)   (complex*16).
C     DBEPS   Tolerance for avoiding singularities at line doublet endpoints.
C     CDDPAD  ln((z-z3)/(z-z1)) added to omega if CZ at i th node (z2).
C     LDBPLT  =.TRUE. if inhomogeneities are to be included in the layout.
c     RDBBTEMP temporary storage of inside aquifer bottom setting
C     RDBKTEMP temporary storage of inside conductivity setting
C     DBKNON  discharge potential or flow across integration segment at first node of line-doublet
C     DBKNOC  discharge potential or flow across integration segment at center of line-doublet
c     dbknoe  discharge potential at last collocation point of line-doublet
C     RDBPOTAVERAGE average potential for string ISTR for use in DBERROR
c     ISTR1OLD initial value for ldominsidedom
C     ISTR2OLD initial value for ldominsidedom
C     LDBMATRIX is true when PHi and Q evaluations at nodes are to be allowed (in the matrix routines).
c     IDBCODE(INOD)   -1  error, node has not been visited
c                     0  node is shared by another earlier visited node and on a common boundary
c                     1  node is on common boundary between 2 or more domains (coinciding nodes on other
c                                                                         domains are flagged 0)
c                     2  node is shared by one or two domains, but not on a common boundary.
c
c     ndbstrengths    number of doublet strength parameters in the matrix
C     RDBKI   inside K-value opposite center of line doublet
C     RDBKO   outside K-value opposite center of line doublet
C     RDBBI   inside b-value opposite center of line doublet
C     RDBBO   outside b-value opposite center of line doublet
c     rdboffsetinhom  offset of collocation point from node on inhomogeneity domain
c
C     -------------------------
C     COMMON /DBRECH/
C     -------------------------
C     RDBGAM  added areal exfiltration rate (>0) or 
C             added areal recharge rate (<0).
C     DBAREA  area of i th inhomogeneity.
C     CDBZ0   center of inhomogeneity (average of all vertices).
C     RDBRAD  average distance from center to vertices.
C     CDDBSS  complex (double precision) strength at vertices of polygon due to
C             RDBGAM only.
C     DBQSTS  real (double precision) strength at midpoints of polygon sides
C             due to RDBGAM only.
C
C             NOTE: CDDBSS and DBSTR are added together before multipication
C             with doublet/dipole functions. DBQSTS and DBQSTR are added together
C             before multiplication with doublet/dipole functions.

