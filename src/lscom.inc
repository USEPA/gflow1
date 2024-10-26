      COMPLEX(8) CLSZS,CLSZE,CLSCONST
      LOGICAL LSGIV,LSDRAIN,LSHIRE,LSHIPE,LSEND,LSBASE,LSFIRST,
     &        LSGALLERY,LSINLET,LSOUTLET,LSLAKE,LSCONNECT,
     &        LSINLETSTRING,LSOUTLETSTRING,LSLAKESTRING
      REAL(8) RLSPOT,RLSIG,RLSH,RLSRES,RLSWID,RLSW0,RLSC0,RLSD0,RLSDEP,
     &RLSBF,RLSOF,RLSOFSIG,RLSOFST,RO2PI,RO8PI,RLSBFMX,RLSOFMX,RLSOVHW,
     &RLSERROR,RLSQ,RLSHMIN,RLSTAGE,RLSTABLE,
     &RLSH0,RLSH1,RLSH2,RLSQ1,RLSQ2,RLSEVAP
      INTEGER(4) NLS,NLSIG,NLSH,NLSTRING,ILSPLOT,NLSMAX,NLSTMX,
     &NLSTBLMAX,NLSTBLSIZEMAX,NLSTAB,NLSTABLENGTH,nlakeiterations,
     &ilsbound,ilsbound0
      CHARACTER(16) ALSLAB,ALSTBLFILENAME,ATEMP
c
      PARAMETER (NLSMAX=5000,NLSTMX=3000) ! max. # line-sinks, max. # strings
      PARAMETER (NLSTBLMAX=30,NLSTBLSIZEMAX=70) ! max. # tables, max. table size
c                                                 (river and lake stage tables)
c
      INTEGER(4) KLSPTS,KLSPTH,KLSTRING,KLSTREND,KLSUP,KLSDN
      COMMON /LSCMN/ NLS,CLSZS(NLSMAX),CLSZE(NLSMAX),RLSPOT(NLSMAX),
     .               RLSIG(NLSMAX),RLSH(NLSMAX),RLSRES(NLSMAX),
     .               KLSPTS(NLSMAX),KLSPTH(NLSMAX),NLSIG,NLSH,
     .               RLSWID(NLSMAX),ALSLAB(NLSMAX),RLSERROR(NLSMAX),
     .               RLSW0,RLSC0,RLSD0,RLSDEP(NLSMAX),LSGIV(NLSMAX),
     .               LSDRAIN(NLSMAX),LSGALLERY(NLSMAX),LSHIRE,LSHIPE,
     .               NLSTRING,KLSTRING(NLSTMX),KLSUP(NLSMAX),
     .               KLSDN(NLSMAX),RLSBF(NLSMAX),RLSOF(NLSMAX),
     .               RLSOFSIG(NLSMAX),KLSTREND(NLSTMX),RLSOFST(NLSTMX),
     .               RO2PI,RO8PI,RLSBFMX,RLSOFMX,LSEND,LSBASE,LSFIRST,
     .               ILSPLOT,RLSOVHW(NLSTMX),CLSCONST(NLSMAX),
     .               RLSQ(NLSMAX),RLSHMIN(NLSMAX),LSINLETSTRING,
     .               LSOUTLETSTRING,LSLAKESTRING,
     .               LSINLET(NLSMAX),LSOUTLET(NLSMAX),LSCONNECT(NLSMAX),
     .               LSLAKE(NLSMAX),RLSTAGE(NLSTBLSIZEMAX,NLSTBLMAX),
     .               RLSTABLE(NLSTBLSIZEMAX,NLSTBLMAX),
     .               RLSH1(NLSTMX),RLSH2(NLSTMX),
     .               RLSQ1(NLSTMX),RLSQ2(NLSTMX),
     .               ALSTBLFILENAME(NLSTMX),NLSTAB(NLSTMX),
     .               RLSH0(NLSTMX),RLSEVAP(NLSTMX),
     .               NLSTABLENGTH(NLSTMX),ATEMP,nlakeiterations,
     .               ilsbound(nlsmax),ilsbound0
C
C     NLSMAX      maximum number of line sinks
C     NLSTMAX     maximum number of streams sections (stream commands)
C     NLS         current number of line sinks
C     CLSZS       starting point of line sink (downstream end)
C     CLSZE       end point of line sink (upstream end)
C     RLSPOT      calculated potential at center of line sink
C     RLSIG       strength of line sink (exfiltration rate) [L^2/T]
C     RLSH        specified head at center of line sink
C     RLSRES      resistance between stream and aquifer
C     KLSPTS      pointer array to discharge specified line sinks
C     KLSPTH      pointer array to head specified line sinks
C     NLSIG       number of strength specified line sinks
C     NLSH        number of head specified line sinks
C     RLSWID      width of stream bed
C     ALSLAB      label of line sink
C     RLSERROR    % error in boundary condition
C     RLSW0       default width of stream bed
C     RLSC0       default resistance of stream bed
C     RLSD0       defaults depth of resistance layer below the stream
C     RLSDEP      depth of resistance layer below the stream
C     LSGIV       .TRUE. if taken out of the matrix solution procedure (drain, etc.)
C     LSDRAIN     .TRUE. if line sink is a drain feature.
C     LSGALLERY   .TRUE. if linesink is part of a gallery.
C     LSHIRE      .TRUE. recharging line sinks are highlighted in the layout
C     LSHIPE      .TRUE. percolating line sinks are highlighted in the layout
C     NLSTRING    number of stream sections (strings)
C     KLSTRING    address of last upstream line sink in string
C     KLSUP       address of next upstream line sink
C     KLSDN       address of next downstream line sink
C     RLSBF       baseflow at downstream end of line sink [L^3/T]
C     RLSOF       overland flow at downstream end of line sink [L^3/T] or
C     RLSOFSIG    infiltration rate for line sink due to overland flow [L^2/T]
C                 (not used directly in GW calculations)
C     KLSTREND    end address of line sink string
C     RLSOFST     total overland flow for line sink string [L^3/T]
C     RLSBFMX     maximum baseflow in the network
C     RLSOFMX     maximum overland flow in the network
C     LSEND       .TRUE. if "end" command encountered when reading a string
C     LSBASE      .TRUE. if baseflow solution generated
C     LSFIRST     .TRUE. when first time in LSCZC (control point routine)
C     ILSPLOT     plot code: 0 - layout ; 1 - baseflow ; 2 - overlandflow ;
C                            3 - streamflow
C     RLSOVHW     overland inflow at head water of string (l**3/t) or
C     CLSCONST    line-sink specific constant used in potential function
C     RLSQ        total discharge for a gallery
C     RLSHMIN     minimum head allowed in gallery
C     LSINLETSTRING .TRUE. if endstream that is inlet to a lake is encountered
C     LSOUTLETSTRING .TRUE. if stream that is outlet to for a lake is encountered
C     LSINLET     line-sink that serves as inlet for a lake (only for "high-k" lakes)
C     LSOUTLET    line-sink that serves as outlet for a lake
C     LSCONNECT   true if outlet stream connects to a line-sink lake
C     LSLAKE      line-sink that is part of a lake
C     RLSTAGE(I,J) stage level i of j-th outlet stream
C     RLSTABLE(I,J) outflow rate i of j-th outlet stream or
C                   area i of the j-th lake
C     RLSH1,RLSH2   successive estimated heads for a lake
c     RLSQ1,RLSQ2   successive water balance deficiencies for a lake
C     ALSTBLFILENAME filename containing table for j-th outlet stream or lake
C     RLSH0        lake bottom elevation
C     RLSEVAP      evapotranspiration rate (l/t) when string is a lake
C     NLSTAB       number of outlet stream or lake tables
C     NLSTABLENGTH number of elements in table
C     NLAKEITERATIONS current number of lake iterations
C     ATEMP        temporary filename variable for table files
c     ilsbound     code: 0 no action, 1 line-sink at center of stream, 2 line-sink at surface water boundary
c     ilsbound0    temporary storage for width parameter code: in LSIN
