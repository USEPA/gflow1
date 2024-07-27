C     Last change:  HMH  16 Jan 2001    8:36 pm
c     This file contians the following routines and functions
c
c     REAL FUNCTION RFPDS       returns exfiltration rate due to all (2D) sink discs
c     REAL FUNCTION RFPDSBOTTOM returns bottom exfiltration due to discharge specified discs
c     REAL FUNCTION RFPDH       returns exfiltration rate due to all head specified discs
c
c
c
c
c
c
c
C
C ---------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFPDS(CZ)
C
C ---------------------------------------------------------------------------
C
C
C     Function returns exfiltration rate at CZ due to all discharge specified
C     sinkdiscs. Note, in case of aquifer bottom exfiltration the function
C     returns the net exfiltration rate: top + bottom
C     Called in CURDAT
C
      IMPLICIT NONE
      INTEGER(4) I,IP
      REAL(8) RDIS
      COMPLEX(8) CZ
      INCLUDE 'pdcom.inc'
      INCLUDE 'lusys.inc'
      RFPDS=0.0
      IF (NPDRC.LE.0) RETURN
      DO 10 I=1,NPDRC
      IP=IPDRC(I)
      RDIS=ABS(CZ-CPDZ(IP))
      IF (RDIS.LE.RPDR(IP)) RFPDS=RFPDS+RPDS(IP)
  10  CONTINUE
      RETURN
      END
C
C ---------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFPDSBOTTOM(CZ)
C
C ---------------------------------------------------------------------------
C
C
C     Function returns bottom exfiltration rate at CZ due to all discharge
C     specified  sinkdiscs.
C     Called in EXTRACT
C
      IMPLICIT NONE
      INTEGER(4) I,IP
      REAL(8) RDIS
      COMPLEX(8) CZ
      INCLUDE 'pdcom.inc'
      INCLUDE 'lusys.inc'
      RFPDSBOTTOM=0.0
      IF (NPDRC.LE.0) RETURN
      DO 10 I=1,NPDRC
      IP=IPDRC(I)
      RDIS=ABS(CZ-CPDZ(IP))
      IF (RDIS.LE.RPDR(IP)) RFPDSBOTTOM=RFPDSBOTTOM+RPDSB(IP)
  10  CONTINUE
      RETURN
      END
C
C ---------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFPDH(CZ)
C
C ---------------------------------------------------------------------------
C
C
C     Function returns exfiltration rate at CZ due to all head specified
C     sinkdiscs. Called in CURDAT
C
      IMPLICIT NONE
      INTEGER(4) I,IP
      REAL(8) RDIS
      COMPLEX(8) CZ
      INCLUDE 'pdcom.inc'
      INCLUDE 'lusys.inc'
      RFPDH=0.0
      IF (NPDHD.LE.0) RETURN
      DO 10 I=1,NPDHD
      IP=IPDHD(I)
      RDIS=ABS(CZ-CPDZ(IP))
      IF (RDIS.LE.RPDR(IP)) RFPDH=RFPDH+RPDS(IP)
  10  CONTINUE
      RETURN
      END
