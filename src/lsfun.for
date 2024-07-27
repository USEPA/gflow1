c       This file contains the following routines and functions:
c
C	SUBROUTINE LSQI         calcuates the discharge vector for all line-sinks
C	COMPLEX FUNCTION COMLS  calculates the complex potential for line-sink between CZ1 and CZ2 
c                             with sink density 1
C	COMPLEX FUNCTION CFLSG  complex potential for all line-sinks with given sink density
C	COMPLEX FUNCTION CFLSU  complex potential for all head specified line-sinks
C	COMPLEX FUNCTION CDLS   derivative of COMLS
C	COMPLEX FUNCTION CDLSG  derivative of CFLSG
C	COMPLEX FUNCTION CDLSU  derivative of CFLSU
C	REAL FUNCTION RFLSPT    real discharge potential due to all line-sinks
c
C
C-----------------------------------------------------------------------------------------------------
C
      SUBROUTINE LSQI (CZ,RQI)
C
C-----------------------------------------------------------------------------------------------------
C
C     Discharge function for constant strength line sinks
C
      IMPLICIT NONE
      INTEGER(4) I
      REAL(8) RQI,RS
      COMPLEX(8) CZ,CZ1,CZ2,CLSDIS,CDLS
      INCLUDE 'LSCOM.INC'
      INCLUDE 'TRACOM.INC'      
      INCLUDE 'LUSYS.INC'
      DIMENSION RQI(3)
      IF (NLS.EQ.0) RETURN
      CLSDIS=(0.0,0.0)
      DO 10 I=1,NLS
      RS=RLSIG(I)      
      IF (RS.EQ.0.0) GOTO 10
      CZ1=CLSZS(I)
      CZ2=CLSZE(I)
      CLSDIS=CLSDIS-RS*CDLS(CZ,CZ1,CZ2)
  10  CONTINUE
      RQI(1)=RQI(1)+REAL(CLSDIS)
      RQI(2)=RQI(2)-AIMAG(CLSDIS)
      RETURN
      END
C
C-----------------------------------------------------------------------------------------------------
C
	COMPLEX(8) FUNCTION COMLS(CZ,CZS,CZE)
C
C-----------------------------------------------------------------------------------------------------
C
C     Last change:  HMH  18 Dec 2006    2:55 pm
c     Complex potential at CZ for line sink with endpoints CZS and CZE and sink density 1.
C     COMLS lacks a constant calculated by LSPREP and added to COMLS where called.
c
      IMPLICIT NONE
      COMPLEX(8) COM1,COM2,CBZ,CZ,CZS,CZE,CFBIGZ
      INCLUDE 'LSCOM.INC'
      INCLUDE 'TRACOM.INC'
      include 'lusys.inc'
      CBZ=CFBIGZ(CZ,CZS,CZE)
      COM1=(.0,.0)
      COM2=(.0,.0)
      IF(ABS(CBZ+1.).GT..000001) COM1=(CBZ+1.)*LOG(CBZ+1.)
      IF(ABS(CBZ-1.).GT..000001) COM2=(CBZ-1.)*LOG(CBZ-1.)
      COMLS=COM1-COM2
      COMLS=.5*RO2PI*ABS(CZE-CZS)*COMLS
C     WARNING: a constant CLSCONST(iad) must be added, where iad is the line-sink address
C              with czs and cze as end points. This constant is added in CFLSG and CFLSU and in LSMAT
      RETURN
      END
C
C-----------------------------------------------------------------------------------------------------
C
      COMPLEX(8) FUNCTION CFLSG(CZ)
C
C-----------------------------------------------------------------------------------------------------
C
c     Complex potential for all line sinks with given strength
c
      IMPLICIT NONE
      INTEGER(4) K,I
      COMPLEX(8) CZ,COMLS
      INCLUDE 'LSCOM.INC'
      INCLUDE 'TRACOM.INC'
      CFLSG=(.0,.0)
      IF (NLSIG.EQ.0) RETURN
      DO 100 I=1,NLSIG
      K=KLSPTS(I)
      CFLSG=CFLSG+RLSIG(K)*(COMLS(CZ,CLSZS(K),CLSZE(K))+CLSCONST(K))
 100  CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------------------------------------
C
      COMPLEX(8) FUNCTION CFLSU(CZ)
C
C-----------------------------------------------------------------------------------------------------
C
c     Complex potential for all head specified line sinks
c
      IMPLICIT NONE
      INTEGER(4) I,K
      REAL(8) RS
      COMPLEX(8) CZ,COMLS
      INCLUDE 'LSCOM.INC'
      INCLUDE 'TRACOM.INC'
      CFLSU=(.0,.0)
      IF (NLSH.EQ.0) RETURN
      DO 100 I=1,NLSH
      K=KLSPTH(I)      
      RS=RLSIG(K)
      IF (RS.EQ.0.0) GOTO 100
      CFLSU=CFLSU+RS*(COMLS(CZ,CLSZS(K),CLSZE(K))+CLSCONST(K))
 100  CONTINUE
      RETURN
      END
C----------------------DERIVATIVES--------------------------
C
C-----------------------------------------------------------------------------------------------------
C
      COMPLEX(8) FUNCTION CDLS(CZ,CZS,CZE)
C
C-----------------------------------------------------------------------------------------------------
C
C     Derivative of COMLS (omega for line sink with endpoints czs and cze) with respect to CZ
C
      IMPLICIT NONE
      COMPLEX(8) CZ,CZS,CZE,CBZ,CBZP1,CBZM1,CFBIGZ
      INCLUDE 'LSCOM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      CBZ=CFBIGZ(CZ,CZS,CZE)
      CBZP1=CBZ+1.0
      CBZM1=CBZ-1.0
      IF (ABS(CBZP1).LT.0.000001) CBZP1=(-0.000001,0.0)
      IF (ABS(CBZM1).LT.0.000001) CBZM1=(0.000001,0.0)
      CDLS=LOG(CBZP1)-LOG(CBZM1)
      CDLS=RO2PI*ABS(CZE-CZS)*CDLS/(CZE-CZS)
      RETURN
      END
C
C-----------------------------------------------------------------------------------------------------
C
      COMPLEX(8) FUNCTION CDLSG(CZ)
C
C-----------------------------------------------------------------------------------------------------
C
c     Derivative of omega for all strength specified line sinks
c
      IMPLICIT NONE
      INTEGER(4) I,K
      COMPLEX(8) CZ,CDLS
      INCLUDE 'LSCOM.INC'
      INCLUDE 'TRACOM.INC'
      CDLSG=(.0,.0)
      IF (NLSIG.EQ.0) RETURN
      DO 100 I=1,NLSIG
      K=KLSPTS(I)
      CDLSG=CDLSG+RLSIG(K)*CDLS(CZ,CLSZS(K),CLSZE(K))
 100  CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------------------------------------
C
      COMPLEX(8) FUNCTION CDLSU(CZ)
C
C-----------------------------------------------------------------------------------------------------
C
c     Derivative of omega for all head specified line sinks
c
      IMPLICIT NONE
      INTEGER(4) I,K
      REAL(8) RS
      COMPLEX(8) CZ,CDLS
      INCLUDE 'LSCOM.INC'
      INCLUDE 'TRACOM.INC'
      CDLSU=(.0,.0)
      IF (NLSH.EQ.0) RETURN
      DO 100 I=1,NLSH
      K=KLSPTH(I)
      RS=RLSIG(K)
      IF (RS.EQ.0.0) GOTO 100
      CDLSU=CDLSU+RS*CDLS(CZ,CLSZS(K),CLSZE(K))
 100  CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFLSPT (CZ)
C
C-----------------------------------------------------------------------------------------------------
C
c     Discharge potential due to all line sinks
c
      IMPLICIT NONE
      COMPLEX(8) CZ,CFLSG,CFLSU
      INCLUDE 'LSCOM.INC'
      INCLUDE 'TRACOM.INC'
      RFLSPT=-REAL(CFLSG(CZ)+CFLSU(CZ))
      RETURN
      END
