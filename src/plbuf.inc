      INTEGER(4) IRBUFSIZE,IKBUFSIZE,NPLBYTES
      INTEGER(2) KBUF
      REAL RBUF
      PARAMETER (IKBUFSIZE=2000,IRBUFSIZE=10)
      COMMON /PLBUF/ RBUF(IRBUFSIZE),KBUF(IKBUFSIZE),NPLBYTES
C
C     IRBUFSIZE   size of real buffer RBUF
C     IKBUFSIZE   size if integer*2 buffer KBUF
C     NPLBYTES    number of bytes stored on disk
C     KBUF        integer*2 plot buffer
C     RBUF        real plot buffer
C
