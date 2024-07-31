      INTEGER(4) ILUIN,ILUOUT,ILUME,ILUER,ILUTMP,ILUPL,ILUECH,ILUDUM
      LOGICAL LECHO,LUCON,LUOUTFILE
      CHARACTER*60 AMESS
      COMMON /LUSYS/ ILUIN,ILUOUT,ILUME,ILUER,ILUTMP,
     &               ILUPL,ILUECH,ILUDUM(2)
      COMMON /LUEXT/ LECHO,LUCON,LUOUTFILE,AMESS(20)
C
C     ILUIN       input of commands from keyboard or file   (*.dat)
C     ILUOUT      output for debugging to screen or file    (*.out)
C     ILUME       output of prompts and messages to screen or file (message.log)
C     ILUER       output of error messages to screen or file (error.log)
C     ILUTMP      for temporary logical unit assignment
C     ILUPL       for plotting instructions to file or plotting device
C     ILUECH      echo of input commands and parameters to file (echo.log)
C     ILUDUM      spare logical units
C
C     LECHO      true if input statements are echoed to a file
C     LUCON      true if input from keyboard (console)
C     LUOUTFILE  true if ILUOUT assigned to a file
C     AMESS      messages string for HALT routine
C


