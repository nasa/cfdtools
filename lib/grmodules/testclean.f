      PROGRAM TESTCLEAN
C     Why won't SCLEAN work?
      WRITE (6,'(A,$)') ' Enter SCLEAN argument: '
      READ  (5, '(I1)') I
      CALL SCLEAN (I)
      STOP ' '
      END
