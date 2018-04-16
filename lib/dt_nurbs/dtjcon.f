      INTEGER FUNCTION DTJCON(I)
C ******************************************
C *                                        *
C *  COPYRIGHT  1987  THE BOEING COMPANY   *
C *                                        *
C *  ALL  RIGHTS  RESERVED                 *
C *                                        *
C ******************************************
C
C MACHINE CONSTANTS FOR PORTABLE CODE SUPPORT ON:
C     IRIS-4D
C
C INPUT  - I (IN 1..15)  IS THE INDEX OF THE DESIRED CONSTANT
C OUTPUT - DTJCON IS THE DESIRED CONSTANT
C ERRORS - IF I IS OUT OF RANGE, DTJCON(1) IS RETURNED
 
      INTEGER K(15)
 
      DATA K(1)/2147483647/
      DATA K(2)/2147483647/
      DATA K(3)/2147483647/
      DATA K(4)/6/
      DATA K(5)/5/
      DATA K(6)/6/
      DATA K(7)/80/
      DATA K(8)/133/
      DATA K(9)/9/
      DATA K(10)/7/
      DATA K(11)/16/
      DATA K(12)/512/
      DATA K(13)/4/
      DATA K(14)/1/
      DATA K(15)/2/
 
      J = I
      IF ( J .LE. 0  .OR.  J .GT. 15 ) J = 1
      DTJCON = K(J)
 
      END
