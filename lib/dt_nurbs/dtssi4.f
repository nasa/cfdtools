        subroutine dtssi4(p, t, h, q)
C  **
C  **           Given base point p and direction t in three space
C  **           and step h, compute q = p + h * t.
C  **
C  **   input:
C  **
C  **           p       base point
C  **
C  **           t       direction vector
C  **
C  **           h       step-size
C  **
C  **   output:
C  **
C  **           q       = p + h * t
C  **
        double precision p(*), t(*), q(*)
        double precision h
        integer i
        do 10 i = 1, 3
        q(i) = p(i) + h * t(i)
10      continue
        return
        end
