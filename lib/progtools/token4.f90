!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine token4 (string, seps, list_delimiter, number, tokens)

!  Description and usage:
!
!        TOKEN4 is an extension of TOKEN2 capable of handling (as one token)
!     a string of the type for which SCAN4 was derived from SCAN2.  Such a
!     string may contain blanks or even other standard delimiters, as used
!     for character data representing, for instance, variable name and units
!     in a plottable file, or for lists of integers or reals as handled by
!     list-directed reads, which can be demarked with (say) a pair of
!     parentheses.  As in SCAN4, such a generalized token is defined either
!     by a pair of the input argument character list_delimiter, or by one
!     of the pairs (), [], or {} where only the leading character needs to
!     be input as list_delimiter.
!
!        Unlike the lower-level SCAN4, this utility retains (for such a
!     generalized token) the paired delimiters as the first and last
!     meaningful characters of the output token string.  This allows the
!     application to identify such a token, which may represent a string
!     with embedded blanks (say) or a list of either integers or reals that
!     is application-specific (with a clue probably gleaned from the token
!     that precedes this generalized token in the output tokens array.
!     That array may or may not be an alternating list of keyword-name/
!     keyword-value pairs.  TOKEN4 effectively allows "keyword value"
!     to be an array of values, except that such an array is returned in
!     character form for decoding as array elements at the higher level.
!
!        Processing continues until the requested number of tokens have
!     been found or the end of the input string is reached.  Note that
!     simple keywords and values are returned as upper case, but the case
!     of each character in a generalized token is untouched.
!
!        CORRECTION:  Unlike TOKEN2, TOKEN4 does NOT convert any token to upper
!     case.  The application can do that as the context indicates.  This eases
!     the handling of character-type keyword values.
!
!  Example of Usage Type 1:
!
!     Inputs:
!
!        string:           variables = "X, m" "Y, m" "Z, m" "qw, W/cm^2"
!        seps:             '=, '  (i.e., equal sign, comma, and blank)
!        list_delimiter:   "
!        number:           at least 5
!
!     Outputs:
!
!        number:           5
!        tokens(1):        variables
!        tokens(2):        "X, m"
!        tokens(3):        "Y, m"
!        tokens(4):        "Z, m"
!        tokens(5):        "qw, W/cm^2"
!
!  Example of Usage Type 2:
!
!     Inputs:
!
!        string:           @block 1    ifree 2 ibc=(8 20 14 8 8 8) kmax=2
!        seps:             '=, '
!        list_delimiter:   (
!        number:           at least 8
!
!     Outputs:
!
!        number:           8
!        tokens(1):        @block
!        tokens(2):        1
!        tokens(3):        ifree
!        tokens(4):        2
!        tokens(5):        ibc
!        tokens(6):        (8 20 14 8 8 8)
!        tokens(7):        kmax
!        tokens(8):        2
!
!  History:
!
!     12 Mar. 1986  R. Kennelly  Original TOKEN2 variant of TOKENS.
!     15 Jan. 2009  D. Saunders  TOKEN4 derived from TOKEN2 and SCAN4.
!     24 Jan. 2009   "      "    Do NOT convert tokens to upper case.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   character, intent (in)  :: string * (*)   ! String to be parsed for tokens
   character, intent (in)  :: seps   * (*)   ! Delimiter(s) for simple tokens;
                                             ! must not include list_delimiter
   character, intent (in)  :: list_delimiter * 1
                                             ! Delimiter for more general tokens
   integer, intent (inout) :: number         ! Input with the desired or maximum
                                             ! allowed for number of tokens;
                                             ! output with the number found
   character, intent (out) :: tokens (number) * (*)
                                             ! Array of simple/upper case tokens
                                             ! and/or untouched general token(s)
!  Local constants:

!!!character, parameter :: blank * 1 = ' '

!  Local variables:

   integer  :: count, first, last, mark
!!!logical  :: generalized

!  Procedures:

   external :: scan2, scan4

!  Execution:

   first = 1;  last = len_trim (string)
   count = 0

   do while (count < number)

!     Find the delimiting indices of the next token, if any.
!     Treat it as a simple token initially.

!!!   generalized = .false.

      call scan2 (string, seps, first, last, mark)

      if (mark == 0) exit

!     Else check for a generalized token:

      if (string (first : first) == list_delimiter) then

         call scan4 (string, list_delimiter, first, last, mark)

         if (mark == 0) exit  ! No valid token - must be erroneous input

!!!      generalized = .true.
         first = first - 1    ! Retain the delimiters omitted by SCAN4
         mark  = mark  + 1

      end if

      count = count + 1

      tokens (count) = string (first : mark)

!!!   if (.not. generalized) call upcase (tokens (count))

      first = mark + 2

   end do  ! Look for another token

!  Fill the rest of TOKENS with blanks and set NUMBER for output.

!!!tokens (count + 1 : number) = blank  ! Why bother?
   number = count

   end subroutine token4
