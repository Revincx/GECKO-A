MODULE tempoci
IMPLICIT NONE

REAL :: primoOH,secoOH,primosci,ycarbo
INTEGER,PARAMETER :: mxsci=15
INTEGER :: nnsci
REAL :: yysci(mxsci)
CHARACTER(LEN=120) :: chemsci(mxsci)

INTEGER,PARAMETER :: mxcarb=6
CHARACTER(LEN=120):: focarb(mxcarb)
INTEGER           :: nc_carb(mxcarb)
REAL              :: ycarb(mxcarb)

INTEGER :: c1ciflag

END MODULE tempoci
