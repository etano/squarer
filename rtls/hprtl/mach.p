c-------------------------------------------
c
c   cpierleo changes
c
c    define all variables
c
c----------------------------------------
c
c these are the machine-dependant parameters for the HP-UX
c     implicit real*8 (a-h,o-z)
cpierleo
      implicit none
      integer mcpw
      real*8 epsmach
cpierleo
      parameter (mcpw=2,epsmach=1.e-6)

