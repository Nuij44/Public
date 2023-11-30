module routine

use iso_fortran_env

IMPLICIT NONE

CONTAINS


subroutine lecture(xdeb,xfin,Ns,CFL,Tfin)
INTEGER, PARAMETER :: rp = REAL64
INTEGER :: Ns
REAL(rp) :: xdeb, xfin, CFL, Tfin


open(unit=50, file='donnees.dat', form='formatted')
read(50,*) xdeb, xfin
read(50, *) Ns
read(50, *) CFL
read(50, *) Tfin
 close(50)

end subroutine lecture



subroutine initialisation(U,xdeb,xfin,NS)
INTEGER, PARAMETER :: rp = REAL64
INTEGER :: Ns, i
REAL(rp), DIMENSION(Ns) :: u
REAL(rp) :: xdeb, xfin, x
REAL(rp) :: pi

pi = acos(1._REAL64)
DO i=1,Ns
	x = xdeb + i*(xfin - xdeb)/Ns
	U(i) = sin(pi*x)
END DO

end subroutine initialisation

REAL(rp) FUNCTION Lax_Friedrich(Ui, Uinext, dx, dt, v)
REAL(rp), INTENT(in) :: Ui, Uinext, dx, dt, v

Lax_Friedrich = v - (dx/2*dt)*(Uinext - Ui)

END FUNCTION


end module routine
