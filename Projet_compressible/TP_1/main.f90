program main

!module utilisé
use iso_fortran_env


IMPLICIT NONE

!déclaration de variable
INTEGER, PARAMETER :: rp = REAL64
INTEGER :: i, Ns
REAL(rp) :: xdeb, xfin, CFL, Tfin, dx, date, dt, a
REAL(rp), DIMENSION(:), ALLOCATABLE :: F, U, Unext		!U stocke les u_i,n et Unext les u_i,n+1


!Lecture du fichier de données
call lecture(xdeb,xfin,Ns,CFL,Tfin,a)


!allocation de tableau
allocate(U(1:Ns), Unext(1:Ns))
allocate(F(1:Ns - 1))

!On initialise U
call initialisation(U,xdeb,xfin,NS)
call sauvegarde_init(U,Ns,xdeb,xfin)
!on sauvegarde la solution exacte
call solex(Ns,xdeb,xfin, Tfin, a)

!On commence la boucle en temps
date = 0_REAL64
dx = (xfin - xdeb)/Ns
DO WHILE (date < Tfin)
	!On utilise la CFL pour choisir un dt pour cette itération
	dt = dx/a
	dt = min(dt, (Tfin - date))
	date = date + dt
	write(6,*) 'dt = ', dt
	
	!calcul du flux 
	DO i = 1, Ns - 1
		F(i) = Lax_Friedrich(U(i), U(i+1), dx, dt, a)
		!write(6,*) 'F(',i,') = ', F(i)
	END DO
	
	DO i = 2, Ns-1
		Unext(i) = U(i) -(dt/dx)*(F(i) - F(i-1))
	END DO
	
	
	!Conditions limites
	
	!Dirichlet
	!Unext(1) = U(1)
	!Unext(Ns) = U(Ns)
	
	!Neumann
	!Unext(1) = Unext(2)
	!Unext(Ns) = Unext(Ns - 1)
	
	!periodique
	Unext(1) = u(1) - (dt/dx)*(F(1) - F(Ns-1))
	Unext(Ns) = U(1)
	!Mise à jour
	DO i = 1, Ns
		U(i) = Unext(i)
	END DO
	
END DO

!on sauvegarde le resultat dans un fichier result.txt
call sauvegarde(U,Ns,xdeb,xfin)

!Desatribution des tableaux
DEALLOCATE(U,Unext,F)
CONTAINS


subroutine lecture(xdeb,xfin,Ns,CFL,Tfin,a)
INTEGER, PARAMETER :: rp = REAL64
INTEGER :: Ns
REAL(rp) :: xdeb, xfin, CFL, Tfin, a


open(unit=50, file='donnees.dat', form='formatted')
read(50,*) xdeb, xfin
read(50, *) Ns
read(50, *) CFL
read(50, *) Tfin
read(50,*) a
 close(50)

end subroutine lecture



subroutine initialisation(U,xdeb,xfin,NS)
INTEGER, PARAMETER :: rp = REAL64
INTEGER :: Ns, i
REAL(rp), DIMENSION(Ns) :: u
REAL(rp) :: xdeb, xfin, x
REAL(rp) :: pi

pi = acos(-1._REAL64)
DO i=1,int(Ns/2)
	U(i) = 2_rp
END DO
DO i = (int(Ns/2) + 1), Ns
	U(i) = 1_rp
END DO

end subroutine initialisation

subroutine sauvegarde(U,Ns,xdeb,xfin)
INTEGER, PARAMETER :: rp = REAL64
INTEGER :: Ns, i
REAL(rp), DIMENSION(Ns) :: u
REAL(rp) :: xdeb, xfin, x

open(unit=50,file='result.txt')
do i = 1,Ns
	x = xdeb + i*(xfin - xdeb)/Ns
	write(50,*) x, U(i)
end do
 close(50)
end subroutine sauvegarde


subroutine sauvegarde_init(U,Ns,xdeb,xfin)
INTEGER, PARAMETER :: rp = REAL64
INTEGER :: Ns, i
REAL(rp), DIMENSION(Ns) :: u
REAL(rp) :: xdeb, xfin, x
	
open(unit=50,file='result_deb.txt')
do i = 1,Ns
	x = xdeb + i*(xfin - xdeb)/Ns
	write(50,*) x, U(i)
end do
close(50)
end subroutine sauvegarde_init

subroutine solex(Ns,xdeb,xfin, Tfin, v)
INTEGER, PARAMETER :: rp = REAL64
INTEGER :: Ns, i
REAL(rp), DIMENSION(Ns) :: Stock
REAL(rp) :: xdeb, xfin, x, Tfin, v

DO i=1, int(Ns*(0.5_rp + Tfin*v)) 
	Stock(i) = 2_rp
END DO
DO i = (int(Ns*(0.5_rp + Tfin*v) + 1)), Ns
	Stock(i) = 1_rp
END DO

open(unit=50,file='solex.txt')
do i = 1,Ns
	x = xdeb + i*(xfin - xdeb)/Ns
	write(50,*) x, Stock(i)
end do
 close(50)
end subroutine solex

REAL(rp) FUNCTION Lax_Friedrich(Ui, Uinext, dx, dt, v)
REAL(rp), INTENT(in) :: Ui, Uinext, dx, dt, v

Lax_Friedrich = 0.5_rp*(v*Ui + v*Uinext) - (dx/(2*dt))*(Uinext - Ui)

END FUNCTION

end program main
