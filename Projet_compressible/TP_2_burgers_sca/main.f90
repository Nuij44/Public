 program main

!module utilisé
use iso_fortran_env
use fonction


IMPLICIT NONE

!déclaration de variable
INTEGER, PARAMETER :: rp = REAL64
INTEGER :: i, Ns
REAL(rp) :: xdeb, xfin, CFL, Tfin, dx, date, dt, a
REAL(rp), DIMENSION(:), ALLOCATABLE :: F, U, Unext		!U stocke les u_i,n et Unext les u_i,n+1


!Lecture du fichier de données
call lecture(xdeb,xfin,Ns,CFL,Tfin)



!allocation de tableau
allocate(U(1:Ns), Unext(1:Ns))
allocate(F(1:Ns - 1))

!On initialise U
call initialisation(U,xdeb,xfin,NS)

!call sauvegarde(U,Ns,xdeb,xfin)
!On commence la boucle en temps
date = 0_REAL64
dx = (xfin - xdeb)/Ns
DO WHILE (date < Tfin)
	!On utilise la CFL pour choisir un dt pour cette itération
	
	!On cherche le min de U pour trouver dt
	a = abs(u(1))
	DO i = 2, Ns
		a = max(abs(U(i)), a)
	END DO
	dt = CFL*dx/a
	dt = min(dt, (Tfin - date))
	date = date + dt
	write(6,*) 'dt_final = ', dt
	
	!calcul du flux 
	DO i = 1, Ns - 1
		F(i) = Lax_Friedrich(U(i), U(i+1), dx, dt, fonction_flux)
		!F(i) = Godunov(U(i), U(i+1), fonction_flux)
	END DO
	
	DO i = 2, Ns-1
		Unext(i) = U(i) -(dt/dx)*(F(i) - F(i-1))
	END DO
	
	
	!Conditions limites
	
	!Dirichlet
	!Unext(1) = U(1)
	!Unext(Ns) = U(Ns)
	
	!Neumann
	Unext(1) = Unext(2)
	Unext(Ns) = Unext(Ns - 1)
	
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

REAL(rp) FUNCTION f_burgers(x)
	REAL(rp), INTENT(in) :: x
	
	f_burgers = (x**2)/2_rp
	
END FUNCTION f_burgers

REAL(rp) FUNCTION U_0(x)
	REAL(rp) ,INTENT(in) :: x
	
	IF (x <= 1) THEN
		U_0 = 0_rp
	ELSE IF (2 <= x) THEN
			U_0 = 2_rp
		ELSE 
			U_0 = 1_rp
	END IF
END FUNCTION U_0


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

DO i = 1, Ns
	x = xdeb + i*(xfin - xdeb)/Ns
	U(i) = U_0(x)
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



end program main
