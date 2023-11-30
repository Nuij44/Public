MODULE fonction

USE iso_fortran_env


CONTAINS

REAL(REAL64) FUNCTION P(rho)
	REAL(REAL64) :: rho

REAL(REAL64) FUNCTION fonction_flux(x)
	REAL(REAL64) :: x
	fonction_flux = 0.5_REAL64 * x**2
END FUNCTION fonction_flux

REAL(REAL64) FUNCTION Lax_Friedrich(UL, UR, dx, dt, flux)
REAL(REAL64), INTENT(in) :: UL, UR, dx, dt
REAL(REAL64), EXTERNAL :: flux

Lax_Friedrich = 0.5_REAL64*(flux(UL) + flux(UR)) - (dx/(2*dt))*(UR - UL)

END FUNCTION

REAL(REAL64) FUNCTION Godunov(U, Uplus, flux)
REAL(REAL64), INTENT(in) :: U, Uplus
REAL(REAL64), EXTERNAL :: flux
REAL(REAL64) :: sigma

IF (U < Uplus) THEN
	IF (0 < U) THEN
		Godunov = flux(U)
	ELSE IF (0 > Uplus) THEN
		Godunov = flux(Uplus)
	ELSE IF ((U < 0) .AND. (0 < Uplus)) THEN
		Godunov = flux(0)
	END IF
ELSE IF (Uplus <= U) THEN
		sigma = (flux(Uplus) - flux(U))/(Uplus - U)
		IF (sigma >= 0) THEN
			Godunov = flux(U)
		ELSE 
			Godunov = flux(Uplus)
		END IF
END IF
END FUNCTION Godunov

END MODULE fonction
