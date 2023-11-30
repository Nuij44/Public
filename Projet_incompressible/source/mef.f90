module MEF
    use iso_fortran_env
    
    implicit none
    
    integer,parameter :: rp = REAL64
contains
    subroutine F(x,y,F1,F2,K,nb_triangle,nb_sommet,coord,connect)
        REAL(rp) :: x,y,F1,F2
        INTEGER :: K,j,nb_triangle,nb_sommet
        REAL(rp), DIMENSION(2,2) :: Bk
        INTEGER, DIMENSION(3) :: n_global
        INTEGER, DIMENSION(3,nb_triangle) :: connect
        REAL(rp), DIMENSION(nb_sommet,2) :: coord
        
        n_global = connect(:,K)
        DO j = 1,2
            Bk(j,1) = coord(n_global(2),j) - coord(n_global(1),j)
            Bk(j,2) = coord(n_global(3),j) - coord(n_global(1),j)
        END DO
        F1 = Bk(1,1)*x + Bk(1,2)*y + coord(n_global(1),1)
        F2 = Bk(2,1)*x + Bk(2,2)*y + coord(n_global(1),2)
    END subroutine F


    REAL(rp) FUNCTION quad_O1(integrande)
        REAL(rp), external :: integrande

        quad_O1 = (1._rp/2._rp)*(integrande(0._rp,0._rp)+integrande(1._rp,0._rp)+integrande(0._rp,1._rp))
    END FUNCTION quad_O1
end module MEF