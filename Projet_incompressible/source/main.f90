program main
    use decla_type
    use recup_affichage
    use MEF
    use iso_fortran_env

    implicit none

    INTEGER, parameter :: pr = REAL64
    CHARACTER(len=70) :: nom
    INTEGER :: nb_element, nb_triangle, i, Nk, j
    REAL(pr) :: F1,F2
    real(pr), Dimension(:,:), ALLOCATABLE :: coord
    INTEGER, Dimension(:), ALLOCATABLE :: position
    integer , DIMENSION(:,:), ALLOCATABLE :: connect
    REAL(pr), DIMENSION(:) :: l
    REAL(pr), DIMENSION(:,:) :: A
    type(MeshDef) :: maillage
    type(Donnees) :: result
    type(variables) :: Var


    
    write(6,*) ("Entrez le noms du fichier : ")
    read(5,*) nom

    !On lance la subroutine pour récuperer le nombre de points
    call nb_points(nom, nb_element, nb_triangle)
    
    !write(6,*) nb_element

    !on récupère les coordonnées des points et leur position
    ALLOCATE(coord(nb_element,2), position(nb_element), connect(3,nb_triangle))
    call recup_point(nom, coord, position, nb_element, connect, nb_triangle)
    
    !on construit maintenant les variables pour la subroutine CellVertexVtk
    maillage%Npoint = nb_element
    maillage%Nelemt = nb_triangle
    ALLOCATE(maillage%coor(2,maillage%Npoint), maillage%Nu(3,maillage%Nelemt))
    
    do i = 1, nb_element
        maillage%coor(1,i) = coord(i,1)
        maillage%coor(2,i) = coord(i,2)
    end do
    
    do i = 1, nb_triangle
        maillage%Nu(1,i) = connect(1,i)
        maillage%Nu(2,i) = connect(2,i)
        maillage%Nu(3,i) = connect(3,i)
    end do
    
    ALLOCATE(result%Z(nb_element),Var%Ua(1,nb_element))
    do i = 1, nb_element
        result%Z(i) = func(coord(i,1),coord(i,2))
        Var%Ua(1,i) = result%Z(i)
    end do
    
    call CellVertexVtk(result, maillage, Var, nom)

    !test de la transformation affine

    !call F(0._pr,1._pr,F1,F2,10,nb_triangle,nb_element,coord,connect)
    !write(6,*) 'Les coordonnées du  1er point du trinagle 10 est (',F1,F2,') et la solution exacte est (0.25 ; 0.96875)'
    
    Nk = 3*nb_triangle
    ALLOCATE(A(Nk,Nk), l(Nk))

    !On remplit le vecteur l
    DO i = 1, nb_triangle
        DO j = 1,3
            l(3*i + j) = l_h(phi(i,j))
        END DO
    END DO


    deallocate(coord, position, connect, result%Z, Var%Ua, maillage%Nu, maillage%coor, A, l)
end program main