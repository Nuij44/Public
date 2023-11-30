module recup_affichage

    use iso_fortran_env

    implicit none
    INTEGER, parameter :: rp = REAL64
contains
    real(rp) function func(x,y)
        real(rp), intent(in) :: x, y
        func = (x - 0.5)**2 + (y - 0.5)**2
        return 
    end function

    subroutine nb_points(nomfic, nb_pt, nb_tri)
        CHARACTER(len=50) :: nomfic, nomficnode, str, nomficele
        INTEGER :: nb_pt, lens, lens2, nb_tri
        
        lens   = INDEX(nomfic,' ')-1   
        str=''
        lens2  = INDEX(str,' ') 
    
        nomficnode = str
        nomficele = str
        nomficnode(lens2:lens2+lens) = nomfic
        nomficnode(lens2+lens:lens2+lens+7)= '.1.node'
        nomficele(lens2:lens2+lens) = nomfic
        nomficele(lens2+lens:lens2+lens+6) = '.1.ele'
    
        OPEN(unit=50,file=nomficnode)
        READ(50,*)nb_pt
        CLOSE(50)
    
        OPEN(unit=50,file=nomficele)
        READ(50,*)nb_tri
        CLOSE(50)
    
    end subroutine nb_points

    subroutine recup_point(nomfic, coord, p, nb_pt, t, nb_tri) 
   
        implicit none
    
        INTEGER          :: n_node, n_ele, Dim, n1, n2, n3, n_per_triangle
        INTEGER          :: i, lens, lens2, stockage
        INTEGER :: nb_pt, nb_tri
        INTEGER, DIMENSION(3,nb_tri) :: t
        INTEGER, DIMENSION(nb_pt) :: p
        REAL(rp), DIMENSION(nb_pt,2) :: coord
        CHARACTER(len=50)                   :: nomfic, nomficnode, nomficele, str
        
        !--------------------------------------------------------!
        !    LECTURE DES FICHIERS CARRE.1.NODE ET CARRE.1.ELE    !
        !--------------------------------------------------------!
    
        lens   = INDEX(nomfic,' ')-1   
        str=''
        lens2  = INDEX(str,' ') 
      
        nomficnode = str
        nomficele = str
        nomficnode(lens2:lens2+lens) = nomfic
        nomficnode(lens2+lens:lens2+lens+7)= '.1.node'
        nomficele(lens2:lens2+lens) = nomfic
        nomficele(lens2+lens:lens2+lens+6) = '.1.ele'
      
        OPEN(unit=50,file=nomficnode)
        READ(50,*)n_node,Dim,n1,n2
        DO i=1,n_node
          READ(50,*) stockage,coord(i,1),coord(i,2),p(i)
        END DO
      
        CLOSE(50)
        
        OPEN(unit=50,file=nomficele)
        READ(50,*)n_ele,n_per_triangle,n3
        DO i=1,n_ele
           READ(50,*)stockage,t(1,i),t(2,i),t(3,i)
        ENDDO
      
        CLOSE(50)
    
      
      END subroutine recup_point

      SUBROUTINE CellVertexVtk(DATA, Mesh, Var, PbName)

        use decla_type
        TYPE(Donnees)   , INTENT(in)     :: DATA
        TYPE(MeshDef)   , INTENT(in)     :: Mesh
        TYPE(variables), INTENT(in) :: Var
        CHARACTER(LEN=70),INTENT(in)     :: PbName
    
    
        CHARACTER(LEN=75) :: vtk=" "
        INTEGER           :: is, jt
        INTEGER           :: lensd3, lPbName
    
    
        WRITE (6,FMT = *) " --------------------------------------------"
        WRITE (6,FMT = *) " ---------      passage dans vtk   ----------"
        WRITE (6,FMT = *) " --------------------------------------------"
    
        ! ECRITURE sous FICHIER vtk !
    
        lPbName                    = INDEX(PbName,' ') - 1
        vtk(1:lPbName)             = PbName(1:lPbName)
        vtk(lPbName+1:lPbName+5)   = '.vtk '
        lensd3                     = lPbName+4
        OPEN(UNIT=61,FILE=vtk(1:lensd3) )
    
        WRITE(61,'(A)')'# vtk DataFile Version 3.0'
        WRITE(61,'(A)')'# Solution du M1 TR'
        WRITE(61,'(A)')'ASCII'
        WRITE(61,'(A)')'DATASET UNSTRUCTURED_GRID'
        WRITE(61,'(A,I7,A)')'POINTS', Mesh%Npoint,'  float'
    
        DO is = 1,Mesh%Npoint
          WRITE(61,'(E13.7,2x,E13.7,2x,E13.7)') Mesh%coor(1,is), Mesh%coor(2,is),  DATA%Z(is)
       END DO
    
        WRITE(61,'(A,1x,I7,1x,I6)')'CELLS',Mesh%Nelemt, 4*Mesh%Nelemt
    
        DO jt = 1,Mesh%Nelemt
           WRITE(61,'(I1,1x,I7,1x,I7,1x,I7)') 3, Mesh%Nu(1,jt)-1, Mesh%Nu(2,jt)-1, Mesh%Nu(3,jt)-1
        END DO
    
        WRITE(61,'(A,1x,I7)')'CELL_TYPES', Mesh%Nelemt
        DO jt = 1,Mesh%Nelemt
           WRITE(61,'(I1)') 5
        ENDDO
    
        WRITE(61,'(A,1x,I7)')'POINT_DATA',Mesh%Npoint
        WRITE(61,'(A)')'SCALARS hauteur float'
        WRITE(61,'(A)')'LOOKUP_TABLE DEFAULT'

        DO is = 1,Mesh%Npoint
          WRITE(61,'(ES20.7)') max(1.E-08,Var%Ua(1,is))
        END DO
    
        WRITE(61,'(A)')'SCALARS topographie float'
        WRITE(61,'(A)')'LOOKUP_TABLE DEFAULT'
    
        DO is = 1,Mesh%Npoint
          WRITE(61,'(ES20.7)') max(1.E-08,DATA%Z(is))
       END DO
    
        CLOSE(61)
    
      END SUBROUTINE CellVertexVtk
      
end module recup_affichage