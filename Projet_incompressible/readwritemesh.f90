subroutine readwritemesh(nomfic) 
   
  INTEGER          :: n_node, n_ele, Dim, n1, n2, n3, n_per_triangle
  INTEGER          :: int, i, lens, lens2, lens3
  INTEGER, DIMENSION(:,:),ALLOCATABLE :: p, t
  REAL(8), DIMENSION(:,:),ALLOCATABLE :: coord
  CHARACTER(len=50)                   :: nomfic, nomficnode, nomficele, str
  
  !--------------------------------------------------------!
  !    LECTURE DES FICHIERS CARRE.1.NODE ET CARRE.1.ELE    !
  !--------------------------------------------------------!
  print*,'nom du fichier a lire'
  read*,nomfic
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
  ALLOCATE(coord(n_node,Dim),p(n_node,n1+n2))
  DO i=1,n_node
     READ(50,*)int,coord(i,:),p(i,:)
  ENDDO

  CLOSE(50)

  OPEN(unit=50,file=nomficele)
  READ(50,*)n_ele,n_per_triangle,n3
  ALLOCATE(t(n_per_triangle+n3,n_ele))
  DO i=1,n_ele
     READ(50,*)int,t(:,i)
  ENDDO

  CLOSE(50)

  !--------------------------------------------------------!
  !    CREATION DU FICHIER DE MAILLAGE CARRE.MESH          !
  !--------------------------------------------------------!
  str = ''
  lens2  = INDEX(str,' ') 
  nomficnode = str
  i = floor(log(n_ele+0._8)/log(10._8)) + 1

  write(str,*)n_ele
  
  lens3 = LEN_TRIM(str)
  nomficnode(lens2:lens2+lens) = nomfic
  nomficnode(lens2+lens:lens2+lens+1) = '.'
  nomficnode(lens2+lens+1:lens2+lens+i+1) = str(lens3-i+1:lens3)
  nomficnode(lens2+lens+i+1:lens2+lens+i+6) = '.mesh'

  print*,'Creation du fichier ',nomficnode

  OPEN(unit=50,file=nomficnode)
  WRITE(50,'(A)')'MeshVersionFormatted 1'
  WRITE(50,'(A)')'Dimension'
  WRITE(50,'(I2)')Dim
  WRITE(50,'(A)')'Vertices'
  WRITE(50,'(I7)')n_node
  DO i=1,n_node
     WRITE(50,'(E19.12,2x,E19.12,2x,I3)')coord(i,1),coord(i,2),p(i,n1+n2)
  ENDDO
  WRITE(50,'(A)')'Triangles'
  WRITE(50,'(I7)')n_ele
  DO i=1,n_ele
     WRITE(50,'(I7,2x,I7,2x,I7,2x,I7)')t(1,i),t(2,i),t(3,i),1
  ENDDO

  CLOSE(50)

END subroutine readwritemesh
