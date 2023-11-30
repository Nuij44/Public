module decla_type
    implicit none



type :: MeshDef
    integer :: Npoint, Nelemt
    integer, Dimension(:,:), POINTER :: Nu
    REAL(8), Dimension(:,:), POINTER :: coor

end type MeshDef

type :: Variables
    REAL(8), Dimension(:,:), POINTER :: Ua 
end type Variables

type :: Donnees
    REAL(8), DIMENSION(:), POINTER :: Z
end type Donnees

end module decla_type