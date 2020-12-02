module porosity

! porE
! Author Kai Trepte
! Version January 23, 2019
! Version April 29th, 2019    -- implement sub-grid division
! Version August 23rd, 2019   -- add many more elements, vdW and covalent radii
! Version September 4th, 2019 -- read initial information from input file
! Version October 15th, 2019  -- start thinking about pore windows and how to get them -- 
!         November 21th, 2019 -- pore windows: read pore centers and sizes from PSD analysis. 
!                                Then, evaluate vectors between the centers. Minimum distance to 
!                                the vdW surface of these vectors is the pore window!
!         January 08th, 2020  -- Start final implementation of pore windows, include in evaluation
! May 28th, 2020              -- restructuring to make it a python module!
! Sep 14th, 2020              -- added more elements (up to Rn). 
!                             -- Introduce return variable for easier handling with python


contains
subroutine OSA(struct,&
                porosity,density,poreV,V_total,V_occupied,V_overlap) ! Return values
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Evaluate the porosity of a crystal structure. Taking information from xyz file and cell parameters.    
  ! 1. evaluation: calculate total volume, occupied volume (excluding overlap) and void volume -> P = V_void/V_total 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  character(len=*), intent(in)        :: struct                    ! full path to input structure (xyz type)
  character(len=200)                  :: name_struct               ! derive name of the structure from the input file
  ! Evaluation method: Overlapping sphere approach (OSA)
  real(8)                             :: sub_overlap               ! overlap volume as evaluated by a subroutine
  real(8)                             :: distance_ab, new_distance ! distance between two atoms, distance evaluated due to PBC
  integer(8)                          :: number_of_atoms
  real(8)                             :: cell_a(3)           ! array for the cell vector in the a direction. Vector.
  real(8)                             :: cell_b(3)           ! array for the cell vector in the b direction. Vector.
  real(8)                             :: cell_c(3)           ! array for the cell vector in the c direction. Vector.
  real(8), allocatable                :: coordinates(:,:)    ! array for the coordinates. Matrix.
  character(2), allocatable           :: elements(:)         ! array for the elements. Vector.
  real(8)                             :: m_total             ! total mass of unit cell
  real(8)                             :: start, finish       ! evaluate the time
  
  real(8), parameter                  :: pi = 4.0D0*atan(1.0D0)      ! define pi
  real(8), parameter                  ::  u = 1.660539D0             ! define atomic mass unit (in 10**-27 kg)
  ! Evaluation
  integer(8)  :: a,b,c,d,e,n,t                                       ! loop parameter

  ! Return values
  real(8),intent(out) :: porosity, density, poreV
  real(8),intent(out) :: V_total, V_occupied, V_overlap 

  ! Initial value
  name_struct = 'User-defined system'
  open(unit=15,file=struct,status='old',action='read')               ! read in the xyz file
  ! Extract name of file name
  do a = 1, len(struct)-3
    if ((struct(a+1:a+3).eq.'xyz').and.(struct(a:a).eq.'.')) then
      do b = 1, len(struct)-3
        if (struct(b:b)=='/') then
          name_struct = struct(b+1:a-1)
        end if
      end do
    end if
  end do
  !
  ! Read in the number of atoms and cell vectors
  !
  read(unit=15,fmt='(I13.0)') number_of_atoms                        ! first entry is the number of atoms
  read(unit=15,fmt=*) cell_a(1:3), cell_b(1:3), cell_c(1:3)          ! second entry contains the cell vectors. Read them in individually (makes it easier later on)
  !
  ! Store elements and coordinates
  !
  allocate(elements(number_of_atoms))                                ! allocate (number_of_atoms) fields for elements. There is one elements each.
  allocate(coordinates(number_of_atoms,3))                           ! allocate (number_of_atoms) fields for coordinates. There are 3 coordinates per entry. 
  do n = 1,number_of_atoms                                           ! go through all atoms 
    read(unit=15,fmt=*) elements(n), coordinates(n,1:3)              ! storing element and coordinates
  end do
  close(unit=15)                                                     ! close the file (no longer necessary)
  !
  ! Output file
  ! 
  open(unit=19,file='output_OSA',status='unknown',action='write')
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calculate the total unit cell volume using the triple product V = a . (b x c) . In A^3 !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call cpu_time(start)                                                                         ! initialize time measurement

  V_total = (cell_a(1)*cell_b(2)*cell_c(3) + cell_a(2)*cell_b(3)*cell_c(1) + cell_a(3)*cell_b(1)*cell_c(2) - &
             cell_a(3)*cell_b(2)*cell_c(1) - cell_a(1)*cell_b(3)*cell_c(2) - cell_a(2)*cell_b(1)*cell_c(3))
  !
  ! Initialization
  !
  V_occupied  = 0.0D0                                                     ! initial value for the occupied volume
  m_total     = 0.0D0                                                     ! and for the total mass
  do t = 1, number_of_atoms                                               ! go through all atoms and evaluate the total occupied volume and the total mass
    call eval_vol_mass(elements(t),V_occupied,m_total)
  end do
  !
  ! Calculate overlap volume. Include PBCs
  !
  V_overlap = 0.0D0                                                       ! initial value for the overlap volume
  do a = 1,number_of_atoms-1                                              ! go through all pairs of atoms (first loop til number_of_atoms - 1)
    do b = a+1,number_of_atoms                                            ! second loop from next atom til the end. Evaluate overlap. Consider periodic boundary conditions (PBCs)! No double counting
      distance_ab = sqrt(sum((coordinates(a,:) - coordinates(b,:))**2))   ! Initial distance between two atoms. No PBCs yet.
      do c = 1,3                                                          ! PBCs in all direction. Here for cell_a (-1,0,+1)
        do d = 1,3                                                        ! here for cell_b
          do e = 1,3                                                      ! here for cell_c. Taking all surrounding unit cells into account
            new_distance = sqrt(sum((coordinates(a,:)-coordinates(b,:)+(c-2)*cell_a(:)+(d-2)*cell_b(:)+(e-2)*cell_c(:))**2))
            if (new_distance < distance_ab) then
             distance_ab = new_distance
            end if
          end do
        end do
      end do

      sub_overlap = 0.0D0
      call eval_overlap(elements(a), elements(b), distance_ab, sub_overlap)                     ! use subroutine to analyze the overlap volume
      V_overlap = V_overlap + sub_overlap

      if (sub_overlap > 0.0D0) then
        write(6,666) a,elements(a),b,elements(b),'  d =', distance_ab,' A    with V_overlap = ', sub_overlap,' A^3'
        write(19,666) a,elements(a),b,elements(b),'  d =', distance_ab,' A    with V_overlap = ', sub_overlap,' A^3'
      end if
    end do
  end do
  !  
  ! Write to screen
  !
  write(6,*) name_struct
  write(6,fmt='(1X a,f10.3,1X a)') 'V_total     = ',V_total,'A^3'
  write(6,fmt='(1X a,f10.3,1X a)') 'V_vdW,atoms = ',V_occupied,'A^3'
  write(6,fmt='(1X a,f10.3,1X a)') 'V_overlap   = ',V_overlap,'A^3'
  write(6,fmt='(1X a,f10.3,1X a)') 'V_occupied  = ',V_occupied - V_overlap,'A^3'
  write(6,fmt='(1X a,f10.3,1X a)') 'V_void      = ',V_total - (V_occupied - V_overlap),'A^3'
  write(6,fmt='(1X a,f10.3,1X a)') 'Porosity    = ',(V_total - (V_occupied - V_overlap))/V_total*100,' %'
  write(6,fmt='(1X a,f10.3,1X a)') 'Density of the structure is (m_total/V_total)   : ',m_total*u/V_total*10**3,'kg/m^3'
  write(6,fmt='(1X a,f10.3,1X a)') 'Pore volume density is (V_void/m_total)         : ' &
             ,(V_total - (V_occupied - V_overlap))/(m_total*u)*10**(0),'cm^3/g'
  write(6,fmt='(1X a,f10.3,1X a)') 'Mass of the unit cell is                        : ',m_total*u,'10**-27 kg'
  call cpu_time(finish)
  write(6,fmt='(A,2X,F12.3,1X,A)') 'Total CPU time: ',finish-start,'s'
  !
  ! Write to file
  !
  write(19,*) name_struct
  write(19,fmt='(1X a,f10.3,1X a)') 'V_total     = ',V_total,'A^3'
  write(19,fmt='(1X a,f10.3,1X a)') 'V_vdW,atoms = ',V_occupied,'A^3'
  write(19,fmt='(1X a,f10.3,1X a)') 'V_overlap   = ',V_overlap,'A^3'
  write(19,fmt='(1X a,f10.3,1X a)') 'V_occupied  = ',V_occupied - V_overlap,'A^3'
  write(19,fmt='(1X a,f10.3,1X a)') 'V_void      = ',V_total - (V_occupied - V_overlap),'A^3'
  write(19,fmt='(1X a,f10.3,1X a)') 'Porosity    = ',(V_total - (V_occupied - V_overlap))/V_total*100,' %'
  write(19,fmt='(1X a,f10.3,1X a)') 'Density of the structure is (m_total/V_total)   : ',m_total*u/V_total*10**3,'kg/m^3'
  write(19,fmt='(1X a,f10.3,1X a)') 'Pore volume density is (V_void/m_total)         : ' &
             ,(V_total - (V_occupied - V_overlap))/(m_total*u)*10**(0),'cm^3/g'
  write(19,fmt='(1X a,f10.3,1X a)') 'Mass of the unit cell is                        : ',m_total*u,'10**-27 kg'
  write(19,fmt='(A,2X,F12.3,1X,A)') 'Total CPU time: ',finish-start,'s'
  close(19)

  deallocate(elements)
  deallocate(coordinates)

  ! 
  ! Define return values
  !
  porosity = (V_total - (V_occupied - V_overlap))/V_total*100.0D0
  density  = m_total*u/V_total*10.0D0**3
  poreV    = (V_total - (V_occupied - V_overlap))/(m_total*u)*10.0D0**(0)
  return

666 format(I5.0,1X,A,I5.0,1X,A,A,F10.5,A,F10.5,A)
end subroutine OSA  



subroutine GPA_FullGrid(struct,probe_r,grid_a,grid_b,grid_c,&
                poro_void,poro_acc,density,poreV_void,poreV_acc) ! return variables
  ! Use the GPA, providing the total grid points per cell vector
  character(len=*), intent(in)   :: struct                    ! full path to input structure (xyz type)

  real(8), intent(in)            :: probe_r                   ! probe radius
  integer, intent(in)            :: grid_a, grid_b, grid_c

  ! dummy arguments for return values
  real(8), intent(out)           :: poro_void,poro_acc,density,poreV_void,poreV_acc

  ! Call actual GPA, with the given grid
  call do_GPA(struct,probe_r,grid_a,grid_b,grid_c,&
          poro_void,poro_acc,density,poreV_void,poreV_acc)
  return
end subroutine GPA_FullGrid


subroutine GPA_GridPerA(struct,probe_r,g,&
                poro_void,poro_acc,density,poreV_void,poreV_acc) ! return values
  ! Use the GPA, providing the approximate number of points per angstrom
  character(len=*), intent(in)   :: struct                    ! full path to input structure (xyz type)

  real(8), intent(in)            :: probe_r, g        ! probe radius
  integer                        :: grid_a, grid_b, grid_c
  real(8)                        :: cell_a(3), cell_b(3), cell_c(3)

  ! dummy arguments for return values
  real(8), intent(out)           :: poro_void,poro_acc,density,poreV_void,poreV_acc

  open(unit=15,file=struct,status='old',action='read')                     ! read in the xyz file
  read(unit=15,fmt=*) 
  read(unit=15,fmt=*) cell_a(1:3), cell_b(1:3), cell_c(1:3)                ! second entry contains the cell vectors. Read them in individually (makes it easier later on)
  close(unit=15)
  !
  ! For number of grid point per Angtrom
  grid_a = ceiling(g*sqrt(cell_a(1)**2 + cell_a(2)**2 + cell_a(3)**2))     ! ceiling -> round to the next higher integer. Use g * len_unit_cell_vector as the number of grid points
  grid_b = ceiling(g*sqrt(cell_b(1)**2 + cell_b(2)**2 + cell_b(3)**2))
  grid_c = ceiling(g*sqrt(cell_c(1)**2 + cell_c(2)**2 + cell_c(3)**2))

  ! Call actual GPA, with the given grid
  call do_GPA(struct,probe_r,grid_a,grid_b,grid_c,&
          poro_void,poro_acc,density,poreV_void,poreV_acc)
  return
end subroutine GPA_GridPerA



subroutine do_GPA(struct,probe_r,grid_a,grid_b,grid_c,&
                poro_void,poro_acc,density,poreV_void,poreV_acc) ! return values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Script to evaluate the porosity of a crystal structure. Taking information from xyz file and cell parameters.                                 !
! 2. evaluation: place grid inside the unit cell, count each point which is in an occupied region, compare to total points -> P = N_occ/N_tot   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  character(len=*), intent(in)                :: struct                    ! full path to input structure (xyz type)
  character(len=200)                          :: name_struct               ! derive name of the structure from the input file

  real(8), intent(in)     :: probe_r                                       ! probe radius
  integer, intent(in)     :: grid_a, grid_b, grid_c                        ! number of grid points along cell vectors (a,b,c), loop variables to write grid points, 

  real(8)     :: grid_point_x, grid_point_y, grid_point_z, factor          ! grid point coordinates for any specific grid point (x,y,z), grid density (grid points per A^3)
  integer     :: aa, bb, cc, running_n                                     ! loop variables to write grid points, 
                                                                           ! running variable for the loop (assign grid_points correctly)
  integer(8)  :: n_coords                                                  ! loop counter for the coordinates (coords)
  integer(8)  :: counter_access, counter_noOccu                            ! counter to evaluate whether a point is accessible, or NOT occupied
  integer(8)  :: n_access, n_occ, n_check_acc, n_noOccu                    ! counter for list assignment -> accessible, occupied, counter for accessibility check list, not occupied
  real(8)     :: dist_point_atom, new_point_atom, dist_point_point         ! distance from a grid point to an atom, distance evaluated due to PBC, distance between grid points
  real(8)     :: V_void, V_accessible                                      ! void and accessible volume
  real(8)     :: grid_per_A_x, grid_per_A_y, grid_per_A_z                  ! grid per angstrom, in the cell vector directions (x == a, y == b, z == c)
  real(8), allocatable, dimension(3)  :: grid_points(:,:)                  ! array for the grid_points. Matrix.
  real(8), allocatable, dimension(3)  :: list_noOccu(:,:)                  ! empty list for all NOT occupied points
  real(8), allocatable, dimension(3)  :: list_access(:,:)                  ! empty list for all accessible points. Initial evaluation
  real(8), allocatable, dimension(3)  :: list_check_acc(:,:)               ! empty list for the check of accessibility. Will be smaller than list_access and thus easier/faster to evaluate
  !! For pore windows
  logical                             :: file_exist                        ! Check whether output_PSD file exists
  real(8), allocatable                :: pore_center(:,:)                  ! empty list for the coordinates of the pore centers
  real(8), allocatable                :: pore_size(:)                      ! empty list for the pore sizes
  real(8), allocatable                :: pore_windows(:)                   ! empty list for the pore windows
  integer(8)  :: counter_1, n_pore, ios, n_points                          ! counter for evaluation of pore window/ pore size, counter for amount of pore sizes, iostat for reading, number of points
  real(8)                             :: junk2                             ! junk for reading
  real(8)                             :: dist1                             ! distance for evaluation
  real(8)                             :: tmp_vec(3),tmp_vec2(3)            ! temporary vector 
  ! General parameters
  integer(8)                          :: number_of_atoms
  real(8)                             :: cell_a(3)           ! array for the cell vector in the a direction. Vector.
  real(8)                             :: cell_b(3)           ! array for the cell vector in the b direction. Vector.
  real(8)                             :: cell_c(3)           ! array for the cell vector in the c direction. Vector.
  real(8), allocatable                :: coordinates(:,:)    ! array for the coordinates. Matrix.
  character(2), allocatable           :: elements(:)         ! array for the elements. Vector.
  real(8)                             :: V_total, m_total    ! total volume of the cell, total mass of unit cell
  real(8)                             :: V_occupied          ! total occupied volume of the atoms in the cell
  
  real(8)                             :: start, finish       ! evaluate the time
  
  real(8), parameter                  :: pi = 4.0D0*atan(1.0D0)      ! define pi
  real(8), parameter                  ::  u = 1.660539D0             ! define atomic mass unit (in 10**-27 kg)
  
  ! Evaluation
  integer(8)  :: a,b,c,d,e,f,n,t,v,w,x                                     ! loop parameter
  
  ! Some arrays for easier handling of vdW radii
  character(2)              :: all_pse(86)                        ! all currently available elements
  real(8)                   :: all_vdW_radii(86)                  ! the respective vdW radii
  integer                   :: all_elements                       ! number of all different elements
  character(2), allocatable :: tmp_pse(:)                         ! temporary list to evaluate the used elements
  character(2), allocatable :: pse(:)                             ! used elements
  real(8), allocatable      :: vdW_radii(:)                       ! used vdW radii
  integer(8)                :: no_elements                        ! number of different elements
  !
  ! IDEA: Subdivide grid into sub-grids per atom. Then, only evaluate the necessary grid points.
  ! Thus, exclude grid points which are very far away. 
  ! - Introduce new array, dimensions are number_of_atoms AND a specific number according to the specific atom (needs to be counted, need two loops over the grid points. Should still be cheap)
  ! - When looping, choose grid point for an atom if it is close to this grid point -> store grid point for the sub-grid for the atom 
  !
  ! For sub-grid generation, one needs a variable array which is allocatable
  type global_array
    integer(8)           :: sub_grid_points                       ! dimensions: n_atoms
    real(8), allocatable :: sub_grids(:,:)                        ! dimensions: n_atoms, sub_grid_points(atom), 3
  end type global_array
  type(global_array), dimension(:), allocatable :: sub_division

  ! Return values
  real(8), intent(out)  :: poro_void,poro_acc,density,poreV_void,poreV_acc

  ! Define some elements with their vdW_radii
  all_elements = 86
  all_pse = (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',&
               'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',&
               'Co', 'Ni', 'Cu', 'Zn', 'Zr', &
               'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', &
               'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', &
               'Rb', 'Sr', 'Y ',       'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', &
               'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', &
               'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', & 
               'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', &
               'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn' /)
  all_vdW_radii = (/ 1.20D0, 1.40D0, 1.82D0, 1.53D0, 1.92D0, 1.70D0, 1.55D0, 1.52D0, 1.47D0, 1.54D0,&
                     2.27D0, 1.73D0, 1.84D0, 2.10D0, 1.80D0, 1.80D0, 1.75D0, 1.88D0, 2.75D0, 2.31D0,&
                     1.92D0, 1.63D0, 1.40D0, 1.39D0, 2.36D0, &
                     2.30D0, 2.15D0, 2.05D0, 2.05D0, 2.05D0, 2.05D0, &
                     2.10D0, 2.10D0, 2.05D0, 1.90D0, 1.90D0, 2.02D0, &
                     2.90D0, 2.55D0, 2.40D0,         2.15D0, 2.10D0, 2.05D0, 2.05D0, 2.00D0, 2.05D0, 2.10D0, 2.20D0, &
                     2.20D0, 2.25D0, 2.20D0, 2.10D0, 2.10D0, 2.16D0, &
                     3.00D0, 2.70D0, 2.50D0, 2.25D0, 2.20D0, 2.10D0, 2.05D0, 2.00D0, 2.00D0, 2.05D0, 2.10D0, 2.05D0, &
                     2.48D0,2.47D0,2.45D0,2.43D0,2.42D0,2.40D0,2.38D0,2.37D0,2.35D0,2.33D0,2.32D0,2.30D0,2.28D0,2.27D0, &
                     2.20D0, 2.30D0, 2.30D0, 2.00D0, 2.00D0, 2.00D0 /)
 ! for > Ca (besides Co, Ni, Cu, Zn, Zr)  : https://www.re3data.org/repository/r3d100000013


  ! Initial value
  name_struct = 'User-defined system'
  ! Read coordinate file
  open(unit=15,file=struct,status='old',action='read')
  ! Extract name of file name
  do a = 1, len(struct)-3
    if ((struct(a+1:a+3).eq.'xyz').and.(struct(a:a).eq.'.')) then
      do b = 1, len(struct)-3
        if (struct(b:b)=='/') then
          name_struct = struct(b+1:a-1)
        end if
      end do
    end if
  end do
  !
  ! Read in the number of atoms and cell vectors
  !
  read(unit=15,fmt='(I13.0)') number_of_atoms                               ! first entry is the number of atoms
  read(unit=15,fmt=*) cell_a(1:3), cell_b(1:3), cell_c(1:3)                 ! second entry contains the cell vectors. Read them in individually (makes it easier later on)
  !
  ! Store elements and coordinates
  !
  allocate(elements(number_of_atoms))                                                         ! allocate (number_of_atoms) fields for elements. There is one elements each.
  allocate(coordinates(number_of_atoms,3))                                                    ! allocate (number_of_atoms) fields for coordinates. There are 3 coordinates per entry. 
  allocate(tmp_pse(number_of_atoms))                                                          ! allocate tmp_pse, dummy
  tmp_pse(:) = 'X'
  no_elements = 0                                                                             ! number of different atoms
  
  do n = 1,number_of_atoms                                                                      ! go through all atoms 
    read(unit=15,fmt=*) elements(n), coordinates(n,1:3)                                         ! storing element and coordinates
    !
    ! Determine what kind of different atoms there are -> use later for the evaluation of the vdW radii
    !
    c = 0                                                                                       ! counter for each entry 
    do a = 1, number_of_atoms
      if (elements(n) == tmp_pse(a)) then                                                       ! if the elements is already in the tmp_pse list -> ignore
        c = 1
      end if
    end do
    if (c == 0) then                                                                            ! if the element is new -> put it into the new list
      no_elements = no_elements + 1
      tmp_pse(no_elements) = elements(n)
    end if
  end do
  close(unit=15)                                                                                ! close the file (no longer necessary)
  !
  ! Write the correct number of atoms and their vdW radii for alter evaluations
  !
  allocate(pse(no_elements))
  allocate(vdW_radii(no_elements))
  do a = 1, all_elements                                                                        ! go through all possible elements
    do b = 1, no_elements                                                                       ! go through the ones that are actually there
      if (all_pse(a) == tmp_pse(b)) then                                                        ! once a new elements has been found
        pse(b)       = all_pse(a)
        vdW_radii(b) = all_vdW_radii(a)
      end if
    end do
  end do
  deallocate(tmp_pse)
  !
  ! Output file
  ! 
  open(unit=19,file='output_GPA',status='unknown',action='write')
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calculate the total unit cell volume using the triple product V = a . (b x c) . In A^3 !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  V_total = (cell_a(1)*cell_b(2)*cell_c(3) + cell_a(2)*cell_b(3)*cell_c(1) + cell_a(3)*cell_b(1)*cell_c(2) - &
             cell_a(3)*cell_b(2)*cell_c(1) - cell_a(1)*cell_b(3)*cell_c(2) - cell_a(2)*cell_b(1)*cell_c(3))
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!! Evaluation method 2 : Uniform grid. Evaluate each point as to whether it is occupied, void !!
  !!!!!                       or accessible                                                        !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call cpu_time(start)                                                                         ! initialize time measurement
  !
  ! Evaluate the total mass (for later evaluation)
  !
  V_occupied  = 0.0D0                                                                          ! Dummy value. Necessary to use subroutine
  m_total     = 0.0D0                                                                          ! initial value for the total mass
  do t = 1,number_of_atoms                                                                     ! go through all atoms and evaluate the total mass
    call eval_vol_mass(elements(t),V_occupied,m_total)
  end do
  !
  ! allocate all lists. Use maximum grid points for each list, as it is not clear how much is needed
  !
  allocate(grid_points(grid_a*grid_b*grid_c,3))                                                ! allocate (grid_size) fields for the grid points. There are 3 coordinates per entry.
  allocate(list_access(grid_a*grid_b*grid_c,3))                                                ! allocate (grid_size) fields for the accessible list. There are 3 coordinates per entry.
  allocate(list_check_acc(grid_a*grid_b*grid_c,3))                                             ! allocate (grid_size) fields for the 'check accessibility' list. There are 3 coordinates per entry.
  allocate(list_noOccu(grid_a*grid_b*grid_c,3))                                                ! allocate (grid_size) fields for the NOT occupied list. There are 3 coordinates per entry.
  !
  ! use these to assign at which point of the respective lists something shall be stored
  !
  n_access = 0
  n_occ = 0
  n_noOccu = 0
  n_check_acc = 0
  !
  !
  ! Write the array of the grid points
  !
  running_n = 0                                                                                ! initialize running variable (for the assignment of the grid_points array)
  do aa = 1,grid_a                                                                             ! go through the grid points and write grid points according to the cell vectors
    do bb = 1,grid_b
      do cc = 1,grid_c
        running_n = running_n + 1                                                                  ! increase running variable at each step
        grid_point_x = cell_a(1)/grid_a*(aa-1) + cell_b(1)/grid_b*(bb-1) + cell_c(1)/grid_c*(cc-1) ! x coordinate of grid point. Choose e.g. aa-1 to include the origin (0,0,0)
        grid_point_y = cell_a(2)/grid_a*(aa-1) + cell_b(2)/grid_b*(bb-1) + cell_c(2)/grid_c*(cc-1) ! y coordinate
        grid_point_z = cell_a(3)/grid_a*(aa-1) + cell_b(3)/grid_b*(bb-1) + cell_c(3)/grid_c*(cc-1) ! z coordinate
        grid_points(running_n,:) = (/ grid_point_x,grid_point_y,grid_point_z /)                    ! assign the respective values to the array
        !
        ! Include evaluation of Occupied or Accessible here, in the grid generation -> no need to go through the entire grid again!
        !
        counter_access = 0                                                                         ! initialize variable for each grid point
        counter_noOccu = 0                                                                         ! initialize variable for each grid point
        loop14: do n_coords = 1,number_of_atoms                                                    ! go through all atoms and evaluate grid points
          dist_point_atom = sqrt(sum((grid_points(running_n,:) - coordinates(n_coords,:))**2))    ! initial distance between grid point and atom
          do c = 1,3                                                                               ! PBCs in all direction. Here for cell_a (-1,0,+1)
            do d = 1,3                                                                             ! here for cell_b
              do e = 1,3                                                                           ! here for cell_c. Taking all surrounding unit cells into account
                new_point_atom = sqrt(sum((grid_points(running_n,:) - coordinates(n_coords,:) + &
                                    (c-2)*cell_a(:) + (d-2)*cell_b(:) + (e-2)*cell_c(:))**2))      ! evaluate new distance due to PBC
                if (new_point_atom < dist_point_atom) then                                         ! if distance is smaller -> use this one !
                  dist_point_atom = new_point_atom
                end if
              end do
            end do
          end do
          !
          ! Evaluate whether point is occupied or accessible
          !
          loop66: do n = 1, no_elements
            if (elements(n_coords) == pse(n)) then
              if (dist_point_atom <= vdW_radii(n)) then                            ! if the grid point is inside any atom (distance is smaller than the vdW radius of the respective atom)
                n_occ = n_occ + 1                                                  ! increase assignemnt counter for the occupied list
                exit loop14                                                        ! stop looping through the atoms at this point. The point is already determined as occupied. Avoid double counting!
              else if (dist_point_atom >= vdW_radii(n) + probe_r) then             ! if grid point is outside an atom + the probe radius -> Immediately accessible
                counter_access = counter_access + 1                                ! add +1 to the counter 'counter_acc' AND to counter_noOccu
                counter_noOccu = counter_noOccu + 1
                exit loop66
              else                                                                 ! If outside an atom -> not occupied
                counter_noOccu = counter_noOccu + 1
                exit loop66
              end if
            end if
          end do loop66
        end do loop14
        !
        !
        !
        if (counter_access == number_of_atoms) then                                                ! if the counter for the accessible points increased for all atoms -> add to list
          n_access = n_access + 1                                                                  ! increase assignment counter for the accessible list
          list_access(n_access,:) = grid_points(running_n,:)
        end if
        if ((counter_noOccu == number_of_atoms) .and. (counter_access .ne. number_of_atoms)) then  ! if the counter for the 'NOT occupied' points increased for all atoms -> add to list
          n_noOccu = n_noOccu + 1                                                                  ! increase assignment counter for the list
          list_noOccu(n_noOccu,:) = grid_points(running_n,:)
        end if
      end do
    end do
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!! PORE WINDOW EVALUATION !!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 
  ! Check whether there is an output_PSD file. If so -> get pore windows
  !
  inquire(file='output_PSD',exist=file_exist)
  if (file_exist) then
    ! 
    ! Read pore sizes from output_PSD
    ! Go through the file. See how many pores there are.  
    ! Allocate pore_center and pore_size accordingly
    !
    n_pore = 0
    open(unit=20,file='output_PSD',status='old',action='read')  
    do a = 1, 100  ! no more than 100 lines in the file.
      read(20,*,iostat=ios) junk2
      !
      ! reading a floating point number was successful
      !
      if (ios==0) then
        !
        ! increase counter for pores
        !
        n_pore = n_pore + 1
      end if
    end do
    !
    ! Go back to the top of the file
    !
    rewind(20)
    !
    ! allocate arrays
    !
    allocate(pore_size(n_pore))
    allocate(pore_center(n_pore,3))
    if (n_pore == 1) then
      allocate(pore_windows(1)) ! for later
    else
      allocate(pore_windows(n_pore*n_pore)) ! for later
    end if
    !
    ! Go through the file again
    !
    n_pore = 0
    do a = 1, 100
      read(20,*,iostat=ios) junk2                          
      !
      ! reading a floating point number was successful
      !
      if (ios==0) then
        !
        ! go up one line
        !
        backspace(20)
        !
        ! increase counter for pores
        !
        n_pore = n_pore + 1
        !
        ! read pore size and pore center
        !
        read(20,*) pore_size(n_pore), junk2, pore_center(n_pore,:)
        !
        ! store the radius, not the diameter
        !
        pore_size(n_pore) = pore_size(n_pore)/2.0D0
      end if
    end do
    close(20)

    write(6,*) 'Pore  Pore radius  Pore diameter       coordinates of pore center'
    write(19,*) 'Pore  Pore radius  Pore diameter       coordinates of pore center'
    do a = 1, n_pore
      write(6,fmt='(I5,F13.4,F15.4,7X,3F11.6)') a, pore_size(a), pore_size(a)*2.0D0, pore_center(a,:)
      write(19,fmt='(I5,F13.4,F15.4,7X,3F11.6)') a, pore_size(a), pore_size(a)*2.0D0, pore_center(a,:)
    end do

    !
    ! initialize all with 5000.0
    ! 
    pore_windows(:) = 5000.0D0
    counter_1 = 0

    !!!!!!!!!!!!!!!!!!!!!!!!!
    ! Pore window evaluation
    !
    ! Use the pore centers as a starting point. Construct vectors between the pore sizes (periodically)
    ! along these vectors, analyze the minimum distance the the vdW surface. 
    ! the smallest value is the pore window!
    ! Use only the vector corresponding to the shortestdistance between two pore centers (take pbcs into account) -> that should do!
    !
    !
    ! Do for all pores, i.e. for b == a as well 
    !                   ! ORG:
    do a = 1, n_pore    ! 1, npore-1
      do b = a, n_pore  ! a+1,n_pore
        dist1 = 1000.0D0 
        do c = 1,3                                                                               ! PBCs in all direction. Here for cell_a (-1,0,+1)
          do d = 1,3                                                                             ! here for cell_b
            do e = 1,3                                                                           ! here for cell_c. Taking all surrounding unit cells into account
              tmp_vec(:) = pore_center(b,:) - (pore_center(a,:) + (c-2)*cell_a(:) + (d-2)*cell_b(:) + (e-2)*cell_c(:))
              if ((sqrt(sum(tmp_vec(:)**2)) < dist1).and.(sqrt(sum(tmp_vec(:)**2)).ne.0.0D0)) then  ! when evaluating the same pore -> exclude zero distances
                !
                ! store new minimum distance and the corresponding values for the cell vectors (periodicity)
                !
                dist1 = sqrt(sum(tmp_vec(:)**2))
                v = c
                w = d
                x = e
              end if
            end do
          end do
        end do
        !
        ! Evaluate pore window
        ! If it is zero -> disregard (pores are not actually adjacent to each other), thus the vector goes through some atoms
        !   can we exclude pores like that before? Would be good
        !
        !
        ! Move along the connecting vector between pore a and pore b. Analyze the distance to the vdW surface. Minimum distance is the
        ! pore window!
        !
        ! Use length of vector between pores to determine number of points to evaluate. e.g. 10 points/A
        ! 
        n_points = int(dist1*10.0D0)  ! 10 points/A to analyze. Should be enough
        dist1 = 100000.0D0            ! initial value for each pore 
        do f = 1, n_points            ! steps along the vector
          !
          ! point to evaluate. Start at pore a. Move step by step further towards pore b.
          !
          tmp_vec(:) = pore_center(a,:) + real(f/real(n_points,8),8)*&
          &(pore_center(b,:) - (pore_center(a,:) + (v-2)*cell_a(:) + (w-2)*cell_b(:) + (x-2)*cell_c(:)))
          !
          ! go through all atoms
          !
          do n_coords = 1, number_of_atoms
            !
            ! For vdW radii
            !
            loop321: do n = 1, no_elements
              if (elements(n_coords) == pse(n)) then
                !
                ! evaluate minimum distance to vdW surface
                !
                dist_point_atom = sqrt(sum((tmp_vec(:) - coordinates(n_coords,:))**2)) - vdW_radii(n)    ! initial distance between grid point and atom
                do c = 1,3                                                                               ! PBCs in all direction. Here for cell_a (-1,0,+1)
                  do d = 1,3                                                                             ! here for cell_b
                    do e = 1,3                                                                           ! here for cell_c. Taking all surrounding unit cells into account
                      new_point_atom = sqrt(sum((tmp_vec(:) - coordinates(n_coords,:) + &
                                          (c-2)*cell_a(:) + (d-2)*cell_b(:) + (e-2)*cell_c(:))**2)) - vdW_radii(n)      ! evaluate new distance due to PBC
                      if (new_point_atom < dist_point_atom) then                                         ! if distance is smaller -> use this one !
                        dist_point_atom = new_point_atom
                      end if
                    end do
                  end do
                end do
                exit loop321
              end if
            end do loop321
            !
            ! if distance is smaller than previous distance
            ! -> store distance and position
            !
            if (dist_point_atom < dist1) then
              tmp_vec2(:) = tmp_vec(:)
              dist1 = dist_point_atom
            end if
          end do
        end do
        !
        ! If dist1 > 0 -> evaluation here
        !
        if (dist1 > 0.0D0) then
          !
          ! Check: ->
          ! Is the distance to the pore centers is very different from the 
          ! pore radii ->> this is not a pore window
          ! Thus, do not take this into account
          !
          ! Distances to the two pore centers, a and b
          !
          dist_point_atom = sqrt(sum((tmp_vec2(:) - pore_center(a,:))**2))
          new_point_atom  = sqrt(sum((tmp_vec2(:) - (pore_center(a,:)+&
       &(pore_center(b,:)-(pore_center(a,:)+(v-2)*cell_a(:)+(w-2)*cell_b(:)+(x-2)*cell_c(:)))))**2))
          !
          ! Compare to pore sizes. If they are very different -> NOT A PORE WINDOW
          ! Check whether distance is within 30% of the pore size
          ! To be checked whether thsi is accurate
          !
          if ((abs((dist_point_atom-pore_size(a))/pore_size(a)) < 0.30D0).and.&
          &    (abs((new_point_atom-pore_size(b))/pore_size(b)) < 0.30D0)) then
            !
            ! If this is all true -> pore window!
            !
            write(6,fmt='(A,I3,A,I3,A,F10.5,A)') 'PORE WINDOW between pore                     ',a,' and pore ',b,' is ',dist1,' A'
            write(19,fmt='(A,I3,A,I3,A,F10.5,A)') 'PORE WINDOW between pore                     ',a,' and pore ',b,' is ',dist1,' A'
            !
            ! Only store values which are physically meaningful (everything larger than H)
            !
            if (dist1 > 1.20D0) then
              counter_1 = counter_1 + 1
              pore_windows(counter_1) = dist1
            end if
          !
          ! In any other case -> not a pore window. Just print window size in between the pores 
          !
          else
            write(6,fmt='(A,I3,A,I3,A,F10.5,A)') 'Minimum radius (NO pore window) between pore ',a,' and pore ',b,' is ',dist1,' A'
            write(19,fmt='(A,I3,A,I3,A,F10.5,A)') 'Minimum radius (NO pore window) between pore ',a,' and pore ',b,' is ',dist1,' A'
          end if
        end if
      end do
    end do
    !
    ! END pore windows
    !
  else
    write(6,*) 'No output_PSD found => no evaluation of pore windows'
    write(19,*) 'No output_PSD found => no evaluation of pore windows'
  end if




  !
  ! for subdivision of the grid
  !
  allocate(sub_division(number_of_atoms))
  !
  ! Intitialize all point as zero
  !
  sub_division(:)%sub_grid_points = 0
  !
  ! first loop -> count points. Then allocate sub_grids and write points into them
  !

  !
  ! Get all points to be double check for accessibility -> don't go through all N_check_acc points!
  ! Get a factor to determine which points shall be taken -> the denser the grid, the smaller the additional distance (epsilon)
  !
  grid_per_A_x = real(grid_a)/sqrt(cell_a(1)**2 + cell_a(2)**2 + cell_a(3)**2)
  grid_per_A_y = real(grid_b)/sqrt(cell_b(1)**2 + cell_b(2)**2 + cell_b(3)**2)
  grid_per_A_z = real(grid_c)/sqrt(cell_c(1)**2 + cell_c(2)**2 + cell_c(3)**2)
  factor = 1.0D0 + 1.0D0/((grid_per_A_x+grid_per_A_y+grid_per_A_z)/3.0D0)                            ! divide by average grid points per A

  !
  ! Get number of sub-grid points per atom
  !
  !
  ! HERE: Idea -> include the firstloop here into the very first loop over the grid points -> avoid looping over all accessible points twice
  !! I.e. do the counting part in the very first loop. Then, in the second loop write grid points
  !
  do n = 1, n_access                                                                           ! go through all accessible points
    loop13: do n_coords = 1,number_of_atoms                                                    ! go through all atoms and evaluate grid points
      dist_point_atom = sqrt(sum((list_access(n,:) - coordinates(n_coords,:))**2))             ! initial distance
      do c = 1,3                                                                               ! PBCs in all direction. Here for cell_a (-1,0,+1)
        do d = 1,3                                                                             ! here for cell_b
          do e = 1,3                                                                           ! here for cell_c. Taking all surrounding unit cells into account
            new_point_atom = sqrt(sum((list_access(n,:) - coordinates(n_coords,:) + &
                                (c-2)*cell_a(:) + (d-2)*cell_b(:) + (e-2)*cell_c(:))**2))      ! evaluate new distance due to PBC
            if (new_point_atom < dist_point_atom) then                                         ! if distance is smaller -> use this one !
              dist_point_atom = new_point_atom
            end if
          end do
        end do
      end do

      loop67: do a = 1, no_elements
        if (elements(n_coords) == pse(a)) then
          if (dist_point_atom < vdW_radii(a) + probe_r*factor) then                             ! if grid point is outside an atom + the probe radius, but close by -> get this point
            sub_division(n_coords)%sub_grid_points = sub_division(n_coords)%sub_grid_points + 1 ! to check accessibility later on
            exit loop67
          end if                                                                                ! take only points in between vdW+probe_r AND vdw+probe_r*factor
        end if
      end do loop67
    end do loop13
  end do

  !
  ! Allocate the respective number of grid points for each sub grid
  !
  do a = 1, number_of_atoms
    allocate(sub_division(a)%sub_grids(sub_division(a)%sub_grid_points,3))
  end do
  !
  ! Reset number to count once more in the next loop
  !
  sub_division(:)%sub_grid_points = 0
  !
  ! Write sub-grid points per atom
  !
  do n = 1, n_access                                                                           ! go through all accessible points
   loop15: do n_coords = 1,number_of_atoms                                                     ! go through all atoms and evaluate grid points
      dist_point_atom = sqrt(sum((list_access(n,:) - coordinates(n_coords,:))**2))             ! initial distance
      do c = 1,3                                                                               ! PBCs in all direction. Here for cell_a (-1,0,+1)
        do d = 1,3                                                                             ! here for cell_b
          do e = 1,3                                                                           ! here for cell_c. Taking all surrounding unit cells into account
            new_point_atom = sqrt(sum((list_access(n,:) - coordinates(n_coords,:) + &
                                (c-2)*cell_a(:) + (d-2)*cell_b(:) + (e-2)*cell_c(:))**2))      ! evaluate new distance due to PBC
            if (new_point_atom < dist_point_atom) then                                         ! if distance is smaller -> use this one !
              dist_point_atom = new_point_atom
            end if
          end do
        end do
      end do

      loop68: do a = 1, no_elements
        if (elements(n_coords) == pse(a)) then
          if (dist_point_atom < vdW_radii(a) + probe_r*factor) then                             ! if grid point is outside an atom + the probe radius, but close by -> get this point
            sub_division(n_coords)%sub_grid_points = sub_division(n_coords)%sub_grid_points + 1
            sub_division(n_coords)%sub_grids(sub_division(n_coords)%sub_grid_points,:) = list_access(n,:)
            exit loop68
          end if                                                                                ! take only points in between vdW+probe_r AND vdw+probe_r*factor
        end if
      end do loop68
    end do loop15
  end do
  !
  ! Write to screen
  !
  write(6,*) ' '
  write(6,*) 'N_acc after 1st loop                ',n_access
  write(6,*) '  Total number of unoccupied points ',n_noOccu
  write(6,*) '  Sum of sub-grid points            ',sum(sub_division(:)%sub_grid_points)
  !
  ! Write to file
  !
  write(19,*) ' '
  write(19,*) 'N_acc after 1st loop                ',n_access
  write(19,*) '  Total number of unoccupied points ',n_noOccu
  write(19,*) '  Sum of sub-grid points            ',sum(sub_division(:)%sub_grid_points)

! HERE -> restructure like
! 1. loop: atoms
! 2. loop: corresponding grid points per atoms
! 3. loop: atoms => make sure that everything is sampled correctly!! evaluate grid points with repsect to all atoms !

! IDEA:
! - loop over noOccu points
! - check distance of the noOccu point to any atoms (thus, only loop over atoms here!)
! - if distance is within a certain distance (some cutoff radius is needed here)
! -> Only evaluate the grid points corresponding to these atoms, and NOT all n_check_acc points !!!!!
! This should be much more efficient, and much faster (if done correctly)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! 2nd looping through the grid points !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Identify last accessible points (the ones which are inside the probe radius sphere) !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do a = 1, n_noOccu

    loop1: do n_coords = 1,number_of_atoms                                                   ! go through all atoms and evaluate grid points
      dist_point_atom = sqrt(sum((list_noOccu(a,:) - coordinates(n_coords,:))**2))           ! initial distance between grid point and atom
      do c = 1,3                                                                             ! PBCs in all direction. Here for cell_a (-1,0,+1)
        do d = 1,3                                                                           ! here for cell_b
          do e = 1,3                                                                         ! here for cell_c. Taking all surrounding unit cells into account
            new_point_atom = sqrt(sum((list_noOccu(a,:) - coordinates(n_coords,:) + &
                                (c-2)*cell_a(:) + (d-2)*cell_b(:) + (e-2)*cell_c(:))**2))    ! evaluate new distance due to PBC
            if (new_point_atom < dist_point_atom) then                                       ! if distance is smaller -> use this one !
              dist_point_atom = new_point_atom
            end if
          end do
        end do
      end do

      ! Get the points which are initially unoccupied, but are accessible
      loop69: do n = 1, no_elements
        if (elements(n_coords) == pse(n)) then
          if (dist_point_atom < vdW_radii(n) + probe_r) then                                ! If grid point is outside an atom (see last if statement), but within a length of the probe radius
            do f = 1, sub_division(n_coords)%sub_grid_points                                ! Go through the 'check accessibility' points (which have been determined before)
              dist_point_point = sqrt(sum((sub_division(n_coords)%sub_grids(f,:) - &
                                           list_noOccu(a,:))**2))                           ! determine the distance to any accessible point
              if (dist_point_point < probe_r) then                                          ! if the distance to any accessible point is smaller than the probe radius -> NEW ACCESSIBLE POINT
                n_access = n_access + 1                                                     ! increase accessible counter
                exit loop1                                                                  ! Stop looping once this is confirmed (i.e. stop looping over the atoms)
              end if
            end do
            exit loop69
          end if
        end if
      end do loop69
    end do loop1                                                                             ! end do atoms
  end do                                                                                     ! end do full grid

  write(6,*) 'N_acc after 2nd loop                ', n_access
  write(19,*) 'N_acc after 2nd loop                ', n_access
  !
  ! Write to screen
  !
  write(6,*) ' '
  write(6,*) '#######################################################################################################'
  write(6,*) '###############################################  OUTPUT ###############################################'
  write(6,*) '#######################################################################################################'
  write(6,*) 'Structure : ',name_struct
  write(6,fmt='(1X a,7X I4,1X a,1X I4,1X a,1x I4)') 'Grid which was used:         ',grid_a,'x',grid_b,'x',grid_c                              ! 1X -> 1 space, 7X -> 7 spaces
  write(6,fmt='(1X,a,I15)') 'Total number of grid points: ',grid_a*grid_b*grid_c
  write(6,fmt='(1X,a,f15.3,a,3(f10.3,2X,a),a)') 'Grid point density:         ',grid_a*grid_b*grid_c/V_total, &
              ' grid points per A^3, with ',grid_per_A_x,' x ',grid_per_A_y,' x ',grid_per_A_z,' grid points per A'                           ! (grid_a*grid_b*grid_c/V_total)**(1./3.)
  write(6,fmt='(1X a,I15)') 'Points OCCUPIED:             ',n_occ
  write(6,fmt='(1X a,I15)') 'Points NOT OCCUPIED (void):  ',grid_a*grid_b*grid_c - n_occ
  !
  ! Pore window evaluation. If smallest pore window is smaller than the probe radius -> unoccupied, but inaccessible!!
  !
  if (file_exist) then  ! make sure the pore windows were evaluated
    if (probe_r > minval(pore_windows)) then
      write(6,fmt='(1X a,I15)') 'Points IN-ACCESSIBLE:        ',n_access
    else
      write(6,fmt='(1X a,I15)') 'Points ACCESSIBLE:           ',n_access
    end if
  else
    write(6,fmt='(1X a,I15)') 'Points ACCESSIBLE:           ',n_access
  end if

  write(6,fmt='(1X a,7X f7.3,1X a)') 'Probe radius:                ',probe_r,'A'
  if (file_exist) then  ! make sure the pore windows were evaluated
    write(6,fmt='(1X a,7X f7.3,1X a)') 'Smallest pore window:        ',minval(pore_windows),'A'
  end if
  write(6,*) ' '
 
  write(6,777) 'Porosity (void):          ',(real(grid_a*grid_b*grid_c) - real(n_occ))/(real(grid_a*grid_b*grid_c))*100D0,'%'
  !
  ! Pore window evaluation. If smallest pore window is smaller than the probe radius -> unoccupied, but inaccessible!!
  !
  if (file_exist) then  ! make sure the pore windows were evaluated
    if (probe_r > minval(pore_windows)) then
      write(6,777) 'Porosity (in-accessible): ',real(n_access)/(real(grid_a*grid_b*grid_c))*100D0,'%'
    else
      write(6,777) 'Porosity (accessible):    ',real(n_access)/(real(grid_a*grid_b*grid_c))*100D0,'%'
    end if
  else
    write(6,777) 'Porosity (accessible):    ',real(n_access)/(real(grid_a*grid_b*grid_c))*100D0,'%'
  end if
 
  V_void       = (real(grid_a*grid_b*grid_c) - real(n_occ))/(real(grid_a*grid_b*grid_c))*V_total
  V_accessible = real(n_access)/(real(grid_a*grid_b*grid_c))*V_total 
 
  write(6,777) 'Volume (void):            ',V_void,'A^3'
  !
  ! Pore window evaluation. If smallest pore window is smaller than the probe radius -> unoccupied, but inaccessible!!
  !
  if (file_exist) then  ! make sure the pore windows were evaluated
    if (probe_r > minval(pore_windows)) then
      write(6,777) 'Volume (in-accessible):   ',V_accessible,'A^3'
    else
      write(6,777) 'Volume (accessible):      ',V_accessible,'A^3'
    end if
  else
    write(6,777) 'Volume (accessible):      ',V_accessible,'A^3'
  end if
  write(6,*) ' '
 
  write(6,fmt='(1X a,f10.3,a)') 'Unit cell volume (V_total):                   ',V_total,' A^3'
  write(6,fmt='(1X a,f10.3,a)') 'Mass of unit cell (m_total):                  ',m_total*u,' 10**-27 kg'
  write(6,fmt='(1X a,f10.3,a,f10.3,a)') 'Density of the structure (m_total/V_total):   ',m_total*u/V_total*10**3,' kg/m^3 = ',&
                                                                                         m_total*u/V_total,' g/cm^3'
  write(6,fmt='(1X a,f10.3,a)') 'Pore volume density (V_void/m_total):         ',V_void/(m_total*u)*10D0**(0),' cm^3/g'
  write(6,fmt='(1X a,f10.3,a)') 'Pore volume density (V_acc/m_total):          ',V_accessible/(m_total*u)*10D0**(0),' cm^3/g'
 
  call cpu_time(finish)
  write(6,*) ' '
  write(6,fmt='(A,2X,F12.3,1X,A)') 'Total CPU time: ',finish-start,'s'
  ! 
  ! Write to file
  ! 
  write(19,*) ' '
  write(19,*) '#######################################################################################################'
  write(19,*) '###############################################  OUTPUT ###############################################'
  write(19,*) '#######################################################################################################'
  write(19,*) 'Structure : ',name_struct
  write(19,fmt='(1X a,7X I4,1X a,1X I4,1X a,1x I4)') 'Grid which was used:         ',grid_a,'x',grid_b,'x',grid_c                              ! 1X -> 1 space, 7X -> 7 spaces
  write(19,fmt='(1X,a,I15)') 'Total number of grid points: ',grid_a*grid_b*grid_c
  write(19,fmt='(1X,a,f15.3,a,3(f10.3,2X,a),a)') 'Grid point density:         ',grid_a*grid_b*grid_c/V_total, &
              ' grid points per A^3, with ',grid_per_A_x,' x ',grid_per_A_y,' x ',grid_per_A_z,' grid points per A'                           ! (grid_a*grid_b*grid_c/V_total)**(1./3.)
  write(19,fmt='(1X a,I15)') 'Points OCCUPIED:             ',n_occ
  write(19,fmt='(1X a,I15)') 'Points NOT OCCUPIED (void):  ',grid_a*grid_b*grid_c - n_occ
  !
  ! Pore window evaluation. If smallest pore window is smaller than the probe radius -> unoccupied, but inaccessible!!
  !
  if (file_exist) then  ! make sure the pore windows were evaluated
    if (probe_r > minval(pore_windows)) then
      write(19,fmt='(1X a,I15)') 'Points IN-ACCESSIBLE:        ',n_access
    else
      write(19,fmt='(1X a,I15)') 'Points ACCESSIBLE:           ',n_access
    end if
  else
    write(19,fmt='(1X a,I15)') 'Points ACCESSIBLE:           ',n_access
  end if

  write(19,fmt='(1X a,7X f7.3,1X a)') 'Probe radius:                ',probe_r,'A'
  if (file_exist) then  ! make sure the pore windows were evaluated
    write(19,fmt='(1X a,7X f7.3,1X a)') 'Smallest pore window:        ',minval(pore_windows),'A'
  end if
  write(19,*) ' '

  write(19,777) 'Porosity (void):          ',(real(grid_a*grid_b*grid_c) - real(n_occ))/(real(grid_a*grid_b*grid_c))*100D0,'%'
  !
  ! Pore window evaluation. If smallest pore window is smaller than the probe radius -> unoccupied, but inaccessible!!
  !
  if (file_exist) then  ! make sure the pore windows were evaluated
    if (probe_r > minval(pore_windows)) then
      write(19,777) 'Porosity (in-accessible): ',real(n_access)/(real(grid_a*grid_b*grid_c))*100D0,'%'
    else
      write(19,777) 'Porosity (accessible):    ',real(n_access)/(real(grid_a*grid_b*grid_c))*100D0,'%'
    end if
  else
    write(19,777) 'Porosity (accessible):    ',real(n_access)/(real(grid_a*grid_b*grid_c))*100D0,'%'
  end if

  V_void       = (real(grid_a*grid_b*grid_c) - real(n_occ))/(real(grid_a*grid_b*grid_c))*V_total
  V_accessible = real(n_access)/(real(grid_a*grid_b*grid_c))*V_total

  write(19,777) 'Volume (void):            ',V_void,'A^3'
  !
  ! Pore window evaluation. If smallest pore window is smaller than the probe radius -> unoccupied, but inaccessible!!
  !
  if (file_exist) then  ! make sure the pore windows were evaluated
    if (probe_r > minval(pore_windows)) then
      write(19,777) 'Volume (in-accessible):   ',V_accessible,'A^3'
    else
      write(19,777) 'Volume (accessible):      ',V_accessible,'A^3'
    end if
  else
    write(19,777) 'Volume (accessible):      ',V_accessible,'A^3'
  end if
  write(19,*) ' '

  write(19,fmt='(1X a,f10.3,a)') 'Unit cell volume (V_total):                   ',V_total,' A^3'
  write(19,fmt='(1X a,f10.3,a)') 'Mass of unit cell (m_total):                  ',m_total*u,' 10**-27 kg'
  write(19,fmt='(1X a,f10.3,a,f10.3,a)') 'Density of the structure (m_total/V_total):   ',m_total*u/V_total*10D0**3,' kg/m^3 = ',&
                                                                                         m_total*u/V_total,' g/cm^3'
  write(19,fmt='(1X a,f10.3,a)') 'Pore volume density (V_void/m_total):         ',V_void/(m_total*u)*10D0**(0),' cm^3/g'
  write(19,fmt='(1X a,f10.3,a)') 'Pore volume density (V_acc/m_total):          ',V_accessible/(m_total*u)*10D0**(0),' cm^3/g'

  write(19,*) ' '
  write(19,fmt='(A,2X,F12.3,1X,A)') 'Total CPU time: ',finish-start,'s'
  close(19)

  ! Define return values
  poro_void  = (real(grid_a*grid_b*grid_c) - real(n_occ))/(real(grid_a*grid_b*grid_c))*100.0D0 
  poro_acc   = real(n_access)/(real(grid_a*grid_b*grid_c))*100.0D0
  density    = m_total*u/V_total*10.0D0**3
  poreV_void = V_void/(m_total*u)*10D0**(0)
  poreV_acc  = V_accessible/(m_total*u)*10D0**(0)
  return

  deallocate(grid_points)
!  deallocate(list_occupi)
  deallocate(list_access)
  deallocate(list_check_acc)
  deallocate(list_noOccu)
  deallocate(sub_division)

777 format(1X,A,F20.5,1X,A)

  deallocate(elements)
  deallocate(coordinates)
  deallocate(pse)
  deallocate(vdW_radii)

end subroutine do_GPA



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine evaluating the occupied volume and the masses of the atoms !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eval_vol_mass(element,vocc,m)           ! element as input, V_occ and mass as output
  character(2), intent(in)      :: element         ! element symbol
  real(8), intent(inout)        :: vocc            ! occupied volume of the atoms (according to their vdW radii)
  real(8), intent(inout)        :: m               ! mass of all atoms
  real(8), parameter            :: pi = 4.0D0*atan(1.0D0)  ! define pi

  character(2)                  :: pse(86)         ! Elements
  real(8)                       :: vdW_radii(86)   ! vdW radii
  real(8)                       :: mass(86)        ! atomic mass
  integer                       :: a               ! loop variables

  pse = (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',&
           'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',&
           'Co', 'Ni', 'Cu', 'Zn', 'Zr', &
           'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', &
           'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', &
           'Rb', 'Sr', 'Y ',       'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', &
           'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', &
           'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', & 
           'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', &
           'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn' /)
  vdW_radii = (/ 1.20D0, 1.40D0, 1.82D0, 1.53D0, 1.92D0, 1.70D0, 1.55D0, 1.52D0, 1.47D0, 1.54D0,&
                 2.27D0, 1.73D0, 1.84D0, 2.10D0, 1.80D0, 1.80D0, 1.75D0, 1.88D0, 2.75D0, 2.31D0,&
                 1.92D0, 1.63D0, 1.40D0, 1.39D0, 2.36D0, &
                 2.30D0, 2.15D0, 2.05D0, 2.05D0, 2.05D0, 2.05D0, &
                 2.10D0, 2.10D0, 2.05D0, 1.90D0, 1.90D0, 2.02D0, &
                 2.90D0, 2.55D0, 2.40D0,         2.15D0, 2.10D0, 2.05D0, 2.05D0, 2.00D0, 2.05D0, 2.10D0, 2.20D0, &
                 2.20D0, 2.25D0, 2.20D0, 2.10D0, 2.10D0, 2.16D0, &
                 3.00D0, 2.70D0, 2.50D0, 2.25D0, 2.20D0, 2.10D0, 2.05D0, 2.00D0, 2.00D0, 2.05D0, 2.10D0, 2.05D0, &
                 2.48D0,2.47D0,2.45D0,2.43D0,2.42D0,2.40D0,2.38D0,2.37D0,2.35D0,2.33D0,2.32D0,2.30D0,2.28D0,2.27D0, &
                 2.20D0, 2.30D0, 2.30D0, 2.00D0, 2.00D0, 2.00D0 /)
  mass = (/  1.0079D0,   4.003D0,   6.941D0,   9.012D0,  10.811D0,  12.011D0,  14.007D0,  15.999D0,  18.998D0,  20.180D0,&
             22.990D0,  24.305D0,  26.982D0,  28.086D0,  30.974D0,  32.066D0,  35.453D0,  39.948D0,  39.099D0,  40.078D0,&
             58.933D0,  58.693D0,  63.546D0,  65.390D0,  91.224D0, &
             44.956D0,  47.867D0,  50.942D0,  51.996D0,  54.938D0,  55.845D0, &
             69.723D0,  72.610D0,  74.922D0,  78.960D0,  79.904D0,  83.800D0, &
             85.468D0,  87.620D0,  88.906D0,  92.906D0,  95.940D0,  98.000D0, 101.070D0, 102.906D0, 106.42D0, 107.868D0, 112.412D0,&
            114.818D0, 118.711D0, 121.760D0, 127.600D0, 126.904D0, 131.290D0, &
            132.905D0, 137.328D0, 138.906D0, 178.490D0, 180.948D0, 183.840D0, &
            186.207D0, 190.230D0, 192.217D0, 195.078D0, 196.967D0, 200.590D0, &
            140.116D0, 140.908D0, 144.240D0, 145.000D0, 150.360D0, 151.964D0, 157.250D0, &
            158.925D0, 162.500D0, 164.930D0, 167.260D0, 168.934D0, 173.040D0, 174.967D0, &
            204.383D0, 207.200D0, 208.980D0, 209.000D0, 210.000D0, 222.000D0 /)
!  pse = (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',&
!           'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',&
!           'Co', 'Ni', 'Cu', 'Zn', 'Zr' /)
!  vdW_radii = (/ 1.20D0, 1.40D0, 1.82D0, 1.53D0, 1.92D0, 1.70D0, 1.55D0, 1.52D0, 1.47D0, 1.54D0,&
!                 2.27D0, 1.73D0, 1.84D0, 2.10D0, 1.80D0, 1.80D0, 1.75D0, 1.88D0, 2.75D0, 2.31D0,&
!                 1.92D0, 1.63D0, 1.40D0, 1.39D0, 2.36D0 /)
!  mass = (/ 1.0079D0,  4.003D0,  6.941D0,  9.012D0, 10.811D0, 12.011D0, 14.007D0, 15.999D0, 18.998D0, 20.180D0,&
!            22.990D0, 24.305D0, 26.982D0, 28.086D0, 30.974D0, 32.066D0, 35.453D0, 39.948D0, 39.099D0, 40.078D0,&
!            58.933D0, 58.693D0, 63.546D0, 65.390D0, 91.224D0 /)
  loop99: do a = 1, 86
    if (element == pse(a)) then
      vocc = vocc + 4.0D0/3.0D0*pi*vdW_radii(a)**(3.0D0)
      m    = m    + mass(a)
      exit loop99
    end if
  end do loop99

  return
end subroutine eval_vol_mass
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to evaluate the overlap volume between two spheres !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine overlap(r_1, r_2, d_12, V_over)             ! radius of sphere 1, radius of sphere 2, distance between them, overlap volume (to be evaluated)
  real(8), intent(in)  :: r_1, r_2, d_12               ! input for the calculation of the overlap of two spheres
  real(8), intent(out) :: V_over                       ! output
  real(8), parameter   :: pi = 4.0D0*atan(1.0D0)       ! define pi
  V_over = (pi*(d_12**4D0-6D0*d_12**2D0*(r_1**2D0+r_2**2D0)+8D0*d_12*(r_1**3D0+r_2**3D0)-3D0*(r_1**2D0-r_2**2D0)**2D0))/(12D0*d_12) 
                                                       ! calculate the overlap of two spheres with radii r_1 and r_2 at distance d_12 analytically
  return                                               ! return value of overlap
end subroutine overlap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to choose the vdW radii according to the elements. Use the first subroutine to calculate the overlap volume. Return that to the program !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eval_overlap(element_a, element_b, dist_ab, sub_over)   ! evaluate the overlap according to the elements and their distance
  character(2), intent(in) :: element_a, element_b                 ! elements involved
  real(8), intent(in)      :: dist_ab                              ! distance between the elements
  real(8)                  :: r_vdw1, r_vdw2                       ! vdW radii of the atoms
  real(8), intent(out)     :: sub_over                             ! overlap volume (evaluated within the subroutine 'overlap')
 
  ! Some arrays for easier handling of vdW and covalent radii
  character(2)           :: pse(86)                               ! Elements
  real(8)                :: vdW_radii(86)                         ! vdW radii
  real(8)                :: cov_radii(86)                         ! cov radii
  integer                :: a, b                                  ! loop variables

  pse = (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',&
           'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',&
           'Co', 'Ni', 'Cu', 'Zn', 'Zr', &
           'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', &
           'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', &
           'Rb', 'Sr', 'Y ',       'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', &
           'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', &
           'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', & 
           'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', &
           'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn' /)
  vdW_radii = (/ 1.20D0, 1.40D0, 1.82D0, 1.53D0, 1.92D0, 1.70D0, 1.55D0, 1.52D0, 1.47D0, 1.54D0,&
                 2.27D0, 1.73D0, 1.84D0, 2.10D0, 1.80D0, 1.80D0, 1.75D0, 1.88D0, 2.75D0, 2.31D0,&
                 1.92D0, 1.63D0, 1.40D0, 1.39D0, 2.36D0, &
                 2.30D0, 2.15D0, 2.05D0, 2.05D0, 2.05D0, 2.05D0, &
                 2.10D0, 2.10D0, 2.05D0, 1.90D0, 1.90D0, 2.02D0, &
                 2.90D0, 2.55D0, 2.40D0,         2.15D0, 2.10D0, 2.05D0, 2.05D0, 2.00D0, 2.05D0, 2.10D0, 2.20D0, &
                 2.20D0, 2.25D0, 2.20D0, 2.10D0, 2.10D0, 2.16D0, &
                 3.00D0, 2.70D0, 2.50D0, 2.25D0, 2.20D0, 2.10D0, 2.05D0, 2.00D0, 2.00D0, 2.05D0, 2.10D0, 2.05D0, &
                 2.48D0,2.47D0,2.45D0,2.43D0,2.42D0,2.40D0,2.38D0,2.37D0,2.35D0,2.33D0,2.32D0,2.30D0,2.28D0,2.27D0, &
                 2.20D0, 2.30D0, 2.30D0, 2.00D0, 2.00D0, 2.00D0 /)
  cov_radii = (/ 0.33D0, 0.28D0, 1.28D0, 0.96D0, 0.84D0, 0.76D0, 0.71D0, 0.66D0, 0.57D0, 0.58D0,&
                 1.67D0, 1.42D0, 1.21D0, 1.11D0, 1.07D0, 1.05D0, 1.02D0, 1.06D0, 2.03D0, 1.76D0,&
                 1.26D0, 1.24D0, 1.30D0, 1.33D0, 1.48D0, &
                 1.44D0, 1.47D0, 1.33D0, 1.35D0, 1.35D0, 1.34D0, &
                 1.22D0, 1.17D0, 1.21D0, 1.22D0, 1.21D0, 1.91D0, &
                 1.47D0, 1.12D0, 1.78D0,         1.48D0, 1.47D0, 1.35D0, 1.40D0, 1.45D0, 1.50D0, 1.59D0, 1.69D0, &
                 1.63D0, 1.46D0, 1.46D0, 1.47D0, 1.40D0, 1.98D0, &
                 1.67D0, 1.34D0, 1.87D0, 1.57D0, 1.43D0, 1.37D0, 1.35D0, 1.37D0, 1.32D0, 1.50D0, 1.50D0, 1.70D0, &
                 1.83D0, 1.82D0, 1.81D0, 1.80D0, 1.80D0, 1.99D0, 1.79D0, 1.76D0, 1.75D0, 1.74D0, 1.73D0, 1.72D0, 1.94D0, 1.72D0, &
                 1.55D0, 1.54D0, 1.54D0, 1.68D0, 1.70D0, 2.40D0 /)

!  pse = (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',&
!           'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',&
!           'Co', 'Ni', 'Cu', 'Zn', 'Zr' /)
!  vdW_radii = (/ 1.20D0, 1.40D0, 1.82D0, 1.53D0, 1.92D0, 1.70D0, 1.55D0, 1.52D0, 1.47D0, 1.54D0,&
!                 2.27D0, 1.73D0, 1.84D0, 2.10D0, 1.80D0, 1.80D0, 1.75D0, 1.88D0, 2.75D0, 2.31D0,&
!                 1.92D0, 1.63D0, 1.40D0, 1.39D0, 2.36D0 /)
!  cov_radii = (/ 0.33D0, 0.28D0, 1.28D0, 0.96D0, 0.84D0, 0.76D0, 0.71D0, 0.66D0, 0.57D0, 0.58D0,&
!                 1.67D0, 1.42D0, 1.21D0, 1.11D0, 1.07D0, 1.05D0, 1.02D0, 1.06D0, 2.03D0, 1.76D0,&
!                 1.26D0, 1.24D0, 1.30D0, 1.33D0, 1.48D0 /)

  do a = 1, 86
    if (element_a == pse(a)) then
      do b = 1, 86
        if (element_b == pse(b)) then
          if (dist_ab < cov_radii(a)+cov_radii(b)) then                                       ! Evaluate, if distance between the atoms is smaller than the sum of the covalent radii
            r_vdw1 = vdW_radii(a)
            r_vdw2 = vdW_radii(b)
            call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)                                   ! evaluate overlap between the vdW spheres
          end if
        end if
      end do
    end if
  end do
  return
end subroutine eval_overlap

! End module definition
end module porosity








module PSD
contains
!
! Calculation pore size distribution
!
subroutine get_PSD(struct,start_points,cycles,&  ! structure, number of different starting points, number of MC steps
                count_pore,pore_sizes,pore_distribution,pore_pos_cart,pore_pos_frac) ! return values: Lists
  implicit none
  ! get PSD
  ! pore_finder
  ! January 13th, 2020: use adaptive step size for MC.
  !                     use larger a_step in the beginning of each MC cycle, reduce over time
  ! May 29th, 2020 : Restructuring to make it available as a module
  
  character(len=*)                    :: struct
  character(len=100)                  :: name_struct
  
  real(8)                             :: start, finish       ! timing
  
  integer(8)                          :: number_of_atoms
  real(8)                             :: cell_a(3)           ! array for the cell vector in the a direction. Vector.
  real(8)                             :: cell_b(3)           ! array for the cell vector in the b direction. Vector.
  real(8)                             :: cell_c(3)           ! array for the cell vector in the c direction. Vector.
  real(8)                             :: len_vec(3)          ! length of cell vectors
  real(8)                             :: len_max             ! length of largest cell vectors
  real(8), allocatable                :: coordinates(:,:)    ! array for the coordinates. Matrix.
  character(2), allocatable           :: elements(:)         ! array for the elements. Vector.
  
  integer(8)                          :: a,b,c,d,e,n             ! loop parameter
  real(8)                             :: rand1, rand2, rand3     ! random numbers to get new coordinates
  
  real(8)                             :: coords1(3), coords2(3)  ! coordinates of point before and after MC step
  real(8), allocatable                :: coords_all_cart(:,:)    ! all coordinates of point after MC, cartesian
  real(8), allocatable                :: coords_all_frac(:,:)    ! all coordinates of point after MC, fractional
  real(8)                             :: distance1, distance2    ! corresponding minimum distance to any atom
  real(8)                             :: tmp_dist, vdw           ! temporary distance, vdW radius of the atom
  real(8)                             :: distribution            ! final PSD distribution for a given radius
  
  integer,intent(in)                  :: start_points, cycles    ! number of starting points, number of MC cycles
  real(8)                             :: stepsize                ! step size for MC steps
  real(8), allocatable                :: all_distances(:)        ! store all distances (maybe need another list to separate different distances which occur more often.. PSD and stuff)
  real(8), allocatable                :: all_distances2(:)       ! store all distances, to double check
  real(8), allocatable                :: final_eval(:,:)         ! store final evaluated results. Use for sorting

  ! return values
  real(8), intent(out)                :: pore_sizes(100), pore_distribution(100) ! output values
  real(8), intent(out)                :: pore_pos_cart(100,3), pore_pos_frac(100,3) ! output values, coordinates of the pores
  integer(8), intent(out)             :: count_pore              ! count how many pores there are

  ! for random seed
  integer                             :: values(1:8), k
  integer, dimension(:), allocatable  :: seed
  !real(8)                             :: vv
  call date_and_time(values=values)
  call random_seed(size=k)
  allocate(seed(1:k))
  seed(:) = values(8)
  call random_seed(put=seed)
  ! end for random seed
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read in the xyz coordinates and the cell vectors !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!
  ! Initialization !
  !!!!!!!!!!!!!!!!!!
  ! Initial value
  name_struct = 'User-defined system'
  ! Read coordinate file
  open(unit=15,file=struct,status='old',action='read')
  ! Extract name of file name
  do a = 1, len(struct)-3
    if ((struct(a+1:a+3).eq.'xyz').and.(struct(a:a).eq.'.')) then
      do b = 1, len(struct)-3
        if (struct(b:b)=='/') then
          name_struct = struct(b+1:a-1)
        end if
      end do
    end if
  end do
  !
  ! Read in the number of atoms and cell vectors
  !
  read(unit=15,fmt='(I13.0)') number_of_atoms                               ! first entry is the number of atoms
  read(unit=15,fmt=*) cell_a(1:3), cell_b(1:3), cell_c(1:3)                 ! second entry contains the cell vectors. Read them in individually (makes it easier later on)
  !
  ! Get length of cell vectors. Use largest one as a starting point for step size in MC cycle
  !
  len_vec(1) = sqrt(sum(cell_a(:)**2))
  len_vec(2) = sqrt(sum(cell_b(:)**2))
  len_vec(3) = sqrt(sum(cell_c(:)**2))
  !
  ! Use 1/10 of the largest cell vector as a starting point for a_step
  !
  len_max = maxval(len_vec)*0.1D0


  allocate(elements(number_of_atoms))                             ! allocate (number_of_atoms) fields for elements. There is one elements each. As many elements as number_of_atoms (makes sense :))
  allocate(coordinates(number_of_atoms,3))                        ! allocate (number_of_atoms) fields for coordinates. There are 3 coordinates per entry. 
  do n = 1,number_of_atoms                                        ! go through all atoms 
    read(unit=15,fmt=*) elements(n), coordinates(n,1:3)           ! storing element and coordinates
  end do
  close(unit=15)                                                  ! close the file (no longer necessary)
  
  ! Output file
  open(unit=19,file='output_PSD',status='unknown',action='write')
  write(6,*) 'Starting points   (recommended: >= 200)   : ',start_points
  write(6,*) 'Monte-Carlo steps (recommended: >= 1000)  : ',cycles
  write(6,*) ' '
  write(19,*) 'Starting points   (recommended: >= 200)   : ',start_points
  write(19,*) 'Monte-Carlo steps (recommended: >= 1000)  : ',cycles
  write(19,*) ' '

  call cpu_time(start)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Monte-Carlo to get pore sizes !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(all_distances(start_points))
  allocate(all_distances2(start_points))
  allocate(coords_all_cart(start_points,3))
  allocate(coords_all_frac(start_points,3))
  !
  ! write stuff about the used step sizes
  !
  write(6,fmt='(A,I8.0,A,F10.5)') 'Until step ',ceiling(0.25D0*cycles),', step size a_step = ',1.0D0*len_max
  write(6,fmt='(A,I8.0,A,F10.5)') 'Until step ',ceiling(0.50D0*cycles),', step size a_step = ',0.1D0*len_max
  write(6,fmt='(A,I8.0,A,F10.5)') 'Until step ',ceiling(0.75D0*cycles),', step size a_step = ',0.01D0*len_max
  write(6,fmt='(A,I8.0,A,F10.5)') 'Until step ',ceiling(1.00D0*cycles),', step size a_step = ',0.001D0*len_max
  write(6,*) ' '

  write(19,fmt='(A,I8.0,A,F10.5)') 'Until step ',ceiling(0.25D0*cycles),', step size a_step = ',1.0D0*len_max
  write(19,fmt='(A,I8.0,A,F10.5)') 'Until step ',ceiling(0.50D0*cycles),', step size a_step = ',0.1D0*len_max
  write(19,fmt='(A,I8.0,A,F10.5)') 'Until step ',ceiling(0.75D0*cycles),', step size a_step = ',0.01D0*len_max
  write(19,fmt='(A,I8.0,A,F10.5)') 'Until step ',ceiling(1.00D0*cycles),', step size a_step = ',0.001D0*len_max
  write(19,*) ' '
  
  ! LOOP OVER START POINTS
  do a = 1, start_points
    ! make random number between 0.1 and 0.9. Points are inside the unit cell and not at a boundary
    call random_number(rand1)
    call random_number(rand2)
    call random_number(rand3)
    ! make random number between 0.1 and 0.9. Points are inside the unit cell and not at a boundary
    !!!rand1 = 0.1 + 0.8*rand1
    !!!rand2 = 0.1 + 0.8*rand2
    !!!rand3 = 0.1 + 0.8*rand3
    coords1(:) = cell_a(:)*rand1 + cell_b(:)*rand2 + cell_c(:)*rand3
  
  ! LOOP MC
    do b = 1, cycles
    !
    ! Adaptive step size. Check where we are in the MC cycle.
    ! Adjust a_step accordingly
    !
    ! KT: new January 2020. 
    ! start at a_step = len_of_largest_cell_vector/10.0
    ! Then, reduce by factor of 10 every1/4 of MC steps
    !
      !
      ! The very beginning -> use large a_step
      !
      if (b <= ceiling(0.25D0*cycles)) stepsize = 1.0D0*len_max
      !
      ! First half -> reduce step size
      !
      if ((b > ceiling(0.25D0*cycles)).and.(b <= ceiling(0.5D0*cycles))) stepsize = 0.1D0*len_max
      !
      ! Beginning of second half -> reduce step size
      !
      if ((b > ceiling(0.5D0*cycles)).and.(b <= ceiling(0.75D0*cycles))) stepsize = 0.01D0*len_max
      !
      ! Last quarter -> reduce step size
      !
      if (b > ceiling(0.75D0*cycles)) stepsize = 0.001D0*len_max
      
      !
      ! Now, do MC
      !
      call random_number(rand1)
      call random_number(rand2)
      call random_number(rand3)
      coords2(1) = coords1(1) + (2*rand1 - 1)*stepsize  
      coords2(2) = coords1(2) + (2*rand2 - 1)*stepsize  
      coords2(3) = coords1(3) + (2*rand3 - 1)*stepsize  
   
  ! Check distances
      distance1 = 100.0D0 ! initial values
      distance2 = 100.0D0
      do n = 1,number_of_atoms                                                                   ! go through all atoms
        if (elements(n) == 'H')  vdw = 1.20D0
        if (elements(n) == 'He') vdw = 1.40D0
        if (elements(n) == 'Li') vdw = 1.82D0
        if (elements(n) == 'Be') vdw = 1.53D0
        if (elements(n) == 'B')  vdw = 1.92D0
        if (elements(n) == 'C')  vdw = 1.70D0
        if (elements(n) == 'N')  vdw = 1.55D0
        if (elements(n) == 'O')  vdw = 1.52D0
        if (elements(n) == 'F')  vdw = 1.47D0
        if (elements(n) == 'Ne') vdw = 1.54D0
        if (elements(n) == 'Na') vdw = 2.27D0
        if (elements(n) == 'Mg') vdw = 1.73D0
        if (elements(n) == 'Al') vdw = 1.84D0
        if (elements(n) == 'Si') vdw = 2.10D0
        if (elements(n) == 'P')  vdw = 1.80D0
        if (elements(n) == 'S')  vdw = 1.80D0
        if (elements(n) == 'Cl') vdw = 1.75D0
        if (elements(n) == 'Ar') vdw = 1.88D0
        if (elements(n) == 'K')  vdw = 2.75D0
        if (elements(n) == 'Ca') vdw = 2.31D0
        if (elements(n) == 'Co') vdw = 1.92D0 ! Los Alamos value
        if (elements(n) == 'Ni') vdw = 1.63D0
        if (elements(n) == 'Cu') vdw = 1.40D0
        if (elements(n) == 'Zn') vdw = 1.39D0
        if (elements(n) == 'Zr') vdw = 2.36D0
        ! from https://www.re3data.org/repository/r3d100000013
        if (elements(n) == 'Sc') vdw = 2.30D0
        if (elements(n) == 'Ti') vdw = 2.15D0
        if (elements(n) == 'V')  vdw = 2.05D0
        if (elements(n) == 'Cr') vdw = 2.05D0
        if (elements(n) == 'Mn') vdw = 2.05D0
        if (elements(n) == 'Fe') vdw = 2.05D0
        if (elements(n) == 'Ga') vdw = 2.10D0
        if (elements(n) == 'Ge') vdw = 2.10D0
        if (elements(n) == 'As') vdw = 2.05D0
        if (elements(n) == 'Se') vdw = 1.90D0
        if (elements(n) == 'Br') vdw = 1.90D0
        if (elements(n) == 'Kr') vdw = 2.02D0
        if (elements(n) == 'Rb') vdw = 2.90D0
        if (elements(n) == 'Sr') vdw = 2.55D0
        if (elements(n) == 'Y')  vdw = 2.40D0
        if (elements(n) == 'Nb') vdw = 2.15D0
        if (elements(n) == 'Mo') vdw = 2.10D0
        if (elements(n) == 'Tc') vdw = 2.05D0
        if (elements(n) == 'Ru') vdw = 2.05D0
        if (elements(n) == 'Rh') vdw = 2.00D0
        if (elements(n) == 'Pd') vdw = 2.05D0
        if (elements(n) == 'Ag') vdw = 2.10D0
        if (elements(n) == 'Cd') vdw = 2.20D0
        if (elements(n) == 'In') vdw = 2.20D0
        if (elements(n) == 'Sn') vdw = 2.25D0
        if (elements(n) == 'Sb') vdw = 2.20D0
        if (elements(n) == 'Te') vdw = 2.10D0
        if (elements(n) == 'I')  vdw = 2.10D0
        if (elements(n) == 'Xe') vdw = 2.16D0
        if (elements(n) == 'Cs') vdw = 3.00D0
        if (elements(n) == 'Ba') vdw = 2.70D0
        if (elements(n) == 'La') vdw = 2.50D0
        if (elements(n) == 'Hf') vdw = 2.25D0
        if (elements(n) == 'Ta') vdw = 2.20D0
        if (elements(n) == 'W')  vdw = 2.10D0
        if (elements(n) == 'Re') vdw = 2.05D0
        if (elements(n) == 'Os') vdw = 2.00D0
        if (elements(n) == 'Ir') vdw = 2.00D0
        if (elements(n) == 'Pt') vdw = 2.05D0
        if (elements(n) == 'Au') vdw = 2.10D0
        if (elements(n) == 'Hg') vdw = 2.05D0
        if (elements(n) == 'Ce') vdw = 2.48D0
        if (elements(n) == 'Pr') vdw = 2.47D0
        if (elements(n) == 'Nd') vdw = 2.45D0
        if (elements(n) == 'Pm') vdw = 2.43D0
        if (elements(n) == 'Sm') vdw = 2.42D0
        if (elements(n) == 'Eu') vdw = 2.40D0
        if (elements(n) == 'Gd') vdw = 2.38D0
        if (elements(n) == 'Tb') vdw = 2.37D0
        if (elements(n) == 'Dy') vdw = 2.35D0
        if (elements(n) == 'Ho') vdw = 2.33D0
        if (elements(n) == 'Er') vdw = 2.32D0
        if (elements(n) == 'Tm') vdw = 2.30D0
        if (elements(n) == 'Yb') vdw = 2.28D0
        if (elements(n) == 'Lu') vdw = 2.27D0
        if (elements(n) == 'Tl') vdw = 2.20D0
        if (elements(n) == 'Pb') vdw = 2.30D0
        if (elements(n) == 'Bi') vdw = 2.30D0
        if (elements(n) == 'Po') vdw = 2.00D0
        if (elements(n) == 'At') vdw = 2.00D0
        if (elements(n) == 'Rn') vdw = 2.00D0

        do c = 1,3                                                                               ! PBCs in all direction. Here for cell_a (-1,0,+1)
          do d = 1,3                                                                             ! here for cell_b
            do e = 1,3                                                                           ! here for cell_c. Taking all surrounding unit cells into account
              tmp_dist = sqrt(sum((coords1(:)-coordinates(n,:)+(c-2)*cell_a(:)+(d-2)*cell_b(:)+(e-2)*cell_c(:))**2))-vdw      ! evaluate new distance due to PBC
              if (tmp_dist < distance1) then                                                     ! if distance is smaller -> use this one !
                distance1 = tmp_dist
              end if
              tmp_dist = sqrt(sum((coords2(:)-coordinates(n,:)+(c-2)*cell_a(:)+(d-2)*cell_b(:)+(e-2)*cell_c(:))**2))-vdw      ! evaluate new distance due to PBC
              if (tmp_dist < distance2) then                                                     ! if distance is smaller -> use this one !
                distance2 = tmp_dist
              end if
            end do
          end do
        end do
      end do
  
  ! Evaluate which distance is larger
      if (distance1 > distance2) then   ! if initial distance is larger
        coords2(:) = coords1(:)         ! reset second set of coordinates to be the first one
      else
        coords1(:) = coords2(:)         ! otherwise, keep new distances
      end if
    end do    ! end MC
  
  ! Get probe diameter
    if (distance1 > distance2) then
  !    write(6,*) "Start ",a," Final distance ",distance1*2.0   ! write diameter, not radius. Store position as well
      all_distances(a) = distance1*2.0D0
      coords_all_cart(a,:) = coords2(:)       
    else if (distance1 < distance2) then
  !    write(6,*) "Start ",a," Final distance ",distance2*2.0   ! write diameter, not radius. Store position as well
      all_distances(a) = distance2*2.0D0
      coords_all_cart(a,:) = coords1(:)
    ! If both are the exact same -> both are still 100.0D0 -> store 0.0D0
    else
      all_distances(a) = 0.0D0
      coords_all_cart(a,:) = (/ 0.0D0, 0.0D0, 0.0D0 /)
    end if
  end do      ! end starting points
  
  !
  ! Get proper cartesian coordinates (within the unit cell, not outside) and fractional coordinates for the pore centers
  ! 
  do a = 1, start_points
    call cart_to_frac(cell_a,cell_b,cell_c,coords_all_cart(a,:),coords_all_frac(a,:))
  end do
  !
  ! store all distances in a second array -> use for double-checking
  !
  all_distances2(:) = all_distances(:)
  
  ! Get distribution
  write(6,*) ' '
  write(6,*) 'Pore size distribution for ',name_struct
  write(6,*) ' '
  write(6,*) '  Pore size   Distribution [%]   coordinate (cartesian)                 coordinate (fractional)'
  
  write(19,*) ' '
  write(19,*) 'Pore size distribution for ',name_struct
  write(19,*) ' '
  write(19,*) '  Pore size   Distribution [%]   coordinate (cartesian)                 coordinate (fractional)'
  !
  ! Allocate final evaluation array
  !
  allocate(final_eval(100,8)) ! diameter, distribution, 3x cartesian coords, 3x fractional coords.
  final_eval(:,:) = 0.0D0
  !
  ! Go through all points. Write final results into array
  !
  count_pore = 0
  do a = 1, start_points
    distance1    = all_distances(a)
    distribution = 0.0D0
    do b = 1, start_points
      if ((abs(all_distances(a) - all_distances2(b)) < 0.10D0).and.(all_distances(a).ne.0.0D0)) then    ! collect data which is within this range of the value
        distribution = distribution + 1.0D0
        !
        ! If all_distance2 larger than all_distances -> wait til this distance is evaluated
        ! -> ensure that always the largest values is taken
        !
        if (all_distances2(b) > distance1) then
          distance1            = all_distances2(b)
          coords_all_cart(a,:) = coords_all_cart(b,:)
          coords_all_frac(a,:) = coords_all_frac(b,:)
        end if
        all_distances2(b) = 1000.0D0                  ! do not evaluate this point again
      end if
    end do
    distribution = distribution/start_points*100.D0  ! distribution in %
    if (distribution > 5.0D0) then      ! if less than 5 % -> do not evaluate
      count_pore = count_pore + 1
      final_eval(count_pore,1) = distance1
      final_eval(count_pore,2) = distribution
      final_eval(count_pore,3) = coords_all_cart(a,1)
      final_eval(count_pore,4) = coords_all_cart(a,2)
      final_eval(count_pore,5) = coords_all_cart(a,3)
      final_eval(count_pore,6) = coords_all_frac(a,1)
      final_eval(count_pore,7) = coords_all_frac(a,2)
      final_eval(count_pore,8) = coords_all_frac(a,3)
    end if
  end do
  !
  ! If any sizes are very close to each other -> add their values up 
  !
  !
  ! Sort the output, from smallest to largest pore.
  !
  do a = 1, count_pore
    do b = 1, count_pore
      !
      ! If pore a is smaller than pore b
      !
      if (final_eval(a,1) <= final_eval(b,1)) then
        !
        ! If entry b comes before entry a
        !
        if (b < a) then
          !
          ! Change entries
          !
          ! set a temporary space to the entries of b
          final_eval(count_pore+1,:) = final_eval(b,:)
          ! set the values of b to the values of a
          final_eval(b,:) = final_eval(a,:)
          ! Set the values of a to the ones of the temporary space
          final_eval(a,:) = final_eval(count_pore+1,:)
        end if
      end if
    end do
  end do

!!!  !
!!!  ! Go thorugh all points once. Evaluate how many pores there are
!!!  !
!!!  count_pore = 0
!!!  do a = 1, start_points
!!!    distribution = 0.0D0
!!!    do b = 1, start_points
!!!      if (abs(all_distances(a) - all_distances2(b)) < 0.10D0) then    ! collect data which is within this range of the value
!!!        distribution = distribution + 1.0D0
!!!        all_distances2(b) = 1000.0D0                  ! do not evaluate this point again
!!!      end if
!!!    end do
!!!    distribution = distribution/start_points*100.D0  ! distribution in %
!!!    if (distribution > 5.0D0) then      ! if less than 5 % -> do not evaluate
!!!      count_pore = count_pore + 1
!!!    end if
!!!  end do
!!!  !
!!!  ! Restore all distances in a second array -> use for final evaluation
!!!  !
!!!  all_distances2(:) = all_distances(:)
!!!  !
!!!  ! Allocate final evaluation array
!!!  !
!!!  allocate(final_eval(count_pore+1,8)) ! diameter, distribution, 3x cartesian coords, 3x fractional coords. Count_pore + 1 for sorting (last entry can be used as temporary storage)
!!!  !
!!!  ! Go through all point again. Write final results into array
!!!  !
!!!  count_pore = 0
!!!  do a = 1, start_points
!!!    distribution = 0.0D0
!!!    do b = 1, start_points
!!!      if (abs(all_distances(a) - all_distances2(b)) < 0.10D0) then    ! collect data which is within this range of the value
!!!        distribution = distribution + 1.0D0
!!!        !
!!!        ! If all_distance2 larger than all_distances -> store the former
!!!        ! -> ensure that always the largest values is taken
!!!        !
!!!        if (all_distances2(b) > all_distances(a)) then
!!!          all_distances(a) = all_distances2(b)
!!!          coords_all_cart(a,:) = coords_all_cart(b,:)
!!!          coords_all_frac(a,:) = coords_all_frac(b,:)
!!!        end if
!!!        all_distances2(b) = 1000.0D0                  ! do not evaluate this point again
!!!      end if
!!!    end do
!!!    distribution = distribution/start_points*100.D0  ! distribution in %
!!!    if (distribution > 5.0D0) then      ! if less than 5 % -> do not evaluate
!!!      count_pore = count_pore + 1
!!!      final_eval(count_pore,1) = all_distances(a)
!!!      final_eval(count_pore,2) = distribution
!!!      final_eval(count_pore,3) = coords_all_cart(a,1)
!!!      final_eval(count_pore,4) = coords_all_cart(a,2)
!!!      final_eval(count_pore,5) = coords_all_cart(a,3)
!!!      final_eval(count_pore,6) = coords_all_frac(a,1)
!!!      final_eval(count_pore,7) = coords_all_frac(a,2)
!!!      final_eval(count_pore,8) = coords_all_frac(a,3)
!!!    end if
!!!  end do
!!!  
!!!  !
!!!  ! Sort the output, from smallest to largest pore.
!!!  !
!!!  do a = 1, count_pore
!!!    do b = 1, count_pore
!!!      !
!!!      ! If pore a is smaller than pore b
!!!      !
!!!      if (final_eval(a,1) <= final_eval(b,1)) then
!!!        !
!!!        ! If entry b comes before entry a
!!!        !
!!!        if (b < a) then
!!!          !
!!!          ! Change entries
!!!          !
!!!          final_eval(count_pore+1,:) = final_eval(b,:)
!!!          final_eval(b,:) = final_eval(a,:)
!!!          final_eval(a,:) = final_eval(count_pore+1,:) 
!!!        end if
!!!      end if
!!!    end do
!!!  end do
  
  !
  ! Write output
  !
  do a = 1, count_pore
    write(6,fmt='(F12.6,F8.2,11X,3F12.6,2X,3F12.6)') final_eval(a,:)
    write(19,fmt='(F12.6,F8.2,11X,3F12.6,2X,3F12.6)') final_eval(a,:)
  end do
  
  call cpu_time(finish)
  write(6,*) ' '
  write(6,fmt='(A,2X,F12.3,1X,A)') 'Total CPU time: ',finish-start,'s'
  
  write(19,*) ' '
  write(19,fmt='(A,2X,F12.3,1X,A)') 'Total CPU time: ',finish-start,'s'
  close(19)
 
  ! return values
  !allocate(pore_sizes(count_pore))
  !allocate(pore_distribution(count_pore))
  pore_sizes(:)         = 0.0D0
  pore_distribution(:)  = 0.0D0
  pore_pos_cart(:,:)    = 0.0D0
  pore_pos_frac(:,:)    = 0.0D0
  do a = 1, count_pore
    pore_sizes(a)        = final_eval(a,1)
    pore_distribution(a) = final_eval(a,2)
    pore_pos_cart(a,1)   = final_eval(a,3)
    pore_pos_cart(a,2)   = final_eval(a,4)
    pore_pos_cart(a,3)   = final_eval(a,5)
    pore_pos_frac(a,1)   = final_eval(a,6)
    pore_pos_frac(a,2)   = final_eval(a,7)
    pore_pos_frac(a,3)   = final_eval(a,8)
  end do
  return 

!  deallocate(pore_sizes)
!  deallocate(pore_distribution)
  deallocate(coordinates)
  deallocate(elements)
  deallocate(all_distances)
  deallocate(all_distances2)
  deallocate(coords_all_cart)
  deallocate(coords_all_frac)

end subroutine get_PSD


subroutine cart_to_frac(vecA,vecB,vecC,pos_cart,pos_frac)
! Transform cart to fractional. Make sure that point are within the unit cell. Transform back to cartesian. Return both cart and frac
real(8), intent(in)    :: vecA(3), vecB(3), vecC(3)
real(8), intent(inout) :: pos_cart(3)
real(8), intent(inout) :: pos_frac(3)
real(8)                :: trans_matrix(3,3)
real(8)                :: determinant
real(8)                :: lenA, lenB, lenC, angleBC, angleAC, angleAB
integer                :: t,f
real(8), parameter     :: pi = 4.0D0*atan(1.0D0)

lenA    = sqrt(vecA(1)**2+vecA(2)**2+vecA(3)**2)
lenB    = sqrt(vecB(1)**2+vecB(2)**2+vecB(3)**2)
lenC    = sqrt(vecC(1)**2+vecC(2)**2+vecC(3)**2)
angleBC = ACOS((vecB(1)*vecC(1) + vecB(2)*vecC(2) + vecB(3)*vecC(3))/(lenB*lenC))*180.0D0/pi
angleAC = ACOS((vecA(1)*vecC(1) + vecA(2)*vecC(2) + vecA(3)*vecC(3))/(lenA*lenC))*180.0D0/pi
angleAB = ACOS((vecA(1)*vecB(1) + vecA(2)*vecB(2) + vecA(3)*vecB(3))/(lenA*lenB))*180.0D0/pi
! deteminant of the matrix comntaining the cell vectors
determinant = vecA(1)*vecB(2)*vecC(3)+vecB(1)*vecC(2)*vecA(3)+vecC(1)*vecA(2)*vecB(3) - &
              vecA(3)*vecB(2)*vecC(1)-vecB(3)*vecC(2)*vecA(1)-vecC(3)*vecA(2)*vecB(1)
! transformation matrix to get fractional coordinates. It is the inverse of the matrix containing the cell vectors
trans_matrix(1,1) = (vecB(2)*vecC(3)-vecB(3)*vecC(2))/determinant
trans_matrix(1,2) = (vecA(3)*vecC(2)-vecA(2)*vecC(3))/determinant
trans_matrix(1,3) = (vecA(2)*vecB(3)-vecA(3)*vecB(2))/determinant
trans_matrix(2,1) = (vecB(3)*vecC(1)-vecB(1)*vecC(3))/determinant
trans_matrix(2,2) = (vecA(1)*vecC(3)-vecA(3)*vecC(1))/determinant
trans_matrix(2,3) = (vecA(3)*vecB(1)-vecA(1)*vecB(3))/determinant
trans_matrix(3,1) = (vecB(1)*vecC(2)-vecB(2)*vecC(1))/determinant
trans_matrix(3,2) = (vecA(2)*vecC(1)-vecA(1)*vecC(2))/determinant
trans_matrix(3,3) = (vecA(1)*vecB(2)-vecA(2)*vecB(1))/determinant
! frac = cart*trans_matrix
pos_frac(1) = pos_cart(1)*trans_matrix(1,1) + pos_cart(2)*trans_matrix(2,1) + pos_cart(3)*trans_matrix(3,1)
pos_frac(2) = pos_cart(1)*trans_matrix(1,2) + pos_cart(2)*trans_matrix(2,2) + pos_cart(3)*trans_matrix(3,2)
pos_frac(3) = pos_cart(1)*trans_matrix(1,3) + pos_cart(2)*trans_matrix(2,3) + pos_cart(3)*trans_matrix(3,3)
! make sure that all fractional coordinates are within 0 and 1
do f = 1, 3
  t = 0
  do while (t < 1)
    if (pos_frac(f) > 1) then
      pos_frac(f) = pos_frac(f) - 1
    end if
    if (pos_frac(f) < 0) then
      pos_frac(f) = pos_frac(f) + 1
    end if
    if ((0 <= pos_frac(f)) .and. (pos_frac(f) <= 1)) then
      t = 1
    end if
  end do
end do
! Transfrom back to cartesian. This ensures that all FODs are inside the unit cell
pos_cart(1) = pos_frac(1)*vecA(1) + pos_frac(2)*vecB(1) + pos_frac(3)*vecC(1)
pos_cart(2) = pos_frac(1)*vecA(2) + pos_frac(2)*vecB(2) + pos_frac(3)*vecC(2)
pos_cart(3) = pos_frac(1)*vecA(3) + pos_frac(2)*vecB(3) + pos_frac(3)*vecC(3)

return
end subroutine cart_to_frac

end module PSD





! Make the fortran routine executable
program porE_all
  USE porosity
  USE PSD
  real(8) :: poro,dens,porV,V_t, V_o,V_oc
  real(8) :: poro_void,poro_acc,density,poreV_void,poreV_acc
  real(8) ::pores(100), distr(100), center_cart(100,3), center_frac(100,3)
  integer(8) :: no_pores
  call get_PSD('pypore.xyz',100,1000,no_pores,pores,distr,center_cart,center_frac)
end program porE_all

