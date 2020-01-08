program porosity

implicit none

! POROsity Whole Analysis Tool (porowat)
! Author Kai Trepte
! Version January 23, 2019
! Version April 29th, 2019 -- implement sub-grid division
! Version August 23rd, 2019 -- add many more elements, vdW and covalent radii
! Version September 4th, 2019 -- read initial information from input file

character(2)                        :: struct
character(len=100)                  :: name_struct
integer                             :: eval_method

integer(8)                          :: number_of_atoms
real(8)                             :: cell_a(3)           ! array for the cell vector in the a direction. Vector.
real(8)                             :: cell_b(3)           ! array for the cell vector in the b direction. Vector.
real(8)                             :: cell_c(3)           ! array for the cell vector in the c direction. Vector.
real(8), allocatable, dimension(3)  :: coordinates(:,:)    ! array for the coordinates. Matrix.
character(2), allocatable           :: elements(:)         ! array for the elements. Vector.
real(8)                             :: V_total             ! total volume of the cell

real(8)                             :: start, finish       ! evaluate the time

real(8), parameter                  :: pi = 3.14159265358979323846 ! define pi
real(8), parameter                  ::  u = 1.660539               ! define atomic mass unit (in 10**-27 kg)

! Evaluation method 1
real(8)     :: sub_overlap                                               ! overlap volume as evaluated by the subroutine
real(8)     :: V_occupied, V_overlap, m_total, distance_ab, new_distance ! volume occupied by vdw spheres, total overlap volume, total mass of unit cell, distance between two atoms, distance evaluated due to PBC
real(8)     :: r_vdw1, r_vdw2                                            ! vdW radii of two species. Needed for evaluation (makes it easier)
integer     :: a,b,c,d,e,f,n,t                                           ! loop parameter

! Evaluation method 2
real(8)     :: probe_r, grid_point_x, grid_point_y, grid_point_z, factor ! probe radius, grid point coordinates for any specific grid point (x,y,z), grid density (grid points per A^3)
integer(8)  :: grid_a, grid_b, grid_c, aa, bb, cc, running_n             ! number of grid points along cell vectors (a,b,c), loop variables to write grid points, running variable for the loop (assign grid_points correctly)
real(8)     :: g                                                         ! if grid size is suppossed to be determined automatically -> use real, not integer
integer     :: check_grid, n_coords                                      ! loop counter for the actual evaluation (grid) and for the coordinates (coords)
integer     :: counter_access, counter_check_acc, counter_noOccu         ! counter to evaluate whether a point is accessible, check accessible, or NOT occupied
integer     :: n_access, n_occ, n_check, n_check_acc, n_noOccu           ! counter for list assignment -> accessible, occupied, used loop variable (store accessible points), counter for accessibility check list, not occupied
integer     :: pbc_a, pbc_b, pbc_c                                       ! loop variables for the check of PBCs (grid point - atom)
real(8)     :: dist_point_atom, new_point_atom, dist_point_point         ! distance from a grid point to an atom, distance evaluated due to PBC, distance between grid points
real(8)     :: V_void, V_accessible                                      ! void and accessible volume
real(8)     :: grid_per_A_x, grid_per_A_y, grid_per_A_z                  ! grid per angstrom
real(8), allocatable, dimension(3)  :: grid_points(:,:)                  ! array for the grid_points. Matrix.
real(8), allocatable, dimension(3)  :: list_noOccu(:,:)                  ! empty list for all NOT occupied points
real(8), allocatable, dimension(3)  :: list_access(:,:)                  ! empty list for all accessible points. Initial evaluation
real(8), allocatable, dimension(3)  :: list_check_acc(:,:)               ! empty list for the check of accessibility. Will be smaller than list_access and thus easier/faster to evaluate
!real(8), allocatable, dimension(3)  :: list_all_access(:,:)              ! empty list for all accessible points. Final evaluation (take points inside the probe radius sphere into account)

!
! IDEA: Subdivide grid into sub-grids per atom. Then, only evaluate the necessary grid points.
! Thus, exclude grid points which are very far away. 
! - Introduce new array, dimensions are number_of_atoms AND a specific number according to the specific atom (needs to be counted, need two loops over the grid points. Should still be cheap)
! - When looping, choose grid point for an atom if it is close to this grid point -> store grid point for the sub-grid for the atom 
!
! For sub-grid generation, one needs a variable array which is allocatable
type global_array
  integer              :: sub_grid_points                       ! dimensions: n_atoms
  real(8), allocatable :: sub_grids(:,:)                        ! dimensions: n_atoms, sub_grid_points(atom), 3
end type global_array
type(global_array), dimension(:), allocatable :: sub_division

! Some arrays for easier handling of vdW radii
character(2)              :: all_pse(25)                        ! all currently available elements
real(8)                   :: all_vdW_radii(25)                  ! the respective vdW radii
integer                   :: all_elements                       ! number of all different elements
character(2), allocatable :: tmp_pse(:)                         ! temporary list to evaluate the used elements
character(2), allocatable :: pse(:)                             ! used elements
real(8), allocatable      :: vdW_radii(:)                       ! used vdW radii
integer                   :: no_elements                        ! number of different elements

all_elements = 25
all_pse = (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',&
             'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',&
             'Co', 'Ni', 'Cu', 'Zn', 'Zr' /)
all_vdW_radii = (/ 1.20D0, 1.40D0, 1.82D0, 1.53D0, 1.92D0, 1.70D0, 1.55D0, 1.52D0, 1.47D0, 1.54D0,&
                   2.27D0, 1.73D0, 1.84D0, 2.10D0, 1.80D0, 1.80D0, 1.75D0, 1.88D0, 2.75D0, 2.31D0,& 
                   1.92D0, 1.63D0, 1.40D0, 1.39D0, 2.36D0 /)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Script to evaluate the porosity of a crystal structure. Taking information from xyz file and cell parameters.                                 !
!                                                                                                                                               !
! 1. evaluation: calculate total volume, occupied volume (excluding overlap) and void volume -> P = V_void/V_total                              !
!                                                                                                                                               !
! 2. evaluation: place grid inside the unit cell, count each point which is in an occupied region, compare to total points -> P = N_occ/N_tot   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!
! Initialization !
!!!!!!!!!!!!!!!!!!
open(unit=16,file='input_porosity',status='old',action='read')      ! Read the input file
read(16,*)                                                          ! ignore the first line
read(16,*) struct
if (struct == 'ud') then                                            ! if user-defined structure is choosen
  read(16,*) name_struct                                            ! read the path to the xyz file
  open(unit=15,file=name_struct,status='old',action='read')         ! read in the xyz file
  read(16,*) name_struct                                            ! read name of the structure
else
  read(16,*)                                                        ! For any other input, ignore the next two lines
  read(16,*)
  !
  ! Read in the xyz coordinates and the cell vectors
  !
  if (struct == 'do') then                                                                      ! if the initial DUT-8(Ni) open structure is choosen
    open(unit=15,file='../structures_xyz/dut_8_open.xyz',status='old',action='read')               ! read in the xyz file
    name_struct = 'DUT-8(Ni) open, exp'
  else if (struct == 'vo') then                                                                 ! if the relaxed DUT-8(Ni) open structure is choosen 
    open(unit=15,file='../structures_xyz/dut_8_open_vcrelax.xyz',status='old',action='read')       ! read in the xyz file
    name_struct = 'DUT-8(Ni) open, vc'
  else if (struct == 'dc') then                                                                 ! if the initial DUT-8(Ni) closed structure is choosen
    open(unit=15,file='../structures_xyz/dut_8_closed.xyz',status='old',action='read')             ! read in the xyz file
    name_struct = 'DUT-8(Ni) closed, exp'
  else if (struct == 'vc') then                                                                 ! if the relaxed DUT-8(Ni) closed structure is choosen
    open(unit=15,file='../structures_xyz/dut_8_closed_vcrelax.xyz',status='old',action='read')     ! read in the xyz file
    name_struct = 'DUT-8(Ni) closed, vc'
  else if (struct == 'u6') then                                                                 ! if UiO-66 (primitive cell) is choosen
    open(unit=15,file='../structures_xyz/uio66.xyz',status='old',action='read')                    ! read in the xyz file
    name_struct = 'UiO-66'
  else if (struct == 'u7') then                                                                 ! if UiO-67 (primitive cell) is choosen
    open(unit=15,file='../structures_xyz/uio67.xyz',status='old',action='read')                    ! read in the xyz file
    name_struct = 'UiO-67'
  else if (struct == 'u8') then                                                                 ! if UiO-68 (primitive cell) is choosen
    open(unit=15,file='../structures_xyz/uio68.xyz',status='old',action='read')                    ! read in the xyz file
    name_struct = 'UiO-68'
  else if (struct == 'm5') then                                                                 ! if MOF-5 (unit cell) is choosen
    open(unit=15,file='../structures_xyz/mof5.xyz',status='old',action='read')                     ! read in the xyz file
    name_struct = 'MOF-5'
  else if (struct == 'ir') then                                                                 ! if IRMOF-10 (unit cell) is choosen
    open(unit=15,file='../structures_xyz/irmof10.xyz',status='old',action='read')                  ! read in the xyz file
    name_struct = 'IRMOF10'
  else if (struct == 'm2') then                                                                 ! if MOF210 (primitive cell) is choosen
    open(unit=15,file='../structures_xyz/mof210.xyz',status='old',action='read')                   ! read in the xyz file
    name_struct = 'MOF-210'
  else if (struct == 'h1') then                                                                 ! if HKUST-1 (primitive cell) is choosen
    open(unit=15,file='../structures_xyz/hkust1.xyz',status='old',action='read')                   ! read in the xyz file
    name_struct = 'HKUST-1'
  else if (struct == 'ho') then                                                                 ! if HKUST-1 (primitive cell) is choosen
    open(unit=15,file='../structures_xyz/hkust1_with_O.xyz',status='old',action='read')             ! read in the xyz file
    name_struct = 'HKUST-1'
  else if (struct == 'c6') then                                                                 ! if C60@MOF (primitive cell) is choosen
    open(unit=15,file='../structures_xyz/c60_MOF.xyz',status='old',action='read')                   ! read in the xyz file
    name_struct = 'C60@MOF'
  else if (struct == 'be') then                                                                 ! if benzene (arbitrary cell) is choosen
    open(unit=15,file='../structures_xyz/benzene.xyz',status='old',action='read')                  ! read in the xyz file
    name_struct = 'Benzene, opt'
  else if (struct == 'b2') then                                                                 ! if benzene, experimental structure (arbitrary cell) is choosen
    open(unit=15,file='../structures_xyz/benzene_exp.xyz',status='old',action='read')              ! read in the xyz file
    name_struct = 'Benzene, exp'
  else if (struct == 'bc') then                                                                 ! if benzene, only C atoms (arbitrary cell) is choosen
    open(unit=15,file='../structures_xyz/benzene_Conly.xyz',status='old',action='read')            ! read in the xyz file
    name_struct = 'Benzene, C only'
  else if (struct == 'ha') then                                                                 ! if H atom (cubic cell) is choosen
    open(unit=15,file='../structures_xyz/h_atom.xyz',status='old',action='read')                   ! read in the xyz file
    name_struct = 'H atom'
  else
    write(6,*) '! Structure invalid !'
    stop
  end if
end if
!
! Read in the number of atoms and cell vectors
!
read(unit=15,fmt='(I13.0)') number_of_atoms                                                 ! first entry is the number of atoms
read(unit=15,fmt=*) cell_a(1:3), cell_b(1:3), cell_c(1:3)                                   ! second entry contains the cell vectors. Read them in individually (makes it easier later on)
!
! Read evaluation method
!
read(16,*)
read(16,*)
read(16,*) eval_method
if ((eval_method /= 1) .and. (eval_method /= 2)) then
  write(6,*) '! Evaluation method invalid !'
  stop
end if
!
! For GPA, read additional values
!
if (eval_method == 2) then
  read(16,*)
  read(16,*)
  read(16,*) probe_r                            ! probe radius
  read(16,*) t                                  ! define which way the grid shall be initialized
  if (t == 1) then                              ! if grid points per unit cell vector -> read number of grid points (3 integers)
    read(16,*) grid_a, grid_b, grid_c
  else if (t == 2) then
    read(16,*) g                                ! if grid point density per angstrom (1 real) -> define grid here
      grid_a = ceiling(g*sqrt(cell_a(1)**2 + cell_a(2)**2 + cell_a(3)**2))                       ! ceiling -> round to the next higher integer. Use g * len_unit_cell_vector as the number of grid points
      grid_b = ceiling(g*sqrt(cell_b(1)**2 + cell_b(2)**2 + cell_b(3)**2))
      grid_c = ceiling(g*sqrt(cell_c(1)**2 + cell_c(2)**2 + cell_c(3)**2))
  else
    write(6,*) '! Grid initialization invalid !'
    stop
  end if
end if
close(16)                                       ! close input file
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
open(unit=19,file='output_porosity',status='unknown',action='write')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the total unit cell volume using the triple product V = a . (b x c) . In A^3 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

V_total = (cell_a(1)*cell_b(2)*cell_c(3) + cell_a(2)*cell_b(3)*cell_c(1) + cell_a(3)*cell_b(1)*cell_c(2) - &
           cell_a(3)*cell_b(2)*cell_c(1) - cell_a(1)*cell_b(3)*cell_c(2) - cell_a(2)*cell_b(1)*cell_c(3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Evaluation method 1 : Calculate overlapping spheres and evaluate the total occupied volume !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (eval_method == 1) then
  call cpu_time(start)                                                                         ! starting time for the evaluation
  !
  ! Initialization
  !
  V_occupied  = 0.0                                                                            ! initial value for the occupied volume
  m_total     = 0.0                                                                            ! and for the total mass
  do t = 1, number_of_atoms                                                                    ! go through all atoms and evaluate the total occupied volume and the total mass
    call eval_vol_mass(elements(t),V_occupied,m_total)
  end do
  !
  ! Calculate overlap volume. Include PBCs
  !
  V_overlap = 0.0                                                                               ! initial value for the overlap volume
  do a = 1,number_of_atoms-1                                                                    ! go through all pairs of atoms (first loop til number_of_atoms - 1)
    do b = a+1,number_of_atoms                                                                  ! second loop from next atom til the end. Evaluate overlap. Consider periodic boundary conditions (PBCs)! No double counting
      distance_ab = sqrt(sum((coordinates(a,:) - coordinates(b,:))**2))                         ! Initial distance between two atoms. No PBCs yet.
      do c = 1,3                                                                                ! PBCs in all direction. Here for cell_a (-1,0,+1)
        do d = 1,3                                                                              ! here for cell_b
          do e = 1,3                                                                            ! here for cell_c. Taking all surrounding unit cells into account
            new_distance = sqrt(sum((coordinates(a,:)-coordinates(b,:)+(c-2)*cell_a(:)+(d-2)*cell_b(:)+(e-2)*cell_c(:))**2))
            if (new_distance < distance_ab) then
             distance_ab = new_distance
            end if
          end do
        end do
      end do 

      sub_overlap = 0.0
      call eval_overlap(elements(a), elements(b), distance_ab, sub_overlap)                     ! use subroutine to analyze the overlap volume
      V_overlap = V_overlap + sub_overlap

      if (sub_overlap > 0.0) then
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

666 format(I5.0,1X,A,I5.0,1X,A,A,F10.5,A,F10.5,A)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Evaluation method 2 : Uniform grid. Evaluate each point as to whether it is occupied, void !!
!!                       or accessible                                                        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

else if (eval_method == 2) then                                                                ! if GPA is chosen
!  write(6,fmt='(A,2X,I15,A)') 'The calculation should take roughly ',&
!                               number_of_atoms*grid_a*grid_b*grid_c/(10**6),' s'               ! only for GPA_subgrid. 

  call cpu_time(start)                                                                         ! initialize time measurement
  !
  ! Evaluate the total mass (for later evaluation)
  !
  V_occupied  = 0.0                                                                            ! Dummy value. Necessary to use subroutine
  m_total     = 0.0                                                                            ! initial value for the total mass
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
  factor = 1.0D0 + 1.0D0/((grid_per_A_x+grid_per_A_y+grid_per_A_z)/3.0D0)                      ! divide by average grid points per A

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
  write(6,fmt='(1X a,I15)') 'Points ACCESSIBLE:           ',n_access
  write(6,fmt='(1X a,7X f7.3,1X a)') 'Probe radius:                ',probe_r,'A'
  write(6,*) ' '
 
  write(6,777) 'Porosity (void):       ',(real(grid_a*grid_b*grid_c) - real(n_occ))/(real(grid_a*grid_b*grid_c))*100,'%'
  write(6,777) 'Porosity (accessible): ',real(n_access)/(real(grid_a*grid_b*grid_c))*100,'%'
 
  V_void       = (real(grid_a*grid_b*grid_c) - real(n_occ))/(real(grid_a*grid_b*grid_c))*V_total
  V_accessible = real(n_access)/(real(grid_a*grid_b*grid_c))*V_total 
 
  write(6,777) 'Volume (void):         ',V_void,'A^3'
  write(6,777) 'Volume (accessible):   ',V_accessible,'A^3'
  write(6,*) ' '
 
  write(6,fmt='(1X a,f10.3,a)') 'Unit cell volume (V_total):                   ',V_total,' A^3'
  write(6,fmt='(1X a,f10.3,a)') 'Mass of unit cell (m_total):                  ',m_total*u,' 10**-27 kg'
  write(6,fmt='(1X a,f10.3,a)') 'Density of the structure (m_total/V_total):   ',m_total*u/V_total*10**3,' kg/m^3'
  write(6,fmt='(1X a,f10.3,a)') 'Pore volume density (V_void/m_total):         ',V_void/(m_total*u)*10**(0),' cm^3/g'
  write(6,fmt='(1X a,f10.3,a)') 'Pore volume density (V_acc/m_total):          ',V_accessible/(m_total*u)*10**(0),' cm^3/g'
 
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
  write(19,fmt='(1X a,I15)') 'Points ACCESSIBLE:           ',n_access
  write(19,fmt='(1X a,7X f7.3,1X a)') 'Probe radius:                ',probe_r,'A'
  write(19,*) ' '

  write(19,777) 'Porosity (void):       ',(real(grid_a*grid_b*grid_c) - real(n_occ))/(real(grid_a*grid_b*grid_c))*100,'%'
  write(19,777) 'Porosity (accessible): ',real(n_access)/(real(grid_a*grid_b*grid_c))*100,'%'

  V_void       = (real(grid_a*grid_b*grid_c) - real(n_occ))/(real(grid_a*grid_b*grid_c))*V_total
  V_accessible = real(n_access)/(real(grid_a*grid_b*grid_c))*V_total

  write(19,777) 'Volume (void):         ',V_void,'A^3'
  write(19,777) 'Volume (accessible):   ',V_accessible,'A^3'
  write(19,*) ' '

  write(19,fmt='(1X a,f10.3,a)') 'Unit cell volume (V_total):                   ',V_total,' A^3'
  write(19,fmt='(1X a,f10.3,a)') 'Mass of unit cell (m_total):                  ',m_total*u,' 10**-27 kg'
  write(19,fmt='(1X a,f10.3,a)') 'Density of the structure (m_total/V_total):   ',m_total*u/V_total*10**3,' kg/m^3'
  write(19,fmt='(1X a,f10.3,a)') 'Pore volume density (V_void/m_total):         ',V_void/(m_total*u)*10**(0),' cm^3/g'
  write(19,fmt='(1X a,f10.3,a)') 'Pore volume density (V_acc/m_total):          ',V_accessible/(m_total*u)*10**(0),' cm^3/g'

  write(19,*) ' '
  write(19,fmt='(A,2X,F12.3,1X,A)') 'Total CPU time: ',finish-start,'s'
  close(19)

  deallocate(grid_points)
!  deallocate(list_occupi)
  deallocate(list_access)
  deallocate(list_check_acc)
  deallocate(list_noOccu)
  deallocate(sub_division)
end if           ! global end if (second method ends here)

777 format(1X,A,F20.5,1X,A)


deallocate(elements)
deallocate(coordinates)
deallocate(pse)
deallocate(vdW_radii)

end program porosity





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine evaluation the occupied volume and the masses of the atoms !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eval_vol_mass(element,vocc,m)           ! element as input, V_occ and mass as output
  character(2), intent(in)      :: element         ! element symbol
  real(8), intent(inout)        :: vocc            ! occupied volume of the atoms (according to their vdW radii)
  real(8), intent(inout)        :: m               ! mass of all atoms
  real(8), parameter            :: pi = 3.14159265358979323846  ! define pi

  character(2)           :: pse(25)                               ! Elements
  real(8)                :: vdW_radii(25)                         ! vdW radii
  real(8)                :: mass(25)                              ! atomic mass
  integer                :: a                                     ! loop variables

  pse = (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',&
           'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',&
           'Co', 'Ni', 'Cu', 'Zn', 'Zr' /)
  vdW_radii = (/ 1.20D0, 1.40D0, 1.82D0, 1.53D0, 1.92D0, 1.70D0, 1.55D0, 1.52D0, 1.47D0, 1.54D0,&
                 2.27D0, 1.73D0, 1.84D0, 2.10D0, 1.80D0, 1.80D0, 1.75D0, 1.88D0, 2.75D0, 2.31D0,&
                 1.92D0, 1.63D0, 1.40D0, 1.39D0, 2.36D0 /)
  mass = (/ 1.0079D0,  4.003D0,  6.941D0,  9.012D0, 10.811D0, 12.011D0, 14.007D0, 15.999D0, 18.998D0, 20.180D0,&
            22.990D0, 24.305D0, 26.982D0, 28.086D0, 30.974D0, 32.066D0, 35.453D0, 39.948D0, 39.099D0, 40.078D0,&
            58.933D0, 58.693D0, 63.546D0, 65.390D0, 91.224D0 /)

  loop99: do a = 1, 25
    if (element == pse(a)) then
      vocc = vocc + 4.0/3.0*pi*vdW_radii(a)**(3.0)
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
  real(8), parameter   :: pi = 3.14159265358979323846  ! define pi
  V_over = (pi*(d_12**4-6*d_12**2*(r_1**2+r_2**2)+8*d_12*(r_1**3+r_2**3)-3*(r_1**2-r_2**2)**2))/(12*d_12) ! calculate the overlap of two spheres with radii r_1 and r_2 at distance d_12 analytically
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
  character(2)           :: pse(25)                               ! Elements
  real(8)                :: vdW_radii(25)                         ! vdW radii
  real(8)                :: cov_radii(25)                         ! cov radii
  integer                :: a, b                                  ! loop variables
  
  pse = (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',&
           'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',&
           'Co', 'Ni', 'Cu', 'Zn', 'Zr' /)
  vdW_radii = (/ 1.20D0, 1.40D0, 1.82D0, 1.53D0, 1.92D0, 1.70D0, 1.55D0, 1.52D0, 1.47D0, 1.54D0,&
                 2.27D0, 1.73D0, 1.84D0, 2.10D0, 1.80D0, 1.80D0, 1.75D0, 1.88D0, 2.75D0, 2.31D0,&
                 1.92D0, 1.63D0, 1.40D0, 1.39D0, 2.36D0 /)
  cov_radii = (/ 0.33D0, 0.28D0, 1.28D0, 0.96D0, 0.84D0, 0.76D0, 0.71D0, 0.66D0, 0.57D0, 0.58D0,&
                 1.67D0, 1.42D0, 1.21D0, 1.11D0, 1.07D0, 1.05D0, 1.02D0, 1.06D0, 2.03D0, 1.76D0,&
                 1.26D0, 1.24D0, 1.30D0, 1.33D0, 1.48D0 /)

  do a = 1, 25
    if (element_a == pse(a)) then
      do b = 1, 25
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
