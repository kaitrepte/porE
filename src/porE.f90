program porosity

implicit none

! porE
! Author Kai Trepte
! Version January 23, 2019
! Version August 23, 2019

character(2)                        :: struct
character(len=100)                  :: name_struct
integer                             :: eval_method

integer                             :: number_of_atoms
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
integer     :: grid_a, grid_b, grid_c, aa, bb, cc, running_n             ! number of grid points along cell vectors (a,b,c), loop variables to write grid points, running variable for the loop (assign grid_points correctly)
real(8)     :: g                                                         ! if grid size is suppossed to be determined automatically -> use real, not integer
integer     :: check_grid, n_coords                                      ! loop counter for the actual evaluation (grid) and for the coordinates (coords)
integer     :: counter_access, counter_check_acc, counter_noOccu         ! counter to evaluate whether a point is accessible, check accessible, or NOT occupied
integer     :: n_access, n_occ, n_check, n_check_acc, n_noOccu           ! counter for list assignment -> accessible, occupied, used loop variable (store accessible points), counter for accessibility check list, not occupied
integer     :: pbc_a, pbc_b, pbc_c                                       ! loop variables for the check of PBCs (grid point - atom)
real(8)     :: dist_point_atom, new_point_atom, dist_point_point         ! distance from a grid point to an atom, distance evaluated due to PBC, distance between grid points
real(8)     :: V_void, V_accessible                                      ! void and accessible volume
real(8)     :: grid_per_A_x, grid_per_A_y, grid_per_A_z                  ! grid per angstrom
real(8), allocatable, dimension(3)  :: grid_points(:,:)                  ! array for the grid_points. Matrix.
!real(8), allocatable, dimension(3)  :: list_occupi(:,:)                  ! empty list to save all the occupied points, thus the ones which are inside atoms
real(8), allocatable, dimension(3)  :: list_noOccu(:,:)                  ! empty list for all NOT occupied points
real(8), allocatable, dimension(3)  :: list_access(:,:)                  ! empty list for all accessible points. Initial evaluation
real(8), allocatable, dimension(3)  :: list_check_acc(:,:)               ! empty list for the check of accessibility. Will be smaller than list_access and thus easier/faster to evaluate
!real(8), allocatable, dimension(3)  :: list_all_access(:,:)              ! empty list for all accessible points. Final evaluation (take points inside the probe radius sphere into account)

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
! Initialization
  V_occupied  = 0.0                                                                            ! initial value for the occupied volume
  m_total     = 0.0                                                                            ! and for the total mass
  do t = 1, number_of_atoms                                                                    ! go through all atoms and evaluate the total occupied volume and the total mass
    call eval_vol_mass(elements(t),V_occupied,m_total)
  end do

! Calculate overlap volume. Include PBCs
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

else if (eval_method == 2) then                                                                ! choosing second evaluation method
  call cpu_time(start)                                                                         ! initialize time measurement

! Evaluate the total mass (for later evaluation)
  V_occupied  = 0.0                                                                            ! Dummy value. Necessary to use subroutine
  m_total     = 0.0                                                                            ! initial value for the total mass
  do t = 1,number_of_atoms                                                                     ! go through all atoms and evaluate the total mass
    call eval_vol_mass(elements(t),V_occupied,m_total)
  end do

! allocate all lists. Use maximum grid points for each list, as it is not clear how much is needed
  allocate(grid_points(grid_a*grid_b*grid_c,3))                                                ! allocate (grid_size) fields for the grid points. There are 3 coordinates per entry.
  allocate(list_access(grid_a*grid_b*grid_c,3))                                                ! allocate (grid_size) fields for the accessible list. There are 3 coordinates per entry.
  allocate(list_check_acc(grid_a*grid_b*grid_c,3))                                             ! allocate (grid_size) fields for the 'check accessibility' list. There are 3 coordinates per entry.
  allocate(list_noOccu(grid_a*grid_b*grid_c,3))                                                ! allocate (grid_size) fields for the NOT occupied list. There are 3 coordinates per entry.
! use these to assign at which point of the respective lists something shall be stored
  n_access = 0
  n_occ = 0
  n_noOccu = 0
  n_check_acc = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! 1st looping through the grid points !!!!!!!!!!!!!
! Get initial list of occupied points and accessible points !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write the array of the grid points
  running_n = 0                                                                                ! initialize running variable (for the assignment of the grid_points array)
  do aa = 1,grid_a                                                                             ! go through the grid points and write grid points according to the cell vectors
    do bb = 1,grid_b
      do cc = 1,grid_c
        running_n = running_n + 1                                                                  ! increase running variable at each step
        grid_point_x = cell_a(1)/grid_a*(aa-1) + cell_b(1)/grid_b*(bb-1) + cell_c(1)/grid_c*(cc-1) ! x coordinate of grid point. Choose e.g. aa-1 to include the origin (0,0,0)
        grid_point_y = cell_a(2)/grid_a*(aa-1) + cell_b(2)/grid_b*(bb-1) + cell_c(2)/grid_c*(cc-1) ! y coordinate
        grid_point_z = cell_a(3)/grid_a*(aa-1) + cell_b(3)/grid_b*(bb-1) + cell_c(3)/grid_c*(cc-1) ! z coordinate
        grid_points(running_n,:) = (/ grid_point_x,grid_point_y,grid_point_z /)                    ! assign the respective values to the array


! Include evaluation of Occupied or Accessible here, in the grid generation -> no need to go through the entire grid again!
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
        if (counter_noOccu == number_of_atoms .and. counter_access .ne. number_of_atoms) then      ! if the counter for the 'NOT occupied' points increased for all atoms -> add to list if it is not immediately accessible
          n_noOccu = n_noOccu + 1                                                                  ! increase assignment counter for the list
          list_noOccu(n_noOccu,:) = grid_points(running_n,:)
        end if
      end do
    end do
  end do













! This can probably be simplified
! Get all points to be double check for accessibility -> don't go through all accessible points again!
! Get a factor to determine which points shall be taken -> the denser the grid, the smaller the additional distance (epsilon)
  grid_per_A_x = real(grid_a)/sqrt(cell_a(1)**2 + cell_a(2)**2 + cell_a(3)**2)
  grid_per_A_y = real(grid_b)/sqrt(cell_b(1)**2 + cell_b(2)**2 + cell_b(3)**2)
  grid_per_A_z = real(grid_c)/sqrt(cell_c(1)**2 + cell_c(2)**2 + cell_c(3)**2)
  factor = 1.0D0 + 1.0D0/((grid_per_A_x+grid_per_A_y+grid_per_A_z)/3.0D0)                      ! divide by average grid points per A

  do n = 1, n_access                                                                           ! go through all accessible points
    loop3: do n_coords = 1,number_of_atoms                                                     ! go through all atoms and evaluate grid points
      dist_point_atom = sqrt(sum((list_access(n,:) - coordinates(n_coords,:))**2))             ! initial distance between grid point and atom
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
            n_check_acc = n_check_acc + 1
            list_check_acc(n_check_acc,:) = list_access(n,:)                                       ! add this point to the check_acc list. Stop looping
            exit loop3
          end if                                                                                ! take only points in between vdW+probe_r AND vdw+probe_r*factor
          exit loop67
        end if
      end do loop67
    end do loop3
  end do




















  write(6,*) ' '
  write(6,*) 'N_acc after 1st loop      ', n_access
  write(6,*) 'N_check_acc after 1st loop', n_check_acc
  write(19,*) ' '
  write(19,*) 'N_acc after 1st loop      ', n_access
  write(19,*) 'N_check_acc after 1st loop', n_check_acc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! 2nd looping through the grid points !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Identify last accessible points (the ones which are inside the probe radius sphere) !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do check_grid = 1, n_noOccu                                                                ! go through all grid points which are NOT occupied
    loop1: do n_coords = 1,number_of_atoms                                                   ! go through all atoms and evaluate grid points
      dist_point_atom = sqrt(sum((list_noOccu(check_grid,:) - coordinates(n_coords,:))**2))  ! initial distance between grid point and atom
      do c = 1,3                                                                             ! PBCs in all direction. Here for cell_a (-1,0,+1)
        do d = 1,3                                                                           ! here for cell_b
          do e = 1,3                                                                         ! here for cell_c. Taking all surrounding unit cells into account
            new_point_atom = sqrt(sum((list_noOccu(check_grid,:) - coordinates(n_coords,:) + &
                                (c-2)*cell_a(:) + (d-2)*cell_b(:) + (e-2)*cell_c(:))**2))    ! evaluate new distance due to PBC
            if (new_point_atom < dist_point_atom) then                                       ! if distance is smaller -> use this one !
              dist_point_atom = new_point_atom
            end if
          end do
        end do
      end do
    

      loop68: do n = 1, no_elements
        if (elements(n_coords) == pse(n)) then
          if (dist_point_atom < vdW_radii(n) + probe_r) then                                ! If grid point is outside an atom (see last if statement), but within a length of the probe radius
            do f = 1, n_check_acc                                                           ! Go through the 'check accessibility' points (which have been determined before)
              dist_point_point = sqrt(sum((list_noOccu(check_grid,:) - list_check_acc(f,:))**2))  ! determine the distance to any accessible point
              if (dist_point_point < probe_r) then                                          ! if the distance to any accessible point is smaller than the probe radius -> NEW ACCESSIBLE POINT
                n_access = n_access + 1                                                     ! increase accessible counter
                exit loop1                                                                  ! Stop looping once this is confirmed (i.e. stop looping over the atoms)
              end if
            end do
            exit loop68
          end if
        end if
      end do loop68
    end do loop1 
  end do                                                                                     ! end do full grid

  write(6,*) 'N_acc after 2nd loop      ', n_access
  write(19,*) 'N_acc after 2nd loop      ', n_access

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
  deallocate(list_access)
  deallocate(list_check_acc)
  deallocate(list_noOccu)
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


subroutine porefinder(structure)      ! subroutine porefinder(structure,all_pore_center,all_pore_radius)
  character(2)                        :: structure
  integer                             :: number_of_atoms
  real(8)                             :: cell_a(3)           ! array for the cell vector in the a direction. Vector.
  real(8)                             :: cell_b(3)           ! array for the cell vector in the b direction. Vector.
  real(8)                             :: cell_c(3)           ! array for the cell vector in the c direction. Vector.
  real(8), allocatable, dimension(3)  :: coordinates(:,:)    ! array for the coordinates. Matrix.
  character(2), allocatable           :: elements(:)         ! array for the elements. Vector.

  integer                             :: a,b,c,d,e,f,n,t         ! loop parameter
  real(8)                             :: rand1, rand2, rand3     ! random numbers to get new coordinates

  real(8)                             :: coords1(3), coords2(3)  ! coordinates of point before and after MC step
  real(8), allocatable                :: all_coords(:,:)         ! keep all coordinates for later usage
  real(8)                             :: distance1, distance2    ! corresponding minimum distance to any atom
  real(8)                             :: tmp_dist, vdw           ! temporary distance, vdW radius of the atom

  integer                             :: start_points, cycles    ! number of starting points, number of MC cycles
  real(8)                             :: stepsize                ! step size for MC steps
  real(8), allocatable                :: all_distances(:)        ! store all distances (maybe need another list to separate different distances which occur more often.. PSD and stuff)
  real(8), allocatable                :: all_distances2(:)       ! store all distances, to double check

  real(8)                             :: all_pore_center(10,3)   ! coordinates of the pore centers          ! intent(inout)
  real(8)                             :: all_pore_radius(10)     ! corresponding radii                      ! intent(inout)

  ! for random seed
  integer                             :: values(1:8), k
  integer, dimension(:), allocatable  :: seed
  real(8)                             :: vv
  call date_and_time(values=values)
  call random_seed(size=k)
  allocate(seed(1:k))
  seed(:) = values(8)
  call random_seed(put=seed)
  ! end for random seed

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read in the xyz coordinates and the cell vectors !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (structure == 'do') then                                                                      ! if the initial DUT-8(Ni) open structure is choosen
    open(unit=15,file='../structures_xyz/dut_8_open.xyz',status='old',action='read')               ! read in the xyz file
  else if (structure == 'vo') then                                                                 ! if the relaxed DUT-8(Ni) open structure is choosen 
    open(unit=15,file='../structures_xyz/dut_8_open_vcrelax.xyz',status='old',action='read')       ! read in the xyz file
  else if (structure == 'dc') then                                                                 ! if the initial DUT-8(Ni) closed structure is choosen
    open(unit=15,file='../structures_xyz/dut_8_closed.xyz',status='old',action='read')             ! read in the xyz file
  else if (structure == 'vc') then                                                                 ! if the relaxed DUT-8(Ni) closed structure is choosen
    open(unit=15,file='../structures_xyz/dut_8_closed_vcrelax.xyz',status='old',action='read')     ! read in the xyz file
  else if (structure == 'u6') then                                                                 ! if UiO-66 (primitive cell) is choosen
    open(unit=15,file='../structures_xyz/uio66.xyz',status='old',action='read')                    ! read in the xyz file
  else if (structure == 'u7') then                                                                 ! if UiO-67 (primitive cell) is choosen
    open(unit=15,file='../structures_xyz/uio67.xyz',status='old',action='read')                    ! read in the xyz file
  else if (structure == 'm5') then                                                                 ! if MOF-5 (unit cell) is choosen
    open(unit=15,file='../structures_xyz/mof5.xyz',status='old',action='read')                     ! read in the xyz file
  else if (structure == 'ir') then                                                                 ! if IRMOF-10 (unit cell) is choosen
    open(unit=15,file='../structures_xyz/irmof10.xyz',status='old',action='read')                  ! read in the xyz file
  else if (structure == 'm2') then                                                                 ! if MOF210 (primitive cell) is choosen
    open(unit=15,file='../structures_xyz/mof210.xyz',status='old',action='read')                   ! read in the xyz file
  else if (structure == 'h1') then                                                                 ! if HKUST-1 (primitive cell) is choosen
    open(unit=15,file='../structures_xyz/hkust1.xyz',status='old',action='read')                   ! read in the xyz file
  else if (structure == 'be') then                                                                 ! if benzene (arbitrary cell) is choosen
    open(unit=15,file='../structures_xyz/benzene.xyz',status='old',action='read')                  ! read in the xyz file
  else if (structure == 'b2') then                                                                 ! if benzene, experimental structure (arbitrary cell) is choosen
    open(unit=15,file='../structures_xyz/benzene_exp.xyz',status='old',action='read')              ! read in the xyz file
  else if (structure == 'bc') then                                                                 ! if benzene, only C atoms (arbitrary cell) is choosen
    open(unit=15,file='../structures_xyz/benzene_Conly.xyz',status='old',action='read')            ! read in the xyz file
  else if (structure == 'ha') then                                                                 ! if H atom (cubic cell) is choosen
    open(unit=15,file='../structures_xyz/h_atom.xyz',status='old',action='read')                   ! read in the xyz file
  end if
  ! Read in the corresponding values
  read(unit=15,fmt='(I13.0)') number_of_atoms                     ! first entry is the number of atoms
  read(unit=15,fmt=*) cell_a(1:3), cell_b(1:3), cell_c(1:3)       ! second entry contains the cell vectors. Read them in individually (makes it easier later on)

  allocate(elements(number_of_atoms))                             ! allocate (number_of_atoms) fields for elements. There is one elements each. As many elements as number_of_atoms (makes sense :))
  allocate(coordinates(number_of_atoms,3))                        ! allocate (number_of_atoms) fields for coordinates. There are 3 coordinates per entry. 
  do n = 1,number_of_atoms                                        ! go through all atoms 
    read(unit=15,fmt=*) elements(n), coordinates(n,1:3)           ! storing element and coordinates
  end do
  close(unit=15)

  !! initialize stuff
  do a = 1, 10
    all_pore_center(a,:) = (/ real(0.0,8), real(0.0,8), real(0.0,8) /)
    all_pore_radius(a)   = real(0.0,8)
  end do


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Monte-Carlo to get pore sizes !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize stuff
  start_points = 100   ! 100
  cycles       = 10000  ! 10000
  stepsize     = 0.01
  allocate(all_distances(start_points))
  allocate(all_distances2(start_points))
  allocate(all_coords(start_points,3))

  ! LOOP OVER START POINTS
  do a = 1, start_points
    call random_number(rand1)
    call random_number(rand2)
    call random_number(rand3)
    ! make random number between 0.1 and 0.9. Points are inside the unit cell and not at a boundary
    rand1 = 0.1 + 0.8*rand1
    rand2 = 0.1 + 0.8*rand2
    rand3 = 0.1 + 0.8*rand3
    coords1(:) = cell_a(:)*rand1 + cell_b(:)*rand2 + cell_c(:)*rand3

  ! LOOP MC
    do b = 1, cycles
      call random_number(rand1)
      call random_number(rand2)
      call random_number(rand3)
      coords2(1) = coords1(1) + (2*rand1 - 1)*stepsize
      coords2(2) = coords1(2) + (2*rand2 - 1)*stepsize
      coords2(3) = coords1(3) + (2*rand3 - 1)*stepsize
  ! Check distances
      distance1 = 100.0 ! initial values
      distance2 = 100.0
      do n = 1,number_of_atoms                                                                   ! go through all atoms
        if (elements(n) == 'H')  vdw = 1.20
        if (elements(n) == 'He') vdw = 1.40
        if (elements(n) == 'Li') vdw = 1.82
        if (elements(n) == 'Be') vdw = 1.53
        if (elements(n) == 'B')  vdw = 1.92
        if (elements(n) == 'C')  vdw = 1.70
        if (elements(n) == 'N')  vdw = 1.55
        if (elements(n) == 'O')  vdw = 1.52
        if (elements(n) == 'F')  vdw = 1.47
        if (elements(n) == 'Ne') vdw = 1.54
        if (elements(n) == 'Na') vdw = 2.27
        if (elements(n) == 'Mg') vdw = 1.73
        if (elements(n) == 'Al') vdw = 1.84
        if (elements(n) == 'Si') vdw = 2.10
        if (elements(n) == 'P')  vdw = 1.80
        if (elements(n) == 'S')  vdw = 1.80
        if (elements(n) == 'Cl') vdw = 1.75
        if (elements(n) == 'Ar') vdw = 1.88
        if (elements(n) == 'K')  vdw = 2.75
        if (elements(n) == 'Ca') vdw = 2.31
        if (elements(n) == 'Co') vdw = 1.92 ! Los Alamos value
        if (elements(n) == 'Ni') vdw = 1.63
        if (elements(n) == 'Cu') vdw = 1.40
        if (elements(n) == 'Zn') vdw = 1.39
        if (elements(n) == 'Zr') vdw = 2.36
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

    ! store all coordinates of the pore centers
    all_coords(a,:) = coords1(:)

  ! Get probe diameter radius
    if (distance1 > distance2) then
      all_distances(a) = distance1
    else
      all_distances(a) = distance2
    end if
  end do      ! end starting points

!!!!!
  ! tranfer all coordinates back into the unit cell!!!!
  do a = 1, start_points
    call frac_cart(cell_a,cell_b,cell_c,all_coords(a,:))
  end do
!!!!!

  all_distances2(:) = all_distances(:)

  ! Get distribution
  write(6,*) ' '
  write(6,*) 'Pore size distribution (diameter in angstrom, PSD in % - excluding too small contributions)'

  do a = 1, start_points
    c = 0
    do b = 1, start_points
      if (abs(all_distances(a) - all_distances2(b)) < 0.10) then    ! collect data which is within this range of the value
        c = c + 1
        all_distances2(b) = 1000.0                  ! do not evaluate this point again
      end if
    end do
    if (c <= 0.05*start_points) then                ! only pores with more than 5 % contribution will be handled
    else
      write(6,*) all_distances(a)*2.0, c   ! print diameter instead of the radius

 ! NEW
      do b = 1, 10
        f = 0
        rand1 = sum(sqrt(all_coords(a,:) - all_pore_center(b,:)))                                ! distance between new coordinate and any other that has been stored

        do c = 1,3                                                                               ! PBCs in all direction. Here for cell_a (-1,0,+1)
          do d = 1,3                                                                             ! here for cell_b
            do e = 1,3                                                                           ! here for cell_c. Taking all surrounding unit cells into account
              rand2 = sqrt(sum((all_coords(a,:)-all_pore_center(b,:)+(c-2)*cell_a(:)+(d-2)*cell_b(:)+(e-2)*cell_c(:))**2))
              if (rand2 < rand1) then                                                     ! if distance is smaller -> use this one !
                rand1 = rand2
              end if
            end do
          end do
        end do

        ! if the new coordinate is within an already existing pore -> see which one is bigger and keep that
        if ((rand1) < (all_pore_radius(b) + all_distances(a))) then
          f = f + 1
          if (all_distances(a) > all_pore_radius(b)) then                                          ! if the new radius is bigger -> overwrite old one
            all_pore_center(b,:) = all_coords(a,:)
            all_pore_radius(b) = all_distances(a)
          end if
        end if

        if (f == 0) then                                                                           ! if no other points is close to this one -> store
          if (all_pore_center(b,1) == 0.0) then                                                    ! for an entry that has not been written yet
            all_pore_center(b,:) = all_coords(a,:)
            all_pore_radius(b) = all_distances(a)
            exit
          end if
        end if
      end do
  ! END NEW
    end if
  end do

  deallocate(coordinates)
  deallocate(elements)
  deallocate(all_distances)
  deallocate(all_distances2)

  close(20)

  return
end subroutine porefinder



subroutine frac_cart(vecA,vecB,vecC,pos_cart)
real(8), intent(in)    :: vecA(3), vecB(3), vecC(3)
real(8), intent(inout) :: pos_cart(3)
real(8)                :: pos_frac(3)
real(8)                :: trans_matrix(3,3)
real(8)                :: determinant
real(8)                :: lenA, lenB, lenC, angleBC, angleAC, angleAB, vol
integer                :: t,f
real(8), parameter     :: pi = 3.14159265358979323846

lenA    = sqrt(vecA(1)**2+vecA(2)**2+vecA(3)**2)
lenB    = sqrt(vecB(1)**2+vecB(2)**2+vecB(3)**2)
lenC    = sqrt(vecC(1)**2+vecC(2)**2+vecC(3)**2)
angleBC = ACOS((vecB(1)*vecC(1) + vecB(2)*vecC(2) + vecB(3)*vecC(3))/(lenB*lenC))*180.0/pi
angleAC = ACOS((vecA(1)*vecC(1) + vecA(2)*vecC(2) + vecA(3)*vecC(3))/(lenA*lenC))*180.0/pi
angleAB = ACOS((vecA(1)*vecB(1) + vecA(2)*vecB(2) + vecA(3)*vecB(3))/(lenA*lenB))*180.0/pi
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
end subroutine frac_cart
