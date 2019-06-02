program porosity

implicit none

! POROsity Whole Analysis Tool (porowat)
! Author Kai Trepte
! Version January 23, 2019
! Version April 29th, 2019 -- try to implement sub-division of grid

character(2)                        :: struct
character(len=25)                   :: name_struct
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

! Define which structure shall be evaluated
write(6,*) '########## STRUCTURE ########'
write(6,*) 'DUT-8(Ni) open           - do'
write(6,*) 'DUT-8(Ni) open vcrelax   - vo'
write(6,*) 'DUT-8(Ni) closed         - dc'
write(6,*) 'DUT-8(Ni) closed vcrelax - vc'
write(6,*) 'UiO-66                   - u6'
write(6,*) 'UiO-67                   - u7'
write(6,*) 'MOF-5                    - m5'
write(6,*) 'IRMOF-10                 - ir'
write(6,*) 'MOF210                   - m2'
write(6,*) 'HKUST-1                  - h1'
write(6,*) 'Benzene, opt             - b1'
write(6,*) 'Benzene, exp             - b2'
write(6,*) 'Benzene, C only          - bc'
write(6,*) 'H atom                   - ha'
write(6,*) '#############################'
read(5,*) struct
! Evaluation stategy
write(6,*) ' '
write(6,*) '########## METHOD OF EVALUATION #################'
write(6,*) '(1) Overlapping sphere approach (OSA, analytical)'
write(6,*) '(2) Grid point approach         (GPA, numerical) '
write(6,*) '#################################################'
read(5,*) eval_method

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read in the xyz coordinates and the cell vectors !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Define the structure (cell vectors are the second line of the given xyz file)
! Define the structure (cell vectors are the second line of the given xyz file)
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
end if
! Read in the corresponding values
read(unit=15,fmt='(I13.0)') number_of_atoms                                                   ! first entry is the number of atoms
read(unit=15,fmt=*) cell_a(1:3), cell_b(1:3), cell_c(1:3)                                     ! second entry contains the cell vectors. Read them in individually (makes it easier later on)

allocate(elements(number_of_atoms))                                                           ! allocate (number_of_atoms) fields for elements. There is one elements each. As many elements as number_of_atoms (makes sense :))
allocate(coordinates(number_of_atoms,3))                                                      ! allocate (number_of_atoms) fields for coordinates. There are 3 coordinates per entry. 
do n = 1,number_of_atoms                                                                      ! go through all atoms 
  read(unit=15,fmt=*) elements(n), coordinates(n,1:3)                                         ! storing element and coordinates
end do
close(unit=15)                                                                                ! close the file (no longer necessary)

! Output file
open(unit=19,file='output',status='unknown',action='write')

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
       write(6,*) a,elements(a),b,elements(b),' d =', distance_ab,' V_overlap = ', sub_overlap,' A^3'
       write(19,*) a,elements(a),b,elements(b),' d =', distance_ab,' V_overlap = ', sub_overlap,' A^3'
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
  write(6,fmt='(1X a,f10.3,1X a)') 'Inverse pore volume density is (V_void/m_total) : ' &
             ,(V_total - (V_occupied - V_overlap))/(m_total*u)*10**(0),'cm^3/g'
  write(6,fmt='(1X a,f10.3,1X a)') 'The mass of unit cell is                        : ',m_total*u,'10**-27 kg'
  call cpu_time(finish)
  write(6,*) 'Total CPU time: ',finish-start,'s'

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
  write(19,fmt='(1X a,f10.3,1X a)') 'Inverse pore volume density is (V_void/m_total) : ' &
             ,(V_total - (V_occupied - V_overlap))/(m_total*u)*10**(0),'cm^3/g'
  write(19,fmt='(1X a,f10.3,1X a)') 'The mass of unit cell is                        : ',m_total*u,'10**-27 kg'
  write(19,*) 'Total CPU time: ',finish-start,'s'
  close(19)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Evaluation method 2 : Uniform grid. Evaluate each point as to whether it is occupied, void !!
!!                       or accessible                                                        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

else if (eval_method == 2) then                                                                ! choosing second evaluation method
! Initialization
  write(6,*) ' '
  write(6,*) '########## PROBE RADIUS #####################################'
  write(6,*) 'Probe radius to evaluate the accessible volume (in angstrom):'                    ! ask for the probe radius
  write(6,*) '#############################################################'
  read(5,*) probe_r                                                                             ! define the probe radius
  write(6,*) ' '
  write(6,*) '########## CELL VECTORS ARE: ################################################'    ! display the cell vectors. Allowing an initial assessment of the number of grid points.
  write(6,*) 'a = ',cell_a,'with length ',sqrt(cell_a(1)**2 + cell_a(2)**2 + cell_a(3)**2)
  write(6,*) 'b = ',cell_b,'with length ',sqrt(cell_b(1)**2 + cell_b(2)**2 + cell_b(3)**2)
  write(6,*) 'c = ',cell_c,'with length ',sqrt(cell_c(1)**2 + cell_c(2)**2 + cell_c(3)**2)
  write(6,*) '#############################################################################'
  write(6,*) ' '
  write(6,*) '########## METHOD TO DEFINE GRID POINTS ##########################'
  write(6,*) '(1) Grid points per unit cell vector                  (3 integers)'
  write(6,*) '(2) Grid points per angstroem (uniform grid)              (1 real)'
  write(6,*) '##################################################################'
  read(5,*) t                                                                                  ! t serves as a temporary variable
  if (t == 1) then
    write(6,*) ' '
    write(6,*) '########## CHOOSE NUMBER OF GRID POINTS ################'
    write(6,*) 'How many grid points per unit cell vector (3 integers) ?'                      ! ask for the number of grid points along each cell vector
    write(6,*) '########################################################'
    read(5,*) grid_a, grid_b, grid_c
  else if (t == 2) then
    write(6,*) ' '
    write(6,*) '########## CHOOSE NUMBER OF GRID POINTS #####'
    write(6,*) 'How many grid points per angstroem (1 real) ?'                              ! ask for the number of grid points per angstroem for each unit cell vector
    write(6,*) '#############################################'
    read(5,*) g
    grid_a = ceiling(g*sqrt(cell_a(1)**2 + cell_a(2)**2 + cell_a(3)**2))                       ! ceiling -> round to the next higher integer. Use g * len_unit_cell_vector as the number of grid points
    grid_b = ceiling(g*sqrt(cell_b(1)**2 + cell_b(2)**2 + cell_b(3)**2))
    grid_c = ceiling(g*sqrt(cell_c(1)**2 + cell_c(2)**2 + cell_c(3)**2))
  end if

  call cpu_time(start)                                                                         ! initialize time measurement

! Evaluate the total mass (for later evaluation)
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
          if ((elements(n_coords) == 'H'  .and. dist_point_atom <= 1.20) .or. &                    ! if the grid point is inside any atom (distance is smaller than the vdW radius of the respective atom)
              (elements(n_coords) == 'C'  .and. dist_point_atom <= 1.70) .or. &                    ! add this point to the occupied list
              (elements(n_coords) == 'N'  .and. dist_point_atom <= 1.55) .or. &
              (elements(n_coords) == 'O'  .and. dist_point_atom <= 1.52) .or. &
              (elements(n_coords) == 'Ni' .and. dist_point_atom <= 1.63) .or. &
              (elements(n_coords) == 'Cu' .and. dist_point_atom <= 1.40) .or. &
              (elements(n_coords) == 'Zn' .and. dist_point_atom <= 1.39) .or. &
              (elements(n_coords) == 'Zr' .and. dist_point_atom <= 2.36)) then
            n_occ = n_occ + 1                                                                      ! increase assignemnt counter for the occupied list
            exit loop14                                                                            ! stop looping through the atoms at this point. The point is already determined as occupied. Avoid double counting!

          else if ((elements(n_coords) == 'H'  .and. dist_point_atom >= 1.20 + probe_r) .or. &     ! if grid point is outside an atom + the probe radius -> Clearly accessible
                   (elements(n_coords) == 'C'  .and. dist_point_atom >= 1.70 + probe_r) .or. &     ! add +1 to the counter 'counter_acc' AND to counter_noOccu
                   (elements(n_coords) == 'N'  .and. dist_point_atom >= 1.55 + probe_r) .or. &
                   (elements(n_coords) == 'O'  .and. dist_point_atom >= 1.52 + probe_r) .or. &
                   (elements(n_coords) == 'Ni' .and. dist_point_atom >= 1.63 + probe_r) .or. &
                   (elements(n_coords) == 'Cu' .and. dist_point_atom >= 1.40 + probe_r) .or. &
                   (elements(n_coords) == 'Zn' .and. dist_point_atom >= 1.39 + probe_r) .or. &
                   (elements(n_coords) == 'Zr' .and. dist_point_atom >= 2.36 + probe_r)) then
            counter_access = counter_access + 1
            counter_noOccu = counter_noOccu + 1
          else                                                                                     ! If outside an atom -> not occupied
            counter_noOccu = counter_noOccu + 1
          end if
        end do loop14                                                                                    ! end do for all atoms

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
  factor = 1.0 + 1.0/((grid_per_A_x+grid_per_A_y+grid_per_A_z)/3.0)                            ! divide by average grid points per A

  !
  ! Get number of sub-grid points per atom
  !
!
! HERE: Idea -> include the firstloop here into the very first loop over the grid points -> avoid looping over all accessible points twice
!! I.e. do the counting part in the very first loop. Then, in the second loop write grid points

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
      if ((elements(n_coords) == 'H'  .and. dist_point_atom < 1.20 + probe_r*factor) .or. &             ! if grid point is outside an atom + the probe radius, but close by -> get this point
          (elements(n_coords) == 'C'  .and. dist_point_atom < 1.70 + probe_r*factor) .or. &             ! to check accessibility later on
          (elements(n_coords) == 'N'  .and. dist_point_atom < 1.55 + probe_r*factor) .or. &             ! take only points in between vdW+probe_r AND vdw+probe_r*factor
          (elements(n_coords) == 'O'  .and. dist_point_atom < 1.52 + probe_r*factor) .or. &
          (elements(n_coords) == 'Ni' .and. dist_point_atom < 1.63 + probe_r*factor) .or. &
          (elements(n_coords) == 'Cu' .and. dist_point_atom < 1.40 + probe_r*factor) .or. &
          (elements(n_coords) == 'Zn' .and. dist_point_atom < 1.39 + probe_r*factor) .or. &
          (elements(n_coords) == 'Zr' .and. dist_point_atom < 2.36 + probe_r*factor)) then
        sub_division(n_coords)%sub_grid_points = sub_division(n_coords)%sub_grid_points + 1
!        exit loop13  ! TEST: include double counting of points. Should matter for this
      end if
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
      if ((elements(n_coords) == 'H'  .and. dist_point_atom < 1.20 + probe_r*factor) .or. &             ! if grid point is outside an atom + the probe radius, but close by -> get this point
          (elements(n_coords) == 'C'  .and. dist_point_atom < 1.70 + probe_r*factor) .or. &             ! to check accessibility later on
          (elements(n_coords) == 'N'  .and. dist_point_atom < 1.55 + probe_r*factor) .or. &             ! take only points in between vdW+probe_r AND vdw+probe_r*factor
          (elements(n_coords) == 'O'  .and. dist_point_atom < 1.52 + probe_r*factor) .or. &
          (elements(n_coords) == 'Ni' .and. dist_point_atom < 1.63 + probe_r*factor) .or. &
          (elements(n_coords) == 'Cu' .and. dist_point_atom < 1.40 + probe_r*factor) .or. &
          (elements(n_coords) == 'Zn' .and. dist_point_atom < 1.39 + probe_r*factor) .or. &
          (elements(n_coords) == 'Zr' .and. dist_point_atom < 2.36 + probe_r*factor)) then
        sub_division(n_coords)%sub_grid_points = sub_division(n_coords)%sub_grid_points + 1
        sub_division(n_coords)%sub_grids(sub_division(n_coords)%sub_grid_points,:) = list_access(n,:)
!        exit loop15   ! TEST
      end if
    end do loop15
  end do

! Write to screen
  write(6,*) ' '
  write(6,*) 'N_acc after 1st loop                ',n_access
  write(6,*) '  Total number of unoccupied points ',n_noOccu
  write(6,*) '  Sum of sub-grid points            ',sum(sub_division(:)%sub_grid_points)

! Write to file
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

      ! condition for NOT occpupied, but not necessarily accessible
      if ((elements(n_coords) == 'H'  .and. dist_point_atom < 1.20 + probe_r) .or. &         ! If grid point is outside an atom (see last if statement), but within a length of the probe radius
          (elements(n_coords) == 'C'  .and. dist_point_atom < 1.70 + probe_r) .or. &    
          (elements(n_coords) == 'N'  .and. dist_point_atom < 1.55 + probe_r) .or. &
          (elements(n_coords) == 'O'  .and. dist_point_atom < 1.52 + probe_r) .or. &
          (elements(n_coords) == 'Ni' .and. dist_point_atom < 1.63 + probe_r) .or. &
          (elements(n_coords) == 'Cu' .and. dist_point_atom < 1.40 + probe_r) .or. &
          (elements(n_coords) == 'Zn' .and. dist_point_atom < 1.39 + probe_r) .or. &
          (elements(n_coords) == 'Zr' .and. dist_point_atom < 2.36 + probe_r)) then
        do f = 1, sub_division(n_coords)%sub_grid_points                                     ! Go through the 'check accessibility' points (which have been determined before)
          dist_point_point = sqrt(sum((sub_division(n_coords)%sub_grids(f,:) - &
                                       list_noOccu(a,:))**2)) ! determine the distance to any accessible point
          if (dist_point_point < probe_r) then                                               ! if the distance to any accessible point is smaller than the probe radius -> NEW ACCESSIBLE POINT
            n_access = n_access + 1                                                          ! increase accessible counter
            exit loop1                                                                       ! Stop looping once this is confirmed (i.e. stop looping over the atoms)
          end if
        end do
      end if

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
  write(6,fmt='(1X,a,2X,I15)') 'Total number of grid points: ',grid_a*grid_b*grid_c
  write(6,fmt='(1X,a,f15.3,a,3(f10.3,2X,a),a)') 'Grid point density:         ',grid_a*grid_b*grid_c/V_total, &
              ' grid points per A^3, with ',grid_per_A_x,' x ',grid_per_A_y,' x ',grid_per_A_z,' grid points per A'                           ! (grid_a*grid_b*grid_c/V_total)**(1./3.)
  write(6,*) 'Points OCCUPIED:             ',n_occ
  write(6,*) 'Points NOT OCCUPIED (void):  ',grid_a*grid_b*grid_c - n_occ
  write(6,*) 'Points ACCESSIBLE:           ',n_access
  write(6,fmt='(1X a,7X f7.3,1X a)') 'probe radius:                ',probe_r,'A'
  write(6,*) ' '
 
  write(6,*) 'Porosity (void):             ',(real(grid_a*grid_b*grid_c) - real(n_occ))/(real(grid_a*grid_b*grid_c))*100,'%'
  write(6,*) 'Porosity (accessible):       ',real(n_access)/(real(grid_a*grid_b*grid_c))*100,'%'
 
  V_void       = (real(grid_a*grid_b*grid_c) - real(n_occ))/(real(grid_a*grid_b*grid_c))*V_total
  V_accessible = real(n_access)/(real(grid_a*grid_b*grid_c))*V_total 
 
  write(6,*) 'Volume (void):               ',V_void,'A^3'
  write(6,*) 'Volume (accessible):         ',V_accessible,'A^3'
  write(6,*) ' '
 
  write(6,fmt='(1X a,f10.3,a)') 'Unit cell volume (V_total)                   : ',V_total,' A^3'
  write(6,fmt='(1X a,f10.3,a)') 'Mass of unit cell (m_total)                  : ',m_total*u,' 10**-27 kg'
  write(6,fmt='(1X a,f10.3,a)') 'Density of the structure (m_total/V_total)   : ',m_total*u/V_total*10**3,' kg/m^3'
  write(6,fmt='(1X a,f10.3,a)') 'Inverse pore volume density (V_void/m_total) : ',V_void/(m_total*u)*10**(0),' cm^3/g'
  write(6,fmt='(1X a,f10.3,a)') 'Inverse pore volume density (V_acc/m_total)  : ',V_accessible/(m_total*u)*10**(0),' cm^3/g'
 
  call cpu_time(finish)
  write(6,*) ' '
  write(6,*) 'Total CPU time was ',finish-start,'s'
 
  ! 
  ! Write to file
  ! 
 ! WRITE OUTPUT
 !
  write(19,*) ' '
  write(19,*) '#######################################################################################################'
  write(19,*) '###############################################  OUTPUT ###############################################'
  write(19,*) '#######################################################################################################'
  write(19,*) 'Structure : ',name_struct
  write(19,fmt='(1X a,7X I4,1X a,1X I4,1X a,1x I4)') 'Grid which was used:         ',grid_a,'x',grid_b,'x',grid_c                              ! 1X -> 1 space, 7X -> 7 spaces
  write(19,fmt='(1X,a,2X,I15)') 'Total number of grid points: ',grid_a*grid_b*grid_c
  write(19,fmt='(1X,a,f15.3,a,3(f10.3,2X,a),a)') 'Grid point density:         ',grid_a*grid_b*grid_c/V_total, &
              ' grid points per A^3, with ',grid_per_A_x,' x ',grid_per_A_y,' x ',grid_per_A_z,' grid points per A'                           ! (grid_a*grid_b*grid_c/V_total)**(1./3.)
  write(19,*) 'Points OCCUPIED:             ',n_occ
  write(19,*) 'Points NOT OCCUPIED (void):  ',grid_a*grid_b*grid_c - n_occ
  write(19,*) 'Points ACCESSIBLE:           ',n_access
  write(19,fmt='(1X a,7X f7.3,1X a)') 'probe radius:                ',probe_r,'A'
  write(19,*) ' '

  write(19,*) 'Porosity (void):             ',(real(grid_a*grid_b*grid_c) - real(n_occ))/(real(grid_a*grid_b*grid_c))*100,'%'
  write(19,*) 'Porosity (accessible):       ',real(n_access)/(real(grid_a*grid_b*grid_c))*100,'%'

  V_void       = (real(grid_a*grid_b*grid_c) - real(n_occ))/(real(grid_a*grid_b*grid_c))*V_total
  V_accessible = real(n_access)/(real(grid_a*grid_b*grid_c))*V_total

  write(19,*) 'Volume (void):               ',V_void,'A^3'
  write(19,*) 'Volume (accessible):         ',V_accessible,'A^3'
  write(19,*) ' '

  write(19,fmt='(1X a,f10.3,a)') 'Unit cell volume (V_total)                   : ',V_total,' A^3'
  write(19,fmt='(1X a,f10.3,a)') 'Mass of unit cell (m_total)                  : ',m_total*u,' 10**-27 kg'
  write(19,fmt='(1X a,f10.3,a)') 'Density of the structure (m_total/V_total)   : ',m_total*u/V_total*10**3,' kg/m^3'
  write(19,fmt='(1X a,f10.3,a)') 'Inverse pore volume density (V_void/m_total) : ',V_void/(m_total*u)*10**(0),' cm^3/g'
  write(19,fmt='(1X a,f10.3,a)') 'Inverse pore volume density (V_acc/m_total)  : ',V_accessible/(m_total*u)*10**(0),' cm^3/g'

  write(19,*) ' '
  write(19,*) 'Total CPU time was ',finish-start,'s'
  close(19)

  deallocate(grid_points)
!  deallocate(list_occupi)
  deallocate(list_access)
  deallocate(list_check_acc)
  deallocate(list_noOccu)
  deallocate(sub_division)
end if           ! global end if (second method ends here)

deallocate(elements)
deallocate(coordinates)

end program porosity





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine evaluation the occupied volume and the masses of the atoms !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eval_vol_mass(element,vocc,m)           ! element as input, V_occ and mass as output
  character(2), intent(in)      :: element         ! element symbol
  real(8), intent(inout)        :: vocc            ! occupied volume of the atoms (according to their vdW radii)
  real(8), intent(inout)        :: m               ! mass of all atoms
  real(8), parameter            :: pi = 3.14159265358979323846  ! define pi
  if (element == 'H') then                         ! if a hydrogen is found,  use r_vdw = 1.20 A and m_atom = 1.0079 (in u)
   vocc = vocc + 4.0/3.0*pi*(1.20)**3
   m    = m    + 1.0079
  else if (element == 'C') then                    ! if a carbon is found,    use r_vdw = 1.70 A and m_atom = 12.011 (in u)
   vocc = vocc + 4.0/3.0*pi*(1.70)**3
   m    = m    + 12.011
  else if (element == 'N') then                    ! if a nitrogen is found,  use r_vdw = 1.55 A and m_atom = 14.007 (in u)
   vocc = vocc + 4.0/3.0*pi*(1.55)**3
   m    = m    + 14.007
  else if (element == 'O') then                    ! if a oxygen is found,    use r_vdw = 1.52 A and m_atom = 15.999 (in u)
   vocc = vocc + 4.0/3.0*pi*(1.52)**3
   m    = m    + 15.999
  else if (element == 'Ni') then                   ! if a nickel is found,    use r_vdw = 1.63 A and m_atom = 58.693 (in u)
   vocc = vocc + 4.0/3.0*pi*(1.63)**3
   m    = m    + 58.693
  else if (element == 'Cu') then                   ! if a copper is found,    use r_vdw = 1.40 A and m_atom = 63.546 (in u)
   vocc = vocc + 4.0/3.0*pi*(1.40)**3
   m    = m    + 63.546
  else if (element == 'Zr') then                   ! if a zirconium is found, use r_vdw = 2.36 A and m_atom = 91.224 (in u)
   vocc = vocc + 4.0/3.0*pi*(2.36)**3
   m    = m    + 91.224
  else if (element == 'Zn') then                   ! if a zinc is found,      use r_vdw = 1.39 A and m_atom = 65.390 (in u)
   vocc = vocc + 4.0/3.0*pi*(1.39)**3
   m    = m    + 65.390
  end if
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
 
                                                                                              ! Evaluate, if distance between the atoms is smaller than the sum of the covalent radii
  if ((element_a == 'H' .and. element_b == 'H') .and. (dist_ab < 0.33 + 0.33)) then           ! two H
   r_vdw1 = 1.20
   r_vdw2 = 1.20
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if ((element_a == 'C'  .and. element_b == 'C')  .and. (dist_ab < 0.76 + 0.76)) then    ! two C
   r_vdw1 = 1.70
   r_vdw2 = 1.70
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if ((element_a == 'N'  .and. element_b == 'N')  .and. (dist_ab < 0.71 + 0.71)) then    ! two N
   r_vdw1 = 1.55
   r_vdw2 = 1.55
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if ((element_a == 'O'  .and. element_b == 'O')  .and. (dist_ab < 0.66 + 0.66)) then    ! two O
   r_vdw1 = 1.52
   r_vdw2 = 1.52
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if ((element_a == 'Ni' .and. element_b == 'Ni') .and. (dist_ab < 1.24 + 1.24)) then    ! two Ni
   r_vdw1 = 1.63
   r_vdw2 = 1.63
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if ((element_a == 'Cu' .and. element_b == 'Cu') .and. (dist_ab < 1.30 + 1.30)) then    ! two Cu
   r_vdw1 = 1.40
   r_vdw2 = 1.40
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if ((element_a == 'Zn' .and. element_b == 'Zn') .and. (dist_ab < 1.33 + 1.33)) then    ! two Zn
   r_vdw1 = 1.39
   r_vdw2 = 1.39
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if ((element_a == 'Zr' .and. element_b == 'Zr') .and. (dist_ab < 1.48 + 1.48)) then    ! two Zr
   r_vdw1 = 2.36
   r_vdw2 = 2.36
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
 ! ALL PAIRS WITH H
  else if (((element_a == 'H'  .and. element_b == 'C')  .or. &
            (element_a == 'C'  .and. element_b == 'H')) .and. (dist_ab < 0.33 + 0.76)) then   ! one H, one C
   r_vdw1 = 1.20
   r_vdw2 = 1.70
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if (((element_a == 'H'  .and. element_b == 'N')  .or. &
            (element_a == 'N'  .and. element_b == 'H')) .and. (dist_ab < 0.33 + 0.71)) then   ! one H, one N
   r_vdw1 = 1.20
   r_vdw2 = 1.55
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if (((element_a == 'H'  .and. element_b == 'O')  .or. &
            (element_a == 'O'  .and. element_b == 'H')) .and. (dist_ab < 0.33 + 0.66)) then   ! one H, one O
   r_vdw1 = 1.20
   r_vdw2 = 1.52
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if (((element_a == 'H'  .and. element_b == 'Ni') .or. &
            (element_a == 'Ni' .and. element_b == 'H')) .and. (dist_ab < 0.33 + 1.24)) then   ! one H, one Ni
   r_vdw1 = 1.20
   r_vdw2 = 1.63
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if (((element_a == 'H'  .and. element_b == 'Cu') .or. &
            (element_a == 'Cu' .and. element_b == 'H')) .and. (dist_ab < 0.33 + 1.30)) then   ! one H, one Cu
   r_vdw1 = 1.20
   r_vdw2 = 1.40
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if (((element_a == 'H'  .and. element_b == 'Zn') .or. &
            (element_a == 'Zn' .and. element_b == 'H')) .and. (dist_ab < 0.33 + 1.33)) then   ! one H, one Zn
   r_vdw1 = 1.20
   r_vdw2 = 1.39
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if (((element_a == 'H'  .and. element_b == 'Zr') .or. &
            (element_a == 'Zr' .and. element_b == 'H')) .and. (dist_ab < 0.33 + 1.48)) then   ! one H, one Zr
   r_vdw1 = 1.20
   r_vdw2 = 2.36
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
 ! ALL PAIRS WITH C
  else if (((element_a == 'C'  .and. element_b == 'N')  .or. &
            (element_a == 'N'  .and. element_b == 'C')) .and. (dist_ab < 0.76 + 0.71)) then   ! one C, one N
   r_vdw1 = 1.70
   r_vdw2 = 1.55
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if (((element_a == 'C'  .and. element_b == 'O')  .or. &
            (element_a == 'O'  .and. element_b == 'C')) .and. (dist_ab < 0.76 + 0.66)) then   ! one C, one O
   r_vdw1 = 1.70
   r_vdw2 = 1.52
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if (((element_a == 'C'  .and. element_b == 'Ni') .or. &
            (element_a == 'Ni' .and. element_b == 'C')) .and. (dist_ab < 0.76 + 1.24)) then   ! one C, one Ni
   r_vdw1 = 1.70
   r_vdw2 = 1.63
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if (((element_a == 'C'  .and. element_b == 'Cu') .or. &
            (element_a == 'Cu' .and. element_b == 'C')) .and. (dist_ab < 0.76 + 1.24)) then   ! one C, one Cu
   r_vdw1 = 1.70
   r_vdw2 = 1.40
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if (((element_a == 'C'  .and. element_b == 'Zn') .or. &
            (element_a == 'Zn' .and. element_b == 'C')) .and. (dist_ab < 0.76 + 1.33)) then   ! one C, one Zn
   r_vdw1 = 1.70
   r_vdw2 = 1.39
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if (((element_a == 'C'  .and. element_b == 'Zr') .or. &
            (element_a == 'Zr' .and. element_b == 'C')) .and. (dist_ab < 0.76 + 1.48)) then   ! one C, one Zr
   r_vdw1 = 1.70
   r_vdw2 = 2.36
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
 ! ALL PAIRS WITH N
  else if (((element_a == 'N'  .and. element_b == 'O')  .or. &
            (element_a == 'O'  .and. element_b == 'N')) .and. (dist_ab < 0.71 + 0.66)) then   ! one N, one O
   r_vdw1 = 1.55
   r_vdw2 = 1.52
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if (((element_a == 'N'  .and. element_b == 'Ni') .or. &
            (element_a == 'Ni' .and. element_b == 'N')) .and. (dist_ab < 0.71 + 1.24)) then   ! one N, one Ni
   r_vdw1 = 1.55
   r_vdw2 = 1.63
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if (((element_a == 'N'  .and. element_b == 'Cu') .or. &
            (element_a == 'Cu' .and. element_b == 'N')) .and. (dist_ab < 0.71 + 1.30)) then   ! one N, one Cu
   r_vdw1 = 1.55
   r_vdw2 = 1.40
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if (((element_a == 'N'  .and. element_b == 'Zn') .or. &
            (element_a == 'Zn' .and. element_b == 'N')) .and. (dist_ab < 0.71 + 1.33)) then   ! one C, one Zn
   r_vdw1 = 1.55
   r_vdw2 = 1.39
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if (((element_a == 'N' .and. element_b == 'Zr') .or. &
            (element_a == 'Zr' .and. element_b == 'N')) .and. (dist_ab < 0.71 + 1.48)) then   ! one C, one Zr
   r_vdw1 = 1.55
   r_vdw2 = 2.36
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
 ! ALL PAIRS WITH O
  else if (((element_a == 'O'  .and. element_b == 'Ni') .or. &
            (element_a == 'Ni' .and. element_b == 'O')) .and. (dist_ab < 0.66 + 1.24)) then   ! one O, one Ni
   r_vdw1 = 1.52
   r_vdw2 = 1.63
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if (((element_a == 'O'  .and. element_b == 'Cu') .or. &
            (element_a == 'Cu' .and. element_b == 'O')) .and. (dist_ab < 0.66 + 1.30)) then   ! one O, one Cu
   r_vdw1 = 1.52
   r_vdw2 = 1.40
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if (((element_a == 'O'  .and. element_b == 'Zn') .or. &
            (element_a == 'Zn' .and. element_b == 'O')) .and. (dist_ab < 0.66 + 1.33)) then   ! one O, one Zn
   r_vdw1 = 1.52
   r_vdw2 = 1.39
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  else if (((element_a == 'O'  .and. element_b == 'Zr') .or. &
            (element_a == 'Zr' .and. element_b == 'O')) .and. (dist_ab < 0.66 + 1.48)) then   ! one O, one Zr
   r_vdw1 = 1.52
   r_vdw2 = 2.36
   call overlap(r_vdw1, r_vdw2, dist_ab, sub_over)
  end if
  return
end subroutine eval_overlap
