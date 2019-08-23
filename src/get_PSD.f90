program porefinder

implicit none

! pore_finder

character(2)                        :: struct
character(len=100)                  :: name_struct

real(8)                             :: start, finish       ! timing

integer                             :: number_of_atoms
real(8)                             :: cell_a(3)           ! array for the cell vector in the a direction. Vector.
real(8)                             :: cell_b(3)           ! array for the cell vector in the b direction. Vector.
real(8)                             :: cell_c(3)           ! array for the cell vector in the c direction. Vector.
real(8), allocatable, dimension(3)  :: coordinates(:,:)    ! array for the coordinates. Matrix.
character(2), allocatable           :: elements(:)         ! array for the elements. Vector.

integer                             :: a,b,c,d,e,f,n,t         ! loop parameter
real(8)                             :: rand1, rand2, rand3     ! random numbers to get new coordinates

real(8)                             :: coords1(3), coords2(3)  ! coordinates of point before and after MC step
real(8), allocatable                :: coords_all_cart(:,:)    ! all coordinates of point after MC, cartesian
real(8), allocatable                :: coords_all_frac(:,:)    ! all coordinates of point after MC, fractional
real(8)                             :: distance1, distance2    ! corresponding minimum distance to any atom
real(8)                             :: tmp_dist, vdw           ! temporary distance, vdW radius of the atom
real(8)                             :: distribution            ! final PSD distribution for a given radius

integer                             :: start_points, cycles    ! number of starting points, number of MC cycles
real(8)                             :: stepsize                ! step size for MC steps
real(8), allocatable                :: all_distances(:)        ! store all distances (maybe need another list to separate different distances which occur more often.. PSD and stuff)
real(8), allocatable                :: all_distances2(:)       ! store all distances, to double check

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
write(6,*) 'HKUST-1, open Cu sites   - h1'
write(6,*) 'HKUST-1, O-Cu-Cu-O       - ho'
write(6,*) 'C60@MOF                  - c6'
write(6,*) 'Benzene, opt             - be'
write(6,*) 'Benzene, exp             - b2'
write(6,*) 'Benzene, C only          - bc'
write(6,*) 'H atom                   - ha'
write(6,*) 'User-defined xyz         - ud'
write(6,*) '#############################'
read(5,*) struct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read in the xyz coordinates and the cell vectors !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
else if (struct == 'ho') then                                                                 ! if HKUST-1 (primitive cell) is choosen
  open(unit=15,file='../structures_xyz/hkust1_with_O.xyz',status='old',action='read')             ! read in the xyz file
  name_struct = 'HKUST-1'
else if (struct == 'c6') then                                                                 ! if HKUST-1 (primitive cell) is choosen
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
else if (struct == 'ud') then                                                                 ! if user-defined cell is choosen
  write(6,*) 'Provide the path and the name of the xyz file (e.g. "../structures_xyz/test.xyz")'
  read(5,*) name_struct 
  open(unit=15,file=name_struct,status='old',action='read')                                  ! read in the xyz file
  write(6,*) 'Provide a name for you structure'
  read(5,*) name_struct
end if

! Initialize stuff
write(6,*) 'How many starting points (recommended: >= 100)?'
read(5,*) start_points
write(6,*) 'How many Monte-Carlo steps (recommended: >= 10000)?'
read(5,*) cycles


call cpu_time(start)

! Read in the corresponding values
read(unit=15,fmt='(I13.0)') number_of_atoms                     ! first entry is the number of atoms
read(unit=15,fmt=*) cell_a(1:3), cell_b(1:3), cell_c(1:3)       ! second entry contains the cell vectors. Read them in individually (makes it easier later on)

allocate(elements(number_of_atoms))                             ! allocate (number_of_atoms) fields for elements. There is one elements each. As many elements as number_of_atoms (makes sense :))
allocate(coordinates(number_of_atoms,3))                        ! allocate (number_of_atoms) fields for coordinates. There are 3 coordinates per entry. 
do n = 1,number_of_atoms                                        ! go through all atoms 
  read(unit=15,fmt=*) elements(n), coordinates(n,1:3)           ! storing element and coordinates
end do
close(unit=15)                                                  ! close the file (no longer necessary)

! Output file
open(unit=19,file='output',status='unknown',action='write')
write(19,*) 'Starting points   (recommended: >= 100)   : ',start_points
write(19,*) 'Monte-Carlo steps (recommended: >= 10000) : ',cycles
write(19,*) ' '

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Monte-Carlo to get pore sizes !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
stepsize     = 0.01
allocate(all_distances(start_points))
allocate(all_distances2(start_points))
allocate(coords_all_cart(start_points,3))
allocate(coords_all_frac(start_points,3))

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

! Get probe diameter
  if (distance1 > distance2) then
!    write(6,*) "Start ",a," Final distance ",distance1*2.0   ! write diameter, not radius. Store position as well
    all_distances(a) = distance1*2.0
    coords_all_cart(a,:) = coords2(:)       
  else
!    write(6,*) "Start ",a," Final distance ",distance2*2.0   ! write diameter, not radius. Store position as well
    all_distances(a) = distance2*2.0
    coords_all_cart(a,:) = coords1(:)
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

do a = 1, start_points
  distribution = 0.0
  do b = 1, start_points
    if (abs(all_distances(a) - all_distances2(b)) < 0.10) then    ! collect data which is within this range of the value
      distribution = distribution + 1.0
      all_distances2(b) = 1000.0                  ! do not evaluate this point again
    end if
  end do
  distribution = distribution/start_points*100.D0  ! distribution in %
  if (distribution < 5.0) then      ! if less than 5 % -> do not evaluate
  else
    write(6,fmt='(F12.6,F8.2,11X,3F12.6,2X,3F12.6)') all_distances(a), distribution, coords_all_cart(a,:), coords_all_frac(a,:)
    write(19,fmt='(F12.6,F8.2,11X,3F12.6,2X,3F12.6)') all_distances(a), distribution, coords_all_cart(a,:), coords_all_frac(a,:)
  end if
end do

call cpu_time(finish)
write(6,*) ' '
write(6,*) 'Total time ',finish-start,' s'

write(19,*) ' '
write(19,*) 'Total time ',finish-start,' s'

close(19)

deallocate(coordinates)
deallocate(elements)
deallocate(all_distances)
deallocate(all_distances2)
deallocate(coords_all_cart)
deallocate(coords_all_frac)

end program porefinder


subroutine cart_to_frac(vecA,vecB,vecC,pos_cart,pos_frac)
! Transform cart to fractional. Make sure that point are within the unit cell. Transform back to cartesian. Return both cart and frac
real(8), intent(in)    :: vecA(3), vecB(3), vecC(3)
real(8), intent(inout) :: pos_cart(3)
real(8), intent(inout) :: pos_frac(3)
real(8)                :: trans_matrix(3,3)
real(8)                :: determinant
real(8)                :: lenA, lenB, lenC, angleBC, angleAC, angleAB
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
end subroutine cart_to_frac
