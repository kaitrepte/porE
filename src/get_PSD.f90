program porefinder

implicit none

! pore_finder

character(2)                        :: struct
character(len=25)                   :: name_struct

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
real(8)                             :: distance1, distance2    ! corresponding minimum distance to any atom
real(8)                             :: tmp_dist, vdw           ! temporary distance, vdW radius of the atom

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
write(6,*) 'HKUST-1                  - h1'
write(6,*) 'Benzene, opt             - be'
write(6,*) 'Benzene, exp             - b2'
write(6,*) 'Benzene, C only          - bc'
write(6,*) 'H atom                   - ha'
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Monte-Carlo to get pore sizes !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
stepsize     = 0.01
allocate(all_distances(start_points))
allocate(all_distances2(start_points))

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
      if (elements(n) == 'C')  vdw = 1.70
      if (elements(n) == 'N')  vdw = 1.55
      if (elements(n) == 'O')  vdw = 1.52
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
!    write(6,*) "Start ",a," Final distance ",distance1*2.0   ! write diameter, not radius
    all_distances(a) = distance1*2.0
  else
!    write(6,*) "Start ",a," Final distance ",distance2*2.0   ! write diameter, not radius
    all_distances(a) = distance2*2.0
  end if
end do      ! end starting points


all_distances2(:) = all_distances(:)

! Get distribution
write(6,*) ' '
write(6,*) 'Pore size distribution for ',name_struct

do a = 1, start_points
  c = 0
  do b = 1, start_points
    if (abs(all_distances(a) - all_distances2(b)) < 0.10) then    ! collect data which is within this range of the value
      c = c + 1
      all_distances2(b) = 1000.0                  ! do not evaluate this point again
    end if
  end do
  if (c == 0) then
  else
    write(6,*) all_distances(a), c
  end if
end do

call cpu_time(finish)
write(6,*) 'Total time ',finish-start,' s'


deallocate(coordinates)
deallocate(elements)
deallocate(all_distances)
deallocate(all_distances2)

end program porefinder
