program porefinder

implicit none

! pore_finder

character(2)                        :: struct

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

call cpu_time(start)

open(unit=15,file='structures_xyz/mof5.xyz',status='old',action='read')

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
! Initialize stuff
start_points = 100
cycles       = 10000
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
    write(6,*) "Start ",a," Final distance ",distance1*2.0   ! write diameter, not radius
    all_distances(a) = distance1*2.0
  else
    write(6,*) "Start ",a," Final distance ",distance2*2.0   ! write diameter, not radius
    all_distances(a) = distance2*2.0
  end if
end do      ! end starting points


all_distances2(:) = all_distances(:)

! Get distribution
write(6,*) ' '
write(6,*) 'Pore size distribution'

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
