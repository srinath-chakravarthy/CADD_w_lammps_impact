! ---- Module to map lammps atoms to local atoms
module mod_lammps
  use mod_global
  use lammps

  implicit none
  ! ----- Predefined types for retrieval from lammps
  ! ---- Compute and fix values
  integer, parameter :: global_data_style = 0, peratom_style = 1, local_data_style = 2
  integer, parameter :: scalar_type = 0, vector_type = 1, array_type = 2

  integer :: kind4int
  integer :: n_lammps_atoms
  !< contains map of lammps atom_index to cadd atom index
  !<     lammps_to cadd --> any array x(i) = x_lammps(lammps_cadd_map(i))
  !< These arrays contain mapping to lammps_extract_atom and extract_fixes etc. 
  integer, dimension(:), allocatable :: lammps_cadd_map, cadd_lammps_map
  integer, dimension(:), allocatable :: lammps_pad_map, pad_lamms_map

  
  !< Gmap arrays contain mapping to lammps_gather and lammps_scatter 
  integer, dimension(:), allocatable :: lammps_cadd_gmap, cadd_lammps_gmap
  
  !< Main coordinate, force and velocity arrays extracted using extract_atom
  real (C_double), dimension(:,:), pointer :: lammps_coord => NULL()
  real (C_double), dimension(:,:), pointer :: lammps_velocity => NULL()
  real (C_double), dimension(:,:), pointer :: lammps_force => NULL()

  
  !< Other main items such as energies
  !< I assume these are coming from computes
  !< To be defined 
  
  !< Various computes values
  !< All computes are extracted using per-atom quantities except where listed explicitly
  !< 1) Current differential displacement from original coordinates lammps_dx ---> mapped to CADD Atomdispl array
  !< 2) Current Averaged displacement over specified interval lammps_avg_dx ---> mapped to CADD avedispl array

  real (C_double), dimension(:,:), pointer :: compute_lammps_dx => NULL()
  

  
  real (C_double), dimension(:,:), pointer :: compute_lammps_stress => NULL()

  
  
contains
  subroutine map_from_lammps(x,id,lmp)
    implicit none
    !<--- Input parameters ----
    integer :: id(NDF,*)
    double precision :: x(nxdm,numnp)
    type(C_ptr) :: lmp
    !< -------------------------

    ! ---- Local Variables 
    character (len=1024) :: command_line

    integer :: error, narg, nsize

    double precision :: tol
    integer :: iatom, catom

    double precision, dimension(:), allocatable :: r
    double precision, dimension(:,:), allocatable :: rcoords

    
    tol = 1.e-5
    

    if (C_ASSOCIATED(lmp)) then
       ! ---- Obtain pointers to internal lammps data ---- 
       call lammps_extract_atom(lammps_coord, lmp, 'x')
!!$       print *, 'Size of CADD domain = ', size(x), shape(x)

       ! ---- Create mapping to CADD -----------
       n_lammps_atoms = size(lammps_coord,2)
       allocate(cadd_lammps_map(n_lammps_atoms))
       allocate(lammps_cadd_map(numnp))
       
       ! Initialize the arrays to large, known, negative value for error checking
       ! default integer is of kind4, max value of 2147483647...so use intrisic function
       ! huge of a kind4 integer (value doesn't 
       kind4int = 0;
       cadd_lammps_map = -HUGE(kind4int)
       lammps_cadd_map = -HUGE(kind4int)
       
!!$       print *, 'Size of arrays = ', size(cadd_lammps_map), size(lammps_cadd_map)
!!$
!!$       print *, 'Number of Lammps Atoms = ', n_lammps_atoms
       do catom = 1, numnp
!!$          print *,"Lammps coord = ", lammps_coord(1,iatom), lammps_coord(2,iatom)
	  if (isRelaxed(catom) /= 0) then 
	    do iatom = 1, n_lammps_atoms
                if (abs(lammps_coord(1,iatom)-x(1, catom)) < tol) then
                   if (abs(lammps_coord(2,iatom)-x(2, catom)) < tol) then
                      cadd_lammps_map(iatom) = catom
                      lammps_cadd_map(catom) = iatom
                   end if
                end if
	    end do
!!$          print *, 'Done mapping atom ', iatom, cadd_lammps_map(iatom)
	end if
       end do
    end if

    
    call lammps_gather_atoms(lmp,'x',3, r)
    nsize = size(r)/3
    allocate(rcoords(3,nsize))
    rcoords = reshape(r,shape(rcoords))

    ! Allocate and initialize arrays
    allocate(lammps_cadd_gmap(numnp))
    allocate(cadd_lammps_gmap(nsize))
    
    cadd_lammps_gmap = -HUGE(kind4int)
    lammps_cadd_gmap = -HUGE(kind4int)

!!$    print *, nsize, n_lammps_atoms
    
    do iAtom = 1, numnp
       if (isRelaxed(iAtom) /= 0) then
          do catom = 1, nsize
             if (abs(x(1,iAtom) - rcoords(1,catom)) < tol) then
                if (abs(x(2,iAtom) - rcoords(2,catom)) < tol) then
                   lammps_cadd_gmap(iAtom) = catom
                   cadd_lammps_gmap(catom) = iAtom
                end if
             end if
          end do
       end if
    end do
    if (allocated(r)) then
       deallocate(r)
    end if
    if (allocated(rcoords)) then
       deallocate(rcoords)
    end if
    print *, 'Finished mapping lammps to CADD atoms'
  end subroutine map_from_lammps

!!$ We need another subroutine to change positions in lammps
!!$    It is not wise manipulating the pointer to the main array directly i.e extract_atom
!!$     So gather_atoms and scatter_atoms is used.
!!$     Since gather atoms is ordered by atomid, we need to manipulate it accordingly
!!$     As with the previous mapping this needs to be called only once in the beginning
!!$     TODO check to make sure this subroutine is called before anything

   subroutine update_lammps_coords(atomcoord, atomdispl, update_pad, update_all, lmp)
    implicit none
    !<--- Input parameters ----
    double precision :: atomcoord(nxdm,numnp), atomdispl(nxdm, numnp)
    type(C_ptr) :: lmp
    !< -------------------------

    ! ---- Local Variables 
    character (len=1024) :: command_line

    integer :: error, narg, natoms

    double precision :: tol, max_coordx, min_coordx, max_coordy, min_coordy
    integer :: iatom, catom, lmpatom, j, minlmpatom, maxlmpatom
    double precision, dimension(:), allocatable :: r
    double precision, dimension(:,:), allocatable :: rcoords
    double precision, dimension(:), allocatable :: types

    real (C_double), pointer :: xlo => NULL()
    real (C_double), pointer :: xhi => NULL()
    real (C_double), pointer :: ylo => NULL()
    real (C_double), pointer :: yhi => NULL()

    character(len=10) :: atom_type

    logical :: update_all, update_pad, update_def
    
    call lammps_extract_global(xlo, lmp, 'boxxlo')
    call lammps_extract_global(xhi, lmp, 'boxxhi')
    call lammps_extract_global(ylo, lmp, 'boxylo')
    call lammps_extract_global(yhi, lmp, 'boxyhi')

!!$    print *, 'Lammps box size = ', xlo, xhi, ylo, yhi

    
    call lammps_gather_atoms(lmp,'x',3, r)
    natoms = size(r)/3
    allocate(rcoords(3,natoms))
    rcoords = reshape(r,shape(rcoords))

    ! --- collecting the types of each atom in lammps
    ! --- so we can ignore the indenter atoms
    !call lammps_gather_atoms(lmp,'type',1,types)

!!$    print*,'--------lammps_cadd_gmap--------'
!!$    print*,'size lammps_cadd_gmap', size(lammps_cadd_gmap)

    do iatom = 1, numnp
       if (update_all) then
          update_def = (isrelaxed(iAtom) /= 0) 
       else
          update_def = (isrelaxed(iAtom) == -1)
       end if

       if (update_def) then 
          lmpatom = lammps_cadd_gmap(iatom)
          rcoords(1,lmpatom) = atomcoord(1,iatom) + atomdispl(1,iatom)
          rcoords(2,lmpatom) = atomcoord(2,iatom) + atomdispl(2,iatom)
       end if
    end do

    do iatom = 1, natoms
       do j = 1,3
          r((iatom-1)*3 + j) = rcoords(j,iatom)
       end do
    end do
    
    call lammps_scatter_atoms(lmp,'x', r)

    return
    
    if (allocated(r)) then
       deallocate(r)
    end if
    if (allocated(rcoords)) then
       deallocate(rcoords)
    end if

    
  end subroutine update_lammps_coords

  subroutine update_from_lammps(AtomDispl, AtomCoord, AtomForce,  AveDispl, Velocity, Virst, AveVirst, lmp)
!!$ Updates local CADD arrays with the coordinates from lammps
!!$ The procedure is as follows
!!$     Lammps does not report displacement, but reports the actual current atom position
!!$     atomcoord contains the original position of the atom
!!$     atomdispl(i) = lammps_coord(
    USE MOD_OUTPUT
    USE MOD_MATERIAL
    USE MOD_PARALLEL
    USE MOD_TIMMING
    USE MOD_COMMON
    USE MOD_FILE
    USE MOD_DISL_PARAMETERS
    implicit none
    ! ----- Input variables -------
    DOUBLE PRECISION :: Atomdispl(NDF,NUMnp), Atomforce(3,*) , Atomcoord(NDF,*) 
    DOUBLE PRECISION :: AveDispl(3,*) , Velocity(3,*)
    DOUBLE PRECISION :: Virst(3,3,*), AveVirst(3,3,*)
    type(c_ptr)::lmp
    ! -----------------------------------


    ! ----- Averaged Displacement Components
    real (C_double), dimension(:), pointer :: compute_lammps_avg_dx => NULL()
    real (C_double), dimension(:), pointer :: compute_lammps_avg_dy => NULL()

  ! --------- Averaged Stress components -------
  real (C_double), dimension(:), pointer :: compute_lammps_avg_stress_xx => NULL()
  real (C_double), dimension(:), pointer :: compute_lammps_avg_stress_yy => NULL()
  real (C_double), dimension(:), pointer :: compute_lammps_avg_stress_zz => NULL()
  real (C_double), dimension(:), pointer :: compute_lammps_avg_stress_xy => NULL()
  real (C_double), dimension(:), pointer :: compute_lammps_avg_stress_zx => NULL()
  real (C_double), dimension(:), pointer :: compute_lammps_avg_stress_yz => NULL()

    
    
    integer :: iatom, i, j, nsize, lmpatom

    call lammps_extract_atom(lammps_force, lmp, 'f')
    
    call lammps_extract_atom(lammps_velocity, lmp, 'v')
    
    call lammps_extract_compute(compute_lammps_dx, lmp, 'dx_free', peratom_style, array_type) 

    call lammps_extract_fix(compute_lammps_avg_dx, lmp, 'dx_ave', peratom_style, vector_type, 1, 1)
    call lammps_extract_fix(compute_lammps_avg_dy, lmp, 'dy_ave', peratom_style, vector_type, 1, 1)
   
    call lammps_extract_compute(compute_lammps_stress, lmp, 'compute_stress', peratom_style, array_type)

    call lammps_extract_fix(compute_lammps_avg_stress_xx, lmp, 'stress_ave_xx', peratom_style, vector_type, 1, 1)
    call lammps_extract_fix(compute_lammps_avg_stress_yy, lmp, 'stress_ave_yy', peratom_style, vector_type, 1, 1)
    call lammps_extract_fix(compute_lammps_avg_stress_zz, lmp, 'stress_ave_zz', peratom_style, vector_type, 1, 1)
    call lammps_extract_fix(compute_lammps_avg_stress_xy, lmp, 'stress_ave_xy', peratom_style, vector_type, 1, 1)
    call lammps_extract_fix(compute_lammps_avg_stress_zx, lmp, 'stress_ave_zx', peratom_style, vector_type, 1, 1)
    call lammps_extract_fix(compute_lammps_avg_stress_yz, lmp, 'stress_ave_yz', peratom_style, vector_type, 1, 1)

    do iatom = 1, numnp
       if (isRelaxed(iatom) /= 0) then
	        AveDispl(1:2, iAtom) = 0.0d0
          if (isRelaxed(iatom) /= -1) then

             lmpatom = lammps_cadd_map(iatom)
			 
             Atomforce(1:2,iatom) = lammps_force(1:2,lmpatom)
             
             Velocity(1:2,iatom) = lammps_velocity(1:2,lmpatom)

             AtomDispl(1:2,iatom) = compute_lammps_dx(1:2, lmpatom)
             
             if (isRelaxed(iAtom) == 2) then 
                AveDispl(1, iatom) = compute_lammps_avg_dx(lmpatom)
                AveDispl(2, iatom) = compute_lammps_avg_dy(lmpatom)
             else
                AveDispl(1:2, iatom) = AtomDispl(1:2,iatom)

             end if

             do i = 1, 3
                do j = 1, 3
                   if (i == j) then 
                      Virst(i, i, iAtom) = compute_lammps_stress(i, lmpatom)
                   else
                      Virst(i, j, iAtom) = compute_lammps_stress(i+j+1, lmpatom)
                   end if
                end do
             end do

             AveVirst(1,1,iAtom) = compute_lammps_avg_stress_xx(lmpatom)
             AveVirst(2,2,iAtom) = compute_lammps_avg_stress_yy(lmpatom)
             AveVirst(3,3,iAtom) = compute_lammps_avg_stress_zz(lmpatom)
             AveVirst(1,2,iAtom) = compute_lammps_avg_stress_xy(lmpatom)
             AveVirst(1,3,iAtom) = compute_lammps_avg_stress_zx(lmpatom)
             AveVirst(2,3,iAtom) = compute_lammps_avg_stress_yz(lmpatom)
!!$             ! ---- Symmetric Stress Tensor
             AveVirst(2,1,iAtom) = AveVirst(1,2,iAtom)
             AveVirst(3,1,iAtom) = AveVirst(1,3,iAtom)
             AveVirst(3,2,iAtom) = AveVirst(2,3,iAtom)

             !end if
          end if
       end if
    end do

!!$    AveVirst(:,:,:) = Virst(:,:,:)

  end subroutine update_from_lammps
end module mod_lammps
