subroutine write_lammps_data(Id, X, Ix, F, B, Itx, xmin, xmax, ymin, ymax)
      !!> Subroutine to write initial atomic config into lammps format for lammps reading
      !!>      Identifies atom types, free atoms --> type 1, pad atoms --> 2
      !!> @TODO (Srinath#1#): Check to ensure data format can be read by lammps ...
      !!>   1) write the data file
      !!>    2) Check with independent lammps script.
      !!>    3) Include type for interface atoms
      !!>    4) All atom types from CADD are identified by a max number and the indenter etc will
      !!>        other atom types.

      USE MOD_GRAIN
      USE MOD_GLOBAL
      USE MOD_FILE
      USE MOD_BOUNDARY
      USE MOD_CRACK
      USE MOD_MATERIAL
      USE MOD_DD_SLIP
      USE MOD_DISL_PARAMETERS
      implicit none
                                !     Input variables
      INTEGER Id , Ix , Itx(3,*)
      DOUBLE PRECISION X , F , B
      DIMENSION Id(NDF,1) , X(NXDm,1) , Ix(NEN1,1) , F(NDF,1) , B(NDF,1)
      double precision, intent(in) :: xmin, xmax, ymin, ymax
      double precision :: zmin, zmax, yparticle
      double precision :: atom1_xmin, atom1_xmax, atom1_ymin, atom1_ymax
      integer :: nsteps
      integer :: atomType
      logical :: top, bot, left, right, itop, ibot, ileft, iright

      integer :: i, j, k, l, natoms, npad, n

      OPEN (UNIT=200,FILE='md.inp',STATUS='old')
      READ (200,*) stadium_width, exclude_top, exclude_bot, exclude_left, exclude_right
      READ (200,*) damp_coeff , damp_ratio
      READ (200,*) lammps_temperature
      READ (200,*) lammps_timestep
      READ (200,*) fem_update_steps
      READ (200,*) num_md_steps
      READ (200,*) lammps_output_steps
      READ (200,*) num_restart_steps
      READ (200,*) num_initial_equil
      READ (200,*) particle_velocity
      READ (200,*) particle_radius
      READ (200,*) particle_height
      READ (200,*) particle_rotation

!       write(*,'(A)') 'Lammps md paramters'
!       write(*,'(A, F15.6, 4(1X,I1))') 'Stadium parameter = ', stadium_width, exclude_top, exclude_bot, exclude_left, exclude_right
!       write(*,'(A, 2E15.6)') 'Damping Coefficeitn =', damp_coeff, damp_ratio
!       write(*,'(A, F15.3)') 'Temperature = ', lammps_temperature
!       write(*,'(A, F15.3)') 'Time Step = ', lammps_timestep/1.e-12
!       write(*,'(A, 5I7)') 'Steps (Update, mdsteps, output, restart, initial) = ', fem_update_steps,  & 
! 		num_md_steps, lammps_output_steps, num_restart_steps, num_initial_equil
!       write(*, '(A, 3F15.3)') 'Particle paramters (velocity, radius, height) = ', &
! 		particle_velocity, particle_radius, particle_height
      
      
      
      CLOSE (200)
      top = .true.
      bot = .true.
      left = .true.
      right = .true. 

      if (exclude_top > 0) then
         top = .false.
      end if

      if (exclude_bot > 0) then
         bot = .false.
      end if

      if (exclude_left > 0) then
         left = .false.
      end if

      if (exclude_right > 0) then
         right = .false.
      end if


      !! --- Assume units metal in lammps
      lammps_timestep = lammps_timestep/1.0d-12
      tstart = lammps_temperature
      tstop = lammps_temperature

      !!> @TODO add exclusion zones for stadium thermostat such as top surface, bottom surface, left, right
      !! --- mainly for impact problem

      if (damp_coeff > 0.0) then 
         damp_coeff = 1.0/(damp_coeff)
      ! TODO error handler 
      end if
      !! Once again assuming metal units in lammps
      damp_coeff = damp_coeff / 1.0d-12

      !!> TODO get zmin and zmax from mat file automatically
      zmin = 0.0d0
      zmax = 2.9573845299d0

      !!		JM: ***Must be careful to allow enough y space to add particle***
      !!		JM: ***over free surface***        
      yparticle = -ymin

      open(unit=1010, file='cadd_atoms.dat', status='UNKNOWN')
      natoms = 0
      atom1_xmin = 1.0d20
      atom1_xmax = -1.0d20
      atom1_ymin = 1.0d20
      atom1_ymax = -1.0d20

      do i = 1, numnp
         if (isRelaxed(i) /= 0) then
            natoms = natoms + 1
            if (isRelaxed(i) /= -1) then 
               if (X(1,i) < atom1_xmin) then 
                  atom1_xmin = X(1,i)
               end if
               if (X(2,i) < atom1_ymin) then 
                  atom1_ymin = X(2,i)
               end if
               if (X(1,i) > atom1_xmax) then 
                  atom1_xmax = X(1,i)
               end if
               if (X(2,i) > atom1_ymax) then 
                  atom1_ymax = X(2,i)
               end if
            end if
         end if
         if (isRelaxed(i) == -1) then
            npad = npad + 1
         end if
      end do
      stadium_xmax = atom1_xmax
      stadium_xmin = atom1_xmin
      stadium_ymin = atom1_ymin
      stadium_ymax = atom1_ymax


      if (particle_radius > abs(atom1_xmax - atom1_xmin)/2.0) then
         call error_handler('Particle radius must be smaller than substrate dimensions')
      end if

      if (particle_radius > abs(atom1_ymax - atom1_ymin)/2.0) then
         call error_handler('Particle radius must be smaller than substrate dimensions')
      end if


      write(1010,*) "CADD input atoms"
      write(1010,*)
      write(1010, fmt='(I7,1X,A10)')  natoms, 'atoms'

!!$		JM Changed number of atoms types to 5 to account for particle
!!$		if impact is being done
      write(1010, fmt='(I3,1X,A15)')  5, 'atom types'
      write(1010, fmt='(2(1X,F15.8),1X,A15)')  xmin, xmax, 'xlo xhi '
      write(1010, fmt='(2(1X,F15.8),1X,A15)')  ymin, yparticle, 'ylo yhi '
      write(1010, fmt='(2(1X,F15.8),1X,A15)')  zmin, zmax, 'zlo zhi '
      write(1010, *)
      write(1010,*) 'Atoms'
      write(1010,*)
      n = 0
      do i = 1, numnp
         if (isRelaxed(i) /=0 ) then
            n = n + 1
            if (isRelaxed(i) == -1) then
               atomType = 2
            elseif (isRelaxed(i) == 2) then
               atomType = 3
            else
               itop = (x(2,i) > atom1_ymax - stadium_width)
               ibot = (x(2,i) < atom1_ymin + stadium_width)
               ileft = (x(1,i) < atom1_xmin + stadium_width)
               iright = (x(1,i) > atom1_xmax - stadium_width)
               if ((itop .and. top) .or. (ibot .and. bot) .or. (ileft .and. left) .or. (iright .and. right)) then 
                  atomType = 4
               else 
                  atomType = 1
               end if
            end if
            write(1010,fmt='(I7,1X,I3,1X,3(1X,F15.8))') n, atomType, X(1,i), X(2,i), 0.0
         end if
      end do
      stadium_ymax = 2.0d0*stadium_width
      write(1010,*)
      close(1010)
end subroutine write_lammps_data


subroutine initialize_lammps(Id,X,Ix,F,B,Itx,lmp)
      !!> Subroutine to initialize lammps
      !!>  This performs the following functions
      !!>     Initialize lammps pointer
      !!>     Initialize potential for use
      !!>     Reads the initial mesh from mesh.f90 through a data file
      !!> @TODO (Srinath#1#): 1) Make data file name avaialble to this routine  ...
      !!>     2) Fixed atoms identified with atom type
      !!>     3) Compare initial configs between this and CADD

      use lammps
      USE MOD_GRAIN
      USE MOD_GLOBAL
      USE MOD_FILE
      USE MOD_BOUNDARY
      USE MOD_CRACK
      USE MOD_MATERIAL
      USE MOD_DD_SLIP
      USE MOD_DISL_PARAMETERS

      implicit none

      type (c_ptr) :: lmp
!!$     Input variables
      INTEGER Id , Ix , Itx(3,*)
      DOUBLE PRECISION X , F , B
      DIMENSION Id(NDF,1) , X(NXDm,1) , Ix(NEN1,1) , F(NDF,1) , B(NDF,1)
!!$        double precision, intent(in) :: xmin, xmax, ymin, ymax, padwidth

      double precision :: zmin, zmax, theta
      character(1024):: command_line
      integer :: iatom, i,j,k,l
      integer :: num_threads
      integer :: omp_get_num_threads

      real (C_double), pointer :: xlo => NULL()
      real (C_double), pointer :: xhi => NULL()
      real (C_double), pointer :: ylo => NULL()
      real (C_double), pointer :: yhi => NULL()

      !!JM: Changed z-bounds from before
      zmin = 0.0d0
      zmax = 2.9573845299d0
      write(*,'(A)') 'Lammps md paramters'
      write(*,'(A, F15.6, 4(1X,I1))') 'Stadium parameter = ', stadium_width, exclude_top, exclude_bot, exclude_left, exclude_right
      write(*,'(A, 2E15.6)') 'Damping Coefficeitn =', damp_coeff, damp_ratio
      write(*,'(A, F15.3)') 'Temperature = ', lammps_temperature
      write(*,'(A, F15.3)') 'Time Step = ', lammps_timestep/1.e-12
      write(*,'(A, 5I7)') 'Steps (Update, mdsteps, output, restart, initial) = ', fem_update_steps,  & 
		num_md_steps, lammps_output_steps, num_restart_steps, num_initial_equil
      write(*, '(A, 3F15.3)') 'Particle paramters (velocity, radius, height, rotation) = ', &
		particle_velocity, particle_radius, particle_height, particle_rotation

!!$        call lammps_open_no_mpi('lmg -log log.CADD', lmp)
!!$     If replacing this entire subroutine by reading input file
!!$        delete all lines in this file and replace by a lammps input file

!!$     call lammps_file(lmp, 'filename')

      call lammps_command(lmp, 'units metal')
      call lammps_command(lmp, 'atom_style atomic')
      call lammps_command(lmp, 'dimension 3')
#ifdef _OPENMP
!$OMP PARALLEL SHARED(num_threads)
      num_threads = omp_get_num_threads()
!$OMP END PARALLEL
      write(command_line,'(A,I3)') 'package omp ', num_threads
      call lammps_command(lmp, command_line)
      call lammps_command(lmp, 'suffix omp')
#endif
      !!		JM: Need a minimum fix distance in upper y to place particle in
      !!		JM: otherwise lammps shirnk wraps to the free surface and one cannot
      !!		JM: place a particle in the system
      call lammps_command(lmp, 'boundary ss sm pp')
      call lammps_command(lmp,'atom_modify sort 0 0.0 map array')

      !!JM: hcp lattice with short z cylinder region length gives a planar hex lattice
      !!JM: with a little hack using the burger's vector as the lattice constant a
      call lammps_command(lmp, 'lattice hcp 2.83087')

      call lammps_command(lmp, 'read_data cadd_atoms.dat')

      !!JM *** Create particle above free surface ***
      !!JM *** must be within simulation box as defined by yparticle ***
      !!JM *** in write_lammps_data subroutine above ***
      write(command_line,'(A,2F15.3,A)') 'region 1 cylinder z 0.0 ', &
           particle_height + particle_radius, particle_radius, ' 0.0 0.1 units box'
      call lammps_command(lmp, command_line)
      
!!$      call lammps_command(lmp, 'region 1 cylinder z 0.0 120.0 100.0 0.0 0.1 units box')

      !! --- creating particle atoms, type 5---
      call lammps_command(lmp, 'create_atoms 5 region 1')   
      call lammps_command(lmp, 'group particle_atoms type 5')
      
      !! --- add in for particle rotation by theta degrees ---
      theta = particle_rotation
      write(command_line,'(A,1F9.3,A,1F9.3,A)') 'displace_atoms particle_atoms rotate 0.0 ', particle_height + particle_radius, ' 0.0 0 0 1 ', theta, ' units box'
      call lammps_command(lmp, command_line)      
      
      !! --- Create groups of atoms for fixes and computes ----
      
      !!call lammps_command(lmp, "group md_atoms type 1 3 4")

      call lammps_command(lmp, "group md_atoms type 1 3 4 5")
      call lammps_command(lmp, "group free_atoms type 1")
      call lammps_command(lmp, "group pad_atoms type 2")
      call lammps_command(lmp, "group interface_atoms type 3")
      call lammps_command(lmp, "group langevin_atoms type 3 4")

      !! ------- EAM potentials
      call lammps_command(lmp, "pair_style eam/alloy")
      call lammps_command(lmp, "pair_coeff * * Al_adams_hex.eam.alloy Al Al Al Al Al")

      call lammps_command(lmp, "neighbor 0.1 bin ")
      call lammps_command(lmp, "neigh_modify delay 0 every 1 check yes")

                                ! ---------- Various Fixes ----------------------------------------------
      write(command_line,*) "variable mytemp equal", lammps_temperature
      call lammps_command(lmp, command_line)
      call lammps_command(lmp, "velocity md_atoms create $(2.0*v_mytemp) 426789 dist uniform")
      call lammps_command(lmp, "velocity md_atoms set NULL NULL 0.0 units box")

      call lammps_command(lmp, "fix fix_integ md_atoms nve")
      write(command_line, fmt='(A38,3(1X,F15.6),I7, A10, 5(1X,F15.6))') "fix fix_temp langevin_atoms langevin ", &
           tstart, tstop, damp_coeff, 699483, " stadium ", stadium_xmin, stadium_xmax, stadium_ymin, stadium_ymax, stadium_width
      call lammps_command(lmp, command_line)

      !! ----------------------------------------------------------------------------------

      !! ---- Temperature Computes and temperature variance for testing  -------------------------------------
      !! ---- This is a 2d Problem so the temperature compute is restricted to partial in the xy plane
      
      call lammps_command(lmp, "compute free_temp free_atoms temp/partial 1 1 0")
      call lammps_command(lmp, "compute stadium_temp langevin_atoms temp/partial 1 1 0")

      !!   --- Variables for actual temperature 
      call lammps_command(lmp, "variable free_tempv equal c_free_temp")
      call lammps_command(lmp, "variable stadium_tempv equal c_stadium_temp")

      !! ---- Canonical Temperature variation 
      call lammps_command(lmp, "variable nfree_atoms equal count(free_atoms)")
      call lammps_command(lmp, "variable nstadium_atoms equal count(langevin_atoms)")

      call lammps_command(lmp, "variable canonical_free equal $((v_mytemp^2)/v_nfree_atoms)")
      call lammps_command(lmp, "variable canonical_stadium equal $((v_mytemp^2)/v_nstadium_atoms)")
      call lammps_command(lmp, "variable delt_free equal (c_free_temp-v_mytemp)^2")
      call lammps_command(lmp, "variable delt_stadium equal (c_stadium_temp-v_mytemp)^2")

      !! ------------ Energies -------------------------
      call lammps_command(lmp, "compute com_pe free_atoms pe/atom")
      call lammps_command(lmp, "compute pe free_atoms reduce sum c_com_pe")
      call lammps_command(lmp, "compute ke free_atoms ke")
      call lammps_command(lmp, "variable tot_energy equal c_pe+c_ke")
      !! ------------------------------------------------------------------------

      !!---- Pad atoms always have zero force so this is fixed here to 0 
      call lammps_command(lmp, "fix fix_zeroforce pad_atoms setforce 0.0 0.0 0.0")
      call lammps_command(lmp, "fix fix_2d all setforce NULL NULL 0.0")


      !! --------- Compute differential displacement from original position
      call lammps_command(lmp, "compute dx_free md_atoms displace/atom")

      call lammps_command(lmp, "compute dx_all all displace/atom")

      !! ---- Compute used for average displacement of interface atoms  -----
      call lammps_command(lmp, "compute dx_inter interface_atoms displace/atom")


      !! ----- Now define a fix to actually calculate the average for interface atoms
      write(command_line, '(A38, 2(1X,I3), A15)') "fix dx_ave interface_atoms ave/atom 1 ", fem_update_steps, fem_update_steps, " c_dx_inter[1]"
      call lammps_command(lmp, command_line)
      write(command_line,  '(A38, 2(1X,I3), A15)') "fix dy_ave interface_atoms ave/atom 1 ", fem_update_steps, fem_update_steps, " c_dx_inter[2]"
      call lammps_command(lmp, command_line)


      !! ---- Compute used for virial stress on atoms
      call lammps_command(lmp, "compute compute_stress md_atoms stress/atom NULL virial")


      !! ----- Now define a fix to actually calculate the stress average for all atoms
      call lammps_command(lmp, "fix stress_ave_xx md_atoms ave/atom 1 25 25 c_compute_stress[1]")
      call lammps_command(lmp, "fix stress_ave_yy md_atoms ave/atom 1 25 25 c_compute_stress[2]")
      call lammps_command(lmp, "fix stress_ave_zz md_atoms ave/atom 1 25 25 c_compute_stress[3]")
      call lammps_command(lmp, "fix stress_ave_xy md_atoms ave/atom 1 25 25 c_compute_stress[4]")
      call lammps_command(lmp, "fix stress_ave_zx md_atoms ave/atom 1 25 25 c_compute_stress[5]")
      call lammps_command(lmp, "fix stress_ave_yz md_atoms ave/atom 1 25 25 c_compute_stress[6]")

      call lammps_extract_global(xlo, lmp, 'boxxlo')
      call lammps_extract_global(xhi, lmp, 'boxxhi')
      call lammps_extract_global(ylo, lmp, 'boxylo')
      call lammps_extract_global(yhi, lmp, 'boxyhi')


      !! ---- Dump data file 
      write(command_line, '(A18,I3,A186)') "dump 1 all custom ", lammps_output_steps, " atom_lmp*.cfg id type x y z c_dx_all[1] c_dx_all[2] vx vy vz c_dx_all[4] f_fix_temp"
		
      call lammps_command(lmp, command_line)
!!$      write(command_line, '(A22,I3,A81)') "dump 2 all custom/vtk ", lammps_output_steps, " atom_lmp*.vtu id type x y z c_dx_all[1] c_dx_all[2] vx vy vz c_dx_all[4]"
!!$      call lammps_command(lmp, command_line)
      
      !!let lammps equilibrate for 10 picoseconds, then call indenter
      !!place this call in another subroutine called equilibrate lammps
      call lammps_command(lmp, "run 0")

      call lammps_command(lmp, "run 0")
end subroutine initialize_lammps


subroutine equilibrate_lammps(lmp, equil_steps)
      !!> Subroutine to let atoms equilibrate to desired temperature
      !!> in our NVE ensemble in lammps (init to 2x desired temp)
      !!> CADD-lammps mapping is done via positional tolerances
      !!> that get screwed up if an initial run # is given in
      !!> initialize lammps

      use lammps

      implicit none
      integer :: equil_steps
      type (c_ptr) :: lmp

      character(1024):: command_line

      write(command_line, '(A,I7)') 'run ', equil_steps
!!$     equilibrate lammps for 10 picoseconds
      call lammps_command(lmp,'run 0')

end subroutine equilibrate_lammps

subroutine add_fix_lammps(lmp, velocity)
      !!> Subroutine to add fixes to lammps post initialization
      !!> For now this will add a particle velocity for cold-spray
      !!> impact studies post equilibrate_lammps

      use lammps

      implicit none

      double precision :: velocity
      type (c_ptr) :: lmp

      character(1024):: command_line

!!$     assign particle a velocity after letting it and the substrate equilibrate

!!$         100 m/s in A/ps
!!$        call lammps_command(lmp, "velocity particle_atoms set NULL -1.0 NULL sum yes units box")

!!$     500 m/s for 15 um radius particle...equivalent velocity for 100 nm particle
!!$        call lammps_command(lmp, "velocity particle_atoms set NULL -1842.98 NULL sum yes units box")

!!$     using KSR: 500 m/s for 15 um radius particle...equivalent velocity for 100 nm particle
      write(command_line, '(A,F15.6,A)') 'velocity particle_atoms set NULL ', velocity, ' NULL sum yes units box'
      call lammps_command(lmp, command_line)
!!$      call lammps_command(lmp, "velocity particle_atoms set NULL -50.0 NULL sum yes units box")


end subroutine add_fix_lammps
