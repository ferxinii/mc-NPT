! ****************************************************************************
!          Molecular Dynamics code to simulate at NVP collectivity 
!          a system of N atoms of Ar inside a cubic box, interacting
!          through a truncated Lennard-Jones force field at 2.5 sigma. 
!          Leap-frog Verlet integration algorithm.
!****************************************************************************
PROGRAM mc_npt
    IMPLICIT NONE
    ! PARAMETERS
    INTEGER, PARAMETER          :: nbins_gdr = 500
    DOUBLE PRECISION, PARAMETER :: r_crit = 2.5  ! Cutoff for LJ interactions
    DOUBLE PRECISION, PARAMETER :: r_crit_gdr = 5
    DOUBLE PRECISION, PARAMETER :: dr_gdr = r_crit_gdr / nbins_gdr
    DOUBLE PRECISION, PARAMETER :: PI = 4.d0*ATAN(1.d0)
    !
    INTEGER :: ii
    ! CONFIGURATION
    INTEGER :: n_steps, freq_sample, n_atoms
    DOUBLE PRECISION :: mass, t_ref, sigma, eps, dr, dv, p_ref, l, beta
    ! 
    DOUBLE PRECISION, ALLOCATABLE :: r(:, :)
    DOUBLE PRECISION, ALLOCATABLE :: gdr(:)
    !
    DOUBLE PRECISION :: Ptail, Ptail_rho2, Utail, Utail_rho
    DOUBLE PRECISION :: denom_nid
    !
    !
    INTEGER :: sample_counter = 0, attempts_r = 0, attempts_v = 0
    INTEGER :: accepted_moves_v = 0, accepted_moves_r = 0
    INTEGER :: rand_id, sample_counter
    DOUBLE PRECISION :: Ptail, Utail

    CALL RANDOM_SEED()

    CALL read_config(n_steps, freq_sample, n_atoms, mass, t_ref, &
                     sigma, eps, dr, dv, l, p_ref)
    ALLOCATE(r(3, n_atoms))

    CALL read_initial_positions(r, n_atoms, l)

    CALL reduce_units(n_atoms, r, l, t_ref, eps, sigma)
    beta = 1/t_ref

    WRITE(*, *) "l: ", l
    WRITE(*, *) "t_ref: ", t_ref
    WRITE(*, *) "r(1,1): ", r(1,1)
     
    ALLOCATE(gdr(nbins_gdr))
    gdr = 0 

    Ptail_rho2 = 16.d0/3.d0*PI*(2.d0/3.d0*(1/r_crit)**9-(1/r_crit)**3)
    Utail_rho = 8.d0/3.d0*PI*n_atoms*(1.d0/3.d0*(1/r_crit)**9-(1/r_crit)**3)

    denom_nid = ((1.d0*n_steps/freq_sample)*n_atoms)

!-----5.1 MAIN LOOP TO GENERATE CONFIGURATIONS
    CALL initialise_output()
    
    accepted_moves_r = 0
    accepted_moves_v = 0
    sample_counter = 1
    attempts_r = 0
    attempts_v = 0
    DO ii = 1, n_steps
        !TODO 
        !Attempt to move npart before attempting to change V:
        CALL RANDOM_NUMBER(rand)
        rand_id =  rand * (n_atoms + 1) + 1
        IF (rand_id .LT. n_atoms) THEN
            CALL mc_pos(n_atoms, r, dr, beta, r_crit, L, &
                        accepted_moves_r)
            attempts_r = attempts_r + 1
        ELSE 
            CALL mc_vol(p_ref, n_atoms, r, L, dv, beta, r_crit, &
                        accepted_moves_v)
            attempts_v = attempts_v + 1
        ENDIF
      
        IF (MOD(ii, freq_sample) .EQ. 0) THEN
            rho = n_atoms / (L**3)
            nid = 4.0 / 3.0 * PI * rho
            CALL sample(r, n_atoms, L, r_crit, u, p_pot, gdr, rho, &
                        nid, denom_nid)
            !we compute the pressure from the virial to make sure its ok
            !(we need rho with new boxlength)
            
            Ptail = Ptail_rho2 * rho**2
            Utail = Utail_rho * rho
            WRITE(5, *) samp_count, p_pot, Ptail, rho*tref, &
                       p_pot+Ptail+rho*tref
            WRITE(16, *) samp_count, u+Utail
            WRITE(17, *) L
            
            sample_count = sample_count + 1
        ENDIF
    ENDDO

    close(5)
    close(16)
    close(17)
    

    ! Print and write acceptance rates
    write(*,*) "samples=", sample_count
    write(*,*) "fraction of accepted moves in r =", &
                (1.0*accepted_moves_r) / attempts_r
    write(*,*) "fraction of accepted moves in v =", &
                (1.0*accepted_moves_v) / attempts_v
    
    OPEN(12, file='acceptance_rate.dat')
    write(12,*) (1.0*accepted_moves_r) / attempts_r
    write(12,*) (1.0*accepted_moves_v) / attempts_v
    CLOSE(12)

END PROGRAM mc_npt

!-----5.2 Normalise gdr & save, simil for vacf and r2 
!     write(filename, "('gdr',i0,'.dat')") actual_it
!     open(12, file= trim(adjustl(filename)))
!     DO is = 1,nb_gdr
!       vb=((is+1)**3-is**3)*delg**3
!       r_gdr=delg*(is+1.0/2.0)
!       WRITE(12,*) r_gdr, gdr(is)/vb
!     ENDDO
!     CLOSE(12)


!-----6. Saving last configuration in A and A/ps 
!     open(11,file='conf_MC.data',status='unknown')
!     do ii = 1,natoms
!        write(11,*) (r(jj,ii)*sigma,jj=1,3)
!     end do         
!     write(11,*) boxlength*sigma
!     close(11)

!     stop
!     end program



SUBROUTINE read_config(n_steps, &
                       freq_sample, &
                       n_atoms, &
                       mass, &
                       t_ref, &
                       sigma, &
                       eps, &
                       dr, &
                       dv, &
                       l, &
                       p_ref)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: n_steps, freq_sample, n_atoms
    DOUBLE PRECISION, INTENT(OUT) :: mass, t_ref, sigma, eps, &
                                     dr, dv, l, p_ref
    CHARACTER(LEN = 256) :: dummy
    INTEGER ios
    OPEN(1, FILE='config.data', STATUS='old', IOSTAT=ios)
    IF (ios /= 0) THEN
        WRITE(*,*) "Error opening file, ios = ", ios
        STOP
    ENDIF
    READ(1, '(A)') dummy
    READ(1, *) dummy, n_steps
    READ(1, *) dummy, freq_sample
    READ(1, '(A)') dummy
    READ(1, *) dummy, mass
    READ(1, *) dummy, sigma
    READ(1, *) dummy, eps  ! Character '/' divides strings
    READ(1, '(A)') dummy
    READ(1, *) dummy, n_atoms
    READ(1, *) dummy, p_ref
    READ(1, *) dummy, t_ref
    READ(1, '(A)') dummy 
    READ(1, *) dummy, l 
    READ(1, *) dummy, dr
    READ(1, *) dummy, dv
    CLOSE(1)
    
    ! Change n_atoms to be a perfect cube if necessary
    IF (NINT(n_atoms**(1.0/3.0))**3 /= n_atoms) THEN
        WRITE(*, *) "ATTENTION! The number of atoms should be a &
                    &perfect cube..."
        WRITE(*, '(A, I0, A, I0)') " Changing n_atoms from ", n_atoms, &
                                   " to ", NINT(n_atoms**(1.0/3.0))**3
        WRITE(*, *)
        n_atoms = NINT(n_atoms**(1.0/3.0))**3
    ENDIF
END SUBROUTINE read_config



SUBROUTINE read_initial_positions(r, n_atoms, l)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n_atoms
    DOUBLE PRECISION, INTENT(IN) :: l
    DOUBLE PRECISION, INTENT(OUT) :: r(3, n_atoms)
    LOGICAL :: file_exists
    CHARACTER(LEN = 256) :: dummy
    INTEGER :: n_lines, ios, divisions
    INTEGER :: ii, jj, kk, ls
    DOUBLE PRECISION :: delta_r

    ! Check if a file with previous positions exists
    INQUIRE(FILE='initial_pos.dat', EXIST=file_exists)
    IF (file_exists) THEN
        OPEN(1, FILE='initial_pos.dat', STATUS='OLD')
        READ(1, '(A)')  ! Ignore first line
        n_lines = 0
        ii = 1
        DO 
            READ(1, *, IOSTAT=ios) (r(jj,ii), jj=1,3)
            ! Check for EOF or out of bounds
            IF (ios /= 0 .or. n_lines - 1 >= n_atoms) EXIT  
            n_lines = n_lines + 1
        ENDDO
        CLOSE(1)
        
        IF (n_lines /= n_atoms) THEN
            WRITE(*, '(A, I0, A, I0, A)') " Read ", n_lines, " initial &
                  &positions, but there are ", n_atoms, " atoms."
            WRITE(*, *) "Creating uniform distribution of initial positions"
            WRITE(*, *)
            file_exists = .FALSE.
        ENDIF
    ENDIF
    
    IF (.NOT. file_exists) THEN
        ! Uniform distribution of particles in the box
        divisions = INT(n_atoms**(1.0/3.0))
        delta_r = l/ divisions
        ls = 1
        DO ii = 0, divisions - 1
            DO jj = 0, divisions - 1
                DO kk = 0, divisions - 1
                    r(1, ls) = -l/2 + delta_r/2 + ii*delta_r
                    r(2, ls) = -l/2 + delta_r/2 + jj*delta_r
                    r(3, ls) = -l/2 + delta_r/2 + kk*delta_r
                    ls = ls + 1
               ENDDO
            ENDDO
        ENDDO
    
        !WRITE(*, *)
        !WRITE(*,*) l, divisions, delta_r
        !DO ii = 1, n_atoms
        !    WRITE(*, *) (r(jj, ii), jj=1,3)
        !END DO
    ENDIF

END SUBROUTINE read_initial_positions




SUBROUTINE initialise_output()
    IMPLICIT NONE
    CHARACTER(LEN = 256) :: filename
    INTEGER :: id
    
    OPEN(1, FILE="current_id.tmp", STATUS="OLD")
    READ(1, *) id
    CLOSE(1)

    WRITE(filename, "('pressure_',I0,'.dat')") id
    OPEN(11, FILE = TRIM(ADJUSTL(filename)))
    WRITE(filename, "('energy_',I0,'.dat')") id
    OPEN(12, FILE = TRIM(ADJUSTL(filename)))
    WRITE(filename, "('boxlength_',I0,'.dat')") id
    OPEN(13, FILE = TRIM(ADJUSTL(filename)))

END SUBROUTINE initialise_output




!********************************************************
!********************************************************
!              subroutine reduced
!********************************************************
!********************************************************

SUBROUTINE reduce_units(n_atoms, &
                              r, &
                              L, &
                              t_ref, &
                              eps, &
                              sigma)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n_atoms
    DOUBLE PRECISION, INTENT(IN) :: eps, sigma
    DOUBLE PRECISION, INTENT(INOUT) :: r(3, n_atoms), L, t_ref
    INTEGER :: ii, jj

    t_ref = t_ref / eps
    L = L / sigma
    DO ii = 1, n_atoms
       DO jj = 1, 3
          r(jj, ii) = r(jj, ii) / sigma 
       ENDDO
    ENDDO

END SUBROUTINE reduce_units

!********************************************************
!********************************************************
!              subroutine mcvol
!********************************************************
!********************************************************
!     SUBROUTINE mcvol(P, npart, r, L, delv, beta, rc,
!    &                 acc_moves_v)

!     implicit double precision (a-h,o-z)
!     dimension r(3,1000)
!     double precision L, LN, lnVN, VN
!     integer ii, jj, acc_moves_v
!     
!     call total_energy(r, rc, L, entot_old, npart)

!     V0 = L**3
!     !random walk on ln(V):
!     lnVN = log(V0) + (rand()-0.5)*delv
!     VN = exp(lnVN)
!     LN = VN**(1.0/3.0)
!     
!     !reescale particles:
!     do ii=1,npart 
!       do jj=1,3
!         r(jj,ii) = r(jj,ii)*LN/L
!       enddo
!     enddo

!     call total_energy(r, rc, LN, entot_new, npart)

!     arg = -beta*( entot_new - entot_old + P*(VN-V0) -
!    &      (npart+1)*log(VN/V0)/beta )

!     !If rejected, return particles to their original scaling+
!     if (rand() .gt. exp(arg)) then
!       do ii=1,npart
!         do jj=1,3
!           r(jj,ii) = r(jj,ii)*L/LN
!         enddo
!       enddo
!     else !if accepted, r already changed, update L and count move
!       acc_moves_v = acc_moves_v + 1
!       L=LN
!     endif
!     
!     RETURN
!     END



!********************************************************
!********************************************************
!              subroutine mcmove
!********************************************************
!********************************************************
!     SUBROUTINE mcmove(npart, r, delr, beta, rc, boxlength, 
!    &                  accepted_moves)
!       implicit double precision (a-h, o-z)
!       dimension r(3,1000), ri_new(3), ri_old(3)
!       integer npart, jj, o, accepted_moves

!       o = int(rand()*npart) + 1 !select random particle
!       
!       do jj = 1,3 
!         ri_old(jj) = r(jj,o)
!         ri_new(jj) = r(jj,o) + (rand()-0.5)*delr
!       enddo
!       
!       call energy(npart, o, ri_new, r, en_new, rc, boxlength)
!       call energy(npart, o, ri_old, r, en_old, rc, boxlength)

!       !Acceptance rule
!       if (rand() .lt. exp(-beta*(en_new-en_old)) ) then
!         do jj=1,3
!           r(jj,o) = ri_new(jj)
!         enddo
!         accepted_moves = accepted_moves + 1
!       endif

!     end
!********************************************************
!********************************************************
!              subroutine total energy
!********************************************************
!********************************************************
!     SUBROUTINE total_energy(r, rc, boxlength, entot, npart)
!     implicit double precision (a-h, o-z)
!     dimension r(3, 1000)
!     integer ii, jj

!     entot = 0;
!     do ii = 1,npart-1
!       do jj = ii+1, npart
!         call lj(ii,jj,r,boxlength,rc,pot,rij_fij)
!         entot = entot + pot
!       enddo
!     enddo
!     
!     RETURN
!     END

!********************************************************
!********************************************************
!              subroutine energy
!********************************************************
!********************************************************
!     SUBROUTINE energy(npart, o, ri, r, en, rc, boxlength)
!       implicit double precision (a-h, o-z)
!       dimension r(3,1000), ri(3)
!       integer npart, jj, o
!       
!       do jj = 1,3
!         r(jj,o)=ri(jj)
!       enddo

!       en = 0
!       do jj=1,npart
!         if (jj .ne. o) then
!           call lj(o,jj,r,boxlength,rc,pot,rij_fij)
!           en = en + pot
!         endif
!       enddo
!       
!     RETURN
!     END

!********************************************************
!********************************************************
!              subroutine Lennard-Jones
!********************************************************
!********************************************************

!     subroutine lj(is,js,r,boxlength,rc,pot,rij_fij)
!     implicit double precision(a-h,o-z)
!     dimension r(3,1000), rij(3)
!     integer is, js

!     rr2 = 0.d0
!     pot = 0.d0
!     rij_fij = 0.d0
!     do l = 1,3
!        rijl = r(l,js) - r(l,is)
!        rij(l) = rijl - boxlength*dnint(rijl/boxlength)
!        rr2 = rr2 + rij(l)*rij(l)
!     end do

!     rr = dsqrt(rr2)

!     if (rr.lt.rc) then
!       ynvrr2 = 1.d0/rr2
!       ynvrr6 = ynvrr2*ynvrr2*ynvrr2
!       ynvrr12 = ynvrr6*ynvrr6
!       rij_fij = rr2*24.d0*(2.d0*ynvrr12-ynvrr6)*ynvrr2
!       pot = 4.d0*(ynvrr12-ynvrr6) 
!     else
!       pot = 0
!       rij_fij = 0
!     endif
!     
!     return
!     end

!********************************************************
!********************************************************
!              subroutine sample
!********************************************************
!********************************************************
!     SUBROUTINE sample(r,npart,boxlength,rc,u,p_pot,gdr, rho, nid,
!    &                  denom_nid)
!       implicit double precision (a-h,o-z)
!       dimension r(3,1000), gdr(500)
!       integer ii, jj, npart, ncycle, nsamp
!       double precision nid
!       u = 0
!       p_pot = 0

!       do ii = 1, npart-1
!         do jj = ii+1, npart
!             call lj(ii,jj,r,boxlength,rc,pot,rij_fij)
!             u = u + pot
!             p_pot = p_pot + rij_fij
!             call function_gdr(ii,jj,r,boxlength,rc,gdr, rho, nid,
!    &                          denom_nid)
!         enddo
!       enddo
!       p_pot = p_pot/(3*boxlength**3)
!     END

!********************************************************
!********************************************************
!              subroutine gdr
!********************************************************
!********************************************************
!     SUBROUTINE function_gdr(is,js,r,boxlength,rc,gdr, rho, nid,
!    &                        denom_nid)
!       implicit double precision (a-h,o-z)
!       dimension r(3,1000), gdr(500), rij(3)
!       integer is, l, ncycle, nsamp, natoms
!       double precision nid

!       rr2 = 0.d0
!       pot = 0.d0
!       do l = 1,3
!         rijl = r(l,js) - r(l,is)
!         rij(l) = rijl - boxlength*dnint(rijl/boxlength)
!         rr2 = rr2 + rij(l)*rij(l)
!       end do

!       rr = dsqrt(rr2)
!     
!       nb_gdr=500
!       rc_gdr=5
!       
!       !we have to initialize all bins to 0...

!       delg=rc_gdr/nb_gdr

!       IF (rr.lt.rc_gdr) THEN
!         ni=nint(rr/delg)
!         gdr(ni)=gdr(ni)+2.d0/(denom_nid*nid)
!         !write(*,*) ni, gdr(ni)
!       END IF
!       
!     END
