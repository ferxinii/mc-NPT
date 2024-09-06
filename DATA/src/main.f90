! Monte Carlo simulation of an isobaric-isothermic (N, P, T) collectivity
! Updated in 30/08/2024 Fernando Muñoz

PROGRAM mc_npt
    IMPLICIT NONE
    ! PARAMETERS
    INTEGER, PARAMETER          :: nbins_gdr = 500
    DOUBLE PRECISION, PARAMETER :: r_crit = 2.5  ! Cutoff for LJ interactions
    DOUBLE PRECISION, PARAMETER :: r_crit_gdr = 5
    DOUBLE PRECISION, PARAMETER :: dr_gdr = r_crit_gdr / nbins_gdr
    DOUBLE PRECISION, PARAMETER :: PI = 4.d0*ATAN(1.d0)
    ! CONFIGURATION
    INTEGER :: n_steps, freq_sample, n_atoms
    DOUBLE PRECISION :: mass, t_ref, sigma, eps, dr, dv, p_ref, l, beta
    ! ARRAYS
    DOUBLE PRECISION, ALLOCATABLE :: r(:, :)
    DOUBLE PRECISION, ALLOCATABLE :: gdr(:)
    ! ENERGY AND PRESSURE
    DOUBLE PRECISION :: u_tail_rho, p_tail_rho2
    DOUBLE PRECISION :: u, p_pot, rho
    DOUBLE PRECISION :: denom_nid, nid
    ! COUNTERS
    INTEGER :: sample_counter = 1, attempts_r = 0, attempts_v = 0
    INTEGER :: accepted_moves_v = 0, accepted_moves_r = 0
    ! AUX
    INTEGER :: ii
    DOUBLE PRECISION :: rand
    INTEGER :: rand_id

    CALL RANDOM_INIT(REPEATABLE = .FALSE., IMAGE_DISTINCT=.TRUE.)

    CALL read_config(n_steps, freq_sample, n_atoms, mass, t_ref, &
                     sigma, eps, dr, dv, l, p_ref)
    ALLOCATE(r(3, n_atoms))
    CALL read_initial_positions(r, n_atoms, l)

    CALL reduce_units(n_atoms, r, l, t_ref, eps, sigma)
    beta = 1 / t_ref

    ALLOCATE(gdr(nbins_gdr))
    gdr = 0 

    p_tail_rho2 = 16.d0/3.d0*PI*(2.d0/3.d0*(1/r_crit)**9-(1/r_crit)**3)
    u_tail_rho = 8.d0/3.d0*PI*n_atoms*(1.d0/3.d0*(1/r_crit)**9-(1/r_crit)**3)

    ! Aux value used for calculation of gdr
    denom_nid = ((1.d0*n_steps / freq_sample) * n_atoms)  

    CALL initialise_output_samples()
    
    DO ii = 1, n_steps
        ! Attempt to move r of n_part (average) before attempting to change V
        CALL RANDOM_NUMBER(rand)
        rand_id =  INT(rand * (n_atoms + 1) + 1)
        IF (rand_id .LT. n_atoms) THEN
            CALL mc_pos(n_atoms, r, dr, beta, r_crit, L, &
                        accepted_moves_r)
            attempts_r = attempts_r + 1
        ELSE 
            CALL mc_vol(n_atoms, p_ref, beta, dv, r_crit, r, L, &
                        accepted_moves_v)
            attempts_v = attempts_v + 1
        ENDIF
      
        IF (MOD(ii, freq_sample) .EQ. 0) THEN
            rho = n_atoms / (L**3)
            nid = 4.0 / 3.0 * PI * rho  ! Aux value used for calculation of gdr

            CALL sample(n_atoms, r, L, r_crit, nid, denom_nid, &
                  r_crit_gdr, nbins_gdr, u, p_pot, gdr)
            
            CALL write_sample(sample_counter, p_tail_rho2, u_tail_rho, &
                              rho, p_pot, t_ref, u, L)
                        
            sample_counter = sample_counter + 1
        ENDIF
    ENDDO
    
    CALL terminate_output_samples()
    
    CALL write_acceptance_rates(sample_counter, accepted_moves_r, attempts_r, &
                                accepted_moves_v, attempts_v)
    
    CALL normalise_and_write_gdr(nbins_gdr, gdr, dr_gdr)
    
    CALL write_last_positions(n_atoms, r, sigma)

END PROGRAM mc_npt



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
    OPEN(1, FILE='config.dat', STATUS='old', IOSTAT=ios)
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
    INTEGER :: n_lines, ios, divisions
    INTEGER :: ii, jj, kk, ls
    DOUBLE PRECISION :: delta_r

    ! Check if a file with previous positions exists
    INQUIRE(FILE='results/initial_pos.dat', EXIST=file_exists)
    IF (file_exists) THEN
        OPEN(1, FILE='results/initial_pos.dat', STATUS='OLD')
        READ(1, '(A)')  ! Ignore first line
        n_lines = 1
        ii = 1
        DO 
            READ(1, *, IOSTAT=ios) (r(jj,ii), jj=1,3)
            ! Check for EOF or out of bounds
            IF (ios /= 0 .or. n_lines - 1 >= n_atoms) EXIT  
            n_lines = n_lines + 1
            ii = ii + 1
        ENDDO
        CLOSE(1)
        
        IF (n_lines /= n_atoms) THEN
            WRITE(*, '(A, I0, A, I0, A)') " Read ", n_lines, " initial &
                  &positions, but there are ", n_atoms, " atoms."
            WRITE(*, *)
            file_exists = .FALSE.
        ELSE 
            WRITE(*, '(A, I0, A)') "Succesfully read ", n_atoms, &
                                   " initial positions"
        ENDIF
    ENDIF
    
    IF (.NOT. file_exists) THEN
        WRITE(*, *) "Creating uniform distribution of initial positions"
        ! Uniform distribution of particles in the box
        divisions = INT(n_atoms**(1.0/3.0))
        delta_r = l / divisions
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




SUBROUTINE initialise_output_samples()
    IMPLICIT NONE
    CHARACTER(LEN = 256) :: filename
    INTEGER :: id
    
    OPEN(1, FILE="tmp/current_id.tmp", STATUS="OLD")
    READ(1, *) id
    CLOSE(1)

    WRITE(filename, "('results/pressure_',I0,'.dat')") id
    OPEN(11, FILE=TRIM(ADJUSTL(filename)))
    WRITE(filename, "('results/energy_',I0,'.dat')") id
    OPEN(12, FILE=TRIM(ADJUSTL(filename)))
    WRITE(filename, "('results/boxlength_',I0,'.dat')") id
    OPEN(13, FILE=TRIM(ADJUSTL(filename)))

END SUBROUTINE initialise_output_samples



SUBROUTINE terminate_output_samples()
    CLOSE(11)
    CLOSE(12)
    CLOSE(13)
END SUBROUTINE



SUBROUTINE write_sample(sample_counter, p_tail_rho2, u_tail_rho, rho, &
                        p_pot, t_ref, u, L)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: sample_counter
    DOUBLE PRECISION, INTENT(IN) :: p_tail_rho2, u_tail_rho, rho, p_pot, &
                                    t_ref, u, L
    DOUBLE PRECISION :: p_tail, u_tail

    p_tail = p_tail_rho2 * rho**2
    u_tail = u_tail_rho * rho
    
    ! Compute virial pressure to compare with imposed pressure
    WRITE(11, *) sample_counter, p_pot, p_tail, rho * t_ref, &
                 p_pot + p_tail + rho * t_ref
    WRITE(12, *) sample_counter, u + u_tail
    WRITE(13, *) sample_counter, L

END SUBROUTINE write_sample



SUBROUTINE write_acceptance_rates(sample_counter, accepted_moves_r, &
                                  attempts_r, accepted_moves_v, attempts_v)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: sample_counter, accepted_moves_r, attempts_r, &
                           accepted_moves_v, attempts_v

    WRITE(*, *) "........ RATES: ........"
    WRITE(*, '(A, I0)') " Total samples = ", sample_counter
    WRITE(*, *) "Fraction of accepted moves in r = ", &
                (1.0 * accepted_moves_r) / attempts_r
    WRITE(*, *) "Fraction of accepted moves in v = ", &
                (1.0 * accepted_moves_v) / attempts_v
    
    OPEN(14, FILE='tmp/acceptance_rate.tmp')
    WRITE(14, *) (1.0 * accepted_moves_r) / attempts_r
    WRITE(14, *) (1.0 * accepted_moves_v) / attempts_v
    CLOSE(14)

END SUBROUTINE write_acceptance_rates



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



SUBROUTINE mc_vol(n_atoms, p_ref, beta, dv, r_crit, r, L, accepted_moves_v)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n_atoms
    DOUBLE PRECISION, INTENT(IN) :: p_ref, beta, dv, r_crit
    DOUBLE PRECISION, INTENT(INOUT) :: r(3, n_atoms), L
    INTEGER, INTENT(INOUT) :: accepted_moves_v
    DOUBLE PRECISION :: Etot_old, Etot_new, v0, LN, ln_vN, vN, rand, arg
    INTEGER :: ii, jj
    
    CALL total_energy(n_atoms, r, r_crit, L, Etot_old)

    v0 = L**3
    ! Random walk on ln(V):
    CALL RANDOM_NUMBER(rand)
    ln_vN = LOG(v0) + (rand - 0.5) * dv
    vN = EXP(ln_vN)
    LN = vN**(1.0/3.0)
    
    ! Rescale particles:
    DO ii = 1, n_atoms 
         DO jj=1,3
              r(jj,ii) = r(jj,ii) * LN / L
         ENDDO
    ENDDO

    CALL total_energy(n_atoms, r, r_crit, LN, Etot_new)

    arg = - beta * ( Etot_new - Etot_old + p_ref * (vN - v0) - &
                     (n_atoms + 1) * log(vN/v0) / beta )

    ! If rejected, return particles to their original scaling+
    CALL RANDOM_NUMBER(rand)
    IF (rand .GT. exp(arg)) THEN
        DO ii=1, n_atoms
            DO jj=1,3
                r(jj,ii) = r(jj,ii) * L / LN
            ENDDO
        ENDDO
    ELSE ! If accepted, r already changed, update L and count move
        accepted_moves_v = accepted_moves_v + 1
        L = LN
    ENDIF

END



SUBROUTINE mc_pos(n_atoms, r, dr, beta, r_crit, L, accepted_moves_r)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n_atoms
      DOUBLE PRECISION, INTENT(INOUT) :: r(3, n_atoms)
      DOUBLE PRECISION, INTENT(IN) :: dr, beta, r_crit, L
      INTEGER, INTENT(INOUT) :: accepted_moves_r
      DOUBLE PRECISION :: rand, ri_new(3), ri_old(3), E_new, E_old
      INTEGER :: jj, o
      
      CALL RANDOM_NUMBER(rand)
      o = INT(rand * n_atoms) + 1  ! Select random particle
      
      DO jj = 1, 3 
          CALL RANDOM_NUMBER(rand)
          ri_old(jj) = r(jj, o)
          ri_new(jj) = r(jj, o) + (rand - 0.5) * dr
      ENDDO
      
      DO jj = 1, 3
          r(jj, o) = ri_new(jj)
      ENDDO
      CALL energy_ii(n_atoms, o, r, r_crit, L, E_new)
      
      DO jj = 1, 3
          r(jj, o) = ri_old(jj)  ! We return r to its old state!
      ENDDO
      CALL energy_ii(n_atoms, o, r, r_crit, L, E_old)
      
      ! Acceptance rule
      CALL RANDOM_NUMBER(rand)
      IF (rand .LT. EXP(-beta * (E_new - E_old))) THEN
          DO jj = 1, 3
              r(jj, o) = ri_new(jj)
            ENDDO
          accepted_moves_r = accepted_moves_r + 1
      ENDIF

END



SUBROUTINE total_energy(n_atoms, r, r_crit, L, E_tot)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n_atoms
    DOUBLE PRECISION, INTENT(IN) :: r(3, n_atoms), r_crit, L
    DOUBLE PRECISION, INTENT(OUT) :: E_tot
    DOUBLE PRECISION :: pot, dummy 
    INTEGER ii, jj

    E_tot = 0
    DO ii = 1, n_atoms - 1
        DO jj = ii + 1, n_atoms
            CALL lj(ii, jj, n_atoms, r, L, r_crit, pot, dummy)
            E_tot = E_tot + pot
        ENDDO
    ENDDO
    
END SUBROUTINE total_energy



SUBROUTINE energy_ii(n_atoms, ii, r, r_crit, L, E)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n_atoms, ii
      DOUBLE PRECISION, INTENT(IN) :: r(3, n_atoms), r_crit, L
      DOUBLE PRECISION, INTENT(OUT) :: E
      DOUBLE PRECISION :: pot, dummy
      INTEGER :: jj
      
      E = 0
      DO jj = 1, n_atoms
          IF (jj .NE. ii) THEN
              CALL lj(ii, jj, n_atoms, r, L, r_crit, pot, dummy)
              E = E + pot
          ENDIF
      ENDDO
      
END SUBROUTINE energy_ii



SUBROUTINE lj(is, js, n_atoms, r, L, r_crit, pot, rij_fij)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, js, n_atoms
    DOUBLE PRECISION, INTENT(IN) :: r(3, n_atoms), L, r_crit
    DOUBLE PRECISION, INTENT(OUT) :: pot, rij_fij
    DOUBLE PRECISION :: rij(3), rijl, rr2, rr
    DOUBLE PRECISION :: ynvrr2, ynvrr6, ynvrr12
    INTEGER :: kk

    rr2 = 0.d0
    pot = 0.d0
    rij_fij = 0.d0
    DO kk = 1, 3
        rijl = r(kk, js) - r(kk, is)
        rij(kk) = rijl - L * DNINT(rijl / L)
        rr2 = rr2 + rij(kk) * rij(kk)
    ENDDO

    rr = dsqrt(rr2)

    IF (rr .LT. r_crit) THEN
        ynvrr2 = 1.d0/rr2
        ynvrr6 = ynvrr2*ynvrr2*ynvrr2
        ynvrr12 = ynvrr6*ynvrr6
        rij_fij = rr2*24.d0*(2.d0*ynvrr12-ynvrr6)*ynvrr2
        pot = 4.d0*(ynvrr12-ynvrr6) 
    ELSE
        pot = 0
        rij_fij = 0
    ENDIF
    
END SUBROUTINE lj



SUBROUTINE sample(n_atoms, r, L, r_crit, nid, denom_nid, &
                  r_crit_gdr, nbins_gdr, u, p_pot, gdr)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n_atoms, nbins_gdr
    DOUBLE PRECISION, INTENT(IN) :: r(3, n_atoms), L, r_crit, nid, &
                                    denom_nid, r_crit_gdr
    DOUBLE PRECISION, INTENT(OUT) :: u, p_pot
    DOUBLE PRECISION, INTENT(INOUT) :: gdr(nbins_gdr)
    INTEGER :: ii, jj
    DOUBLE PRECISION :: pot, rij_fij

    u = 0
    p_pot = 0
    DO ii = 1, n_atoms - 1
        DO jj = ii + 1, n_atoms
            CALL lj(ii, jj, n_atoms, r, L, r_crit, pot, rij_fij)
            u = u + pot
            p_pot = p_pot + rij_fij

            CALL function_gdr(ii, jj, n_atoms, r, L, r_crit_gdr, nid, &
                              denom_nid, nbins_gdr, gdr)
        ENDDO
    ENDDO
    p_pot = p_pot / (3 * L**3)

END SUBROUTINE sample



SUBROUTINE function_gdr(is, js, n_atoms, r, L, r_crit_gdr, nid, denom_nid, &
                        nbins_gdr, gdr)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: is, js, n_atoms, nbins_gdr
    DOUBLE PRECISION, INTENT(IN) :: r(3, n_atoms), L, r_crit_gdr, nid, denom_nid
    DOUBLE PRECISION, INTENT(INOUT) :: gdr(nbins_gdr)
    INTEGER :: kk, ni
    DOUBLE PRECISION :: pot, rr2, rijl, rij(3), rr, delg

    pot = 0.d0
    rr2 = 0.d0
    DO kk = 1, 3
        rijl = r(kk, js) - r(kk, is)
        rij(kk) = rijl - L * DNINT(rijl / L)
        rr2 = rr2 + rij(kk) * rij(kk)
    ENDDO

    rr = DSQRT(rr2)
  
    delg = r_crit_gdr / nbins_gdr

    IF (rr .LT. r_crit_gdr) THEN
        ni = NINT(rr / delg)
        gdr(ni) = gdr(ni) + 2.d0 / (denom_nid * nid)
        !write(*,*) ni, gdr(ni)
    END IF
    
END



SUBROUTINE normalise_and_write_gdr(nbins_gdr, gdr, dr_gdr)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbins_gdr
    DOUBLE PRECISION, INTENT(IN) :: gdr(nbins_gdr), dr_gdr
    CHARACTER(LEN = 256) :: filename
    INTEGER :: id, ii
    DOUBLE PRECISION :: vb, r_gdr
    
    OPEN(1, FILE="tmp/current_id.tmp", STATUS="OLD")
    READ(1, *) id
    CLOSE(1)

    WRITE(filename, "('results/gdr_',i0,'.dat')") id
    OPEN(15, FILE = TRIM(ADJUSTL(filename)))
    DO ii = 1, nbins_gdr
        vb = ((ii + 1)**3 - ii**3) * dr_gdr**3
        r_gdr = dr_gdr * (ii + 1.0/2.0)
        WRITE(15, *) r_gdr, gdr(ii)/vb
    ENDDO
    CLOSE(15)

END SUBROUTINE normalise_and_write_gdr



SUBROUTINE write_last_positions(n_atoms, r, sigma)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n_atoms
    DOUBLE PRECISION, INTENT(IN) :: r(3, n_atoms), sigma
    INTEGER :: ii, jj

    OPEN(16, FILE='results/initial_pos.dat', STATUS='unknown')
    ! Write in Amstrongs!!! Un-reduced positions!!
    DO ii = 1, n_atoms
       WRITE(16, *) (r(jj,ii) * sigma, jj=1,3)
    ENDDO         
    CLOSE(16)

END SUBROUTINE write_last_positions
