! Schedules different runs of mc_sampling and updates parameters automatically

PROGRAM main
    IMPLICIT NONE
    ! -- MC config --
    INTEGER :: n_steps, freq_sample, n_atoms
    DOUBLE PRECISION :: mass, sigma, eps, p_ref, t_ref, l, dr=0, dv=0
    DOUBLE PRECISION :: acceptance_r, acceptance_v
    ! -- Varying pressure --
    INTEGER :: steps_pressure, indicator_acc_r, indicator_acc_v
    DOUBLE PRECISION :: pressure_vec(20)

    INTEGER :: indicator_continue, ii, nsteps_reach_equil, nsteps_sampling  
    DOUBLE PRECISION :: p0, pf, dr_0, dv_0, l_0
    

    ! CONFIG
    steps_pressure = 10
    p0 = 0.15
    pf = 15

    nsteps_reach_equil = 1000000 !10000000

    nsteps_sampling =    5000000 !5000000
    
    freq_sample = 500

    mass = 39.948
    sigma = 3.405
    eps = 119.8

    n_atoms = 512
    t_ref = 239.6

    l_0 = 61.3
    dr_0 = 28
    dv_0 = 0.15


    ! Define pressure_vec
    DO ii = 1, steps_pressure
        pressure_vec(ii) = p0 + (ii - 1) * (pf - p0) / (steps_pressure - 1)
    ENDDO


    dr = dr_0
    dv = dv_0
    l = l_0
    DO ii = 1, steps_pressure
        p_ref = pressure_vec(ii)
        WRITE(*, *)
        WRITE(*, *)
        WRITE(*, '(A, A, F8.4)') " --- MAIN ---   ", "PRESSURE: ", p_ref

        OPEN(1, file='tmp/current_id.tmp')
        WRITE(1,*) ii
        CLOSE(1)

        indicator_continue = 0
        n_steps = nsteps_reach_equil 
        DO WHILE (indicator_continue .le. 1)
            IF (indicator_continue .eq. 1) THEN
                WRITE(*, *) "--- MAIN ---   ", "    Assume we have reached &
                            &equilibrium, now sampling..."
                n_steps = nsteps_sampling
            ELSE 
                WRITE(*, *) "--- MAIN ---   ", "    Now trying to overcome &
                            &transient state and reach equilibrium..."
            ENDIF

            ! UPDATE MC PARAMETERS
            CALL write_mc_config(n_steps, freq_sample, n_atoms, mass, t_ref, &
                                 sigma, eps, dr, dv, l, p_ref)

            ! RUN MC SIMULATION
            CALL EXECUTE_COMMAND_LINE("./mc_sampling > /dev/null")

            ! CHECK ACCEPTANCE RATE  
            OPEN(2, file='tmp/acceptance_rate.tmp')
            READ(2,*) acceptance_r
            READ(2,*) acceptance_v
            CLOSE(2)

            ! ADAPT delr, delv ACCORDINGLY
            indicator_acc_r = 1
            IF (acceptance_r .ge. 0.61) THEN
                dr = dr * (1.2)
                indicator_acc_r = 0
            ELSEIF (acceptance_r .le. 0.39) THEN
                dr = dr * (0.8)
                indicator_acc_r = 0
            ENDIF

            indicator_acc_v = 1
            IF (acceptance_v .ge. 0.61) THEN
                dv = dv * (1.2)
                indicator_acc_v = 0
            ELSEIF (acceptance_r .le. 0.39) THEN
                dv = dv * (0.8)
                indicator_acc_v = 0
            ENDIF

            IF ((indicator_acc_r .eq. 1) .and. (indicator_acc_v .eq. 1)) THEN
                indicator_continue = indicator_continue + 1
            ENDIF

            WRITE(*, '(A, I2, F8.3, A, F8.3, A, F8.4, A, F8.4, A, F8.4, I0)') &
                       " --- MAIN ---           ", ii, p_ref, ", dr=", dr, &
                       ", acc_r=", acceptance_r, ", dv=", dv, ", acc_v=", &
                       acceptance_v, indicator_continue
            ENDDO
    ENDDO

END PROGRAM main



SUBROUTINE write_mc_config(n_steps, freq_sample, n_atoms, mass, t_ref, &
                           sigma, eps, dr, dv, l, p_ref)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: n_steps, freq_sample, n_atoms
    DOUBLE PRECISION, INTENT(OUT) :: mass, t_ref, sigma, eps, &
                                     dr, dv, l, p_ref

    OPEN(1, FILE='config.dat', STATUS='unknown')
    WRITE(1, '(A)') "-------- MONTE CARLO CONFIG --------"
    WRITE(1, *) "N_STEPS", n_steps
    WRITE(1, *) "FREQ_SAMPLE", freq_sample
    WRITE(1, '(A)') "---------- PHYSICAL CTS ------------"
    WRITE(1, *) "MASS(UMA)", mass
    WRITE(1, *) "SIGMA(Å)", sigma
    WRITE(1, *) "EPSILON(K)", eps  ! Character '/' divides strings
    WRITE(1, '(A)') "-------- ENSEMBLE VARIABLES --------"
    WRITE(1, *) "N_ATOMS", n_atoms
    WRITE(1, *) "PRESSURE_REF", p_ref
    WRITE(1, *) "T_REF(K)", t_ref
    WRITE(1, '(A)') "--- VARIABLES / INITIAL GUESSES ----"
    WRITE(1, *) "L(Å)", l 
    WRITE(1, *) "DELTAR", dr
    WRITE(1, *) "DELTAV", dv
    WRITE(1, *)
    WRITE(1, *)
    WRITE(1, *) "***************************************"
    WRITE(1, *) "**  ! DO NOT CHANGE FORMAT OF FILE ! **"
    WRITE(1, *) "** ONLY EDIT THE VALUES OF VARIABLES **"
    WRITE(1, *) "***************************************"
    CLOSE(1)

END SUBROUTINE write_mc_config
