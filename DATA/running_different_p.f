      PROGRAM running_different_p
        IMPLICIT NONE
        INTEGER :: ncycle, nsamp, natoms
        INTEGER ::steps_pressure, indicator_acc_r, indicator_acc_v
        INTEGER ::indicator_continue, ii       
        DOUBLE PRECISION :: mass, tref, sigma, epsil, delr, delv
        DOUBLE PRECISION :: pressure_ref, acceptance_r, acceptance_v
        DOUBLE PRECISION :: p0, pf
        DOUBLE PRECISION :: pressure_vec(20)
        DOUBLE PRECISION :: mean_en_vec(20), std_en_vec(20)
        DOUBLE PRECISION :: mean_pvirial_vec(20), std_pvirial_vec(20)
        DOUBLE PRECISION :: r_gdr_vec(40), gdr_vec(20)
      
        !config pressure range
        steps_pressure = 10
        p0 = 0.15
        pf = 15

        !define pressure_vec
        DO ii = 1, steps_pressure
          pressure_vec(ii) = p0 + (ii-1)*(pf-p0)/(steps_pressure-1)
          write(*,*) pressure_vec(ii)
        ENDDO

        !define constants
        ncycle = 10000000
        nsamp = 500
        natoms = 500
        mass = 39.948d0
        tref = 2
        sigma = 3.405d0
        epsil = 119.8d0
        
        
        delr = 28
        delv = 0.15
        
        DO ii = 1, steps_pressure
          pressure_ref = pressure_vec(ii)

          OPEN(1, file='actual_counter.dat')
          WRITE(1,*) ii
          CLOSE(1)

          indicator_continue = 0
          ncycle = 1000000
          DO WHILE (indicator_continue .le. 1)

          IF (indicator_continue .eq. 1) THEN
            ncycle = 5000000
          ENDIF

          !WRITE LEAP-LJ.DATA FILE
          OPEN(1,file='leap-lj.data')
          WRITE(1,*) ncycle, nsamp
          WRITE(1,*) natoms
          WRITE(1,*) mass, tref
          WRITE(1,*) sigma, epsil
          WRITE(1,*) delr, delv
          WRITE(1,*) pressure_ref
          WRITE(1,*) '**************************************************'
          WRITE(1,*) '1.ncycl, nsamp (number of cycles between samples)'
          WRITE(1,*) '2.natoms'
          WRITE(1,*) '3.mass, Tref (Tref already in reduced units!)'
          WRITE(1,*) '4.sigma,epsil (LJ parameters in A and K)'
          WRITE(1,*) '5.deltar, deltav (time step in ps), ~0.05,20'
          WRITE(1,*) '6.pressure (reference, should be cte!)'
          CLOSE(1)

          !RUN MC SIMULATION
          CALL EXECUTE_COMMAND_LINE ("./leapfroglj")

          !CHECK ACCEPTANCE RATE  
          OPEN(2, file='acceptance_rate.dat')
          READ(2,*) acceptance_r
          READ(2,*) acceptance_v
          CLOSE(2)
          
          !ADAPT delr, delv ACCORDINGLY
          indicator_acc_r = 1
          IF ( acceptance_r .ge. 0.61 ) THEN
            delr = delr * (1.2)
            indicator_acc_r = 0
          ELSEIF ( acceptance_r .le. 0.39 ) THEN
            delr = delr * (0.8)
            indicator_acc_r = 0
          ENDIF
          
          indicator_acc_v = 1
          IF ( acceptance_v .ge. 0.61 ) THEN
            delv = delv * (1.2)
            indicator_acc_v = 0
          ELSEIF ( acceptance_r .le. 0.39 ) THEN
            delv = delv * (0.8)
            indicator_acc_v = 0
          ENDIF
          
          IF ( (indicator_acc_r .eq. 1) .and. (indicator_acc_v .eq. 1) )
       &      THEN
            indicator_continue = indicator_continue + 1
          ENDIF
          
          write(*,*) ii, pressure_ref, "delr",delr, "acc_r=",acceptance_r,
       &          "delv", delv, "acc_v=", acceptance_v, indicator_continue
          ENDDO
        ENDDO

      END PROGRAM running_different_p
