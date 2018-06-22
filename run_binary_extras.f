! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation;either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! *********************************************************************** 
      module run_binary_extras 

      use star_lib
      use star_def
      use const_def
      use const_def
      use chem_def
      use num_lib
      use binary_def
      use crlibm_lib
      
      implicit none
      
      contains
      
      subroutine extras_binary_controls(binary_id, ierr)
         integer :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr .ne. 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         ! write(*,*) 'hello from extra_binary_controls'
         b% other_mdot_edd => mdot_edd_routine
         b% other_jdot_mb => jdot_mb_routine

         b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
         b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

         b% extras_binary_startup => extras_binary_startup
         b% extras_binary_check_model => extras_binary_check_model
         b% extras_binary_finish_step => extras_binary_finish_step
         b% extras_binary_after_evolve => extras_binary_after_evolve
         b% warn_binary_extra =.false.
      end subroutine extras_binary_controls

      subroutine mdot_edd_routine(binary_id, mdot_edd, ierr)
         use const_def, only: dp
         integer, intent(in) :: binary_id
         real(dp), intent(out) :: mdot_edd
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr .ne. 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         ! changing mdot_edd from default MESA from BH to NS

         ! mdot_edd = 4*pi*clight*b% s1% x_ctrl(5)/(0.2*(1+b% s1% surface_h1))
         ! hard coding in the radius of 11.5km into the equation results in the next line.
         mdot_edd = 2.1666d18 / ((1.d0 + b% s1% surface_h1)) 
         ! write (*,*) "Modified Mdot_edd = ", mdot_edd
         ! write (*,*) "mdot_system_transfer = ", b% mdot_system_transfer(b% a_i)
         ! write (*,*) " "

      end subroutine mdot_edd_routine

      subroutine check_radiative_core(b)
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         
         real(dp) :: sum_conv, q_loc, sum_div_qloc 
         integer :: i, k, id

         include 'formats.inc'

         do i=1,2
            if (i == 1) then
               s => b% s_donor
               id = b% d_i
            else if (b% point_mass_i == 0 .and. b% include_accretor_mb) then
               s => b% s_accretor
               id = b% a_i
            else
               exit
            end if

            ! calculate how much of inner region is convective
            sum_conv = 0; q_loc = 0
            do k = s% nz, 1, -1
               q_loc = s% q(k)
               if (q_loc > 0.5d0) exit 
               if (s% mixing_type(k) == convective_mixing) &
                  sum_conv = sum_conv + s% dq(k)
            end do
            
            sum_div_qloc = (b% sum_div_qloc(id) + sum_conv/q_loc)/2
            b% sum_div_qloc(id) = sum_div_qloc
            
            if (b% have_radiative_core(id)) then ! check if still have rad core
               if (sum_div_qloc > 0.75d0) then
                  b% have_radiative_core(id) = .false.
                  write(*,*)
                  write(*,*) 'turn off magnetic braking because radiative core has gone away'
                  write(*,*)
                  ! required mdot for the implicit scheme may drop drastically,
                  ! so its neccesary to increase change factor to avoid implicit 
                  ! scheme from getting stuck
                  b% change_factor = b% max_change_factor
               end if
            else if (sum_div_qloc < 0.25d0) then ! check if now have rad core
               if (.not. b% have_radiative_core(id)) then
                  write(*,*)
                  write(*,*) 'turn on magnetic braking'
                  write(*,*)
               end if
               b% have_radiative_core(id) = .true.
            end if
         end do
            
      end subroutine check_radiative_core

      subroutine jdot_mb_routine(binary_id, ierr)
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         integer :: k, nz
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         real(dp) :: turnover_time, tt_temp, tt_old, tt_diff, tt_mesa
         real(dp) :: dr, tot_r, mb, jdot_mb
         real(dp) :: wind_factor, tt_factor, rot_factor, saturation_factor
         real(dp) :: wind_boost, tt_boost, rotation_scaling
         real(dp) :: vel, vel_ratio, upper_lim, lower_lim, tau_lim
         real(dp) :: vel_diff
         real(dp) :: rsun4, two_pi_div_p3, rad4
         real(dp) :: eps_nuc_lim, eps_nuc
         real(dp) :: mag_field, mag_temp, mag_old, mag_diff, delta_mag_chk
         common/ old_var/ mag_old, tt_old, tt_mesa
         logical :: conv_env_found
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr .ne. 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         ! write (*,*) " "
         ! write (*,*) "====================================================="
         ! write (*,*) 'doing jdot'
         ! write (*,*) "====================================================="
         ! write (*,*) " "

         s => b% s_donor
         nz = s% nz
         vel_ratio = s% x_ctrl(1) ! originally x_ctrl(3)
         tau_lim = s% x_ctrl(2) ! originally x_ctrl(4)
         wind_factor = s% x_ctrl(3)
         tt_factor = s% x_ctrl(4)
         rot_factor = s% x_ctrl(5)
         saturation_factor = s% x_ctrl(6)

         conv_env_found = .false.

         tot_r = 0.0
         turnover_time = 0.0
         tt_temp = 0.0
         eps_nuc_lim = 1.0d-2
         vel_diff = 0.0

         mag_temp = 0.0
         ! mag_old = 0.0

         do k = nz, 1, -1

            eps_nuc = s% eps_nuc(k)
            ! write(*,*) "conv_vel_old = ", s% conv_vel_old(k)
            if ((s% gradr(k) .gt. s% grada(k)) .and. (eps_nuc .lt. eps_nuc_lim)) then
               conv_env_found = .true.
            end if
            ! write(*,*) "k = ", k

            ! if (conv_env_found) then
            !    write(*,*) "conv ", s% eps_nuc(k)
            ! else
            !    write(*,*) "not_conv ", s% eps_nuc(k)
            ! end if

            if (conv_env_found) then

               if (k .lt. s% nz) then
                  dr = (s% r(k) - s% r(k + 1))
               else
                  dr = (s% r(k) - s% R_center)
               end if
               
               if (s% mixing_type(k) == convective_mixing) then
                  vel = s% conv_vel(k)
                  lower_lim = vel_ratio * s% csound(k)
                  upper_lim = 1.0 * s% csound(k)
                  if (vel .lt. lower_lim) then
                     vel = lower_lim
                  else if (vel .gt. upper_lim) then
                     vel = upper_lim
                  end if
               else
                  vel = s% csound(k)
               end if

               if (s% tau(k) .gt. tau_lim) then
                  tt_temp = tt_temp + (dr / vel)
                  ! turnover_time = turnover_time + (dr / vel)
                  tot_r = tot_r + dr
               end if
            end if

         end do

         write (*,*) "model_num = ", s% model_number
         write (*,*) "turnover_time = ", tt_temp

         if (s% model_number == 1) then
            turnover_time = tt_temp
            mag_field = (turnover_time / 2.8d6) * (2073600. / b% period)
         else
            mag_temp = (tt_temp / 2.8d6) * (2073600. / b% period)
            mag_diff = abs(mag_old - mag_temp)
            ! write (*,*) "mag_old = ", mag_old
            ! write (*,*) "mag_temp = ", mag_temp
            ! write (*,*) "mag_diff = ", mag_diff
            delta_mag_chk = s% dt / tt_temp
            ! write (*,*) "delta_mag = ", delta_mag_chk
            if ((mag_diff > delta_mag_chk) .and. (s% dt < tt_old)) then
               write (*,*) "small timestep, adjusting"
               write (*,*) "tt_temp = ", tt_temp
               tt_diff = tt_temp - tt_mesa
               write (*,*) "tt_diff = ", tt_diff
               ! if (tt_temp .lt. tt_old) then
               !    turnover_time = tt_old - tt_diff * (s% dt / tt_old)
               ! else if (tt_temp .gt. tt_old) then
               !    turnover_time = tt_old + tt_diff * (s% dt / tt_old)
               ! end if
               turnover_time = tt_old + tt_diff * (s% dt / tt_old)
               write (*,*) "turnover_time = ", turnover_time
               write (*,*) "adjusting convective velocity"

               do k = nz, 1, -1
                  vel_diff = s% conv_vel_old(k) - s% conv_vel(k)
                  s% conv_vel(k) = s% conv_vel_old(k) + vel_diff * (s% dt / tt_old)
               end do

            else 
               turnover_time = tt_temp
            end if
            mag_field = mag_temp
         end if

         tt_mesa = tt_temp
         tt_old = turnover_time
         mag_old = mag_field

         b% jdot_mb = 0
         rsun4 = pow4(rsun)
         call check_radiative_core(b)
         two_pi_div_p3 = (2.0*pi/b% period)*(2.0*pi/b% period)*(2.0*pi/b% period)

         ! use the formula from rappaport, verbunt, and joss.  apj, 275, 713-731. 1983.
         if (b% have_radiative_core(b% d_i) .or. b% keep_mb_on) then

            mb = -3.8d-30*b% m(b% d_i)*rsun4* &         
               pow_cr(min(b% r(b% d_i),b% rl(b% d_i))/rsun,b% magnetic_braking_gamma)* &
               two_pi_div_p3

            if (wind_factor .ne. 0.0) then
               ! write (*,*) "current donor wind = ", b% mdot_system_wind(b% d_i)
               wind_boost = (b% mdot_system_wind(b% d_i) / (-1.6d12)) ** wind_factor ! Grams per second
               ! write (*,*) "MB wind boost by = ", wind_boost
            else
               wind_boost = 1.0
            end if

            if (tt_factor .ne. 0.0) then
               ! 2.8d6 is turnover time in seconds for a MESA model using initial mass of 
               ! 1.0 solar masses, solar metalicity, at age 4.6 Gyr. This is approximately
               ! turnover time of the Sun.
               ! write (*,*) "turnover time = ", turnover_time
               tt_boost = (turnover_time / 2.8d6 ) ** tt_factor
               ! write (*,*) "MB tt boost by = ", tt_boost
            else
               tt_boost = 1.0
            end if

            if (rot_factor .ne. 0.0) then
               ! write (*,*) "rotation rate = ", 1 / b% period
               rotation_scaling = (2073600. / b% period ) ** rot_factor
               ! rotation_scaling = ( b% period / 2073600. ) ** rot_factor
               ! write (*,*) "MB rotation scaled by = ", rotation_scaling
            else
               rotation_scaling = 1.0
            end if

            jdot_mb = (wind_boost) * (tt_boost) * (rotation_scaling) * mb

            ! if ((wind_factor .ne. 0.0) .or. (tt_factor .ne. 0.0) .or. (rot_factor .ne. 0.0)) then
            !    write (*,*) "Magnetic Braking pre boost = ", mb
            !    write (*,*) "Magnetic Braking post scaling = ", jdot_mb
            ! end if

            ! taking the period of the sun to be 24days => 10 * P < Psun
            ! 2.4 days * 24hours/day * 60 minutes/hour * 60sec/min = 207360
            if (b% period < 207360) then
               ! write (*,*) "Rotation saturated"
               ! use the formula from Ivanova & Taam 2003 for quickly rotating stars
               b% jdot_mb = jdot_mb * ( b% period / 207360. ) ** saturation_factor
               ! write (*,*) "Magnetic Braking post saturation = ", b% jdot_mb
            else
               b% jdot_mb = jdot_mb
            end if
            
         end if

         if (b% point_mass_i == 0 .and. b% include_accretor_mb .and. &
            (b% have_radiative_core(b% a_i) .or. b% keep_mb_on)) then
            b% jdot_mb = b% jdot_mb - &
                           3.8d-30*b% m(b% a_i)*rsun4* &
                           pow_cr(min(b% r(b% a_i),b% rl(b% a_i))/rsun,b% magnetic_braking_gamma)* &
                           two_pi_div_p3
         end if

      end subroutine jdot_mb_routine

      integer function how_many_extra_binary_history_columns(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_columns = 2
      end function how_many_extra_binary_history_columns
      
      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         use const_def, only: dp
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         integer :: k, nz
         real(dp) :: turnover_time, tt_temp, tt_old, tt_diff, tt_mesa
         real(dp) :: dr, tot_r, mb, jdot_mb
         real(dp) :: vel, vel_ratio, upper_lim, lower_lim, tau_lim
         real(dp) :: vel_diff
         real(dp) :: eps_nuc_lim, eps_nuc
         real(dp) :: mag_field, mag_temp, mag_old, mag_diff, delta_mag_chk
         real(dp) :: wind_factor, tt_factor, rot_factor, saturation_factor
         real(dp) :: wind_boost, tt_boost, rotation_scaling
         common/ old_var/ mag_old, tt_old, tt_mesa
         logical :: conv_env_found
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr .ne. 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         s => b% s_donor
         nz = s% nz
         vel_ratio = s% x_ctrl(1) ! originally x_ctrl(3)
         tau_lim = s% x_ctrl(2) ! originally x_ctrl(4)
         wind_factor = s% x_ctrl(3)
         tt_factor = s% x_ctrl(4)
         rot_factor = s% x_ctrl(5)
         saturation_factor = s% x_ctrl(6)

         conv_env_found = .false.

         tot_r = 0.0
         turnover_time = 0.0
         tt_temp = 0.0
         eps_nuc_lim = 1.0d-2

         mag_temp = 0.0

         do k = nz, 1, -1
            eps_nuc = s% eps_nuc(k)

            if ((s% gradr(k) .gt. s% grada(k)) .and. (eps_nuc .lt. eps_nuc_lim)) then
               conv_env_found = .true.
            end if

            if (conv_env_found) then
               if (k .lt. s% nz) then
                  dr = (s% r(k) - s% r(k + 1))
               else
                  dr = (s% r(k) - s% R_center)
               end if
               
               if (s% mixing_type(k) == convective_mixing) then
                  vel = s% conv_vel(k)
                  lower_lim = vel_ratio * s% csound(k)
                  upper_lim = 1.0 * s% csound(k)
                  if (vel .lt. lower_lim) then
                     vel = lower_lim
                  else if (vel .gt. upper_lim) then
                     vel = upper_lim
                  end if
               else
                  vel = s% csound(k)
               end if

               if (s% tau(k) .gt. tau_lim) then
                  tt_temp = tt_temp + (dr / vel)
                  tot_r = tot_r + dr
               end if
            end if

         end do
         
         if (s% model_number == 1) then
            turnover_time = tt_temp
            mag_field = (turnover_time / 2.8d6) * (2073600. / b% period)
         else
            mag_temp = (tt_temp / 2.8d6) * (2073600. / b% period)
            mag_diff = abs(mag_old - mag_temp)
            delta_mag_chk = s% dt / tt_temp
            if ((mag_diff > delta_mag_chk) .and. (s% dt < tt_old)) then
               tt_diff = tt_temp - tt_mesa 
               ! if (tt_temp .lt. tt_old) then
               !    turnover_time = tt_old - tt_diff * (s% dt / tt_old)
               ! else if (tt_temp .gt. tt_old) then
               !    turnover_time = tt_old + tt_diff * (s% dt / tt_old)
               ! end if
               turnover_time = tt_temp + tt_diff * (s% dt / tt_old)

               do k = nz, 1, -1
                  vel_diff = s% conv_vel_old(k) - s% conv_vel(k)
                  s% conv_vel(k) = s% conv_vel_old(k) + vel_diff * (s% dt / tt_old)
               end do
            else 
               turnover_time = tt_temp
            end if
            mag_field = mag_temp
         end if

         tt_mesa = tt_temp
         tt_old = turnover_time
         mag_old = mag_field

         names(1) = "turnover_time"
         vals(1) = turnover_time

         names(2) = "mag_field"
         vals(2) = mag_field

      end subroutine data_for_extra_binary_history_columns
      
      
      integer function extras_binary_startup(binary_id,restart,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         logical, intent(in) :: restart
         call binary_ptr(binary_id, b, ierr)
         if (ierr .ne. 0) then ! failure in  binary_ptr
            return
         end if
         
!          b% s1% job% warn_run_star_extras = .false.
         ! extras_binary_startup = keep_going

         write(*,*) "starting modified MB"
         if (b% angular_momentum_j <= 0) then
            write(*,*) "angular_momentum_j <= 0", b% angular_momentum_j
            b% s1% dt_next = min(b% s1% dt * 0.50, b% s1% dt_next)
            extras_binary_startup = retry
         end if
      end function  extras_binary_startup
      
      !Return either retry,backup,keep_going or terminate
      integer function extras_binary_check_model(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         real(dp) :: j_check1
         call binary_ptr(binary_id, b, ierr)
         if (ierr .ne. 0) then ! failure in  binary_ptr
            return
         end if

         ! write (*,*) " "
         ! write (*,*) "====================================================="
         ! write (*,*) "checking model"
         ! write (*,*) "====================================================="
         ! write (*,*) " "

         ! write (*,*) "angular momentum before jdot = ", b% angular_momentum_j

         if (b% angular_momentum_j <= 0) then
            write(*,*) "bad angular momentum"
            b% s1% dt_next = min(b% s1% dt * 0.50, b% s1% dt_next)
            extras_binary_check_model = retry
         end if
        
         extras_binary_check_model = keep_going
      end function extras_binary_check_model
      
      
      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_binary_finish_step(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr .ne. 0) then ! failure in binary_ptr
            return
         end if  
         extras_binary_finish_step = keep_going
         
      end function extras_binary_finish_step
      
      subroutine extras_binary_after_evolve(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr .ne. 0) then ! failure in binary_ptr
            return
         end if      
         
 
      end subroutine extras_binary_after_evolve     
      
      end module run_binary_extras
