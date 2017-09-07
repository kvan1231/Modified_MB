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
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         ! write(*,*) 'hello from extra_binary_controls'
         b% other_jdot_mb => jdot_mb_routine
         b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
         b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

         b% other_mdot_edd => mdot_edd_routine
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
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         ! changing mdot_edd from default MESA from BH to NS
         ! x_ctrl(5) is the radius of the NS

         ! mdot_edd = 4*pi*clight*b% s1% x_ctrl(5)/(0.2*(1+b% s1% surface_h1))
         ! hard coding in te radius of 11.5km into the equation results in the next line.
         mdot_edd = 2.1666d18 / ((1.d0 + b% s1% surface_h1)) 
         write (*,*), "Modified Mdot_edd = ", mdot_edd
         write (*,*), "mdot_system_transfer = ", b% mdot_system_transfer(b% a_i)
         write (*,*), " "

      end subroutine mdot_edd_routine

      subroutine jdot_mb_routine(binary_id, ierr)
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         integer :: k, nz
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         real(dp) :: turnover_time, envelope_edge
         real(dp) :: dr, tot_r, mb, jdot_mb
         real(dp) :: eta, wind_fac, saturate_fac
         real(dp) :: tt_boost, wind_boost
         real(dp) :: vel_ratio, tau_lim
         real(dp) :: rsun4, two_pi_div_p3, rad4
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         write (*,*), " "
         write (*,*), "====================================================="
         write (*,*), 'doing jdot'
         write (*,*), "====================================================="
         write (*,*), " "

         s => b% s_donor
         nz = s% nz
         eta = s% x_ctrl(1)
         wind_fac = s% x_ctrl(2)
         vel_ratio = s% x_ctrl(3)
         tau_lim = s% x_ctrl(4)
         saturate_fac = s% x_ctrl(5)

         tot_r = 0.0
         turnover_time = 0.0
         envelope_edge = 0.0
         envelope_edge = max(s% conv_mx1_bot_r, s% conv_mx2_bot_r)

         do k = nz, 1, -1
            if (s% mixing_type(k) == convective_mixing) then
               if ( s% r(k) .gt. envelope_edge) then
                  if (k < s% nz) then
                     dr = (s% r(k) - s% r(k + 1))
                  else
                     dr = (s% r(k) - s% R_center)
                  end if
                  if (s% conv_vel(k) .gt. vel_ratio * s% csound(k) .and. s% tau(k) .gt. tau_lim) then
                     turnover_time = turnover_time + (dr/ s% conv_vel(k))
                     tot_r = tot_r + dr
                  end if
               else
                  turnover_time = turnover_time
                  tot_r = tot_r + dr
               end if
            end if
         end do

         ! b% jdot_mb = 0
         rsun4 = rsun*rsun*rsun*rsun 
         two_pi_div_p3 = (2.0*pi/b% period)*(2.0*pi/b% period)*(2.0*pi/b% period)

         mb = -3.8d-30*b% m(b% d_i)*rsun4* &         
                        pow_cr(min(b% r(b% d_i),b% rl(b% d_i))/rsun,b% magnetic_braking_gamma)* &
                        two_pi_div_p3

         ! use the formula from rappaport, verbunt, and joss.  apj, 275, 713-731. 1983.
         if (b% have_radiative_core(b% d_i) .or. b% keep_mb_on) &

            write (*,*), "current donor wind = ", b% mdot_system_wind(b% d_i)
            wind_boost = (b% mdot_system_wind(b% d_i) / (-1.6d12)) ** wind_fac ! Grams per second
            write (*,*), "MB wind boost by = ", wind_boost

            write (*,*), "turnover time = ", turnover_time
            tt_boost = (turnover_time / 2.8d6 ) ** eta
            write (*,*), "MB tt boost by = ", tt_boost

            ! 2.8d6 is turnover time in seconds for a MESA model using initial mass of 
            ! 1.0 solar masses, solar metalicity, at age 4.6 Gyr. This is approximately
            ! turnover time of the Sun.

            jdot_mb = (wind_boost) * (tt_boost) * mb
            write (*,*), "Magnetic Braking pre boost = ", mb
            write (*,*), "Magnetic Braking post scaling = ", jdot_mb
            ! taking the period of the sun to be 24days => 10 * P < Psun
            ! 2.4 days * 24hours/day * 60 minutes/hour * 60sec/min = 207360
            if (b% period < 207360) then
               write (*,*), "Rotation saturated"
               ! use the formula from Ivanova & Taam 2003 for quickly rotating stars
               rad4 = b% r(b% d_i) * b% r(b% d_i) * b% r(b% d_i) * b% r(b% d_i)
               b% jdot_mb =  (-6.0d30 * rad4 / rsun4 * 10 ** (1.7) *&
                             (2073600 / b% period) ** saturate_fac) * tt_boost * wind_boost
            else
               b% jdot_mb = jdot_mb
            end if
            write (*,*), "Magnetic Braking post saturation = ", b% jdot_mb

         if (b% evolve_both_stars .and. b% include_accretor_mb .and. &
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
         how_many_extra_binary_history_columns = 3
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
         real(dp) :: turnover_time, envelope_edge
         real(dp) :: dr, tot_r
         real(dp) :: eta, wind_fac
         real(dp) :: vel_ratio, tau_lim

         ierr = 0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         s => b% s_donor
         nz = s% nz
         eta = s% x_ctrl(1)
         wind_fac = s% x_ctrl(2)
         vel_ratio = s% x_ctrl(3)
         tau_lim = s% x_ctrl(4)

         tot_r = 0.0
         turnover_time = 0.0
         envelope_edge = 0.0
         envelope_edge = max(s% conv_mx1_bot_r, s% conv_mx2_bot_r)

         do k = nz, 1, -1
            if (s% mixing_type(k) == convective_mixing) then
               if ( s% r(k) .gt. envelope_edge) then
                  if (k < s% nz) then
                     dr = s% r(k) - s% r(k + 1)
                  else
                     dr = s% r(k) - s% R_center
                  end if
                  if (s% conv_vel(k) .gt. vel_ratio * s% csound(k) .and. s% tau(k) .gt. tau_lim) then
                     turnover_time = turnover_time + (dr/ s% conv_vel(k))
                     tot_r = tot_r + dr
                  end if
               end if
            end if
         end do

         names(1) = "turnover_time"
         vals(1) = turnover_time

         names(2) = "total_r"
         vals(2) = tot_r

         names(3) = "mdot_edd"
         vals(3) = 4*pi*clight*b% s1% x_ctrl(5)/(0.2*(1+b% s1% surface_h1))

      end subroutine data_for_extra_binary_history_columns
      
      
      integer function extras_binary_startup(binary_id,restart,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         logical, intent(in) :: restart
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
         
!          b% s1% job% warn_run_star_extras = .false.
         ! extras_binary_startup = keep_going

         write(*,*) "starting modified MB"
         if (b% angular_momentum_j <= 0) then
            write(*,*), "angular_momentum_j <= 0", b% angular_momentum_j
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
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if

         write (*,*), " "
         write (*,*), "====================================================="
         write (*,*), "checking model"
         write (*,*), "====================================================="
         write (*,*), " "

         ! write (*,*), "angular momentum before jdot = ", b% angular_momentum_j

         if (b% angular_momentum_j <= 0) then
            write(*,*), "bad angular momentum"
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
         if (ierr /= 0) then ! failure in binary_ptr
            return
         end if  
         extras_binary_finish_step = keep_going
         
      end function extras_binary_finish_step
      
      subroutine extras_binary_after_evolve(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in binary_ptr
            return
         end if      
         
 
      end subroutine extras_binary_after_evolve     
      
      end module run_binary_extras
