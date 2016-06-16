! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
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
         b% warn_binary_extra =.false.
      end subroutine extras_binary_controls

      subroutine jdot_mb_routine(binary_id, ierr)
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         integer :: k, nz
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         real(dp) :: turnover_time, dr, eta
         real(dp) :: rsun4, two_pi_div_p3
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         ! write(*,*) 'doing jdot'
         b% jdot_mb = 0
         s => b% s_donor
         nz = s% nz
         eta = s% x_ctrl(1)

         rsun4 = rsun*rsun*rsun*rsun
         two_pi_div_p3 = (2.0*pi/b% period)*(2.0*pi/b% period)*(2.0*pi/b% period)

         turnover_time = 0.0
         do k = 1, nz
            if (s% mixing_type(k) == convective_mixing) then
               if (k == nz) then
                  dr = s% r(k) - s% R_center
               else
                  dr = s% r(k) - s% r(k+1)
               end if
               ! s% conv_vel in cm/s
               ! dr in cm
               turnover_time = turnover_time + (dr/ s% conv_vel(k))
            end if
         end do

         ! use the formula from rappaport, verbunt, and joss.  apj, 275, 713-731. 1983.
         if (b% have_radiative_core(b% d_i) .or. b% keep_mb_on) &
            b% jdot_mb = -3.8d-30*b% m(b% d_i)*rsun4* &         
                           pow_cr(min(b% r(b% d_i),b% rl(b% d_i))/rsun,b% magnetic_braking_gamma)* &
                           two_pi_div_p3

            write (*,*), "current donor wind = ", b% mdot_system_wind(b% d_i)
            write (*,*), "MB wind boost by ", &
            b% mdot_system_wind(b% d_i) / (-1.6d12) ! Grams per second
            write (*,*), "turnover time = ", turnover_time
            write (*,*), "MB tt boost by ", &
            (turnover_time / 2.8d6 ) ** eta
            ! 2.8d6 is turnover time in seconds for a MESA model using initial mass of 
            ! 1.0 solar masses, solar metalicity, at age 4.6 Gyr

            b% jdot_mb = (b% mdot_system_wind(b% d_i) / (-1.6d12)) * ((turnover_time / 2.8d6) ** eta) * b% jdot_mb
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
         how_many_extra_binary_history_columns = 0
      end function how_many_extra_binary_history_columns
      
      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         use const_def, only: dp
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         real(dp) :: beta
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         
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
          extras_binary_startup = keep_going
      end function  extras_binary_startup
      
      !Return either rety,backup,keep_going or terminate
      integer function extras_binary_check_model(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
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
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  
         extras_binary_finish_step = keep_going
         
      end function extras_binary_finish_step
      
      subroutine extras_binary_after_evolve(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if      
         
 
      end subroutine extras_binary_after_evolve     
      
      end module run_binary_extras
