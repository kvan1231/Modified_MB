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
 
		module run_star_extras 

		use star_lib
		use star_def
		use const_def
		use const_def
		use chem_def
		use binary_def
		
		implicit none
		
		integer :: time0, time1, clock_rate
		real(dp), parameter :: expected_runtime = 1 ! minutes

		integer, parameter :: restart_info_alloc = 1
		integer, parameter :: restart_info_get = 2
		integer, parameter :: restart_info_put = 3
		
		
		contains
		
		subroutine extras_controls(id, ierr)
			integer, intent(in) :: id
			integer, intent(out) :: ierr
			type (star_info), pointer :: s

			ierr = 0
			call star_ptr(id, s, ierr)
			if (ierr /= 0) return

			s% how_many_extra_history_columns => how_many_extra_history_columns
			s% data_for_extra_history_columns => data_for_extra_history_columns

			s% how_many_extra_profile_columns => how_many_extra_profile_columns
			s% data_for_extra_profile_columns => data_for_extra_profile_columns
		end subroutine extras_controls
		
		
		integer function extras_startup(id, restart, ierr)
			integer, intent(in) :: id
			logical, intent(in) :: restart
			integer, intent(out) :: ierr
			type (star_info), pointer :: s
			integer :: restart_time, prev_time_used
			ierr = 0
			call star_ptr(id, s, ierr)
			if (ierr /= 0) return
			if (.not. restart) then
				call system_clock(time0,clock_rate)
				call alloc_restart_info(s)
			else
				call unpack_restart_info(s)
				call system_clock(restart_time,clock_rate)
				prev_time_used = time1 - time0
				time1 = restart_time
				time0 = time1 - prev_time_used
			end if
			extras_startup = keep_going
		end function extras_startup
		

		integer function extras_check_model(id, id_extra)
			integer, intent(in) :: id, id_extra
			extras_check_model = keep_going
		end function extras_check_model


		integer function how_many_extra_history_columns(id, id_extra)
			integer, intent(in) :: id, id_extra
			how_many_extra_history_columns = 4
		end function how_many_extra_history_columns
		
		
		subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
			integer, intent(in) :: id, id_extra, n
			character (len=maxlen_history_column_name) :: names(n)
			real(dp) :: vals(n)
			type (star_info), pointer :: s
			integer, intent(out) :: ierr
            integer :: k, nz
            real(dp) :: turnover_time, tt_temp, tt_temp_scaled, tt_old, tt_diff
            real(dp) :: vel, vel_ratio, vel_diff, upper_lim, lower_lim, scaled_vel
            real(dp) :: dr, conv_env_r, conv_env_m, sonic_cross_time, tau_lim
            real(dp) :: eps_nuc_lim, eps_nuc, delta_mag_chk
            logical :: conv_env_found
            common/ old_var/ tt_old
			ierr = 0
			call star_ptr(id, s, ierr)
			if (ierr /= 0) return
			
			! if (n /= 1) then
			!    stop 'bad n for data_for_extra_history_columns'
			! end if

			nz = s% nz
			vel_ratio = s% x_ctrl(1)            ! originally x_ctrl(3)
			tau_lim = s% x_ctrl(2)              ! originally x_ctrl(4)

			conv_env_found = .false.

			turnover_time = 0.0
			tt_temp = 0.0
			tt_temp_scaled = 0.0

			eps_nuc_lim = 1.0d-2
			vel_diff = 0.0
			scaled_vel = 0.0

            conv_env_r = 0.0
            conv_env_m = 0.0
            sonic_cross_time = 0.0

		  ! INITIAL TURNOVER TIME CALCULATION

				do k = nz, 1, -1 ! beginning of do loop to calculate convective turnover time

					 eps_nuc = s% eps_nuc(k)
					 ! check if the cell we are looping through satisfies our convection criteria
					 if ((s% gradr(k) .gt. s% grada(k)) .and. (eps_nuc .lt. eps_nuc_lim)) then
						  ! toggle the boolean to begin integration
						  conv_env_found = .true.
					 end if

					 ! only enter this portion if the convective boolean is true
					 ! this loop will go from the innermost cell that is convective to 
					 ! the surface. This is to try and smooth through any numeric issues
					 ! with convective zones appearing and disappearing in MESA.
					 if (conv_env_found) then

						  ! loop to calculate the size of the cell, the innermost cell
						  ! needs special consideration as it is above the core
						  if (k .lt. s% nz) then
								dr = (s% r(k) - s% r(k + 1))
						  else
								dr = (s% r(k) - s% R_center)
						  end if
						  
						  ! determine the convective velocity inside each given cell
						  ! cells that 
						  if (s% gradr(k) .gt. s% grada(k)) then

								! need to ensure that the convective velocity is within
								! our defined limits, if they are outside of these limits
								! set them to be the max/min value allowed.
								vel = s% conv_vel(k)
								lower_lim = vel_ratio * s% csound(k)
								upper_lim = 1.0 * s% csound(k)

								if (vel .lt. lower_lim) then
									 vel = lower_lim
								else if (vel .gt. upper_lim) then
									 vel = upper_lim
								end if
						  
						  ! if the cell isnt defined by MESA to be convective take the
						  ! convective velocity to be equal to sound speed
						  else
								vel = s% csound(k)
						  end if

						  ! Final check involving the opacity of the given cell. If the 
						  ! cell isn't near the surface (low tau) then include it in our integration
						  if (s% tau(k) .gt. tau_lim) then
								tt_temp = tt_temp + (dr / vel)
								sonic_cross_time = sonic_cross_time + (dr / s% csound(k))
								conv_env_r = conv_env_r + dr
								conv_env_m = conv_env_m + s% dm(k)
						  end if
					 end if

				end do ! end of do loop to calculate convective turnover time

				! reset the boolean just in case
				conv_env_found = .false.


		  ! TURNOVER TIME CHECK, THIS IS TO TRY AND AVOID LARGE CHANGES

				! simply set the turnover time to the internal variable calculated above
				turnover_time = tt_temp

				if (s% model_number .gt. 1) then
					 ! calculate the variables used to check if our system is rapidly evolving
					 tt_diff = abs(tt_old - tt_temp) / tt_old
					 delta_mag_chk = s% dt / tt_old

					 ! write (*,*) "tt_diff = ", tt_diff
					 ! write (*,*) "delta_mag = ", delta_mag_chk
					 ! write (*,*) "turnover_time = ", turnover_time
					 ! write (*,*) "tt_old = ", tt_old

					 ! check if timesteps are very small or if the relative change is very large
					 if (tt_diff .gt. delta_mag_chk) then 
						  ! write (*,*) "large change, adjusting accordingly"
						  turnover_time = tt_old + (tt_temp - tt_old) * min((s% dt / tt_old), 0.5)
 
					 end if ! end of timestep/relative change check
				
				end if

				! remember the current values to be used as comparison in the next step
				tt_old = turnover_time

				! write (*,*) "outputting values"

				names(1) = "turnover_time"
				vals(1) = turnover_time

				names(2) = "conv_env_r"
				vals(2) = conv_env_r

				names(3) = "conv_env_m"
				vals(3) = conv_env_m

				names(4) = "sonic_cross_time"
				vals(4) = sonic_cross_time

		end subroutine data_for_extra_history_columns

		
		integer function how_many_extra_profile_columns(id, id_extra)
			integer, intent(in) :: id, id_extra
			integer :: ierr
			type (star_info), pointer :: s
			ierr = 0
			call star_ptr(id, s, ierr)
			if (ierr /= 0) return
			how_many_extra_profile_columns = 2
		end function how_many_extra_profile_columns
		
		
		subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
			integer, intent(in) :: id, id_extra, n, nz
			character (len=maxlen_profile_column_name) :: names(n)
			real(dp) :: vals(nz,n)
			integer, intent(out) :: ierr
			type (star_info), pointer :: s
			integer :: k
			ierr = 0
			call star_ptr(id, s, ierr)
			if (ierr /= 0) return

			names(1) = "dr"
			names(2) = "turnover_time"

			do k = 1, s% nz  
				if (k < s% nz) then
					vals(k, 1) = (s% r(k) - s% r(k + 1))
					vals(k, 2) = (s% r(k) - s% r(k + 1)) / s% conv_vel(k)
				else
					vals(k, 1) = (s% r(k) - s% R_center)
					vals(k, 2) = (s% r(k) - s% R_center) / s% conv_vel(k)
				end if
			end do

		end subroutine data_for_extra_profile_columns
		

		integer function extras_finish_step(id, id_extra)
			integer, intent(in) :: id, id_extra
			integer :: ierr
			type (star_info), pointer :: s
			ierr = 0
			call star_ptr(id, s, ierr)
			if (ierr /= 0) return
			extras_finish_step = keep_going
			call system_clock(time1,clock_rate)
			call store_restart_info(s)
		end function extras_finish_step
		
		
		subroutine extras_after_evolve(id, id_extra, ierr)
			integer, intent(in) :: id, id_extra
			integer, intent(out) :: ierr
			type (star_info), pointer :: s
			real(dp) :: dt
			ierr = 0
			call star_ptr(id, s, ierr)
			if (ierr /= 0) return
			dt = dble(time1 - time0) / clock_rate / 60
		end subroutine extras_after_evolve
		
		
		! routines for saving and restoring data so can do restarts

		
		subroutine alloc_restart_info(s)
			type (star_info), pointer :: s
			call move_restart_info(s,restart_info_alloc)
		end subroutine alloc_restart_info
		
		
		subroutine unpack_restart_info(s)
			type (star_info), pointer :: s
			call move_restart_info(s,restart_info_get)
		end subroutine unpack_restart_info
		
		
		subroutine store_restart_info(s)
			type (star_info), pointer :: s
			call move_restart_info(s,restart_info_put)
		end subroutine store_restart_info
		
		
		subroutine move_restart_info(s,op)
			type (star_info), pointer :: s
			integer, intent(in) :: op
			
			integer :: i, j, num_ints, num_dbls, ierr
			
			i = 0
			! call move_int or move_flg 
			call move_int(time0)
			call move_int(time1)
			
			num_ints = i
			
			i = 0
			! call move_dbl 
			
			num_dbls = i
			
			if (op /= restart_info_alloc) return
			if (num_ints == 0 .and. num_dbls == 0) return
			
			ierr = 0
			call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
			if (ierr /= 0) then
				write(*,*) 'failed in star_alloc_extras'
				write(*,*) 'alloc_extras num_ints', num_ints
				write(*,*) 'alloc_extras num_dbls', num_dbls
				stop 1
			end if
			
			contains
			
			subroutine move_dbl(dbl)
				real(dp) :: dbl
				i = i+1
				select case (op)
				case (restart_info_get)
					dbl = s% extra_work(i)
				case (restart_info_put)
					s% extra_work(i) = dbl
				end select
			end subroutine move_dbl
			
			subroutine move_int(int)
				integer :: int
				include 'formats'
				i = i+1
				select case (op)
				case (restart_info_get)
					!write(*,3) 'restore int', i, s% extra_iwork(i)
					int = s% extra_iwork(i)
				case (restart_info_put)
					!write(*,3) 'save int', i, int
					s% extra_iwork(i) = int
				end select
			end subroutine move_int
			
			subroutine move_flg(flg)
				logical :: flg
				i = i+1
				select case (op)
				case (restart_info_get)
					flg = (s% extra_iwork(i) /= 0)
				case (restart_info_put)
					if (flg) then
						s% extra_iwork(i) = 1
					else
						s% extra_iwork(i) = 0
					end if
				end select
			end subroutine move_flg

		end subroutine move_restart_info
		
		


		end module run_star_extras
		
