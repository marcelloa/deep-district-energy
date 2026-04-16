! DEEP : District Energy & Environmental Planning Framework
!
! Copyright 2026 Marcello Aprile
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! 
program DEEP
use version_module
use district_module
use building_complex_module
use thermal_grid_module
use thermal_plants_module
use thermal_stations_module
implicit none

character(500)::input_file_name_str, &
				data_file_path_str, &
                timestamp, hours_str
type(district),target::current_district
integer::n_arguments,i,j,k,m,h,job,t,n_hours,stat,y
double precision::time,min_timestep,max_timestep,old_time ! integration time within the main TIMESTEP, [0,TIMESTEP]

! command line arguments
n_arguments=command_argument_count()
write(*,*) 'DEEP ',app_version
if (n_arguments<3)  then
	write(*,*) 'Usage: DEEP <input_file_name> <data_file_path> <timestamp> [<hours>]'
	stop ' Error: command line arguments'
end if
call get_command_argument(1,input_file_name_str)
call get_command_argument(2,data_file_path_str)
call get_command_argument(3,timestamp)
if (n_arguments>4-1) then 
    call get_command_argument(4,hours_str)
else
    hours_str='8760'
end if    
n_hours=text_to_integer(trim(hours_str))

! files names
fn%input_file_path_str=getPath(input_file_name_str)
fn%data_path_name=trim(data_file_path_str)
fn%input_file_name_ext=trim(input_file_name_str)//'.json'
fn%output_file_name_ext=trim(input_file_name_str)//'_output.json'
fn%output_file_name_bin=trim(input_file_name_str)//'_output.bin'
fn%scheduler_file_name_ext=trim(input_file_name_str)//'_scheduler.json'

simlog%logfile=trim(input_file_name_str)//'.log.'//trim(timestamp)
simstat%logfile=trim(input_file_name_str)//'.stat.'//trim(timestamp)
current_district%timestamp=trim(timestamp)

! read inputs
!call simlog%output('Reading project files '//fn%input_file_name_ext//', '//fn%output_file_name_ext,'N')
call simlog%output('Reading project file '//fn%input_file_name_ext,'N')
call current_district%read_json(fn%input_file_name_ext,fn%output_file_name_ext,fn%data_path_name,fn%input_file_path_str)
call simstat%output('10.000','N')

if (current_district%multi_year) then
    call simlog%output('Multi-year simulation.')
    simulation_years=current_district%years
    ! make a backup copy of meteo file and schedules file
    call copy_file(fn%data_path_name//'/'//trim(current_district%meteo_file_name),fn%data_path_name//'/'//trim(current_district%meteo_file_name)//'.bkp',stat)
    if (stat/=0) then
        call simlog%output('Error: unable to make a backup copy of the meteo file.')
        stop        
    end if    
    call copy_file(fn%input_file_path_str//'/'//trim(current_district%schedules_file_name),fn%input_file_path_str//'/'//trim(current_district%schedules_file_name)//'.bkp',stat)
    if (stat/=0) then
        call simlog%output('Error: unable to make a backup copy of the schedules file.')
        stop        
    end if    
else
    simulation_years=[0]
end if    
do y=1,size(simulation_years)
    current_year=y
    ! prepare input files for multi-year simulation
    if (current_district%multi_year) then
        call simlog%output('Processing year '//integer_to_text(simulation_years(y))//'...')
        call copy_file(fn%data_path_name//'/'//trim(current_district%meteo_file_name)//'.'//integer_to_text(simulation_years(y)),fn%data_path_name//'/'//trim(current_district%meteo_file_name),stat)
        if (stat/=0) then
            call simlog%output('Error: unable to copy meteo file for multi-year simulation (year '//integer_to_text(simulation_years(y))//').')
            stop        
        end if    
        call copy_file(fn%input_file_path_str//'/'//trim(current_district%schedules_file_name)//'.'//integer_to_text(simulation_years(y)),fn%input_file_path_str//'/'//trim(current_district%schedules_file_name),stat)
        if (stat/=0) then
            call simlog%output('Error: unable to copy schedules file for multi-year simulation (year '//integer_to_text(simulation_years(y))//').')
            stop        
        end if    
    end if
    ! meteo data and schedules initialization
    call current_district%building_complex%meteodata%load_file(fn%data_path_name//'/'//trim(current_district%meteo_file_name))
    call current_district%read_schedules_json(current_district%schedules_file_name,fn%input_file_path_str)
    T_ground=current_district%building_complex%meteodata%get_T_ground() ! set T_ground to the undisturbed ground temperature value.
    ! execute jobs
    do job=1,size(current_district%scheduler)
        if (current_district%scheduler(job)%execute) then    
            select case(current_district%scheduler(job)%id)
                case(1)
                if (y==1) then
                    call simlog%output('Calculate geometry and thermal zone of building units')
                    do i=1,size(current_district%building_complex%buildings)
                        call current_district%building_complex%buildings(i)%calculate_geometry()
                        call current_district%building_complex%buildings(i)%create_thermal_zones()
                    end do
                    call simlog%output('Calculate shading on external surfaces')
                    do i=1,size(current_district%building_complex%neighbors)
                        call current_district%building_complex%neighbors(i)%calculate_geometry()
                    end do
                    do m=1,12
                        do h=1,24
                            call current_district%building_complex%calculate_fsol(m,h)
                        end do
                    end do
                end if   
                case(2)
                    call simlog%output('Calculate district energy needs')
                    call current_district%building_complex%calculate_energy_needs()
                case(3)
                    call simlog%output('Calculate hvac thermal loads')
                    do i=1,size(current_district%building_complex%buildings)
                        do j=1,size(current_district%building_complex%buildings(i)%units)
                            call current_district%building_complex%buildings(i)%units(j)%calculate_hvac()
                        end do    
                    end do
                case(4)
                    if (y==1) then   
                        call simlog%output('Create thermo-hydraulic network')
                        call current_district%thermal_grids(1)%create_networks()
                    end if    
                    call simlog%output('Calculate thermo-hydraulic network energy demand')
                    call current_district%thermal_grids(1)%calculate_energy_demand(n_hours)
            end select    
        end if    
    end do    
    ! write output
    call simlog%output('Print output file')
    call current_district%output_json(fn%scheduler_file_name_ext,fn%output_file_name_ext)
    ! save output files for multi-year simulation
    if (current_district%multi_year) then
        call copy_file(fn%output_file_name_ext,fn%output_file_name_ext//'.'//integer_to_text(simulation_years(y)),stat)
        if (stat/=0) then
            call simlog%output('Error: unable to copy output file for multi-year simulation (year '//integer_to_text(simulation_years(y))//').')
            stop        
        end if    
        call copy_file(fn%output_file_name_bin,fn%output_file_name_bin//'.'//integer_to_text(simulation_years(y)),stat)
        if (stat/=0) then
            call simlog%output('Error: unable to copy output binary file for multi-year simulation (year '//integer_to_text(simulation_years(y))//').')
            stop        
        end if    
    end if    
end do
if (current_district%multi_year) then
    ! revert original meteo file and schedules file
    call copy_file(fn%data_path_name//'/'//trim(current_district%meteo_file_name)//'.bkp',&
                   fn%data_path_name//'/'//trim(current_district%meteo_file_name),stat)
    if (stat/=0) then
        call simlog%output('Error: unable to revert original meteo file.')
        stop        
    end if    
    call copy_file(fn%input_file_path_str//'/'//trim(current_district%schedules_file_name)//'.bkp',&
                   fn%input_file_path_str//'/'//trim(current_district%schedules_file_name),stat)
    if (stat/=0) then
        call simlog%output('Error: unable to revert original schedules file.')
        stop        
    end if 
end if
call simstat%output('100.000','N')
call simlog%output('All done')

end program

