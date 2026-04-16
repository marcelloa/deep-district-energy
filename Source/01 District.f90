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
module district_module

	use building_complex_module
    use thermal_stations_module
    use json_module
    implicit none

    type calculation_step
         integer::id
         character(:),allocatable::step
         character(:),allocatable::timestamp
         logical::execute
    end type

   	type district
		type(building_complex)::building_complex
		type(thermal_grid),dimension(:),allocatable::thermal_grids
        character(:),allocatable::timestamp        
        character(:),allocatable::schedules_file_name
        character(:),allocatable::meteo_file_name
        type(calculation_step),dimension(:),allocatable::scheduler
        logical::multi_year
        integer,dimension(:),allocatable::years
	contains
		procedure::read_json=>district_read_json
		procedure::read_materials_json=>district_read_materials_json
		procedure::read_structures_json=>district_read_structures_json
		procedure::read_schedules_json=>district_read_schedules_json
        procedure::output_json=>district_output_json
	end type	

contains


  subroutine district_read_materials_json(this,file_name,data_path_name)
        class(district),target::this
        character(*)::file_name,data_path_name
        type(json_file) :: jsonf_input
        type(json_core) :: jsonc_input
        type(json_value), pointer:: p_input, child_input
        integer::n_input_elements,n_input_layers
        logical::found
        integer::i,j,k
        character(len=80) :: str_name, str_layer_name
        character(len=20) :: str_id,layer_id   
        call jsonf_input%initialize()
	    call jsonf_input%load(trim(data_path_name)//'/'//trim(file_name))
        if (jsonf_input%failed()) then    
            call jsonf_input%print_error_message()
            call simlog%output('Error: input file for materials not found.')
            stop
        end if          
        call jsonc_input%initialize()
        ! materials
        call jsonf_input%info('', found, n_children=n_input_elements)
        if (.not.found) then   
        	call simlog%output('Error: array [materials] not found in input json file.')
            stop
        end if
        allocate(this%building_complex%materials(n_input_elements))
        do i=1,n_input_elements
            write (str_id, *) i
            str_name='('//trim(adjustl(str_id))//')' ! create name of the structure element to find, i.e. (<i>)
            call jsonf_input%get(trim(str_name), p_input, found) ! now p point to the i-th element of the structure that we want to read
            call jsonc_input%get(p_input,'name',this%building_complex%materials(i)%name)  ! get the name of the structure element
            call jsonc_input%get(p_input,'description',this%building_complex%materials(i)%description)  ! get the name of the structure element
            call jsonc_input%get(p_input,'lambda',this%building_complex%materials(i)%properties%lambda)  
            call jsonc_input%get(p_input,'rho',this%building_complex%materials(i)%properties%rho)  
            call jsonc_input%get(p_input,'cp',this%building_complex%materials(i)%properties%cp)  
            call jsonc_input%get(p_input,'vrf',this%building_complex%materials(i)%properties%vrf)  
            call jsonc_input%get(p_input,'resistance',this%building_complex%materials(i)%properties%resistance)  
            if (jsonc_input%failed()) then    
                call jsonc_input%print_error_message()
                call simlog%output('Error: materials data are corrupted.')
                stop
            end if    
        end do        
        call jsonf_input%destroy()
        if (jsonf_input%failed()) then
        	call simlog%output('Error: unable to destroy input json file')
            stop
        end if
        call jsonc_input%destroy()
        if (jsonc_input%failed()) then
        	call simlog%output('Error: unable to destroy json object')
            stop
        end if        
    end subroutine

    subroutine district_read_structures_json(this,file_name,data_path_name)
        class(district),target::this
        character(*)::file_name,data_path_name
        type(json_file) :: jsonf_input
        type(json_core) :: jsonc_input
        type(json_value), pointer:: p_input, child_input
        integer::n_input_elements,n_input_layers
        logical::found
        integer::i,j,k
        character(len=80) :: str_name, str_layer_name
        character(len=20) :: str_id,layer_id   
        call jsonf_input%initialize()
	    call jsonf_input%load(trim(data_path_name)//'/'//trim(file_name))
        if (jsonf_input%failed()) then    
            call jsonf_input%print_error_message()
            call simlog%output('Error: input file for structures not found.')
            stop
        end if          
        call jsonc_input%initialize()
        ! wall_slab_structures
        call jsonf_input%info('wall_slab_structures', found, n_children=n_input_elements)
        if (.not.found) then
        	call simlog%output('Error: structure named [wall_slab_structures] not found in input json file.')
            stop
        end if
        allocate(this%building_complex%wall_slab_structures(n_input_elements))
        do i=1,n_input_elements
            write (str_id, *) i
            str_name='wall_slab_structures('//trim(adjustl(str_id))//')' ! create name of the structure element to find, i.e. wall_slab_structures(<i>)
            call jsonf_input%get(trim(str_name), p_input, found) ! now p point to the i-th element of the structure that we want to read
            call jsonc_input%get(p_input,'name',this%building_complex%wall_slab_structures(i)%name)  ! get the name of the structure element        
            call jsonc_input%get(p_input,'ISOclass',this%building_complex%wall_slab_structures(i)%ISOclass)
            call jsonc_input%get(p_input,'internal_resistance',this%building_complex%wall_slab_structures(i)%internal_resistance)  
            call jsonc_input%get(p_input,'external_resistance',this%building_complex%wall_slab_structures(i)%external_resistance)  
            ! get number of layers
            call jsonc_input%info(p_input,'layers', found, n_children=n_input_layers)    
            allocate(this%building_complex%wall_slab_structures(i)%layers(n_input_layers))
            do j=1,n_input_layers
                write (layer_id, *) j
                str_layer_name=trim(str_name)//'.layers('//trim(adjustl(layer_id))//')' ! create name of the structure element to find, i.e. layers(<j>)                
                call jsonf_input%get(trim(str_layer_name), p_input, found) ! now p point to the i-th element of the structure that we want to read
                call jsonc_input%get(p_input,'material_name',this%building_complex%wall_slab_structures(i)%layers(j)%material_name)  ! get the name of the material element
                call jsonc_input%get(p_input,'thickness',this%building_complex%wall_slab_structures(i)%layers(j)%thickness)  
                if (jsonc_input%failed()) then    
                    call jsonc_input%print_error_message()
                    call simlog%output('Error: layers data are corrupted.')
                    stop
                end if    
            end do    
            if (jsonc_input%failed()) then    
                call jsonc_input%print_error_message()
                call simlog%output('Error: structures data are corrupted.')
                stop
            end if    
        end do        
        ! windows
        call jsonf_input%info('window_structures', found, n_children=n_input_elements)
        if (.not.found) then
        	call simlog%output('Error: structure named [window_structures] not found in input json file.')
            stop
        end if
        allocate(this%building_complex%window_structures(n_input_elements))
        do i=1,n_input_elements
            write (str_id, *) i
            str_name='window_structures('//trim(adjustl(str_id))//')' ! create name of the structure element to find, i.e. windows(<i>)
            call jsonf_input%get(trim(str_name), p_input, found) ! now p point to the i-th element of the structure that we want to read
            call jsonc_input%get(p_input,'name',this%building_complex%window_structures(i)%name)  ! get the name of the structure element
            call jsonc_input%get(p_input,'description',this%building_complex%window_structures(i)%description)  ! get the name of the structure element
            call jsonc_input%get(p_input,'Ug',this%building_complex%window_structures(i)%Ug)  
            call jsonc_input%get(p_input,'gn',this%building_complex%window_structures(i)%gn)  
            call jsonc_input%get(p_input,'Uf',this%building_complex%window_structures(i)%Uf)  
            call jsonc_input%get(p_input,'Htb',this%building_complex%window_structures(i)%Htb)  
            call jsonc_input%get(p_input,'Ff',this%building_complex%window_structures(i)%Ff)  
            call jsonc_input%get(p_input,'Fw',this%building_complex%window_structures(i)%Fw)  
            if (jsonc_input%failed()) then    
                call jsonc_input%print_error_message()
                call simlog%output('Error: windows data are corrupted.')
                stop
            end if    
        end do      
        call jsonf_input%destroy()
        if (jsonf_input%failed()) then
        	call simlog%output('Error: unable to destroy input json file')
            stop
        end if
        call jsonc_input%destroy()
        if (jsonc_input%failed()) then
        	call simlog%output('Error: unable to destroy json object')
            stop
        end if        
    end subroutine
      
    subroutine district_read_schedules_json(this,file_name,data_path_name)
        class(district)::this
        character(*)::file_name,data_path_name
        type(json_file) :: jsonf_input
        type(json_core) :: jsonc_input
        type(json_value), pointer:: p_input, child_input
        integer::n_input_elements
        logical::found,first_allocation
        integer::i,j,k
        character(len=80) :: str_name
        character(len=20) :: str_id 
        character(:),allocatable::str
        call jsonf_input%initialize()
	    call jsonf_input%load(trim(data_path_name)//'/'//trim(file_name))
        if (jsonf_input%failed()) then    
            call jsonf_input%print_error_message()
            call simlog%output('Error: input file for schedules not found.')
            stop
        end if          
        call jsonc_input%initialize()
        ! schedules
        call jsonf_input%info('schedules', found, n_children=n_input_elements)
        if (.not.found) then
        	call simlog%output('Error: structure named [schedules] not found in input json file.')
            stop
        end if
        if (.not.associated(this%building_complex%schedules)) then
           allocate(this%building_complex%schedules(n_input_elements))
           first_allocation=.true.
        else
           first_allocation=.false.            
        end if
        ! integrity check: if it is not first allocation, check the content of the schedules file is the same as the one already in memory
        if (.not.first_allocation) then
            if (n_input_elements/=size(this%building_complex%schedules)) then
            	call simlog%output('Error: the schedules file being loaded does not have the expected number of elements.')
                stop
            end if    
        end if    
        do i=1,n_input_elements
            write (str_id, *) i
            str_name='schedules('//trim(adjustl(str_id))//')' ! create name of the structure element to find, i.e. schedules(<i>)
            call jsonf_input%get(trim(str_name), p_input, found) ! now p point to the i-th element of the structure that we want to read
            call jsonc_input%get(p_input,'name',str)  ! get the name of the structure element
            if (.not.first_allocation) then
                if (str/=this%building_complex%schedules(i)%name) then
            	    call simlog%output('Error: the schedules file being loaded does not have the expected schedule name ('//this%building_complex%schedules(i)%name//') in position '//integer_to_text(i)//'.')
                    stop
                end if
            else
                this%building_complex%schedules(i)%name=str
            end if    
            call jsonc_input%get(p_input,'timestep',this%building_complex%schedules(i)%timestep)  
            call jsonc_input%get(p_input,'size',this%building_complex%schedules(i)%size)  
            call jsonc_input%get(p_input,'values',this%building_complex%schedules(i)%values)  
            if (jsonc_input%failed()) then    
                call jsonc_input%print_error_message()
            	call simlog%output('Error: schedule data are corrupted.')
                stop
            end if    
        end do        
        call jsonf_input%destroy()
        if (jsonf_input%failed()) then
        	call simlog%output('Error: unable to destroy input json file')
            stop
        end if
        call jsonc_input%destroy()
        if (jsonc_input%failed()) then
        	call simlog%output('Error: unable to destroy json object')
            stop
        end if        
    end subroutine

    subroutine district_read_json(this,file_name,file_output_name,data_path_name,input_path_name)
        class(district),target::this
        character(*)::file_name,file_output_name,data_path_name,input_path_name
        character(:),allocatable::materials_file_name,structures_file_name,tmp_str
		logical::is_binary
        type(json_file) :: jsonf_input
        type(json_core) :: jsonc_input
        type(json_value), pointer:: p_input, child_input,jmat,jrow,jnext
        integer::n_input_elements,n_input_units,n_input_units_materials,n_input_external_structures,n_input_fsol_1,n_input_fsol_2,n_input_fsol_3
        integer::n_input_vertices,n_input_grid,n_input_segments,n_input_windows,n_input_childs
        logical::found
        integer::i,j,k,l,m,nrow,ncol
        character(len=80) :: str_name,str_unit_name,str_unit_material_name,str_structure_name,str_par_name
        character(len=20) :: str_id,unit_id,material_id,structure_id,l_id,m_id,segment_id 
        character(:),allocatable::material_name,structure_name,window_name
        double precision::material_mass
        double precision,dimension(:),allocatable::vec_dbl
		type(textbuffer)::buffer
        call jsonf_input%initialize()
	    call jsonf_input%load(file_name)
        if (jsonf_input%failed()) then    
            call jsonf_input%print_error_message()
            call simlog%output('Error: input file not found or corrupted.') 
            stop
        end if          
        call jsonc_input%initialize()
        ! district main data
        call jsonf_input%get('name',this%building_complex%name,found)
        call jsonf_input%info('scheduler', found, n_children=n_input_elements)
        if (.not.found) then
        	call simlog%output('Error: structure named [scheduler] not found in input json file.')
            stop
        end if
        allocate(this%scheduler(n_input_elements))
        do i=1,n_input_elements
            write(str_id,*) i
            call jsonf_input%get('scheduler['//trim(adjustl(str_id))//'].id',this%scheduler(i)%id,found); if (.not.found) stop
            call jsonf_input%get('scheduler['//trim(adjustl(str_id))//'].step',this%scheduler(i)%step,found); if (.not.found) stop
            call jsonf_input%get('scheduler['//trim(adjustl(str_id))//'].timestamp',this%scheduler(i)%timestamp,found); if (.not.found) stop
            call jsonf_input%get('scheduler['//trim(adjustl(str_id))//'].execute',this%scheduler(i)%execute,found); if (.not.found) stop
        end do
        call jsonf_input%get('position.latitude',this%building_complex%position%latitude,found)
        call jsonf_input%get('position.longitude',this%building_complex%position%longitude,found)
        ! read multi-year simulation inputs
        call jsonf_input%get('multi_year',this%multi_year,found)
        call jsonf_input%get('years',this%years,found)        
        ! read meteo file name
        call jsonf_input%get('meteofile',this%meteo_file_name,found)
        ! read materials
        call jsonf_input%get('materials',materials_file_name,found)
        call this%read_materials_json(materials_file_name,data_path_name)
        ! read building structures
        call jsonf_input%get('structures',structures_file_name,found)
        call this%read_structures_json(structures_file_name,input_path_name)
        ! read building schedules
        call jsonf_input%get('schedules',this%schedules_file_name,found)
        ! calculate windows transmittance and resistance
        do i=1,size(this%building_complex%window_structures)
            call this%building_complex%window_structures(i)%calculate()
        end do        
        ! copy materials properties in building structures data and calculate transmittance values
        do i=1,size(this%building_complex%wall_slab_structures)
            do j=1,size(this%building_complex%wall_slab_structures(i)%layers)
                do k=1,size(this%building_complex%materials)
                    if (this%building_complex%wall_slab_structures(i)%layers(j)%material_name==this%building_complex%materials(k)%name) then
                        this%building_complex%wall_slab_structures(i)%layers(j)%properties=this%building_complex%materials(k)%properties
                        ! calculate layer resistance
                        if (this%building_complex%wall_slab_structures(i)%layers(j)%properties%lambda>0.d0 .and. this%building_complex%wall_slab_structures(i)%layers(j)%properties%resistance==0d0) then
                            if (this%building_complex%wall_slab_structures(i)%layers(j)%thickness<=0d0) then 
                                call simlog%output('Error: thickness of layer '//int2str(j)//' in wall/slab '//this%building_complex%wall_slab_structures(i)%name//' does not have a positive value','A')
                                stop
                            end if    
                            if (this%building_complex%wall_slab_structures(i)%layers(j)%properties%lambda<=0d0) then 
                                call simlog%output('Error: lambda of layer '//int2str(j)//' in wall/slab '//this%building_complex%wall_slab_structures(i)%name//' does not have a positive value','A')
                                stop
                            end if    
                            this%building_complex%wall_slab_structures(i)%layers(j)%properties%resistance=this%building_complex%wall_slab_structures(i)%layers(j)%thickness/this%building_complex%wall_slab_structures(i)%layers(j)%properties%lambda
                            else
                        end if  
                    end if    
                end do    
            end do
            ! calculate transmittance and resistance
            call this%building_complex%wall_slab_structures(i)%calculate()
        end do    
        ! buildings
        call jsonf_input%info('buildings', found, n_children=n_input_elements)
        if (.not.found) then
        	call simlog%output('Error: structure named [buildings] not found in input json file.')
            stop
        end if
        allocate(this%building_complex%buildings(n_input_elements))
        do i=1,n_input_elements
            write (str_id, *) i
            str_name='buildings('//trim(adjustl(str_id))//')' ! create name of the structure element to find, i.e. buildings(<i>)
            call jsonf_input%get(trim(str_name), p_input, found) ! now p point to the i-th element of the structure that we want to read
            this%building_complex%buildings(i)%id=i
            call jsonc_input%get(p_input,'name',this%building_complex%buildings(i)%name)  ! get the name of the structure element
            call jsonc_input%get(p_input, 'X', this%building_complex%buildings(i)%X)
            call jsonc_input%get(p_input, 'Y', this%building_complex%buildings(i)%Y)
            call jsonc_input%get(p_input,'floors_ag',this%building_complex%buildings(i)%floors_ag)  
            call jsonc_input%get(p_input,'floors_bg',this%building_complex%buildings(i)%floors_bg)  
            call jsonc_input%get(p_input,'height_ag',this%building_complex%buildings(i)%height_ag)  
            call jsonc_input%get(p_input,'height_bg',this%building_complex%buildings(i)%height_bg)  
            ! get number of units
            call jsonc_input%info(p_input,'units', found, n_children=n_input_units)    
            allocate(this%building_complex%buildings(i)%units(n_input_units))
            do j=1,n_input_units
                this%building_complex%buildings(i)%units(j)%building_id=i
                this%building_complex%buildings(i)%units(j)%p_buildings=>this%building_complex%buildings
                this%building_complex%buildings(i)%units(j)%p_building_complex=>this%building_complex
                write (unit_id, *) j
                str_unit_name=trim(str_name)//'.units('//trim(adjustl(unit_id))//')' ! create name of the structure element to find, i.e. units(<j>)                
                call jsonf_input%get(trim(str_unit_name), p_input, found) ! now p point to the i-th element of the structure that we want to read
                call jsonc_input%get(p_input,'name',this%building_complex%buildings(i)%units(j)%name)  
                call jsonc_input%get(p_input,'footprint_to_gross_area',this%building_complex%buildings(i)%units(j)%footprint_to_gross_area)
                call jsonc_input%get(p_input,'gross_to_heated_area',this%building_complex%buildings(i)%units(j)%gross_to_heated_area)
                call jsonc_input%get(p_input,'heated_to_useful_area',this%building_complex%buildings(i)%units(j)%heated_to_useful_area)
                this%building_complex%buildings(i)%units(j)%useful_floor_height=0.9 ! default value, can be improved considering interfloor slab thickness...
                call jsonc_input%get(p_input,'window_to_useful_area',this%building_complex%buildings(i)%units(j)%window_to_useful_area)
                call jsonc_input%get(p_input,'window_to_segment',this%building_complex%buildings(i)%units(j)%window_to_segment)
                call jsonc_input%get(p_input,'wall_is_adiabatic',this%building_complex%buildings(i)%units(j)%wall_is_adiabatic)
                call jsonc_input%get(p_input,'floor_start_type',tmp_str); this%building_complex%buildings(i)%units(j)%floor_start_type=decode_floor_type(tmp_str)
                call jsonc_input%get(p_input,'floor_end_type', tmp_str);  this%building_complex%buildings(i)%units(j)%floor_end_type=decode_floor_type(tmp_str)
                call jsonc_input%get(p_input,'floor_start',this%building_complex%buildings(i)%units(j)%floor_start)
                call jsonc_input%get(p_input,'floor_end',this%building_complex%buildings(i)%units(j)%floor_end)
                call jsonc_input%get(p_input,'name_slab_basement',this%building_complex%buildings(i)%units(j)%name_slab_basement)
                call jsonc_input%get(p_input,'name_slabs_interfloor',this%building_complex%buildings(i)%units(j)%name_slab_interfloor)
                call jsonc_input%get(p_input,'name_roof',this%building_complex%buildings(i)%units(j)%name_roof)
                call jsonc_input%get(p_input,'name_walls',this%building_complex%buildings(i)%units(j)%name_walls)
                call jsonc_input%get(p_input,'name_windows',this%building_complex%buildings(i)%units(j)%name_windows)
                call jsonc_input%get(p_input,'linear_thermal_bridge',this%building_complex%buildings(i)%units(j)%linear_thermal_bridge)
                call jsonc_input%get(p_input,'min_RH',this%building_complex%buildings(i)%units(j)%min_RH)
                call jsonc_input%get(p_input,'max_RH',this%building_complex%buildings(i)%units(j)%max_RH)
                call jsonc_input%get(p_input,'schedule_ventilation',this%building_complex%buildings(i)%units(j)%schedule_ventilation)
                call jsonc_input%get(p_input,'schedule_occupancy',this%building_complex%buildings(i)%units(j)%schedule_occupancy)
                call jsonc_input%get(p_input,'schedule_appliances',this%building_complex%buildings(i)%units(j)%schedule_appliances)
                call jsonc_input%get(p_input,'schedule_heating_setpoint',this%building_complex%buildings(i)%units(j)%schedule_heating_setpoint)
                call jsonc_input%get(p_input,'schedule_cooling_setpoint',this%building_complex%buildings(i)%units(j)%schedule_cooling_setpoint)
                call jsonc_input%get(p_input,'schedule_dhw',this%building_complex%buildings(i)%units(j)%schedule_dhw)
                call jsonc_input%get(p_input,'ventilation_rate',this%building_complex%buildings(i)%units(j)%ventilation_rate)
                call jsonc_input%get(p_input,'crowd_density',this%building_complex%buildings(i)%units(j)%crowd_density)
                call jsonc_input%get(p_input,'power_density',this%building_complex%buildings(i)%units(j)%power_density)
                call jsonc_input%get(p_input,'dhw_use',this%building_complex%buildings(i)%units(j)%dhw_use)
                call jsonc_input%get(p_input,'blinds_G',this%building_complex%buildings(i)%units(j)%blinds_G)
                call jsonc_input%get(p_input,'blinds_b',this%building_complex%buildings(i)%units(j)%blinds_b)
                call jsonc_input%get(p_input,'hvac_type',this%building_complex%buildings(i)%units(j)%hvac_type)
                call jsonc_input%get(p_input,'dhw_type',this%building_complex%buildings(i)%units(j)%dhw_type)
                call jsonc_input%get(p_input,'emission_efficiency',this%building_complex%buildings(i)%units(j)%emission_efficiency)
                call jsonc_input%get(p_input,'control_efficiency',this%building_complex%buildings(i)%units(j)%control_efficiency)
                call jsonc_input%get(p_input,'distribution_efficiency',this%building_complex%buildings(i)%units(j)%distribution_efficiency)
                call jsonc_input%get(p_input,'distribution_volumetric_efficiency',this%building_complex%buildings(i)%units(j)%distribution_volumetric_efficiency)
                call jsonc_input%get(p_input,'heat_recovery_effectiveness',this%building_complex%buildings(i)%units(j)%heat_recovery_effectiveness)
                call jsonc_input%get(p_input,'maximum_airvolume_changes',this%building_complex%buildings(i)%units(j)%maximum_airvolume_changes)
                call jsonc_input%get(p_input,'average_flat_area',this%building_complex%buildings(i)%units(j)%average_flat_area)
                call jsonc_input%get(p_input,'dhw_peak_flow',this%building_complex%buildings(i)%units(j)%dhw_peak_flow)
                call jsonc_input%get(p_input,'dhw_peak_duration',this%building_complex%buildings(i)%units(j)%dhw_peak_duration)
                call jsonc_input%get(p_input,'dhw_preheating',this%building_complex%buildings(i)%units(j)%dhw_preheating)
                call jsonc_input%get(p_input,'water_main_temperature',this%building_complex%buildings(i)%units(j)%water_main_temperature)
                call jsonc_input%get(p_input,'dhw_storage_temperature',this%building_complex%buildings(i)%units(j)%dhw_storage_temperature)
                ! optional fields
                call jsonc_input%info(p_input,'Q_h_specific', found)
                if (found) then 
                    call jsonc_input%get(p_input,'Q_h_specific',this%building_complex%buildings(i)%units(j)%Q_h_specific)
                else
                    this%building_complex%buildings(i)%units(j)%Q_h_specific=50d0
                end if
                call jsonc_input%info(p_input,'Q_c_specific', found)
                if (found) then 
                    call jsonc_input%get(p_input,'Q_c_specific',this%building_complex%buildings(i)%units(j)%Q_c_specific)
                else
                    this%building_complex%buildings(i)%units(j)%Q_c_specific=80d0
                end if
                call jsonc_input%info(p_input,'schedule_power_available', found)
                if (found) then 
                    call jsonc_input%get(p_input,'schedule_power_available',this%building_complex%buildings(i)%units(j)%schedule_power_available)
                else
                    this%building_complex%buildings(i)%units(j)%schedule_power_available='DEFAULT'
                end if    
                call jsonc_input%info(p_input,'schedule_internal_gain', found)
                if (found) then 
                    call jsonc_input%get(p_input,'schedule_internal_gain',this%building_complex%buildings(i)%units(j)%schedule_internal_gain)
                else
                    this%building_complex%buildings(i)%units(j)%schedule_internal_gain='DEFAULT'
                end if    
                if (jsonc_input%failed()) then    
                    call jsonc_input%print_error_message()
                    call simlog%output('Error: building unit data are corrupted')
                    stop
                end if    
            end do             
            if (jsonc_input%failed()) then    
                call jsonc_input%print_error_message()
                call simlog%output('Error: building data are corrupted')
                stop
            end if    
        end do    
        ! neighbors
        call jsonf_input%info('neighbors', found, n_children=n_input_elements)
        if (.not.found) then
        	call simlog%output('Warning: structure named [neighbors] not found in input json file.')
        end if
        allocate(this%building_complex%neighbors(n_input_elements))
        do i=1,n_input_elements
            write (str_id, *) i
            str_name='neighbors('//trim(adjustl(str_id))//')' ! create name of the structure element to find, i.e. buildings(<i>)
            call jsonf_input%get(trim(str_name), p_input, found) ! now p point to the i-th element of the structure that we want to read
            this%building_complex%neighbors(i)%id=i
            call jsonc_input%get(p_input,'name',this%building_complex%neighbors(i)%name)  ! get the name of the structure element
            call jsonc_input%get(p_input, 'X', this%building_complex%neighbors(i)%X)
            call jsonc_input%get(p_input, 'Y', this%building_complex%neighbors(i)%Y)
            call jsonc_input%get(p_input,'height_ag',this%building_complex%neighbors(i)%height_ag)  
            if (jsonc_input%failed()) then    
                call jsonc_input%print_error_message()
                call simlog%output('Error: neighbors data are corrupted.')
                stop
            end if    
        end do
        ! thermal_grids
        call jsonf_input%info('thermal_grids', found, n_children=n_input_elements)
        if (.not.found) then
            call simlog%output('Warning: structure named [thermal_grids] not found in input json file.')
        end if
        allocate(this%thermal_grids(n_input_elements))
        do i=1,n_input_elements
            write (str_id, *) i
            str_name='thermal_grids('//trim(adjustl(str_id))//')' ! create name of the structure element to find, i.e. buildings(<i>)
            call jsonf_input%get(trim(str_name), p_input, found) ! now p point to the i-th element of the structure that we want to read
            ! get number of double_connections
            call jsonc_input%info(p_input,'double_connections', found, n_children=n_input_units)    
            allocate(this%thermal_grids(i)%double_connections(n_input_units))
            do j=1,n_input_units
                write (unit_id, *) j
                str_unit_name=trim(str_name)//'.double_connections('//trim(adjustl(unit_id))//')' ! create name of the structure element to find, i.e. units(<j>)                
				this%thermal_grids(i)%double_connections(j)%p_thermal_grid=>this%thermal_grids(i)
                call jsonf_input%get(trim(str_unit_name), child_input, found) ! now p point to the i-th element of the structure that we want to read
                call jsonc_input%get(child_input,'name',this%thermal_grids(i)%double_connections(j)%name)  
                call jsonc_input%get(child_input,'P0',vec_dbl)
                this%thermal_grids(i)%double_connections(j)%P0=point(x=vec_dbl(1),y=vec_dbl(2),z=vec_dbl(3))  			
			end do		
            ! get number of double_pipes
            call jsonc_input%info(p_input,'double_pipes', found, n_children=n_input_units)    
            allocate(this%thermal_grids(i)%double_pipes(n_input_units))
            do j=1,n_input_units
                write (unit_id, *) j
                str_unit_name=trim(str_name)//'.double_pipes('//trim(adjustl(unit_id))//')' ! create name of the structure element to find, i.e. units(<j>)                
				this%thermal_grids(i)%double_pipes(j)%p_thermal_grid=>this%thermal_grids(i)
                call jsonf_input%get(trim(str_unit_name), child_input, found) ! now p point to the i-th element of the structure that we want to read
                call jsonc_input%get(child_input,'name',this%thermal_grids(i)%double_pipes(j)%name) 
				call jsonc_input%get(child_input,'D_int',this%thermal_grids(i)%double_pipes(j)%D_int)
				call jsonc_input%get(child_input,'ins_lambda',this%thermal_grids(i)%double_pipes(j)%ins_lambda)
				call jsonc_input%get(child_input,'ins_delta',this%thermal_grids(i)%double_pipes(j)%ins_delta)
				call jsonc_input%get(child_input,'eps',this%thermal_grids(i)%double_pipes(j)%eps)   		 	
				call jsonc_input%get(child_input,'K',this%thermal_grids(i)%double_pipes(j)%K)  	
				! get number of path points
				call jsonc_input%info(child_input,'path', found, n_children=n_input_vertices)
				allocate(this%thermal_grids(i)%double_pipes(j)%path%points(n_input_vertices))
				do k=1,n_input_vertices
					write (unit_id, *) k
					str_unit_name='path('//trim(adjustl(unit_id))//')' ! create name of the structure element to find, i.e. units(<j>)                
					call jsonc_input%get(child_input,str_unit_name,vec_dbl)
					this%thermal_grids(i)%double_pipes(j)%path%points(k)=point(x=vec_dbl(1),y=vec_dbl(2),z=vec_dbl(3))
                end do	
                call this%thermal_grids(i)%double_pipes(j)%path%set_properties()
			end do		
            if (jsonc_input%failed()) then    
                call jsonc_input%print_error_message()
                call simlog%output('Error: pipes or connections data are corrupted.')
                stop
            end if    
            ! get number of substations
			call jsonc_input%info(p_input,'substations', found, n_children=n_input_units)    
			allocate(this%thermal_grids(i)%substations(n_input_units))
            do j=1,n_input_units
                write (unit_id, *) j
                str_unit_name=trim(str_name)//'.substations('//trim(adjustl(unit_id))//')' ! create name of the structure element to find, i.e. units(<j>)                
                call jsonf_input%get(trim(str_unit_name), child_input, found) ! now p point to the i-th element of the structure that we want to read
                call jsonc_input%get(child_input,'type_sub',tmp_str) 
                call substation_wrapper_set(this%thermal_grids(i)%substations(j),tmp_str) ! allocate item based on type_sub
				this%thermal_grids(i)%substations(j)%item%type_sub=tmp_str
				this%thermal_grids(i)%substations(j)%item%p_thermal_grid=>this%thermal_grids(i)
				this%thermal_grids(i)%substations(j)%item%p_building_complex=>this%building_complex
                call jsonc_input%get(child_input,'name',this%thermal_grids(i)%substations(j)%item%name) 
                call jsonc_input%get(child_input,'building_units',tmp_str) 
				call buffer%splitsep(string(tmp_str),",")
				allocate(this%thermal_grids(i)%substations(j)%item%building_units(buffer%num_tokens))
				do k=1,buffer%num_tokens
					this%thermal_grids(i)%substations(j)%item%building_units(k)=buffer%tokens(k)
				end do
                call jsonc_input%get(child_input,'P0',vec_dbl)
                this%thermal_grids(i)%substations(j)%item%P0=point(x=vec_dbl(1),y=vec_dbl(2),z=vec_dbl(3))  			
				this%thermal_grids(i)%substations(j)%item%DHW_tank=dhw_storage()
				! get number of parameters
				call jsonc_input%info(child_input,'parameters', found, n_children=n_input_childs)
				allocate(this%thermal_grids(i)%substations(j)%item%parameters(n_input_childs))
				do k=1,n_input_childs
					write (unit_id, *) k
					str_par_name=trim(str_unit_name)//'.parameters('//trim(adjustl(unit_id))//')' ! create name of the structure element to find, i.e. units(<j>)                
					call jsonf_input%get(trim(str_par_name), child_input, found) ! now p point to the i-th element of the structure that we want to read
					call jsonc_input%get(child_input,'name',this%thermal_grids(i)%substations(j)%item%parameters(k)%name%text)
                    call jsonc_input%info(child_input,'val_int',found)
                    if (found) call jsonc_input%get(child_input,'val_int',this%thermal_grids(i)%substations(j)%item%parameters(k)%val_int)
                    call jsonc_input%info(child_input,'val_dbl',found)
					if (found) call jsonc_input%get(child_input,'val_dbl',this%thermal_grids(i)%substations(j)%item%parameters(k)%val_dbl)
                    call jsonc_input%info(child_input,'val_str',found)
					if (found) call jsonc_input%get(child_input,'val_str',this%thermal_grids(i)%substations(j)%item%parameters(k)%val_str)
                    call jsonc_input%info(child_input,'arr_dbl',found)
					if (found) call jsonc_input%get(child_input,'arr_dbl',this%thermal_grids(i)%substations(j)%item%parameters(k)%arr_dbl)
                end do	
				call this%thermal_grids(i)%substations(j)%item%load() ! call the type-specific method to load parameters in the substation's objects
			end do	
            ! get number of energy_centres
			call jsonc_input%info(p_input,'energy_centres', found, n_children=n_input_units)    
			allocate(this%thermal_grids(i)%energy_centres(n_input_units))
            do j=1,n_input_units
                write (unit_id, *) j
                str_unit_name=trim(str_name)//'.energy_centres('//trim(adjustl(unit_id))//')' ! create name of the structure element to find, i.e. units(<j>)                
                call jsonf_input%get(trim(str_unit_name), child_input, found) ! now p point to the i-th element of the structure that we want to read
                call jsonc_input%get(child_input,'type_centre',tmp_str) 
                call energy_centre_wrapper_set(this%thermal_grids(i)%energy_centres(j),tmp_str) ! allocate item based on type_centre
				this%thermal_grids(i)%energy_centres(j)%item%type_centre=tmp_str
				this%thermal_grids(i)%energy_centres(j)%item%p_thermal_grid=>this%thermal_grids(i)
				this%thermal_grids(i)%energy_centres(j)%item%p_building_complex=>this%building_complex                
                call jsonc_input%get(child_input,'name',this%thermal_grids(i)%energy_centres(j)%item%name) 
                call jsonc_input%get(child_input,'P0',vec_dbl)
                this%thermal_grids(i)%energy_centres(j)%item%P0=point(x=vec_dbl(1),y=vec_dbl(2),z=vec_dbl(3))  			
				! get number of parameters
				call jsonc_input%info(child_input,'parameters', found, n_children=n_input_childs)
				allocate(this%thermal_grids(i)%energy_centres(j)%item%parameters(n_input_childs))
				do k=1,n_input_childs
					write (unit_id, *) k
					str_par_name=trim(str_unit_name)//'.parameters('//trim(adjustl(unit_id))//')' ! create name of the structure element to find, i.e. units(<j>)                
					call jsonf_input%get(trim(str_par_name), child_input, found) ! now p point to the i-th element of the structure that we want to read
					call jsonc_input%get(child_input,'name',this%thermal_grids(i)%energy_centres(j)%item%parameters(k)%name%text)
                    call jsonc_input%info(child_input,'val_int',found)
                    if (found) call jsonc_input%get(child_input,'val_int',this%thermal_grids(i)%energy_centres(j)%item%parameters(k)%val_int)
                    call jsonc_input%info(child_input,'val_dbl',found)
					if (found) call jsonc_input%get(child_input,'val_dbl',this%thermal_grids(i)%energy_centres(j)%item%parameters(k)%val_dbl)
                    call jsonc_input%info(child_input,'val_str',found)
					if (found) call jsonc_input%get(child_input,'val_str',this%thermal_grids(i)%energy_centres(j)%item%parameters(k)%val_str)
                    call jsonc_input%info(child_input,'arr_dbl',found)
					if (found) call jsonc_input%get(child_input,'arr_dbl',this%thermal_grids(i)%energy_centres(j)%item%parameters(k)%arr_dbl)
                    call jsonc_input%info(child_input,'mat_dbl',found)
					if (found) then
						call jsonc_input%get(child_input,'mat_dbl',jmat)
                        nrow=jsonc_input%count(jmat)
						if (nrow<=0) then
							allocate(this%thermal_grids(i)%energy_centres(j)%item%parameters(k)%mat_dbl(0, 0))							
						else
							call jsonc_input%get_child(jmat,jrow)
							ncol=jsonc_input%count(jrow)
							allocate(this%thermal_grids(i)%energy_centres(j)%item%parameters(k)%mat_dbl(nrow, ncol))
							do l = 1, nrow
							    call jsonc_input%get(jrow, vec_dbl)
								this%thermal_grids(i)%energy_centres(j)%item%parameters(k)%mat_dbl(l, :)=vec_dbl
								call jsonc_input%get_next(jrow, jnext)
								if (.not. associated(jnext)) exit
								jrow => jnext							
							end do
						end if	
					end if	
                end do	
				call this%thermal_grids(i)%energy_centres(j)%item%load() ! call the type-specific method to load parameters in the substation's objects
			end do	
        end do
        ! cleaning
        call jsonf_input%destroy()
        if (jsonf_input%failed()) then
            call simlog%output('Error: unable to destroy input json file')
            stop
        end if
        call jsonc_input%destroy()
        if (jsonc_input%failed()) then
        	call simlog%output('Error: unable to destroy json object')
            stop
        end if
  
        ! =============================================================================================
        ! load output file in memory
        ! =============================================================================================
  !      call jsonf_input%initialize()
		!inquire(file=file_output_name,exist=found)
		!if (found) then
		!	call jsonf_input%load(file_output_name)
		!	if (jsonf_input%failed()) then    
		!		call jsonf_input%print_error_message()
		!		stop
		!	end if
		!else
		!	return
		!end if	
  !      call jsonc_input%initialize()
		!! get file format
  !      call jsonf_input%get('binary_file', p_input,found)
		!if (found) then
		!	call jsonc_input%get(p_input,fn%output_file_name_bin)
		!	is_binary=.true.
		!else
		!	is_binary=.false.
		!end if
		!
  !      ! thermal_grids
  !      call jsonf_input%info('thermal_grids', found, n_children=n_input_elements)
  !      if (.not.found) then
  !      	call simlog%output('Warning: structure named [thermal_grids] not found in output json file.')
		!else
  !          n_input_elements=min(size(this%thermal_grids),n_input_elements)
  !          do i=1,n_input_elements
  !              write (str_id, *) i
  !              str_name='thermal_grids('//trim(adjustl(str_id))//')' ! create name of the structure element to find, i.e. buildings(<i>)
  !              call jsonf_input%get(trim(str_name), p_input, found) ! now p point to the i-th element of the structure that we want to read
		!		! to do ... restart not working for thermal grids, for now it is skipped
		!		!call this%thermal_grids(i)%from_json(jsonc_input,p_input)
		!	end do
  !      end if
	 !
  !      ! buildings
  !      call jsonf_input%info('buildings', found, n_children=n_input_elements)
  !      if (.not.found) then
  !      	call simlog%output('Error: structure named [buildings] not found in output json file.')
  !          stop
  !      else    
  !          n_input_elements=min(size(this%building_complex%buildings),n_input_elements)
  !          do i=1,n_input_elements
  !              write (str_id, *) i
  !              str_name='buildings('//trim(adjustl(str_id))//')' ! create name of the structure element to find, i.e. buildings(<i>)
  !              call jsonf_input%get(trim(str_name), p_input, found) ! now p point to the i-th element of the structure that we want to read
  !              ! we presume no changes in building data, otherwise the user must set in the scheduler the building calculations from scratch
  !              ! get number of footprint_segments
  !              call jsonc_input%info(p_input,'footprint_segments', found, n_children=n_input_segments)  
  !              if (n_input_segments>0) then
  !                  allocate(this%building_complex%buildings(i)%footprint_segments(n_input_segments))
  !                  do j=1,n_input_segments
  !                      write (segment_id, *) j
  !                      call jsonc_input%get(p_input,'footprint_segments.['//trim(adjustl(segment_id))//'].X1',this%building_complex%buildings(i)%footprint_segments(j)%X1)
  !                      call jsonc_input%get(p_input,'footprint_segments.['//trim(adjustl(segment_id))//'].Y1',this%building_complex%buildings(i)%footprint_segments(j)%Y1)
  !                      call jsonc_input%get(p_input,'footprint_segments.['//trim(adjustl(segment_id))//'].X2',this%building_complex%buildings(i)%footprint_segments(j)%X2)
  !                      call jsonc_input%get(p_input,'footprint_segments.['//trim(adjustl(segment_id))//'].Y2',this%building_complex%buildings(i)%footprint_segments(j)%Y2)
  !                      call jsonc_input%get(p_input,'footprint_segments.['//trim(adjustl(segment_id))//'].azimuth',this%building_complex%buildings(i)%footprint_segments(j)%azimuth)
  !                      call jsonc_input%get(p_input,'footprint_segments.['//trim(adjustl(segment_id))//'].length',this%building_complex%buildings(i)%footprint_segments(j)%length)
  !                  end do    
  !              end if
  !              ! get number of units
  !              call jsonc_input%info(p_input,'units', found, n_children=n_input_units)    
  !              n_input_units=min(size(this%building_complex%buildings(i)%units),n_input_units)
  !              do j=1,n_input_units
  !                  write (unit_id, *) j
  !                  str_unit_name=trim(str_name)//'.units('//trim(adjustl(unit_id))//')' ! create name of the structure element to find, i.e. units(<j>)                
  !                  call jsonf_input%get(trim(str_unit_name), p_input, found) ! now p point to the i-th element of the structure that we want to read
  !                  ! get materials
  !                  call jsonc_input%info(p_input,'materials', found, n_children=n_input_units_materials)  
  !                  call this%building_complex%buildings(i)%units(j)%materials_inventory%initialize()
  !                  do k=1,n_input_units_materials
  !                      write (material_id, *) k
  !                      str_unit_material_name=trim(str_name)//'.units('//trim(adjustl(unit_id))//')'//'.materials('//trim(adjustl(material_id))//')' ! create name of the structure element to find, i.e. units(<j>)                
  !                      call jsonf_input%get(trim(str_unit_material_name), p_input, found) ! now p point to the i-th element of the structure     
  !                      call jsonc_input%get(p_input,'material_id',material_name)
  !                      call jsonc_input%get(p_input,'material_mass',material_mass)
  !                      !print *, 'add material to inventory list...',material_name,material_mass
  !                      call this%building_complex%buildings(i)%units(j)%materials_inventory%update(material_name,material_mass)
  !                      if (jsonc_input%failed()) then    
  !                          call jsonc_input%print_error_message()
  !                          call simlog%output('Error: materials not found in output file.')
  !                          stop
  !                      end if    
  !                  end do    
  !                  call jsonf_input%get(trim(str_unit_name), p_input, found) ! now p point to the i-th element of the structure that we want to read
  !                  ! get dhw 
  !                  call jsonc_input%get(p_input,'dhw.V_storage',this%building_complex%buildings(i)%units(j)%dhw%V_storage)
  !                  call jsonc_input%get(p_input,'dhw.UA_storage',this%building_complex%buildings(i)%units(j)%dhw%UA_storage)
  !                  call jsonc_input%get(p_input,'dhw.Q_heater',this%building_complex%buildings(i)%units(j)%dhw%Q_heater)
		!			if (is_binary) then
		!				call read_series_bin(jsonc_input,p_input,'dhw.Q_dhw',this%building_complex%buildings(i)%units(j)%dhw%Q_dhw)
		!				call read_series_bin(jsonc_input,p_input,'dhw.mdotcp_dhw',this%building_complex%buildings(i)%units(j)%dhw%mdotcp_dhw)
		!				call read_series_bin(jsonc_input,p_input,'dhw.T_main',this%building_complex%buildings(i)%units(j)%dhw%T_main)
		!			else
		!				call jsonc_input%get(p_input,'dhw.Q_dhw',vec_dbl);this%building_complex%buildings(i)%units(j)%dhw%Q_dhw(1:8760)=vec_dbl(1:8760)
		!				call jsonc_input%get(p_input,'dhw.mdotcp_dhw',vec_dbl);this%building_complex%buildings(i)%units(j)%dhw%mdotcp_dhw(1:8760)=vec_dbl(1:8760)
		!				call jsonc_input%get(p_input,'dhw.T_main',vec_dbl);this%building_complex%buildings(i)%units(j)%dhw%T_main(1:8760)=vec_dbl(1:8760)
		!			end if	
  !                  ! get hvac
  !                  call jsonc_input%get(p_input,'hvac.T_supply_hw',this%building_complex%buildings(i)%units(j)%hvac%T_supply_hw)
  !                  call jsonc_input%get(p_input,'hvac.T_supply_cw',this%building_complex%buildings(i)%units(j)%hvac%T_supply_cw)
		!			if (is_binary) then
		!				call read_series_bin(jsonc_input,p_input,'hvac.Q_heat_hyd',this%building_complex%buildings(i)%units(j)%hvac%Q_heat_hyd)
		!				call read_series_bin(jsonc_input,p_input,'hvac.Q_cool_hyd',this%building_complex%buildings(i)%units(j)%hvac%Q_cool_hyd)
		!				call read_series_bin(jsonc_input,p_input,'hvac.Q_heat_AHU',this%building_complex%buildings(i)%units(j)%hvac%Q_heat_AHU)
		!				call read_series_bin(jsonc_input,p_input,'hvac.Q_cool_AHU',this%building_complex%buildings(i)%units(j)%hvac%Q_cool_AHU)
		!				call read_series_bin(jsonc_input,p_input,'hvac.W_par_hyd',this%building_complex%buildings(i)%units(j)%hvac%W_par_hyd)
		!				call read_series_bin(jsonc_input,p_input,'hvac.W_par_AHU',this%building_complex%buildings(i)%units(j)%hvac%W_par_AHU)
		!			else
		!				call jsonc_input%get(p_input,'hvac.Q_heat_hyd',vec_dbl);this%building_complex%buildings(i)%units(j)%hvac%Q_heat_hyd(1:8760)=vec_dbl(1:8760)
		!				call jsonc_input%get(p_input,'hvac.Q_cool_hyd',vec_dbl);this%building_complex%buildings(i)%units(j)%hvac%Q_cool_hyd(1:8760)=vec_dbl(1:8760)
		!				call jsonc_input%get(p_input,'hvac.Q_heat_AHU',vec_dbl);this%building_complex%buildings(i)%units(j)%hvac%Q_heat_AHU(1:8760)=vec_dbl(1:8760)
		!				call jsonc_input%get(p_input,'hvac.Q_cool_AHU',vec_dbl);this%building_complex%buildings(i)%units(j)%hvac%Q_cool_AHU(1:8760)=vec_dbl(1:8760)
		!				call jsonc_input%get(p_input,'hvac.W_par_hyd',vec_dbl);this%building_complex%buildings(i)%units(j)%hvac%W_par_hyd(1:8760)=vec_dbl(1:8760)
		!				call jsonc_input%get(p_input,'hvac.W_par_AHU',vec_dbl);this%building_complex%buildings(i)%units(j)%hvac%W_par_AHU(1:8760)=vec_dbl(1:8760)      
		!			end if		
  !                  ! get geometry
  !                  call jsonc_input%get(p_input,'geometry.footprint_area',this%building_complex%buildings(i)%units(j)%geometry%footprint_area)
  !                  call jsonc_input%get(p_input,'geometry.footprint_azimuth',this%building_complex%buildings(i)%units(j)%geometry%footprint_azimuth)
  !                  call jsonc_input%get(p_input,'geometry.floor_area_gross',this%building_complex%buildings(i)%units(j)%geometry%floor_area_gross)
  !                  call jsonc_input%get(p_input,'geometry.floor_area_heated',this%building_complex%buildings(i)%units(j)%geometry%floor_area_heated)
  !                  call jsonc_input%get(p_input,'geometry.floor_area_useful',this%building_complex%buildings(i)%units(j)%geometry%floor_area_useful)
  !                  call jsonc_input%get(p_input,'geometry.volume_gross',this%building_complex%buildings(i)%units(j)%geometry%volume_gross)
  !                  call jsonc_input%get(p_input,'geometry.volume_useful',this%building_complex%buildings(i)%units(j)%geometry%volume_useful)
  !                  call jsonc_input%get(p_input,'geometry.window_area',this%building_complex%buildings(i)%units(j)%geometry%window_area) 
  !                  call jsonc_input%get(p_input,'geometry.number_of_storeys_above_ground',this%building_complex%buildings(i)%units(j)%geometry%number_of_storeys_above_ground) 
  !                  call jsonc_input%get(p_input,'geometry.number_of_storeys_below_ground',this%building_complex%buildings(i)%units(j)%geometry%number_of_storeys_below_ground) 
  !                  call jsonc_input%get(p_input,'geometry.number_of_storeys',this%building_complex%buildings(i)%units(j)%geometry%number_of_storeys) 
  !                  call jsonc_input%get(p_input,'geometry.height_above_ground',this%building_complex%buildings(i)%units(j)%geometry%height_above_ground) 
  !                  call jsonc_input%get(p_input,'geometry.height_below_ground',this%building_complex%buildings(i)%units(j)%geometry%height_below_ground) 
  !                  call jsonc_input%get(p_input,'geometry.height',this%building_complex%buildings(i)%units(j)%geometry%height) 
  !                  call jsonc_input%get(p_input,'geometry.height_floor_start',this%building_complex%buildings(i)%units(j)%geometry%height_floor_start)
  !                  ! get thermal_zone
  !                  call jsonc_input%get(p_input,'thermal_zone.air_mass',this%building_complex%buildings(i)%units(j)%thermal_zone%air_mass)
  !                  call jsonc_input%get(p_input,'thermal_zone.internal_capacity',this%building_complex%buildings(i)%units(j)%thermal_zone%internal_capacity) 
  !                  call jsonc_input%get(p_input,'thermal_zone.envelope_area',this%building_complex%buildings(i)%units(j)%thermal_zone%envelope_area) 
  !                  call jsonc_input%get(p_input,'thermal_zone.envelope_HT',this%building_complex%buildings(i)%units(j)%thermal_zone%envelope_HT) 
		!			if (is_binary) then
		!				call read_series_bin(jsonc_input,p_input,'thermal_zone.T_air',this%building_complex%buildings(i)%units(j)%thermal_zone%T_air)
		!				call read_series_bin(jsonc_input,p_input,'thermal_zone.X_air',this%building_complex%buildings(i)%units(j)%thermal_zone%X_air)
		!				call read_series_bin(jsonc_input,p_input,'thermal_zone.T_op',this%building_complex%buildings(i)%units(j)%thermal_zone%T_op)
		!				call read_series_bin(jsonc_input,p_input,'thermal_zone.Q_sh_sen',this%building_complex%buildings(i)%units(j)%thermal_zone%Q_sh_sen)
		!				call read_series_bin(jsonc_input,p_input,'thermal_zone.Q_sc_sen',this%building_complex%buildings(i)%units(j)%thermal_zone%Q_sc_sen)
		!				call read_series_bin(jsonc_input,p_input,'thermal_zone.Q_sh_lat',this%building_complex%buildings(i)%units(j)%thermal_zone%Q_sh_lat)
		!				call read_series_bin(jsonc_input,p_input,'thermal_zone.Q_sc_lat',this%building_complex%buildings(i)%units(j)%thermal_zone%Q_sc_lat)
		!			else
		!				call jsonc_input%get(p_input,'thermal_zone.T_air',vec_dbl);this%building_complex%buildings(i)%units(j)%thermal_zone%T_air(1:8760)=vec_dbl(1:8760) 
		!				call jsonc_input%get(p_input,'thermal_zone.X_air',vec_dbl);this%building_complex%buildings(i)%units(j)%thermal_zone%X_air(1:8760)=vec_dbl(1:8760) 
		!				call jsonc_input%get(p_input,'thermal_zone.T_op',vec_dbl);this%building_complex%buildings(i)%units(j)%thermal_zone%T_op(1:8760)=vec_dbl(1:8760) 
		!				call jsonc_input%get(p_input,'thermal_zone.Q_sh_sen',vec_dbl);this%building_complex%buildings(i)%units(j)%thermal_zone%Q_sh_sen(1:8760)=vec_dbl(1:8760) 
		!				call jsonc_input%get(p_input,'thermal_zone.Q_sc_sen',vec_dbl);this%building_complex%buildings(i)%units(j)%thermal_zone%Q_sc_sen(1:8760)=vec_dbl(1:8760) 
		!				call jsonc_input%get(p_input,'thermal_zone.Q_sh_lat',vec_dbl);this%building_complex%buildings(i)%units(j)%thermal_zone%Q_sh_lat(1:8760)=vec_dbl(1:8760)
		!				call jsonc_input%get(p_input,'thermal_zone.Q_sc_lat',vec_dbl);this%building_complex%buildings(i)%units(j)%thermal_zone%Q_sc_lat(1:8760)=vec_dbl(1:8760)						 
		!			end if	
  !                  ! get windows
  !                  call jsonc_input%info(p_input,'thermal_zone.windows', found, n_children=n_input_windows)
  !                  if (n_input_windows>0) then
  !                      allocate(this%building_complex%buildings(i)%units(j)%thermal_zone%windows(n_input_windows))
  !                      do l=1,n_input_windows
  !                          write (l_id, *) l
  !                          call jsonc_input%get(p_input,'thermal_zone.windows['//trim(adjustl(l_id))//'].window_name',window_name)
  !                          this%building_complex%buildings(i)%units(j)%thermal_zone%windows(l)%struct=>this%building_complex%get_window(window_name)
  !                          call jsonc_input%get(p_input,'thermal_zone.windows['//trim(adjustl(l_id))//'].area',this%building_complex%buildings(i)%units(j)%thermal_zone%windows(l)%area)
  !                          call jsonc_input%get(p_input,'thermal_zone.windows['//trim(adjustl(l_id))//'].azimuth',this%building_complex%buildings(i)%units(j)%thermal_zone%windows(l)%azimuth)
  !                          call jsonc_input%get(p_input,'thermal_zone.windows['//trim(adjustl(l_id))//'].tilt',this%building_complex%buildings(i)%units(j)%thermal_zone%windows(l)%tilt)
  !                      end do    
  !                  end if                        
  !                  ! get external_structures
  !                  call jsonc_input%info(p_input,'thermal_zone.external_structures', found, n_children=n_input_external_structures)  
  !                  allocate(this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(n_input_external_structures))
  !                  do k=1,n_input_external_structures
  !                      write (structure_id, *) k
  !                      str_structure_name=trim(str_unit_name)//'.thermal_zone.external_structures('//trim(adjustl(structure_id))//')' ! create name of the structure element to find, i.e. units(<j>)                
  !                      call jsonf_input%get(trim(str_structure_name), p_input, found) ! now p point to the i-th element of the structure 
  !                      call jsonc_input%get(p_input,'structure_name',structure_name)
  !                      this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%struct=>this%building_complex%get_wall_slab(structure_name)
  !                      call jsonc_input%get(p_input,'area',this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%area)
  !                      call jsonc_input%get(p_input,'area_above_ground',this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%area_above_ground)
  !                      call jsonc_input%get(p_input,'area_below_ground',this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%area_below_ground)
  !                      call jsonc_input%get(p_input,'azimuth',this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%azimuth)
  !                      call jsonc_input%get(p_input,'tilt',this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%tilt)
  !                      call jsonc_input%get(p_input,'is_adiabatic',this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%is_adiabatic)
  !                      call jsonc_input%get(p_input,'Htb',this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%Htb)
  !                      call jsonc_input%get(p_input,'resistance',this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%resistance)
  !                      call jsonc_input%get(p_input,'transmittance',this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%transmittance)
  !                      call jsonc_input%info(p_input,'vertices', found, n_children=n_input_vertices)
  !                      if (n_input_vertices>0) then
  !                          allocate(this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%vertices(n_input_vertices))
  !                          do l=1,n_input_vertices
  !                              write (l_id, *) l
  !                              call jsonc_input%get(p_input,'vertices['//trim(adjustl(l_id))//'].x',this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%vertices(l)%x)
  !                              call jsonc_input%get(p_input,'vertices['//trim(adjustl(l_id))//'].y',this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%vertices(l)%y)
  !                              call jsonc_input%get(p_input,'vertices['//trim(adjustl(l_id))//'].z',this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%vertices(l)%z)
  !                          end do    
  !                      end if    
  !                      call jsonc_input%info(p_input,'grid', found, n_children=n_input_grid)
  !                      if (n_input_grid>0) then
  !                          allocate(this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%grid(n_input_grid))
  !                          do l=1,n_input_grid
  !                              write (l_id, *) l
  !                              call jsonc_input%get(p_input,'grid['//trim(adjustl(l_id))//'].x',this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%grid(l)%x)
  !                              call jsonc_input%get(p_input,'grid['//trim(adjustl(l_id))//'].y',this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%grid(l)%y)
  !                              call jsonc_input%get(p_input,'grid['//trim(adjustl(l_id))//'].z',this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%grid(l)%z)
  !                          end do    
  !                      end if    
  !                      call jsonc_input%info(p_input,'fsol', found, n_children=n_input_fsol_1)
  !                      if (n_input_fsol_1>0) then
  !                          call jsonc_input%info(p_input,'fsol[1]', found, n_children=n_input_fsol_2)
  !                          call jsonc_input%info(p_input,'fsol[1][1]', found, n_children=n_input_fsol_3)
  !                          allocate(this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%fsol(n_input_fsol_1,n_input_fsol_2,n_input_fsol_3))
  !                          do l=1,n_input_fsol_1
  !                              do m=1,n_input_fsol_2
  !                                  write (l_id, *) l
  !                                  write (m_id, *) m
  !                                  call jsonc_input%get(p_input,'fsol['//trim(adjustl(l_id))//']['//trim(adjustl(m_id))//']',vec_dbl)
  !                                  this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%fsol(l,m,1:n_input_fsol_3)=vec_dbl(1:n_input_fsol_3)
  !                              end do
  !                          end do    
  !                      end if    
  !                  end do                    
  !                  if (jsonc_input%failed()) then    
  !                      call jsonc_input%print_error_message()
  !                      call simlog%output('Error: thermal zones data are corrupted in output file.')
  !                      stop
  !                  end if    
  !              end do             
  !              if (jsonc_input%failed()) then    
  !                  call jsonc_input%print_error_message()
  !                  call simlog%output('Error: building data are corrupted in output file.')
  !                  stop
  !              end if    
  !          end do    
  !         
  !      end if
		                
    end subroutine    
    
!buildings[{}]
!   name    
!	footprint_segments[{}]
!		X1
!		Y1
!		X2
!		Y2
!       azimuth       
!       length        
!       is_adiabatic           
!	units[{}]
!       name
!		geometry{}
!			footprint_area    
!			footprint_azimuth 
!			floor_area_gross   
!			floor_area_heated 
!			floor_area_useful 
!			volume_gross      
!			volume_useful    
!			window_area      
!			number_of_storeys_above_ground 
!			number_of_storeys_below_ground 
!			number_of_storeys          
!			height_above_ground 
!			height_below_ground 
!			height            
!			height_floor_start				
!		thermal_zone{}
!			air_mass 
!			internal_capacity  
!			envelope_area 
!			envelope_HT 
!			external_structures[] 
!				struct->name
!				struct->transmittance
!               struct->resistance
!				area
!				area_above_ground
!				area_below_ground
!				azimuth
!				tilt
!        		is_adiabatic
!				Htb  
!				resistance 
!				transmittance        
!				vertices[{}] 
!					x
!					y
!					z
!				grid [{}]
!					x
!					y
!					z
!				fsol[] 
!			windows[{}]
!				area      
!				azimuth   
!				tilt  	
!				struct->name
!				struct->transmittance
!				struct->gw
!thermal_grids[{}]
!	fluid{}
!		rho
!		cp
!	double_pipes[{}]
!		pipe_A{}
!			mdot[8760]
!			thermal_nodes[{}]
!				T[8760]
!				Qin[8760]
!		pipe_B{}
!			mdot[8760]
!			thermal_nodes[{}]
!				T[8760]
!				Qin[8760]
!	double_connections[{}]
!		
!	substations[{}]
!		
!	energy_centres[{}]
!		

    subroutine district_output_json(this,scheduler_filename,output_filename)
        class(district)::this
        character(*)::scheduler_filename,output_filename
        integer::i,j,k,l,m
        integer::iunit
        type(material_list_element),pointer::material_start
        type(json_value),pointer::p_root,p_buildings,p_segments,p_units,p_structures,p_windows,p_grid,p_vertices,p_fsol_1,p_fsol_2,p_materials,p_scheduler,p_thermal_grids,p_obj
        type(json_core)::json  ! library for manipulating json_value pointers
        ! scheduler file
        call json%initialize()
        nullify(p_root)
        call json%create_object(p_root,scheduler_filename)
        nullify(p_scheduler)
        call json%create_array(p_scheduler,'scheduler')
        call json%add(p_root, p_scheduler)
        do i=1,size(this%scheduler)
            call add_data_to_scheduler(this%scheduler(i))
        end do    
        ! save output to file
        open(newunit=iunit, file=scheduler_filename, status='REPLACE')
        call json%print(p_root,iunit)
        if (json%failed()) then
            call json%print_error_message()
        end if
        close(iunit)        
        ! output file
        call json%initialize()
        nullify(p_root)
        call json%create_object(p_root,output_filename)
        nullify(p_buildings)
        ! time series data format
        open(newunit=iunit, file=fn%output_file_name_bin, status='REPLACE')
        close(iunit)
		call json%add(p_root, 'format_version', 1)
		call json%add(p_root, 'binary_file', fn%output_file_name_bin)
		call json%add(p_root, 'datatype', 'float64')
		call json%add(p_root, 'endianness', 'little')
        call this%building_complex%meteodata%to_json(json,p_root,'meteo')
        !call json%add(p_root,p_meteo)
        call json%create_array(p_buildings,'buildings')
        call json%add(p_root, p_buildings)
        do i=1,size(this%building_complex%buildings)
            nullify(p_segments)
            call json%create_array(p_segments,'footprint_segments')
            do j=1,size(this%building_complex%buildings(i)%footprint_segments)
                call add_data_to_segment(this%building_complex%buildings(i)%footprint_segments(j))
            end do  
            nullify(p_units)
            call json%create_array(p_units,'units')    
            do j=1,size(this%building_complex%buildings(i)%units)
                nullify(p_materials)
                call json%create_array(p_materials,'materials')
                nullify(p_structures)
                call json%create_array(p_structures,'external_structures')
                do k=1,size(this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures)
                    nullify(p_grid)
                    call json%create_array(p_grid,'grid')
                    do l=1,size(this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%grid)
                        call add_data_to_grid(this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%grid(l))
                    end do
                    nullify(p_vertices)
                    call json%create_array(p_vertices,'vertices')
                    do l=1,size(this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%vertices)
                        call add_data_to_vertices(this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%vertices(l))
                    end do
                    nullify(p_fsol_1)
                    call json%create_array(p_fsol_1,'fsol')
                    do l=1,size(this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%fsol,1)
                        nullify(p_fsol_2)
                        call json%create_array(p_fsol_2,'')
                        do m=1,size(this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%fsol,2)
                            call json%add(p_fsol_2,'',this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k)%fsol(l,m,:))
                        end do    
                        call json%add(p_fsol_1,p_fsol_2)
                        nullify(p_fsol_2)
                    end do                    
                    call add_data_to_structure(this%building_complex%buildings(i)%units(j)%thermal_zone%external_structures(k))
                    nullify(p_grid)
                    nullify(p_vertices)
                    nullify(p_fsol_1)
                end do 
                nullify(p_windows)
                call json%create_array(p_windows,'windows')
                do k=1,size(this%building_complex%buildings(i)%units(j)%thermal_zone%windows)
                    call add_data_to_window(this%building_complex%buildings(i)%units(j)%thermal_zone%windows(k))
                end do 
                material_start=>this%building_complex%buildings(i)%units(j)%materials_inventory%first
                if (associated(material_start)) then
                    do
                        call add_data_to_material(material_start%value)
                        material_start=>material_start%next
                        if (.not.associated(material_start)) exit
                    end do  
                end if
                call add_data_to_unit(this%building_complex%buildings(i)%units(j))
                nullify(p_materials)
                nullify(p_structures)
                nullify(p_windows)
            end do  
            call add_data_to_building(this%building_complex%buildings(i))
            nullify(p_segments)
            nullify(p_units)
        end do  
        nullify(p_buildings)
        if (size(this%thermal_grids)>0) then
            call json%create_array(p_thermal_grids,'thermal_grids')
		    call this%thermal_grids(1)%to_json(json, p_obj)
		    call json%add(p_thermal_grids,p_obj)
            call json%add(p_root, p_thermal_grids)
            nullify(p_thermal_grids)
        end if    
        ! save output to file
        open(newunit=iunit, file=output_filename, status='REPLACE')
        call json%print(p_root,iunit)
        if (json%failed()) then
            call json%print_error_message()
        end if
        close(iunit)
        nullify(p_root)
    
    contains
    	
       subroutine add_data_to_scheduler(s)
            type(json_value),pointer :: var
            type(calculation_step)::s
            !initialize:
            nullify(var)
            !create the object before data can be added:
            call json%create_object(var,'')    !name does not matter
            !variable info:
            call json%add(var,'id',s%id)
            call json%add(var,'step',s%step)
            if (s%execute) then
                call json%add(var,'timestamp',this%timestamp)
            else
                call json%add(var,'timestamp',s%timestamp)
            end if    
            call json%add(var,'execute',s%execute)
            call json%add(p_scheduler, var)
            nullify(var)
        end subroutine 
          
        subroutine add_data_to_building(b)
            type(json_value),pointer :: var
            type(building)::b
            !initialize:
            nullify(var)
            !create the object before data can be added:
            call json%create_object(var,'')    !name does not matter
            !variable info:
            call json%add(var,'name',b%name)
            call json%add(var,p_segments)
            call json%add(var,p_units)
            call json%add(p_buildings, var)
            nullify(var)
        end subroutine 
        
        subroutine add_data_to_segment(fseg)
            type(json_value),pointer :: var
            type(footprint_edge)::fseg
            !initialize:
            nullify(var)
            !create the object before data can be added:
            call json%create_object(var,'')    !name does not matter
            !variable info:
            call json%add(var,'X1',fseg%X1)
            call json%add(var,'Y1',fseg%Y1)
            call json%add(var,'X2',fseg%X2)
            call json%add(var,'Y2',fseg%Y2)
            call json%add(var,'azimuth',fseg%azimuth)
            call json%add(var,'length',fseg%length)
            call json%add(p_segments, var)
            nullify(var)
        end subroutine 
        
        subroutine add_data_to_unit(bunit)
            type(json_value),pointer :: var,geom,tzone,dhw,hvac
            type(building_unit)::bunit
            !initialize:
            nullify(var)
            !create the object before data can be added:
            call json%create_object(var,'')
            call json%create_object(geom,'geometry')
            call json%create_object(tzone,'thermal_zone')
            call json%create_object(dhw,'dhw')
            call json%create_object(hvac,'hvac')
            !variable info:		            
            call json%add(geom,'footprint_area',bunit%geometry%footprint_area)
            call json%add(geom,'footprint_azimuth',bunit%geometry%footprint_azimuth)
            call json%add(geom,'floor_area_gross',bunit%geometry%floor_area_gross)
            call json%add(geom,'floor_area_heated',bunit%geometry%floor_area_heated)
            call json%add(geom,'floor_area_useful',bunit%geometry%floor_area_useful)
            call json%add(geom,'volume_gross',bunit%geometry%volume_gross)
            call json%add(geom,'volume_useful',bunit%geometry%volume_useful)
            call json%add(geom,'window_area',bunit%geometry%window_area)
            call json%add(geom,'number_of_storeys_above_ground',bunit%geometry%number_of_storeys_above_ground)
            call json%add(geom,'number_of_storeys_below_ground',bunit%geometry%number_of_storeys_below_ground)
            call json%add(geom,'number_of_storeys',bunit%geometry%number_of_storeys)
            call json%add(geom,'height_above_ground',bunit%geometry%height_above_ground)
            call json%add(geom,'height_below_ground',bunit%geometry%height_below_ground)
            call json%add(geom,'height',bunit%geometry%height)
            call json%add(geom,'height_floor_start',bunit%geometry%height_floor_start)
            call json%add(tzone,'air_mass',bunit%thermal_zone%air_mass)
            call json%add(tzone,'internal_capacity',bunit%thermal_zone%internal_capacity)
            call json%add(tzone,'envelope_area',bunit%thermal_zone%envelope_area)
            call json%add(tzone,'envelope_HT',bunit%thermal_zone%envelope_HT) 
			call save_series_bin(json,tzone,'T_air',bunit%thermal_zone%T_air)			
            call save_series_bin(json,tzone,'X_air',bunit%thermal_zone%X_air)
            call save_series_bin(json,tzone,'T_op',bunit%thermal_zone%T_op)
            call save_series_bin(json,tzone,'Q_sh_sen',bunit%thermal_zone%Q_sh_sen)
            call save_series_bin(json,tzone,'Q_sc_sen',bunit%thermal_zone%Q_sc_sen)
            call save_series_bin(json,tzone,'Q_sh_lat',bunit%thermal_zone%Q_sh_lat)
            call save_series_bin(json,tzone,'Q_sc_lat',bunit%thermal_zone%Q_sc_lat)            
            call json%add(tzone,p_structures)
            call json%add(tzone,p_windows)
            call json%add(dhw,'V_storage',bunit%dhw%V_storage)
            call json%add(dhw,'UA_storage',bunit%dhw%UA_storage)
            call json%add(dhw,'Q_heater',bunit%dhw%Q_heater)
            call save_series_bin(json,dhw,'Q_dhw',bunit%dhw%Q_dhw)
            call save_series_bin(json,dhw,'mdotcp_dhw',bunit%dhw%mdotcp_dhw)
            call save_series_bin(json,dhw,'T_main',bunit%dhw%T_main)
            call json%add(var,'name',bunit%name)
            call json%add(hvac,'T_supply_hw',bunit%hvac%T_supply_hw)
            call json%add(hvac,'T_supply_cw',bunit%hvac%T_supply_cw)
            call save_series_bin(json,hvac,'Q_heat_hyd',bunit%hvac%Q_heat_hyd)
            call save_series_bin(json,hvac,'Q_cool_hyd',bunit%hvac%Q_cool_hyd)
            call save_series_bin(json,hvac,'Q_heat_AHU',bunit%hvac%Q_heat_AHU)
            call save_series_bin(json,hvac,'Q_cool_AHU',bunit%hvac%Q_cool_AHU)
            call save_series_bin(json,hvac,'W_par_hyd',bunit%hvac%W_par_hyd)
            call save_series_bin(json,hvac,'W_par_AHU',bunit%hvac%W_par_AHU)
            call json%add(var,p_materials)
            call json%add(var,geom)
            call json%add(var,tzone)
            call json%add(var,dhw)
            call json%add(var,hvac)
            call json%add(p_units, var)
            nullify(var)
        end subroutine   

       subroutine add_data_to_material(info)
            type(json_value),pointer :: var
            type(material_mass)::info
            !initialize:
            nullify(var)
            !create the object before data can be added:
            call json%create_object(var,'')    !name does not matter
            !variable info:
            call json%add(var,'material_id',info%name)
            call json%add(var,'material_mass',info%mass)
            call json%add(p_materials, var)
            nullify(var)
       end subroutine
                
       subroutine add_data_to_structure(info)
            type(json_value),pointer :: var
            type(wall_slab)::info
            !initialize:
            nullify(var)
            !create the object before data can be added:
            call json%create_object(var,'')    !name does not matter
            !variable info:
            call json%add(var,'structure_name',info%struct%name)
            call json%add(var,'structure_resistance',info%struct%resistance)
            call json%add(var,'structure_transmittance',info%struct%transmittance)
            call json%add(var,'area',info%area)
            call json%add(var,'area_above_ground',info%area_above_ground)
            call json%add(var,'area_below_ground',info%area_below_ground)
            call json%add(var,'azimuth',info%azimuth)
            call json%add(var,'tilt',info%tilt)
            call json%add(var,'is_adiabatic',info%is_adiabatic)
            call json%add(var,'Htb',info%Htb)
            call json%add(var,'resistance',info%resistance)
            call json%add(var,'transmittance',info%transmittance)
            call json%add(var,p_vertices)
            call json%add(var,p_grid)
            call json%add(var,p_fsol_1)
            call json%add(p_structures, var)
            nullify(var)
       end subroutine

       subroutine add_data_to_window(info)
            type(json_value),pointer :: var
            type(window)::info
            !initialize:
            nullify(var)
            !create the object before data can be added:
            call json%create_object(var,'')    !name does not matter
            !variable info:
            call json%add(var,'window_name',info%struct%name)
            call json%add(var,'window_resistance',info%struct%resistance)
            call json%add(var,'window_transmittance',info%struct%transmittance)
            call json%add(var,'area',info%area)
            call json%add(var,'azimuth',info%azimuth)
            call json%add(var,'tilt',info%tilt)
            call json%add(p_windows, var)
            nullify(var)
       end subroutine
       
        subroutine add_data_to_grid(p)
            type(point)::p
            type(json_value),pointer :: var
            !initialize:
            nullify(var)
            !create the object before data can be added:
            call json%create_object(var,'')    !name does not matter
            !variable info:
            call json%add(var, 'x',p%x)
            call json%add(var, 'y',p%y)
            call json%add(var, 'z',p%z)
            call json%add(p_grid, var)
            nullify(var)
        end subroutine

        subroutine add_data_to_vertices(p)
            type(point)::p
            type(json_value),pointer :: var
            !initialize:
            nullify(var)
            !create the object before data can be added:
            call json%create_object(var,'')    !name does not matter
            !variable info:
            call json%add(var, 'x',p%x)
            call json%add(var, 'y',p%y)
            call json%add(var, 'z',p%z)
            call json%add(p_vertices, var)
            nullify(var)
        end subroutine
        
    end subroutine

end module