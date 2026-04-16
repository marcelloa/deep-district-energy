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
! The module contains classes and methods related to substations and energy_centres as polymorphic objects derived from 
! respective core progenitor classes (substation, energy_centre).
! The developer can easily add new configurations for substations and energy centres by extending the progenitor classes.
!
module thermal_stations_module
	use thermal_grid_module

	! =============================================================================================================
	! substation classes 
	! =============================================================================================================

	! substation 001: HP(SH&DHW)+HX(SC)	
    type, extends(substation)::substation_001
		type(ww_heat_pump)::HP											! the water-water heat pump 
        type(pump)::pump_heat_pump	        							! heat pump's pump
        type(heat_exchanger)::HX_heat_pump              				! heat pump heat exchanger (evaporator)
		type(pipe)::pipe_heat_pump					    				! pipe from heat pump heat exchanger to network node B
        type(pump)::pump_cooling              							! cooling pump
        type(heat_exchanger)::HX_cooling                				! cooling heat exchanger 
		type(pipe)::pipe_cooling					    				! pipe from cooling heat exchanger to network node A
        double precision, dimension(8760)::Q_aux_h						! auxiliary heat at the heat pump heat exchanger
        double precision, dimension(8760)::Q_aux_c						! auxiliary cooling at the cooling heat exchanger
        double precision::T_grid_h_min=2								! lower grid tempertaure limit (°C) at the the heat pump heat exchanger
        double precision::T_grid_c_max=18								! upper grid tempertaure limit (°C) at the the cooling heat exchanger
    contains
        procedure::load=>substation_001_load 							! load parameters in the specific type objects
        procedure::create_networks=>substation_001_create_networks 		! create connections, pumps, heat exchangers, pipes 
        procedure::num_elements=>substation_001_num_elements 			! return the number of elements (see num_elements type) 
        procedure::initialize=>substation_001_initialize 				! get heating, cooling and dhw demands from the associated building units
        procedure::solve=>substation_001_solve 							! calculate the energy flows in the time step
		procedure::to_json=>substation_001_to_json						! generate json output
		procedure::from_json=>substation_001_from_json					! read json output
    end type    

	! substation 002: HP(SH&DHW)+HX(SC)+TANK for SH	(NOTE: for Linus case study)
    type, extends(substation)::substation_002
		type(ww_heat_pump)::HP											! the water-water heat pump 
		!type(water_storage)::SHS										! the storage tank for space heating (NOTE for Linus: create water_storage class, not existing now, from dhw_storage class, possible common class)
        type(pump)::pump_heat_pump	        							! heat pump's pump
        type(heat_exchanger)::HX_heat_pump              				! heat pump heat exchanger (evaporator)
		type(pipe)::pipe_heat_pump					    				! pipe from heat pump heat exchanger to network node B
        type(pump)::pump_cooling              							! cooling pump
        type(heat_exchanger)::HX_cooling                				! cooling heat exchanger 
		type(pipe)::pipe_cooling					    				! pipe from cooling heat exchanger to network node A
    !contains
        !procedure::load=>substation_002_load 							! load parameters in the specific type objects
        !procedure::create_networks=>substation_002_create_networks 	! create connections, pumps, heat exchangers, pipes 
        !procedure::num_elements=>substation_002_num_elements 			! return the number of elements (see num_elements type) 
        !procedure::initialize=>substation_002_initialize 				! get heating, cooling and dhw demands from the associated building units
        !procedure::solve=>substation_002_solve 						! calculate the energy flows in the time step
		!procedure::to_json=>substation_001_to_json						! generate json output
		!procedure::from_json=>substation_002_from_json					! read json output
    end type    

! =============================================================================================================
! energy_centre classes 
! =============================================================================================================

	! energy centre 000: groundwater reversible EHP+TES
    type, extends(energy_centre)::energy_centre_000
		type(ww_heat_pump_chiller)::WWHPC			 ! ww heat pump chiller
		type(pipe)::pipe_A							 ! pipe from buffer_tank to network's pipe A
		type(pipe)::pipe_B							 ! pipe from buffer_tank to network's pipe B
        type(buffer_tank)::TES  					 ! buffer_tank element
    contains
        procedure::load=>energy_centre_000_load ! load parameters in the specific type objects
        procedure::create_networks=>energy_centre_000_create_networks  ! create connections, pipes, tank
        procedure::num_elements=>energy_centre_000_num_elements ! return the number of elements (see num_elements type)
        procedure::initialize=>energy_centre_000_initialize ! set the initial parameters
        procedure::solve=>energy_centre_000_solve ! calculate the energy flows in the time step
        procedure::update=>energy_centre_000_update ! calculate the outputs
        procedure::to_json=>energy_centre_000_to_json
        procedure::from_json=>energy_centre_000_from_json
    end type    

    ! energy centre 001: reversible air-source EHP+TES
    type, extends(energy_centre)::energy_centre_001
		type(as_heat_pump_chiller)::ASHPC			 ! air-source heat pump chiller
		type(pipe)::pipe_A							 ! pipe from buffer_tank to network's pipe A
		type(pipe)::pipe_B							 ! pipe from buffer_tank to network's pipe B
        type(buffer_tank)::TES  					 ! buffer_tank element
    contains
        procedure::load=>energy_centre_001_load ! load parameters in the specific type objects
        procedure::create_networks=>energy_centre_001_create_networks  ! create connections, pipes, tank
        procedure::num_elements=>energy_centre_001_num_elements ! return the number of elements (see num_elements type)
        procedure::initialize=>energy_centre_001_initialize ! set the initial parameters
        procedure::solve=>energy_centre_001_solve ! calculate the energy flows in the time step
        procedure::update=>energy_centre_001_update ! calculate the outputs
        procedure::to_json=>energy_centre_001_to_json
        procedure::from_json=>energy_centre_001_from_json
    end type    

	! energy centre 002: borehole field
    type, extends(energy_centre)::energy_centre_002
		type(borehole_field)::bhe_field			! borehole field
    contains
        procedure::load=>energy_centre_002_load ! load parameters in the specific type objects
        procedure::create_networks=>energy_centre_002_create_networks  ! create connections and borehole_field
        procedure::num_elements=>energy_centre_002_num_elements ! return the number of elements (see num_elements type)
        procedure::initialize=>energy_centre_002_initialize ! set the initial parameters
        procedure::solve=>energy_centre_002_solve ! calculate the energy flows in the time step
        procedure::update=>energy_centre_002_update ! calculate the outputs
        procedure::to_json=>energy_centre_002_to_json
        procedure::from_json=>energy_centre_002_from_json
    end type    
       
 contains
 
    subroutine substation_wrapper_set(this,type_name)
		class(substation_wrapper)::this
		character(*)::type_name
		select case(type_name)
        case('S001')
            allocate(substation_001::this%item)
        !case('S00X')
        !    allocate(substation_00X::this%item)
        case default
			call simlog%output( 'Error: substation type not defined ('//type_name//').')
			stop			
		end select		
	end subroutine

    subroutine energy_centre_wrapper_set(this,type_name)
		class(energy_centre_wrapper)::this
		character(*)::type_name
		select case(type_name)
        case('C000')
            allocate(energy_centre_000::this%item)
        case('C001')
            allocate(energy_centre_001::this%item)
        case('C002')
            allocate(energy_centre_002::this%item)
        !case('C00X')
        !    allocate(energy_centre_00X::this%item)
        case default
			call simlog%output( 'Error: energy centre type not defined ('//type_name//').')
			stop			
		end select		
	end subroutine


! =============================================================================================================
! substation methods 
! =============================================================================================================

	 ! substation 001: HP(SH&DHW)+HX(SC)
	 subroutine substation_001_load(this)
		class(substation_001)::this
		integer::i
		do i=1,size(this%parameters)
            call this%parameters(i)%get_keys()
			select case(this%parameters(i)%key(1))
				case('HP')
					call this%HP%get_par(this%parameters(i))
				case('pump_heat_pump')
					this%pump_heat_pump%name=this%get_name('pump_heat_pump')
					call this%pump_heat_pump%get_par(this%parameters(i))
				case('HX_heat_pump')
					this%HX_heat_pump%name=this%get_name('HX_heat_pump')
					call this%HX_heat_pump%get_par(this%parameters(i))
				case('pipe_heat_pump')
					this%pipe_heat_pump%name=this%get_name('pipe_heat_pump')
					call this%pipe_heat_pump%get_par(this%parameters(i))
				case('pump_cooling')
					this%pump_cooling%name=this%get_name('pump_cooling')
					call this%pump_cooling%get_par(this%parameters(i))
				case('HX_cooling')
					this%HX_cooling%name=this%get_name('HX_cooling')
					call this%HX_cooling%get_par(this%parameters(i))
				case('pipe_cooling')
					this%pipe_cooling%name=this%get_name('pipe_cooling')
					call this%pipe_cooling%get_par(this%parameters(i))
                case ('T_grid_h_min')
					this%T_grid_h_min=this%parameters(i)%val_dbl    
                case ('T_grid_c_max')
					this%T_grid_c_max=this%parameters(i)%val_dbl    
			end select
		end do
	 end subroutine
	 
     ! create hydraulics nodes, pumps, heat exchangers for the substation_001
     subroutine substation_001_create_networks(this)
	    class(substation_001)::this
		integer::idx_node_A,idx_node_B 		! indexes of thermal grid nodes A and B
		integer::idx_node_HP1,idx_node_HP2 	! indexes of heat pump evaporator inlet and outlet nodes
		integer::idx_node_HX1,idx_node_HX2 	! indexes of cooling heat exchanger inlet and outlet nodes
        associate(thermal_grid=>this%p_thermal_grid)
		! create hydraulic nodes with unique coordinates based on substation position P0 
        idx_node_A=thermal_grid%add_node(name=this%get_name('node_A'),P0=this%P0)
        idx_node_B=thermal_grid%add_node(name=this%get_name('node_A'),P0=this%P0%translate((/0d0,0d0,1d0/)))
        idx_node_HP1=thermal_grid%add_node(name=this%get_name('node_HP1'),P0=this%P0%translate((/-1d0,1d0,0d0/)))
        idx_node_HP2=thermal_grid%add_node(name=this%get_name('node_HP2'),P0=this%P0%translate((/-1d0,1d0,1d0/)))
        idx_node_HX1=thermal_grid%add_node(name=this%get_name('node_HX1'),P0=this%P0%translate((/+1d0,1d0,1d0/)))
        idx_node_HX2=thermal_grid%add_node(name=this%get_name('node_HX2'),P0=this%P0%translate((/+1d0,1d0,0d0/)))
		! connect hydraulic elements to the hydraulic nodes	
		call this%pump_heat_pump%connect(node_in=idx_node_A, node_out=idx_node_HP1)
		call this%HX_heat_pump%connect(node_in=idx_node_HP1, node_out=idx_node_HP2)
		call this%pipe_heat_pump%connect(node_in=idx_node_HP2, node_out=idx_node_B)
		call this%pump_cooling%connect(node_in=idx_node_B, node_out=idx_node_HX1)
		call this%HX_cooling%connect(node_in=idx_node_HX1, node_out=idx_node_HX2)
		call this%pipe_cooling%connect(node_in=idx_node_HX2, node_out=idx_node_A)	
		! add elements to the grid
        call thermal_grid%add_element(this%pump_heat_pump)
        call thermal_grid%add_element(this%HX_heat_pump)
        call thermal_grid%add_element(this%pipe_heat_pump)
        call thermal_grid%add_element(this%pump_cooling)
        call thermal_grid%add_element(this%HX_cooling)
        call thermal_grid%add_element(this%pipe_cooling)		
	 end associate
     end subroutine
	 
	 function substation_001_num_elements(this) result(num)
		class(substation_001)::this
		type(num_elements)::num
		num%num_pipes=2
		num%num_heat_exchangers=2
		num%num_pumps=2
		call num%set_num_thermal_nodes()
	 end function
	 
     subroutine substation_001_initialize(this)
		class(substation_001),intent(inout)::this
		this%HP%component_feature=ALGEBRAIC
        this%HX_heat_pump%T_aux_h=this%T_grid_h_min
        this%HX_cooling%T_aux_c=this%T_grid_c_max
	 end subroutine
	 
     subroutine substation_001_solve(this,hour)
		class(substation_001)::this
		integer::hour,previous_hour
		double precision::DT_HX_heat_pump,DT_HX_cooling,Q_HP_condenser,Q_HX_heat_pump,Q_HX_cooling
        ! calculate auxiliary heating and cooling that, if different than zero, reval a capacity adequacy problem of the network
        this%Q_aux_h(hour)=this%HX_heat_pump%get_aux_h(hour)
        this%Q_aux_c(hour)=this%HX_cooling%get_aux_c(hour)
		! calculate energy flows, DHW storage T (assumption: one DHW tank for each associated building unit), mass flow at network heat exchangers based on fixed DT
		DT_HX_heat_pump=-5d0 ! set the delta T of the network stream across the heat pump heat exchanger (water is cooled, negative value)
		DT_HX_cooling=5d0 ! set the delta T of the network stream across the cooling heat exchanger (water is heated, positive value)
		call this%DHW_tank%calculate(hour) ! update state and outputs in DHW tank
		! overall heat requested to the heat pump for dhw storage charging and space heating, Wh
		this%HP%Q_cond(hour)=this%DHW_tank%Q_supply(hour)+this%Q_supply_heat(hour)
		this%HP%T_hwo=50d0 ! to do ... set this values based on heat distribution temperature and DHW tank temperature
		this%HP%T_cwo=this%HX_heat_pump%thermal_nodes(1)%T(hour)
		call this%HP%calculate(hour) ! NOTE: Q_cond can be downgraded internally if HP capacity is not sufficient
		! heat transfers to the grid
		Q_HX_heat_pump=min(0d0,-this%HP%Q_evap(hour)) ! heat positive when supplied to the network, thus we must ensure a negative value
		Q_HX_cooling=max(0d0,-this%Q_supply_cool(hour)) ! heat positive when supplied to the network, thus we must ensure a positive value
		! save Q in the thermal nodes related to the elements HX heat pump and the HX cooling 
		this%HX_heat_pump%thermal_nodes(1)%Qin(hour)=Q_HX_heat_pump
		this%HX_cooling%thermal_nodes(1)%Qin(hour)=Q_HX_cooling
		! impose the mass flow rates of the two pumps
		this%pump_heat_pump%fix_flow(hour)=1
		this%pump_heat_pump%mdot(hour)=Q_HX_heat_pump/(this%HX_heat_pump%fluid%cp*DT_HX_heat_pump)
		this%pump_cooling%fix_flow(hour)=1
		this%pump_cooling%mdot(hour)=Q_HX_cooling/(this%HX_cooling%fluid%cp*DT_HX_cooling)
	 end subroutine	
	 
	 subroutine substation_001_to_json(this,json,obj)
		class(substation_001)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(out)::obj
		call substation_to_json(this,json,obj)
		call this%HP%to_json(json,obj,'HP')
		call this%pump_heat_pump%to_json(json,obj,'pump_heat_pump')
		call this%HX_heat_pump%to_json(json,obj,'HX_heat_pump')
		call this%pipe_heat_pump%to_json(json,obj,'pipe_heat_pump')
		call this%pump_cooling%to_json(json,obj,'pump_cooling')
		call this%HX_cooling%to_json(json,obj,'HX_cooling')
		call this%pipe_cooling%to_json(json,obj,'pipe_cooling')
   		call save_series_bin(json,obj,'Q_aux_h',this%Q_aux_h)
		call save_series_bin(json,obj,'Q_aux_c',this%Q_aux_c)
	 end subroutine
	 
	 subroutine substation_001_from_json(this,json,obj,path)
		class(substation_001)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(in)::obj
		character(*)::path
		call substation_from_json(this,json,obj,path)
		call this%HP%from_json(json,obj,path//'.HP')
		call this%pump_heat_pump%from_json(json,obj,path//'.pump_heat_pump')
		call this%HX_heat_pump%from_json(json,obj,path//'.HX_heat_pump')
		call this%pipe_heat_pump%from_json(json,obj,path//'.pipe_heat_pump')
		call this%pump_cooling%from_json(json,obj,path//'.pump_cooling')
		call this%HX_cooling%from_json(json,obj,path//'.HX_cooling')
		call this%pipe_cooling%from_json(json,obj,path//'.pipe_cooling')
   		call read_series_bin(json,obj,path//'.Q_aux_h',this%Q_aux_h)
   		call read_series_bin(json,obj,path//'.Q_aux_c',this%Q_aux_c)
     end subroutine
	 
	 ! substation 001: end

! =============================================================================================================
! energy_centre methods 
! =============================================================================================================

	 ! energy centre 000: groundwater reversible heat pump + TES
	 subroutine energy_centre_000_load(this)
		class(energy_centre_000)::this
		integer::i
		do i=1,size(this%parameters)
            call this%parameters(i)%get_keys()
			select case(this%parameters(i)%key(1))
				case('WWHPC')
					call this%WWHPC%get_par(this%parameters(i))                
				case('TES')
					this%TES%name=this%get_name('TES')
                    call this%TES%get_par(this%parameters(i))
				case('pipe_A')
					this%pipe_A%name=this%get_name('pipe_A')
                    call this%pipe_A%get_par(this%parameters(i))
				case('pipe_B')
					this%pipe_B%name=this%get_name('pipe_B')
                    call this%pipe_B%get_par(this%parameters(i))
			end select
		end do
	 end subroutine
	 
     ! create hydraulic nodes (storage) and pipes for the energy_centre_001
     subroutine energy_centre_000_create_networks(this)
	    class(energy_centre_000)::this
        integer:: idx_node_TES 					! index for the central node of the TES
		integer:: idx_node_A, idx_node_B 		! indexes for connections to thermal grid pipes A and B
		double precision::D_pipe 				! connection pipe diameter, mm
		double precision::M,UA,MC_wall 			! tank mass of water (M), UA coefficient (W/K) and thermal capacity of the wall (J/K)
        associate(thermal_grid=>this%p_thermal_grid)
		! create hydraulic nodes with unique coordinates based on energy_centre position P0 
        idx_node_A=thermal_grid%add_node(name=this%get_name('node_A'),P0=this%P0)
        idx_node_B=thermal_grid%add_node(name=this%get_name('node_B'),P0=this%P0%translate((/0d0,0d0,1d0/)))
        idx_node_TES=thermal_grid%add_node(name=this%get_name('node_TES'),P0=this%P0%translate((/-1d0,0d0,0.5d0/)))
		! connect hydraulic elements to the hydraulic nodes	
		call this%pipe_A%connect(node_in=idx_node_A,node_out=idx_node_TES)
		call this%pipe_B%connect(node_in=idx_node_TES,node_out=idx_node_B)
		call this%TES%connect(hydraulic_node=idx_node_TES, element_to_top=this%pipe_A, element_to_base=this%pipe_B)
		! add hydraulic elements to the thermal grid
		call thermal_grid%add_element(this%pipe_A)
		call thermal_grid%add_element(this%pipe_B)
		call thermal_grid%add_element(this%TES)
		end associate		
     end subroutine
     
	 function energy_centre_000_num_elements(this) result(num)
		class(energy_centre_000)::this
		type(num_elements)::num
		num%num_pipes=2
		num%num_tanks=1
		call num%set_num_thermal_nodes()
	 end function
	 
     subroutine energy_centre_000_initialize(this)
		class(energy_centre_000)::this
		! to do ... use this method to perform automatic sizing
		this%WWHPC%component_feature=DYNAMICAL
	 end subroutine		 
	 
     subroutine energy_centre_000_solve(this,hour,time)
		class(energy_centre_000)::this
		integer::hour
		double precision,dimension(2)::time
		! invoke the buffer controller to test its status and update its mode (signal)
		call this%TES%controller%update(hour,this%TES%T_middle(hour))
		! test the signal of the controller 
		select case(this%TES%controller%mode)
		case (COOLING_ON)
			! WWHPC settings for cooling mode
			this%WWHPC%operation_mode=COOLING
			this%WWHPC%T_wo_set=this%TES%controller%T_storage_max(hour)-this%WWHPC%DT_w  ! chilled water set point, DT_w below the maximum storage temperature (fast charge)
			this%WWHPC%T_wi=this%TES%T_top(hour)					 ! heat pump inlet water is taken from the TES top in cooling mode			
			this%WWHPC%T_source_sink=max(T_ground(hour)+10d0,25d0)         ! minimum operating temperature of the water leaving the condenser			
			call this%WWHPC%calculate(hour)
			call this%TES%charge_base(hour,this%WWHPC%T_wo,this%WWHPC%mdot_w)
		case (HEATING_ON)
			! WWHPC settings for heating mode
			this%WWHPC%operation_mode=HEATING
			this%WWHPC%T_wo_set=this%TES%controller%T_storage_min(hour)+this%WWHPC%DT_w  ! hot water set point, DT_w above the minimum storage temperature (fast charge) 
			this%WWHPC%T_wi=this%TES%T_base(hour)					 ! heat pump inlet water is taken from the TES base in heating mode	            
			this%WWHPC%T_source_sink=T_ground(hour)-7d0                    ! operating temperature of the water leaving the evaporator			
			call this%WWHPC%calculate(hour)
			call this%TES%charge_top(hour,this%WWHPC%T_wo,this%WWHPC%mdot_w)
		case (INACTIVE)
			! WWHPC in standby mode 
			this%WWHPC%operation_mode=STANDBY				
			call this%WWHPC%calculate(hour)
			call this%TES%charge_zero(hour)
		end select	
	 end subroutine		 

     subroutine energy_centre_000_update(this,hour,time) ! update output data
		class(energy_centre_000)::this
		integer::hour
		double precision,dimension(2)::time
		call this%WWHPC%update(hour,time)
	 end subroutine	 
	 	 
 	 subroutine energy_centre_000_to_json(this,json,obj)
		class(energy_centre_000)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(out)::obj
		call energy_centre_to_json(this,json,obj)
		call this%WWHPC%to_json(json,obj,'WWHPC')
		call this%pipe_A%to_json(json,obj,'pipe_A')
		call this%pipe_B%to_json(json,obj,'pipe_B')
		call this%TES%to_json(json,obj,'TES')
	 end subroutine
	 
 	 subroutine energy_centre_000_from_json(this,json,obj,path)
		class(energy_centre_000)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(in)::obj
		character(*)::path
		call energy_centre_from_json(this,json,obj,path)
		call this%WWHPC%from_json(json,obj,path//'.WWHPC')
		call this%pipe_A%from_json(json,obj,path//'.pipe_A')
		call this%pipe_B%from_json(json,obj,path//'.pipe_B')
		call this%TES%from_json(json,obj,path//'.TES')
	 end subroutine	 
	 
	 ! energy centre 000: end
     
     
	 ! energy centre 001: Air-source reversible heat pump + TES
	 subroutine energy_centre_001_load(this)
		class(energy_centre_001)::this
		integer::i
		do i=1,size(this%parameters)
            call this%parameters(i)%get_keys()
			select case(this%parameters(i)%key(1))
				case('ASHPC')
					call this%ASHPC%get_par(this%parameters(i))                
				case('TES')
					this%TES%name=this%get_name('TES')
                    call this%TES%get_par(this%parameters(i))
				case('pipe_A')
					this%pipe_A%name=this%get_name('pipe_A')
                    call this%pipe_A%get_par(this%parameters(i))
				case('pipe_B')
					this%pipe_B%name=this%get_name('pipe_B')
                    call this%pipe_B%get_par(this%parameters(i))
			end select
		end do
	 end subroutine
	 
     ! create hydraulic nodes (storage) and pipes for the energy_centre_001
     subroutine energy_centre_001_create_networks(this)
	    class(energy_centre_001)::this
        integer:: idx_node_TES 					! index for the central node of the TES
		integer:: idx_node_A, idx_node_B 		! indexes for connections to thermal grid pipes A and B
		double precision::D_pipe 				! connection pipe diameter, mm
		double precision::M,UA,MC_wall 			! tank mass of water (M), UA coefficient (W/K) and thermal capacity of the wall (J/K)
        associate(thermal_grid=>this%p_thermal_grid)
		! create hydraulic nodes with unique coordinates based on energy_centre position P0 
        idx_node_A=thermal_grid%add_node(name=this%get_name('node_A'),P0=this%P0)
        idx_node_B=thermal_grid%add_node(name=this%get_name('node_B'),P0=this%P0%translate((/0d0,0d0,1d0/)))
        idx_node_TES=thermal_grid%add_node(name=this%get_name('node_TES'),P0=this%P0%translate((/-1d0,0d0,0.5d0/)))
		! connect hydraulic elements to the hydraulic nodes	
		call this%pipe_A%connect(node_in=idx_node_A,node_out=idx_node_TES)
		call this%pipe_B%connect(node_in=idx_node_TES,node_out=idx_node_B)
		call this%TES%connect(hydraulic_node=idx_node_TES, element_to_top=this%pipe_A, element_to_base=this%pipe_B)
		! add hydraulic elements to the thermal grid
		call thermal_grid%add_element(this%pipe_A)
		call thermal_grid%add_element(this%pipe_B)
		call thermal_grid%add_element(this%TES)
		end associate		
     end subroutine
     
	 function energy_centre_001_num_elements(this) result(num)
		class(energy_centre_001)::this
		type(num_elements)::num
		num%num_pipes=2
		num%num_tanks=1
		call num%set_num_thermal_nodes()
	 end function
	 
     subroutine energy_centre_001_initialize(this)
		class(energy_centre_001)::this
		! to do ... use this method to perform automatic sizing
		this%ASHPC%component_feature=DYNAMICAL
     end subroutine	
     	 
     subroutine energy_centre_001_solve(this,hour,time)
		class(energy_centre_001)::this
		integer::hour
		double precision,dimension(2)::time
		! invoke the buffer controller to test its status and update its mode (signal)
		call this%TES%controller%update(hour,this%TES%T_middle(hour))
		! test the signal of the controller 
		select case(this%TES%controller%mode)
		case (COOLING_ON)
			! ASHPC settings for cooling mode
			this%ASHPC%operation_mode=COOLING
			this%ASHPC%T_wo_set=this%TES%controller%T_storage_max(hour)-this%ASHPC%DT_w  ! chilled water set point, DT_w below the maximum storage temperature (fast charge)
			this%ASHPC%T_wi=this%TES%T_top(hour)					 ! heat pump inlet water is taken from the TES top in cooling mode			
			this%ASHPC%T_air=this%p_building_complex%meteodata%T_amb(hour)
			call this%ASHPC%calculate(hour)
			call this%TES%charge_base(hour,this%ASHPC%T_wo,this%ASHPC%mdot_w)
		case (HEATING_ON)
			! ASHPC settings for heating mode
			this%ASHPC%operation_mode=HEATING
			this%ASHPC%T_wo_set=this%TES%controller%T_storage_min(hour)+this%ASHPC%DT_w  ! hot water set point, DT_w above the minimum storage temperature (fast charge)
			this%ASHPC%T_wi=this%TES%T_base(hour)					 ! heat pump inlet water is taken from the TES base in heating mode	            
			this%ASHPC%T_air=this%p_building_complex%meteodata%T_amb(hour)
			call this%ASHPC%calculate(hour)
			call this%TES%charge_top(hour,this%ASHPC%T_wo,this%ASHPC%mdot_w)
		case (INACTIVE)
			! ASHPC in standby mode 
			this%ASHPC%operation_mode=STANDBY				
			call this%ASHPC%calculate(hour)
			call this%TES%charge_zero(hour)
		end select	
	 end subroutine		 

     subroutine energy_centre_001_update(this,hour,time) ! update output data
		class(energy_centre_001)::this
		integer::hour
		double precision,dimension(2)::time
		call this%ASHPC%update(hour,time)
	 end subroutine	 
	 	 
 	 subroutine energy_centre_001_to_json(this,json,obj)
		class(energy_centre_001)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(out)::obj
		call energy_centre_to_json(this,json,obj)
		call this%ASHPC%to_json(json,obj,'ASHPC')
		call this%pipe_A%to_json(json,obj,'pipe_A')
		call this%pipe_B%to_json(json,obj,'pipe_B')
		call this%TES%to_json(json,obj,'TES')
	 end subroutine
	 
 	 subroutine energy_centre_001_from_json(this,json,obj,path)
		class(energy_centre_001)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(in)::obj
		character(*)::path
		call energy_centre_from_json(this,json,obj,path)
		call this%ASHPC%from_json(json,obj,path//'.ASHPC')
		call this%pipe_A%from_json(json,obj,path//'.pipe_A')
		call this%pipe_B%from_json(json,obj,path//'.pipe_B')
		call this%TES%from_json(json,obj,path//'.TES')
	 end subroutine	 
	 
	 ! energy centre 001: end

	 ! energy centre 002: Borehole field
 	 subroutine energy_centre_002_load(this)
		class(energy_centre_002)::this
		integer::i
		do i=1,size(this%parameters)
            call this%parameters(i)%get_keys()
			select case(this%parameters(i)%key(1))
				case('bhe_field')
				this%bhe_field%name=this%get_name('bhe_field')		
                call this%bhe_field%get_par(this%parameters(i))
			end select
		end do
		this%bhe_field%bhe%g_function=g_function()
     end subroutine	
	 
     ! create hydraulic nodes and elements for the energy_centre_002
     subroutine energy_centre_002_create_networks(this)
	    class(energy_centre_002)::this
		integer:: idx_node_A, idx_node_B ! indexes for connections to thermal grid pipes A and B
        associate(thermal_grid=>this%p_thermal_grid)
		! create hydraulic nodes with unique coordinates based on energy_centre position P0 
        idx_node_A=thermal_grid%add_node(name=this%get_name('node_A'),P0=this%P0)
        idx_node_B=thermal_grid%add_node(name=this%get_name('node_B'),P0=this%P0%translate((/0d0,0d0,1d0/)))
		! connect hydraulic elements to the hydraulic nodes	        
		call this%bhe_field%connect(node_in=idx_node_A,node_out=idx_node_B)
		! add element to the thermal grid
		call thermal_grid%add_element(this%bhe_field)
		end associate		
     end subroutine
	 
     function energy_centre_002_num_elements(this) result(num)
		class(energy_centre_002)::this
		type(num_elements)::num
		num%num_heat_exchangers=1
		call num%set_num_thermal_nodes()
	 end function
	 
     subroutine energy_centre_002_initialize(this)
		class(energy_centre_002)::this
        this%bhe_field%bhe%Tg=this%p_building_complex%meteodata%get_T_ground()
		call this%bhe_field%bhe%initialize()
	 end subroutine		 
     
	 subroutine energy_centre_002_solve(this,hour,time)
		class(energy_centre_002)::this
		integer::hour
		double precision,dimension(2)::time
		associate(bhe=>this%bhe_field%bhe, bhe_field=>this%bhe_field)
		! calculate Qin as the Qex extracted from the borehole field
		bhe_field%thermal_nodes(1)%Qin(hour)=bhe%Qex(&
				hour,&								! the current hour
				bhe_field%mdot(hour),&				! fluid mass flow rate
				bhe_field%inflow_T(hour),&      	! fluid inlet temperature (from adjacent element)
				bhe_field%thermal_nodes(1)%T(hour)) ! fluid temperature (in the control volume)
		end associate
	 end subroutine		 

     subroutine energy_centre_002_update(this,hour,time) ! update output data
		class(energy_centre_002)::this
		integer::hour,i
		double precision,dimension(2)::time
		call this%bhe_field%bhe%update(hour,time)
	 end subroutine

 	 subroutine energy_centre_002_to_json(this,json,obj)
		class(energy_centre_002)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(out)::obj
		type(json_value),pointer::bhe_field_obj
		call energy_centre_to_json(this,json,obj)
		call this%bhe_field%to_json(json,obj,'bhe_field')
	 end subroutine

 	 subroutine energy_centre_002_from_json(this,json,obj,path)
		class(energy_centre_002)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(in)::obj
		character(*)::path
		call energy_centre_from_json(this,json,obj,path)
		call this%bhe_field%from_json(json,obj,path//'.bhe_field')
	 end subroutine

	 ! energy centre 002: end


end module