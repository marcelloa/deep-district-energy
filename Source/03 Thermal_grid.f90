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
module thermal_grid_module
	use constants_module
	use building_complex_module
	use thermal_plants_module
    use datalog
    use meteolib
    use maths
    use utility
	use json_module
    implicit none

	double precision,parameter::TOL_T=0.1d0 ! K, maximum temperature variation during integration
    double precision,parameter::MIN_FLOW_RATE=1d-3 ! kg/s, minimum mass flow rate 
    double precision,parameter::MIN_MC=TIMESTEP*MIN_FLOW_RATE*CP_WATER ! J/K, minimum thermal capacity
	! thermal solver modes
    enum, bind(c)
        enumerator::CAPACITIVE_ADVECTIVE=1,NON_CAPACITIVE_ADVECTIVE,CAPACITIVE_NON_ADVECTIVE,ADVECTIVE
    end enum 
    
! =============================================================================================================
! core classes 
! =============================================================================================================
    
	! a path of segments interconnecting an ordered set of points	
    type polyline
        type(point),dimension(:),allocatable::points	! ordered set of points
        double precision::length						! path length, a calculated property (m)
        type(point)::P1									! the start point, a calculate property (x,y,z in m)
        type(point)::P2									! the end point, a calculate property (x,y,z in m)
    contains
        procedure::set_properties=>polyline_set_properties ! calculate properties starting from the set of points
    end type    
        	
	double precision,dimension(8760),target::T_ground=12d0		! ground temperature (°C)
	double precision,dimension(8760),target::T_unheated=18d0	! temperature of unheated spaces, e.g. technical room (°C)
	
	! record containing the number of hydraulic elements and thermal nodes
	type num_elements
		integer::num_pipes=0
		integer::num_heat_exchangers=0
		integer::num_pumps=0
		integer::num_tanks=0
		integer::num_thermal_nodes=0
	contains
		procedure::hydraulic_channels_like=>num_elements_hydraulic_channels_like
		procedure::non_hydraulic_channels_like=>num_elements_non_hydraulic_channels_like
		procedure::set_num_thermal_nodes=>num_elements_set_num_thermal_nodes
	end type	
			
	! operator that sums two num_elements objects
    interface operator (+)
        procedure add_num_elements
    end interface operator  (+)	
	
	! the progenitor of all hydraulic elements, both the channel-like and the non-channel-like ones	
	type hydraulic_element
		character(:),allocatable::name
        type(thermal_node),dimension(:),pointer::thermal_nodes=>null()	! pointer to the thermal node subarray
		type(incompressible_fluid),pointer::fluid=>null()				! pointer to the thermal grid fluid
        double precision::UA=0d0										! overall heat loss coefficient (W/K)
        double precision::MC=MIN_MC										! overall thermal capacity (J/K)
	contains
		procedure::to_json=>hydraulic_element_to_json
		procedure::from_json=>hydraulic_element_from_json
    end type

	! the progenitor of all channel-like hydraulic elements
	type,extends(hydraulic_element)::hydraulic_channel
        integer::node_in=0												! index of inlet hydraulic node
        integer::node_out=0												! index of outlet hydraulic node
        integer,dimension(8760)::fix_flow=0								! 1 if flow is fixed otherwise 0
		double precision,dimension(8760)::mdot=0d0						! mass flow rate, kg/s, positive from node_in to node_out
		double precision::a_hyd=0d0,E_hyd=0d0,dEdm_hyd=0d0,Em_hyd=0d0	! work variables (hydraulic solver) => to remove ...
	contains
		procedure::connect=>hydraulic_channel_connect					! assign node_in and node_out
		procedure::initialize=>hydraulic_channel_initialize              ! initialization, e.g. MC calculation
		procedure::to_json=>hydraulic_channel_to_json
		procedure::from_json=>hydraulic_channel_from_json
    end type
    	
	! wrapper to create arrays of polymorphic hydraulic elements 	
	type hydraulic_element_wrapper
		class(hydraulic_element),pointer::item
	end type

	! wrapper to create arrays of polymorphic channel-like hydraulic elements
	type hydraulic_channel_wrapper
		class(hydraulic_channel),pointer::item
	end type
	
	! pipe dimension and hydraulic and thermal characteristic
    type,extends(hydraulic_channel)::pipe
        double precision::D                         ! internal diameter (m)
        double precision::L                         ! length (m)
        double precision::K                         ! minor pressure drops coefficient (-)
        double precision::eps                   	! relative roughness, pipe roughness/pipe diameter (-)
		double precision::ins_delta					! insulation thickness (m)
		double precision::ins_lambda				! insulation thermal conductivity (W/(m,K))
	contains	
		procedure::set=>pipe_set 					! assign node_in and node_out
		procedure::initialize=>pipe_initialize       ! calculate MC and UA
		procedure::to_json=>pipe_to_json
		procedure::from_json=>pipe_from_json
		procedure::get_par=>pipe_get_par
    end type
	interface pipe
        module procedure init_pipe					! constructor
    end interface
	
	type inflow
		double precision,dimension(:),pointer::T 	! the array of temperatures 1..8760
		double precision::mdot 						! the mass flow rate last computed
	end type

	! heat exchanger stream and its hydraulic characteristic
    type,extends(hydraulic_channel)::heat_exchanger                             
        double precision::DP_nom                    	! nominal pressure drop (Pa)
        double precision::m_nom                     	! nominal mass flow rate (kg/s)
        double precision::Q_nom                     	! nominal heat duty (W)
        type(inflow),dimension(:),allocatable::inflows 	! array of inflow types 
		integer::num_inflows=0							! number of inflows
        double precision::T_aux_h=-1d12                 ! limiting temperature below which a virtual auxiliary heater is activated 
        double precision::T_aux_c=+1d12					! limiting temperature above which a virtual auxiliary cooler is activated
	contains
		procedure::inflow_T=>heat_exchanger_inflow_T 	! get the temperature of the adjacent thermal node that supplies a positive mass flow 
		procedure::set=>heat_exchanger_set 				! assign node_in and node_out
		procedure::initialize=>heat_exchanger_initialize
		procedure::to_json=>heat_exchanger_to_json
		procedure::from_json=>heat_exchanger_from_json
   		procedure::get_par=>heat_exchanger_get_par
        procedure::get_aux_h=>heat_exchanger_get_aux_h
        procedure::get_aux_c=>heat_exchanger_get_aux_c
    end type
       
	! pump and its hydraulic characteristic
    type,extends(hydraulic_channel)::pump                             
        double precision::a1                        ! coefficient a1 in H=a1+a2/m+a3/m^2, (m.w.c.)
        double precision::a2                        ! coefficient a2 in H=a1+a2/m+a3/m^2, (m.w.c.*(kg/s))
        double precision::a3                        ! coefficient a3 in H=a1+a2/m+a3/m^2, (m.w.c.*(kg/s)^2)
		double precision,dimension(8760)::W=0d0		! pump electricity consumption, W
	contains	
		procedure::set=>pump_set 					! assign node_in and node_out
		procedure::to_json=>pump_to_json
		procedure::from_json=>pump_from_json
   		procedure::get_par=>pump_get_par
    end type

	! borehole field water stream and its hydraulic and thermal characteristic
    type,extends(heat_exchanger)::borehole_field                             
		type(borehole_heat_exchanger)::bhe          ! the borehole heat exchanger object with thermal data and methods
	contains	
		procedure::set=>borehole_field_set			! assign node_in and node_out
		procedure::initialize=>borehole_field_initialize
		procedure::to_json=>borehole_field_to_json
		procedure::from_json=>borehole_field_from_json
		procedure::get_par=>borehole_field_get_par
    end type
	interface borehole_field
        module procedure init_borehole_field		! constructor
    end interface
	
	integer,parameter::TANK_NODES=3					! number of thermal nodes (i.e., horizontal isothermal layers) in a tank 
	
	! tank with hydraulic and thermal characteristics, and topological data (refencing index and pointers to hydraulic node and inlet / outlet channels) 
    type,extends(hydraulic_element)::tank 
		! hydraulic modelling
        double precision::K_in=1d0                   ! minor pressure drops coefficient at tank inlet (-)
        double precision::K_out=1d0                  ! minor pressure drops coefficient at tank outlet (-)
		! thermal modelling 
        integer::num_thermal_nodes=TANK_NODES       ! number of thermal nodes (i.e., horizontal isothermal layers) 
		! topological information
		integer::hydraulic_node=0					! the hydraulic node representing the tank central point in the hydraulic network
		class(hydraulic_channel),pointer::element_to_top=>null()	! the element (e.g., pipe) connected to the top layer of the tank
		class(hydraulic_channel),pointer::element_to_base=>null()	! the element (e.g., pipe) connected to the base layer of the tank
	contains
		procedure::connect=>tank_connect ! assign hydraulic node and pointers to elements _to_top and _to_base
		procedure::mass_flow_in_top=>tank_mass_flow_in_top ! get the mass flow rate, positive if entering in the top layer (and leaving from the base layer) 
		procedure::Q_advection=>tank_Q_advection ! advection heat transfer, W  
		procedure::Q_loss=>tank_Q_loss ! heat losses to the environment, W  
		procedure::T_middle=>tank_T_middle ! return the temperature in the middle of the tank
		procedure::T_top=>tank_T_top ! return the temperature at the top of the tank
		procedure::T_base=>tank_T_base ! return the temperature at the base of the tank
		procedure::charge_top=>tank_charge_top ! charge the tank from the top
		procedure::charge_base=>tank_charge_base ! charge the tank from the base
		procedure::charge_zero=>tank_charge_zero ! set to zero the charging heating/cooling power
		procedure::to_json=>tank_to_json
		procedure::from_json=>tank_from_json
    end type
			
	! buffer tank with inlet/outlet pressure drops and control signals, charged by flow rate at given T from top or base 
    type,extends(tank)::buffer_tank 
		type(buffer_controller)::controller 
	contains
		procedure::set=>buffer_tank_set
		procedure::to_json=>buffer_tank_to_json
		procedure::from_json=>buffer_tank_from_json
   		procedure::get_par=>buffer_tank_get_par        
    end type
    
	! connection (hydraulic node) with its hydraulic and topological information
    type connection 								
		character(:),allocatable::name				! name
        double precision,dimension(8760)::H         ! hydraulic head (m.w.c)
        integer,dimension(8760)::fix_height         ! 1 if head is fixed (e.g., expansion vessel) otherwise 0   
		type(point)::P0 							! space coordinates(m) of connection in the local reference system
		! topological information of the hydraulic network
		integer,dimension(:),allocatable::channels  ! index of the channel-like elements connected to the node
		double precision,dimension(:),allocatable::signs  ! flow convention: 1 = positive if entering the node, -1 = positive if leaving the node
		double precision,dimension(:),allocatable::a,b ! coefficients for calculating positive mass flow rates
	contains
		procedure::to_json=>connection_to_json
		procedure::from_json=>connection_from_json
    end type
        
  	! the hydraulic network as group of connections (hydraulic nodes) interconnected by channel-like hydraulic elements
    type hydraulic_network
        type(connection),dimension(:),allocatable::node ! the array of hydraulic nodes (i.e. known hydraulic heads)
		type(hydraulic_channel_wrapper),dimension(:),allocatable::elements ! the array of channel-like elements, i.e., pipes, heat exchangers, and pumps
        ! counters for each element type
        integer::num_pipes=0
		integer::num_heat_exchangers=0
		integer::num_pumps=0
		! network output
        double precision,dimension(:,:),allocatable::m_new_tot
        double precision,dimension(:,:),allocatable::H       
       ! work matrices
        double precision,dimension(:,:),allocatable::Aa_hyd
        double precision,dimension(:,:),allocatable::H_fix
        double precision,dimension(:,:),allocatable::q_hyd
        double precision,dimension(:,:),allocatable::A10
        double precision,dimension(:,:),allocatable::A21
        double precision,dimension(:,:),allocatable::A10_tot
        double precision,dimension(:,:),allocatable::A21_tot   
        double precision,dimension(:,:),allocatable::E_hyd
        double precision,dimension(:,:),allocatable::E_hyd_tot
        double precision,dimension(:,:),allocatable::dEdm_hyd
        double precision,dimension(:,:),allocatable::dEdm_hyd_tot
        double precision,dimension(:,:),allocatable::Em_hyd
        double precision,dimension(:,:),allocatable::Em_hyd_tot
        double precision,dimension(:,:),allocatable::m
        double precision,dimension(:,:),allocatable::m_new
        double precision,dimension(:,:),allocatable::a_hyd
        double precision,dimension(:,:),allocatable::a_hyd_tot
        double precision,dimension(:,:),allocatable::b_hyd
        double precision,dimension(:,:),allocatable::A12
        double precision,dimension(:,:),allocatable::Jacob
        double precision,dimension(:,:),allocatable::Jacob_inv
        double precision,dimension(:,:),allocatable::A11
        double precision,dimension(:,:),allocatable::A_hyd_sys
    end type
           
	! thermal node, concentrated thermal capacity with UA for heat loss to the environment and heat input
	type thermal_node
		integer::idx=0 ! the index in the array of thermal nodes within a thermal network
		double precision::MC  ! thermal capacity, J/K
		double precision::UA  ! overall heat loss coefficient, W/K
		double precision,dimension(8760)::Qin=0d0 ! heat supplied in input to the thermal node
		double precision,dimension(8760)::T=12d0  ! temperature of the thermal node 
		double precision,dimension(:),pointer::T_env ! the temperature of the environment (ground, unheated space)
	contains
		procedure::to_json=>thermal_node_to_json
		procedure::from_json=>thermal_node_from_json
	end type	
  
	! positive defined mass flow rate for mass exchange between two thermal nodes
	! NOTE the mass flow is a virtual pipe between two thermal nodes, node_in and node_out, positive from node_in to node_out
	type mass_flow
		integer::node_in							! index of thermal node that is the inlet for the mass flow virtual pipe (so the flow leaves from node_in thermal node) 
		integer::node_out							! index of thermal node that is the outlet for the mass flow virtual pipe (so the flow enters in node_out thermal node) 
		double precision::mass_flow_rate			! mass flow rate, kg/s
	end type
	
	! the thermal network as group of thermal nodes interconnected by mass flows
	type thermal_network
		type(thermal_node),dimension(:),allocatable::thermal_nodes ! array of thermal nodes
		type(mass_flow),dimension(:),allocatable::mass_flows    ! array of positive defined mass flows, each connecting two thermal nodes
		integer::num_thermal_nodes ! counter of thermal nodes
		! work matrices
		double precision,dimension(:,:),allocatable::A,B,T
		double precision,dimension(:),allocatable::C,D,mdotcp
		integer,dimension(:),allocatable::solver_mode
		integer,dimension(:),allocatable::algebraic_idx
	end type	

    ! connection between two or more double pipes
    type double_connection
        type(thermal_grid),pointer::p_thermal_grid
        type(point)::P0							! coordinates must coincide with start or end point coordinates of double pipe
        character(:),allocatable::name			! name of the connection
        integer::idx_node_A						! index of node_A (hot, upper level) in hydraulic network's nodes array
        integer::idx_node_B						! index of node_B (cold, lower level) in hydraulic network's nodes array
    contains
        procedure::create_networks=>double_connection_create_networks ! create hydraulic connections
		procedure::get_name=>double_connection_get_name
		procedure::to_json=>double_connection_to_json
		procedure::from_json=>double_connection_from_json
    end type  
    
    ! double pipe (forward and return)
    type double_pipe
        type(thermal_grid),pointer::p_thermal_grid
        ! inputs
        character(:),allocatable::name 		! name of the pipe
        type(polyline)::path 				! sequence of points (pipe's path can be a broken line)
        double precision::D_int  		 	! internal diameter, mm
        double precision::ins_lambda     	! insulation thermal conductivity, W/(m,K)
        double precision::ins_delta  	 	! insulation thickness, mm
        double precision::eps   		 	! roughness, mm
        double precision::K     		 	! minor loss coefficient, -
        ! network elements
		type(pipe)::pipe_A					! network's pipe A
		type(pipe)::pipe_B					! network's pipe B
    contains
        procedure::create_networks=>double_pipe_create_networks ! create hydraulic pipes
		procedure::get_name=>double_pipe_get_name
        procedure::num_elements=>double_pipe_num_elements ! count the number of hydraulic pipes 
		procedure::to_json=>double_pipe_to_json
		procedure::from_json=>double_pipe_from_json
    end type   
		
	! substation wrapper to create arrays of polymorphic substation objects
	type substation_wrapper
		class(substation),allocatable::item
	end type	
	
	! the progenitor of all substation-like types
    type substation
		type(building_complex),pointer::p_building_complex				! pointer to the building complex (see building_complex_module)
		type(thermal_grid),pointer::p_thermal_grid=>null()				! pointer to the thermal grid object that contains the substation
        type(point)::P0=point(x=0d0,y=0d0,z=0d0)						! coordinates must coincide with start or end point coordinates of double pipes
        character(:),allocatable::name 									! name of the substation
        character(:),allocatable::type_sub 								! type of substation
		type(string),dimension(:),allocatable::building_units       	! array of names of the associated building units
		double precision,dimension(8760)::Q_supply_heat=0d0         	! heat supplied by the substation to the associated building units, Wh
		double precision,dimension(8760)::Q_supply_cool=0d0				! cool supplied by the substation to the associated building units, Wh
		type(dhw_storage)::DHW_tank=dhw_storage()						! the DHW storage tank, with inputs, state of charge and outputs
        type(param),dimension(:),allocatable::parameters				! the list of type-specific parameters
    contains
        procedure::load=>substation_load 								! load parameters in the specific type objects
        procedure::create_networks=>substation_create_networks 			! create connections, pumps, heat exchangers, pipes 
		procedure::get_name=>substation_get_name
        procedure::num_elements=>substation_num_elements 				! return the number of elements (see num_elements type)         
        procedure::initialize=>substation_initialize 					! get heating, cooling and dhw demands from the associated building units
        procedure::solve=>substation_solve 								! calculate energy flows in the time step
		procedure::to_json=>substation_to_json
		procedure::from_json=>substation_from_json
    end type    
	
	! energy centre wrapper to create array of polymorphic energy_centre objects
	type energy_centre_wrapper
		class(energy_centre),allocatable::item
	end type	
	
    ! the progenitor of all energy centre-like types
    type energy_centre
		type(building_complex),pointer::p_building_complex=>null() 		! pointer to the building complex (see building_complex_module)
		type(thermal_grid),pointer::p_thermal_grid=>null() 				! pointer to the thermal grid object that contains the energy centre
        type(point)::P0=point(x=0d0,y=0d0,z=0d0)						! coordinates must coincide with start or end point coordinates of double pipes
        character(:),allocatable::name  								! name of the energy centre
        character(:),allocatable::type_centre 							! type of energy centre
		double precision::last_update=0d0								! the time (in seconds) within the TIMESTEP at which the solve method was last invoked 
        type(param),dimension(:),allocatable::parameters				! the list of type-specific parameters
    contains
        procedure::load=>energy_centre_load								! load parameters in the specific type objects
        procedure::create_networks=>energy_centre_create_networks  		! create connections, pipes, tank
		procedure::get_name=>energy_centre_get_name
        procedure::num_elements=>energy_centre_num_elements 			! return the number of elements (see num_elements type) 
        procedure::initialize=>energy_centre_initialize 				! set the initial parameters
        procedure::solve=>energy_centre_solve 							! calculate energy flows in the time step
        procedure::update=>energy_centre_update 						! calculate outputs in the time step
		procedure::update_output=>energy_centre_update_output			! calculate running average of an output value in the sub-interval [0:time], time<=TIMESTEP
		procedure::to_json=>energy_centre_to_json
		procedure::from_json=>energy_centre_from_json
    end type    		
    
    ! thermal grid 
    type thermal_grid
		! data objects provided in input
		type(incompressible_fluid)::fluid=incompressible_fluid(rho=RHO_WATER,cp=CP_WATER)
        type(double_pipe),dimension(:),allocatable::double_pipes
        type(double_connection),dimension(:),allocatable::double_connections
        type(substation_wrapper),dimension(:),allocatable::substations
        type(energy_centre_wrapper),dimension(:),allocatable::energy_centres
		! data objects internally managed by the method create_networks 
		type(hydraulic_element_wrapper),dimension(:),allocatable::elements ! the array of non channel-like elements		
        type(hydraulic_network)::hydraulic_network ! hydraulic nodes interconnected by channel-like hydraulic elements (pipes, pumps, heat_exchangers) with mass flow rate Vs head loss relationship 
        type(thermal_network)::thermal_network ! fluid thermal nodes exchanging mass through positive defined mass flow rates (advection heat transfer)
		integer::num_elements=0 ! counter of non channel-like elements (all types)
		integer::num_tanks ! number of tank elements (NOTE for now, the tank is the only non channel-like element)
		!double precision::err_T ! the error of the thermal solver as absolute temperature difference between the last two iterations
    contains
        procedure::create_networks=>thermal_grid_create_networks
		procedure::initialize=>thermal_grid_initialize ! calculate incidence matrix for the network elements and nodes
        procedure::add_node=>thermal_grid_add_node ! add a connection node (to the hydraulic network) and return an index 
        procedure::thermal_grid_add_pipe ! add a pipe element (to hydraulic and thermal networks) and return a pointer		
        procedure::thermal_grid_add_pump ! add a pump element (to hydraulic and thermal networks) and return a pointer		
        procedure::thermal_grid_add_heat_exchanger ! add a heat exchanger element (to hydraulic and thermal networks) and return a pointer 	
		procedure::thermal_grid_add_borehole_field ! add a borehole field element (to hydraulic and thermal networks) and return a pointer		
        procedure::thermal_grid_add_buffer_tank ! add a buffer tank element to both the hydraulic and the thermal network 		
        generic::add_element=>thermal_grid_add_pipe,&
							  thermal_grid_add_pump,&
							  thermal_grid_add_heat_exchanger,&
							  thermal_grid_add_borehole_field,&
							  thermal_grid_add_buffer_tank 
        procedure::solve_substations=>thermal_grid_solve_substations ! algebraic methods that calculate energy flows of the substations
        procedure::solve_energy_centres=>thermal_grid_solve_energy_centres ! algebraic methods that calculate energy flows of the energy_centres
        procedure::update_energy_centres=>thermal_grid_update_energy_centres ! algebraic methods that calculate outputs of the energy_centres
        procedure::solve_Hm=>thermal_grid_solve_Hm ! solve the hydraulic problem (H,m) at a given hour using Todini's method
		procedure::solve_Tq=>thermal_grid_solve_Tq ! solve the thermal problem using Euler explicit variable time step
		procedure::calculate_energy_demand=>thermal_grid_calculate_energy_demand ! simulate the thermal grid and calculate yearly energy demand
		procedure::to_json=>thermal_grid_to_json
		procedure::from_json=>thermal_grid_from_json
    end type
    
      
 contains
 
 
! =============================================================================================================
! output (to_json) methods 
! =============================================================================================================

	subroutine connection_to_json(this, json, obj)
		class(connection), intent(in) :: this
		type(json_core), intent(inout)  :: json
		type(json_value), pointer, intent(out) :: obj
		type(json_value), pointer :: obj_P
		call json%create_object(obj,'')
		if (allocated(this%name)) then
		  call json%add(obj,'name',this%name)
		end if
		call json%create_object(obj_P,'P0')
  	    call json%add(obj_P,'x',this%P0%x)
  	    call json%add(obj_P,'y',this%P0%y)
  	    call json%add(obj_P,'z',this%P0%z)
		call json%add(obj,obj_P)
		call save_series_bin(json,obj,'H', this%H)
	end subroutine 
	
	subroutine connection_from_json(this, json, obj, path)
		class(connection) :: this
		type(json_core), intent(inout)  :: json
		type(json_value), pointer, intent(in) :: obj
		character(*)::path
		call json%get(obj,path//'.name',this%name)
		call json%get(obj,path//'.P0.x',this%P0%x)
		call json%get(obj,path//'.P0.y',this%P0%y)
		call json%get(obj,path//'.P0.z',this%P0%z)
		call read_series_bin(json,obj,path//'.H', this%H)
	end subroutine 

	subroutine thermal_node_to_json(this, json, obj)
		class(thermal_node), intent(in) :: this
		type(json_core), intent(inout)  :: json
		type(json_value), pointer, intent(out) :: obj
		call json%create_object(obj,'')
		call save_series_bin(json,obj,'Qin',this%Qin)
		call save_series_bin(json,obj,'T', this%T)
		call json%add(obj,'UA',this%UA)
		call json%add(obj,'MC',this%MC)
	end subroutine 
	
	subroutine thermal_node_from_json(this, json, obj, path)
		class(thermal_node) :: this
		type(json_core), intent(inout)  :: json
		type(json_value), pointer, intent(in) :: obj
		character(*)::path
		call read_series_bin(json,obj,path//'.Qin',this%Qin)
		call read_series_bin(json,obj,path//'.T', this%T)
		call json%get(obj,path//'.MC',this%MC)
		call json%get(obj,path//'.UA',this%UA)
	end subroutine 
	
	subroutine hydraulic_element_to_json(this, json, obj,obj_name)
		class(hydraulic_element),intent(in)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(inout)::obj
		character(*)::obj_name
		type(json_value),pointer::arr,node_obj
		integer::i
		call json%create_object(obj,obj_name)
		if (allocated(this%name)) then
		  call json%add(obj,'name',this%name)
		end if
		call json%add(obj,'MC', this%MC)
		call json%add(obj,'UA',this%UA)
		if (associated(this%thermal_nodes)) then
		  call json%create_array(arr,'thermal_nodes')
		  do i = 1, size(this%thermal_nodes)
			 call this%thermal_nodes(i)%to_json(json,node_obj)
			 call json%add(arr,node_obj)
		  end do
		  call json%add(obj,arr)
		end if
	end subroutine 

	subroutine hydraulic_element_from_json(this, json, obj,path)
		class(hydraulic_element)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(in)::obj
		character(*)::path
		type(json_value),pointer::arr,node_obj
		integer::i,n_elements
		logical::found
		call json%get(obj,path//'.name',this%name)	
		call json%get(obj,path//'.MC',this%MC)	
		call json%get(obj,path//'.UA',this%UA)	
		call json%info(obj,path//'.thermal_nodes',found,n_children=n_elements)
		allocate(this%thermal_nodes(n_elements))
	    do i = 1, size(this%thermal_nodes)
			call this%thermal_nodes(i)%from_json(json,obj,path//'.thermal_nodes('//integer_to_text(i)//')')
	    end do
	end subroutine 

	subroutine hydraulic_channel_to_json(this, json, obj,obj_name)
		class(hydraulic_channel),intent(in)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(inout)::obj
		character(*)::obj_name
		call hydraulic_element_to_json(this,json,obj,obj_name)
		call json%add(obj,'node_in', this%node_in)
		call json%add(obj,'node_out',this%node_out)
		call save_series_bin(json,obj,'mdot',this%mdot)
	end subroutine 

	subroutine hydraulic_channel_from_json(this, json, obj,path)
		class(hydraulic_channel)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(in)::obj
		character(*)::path
		call hydraulic_element_from_json(this, json, obj,path)
		call json%get(obj,path//'.node_in',this%node_in)	
		call json%get(obj,path//'.node_out',this%node_out)	
		call read_series_bin(json,obj,path//'.mdot',this%mdot)
	end subroutine 

	subroutine pipe_to_json(this, json, obj,obj_name)
		class(pipe),intent(in)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(inout)::obj
		type(json_value),pointer::obj_pipe
		character(*)::obj_name
		call hydraulic_channel_to_json(this,json, obj_pipe,obj_name)
		call json%add(obj_pipe,'D',this%D)	
		call json%add(obj_pipe,'L',this%L)	
		call json%add(obj_pipe,'K',this%K)	
		call json%add(obj_pipe,'eps',this%eps)	
		call json%add(obj_pipe,'ins_delta',this%ins_delta)	
		call json%add(obj_pipe,'ins_lambda',this%ins_lambda)
		call json%add(obj,obj_pipe)
	end subroutine 

	subroutine pipe_from_json(this, json, obj,path)
		class(pipe)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(in)::obj
		character(*)::path
		type(json_value),pointer::arr,node_obj
		integer::i
		call hydraulic_channel_from_json(this,json, obj,path)
		call json%get(obj,path//'.D',this%D)	
		call json%get(obj,path//'.L',this%L)	
		call json%get(obj,path//'.K',this%K)	
		call json%get(obj,path//'.eps',this%eps)	
		call json%get(obj,path//'.ins_delta',this%ins_delta)	
		call json%get(obj,path//'.ins_lambda',this%ins_lambda)	
    end subroutine 

	subroutine pipe_get_par(this, par)
		class(pipe)::this
		type(param)::par
  		select case(par%key(2))
			case ('D')
				this%D=par%val_dbl
			case ('L')
				this%L=par%val_dbl
			case ('K')
				this%K=par%val_dbl
			case ('eps')	
				this%eps=par%val_dbl
			case ('ins_delta')	
				this%ins_delta=par%val_dbl
			case ('ins_lambda')	
				this%ins_lambda=par%val_dbl
		end select      
	end subroutine 


    subroutine pump_to_json(this, json, obj,obj_name)
		class(pump),intent(in)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(inout)::obj
		type(json_value),pointer::obj_pump
		character(*)::obj_name
		call hydraulic_channel_to_json(this,json, obj_pump,obj_name)
		call save_series_bin(json,obj_pump,'W',this%W)		
		call json%add(obj_pump,'a1',this%a1)	
		call json%add(obj_pump,'a2',this%a2)	
		call json%add(obj_pump,'a3',this%a3)
		call json%add(obj,obj_pump)		
	end subroutine 

	subroutine pump_from_json(this,json,obj,path)
		class(pump)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(in)::obj
		character(*)::path
		call hydraulic_channel_from_json(this,json, obj,path)
		call read_series_bin(json,obj,path//'.W',this%W)
		call json%get(obj,path//'.a1',this%a1)	
		call json%get(obj,path//'.a2',this%a2)	
		call json%get(obj,path//'.a3',this%a3)	
    end subroutine 
    
    subroutine pump_get_par(this, par)
		class(pump)::this
		type(param)::par
		select case(par%key(2))
			case ('a1')
				this%a1 = par%val_dbl
			case ('a2')
				this%a2 = par%val_dbl
			case ('a3')
				this%a3 = par%val_dbl
		end select
	end subroutine
	
	subroutine heat_exchanger_to_json(this, json, obj,obj_name)
		class(heat_exchanger),intent(in)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(inout)::obj
		type(json_value),pointer::obj_hx
		character(*)::obj_name
		call hydraulic_channel_to_json(this,json, obj_hx,obj_name)
		call json%add(obj_hx,'DP_nom',this%DP_nom)
		call json%add(obj_hx,'m_nom',this%m_nom)
		call json%add(obj_hx,'Q_nom',this%Q_nom)
		call json%add(obj,obj_hx)
	end subroutine 

	subroutine heat_exchanger_from_json(this, json,obj,path)
		class(heat_exchanger)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(in)::obj
		character(*)::path
		call hydraulic_channel_from_json(this,json, obj,path)
		call json%get(obj,path//'.DP_nom',this%DP_nom)			
		call json%get(obj,path//'.m_nom',this%m_nom)			
		call json%get(obj,path//'.Q_nom',this%Q_nom)			
    end subroutine 
    
    subroutine heat_exchanger_get_par(this, par)
		class(heat_exchanger)::this
		type(param)::par
		select case(par%key(2))
			case ('DP_nom')
				this%DP_nom = par%val_dbl
			case ('m_nom')
				this%m_nom = par%val_dbl
			case ('Q_nom')
				this%Q_nom = par%val_dbl
		end select
    end subroutine
    
    function heat_exchanger_get_aux_h(this,hour) result(Q_aux)
		class(heat_exchanger)::this
        integer::hour,prev_hour
        double precision::Q_aux
        prev_hour=hour-1
        if (prev_hour==0) prev_hour=8760 
        if (this%thermal_nodes(1)%T(hour)<this%T_aux_h) then
            Q_aux=this%mdot(prev_hour)*this%fluid%cp*(this%T_aux_h-this%thermal_nodes(1)%T(hour))
            this%thermal_nodes(1)%T(hour)=this%T_aux_h
            this%thermal_nodes(1)%T(prev_hour)=this%T_aux_h
        else
            Q_aux=0d0
        end if
    end function    

    function heat_exchanger_get_aux_c(this,hour) result(Q_aux)
		class(heat_exchanger)::this
        integer::hour,prev_hour
        double precision::Q_aux
        prev_hour=hour-1
        if (prev_hour==0) prev_hour=8760 
        if (this%thermal_nodes(1)%T(hour)>this%T_aux_c) then
            Q_aux=this%mdot(prev_hour)*this%fluid%cp*(this%T_aux_c-this%thermal_nodes(1)%T(hour))
            this%thermal_nodes(1)%T(hour)=this%T_aux_c
            this%thermal_nodes(1)%T(prev_hour)=this%T_aux_c
        else
            Q_aux=0d0
        end if
    end function     

	subroutine borehole_field_to_json(this, json, obj,obj_name)
		class(borehole_field),intent(in)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(inout)::obj
		type(json_value),pointer::obj_bf
		character(*)::obj_name
		type(json_value),pointer::obj_bhe
		logical::found
		call heat_exchanger_to_json(this,json, obj,obj_name)
		call json%get(obj,obj_name,obj_bf,found)
		call this%bhe%to_json(json,obj_bf,'bhe')
	end subroutine 

	subroutine borehole_field_from_json(this, json, obj,path)
		class(borehole_field)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(in)::obj
		character(*)::path
		call heat_exchanger_from_json(this,json, obj,path)
		call this%bhe%from_json(json,obj,path//'.bhe')		
    end subroutine 

	subroutine borehole_field_get_par(this,par)
		class(borehole_field)::this
		type(param)::par
		select case(par%key(2))
		case('DP_nom')
			this%DP_nom=par%val_dbl			
		case('m_nom')
			this%m_nom=par%val_dbl			
		case('bhe')
			select case(par%key(3))
				case('bhe_XY_coordinates')
					this%bhe%bhe_XY_coordinates=par%mat_dbl
				case('N_bhe')
					this%bhe%N_bhe=par%val_int
				case('layout_mode')
					this%bhe%layout_mode=par%val_int
				case('bhe_type')
					this%bhe%bhe_type=par%val_int
				case('D_bhe')
					this%bhe%D_bhe=par%val_dbl
				case('H')
					this%bhe%H=par%val_dbl
				case('D_p')
					this%bhe%D_p=par%val_dbl
				case('s_p')
					this%bhe%s_p=par%val_dbl
				case('D')
					this%bhe%D=par%val_dbl
				case('lambda_p')
					this%bhe%lambda_p=par%val_dbl
				case('lambda_grout')
					this%bhe%lambda_grout=par%val_dbl
				case('lambda')
					this%bhe%lambda=par%val_dbl
				case('C')
					this%bhe%C=par%val_dbl
				case('Q_ghd')
					this%bhe%Q_ghd=par%val_dbl
				case('Q_gcd')
					this%bhe%Q_gcd=par%val_dbl
				case('Q_ghm')
					this%bhe%Q_ghm=par%val_dbl
				case('Q_gcm')
					this%bhe%Q_gcm=par%val_dbl
				case('Q_ga')
					this%bhe%Q_ga=par%val_dbl
				case('B_des')
					this%bhe%B_des=par%val_dbl
				case('H_des')
					this%bhe%H_des=par%val_dbl
				case('T_fih_design')
					this%bhe%T_fih_design=par%val_dbl
				case('T_foh_design')
					this%bhe%T_foh_design=par%val_dbl
				case('T_fic_design')
					this%bhe%T_fic_design=par%val_dbl
				case('T_foc_design')
					this%bhe%T_foc_design=par%val_dbl
			end select     
		end select				
	end subroutine    
	    
	subroutine tank_to_json(this, json, obj,obj_name)
		class(tank),intent(in)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(inout)::obj
		character(*)::obj_name
		type(json_value),pointer::obj_tank
		call hydraulic_element_to_json(this,json, obj_tank,obj_name)
		call json%add(obj_tank,'K_in',this%K_in)
		call json%add(obj_tank,'K_out',this%K_out)
		call json%add(obj_tank,'num_thermal_nodes',this%num_thermal_nodes)
		call json%add(obj,obj_tank)
	end subroutine 

	subroutine tank_from_json(this, json, obj,path)
		class(tank)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(in)::obj
		character(*)::path
		call hydraulic_element_from_json(this,json, obj,path)
		call json%get(obj,path//'.K_in',this%K_in)			
		call json%get(obj,path//'.K_out',this%K_out)			
		call json%get(obj,path//'.num_thermal_nodes',this%num_thermal_nodes)
    end subroutine 

   	subroutine buffer_tank_to_json(this, json, obj,obj_name)
		class(buffer_tank),intent(in)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(inout)::obj
		character(*)::obj_name
		type(json_value),pointer::obj_tank
		call hydraulic_element_to_json(this,json, obj_tank,obj_name)
		call json%add(obj_tank,'K_in',this%K_in)
		call json%add(obj_tank,'K_out',this%K_out)
		call json%add(obj_tank,'num_thermal_nodes',this%num_thermal_nodes)
        call this%controller%to_json(json,obj_tank,'controller')
 		call json%add(obj,obj_tank)
	end subroutine 

	subroutine buffer_tank_from_json(this, json, obj,path)
		class(buffer_tank)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(in)::obj
		character(*)::path
		call hydraulic_element_from_json(this,json, obj,path)
		call json%get(obj,path//'.K_in',this%K_in)			
		call json%get(obj,path//'.K_out',this%K_out)			
		call json%get(obj,path//'.num_thermal_nodes',this%num_thermal_nodes)					
        call this%controller%from_json(json,obj,path//'.controller')	
    end subroutine 
    
     subroutine buffer_tank_get_par(this,par)
		class(buffer_tank)::this
        type(param)::par
		select case(par%key(2))
			case ('volume')
				this%MC=par%val_dbl*RHO_WATER*CP_WATER
			case ('UA')
				this%UA=par%val_dbl
			case ('K_in')
				this%K_in=par%val_dbl
			case ('K_out')
				this%K_out=par%val_dbl
			case ('controller')
			select case(par%key(3))
				case('T_storage_max')
					this%controller%T_storage_max=par%val_dbl
				case('T_storage_min')
					this%controller%T_storage_min=par%val_dbl
            end select        
		end select	
	end subroutine    
   
	subroutine double_connection_to_json(this, json, obj)
		class(double_connection), intent(in) :: this
		type(json_core), intent(inout)  :: json
		type(json_value), pointer, intent(out) :: obj
		call json%create_object(obj,'')
		if (allocated(this%name)) then
		  call json%add(obj,'name',this%name)
		end if
		call json%add(obj,'idx_node_A',this%idx_node_A)
		call json%add(obj,'idx_node_B',this%idx_node_B)
	end subroutine 

	subroutine double_connection_from_json(this, json, obj,path)
		class(double_connection)::this
		type(json_core), intent(inout)  :: json
		type(json_value), pointer, intent(in) :: obj
		character(*)::path
		call json%get(obj,path//'.idx_node_A',this%idx_node_A)
		call json%get(obj,path//'.idx_node_B',this%idx_node_B)
	end subroutine 
	
	subroutine double_pipe_to_json(this, json, obj)
		class(double_pipe), intent(in) :: this
		type(json_core), intent(inout)  :: json
		type(json_value), pointer, intent(out) :: obj
		type(json_value), pointer :: pipe_obj
		call json%create_object(obj,'')
		if (allocated(this%name)) then
		  call json%add(obj,'name',this%name)
		end if
		call this%pipe_A%to_json(json,obj,'pipe_A')
		call this%pipe_B%to_json(json,obj,'pipe_B')
	end subroutine 

	subroutine double_pipe_from_json(this, json, obj,path)
		class(double_pipe)::this
		type(json_core), intent(inout)  :: json
		type(json_value), pointer, intent(in) :: obj
		character(*)::path
		call this%pipe_A%from_json(json,obj,path//'.pipe_A')
		call this%pipe_B%from_json(json,obj,path//'.pipe_B')
	end subroutine 
	
	function double_pipe_get_name(this,str) result (txt)
		class(double_pipe)::this
		character(*)::str
		character(:),allocatable::txt
		txt=this%name//'.'//str
	end function

	function double_connection_get_name(this,str) result (txt)
		class(double_connection)::this
		character(*)::str
		character(:),allocatable::txt
		txt=this%name//'.'//str
	end function

	function substation_get_name(this,str) result (txt)
		class(substation)::this
		character(*)::str
		character(:),allocatable::txt
		txt=this%name//'.'//str
	end function

	function energy_centre_get_name(this,str) result (txt)
		class(energy_centre)::this
		character(*)::str
		character(:),allocatable::txt
		txt=this%name//'.'//str
	end function

	subroutine substation_to_json(this,json,obj)
		class(substation)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(out)::obj
		call json%create_object(obj,'')
		if (allocated(this%name)) then
		  call json%add(obj,'name',this%name)
		end if
		if (allocated(this%type_sub)) then
		  call json%add(obj,'type_sub',this%type_sub)
		end if
		call save_series_bin(json,obj,'Q_supply_heat',this%Q_supply_heat)
		call save_series_bin(json,obj,'Q_supply_cool',this%Q_supply_cool)
		call this%DHW_tank%to_json(json,obj,'DHW_tank') 
	end subroutine	
	
	subroutine substation_from_json(this,json,obj,path)
		class(substation)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(in)::obj
		character(*)::path		
		call json%get(obj,path//'.name',this%name)
		call json%get(obj,path//'.type_sub',this%type_sub)
		call read_series_bin(json,obj,path//'.Q_supply_heat',this%Q_supply_heat)
		call read_series_bin(json,obj,path//'.Q_supply_cool',this%Q_supply_cool)
		call this%DHW_tank%from_json(json,obj,path//'.DHW_tank') 
	end subroutine	
		
	subroutine energy_centre_to_json(this,json,obj)
		class(energy_centre)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(out)::obj
		call json%create_object(obj,'')		
		if (allocated(this%name)) then
		  call json%add(obj,'name',this%name)
		end if		
		if (allocated(this%type_centre)) then
		  call json%add(obj,'type_centre',this%type_centre)
		end if
	end subroutine	

	subroutine energy_centre_from_json(this,json,obj,path)
		class(energy_centre)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(in)::obj
		character(*)::path
		call json%get(obj,path//'.name',this%name)
		call json%get(obj,path//'.type_centre',this%type_centre)
	end subroutine	

	subroutine thermal_grid_to_json(this, json, obj)
		class(thermal_grid),intent(in)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(out)::obj
		type(json_value),pointer::connection_arr,connection_obj
		type(json_value),pointer::double_connection_arr,double_connection_obj
		type(json_value),pointer::double_pipe_arr,double_pipe_obj
		type(json_value),pointer::substation_arr,substation_obj
		type(json_value),pointer::energy_centre_arr,energy_centre_obj
		integer::i
		call json%create_object(obj,'')
		call json%create_array(connection_arr,'connections')
		do i=1,size(this%hydraulic_network%node)
			call this%hydraulic_network%node(i)%to_json(json,connection_obj)
			call json%add(connection_arr,connection_obj)
		end do	
		call json%create_array(double_connection_arr,'double_connections')
		do i=1,size(this%double_connections)
			call this%double_connections(i)%to_json(json,double_connection_obj)
			call json%add(double_connection_arr,double_connection_obj)
		end do	
		call json%create_array(double_pipe_arr,'double_pipes')
		do i=1,size(this%double_pipes)
			call this%double_pipes(i)%to_json(json,double_pipe_obj)
			call json%add(double_pipe_arr,double_pipe_obj)
		end do	
		call json%create_array(substation_arr,'substations')
		do i=1,size(this%substations)
			call this%substations(i)%item%to_json(json,substation_obj)
			call json%add(substation_arr,substation_obj)
		end do	
		call json%create_array(energy_centre_arr,'energy_centres')
		do i=1,size(this%energy_centres)
			call this%energy_centres(i)%item%to_json(json,energy_centre_obj)
			call json%add(energy_centre_arr,energy_centre_obj)
		end do	
		call this%fluid%to_json(json,obj,'fluid')
		call json%add(obj,connection_arr)
		call json%add(obj,double_connection_arr)
		call json%add(obj,double_pipe_arr)
		call json%add(obj,substation_arr)
		call json%add(obj,energy_centre_arr)		
	end subroutine 

	subroutine thermal_grid_from_json(this, json,obj)
		class(thermal_grid)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(inout)::obj
		integer::i,n_elements
		logical::found
		call json%info(obj,'connections',found, n_children=n_elements)
		allocate(this%hydraulic_network%node(n_elements))
		do i=1,n_elements
			call this%hydraulic_network%node(i)%from_json(json,obj,'connections('//integer_to_text(i)//')')
		end do	
		call json%info(obj,'double_connections',found, n_children=n_elements)
		do i=1,n_elements
			call this%double_connections(i)%from_json(json,obj,'double_connections('//integer_to_text(i)//')')
		end do	
		call json%info(obj,'double_pipes',found, n_children=n_elements)
		do i=1,n_elements
			call this%double_pipes(i)%from_json(json,obj,'double_pipes('//integer_to_text(i)//')')
		end do	
		call json%info(obj,'substations',found, n_children=n_elements)
		do i=1,n_elements
			call this%substations(i)%item%from_json(json,obj,'substations('//integer_to_text(i)//')')
		end do	
		call json%info(obj,'energy_centres',found, n_children=n_elements)
		do i=1,n_elements
			call this%energy_centres(i)%item%from_json(json,obj,'energy_centres('//integer_to_text(i)//')')
		end do	
		call this%fluid%from_json(json,obj,'fluid')
	end subroutine 

! =============================================================================================================
! core methods 
! =============================================================================================================
        				
    subroutine polyline_set_properties(this)
        class(polyline)::this
        integer::i
        this%length=0d0
        do i=2,size(this%points)
            this%length=this%length+((this%points(i)%X-this%points(i-1)%X)**2+ &
                                     (this%points(i)%Y-this%points(i-1)%Y)**2+ &
                                     (this%points(i)%Z-this%points(i-1)%Z)**2)**0.5
        end do    
        this%P1=this%points(1)
        this%P2=this%points(size(this%points))
    end subroutine
    	
   	function add_num_elements(num_1,num_2) 
		type(num_elements),intent(in)::num_1,num_2
		type(num_elements)::add_num_elements
		add_num_elements%num_pipes=num_1%num_pipes+num_2%num_pipes
		add_num_elements%num_heat_exchangers=num_1%num_heat_exchangers+num_2%num_heat_exchangers
		add_num_elements%num_pumps=num_1%num_pumps+num_2%num_pumps
		add_num_elements%num_tanks=num_1%num_tanks+num_2%num_tanks
		add_num_elements%num_thermal_nodes=num_1%num_thermal_nodes+num_2%num_thermal_nodes
	end function
	
	function num_elements_hydraulic_channels_like(this) result(num)
		class(num_elements)::this
		integer::num
		num=this%num_pipes+this%num_heat_exchangers+this%num_pumps
	end function	

	function num_elements_non_hydraulic_channels_like(this) result(num)
		class(num_elements)::this
		integer::num
		num=this%num_tanks
	end function	
	
	subroutine num_elements_set_num_thermal_nodes(this)
		class(num_elements)::this
		this%num_thermal_nodes=TANK_NODES*this%num_tanks+this%num_pipes+this%num_heat_exchangers+this%num_pumps 		
	end subroutine

	subroutine hydraulic_channel_connect(this,node_in,node_out)
		class(hydraulic_channel)::this
		integer,optional::node_in,node_out
		if (present(node_in)) this%node_in=node_in
		if (present(node_out))this%node_out=node_out
	end subroutine

	subroutine hydraulic_channel_initialize(this)
		class(hydraulic_channel)::this
	end subroutine

	subroutine tank_connect(this,hydraulic_node,element_to_top,element_to_base)
		class(tank)::this
		class(hydraulic_channel),target,optional::element_to_top
		class(hydraulic_channel),target,optional::element_to_base
		integer,optional::hydraulic_node
		if (present(hydraulic_node)) this%hydraulic_node=hydraulic_node
		if (present(element_to_top))this%element_to_top=>element_to_top
		if (present(element_to_base))this%element_to_base=>element_to_base
	end subroutine

	subroutine pipe_set(this,UA,MC,D,L,K,eps,ins_delta,ins_lambda)
		class(pipe)::this
        double precision,optional::UA                     ! overall heat loss coefficient (W/K)
        double precision,optional::MC                     ! overall thermal capacity (J/K)
        double precision,optional::D                      ! internal diameter (m)
        double precision,optional::L                      ! length (m)
        double precision,optional::K                      ! minor pressure drops coefficient (-)
        double precision,optional::eps                    ! relative roughness, pipe roughness/pipe diameter (-)
		double precision,optional::ins_delta			  ! insulation thickness (m)
		double precision,optional::ins_lambda			  ! insulation thermal conductivity (W/(m,K))
		if (present(UA)) this%UA=UA
		if (present(MC)) this%MC=MC
		if (present(D)) this%D=D
		if (present(L)) this%L=L
		if (present(K)) this%K=K
		if (present(eps)) this%eps=eps
		if (present(ins_delta)) this%ins_delta=ins_delta
		if (present(ins_lambda)) this%ins_lambda=ins_lambda
 	end subroutine

	subroutine heat_exchanger_set(this,UA,MC,DP_nom,Q_nom,m_nom)
		class(heat_exchanger)::this
        double precision,optional::UA                     ! overall heat loss coefficient (W/K)
        double precision,optional::MC                     ! overall thermal capacity (J/K)
	    double precision,optional::DP_nom                 ! nominal pressure drop (Pa)
        double precision,optional::m_nom                  ! nominal mass flow rate (kg/s)
        double precision,optional::Q_nom                  ! nominal heat duty (W)
		if (present(UA)) this%UA=UA
		if (present(MC)) this%MC=MC
		if (present(DP_nom)) this%DP_nom=DP_nom
		if (present(m_nom)) this%m_nom=m_nom
		if (present(Q_nom)) this%Q_nom=Q_nom
	end subroutine

	subroutine pump_set(this,UA,MC,a1,a2,a3)
		class(pump)::this
        double precision,optional::UA                     ! overall heat loss coefficient (W/K)
        double precision,optional::MC                     ! overall thermal capacity (J/K)
        double precision,optional::a1                     ! coefficient a1 in H=a1+a2/m+a3/m^2, (m.w.c.)
        double precision,optional::a2                     ! coefficient a2 in H=a1+a2/m+a3/m^2, (m.w.c.*(kg/s))
        double precision,optional::a3                     ! coefficient a3 in H=a1+a2/m+a3/m^2, (m.w.c.*(kg/s)^2)
		if (present(UA)) this%UA=UA
		if (present(MC)) this%MC=MC
		if (present(a1)) this%a1=a1
		if (present(a2)) this%a2=a2
		if (present(a3)) this%a3=a3
	end subroutine

	subroutine borehole_field_set(this,UA,MC,bhe)
		class(borehole_field)::this
        double precision,optional::UA                     ! overall heat loss coefficient (W/K)
        double precision,optional::MC                     ! overall thermal capacity (J/K)
		type(borehole_heat_exchanger),optional::bhe       ! the borehole heat exchanger object with thermal data and methods
		if (present(UA)) this%UA=UA
		if (present(MC)) this%MC=MC
		if (present(bhe)) this%bhe=bhe
	end subroutine
	
	subroutine buffer_tank_set(this,UA,MC,K_in,K_out,controller)
		class(buffer_tank)::this
        double precision,optional::UA                     ! overall heat loss coefficient (W/K)
        double precision,optional::MC                     ! overall thermal capacity (J/K)
        double precision,optional::K_in                   ! minor pressure drops coefficient at tank inlet (-)
        double precision,optional::K_out                  ! minor pressure drops coefficient at tank outlet (-)
        type(buffer_controller),optional::controller      ! controller
		if (present(UA)) this%UA=UA
		if (present(MC)) this%MC=MC
		if (present(K_in)) this%K_in=K_in
		if (present(K_out)) this%K_out=K_out
		if (present(controller)) this%controller=controller
	end subroutine
		
	function heat_exchanger_inflow_T(this,hour) result(T)
		class(heat_exchanger)::this
		integer::hour
		double precision::T,m
        integer::i
		T=0d0
		m=0d0
		do i=1,size(this%inflows)
			T=T+this%inflows(i)%T(hour)*this%inflows(i)%mdot
			m=m+this%inflows(i)%mdot
        end do
        if (m>0) T=T/m	
	end function	
	
	function tank_mass_flow_in_top(this,hour) result(mdot)
		class(tank)::this
		integer::hour
		double precision::mdot
		mdot=this%element_to_top%mdot(hour)
		if (this%hydraulic_node==this%element_to_top%node_in) mdot=-1d0*mdot
	end function

	function tank_T_middle(this,hour) result (T)
		class(tank)::this
		integer::hour
		double precision::T
		integer::central_node
		central_node=int((this%num_thermal_nodes+1d0)/2d0)
		T=this%thermal_nodes(central_node)%T(hour)
	end function

	function tank_T_top(this,hour) result (T)
		class(tank)::this
		integer::hour
		double precision::T
		T=this%thermal_nodes(1)%T(hour)
	end function

	function tank_T_base(this,hour) result (T)
		class(tank)::this
		integer::hour
		double precision::T
		T=this%thermal_nodes(this%num_thermal_nodes)%T(hour)
	end function

	subroutine tank_charge_base(this,hour,T_fluid,mdot)
		class(tank)::this
		integer::hour
		double precision::T_fluid,mdot
		double precision::T_in
		integer::i
		do i=this%num_thermal_nodes,1,-1
			if (i==this%num_thermal_nodes) then
				T_in=T_fluid
			else
				T_in=this%thermal_nodes(i+1)%T(hour)
			end if
			this%thermal_nodes(i)%Qin(hour)=mdot*this%fluid%cp*(T_in-this%thermal_nodes(i)%T(hour))
		end do	
	end subroutine	
	
	subroutine tank_charge_top(this,hour,T_fluid,mdot)
		class(tank)::this
		integer::hour
		double precision::T_fluid,mdot
		double precision::T_in
		integer::i
		do i=1,this%num_thermal_nodes
			if (i==1) then
				T_in=T_fluid
			else
				T_in=this%thermal_nodes(i-1)%T(hour)
			end if
			this%thermal_nodes(i)%Qin(hour)=mdot*this%fluid%cp*(T_in-this%thermal_nodes(i)%T(hour))
		end do							
    end subroutine	

	subroutine tank_charge_zero(this,hour)
		class(tank)::this
		integer::hour
		integer::i
		do i=1,this%num_thermal_nodes
			this%thermal_nodes(i)%Qin(hour)=0d0
		end do	
	end subroutine	

    function tank_Q_advection(this,i,hour_mdot,hour_T) result(Qdot)
		class(tank)::this
		integer::hour_mdot,hour_T
		double precision::Qdot
		integer::i
		double precision::mdot,Tin
		mdot=this%mass_flow_in_top(hour_mdot)
		if (mdot>0) then
			if (i==1) then
				Tin=this%element_to_top%thermal_nodes(1)%T(hour_T)
			else
				Tin=this%thermal_nodes(i-1)%T(hour_T)
			end if	
			Qdot=mdot*this%fluid%cp*(Tin-this%thermal_nodes(i)%T(hour_T))
		else
			if (i==this%num_thermal_nodes) then
				Tin=this%element_to_base%thermal_nodes(1)%T(hour_T)
			else
				Tin=this%thermal_nodes(i+1)%T(hour_T)
			end if	
			Qdot=abs(mdot)*this%fluid%cp*(Tin-this%thermal_nodes(i)%T(hour_T))
		end if
	end function

	function tank_Q_loss(this,i,hour) result(Qdot)
		class(tank)::this
		integer::hour
		double precision::Qdot
		integer::i
		Qdot=this%UA*(this%thermal_nodes(i)%T(hour)-this%thermal_nodes(i)%T_env(hour))
	end function
	
	
    subroutine thermal_grid_initialize(this)
        class(thermal_grid),target::this
        integer::i,j,k,num_mass_flows,num_mass_flows_per_node
		associate(h_n=>this%hydraulic_network,t_n=>this%thermal_network)	
		! calculate incidence matrix for the hydraulic network
		allocate(h_n%Aa_hyd(size(h_n%node),size(h_n%elements)))
		do i=1,size(h_n%node)      
			do j=1,size(h_n%elements)
				if(h_n%elements(j)%item%node_in==i) then
					h_n%Aa_hyd(i,j)=1d0
				else if(h_n%elements(j)%item%node_out==i) then
					h_n%Aa_hyd(i,j)=-1d0
				else
					h_n%Aa_hyd(i,j)=0d0
				end if
			end do
		end do
	    ! settings for the reference node (expansion vessel), by convention node 1 is selected
        h_n%node(1)%fix_height(1:8760)=1 ! flag
		h_n%node(1)%H(1:8760)=30d0 ! m.w.c.
		! mapping between hydraulic network mass flows and thermal network mass flows (defined positive)
		! for each node, the array of channels connected to the node and their signs are set
        num_mass_flows=0
		do i=1,size(h_n%node)   
			num_mass_flows_per_node=sum(abs(h_n%Aa_hyd(i,:)))
			allocate(h_n%node(i)%channels(num_mass_flows_per_node))
			allocate(h_n%node(i)%signs(num_mass_flows_per_node))
			allocate(h_n%node(i)%a(num_mass_flows_per_node))
			allocate(h_n%node(i)%b(num_mass_flows_per_node))
			! two mass flows from channels 1 and 2 => two positive flows (1->2, 2->1)
			! three mass flows from channels 1, 2, and 3 => six positive flows (1->2, 1->3,2->1,2->3,3->1,3->2)
			! and so on, thus if n is the number of channels, n*(n-1) is the number of positive mass flows among them
			num_mass_flows=num_mass_flows+num_mass_flows_per_node*(num_mass_flows_per_node-1)
			k=0
			do j=1,size(h_n%elements)
				if(h_n%Aa_hyd(i,j)==1d0) then
					k=k+1
					h_n%node(i)%channels(k)=j
					h_n%node(i)%signs(k)=-1d0
				else if(h_n%Aa_hyd(i,j)==-1d0) then
					k=k+1
					h_n%node(i)%channels(k)=j
					h_n%node(i)%signs(k)=+1d0
                end if
            end do    
		end do
		! allocate mass_flows of the thermal_network, adding mass_flow inside tanks + mass flows at tank inlet and outlet 
		! (only 2 because the other 2 are taken from the mass flows passing through the tank's hydraulic node)
		num_mass_flows=num_mass_flows+this%num_tanks*(2*(TANK_NODES-1)+2)
		allocate(t_n%mass_flows(num_mass_flows))
		! allocate A,B,T,C,mdotcp,D matrices
		allocate(t_n%A(t_n%num_thermal_nodes,t_n%num_thermal_nodes))
		allocate(t_n%B(t_n%num_thermal_nodes,1))
		allocate(t_n%T(1,t_n%num_thermal_nodes))
		allocate(t_n%C(t_n%num_thermal_nodes))
		allocate(t_n%mdotcp(t_n%num_thermal_nodes))
		allocate(t_n%D(t_n%num_thermal_nodes))
		allocate(t_n%solver_mode(t_n%num_thermal_nodes))
		allocate(t_n%algebraic_idx(t_n%num_thermal_nodes))
		do i=1,t_n%num_thermal_nodes
			t_n%T(1,i)=t_n%thermal_nodes(i)%T(1)
		end do
		end associate	
     end subroutine	
            
     ! create hydraulics nodes for the double_connection
     subroutine double_connection_create_networks(this)
	    class(double_connection)::this
        associate(thermal_grid=>this%p_thermal_grid)
		this%idx_node_A=thermal_grid%add_node(name=this%get_name('node_A'),P0=this%P0)
		this%idx_node_B=thermal_grid%add_node(name=this%get_name('node_B'),P0=this%P0%translate((/0d0,0d0,1d0/))) 
		end associate		
     end subroutine

     ! create hydraulics pipes for the double_pipe
     subroutine double_pipe_create_networks(this)
	    class(double_pipe)::this
        type(hydraulic_network),pointer::hydraulic_network
        integer::idx_node_A1,idx_node_A2 ! index of conventional start node (A1) and end node (A2) of pipe A
        integer::idx_node_B1,idx_node_B2 ! index of conventional start node (B1) and end node (B2) of pipe B
		integer::idx
        associate(thermal_grid=>this%p_thermal_grid)
        ! set nodes
		idx_node_A1=thermal_grid%add_node(name=this%get_name('pipe_A1'),P0=this%path%P1)
		idx_node_A2=thermal_grid%add_node(name=this%get_name('pipe_A2'),P0=this%path%P2)
		idx_node_B1=thermal_grid%add_node(name=this%get_name('pipe_B1'),P0=this%path%P1%translate((/0d0,0d0,1d0/)))
		idx_node_B2=thermal_grid%add_node(name=this%get_name('pipe_B2'),P0=this%path%P2%translate((/0d0,0d0,1d0/)))
		! set pipes parameters
		if (associated(this%pipe_A%thermal_nodes)) then
			deallocate(this%pipe_A%thermal_nodes)
		end if	
		this%pipe_A=pipe(name=this%get_name('pipe_A'),node_in=0,node_out=0,D=1d-3*this%D_int,L=this%path%length,K=this%k,eps=this%eps/this%D_int, &
			ins_delta=1d-3*this%ins_delta,ins_lambda=this%ins_lambda)
		if (associated(this%pipe_B%thermal_nodes)) deallocate(this%pipe_B%thermal_nodes)
		this%pipe_B=pipe(name=this%get_name('pipe_B'),node_in=0,node_out=0,D=1d-3*this%D_int,L=this%path%length,K=this%k,eps=this%eps/this%D_int, &
			ins_delta=1d-3*this%ins_delta,ins_lambda=this%ins_lambda)
        ! connect pipes
        call this%pipe_A%connect(node_in=idx_node_A1,node_out=idx_node_A2)
		call this%pipe_B%connect(node_in=idx_node_B1,node_out=idx_node_B2)
        ! add pipes to the thermal grid
		call thermal_grid%add_element(this%pipe_A)	
		call thermal_grid%add_element(this%pipe_B)	
		end associate
     end subroutine
	 
     ! count number of hydraulic elements
     function double_pipe_num_elements(this) result(num)
	        class(double_pipe)::this
	        type(num_elements)::num
	        num%num_pipes=2
			call num%set_num_thermal_nodes()
     end function
	 
	 ! substation 000 (ancestor) dummy methods
     subroutine substation_load(this)
	        class(substation)::this
	 end subroutine		
     subroutine substation_create_networks(this)
	    class(substation)::this
     end subroutine
	 function substation_num_elements(this) result(num)
	        class(substation)::this
			type(num_elements)::num
     end function
     subroutine substation_initialize(this)
	        class(substation),intent(inout)::this
			integer::bu,un,j,num
			! get heating, cooling, and DHW demands from the associated building units
			! NOTE DHW systems in parallel, so Volume, UA and Q_heater are the aggregated (summed) values
			! taken from the connected building units
			associate(buildings=>this%p_building_complex%buildings,&
			    Q_supply_heat=>this%Q_supply_heat,&
				Q_supply_cool=>this%Q_supply_cool,&
				mdotcp_dhw=>this%DHW_tank%mdotcp_dhw,&
				V_storage_dhw=>this%DHW_tank%volume,&
				UA_storage_dhw=>this%DHW_tank%UA,&
				Q_heater_dhw=>this%DHW_tank%Q_heater,&
				T_main=>this%DHW_tank%T_main)
                ! reset totals for the substation
                Q_supply_heat=0d0
                Q_supply_cool=0d0
                mdotcp_dhw=0d0
                V_storage_dhw=0d0
                UA_storage_dhw=0d0
                Q_heater_dhw=0d0
                T_main=0d0
                ! calculate totals for the substation based on demand of the related building units
				num=0
				do j=1,size(this%building_units)
					do bu=1,size(buildings)
						do un=1,size(buildings(bu)%units)
							if(trim(buildings(bu)%name)//'('//integer_to_text(un)//')'==trim(this%building_units(j)%text)) then
							  num=num+1
							  Q_supply_heat=Q_supply_heat+buildings(bu)%units(un)%hvac%Q_heat_hyd+buildings(bu)%units(un)%hvac%Q_heat_AHU
							  Q_supply_cool=Q_supply_cool+buildings(bu)%units(un)%hvac%Q_cool_hyd+buildings(bu)%units(un)%hvac%Q_cool_AHU
							  mdotcp_dhw=mdotcp_dhw+buildings(bu)%units(un)%dhw%mdotcp_dhw
							  V_storage_dhw=V_storage_dhw+buildings(bu)%units(un)%dhw%V_storage 
							  UA_storage_dhw=UA_storage_dhw+buildings(bu)%units(un)%dhw%UA_storage 
							  Q_heater_dhw=Q_heater_dhw+buildings(bu)%units(un)%dhw%Q_heater
							  T_main=T_main+buildings(bu)%units(un)%dhw%T_main
							end if
						end do
					end do
				end do
				if (num>0) T_main=T_main/num			
                end associate	
                this%DHW_tank%T_loss=>T_unheated
	 end subroutine     
     subroutine substation_solve(this,hour)
	        class(substation)::this
	        integer::hour
	 end subroutine		
	 	 
	 ! energy centre (ancestor) dummy methods
	 subroutine energy_centre_load(this)
		class(energy_centre)::this
	 end subroutine
	 subroutine energy_centre_create_networks(this)
		class(energy_centre)::this
	 end subroutine
	 function energy_centre_num_elements(this) result(num)
		class(energy_centre)::this
		type(num_elements)::num
	 end function
     subroutine energy_centre_initialize(this)
	        class(energy_centre)::this
	 end subroutine     
     subroutine energy_centre_solve(this,hour,time)
	        class(energy_centre)::this
	        integer::hour
			double precision,dimension(2)::time
	 end subroutine		 
     subroutine energy_centre_update(this,hour,time)
	        class(energy_centre)::this
	        integer::hour
			double precision,dimension(2)::time
	 end subroutine		 
	 ! update incrementally an energy centre output during the integration inside the principal timestep of duration TIMESTEP
	 ! maintaining always the output value to the average value within the interval [0,time] 
	 subroutine energy_centre_update_output(this,arr,current_val,hour,time) 
		class(energy_centre)::this
		double precision,dimension(:),intent(inout)::arr
        double precision::current_val
		integer::hour
		double precision::time
		if (this%last_update>0d0) then
			arr(hour)=(arr(hour)*this%last_update+current_val*(time-this%last_update))/time
		else
			arr(hour)=current_val
		end if	
	 end subroutine	 

     ! create hydraulics nodes, pumps, heat exchangers, pipes for all objects comprising the thermal grid
	 ! performing initializations
     subroutine thermal_grid_create_networks(this)
	    class(thermal_grid)::this
	    integer::i
		type(num_elements)::num
		! for safe restart, deallocate all hydraulic nodes
		if (allocated(this%hydraulic_network%node)) deallocate(this%hydraulic_network%node)
        ! calculate number of elements of hydraulic network and thermal network
        num=num_elements(num_pipes=0,num_heat_exchangers=0,num_pumps=0,num_thermal_nodes=0)
	    do i=1,size(this%double_pipes)
		    num=num+this%double_pipes(i)%num_elements()
	    end do
	    do i=1,size(this%substations)
		    num=num+this%substations(i)%item%num_elements()
	    end do
	    do i=1,size(this%energy_centres)
		    num=num+this%energy_centres(i)%item%num_elements()
	    end do	
		! allocate elements of hydraulic network
        allocate(this%hydraulic_network%elements(num%hydraulic_channels_like()))
		! allocate thermal nodes of thermal network
        allocate(this%thermal_network%thermal_nodes(num%num_thermal_nodes))
		! allocate non-channel elements of thermal grid
        allocate(this%elements(num%non_hydraulic_channels_like()))
        this%num_tanks=num%num_tanks
        ! create connections
	    do i=1,size(this%double_connections)
		    call this%double_connections(i)%create_networks()
        end do
        ! create pipes, substations, energy_centres
	    do i=1,size(this%double_pipes)
		    call this%double_pipes(i)%create_networks()
	    end do
	    do i=1,size(this%substations)
		    call this%substations(i)%item%create_networks()
	    end do
	    do i=1,size(this%energy_centres)
		    call this%energy_centres(i)%item%create_networks()
	    end do	
	    ! initialize network (matrixes of incidence)
        call this%initialize()
	    ! initialization of substations and energy_centres
		do i=1,size(this%substations)
			call this%substations(i)%item%initialize()
		end do
		do i=1,size(this%energy_centres)
			call this%energy_centres(i)%item%initialize()
		end do
     end subroutine
	 
	 ! call the solve method of each substation
	 subroutine thermal_grid_solve_substations(this,hour)
		class(thermal_grid)::this
		integer::i,hour
	    do i=1,size(this%substations)
		    call this%substations(i)%item%solve(hour)
	    end do	 
	 end subroutine

	 ! call the solve method of each energy_centre
	 subroutine thermal_grid_solve_energy_centres(this,hour,time)
		class(thermal_grid)::this
		integer::i,hour
		double precision,dimension(2)::time
	    do i=1,size(this%energy_centres)
		    call this%energy_centres(i)%item%solve(hour,time)
	    end do	 
	 end subroutine

	 ! call the update method of each energy_centre
	 subroutine thermal_grid_update_energy_centres(this,hour,time)
		class(thermal_grid)::this
		integer::i,hour
		double precision,dimension(2)::time
	    do i=1,size(this%energy_centres)
		    call this%energy_centres(i)%item%update(hour,time)
	    end do	 
	 end subroutine
	      
     ! add a connection node to the hydraulic network and return its index
     function thermal_grid_add_node(this,name,P0) result(idx)
        class(thermal_grid),target::this
		character(*)::name
        type(point)::P0
        type(connection),dimension(:),allocatable::tmp_node
        integer::idx
        integer::i
        logical::found
		associate(h_n=>this%hydraulic_network)
        if (.not.allocated(h_n%node)) then
            idx=1
            allocate(h_n%node(1))
            h_n%node(1)%P0=P0
			h_n%node(1)%name=name
        else  
            ! test if node already exists
            found=.false.
            do i=1,size(h_n%node)
                if (h_n%node(i)%P0%equal(P0,0.1d0)) then
                    found=.true.
                    exit 
                end if
            end do    
            if (found) then
                idx=i
				h_n%node(idx)%name=name
            else    
                ! increase node array size and add a new node at the end
                allocate(tmp_node(size(h_n%node)))
                do i=1,size(h_n%node)
                    tmp_node(i)=h_n%node(i)
                end do    
                deallocate(h_n%node)
                allocate(h_n%node(size(tmp_node)+1))
                do i=1,size(h_n%node)-1
                    h_n%node(i)=tmp_node(i)
                end do    
                idx=size(h_n%node)
                h_n%node(idx)%P0=P0
				h_n%node(idx)%name=name
            end if    
        end if  
		end associate	
     end function

!=============================================================================================	 
!	hydraulic elements methods 
!=============================================================================================	 

 	 function init_pipe(name,node_in,node_out,D,L,K,eps,ins_delta,ins_lambda) result(new_pipe)
		double precision, intent(in)::D,L,K,eps,ins_delta,ins_lambda
		type(pipe)::new_pipe
		character(*)::name
		integer::node_in,node_out
		new_pipe%name=name
		new_pipe%node_in=node_in
		new_pipe%node_out=node_out
		! all inputs at this stage must be in SI units
		new_pipe%D=D
		new_pipe%L=L
		new_pipe%K=K
		new_pipe%eps=eps
		new_pipe%ins_delta=ins_delta
		new_pipe%ins_lambda=ins_lambda
		new_pipe%MC=L*(D**2)/4*pi*RHO_WATER*CP_WATER  ! thermal capacity (J/K)
		new_pipe%UA=pi*L*2*ins_lambda/log((D+2*ins_delta)/D) ! UA (W/K)
	 end function

	 subroutine pipe_initialize(this)
		class(pipe)::this
		this%MC=this%L*(this%D**2)/4*pi*RHO_WATER*CP_WATER  ! thermal capacity (J/K)
		this%UA=pi*this%L*2*this%ins_lambda/log((this%D+2*this%ins_delta)/this%D) ! UA (W/K)
	 end subroutine

	 subroutine heat_exchanger_initialize(this)
		class(heat_exchanger)::this	
		this%MC=this%Q_nom*1d-3*(0.1d-3*RHO_WATER*CP_WATER + 0.5*500)! thermal capacity (J/K)
		this%UA=0d0
	 end subroutine

	 subroutine thermal_grid_add_pipe(this,element) 
		class(thermal_grid),target::this
        type(pipe),target::element
 		associate(h_n=>this%hydraulic_network,t_n=>this%thermal_network)
		! increment networks counters
		h_n%num_pipes=h_n%num_pipes+1 
		t_n%num_thermal_nodes=t_n%num_thermal_nodes+1
		! set the fluid 
		element%fluid=>this%fluid
		! initialize the pipe
		call element%initialize()
		! add  element to the hydraulic network
		h_n%elements(h_n%num_pipes+h_n%num_pumps+h_n%num_heat_exchangers)%item=>element
		! create and add the associated thermal nodes to the thermal network
		t_n%thermal_nodes(t_n%num_thermal_nodes)=thermal_node(idx=t_n%num_thermal_nodes,MC=element%MC,UA=element%UA,T_env=T_ground)
		! save thermal node pointer
		if (associated(element%thermal_nodes)) then
			deallocate(element%thermal_nodes)
		end if	
        element%thermal_nodes=>t_n%thermal_nodes(t_n%num_thermal_nodes:t_n%num_thermal_nodes)
		end associate
	 end subroutine

	 subroutine thermal_grid_add_pump(this,element) 
		class(thermal_grid),target::this
        type(pump),target::element
 		associate(h_n=>this%hydraulic_network,t_n=>this%thermal_network)
		! increment networks counters
		h_n%num_pumps=h_n%num_pumps+1 
		t_n%num_thermal_nodes=t_n%num_thermal_nodes+1
		! set the fluid 
		element%fluid=>this%fluid
		! initialize
		call element%initialize()
		! add  element to the hydraulic network
		h_n%elements(h_n%num_pipes+h_n%num_pumps+h_n%num_heat_exchangers)%item=>element
		! create and add the associated thermal nodes to the thermal network
		t_n%thermal_nodes(t_n%num_thermal_nodes)=thermal_node(idx=t_n%num_thermal_nodes,MC=element%MC,UA=element%UA,T_env=T_unheated)
		! save thermal node pointer
		if (associated(element%thermal_nodes)) deallocate(element%thermal_nodes)
        element%thermal_nodes=>t_n%thermal_nodes(t_n%num_thermal_nodes:t_n%num_thermal_nodes)
		end associate
	 end subroutine
	 
	 subroutine thermal_grid_add_heat_exchanger(this,element) 
		class(thermal_grid),target::this
        type(heat_exchanger),target::element
 		associate(h_n=>this%hydraulic_network,t_n=>this%thermal_network)
		! increment networks counters
		h_n%num_heat_exchangers=h_n%num_heat_exchangers+1 
		t_n%num_thermal_nodes=t_n%num_thermal_nodes+1
		! set the fluid 
		element%fluid=>this%fluid
		! initialize
		call element%initialize()
		! add  element to the hydraulic network
		h_n%elements(h_n%num_pipes+h_n%num_pumps+h_n%num_heat_exchangers)%item=>element
		! create and add the associated thermal nodes to the thermal network
		t_n%thermal_nodes(t_n%num_thermal_nodes)=thermal_node(idx=t_n%num_thermal_nodes,MC=element%MC,UA=element%UA,T_env=T_unheated)
		! save thermal node pointer
		if (associated(element%thermal_nodes)) deallocate(element%thermal_nodes)
        element%thermal_nodes=>t_n%thermal_nodes(t_n%num_thermal_nodes:t_n%num_thermal_nodes)
		end associate
	 end subroutine

	 function init_borehole_field(name,node_in,node_out,DP_nom,m_nom,bhe) result(new_bhe)
		!integer,intent(in)::bhe_type
		!double precision,intent(in)::D_bhe,H,D_p,s_p,D,lambda_p,lambda_grout
		character(*)::name
		integer::node_in,node_out
		double precision,intent(in)::DP_nom,m_nom
		type(borehole_heat_exchanger)::bhe
		type(borehole_field),allocatable::new_bhe	
		allocate(new_bhe)
		new_bhe%name=name
		new_bhe%node_in=node_in
		new_bhe%node_out=node_out
		! all inputs at this stage must in SI units
		new_bhe%DP_nom=DP_nom
		new_bhe%m_nom=m_nom
		new_bhe%MC=bhe%N_bhe*bhe%H*((bhe%D_bhe)**2)/4*pi*RHO_WATER*CP_WATER ! thermal capacity (J/K)
		new_bhe%UA=0d0
		new_bhe%bhe=bhe
        !allocate(new_bhe%bhe%g_function)
	 end function

	 subroutine borehole_field_initialize(this)
		class(borehole_field)::this	
		this%MC=this%bhe%N_bhe*this%bhe%H*((this%bhe%D_bhe)**2)/4*pi*RHO_WATER*CP_WATER ! thermal capacity (J/K)
		this%UA=0d0
	 end subroutine

	 subroutine thermal_grid_add_borehole_field(this,element) 
		class(thermal_grid),target::this
        type(borehole_field),target::element
		associate(h_n=>this%hydraulic_network,t_n=>this%thermal_network)
		! increment networks counters
		!h_n%num_borehole_fields=h_n%num_borehole_fields+1 
		h_n%num_heat_exchangers=h_n%num_heat_exchangers+1 
		t_n%num_thermal_nodes=t_n%num_thermal_nodes+1
		! set the fluid 
		element%fluid=>this%fluid
		! initialize
		call element%initialize()
		! add  element to the hydraulic network
		h_n%elements(h_n%num_pipes+h_n%num_pumps+h_n%num_heat_exchangers)%item=>element
		! set the returned pointer
		! create and add  the associated thermal nodes to the thermal network
		t_n%thermal_nodes(t_n%num_thermal_nodes)=thermal_node(idx=t_n%num_thermal_nodes,MC=element%MC,UA=0d0,T_env=T_ground)
		! save thermal node pointer(s) 
		if (associated(element%thermal_nodes)) deallocate(element%thermal_nodes)
        element%thermal_nodes=>t_n%thermal_nodes(t_n%num_thermal_nodes:t_n%num_thermal_nodes)		
		end associate
	 end subroutine
	 
	 subroutine thermal_grid_add_buffer_tank(this,element)
		class(thermal_grid),target::this
        type(buffer_tank),target::element
        integer::i
		associate(h_n=>this%hydraulic_network,t_n=>this%thermal_network)
		! increment elements counter
		this%num_elements=this%num_elements+1
		! set the fluid 
		element%fluid=>this%fluid
		! add  element to the thermal grid
		this%elements(this%num_elements)%item=>element
		! create and add the associated thermal nodes to the thermal network, saving thermal node pointers
		if (associated(element%thermal_nodes)) deallocate(element%thermal_nodes)
        element%thermal_nodes=>t_n%thermal_nodes(t_n%num_thermal_nodes+1:t_n%num_thermal_nodes+element%num_thermal_nodes)
		do i=1,element%num_thermal_nodes
			t_n%num_thermal_nodes=t_n%num_thermal_nodes+1
			t_n%thermal_nodes(t_n%num_thermal_nodes)=thermal_node(idx=t_n%num_thermal_nodes,MC=element%MC/element%num_thermal_nodes,UA=element%UA/element%num_thermal_nodes,T_env=T_ground)
		end do	
		end associate
	 end subroutine
	 
	 ! solve the thermal problem using Euler's Explicit method
	 ! - variable time step, controlled by Courant < 1 and dT max=0.1 K  
	 subroutine thermal_grid_solve_Tq(this,hour,time)
        class(thermal_grid),target::this
        integer,intent(in)::hour ! the current hour
		double precision,intent(inout)::time ! seconds from beginning of current hourly time step
		double precision::dtime,dtime_c,dtime_T,dtime_c_min ! inner, variable, time step
		integer::next_hour,i,j,k,l,j_top,j_base,flag_T
		double precision::m0,mdot_ki,mdot_ij,mdot_down,mdot_up
        type(tank),pointer::p_test
		associate(h_n=>this%hydraulic_network,t_n=>this%thermal_network)
		! energy conservation equation for each thermal node:
		! (dT/dt)=1/MC*(-UA*(T-T_env)+Qin(T)+mdot*cp*(Tin-T))
		t_n%C=0d0 ! the Courant time step for each node = MC/sum(cp*mdot_in)	
		t_n%mdotcp=0d0 ! sum(cp*mdot_in)	
		! advection terms for Courant time step calculation
		do i=1,size(t_n%mass_flows)
			if (t_n%mass_flows(i)%mass_flow_rate>0) then
				! add heat capacity rate to node_out 
				t_n%mdotcp(t_n%mass_flows(i)%node_out)=t_n%mdotcp(t_n%mass_flows(i)%node_out)+t_n%mass_flows(i)%mass_flow_rate*this%fluid%cp
			end if	
        end do  
		! calculate Courant time steps
		! decide if node is domindated by advection or not
        dtime_c_min=0d0
		do i=1,size(t_n%thermal_nodes)
			if (t_n%mdotcp(i)>1.5d0*MIN_FLOW_RATE*this%fluid%cp) then
				t_n%solver_mode(i)=ADVECTIVE 				
				t_n%C(i)=min(0.5d0*TIMESTEP,t_n%thermal_nodes(i)%MC/t_n%mdotcp(i))
                dtime_c_min=max(dtime_c_min,t_n%C(i))
			else
				t_n%solver_mode(i)=CAPACITIVE_NON_ADVECTIVE ! NOTE all nodes must have MC>=MIN_MC, some can have UA=0				
				t_n%C(i)=TIMESTEP
			end if	
        end do	
		! decide if ADVECTIVE nodes must be treated as capacitive based on Courant number
		! set capacitive flag (Courant time < 60 seconds is low-capacity if 60 seconds > 20% maximum Courant time)
        ! at least one node must be capacitive
		dtime_c_min=min(0.2d0*dtime_c_min,60d0)
		k=0 ! counter of linear algebraic equations
		do i=1,size(t_n%thermal_nodes)
		  if (t_n%solver_mode(i)==ADVECTIVE) then
		    if (t_n%C(i)<dtime_c_min) then
				t_n%solver_mode(i)=NON_CAPACITIVE_ADVECTIVE
				k=k+1
				t_n%algebraic_idx(i)=k
				t_n%C(i)=dtime_c_min
            else
				t_n%solver_mode(i)=CAPACITIVE_ADVECTIVE
			end	if
		  end if	
        end do	
        dtime_c=minval(t_n%C)
		! energy balance at the nodes (besides advection terms)
		! NOTE thermal_node%T(hour) contains the solution at hour-1 at the beginning of the time step, 
		! and during integration is updated regularly to the value at the end of the last computed inner time step
		! NOTE: integration scheme
		! 1. calculate T in non-capacitive advective nodes (algebraic)
		! 2. calculate derivatives for capacitive advective and non-advective
		!    a) include advection terms for advective (capacitive or non-capacitive) only
		!    b) include generation and heat loss 
		! 3. set variable integration timestep
		! 4. integrate T for capacitive advective and non-advective 
		! use t_n%D (derivatives of T) to store temporarily the RHS of energy conservation equation
		t_n%D=0d0
		k=0 ! counter of linear algebraic equations
		t_n%A=0d0
		t_n%B=0d0
		do i=1,size(t_n%thermal_nodes)
		    if (t_n%solver_mode(i)==NON_CAPACITIVE_ADVECTIVE) then
				! 1. calculate temperatures in non-capacitive nodes
				k=k+1
				t_n%A(k,k)=t_n%thermal_nodes(i)%UA+t_n%mdotcp(i)
				t_n%B(k,1)=t_n%thermal_nodes(i)%UA*t_n%thermal_nodes(i)%T_env(hour)+t_n%thermal_nodes(i)%Qin(hour)
				! include incrementally advection terms 
				do j=1,size(t_n%mass_flows)
					! consider only positive advection terms entering in node i (node_out==i)
					if (t_n%mass_flows(j)%mass_flow_rate>0 .and. t_n%mass_flows(j)%node_out==i) then
						! distributing them on LHS (A) if node_in is non-capacitive, else on the RHS (B)
						if (t_n%solver_mode(t_n%mass_flows(j)%node_in)==NON_CAPACITIVE_ADVECTIVE) then
						! remove advection heat from node_in with temperature of node_in 
							t_n%A(k,t_n%algebraic_idx(t_n%mass_flows(j)%node_in))=-t_n%mass_flows(j)%mass_flow_rate*this%fluid%cp
						else
							t_n%B(k,1)=t_n%B(k,1)+t_n%mass_flows(j)%mass_flow_rate*this%fluid%cp*t_n%thermal_nodes(t_n%mass_flows(j)%node_in)%T(hour)
						end if	
					end if	
				end do  
			else if (t_n%solver_mode(i)==CAPACITIVE_ADVECTIVE .or. t_n%solver_mode(i)==CAPACITIVE_NON_ADVECTIVE) then
				! 2.b. include generation and heat loss for capacitive nodes
				t_n%D(i)=t_n%thermal_nodes(i)%UA*(t_n%thermal_nodes(i)%T_env(hour)-t_n%thermal_nodes(i)%T(hour))+t_n%thermal_nodes(i)%Qin(hour)
			end if	
        end do
        if (k>0) then
			call linsolve(t_n%A(1:k,1:k),t_n%B(1:k,1),t_n%T(1,1:k),k,flag_T)
			if (flag_T/=0) then 
				call simlog%output( 'Error: solution of thermal problem failed at hour '//integer_to_text(hour)//'.')
				stop
			end if	
			k=0 
			do i=1,size(t_n%thermal_nodes)
				if (t_n%solver_mode(i)==NON_CAPACITIVE_ADVECTIVE) then
					k=k+1
					t_n%thermal_nodes(i)%T(hour)=t_n%T(1,k) 
				end if
			end do		
		end if	
		! 2.a. include advection terms
		do i=1,size(t_n%mass_flows)
			if (t_n%mass_flows(i)%mass_flow_rate>0) then
				!if (t_n%solver_mode(t_n%mass_flows(i)%node_in)==CAPACITIVE_ADVECTIVE .and. t_n%solver_mode(t_n%mass_flows(i)%node_out)==CAPACITIVE_ADVECTIVE) then
				if ((t_n%solver_mode(t_n%mass_flows(i)%node_in)==CAPACITIVE_ADVECTIVE .or. t_n%solver_mode(t_n%mass_flows(i)%node_in)==NON_CAPACITIVE_ADVECTIVE).and.&
				    (t_n%solver_mode(t_n%mass_flows(i)%node_out)==CAPACITIVE_ADVECTIVE .or. t_n%solver_mode(t_n%mass_flows(i)%node_out)==NON_CAPACITIVE_ADVECTIVE)) then
					! remove advection heat from node_in with temperature of node_in
					t_n%D(t_n%mass_flows(i)%node_in)=t_n%D(t_n%mass_flows(i)%node_in)-t_n%mass_flows(i)%mass_flow_rate*this%fluid%cp*t_n%thermal_nodes(t_n%mass_flows(i)%node_in)%T(hour)
					! add advection heat to node_out with temperature of node_in 
					t_n%D(t_n%mass_flows(i)%node_out)=t_n%D(t_n%mass_flows(i)%node_out)+t_n%mass_flows(i)%mass_flow_rate*this%fluid%cp*t_n%thermal_nodes(t_n%mass_flows(i)%node_in)%T(hour)
				end if	
			end if	
        end do    
		! 2. calculate derivatives for capacitive (advective and non-advective) nodes, and set derivative zero for non-capacitive nodes
		do i=1,size(t_n%thermal_nodes)
			 if (t_n%solver_mode(i)==CAPACITIVE_ADVECTIVE .or. t_n%solver_mode(i)==CAPACITIVE_NON_ADVECTIVE) then
				 t_n%D(i)=t_n%D(i)/t_n%thermal_nodes(i)%MC
			 else if (t_n%solver_mode(i)==NON_CAPACITIVE_ADVECTIVE) then
				 t_n%D(i)=0d0 
		     end if
		end do		
		! 3. calculate the inner timestep 
		dtime_T=minval(TOL_T/abs(t_n%D), mask = t_n%D/=0d0)
	    dtime=min(dtime_T,TIMESTEP-time,dtime_c)
		! 4. calculate new temperatures in capacitive nodes 
		do i=1,size(t_n%thermal_nodes)
		    if (t_n%solver_mode(i)==CAPACITIVE_ADVECTIVE .or. t_n%solver_mode(i)==CAPACITIVE_NON_ADVECTIVE) then
				t_n%thermal_nodes(i)%T(hour)=t_n%thermal_nodes(i)%T(hour)+t_n%D(i)*dtime
			end if
		end do
		! advance time	
		time=time+dtime
		if (time==TIMESTEP) then
			! propagate capacitive node temperatures to non-capacitive ones (in case the next hour flow rates are zero)
			t_n%A=0d0
			t_n%B=0d0
			k=0
			do i=1,size(t_n%thermal_nodes)
				if (t_n%solver_mode(i)==NON_CAPACITIVE_ADVECTIVE) then
					k=k+1
					t_n%A(k,k)=t_n%thermal_nodes(i)%UA+t_n%mdotcp(i)
					t_n%B(k,1)=t_n%thermal_nodes(i)%UA*t_n%thermal_nodes(i)%T_env(hour)+t_n%thermal_nodes(i)%Qin(hour)
					! include incrementally advection terms 
					do j=1,size(t_n%mass_flows)
						! consider only positive advection terms entering in node i (node_out==i)
						if (t_n%mass_flows(j)%mass_flow_rate>0 .and. t_n%mass_flows(j)%node_out==i) then
							! distributing them on LHS (A) if node_in is low-capacitive, else on the RHS (B)
							if (t_n%solver_mode(t_n%mass_flows(j)%node_in)==NON_CAPACITIVE_ADVECTIVE) then
								! remove advection heat from node_in with temperature of node_in 
								t_n%A(k,t_n%algebraic_idx(t_n%mass_flows(j)%node_in))=-t_n%mass_flows(j)%mass_flow_rate*this%fluid%cp
							else
								t_n%B(k,1)=t_n%B(k,1)+t_n%mass_flows(j)%mass_flow_rate*this%fluid%cp*t_n%thermal_nodes(t_n%mass_flows(j)%node_in)%T(hour)
							end if	
						end if	
					end do    
				end if	
			end do
			if (k>0) then
				! recalculate algebraic temperatures for low-capacitive nodes
				call linsolve(t_n%A(1:k,1:k),t_n%B(1:k,1),t_n%T(1,1:k),k,flag_T)
				if (flag_T/=0) then 
					call simlog%output( 'Error: solution of thermal problem failed at hour '//integer_to_text(hour)//'.')
					stop
				end if	
				k=0 
				do i=1,size(t_n%thermal_nodes)
					if (t_n%solver_mode(i)==NON_CAPACITIVE_ADVECTIVE) then
						k=k+1
						t_n%thermal_nodes(i)%T(hour)=t_n%T(1,k)
					end if
				end do	
			end if
			! carry forward solution to the next hour if the end of the TIMESTEP is reached
			next_hour=hour+1
			if (next_hour>8760) next_hour=1
			do i=1,size(t_n%thermal_nodes)
				t_n%thermal_nodes(i)%T(next_hour)=t_n%thermal_nodes(i)%T(hour)
			end do				
		end if
 		end associate
	 end subroutine
	 	 
     ! solve the hydraulic problem using Todini's method
     subroutine thermal_grid_solve_Hm(this,hour)
        class(thermal_grid),target::this
        double precision,parameter:: tol=1d-3 ! set to lower value? e.g. 1d-5
        integer::hour
        integer::i
        integer::j
        integer::count_fix_height
        integer::count_fix_flow
        integer::node_under_analysis
        integer::error_flag_inv
        integer::flag_H
        integer::imposed_flows,active_pumps
        double precision::lambda
        double precision::error
		integer::k,l,j_top,j_base,flag_T
		double precision::m0,mdot_ki,mdot_ij,mdot_down,mdot_up
        type(tank),pointer::p_test
		type(inflow),dimension(:),allocatable::tmp_inflows
		
		associate(h_n=>this%hydraulic_network,t_n=>this%thermal_network)		
		
        count_fix_height=0
        do i=1,size(h_n%node)
            if(h_n%node(i)%fix_height(hour)==1) then
                count_fix_height=count_fix_height+1
            end if
        end do
                
        count_fix_flow=0               
 
        do i=1,size(h_n%elements)
            if(h_n%elements(i)%item%fix_flow(hour)==1) then
                count_fix_flow=count_fix_flow+1
            end if
        end do

                
        call alloc_array_2(h_n%H_fix,count_fix_height,1)
        call alloc_array_2(h_n%q_hyd,size(h_n%node)-count_fix_height,1)
        call alloc_array_2(h_n%A21_tot,size(h_n%node)-count_fix_height,size(h_n%elements))
        call alloc_array_2(h_n%A10_tot,size(h_n%elements),count_fix_height)
        call alloc_array_2(h_n%A21,size(h_n%node)-count_fix_height,size(h_n%elements)-count_fix_flow)
        call alloc_array_2(h_n%A10,size(h_n%elements)-count_fix_flow,count_fix_height)
                
        associate(&                
        H_fix=>h_n%H_fix,&
        q_hyd=>h_n%q_hyd,&
        A21_tot=>h_n%A21_tot,&
        A10_tot=>h_n%A10_tot,&
        A21=>h_n%A21,&
        A10=>h_n%A10)
                
        count_fix_height=0
        do i=1,size(h_n%node)
            if(h_n%node(i)%fix_height(hour)==1) then
                count_fix_height=count_fix_height+1
                H_fix(count_fix_height,1)=h_n%node(i)%H(hour)
            end if
        end do
                
        q_hyd=0d0
        node_under_analysis=0
        do i=1,size(h_n%elements)
            if(abs(h_n%elements(i)%item%mdot(hour))<MIN_FLOW_RATE) h_n%elements(i)%item%mdot(hour)=sign(MIN_FLOW_RATE,h_n%elements(i)%item%mdot(hour))
        end do        
                
        do i=1,size(h_n%node)-count_fix_height              
            node_under_analysis=node_under_analysis+1
            if(h_n%node(i)%fix_height(hour)==1) node_under_analysis=i+1
            do j=1,size(h_n%elements)
                if(h_n%elements(j)%item%fix_flow(hour)==1 .and. (h_n%elements(j)%item%node_in==node_under_analysis .or. h_n%elements(j)%item%node_out==node_under_analysis)) then
                    if(h_n%elements(j)%item%node_out==node_under_analysis) q_hyd(i,1)=q_hyd(i,1)+h_n%elements(j)%item%mdot(hour)
                    if(h_n%elements(j)%item%node_in==node_under_analysis) q_hyd(i,1)=q_hyd(i,1)-h_n%elements(j)%item%mdot(hour)
                end if
            end do			
        end do
    
        A10_tot=0d0
        count_fix_height=1
        do i=1,size(h_n%node)
            if(h_n%node(i)%fix_height(hour)==1) then  
                do j=1,size(h_n%elements)
                    if(h_n%elements(j)%item%node_in==i) then
                        A10_tot(j,count_fix_height)=1d0
                    else if(h_n%elements(j)%item%node_out==i) then
                        A10_tot(+j,count_fix_height)=-1d0
                    end if
                end do

                count_fix_height=count_fix_height+1
            end if
        end do
                
        j=1
        do i=1,size(h_n%elements)
            if(h_n%elements(i)%item%fix_flow(hour)/=1) then
                A10(j,:)=A10_tot(i,:)
                j=j+1 
            end if
        end do		
                
        j=1
        do i=1,size(h_n%node)
            if(h_n%node(i)%fix_height(hour)/=1) then
                A21_tot(j,:)=h_n%Aa_hyd(i,1:(size(h_n%Aa_hyd,2)))
                j=j+1
            end if
        end do
                
        j=1
		do i=1,size(h_n%elements)
            if(h_n%elements(i)%item%fix_flow(hour)/=1) then
                A21(:,j)=A21_tot(:,i)
                j=j+1
            end if
        end do
                
        call alloc_array_2 (h_n%E_hyd,size(A21,2),1)
        call alloc_array_2 (h_n%E_hyd_tot,size(A21_tot,2),1)
        call alloc_array_2 (h_n%dEdm_hyd,size(A21,2),1)
        call alloc_array_2 (h_n%dEdm_hyd_tot,size(A21_tot,2),1)
        call alloc_array_2 (h_n%Em_hyd,size(A21,2),1)
        call alloc_array_2 (h_n%Em_hyd_tot,size(A21_tot,2),1)
        call alloc_array_2 (h_n%m,size(A21,2),1)
        call alloc_array_2 (h_n%m_new,size(A21,2),1)
        call alloc_array_2 (h_n%m_new_tot,size(A21_tot,2),1)
        call alloc_array_2 (h_n%a_hyd,size(A21,2),1)
        call alloc_array_2 (h_n%a_hyd_tot,size(A21_tot,2),1)
        call alloc_array_2 (h_n%b_hyd,size(A21,1),1)
        call alloc_array_2 (h_n%A12,size(A21,2),size(A21,1))
        call alloc_array_2 (h_n%Jacob,size(A21,2),size(A21,2))
        call alloc_array_2 (h_n%Jacob_inv,size(A21,2),size(A21,2))
        call alloc_array_2 (h_n%A11,size(A21,2),size(A21,2))
        call alloc_array_2 (h_n%A_hyd_sys,size(A21,1),size(A21,1))
        call alloc_array_2 (h_n%H,1,size(A21,1))
                
        associate(&                
        E_hyd=>h_n%E_hyd,&
        E_hyd_tot=>h_n%E_hyd_tot,&
        dEdm_hyd=>h_n%dEdm_hyd,&
        dEdm_hyd_tot=>h_n%dEdm_hyd_tot,&
        Em_hyd=>h_n%Em_hyd,&
        Em_hyd_tot=>h_n%Em_hyd_tot,&
        m=>h_n%m,&
        m_new=>h_n%m_new,&
        a_hyd=>h_n%a_hyd,&
        a_hyd_tot=>h_n%a_hyd_tot,&
        b_hyd=>h_n%b_hyd,&
        A12=>h_n%A12,&
        Jacob=>h_n%Jacob,&
        Jacob_inv=>h_n%Jacob_inv,&
        A11=>h_n%A11,&
        A_hyd_sys=>h_n%A_hyd_sys,&
        H=>h_n%H)
                
        ! Hydraulic problem
        do i=1,size(h_n%elements)
			select type(element=>h_n%elements(i)%item)
			class is (pipe)
				lambda=1d0/4d0*(log10(element%eps/element%fluid%rho/(3.71d0*element%D)))**(-2d0)
				element%a_hyd=(8d0/(element%fluid%rho**2d0*9.81d0*pi**2d0*element%D**4d0)*(lambda*element%L/element%D+element%K))
			class is (heat_exchanger)	
				element%a_hyd=element%DP_nom/(element%fluid%rho*9.81d0*element%m_nom**2d0)
			class is (pump)
				! missing a_hyd value ?
			end select
            h_n%m_new_tot(i,1)=h_n%elements(i)%item%mdot(hour)
        end do
		            
        j=1

        do i=1,size(h_n%elements)
            if(h_n%elements(i)%item%fix_flow(hour)/=1) then
                m_new(j,:)=h_n%m_new_tot(i,:)
                j=j+1 
            end if
        end do
                
        m=m_new
        A12=transpose(A21)
        error=tol*10d0
                
        imposed_flows=0 ! number of elements imposing a flow constraints
		active_pumps=0 ! number of pumps that are running at fixed speed
        do i=1,size(h_n%elements)
            if(abs(h_n%elements(i)%item%mdot(hour))>MIN_FLOW_RATE .and. h_n%elements(i)%item%fix_flow(hour)==1) then
                imposed_flows=imposed_flows+1
            end if
			select type(element=>h_n%elements(i)%item)
			class is (pump)
				if (element%fix_flow(hour)==0) then
					active_pumps=active_pumps+1
				end if	
			end select
        end do		
              
        if(imposed_flows>0 .or. active_pumps>0) then
                
            do while (error>tol)
                    
                if(minval(abs(h_n%m_new_tot))<MIN_FLOW_RATE) then
                    do i=1,size(h_n%elements)
                        if(abs(h_n%elements(i)%item%mdot(hour))<MIN_FLOW_RATE) then
                            h_n%elements(i)%item%mdot(hour)=sign(MIN_FLOW_RATE,h_n%elements(i)%item%mdot(hour))
                            h_n%m_new_tot(i,1)=h_n%elements(i)%item%mdot(hour)
                        end if
                    end do

                end if
                    
                j=1
				do i=1,size(h_n%elements)
					select type(element=>h_n%elements(i)%item)
					class is (pipe)
						element%E_hyd=element%mdot(hour)*abs(element%mdot(hour))*element%a_hyd
						element%dEdm_hyd=2d0*abs(element%mdot(hour))*element%a_hyd
					class is (heat_exchanger)	
						element%E_hyd=element%mdot(hour)*abs(element%mdot(hour))*element%a_hyd
						element%dEdm_hyd=2d0*abs(element%mdot(hour))*element%a_hyd
					class is (pump)
						element%a_hyd=-(element%a1+element%a2/abs(element%mdot(hour))+element%a3/abs(element%mdot(hour))**2d0)
						element%dEdm_hyd=-(2d0*abs(element%mdot(hour))*element%a1+element%a2)
						element%E_hyd=abs(element%mdot(hour))*abs(element%mdot(hour))*element%a_hyd
					end select					
                    h_n%elements(i)%item%E_hyd=h_n%elements(i)%item%mdot(hour)*abs(h_n%elements(i)%item%mdot(hour))*h_n%elements(i)%item%a_hyd
                    h_n%elements(i)%item%dEdm_hyd=2d0*abs(h_n%elements(i)%item%mdot(hour))*h_n%elements(i)%item%a_hyd
                    h_n%elements(i)%item%Em_hyd=abs(h_n%elements(i)%item%mdot(hour))*h_n%elements(i)%item%a_hyd
                    if(abs(h_n%elements(i)%item%mdot(hour))<MIN_FLOW_RATE) h_n%elements(i)%item%Em_hyd=sign(MIN_FLOW_RATE,h_n%elements(i)%item%mdot(hour))*h_n%elements(i)%item%a_hyd
                    a_hyd_tot(i,1)=h_n%elements(i)%item%a_hyd
                    E_hyd_tot(i,1)=h_n%elements(i)%item%E_hyd
                    dEdm_hyd_tot(i,1)=h_n%elements(i)%item%dEdm_hyd
                    Em_hyd_tot(i,1)=h_n%elements(i)%item%Em_hyd
                            
                    if(h_n%elements(i)%item%fix_flow(hour)/=1) then
                        a_hyd(j,:)=a_hyd_tot(i,:)
                        E_hyd(j,:)=E_hyd_tot(i,:)
                        dEdm_hyd(j,:)=dEdm_hyd_tot(i,:)
                        Em_hyd(j,:)=Em_hyd_tot(i,:)
                        j=j+1 
                    end if   
                end do
                    
                Jacob=0d0
                do i=1,size(A21,2)
                    Jacob(i,i)=dEdm_hyd(i,1)
                end do
        
                A11=0d0
                do i=1,size(A21,2)
                    A11(i,i)=Em_hyd(i,1)
                end do
                    
                call inv(Jacob,Jacob_inv,size(Jacob,1),error_flag_inv)
                                            
                A_hyd_sys=matmul(matmul(A21,Jacob_inv),A12)
                b_hyd=matmul(A21,m)-q_hyd-matmul(matmul(A21,Jacob_inv),matmul(A11,m)+matmul(A10,H_fix))
                                                
                call linsolve(A_hyd_sys,b_hyd,H,size(b_hyd,1),flag_H)
                if(flag_H/=0) then
					call simlog%output( 'Error: solution of hydraulic problem failed at hour '//integer_to_text(hour)//'.')
					stop
                end if        
                m_new=m-matmul(Jacob_inv,matmul(A11,m)+matmul(A12,transpose(H))+matmul(A10,H_fix))
                                   
                j=1
                do i=1,size(h_n%elements)
                    if(h_n%elements(i)%item%fix_flow(hour)/=1) then
                        h_n%m_new_tot(i,:)=m_new(j,:)
                        j=j+1 
                    end if
                    h_n%elements(i)%item%mdot(hour)=h_n%m_new_tot(i,1)
                end do				
                    
                error=maxval(abs(m_new-m))
                    
                m=m_new                    
            end do
                    
        else
            H(1,:)=H_fix(1,1)
            h_n%m_new_tot=1d-8
			do i=1,size(h_n%elements)
                h_n%elements(i)%item%mdot(hour)=h_n%m_new_tot(i,1)
            end do
			
        end if
                    
        j=1
        do i=1,size(h_n%node)
            if(h_n%node(i)%fix_height(hour)/=1) then
                h_n%node(i)%H(hour)=H(1,j)
                j=j+1
            end if
        end do
					        
        ! clean memory
        H_fix=0d0
        q_hyd=0d0
        A21_tot=0d0
        A10_tot=0d0
        A21=0d0
        A10=0d0
        E_hyd=0d0
        E_hyd_tot=0d0
        dEdm_hyd=0d0
        dEdm_hyd_tot=0d0
        Em_hyd=0d0
        Em_hyd_tot=0d0
        m=0d0
        m_new=0d0
        a_hyd=0d0
        a_hyd_tot=0d0
        b_hyd=0d0
        A12=0d0
        Jacob=0d0
        Jacob_inv=0d0
        A11=0d0
        A_hyd_sys=0d0
        H=0d0
        
        end associate
        end associate
		
		! generate positive mass flow rates
	    l=0
		do k=1,size(h_n%node)
			m0=0
			do i=1,size(h_n%node(k)%channels)
				mdot_ki=h_n%elements(h_n%node(k)%channels(i))%item%mdot(hour)*h_n%node(k)%signs(i) ! positive if entering the node k
				if (mdot_ki>0) then	
					m0=m0+mdot_ki
                end if
            end do
			if (m0>0) then
				do i=1,size(h_n%node(k)%channels)
					mdot_ki=h_n%elements(h_n%node(k)%channels(i))%item%mdot(hour)*h_n%node(k)%signs(i) ! positive if entering the node k
					if (mdot_ki>0) then	
						h_n%node(k)%a(i)=mdot_ki/m0
						h_n%node(k)%b(i)=0
					else if (mdot_ki<0) then
						h_n%node(k)%a(i)=0
						h_n%node(k)%b(i)=-mdot_ki/m0
					else
						h_n%node(k)%a(i)=0
						h_n%node(k)%b(i)=0
					end if	
				end do  
			end if	
			! generate positive flow rates i->j associated to node k
			do i=1,size(h_n%node(k)%channels)
				do j=1,size(h_n%node(k)%channels)
					mdot_ij=m0*h_n%node(k)%a(i)*h_n%node(k)%b(j)
					if (i/=j) then
						l=l+1
						! NOTE the two lists thermal_nodes and elements are not aligned. 
						! Thus, the element index must be converted into a thermal_node index using element%thermal_nodes(1)%idx
						t_n%mass_flows(l)=mass_flow(node_in=h_n%elements(h_n%node(k)%channels(i))%item%thermal_nodes(1)%idx,&
													node_out=h_n%elements(h_n%node(k)%channels(j))%item%thermal_nodes(1)%idx,&
													mass_flow_rate=mdot_ij)
					end if
				end do
			end do
        end do
        ! generate positive mass flow rates for tank elements
		do k=1,size(this%elements)
			select type(item=>this%elements(k)%item)
            class is(tank)
                p_test=>item
				! alter the connection of already existing mass flow rates between the channel-like elements connected to the tank's hydraulic node
				do j=1,l
					! replace mass_flow's node_out with the thermal node at the bottom of the tank 
					! when mass flow is from element linked to tank's base to element linked to tank's top
					if (t_n%mass_flows(j)%node_in==item%element_to_base%thermal_nodes(1)%idx .and. &
					    t_n%mass_flows(j)%node_out == item%element_to_top%thermal_nodes(1)%idx) then
						t_n%mass_flows(j)%node_out=item%thermal_nodes(item%num_thermal_nodes)%idx ! the last thermal node of the tank is the base
						mdot_up=t_n%mass_flows(j)%mass_flow_rate	
						j_base=j
					end if	
					! replace mass_flow's node_out with the thermal node at the top of the tank 
					! when mass flow is from element linked to tank's top to element linked to tank's base
					if (t_n%mass_flows(j)%node_in==item%element_to_top%thermal_nodes(1)%idx .and. &
					    t_n%mass_flows(j)%node_out == item%element_to_base%thermal_nodes(1)%idx) then
						t_n%mass_flows(j)%node_out=item%thermal_nodes(1)%idx ! the first thermal node of the tank is the top
						mdot_down=t_n%mass_flows(j)%mass_flow_rate	
						j_top=j
					end if	
				end do
				! having calculated mdot_down and mdot_up, add the mass_flows dual of j_top and j_base 
				l=l+1
				! dual of j_top mass flow
				t_n%mass_flows(l)%mass_flow_rate=mdot_up
				t_n%mass_flows(l)%node_in=t_n%mass_flows(j_top)%node_out
				t_n%mass_flows(l)%node_out=t_n%mass_flows(j_top)%node_in
				! dual of j_base mass flow
				l=l+1
				t_n%mass_flows(l)%mass_flow_rate=mdot_down
				t_n%mass_flows(l)%node_in=t_n%mass_flows(j_base)%node_out
				t_n%mass_flows(l)%node_out=t_n%mass_flows(j_base)%node_in
				! internal mass flows between layers
				do i=1,item%num_thermal_nodes-1
					l=l+1
					t_n%mass_flows(l)=mass_flow(node_in=item%thermal_nodes(i)%idx,&
												node_out=item%thermal_nodes(i+1)%idx,&
												mass_flow_rate=mdot_down)
					l=l+1
					t_n%mass_flows(l)=mass_flow(node_in=item%thermal_nodes(i+1)%idx,&
												node_out=item%thermal_nodes(i)%idx,&
												mass_flow_rate=mdot_up)
				end do								
			end select
		end do
		
		! following positive mass_flow_rates, find thermal node associated to the node_out and save links in the corresponding element if heat_exchanger-like
		do l=1,size(t_n%mass_flows)
			!write(*,*) l,t_n%mass_flows(l)%node_out
			do i=1,size(h_n%elements)
                !write(*,*) h_n%elements(i)%item%name
				select type(item=>h_n%elements(i)%item)
				! NOTE with "class is" it looks for heat exchangers and all its extensions, e.g., boreholes
				class is(heat_exchanger) 
					!write(*,*) '->', i, item%thermal_nodes(1)%idx
					if (l==1) then
                        item%num_inflows=0
						if (allocated(item%inflows)) deallocate(item%inflows)
						allocate(item%inflows(10)) ! make enough space (10 flows in one node!)
                    end if   
					if (item%thermal_nodes(1)%idx==t_n%mass_flows(l)%node_out) then
						item%num_inflows=item%num_inflows+1
						item%inflows(item%num_inflows)%T=>t_n%thermal_nodes(t_n%mass_flows(l)%node_in)%T
						item%inflows(item%num_inflows)%mdot=t_n%mass_flows(l)%mass_flow_rate						
					end if
					! at the end, remove space not used in each element
					if (l==size(t_n%mass_flows)) then
						allocate(tmp_inflows(item%num_inflows))
						tmp_inflows=item%inflows(1:item%num_inflows)
						deallocate(item%inflows)
						allocate(item%inflows(item%num_inflows))
						item%inflows=tmp_inflows
						deallocate(tmp_inflows)
					end if
				end select	
			end do
		end do
		
		end associate
                
    end subroutine    
	
	subroutine thermal_grid_calculate_energy_demand(this,input_hours)
		class(thermal_grid)::this
		integer,optional::input_hours
		integer::t,num_hours ! hour time step
		double precision::max_timestep,min_timestep
		double precision,dimension(2)::time ! current time and last time within the hourly time step
		integer::k,i
        ! connect the buildings energy demand to the related substations
        do i=1,size(this%substations)
            call substation_initialize(this%substations(i)%item)
        end do
        ! calculate the energy performance of the network
		num_hours=8760
		if (present(input_hours)) num_hours=input_hours
		do t=1,num_hours ! one year
			! solve thermal problem in the substations (thermal energy demands and DHW storage)
			call this%solve_substations(t) 
			! solve the hydraulic problem in the network
			call this%solve_Hm(t)
			if (debugging) then
			 write (*,*) ' hydraulic network solved at hour ',t
			 do i=1,size(this%hydraulic_network%node)
				 write (*,'(a,i,f7.2,f7.2,f7.2,f7.2)') 'node',i,&
					 this%hydraulic_network%node(i)%P0%X,&
					 this%hydraulic_network%node(i)%P0%Y,&
					 this%hydraulic_network%node(i)%P0%Z,&
					 this%hydraulic_network%node(i)%H(t)
			 end do
			 do i=1,size(this%hydraulic_network%elements)
				 associate (t_g=>this%thermal_network)
				 select type (element=>this%hydraulic_network%elements(i)%item)
				 class is (pipe)    
					 write (*,'(a,i,i,i,f7.2,f14.2)') 'pipe',i,element%node_in,element%node_out,element%mdot(t),element%thermal_nodes(1)%T(t)
				 class is (pump)    
					 write (*,'(a,i,i,i,f7.2,f14.2)') 'pump',i,element%node_in,element%node_out,element%mdot(t),element%thermal_nodes(1)%T(t)
				 class is (heat_exchanger)    
					 write (*,'(a,i,i,i,f7.2,f14.2)') 'HX  ',i,element%node_in,element%node_out,element%mdot(t),element%thermal_nodes(1)%T(t)
				 end select
				 end associate
			 end do
			end if 
			time=0d0
			k=0
			max_timestep=0d0
			min_timestep=TIMESTEP
			do while (time(1)<TIMESTEP)
				! solve thermal problem in the energy centres 
				call this%solve_energy_centres(t,time)
				! solve the thermal problem in the network
				time(2)=time(1)
				call this%solve_Tq(t,time(1))
				k=k+1
				max_timestep=max(max_timestep,time(1)-time(2))
				min_timestep=min(min_timestep,time(1)-time(2))
				! incremental update of TIMESTEP-averaged outputs
				call this%update_energy_centres(t,time)
			end do    
			if (debugging) then
			 write (*,*) ' thermal network solved at hour ',t
			 do i=1,size(this%hydraulic_network%elements)
				 associate (t_g=>this%thermal_network)
				 select type (element=>this%hydraulic_network%elements(i)%item)
				 class is (pipe)    
					 write (*,'(a,i,i,i,f7.2,f14.2)') 'pipe',i,element%node_in,element%node_out,element%mdot(t),element%thermal_nodes(1)%T(t)
				 class is (pump)    
					 write (*,'(a,i,i,i,f7.2,f14.2)') 'pump',i,element%node_in,element%node_out,element%mdot(t),element%thermal_nodes(1)%T(t)
				 class is (heat_exchanger)    
					 write (*,'(a,i,i,i,f7.2,f14.2)') 'HX  ',i,element%node_in,element%node_out,element%mdot(t),element%thermal_nodes(1)%T(t)
				 end select
				 end associate
			 end do
			 do i=1,size(this%elements)
				 associate (t_g=>this%thermal_network)
				 select type (element=>this%elements(i)%item)
				 class is (tank)    
					 write (*,'(a,i,f14.2,f14.2,f14.2)') 'tank',i,element%thermal_nodes(1)%T(t),&
						 element%thermal_nodes(2)%T(t),element%thermal_nodes(3)%T(t)
				 end select
				 end associate
			 end do
			 write (*,*) ' hour, timesteps count, min, max',t,k,min_timestep,max_timestep  
            end if 
            if (mod(t, 168) == 0) then 
                call simlog%output('.','a')
            end if    
		end do
		call simlog%output('')
	end subroutine
    
end module