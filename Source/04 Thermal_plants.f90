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
! The module contains classes and methods related to the main thermal plants components, including:
! - algebraic components (A), 
! - dynamical components (D),
! - controllers (C).
! The requirement of algebraic components is to calculate average output values within the main time step of the simulation (one hour).
! They can have internal states, which are updated at the end of the main time step, but their calculation is carried out internally, without
! interacting with the dynamical solver (e.g. DHW storage and its state of charge). 
! The load of an algebraic generator is pre-defined. Virtual auxiliary are used to compensate insufficient capacity.     
! Dynamical components are solved simultaneously with the thermal grid, thus interacting with the dynamical solver (e.g., ground-coupled heat exchanger)
! The load of a dynamical generator is calculated based on its controller (e.g. temperature setpoint). Virtual auxiliaries are useless in this case.   
! Algebraic components require two main methods: 
! - initialize() 		  : set internal parameters, allocate memory, pre-calculations (e.g. automatic sizing)
! - calculate(hour)  	  : calculate energy flows (and internal states if present) at each main (hourly) time step  
! Dynamical components require three main methods: 
! - initialize() 		  : set internal parameters, allocate memory, pre-calculations (e.g. automatic sizing)
! - calculate(hour,time)  : calculate energy flows at each inner (sub-hourly) time step  
! - update(hour,time)	  : calculate internal states and outputs at the end of the main (hourly) time step
! NOTE One component can be defined with both Algebraic and Dynamical features (the calculate method must differentiate the treatment of the load and the virtual auxiliary)
    
module thermal_plants_module  
	use constants_module
	use maths
    use datalog
	use json_module
    use utility
    implicit none
    
 	! plant component feature
    enum, bind(c)
        enumerator::ALGEBRAIC=1,DYNAMICAL,CONTROLLER
    end enum 

    ! borehole heat exchanger settings
    enum, bind(c)
        enumerator::MANUAL_INPUT=1,AUTO_DESIGN,REBUILD_ALWAYS,BUILD_IF_EMPTY
    end enum 
    
	! operation modes
    enum, bind(c)
        enumerator::OFF=1,STANDBY,HEATING,COOLING,COOLING_PLUS_HEAT_RECOVERY,COOLING_ON,HEATING_ON,INACTIVE
    end enum 
		
	! simple incompressible fluid with constant properties
    type incompressible_fluid
        double precision::cp        ! specific heat, J/(kg,K)
        double precision::rho       ! density, kg/m3
	contains	
		procedure::to_json=>incompressible_fluid_to_json
		procedure::from_json=>incompressible_fluid_from_json
    end type

	! air-source heat pump / chiller (A or D)
	type as_heat_pump_chiller
        integer::component_feature=DYNAMICAL			! Dynamical by default, it can be turned into ALGEBRAIC (i.e., for use inside subsations)
		! nominal data
		double precision::Q_heat_nom=100d3				! heating capacity at nominal condition, W
		double precision::COP_nom=3.8d0					! COP at nominal condition
		double precision::T_hwo_nom=50d0				! hot water leaving temperature at nominal condition, °C 
		double precision::T_air_h_nom=7d0				! air temperature at nominal condition in heating, °C
		double precision::Q_cool_nom=100d3				! cooling capacity at nominal condition, W (>0)
		double precision::EER_nom=3.8d0					! EER at nominal condition
		double precision::T_cwo_nom=7d0					! chilled water leaving temperature at nominal condition, °C 
		double precision::T_air_c_nom=32d0				! air temperature at nominal condition in cooling, °C
		double precision::W_standby=10					! electricity consumption in standby, W
		type(incompressible_fluid)::fluid=>incompressible_fluid(cp=CP_WATER,rho=RHO_WATER)	! fluid on demand side (water with or without glycol)
		! data for capacity and COP/EER variations (based on the predefined model)
        integer::model_type                             ! 1) Carnot,PLR quadratic,temperature-based frost efficiency (F=Carnot efficiency,Q=PLR*Qmax, Qmax = Q at full load at same operating conditions)
														!    The parameters of the model are 16, 1-10 for heating and 11-16 for cooling:
														!    (heating: 1-10)
														!	 DT_e (water-refrigerant), DT_c (water-refrigerant), p (exponent for F/Fnom), r (exponent for DTlift/DTlift,nom), 
                                                        !    a0, a1 (coefficients for W/Wmax=a0+a1*PLR+(1-a0-a1)*PLR^2)
                                                        !    b,T0,alpha_COP,alpha_Q (coefficients for frosting degradation)
														!    (cooling: 11-16)
														!	 DT_e (water-refrigerant), DT_c (water-refrigerant), p (exponent for G/Gnom), r (exponent for DTlift/DTlift,nom), 
                                                        !    a0, a1 (coefficients for W/Wmax=a0+a1*PLR+(1-a0-a1)*PLR^2)
        double precision,dimension(:),allocatable::model_par ! array of the model parameters
		! input data
		double precision::T_wo_set=50d0					! set-point for the water leaving temperature 
		double precision::T_air=10d0					! air temperature at work condition
		double precision::DT_w=5d0						! absolute inlet/outlet temperature difference of the water 
		integer::operation_mode=OFF						! operation mode (OFF, STANDBY, HEATING or COOLING)
		double precision::T_wi=0d0						! temperature of the inlet water
		! output data	
		double precision::mdot_w=0d0					! the water mass flow rate such that nominal capacity is delivered at the imposed DT_w 
		double precision::T_wo=0d0						! the water leaving temperature at work condition 		
		double precision::Qdot_heat=0d0,Qdot_cool=0d0,Wdot_hpc=0d0,PLR_now=0d0  ! instantaneous heating (W), cooling (W), power input (W) and PLR (-)
		double precision,dimension(8760)::Q_heat=0d0	! average heating capacity in the timestep, W (output if Dynamical, input if Algebraic)
		double precision,dimension(8760)::Q_cool=0d0	! average cooling capacity in the timestep, W (output if Dynamical, input if Algebraic)
		double precision,dimension(8760)::Q_aux_h=0d0	! average heating capacity of virtual auxiliary heater, W
		double precision,dimension(8760)::Q_aux_c=0d0	! average cooling capacity of virtual axuiliary cooler, W
		double precision,dimension(8760)::PLR=0d0		! partial load ratio at work condition
		double precision,dimension(8760)::W_hpc=0d0 	! average electricity consumption of compressor and fan, W
	contains
		procedure::initialize=>as_heat_pump_chiller_initialize
		procedure::calculate=>as_heat_pump_chiller_calculate
		procedure::update=>as_heat_pump_chiller_update
		procedure::to_json=>as_heat_pump_chiller_to_json
		procedure::from_json=>as_heat_pump_chiller_from_json
		procedure::get_par=>as_heat_pump_chiller_get_par
    end type

	! water-water heat pump / chiller (A or D)
	type ww_heat_pump_chiller
        integer::component_feature=DYNAMICAL			! Dynamical by default, it can be turned into ALGEBRAIC (i.e., for use inside subsations)
		! nominal data
		double precision::Q_heat_nom=100d3				! heating capacity at nominal condition, W
		double precision::COP_nom=3.8d0					! COP at nominal condition
		double precision::T_hwo_nom=50d0				! hot water leaving temperature at nominal condition (heating mode), °C 
		double precision::T_cwo_h_nom=0d0				! chilled water leaving temperature at nominal condition in heating mode, °C
		double precision::Q_cool_nom=100d3				! cooling capacity at nominal condition, W (>0)
		double precision::EER_nom=3.8d0					! EER at nominal condition
		double precision::T_cwo_nom=0d0				    ! chilled water leaving temperature at nominal condition (cooling mode), °C
		double precision::T_hwo_c_nom=30d0				! hot water leaving temperature at nominal condition in cooling mode, °C 
		double precision::W_standby=10					! electricity consumption in standby, W
		type(incompressible_fluid)::fluid=>incompressible_fluid(cp=CP_WATER,rho=RHO_WATER)	! fluid on demand side (water with or without glycol) 
		! data for capacity and COP/EER variations (based on the predefined model)
        integer::model_type                             ! 1) Carnot,PLR quadratic (F=Carnot efficiency,Q=PLR*Qmax, Qmax = Q at full load at same operating conditions)
														!    The parameters of the model are 12, 1-6 for heating and 7-12 for cooling:
														!    (heating: 1-6)
														!	 DT_e (water-refrigerant), DT_c (water-refrigerant), p (exponent for F/Fnom), r (exponent for DTlift/DTlift,nom), 
                                                        !    a0, a1 (coefficients for W/Wmax=a0+a1*PLR+(1-a0-a1)*PLR^2)
                                                		!    (cooling: 7-12)
														!	 DT_e (water-refrigerant), DT_c (water-refrigerant), p (exponent for G/Gnom), r (exponent for DTlift/DTlift,nom), 
                                                        !    a0, a1 (coefficients for W/Wmax=a0+a1*PLR+(1-a0-a1)*PLR^2)
        double precision,dimension(:),allocatable::model_par ! array of the model parameters
		! input data
		double precision::T_wo_set=50d0					! set-point for the water leaving temperature 
		double precision::T_source_sink=10d0            ! the temperature of the water leaving the source (for heating) or the sink (for cooling)
		double precision::DT_w=5d0						! absolute inlet/outlet temperature difference of the water 
		integer::operation_mode=OFF						! operation mode (OFF, STANDBY, HEATING or COOLING)
		double precision::T_wi=0d0						! temperature of the inlet water
		! output data	
		double precision::mdot_w=0d0					! the water mass flow rate such that nominal capacity is delivered at the imposed DT_w 
		double precision::T_wo=0d0						! the water leaving temperature at work condition 		
		double precision::Qdot_heat=0d0,Qdot_cool=0d0,Wdot_hpc=0d0,PLR_now=0d0  ! instantaneous heating (W), cooling (W), power input (W) and PLR (-)
		double precision,dimension(8760)::Q_heat=0d0	! average heating capacity in the timestep, W (output if Dynamical, input if Algebraic)
		double precision,dimension(8760)::Q_cool=0d0	! average cooling capacity in the timestep, W (output if Dynamical, input if Algebraic)
		double precision,dimension(8760)::Q_aux_h=0d0	! average heating capacity of virtual auxiliary heater, W
		double precision,dimension(8760)::Q_aux_c=0d0	! average cooling capacity of virtual axuiliary cooler, W
		double precision,dimension(8760)::PLR=0d0		! partial load ratio at work condition
		double precision,dimension(8760)::W_hpc=0d0 	! average electricity consumption of compressor and fan, W
	contains
		procedure::initialize=>ww_heat_pump_chiller_initialize
		procedure::calculate=>ww_heat_pump_chiller_calculate
		procedure::update=>ww_heat_pump_chiller_update
		procedure::to_json=>ww_heat_pump_chiller_to_json
		procedure::from_json=>ww_heat_pump_chiller_from_json
		procedure::get_par=>ww_heat_pump_chiller_get_par
    end type
    
 	! water-water heat pump (A or D)
	type ww_heat_pump
        integer::component_feature=ALGEBRAIC              ! ALGEBRAIC by default, can be turned into a DYNAMICAL model
		! nominal data
		double precision::Q_heat_nom=100d3				  ! heating capacity at nominal condition, W
		double precision::COP_nom=4.9d0					  ! COP at nominal condition
		double precision::T_hwo_nom=35d0				  ! hot water leaving temperature at nominal condition, °C 
		double precision::T_cwo_nom=0d0				      ! chilled water leaving temperature at nominal condition, °C
		double precision::W_standby=10					  ! electricity consumption in standby, W        
		type(incompressible_fluid)::fluid=>incompressible_fluid(cp=CP_WATER,rho=RHO_WATER)	! fluid on demand side (water with or without glycol)
		! data for capacity and COP variations (based on the predefined model)
        integer::model_type                               ! 1) Carnot,PLR quadratic (F=Carnot efficiency,Q=PLR*Qmax, Qmax = Q at full load at same operating conditions)
														  !   The parameters of the model are 6:
														  !	  DT_e (water-refrigerant), DT_c (water-refrigerant), p (exponent for F/Fnom), r (exponent for DTlift/DTlift,nom), 
                                                          !   a0, a1 (coefficients for W/Wmax=a0+a1*PLR+(1-a0-a1)*PLR^2)
        double precision,dimension(:),allocatable::model_par ! array of model parameters
		! input data
		double precision::T_hwo_set=35d0				  ! hot water leaving temperature at work condition (set-point)
		double precision::T_cwo=0d0					      ! chilled water leaving temperature at work condition (assigned)
		double precision::DT_hw=5d0						  ! absolute inlet/outlet temperature difference of the water 
		integer::operation_mode=OFF						  ! operation mode (OFF, STANDBY, HEATING or COOLING)
		double precision::T_hwi=0d0						  ! temperature of hot water at the inlet, °C
		! output data	
		double precision::mdot_hw=0d0                     ! the hot water mass flow rate such that nominal capacity is delivered at the imposed DT_w 
		double precision::T_hwo=0d0                       ! the hot water leaving temperature at work condition 		
		double precision::Qdot_cond=0d0,Qdot_evap=0d0,Wdot_comp=0d0,PLR_now=0d0  ! instantaneous heating (W), cooling (W), power input (W) and PLR (-)
		double precision,dimension(8760)::Q_cond=0d0	  ! average heating capacity at condenser, W (output if Dynamical, input if Algebraic)
		double precision,dimension(8760)::Q_aux_h=0d0	  ! average heating capacity of auxiliary heater, W
		double precision,dimension(8760)::PLR=0d0		  ! partial load ratio at work condition
		double precision,dimension(8760)::Q_evap=0d0  	  ! average cooling capacity, W
		double precision,dimension(8760)::W_comp=0d0 	  ! average electricity consumption of the compressor, W
	contains
		procedure::initialize=>ww_heat_pump_initialize
		procedure::calculate=>ww_heat_pump_calculate
		procedure::update=>ww_heat_pump_update
		procedure::to_json=>ww_heat_pump_to_json
		procedure::from_json=>ww_heat_pump_from_json
        procedure::get_par=>ww_heat_pump_get_par
	end type
   
	! the g-function class
	type g_function
		double precision::t_s=0d0 ! characteristic time, s
		double precision,dimension(:),allocatable::tau  ! non-dimensional time (t/t_s)
		double precision,dimension(:),allocatable::gval ! g-function value
		type(linear_interp_1d),allocatable::interp 		! interpolation class 
	contains
		procedure::load=>g_function_load ! load data in the the interpolation class
		procedure::get=>g_function_get ! get the g-function value at a given time (seconds)
	end type	
    
	! borehole heat exchanger field (D)
    type borehole_heat_exchanger
        integer::component_feature=DYNAMICAL
        ! borehole data 
        integer::bhe_type 								! 1 for single U-pipe, 2 for double U-pipe
        double precision::D_bhe 						! borehole diameter [m]
        double precision::D_p 							! pipe outer diameter [m]
        double precision::s_p 							! pipe thickness [m]
        double precision::D 							! distance between the legs of a pipe [m]
        double precision::lambda_p 						! pipe thermal conductivity [W/(m K)]
        double precision::lambda_grout 					! grout thermal conductivity [W/(m K)]
        double precision,dimension(:,:),allocatable::bhe_XY_coordinates ! X,Y coordinates of BHEs
        ! fluid data
        double precision::rho_f=1000d0					! fluid density [kg/(m3)]
        double precision::cp_f=4186d0					! fluid specific heat [J/(kg K)]
        ! ground data
        double precision::lambda 						! ground thermal conductivity [W/(m K)]
        double precision::C 							! ground volumetric thermal capacity [J/(m3 K]
        double precision::Tg 							! undisturbed ground temperature [°C]
 		! input design data (only for automatic sizing option)
		double precision::B_des=7d0 					! design bhe to bhe distance [m]
        double precision::H_des=200d0 					! design bhe depth [m]
        double precision::Q_ghd=44211d0 				! heating design peak load, positive [W]
        double precision::Q_gcd=-104242d0				! cooling design peak load, negative [W]
        double precision::Q_ghm=15960d0 				! heating design month load, positive [W]
        double precision::Q_gcm=-30126d0 				! cooling design month load, negative [W]
        double precision::Q_ga=-3626d0 					! annual design load, balance between (positive) heating and (negative) cooling [W]
        double precision::T_fih_design=2d0 				! fluid temperature entering the BHE in heating operation (design) [°C]
        double precision::T_foh_design=6d0 				! fluid temperature exiting the BHE in heating operation (design) [°C]
        double precision::T_fic_design=32d0 			! fluid temperature entering the BHE in cooling operation (design) [°C]
        double precision::T_foc_design=28d0 			! fluid temperature exiting the BHE in cooling operation (design) [°C]
        ! output design data (only for automatic sizing option)
        double precision::L_new=0d0						! total length, m
        character(100)::mode=''                         ! 'heating' or 'cooling' is the predominant mode
        integer::N_bhe=0								! total number of boreholes (design or input)
        integer::Nx=0									! number of boreholes along X direction
        integer::Ny=0									! number of boreholes along X direction
        integer::N_odd=0								! number of boreholes along X direction
        double precision::H=0d0							! length of boreholes (m)
        double precision::R_bhe_cond=0d0				! conductive thermal resistance of the borehole (grout and pipe)
        ! g-function data
		integer::layout_mode							! MANUAL_INPUT or AUTO_DESIGN
		type(g_function)::g_function    				! the g-function object with data and methods for interpolation
		! output data	
        double precision::Qdot_ground=0d0 				! heat load to the ground, instantaneous [W]
        double precision,dimension(8760)::Q_ground=0d0 	! heat load to the ground, TIMESTEP-averaged [W]
        double precision,dimension(:),allocatable::DQ_g ! variation w.r.t. previous time step of heat transfer rate per unit length to the ground  [W/m]
        double precision,dimension(8760)::DT_g=0d0 		! variation w.r.t. the undisturbed ground temperature at borehole boundary [K]
		double precision::DT_g_next=0d0    				! the DT_g at the end of next hour due to previous heat interactions with the ground [K]
	contains	
		procedure::initialize=>borehole_heat_exchanger_initialize  				! set the main characteristics, coefficients and design conditions
		procedure::Reynolds=>borehole_heat_exchanger_Reynolds ! Reynolds number inside borehole pipes
		procedure::convective_resistance=>borehole_heat_exchanger_convective_resistance ! convective resistance (m*K/W) inside borehole pipes
		procedure::design=>borehole_heat_exchanger_design  						! automatically design the bhe field (compact rectangular) 
		procedure::create_g_function=>borehole_heat_exchanger_create_g_function ! calculate the g-function and store its values
		procedure::Qex=>borehole_heat_exchanger_Qex ! calculate the heat transfer from the borehole to the thermal medium
		procedure::update_DT_ground=>borehole_heat_exchanger_update_DT_ground ! update the temperature variation of the ground at the borehole perimeter
		procedure::update=>borehole_heat_exchanger_update ! update ground heat exchange and borehole temperature
		procedure::to_json=>borehole_heat_exchanger_to_json
		procedure::from_json=>borehole_heat_exchanger_from_json
    end type

	! domestic hot water storage (A)
	type dhw_storage
        integer::component_feature=ALGEBRAIC
		! inputs
		double precision::volume=0d0 						  ! water volume, m3
		double precision::UA=0d0 							  ! heat loss coefficient, W/K
		double precision::Q_heater=0d0 						  ! heating capacity of the DHW heat exchanger, W
		double precision,dimension(8760)::mdotcp_dhw=0d0	  ! the hourly demand of hot water at 40°C
		double precision,dimension(8760)::T_main=0d0 		  ! the temperature of water main
		double precision,dimension(:),pointer::T_loss=>null() ! the temperature of the environment used for the calculation of heat losses 
		! outputs
		double precision,dimension(8760)::T_storage=50d0 	  ! the hourly state of charge of the dhw tank
		double precision,dimension(8760)::Q_supply=0d0  	  ! the heating supplied to the dhw on average in the time step (W)
		! control parameters
		double precision::f_dhw_on=0d0                   	  ! the charging state (on/off) updated at the end of each main time step
		double precision::T_lower_dhw=45d0					  ! lowest temperature in the dhw tank (heater turns on)
		double precision::T_upper_dhw=60d0		        	  ! highest temperature in the dhw tank (heater turns off)
	contains
		procedure::initialize=>dhw_storage_initialize
		procedure::calculate=>dhw_storage_calculate
		procedure::to_json=>dhw_storage_to_json
		procedure::from_json=>dhw_storage_from_json
	end type

	! buffer controller data and methods (C)
	type buffer_controller 	
        integer::component_feature=CONTROLLER
		double precision,dimension(8760)::T_storage_max=18d0	! central node storage temperature above which cooling is activated, °C 
		double precision,dimension(8760)::T_storage_min=9d0 	! central node storage temperature below which heating is activated, °C
		! output
		integer::mode=INACTIVE									! signal communicating the operation status to the external heating/cooling generator: COOLING_OFF, COOLING_ON, HEATING_OFF, HEATING_ON
	contains
		procedure::update=>buffer_controller_update
		procedure::to_json=>buffer_controller_to_json
		procedure::from_json=>buffer_controller_from_json
	end	type
	
contains	

	subroutine incompressible_fluid_to_json(this,json,obj,obj_name)
		class(incompressible_fluid)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(inout)::obj
		type(json_value),pointer::obj_fluid
		character(*)::obj_name
		call json%create_object(obj_fluid,obj_name)
		call json%add(obj_fluid,'rho',this%rho)
		call json%add(obj_fluid,'cp',this%cp)
		call json%add(obj,obj_fluid)
    end subroutine		

	subroutine incompressible_fluid_from_json(this,json,obj,path)
		class(incompressible_fluid)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(in)::obj
		character(*)::path
		call json%get(obj,path//'.rho',this%rho)
		call json%get(obj,path//'.cp',this%cp)
    end subroutine		

	subroutine buffer_controller_update(this,hour,T_tank)
		class(buffer_controller)::this
		integer::hour
		double precision::T_tank
		if (T_tank<=this%T_storage_min(hour) .and. this%mode/=HEATING_ON) then
			this%mode=HEATING_ON
        else if (T_tank>this%T_storage_min(hour)+2d0 .and.this%mode==HEATING_ON ) then ! 2 K hysteresis
			this%mode=INACTIVE
        else if (T_tank>=this%T_storage_max(hour) .and. this%mode/=COOLING_ON)  then
			this%mode=COOLING_ON
        else if (T_tank<this%T_storage_max(hour)-2d0 .and.this%mode==COOLING_ON ) then ! 2 K hysteresis
			this%mode=INACTIVE
        end if    
    end subroutine

	! core subroutines for heat pump & chiller models	
    
    subroutine heat_pump_model_1_heating (T_e,T_c,T_e_nom,T_c_nom,Q_heat_nom,COP_nom,p,r,a0,a1,Q_load,COP,Qfull)
        double precision,intent(in)::Q_load,T_e,T_c,T_e_nom,T_c_nom,COP_nom,Q_heat_nom,r,p,a0,a1
        double precision,intent(out)::COP
        double precision,intent(out)::Qfull
        double precision::PLF,PLR,COPfull,T_e0
        ! adjust Te so to respect the typical compressor operating envelope limit (Tc-Te>15 K)
        T_e0=min(T_e,T_c-15d0)
		COPfull = COP_nom * (Carnot_heating(T_e0,T_c)/Carnot_heating(T_e_nom,T_c_nom))**p
		Qfull = Q_heat_nom/COP_nom * COPfull * ((T_c-T_e0)/(T_c_nom-T_e_nom))**r
        PLR = max(0.15d0,min(1d0,Q_load/Qfull))
		PLF = PartLoadFactor(PLR,a0,a1)
		COP = PLF * COPfull 
    end subroutine

    subroutine heat_pump_model_1_cooling (T_e,T_c,T_e_nom,T_c_nom,Q_cool_nom,EER_nom,p,r,a0,a1,Q_load,EER,Qfull)
        double precision,intent(in)::Q_load,T_e,T_c,T_e_nom,T_c_nom,EER_nom,Q_cool_nom,r,p,a0,a1
        double precision,intent(out)::EER
        double precision,intent(out)::Qfull
        double precision::PLF,PLR,EERfull,T_c0
        ! adjust Tc so to respect the typical compressor operating envelope limit (Tc-Te>15 K)
        T_c0=max(T_c,T_e+15d0)
		EERfull = EER_nom * (Carnot_cooling(T_e,T_c0)/Carnot_cooling(T_e_nom,T_c_nom))**p
		Qfull = Q_cool_nom/EER_nom * EERfull * ((T_c0-T_e)/(T_c_nom-T_e_nom))**r
        PLR = max(0.15d0,min(1d0,Q_load/Qfull))
		PLF = PartLoadFactor(PLR,a0,a1)
		EER = PLF *EERfull 
    end subroutine

    function Carnot_heating(T_e,T_c) result (F)
		double precision::T_e,T_c,T_eK,T_cK,Delta_T,F
		T_eK = T_e + 273.15d0
		T_cK = T_c + 273.15d0
		Delta_T = max(3d0,T_cK - T_eK)
		F = T_cK / Delta_T
    end function

    function Carnot_cooling(T_e,T_c) result (G)
		double precision::T_e,T_c,T_eK,T_cK,Delta_T,G
		T_eK = T_e + 273.15d0
		T_cK = T_c + 273.15d0
		Delta_T = max(3d0,T_cK - T_eK)
		G = T_eK / Delta_T
    end function

    function PartLoadFactor(PLR,a0,a1) result (PLF)
		double precision,intent(in)::PLR,a0,a1
		double precision::PLF
        PLF=PLR/(a0 + a1*PLR + (1d0-a0-a1)*PLR**2) ! PLF = COP / COPfull or EER / EERfull
    end function
    
    ! as_heat_pump_chiller
    
	subroutine as_heat_pump_chiller_initialize(this)
		class(as_heat_pump_chiller)::this
    end subroutine	
    
    ! NOTE COOLING_PLUS_HEAT_RECOVERY mode is foreseen but not managed
	subroutine as_heat_pump_chiller_calculate(this,hour)
		class(as_heat_pump_chiller)::this
		integer::hour
		double precision::EER,COP,Qfull_heat,Qfull_cool
        double precision::ST,b,T0,alpha_COP,alpha_Q
        ! 1. compute the load 
        if (this%component_feature==DYNAMICAL) then
            ! load based on controller logic
			if (this%operation_mode==OFF) then
				this%Qdot_cool=0d0
				this%Qdot_heat=0d0
				this%Wdot_hpc=0d0
                this%PLR_now=0d0
			else if (this%operation_mode==STANDBY) then
				this%Qdot_cool=0d0
				this%Qdot_heat=0d0
				this%Wdot_hpc=this%W_standby			
                this%PLR_now=0d0
			else if (this%operation_mode==COOLING .or. this%operation_mode==COOLING_PLUS_HEAT_RECOVERY) then
				! calculate the cooling load based on outlet water temperature setpoint and mass flow rate control
				this%mdot_w=this%Q_cool_nom/(this%fluid%cp*this%DT_w)             ! for now, only fixed flow rate
				this%T_wo=max(this%T_wi-this%DT_w,this%T_wo_set)   				  ! temperature of outlet water is adjusted based on setpoint 
				this%Qdot_cool=min(0d0,this%mdot_w*this%fluid%cp*(this%T_wo-this%T_wi))    ! the requested cooling load (<0), that can be larger in modulus than actual cooling capacity at this stage
			else if (this%operation_mode==HEATING) then
				! calculate the heating load based on outlet water temperature setpoint and mass flow rate control
				this%mdot_w=this%Q_heat_nom/(this%fluid%cp*this%DT_w)             ! for now, only fixed flow rate
				this%T_wo=min(this%T_wi+this%DT_w,this%T_wo_set)   				  ! temperature of outlet water is adjusted based on setpoint
				this%Qdot_heat=max(0d0,this%mdot_w*this%fluid%cp*(this%T_wo-this%T_wi))    ! the requested heating load, that can be larger than actual heating capacity at this stage
            end if    
        else if (this%component_feature==ALGEBRAIC) then
            ! load is pre-defined (assigned externally)
            this%Qdot_heat=this%Q_heat(hour)
            this%Qdot_cool=this%Q_cool(hour)
            if (this%Qdot_heat>0d0 .AND. this%Qdot_cool==0d0) then
                this%operation_mode=HEATING
            else if (this%Qdot_heat==0d0 .AND. this%Qdot_cool<0d0) then
                this%operation_mode=COOLING
            else if (this%Qdot_heat>0d0 .AND. this%Qdot_cool<0d0) then
                this%operation_mode=COOLING_PLUS_HEAT_RECOVERY
            else     
				this%operation_mode=STANDBY
				this%Q_heat(hour)=0d0
				this%Q_cool(hour)=0d0
				this%Q_aux_h(hour)=0d0
				this%Q_aux_c(hour)=0d0
				this%PLR(hour)=0d0
            end if    
        end if    
		! 2. evaluate thermodynamic model performance  
		select case(this%model_type)
        case(1)
            if (this%operation_mode==HEATING) then
				call heat_pump_model_1_heating (T_e=this%T_air-this%model_par(1),&
												T_c=this%T_wo+this%model_par(2),&
												T_e_nom=this%T_air_h_nom-this%model_par(1),&
												T_c_nom=this%T_hwo_nom+this%model_par(2),&
												Q_heat_nom=this%Q_heat_nom,&
												COP_nom=this%COP_nom,&
												p=this%model_par(3),&
												r=this%model_par(4),&
												a0=this%model_par(5),&
												a1=this%model_par(6),&
												Q_load=this%Qdot_heat,&
												COP=COP,Qfull=Qfull_heat)
				! frost efficiency: COP and Qfull degradation
                b=this%model_par(7)
                T0=this%model_par(8)
                alpha_COP=this%model_par(9)
                alpha_Q=this%model_par(10)
                ST=1d0/(1d0+exp(b*(this%T_air-T0)))
                Qfull_heat=Qfull_heat*(1d0-alpha_Q*ST)
                COP=COP*(1d0-alpha_COP*ST)   
            else if (this%operation_mode==COOLING .or. this%operation_mode==COOLING_PLUS_HEAT_RECOVERY) then
				call heat_pump_model_1_cooling (T_e=this%T_wo-this%model_par(11),&
												T_c=this%T_air+this%model_par(12),&
												T_e_nom=this%T_cwo_nom-this%model_par(11),&
												T_c_nom=this%T_air_c_nom+this%model_par(12),&
												Q_cool_nom=-this%Q_cool_nom,&
												EER_nom=this%EER_nom,&
												p=this%model_par(13),&
												r=this%model_par(14),&
												a0=this%model_par(15),&
												a1=this%model_par(16),&
												Q_load=this%Qdot_cool,&
												EER=EER,Qfull=Qfull_cool)
			end if                
		case default
			call simlog%output('Error: unable to find heat_pump model')
			stop
        end select        
        ! 3. calculate outputs
        if (this%component_feature==DYNAMICAL) then
            ! limit load to available capacity and calculate electricity consumption
            if (this%operation_mode==HEATING) then
				this%Qdot_heat=min(this%Qdot_heat,Qfull_heat)
                this%T_wo=this%T_wi+this%Qdot_heat/(this%mdot_w*this%fluid%cp)
				this%PLR_now=this%Qdot_heat/Qfull_heat
				this%Wdot_hpc=this%Qdot_heat/COP
            else if (this%operation_mode==COOLING .or. this%operation_mode==COOLING_PLUS_HEAT_RECOVERY) then
				this%Qdot_cool=max(this%Qdot_cool,Qfull_cool)
                this%T_wo=this%T_wi+this%Qdot_cool/(this%mdot_w*this%fluid%cp)
				this%PLR_now=this%Qdot_cool/Qfull_cool
				this%Wdot_hpc=-this%Qdot_cool/EER
            end if    
        else if (this%component_feature==ALGEBRAIC) then
			! evaluate partial load ratio and calculate auxiliary heat
            if (this%operation_mode==HEATING) then
				if (this%Qdot_heat<=Qfull_heat) then
					this%Q_aux_h(hour)=0d0
					this%PLR(hour)=this%Q_heat(hour)/Qfull_heat            
				else
					this%Q_aux_h(hour)=this%Q_heat(hour)-Qfull_heat
					this%Q_heat(hour)=Qfull_heat
					this%PLR(hour)=1d0
				end if    
				! calculate electricity consumption
				this%W_hpc(hour)=this%Q_heat(hour)/COP
            else if (this%operation_mode==COOLING .or. this%operation_mode==COOLING_PLUS_HEAT_RECOVERY) then
				if (this%Qdot_cool>=Qfull_cool) then
					this%Q_aux_c(hour)=0d0
					this%PLR(hour)=this%Q_cool(hour)/Qfull_cool            
				else
					this%Q_aux_c(hour)=this%Q_cool(hour)-Qfull_cool
					this%Q_cool(hour)=Qfull_cool
					this%PLR(hour)=1d0
				end if    
				! calculate electricity consumption
				this%W_hpc(hour)=-this%Q_cool(hour)/EER
            end if    
		end if            
	end subroutine	
	
	subroutine as_heat_pump_chiller_update(this,hour,time)
		class(as_heat_pump_chiller)::this
		integer::hour
		double precision,dimension(2)::time ! time(1) new time, time(2) old time
		! update incrementally TIMESTEP-averaged outputs
		this%Q_heat(hour)=(this%Q_heat(hour)*time(2)+this%Qdot_heat*(time(1)-time(2)))/time(1)
		this%Q_cool(hour)=(this%Q_cool(hour)*time(2)+this%Qdot_cool*(time(1)-time(2)))/time(1)
		this%W_hpc(hour)=(this%W_hpc(hour)*time(2)+this%Wdot_hpc*(time(1)-time(2)))/time(1)
		this%PLR(hour)=(this%PLR(hour)*time(2)+this%PLR_now*(time(1)-time(2)))/time(1)
	end subroutine

	subroutine as_heat_pump_chiller_to_json(this,json,obj,obj_name)
		class(as_heat_pump_chiller)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(inout)::obj
		type(json_value),pointer::obj_hp
		character(*)::obj_name
		call json%create_object(obj_hp,obj_name)
		call json%add(obj_hp,'Q_heat_nom',this%Q_heat_nom)	
		call json%add(obj_hp,'COP_nom',this%COP_nom)	
		call json%add(obj_hp,'T_hwo_nom',this%T_hwo_nom)	
		call json%add(obj_hp,'T_air_h_nom',this%T_air_h_nom)	
		call json%add(obj_hp,'Q_cool_nom',this%Q_cool_nom)	
		call json%add(obj_hp,'EER_nom',this%EER_nom)	
		call json%add(obj_hp,'T_cwo_nom',this%T_cwo_nom)	
		call json%add(obj_hp,'T_air_c_nom',this%T_air_c_nom)	
		call json%add(obj_hp,'W_standby',this%W_standby)	
		call this%fluid%to_json(json,obj_hp,'fluid')	
		call json%add(obj_hp,'model_type',this%model_type)					
		call json%add(obj_hp,'model_par',this%model_par)					
		call save_series_bin(json,obj_hp,'Q_heat',this%Q_heat)
		call save_series_bin(json,obj_hp,'Q_cool',this%Q_cool)
		call save_series_bin(json,obj_hp,'PLR',this%PLR)
		call save_series_bin(json,obj_hp,'W_hpc',this%W_hpc)
		call json%add(obj,obj_hp)
    end subroutine		

	subroutine as_heat_pump_chiller_from_json(this,json,obj,path)
		class(as_heat_pump_chiller)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(in)::obj
		character(*)::path
		call json%get(obj,path//'.Q_heat_nom',this%Q_heat_nom)	
		call json%get(obj,path//'.COP_nom',this%COP_nom)	
		call json%get(obj,path//'.T_hwo_nom',this%T_hwo_nom)	
		call json%get(obj,path//'.T_air_h_nom',this%T_air_h_nom)	
		call json%get(obj,path//'.Q_cool_nom',this%Q_cool_nom)	
		call json%get(obj,path//'.EER_nom',this%EER_nom)	
		call json%get(obj,path//'.T_cwo_nom',this%T_cwo_nom)	
		call json%get(obj,path//'.T_air_c_nom',this%T_air_c_nom)	
		call json%get(obj,path//'.W_standby',this%W_standby)	
		call json%get(obj,path//'.model_type',this%model_type)					
		call json%get(obj,path//'.model_par',this%model_par)					
		call this%fluid%from_json(json,obj,path//'.fluid')	
		call read_series_bin(json,obj,path//'.Q_heat',this%Q_heat)
		call read_series_bin(json,obj,path//'.Q_cool',this%Q_cool)
		call read_series_bin(json,obj,path//'.PLR',this%PLR)
		call read_series_bin(json,obj,path//'.W_hpc',this%W_hpc)
    end subroutine		
    
    subroutine as_heat_pump_chiller_get_par(this,par)
		class(as_heat_pump_chiller)::this
        type(param)::par
		select case(par%key(2))
			case ('Q_heat_nom')
				this%Q_heat_nom=par%val_dbl
			case ('COP_nom')
				this%COP_nom=par%val_dbl
			case ('T_hwo_nom')
				this%T_hwo_nom=par%val_dbl
			case ('T_air_h_nom')
				this%T_air_h_nom=par%val_dbl
			case ('Q_cool_nom')
				this%Q_cool_nom=par%val_dbl
			case ('EER_nom')
				this%EER_nom=par%val_dbl
			case ('T_cwo_nom')
				this%T_cwo_nom=par%val_dbl
			case ('T_air_c_nom')
				this%T_air_c_nom=par%val_dbl
			case ('W_standby')
				this%W_standby=par%val_dbl
			case ('model_type')
				this%model_type=par%val_int
			case ('model_par')
				this%model_par=par%arr_dbl
		end select
    end subroutine
    

   ! ww_heat_pump_chiller
    
	subroutine ww_heat_pump_chiller_initialize(this)
		class(ww_heat_pump_chiller)::this
    end subroutine	
    
    ! NOTE COOLING_PLUS_HEAT_RECOVERY mode is foreseen but not managed
	subroutine ww_heat_pump_chiller_calculate(this,hour)
		class(ww_heat_pump_chiller)::this
		integer::hour
		double precision::EER,COP,Qfull_heat,Qfull_cool
        double precision::ST,b,T0,alpha_COP,alpha_Q
        ! 1. compute the load 
        if (this%component_feature==DYNAMICAL) then
            ! load based on controller logic
			if (this%operation_mode==OFF) then
				this%Qdot_cool=0d0
				this%Qdot_heat=0d0
				this%Wdot_hpc=0d0
                this%PLR_now=0d0
			else if (this%operation_mode==STANDBY) then
				this%Qdot_cool=0d0
				this%Qdot_heat=0d0
				this%Wdot_hpc=this%W_standby			
                this%PLR_now=0d0
			else if (this%operation_mode==COOLING .or. this%operation_mode==COOLING_PLUS_HEAT_RECOVERY) then
				! calculate the cooling load based on outlet water temperature setpoint and mass flow rate control
				this%mdot_w=this%Q_cool_nom/(this%fluid%cp*this%DT_w)             ! for now, only fixed flow rate
				this%T_wo=max(this%T_wi-this%DT_w,this%T_wo_set)   				  ! temperature of outlet water is adjusted based on setpoint 
				this%Qdot_cool=min(0d0,this%mdot_w*this%fluid%cp*(this%T_wo-this%T_wi))    ! the requested cooling load (<0), that can be larger in modulus than actual cooling capacity at this stage
			else if (this%operation_mode==HEATING) then
				! calculate the heating load based on outlet water temperature setpoint and mass flow rate control
				this%mdot_w=this%Q_heat_nom/(this%fluid%cp*this%DT_w)             ! for now, only fixed flow rate
				this%T_wo=min(this%T_wi+this%DT_w,this%T_wo_set)   				  ! temperature of outlet water is adjusted based on setpoint
				this%Qdot_heat=max(0d0,this%mdot_w*this%fluid%cp*(this%T_wo-this%T_wi))    ! the requested heating load, that can be larger than actual heating capacity at this stage
            end if    
        else if (this%component_feature==ALGEBRAIC) then
            ! load is pre-defined (assigned externally)
            this%Qdot_heat=this%Q_heat(hour)
            this%Qdot_cool=this%Q_cool(hour)
            if (this%Qdot_heat>0d0 .AND. this%Qdot_cool==0d0) then
                this%operation_mode=HEATING
            else if (this%Qdot_heat==0d0 .AND. this%Qdot_cool<0d0) then
                this%operation_mode=COOLING
            else if (this%Qdot_heat>0d0 .AND. this%Qdot_cool<0d0) then
                this%operation_mode=COOLING_PLUS_HEAT_RECOVERY
            else     
				this%operation_mode=STANDBY
				this%Q_heat(hour)=0d0
				this%Q_cool(hour)=0d0
				this%Q_aux_h(hour)=0d0
				this%Q_aux_c(hour)=0d0
				this%PLR(hour)=0d0
            end if    
        end if    
		! 2. evaluate thermodynamic model performance  
		select case(this%model_type)
        case(1)
            if (this%operation_mode==HEATING) then
				call heat_pump_model_1_heating (T_e=this%T_source_sink-this%model_par(1),&
												T_c=this%T_wo+this%model_par(2),&
												T_e_nom=this%T_cwo_h_nom-this%model_par(1),&
												T_c_nom=this%T_hwo_nom+this%model_par(2),&
												Q_heat_nom=this%Q_heat_nom,&
												COP_nom=this%COP_nom,&
												p=this%model_par(3),&
												r=this%model_par(4),&
												a0=this%model_par(5),&
												a1=this%model_par(6),&
												Q_load=this%Qdot_heat,&
												COP=COP,Qfull=Qfull_heat)
            else if (this%operation_mode==COOLING .or. this%operation_mode==COOLING_PLUS_HEAT_RECOVERY) then
				call heat_pump_model_1_cooling (T_e=this%T_wo-this%model_par(7),&
												T_c=this%T_source_sink+this%model_par(8),&
												T_e_nom=this%T_cwo_nom-this%model_par(7),&
												T_c_nom=this%T_hwo_c_nom+this%model_par(8),&
												Q_cool_nom=-this%Q_cool_nom,&
												EER_nom=this%EER_nom,&
												p=this%model_par(9),&
												r=this%model_par(10),&
												a0=this%model_par(11),&
												a1=this%model_par(12),&
												Q_load=this%Qdot_cool,&
												EER=EER,Qfull=Qfull_cool)
			end if                
		case default
			call simlog%output('Error: unable to find heat_pump model')
			stop
        end select        
        ! 3. calculate outputs
        if (this%component_feature==DYNAMICAL) then
            ! limit load to available capacity and calculate electricity consumption
            if (this%operation_mode==HEATING) then
				this%Qdot_heat=min(this%Qdot_heat,Qfull_heat)
                this%T_wo=this%T_wi+this%Qdot_heat/(this%mdot_w*this%fluid%cp)
				this%PLR_now=this%Qdot_heat/Qfull_heat
				this%Wdot_hpc=this%Qdot_heat/COP
            else if (this%operation_mode==COOLING .or. this%operation_mode==COOLING_PLUS_HEAT_RECOVERY) then
				this%Qdot_cool=max(this%Qdot_cool,Qfull_cool)
                this%T_wo=this%T_wi+this%Qdot_cool/(this%mdot_w*this%fluid%cp)
				this%PLR_now=this%Qdot_cool/Qfull_cool
				this%Wdot_hpc=-this%Qdot_cool/EER
            end if    
        else if (this%component_feature==ALGEBRAIC) then
			! evaluate partial load ratio and calculate auxiliary heat
            if (this%operation_mode==HEATING) then
				if (this%Qdot_heat<=Qfull_heat) then
					this%Q_aux_h(hour)=0d0
					this%PLR(hour)=this%Q_heat(hour)/Qfull_heat            
				else
					this%Q_aux_h(hour)=this%Q_heat(hour)-Qfull_heat
					this%Q_heat(hour)=Qfull_heat
					this%PLR(hour)=1d0
				end if    
				! calculate electricity consumption
				this%W_hpc(hour)=this%Q_heat(hour)/COP
            else if (this%operation_mode==COOLING .or. this%operation_mode==COOLING_PLUS_HEAT_RECOVERY) then
				if (this%Qdot_cool>=Qfull_cool) then
					this%Q_aux_c(hour)=0d0
					this%PLR(hour)=this%Q_cool(hour)/Qfull_cool            
				else
					this%Q_aux_c(hour)=this%Q_cool(hour)-Qfull_cool
					this%Q_cool(hour)=Qfull_cool
					this%PLR(hour)=1d0
				end if    
				! calculate electricity consumption
				this%W_hpc(hour)=-this%Q_cool(hour)/EER
            end if    
		end if            
	end subroutine	
	
	subroutine ww_heat_pump_chiller_update(this,hour,time)
		class(ww_heat_pump_chiller)::this
		integer::hour
		double precision,dimension(2)::time ! time(1) new time, time(2) old time
		! update incrementally TIMESTEP-averaged outputs
		this%Q_heat(hour)=(this%Q_heat(hour)*time(2)+this%Qdot_heat*(time(1)-time(2)))/time(1)
		this%Q_cool(hour)=(this%Q_cool(hour)*time(2)+this%Qdot_cool*(time(1)-time(2)))/time(1)
		this%W_hpc(hour)=(this%W_hpc(hour)*time(2)+this%Wdot_hpc*(time(1)-time(2)))/time(1)
		this%PLR(hour)=(this%PLR(hour)*time(2)+this%PLR_now*(time(1)-time(2)))/time(1)
	end subroutine

	subroutine ww_heat_pump_chiller_to_json(this,json,obj,obj_name)
		class(ww_heat_pump_chiller)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(inout)::obj
		type(json_value),pointer::obj_hp
		character(*)::obj_name
		call json%create_object(obj_hp,obj_name)
		call json%add(obj_hp,'Q_heat_nom',this%Q_heat_nom)	
		call json%add(obj_hp,'COP_nom',this%COP_nom)	
		call json%add(obj_hp,'T_hwo_nom',this%T_hwo_nom)	
		call json%add(obj_hp,'T_cwo_h_nom',this%T_cwo_h_nom)	
		call json%add(obj_hp,'Q_cool_nom',this%Q_cool_nom)	
		call json%add(obj_hp,'EER_nom',this%EER_nom)	
		call json%add(obj_hp,'T_cwo_nom',this%T_cwo_nom)	
		call json%add(obj_hp,'T_hwo_c_nom',this%T_hwo_c_nom)	
		call json%add(obj_hp,'W_standby',this%W_standby)	
		call this%fluid%to_json(json,obj_hp,'fluid')	
		call json%add(obj_hp,'model_type',this%model_type)					
		call json%add(obj_hp,'model_par',this%model_par)					
		call save_series_bin(json,obj_hp,'Q_heat',this%Q_heat)
		call save_series_bin(json,obj_hp,'Q_cool',this%Q_cool)
		call save_series_bin(json,obj_hp,'PLR',this%PLR)
		call save_series_bin(json,obj_hp,'W_hpc',this%W_hpc)
		call json%add(obj,obj_hp)
    end subroutine		

	subroutine ww_heat_pump_chiller_from_json(this,json,obj,path)
		class(ww_heat_pump_chiller)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(in)::obj
		character(*)::path
		call json%get(obj,path//'.Q_heat_nom',this%Q_heat_nom)	
		call json%get(obj,path//'.COP_nom',this%COP_nom)	
		call json%get(obj,path//'.T_hwo_nom',this%T_hwo_nom)	
		call json%get(obj,path//'.T_cwo_h_nom',this%T_cwo_h_nom)	
		call json%get(obj,path//'.Q_cool_nom',this%Q_cool_nom)	
		call json%get(obj,path//'.EER_nom',this%EER_nom)	
		call json%get(obj,path//'.T_cwo_nom',this%T_cwo_nom)	
		call json%get(obj,path//'.T_hwo_c_nom',this%T_hwo_c_nom)	
		call json%get(obj,path//'.W_standby',this%W_standby)	
		call json%get(obj,path//'.model_type',this%model_type)					
		call json%get(obj,path//'.model_par',this%model_par)					
		call this%fluid%from_json(json,obj,path//'.fluid')	
		call read_series_bin(json,obj,path//'.Q_heat',this%Q_heat)
		call read_series_bin(json,obj,path//'.Q_cool',this%Q_cool)
		call read_series_bin(json,obj,path//'.PLR',this%PLR)
		call read_series_bin(json,obj,path//'.W_hpc',this%W_hpc)
    end subroutine		
    
    subroutine ww_heat_pump_chiller_get_par(this,par)
		class(ww_heat_pump_chiller)::this
        type(param)::par
		select case(par%key(2))
			case ('Q_heat_nom')
				this%Q_heat_nom=par%val_dbl
			case ('COP_nom')
				this%COP_nom=par%val_dbl
			case ('T_hwo_nom')
				this%T_hwo_nom=par%val_dbl
			case ('T_cwo_h_nom')
				this%T_cwo_h_nom=par%val_dbl
			case ('Q_cool_nom')
				this%Q_cool_nom=par%val_dbl
			case ('EER_nom')
				this%EER_nom=par%val_dbl
			case ('T_cwo_nom')
				this%T_cwo_nom=par%val_dbl
			case ('T_hwo_c_nom')
				this%T_hwo_c_nom=par%val_dbl
			case ('W_standby')
				this%W_standby=par%val_dbl
			case ('model_type')
				this%model_type=par%val_int
			case ('model_par')
				this%model_par=par%arr_dbl
		end select
    end subroutine
     
   ! ww_heat_pump 
    
	subroutine ww_heat_pump_initialize(this)
		class(ww_heat_pump)::this
    end subroutine		
    
	subroutine ww_heat_pump_calculate(this,hour)
		class(ww_heat_pump)::this
		double precision::Q_cond_full,COP
		integer::hour
        ! 1. compute the load 
        if (this%component_feature==DYNAMICAL) then
            ! load based on controller logic
			if (this%operation_mode==OFF) then
				this%Qdot_cond=0d0
				this%Qdot_evap=0d0
				this%Wdot_comp=0d0
                this%PLR_now=0d0
			else if (this%operation_mode==STANDBY) then
				this%Qdot_cond=0d0
				this%Qdot_evap=0d0
				this%Wdot_comp=this%W_standby			
                this%PLR_now=0d0
			else if (this%operation_mode==HEATING) then
				! calculate the heating load based on outlet water temperature setpoint and mass flow rate control
				this%mdot_hw=this%Q_heat_nom/(this%fluid%cp*this%DT_hw)           ! for now, only fixed flow rate
				this%T_hwo=min(this%T_hwi+this%DT_hw,this%T_hwo_set)                ! temperature of outlet water is adjusted based on setpoint
				this%Qdot_cond=max(0d0,this%mdot_hw*this%fluid%cp*(this%T_hwo-this%T_hwi))  ! the requested heating load, that can be larger than actual heating capacity at this stage
            end if    
        else if (this%component_feature==ALGEBRAIC) then
            ! load is pre-defined (assigned externally)
            this%Qdot_cond=this%Q_cond(hour)
            if (this%Qdot_cond>0d0) then
                this%operation_mode=HEATING
            else                  
   				this%operation_mode=STANDBY
				this%Q_evap(hour)=0d0
				this%W_comp(hour)=this%W_standby
            	this%Q_aux_h(hour)=0d0
				this%PLR(hour)=0d0
            end if    
        end if    
		! 2. evaluate thermodynamic model performance          
		select case(this%model_type)
        case(1)
            if (this%operation_mode==HEATING) then           
				call heat_pump_model_1_heating (T_e=this%T_cwo-this%model_par(1),&
												T_c=this%T_hwo+this%model_par(2),&
												T_e_nom=this%T_cwo_nom-this%model_par(1),&
												T_c_nom=this%T_hwo_nom+this%model_par(2),&
												Q_heat_nom=this%Q_heat_nom,&
												COP_nom=this%COP_nom,&
												p=this%model_par(3),&
												r=this%model_par(4),&
												a0=this%model_par(5),&
												a1=this%model_par(6),&
												Q_load=this%Qdot_cond,&
												COP=COP,Qfull=Q_cond_full)
             end if   
        case default
            call simlog%output('Error: unable to find heat_pump model')
            stop
        end select  
        ! 3. calculate outputs
        if (this%component_feature==DYNAMICAL) then
            ! limit load to available capacity and calculate electricity consumption
            if (this%operation_mode==HEATING) then
				this%Qdot_cond=min(this%Qdot_cond,Q_cond_full)
                this%T_hwo=this%T_hwi+this%Qdot_cond/(this%mdot_hw*this%fluid%cp)
				this%PLR_now=this%Qdot_cond/Q_cond_full
				this%Wdot_comp=this%Qdot_cond/COP
            end if
		else if (this%component_feature==ALGEBRAIC) then
			! evaluate partial load ratio and calculate auxiliary heat
			if (this%operation_mode==HEATING) then
				if (this%Qdot_cond<=Q_cond_full) then
					this%Q_aux_h(hour)=0d0
					this%PLR(hour)=this%Q_cond(hour)/Q_cond_full            
				else
					this%Q_aux_h(hour)=this%Q_cond(hour)-Q_cond_full
					this%Q_cond(hour)=Q_cond_full
					this%PLR(hour)=1d0
				end if    
				! calculate electricity consumption and evaporator cooling power
				this%W_comp(hour)=this%Q_cond(hour)/COP
				this%Q_evap(hour)=this%Q_cond(hour)-this%W_comp(hour)	
			end if                        
		end if                     
    end subroutine
    
	subroutine ww_heat_pump_update(this,hour,time)
		class(ww_heat_pump)::this
		integer::hour
		double precision,dimension(2)::time ! time(1) new time, time(2) old time
		! update incrementally TIMESTEP-averaged outputs
		this%Q_cond(hour)=(this%Q_cond(hour)*time(2)+this%Qdot_cond*(time(1)-time(2)))/time(1)
		this%Q_evap(hour)=(this%Q_evap(hour)*time(2)+this%Qdot_evap*(time(1)-time(2)))/time(1)
		this%W_comp(hour)=(this%W_comp(hour)*time(2)+this%Wdot_comp*(time(1)-time(2)))/time(1)
		this%PLR(hour)=(this%PLR(hour)*time(2)+this%PLR_now*(time(1)-time(2)))/time(1)
	end subroutine
    
	subroutine ww_heat_pump_to_json(this,json,obj,obj_name)
		class(ww_heat_pump)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(inout)::obj
		type(json_value),pointer::obj_ww
		character(*)::obj_name
		call json%create_object(obj_ww,obj_name)
		call json%add(obj_ww,'Q_heat_nom',this%Q_heat_nom)	
		call json%add(obj_ww,'COP_nom',this%COP_nom)	
		call json%add(obj_ww,'T_hwo_nom',this%T_hwo_nom)	
		call json%add(obj_ww,'T_cwo_nom',this%T_cwo_nom)
        call json%add(obj_ww,'W_standby',this%W_standby)	        
		call json%add(obj_ww,'model_type',this%model_type)					
		call json%add(obj_ww,'model_par',this%model_par)					
		call save_series_bin(json,obj_ww,'Q_cond',this%Q_cond)
		call save_series_bin(json,obj_ww,'PLR',this%PLR)
		call save_series_bin(json,obj_ww,'Q_evap',this%Q_evap)
		call save_series_bin(json,obj_ww,'W_comp',this%W_comp)
		call save_series_bin(json,obj_ww,'Q_aux_h',this%Q_aux_h)
		call json%add(obj,obj_ww)
    end subroutine		

	subroutine ww_heat_pump_from_json(this,json,obj,path)
		class(ww_heat_pump)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(in)::obj
		character(*)::path
		call json%get(obj,path//'.Q_heat_nom',this%Q_heat_nom)
		call json%get(obj,path//'.COP_nom',this%COP_nom)
		call json%get(obj,path//'.T_hwo_nom',this%T_hwo_nom)
		call json%get(obj,path//'.T_cwo_nom',this%T_cwo_nom)
		call json%get(obj,path//'.W_standby',this%W_standby)	
		call json%get(obj,path//'.model_type',this%model_type)					
		call json%get(obj,path//'.model_par',this%model_par)					
		call read_series_bin(json,obj,path//'.Q_cond',this%Q_cond)
		call read_series_bin(json,obj,path//'.PLR',this%PLR)
		call read_series_bin(json,obj,path//'.Q_evap',this%Q_evap)
		call read_series_bin(json,obj,path//'.W_comp',this%W_comp)
		call read_series_bin(json,obj,path//'.Q_aux_h',this%Q_aux_h)
    end subroutine		
    
    subroutine ww_heat_pump_get_par(this,par)
		class(ww_heat_pump)::this
        type(param)::par
		select case(par%key(2))
			case ('Q_heat_nom')
				this%Q_heat_nom=par%val_dbl
			case ('COP_nom')
				this%COP_nom=par%val_dbl
			case ('T_hwo_nom')
				this%T_hwo_nom=par%val_dbl
			case ('T_cwo_nom')
				this%T_cwo_nom=par%val_dbl
			case ('W_standby')
				this%W_standby=par%val_dbl
			case ('model_type')
				this%model_type=par%val_int
			case ('model_par')
				this%model_par=par%arr_dbl
		end select
    end subroutine
    
	subroutine dhw_storage_initialize(this)
		class(dhw_storage)::this
    end subroutine

    subroutine dhw_storage_calculate(this,hour)
		class(dhw_storage)::this
		integer::hour
		double precision::f_dhw_on_new,Q_dhw_out,T_tank_dhw_new,Q_heater_dhw_new
		integer::previous_hour
		previous_hour=previous_step(hour)
		if (this%volume>0) then
			f_dhw_on_new=0d0
			Q_dhw_out=this%mdotcp_dhw(hour)*(40d0-this%T_main(hour))*TIMESTEP+this%UA*(this%T_storage(previous_hour)-this%T_loss(hour))*TIMESTEP
			! calculate new temperature in the DHW storage tank at the end of time step using current control signal
			Q_heater_dhw_new=this%f_dhw_on*this%Q_heater
			T_tank_dhw_new=this%T_storage(previous_hour)+(Q_heater_dhw_new*TIMESTEP-Q_dhw_out)/(this%volume*RHO_WATER*CP_WATER)
			! calculate fraction of timestep during which the heater is on (f_dhw_on_new) so that temperature is within bounds 
			! decreasing temperatures and new temperature below lower setpoint 
			if (T_tank_dhw_new<this%T_lower_dhw .and. this%T_storage(previous_hour)>T_tank_dhw_new) f_dhw_on_new=(this%T_lower_dhw-T_tank_dhw_new)/(this%T_storage(previous_hour)-T_tank_dhw_new)
			! increasing temperatures and new temperature above upper setpoint 
			if (T_tank_dhw_new>this%T_upper_dhw .and. this%T_storage(previous_hour)<T_tank_dhw_new) f_dhw_on_new=1d0-(T_tank_dhw_new-this%T_upper_dhw)/(T_tank_dhw_new-this%T_storage(previous_hour))
			! recalculate storage temperature if fraction of time during which the dhw tank is charged (f_dhw_on) has changed
			f_dhw_on_new=max(0d0,min(1d0,f_dhw_on_new))
			if (this%f_dhw_on==0 .and. f_dhw_on_new > 0) then
				 ! set heating energy to the minimum between the current value and the value giving T upper at the end of time step
				Q_heater_dhw_new=min(f_dhw_on_new*this%Q_heater,((this%volume*RHO_WATER*CP_WATER)*(this%T_upper_dhw-this%T_storage(previous_hour))+Q_dhw_out)/TIMESTEP)
				this%T_storage(hour)=this%T_storage(previous_hour)+(Q_heater_dhw_new*TIMESTEP-Q_dhw_out)/(this%volume*RHO_WATER*CP_WATER)
				! set control signal for next time step
				this%f_dhw_on=1d0
			else if (this%f_dhw_on==1 .and. f_dhw_on_new < 1) then
				 ! set heating energy to the maximum between the current value and the value giving T lower at the end of time step
				Q_heater_dhw_new=max(f_dhw_on_new*this%Q_heater,((this%volume*RHO_WATER*CP_WATER)*(this%T_lower_dhw-this%T_storage(previous_hour))+Q_dhw_out)/TIMESTEP)
				this%T_storage(hour)=this%T_storage(previous_hour)+(Q_heater_dhw_new*TIMESTEP-Q_dhw_out)/(this%volume*RHO_WATER*CP_WATER)
				! set control signal for next time step
				this%f_dhw_on=0d0   
			else
				this%T_storage(hour)=T_tank_dhw_new
			end if
		else
			! tankless DHW preparation
			Q_heater_dhw_new=this%mdotcp_dhw(hour)*(40d0-this%T_main(hour))
        end if 	
        this%Q_supply(hour)=Q_heater_dhw_new
	end subroutine

    subroutine dhw_storage_to_json(this,json,obj,obj_name)
		class(dhw_storage)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(inout)::obj
		type(json_value),pointer::obj_dhw
		character(*)::obj_name
		call json%create_object(obj_dhw,obj_name)
		call json%add(obj_dhw,'volume',this%volume)
		call json%add(obj_dhw,'UA',this%UA)
		call json%add(obj_dhw,'Q_heater',this%Q_heater)
		call save_series_bin(json,obj_dhw,'mdotcp_dhw',this%mdotcp_dhw)
		call save_series_bin(json,obj_dhw,'T_main',this%T_main)
		call save_series_bin(json,obj_dhw,'T_storage',this%T_storage)
		call save_series_bin(json,obj_dhw,'Q_supply',this%Q_supply)
		call json%add(obj,obj_dhw)
	end subroutine

    subroutine dhw_storage_from_json(this,json,obj,path)
		class(dhw_storage)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(in)::obj
		character(*)::path
		call json%get(obj,path//'.volume',this%volume)
		call json%get(obj,path//'.UA',this%UA)
		call json%get(obj,path//'.Q_heater',this%Q_heater)
		call read_series_bin(json,obj,path//'.mdotcp_dhw',this%mdotcp_dhw)
		call read_series_bin(json,obj,path//'.T_main',this%T_main)
		call read_series_bin(json,obj,path//'.T_storage',this%T_storage)
		call read_series_bin(json,obj,path//'.Q_supply',this%Q_supply)
    end subroutine

    subroutine buffer_controller_to_json(this,json,obj,obj_name)
		class(buffer_controller)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(inout)::obj
		type(json_value),pointer::obj_bcontol
		character(*)::obj_name
		call json%create_object(obj_bcontol,obj_name)
		call save_series_bin(json,obj_bcontol,'T_storage_min',this%T_storage_min)
		call save_series_bin(json,obj_bcontol,'T_storage_max',this%T_storage_max)
		call json%add(obj,obj_bcontol)
	end subroutine
    
    subroutine buffer_controller_from_json(this,json,obj,path)
		class(buffer_controller)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(in)::obj
		character(*)::path
		call read_series_bin(json,obj,path//'.T_storage_min',this%T_storage_min)
		call read_series_bin(json,obj,path//'.T_storage_max',this%T_storage_max)
    end subroutine
        
	function borehole_heat_exchanger_Qex(this,hour,mdot_f,T_fi,T_fo) result (Qex)
		class(borehole_heat_exchanger)::this
		double precision::mdot_f,T_fi,T_fo
		integer::hour
		double precision::Qex	
		double precision::eps,NTU,UA,Cmin,Tb
		Tb=this%Tg+this%DT_g(hour)! borehole temperature
        UA=this%N_bhe*this%H/(this%R_bhe_cond+this%convective_resistance(this%Reynolds(mdot_f,T_fo),T_fo))	
		Cmin=abs(mdot_f)*this%cp_f	
		if (Cmin>0.2d0*UA) then
			NTU=UA/Cmin
			eps=1d0-exp(-NTU)
			Qex=eps*Cmin*(Tb-T_fi)
		else
			Qex=UA*(Tb-T_fo)
		end if
		this%Qdot_ground=-Qex
	end function
	
	subroutine borehole_heat_exchanger_update_DT_ground(this,hour,time,DQ_out) 
		class(borehole_heat_exchanger)::this
		double precision::DQ_out,time,alpha
		integer::hour,next_hour,prev_hour
        integer::i,mode
        ! heat in input to the ground per unit length (W/m)
		this%DQ_g(cumulative_step(hour))=DQ_out/(this%N_bhe*this%H) 
		! update ground temperature at the borehole boundary using the g-function
		next_hour=next_step(hour)
		prev_hour=previous_step(hour)
		this%DT_g(hour)=this%g_function%get(time)*this%DQ_g(cumulative_step(hour))/(2d0*pi*this%lambda) ! contribution of current timestep
		if (time==TIMESTEP) then
			this%DT_g(hour)=this%DT_g(hour)+this%DT_g_next			
			! contribution of previous timesteps at the end of next timestep
			this%DT_g_next=0d0
			do i=1,cumulative_step(hour)
				this%DT_g_next=this%DT_g_next+this%g_function%get(i*3600d0)*this%DQ_g(cumulative_step(hour)+1-i)/(2d0*pi*this%lambda)
			end do	
		else	
			! add contribution of previous timesteps as linearized term of the respective DT evaluated at the beginning and the end of current timestep
			alpha=time/TIMESTEP
			this%DT_g(hour)=this%DT_g(hour)+this%DT_g(prev_hour)*(1d0-alpha)+this%DT_g_next*alpha
		end if
		if (time==TIMESTEP) then
			! carry forward DT_g to next hour
			this%DT_g(next_hour)=this%DT_g(hour)
		end if
	end subroutine
	
    subroutine borehole_heat_exchanger_update(this,hour,time)
        class(borehole_heat_exchanger)::this
        integer::hour
		double precision,dimension(2)::time
		! timestep averaged Q_ground
		this%Q_ground(hour)=(this%Q_ground(hour)*time(2)+this%Qdot_ground*(time(1)-time(2)))/time(1)
        ! calculate the temperature variation at the borehole perimeter
		call this%update_DT_ground(hour,&   							 ! the current hour
				time(1),&						   						 ! the current time
				this%Q_ground(hour)-this%Q_ground(previous_step(hour))&  ! the variation in Q injected in the ground between hour and hour-1
				) 
	end subroutine

	subroutine borehole_heat_exchanger_to_json(this,json,obj,obj_name)
		class(borehole_heat_exchanger)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(inout)::obj
		type(json_value),pointer::obj_bhe
		type(json_value),pointer::g_function_obj,mat,row	
		character(*)::obj_name
        integer::i
		call json%create_object(obj_bhe,obj_name)
		call json%add(obj_bhe,'H',this%H)
        call json%create_array(mat,'bhe_XY_coordinates')
		do i = 1, size(this%bhe_XY_coordinates,1)
			call json%add(mat, '',this%bhe_XY_coordinates(i,1:2))
		end do
		call json%add(obj_bhe,mat)        
		call json%add(obj_bhe,'Tg',this%Tg)
		call save_series_bin(json,obj_bhe,'Q_ground',this%Q_ground)
		call save_series_bin(json,obj_bhe,'DT_g',this%DT_g)
		call json%create_object(g_function_obj,'g_function')
		if (allocated(this%g_function%tau).and.allocated(this%g_function%gval)) then
			call json%add(g_function_obj,'t_s',this%g_function%t_s)
			call json%add(g_function_obj,'tau',this%g_function%tau)
			call json%add(g_function_obj,'gval',this%g_function%gval)
		end if	
		call json%add(obj_bhe,g_function_obj)
		call json%add(obj,obj_bhe)		
    end subroutine			

	subroutine borehole_heat_exchanger_from_json(this,json,obj,path)
		class(borehole_heat_exchanger)::this
		type(json_core),intent(inout)::json
		type(json_value),pointer,intent(in)::obj
		character(*)::path
		call read_series_bin(json,obj,path//'.Q_ground',this%Q_ground)
		call read_series_bin(json,obj,path//'.DT_g',this%DT_g)
		call json%get(obj,path//'.g_function.t_s',this%g_function%t_s)
		call json%get(obj,path//'.g_function.tau',this%g_function%tau)
		call json%get(obj,path//'.g_function.gval',this%g_function%gval)
    end subroutine			
        
	! load data in the the interpolation class
	subroutine g_function_load(this,tau_array,gval_array)
		class(g_function)::this
		double precision,dimension(:)::tau_array,gval_array
		integer::flag
        allocate(this%interp)
        if (allocated(this%tau)) deallocate(this%tau)
        allocate(this%tau(size(tau_array)))
        if (allocated(this%gval)) deallocate(this%gval)
        allocate(this%gval(size(tau_array)))
		this%tau=tau_array
		this%gval=gval_array
		call this%interp%initialize(this%tau,this%gval,flag)
		if (flag>0) then
			call simlog%output('Error: g_function_load failed.')
			stop
		end if	
	end subroutine
	
	! get the g-function value at a given time (seconds)
	function g_function_get(this,time) result(gval)
		class(g_function)::this
		double precision::time
		double precision::gval
		call this%interp%evaluate(time/this%t_s,gval)
	end function	
	
    subroutine borehole_heat_exchanger_initialize(this)
        class(borehole_heat_exchanger),target::this
		double precision::R_p,R_grout		
		this%DT_g=0d0
        allocate(this%DQ_g(8760*size(simulation_years))) ! NOTE the full history of DQ_g's must be stored in multi-year simulation since hour 1 of year 1.
		this%DQ_g=0d0
		! pipe conductive resistance (m*K/W)
		R_p=log(this%D_p/(this%D_p-2d0*this%s_p))/(2d0*pi*this%lambda_p*this%bhe_type*2d0) 	
		! grout thermal resistance
		if (this%bhe_type==1) then
			! single-U pipe (m*K/W)
			R_grout=1d0/(4d0*pi*this%lambda_grout)*(log(this%D_bhe/this%D_p)+log(this%D_bhe/this%D)+(this%lambda_grout-this%lambda)/(this%lambda_grout+this%lambda)*log(1/(1-(this%D/this%D_bhe)**4d0))) 
		else
			! double-U pipe (m*K/W)
			R_grout=1d0/(2d0*pi*this%lambda_grout)*(log(this%D_bhe/this%D_p)-3d0/4d0+(this%D/this%D_bhe)**2d0-1d0/4d0*log(1d0-(this%D/this%D_bhe)**8d0)-1d0/2d0*log(2d0**0.5d0*this%D/this%D_p)-1d0/4d0*log(2d0*this%D/this%D_p)) 
		end if
		! calculate overall conductive thermal resistance
		this%R_bhe_cond=R_p+R_grout		
        if (this%layout_mode==AUTO_DESIGN) then
			call this%design()
        end if    
        call this%create_g_function(size(simulation_years))
    end subroutine
	
	function dynamic_viscosity(T) result (mu)
		double precision::T  ! °C
		double precision::mu ! Pa*s
		double precision::T0
		T0=min(max(T,0d0),50d0) ! validity range of correlations 0 - 50°C
		mu=2.414d-5*10**(247.8d0/(T0-140d0))
	end function
	
	function borehole_heat_exchanger_Reynolds(this,mdot,Tf) result (Re)
		class(borehole_heat_exchanger)::this
		double precision::mdot ! overall mass flow rate in the borehole field
		double precision::Tf ! mean fluid temperature 
		double precision::Re ! Reynolds number
		double precision::mu ! dynamic viscosity
		mu=dynamic_viscosity(Tf)
		Re=4d0*(abs(mdot)/this%N_bhe)/(this%bhe_type*pi*mu*(this%D_p-2d0*this%s_p))
 	end function
	
	function borehole_heat_exchanger_convective_resistance(this,Re,Tf) result(R_cv)
		class(borehole_heat_exchanger)::this
		double precision::Re ! Reynolds number in the pipes
		double precision::Tf,Tf0 ! mean fluid temperature 
		double precision::Nu_lam,f_3000,Nu_3000,Nu,w,h,Pr,Lambda,f,R_cv
		Tf0=min(max(Tf,0d0),50d0) ! validity range of correlations 0 - 50°C
		Pr = 113.471d0-0.4729d0*Tf0+0.008946d0*Tf0**2-5.30d-5*Tf0**3 ! Prandtl of water in °C
		Lambda=0.561d0+1.90d-3*Tf0-1.60d-5*Tf0**2 ! thermal conductivity of water
		Nu_lam=4.36d0
		f_3000=(0.79d0*log(3000d0)-1.64d0)**(-2d0)
		Nu_3000=(f_3000/8d0)*(3000d0-1000d0)*Pr/(1d0+12.7d0*(f_3000/8d0)**0.5d0*(Pr**(2d0/3d0)-1d0))
		if (Re<=2300d0) then
			Nu=Nu_lam
		else if (Re>=3000d0)then
			f=(0.79d0*log(Re)-1.64d0)**(-2d0) ! Darcy friction coefficient
			Nu=(f/8d0)*(Re-1000d0)*Pr/(1d0+12.7d0*(f/8d0)**0.5*(Pr**(2d0/3d0)-1d0)) ! Gnielinsky correlation 3000<Re<5E6
		else if ((2300d0<Re) .and. (Re<3000d0)) then
			w=(Re-2300d0)/(3000d0-2300d0);
			Nu=Nu_lam*(1-w)+Nu_3000*w;
		end if
		h=Nu*Lambda/(this%D_p-2d0*this%s_p) ! internal convective coefficient [W/(m2.K)]
		R_cv=1/(pi*(this%D_p-2*this%s_p)*h*this%bhe_type*2d0) ! convective resistance [m.K/W]
	end function 
	
    subroutine borehole_heat_exchanger_design(this)
        ! Assumptions:
        ! 1) groundwater flow absent
        ! 2) water as fluid in the bhe
        ! 3) single or double U pipes
        ! 4) BHE layout is square grid with compact rectangular layout
        ! 5) boreholes in parallel
    
        class(borehole_heat_exchanger)::this
        integer::ih
        integer::ic
        integer::i,j

        double precision,dimension(:),allocatable::Xs
        double precision,dimension(:),allocatable::Ys
        double precision,dimension(:),allocatable::YY 
        double precision,dimension(10,3)::a 		! coefficients for ICS solution

        
        double precision::Re
        double precision::R_cv
        double precision::R_6H
        double precision::R_1M
        double precision::R_10Y
        
        double precision::L_h0
        double precision::L_c0

        integer::Nx0
        integer::Ny0
        integer::N_odd0
        double precision::H_0
        
        double precision::m_h
        double precision::m_c
        
        double precision::R_cv_h
        double precision::R_cv_c
        double precision::R_bhe_h
        double precision::R_bhe_c
        
        double precision::R_cv_h_old
        double precision::R_cv_c_old
        double precision::R_bhe_h_old
        double precision::R_bhe_c_old
        
        double precision::r_max 
        integer::N_cc
        double precision::Dr 
        double precision,dimension(:),allocatable::r_cc
        double precision,dimension(:),allocatable::r_icc
        integer::t_10y
        double precision,dimension(:),allocatable::X
        double precision,dimension(:),allocatable::I_temp
        double precision::Q_stored
        double precision::Tp1
        
        double precision::tol
        double precision::L_old
        integer::iteration
        
        double precision::L_h
        double precision::L_c
        
        double precision::Tp
        integer::N1
        integer::N2
        integer::N3
        integer::N4
        
        double precision::alpha
		double precision::Th,Tc,mu_h,mu_c,Re_h,Re_c,m_hbhe,m_cbhe,R_bhe_res
                
        
        ! matrix of the coefficients for ground resistances correlation based on
		! ICS solution, [f6h, f1m, f10y]		
		a(1,:)=[0.6619352d0,0.4132728d0,0.3057646d0]
		a(2,:)=[-4.815693d0,0.2912981d0,0.08987446d0]
		a(3,:)=[15.03571d0,0.07589286d0,-0.09151786d0]
		a(4,:)=[-0.09879421d0,0.1563978d0,-0.03872451d0]
		a(5,:)=[0.02917889d0,-0.2289355d0,0.1690853d0]
		a(6,:)=[0.1138498d0,-0.004927554d0,-0.02881681d0]
		a(7,:)=[0.005610933d0,-0.002694979d0,-0.002886584d0]
		a(8,:)=[0.7796329d0,-0.6380360d0,-0.1723169d0]
		a(9,:)=[-0.3243880d0,0.2950815d0,0.03112034d0]
		a(10,:)=[-0.01824101d0,0.1493320d0,-0.1188438d0]	

		Th=(this%T_fih_design+this%T_foh_design)/2d0
		Tc=(this%T_fic_design+this%T_foc_design)/2d0
		
		mu_h=dynamic_viscosity(Th)
		mu_c=dynamic_viscosity(Tc)
				
		Re=4000d0
		R_cv=this%convective_resistance(Re,Th)		
		R_bhe_res=this%R_bhe_cond+R_cv
		
		alpha=this%lambda/this%C*24d0*3600d0 ! thermal diffusivity [m2/day]
		R_6h=1d0/this%lambda*(a(1,1)+a(2,1)*(this%D_bhe/2d0)+a(3,1)*(this%D_bhe/2d0)**2d0+a(4,1)*alpha+a(5,1)*alpha**2d0+a(6,1)*log(alpha)+a(7,1)*(log(alpha))**2d0+a(8,1)*(this%D_bhe/2d0)*alpha+a(9,1)*(this%D_bhe/2d0)*log(alpha)+a(10,1)*alpha*log(alpha)) ! 6h pulse resistance
		R_1m=1d0/this%lambda*(a(1,2)+a(2,2)*(this%D_bhe/2d0)+a(3,2)*(this%D_bhe/2d0)**2d0+a(4,2)*alpha+a(5,2)*alpha**2d0+a(6,2)*log(alpha)+a(7,2)*(log(alpha))**2d0+a(8,2)*(this%D_bhe/2d0)*alpha+a(9,2)*(this%D_bhe/2d0)*log(alpha)+a(10,2)*alpha*log(alpha)) ! monthly pulse resistance 
		R_10y=1d0/this%lambda*(a(1,3)+a(2,3)*(this%D_bhe/2d0)+a(3,3)*(this%D_bhe/2d0)**2d0+a(4,3)*alpha+a(5,3)*alpha**2d0+a(6,3)*log(alpha)+a(7,3)*(log(alpha))**2d0+a(8,3)*(this%D_bhe/2d0)*alpha+a(9,3)*(this%D_bhe/2d0)*log(alpha)+a(10,3)*alpha*log(alpha)) ! 10y pulse resistance

		L_h=(this%Q_ga*R_10y+this%Q_ghm*R_1m+this%Q_ghd*(R_bhe_res+R_6h))/(this%Tg-((this%T_fih_design+this%T_foh_design)/2d0+273.15d0)+273.15d0) ! total bhe length heating mode, first guess [m]
		L_c=(this%Q_ga*R_10y+this%Q_gcm*R_1m+this%Q_gcd*(R_bhe_res+R_6h))/(this%Tg-((this%T_fic_design+this%T_foc_design)/2d0+273.15d0)+273.15d0) ! total bhe length cooling mode, first guess [m]
		this%L_new=max(L_h,L_c) ! total bhe length, first guess [m]

		call layout(this%N_bhe,this%Nx,this%Ny,this%N_odd,this%H,this%L_new,this%H_des) ! preliminary sizing, neglecting interaction penalty

		m_h=this%Q_ghd/(this%cp_f*(this%T_foh_design-this%T_fih_design)) ! total mass flow rate heating mode [kg/s]
		m_c=this%Q_gcd/(this%cp_f*(this%T_foc_design-this%T_fic_design)) ! total mass flow rate cooling mode [kg/s]
        !m_hbhe=m_h/this%N_bhe
		!m_cbhe=m_c/this%N_bhe
        
		!Re_h=4d0*m_hbhe/(this%bhe_type*pi*mu_h*(this%D_p-2d0*this%s_p))
		!Re_c=4d0*m_cbhe/(this%bhe_type*pi*mu_c*(this%D_p-2d0*this%s_p))
		R_cv_h = this%convective_resistance(this%Reynolds(m_h,Th),Th)
		R_cv_c = this%convective_resistance(this%Reynolds(m_c,Tc),Tc)
		R_bhe_h=this%R_bhe_cond+R_cv_h
		R_bhe_c=this%R_bhe_cond+R_cv_c
		
		r_max=10d0 ! distance from each bhe where temperature perturbation arrives in 10y [m] 
		N_cc=10 ! number of circular sectors
		Dr=(r_max-this%B_des/2d0)/N_cc ! circular crowns thickness [m]
		allocate(r_cc(floor(abs(((this%B_des/2d0+Dr/2d0)-(r_max-Dr/2d0))/Dr))+1))
		allocate(r_icc(floor(abs((this%B_des/2d0-(r_max-Dr))/Dr))+1))
		r_cc(1)=(this%B_des/2+Dr/2) ! average radii of the circular crowns [m]
		do i=2,size(r_cc)
			r_cc(i)=r_cc(i-1)+Dr
		end do
		
		r_icc(1)=this%B_des/2d0 ! inner radii of the circular crowns [m]
		do i=2,size(r_icc)
			r_icc(i)=r_icc(i-1)+Dr
		end do
		
		t_10y=3650 ! [days]
		allocate(X(size(r_cc)))
		allocate(I_temp(size(r_icc)))
        I_temp=0d0
		X=r_cc/(2d0*(alpha*t_10y)**0.5d0) ! non-dimensional variable for ILS solution (correlation from UNI 11466)
		do i=1,size(X)
		    if (0.01d0<=X(i) .and. X(i)<0.5d0) then
		        I_temp(i)=-0.932992d0*log(X(i))-0.14601d0
		    else if (0.5d0<=X(i) .and. X(i)<=1d0) then
				I_temp(i)=-0.577078d0*log(X(i))+0.1d0
		    end if
		end do
		
		Q_stored=this%C*sum((this%L_new*pi*((r_icc+Dr)**2d0-r_icc**2d0))*(this%Q_ga*I_temp/(2d0*pi*this%lambda*this%L_new))) ! heat stored in the annulus between B_des/2 and r_max in 10 ys [J]
		Tp1=Q_stored/(this%C*this%B_des**2d0*this%L_new) ! temperature penalty for a bhe surrounded by others on the 4 sides [K]
		
		call classBHE(N1,N2,N3,N4,this%N_bhe,this%Nx,this%Ny,this%N_odd)
		Tp=penalty(N1,N2,N3,N4,this%N_bhe,Tp1)
		
		L_h=(this%Q_ga*R_10y+this%Q_ghm*R_1m+this%Q_ghd*(R_bhe_h+R_6h))/(this%Tg+273.15d0-((this%T_fih_design+this%T_foh_design)/2d0+273.15d0)-Tp) ! total bhe length heating mode, first guess [m]
		L_c=(this%Q_ga*R_10y+this%Q_gcm*R_1m+this%Q_gcd*(R_bhe_c+R_6h))/(this%Tg+273.15d0-((this%T_fic_design+this%T_foc_design)/2d0+273.15d0)-Tp) ! total bhe length cooling mode, first guess [m]
		this%L_new=maxval([L_h,L_c])! total bhe length, including temperature penalty effect [m]
		L_old=this%L_new
		
		call layout(this%N_bhe,this%Nx,this%Ny,this%N_odd,this%H,this%L_new,this%H_des) ! 1st sizing, including penalty

		call classBHE(N1,N2,N3,N4,this%N_bhe,this%Nx,this%Ny,this%N_odd)
		Tp = penalty(N1,N2,N3,N4,this%N_bhe,Tp1) ! penalty is recalculated

		L_h=(this%Q_ga*R_10y+this%Q_ghm*R_1m+this%Q_ghd*(R_bhe_h+R_6h))/(this%Tg+273.15d0-((this%T_fih_design+this%T_foh_design)/2d0+273.15d0)-Tp) ! total bhe length heating mode, 2nd guess [m]
		L_c=(this%Q_ga*R_10y+this%Q_gcm*R_1m+this%Q_gcd*(R_bhe_c+R_6h))/(this%Tg+273.15d0-((this%T_fic_design+this%T_foc_design)/2d0+273.15d0)-Tp) ! total bhe length cooling mode, 2nd guess [m]
		this%L_new=maxval([L_h,L_c]) ! total bhe length [m]
		if (this%L_new==L_h) then
			this%mode='heating'
		else
			this%mode='cooling'
		end if

		tol=0.001d0
		iteration=0
		
		do while (abs(this%L_new-L_old)>tol*this%L_new .and. iteration<100)
			L_old=this%L_new
			R_cv_h_old=R_cv_h
			R_cv_c_old=R_cv_c
			R_bhe_h_old=R_bhe_h
			R_bhe_c_old=R_bhe_c
			call layout(this%N_bhe,this%Nx,this%Ny,this%N_odd,this%H,this%L_new,this%H_des)
			call classBHE(N1,N2,N3,N4,this%N_bhe,this%Nx,this%Ny,this%N_odd)
			Tp = penalty(N1,N2,N3,N4,this%N_bhe,Tp1)
			R_cv_h = this%convective_resistance(this%Reynolds(m_h,Th),Th)
			R_cv_c = this%convective_resistance(this%Reynolds(m_c,Tc),Tc)
			R_bhe_h=R_bhe_h_old-R_cv_h_old+R_cv_h
			R_bhe_c=R_bhe_c_old-R_cv_c_old+R_cv_c
			L_h=(this%Q_ga*R_10y+this%Q_ghm*R_1m+this%Q_ghd*(R_bhe_h+R_6h))/(this%Tg+273.15-((this%T_fih_design+this%T_foh_design)/2d0+273.15d0)-Tp) ! total bhe length heating mode [m]
			L_c=(this%Q_ga*R_10y+this%Q_gcm*R_1m+this%Q_gcd*(R_bhe_c+R_6h))/(this%Tg+273.15-((this%T_fic_design+this%T_foc_design)/2d0+273.15d0)-Tp) ! total bhe length cooling mode [m]
			this%L_new=maxval([L_h,L_c])
				if (this%L_new==L_h) then
					this%mode='heating'
				else
					this%mode='cooling'
				end if
			iteration=iteration+1
        end do
  
		! calculate XY coordinates 	
        if (allocated(this%bhe_XY_coordinates)) deallocate(this%bhe_XY_coordinates)
		allocate(this%bhe_XY_coordinates(this%N_bhe,2))
		allocate(Xs(this%Nx))
		allocate(Ys(this%Ny))
		allocate(YY(this%Nx*this%Ny))

		associate(coord=>this%bhe_XY_coordinates)
		coord=0d0

		Xs(1)=0d0
		do i=2,size(Xs)
			Xs(i)=Xs(i-1)+this%B_des    
		end do
	
		Ys(1)=0d0
		do i=2,size(Ys)
			Ys(i)=Ys(i-1)+this%B_des    
		end do      
		j=1
		do i=1,size(YY)
			YY(i)=Ys(j)
			j=j+1
			if(j>size(Ys)) j=1
		end do
	
		do i=1,this%N_bhe-this%N_odd
			coord(i,2)=YY(i)
		end do
	
		do i=1,this%Nx
			coord((i-1)*this%Ny+1:i*this%Ny,1)=Xs(i)
		end do
	
		if (this%N_odd>0) then
			coord(this%N_bhe-this%N_odd+1:this%N_bhe,1)=this%B_des*this%Nx
            do j=1,this%N_odd
				coord(this%N_bhe-this%N_odd+j,2)=this%B_des*(j-1)
            end do    
		end if	
		end associate
		
		deallocate(Xs)
		deallocate(Ys)
		deallocate(YY)
		deallocate(r_cc)
		deallocate(r_icc)
		deallocate(X)
		deallocate(I_temp)

    contains
    
            subroutine layout(N_bhe,Nx,Ny,N_odd,H,L,H_des)
                integer::N_bhe
                integer::Nx
                integer::Ny
                integer::N_odd
                double precision::H
                double precision::L
                double precision::H_des
                
                integer::j
                integer::count_factor
                integer::N_bhe_even
                integer::Nx1
                integer::NX2
                integer::Ny1
                integer::Ny2
                double precision::AR1
                double precision::AR2
                integer,dimension(:),allocatable::Nyy1
                integer,dimension(:),allocatable::Nyy2
                integer,dimension(:),allocatable::f
                integer,dimension(2)::Ntry
                double precision,dimension(2)::Htry
                double precision,dimension(2)::DH
                
                
				if (L<=H_des) then
				   H=L
				   N_bhe=1
				   Nx=1
				   Ny=1
				else
				   Ntry=[floor(L/H_des), ceiling(L/H_des)] ! rounding to next integers
				   Htry=L/Ntry
				   DH=abs(H_des-Htry)
				   if (Ntry(1)==Ntry(2) .or. DH(1)==DH(2)) then
					   j=1
				   else
					   j=1
					   if(abs(DH(2))<abs(DH(1))) j=2
				   end if
				   N_bhe=Ntry(j) ! choosing the nearest integer
				   H=Htry(j)
				   if (N_bhe==2 .or. N_bhe==3) then  ! line arrangement for N_bhe=2 or 3
					  Nx=N_bhe
					  Ny=1
				   else if (nint(N_bhe**0.5d0)-N_bhe**0.5d0==0) then ! squared arrangement 
						   Nx=N_bhe**0.5d0
						   Ny=Nx
				   else
					   
					   call factors_count(N_bhe,count_factor)
					   if (count_factor==0) then
						   count_factor=1
						   allocate(f(count_factor))
						   f(1)=1
					   else 
							allocate(f(count_factor))
							call factors(N_bhe,f,count_factor)
					   end if
				 
					   if (size(f)==1) then   ! N_bhe is prime
							   deallocate(f)
							   N_bhe_even=N_bhe-1 ! the nearest lower even number is defined
							   call factors_count(N_bhe_even,count_factor)
							   !allocate(f(count_factor))
							   if (count_factor==0) then
								   count_factor=1
								   allocate(f(count_factor))
								   f(1)=1
							   else 
								   allocate(f(count_factor))
								   call factors(N_bhe_even,f,count_factor)
							   end if
							   
					   end if
				   
					   if (nint(N_bhe_even**0.5d0)-N_bhe_even**0.5d0==0) then ! squared arrangement 
						   Nx=N_bhe**0.5d0
						   Ny=Nx
					   
					   else if (size(f)==2) then ! if there are only 2 prime factors
							  Nx=f(1)
							  Ny=f(2)
					   else if(size(f)>=2) then! if there are more than 2 prime factors, 2 layouts are created and compared
						   allocate(Nyy1(size(f)-1))
						   allocate(Nyy2(size(f)-1))
						   Nx1=f(1)*f(size(f)) ! Nx1 is obtained by multiplying the lowest and the highest factors
						   Nyy1(1)=1
						   do i=2,(size(f)-1)
							  Nyy1(i)=f(i)*Nyy1(i-1)
						   end do
						   Ny1=Nyy1(size(Nyy1))
						   AR1=maxval([Nx1,Ny1])/minval([Nx1,Ny1])
						   Nx2=f(1)*f(2) ! Nx2 is obtained by multiplying the first two factors
						   Nyy2(1)=1
						   do i=2,(size(f)-1)
							   Nyy2(i)=f(i+1)*Nyy2(i-1)
						   end do
						   Ny2=Nyy2(size(Nyy2))
						   AR2=maxval([Nx2,Ny2])/minval([Nx2,Ny2])
						   if (AR1<AR2) then
							   Nx=Nx1
							   Ny=Ny1
						   else
							   Nx=Nx2
							   Ny=Ny2
						   end if
						   deallocate(Nyy1)
						   deallocate(Nyy2)
					   end if
				   end if
				end if

				N_odd=N_bhe-Nx*Ny ! if N_bhe is prime, the configuration Nx X Ny misses 1 bhe, i.e. N_odd=1
				                    
                if (allocated(f)) deallocate(f)     
            end subroutine
            
            subroutine classBHE(N1,N2,N3,N4,N_bhe,Nx,Ny,N_odd)
                integer::N1
                integer::N2
                integer::N3
                integer::N4
                integer::N_bhe
                integer::Nx
                integer::Ny
                integer::N_odd
                integer::Nper
                if (N_odd<=1) then    
					N1=N_odd ! number of bhe with only 1 neighbor
					N2=4-N_odd ! number of bhe with 2 neighbors
					Nper=2*Nx+2*Ny-4 ! number of bhe on the perimeter
					N3=Nper-N2 ! number of bhe with 3 neighbors
					N4=N_bhe-N1-N2-N3 ! number of bhe with 4 neighbors
                else ! e.g. 54 boreholes, 7x7 + 5 on the last column
					N1=0 ! number of bhe with only 1 neighbor
					N2=5 ! number of bhe with 2 neighbors
					Nper=2*Nx+2*Ny+1-4 ! number of bhe on the perimeter
					N3=Nper-N2 ! number of bhe with 3 neighbors
					N4=N_bhe-N1-N2-N3 ! number of bhe with 4 neighbors
                end if    
            end subroutine
    
            function penalty(N1,N2,N3,N4,N_bhe,Tp1) result(Tp)
                integer::N1
                integer::N2
                integer::N3
                integer::N4
                integer::N_bhe
                double precision::Tp1
                double precision::Tp
                
                Tp=Tp1*(1*N4+0.5*N3+0.25*N2+0.1*N1)/N_bhe ! temperature penalty for the bhe field [K]
            end function penalty
                
            subroutine factors_count(N,count)
                integer::N
                integer::count   
                integer::i,j
                integer::N1
                logical::is_prime
				count=0
				N1=N
				i=2
				do while(i*i<=N1)
					do while(mod(N1,i)==0)
						count=count+1
						N1=N1/i
					end do
					i=i+1
				end do			
				if(N1>1) count=count+1    
            end subroutine
  
            subroutine factors(N,f,count)
            implicit none                
                integer::N
                integer::count
                integer,dimension(count)::f
                integer::i,j,k
                logical::is_prime  
				N1=N
				j=1
				i=2
				do while(i*i<=N1)
					do while(mod(N1,i)==0)
						f(j)=i
						N1=N1/i
						j=j+1
					end do
					i=i+1
				end do         
				if(N1>1) f(j)=N1 
            end subroutine        
    end subroutine	
    
    subroutine borehole_heat_exchanger_create_g_function(this,years)
        class(borehole_heat_exchanger),target::this
        integer::r_short_size
        integer::i,j,k
        integer::tstep
        integer::endtime
        integer::years
        double precision::temp
        double precision::h,n
        double precision::beta1
        double precision::alpha
        double precision,parameter::tol=1d-6
        
        double precision,dimension(this%N_bhe,this%N_bhe)::r
        double precision,dimension(this%N_bhe**2-this%N_bhe)::AA
        
        double precision,dimension(:),allocatable::A
        double precision,dimension(:),allocatable::r_short
        double precision,dimension(:),allocatable::N_short
        double precision,dimension(:),allocatable::beta
        double precision,dimension(:),allocatable::t
        double precision,dimension(:),allocatable::t_star
        double precision,dimension(:),allocatable::gamma1
        double precision,dimension(:),allocatable::g_interact
        double precision,dimension(:),allocatable::g_1
        double precision,dimension(:,:),allocatable::Da
        double precision,dimension(:,:),allocatable::Db
        double precision,dimension(:,:),allocatable::I1
        double precision,dimension(:,:),allocatable::I2
        double precision,dimension(:,:),allocatable::gg
        double precision,dimension(:,:),allocatable::W
         
        double precision,dimension(:),allocatable::Da1
        double precision,dimension(:),allocatable::Db1
        double precision,dimension(:),allocatable::I11
        double precision,dimension(:),allocatable::I21
        
        double precision::fun_gamma,fun_beta
        
        call simlog%output( 'Calculate the g-function for '//integer_to_text(years)//' years...')
        
		associate (coord=>this%bhe_XY_coordinates)
		alpha=this%lambda/this%C	
		r=0d0
		do i=1,this%N_bhe
			do j=1,this%N_bhe
				r(i,j)=(((coord(i,1)-coord(j,1))**2d0+(coord(i,2)-coord(j,2))**2d0))**0.5d0
			end do
        end do
        end associate
	
		allocate(A(this%N_bhe**2))
	
		j=1
		k=1
		do i=1,size(A)
			A(i)=r(j,k)
			j=j+1
			if(j>size(r,1)) then
				k=k+1
				j=1
			end if
		end do
	
		do i=1,size(A)-1
			do j=1,size(A)-i
				if(A(j)>A(j+1)) then
					temp=A(j)
					A(j)=A(j+1)
					A(j+1)=temp
				end if
			end do
		end do
	
		AA=A(this%N_bhe+1:this%N_bhe**2)
	
		h=1d0
		n=h
		r_short_size=1
		do k=2,size(AA,1)
			if (AA(k)==AA(k-1)) then
				n=n+1
				h=h+1
			else
				r_short_size=k-h+1
				n=1
			end if
		end do
	
		allocate(r_short(r_short_size))
		allocate(N_short(r_short_size))
		r_short(1)=AA(1)
		h=1d0
		n=h
		do k=2,size(AA,1)
			if (AA(k)==AA(k-1)) then
				n=n+1
				h=h+1
			else
				r_short(k-h+1)=AA(k)
				N_short(k-h)=n
				n=1
			end if
		end do
	
		N_short(size(r_short))=n
		beta=r_short/this%H

		tstep=3600*24
		endtime=years*365*24*3600

		allocate(t((endtime-tstep)/tstep+1))
		t(1)=tstep
		do i=2,size(t)
			t(i)=t(i-1)+tstep
		end do
	
		allocate(t_star(size(t)))
		allocate(gamma1(size(t)))
		this%g_function%t_s=this%H**2d0/(9d0*alpha)
		t_star=t/this%g_function%t_s
		gamma1=3d0/(2d0*(t_star**0.5d0))
		
	
		allocate(Da(size(gamma1),size(beta)))
		allocate(Db(size(gamma1),size(beta)))
		allocate(I1(size(gamma1),size(beta)))
		allocate(I2(size(gamma1),size(beta)))
		allocate(gg(size(gamma1),size(beta)))
		allocate(W(size(gamma1),size(beta)))
	
		do i=1,size(beta)
			do j=1,size(gamma1)
				Da(j,i)=(beta(i)**2d0+1d0)**0.5d0*erfc(gamma1(j)*(beta(i)**2d0+1d0)**0.5d0)-beta(i)*erfc(gamma1(j)*beta(i))-(exp(-gamma1(j)*2*(beta(i)**2+1d0))-exp(-gamma1(j)**2*beta(i)**2))/(gamma1(j)*pi**0.5d0) 
				Db(j,i)=(beta(i)**2d0+1d0)**0.5d0*erfc(gamma1(j)*(beta(i)**2d0+1d0)**0.5d0)-0.5d0*(beta(i)*erfc(gamma1(j)*beta(i))+(beta(i)**2d0+4d0)**0.5d0*erfc(gamma1(j)*(beta(i)**2d0+4d0)**0.5d0))-(exp(-gamma1(j)**2d0*(beta(i)**2d0+1d0))-0.5d0*(exp(-gamma1(j)**2d0*beta(i)**2d0)+exp(-gamma1(j)**2d0*(beta(i)**2d0+4d0))))/(gamma1(j)*pi**0.5d0)                    
				fun_gamma=gamma1(j)
                fun_beta=beta(i)
                call Integration_adaptive_simpson(fun1,beta(i),(beta(i)**2d0+1d0)**0.5d0, tol, I1(j,i))
				call Integration_adaptive_simpson(fun1,(beta(i)**2d0+1d0)**0.5d0,(beta(i)**2d0+4d0)**0.5d0, tol, I2(j,i))         
			end do
		end do
	
		gg=I1-Da-I2-Db
		do i=1,size(r_short)
			W(:,i)=gg(:,i)*N_short(i)
		end do
	
		allocate(g_interact(size(gamma1)))
		do i=1,size(g_interact)
			g_interact(i)=sum(W(i,:))/this%N_bhe
		end do
	
		allocate(Da1(size(gamma1)))
		allocate(Db1(size(gamma1)))
		allocate(I11(size(gamma1)))
		allocate(I21(size(gamma1)))
	
		beta1=this%D_bhe/2d0/this%H    
		do j=1,size(gamma1)
			 Da1(j)=(beta1**2d0+1d0)**0.5d0*erfc(gamma1(j)*(beta1**2d0+1d0)**0.5d0)-beta1*erfc(gamma1(j)*beta1)-(exp(-gamma1(j)**2d0*(beta1**2d0+1d0))-exp(-gamma1(j)**2d0*beta1**2d0))/(gamma1(j)*pi**0.5d0)
			 Db1(j)=(beta1**2d0+1d0)**0.5d0*erfc(gamma1(j)*(beta1**2d0+1d0)**0.5d0)-0.5d0*(beta1*erfc(gamma1(j)*beta1)+(beta1**2d0+4d0)**0.5d0*erfc(gamma1(j)*sqrt(beta1**2d0+4d0)))-(exp(-gamma1(j)**2d0*(beta1**2d0+1d0))-0.5d0*(exp(-gamma1(j)**2d0*beta1**2d0)+exp(-gamma1(j)**2d0*(beta1**2d0+4d0))))/(gamma1(j)*pi**0.5d0)     
				fun_gamma=gamma1(j)
                fun_beta=beta1
			 call Integration_adaptive_simpson(fun1,beta1, (beta1**2d0+1d0)**0.5d0, tol, I11(j))
			 call Integration_adaptive_simpson(fun1,(beta1**2d0+1d0)**0.5d0, (beta1**2d0+4d0)**0.5d0, tol, I21(j))
		end do
	
		allocate(g_1(size(gamma1)))
		g_1=I11-Da1-I21-Db1
		
		call this%g_function%load(t_star,g_1+g_interact)
			
		deallocate(A)
		deallocate(r_short)
		deallocate(N_short)
		deallocate(beta)
		deallocate(t)
		deallocate(t_star)
		deallocate(gamma1)
		deallocate(g_interact)
		deallocate(g_1)
		deallocate(Da)
		deallocate(Db)
		deallocate(I1)
		deallocate(I2)
		deallocate(gg)
		deallocate(W)
		deallocate(Da1)
		deallocate(Db1)
		deallocate(I11)
		deallocate(I21)
			 
	contains
	
		function fun1(x) result(y)
			double precision :: y
			double precision, intent(in) :: x
			if(abs(x**2d0-fun_beta**2d0)/=0d0) then
				y=erfc(fun_gamma*x)/(x**2d0-fun_beta**2d0)**0.5d0
			else
				y=erfc(fun_gamma*(x+1d-10))/((x+1d-10)**2d0-fun_beta**2d0)**0.5d0
			end if
		end function fun1
          
    end subroutine
	
end module