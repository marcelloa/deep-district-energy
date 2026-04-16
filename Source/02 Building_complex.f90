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
module building_complex_module
    use constants_module
    use datalog
    use meteolib
    use maths
    use utility
    implicit none
        
    enum, bind(c)
        enumerator::WALL=1,SLAB,SLAB_ROOF
    end enum    
    
    enum, bind(c)
        enumerator::GROUND=1,UNHEATED,HEATED,ROOF
    end enum    

    enum, bind(c)
        enumerator::S_PEOPLE_SEN=1,S_PEOPLE_LAT,S_APPLIANCES,S_VENTILATION,S_TH_SET,S_TC_SET,S_DHW,S_INTERNAL_GAIN,S_POWER_AVAILABLE
    end enum    
	
	type point
        double precision::x
        double precision::y
        double precision::z
    contains
        procedure::equal=>point_equal ! test if X,Y,Z of this point and point in argument are equal within a given tolerance
        procedure::translate=>point_translate ! return a new point translated by a given array of DX,DY,DZ
    end type  
    
    type material_mass
        character(:),allocatable::name ! name (code) of material
        double precision::mass=0d0 ! mass in kg
    end type    
    
    type material_list_element
        type(material_mass)::value
        type(material_list_element),pointer::next=>null()
    end type    

    type material_list
        type(material_list_element),pointer::first=>null()
    contains
        procedure::update=>material_list_update
        procedure::initialize=>material_list_initialize
    end type    

    type material_properties
        double precision::lambda !
        double precision::rho ! 
        double precision::cp !
        double precision::vrf !
        double precision::resistance ! (m2,K)/W
    end type    

    type material_type
        character(:),allocatable::name
        character(:),allocatable::description
        type(material_properties)::properties
    end type    
           
	type layer
		double precision::thickness
        character(:),allocatable::material_name        
		type(material_properties)::properties
    end type

    type window_structure
        character(:),allocatable::name
        character(:),allocatable::description
        double precision::Ug    ! U value for glass, W/m2,K
        double precision::gn    ! g value for glass, normal direction, -
        double precision::Ff    ! frame fraction, e.g 0.2 = 20% is frame
        double precision::Uf    ! U value for the frame , W/m2,K
        double precision::Htb   ! thermal bridge, W/m2,K 
        double precision::Fw    ! solar transmittance angle modifier (mean value), -
        double precision::internal_resistance ! (m2,K)/W
        double precision::external_resistance ! W/(m2,K)
        double precision::gw    ! mean g value of window, releated to window gross area, -
        double precision::resistance ! (m2,K)/W, not including external and internal resistances 
        double precision::transmittance ! W/(m2,K), already including external and internal resistances
    contains
        procedure::calculate=>window_structure_calculate
    end type

    type window
		double precision::area      ! surface, m2
		double precision::azimuth   ! azimuth angle, rad
		double precision::tilt      ! tilt angle, rad
		type(window_structure),pointer::struct
    end type
     
    type wall_slab_structure
        character(:),allocatable::name
        character(:),allocatable::ISOclass
        double precision::internal_resistance ! (m2,K)/W
        type(layer),dimension(:),allocatable::layers ! 1=internal, N=external
        double precision::external_resistance ! (m2,K)/W
        double precision::resistance ! (m2,K)/W, not including external and internal resistances
        double precision::transmittance ! W/(m2,K), already including external and internal resistances
        double precision::specific_capacity ! J/(m2,K)
        double precision,dimension(5)::ISO_C_coeff
        double precision,dimension(4)::ISO_R_coeff
    contains
        procedure::calculate=>wall_slab_structure_calculate
    end type
    
	type wall_slab
		integer::type ! see wall_slab_type
		double precision::area
		double precision::area_above_ground
		double precision::area_below_ground
		double precision::azimuth
		double precision::tilt
		type(wall_slab_structure),pointer::struct
        logical::is_adiabatic
        double precision::Htb   ! thermal bridge, W/K
        double precision::resistance ! (m2,K)/W, not including external and internal resistances, including thermal bridge effect
        double precision::transmittance ! W/(m2,K), already including external, internal resistances and thermal bridge effect       
        type(point),dimension(:),allocatable::vertices ! the points representing vertices of the external gross area (including windows) above ground
        type(point),dimension(:),allocatable::grid ! grid of points in which the external gross area is discretized
        double precision,dimension(:,:,:),allocatable::fsol  ! solar factor 0 (the point is in shadow) or 1 (the point can see the sun), size(grid) x 12 x 24
    contains
        procedure::make_horizontal_grid=>wall_slab_make_horizontal_grid
        procedure::make_vertical_grid=>wall_slab_make_vertical_grid
        procedure::calculate_materials=>wall_slab_calculate_materials
    end type
    
    type dhw_system
        double precision::V_storage     ! m3 Volume of the tank
        double precision::UA_storage    ! W/K UA value of the tank
        double precision::Q_heater      ! W, heating capacity requested by the DHW tank
        double precision,dimension(8760)::Q_dhw         ! Wh, dhw load (instantaneous demand)
        double precision,dimension(8760)::mdotcp_dhw    ! Wh/K, dhw heat capacity rate (instantaneous demand)
        double precision,dimension(8760)::T_main        ! °C Temperature of the main
    end type    

    type hvac_system
        double precision::T_supply_hw ! °C supply temperature of hot water loop
        double precision::T_supply_cw ! °C supply temperature of chilled water loop
        double precision,dimension(8760)::Q_heat_hyd ! Wh, heat supplied to hydronic system (fan-coil, radiant floor)    
        double precision,dimension(8760)::Q_cool_hyd ! Wh, cool supplied to hydronic system (fan-coil, radiant floor)    
        double precision,dimension(8760)::Q_heat_AHU ! Wh, heat supplied to AHU heating coils    
        double precision,dimension(8760)::Q_cool_AHU ! Wh, cool supplied to AHU cooling coil
        double precision,dimension(8760)::W_par_hyd  ! Wh, parasitic consumption for hydronic system, heat & cool combined
        double precision,dimension(8760)::W_par_AHU  ! Wh, parasitic consumption for AHU, heat & cool combined
    end type
        
	type building_thermal_zone
		double precision::air_mass ! kg
		double precision::internal_capacity ! J/K, includes forniture 
		type(wall_slab),dimension(:),allocatable::external_structures ! walls, floor, roof
		type(wall_slab),dimension(:),allocatable::internal_structures ! interfloor slabs and partition walls (not used now)
		type(window),dimension(:),allocatable::windows
        double precision::envelope_area ! external heat transfer area, excluding adiabatic walls, m2
        double precision::radiant_area ! external heat transfer area, including adiabatic walls, m2
        double precision::envelope_HT ! transmittance-area product, W/K
        double precision::Htb ! thermal bridge, W/K
        double precision,dimension(8760,9)::schedules
        double precision,dimension(8760)::T_air
        double precision,dimension(8760)::X_air
        double precision,dimension(8760)::T_op
        double precision,dimension(8760)::Q_sh_sen ! Wh, space heating sensible load
        double precision,dimension(8760)::Q_sc_sen ! Wh, space cooling sensible load
        double precision,dimension(8760)::Q_sh_lat ! Wh, space heating latent load
        double precision,dimension(8760)::Q_sc_lat ! Wh, space cooling latent load
        double precision,dimension(8760)::W_app ! Wh, plug loads
    end type
    
    ! edge of the buidling footprint
    type footprint_edge
        double precision::X1,Y1,X2,Y2   ! 1 (start) and 2 (end) X,Y coordinate
        double precision::azimuth       ! azimuth for the side of the edge facing the outdoor air
        double precision::length        ! length
    contains
        procedure::update=>footprint_edge_update ! calculate length and azimuth
    end type    
    
    type building_unit_geometry
        double precision::footprint_area    ! gross area of the building (unit) footprint, net of correction factor footprint to gross 
        double precision::footprint_azimuth ! the azimuth of segment having miniumum absolute azimuth (relevant for solar device orientation)
        double precision::floor_area_gross  ! gross floor area, including non heated indoor and outdoor spaces (e.g. stairwell, balconies) 
        double precision::floor_area_heated ! floor area heated, including external and internal walls
        double precision::floor_area_useful ! floor area heated, excluding external and internal walls
        double precision::volume_gross      ! gross volume, including non heated indoor and outdoor spaces (e.g. stairwell, balconies), walls, slabs, roof
        double precision::volume_useful     ! useful volume, excluding external, internal walls, slabs and roof
        double precision::window_area       ! window area
        integer::number_of_storeys_above_ground ! in case the building unit is partly above and partly below ground
        integer::number_of_storeys_below_ground ! in case the building unit is partly above and partly below ground
        integer::number_of_storeys          ! total storeys above and below ground
        double precision::height_above_ground ! in case the building unit is partly above and partly below ground
        double precision::height_below_ground ! in case the building unit is partly above and partly below ground
        double precision::height            ! total height above and below ground
        double precision::height_floor_start! height from ground of the starting floor, negative below ground, positive above ground, m
    end type

    type building_unit
        integer::building_id
        type(building),dimension(:),pointer::p_buildings
        type(building_complex),pointer::p_building_complex
	    character(:),allocatable::name
        type(material_list)::materials_inventory ! the list of materials names and masses for LCA
	    double precision::footprint_to_gross_area ! related to a generic floor of the 2.5D unit
	    double precision::gross_to_heated_area ! related to the whole unit
	    double precision::heated_to_useful_area ! related to the whole unit
	    double precision::useful_floor_height ! (for 1 floor) the fraction of floor height that is useful (occupied by the air of the thermal zone)
        double precision::window_to_useful_area ! unit window area = window_to_useful_area * unit useful_area * number of floors above ground / number of floors
        double precision,dimension(:),allocatable::window_to_segment ! window area on side j = total window area * segment(j)*window_to_segment(j)/sum_j=1:N(segment(j)*window_to_segment(j)), in this way 0 window area on one side can be managed, e.g. [1,1,1,0] indicates 0 window area on segment 4.
        integer,dimension(:),allocatable::wall_is_adiabatic ! flag 1 indicates adiabatic external wall (i.e., separation between adjoning building units)
	    integer::floor_start_type ! type of space in communication with start floor: ground, unheated, heated, roof
	    integer::floor_end_type   ! type of space in communication with end floor: ground, unheated, heated, roof     
        integer::floor_start      ! the starting floor w.r.t. the building storeys, above ground storeys: +1,+2,.. / below ground storeys: -1,-2,..
	    integer::floor_end        ! the ending floor w.r.t. the building storeys,above ground storeys: +1,+2,.. / below ground storeys: -1,-2,..  
        type(building_unit_geometry)::geometry
        type(building_thermal_zone)::thermal_zone    
        type(dhw_system)::dhw
        type(hvac_system)::hvac
	    character(:),allocatable::name_slab_basement
	    character(:),allocatable::name_slab_interfloor
	    character(:),allocatable::name_roof
	    character(:),allocatable::name_walls
	    character(:),allocatable::name_windows       
        double precision::linear_thermal_bridge ! W/(m,K), thermal bridge associated to the base of external walls in each floor
	    double precision::min_RH
	    double precision::max_RH
        character(:),allocatable::schedule_ventilation
        character(:),allocatable::schedule_occupancy
        character(:),allocatable::schedule_appliances
        character(:),allocatable::schedule_dhw
        character(:),allocatable::schedule_heating_setpoint
        character(:),allocatable::schedule_cooling_setpoint
        ! optional fields
        character(:),allocatable::schedule_power_available !(-), 'default' = 1, 0 in case of electric outage (lights turned off, mechanical ventilation off, heating/cooling off)
        character(:),allocatable::schedule_internal_gain !  (W), 'default' = 0, used to include the effect of passive solutions not modeled by EN 52016 
        double precision::Q_h_specific !(W/m2), default = 50
        double precision::Q_c_specific !(W/m2), default = 80
        ! end optional fields
	    double precision::ventilation_rate
	    double precision::crowd_density
	    double precision::power_density
	    double precision::dhw_use
	    double precision::blinds_G
	    double precision::blinds_b
	    integer::hvac_type
	    integer::dhw_type
	    double precision::emission_efficiency
	    double precision::control_efficiency
	    double precision::distribution_efficiency
	    double precision::distribution_volumetric_efficiency
	    double precision::heat_recovery_effectiveness
	    double precision::maximum_airvolume_changes
	    double precision::average_flat_area
	    double precision::dhw_peak_flow
	    double precision::dhw_peak_duration
	    double precision::dhw_preheating
	    double precision::water_main_temperature
	    double precision::dhw_storage_temperature	
        double precision,dimension(:),allocatable::T_old
        double precision,dimension(:),allocatable::T_new
        double precision,dimension(:,:),allocatable::A,invAww
        double precision,dimension(:),allocatable::B   
    contains
        procedure::create_thermal_zone=>building_unit_create_thermal_zone
        procedure::create_dhw=>building_unit_create_dhw
        procedure::calculate_T_op=>building_unit_calculate_T_op
        procedure::calculate_AB=>building_unit_calculate_AB       
        procedure::calculate_hvac=>building_unit_calculate_hvac
    end type	
    
    ! 2.5D building and its associated building units.
    ! A complex building can be split in a number of 2.5D buildings.
    ! Building Units can be used to split a building vertically. 
    ! Energy needs are calculated at building unit level.
    ! Building Units belonging to different (adjoining) buildings can be merged if the resulting building unit is a 2.5D geometry:
    ! the merged building unit is the single entity used for energy calculations. 
    type building
        integer::id
        character(:),allocatable::name
        double precision,dimension(:),allocatable::X ! X coordinates of the X,Y sequence describing the footprint polygon, ordered in clock-wise direction
        double precision,dimension(:),allocatable::Y ! Y coordinates of the X,Y sequence describing the footprint polygon, ordered in clock-wise direction
        type(footprint_edge),dimension(:),allocatable:: footprint_segments ! segments describing the footprint polygon, ordered in clock-wise direction
		integer::floors_ag  ! number of storeys above ground
		integer::floors_bg  ! number of storeys below ground
		double precision::height_ag  ! height above ground   
		double precision::height_bg  ! height below ground      
        type(building_unit),dimension(:),allocatable::units
        !type(substation),dimension(:),allocatable::substations
    contains
        procedure::calculate_geometry=>building_calculate_geometry
        procedure::create_thermal_zones=>building_create_thermal_zones
        procedure::create_dhw_systems=>building_create_dhw_systems
        procedure::calculate_energy_needs=>building_calculate_energy_needs 
    end type
    
    type position_deg
        double precision::latitude  ! deg
        double precision::longitude ! deg
    end type
    
    type schedule
        character(:),allocatable::name
        character(:),allocatable::timestep
        integer::size
        double precision,dimension(:),allocatable::values        
    end type
        	
	type building_complex
        type(meteo)::meteodata
        character(:),allocatable::name
        type(position_deg)::position
        type(material_type),dimension(:),allocatable::materials
        type(wall_slab_structure),dimension(:),pointer::wall_slab_structures=>null()
        type(window_structure),dimension(:),pointer::window_structures=>null()
        type(building),dimension(:),allocatable::buildings
        type(building),dimension(:),allocatable::neighbors
        type(schedule),dimension(:),pointer::schedules=>null()
     contains
        procedure::get_wall_slab=>building_complex_get_wall_slab
        procedure::get_window=>building_complex_get_window
        procedure::get_schedule=>building_complex_get_schedule
        procedure::calculate_fsol=>building_complex_calculate_fsol
        procedure::calculate_energy_needs=>building_complex_calculate_energy_needs
        !procedure::calculate_substations=>calculate_substation_system
    end type    
       
    contains
	
	function point_translate(this,shift) result(new_point)
	    class(point)::this
	    type(point)::new_point
	    double precision,dimension(3)::shift
	    new_point%X=this%X+shift(1)
	    new_point%Y=this%Y+shift(2)
	    new_point%Z=this%Z+shift(3)
    end function

    function point_equal(this,test_point,eps) result(test)
        class(point)::this
        type(point)::test_point
        double precision::eps
        logical::test
        if ((abs(this%X-test_point%X)<=eps).and.(abs(this%Y-test_point%Y)<=eps).and.(abs(this%Z-test_point%Z)<=eps)) then
            test=.true.
        else
            test=.false.
        end if    
    end function
            
    subroutine material_list_initialize(this)
        class(material_list)::this
        type(material_list_element),pointer::material_first,material_new
        if (.not.associated(this%first)) then
            return
        else    
            do ! loop on material list elements until next element is null
                material_first=>this%first
                if (.not.associated(material_first%next)) then
                    this%first=>null()
                    deallocate(material_first)
                    exit
                else
                    this%first=>material_first%next
                    deallocate(material_first)
                end if    
            end do
        end if
    end subroutine    

    ! serach in the list of materials the given material name and increment its mass
    subroutine material_list_update(this,material_name,material_mass)
        class(material_list)::this
        type(material_list_element),pointer::material_start,material_new
        logical::found
        character(*)::material_name
        double precision::material_mass
        if (.not.associated(this%first)) then
            allocate(material_new)
            material_new%value%name=material_name
            material_new%value%mass=material_mass
            this%first=>material_new
        else    
            material_start=>this%first
            found=.false.
            do ! loop on material list elements until next element is null
                if (material_start%value%name==material_name) then
                    material_start%value%mass=material_start%value%mass+material_mass
                    found=.true.
                end if    
                if (.not.associated(material_start%next)) then
                    exit
                else
                    material_start=>material_start%next
                end if    
            end do
            if (.not.found) then
                allocate(material_new)
                material_new%value%name=material_name
                material_new%value%mass=material_mass
                material_start%next=>material_new
            end if    
        end if
    end subroutine    
    
    ! provides the intersection between a segment and a segment
    ! segment A from (X1,Y1) to (X2,Y2)
    ! segment B from (X3,Y3) to (X4,Y4)
    ! X,Y are the coordinates of the intersection point
    ! stat=1 if X,Y are found
    ! stat=0 if X,Y are not found (parallel lines or (X,Y) not in segment A and B) 
    subroutine intersect_ss(X1,Y1,X2,Y2,X3,Y3,X4,Y4,X,Y,stat)
        double precision::X1,Y1,X2,Y2,X3,Y3,X4,Y4,X,Y,aA,bA,cA,aB,bB,cB
        integer::stat
        ! consider two lines in the form a*x+b*y=c
        ! a) x2<>x1
        ! y=m*(x-x1)+y1
        ! m=(y2-y1)/(x2-x1)
        ! (x2-x1)*y=(y2-y1)*(x-x1)+y1*(x2-x1)
        ! (y1-y2)*x +(x2-x1)*y=(y1-y2)*x1+y1*(x2-x1)
        ! a=(y1-y2)
        ! b=(x2-x1)
        ! c=(y1-y2)*x1+y1*(x2-x1)
        ! b) y2<>y1
        ! x=m*(y-y1)+x1
        ! m=(x2-x1)/(y2-y1)
        ! (y2-y1)*x=(x2-x1)*(y-y1)+x1*(y2-y1)
        ! (x1-x2)*y +(y2-y1)*x=(x1-x2)*y1+x1*(y2-y1)
        ! a=(y1-y2)
        ! b=(x2-x1)
        ! c=(x2-x1)*y1+x1*(y1-y2)
        aA=Y1-Y2
        bA=X2-X1
        cA=(Y1-Y2)*X1+Y1*(X2-X1)
        aB=Y3-Y4
        bB=X4-X3
        cB=(Y3-Y4)*X3+Y3*(X4-X3)
        if ((aA*bB-aB*bA).ne.0) then
            X=( bB*cA-bA*cB)/(aA*bB-aB*bA)
            Y=(-aB*cA+aA*cB)/(aA*bB-aB*bA)
            if (X>=min(X3,X4) .and. X<=max(X3,X4) .and. Y>=min(Y3,Y4) .and. Y<=max(Y3,Y4) .and. &
               X>=min(X1,X2) .and. X<=max(X1,X2) .and. Y>=min(Y1,Y2) .and. Y<=max(Y1,Y2)) then     
                stat=1
                return
            end if
        end if
        X=0d0
        Y=0d0
        stat=0
    end subroutine   
    
    ! return the point along segment (X1,Y1)-(X2,Y2) at dimensional distance
    ! d from (X1,Y1)
    subroutine segment_pd(X1,Y1,X2,Y2,d,X,Y)
        double precision::X1,Y1,X2,Y2,d,X,Y,l,s
        l=((X2-X1)**2+(Y2-Y1)**2)**0.5
        s=d/l
        X=X1+(X2-X1)*s
        Y=Y1+(Y2-Y1)*s
    end    

    ! return the point along segment (X1,Y1)-(X2,Y2) at nondimensional distance
    ! s from (X1,Y1)
    subroutine segment_ps(X1,Y1,X2,Y2,s,X,Y)
        double precision::X1,Y1,X2,Y2,s,X,Y
        X=X1+(X2-X1)*s
        Y=Y1+(Y2-Y1)*s
    end subroutine
    
    function decode_floor_type(name) result(i)
        integer::i
        character(*)::name
        i=0;
        if (name=="ground") i=GROUND
        if (name=="unheated") i=UNHEATED
        if (name=="heated") i=HEATED
        if (name=="roof") i=ROOF
        if (i==0) then
            call simlog%output('Error: floor type value '//name//' not valid.')
            stop
        end if    
    end function
    
    function int_to_logical(i) result(l)
        integer::i
        logical::l
        if (i<=0) then
            l=.false.
        else
            l=.true.
        end if    
    end function

    function building_complex_get_schedule(this,name) result(sched)
        class(building_complex)::this
        type(schedule),pointer::sched
        character(*)::name
        integer::i
        sched=>null()
        do i=1,size(this%schedules)
            if (this%schedules(i)%name==name) then
                sched=>this%schedules(i)
                exit
            end if    
        end do    
        if (.not.associated(sched)) then
            call simlog%output('Error: schedule name '//name//' not a valid schedule.')
            stop
        end if    
    end function
        
    function building_complex_get_wall_slab(this,name) result(struct)
        class(building_complex)::this
        type(wall_slab_structure),pointer::struct
        character(*)::name
        integer::i
        struct=>null()
        do i=1,size(this%wall_slab_structures)
            if (this%wall_slab_structures(i)%name==name) then
                struct=>this%wall_slab_structures(i)
                exit
            end if    
        end do    
        if (.not.associated(struct)) then
            call simlog%output('Error: structure name '//name//' not a valid slab.')
            stop
        end if    
    end function
    
    function building_complex_get_window(this,name) result(struct)
        class(building_complex)::this
        type(window_structure),pointer::struct
        character(*)::name
        integer::i
        struct=>null()
        do i=1,size(this%window_structures)
            if (this%window_structures(i)%name==name) then
                struct=>this%window_structures(i)
                exit
            end if    
        end do    
        if (.not.associated(struct)) then
            call simlog%output('Error: structure name '//name//' not a valid window.')
            stop
        end if    
    end function
    
    subroutine wall_slab_structure_calculate(this)
        class(wall_slab_structure)::this
        character(len=2)::ISOclass
        integer::i
        this%resistance=0d0
        do i=1,size(this%layers)
            this%resistance=this%resistance+this%layers(i)%properties%resistance
        end do    
        this%transmittance=1d0/(this%internal_resistance+this%resistance+this%external_resistance)
        this%specific_capacity=0d0
        do i=1,size(this%layers)
            this%specific_capacity=this%specific_capacity+this%layers(i)%properties%rho*this%layers(i)%properties%cp*this%layers(i)%thickness
        end do    
        this%ISO_R_coeff=(/1d0/6,1d0/3,1d0/3,1d0/6/)
        ISOclass=this%ISOclass
        select case (ISOclass)
            case ('D')
                this%ISO_C_coeff=(/1d0/8,1d0/4,1d0/4,1d0/4,1d0/8/)
            case ('I')
                this%ISO_C_coeff=(/0d0,0d0,0d0,0d0,1d0/)
            case ('E')
                this%ISO_C_coeff=(/1d0,0d0,0d0,0d0,0d0/)
            case ('IE')
                this%ISO_C_coeff=(/1d0/2,0d0,0d0,0d0,1d0/2/)
            case ('M')
                this%ISO_C_coeff=(/0d0,0d0,1d0,0d0,0d0/)
            case default
                call simlog%output('Error: ISO class '//this%ISOclass//' not found.')   
                stop
        end select
         
    end subroutine

    subroutine window_structure_calculate(this)
        class(window_structure)::this
        this%gw=this%gn*(1-this%Ff)*this%Fw
        this%resistance=1d0/(this%Ug*(1-this%Ff)+this%Uf*(this%Ff)+this%Htb)
        this%transmittance=1d0/(this%internal_resistance+this%resistance+this%external_resistance)
    end subroutine
    
    
    subroutine footprint_edge_update(this)
        class(footprint_edge)::this
        double precision::DX,DY,az,tmp
        ! calculate length
        this%length=((this%X1-this%X2)**2+(this%Y1-this%Y2)**2)**0.5
        ! calculate azimuth
        DX=this%X2-this%X1
        DY=this%Y2-this%Y1
        if (DX.ne.0) then
            tmp=atan(DY/DX)
            if (DY>=0 .and. DX>0) then
                az=1.5*pi-tmp-pi/2
            else if (DY<=0 .and. DX<0) then 
                az=-0.5*pi-tmp+pi/2
            else if (DY>=0 .and. DX<0) then
                az=-tmp
            else ! DY<=0 .and. DX>0
                az=-pi-tmp
            end if
        else ! DX==0
            if (DY>0) then
                az=pi/2
            else
                az=-pi/2
            end if
        end if  
        this%azimuth=az
    end subroutine
    
    subroutine wall_slab_make_horizontal_grid(this,ds)
        class(wall_slab)::this
        double precision::ds
        integer::num_points,num_vertices,i,ix,iy,id,NX,NY,count,stat,s
        double precision::xMin,xMax,yMin,yMax,Xa,Ya,Xb,Yb,xInt,yInt,xTmp,yTmp,Xold,Yold,Xnew,Ynew,X1,Y1,X2,Y2  
        integer,dimension(:),allocatable::idx
        num_vertices=size(this%vertices)
        xMin=1d12
        yMin=1d12
        xMax=-1d12
        yMax=-1d12
        do i=1,num_vertices
            xMin=min(xMin,this%vertices(i)%x)
            xMax=max(xMax,this%vertices(i)%x)
            yMin=min(yMin,this%vertices(i)%y)
            yMax=max(yMax,this%vertices(i)%y)            
        end do    
        ! internal point near first segment (0.5 m distance)
        call segment_pd(this%vertices(1)%x,this%vertices(1)%y,this%vertices(2)%x,this%vertices(2)%y,0.5d0,Xa,Ya);
        call segment_pd(this%vertices(1)%x,this%vertices(1)%y,this%vertices(num_vertices)%x,this%vertices(num_vertices)%y,0.5d0,Xb,Yb);
        xInt=0.5d0*(Xa+Xb)
        yInt=0.5d0*(Ya+Yb)
        NX=nint((xMax-xMin)/ds)
        NY=nint((yMax-yMin)/ds)
        allocate(idx(NX*NY))
        id=0
        do ix=1,NX
            do iy=1,NY
                id=id+1
                xTmp=xMin+1+(ix-1)*ds
                yTmp=yMin+1+(iy-1)*ds
                ! test intersections between 'internal point to point' line and each segment of the polygon
                count=0
                Xold=-2*xMin
                Yold=-2*yMin
                do s=1,num_vertices
                    X1=this%vertices(s)%x
                    Y1=this%vertices(s)%y
                    if (s<num_vertices) then
                        X2=this%vertices(s+1)%x
                        Y2=this%vertices(s+1)%y
                    else
                        X2=this%vertices(1)%x
                        Y2=this%vertices(1)%y
                    end if                        
                    call intersect_ss(xInt,yInt,xTmp,yTmp, &
                        X1,Y1,X2,Y2,Xnew,Ynew,stat)
                    if ((stat==1) .and. (Xnew.ne.Xold) .and. (Ynew.ne.Yold)) then ! avoid double counting of vertices
                        count=count+1
                        Xold=Xnew
                        Yold=Ynew
                    end if  
                end do
                if (mod(count,2)==0) then
                    idx(id)=1
                else
                    idx(id)=0
                end if
            end do
        end do 
        num_points=sum(idx)
        allocate(this%grid(num_points))
        id=0
        i=0
        do ix=1,NX
            do iy=1,NY
                id=id+1
                if (idx(id)==1) then
                    i=i+1
                    xTmp=xMin+ds/2+(ix-1)*ds
                    yTmp=yMin+ds/2+(iy-1)*ds
                    this%grid(i)=point(xTmp,yTmp,this%vertices(1)%z)
                    !print *,this%grid(i)%x,this%grid(i)%y,this%grid(i)%z
                end if    
            end do
        end do    
        ! allocate fsol
        allocate(this%fsol(num_points,12,24))
    end subroutine
    
    subroutine wall_slab_make_vertical_grid(this,ds)
        class(wall_slab)::this
        double precision::ds
        integer::num_points,is,iz,id,NS,NZ
        double precision::height,d,X,Y
        id=0
        height=this%vertices(3)%z-this%vertices(1)%z
        NZ=floor(height/ds)
        d=((this%vertices(1)%x-this%vertices(2)%x)**2+(this%vertices(1)%y-this%vertices(2)%y)**2)**0.5
        NS=nint(d/ds) 
        num_points=NS*NZ
        allocate(this%grid(num_points))
        do is=1,NS
            call segment_ps(this%vertices(1)%x,this%vertices(1)%y,this%vertices(2)%x,this%vertices(2)%y,(is-0.5d0)/NS,X,Y)
            do iz=NZ,1,-1
                id=id+1
                this%grid(id)=point(X,Y,this%vertices(1)%z+ds/2+(iz-1)*ds)
                !print *,this%grid(id)%x,this%grid(id)%y,this%grid(id)%z
            end do   
        end do
        ! allocate fsol
        allocate(this%fsol(num_points,12,24))        
    end subroutine    

    ! for each layer, calculate the mass and update the table of materials
    subroutine wall_slab_calculate_materials(this,materials_inventory)
        class(wall_slab)::this
        type(material_list)::materials_inventory
        integer::i
        double precision::mass
        do i=1,size(this%struct%layers)
            mass=this%area*this%struct%layers(i)%thickness*this%struct%layers(i)%properties%rho
            call materials_inventory%update(this%struct%layers(i)%material_name,mass)    
        end do    
    end subroutine    
    
    subroutine building_complex_calculate_energy_needs(this)
        class(building_complex)::this
        type(schedule),pointer::p_schedule
        integer::b,bu,h,hd
        double precision::people
        ! calculate schedules of building unit thermal zone
        do b=1,size(this%buildings)
            do bu=1,size(this%buildings(b)%units)
                ! occupancy, sensible and latent
                people=this%buildings(b)%units(bu)%geometry%floor_area_useful*this%buildings(b)%units(bu)%crowd_density
                p_schedule=>this%get_schedule(this%buildings(b)%units(bu)%schedule_occupancy)
                this%buildings(b)%units(bu)%thermal_zone%schedules(1:8760,S_PEOPLE_SEN)=p_schedule%values(1:8760)*people*80d0 ! W
                this%buildings(b)%units(bu)%thermal_zone%schedules(1:8760,S_PEOPLE_LAT)=p_schedule%values(1:8760)*people*38d0 ! W
                ! appliances
                p_schedule=>this%get_schedule(this%buildings(b)%units(bu)%schedule_appliances)
                this%buildings(b)%units(bu)%thermal_zone%schedules(1:8760,S_APPLIANCES)=p_schedule%values(1:8760)*this%buildings(b)%units(bu)%power_density*this%buildings(b)%units(bu)%geometry%floor_area_useful ! W
                ! ventilation
                p_schedule=>this%get_schedule(this%buildings(b)%units(bu)%schedule_ventilation)
                this%buildings(b)%units(bu)%thermal_zone%schedules(1:8760,S_VENTILATION)=p_schedule%values(1:8760)*this%buildings(b)%units(bu)%ventilation_rate*this%buildings(b)%units(bu)%geometry%volume_useful/3600d0 ! m3/s
                ! heating set points    
                p_schedule=>this%get_schedule(this%buildings(b)%units(bu)%schedule_heating_setpoint)
                this%buildings(b)%units(bu)%thermal_zone%schedules(1:8760,S_TH_SET)=p_schedule%values(1:8760)! °C
                ! cooling set points                
                p_schedule=>this%get_schedule(this%buildings(b)%units(bu)%schedule_cooling_setpoint)
                this%buildings(b)%units(bu)%thermal_zone%schedules(1:8760,S_TC_SET)=p_schedule%values(1:8760)! °C 
                ! dhw use
                p_schedule=>this%get_schedule(this%buildings(b)%units(bu)%schedule_dhw)
                this%buildings(b)%units(bu)%thermal_zone%schedules(1:8760,S_DHW)=p_schedule%values(1:8760)*people*this%buildings(b)%units(bu)%dhw_use/3600 ! kg/s
                ! internal gain
                if (this%buildings(b)%units(bu)%schedule_internal_gain=='DEFAULT') then
                    this%buildings(b)%units(bu)%thermal_zone%schedules(1:8760,S_INTERNAL_GAIN)=0d0
                else
                    p_schedule=>this%get_schedule(this%buildings(b)%units(bu)%schedule_internal_gain)
                    this%buildings(b)%units(bu)%thermal_zone%schedules(1:8760,S_INTERNAL_GAIN)=p_schedule%values(1:8760) ! W
                end if    
                ! power outage 
                if (this%buildings(b)%units(bu)%schedule_power_available=='DEFAULT') then
                    this%buildings(b)%units(bu)%thermal_zone%schedules(1:8760,S_POWER_AVAILABLE)=1d0
                else
                    p_schedule=>this%get_schedule(this%buildings(b)%units(bu)%schedule_power_available)
                    this%buildings(b)%units(bu)%thermal_zone%schedules(1:8760,S_POWER_AVAILABLE)=p_schedule%values(1:8760) ! W
                end if    
            end do
        end do  
        ! calculate energy needs
        do b=1,size(this%buildings)
            call this%buildings(b)%create_dhw_systems()
            call simlog%output('Building '//Int2Str(b))
            call this%buildings(b)%calculate_energy_needs()
        end do    
    end subroutine
    
    ! psychrometric functions
    function get_X(T,Patm,RH) result(X)
        double precision::X,T,Patm,RH,Psat ! RH 0..1, T °C
        Psat=611.2*exp((17.62*T)/(243.12+T))
        X=0.622*RH*Psat/(Patm-RH*Psat)
    end function
    
    function get_H(T,X) result(H)
        double precision::H,T,X    
        H=CP_AIR*T+X*(CP_VAP*T+HFG_WATER)
    end function
    
    function get_T(H,X) result(T)
        double precision::T,H,X  
        T=(H-HFG_WATER*X)/(CP_AIR+X*CP_VAP)
    end function

    function get_Tdew(X,Patm) result(Tdew)
        double precision::Tdew,Tmax,Tavg,X,Patm 
        integer::iflag
        Tdew=-20d0
        Tmax=100d0
        Tavg=0.5*(Tdew+Tmax)
        call fzero(errTdew,Tdew,Tmax,Tavg,1d-3,1d-6,iflag)
        if (iflag>2) then
            call simlog%output('Error: function get_Tdew failed')
            stop
        end if    
    contains
        subroutine errTdew(err,T)
            double precision::err,T
            err=get_X(T,Patm,1d0)-X
        end subroutine
    end function
    
    subroutine building_calculate_energy_needs(this)
        class(building)::this
        double precision::Q_h_max,Q_c_max,Q_h_need,Q_c_need,T_op_0,T_op_max,Th_set,Tc_set,Patm,RH,m_ven,G_hum,G_deh,G_int,G_room,X_min,X_max,X_air,X_amb,T_amb,T_mean,T_main,T_max,T_min,omega,z,to,Ao,gam,alpha
        double precision::power_available
        double precision,dimension(365)::T_day
        integer::bu,el,num_nodes,hour,step,info,i,j,d,calc_mode
        logical::recalc
		calc_mode=1
        ! calculate size of vector T (temperatures of air node and all other nodes)
        do bu=1,size(this%units)
            call simlog%output('- Unit '//Int2Str(bu))
            this%units(bu)%thermal_zone%X_air=0d0
            Q_h_max=this%units(bu)%geometry%floor_area_useful*this%units(bu)%Q_h_specific ! max heating power, W
            Q_c_max=this%units(bu)%geometry%floor_area_useful*this%units(bu)%Q_c_specific ! max cooling power, W
            num_nodes=1 ! air node
            do el=1,size(this%units(bu)%thermal_zone%external_structures)
                num_nodes=num_nodes+5 ! five nodes in each opaque elements, 1=external, 5=internal  
            end do    
            do el=1,size(this%units(bu)%thermal_zone%windows)
                num_nodes=num_nodes+2 ! two nodes in each opaque elements, 1=external, 2=internal  
            end do            
            ! allocate vectors (if not yet allocated in a multi-year simulation 
            if (.not.allocated(this%units(bu)%T_old)) allocate(this%units(bu)%T_old(num_nodes))
            if (.not.allocated(this%units(bu)%T_new)) allocate(this%units(bu)%T_new(num_nodes))
            if (.not.allocated(this%units(bu)%B)) allocate(this%units(bu)%B(num_nodes))
            if (.not.allocated(this%units(bu)%A)) allocate(this%units(bu)%A(num_nodes,num_nodes))
            if (.not.allocated(this%units(bu)%invAww)) allocate(this%units(bu)%invAww(num_nodes-1,num_nodes-1))
            ! integration along 1 year + warm-up period
            this%units(bu)%T_old=15d0
            this%units(bu)%T_new=0d0
            this%units(bu)%B=0d0
            this%units(bu)%A=0d0
            this%units(bu)%invAww=0d0
			d=0
			T_day=0d0
            do step=-800,8760
                hour=step
                if (step<=0) hour=hour+8760
                T_amb=this%units(bu)%p_building_complex%meteodata%get_value(hour*3600d0-1800d0,meteo_Tdb)
				if (step>0) then
					if (mod(hour,24)==1) d=d+1
					T_day(d)=T_day(d)+T_amb/24
                    power_available=this%units(bu)%thermal_zone%schedules(hour,S_POWER_AVAILABLE)
                else
                    power_available=1d0
				end if   
                Patm=this%units(bu)%p_building_complex%meteodata%get_value(hour*3600d0-1800d0,meteo_Patm)
                RH=this%units(bu)%p_building_complex%meteodata%get_value(hour*3600d0-1800d0,meteo_RH)
                ! free floating    
                call this%units(bu)%calculate_AB(power_available,hour,0d0)
				if (calc_mode==0) then
					call linsolve(this%units(bu)%A,this%units(bu)%B,this%units(bu)%T_new,num_nodes,info)
				else if (calc_mode==1) then
					call linsolve_partitions(this%units(bu)%A,this%units(bu)%B,this%units(bu)%T_new,this%units(bu)%invAww,this%units(bu)%T_old(1),num_nodes,info)
				end if	
                if (info.ne.0) then
                    call simlog%output('Error: solution not found (building energy needs).')
                    stop
                end if    
                T_op_0=this%units(bu)%calculate_T_op()  
                this%units(bu)%thermal_zone%T_op(hour)=T_op_0
                this%units(bu)%thermal_zone%Q_sh_sen(hour)=0d0
                this%units(bu)%thermal_zone%Q_sc_sen(hour)=0d0
                ! control set-points and activate heating or cooling
                Th_set=this%units(bu)%thermal_zone%schedules(hour,S_TH_SET)
                Tc_set=this%units(bu)%thermal_zone%schedules(hour,S_TC_SET)
                if (power_available>0d0) then
                    if (T_op_0<Th_set .or. &
                        T_op_0>Tc_set) then
                        if (T_op_0<Th_set) then
                            call this%units(bu)%calculate_AB(power_available,hour,Q_h_max)
                            this%units(bu)%thermal_zone%Q_sh_sen(hour)=Q_h_max
                            this%units(bu)%thermal_zone%Q_sc_sen(hour)=0d0
                        end if    
                        if (T_op_0>Tc_set) then
                            call this%units(bu)%calculate_AB(power_available,hour,-Q_c_max)
                            this%units(bu)%thermal_zone%Q_sh_sen(hour)=0d0
                            this%units(bu)%thermal_zone%Q_sc_sen(hour)=-Q_c_max
                        end if   
					    if (calc_mode==0) then
						    call linsolve(this%units(bu)%A,this%units(bu)%B,this%units(bu)%T_new,num_nodes,info)
					    else if (calc_mode==1) then
						    call linsolve_partitions(this%units(bu)%A,this%units(bu)%B,this%units(bu)%T_new,this%units(bu)%invAww,this%units(bu)%T_old(1),num_nodes,info)
					    end if	
                        if (info.ne.0) then
                            call simlog%output('Error: solution not found (building energy needs).')
                            stop
                        end if      
                        T_op_max=this%units(bu)%calculate_T_op() 
                        this%units(bu)%thermal_zone%T_op(hour)=T_op_max
                        recalc=.false.
                        if (T_op_0<Th_set) then
                            Q_h_need=Q_h_max*(Th_set-T_op_0)/(T_op_max-T_op_0)
                            if (Q_h_need < Q_h_max) then
                                call this%units(bu)%calculate_AB(power_available,hour,Q_h_need)    
                                this%units(bu)%thermal_zone%Q_sh_sen(hour)=Q_h_need
                                this%units(bu)%thermal_zone%Q_sc_sen(hour)=0d0                            
                                recalc=.true.
                            end if    
                        end if    
                        if (T_op_0>Tc_set) then
                            Q_c_need=Q_c_max*(Tc_set-T_op_0)/(T_op_max-T_op_0)
                            if (Q_c_need < Q_c_max) then
                                call this%units(bu)%calculate_AB(power_available,hour,-Q_c_need) 
                                recalc=.true.
                                this%units(bu)%thermal_zone%Q_sh_sen(hour)=0d0
                                this%units(bu)%thermal_zone%Q_sc_sen(hour)=-Q_c_need
                            end if    
                        end if  
                        if (recalc) then
						    if (calc_mode==0) then
							    call linsolve(this%units(bu)%A,this%units(bu)%B,this%units(bu)%T_new,num_nodes,info)
						    else if (calc_mode==1) then
							    call linsolve_partitions(this%units(bu)%A,this%units(bu)%B,this%units(bu)%T_new,this%units(bu)%invAww,this%units(bu)%T_old(1),num_nodes,info)
						    end if		
                            if (info.ne.0) then
                                call simlog%output('Error: solution not found (building energy needs).')
                                stop
                            end if        
                            this%units(bu)%thermal_zone%T_op(hour)=this%units(bu)%calculate_T_op()
                        end if    
                    end if
                end if    
                this%units(bu)%thermal_zone%T_air(hour)=this%units(bu)%T_new(1)    
                this%units(bu)%T_old=this%units(bu)%T_new
                ! calculate latent load
                X_min=get_X(this%units(bu)%thermal_zone%T_air(hour),Patm,this%units(bu)%min_RH)
                X_max=get_X(this%units(bu)%thermal_zone%T_air(hour),Patm,this%units(bu)%max_RH)
                X_amb=get_X(T_amb,Patm,RH)
                if (hour==1) then
                    X_air=this%units(bu)%thermal_zone%X_air(8760)
                else
                    X_air=this%units(bu)%thermal_zone%X_air(hour-1)
                end if    
                m_ven=RHO_AIR*this%units(bu)%thermal_zone%schedules(hour,S_VENTILATION)
                G_int=this%units(bu)%thermal_zone%schedules(hour,S_PEOPLE_LAT)/HFG_WATER
                G_hum=m_ven*(X_min-X_amb)-G_int &
                    +this%units(bu)%thermal_zone%air_mass/TIMESTEP*(X_min-X_air)
                G_hum=max(0.d0,G_hum)*power_available
                G_deh=m_ven*(X_max-X_amb)-G_int &
                    +this%units(bu)%thermal_zone%air_mass/TIMESTEP*(X_max-X_air)
                G_deh=min(G_deh,0.d0)*power_available
                ! calculate latent load
                G_room=0d0
                this%units(bu)%thermal_zone%Q_sh_lat(hour)=0d0
                this%units(bu)%thermal_zone%Q_sc_lat(hour)=0d0
                if (G_hum>0) then
                    G_room=G_hum
                    this%units(bu)%thermal_zone%Q_sh_lat(hour)=G_room*HFG_WATER
                end if
                if (G_deh<0) then
                    G_room=G_deh
                    this%units(bu)%thermal_zone%Q_sc_lat(hour)=G_room*HFG_WATER
                end if
                ! update X_air
                this%units(bu)%thermal_zone%X_air(hour)=(m_ven*X_amb+G_room+G_int &
                    +this%units(bu)%thermal_zone%air_mass*X_air/TIMESTEP)/(m_ven+this%units(bu)%thermal_zone%air_mass/TIMESTEP)
                ! electricity consumption (plug loads)
                if (step>0) then
                    this%units(bu)%thermal_zone%W_app(hour)=this%units(bu)%thermal_zone%schedules(hour,S_APPLIANCES)*power_available
                end if    
            end do 
            
			! calculate dhw load 
			T_mean=sum(T_day)/365
			T_max=maxval(T_day)
			T_min=minval(T_day)
			! calculate T_main as ground temperature at 1.5 m depth
			omega=2*pi/(365*24*3600)
			z=1.5 ! depth
			Ao=(T_max-T_min)/2
			to=15d0/2*24*3600
            alpha=5d-7
			gam=(omega/2/alpha)**0.5  
    		do hour=1,8760
				this%units(bu)%dhw%T_main(hour)=T_mean-Ao*exp(-gam*z)*cos(omega*(hour*3600-to)-gam*z)
				this%units(bu)%dhw%mdotcp_dhw(hour)=this%units(bu)%thermal_zone%schedules(hour,S_DHW)*CP_WATER ! W/K
				this%units(bu)%dhw%Q_dhw(hour)=this%units(bu)%dhw%mdotcp_dhw(hour)*(40-this%units(bu)%dhw%T_main(hour)) ! W
            end do	
        end do                    
    end subroutine
    
    subroutine building_unit_calculate_AB(this,power_available,hour,Q_hc)
        class(building_unit)::this
        double precision::power_available ! 1 (electric grid is operational) or 0 (electric outage)
        integer::hour ! hour of the year, 1:8760
        double precision::Q_hc ! power of the heating (>0) and cooling (<0) system, W
        double precision::GbT ! solar beam radiation on the tilt surface. W/m2
        double precision::GdT ! solar diffuse radiation on the tilt surface, W/m2
        double precision::Tamb ! outdoor ambient temperature, °C
        double precision::theta ! solar radiation incidence angle, rad
        double precision::Q_int ! internal sensible heat from people and appliances, W   
        double precision::Q_sol ! solar heat through windows, W
        double precision::H_ven ! heat capacity rate (mdot*cp) of ventilation air, W/K
        double precision::fsol ! solar factor due to shading (0 = shading, 1 = no shading)
        double precision::blinds_b ! shading of window blinds
        double precision::hc_ve,hr_ve,hc_vi,hr_vi,hc_ue,hr_ue,hc_de,hr_de ! heat transfer coefficients, W/(m2,K)
        double precision::q_sky ! heat flux to sky (difference w.r.t to sky temperature) 
        double precision::fc_Qhc,fc_Qint,fc_Qsol ! fraction exchanged by convection
        double precision::h_opaque ! thermal conductance between opaque elements nodes (1 to 5 nodes, thus 4 conductances), left of right of current node
        double precision::h_opaque_r ! thermal conductance between opaque elements nodes (1 to 5 nodes, thus 4 conductances), right of current node
        double precision::k_opaque ! specific capacity for opaque elements nodes (1 to 5 nodes, thus 5 capacities)
        double precision::alpha_sol ! solar absorptance of opaque elements,-
        double precision::flag_adiabatic ! 0 if wall is adiabatic
        integer::month,day_hour,w,r,j,k,n_ext
        ! default values 
        hc_ve=10.4d0
        hr_ve=5d0
        hc_vi=2.7d0
        hr_vi=5d0
        fc_Qhc=0.5d0
        fc_Qint=0.5d0
        fc_Qsol=0.5d0
        alpha_sol=0.6d0
        q_sky=0.5d0*hr_ve*11d0 ! 11 K difference between ambient and sky temperature
        ! get outdoor temperature
        Tamb=this%p_building_complex%meteodata%get_value(hour*3600d0-1800d0,meteo_Tdb)
        ! prepare matrix A and B
        this%A=0d0
        this%B=0d0
        ! calculate month (1:12) and hour of the day (1:24) from hour of the year (1:8760)
        month=1
        do while (24*sum(days_per_month(1:month))-hour<0)
            month=month+1
        end do    
        day_hour=mod(hour,24)
        if (day_hour==0) day_hour=24
        ! calculate energy balance 
        n_ext=size(this%thermal_zone%external_structures)
        Q_sol=0d0
        Q_int=this%thermal_zone%schedules(hour,S_APPLIANCES)*power_available+this%thermal_zone%schedules(hour,S_PEOPLE_SEN)+this%thermal_zone%schedules(hour,S_INTERNAL_GAIN)
        H_ven=this%thermal_zone%schedules(hour,S_VENTILATION)*RHO_AIR*CP_AIR
        do w=1,n_ext  
            ! tilt radiation on external elements' area
            call this%p_building_complex%meteodata%tilt_radiation(hour*3600d0-1800d0,this%thermal_zone%external_structures(w)%tilt,this%thermal_zone%external_structures(w)%azimuth,0.2d0,GbT,GdT,theta)
            ! calculate fsol for the external elements 
            fsol=0d0
            if (allocated(this%thermal_zone%external_structures(w)%fsol)) then
				if (size(this%thermal_zone%external_structures(w)%fsol(:,month,day_hour))>0) then
					fsol=sum(this%thermal_zone%external_structures(w)%fsol(:,month,day_hour))/size(this%thermal_zone%external_structures(w)%fsol(:,month,day_hour))
				end if	
            end if    
            ! calculate Q_sol for the window and add the result to the total Q_sol of the thermal zone
            if (w<=size(this%thermal_zone%windows)) then ! windows and external walls are 1:1
                ! set blinds shading factor
                blinds_b=1d0
                if (GbT+GdT>this%blinds_G) blinds_b=this%blinds_b
                ! GdT+GbT*fsol is the radiation on the external elements, for external walls correspond to that on the related windows
                Q_sol=Q_sol+this%thermal_zone%windows(w)%area*this%thermal_zone%windows(w)%struct%gw*(GdT+GbT*fsol)*blinds_b
                ! Window w, j=1 external node, j=2 internal node
                do j=1,2
                    r=1+n_ext*5+(w-1)*2+j
                    if (j==1) then 
                        ! external node (node 1=T1, node 2=T2):  
                        ! (1/R)*(T1-T2)+(hc+hr)*(T1-Tamb)+q_sky=0 
                        ! [(1/R) + (hc+hr)]*T1 - (1/R)*T2 = (hc+hr)*Tamb - q_sky
                        this%A(r,r)=this%A(r,r)+(hc_ve+hr_ve)+1d0/this%thermal_zone%windows(w)%struct%resistance                
                        this%A(r,r+1)=this%A(r,r+1)-1d0/this%thermal_zone%windows(w)%struct%resistance
                        this%B(r)=this%B(r)+(hc_ve+hr_ve)*Tamb-q_sky
                    else 
                        ! internal node (node 1=T1, node 2=T2):
                        ! (1/R)*(T2-T1)+hc*(T2-T_air_zone)+hr*[T2 - sum(Aj/Atot*Tj)]=[(1-fc_Qhc)*Qhc+(1-fc_Qint)*Qint+(1-fc_Qsol)*Qsol)/Atot 
                        ! - (1/R)*T1+[(1/R) + (hc+hr)]*T2 - hc*T_air_zone - hr*sum(Aj/Atot*Tj) = [(1-fc_Qhc)*Qhc+(1-fc_Qint)*Qint+(1-fc_Qsol)*Qsol)/Atot
                        this%A(r,r)=this%A(r,r)+(hc_vi+hr_vi)+1d0/this%thermal_zone%windows(w)%struct%resistance                
                        this%A(r,r-1)=this%A(r,r-1)-1d0/this%thermal_zone%windows(w)%struct%resistance
                        this%A(r,1)=this%A(r,1)-hc_vi
                        do k=1,n_ext
                            this%A(r,1+(k-1)*5+5)=this%A(r,1+(k-1)*5+5)-hr_vi*this%thermal_zone%external_structures(k)%area/this%thermal_zone%radiant_area
                        end do    
                        do k=1,size(this%thermal_zone%windows)
                            this%A(r,1+n_ext*5+(k-1)*2+2)=this%A(r,1+n_ext*5+(k-1)*2+2)-hr_vi*this%thermal_zone%windows(k)%area/this%thermal_zone%radiant_area         
                        end do    
                        this%B(r)=this%B(r)+((1-fc_Qint)*Q_int+(1-fc_Qsol)*Q_sol+(1-fc_Qhc)*Q_hc)/this%thermal_zone%radiant_area
                    end if 
                end do    
            end if    
            ! External element (wall or slab) w, j=1 external node, j=5 internal node
            do j=1,5
                r=1+(w-1)*5+j
                if (j==1) then
                    ! external node (node 1=T1, node 2=T2):  
                    ! C*(T1-T1_old)/dt+(1/R)*(T1-T2)+(hc+hr)*(T1-Tamb)+q_sky=GT*alpha_sol 
                    ! [C/dt+(1/R) + (hc+hr)]*T1 - (1/R)*T2 = (hc+hr)*Tamb - q_sky + C/dt*T1_old+GT*alpha_sol
                    flag_adiabatic=1d0
                    if (this%thermal_zone%external_structures(w)%is_adiabatic) flag_adiabatic=0d0
                    h_opaque=(this%thermal_zone%external_structures(w)%struct%resistance*this%thermal_zone%external_structures(w)%struct%ISO_R_coeff(1))**-1
                    k_opaque=(this%thermal_zone%external_structures(w)%struct%specific_capacity*this%thermal_zone%external_structures(w)%struct%ISO_C_coeff(1))
                    this%A(r,r)=this%A(r,r)+flag_adiabatic*(hc_ve+hr_ve)+h_opaque+k_opaque/TIMESTEP                
                    this%A(r,r+1)=this%A(r,r+1)-h_opaque 
                    this%B(r)=this%B(r)+flag_adiabatic*(hc_ve+hr_ve)*Tamb-flag_adiabatic*q_sky+k_opaque/TIMESTEP*this%T_old(r)+flag_adiabatic*(GdT+GbT*fsol)*alpha_sol
                 else if (j==5) then
                    ! internal node (node 4=T1, node 5=T2):
                    ! C/dt*(T2-T2_old)+(1/R)*(T2-T1)+hc*(T2-T_air_zone)+hr*[T2 - sum(Aj/Atot*Tj)]=[(1-fc_Qhc)*Qhc+(1-fc_Qint)*Qint+(1-fc_Qsol)*Qsol)/Atot 
                    ! - (1/R)*T1+[C/dt+(1/R) + (hc+hr)]*T2 - hc*T_air_zone - hr*sum(Aj/Atot*Tj) = [(1-fc_Qhc)*Qhc+(1-fc_Qint)*Qint+(1-fc_Qsol)*Qsol)/Atot+C/dt*T2_old
                    h_opaque=(this%thermal_zone%external_structures(w)%struct%resistance*this%thermal_zone%external_structures(w)%struct%ISO_R_coeff(4))**-1
                    k_opaque=(this%thermal_zone%external_structures(w)%struct%specific_capacity*this%thermal_zone%external_structures(w)%struct%ISO_C_coeff(5))
                    this%A(r,r)=this%A(r,r)+(hc_vi+hr_vi)+h_opaque+k_opaque/TIMESTEP                
                    this%A(r,r-1)=this%A(r,r-1)-h_opaque
                    this%A(r,1)=this%A(r,1)-hc_vi
                    do k=1,n_ext
                        this%A(r,1+(k-1)*5+5)=this%A(r,1+(k-1)*5+5)-hr_vi*this%thermal_zone%external_structures(k)%area/this%thermal_zone%radiant_area
                    end do    
                    do k=1,size(this%thermal_zone%windows)
                        this%A(r,1+n_ext*5+(k-1)*2+2)=this%A(r,1+n_ext*5+(k-1)*2+2)-hr_vi*this%thermal_zone%windows(k)%area/this%thermal_zone%radiant_area         
                    end do    
                    this%B(r)=this%B(r)+((1-fc_Qint)*Q_int+(1-fc_Qsol)*Q_sol+(1-fc_Qhc)*Q_hc)/this%thermal_zone%radiant_area+k_opaque/TIMESTEP*this%T_old(r)           
                else 
                    ! intermediate node (e.g., for node 2, node 2=T1, node 3=T2, node 1=T0):
                    ! C/dt*(T1-T1_old)+(1/R)right*(T1-T2)+(1/R)left*(T1-T0)=0 
                    ! - (1/R)left*T0+[C/dt+(1/R)right+(1/R)left]*T2 -(1/R)right*T2 = C/dt*T1_old
                    h_opaque=(this%thermal_zone%external_structures(w)%struct%resistance*this%thermal_zone%external_structures(w)%struct%ISO_R_coeff(j-1))**-1
                    h_opaque_r=(this%thermal_zone%external_structures(w)%struct%resistance*this%thermal_zone%external_structures(w)%struct%ISO_R_coeff(j))**-1
                    k_opaque=(this%thermal_zone%external_structures(w)%struct%specific_capacity*this%thermal_zone%external_structures(w)%struct%ISO_C_coeff(j))
                    this%A(r,r-1)=this%A(r,r-1)-h_opaque
                    this%A(r,r)=this%A(r,r)+h_opaque+h_opaque_r+k_opaque/TIMESTEP
                    this%A(r,r+1)=this%A(r,r+1)-h_opaque_r
                    this%B(r)=this%B(r)+k_opaque/TIMESTEP*this%T_old(r)          
                end if  
            end do    
        end do    
        ! calculate energy balance for the air node (Q_sol is ready at this stage)
        ! MC/dt*(T1-T1_old)+sum(Aj*hc*(T1-Tj)+Hven*(T1-Tamb)=fc_Qhc*Q_hc+fc_Qsol*Q_sol+fc_Qint*Q_int
        ! [MC/dt+sum(Aj*hc)]*T1 - sum(Aj*hc*Tj+Hven) = fc_Qhc*Q_hc+fc_Qsol*Q_sol+fc_Qint*Q_int+MC/dt*T1_old+Hven*Tamb
        r=1 
        this%A(r,r)=this%A(r,r)+this%thermal_zone%internal_capacity/TIMESTEP+H_ven+this%thermal_zone%Htb
        do k=1,n_ext
            this%A(r,r)=this%A(r,r)+hc_vi*this%thermal_zone%external_structures(k)%area
            this%A(r,1+(k-1)*5+5)=this%A(r,1+(k-1)*5+5)-hc_vi*this%thermal_zone%external_structures(k)%area
        end do    
        do k=1,size(this%thermal_zone%windows)
            this%A(r,r)=this%A(r,r)+hc_vi*this%thermal_zone%windows(k)%area
            this%A(r,1+n_ext*5+(k-1)*2+2)=this%A(r,1+n_ext*5+(k-1)*2+2)-hc_vi*this%thermal_zone%windows(k)%area         
        end do 
        this%B(r)=fc_Qhc*Q_hc+fc_Qsol*Q_sol+fc_Qint*Q_int+this%thermal_zone%internal_capacity/TIMESTEP*this%T_old(r)+H_ven*Tamb+this%thermal_zone%Htb*Tamb  
    end subroutine
    
    function building_unit_calculate_T_op(this) result(T_op)
        class(building_unit)::this
        double precision::T_op
        integer::i,idx
        double precision::Trm,Ar,T_air
        idx=1
        T_air=this%T_new(idx)
        Trm=0d0
        Ar=0d0
        do i=1,size(this%thermal_zone%external_structures)
            idx=idx+5
            Trm=Trm+this%T_new(idx)*this%thermal_zone%external_structures(i)%Area
            Ar=Ar+this%thermal_zone%external_structures(i)%Area
        end do    
        do i=1,size(this%thermal_zone%windows)
            idx=idx+2
            Trm=Trm+this%T_new(idx)*this%thermal_zone%windows(i)%Area
            Ar=Ar+this%thermal_zone%windows(i)%Area
        end do    
        Trm=Trm/Ar
        T_op=(Trm+T_air)/2d0
    end function
 
    
    ! create the dhw distribution system for the building
    ! must be called after energy needs calculation
    ! automatically size the dhw system (tank volume, heating capacity)
    ! possible types: 
    ! - 'tank' tank for centralized hot water production
    ! - 'tankless' hot distribution for localized hot water production
    subroutine building_unit_create_dhw(this)
        class(building_unit)::this
        double precision::N_flats,sim_factor,Tdhw,qmax,t1,t2,Tc,Th,VL
        Tdhw=40d0
        N_flats=nint(this%geometry%floor_area_heated/this%average_flat_area)
        sim_factor= N_flats**(-0.5) ! simultaneity
        qmax=this%dhw_peak_flow*N_flats*sim_factor ! l/h, peak flow
        qmax=qmax/3600/RHO_WATER ! m3/s
        t1=this%dhw_peak_duration*3600 ! s
        t2=this%dhw_preheating*3600 ! s
        Tc=this%water_main_temperature
        Th=this%dhw_storage_temperature
        select case (this%dhw_type)
            case (1)
                ! size tank volume (m3) and heater capacity (W)
                VL=qmax*t1
                this%dhw%V_storage=VL*((Tdhw-Tc)/(Th-Tc))*(t2/(t1+t2))
                this%dhw%UA_storage=2.6d0+3.9d0*this%dhw%V_storage**(3d0/4d0)
                this%dhw%Q_heater=RHO_WATER*this%dhw%V_storage*CP_WATER*(Th-Tc)/t2
            case (2)
                this%dhw%V_storage=0d0
                this%dhw%UA_storage=0d0
                ! size heating capacity
                this%dhw%Q_heater=qmax*RHO_WATER*CP_WATER*(Tdhw-Tc)
        end select
    end subroutine
    
    
    subroutine building_unit_create_thermal_zone(this)
        class(building_unit)::this
        integer::i,j,num_segments,num_interfloors
        double precision,dimension(:),allocatable::base
        double precision::envelope_area ! m2
        double precision::radiant_area ! m2
        double precision::envelope_HT ! W/K
        call this%materials_inventory%initialize()
        num_segments=size(this%p_buildings(this%building_id)%footprint_segments)
        allocate(base(num_segments))
        ! check data integrity (correct size of vectors)
        if (num_segments.ne.size(this%window_to_segment)) then
            call simlog%output('Error: number of elements in window_to_segment array is incorrect for building '//this%p_buildings(this%building_id)%name //' unit '//this%name)
            call simlog%output('- number of elements in window_to_segment array: '//int2str(size(this%window_to_segment)))
            call simlog%output('- number of segments in building footprint: '//int2str(num_segments))
            stop
        end if    
        if (num_segments.ne.size(this%wall_is_adiabatic)) then
            call simlog%output('Error: number of elements in wall_is_adiabatic array is incorrect for building '//this%p_buildings(this%building_id)%name //' unit '//this%name)
            call simlog%output('- number of elements in wall_is_adiabatic array: '//int2str(size(this%wall_is_adiabatic)))
            call simlog%output('- number of segments in building footprint: '//int2str(num_segments))
            stop
        end if    
        ! initialize base vector 
        do i=1,num_segments
            base(i)=this%p_buildings(this%building_id)%footprint_segments(i)%length
        end do     
        envelope_area=0d0
        radiant_area=0d0
        envelope_HT=0d0
        ! for each footprint segment, add a window 
        if (allocated(this%thermal_zone%windows)) deallocate(this%thermal_zone%windows)
        allocate(this%thermal_zone%windows(num_segments))
        do i=1,num_segments
            ! window_area = total window area * segment(j)*window_to_segment(j)/sum_j=1:N(segment(j)*window_to_segment(j))
            this%thermal_zone%windows(i)%area=this%geometry%window_area*base(i)*this%window_to_segment(i)/dot_product(base,this%window_to_segment)
            this%thermal_zone%windows(i)%azimuth=this%p_buildings(this%building_id)%footprint_segments(i)%azimuth
            this%thermal_zone%windows(i)%tilt=pi/2
            this%thermal_zone%windows(i)%struct=>this%p_building_complex%get_window(this%name_windows)
            ! TO DO: materials inventory for windows
            envelope_area=envelope_area+this%thermal_zone%windows(i)%area
            radiant_area=radiant_area+this%thermal_zone%windows(i)%area
            envelope_HT=envelope_HT+this%thermal_zone%windows(i)%struct%transmittance*this%thermal_zone%windows(i)%area
        end do  
        ! allocate space for exernal structures walls, floor and roof
        ! NOTE: slabs adjacent to heated spaces are excluded
        i=0
        if (this%floor_end_type==ROOF) i=i+1
        if (this%floor_start_type==UNHEATED .or. this%floor_start_type==GROUND) i=i+1
        if (allocated(this%thermal_zone%external_structures)) deallocate(this%thermal_zone%external_structures)
        allocate(this%thermal_zone%external_structures(num_segments+i))
        ! for each footprint segment, add a wall
        do i=1,num_segments
            this%thermal_zone%external_structures(i)%type=WALL
            ! wall area and area above ground are calculated net of the respective window area 
            this%thermal_zone%external_structures(i)%area=base(i)*this%geometry%height-this%thermal_zone%windows(i)%area
            this%thermal_zone%external_structures(i)%area_above_ground=base(i)*this%geometry%height_above_ground-this%thermal_zone%windows(i)%area
            this%thermal_zone%external_structures(i)%area_below_ground=base(i)*this%geometry%height_below_ground
            this%thermal_zone%external_structures(i)%azimuth=this%p_buildings(this%building_id)%footprint_segments(i)%azimuth
            this%thermal_zone%external_structures(i)%tilt=pi/2
            this%thermal_zone%external_structures(i)%struct=>this%p_building_complex%get_wall_slab(this%name_walls)
            this%thermal_zone%external_structures(i)%is_adiabatic=int_to_logical(this%wall_is_adiabatic(i))
            call this%thermal_zone%external_structures(i)%calculate_materials(this%materials_inventory)            
            this%thermal_zone%external_structures(i)%Htb=this%linear_thermal_bridge*base(i)*this%geometry%number_of_storeys
            this%thermal_zone%Htb=this%thermal_zone%Htb+this%thermal_zone%external_structures(i)%Htb
            ! calculate wall transmittance and resistance to include thermal bdrige effect
            this%thermal_zone%external_structures(i)%transmittance=this%thermal_zone%external_structures(i)%struct%transmittance+ &
                this%thermal_zone%external_structures(i)%Htb/this%thermal_zone%external_structures(i)%area
            this%thermal_zone%external_structures(i)%resistance=1d0/(1d0/this%thermal_zone%external_structures(i)%struct%resistance+ &
                this%thermal_zone%external_structures(i)%Htb/this%thermal_zone%external_structures(i)%area)
            ! update envelope area and HT if wall is not adiabatic
            radiant_area=radiant_area+this%thermal_zone%external_structures(i)%area
            if (.not.this%thermal_zone%external_structures(i)%is_adiabatic) then
                envelope_area=envelope_area+this%thermal_zone%external_structures(i)%area
                envelope_HT=envelope_HT+this%thermal_zone%external_structures(i)%transmittance*this%thermal_zone%external_structures(i)%area
                ! create vertices of external area, the portion above ground
                allocate(this%thermal_zone%external_structures(i)%vertices(4))
                this%thermal_zone%external_structures(i)%vertices(1)=point(this%p_buildings(this%building_id)%footprint_segments(i)%X1,&
                                                                           this%p_buildings(this%building_id)%footprint_segments(i)%Y1,&
                                                                           max(this%geometry%height_floor_start,0d0))
                this%thermal_zone%external_structures(i)%vertices(2)=point(this%p_buildings(this%building_id)%footprint_segments(i)%X2,&
                                                                           this%p_buildings(this%building_id)%footprint_segments(i)%Y2,&
                                                                           max(this%geometry%height_floor_start,0d0))
                this%thermal_zone%external_structures(i)%vertices(3)=point(this%p_buildings(this%building_id)%footprint_segments(i)%X1,&
                                                                           this%p_buildings(this%building_id)%footprint_segments(i)%Y1,&
                                                                           max(this%geometry%height_floor_start+this%geometry%height_above_ground,0d0))
                this%thermal_zone%external_structures(i)%vertices(4)=point(this%p_buildings(this%building_id)%footprint_segments(i)%X2,&
                                                                           this%p_buildings(this%building_id)%footprint_segments(i)%Y2,&
                                                                           max(this%geometry%height_floor_start+this%geometry%height_above_ground,0d0))
                call this%thermal_zone%external_structures(i)%make_vertical_grid(GRID_SIZE)
            end if    
        end do         
        i=i-1
        ! add the floor slab
        if (this%floor_start_type==UNHEATED .or. this%floor_start_type==GROUND) then
            i=i+1
            this%thermal_zone%external_structures(i)%type=SLAB
            this%thermal_zone%external_structures(i)%area=this%geometry%footprint_area
            this%thermal_zone%external_structures(i)%area_above_ground=0d0
            this%thermal_zone%external_structures(i)%area_below_ground=0d0
            this%thermal_zone%external_structures(i)%azimuth=this%geometry%footprint_azimuth
            this%thermal_zone%external_structures(i)%tilt=pi ! NOTE: tilt 0 is valid for the roof, since the floor points downwards the correct tilt is 180 deg.
            this%thermal_zone%external_structures(i)%struct=>this%p_building_complex%get_wall_slab(this%name_slab_basement)
            this%thermal_zone%external_structures(i)%is_adiabatic=.false.
            call this%thermal_zone%external_structures(i)%calculate_materials(this%materials_inventory)
            this%thermal_zone%external_structures(i)%Htb=0d0
            this%thermal_zone%Htb=this%thermal_zone%Htb+this%thermal_zone%external_structures(i)%Htb
            ! calculate wall transmittance and resistance to include thermal bdrige effect
            this%thermal_zone%external_structures(i)%transmittance=this%thermal_zone%external_structures(i)%struct%transmittance+ &
                this%thermal_zone%external_structures(i)%Htb/this%thermal_zone%external_structures(i)%area
            this%thermal_zone%external_structures(i)%resistance=1d0/(1d0/this%thermal_zone%external_structures(i)%struct%resistance+ &
                this%thermal_zone%external_structures(i)%Htb/this%thermal_zone%external_structures(i)%area)            
            envelope_area=envelope_area+this%thermal_zone%external_structures(i)%area
            radiant_area=radiant_area+this%thermal_zone%external_structures(i)%area
            envelope_HT=envelope_HT+this%thermal_zone%external_structures(i)%transmittance*this%thermal_zone%external_structures(i)%area
        end if   
        ! add the roof slab
        if (this%floor_end_type==ROOF) then
            i=i+1
            this%thermal_zone%external_structures(i)%type=SLAB_ROOF
            this%thermal_zone%external_structures(i)%area=this%geometry%footprint_area
            this%thermal_zone%external_structures(i)%area_above_ground=0d0
            this%thermal_zone%external_structures(i)%area_below_ground=0d0
            this%thermal_zone%external_structures(i)%azimuth=this%geometry%footprint_azimuth
            this%thermal_zone%external_structures(i)%tilt=0d0
            this%thermal_zone%external_structures(i)%struct=>this%p_building_complex%get_wall_slab(this%name_roof)
            this%thermal_zone%external_structures(i)%is_adiabatic=.false.
            call this%thermal_zone%external_structures(i)%calculate_materials(this%materials_inventory)
            this%thermal_zone%external_structures(i)%Htb=0d0
            this%thermal_zone%Htb=this%thermal_zone%Htb+this%thermal_zone%external_structures(i)%Htb
            ! calculate wall transmittance and resistance to include thermal bdrige effect
            this%thermal_zone%external_structures(i)%transmittance=this%thermal_zone%external_structures(i)%struct%transmittance+ &
                this%thermal_zone%external_structures(i)%Htb/this%thermal_zone%external_structures(i)%area
            this%thermal_zone%external_structures(i)%resistance=1d0/(1d0/this%thermal_zone%external_structures(i)%struct%resistance+ &
                this%thermal_zone%external_structures(i)%Htb/this%thermal_zone%external_structures(i)%area)            
            envelope_area=envelope_area+this%thermal_zone%external_structures(i)%area
            radiant_area=radiant_area+this%thermal_zone%external_structures(i)%area
            envelope_HT=envelope_HT+this%thermal_zone%external_structures(i)%transmittance*this%thermal_zone%external_structures(i)%area
            ! create vertices of external area, the portion above ground
            allocate(this%thermal_zone%external_structures(i)%vertices(num_segments))
            do j=1,num_segments
                this%thermal_zone%external_structures(i)%vertices(j)=point(this%p_buildings(this%building_id)%footprint_segments(j)%X1,&
                                                                           this%p_buildings(this%building_id)%footprint_segments(j)%Y1,&
                                                                           max(this%geometry%height_floor_start+this%geometry%height_above_ground,0d0))
            end do    
            call this%thermal_zone%external_structures(i)%make_horizontal_grid(GRID_SIZE)
        end if   
        ! allocate space for internal structures 
        ! NOTE: for now, only interfloor slabs are included, while internal partitions are considered negligible, i.e. active capacity included in air node
        num_interfloors=this%floor_end-this%floor_start+2 ! maximum number of slabs including external structures
        if (this%floor_start*this%floor_end<0) num_interfloors=num_interfloors-1 ! correction due to floor numbering convention ...,-2, -1, +1, +2, ...
        if (this%floor_start_type==UNHEATED .or. this%floor_start_type==GROUND) num_interfloors=num_interfloors-1 ! remove lower external slab
        if (this%floor_end_type==ROOF) num_interfloors=num_interfloors-1 ! remove upper external slab
        if (allocated(this%thermal_zone%internal_structures)) deallocate(this%thermal_zone%internal_structures)
        allocate(this%thermal_zone%internal_structures(num_interfloors))
        do i=1,num_interfloors
            this%thermal_zone%internal_structures(i)%type=SLAB
            this%thermal_zone%internal_structures(i)%area=this%geometry%footprint_area
            this%thermal_zone%internal_structures(i)%struct=>this%p_building_complex%get_wall_slab(this%name_slab_interfloor)
            call this%thermal_zone%internal_structures(i)%calculate_materials(this%materials_inventory)
        end do    
        this%thermal_zone%envelope_area=envelope_area
        this%thermal_zone%radiant_area=radiant_area
        this%thermal_zone%envelope_HT=envelope_HT+this%thermal_zone%Htb
        this%thermal_zone%air_mass=this%geometry%volume_useful*1.2 ! 1.2 density of air
        this%thermal_zone%internal_capacity=this%geometry%volume_useful*1.2*1.005*1.5 ! 1.5 factor to include active thermal capacity of forniture and internal partitions
    end subroutine
                   
    subroutine building_complex_calculate_fsol(this,month,hour)
        class(building_complex),target::this
        integer::b,u,s,g,ss,bb,nb,month,hour
        double precision::time_month_hour ! time in seconds for the representative day in the month and the hour 
        type(point),pointer::p
        double precision::r,DX,DY,gamma_s,alpha_s,beta_p,gamma_p,h_p,cos_theta
        double precision,dimension(2)::angles
        ! calculate azimuth and altitude of the Sun for given month and hour
        time_month_hour=3600d0*((sum(days_per_month(1:month-1))+15)*24+hour-0.5d0)
        angles=this%meteodata%get_values(time_month_hour,(/meteo_gammas,meteo_thetaz/))
        gamma_s=angles(1)
        alpha_s=pi/2-angles(2)
        r=1d5
        DX=-r*sin(gamma_s)
        DY=-r*cos(gamma_s)
        do b=1,size(this%buildings)
            do u=1,size(this%buildings(b)%units)
                ! first external structures are walls, 1:1 with footprint_segments
                do s=1,size(this%buildings(b)%footprint_segments)
                    ! calculate angle of incidence
                    beta_p=this%buildings(b)%units(u)%thermal_zone%external_structures(s)%tilt
                    gamma_p=this%buildings(b)%units(u)%thermal_zone%external_structures(s)%azimuth
                    cos_theta=cos(pi/2-alpha_s)*cos(beta_p)+sin(pi/2-alpha_s)*sin(beta_p)*cos(gamma_p-gamma_s)
                    if (alpha_s>1d-2 .and. cos_theta>0) then
                        ! assume surface is in shadow 
                        this%buildings(b)%units(u)%thermal_zone%external_structures(s)%fsol(:,month,hour)=0d0
                        ! for each grid point p
                        do g=1,size(this%buildings(b)%units(u)%thermal_zone%external_structures(s)%grid)
                            p=>this%buildings(b)%units(u)%thermal_zone%external_structures(s)%grid(g)
                            h_p=0d0
                            ! test shadow from all segments != s of same building b
                            do ss=1,size(this%buildings(b)%footprint_segments)
                                if (ss.ne.s) then
                                    h_p=get_projection(this%buildings(b)%footprint_segments(ss),this%buildings(b)%height_ag)
                                end if    
                            end do    
                            ! test shadow from all segments of all buildings !=b
                            do bb=1,size(this%buildings)
                                if (bb.ne.b) then
                                    do ss=1,size(this%buildings(bb)%footprint_segments)
                                        h_p=get_projection(this%buildings(bb)%footprint_segments(ss),this%buildings(bb)%height_ag)
                                    end do    
                                end if    
                            end do    
                            ! test shadow from all segments of all neighboring buildings 
                            do nb=1,size(this%neighbors)
                                do ss=1,size(this%neighbors(nb)%footprint_segments)
                                    h_p=get_projection(this%neighbors(nb)%footprint_segments(ss),this%neighbors(nb)%height_ag)
                                end do    
                            end do   
                            ! if p not in shadow set fsol to 1
                            if (p%z>h_p) this%buildings(b)%units(u)%thermal_zone%external_structures(s)%fsol(g,month,hour)=1d0
                        end do
                    else ! dark hours
                        this%buildings(b)%units(u)%thermal_zone%external_structures(s)%fsol(:,month,hour)=0d0
                    end if
                end do
                ! building's roof
                do s=size(this%buildings(b)%footprint_segments)+1,size(this%buildings(b)%units(u)%thermal_zone%external_structures) 
                    if (this%buildings(b)%units(u)%thermal_zone%external_structures(s)%type==SLAB_ROOF) then
                        ! calculate angle of incidence
                        beta_p=this%buildings(b)%units(u)%thermal_zone%external_structures(s)%tilt
                        gamma_p=this%buildings(b)%units(u)%thermal_zone%external_structures(s)%azimuth
                        cos_theta=cos(pi/2-alpha_s)*cos(beta_p)+sin(pi/2-alpha_s)*sin(beta_p)*cos(gamma_p-gamma_s)
                        if (alpha_s>1d-2 .and. cos_theta>0) then
                            ! assume surface is in shadow 
                            this%buildings(b)%units(u)%thermal_zone%external_structures(s)%fsol(:,month,hour)=0d0
                            ! for each grid point p
                            do g=1,size(this%buildings(b)%units(u)%thermal_zone%external_structures(s)%grid)
                                p=>this%buildings(b)%units(u)%thermal_zone%external_structures(s)%grid(g)
                                h_p=0d0
                                ! test shadow from all segments of all buildings !=b
                                do bb=1,size(this%buildings)
                                    if (bb.ne.b) then
                                        do ss=1,size(this%buildings(bb)%footprint_segments)
                                            h_p=get_projection(this%buildings(bb)%footprint_segments(ss),this%buildings(bb)%height_ag)
                                        end do    
                                    end if    
                                end do    
                                ! test shadow from all segments of all neighboring buildings 
                                do nb=1,size(this%neighbors)
                                    do ss=1,size(this%neighbors(nb)%footprint_segments)
                                        h_p=get_projection(this%neighbors(nb)%footprint_segments(ss),this%neighbors(nb)%height_ag)
                                    end do    
                                end do  
                                ! if p not in shadow set fsol to 1
                                if (p%z>h_p) this%buildings(b)%units(u)%thermal_zone%external_structures(s)%fsol(g,month,hour)=1d0
                            end do    
                        else ! dark hours
                            this%buildings(b)%units(u)%thermal_zone%external_structures(s)%fsol(:,month,hour)=0d0
                        end if                            
                    end if
                end do
            end do
        end do  
        
    contains
    
        function get_projection(segment,height_ag) result(h)
            double precision::X,Y,d,h,height_ag
            type(footprint_edge)::segment
            integer::stat
            call intersect_ss(p%x,p%y,p%x+DX,p%y+DY,segment%X1,segment%Y1,segment%X2,segment%Y2,X,Y,stat)    
            if (stat==1) then 
    	        d=((X-p%x)**2+(Y-p%y)**2)**0.5
                h=max(h_p,height_ag-d*tan(alpha_s))
            else
                h=h_p
            end if
        end function
    
    end subroutine
    
    subroutine building_calculate_geometry(this)
        class(building)::this
        integer::i
        double precision::footprint_area,footprint_azimuth
        ! building geometry data
        footprint_area=0d0
        footprint_azimuth=1d9
        if (allocated(this%footprint_segments)) deallocate(this%footprint_segments)
        allocate(this%footprint_segments(size(this%X)-1))
        do i=1,size(this%X)-1
            footprint_area=footprint_area+this%X(i)*this%Y(i+1)-this%X(i+1)*this%Y(i)
            this%footprint_segments(i)%X1=this%X(i)
            this%footprint_segments(i)%X2=this%X(i+1)
            this%footprint_segments(i)%Y1=this%Y(i)
            this%footprint_segments(i)%Y2=this%Y(i+1)
            call this%footprint_segments(i)%update()
            if (abs(this%footprint_segments(i)%azimuth)<footprint_azimuth) footprint_azimuth=this%footprint_segments(i)%azimuth
        end do    
        footprint_area=abs(footprint_area)/2d0
        ! building unit geometry
        do i=1,size(this%units)
            this%units(i)%geometry%footprint_area=this%units(i)%footprint_to_gross_area*footprint_area
            this%units(i)%geometry%footprint_azimuth=footprint_azimuth
            if (this%units(i)%floor_end>0) then
                this%units(i)%geometry%number_of_storeys_above_ground=this%units(i)%floor_end-max(this%units(i)%floor_start,1)+1
            else
                this%units(i)%geometry%number_of_storeys_above_ground=0
            end if 
            if (this%units(i)%floor_start<0) then
                this%units(i)%geometry%number_of_storeys_below_ground=min(this%units(i)%floor_end,-1)-this%units(i)%floor_start+1
            else
                this%units(i)%geometry%number_of_storeys_below_ground=0
            end if    
            this%units(i)%geometry%number_of_storeys=this%units(i)%geometry%number_of_storeys_above_ground+this%units(i)%geometry%number_of_storeys_below_ground
            this%units(i)%geometry%floor_area_gross=this%units(i)%geometry%footprint_area*this%units(i)%geometry%number_of_storeys
            this%units(i)%geometry%floor_area_heated=this%units(i)%geometry%floor_area_gross*this%units(i)%gross_to_heated_area
            this%units(i)%geometry%floor_area_useful=this%units(i)%geometry%floor_area_heated*this%units(i)%heated_to_useful_area 
            this%units(i)%geometry%window_area=this%units(i)%window_to_useful_area*this%units(i)%geometry%floor_area_useful*this%units(i)%geometry%number_of_storeys_above_ground/this%units(i)%geometry%number_of_storeys
            if (this%floors_ag>0) then
                this%units(i)%geometry%height_above_ground=(this%height_ag*this%units(i)%geometry%number_of_storeys_above_ground)/this%floors_ag
            else
                this%units(i)%geometry%height_above_ground=0d0
            end if    
            if (this%floors_bg>0) then
                this%units(i)%geometry%height_below_ground=(this%height_bg*this%units(i)%geometry%number_of_storeys_below_ground)/this%floors_bg
            else
                this%units(i)%geometry%height_below_ground=0d0
            end if
            this%units(i)%geometry%height=this%units(i)%geometry%height_above_ground+this%units(i)%geometry%height_below_ground
            this%units(i)%geometry%height_floor_start=max(0,this%units(i)%floor_start-1)*this%height_ag/max(this%floors_ag,1)+&
                                                      min(0,this%units(i)%floor_start+1)*this%height_bg/max(this%floors_bg,1)
            this%units(i)%geometry%volume_gross=this%units(i)%geometry%floor_area_gross/this%units(i)%geometry%number_of_storeys*this%units(i)%geometry%height
            this%units(i)%geometry%volume_useful=this%units(i)%geometry%floor_area_useful/this%units(i)%geometry%number_of_storeys*this%units(i)%geometry%height*this%units(i)%useful_floor_height
        end do    
   end subroutine

    subroutine building_create_thermal_zones(this)
        class(building)::this
        integer::i
        ! iterate for each building unit 
        do i=1,size(this%units)
            call this%units(i)%create_thermal_zone()
        end do    
    end subroutine
   
     subroutine building_create_dhw_systems(this)
        class(building)::this
        integer::i
        ! iterate for each building unit 
        do i=1,size(this%units)
            call this%units(i)%create_dhw()
        end do    
     end subroutine
     
    subroutine building_unit_calculate_hvac(this)
        class(building_unit)::this
        integer::t
        double precision::aef_hyd,aef_AHU,lat_to_sen,T_Op,Q_fan,Q_air,Q_hr
        double precision::m_ven,m_AHU,m_AHU_0,Patm,RH_OA,T_OA,X_OA,T_R,X_R,H_R,H_OA,Qlat_ven,Qsen_ven,Qlat_noven,Qsen_noven,X_S,DX
        double precision::T_dew,H_dew,T_S,H_S,Q_hyd,T_HR,H_HR,mix_ratio,X_mix,H_mix,T_mix,DH_ahu_lat,DH_ahu_sen,Q_ahu_cool,Q_ahu_heat
        logical::plant_is_on
        this%hvac%Q_heat_hyd=0d0
        this%hvac%Q_cool_hyd=0d0
        this%hvac%Q_heat_AHU=0d0
        this%hvac%Q_cool_AHU=0d0
        this%hvac%W_par_hyd=0d0
        this%hvac%W_par_AHU=0d0
        select case(this%hvac_type)
            case(1)
                aef_hyd=0.025d0 ! auxilary energy for hydraulic loop, fraction of thermal energy
                aef_AHU=0d0 ! auxilary energy for AHU, specific for m3/h
                lat_to_sen=0.3d0 ! latent to sensibile cooling ratio for the fan-coil, fraction
                this%hvac%T_supply_hw=45
                this%hvac%T_supply_cw=7
                do t=1,8760
                    if (this%thermal_zone%Q_sh_sen(t)>0) then
                        ! heating
                            this%hvac%Q_heat_hyd(t)=this%thermal_zone%Q_sh_sen(t)/(this%emission_efficiency*this%control_efficiency*this%distribution_efficiency)
                            this%hvac%W_par_hyd(t)=aef_hyd*this%hvac%Q_heat_hyd(t) 
                    end if    
                    if (this%thermal_zone%Q_sc_sen(t)<0) then
                        ! cooling
                            this%hvac%Q_cool_hyd(t)=this%thermal_zone%Q_sc_sen(t)*(1+lat_to_sen)/(this%emission_efficiency*this%control_efficiency*this%distribution_efficiency)
                            this%hvac%W_par_hyd(t)=this%hvac%W_par_hyd(t)-aef_hyd*this%hvac%Q_cool_hyd(t) 
                    end if
                end do               
            case(2)
                aef_hyd=0.025d0 ! auxilary energy for hydraulic loop, fraction of thermal energy
                aef_AHU=0.5d0 ! auxilary energy for AHU, W specific for m3/h
                lat_to_sen=0.3d0 ! latent to sensibile cooling ratio for the fan-coil, fraction
                this%hvac%T_supply_hw=45
                this%hvac%T_supply_cw=7
                do t=1,8760
                    m_ven=RHO_AIR*this%thermal_zone%schedules(t,S_VENTILATION)
                    T_OA=this%p_building_complex%meteodata%get_value(t*3600d0-1800d0,meteo_Tdb)
                    T_R=this%thermal_zone%T_air(t)
                    T_Op=T_OA-this%heat_recovery_effectiveness*(T_OA-T_R)  ! heat recovery
                    Q_fan=aef_AHU*(m_ven*3600/RHO_AIR/this%distribution_volumetric_efficiency)
                    this%hvac%W_par_AHU(t)=Q_fan
                    if (this%thermal_zone%Q_sh_sen(t)>0d0) then
                        ! heating, CMV (controlled mechanical ventilation)
                        Q_hr=max(m_ven/this%distribution_volumetric_efficiency*CP_AIR*(T_Op-T_OA),0d0)
                        ! heat recovery and 1/2 fan heat 
                        Q_air=min((Q_hr+0.5*Q_fan)*this%distribution_volumetric_efficiency,this%thermal_zone%Q_sh_sen(t)) ! simulate HR control at very low heating loads
                        ! heating, fan-coils
                        this%hvac%Q_heat_hyd(t)=(this%thermal_zone%Q_sh_sen(t)-Q_air)/(this%emission_efficiency*this%control_efficiency*this%distribution_efficiency)
                        this%hvac%W_par_hyd(t)=aef_hyd*this%hvac%Q_heat_hyd(t)        
                    end if    
                    if (this%thermal_zone%Q_sc_sen(t)<0d0) then
                        ! cooling, CMV
                        Q_hr=min(m_ven/this%distribution_volumetric_efficiency*CP_AIR*(T_Op-T_OA),0d0)
                        ! heat recovery and 1/2 fan heat
                        Q_air=max((Q_hr+0.5*Q_fan)*this%distribution_volumetric_efficiency,this%thermal_zone%Q_sc_sen(t)) ! simulate HR control at very low cooling loads
                        ! cooling, fan-coils
                        this%hvac%Q_cool_hyd(t)=(this%thermal_zone%Q_sc_sen(t)-Q_air)*(1+lat_to_sen)/(this%emission_efficiency*this%control_efficiency*this%distribution_efficiency)
                        this%hvac%W_par_hyd=-aef_hyd*this%hvac%Q_cool_hyd(t)       
                    end if
                end do               
            case(3,4)
                ! sensible and latent heat due to ventilation only
                aef_hyd=0.025d0 ! auxilary energy for hydraulic loop, fraction of thermal energy
                aef_AHU=0.7d0 ! auxilary energy for AHU, W specific for m3/h
                select case(this%hvac_type)
                    case(3)
                        this%hvac%T_supply_hw=45
                        this%hvac%T_supply_cw=7
                    case(4)
                        this%hvac%T_supply_hw=45 ! same generation system for radiant floor and AHU
                        this%hvac%T_supply_cw=7 ! same generation system for radiant floor and AHU
                end select    
                ! primary air AHU works at constant air flow rate, modulating ventilation through mixing damper with return air    
                ! calculate mass flow rate of AHU as the maximum between the hourly maximum ventilation and the fixed flow rate set by air volume changes
                ! in principle, user inputs should be such that the second term is higher than the first term
                ! this approach is safer since it automatically corrects incoherent inputs on ventilation (fresh air) and AHU flow rate
                m_AHU_0=RHO_AIR*max(maxval(this%thermal_zone%schedules(1:8760,S_VENTILATION)),this%maximum_airvolume_changes*this%geometry%volume_useful/3600d0)
                Q_fan=aef_AHU*(m_AHU*3600/RHO_AIR/this%distribution_volumetric_efficiency)   
                ! TODO: Q_fan not included in the energy balance of the AHU (lower heating, higher cooling!)
                do t=1,8760
                    if (this%thermal_zone%Q_sh_sen(t)==0d0 .and. this%thermal_zone%Q_sc_sen(t)==0d0 .and. this%thermal_zone%Q_sh_lat(t)==0d0 .and. this%thermal_zone%Q_sc_lat(t)==0d0) then
                        plant_is_on=.false.
                    else
                        plant_is_on=.true.
                    end if
                    if (plant_is_on) then
                        !if (t==17) then
                        !    m_ven=0d0
                        !end if    
                        m_AHU=m_AHU_0
                        m_ven=RHO_AIR*this%thermal_zone%schedules(t,S_VENTILATION)
                        Patm=this%p_building_complex%meteodata%get_value(t*3600d0-1800d0,meteo_Patm)
                        RH_OA=this%p_building_complex%meteodata%get_value(t*3600d0-1800d0,meteo_RH)
                        T_OA=this%p_building_complex%meteodata%get_value(t*3600d0-1800d0,meteo_Tdb)
                        X_OA=get_X(T_OA,Patm,RH_OA)
                        T_R=this%thermal_zone%T_air(t)
                        X_R=this%thermal_zone%X_air(t)
                        H_R=get_H(T_R,X_R)
                        H_OA=get_H(T_OA,X_OA)
                        Qlat_ven=m_ven*HFG_WATER*(X_R-X_OA)
                        Qsen_ven=m_ven*(H_R-H_OA)-Qlat_ven
                        Qlat_noven=this%thermal_zone%Q_sh_lat(t)+this%thermal_zone%Q_sc_lat(t)-Qlat_ven
                        Qsen_noven=this%thermal_zone%Q_sh_sen(t)+this%thermal_zone%Q_sc_sen(t)-Qsen_ven
                        ! AHU outdoor air after heat recovery (HR)
                        T_HR=T_OA-this%heat_recovery_effectiveness*(T_OA-T_R)
                        H_HR=get_H(T_HR,X_OA)
                        ! AHU mixing of return air with outdoor air after heat recovery (HR)
                        mix_ratio=m_ven/m_AHU
                        X_mix=X_OA*mix_ratio+X_R*(1-mix_ratio)
                        H_mix=H_HR*mix_ratio+H_R*(1-mix_ratio)
                        T_mix=get_T(H_mix,X_mix)                        
                        ! if dehumidification (overall, including ventilation air) is needed, check that X_S>=0.0085 (Tdew cannot be lower than 12 °C)
                        if (this%thermal_zone%Q_sh_lat(t)+this%thermal_zone%Q_sc_lat(t)<0d0) then
                            ! calculate X_S, humidity ratio (kg/kg) of supply air necessary to abate internal latent load
                            ! m_AHU*(X_R-X_S)=-Qlat_noven/HFG => X_S=X_R+Qlat_noven/(m_AHU*HFG)
                            X_S=X_R+Qlat_noven/(m_AHU*HFG_WATER) 
                            if (X_S<0.0085) then
                                X_S=0.0085
                                ! unusual (nearly impossible) conditions, treated only to avoid numerical problems
                                if (X_R > X_S) then
                                    if (m_AHU < -Qlat_noven/HFG_WATER/(X_R-X_S)) then
                                        ! since internal latent load cannot be covered, an automatic correction to the internal latent load is applied
                                        Qlat_noven=-m_AHU*HFG_WATER/(X_R-X_S)
                                    end if  
                                else
                                    ! if the room is already dry, internal latent load is not relevant for the AHU
                                    Qlat_noven=0d0
                                end if  
                            end if 
                        ! humidification, bring supply air right below room conditions    
                        else if (this%thermal_zone%Q_sh_lat(t)+this%thermal_zone%Q_sc_lat(t)>0d0) then
                            X_S=X_R+Qlat_noven/(m_AHU*HFG_WATER) 
                        else    
                            X_S=X_mix
                        end if     
                        ! calculate T_S so that the difference H_R-H_S at X_S can cover the sensible load net of ventilation Qsen_noven (partialization)
                        ! primary air T_S = 30 °C for heating, Tdew + 4 K for cooling, m_AHU fixed
                        if (Qsen_noven<0) then
                            ! calculate air conditions at the cooling coil outlet, considering 2 K difference with dew point due to bypass of air 
                            T_dew=get_Tdew(X_S,Patm)
                            H_dew=get_H(T_dew+2d0,X_S)
                            T_S=max(get_T(H_R+(Qsen_noven+Qlat_noven)/m_AHU,X_S),T_dew+4d0)
                        else if (Qsen_noven>0) then
                            T_S=min(get_T(H_R+(Qsen_noven+Qlat_noven)/m_AHU,X_S),30d0)
                        else
                            T_S=T_R
                        end if   
                        H_S=get_H(T_S,X_S)
                        ! calculate contribution of hydronic system (fan-coil or radiant floor)
                        Q_hyd=(Qsen_noven+Qlat_noven)-m_AHU*(H_S-H_R)
                        ! the real AHU and ventilation flow rates (in the AHU) must account for air leakages in air ducts, => /distribution_volumetric_efficiency
                        m_AHU=m_AHU/this%distribution_volumetric_efficiency
                        m_ven=m_ven/this%distribution_volumetric_efficiency                    
                        ! calculate enthalpy change of supply air between S and mix (the task of the AHU heating/cooling coils and humidifier) 
                        DH_ahu_lat=(X_S-X_mix)*HFG_WATER ! latent heating/cooling
                        DH_ahu_sen=H_S-H_mix-DH_ahu_lat ! sensible heating/cooling
                        ! mode 1: heating and dehumidification
                        if (DH_ahu_sen>=0 .and. DH_ahu_lat<0) then
                            Q_ahu_cool=(H_dew-H_mix)*m_AHU
                            Q_ahu_heat=(H_S-H_dew)*m_AHU
                        ! mode 2: heating and humidification
                        else if (DH_ahu_sen>=0 .and. DH_ahu_lat>=0) then   
                            Q_ahu_cool=0d0
                            Q_ahu_heat=(DH_ahu_sen+DH_ahu_lat)*m_AHU
                        ! mode 3: cooling and dehumidification
                        else if (DH_ahu_sen<0 .and. DH_ahu_lat<0) then    
                            Q_ahu_cool=(H_dew-H_mix)*m_AHU
                            Q_ahu_heat=(H_S-H_dew)*m_AHU ! reheat
                        ! mode 4: cooling and humidification
                        else if (DH_ahu_sen<0 .and. DH_ahu_lat>=0) then   
                            if (abs(DH_ahu_sen)>DH_ahu_lat) then
                                ! evaporative cooling + sensible cooling
                                Q_ahu_cool=(DH_ahu_sen+DH_ahu_lat)*m_AHU
                                Q_ahu_heat=0d0
                            else 
                                ! evaporative cooling + sensible heating   
                                Q_ahu_cool=0d0
                                Q_ahu_heat=(DH_ahu_sen+DH_ahu_lat)*m_AHU
                            end if   
                        end if   
                    else
                       Q_ahu_heat=0d0
                       Q_ahu_cool=0d0
                       Q_fan=0d0
                       Q_hyd=0d0
                    end if                
                    if (Q_hyd<0) then
                        this%hvac%Q_cool_hyd(t)=Q_hyd/(this%emission_efficiency*this%control_efficiency*this%distribution_efficiency)
                        this%hvac%W_par_hyd=-aef_hyd*this%hvac%Q_cool_hyd(t)                        
                    else
                        this%hvac%Q_heat_hyd(t)=Q_hyd/(this%emission_efficiency*this%control_efficiency*this%distribution_efficiency)
                        this%hvac%W_par_hyd=aef_hyd*this%hvac%Q_heat_hyd(t)                        
                    end if 
                    this%hvac%Q_heat_AHU(t)=Q_ahu_heat/this%distribution_efficiency
                    this%hvac%Q_cool_AHU(t)=Q_ahu_cool/this%distribution_efficiency
                    ! parasitic energy update 
                    ! NOTE pump work to feed AHU heating/cooling coils adds to hydronic parasitic consumption
                    this%hvac%W_par_AHU(t)=Q_fan
                    this%hvac%W_par_hyd(t)=this%hvac%W_par_AHU(t)+aef_hyd*(this%hvac%Q_heat_AHU(t)-this%hvac%Q_cool_AHU(t))
                end do    
            case default
                call simlog%output('Error: hvac type not found ('//int2str(this%hvac_type)//')')
                stop
        end select    
    end subroutine
    
end module