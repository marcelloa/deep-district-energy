module constants_module
    implicit none

    double precision, parameter::TIMESTEP=3600d0 ! s, 1 hour
    double precision, parameter::GRID_SIZE=2d0 ! m
    double precision, parameter::RHO_AIR=1.2 ! kg/m3
    double precision, parameter::CP_AIR=1005 ! J/(kg,K)
    double precision, parameter::CP_VAP=1806 ! J/(kg,K)
    double precision, parameter::RHO_WATER=1000 ! kg/m3
    double precision, parameter::CP_WATER=4186 ! J/(kg,K)
    double precision, parameter::HFG_WATER=2466000 ! J/kg
    
    integer,dimension(:),allocatable::simulation_years
    integer::current_year
    
contains  

    function previous_step(hour) result (prev_hour)
		integer::hour,prev_hour
		prev_hour=hour-1
		if (prev_hour==0) prev_hour=8760
	end function	

	function next_step(hour) result (next_hour)
		integer::hour,next_hour
		next_hour=hour+1
		if (next_hour>8760) next_hour=1
    end function
    
    function cumulative_step(hour) result (cumulative_hour)
		integer::hour,cumulative_hour
        cumulative_hour=(current_year-1)*8760+hour
    end function
    
end module