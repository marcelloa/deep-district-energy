module meteolib
	use utility
    use datalog
	use json_module
    use maths, only:pi
	implicit none

	integer,parameter::meteo_G=1
	integer,parameter::meteo_Gbn=2
	integer,parameter::meteo_Gd=3
	integer,parameter::meteo_Tdb=4
	integer,parameter::meteo_Tdp=5
	integer,parameter::meteo_RH=6
	integer,parameter::meteo_Patm=7
	integer,parameter::meteo_Wang=8
	integer,parameter::meteo_Wspd=9
	integer,parameter::meteo_Thetaz=10
	integer,parameter::meteo_Gammas=11
    
    double precision,dimension(12),parameter::days_per_month=(/31,28,31,30,31,30,31,31,30,31,30,31/)
    
    type meteo
        type(string)::file_type
	    type(timeseries)::values
    contains
        procedure::update_solar_angles=>meteo_update_solar_angles
        procedure::load_file=>meteo_load_file
        procedure::get_values=>meteo_get_values
        procedure::get_value=>meteo_get_value
        procedure::T_amb=>meteo_T_amb
        procedure::get_T_ground=>meteo_get_T_ground
        procedure::tilt_radiation=>meteo_tilt_radiation
        procedure::to_json=>meteo_to_json
    end type

    contains
    
    subroutine meteo_to_json(this,json,obj,obj_name)
        class(meteo)::this
    	type(json_core),intent(inout)::json
		type(json_value),pointer,intent(inout)::obj
		type(json_value),pointer::obj_meteo
        double precision,dimension(1:8760)::tmp
		character(*)::obj_name
		call json%create_object(obj_meteo,obj_name)
        tmp=this%values%y(1:8760,meteo_Tdb)
		call save_series_bin(json,obj_meteo,'T_amb',tmp)
        call json%add(obj,obj_meteo)
    end subroutine
    
    function meteo_T_amb(this,hour) result(T)
        class(meteo)::this
        integer::hour
        double precision::T
        T=this%get_value(hour*3600d0-1800d0,meteo_Tdb)
    end function

    function meteo_get_T_ground(this) result(T)
        class(meteo)::this
        integer::hour
        double precision::T
        T=0d0
        do hour=1,8760
            T=T+this%get_value(hour*3600d0-1800d0,meteo_Tdb)
        end do
        T=T/8760
    end function

    subroutine meteo_tilt_radiation(this,time,tilt,azimuth,rhoground,RadGbT,RadGdT,theta)
        class(meteo)::this
        double precision::time,tilt,azimuth,rhoground,RadGH,RadDH,RadBN,tetaSun,gammaSun,RadGbT,RadGdT,theta
        double precision,dimension(5)::meteo_data
        meteo_data=this%get_values(time,(/meteo_G,meteo_Gd,meteo_Gbn,meteo_Thetaz,meteo_Gammas/))
        call TiltSurfaceRad_HDKR(time,tilt,azimuth,rhoground,meteo_data(1),meteo_data(2),meteo_data(3),meteo_data(4),meteo_data(5),RadGbT,RadGdT,theta)
     end subroutine    

	! update solar angles based on latitude and longitude (when not provided by the meteo file)
	subroutine meteo_update_solar_angles(this,latitude, longitude)
        class(meteo)::this
		double precision::latitude, longitude ! rad
		double precision::timezone
		integer::r
		timezone=nint(longitude/(15*pi/180)) ! assumed timezone
		do r=1,this%values%num_rows
			call SunAngles(latitude, longitude, timezone,dble(r-1)-0.5d0,this%values%y(r,10),this%values%y(r,11))
		end do
    end subroutine
    
	! create meteo data with the following structure
    ! t		time (s)
    ! y(1) 	G (W/m2)
    ! y(2) 	Gbn (W/m2)
    ! y(3) 	Gd (W/m2)
    ! y(4) 	Tdb (°C)
    ! y(5) 	Tdp (°C)
    ! y(6) 	RH (0..1)
    ! y(7) 	Pressure (Pa)
    ! y(8) 	Wind angle (rad)
    ! y(9) 	Wind speed (m/s)
    ! y(10)	Solar zenith, rad
    ! y(11)	Solar azimuth, rad
	subroutine meteo_load_file(this,file_name)
        class(meteo)::this
        character(*)::file_name
		type(string)::file_name_str,extension,location,latC,longC
		type(textbuffer)::buffer
	    type(input_data)::meteo_input
        character(:),allocatable::substr
		integer::r,j,ierr
		double precision::timezone,latDeg,latMin,latitude,longDeg,longMin,longitude
        file_name_str%text=file_name
		call buffer%splitsep(file_name_str,'.')
		extension=buffer%tokens(buffer%num_tokens)
		this%file_type=extension
		call meteo_input%read_from_file(file_name)
	    select case (extension%text)
	        case('tm2') ! TMY2 meteo file
	        	if (meteo_input%num_rows /= 8761) then
	        		write(*,*) 'Error: tm2 file ',file_name,' is corrupted.'
	        		stop 1061
	        	end if
				call this%values%initialize(meteo_input%num_rows+1,11)
				! read header line
				location%text=meteo_input%rows(1)%substring(8,22)//' '//meteo_input%rows(1)%substring(31,2)
                ! location timezone
                timezone = to_double(meteo_input%rows(1)%substring(34, 3))
                ! latitude data
                latC%text = meteo_input%rows(1)%substring(38, 1)
                latDeg = to_double(meteo_input%rows(1)%substring(40, 2))
                latMin = to_double(meteo_input%rows(1)%substring(43, 2))
                latitude = (latDeg + (1d0 / 60d0) * latMin) / 180d0 * pi
                if (.not. (latC%text=='N')) then
                    latitude = -latitude
                end if
                ! longitude data
                longC%text = meteo_input%rows(1)%substring(46, 1)
                longDeg = to_double(meteo_input%rows(1)%substring(48, 3))
                longMin = to_double(meteo_input%rows(1)%substring(52, 2))
                longitude = (longDeg + (1d0 / 60d0) * longMin) / 180d0 * pi
                if (.not. (longC%text=='W')) then
                    longitude = -longitude
                end if
	            do r=2,this%values%num_rows-1
               		this%values%t(r)=(dble(r-1)-0.5d0)*3600d0
	                this%values%y(r,1) = to_double(meteo_input%rows(r)%substring(18, 4))
	                this%values%y(r,2) = to_double(meteo_input%rows(r)%substring(24, 4))
	                this%values%y(r,3) = to_double(meteo_input%rows(r)%substring(30, 4))
	                this%values%y(r,4) = to_double(meteo_input%rows(r)%substring(68, 4))/10d0
	                this%values%y(r,5) = to_double(meteo_input%rows(r)%substring(74, 4))/10d0
	                this%values%y(r,6) = to_double(meteo_input%rows(r)%substring(80, 3))/100d0
	                this%values%y(r,7) = to_double(meteo_input%rows(r)%substring(85, 4))*100d0 ! Pa
	                this%values%y(r,8) = to_double(meteo_input%rows(r)%substring(91, 3))*pi/180d0 ! rad
	                this%values%y(r,9) = to_double(meteo_input%rows(r)%substring(96, 3))/10d0
					call SunAngles(latitude, longitude, timezone,dble(r-1)-0.5d0,this%values%y(r,10),this%values%y(r,11))
	            end do
	            this%values%t(1)=-this%values%t(2)
	            this%values%y(1,:)=this%values%y(this%values%num_rows-1,:)
	            this%values%t(this%values%num_rows)=this%values%t(this%values%num_rows-1)+1800d0
	            this%values%y(this%values%num_rows,:)=this%values%y(2,:)
	        case('csv') ! EnergyPlus meteo file, csv format, without solar angles
	        	if (meteo_input%num_rows /= 8761) then
	        		write(*,*) 'Error: csv file ',file_name,' is corrupted.'
	        		stop 1062
	        	end if
				call this%values%initialize(meteo_input%num_rows+1,11)
	            do r=2,this%values%num_rows-1
               		this%values%t(r)=(dble(r-1)-0.5d0)*3600d0
               		call buffer%splitsep(meteo_input%rows(r),',')
	                this%values%y(r,1) = to_double(buffer%tokens(12)%text)   		! G, W/m2
	                this%values%y(r,2) = to_double(buffer%tokens(13)%text)   		! Gbn, W/m2
	                this%values%y(r,3) = to_double(buffer%tokens(14)%text)   		! Gd, W/m2
	                this%values%y(r,4) = to_double(buffer%tokens(5)%text)   		! Tdb, °C
	                this%values%y(r,5) = to_double(buffer%tokens(6)%text) 		! Tdp, °C
	                this%values%y(r,6) = to_double(buffer%tokens(7)%text)/100d0 	! rh, 0..1
	                this%values%y(r,7) = to_double(buffer%tokens(8)%text)	 		! Pressure, Pa
	            end do
	            this%values%t(1)=-this%values%t(2)
	            this%values%y(1,:)=this%values%y(this%values%num_rows-1,:)
	            this%values%t(this%values%num_rows)=this%values%t(this%values%num_rows-1)+1800d0
	            this%values%y(this%values%num_rows,:)=this%values%y(2,:)
	        case('met') ! proprietary format meteo file
	        	if (meteo_input%num_rows /= 8761) then
	        		write(*,*) 'Error: met file ',file_name,' is corrupted.'
	        		stop 1063
	        	end if
				call this%values%initialize(meteo_input%num_rows,11)
	            do r=1,this%values%num_rows
               		this%values%t(r)=(dble(r-1)-0.5d0)*3600d0
               		call buffer%split(meteo_input%rows(r))
                	this%values%y(r,1)=to_double(buffer%tokens(2)%text) ! G(W/m2)
                	this%values%y(r,2)=to_double(buffer%tokens(3)%text) ! Gbn (W/m2)
                	this%values%y(r,3)=to_double(buffer%tokens(4)%text) ! Gd (W/m2)
                	this%values%y(r,4)=to_double(buffer%tokens(5)%text) ! Ta (°C)
                	this%values%y(r,6)=to_double(buffer%tokens(6)%text) ! rh (0..1)
                	this%values%y(r,10)=to_double(buffer%tokens(7)%text)! theta_z (rad)
                	this%values%y(r,11)=to_double(buffer%tokens(8)%text)! gamma_s (rad)
	            end do
	    end select
    end subroutine

    subroutine TiltSurfaceRad_HDKR(t,tilt,azimuth,rhoground,RadGH,RadDH,RadBN,tetaSun,gammaSun,RadGbT,RadGdT,theta)
        ! Anisotropic sky model of Hay, Davies, Klucher & Reindl
        ! NOTE, HDKR gives strange peaks when thetaSun is near 90° and a large value of beam radiation exists
        ! To avoid such a strange behavior, the LJ model is used whenever thetaSun > 85°
        double precision::t,tilt,azimuth,rhoground,RadGH,RadDH,RadBN,tetaSun,gammaSun,RadGbT,RadGdT,theta,Gon,nDay,Ib,I,Id,Ai,Rb,cosTheta,f
        if (tetaSun>85d0/180d0*pi) then
            call TiltSurfaceRad_LJ(t,tilt,azimuth,rhoground,RadGH,RadDH,RadBN,tetaSun,gammaSun,RadGbT,RadGdT,theta)
            return
        end if
        nDay = dble(int((t / 3600d0 - 1e-15) / 24d0)) + 1d0
        Gon = 1367d0 * (1d0 + 0.033d0 * cos(2d0 * pi * nDay / 365d0))
        Ib = RadBN * cos(tetaSun)
        I = RadGH
        Id = RadDH
        Ai = RadBN / Gon
        Rb = 1d0
        cosTheta = CosIncidenceAngle(tilt, azimuth,tetaSun, gammaSun)
        if (cos(tetaSun) > 0d0) then
            Rb = cosTheta / cos(tetaSun)
        end if
        f = 0d0
        if (I > 0 .and. Ib > 0) then
            f = (Ib / I)**0.5d0
        end if
        RadGbT = (Ib + Id * Ai) * Rb 
        RadGdT = Id * ((1d0 - Ai) * (1d0 + cos(tilt)) / 2d0 * (1d0 + f * (sin(tilt / 2d0))**3)) + I * rhoground * (1d0 - cos(tilt)) / 2d0
        theta = acos(cosTheta)
    end

    subroutine TiltSurfaceRad_LJ( t, tilt,  azimuth,  rhoground, RadGH, RadDH, RadBN,tetaSun, gammaSun,RadGbT,RadGdT, theta)
    	double precision::t, tilt,  azimuth,  rhoground, RadGH, RadDH, RadBN,tetaSun, gammaSun,RadGbT,RadGdT,theta,cosTheta
        ! Isotropic sky model of Liu & Jordan
        ! See Duffie Beckman, p% 95 / p%25
        ! Remark, RadBN = Ib/cos(tetaz)
        cosTheta=CosIncidenceAngle(tilt, azimuth,tetaSun, gammaSun)
        RadGbT = RadBN * cosTheta  
        RadGdT = RadDH * (1d0 + cos(tilt)) / 2d0 + rhoground * RadGH * (1d0 - cos(tilt)) / 2d0
        theta = acos(cosTheta)
    end

    function CosIncidenceAngle(tilt,  azimuth, tetaSun, gammaSun) result(y)
    	double precision::tilt,  azimuth, tetaSun, gammaSun,y
    	double precision::nZs,nXs,nYs,nZt,nXt,nYt
        ! XYZ frame, X heading south, Y east, Z upward
        nZs = cos(tetaSun)
        nXs = sin(tetaSun) * cos(gammaSun)
        nYs = sin(tetaSun) * sin(-gammaSun)
        nZt = cos(tilt)
        nXt = sin(tilt) * cos(azimuth)
        nYt = sin(tilt) * sin(-azimuth)
        y=max(nXs * nXt + nYs * nYt + nZs * nZt, 0d0)
    end

    subroutine SunAngles(latRad, longRad, timezone, sTimeHrs,sunZenithRad, sunAzimuthRad)
    	double precision::latRad, longRad, timezone, sTimeHrs,sunZenithRad, sunAzimuthRad,latRad0,nDay,B,delta,solarTime,omega,cosTetaSun,E,C1,C2,C3,omegaew,gammap
    	! See Duffie Beckman, pp. 11-16
        nDay = dble(int((sTimeHrs - 1e-15) / 24)) + 1d0
        B = dble(int((sTimeHrs - 1e-15) / 24)) * 2d0 * pi / 365d0 ! Rad
        delta = pi / 180d0 * 23.45d0 * sin(2d0 * pi / 365d0 *(284d0 + nDay))
        E = 229.2d0 * (0.000075d0 + 0.001868d0 * cos(B) - 0.032077d0 * sin(B) &
                       - 0.014615d0 * cos(2d0 * B) - 0.04089d0 * sin(2d0 * B))   ! Min
        solarTime = sTimeHrs * 60d0 + E + 4d0 * (-timezone * 15d0 -longRad / pi * 180d0)  ! Min
        omega = (solarTime / 60d0 - (12d0 + 24d0 * (nDay - 1))) * 15d0 / 180d0 * pi ! Rad
        cosTetaSun = max(cos(latRad) * cos(delta) * cos(omega) + sin(latRad) * sin(delta),0d0)
        sunZenithRad = acos(cosTetaSun)
        latRad0 = latRad
        if (latRad0 == 0d0) then
            latRad0 = 0.001d0 ! nearly zero, to avoid singualarity
        end if
        if (abs(tan(delta) / tan(latRad0)) > 1d0) then
            C1 = 1d0
        else
            omegaew = acos(tan(delta) / tan(latRad0))
            if (abs(omega) < omegaew) then
                C1 = 1d0
            else
                C1 = -1d0
            end if
        end if
        if (latRad0 * (latRad0 - delta) >= 0) then
            C2 = 1d0
        else
            C2 = -1d0
        end if
        if (omega >= 0) then
            C3 = 1d0
        else
            C3 = -1d0
        end if
        gammap = asin(max(min(sin(omega) * cos(delta) / sin(sunZenithRad),1d0),-1d0))
        sunAzimuthRad = C1 * C2 * gammap + C3 * (1d0 - C1 * C2) / 2d0 * pi
    end

    ! interpolate meteo data
   	function meteo_get_values(this,time,codes) result(y)
        class(meteo)::this
		double precision::time,x
		integer,dimension(:)::codes
		double precision,dimension(size(codes))::y
		integer::i,h
	    h=int((time+1800d0)/3600d0)+1
	    x=time-(h-1.5d0)*3600d0
	    do i=1,size(codes)
	    	y(i)=this%values%y(h+1,codes(i))*x/3600d0
	    	y(i)=y(i)+this%values%y(h,codes(i))*(1-x/3600d0)
	    end do
    end function

    ! interpolate meteo data
   	function meteo_get_value(this,time,code) result(y)
        class(meteo)::this
		double precision::time,x
		integer::code
		double precision::y
		integer::i,h
	    h=int((time+1800d0)/3600d0)+1
	    x=time-(h-1.5d0)*3600d0
       	y=this%values%y(h+1,code)*x/3600d0
	   	y=y+this%values%y(h,code)*(1-x/3600d0)
	end function
    
end module
