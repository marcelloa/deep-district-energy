module datalog
    use json_module
    use iso_fortran_env, only: output_unit
    implicit none
    type logger
        character(500)::logfile
        logical::is_open=.false.
        integer::file_id=0
    contains
        procedure::output=>logger_output
    end type

    type(logger)::simlog,simstat
	
	logical::debugging=.false.
    
    type filenames
		character(:),allocatable::input_file_path_str
		character(:),allocatable::input_file_name_ext
		character(:),allocatable::data_path_name
        character(:),allocatable::output_file_name_ext
        character(:),allocatable::output_file_name_bin
        character(:),allocatable::scheduler_file_name_ext
	end type

	type(filenames)::fn

contains

    subroutine save_series_bin(json,json_var,var_name,var_data)
		type(json_core),intent(inout)::json
		type(json_value),pointer::json_var,struct
		character(*)::var_name
		double precision,dimension(:),intent(in)::var_data
		double precision,dimension(2)::meta_data
		integer(kind=8)::offset_bytes
		integer::u,num_data
		inquire(file=fn%output_file_name_bin, size=offset_bytes)
		open(newunit=u, file=fn%output_file_name_bin, access="stream", &
			form="unformatted", status="unknown", action="write")
		write(u,pos=offset_bytes+1) var_data
		close(u)
		num_data=size(var_data)
		meta_data=[dble(offset_bytes),dble(num_data)]
		call json%create_object(struct,var_name) 
		call json%add(struct,'format','binary')
		call json%add(struct,'info',meta_data)
		call json%add(json_var,struct)
    end subroutine
	
	subroutine read_series_bin(json,p_var,var_name,var_data)
       type(json_core),intent(inout)::json
	   type(json_value),pointer::p_var
	   character(*),intent(in)::var_name
	   double precision,dimension(:),allocatable::meta_data
	   double precision,dimension(8760),intent(out)::var_data
	   integer::u
	   call json%get(p_var,var_name//'.info',meta_data)
	   if (size(meta_data)==2) then
		   open(newunit=u, file=fn%output_file_name_bin, access="stream", &
				form="unformatted", status="old", action="read")
		   read(u, pos=dint(meta_data(1))+1) var_data(1:dint(meta_data(2)))
		   close(u)
		else
			call simlog%output('Error: metadata of binary series are incorrect, check variable '//var_name//' in the json file '//fn%output_file_name_bin)
			stop
		end if	
	end subroutine 

    subroutine logger_output(this,message,input_mode)
        class(logger)::this
        integer :: io
        character,optional::input_mode
        character::mode
        character(*)::message
        logical :: exists
		if (present(input_mode)) then
			mode=input_mode
		else
			mode='A'
		end if	
        inquire(file=this%logfile, exist=exists)
        if (.not.this%is_open) then
            if (exists) then
                if (mode=='A'.OR.mode=='a') then
                    open(newunit=io, file=this%logfile, position="append", status="old", action="write")
                else    
                    open(newunit=io, file=this%logfile, status="replace", action="write")
                end if        
            else
                open(newunit=io, file=this%logfile, status="new", action="write")    
            end if        
            this%is_open=.true.
            this%file_id=io
        else
            io=this%file_id
        end if    
        if (mode=='a') then
            write(io,'(A)', advance="no") message
            flush(io)
            write(*,'(A)', advance="no") message
            flush(output_unit)
        else
            write(io,'(A)') message
            write(*,'(A)'), message 
            close(io)
            this%is_open=.false.
            this%file_id=0
        end if    
    end subroutine
    
    function Int2Str(i)
        character(100)::tmp
        character(:),allocatable::Int2Str
        integer::i
        write(tmp,"(I0)") i
        Int2Str=trim(tmp)
    end function
    
    function getPath(filename)
        integer::length,i
        character(*)::filename
        character(:),allocatable::getPath
        character(:),allocatable::tmp
        tmp=trim(filename)
        length=len(tmp)
        do i=length,1,-1
            if (tmp(i:i)=='\' .OR. tmp(i:i)=='/') then
                tmp(i:i)=' '
                exit
            else
              tmp(i:i)=' '  
            end if    
        end do
        getPath=trim(tmp)
    end function
    
end module