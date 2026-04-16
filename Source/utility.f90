module utility
	use json_module
    implicit none
    
	type string
		character(:),allocatable::text
	contains
		procedure::substring=>string_substring
    end type
    	
	
    type wrapper
        class(*),pointer::item
    end type
    
    interface
		function equal(a,b) result(test)
			class(*)::a,b
            logical::test
        end function
    end interface    
    
    type flexible
        type(wrapper),dimension(:),allocatable::array
        procedure(equal),pointer,nopass::equal
    contains
		procedure::push=>flexible_push     ! add a new item at the end of the array, increasing its size, and return its index
		procedure::index=>flexible_index   ! find the index of an item object, based on the "equal" function
		procedure::unique=>flexible_unique ! add an item at the end of the array if not found, else return its index
    end type    
    
    type timeseries
    	!" The class for time series t(:),y(:,:)
    	!" that can be used to hold disturbances like
    	!" weather data (e.g., outdoor temperature, solar radiation)
    	!" and scheduled inputs (e.g., hot water load profile).
    	integer::num_columns
        !" is the number of data columns (col_1,col_2,...)
        !" which are used to hold the time series.
        integer::num_rows
        !" is the number of time steps (ts_1,ts_2,...).
        double precision,dimension(:),allocatable::t
        !" is the time vector (1:num_rows).
        double precision,dimension(:,:),allocatable::y
        !" is the data matrix (1:num_rows,1:num_columns).
	contains
        procedure::initialize=>timeseries_initialize
        !" allocates time vector and data matrix.
    end type
	    
	type textbuffer	
		integer::num_tokens
		type(string),dimension(:),allocatable::tokens
		contains
			procedure::split => textbuffer_split
			procedure::splitsep => textbuffer_splitsep
	end type    
	
	type,public:: input_data
		type(string),dimension(:),allocatable::rows
		integer::num_rows
		integer::section_num_rows
		type(string),dimension(:),allocatable::section_rows
	contains
		procedure::read_from_file => input_data_read_from_file
		!procedure::get_string => input_data_get_string
    end type	
    
    type param
        type(string)::name
        integer::val_int
        double precision::val_dbl
        character(:),allocatable::val_str
        double precision,dimension(:),allocatable::arr_dbl
        double precision,dimension(:,:),allocatable::mat_dbl
		type(textbuffer)::buffer
    contains
		procedure::get_keys=>param_get_keys
		procedure::key=>param_key
    end type    
	
contains

    function flexible_index(this,item) result(idx)
        class(flexible),target::this
        class(*),target::item
        integer::idx
        integer::i
        idx=0
        if (.not.associated(this%equal)) then
            write(*,*) 'Error: equal function not associated in flexible array.'
            stop
        end if
        do i=1,size(this%array)
            if (this%equal(this%array(i)%item,item)) then
				idx=i
				exit
            end if    
        end do         
    end function    

    function flexible_unique(this,item) result(idx)
        class(flexible),target::this
        class(*),target::item
        integer::idx
        idx=this%index(item)
        if (idx==0) idx=this%push(item)
    end function    
    
    function flexible_push(this,item) result(idx)
        class(flexible),target::this
        class(*),target::item
        type(wrapper),dimension(:),allocatable::tmp
        integer::idx
        integer::i
        if (.not.allocated(this%array)) then
            idx=1
            allocate(this%array(1))
            this%array(1)%item=>item
        else  
            ! increase array size and add a new element at the end
            allocate(tmp(size(this%array)))
            do i=1,size(this%array)
                tmp(i)=this%array(i)   ! tmp(i)%item=>this%array(i)%item ?
            end do    
            deallocate(this%array)
            allocate(this%array(size(tmp)+1))
            do i=1,size(this%array)-1
                this%array(i)=tmp(i)
            end do    
            idx=size(this%array)
            this%array(idx)%item=>item
        end if    
    end function    
    

    ! allocate a 2-dimensional array only if not allocated or allocated and of different size
    subroutine alloc_array_2(arr,n1,n2)
		double precision,dimension(:,:),allocatable::arr
        integer::n1,n2
		if (allocated(arr)) then
            if (size(arr,1)==n1 .and. size(arr,2)==n2) then
                ! do nothing: size has not changed
            else
				! destroy and reallocate
                deallocate(arr)
                allocate(arr(n1,n2))
                arr=0d0
            end if    
        else
            allocate(arr(n1,n2))
        end if    
    end subroutine
    
    ! true if test_char is within array_char
	function within(test_char,array_char) result(test)
		character(1),intent(in)::test_char
		character,dimension(:),intent(in)::array_char
		logical::test
		integer::i,n
		test=.false.
		n=size(array_char,1)
		do i=1,n
			if (test_char==array_char(i)) then
				test=.true.
			end if
		end do
    end function
    
    subroutine textbuffer_split(this,inputstr)
		class(textbuffer)::this
		type(string)::inputstr
		integer::i,stato,lenstr
		double precision::num
		integer,dimension(50000)::idx_begin,idx_end
		logical::previous_separator
		lenstr=len(inputstr%text)
		this%num_tokens=0
		previous_separator=.false.
		idx_begin(1)=1
		do i=1,lenstr
		    if (within(inputstr%text(i:i),(/' ',char(9),'%',char(10),',',';'/)) .or. i==lenstr) then
				if (.not. previous_separator) then
					this%num_tokens=this%num_tokens+1
					if (i<lenstr) then
						idx_end(this%num_tokens)=i-1
					else
						idx_end(this%num_tokens)=i
					end if
				else
					if (i==lenstr .and. .not. within(inputstr%text(i:i),(/' ',char(9),'%',char(10),',',';'/)) ) then
						this%num_tokens=this%num_tokens+1
						idx_begin(this%num_tokens)=i
						idx_end(this%num_tokens)=i
					end if
				end if
				previous_separator = .true.
				if (inputstr%text(i:i)=='%') exit
			else
				if (previous_separator) then
					idx_begin(this%num_tokens+1)=i
				end if
				previous_separator = .false.
			end if
		end do
		if (allocated(this%tokens)) deallocate(this%tokens)
		allocate(this%tokens(this%num_tokens))
		do i=1,this%num_tokens
			this%tokens(i)%text=inputstr%text(idx_begin(i):idx_end(i))
		end do
	end subroutine


	subroutine textbuffer_splitsep(this,inputstr,sep)
		class(textbuffer)::this
		type(string)::inputstr
		character::sep
		integer::i,stato,lenstr
		double precision::num
		integer,dimension(50000)::idx_begin,idx_end
		lenstr=len(inputstr%text)
		this%num_tokens=0
		idx_begin(1)=1
		do i=1,lenstr
			if (i==lenstr) then
				this%num_tokens=this%num_tokens+1
				idx_end(this%num_tokens)=i
			elseif (inputstr%text(i:i)==sep) then
				this%num_tokens=this%num_tokens+1
				idx_end(this%num_tokens)=i-1
				idx_begin(this%num_tokens+1)=i+1
			end if
		end do
		if (allocated(this%tokens)) deallocate(this%tokens)
		allocate(this%tokens(this%num_tokens))
		do i=1,this%num_tokens
			this%tokens(i)%text=inputstr%text(idx_begin(i):idx_end(i))
		end do
	end subroutine

    subroutine timeseries_initialize(this,num_rows,num_columns)
    	! allocate timeseries data
        class(timeseries)::this
        integer::num_rows,num_columns
        this%num_rows=num_rows
        this%num_columns=num_columns
        if (allocated(this%t)) deallocate(this%t)
        if (allocated(this%y)) deallocate(this%y)
        allocate(this%t(num_rows))
        allocate(this%y(num_rows,num_columns))
    end subroutine
    
	function string_substring(str_source,begin_or_length,length) result(str)
		class(string),intent(in)::str_source
		character(:),allocatable::str
		integer::begin_or_length
		integer,optional::length
		integer::nend,nstart
		integer::i1,i2
		if (.not.present(length)) then
			nend=begin_or_length
			nstart=1
		else
			nend=begin_or_length+length-1
			nstart=begin_or_length
		end if
		i1=min(len(str_source%text),nstart)
		i2=min(len(str_source%text),nend)
		if (i1>0) then 
            str=str_source%text(i1:i2)
        else
            str=''
        end if    
    end function	
    
	function integer_to_text(i) result(text)
		integer::i
		character(len=20)::tmp
		character(:),allocatable::text
		write (tmp, *) i
		text=trim(adjustl(tmp))
    end function
    
    function text_to_integer(str) result(i)
		character(len=*),intent(in) :: str
		integer:: i
		integer:: stat
		read(str,*,iostat=stat)  i
	end function
    
    function to_double(substr) result(dbl)
		character(*)::substr
        integer::ierr
        double precision::dbl
        read(substr,fmt=*,iostat=ierr) dbl
    end function
    
	! read data from input file and save it internally as a list of strings
	subroutine input_data_read_from_file(this,str_file_name)
		class(input_data)::this
		integer::i,istatus
		character(*)::str_file_name
		character(len=1024)::text_row
		! count number of rows
		open (unit=15,file=str_file_name,status='old',action='read',iostat=istatus)
  		i = 0
		if (istatus==0) then
	   		do
	   			i = i +1
	   			read (15,*,iostat=istatus)
		    	if (istatus /=0) exit
	   		end do
		else
			print *, 'ERROR: file ',str_file_name,' not found.'
			stop
		end if
  		close(unit=15)
  		! allocate text rows
  		this%num_rows=i
  		if (allocated(this%rows)) deallocate(this%rows)
  		allocate(this%rows(this%num_rows))
  		! read text rows
		open (unit=15,file=str_file_name,status='old',action='read',iostat=istatus)
  		i = 0
		if (istatus==0) then
	   		do
	   			i = i +1
	   			text_row=''
	   			read (15,'(a)',iostat=istatus) text_row
	   			this%rows(i)%text=trim(text_row)
		    	if (istatus /=0) exit
	   		end do
		else
			print *, 'ERROR: file ', str_file_name, ' not found.'
			stop 
		end if
		! eliminate last row if empty
		if (len(this%rows(i)%text)==0) then
			this%num_rows=this%num_rows-1
		end if
		close(unit=15)
    end subroutine    
    
    subroutine param_get_keys(this)
		class(param)::this
		call this%buffer%splitsep(this%name,'.')
    end subroutine
	
	function param_key(this,i) result(str)
		class(param)::this
		integer::i
		character(:),allocatable::str
		if (i<=this%buffer%num_tokens) then
			str=this%buffer%tokens(i)%text
		else
			str=''
		end if
    end function
    
    subroutine copy_file(input_file, output_file, stat)
        use iso_fortran_env, only : int64
        implicit none

        character(len=*), intent(in)  :: input_file
        character(len=*), intent(in)  :: output_file
        integer,          intent(out) :: stat

        integer :: in, out, ios, i
        integer, parameter :: chunk = 4096
        integer(int64) :: filesize, nfull, rem
        character(len=1) :: buffer(chunk)

        stat = 0
        in  = -1
        out = -1

        ! Get input file size
        inquire(file=input_file, size=filesize, iostat=ios)
        if (ios /= 0) then
            stat = 1    ! cannot inquire input file
            return
        end if

        ! Open input file
        open(newunit=in, file=input_file, access="stream", form="unformatted", &
             status="old", action="read", iostat=ios)
        if (ios /= 0) then
            stat = 2    ! cannot open input file
            return
        end if

        ! Open output file
        open(newunit=out, file=output_file, access="stream", form="unformatted", &
             status="replace", action="write", iostat=ios)
        if (ios /= 0) then
            stat = 3    ! cannot open output file
            close(in)
            return
        end if

        nfull = filesize / chunk
        rem   = mod(filesize, int(chunk, int64))

        ! Copy full chunks
        do i = 1, int(nfull)
            read(in, iostat=ios) buffer
            if (ios /= 0) then
                stat = 4    ! read error
                exit
            end if

            write(out, iostat=ios) buffer
            if (ios /= 0) then
                stat = 5    ! write error
                exit
            end if
        end do

        ! Copy remainder
        if (stat == 0 .and. rem > 0) then
            read(in, iostat=ios) buffer(1:int(rem))
            if (ios /= 0) then
                stat = 4
            else
                write(out, iostat=ios) buffer(1:int(rem))
                if (ios /= 0) stat = 5
            end if
        end if

        if (in  /= -1) close(in)
        if (out /= -1) close(out)

    end subroutine copy_file
    
end module