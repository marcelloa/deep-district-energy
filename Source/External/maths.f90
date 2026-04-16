module maths

    implicit none
    
    double precision::NaV ! undetermined number, used when a function requires a dummy input
    double precision,parameter::Inf=1d32
    double precision,parameter::pi=3.141592653589793d0 !! pi number
    double precision, parameter, private :: HALF=0.5d0,ONE=1.d0,TWO=2.d0,THREE=3.d0

    type,public,abstract :: linear_interp_class
        !! Base class for the linear interpolation types
        private
    contains
        private
        procedure(destroy_func),deferred,public :: destroy  !! destructor
    end type linear_interp_class

    abstract interface
        pure elemental subroutine destroy_func(me)  !! interface for bspline destructor routines
        import :: linear_interp_class
        implicit none
        class(linear_interp_class),intent(inout) :: me
        end subroutine destroy_func
    end interface

    type,extends(linear_interp_class),public :: linear_interp_1d
        !! Class for 1d linear interpolation.
        private
        double precision,dimension(:),allocatable :: f
        double precision,dimension(:),allocatable :: x
        integer :: ilox = 1
        contains
        private
        procedure,public :: initialize => initialize_1d
        procedure,public :: evaluate   => interp_1d
        procedure,public :: destroy    => destroy_1d
        final :: finalize_1d
    end type linear_interp_1d

    type,extends(linear_interp_class),public :: linear_interp_2d
        !! Class for 2d linear interpolation.
        !private
        double precision,dimension(:,:),allocatable :: f
        double precision,dimension(:),allocatable :: x
        double precision,dimension(:),allocatable :: y
        integer :: ilox = 1
        integer :: iloy = 1
        contains
        private
        procedure,public :: initialize => initialize_2d
        procedure,public :: evaluate   => interp_2d
        procedure,public :: destroy    => destroy_2d
        final :: finalize_2d
    end type linear_interp_2d

    type,extends(linear_interp_class),public :: linear_interp_3d
        !! Class for 3d linear interpolation.
        !private
        double precision,dimension(:,:,:),allocatable :: f
        double precision,dimension(:),allocatable :: x
        double precision,dimension(:),allocatable :: y
        double precision,dimension(:),allocatable :: z
        integer :: ilox  = 1
        integer :: iloy  = 1
        integer :: iloz  = 1
        contains
        private
        procedure,public :: initialize => initialize_3d
        procedure,public :: evaluate   => interp_3d
        procedure,public :: destroy    => destroy_3d
        final :: finalize_3d
    end type linear_interp_3d

    type,extends(linear_interp_class),public :: linear_interp_4d
        !! Class for 4d linear interpolation.
        !private
        double precision,dimension(:,:,:,:),allocatable :: f
        double precision,dimension(:),allocatable :: x
        double precision,dimension(:),allocatable :: y
        double precision,dimension(:),allocatable :: z
        double precision,dimension(:),allocatable :: q
        integer :: ilox  = 1
        integer :: iloy  = 1
        integer :: iloz  = 1
        integer :: iloq  = 1
        contains
        private
        procedure,public :: initialize => initialize_4d
        procedure,public :: evaluate   => interp_4d
        procedure,public :: destroy    => destroy_4d
        final :: finalize_4d
    end type linear_interp_4d

    type,extends(linear_interp_class),public :: linear_interp_5d
        !! Class for 5d linear interpolation.
        !private
        double precision,dimension(:,:,:,:,:),allocatable :: f
        double precision,dimension(:),allocatable :: x
        double precision,dimension(:),allocatable :: y
        double precision,dimension(:),allocatable :: z
        double precision,dimension(:),allocatable :: q
        double precision,dimension(:),allocatable :: r
        integer :: ilox  = 1
        integer :: iloy  = 1
        integer :: iloz  = 1
        integer :: iloq  = 1
        integer :: ilor  = 1
        contains
        private
        procedure,public :: initialize => initialize_5d
        procedure,public :: evaluate   => interp_5d
        procedure,public :: destroy    => destroy_5d
        final :: finalize_5d
    end type linear_interp_5d

    type,extends(linear_interp_class),public :: linear_interp_6d
        !! Class for 6d linear interpolation.
        !private
        double precision,dimension(:,:,:,:,:,:),allocatable :: f
        double precision,dimension(:),allocatable :: x
        double precision,dimension(:),allocatable :: y
        double precision,dimension(:),allocatable :: z
        double precision,dimension(:),allocatable :: q
        double precision,dimension(:),allocatable :: r
        double precision,dimension(:),allocatable :: s
        integer :: ilox  = 1
        integer :: iloy  = 1
        integer :: iloz  = 1
        integer :: iloq  = 1
        integer :: ilor  = 1
        integer :: ilos  = 1
        contains
        private
        procedure,public :: initialize => initialize_6d
        procedure,public :: evaluate   => interp_6d
        procedure,public :: destroy    => destroy_6d
        final :: finalize_6d
    end type linear_interp_6d

    contains
    
    
    subroutine find_position(x,x_array,N_step_x,position,x_inf,x_sup)
        !Variable declaration
        !Input
        integer::N_step_x,j
        double precision::x
        double precision,dimension(N_step_x+1)::x_array
        !Internal varibles
        !Output
        integer::position
        double precision::x_inf,x_sup
        position=1
        x_inf=x_array(1)
        x_sup=x_array(1)
        do j=1,N_step_x
            if (x>=x_array(j) .and. x<x_array(j+1)) then
                x_inf=x_array(j)
                x_sup=x_array(j+1)
                position=j
                return
            end if
        end do
    end subroutine find_position
    
    subroutine double_interpolation(y,x,y_inf,y_sup,x_inf,x_sup,val_x_inf_y_inf,val_x_sup_y_inf,val_x_inf_y_sup,val_x_sup_y_sup,val)
        !Input
        double precision::y,x,y_inf,y_sup,x_inf,x_sup,val_x_inf_y_inf,val_x_sup_y_inf,val_x_inf_y_sup,val_x_sup_y_sup
        !Internal variables
        double precision::val_y_inf,val_y_sup
        !Output
        double precision::val
        val_y_inf=val_x_inf_y_inf+(x-x_inf)/(x_sup-x_inf)*(val_x_sup_y_inf-val_x_inf_y_inf)
        val_y_sup=val_x_inf_y_sup+(x-x_inf)/(x_sup-x_inf)*(val_x_sup_y_sup-val_x_inf_y_sup)
        val=val_y_inf+(val_y_sup-val_y_inf)*(y-y_inf)/(y_sup-y_inf)
    end subroutine double_interpolation
    
    subroutine single_interpolation(x,x_inf,x_sup,val_inf,val_sup,val)
        !Input
        double precision::x,x_sup,x_inf,val_inf,val_sup
        !Output
        double precision::val
        val=val_inf+(val_sup-val_inf)*(x-x_inf)/(x_sup-x_inf)        
    end subroutine single_interpolation
    
    

    ! find the zero of single valued function using the secant method
    ! F = subroutine F(err,X) where err is the error to zero, X is the independent variable
    ! B = minimum value and also output variable (X s.t. F(X)=0)
    ! C = maximum value (must be variable)
    ! R = guessed starting point (must be variable)
    ! re = relative tolerance
    ! AE = absolute tolerance
    ! IFLAG = 1 all fine, 2 acceptable, >= 3 some error occurred
    ! NOTE: maximum number of iterations is internally set to 500
    !       when exceeded, IFLAG=5
    recursive subroutine fzero (F, B, C, R, re, AE, IFLAG)
        double precision A,ACBS,ACMB,AE,AW,B,C,CMB,ER,FA,FB,FC,FX,FZ,P,Q,R,re,RW,T,TOL,Z
        integer IC,IFLAG,KOUNT
        external F
        ER = 1d-9
        !   Initialize.
        Z = R
        if (R .LE. min(B,C)  .OR.  R .GE. max(B,C)) Z = C
        RW = max(re,ER)
        AW = max(AE,0d0)
        IC = 0
        T = Z
        call F(FZ,T)
        FC = FZ
        T = B
        call F(FB,T)
        KOUNT = 2
        if (sign(1.0D0,FZ) .EQ. sign(1.0D0,FB)) go to 1
        C = Z
        go to 2
        1 if (Z .EQ. C) go to 2
        T = C
        call F(FC,T)
        KOUNT = 3
        if (sign(1.0D0,FZ) .EQ. sign(1.0D0,FC)) go to 2
        B = Z
        FB = FZ
        2 A = C
        FA = FC
        ACBS = abs(B-C)
        FX = max(abs(FB),abs(FC))
        3 if (abs(FC) .GE. abs(FB)) go to 4
        !   Perform interchange.
        A = B
        FA = FB
        B = C
        FB = FC
        C = A
        FC = FA
        4 CMB = 0.5D0*(C-B)
        ACMB = abs(CMB)
        TOL = RW*abs(B) + AW
        !   Test stopping criterion and function count.
        if (ACMB .LE. TOL) go to 10
        if (FB .EQ. 0.D0) go to 11
        if (KOUNT .GE. 500) go to 14
        !   Calculate new iterate implicitly as B+P/Q, where we arrange
        !   P .GE. 0.  The implicit form is used to prevent overflow.
        P = (B-A)*FB
        Q = FA - FB
        if (P .GE. 0.D0) go to 5
        P = -P
        Q = -Q
        !   Update A and check for satisfactory reduction in the size of the
        !   bracketing interval.  If not, perform bisection.
        5 A = B
        FA = FB
        IC = IC + 1
        if (IC .LT. 4) go to 6
        if (8.0D0*ACMB .GE. ACBS) go to 8
        IC = 0
        ACBS = ACMB
        !   Test for too small a change.
        6 if (P .GT. abs(Q)*TOL) go to 7
        !   Increment by TOLerance.
        B = B + sign(TOL,CMB)
        go to 9
        !   Root ought to be between B and (C+B)/2.
        7 if (P .GE. CMB*Q) go to 8
        !   Use secant rule.
        B = B + P/Q
        go to 9
        !   Use bisection (C+B)/2.
        8 B = B + CMB
        !   Have completed computation for new iterate B.
        9 T = B
        call F(FB,T)
        KOUNT = KOUNT + 1
        !   Decide whether next step is interpolation or extrapolation.
        if (sign(1.0D0,FB) .NE. sign(1.0D0,FC)) go to 3
        C = A
        FC = FA
        go to 3
        !   Finished.  Process results for proper setting of IFLAG.
        10 if (sign(1.0D0,FB) .EQ. sign(1.0D0,FC)) go to 13
        if (abs(FB) .GT. FX) go to 12
        IFLAG = 1
        return
        11 IFLAG = 2
        return
        12 IFLAG = 3
        return
        13 IFLAG = 4
        return
        14 IFLAG = 5
        return
    end

    ! <TODO> eliminate fzero2, declaring fzero recursive
    subroutine fzero2 (F, B, C, R, re, AE, IFLAG)
        double precision A,ACBS,ACMB,AE,AW,B,C,CMB,ER,FA,FB,FC,FX,FZ,P,Q,R,re,RW,T,TOL,Z
        integer IC,IFLAG,KOUNT
        external F
        ER = 1d-9
        !   Initialize.
        Z = R
        if (R .LE. min(B,C)  .OR.  R .GE. max(B,C)) Z = C
        RW = max(re,ER)
        AW = max(AE,0d0)
        IC = 0
        T = Z
        call F(FZ,T)
        FC = FZ
        T = B
        call F(FB,T)
        KOUNT = 2
        if (sign(1.0D0,FZ) .EQ. sign(1.0D0,FB)) go to 1
        C = Z
        go to 2
        1 if (Z .EQ. C) go to 2
        T = C
        call F(FC,T)
        KOUNT = 3
        if (sign(1.0D0,FZ) .EQ. sign(1.0D0,FC)) go to 2
        B = Z
        FB = FZ
        2 A = C
        FA = FC
        ACBS = abs(B-C)
        FX = max(abs(FB),abs(FC))
        3 if (abs(FC) .GE. abs(FB)) go to 4
        !   Perform interchange.
        A = B
        FA = FB
        B = C
        FB = FC
        C = A
        FC = FA
        4 CMB = 0.5D0*(C-B)
        ACMB = abs(CMB)
        TOL = RW*abs(B) + AW
        !   Test stopping criterion and function count.
        if (ACMB .LE. TOL) go to 10
        if (FB .EQ. 0.D0) go to 11
        if (KOUNT .GE. 500) go to 14
        !   Calculate new iterate implicitly as B+P/Q, where we arrange
        !   P .GE. 0.  The implicit form is used to prevent overflow.
        P = (B-A)*FB
        Q = FA - FB
        if (P .GE. 0.D0) go to 5
        P = -P
        Q = -Q
        !   Update A and check for satisfactory reduction in the size of the
        !   bracketing interval.  If not, perform bisection.
        5 A = B
        FA = FB
        IC = IC + 1
        if (IC .LT. 4) go to 6
        if (8.0D0*ACMB .GE. ACBS) go to 8
        IC = 0
        ACBS = ACMB
        !   Test for too small a change.
        6 if (P .GT. abs(Q)*TOL) go to 7
        !   Increment by TOLerance.
        B = B + sign(TOL,CMB)
        go to 9
        !   Root ought to be between B and (C+B)/2.
        7 if (P .GE. CMB*Q) go to 8
        !   Use secant rule.
        B = B + P/Q
        go to 9
        !   Use bisection (C+B)/2.
        8 B = B + CMB
        !   Have completed computation for new iterate B.
        9 T = B
        call F(FB,T)
        KOUNT = KOUNT + 1
        !   Decide whether next step is interpolation or extrapolation.
        if (sign(1.0D0,FB) .NE. sign(1.0D0,FC)) go to 3
        C = A
        FC = FA
        go to 3
        !   Finished.  Process results for proper setting of IFLAG.
        10 if (sign(1.0D0,FB) .EQ. sign(1.0D0,FC)) go to 13
        if (abs(FB) .GT. FX) go to 12
        IFLAG = 1
        return
        11 IFLAG = 2
        return
        12 IFLAG = 3
        return
        13 IFLAG = 4
        return
        14 IFLAG = 5
        return
    end

    ! guess x0 for a 2x2 equations system
    function guess2x2(F,xBounds,nSteps) result(x0)
        double precision,dimension(2,1)::x0,F0,x1,F1
        double precision,dimension(2,2)::xBounds
        integer,dimension(2)::nSteps
        integer::i,j
        external F
        F0(1,1)=Inf
        F0(2,1)=Inf
        do i=1,nSteps(1)
        	do j=1,nSteps(2)
        		x1(1,1)=xBounds(1,1)+(xBounds(1,2)-xBounds(1,1))/(nSteps(1)-1)*(i-1)
        		x1(2,1)=xBounds(2,1)+(xBounds(2,2)-xBounds(2,1))/(nSteps(2)-1)*(j-1)
            	call F(x1,F1)
	            if (maxval(abs(F1))<maxval(abs(F0))) then
                	F0=F1
                	x0=x1
                end if
            end do
        end do
	end function

    ! Newton-Raphson method 2x2
    function NR2x2(x0,F,xBounds,eps,info) result(x)
        double precision,dimension(2)::x0,x,F0,aerr
        double precision,dimension(2,2)::xBounds
        double precision,dimension(2)::eps,h
        integer::info,i,j
        external F
        info=1
        h=(/eps(1)*(xBounds(1,2)-xBounds(1,1)),eps(2)*(xBounds(2,2)-xBounds(2,1))/)
        do i=1,500
            call F(x0,F0)
            x=x0-0.5*matmul(inv2x2(Jac(F,F0,x0)),F0)
            aerr=abs(x-x0)
            if (aerr(1)<h(1).and.aerr(2)<h(2)) then
                info=0
                exit
            else
	            x(1)=max(min(x(1),xBounds(1,2)),xBounds(1,1))
    	        x(2)=max(min(x(2),xBounds(2,2)),xBounds(2,1))
                x0=x
            end if
        end do
        if (i>500) then
            print *, 'Error: max iterations exceeded'
            write(*,*) aerr(1), ' ', aerr(2)
        end if
        contains
        function Jac(F,F0,x) result(M)
            double precision,dimension(2,2)::M
            double precision,dimension(2)::x,x1,x2,F1,F2,F0,dF1,dF2
            external F
            x1(1)=x(1)+h(1)
            x1(2)=x(2)
            x2(1)=x(1)
            x2(2)=x(2)+h(2)
            call F(x1,F1)
            call F(x2,F2)
            dF1=(F1-F0)/h(1)
            dF2=(F2-F0)/h(2)
            M(1,1)=dF1(1)
            M(2,1)=dF1(2)
            M(1,2)=dF2(1)
            M(2,2)=dF2(2)
        end function
    end function

    function inv2x2(M) result(invM)
        double precision,dimension(2,2)::M,invM
        double precision::delta
        delta=M(1,1)*M(2,2)-M(1,2)*M(2,1)
        invM(1,1)=M(2,2)
        invM(1,2)=-M(1,2)
        invM(2,1)=-M(2,1)
        invM(2,2)=M(1,1)
        if (delta/=0) then
            invM=invM/delta
        else
            invM=0
            print *, 'Error: singular matrix'
            !invM=invM*abs(M(1,1)*M(2,2))
        end if
    end function

    ! Newton-Raphson method NxN symmetric Jacobian
    function NRNxNsym(x0,F,N,rtol,atol,maxi,info,xmin) result(x)
    	integer::N
        double precision,dimension(N)::x00,x0,x,Dx,F0
        double precision,dimension(N,N)::Jac
        double precision::aerr,rerr,rden,rtol,atol,aerr0,r,rmax,rmin,xmin
        integer,intent(out)::info
        integer::maxi,i,j,try,info2
        external F
        info=-1
        r=0.1 ! initial radius
        rmin=0.0001 ! minimum radius
        rmax=0.5 ! maximum radius
        do try=1,2
          if (info==-1) then
	        do i=1,maxi
	            call F(N,x0,F0,Jac)
				aerr=norm2(F0)/N
                !write(*,*) 'iter',i,'err',aerr
	            !if (mod(i,2)==0) then
	            !	x=0.5*(x00+x0)
	            !else
		            call linsolvesym(Jac,F0,Dx,N,info2)
	            	x=x0-r*Dx
	            !end if	
	            if (i>1) then
	            	! update radius
	            	if (aerr0>aerr) then
	            		if (try==1)	r=min(1.2*r,rmax)
	            	else
	            		r=max(0.8*r,rmin)
	            	end if		
	            end if
	            rden=maxval(x)-minval(x)
	            if (rden>0) then
	            	rerr=maxval(abs(x-x0))/rden
	            else
	            	rerr=0d0
	            end if		
	            if (aerr<atol .and. rerr<rtol) then
	                info=0
	                exit
	            else
	            	aerr0=aerr
	            	x00=x0
	                x0=x
	            end if
	            !write(*,*) i,r,aerr,rerr
	        end do
	      end if  
	    end do    
!        if (i>maxi .and. try>2) then
        if (info==-1) then
            print *, 'Error: max iterations exceeded'
            write(*,*) 'absolute error', aerr
            write(*,*) 'relative error', rerr
            write(*,*) 'radius',r
            write(*,*) 'vector x'
            write(*,'(f32.16)') x
        end if
    end function

    ! smooth decay from 1 to 0 between x=0 and x=1
    function decay(x)
    	double precision::decay,x,t
    	t=between(x,0d0,1d0)
     	decay=2d0*t**3-3d0*t**2+1d0
    end function

    ! limit x between x0 and x1>x0
    function between(x,x0,x1)
    	double precision::between,x,x0,x1
    	between=max(min(x,x1),x0)
    end function

    ! linear interpolation for x in (0,1) and y in (y0,y1)
    function linear(x,y0,y1)
    	double precision::linear,x,y0,y1
    	linear=y0+(y1-y0)*x
	end function


!*****************************************************************************************
!  Multidimensional linear interpolation/extrapolation.
!*****************************************************************************************
!  Uses repeated linear interpolation to evaluate
!  functions \(f(x), f(x,y), f(x,y,z), f(x,y,z,q), f(x,y,z,q,r), f(x,y,z,q,r,s) \)
!  which have been tabulated at the nodes of an n-dimensional rectangular grid.
!  If any coordinate \( (x_i, y_i, ...) \) lies outside the range of the corresponding
!  variable, then extrapolation is performed using the two nearest points.
!>
!  Finalizer for a [[linear_interp_1d]] type.

    pure elemental subroutine finalize_1d(me)

    implicit none

    type(linear_interp_1d),intent(inout) :: me
    call me%destroy()

    end subroutine finalize_1d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Finalizer for a [[linear_interp_2d]] type.

    pure elemental subroutine finalize_2d(me)

    implicit none

    type(linear_interp_2d),intent(inout) :: me
    call me%destroy()

    end subroutine finalize_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Finalizer for a [[linear_interp_3d]] type.

    pure elemental subroutine finalize_3d(me)

    implicit none

    type(linear_interp_3d),intent(inout) :: me
    call me%destroy()

    end subroutine finalize_3d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Finalizer for a [[linear_interp_4d]] type.

    pure elemental subroutine finalize_4d(me)

    implicit none

    type(linear_interp_4d),intent(inout) :: me
    call me%destroy()

    end subroutine finalize_4d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Finalizer for a [[linear_interp_5d]] type.

    pure elemental subroutine finalize_5d(me)

    implicit none

    type(linear_interp_5d),intent(inout) :: me
    call me%destroy()

    end subroutine finalize_5d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Finalizer for a [[linear_interp_6d]] type.

    pure elemental subroutine finalize_6d(me)

    implicit none

    type(linear_interp_6d),intent(inout) :: me
    call me%destroy()

    end subroutine finalize_6d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for a [[linear_interp_1d]] class.

    pure elemental subroutine destroy_1d(me)

    implicit none

    class(linear_interp_1d),intent(inout) :: me

    if (allocated(me%f)) deallocate(me%f)
    if (allocated(me%x)) deallocate(me%x)
    me%ilox = 1

    end subroutine destroy_1d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for a [[linear_interp_2d]] class.

    pure elemental subroutine destroy_2d(me)

    implicit none

    class(linear_interp_2d),intent(inout) :: me

    if (allocated(me%f)) deallocate(me%f)
    if (allocated(me%x)) deallocate(me%x)
    if (allocated(me%y)) deallocate(me%y)
    me%ilox = 1
    me%iloy = 1

    end subroutine destroy_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for a [[linear_interp_3d]] class.

    pure elemental subroutine destroy_3d(me)

    implicit none

    class(linear_interp_3d),intent(inout) :: me

    if (allocated(me%f)) deallocate(me%f)
    if (allocated(me%x)) deallocate(me%x)
    if (allocated(me%y)) deallocate(me%y)
    if (allocated(me%z)) deallocate(me%z)
    me%ilox = 1
    me%iloy = 1
    me%iloz = 1

    end subroutine destroy_3d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for a [[linear_interp_4d]] class.

    pure elemental subroutine destroy_4d(me)

    implicit none

    class(linear_interp_4d),intent(inout) :: me

    if (allocated(me%f)) deallocate(me%f)
    if (allocated(me%x)) deallocate(me%x)
    if (allocated(me%y)) deallocate(me%y)
    if (allocated(me%z)) deallocate(me%z)
    if (allocated(me%q)) deallocate(me%q)
    me%ilox = 1
    me%iloy = 1
    me%iloz = 1
    me%iloq = 1

    end subroutine destroy_4d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for a [[linear_interp_5d]] class.

    pure elemental subroutine destroy_5d(me)

    implicit none

    class(linear_interp_5d),intent(inout) :: me

    if (allocated(me%f)) deallocate(me%f)
    if (allocated(me%x)) deallocate(me%x)
    if (allocated(me%y)) deallocate(me%y)
    if (allocated(me%z)) deallocate(me%z)
    if (allocated(me%q)) deallocate(me%q)
    if (allocated(me%r)) deallocate(me%r)
    me%ilox = 1
    me%iloy = 1
    me%iloz = 1
    me%iloq = 1
    me%ilor = 1

    end subroutine destroy_5d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for a [[linear_interp_6d]] class.

    pure elemental subroutine destroy_6d(me)

    implicit none

    class(linear_interp_6d),intent(inout) :: me

    if (allocated(me%f)) deallocate(me%f)
    if (allocated(me%x)) deallocate(me%x)
    if (allocated(me%y)) deallocate(me%y)
    if (allocated(me%z)) deallocate(me%z)
    if (allocated(me%q)) deallocate(me%q)
    if (allocated(me%r)) deallocate(me%r)
    if (allocated(me%s)) deallocate(me%s)
    me%ilox = 1
    me%iloy = 1
    me%iloz = 1
    me%iloq = 1
    me%ilor = 1
    me%ilos = 1

    end subroutine destroy_6d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[linear_interp_1d]] class.

    pure subroutine initialize_1d(me,x,f,istat)

    implicit none

    class(linear_interp_1d),intent(inout) :: me
    double precision,dimension(:),intent(in)      :: x
    double precision,dimension(:),intent(in)      :: f
    integer,intent(out)                   :: istat  !! `0`  : no problems,
                                                    !! `1`  : `x` is not strictly increasing,
                                                    !! `10` : `x` is not equal to size(f,1).

    call me%destroy()

    istat = 0

    if (istat==0 .and. size(x)/=size(f,1)) istat = 10

    if (istat==0) then
        call check_inputs(x=x,ierr=istat)
        if (istat==0) then
            allocate(me%f(size(x))); me%f = f
            allocate(me%x(size(x))); me%x = x
        end if
    end if

    end subroutine initialize_1d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[linear_interp_2d]] class.

    pure subroutine initialize_2d(me,x,y,f,istat)

    implicit none

    class(linear_interp_2d),intent(inout) :: me
    double precision,dimension(:),intent(in)      :: x
    double precision,dimension(:),intent(in)      :: y
    double precision,dimension(:,:),intent(in)    :: f
    integer,intent(out)                   :: istat  !! `0`  : no problems,
                                                    !! `1`  : `x` is not strictly increasing,
                                                    !! `2`  : `y` is not strictly increasing,
                                                    !! `10` : `x` is not equal to size(f,1),
                                                    !! `20` : `y` is not equal to size(f,2).

    call me%destroy()

    istat = 0

    if (istat==0 .and. size(x)/=size(f,1)) istat = 10
    if (istat==0 .and. size(y)/=size(f,2)) istat = 20

    if (istat==0) then
        call check_inputs(x=x,y=y,ierr=istat)
        if (istat==0) then
            allocate(me%f(size(x),size(y))); me%f = f
            allocate(me%x(size(x))); me%x = x
            allocate(me%y(size(y))); me%y = y
        end if
    end if

    end subroutine initialize_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[linear_interp_3d]] class.

    pure subroutine initialize_3d(me,x,y,z,f,istat)

    implicit none

    class(linear_interp_3d),intent(inout) :: me
    double precision,dimension(:),intent(in)      :: x
    double precision,dimension(:),intent(in)      :: y
    double precision,dimension(:),intent(in)      :: z
    double precision,dimension(:,:,:),intent(in)  :: f
    integer,intent(out)                   :: istat  !! `0`  : no problems,
                                                    !! `1`  : `x` is not strictly increasing,
                                                    !! `2`  : `y` is not strictly increasing,
                                                    !! `3`  : `z` is not strictly increasing,
                                                    !! `10` : `x` is not equal to size(f,1),
                                                    !! `20` : `y` is not equal to size(f,2),
                                                    !! `30` : `z` is not equal to size(f,3).

    call me%destroy()

    istat = 0

    if (istat==0 .and. size(x)/=size(f,1)) istat = 10
    if (istat==0 .and. size(y)/=size(f,2)) istat = 20
    if (istat==0 .and. size(z)/=size(f,3)) istat = 30

    if (istat==0) then
        call check_inputs(x=x,y=y,z=z,ierr=istat)
        if (istat==0) then
            allocate(me%f(size(x),size(y),size(z))); me%f = f
            allocate(me%x(size(x))); me%x = x
            allocate(me%y(size(y))); me%y = y
            allocate(me%z(size(z))); me%z = z
        end if
    end if

    end subroutine initialize_3d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[linear_interp_4d]] class.

    pure subroutine initialize_4d(me,x,y,z,q,f,istat)

    implicit none

    class(linear_interp_4d),intent(inout)  :: me
    double precision,dimension(:),intent(in)       :: x
    double precision,dimension(:),intent(in)       :: y
    double precision,dimension(:),intent(in)       :: z
    double precision,dimension(:),intent(in)       :: q
    double precision,dimension(:,:,:,:),intent(in) :: f
    integer,intent(out)                    :: istat !! `0`  : no problems,
                                                    !! `1`  : `x` is not strictly increasing,
                                                    !! `2`  : `y` is not strictly increasing,
                                                    !! `3`  : `z` is not strictly increasing,
                                                    !! `4`  : `q` is not strictly increasing,
                                                    !! `10` : `x` is not equal to size(f,1),
                                                    !! `20` : `y` is not equal to size(f,2),
                                                    !! `30` : `z` is not equal to size(f,3),
                                                    !! `40` : `q` is not equal to size(f,4).

    call me%destroy()

    istat = 0

    if (istat==0 .and. size(x)/=size(f,1)) istat = 10
    if (istat==0 .and. size(y)/=size(f,2)) istat = 20
    if (istat==0 .and. size(z)/=size(f,3)) istat = 30
    if (istat==0 .and. size(q)/=size(f,4)) istat = 40

    if (istat==0) then
        call check_inputs(x=x,y=y,z=z,q=q,ierr=istat)
        if (istat==0) then
            allocate(me%f(size(x),size(y),size(z),size(q))); me%f = f
            allocate(me%x(size(x))); me%x = x
            allocate(me%y(size(y))); me%y = y
            allocate(me%z(size(z))); me%z = z
            allocate(me%q(size(q))); me%q = q
        end if
    end if

    end subroutine initialize_4d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[linear_interp_5d]] class.

    pure subroutine initialize_5d(me,x,y,z,q,r,f,istat)

    implicit none

    class(linear_interp_5d),intent(inout)    :: me
    double precision,dimension(:),intent(in)         :: x
    double precision,dimension(:),intent(in)         :: y
    double precision,dimension(:),intent(in)         :: z
    double precision,dimension(:),intent(in)         :: q
    double precision,dimension(:),intent(in)         :: r
    double precision,dimension(:,:,:,:,:),intent(in) :: f
    integer,intent(out)                      :: istat   !! `0`  : no problems,
                                                        !! `1`  : `x` is not strictly increasing,
                                                        !! `2`  : `y` is not strictly increasing,
                                                        !! `3`  : `z` is not strictly increasing,
                                                        !! `4`  : `q` is not strictly increasing,
                                                        !! `5`  : `r` is not strictly increasing,
                                                        !! `10` : `x` is not equal to size(f,1),
                                                        !! `20` : `y` is not equal to size(f,2),
                                                        !! `30` : `z` is not equal to size(f,3),
                                                        !! `40` : `q` is not equal to size(f,4),
                                                        !! `50` : `r` is not equal to size(f,5).

    call me%destroy()

    istat = 0

    if (istat==0 .and. size(x)/=size(f,1)) istat = 10
    if (istat==0 .and. size(y)/=size(f,2)) istat = 20
    if (istat==0 .and. size(z)/=size(f,3)) istat = 30
    if (istat==0 .and. size(q)/=size(f,4)) istat = 40
    if (istat==0 .and. size(r)/=size(f,5)) istat = 50

    if (istat==0) then
        call check_inputs(x=x,y=y,z=z,q=q,r=r,ierr=istat)
        if (istat==0) then
            allocate(me%f(size(x),size(y),size(z),size(q),size(r))); me%f = f
            allocate(me%x(size(x))); me%x = x
            allocate(me%y(size(y))); me%y = y
            allocate(me%z(size(z))); me%z = z
            allocate(me%q(size(q))); me%q = q
            allocate(me%r(size(r))); me%r = r
        end if
    end if

    end subroutine initialize_5d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[linear_interp_6d]] class.

    pure subroutine initialize_6d(me,x,y,z,q,r,s,f,istat)

    implicit none

    class(linear_interp_6d),intent(inout)      :: me
    double precision,dimension(:),intent(in)           :: x
    double precision,dimension(:),intent(in)           :: y
    double precision,dimension(:),intent(in)           :: z
    double precision,dimension(:),intent(in)           :: q
    double precision,dimension(:),intent(in)           :: r
    double precision,dimension(:),intent(in)           :: s
    double precision,dimension(:,:,:,:,:,:),intent(in) :: f
    integer,intent(out)                        :: istat !! `0`  : no problems,
                                                        !! `1`  : `x` is not strictly increasing,
                                                        !! `2`  : `y` is not strictly increasing,
                                                        !! `3`  : `z` is not strictly increasing,
                                                        !! `4`  : `q` is not strictly increasing,
                                                        !! `5`  : `r` is not strictly increasing,
                                                        !! `6`  : `s` is not strictly increasing,
                                                        !! `10` : `x` is not equal to size(f,1),
                                                        !! `20` : `y` is not equal to size(f,2),
                                                        !! `30` : `z` is not equal to size(f,3),
                                                        !! `40` : `q` is not equal to size(f,4),
                                                        !! `50` : `r` is not equal to size(f,5),
                                                        !! `60` : `s` is not equal to size(f,6).

    call me%destroy()

    istat = 0

    if (istat==0 .and. size(x)/=size(f,1)) istat = 10
    if (istat==0 .and. size(y)/=size(f,2)) istat = 20
    if (istat==0 .and. size(z)/=size(f,3)) istat = 30
    if (istat==0 .and. size(q)/=size(f,4)) istat = 40
    if (istat==0 .and. size(r)/=size(f,5)) istat = 50
    if (istat==0 .and. size(s)/=size(f,6)) istat = 60

    if (istat==0) then
        call check_inputs(x=x,y=y,z=z,q=q,r=r,s=s,ierr=istat)
        if (istat==0) then
            allocate(me%f(size(x),size(y),size(z),size(q),size(r),size(s))); me%f = f
            allocate(me%x(size(x))); me%x = x
            allocate(me%y(size(y))); me%y = y
            allocate(me%z(size(z))); me%z = z
            allocate(me%q(size(q))); me%q = q
            allocate(me%r(size(r))); me%r = r
            allocate(me%s(size(s))); me%s = s
        end if
    end if

    end subroutine initialize_6d
!*****************************************************************************************

!*****************************************************************************************
!>
!  1D linear interpolation routine.

    pure subroutine interp_1d(me,x,fx)

    implicit none

    class(linear_interp_1d),intent(inout) :: me
    double precision,intent(in)                   :: x
    double precision,intent(out)                  :: fx  !! Interpolated \( f(x) \)

    integer,dimension(2) :: ix
    double precision :: p1
    double precision :: q1
    integer :: mflag

    call dintrv(me%x,x,me%ilox,ix(1),ix(2),mflag)

    q1 = (x-me%x(ix(1)))/(me%x(ix(2))-me%x(ix(1)))
    p1 = one-q1

    fx = p1*me%f(ix(1)) + q1*me%f(ix(2))

    end subroutine interp_1d
!*****************************************************************************************

!*****************************************************************************************
!>
!  2D linear interpolation routine.

    pure subroutine interp_2d(me,x,y,fxy)

    implicit none

    class(linear_interp_2d),intent(inout) :: me
    double precision,intent(in)                   :: x
    double precision,intent(in)                   :: y
    double precision,intent(out)                  :: fxy  !! Interpolated \( f(x,y) \)

    integer,dimension(2) :: ix, iy
    double precision :: p1, p2
    double precision :: q1, q2
    integer :: mflag
    double precision :: fx1, fx2

    call dintrv(me%x,x,me%ilox,ix(1),ix(2),mflag)
    call dintrv(me%y,y,me%iloy,iy(1),iy(2),mflag)

    q1 = (x-me%x(ix(1)))/(me%x(ix(2))-me%x(ix(1)))
    q2 = (y-me%y(iy(1)))/(me%y(iy(2))-me%y(iy(1)))
    p1 = one-q1
    p2 = one-q2

    fx1 = p1*me%f(ix(1),iy(1)) + q1*me%f(ix(2),iy(1))
    fx2 = p1*me%f(ix(1),iy(2)) + q1*me%f(ix(2),iy(2))

    fxy = p2*( fx1 ) + q2*( fx2 )

    end subroutine interp_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  3D linear interpolation routine.

    pure subroutine interp_3d(me,x,y,z,fxyz)

    implicit none

    class(linear_interp_3d),intent(inout) :: me
    double precision,intent(in)                   :: x
    double precision,intent(in)                   :: y
    double precision,intent(in)                   :: z
    double precision,intent(out)                  :: fxyz  !! Interpolated \( f(x,y,z) \)

    integer,dimension(2) :: ix, iy, iz
    double precision :: p1, p2, p3
    double precision :: q1, q2, q3
    integer :: mflag
    double precision :: fx11, fx21, fx12, fx22, fxy1, fxy2

    call dintrv(me%x,x,me%ilox,ix(1),ix(2),mflag)
    call dintrv(me%y,y,me%iloy,iy(1),iy(2),mflag)
    call dintrv(me%z,z,me%iloz,iz(1),iz(2),mflag)

    q1 = (x-me%x(ix(1)))/(me%x(ix(2))-me%x(ix(1)))
    q2 = (y-me%y(iy(1)))/(me%y(iy(2))-me%y(iy(1)))
    q3 = (z-me%z(iz(1)))/(me%z(iz(2))-me%z(iz(1)))
    p1 = one-q1
    p2 = one-q2
    p3 = one-q3

    fx11 = p1*me%f(ix(1),iy(1),iz(1)) + q1*me%f(ix(2),iy(1),iz(1))
    fx21 = p1*me%f(ix(1),iy(2),iz(1)) + q1*me%f(ix(2),iy(2),iz(1))
    fx12 = p1*me%f(ix(1),iy(1),iz(2)) + q1*me%f(ix(2),iy(1),iz(2))
    fx22 = p1*me%f(ix(1),iy(2),iz(2)) + q1*me%f(ix(2),iy(2),iz(2))
    fxy1 = p2*( fx11 ) + q2*( fx21 )
    fxy2 = p2*( fx12 ) + q2*( fx22 )

    fxyz = p3*( fxy1 ) + q3*( fxy2 )

    end subroutine interp_3d
!*****************************************************************************************

!*****************************************************************************************
!>
!  4D linear interpolation routine.

    pure subroutine interp_4d(me,x,y,z,q,fxyzq)

    implicit none

    class(linear_interp_4d),intent(inout) :: me
    double precision,intent(in)                   :: x
    double precision,intent(in)                   :: y
    double precision,intent(in)                   :: z
    double precision,intent(in)                   :: q
    double precision,intent(out)                  :: fxyzq  !! Interpolated \( f(x,y,z,q) \)

    integer,dimension(2) :: ix, iy, iz, iq
    double precision :: p1, p2, p3, p4
    double precision :: q1, q2, q3, q4
    integer :: mflag
    double precision :: fx111,fx211,fx121,fx221,fxy11,fxy21,fxyz1,&
                fx112,fx212,fx122,fx222,fxy12,fxy22,fxyz2

    call dintrv(me%x,x,me%ilox,ix(1),ix(2),mflag)
    call dintrv(me%y,y,me%iloy,iy(1),iy(2),mflag)
    call dintrv(me%z,z,me%iloz,iz(1),iz(2),mflag)
    call dintrv(me%q,q,me%iloq,iq(1),iq(2),mflag)

    q1 = (x-me%x(ix(1)))/(me%x(ix(2))-me%x(ix(1)))
    q2 = (y-me%y(iy(1)))/(me%y(iy(2))-me%y(iy(1)))
    q3 = (z-me%z(iz(1)))/(me%z(iz(2))-me%z(iz(1)))
    q4 = (q-me%q(iq(1)))/(me%q(iq(2))-me%q(iq(1)))
    p1 = one-q1
    p2 = one-q2
    p3 = one-q3
    p4 = one-q4

    fx111 = p1*me%f(ix(1),iy(1),iz(1),iq(1)) + q1*me%f(ix(2),iy(1),iz(1),iq(1))
    fx211 = p1*me%f(ix(1),iy(2),iz(1),iq(1)) + q1*me%f(ix(2),iy(2),iz(1),iq(1))
    fx121 = p1*me%f(ix(1),iy(1),iz(2),iq(1)) + q1*me%f(ix(2),iy(1),iz(2),iq(1))
    fx221 = p1*me%f(ix(1),iy(2),iz(2),iq(1)) + q1*me%f(ix(2),iy(2),iz(2),iq(1))
    fx112 = p1*me%f(ix(1),iy(1),iz(1),iq(2)) + q1*me%f(ix(2),iy(1),iz(1),iq(2))
    fx212 = p1*me%f(ix(1),iy(2),iz(1),iq(2)) + q1*me%f(ix(2),iy(2),iz(1),iq(2))
    fx122 = p1*me%f(ix(1),iy(1),iz(2),iq(2)) + q1*me%f(ix(2),iy(1),iz(2),iq(2))
    fx222 = p1*me%f(ix(1),iy(2),iz(2),iq(2)) + q1*me%f(ix(2),iy(2),iz(2),iq(2))

    fxy11 = p2*fx111 + q2*fx211
    fxy21 = p2*fx121 + q2*fx221
    fxy12 = p2*fx112 + q2*fx212
    fxy22 = p2*fx122 + q2*fx222

    fxyz1 = p3*fxy11 + q3*fxy21
    fxyz2 = p3*fxy12 + q3*fxy22

    fxyzq = p4*fxyz1 + q4*fxyz2

    end subroutine interp_4d
!*****************************************************************************************

!*****************************************************************************************
!>
!  5D linear interpolation routine.

    pure subroutine interp_5d(me,x,y,z,q,r,fxyzqr)

    implicit none

    class(linear_interp_5d),intent(inout) :: me
    double precision,intent(in)                   :: x
    double precision,intent(in)                   :: y
    double precision,intent(in)                   :: z
    double precision,intent(in)                   :: q
    double precision,intent(in)                   :: r
    double precision,intent(out)                  :: fxyzqr  !! Interpolated \( f(x,y,z,q,r) \)

    integer,dimension(2) :: ix, iy, iz, iq, ir
    double precision :: p1, p2, p3, p4, p5
    double precision :: q1, q2, q3, q4, q5
    integer :: mflag
    double precision :: fx1111, fx2111, fx1211, fx2211, fx1121, fx2121, fx1221, fx2221, &
                fxy111, fxy211, fxy121, fxy221, fxyz11, fxyz21, fxyzq1, fx1112, &
                fx2112, fx1212, fx2212, fx1122, fx2122, fx1222, fx2222, fxy112, &
                fxy212, fxy122, fxy222, fxyz12, fxyz22, fxyzq2

    call dintrv(me%x,x,me%ilox,ix(1),ix(2),mflag)
    call dintrv(me%y,y,me%iloy,iy(1),iy(2),mflag)
    call dintrv(me%z,z,me%iloz,iz(1),iz(2),mflag)
    call dintrv(me%q,q,me%iloq,iq(1),iq(2),mflag)
    call dintrv(me%r,r,me%ilor,ir(1),ir(2),mflag)

    q1 = (x-me%x(ix(1)))/(me%x(ix(2))-me%x(ix(1)))
    q2 = (y-me%y(iy(1)))/(me%y(iy(2))-me%y(iy(1)))
    q3 = (z-me%z(iz(1)))/(me%z(iz(2))-me%z(iz(1)))
    q4 = (q-me%q(iq(1)))/(me%q(iq(2))-me%q(iq(1)))
    q5 = (r-me%r(ir(1)))/(me%r(ir(2))-me%r(ir(1)))
    p1 = one-q1
    p2 = one-q2
    p3 = one-q3
    p4 = one-q4
    p5 = one-q5

    fx1111 = p1*me%f(ix(1),iy(1),iz(1),iq(1),ir(1)) + q1*me%f(ix(2),iy(1),iz(1),iq(1),ir(1))
    fx2111 = p1*me%f(ix(1),iy(2),iz(1),iq(1),ir(1)) + q1*me%f(ix(2),iy(2),iz(1),iq(1),ir(1))
    fx1211 = p1*me%f(ix(1),iy(1),iz(2),iq(1),ir(1)) + q1*me%f(ix(2),iy(1),iz(2),iq(1),ir(1))
    fx2211 = p1*me%f(ix(1),iy(2),iz(2),iq(1),ir(1)) + q1*me%f(ix(2),iy(2),iz(2),iq(1),ir(1))
    fx1121 = p1*me%f(ix(1),iy(1),iz(1),iq(2),ir(1)) + q1*me%f(ix(2),iy(1),iz(1),iq(2),ir(1))
    fx2121 = p1*me%f(ix(1),iy(2),iz(1),iq(2),ir(1)) + q1*me%f(ix(2),iy(2),iz(1),iq(2),ir(1))
    fx1221 = p1*me%f(ix(1),iy(1),iz(2),iq(2),ir(1)) + q1*me%f(ix(2),iy(1),iz(2),iq(2),ir(1))
    fx2221 = p1*me%f(ix(1),iy(2),iz(2),iq(2),ir(1)) + q1*me%f(ix(2),iy(2),iz(2),iq(2),ir(1))
    fx1112 = p1*me%f(ix(1),iy(1),iz(1),iq(1),ir(2)) + q1*me%f(ix(2),iy(1),iz(1),iq(1),ir(2))
    fx2112 = p1*me%f(ix(1),iy(2),iz(1),iq(1),ir(2)) + q1*me%f(ix(2),iy(2),iz(1),iq(1),ir(2))
    fx1212 = p1*me%f(ix(1),iy(1),iz(2),iq(1),ir(2)) + q1*me%f(ix(2),iy(1),iz(2),iq(1),ir(2))
    fx2212 = p1*me%f(ix(1),iy(2),iz(2),iq(1),ir(2)) + q1*me%f(ix(2),iy(2),iz(2),iq(1),ir(2))
    fx1122 = p1*me%f(ix(1),iy(1),iz(1),iq(2),ir(2)) + q1*me%f(ix(2),iy(1),iz(1),iq(2),ir(2))
    fx2122 = p1*me%f(ix(1),iy(2),iz(1),iq(2),ir(2)) + q1*me%f(ix(2),iy(2),iz(1),iq(2),ir(2))
    fx1222 = p1*me%f(ix(1),iy(1),iz(2),iq(2),ir(2)) + q1*me%f(ix(2),iy(1),iz(2),iq(2),ir(2))
    fx2222 = p1*me%f(ix(1),iy(2),iz(2),iq(2),ir(2)) + q1*me%f(ix(2),iy(2),iz(2),iq(2),ir(2))

    fxy111 = p2*( fx1111 ) + q2*( fx2111 )
    fxy211 = p2*( fx1211 ) + q2*( fx2211 )
    fxy121 = p2*( fx1121 ) + q2*( fx2121 )
    fxy221 = p2*( fx1221 ) + q2*( fx2221 )
    fxy112 = p2*( fx1112 ) + q2*( fx2112 )
    fxy212 = p2*( fx1212 ) + q2*( fx2212 )
    fxy122 = p2*( fx1122 ) + q2*( fx2122 )
    fxy222 = p2*( fx1222 ) + q2*( fx2222 )

    fxyz11 = p3*( fxy111 ) + q3*( fxy211 )
    fxyz21 = p3*( fxy121 ) + q3*( fxy221 )
    fxyz12 = p3*( fxy112 ) + q3*( fxy212 )
    fxyz22 = p3*( fxy122 ) + q3*( fxy222 )

    fxyzq1 = p4*( fxyz11 ) + q4*( fxyz21 )
    fxyzq2 = p4*( fxyz12 ) + q4*( fxyz22 )

    fxyzqr = p5*fxyzq1 + q5*fxyzq2

    end subroutine interp_5d
!*****************************************************************************************

!*****************************************************************************************
!>
!  6D linear interpolation routine.

    pure subroutine interp_6d(me,x,y,z,q,r,s,fxyzqrs)

    implicit none

    class(linear_interp_6d),intent(inout) :: me
    double precision,intent(in)                   :: x
    double precision,intent(in)                   :: y
    double precision,intent(in)                   :: z
    double precision,intent(in)                   :: q
    double precision,intent(in)                   :: r
    double precision,intent(in)                   :: s
    double precision,intent(out)                  :: fxyzqrs  !! Interpolated \( f(x,y,z,q,r,s) \)

    integer,dimension(2) :: ix, iy, iz, iq, ir, is
    double precision :: p1, p2, p3, p4, p5, p6
    double precision :: q1, q2, q3, q4, q5, q6
    integer :: mflag
    double precision :: fx11111, fx21111, fx12111, fx22111, fx11211, fx21211, fx12211, &
                fx22211, fxy1111, fxy2111, fxy1211, fxy2211, fxyz111, fxyz211, &
                fxyzq11, fx11121, fx21121, fx12121, fx22121, fx11221, fx21221, &
                fx12221, fx22221, fxy1121, fxy2121, fxy1221, fxy2221, fxyz121, &
                fxyz221, fxyzq21, fx11112, fx21112, fx12112, fx22112, fx11212, &
                fx21212, fx12212, fx22212, fxy1112, fxy2112, fxy1212, fxy2212, &
                fxyz112, fxyz212, fxyzq12, fx11122, fx21122, fx12122, fx22122, &
                fx11222, fx21222, fx12222, fx22222, fxy1122, fxy2122, fxy1222, &
                fxy2222, fxyz122, fxyz222, fxyzq22, fxyzqr1, fxyzqr2

    call dintrv(me%x,x,me%ilox,ix(1),ix(2),mflag)
    call dintrv(me%y,y,me%iloy,iy(1),iy(2),mflag)
    call dintrv(me%z,z,me%iloz,iz(1),iz(2),mflag)
    call dintrv(me%q,q,me%iloq,iq(1),iq(2),mflag)
    call dintrv(me%r,r,me%ilor,ir(1),ir(2),mflag)
    call dintrv(me%s,s,me%ilos,is(1),is(2),mflag)

    q1 = (x-me%x(ix(1)))/(me%x(ix(2))-me%x(ix(1)))
    q2 = (y-me%y(iy(1)))/(me%y(iy(2))-me%y(iy(1)))
    q3 = (z-me%z(iz(1)))/(me%z(iz(2))-me%z(iz(1)))
    q4 = (q-me%q(iq(1)))/(me%q(iq(2))-me%q(iq(1)))
    q5 = (r-me%r(ir(1)))/(me%r(ir(2))-me%r(ir(1)))
    q6 = (s-me%s(is(1)))/(me%s(is(2))-me%s(is(1)))
    p1 = one-q1
    p2 = one-q2
    p3 = one-q3
    p4 = one-q4
    p5 = one-q5
    p6 = one-q6

    fx11111 = p1*me%f(ix(1),iy(1),iz(1),iq(1),ir(1),is(1)) + q1*me%f(ix(2),iy(1),iz(1),iq(1),ir(1),is(1))
    fx21111 = p1*me%f(ix(1),iy(2),iz(1),iq(1),ir(1),is(1)) + q1*me%f(ix(2),iy(2),iz(1),iq(1),ir(1),is(1))
    fx12111 = p1*me%f(ix(1),iy(1),iz(2),iq(1),ir(1),is(1)) + q1*me%f(ix(2),iy(1),iz(2),iq(1),ir(1),is(1))
    fx22111 = p1*me%f(ix(1),iy(2),iz(2),iq(1),ir(1),is(1)) + q1*me%f(ix(2),iy(2),iz(2),iq(1),ir(1),is(1))
    fx11211 = p1*me%f(ix(1),iy(1),iz(1),iq(2),ir(1),is(1)) + q1*me%f(ix(2),iy(1),iz(1),iq(2),ir(1),is(1))
    fx21211 = p1*me%f(ix(1),iy(2),iz(1),iq(2),ir(1),is(1)) + q1*me%f(ix(2),iy(2),iz(1),iq(2),ir(1),is(1))
    fx12211 = p1*me%f(ix(1),iy(1),iz(2),iq(2),ir(1),is(1)) + q1*me%f(ix(2),iy(1),iz(2),iq(2),ir(1),is(1))
    fx22211 = p1*me%f(ix(1),iy(2),iz(2),iq(2),ir(1),is(1)) + q1*me%f(ix(2),iy(2),iz(2),iq(2),ir(1),is(1))
    fx11121 = p1*me%f(ix(1),iy(1),iz(1),iq(1),ir(2),is(1)) + q1*me%f(ix(2),iy(1),iz(1),iq(1),ir(2),is(1))
    fx21121 = p1*me%f(ix(1),iy(2),iz(1),iq(1),ir(2),is(1)) + q1*me%f(ix(2),iy(2),iz(1),iq(1),ir(2),is(1))
    fx12121 = p1*me%f(ix(1),iy(1),iz(2),iq(1),ir(2),is(1)) + q1*me%f(ix(2),iy(1),iz(2),iq(1),ir(2),is(1))
    fx22121 = p1*me%f(ix(1),iy(2),iz(2),iq(1),ir(2),is(1)) + q1*me%f(ix(2),iy(2),iz(2),iq(1),ir(2),is(1))
    fx11221 = p1*me%f(ix(1),iy(1),iz(1),iq(2),ir(2),is(1)) + q1*me%f(ix(2),iy(1),iz(1),iq(2),ir(2),is(1))
    fx21221 = p1*me%f(ix(1),iy(2),iz(1),iq(2),ir(2),is(1)) + q1*me%f(ix(2),iy(2),iz(1),iq(2),ir(2),is(1))
    fx12221 = p1*me%f(ix(1),iy(1),iz(2),iq(2),ir(2),is(1)) + q1*me%f(ix(2),iy(1),iz(2),iq(2),ir(2),is(1))
    fx22221 = p1*me%f(ix(1),iy(2),iz(2),iq(2),ir(2),is(1)) + q1*me%f(ix(2),iy(2),iz(2),iq(2),ir(2),is(1))
    fx11112 = p1*me%f(ix(1),iy(1),iz(1),iq(1),ir(1),is(2)) + q1*me%f(ix(2),iy(1),iz(1),iq(1),ir(1),is(2))
    fx21112 = p1*me%f(ix(1),iy(2),iz(1),iq(1),ir(1),is(2)) + q1*me%f(ix(2),iy(2),iz(1),iq(1),ir(1),is(2))
    fx12112 = p1*me%f(ix(1),iy(1),iz(2),iq(1),ir(1),is(2)) + q1*me%f(ix(2),iy(1),iz(2),iq(1),ir(1),is(2))
    fx22112 = p1*me%f(ix(1),iy(2),iz(2),iq(1),ir(1),is(2)) + q1*me%f(ix(2),iy(2),iz(2),iq(1),ir(1),is(2))
    fx11212 = p1*me%f(ix(1),iy(1),iz(1),iq(2),ir(1),is(2)) + q1*me%f(ix(2),iy(1),iz(1),iq(2),ir(1),is(2))
    fx21212 = p1*me%f(ix(1),iy(2),iz(1),iq(2),ir(1),is(2)) + q1*me%f(ix(2),iy(2),iz(1),iq(2),ir(1),is(2))
    fx12212 = p1*me%f(ix(1),iy(1),iz(2),iq(2),ir(1),is(2)) + q1*me%f(ix(2),iy(1),iz(2),iq(2),ir(1),is(2))
    fx22212 = p1*me%f(ix(1),iy(2),iz(2),iq(2),ir(1),is(2)) + q1*me%f(ix(2),iy(2),iz(2),iq(2),ir(1),is(2))
    fx11122 = p1*me%f(ix(1),iy(1),iz(1),iq(1),ir(2),is(2)) + q1*me%f(ix(2),iy(1),iz(1),iq(1),ir(2),is(2))
    fx21122 = p1*me%f(ix(1),iy(2),iz(1),iq(1),ir(2),is(2)) + q1*me%f(ix(2),iy(2),iz(1),iq(1),ir(2),is(2))
    fx12122 = p1*me%f(ix(1),iy(1),iz(2),iq(1),ir(2),is(2)) + q1*me%f(ix(2),iy(1),iz(2),iq(1),ir(2),is(2))
    fx22122 = p1*me%f(ix(1),iy(2),iz(2),iq(1),ir(2),is(2)) + q1*me%f(ix(2),iy(2),iz(2),iq(1),ir(2),is(2))
    fx11222 = p1*me%f(ix(1),iy(1),iz(1),iq(2),ir(2),is(2)) + q1*me%f(ix(2),iy(1),iz(1),iq(2),ir(2),is(2))
    fx21222 = p1*me%f(ix(1),iy(2),iz(1),iq(2),ir(2),is(2)) + q1*me%f(ix(2),iy(2),iz(1),iq(2),ir(2),is(2))
    fx12222 = p1*me%f(ix(1),iy(1),iz(2),iq(2),ir(2),is(2)) + q1*me%f(ix(2),iy(1),iz(2),iq(2),ir(2),is(2))
    fx22222 = p1*me%f(ix(1),iy(2),iz(2),iq(2),ir(2),is(2)) + q1*me%f(ix(2),iy(2),iz(2),iq(2),ir(2),is(2))

    fxy1111 = p2*( fx11111 ) + q2*( fx21111 )
    fxy2111 = p2*( fx12111 ) + q2*( fx22111 )
    fxy1211 = p2*( fx11211 ) + q2*( fx21211 )
    fxy2211 = p2*( fx12211 ) + q2*( fx22211 )
    fxy1121 = p2*( fx11121 ) + q2*( fx21121 )
    fxy2121 = p2*( fx12121 ) + q2*( fx22121 )
    fxy1221 = p2*( fx11221 ) + q2*( fx21221 )
    fxy2221 = p2*( fx12221 ) + q2*( fx22221 )
    fxy1112 = p2*( fx11112 ) + q2*( fx21112 )
    fxy2112 = p2*( fx12112 ) + q2*( fx22112 )
    fxy1212 = p2*( fx11212 ) + q2*( fx21212 )
    fxy2212 = p2*( fx12212 ) + q2*( fx22212 )
    fxy1122 = p2*( fx11122 ) + q2*( fx21122 )
    fxy2122 = p2*( fx12122 ) + q2*( fx22122 )
    fxy1222 = p2*( fx11222 ) + q2*( fx21222 )
    fxy2222 = p2*( fx12222 ) + q2*( fx22222 )

    fxyz111 = p3*( fxy1111 ) + q3*( fxy2111 )
    fxyz211 = p3*( fxy1211 ) + q3*( fxy2211 )
    fxyz121 = p3*( fxy1121 ) + q3*( fxy2121 )
    fxyz221 = p3*( fxy1221 ) + q3*( fxy2221 )
    fxyz112 = p3*( fxy1112 ) + q3*( fxy2112 )
    fxyz212 = p3*( fxy1212 ) + q3*( fxy2212 )
    fxyz122 = p3*( fxy1122 ) + q3*( fxy2122 )
    fxyz222 = p3*( fxy1222 ) + q3*( fxy2222 )

    fxyzq11 = p4*( fxyz111 ) + q4*( fxyz211 )
    fxyzq21 = p4*( fxyz121 ) + q4*( fxyz221 )
    fxyzq12 = p4*( fxyz112 ) + q4*( fxyz212 )
    fxyzq22 = p4*( fxyz122 ) + q4*( fxyz222 )

    fxyzqr1 = p5*fxyzq11 + q5*fxyzq21
    fxyzqr2 = p5*fxyzq12 + q5*fxyzq22

    fxyzqrs = p6*fxyzqr1 + q6*fxyzqr2

    end subroutine interp_6d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns the indices in `xt` that bound `x`, to use for interpolation.
!  If outside the range, then the indices are returned that can
!  be used for extrapolation.
!  Precisely,
!
!```fortran
!         if            x < xt(1)   then ileft=1,   iright=2,    mflag=-1
!         if   xt(i) <= x < xt(i+1) then ileft=i,   iright=i+1,  mflag=0
!         if   xt(n) <= x           then ileft=n-1, iright=n,    mflag=1
!```
!
!### History
!
!  * interv written by carl de boor [5]
!  * dintrv author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * Jacob Williams, 2/24/2015 : updated to free-form Fortran.
!  * Jacob Williams, 2/17/2016 : additional refactoring (eliminated GOTOs).
!  * Jacob Williams, 2/22/2016 : modified bspline-fortran `dintrv` routine for
!    linear interpolation/extrapolation use.

    pure subroutine dintrv(xt,x,ilo,ileft,iright,mflag)

    implicit none

    double precision,dimension(:),intent(in) :: xt     !! a knot or break point vector
    double precision,intent(in)              :: x      !! argument
    integer,intent(inout)            :: ilo    !! an initialization parameter which must be set
                                               !! to 1 the first time the array `xt` is
                                               !! processed by dintrv. `ilo` contains information for
                                               !! efficient processing after the initial call and `ilo`
                                               !! must not be changed by the user.  each dimension
                                               !! requires a distinct `ilo` parameter.
    integer,intent(out)              :: ileft  !! left index
    integer,intent(out)              :: iright !! right index
    integer,intent(out)              :: mflag  !! signals when `x` lies out of bounds

    integer :: ihi, istep, imid, n

    n = size(xt)

    ihi = ilo + 1
    if ( ihi>=n ) then
        if ( x>=xt(n) ) then
            mflag = 1
            ileft = n-1
            iright= n
            return
        end if
        if ( n<=1 ) then
            mflag = -1
            ileft = 1
            iright= 2
            return
        end if
        ilo = n - 1
        ihi = n
    endif

    if ( x>=xt(ihi) ) then

        ! now x >= xt(ilo). find upper bound
        istep = 1
        do
            ilo = ihi
            ihi = ilo + istep
            if ( ihi>=n ) then
                if ( x>=xt(n) ) then
                    mflag = 1
                    ileft = n-1
                    iright= n
                    return
                end if
                ihi = n
            elseif ( x>=xt(ihi) ) then
                istep = istep*2
                cycle
            endif
            exit
        end do

    else

        if ( x>=xt(ilo) ) then
            mflag = 0
            ileft = ilo
            iright= ilo+1
            return
        end if
        ! now x <= xt(ihi). find lower bound
        istep = 1
        do
            ihi = ilo
            ilo = ihi - istep
            if ( ilo<=1 ) then
                ilo = 1
                if ( x<xt(1) ) then
                    mflag = -1
                    ileft = 1
                    iright= 2
                    return
                end if
            elseif ( x<xt(ilo) ) then
                istep = istep*2
                cycle
            endif
            exit
        end do

    endif

    ! now xt(ilo) <= x < xt(ihi). narrow the interval
    do
        imid = (ilo+ihi)/2
        if ( imid==ilo ) then
            mflag = 0
            ileft = ilo
            iright= ilo+1
            return
        end if
        ! note. it is assumed that imid = ilo in case ihi = ilo+1
        if ( x<xt(imid) ) then
            ihi = imid
        else
            ilo = imid
        endif
    end do

    end subroutine dintrv
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns true if all the elements in the array `x` are unique.
!  Note: the array must be sorted.
!
!@note This routine is not currently used in the module.

    pure function check_if_unique(x) result(unique)

    implicit none

    double precision,dimension(:),intent(in) :: x       !! a sorted array
    logical                          :: unique  !! true if all elements are unique

    integer :: i  !! counter

    unique = .true. ! initialize

    do i=1,size(x)-1
        if (x(i)==x(i+1)) then
            unique = .false.
            exit
        end if
    end do

    end function check_if_unique
!*****************************************************************************************

!*****************************************************************************************
!>
!  Sorts an array `dx` in increasing order,
!  carrying along an additional array `dy`.
!
!  Uses a non-recursive quicksort, reverting to insertion sort on arrays of
!  size <= 20. Dimension of `stack` limits array size to about \(2^32\).
!
!### License
!  * [Original LAPACK license](http://www.netlib.org/lapack/LICENSE.txt)
!
!### History
!  * Based on the LAPACK routine [DLASRT](http://www.netlib.org/lapack/explore-html/df/ddf/dlasrt_8f.html).
!  * Extensively modified by Jacob Williams, Feb. 2016. Converted to
!    modern Fortran and added the `dy` output. Removed the descending sort option.
!
!@note This routine is not currently used in the module.

    pure subroutine sort(dx,dy)

    implicit none

    double precision,dimension(:),intent(inout) :: dx  !! on entry, the array to be sorted.
                                               !! on exit, `dx` has been sorted into increasing order
                                               !! (`dx(1) <= ... <= dx(n)`) or into decreasing order
                                               !! (`dx(1) >= ... >= dx(n)`), depending on `id`.
    double precision,dimension(:),intent(inout) :: dy  !! array carried along with `dx`.

    integer,parameter :: select = 20  !! max size for using insertion sort.

    integer :: endd, i, j, n, start, stkpnt
    double precision :: d1, d2, d3, dmnmx, dmnmy, tmp
    integer,dimension(2,32) :: stack

    ! number of elements to sort:
    n = size(dx)

    if ( n>1 ) then

        stkpnt     = 1
        stack(1,1) = 1
        stack(2,1) = n

        do

            start  = stack(1,stkpnt)
            endd   = stack(2,stkpnt)
            stkpnt = stkpnt - 1
            if ( endd-start<=select .and. endd>start ) then

                ! do insertion sort on dx( start:endd )
                insertion: do i = start + 1, endd
                    do j = i, start + 1, -1
                        if ( dx(j)>=dx(j-1) ) cycle insertion
                        dmnmx   = dx(j)
                        dx(j)   = dx(j-1)
                        dx(j-1) = dmnmx
                        dmnmy   = dy(j)
                        dy(j)   = dy(j-1)
                        dy(j-1) = dmnmy
                    enddo
                enddo insertion

            elseif ( endd-start>select ) then

                ! partition dx( start:endd ) and stack parts, largest one first
                ! choose partition entry as median of 3

                d1 = dx(start)
                d2 = dx(endd)
                i  = (start+endd)/2
                d3 = dx(i)
                if ( d1<d2 ) then
                    if ( d3<d1 ) then
                        dmnmx = d1
                    elseif ( d3<d2 ) then
                        dmnmx = d3
                    else
                        dmnmx = d2
                    endif
                elseif ( d3<d2 ) then
                    dmnmx = d2
                elseif ( d3<d1 ) then
                    dmnmx = d3
                else
                    dmnmx = d1
                endif

                i = start - 1
                j = endd + 1
                do
                    do
                        j = j - 1
                        if ( dx(j)<=dmnmx ) exit
                    end do
                    do
                        i = i + 1
                        if ( dx(i)>=dmnmx ) exit
                    end do
                    if ( i<j ) then
                        tmp   = dx(i)
                        dx(i) = dx(j)
                        dx(j) = tmp
                        tmp   = dy(i)
                        dy(i) = dy(j)
                        dy(j) = tmp
                    else
                        exit
                    endif
                end do
                if ( j-start>endd-j-1 ) then
                    stkpnt          = stkpnt + 1
                    stack(1,stkpnt) = start
                    stack(2,stkpnt) = j
                    stkpnt          = stkpnt + 1
                    stack(1,stkpnt) = j + 1
                    stack(2,stkpnt) = endd
                else
                    stkpnt          = stkpnt + 1
                    stack(1,stkpnt) = j + 1
                    stack(2,stkpnt) = endd
                    stkpnt          = stkpnt + 1
                    stack(1,stkpnt) = start
                    stack(2,stkpnt) = j
                endif

            endif

            if ( stkpnt<=0 ) exit

        end do

    end if

    end subroutine sort
!*****************************************************************************************

!*****************************************************************************************
!>
!  Check the validity of the inputs to the initialize routines.
!  Prints warning message if there is an error,
!  and also sets `ierr` (/=0 if there were any errors).
!
!  Supports up to 6D: x,y,z,q,r,s
!
!# History
!  * Jacob Williams, 2/24/2015 : Created this routine.
!  * Jacob Williams, 2/23/2016 : modified for linear interp module.

    pure subroutine check_inputs(x,y,z,q,r,s,ierr)

    use iso_fortran_env,    only: error_unit

    implicit none

    double precision,dimension(:),intent(in),optional  :: x     !! `x` abscissa vector
    double precision,dimension(:),intent(in),optional  :: y     !! `y` abscissa vector
    double precision,dimension(:),intent(in),optional  :: z     !! `z` abscissa vector
    double precision,dimension(:),intent(in),optional  :: q     !! `q` abscissa vector
    double precision,dimension(:),intent(in),optional  :: r     !! `r` abscissa vector
    double precision,dimension(:),intent(in),optional  :: s     !! `s` abscissa vector
    integer,intent(out)                        :: ierr  !! `0` : no problems,
                                                        !! `1` : `x` is not strictly increasing,
                                                        !! `2` : `y` is not strictly increasing,
                                                        !! `3` : `z` is not strictly increasing,
                                                        !! `4` : `q` is not strictly increasing,
                                                        !! `5` : `r` is not strictly increasing,
                                                        !! `6` : `s` is not strictly increasing,

    ierr = 0  ! initialize

    if (present(x)) call check(x,1,ierr); if (ierr/=0) return
    if (present(y)) call check(y,2,ierr); if (ierr/=0) return
    if (present(z)) call check(z,3,ierr); if (ierr/=0) return
    if (present(q)) call check(q,4,ierr); if (ierr/=0) return
    if (present(r)) call check(r,5,ierr); if (ierr/=0) return
    if (present(s)) call check(s,6,ierr); if (ierr/=0) return

    contains
!*****************************************************************************************

        pure subroutine check(v,error_code,ierr)

        implicit none

        double precision,dimension(:),intent(in) :: v          !! abcissae vector
        integer,intent(in)               :: error_code !! error code for check
        integer,intent(inout)            :: ierr       !! will be set to `error_code` if there is a problem

        integer :: i  !! counter
        integer :: n  !! size of the input `v` array

        n = size(v)
        do i=2,n
            if (v(i) <= v(i-1)) then
                ierr = error_code
                exit
            end if
        end do

        end subroutine check

    end subroutine check_inputs

!*****************************************************************************************
! Linear Algebra
!*****************************************************************************************

	! Calculate Rank of Matrix A (N x M)
	function get_rank(A,N,M,info) result(r)
	  integer::N,M,NM,info,r,i,j
	  double precision,dimension(N,M)::A
	  double precision,dimension(:,:),allocatable::B
	  NM=max(M,N)
	  ! copy A in B with allocation of B, because get_rank destroys B
	  B=A
	  call calc_rank(B,N,M,NM,info,r)
	end function


subroutine drotg(da, db, dc, ds)

!     designed by c.l.lawson, jpl, 1977 sept 08

!     construct the givens transformation

!         ( dc  ds )
!     g = (        ) ,    dc**2 + ds**2 = 1 ,
!         (-ds  dc )

!     which zeros the second entry of the 2-vector  (da,db)**t .

!     the quantity r = (+/-)sqrt(da**2 + db**2) overwrites da in
!     storage.  the value of db is overwritten by a value z which
!     allows dc and ds to be recovered by the following algorithm:
!           if z=1  set  dc=0.d0  and  ds=1.d0
!           if dabs(z) < 1  set  dc=sqrt(1-z**2)  and  ds=z
!           if dabs(z) > 1  set  dc=1/z  and  ds=sqrt(1-dc**2)

!     normally, the subprogram drot(n,dx,incx,dy,incy,dc,ds) will
!     next be called to apply the transformation to a 2 by n matrix.

! ------------------------------------------------------------------

double precision, intent(in out)  :: da
double precision, intent(in out)  :: db
double precision, intent(out)     :: dc
double precision, intent(out)     :: ds

double precision  :: u, v, r
if (abs(da) <= abs(db)) go to 10

! *** here abs(da) > abs(db) ***

u = da + da
v = db / u

!     note that u and r have the sign of da

r = sqrt(.25d0 + v**2) * u

!     note that dc is positive

dc = da / r
ds = v * (dc + dc)
db = ds
da = r
return

! *** here abs(da) <= abs(db) ***

10 if (db == 0.d0) go to 20
u = db + db
v = da / u

!     note that u and r have the sign of db
!     (r is immediately stored in da)

da = sqrt(.25d0 + v**2) * u

!     note that ds is positive

ds = db / da
dc = v * (ds + ds)
if (dc == 0.d0) go to 15
db = 1.d0 / dc
return
15 db = 1.d0
return

! *** here da = db = 0.d0 ***

20 dc = 1.d0
ds = 0.d0
return

end subroutine drotg
!
!
subroutine dswap1 (n, dx, dy)

!     interchanges two vectors.
!     uses unrolled loops for increments equal one.
!     jack dongarra, linpack, 3/11/78.
!     this version is for increments = 1.

integer, intent(in)        :: n
double precision, intent(in out)  :: dx(*)
double precision, intent(in out)  :: dy(*)

double precision  :: dtemp
integer    :: i, m, mp1

if(n <= 0) return

!       code for both increments equal to 1

!       clean-up loop

m = mod(n,3)
if( m == 0 ) go to 40
do  i = 1,m
  dtemp = dx(i)
  dx(i) = dy(i)
  dy(i) = dtemp
end do
if( n < 3 ) return
40 mp1 = m + 1
do  i = mp1,n,3
  dtemp = dx(i)
  dx(i) = dy(i)
  dy(i) = dtemp
  dtemp = dx(i + 1)
  dx(i + 1) = dy(i + 1)
  dy(i + 1) = dtemp
  dtemp = dx(i + 2)
  dx(i + 2) = dy(i + 2)
  dy(i + 2) = dtemp
end do
return
end subroutine  dswap1


subroutine  drot1 (n, dx, dy, c, s)

!     applies a plane rotation.
!     jack dongarra, linpack, 3/11/78.
!     this version is for increments = 1.

integer, intent(in)        :: n
double precision, intent(in out)  :: dx(*)
double precision, intent(in out)  :: dy(*)
double precision, intent(in)      :: c
double precision, intent(in)      :: s

double precision  :: dtemp
integer    :: i

if(n <= 0) return
!       code for both increments equal to 1

do  i = 1,n
  dtemp = c*dx(i) + s*dy(i)
  dy(i) = c*dy(i) - s*dx(i)
  dx(i) = dtemp
end do
return
end subroutine  drot1


subroutine calc_rank(a,n,p,m,info,rank)
	double precision, intent(in out)  :: a(:,:)
	integer, intent(in)        :: n  ! number of rows
	integer, intent(in)        :: p  ! number of cols
	integer, intent(in)        :: m  ! max(n,p)
	integer,intent(out)		   :: rank, info
	double precision:: s(m), u(n,n), v(p,p), e(p), eps
	integer :: i, imax

	call dsvdc(a, n, p, s, e, u, v, 0, info)
  	eps = 5e-10
  	rank = 0
  	imax=min(n,p)
  	do i=1, imax
    	if (s(i) .gt. eps) then
      		rank = rank + 1
    	end if
  	end do

end subroutine

subroutine dsvdc(x, n, p, s, e, u, v, job, info)

integer, intent(in)        :: n
integer, intent(in)        :: p
double precision, intent(in out)  :: x(:,:)
double precision, intent(out)     :: s(:)
double precision, intent(out)     :: e(:)
double precision, intent(out)     :: u(:,:)
double precision, intent(out)     :: v(:,:)
integer, intent(in)        :: job
integer, intent(out)       :: info

!     dsvdc is a subroutine to reduce a double precision nxp matrix x
!     by orthogonal transformations u and v to diagonal form.  the
!     diagonal elements s(i) are the singular values of x.  the
!     columns of u are the corresponding left singular vectors,
!     and the columns of v the right singular vectors.

!     on entry

!         x         double precision(ldx,p), where ldx.ge.n.
!                   x contains the matrix whose singular value
!                   decomposition is to be computed.  x is
!                   destroyed by dsvdc.

!         ldx       integer.
!                   ldx is the leading dimension of the array x.

!         n         integer.
!                   n is the number of rows of the matrix x.

!         p         integer.
!                   p is the number of columns of the matrix x.

!         ldu       integer.
!                   ldu is the leading dimension of the array u.
!                   (see below).

!         ldv       integer.
!                   ldv is the leading dimension of the array v.
!                   (see below).

!         job       integer.
!                   job controls the computation of the singular
!                   vectors.  it has the decimal expansion ab
!                   with the following meaning

!                        a.eq.0    do not compute the left singular vectors.
!                        a.eq.1    return the n left singular vectors in u.
!                        a.ge.2    return the first min(n,p) singular
!                                  vectors in u.
!                        b.eq.0    do not compute the right singular vectors.
!                        b.eq.1    return the right singular vectors in v.

!     on return

!         s         double precision(mm), where mm=min(n+1,p).
!                   the first min(n,p) entries of s contain the singular
!                   values of x arranged in descending order of magnitude.

!         e         double precision(p).
!                   e ordinarily contains zeros.  however see the
!                   discussion of info for exceptions.

!         u         double precision(ldu,k), where ldu.ge.n.  if
!                                   joba.eq.1 then k.eq.n, if joba.ge.2
!                                   then k.eq.min(n,p).
!                   u contains the matrix of left singular vectors.
!                   u is not referenced if joba.eq.0.  if n.le.p
!                   or if joba.eq.2, then u may be identified with x
!                   in the subroutine call.

!         v         double precision(ldv,p), where ldv.ge.p.
!                   v contains the matrix of right singular vectors.
!                   v is not referenced if job.eq.0.  if p.le.n,
!                   then v may be identified with x in the
!                   subroutine call.

!         info      integer.
!                   the singular values (and their corresponding singular
!                   vectors) s(info+1),s(info+2),...,s(m) are correct
!                   (here m=min(n,p)).  thus if info.eq.0, all the
!                   singular values and their vectors are correct.
!                   in any event, the matrix b = trans(u)*x*v is the
!                   bidiagonal matrix with the elements of s on its diagonal
!                   and the elements of e on its super-diagonal (trans(u)
!                   is the transpose of u).  thus the singular values
!                   of x and b are the same.

!     linpack. this version dated 03/19/79 .
!     g.w. stewart, university of maryland, argonne national lab.

!     dsvdc uses the following functions and subprograms.

!     external drot
!     blas daxpy,ddot,dscal,dswap,dnrm2,drotg
!     fortran dabs,dmax1,max0,min0,mod,dsqrt

!     internal variables

integer :: iter, j, jobu, k, kase, kk, l, ll, lls, lm1, lp1, ls, lu, m, maxit,  &
    mm, mm1, mp1, nct, nctp1, ncu, nrt, nrtp1
double precision :: t, work(n)
double precision :: b, c, cs, el, emm1, f, g, scale, shift, sl, sm, sn,  &
    smm1, t1, test, ztest
logical :: wantu, wantv

!     set the maximum number of iterations.

maxit = 30

!     determine what is to be computed.

wantu = .false.
wantv = .false.
jobu = mod(job,100)/10
ncu = n
if (jobu > 1) ncu = min(n,p)
if (jobu /= 0) wantu = .true.
if (mod(job,10) /= 0) wantv = .true.

!     reduce x to bidiagonal form, storing the diagonal elements
!     in s and the super-diagonal elements in e.

info = 0
nct = min(n-1, p)
s(1:nct+1) = 0.0
nrt = max(0, min(p-2,n))
lu = max(nct,nrt)
if (lu < 1) go to 170
do  l = 1, lu
  lp1 = l + 1
  if (l > nct) go to 20

!           compute the transformation for the l-th column and
!           place the l-th diagonal in s(l).

  s(l) = sqrt( sum( x(l:n,l)**2 ) )
  if (s(l) == 0.0d0) go to 10
  if (x(l,l) /= 0.0d0) s(l) = sign(s(l), x(l,l))
  x(l:n,l) = x(l:n,l) / s(l)
  x(l,l) = 1.0d0 + x(l,l)

  10 s(l) = -s(l)

  20 if (p < lp1) go to 50
  do  j = lp1, p
    if (l > nct) go to 30
    if (s(l) == 0.0d0) go to 30

!              apply the transformation.

    t = -dot_product(x(l:n,l), x(l:n,j)) / x(l,l)
    x(l:n,j) = x(l:n,j) + t * x(l:n,l)

!           place the l-th row of x into  e for the
!           subsequent calculation of the row transformation.

    30 e(j) = x(l,j)
  end do

  50 if (.not.wantu .or. l > nct) go to 70

!           place the transformation in u for subsequent back multiplication.

  u(l:n,l) = x(l:n,l)

  70 if (l > nrt) cycle

!           compute the l-th row transformation and place the
!           l-th super-diagonal in e(l).

  e(l) = sqrt( sum( e(lp1:p)**2 ) )
  if (e(l) == 0.0d0) go to 80
  if (e(lp1) /= 0.0d0) e(l) = sign(e(l), e(lp1))
  e(lp1:lp1+p-l-1) = e(lp1:p) / e(l)
  e(lp1) = 1.0d0 + e(lp1)

  80 e(l) = -e(l)
  if (lp1 > n .or. e(l) == 0.0d0) go to 120

!              apply the transformation.

  work(lp1:n) = 0.0d0
  do  j = lp1, p
    work(lp1:lp1+n-l-1) = work(lp1:lp1+n-l-1) + e(j) * x(lp1:lp1+n-l-1,j)
  end do
  do  j = lp1, p
    x(lp1:lp1+n-l-1,j) = x(lp1:lp1+n-l-1,j) - (e(j)/e(lp1)) * work(lp1:lp1+n-l-1)
  end do

  120 if (.not.wantv) cycle

!              place the transformation in v for subsequent
!              back multiplication.

  v(lp1:p,l) = e(lp1:p)
end do

!     set up the final bidiagonal matrix of order m.

170 m = min(p,n+1)
nctp1 = nct + 1
nrtp1 = nrt + 1
if (nct < p) s(nctp1) = x(nctp1,nctp1)
if (n < m) s(m) = 0.0d0
if (nrtp1 < m) e(nrtp1) = x(nrtp1,m)
e(m) = 0.0d0

!     if required, generate u.

if (.not.wantu) go to 300
if (ncu < nctp1) go to 200
do  j = nctp1, ncu
  u(1:n,j) = 0.0
  u(j,j) = 1.0
end do

200 do  ll = 1, nct
  l = nct - ll + 1
  if (s(l) == 0.0d0) go to 250
  lp1 = l + 1
  if (ncu < lp1) go to 220
  do  j = lp1, ncu
    t = -dot_product(u(l:n,l), u(l:n,j)) / u(l,l)
    u(l:n,j) = u(l:n,j) + t * u(l:n,l)
  end do

  220 u(l:n,l) = -u(l:n,l)
  u(l,l) = 1.0d0 + u(l,l)
  lm1 = l - 1
  if (lm1 < 1) cycle
  u(1:lm1,l) = 0.0
  cycle

  250 u(1:n,l) = 0.0
  u(l,l) = 1.0
end do

!     if it is required, generate v.

300 if (.not.wantv) go to 350
do  ll = 1, p
  l = p - ll + 1
  lp1 = l + 1
  if (l > nrt) go to 320
  if (e(l) == 0.0d0) go to 320
  do  j = lp1, p
    t = -dot_product(v(lp1:lp1+p-l-1,l), v(lp1:lp1+p-l-1,j)) / v(lp1,l)
    v(lp1:lp1+p-l-1,j) = v(lp1:lp1+p-l-1,j) + t * v(lp1:lp1+p-l-1,l)
  end do

  320 v(1:p,l) = 0.0d0
  v(l,l) = 1.0d0
end do

!     main iteration loop for the singular values.

350 mm = m
iter = 0

!        quit if all the singular values have been found.

!     ...exit
360 if (m == 0) go to 620

!        if too many iterations have been performed, set flag and return.

if (iter < maxit) go to 370
info = m
!     ......exit
go to 620

!        this section of the program inspects for negligible elements
!        in the s and e arrays.  on completion
!        the variables kase and l are set as follows.

!           kase = 1     if s(m) and e(l-1) are negligible and l < m
!           kase = 2     if s(l) is negligible and l < m
!           kase = 3     if e(l-1) is negligible, l < m, and
!                        s(l), ..., s(m) are not negligible (qr step).
!           kase = 4     if e(m-1) is negligible (convergence).

370 do  ll = 1, m
  l = m - ll
!        ...exit
  if (l == 0) exit
  test = abs(s(l)) + abs(s(l+1))
  ztest = test + abs(e(l))
  if (ztest /= test) cycle
  e(l) = 0.0d0
!        ......exit
  exit
end do

if (l /= m - 1) go to 410
kase = 4
go to 480

410 lp1 = l + 1
mp1 = m + 1
do  lls = lp1, mp1
  ls = m - lls + lp1
!           ...exit
  if (ls == l) exit
  test = 0.0d0
  if (ls /= m) test = test + abs(e(ls))
  if (ls /= l + 1) test = test + abs(e(ls-1))
  ztest = test + abs(s(ls))
  if (ztest /= test) cycle
  s(ls) = 0.0d0
!           ......exit
  exit
end do

if (ls /= l) go to 450
kase = 3
go to 480

450 if (ls /= m) go to 460
kase = 1
go to 480

460 kase = 2
l = ls
480 l = l + 1

!        perform the task indicated by kase.

select case ( kase )
  case (    1)
    go to 490
  case (    2)
    go to 520
  case (    3)
    go to 540
  case (    4)
    go to 570
end select

!        deflate negligible s(m).

490 mm1 = m - 1
f = e(m-1)
e(m-1) = 0.0d0
do  kk = l, mm1
  k = mm1 - kk + l
  t1 = s(k)
  call drotg(t1, f, cs, sn)
  s(k) = t1
  if (k == l) go to 500
  f = -sn*e(k-1)
  e(k-1) = cs*e(k-1)

  500 if (wantv) call drot1(p, v(1:,k), v(1:,m), cs, sn)
end do
go to 610

!        split at negligible s(l).

520 f = e(l-1)
e(l-1) = 0.0d0
do  k = l, m
  t1 = s(k)
  call drotg(t1, f, cs, sn)
  s(k) = t1
  f = -sn*e(k)
  e(k) = cs*e(k)
  if (wantu) call drot1(n, u(1:,k), u(1:,l-1), cs, sn)
end do
go to 610

!        perform one qr step.

!           calculate the shift.

540 scale = max(abs(s(m)), abs(s(m-1)), abs(e(m-1)), abs(s(l)), abs(e(l)))
sm = s(m)/scale
smm1 = s(m-1)/scale
emm1 = e(m-1)/scale
sl = s(l)/scale
el = e(l)/scale
b = ((smm1 + sm)*(smm1 - sm) + emm1**2)/2.0d0
c = (sm*emm1)**2
shift = 0.0d0
if (b == 0.0d0 .and. c == 0.0d0) go to 550
shift = sqrt(b**2+c)
if (b < 0.0d0) shift = -shift
shift = c/(b + shift)

550 f = (sl + sm)*(sl - sm) - shift
g = sl*el

!           chase zeros.

mm1 = m - 1
do  k = l, mm1
  call drotg(f, g, cs, sn)
  if (k /= l) e(k-1) = f
  f = cs*s(k) + sn*e(k)
  e(k) = cs*e(k) - sn*s(k)
  g = sn*s(k+1)
  s(k+1) = cs*s(k+1)
  if (wantv) call drot1(p, v(1:,k), v(1:,k+1), cs, sn)
  call drotg(f, g, cs, sn)
  s(k) = f
  f = cs*e(k) + sn*s(k+1)
  s(k+1) = -sn*e(k) + cs*s(k+1)
  g = sn*e(k+1)
  e(k+1) = cs*e(k+1)
  if (wantu .and. k < n) call drot1(n, u(1:,k), u(1:,k+1), cs, sn)
end do
e(m-1) = f
iter = iter + 1
go to 610

!        convergence.

!           make the singular value  positive.

570 if (s(l) >= 0.0d0) go to 590
s(l) = -s(l)
if (wantv) v(1:p,l) = -v(1:p,l)

!           order the singular value.

590 if (l == mm) go to 600
!           ...exit
if (s(l) >= s(l+1)) go to 600
t = s(l)
s(l) = s(l+1)
s(l+1) = t
if (wantv .and. l < p) call dswap1(p, v(1:,l), v(1:,l+1))
if (wantu .and. l < n) call dswap1(n, u(1:,l), u(1:,l+1))
l = l + 1
go to 590

600 iter = 0
m = m - 1

610 go to 360

620 return
end subroutine dsvdc


subroutine inv(matrix, inverse, n, errorflag)
	implicit none
	!declarations
	integer, intent(in) :: n
	integer, intent(out) :: errorflag  !return error status. -1 for error, 0 for normal
	double precision, intent(in), dimension(n,n) :: matrix  !input matrix
	double precision, intent(out), dimension(n,n) :: inverse !inverted matrix

	logical :: flag = .true.
	integer :: i, j, k, l
	double precision :: m
	double precision, dimension(n,2*n) :: augmatrix !augmented matrix

	!augment input matrix with an identity matrix
	do i = 1, n
		do j = 1, 2*n
			if (j <= n ) then
				augmatrix(i,j) = matrix(i,j)
			else if ((i+n) == j) then
				augmatrix(i,j) = 1
			else
				augmatrix(i,j) = 0
			endif
		end do
	end do

	!reduce augmented matrix to upper traingular form
	do k =1, n-1
		if (augmatrix(k,k) == 0) then
			flag = .false.
			do i = k+1, n
				if (augmatrix(i,k) /= 0) then
					do j = 1,2*n
						augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
					end do
					flag = .true.
					exit
				endif
				if (flag .eqv. .false.) then
					print*, "matrix is non - invertible"
					inverse = 0
					errorflag = -1
					return
				endif
			end do
		endif
		do j = k+1, n
			m = augmatrix(j,k)/augmatrix(k,k)
			do i = k, 2*n
				augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
			end do
		end do
	end do

	!test for invertibility
	do i = 1, n
		if (augmatrix(i,i) == 0) then
			print*, "matrix is non - invertible"
			inverse = 0
			errorflag = -1
			return
		endif
	end do

	!make diagonal elements as 1
	do i = 1 , n
		m = augmatrix(i,i)
		do j = i , (2 * n)
			   augmatrix(i,j) = (augmatrix(i,j) / m)
		end do
	end do

	!reduced right side half of augmented matrix to identity matrix
	do k = n-1, 1, -1
		do i =1, k
		m = augmatrix(i,k+1)
			do j = k, (2*n)
				augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
			end do
		end do
	end do

	!store answer
	do i =1, n
		do j = 1, n
			inverse(i,j) = augmatrix(i,j+n)
		end do
	end do
	errorflag = 0
end subroutine inv

!  ***************************************************************
!  * given an n x n matrix a, this routine replaces it by the lu *
!  * decomposition of a rowwise permutation of itself. a and n   *
!  * are input. indx is an output vector which records the row   *
!  * permutation effected by the partial pivoting; d is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. this routine is used *
!  * in combination with lubksb to solve linear equations or to  *
!  * invert a matrix. return code is 1, if matrix is singular.   *
!  ***************************************************************
 subroutine ludcmp(a,n,indx,d,code)
 integer n,i,j,k,imax
 integer, parameter :: nmax=100
 double precision, parameter :: tiny=1.5d-16
 double precision  amax,dum, sum, a(n,n),vv(nmax)
 integer code, d, indx(n)

 d=1; code=0

 do i=1,n
   amax=0.d0
   do j=1,n
     if (abs(a(i,j)).gt.amax) amax=abs(a(i,j))
   end do ! j loop
   if(amax.lt.tiny) then
     code = 1
     return
   end if
   vv(i) = 1.d0 / amax
 end do ! i loop

 do j=1,n
   do i=1,j-1
     sum = a(i,j)
     do k=1,i-1
       sum = sum - a(i,k)*a(k,j)
     end do ! k loop
     a(i,j) = sum
   end do ! i loop
   amax = 0.d0
   do i=j,n
     sum = a(i,j)
     do k=1,j-1
       sum = sum - a(i,k)*a(k,j)
     end do ! k loop
     a(i,j) = sum
     dum = vv(i)*abs(sum)
     if(dum.ge.amax) then
       imax = i
       amax = dum
     end if
   end do ! i loop

   if(j.ne.imax) then
     do k=1,n
       dum = a(imax,k)
       a(imax,k) = a(j,k)
       a(j,k) = dum
     end do ! k loop
     d = -d
     vv(imax) = vv(j)
   end if

   indx(j) = imax
   if(abs(a(j,j)) < tiny) a(j,j) = tiny

   if(j.ne.n) then
     dum = 1.d0 / a(j,j)
     do i=j+1,n
       a(i,j) = a(i,j)*dum
     end do ! i loop
   end if
 end do ! j loop

 return
 end subroutine ludcmp


!  ******************************************************************
!  * solves the set of n linear equations a . x = b.  here a is     *
!  * input, not as the matrix a but rather as its lu decomposition, *
!  * determined by the routine ludcmp. indx is input as the permuta-*
!  * tion vector returned by ludcmp. b is input as the right-hand   *
!  * side vector b, and returns with the solution vector x. a, n and*
!  * indx are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. this routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
 subroutine lubksb(a,n,indx,b)
 integer n,ii,i,ll,j
 double precision  sum, a(n,n),b(n)
 integer indx(n)

 ii = 0

 do i=1,n
   ll = indx(i)
   sum = b(ll)
   b(ll) = b(i)
   if(ii.ne.0) then
     do j=ii,i-1
       sum = sum - a(i,j)*b(j)
     end do ! j loop
   else if(sum.ne.0.d0) then
     ii = i
   end if
   b(i) = sum
 end do ! i loop

 do i=n,1,-1
   sum = b(i)
   if(i < n) then
     do j=i+1,n
       sum = sum - a(i,j)*b(j)
     end do ! j loop
   end if
   b(i) = sum / a(i,i)
 end do ! i loop

 return
 end subroutine lubksb


 subroutine inv2(a, y, n, rc)
	implicit none
	!declarations
	integer, intent(in) :: n
	integer, intent(out) :: rc  !return error status. -1 for error, 0 for normal
	double precision, intent(in), dimension(n,n) :: a  !input matrix
	double precision, intent(out), dimension(n,n) :: y !inverted matrix
	integer, dimension(n) :: indx
	logical :: flag = .true.
	integer :: i,j,d

	y=0
	do i=1,n
		y(i,i)=1.0
	end do

	!call lu decomposition routine (only once)
    call ludcmp(a,n,indx,d,rc)

    !call solver if previous return code is ok
    !to obtain inverse of a one column at a time
    if (rc.eq.0) then
      do j=1, n
        call lubksb(a,n,indx,y(1,j))
      end do
    else
    	rc=-1
    endif

    !the inverse matrix is now in matrix y
    !the original matrix a is destroyed

 end subroutine

    ! linear system A*X=B with symmetric matrix
  	subroutine linsolvesym(A,B,X,N,info)
  		integer::N
        double precision,dimension(N,1)::B
        double precision,dimension(N,N)::A
        double precision,dimension(N)::WORK
        double precision,dimension(N)::X
        integer,dimension(N)::IPIV
        integer::info
  		call DSYSV('U',N,1,A,N,IPIV,B,N,WORK,N,info)
  		X=B(:,1)
    end subroutine
!
    ! linear system A*X=B with generic matrix
  	subroutine linsolve(A,B,X,N,info)
  		integer::N
        double precision,dimension(N)::B
        double precision,dimension(N,N)::A
        double precision,dimension(N)::X
        double precision,dimension(N,1)::B1
        integer,dimension(N)::IPIV
        integer::info
        B1(:,1)=B
  		call DGESV(N,1,A,N,IPIV,B1,N,info)
  		X=B1(:,1)
    end subroutine
 
    ! linear system A*X=B with generic matrix
	! Aaa Aaw Xa  Ba
	! Awa Aww Xw  Bw
	! where Aww is constant and invertible
	! Aaa(1:1) /= 0
  	subroutine linsolve_partitions(A,B,X,invAww,Xa0,N,info)
  		integer::N
        double precision,dimension(N)::B
        double precision,dimension(N,N)::A
        double precision,dimension(N)::X
        double precision,dimension(N-1,N-1),intent(inout)::invAww
        double precision,dimension(N,1)::B1,X1
		double precision::Xa0,Xa1
        integer::info,iter
		B1(:,1)=B
		! test A(1,1)
		if (abs(A(1,1))<1d-16) then
			info=-2
			return
		end if	
		! perform inversion of Aww only once
		if (all(invAww == 0)) then
			call inv2(A(2:N,2:N), invAww, N-1, info)
			if (info.ne.0) return
		end if
		!    Xw=Aww^-1(Bw-Awa*Xa)
		!    Xa=(Ba-Aaw*Xw)/Aaa
		Xa1=Xa0
		do iter=1,500
			X1(2:N,1)=matmul(invAww,B1(2:N,1)-Xa1*A(2:N,1))
			X1(1,1)=(B1(1,1)-dot_product(A(1,2:N),X1(2:N,1)))/A(1,1)
			if (abs(Xa1-X1(1,1))<1d-3) then	
				exit
			else
				Xa1=0.1*Xa1+0.9*X1(1,1)
			end if	
		end do
		if (iter.ge.500) then
			info=-3
			return
		end if	
  		X=X1(:,1)
    end subroutine
    
	! Adaptive Simpson Integration function
	subroutine Integration_adaptive_simpson(f, a, b, tol, result)
		implicit none
        
        interface
          function f(x) result(y)
            double precision :: y
            double precision, intent(in) :: x
          end function f
        end interface
        
        double precision, intent(in) :: a, b, tol
        double precision, intent(out) :: result
        integer :: max_recursions
        logical :: converged
        
            max_recursions = 2000
            converged = .true.        
        
            result = simpson_recursive_function(f, a, b, tol, max_recursions, converged)

            if (.not. converged) then
              write(*,*) "Numerical Integration does not converge"
            end if
	contains

        recursive function simpson_recursive_function(f,a, b, tol, max_recursions, converged) result(integral)
          double precision :: a, b, tol, c, integral
          double precision :: fa, fb, fc, simpson1, simpson2, error
          integer, intent(inout) :: max_recursions
          logical, intent(out) :: converged
          
          interface
            function f(x) result(y)
              double precision :: y
              double precision, intent(in) :: x
            end function f
          end interface
          
          if (max_recursions <= 0) then
             integral = 0d0
             converged = .false.
             write(*,*) 'Error: ',error
             return
          end if

          c = (a + b)/2d0

          fa = f(a)
          fb = f(b)
          fc = f(c)

          simpson1 = (b - a) / 6d0 * (fa + 4d0 * fc + fb)
          simpson2 = (b - a) / 12d0 * (fa + 4d0 * f((a + c) / 2d0) + 2d0 * fc + 4d0 * f((c + b) / 2d0) + fb)

          error = abs(simpson2 - simpson1)

          if (error < 15d0 * tol) then
            integral = simpson2 + (simpson2 - simpson1) / 15d0
          else
            max_recursions = max_recursions - 1
            ! Ricorsione su ciascun intervallo
            !integral = simpson_recursive_function(f,gamma,beta, a, c, tol / 2d0,max_recursions, converged) + simpson_recursive_function(f,gamma,beta, c, b, tol / 2d0,max_recursions, converged)
            integral = simpson_recursive_function(f, a, c, tol,max_recursions, converged) + simpson_recursive_function(f, c, b, tol,max_recursions, converged)
          end if
        end function simpson_recursive_function
        
	end subroutine Integration_adaptive_simpson

end module
