! ----------------------------------------------------------------------------------------------	
module property_ftns

   ! Variables Declaration----------------------------------------------------------------------	
	
   implicit none
        
   ! For loop
   integer :: i, j, k, l, nfile
   
   ! File I/O system
   character(80)                            :: infile, outfile, invelfile, options, sym_file
   real                                     :: dummyr
   integer, dimension(20)                   :: dummyi(20)
   character(4)                             :: dummyc
   character(80), allocatable, dimension(:) :: filename
   integer                                  :: num_of_args

   ! Checking Trajectory variables
   integer 	 :: init_check, dt_check, end_check, num_check
   double precision :: max_r

   ! Data variables
   real, allocatable, dimension(:,:,:)    :: x, y, z
   real(8), allocatable, dimension(:,:,:) :: box 	   ! Is 3-D array best ?
   integer                                :: nmedia, total_traj
   integer, allocatable, dimension(:)     :: ntraj
   				         !nmedia : number of particles
   				         !total_traj : number of trajectory

! ----------------------------------------------------------------------------------------------	

   contains

   subroutine dcd_reader() ! for ensembles

   !------------------------------------------------!
   !                                                !
   !                Main Link Program               !
   !              made by Park HyungShick           !
   !                modified by EunHyeok            !
   !                modified by Jina                !
   !                                                !
   !               Ver.1 2020.06.26                 !
   !               Ver.2 2022.01.03                 !
   !   		   Ver.3 2025.04.27                 !
   !                                                !
   !------------------------------------------------!

   ! Refer to JENOGMIN KIM READ_DCD program
   !          Mar. 18, 2013

   
   ! Argument-----------------------------------------------------------------	

        num_of_args=iargc() ! # of files

        if (num_of_args .eq. 0) then
           print*, "------------ERROR------------"
           print*, "USAGE: ./run.x trajectory ..."
           print*, "-----------------------------"
           stop
        endif

	allocate(filename(num_of_args), ntraj(num_of_args))

	print*
        print*, "-----------------------------"
	print*, "(  Let's start!  )> \(>.< ;;)"
        print*, "-----------------------------"
	print*
	print*, num_of_args, " files are loaded"
	print*
	do i=1,num_of_args
           call getarg(i, filename(i)) ! ith arg is saved at filename(i)
	enddo
	
   ! Read DCD file------------------------------------------------------------------------------	
    
	! Checking ensembles to have same info--------------------------------------------------
	do nfile=1,num_of_args
	
           open(10+nfile,file=trim(filename(nfile)),status='old',form='unformatted')

           read(10+nfile) dummyc, ntraj(nfile), (dummyi(i),i=1,8), dummyr, (dummyi(i),i=9,17)
           read(10+nfile) dummyi(18), dummyr
           read(10+nfile) nmedia

           ! dummyi(1); trajectory starting point
           ! dummyi(2); means dt
           ! dummyi(3); trajectory end point

	   if (nfile .eq. 1) then
	      init_check=dummyi(1)
	      dt_check  =dummyi(2)
	      end_check =dummyi(3)
	      num_check =nmedia
	   else
	      ! Initial point check
	      if (init_check .ne. dummyi(1)) then
                 print*, "---------------------ERROR---------------------"
                 print*, "Trajectories should have the same initial point"
                 print*, "-----------------------------------------------"
                 stop
              endif

	      ! Trajectory interval check
	      if (dt_check .ne. dummyi(2)) then
                 print*, "---------------ERROR---------------"
                 print*, "Trajectory interval should be equal"
                 print*, "-----------------------------------"
                 stop
              endif

	      ! Initial point check
	      if (end_check .ne. dummyi(3)) then
                 print*, "-------------------ERROR-------------------"
                 print*, "Trajectories should have the same end point"
                 print*, "-------------------------------------------"
                 stop
              endif

	      ! Nmedia check
	      if (num_check .ne. nmedia) then
                 print*, "-------------ERROR-------------"
                 print*, "Number of atoms should be equal"
                 print*, "-------------------------------"
                 stop
              endif
	   endif
	enddo
	!---------------------------------------------------------------------------------------

	total_traj=int((dummyi(3)-dummyi(1))/dummyi(2))+1
	
	print*, "------------INFO-------------"
	print*, "Initial time    |", dummyi(1)
	print*, "End time        |", dummyi(3)
	print*, "Time interval   |", dummyi(2)
	print*, "Number of trajs |", total_traj
	print*, "Number of atoms |", nmedia
        print*, "-----------------------------"
	print*

        allocate(x(num_of_args,ntraj(num_of_args),nmedia),&
		 y(num_of_args,ntraj(num_of_args),nmedia),&
		 z(num_of_args,ntraj(num_of_args),nmedia),&
		 box(num_of_args,ntraj(num_of_args),6))

	! Reading Trajectory Data---------------------------------------------------------------

	do nfile=1,num_of_args
	   
	   print*, nfile, "th file READING.."
	
	   do i=1,ntraj(nfile)
	      read(10+nfile) (box(nfile,i,j), j=1,6) ! xx xy yy yz zx zz
	      read(10+nfile) (x(nfile,i,j), j=1,nmedia)
	      read(10+nfile) (y(nfile,i,j), j=1,nmedia)
	      read(10+nfile) (z(nfile,i,j), j=1,nmedia)
	   enddo

	enddo
	print*
        print*, "-----------------------------"
	print*, "(   COMPLETE!   )>   \(^O^)/"
        print*, "-----------------------------"
	print*
	!---------------------------------------------------------------------------------------

   ! -------------------------------------------------------------------------------------------	
   end subroutine dcd_reader

! ----------------------------------------------------------------------------------------------	
end module property_ftns
