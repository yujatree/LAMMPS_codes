! ----------------------------------------------------------------------------------------------	
module property_ftns

   ! Variables Declaration----------------------------------------------------------------------	
	
   implicit none
        
   ! For loop
   integer :: i, j, k, l, nfile
   
   ! File I/O system
   real                                     :: dummyr
   integer, dimension(20)                   :: dummyi(20)
   character(4)                             :: dummyc
   character(80), allocatable, dimension(:) :: filename
   integer                                  :: num_of_args

   ! Checking Trajectory variables
   integer 	    :: init_check, dt_check, end_check, num_check

   ! Data variables
   type 				     :: traj_ptr
      real, pointer, dimension(:,:)	     :: p
   end type
   type 				     :: box_ptr
      real(8), pointer, dimension(:,:)	     :: p
   end type

   type(traj_ptr), allocatable, dimension(:) :: x, y, z
   type(box_ptr), allocatable, dimension(:)  :: box
   integer                                   :: nmedia, total_traj
   integer, allocatable, dimension(:)        :: ntraj
					       !num_of_args : number of ensembles
   				               !nmedia : number of particles
   				               !total_traj : number of trajectory

   contains

! ----------------------------------------------------------------------------------------------	

   subroutine dcd_reader() ! for ensembles

   !------------------------------------------------!
   !                                                !
   !                Main Link Program               !
   !              made by Park HyungShick           !
   !                modified by EunHyeok            !
   !                modified by Jina    	    !
   !				                    !
   !               Ver.1 2020.06.26                 !
   !               Ver.2 2022.01.03 		    !
   !   		   Ver.3 2025.04.27                 !
   !				                    !
   !------------------------------------------------!

   ! Refer to JENOGMIN KIM READ_DCD program
   !          Mar. 18, 2013

   
   ! Argument-----------------------------------------------------------------------------------	

        num_of_args=iargc() ! # of files

        if (num_of_args .eq. 0) then
           print*, "------------ERROR------------"
           print*, "USAGE : ./run.x dcd1 dcd2 ..."
           print*, "-----------------------------"
           stop
        endif

	allocate(filename(num_of_args), ntraj(num_of_args))
	allocate(x(num_of_args), y(num_of_args), z(num_of_args), box(num_of_args))

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
	
	do i=1,num_of_args
           allocate(x(i)%p(total_traj,nmedia),&
		    y(i)%p(total_traj,nmedia),&
		    z(i)%p(total_traj,nmedia),&
		    box(i)%p(total_traj,6))
	enddo

	! Reading Trajectory Data---------------------------------------------------------------

	do nfile=1,num_of_args
	   
	   print*, nfile, "th file READING.."
	
	   do i=1,total_traj
	      read(10+nfile) (box(nfile)%p(i,j), j=1,6) ! xx xy yy yz zx zz
	      read(10+nfile) (x(nfile)%p(i,j), j=1,nmedia)
	      read(10+nfile) (y(nfile)%p(i,j), j=1,nmedia)
	      read(10+nfile) (z(nfile)%p(i,j), j=1,nmedia)
	   enddo

	enddo

	print*
        print*, "-----------------------------"
	print*, "( READ COMPLETE! )>   \(^O^)/"
        print*, "-----------------------------"
	print*

	!---------------------------------------------------------------------------------------

   end subroutine dcd_reader

! ----------------------------------------------------------------------------------------------	

   subroutine mean_squared_displacement() ! with Non-Gaussian Parameter

	! Variables Declaration-----------------------------------------------------------------
        
        implicit none 
        integer :: num_Li, num_P, num_S, atom_idx
        double precision :: dx, dy, dz

	double precision, dimension(4) :: displ, msd, square_msd, msd_avg, msd_std, &
					  square_displ, mqd, square_mqd, mqd_avg, mqd_std, &
					  ngp, ngp_std
					  ! quartic displacement for NGP

	! Setting-------------------------------------------------------------------------------

	call dcd_reader()

	open(101, file='MSD.out', status='unknown')
	open(102, file='NGP.out', status='unknown')

	write(101,*) "# time(ps) Li Li_std P P_std S S_std Total Total_std"	
999	format(I10,F22.16,F22.16,F22.16,F22.16,F22.16,F22.16,F22.16,F22.16/)
998	format(A,I4,A)

	num_Li=nmedia*3/8
	num_P=nmedia/8
	num_S=nmedia/2

	print*, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
	print*, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	print*
        print*, "-----------------------------"
	print*, "(   MSD start!   )> \(o_o ;;)"
        print*, "-----------------------------"
	print*
	print*, "------------INFO-------------"
	print*, " Number of Li+ | ", num_Li
	print*, " Number of P5+ | ", num_P
	print*, " Number of S2- | ", num_S
        print*, "-----------------------------"
	print*

	!---------------------------------------------------------------------------------------
	
	do i=1, int(total_traj/3)
	  
	   if (mod(i,10)==0) write(*,998) " Processing ", i, " th trajectory"

	   msd=0.d0 ; square_msd=0.d0
	   mqd=0.d0 ; square_mqd=0.d0

	   do j=1, num_of_args
	   
	      displ=0.d0
 	      square_displ=0.d0

	      do k=1,total_traj-i

		 ! MSD of Li
		 do l=1, num_Li
		    dx=x(j)%p(i+k,l)-x(j)%p(k,l)
		    dy=y(j)%p(i+k,l)-y(j)%p(k,l)
		    dz=z(j)%p(i+k,l)-z(j)%p(k,l)
		    displ(1)=displ(1)+dx*dx+dy*dy+dz*dz
		    square_displ(1)=square_displ(1)&
		                   +(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)
		 enddo

		 ! MSD of P
		 do l=1, num_P
		    atom_idx=num_Li+l
		    dx=x(j)%p(i+k,atom_idx)-x(j)%p(k,atom_idx)
		    dy=y(j)%p(i+k,atom_idx)-y(j)%p(k,atom_idx)
		    dz=z(j)%p(i+k,atom_idx)-z(j)%p(k,atom_idx)
		    displ(2)=displ(2)+dx*dx+dy*dy+dz*dz
		    square_displ(2)=square_displ(2)&
		                   +(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)
		 enddo

		 ! MSD of S
		 do l=1, num_S
		    atom_idx=num_Li+num_P+l
		    dx=x(j)%p(i+k,atom_idx)-x(j)%p(k,atom_idx)
		    dy=y(j)%p(i+k,atom_idx)-y(j)%p(k,atom_idx)
		    dz=z(j)%p(i+k,atom_idx)-z(j)%p(k,atom_idx)
		    displ(3)=displ(3)+dx*dx+dy*dy+dz*dz
		    square_displ(3)=square_displ(3)&
		                   +(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)
		 enddo

		 ! MSD of total
		 displ(4)=displ(1)+displ(2)+displ(3)
		 square_displ(4)=square_displ(1)+square_displ(2)+square_displ(3)

	      enddo
	       
	      displ(1)=displ(1)/dble(total_traj-i)/dble(num_Li)
	      displ(2)=displ(2)/dble(total_traj-i)/dble(num_P)
	      displ(3)=displ(3)/dble(total_traj-i)/dble(num_S)
	      displ(4)=displ(4)/dble(total_traj-i)/dble(nmedia)
	      square_displ(1)=square_displ(1)/dble(total_traj-i)/dble(num_Li)
	      square_displ(2)=square_displ(2)/dble(total_traj-i)/dble(num_P)
	      square_displ(3)=square_displ(3)/dble(total_traj-i)/dble(num_S)
	      square_displ(4)=square_displ(4)/dble(total_traj-i)/dble(nmedia)

	      do k=1,4
	         msd(k)=msd(k)+displ(k)
		 square_msd(k)=square_msd(k)+displ(k)*displ(k)
		 mqd(k)=mqd(k)+square_displ(k)
		 square_mqd(k)=square_mqd(k)+square_displ(k)*square_displ(k)
	      enddo

	   enddo
	   
	   do j=1,4
	      msd_avg(j)=msd(j)/dble(num_of_args)
	      msd_std(j)=dsqrt(square_msd(j)/dble(num_of_args)-msd_avg(j)*msd_avg(j))
	      mqd_avg(j)=mqd(j)/dble(num_of_args)
	      mqd_std(j)=dsqrt(square_mqd(j)/dble(num_of_args)-mqd_avg(j)*mqd_avg(j))
	   enddo
				      
	   write(101,999,advance='no') i, msd_avg(1), msd_std(1), &
 				       msd_avg(2), msd_std(2), &
				       msd_avg(3), msd_std(3), &
				       msd_avg(4), msd_std(4)
	  
	   do j=1,4
	      ngp(j)=0.6d0*mqd_avg(j)/msd_avg(j)/msd_avg(j)-1.d0
	      ngp_std(j)=dsqrt((-1.2d0*mqd_avg(j)*msd_std(j)/msd_avg(j)/msd_avg(j)/msd_avg(j))**2.d0&
			 +(0.6d0*mqd_std(j)/msd_avg(j)/msd_avg(j))**2.d0)
	   enddo

	   write(102,999,advance='no') i, ngp(1), ngp_std(1), ngp(2), ngp_std(2), &
					  ngp(3), ngp_std(3), ngp(4), ngp_std(4)

	enddo      

	print*
	print*, "   Saved as 'MSD.out' file   "
	print*, "            'NGP.out' file   "
        print*, "-----------------------------"
	print*, "(   COMPLETE!   )>    \(^O^)/"
        print*, "-----------------------------"
	print*
	close(101)
	close(102)

   end subroutine mean_squared_displacement

! ----------------------------------------------------------------------------------------------	
end module property_ftns
