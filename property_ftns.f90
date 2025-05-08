! ----------------------------------------------------------------------------------------------	
module property_ftns

   ! Variables Declaration----------------------------------------------------------------------	
	
   implicit none
        
   ! For loop
   integer :: i, j, k, l, nfile
   
   ! File I/O system
   real                                      :: dummyr
   integer, dimension(20)                    :: dummyi(20)
   character(4)                              :: dummyc
   character(80), allocatable, dimension(:)  :: filename
   integer                                   :: num_of_args

   ! Checking Trajectory variables
   integer 	    			     :: init_check, dt_check, end_check, num_check

   ! Data variables
   type 				     :: traj_ptr
      real, pointer, dimension(:,:)	     :: p
   end type
   type 				     :: box_ptr
      real(8), pointer, dimension(:,:)	     :: p
      real(8), pointer, dimension(:)	     :: v
   end type

   type(traj_ptr), allocatable, dimension(:) :: x, y, z
   type(box_ptr), allocatable, dimension(:)  :: box, volume
   integer                                   :: nmedia, total_traj, atom_idx
   integer, allocatable, dimension(:)        :: ntraj
					       !num_of_args : number of ensembles
   				               !nmedia : number of particles
   				               !total_traj : number of trajectory

   ! For LPS system
   integer				     :: num_Li, num_P, num_S	

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
	allocate(x(num_of_args), y(num_of_args), z(num_of_args), box(num_of_args), volume(num_of_args))

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
		    box(i)%p(total_traj,6),&
		    volume(i)%v(total_traj))
	enddo

	num_Li=nmedia*3/8
	num_P=nmedia/8
	num_S=nmedia/2

	print*
	print*, "---------LPS System----------"
	print*, " Number of Li+ | ", num_Li
	print*, " Number of P5+ | ", num_P
	print*, " Number of S2- | ", num_S
        print*, "-----------------------------"

	! Reading Trajectory Data---------------------------------------------------------------

	do nfile=1,num_of_args
	   
	   print*, nfile, "th file READING.."
	
	   do i=1,total_traj
	      read(10+nfile) (box(nfile)%p(i,j), j=1,6) ! xx xy yy yz zx zz
	      read(10+nfile) (x(nfile)%p(i,j), j=1,nmedia)
	      read(10+nfile) (y(nfile)%p(i,j), j=1,nmedia)
	      read(10+nfile) (z(nfile)%p(i,j), j=1,nmedia)
	      volume(nfile)%v(i)=box(nfile)%p(i,1)*box(nfile)%p(i,3)*box(nfile)%p(i,6)
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

   subroutine mean_squared_displacement()

	! Variables Declaration-----------------------------------------------------------------
        
        implicit none
 
        double precision	       :: dx, dy, dz
	double precision, dimension(4) :: displ, msd, square_msd, msd_avg, msd_std

	! Setting-------------------------------------------------------------------------------

	call dcd_reader()

	open(101, file='MSD.out', status='unknown')

	write(101,*) "# time(ps) Li Li_std P P_std S S_std Total Total_std"	

999	format(I10,F22.16,F22.16,F22.16,F22.16,F22.16,F22.16,F22.16,F22.16/)
998	format(A,I4,A)

	print*, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
	print*, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	print*
        print*, "-----------------------------"
	print*, "(   MSD start!   )> \(o_o ;;)"
        print*, "-----------------------------"
	print*

	!---------------------------------------------------------------------------------------
	
	do i=1, int(total_traj*0.7)
	  
	   if (mod(i,10)==0) write(*,998) " Processing ", i, " th trajectory"

	   msd=0.d0 ; square_msd=0.d0

	   do j=1, num_of_args
	   
	      displ=0.d0	     

	      do k=1,total_traj-i

		 ! MSD of Li
		 do l=1, num_Li
		    dx=x(j)%p(i+k,l)-x(j)%p(k,l)
		    dy=y(j)%p(i+k,l)-y(j)%p(k,l)
		    dz=z(j)%p(i+k,l)-z(j)%p(k,l)
		    displ(1)=displ(1)+dx*dx+dy*dy+dz*dz
		 enddo

		 ! MSD of P
		 do l=1, num_P
		    atom_idx=num_Li+l
		    dx=x(j)%p(i+k,atom_idx)-x(j)%p(k,atom_idx)
		    dy=y(j)%p(i+k,atom_idx)-y(j)%p(k,atom_idx)
		    dz=z(j)%p(i+k,atom_idx)-z(j)%p(k,atom_idx)
		    displ(2)=displ(2)+dx*dx+dy*dy+dz*dz
		 enddo

		 ! MSD of S
		 do l=1, num_S
		    atom_idx=num_Li+num_P+l
		    dx=x(j)%p(i+k,atom_idx)-x(j)%p(k,atom_idx)
		    dy=y(j)%p(i+k,atom_idx)-y(j)%p(k,atom_idx)
		    dz=z(j)%p(i+k,atom_idx)-z(j)%p(k,atom_idx)
		    displ(3)=displ(3)+dx*dx+dy*dy+dz*dz
		 enddo

		 ! MSD of total
		 displ(4)=displ(1)+displ(2)+displ(3)

	      enddo
	       
	      displ(1)=displ(1)/dble(total_traj-i)/dble(num_Li)
	      displ(2)=displ(2)/dble(total_traj-i)/dble(num_P)
	      displ(3)=displ(3)/dble(total_traj-i)/dble(num_S)
	      displ(4)=displ(4)/dble(total_traj-i)/dble(nmedia)

	      do k=1,4
	         msd(k)=msd(k)+displ(k)
		 square_msd(k)=square_msd(k)+displ(k)*displ(k)
	      enddo

	   enddo
	   
	   do j=1,4
	      msd_avg(j)=msd(j)/dble(num_of_args)
	      msd_std(j)=dsqrt(square_msd(j)/dble(num_of_args)-msd_avg(j)*msd_avg(j))
	   enddo
				      
	   write(101,999,advance='no') i, msd_avg(1), msd_std(1), &
 				       msd_avg(2), msd_std(2), &
				       msd_avg(3), msd_std(3), &
				       msd_avg(4), msd_std(4)

	enddo      

	print*
	print*, "   Saved as 'MSD.out' file   "
        print*, "-----------------------------"
	print*, "(   COMPLETE!   )>    \(^O^)/"
        print*, "-----------------------------"
	print*
	close(101)

   end subroutine mean_squared_displacement

! ----------------------------------------------------------------------------------------------	

   subroutine non_gaussian_parameter()

	! Variables Declaration-----------------------------------------------------------------
        
        implicit none
 
        double precision	       :: dx, dy, dz
	double precision, dimension(4) :: displ, msd, square_msd, msd_avg, msd_std, &
					  square_displ, mqd, square_mqd, mqd_avg, mqd_std, &
					  ngp, ngp_std
					  ! quartic displacement for NGP

	! Setting-------------------------------------------------------------------------------

	call dcd_reader()

999	format(I10,F22.16,F22.16,F22.16,F22.16,F22.16,F22.16,F22.16,F22.16/)
998	format(A,I4,A)

	open(102, file='NGP.out', status='unknown')

	write(102,*) "# time(ps) Li Li_std P P_std S S_std Total Total_std"	

	print*, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
	print*, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	print*
        print*, "-----------------------------"
	print*, "(   NGP start!   )> \(o_o ;;)"
        print*, "-----------------------------"
	print*

	!---------------------------------------------------------------------------------------
	
	do i=1, int(total_traj/3)
	  
	   if (mod(i,10)==0) write(*,998) " Processing ", i, " th trajectory"

	   msd=0.d0 ; square_msd=0.d0 ! for MSD
	   mqd=0.d0 ; square_mqd=0.d0 ! for NGP

	   do j=1, num_of_args
	   
	      displ=0.d0	      ! for MSD
 	      square_displ=0.d0	      ! for NGP

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
				      
	   do j=1,4
	      ngp(j)=0.6d0*mqd_avg(j)/msd_avg(j)/msd_avg(j)-1.d0
	      ngp_std(j)=dsqrt((-1.2d0*mqd_avg(j)*msd_std(j)/msd_avg(j)/msd_avg(j)/msd_avg(j))**2.d0&
			 +(0.6d0*mqd_std(j)/msd_avg(j)/msd_avg(j))**2.d0)
	   enddo

	   write(102,999,advance='no') i, ngp(1), ngp_std(1), ngp(2), ngp_std(2), &
					  ngp(3), ngp_std(3), ngp(4), ngp_std(4)

	enddo      

	print*
	print*, "   Saved as 'NGP.out' file   "
        print*, "-----------------------------"
	print*, "(   COMPLETE!   )>    \(^O^)/"
        print*, "-----------------------------"
	print*
	close(102)

   end subroutine non_gaussian_parameter

! ----------------------------------------------------------------------------------------------
	
   subroutine ionic_conductivity()

	! Variables Declaration-----------------------------------------------------------------
        
        implicit none
 
        logical 		       :: file_exists

        integer 		       :: ios, temp, last_step
	double precision	       :: msd, std, icd_avg, icd_std, vol, var, &
					  last_msd, last_std		   ! Only for Li+ ions
	double precision, parameter    :: const=1.6d0*1.6d0/1.38d2/6d0     ! Constant of N-E eq.
									   ! Optimized w/ 'ps'

	! Setting-------------------------------------------------------------------------------

999	format(I10,F22.16,F22.16,F22.16,F22.16,F22.16,F22.16,F22.16,F22.16/)
997	format(I10,F22.16,F22.16)

	print*, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
	print*, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	print*
        print*, "-----------------------------"
	print*, "(   ICD start!   )> \(o_o ;;)"
        print*, "-----------------------------"
	print*

	! Reading MSD.out
	inquire(file="MSD.out", exist=file_exists)

        if (file_exists) then
	   call dcd_reader() ! for volume
	else
	   print*, "There is no 'MSD.out' file."
	   print*, "Calling MSD program..."
	   call mean_squared_displacement()
	endif

	! Input file
	open(101,file="MSD.out",status='old',action='read')
	! Output file
	open(103,file="ICD.out",status='unknown')

	! --------------------------------------------------------------------------------------
	print*, "Default = 300 K"
	write(*,'(A)',advance='no') " Enter the temperature : "
	read(*,'(A)',iostat=ios) temp

	if (ios .ne. 0) temp=300

	do 
 	   read(101,997,iostat=ios) last_step, msd, std
	   if (ios .ne. 0) exit
	enddo

	print*, last_step
	print*, last_msd
	print*, last_std

	! Volume ensemble average
	vol=0.d0
	do j=1, num_of_args
	   vol=vol+volume(j)%v(last_step) ! for last step
	enddo
	vol=vol/dble(num_of_args)

	print*, vol	
	var=const*dble(num_Li)/dble(temp)/vol/dble(last_step)
	print*, var
	icd_avg=var*last_msd ! S/cm
	icd_std=var*last_std
	
	write(103,'(A)') "# Ionic Conductivity     Standard deviation"
	write(103,*) icd_avg, icd_std

	print*
	print*, "   Saved as 'ICD.out' file   "
        print*, "-----------------------------"
	print*, "(   COMPLETE!   )>    \(^O^)/"
        print*, "-----------------------------"
	print*
	close(103)

	! --------------------------------------------------------------------------------------

   end subroutine ionic_conductivity

! ----------------------------------------------------------------------------------------------	

   subroutine radial_distribution_function()

	! Variables Declaration-----------------------------------------------------------------
        
        implicit none 
        integer :: bin
        double precision :: dx, dy, dz, distance, bin_size, median, r
	double precision, parameter :: pi=datan(1.d0)*4.d0
	double precision, dimension(7) :: rdf_avg, rdf_std
	double precision, dimension(:,:), allocatable :: g, rdf, square_rdf

	! Setting-------------------------------------------------------------------------------

	call dcd_reader()

	open(103, file='RDF.out', status='unknown')

	write(103,*) "# r(A) Li_Li P_P S_S Li_P Li_S P_S Total"	
999	format(F10.3,F30.16,F30.16,F30.16,F30.16,F30.16,F30.16,F30.16,F30.16,&
		     F30.16,F30.16,F30.16,F30.16,F30.16,F30.16/)

	allocate(g(7,bin), rdf(7,bin), square_rdf(7,bin)) ! 수정 필요!!!!!!!!!!!!!!!!!!!

	! Default values
	bin_size=0.05d0
	bin=1000

	print*, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
	print*, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	print*
        print*, "-----------------------------"
	print*, "(   RDF start!   )>  (0o0 ;;)"
        print*, "-----------------------------"
	print*
	print*, "Enter the size/number of bins"
	print*, "        (Default : 0.05/1000)"
	read(*,*) bin_size, bin
	print*, bin_size, bin
	
	!---------------------------------------------------------------------------------------
	
	!do i=1, total_traj
	
	 
	   
	
	!---------------------------------------------------------------------------------------
	print*
	print*, "   Saved as 'RDF.out' file   "
        print*, "-----------------------------"
	print*, "(   COMPLETE!   )>    V(0_0)V"
        print*, "-----------------------------"
	print*
	close(103)

   end subroutine radial_distribution_function

! ----------------------------------------------------------------------------------------------	
end module property_ftns
