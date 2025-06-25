! ----------------------------------------------------------------------------------------------	
module property_ftns

   !------------------------------------------------!
   !                                                !
   !              LAMMPS ANALYSIS PROGRAM           !
   !                   made by Jina                 !
   !						    !
   !          	       DCD  : 60 		    !
   !          	       MSD  : 225                   !
   !          	       NGP  : 332    	            !
   !          	       ICD  : 467                   !
   !          	       RDF  : 569		    !
   !		       GSRT :			    !
   !						    !
   !------------------------------------------------!

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
   integer, parameter			     :: crst_nmedia=4032
   integer, allocatable, dimension(:)        :: ntraj
					       !num_of_args : number of ensembles
   				               !nmedia : number of particles
   				               !total_traj : number of trajectory

   ! For LPS system
   integer				     :: num_Li, num_P, num_S	

   contains

! ----------------------------------------------------------------------------------------------	

   subroutine dcd_reader ! for ensembles

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

	if (mod(nmedia,8)==0) then
	   num_Li=nmedia*3/8
	   num_P=nmedia/8
	   num_S=nmedia/2
	else
	   print*, "There may be some defects in the system"
	   num_Li=crst_nmedia*3/8-(crst_nmedia-nmedia)/6
	   num_P=crst_nmedia/8-(crst_nmedia-nmedia)/6
	   num_S=crst_nmedia/2-(crst_nmedia-nmedia)*2/3
	endif

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

   subroutine mean_squared_displacement

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
	
	do i=1, int(dble(total_traj*0.7d0))
	  
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

   subroutine non_gaussian_parameter

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
	
   subroutine ionic_conductivity

	! Variables Declaration-----------------------------------------------------------------
        
        implicit none
 
        logical 		       :: file_exists
        character(20)		       :: input_line

        integer 		       :: ios, temp, last_step
	double precision	       :: msd, std, vol, const, diff_coeff, icd_avg, icd_std, &
					  last_msd, last_std		   ! Only for Li+ ions
	double precision, parameter    :: q=1.60217663d-19, &
					  kb=1.380649d-23,&   ! Constant of N-E eq.
					  A_to_cm=1.d-8, ps_to_s=1.d-12	   ! Optimized w/ 'ps'

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
	   call dcd_reader ! for volume
	else
	   print*, "There is no 'MSD.out' file."
	   print*, "Calling MSD program..."
	   call mean_squared_displacement
	endif

	! Input file
	open(101,file="MSD.out",status='old',action='read')
	! Output file
	open(103,file="ICD.out",status='unknown')

	! --------------------------------------------------------------------------------------
	print*, "Default = 300 K"
	write(*,'(A)',advance='no') " Enter the temperature : "
	read(*,'(A)',iostat=ios) input_line

	if (trim(input_line)=="") then
	   temp=300
	else
	   read(input_line, *, iostat=ios) temp
           if (ios /= 0) then
              print*, "Invalid input. Using default temperature."
              temp = 300
           endif
	endif

	print*, temp
	read(101,*)
	do
	   last_step=i ; last_msd=msd ; last_std=std 
 	   read(101,997) i, msd, std
	   if (i == 0) exit
	enddo

	! Volume ensemble average
	vol=0.d0
	do j=1, num_of_args
	   vol=vol+volume(j)%v(last_step) ! for last step
	enddo
	vol=vol/dble(num_of_args)*A_to_cm*A_to_cm*A_to_cm

	const=(dble(num_Li)/vol)*q*q/kb/dble(temp)     ! (C^2/J.cm^3)
	diff_coeff=last_msd/6d0/dble(last_step)        ! (A^2/ps)
	diff_coeff=diff_coeff*A_to_cm*A_to_cm/ps_to_s  ! (cm^2/s)
	icd_avg=const*diff_coeff		       ! (A/Jcm=S/cm)
	icd_std=const*diff_coeff/last_msd*last_std
	print*, const	
	write(103,'(A)') "# Ionic Conductivity     Standard deviation"
	write(103,*) icd_avg, icd_std
	
	print*, "Diffusion Coefficient : "
	print*, diff_coeff, "(cm2/s)"
	print*, "Ionic conductivity : "
	print*, icd_avg, "(S/m)"
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

   subroutine radial_distribution_function

	! Variables Declaration-----------------------------------------------------------------
        
        implicit none 
        integer :: bin, bin_num, idx1, idx2
        double precision :: dx, dy, dz, dist, bin_size, r
	double precision, parameter :: pi=datan(1.d0)*4.d0
	double precision, dimension(7) :: rdf_avg, rdf_std
	double precision, dimension(:), allocatable   :: median
	double precision, dimension(:,:), allocatable :: g, g_avg, rdf, square_rdf

	! Setting-------------------------------------------------------------------------------

	call dcd_reader

	open(104, file='RDF.out', status='unknown')

	write(104,*) "# r(A) Li_Li P_P S_S Li_P Li_S P_S Total"	
999	format(F10.3,14(F30.16)/)
998	format(A,I4,A)
	
	print*, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
	print*, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	print*
        print*, "-----------------------------"
	print*, "(   RDF start!   )>  (0o0 ;;)"
        print*, "-----------------------------"
	print*
	print*, "Enter the size/number of bins"
	print*, "        (Default : 0.05/1000)"
	
	! Default values
	bin_size=0.05d0
	bin=1000
	
	allocate(median(bin))
	allocate(g(7,bin), g_avg(7,bin), rdf(7,bin), square_rdf(7,bin))

	!---------------------------------------------------------------------------------------
	
	rdf=0.d0 ; square_rdf=0.d0

	do i=1, num_of_args
	   
	   write(*,998) " Processing ", i, " th ensemble"

	   g_avg=0.d0

	   do j=1, total_traj

	      if (mod(j,10)==0) write(*,998) " Processing ", j, " th trajectory"

	      g=0.d0
	
	      ! Li-Li RDF
	      do k=1, num_Li
		 do l=1, k-1

		    if (k==1) goto 1001
		    dx=x(i)%p(j,k)-x(i)%p(j,l)
		    dy=y(i)%p(j,k)-y(i)%p(j,l)
		    dz=z(i)%p(j,k)-z(i)%p(j,l)
		    dx=dx-nint(dx/box(i)%p(j,1))*box(i)%p(j,1)
		    dy=dy-nint(dy/box(i)%p(j,3))*box(i)%p(j,3)
		    dz=dz-nint(dz/box(i)%p(j,6))*box(i)%p(j,6)
		    dist=sqrt(dx*dx+dy*dy+dz*dz)

		    bin_num=int(dist/bin_size)+1
		    g(1,bin_num)=g(1,bin_num)+2.d0

		 enddo
1001	      enddo
	     
	      ! P-P RDF
	      do k=1, num_P
		 idx1=num_Li+k    ! To match P idx

		 do l=1, k-1

		    if (k==1) goto 1002
		    idx2=num_Li+l ! To match P idx
		    dx=x(i)%p(j,idx1)-x(i)%p(j,idx2)
		    dy=y(i)%p(j,idx1)-y(i)%p(j,idx2)
		    dz=z(i)%p(j,idx1)-z(i)%p(j,idx2)
		    dx=dx-nint(dx/box(i)%p(j,1))*box(i)%p(j,1)
		    dy=dy-nint(dy/box(i)%p(j,3))*box(i)%p(j,3)
		    dz=dz-nint(dz/box(i)%p(j,6))*box(i)%p(j,6)
		    dist=sqrt(dx*dx+dy*dy+dz*dz)

		    bin_num=int(dist/bin_size)+1
		    g(2,bin_num)=g(2,bin_num)+2.d0

		 enddo
1002	      enddo

	
	      ! S-S RDF
	      do k=1, num_S
		 idx1=num_Li+num_P+k    ! To match S idx
		 do l=1, k-1

		    if (k==1) goto 1003
		    idx2=num_Li+num_P+l ! To match S idx
		    dx=x(i)%p(j,idx1)-x(i)%p(j,idx2)
		    dy=y(i)%p(j,idx1)-y(i)%p(j,idx2)
		    dz=z(i)%p(j,idx1)-z(i)%p(j,idx2)
		    dx=dx-nint(dx/box(i)%p(j,1))*box(i)%p(j,1)
		    dy=dy-nint(dy/box(i)%p(j,3))*box(i)%p(j,3)
		    dz=dz-nint(dz/box(i)%p(j,6))*box(i)%p(j,6)
		    dist=sqrt(dx*dx+dy*dy+dz*dz)

		    bin_num=int(dist/bin_size)+1
		    g(3,bin_num)=g(3,bin_num)+2.d0

		 enddo
1003	      enddo

	      ! Li-P RDF
	      do k=1, num_Li
		 do l=1, num_P

		    idx1=num_Li+l ! To match P idx
		    dx=x(i)%p(j,k)-x(i)%p(j,idx1)
		    dy=y(i)%p(j,k)-y(i)%p(j,idx1)
		    dz=z(i)%p(j,k)-z(i)%p(j,idx1)
		    dx=dx-nint(dx/box(i)%p(j,1))*box(i)%p(j,1)
		    dy=dy-nint(dy/box(i)%p(j,3))*box(i)%p(j,3)
		    dz=dz-nint(dz/box(i)%p(j,6))*box(i)%p(j,6)
		    dist=sqrt(dx*dx+dy*dy+dz*dz)

		    bin_num=int(dist/bin_size)+1
		    g(4,bin_num)=g(4,bin_num)+1.d0

		 enddo
	      enddo

	      ! Li-S RDF
	      do k=1, num_Li
		 do l=1, num_S

		    idx1=num_Li+num_P+l ! To match S idx
		    dx=x(i)%p(j,k)-x(i)%p(j,idx1)
		    dy=y(i)%p(j,k)-y(i)%p(j,idx1)
		    dz=z(i)%p(j,k)-z(i)%p(j,idx1)
		    dx=dx-nint(dx/box(i)%p(j,1))*box(i)%p(j,1)
		    dy=dy-nint(dy/box(i)%p(j,3))*box(i)%p(j,3)
		    dz=dz-nint(dz/box(i)%p(j,6))*box(i)%p(j,6)
		    dist=sqrt(dx*dx+dy*dy+dz*dz)

		    bin_num=int(dist/bin_size)+1
		    g(5,bin_num)=g(5,bin_num)+1.d0

		 enddo
	      enddo

	
	      ! P-S RDF
	      do k=1, num_P
		 idx1=num_Li+k
		 do l=1, num_S

		    idx2=num_Li+num_P+l ! To match P idx
		    dx=x(i)%p(j,idx1)-x(i)%p(j,idx2)
		    dy=y(i)%p(j,idx1)-y(i)%p(j,idx2)
		    dz=z(i)%p(j,idx1)-z(i)%p(j,idx2)
		    dx=dx-nint(dx/box(i)%p(j,1))*box(i)%p(j,1)
		    dy=dy-nint(dy/box(i)%p(j,3))*box(i)%p(j,3)
		    dz=dz-nint(dz/box(i)%p(j,6))*box(i)%p(j,6)
		    dist=sqrt(dx*dx+dy*dy+dz*dz)

		    bin_num=int(dist/bin_size)+1
		    g(6,bin_num)=g(6,bin_num)+1.d0

		 enddo
	      enddo

	      do k=1, bin
	         g(1,k)=g(1,k)/dble(num_Li*num_Li)
	         g(2,k)=g(2,k)/dble(num_P*num_P)
	         g(3,k)=g(3,k)/dble(num_S*num_S)
	         g(4,k)=g(4,k)/dble(num_Li*num_P)
	         g(5,k)=g(5,k)/dble(num_Li*num_S)
	         g(6,k)=g(6,k)/dble(num_P*num_S)
		 g(7,k) = ( &				! TOTAL RDF
		    dble(num_Li*num_Li)*g(1,k) + &
		    dble(num_P*num_P)*g(2,k) + &
		    dble(num_S*num_S)*g(3,k) + &
		    dble(num_Li*num_P)*g(4,k) + &
		    dble(num_Li*num_S)*g(5,k) + &
 		   dble(num_P*num_S)*g(6,k) ) / dble(nmedia*nmedia)
	      enddo

	      g=g*volume(i)%v(j) ! division by number of atoms & number density = x volume

	      do k=1, bin
	         do l=1, 7
		     g_avg(l,k)=g_avg(l,k)+g(l,k)
		 enddo
	      enddo

	   enddo

	   do j=1, bin
	      median(j)=dble(j-0.5d0)*bin_size
	      r=dble(j-1)*bin_size
	
	      do k=1,7
	 	 if (j==1) then
		    g_avg(k,j)=0.d0
		 else
		    g_avg(k,j)=g_avg(k,j)/dble(total_traj)/(4.d0*pi*r*r*bin_size)
		 endif
	         rdf(k,j)=rdf(k,j)+g_avg(k,j)
		 square_rdf(k,j)=square_rdf(k,j)+g_avg(k,j)*g_avg(k,j)
	      enddo

	   enddo

	enddo

	do i=1, bin
	   do j=1,7
	      rdf_avg(j)=rdf(j,i)/dble(num_of_args)
	      rdf_std(j)=square_rdf(j,i)/dble(num_of_args)-rdf_avg(j)*rdf_avg(j)
	   enddo
	   write(104,999,advance='no') median(i), rdf_avg(1), rdf_std(1), rdf_avg(2), rdf_std(2),&
	      			      rdf_avg(3), rdf_std(3), rdf_avg(4), rdf_std(4), rdf_avg(5),&
	      			      rdf_std(5), rdf_avg(6), rdf_std(6), rdf_avg(7), rdf_std(7)
	enddo
		 	
	!---------------------------------------------------------------------------------------
	print*
	print*, "   Saved as 'RDF.out' file   "
        print*, "-----------------------------"
	print*, "(   COMPLETE!   )>    V(0_0)V"
        print*, "-----------------------------"
	print*

	deallocate(g,rdf,square_rdf)
	close(104)

   end subroutine radial_distribution_function

! ----------------------------------------------------------------------------------------------	
end module property_ftns
! ----------------------------------------------------------------------------------------------

program main

   use property_ftns

   implicit none
   integer :: input
   character(1) :: answer

123 format(A/)
   print*
   write(*,123,advance='no') "======================================================================================================"
   write(*,123,advance='no') "   __    ___    __  ___ __  ___ ___   ____    ___   _____ ___      ____ ____ __    ____ ____          " 
   write(*,123,advance='no') "  / /   / _ |  /  |/  //  |/  // _ \ / __/   / _ \ / ___// _ \    / __//  _// /   / __// __/          " 
   write(*,123,advance='no') " / /__ / __ | / /|_/ // /|_/ // ___/_\ \    / // // /__ / // /   / _/ _/ / / /__ / _/ _\ \            " 
   write(*,123,advance='no') "/____//_/ |_|/_/  /_//_/  /_//_/   /___/   /____/ \___//____/   /_/  /___//____//___//___/            " 
   write(*,123,advance='no') "      ___   ___   ____   ___   ____ ___  ________  __    ___   ___   ____   _____ ___   ___    __  ___"
   write(*,123,advance='no') "     / _ \ / _ \ / __ \ / _ \ / __// _ \/_  __/\ \/ /   / _ \ / _ \ / __ \ / ___// _ \ / _ |  /  |/  /"
   write(*,123,advance='no') "    / ___// , _// /_/ // ___// _/ / , _/ / /    \  /   / ___// , _// /_/ // (_ // , _// __ | / /|_/ / " 
   write(*,123,advance='no') "   /_/   /_/|_| \____//_/   /___//_/|_| /_/     /_/   /_/   /_/|_| \____/ \___//_/|_|/_/ |_|/_/  /_/  " 
   write(*,123,advance='no') " "
1  write(*,123,advance='no') "======================================================================================================"
                                                                                             
   print*
   print*, "          [1] Mean Squared Displacement"
   print*, "          [2] Non-Gaussian Parameter"
   print*, "          [3] Ionic Conductivity"
   print*, "          [4] Radial Distribution Function"
   print*, "          [5] Self-part van Hove Correlation Function"
   print*, "          [6] Distinct-part van Have Correlation Function"
   print*, "          [7] Hop Function"
   print*, "          [8] EXIT"
   print*

   write(*,'(A)',advance='no') "           [Enter the property index] : "
!   read(*,*) input
   input=4
   print*
   write(*,123,advance='no') "======================================================================================================"
   print*

   if (input==1) then
      call mean_squared_displacement
   elseif (input==2) then
      call non_gaussian_parameter
   elseif (input==3) then
      call ionic_conductivity
   elseif (input==4) then
      call radial_distribution_function
   elseif (input==5) then
      print*, "   Sorry, It's Updating... :( "
      print*
      goto 1
      !call self_part_van_hove_correlation_function()
   elseif (input==6) then
      print*, "   Sorry, It's Updating... :( "
      print*
      goto 1
      !call distinct_part_van_hove_correlation_function()
   elseif (input==7) then
      print*, "   Sorry, It's Updating... :( "
      print*
      goto 1
   elseif (input==8) then
      print*, "   PROGRAM ENDS ..."
      print*
      stop
   else
      print*, "   Enter the right index"
      goto 1
   end if
  
   print*
   write(*,123,advance='no') "======================================================================================================"
   print*
   write(*,'(A)',advance='no') "    Do you want to get another property? (Y/N) : "
  ! read(*,'(A)') answer
answer='n'
   if ((answer=='y') .or. (answer=='Y')) goto 1
  
   print*
   print*, "   PROGRAM ENDS ..."
   print*
   write(*,123,advance='no') "======================================================================================================"
   print*

end program

! ----------------------------------------------------------------------------------------------	
