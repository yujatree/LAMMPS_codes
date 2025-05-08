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
   print*, "   [1] Mean Squared Displacement"
   print*, "   [2] Non-Gaussian Parameter"
   print*, "   [3] Ionic Conductivity"
   print*, "   [4] Radial Distribution Function"
   print*, "   [5] Self-part van Hove Correlation Function"
   print*, "   [6] Distinct-part van Have Correlation Function"
   print*, "   [7] EXIT"
   print*

   write(*,'(A)',advance='no') "    [Enter the index of the property] : "
   read(*,*) input
   
   print*
   write(*,123,advance='no') "======================================================================================================"
   print*

   if (input==1) then
      call mean_squared_displacement()
   elseif (input==2) then
      call non_gaussian_parameter() 
   elseif (input==3) then
      call ionic_conductivity()
   elseif (input==4) then
      print*, "   Sorry, It's Updating... :( "
      print*
      !call radial_distribution_function()
      goto 1
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
   read(*,'(A)') answer

   if ((answer=='y') .or. (answer=='Y')) goto 1
  
   print*, "   PROGRAM ENDS ..."
   print*
   write(*,123,advance='no') "======================================================================================================"
   print*

end program
