module pair_correlation
    use iso_c_binding, only: c_int, c_float, c_double 
    implicit none 
    private 

    real(8), parameter :: pi = 3.14159265358979323d0

    public :: build_homo_pair_distance_histogram, & 
            & normalize_histogram 

contains

    subroutine build_homo_pair_distance_histogram(total_atoms, &
                                                 & cutoff, &
                                                 & num_bins, &
                                                 & XYZ, &
                                                 & L_xyz, &
                                                 & rdf_histogram, &
                                                 & volume) bind(c,name="call_homo_pair_dist_his")
        implicit none
        
        ! Passed 
        integer(c_int), intent(in) :: total_atoms,num_bins 
        real(c_double), intent(in) :: cutoff
        real(c_double), intent(in),dimension(1:3,1:total_atoms) :: XYZ
        real(c_double), intent(in),dimension(1:3) :: L_xyz

        ! Local 
        real(8),dimension(1:size(L_xyz)) :: xyz_temp,xyz_separate
        integer :: i,j,bin_index
        real(8) :: distance_sqr,distance,cutoff_sqr,r_interval 

        ! Output 
        real(c_double), intent(out),dimension(1:num_bins) :: rdf_histogram
        real(c_double), intent(out) :: volume 

        r_interval = cutoff/num_bins 

        cutoff_sqr = cutoff*cutoff

        rdf_histogram = 0.0d0 

        volume = product(L_xyz)

        do i = 1, total_atoms-1

            xyz_temp = xyz(1:3,i)

            do j = i + 1,total_atoms

                ! compute the separation vector between i,j  

                xyz_separate = xyz_temp - xyz(1:3,j)

                ! apply minimum image convention 

                xyz_separate = xyz_separate - L_xyz*dnint(xyz_separate/L_xyz)

                ! sum of squared separation vector  

                distance_sqr = xyz_separate(1)*xyz_separate(1) & 
                           & + xyz_separate(2)*xyz_separate(2) & 
                           & + xyz_separate(3)*xyz_separate(3)

                if (distance_sqr < cutoff_sqr ) then 
                    
                    bin_index = int(dsqrt(distance_sqr)/r_interval) + 1
                    
                    rdf_histogram(bin_index)  =  rdf_histogram(bin_index) + 2

                end if

            end do 

        end do

        end subroutine 

    subroutine normalize_histogram(rdf_histogram, &
                                 & num_bins,&
                                 & cutoff,&
                                 & natoms,&
                                 & num_configs,&
                                 & bulk_density,&
                                 & gr) bind(c,name="call_normalize_histogram")
        implicit none  

        ! Passed 
        integer(c_int), intent(in) :: num_bins,natoms,num_configs
        real(c_double), intent(in),dimension(1:num_bins) :: rdf_histogram 
        real(c_double), intent(in) :: cutoff,bulk_density  

        ! Local
        real(8) :: preconst,vshell_i,half_interval,upper_r,lower_r,center_r,r_interval  
        integer :: i 

        ! Output
        real(c_double),dimension(1:num_bins),intent(out) :: gr  

        preconst = 4.0d0/3.0d0*pi

        half_interval = r_interval/2.0d0*natoms 

        r_interval = cutoff/num_bins

        do i = 1,num_bins 

            upper_r = i*r_interval 

            lower_r = (i-1)*r_interval

            center_r = half_interval + r_interval*(i-1)

            vshell_i = preconst*( ( upper_r*upper_r*upper_r ) - (lower_r*lower_r*lower_r) )  

            gr(i) = rdf_histogram(i)/( vshell_i*bulk_density*natoms*num_configs)

        end do 				

        end subroutine 

end module 
