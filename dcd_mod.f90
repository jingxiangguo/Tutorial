module dcd_reader
    use iso_c_binding, only: c_int, c_float, c_double, c_char
    implicit none 

    private

    ! accessible by other Fortran module or routines
    public :: readdcdheader, read_xyz_box

contains

    subroutine readdcdheader(dcdfile,strleng,total_atoms,total_frames) bind(c, name="call_dcd_header")
        implicit none

        ! Passed
        integer(c_int),intent(in) :: strleng
        character(kind=c_char,len=1),intent(in),dimension(1:strleng) :: dcdfile

        ! Local
        character(len=:),allocatable :: fortran_dcdfile_name 
        integer,dimension(20) :: ICNTRL
        character(len=3) :: fext
        character(len=4) :: HDR
        real(4) :: delta
        integer :: i,unitnum

        ! Return 
        integer(c_int), intent(out) :: total_atoms, total_frames

        unitnum = 12345

        call convert_c_string_f_string(dcdfile,strleng,fortran_dcdfile_name)

        open(unit=unitnum,file=trim(fortran_dcdfile_name),form='unformatted',status='old')

            read(unitnum) HDR, (ICNTRL(i),i=1,9),delta,(ICNTRL(i),i=11,20)
            read(unitnum)
            read(unitnum) total_atoms

            total_frames =  ICNTRL(1)

        close(unitnum)

        end subroutine

    subroutine read_xyz_box(dcdfile, strleng, total_atoms, current_frame, xyz, box) bind(c, name="call_dcd_traj")
        implicit none
        ! Passed
        integer(c_int),intent(in) :: strleng,total_atoms,current_frame
        character(kind=c_char,len=1),intent(in),dimension(1:strleng) :: dcdfile
    
        ! Local
        integer ::i,iframe,imol,counter,unitnum
        real(8),dimension(6) :: XTLABC
        character(len=:),allocatable :: fortran_dcdfile_name

        ! Return
        real(c_double), intent(out), dimension(1:3) :: box
        real(c_float), intent(out), dimension(1:3,1:total_atoms) :: xyz

        call convert_c_string_f_string(dcdfile, strleng, fortran_dcdfile_name)
        
        unitnum = 12345

        open(unit=unitnum,file=trim(fortran_dcdfile_name),form='unformatted',status='old')

            read(unitnum)
            read(unitnum)
            read(unitnum)

            do i = 1, current_frame - 1

                read(unitnum)
                read(unitnum)
                read(unitnum)
                read(unitnum)

            end do

            read(unitnum) (XTLABC(i),i=1,6)

            read(unitnum) (xyz(1,imol),imol=1,total_atoms)
            read(unitnum) (xyz(2,imol),imol=1,total_atoms)
            read(unitnum) (xyz(3,imol),imol=1,total_atoms)
           
            box =  [ XTLABC(1),XTLABC(3),XTLABC(6) ]
            
        close(unitnum)

        end subroutine

    subroutine convert_c_string_f_string(str,strlen,f_string)
        implicit none
        integer ,intent(in):: strlen
        character(len=1),intent(in),dimension(1:strlen) :: str

        character(len=:),allocatable,intent(out) :: f_string
        integer :: i

        f_string = ""

        do i = 1,strlen

            f_string = f_string//str(i)

        end do

        end subroutine

    
end module 
