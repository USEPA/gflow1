!! The original gflow1 code was exclusively compiled with the Lahey Fortran compiler.
!! This compiler provided a number of utilities which aren't supported by other
!! compilers. Fortunately, these utilities are generally available in Fortran90,
!! albeit with a different name or a slightly different form.
!! 
!! The utilities in lahey_utils.f90 wrap the standard Fortran90 intrinsics to provide
!! the same logic as the Lahey-specific intrinsics.
module lahey_utils

    implicit none
    private  ! Make everything private by default
    ! Explicitly export only what we want to be public
    public :: nblank
    public :: timer
    public :: getcl
    public :: iostat_msg

contains

    !> Returns length of string excluding trailing blanks (Lahey NBLANK compatible)
    function nblank(str) result(length)
        implicit none
        character(*), intent(in) :: str
        integer :: length

        length = len_trim(str)
    end

    !> Gets command line arguments (Lahey GETCL compatible)
    !!
    !! Arguments:
    !!     buffer (character(*)): Output buffer to store command line (256 char max)
    !!
    !! Combines all command line arguments into a single space-separated string.
    !! Stops with error if command line exceeds buffer size.
    subroutine getcl(buffer)
        implicit none
        character(*), intent(out) :: buffer
        integer :: i, nargs, stat, pos
        character(256) :: arg

        if (len(buffer) > 256) then
            write(*,*) 'GETCL: Buffer length exceeds 256 characters'
            stop
        end if

        buffer = ' '  ! Initialize to blanks
        pos = 1

        nargs = command_argument_count()
        do i = 1, nargs
            ! Get the argument
            call get_command_argument(i, arg, status=stat)
            if (stat /= 0) cycle

            ! Check if adding this argument would exceed buffer
            if (pos + len_trim(arg) + merge(1, 0, i > 1) > 256) then
                write(*,*) 'GETCL: Command line too long (max 256 characters)'
                stop
            end if

            ! Add space before argument (except for first one)
            if (i > 1) then
                buffer(pos:pos) = ' '
                pos = pos + 1
            end if

            ! Add argument
            buffer(pos:) = trim(arg)
            pos = pos + len_trim(arg)
        end do
    end

    !> Provides system timing functionality compatible with Lahey Fortran's timer
    !!
    !! This subroutine wraps the standard SYSTEM_CLOCK intrinsic to provide
    !! timing functionality compatible with legacy Lahey Fortran code. Returns
    !! time in centiseconds (1/100ths of a second) since system startup.
    subroutine timer(count_out)
        implicit none
        integer, intent(out) :: count_out
        integer :: count, count_rate
        call system_clock(count, count_rate)
        count_out = (count * 100) / count_rate
    end

    !> Get runtime I/O error message (Lahey IOSTAT_MSG compatible)
    !!
    !! Arguments:
    !!     iostat (integer): IOSTAT value from I/O statement
    !!     message (character(*)): Output buffer for error message (min 256 chars)
    !!
    !! Provides descriptive messages for common I/O errors.
    !! Based on standard Fortran IOSTAT values.
    subroutine iostat_msg(iostat, message)
        implicit none
        integer, intent(in) :: iostat
        character(*), intent(out) :: message

        message = ' '  ! Initialize to blanks
        select case (iostat)
            case (0)
                message = "No error"
            case (-1)
                message = "End of file"
            case (-2)
                message = "End of record"
            case (5000:5999)  ! Common range for OS-specific errors
                message = "Operating system error"
            case (1)
                message = "Invalid file specification"
            case (2)
                message = "File not found"
            case (3)
                message = "Invalid path specification"
            case (4)
                message = "Too many files open"
            case (5)
                message = "File access denied"
            case (6)
                message = "Invalid file handle"
            case (7)
                message = "Memory allocation error"
            case (8)
                message = "Invalid file format"
            case (9)
                message = "Invalid file access"
            case (10)
                message = "Invalid record specification"
            case (11)
                message = "Record length mismatch"
            case (12)
                message = "Record not found"
            case (13)
                message = "Invalid record length"
            case (14)
                message = "Write protection violation"
            case (15)
                message = "Device I/O error"
            case default
                write(message, '(A,I0)') "Unknown I/O error. IOSTAT = ", iostat
        end select
    end

end module