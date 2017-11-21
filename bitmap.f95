module bitmap
use params
implicit none

type :: Fileheader
    integer(8) :: Size
    integer(8) :: OffBits
    integer(8) :: Width
    integer(8) :: Height
    integer(4) :: BitCount
end type

type(FileHeader) :: Header
character, dimension(:,:), ALLOCATABLE :: Bmp
integer(2), dimension(:,:), ALLOCATABLE :: Img

Contains

  Subroutine load_bmp
    use params
    implicit none
    integer :: i,j,nw
    integer(4) :: read_short
    integer(8) :: read_long
    if (BMPLOADED .eq. 1 ) then
        return
    end if
    
    !write(6,*) "Reading: ",pot_filename
    OPEN(UNIT=11, FILE=pot_filename, STATUS="OLD", ACCESS="STREAM")

    READ(11, POS=3) read_long
    Header%Size = read_long
    READ(11, POS=11) read_long
    Header%OffBits = read_long
    READ(11, POS=19) read_long
    Header%Width = read_long
    READ(11, POS=23) read_long
    Header%Height = read_long
    READ(11, POS=29) read_short
    Header%BitCount = read_short

    if (Header%Width .ne. NX+1 .or. Header%Height .ne. NY+1) then
        write(6,*) "ERROR WRONG BMP SIZE"
        stop
    end if

    nw = INT(Header%Width*Header%BitCount/8)
    if (MODULO(nw,4) .eq. 0) then
        ALLOCATE(Bmp(0:nw-1,0:NY))
    else
        ALLOCATE(Bmp(0:nw+4-MODULO(nw,4)-1,0:NY))
    end if
    READ(11, POS=1+Header%OffBits) Bmp

    ALLOCATE(Img(-NX/2:NX/2,-NY/2:NY/2))

    do i = 0,NX
    do j = 0,NY
        Img(i-NX/2,j-NY/2) = INT((ICHAR(Bmp(i*Header%BitCount/8,j))&
                             +ICHAR(Bmp(i*Header%BitCount/8+1,j))&
                             +ICHAR(Bmp(i*Header%BitCount/8+2,j)))/3,kind=2)
    end do
    end do

    CLOSE(UNIT=11)
    BMPLOADED = 1
  End Subroutine 

 end module
