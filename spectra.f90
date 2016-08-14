program spectra
!To calculate the emergent spectra of a planet. 
!The structure of the planet is described by "atmocompo.dat" and  "structure.dat"
!Sheng.Jin 
 
implicit none
include                'commonnatcon.h'

double precision             ::  press_bottom,press_top 
character(LEN=200)           ::  filename(2)
character(LEN=200)           ::  opacitydir
 
double precision             :: lamdarange(2)
integer                      :: n,lamdadivide ! Divide grid-space to sub-grid
double precision             :: lamdarangeonce(2) ! The lamdarange in each sub-grid

double precision,dimension(:),allocatable ::  press,temp,gacc,fCH4,fCO,fH2O,fCO2

double precision                                :: Tstar,Rstar
double precision                                :: Rplanet
double precision                                :: Tempcomparison

double precision,dimension(:),allocatable       :: wavearray
double precision,dimension(:,:),allocatable     :: wavearyspskap 
double precision,allocatable,dimension(:,:,:,:) :: kappa
 
logical,parameter            ::  DEBUG=.TRUE.
!logical,parameter           ::  DEBUG=.FALSE.
 

interface
    subroutine readstructure(filename,press_bottom,press_top,press,temp,gacc,fCH4,fCO,fH2O,fCO2)
    character(LEN=200)           ::  filename(2)
    double precision             ::  press_bottom,press_top 
    double precision,dimension(:),allocatable ::  press,temp,gacc,fCH4,fCO,fH2O,fCO2
    end subroutine readstructure
end interface


interface
    subroutine readopacity(opacitydir,lamdarangeonce,wavearray,kappa)
    character(LEN=200)           ::  opacitydir
    double precision             ::  lamdarangeonce(2)
    double precision,dimension(:),allocatable       :: wavearray
    double precision,allocatable,dimension(:,:,:,:) :: kappa
    end subroutine readopacity
end interface


interface
    subroutine calcspectrum(press,temp,gacc,fCH4,fCO,fH2O,fCO2,wavearray,kappa,Tstar,Rstar,Rplanet,Tempcomparison)
    double precision                                :: Tstar,Rstar
    double precision                                :: Rplanet
    double precision                                :: Tempcomparison
    double precision,dimension(:),allocatable ::  press,temp,gacc,fCH4,fCO,fH2O,fCO2
    double precision,dimension(:),allocatable       :: wavearray
    double precision,allocatable,dimension(:,:,:,:) :: kappa
    end subroutine calcspectrum
end interface


call natconst

!!!!!!!!!!!!!!!!!!!    SET INITIAL PARAMETER

! Set the stellar and planet parameter
Tstar=TSun
Rplanet=1.0*Rjup              ! Jupiter
!Rplanet=(9.4492/11.209)*Rjup   ! Saturn
Rstar=1.0*RSun

! Set Pressure range
press_bottom=10.0  !bar
!press_top=0.0001      !bar
press_top=9.8693e-4      !bar

! Set Spectra(lamda) range, short, long
lamdarange(1)=0.4*micron  !cm
lamdarange(2)=28.0125*micron   !cm

! Divide the whole spetra range to N (N=lamdadivide) section, in order to
! reduce the memory consumption. Could simply set to 1 for wavelength grid
! less than 100,000
lamdadivide=3

! Atmos structure file
filename(1)='atmocompo.dat'
filename(2)='structure.dat'

! Opacity dir
!opacitydir='/home/jin/yamila_opacAll/'
opacitydir='/home/jin/yamila_opac1of10/'
!opacitydir='/home/jin/yamila_opac10mean/'

! For comparison, Black body intensity (in output file) at a certain temperature
Tempcomparison=2000.0

!!!!!!!!!!!!!!!       SET INITIAL PARAMETER END
 

! Clear existing "emergentflux.dat"
open(unit=90,file='emergentflux.dat')
write(90,*)
close(90)
!Nine columns in "emergentflux.dat": wavelength(micro meter), intensity(emergent flux), Fplanet/Fstar(emergent flux)
!Bnu(bottom layer), Bnu(top layer), Bnu(Tempcomparison))
!Fplanet(only_bottom_layer)/Fstar, Fplanet(only_top_layer/Fstar, Flux(Tempcomparison)/Fstar



do n=1,lamdadivide

    ! Get the spectra range in each calculation
    lamdarangeonce(1)=exp(log(lamdarange(1))-real(n-1)*(log(lamdarange(1))-log(lamdarange(2)))/real(lamdadivide))
    lamdarangeonce(2)=exp(log(lamdarange(1))-real(n)*(log(lamdarange(1))-log(lamdarange(2)))/real(lamdadivide))

    if(DEBUG)then
        write(*,*)'lamdadivide,n',lamdadivide,n
        write(*,*)'lamdarangeonce(1)',lamdarangeonce(1)
        write(*,*)'lamdarangeonce(2)',lamdarangeonce(2)
    end if

    ! Read the Atmos Structure that in the Pressure range
    call readstructure(filename,press_bottom,press_top,press,temp,gacc,fCH4,fCO,fH2O,fCO2)
 

    ! Read Yamila's opacity
    call readopacity(opacitydir,lamdarangeonce,wavearray,kappa)

    if(DEBUG)then
        write(*,*)'End of reading files'
        write(*,*)'Number of Layers:', size(press)
        write(*,*)'size(wavearray)',size(wavearray)
        write(*,*)'size(kappa)/size(wavearray)',size(kappa)/size(wavearray)
    end if

    call calcspectrum(press,temp,gacc,fCH4,fCO,fH2O,fCO2,wavearray,kappa,Tstar,Rstar,Rplanet,Tempcomparison)

end do


stop ! Normal End
 
end program spectra

