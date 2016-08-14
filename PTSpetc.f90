subroutine PTSpetc
!read Opacity tables and creat a new one
 
implicit none
include          'commonPTSpetc.h'

logical,parameter             ::  DEBUG=.FALSE.
!logical,parameter             ::  DEBUG=.TRUE.


temperature(1)=500.0
temperature(2)=1000.0
temperature(3)=1500.0
temperature(4)=2000.0

strTemp(1)='500'
strTemp(2)='1000'
strTemp(3)='1500'
strTemp(4)='2000'

pressure(1)=0.0001
pressure(2)=0.001
pressure(3)=0.01
pressure(4)=0.1
pressure(5)=1.0
pressure(6)=10.0

strPress(1)='00001'
strPress(2)='0001'
strPress(3)='001'
strPress(4)='01'
strPress(5)='1'
strPress(6)='10'

! CAUTION!!! if change the order of each Species,
! Must modefy the interpolopacity.f90 also!!!
species(1)='H2O'
species(2)='CO2'
species(3)='CO'
species(4)='CH4'

! 1.6726d-24   ! Mass of proton [g]
moleweight(1)=1.6726d-24*18.0 ! [g/molecule]
moleweight(2)=1.6726d-24*44.0
moleweight(3)=1.6726d-24*28.0
moleweight(4)=1.6726d-24*16.0

IF(DEBUG)THEN
        write(*,*),species(1)
        write(*,*),species(2)
        write(*,*),species(3)
        write(*,*),species(4)

        write(*,*),pressure(1)
        write(*,*),pressure(2)
        write(*,*),pressure(3)
        write(*,*),pressure(4)
        write(*,*),pressure(5)
        write(*,*),pressure(6)

        write(*,*),strPress(1)
        write(*,*),strPress(2)
        write(*,*),strPress(3)
        write(*,*),strPress(4)
        write(*,*),strPress(5)
        write(*,*),strPress(6)
        
        write(*,*),moleweight(1)
        write(*,*),moleweight(2)
        write(*,*),moleweight(3)
        write(*,*),moleweight(4)
        
        write(*,*),temperature(1)
        write(*,*),temperature(2)
        write(*,*),temperature(3)
        write(*,*),temperature(4)
        
        write(*,*),strTemp(1)
        write(*,*),strTemp(2)
        write(*,*),strTemp(3)
        write(*,*),strTemp(4)
END IF

nTemp=size(temperature)
nPress=size(pressure)
nSpecies=size(species)

return 

end subroutine PTSpetc
