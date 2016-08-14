subroutine natconst
!     Natural constants in CGS
 
implicit none
include         'commonnatcon.h'

GG  = 6.672d-8       ! Gravitational constant
mp  = 1.6726d-24     ! Mass of proton          [g]
me  = 9.1095d-28     ! Mass of electron        [g]
kk  = 1.3807d-16     ! Bolzmann's constant     [erg/K]
hh  = 6.6262d-27     ! Planck's constant       [erg.s]
ee  = 4.8032d-10     ! Unit charge             
cc  = 2.9979d10      ! Light speed             [cm/s]
st  = 6.6524d-25     ! Thompson cross-section  [cm^2]
ss  = 5.6703d-5      ! Stefan-Boltzmann const  [erg/cm^2/K^4/s]
aa  = 7.5657d-15     ! 4 ss / cc               [erg/cm^3/K^4]
                                                               
!     Gas constants                                            
muh2= 2.3000d0       ! Mean molec weight H2+He+Metals
                                                               
!     Alternative units                                        
ev  = 1.6022d-12     ! Electronvolt            [erg]
kev = 1.6022d-9      ! Kilo electronvolt       [erg]
micron= 1.d-4        ! Micron                  [cm]
km  = 1.d5           ! Kilometer               [cm]
angstrom= 1.d-8      ! Angstroem               [cm]
                                                               
!     Astronomy constants                                      
LSun  = 3.8525d33      ! Solar luminosity        [erg/s]
!RSun  = 6.96d10       ! Solar radius            [cm]
RSun  = 6.96342*1e5*1e5 ! Solar radius            [cm]
Rjup  = 6.9911*1e4*1e5 ! Jupiter radius            [cm]
Mjup  = 1.8986e30     ! [g]
MSun  = 1.9891d33    ! Solar mass              [g]
!MSun  = 1.99d33       ! Solar mass              [g]
TSun  = 5780.d0        ! Solar temperature       [K]
AU    = 1.496d13       ! Astronomical Unit       [cm]
pc    = 3.08572d18     ! Parsec                  [cm]
                                                               
!     Time units                                               
year= 3.1536d7       ! Year                    [s]
hour= 3.6000d3       ! Hour                    [s]
day = 8.64d4         ! Day                     [s]
 
!     Math constants
pi  = 3.1415926535897932385d0 
 
return
 
end subroutine natconst
