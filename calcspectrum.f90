subroutine calcspectrum(press,temp,gacc,fCH4,fCO,fH2O,fCO2,wavearray,kappa,Tstar,Rstar,Rplanet,Tempcomparison)

implicit none
!include                'commonPTSpetc.h'
include                'commonnatcon.h'

double precision,dimension(:),allocatable       :: press,temp,gacc,fCH4,fCO,fH2O,fCO2
double precision,dimension(:),allocatable       :: wavearray
double precision,allocatable,dimension(:,:,:,:) :: kappa !Sp,P,T,Wave:kappa

double precision                                :: Tstar,Rstar
double precision                                :: Rplanet
double precision,dimension(:),allocatable       :: intStar !N_nu 1/wavearray
double precision,dimension(:),allocatable       :: ratio   !N_nu 1/wavearray

double precision                                :: mu
integer                                         :: N_nu,N_layer
integer                                         :: i,j,k
double precision,dimension(:),allocatable       :: dP  !N_layer
double precision,dimension(:),allocatable       :: dz  !N_layer
double precision,dimension(:),allocatable       :: z  !N_layer
double precision,dimension(:),allocatable       :: rho  !N_layer, density
double precision,dimension(:),allocatable       :: sigma !N_layer, surface density

double precision,dimension(:),allocatable       :: nu !N_nu 1/wavearray

double precision,dimension(:,:),allocatable     :: intensity  !N_Layer*N_nu: I
double precision,dimension(:,:),allocatable     :: dtau    !N_layer*N_nu Tau
double precision,dimension(:,:),allocatable     :: Srcfunc !N_layer*N_nu: Source function
double precision,dimension(:,:),allocatable     :: Bnu     !N_layer*N_nu, Planck function

! Black body intensity at a certain T, for comparison
double precision                                :: Tempcomparison

double precision                                :: layerP,layerT
double precision                                :: bplanck
double precision,dimension(:,:),allocatable     :: wavearyspskap !Sp,Nu
double precision                                :: kappatmp
double precision,dimension(:,:),allocatable     :: kappaLayerNu !kapp at each grid, N_layer*N_nu

logical,parameter                               :: DEBUG=.TRUE.

interface
    subroutine interpolopacity(press,temp,wavearray,kappa,kappaLayerNu,fCH4,fCO,fH2O,fCO2)
    double precision,dimension(:),allocatable       :: press,temp,fCH4,fCO,fH2O,fCO2
    double precision,dimension(:),allocatable       :: wavearray
    double precision,dimension(:,:),allocatable     :: kappaLayerNu
    double precision,allocatable,dimension(:,:,:,:) :: kappa
    end subroutine interpolopacity
end interface


! Get nature constant
call natconst


N_layer=size(press)
N_nu=size(wavearray)


if(ALLOCATED(dP))then
    deallocate(dP)
end if
if(ALLOCATED(dz))then
    deallocate(dz)
end if
if(ALLOCATED(z))then
    deallocate(z)
end if
if(ALLOCATED(rho))then
    deallocate(rho)
end if
if(ALLOCATED(sigma))then
    deallocate(sigma)
end if
allocate(dP(N_layer))
allocate(dz(N_layer))
allocate(z(N_layer))
allocate(rho(N_layer))
allocate(sigma(N_layer))

if(ALLOCATED(nu))then
    deallocate(nu)
end if
if(ALLOCATED(intStar))then
    deallocate(intStar)
end if
if(ALLOCATED(ratio))then
    deallocate(ratio)
end if
if(ALLOCATED(intensity))then
    deallocate(intensity)
end if
if(ALLOCATED(dtau))then
    deallocate(dtau)
end if
if(ALLOCATED(Srcfunc))then
    deallocate(Srcfunc)
end if
if(ALLOCATED(Bnu))then
    deallocate(Bnu)
end if
if(ALLOCATED(kappaLayerNu))then
    deallocate(kappaLayerNu)
end if
allocate(nu(N_nu))
allocate(intStar(N_nu))
allocate(ratio(N_nu))
!intensity(N_layer+1,N_nu), intensity(2,Nu) corresponding to dtua(1,N_nu),Srcfunc(1,N_nu
allocate(intensity(N_layer+1,N_nu))
allocate(dtau(N_layer,N_nu))
allocate(Srcfunc(N_layer,N_nu))
allocate(Bnu(N_layer,N_nu))
allocate(kappaLayerNu(N_layer,N_nu))


!!!!! Get kappaLayerNu, kapp at each grid, N_layer*N_nu
!
call interpolopacity(press,temp,wavearray,kappa,kappaLayerNu,fCH4,fCO,fH2O,fCO2)


!!!!!! Get dP, press, N_layer
! the first layer is at the bottom
do i=1,N_layer-1
    dP(i)=press(i)-press(i+1)
end do
    dP(N_layer)=dP(N_layer-1)
! convert bar to bayre
do i=1,N_layer
    dP(i)=dP(i)*1.01325E6  
    !write(*,*)'dP(i),i',dP(i),i
    press(i)=press(i)*1.01325E6
end do

!!!!! Get dz, N_layer
mu=2.3*mp !mean molecular weight, for mixed gases in exoplanetary atmosphere
z(N_layer)=0.0
do j=1,N_layer-1
    i=N_layer-j
    dz(i)=(dP(i+1)/press(i+1))*(kk*temp(i+1))/(mu*gacc(i+1)) ! [cm]
    !write(*,*)'i,dz(i)',i,dz(i)
    z(i)=z(i+1)-dz(i)
end do

do i=1,N_layer-1
    dz(i)=z(i)-z(i+1)
    !write(*,*)'i,dz(i)',i,dz(i)
end do
    dz(N_layer)=dz(N_layer-1)
    !write(*,*)'dz(N_layer),N_layer',dz(N_layer),N_layer

!!!!! Get rho,  N_layer
do i=1,N_layer
    rho(i)=dP(i)/(gacc(i)*dz(i))  ! [g/cm3]  density
    !write(*,*)'i,rho(i)',i,rho(i)  
end do

!!!!! Get Sigma, N_layer
do i=1,N_layer
    sigma(i)=rho(i)*dz(i)     ! [g/cm2]  surface density
    !write(*,*)'i,sigma(i)',i,sigma(i)  
end do

!get nu(N_nu)
do i=1,N_nu
    nu(i)=cc*wavearray(i)
    !write(*,*)'nu(i),wavearray(i)',nu(i),wavearray(i)
end do

! Get dtau, N_layer*N_Nu
do i=1,N_layer
    do j=1,N_nu
        dtau(i,j)=sigma(i)*kappaLayerNu(i,j)
        !write(*,*)'i,j,sigma(i),kappaLayerNu(i,j),dtau(i,j)',i,j,sigma(i),kappaLayerNu(i,j),dtau(i,j)  
    end do
end do
 
!get Bnu(N_layer,N_nu)
do i=1,N_layer
    do j=1,N_nu
        Bnu(i,j)=bplanck(temp(i),nu(j))
        !write(*,*)'Bnu(i,j)',Bnu(i,j)
    end do
end do

!get Source function
do i=1,N_layer
    do j=1,N_nu
        Srcfunc(i,j)=Bnu(i,j)*(1.0-exp(-dtau(i,j)))
        !write(*,*)'i,j,Bnu(i,j),dtau(i,j),Srcfunc(i,j)',i,j,Bnu(i,j),dtau(i,j),Srcfunc(i,j)
    end do
end do

! Calc the intensity at bottom 
do j=1,N_nu
    !intensity(1,j),bottom; intensity(N_layer+1,j),top
    intensity(1,j)=Bnu(1,j)
end do

do i=1,N_layer
    do j=1,N_nu
        !intensity(N_layer+1,N_nu), intensity(2,Nu) corresponding to dtua(1,N_nu),Srcfunc(1,N_nu
        intensity(i+1,j)=intensity(i,j)*exp(-dtau(i,j))+Srcfunc(i,j)
    end do
end do

! get the stellar intensity and flux ratio
do j=1,N_nu
    intStar(j)=bplanck(Tstar,nu(j))
    ratio(j)=(intensity(N_layer+1,j)/intStar(j))*(Rplanet**2.0/Rstar**2.0)
end do

write(*,*)'End of calc, begin writting ...'

open(unit=90,file='emergentflux.dat',ACCESS='APPEND')

! write the intensity at the top layer
do j=1,N_nu
    write(90,*),10000/wavearray(j),intensity(N_layer+1,j),ratio(j),   &
                Bnu(1,j),Bnu(N_layer,j),bplanck(Tempcomparison,nu(j)),        &
                Bnu(1,j)/intStar(j)*(Rplanet**2.0/Rstar**2.0),        &
                Bnu(N_layer,j)/intStar(j)*(Rplanet**2.0/Rstar**2.0),  &
                bplanck(Tempcomparison,nu(j))/intStar(j)*(Rplanet**2.0/Rstar**2.0)
                !wavelength(micro meter), intensity(emergent flux), Fplanet/Fstar(emergent flux)  
                !Bnu(bottom layer), Bnu(top layer), Bnu(Tempcomparison))
                !Fplanet(only_bottom_layer)/Fstar, Fplanet(only_top_layer/Fstar, Flux(Tempcomparison)/Fstar
end do

close(90)

if(DEBUG)then
    write(*,*)'Layer number, Nu number',N_layer,N_nu
    write(*,*)'Size of intensity(size(wavearray))', size(intensity)
end if



end subroutine calcspectrum




function bplanck(temp,nu)
  implicit none
  double precision :: temp
  double precision :: nu
  double precision :: bplanck
  if(temp.eq.0.d0) then 
     bplanck = 0.d0
     return
  endif
  bplanck = 1.47455d-47 * nu * nu * nu /                    &
            (exp(4.7989d-11 * nu / temp)-1.d0) + 1.d-290
  return
end function bplanck

