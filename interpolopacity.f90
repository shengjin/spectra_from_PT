subroutine interpolopacity(press,temp,wavearray,kappa,kappaLayerNu,fCH4,fCO,fH2O,fCO2)
! Interpolate the opacities onto the wavearray grid
! kappaLayerNu, N_layer*N_nu

! CAUTION !!! HERE THE UNITY OF P IS BAR!!!
! CAUTION !!! HERE THE UNITY OF P IS BAR!!!
! CAUTION !!! HERE THE UNITY OF P IS BAR!!!
 
implicit none
include          'commonPTSpetc.h'

double precision,dimension(:),allocatable   :: press,temp,fCH4,fCO,fH2O,fCO2
integer                       ::  i,j,k
integer                       ::  N_layer,N_nu

double precision              ::  Tloc
integer                       ::  iT,iTup,iTlow,iP,iPup,iPlow
double precision              ::  kappaiTlowiPlow,tempiTlow,kappaiTupiPlow,tempiTup,kappaPlow
double precision              ::  kappaiTlowiPup,kappaiTupiPup,kappaPup
double precision              ::  pressiPlow,pressiPup,kappaSingleSP

double precision              ::  kappatemp

double precision,dimension(:),allocatable       ::  wavearray
double precision,dimension(:,:),allocatable     ::  kappaLayerNu
double precision,allocatable,dimension(:,:,:,:) ::  kappa ! nSpecies,nT,nP,nwave,kappa

logical,parameter             ::  DEBUG=.FALSE.
!logical,parameter            ::  DEBUG=.TRUE.


! load Press,Temp,Species,nP,nT,nSp,strP,strT,MW parameter
call PTSpetc

N_layer=size(press)
N_nu=size(wavearray)
        
        
!BILINEAR INTERPOLATION

!!!!! Get kappaLayerNu (kapp at each grid, N_layer*N_nu)
! using, ! press(N_layer),temp(N_layer),fH2O(N_layer),fCO2(N_layer),fCH4(N_layer),kappa(nSpecies,P,T,Wave)
!LATER: using, press(N_layer),temp(N_layer),fH2O(N_layer),fCO2(N_layer),fCH4(N_layer),fCO(N_layer),kappa(nSpecies,P,T,Wave)
!
!
! 1, For each layer
do i=1,N_layer
    !
    !case P outside tabulated area
    if(press(i)<=pressure(1))then
        write(*,*)'layerP < Pmin'
        write(*,*)'layerP,Pmin',press(i),pressure(1)
        stop
    end if
    !
    if(press(i)>=pressure(size(pressure)))then
        write(*,*)'layerP > Pmax'
        write(*,*)'layerP,Pmax',press(i),pressure(size(pressure))
        stop
    end if
    !
    !set Tloc = temp(i)
    Tloc=temp(i)
    !case T outside tabulated area
    if(Tloc<=temperature(1))then
        write(*,*)'layerT < Tmin'
        stop
    end if
    !
    if( Tloc > 2100.0) then
        write(*,*)'layerT > Tmax and layerT > 2100'
        stop
    else if( (Tloc>=temperature(size(temperature))) .and. (Tloc <= 2100.0))then
        Tloc=2000.0
        write(*,*)'iTlow,layerTemp(i),layerTloc,iTup',iTlow,temp(i),Tloc,iTup
    end if
    !
    ! get the P,T grid
    !
    iT=1
    do while (temperature(iT).lt.Tloc)
        iT=iT + 1
    end do
    iTup=iT
    iTlow=iT-1
    if(DEBUG)write(*,*)'iTlow,layerTemp(i),layerTloc,iTup',iTlow,temp(i),Tloc,iTup
    if(DEBUG)write(*,*)'temperature(iTlow),temprature(iTup)',temperature(iTlow),temperature(iTup)
    !
    iP=1
    do while (pressure(iP).lt.press(i))
        iP=iP + 1
    end do
    iPup=iP
    iPlow=iP-1
    if(DEBUG)write(*,*)'iPlow,layerP,iPup',iPlow,press(i),iPup
    if(DEBUG)write(*,*)'pressure(iPlow),pressure(iPup)',pressure(iPlow),pressure(iPup)
    !
    !
    !2,interpolate, loop for N_nu 
    !
    do j=1,N_nu
        !
        kappatemp=0D0
        !
        do k=1,nSpecies
            !
            !
            !3 For Each Species
            !
            !3(1) interpolate in Temperature at constant Pressure
            tempiTlow=temperature(iTlow)
            tempiTup=temperature(iTup)
            ! at iPlow
            kappaiTlowiPlow=kappa(k,iTlow,iPlow,j)
            kappaiTupiPlow=kappa(k,iTup,iPlow,j)
            kappaPlow=0D0
            call LinIntPol(1,Tloc,kappaiTlowiPlow,tempiTlow,kappaiTupiPlow,tempiTup,kappaPlow)
            ! at iPup
            kappaiTlowiPup=kappa(k,iTlow,iPup,j)
            kappaiTupiPup=kappa(k,iTup,iPup,j)
            kappaPup=0D0
            call LinIntPol(1,Tloc,kappaiTlowiPup,tempiTlow,kappaiTupiPup,tempiTup,kappaPup)
            !
            !3(2) interpolate in Pressure at constant Temperature
            pressiPlow=pressure(iPlow)
            pressiPup=pressure(iPup)
            call LinIntPol(1,press(i),kappaPlow,pressiPlow,kappaPup,pressiPup,kappaSingleSP)
            !
            !
            if ( k .EQ. 1) then
                kappatemp=kappatemp+kappaSingleSP*fH2O(i)
            else if( k .EQ. 2 ) then
                kappatemp=kappatemp+kappaSingleSP*fCO2(i)
            else if( k .EQ. 3 ) then
                kappatemp=kappatemp+kappaSingleSP*fCH4(i)
            else if( k .EQ. 4 ) then
                kappatemp=kappatemp+kappaSingleSP*fCO(i)
            else 
                write(*,*)'Error: Species number does not match'
                stop
            end if
            !
        end do
        !
        kappaLayerNu(i,j)=kappatemp
        !write(*,*)'i,j,kappaLayerNu(i,j)',i,j,kappaLayerNu(i,j)
        !
    end do
    !    
end do


if(DEBUG)then
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'nSpecies',nSpecies
    write(*,*)'N_nu',N_nu
    write(*,*)'N_layer',N_layer
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
end if


return 

end subroutine interpolopacity






subroutine LinIntPol(n,x,y1,x1,y2,x2,y)
!makes linear interpolation for f(x) (Dim n) at the point x, if (x1,f(x1)=y1) and (x2,f(x2)=y2)
!are given (dimension y1, y2 also n of course) using Lagrange's formula

IMPLICIT NONE
INTEGER                         ::n
DOUBLE PRECISION,DIMENSION(n)   ::y1,y2,y
DOUBLE PRECISION                ::x,x1,x2
LOGICAL,PARAMETER               ::debug=.FALSE.
IF(debug)THEN
   WRITE(*,*) 'Linintpol n,x,y1,y2,x1,x2,y',n,x,y1,y2,x1,x2,y
END IF

y=((x-x2)/(x1-x2))*y1+((x-x1)/(x2-x1))*y2

IF(debug)WRITE(*,*) 'y after',y

end subroutine LinIntPol


