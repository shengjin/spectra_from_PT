subroutine readstructure(filename,press_bottom,press_top,press,temp,gacc,fCH4,fCO,fH2O,fCO2)
!read P,T,G,and fraction of wanted species in the P range
 
implicit none
character(LEN=200),INTENT(IN) ::  filename(2)
character(LEN=200)           ::  filenamenew
double precision             ::  press_bottom,press_top
integer,parameter            ::  MAXLINE=50000        ! maximum layer number
logical,parameter            ::  DEBUG=.FALSE.


integer                      ::  i,j
integer                      ::  nLayers,ipresstopbelow,ipressbottomabove
double precision             ::  bog
!logical,parameter            ::  WRITENEWFILE=.FALSE.
logical,parameter            ::  WRITENEWFILE=.TRUE.
double precision,dimension(:),allocatable  ::  fCH4tmp,fCOtmp,fH2Otmp,fCO2tmp,presstmp,temptmp,gacctmp
double precision,dimension(:),allocatable       ::  fCH4,fCO,fH2O,fCO2,press,temp,gacc
 
filenamenew='atmocomponew.dat'


if(ALLOCATED(presstmp))then
    deallocate(presstmp)
end if
if(ALLOCATED(temptmp))then
    deallocate(temptmp)
end if
if(ALLOCATED(gacctmp))then
    deallocate(gacctmp)
end if
if(ALLOCATED(fCH4tmp))then
    deallocate(fCH4tmp)
end if
if(ALLOCATED(fCOtmp))then
    deallocate(fCOtmp)
end if
if(ALLOCATED(fH2Otmp))then
    deallocate(fH2Otmp)
end if
if(ALLOCATED(fCO2tmp))then
    deallocate(fCO2tmp)
end if
allocate(presstmp(MAXLINE))
allocate(temptmp(MAXLINE))
allocate(gacctmp(MAXLINE))
allocate(fCH4tmp(MAXLINE))
allocate(fCOtmp(MAXLINE))
allocate(fH2Otmp(MAXLINE))
allocate(fCO2tmp(MAXLINE))

open(unit=20,file=trim(filename(1)),form='formatted',access='sequential',err=1020)
open(unit=19,file=trim(filename(2)),form='formatted',access='sequential',err=1030)
!press(bar),temp,gacc,fH2,fH,fHe,fCH4,fCH3,fCO,fH2O,fCO2,fC2H2,fN2,fNH3,fHCN,fOH,fO,fHCHO,fTiO,fTiO2,fNO   
!
if(WRITENEWFILE)then
    open(unit=21,file=trim(filenamenew))
end if
! 
do i=1,MAXLINE       ! read loop
    !
    !set the the line begin,line end that in pressure range
    read(20,*)presstmp(i),temptmp(i),bog,bog,bog,fCH4tmp(i),bog,fCOtmp(i),fH2Otmp(i),fCO2tmp(i)
    read(19,*)bog,bog,bog,bog,bog,bog,bog,bog,bog,bog,bog,bog,bog,gacctmp(i)
    if( presstmp(i) > press_bottom )then
        ipressbottomabove=i
    else if( presstmp(i) <= press_top )then
        ipresstopbelow=i
        goto 1050
    end if
    !
end do      
!
1050  close(20) !1050  rewind(20)
      close(19) !1050  rewind(20)
!
!nLayers=(ipresstopbelow-1)-(ipressbottomabove+1)+1
nLayers=ipresstopbelow-1-ipressbottomabove
!
if(ALLOCATED(press))then
    deallocate(press)
end if
if(ALLOCATED(temp))then
    deallocate(temp)
end if
if(ALLOCATED(gacc))then
    deallocate(gacc)
end if
if(ALLOCATED(fCH4))then
    deallocate(fCH4)
end if
if(ALLOCATED(fCO))then
    deallocate(fCO)
end if
if(ALLOCATED(fH2O))then
    deallocate(fH2O)
end if
if(ALLOCATED(fCO2))then
    deallocate(fCO2)
end if
allocate(press(nLayers))
allocate(temp(nLayers))
allocate(gacc(nLayers))
allocate(fCH4(nLayers))
allocate(fCO(nLayers))
allocate(fH2O(nLayers))
allocate(fCO2(nLayers))
!
do j=1, nLayers
    !
    press(j)=presstmp(ipressbottomabove+j)
    temp(j)=temptmp(ipressbottomabove+j)
    gacc(j)=gacctmp(ipressbottomabove+j)
    fCH4(j)=fCH4tmp(ipressbottomabove+j)
    fCO(j)=fCOtmp(ipressbottomabove+j)
    fH2O(j)=fH2Otmp(ipressbottomabove+j)
    fCO2(j)=fCO2tmp(ipressbottomabove+j)
    !
    if(WRITENEWFILE)then
        write(21,'(200ES12.4)')press(j),temp(j),gacc(j),fCH4(j),fCO(j),fH2O(j),fCO2(j)
    end if
    !
end do
 

if(WRITENEWFILE)then
    close(21)
end if
 

return

1020 write(*,*)'Error: file not find:',trim(filename(1))
1030 write(*,*)'Error: file not find:',trim(filename(2))
stop

!1030 write(*,*)'Warning: $Filename could not cover the pressure range',trim(filename)
!stop

end subroutine readstructure

