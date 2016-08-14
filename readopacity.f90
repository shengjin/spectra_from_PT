subroutine readopacity(opacitydir,lamdarange,wavearray,kappa)
!read Opacity tables and creat a new one
 
implicit none
include          'commonPTSpetc.h'

character(LEN=200),INTENT(IN) ::  opacitydir
character(LEN=300)            ::  opacityfile
 
double precision              ::  lamdarange(2) !spectra range

integer                       ::  i,j,k,l ! do loop for nSp,nT,nP,nW
integer                       ::  totline ! total line in opacity table
integer                       ::  iwavebottomabove,iwavetopabove
double precision              ::  wavebottom,wavetop
double precision              ::  opactmp,wavetmp,bog
integer                       ::  n,m      ! do loop for nwave nkappa
integer                       ::  nWave,nFile


double precision,dimension(:),allocatable      ::  wavearray
double precision,allocatable,dimension(:,:,:,:)  :: kappa ! nSpecies,nT,nP,nwave,kappa

!logical,parameter             ::  DEBUG=.FALSE.
logical,parameter             ::  DEBUG=.TRUE.

!lamdarange(1)=C'H2O'
!lamdarange(2)=C'H2O'
! top: large and also in the top part of the opacity files
wavebottom=1.0D0/lamdarange(2)
wavetop=1.0D0/lamdarange(1)
if(DEBUG)then
    write(*,*)'wavebottom,lamdarange(2)',wavebottom,lamdarange(2)
    write(*,*)'wavetop,lamdarange(1)',wavetop,lamdarange(1)
end if


! load Press,Temp,Species,nP,nT,nSp,strP,strT,MW parameter
call PTSpetc


nFile=70 ! 70+ set file number of opacity tables

do i=1,nSpecies
    !
    do j=1,nTemp
        !
        do k=1,nPress
        !
        ! get the name of opacity file
        opacityfile=trim(opacitydir) // 'k_vs_nu-' // trim(species(i)) // &
                  '-T' // trim(strTemp(j))  // '-P' // trim(strPress(k)) // '.dat'
        !
        if(DEBUG)then
            write(*,*),pressure(k),strPress(k)
            write(*,*),temperature(j),strTemp(j)
            write(*,*),opacityfile
        end if
        !
        ! get the number of lines(opacity table) that are inside the spectra range 
        ! and allocate kappa, only for the first i,j,k=1,1,1 case
        if (i.EQ.1 .and. j.EQ.1 .and. k.EQ.1)then
            call getline(opacityfile,totline)
            if(ALLOCATED(wavearray))then
                deallocate(wavearray)
            end if
            allocate(wavearray(totline))
            open(unit=70,file=trim(opacityfile),form='FORMATTED',access='sequential',err=7000)
            do l=1,totline
                read(70,*),bog,wavearray(l)
                    if(wavearray(l)>wavetop)then
                        iwavetopabove=l
                    else if(wavearray(l)>wavebottom)then
                        iwavebottomabove=l
                    else 
                        goto 7010
                    end if
7010        end do
            close(70)
            nWave=iwavebottomabove - iwavetopabove + 1
            deallocate(wavearray)
            allocate(wavearray(nWave))
            if(ALLOCATED(kappa))then
                deallocate(kappa)
            end if
            allocate(kappa(nSpecies,nTemp,nPress,nWave))
        end if
        !
        if(DEBUG)then
            write(*,*)'Reading opacityfile',trim(opacityfile)
            write(*,*)'total line',totline
            write(*,*)'nSpecies,nTemp,nPress,nWave',nSpecies,nTemp,nPress,nWave
            write(*,*)'iwavetopabove,iwavebottomabove',iwavetopabove,iwavebottomabove
            write(*,*)'Size of wavearray',size(wavearray)
            write(*,*)'Size of kappa',size(kappa)
        end if
        !
        ! read the tables, read wave
        nFile=nFile+1
        open(unit=nFile,file=trim(opacityfile),form='formatted',access='sequential',err=7000)
        m=1
        do n=1,totline
            read(nFile,*),opactmp,wavetmp
            if( n >= iwavetopabove .and. n <= iwavebottomabove)then
                !kappa(i,j,k,m)=opactmp         [cm2/molecule]
                kappa(i,j,k,m)=opactmp/moleweight(i) ! [cm2/g]
                if(i.EQ.1 .and. j.EQ.1 .and. k.EQ.1)then
                    wavearray(m)=wavetmp ! read wavearray in the first
                end if
                m=m+1
            else if( n > iwavebottomabove)then
                goto 7020
            end if
7020    end do
        close(nFile)
        !
        if(DEBUG)then
        end if
        !
        end do 
        !
    end do
    !
end do


return 

7000  write(*,*)'ERROR opening opacity table:',trim(opacityfile)
      stop

end subroutine readopacity
