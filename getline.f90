subroutine getline(filename,line)
implicit none
character(LEN=200),INTENT(IN)     ::  filename
integer,INTENT(OUT)              ::  line           
integer,parameter                ::  MAXLINE=2000000  
integer                          ::  i,j   
 
open(unit=30,file=trim(filename),form='formatted',access='sequential')
j=0
do i=1,MAXLINE
read(unit=30,fmt=*,end=1030)
j=j+1
end do
1030 close(unit=30)
 
line=j
 
end subroutine getline
