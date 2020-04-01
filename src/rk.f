module rungekutta

use parameters
use constants
use derivatives
use analysis


implicit none

private


public rk

contains


subroutine rk(y0)

 
!Arguments
real(kind=dp),intent(IN),dimension(entries) :: y0 !initial conditions

!Other
real(kind=dp), dimension(size(y0)) :: y, y1,dy !Some useful vectors used during integration
real(kind=dp), dimension(nrows,ncols) :: AllData !Big array to save all data. 12 coordinates + tau +
real(kind=dp), dimension(:,:),allocatable :: output !smaller array which will be outout
integer(kind=dp) :: i,j,NSteps !,nsteps !index for saving to array


real(kind=dp) :: tau
real(kind=dp) :: mm, xC, yC, zC !Cartesian components



!Set the integration tolerance
if (dp .EQ. 8) then
escal = 1.0d15
else if (dp .EQ. 16) then
escal = 1.0d19
endif


!Assign y0 to vecor y
y  = y0
tau = 0.0_dp

!Save the first row to array
i = 1
AllData(i,1:12) = y
AllData(i,13) = tau !tau


!Integrate
do while ( y(1) .LT. time_cutoff )
!do while (i .LT. 5)  

    !Update
    call RKF(y,y1)
    y = y1
 


    !Print statements
!    print *, y(1)/ time_cutoff, h, y(2), i


   
    !Save the output
    i = i + 1
    AllData(i,1:12) = y
    tau = tau + h
    AllData(i,13) = tau




enddo
NSteps = i




print *, 'Total number of steps is = ', i
print *, 'Runge Kutta completed. Start data I/O'
!!!!!!!!!! - Save the output for analysis and plotting - !!!!!!!
!!!!!!!!!! - Save the output for analysis and plotting - !!!!!!!
!!!!!!!!!! - Save the output for analysis and plotting - !!!!!!!



!Binary format. See discussion at https://stackoverflow.com/questions/24395686/best-way-to-write-a-large-array-to-file-in-fortran-text-vs-other



!First reallocate to create a smaller array
allocate(output(NSteps,ncols))
output = AllData(1:i, :)




print *, 'savefile1'
!Savefile 1
open(unit=30,file=savefile1,status='replace',form='formatted')
do j=1,NSteps

!Convert to Cartesian first

        mm = sqrt(output(j,2)**2 + a**2)
        xC = mm * sin(output(j,3)) * cos(output(j,4))
        yC = mm * sin(output(j,3)) * sin(output(j,4))
        zC = mm * cos(output(j,3)) 

write(30,*) output(j,1),xC,yC,zC
enddo
close(30)


print *, 'savefile2'
!Savefile 2
open(unit=30,file=savefile2,status='replace',form='formatted')
do j=1,i
write(30,*) output(j,1)
enddo
close(30)


    







end subroutine rk


end module rungekutta
