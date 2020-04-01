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
real(kind=dp), dimension(size(y0)) :: y, y1,dy
!an array to store all the data. !Should be faster than dynamically allocating
!However this does require more exploration
!see https://stackoverflow.com/questions/8384406/how-to-increase-array-size-on-the-fly-in-fortran
real(kind=dp), dimension(nrows,size(y0)+2) :: AllData !Big array to save all data. 12 coordinates + tau + dt/dtau
real(kind=dp), dimension(nrows,4) :: DerivativesStore !Big array to save all data. 12 coordinates + tau + dt/dtau
real(kind=dp), dimension(:,:),allocatable :: output !smaller array which will be outout
integer(kind=dp) :: i,j !,nsteps !index for saving to array
real(kind=dp) :: mm, xC, yC, zC !Cartesian components
real(kind=dp) :: tau, roemer 
real(kind=dp) :: r,theta,phi,S1,S2,S3,Sx,Sy,Sz, thetaSL, phiSL
real(kind=dp) :: time_cutoff
real(kind=dp),dimension(4) :: uVector, xVector, OVector, OVector_comoving
real(kind=dp) :: top,bot, critical_phase







y  = y0


!Save the first row to array
i = 1
AllData(i,1:12) = y
AllData(i,13) = tau
 

call derivs(y,dy)
DerivativesStore(i,:) = dy(1:4)


AllData(i,14) = dy(1)





time_cutoff = 2.0_dp*PI*semi_major**(3.0_dp/2.0_dp) * N_orbit
!Integrate
!do while ( abs( y(4) - y0(4)) .LT. 1.10_dp*N_orbit*2.0_dp*PI)    
do while ( y(1) .LT. time_cutoff )
  

    !Update
    call RKF(y,y1)
    y = y1
 


    !Print statements
 !   print *, y(1)/ time_cutoff, h


   
    !Save the output
    i = i + 1
    AllData(i,1:12) = y
    tau = tau + h
    AllData(i,13) = tau

    !And calculate some derivative info - probably not the most effective way to do this

    call derivs(y,dy)

    
    DerivativesStore(i,:) = dy(1:4)


    



    AllData(i,14) = dy(1)

 !   print *, y


enddo






print *, 'Runge Kutta completed. Start data I/O'
!!!!!!!!!! - Save the output for analysis and plotting - !!!!!!!
!!!!!!!!!! - Save the output for analysis and plotting - !!!!!!!
!!!!!!!!!! - Save the output for analysis and plotting - !!!!!!!



!Binary format. See discussion at https://stackoverflow.com/questions/24395686/best-way-to-write-a-large-array-to-file-in-fortran-text-vs-other



!First reallocate to create a smaller array
allocate(output(i,entries))
output = AllData(1:i, :)

!Now save
!open(unit=10,file=BinaryData,status='replace',form='unformatted')
!write(10) output
!close(10)




!Spatial



!Time
print *, 'Writing Time'
open(unit=30,file=TimeFile,status='replace',form='formatted')
do j=1,i
write(30,*) output(j,13), output(j,1), output(j,14), output(j,4)
enddo

close(30)


!Spin
print *, 'Writing Spin'
open(unit=40,file=SpinFile,status='replace',form='formatted')
do j=1,i

    !This should be vectorised for speed
    r = output(j,2) ; theta = output(j,3) ; phi = output(j,4)
    s1 = output(j,10) ; s2 = output(j,11) ; s3 = output(j,12)


    !Proper time
    tau = output(j,13)


    !Some calculations with the spin components
    Sx = s1*sin(theta)*cos(phi) + s2*r*cos(theta)*cos(phi) - s3*r*sin(theta)*sin(phi)
    Sy = s1*sin(theta)*sin(phi) + s2*r*cos(theta)*sin(phi) + s3*r*sin(theta)*cos(phi)
    Sz = s1*cos(theta) - s2*r*sin(theta)


    thetaSL = atan2(sqrt(Sx**2 + Sy**2),Sz)
    phiSL = atan2(Sy,Sx)

    !critical phase angle
    top = cos(phiSL)*cos(thetaSL) - sin(thetaSL)
    bot = Sqrt(cos(phiSL)**2 * cos(thetaSL)**2 + sin(phiSL)**2 + sin(thetaSL)**2 - cos(phiSL)*sin(2.0_dp*thetaSL))
    
    critical_phase=acos(top/bot)

    !Relativistic aberration

!    uVector(1:4) = DerivativesStore(j,1:4)
!    xVector(1:4) = output(j,1:4)
    
    !Contravariant components of the observaion vector in the coordinate frame
!    OVector(1) = 0.0_dp
!    Ovector(2) = sin(theta)*cos(phi)*ObsX + cos(theta)*ObsZ !r_cpt
!    OVector(3) = cos(theta)*cos(phi)*ObsX/r - sin(theta)*ObsZ/r !theta cpt
!    OVector(4) = -sin(phi) / (r*sin(theta)) * ObsX


!    call transform_to_comoving_frame(OVector,uVector,xVector,OVector_comoving)


    write(40,*) tau/convert_s, phi, thetaSL, phiSL, output(j,9), critical_phase, tau/PeriodEst

enddo



!Roemer

print *, 'Writing Roemer'
open(unit=30,file=RoemerFile,status='replace',form='formatted')
do j=1,i

        tau = output(j,13)
        mm = sqrt(output(j,2)**2 + a**2)
        xC = mm * sin(output(j,3)) * cos(output(j,4))
        yC = mm * sin(output(j,3)) * sin(output(j,4))
        zC = mm * cos(output(j,3)) 

write(30,*) tau/PeriodEst, xC, yC, zC, output(j,4),tau
enddo

close(30)

print *, 'Writing Einstein'
open(unit=40,file=EinsteinFile,status='replace',form='formatted')
do j=1,i

        tau = output(j,13)
        write(40,*) tau/PeriodEst, (output(j,1) - tau)/convert_s 
enddo

close(40)




if (plot .EQ. 1) then
!Save formatted data for plotting
    open(unit=20,file=PlotData,status='replace',form='formatted')
    do j = 1,i
        if (mod(real(j), coarse) .EQ. 0.0_dp) then
        mm = sqrt(output(j,2)**2 + a**2)
        xC = mm * sin(output(j,3)) * cos(output(j,4))
        yC = mm * sin(output(j,3)) * sin(output(j,4))
        zC = mm * cos(output(j,3)) 
    



        write(20, *) output(j,1), xC, yC, zC,&!t,x,y,z
                     output(j,2), output(j,3),output(j,4),& !r,theta,phi
                     output(j,5), & !p0 i.e. time component of momentum
                     output(j,9), output(j,10),output(j,11), output(j,12), & !spin
                     output(j,13) !tau

        endif 
    

    enddo
    close(20)



endif




end subroutine rk


end module rungekutta
