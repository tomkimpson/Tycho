module derivatives

use parameters
use constants
use tensors
use quadrupole_expressions 

implicit none

private 


public RKF,derivs

contains


subroutine derivs(y, dy)

!Arguments
real(kind=dp), intent(IN), dimension(entries) :: y
real(kind=dp), intent(OUT), dimension(entries) :: dy !deriatives

!Other
real(kind=dp), dimension(4,4) :: metric,metricCONTRA !covariant and contravariant metric components
real(kind=dp), dimension(4) ::  SVector, PVector,perts
real(kind=dp), dimension(4) :: Xprime, Sprime, Pprime !derivative vectors

!Read in the data
PVector = y(5:8)
SVector = y(9:12)





!Calculate the metric components
call calculate_covariant_metric(y(2), y(3), metric)
call calculate_contravariant_metric(metric, metricCONTRA)



!Calculate Christoffel symbols - these are saved globally
call ChristoffelQuad(y(2),y(3)) 
!Calculate Riemann tensor - components are saved globally
call RiemannQuad(y(2), y(3))
!Calculate components of antisymmetric spin tensor - cpts saved globally

call calculate_spintensor(y(2), y(3), &
                            SVector, PVector, metricCONTRA)

!Calculate 4-velocity
call calculate_FourVelocity(PVector, metric,Xprime)


!Calculate 4-momentum
call calculate_FourMom(Xprime,PVector,Pprime)


!Calculate 4-spin
call calculate_FourSpin(Xprime,PVector,SVector,Sprime)



dy(1:4) = Xprime
dy(5:8) = Pprime
dy(9:12) = Sprime



end subroutine derivs




subroutine RKF(yIN, yOUT)
!Integrate using the Runge-Kutta-Fehlberg algorithm

!Arguments
real(kind=dp), dimension(entries), intent(in) :: yIN
real(kind=dp), dimension(entries), intent(out) :: yOUT

!Other
real(kind=dp), dimension(size(yIN)) :: y1,y2,y3,y4,y5,y6
real(kind=dp), dimension(size(yIN)) :: k1,k2,k3,k4,k5,k6
real(kind=dp), dimension(size(yIN)) :: dy1, dy2, dy3, dy4, dy5,dy6
real(kind=dp), dimension(size(yIN)) :: ynew, yerr
real(kind=dp), dimension(size(yIN)) :: deltaErr, yscal, ratio
real(kind=dp) :: errmax



11 continue

! Y1
y1 = yIN
call derivs(y1,dy1)
k1 = h * dy1



!Y2
y2 = y1 + B21*k1
call derivs(y2,dy2)
k2 = h * dy2



!Y3
y3 = y1 + B31*k1 + B32*k2
call derivs(y3,dy3)
k3 = h * dy3


!Y4
y4 = y1 + B41*k1 + B42*k2 + B43*k3
call derivs(y4,dy4)
k4 = h * dy4


!Y5
y5 = y1 + B51*k1 + B52*k2 + B53*k3 + B54*k4 
call derivs(y5,dy5)
k5 = h * dy5


!Y6
y6 = y1 + B61*k1 + B62*k2 + B63*k3 + B64*k4 + B65*k5
call derivs(y6,dy6)
k6 = h * dy6


!Update
ynew = y1 + c1*k1  + c3*k3 + c4*k4  +c6*k6 
yerr = y1 + cbar1*k1 + cbar3*k3 + cbar4*k4 + cbar5*k5 + cbar6*k6


deltaErr = abs(ynew - yerr)
yscal = abs(y1) + abs(k1) + 1.0d-3
ratio = deltaErr/yscal
errmax = escal * maxval(ratio)





if (adaptive .EQ. 1) then


if (errmax .GT. 1.0_dp) then
!This is not good. Do not update yOUT and reduce the stepsize
call ShrinkStepsize(errmax)
yOUT = yIN
goto 11
else
!This is good. Update yOUT and try to increase the stepsize a little bit
call GrowStepsize(errmax)
yOUT = ynew
endif



else


!Just always update steps

yOUT = ynew

endif









end subroutine RKF



subroutine GrowStepsize(errmax)
real(kind=dp) :: errmax


if (errmax .GT. errcon) then
h = S*h*errmax**Pgrow
else
h = h * 5.0_dp
endif


end subroutine GrowStepsize



subroutine ShrinkStepsize(errmax)
real(kind=dp) :: errmax
real(kind=dp) :: htemp

htemp = S*h*errmax**Pshrink
h = sign(max(abs(htemp),0.10_dp*abs(h)),h)



end subroutine ShrinkStepsize




end module derivatives
