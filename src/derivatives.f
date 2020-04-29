module derivatives

use parameters
use constants
use tensors

implicit none

private PK_tdot


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
real(kind=dp) :: r_i, theta_i
!Read in the data
PVector = y(5:8)
SVector = y(9:12)

r_i = y(2)
theta_i = y(3)




!Calculate the metric components
call calculate_covariant_metric(r_i, theta_i, metric)
call calculate_contravariant_metric(r_i,theta_i, metricCONTRA)



!Calculate Christoffel symbols - these are saved globally
call calculate_christoffel(r_i,theta_i) 
!Calculate Riemann tensor - components are saved globally
call calculate_riemann(r_i, theta_i)

!Calculate components of antisymmetric spin tensor - cpts saved globally
call calculate_spintensor(r_i, theta_i, &
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


!print *, Xprime
!print *, Pprime
!print *, Sprime
!stop


end subroutine derivs




subroutine RKF(yIN, yOUT,tPK)
!Integrate using the Runge-Kutta-Fehlberg algorithm

!Arguments
real(kind=dp), dimension(entries), intent(in) :: yIN
real(kind=dp), dimension(entries), intent(out) :: yOUT
real(kind=dp), intent(inout) :: tPK

!Other
real(kind=dp), dimension(size(yIN)) :: y1,y2,y3,y4,y5,y6
real(kind=dp), dimension(size(yIN)) :: k1,k2,k3,k4,k5,k6
real(kind=dp), dimension(size(yIN)) :: dy1, dy2, dy3, dy4, dy5,dy6
real(kind=dp), dimension(size(yIN)) :: ynew, yerr
real(kind=dp), dimension(size(yIN)) :: deltaErr, yscal, ratio
real(kind=dp) :: errmax

!Post-keplerian calculations
real(kind=dp) :: dt1,dt2,dt3,dt4,dt5,dt6
real(kind=dp) :: kt1,kt2,kt3,kt4,kt5,kt6

11 continue

! ---------------Y1 ------------
y1 = yIN
call derivs(y1,dy1)
k1 = h * dy1


!Y1 - additional PK
call PK_tdot(y1,dy1,dt1)
kt1 = h*dt1



! ---------------Y2 ------------
y2 = y1 + B21*k1
call derivs(y2,dy2)
k2 = h * dy2


!Y2 - additional PK
call PK_tdot(y2,dy2,dt2)
kt2 = h*dt2



! ---------------Y3 ------------
y3 = y1 + B31*k1 + B32*k2
call derivs(y3,dy3)
k3 = h * dy3


!Y3 - additional PK
call PK_tdot(y3,dy3,dt3)
kt3 = h*dt3



! ---------------Y4 ------------
y4 = y1 + B41*k1 + B42*k2 + B43*k3
call derivs(y4,dy4)
k4 = h * dy4



!Y4 - additional PK
call PK_tdot(y4,dy4,dt4)
kt4 = h*dt4


! ---------------Y5 ------------
y5 = y1 + B51*k1 + B52*k2 + B53*k3 + B54*k4 
call derivs(y5,dy5)
k5 = h * dy5



!Y5 - additional PK
call PK_tdot(y5,dy5,dt5)
kt5 = h*dt5



! ---------------Y6 ------------
y6 = y1 + B61*k1 + B62*k2 + B63*k3 + B64*k4 + B65*k5
call derivs(y6,dy6)
k6 = h * dy6



!Y6 - additional PK
call PK_tdot(y6,dy6,dt6)
kt6 = h*dt6



!Update
ynew = y1 + c1*k1  + c3*k3 + c4*k4  +c6*k6 
yerr = y1 + cbar1*k1 + cbar3*k3 + cbar4*k4 + cbar5*k5 + cbar6*k6


!And PK calculations
tPK = tPK + c1*kt1 + c3*kt3 + c4*kt4 + c6*kt6



!Errors for adaptive stepsize

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



subroutine PK_tdot(y,dy,dt)

!What is the derivative dt/dtau according to the first order PK metric

!Arguments
real(kind=dp), intent(IN), dimension(entries) :: y
real(kind=dp), intent(IN), dimension(entries) :: dy !deriatives
real(kind=dp), intent(out) :: dt


real(kind=dp) :: r,theta,phi,rh,pot
real(kind=dp) :: ur, utheta, uphi
real(kind=dp) :: vmag, vx, vy,vz,top,bot

r = y(2) ; theta = y(3) ; phi = y(4) 
ur = dy(2) ; utheta = dy(3) ; uphi = dy(4)


!Define the harmonic radial coordinate
rh = r-1.0_dp

!And the general potential
pot = -1.0_dp / rh






vx = sin(theta)*cos(phi)*ur + (r-1.0_dp)*cos(theta)*cos(phi) * utheta - (r-1.0_dp)*sin(theta)*sin(phi)*uphi

vy = sin(theta)*sin(phi)*ur + (r-1.0_dp)*cos(theta)*sin(phi) * utheta + (r-1.0_dp)*sin(theta)*cos(phi)*uphi

vz = cos(theta)*ur - (r-1.0_dp)*sin(theta)*utheta

vmag = vx**2 + vy**2 + vz**2



dt = (1.0_dp + (1.0_dp + 2.0_dp/rh)*vmag)/(1-2.0_dp/rh)
dt = sqrt(dt)

end subroutine PK_tdot




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
