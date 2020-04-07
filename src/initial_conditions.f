module IC

use parameters
use constants
use tensors



implicit none

private function_f, function_g, function_h, function_d, checks

public calculate_EQL, calculate_IC, calculate_EQL_circular, setup, update_class

contains 




subroutine update_class()
!Arguments
!character(len=3) :: mode



if (mode .EQ. 'MPD') then

print *, mode, ' mode' 
a = a_fixed
p0 = 1e-3
lambda=1.0_dp
s0 = convert_spin*2.0_dp*PI*inertia/p0


else if (mode .EQ. 'KER') then
print *, mode, ' mode' 
a = a_fixed
lambda = 0.0_dp
s0 = 0.0_dp


else if (mode .EQ. 'SCH') then
a = 0.0_dp
lambda = 0.0_dp
s0=0.0_dp


endif


end subroutine update_class


subroutine setup()

character(len=10) :: LamStr
character(len=10) :: QStr, EStr, Astr


print *, 'System Info----- '



if (adaptive .EQ. 0) then
print *, 'Adaptive stepsize is off, with h = ', h
else
print *, 'Adaptive stepsize is on' 
endif


print *, 'SMA = ', semi_major
print *, 'Eccentricity = ', eccentricity
print *, 'Estimated Orbital Period = ', PeriodEst/(convert_s*3600.0_dp * 24.0_dp*365.0_dp), ' years'


savefile1 = trim(adjustl(IO_path))//'trajectory_'//mode//'.txt'
savefile2 = trim(adjustl(IO_path))//'TimeDelay_'//mode//'.txt'
savefile3 = trim(adjustl(IO_path))//'TimeDelay_PK.txt' !Only used for mode = 'SCH'

print *, 'Outfiles:'
print *, savefile1
print *, savefile2


end subroutine setup


SUBROUTINE checks(RR,TT)
real(kind = dp) :: RR, TT

if (RR .LT. 0 .and. abs(RR) < precision_limit) then
print *, ' RR is negative but small. Correction applied :', RR, 0.00_dp
RR = 0.00_dp
else if (RR .LT. 0 .and. abs(RR) > 1d-16) then
print *, 'RR is negative and big. Something not right :', RR
STOP
endif


if (TT .LT. 0 .and. abs(TT) < precision_limit) then
print *, ' TT is negative but small. Correction applied :', TT, 0.00_dp
TT = 0.00_dp
else if (TT .LT. 0 .and. abs(TT) > 1d-16) then
print *, 'TT is negative and big. Something not right :', TT
STOP
endif


END SUBROUTINE checks


subroutine calculate_IC(E,L,Q, &
                        SVector, PVector)

!Arguments
real(kind=dp), intent(IN) :: E,L,Q !Energy, AngMom, Carter
real(kind=dp), intent(INOUT), dimension(4) :: SVector !spin
real(kind=dp), intent(OUT),dimension(4) :: PVector !momentum vectors for IC

!Internals
real(kind=dp) :: r, theta, phi !initial location of particle
real(kind=dp) :: sigma, delta,PP,RR,TT !Some useful functions
real(kind=dp) :: tdot, rdot, thetadot, phidot !Kerr differential equations
real(kind=dp), dimension(4,4) :: metric
real(kind=dp) :: f3, h3, ur2,ut,uphi,utheta2

!Load the variables
r = r_init
theta = theta_init
phi = phi_init

!Define some useful stuff
sigma = r**2.0_dp +a**2 * cos(theta)
delta = r**2.0_dp +a**2 - 2.0_dp*r



!Calculate the initial 4 velocity
PP = E* (r**2.0_dp + a**2) - a*L




ut = ((r**2+a**2)*PP/delta -a*(a*E-L))/sigma
ur2 = (PP**2.0_dp - delta*(r**2.0_dp + (L - a*E)**2.0_dp) )/sigma**2.0_dp
utheta2 = (Q - cos(theta)**2.0_dp*(a**2.0_dp * (1.0_dp - E**2.0_dp)+L**2.0_dp/sin(theta)**2.0_dp) ) / sigma**2
uphi = (a*PP/delta -a*E + L/sin(theta)**2)/sigma




call checks(ur2,utheta2)




!Convert to a momentum
PVector(1) = m0 * ut
PVector(2) = - m0*sqrt(ur2)  !made this negative
PVector(3) = m0*sqrt(utheta2)
PVector(4) = m0*uphi



!Do we need to account for the oreintation of the initial momentum here?
!PVector(2) = PVector(2) * cos(PI/2.0_dp)
!PVector(4) = PVector(4) * sin(PI/2.0_dp)

call calculate_covariant_metric(r,theta,metric)




!Now calculate some extras using the covariant metric



SVector(1) = -( &
             (metric(2,1) * PVector(1) + metric(2,2)*PVector(2) + metric(2,3)*PVector(3) + metric(2,4) * PVector(4) ) * SVector(2)+&
             (metric(3,1) * PVector(1) + metric(3,2)*PVector(2) + metric(3,3)*PVector(3) + metric(3,4) * PVector(4) ) * SVector(3)+&
             (metric(4,1) * PVector(1) + metric(4,2)*PVector(2) + metric(4,3)*PVector(3) + metric(4,4) * PVector(4) ) * SVector(4) &
             ) / &
             (metric(1,1)*PVector(1) + metric(2,1) *PVector(2) + metric(3,1)*PVector(3) + metric(4,1) * PVector(4))



     !do we really need to calculate these here?

call magnitude(metric, PVector, m_sq)


call magnitude(metric, SVector, s_sq)

m_sq = -m_sq




end subroutine calculate_IC








subroutine calculate_EQL_circular(E,Q,L)
real(kind=dp) :: E, Q, L
real(kind=dp) :: N, dL


N = (1.00_dp - 3.00_dp/r_init + 2.00_dp*a*r_init**(-1.50_dp))**0.50_dp
E = (1.00_dp - 2.00_dp/r_init + a*r_init**(-1.50_dp))/N
L = r_init**0.50_dp * (1+(a/r_init)**2.00_dp - 2*a*r_init**(-1.50_dp))/N
dL = 0.00_dp
L = L +dL
Q = 0.00_dp




end subroutine calculate_EQL_circular




subroutine calculate_EQL(E, Q, L)
real(kind=dp) :: f1,g1,h1,d1 !f_functions used in defining the determinants
real(kind=dp) :: f2,g2,h2,d2 !f_functions used in defining the determinants
real(kind=dp) :: kappa, epsil, rho, eta, sigma !determinants   
real(kind=dp) :: DD ! labels prograge or retrograde orbits

real(kind=dp) :: E1, E2, E3, E !Different parts of the energy expression
real(kind=dp) :: L1, L2, L !Different parts of the momentum expression
real(kind=dp) :: Q !Carter constant


if (a .LT. 0.0_dp) then
  DD = -1.0_dp

else if (a .GT. 0.0_dp) then
  DD = 1.0_dp
else if (A .EQ. 0.0_dp) then
  DD = 1.0_dp
  print *, 'Spin parameter a = 0 (Schwarzchild). Setting DD = +1 in initial_EQL module'
  endif




call function_f(rp,f1)
call function_g(rp,g1)
call function_h(rp,h1)
call function_d(rp,d1)


call function_f(ra,f2)
call function_g(ra,g2)
call function_h(ra,h2)
call function_d(ra,d2)


kappa = d1*h2 - d2*h1
epsil = d1*g2 - d2*g1
rho = f1*h2 - f2*h1
eta = f1*g2 - f2*g1
sigma = g1*h2 - g2*h1


!Now calculate the energy

E1 = sigma*(sigma*epsil**2.0_dp + rho*epsil*kappa - eta*kappa**2.0_dp)
E2 = kappa*rho + 2.0_dp*epsil*sigma - 2.0_dp*DD*sqrt(E1)
E3 = rho**2.0_dp + 4.0_dp*eta*sigma

E = sqrt(E2/E3)




!And the angular momentum
L1 = -g1*E/h1
L2 = g1**2.0_dp * E**2.0_dp + (f1*E**2.0_dp - d1)*h1

L = L1 + DD*sqrt(L2)/h1



!And finally the Carter constant

Q = zMinus * (a**2 * (1.0_dp - E**2.0_dp) + L**2.0_dp/(1.0_dp-zMinus))

end subroutine calculate_EQL








subroutine function_f(r,f)
real(kind=dp) :: r,f ! in , out
real(kind=dp) :: delta
real(kind=dp) :: ans1, ans2, rr
integer(kind=dp) :: pow





delta = r**2 - 2.0_dp*r + a**2
f = r**4 + a**2 * (r*(r+2.0_dp))


end subroutine function_f


subroutine function_g(r,g)
real(kind=dp) :: r,g
g = 2.0_dp*a*r
end subroutine function_g


subroutine function_h(r,h)
real(kind=dp) :: r,h
real(kind=dp) :: f3,h3


h = r*(r-2.0_dp) 

end subroutine function_h




subroutine function_d(r,d)
real(kind=dp) :: r, d
real(kind=dp) :: delta, f3

delta = r**2.0_dp - 2.0_dp*r + a**2

d = (r**2.0_dp)*delta 

end subroutine function_d
























end module IC

