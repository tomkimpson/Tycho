module IC

use parameters
use constants
use tensors



implicit none

private function_f, function_g, function_h, function_d, checks, function_f3, function_h3

public calculate_EQL, calculate_IC, calculate_EQL_circular, setup

contains 


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






!Declare savefiles
write(LamStr,'(F10.2)') lambda
write(QStr,'(F10.2)') epsQ
write(EStr,'(F10.2)') eccentricity
write(AStr,'(F10.2)') a

call get_environment_variable("QuadDir", PathOut)

Fname = 'data_eps='//trim(adjustl(QStr))//'lambda='//trim(adjustl(LamStr))


!Fname = 'data_lambda='//trim(adjustl(LamStr))//'_epsQ='//trim(adjustl(Qstr))//'_ecc='//trim(adjustl(EStr))//'_a='//trim(adjustl(AStr))
BinaryData = trim(adjustl(PathOut))//&
             trim(adjustl(Fname))//&
             '.dat'

PlotData = trim(adjustl(PathOut))//&
             trim(adjustl(Fname))//&
             '.txt'


TimeFile = trim(adjustl(PathOut))//&
           'V2/TimeEvolution/'//&
           trim(adjustl(Fname))//'.txt'



SpinFile = trim(adjustl(PathOut))//&
           'V2/SpinEvolution/'//&
           trim(adjustl(FileID))//'.txt'


RoemerFile = trim(adjustl(PathOut))//&
           'V2/Roemer/'//&
           trim(adjustl(FileID))//'.txt'


EinsteinFile = trim(adjustl(PathOut))//&
           'V2/Einstein/'//&
           trim(adjustl(FileID))//'.txt'

end subroutine setup


SUBROUTINE checks(RR,TT)
real(kind = dp) :: RR, TT

if (RR .LT. 0 .and. abs(RR) < 1d-16) then
print *, ' RR is negative but small. Correction applied :', RR, 0.00_dp
RR = 0.00_dp
else if (RR .LT. 0 .and. abs(RR) > 1d-16) then
print *, 'RR is negative and big. Something not right :', RR
STOP
endif


if (TT .LT. 0 .and. abs(TT) < 1d-16) then
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
real(kind=dp) :: f3, h3, ur,ut,uphi,utheta


r = r_init
theta = theta_init
phi = phi_init

sigma = r**2.0_dp +a**2 * cos(theta)
delta = r**2.0_dp +a**2 - 2.0_dp*r


call function_f3(r,f3)
call function_h3(r,h3)

PP = E* (r**2.0_dp + a**2) - a*L



ut = ((r**2+a**2)*PP/delta -a*(a*E-L))/r**2.0_dp - epsQ*(1.0_dp - 2.0_dp/r)**(-1.0_dp) * f3*E

ur = (PP**2.0_dp - delta*(r**2.0_dp + (L - a*E)**2.0_dp) - epsQ*r**4.0_dp * (1.0_dp -2.0_dp/r)*((f3-h3)*L**2.0_dp/r**2.0_dp + f3) )&
/r**4.0_dp

utheta = 0.0_dp


uphi = (a*PP/delta -a*E + L)/r**2.0_dp - epsQ*h3*L/r**2.0_dp


call checks(RR,TT)


if (ur .LT. 0.0_dp .and. abs(ur) .LT. 1e-33) then
ur = 0.0_dp
endif
!fixing float bug? check this





PVector(1) = m0 * ut
PVector(2) = - m0*sqrt(ur)  !made this negative
PVector(3) = 0.0_dp
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




if (iota .NE. 0.0_dp) then
print *, 'Error! Quadrupole orbital parameter mapping is only defined for equatorial plane'
stop
endif

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




call function_f3(r,f3)
call function_h3(r,h3)


h = r*(r-2.0_dp) -epsQ *(2.0_dp*f3*r - 2.0_dp*h3*r - f3*r**2 + h3*r**2)

end subroutine function_h




subroutine function_d(r,d)
real(kind=dp) :: r, d
real(kind=dp) :: delta, f3

delta = r**2.0_dp - 2.0_dp*r + a**2


call function_f3(r,f3)

d = (r**2.0_dp)*delta - epsQ*(2.0_dp*f3*r**3 - f3*r**4)

end subroutine function_d




subroutine function_f3(r,f3)
real(kind=dp) :: r,f3
real(kind=dp) :: F2,A1, F2_upper, F2_lower

F2_upper = -5.0_dp*(r-1.0_dp)*(2.0_dp + 6.0_dp*r - 3.0_dp*r**2)
F2_lower = 8.0_dp*r*(r-2.0_dp)
F2 = F2_upper/ F2_lower

A1= 15.0_dp*r*(r-2.0_dp)*log(r/(r-2.0_dp)) / 16.0_dp


f3 = F2 - A1


end subroutine function_f3





subroutine function_h3(r,h3)
real(kind=dp) :: r,h3
real(kind=dp) :: H2,A2



H2 = 5.0_dp*(2.0_dp -3.0_dp*r - 3.0_dp*r**2)/(8.0_dp*r)

A2 = -15.0_dp*(r**2 - 2.0_dp)*log(r/(r-2.0_dp))/16.0_dp


h3 = H2 - A2


end subroutine function_h3







end module IC

