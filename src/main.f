program main



use parameters
use constants
implicit none



PeriodEst = 2.0_dp * PI * semi_major**(3.0_dp / 2.0_dp) 



print *, 'Estimated Orbital Period = ', PeriodEst/convert_s, ' seconds'
print *, 'Estimated Orbital Period = ', PeriodEst/(convert_s*3600.0_dp * 24.0_dp), ' days'
print *, 'Estimated Orbital Period = ', PeriodEst/(convert_s*3600.0_dp * 24.0_dp*365.0_dp), ' years'



print *, 'Convert_s =', convert_s


!stop









!Run A
h = 1000.0_dp
FileID = 'A'
p0 = 1e-3
lambda=1.0_dp
s0 = convert_spin*2.0_dp*PI*inertia/p0
call run()

print *, 'Completed OK'






end program main



subroutine run()

use parameters
use constants
use IC
use rungekutta

implicit none

real(kind=dp) :: E,L, Q !Energy, AngMom, Carter
real(kind=dp), dimension(4) :: SVector, PVector !spin and momentum vectors for IC. Contravariant
real(kind=dp), dimension(entries) :: Y_init !the array of initial conditions to pass to the rk solver


!Setup sma, periapsis etc.
call setup()

!Calculate the initial E,L,Q from the Keplerian orbital parameters
r_init = semi_major
r_init = rp

if (eccentricity .EQ. 0.0_dp) then

call calculate_EQL_circular(E,Q,L)

else

call calculate_EQL(E,Q,L)

endif


!Set the spatial components of the spin vector
SVector(2) = s0 * sin(stheta) * cos(sphi) !S1
SVector(3) = -s0 *cos(stheta)/r_init !S2 
SVector(4) = s0*sin(stheta)*sin(sphi)/(r_init*sin(theta_init)) !S3





!Calculate the remaining initial conditions

call calculate_IC(E,L,Q, SVector,PVector)



!Now do the numerical integration using a runge-kutta

Y_init(1) = t_init
Y_init(2) = r_init
Y_init(3) = theta_init
Y_init(4) = phi_init

Y_init(5:8) = PVector
Y_init(9:12) = SVector



print *, 'start RK'
call rk(Y_init)


end subroutine run


