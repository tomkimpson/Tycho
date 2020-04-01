module analysis

use parameters
use constants
use tensors

implicit none

private 


public transform_to_comoving_frame

contains


subroutine transform_to_comoving_frame(vector,u,x,vector_comoving)

!Arguments
real(kind=dp), intent(IN), dimension(4) :: vector,u,x !vector, 4-velocity, position
real(kind=dp), intent(OUT), dimension(4) :: vector_comoving

!Other
real(kind=dp), dimension(4) :: vector_covar
real(kind=dp), dimension(4,4) :: metric,metricCONTRA !covariant and contravariant metric components
real(kind=dp),dimension(4) :: u_covar
real(kind=dp) :: r, theta,phi,u0,u1,u2,u3,u0_covar,u1_covar,u2_covar,u3_covar
real(kind=dp) :: delta, grr,gthth, N1, N2, N3, pmag, mag1, mag2
real(kind=dp), dimension(4,4) :: transformation_matrix_inverse,transformation_matrix


!Calculate the metric components
call calculate_covariant_metric(x(2), x(3), metric)
call calculate_contravariant_metric(metric, metricCONTRA)


!Transform the contravariant components of the MSP 4-velocity
u_covar = MATMUL(metric,u)




!Extract variables from the vector format
r = x(2) ; theta = x(3); phi = x(4)
u0 = u(1) ; u1=u(2) ; u2 = u(3) ; u3 = u(4)
u0_covar = u_covar(1) ; u1_covar=u_covar(2) ; u2_covar = u_covar(3) ; u3_covar = u_covar(4)
grr = metric(2,2) ; gthth = metric(3,3)





!print *, '4 velocity magnitude =', u0*u0_covar + u1*u1_covar + u2*u2_covar + u3*u3_covar

!Construct transformation matrix
delta = r**2 + a**2 - 2.0_dp*r
N1 = sqrt(-grr * (u0_covar * u0 + u3_covar*u3)*(1.0_dp + u2_covar*u2))
N2 = sqrt(gthth*(1.0_dp + u2_covar*u2))
N3 = sqrt(-(u0_covar*u0 + u3_covar*u3)*delta*sin(theta)**2)



transformation_matrix_inverse(1,1) = -u0_covar
transformation_matrix_inverse(1,2) = -u1_covar
transformation_matrix_inverse(1,3) = -u2_covar
transformation_matrix_inverse(1,4) = -u3_covar

transformation_matrix_inverse(2,1) = u1_covar*u0_covar/N1
transformation_matrix_inverse(2,2) = -grr*(u0_covar*u0 + u3_covar*u3)/N1
transformation_matrix_inverse(2,3) = 0.0_dp
transformation_matrix_inverse(2,4) = u1_covar*u3_covar/N1

transformation_matrix_inverse(3,1) = u2_covar*u0_covar/N2
transformation_matrix_inverse(3,2) = u2_covar*u1_covar/N2
transformation_matrix_inverse(3,3) = gthth*(1.0_dp + u2_covar*u2)/N2
transformation_matrix_inverse(3,4) = u2_covar*u3_covar/N2

transformation_matrix_inverse(4,1) = -delta*sin(theta)**2*u3/N3
transformation_matrix_inverse(4,2) = 0.0_dp
transformation_matrix_inverse(4,3) = 0.0_dp
transformation_matrix_inverse(4,4) = delta*sin(theta)**2*u0/N3






transformation_matrix(1,1) = u0
transformation_matrix(1,2) = u1
transformation_matrix(1,3) = u2
transformation_matrix(1,4) = u3

transformation_matrix(2,1) = (u1_covar*u0) / N1
transformation_matrix(2,2) = (-(u0_covar * u0 + u3_covar * u3)) / N1
transformation_matrix(2,3) = (0.0_dp)/N1
transformation_matrix(2,4) = (u1_covar * u3)/N1

transformation_matrix(3,1) = (u2_covar*u0) / N2
transformation_matrix(3,2) = (u2_covar*u1) / N2
transformation_matrix(3,3) = (1.0_dp + u2_covar*u2)/N2
transformation_matrix(3,4) = (u2_covar*u3)/N2

transformation_matrix(4,1) = (u3_covar) / N3
transformation_matrix(4,2) = (0.0_dp) / N3
transformation_matrix(4,3) = (0.0_dp)/N3
transformation_matrix(4,4) = (-u0_covar)/N3


!Swithc the input vector to covariant form

vector_covar = matmul(metric,vector)

!Get the magnitude before the transform
mag1 = vector(1)*vector_covar(1) + vector(2)*vector_covar(2) + vector(3)*vector_covar(3) + vector(4)*vector_covar(4)
!print *, pmag



!Perform the transformation of the p vector
vector_comoving = matmul(transformation_matrix_inverse,vector)

!Get the magnitude after the transform - should be the same
mag2 = -vector_comoving(1)**2 + vector_comoving(2)**2 + vector_comoving(3)**2 + vector_comoving(4)**2


!print *, mag1
!print *, mag2
!stop


end subroutine transform_to_comoving_frame





end module analysis
