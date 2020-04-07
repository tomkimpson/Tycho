module parameters
implicit none


!Define float precision
!integer, parameter :: dp = selected_real_kind(33,4931)
integer, parameter :: dp = selected_real_kind(15, 307)
!Define useful stuff
real(kind=dp), parameter :: PI = 4.D0*ATAN(1.D0) 


!Orbital parameters

real(kind=dp), parameter :: semi_major = 820.0_dp 
real(kind=dp), parameter :: eccentricity = 0.70_dp !Orbital eccentricity
real(kind=dp), parameter :: iota = 0.0_dp !Inclination w.r.t equatorial plane in degrees
real(kind=dp), parameter :: N_orbit = 5.50_dp !Number of orbits to integrate


!BH intrinsic parameters
real(kind=dp), parameter :: MBH = 4.31d6!BH mass in solar masses
!real(kind=dp), parameter :: a= +0.10_dp !BH spin parameter Now set later
real(kind=dp), parameter :: a_fixed = +0.60_dp

!PSR intrinsic parameters
real(kind=dp), parameter :: MPSR = 1.40_dp !pulsar mass in solar masses
real(kind=dp), parameter :: RPSR = 10.0_dp !pulsar radius in km
real(kind=dp), parameter :: stheta = PI/6.0_dp, sphi = 0.0_dp 


!Integration parameters

integer(kind=dp), parameter :: adaptive = 0 !Turn on/off adaptive steps


!IO location
character(len=200) :: IO_path = '/Users/tomkimpson/Data/Tycho/'


!Debugging
integer(kind=dp), parameter :: print_status = 0 !Turns on/off 1/0 print commands 



end module parameters
