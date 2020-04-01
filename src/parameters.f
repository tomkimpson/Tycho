module parameters
implicit none


!Define float precision
!integer, parameter :: dp = selected_real_kind(33,4931)
integer, parameter :: dp = selected_real_kind(15, 307)
!Define useful stuff
real(kind=dp), parameter :: PI = 4.D0*ATAN(1.D0) 


!Orbital parameters

real(kind=dp), parameter :: semi_major = 8.0_dp 
real(kind=dp), parameter :: eccentricity = 0.01_dp !Orbital eccentricity
real(kind=dp), parameter :: iota = 0.0_dp !Inclination w.r.t equatorial plane in degrees
real(kind=dp), parameter :: N_orbit = 100.50_dp !Number of orbits to integrate


!BH intrinsic parameters
real(kind=dp), parameter :: MBH = 2.200d3!BH mass in solar masses
real(kind=dp), parameter :: a= -0.60_dp !BH spin parameter

!PSR intrinsic parameters
real(kind=dp), parameter :: MPSR = 1.40_dp !pulsar mass in solar masses
real(kind=dp), parameter :: RPSR = 10.0_dp !pulsar radius in km
real(kind=dp), parameter :: stheta = PI/6.0_dp, sphi = 0.0_dp 


!Integration parameters

integer(kind=dp), parameter :: adaptive = 1 !Turn on/off adaptive steps


!IO location
character(len=200) :: IO_path = '/Users/tomkimpson/Data/Tycho/'



end module parameters
