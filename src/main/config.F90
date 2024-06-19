!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module dim
!
! Module to determine storage based on compile-time configuration
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!
<<<<<<< HEAD
   implicit none
#include "../../build/phantom-version.h"
   integer, parameter, public :: phantom_version_major = PHANTOM_VERSION_MAJOR
   integer, parameter, public :: phantom_version_minor = PHANTOM_VERSION_MINOR
   integer, parameter, public :: phantom_version_micro = PHANTOM_VERSION_MICRO
   character(len=*), parameter, public :: phantom_version_string = PHANTOM_VERSION_STRING
   character(len=80), parameter :: &  ! module version
      modid="$Id: 4ec8327f7d2ee81e60ce8c7b21bea0835237ae31 $"

   public

   character(len=80), parameter :: &
      tagline='Phantom v'//phantom_version_string//' (c) 2007-2023 The Authors'

   ! maximum number of particles
   integer :: maxp = 0 ! memory not allocated initially
#ifdef MAXP
   integer, parameter :: maxp_hard = MAXP
#else
   integer, parameter :: maxp_hard = 5200000
#endif

   ! maximum number of point masses
#ifdef MAXPTMASS
   integer, parameter :: maxptmass = MAXPTMASS
#else
   integer, parameter :: maxptmass = 1000
#endif
   integer, parameter :: nsinkproperties = 18

   ! storage of thermal energy or not
#ifdef ISOTHERMAL
   integer, parameter :: maxvxyzu = 3
#else
   integer, parameter :: maxvxyzu = 4
#endif

   integer :: maxTdust = 0
   logical :: store_dust_temperature = .false.
#ifdef SINK_RADIATION
   logical, parameter :: sink_radiation = .true.
#else
   logical, parameter :: sink_radiation = .false.
#endif

   ! maximum allowable number of neighbours (safest=maxp)
#ifdef MAXNEIGH
   integer, parameter :: maxneigh = MAXNEIGH
#else
   integer, parameter :: maxneigh = maxp_hard
#endif

! maxmimum storage in linklist
   integer         :: ncellsmax
   integer(kind=8) :: ncellsmaxglobal
=======
 implicit none
#include "../../build/phantom-version.h"
 integer, parameter, public :: phantom_version_major = PHANTOM_VERSION_MAJOR
 integer, parameter, public :: phantom_version_minor = PHANTOM_VERSION_MINOR
 integer, parameter, public :: phantom_version_micro = PHANTOM_VERSION_MICRO
 character(len=*), parameter, public :: phantom_version_string = PHANTOM_VERSION_STRING
 character(len=80), parameter :: &  ! module version
    modid="$Id: 4ec8327f7d2ee81e60ce8c7b21bea0835237ae31 $"

 public

 character(len=80), parameter :: &
    tagline='Phantom v'//phantom_version_string//' (c) 2007-2023 The Authors'

 ! maximum number of particles
 integer :: maxp = 0 ! memory not allocated initially
#ifdef MAXP
 integer, parameter :: maxp_hard = MAXP
#else
 integer, parameter :: maxp_hard = 5200000
#endif

 ! maximum number of point masses
#ifdef MAXPTMASS
 integer, parameter :: maxptmass = MAXPTMASS
#else
 integer, parameter :: maxptmass = 1000
#endif
 integer, parameter :: nsinkproperties = 18

 ! storage of thermal energy or not
#ifdef ISOTHERMAL
 integer, parameter :: maxvxyzu = 3
#else
 integer, parameter :: maxvxyzu = 4
#endif

 integer :: maxTdust = 0
 logical :: store_dust_temperature = .false.
#ifdef SINK_RADIATION
 logical, parameter :: sink_radiation = .true.
#else
 logical, parameter :: sink_radiation = .false.
#endif

 ! maximum allowable number of neighbours (safest=maxp)
#ifdef MAXNEIGH
 integer, parameter :: maxneigh = MAXNEIGH
#else
 integer, parameter :: maxneigh = maxp_hard
#endif

! maxmimum storage in linklist
 integer         :: ncellsmax
 integer(kind=8) :: ncellsmaxglobal
>>>>>>> bb654c34b755fa9698556d3593fbe0482b742aba

!------
! Dust
!------
<<<<<<< HEAD
   integer :: maxp_dustfrac = 0
   integer :: maxp_growth = 0
#ifdef DUST
   logical, parameter :: use_dust = .true.

#ifdef MAXDUSTLARGE
   integer, parameter :: maxdustlarge = MAXDUSTLARGE
#else
   integer, parameter :: maxdustlarge = 11
#endif
#ifdef MAXDUSTSMALL
   integer, parameter :: maxdustsmall = MAXDUSTSMALL
#else
   integer, parameter :: maxdustsmall = 11
#endif

#ifdef DUSTGROWTH
   logical, parameter :: use_dustgrowth = .true.
#else
   logical, parameter :: use_dustgrowth = .false.
#endif
#else
   logical, parameter :: use_dust = .false.
   ! integer, parameter :: ndustfluids = 0
   ! integer, parameter :: ndusttypes = 1 ! to avoid seg faults
   integer, parameter :: maxdustlarge = 1
   integer, parameter :: maxdustsmall = 1
   logical, parameter :: use_dustgrowth = .false.
#endif
   integer, parameter :: maxdusttypes = maxdustsmall + maxdustlarge

   ! kdtree
   integer, parameter :: minpart = 10

   integer :: maxprad = 0

   integer, parameter :: &
   radensumforce      = 1,&
   radenxpartvecforce = 7,&
   radensumden        = 3,&
   radenxpartvetden   = 1
#ifdef RADIATION
   logical, parameter :: do_radiation = .true.
#else
   logical, parameter :: do_radiation = .false.
#endif

   ! rhosum
   integer, parameter :: maxrhosum = 39 + &
                                    maxdustlarge - 1 + &
                                    radensumden

   ! fsum
   integer, parameter :: fsumvars = 20 ! Number of scalars in fsum
   integer, parameter :: fsumarrs = 5  ! Number of arrays  in fsum
   integer, parameter :: maxfsum  = fsumvars + &                  ! Total number of values
                                    fsumarrs*(maxdusttypes-1) + &
                                    radensumforce

! xpartveci
   integer, parameter :: maxxpartvecidens = 14 + radenxpartvetden

   integer, parameter :: maxxpartvecvars = 57 ! Number of scalars in xpartvec
   integer, parameter :: maxxpartvecarrs = 2  ! Number of arrays in xpartvec
   integer, parameter :: maxxpartvecGR   = 33 ! Number of GR values in xpartvec (1 for dens, 16 for gcov, 16 for gcon)
   integer, parameter :: maxxpartveciforce = maxxpartvecvars + &              ! Total number of values
                                             maxxpartvecarrs*(maxdusttypes-1) + &
                                             radenxpartvecforce + &
                                             maxxpartvecGR

#ifdef MPI
   logical, parameter :: mpi = .true.
#else
   logical, parameter :: mpi = .false.
#endif

   ! storage for artificial viscosity switch
   integer :: maxalpha = 0
#ifdef DISC_VISCOSITY
   integer, parameter :: nalpha = 0
#else
#ifdef CONST_AV
   integer, parameter :: nalpha = 0
#else
#ifdef USE_MORRIS_MONAGHAN
   integer, parameter :: nalpha = 1
#else
   integer, parameter :: nalpha = 3
=======
 integer :: maxp_dustfrac = 0
 integer :: maxp_growth = 0
#ifdef DUST
 logical, parameter :: use_dust = .true.

#ifdef MAXDUSTLARGE
 integer, parameter :: maxdustlarge = MAXDUSTLARGE
#else
 integer, parameter :: maxdustlarge = 11
#endif
#ifdef MAXDUSTSMALL
 integer, parameter :: maxdustsmall = MAXDUSTSMALL
#else
 integer, parameter :: maxdustsmall = 11
#endif

#ifdef DUSTGROWTH
 logical, parameter :: use_dustgrowth = .true.
#else
 logical, parameter :: use_dustgrowth = .false.
#endif
#else
 logical, parameter :: use_dust = .false.
 ! integer, parameter :: ndustfluids = 0
 ! integer, parameter :: ndusttypes = 1 ! to avoid seg faults
 integer, parameter :: maxdustlarge = 1
 integer, parameter :: maxdustsmall = 1
 logical, parameter :: use_dustgrowth = .false.
#endif
 integer, parameter :: maxdusttypes = maxdustsmall + maxdustlarge

 ! kdtree
 integer, parameter :: minpart = 10

 integer :: maxprad = 0

 integer, parameter :: &
 radensumforce      = 1,&
 radenxpartvecforce = 7,&
 radensumden        = 3,&
 radenxpartvetden   = 1
#ifdef RADIATION
 logical, parameter :: do_radiation = .true.
#else
 logical, parameter :: do_radiation = .false.
#endif

 ! rhosum
 integer, parameter :: maxrhosum = 39 + &
                                   maxdustlarge - 1 + &
                                   radensumden

 ! fsum
 integer, parameter :: fsumvars = 20 ! Number of scalars in fsum
 integer, parameter :: fsumarrs = 5  ! Number of arrays  in fsum
 integer, parameter :: maxfsum  = fsumvars + &                  ! Total number of values
                                  fsumarrs*(maxdusttypes-1) + &
                                  radensumforce

! xpartveci
 integer, parameter :: maxxpartvecidens = 14 + radenxpartvetden

 integer, parameter :: maxxpartvecvars = 57 ! Number of scalars in xpartvec
 integer, parameter :: maxxpartvecarrs = 2  ! Number of arrays in xpartvec
 integer, parameter :: maxxpartvecGR   = 33 ! Number of GR values in xpartvec (1 for dens, 16 for gcov, 16 for gcon)
 integer, parameter :: maxxpartveciforce = maxxpartvecvars + &              ! Total number of values
                                           maxxpartvecarrs*(maxdusttypes-1) + &
                                           radenxpartvecforce + &
                                           maxxpartvecGR

#ifdef MPI
 logical, parameter :: mpi = .true.
#else
 logical, parameter :: mpi = .false.
#endif

 ! storage for artificial viscosity switch
 integer :: maxalpha = 0
#ifdef DISC_VISCOSITY
 integer, parameter :: nalpha = 0
#else
#ifdef CONST_AV
 integer, parameter :: nalpha = 0
#else
#ifdef USE_MORRIS_MONAGHAN
 integer, parameter :: nalpha = 1
#else
 integer, parameter :: nalpha = 3
>>>>>>> bb654c34b755fa9698556d3593fbe0482b742aba
#endif
#endif
#endif

<<<<<<< HEAD
   ! optional storage of curl v
#ifdef CURLV
   integer, parameter :: ndivcurlv = 4
#else
   integer, parameter :: ndivcurlv = 1
#endif
   ! storage of velocity derivatives
   integer :: maxdvdx = 0  ! set to maxp when memory allocated

   ! periodic boundaries
#ifdef PERIODIC
   logical, parameter :: periodic = .true.
#else
   logical, parameter :: periodic = .false.
#endif

   ! Maximum number of particle types
   !
   integer, parameter :: maxtypes = 7 + 2*maxdustlarge - 1

   !
   ! Number of dimensions, where it is needed
   ! (Phantom is hard wired to ndim=3 in a lot of
   !  places; changing this does NOT change the
   !  code dimensionality, it just allows routines
   !  to be written in a way that are agnostic to
   !  the number of dimensions)
   !
   integer, parameter :: ndim = 3
=======
 ! optional storage of curl v
#ifdef CURLV
 integer, parameter :: ndivcurlv = 4
#else
 integer, parameter :: ndivcurlv = 1
#endif
 ! storage of velocity derivatives
 integer :: maxdvdx = 0  ! set to maxp when memory allocated

 ! periodic boundaries
#ifdef PERIODIC
 logical, parameter :: periodic = .true.
#else
 logical, parameter :: periodic = .false.
#endif

 ! Maximum number of particle types
 !
 integer, parameter :: maxtypes = 7 + 2*maxdustlarge - 1

 !
 ! Number of dimensions, where it is needed
 ! (Phantom is hard wired to ndim=3 in a lot of
 !  places; changing this does NOT change the
 !  code dimensionality, it just allows routines
 !  to be written in a way that are agnostic to
 !  the number of dimensions)
 !
 integer, parameter :: ndim = 3
>>>>>>> bb654c34b755fa9698556d3593fbe0482b742aba


!-----------------
! KROME chemistry
!-----------------
<<<<<<< HEAD
   integer :: maxp_krome = 0
#ifdef KROME
   logical, parameter :: use_krome = .true.
#else
   logical, parameter :: use_krome = .false.
=======
 integer :: maxp_krome = 0
#ifdef KROME
 logical, parameter :: use_krome = .true.
#else
 logical, parameter :: use_krome = .false.
>>>>>>> bb654c34b755fa9698556d3593fbe0482b742aba
#endif

!-----------------
! Magnetic fields
!-----------------
<<<<<<< HEAD
   integer :: maxmhd = 0
#ifdef MHD
   logical, parameter :: mhd = .true.
#else
   logical, parameter :: mhd = .false.
#endif
   integer, parameter :: maxBevol  = 4  ! size of B-arrays (Bx,By,Bz,psi)
   integer, parameter :: ndivcurlB = 4
=======
 integer :: maxmhd = 0
#ifdef MHD
 logical, parameter :: mhd = .true.
#else
 logical, parameter :: mhd = .false.
#endif
 integer, parameter :: maxBevol  = 4  ! size of B-arrays (Bx,By,Bz,psi)
 integer, parameter :: ndivcurlB = 4
>>>>>>> bb654c34b755fa9698556d3593fbe0482b742aba

! Non-ideal MHD
! if fast_divcurlB=true, then divcurlB is calculated simultaneous with density which leads to a race condition and errors (typically less than a percent)
! divcurlB is only used as diagnostics & divergence cleaning in ideal MHD, so fast_divcurlB=true is reasonable
! divcurlB is used to update the non-ideal terms, so fast_divcurlB=false is required for accuracy (especially if there will be jumps in density)
<<<<<<< HEAD
   integer :: maxmhdni = 0
#ifdef NONIDEALMHD
   logical, parameter :: mhd_nonideal    = .true.
   logical, parameter :: fast_divcurlB   = .false.
   integer, parameter :: n_nden_phantom  = 13      ! number density of chemical species, electrons & n_grains; defined in nicil == 11+2*na
#else
   logical, parameter :: mhd_nonideal    = .false.
   logical, parameter :: fast_divcurlB   = .true.
   integer, parameter :: n_nden_phantom  = 0
#endif
   logical            :: calculate_density  = .true.  ! do not toggle; initialised for efficiency
   logical            :: calculate_divcurlB = .true.  ! do not toggle; initialised for efficiency
=======
 integer :: maxmhdni = 0
#ifdef NONIDEALMHD
 logical, parameter :: mhd_nonideal    = .true.
 logical, parameter :: fast_divcurlB   = .false.
 integer, parameter :: n_nden_phantom  = 13      ! number density of chemical species, electrons & n_grains; defined in nicil == 11+2*na
#else
 logical, parameter :: mhd_nonideal    = .false.
 logical, parameter :: fast_divcurlB   = .true.
 integer, parameter :: n_nden_phantom  = 0
#endif
 logical            :: calculate_density  = .true.  ! do not toggle; initialised for efficiency
 logical            :: calculate_divcurlB = .true.  ! do not toggle; initialised for efficiency
>>>>>>> bb654c34b755fa9698556d3593fbe0482b742aba

!--------------------
! H2 Chemistry
!--------------------
<<<<<<< HEAD
   integer :: maxp_h2 = 0
#ifdef H2CHEM
   logical, parameter :: h2chemistry = .true.
#else
   logical, parameter :: h2chemistry = .false.
#endif
   integer, parameter :: nabundances = 5
=======
 integer :: maxp_h2 = 0
#ifdef H2CHEM
 logical, parameter :: h2chemistry = .true.
#else
 logical, parameter :: h2chemistry = .false.
#endif
 integer, parameter :: nabundances = 5
>>>>>>> bb654c34b755fa9698556d3593fbe0482b742aba

!--------------------
! Self-gravity
!--------------------
<<<<<<< HEAD
   integer :: maxgrav = 0
#ifdef GRAVITY
   logical, parameter :: gravity = .true.
   integer, parameter :: ngradh = 2
#else
   logical, parameter :: gravity = .false.
   integer, parameter :: ngradh = 1
=======
 integer :: maxgrav = 0
#ifdef GRAVITY
 logical, parameter :: gravity = .true.
 integer, parameter :: ngradh = 2
#else
 logical, parameter :: gravity = .false.
 integer, parameter :: ngradh = 1
>>>>>>> bb654c34b755fa9698556d3593fbe0482b742aba
#endif

!--------------------
! General relativity
!--------------------
<<<<<<< HEAD
   integer :: maxgr = 0
#ifdef GR
   logical, parameter :: gr = .true.
#else
   logical, parameter :: gr = .false.
=======
 integer :: maxgr = 0
#ifdef GR
 logical, parameter :: gr = .true.
#else
 logical, parameter :: gr = .false.
>>>>>>> bb654c34b755fa9698556d3593fbe0482b742aba
#endif

!--------------------
! Supertimestepping
!--------------------
<<<<<<< HEAD
   integer :: maxsts = 1
=======
 integer :: maxsts = 1
>>>>>>> bb654c34b755fa9698556d3593fbe0482b742aba

!--------------------
! Dust formation
!--------------------
<<<<<<< HEAD
   logical :: do_nucleation = .false.
   integer :: itau_alloc    = 0
   integer :: inucleation   = 0
   !number of elements considered in the nucleation chemical network
   integer, parameter :: nElements = 10
#ifdef DUST_NUCLEATION
   logical :: nucleation = .true.
#else
   logical :: nucleation = .false.
#endif
   integer :: maxp_nucleation = 0
=======
 logical :: do_nucleation = .false.
 integer :: itau_alloc    = 0
 integer :: inucleation   = 0
 !number of elements considered in the nucleation chemical network
 integer, parameter :: nElements = 10
#ifdef DUST_NUCLEATION
 logical :: nucleation = .true.
#else
 logical :: nucleation = .false.
#endif
 integer :: maxp_nucleation = 0
>>>>>>> bb654c34b755fa9698556d3593fbe0482b742aba

!--------------------
! MCFOST library
!--------------------
#ifdef MCFOST
<<<<<<< HEAD
   logical, parameter :: compiled_with_mcfost = .true.
#else
   logical, parameter :: compiled_with_mcfost = .false.
=======
 logical, parameter :: compiled_with_mcfost = .true.
#else
 logical, parameter :: compiled_with_mcfost = .false.
>>>>>>> bb654c34b755fa9698556d3593fbe0482b742aba
#endif

!--------------------
! Light curve stuff
!--------------------
<<<<<<< HEAD
   integer :: maxlum = 0
#ifdef LIGHTCURVE
   logical, parameter :: lightcurve = .true.
#else
   logical, parameter :: lightcurve = .false.
=======
 integer :: maxlum = 0
#ifdef LIGHTCURVE
 logical, parameter :: lightcurve = .true.
#else
 logical, parameter :: lightcurve = .false.
>>>>>>> bb654c34b755fa9698556d3593fbe0482b742aba
#endif

!--------------------
! logical for bookkeeping
!--------------------
#ifdef INJECT_PARTICLES
<<<<<<< HEAD
   logical, parameter :: particles_are_injected = .true.
#else
   logical, parameter :: particles_are_injected = .false.
=======
 logical, parameter :: particles_are_injected = .true.
#else
 logical, parameter :: particles_are_injected = .false.
>>>>>>> bb654c34b755fa9698556d3593fbe0482b742aba
#endif

!--------------------
! individual timesteps
!--------------------
#ifdef IND_TIMESTEPS
<<<<<<< HEAD
   logical, parameter :: ind_timesteps = .true.
#else
   logical, parameter :: ind_timesteps = .false.
#endif

   !--------------------
   ! Analysis array sizes
   !--------------------
   integer :: maxan = 0
   integer :: maxmhdan = 0
   integer :: maxdustan = 0
   integer :: maxgran = 0

   !--------------------
   ! Phase and gradh sizes - inconsistent with everything else, but keeping to original logic
   !--------------------
   integer :: maxphase = 0
   integer :: maxgradh = 0

   !--------------------
   ! a place to store the number of the dumpfile; required for restart dumps
   !--------------------
   integer :: idumpfile = 0
=======
 logical, parameter :: ind_timesteps = .true.
#else
 logical, parameter :: ind_timesteps = .false.
#endif

 !--------------------
 ! Analysis array sizes
 !--------------------
 integer :: maxan = 0
 integer :: maxmhdan = 0
 integer :: maxdustan = 0
 integer :: maxgran = 0

 !--------------------
 ! Phase and gradh sizes - inconsistent with everything else, but keeping to original logic
 !--------------------
 integer :: maxphase = 0
 integer :: maxgradh = 0

 !--------------------
 ! a place to store the number of the dumpfile; required for restart dumps
 !--------------------
 integer :: idumpfile = 0
>>>>>>> bb654c34b755fa9698556d3593fbe0482b742aba

contains

subroutine update_max_sizes(n,ntot)
<<<<<<< HEAD
   integer,                   intent(in) :: n
   integer(kind=8), optional, intent(in) :: ntot

   maxp = n

#ifdef KROME
   maxp_krome = maxp
#endif

#ifdef SINK_RADIATION
   store_dust_temperature = .true.
#endif

   if (store_dust_temperature) maxTdust = maxp
   if (do_nucleation) maxp_nucleation = maxp

#ifdef NCELLSMAX
   ncellsmax       = NCELLSMAX
   ncellsmaxglobal = NCELLSMAX
#else
   ncellsmax = 2*maxp
   if (present(ntot)) then
      ncellsmaxglobal = 2*ntot
   else
      ncellsmaxglobal = ncellsmax
   endif
#endif

#ifdef DUST
   maxp_dustfrac = maxp
#ifdef DUSTGROWTH
   maxp_growth = maxp
=======
 integer,                   intent(in) :: n
 integer(kind=8), optional, intent(in) :: ntot

 maxp = n

#ifdef KROME
 maxp_krome = maxp
#endif

#ifdef SINK_RADIATION
 store_dust_temperature = .true.
#endif

 if (store_dust_temperature) maxTdust = maxp
 if (do_nucleation) maxp_nucleation = maxp

#ifdef NCELLSMAX
 ncellsmax       = NCELLSMAX
 ncellsmaxglobal = NCELLSMAX
#else
 ncellsmax = 2*maxp
 if (present(ntot)) then
    ncellsmaxglobal = 2*ntot
 else
    ncellsmaxglobal = ncellsmax
 endif
#endif

#ifdef DUST
 maxp_dustfrac = maxp
#ifdef DUSTGROWTH
 maxp_growth = maxp
>>>>>>> bb654c34b755fa9698556d3593fbe0482b742aba
#endif
#endif

#ifdef DISC_VISCOSITY
<<<<<<< HEAD
   maxalpha = 0
#else
#ifdef CONST_AV
   maxalpha = 0
#else
#ifdef USE_MORRIS_MONAGHAN
   maxalpha = maxp
#else
   maxalpha = maxp
=======
 maxalpha = 0
#else
#ifdef CONST_AV
 maxalpha = 0
#else
#ifdef USE_MORRIS_MONAGHAN
 maxalpha = maxp
#else
 maxalpha = maxp
>>>>>>> bb654c34b755fa9698556d3593fbe0482b742aba
#endif
#endif
#endif

#ifdef MHD
<<<<<<< HEAD
   maxmhd = maxp
#ifdef NONIDEALMHD
   maxmhdni = maxp
=======
 maxmhd = maxp
#ifdef NONIDEALMHD
 maxmhdni = maxp
>>>>>>> bb654c34b755fa9698556d3593fbe0482b742aba
#endif
#endif

#ifdef H2CHEM
<<<<<<< HEAD
   maxp_h2 = maxp
#endif

#ifdef GRAVITY
   maxgrav = maxp
#endif

#ifdef GR
   maxgr = maxp
=======
 maxp_h2 = maxp
#endif

#ifdef GRAVITY
 maxgrav = maxp
#endif

#ifdef GR
 maxgr = maxp
>>>>>>> bb654c34b755fa9698556d3593fbe0482b742aba
#endif

#ifdef STS_TIMESTEPS
#ifdef IND_TIMESTEPS
<<<<<<< HEAD
   maxsts = maxp
=======
 maxsts = maxp
>>>>>>> bb654c34b755fa9698556d3593fbe0482b742aba
#endif
#endif

#if LIGHTCURVE
<<<<<<< HEAD
   maxlum = maxp
#endif

#ifndef ANALYSIS
   maxan = maxp
   maxmhdan = maxmhd
   maxdustan = maxp_dustfrac
   maxgran = maxgr
#endif

#ifdef RADIATION
   maxprad = maxp
   maxlum = maxp
#endif
! Very convoluted, but follows original logic...
   maxphase = maxan
   maxgradh = maxan
   maxdvdx = maxan
=======
 maxlum = maxp
#endif

#ifndef ANALYSIS
 maxan = maxp
 maxmhdan = maxmhd
 maxdustan = maxp_dustfrac
 maxgran = maxgr
#endif

#ifdef RADIATION
 maxprad = maxp
 maxlum = maxp
#endif
! Very convoluted, but follows original logic...
 maxphase = maxan
 maxgradh = maxan
 maxdvdx = maxan
>>>>>>> bb654c34b755fa9698556d3593fbe0482b742aba

end subroutine update_max_sizes

end module dim
