!vw ------------------------------------------------------
!vw Compiler

!vw whether the compiler supports the attribute
#if defined(__HASNT_F08_CONTIGUOUS)
#define __CONTIGUOUS
#else
#define __CONTIGUOUS , CONTIGUOUS
#endif

!vw whether the compiler supports the attribute
#if defined(__HASNT_F08_ASYNCHRONOUS)
#define __ASYNCHRONOUS
#else
#define __ASYNCHRONOUS , ASYNCHRONOUS
#endif

!vw whether the compiler supports the attribute
#if defined(__HASNT_F03_ISO_FORTRAN_ENV)
#define _HASNT_F03_ISO_FORTRAN_ENV
#define _HASNT_F08_ISO_FORTRAN_ENV
#endif

!vw whether the compiler supports the attribute
#if defined(__HASNT_F08_ISO_FORTRAN_ENV)
#define _HASNT_F08_ISO_FORTRAN_ENV
#endif

!vw whether the compiler supports command line
#if defined(__HASNT_F03_EXECUTE_COMMAND_LINE)
#define _HASNT_F03_EXECUTE_COMMAND_LINE
#endif

#if defined(__HASNT_F08_POINTER_REMAPPING)
#define _HASNT_F08_POINTER_REMAPPING
#endif

!vw whether the mpi library support standard 3.0
#if defined(__HASNT_MPI_30)
#define _HASNT_MPI_30
#endif

!vw whether the compiler supports the attribute
#if defined(__HASNT_OMP_COLLAPSE)
#define __COLLAPSE2
#define __COLLAPSE3
#else
#define __COLLAPSE2 collapse(2)
#define __COLLAPSE3 collapse(3)
#endif

!vw whether the compiler supports openmp 4.5
#if defined(__HASNT_OMP_45)
#define _HASNT_OMP_45
#endif

!vw whether the compiler supports omp_set/get_nested
#if defined(__HASNT_OMP_SET_NESTED)
#define _HASNT_OMP_SET_NESTED
#endif

!vw whether the compiler supports stream io
#if defined(__HASNT_BF_STREAM_IO)
#define _HASNT_BF_STREAM_IO
#endif

!vw bypass multithreading if not supported
#if defined(__HASNT_MULTITHREAD_MPI_SUPPORT)
#define _HASNT_MULTITHREAD_MPI_SUPPORT
#endif

!vw whether to use the external c math library
#if defined(__HAS_EXTERNAL_C_ERF)
#define _HAS_EXTERNAL_C_ERF
#endif

!vw ------------------------------------------------------
!vw EXTERNAL LIBS
#if defined(__HAS_LIBXC)
#define _HAS_LIBXC
#endif

!tk DGEMMT not yet in standard BLAS
#if defined(__HAS_DGEMMT) || defined(__INTEL_MKL)
#define _HAS_DGEMMT
#endif

!tk grimme lib for vdw
#if defined(__HAS_LIBGRIMMEVDW)
#define _HAS_LIBGRIMMEVDW
#endif

!tk use scratchmodule
#if defined(__USE_SCRATCHLIBRARY)
#define _USE_SCRATCHLIBRARY
#endif

!tk Disable Scalapack
#if defined(__HAS_SCALAPACK) || defined(__INTEL_MKL)
#define _HAS_SCALAPACK
#endif

!------------------------------------------------------
! FFT library
!
! The configuration selects the FFT library with one of
!   __HAS_FFT_FFTW3    FFTW3 (or the FFTW3 interface of Intel MKL)
!   __HAS_FFT_ESSL     IBM ESSL
!   __HAS_FFT_HP       HP MLIB
!   __HAS_FFT_DEFAULT  built-in FFT (from S. Goedecker)
! and the sources only test the _HAS_FFT_* macro defined here, so
! that exactly one library is used. An explicitly requested library
! takes precedence over the built-in FFT. __INTEL_MKL selects FFTW3
! unless another FFT is requested, and without any selection the
! built-in FFT is used.
!------------------------------------------------------
#if defined(__HAS_FFT_FFTW3)
#define _HAS_FFT_FFTW3
#elif defined(__HAS_FFT_ESSL)
#define _HAS_FFT_ESSL
#elif defined(__HAS_FFT_HP)
#define _HAS_FFT_HP
#elif defined(__HAS_FFT_DEFAULT) || !defined(__INTEL_MKL)
#define _HAS_FFT_DEFAULT
#else
#define _HAS_FFT_FFTW3
#endif

! FFT tuning switches
!   __FFT_HAS_LOW_LEVEL_TIMERS    time the individual steps of the FFT
!   __FFT_HASNT_THREADED_COPIES   run the pack/unpack copies on one thread
!   __FFT_HAS_SPECIAL_COPY=<n>    variant of the pack/unpack copies
!   __FFT_NCACHE=<n>              cache size used by the built-in FFT
#if defined(__FFT_HAS_LOW_LEVEL_TIMERS)
#define _FFT_LOW_LEVEL_TIMERS .TRUE.
#else
#define _FFT_LOW_LEVEL_TIMERS .FALSE.
#endif

#if defined(__FFT_HASNT_THREADED_COPIES)
#define _FFT_NO_THREADED_COPIES .TRUE.
#else
#define _FFT_NO_THREADED_COPIES .FALSE.
#endif

#if defined(__FFT_HAS_SPECIAL_COPY)
#define _FFT_SPECIAL_COPY __FFT_HAS_SPECIAL_COPY
#else
#define _FFT_SPECIAL_COPY 0
#endif

#if defined(__FFT_NCACHE)
#define _FFT_NCACHE __FFT_NCACHE
#endif

!gm elpa lib
#if defined(__HAS_LIBELPA)
#define _HAS_LIBELPA
#endif
!tk ------------------------------------------------------
!tk use scratch module
!vw ------------------------------------------------------
!vw CUDA 

!vw whether CUDA is supported
#if defined(__HAS_CUDA)
#define _HAS_CUDA
#endif

!vw whether timer are used
#if defined(__HAS_NVTX_TIMER)
#define __NVTX_TIMER_START( name )  CALL nvtx_range_push_a ( name )
#define __NVTX_TIMER_STOP  CALL nvtx_range_pop ( )
#else
#define __NVTX_TIMER_START( name )
#define __NVTX_TIMER_STOP
#endif

!tk ------------------------------------------------------
!tk INTEL MKL needs special setting for nested calls
#if defined(__INTEL_MKL)
#define _INTEL_MKL
#endif
