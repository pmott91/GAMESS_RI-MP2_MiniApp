!  Copyright (C) 2020, Argonne National Laboratory. All Rights Reserved.
!  Licensed under the NCSA open source license

module hipblasf
  !-------------------------------------------------------------------------------------------------
  ! Interface to cuBLAS routines
  !-------------------------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  
  enum, bind(c) !:: hipblasStatus_t
    enumerator :: HIPBLAS_STATUS_SUCCESS = 0
    enumerator :: HIPBLAS_STATUS_NOT_INITIALIZED = 1
    enumerator :: HIPBLAS_STATUS_ALLOC_FAILED = 3
    enumerator :: HIPBLAS_STATUS_INVALID_VALUE = 7
    enumerator :: HIPBLAS_STATUS_ARCH_MISMATCH = 8
    enumerator :: HIPBLAS_STATUS_MAPPING_ERROR = 11
    enumerator :: HIPBLAS_STATUS_EXECUTION_FAILED = 13
    enumerator :: HIPBLAS_STATUS_INTERNAL_ERROR = 14
  end enum !hipblasStatus_t

  enum, bind(c) !:: hipblasFillMode_t
    enumerator :: HIPBLAS_FILL_MODE_LOWER = 0
    enumerator :: HIPBLAS_FILL_MODE_UPPER = 1
  end enum !hipblasFillMode_t

  enum, bind(c) !:: hipblasDiag type_t
    enumerator :: HIPBLAS_DIAG_NON_UNIT = 0
    enumerator :: HIPBLAS_DIAG_UNIT = 1
  end enum !hipblasDiag    type_t

  enum, bind(c) !:: hipblasSideMode_t
    enumerator :: HIPBLAS_SIDE_LEFT = 0
    enumerator :: HIPBLAS_SIDE_RIGHT = 1
  end enum !hipblasSideMode_t

  enum, bind(c) !:: hipblasOperation_t
    enumerator :: HIPBLAS_OP_N = 0
    enumerator :: HIPBLAS_OP_T = 1
    enumerator :: HIPBLAS_OP_C = 2
  end enum !hipblasOperation_t

  interface

    integer(c_int) function &
         hipblasInit() &
         bind(c, name="hipblasInit")
      use, intrinsic :: iso_c_binding
      
    end function hipblasInit

    integer(c_int) function &
         hipblasShutdown() &
         bind(c, name="hipblasShutdown")
      use, intrinsic :: iso_c_binding
      
    end function hipblasShutdown

    integer(c_int) function &
         hipblasCreate(handle) &
         bind(c, name="hipblasCreate")
      use, intrinsic :: iso_c_binding
      
      type(c_ptr) :: handle
    end function hipblasCreate

    integer(c_int) function &
         hipblasDestroy(handle) &
         bind(c, name="hipblasDestroy")
      use, intrinsic :: iso_c_binding
      
      type(c_ptr), value :: handle
    end function hipblasDestroy

    integer(c_int) function &
         hipblasGetStream(handle, stream) &
         bind(c, name="hipblasGetStream")
      use, intrinsic :: iso_c_binding
      
      type(c_ptr), value :: handle
      type(c_ptr) :: stream
    end function hipblasGetStream

    integer(c_int) function &
         hipblasSetStream(handle, stream) &
         bind(c, name="hipblasSetStream")
      use, intrinsic :: iso_c_binding
      
      type(c_ptr), value :: handle
      type(c_ptr), value :: stream
    end function hipblasSetStream

    integer(c_int) function &
         hipDeviceSynchronize() &
         bind(c, name="hipDeviceSynchronize")
      use, intrinsic :: iso_c_binding
      
    end function hipDeviceSynchronize

    integer(c_int) function &
         hipblasDgemmBatched(handle, transa, transb, m, n, k, alpha, dA,&
         ldda, dB, lddb, beta, dC, lddc, nbatch) &
         bind(c, name="hipblasDgemmBatched")
      use, intrinsic :: iso_c_binding
      
      type(c_ptr), value :: handle
      integer(c_int), value :: transa
      integer(c_int), value :: transb
      integer(c_int), value :: m
      integer(c_int), value :: n
      integer(c_int), value :: k
      real(c_double) :: alpha
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dB
      integer(c_int), value :: lddb
      real(c_double) :: beta
      type(c_ptr), value :: dC
      integer(c_int), value :: lddc
      integer(c_int), value :: nbatch
    end function hipblasDgemmBatched

    integer(c_int) function &
         hipblasDgemmStridedBatched(handle, transa, transb, m, n, k, &
         alpha, &
         dA, ldda, strideA, dB, lddb, strideB, beta, dC, lddc, &
         strideC,nbatch) &
         bind(c, name="hipblasDgemmStridedBatched")
      use, intrinsic :: iso_c_binding
      
      type(c_ptr), value :: handle
      integer(c_int), value :: transa
      integer(c_int), value :: transb
      integer(c_int), value :: m
      integer(c_int), value :: n
      integer(c_int), value :: k
      real(c_double) :: alpha
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      integer(c_int), value :: strideA
      type(c_ptr), value :: dB
      integer(c_int), value :: lddb
      integer(c_int), value :: strideB
      real(c_double) :: beta
      type(c_ptr), value :: dC
      integer(c_int), value :: strideC
      integer(c_int), value :: lddc
      integer(c_int), value :: nbatch
    end function hipblasDgemmStridedBatched

    integer(c_int) function &
         hipblasDgemm(transa, transb, m, n, k, alpha, dA, ldda, dB,&
         lddb, beta, dC, lddc) &
         bind(c, name="hipblasDgemm")
      use, intrinsic :: iso_c_binding
      
      character(c_char), value :: transa
      character(c_char), value :: transb
      integer(c_int), value :: m
      integer(c_int), value :: n
      integer(c_int), value :: k
      real(c_double), value :: alpha
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dB
      integer(c_int), value :: lddb
      real(c_double), value :: beta
      type(c_ptr), value :: dC
      integer(c_int), value :: lddc
    end function hipblasDgemm

    integer(c_int) function &
         hipblasDgemm(handle, transa, transb, m, n, k, alpha, dA,&
         ldda, dB, lddb, beta, dC, lddc) &
         bind(c, name="hipblasDgemm")
      use, intrinsic :: iso_c_binding
      
      type(c_ptr), value :: handle
      integer(c_int), value :: transa
      integer(c_int), value :: transb
      integer(c_int), value :: m
      integer(c_int), value :: n
      integer(c_int), value :: k
      real(c_double) :: alpha
      real(c_double),dimension(*) :: dA, dB, dC
      !type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      !type(c_ptr), value :: dB
      integer(c_int), value :: lddb
      real(c_double) :: beta
      !type(c_ptr), value :: dC
      integer(c_int), value :: lddc
    end function hipblasDgemm

  !  integer(c_int) function &
  !       cublasxtcreate(handle) &
  !       bind(c, name="cublasXtCreate")
  !    use, intrinsic :: iso_c_binding
  !    type(c_ptr) :: handle
  !  end function cublasxtcreate

  !  integer(c_int) function &
  !       cublasxtdestroy(handle) &
  !       bind(c, name="cublasXtDestroy")
  !    use, intrinsic :: iso_c_binding
  !    type(c_ptr), value :: handle
  !  end function cublasxtdestroy

  !  integer(c_int) function &
  !       cublasXtDeviceSelect(handle, nbDevices, deviceId) &
  !       bind(c, name="cublasXtDeviceSelect")
  !    use, intrinsic :: iso_c_binding
  !    type(c_ptr), value :: handle
  !    integer(c_int), value :: nbDevices
  !    integer(c_int),dimension(*) :: deviceId
  !  end function cublasXtDeviceSelect

  !  integer(c_int) function &
  !       cublasXtSetBlockDim(handle, blockDim) &
  !       bind(c, name="cublasXtSetBlockDim")
  !    use, intrinsic :: iso_c_binding
  !    type(c_ptr), value :: handle
  !    integer(c_int), value :: blockDim
  !  end function cublasXtSetBlockDim

  !  integer(c_int) function &
  !       cublasXtDgemm(handle, transa, transb, m, n, k, alpha, dA,&
  !       ldda, dB, lddb, beta, dC, lddc) &
  !       bind(c, name="cublasXtDgemm")
  !    use, intrinsic :: iso_c_binding
  !    type(c_ptr), value :: handle
  !    integer(c_int), value :: transa
  !    integer(c_int), value :: transb
  !    integer(c_int), value :: m
  !    integer(c_int), value :: n
  !    integer(c_int), value :: k
  !    real(c_double) :: alpha
  !    real(c_double),dimension(*) :: dA, dB, dC
  !    !type(c_ptr), value :: dA
  !    integer(c_int), value :: ldda
  !    !type(c_ptr), value :: dB
  !    integer(c_int), value :: lddb
  !    real(c_double) :: beta
  !    !type(c_ptr), value :: dC
  !    integer(c_int), value :: lddc
  !  end function cublasXtDgemm

  end interface

end module hipblasf
