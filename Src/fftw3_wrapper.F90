  module fftw3_siesta_wrapper
# ifdef SIESTA__FFTW
  use, intrinsic :: iso_c_binding
  include 'fftw3.f03'
#	endif
  end module fftw3_siesta_wrapper
