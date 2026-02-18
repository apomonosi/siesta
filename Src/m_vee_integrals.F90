
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!> \brief General purpose of the m_vee_integrals module:
!! Compute the electron–electron interactions,
!! that are expressed as the integrals of the Coulomb kernel 
!! on the wave functions of the localized basis set (e.g.\f$d\f$ atomic states),
!! labeled by the index \f$m\f$.
!!
!! Following Ref. \cite Himmetoglu:2014:RHC,  
!! \f{eqnarray*}{
!! {\tt vee\_integral\_real(m,m^{\prime},m^{\prime\prime},m^{\prime\prime\prime})}
!!      = \langle m, m^{\prime \prime} \vert V_{\rm ee} 
!!      \vert m^{\prime}, m^{\prime \prime \prime}  \rangle 
!! = \int d^{3} r \int d^{3} r^\prime 
!!   \psi_{lm}^{\ast} ({\bf r}) 
!!   \psi_{lm^{\prime\prime}}^{\ast} ({\bf r}^\prime) 
!!   \frac{e^{2}}{\vert {\bf r} - {\bf r}^{\prime} \vert}
!!   \psi_{lm^{\prime}} ({\bf r}) 
!!   \psi_{lm^{\prime \prime \prime}} ({\bf r}^{\prime}). 
!! \f}
!! Assuming that atomic (e.g. \f$d\f$ or \f$f\f$) states are chosen as the 
!! localized basis, these integrals can be factorized in a radial and 
!! an angular contributions. 
!!
!! First, these four center integrals are computed assuming that the
!! spherical harmonics are written in complex form,
!! in the function \f$ee\_integrals\_cmplx\f$, and then transformed
!! into a basis of real spherical harmonics using the transformation
!! matrices computed in \f$complex\_to\_real\_harmonics\f$.

module m_vee_integrals

  use precision, only: dp    ! Double precision

  implicit none

  private

  public :: ee_4index_int_real

  CONTAINS 

!> \brief General purpose of the function ee_integrals_cmplx:
!! compute the four center integrals assuming that the spherical harmonics
!! are complex functions
!!
!! Here, we shall compute integrals of the type
!!
!! \f{eqnarray*}{
!! {\tt vee\_integral\_cmplx (m, m^{\prime}, m^{\prime\prime}, 
!!    m^{\prime\prime\prime}) } \equiv
!!    \langle m, m^{\prime \prime} \vert V_{\rm ee} 
!!      \vert m^{\prime}, m^{\prime \prime \prime}  \rangle 
!! = \int d^{3} r \int d^{3} r^\prime 
!!   \psi_{lm}^{\ast} ({\bf r}) 
!!   \psi_{lm^{\prime\prime}}^{\ast} ({\bf r}^\prime) 
!!   \frac{e^{2}}{\vert {\bf r} - {\bf r}^{\prime} \vert}
!!   \psi_{lm^{\prime}} ({\bf r}) 
!!   \psi_{lm^{\prime \prime \prime}} ({\bf r}^{\prime}). 
!! \f}
!! 
!! asumming the functions \f$ \psi_{lm} ({\bf r}) \f$
!! are the product of a radial function times a complex spherical harmonic.
!!
!! These integrals can be factorized in a radial and
!! an angular contributions. This factorization stems from the expansion
!! of the Coulomb kernel in spherical harmonics and yields
!! \f{eqnarray*}{
!! \langle m, m^{\prime \prime} \vert V_{\rm ee}
!!      \vert m^{\prime}, m^{\prime \prime \prime}  \rangle
!! = \sum_{k} a_{k}(m,m^{\prime},m^{\prime \prime},m^{\prime\prime\prime})F^{k},
!! \f}
!!  where \f$ 0 \le k \le 2l \f$
!! (\f$l\f$ being the angular quantum number of the localized manifold with
!! \f$-l \le m \le l\f$).
!!
!! The \f$a_{k}\f$ represents the angular factors and corresponds to products
!! of Gaunt coefficients
!!
!! \f{eqnarray*}{
!!   a_{k}(m,m^\prime,m^{\prime\prime},m^{\prime\prime\prime}) =
!!   \frac{4\pi}{2k+1} \sum_{q=-k}^{k} \langle lm \vert Y_{kq}
!!                                                \vert lm^{\prime} \rangle
!!   \langle lm^{\prime \prime} \vert Y^\ast_{kq} \vert
!!                                             lm^{\prime\prime\prime} \rangle.
!! \f}     
!! The quantities \f$ F^{k} \f$ are the Slater integrals involving
!! the radial part of the
!! atomic wave functions \f$ R_{nl} \f$ (\f$ n \f$ indicating the atomic shell
!! they belong to). They have the following expression
!! \f{eqnarray*}{
!!   F^{k} = e^{2} \int d^{3} r \int d^{3} r^\prime r^{2} r^{\prime 2}
!!           R^{2}_{nl} ({\bf r})
!!           \frac{r^{k}_{<}}{r^{k+1}_{>}} R^{2}_{nl}({\bf r}^\prime),
!! \f}
!! where \f$ r_{<} \f$  and \f$ r_{>} \f$ indicate, respectively,
!! the shorter and the
!! larger radial distances between \f$ r \f$  and \f$ r^{\prime} \f$.
!! For \f$ d \f$ electrons, only \f$ F^{0} \f$, \f$ F^{2} \f$, and \f$ F^{4} \f$
!! are needed to compute the \f$ V_{\rm ee} \f$ matrix elements
!! (for higher \f$ k \f$ values the corresponding \f$ a_{k} \f$ vanish)
!! while the \f$ f \f$ electrons also require \f$ F^{6} \f$.
!!
!! The integration over the product of three spherical harmonics 
!! in the Gaunt coefficients can be 
!! simplified as a simple expression involving only a normalization 
!! constant and two Wigner 3\f$ j \f$ symbols, that will
!! be computed in the function wigner\_3j below.
!! 
!! The following expressions are borrowed from [Eq.(1.23)] of
!!
!! https://sahussaintu.files.wordpress.com/2014/03/spherical_harmonics.pdf
!!
!! So, the two products that appear in the expression needed to 
!! compute the four center matrix elements,
!! will be referred to in the function as \f$ {\tt Gab} \f$ and
!! \f$ {\tt Gapbp} \f$ respectively, and can be expressed as
!!
!! \f{eqnarray*}{
!!   {\tt Gab } \equiv 
!!   \langle lm \vert Y_{kq} \vert lm^\prime \rangle & = 
!!        \int_{0}^{2 \pi} \int_{0}^{\pi} Y^{\ast}_{l,m} (\theta, \phi) 
!!                                     Y_{k,q} (\theta, \phi) 
!!                                      Y_{l,m^{\prime}} (\theta, \phi) 
!!                                      \sin \theta d \theta d \phi
!!   \nonumber \\
!!       & = \int_{0}^{2 \pi} \int_{0}^{\pi} (-1)^{m} Y_{l,-m} (\theta, \phi) 
!!                                      Y_{k,m-m^{\prime}} (\theta, \phi) 
!!                                      Y_{l,m^{\prime}} (\theta, \phi) 
!!                                      \sin \theta d \theta d \phi
!!   \nonumber \\
!!       & = (-1)^{m} \sqrt{\frac{(2l+1)\times(2l+1)\times(2k+1)}{4 \pi}} 
!!            \left( \begin{array}{ccc}
!!                     l & k & l \\\
!!                     0 & 0 & 0 \end{array} \right)
!!            \left( \begin{array}{ccc}
!!                     l & k & l \\\
!!                     -m & m-m^{\prime} & m^{\prime} \end{array} \right)
!! \f}
!! \f{eqnarray*}{
!!   {\tt Gapbp } \equiv 
!!   \langle lm^{\prime\prime} \vert Y^{\ast}_{kq} \vert 
!!           lm^{\prime\prime\prime} \rangle & = 
!!       \int_{0}^{2 \pi} \int_{0}^{\pi} 
!!                                 Y^{\ast}_{l,m^{\prime\prime}} (\theta, \phi) 
!!                                 Y^{\ast}_{k,q} (\theta, \phi) 
!!                                 Y_{l,m^{\prime\prime\prime}} (\theta, \phi) 
!!                                 \sin \theta d \theta d \phi
!!   \nonumber \\
!!       & = \int_{0}^{2 \pi} \int_{0}^{\pi} 
!!                                 Y^{\ast}_{l,m^{\prime\prime}} (\theta, \phi) 
!!                          (-1)^q Y_{k,-q} (\theta, \phi) 
!!                                 Y_{l,m^{\prime\prime\prime}} (\theta, \phi) 
!!                                 \sin \theta d \theta d \phi
!!   \nonumber \\
!!       & = \int_{0}^{2 \pi} \int_{0}^{\pi} 
!!              (-1)^{m^{\prime\prime}} Y_{l,-m^{\prime\prime}} (\theta, \phi) 
!!                          (-1)^q Y_{k,-q} (\theta, \phi) 
!!                                 Y_{l,m^{\prime\prime\prime}} (\theta, \phi) 
!!                                 \sin \theta d \theta d \phi
!!   \nonumber \\
!!       & = \int_{0}^{2 \pi} \int_{0}^{\pi} 
!!              (-1)^{m^{\prime\prime}} Y_{l,-m^{\prime\prime}} (\theta, \phi) 
!!          (-1)^q Y_{k,m^{\prime\prime}-m^{\prime\prime\prime}} (\theta, \phi) 
!!                                 Y_{l,m^{\prime\prime\prime}} (\theta, \phi) 
!!                                 \sin \theta d \theta d \phi
!!   \nonumber \\
!!       & = (-1)^{m^{\prime\prime} + q} 
!!           \sqrt{\frac{(2l+1)\times(2l+1)\times(2k+1)}{4 \pi}} 
!!            \left( \begin{array}{ccc}
!!                     l & k & l \\\
!!                     0 & 0 & 0 \end{array} \right)
!!            \left( \begin{array}{ccc}
!!                     l & k & l \\\
!!                     -m^{\prime\prime} & 
!!                      m^{\prime\prime}-m^{\prime\prime\prime} & 
!!                      m^{\prime\prime\prime}
!!                     \end{array} \right)
!! \f}
!! 
!! In these expressions we have made use of the fact that:
!! 
!! \f$ {\tt Gab} \ne 0 \f$ if and only if 
!! \f$ q + m^{\prime} = m \Rightarrow q = m - m^{\prime}, \f$
!!
!! and
!! 
!! \f$ {\tt Gapbp} \ne 0 \f$ if and only if 
!! \f$ -q + m^{\prime\prime\prime} = m^{\prime\prime} 
!!     \Rightarrow q = m^{\prime\prime\prime} - m^{\prime\prime}. \f$
!!
!! The calculations of the Gaunt integrals were tested against a python code,
!! available in this
!! <a href="https://docs.sympy.org/latest/modules/physics/wigner.html">link</a>.
!! An example is given by 
!!
!! \f$ >>> \f$  from sympy.physics.wigner import gaunt
!!
!! \f$ >>> \f$  gaunt(1,0,1,1,0,-1)
!! 
!! Simply take into account that in the python code the Gaunt integrals
!! were defined as 
!! 
!! \f{eqnarray*}{
!!   {\tt Gaunt } (l_{1}, l_{2}, l_{3}, m_{1}, m_{2}, m_{3}) & \equiv 
!!   \int Y_{l_{1},m_{1}}(\Omega) Y_{l_{2},m_{2}}(\Omega) 
!!        Y_{l_{3},m_{3}}(\Omega) d \Omega
!!        \nonumber \\
!!        & = \sqrt{\frac{(2l_{1} +1) (2l_{2} +1) (2l_{3} +1)}{4 \pi}}
!!        {\tt Wigner3j }   (l_{1}, l_{2}, l_{3}, 0, 0, 0) 
!!        {\tt Wigner3j }   (l_{1}, l_{2}, l_{3}, m_{1}, m_{2}, m_{3}) 
!! \f}
!! So a proper change in the sign of the \f$ {\tt Gab} \ne 0 \f$
!! might have to be consider in the comparison.

  logical function ee_integrals_cmplx( l_value, Slater_F, vee_integral_cmplx )     
    implicit none
     
    integer,  intent(in)    :: l_value     ! Angular momentum of the 
                                           !    localized basis of atomic states
                                           !    (l = 2 => d-shell)
                                           !    (l = 3 => f-shell)
    real(dp), intent(in)    :: Slater_F(0:2*l_value) 
                                           ! Value of the Slater integrals
    real(dp), intent(inout) ::             &
 &      vee_integral_cmplx(2*l_value+1,2*l_value+1,2*l_value+1,2*l_value+1)
                                           ! Value of the four center integral
                                           !    <m m'' | V_ee | m' m'''>
                                           !    assuming complex spherical 
                                           !    harmonics.
                                           !    Since -l <= m <= +l,
                                           !    we need to compute and store
                                           !    (2l+1)^4 integrals.

! 
!   Internal variables
!  

    integer  :: m       ! Magnetic quantum number of the first  orbital
    integer  :: mprime  ! Magnetic quantum number of the second orbital
    integer  :: m2prime ! Magnetic quantum number of the third  orbital
    integer  :: m3prime ! Magnetic quantum number of the fourth orbital
   
    integer  :: k       ! Counter for the loop on the a_k coefficients
                        !    in the equations above
                        !    (0 <= k <= 2l_value)
    integer  :: q       ! Counter for the loop on spherical harmonic matrix
                        !    elements
                        !    (-k <= q <= k)
    real(dp) :: pi      ! Value of pi
    real(dp) :: myint   ! Partial value of the four center integral 
    real(dp) :: Gab     ! Gaunt coefficient: < lm | Y_{kq} | lm' > 
    real(dp) :: Gapbp   ! Gaunt coefficient: < lm'' | Y*_{kq} | lm''' >

    ee_integrals_cmplx = .true.
    
!
!   Initialize pi
!
    pi=acos(-1.0_dp)    

!
!   Initialize the value of the integrals
!
    vee_integral_cmplx = 0.0_dp

!   The four loops corresponding to the indexes
    do m = -l_value, l_value
      do mprime = -l_value, l_value
        do m2prime = -l_value, l_value
          do m3prime = -l_value, l_value
            myint = 0.0_dp

!           Sum on k to compute the four center matrix elements
!           < m, m'' | V_ee | m' m''' > = \sum_{k} a_{k} (m, m', m'', m''') F^k
            do k = 0, 2*l_value, 2

!             For each k, we have to compute the corresponding Gaunt
!             coefficient
!             a_k(m, m', m'', m''') = (4pi/(2k+1)) * 
!                           sum_{q=-k}^{k}  <lm|Y_{kq}|lm'> <lm''|Y*_{kq}|lm'''>
              do q = -k, k

!               The integral <lm|Y_{kq}|lm'> does not vanish if and only if
!               (m-m') = q
!               The integral <lm''|Y*_{kq}|lm'''> does not vanish if and only if
!               (m''-m''') = -q
!               Only in these two cases we need to compute the product
                if ( ( (m - mprime) .eq. q ) .and.               &
 &                   ( (m2prime - m3prime) .eq. -q) ) then

!                 Get the Gaunt coefficients 
!                 (see the expressions in the documentation above)

!
!                 Gab = <lm|Y_{kq}|lm'>
!
                  Gab = dsqrt(( (2*l_value+1)**2*(2*k+1) )/(4.0_dp*pi)) *    &
 &                      (-1)**m * wigner_3j(l_value,k,l_value,0,0,0)    *    &
 &                      wigner_3j(l_value,k,l_value,-m,m-mprime,mprime)

!!                 For debugging
!                  write(6,'(a,4i5,f12.5)')                                  &
! &                  'ee_integrals_cmplx: m,       mprime , k, q, Gab   = ', &
! &                  m, mprime, k, q, Gab
!!                 End debugging

!
!                 Gapbp = <lm''|Y*_{kq}|lm'''>
!

                  Gapbp = dsqrt(( (2*l_value+1)**2*(2*k+1) )/(4.d0*pi)) *    &
 &                       (-1)**(m2prime+q)                              *    &
 &                 wigner_3j(l_value,k,l_value,0,0,0)                   *    &
 &                 wigner_3j(l_value,k,l_value,-m2prime,m2prime-m3prime,m3prime)

!!                 For debugging
!                  write(6,'(a,4i5,f12.5)')                                  &
! &                  'ee_integrals_cmplx: m2prime, m3prime, k, q, Gapbp = ', &
! &                  m2prime, m3prime, k, q, Gapbp
!!                 End debugging

!                 Add values to the four-center integral
                  myint = myint + (4*pi/(2*k+1))        *        &
 &                                Gab * Gapbp * Slater_F(k)
                endif
              enddo     ! Close the loop on the sum over q
            enddo       ! Close the loop on the sum over k
            vee_integral_cmplx( m+l_value+1,                 &
 &                              mprime+l_value+1,            &
 &                              m2prime+l_value+1,           &
 &                              m3prime+l_value+1 ) = myint
          enddo         ! Close the loop on m3prime
        enddo           ! Close the loop on m2prime
      enddo             ! Close the loop on mprime
    enddo               ! Close the loop on m

!!   For debugging
!    do m = 1, 2*l_value + 1
!      do mprime = 1, 2*l_value + 1
!        do m2prime = 1, 2*l_value + 1
!          do m3prime = 1, 2*l_value + 1
!            write(6,'(a,4i5,f12.5)')                                          &
! &  'ee_integrals_cmplx: m, mprime, m2prime, m3prime, vee_integral_cmplx = ', &
! &            m, mprime, m2prime, m3prime,                                    &
! &            vee_integral_cmplx(m, mprime, m2prime, m3prime) 
!          enddo 
!        enddo 
!      enddo 
!    enddo 
!!   End debugging

    ee_integrals_cmplx = .false.
    
  end function ee_integrals_cmplx

!> \brief General purpose of the function wigner_3j: It computes the
!!        Wigner_3j symbols
!!
!! The 3-j symbols are given in terms of the Clebsch–Gordan coefficients by
!! 
!! The following expressions are borrowed from Wikipedia
!!
!! https://en.wikipedia.org/wiki/3-j_symbol
!!
!! \f{eqnarray*}{
!!   {\tt wigner\_3j } & \equiv 
!!            \left( \begin{array}{ccc}
!!                     j_{1} & j_{2} & j_{3} \\\
!!                     m_{1} & m_{2} & m_{3} \end{array} \right)
!!            \equiv \frac{(-1)^{j_{1}-j_{2}-m_{3}}}{\sqrt{2j_{3}+1}}
!!            \langle j_{1} m_{1} j_{2} m_{2} \vert j_{3} (-m_{3}) \rangle,
!! \f}
!! where \f$ \langle j_{1} m_{1} j_{2} m_{2} \vert j_{3} (-m_{3}) \rangle\f$
!! is a Clebsch–Gordan coefficient

  real(dp) function wigner_3j( l1, l2, l3, m1, m2, m3 )
    integer, intent(in):: l1, l2, l3 ! Angular momentum quantum numbers 
                                     !   (upper row of 3j symbol)
    integer,intent(in):: m1, m2, m3  ! Magnetic quantum numbers 
                                     !   (lower row of 3j symbol)

!   The clebsch_gordan_coef can be defined for semiinteger values of j and m,
!   that is why we need to transfor here the integer l and m's to real numbers
    wigner_3j = clebsch_gordan_coef( real(l1,kind=dp),             &
 &                                   real(l2,kind=dp),             &
 &                                   real(l3,kind=dp),             &
 &                                   real(m1,kind=dp),             &
 &                                   real(m2,kind=dp),             &
 &                                   -1.0_dp*real(m3,kind=dp) ) *  &
 &              ((-1.0_dp)**(l1-l2-m3))/dsqrt(2.0_dp*l3+1.0_dp)

  end function wigner_3j

!> \brief General purpose of the function clebsch_gordan_coef: It computes the
!!        Clebsch–Gordan coefficient
!!
!! The explicit expression for the Clebsch–Gordan coefficients, taken from
!! the Racah formula [A. Messiah, Mécanique Quantique (Dunod, Paris, 1960) Appendix C]
!! and written in wikipedia in this 
!! <a href="https://en.wikipedia.org/wiki/Table_of_Clebsch–Gordan_coefficients">link</a> 
!! is
!! \f{eqnarray*}{
!!   {\tt clebsch\_gordan\_coef } & \equiv 
!!      \langle j_{1}, j_{2} ; m_{1}, m_{2} \vert j_{1}, j_{2}, J, M \rangle =
!!      \delta_{M, m_{1} + m_{2}} 
!!      \sqrt{ \frac{(2J+1) (J + j_{1} -j_{2}) !
!!             (J - j_{1} +j_{2})! (j_{1} +j_{2} -J)! }
!!           {( j_{1} + j_{2} +J +1)!}} 
!!      \nonumber \\
!!      & \times \sqrt{ (J+M)! (J-M)! (j_{1}-m_{1})! (j_{1} + m_{1})! 
!!                      (j_{2}-m_{2})! (j_{2}+m_{2})! }
!!      \nonumber \\
!!      & \times \left( \sum_{k} \frac{(-1)^k}{k! (j_{1}+j_{2}-J-k)! 
!!                      (j_{1}-m_{1}-k)! (j_{2}+m_{2}-k)! 
!!                      (J-j_{2} + m_{1} +k)! (J - j_{1} -m_{2} +k)!}  \right)
!! \f}
!!
!! The same expressions are available in
!!
!! G. Rudnicki-Bujnowski, Computer Physics Communications 10, 245-250 (1975)
!!
!! http://phi.phys.nagoya-u.ac.jp/NOP/Papers/Rundnicki-Bujnowski-CPC10(1975)245.pdf


  real(dp) function clebsch_gordan_coef(l1,l2,l3,m1,m2,m3)

    real(dp), intent(in):: l1, l2, l3 ! Angular momentum quantum numbers 
                                      !   (upper row of 3j symbol)
                                      !   For the sake of generality,
                                      !   the Clebsch-Gordan coefficients 
                                      !   can be defined
                                      !   for semi-integer values of j.
                                      !   That is why we define the 
                                      !   l and the m as reals
    real(dp),intent(in):: m1, m2, m3  ! Magnetic quantum numbers 
                                      !   (lower row of 3j symbol)

    real(dp) :: l3max                 ! Sum of l1 plus l2
    real(dp) :: l3min                 ! Absolute value of the difference of 
                                      !    l1 and l2
    integer  :: i                     ! Counter for the loop to compute 
                                      !    the factorial
    integer  :: r                     ! Counter for the sum on the norm
    integer  :: phase(0:40)           ! (-1)^n, where n is the index within the
                                      !    array
    real(dp) :: pi                    ! Value of pi
    real(dp) :: fact(0:40)            ! Factorial of the first fourty integers
    real(dp) :: norm                  ! Normalization factor of the 
                                      !    Clebsch-Gordan coefficient
    real(dp), parameter  :: epsilon = 1.0d-6
                                      ! Small number to make comparisons 
                                      !    between real numbers
    real(dp) :: CG_coeff              ! Clebsch-Gordan coefficient

!   Initialize the output, if input does not make sense we will return 0.0
    clebsch_gordan_coef = 0.0_dp

    pi = acos(-1.0_dp)
    
!   Get some initial values and carry sanity check
    l3max = l1 + l2
    l3min = abs(l1-l2)
    if (m1 > l1 .or. m1 < -l1) return    ! m1 has to take a value between 
                                         !    -l1 < m1 < l1
                                         ! Same for m2 and m3 below   
    if (m2 > l2 .or. m2 < -l2) return
    if (m3 > l3 .or. m3 < -l3) return
    if ( abs(m1 + m2 - m3) .ge. epsilon  ) return    
                                         ! There is a Kronecker delta 
                                         !   delta_{M,m1+m2}
                                         !   in the definition of the 
                                         !   Clebsch-Gordan coefficient
    if (l3 < l3min) return
    if (l3 > l3max) return

!   Initialize some data
!   Between them, the factorial of the first fourty integers
!   and (-1) up to the index of the array
    fact(0)= 1.0_dp
    phase(0)= 1
    do i = 1, 40
      fact(i)  =  fact(i-1) * i
      phase(i) = -phase(i-1)
    enddo

!   Calculate the Clebsch-Gordan coefficient      
    norm = 0.0_dp
    do r = 0, 20
      if ( l3 - l2 + r + m1 < 0 ) cycle
      if ( l3 - l1 + r - m2 < 0 ) cycle   
      if ( l1 + l2 - l3 - r < 0 ) cycle  
      if ( l1 - r - m1      < 0 ) cycle
      if ( l2 - r + m2      < 0 ) cycle
      norm = norm + phase(r) * 1.0_dp /                     &
 &           ( fact(r) * fact( nint(l3 - l2 + r + m1) ) *   &
 &             fact( nint(l3 - l1 + r - m2) )           *   &
 &             fact( nint(l1 + l2 - l3 - r) )           *   &
 &             fact( nint(l1 - r - m1 ) )               *   &
 &             fact( nint(l2 - r + m2 ) ) )   
    end do
    CG_coeff = norm * &
               dsqrt( fact(nint(l1+m1)) * fact(nint(l1-m1))                * &
 &                    fact(nint(l2+m2)) * fact(nint(l2-m2))                * &
 &                    fact(nint(l3+m3)) * fact(nint(l3-m3)) )              * & 
 &             dsqrt( (fact(nint(l1 + l2 - l3)) * fact(nint(l2 + l3 - l1)) * &
 &             fact(nint(l3 + l1 - l2)) ) / fact(nint(l1 + l2 + l3 + 1)) ) * &
 &             dsqrt( 2.0_dp*l3 + 1.0_dp ) 

    clebsch_gordan_coef = CG_coeff
!!   For debugging
!    write(6,'(a,f12.5)')'clebsch_gordan_coef: l1 = ', l1 
!    write(6,'(a,f12.5)')'clebsch_gordan_coef: l2 = ', l2 
!    write(6,'(a,f12.5)')'clebsch_gordan_coef: l3 = ', l3 
!    write(6,'(a,f12.5)')'clebsch_gordan_coef: m1 = ', m1 
!    write(6,'(a,f12.5)')'clebsch_gordan_coef: m2 = ', m2 
!    write(6,'(a,f12.5)')'clebsch_gordan_coef: m3 = ', m3 
!    write(6,'(a,f12.5)')'clebsch_gordan_coef: norm = ', norm
!    write(6,'(a,f12.5)')'clebsch_gordan_coef: factor1 = ',             &
! &             dsqrt( fact(nint(l1+m1)) * fact(nint(l1-m1))          * &
! &                    fact(nint(l2+m2)) * fact(nint(l2-m2))          * &
! &                    fact(nint(l3+m3)) * fact(nint(l3-m3) ))                  
!    write(6,'(a,f12.5)')'clebsch_gordan_coef: factor2 = ',             &
! &             dsqrt( (fact(nint(l1 + l2 - l3)) * fact(nint(l2 + l3 - l1)) * &
! &                     fact(nint(l3 + l1 - l2)) ) / fact(nint(l1 + l2 + l3 + 1)))  * &
! &             dsqrt( 2.0_dp*l3 + 1.0_dp ) 
!    write(6,'(a,f12.5)')'clebsch_gordan_coef: clebsch_gordan_coef =', CG_coeff
!    call die()
!!   End debugging

  end function clebsch_gordan_coef

!> \brief General purpose of the function complex_to_real_harmonics: 
!! It produces the transfomation matrices from imaginary to 
!! real spherical harmonics
!!
!! The explicit expression for the matrices can be found in this
!! <a href="https://en.wikipedia.org/wiki/Spherical_harmonics">link</a> 

  logical function complex_to_real_harmonics( l_value, transmat )

    use sys,       only : die      ! Termination routine

    implicit none
     
    integer,     intent(in)    ::  l_value    
    ! Value of the angular momentum number 
    ! This function assumes the following order 
    ! of the orbitals for the different 
    ! magnetic quantum numbers (l,m)
    ! l = 0: trivial
    ! l = 1: (-1,0,1)          => (py,pz,px)
    ! l = 2: (-2,-1,0,1,2)     => (xy,yz,z2,xz,x2-y2)
    ! l = 3: (-3-2,-1,0,1,2,3) => (y(3x2-y2), xyz, yz2, z3, z(x2-y2), x(x2-3y2))
    ! This order is in agreement with what is implemented in the
    ! subroutine rlylm in spher_harm.F90

    complex(dp), intent(out) ::    transmat( 2*l_value+1, 2*l_value+1 )
    ! Transformation matrix between the imaginary
    !    and real spherical harmonics
    
    complex(dp):: zero       ! Zero                           = (0.0, 0.0)
    complex(dp):: sq2        ! Square root of 2 (real number) = (sqrt(2.0),0.0)
    complex(dp):: one        ! One  (real number)             = (1.0,0.0)
    complex(dp):: im         ! i (imaginary number)           = (0.0,1.0)

    complex_to_real_harmonics=.true.

!
!   Define some of the numbers that will be used in the transformation matrices
!
    sq2  = cmplx( dsqrt(2.0_dp), 0.0_dp, kind=dp )
    zero = cmplx( 0.0_dp,        0.0_dp, kind=dp )
    one  = cmplx( 1.0_dp,        0.0_dp, kind=dp )
    im   = cmplx( 0.0_dp,        1.0_dp, kind=dp )

!
!   Define the transformation matrices:
!   The transformations below are from the 
!   complex to the real spherical spherical harmonics.
!   The first index of transmat refers to the real spherical harmonic
!   The second index of transmat refers to the complex spherical harmonic
!   real_spherical_harmonic(first_index) = 
!       \sum_{second_index} transmat(first_index,second_index) \times
!                           complex_spherical_harmonic(second_index)
!   These transformations are unitary.
!   Due to the way reshape is defined here, we need to take the transpose
!   of the matrix to adapt to the former convention.

    if ( l_value .eq. 0 ) then
!     if l = 0, then the transformation matrix is trivial
      transmat = 1.0_dp

    elseif ( l_value .eq. 1 ) then  
!       if l = 1, then the transformation matrix is a (3x3) matrix
!       from the three complex spherical harmonic (1,-1), (1, 0), and (1,1)
!       to the real spherical harmonic,
!       with the convention that
!       (1, -1) => p_y
!       (1,  0) => p_z
!       (1, +1) => p_x
        transmat = reshape( (/ im/sq2,  zero, im/sq2,      &
                               zero,    one,  zero,        &
                               one/sq2, zero, -one/sq2 /), &
                            (/3,3/))
        transmat = transpose(transmat)

    elseif ( l_value .eq. 2 ) then 
!       if l = 2, then the transformation matrix is a (5x5) matrix
!       from the five complex spherical harmonic 
!       (2,-2), (2, -1), (2,0), (2,1), and (2,2)
!       to the real spherical harmonic,
!       with the convention that
!       (2, -2) => d_xy
!       (2, -1) => d_yz
!       (2,  0) => d_z2
!       (2, +1) => d_xz
!       (2, +2) => d_x2-y2
        transmat = reshape( (/ im/sq2,  zero,    zero,  zero,    -im/sq2,     &
 &                             zero,    im/sq2,  zero,  im/sq2,   zero,       &
 &                             zero,    zero,    one,   zero,     zero,       &
 &                             zero,    one/sq2, zero, -one/sq2,  zero,       &
 &                             one/sq2, zero,    zero,  zero,     one/sq2 /), &
 &                          (/5,5/) )
        transmat = transpose(transmat)

    else if (l_value.eq.3) then
!       if l = 3, then the transformation matrix is a (7x7) matrix
!       from the seven complex spherical harmonic 
!       (3,-3), (3, -2), (3,-1), (3,0), (3,1), (3,2), and (3,3)
!       to the real spherical harmonic,
!       with the convention that
!       (3, -3) => f_y(3x2-y2)
!       (3, -2) => f_xyz
!       (3, -1) => f_yz2
!       (3,  0) => f_z3
!       (3, +1) => f_xz2
!       (3, +2) => f_z(x2-y2)
!       (3, +3) => f_x(x2-3y2)
        transmat = reshape( (/                                                 &
 &           im/sq2,  zero,    zero,     zero,  zero,    zero,     im/sq2,     &
 &           zero,    im/sq2,  zero,     zero,  zero,   -im/sq2,   zero,       &
 &           zero,    zero,    im/sq2,   zero,  im/sq2,  zero,     zero,       &
 &           zero,    zero,    zero,     one,   zero,    zero,     zero,       &
 &           zero,    zero,    one/sq2,  zero, -one/sq2, zero,     zero,       &
 &           zero,    one/sq2, zero,     zero,  zero,    one/sq2,  zero,       &
 &           one/sq2, zero,    zero,     zero,  zero,    zero,    -one/sq2 /), &
 &                          (/7,7/) )
        transmat = transpose(transmat)
    else 
        call die('complex_to_real_harmonics:l larger than 3 not implemented')
    end if

    complex_to_real_harmonics=.false.
    
  end function complex_to_real_harmonics

!> \brief General purpose of the subroutine ee_4index_int_real: Compute the 
!!        four-center integrals of the electron-electron interactions
!!        when the angular part of the orbitals is given as a real spherical 
!!        harmonic
!!
!! The transformation operator between real and complex spherical harmonics
!! is computed in the function complex_to_real_harmonics,
!! and the transformation that is applied below can be summarized as
!!
!! \f{eqnarray*}{
!!     \vert Y_{lm}^{\rm real}(\theta,\phi) \rangle = 
!!     \sum_{a} T_{m,a} \vert Y_{la}^{\rm complex}(\theta,\phi) \rangle
!! \f}
!!
!! \f{eqnarray*}{
!!     \langle Y_{lm}^{\rm real}(\theta,\phi) \vert = 
!!     \sum_{a} \langle Y_{la}^{\rm complex}(\theta,\phi) \vert 
!!              T^{\ast}_{a,m}
!! \f}
!!
!!
  subroutine ee_4index_int_real( l_value, Slater_F, vee_integral_real )

    implicit none
     
    integer,  intent(in)    :: l_value     ! Angular momentum of the 
                                           !    localized basis of atomic states
                                           !    (l = 2 => d-shell)
                                           !    (l = 3 => f-shell)
    real(dp), intent(inout) :: Slater_F(0:2*l_value) 
                                           ! Value of the Slater integrals
    real(dp), intent(out)   :: vee_integral_real( 2*l_value+1, 2*l_value+1, 2*l_value+1, 2*l_value+1 )

!
!   Internal variables
!

    integer  :: nfuns
                                  ! Total number of magnetic quantum numbers
    integer  :: m_a, m_b, m_ap, m_bp
                                  ! Counters for the loop on angular momenta
    integer  :: m, mprime, m2prime, m3prime
                                  ! Counters for the loop on angular momenta
    real(dp) :: imag_ints( 2*l_value+1, 2*l_value+1, 2*l_value+1, 2*l_value+1 )
                                  ! Value of the four center integral
                                  !    <m m'' | V_ee | m' m'''>
                                  !    between complex spherical harmonics
                                  !    Since -l <= m <= +l,
                                  !    we need to compute and store
                                  !    (2l+1)^4 integrals.
    complex(dp) ::  real_ints(2*l_value+1,2*l_value+1,2*l_value+1,2*l_value+1)  
                                  ! Value of the four center integral
                                  !    <m m'' | V_ee | m' m'''>
                                  !    between real spherical harmonics
                                  !    Since -l <= m <= +l,
                                  !    we need to compute and store
                                  !    (2l+1)^4 integrals.
    complex(dp) ::  transmat( 2*l_value+1, 2*l_value+1 )
                                  ! Transformation matrix between the imaginary
                                  !    and real spherical harmonics
    complex(dp) ::  trans_inv( 2*l_value+1, 2*l_value+1 )
                                  ! Inverse of the transformation matrix 
                                  !    between the imaginary
                                  !    and real spherical harmonics

!!   For debugging
!    write(6,'(a,i5)')                             & 
! &    'ee_4index_int_real: l    = ', l_value
!    do a = 0, 2*l_value
!      write(6,'(a,i1,a,f12.5,a)')                 & 
! &      'ee_4index_int_real: F(', a,') = ', Slater_F(a) * 13.6058_dp, ' eV'
!    enddo 
!!   End debugging

!   Get the four center integrals <m m'' | V_ee | m' m'''>
!   computed between complex spherical  harmonics
!   (if returned, something went wrong)
    if( ee_integrals_cmplx( l_value, Slater_F, imag_ints ) ) return 

!!   For debugging
!    do m = 1, 2*l_value + 1
!      do mprime = 1, 2*l_value + 1
!        do m2prime = 1, 2*l_value + 1
!          do m3prime = 1, 2*l_value + 1
!            write(6,'(a,4i5,f12.5)')                                  &
! &            'ee_4index_int_real: m, mp, m2p, m3p, imag_ints = ',    &
! &            m, mprime, m2prime, m3prime,                            &
! &            imag_ints(m, mprime, m2prime, m3prime)
!          enddo 
!        enddo 
!      enddo 
!    enddo 
!!   End debugging

!   Get the transformation matrix from imaginary to real harmonics  
    if( complex_to_real_harmonics( l_value, transmat ) ) return 
                                                             ! something failed
!   It is a unitary matrix, the inverse is the adjoint (conjugate transpose)
    trans_inv = transpose( conjg(transmat) )

!!   For debugging
!    do m = 1, 2*l_value +1 
!      do m_a = 1, 2*l_value +1 
!        write(6,'(a,2i5,2f12.5)')                          &
! &        'ee_4index_int_real: a, b, transmat(a,b) = ',    &
! &        m, m_a, transmat(m,m_a)
!      enddo 
!    enddo 
!
!    do m = 1, 2*l_value +1 
!      do m_a = 1, 2*l_value +1 
!        write(6,'(a,2i5,2f12.5)')                          &
! &        'ee_4index_int_real: a, b, trans_inv(a,b) = ',   &
! &        m, m_a, trans_inv(m,m_a)
!      enddo 
!    enddo 
!    call die
!!   End debugging


!   Prepare the integrals
    nfuns = 2*l_value + 1
!   Initialize the integrals between real spherical harmonics
    real_ints = 0.0_dp

!   Transform the integrals between imaginary spherical harmonics
!   into the integrals between real spherical harmonics 
    do m = 1,nfuns
      do m_a = 1, nfuns
        do mprime = 1, nfuns
          do m_b = 1, nfuns
            do m2prime = 1, nfuns
              do m_ap = 1, nfuns
                do m3prime = 1, nfuns
                  do m_bp = 1, nfuns
                    real_ints(m,mprime,m2prime,m3prime) =                      &
 &                    real_ints(m,mprime,m2prime,m3prime)    +                 &
 &                    trans_inv( m_a, m )        * transmat( mprime, m_b )   * &
 &                    trans_inv( m_ap, m2prime ) * transmat( m3prime, m_bp ) * &
 &                     imag_ints( m_a, m_b, m_ap, m_bp )
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    
    do m = 1, nfuns
      do mprime = 1, nfuns
        do m2prime = 1, nfuns
          do m3prime = 1, nfuns
            vee_integral_real(m, mprime, m2prime, m3prime) =        &
 &               real( real_ints(m, mprime, m2prime, m3prime) )
          enddo
        enddo
      enddo
    enddo

!! For debugging
!    do m = 1, nfuns
!      do mprime = 1, nfuns
!        do m2prime = 1, nfuns
!          do m3prime = 1, nfuns
!            write(6,'(a,4i5,3f12.5)')                                 &
! &            'ee_4index_int_real: m, mp, m2p, m3p, real_ints = ',    &
! &            m, mprime, m2prime, m3prime,                            &
! &            real_ints(m, mprime, m2prime, m3prime),                 &
! &            vee_integral_real(m, mprime, m2prime, m3prime) 
!          enddo
!        enddo
!      enddo
!    enddo
!    call die
!! End debugging

! Check the results are correct
!do a = 1,nfuns
!do b = 1,nfuns
!do ap = 1,nfuns
!do bp = 1,nfuns
!write(*,'("realint ",4i4," = ",2F18.10)') a-1-l_value,b-1-l_value,ap-1-l_value,bp-1-l_value,real_ints(a,b,ap,bp)
!end do
!end do
!end do
!end do
    
  end subroutine ee_4index_int_real 

end module m_vee_integrals
