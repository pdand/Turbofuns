module MODNormal

    public:: bvn, PPND16, d2norm, phi

    private :: normp

   contains

    !------------------------------------------------------------------------
! All the code (except the subroutine normp) in the module MODNormal is software obtained from
! the website http://www.math.wsu.edu/faculty/genz/software/software.html.
! Below please find the license information about Alan Genz's software.
! In addition, these Fortran code are included in the R package mnormt. Its license is GPL-2 | GPL-3, which
! is given at https://www.gnu.org/licenses/gpl-3.0.en.html

! All Alan Millers' code has been released to the public domain. More details can be
! found at https://jblevins.org/mirror/amiller/


!    The software is based on work described in the paper
!     "Numerical Computation of Rectangular Bivariate and Trivariate
!      Normal and t Probabilities", by the code author:
!
!       Alan Genz
!       Department of Mathematics
!       Washington State University
!       Pullman, WA 99164-3113
!       Email : alangenz@wsu.edu
!
!
! Copyright (C) 2013, Alan Genz,  All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided the following conditions are met:
!   1. Redistributions of source code must retain the above copyright
!      notice, this list of conditions and the following disclaimer.
!   2. Redistributions in binary form must reproduce the above copyright
!      notice, this list of conditions and the following disclaimer in
!      the documentation and/or other materials provided with the
!      distribution.
!   3. The contributor name(s) may not be used to endorse or promote
!      products derived from this software without specific prior
!      written permission.
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
! OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
! ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
! TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

FUNCTION bvn ( lower, upper, infin, correl ) RESULT(fn_val)

!     A function for computing bivariate normal probabilities.
!     Extracted from Alan Genz's package for multivariate normal integration.

!  Parameters

!     LOWER  REAL, array of lower integration limits.
!     UPPER  REAL, array of upper integration limits.
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
!            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
!            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
!     CORREL REAL, correlation coefficient.

INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN) :: lower(:), upper(:), correl
INTEGER, INTENT(IN)   :: infin(:)
REAL (dp)             :: fn_val

REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp

fn_val = 0.0_dp ! added to please a gfortran compiler

IF ( infin(1) == 2  .AND. infin(2) == 2 ) THEN
  fn_val =  bvnu ( lower(1), lower(2), correl )    &
            - bvnu ( upper(1), lower(2), correl )  &
            - bvnu ( lower(1), upper(2), correl )  &
            + bvnu ( upper(1), upper(2), correl )
ELSE IF ( infin(1) == 2  .AND. infin(2) == 1 ) THEN
  fn_val =  bvnu ( lower(1), lower(2), correl )  &
            - bvnu ( upper(1), lower(2), correl )
ELSE IF ( infin(1) == 1  .AND. infin(2) == 2 ) THEN
  fn_val =  bvnu ( lower(1), lower(2), correl )  &
            - bvnu ( lower(1), upper(2), correl )
ELSE IF ( infin(1) == 2  .AND. infin(2) == 0 ) THEN
  fn_val =  bvnu ( -upper(1), -upper(2), correl )  &
            - bvnu ( -lower(1), -upper(2), correl )
ELSE IF ( infin(1) == 0  .AND. infin(2) == 2 ) THEN
  fn_val =  bvnu ( -upper(1), -upper(2), correl )  &
            - bvnu ( -upper(1), -lower(2), correl )
ELSE IF ( infin(1) == 1  .AND. infin(2) == 0 ) THEN
  fn_val =  bvnu ( lower(1), -upper(2), -correl )
ELSE IF ( infin(1) == 0  .AND. infin(2) == 1 ) THEN
  fn_val =  bvnu ( -upper(1), lower(2), -correl )
ELSE IF ( infin(1) == 1  .AND. infin(2) == 1 ) THEN
  fn_val =  bvnu ( lower(1), lower(2), correl )
ELSE IF ( infin(1) == 0  .AND. infin(2) == 0 ) THEN
  fn_val =  bvnu ( -upper(1), -upper(2), correl )
END IF

RETURN


 CONTAINS


FUNCTION bvnu( sh, sk, r ) RESULT(fn_val)

!     A function for computing bivariate normal probabilities.

!       Yihong Ge
!       Department of Computer Science and Electrical Engineering
!       Washington State University
!       Pullman, WA 99164-2752
!       Email : yge@eecs.wsu.edu
!     and
!       Alan Genz
!       Department of Mathematics
!       Washington State University
!       Pullman, WA 99164-3113
!       Email : alangenz@wsu.edu

! BVN - calculate the probability that X is larger than SH and Y is
!       larger than SK.

! Parameters

!   SH  REAL, integration limit
!   SK  REAL, integration limit
!   R   REAL, correlation coefficient
!   LG  INTEGER, number of Gauss Rule Points and Weights

REAL (dp), INTENT(IN) :: sh, sk, r
REAL (dp)             :: fn_val

! Local variables
INTEGER              :: i, lg, ng
REAL (dp), PARAMETER :: twopi = 6.283185307179586
REAL (dp)            :: as, a, b, c, d, rs, xs
REAL (dp)            :: bvn, sn, asr, h, k, bs, hs, hk
!     Gauss Legendre Points and Weights, N =  6
! DATA ( w(i,1), x(i,1), i = 1,3) /  &
! 0.1713244923791705D+00,-0.9324695142031522D+00,  &
! 0.3607615730481384D+00,-0.6612093864662647D+00,  &
! 0.4679139345726904D+00,-0.2386191860831970D+00/
!     Gauss Legendre Points and Weights, N = 12
! DATA ( w(i,2), x(i,2), i = 1,6) /  &
! 0.4717533638651177D-01,-0.9815606342467191D+00,  &
! 0.1069393259953183D+00,-0.9041172563704750D+00,  &
! 0.1600783285433464D+00,-0.7699026741943050D+00,  &
! 0.2031674267230659D+00,-0.5873179542866171D+00,  &
! 0.2334925365383547D+00,-0.3678314989981802D+00,  &
! 0.2491470458134029D+00,-0.1252334085114692D+00/
!     Gauss Legendre Points and Weights, N = 20
! DATA ( w(i,3), x(i,3), i = 1,10) /  &
! 0.1761400713915212D-01,-0.9931285991850949D+00,  &
! 0.4060142980038694D-01,-0.9639719272779138D+00,  &
! 0.6267204833410906D-01,-0.9122344282513259D+00,  &
! 0.8327674157670475D-01,-0.8391169718222188D+00,  &
! 0.1019301198172404D+00,-0.7463319064601508D+00,  &
! 0.1181945319615184D+00,-0.6360536807265150D+00,  &
! 0.1316886384491766D+00,-0.5108670019508271D+00,  &
! 0.1420961093183821D+00,-0.3737060887154196D+00,  &
! 0.1491729864726037D+00,-0.2277858511416451D+00,  &
! 0.1527533871307259D+00,-0.7652652113349733D-01/
REAL (dp), PARAMETER :: w(10,3) = RESHAPE( (/      &
      0.1713244923791705D+00, 0.3607615730481384D+00, 0.4679139345726904D+00, &
        0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,                  &
      0.4717533638651177D-01, 0.1069393259953183D+00, 0.1600783285433464D+00, &
      0.2031674267230659D+00, 0.2334925365383547D+00, 0.2491470458134029D+00, &
       0.0D0, 0.0D0, 0.0D0, 0.0D0,                        &
      0.1761400713915212D-01, 0.4060142980038694D-01, 0.6267204833410906D-01, &
      0.8327674157670475D-01, 0.1019301198172404D+00, 0.1181945319615184D+00, &
      0.1316886384491766D+00, 0.1420961093183821D+00, 0.1491729864726037D+00, &
      0.1527533871307259D+00 /), (/ 10, 3 /) )
REAL (dp), PARAMETER :: x(10,3) = RESHAPE( (/      &
      -0.9324695142031522D+00, -0.6612093864662647D+00,   &
      -0.2386191860831970D+00,  0.0D0, 0.0D0, 0.0D0,      &
       0.0D0, 0.0D0, 0.0D0, 0.0D0,                        &
      -0.9815606342467191D+00, -0.9041172563704750D+00,   &
      -0.7699026741943050D+00, -0.5873179542866171D+00,   &
      -0.3678314989981802D+00, -0.1252334085114692D+00,   &
       0.0D0, 0.0D0, 0.0D0, 0.0D0,                        &
      -0.9931285991850949D+00, -0.9639719272779138D+00,   &
      -0.9122344282513259D+00, -0.8391169718222188D+00,   &
      -0.7463319064601508D+00, -0.6360536807265150D+00,   &
      -0.5108670019508271D+00, -0.3737060887154196D+00,   &
      -0.2277858511416451D+00, -0.7652652113349733D-01 /), (/ 10, 3 /) )

IF ( ABS(r) < 0.3 ) THEN
  ng = 1
  lg = 3
ELSE IF ( ABS(r) < 0.75 ) THEN
  ng = 2
  lg = 6
ELSE
  ng = 3
  lg = 10
END IF
h = sh
k = sk
hk = h*k
bvn = zero
IF ( ABS(r) < 0.925 ) THEN
  hs = ( h*h + k*k )/2
  asr = ASIN(r)
  DO  i = 1, lg
    sn = SIN(asr*( x(i,ng)+1 )/2)
    bvn = bvn + w(i,ng)*EXP( ( sn*hk - hs )/(one - sn*sn ) )
    sn = SIN(asr*(-x(i,ng)+1 )/2)
    bvn = bvn + w(i,ng)*EXP( ( sn*hk - hs )/(one - sn*sn ) )
  END DO
  bvn = bvn*asr/(2*twopi) + phi(-h)*phi(-k)
ELSE
  IF ( r < zero ) THEN
    k = -k
    hk = -hk
  END IF
  IF ( ABS(r) < one ) THEN
    as = ( one - r )*( one + r )
    a = SQRT(as)
    bs = ( h - k )**2
    c = ( 4. - hk )/8
    d = ( 12. - hk )/16.
    bvn = a*EXP( -(bs/as + hk)/2. )  &
    *( one - c*(bs - as)*(one - d*bs/5.)/3. + c*d*as*as/5. )
    IF ( hk > -160. ) THEN
      b = SQRT(bs)
      bvn = bvn - EXP(-hk/2)*SQRT(twopi)*phi(-b/a)*b  &
                      *( one - c*bs*( one - d*bs/5. )/3. )
    END IF
    a = a/2
    DO i = 1, lg
      xs = ( a*(x(i,ng) + one) )**2
      rs = SQRT( one - xs )
      bvn = bvn + a*w(i,ng)*( EXP( -bs/(2*xs) - hk/(1+rs) )/rs  &
                - EXP( -(bs/xs+hk)/2. )*( one + c*xs*( one + d*xs ) ) )
      xs = as*(-x(i,ng) + one)**2/4.
      rs = SQRT( 1 - xs )
      bvn = bvn + a*w(i,ng)*EXP( -(bs/xs + hk)/2 )  &
                * ( EXP( -hk*(one - rs)/(2*(one + rs)) )/rs - &
                       ( one + c*xs*( one + d*xs ) ) )
    END DO
    bvn = -bvn/twopi
  END IF
  IF ( r > 0 ) bvn =  bvn + phi( -MAX( h, k ) )
  IF ( r < 0 ) bvn = -bvn + MAX( zero, phi(-h) - phi(-k) )
END IF
fn_val = bvn

RETURN
END FUNCTION bvnu

END FUNCTION bvn


SUBROUTINE normp(z, p, q, pdf)

! Normal distribution probabilities accurate to 1.e-15.
! Z = no. of standard deviations from the mean.
! P, Q = probabilities to the left & right of Z.   P + Q = 1.
! PDF = the probability density.

! Based upon algorithm 5666 for the error function, from:
! Hart, J.F. et al, 'Computer Approximations', Wiley 1968

! Programmer: Alan Miller

! Latest revision of Fortran 77 version - 30 March 1986
! Latest revision of Fortran 90 version - 12 August 1997

INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)


REAL(dp), INTENT(IN)            :: z
REAL(dp), INTENT(OUT), OPTIONAL :: p, q, pdf

! Local variables
REAL(dp), PARAMETER :: p0 = 220.2068679123761D0, p1 = 221.2135961699311D0,  &
                        p2 = 112.0792914978709D0, p3 = 33.91286607838300D0,  &
                        p4 = 6.373962203531650D0, p5 = .7003830644436881D0,  &
                        p6 = .3526249659989109D-01,  &
                        q0 = 440.4137358247522D0, q1 = 793.8265125199484D0,  &
                        q2 = 637.3336333788311D0, q3 = 296.5642487796737D0,  &
                        q4 = 86.78073220294608D0, q5 = 16.06417757920695D0,  &
                        q6 = 1.755667163182642D0, q7 = .8838834764831844D-1, &
                        cutoff = 7.071D0, root2pi = 2.506628274631001D0
REAL(dp), PARAMETER  :: zero = 0.0D0, one = 1.0d0
REAL(dp)            :: zabs, expntl, pp, qq, ppdf

zabs = ABS(z)

! |Z| > 37.

IF (zabs > 37.d0) THEN
  IF (PRESENT(pdf)) pdf = zero
  IF (z > zero) THEN
    IF (PRESENT(p)) p = one
    IF (PRESENT(q)) q = zero
  ELSE
    IF (PRESENT(p)) p = zero
    IF (PRESENT(q)) q = one
  END IF
  RETURN
END IF

! |Z| <= 37.

expntl = EXP(-0.5D0*zabs**2)
ppdf = expntl/root2pi
IF (PRESENT(pdf)) pdf = ppdf

! |Z| < CUTOFF = 10/sqrt(2).

IF (zabs < cutoff) THEN
  pp = expntl*((((((p6*zabs + p5)*zabs + p4)*zabs + p3)*zabs + p2)*zabs     &
                   + p1)*zabs + p0) / (((((((q7*zabs + q6)*zabs + q5)*zabs &
                   + q4)*zabs + q3)*zabs + q2)*zabs + q1)*zabs +q0)

! |Z| >= CUTOFF.

ELSE
  pp = ppdf/(zabs + one/(zabs + 2.d0/(zabs + 3.d0/(zabs + 4.d0/(zabs + 0.65D0)))))
END IF

IF (z < zero) THEN
  qq = one - pp
ELSE
  qq = pp
  pp = one - qq
END IF

IF (PRESENT(p)) p = pp
IF (PRESENT(q)) q = qq

RETURN
END SUBROUTINE normp



FUNCTION phi(z) RESULT(p)

INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL(dp), INTENT(IN) :: z
REAL(dp)             :: p

CALL normp(z, p)

RETURN
END FUNCTION phi


!------------------------------------------------------------------------


DOUBLE PRECISION FUNCTION PPND16 (P, IFAULT)
! ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

! Produces the normal deviate Z corresponding to a given lower
! tail area of P; Z is accurate to about 1 part in 10**16.

! The hash sums below are the sums of the mantissas of the
! coefficients.   They are included for use in checking
! transcription.

DOUBLE PRECISION P
INTEGER IFAULT

DOUBLE PRECISION ZERO, ONE, HALF, SPLIT1, SPLIT2, CONST1,   &
 CONST2, A0, A1, A2, A3, A4, A5, A6, A7, B1, B2, B3, &
        B4, B5, B6, B7,                                     &
    C0, C1, C2, C3, C4, C5, C6, C7, D1, D2, D3, D4, D5, &
    D6, D7, E0, E1, E2, E3, E4, E5, E6, E7, F1, F2, F3, &
    F4, F5, F6, F7, Q, R

PARAMETER (ZERO = 0.D0, ONE = 1.D0, HALF = 0.5D0, &
SPLIT1 = 0.425D0, SPLIT2 = 5.D0,          &
CONST1 = 0.180625D0, CONST2 = 1.6D0)

! Coefficients for P close to 0.5
PARAMETER (A0 = 3.3871328727963666080D0,  &
    A1 = 1.3314166789178437745D+2, &
    A2 = 1.9715909503065514427D+3, &
    A3 = 1.3731693765509461125D+4, &
    A4 = 4.5921953931549871457D+4, &
    A5 = 6.7265770927008700853D+4, &
    A6 = 3.3430575583588128105D+4, &
    A7 = 2.5090809287301226727D+3, &
    B1 = 4.2313330701600911252D+1, &
    B2 = 6.8718700749205790830D+2, &
    B3 = 5.3941960214247511077D+3, &
    B4 = 2.1213794301586595867D+4, &
    B5 = 3.9307895800092710610D+4, &
    B6 = 2.8729085735721942674D+4, &
    B7 = 5.2264952788528545610D+3)
! HASH SUM AB  55.88319 28806 14901 4439

! Coefficients for P not close to 0, 0.5 or 1.
PARAMETER (C0 = 1.42343711074968357734D0,   &
    C1 = 4.63033784615654529590D0,   &
    C2 = 5.76949722146069140550D0,   &
    C3 = 3.64784832476320460504D0,   &
    C4 = 1.27045825245236838258D0,   &
    C5 = 2.41780725177450611770D-1,  &
    C6 = 2.27238449892691845833D-2,  &
    C7 = 7.74545014278341407640D-4,  &
    D1 = 2.05319162663775882187D0,   &
    D2 = 1.67638483018380384940D0,   &
    D3 = 6.89767334985100004550D-1,  &
    D4 = 1.48103976427480074590D-1,  &
    D5 = 1.51986665636164571966D-2,  &
    D6 = 5.47593808499534494600D-4,  &
    D7 = 1.05075007164441684324D-9)
! HASH SUM CD  49.33206 50330 16102 89036

! Coefficients for P near 0 or 1.
PARAMETER (E0 = 6.65790464350110377720D0,   &
    E1 = 5.46378491116411436990D0,   &
    E2 = 1.78482653991729133580D0,   &
    E3 = 2.96560571828504891230D-1,  &
    E4 = 2.65321895265761230930D-2,  &
    E5 = 1.24266094738807843860D-3,  &
    E6 = 2.71155556874348757815D-5,  &
    E7 = 2.01033439929228813265D-7,  &
    F1 = 5.99832206555887937690D-1,  &
    F2 = 1.36929880922735805310D-1,  &
    F3 = 1.48753612908506148525D-2,  &
    F4 = 7.86869131145613259100D-4,  &
    F5 = 1.84631831751005468180D-5,  &
    F6 = 1.42151175831644588870D-7,  &
    F7 = 2.04426310338993978564D-15)
! HASH SUM EF    47.52583 31754 92896 71629

IFAULT = 0
Q = P - HALF
IF (ABS(Q) .LE. SPLIT1) THEN
  R = CONST1 - Q * Q
  PPND16 = Q * (((((((A7 * R + A6) * R + A5) * R + A4) * R + A3) &
                    * R + A2) * R + A1) * R + A0) / &
               (((((((B7 * R + B6) * R + B5) * R + B4) * R + B3) &
            * R + B2) * R + B1) * R + ONE)
  RETURN
ELSE
  IF (Q .LT. ZERO) THEN
     R = P
  ELSE
     R = ONE - P
  END IF
  IF (R .LE. ZERO) THEN
     IFAULT = 1
     PPND16 = ZERO
     RETURN
  END IF
  R = SQRT(-LOG(R))
  IF (R .LE. SPLIT2) THEN
     R = R - CONST2
     PPND16 = (((((((C7 * R + C6) * R + C5) * R + C4) * R + C3) &
        * R + C2) * R + C1) * R + C0) / &
        (((((((D7 * R + D6) * R + D5) * R + D4) * R + D3) &
        * R + D2) * R + D1) * R + ONE)
  ELSE
     R = R - SPLIT2
     PPND16 = (((((((E7 * R + E6) * R + E5) * R + E4) * R + E3) &
              * R + E2) * R + E1) * R + E0) / &
              (((((((F7 * R + F6) * R + F5) * R + F4) * R + F3) &
             * R + F2) * R + F1) * R + ONE)
  END IF
  IF (Q .LT. ZERO) PPND16 = - PPND16
  RETURN
END IF

END FUNCTION PPND16

!------------------------------------------------------------------------

function d2norm(a,b,rho)
! function d2norm calculates the density of bivariate
! normal distribution

INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL(dp),INTENT(in)::a,b,rho
REAL(dp)::d2norm

REAL(dp),PARAMETER::pi=3.141592653589793d0,pi_2=pi*2.0d0,two=2.0d0

d2norm=dexp((a*a-two*a*b*rho+b*b)/(two*rho*rho-two))/(pi_2*dsqrt(1.0d0-rho*rho))

end function d2norm

!------------------------------------------------------------------------



end module MODNormal

!------------------------------------------------------------------------

module MODpolyACM

! The Fortran code in the module MODpolyACM are programmed by Guangjian Zhang.
! His contact information: 390 Corbett Family Hall, Notre Dame, IN 46556.
! Email gzhang3@nd.edu and phone 574-631-3751.
! The code is under the license GPL-2 | GPL-3, which
! can be found at https://www.gnu.org/licenses/gpl-3.0.en.html
! The R code computes polychoric correlations with the two stage method (Olsson, 1979) and the
! asymptotic covariance matrix with an estimating equation method (Christoffersson & Gunsjo, 1996).
! Some further details are also given in Joreskog (1994).
! Christoffersson & Gunsjo (1996). A short note on the estimation of the asymptotic covariance matrix for
! polychoric correlations. Psychometrika, 61, 173-175.
! Joreskog (1994). On the estimation of polychoric correlations and their asymptotic covariance matrix. Psychometrika,
! 59, 381-389.
! Olsson (1979). Maximum likelihood estimation of the polychoric correlation coefficient. Psychometrika, 44, 443-460.

    implicit none

    PUBLIC:: SubEstimatePolyACM

    private ::cor_start, SUBPreprocessing,Threshold1D, Make1Contingency, Estimate1Correlation, &
            SubComputeACM, SubMakeGamma

contains

  subroutine SUBPreprocessing (nvar,DataMatrixIn,DataMatrixOut,Category,T_NCount,IError)
! 2020-05-30, Saturday! Guangjian Zhang
! SUBPreprocessing cleans up data first.
! It will find out the real numbers of categories for each variables and remove unnecessary categories.

  implicit none
  INTEGER, INTENT(IN) ::  nvar
  INTEGER,INTENT(IN) :: DataMatrixIn(:,:)
  INTEGER,INTENT(OUT) :: DataMatrixOut(:,:)
  INTEGER,DIMENSION(:),INTENT(OUT) :: Category
  Integer, DIMENSION(1:10,nvar), INTENT(OUT) :: T_NCount
  INTEGER,DIMENSION(nvar),INTENT(out)::IError

  Integer :: NCount(10)
  INTEGER :: i, j, IMin, nCase, k, kreal


  Category=0
  IError=0
  T_NCount=0
  DataMatrixOut = DataMatrixIn


    DO j = 1, nvar
    IMin = MINVAL(datamatrixIn(:,j))
    Category(j) = MAXVAL(datamatrixIn(:,j))  - IMin + 1
    if (Category(j)==1) IError(j) = -1 ! zero variation
    if (Category(j)> 9) IError(j) = -10 ! More than 10 categories
    DataMatrixOut(:,j) = DataMatrixIn(:,j) - IMin + 1
    END DO


    if (any(IError /= 0)) then
    return
    END IF

    ! remove empty categories

    nCase=SIZE(datamatrixIn(:,1))

    DO j = 1, nvar
        Ncount = 0
       do i = 1, nCase
         Ncount(DataMatrixOut(i,j)) = Ncount(DataMatrixOut(i,j)) + 1
       end do

       if (all(Ncount(1:Category(j))>0)) then
        T_NCount(1:Category(j),j) = Ncount(1:Category(j))
        cycle
       end if


       kreal = 0

       do k = 1, Category(j)

        if (NCount(k)>0) then
          kreal = kreal+1
          T_NCount(kreal,j) = NCount(k)
        else
          where (DataMatrixOut(:,j) > kreal) DataMatrixOut(:,j) = DataMatrixOut(:,j) - 1
        end if

       end do
       Category(j) = kreal

    End do



  end subroutine SUBPreprocessing

!-------------------------------------------------------
  subroutine Threshold1D(NCount, NCase, NCategory, Threshold, Ierror)
    use MODNormal, only: PPND16
    implicit none

    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)


    integer, dimension(:), intent(in) :: NCount
    integer, intent(in) :: Ncase, NCategory
    real(dp), dimension(0:10,2), intent(out) :: Threshold !(:,1): without a constant; (:,2) with a constant
    integer, dimension(1:10,2), intent(out) :: Ierror

    real(dp):: xtemp, xtemp1d(NCategory)
    integer:: i

    xtemp = 0.0d0

    Threshold = 0.0d0
    Ierror = 0

    Threshold(0,1:2) = -1.0d10
    Threshold(NCategory,1:2) = 1.0d10


    xtemp1d = dble(NCount(1:NCategory)) / dble(NCase)

    do i = 1, NCategory - 1
      xtemp = xtemp + xtemp1d(i)
      Threshold(i,1)=ppnd16(xtemp,Ierror(i,1))
    end do

! if a small constant 1/(nr*nc) is added to each entry
    xtemp1d = (dble(NCount(1:NCategory)) + 1.0d0 / dble(NCategory)) / dble(NCase+1)
    xtemp = 0.0d0
    do i = 1, NCategory - 1
      xtemp = xtemp + xtemp1d(i)
      Threshold(i,2)=ppnd16(xtemp,Ierror(i,2))
    end do



  end subroutine Threshold1D


  subroutine Make1Contingency(VariableI,VariableJ,Contingency)
    implicit none
    INTEGER,INTENT(in)::VariableI(:),VariableJ(:)
    INTEGER,INTENT(out)::Contingency(:,:)

    integer:: Ncase
    integer:: l

    ncase=SIZE(VariableI)


    Contingency = 0

    do l = 1, Ncase
       Contingency(VariableI(l),VariableJ(l)) = Contingency(VariableI(l),VariableJ(l)) + 1
    end do


  end subroutine Make1Contingency

!-------------------------------------------------------
subroutine Estimate1Correlation (Contingency, CategoryI, CategoryJ, IAdjust, ThresholdI2, ThresholdJ2, xContingency, &
            CellProbability, CellDerivative, Correlation, information,Ierror)
    use MODNormal, only: PPND16

! more options of dealing with empty cells


    implicit none

    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

    Integer, dimension(:,:), intent(in) :: Contingency
    Integer, intent(in) :: CategoryI, CategoryJ, IAdjust
    real(dp), dimension(0:10,2), intent(in) :: ThresholdI2, ThresholdJ2
    real(dp), dimension(:,:), intent(out) :: xContingency, CellProbability, CellDerivative
    real(dp), dimension(2), intent(out) :: Correlation
    real(dp), intent(out) :: information
    integer, dimension(2), intent(out) :: Ierror


    real(dp) :: xTable(CategoryI,CategoryJ), ThresholdI(0:10), ThresholdJ(0:10)
    real(dp) :: PearsonCorrelation, fx, dfx, dfxold, xscore, xtemp, xcase
    real(dp) :: xtolerance = 1.0D-5

    integer:: i
    integer :: NIteration = 20

    Ierror = 0


    xcase = dble(sum(Contingency))
    xContingency = dble (Contingency)


 !   if ((any(Contingency(1:CategoryI,1:CategoryJ)==0)) .and. (IAdjust /= 0) ) then
 !          xContingency(1:CategoryI,1:CategoryJ) = xContingency(1:CategoryI,1:CategoryJ) + 1.0d0 / dble(CategoryI * CategoryJ)
  !         Ierror(1) = 1
  !         xcase = xcase + 1
  !         ThresholdI(0:10) = ThresholdI2(0:10,2)
  !         ThresholdJ(0:10) = ThresholdJ2(0:10,2)
  !    else
  !         ThresholdI(0:10) = ThresholdI2(0:10,1)
   !        ThresholdJ(0:10) = ThresholdJ2(0:10,1)
   ! end if

if (any(Contingency(1:CategoryI,1:CategoryJ)==0)) then

    select case (IAdjust)

    case (0)

    ThresholdI(0:10) = ThresholdI2(0:10,1)
    ThresholdJ(0:10) = ThresholdJ2(0:10,1)

    case (1)

    xContingency(1:CategoryI,1:CategoryJ) = xContingency(1:CategoryI,1:CategoryJ) + 1.0d0 / dble(CategoryI * CategoryJ)

    Ierror(1) = 1
    xcase = xcase + 1.0d0
    ThresholdI(0:10) = ThresholdI2(0:10,2)
    ThresholdJ(0:10) = ThresholdJ2(0:10,2)

    case (2)

    xContingency(1:CategoryI,1:CategoryJ) = xContingency(1:CategoryI,1:CategoryJ) + 0.1d0
    Ierror(1) = 1
    xcase = xcase + dble(CategoryI * CategoryJ)* 0.1d0

    call ComputeThresholdIJ

    case (3)

    xContingency(1:CategoryI,1:CategoryJ) = xContingency(1:CategoryI,1:CategoryJ) + 0.5d0
    Ierror(1) = 1
    xcase = xcase + dble(CategoryI * CategoryJ)* 0.5d0

    call ComputeThresholdIJ

    case (11)

     where (Contingency(1:CategoryI,1:CategoryJ)==0) xContingency(1:CategoryI,1:CategoryJ) = 1.0d0 / dble(CategoryI * CategoryJ)

    Ierror(1) = 1
    xcase = xcase + dble(count(Contingency(1:CategoryI,1:CategoryJ)==0))* 1.0d0 / dble(CategoryI * CategoryJ)

    call ComputeThresholdIJ

    case (12)

     where (Contingency(1:CategoryI,1:CategoryJ)==0) xContingency(1:CategoryI,1:CategoryJ) = 0.1d0

    Ierror(1) = 1
    xcase = xcase + dble(count(Contingency(1:CategoryI,1:CategoryJ)==0)) * 0.1d0

    call ComputeThresholdIJ

    case (13)

     where (Contingency(1:CategoryI,1:CategoryJ)==0) xContingency(1:CategoryI,1:CategoryJ) = 0.5d0

    Ierror(1) = 1
    xcase = xcase + dble(count(Contingency(1:CategoryI,1:CategoryJ)==0)) * 0.5d0

    call ComputeThresholdIJ


    case default

    end select

    else
    ThresholdI(0:10) = ThresholdI2(0:10,1)
    ThresholdJ(0:10) = ThresholdJ2(0:10,1)
end if

    xTable = xContingency(1:CategoryI,1:CategoryJ)
    call cor_start(xTable,CategoryI, CategoryJ,PearsonCorrelation,Ierror(2))


xtemp = PearsonCorrelation
dfxold = 1.0d0

do i = 1, NIteration

call updateCorrelation(xtemp,CategoryI, CategoryJ,ThresholdI, ThresholdJ, CellProbability, CellDerivative)

fx =  SUM(xTable*dlog(CellProbability(1:CategoryI,1:CategoryJ)))/xcase
dfx = SUM(xTable / CellProbability(1:CategoryI,1:CategoryJ)*CellDerivative(1:CategoryI,1:CategoryJ))/xcase
xscore = -SUM( CellDerivative(1:CategoryI,1:CategoryJ)**2 / CellProbability(1:CategoryI,1:CategoryJ))

! write (unit=3,fmt='("NIteration, Correlation, fx, dfx, and xscore = ",i4,4E15.4)') i, xtemp, fx, dfx, xscore

xtemp = xtemp - dfx / xscore

if ((abs(dfx)<xtolerance) .and. (abs(dfxold)<xtolerance)) exit

dfxold = dfx

end do


Correlation(1) = PearsonCorrelation
Correlation(2) = xtemp
information =  - xscore ! possible error? 2020-07-12
Ierror(2) = i

contains
!-------------------------------
subroutine ComputeThresholdIJ
    integer :: i,j,DummyError
    real(dp):: xtemp, xtemp1d(10)


     ThresholdI(0) = -1.0d10
     ThresholdJ(0) = -1.0d10
     ThresholdI(CategoryI) = 1.0d10
     ThresholdJ(CategoryJ) = 1.0d10

    xtemp=0.0d0
    xtemp1d = 0.0d0
    xtemp1d(1:CategoryI) = sum(xContingency(1:CategoryI,1:CategoryJ), dim=2) /xcase


    do i = 1, CategoryI - 1
      xtemp = xtemp + xtemp1d(i)
      ThresholdI(i)=ppnd16(xtemp,DummyError)
    end do



    xtemp=0.0d0
    xtemp1d = 0.0d0
    xtemp1d(1:CategoryJ) = sum(xContingency(1:CategoryI,1:CategoryJ), dim = 1) /xcase


    do j = 1, CategoryJ - 1
      xtemp = xtemp + xtemp1d(j)
      ThresholdJ(j)=ppnd16(xtemp,DummyError)
    end do

end subroutine ComputeThresholdIJ


end subroutine Estimate1Correlation
!-------------------------------------------------------
subroutine cor_start(xtable,CategoryI, CategoryJ,correlation,nerror)
! Modified on 2020-06-02, GZ, to replace allocatable arrays so it is compatible with openMP
! modified on 2004-12-03, by change table from integer to real number.
! subroutine cor_start estimates correlations from contingcy tables
! Guangjian Zhang, 2004-10-10
implicit none

INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)


REAL(dp),INTENT(in)::xtable(:,:)
Integer, intent(in) :: CategoryI, CategoryJ
REAL(dp),INTENT(out)::correlation
INTEGER,INTENT(out)::nerror
! No error                 -> nerror = 0
! either s.d.'s is zero    -> nerror = -10000
! correlation is -1        -> nerror = -20000
! correlation is +1        -> nerror = -30000

REAL(dp) ::xproduct(CategoryI, CategoryJ)
REAL(dp)::xtotal,xtemp1,xtemp2,srow,ssrow,scolumn,sscolumn,sproduct
INTEGER::nr,nc
INTEGER::i,j
LOGICAL::Ltable(CategoryI, CategoryJ)

! nr=SIZE(xtable,1)
! nc=SIZE(xtable,2)

nr = CategoryI
nc = CategoryJ

correlation = 0.0d0


! ALLOCATE (xproduct(nr,nc),ltable(nr,nc))


ltable=(xtable>0.0d0)

xtemp1=0.0d0 ! for sum
xtemp2=0.0d0 ! for sum of squares

do i=1,nr
xtemp1=xtemp1+SUM(xtable(i,:),ltable(i,:))*DBLE(i)
xtemp2=xtemp2+SUM(xtable(i,:),ltable(i,:))*(DBLE(i))**2
end do

srow=xtemp1
ssrow=xtemp2


xtemp1=0.0d0
xtemp2=0.0d0
do j=1,nc
xtemp1=xtemp1+SUM(xtable(:,j),ltable(:,j))*DBLE(j)
xtemp2=xtemp2+SUM(xtable(:,j),ltable(:,j))*(DBLE(j))**2
end do
scolumn=xtemp1
sscolumn=xtemp2

do j=1,nc
 do i=1,nr
 xproduct(i,j)=DBLE(i*j)
 end do
end do

sproduct=SUM(xproduct*xtable,ltable)


xtemp1=0 ! SSE for row
xtemp2=0 ! SSE for column

xtotal=SUM(xtable,ltable)

xtemp1=ssrow-srow**2/xtotal
xtemp2=sscolumn-scolumn**2/xtotal


if (xtemp1<1.0d-10.or.xtemp2<1.0d-10) then
nerror=-10000
return
else
nerror=0
correlation=(sproduct-srow*scolumn/xtotal)/dsqrt(xtemp1*xtemp2)
end if

IF(dabs(correlation+1.0d0)<1.0d-8) nerror=-20000
IF(dabs(correlation-1.0d0)<1.0d-8) nerror=-30000

! deallocate (xproduct,ltable)
end subroutine cor_start
!-------------------------------------------------------------------------
subroutine updateCorrelation(rho,nrow,ncolumn,Thresholdr, Thresholdc, probabilities, derivatives)
! Guangjian Zhang, 2004-10-08
! It updates probabilities and derivatives based on procedures described by OLSSON (1979)
! in the Psychometrika paper, vol 44, 443-460
use MODNormal, only: bvn, d2norm

implicit none

INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)


REAL(dp),INTENT(in):: rho
integer, intent(in):: nrow,ncolumn
real(dp), dimension(0:10), intent(in) :: Thresholdr, Thresholdc
real(dp), dimension(:,:), intent(out):: probabilities, derivatives

REAL(dp)::xlow(2),xupper(2),temp(2,2)
INTEGER::inx(2)
INTEGER::i,j


do j=1,ncolumn ! ncolumn is declared in polychoric_data
   do i=1,nrow ! nrow is declared in polychoric_data

   ! update probability first
   inx=2
   IF(i==1) inx(1)=0
   IF(j==1) inx(2)=0
   IF(i==nrow) inx(1)=1
   IF(j==ncolumn) inx(2)=1
   xlow(1)=thresholdr(i-1) ! xlow(1) -> row ; xlow(2) -> column
   xlow(2)=thresholdc(j-1)
   xupper(1)=thresholdr(i)
   xupper(2)=thresholdc(j)
   probabilities(i,j)=bvn(xlow,xupper,inx,rho) ! probabilities declared in polychoric_data, bvn is in module distribution



   temp(1,1)=d2norm(thresholdr(i-1),thresholdc(j-1),rho)
   temp(1,2)=d2norm(thresholdr(i-1),thresholdc(j),rho)
   temp(2,1)=d2norm(thresholdr(i),thresholdc(j-1),rho)
   temp(2,2)=d2norm(thresholdr(i),thresholdc(j),rho)
   derivatives(i,j)=temp(2,2)+temp(1,1)-temp(1,2)-temp(2,1)

   end do ! nrow
end do ! ncolumn

end subroutine updateCorrelation

!---------------------------------------------------------------
subroutine SubMakeGamma(NCategoryI, NCategoryJ, ThresholdI, ThresholdJ, &
            Probability, Contigency, Derivative, Correlation, Information, GammaMatrix, Omega)

use MODNormal, only: phi

implicit none

INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)


integer, intent(in) :: NCategoryI, NCategoryJ
integer, dimension(:,:), intent(in) :: Contigency
real(dp), dimension(0:10,2), intent(in) :: ThresholdI, ThresholdJ
real(dp), dimension(:,:), intent(in) :: Probability, Derivative
real(dp), intent(in) :: Correlation, Information
real(dp), dimension(:,:), intent(out) :: GammaMatrix
real(dp), intent(out) :: Omega

REAL(dp),PARAMETER::rootpi2=sqrt(3.141592653589793d0*2.0d0)

real(dp)::  ThresholdI1D(0:10), ThresholdJ1D(0:10), BetaI(NCategoryI-1), BetaJ(NCategoryJ-1), &
        AlphaMatrix(NCategoryI, NCategoryJ), xtempI1(NCategoryI), xtempJ1(NCategoryJ), &
        DenThresholdI1D(0:10), DenThresholdJ1D(0:10), xContigency(NCategoryI, NCategoryJ)
integer::i,j

! step 1, choose proper thresholds

    ThresholdI1D(0:10) = ThresholdI(0:10,1)
    ThresholdJ1D(0:10) = ThresholdJ(0:10,1)

! Step 2, Compute AlphaMatrix

AlphaMatrix(1:NCategoryI, 1:NCategoryJ) = Derivative(1:NCategoryI, 1:NCategoryJ)/ Probability(1:NCategoryI, 1:NCategoryJ)

! Step 3, Compute BetaJ and BetaI
   ! Step 3.1, compute BetaJ
    call ComputeBeta (NCategoryI, NCategoryJ, ThresholdI1D, ThresholdJ1D, AlphaMatrix, BetaJ)
   ! Step 3.2, compute BetaI
    call ComputeBeta (NCategoryJ, NCategoryI, ThresholdJ1D, ThresholdI1D, transpose(AlphaMatrix), BetaI)


   ! Step 4, Compute GammaMatrix and Omega
   DenThresholdI1D = 0.0d0
   DenThresholdJ1D = 0.0d0

   DenThresholdI1D(1:(NCategoryI-1)) = exp(- 0.5d0 * ThresholdI1D(1:(NCategoryI-1))**2 ) / rootpi2
   DenThresholdJ1D(1:(NCategoryJ-1)) = exp(- 0.5d0 * ThresholdJ1D(1:(NCategoryJ-1))**2 ) / rootpi2

   GammaMatrix = 0.0d0
   xtempI1 = 0.0d0
   xtempJ1 = 0.0d0

   do j =1, NCategoryI -1
     xtempI1(j) =  - sum(  BetaI(j:(NCategoryI - 1))/DenThresholdI1D(j:(NCategoryI - 1)))
   end do

   do j =1, NCategoryJ -1
     xtempJ1(j) = - sum( BetaJ(j:(NCategoryJ - 1))/DenThresholdJ1D(j:(NCategoryJ - 1)))
   end do


GammaMatrix(1:NCategoryI, 1:NCategoryJ) = AlphaMatrix(1:NCategoryI, 1:NCategoryJ)

do j = 1, NCategoryJ
    do i=1, NCategoryI
      GammaMatrix(i,j) = GammaMatrix(i,j) +  xtempI1(i) + xtempJ1(j)
    end do
end do

GammaMatrix(1:NCategoryI, 1:NCategoryJ) = GammaMatrix(1:NCategoryI, 1:NCategoryJ) / Information

 xContigency = 0.0d0
 xContigency(1:NCategoryI, 1:NCategoryJ) = dble(Contigency(1:NCategoryI, 1:NCategoryJ))
 xContigency(1:NCategoryI, 1:NCategoryJ) = xContigency(1:NCategoryI, 1:NCategoryJ) / sum(xContigency(1:NCategoryI, 1:NCategoryJ))

 Omega = sum(GammaMatrix(1:NCategoryI, 1:NCategoryJ) * xContigency(1:NCategoryI, 1:NCategoryJ))


return

contains

subroutine ComputeBeta (nrow, ncolumn, ThresholdRow, ThresholdColumn, Alpha, BetaColumn)
! It is an internal subroutine to compute BetaJ (and BetaI). Operations are done mostly column-wise.
integer, intent(in):: nrow, ncolumn
real(dp), dimension(0:10), intent(in) :: ThresholdRow, ThresholdColumn
real(dp), dimension(:,:), intent(in) :: Alpha
real(dp), dimension(:), intent(out) :: BetaColumn


REAL(dp),PARAMETER::rootpi2=sqrt(3.141592653589793d0*2.0d0)


integer :: i, j
real(dp) :: DerivativePhi2Tau(0:nrow, ncolumn-1), DensitycColumn(ncolumn-1), xtemp2D(nrow, ncolumn-1)
real(dp) :: root1mrhosquare

! Step 1, Compute the derivative of Phi_2 WRT Tau_j

root1mrhosquare = sqrt(1.0d0-correlation*correlation) ! typo? 2020-07-12

DerivativePhi2Tau = 0.0d0
DensitycColumn = 0.0d0

DensitycColumn(1:ncolumn-1) = exp(- 0.5d0 * ThresholdColumn(1:ncolumn-1)**2 ) / rootpi2
DerivativePhi2Tau (nrow,1:ncolumn-1) = DensitycColumn(1:ncolumn-1)
do j = 1, ncolumn - 1
    do i = 1, nrow -1
       DerivativePhi2Tau(i,j) = DensitycColumn(j) * phi((ThresholdRow(i)-correlation * ThresholdColumn(j))/ root1mrhosquare )
    end do
end do

! Step 2, Compute BetaColumn

xtemp2D = 0.0d0

 xtemp2D(1:nrow,1:ncolumn-1) = (DerivativePhi2Tau(1:nrow,1:ncolumn-1) - DerivativePhi2Tau(0:nrow-1,1:ncolumn-1)) &
  * (Alpha(1:nrow,1:ncolumn-1)-Alpha(1:nrow,2:ncolumn))



BetaColumn(1:ncolumn-1) = sum(xtemp2D,dim=1)



end subroutine ComputeBeta

end subroutine SubMakeGamma
!---------------------------------------------------------------

subroutine SubComputeACM(Ncase, DataI, DataJ, DataK, DataL, GammaIJ, GammaKL, OmegaIJ, OmegaKL, acmIJ)
    implicit none

    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)


    integer, intent(in) :: Ncase
    integer, dimension(Ncase), intent(in) :: DataI, DataJ, DataK, DataL
    real(dp), dimension(:,:), intent(in) :: GammaIJ, GammaKL
    real(dp), intent(in) :: OmegaIJ, OmegaKL
    real(dp), intent(out) :: acmIJ

    integer :: i

    acmIJ = 0.0d0

    do i = 1, Ncase
        acmIJ = acmIJ + GammaIJ(DataI(i),DataJ(i)) * GammaKL(DataK(i),DataL(i))
    end do
      acmIJ = acmIJ / dble(Ncase) - OmegaIJ * OmegaKL

end subroutine SubComputeACM
!-----------------------------------------------------------------------------------------


subroutine SubEstimatePolyACM (Ncase, nvar, IAdjustIn, NCore, DataIn, ThresholdOut, PolyR,IError, ACM1d)

    use omp_lib

    implicit none

    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)


    integer, intent(in) :: Ncase, nvar, IAdjustIn, NCore
    integer, dimension(Ncase, nvar),intent(in):: DataIn
    real(dp), dimension(0:10,nvar), intent(out) :: ThresholdOut
    real(dp), dimension(nvar*(nvar-1)/2), intent(out) :: PolyR
    integer, dimension(2,nvar*(nvar-1)/2), intent(out) :: IError
    real(dp), dimension((nvar*(nvar-1))*(nvar*(nvar-1)+2)/8), intent(out), optional :: ACM1D


    real(dp), allocatable:: Total_Threshold(:,:,:), Total_xContingency(:,:,:), Total_CellProbability(:,:,:), &
                        Total_CellDerivative(:,:,:), Total_Correlation(:,:), &
                        Total_Gamma(:,:,:),Total_information(:), Total_Omega(:),&
                        CorrelationT (:,:)

    integer, allocatable :: Total_Table (:,:,:), Table_Temp(:,:), DataMatrixIn(:,:), DataMatrixOut(:,:), &
                            Category(:), T_Ncount(:,:), IJ_Index(:,:), Total_Error(:), IJKL_index(:,:), &
                            IError1d(:)
    integer ::  IerrorTh(10,2),nmaxc, IAdjustDummy, NCoreDummy
    integer:: i,j,ij,kl, ijkl

    ! Step 0, take care of local variables
    ! mthread = omp_get_max_threads()

     if (NCore < 1) then
         ! NCoreDummy = mthread - 1
         NCoreDummy = 2
        else
         NCoreDummy = NCore
     end if



allocate (DataMatrixIn(ncase,nvar), DataMatrixOut(ncase,nvar), Category(nvar), IError1d(nvar), &
          T_Ncount(10,nvar),Total_Threshold(0:10,2,nvar))

     DataMatrixIn = DataIn




    ! Step 1, Pre-processing the matrix of integers

    call SUBPreprocessing (nvar,DataMatrixIn,DataMatrixOut,Category,T_NCount,IError1d)

    if (any(IError1d /= 0)) then
       Ierror(1,1:nvar) = IError1d(1:nvar)
       ThresholdOut = 0.0d0
       PolyR = 0.0d0
       return
    end if


    ! Step 2, Univariate Compute threshold estimates in two situations (IAdjust=0 or IAdjust=1)

 !  !$omp parallel do

     do j = 1, nvar
       call Threshold1D(T_NCount(:,j), NCase, Category(j), Total_Threshold(0:10,1:2, j), IerrorTh)
     end do
! !$omp end parallel do


    ! Step 3, Bivariate case

        nmaxc = maxval(Category)
       allocate (IJ_Index(2,nvar*(nvar-1)/2), Total_Table (nmaxc,nmaxc,nvar*(nvar-1)/2),Table_Temp(nmaxc,nmaxc))
i=0 ! added to please a gfortran compiler
ij = 0
IJ_Index = 0
Total_Table = 0
Table_Temp = 0
do j = 2, nvar
    do i = 1, j-1
      ij = ij + 1
       IJ_Index(1,ij) = i
       IJ_Index(2,ij) = j
    end do
end do


         ! Step 3.1 Make contingency tables

 ! !$omp parallel do


do ij = 1, nvar*(nvar-1)/2
    call Make1Contingency(DataMatrixOut(:,IJ_Index(1,ij)),DataMatrixOut(:,IJ_Index(2,ij)),Table_Temp)
    Total_Table(:,:,ij) = Table_Temp
end do

! !$omp end parallel do


         ! Step 3.2 Estimate polychoric correlations


allocate (Total_xContingency(nmaxc,nmaxc,nvar*(nvar-1)/2), Total_CellProbability(nmaxc,nmaxc,nvar*(nvar-1)/2), &
                        Total_CellDerivative(nmaxc,nmaxc,nvar*(nvar-1)/2), Total_Correlation(2,nvar*(nvar-1)/2),&
                        Total_information(nvar*(nvar-1)/2))



IAdjustDummy = IAdjustIn


  call omp_set_num_threads(NCoreDummy)


  !$omp parallel do


 do ij = 1, nvar*(nvar-1)/2
!   ij = 302

    call Estimate1Correlation (Total_Table(:,:,ij), Category(IJ_Index(1,ij)), Category(IJ_Index(2,ij)), IAdjustDummy,&
    Total_Threshold(:,:,IJ_Index(1,ij)), Total_Threshold(:,:,IJ_Index(2,ij)), &
    Total_xContingency(:,:,ij), Total_CellProbability(:,:,ij),Total_CellDerivative(:,:,ij), &
    Total_Correlation(:,ij), Total_information(ij), IError(:,ij))

 end do

 !$omp end parallel do


         ! Step 3.3 if not computing ACM, deallocate arrays, and return
            ThresholdOut = 0.0d0
         do j = 1, nvar
            ThresholdOut(0:10,j) = Total_Threshold(0:10,1,j)
         end do

            allocate (CorrelationT(nvar*(nvar-1)/2,2))
            CorrelationT = transpose(Total_Correlation)

            PolyR(1:nvar*(nvar-1)/2) =  CorrelationT(1:nvar*(nvar-1)/2,2)


        if (.not.present(ACM1d)) then

            deallocate (DataMatrixIn, DataMatrixOut, Category, IError1d,T_Ncount,Total_Threshold)
            deallocate (IJ_Index, Total_Table,Table_Temp)
            deallocate (Total_xContingency, Total_CellProbability,Total_CellDerivative, Total_Correlation,&
                        Total_information, CorrelationT)

            return

         else   ! (.not.present(ACM1d))

    ! Step 4, Compute ACM

         ! Step 4.1, compute gamma and omega

allocate (Total_Gamma(nmaxc,nmaxc,nvar*(nvar-1)/2), Total_Omega(nvar*(nvar-1)/2), Total_Error(nvar*(nvar-1)/2))

IAdjustDummy = IAdjustIn

do ij = 1, nvar*(nvar-1)/2
!

call SubMakeGamma(Category(IJ_Index(1,ij)), Category(IJ_Index(2,ij)), &
Total_Threshold(:,:,IJ_Index(1,ij)),Total_Threshold(:,:,IJ_Index(2,ij)), &
Total_CellProbability(:,:,ij), Total_Table(:,:,ij), &
Total_CellDerivative(:,:,ij), Total_Correlation(2,ij), Total_information(ij), Total_Gamma(:,:,ij), &
Total_Omega(ij))

end do

allocate (IJKL_index(2,(nvar*(nvar-1))*(nvar*(nvar-1)+2)/8))
ACM1d = 0.0d0
IJKL_index = 0

kl = 0
ij =0
ijkl = 0
do j=1, nvar*(nvar-1)/2
 kl = kl + 1
 ij = 0
  do i = 1, j
    ij = ij + 1
    ijkl = ijkl + 1
    IJKL_index(1,ijkl) = ij
    IJKL_index(2,ijkl) = kl
  end do
end do


         ! Step 4.3 Compute ACM


   call omp_set_num_threads(NCoreDummy)

  !$omp parallel do

do i = 1, nvar*(nvar-1)*(nvar*(nvar-1)+2)/8

call SubComputeACM(Ncase, DataMatrixOut(1:Ncase, IJ_Index(1,IJKL_index(1,i))), &
    DataMatrixOut(1:Ncase, IJ_Index(2,IJKL_index(1,i))),&
    DataMatrixOut(1:Ncase, IJ_Index(1,IJKL_index(2,i))),&
    DataMatrixOut(1:Ncase, IJ_Index(2,IJKL_index(2,i))),&
    Total_Gamma(:,:,IJKL_index(1,i)), Total_Gamma(:,:,IJKL_index(2,i)),&
    Total_Omega(IJKL_index(1,i)), Total_Omega(IJKL_index(2,i)),ACM1d(i)) ! found a bug in 2020-07-24, GZ

end do

 !$omp end parallel do


    ! Step 5, deallocate arrays


            deallocate (DataMatrixIn, DataMatrixOut, Category, IError1D,T_Ncount,Total_Threshold)
            deallocate (IJ_Index, Total_Table,Table_Temp)
            deallocate (Total_xContingency, Total_CellProbability,Total_CellDerivative, Total_Correlation,&
                        Total_information)

      !      deallocate (Total_MB_NA,Total_MB_A)
            deallocate (Total_Gamma, Total_Omega, Total_Error)
            deallocate (IJKL_index)



        end if  ! (.not.present(ACM1d))

end subroutine SubEstimatePolyACM
!--------------------------------------------------------------------------------
end module MODpolyACM
!-------------------------------------------------------------------------------------

module ModPolyEnvelop
    use, intrinsic :: iso_c_binding

    implicit none
    private
    public :: polyR_f, polyACM_f

contains

!-------------------------------------------------------------------------------------
    subroutine polyACM_f(ncase, nvar, IAdjust, NCore, iRaw,  xThreshold, xPoly, IError, xACM) bind(C, name = "polyACM_f_")

 use MODpolyACM, only: SubEstimatePolyACM


        integer(kind = c_int), intent(in)        :: ncase    ! The # of participants
        integer(kind = c_int), intent(in)        :: nvar    ! # of variables,
        integer(kind = c_int), intent(in)        :: IAdjust    ! whether to add a small constant to empty entries
        integer(kind = c_int), intent(in)        :: NCore    ! whether to add a small constant to empty entries
        integer(kind = c_int), intent(in), dimension(ncase,nvar)               :: iRaw ! Likert variables
        real(kind = c_double), intent(out), dimension(nvar,nvar) :: xPoly  ! polychoric correlations
        real(kind = c_double), intent(out), dimension(11,nvar) :: xThreshold  ! threshold estimates for each Likert variables
        Integer(kind = c_int), intent(out),dimension(2,nvar*(nvar-1)/2)        :: iError   ! Errors for each contingency tables
        real(kind = c_double), intent(out), dimension(nvar*(nvar-1)/2,nvar*(nvar-1)/2) :: xACM  ! polychoric correlations



        real*8, dimension(nvar*(nvar-1)/2):: Correlation1D
        real*8, dimension(nvar*(nvar-1)*(nvar*(nvar-1)+2)/8):: ACM1D

        integer:: i,j, ij


        call SubEstimatePolyACM (Ncase, nvar, IAdjust, NCore, IRaw, xThreshold, Correlation1D,IError,ACM1D)

          xPoly = 0.0d0
          ij=0
        do j=2,nvar
            do i=1,j-1
              ij = ij + 1
              xPoly(i,j) = Correlation1D(ij)
            end do
        end do

        xPoly = xPoly + transpose(xPoly)

        do j=1,nvar
            Xpoly(j,j) = 1.0d0
        end do

        xACM=0.0d0
        ij=0

        do j=1,nvar*(nvar-1)/2
            do i=1,j
                ij=ij+1
                xACM(i,j) = ACM1D(ij)
                xACM(j,i) = ACM1D(ij)
            end do
        end do


    end subroutine polyACM_f

!--------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------
    subroutine polyR_f(ncase, nvar, IAdjust, NCore, iRaw,  xThreshold, xPoly, IError) bind(C, name = "polyR_f_")

 use MODpolyACM, only: SubEstimatePolyACM


        integer(kind = c_int), intent(in)        :: ncase    ! The # of participants
        integer(kind = c_int), intent(in)        :: nvar    ! # of variables,
        integer(kind = c_int), intent(in)        :: IAdjust    ! whether to add a small constant to empty entries
        integer(kind = c_int), intent(in)        :: NCore    ! whether to add a small constant to empty entries
        integer(kind = c_int), intent(in), dimension(ncase,nvar)               :: iRaw ! Likert variables
        real(kind = c_double), intent(out), dimension(nvar,nvar) :: xPoly  ! polychoric correlations
        real(kind = c_double), intent(out), dimension(11,nvar) :: xThreshold  ! threshold estimates for each Likert variables
        Integer(kind = c_int), intent(out),dimension(2,nvar*(nvar-1)/2)        :: iError   ! Errors for each contingency tables



        real*8, dimension(nvar*(nvar-1)/2):: Correlation1D


        integer:: i,j, ij


        call SubEstimatePolyACM (Ncase, nvar, IAdjust, NCore, IRaw, xThreshold, Correlation1D,IError)

          xPoly = 0.0d0
          ij=0
        do j=2,nvar
            do i=1,j-1
              ij = ij + 1
              xPoly(i,j) = Correlation1D(ij)
            end do
        end do

        xPoly = xPoly + transpose(xPoly)

        do j=1,nvar
            Xpoly(j,j) = 1.0d0
        end do

    end subroutine polyR_f

!--------------------------------------------------------------------------------------------


end module ModPolyEnvelop


