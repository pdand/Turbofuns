
module MODNormal

    public:: bvn, PPND16, d2norm, phi

    private :: normp

   contains

    !------------------------------------------------------------------------

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

REAL*8, INTENT(IN)            :: z
REAL*8, INTENT(OUT), OPTIONAL :: p, q, pdf

! Local variables
REAL*8, PARAMETER :: p0 = 220.2068679123761D0, p1 = 221.2135961699311D0,  &
                        p2 = 112.0792914978709D0, p3 = 33.91286607838300D0,  &
                        p4 = 6.373962203531650D0, p5 = .7003830644436881D0,  &
                        p6 = .3526249659989109D-01,  &
                        q0 = 440.4137358247522D0, q1 = 793.8265125199484D0,  &
                        q2 = 637.3336333788311D0, q3 = 296.5642487796737D0,  &
                        q4 = 86.78073220294608D0, q5 = 16.06417757920695D0,  &
                        q6 = 1.755667163182642D0, q7 = .8838834764831844D-1, &
                        cutoff = 7.071D0, root2pi = 2.506628274631001D0
REAL*8, PARAMETER  :: zero = 0.0D0, one = 1.0d0
REAL*8            :: zabs, expntl, pp, qq, ppdf

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
REAL*8, INTENT(IN) :: z
REAL*8             :: p

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
REAL*8,INTENT(in)::a,b,rho
REAL*8::d2norm

REAL*8,PARAMETER::pi=3.141592653589793d0,pi_2=pi*2.0d0,two=2.0d0

d2norm=dexp((a*a-two*a*b*rho+b*b)/(two*rho*rho-two))/(pi_2*dsqrt(1.0d0-rho*rho))

end function d2norm

!------------------------------------------------------------------------



end module MODNormal

!------------------------------------------------------------------------

module MODpolyACM
! Make sure to include the license information before release.

    implicit none

    PUBLIC:: SubEstimatePolyACM

    private ::cor_start, SUBPreprocessing,Threshold1D, Make1Contingency, Estimate1Correlation, &
            SubMakeB, SubComputeACM, SubMakeGamma

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
    if (Category(j)==1) IError(j) = 1 ! zero variation
    if (Category(j)> 9) IError(j) = 10 ! More than 10 categories
    DataMatrixOut(:,j) = DataMatrixIn(:,j) - IMin + 1
    END DO


    if (any(IError>0)) then
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

    integer, dimension(:), intent(in) :: NCount
    integer, intent(in) :: Ncase, NCategory
    real*8, dimension(0:10,2), intent(out) :: Threshold !(:,1): without a constant; (:,2) with a constant
    integer, dimension(1:10,2), intent(out) :: Ierror

    real*8:: xtemp, xtemp1d(NCategory)
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

    implicit none
    Integer, dimension(:,:), intent(in) :: Contingency
    Integer, intent(in) :: CategoryI, CategoryJ, IAdjust
    real*8, dimension(0:10,2), intent(in) :: ThresholdI2, ThresholdJ2
    real*8, dimension(:,:), intent(out) :: xContingency, CellProbability, CellDerivative
    real*8, dimension(2), intent(out) :: Correlation
    real*8, intent(out) :: information
    integer, dimension(2), intent(out) :: Ierror


    real*8 :: xTable(CategoryI,CategoryJ), ThresholdI(0:10), ThresholdJ(0:10)
    real*8 :: PearsonCorrelation, fx, dfx, dfxold, xscore, xtemp, xcase
    real*8 :: xtolerance = 1.0D-5

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

    xContingency(1:CategoryI,1:CategoryJ) = xContingency(1:CategoryI,1:CategoryJ) + 0.5d0
    Ierror(1) = 1
    xcase = xcase + dble(CategoryI * CategoryJ)* 0.5d0

    call ComputeThresholdIJ

    case default

    end select

    else
    ThresholdI(0:10) = ThresholdI2(0:10,1)
    ThresholdJ(0:10) = ThresholdJ2(0:10,1)
end if

    xTable = xContingency(1:CategoryI,1:CategoryJ)
    call cor_start(xTable,CategoryI, CategoryJ,PearsonCorrelation,Ierror(2))




 !   write (3, '(3i8, f15.4) ') CategoryI, CategoryJ, Ierror1, PearsonCorrelation




!    call updateCorrelation(PearsonCorrelation,CategoryI, CategoryJ,ThresholdI, ThresholdJ, CellProbability, CellDerivative)


!  write (unit=3,fmt='("xContingency = ")')
!  do i = 1, CategoryI
!    write (unit=3,fmt='(10E15.4)') (xContingency(i,j),j=1,CategoryJ)
! end do


! write (unit=3,fmt='("CellProbability = ")')
! do i = 1, CategoryI
!    write (unit=3,fmt='(10E15.4)') (CellProbability(i,j),j=1,CategoryJ)
! end do

! write (unit=3,fmt='("CellDerivative = ")')
! do i = 1, CategoryI
!    write (unit=3,fmt='(10E15.4)') (CellDerivative(i,j),j=1,CategoryJ)
! end do

! fx = SUM(xTable*dlog(CellProbability(1:CategoryI,1:CategoryJ)))/xcase
! dfx = SUM(xTable / CellProbability(1:CategoryI,1:CategoryJ)*CellDerivative(1:CategoryI,1:CategoryJ))/xcase
! xscore = -SUM( CellDerivative(1:CategoryI,1:CategoryJ)**2 / CellProbability(1:CategoryI,1:CategoryJ))

! write (unit=3,fmt='("Correlation, fx, dfx, and xscore = ",4E15.4)') PearsonCorrelation, fx, dfx, xscore

! xtemp = sign(5.0d-1 + abs(PearsonCorrelation)/2.0d0, PearsonCorrelation)

! call updateCorrelation(xtemp,CategoryI, CategoryJ,ThresholdI, ThresholdJ, CellProbability, CellDerivative)

! fx = SUM(xTable*dlog(CellProbability(1:CategoryI,1:CategoryJ)))/xcase
! dfx = SUM(xTable / CellProbability(1:CategoryI,1:CategoryJ)*CellDerivative(1:CategoryI,1:CategoryJ))/xcase
! xscore = -SUM( CellDerivative(1:CategoryI,1:CategoryJ)**2 / CellProbability(1:CategoryI,1:CategoryJ))

! write (unit=3,fmt='("Correlation, fx, dfx, and xscore = ",4E15.4)') xtemp, fx, dfx, xscore

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
information = xscore
Ierror(2) = i

contains
!-------------------------------
subroutine ComputeThresholdIJ
    integer :: i,j,DummyError
    real*8:: xtemp, xtemp1d(10)


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
REAL*8,INTENT(in)::xtable(:,:)
Integer, intent(in) :: CategoryI, CategoryJ
REAL*8,INTENT(out)::correlation
INTEGER,INTENT(out)::nerror
! No error                 -> nerror = 0
! either s.d.'s is zero    -> nerror = -10000
! correlation is -1        -> nerror = -20000
! correlation is +1        -> nerror = -30000

REAL*8 ::xproduct(CategoryI, CategoryJ)
REAL*8::xtotal,xtemp1,xtemp2,srow,ssrow,scolumn,sscolumn,sproduct
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
! USE polychoric_data
use MODNormal, only: bvn, d2norm

implicit none

REAL*8,INTENT(in):: rho
integer, intent(in):: nrow,ncolumn
real*8, dimension(0:10), intent(in) :: Thresholdr, Thresholdc
real*8, dimension(:,:), intent(out):: probabilities, derivatives

REAL*8::xlow(2),xupper(2),temp(2,2)
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



   ! upate derivatives next
   ! I will use -1.0d10 to denote negative infinity and 1.0d10 to denote postive infinity
   ! It slows the program, if it is an issue, I will rewrite the program
   ! It will cause the underflow, but it will gives me an zero
   temp(1,1)=d2norm(thresholdr(i-1),thresholdc(j-1),rho)
   temp(1,2)=d2norm(thresholdr(i-1),thresholdc(j),rho)
   temp(2,1)=d2norm(thresholdr(i),thresholdc(j-1),rho)
   temp(2,2)=d2norm(thresholdr(i),thresholdc(j),rho)
   derivatives(i,j)=temp(2,2)+temp(1,1)-temp(1,2)-temp(2,1)

   end do ! nrow
end do ! ncolumn

end subroutine updateCorrelation

!----------------------------------------------------------------

subroutine SubMakeB (NCountJ, NCategoryJ, IAdjust, ThresholdJ, MBJ, Ierror)

    implicit none

    integer, dimension(:), intent(in) :: NCountJ
    integer, intent(in) :: NCategoryJ, IAdjust
    real*8, dimension(0:10), intent(in) :: ThresholdJ
    real*8, dimension(:,:), intent(out) :: MBJ
    integer, intent(out) :: Ierror


    REAL*8,PARAMETER::rootpi2=sqrt(3.141592653589793d0*2.0d0)

    real*8 :: xa(NCategoryJ-1), xb(NCategoryJ-2), xprobJ(NCategoryJ), density1norm(NCategoryJ-1), &
              xinverse2d(NCategoryJ-1, NCategoryJ-1),  xcase, MBJT(NCategoryJ-1, NCategoryJ)
    integer :: j !, i

    Ierror = 0
    MBJT = 0.0d0

   density1norm = 0.0d0
   density1norm = exp(- 0.5d0 * ThresholdJ(1:(NCategoryJ-1))**2 ) / rootpi2

   xcase = 0.0d0
   xprobJ = 0.0d0

   xcase = dble (sum(NCountJ(1:NCategoryJ)))
   xprobJ(1:NCategoryJ) = dble(NCountJ(1:NCategoryJ))

   if (IAdjust /= 0) then
    xcase = xcase + 1.0d0
    xprobJ = xprobJ + 1.0d0 / dble(NCategoryJ)
   end if

   xprobJ = xprobJ / xcase

   xa (1: (NCategoryJ -1 )) = ( 1.0d0 / xprobJ(1:(NCategoryJ-1)) + 1.0d0 / xprobJ(2 : NCategoryJ)) &
       * density1norm (1: (NCategoryJ-1)) * density1norm (1: (NCategoryJ-1))

   if (NCategoryJ == 2) then
    MBJT(1,1:2) =  1.0d0 / xprobJ(1:2)
    MBJT(1,1:2) = density1norm(1) * MBJ(1,1:2) / xa(1)
    MBJT(1,2) = - MBJ(1,2)
    MBJ(1:NCategoryJ, 1: NCategoryJ-1) = transpose(MBJT(1:NCategoryJ-1, 1: NCategoryJ))

    return

   end if ! (NCategoryJ == 2)

    xb(1: (NCategoryJ - 2)) = - density1norm (1: (NCategoryJ-2)) * density1norm (2: (NCategoryJ-1)) &
                              / xprobJ(2 : (NCategoryJ-1))

    xinverse2d = 0.0d0

    call InvSymTridiag(xa, xb, NCategoryJ - 1, xinverse2d, Ierror)

    MBJT(:,1) = xinverse2d(:,1) * density1norm(1) / xprobJ(1)
    MBJT(:,NCategoryJ) = - xinverse2d(:,NCategoryJ-1) * density1norm(NCategoryJ-1) / xprobJ(NCategoryJ)

    do j = 2, NCategoryJ - 1
    MBJT(:,j) = ( - xinverse2d(:, j - 1) * density1norm(j-1) + xinverse2d(:, j) * density1norm(j)) / xprobJ(j)
    end do

    MBJ(1:NCategoryJ, 1: NCategoryJ-1) = transpose(MBJT(1:NCategoryJ-1, 1: NCategoryJ))


!    write (unit=3,fmt='("xa",10E15.4)') xa
!    write (unit=3,fmt='("xb",10E15.4)') xb
!    write (unit=3,fmt='(i8)') Ierror

!    write (unit=3,fmt='("xinverse")')
!    do i = 1, NCategoryJ - 1
!    write (unit=3,fmt='(10E15.4)') xinverse2d (i,:)
!    end do

 !   write (unit=3,fmt='("densit1norm")')
 !   write (unit=3,fmt='(10E15.4)') density1norm

 !   write (unit=3,fmt='("MBJ")')
 !   do i = 1, NCategoryJ - 1
 !   write (unit=3,fmt='(10E15.4)') MBJ (i,:)
 !   end do


    contains

    subroutine InvSymTridiag(xa, xb, n, xinverse, Ierror)
        ! This is an internal subroutine
        real*8, intent(in) :: xa(1:n), xb(1:n-1)
        integer, intent(in) :: n
        real*8, intent(out) :: xinverse(n,n)
        integer, intent(out) :: Ierror

        real*8 :: theta(0:n), phi(1:(n+1)), xtemp
        integer:: i,j,k

        xinverse = 0.0d0
        Ierror = 0

        theta= 0.0d0
        phi = 0.0d0

        theta(0) = 1.0d0
        theta(1) = xa(1)

        do i=2,n
           theta(i) = xa(i) * theta(i-1) - xb(i-1) * xb(i-1) * theta(i-2)
!           write (unit=3, fmt='(i8, 5f15.4)') i, theta(i-2:i), xa(i), xb(i-1)
        end do

        if (abs(theta(n))<1.0d-10) then ! The determinant of the symmetric tri-diagonal matrix is zero.
            Ierror = 1
            return
        end if

        phi(n+1) = 1.d0
        phi(n) = xa(n)

        do i = n-1, 1, -1
            phi(i) = xa(i) * phi(i+1) - xb(i) * xb(i) * phi(i+2)
        end do


        do j = 1, n
            do i = 1, j
                 if (i == j) then
                      xinverse(i,j) = theta(i-1) * phi(j+1)
                    else

                        xtemp = 1.0d0
                      do k = i, j -1
                        xtemp = xtemp * xb(k)
                      end do ! k

                    if (mod(i+j,2) == 0) then
                          xinverse(i,j) = xtemp * theta(i-1) * phi(j+1)
                        else
                         xinverse(i,j) = - xtemp * theta(i-1) * phi(j+1)
                    end if ! (mod(i+j,2) == 0)

                       xinverse(j,i) = xinverse(i,j)

                 end if ! (i == j)

            end do ! i
        end do ! j

        xinverse = xinverse / theta(n)

    end subroutine InvSymTridiag

end subroutine SubMakeB

!---------------------------------------------------------------
subroutine SubMakeGamma(NCategoryI, NCategoryJ, IAdjust, IEmpty, ThresholdI, ThresholdJ, BMatrix0I, BMatrix0J, &
            BMatrix1I, BMatrix1J, Probability, Derivative, Correlation, Information, GammaMatrix, Omega)

use MODNormal, only: phi

implicit none

integer, intent(in) :: NCategoryI, NCategoryJ, IAdjust, IEmpty
real*8, dimension(0:10,2), intent(in) :: ThresholdI, ThresholdJ
real*8, dimension(:,:), intent(in) :: BMatrix0I, BMatrix0J, BMatrix1I, BMatrix1J, Probability, Derivative
real*8, intent(in) :: Correlation, Information
real*8, dimension(:,:), intent(out) :: GammaMatrix
real*8, intent(out) :: Omega

real*8::  ThresholdI1D(0:10), ThresholdJ1D(0:10), BetaI(NCategoryI-1), BetaJ(NCategoryJ-1), &
        AlphaMatrix(NCategoryI, NCategoryJ), xtempI1(NCategoryI), xtempJ1(NCategoryJ), &
        BMatrixI(NCategoryI,NCategoryI-1), BMatrixJ(NCategoryJ,NCategoryJ-1)
integer::i,j

! step 1, choose proper thresholds according to whether a small constant is added

if ((IAdjust /=0) .and. (IEmpty /=0) ) then
    ThresholdI1D(0:10) = ThresholdI(0:10,2)
    ThresholdJ1D(0:10) = ThresholdJ(0:10,2)
    BMatrixI(1:NCategoryI,1:NCategoryI-1) = BMatrix1I(NCategoryI,NCategoryI-1)
    BMatrixJ(1:NCategoryJ,1:NCategoryJ-1) = BMatrix1J(NCategoryJ,NCategoryJ-1)
    else
    ThresholdI1D(0:10) = ThresholdI(0:10,1)
    ThresholdJ1D(0:10) = ThresholdJ(0:10,1)
    BMatrixI(1:NCategoryI,1:NCategoryI-1) = BMatrix0I(NCategoryI,NCategoryI-1)
    BMatrixJ(1:NCategoryJ,1:NCategoryJ-1) = BMatrix0J(NCategoryJ,NCategoryJ-1)

end if

! Step 2, Compute AlphaMatrix

AlphaMatrix(1:NCategoryI, 1:NCategoryJ) = Derivative(1:NCategoryI, 1:NCategoryJ)/ Probability(1:NCategoryI, 1:NCategoryJ)

! Step 3, Compute BetaJ and BetaI
   ! Step 3.1, compute BetaJ
    call ComputeBeta (NCategoryI, NCategoryJ, ThresholdI1D, ThresholdJ1D, AlphaMatrix, BetaJ)
   ! Step 3.2, compute BetaI
    call ComputeBeta (NCategoryJ, NCategoryI, ThresholdJ1D, ThresholdI1D, transpose(AlphaMatrix), BetaI)

! write (unit=3, fmt='("BetaJ", 10E15.4)') BetaJ
! write (unit=3, fmt='("BetaI", 10E15.4)') BetaI


   ! Step 4, Compute GammaMatrix and Omega
   GammaMatrix = 0.0d0
   xtempI1 = 0.0d0
   xtempJ1 = 0.0d0


   do j =1, NCategoryI -1
     xtempI1(1:NcategoryI) = xtempI1(1:NcategoryI) + BMatrixI(1:NcategoryI, j) * BetaI(j)
   end do

   do j =1, NCategoryJ -1
     xtempJ1(1:NcategoryJ) = xtempJ1(1:NcategoryJ) + BMatrixJ(1:NcategoryJ, j) * BetaJ(j)
   end do

GammaMatrix(1:NCategoryI, 1:NCategoryJ) = AlphaMatrix(1:NCategoryI, 1:NCategoryJ)

do j = 1, NCategoryJ
    do i=1, NCategoryI
      GammaMatrix(i,j) = GammaMatrix(i,j) +  xtempI1(i) + xtempJ1(j)
    end do
end do

GammaMatrix(1:NCategoryI, 1:NCategoryJ) = GammaMatrix(1:NCategoryI, 1:NCategoryJ) / Information

Omega = sum(GammaMatrix(1:NCategoryI, 1:NCategoryJ) * Probability(1:NCategoryI, 1:NCategoryJ))

return

contains

subroutine ComputeBeta (nrow, ncolumn, ThresholdRow, ThresholdColumn, Alpha, BetaColumn)
! It is an internal subroutine to compute BetaJ (and BetaI). Operations are done mostly column-wise.
integer, intent(in):: nrow, ncolumn
real*8, dimension(0:10), intent(in) :: ThresholdRow, ThresholdColumn
real*8, dimension(:,:), intent(in) :: Alpha
real*8, dimension(:), intent(out) :: BetaColumn


REAL*8,PARAMETER::rootpi2=sqrt(3.141592653589793d0*2.0d0)


integer :: i, j
real*8 :: DerivativePhi2Tau(0:nrow, ncolumn-1), DensitycColumn(ncolumn-1), xtemp2D(nrow, ncolumn-1)
real*8 :: root1mrhosquare

! Step 1, Compute the derivative of Phi_2 WRT Tau_j

root1mrhosquare = sqrt(1-correlation*correlation)

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

! do i=1,nrow
!    write (unit=3, fmt='(10E15.4)') xtemp2D(i,:)
! end do


end subroutine ComputeBeta

end subroutine SubMakeGamma
!---------------------------------------------------------------

subroutine SubComputeACM(Ncase, DataI, DataJ, DataK, DataL, GammaIJ, GammaKL, OmegaIJ, OmegaKL, acmIJ)
    implicit none
    integer, intent(in) :: Ncase
    integer, dimension(Ncase), intent(in) :: DataI, DataJ, DataK, DataL
    real*8, dimension(:,:), intent(in) :: GammaIJ, GammaKL
    real*8, intent(in) :: OmegaIJ, OmegaKL
    real*8, intent(out) :: acmIJ

    integer :: i

    acmIJ = 0.0d0

    do i = 1, Ncase
        acmIJ = acmIJ + GammaIJ(DataI(i),DataJ(i)) * GammaKL(DataK(i),DataL(i))
    end do
        acmIJ = acmIJ / dble(Ncase) + OmegaIJ * OmegaKL

end subroutine SubComputeACM
!-----------------------------------------------------------------------------------------


subroutine SubEstimatePolyACM (Ncase, nvar, IAdjustIn, NCore, DataIn, ThresholdOut, PolyR,IError, ACM1d)

    use omp_lib

    implicit none

    integer, intent(in) :: Ncase, nvar, IAdjustIn, NCore
    integer, dimension(Ncase, nvar),intent(in):: DataIn
    real*8, dimension(0:10,nvar), intent(out) :: ThresholdOut
    real*8, dimension(nvar*(nvar-1)/2), intent(out) :: PolyR
    integer, dimension(2,nvar*(nvar-1)/2), intent(out) :: IError
    real*8, dimension((nvar*(nvar-1))*(nvar*(nvar-1)+2)/8), intent(out), optional :: ACM1D


    real*8, allocatable:: Total_Threshold(:,:,:), Total_xContingency(:,:,:), Total_CellProbability(:,:,:), &
                        Total_CellDerivative(:,:,:), Total_Correlation(:,:), Total_MB_NA(:,:,:), &
                        Total_MB_A(:,:,:), Total_Gamma(:,:,:),Total_information(:), Total_Omega(:),&
                        CorrelationT (:,:)

    integer, allocatable :: Total_Table (:,:,:), Table_Temp(:,:), DataMatrixIn(:,:), DataMatrixOut(:,:), &
                            Category(:), T_Ncount(:,:), IJ_Index(:,:), Total_Error(:), IJKL_index(:,:), &
                            IError1d(:)
    integer ::  IerrorTh(10,2),nmaxc, mthread, IAdjustDummy, Ierror1, NCoreDummy
    integer:: i,j,ij,kl, ijkl

    ! Step 0, take care of local variables
     mthread = omp_get_max_threads()

     if (NCore < 1) then
          NCoreDummy = mthread - 1
        else
         NCoreDummy = NCore
     end if



allocate (DataMatrixIn(ncase,nvar), DataMatrixOut(ncase,nvar), Category(nvar), IError1d(nvar), &
          T_Ncount(10,nvar),Total_Threshold(0:10,2,nvar))

     DataMatrixIn = DataIn




    ! Step 1, Pre-processing the matrix of integers

    call SUBPreprocessing (nvar,DataMatrixIn,DataMatrixOut,Category,T_NCount,IError1d)


    ! Step 2, Univariate Compute threshold estimates in two situations (IAdjust=0 or IAdjust=1)

     do j = 1, nvar
       call Threshold1D(T_NCount(:,j), NCase, Category(j), Total_Threshold(0:10,1:2, j), IerrorTh)
     end do


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


do ij = 1, nvar*(nvar-1)/2
    call Make1Contingency(DataMatrixOut(:,IJ_Index(1,ij)),DataMatrixOut(:,IJ_Index(2,ij)),Table_Temp)
    Total_Table(:,:,ij) = Table_Temp
end do



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


         ! Step 4.1, Compute B matrices

         allocate (Total_MB_NA(nmaxc,nmaxc-1,nvar),Total_MB_A(nmaxc,nmaxc-1,nvar))
IAdjustDummy = 0

do j = 1, nvar
  call SubMakeB (T_NCount(:,j), Category(j), IAdjustDummy, Total_Threshold(0:10,1,j), Total_MB_NA(:,:,j), Ierror1)
end do

IAdjustDummy = 1

do j = 1, nvar
  call SubMakeB (T_NCount(:,j), Category(j), IAdjustDummy, Total_Threshold(0:10,2,j), Total_MB_A(:,:,j), Ierror1)
end do


         ! Step 4.2, compute gamma and omega

allocate (Total_Gamma(nmaxc,nmaxc,nvar*(nvar-1)/2), Total_Omega(nvar*(nvar-1)/2), Total_Error(nvar*(nvar-1)/2))

IAdjustDummy = IAdjustIn

do ij = 1, nvar*(nvar-1)/2
!

call SubMakeGamma(Category(IJ_Index(1,ij)), Category(IJ_Index(2,ij)), IAdjustDummy, Total_Error(ij), &
Total_Threshold(:,:,IJ_Index(1,ij)),Total_Threshold(:,:,IJ_Index(2,ij)), Total_MB_NA(:,:,i), &
Total_MB_NA(:,:,j), Total_MB_A(:,:,i), Total_MB_A(:,:,j), Total_CellProbability(:,:,ij),&
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
    Total_Omega(IJKL_index(1,i)), Total_Omega(IJKL_index(1,2)),ACM1d(i))

end do

 !$omp end parallel do


    ! Step 5, deallocate arrays


            deallocate (DataMatrixIn, DataMatrixOut, Category, IError1D,T_Ncount,Total_Threshold)
            deallocate (IJ_Index, Total_Table,Table_Temp)
            deallocate (Total_xContingency, Total_CellProbability,Total_CellDerivative, Total_Correlation,&
                        Total_information)

            deallocate (Total_MB_NA,Total_MB_A)
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


