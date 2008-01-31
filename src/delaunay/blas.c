#include "mf2c.h"
integer isamax_(n, sx, incx)
integer *n;
real *sx;
integer *incx;
{
    /* System generated locals */
    integer ret_val, i__1;
    real r__1;

    /* Local variables */
    static real smax;
    static integer i, ix;


/*     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --sx;

    /* Function Body */
    ret_val = 0;
    if (*n < 1) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        CODE FOR INCREMENT NOT EQUAL TO 1 */

    ix = 1;
    smax = dabs(sx[1]);
    ix += *incx;
    i__1 = *n;
    for (i = 2; i <= i__1; ++i) {
	if ((r__1 = sx[ix], dabs(r__1)) <= smax) {
	    goto L5;
	}
	ret_val = i;
	smax = (r__1 = sx[ix], dabs(r__1));
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*        CODE FOR INCREMENT EQUAL TO 1 */

L20:
    smax = dabs(sx[1]);
    i__1 = *n;
    for (i = 2; i <= i__1; ++i) {
	if ((r__1 = sx[i], dabs(r__1)) <= smax) {
	    goto L30;
	}
	ret_val = i;
	smax = (r__1 = sx[i], dabs(r__1));
L30:
	;
    }
    return ret_val;
} /* isamax_ */

doublereal sasum_(n, sx, incx)
integer *n;
real *sx;
integer *incx;
{
    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1, r__2, r__3, r__4, r__5, r__6;

    /* Local variables */
    static integer i, m, nincx;
    static real stemp;
    static integer mp1;


/*     TAKES THE SUM OF THE ABSOLUTE VALUES. */
/*     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --sx;

    /* Function Body */
    ret_val = (float)0.;
    stemp = (float)0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        CODE FOR INCREMENT NOT EQUAL TO 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i = 1; i__2 < 0 ? i >= i__1 : i <= i__1; i += i__2) {
	stemp += (r__1 = sx[i], dabs(r__1));
/* L10: */
    }
    ret_val = stemp;
    return ret_val;

/*        CODE FOR INCREMENT EQUAL TO 1 */


/*        CLEAN-UP LOOP */

L20:
    m = *n % 6;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i = 1; i <= i__2; ++i) {
	stemp += (r__1 = sx[i], dabs(r__1));
/* L30: */
    }
    if (*n < 6) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i = mp1; i <= i__2; i += 6) {
	stemp = stemp + (r__1 = sx[i], dabs(r__1)) + (r__2 = sx[i + 1], dabs(
		r__2)) + (r__3 = sx[i + 2], dabs(r__3)) + (r__4 = sx[i + 3], 
		dabs(r__4)) + (r__5 = sx[i + 4], dabs(r__5)) + (r__6 = sx[i + 
		5], dabs(r__6));
/* L50: */
    }
L60:
    ret_val = stemp;
    return ret_val;
} /* sasum_ */

/* Subroutine */ int saxpy_(n, sa, sx, incx, sy, incy)
integer *n;
real *sa, *sx;
integer *incx;
real *sy;
integer *incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*     CONSTANT TIMES A VECTOR PLUS A VECTOR. */
/*     USES UNROLLED LOOP FOR INCREMENTS EQUAL TO ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*sa == (float)0.) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
/*          NOT EQUAL TO 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	sy[iy] += *sa * sx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*        CLEAN-UP LOOP */

L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
	sy[i] += *sa * sx[i];
/* L30: */
    }
    if (*n < 4) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 4) {
	sy[i] += *sa * sx[i];
	sy[i + 1] += *sa * sx[i + 1];
	sy[i + 2] += *sa * sx[i + 2];
	sy[i + 3] += *sa * sx[i + 3];
/* L50: */
    }
    return 0;
} /* saxpy_ */

/* Subroutine */ int scopy_(n, sx, incx, sy, incy)
integer *n;
real *sx;
integer *incx;
real *sy;
integer *incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*     COPIES A VECTOR, X, TO A VECTOR, Y. */
/*     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO 1. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
/*          NOT EQUAL TO 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	sy[iy] = sx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*        CLEAN-UP LOOP */

L20:
    m = *n % 7;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
	sy[i] = sx[i];
/* L30: */
    }
    if (*n < 7) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 7) {
	sy[i] = sx[i];
	sy[i + 1] = sx[i + 1];
	sy[i + 2] = sx[i + 2];
	sy[i + 3] = sx[i + 3];
	sy[i + 4] = sx[i + 4];
	sy[i + 5] = sx[i + 5];
	sy[i + 6] = sx[i + 6];
/* L50: */
    }
    return 0;
} /* scopy_ */

doublereal sdot_(n, sx, incx, sy, incy)
integer *n;
real *sx;
integer *incx;
real *sy;
integer *incy;
{
    /* System generated locals */
    integer i__1;
    real ret_val;

    /* Local variables */
    static integer i, m;
    static real stemp;
    static integer ix, iy, mp1;


/*     FORMS THE DOT PRODUCT OF TWO VECTORS. */
/*     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    stemp = (float)0.;
    ret_val = (float)0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
/*          NOT EQUAL TO 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	stemp += sx[ix] * sy[iy];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    ret_val = stemp;
    return ret_val;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*        CLEAN-UP LOOP */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
	stemp += sx[i] * sy[i];
/* L30: */
    }
    if (*n < 5) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 5) {
	stemp = stemp + sx[i] * sy[i] + sx[i + 1] * sy[i + 1] + sx[i + 2] * 
		sy[i + 2] + sx[i + 3] * sy[i + 3] + sx[i + 4] * sy[i + 4];
/* L50: */
    }
L60:
    ret_val = stemp;
    return ret_val;
} /* sdot_ */

doublereal smach_(job)
integer *job;
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    static real huge, tiny, s, eps;


/*     SMACH COMPUTES MACHINE PARAMETERS OF FLOATING POINT */
/*     ARITHMETIC FOR USE IN TESTING ONLY.  NOT REQUIRED BY */
/*     LINPACK PROPER. */

/*     IF TROUBLE WITH AUTOMATIC COMPUTATION OF THESE QUANTITIES, */
/*     THEY CAN BE SET BY DIRECT ASSIGNMENT STATEMENTS. */
/*     ASSUME THE COMPUTER HAS */

/*        B = BASE OF ARITHMETIC */
/*        T = NUMBER OF BASE  B  DIGITS */
/*        L = SMALLEST POSSIBLE EXPONENT */
/*        U = LARGEST POSSIBLE EXPONENT */

/*     THEN */

/*        EPS = B**(1-T) */
/*        TINY = 100.0*B**(-L+T) */
/*        HUGE = 0.01*B**(U-T) */

/*     DMACH SAME AS SMACH EXCEPT T, L, U APPLY TO */
/*     DOUBLE PRECISION. */

/*     CMACH SAME AS SMACH EXCEPT IF COMPLEX DIVISION */
/*     IS DONE BY */

/*        1/(X+I*Y) = (X-I*Y)/(X**2+Y**2) */

/*     THEN */

/*        TINY = SQRT(TINY) */
/*        HUGE = SQRT(HUGE) */


/*     JOB IS 1, 2 OR 3 FOR EPSILON, TINY AND HUGE, RESPECTIVELY. */



    eps = (float)1.;
L10:
    eps /= (float)2.;
    s = eps + (float)1.;
    if (s > (float)1.) {
	goto L10;
    }
    eps *= (float)2.;

    s = (float)1.;
L20:
    tiny = s;
    s /= (float)16.;
    if (s * (float)100. != (float)0.) {
	goto L20;
    }
    tiny = tiny / eps * (float)100.;
    huge = (float)1. / tiny;

    if (*job == 1) {
	ret_val = eps;
    }
    if (*job == 2) {
	ret_val = tiny;
    }
    if (*job == 3) {
	ret_val = huge;
    }
    return ret_val;
} /* smach_ */

doublereal snrm2_(n, sx, incx)
integer *n;
real *sx;
integer *incx;
{
    /* Initialized data */

    static real zero = (float)0.;
    static real one = (float)1.;
    static real cutlo = (float)4.441e-16;
    static real cuthi = (float)1.304e19;

    /* Format strings */
    static char fmt_30[] = "";
    static char fmt_50[] = "";
    static char fmt_70[] = "";
    static char fmt_110[] = "";

    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static real xmax;
    static integer next, i, j, nn;
    static real hitest, sum;

    /* Assigned format variables */
    char *next_fmt;

    /* Parameter adjustments */
    --sx;

    /* Function Body */

/*     EUCLIDEAN NORM OF THE N-VECTOR STORED IN SX() WITH STORAGE */
/*     INCREMENT INCX . */
/*     IF    N .LE. 0 RETURN WITH RESULT = 0. */
/*     IF N .GE. 1 THEN INCX MUST BE .GE. 1 */

/*           C.L.LAWSON, 1978 JAN 08 */

/*     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE */
/*     HOPEFULLY APPLICABLE TO ALL MACHINES. */
/*         CUTLO = MAXIMUM OF  SQRT(U/EPS)  OVER ALL KNOWN MACHINES. */
/*         CUTHI = MINIMUM OF  SQRT(V)      OVER ALL KNOWN MACHINES. */
/*     WHERE */
/*         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1. */
/*         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT) */
/*         V   = LARGEST  NO.            (OVERFLOW  LIMIT) */

/*     BRIEF OUTLINE OF ALGORITHM.. */

/*     PHASE 1    SCANS ZERO COMPONENTS. */
/*     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO */
/*     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO */
/*     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M */
/*     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX. */

/*     VALUES FOR CUTLO AND CUTHI.. */
/*     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER */
/*     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS.. */
/*     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE 
*/
/*                   UNIVAC AND DEC AT 2**(-103) */
/*                   THUS CUTLO = 2**(-51) = 4.44089E-16 */
/*     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC. */
/*                   THUS CUTHI = 2**(63.5) = 1.30438E19 */
/*     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC. */
/*                   THUS CUTLO = 2**(-33.5) = 8.23181D-11 */
/*     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19 */
/*     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 / */
/*     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 / */

    if (*n > 0) {
	goto L10;
    }
    ret_val = zero;
    goto L300;

L10:
    next = 0;
    next_fmt = fmt_30;
    sum = zero;
    nn = *n * *incx;
/*                                                 BEGIN MAIN LOOP */
    i = 1;
L20:
    switch ((int)next) {
	case 0: goto L30;
	case 1: goto L50;
	case 2: goto L70;
	case 3: goto L110;
    }
L30:
    if ((r__1 = sx[i], dabs(r__1)) > cutlo) {
	goto L85;
    }
    next = 1;
    next_fmt = fmt_50;
    xmax = zero;

/*                        PHASE 1.  SUM IS ZERO */

L50:
    if (sx[i] == zero) {
	goto L200;
    }
    if ((r__1 = sx[i], dabs(r__1)) > cutlo) {
	goto L85;
    }

/*                                PREPARE FOR PHASE 2. */
    next = 2;
    next_fmt = fmt_70;
    goto L105;

/*                                PREPARE FOR PHASE 4. */

L100:
    i = j;
    next = 3;
    next_fmt = fmt_110;
    sum = sum / sx[i] / sx[i];
L105:
    xmax = (r__1 = sx[i], dabs(r__1));
    goto L115;

/*                   PHASE 2.  SUM IS SMALL. */
/*                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW. */

L70:
    if ((r__1 = sx[i], dabs(r__1)) > cutlo) {
	goto L75;
    }

/*                     COMMON CODE FOR PHASES 2 AND 4. */
/*                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW. 
*/

L110:
    if ((r__1 = sx[i], dabs(r__1)) <= xmax) {
	goto L115;
    }
/* Computing 2nd power */
    r__1 = xmax / sx[i];
    sum = one + sum * (r__1 * r__1);
    xmax = (r__1 = sx[i], dabs(r__1));
    goto L200;

L115:
/* Computing 2nd power */
    r__1 = sx[i] / xmax;
    sum += r__1 * r__1;
    goto L200;


/*                  PREPARE FOR PHASE 3. */

L75:
    sum = sum * xmax * xmax;


/*     FOR REAL OR D.P. SET HITEST = CUTHI/N */
/*     FOR COMPLEX      SET HITEST = CUTHI/(2*N) */

L85:
    hitest = cuthi / (real) (*n);

/*                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING. */

    i__1 = nn;
    i__2 = *incx;
    for (j = i; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
	if ((r__1 = sx[j], dabs(r__1)) >= hitest) {
	    goto L100;
	}
/* L95: */
/* Computing 2nd power */
	r__1 = sx[j];
	sum += r__1 * r__1;
    }
    ret_val = sqrt(sum);
    goto L300;

L200:
    i += *incx;
    if (i <= nn) {
	goto L20;
    }

/*              END OF MAIN LOOP. */

/*              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING. */

    ret_val = xmax * sqrt(sum);
L300:
    return ret_val;
} /* snrm2_ */

/* Subroutine */ int srot_(n, sx, incx, sy, incy, c, s)
integer *n;
real *sx;
integer *incx;
real *sy;
integer *incy;
real *c, *s;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i;
    static real stemp;
    static integer ix, iy;


/*     APPLIES A PLANE ROTATION. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL */
/*         TO 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	stemp = *c * sx[ix] + *s * sy[iy];
	sy[iy] = *c * sy[iy] - *s * sx[ix];
	sx[ix] = stemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       CODE FOR BOTH INCREMENTS EQUAL TO 1 */

L20:
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	stemp = *c * sx[i] + *s * sy[i];
	sy[i] = *c * sy[i] - *s * sx[i];
	sx[i] = stemp;
/* L30: */
    }
    return 0;
} /* srot_ */


/* Subroutine */ int sscal_(n, sa, sx, incx)
integer *n;
real *sa, *sx;
integer *incx;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i, m, nincx, mp1;


/*     SCALES A VECTOR BY A CONSTANT. */
/*     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO 1. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --sx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        CODE FOR INCREMENT NOT EQUAL TO 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i = 1; i__2 < 0 ? i >= i__1 : i <= i__1; i += i__2) {
	sx[i] = *sa * sx[i];
/* L10: */
    }
    return 0;

/*        CODE FOR INCREMENT EQUAL TO 1 */


/*        CLEAN-UP LOOP */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i = 1; i <= i__2; ++i) {
	sx[i] = *sa * sx[i];
/* L30: */
    }
    if (*n < 5) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i = mp1; i <= i__2; i += 5) {
	sx[i] = *sa * sx[i];
	sx[i + 1] = *sa * sx[i + 1];
	sx[i + 2] = *sa * sx[i + 2];
	sx[i + 3] = *sa * sx[i + 3];
	sx[i + 4] = *sa * sx[i + 4];
/* L50: */
    }
    return 0;
} /* sscal_ */

/* Subroutine */ int sswap_(n, sx, incx, sy, incy)
integer *n;
real *sx;
integer *incx;
real *sy;
integer *incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m;
    static real stemp;
    static integer ix, iy, mp1;


/*     INTERCHANGES TWO VECTORS. */
/*     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO 1. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL */
/*         TO 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	stemp = sx[ix];
	sx[ix] = sy[iy];
	sy[iy] = stemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*       CLEAN-UP LOOP */

L20:
    m = *n % 3;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
	stemp = sx[i];
	sx[i] = sy[i];
	sy[i] = stemp;
/* L30: */
    }
    if (*n < 3) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 3) {
	stemp = sx[i];
	sx[i] = sy[i];
	sy[i] = stemp;
	stemp = sx[i + 1];
	sx[i + 1] = sy[i + 1];
	sy[i + 1] = stemp;
	stemp = sx[i + 2];
	sx[i + 2] = sy[i + 2];
	sy[i + 2] = stemp;
/* L50: */
    }
    return 0;
} /* sswap_ */

