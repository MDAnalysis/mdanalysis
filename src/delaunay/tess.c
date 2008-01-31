
#include "tess.h"

/* Subroutine */ int tess_(mrowp, mp, np, p, mrows, ncols, ms, ns, nsim, nadj,
	 work, llfact, iwork, ierr)
integer *mrowp, *mp, *np;
real *p;
integer *mrows, *ncols, *ms, *ns, *nsim, *nadj;
real *work;
integer *llfact, *iwork, *ierr;
{
    /* System generated locals */
    integer nsim_dim1, nsim_offset, nadj_dim1, nadj_offset, p_dim1, p_offset;

    /* Local variables */
    static integer klen, llen, kfwd, lmin;
    extern /* Subroutine */ int init_();
    static integer lnew, khead, lhead, lvail, lista, lnext, listv, lmpty, ls, 
	    lv, lkount, iv1, iv2, iv3, iv4, iv5, iv6;
    extern /* Subroutine */ int nxtsim_();
    static integer kbk, len;


/* ***begin prologue tess                                      */
/* ***date written 950115   (yymmdd)                           */
/* ***keywords  Delaunay tetrahedrization tessellation         */
/* ***author Hazlewood, Carol                                  */
/*           ch04@academia.swt.edu                             */
/*           Department of Computer Science                    */
/*           Southwest Texas State University                  */
/*           San Marcos, TX 78666                              */
/*                                                             */
/* ***purpose  this routine computes a Delaunay                */
/*             tetrahedrization of the points in p.            */

/* ***description                                              */

/*  argument description ***                                   */

/*  on input                                                   */

/*   *mrowp  integer                                           */
/*             the row dimension of the array p as declared    */
/*             by the calling routine. mrowp .ge. mp           */
/*   *mp     integer                                           */
/*             the dimension of the space that contains the    */
/*             points of p.  mrowp .ge mp                      */
/*   *np     integer                                           */
/*             the number of points in p                       */
/*    p      real(mrowp*np)                                    */
/*             a one dimensional real array that contains the  */
/*             coordinates of the points to be tessellated.    */
/*             Conceptually the array p is a two dimensional   */
/*             array with mrowp rows and np columns, where     */
/*             each column of p contains the coordinates       */
/*             of one point.                                   */
/*   *mrows  integer                                           */
/*             the row dimension of arrays nsim and nadj as    */
/*             declared by the calling routine.  mrows .ge. ms */
/*   *ncols  integer                                           */
/*             the column dimension of arrays nsim and nadj as */
/*             declared by the calling routine.  this is the   */
/*             maximum number of simplices that can be         */
/*             represented.  the expected number of simplices  */
/*             in a three dimensional tessellation of np       */
/*             points is 7*np (see stoyan, kendall, and mecke, */
/*             stochastic geometry)                            */
/*    iwork  integer(*)                                        */
/*             integer array to be used as workspace.          */
/*             It should be declared by the calling routine    */
/*             to be of size 140*ncols.                        */
/*    work   real(*)                                           */
/*             real array to be used as workspace.  It should  */
/*             be declared by the calling routine to be of     */
/*             size at least 5*mrowp + mrowp*mrowp.            */

/*  on return                                                  */

/*   *ms     integer                                           */
/*             the number of vertices in a simplex.            */
/*             normally ms = mp + 1                            */
/*   *ns     integer                                           */
/*             the number of simplices in the tessellation     */
/*    nsim   integer(mrows*ncols) 
/*             a one-dimensional real array that contains      */
/*             a representation of the vertices of the         */
/*             simplices in the tessellation.   Conceptually   */
/*             it is a two dimensional integer array of mrows  */
/*             rows and ncols columns.  Each column of nsim    */
/*             represents one simplex of the tessellation and  */
/*             contains the indices in p of the vertices of    */
/*             the simplex                                     */
/*    nadj   integer(mrows,ncols) */
/*             a one-dimensional real array that contains      */
/*             a representation of the adjacent  simplices     */
/*             in the tessellation.   Conceptually it is a     */
/*             two dimensional integer array of mrows rows and */
/*             ncols columns.  Each column of nadj contains    */
/*             information about the neighbors of one simplex. */
/*             nadj(i,j) is the index in nsim of the simplex   */
/*             that is adjacent to simplex j across the facet  */
/*             that does into contain vertex i.                */
/*   *ierr   integer                                           */
/*             an error flag.  if ierr = 0, the tessellation   */
/*             is complete.  if ierr .ne. 0, an error has      */
/*             occurred                                        */

/*    mrowp,mp,p,mrows,and ncols are unchanged.                */
/*    work and iwork may be discarded.                         */

/* ***references                                               */
/* ***routines called  init,nxtsim                             */
/* ***end prologue tess                                        */





/* decompose workspace */

    /* Parameter adjustments */
    p_dim1 = *mrowp;
    p_offset = p_dim1 + 1;
    p -= p_offset;
    nadj_dim1 = *mrows;
    nadj_offset = nadj_dim1 + 1;
    nadj -= nadj_offset;
    nsim_dim1 = *mrows;
    nsim_offset = nsim_dim1 + 1;
    nsim -= nsim_offset;
    --work;
    --iwork;

    /* Function Body */
    len = *np * *llfact;
    listv = 1;
    lista = listv + *mrows;
    llen = lista + *mrows;
    lhead = llen + 1;
    lmin = lhead + 1;
    lvail = lmin + 1;
    lnew = lvail + 1;
    lnext = lnew + 1;
    ls = lnext + len;
    lv = ls + len;
    lkount = lv + len;
    klen = lkount + *np;
    khead = klen + 1;
    kfwd = khead + 1;
    kbk = kfwd + *np;
    iv1 = 1;
    iv2 = iv1 + *mp;
    iv3 = iv2 + *mp;
    iv4 = iv3 + *mp;
    iv5 = iv4 + *mp;
    iv6 = iv5 + *mp;

    init_(mrowp, mp, np, &p[p_offset], mrows, ncols, ms, ns, &nsim[
	    nsim_offset], &nadj[nadj_offset], &len, &iwork[listv], &iwork[
	    lista], &iwork[llen], &iwork[lhead], &iwork[lmin], &iwork[lvail], 
	    &iwork[lnew], &iwork[lnext], &iwork[ls], &iwork[lv], &iwork[
	    lkount], &iwork[klen], &iwork[khead], &iwork[kfwd], &iwork[kbk], &
	    work[iv1], &work[iv2], &work[iv3], &work[iv4], &work[iv5], &work[
	    iv6], &lmpty, ierr);
L1:
    if (lmpty != 0 && *ierr == 0) {
	nxtsim_(mrowp, mp, np, &p[p_offset], mrows, ncols, ms, ns, &nsim[
		nsim_offset], &nadj[nadj_offset], &iwork[listv], &iwork[lista]
		, &iwork[llen], &iwork[lhead], &iwork[lmin], &iwork[lnew], &
		iwork[lvail], &iwork[lnext], &iwork[ls], &iwork[lv], &iwork[
		lkount], &iwork[klen], &iwork[khead], &iwork[kfwd], &iwork[
		kbk], &work[iv1], &work[iv2], &work[iv3], &work[iv4], &work[
		iv5], &work[iv6], &lmpty, ierr);
	goto L1;
    }
    return 0;
} /* tess_ */




int init_(mrowp, mp, np, p, mrows, ncols, ms, ns, nsim, nadj,
          len, listv, lista, llen, lhead, lmin, lvail, lnew,
          lnext, ls, lv, lkount, klen, khead, kfwd, kbk, 
          v1, v2, v3, v4, v5, v6, lmpty, ierr)

  integer *mrowp, *mp, *np, *mrows, *ncols, *ms, *ns, *nsim, 
          *nadj, *len, *listv, *lista, *llen, *lhead, *lmin, 
          *lvail, *lnew, *lnext, *ls, *lv, *lkount, *klen, 
          *khead, *kfwd, *kbk, *lmpty, *ierr;
  real    *p, *v1, *v2, *v3, *v4, *v5, *v6;


/* This routine finds a facet in the boundary of the      */
/* convex hull of p (necessarily a facet of a simplex     */
/* in the simplicial tessellation of p) and initializes   */
/* data structures.                                       */


{
    /* System generated locals */
    integer nsim_dim1, nsim_offset, nadj_dim1, nadj_offset, 
            p_dim1, p_offset, v6_dim1, v6_offset, i__1;

    /* Local variables */
    extern  int addl_(), adds_(), initl_(), initk_(), inits_(),
	               facet1_(), onesim_();

    /* Parameter adjustments */
    v6_dim1 = *mp;
    v6_offset = v6_dim1 + 1;
    v6 -= v6_offset;
    --v5;
    --v4;
    --v3;
    --v2;
    --v1;
    --kbk;
    --kfwd;
    --lkount;
    p_dim1 = *mrowp;
    p_offset = p_dim1 + 1;
    p -= p_offset;
    --lista;
    --listv;
    nadj_dim1 = *mrows;
    nadj_offset = nadj_dim1 + 1;
    nadj -= nadj_offset;
    nsim_dim1 = *mrows;
    nsim_offset = nsim_dim1 + 1;
    nsim -= nsim_offset;
    --lv;
    --ls;
    --lnext;

    /* Function Body */

    initl_(np, len, llen, lhead, lmin, lvail, lnew, 
           &lnext[1], &ls[1], &lv[1], &lkount[1]);
    inits_(mrows, ncols, ms, ns, &nsim[nsim_offset], 
           &nadj[nadj_offset]);
    initk_(np, klen, khead, &kfwd[1], &kbk[1]);

    if (*np > *ms) 
    {

       /* more than one simplex in triangulation */

       i__1 = *ms - 1;
	      facet1_(&i__1, mrowp, mp, np, &p[p_offset], 
               &listv[1], &v1[1], &v2[1], &v3[1], 
               &v4[1], &v5[1], &v6[v6_offset], ierr);
	      if (*ierr == 0) 
       {
	          i__1 = *ms - 1;
	          adds_(&i__1, &listv[1], mrows, ncols, ms, ns, 
                 &nsim[nsim_offset], ierr);
      	}
      	if (*ierr == 0) 
       {
	          i__1 = *ms - 1;
	          addl_(&listv[1], &listv[1], np, &c__1, ms, 
                 &i__1, &listv[1], len, 
	          	    lhead, lmin, lvail, lnew, &lnext[1], 
                &ls[1], &lv[1], &lkount[1], ierr);
      	}
    } 

    else 
    {

      /* exactly one simplex in triangulation */

     	onesim_(np, &listv[1], &lista[1]);
     	adds_(np, &listv[1], mrows, ncols, ms, ns, 
            &nsim[nsim_offset], ierr);
    }

    *lmpty = *lhead;

    return 0;
} /* init_ */






/* Subroutine */ int onesim_(n, listv, lista)
integer *n, *listv, *lista;
{

/* ***begin prologue onesim                                    */
/* ***date written 950115   (yymmdd)                           */
/* ***keywords  Delaunay tetrahedrization tessellation         */
/* ***author Hazlewood, Carol                                  */
/*           ch04@academia.swt.edu                             */
/*           Department of Computer Science                    */
/*           Southwest Texas State University                  */
/*           San Marcos, TX 78666                              */
/*                                                             */
/* ***purpose  this routine sets up vertices and adjacency */
/*            information for a tessellation with exactly */
/*            one simplex */
/* %%% */
/* ***end prologue onesim */


    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i;



    /* Parameter adjustments */
    --lista;
    --listv;

    /* Function Body */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
       	listv[i] = i;
	       lista[i] = 0;
/* L1: */
    }
    return 0;
} /* onesim_ */






/* Subroutine */ int facet1_(md, mrowp, mp, np, p, listv, xn, e, v, w, rk, a, 
	ierr)
integer *md, *mrowp, *mp, *np;
real *p;
integer *listv;
real *xn, *e, *v, *w, *rk, *a;
integer *ierr;
{
    /* System generated locals */
    integer p_dim1, p_offset, a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    extern doublereal rhof_();
    static integer lmax;
    static real rmax;
    extern integer kont_();
    static integer i, j;
    static real r;
    extern /* Subroutine */ int scopy_(), saxpy_(), vsort_();
    extern integer minfnd_();
    extern /* Subroutine */ int getvec_(), brktie_();
    static integer loc;


/* this subroutine giftwraps a single facet in the boundary */
/* of the convex hull of the points in p, starting with a */
/* point of minimum first coordinate */



    /* Parameter adjustments */
    --listv;
    a_dim1 = *mp;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --rk;
    --w;
    --v;
    --e;
    --xn;
    p_dim1 = *mrowp;
    p_offset = p_dim1 + 1;
    p -= p_offset;

    /* Function Body */
    loc = minfnd_(mrowp, np, &p[p_offset]);
    i__1 = *md;
    for (i = 1; i <= i__1; ++i) {
	listv[i] = 0;
/* L1: */
    }
    listv[1] = loc;

    if (*md >= 2) {
	i__1 = *md;
	for (i = 2; i <= i__1; ++i) {
	    getvec_(&i, mrowp, mp, np, &p[p_offset], &listv[1], &loc, &rmax, &
		    xn[1], &e[1], &w[1], &rk[1], &a[a_offset], ierr);
	    rmax = (float)-1e32;
	    i__2 = *np;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = i - 1;
		if (kont_(&i__3, &listv[1], &j) == 0) {
		    scopy_(mp, &p[j * p_dim1 + 1], &c__1, &v[1], &c__1);
		    saxpy_(mp, &c_b100, &p[loc * p_dim1 + 1], &c__1, &v[1], &
			    c__1);
		    r = rhof_(mp, &v[1], &e[1], &xn[1]);
/* 	WRITE(*,*)'  J = ',J,'  RHO(J) = ',R */
		    if (r > rmax + TOL) {
			rmax = r;
			lmax = j;
		    } else if (r > rmax - TOL) {
			i__3 = i - 1;
			brktie_(&rmax, &lmax, &r, &j, mrowp, mp, np, &p[
				p_offset], &i__3, &listv[1]);
		    }
		}
/* L2: */
	    }
	    listv[i] = lmax;
/*      WRITE(*,*)'VERTEX = ',LMAX */
/* L3: */
	}
    }
/* 	WRITE(*,*)'  MD = ',MD */
/*      WRITE(*,*)'FACET = ',(LISTV(ICK),ICK=1,MD) */
    vsort_(md, &listv[1]);
/* C      WRITE(*,*)'SORTED FACET = ',(LISTV(ICK),ICK=1,MD) */
    return 0;
} /* facet1_ */






integer minfnd_(mrowp, np, p)
integer *mrowp, *np;
real *p;
{
    /* System generated locals */
    integer p_dim1, p_offset, ret_val, i__1;

    /* Local variables */
    static integer j, loc;




    /* Parameter adjustments */
    p_dim1 = *mrowp;
    p_offset = p_dim1 + 1;
    p -= p_offset;

    /* Function Body */
    loc = 1;
    i__1 = *np;
    for (j = 2; j <= i__1; ++j) {
	if (p[j * p_dim1 + 1] < p[loc * p_dim1 + 1]) {
	    loc = j;
	}
/* L1: */
    }
    ret_val = loc;
    return ret_val;
} /* minfnd_ */






integer kont_(n, list, k)
integer *n, *list, *k;
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer i;


/* linear search */
/* returns 0 if k is not in the first n elements of list */
/* else it returns the position of the first occurrence of */
/* k in the list. */



    /* Parameter adjustments */
    --list;

    /* Function Body */
    i = 1;
L1:
    if (i <= *n && list[i] != *k) {
	++i;
	goto L1;
    }
    if (i > *n) {
	i = 0;
    }
    ret_val = i;
    return ret_val;
} /* kont_ */






doublereal rhof_(mp, v, e, xn)
integer *mp;
real *v, *e, *xn;
{
    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Builtin functions */
    /* double r_sign();*/

    /* Local variables */
    extern doublereal sdot_();
    static real d1, d2;


/* giftwrapping primitive */



    /* Parameter adjustments */
    --xn;
    --e;
    --v;

    /* Function Body */
    d1 = sdot_(mp, &xn[1], &c__1, &v[1], &c__1);
    d2 = -(doublereal)sdot_(mp, &e[1], &c__1, &v[1], &c__1);
    if (d1 != (float)0.) {
	d2 /= d1;
/* Computing MIN */
	r__2 = dabs(d2);
	r__1 = dmin(r__2,(float)1e32);
        ret_val = (d2 >= 0.0) ? r__1 : -r__1;
    } else {
        ret_val = (d2 >= 0.0) ? c_b112 : c_b112;
    }
    return ret_val;
} /* rhof_ */






/* Subroutine */ int getvec_(n, mrowp, mp, np, p, listv, loc, xmax, xn, e, w, 
	rk, a, ierr)
integer *n, *mrowp, *mp, *np;
real *p;
integer *listv, *loc;
real *xmax, *xn, *e, *w, *rk, *a;
integer *ierr;
{
    /* System generated locals */
    integer p_dim1, p_offset, a_dim1, a_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2;

    /* Builtin functions */
    /*double r_sign();*/


    /* Local variables */
    extern doublereal sdot_(), snrm2_();
    static integer i, j, k;
    static real s, t;
    extern /* Subroutine */ int sscal_(), scopy_(), saxpy_(), hh_();
    static integer kk;


/* updates normal and tangent vectors for giftwrapping */



    /* Parameter adjustments */
    --listv;
    a_dim1 = *mp;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --rk;
    --w;
    --e;
    --xn;
    p_dim1 = *mrowp;
    p_offset = p_dim1 + 1;
    p -= p_offset;

    /* Function Body */
    if (*n <= 2) {

/* initialization */

	i__1 = *mp;
	for (k = 1; k <= i__1; ++k) {
	    xn[k] = (float)0.;
	    e[k] = (float)0.;
/* L1: */
	}
	xn[1] = (float)1.;
	e[2] = (float)1.;
	scopy_(mp, &e[1], &c__1, &a[*mp * a_dim1 + 1], &c__1);
    } else {

/* i is the number of vertices already found */
/* 2 .le. i .le mp-1 */

	i = *n - 1;

/* update xn and put in column i of a */

	if (*xmax == (float)0.) {

/* new normal is e.  swap e and xn. */

	    i__1 = *mp;
	    for (k = 1; k <= i__1; ++k) {
		t = e[k];
		e[k] = xn[k];
		xn[k] = e[k];
/* L2: */
	    }
	} else {
	    i__1 = *mp;
	    for (k = 1; k <= i__1; ++k) {
		xn[k] = *xmax * xn[k] + e[k];
/* L3: */
	    }
	    r__2 = sdot_(mp, &xn[1], &c__1, &e[1], &c__1);
            r__1 = ((r__2 >= 0.0) ? c_b121 : -c_b121)  
                   /snrm2_(mp, &xn[1], &c__1);
	    /*r__1 = r_sign(&c_b121, &r__2) / snrm2_(mp, &xn[1], &c__1);*/
	    sscal_(mp, &r__1, &xn[1], &c__1);
	    scopy_(mp, &xn[1], &c__1, &a[i * a_dim1 + 1], &c__1);
	}
/* 	WRITE(*,*) 'XMAX = ',XMAX, ' XN = ',XN */

/* construct vi in column i-1 of a */

	k = listv[i];
	scopy_(mp, &p[k * p_dim1 + 1], &c__1, &a[(i - 1) * a_dim1 + 1], &c__1)
		;
	saxpy_(mp, &c_b100, &p[*loc * p_dim1 + 1], &c__1, &a[(i - 1) * a_dim1 
		+ 1], &c__1);

/* apply Householder trans 1,...,i-2 to col i-1 and i */

	if (i > 2) {
	    i__1 = i;
	    for (j = i - 1; j <= i__1; ++j) {
		i__2 = i - 2;
		for (k = 1; k <= i__2; ++k) {
		    i__3 = *mp - k + 1;
		    i__4 = *mp - k + 1;
		    r__1 = -(doublereal)w[k] * sdot_(&i__4, &a[k + k * a_dim1]
			    , &c__1, &a[k + j * a_dim1], &c__1);
		    saxpy_(&i__3, &r__1, &a[k + k * a_dim1], &c__1, &a[k + j *
			     a_dim1], &c__1);
/* L4: */
		}
/* L5: */
	    }
	}

/* compute Householder trans for cols i-1 and i, */
/* apply col i-1 to col mp */
/* apply col i to col mp after it is moved to e */
	i__1 = *mp - i + 2;
	hh_(&i__1, &c__2, &a[i - 1 + (i - 1) * a_dim1], &w[i - 1], &rk[i - 1],
		 ierr);
	i__1 = *mp - i + 2;
	i__2 = *mp - i + 2;
	r__1 = -(doublereal)w[i - 1] * sdot_(&i__2, &a[i - 1 + (i - 1) * 
		a_dim1], &c__1, &a[i - 1 + *mp * a_dim1], &c__1);
	saxpy_(&i__1, &r__1, &a[i - 1 + (i - 1) * a_dim1], &c__1, &a[i - 1 + *
		mp * a_dim1], &c__1);
	scopy_(mp, &a[*mp * a_dim1 + 1], &c__1, &e[1], &c__1);
	i__1 = *mp - i + 1;
	i__2 = *mp - i + 1;
	r__1 = -(doublereal)w[i] * sdot_(&i__2, &a[i + i * a_dim1], &c__1, &e[
		i], &c__1);
	saxpy_(&i__1, &r__1, &a[i + i * a_dim1], &c__1, &e[i], &c__1);

/* apply Householder transformations i,...,1 to */
/* (0,...,0,e(i+1),...,e(mp)) to obtain projection */
/* of 'old' e onto the orthogonal complement of */
/* the column space of the first i columns of a. */
/* See Stewart, p.241. */

/* 	WRITE(*,*)' ***E = ',E */
	i__1 = *mp - i;
	s = snrm2_(&i__1, &e[i + 1], &c__1);
/* 	WRITE(*,*) ' I = ',I, ' S = ',S,' E = ', (E(ICK),ICK=1,MP) */
	scopy_(&i, &c_b149, &c__0, &e[1], &c__1);
/* 	WRITE(*,*) 'E = ', (E(ICK),ICK=1,MP) */
	i__1 = i;
	for (kk = 1; kk <= i__1; ++kk) {
	    k = i - kk + 1;
	    i__2 = *mp - k + 1;
	    i__3 = *mp - k + 1;
	    r__1 = -(doublereal)w[k] * sdot_(&i__3, &a[k + k * a_dim1], &c__1,
		     &e[k], &c__1);
	    saxpy_(&i__2, &r__1, &a[k + k * a_dim1], &c__1, &e[k], &c__1);
/* L7: */
	}

	s = snrm2_(mp, &e[1], &c__1);
/* 	WRITE(*,*) ' S = ',S,' E = ', (E(ICK),ICK=1,MP) */
	sscal_(mp, &s, &e[1], &c__1);
/*      do 8 k = i+1,mp */
/*        e(k) = e(k)/s */
/*    8 continue */
/* 	WRITE(*,*)' ---E = ',E */
    }
/* 	WRITE(*,*)'  XN = ',(XN(ICK),ICK=1,MP) */
/* 	WRITE(*,*)'   E = ',(E(ICK),ICK=1,MP) */
/* 	WRITE(*,*)'norm(xn) = ',snrm2(mp,xn,1) */
/* 	WRITE(*,*)'norm(e) = ',snrm2(mp,e,1) */
/* 	WRITE(*,*)'XN DOT E = ',SDOT(MP,XN,1,E,1) */
/* 	DO 100 K = 2,I */
/* 	KK = LISTV(K) */
/* 	WRITE(*,*)' VI = ',KK */
/* 	IF (KK .NE. 0) THEN */
/* 	WRITE(*,*)' I = ',I,' E DOT VI = ',SDOT(MP,E,1,P(1,KK),1)- */
/*     $             SDOT(MP,E,1,P(1,LOC),1) */
/* 	WRITE(*,*)'XN DOT VI = ',SDOT(MP,XN,1,P(1,KK),1)- */
/*     $             SDOT(MP,XN,1,P(1,LOC),1) */
/* 	ENDIF */
/*  100   CONTINUE */
    return 0;
} /* getvec_ */






/* Subroutine */ int brktie_(rmax, lmax, r, j, mrowp, mp, np, p, n, listv)
real *rmax;
integer *lmax;
real *r;
integer *j, *mrowp, *mp, *np;
real *p;
integer *n, *listv;


{
    /* System generated locals */
    integer p_dim1, p_offset;

    integer i;


    /* Parameter adjustments */
    p_dim1 = *mrowp;
    p_offset = p_dim1 + 1;
    p -= p_offset;
    --listv;

    /* Function Body */
  return(0);
} /* brktie_ */

/* Subroutine */ int hh_(m, n, a, w, rk, ierr)
integer *m, *n;
real *a, *w, *rk;
integer *ierr;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    real r__1;

    /* Builtin functions */
    /* double r_sign();*/

    /* Local variables */
    extern doublereal sdot_(), snrm2_();
    static integer j, k;
    static real s;
    extern /* Subroutine */ int saxpy_();
    static integer mm;
    static real akk, wkk;


/* computer householder qr factorization of m by n matrix a */
/* and places result in a, w, and rk */



    /* Parameter adjustments */
    --rk;
    --w;
    a_dim1 = *m;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	mm = *m - k + 1;
	s = snrm2_(&mm, &a[k + k * a_dim1], &c__1);
	if (s == (float)0.) {
	    goto L3;
	}
	/* akk = -(doublereal)r_sign(&s, &a[k + k * a_dim1]);*/
        akk = -(doublereal)(a[k + k * a_dim1] >= 0.0) ? s : -s;
	rk[k] = akk;
	wkk = a[k + k * a_dim1] - akk;
	w[k] = (float)-1. / (wkk * akk);
	a[k + k * a_dim1] = wkk;
	if (k < *n) {
	    i__2 = *n;
	    for (j = k + 1; j <= i__2; ++j) {
		r__1 = -(doublereal)w[k] * sdot_(&mm, &a[k + k * a_dim1], &
			c__1, &a[k + j * a_dim1], &c__1);
		saxpy_(&mm, &r__1, &a[k + k * a_dim1], &c__1, &a[k + j * 
			a_dim1], &c__1);
/* L1: */
	    }
	}
/* L2: */
    }
    goto L4;

L3:
    *ierr = 301;
/*    4  	WRITE(*,*)' A = ' */
/* 	DO 5 ICK=1,M */
/* 	  WRITE(*,*)(A(ICK,JCK),JCK=1,N) */
/*    5	CONTINUE */
L4:
    return 0;
} /* hh_ */






/* Subroutine */ int qv_(m, v, a, w)
integer *m;
real *v, *a, *w;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    real r__1;

    /* Local variables */
    extern doublereal sdot_();
    static integer k;
    extern /* Subroutine */ int scopy_(), saxpy_();
    static integer kk, mm;


/* given householder qr factorization of a matrix in a and */
/* w, this routine computes q*e(m) and returns the result */
/* in v. */



    /* Parameter adjustments */
    --w;
    a_dim1 = *m;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --v;

    /* Function Body */
    scopy_(m, &c_b149, &c__0, &v[1], &c__1);
    v[*m] = (float)1.;
    i__1 = *m - 1;
    for (kk = 1; kk <= i__1; ++kk) {
	k = *m - kk;
	mm = *m - k + 1;
	r__1 = -(doublereal)w[k] * sdot_(&mm, &a[k + k * a_dim1], &c__1, &v[k]
		, &c__1);
	saxpy_(&mm, &r__1, &a[k + k * a_dim1], &c__1, &v[k], &c__1);
/* L1: */
    }
/* 	WRITE(*,*)' V = ', (V(Ick),ICK=1,M) */
    return 0;
} /* qv_ */






doublereal opp_(mrowp, mp, np, p, listv, ms, k1, k2, init, p1tv, d1, v, w, rk,
	 a, ierr)
integer *mrowp, *mp, *np;
real *p;
integer *listv, *ms, *k1, *k2, *init;
real *p1tv, *d1, *v, *w, *rk, *a;
integer *ierr;
{
    /* System generated locals */
    integer p_dim1, p_offset, a_dim1, a_offset, i__1;
    real ret_val;

    /* Local variables */
    extern doublereal sdot_();
    static integer i, i1;
    extern /* Subroutine */ int scopy_(), saxpy_(), hh_();
    static integer ii;
    extern /* Subroutine */ int qv_();


/* if opp .gt. 0 then p(k1) and p(k2) are on the same side */
/*                    of the hyperplane that contains the */
/*                    points of p indexed by listv */
/*        .eq. 0 then p(k1) and p(k2) are in the hyperplane */
/*        .lt. 0 then p(k1) and p(k2) are on opposite sides */
/*                    of the hyperplane */



/* 	WRITE(*,*)' INIT = ',INIT,' IERR = ',IERR, */
/*     $	' MROWP = ',MROWP,' MP = ',MP,' NP = ',NP */
/* 	WRITE(*,*)' P1TV = ',P1TV,' D1 = ',D1 */
    /* Parameter adjustments */
    a_dim1 = *mp;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --rk;
    --w;
    --v;
    p_dim1 = *mrowp;
    p_offset = p_dim1 + 1;
    p -= p_offset;
    --listv;

    /* Function Body */
    if (*init == 0) {
	*init = 1;
	i1 = listv[1];
	i__1 = *ms - 2;
	for (i = 1; i <= i__1; ++i) {
	    ii = listv[i + 1];
	    scopy_(mp, &p[ii * p_dim1 + 1], &c__1, &a[i * a_dim1 + 1], &c__1);
	    saxpy_(mp, &c_b100, &p[i1 * p_dim1 + 1], &c__1, &a[i * a_dim1 + 1]
		    , &c__1);
/* L1: */
	}
/*    4  	WRITE(*,*)' A = ' */
/* 	DO 5 ICK=1,Mp */
/* 	  WRITE(*,*)(A(ICK,JCK),JCK=1,mp) */
/*    5	CONTINUE */
	i__1 = *ms - 2;
	hh_(mp, &i__1, &a[a_offset], &w[1], &rk[1], ierr);
	if (*ierr == 0) {
	    qv_(mp, &v[1], &a[a_offset], &w[1]);
	    *p1tv = sdot_(mp, &p[i1 * p_dim1 + 1], &c__1, &v[1], &c__1);
	    *d1 = sdot_(mp, &p[*k1 * p_dim1 + 1], &c__1, &v[1], &c__1) - *
		    p1tv;
	}
    }

    if (*ierr == 0) {
	ret_val = *d1 * (sdot_(mp, &p[*k2 * p_dim1 + 1], &c__1, &v[1], &c__1) 
		- *p1tv);
    }
/* 	WRITE(*,*)' INIT = ',INIT,' IERR = ',IERR, */
/*     $	' MROWP = ',MROWP,' MP = ',MP,' NP = ',NP */
/* 	WRITE(*,*)' P1TV = ',P1TV,' D1 = ',D1 */
    return ret_val;
} /* opp_ */






doublereal rhoa_(iap, mrowp, mp, np, p, ms, listv, init, c, x, v, w, rk, a, 
	ierr)
integer *iap, *mrowp, *mp, *np;
real *p;
integer *ms, *listv, *init;
real *c, *x, *v, *w, *rk, *a;
integer *ierr;
{
    /* System generated locals */
    integer p_dim1, p_offset, a_dim1, a_offset, i__1;
    real ret_val, r__1;

    /* Local variables */
    extern doublereal sdot_();
    static integer i, i1;
    extern /* Subroutine */ int scopy_(), saxpy_(), cc_(), hh_();
    static integer ii;
    extern /* Subroutine */ int qv_();




/*      WRITE(*,*)'LISTV = ',(LISTV(ICK),ICK=1,MS) */
    /* Parameter adjustments */
    a_dim1 = *mp;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --rk;
    --w;
    --v;
    --x;
    --c;
    p_dim1 = *mrowp;
    p_offset = p_dim1 + 1;
    p -= p_offset;
    --listv;

    /* Function Body */
    if (*init < 2) {
	if (*init < 1) {
	    i1 = listv[1];
	    i__1 = *ms - 2;
	    for (i = 1; i <= i__1; ++i) {
		ii = listv[i + 1];
		scopy_(mp, &p[ii * p_dim1 + 1], &c__1, &a[i * a_dim1 + 1], &
			c__1);
		saxpy_(mp, &c_b100, &p[i1 * p_dim1 + 1], &c__1, &a[i * a_dim1 
			+ 1], &c__1);
/* L1: */
	    }
/*    4  	WRITE(*,*)' A = ' */
/* 	DO 5 ICK=1,M */
/* 	  WRITE(*,*)(A(ICK,JCK),JCK=1,N) */
/*    5	CONTINUE */
	    i__1 = *ms - 2;
	    hh_(mp, &i__1, &a[a_offset], &w[1], &rk[1], ierr);
	    if (*ierr == 0) {
		qv_(mp, &v[1], &a[a_offset], &w[1]);
	    }
	}
	*init = 2;
	cc_(mp, &c[1], &w[1], &rk[1], &a[a_offset]);
    }

    i1 = listv[1];
    scopy_(mp, &p[*iap * p_dim1 + 1], &c__1, &x[1], &c__1);
    saxpy_(mp, &c_b100, &p[i1 * p_dim1 + 1], &c__1, &x[1], &c__1);
    ret_val = (sdot_(mp, &x[1], &c__1, &x[1], &c__1) - sdot_(mp, &c[1], &c__1,
	     &x[1], &c__1) * (float)2.) / ((r__1 = sdot_(mp, &v[1], &c__1, &x[
	    1], &c__1), dabs(r__1)) * (float)2.);
    return ret_val;
} /* rhoa_ */






/* Subroutine */ int cc_(mp, c, w, rk, a)
integer *mp;
real *c, *w, *rk, *a;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    extern doublereal sdot_();
    static integer i, k;
    extern /* Subroutine */ int saxpy_();
    static integer kk, mm;




    /* Parameter adjustments */
    a_dim1 = *mp;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --rk;
    --w;
    --c;

    /* Function Body */
    c[1] = rk[1] / (float)2.;
    c[*mp] = (float)0.;
    if (*mp > 2) {
	i__1 = *mp - 1;
	for (k = 2; k <= i__1; ++k) {
	    c[k] = (float)0.;
	    i__2 = k - 1;
	    for (i = 1; i <= i__2; ++i) {
		c[k] += a[i + k * a_dim1] * (a[i + k * a_dim1] - c[i] * (
			float)2.);
/* L1: */
	    }
	    c[k] = (rk[k] + c[k] / rk[k]) / (float)2.;
/* L2: */
	}
    }
    i__1 = *mp - 1;
    for (kk = 1; kk <= i__1; ++kk) {
	k = *mp - kk;
	mm = *mp - k + 1;
	r__1 = -(doublereal)w[k] * sdot_(&mm, &a[k + k * a_dim1], &c__1, &c[k]
		, &c__1);
	saxpy_(&mm, &r__1, &a[k + k * a_dim1], &c__1, &c[k], &c__1);
/* L3: */
    }
    return 0;
} /* cc_ */


/* operations on data structures */

/* Subroutine */ int inits_(mrows, ncols, ms, ns, nsim, nadj)
integer *mrows, *ncols, *ms, *ns, *nsim, *nadj;
{
    /* System generated locals */
    integer nsim_dim1, nsim_offset, nadj_dim1, nadj_offset, i__1, i__2;

    /* Local variables */
    static integer i, j;




/* 	WRITE(*,*)'MROWS = ',MROWS,'   NCOLS = ',NCOLS */
    /* Parameter adjustments */
    nadj_dim1 = *mrows;
    nadj_offset = nadj_dim1 + 1;
    nadj -= nadj_offset;
    nsim_dim1 = *mrows;
    nsim_offset = nsim_dim1 + 1;
    nsim -= nsim_offset;

    /* Function Body */
    i__1 = *ncols;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *mrows;
	for (i = 1; i <= i__2; ++i) {
	    nsim[i + j * nsim_dim1] = 0;
	    nadj[i + j * nadj_dim1] = 0;
/* L1: */
	}
/* L2: */
    }
    *ns = -1;
    return 0;
} /* inits_ */






/* Subroutine */ int adds_(nv, listv, mrows, ncols, ms, ns, nsim, ierr)
integer *nv, *listv, *mrows, *ncols, *ms, *ns, *nsim, *ierr;
{
    /* System generated locals */
    integer nsim_dim1, nsim_offset, i__1;

    /* Local variables */
    static integer i, j;



/* inserts listv(vertices) and lista(adjacency info) into */
/* list of simplices */


    /* Parameter adjustments */
    --listv;
    nsim_dim1 = *mrows;
    nsim_offset = nsim_dim1 + 1;
    nsim -= nsim_offset;

    /* Function Body */
    if (*ns < *ncols) {
	if (*ns == -1) {
	    *ns = 0;
	    j = 1;
	} else {
	    ++(*ns);
	    j = *ns;
	}
	i__1 = *nv;
	for (i = 1; i <= i__1; ++i) {
	    nsim[i + j * nsim_dim1] = listv[i];
/* L1: */
	}
    } else {
	*ierr = 111;
    }
/* 	WRITE(*,*)'ADDS:  NS =',NS,'  ',(NSIM(ICK,J),ICK=1,MS) */
    return 0;
} /* adds_ */






/* Subroutine */ int vsort_(nv, listv)
integer *nv, *listv;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer lmin, i, j, min_;



/* sorts first nv items of listv in increasing order */
/* permutes lista to preserve original correspondence */


    /* Parameter adjustments */
    --listv;

    /* Function Body */
    if (*nv >= 2) {
	i__1 = *nv - 1;
	for (j = 1; j <= i__1; ++j) {
	    min_ = listv[j];
	    lmin = j;
	    i__2 = *nv;
	    for (i = j + 1; i <= i__2; ++i) {
		if (listv[i] < min_) {
		    min_ = listv[i];
		    lmin = i;
		}
/* L1: */
	    }
	    listv[lmin] = listv[j];
	    listv[j] = min_;
/* L2: */
	}
    }
    return 0;
} /* vsort_ */






/* Subroutine */ int initl_(np, len, llen, lhead, lmin, lvail, lnew, lnext, 
	ls, lv, lkount)
integer *np, *len, *llen, *lhead, *lmin, *lvail, *lnew, *lnext, *ls, *lv, *
	lkount;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i;



/* this routine initializes a list of open facets */


    /* Parameter adjustments */
    --lkount;
    --lv;
    --ls;
    --lnext;

    /* Function Body */
    i__1 = *len;
    for (i = 1; i <= i__1; ++i) {
	lnext[i] = 0;
	ls[i] = 0;
	lv[i] = 0;
/* L1: */
    }
    i__1 = *np;
    for (i = 1; i <= i__1; ++i) {
	lkount[i] = -1;
/* L2: */
    }
    *llen = *len;
    *lhead = 0;
    *lmin = *np + 1;
    *lvail = 0;
    *lnew = *np + 1;
    return 0;
} /* initl_ */






/* Subroutine */ int addl_(loc, minv, np, ns, mv, nv, listv, len, lhead, lmin,
	 lvail, lnew, lnext, ls, lv, lkount, ierr)
integer *loc, *minv, *np, *ns, *mv, *nv, *listv, *len, *lhead, *lmin, *lvail, 
	*lnew, *lnext, *ls, *lv, *lkount, *ierr;
{
    /* System generated locals */
    integer i__1;


    /* Local variables */
    static integer i, kv;




/* 	WRITE(*,*)' ENTER ADDL: LMIN = ',LMIN,'  LHEAD = ',LHEAD */
/* 	WRITE(*,*)' LOC = ',LOC,' MV = ',MV */
/* 	WRITE(*,*)'LISTV = ', (LISTV(ICK), ICK=1,NV) */
    /* Parameter adjustments */
    --lkount;
    --listv;
    --lv;
    --ls;
    --lnext;

    /* Function Body */
    if (*loc == 0) {
	if (*mv == 1) {
	    *loc = listv[2];
	} else {
	    *loc = listv[1];
	}
    }
    if (*lhead == 0 || *lhead == *np + 1) {
	*lhead = *loc;
    }

    if (ls[*loc] != 0) {
	if (*lvail != 0) {
	    lnext[*loc] = *lvail;
	    *lvail = lnext[*lvail];
	} else {
	    if (*lnew <= *len) {
		lnext[*loc] = *lnew;
		++(*lnew);
	    } else {
		*ierr = 101;
		return(0);
	    }
	}
	*loc = lnext[*loc];
    }

/* 	WRITE(*,*) 'LOC = ',LOC, ' IERR = ',IERR */
    if (*ierr == 0) {
	ls[*loc] = *ns;
	lv[*loc] = *mv;
	lnext[*loc] = 0;
	if (*minv < *lmin) {
	    *lmin = *minv;
	}

	i__1 = *nv;
	for (i = 1; i <= i__1; ++i) {
	    if (i != *mv) {
		kv = listv[i];
		if (lkount[kv] == -1) {
		    lkount[kv] = 1;
/* 	WRITE(*,*) 'KV = ',KV, ' LKOUNT(KV) = ',LKOUNT(KV) */
		} else if (lkount[kv] >= 0) {
		    ++lkount[kv];
/* 	WRITE(*,*) 'KV = ',KV, ' LKOUNT(KV) = ',LKOUNT(KV) */
		} else {
		    *ierr = 102;
/* 	WRITE(*,*) 'KV = ',KV, ' LKOUNT(KV) = ',LKOUNT(KV) */
		}
	    }
/* L1: */
	}
    }
/* 	WRITE(*,*)' EXIT ADDL: LMIN = ',LMIN,'  LHEAD = ',LHEAD */
    return 0;
} /* addl_ */






/* Subroutine */ int reml_(loc, lprev, np, ms, listv, mv, len, lhead, lmin, 
	lvail, lnext, ls, lv, lkount, ierr)
integer *loc, *lprev, *np, *ms, *listv, *mv, *len, *lhead, *lmin, *lvail, *
	lnext, *ls, *lv, *lkount, *ierr;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, k;




/* remove node at loc from list */
/* put on availlist unless head node */

/* 	WRITE(*,*)' LOC = ',LOC */
/* 	WRITE(*,*)' LS(LOC) = ',LS(LOC),'LV(LOC) = ',LV(LOC) */
/* 	WRITE(*,*)' ENTER REML: LMIN = ',LMIN,'  LHEAD = ',LHEAD */
    /* Parameter adjustments */
    --lkount;
    --listv;
    --lv;
    --ls;
    --lnext;

    /* Function Body */
    ls[*loc] = 0;
    lv[*loc] = 0;
    if (*lprev != 0) {
	lnext[*lprev] = lnext[*loc];
	if (*lvail != 0) {
	    lnext[*loc] = lnext[*lvail];
	} else {
	    lnext[*loc] = 0;
	}
	*lvail = *loc;
    }

/* update lmin and lhead */

    if (*mv == 1) {
	k = listv[2];
    } else {
	k = listv[1];
    }
    if (lnext[k] == 0 && ls[k] == 0) {
	if (*lmin == k) {
	    i__1 = *np;
	    for (i = 1; i <= i__1; ++i) {
		if (lnext[i] != 0 || ls[i] != 0) {
		    *lmin = i;
		    goto L2;
		}
/* L1: */
	    }
	    *lmin = *np + 1;
	}
L2:
	if (*lhead == k) {
	    *lhead = *lmin;
	}
    }

/* update kount */

    i__1 = *ms;
    for (i = 1; i <= i__1; ++i) {
	if (i != *mv) {
	    k = listv[i];
	    if (lkount[k] > 0) {
		--lkount[k];
/* 	WRITE(*,*)' REML:  I = ',I,' LISTV(I) = ',K, */
/*     $            ' LKNT(K) = ',LKOUNT(K) */
	    } else {
		*ierr = 203;
	    }
	}
/* L3: */
    }
/* 	WRITE(*,*)' EXIT REML: LMIN = ',LMIN,'  LHEAD = ',LHEAD */
    return 0;
} /* reml_ */






/* Subroutine */ int lookup_(loc, lprev, mv, listv, list, len, lnext, ls, lv, 
	mrows, ncols, ms, ns, nsim)
integer *loc, *lprev, *mv, *listv, *list, *len, *lnext, *ls, *lv, *mrows, *
	ncols, *ms, *ns, *nsim;
{
    /* System generated locals */
    integer nsim_dim1, nsim_offset, i__1;

    /* Local variables */
    static integer ifnd, lsav, i, k;
    extern /* Subroutine */ int match_(), svtol_(), insert_();




/* 	WRITE(*,*)' LISTV1 = ',(LISTV(ICK),ICK=1,MS) */
/* 	WRITE(*,*)' MV = ',MV */
    /* Parameter adjustments */
    --lv;
    --ls;
    --lnext;
    nsim_dim1 = *mrows;
    nsim_offset = nsim_dim1 + 1;
    nsim -= nsim_offset;
    --list;
    --listv;

    /* Function Body */
    ifnd = 0;
    *lprev = 0;
    lsav = listv[*mv];
    if (*mv < *ms) {
	i__1 = *ms - 1;
	for (i = *mv; i <= i__1; ++i) {
	    listv[i] = listv[i + 1];
/* L1: */
	}
    }
/* 	WRITE(*,*)' LISTV2 = ',(LISTV(ICK),ICK=1,MS) */
    k = listv[1];
    if (ls[k] == 0) {
	*lprev = k;
	k = lnext[k];
    }
L2:
    if (k != 0 && ifnd == 0) {
/* 	WRITE(*,*)' LISTV3 = ',(LISTV(ICK),ICK=1,MS) */
/* 	WRITE(*,*)' K = ',K,' LS(K) = ',LS(K),' LV(K) = ',LV(K) */
/* 	WRITE(*,*)'NSIM = ',(NSIM(ICK,LV(K)),ICK=1,MS) */
	svtol_(&ls[k], &lv[k], &list[1], mrows, ncols, ms, ns, &nsim[
		nsim_offset]);
	i__1 = *ms - 1;
	match_(&i__1, &listv[1], &list[1], &ifnd);
/* 	WRITE(*,*) 'IFND = ',IFND */
	if (ifnd == 0) {
	    *lprev = k;
	    k = lnext[k];
	}
	goto L2;
    }
/* 	WRITE(*,*)' LISTV4 = ',(LISTV(ICK),ICK=1,MS) */
    *loc = k;
    if (*mv < *ms) {
	insert_(&lsav, ms, &listv[1]);
    }
/* 	WRITE(*,*)' LISTV5 = ',(LISTV(ICK),ICK=1,MS) */
/* 	WRITE(*,*)'LOOKUP: LOC = ',LOC,' LPREV = ',LPREV */
/* 	IF (LOC .NE. 0) THEN */
/*        WRITE(*,*)'LS(LOC) = ',LS(LOC), */
/*     $  'LV(LOC) = ',LV(LOC) */
/*        WRITE(*,*)'LISTV = ',(LISTV(ICK),ICK=1,MS) */
/* 	JCK=LS(LOC) */
/* 	WRITE(*,*)'NSIM(,) = ',(NSIM(ICK,JCK),ICK=1,MS) */
/*        ENDIF */
    return 0;
} /* lookup_ */




/* Subroutine */ int match_(n, list1, list2, ifnd)
integer *n, *list1, *list2, *ifnd;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i;



/* if list1 = list2 ifnd is returned with value 1 */
/* else ifnd is unchanged */


    /* Parameter adjustments */
    --list2;
    --list1;

    /* Function Body */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	if (list1[i] != list2[i]) {
	    goto L2;
	}
/* L1: */
    }
    *ifnd = 1;
/*        WRITE(*,*)(LIST1(ICK),ICK=1,N) */
/* 	WRITE(*,*)(LIST2(ICK),ICK=1,N) */
/* 	WRITE(*,*)' IFND = ',IFND */
L2:
    return 0;
} /* match_ */




/* Subroutine */ int svtol_(nt, mv, listv, mrows, ncols, ms, ns, nsim)
integer *nt, *mv, *listv, *mrows, *ncols, *ms, *ns, *nsim;
{
    /* System generated locals */
    integer nsim_dim1, nsim_offset, i__1;

    /* Local variables */
    static integer i, k;




/*      WRITE(*,*)'SVTOL: NT = ',NT,' MV = ',MV,'MS = ',MS */
    /* Parameter adjustments */
    --listv;
    nsim_dim1 = *mrows;
    nsim_offset = nsim_dim1 + 1;
    nsim -= nsim_offset;

    /* Function Body */
    if (*nt > 0) {
	i__1 = *ms;
	for (i = 1; i <= i__1; ++i) {
	    listv[i] = 0;
/* L1: */
	}
	k = 1;
	i__1 = *ms;
	for (i = 1; i <= i__1; ++i) {
	    if (i != *mv) {
		listv[k] = nsim[i + *nt * nsim_dim1];
		++k;
	    }
/* L2: */
	}
    }
/* 	WRITE(*,*)'LISTV = ',(LISTV(ICK),ICK=1,MS) */
    return 0;
} /* svtol_ */






/* Subroutine */ int initk_(n, klen, khead, kfwd, kbk)
integer *n, *klen, *khead, *kfwd, *kbk;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i;




    /* Parameter adjustments */
    --kbk;
    --kfwd;

    /* Function Body */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	kfwd[i] = i + 1;
	kbk[i] = i - 1;
/* L1: */
    }
    kfwd[*n] = 0;
    *klen = *n;
    *khead = 1;
    return 0;
} /* initk_ */






/* Subroutine */ int delk_(node, klen, khead, kfwd, kbk)
integer *node, *klen, *khead, *kfwd, *kbk;
{
    static integer inxt, iprev;



/* deletes node from candidate list */


    /* Parameter adjustments */
    --kbk;
    --kfwd;

    /* Function Body */
    iprev = kbk[*node];
    if (iprev != 0) {
	kfwd[iprev] = kfwd[*node];
    }
    inxt = kfwd[*node];
    if (inxt != 0) {
	kbk[inxt] = kbk[*node];
    }
    if (*khead == *node) {
	*khead = kfwd[*node];
    }
    kfwd[*node] = 0;
    kbk[*node] = 0;
    return 0;
} /* delk_ */

/* Subroutine */ int nxtsim_(mrowp, mp, np, p, mrows, ncols, ms, ns, nsim, 
	nadj, listv, lista, len, lhead, lmin, lnew, lvail, lnext, ls, lv, 
	lkount, klen, khead, kfwd, kbk, v1, v2, v3, v4, v5, v6, lmpty, ierr)
integer *mrowp, *mp, *np;
real *p;
integer *mrows, *ncols, *ms, *ns, *nsim, *nadj, *listv, *lista, *len, *lhead, 
	*lmin, *lnew, *lvail, *lnext, *ls, *lv, *lkount, *klen, *khead, *kfwd,
	 *kbk;
real *v1, *v2, *v3, *v4, *v5, *v6;
integer *lmpty, *ierr;
{
    /* System generated locals */
    integer nsim_dim1, nsim_offset, nadj_dim1, nadj_offset, p_dim1, p_offset, 
	    v6_dim1, v6_offset;

    /* Local variables */
    extern /* Subroutine */ int apex_();
    static integer mv;
    extern /* Subroutine */ int update_(), facnxt_();
    static integer iap;




/* 	WRITE(*,*)'ENTERING NXTSIM' */
    /* Parameter adjustments */
    v6_dim1 = *mp;
    v6_offset = v6_dim1 + 1;
    v6 -= v6_offset;
    --v5;
    --v4;
    --v3;
    --v2;
    --v1;
    --kbk;
    --kfwd;
    --lkount;
    p_dim1 = *mrowp;
    p_offset = p_dim1 + 1;
    p -= p_offset;
    --lista;
    --listv;
    nadj_dim1 = *mrows;
    nadj_offset = nadj_dim1 + 1;
    nadj -= nadj_offset;
    nsim_dim1 = *mrows;
    nsim_offset = nsim_dim1 + 1;
    nsim -= nsim_offset;
    --lv;
    --ls;
    --lnext;

    /* Function Body */
    facnxt_(mrows, ncols, ms, ns, &nsim[nsim_offset], &nadj[nadj_offset], len,
	     &listv[1], &mv, lhead, &lnext[1], &ls[1], &lv[1]);
/* C	WRITE(*,*) 'NT = ','  MV = ',MV */
/* C	WRITE(*,*)' LISTV11 = ',(LISTV(ICK),ICK=1,MS) */
/* 	WRITE(*,*)'ENTERING APEX' */
    apex_(&iap, mrowp, mp, np, &p[p_offset], &listv[1], ms, &mv, klen, khead, 
	    &kfwd[1], &kbk[1], &v1[1], &v2[1], &v3[1], &v4[1], &v5[1], &v6[
	    v6_offset], ierr);
/* 	WRITE(*,*)' LISTV12 = ',(LISTV(ICK),ICK=1,MS) */
    update_(&iap, &listv[1], &lista[1], np, mrows, ncols, ms, ns, &nsim[
	    nsim_offset], &nadj[nadj_offset], len, lhead, lmin, lvail, lnew, &
	    lnext[1], &ls[1], &lv[1], &lkount[1], klen, khead, &kfwd[1], &kbk[
	    1], ierr);
    if (*lhead == *np + 1) {
	*lmpty = 0;
    }
    return 0;
} /* nxtsim_ */






/* Subroutine */ int facnxt_(mrows, ncols, ms, ns, nsim, nadj, len, listv, mv,
	 lhead, lnext, ls, lv)
integer *mrows, *ncols, *ms, *ns, *nsim, *nadj, *len, *listv, *mv, *lhead, *
	lnext, *ls, *lv;
{
    /* System generated locals */
    integer nsim_dim1, nsim_offset, nadj_dim1, nadj_offset;

    /* Local variables */
    extern /* Subroutine */ int svtol_();
    static integer lh, nt;




/*      WRITE(*,*)'FACNXT:  LHEAD = ',LHEAD,lnext(lhead) */
    /* Parameter adjustments */
    --listv;
    nadj_dim1 = *mrows;
    nadj_offset = nadj_dim1 + 1;
    nadj -= nadj_offset;
    nsim_dim1 = *mrows;
    nsim_offset = nsim_dim1 + 1;
    nsim -= nsim_offset;
    --lv;
    --ls;
    --lnext;

    /* Function Body */
    lh = *lhead;
    if (ls[lh] == 0) {
	lh = lnext[lh];
    }
    nt = ls[lh];
    *mv = lv[lh];
/*      	WRITE(*,*)'NSIM = ',(NSIM(ICK,NT),ICK=1,MS) */
    svtol_(&nt, mv, &listv[1], mrows, ncols, ms, ns, &nsim[nsim_offset]);
/* 	WRITE(*,*)'LISTV = ',(LISTV(ICK),ICK=1,MS),' MV = ',MV */
    *mv = nsim[*mv + nt * nsim_dim1];
    return 0;
} /* facnxt_ */






/* Subroutine */ int apex_(iap, mrowp, mp, np, p, listv, ms, mv, klen, khead, 
	kfwd, kbk, v1, v2, v3, v4, v5, v6, ierr)
integer *iap, *mrowp, *mp, *np;
real *p;
integer *listv, *ms, *mv, *klen, *khead, *kfwd, *kbk;
real *v1, *v2, *v3, *v4, *v5, *v6;
integer *ierr;
{
    /* System generated locals */
    integer p_dim1, p_offset, v6_dim1, v6_offset;

    /* Local variables */
    extern doublereal rhoa_();
    static integer init;
    static real rmin, r;
    extern /* Subroutine */ int cosph_();
    static real d1;
    extern /* Subroutine */ int kandnx_();
    static integer min_;
    static real p1tv;




    /* Parameter adjustments */
    v6_dim1 = *mp;
    v6_offset = v6_dim1 + 1;
    v6 -= v6_offset;
    --v5;
    --v4;
    --v3;
    --v2;
    --v1;
    p_dim1 = *mrowp;
    p_offset = p_dim1 + 1;
    p -= p_offset;
    --listv;
    --kbk;
    --kfwd;

    /* Function Body */
    init = 0;
    *iap = 0;
    rmin = (float)1e32;
    min_ = 0;
    kandnx_(mrowp, mp, np, &p[p_offset], ms, &listv[1], mv, iap, &init, klen, 
	    khead, &kfwd[1], &v3[1], &v4[1], &v5[1], &v6[v6_offset], &p1tv, &
	    d1, ierr);
L1:
    if (*iap != 0 && *ierr == 0) {
	r = rhoa_(iap, mrowp, mp, np, &p[p_offset], ms, &listv[1], &init, &v1[
		1], &v2[1], &v3[1], &v4[1], &v5[1], &v6[v6_offset], ierr);
/* C	WRITE(*,*)'IAP = ',IAP,' RHO(IAP) = ',R */
	if (r < rmin - TOL) {
	    rmin = r;
	    min_ = *iap;
	} else {
	    if (r <= rmin + TOL) {
		cosph_(mrowp, mp, np, &p[p_offset]);
	    }
	}
	kandnx_(mrowp, mp, np, &p[p_offset], ms, &listv[1], mv, iap, &init, 
		klen, khead, &kfwd[1], &v3[1], &v4[1], &v5[1], &v6[v6_offset],
		 &p1tv, &d1, ierr);
	goto L1;
    }
    *iap = min_;
/*     listv(ms) = iap */
/*      if (min .eq. 0) ierr = 601 */
/* C	WRITE(*,*)'APEX:  IAP = ',IAP */
    return 0;
} /* apex_ */




/* Subroutine */ int kandnx_(mrowp, mp, np, p, ms, listv, mv, iap, init, klen,
	 khead, kfwd, v1, v2, v3, v6, p1tv, d1, ierr)
integer *mrowp, *mp, *np;
real *p;
integer *ms, *listv, *mv, *iap, *init, *klen, *khead, *kfwd;
real *v1, *v2, *v3, *v6, *p1tv, *d1;
integer *ierr;
{
    /* System generated locals */
    integer p_dim1, p_offset, v6_dim1, v6_offset, i__1;

    /* Local variables */
    extern integer kont_();
    extern doublereal opp_();
    static integer knx;




    /* Parameter adjustments */
    v6_dim1 = *mp;
    v6_offset = v6_dim1 + 1;
    v6 -= v6_offset;
    --v3;
    --v2;
    --v1;
    --kfwd;
    p_dim1 = *mrowp;
    p_offset = p_dim1 + 1;
    p -= p_offset;
    --listv;

    /* Function Body */
    if (*iap == 0) {
	knx = *khead;
    } else {
	knx = kfwd[*iap];
    }
L1:
    if (knx != 0) {
	i__1 = *ms - 1;
	if (kont_(&i__1, &listv[1], &knx) == 0 && knx != *mv) {
	    if (*mv == 0) {
		goto L2;
	    }
/* 	WRITE(*,*)'KANDNX:LISTV = ',(LISTV(ICK),ICK=1,MS),' MV = ',MV 
*/
/* 	WRITE(*,*)' KNX = ',KNX */
/*          o = opp(mrowp,mp,np,p,listv,ms,mv,knx,init,p1tv,d1, */
/*     $            v1,v2,v3,v6,ierr) */
/* 	WRITE(*,*) 'OPP(',MV,',',KNX,') = ',O */
	    if (opp_(mrowp, mp, np, &p[p_offset], &listv[1], ms, mv, &knx, 
		    init, p1tv, d1, &v1[1], &v2[1], &v3[1], &v6[v6_offset], 
		    ierr) < (float)0.) {
		goto L2;
	    }
	}
	knx = kfwd[knx];
	goto L1;
    }
L2:
    *iap = knx;
/*      WRITE(*,*)'FROM KANDNX:  IERR = ',IERR */
    return 0;
} /* kandnx_ */




/* Subroutine */ int cosph_(mrowp, mp, np, p)
integer *mrowp, *mp, *np;
real *p;
{
    /* System generated locals */
    integer p_dim1, p_offset;


    /* Parameter adjustments */
    p_dim1 = *mrowp;
    p_offset = p_dim1 + 1;
    p -= p_offset;

    /* Function Body */
    return 0;
} /* cosph_ */

/* Subroutine */ int update_(iap, listv, list, np, mrows, ncols, ms, ns, nsim,
	 nadj, len, lhead, lmin, lvail, lnew, lnext, ls, lv, kount, klen, 
	khead, kfwd, kbk, ierr)
integer *iap, *listv, *list, *np, *mrows, *ncols, *ms, *ns, *nsim, *nadj, *
	len, *lhead, *lmin, *lvail, *lnew, *lnext, *ls, *lv, *kount, *klen, *
	khead, *kfwd, *kbk, *ierr;
{
    /* System generated locals */
    integer nsim_dim1, nsim_offset, nadj_dim1, nadj_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int addl_(), adds_(), delk_(), reml_();
    static integer minv, i;
    extern /* Subroutine */ int fudge_();
    static integer lprev, kv, ks;
    extern /* Subroutine */ int insert_(), lookup_();
    static integer loc;




/* 	WRITE(*,*)'IAP = ',IAP, 'LISTV = ',(LISTV(ICK),ICK=1,MS) */
    /* Parameter adjustments */
    --kount;
    nadj_dim1 = *mrows;
    nadj_offset = nadj_dim1 + 1;
    nadj -= nadj_offset;
    nsim_dim1 = *mrows;
    nsim_offset = nsim_dim1 + 1;
    nsim -= nsim_offset;
    --list;
    --listv;
    --lv;
    --ls;
    --lnext;
    --kbk;
    --kfwd;

    /* Function Body */
    if (*iap != 0) {
	insert_(iap, ms, &listv[1]);
	if (*ns == 0) {
	    fudge_(iap, ms, &listv[1], len, &lnext[1], &ls[1], &lv[1], lhead);
	}
	adds_(ms, &listv[1], mrows, ncols, ms, ns, &nsim[nsim_offset], ierr);
	i__1 = *ms;
	for (i = 1; i <= i__1; ++i) {
/* 	WRITE(*,*)'LISTV = ',(LISTV(ICK),ICK=1,MS) */
	    lookup_(&loc, &lprev, &i, &listv[1], &list[1], len, &lnext[1], &
		    ls[1], &lv[1], mrows, ncols, ms, ns, &nsim[nsim_offset]);
/* 	WRITE(*,*)'LOC = ',LOC,' LPREV = ',LPREV, 'MV = ',I */
	    if (loc != 0) {
		kv = lv[loc];
		ks = ls[loc];
		if (ks > 0 && ks != *ns) {
		    nadj[kv + ks * nadj_dim1] = *ns;
		    nadj[i + *ns * nadj_dim1] = ks;
		} else {
		    nadj[i + *ns * nadj_dim1] = 0;
		}
		reml_(&loc, &lprev, np, ms, &listv[1], &i, len, lhead, lmin, 
			lvail, &lnext[1], &ls[1], &lv[1], &kount[1], ierr);
/* 	WRITE(*,*)'LISTV = ',(LISTV(ICK),ICK=1,MS) */
	    } else {
		if (i == 1) {
		    minv = listv[2];
		} else {
		    minv = listv[1];
		}
		addl_(&lprev, &minv, np, ns, &i, ms, &listv[1], len, lhead, 
			lmin, lvail, lnew, &lnext[1], &ls[1], &lv[1], &kount[
			1], ierr);
/* 	WRITE(*,*)'LISTV = ',(LISTV(ICK),ICK=1,MS) */
	    }
/* L1: */
	}
	i__1 = *ms;
	for (i = 1; i <= i__1; ++i) {
	    kv = listv[i];
	    if (kount[kv] == 0) {
		delk_(&kv, klen, khead, &kfwd[1], &kbk[1]);
	    }
/* L2: */
	}
    } else {
	lookup_(&loc, &lprev, ms, &listv[1], &list[1], len, &lnext[1], &ls[1],
		 &lv[1], mrows, ncols, ms, ns, &nsim[nsim_offset]);
/* 	WRITE(*,*)'LISTV = ',(LISTV(ICK),ICK=1,MS) */
	reml_(&loc, &lprev, np, ms, &listv[1], ms, len, lhead, lmin, lvail, &
		lnext[1], &ls[1], &lv[1], &kount[1], ierr);
    }
/* 	CALL WTESS(MROWS,NCOLS,MS,NS,NSIM,NADJ,IERR) */
/* 	CALL WRDS(LEN,LHEAD,LMIN,LVAIL,LNEW,LNEXT,LS,LV,KOUNT, */
/*     $            KLEN,KHEAD,KFWD,KBK) */
    return 0;
} /* update_ */




/* Subroutine */ int insert_(iap, ms, listv)
integer *iap, *ms, *listv;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, ii;




    /* Parameter adjustments */
    --listv;

    /* Function Body */
    i__1 = *ms - 1;
    for (i = 1; i <= i__1; ++i) {
	ii = *ms - i;
	if (listv[ii] > *iap) {
	    listv[ii + 1] = listv[ii];
	} else {
	    goto L2;
	}
/* L1: */
    }
    ii = 0;
L2:
    listv[ii + 1] = *iap;
    return 0;
} /* insert_ */

/* Subroutine */ int wrds_(len, lhead, lmin, lvail, lnew, lnext, ls, lv, 
	kount, klen, khead, kfwd, kbk)
integer *len, *lhead, *lmin, *lvail, *lnew, *lnext, *ls, *lv, *kount, *klen, *
	khead, *kfwd, *kbk;
{
/* 	WRITE(*,*) ' LEN = ',LEN */
/* 	WRITE(*,*)' LHEAD = ',LHEAD */
/* 	WRITE(*,*)' LMIN = ',LMIN */
/* 	WRITE(*,*)' LVAIL = ',LVAIL */
/* 	WRITE(*,*)' LNEW = ',LNEW */
/* 	WRITE(*,*)' KLEN = ', KLEN */
/* 	WRITE(*,*)' KHEAD = ',KHEAD */
/* 	WRITE(*,*)' LIST = ' */
/* 	WRITE(*,*)(ICK,'. ',LNEXT(ICK),LS(ICK),LV(ICK),ICK=1,LNEW) */
/* 	WRITE(*,*)' KOUNT = ' */
/* 	WRITE(*,*)(ICK,'. ',KOUNT(ICK),ICK=1,KLEN) */
/* 	WRITE(*,*)' KAND = ' */
/* 	WRITE(*,*)(ICK,'. ',KFWD(ICK),KBK(ICK),ICK=1,KLEN) */
    /* Parameter adjustments */
    --lv;
    --ls;
    --lnext;
    --kbk;
    --kfwd;
    --kount;

    /* Function Body */
    return 0;
} /* wrds_ */

/* Subroutine */ int fudge_(iap, ms, listv, len, lnext, ls, lv, lhead)
integer *iap, *ms, *listv, *len, *lnext, *ls, *lv, *lhead;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i;

    /* Parameter adjustments */
    --listv;
    --lv;
    --ls;
    --lnext;

    /* Function Body */
    i__1 = *ms;
    for (i = 1; i <= i__1; ++i) {
	if (listv[i] == *iap) {
	    goto L2;
	}
/* L1: */
    }
L2:
    lv[*lhead] = i;
    return 0;
} /* fudge_ */

