/* -*- mode: c; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- 
 *
 * $Id$
 *
 * Copyright (c) Erik Lindahl, David van der Spoel 2003,2004.
 * Copyright (c) Manuel Melo <manuel.nuno.melo@gmail.com> 2013,2014.
 * Coordinate compression (c) by Frans van Hoesel. 
 * XTC/TRR seeking and indexing (c) Manuel Melo.
 *
 *    This file is part of libxdrfile2.
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "xdrfile.h"
#include "xdrfile_trr.h"

#define BUFSIZE		128
#define GROMACS_MAGIC   1993

typedef struct		/* This struct describes the order and the	*/
/* sizes of the structs in a trjfile, sizes are given in bytes.	*/
{
    mybool  bDouble;        /* Double precision?                            */
    int	ir_size;	/* Backward compatibility		        */
    int	e_size;		/* Backward compatibility		        */
    int	box_size;	/* Non zero if a box is present			*/
    int   vir_size;       /* Backward compatibility		        */
    int   pres_size;      /* Backward compatibility		        */
    int	top_size;	/* Backward compatibility		        */
    int	sym_size;	/* Backward compatibility		        */
    int	x_size;		/* Non zero if coordinates are present		*/
    int	v_size;		/* Non zero if velocities are present		*/
    int	f_size;		/* Non zero if forces are present		*/

    int	natoms;		/* The total number of atoms			*/
    int	step;		/* Current step number				*/
    int	nre;		/* Backward compatibility		        */
    float	tf;		/* Current time					*/
    float	lambdaf;		/* Current value of lambda			*/
    double	td;		/* Current time					*/
    double	lambdad;		/* Current value of lambda			*/
} t_trnheader;

static int nFloatSize(t_trnheader *sh,int *nflsz)
{
    int nflsize=0;
  
    if (sh->box_size)
        nflsize = sh->box_size/(DIM*DIM);
    else if (sh->x_size)
        nflsize = sh->x_size/(sh->natoms*DIM);
    else if (sh->v_size)
        nflsize = sh->v_size/(sh->natoms*DIM);
    else if (sh->f_size)
        nflsize = sh->f_size/(sh->natoms*DIM);
    else 
        return exdrHEADER;
  
    if (((nflsize != sizeof(float)) && (nflsize != sizeof(double))))
        return exdrHEADER;
      
    *nflsz = nflsize;
  
    return exdrOK;
}

static int do_trnheader(XDRFILE *xd,mybool bRead,t_trnheader *sh)
{
	int magic=GROMACS_MAGIC;
	int nflsz,slen,result;
	char *version = "GMX_trn_file";
	char buf[BUFSIZE];
  
	if (xdrfile_read_int(&magic,1,xd) != 1)
		return exdrINT;
  
	if (bRead) 
    {
        if (xdrfile_read_int(&slen,1,xd) != 1)
            return exdrINT;
        if (slen != strlen(version)+1)
            return exdrSTRING;
        if (xdrfile_read_string(buf,BUFSIZE,xd) <= 0)
            return exdrSTRING;
    }
	else 
    {
        slen = strlen(version)+1;
        if (xdrfile_read_int(&slen,1,xd) != 1)
            return exdrINT;
        if (xdrfile_write_string(version,xd) != (strlen(version)+1) )
            return exdrSTRING;
    }
	if (xdrfile_read_int(&sh->ir_size,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->e_size,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->box_size,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->vir_size,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->pres_size,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->top_size,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->sym_size,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->x_size,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->v_size,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->f_size,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->natoms,1,xd) != 1)
		return exdrINT;
	
	if ((result = nFloatSize(sh,&nflsz)) != exdrOK)
		return result;
	sh->bDouble = (nflsz == sizeof(double));
	
	if (xdrfile_read_int(&sh->step,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->nre,1,xd) != 1)
		return exdrINT;
	if (sh->bDouble) 
    {
        if (xdrfile_read_double(&sh->td,1,xd) != 1)
            return exdrDOUBLE;
        sh->tf = sh->td;
        if (xdrfile_read_double(&sh->lambdad,1,xd) != 1)
            return exdrDOUBLE;
        sh->lambdaf = sh->lambdad;
    }
	else 
    {
        if (xdrfile_read_float(&sh->tf,1,xd) != 1)
            return exdrFLOAT;
        sh->td = sh->tf;
        if (xdrfile_read_float(&sh->lambdaf,1,xd) != 1)
            return exdrFLOAT;
        sh->lambdad = sh->lambdaf;
    }
  
    return exdrOK;
}

static int do_htrn(XDRFILE *xd,mybool bRead,t_trnheader *sh,
				   matrix box,rvec *x,rvec *v,rvec *f)
{
	double pvd[DIM*DIM];
	double *dx=NULL;
	float  pvf[DIM*DIM];
	float  *fx=NULL;
	int    i,j;
	
	if (sh->bDouble) 
	{
		if (sh->box_size != 0) 
        {
            if (!bRead) 
            {
                for(i=0; (i<DIM); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != box)
                        {
                            pvd[i*DIM+j] = box[i][j];
                        }
            }
            if (xdrfile_read_double(pvd,DIM*DIM,xd) == DIM*DIM) 
            {
                for(i=0; (i<DIM); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != box)
                        {
                            box[i][j] = pvd[i*DIM+j];
                        }
            }
            else
                return exdrDOUBLE;
        }
			
		if (sh->vir_size != 0) 
        {
            if (xdrfile_read_double(pvd,DIM*DIM,xd) != DIM*DIM) 
                return exdrDOUBLE;
        }
		
		if (sh->pres_size!= 0) 
        {
            if (xdrfile_read_double(pvd,DIM*DIM,xd) != DIM*DIM) 
                return exdrDOUBLE;
        }
		
		if ((sh->x_size != 0) || (sh->v_size != 0) || (sh->f_size != 0)) {
			dx = (double *)calloc(sh->natoms*DIM,sizeof(dx[0]));
			if (NULL == dx)
				return exdrNOMEM;
		}
		if (sh->x_size   != 0) 
        {
            if (!bRead) 
            {
                for(i=0; (i<sh->natoms); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != x)
                        {
                            dx[i*DIM+j] = x[i][j];
                        }
            }
            if (xdrfile_read_double(dx,sh->natoms*DIM,xd) == sh->natoms*DIM)
            {
                if (bRead) 
                {
                    for(i=0; (i<sh->natoms); i++)
                        for(j=0; (j<DIM); j++)
                            if (NULL != x)
                            {
                                x[i][j] = dx[i*DIM+j];
                            }
                }
            }
            else
                return exdrDOUBLE;
        }
		if (sh->v_size   != 0) 
        {
            if (!bRead) 
            {
                for(i=0; (i<sh->natoms); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != x)
                        {
                            dx[i*DIM+j] = v[i][j];
                        }
            }
            if (xdrfile_read_double(dx,sh->natoms*DIM,xd) == sh->natoms*DIM)
            {
                for(i=0; (i<sh->natoms); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != v)
                        {
                            v[i][j] = dx[i*DIM+j];
                        }
            }
            else
                return exdrDOUBLE;
        }
		if (sh->f_size   != 0) 
        {
            if (!bRead) 
            {
                for(i=0; (i<sh->natoms); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != x)
                        {
                            dx[i*DIM+j] = f[i][j];
                        }
            }
            if (xdrfile_read_double(dx,sh->natoms*DIM,xd) == sh->natoms*DIM)
            {
                for(i=0; (i<sh->natoms); i++)
                {
                    for(j=0; (j<DIM); j++)
                    {
                        if (NULL != f)
                        {
                            f[i][j] = dx[i*DIM+j];
                        }
                    }
                }
            }
            else
                return exdrDOUBLE;
        }
		if ((sh->x_size != 0) || (sh->v_size != 0) || (sh->f_size != 0)) {
			free(dx);
		}
	}
	else
		/* Float */
	{
		if (sh->box_size != 0) 
        {
            if (!bRead) 
            {
                for(i=0; (i<DIM); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != box)
                        {
                            pvf[i*DIM+j] = box[i][j];
                        }
            }
            if (xdrfile_read_float(pvf,DIM*DIM,xd) == DIM*DIM) 
            {
                for(i=0; (i<DIM); i++)
                {
                    for(j=0; (j<DIM); j++)
                    {
                        if (NULL != box)
                        {
                            box[i][j] = pvf[i*DIM+j];
                        }
                    }
                }
            }
            else
                return exdrFLOAT;
        }
			
		if (sh->vir_size != 0) 
        {
            if (xdrfile_read_float(pvf,DIM*DIM,xd) != DIM*DIM) 
                return exdrFLOAT;
        }
		
		if (sh->pres_size!= 0) 
        {
            if (xdrfile_read_float(pvf,DIM*DIM,xd) != DIM*DIM) 
                return exdrFLOAT;
        }
		
		if ((sh->x_size != 0) || (sh->v_size != 0) || (sh->f_size != 0)) {
			fx = (float *)calloc(sh->natoms*DIM,sizeof(fx[0]));
			if (NULL == fx)
				return exdrNOMEM;
		}
		if (sh->x_size   != 0) 
        {
            if (!bRead) 
            {
                for(i=0; (i<sh->natoms); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != x)
                        {
                            fx[i*DIM+j] = x[i][j];
                        }
            }
            if (xdrfile_read_float(fx,sh->natoms*DIM,xd) == sh->natoms*DIM)
            {
                if (bRead) 
                {
                    for(i=0; (i<sh->natoms); i++)
                        for(j=0; (j<DIM); j++)
                            if (NULL != x)
                                x[i][j] = fx[i*DIM+j];
                }
            }
            else
                return exdrFLOAT;
        }
		if (sh->v_size   != 0) 
        {
            if (!bRead) 
            {
                for(i=0; (i<sh->natoms); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != x)
                        {
                            fx[i*DIM+j] = v[i][j];
                        }
            }
            if (xdrfile_read_float(fx,sh->natoms*DIM,xd) == sh->natoms*DIM)
            {
                for(i=0; (i<sh->natoms); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != v)
                            v[i][j] = fx[i*DIM+j];
            }
            else
                return exdrFLOAT;
        }
		if (sh->f_size   != 0) 
        {
           if (!bRead) 
            {
                for(i=0; (i<sh->natoms); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != x)
                        {
                            fx[i*DIM+j] = f[i][j];
                        }
            }
             if (xdrfile_read_float(fx,sh->natoms*DIM,xd) == sh->natoms*DIM)
            {
                for(i=0; (i<sh->natoms); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != f)
                            f[i][j] = fx[i*DIM+j];
            }
            else
                return exdrFLOAT;
        }
		if ((sh->x_size != 0) || (sh->v_size != 0) || (sh->f_size != 0)) {
			free(fx);
		}
	}
	return exdrOK;
}

static int do_trn(XDRFILE *xd,mybool bRead,int *step,float *t,float *lambda,
				  matrix box,int *natoms,rvec *x,rvec *v,rvec *f, int *has_prop)
{
    t_trnheader *sh;
    int result;
  
    sh = (t_trnheader *)calloc(1,sizeof(*sh));
  
    if (!bRead) {
        sh->box_size = (NULL != box) ? sizeof(matrix):0;
        sh->x_size   = ((NULL != x) ? (*natoms*sizeof(x[0])):0);
        sh->v_size   = ((NULL != v) ? (*natoms*sizeof(v[0])):0);
        sh->f_size   = ((NULL != f) ? (*natoms*sizeof(f[0])):0);
        sh->natoms = *natoms;
        sh->step   = *step;
        sh->nre    = 0;
        sh->td      = *t;
        sh->lambdad = *lambda;
        sh->tf      = *t;
        sh->lambdaf = *lambda;
    }
    if ((result = do_trnheader(xd,bRead,sh)) != exdrOK)
        return result;
    if (bRead) {
        *natoms = sh->natoms;
        *step   = sh->step;
        *t      = sh->td;
        *lambda = sh->lambdad;
        /* Flag what we read */
        if (sh->x_size)
            *has_prop |= HASX;
        if (sh->v_size)
            *has_prop |= HASV;
        if (sh->f_size)
            *has_prop |= HASF;
    }
    if ((result = do_htrn(xd,bRead,sh,box,x,v,f)) != exdrOK)
        return result;

    free(sh);
  
    return exdrOK;
}

/************************************************************
 *
 *  The following routines are the exported ones
 *
 ************************************************************/
 
int read_trr_natoms(char *fn,int *natoms)
{
	XDRFILE *xd;
	t_trnheader sh;
	int  result;
	
	xd = xdrfile_open(fn,"r");
	if (NULL == xd)
		return exdrFILENOTFOUND;
	if ((result = do_trnheader(xd,1,&sh)) != exdrOK)
		return result;
	xdrfile_close(xd);
	*natoms = sh.natoms;
	
	return exdrOK;
}

int read_trr_n_frames(char *fn, int *n_frames, int64_t **offsets)
{
	XDRFILE *xd;
	t_trnheader sh;
	float time, lambda;
	int result, framebytes, est_nframes, totalframebytes;
    int64_t filesize, frame_offset;
	
	if ((xd = xdrfile_open(fn,"r"))==NULL)
		return exdrFILENOTFOUND;
    if (xdr_seek(xd, 0L, SEEK_END) != exdrOK)
    {
        xdrfile_close(xd);
        return exdrNR;
    }
    filesize = xdr_tell(xd);
    if (xdr_seek(xd, 0L, SEEK_SET) != exdrOK)
    {
        xdrfile_close(xd);
        return exdrNR;
    }

	if ((result = do_trnheader(xd,1,&sh)) != exdrOK)
    {
        xdrfile_close(xd);
		return result;
    }

    framebytes = sh.ir_size + sh.e_size + sh.box_size +
                 sh.vir_size + sh.pres_size + sh.top_size +
                 sh.sym_size + sh.x_size + sh.v_size + sh.f_size;

    est_nframes = (int) (filesize/((int64_t) (framebytes + TRR_MIN_HEADER_SIZE)) + 1); // add one because it'd be easy to underestimate low frame numbers. 
    est_nframes += est_nframes/5;

    /* Allocate memory for the frame index array */
    if ((*offsets=(int64_t *)malloc(sizeof(int64_t)*est_nframes))==NULL)
    {
        xdrfile_close(xd);
        return exdrNOMEM;
    }

    (*offsets)[0] = 0L;
    *n_frames = 1;
    while (1)
    {
        if (xdr_seek(xd, (int64_t) (framebytes), SEEK_CUR) != exdrOK) {
            free(*offsets);
            xdrfile_close(xd);
            return exdrNR;
        }
        frame_offset = xdr_tell(xd); /* Store it now, before we read the header */
        if ((result = do_trnheader(xd,1,&sh)) != exdrOK) /* Interpreting as EOF */
            break;
        /* Read was successful; this is another frame */
        /* Check if we need to enlarge array */
        if (*n_frames == est_nframes){
            est_nframes += est_nframes/5 + 1; // Increase in 20% stretches
            if ((*offsets = realloc(*offsets, sizeof(int64_t)*est_nframes))==NULL)
            {
                xdrfile_close(xd);
                return exdrNOMEM;
            }
        }
        (*offsets)[*n_frames] = frame_offset;
        (*n_frames)++;
        /* Calculate how much to skip this time */
        framebytes = sh.ir_size + sh.e_size + sh.box_size +
                     sh.vir_size + sh.pres_size + sh.top_size +
                     sh.sym_size + sh.x_size + sh.v_size + sh.f_size;
    }
	xdrfile_close(xd);
	return exdrOK;
}


int write_trr(XDRFILE *xd,int natoms,int step,float t,float lambda,
			  matrix box,rvec *x,rvec *v,rvec *f)
{
    int *plcholder;
	return do_trn(xd,0,&step,&t,&lambda,box,&natoms,x,v,f, plcholder);
}

int read_trr(XDRFILE *xd,int natoms,int *step,float *t,float *lambda,
			 matrix box,rvec *x,rvec *v,rvec *f, int *has_prop)
{
	return do_trn(xd,1,step,t,lambda,box,&natoms,x,v,f,has_prop);
}

