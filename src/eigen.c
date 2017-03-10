/*
IgPhyML: a program that computes maximum likelihood phylogenies under
non-reversible codon models designed for antibody lineages.

Copyright (C) Kenneth B Hoehn. Sept 2016 onward.

built upon

codonPHYML: a program that  computes maximum likelihood phylogenies from
CODON homologous sequences.

Copyright (C) Marcelo Serrano Zanetti. Oct 2010 onward.

built upon

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

/***********************************************************
*  This eigen() routine works for eigenvalue/vector analysis
*         for real general square matrix A
*         A will be destroyed
*         rr,ri are vectors containing eigenvalues
*         vr,vi are matrices containing (right) eigenvectors
*
*              A*[vr+vi*i] = [vr+vi*i] * diag{rr+ri*i}
*
*  Algorithm: Handbook for Automatic Computation, vol 2
*             by Wilkinson and Reinsch, 1971
*             most of source codes were taken from a public domain
*             solftware called MATCALC.
*  Credits:   to the authors of MATCALC
*
*  return     -1 not converged
*              0 no phycomplex eigenvalues/vectors
*              1 phycomplex eigenvalues/vectors
*  Tianlin Wang at University of Illinois
*  Thu May  6 15:22:31 CDT 1993
***************************************************************/

#include "eigen.h"


#define BASE        2    /* base of floating point arithmetic */

/* no. of digits to the base BASE in the fraction */
#define DIGITS 40
/* 
#define DIGITS 53 
*/

#define MAXITER    30    /* max2. no. of iterations to converge */

#define pos(i,j,n)      ((i)*(n)+(j))


/* rr/vr : real parts of eigen values/vectors */
/* ri/vi : imaginary part s of eigen values/vectors */

int Eigen(int job, phydbl *A, int n, phydbl *rr, phydbl *ri, 
          phydbl *vr, phydbl *vi, phydbl *work)
{    
  /* job=0: eigen values only
  1: both eigen values and eigen vectors 
  phydbl w[n*2]: work space
  */
  int low,hi,i,j,k, it, istate=0;
  phydbl tiny, t; 
  
  /*     tiny=SQRT(POW((phydbl)BASE,(phydbl)(1-(int)DIGITS))); */
  tiny=FLT_MIN;
  
  balance(A,n,&low,&hi,work);
  elemhess(job,A,n,low,hi,vr,vi, (int*)(work+n));
  if (-1 == realeig(job,A,n,low,hi,rr,ri,vr,vi)) return (-1);
  if (job) unbalance(n,vr,vi,low,hi,work);
  
  /* sort, added by Z. Yang */
  for (i=0; i<n; i++) {
    for (j=i+1,it=i,t=rr[i]; j<n; j++)
      if (t<rr[j]) { t=rr[j]; it=j; }
      rr[it]=rr[i];   rr[i]=t;
    t=ri[it];       ri[it]=ri[i];  ri[i]=t;
    for (k=0; k<n; k++) {
      t=vr[k*n+it];  vr[k*n+it]=vr[k*n+i];  vr[k*n+i]=t;
      t=vi[k*n+it];  vi[k*n+it]=vi[k*n+i];  vi[k*n+i]=t;
    }
    if (FABS(ri[i])>tiny) istate=1;
  }
  
  return (istate);
}

/* phycomplex funcctions
*/

phycomplex compl (phydbl re,phydbl im)
{
    phycomplex r;

    r.re = re;
    r.im = im;
    return(r);
}

phycomplex _conj (phycomplex a)
{
    a.im = -a.im;
    return(a);
}


phycomplex cplus (phycomplex a, phycomplex b)
{
   phycomplex c;
   c.re = a.re+b.re;  
   c.im = a.im+b.im;   
   return (c);
}

phycomplex cminus (phycomplex a, phycomplex b)
{
   phycomplex c;
   c.re = a.re-b.re;  
   c.im = a.im-b.im;   
   return (c);
}

phycomplex cby (phycomplex a, phycomplex b)
{
   phycomplex c;
   c.re = a.re*b.re-a.im*b.im ;
   c.im = a.re*b.im+a.im*b.re ;
   return (c);
}

phycomplex cdiv (phycomplex a,phycomplex b)
{
    phydbl ratio, den;
    phycomplex c;

    if (FABS(b.re) <= FABS(b.im)) {
        ratio = b.re / b.im;
        den = b.im * (1 + ratio * ratio);
        c.re = (a.re * ratio + a.im) / den;
        c.im = (a.im * ratio - a.re) / den;
    }
    else {
        ratio = b.im / b.re;
        den = b.re * (1 + ratio * ratio);
        c.re = (a.re + a.im * ratio) / den;
        c.im = (a.im - a.re * ratio) / den;
    }
    return(c);
}

/* phycomplex local_cexp (phycomplex a) */
/* { */
/*    phycomplex c; */
/*    c.re = EXP(a.re); */
/*    if (FABS(a.im)==0) c.im = 0;  */
/*    else  { c.im = c.re*sin(a.im); c.re*=cos(a.im); } */
/*    return (c); */
/* } */

phycomplex cfactor (phycomplex x, phydbl a)
{
   phycomplex c;
   c.re = a*x.re; 
   c.im = a*x.im;
   return (c);
}

int cxtoy (phycomplex *x, phycomplex *y, int n)
{
   int i;
   For (i,n) y[i]=x[i];
   return (0);
}

int cmatby (phycomplex *a, phycomplex *b, phycomplex *c, int n,int m,int k)
/* a[n*m], b[m*k], c[n*k]  ......  c = a*b
*/
{
   int i,j,i1;
   phycomplex t;

   For (i,n)  For(j,k) {
       for (i1=0,t=compl(0,0); i1<m; i1++)  
           t = cplus (t, cby(a[i*m+i1],b[i1*k+j]));
       c[i*k+j] = t;
   }
   return (0);
}

int cmatinv( phycomplex *x, int n, int m, phydbl *space)
{
/* x[n*m]  ... m>=n
*/
   int i,j,k, *irow=(int*) space;
   phydbl xmaxsize, ee=1e-20;
   phycomplex xmax, t,t1;

   For(i,n)  {
       xmaxsize = 0.;
       for (j=i; j<n; j++) {
          if ( xmaxsize < csize (x[j*m+i]))  {
               xmaxsize = csize (x[j*m+i]);
               xmax = x[j*m+i];
               irow[i] = j;
          }
       }
       if (xmaxsize < ee)   {
           PhyML_Printf("\nDet goes to zero at %8d!\t\n", i+1);
           return(-1);
       }
       if (irow[i] != i) {
           For(j,m) {
                t = x[i*m+j];
                x[i*m+j] = x[irow[i]*m+j];
                x[ irow[i]*m+j] = t;
           }
       }
       t = cdiv (compl(1,0), x[i*m+i]);
       For(j,n) {
           if (j == i) continue;
           t1 = cby (t,x[j*m+i]);
           For(k,m)  x[j*m+k] = cminus (x[j*m+k], cby(t1,x[i*m+k]));
           x[j*m+i] = cfactor (t1, -1);
       }
       For(j,m)   x[i*m+j] = cby (x[i*m+j], t);
       x[i*m+i] = t;
   }                         
   for (i=n-1; i>=0; i--) {
        if (irow[i] == i) continue;
        For(j,n)  {
           t = x[j*m+i];
           x[j*m+i] = x[j*m+irow[i]];
           x[ j*m+irow[i]] = t;
        }
   }
   return (0);
}


void balance(phydbl *mat, int n,int *low, int *hi, phydbl *scale)
{
/* Balance a matrix for calculation of eigenvalues and eigenvectors
*/
    phydbl c,f,g,r,s;
    int i,j,k,l,done;
        /* search for rows isolating an eigenvalue and push them down */
    for (k = n - 1; k >= 0; k--) 
      {
        for (j = k; j >= 0; j--) 
	  {
	    for (i = 0; i <= k; i++) 
	      {
		if (i != j && FABS(mat[pos(j,i,n)]) > SMALL) break;
	      }

            if (i > k) {
                scale[k] = j;

                if (j != k) {
                    for (i = 0; i <= k; i++) {
                       c = mat[pos(i,j,n)];
                       mat[pos(i,j,n)] = mat[pos(i,k,n)];
                       mat[pos(i,k,n)] = c;
                    }

                    for (i = 0; i < n; i++) {
                       c = mat[pos(j,i,n)];
                       mat[pos(j,i,n)] = mat[pos(k,i,n)];
                       mat[pos(k,i,n)] = c;
                    }
                }
                break;
            }
        }
        if (j < 0) break;
    }

    /* search for columns isolating an eigenvalue and push them left */

    for (l = 0; l <= k; l++) {
        for (j = l; j <= k; j++) {
            for (i = l; i <= k; i++) {
	      if (i != j && FABS(mat[pos(i,j,n)]) > SMALL) break;
            }
            if (i > k) {
                scale[l] = j;
                if (j != l) {
                    for (i = 0; i <= k; i++) {
                       c = mat[pos(i,j,n)];
                       mat[pos(i,j,n)] = mat[pos(i,l,n)];
                       mat[pos(i,l,n)] = c;
                    }

                    for (i = l; i < n; i++) {
                       c = mat[pos(j,i,n)];
                       mat[pos(j,i,n)] = mat[pos(l,i,n)];
                       mat[pos(l,i,n)] = c;
                    }
                }

                break;
            }
        }

        if (j > k) break;
    }

    *hi = k;
    *low = l;

    /* balance the submatrix in rows l through k */

    for (i = l; i <= k; i++) {
        scale[i] = 1;
    }

    do {
        for (done = 1,i = l; i <= k; i++) {
            for (c = 0,r = 0,j = l; j <= k; j++) {
                if (j != i) {
                    c += FABS(mat[pos(j,i,n)]);
                    r += FABS(mat[pos(i,j,n)]);
                }
            }

/*             if (c != 0 && r != 0) {  */
            if (FABS(c) > SMALL && FABS(r) > SMALL) {
               g = r / BASE;
                f = 1;
                s = c + r;

                while (c < g) {
                    f *= BASE;
                    c *= BASE * BASE;
                }

                g = r * BASE;

                while (c >= g) {
                    f /= BASE;
                    c /= BASE * BASE;
                }

                if ((c + r) / f < 0.95 * s) {
                    done = 0;
                    g = 1 / f;
                    scale[i] *= f;

                    for (j = l; j < n; j++) {
                        mat[pos(i,j,n)] *= g;
                    }

                    for (j = 0; j <= k; j++) {
                        mat[pos(j,i,n)] *= f;
                    }
                }
            }
        }
    } while (!done);
}


/*
 * Transform back eigenvectors of a balanced matrix
 * into the eigenvectors of the original matrix
 */
void unbalance(int n,phydbl *vr,phydbl *vi, int low, int hi, phydbl *scale)
{
    int i,j,k;
    phydbl tmp;

    for (i = low; i <= hi; i++) {
        for (j = 0; j < n; j++) {
            vr[pos(i,j,n)] *= scale[i];
            vi[pos(i,j,n)] *= scale[i];
        }
    }

    for (i = low - 1; i >= 0; i--) {
        if ((k = (int)scale[i]) != i) {
            for (j = 0; j < n; j++) {
                tmp = vr[pos(i,j,n)];
                vr[pos(i,j,n)] = vr[pos(k,j,n)];
                vr[pos(k,j,n)] = tmp;

                tmp = vi[pos(i,j,n)];
                vi[pos(i,j,n)] = vi[pos(k,j,n)];
                vi[pos(k,j,n)] = tmp;        
            }
        }
    }

    for (i = hi + 1; i < n; i++) {
        if ((k = (int)scale[i]) != i) {
            for (j = 0; j < n; j++) {
                tmp = vr[pos(i,j,n)];
                vr[pos(i,j,n)] = vr[pos(k,j,n)];
                vr[pos(k,j,n)] = tmp;

                tmp = vi[pos(i,j,n)];
                vi[pos(i,j,n)] = vi[pos(k,j,n)];
                vi[pos(k,j,n)] = tmp;        
            }
        }
    }
}

/*
 * Reduce the submatrix in rows and columns low through hi of real matrix mat to
 * Hessenberg form by elementary similarity transformations
 */
void elemhess(int job,phydbl *mat,int n,int low,int hi, phydbl *vr,
              phydbl *vi, int *work)
{
/* work[n] */
    int i,j,m;
    phydbl x,y;

    for (m = low + 1; m < hi; m++) {
        for (x = 0,i = m,j = m; j <= hi; j++) {
            if (FABS(mat[pos(j,m-1,n)]) > FABS(x)) {
                x = mat[pos(j,m-1,n)];
                i = j;
            }
        }

        if ((work[m] = i) != m) {
            for (j = m - 1; j < n; j++) {
               y = mat[pos(i,j,n)];
               mat[pos(i,j,n)] = mat[pos(m,j,n)];
               mat[pos(m,j,n)] = y;
            }

            for (j = 0; j <= hi; j++) {
               y = mat[pos(j,i,n)];
               mat[pos(j,i,n)] = mat[pos(j,m,n)];
               mat[pos(j,m,n)] = y;
            }
        }

        if (FABS(x) > SMALL) {
            for (i = m + 1; i <= hi; i++) {
	      if (FABS(y = mat[pos(i,m-1,n)]) > SMALL) {
                    y = mat[pos(i,m-1,n)] = y / x;

                    for (j = m; j < n; j++) {
                        mat[pos(i,j,n)] -= y * mat[pos(m,j,n)];
                    }

                    for (j = 0; j <= hi; j++) {
                        mat[pos(j,m,n)] += y * mat[pos(j,i,n)];
                    }
                }
            }
        }
    }
    if (job) {
       for (i=0; i<n; i++) {
          for (j=0; j<n; j++) {
             vr[pos(i,j,n)] = 0.0; vi[pos(i,j,n)] = 0.0;
          }
          vr[pos(i,i,n)] = 1.0;
       }

       for (m = hi - 1; m > low; m--) {
          for (i = m + 1; i <= hi; i++) {
             vr[pos(i,m,n)] = mat[pos(i,m-1,n)];
          }

         if ((i = work[m]) != m) {
            for (j = m; j <= hi; j++) {
               vr[pos(m,j,n)] = vr[pos(i,j,n)];
               vr[pos(i,j,n)] = 0.0;
            }
            vr[pos(i,m,n)] = 1.0;
         }
      }
   }
}

/*
 * Calculate eigenvalues and eigenvectors of a real upper Hessenberg matrix
 * Return 1 if converges successfully and 0 otherwise
 */
 
int realeig(int job,phydbl *mat,int n,int low, int hi, phydbl *valr,
      phydbl *vali, phydbl *vr,phydbl *vi)
{
   phycomplex v;
   phydbl p=.0,q=.0,r=.0,s=.0,t,w,x,y,z=0,ra,sa,norm,eps;
   int niter,en,i,j,k,l,m;
   phydbl precision  = POW((phydbl)BASE,(phydbl)(1-(int)DIGITS));

   eps = precision;
   for (i=0; i<n; i++) {
      valr[i]=0.0;
      vali[i]=0.0;
   }
      /* store isolated roots and calculate norm */
   for (norm = 0,i = 0; i < n; i++) {
      for (j = MAX(0,i-1); j < n; j++) {
         norm += FABS(mat[pos(i,j,n)]);
      }
      if (i < low || i > hi) valr[i] = mat[pos(i,i,n)];
   }
   t = 0;
   en = hi;

   while (en >= low) {
      niter = 0;
      for (;;) {

       /* look for single small subdiagonal element */

         for (l = en; l > low; l--) {
            s = FABS(mat[pos(l-1,l-1,n)]) + FABS(mat[pos(l,l,n)]);
            if (FABS(s) < SMALL) s = norm;
            if (FABS(mat[pos(l,l-1,n)]) <= eps * s) break;
         }

         /* form shift */

         x = mat[pos(en,en,n)];

         if (l == en) {             /* one root found */
            valr[en] = x + t;
            if (job) mat[pos(en,en,n)] = x + t;
            en--;
            break;
         }

         y = mat[pos(en-1,en-1,n)];
         w = mat[pos(en,en-1,n)] * mat[pos(en-1,en,n)];

         if (l == en - 1) {                /* two roots found */
            p = (y - x) / 2;
            q = p * p + w;
            z = SQRT(FABS(q));
            x += t;
            if (job) {
               mat[pos(en,en,n)] = x;
               mat[pos(en-1,en-1,n)] = y + t;
            }
            if (q < 0) {                /* phycomplex pair */
               valr[en-1] = x+p;
               vali[en-1] = z;
               valr[en] = x+p;
               vali[en] = -z;
            }
            else {                      /* real pair */
               z = (p < 0) ? p - z : p + z;
               valr[en-1] = x + z;
               valr[en] = (FABS(z) < SMALL) ? x + z : x - w / z;
               if (job) {
                  x = mat[pos(en,en-1,n)];
                  s = FABS(x) + FABS(z);
                  p = x / s;
                  q = z / s;
                  r = SQRT(p*p+q*q);
                  p /= r;
                  q /= r;
                  for (j = en - 1; j < n; j++) {
                     z = mat[pos(en-1,j,n)];
                     mat[pos(en-1,j,n)] = q * z + p *
                     mat[pos(en,j,n)];
                     mat[pos(en,j,n)] = q * mat[pos(en,j,n)] - p*z;
                  }
                  for (i = 0; i <= en; i++) {
                     z = mat[pos(i,en-1,n)];
                     mat[pos(i,en-1,n)] = q * z + p * mat[pos(i,en,n)];
                     mat[pos(i,en,n)] = q * mat[pos(i,en,n)] - p*z;
                  }
                  for (i = low; i <= hi; i++) {
                     z = vr[pos(i,en-1,n)];
                     vr[pos(i,en-1,n)] = q*z + p*vr[pos(i,en,n)];
                     vr[pos(i,en,n)] = q*vr[pos(i,en,n)] - p*z;
                  }
               }
            }
            en -= 2;
            break;
         }
         if (niter == MAXITER) return(-1);
         if (niter != 0 && niter % 10 == 0) {
            t += x;
            for (i = low; i <= en; i++) mat[pos(i,i,n)] -= x;
            s = FABS(mat[pos(en,en-1,n)]) + FABS(mat[pos(en-1,en-2,n)]);
            x = y = 0.75 * s;
            w = -0.4375 * s * s;
         }
         niter++;
           /* look for two consecutive small subdiagonal elements */
         for (m = en - 2; m >= l; m--) {
            z = mat[pos(m,m,n)];
            r = x - z;
            s = y - z;
            p = (r * s - w) / mat[pos(m+1,m,n)] + mat[pos(m,m+1,n)];
            q = mat[pos(m+1,m+1,n)] - z - r - s;
            r = mat[pos(m+2,m+1,n)];
            s = FABS(p) + FABS(q) + FABS(r);
            p /= s;
            q /= s;
            r /= s;
            if (m == l || FABS(mat[pos(m,m-1,n)]) * (FABS(q)+FABS(r)) <=
                eps * (FABS(mat[pos(m-1,m-1,n)]) + FABS(z) +
                FABS(mat[pos(m+1,m+1,n)])) * FABS(p)) break;
         }
         for (i = m + 2; i <= en; i++) mat[pos(i,i-2,n)] = 0;
         for (i = m + 3; i <= en; i++) mat[pos(i,i-3,n)] = 0;
             /* phydbl QR step involving rows l to en and columns m to en */
         for (k = m; k < en; k++) {
            if (k != m) {
               p = mat[pos(k,k-1,n)];
               q = mat[pos(k+1,k-1,n)];
               r = (k == en - 1) ? 0 : mat[pos(k+2,k-1,n)];
               if (FABS(x = FABS(p) + FABS(q) + FABS(r)) < SMALL) continue;
               p /= x;
               q /= x;
               r /= x;
            }
            s = SQRT(p*p+q*q+r*r);
            if (p < 0) s = -s;
            if (k != m) {
               mat[pos(k,k-1,n)] = -s * x;
            }
            else if (l != m) {
               mat[pos(k,k-1,n)] = -mat[pos(k,k-1,n)];
            }
            p += s;
            x = p / s;
            y = q / s;
            z = r / s;
            q /= p;
            r /= p;
                /* row modification */
            for (j = k; j <= (!job ? en : n-1); j++){
               p = mat[pos(k,j,n)] + q * mat[pos(k+1,j,n)];
               if (k != en - 1) {
                  p += r * mat[pos(k+2,j,n)];
                  mat[pos(k+2,j,n)] -= p * z;
               }
               mat[pos(k+1,j,n)] -= p * y;
               mat[pos(k,j,n)] -= p * x;
            }
            j = MIN(en,k+3);
              /* column modification */
            for (i = (!job ? l : 0); i <= j; i++) {
               p = x * mat[pos(i,k,n)] + y * mat[pos(i,k+1,n)];
               if (k != en - 1) {
                  p += z * mat[pos(i,k+2,n)];
                  mat[pos(i,k+2,n)] -= p*r;
               }
               mat[pos(i,k+1,n)] -= p*q;
               mat[pos(i,k,n)] -= p;
            }
            if (job) {             /* accumulate transformations */
               for (i = low; i <= hi; i++) {
                  p = x * vr[pos(i,k,n)] + y * vr[pos(i,k+1,n)];
                  if (k != en - 1) {
                     p += z * vr[pos(i,k+2,n)];
                     vr[pos(i,k+2,n)] -= p*r;
                  }
                  vr[pos(i,k+1,n)] -= p*q;
                  vr[pos(i,k,n)] -= p;
               }
            }
         }
      }
   }

   if (!job) return(0);
   if (FABS(norm) > SMALL) {
       /* back substitute to find vectors of upper triangular form */
      for (en = n-1; en >= 0; en--) {
         p = valr[en];
         if ((q = vali[en]) < 0) {            /* phycomplex vector */
            m = en - 1;
            if (FABS(mat[pos(en,en-1,n)]) > FABS(mat[pos(en-1,en,n)])) {
               mat[pos(en-1,en-1,n)] = q / mat[pos(en,en-1,n)];
               mat[pos(en-1,en,n)] = (p - mat[pos(en,en,n)]) /
                     mat[pos(en,en-1,n)];
            }
            else {
               v = cdiv(compl(0.0,-mat[pos(en-1,en,n)]),
                    compl(mat[pos(en-1,en-1,n)]-p,q));
               mat[pos(en-1,en-1,n)] = v.re;
               mat[pos(en-1,en,n)] = v.im;
            }
            mat[pos(en,en-1,n)] = 0;
            mat[pos(en,en,n)] = 1;
            for (i = en - 2; i >= 0; i--) {
               w = mat[pos(i,i,n)] - p;
               ra = 0;
               sa = mat[pos(i,en,n)];
               for (j = m; j < en; j++) {
                  ra += mat[pos(i,j,n)] * mat[pos(j,en-1,n)];
                  sa += mat[pos(i,j,n)] * mat[pos(j,en,n)];
               }
               if (vali[i] < 0) {
                  z = w;
                  r = ra;
                  s = sa;
               }
               else {
                  m = i;
                  if (FABS(vali[i]) < SMALL) {
                     v = cdiv(compl(-ra,-sa),compl(w,q));
                     mat[pos(i,en-1,n)] = v.re;
                     mat[pos(i,en,n)] = v.im;
                  }
                  else {                      /* solve phycomplex equations */
                     x = mat[pos(i,i+1,n)];
                     y = mat[pos(i+1,i,n)];
                     v.re = (valr[i]- p)*(valr[i]-p) + vali[i]*vali[i] - q*q;
                     v.im = (valr[i] - p)*2*q;
                     if (FABS(v.re) + FABS(v.im) < SMALL) {
                        v.re = eps * norm * (FABS(w) +
                                FABS(q) + FABS(x) + FABS(y) + FABS(z));
                     }
                     v = cdiv(compl(x*r-z*ra+q*sa,x*s-z*sa-q*ra),v);
                     mat[pos(i,en-1,n)] = v.re;
                     mat[pos(i,en,n)] = v.im;
                     if (FABS(x) > FABS(z) + FABS(q)) {
                        mat[pos(i+1,en-1,n)] = 
                             (-ra - w * mat[pos(i,en-1,n)] +
                             q * mat[pos(i,en,n)]) / x;
                        mat[pos(i+1,en,n)] = (-sa - w * mat[pos(i,en,n)] -
                             q * mat[pos(i,en-1,n)]) / x;
                     }
                     else {
                        v = cdiv(compl(-r-y*mat[pos(i,en-1,n)],
                             -s-y*mat[pos(i,en,n)]),compl(z,q));
                        mat[pos(i+1,en-1,n)] = v.re;
                        mat[pos(i+1,en,n)] = v.im;
                     }
                  }
               }
            }
         }
         else if (FABS(q) < SMALL) {                             /* real vector */
            m = en;
            mat[pos(en,en,n)] = 1;
            for (i = en - 1; i >= 0; i--) {
               w = mat[pos(i,i,n)] - p;
               r = mat[pos(i,en,n)];
               for (j = m; j < en; j++) {
                  r += mat[pos(i,j,n)] * mat[pos(j,en,n)];
               }
               if (vali[i] < 0) {
                  z = w;
                  s = r;
               }
               else {
                  m = i;
                  if (FABS(vali[i]) < SMALL) {
                     if (FABS(t = w) < SMALL) t = eps * norm;
                     mat[pos(i,en,n)] = -r / t;
                  }
                  else {            /* solve real equations */
                     x = mat[pos(i,i+1,n)];
                     y = mat[pos(i+1,i,n)];
                     q = (valr[i] - p) * (valr[i] - p) + vali[i]*vali[i];
                     t = (x * s - z * r) / q;
                     mat[pos(i,en,n)] = t;
                     if (FABS(x) <= FABS(z)) {
                        mat[pos(i+1,en,n)] = (-s - y * t) / z;
                     }
                     else {
                        mat[pos(i+1,en,n)] = (-r - w * t) / x;
                     }
                  }
               }
            }
         }
      }
             /* vectors of isolated roots */
      for (i = 0; i < n; i++) {
         if (i < low || i > hi) {
            for (j = i; j < n; j++) {
               vr[pos(i,j,n)] = mat[pos(i,j,n)];
            }
         }
      }
       /* multiply by transformation matrix */

      for (j = n-1; j >= low; j--) {
         m = MIN(j,hi);
         for (i = low; i <= hi; i++) {
            for (z = 0,k = low; k <= m; k++) {
               z += vr[pos(i,k,n)] * mat[pos(k,j,n)];
            }
            vr[pos(i,j,n)] = z;
         }
      }
   }
    /* rearrange phycomplex eigenvectors */
   for (j = 0; j < n; j++) {
     if (FABS(vali[j]) > SMALL) {
         for (i = 0; i < n; i++) {
            vi[pos(i,j,n)] = vr[pos(i,j+1,n)];
            vr[pos(i,j+1,n)] = vr[pos(i,j,n)];
            vi[pos(i,j+1,n)] = -vi[pos(i,j,n)];
         }
         j++;
      }
   }
   return(0);
}


#define LUDCMP_TINY 1.0e-20;

int ludcmp(phydbl **a, int n, phydbl *d)
{
   int i,imax,j,k;
   phydbl big,dum,sum,temp;
   phydbl *vv;
   
   imax = 0;
   vv = (phydbl *)mCalloc(n,sizeof(phydbl));

   *d=1.0;
   for (i=0;i<n;i++) 
     {
       big=0.0;
       for (j=0;j<n;j++)
         if ((temp=FABS(a[i][j])) > big) big=temp;
       if (FABS(big) < SMALL) Exit("\n. Singular matrix in routine LUDCMP");
       vv[i]=1.0/big;
     }
   for (j=0;j<n;j++) 
     {
       for (i=0;i<j;i++) 
	 {
	   sum=a[i][j];
	   for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
	   a[i][j]=sum;
	 }
      big=0.0;
      for (i=j;i<n;i++) {
	sum=a[i][j];
	for (k=0;k<j;k++)
	  sum -= a[i][k]*a[k][j];
	a[i][j]=sum;
	if ((dum=vv[i]*FABS(sum)) >= big) 
	  {
            big=dum;
            imax=i;
	  }
      }
      if (j != imax) 
	{
	  for (k=0;k<n;k++) 
	    {
	      dum=a[imax][k];
	      a[imax][k]=a[j][k];
	      a[j][k]=dum;
	    }
	  *d = -(*d);
	  vv[imax]=vv[j];
	}
      if (FABS(a[j][j]) < SMALL) a[j][j]=LUDCMP_TINY;
      if (j != n) {
	dum=1.0/(a[j][j]);
	for (i=j+1;i<n;i++) a[i][j] *= dum;
      }
     }
   free(vv);
   return(0);
}

void det(phydbl **a, int n, phydbl *d)
{
  int j;
  ludcmp(a,n,d);
  For(j,n) *d *= a[j][j];
}



int ludcmp_1D(phydbl *a, int n, phydbl *d)
{
   int i,imax,j,k;
   phydbl big,dum,sum,temp;
   phydbl *vv;
   
   imax = 0;
   vv = (phydbl *)mCalloc(n,sizeof(phydbl));

   *d=1.0;
   for (i=0;i<n;i++) 
     {
       big=0.0;
       for (j=0;j<n;j++)
         if ((temp=FABS(a[i*n+j])) > big) big=temp;
       if (FABS(big) < SMALL) Exit("\n. Singular matrix in routine LUDCMP");
       vv[i]=1.0/big;
     }
   for (j=0;j<n;j++) 
     {
       for (i=0;i<j;i++) 
	 {
	   sum=a[i*n+j];
	   for (k=0;k<i;k++) sum -= a[i*n+k]*a[k*n+j];
	   a[i*n+j]=sum;
	 }
      big=0.0;
      for (i=j;i<n;i++) {
	sum=a[i*n+j];
	for (k=0;k<j;k++)
	  sum -= a[i*n+k]*a[k*n+j];
	a[i*n+j]=sum;
	if ((dum=vv[i]*FABS(sum)) >= big) 
	  {
            big=dum;
            imax=i;
	  }
      }
      if (j != imax) 
	{
	  for (k=0;k<n;k++) 
	    {
	      dum=a[imax*n+k];
	      a[imax*n+k]=a[j*n+k];
	      a[j*n+k]=dum;
	    }
	  *d = -(*d);
	  vv[imax]=vv[j];
	}
      if (FABS(a[j*n+j]) < SMALL) a[j*n+j]=LUDCMP_TINY;
      if (j != n) {
	dum=1.0/(a[j*n+j]);
	for (i=j+1;i<n;i++) a[i*n+j] *= dum;
      }
     }
   free(vv);
   return(0);
}

void det_1D(phydbl *a, int n, phydbl *d)
{
  int j;
  ludcmp_1D(a,n,d);
  For(j,n) *d *= a[j*n+j];
}

/* Find L such that L.L' = A */
phydbl *Cholesky_Decomp(phydbl *A,  int dim)
{
  int i,j,k;
  phydbl sum;
  phydbl *L;

  L = (phydbl *)mCalloc(dim*dim,sizeof(phydbl));

  For(i,dim)
    {
      for(j=i;j<dim;j++)
	{
	  sum = A[j*dim+i];
	  for(k=0;k<i;k++) sum -= L[i*dim+k] * L[j*dim+k];

	  if(i == j)
	    {
	      if(sum < 1.E-20)
		{
		  PhyML_Printf("\n. sum=%G i=%d j=%d",sum,i,j);
		  PhyML_Printf("\n. Numerical precision issue detected...");
		  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
		  Warn_And_Exit("");
		}
	      L[j*dim+i] = SQRT(sum);
	    }

	  else L[j*dim+i] = sum / (L[i*dim+i]);

 	}
    }

  return L;
}

int EigenQ(phydbl *Q, phydbl *Qbuff, phydbl *pi, int n, phydbl *R, phydbl *Ri, phydbl *U, phydbl *Ui, phydbl *V, phydbl *space)
{
  phydbl scalar=1.0;
  int n_iter=0, result=0, i, nn=n*n;
  
  For(i,nn) Qbuff[i] = Q[i];
  
  if(!Eigen(1, Qbuff, n, R, Ri, U, Ui, space))
  {
    /* compute inverse(Vr) into Vi */
    For (i,nn) V[i] = U[i];
    
    while(!Matinv(V, n, n))
    {
      PhyML_Printf("\n. Trying Q<-Q*scalar and then Root<-Root/scalar to fix this...\n");
      scalar += scalar / 3.;
      For(i,nn) Qbuff[i] = Q[i];
      For(i,nn) Qbuff[i]*= scalar;
      result = Eigen(1, Qbuff, n, R, Ri, U, Ui, space);
      if (result == -1)
	Exit("\n. Eigenvalues/vectors computation did not converge : computation cancelled\n");
      else if (result == 1)
	Exit("\n. phycomplex eigenvalues/vectors : computation cancelled\n");
      For (i,nn) V[i] = U[i];
      n_iter++;
      if(n_iter > 10) Exit("\n. Cannot work out eigen vectors\n");
    };
    
    if(scalar!=1.0) For(i,n) R[i] /= scalar;
  }
  else
  {
    PhyML_Printf("\n. Eigenvalues/vectors computation does not converge : computation cancelled");
    Warn_And_Exit("\n");
  }
  return 0;
}

//!< Eigen routines below were EXTRACTED from PAML 4.4c 2010 |
//!< All the credit given to Prof. Yang                      V
int EigenQREV (double Q[], double pi[], int n, double Root[], double U[], double V[], double spacesqrtpi[])
{
/* 
   This finds the eigen solution of the rate matrix Q for a time-reversible 
   Markov process, using the algorithm for a real symmetric matrix.
   Rate matrix Q = S * diag{pi} = U * diag{Root} * V, 
   where S is symmetrical, all elements of pi are positive, and U*V = I.
   space[n] is for storing sqrt(pi).

   [U 0] [Q_0 0] [U^-1 0]    [Root  0]
   [0 I] [0   0] [0    I]  = [0     0]

   Ziheng Yang, 25 December 2001 (ref is CME/eigenQ.pdf)
*/
   int i,j, inew, jnew, nnew, status;
   double *pi_sqrt=spacesqrtpi, small=1e-100;

   for(j=0,nnew=0; j<n; j++)
      if(pi[j]>small)
         pi_sqrt[nnew++] = sqrt(pi[j]);

   /* store in U the symmetrical matrix S = sqrt(D) * Q * sqrt(-D) */

   if(nnew==n) {
      for(i=0; i<n; i++)
         for(j=0,U[i*n+i] = Q[i*n+i]; j<i; j++)
            U[i*n+j] = U[j*n+i] = (Q[i*n+j] * pi_sqrt[i]/pi_sqrt[j]);

      status = eigenRealSym(U, n, Root, V);
      for(i=0; i<n; i++) for(j=0; j<n; j++)  V[i*n+j] = U[j*n+i] * pi_sqrt[j];
      for(i=0; i<n; i++) for(j=0; j<n; j++)  U[i*n+j] /= pi_sqrt[i];
   }
   else {
      for(i=0,inew=0; i<n; i++) {
         if(pi[i]>small) {
            for(j=0,jnew=0; j<i; j++) 
               if(pi[j]>small) {
                  U[inew*nnew+jnew] = U[jnew*nnew+inew] 
                                    = Q[i*n+j] * pi_sqrt[inew]/pi_sqrt[jnew];
                  jnew++;
               }
            U[inew*nnew+inew] = Q[i*n+i];
            inew++;
         }
      }

      status=eigenRealSym(U, nnew, Root, V);

      for(i=n-1,inew=nnew-1; i>=0; i--)   /* construct Root */
         Root[i] = (pi[i]>small ? Root[inew--] : 0);
      for(i=n-1,inew=nnew-1; i>=0; i--) {  /* construct V */
         if(pi[i]>small) {
            for(j=n-1,jnew=nnew-1; j>=0; j--)
               if(pi[j]>small) {
                  V[i*n+j] = U[jnew*nnew+inew]*pi_sqrt[jnew];
                  jnew--;
               }
               else 
                  V[i*n+j] = (i==j);
            inew--;
         }
         else 
            for(j=0; j<n; j++)  V[i*n+j] = (i==j);
      }
      for(i=n-1,inew=nnew-1; i>=0; i--) {  /* construct U */
         if(pi[i]>small) {
            for(j=n-1,jnew=nnew-1;j>=0;j--)
               if(pi[j]>small) {
                  U[i*n+j] = U[inew*nnew+jnew]/pi_sqrt[inew];
                  jnew--;
               }
               else 
                  U[i*n+j] = (i==j);
            inew--;
         }
         else 
            for(j=0;j<n;j++)
               U[i*n+j] = (i==j);
      }
   }

/*   This routine works on P(t) as well as Q. */
/*
   if(fabs(Root[0])>1e-10 && noisy) printf("Root[0] = %.5e\n",Root[0]);
   Root[0]=0; 
*/
   return(status);
}

int eigenRealSym(double A[], int n, double Root[], double work[])
{
/* This finds the eigen solution of a real symmetrical matrix A[n*n].  In return, 
   A has the right vectors and Root has the eigenvalues. 
   work[n] is the working space.
   The matrix is first reduced to a tridiagonal matrix using HouseholderRealSym(), 
   and then using the QL algorithm with implicit shifts.  

   Adapted from routine tqli in Numerical Recipes in C, with reference to LAPACK
   Ziheng Yang, 23 May 2001
*/
   int status=0;
   HouseholderRealSym(A, n, Root, work);
   status = EigenTridagQLImplicit(Root, work, n, A);
   EigenSort(Root, A, n);

   return(status);
}

void EigenSort(double d[], double U[], int n)
{
/* this sorts the eigen values d[] and rearrange the (right) eigen vectors U[]
*/
   int k,j,i;
   double p;

   for (i=0;i<n-1;i++) {
      p=d[k=i];
      for (j=i+1;j<n;j++)
         if (d[j] >= p) p=d[k=j];
      if (k != i) {
         d[k]=d[i];
         d[i]=p;
         for (j=0;j<n;j++) {
            p=U[j*n+i];
            U[j*n+i]=U[j*n+k];
            U[j*n+k]=p;
         }
      }
   }
}

void HouseholderRealSym(double a[], int n, double d[], double e[])
{
/* This uses HouseholderRealSym transformation to reduce a real symmetrical matrix 
   a[n*n] into a tridiagonal matrix represented by d and e.
   d[] is the diagonal (eigends), and e[] the off-diagonal.
*/
   int m,k,j,i;
   double scale,hh,h,g,f;

   for (i=n-1;i>=1;i--) {
      m=i-1;
      h=scale=0;
      if (m > 0) {
         for (k=0;k<=m;k++)
            scale += fabs(a[i*n+k]);/***********************************************************
*  This eigen() routine works for eigenvalue/vector analysis
*         for real general square matrix A
*         A will be destroyed
*         rr,ri are vectors containing eigenvalues
*         vr,vi are matrices containing (right) eigenvectors
*
*              A*[vr+vi*i] = [vr+vi*i] * diag{rr+ri*i}
*
*  Algorithm: Handbook for Automatic Computation, vol 2
*             by Wilkinson and Reinsch, 1971
*             most of source codes were taken from a public domain
*             solftware called MATCALC.
*  Credits:   to the authors of MATCALC
*
*  return     -1 not converged
*              0 no phycomplex eigenvalues/vectors
*              1 phycomplex eigenvalues/vectors
*  Tianlin Wang at University of Illinois
*  Thu May  6 15:22:31 CDT 1993
***************************************************************/
         if (scale == 0)
            e[i]=a[i*n+m];
         else {
            for (k=0;k<=m;k++) {
               a[i*n+k] /= scale;
               h += a[i*n+k]*a[i*n+k];
            }
            f=a[i*n+m];
            g=(f >= 0 ? -sqrt(h) : sqrt(h));
            e[i]=scale*g;
            h -= f*g;
            a[i*n+m]=f-g;
            f=0;
            for (j=0;j<=m;j++) {
               a[j*n+i]=a[i*n+j]/h;
               g=0;
               for (k=0;k<=j;k++)
                  g += a[j*n+k]*a[i*n+k];
               for (k=j+1;k<=m;k++)
                  g += a[k*n+j]*a[i*n+k];
               e[j]=g/h;
               f += e[j]*a[i*n+j];
            }
            hh=f/(h*2);
            for (j=0;j<=m;j++) {
               f=a[i*n+j];
               e[j]=g=e[j]-hh*f;
               for (k=0;k<=j;k++)
                  a[j*n+k] -= (f*e[k]+g*a[i*n+k]);
            }
         }
      } 
      else
         e[i]=a[i*n+m];
      d[i]=h;
   }
   d[0]=e[0]=0;

   /* Get eigenvectors */
   for (i=0;i<n;i++) {
      m=i-1;
      if (d[i]) {
         for (j=0;j<=m;j++) {
            g=0;
            for (k=0;k<=m;k++)
               g += a[i*n+k]*a[k*n+j];
            for (k=0;k<=m;k++)
               a[k*n+j] -= g*a[k*n+i];
         }
      }
      d[i]=a[i*n+i];
      a[i*n+i]=1;
      for (j=0;j<=m;j++) a[j*n+i]=a[i*n+j]=0;
   }
}

int EigenTridagQLImplicit(double d[], double e[], int n, double z[])
{
/* This finds the eigen solution of a tridiagonal matrix represented by d and e.  
   d[] is the diagonal (eigenvalues), e[] is the off-diagonal
   z[n*n]: as input should have the identity matrix to get the eigen solution of the 
   tridiagonal matrix, or the output from HouseholderRealSym() to get the 
   eigen solution to the original real symmetric matrix.
   z[n*n]: has the orthogonal matrix as output

   Adapted from routine tqli in Numerical Recipes in C, with reference to
   LAPACK fortran code.
   Ziheng Yang, May 2001
*/
   int m,j,iter,niter=30, status=0, i,k;
   double s,r,p,g,f,dd,c,b, aa,bb;

   for (i=1;i<n;i++) e[i-1]=e[i];  e[n-1]=0;
   for (j=0;j<n;j++) {
      iter=0;
      do {
         for (m=j;m<n-1;m++) {
            dd=fabs(d[m])+fabs(d[m+1]);
            if (fabs(e[m])+dd == dd) break;  /* ??? */
         }
         if (m != j) {
            if (iter++ == niter) {
               status=-1;
               break;
            }
            g=(d[j+1]-d[j])/(2*e[j]);

            /* r=pythag(g,1); */

            if((aa=fabs(g))>1)  r=aa*sqrt(1+1/(g*g));
            else                r=sqrt(1+g*g);

            g=d[m]-d[j]+e[j]/(g+SIGN(r,g));
            s=c=1;
            p=0;
            for (i=m-1;i>=j;i--) {
               f=s*e[i];
               b=c*e[i];

               /*  r=pythag(f,g);  */
               aa=fabs(f); bb=fabs(g);
               if(aa>bb)       { bb/=aa;  r=aa*sqrt(1+bb*bb); }
               else if(bb==0)             r=0;
               else            { aa/=bb;  r=bb*sqrt(1+aa*aa); }

               e[i+1]=r;
               if (r == 0) {
                  d[i+1] -= p;
                  e[m]=0;
                  break;
               }
               s=f/r;
               c=g/r;
               g=d[i+1]-p;
               r=(d[i]-g)*s+2*c*b;
               d[i+1]=g+(p=s*r);
               g=c*r-b;
               for (k=0;k<n;k++) {
                  f=z[k*n+i+1];
                  z[k*n+i+1]=s*z[k*n+i]+c*f;
                  z[k*n+i]=c*z[k*n+i]-s*f;
               }
            }
            if (r == 0 && i >= j) continue;
            d[j]-=p; e[j]=g; e[m]=0;
         }
      } while (m != j);
   }
   return(status);
}


