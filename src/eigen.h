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

#ifndef EIGEN_H
#define EIGEN_H

#include "utilities.h"
#include "free.h"

#ifdef RWRAPPER
#include <R.h>
#endif

int ludcmp_1D(phydbl *a, int n, phydbl *d);
void det_1D(phydbl *a, int n, phydbl *d);

int Eigen(int job, phydbl *A, int n, phydbl *rr, phydbl *ri,
          phydbl *vr, phydbl *vi, phydbl *w);
void balance(phydbl *mat, int n, int *low, int *hi, phydbl *scale);
void unbalance(int n, phydbl *vr, phydbl *vi, int low, int hi,
               phydbl *scale);
int realeig(int job, phydbl *mat, int n,int low, int hi, phydbl *valr,
            phydbl *vali, phydbl *vr, phydbl *vi);
void elemhess(int job, phydbl *mat, int n, int low, int hi, 
            phydbl *vr, phydbl *vi, int *work);

int ludcmp(phydbl **a, int n, phydbl *d);
void det(phydbl **a, int n, phydbl *d);

/* phycomplex functions */

typedef struct { phydbl re, im; } phycomplex;
#define csize(a) (FABS(a.re)+FABS(a.im))

phycomplex compl (phydbl re,phydbl im);
phycomplex _conj (phycomplex a);
phycomplex cplus (phycomplex a, phycomplex b);
phycomplex cminus (phycomplex a, phycomplex b);
phycomplex cby (phycomplex a, phycomplex b);
phycomplex cdiv (phycomplex a,phycomplex b);
/* phycomplex local_cexp (phycomplex a); */
phycomplex cfactor (phycomplex x, phydbl a);
int cxtoy (phycomplex *x, phycomplex *y, int n);
int cmatby (phycomplex *a, phycomplex *b, phycomplex *c, int n,int m,int k);
int cmatout (FILE * fout, phycomplex *x, int n, int m);
int cmatinv( phycomplex *x, int n, int m, phydbl *space);
phydbl *Cholesky_Decomp(phydbl *A,  int dim);

int EigenQ(phydbl *Q, phydbl *Qbuff, phydbl *pi, int n, phydbl *R, phydbl *Ri, phydbl *U, phydbl *Ui, phydbl *V, phydbl *space); //!< Added by Marcelo.

//!< Eigen routines below were EXTRACTED from PAML 4.4c 2010 |
//!<                                                         V
void HouseholderRealSym(double a[], int n, double d[], double e[]);
int EigenTridagQLImplicit(double d[], double e[], int n, double z[]);
void EigenSort(double d[], double U[], int n);
int EigenQREV (double Q[], double pi[], int n, double Root[], double U[], double V[], double spacesqrtpi[]);
int eigenRealSym(double A[], int n, double Root[], double work[]);

#endif

