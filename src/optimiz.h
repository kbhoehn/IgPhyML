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

#ifndef OPTIMIZ_H
#define OPTIMIZ_H

#include "utilities.h"
#include "lk.h"
#include "free.h"
#include "models.h"
/* unused ken 5/1
#include "mg.h"
*/
void Set_Ancestors(t_node *a, t_node *d, t_tree *tree); //added by ken
void Set_Ancestors_Root(t_tree *tree);
void      Optimiz_Ext_Br(t_tree *tree);
void      Optimize_Alpha(t_tree *tree);
void      Optimize_Kappa(t_tree *tree);
void      Optimize_Lambda(t_tree *tree);
void      Optimize_Param_Parall(t_tree *tree);
phydbl    Optimize_Branch_Quad(t_tree *tree, calign *cdata, t_edge *b_fcus);
void      Optimize_After_Hide(t_tree *tree, calign *cdata, t_node *h);
void      Round_Optimize(t_tree *tree, calign *data, int n_round_max);
int       Dist_Seq_Brak(phydbl *ax, phydbl *bx, phydbl *cx, 
			phydbl *fa, phydbl *fb, phydbl *fc, 
			calign *data, int num1, int num2, model *mod);
phydbl    Dist_Seq_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
			 phydbl *xmin, calign *data, 
			 int num1, int num2, model *mod);
phydbl    Kappa_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
		       phydbl *xmin, t_tree *tree, calign *cdata);
phydbl    Lambda_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
			phydbl *xmin, t_tree *tree, calign *cdata);
phydbl    Alpha_Golden_Br_Opt(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
			      phydbl *xmin, t_tree *tree, calign *cdata, 
			      int n_opt, phydbl *init_l);
phydbl    Alpha_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol,phydbl *xmin, 
		       t_tree *tree, calign *cdata);
phydbl    Br_Len_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
			phydbl *xmin, t_edge *b_fcus, t_tree *tree);
phydbl    Br_Len_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
		       t_edge *b_fcus, t_tree *tree, int n_iter_max, int quickdirty);
int       Br_Len_Brak(phydbl *ax, phydbl *bx, phydbl *cx, 
		      phydbl *fa, phydbl *fb, phydbl *fc, 
		      t_edge *b_fcus, t_tree *tree);
phydbl    Optimize_Path_Length(model *mod, calign *cdata, t_edge *a, 
			       int lra, t_edge *b, int lrb, phydbl i_len);
void      Optimize_Param_Serie(t_node *a, t_node *d, t_edge *b_fcus, t_tree *tree, 
			       calign *cdata, int n_passes);
phydbl    Optimize_Dist(model *mod, phydbl init, calign *twoseqs);
phydbl    Pinvar_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
			phydbl *xmin, t_tree *tree, calign *cdata, int n_iter_max);
void      Optimize_Pinvar(t_tree *tree);
int       Lambda_Brak(phydbl *ax, phydbl *bx, phydbl *cx, 
		      phydbl *fa, phydbl *fb, phydbl *fc, 
		      t_tree *tree);
int       Kappa_Brak(phydbl *ax, phydbl *bx, phydbl *cx, 
		      phydbl *fa, phydbl *fb, phydbl *fc, 
		      t_tree *tree);
int       Alpha_Brak(phydbl *ax, phydbl *bx, phydbl *cx, 
		      phydbl *fa, phydbl *fb, phydbl *fc, 
		      t_tree *tree);
int       Pinvar_Brak(phydbl *ax, phydbl *bx, phydbl *cx, 
		      phydbl *fa, phydbl *fb, phydbl *fc, 
		      t_tree *tree);
void Optimiz_All_Free_Param(t_tree *tree, int verbose);
void      Optimiz_RRparam_GTR(t_tree *tree, int num_param);
phydbl    RRparam_GTR_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
		   	     phydbl *xmin, t_tree *tree, calign *cdata, phydbl *param, int n_iter_max);

int Powell_GTR_Param(t_tree *tree, phydbl *p, int n, phydbl ftol);
phydbl Linmin_GTR_Param(t_tree *tree,phydbl *p, phydbl *xi, int n);
phydbl F1dim(t_tree *tree, phydbl x, phydbl *p, phydbl *xi, phydbl n);
int Mnbrak_1dim(phydbl *ax, phydbl *bx, phydbl *cx, 
		phydbl *fa, phydbl *fb, phydbl *fc,
		t_tree *tree,
		phydbl *p,  phydbl *xi, phydbl n);
phydbl Brent_1dim(phydbl ax, phydbl bx, phydbl cx, 
		  phydbl tol, phydbl *xmin,
		  t_tree *tree,
		  phydbl *p, phydbl *xi, phydbl n);

int Min_With_Derivatives(t_tree *tree, phydbl *p, int n, phydbl ftol, phydbl step_size, 
			 phydbl (*func) (), void (*dfunc)(), phydbl (*linmin)());
void BFGS(t_tree *tree, 
	  phydbl *p, 
	  int n, 
	  phydbl gtol, 
	  phydbl step_size,
	  phydbl(*func)(t_tree *tree), 
	  int(*dfunc)(t_tree *tree,phydbl *param,int n_param,phydbl stepsize,phydbl(*func)(t_tree *tree),phydbl *derivatives), 
	  int(*lnsrch)(t_tree *tree, int n, phydbl *xold, phydbl fold,phydbl *g, phydbl *p, phydbl *x,phydbl *f, phydbl stpmax, int *check),
	  int *failed);
int Lnsrch_RR_Param(t_tree *tree, int n, phydbl *xold, phydbl fold, phydbl *g, phydbl *p, phydbl *x,
		     phydbl *f, phydbl stpmax, int *check);
void Optimize_Single_Param_Generic(t_tree *tree, phydbl *param, phydbl lim_inf, phydbl lim_sup, phydbl tol, int n_max_iter, int quickdirty);
int Generic_Brak(phydbl *param,
		 phydbl *ax, phydbl *bx, phydbl *cx, 
		 phydbl *fa, phydbl *fb, phydbl *fc,
		 phydbl lim_inf, phydbl lim_sup,
		 t_tree *tree);
phydbl Generic_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
		     phydbl *xmin, t_tree *tree, int n_iter_max,int quickdirty);
void Optimize_Br_Len_Serie(t_node *a, t_node *d, t_edge *b_fcus, t_tree *tree,calign *cdata);
int Lnsrch_Nucleotide_Frequencies(t_tree *tree, int n, phydbl *xold, 
				   phydbl fold, phydbl *g, phydbl *p, phydbl *x,
				   phydbl *f, phydbl stpmax, int *check);

void Optimize_Global_Rate(t_tree *tree);
phydbl Br_Len_Brent_Default(t_edge *b_fcus, t_tree *tree);

void EM_Dist(model *mod, calign *data);
phydbl Dist_F_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, int n_iter_max, 
		    phydbl *param, phydbl *F, model *mod);
int Dist_F_Brak(phydbl *ax, phydbl *bx, phydbl *cx, phydbl *F, phydbl *param, model *mod);
void Opt_Dist_F(phydbl *dist, phydbl *F, model *mod);
phydbl Missing_Dist_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, int n_iter_max, 
			  int x, int y, matrix *mat);
int Missing_Dist_Brak(phydbl *ax, phydbl *bx, phydbl *cx, int x, int y, matrix *mat);
void Opt_Missing_Dist(int x, int y, matrix *mat);
int Optimiz_Alpha_And_Pinv(t_tree *tree);
int Lnsrch_RR_Cov_Param(t_tree *tree, int n, phydbl *xold, phydbl fold, 
			 phydbl *g, phydbl *p, phydbl *x,
			 phydbl *f, phydbl stpmax, int *check);
phydbl Node_Time_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol,
		       t_node *anc, t_node *des, t_tree *tree, int n_iter_max);
phydbl Time_Stamps_Mult_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol,
			      t_tree *tree, int n_iter_max);
phydbl Branch_Rate_Shape_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
			       phydbl *xmin, t_tree *tree, int n_iter_max);
phydbl Node_Time_Brent_Fixed_Br_Len(phydbl ax, phydbl bx, phydbl cx, phydbl tol,
				    t_node *n, t_tree *tree, int n_iter_max);

phydbl Generic_Brent_Lk(phydbl *param, phydbl ax, phydbl cx, phydbl tol, 
			int n_iter_max, int quickdirty,
			phydbl (*obj_func)(t_edge *,t_tree *,supert_tree *), 
			t_edge *branch, t_tree *tree, supert_tree *stree);
phydbl Optwrap_Lk(t_edge *b, t_tree *tree, supert_tree *stree);
phydbl Optwrap_Lk_At_Given_Edge(t_edge *b, t_tree *tree, supert_tree *stree);
phydbl Optwrap_Part_Lk_At_Given_Edge(t_edge *b, t_tree *tree, supert_tree *stree);
phydbl Optwrap_Part_Lk(t_edge *b, t_tree *tree, supert_tree *stree);


void      Round_Optimize_New(t_tree *tree, calign *data, int n_round_max); //!Added by Marcelo.

#define ITMAX 200 
#define EPS   3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0
#define SQR(a) ((sqrarg=(a)) < SMALL ? 0.0 : sqrarg*sqrarg)

phydbl Br_Len_Brent_Codon_Pairwise(phydbl ax, phydbl bx, phydbl cx, phydbl tol, phydbl *b_fcus, phydbl *Pij, phydbl *pi, eigen *eigenStruct, calign *data, int ns, int n_iter_max, int quickdirty, phydbl *uexpt, phydbl *expt); //!< Added by Marcelo.

int gradientB (int n, phydbl x[], phydbl f0, phydbl g[], t_tree * tree, phydbl space[], int xmark[]);
extern phydbl SIZEp;
extern int noisy, Iround,NFunCall;
extern int AlwaysCenter;
extern phydbl Small_Diff;
#define min2(a,b) ((a)<(b)?(a):(b))
#define square(a) ((a)*(a))
#define max2(a,b) ((a)>(b)?(a):(b))
#define square(a) ((a)*(a))
#define CDFGamma(x,alpha,beta) IncompleteGamma((beta)*(x),alpha,LnGamma(alpha))
phydbl innerp (phydbl x[], phydbl y[], int n);
int xtoy (phydbl x[], phydbl y[], int n);
phydbl fun_LineSearch (phydbl t, t_tree *tree, phydbl x0[], phydbl p[], phydbl x[], int n);
phydbl LineSearch2 (t_tree * tree, phydbl *f, phydbl x0[], phydbl p[], phydbl step, phydbl limit, phydbl e, phydbl space[], int n);
int zero (phydbl x[], int n);
int identity (phydbl x[], int n);
phydbl norm (phydbl x[], int n);
int H_end (phydbl x0[], phydbl x1[], phydbl f0, phydbl f1,phydbl e1, phydbl e2, int n);
phydbl distance (phydbl x[], phydbl y[], int n);
int BFGS_from_CODEML (phydbl *f, t_tree *tree, phydbl *x, phydbl xb[120][2], phydbl space[], phydbl e, int n);
phydbl Sum (phydbl x[], int n);
#endif

