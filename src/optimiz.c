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

#include "optimiz.h"

/*********************************************************/

phydbl Br_Len_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol,
		    t_edge *b_fcus, t_tree *tree, int n_iter_max, int quickdirty)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL, init_lnL;

  /*if(tree->mod->s_opt->opt_topo){
	  if(b_fcus->anc_node->num != tree->mod->startnode){
		  Fill_UPP_single(tree,b_fcus);
	  }else{
		  Fill_UPP_root(tree,b_fcus);
	  }
  }*/

  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  old_lnL = UNLIKELY;
  b_fcus->l = FABS(bx);
  fw=fv=fx=fu=-Lk_At_Given_Edge(b_fcus,tree);
  init_lnL = -fw;

/*   PhyML_Printf("\n. INIT BRENT t_edge %3d l=%f lnL=%20f",b_fcus->num,b_fcus->l,fu); */

  for(iter=1;iter<=BRENT_ITMAX;iter++)
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*FABS(x)+BRENT_ZEPS);

      if((tree->c_lnL > init_lnL + tol) && (quickdirty))
	{
	  b_fcus->l = x;
	  Lk_At_Given_Edge(b_fcus,tree);
/* 	  PhyML_Printf("\n> iter=%3d max=%3d v=%f lnL=%f init_lnL=%f tol=%f",iter,n_iter_max,(*xmin),tree->c_lnL,init_lnL,tol); */
/* 	  Exit("\n"); */
	  return tree->c_lnL;	  
	}

/*       if(((FABS(tree->c_lnL-old_lnL) < tol) && (tree->c_lnL > init_lnL - tol)) || (iter > n_iter_max - 1)) */
      if((FABS(tree->c_lnL-old_lnL) < tol) || (iter > n_iter_max - 1))
	{
	  b_fcus->l=x;
	  Lk_At_Given_Edge(b_fcus,tree);
/* 	  PhyML_Printf("\n. iter=%3d max=%3d l=%f lnL=%f init_lnL=%f",iter,n_iter_max,b_fcus->l,tree->c_lnL,init_lnL); */
/* 	  Exit("\n"); */
	  return tree->c_lnL;
	}
      
      if(FABS(e) > tol1)
	{
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if(q > 0.0) p = -p;
	  q=FABS(q);
	  etemp=e;
	  e=d;
	  if(FABS(p) >= FABS(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	  else{
	    d=p/q;
	    u=x+d;
	    if (u-a < tol2 || b-u < tol2)
	      d=SIGN(tol1,xm-x);
	  }
	}
      else
	{
	  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	}
      u=(FABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));
/*       if(u<BL_MIN) u = BL_MIN; */
      b_fcus->l=FABS(u);
      old_lnL = tree->c_lnL;
      fu=-Lk_At_Given_Edge(b_fcus,tree);

/*       PhyML_Printf("\n. BRENT t_edge %3d l=%f lnL=%20f iter=%3d",b_fcus->num,b_fcus->l,fu,iter); */

/*       if(fu <= fx) */
      if(fu < fx)
	{
	  if(u > x) a=x; else b=x;
/* 	  if(u >= x) a=x; else b=x; */
	  SHFT(v,w,x,u)
	  SHFT(fv,fw,fx,fu)
	}
      else
	{
	  if (u < x) a=u; else b=u;
/* 	  if (fu <= fw || w == x) */
	  if (fu < fw || FABS(w-x) < SMALL)
	    {
	      v=w;
	      w=u;
	      fv=fw;
	      fw=fu;
	    }
	  else if (fu < fv || FABS(v-x) < SMALL || FABS(v-w) < SMALL)
/* 	  else if (fu <= fv || v == x || v == w)  */
	    {
	      v=u;
	      fv=fu;
	    }
	}
    }

  if(iter > BRENT_ITMAX) PhyML_Printf("\n. Too many iterations in BRENT (%d) (%f)",iter,b_fcus->l);
  return(-1);
  /* Not Reached ??  *xmin=x;   */
  /* Not Reached ??  return fx; */
}

/*********************************************************/
//Modified by Ken
void Round_Optimize(t_tree *tree, calign *data, int n_round_max)
{
  int n_round,each;
  phydbl lk_old, lk_new, tol;
  t_node *root;

  lk_new = tree->c_lnL;
  lk_old = UNLIKELY;
  n_round = 0;
  each = 0;
  tol = 1.e-2;
  root = tree->noeud[tree->mod->startnode];
  if(tree->mod->whichrealmodel == HLP17){
	  root = tree->noeud[tree->mod->startnode];
  }
  

  while(n_round < n_round_max)
    {
	  if(tree->mod->whichrealmodel != HLP17){
		  (!((n_round+2)%2))?(root=tree->noeud[0]):(root=tree->noeud[tree->n_otu-1]);
	  }
      
      if(tree->has_branch_lengths)//!< Added by Marcelo ... in this case opt parameters first
      {
	if(!each)
	{
	  each = 1;
	  Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));

	}
	
	if(tree->mod->s_opt->opt_bl){
		Lk(tree);
		if(tree->mod->whichrealmodel == HLP17){Get_UPP(root, root->v[0], tree);}
	//	Get_Lhood(tree->noeud[tree->mod->startnode], tree->noeud[tree->mod->startnode]->v[0], tree);
		Optimize_Br_Len_Serie(root,root->v[0],root->b[0],tree,data);
	}


	lk_new = tree->c_lnL;
	if((tree->mod->s_opt->print) && (!tree->io->quiet)) Print_Lk(tree,"[Branch lengths     ]");
      }
      else
      {
	if(tree->mod->s_opt->opt_bl){
		tree->both_sides = 1;
     	Lk(tree);
     	if(tree->mod->whichrealmodel == HLP17){Get_UPP(root, root->v[0], tree);}
	//	Get_Lhood(tree->noeud[tree->mod->startnode], tree->noeud[tree->mod->startnode]->v[0], tree);
		Optimize_Br_Len_Serie(root,root->v[0],root->b[0],tree,data);
	}
	

	if((tree->mod->s_opt->print) && (!tree->io->quiet)) Print_Lk(tree,"[Branch lengths     ]");

	if(!each)
	{
	  each = 1;
	  Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));
	}
      }
      Lk(tree);
      if(tree->mod->whichrealmodel == HLP17){
    	  Get_UPP(root, root->v[0], tree);
    	  Get_Lhood(root,root->v[0],tree);
      }
      tree->both_sides = 1;
      Lk(tree);
      lk_new = tree->c_lnL;
    //  printf("new tree lk2: %lf\n",lk_new);
      if(lk_new < lk_old - tree->mod->s_opt->min_diff_lk_global){
    	  printf("Old: %lf, New: %lf\n",lk_old,lk_new);
    	  Exit("\n. Optimisation failed ! (Round_Optimize)\n");
      }
      if(FABS(lk_new - lk_old) < tree->mod->s_opt->min_diff_lk_global)  break;
      else lk_old  = lk_new;
      n_round++;
      each--;
    }
  
  Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));
}

/*********************************************************/
//Edited by Ken
void Optimize_Br_Len_Serie(t_node *a, t_node *d, t_edge *b_fcus, t_tree *tree, calign *cdata)
{
  int i;
  phydbl l_infa,l_max,l_infb;
  phydbl lk_init;

//removed by Ken 8/9
//  tree->both_sides = 1;
//  Lk(tree);
  lk_init = tree->c_lnL;
  
  l_infa = l_max  = l_infb = BL_MIN;
  
  l_infa = BL_MAX;
  l_max  = b_fcus->l;
  l_infb = BL_MIN;



  /*if(a->num == tree->mod->startnode){
	   printf("doing root adjustment6\n");
	   Update_P_Lk(tree,d->anc_edge,d);
	   Fill_UPP_root(tree,d->anc_edge,d,a);
  }*/


  Br_Len_Brent(l_infa,l_max,l_infb,
	       tree->mod->s_opt->min_diff_lk_local,
	       b_fcus,tree,
	       tree->mod->s_opt->brent_it_max,
	       tree->mod->s_opt->quickdirty);

  	//Added to catch potential issues with branch optimization

  if(tree->mod->whichrealmodel == HLP17){
   if(tree->c_lnL < lk_init - tree->mod->s_opt->min_diff_lk_local){

		  PhyML_Printf("\n. %f %f %f %f %d %d %d\n",l_infa,l_max,l_infb,b_fcus->l,b_fcus->num,a->num,d->num);
		  PhyML_Printf("\n. %f -- %f\n",lk_init,tree->c_lnL);
		  Warn_And_Exit("\n. Err. in Optimize_Br_Len_Serie\n");
   }
  }

    
  if(a->num == tree->mod->startnode){
	//   printf("doing root adjustment7\n");
	//   Update_P_Lk(tree,d->anc_edge,d);
	  if(tree->mod->whichrealmodel == HLP17){Fill_UPP_root(tree,d->anc_edge);}
  }
  //printf("did root adjustment\n");

  if(d->tax) return;
   else For(i,3) if(d->v[i] != a)
   {
	   Update_P_Lk(tree,d->b[i],d);
	   if(tree->mod->whichrealmodel == HLP17){Fill_UPP_single(tree,d->b[i]);}
	   Optimize_Br_Len_Serie(d,d->v[i],d->b[i],tree,cdata);
   }
  
  For(i,3) if((d->v[i] == a) && !(d->v[i]->tax)) Update_P_Lk(tree,d->b[i],d);
}

/*********************************************************/

void Optimiz_Ext_Br(t_tree *tree)
{
  int i;
  t_edge *b;
  phydbl l_infa,l_max,l_infb;
  phydbl lk, lk_init,l_init;
  
  lk_init = tree->c_lnL;
  

  For(i,2*tree->n_otu-3)
    {
      b = tree->t_edges[i];
      if((b->left->tax) || (b->rght->tax))
	{

	  l_init = b->l;

/* 	  Fast_Br_Len(b,tree); */
/* 	  lk = Lk_At_Given_Edge(tree,b); */

	  l_infa = 10.*b->l;
	  l_max  = b->l;
	  l_infb = BL_MIN;

	  int br;
	  	      	  int n_edges=2*tree->n_otu-3;
	  	   //   	  For(br,n_edges) Update_PMat_At_Given_Edge(tree->t_edges[br],tree);
	//  Get_UPP(tree->noeud[tree->mod->startnode], tree->noeud[tree->mod->startnode]->v[0], tree);

	  lk = Br_Len_Brent(l_infa,l_max,l_infb,
			    tree->mod->s_opt->min_diff_lk_local,
			    b,tree,
			    tree->mod->s_opt->brent_it_max,
			    tree->mod->s_opt->quickdirty);

	  b->nni->best_l    = b->l;
	  b->nni->l0        = b->l;
	  b->nni->best_conf = 0;
	  b->l              = l_init;

	  if(tree->mod->whichrealmodel == HLP17){
	  Update_PMat_At_Given_Edge(b,tree);

	  }
	}
    }
  tree->c_lnL = lk_init; 
}

/*********************************************************/
phydbl gemin=1e-6; //!< Added by Marcelo.

//Modified by Ken to update h
void Optimiz_All_Free_Param(t_tree *tree, int verbose)
{

  if(tree->mod->io->datatype==CODON) //!< Added by Marcelo.
  {
    char s[100],r[100];
    int  init_both_sides, numParams = 0, i, n;
      
    if(tree->mod->s_opt->opt_method==optPAML)
    {
      phydbl x2min[120], x2minbound[120][2], fx, intf,newf, *space; //!< 120 is a very pessimist number of parametrs adopted to simplify the selection for optimiuazion.  
      
      init_both_sides  = tree->both_sides;
      tree->both_sides = 0;
      
      if ((tree->mod->whichmodel==GYECMS05  || //! Those models are updated only once ... need to clean up and make this more robust to consider the case where no parameters are optimized
	tree->mod->whichmodel==YAPECMS05 || 
	tree->mod->whichmodel==GYECMK07  || 
	tree->mod->whichmodel==YAPECMK07 ||
	tree->mod->whichmodel==GYECMS05F   ||
	tree->mod->whichmodel==GYECMK07F   ||
	tree->mod->whichmodel==YAPECMS05F  ||
	tree->mod->whichmodel==YAPECMK07F ||
	tree->mod->whichmodel==GYECMUSR ||
	tree->mod->whichmodel==GYECMUSRF ||
	tree->mod->whichmodel==MGECMUSRF ||
	tree->mod->whichmodel==MGECMUSR ||
	tree->mod->whichmodel==YAPECMUSR ||
	tree->mod->whichmodel==YAPECMUSRF ) && 
	tree->mod->s_opt->opt_state_freq==NO &&
	(tree->mod->n_catg>1 && tree->mod->n_w_catg==1)&&
	(tree->mod->pcaModel==0))
      {
	tree->mod->update_eigen = 0;
      }
      else
      {
	tree->mod->update_eigen = 1;
      }
      
      if(tree->io->heuristicExpm){ tree->io->expm=TAYLOR; tree->io->optParam=1; }
      
      if(tree->mod->s_opt->opt_kappa)
      {
	if(tree->mod->whichmodel!=GYECMK07WK  && tree->mod->whichmodel!=GYECMK07WKF && 
	  tree->mod->whichmodel!=GYECMS05WK  && tree->mod->whichmodel!=GYECMS05WKF && 
	  tree->mod->whichmodel!=MGECMK07WK  && tree->mod->whichmodel!=MGECMK07WKF && 
	  tree->mod->whichmodel!=MGECMS05WK  && tree->mod->whichmodel!=MGECMS05WKF && 
	  tree->mod->whichmodel!=YAPECMK07WK && tree->mod->whichmodel!=YAPECMK07WKF&& 
	  tree->mod->whichmodel!=YAPECMS05WK && tree->mod->whichmodel!=YAPECMS05WKF &&
	  tree->mod->whichmodel!=GYECMUSRWK && tree->mod->whichmodel!=GYECMUSRWKF &&
	  tree->mod->whichmodel!=MGECMUSRWKF && tree->mod->whichmodel!=YAPECMUSRWKF &&
	  tree->mod->whichmodel!=MGECMUSRWK && tree->mod->whichmodel!=YAPECMUSRWK)
	{
	  x2min[numParams++]               = tree->mod->kappa;
	  x2minbound[numParams-1][0]       = TREEBOUNDLOW;
	  x2minbound[numParams-1][1]       = TREEBOUNDHIGH;
	}
	else
	{
	  switch(tree->io->kappaECM)
	  {
	    case kap1: 
	      break;
	    case kap2:
	    case kap3: 
	    case kap6:  
	    {
	      x2min[numParams++]           = tree->mod->pkappa[0]; 
	      x2minbound[numParams-1][0]   = TREEBOUNDLOW;
	      x2minbound[numParams-1][1]   = TREEBOUNDHIGH;
	      break;
	    }
	    case kap4:
	    {
	      x2min[numParams++]           = tree->mod->pkappa[0]; 
	      x2minbound[numParams-1][0]   = TREEBOUNDLOW;
	      x2minbound[numParams-1][1]   = TREEBOUNDHIGH;
	      x2min[numParams++]           = tree->mod->pkappa[1];
	      x2minbound[numParams-1][0]   = TREEBOUNDLOW;
	      x2minbound[numParams-1][1]   = TREEBOUNDHIGH;
	      break;
	    }
	    case kap5:
	    {
	      For(i,tree->mod->nkappa-1)
	      {
		x2min[numParams++]         = tree->mod->unspkappa[i];
		x2minbound[numParams-1][0] = -99.00;
		x2minbound[numParams-1][1] = 99.00;
	      }
	      break;
	    }
	    default: Warn_And_Exit("Error in Kappa assignment.");
	    break;
	  }
	}
      }
      
      if(tree->mod->s_opt->opt_omega)
      {
	if(tree->mod->omegaSiteVar==DM0)
	{
		int omegai; //added by Ken 18/8/2016
		 for(omegai=0;omegai<tree->mod->nomega_part;omegai++){
		    x2min[numParams++]           = tree->mod->omega_part[omegai];
		    x2minbound[numParams-1][0]   = TREEBOUNDLOW*100; //changed by Ken 9/2/2017 due to underflow issues with highly polymorphic lineages
		    x2minbound[numParams-1][1]   = TREEBOUNDHIGH;

		  }
	}
	else if(tree->mod->omegaSiteVar==DMODELK)
	{	
	  For(i,tree->mod->n_w_catg)
	  {
	    x2min[numParams++]         = tree->mod->omegas[i];
	    x2minbound[numParams-1][0] = TREEBOUNDLOW;
	    x2minbound[numParams-1][1] = TREEBOUNDHIGH;
	  }
	  For(i,tree->mod->n_w_catg-1)
	  {
	    x2min[numParams++]         = tree->mod->prob_omegas_uns[i];
	    x2minbound[numParams-1][0] = -99.0;
	    x2minbound[numParams-1][1] = 99.0;
	  }
	}
	else if(tree->mod->omegaSiteVar==DGAMMAK)
	{
	  x2min[numParams++]           = tree->mod->alpha; 
	  x2minbound[numParams-1][0]   = 2e-3;
	  x2minbound[numParams-1][1]   = TREEBOUNDHIGH;
	  x2min[numParams++]           = tree->mod->beta; 
	  x2minbound[numParams-1][0]   = 2e-3;
	  x2minbound[numParams-1][1]   = TREEBOUNDHIGH;	
	}
      }
      if(tree->mod->opthotness == 1) //added by Kenneth Hoehn
      {
    		  int c;
    	  	  for(c=0;c<tree->mod->nhotness;c++){
    		  if(tree->mod->hoptindex[c] == 1){
    	  		  x2min[numParams++]           = tree->mod->hotness[c];
    	  	  	  x2minbound[numParams-1][0]   = -1; //new minimum value as of 12/June/2016
    	  	  	  x2minbound[numParams-1][1]   = TREEBOUNDHIGH;
    	  	 	 }
    	  	  }

      }
      
      if(tree->mod->s_opt->opt_state_freq)
      {
	switch(tree->mod->freq_model)
	{
	  case F1XSENSECODONS:
	  {
	    For(i,tree->mod->num_base_freq-1)
	    {
	      x2min[numParams++]         = tree->mod->pi_unscaled[i];
	      x2minbound[numParams-1][0] = -99.0;
	      x2minbound[numParams-1][1] = 99.0;
	    }
	    break;
	  }
	  case F1X4:
	  {
	    For(i,tree->mod->num_base_freq-1)
	    {
	      x2min[numParams++]         = tree->mod->uns_base_freq[i];
	      x2minbound[numParams-1][0] = -99.0;
	      x2minbound[numParams-1][1] = 99.0;
	    }	  
	    break;
	  }
	  case F3X4:
	  case CF3X4:
	  {
	    for(i=0;i<3;i++)
	    {
	      x2min[numParams++]         = tree->mod->uns_base_freq[i];
	      x2minbound[numParams-1][0] = -99.0;
	      x2minbound[numParams-1][1] = 99.0;
	    }
	    for(i=4;i<7;i++)
	    {
	      x2min[numParams++]         = tree->mod->uns_base_freq[i];
	      x2minbound[numParams-1][0] = -99.0;
	      x2minbound[numParams-1][1] = 99.0;
	    }
	    for(i=8;i<11;i++)
	    {
	      x2min[numParams++]         = tree->mod->uns_base_freq[i];
	      x2minbound[numParams-1][0] = -99.0;
	      x2minbound[numParams-1][1] = 99.0;
	    }
	    break;
	  }
	  default:
	    break;
	}
      }
      
      if(tree->mod->s_opt->opt_pinvar)
      {
	x2min[numParams++]               = tree->mod->pinvar;
	x2minbound[numParams-1][0]       = 0.0;
	x2minbound[numParams-1][1]       = 1.0;
      }
      
      if(tree->mod->n_catg>1 && tree->mod->n_w_catg==1 && tree->mod->s_opt->opt_alphaCD)
      { 
	
	x2min[numParams++]           = tree->mod->alpha;
	x2minbound[numParams-1][0]   = 2e-3;
	x2minbound[numParams-1][1]   = TREEBOUNDHIGH;
      }
      
      if(tree->mod->pcaModel==1)
      {
	
	For(i,tree->mod->npcs)
	{
	  x2min[numParams++]         = tree->mod->pcsC[i];
	  x2minbound[numParams-1][0] = (-1) * TREEBOUNDHIGH;
	  x2minbound[numParams-1][1] = TREEBOUNDHIGH ;
	}	
      }
      
      intf=fx=tree->c_lnL;
      
      if(numParams>0)
      {
	space = (phydbl *) mCalloc((20+5*numParams)*numParams, sizeof(phydbl));
	
	BFGS_from_CODEML (&fx, tree, x2min, x2minbound, space, gemin, numParams);
	
	newf=tree->c_lnL;
	gemin/=2;  if(FABS(newf-intf)<1) gemin/=2;
	if(FABS(newf-intf)<0.5)     gemin = min2(gemin,1e-3); 
	else if(FABS(newf-intf)>10) gemin = max2(gemin,0.1); 
	gemin = max2(gemin,1e-6);
	
	free(space);
      }
      
      tree->both_sides = init_both_sides;
      
      if(tree->io->heuristicExpm)
      {
	tree->mod->update_eigen = 1; 
	tree->io->optParam=0;
	tree->io->expm=EIGEN;
	Lk(tree);
	tree->mod->update_eigen = 0;  
      }
      else
      {
	tree->mod->update_eigen = 0;
	if(tree->both_sides) Lk(tree); /* Needed to update all partial likelihoods.*/  
      }
    }
    else
    {
      //!< use the old phyml method of parameter by parameter with brent.
      
      init_both_sides  = tree->both_sides;
      tree->both_sides = 0;
      For(n, tree->mod->s_opt->nBrentCycles)
      {
	if ((tree->mod->whichmodel==GYECMS05  || //! Those models are updated only once ... need to clean up and make this more robust to consider the case where no parameters are optimized
	  tree->mod->whichmodel==YAPECMS05 || 
	  tree->mod->whichmodel==GYECMK07  || 
	  tree->mod->whichmodel==YAPECMK07 ||
	  tree->mod->whichmodel==GYECMS05F   ||
	  tree->mod->whichmodel==GYECMK07F   ||
	  tree->mod->whichmodel==YAPECMS05F  ||
	  tree->mod->whichmodel==YAPECMK07F ||
	  tree->mod->whichmodel==GYECMUSR ||
	  tree->mod->whichmodel==GYECMUSRF ||
	  tree->mod->whichmodel==MGECMUSRF ||
	  tree->mod->whichmodel==MGECMUSR ||
	  tree->mod->whichmodel==YAPECMUSR ||
	  tree->mod->whichmodel==YAPECMUSRF ) && 
	  tree->mod->s_opt->opt_state_freq==NO &&
	  (tree->mod->n_catg>1 && tree->mod->n_w_catg==1)&&
	  (tree->mod->pcaModel==0))
	{
	  tree->mod->update_eigen = 0;
	}
	else
	{
	  tree->mod->update_eigen = 1;
	}
	
	if(tree->io->heuristicExpm){ tree->io->expm=TAYLOR; tree->io->optParam=1; }
	
	if(tree->mod->s_opt->opt_kappa)
	{
	  if(tree->mod->whichmodel!=GYECMK07WK  && tree->mod->whichmodel!=GYECMK07WKF && 
	    tree->mod->whichmodel!=GYECMS05WK  && tree->mod->whichmodel!=GYECMS05WKF && 
	    tree->mod->whichmodel!=MGECMK07WK  && tree->mod->whichmodel!=MGECMK07WKF && 
	    tree->mod->whichmodel!=MGECMS05WK  && tree->mod->whichmodel!=MGECMS05WKF && 
	    tree->mod->whichmodel!=YAPECMK07WK && tree->mod->whichmodel!=YAPECMK07WKF&& 
	    tree->mod->whichmodel!=YAPECMS05WK && tree->mod->whichmodel!=YAPECMS05WKF &&
	    tree->mod->whichmodel!=GYECMUSRWK && tree->mod->whichmodel!=GYECMUSRWKF &&
	    tree->mod->whichmodel!=MGECMUSRWKF && tree->mod->whichmodel!=YAPECMUSRWKF &&
	    tree->mod->whichmodel!=MGECMUSRWK && tree->mod->whichmodel!=YAPECMUSRWK)
	  {
	    Generic_Brent_Lk(&(tree->mod->kappa),
			     TREEBOUNDLOW,TREEBOUNDHIGH,
			     tree->mod->s_opt->min_diff_lk_global,
			     tree->mod->s_opt->brent_it_max,
			     tree->mod->s_opt->quickdirty,
			     Optwrap_Lk,NULL,tree,NULL);
	  }
	  else
	  {
	    switch(tree->io->kappaECM)
	    {
	      case kap1: 
		break;
	      case kap2:
	      case kap3: 
	      case kap6:  
	      {
		Generic_Brent_Lk(&(tree->mod->pkappa[0]),
				 TREEBOUNDLOW,TREEBOUNDHIGH,
				 tree->mod->s_opt->min_diff_lk_global,
				 tree->mod->s_opt->brent_it_max,
				 tree->mod->s_opt->quickdirty,
				 Optwrap_Lk,NULL,tree,NULL);
				 break;
	      }
	      case kap4:
	      {
		Generic_Brent_Lk(&(tree->mod->pkappa[0]),
				 TREEBOUNDLOW,TREEBOUNDHIGH,
				 tree->mod->s_opt->min_diff_lk_global,
				 tree->mod->s_opt->brent_it_max,
				 tree->mod->s_opt->quickdirty,
				 Optwrap_Lk,NULL,tree,NULL);
				 Generic_Brent_Lk(&(tree->mod->pkappa[1]),
						  TREEBOUNDLOW,TREEBOUNDHIGH,
						  tree->mod->s_opt->min_diff_lk_global,
						  tree->mod->s_opt->brent_it_max,
						  tree->mod->s_opt->quickdirty,
						  Optwrap_Lk,NULL,tree,NULL);
						  break;
	      }
	      case kap5:
	      {
		For(i,tree->mod->nkappa-1)
		{
		  Generic_Brent_Lk(&(tree->mod->unspkappa[i]),
				   -99.0,99.0,
				   tree->mod->s_opt->min_diff_lk_global,
				   tree->mod->s_opt->brent_it_max,
				   tree->mod->s_opt->quickdirty,
				   Optwrap_Lk,NULL,tree,NULL);
		}
		break;
	      }
	      default: Warn_And_Exit("Error in Kappa assignment.");
	      break;
	    }
	  }
	}
	
	if(tree->mod->s_opt->opt_omega)
	{
	  if(tree->mod->omegaSiteVar==DM0)
	  {
		  int omegai; //added by Ken 18/8/2016
		  for(omegai=0;omegai<tree->mod->nomega_part;omegai++){ //Ken 18/8

		  Generic_Brent_Lk(&(tree->mod->omega_part[omegai]),
			     TREEBOUNDLOW,TREEBOUNDHIGH,
			     tree->mod->s_opt->min_diff_lk_global,
			     tree->mod->s_opt->brent_it_max,
			     tree->mod->s_opt->quickdirty,
			     Optwrap_Lk,NULL,tree,NULL);
		  }
	  }
	  else if(tree->mod->omegaSiteVar==DMODELK)
	  {
	    For(i,tree->mod->n_w_catg)
	    {
	      Generic_Brent_Lk(&(tree->mod->omegas[i]),
			       TREEBOUNDLOW,TREEBOUNDHIGH,
			       tree->mod->s_opt->min_diff_lk_global,
			       tree->mod->s_opt->brent_it_max,
			       tree->mod->s_opt->quickdirty,
			       Optwrap_Lk,NULL,tree,NULL);
			       
	      if(i<tree->mod->n_w_catg-1)
	      {
		
		Generic_Brent_Lk(&(tree->mod->prob_omegas_uns[i]),
				-99.0,99.0,
				tree->mod->s_opt->min_diff_lk_global,
				tree->mod->s_opt->brent_it_max,
				tree->mod->s_opt->quickdirty,
				Optwrap_Lk,NULL,tree,NULL);
	      }
	    }
	  }
	  else if(tree->mod->omegaSiteVar==DGAMMAK)
	  {
	    Generic_Brent_Lk(&(tree->mod->alpha),
			     2e-3,TREEBOUNDHIGH,
			     tree->mod->s_opt->min_diff_lk_global,
			     tree->mod->s_opt->brent_it_max,
			     tree->mod->s_opt->quickdirty,
			     Optwrap_Lk,NULL,tree,NULL);
			     Generic_Brent_Lk(&(tree->mod->beta),
					      2e-3,999.0,
					      tree->mod->s_opt->min_diff_lk_global,
					      tree->mod->s_opt->brent_it_max,
					      tree->mod->s_opt->quickdirty,
					      Optwrap_Lk,NULL,tree,NULL);
	  }
	}
	printf("opt state freq: %d\n",tree->mod->s_opt->opt_state_freq);
	if(tree->mod->s_opt->opt_state_freq)
	{
	  switch(tree->mod->freq_model)
	  {
	    case F1XSENSECODONS:
	    {
	      For(i,tree->mod->num_base_freq-1)
	      {
		Generic_Brent_Lk(&(tree->mod->pi_unscaled[i]),
				 -99.0,99.0,
				 tree->mod->s_opt->min_diff_lk_global,
				 tree->mod->s_opt->brent_it_max,
				 tree->mod->s_opt->quickdirty,
				 Optwrap_Lk,NULL,tree,NULL);
	      }
	      break;
	    }
	    case F1X4:
	    {
	      For(i,tree->mod->num_base_freq-1)
	      {
		Generic_Brent_Lk(&(tree->mod->uns_base_freq[i]),
				 -99.0,99.0,
				 tree->mod->s_opt->min_diff_lk_global,
				 tree->mod->s_opt->brent_it_max,
				 tree->mod->s_opt->quickdirty,
				 Optwrap_Lk,NULL,tree,NULL);
	      }
	     PhyML_Printf("%lf %lf %lf %lf\n",tree->mod->uns_base_freq[0],tree->mod->uns_base_freq[1],tree->mod->uns_base_freq[2],tree->mod->uns_base_freq[3]);
	      break;
	    }
	    case F3X4:
	    case CF3X4:
	    {
	      for(i=0;i<3;i++)
	      {
	    	  printf("doing optimizations\n");

		Generic_Brent_Lk(&(tree->mod->uns_base_freq[i]),
				 -99.0,99.0,
				 tree->mod->s_opt->min_diff_lk_global,
				 tree->mod->s_opt->brent_it_max,
				 tree->mod->s_opt->quickdirty,
				 Optwrap_Lk,NULL,tree,NULL);
	      }
	      for(i=4;i<7;i++)
	      {
		Generic_Brent_Lk(&(tree->mod->uns_base_freq[i]),
				 -99.0,99.0,
				 tree->mod->s_opt->min_diff_lk_global,
				 tree->mod->s_opt->brent_it_max,
				 tree->mod->s_opt->quickdirty,
				 Optwrap_Lk,NULL,tree,NULL);
	      }
	      for(i=8;i<11;i++)
	      {
		Generic_Brent_Lk(&(tree->mod->uns_base_freq[i]),
				 -99.0,99.0,
				 tree->mod->s_opt->min_diff_lk_global,
				 tree->mod->s_opt->brent_it_max,
				 tree->mod->s_opt->quickdirty,
				 Optwrap_Lk,NULL,tree,NULL);
	      }
	      break;
	    }
	    default:
	      break;
	  }
	}
	
	if(tree->mod->s_opt->opt_pinvar)
	{
	  Generic_Brent_Lk(&(tree->mod->pinvar),
			   0.0,1.0,
			   tree->mod->s_opt->min_diff_lk_global,
			   tree->mod->s_opt->brent_it_max,
			   tree->mod->s_opt->quickdirty,
			   Optwrap_Lk,NULL,tree,NULL);
	}
	
	if(tree->mod->n_catg>1 && tree->mod->n_w_catg==1 && tree->mod->s_opt->opt_alphaCD)
	{ 
	  Generic_Brent_Lk(&(tree->mod->alpha),
			   2e-3,TREEBOUNDHIGH,
			   tree->mod->s_opt->min_diff_lk_global,
			   tree->mod->s_opt->brent_it_max,
			   tree->mod->s_opt->quickdirty,
			   Optwrap_Lk,NULL,tree,NULL);
	}
	
	if(tree->mod->pcaModel==1)
	{
            printf("hier\n");
	  For(i,tree->mod->npcs)
	  {
	    Generic_Brent_Lk(&(tree->mod->pcsC[i]),
			     ((-1) * TREEBOUNDHIGH),TREEBOUNDHIGH,
			     tree->mod->s_opt->min_diff_lk_global,
			     tree->mod->s_opt->brent_it_max,
			     tree->mod->s_opt->quickdirty,
			     Optwrap_Lk,NULL,tree,NULL);
	  }	
	}
	
	tree->both_sides = init_both_sides;
	
	if(tree->io->heuristicExpm)
	{
	  tree->mod->update_eigen = 1; 
	  tree->io->optParam=0;
	  tree->io->expm=EIGEN;
	  Lk(tree);
	  tree->mod->update_eigen = 0;  
	}
	else
	{
	  tree->mod->update_eigen = 0;
	  if(tree->both_sides) Lk(tree); /* Needed to update all partial likelihoods.*/  
	}
      }
    }
    if(verbose) 
    {
      if(tree->mod->s_opt->opt_kappa)
      {
	if(tree->mod->whichmodel!=GYECMK07WK  && tree->mod->whichmodel!=GYECMK07WKF && 
	  tree->mod->whichmodel!=GYECMS05WK  && tree->mod->whichmodel!=GYECMS05WKF &&
	  tree->mod->whichmodel!=MGECMK07WK  && tree->mod->whichmodel!=MGECMK07WKF && 
	  tree->mod->whichmodel!=MGECMS05WK  && tree->mod->whichmodel!=MGECMS05WKF &&
	  tree->mod->whichmodel!=YAPECMK07WK && tree->mod->whichmodel!=YAPECMK07WKF && 
	  tree->mod->whichmodel!=YAPECMS05WK && tree->mod->whichmodel!=YAPECMS05WKF &&
	  tree->mod->whichmodel!=GYECMUSRWK && tree->mod->whichmodel!=GYECMUSRWKF &&
	  tree->mod->whichmodel!=MGECMUSRWKF && tree->mod->whichmodel!=YAPECMUSRWKF &&
	  tree->mod->whichmodel!=MGECMUSRWK && tree->mod->whichmodel!=YAPECMUSRWK ) 
	{
	  Print_Lk(tree,"[ts/tv ratio        ]");
	  PhyML_Printf("[%.2f ]",tree->mod->kappa);
	}
	else
	{
	  switch(tree->io->kappaECM)
	  {
	    case kap1: 
	      break;
	    case kap2:
	    case kap3: 
	    {
	      Print_Lk(tree,"[Emp. ts/tv ratio   ]");
	      PhyML_Printf("[%.2f ]",tree->mod->pkappa[0]);
	      break;
	    }
	    case kap4:
	    {
	      Print_Lk(tree,"[Emp. ts/tv ratio   ]");
	      PhyML_Printf("[ %.2f %.2f ]",tree->mod->pkappa[0], tree->mod->pkappa[1]);
	      break;
	    }
	    case kap5:
	    {
	      Print_Lk(tree,"[Emp. ts/tv ratio   ]");
	      break;
	    }
	    case kap6: 
	    {
	      Print_Lk(tree,"[Multi-NT parameter ]");
	      PhyML_Printf("[%.2f ]",tree->mod->pkappa[0]);
	      break;
	    }
	    default: Warn_And_Exit("Error in Kappa assignment.");
	    break;
	  }
	}
      }
      if(tree->mod->s_opt->opt_omega)
      {
	if(tree->mod->omegaSiteVar==DM0)
	{
	  if(tree->mod->whichmodel!=GYECMK07WK  && tree->mod->whichmodel!=GYECMK07WKF &&
	    tree->mod->whichmodel!=GYECMS05WK  && tree->mod->whichmodel!=GYECMS05WKF &&
	    tree->mod->whichmodel!=MGECMK07WK  && tree->mod->whichmodel!=MGECMK07WKF &&
	    tree->mod->whichmodel!=MGECMS05WK  && tree->mod->whichmodel!=MGECMS05WKF &&
	    tree->mod->whichmodel!=YAPECMK07WK && tree->mod->whichmodel!=YAPECMK07WKF&&
	    tree->mod->whichmodel!=YAPECMS05WK && tree->mod->whichmodel!=YAPECMS05WKF &&
	    tree->mod->whichmodel!=GYECMUSRWK && tree->mod->whichmodel!=GYECMUSRWKF &&
	    tree->mod->whichmodel!=MGECMUSRWKF && tree->mod->whichmodel!=YAPECMUSRWKF &&
	    tree->mod->whichmodel!=MGECMUSRWK &&tree->mod->whichmodel!=YAPECMUSRWK      )
	  { 
		  int omegai; //added by Ken 18/8/2016
		  for(omegai=0;omegai<tree->mod->nomega_part;omegai++){
			  Print_Lk(tree,"[dn/ds ratio        ]");
			  PhyML_Printf("[%.2f ]",tree->mod->omega_part[omegai]);
		  }
	  }
	  else 
	  {                 
	    Print_Lk(tree,"[Emp. dn/ds ratio   ]");
    	if(tree->mod->nparts > 1){printf("options not compatible with partitioned model error 5\n");exit(EXIT_FAILURE);}
	    PhyML_Printf("[%.2f ]",Omega_ECMtoMmechModels(tree->mod->pi, tree->mod->qmat_part[0], tree->mod->qmat_buff_part[0], tree->mod->ns, tree->mod->n_w_catg));
	  }
	}
	else if(tree->mod->omegaSiteVar==DMODELK)
	{	
	  if(tree->mod->n_w_catg<5)
	  {
	    if((tree->mod->whichmodel!=MGECMUSRWK)&&(tree->mod->whichmodel!=GYECMK07WK)&&(tree->mod->whichmodel!=GYECMK07WKF)&&(tree->mod->whichmodel!=GYECMS05WK)&&(tree->mod->whichmodel!=GYECMS05WKF) && (tree->mod->whichmodel!=MGECMK07WK)&&(tree->mod->whichmodel!=MGECMK07WKF)&&(tree->mod->whichmodel!=MGECMS05WK)&&(tree->mod->whichmodel!=MGECMS05WKF) && (tree->mod->whichmodel!=YAPECMK07WK)&&(tree->mod->whichmodel!=YAPECMK07WKF)&&(tree->mod->whichmodel!=YAPECMS05WK)&&(tree->mod->whichmodel!=YAPECMS05WKF) &&(tree->mod->whichmodel!=GYECMUSRWK && tree->mod->whichmodel!=GYECMUSRWKF && tree->mod->whichmodel!=MGECMUSRWKF && tree->mod->whichmodel!=YAPECMUSRWKF && tree->mod->whichmodel!=YAPECMUSRWK))
	    { 
	      int m;         
	      Print_Lk(tree,"[dn/ds ratio        ]");
	      printf("[");
	      For(m,tree->mod->n_w_catg) printf("%.2f ",tree->mod->omegas[m]);
	      printf("]");
	      Print_Lk(tree,"[dn/ds frequencies  ]");
	      printf("[");
	      For(m,tree->mod->n_w_catg) printf("%.2f ",tree->mod->prob_omegas[m]);
	      printf("]");
          fflush(NULL);
	    }
	    else
	    {
	      int m;         
	      Print_Lk(tree,"[Emp. dn/ds ratio   ]");
	      printf("[");
	      if(tree->mod->nparts > 1){printf("options not compatible with partitioned model error 6\n");exit(EXIT_FAILURE);}
	      For(m,tree->mod->n_w_catg) printf("%.2f ",Omega_ECMtoMmechModels(tree->mod->pi, tree->mod->qmat_part[0]+m*tree->mod->ns*tree->mod->ns, tree->mod->qmat_buff_part[0]+m*tree->mod->ns*tree->mod->ns, tree->mod->ns, tree->mod->n_w_catg));
	      printf("]");
	      Print_Lk(tree,"[dn/ds frequencies  ]");
	      printf("[");
	      For(m,tree->mod->n_w_catg) printf("%.2f ",tree->mod->prob_omegas[m]);
	      printf("]");
	    }
	  }
	  else
	  {
	    Print_Lk(tree,"[dn/ds ratio        ]");
	    Print_Lk(tree,"[dn/ds frequencies  ]");
	  }
	}
	else if(tree->mod->omegaSiteVar==DGAMMAK)
	{
	  if(tree->mod->n_w_catg<5)
	  {
	    int m; 
	    if((tree->mod->whichmodel!=MGECMK07WK)&&(tree->mod->whichmodel!=GYECMK07WK)&&(tree->mod->whichmodel!=GYECMK07WKF)&&(tree->mod->whichmodel!=GYECMS05WK)&&(tree->mod->whichmodel!=GYECMS05WKF)&&(tree->mod->whichmodel!=MGECMK07WK)&&(tree->mod->whichmodel!=MGECMK07WKF)&&(tree->mod->whichmodel!=MGECMS05WK)&&(tree->mod->whichmodel!=MGECMS05WKF)&&(tree->mod->whichmodel!=YAPECMK07WK)&&(tree->mod->whichmodel!=YAPECMK07WKF)&&(tree->mod->whichmodel!=YAPECMS05WK)&&(tree->mod->whichmodel!=YAPECMS05WKF)&&(tree->mod->whichmodel!=GYECMUSRWK && tree->mod->whichmodel!=GYECMUSRWKF && tree->mod->whichmodel!=MGECMUSRWKF && tree->mod->whichmodel!=YAPECMUSRWKF && tree->mod->whichmodel!=YAPECMUSRWK))
	    { 
	      Print_Lk(tree,"[dn/ds ratio        ]");
	      printf("[");
	      For(m,tree->mod->n_w_catg) printf("%.2f ",tree->mod->omegas[m]);
	      printf("]");
	      Print_Lk(tree,"[Alpha and Beta     ]");
	      PhyML_Printf("[%.2f %.2f ]",tree->mod->alpha, tree->mod->beta);
	    }
	    else
	    {
	      Print_Lk(tree,"[Emp. dn/ds ratio   ]");
	      printf("[");
      	if(tree->mod->nparts > 1){printf("options not compatible with partitioned model error 7\n");exit(EXIT_FAILURE);}
	      For(m,tree->mod->n_w_catg) printf("%.2f ",Omega_ECMtoMmechModels(tree->mod->pi, tree->mod->qmat_part[0]+(m*tree->mod->ns*tree->mod->ns), tree->mod->qmat_buff_part[0]+(m*tree->mod->ns*tree->mod->ns), tree->mod->ns,tree->mod->n_w_catg ));
	      printf("]");
	      Print_Lk(tree,"[Alpha and Beta     ]");
	      PhyML_Printf("[%.2f %.2f ]",tree->mod->alpha, tree->mod->beta);
	    }
	  }
	  else
	  {
	    Print_Lk(tree,"[dn/ds ratio        ]");
	    Print_Lk(tree,"[Alpha and Beta     ]");
	    PhyML_Printf("[%.2f %.2f ]",tree->mod->alpha, tree->mod->beta);
	  } 
	}
      }
      
      if(tree->mod->opthotness==1){
    	 int c,d;
    	 for(c=0;c<tree->mod->nmotifs;c++){
    	  char *info = malloc(22);
    	  char motifh[10];
    	  sprintf(motifh,"%d",tree->mod->motif_hotness[c]);
    	  strcpy(info, "[");
    	  strcat(info, tree->mod->motifs[c]);
    	  strcat(info, " h");
    	  strcat(info, motifh);
    	  for(d=strlen(info);d<20;d++){
    		  strcat(info, " ");
    	  }
    	  strcat(info,"]");

    	  Print_Lk(tree,info);
    	  PhyML_Printf("[%.3f]",tree->mod->hotness[tree->mod->motif_hotness[c]]);
    	 }

      }

      if(tree->mod->s_opt->opt_state_freq)
      {
	switch(tree->mod->freq_model)
	{
	  case F1XSENSECODONS:
	  {
	    strcpy(s,"[F1x"); 
	    sprintf(r,"%d",tree->mod->ns);
	    strcat(s,r);
	    strcat(s,"        freqs.]");
	    i=-1;
	    while(s[++i]!=']');
	    s[++i]=0;
	    Print_Lk(tree,s);	    
	    break;
	  }
	  case F1X4: Print_Lk(tree,"[F1x4 freqs.        ]"); break;
	  case F3X4: Print_Lk(tree,"[F3x4 freqs.        ]"); break;
	  case CF3X4: Print_Lk(tree,"[CF3x4 freqs.       ]"); break;
	  default:
	    break;
	}
      }
      
      if(tree->mod->s_opt->opt_pinvar)
      {
	Print_Lk(tree,"[P-inv              ]");
	PhyML_Printf("[%.2f]",tree->mod->pinvar);
      }
      
      if(tree->mod->n_catg>1 && tree->mod->n_w_catg==1 && tree->mod->s_opt->opt_alphaCD)
      {
	Print_Lk(tree,"[Alpha              ]");
	PhyML_Printf("[%.2f]",tree->mod->alpha);
      }
      
      if(tree->mod->pcaModel==1)
      {
	Print_Lk(tree,"[PCA linear coeff.  ]");
      }
    }          
  } ///////!<End of CODON MODELS///////////

}
  


static phydbl sqrarg;

void BFGS(t_tree *tree, 
	  phydbl *p, 
	  int n, 
	  phydbl gtol, 
	  phydbl step_size,
	  phydbl(*func)(t_tree *tree), 
	  int(*dfunc)(t_tree *tree,phydbl *param,int n_param,phydbl stepsize,phydbl(*func)(t_tree *tree),phydbl *derivatives), 
	  int(*lnsrch)(t_tree *tree, int n, phydbl *xold, phydbl fold,phydbl *g, phydbl *p, phydbl *x,phydbl *f, phydbl stpmax, int *check),
	  int *failed)
{

  int check,i,its,j;
  phydbl den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test,fret;
  phydbl *dg,*g,*hdg,**hessin,*pnew,*xi;
  
  hessin = (phydbl **)mCalloc(n,sizeof(phydbl *));
  For(i,n) hessin[i] = (phydbl *)mCalloc(n,sizeof(phydbl));
  dg   = (phydbl *)mCalloc(n,sizeof(phydbl ));
  g    = (phydbl *)mCalloc(n,sizeof(phydbl ));
  pnew = (phydbl *)mCalloc(n,sizeof(phydbl ));
  hdg  = (phydbl *)mCalloc(n,sizeof(phydbl ));
  xi   = (phydbl *)mCalloc(n,sizeof(phydbl ));
  

/*   PhyML_Printf("\n. ENTER BFGS WITH: %f\n",Lk(tree)); */

  fp=(*func)(tree);
  (*dfunc)(tree,p,n,step_size,func,g);

  for (i=0;i<n;i++) 
    {
      for (j=0;j<n;j++) hessin[i][j]=0.0;
      hessin[i][i]=1.0;
      xi[i] = -g[i];
      sum += p[i]*p[i];
    }

  stpmax=STPMX*MAX(SQRT(sum),(phydbl)n);

  for(its=1;its<=ITMAX;its++) 
    {
      lnsrch(tree,n,p,fp,g,xi,pnew,&fret,stpmax,&check);

/*       PhyML_Printf("BFGS -> %f\n",tree->c_lnL); */

      fp = fret;
      
      for (i=0;i<n;i++) 
	{
	  xi[i]=pnew[i]-p[i];
	  p[i]=pnew[i];
	}

      test=0.0;
      for (i=0;i<n;i++) 
	{
	  temp=FABS(xi[i])/MAX(FABS(p[i]),1.0);
	  if (temp > test) test=temp;
	}
      if (test < TOLX) 
	{
	  (*func)(tree);
	  For(i,n) free(hessin[i]);
	  free(hessin);
	  free(xi);
	  free(pnew);
	  free(hdg);
	  free(g);
	  free(dg);   

	  if(its == 1) 
	    {
/* 	      PhyML_Printf("\n. WARNING : BFGS failed ! \n"); */
	      *failed = 1;
	    }
	  return;
	}

      for (i=0;i<n;i++) dg[i]=g[i];

      (*dfunc)(tree,p,n,step_size,func,g);

      test=0.0;
      den=MAX(fret,1.0);
      for (i=0;i<n;i++) 
	{
	  temp=FABS(g[i])*MAX(FABS(p[i]),1.0)/den;
	  if (temp > test) test=temp;
	}
      if (test < gtol) 
	{
	  (*func)(tree);
	  For(i,n) free(hessin[i]);
	  free(hessin);
	  free(xi);
	  free(pnew);
	  free(hdg);
	  free(g);
	  free(dg);   
	  return;
	}

    for (i=0;i<n;i++) dg[i]=g[i]-dg[i];

    for (i=0;i<n;i++) 
      {
	hdg[i]=0.0;
	for (j=0;j<n;j++) hdg[i] += hessin[i][j]*dg[j];
      }

    fac=fae=sumdg=sumxi=0.0;
    for (i=0;i<n;i++) 
      {
	fac += dg[i]*xi[i];
	fae += dg[i]*hdg[i];
	sumdg += SQR(dg[i]);
	sumxi += SQR(xi[i]);
      }
    
    if(fac*fac > EPS*sumdg*sumxi) 
      {
	fac=1.0/fac;
	fad=1.0/fae;
	for (i=0;i<n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
	for (i=0;i<n;i++) 
	  {
	    for (j=0;j<n;j++) 
	      {
		hessin[i][j] += fac*xi[i]*xi[j]
		  -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
	      }
	  }
      }
    for (i=0;i<n;i++) 
      {
	xi[i]=0.0;
	for (j=0;j<n;j++) xi[i] -= hessin[i][j]*g[j];
      }
    }
/*   PhyML_Printf("\n. Too many iterations in BFGS...\n"); */
  *failed = 1;
  For(i,n) free(hessin[i]);
  free(hessin);
  free(xi);
  free(pnew);
  free(hdg);
  free(g);
  free(dg);   
}

#undef ITMAX
#undef EPS
#undef TOLX
#undef STPMX

/*********************************************************/


#define ALF 1.0e-4
#define TOLX 1.0e-7
#undef ALF
#undef TOLX
#undef NRANSI

/*********************************************************/

#define ALF 1.0e-4
#define TOLX 1.0e-7

/* void Optimize_Global_Rate(t_tree *tree) */
/* { */
/*     PhyML_Printf("\n. Global rate (%f->)",tree->c_lnL); */
/*     Optimize_Single_Param_Generic(tree,&(tree->tbl),tree->tbl,BL_MIN,1.E+4,100); */
/*     PhyML_Printf("%f)\n",tree->c_lnL); */
/* } */


#undef ALF
#undef TOLX
#undef NRANSI

/*********************************************************/

phydbl Dist_F_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, int n_iter_max, 
		    phydbl *param, phydbl *F, model *mod)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL,init_lnL, curr_lnL;

  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x = w = v = bx;
  old_lnL = UNLIKELY;
  fw = fv = fx = -Lk_Dist(F,FABS(bx),mod);
  curr_lnL = init_lnL = -fw;

  for(iter=1;iter<=BRENT_ITMAX;iter++) 
    {
      xm=0.5*(a+b);

      tol2=2.0*(tol1=tol*FABS(x)+BRENT_ZEPS);

      if(
	 ((FABS(curr_lnL-old_lnL) < mod->s_opt->min_diff_lk_local) && 
	  (curr_lnL > init_lnL - mod->s_opt->min_diff_lk_local)) ||	 
	  (iter > n_iter_max - 1)
	 )	 
	{
	  *param = x;
	  curr_lnL = Lk_Dist(F,*param,mod);
	  return -curr_lnL;
	}
      
      if(FABS(e) > tol1) 
	{
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if(q > 0.0) p = -p;
	  q=FABS(q);
	  etemp=e;
	  e=d;
	  if(FABS(p) >= FABS(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    {
	      d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	      /*                   PhyML_Printf("Golden section step\n"); */
	    }
	  else
	    {
	      d=p/q;
	      u=x+d;
	      if (u-a < tol2 || b-u < tol2)
		d=SIGN(tol1,xm-x);
	      /*                   PhyML_Printf("Parabolic step\n"); */
	    }
        }
      else
	{
	  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	  /*               PhyML_Printf("Golden section step (default)\n"); */
	}
      
      u=(FABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      if(u<BL_MIN) u = BL_MIN;
      (*param) = FABS(u);
      old_lnL = curr_lnL;
      fu = -Lk_Dist(F,FABS(u),mod);
      curr_lnL = -fu;      
/*       PhyML_Printf("param=%f LOGlk=%f\n",*param,fu); */
      
/*       if(fu <= fx)  */
      if(fu < fx) 
	{
	  if(iter > n_iter_max) return -fu;

	  if(u >= x) a=x; else b=x;
	  SHFT(v,w,x,u)
	  SHFT(fv,fw,fx,fu)
	} 
      else
	{
	  if (u < x) a=u; else b=u;
/* 	  if (fu <= fw || w == x)  */
	  if (fu < fw || FABS(w-x) < SMALL)
	    {
	      v=w;
	      w=u;
	      fv=fw;
	      fw=fu;
	    } 
/* 	  else if (fu <= fv || v == x || v == w)  */
	  else if (fu < fv || FABS(v-x) < SMALL || FABS(v-w) < SMALL)
	    {
	      v=u;
	      fv=fu;
	    }
	}
    }
  Exit("\n. Too many iterations in BRENT !");
  return(-1);
}

/*********************************************************/

void Opt_Dist_F(phydbl *dist, phydbl *F, model *mod)
{
  phydbl ax,bx,cx;

  if(*dist < BL_MIN) *dist = BL_MIN;

  ax = BL_MIN;
  bx =  (*dist);
  cx = BL_MAX;

/*   Dist_F_Brak(&ax,&bx,&cx,F,dist,mod); */
  Dist_F_Brent(ax,bx,cx,1.E-10,1000,dist,F,mod);
}

/*********************************************************/

phydbl Optwrap_Lk(t_edge *b, t_tree *tree, supert_tree *stree)
{
  Lk(tree);
  return tree->c_lnL;
}

/*********************************************************/

phydbl Optwrap_Lk_At_Given_Edge(t_edge *b, t_tree *tree, supert_tree *stree)
{
  Lk_At_Given_Edge(b,tree);
  return tree->c_lnL;
}

/*********************************************************/

phydbl Generic_Brent_Lk(phydbl *param, phydbl ax, phydbl cx, phydbl tol, 
			int n_iter_max, int quickdirty,
			phydbl (*obj_func)(t_edge *,t_tree *,supert_tree *), 
			t_edge *branch, t_tree *tree, supert_tree *stree)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL,init_lnL,cur_lnL;
  phydbl bx = *param;

  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  old_lnL = UNLIKELY;
  (*param) = FABS(bx);
  fw=fv=fx=fu=-(*obj_func)(branch,tree,stree);
  init_lnL = -fw;

/*   PhyML_Printf("\n. init_lnL = %f a=%f b=%f c=%f\n",init_lnL,ax,bx,cx); */

  for(iter=1;iter<=BRENT_ITMAX;iter++) 
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*FABS(x)+BRENT_ZEPS);

      cur_lnL = (stree)?(stree->tree->c_lnL):(tree->c_lnL);

      if((cur_lnL > init_lnL + tol) && (quickdirty))
	{
	  (*param) = x;
	  (*obj_func)(branch,tree,stree);
/* 	  Exit("\n"); */
	  return (stree)?(stree->tree->c_lnL):(tree->c_lnL);	  
	}

/*       if(((FABS(cur_lnL-old_lnL) < tol) && (cur_lnL > init_lnL - tol)) || (iter > n_iter_max - 1)) */
      if((FABS(cur_lnL-old_lnL) < tol) || (iter > n_iter_max - 1))
	{
	  (*param) = x;
	  (*obj_func)(branch,tree,stree);
/* 	  Exit("\n"); */
	  return (stree)?(stree->tree->c_lnL):(tree->c_lnL);	  
	}
      
      if(FABS(e) > tol1) 
	{
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if(q > 0.0) p = -p;
	  q=FABS(q);
	  etemp=e;
	  e=d;
	  if(FABS(p) >= FABS(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    {
	      d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
/* 	      PhyML_Printf("Golden section step\n"); */
	    }
	  else
	    {
	      d=p/q;
	      u=x+d;
	      if (u-a < tol2 || b-u < tol2) d=SIGN(tol1,xm-x);
/* 	      PhyML_Printf("Parabolic step [e=%f]\n",e); */
	    }
        }
      else
	{
	  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
/* 	  PhyML_Printf("Golden section step (default) [e=%f tol1=%f a=%f b=%f d=%f]\n",e,tol1,a,b,d); */
	}
      
      u=(FABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      (*param) = FABS(u);
      old_lnL = (stree)?(stree->tree->c_lnL):(tree->c_lnL);
      fu = -(*obj_func)(branch,tree,stree);
      
/*       PhyML_Printf("\n. iter=%d/%d param=%f LOGlk=%f",iter,BRENT_ITMAX,*param,tree->c_lnL); */

      if(fu <= fx)
	{
	  if(u >= x) a=x; else b=x;
	  SHFT(v,w,x,u)
	  SHFT(fv,fw,fx,fu)
	} 
      else
	{
	  if (u < x) a=u; else b=u;
/* 	  if (fu <= fw || w == x) */
	  if (fu < fw || FABS(w-x) < SMALL)
	    {
	      v=w;
	      w=u;
	      fv=fw;
	      fw=fu;
	    } 
/* 	  else if (fu <= fv || v == x || v == w) */
	  else if (fu < fv || FABS(v-x) < SMALL || FABS(v-w) < SMALL)
	    {
	      v=u;
	      fv=fu;
	    }
	}
    }

  Exit("\n. Too many iterations in BRENT !");
  return(-1);
  /* Not Reached ??  *param=x;   */
  /* Not Reached ??  return fx; */
}

/*********************************************************/
phydbl Br_Len_Brent_Codon_Pairwise(phydbl ax, phydbl bx, phydbl cx, phydbl tol, phydbl *b_fcus, phydbl *Pij, phydbl *pi, eigen *eigenStruct, calign *data, int ns, int n_iter_max, int quickdirty, phydbl *uexpt, phydbl *expt) //!< Added by Marcelo.
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL, init_lnL, curr_lnL;

  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  old_lnL = UNLIKELY;
  (*b_fcus) = FABS(bx);
  fw=fv=fx=fu=-LK_Codon_Pairwise(data, Pij, pi, ns, (*b_fcus), eigenStruct, uexpt, expt);
  curr_lnL = init_lnL = -fw;

  for(iter=1;iter<=BRENT_ITMAX;iter++)
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*FABS(x)+BRENT_ZEPS);

      if((curr_lnL > init_lnL + tol) && (quickdirty))
	{
	  (*b_fcus) = x;
	  curr_lnL=LK_Codon_Pairwise(data, Pij, pi, ns, (*b_fcus), eigenStruct, uexpt, expt);
	  return curr_lnL;	  
	}

      if((FABS(curr_lnL-old_lnL) < tol) || (iter > n_iter_max - 1))
	{
	  (*b_fcus)=x;
	  curr_lnL=LK_Codon_Pairwise(data, Pij, pi, ns, (*b_fcus), eigenStruct, uexpt, expt);
	  return curr_lnL;	
	}
      
      if(FABS(e) > tol1)
	{
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if(q > 0.0) p = -p;
	  q=FABS(q);
	  etemp=e;
	  e=d;
	  if(FABS(p) >= FABS(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	  else{
	    d=p/q;
	    u=x+d;
	    if (u-a < tol2 || b-u < tol2)
	      d=SIGN(tol1,xm-x);
	  }
	}
      else
	{
	  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	}
      u=(FABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));

      (*b_fcus)=FABS(u);
      old_lnL = curr_lnL;
      fu=-LK_Codon_Pairwise(data, Pij, pi, ns, (*b_fcus), eigenStruct, uexpt, expt);
      curr_lnL=-fu;
      if(fu < fx)
	{
	  if(u > x) a=x; else b=x;
	  SHFT(v,w,x,u)
	  SHFT(fv,fw,fx,fu)
	}
      else
	{
	  if (u < x) a=u; else b=u;
	  if (fu < fw || FABS(w-x) < SMALL)
	    {
	      v=w;
	      w=u;
	      fv=fw;
	      fw=fu;
	    }
	  else if (fu < fv || FABS(v-x) < SMALL || FABS(v-w) < SMALL)
	    {
	      v=u;
	      fv=fu;
	    }
	}
    }
  if(iter > BRENT_ITMAX) PhyML_Printf("\n. Too many iterations in BRENT (%d) (%f)",iter,(*b_fcus));
  return(-1);
}

/*********************************************************/
#define BFGS
/*
#define SR1
#define DFP
*/
phydbl SIZEp=0;
int noisy=0, Iround=0, NFunCall=0; 
int BFGS_from_CODEML (phydbl *f, t_tree *tree, phydbl *x, phydbl xb[120][2], phydbl space[], phydbl e, int n)
{
/* n-variate minimization with bounds using the BFGS algorithm
     g0[n] g[n] p[n] x0[n] y[n] s[n] z[n] H[n*n] C[n*n] tv[2*n]
     xmark[n],ix[n]
   Size of space should be (check carefully?)
      #define spaceming2(n) ((n)*((n)*2+9+2)*sizeof(phydbl))
   nfree: # free variables
   xmark[i]=0 for inside space; -1 for lower boundary; 1 for upper boundary.
   x[] has initial values at input and returns the estimates in return.
   ix[i] specifies the i-th free parameter
  
   ALL CREDIT GOES TO PROF YANG and CODEML
*/
   int i,j, i1,i2,it, maxround=10000, fail=0, *xmark, *ix, nfree;
   int Ngoodtimes=2, goodtimes=0;
   phydbl small=1.e-30, sizep0=0;     /* small value for checking |w|=0 */
   phydbl f0, *g0, *g, *p, *x0, *y, *s, *z, *H, *C, *tv;
   phydbl w,v, alpha, am, h, maxstep=8;

   if(n==0) return(0);
   g0=space;   g=g0+n;  p=g+n;   x0=p+n;
   y=x0+n;     s=y+n;   z=s+n;   H=z+n;  C=H+n*n, tv=C+n*n;
   xmark=(int*)(tv+2*n);  ix=xmark+n;

   for(i=0; i<n; i++)  { xmark[i]=0; ix[i]=i; }
   
   for(i=0,nfree=0;i<n;i++) {
      if(x[i]<=xb[i][0]) { x[i]=xb[i][0]; xmark[i]=-1; continue; }
      if(x[i]>=xb[i][1]) { x[i]=xb[i][1]; xmark[i]= 1; continue; }
      ix[nfree++]=i;
   }
   f0=*f=LK_BFGS_from_CODEML(tree,x,n);

   xtoy(x,x0,n);
   SIZEp=99;

   gradientB (n, x0, f0, g0, tree, tv, xmark);

   identity (H,nfree);
   for(Iround=0; Iround<maxround; Iround++) {
     

      for (i=0,zero(p,n); i<nfree; i++)  For (j,nfree)
         p[ix[i]] -= H[i*nfree+j]*g0[ix[j]];
      sizep0 = SIZEp; 
      SIZEp  = norm(p,n);      /* check this */

      for (i=0,am=maxstep; i<n; i++) {  /* max step length */
         if (p[i]>0 && (xb[i][1]-x0[i])/p[i]<am) am=(xb[i][1]-x0[i])/p[i];
         else if (p[i]<0 && (xb[i][0]-x0[i])/p[i]<am) am=(xb[i][0]-x0[i])/p[i];
      }

      if (Iround==0) {
         h=fabs(2*f0*.01/innerp(g0,p,n));  /* check this?? */
         h=min2(h,am/2000);

      }
      else {
         h=norm(s,nfree)/SIZEp;
         h=max2(h,am/500);
      }
      h = max2(h,1e-5);   h = min2(h,am/5);
      *f = f0;
      alpha = LineSearch2(tree,f,x0,p,h,am, min2(1e-3,e), tv,n); /* n or nfree? */

      if (alpha<=0) {
         if (fail) {
            if (AlwaysCenter) { Iround=maxround;  break; }
            else { AlwaysCenter=1; identity(H,n); fail=1; }
         }
         else   
            { if(noisy>2) printf(".. ");  identity(H,nfree); fail=1; }
      }
      else  {
         fail=0;
         For(i,n)  x[i]=x0[i]+alpha*p[i];
         w=min2(2,e*1000); if(e<1e-4 && e>1e-6) w=0.01;

         if(Iround==0 || SIZEp<sizep0 || (SIZEp<.001 && sizep0<.001)) goodtimes++;
         else  goodtimes=0;
         if((n==1||goodtimes>=Ngoodtimes) && SIZEp<(e>1e-5?1:.001)
            && H_end(x0,x,f0,*f,e,e,n))
            break;
      }
     
      gradientB (n, x, *f, g, tree, tv, xmark);
/*
for(i=0; i<n; i++) fprintf(frst,"%9.5f", x[i]); fprintf(frst, "%6d",AlwaysCenter);
for(i=0; i<n; i++) fprintf(frst,"%9.2f", g[i]); FPN(frst);
*/
      /* modify the working set */
      for(i=0; i<n; i++) {         /* add constraints, reduce H */
         if (xmark[i]) continue;
         if (fabs(x[i]-xb[i][0])<1e-6 && -g[i]<0)  xmark[i]=-1;
         else if (fabs(x[i]-xb[i][1])<1e-6 && -g[i]>0)  xmark[i]=1;
         if (xmark[i]==0) continue;
         xtoy (H, C, nfree*nfree);
         for(it=0; it<nfree; it++) if (ix[it]==i) break;
         for (i1=it; i1<nfree-1; i1++) ix[i1]=ix[i1+1];
         for (i1=0,nfree--; i1<nfree; i1++) For (i2,nfree)
            H[i1*nfree+i2]=C[(i1+(i1>=it))*(nfree+1) + i2+(i2>=it)];
      }
      for (i=0,it=0,w=0; i<n; i++) {  /* delete a constraint, enlarge H */
         if (xmark[i]==-1 && -g[i]>w)     { it=i; w=-g[i]; }
         else if (xmark[i]==1 && -g[i]<-w) { it=i; w=g[i]; }
      }
      if (w>10*SIZEp/nfree) {          /* *** */
         xtoy (H, C, nfree*nfree);
         For (i1,nfree) For (i2,nfree) H[i1*(nfree+1)+i2]=C[i1*nfree+i2];
         For (i1,nfree+1) H[i1*(nfree+1)+nfree]=H[nfree*(nfree+1)+i1]=0;
         H[(nfree+1)*(nfree+1)-1]=1;
         xmark[it]=0;   ix[nfree++]=it;
      }

     
      for (i=0,f0=*f; i<nfree; i++)
        {  y[i]=g[ix[i]]-g0[ix[i]];  s[i]=x[ix[i]]-x0[ix[i]]; }
      For (i,n) { g0[i]=g[i]; x0[i]=x[i]; }


      /* renewal of H varies with different algorithms   */
#if (defined SR1)
      /*   Symmetrical Rank One (Broyden, C. G., 1967) */
      for (i=0,w=.0; i<nfree; i++) {
         for (j=0,v=.0; j<nfree; j++) v += H[i*nfree+j] * y[j];
         z[i]=s[i] - v;
         w += y[i]*z[i];
      }
      if (fabs(w)<small)   { identity(H,nfree); fail=1; continue; }
      For (i,nfree)  For (j,nfree)  H[i*nfree+j] += z[i]*z[j]/w;
#elif (defined DFP)
      /* Davidon (1959), Fletcher and Powell (1963). */
      for (i=0,w=v=0.; i<nfree; i++) {
         for (j=0,z[i]=0; j<nfree; j++) z[i] += H[i*nfree+j] * y[j];
         w += y[i]*z[i];  v += y[i]*s[i];
      }
      if (fabs(w)<small || fabs(v)<small)  { identity(H,nfree); fail=1; continue;}
      For (i,nfree)  For (j,nfree)  
         H[i*nfree+j] += s[i]*s[j]/v - z[i]*z[j]/w;
#else /* BFGS */
      for (i=0,w=v=0.; i<nfree; i++) {
         for (j=0,z[i]=0.; j<nfree; j++) z[i]+=H[i*nfree+j]*y[j];
         w+=y[i]*z[i];    v+=y[i]*s[i];
      }
      if (fabs(v)<small)   { identity(H,nfree); fail=1; continue; }
      For (i,nfree)  For (j,nfree)
         H[i*nfree+j] += ((1+w/v)*s[i]*s[j]-z[i]*s[j]-s[i]*z[j])/v;
#endif

   }    /* for (Iround,maxround)  */

   /* try to remove this after updating LineSearch2() */
   *f=LK_BFGS_from_CODEML(tree,x,n);
  
   if (Iround==maxround) {

      return(-1);
   }
   if(nfree==n) { 
      xtoy(H, space, n*n);  /* H has variance matrix, or inverse of Hessian */
      return(1);
   }
   return(0);
}
int AlwaysCenter=0;
phydbl Small_Diff=.5e-6;
/*********************************************************/

int gradientB (int n, phydbl x[], phydbl f0, phydbl g[], 
    t_tree * tree, phydbl space[], int xmark[])
{
/* f0=fun(x) is always provided.
   xmark=0: central; 1: upper; -1: down
*/
   int i,j;
   phydbl *x0=space, *x1=space+n, eh0=Small_Diff, eh;  /* eh0=1e-6 || 1e-7 */

   For(i,n) {
      eh=eh0*(fabs(x[i])+1);
      if (xmark[i]==0 && (AlwaysCenter || SIZEp<1)) {   /* central */
         For (j, n)  x0[j]=x1[j]=x[j];
         eh=pow(eh,.67);  x0[i]-=eh; x1[i]+=eh;
         g[i]=(LK_BFGS_from_CODEML(tree,x1,n)-LK_BFGS_from_CODEML(tree,x0,n))/(eh*2.0);
      }
      else  {                         /* forward or backward */
         For (j, n)  x1[j]=x[j];
         if (xmark[i]) eh*=-xmark[i];
         x1[i] += eh;
         g[i]=(LK_BFGS_from_CODEML(tree,x1,n)-f0)/eh;
      }
   }
   return(0);
}

/*********************************************************/
phydbl innerp (phydbl x[], phydbl y[], int n)
{ int i; phydbl t=0;  for(i=0; i<n; i++)  t += x[i]*y[i];  return(t); }
/*********************************************************/
int xtoy (phydbl x[], phydbl y[], int n)
{ int i; for (i=0; i<n; i++) {y[i]=x[i];} return(0); }
/*********************************************************/
phydbl LineSearch2 (t_tree * tree, phydbl *f, phydbl x0[], 
       phydbl p[], phydbl step, phydbl limit, phydbl e, phydbl space[], int n)
{
/* linear search using quadratic interpolation 
   from x0[] in the direction of p[],
                x = x0 + a*p        a ~(0,limit)
   returns (a).    *f: f(x0) for input and f(x) for output

   x0[n] x[n] p[n] space[n]

   adapted from Wolfe M. A.  1978.  Numerical methods for unconstrained
   optimization: An introduction.  Van Nostrand Reinhold Company, New York.
   pp. 62-73.
   step is used to find the bracket and is increased or reduced as necessary, 
   and is not terribly important.
*/
   int ii=0, maxround=10, status;
   // int i, nsymb=0;
   phydbl *x=space, factor=4, small=1e-10, smallgapa=0.2;
   phydbl a0,a1,a2,a3,a4=-1,a5,a6, f0,f1,f2,f3,f4=-1,f5,f6;

/* look for bracket (a1, a2, a3) with function values (f1, f2, f3)
   step length step given, and only in the direction a>=0
*/

  
   if (step<=0 || limit<small || step>=limit) {
     
      return (0);
   }
   a0=a1=0; f1=f0=*f;
   a2=a0+step; f2=fun_LineSearch(a2, tree,x0,p,x,n);
   if (f2>f1) {  /* reduce step length so the algorithm is decreasing */
      for (; ;) {
         step/=factor;
         if (step<small) return (0);
         a3=a2;    f3=f2;
         a2=a0+step;  f2=fun_LineSearch(a2, tree,x0,p,x,n);
         if (f2<=f1) break;
        
      }
   }
   else {       /* step length is too small? */
      for (; ;) {
         step*=factor;
         if (step>limit) step=limit;
         a3=a0+step;  f3=fun_LineSearch(a3, tree,x0,p,x,n);
         if (f3>=f2) break;

        
         a1=a2; f1=f2;    a2=a3; f2=f3;
         if (step>=limit) {
            
            *f=f3; return(a3);
         }
      }
   }

   /* iteration by quadratic interpolation, fig 2.2.9-10 (pp 71-71) */
   for (ii=0; ii<maxround; ii++) {
      /* a4 is the minimum from the parabola over (a1,a2,a3)  */
      a4 = (a2-a3)*f1+(a3-a1)*f2+(a1-a2)*f3;
      if(fabs(a4)>1e-100) 
         a4 = ((a2*a2-a3*a3)*f1+(a3*a3-a1*a1)*f2+(a1*a1-a2*a2)*f3)/(2*a4);
      if (a4>a3 || a4<a1) {   /* out of range */
         a4=(a1+a2)/2;
         status='N';
      }
      else {
         if((a4<=a2 && a2-a4>smallgapa*(a2-a1)) || (a4>a2 && a4-a2>smallgapa*(a3-a2)))
            status='Y';
         else 
            status='C';
      }
      f4 = fun_LineSearch(a4, tree,x0,p,x,n);
      
      if (fabs(f2-f4)<e*(1+fabs(f2))) {
         
         break;
      }

      /* possible multiple local optima during line search */
      
      if (a4<=a2) {    /* fig 2.2.10 */
         if (a2-a4>smallgapa*(a2-a1)) {
            if (f4<=f2) { a3=a2; a2=a4;  f3=f2; f2=f4; }
            else        { a1=a4; f1=f4; }
         }
         else {
            if (f4>f2) {
               a5=(a2+a3)/2; f5=fun_LineSearch(a5, tree,x0,p,x,n);
               if (f5>f2) { a1=a4; a3=a5;  f1=f4; f3=f5; }
               else       { a1=a2; a2=a5;  f1=f2; f2=f5; }
            }
            else {
               a5=(a1+a4)/2; f5=fun_LineSearch(a5, tree,x0,p,x,n);
               if (f5>=f4)
                  { a3=a2; a2=a4; a1=a5;  f3=f2; f2=f4; f1=f5; }
               else {
                  a6=(a1+a5)/2; f6=fun_LineSearch(a6, tree,x0,p,x,n);
                  if (f6>f5)
                       { a1=a6; a2=a5; a3=a4;  f1=f6; f2=f5; f3=f4; }
                  else { a2=a6; a3=a5; f2=f6; f3=f5; }
               }
            }
         }
      }
      else {                     /* fig 2.2.9 */
         if (a4-a2>smallgapa*(a3-a2)) {
            if (f2>=f4) { a1=a2; a2=a4;  f1=f2; f2=f4; }
            else        { a3=a4; f3=f4; }
         }
         else {
            if (f4>f2) {
               a5=(a1+a2)/2; f5=fun_LineSearch(a5, tree,x0,p,x,n);
               if (f5>f2) { a1=a5; a3=a4;  f1=f5; f3=f4; }
               else       { a3=a2; a2=a5;  f3=f2; f2=f5; }
            }
            else {
               a5=(a3+a4)/2; f5=fun_LineSearch(a5, tree,x0,p,x,n);
               if (f5>=f4)
                  { a1=a2; a2=a4; a3=a5;  f1=f2; f2=f4; f3=f5; }
               else {
                  a6=(a3+a5)/2; f6=fun_LineSearch(a6, tree,x0,p,x,n);
                  if (f6>f5)
                      { a1=a4; a2=a5; a3=a6;  f1=f4; f2=f5; f3=f6; }
                  else { a1=a5; a2=a6;  f1=f5; f2=f6; }
               }
            }
         }
      }
   }

   if (f2>f0 && f4>f0)  a4=0;
   if (f2<=f4)  { *f=f2; a4=a2; }
   else         *f=f4;
   

   return (a4);
}

/*********************************************************/
phydbl fun_LineSearch (phydbl t, t_tree *tree, phydbl x0[], phydbl p[], phydbl x[], int n)
{  
  int i;   
  For (i,n) x[i]=x0[i] + t*p[i];   
  return( LK_BFGS_from_CODEML(tree,x,n) ); 
}
/*********************************************************/
int zero (phydbl x[], int n)
{ int i; for(i=0; i<n; i++) x[i]=0; return (0);}
/*********************************************************/
int identity (phydbl x[], int n)
{ int i,j;  for(i=0; i<n; i++)  { for(j=0; j<n; j++)   x[i*n+j]=0;  x[i*n+i]=1; }  return (0); }
/*********************************************************/
phydbl norm (phydbl x[], int n)
{ int i; phydbl t=0;  for(i=0; i<n; i++)  t += x[i]*x[i];  return sqrt(t); }
/*********************************************************/
int H_end (phydbl x0[], phydbl x1[], phydbl f0, phydbl f1,
    phydbl e1, phydbl e2, int n)
/*   Himmelblau termination rule.   return 1 for stop, 0 otherwise.
*/
{
   phydbl r;
   if((r=norm(x0,n))<e2)
      r=1;
   r*=e1;
   if(distance(x1,x0,n)>=r)
      return(0);
   r=fabs(f0);  if(r<e2) r=1;     
   r*=e1;
   if(fabs(f1-f0)>=r) 
      return(0);
   return (1);
}
/*********************************************************/
phydbl distance (phydbl x[], phydbl y[], int n)
{  int i; phydbl t=0;
   for (i=0; i<n; i++) t += square(x[i]-y[i]);
   return(sqrt(t));
}
/*********************************************************/
phydbl Sum (phydbl x[], int n)
{ 
  int i; 
  phydbl t=0;  
  for(i=0; i<n; i++) t += x[i];    
  return(t); 
}
/*********************************************************/
