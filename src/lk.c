/*
IgPhyML: a program that computes maximum likelihood phylogenies under
non-reversible codon models designed for antibody lineages.

Copyright (C) Kenneth B Hoehn. Sept 2016 onward.

built upon

codonPHYML: a program that  computes maximum likelihood phylogenies from
CODON homologous sequences.

Copyright (C) Marcelo Serrano Zanetti. Oct 2010 onward.

built upon

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "lk.h"
extern int senseCodons[64];
extern int indexSenseCodons[64];
extern phydbl ecmK07freq[61];
extern phydbl ecmS05freq[64];
/* int    LIM_SCALE; */
/* phydbl LIM_SCALE_VAL; */
/* phydbl BIG; */
/* phydbl SMALL; */

/*********************************************************/

void Init_Tips_At_One_Site_Nucleotides_Float(char state, int pos, phydbl *p_lk)
{
  switch(state)
    {
    case 'A' : p_lk[pos+0]=1.; p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+3]=.0;
      break;
    case 'C' : p_lk[pos+1]=1.; p_lk[pos+0]=p_lk[pos+2]=p_lk[pos+3]=.0;
      break;
    case 'G' : p_lk[pos+2]=1.; p_lk[pos+1]=p_lk[pos+0]=p_lk[pos+3]=.0;
      break;
    case 'T' : p_lk[pos+3]=1.; p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+0]=.0;
      break;
    case 'U' : p_lk[pos+3]=1.; p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+0]=.0;
      break;
    case 'M' : p_lk[pos+0]=p_lk[pos+1]=1.; p_lk[pos+2]=p_lk[pos+3]=.0;
      break;
    case 'R' : p_lk[pos+0]=p_lk[pos+2]=1.; p_lk[pos+1]=p_lk[pos+3]=.0;
      break;
    case 'W' : p_lk[pos+0]=p_lk[pos+3]=1.; p_lk[pos+1]=p_lk[pos+2]=.0;
      break;
    case 'S' : p_lk[pos+1]=p_lk[pos+2]=1.; p_lk[pos+0]=p_lk[pos+3]=.0;
      break;
    case 'Y' : p_lk[pos+1]=p_lk[pos+3]=1.; p_lk[pos+0]=p_lk[pos+2]=.0;
      break;
    case 'K' : p_lk[pos+2]=p_lk[pos+3]=1.; p_lk[pos+0]=p_lk[pos+1]=.0;
      break;
    case 'B' : p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+3]=1.; p_lk[pos+0]=.0;
      break;
    case 'D' : p_lk[pos+0]=p_lk[pos+2]=p_lk[pos+3]=1.; p_lk[pos+1]=.0;
      break;
    case 'H' : p_lk[pos+0]=p_lk[pos+1]=p_lk[pos+3]=1.; p_lk[pos+2]=.0;
      break;
    case 'V' : p_lk[pos+0]=p_lk[pos+1]=p_lk[pos+2]=1.; p_lk[pos+3]=.0;
      break;
    case 'N' : case 'X' : case '?' : case 'O' : case '-' :
      p_lk[pos+0]=p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+3]=1.;break;
    default :
      {
	PhyML_Printf("\n. Unknown character state : %c\n",state);
	Exit("\n. Init failed (check the data type)\n");
	break;
      }
    }
}

/*********************************************************/

void Init_Tips_At_One_Site_Nucleotides_Int(char state, int pos, short int *p_pars)
{
  switch(state)
    {
    case 'A' : p_pars[pos+0]=1; p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+3]=0;
      break;
    case 'C' : p_pars[pos+1]=1; p_pars[pos+0]=p_pars[pos+2]=p_pars[pos+3]=0;
      break;
    case 'G' : p_pars[pos+2]=1; p_pars[pos+1]=p_pars[pos+0]=p_pars[pos+3]=0;
      break;
    case 'T' : p_pars[pos+3]=1; p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+0]=0;
      break;
    case 'U' : p_pars[pos+3]=1; p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+0]=0;
      break;
    case 'M' : p_pars[pos+0]=p_pars[pos+1]=1; p_pars[pos+2]=p_pars[pos+3]=0;
      break;
    case 'R' : p_pars[pos+0]=p_pars[pos+2]=1; p_pars[pos+1]=p_pars[pos+3]=0;
      break;
    case 'W' : p_pars[pos+0]=p_pars[pos+3]=1; p_pars[pos+1]=p_pars[pos+2]=0;
      break;
    case 'S' : p_pars[pos+1]=p_pars[pos+2]=1; p_pars[pos+0]=p_pars[pos+3]=0;
      break;
    case 'Y' : p_pars[pos+1]=p_pars[pos+3]=1; p_pars[pos+0]=p_pars[pos+2]=0;
      break;
    case 'K' : p_pars[pos+2]=p_pars[pos+3]=1; p_pars[pos+0]=p_pars[pos+1]=0;
      break;
    case 'B' : p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+3]=1; p_pars[pos+0]=0;
      break;
    case 'D' : p_pars[pos+0]=p_pars[pos+2]=p_pars[pos+3]=1; p_pars[pos+1]=0;
      break;
    case 'H' : p_pars[pos+0]=p_pars[pos+1]=p_pars[pos+3]=1; p_pars[pos+2]=0;
      break;
    case 'V' : p_pars[pos+0]=p_pars[pos+1]=p_pars[pos+2]=1; p_pars[pos+3]=0;
      break;
    case 'N' : case 'X' : case '?' : case 'O' : case '-' :
      p_pars[pos+0]=p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+3]=1;break;
    default :
      {
	PhyML_Printf("\n. Unknown character state : %c\n",state);
	Exit("\n. Init failed (check the data type)\n");
	break;
      }
    }
}

/*********************************************************/



void Init_Tips_At_One_Site_AA_Float(char aa, int pos, phydbl *p_lk)
{
  int i;

  For(i,20) p_lk[pos+i] = .0;

  switch(aa){
  case 'A' : p_lk[pos+0]= 1.; break;/* Alanine */
  case 'R' : p_lk[pos+1]= 1.; break;/* Arginine */
  case 'N' : p_lk[pos+2]= 1.; break;/* Asparagine */
  case 'D' : p_lk[pos+3]= 1.; break;/* Aspartic acid */
  case 'C' : p_lk[pos+4]= 1.; break;/* Cysteine */
  case 'Q' : p_lk[pos+5]= 1.; break;/* Glutamine */
  case 'E' : p_lk[pos+6]= 1.; break;/* Glutamic acid */
  case 'G' : p_lk[pos+7]= 1.; break;/* Glycine */
  case 'H' : p_lk[pos+8]= 1.; break;/* Histidine */
  case 'I' : p_lk[pos+9]= 1.; break;/* Isoleucine */
  case 'L' : p_lk[pos+10]=1.; break;/* Leucine */
  case 'K' : p_lk[pos+11]=1.; break;/* Lysine */
  case 'M' : p_lk[pos+12]=1.; break;/* Methionine */
  case 'F' : p_lk[pos+13]=1.; break;/* Phenylalanin */
  case 'P' : p_lk[pos+14]=1.; break;/* Proline */
  case 'S' : p_lk[pos+15]=1.; break;/* Serine */
  case 'T' : p_lk[pos+16]=1.; break;/* Threonine */
  case 'W' : p_lk[pos+17]=1.; break;/* Tryptophan */
  case 'Y' : p_lk[pos+18]=1.; break;/* Tyrosine */
  case 'V' : p_lk[pos+19]=1.; break;/* Valine */

  case 'B' : p_lk[pos+2]= 1.; break;/* Asparagine */
  case 'Z' : p_lk[pos+5]= 1.; break;/* Glutamine */

  case 'X' : case '?' : case '-' : For(i,20) p_lk[pos+i] = 1.; break;
  default :
    {
      PhyML_Printf("\n. Unknown character state : %c\n",aa);
      Exit("\n. Init failed (check the data type)\n");
      break;
    }
  }
}

/*********************************************************/

void Init_Tips_At_One_Site_AA_Int(char aa, int pos, short int *p_pars)
{
  int i;

  For(i,20) p_pars[pos+i] = .0;

  switch(aa){
  case 'A' : p_pars[pos+0]  = 1; break;/* Alanine */
  case 'R' : p_pars[pos+1]  = 1; break;/* Arginine */
  case 'N' : p_pars[pos+2]  = 1; break;/* Asparagine */
  case 'D' : p_pars[pos+3]  = 1; break;/* Aspartic acid */
  case 'C' : p_pars[pos+4]  = 1; break;/* Cysteine */
  case 'Q' : p_pars[pos+5]  = 1; break;/* Glutamine */
  case 'E' : p_pars[pos+6]  = 1; break;/* Glutamic acid */
  case 'G' : p_pars[pos+7]  = 1; break;/* Glycine */
  case 'H' : p_pars[pos+8]  = 1; break;/* Histidine */
  case 'I' : p_pars[pos+9]  = 1; break;/* Isoleucine */
  case 'L' : p_pars[pos+10] = 1; break;/* Leucine */
  case 'K' : p_pars[pos+11] = 1; break;/* Lysine */
  case 'M' : p_pars[pos+12] = 1; break;/* Methionine */
  case 'F' : p_pars[pos+13] = 1; break;/* Phenylalanin */
  case 'P' : p_pars[pos+14] = 1; break;/* Proline */
  case 'S' : p_pars[pos+15] = 1; break;/* Serine */
  case 'T' : p_pars[pos+16] = 1; break;/* Threonine */
  case 'W' : p_pars[pos+17] = 1; break;/* Tryptophan */
  case 'Y' : p_pars[pos+18] = 1; break;/* Tyrosine */
  case 'V' : p_pars[pos+19] = 1; break;/* Valine */

  case 'B' : p_pars[pos+2]  = 1; break;/* Asparagine */
  case 'Z' : p_pars[pos+5]  = 1; break;/* Glutamine */

  case 'X' : case '?' : case '-' : For(i,20) p_pars[pos+i] = 1; break;
  default :
    {
      PhyML_Printf("\n. Unknown character state : %c\n",aa);
      Exit("\n. Init failed (check the data type)\n");
      break;
    }
  }
}

/*********************************************************/

void Init_Tips_At_One_Site_Generic_Float(char *state, int ns, int state_len, int pos, phydbl *p_lk)
{
  int i;
  int state_int;
 
  For(i,ns) p_lk[pos+i] = 0.;

  if(Is_Ambigu(state,GENERIC,state_len)) For(i,ns) p_lk[pos+i] = 1.;
  else
    {
      char format[6];
      sprintf(format,"%%%dd",state_len);
      if(!sscanf(state,format,&state_int))
	{
	  PhyML_Printf("\n. state='%c'",state);
	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      if(state_int > ns)
	{
	  PhyML_Printf("\n. %s %d cstate: %.2s istate: %d state_len: %d.\n",__FILE__,__LINE__,state,state_int,state_len);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");	  
	}
      p_lk[pos+state_int] = 1.;
      /*       PhyML_Printf("\n. %s %d cstate: %.2s istate: %d state_len: %d ns: %d pos: %d",__FILE__,__LINE__,state,state_int,state_len,ns,pos); */
    }
}

/*********************************************************/

void Init_Tips_At_One_Site_Generic_Int(char *state, int ns, int state_len, int pos, short int *p_pars)
{
  int i;
  int state_int;

  For(i,ns) p_pars[pos+i] = 0;
  
  if(Is_Ambigu(state,GENERIC,state_len)) For(i,ns) p_pars[pos+i] = 1;
  else 
    {
      char format[6];
      sprintf(format,"%%%dd",state_len);
      if(!sscanf(state,format,&state_int))
	{
	  PhyML_Printf("\n. state='%c'",state);
	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      if(state_int > ns)
	{
	  PhyML_Printf("\n. %s %d cstate: %.2s istate: %d state_len: %d.\n",__FILE__,__LINE__,state,state_int,state_len);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");	  
	}
      p_pars[pos+state_int] = 1;
    }
}

/*********************************************************/

void Get_All_Partial_Lk_Scale(t_tree *tree, t_edge *b_fcus, t_node *a, t_node *d)
{
  if(d->tax) return;
  else Update_P_Lk(tree,b_fcus,d);
}

/*********************************************************/

void Post_Order_Lk(t_node *a, t_node *d, t_tree *tree)
{
  int i,dir;

  dir = -1;

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    Post_Order_Lk(d,d->v[i],tree);
	  else dir = i;
	}      
      Get_All_Partial_Lk_Scale(tree,d->b[dir],a,d);
    }
}


/*********************************************************/

void Pre_Order_Lk(t_node *a, t_node *d, t_tree *tree)
{
  int i;

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    {
	      Get_All_Partial_Lk_Scale(tree,d->b[i],d->v[i],d);
	      Pre_Order_Lk(d,d->v[i],tree);
	    }
	}
    }
}



/*********************************************************/
phydbl Lk(t_tree *tree)
{
  int br, n_patterns, n_edges;  
  n_edges=2*tree->n_otu-3;
  tree->old_lnL = tree->c_lnL;
 
  if(tree->rates && tree->rates->lk_approx == NORMAL)
  {
    tree->c_lnL = Dnorm_Multi_Given_InvCov_Det(tree->rates->u_cur_l, 
					       tree->rates->u_ml_l,
					       tree->rates->invcov,
					       tree->rates->covdet,
					       2*tree->n_otu-3,YES);
					       return tree->c_lnL;
  }
  
  n_patterns = tree->n_pattern;
  
  Set_Model_Parameters(tree->mod);
  #ifdef TIMES
  if((tree->rates) && (tree->rates->bl_from_rt)) RATES_Update_Cur_Bl(tree);
  if(tree->bl_from_node_stamps) TIMES_Bl_From_T(tree);
  #endif
 
  #if defined OMP || defined BLAS_OMP
   
  #pragma omp parallel for if(tree->io->datatype==CODON)
   
  #endif
   
  For(br,n_edges) Update_PMat_At_Given_Edge(tree->t_edges[br],tree);
  
  //added and modified by Ken to start likelihood calculation at specified node
  int startnode = 0;
  if(tree->mod->whichrealmodel == HLP17){
	  startnode = tree->mod->startnode;
  }
  Post_Order_Lk(tree->noeud[startnode],tree->noeud[startnode]->v[0],tree);

  if(tree->both_sides)
    Pre_Order_Lk(tree->noeud[startnode], tree->noeud[startnode]->v[0], tree);

  tree->c_lnL             = .0;
  tree->sum_min_sum_scale = .0;

  if(tree->mod->whichrealmodel==HLP17){
	  Fill_UPP_root(tree,tree->noeud[startnode]->b[0]);
  }

	tree->c_lnL             = .0;
	tree->sum_min_sum_scale = .0;

  For(tree->curr_site,n_patterns)
  {
	  if(tree->data->wght[tree->curr_site] > SMALL && tree->mod->whichrealmodel==HLP17) Lk_Core_UPP(tree->noeud[startnode]->b[0],tree,tree->noeud[startnode],tree->noeud[startnode]->b[0]->des_node);
	  else if(tree->data->wght[tree->curr_site] > SMALL) Lk_Core(tree->noeud[startnode]->b[0],tree);
  }

  Adjust_Min_Diff_Lk(tree);
  
  return tree->c_lnL;
}


void upAllPmats(t_tree *tree)
{
  int br, n_patterns, n_edges;
  n_edges=2*tree->n_otu-3;

  For(br,n_edges) Update_PMat_At_Given_Edge(tree->t_edges[br],tree);
}


/*********************************************************/
//Added by Ken 20/10/2016
void Get_UPP(t_node *a, t_node *d,  t_tree *tree)
{
  int i;

  if(d->tax) return;
  else{
	  //recurse
	  int desc[2];
	  int count=0;
      For(i,3){//figure out which nodes are the descendants
    	  if(d->v[i]->num != d->anc->num){
			  desc[count]=i;
			  count++;
    	  }
      }
      //fill in upp matrix
      if(a->num == tree->mod->startnode){//if next to the root
    	  Fill_UPP_root(tree,d->anc_edge);
      }
      Fill_UPP_single(tree,d->b[desc[0]]);
      Fill_UPP_single(tree,d->b[desc[1]]);

      Get_UPP(d,d->v[desc[0]],tree);
      Get_UPP(d,d->v[desc[1]],tree);
    }
}

void Fill_UPP_single(t_tree *tree,  t_edge *target)
{
/*	  anc
       |
	   |<- b
	   |
	   d
      / \
u_e->/   \ <-c_e
  	/     \
n_v1   n_v2
*/
//printf("about to do this\n");
t_node *d = target->anc_node;

  t_node *adjn, *targn;
  t_edge *b, *adje;
  phydbl padj_lk,ptar_lk;
  phydbl *p_lk,*p_lk_adj,*p_lk_tar;
  int *sum_scale, *sum_scale_v1, *sum_scale_v2, *sum_scale_adj, *sum_scale_tar, *sum_scale_upp;
  int sum_scale_adj_val, sum_scale_upp_val, sum_scale_tar_val, sum_scale_val;
  int i,j;
  int catg,site;
  int dir1,dir2,adjdir,tardir;
  int n_patterns;
  short int ambiguity_check_adj,ambiguity_check_tar;
  int state_adj,state_tar;
  int dim1, dim2, dim3;
  phydbl curr_scaler;
  int curr_scaler_pow, piecewise_scaler_pow;
  phydbl p_lk_lim_inf;
  phydbl smallest_p_lk;

  p_lk_lim_inf = (phydbl)P_LK_LIM_INF;

  dim1 = tree->mod->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;
  dim3 = tree->mod->ns * tree->mod->ns;

  state_adj = state_tar = -1;
  ambiguity_check_adj = ambiguity_check_tar = NO;
  sum_scale_adj_val = sum_scale_tar_val = 0;

  if(d->tax){
      PhyML_Printf("\n. t_node %d is a leaf...",d->num);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("\n");
    }

  n_patterns = tree->n_pattern;

  b=d->anc_edge;
  dir1=dir2=-1;
  For(i,3){
	  if(d->b[i] != b){
		  (dir1<0)?(dir1=i):(dir2=i);
	  }
  }

  if(target->num == d->b[dir1]->num){
	  adjdir=dir2;
	  tardir=dir1;
  }else if(target->num == d->b[dir2]->num){
	  adjdir=dir1;
	  tardir=dir2;
  }else{
	  printf("didn't work\n");
	  printf("Edge %d error as descendant of node d: %d b: %d target_des: %d | %d %d | %d %d %d\n",target->num,d->num,b->num,target->des_node->num,dir1,dir2, d->b[0]->num, d->b[1]->num, d->b[2]->num);
	  exit(EXIT_FAILURE);
  }

  adje = d->b[adjdir];
  adjn = d->v[adjdir];
  targn = d->v[tardir];

  if(d->num == b->left->num){
      p_lk = b->p_lk_left;
    }
  else{
      p_lk = b->p_lk_rght;
    }

  //upp scaling
  sum_scale = target->sum_scale_upp;
  sum_scale_upp = b->sum_scale_upp;

  if(d->num == d->b[adjdir]->left->num){
        p_lk_adj = d->b[adjdir]->p_lk_rght;
        sum_scale_adj = d->b[adjdir]->sum_scale_rght;
      }else{
        p_lk_adj = d->b[adjdir]->p_lk_left;
        sum_scale_adj = d->b[adjdir]->sum_scale_left;
      }

  if(d->num == d->b[tardir]->left->num){
       p_lk_tar = d->b[tardir]->p_lk_rght;
       sum_scale_tar = d->b[tardir]->sum_scale_rght;
     }else{
       p_lk_tar = d->b[tardir]->p_lk_left;
       sum_scale_tar = d->b[tardir]->sum_scale_left;
     }

  phydbl tarsum=0;

      int interest = 29;
  For(site,n_patterns){
      phydbl tsitesum=0.0;

      state_adj = -1;
      ambiguity_check_adj = NO;
      if(!tree->mod->s_opt->greedy){
	   /* n_v1 and n_v2 are tip nodes */
	      if(adjn->tax){
	      /* Is the state at this tip ambiguous? */
	    	  if(adjn->num == target->left->num){
	    		  printf("left taxa %d %s\n",b->num,d->name);
	    	      	 exit(EXIT_FAILURE);
	    	  }
	       ambiguity_check_adj = tree->data->c_seq[adjn->num]->is_ambigu[site];
	       if(ambiguity_check_adj == NO) state_adj = Get_State_From_P_Pars(adje->p_lk_tip_r,site*dim2,tree);
	      }
      }

      if(tree->mod->use_m4mod){
    	  ambiguity_check_adj = YES;
      }

      /* For all the rate classes */
      For(catg,tree->mod->n_catg){
	 smallest_p_lk  =  BIG;

	  /* For all the state at node d */
	  For(i,tree->mod->ns){
		  target->upp[site][i]=0.0;
	      /* n_v1 is a tip */
	      if((adjn->tax) && (!tree->mod->s_opt->greedy)){
		  if(ambiguity_check_adj == NO){
		      /* For the (non-ambiguous) state at node n_v1 */
			  target->upp[site][i] = d->b[adjdir]->bPmat_part[tree->mod->partIndex[site]][catg*dim3+i*dim2+state_adj]; // Ken 17/8/2016
		    }else{
		      For(j,tree->mod->ns){
		    	phydbl val = d->b[adjdir]->bPmat_part[tree->mod->partIndex[site]][catg*dim3+i*dim2+j]*(phydbl)adje->p_lk_tip_r[site*dim2+j];
		    	target->upp[site][i] += val;
		      }
		    }
		  }
	      /* n_v1 is an internal node */
	   else{
		    For(j,tree->mod->ns){
		    	phydbl val = d->b[adjdir]->bPmat_part[tree->mod->partIndex[site]][catg*dim3+i*dim2+j] * p_lk_adj[site*dim1+catg*dim2+j];
		    	target->upp[site][i] += val;
		    }
	   }
	      //add in upper likelihood of b
	   phydbl sum=0.0;
	   For(j,tree->mod->ns){ /* For all the states at node n_v2 */
          phydbl val = b->bPmat_part[tree->mod->partIndex[site]][catg*dim3+j*dim2+i]*b->upp[site][j];
          sum += val;
        }
	   target->upp[site][i] *= sum;


	     tsitesum += target->upp[site][i];


	   	   if(target->upp[site][i] < smallest_p_lk) smallest_p_lk = target->upp[site][i] ;
	   }//For(i,tree->mod->ns)

	   	  sum_scale_adj_val = (sum_scale_adj)?(sum_scale_adj[catg*n_patterns+site]):(0);
	  	  sum_scale_upp_val = (sum_scale_upp)?(sum_scale_upp[catg*n_patterns+site]):(0);
	  	  sum_scale[catg*n_patterns+site] = sum_scale_adj_val + sum_scale_upp_val;

	  	  if(smallest_p_lk < p_lk_lim_inf){
	  	      curr_scaler_pow = (int)(LOG(p_lk_lim_inf)-LOG(smallest_p_lk))/LOG2;
	  	      curr_scaler     = (phydbl)((unsigned long long)(1) << curr_scaler_pow);
	  	      sum_scale[catg*n_patterns+site] += curr_scaler_pow;
	  	      do{
	  		  piecewise_scaler_pow = MIN(curr_scaler_pow,63);
	  		  curr_scaler = (phydbl)((unsigned long long)(1) << piecewise_scaler_pow);
	  		  For(i,tree->mod->ns){
	  			  target->upp[site][i] *= curr_scaler;
	  		      if(target->upp[site][i] > BIG){
	  			  PhyML_Printf("\n. p_lk_lim_inf = %G smallest_p_lk = %G",p_lk_lim_inf,smallest_p_lk);
	  			  PhyML_Printf("\n. curr_scaler_pow = %d",curr_scaler_pow);
	  			  PhyML_Printf("\n. Err in file %s at line %d.",__FILE__,__LINE__);
	  			  Warn_And_Exit("\n");
	  			}
	  		  }
	  		  curr_scaler_pow -= piecewise_scaler_pow;
	  		}
	  	      while(curr_scaler_pow != 0);
	  	 }

	  }//For(catg,tree->mod->n_catg){
    }//For(n,sites
}

void Fill_UPP_root(t_tree *tree, t_edge *b)
{
/*	  anc
       |
	   |<- b
	   |
	   d
      / \
u_e->/   \ <-c_e
  	/     \
n_v1   n_v2
*/

	t_node *d = b->des_node;
	t_node *anc = b->anc_node;
  t_node *n_v1, *n_v2;
  phydbl p1_lk1,p2_lk2;
  phydbl *p_lk,*p_lk_v1,*p_lk_v2;
  int *sum_scale, *sum_scale_adj, *sum_scale_upp;
  int sum_scale_v1_val, sum_scale_v2_val;
  int i,j;
  int catg,site;
  int dir1,dir2,adjdir,tardir;
  int n_patterns;
  short int ambiguity_check_anc;
  int state_anc;
  int dim1, dim2, dim3;
  phydbl curr_scaler;
  int curr_scaler_pow, piecewise_scaler_pow;
  phydbl p_lk_lim_inf;
  phydbl smallest_p_lk;

  p_lk_lim_inf = (phydbl)P_LK_LIM_INF;

  dim1 = tree->mod->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;
  dim3 = tree->mod->ns * tree->mod->ns;

  state_anc = -1;
  ambiguity_check_anc = NO;

  sum_scale = b->sum_scale_upp;

  int sum_scale_adj_val;


  if(d->tax)
    {
      PhyML_Printf("\n. t_node %d is a leaf...",d->num);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("\n");
    }

  n_patterns = tree->n_pattern;

  /* Get the partial likelihood vectors on edge b and the two pendant
     edges (i.e., the two other edges connected to d) */
  if(d == b->left){
      p_lk = b->p_lk_left;

    }
  else{
      p_lk = b->p_lk_rght;

    }

  For(site,n_patterns){
      state_anc = -1;
      ambiguity_check_anc= NO;
      if(!tree->mod->s_opt->greedy){
	  /* n_v1 and n_v2 are tip nodes */
	  if(anc->tax){
	      /* Is the state at this tip ambiguous? */
	      ambiguity_check_anc = tree->data->c_seq[anc->num]->is_ambigu[site];
	      if(ambiguity_check_anc == NO) state_anc = Get_State_From_P_Pars(anc->b[0]->p_lk_tip_r,site*dim2,tree);
	    }
      }

      if(tree->mod->use_m4mod){
    	  ambiguity_check_anc = YES;
      }

      /* For all the rate classes */
      For(catg,tree->mod->n_catg){
	  smallest_p_lk  =  BIG;

	  /* For all the state at node d */
	  For(i,tree->mod->ns){
		  b->upp[site][i]=0.0;
	      /* n_v1 is a tip */
	      if((anc->tax) && (!tree->mod->s_opt->greedy)){
		    if(ambiguity_check_anc == NO){
		      /* For the (non-ambiguous) state at node n_v1 */
		    	phydbl state = 0.0;
		    	if(state_anc==i){state=1;}
			     b->upp[site][i] = state; // Ken 17/8/2016
		    }else{
		    	 b->upp[site][i] = anc->b[0]->p_lk_tip_r[site*tree->mod->ns+i];
		    }
		  }

	      //if not at root node, add in upper likelihood of b
	      if(anc->num != tree->mod->startnode){
	    	  printf("something weird happened\n");
	    	  exit(EXIT_FAILURE);
	      }
	   	   if(b->upp[site][i] < smallest_p_lk) smallest_p_lk = b->upp[site][i] ;
	    }
	  sum_scale[catg*n_patterns+site]=0;
    }
  }
}


void Get_Lhood(t_node *a, t_node *d,  t_tree *tree)
{
  int i;

  phydbl lk = Lk_At_Given_Edge(d->anc_edge,tree);
  if(d->tax) return;
  else{
	  //recurse
	  int desc[2];
	  int count=0;
      For(i,3){//figure out which nodes are the descendants
    	  if(d->v[i]->num != d->anc->num){
			  desc[count]=i;
			  count++;
    	  }
      }
      Get_Lhood(d,d->v[desc[0]],tree);
      Get_Lhood(d,d->v[desc[1]],tree);
    }
}


/*********************************************************/
//Updated by Ken 3/8/2016
phydbl Lk_At_Given_Edge(t_edge *b_fcus, t_tree *tree)
{
	if(b_fcus->left->tax){
	   PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	   Warn_And_Exit("");
	}

	if(tree->mod->whichrealmodel == HLP17){//Added by Ken 3/8/2016
	  Update_PMat_At_Given_Edge(b_fcus,tree);

	  if(tree->mod->s_opt->opt_topo){
		  if(b_fcus->anc_node->num != tree->mod->startnode){
			  Fill_UPP_single(tree,b_fcus);
		  }else{
			  Fill_UPP_root(tree,b_fcus);
		  }
	  }
	  tree->c_lnL             = .0;
	   tree->sum_min_sum_scale = .0;
	   For(tree->curr_site,tree->n_pattern){
	   	 Lk_Core_UPP(b_fcus,tree, b_fcus->anc_node, b_fcus->des_node);
	   }
	  return tree->c_lnL;
  }else{

  int n_patterns;
  n_patterns = tree->n_pattern;

  #ifdef TIMES
  if((tree->rates) && (tree->rates->bl_from_rt)) RATES_Update_Cur_Bl(tree);
  if(tree->bl_from_node_stamps) TIMES_Bl_From_T(tree);
  #endif
  Update_PMat_At_Given_Edge(b_fcus,tree);

  tree->c_lnL             = .0;
  tree->sum_min_sum_scale = .0;
  For(tree->curr_site,n_patterns)
  {
    if(tree->data->wght[tree->curr_site] > SMALL) Lk_Core(b_fcus,tree);

  }
  //printf("%d tree LK3: %lf \n",b_fcus->num,tree->c_lnL);

  return tree->c_lnL;
  }
}

/*********************************************************/
/* Core of the likelihood calculcation. Assume that the partial likelihoods on both
   sides of t_edge *b are up-to-date. Calculate the log-likelihood at one site.
*/


phydbl Lk_Core_UPP(t_edge *b, t_tree *tree, t_node *anc, t_node *d)
{
  phydbl log_site_lk;
  phydbl site_lk_cat, site_lk;
  int fact_sum_scale;
  phydbl max_sum_scale;
  phydbl sum;
  int ambiguity_check,state;
  int catg,ns,k,l,site;
  int dim1,dim2,dim3;
  int *sum_scale_down_cat,*sum_scale_upp_cat,*sum_scale_down,*sum_scale_upp;
  phydbl multiplier;
  int exponent, piecewise_exponent;
  phydbl tmp;
  phydbl logbig;
  phydbl inv_site_lk;

  logbig = LOG((phydbl)BIG);

  dim1 = tree->mod->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;
  dim3 = tree->mod->ns * tree->mod->ns;

  log_site_lk     = .0;
  site_lk         = .0;
  site_lk_cat     = .0;
  ambiguity_check = -1;
  state           = -1;
  site            = tree->curr_site;
  ns              = tree->mod->ns;

  phydbl *p_lk;

  if(d->num == b->left->num){
      p_lk = b->p_lk_left;
      sum_scale_down = b->sum_scale_left;
      sum_scale_down_cat = b->sum_scale_left_cat;
     if((d->tax) && (!tree->mod->s_opt->greedy)){
    	 if(site==16){printf("left taxa %d %s\n",b->num,d->name);}
    	 exit(EXIT_FAILURE);
        ambiguity_check = tree->data->c_seq[d->num]->is_ambigu[site];
        if(!ambiguity_check)state = Get_State_From_P_Pars(b->p_lk_tip_l,site*dim2,tree);
      }
    }
  else{
      p_lk = b->p_lk_rght;
      sum_scale_down = b->sum_scale_rght;
      sum_scale_down_cat = b->sum_scale_rght_cat;
      if((d->tax) && (!tree->mod->s_opt->greedy)){
        ambiguity_check = tree->data->c_seq[d->num]->is_ambigu[site];
        if(!ambiguity_check)state = Get_State_From_P_Pars(b->p_lk_tip_r,site*dim2,tree);
      }
    }
  sum_scale_upp_cat = b->sum_scale_upp_cat;
  sum_scale_upp = b->sum_scale_upp;

  if(tree->mod->use_m4mod) ambiguity_check = 1;


  #if defined OMP || defined BLAS_OMP

  #pragma omp parallel for if(tree->mod->n_w_catg>1) private(site_lk_cat,sum,k,l)

  #endif


  For(catg,tree->mod->n_catg){
      site_lk_cat = .0;
      /* b is an external edge */
      if((d->tax) && (!tree->mod->s_opt->greedy)){
          /* If the character observed at the tip is NOT ambiguous: ns x 1 terms to consider */
          if(!ambiguity_check){
              sum = .0;
              For(l,ns){
                phydbl val = b->bPmat_part[tree->mod->partIndex[site]][catg*dim3+l*dim2+state]*b->upp[site][l];
                sum+=val;
              }
              site_lk_cat = sum;
            }
          /* If the character observed at the tip is ambiguous: ns x ns terms to consider */
          else{
              For(k,ns){
                  sum = .0;
                //  if(b->p_lk_tip_r[site*dim2+l] > .0){
                      For(l,ns){
                        phydbl val = b->bPmat_part[tree->mod->partIndex[site]][catg*dim3+k*dim2+l]*b->p_lk_tip_r[site*dim2+l]*b->upp[site][k];
                        sum+=val;
                      }
                      site_lk_cat += sum;
                    }
            //    }
            }
        }
      /* b is an internal edge: ns x ns terms to consider */
      else{
          For(k,ns){
              sum = .0;
         //     if(b->p_lk_rght[site*dim1+catg*dim2+k] > .0){
                  For(l,ns){ 
                    	phydbl val = b->bPmat_part[tree->mod->partIndex[site]][catg*dim3+k*dim2+l]*p_lk[site*dim1+catg*dim2+l]*b->upp[site][k];//Modified by Ken 17/8/2016
                      sum+=val;
                    }
                  site_lk_cat+=sum;
                }
        //    }
        }
      tree->site_lk_cat[catg] = site_lk_cat;
    }

  max_sum_scale =  (phydbl)BIG;
  For(catg,tree->mod->n_catg){

      sum_scale_down_cat[catg] =
        (sum_scale_down)?
        (sum_scale_down[catg*tree->n_pattern+site]):
        (0.0);

      sum_scale_upp_cat[catg] =
        (sum_scale_upp)?
        (sum_scale_upp[catg*tree->n_pattern+site]):
        (0.0);

      sum = sum_scale_down_cat[catg] + sum_scale_upp_cat[catg];
      if(sum < .0){
    	  PhyML_Printf("\n. sum = %G",sum);
    	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
    	  Warn_And_Exit("\n");
      }
      tmp = sum + (logbig - LOG(tree->site_lk_cat[catg]))/(phydbl)LOG2;
      if(tmp < max_sum_scale) max_sum_scale = tmp; /* min of the maxs */
    }

  fact_sum_scale = (int)(max_sum_scale / 2);

  // Apply scaling factors
  For(catg,tree->mod->n_catg){
      exponent = -(sum_scale_down_cat[catg]+sum_scale_upp_cat[catg])+fact_sum_scale;
      site_lk_cat = tree->site_lk_cat[catg];
      if(exponent > 0){
    	  do{
    		  piecewise_exponent = MIN(exponent,63);
    		  multiplier = (phydbl)((unsigned long long)(1) << piecewise_exponent);
    		  site_lk_cat *= multiplier;
    		  exponent -= piecewise_exponent;
    	  }
    	  while(exponent != 0);
      }
      else{
    	  do{
    		  piecewise_exponent = MAX(exponent,-63);
    		  multiplier = 1. / (phydbl)((unsigned long long)(1) << -piecewise_exponent);
    		  site_lk_cat *= multiplier;
    		  exponent -= piecewise_exponent;
    	  }
    	  while(exponent != 0);
      }

      if(isinf(site_lk_cat)){
    	  PhyML_Printf("\n+ site=%4d site_lk_cat=%G sum_scale=%d max=%G fact=%d expo=%d dbl=%G",
		       tree->curr_site,
		       tree->site_lk_cat[catg],
		       sum_scale_down_cat[catg]+sum_scale_upp_cat[catg],
		       max_sum_scale,
		       fact_sum_scale,
		       -(sum_scale_down_cat[catg]+sum_scale_upp_cat[catg])+fact_sum_scale,
		       (double)tree->site_lk_cat[catg] * pow(2.,-(sum_scale_down_cat[catg]+sum_scale_upp_cat[catg])+fact_sum_scale));
    	  site_lk_cat = BIG / 10;
      }

      if(site_lk_cat < SMALL){
    	  site_lk_cat = SMALL;
    	  printf("site lk = 0 %d %d %lf %d %d %d.\nSetting underflow to DBL_MIN. See Manual.\n",anc->num,d->num,site_lk_cat,site,sum_scale_upp_cat[catg],sum_scale_down_cat[catg]);
    	  printf("%lf %lf",tree->mod->omega_part[0],tree->mod->kappa);
       	  printf("\n");
      }
      tree->site_lk_cat[catg] = site_lk_cat;
    }


  site_lk = .0;
  if(tree->io->datatype!=CODON ) For(catg,tree->mod->n_catg) site_lk += tree->site_lk_cat[catg]*tree->mod->gamma_r_proba[catg];
  else //!< Added by Marcelo.
  {
    if(tree->mod->n_catg>1 && tree->mod->n_w_catg==1) For(catg,tree->mod->n_catg) site_lk += tree->site_lk_cat[catg]*tree->mod->gamma_r_proba[catg];
    else For(catg,tree->mod->n_catg) site_lk += tree->site_lk_cat[catg]*tree->mod->prob_omegas[catg];
  }

  /* The substitution model does include invariable sites */
   if(tree->mod->invar)
     {
       /* The site is invariant */
       if(tree->data->invar[site] > -0.5)
         {
 	  /* Multiply P(D|r=0) by 2^(fact_sum_scale) */

 	  inv_site_lk = tree->mod->pi[tree->data->invar[site]];
 	  exponent = fact_sum_scale;
 	  do
 	    {
 	      piecewise_exponent = MIN(exponent,63);
 	      multiplier = (phydbl)((unsigned long long)(1) << piecewise_exponent);
 	      inv_site_lk *= multiplier;
 	      exponent -= piecewise_exponent;
 	    }
 	  while(exponent != 0);
 	  inv_site_lk *= tree->mod->pinvar;
 	  /* Update the value of site_lk */
 	  site_lk *= (1. - tree->mod->pinvar);
 	  site_lk += inv_site_lk;
         }
       else
         {
           /* Same formula as above with P(D | subs, rate = 0) = 0 */
 	  site_lk *= (1. - tree->mod->pinvar);
         }
     }

  log_site_lk = LOG(site_lk) - (phydbl)LOG2 * fact_sum_scale;
  For(catg,tree->mod->n_catg) tree->log_site_lk_cat[catg][site] = LOG(tree->site_lk_cat[catg]) - (phydbl)LOG2 * fact_sum_scale;
  if(isinf(log_site_lk) || isnan(log_site_lk))
    {
      PhyML_Printf("\n. site = %d",site);
      PhyML_Printf("\n. invar = %f",tree->data->invar[site]);
      PhyML_Printf("\n. scale_down = %d scale_upp = %d %lf",sum_scale_down[0],sum_scale_upp[0],fact_sum_scale);
      PhyML_Printf("\n. Lk = %G LOG(Lk) = %f < %G",site_lk,log_site_lk,-BIG);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      PhyML_Printf("\n. %d %d %d\n\n",b->num,d->num,anc->num);
      Warn_And_Exit("\n");
    }
  tree->cur_site_lk[site] = log_site_lk;
  /* Multiply log likelihood by the number of times this site pattern is found in the data */
  tree->c_lnL_sorted[site] = tree->data->wght[site]*log_site_lk;
  tree->c_lnL += tree->data->wght[site]*log_site_lk;
  return log_site_lk;
}

#ifndef USE_OLD_LK
phydbl Lk_Core(t_edge *b, t_tree *tree)
{
  phydbl log_site_lk;
  phydbl site_lk_cat, site_lk;
  int sum_scale_left, sum_scale_rght;
  int fact_sum_scale;
  phydbl max_sum_scale;
  phydbl sum;
  int ambiguity_check,state;
  int catg,ns,k,l,site;
  int dim1,dim2,dim3;
  int *sum_scale_left_cat,*sum_scale_rght_cat;
  phydbl multiplier;
  int exponent, piecewise_exponent;
  phydbl tmp;
  phydbl logbig;
  phydbl inv_site_lk;

  logbig = LOG((phydbl)BIG);


  dim1 = tree->mod->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;
  dim3 = tree->mod->ns * tree->mod->ns;

  log_site_lk     = .0;
  site_lk         = .0;
  site_lk_cat     = .0;
  ambiguity_check = -1;
  state           = -1;
  site            = tree->curr_site;
  ns              = tree->mod->ns;



  if((b->rght->tax) && (!tree->mod->s_opt->greedy))
    {
      ambiguity_check = tree->data->c_seq[b->rght->num]->is_ambigu[site];
      if(!ambiguity_check) state = Get_State_From_P_Pars(b->p_lk_tip_r,site*dim2,tree);
    }

  sum_scale_left = .0;
  sum_scale_rght = .0;
  
  if(tree->mod->use_m4mod) ambiguity_check = 1;
  
  sum_scale_left_cat = b->sum_scale_left_cat;
  sum_scale_rght_cat = b->sum_scale_rght_cat;

  /* Actual likelihood calculation */
  /* For all classes of rates */
 
  #if defined OMP || defined BLAS_OMP 
  
  #pragma omp parallel for if(tree->mod->n_w_catg>1) private(site_lk_cat,sum,k,l) 
     
  #endif

 
  For(catg,tree->mod->n_catg)
     {

      site_lk_cat = .0;

      /* b is an external edge */
      if((b->rght->tax) && (!tree->mod->s_opt->greedy))
        {
          /* If the character observed at the tip is NOT ambiguous: ns x 1 terms to consider */
          if(!ambiguity_check)
            {
              sum = .0;
              For(l,ns)
                {
                  sum +=
                    b->bPmat_part[tree->mod->partIndex[site]][catg*dim3+state*dim2+l] * //Modified by Ken 17/8/2016
                    b->p_lk_left[site*dim1+catg*dim2+l];
                }
              site_lk_cat += sum * tree->mod->pi[state];
            }
          /* If the character observed at the tip is ambiguous: ns x ns terms to consider */
          else
            {
              For(k,ns)
                {
                  sum = .0;
                  if(b->p_lk_tip_r[site*dim2+k] > .0)
                    {
                      For(l,ns)
                        {
                          sum +=
                        	b->bPmat_part[tree->mod->partIndex[site]][catg*dim3+k*dim2+l] *			//Modified by Ken 17/8/2016
                            b->p_lk_left[site*dim1+catg*dim2+l];
                        }
                      site_lk_cat +=
                        sum *
                        tree->mod->pi[k] *
                        b->p_lk_tip_r[site*dim2+k];
                    }
                }
            }
        }
      /* b is an internal edge: ns x ns terms to consider */
      else
        {
          For(k,ns)
            {
              sum = .0;
              if(b->p_lk_rght[site*dim1+catg*dim2+k] > .0)
                {
                  For(l,ns)
                    {
                      sum +=
                    	b->bPmat_part[tree->mod->partIndex[site]][catg*dim3+k*dim2+l] * //Modified by Ken 17/8/2016
                        b->p_lk_left[site*dim1+catg*dim2+l];
                      //printf("lk internal!\n");
                    }
                  site_lk_cat +=
                    sum *
                    tree->mod->pi[k] *
                    b->p_lk_rght[site*dim1+catg*dim2+k];
                }
            }
        }
      tree->site_lk_cat[catg] = site_lk_cat;
    }
    
  max_sum_scale =  (phydbl)BIG;
  For(catg,tree->mod->n_catg)
    {
      sum_scale_left_cat[catg] =
        (b->sum_scale_left)?
        (b->sum_scale_left[catg*tree->n_pattern+site]):
        (0.0);

      sum_scale_rght_cat[catg] =
        (b->sum_scale_rght)?
        (b->sum_scale_rght[catg*tree->n_pattern+site]):
        (0.0);

      sum = sum_scale_left_cat[catg] + sum_scale_rght_cat[catg];

      if(sum < .0)
	{
	  PhyML_Printf("\n. sum = %G",sum);
	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  Warn_And_Exit("\n");
	}

      tmp = sum + (logbig - LOG(tree->site_lk_cat[catg]))/(phydbl)LOG2;
      if(tmp < max_sum_scale) max_sum_scale = tmp; /* min of the maxs */
    }

/*   fact_sum_scale = (int)((max_sum_scale + min_sum_scale) / 2); */


  fact_sum_scale = (int)(max_sum_scale / 2);


  /* Apply scaling factors */
  For(catg,tree->mod->n_catg)
    {
      exponent = -(sum_scale_left_cat[catg]+sum_scale_rght_cat[catg])+fact_sum_scale;

      site_lk_cat = tree->site_lk_cat[catg];

      if(exponent > 0)
	{
	  do
	    {
	      piecewise_exponent = MIN(exponent,63);
	      multiplier = (phydbl)((unsigned long long)(1) << piecewise_exponent);
	      site_lk_cat *= multiplier;
	      exponent -= piecewise_exponent;
	    }
	  while(exponent != 0);
	}
      else
	{
	  do
	    {
	      piecewise_exponent = MAX(exponent,-63);
	      multiplier = 1. / (phydbl)((unsigned long long)(1) << -piecewise_exponent);
	      site_lk_cat *= multiplier;
	      exponent -= piecewise_exponent;
	    }
	  while(exponent != 0);
	}

      if(isinf(site_lk_cat))
	{
	  PhyML_Printf("\n+ site=%4d site_lk_cat=%G sum_scale=%d max=%G fact=%d expo=%d dbl=%G",
		       tree->curr_site,
		       tree->site_lk_cat[catg],
		       sum_scale_left_cat[catg]+sum_scale_rght_cat[catg],
		       max_sum_scale,
		       fact_sum_scale,
		       -(sum_scale_left_cat[catg]+sum_scale_rght_cat[catg])+fact_sum_scale,
		       (double)tree->site_lk_cat[catg] * pow(2.,-(sum_scale_left_cat[catg]+sum_scale_rght_cat[catg])+fact_sum_scale));
	  site_lk_cat = BIG / 10;
	}

      //if(site_lk_cat < SMALL) site_lk_cat = 0;
      if(site_lk_cat < SMALL){
    	  site_lk_cat = SMALL;
    	  printf("site lk = 0 %lf %d %d %d.\nSetting underflow to DBL_MIN. Probably not a big problem during initial parameter searching..\n",site_lk_cat,site,sum_scale_left_cat[catg],sum_scale_rght_cat[catg]);
    	  printf("%lf %lf",tree->mod->omega_part[0],tree->mod->kappa);
       	  printf("\n");
      }
      tree->site_lk_cat[catg] = site_lk_cat;
    }


  site_lk = .0;
  if(tree->io->datatype!=CODON ) For(catg,tree->mod->n_catg) site_lk += tree->site_lk_cat[catg]*tree->mod->gamma_r_proba[catg];
  else //!< Added by Marcelo.
  {
    if(tree->mod->n_catg>1 && tree->mod->n_w_catg==1) For(catg,tree->mod->n_catg) site_lk += tree->site_lk_cat[catg]*tree->mod->gamma_r_proba[catg];
    else For(catg,tree->mod->n_catg) site_lk += tree->site_lk_cat[catg]*tree->mod->prob_omegas[catg];
  }

 /* The substitution model does include invariable sites */
  if(tree->mod->invar)
    {
      /* The site is invariant */
      if(tree->data->invar[site] > -0.5)
        {
	  /* Multiply P(D|r=0) by 2^(fact_sum_scale) */
	  inv_site_lk = tree->mod->pi[tree->data->invar[site]];
	  exponent = fact_sum_scale;
	  do
	    {
	      piecewise_exponent = MIN(exponent,63);
	      multiplier = (phydbl)((unsigned long long)(1) << piecewise_exponent);
	      inv_site_lk *= multiplier;
	      exponent -= piecewise_exponent;
	    }
	  while(exponent != 0);
	  inv_site_lk *= tree->mod->pinvar;
	  /* Update the value of site_lk */
	  site_lk *= (1. - tree->mod->pinvar);
	  site_lk += inv_site_lk;
        }
      else
        {
          /* Same formula as above with P(D | subs, rate = 0) = 0 */
	  site_lk *= (1. - tree->mod->pinvar);
        }
    }

  log_site_lk = LOG(site_lk) - (phydbl)LOG2 * fact_sum_scale;

  For(catg,tree->mod->n_catg) tree->log_site_lk_cat[catg][site] = LOG(tree->site_lk_cat[catg]) - (phydbl)LOG2 * fact_sum_scale;
  //HERE!!!
  if(isinf(log_site_lk) || isnan(log_site_lk))
    {
      PhyML_Printf("\n. site = %d",site);
      PhyML_Printf("\n. invar = %f",tree->data->invar[site]);
      PhyML_Printf("\n. scale_left = %d scale_rght = %d",sum_scale_left,sum_scale_rght);
      PhyML_Printf("\n. Lk = %G LOG(Lk) = %f < %G",site_lk,log_site_lk,-BIG);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("\n");
    }

  tree->cur_site_lk[site] = log_site_lk;

  /* Multiply log likelihood by the number of times this site pattern is found in the data */
  tree->c_lnL_sorted[site] = tree->data->wght[site]*log_site_lk;

  tree->c_lnL += tree->data->wght[site]*log_site_lk;

  return log_site_lk;
}


/*********************************************************/

/* Update partial likelihood on edge b on the side of b where
   node d lies.
*/
void Update_P_Lk(t_tree *tree, t_edge *b, t_node *d)
{
/*
       |
	   |<- b
	   |
	   d
      / \
   	 /   \
  	/     \
n_v1   n_v2
*/
  t_node *n_v1, *n_v2;
  phydbl p1_lk1,p2_lk2;
  phydbl *p_lk,*p_lk_v1,*p_lk_v2;
 // phydbl *Pij1,*Pij2;
  int *sum_scale, *sum_scale_v1, *sum_scale_v2;
  int sum_scale_v1_val, sum_scale_v2_val;
  int i,j;
  int catg,site;
  int dir1,dir2;
  int n_patterns;
  short int ambiguity_check_v1,ambiguity_check_v2;
  int state_v1,state_v2;
  int dim1, dim2, dim3;
  phydbl curr_scaler;
  int curr_scaler_pow, piecewise_scaler_pow;
  phydbl p_lk_lim_inf;
  phydbl smallest_p_lk;

  p_lk_lim_inf = (phydbl)P_LK_LIM_INF;
 
  dim1 = tree->mod->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;
  dim3 = tree->mod->ns * tree->mod->ns;

  state_v1 = state_v2 = -1;
  ambiguity_check_v1 = ambiguity_check_v2 = NO;
  sum_scale_v1_val = sum_scale_v2_val = 0;
  p1_lk1 = p2_lk2 = .0;

  if(d->tax)
    {
      PhyML_Printf("\n. t_node %d is a leaf...",d->num);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("\n");
    }

  n_patterns = tree->n_pattern;
  
  /* TO DO: Might be worth keeping these directions in memory instead of
     calculating them every time... */
  dir1=dir2=-1;
  For(i,3) if(d->b[i] != b) (dir1<0)?(dir1=i):(dir2=i);
  
  n_v1 = d->v[dir1];
  n_v2 = d->v[dir2];

  /* Get the partial likelihood vectors on edge b and the two pendant
     edges (i.e., the two other edges connected to d) */
  if(d == b->left)
    {
      p_lk = b->p_lk_left;
      sum_scale = b->sum_scale_left;
    }
  else
    {
      p_lk = b->p_lk_rght;
      sum_scale = b->sum_scale_rght;
    }
      
  if(d == d->b[dir1]->left)
    {
      p_lk_v1 = d->b[dir1]->p_lk_rght;
      sum_scale_v1 = d->b[dir1]->sum_scale_rght;
    }
  else
    {
      p_lk_v1 = d->b[dir1]->p_lk_left;
      sum_scale_v1 = d->b[dir1]->sum_scale_left;
    }
  
  if(d == d->b[dir2]->left)
    {
      p_lk_v2 = d->b[dir2]->p_lk_rght;
      sum_scale_v2 = d->b[dir2]->sum_scale_rght;
    }
  else
    {
      p_lk_v2 = d->b[dir2]->p_lk_left;
      sum_scale_v2 = d->b[dir2]->sum_scale_left;
    }
  

  phydbl **Ppart1 = d->b[dir1]->bPmat_part;
  phydbl **Ppart2 = d->b[dir2]->bPmat_part;
  
  /* For every site in the alignment */
  
  #if defined OMP || defined BLAS_OMP
  
  #pragma omp parallel for if(tree->io->datatype==CODON) private(state_v1, state_v2, ambiguity_check_v1, ambiguity_check_v2, smallest_p_lk, p1_lk1, catg, i, j, p2_lk2, sum_scale_v1_val, sum_scale_v2_val, curr_scaler_pow, curr_scaler, piecewise_scaler_pow )
  
  #endif
  
  For(site,n_patterns)
    {
      state_v1 = state_v2 = -1;
      ambiguity_check_v1 = ambiguity_check_v2 = NO;
      if(!tree->mod->s_opt->greedy)
	{
	  /* n_v1 and n_v2 are tip nodes */
	  if(n_v1->tax)
	    {
	      /* Is the state at this tip ambiguous? */
	      ambiguity_check_v1 = tree->data->c_seq[n_v1->num]->is_ambigu[site];
	      if(ambiguity_check_v1 == NO) state_v1 = Get_State_From_P_Pars(n_v1->b[0]->p_lk_tip_r,site*dim2,tree);
	    }
	      
	  if(n_v2->tax)
	    {
	      /* Is the state at this tip ambiguous? */
	      ambiguity_check_v2 = tree->data->c_seq[n_v2->num]->is_ambigu[site];
	      if(ambiguity_check_v2 == NO) state_v2 = Get_State_From_P_Pars(n_v2->b[0]->p_lk_tip_r,site*dim2,tree);

	    }
	}
      
      if(tree->mod->use_m4mod)
	{
	  ambiguity_check_v1 = YES;
	  ambiguity_check_v2 = YES;
	}

      /* For all the rate classes */
      For(catg,tree->mod->n_catg)
	{
	  smallest_p_lk  =  BIG;

	  /* For all the state at node d */
	  For(i,tree->mod->ns)
	    {
	      p1_lk1 = .0;
	      
	      /* n_v1 is a tip */
	      if((n_v1->tax) && (!tree->mod->s_opt->greedy))
		{
		  if(ambiguity_check_v1 == NO)
		    {
		      /* For the (non-ambiguous) state at node n_v1 */
		      p1_lk1 = d->b[dir1]->bPmat_part[tree->mod->partIndex[site]][catg*dim3+i*dim2+state_v1]; // Ken 17/8/2016
		    }
		  else
		    {
		      /* For all the states at node n_v1 */
		      For(j,tree->mod->ns)
			{
			  p1_lk1 +=d->b[dir1]->bPmat_part[tree->mod->partIndex[site]][catg*dim3+i*dim2+j] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+j]; //!!!!!
			  phydbl temp = (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+j];
			}
		    }
		}
	      /* n_v1 is an internal node */
	      else
		{
		  /* For the states at node n_v1 */
		  For(j,tree->mod->ns)
		    {
		      p1_lk1 += d->b[dir1]->bPmat_part[tree->mod->partIndex[site]][catg*dim3+i*dim2+j] * p_lk_v1[site*dim1+catg*dim2+j];//!!!!!!
		    }
		}
	      
	      p2_lk2 = .0;
	      
	      /* We do exactly the same as for node n_v1 but for node n_v2 this time.*/
	      /* n_v2 is a tip */
	      if((n_v2->tax) && (!tree->mod->s_opt->greedy))
		{
		  if(ambiguity_check_v2 == NO)
		    {
		      /* For the (non-ambiguous) state at node n_v2 */
		      p2_lk2 =  d->b[dir2]->bPmat_part[tree->mod->partIndex[site]][catg*dim3+i*dim2+state_v2];

		    }
		  else
		    {
		      /* For all the states at node n_v2 */
		      For(j,tree->mod->ns)
			{
			  p2_lk2 += d->b[dir2]->bPmat_part[tree->mod->partIndex[site]][catg*dim3+i*dim2+j] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+j]; //!!!!!
			  phydbl temp = (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+j];
			}
		    }
		}
	      /* n_v2 is an internal node */
	      else
		{
		  /* For all the states at node n_v2 */
		  For(j,tree->mod->ns)
		    {
		      p2_lk2 += d->b[dir2]->bPmat_part[tree->mod->partIndex[site]][catg*dim3+i*dim2+j] * p_lk_v2[site*dim1+catg*dim2+j];//!!!!!
		    }
		}
	      
	      p_lk[site*dim1+catg*dim2+i] = p1_lk1 * p2_lk2;	    

	      if(p_lk[site*dim1+catg*dim2+i] < smallest_p_lk) smallest_p_lk = p_lk[site*dim1+catg*dim2+i] ; 
	    }
	      
	  /* Current scaling values at that site */
 	  sum_scale_v1_val = (sum_scale_v1)?(sum_scale_v1[catg*n_patterns+site]):(0);
	  sum_scale_v2_val = (sum_scale_v2)?(sum_scale_v2[catg*n_patterns+site]):(0);
	  sum_scale[catg*n_patterns+site] = sum_scale_v1_val + sum_scale_v2_val;
	  
	  /* Scaling */
	  if(smallest_p_lk < p_lk_lim_inf)
	    {
	      curr_scaler_pow = (int)(LOG(p_lk_lim_inf)-LOG(smallest_p_lk))/LOG2;
	      curr_scaler     = (phydbl)((unsigned long long)(1) << curr_scaler_pow);

	      sum_scale[catg*n_patterns+site] += curr_scaler_pow;

	      do
		{
		  piecewise_scaler_pow = MIN(curr_scaler_pow,63);
		  curr_scaler = (phydbl)((unsigned long long)(1) << piecewise_scaler_pow);
		  For(i,tree->mod->ns)
		    {
		      p_lk[site*dim1+catg*dim2+i] *= curr_scaler;
		      
		      if(p_lk[site*dim1+catg*dim2+i] > BIG)
			{
			  PhyML_Printf("\n. p_lk_lim_inf = %G smallest_p_lk = %G",p_lk_lim_inf,smallest_p_lk);
			  PhyML_Printf("\n. curr_scaler_pow = %d",curr_scaler_pow);
			  PhyML_Printf("\n. Err in file %s at line %d.",__FILE__,__LINE__);
			  Warn_And_Exit("\n");
			}
		    }
		  curr_scaler_pow -= piecewise_scaler_pow;
		}
	      while(curr_scaler_pow != 0);
	    }
	}
    }
  if(tree->mod->s_opt->opt_topo && tree->mod->whichrealmodel == HLP17){
	  if(b->anc_node->num != tree->mod->startnode){
		  Fill_UPP_single(tree,b);
	  }else{
		  Fill_UPP_root(tree,b);
	  }
  }
}
#endif

/*********************************************************/

matrix *ML_Dist(calign *data, model *mod)
{
  int i,j,k,l;
  phydbl init;
  int n_catg;
  phydbl d_max,sum;
  matrix *mat;
  calign *twodata,*tmpdata;
  int state0, state1,len;
  phydbl *F;
  eigen *eigen_struct;
  time_t tbegin, tend; //! Added by Marcelo
  
  tmpdata         = (calign *)mCalloc(1,sizeof(calign));
  tmpdata->c_seq  = (align **)mCalloc(2,sizeof(align *));
  tmpdata->b_frq  = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));
  tmpdata->ambigu = (short int *)mCalloc(data->crunch_len,sizeof(short int));
  F               = (phydbl *)mCalloc(mod->ns*mod->ns,sizeof(phydbl ));
  eigen_struct    = (eigen *)Make_Eigen_Struct(mod);
  
  tmpdata->n_otu  = 2;
  
  tmpdata->crunch_len = data->crunch_len;
  tmpdata->init_len   = data->init_len;
  
  
  mat = NULL;

  For(i,mod->n_catg) /* Don't use the discrete gamma distribution */
  {
    mod->gamma_rr[i]      = 1.0;
    mod->gamma_r_proba[i] = 1.0;
  }
  
  n_catg = mod->n_catg;
  mod->n_catg = 1;
  
  time(&tbegin); //! Added by Marcelo.
  
  int progress = 0; //!Added by Marcelo.
  double loops = ((data->n_otu)*(data->n_otu)-(data->n_otu))/2.00;//!Added by Marcelo.
  For(j,data->n_otu-1)
  {
    tmpdata->c_seq[0]       = data->c_seq[j];
    tmpdata->c_seq[0]->name = data->c_seq[j]->name;
    tmpdata->wght           = data->wght;
    
    for(k=j+1;k<data->n_otu;k++)
    {
      
      tmpdata->c_seq[1]       = data->c_seq[k];
      tmpdata->c_seq[1]->name = data->c_seq[k]->name;
      twodata = Compact_Cdata(tmpdata,mod->io);
      For(l,mod->ns) twodata->b_frq[l] = data->b_frq[l];
      Check_Ambiguities(twodata,mod->io->datatype,mod->io->mod->state_len);
      
      Hide_Ambiguities(twodata);
      
      init = mat->dist[j][k];
      
      if((init > DIST_MAX-SMALL) || (init < .0)) init = 0.1;
      
      d_max = init;
      
      For(i,mod->ns*mod->ns) F[i]=.0;
      len = 0;
      For(l,twodata->c_seq[0]->len)
      {
	state0 = Assign_State(twodata->c_seq[0]->state+l*mod->io->mod->state_len,mod->io->datatype,mod->io->mod->state_len);
	state1 = Assign_State(twodata->c_seq[1]->state+l*mod->io->mod->state_len,mod->io->datatype,mod->io->mod->state_len);
	
	if((state0 > -1) && (state1 > -1))
	{
	  if(mod->io->datatype == CODON) //! Added by Marcelo.
	  { 
	    state0=indexSenseCodons[state0];
	    state1=indexSenseCodons[state1];
	  }
	  F[mod->ns*state0+state1] += (int)twodata->wght[l];
	  len += (int)twodata->wght[l];
	}
      }
      sum=0;
      For(i,mod->ns*mod->ns) sum += F[i];
           
      if(len > .0) {For(i,mod->ns*mod->ns) F[i] /= (phydbl)len;}
      
      sum = 0.;
      For(i,mod->ns*mod->ns) sum += F[i];
      
      
      if(sum < .001) d_max = -1.;
      else if((sum > 1. - .001) && (sum < 1. + .001)) Opt_Dist_F(&(d_max),F,mod);
      else
      {
	PhyML_Printf("\n. sum = %f\n",sum);
	PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	Exit("");
      }
      
      if(d_max >= DIST_MAX) d_max = DIST_MAX;
      
      
      /* Do not correct for dist < BL_MIN, otherwise Fill_Missing_Dist
      *  will not be called
      */
      
      mat->dist[j][k] = d_max;
      mat->dist[k][j] = mat->dist[j][k];
      Free_Cseq(twodata);
      progress++;
      if(!mod->io->quiet) PhyML_Printf("\r. Computing pairwise distances...%3.2f%c concluded.", progress*100.00/loops, '%');
    }
  }
  if(!mod->io->quiet) PhyML_Printf("\n");
  time(&tend); //! Added by Marcelo.
  
  if((tend-tbegin)>=1) //! Added by Marcelo.
  {
    PhyML_Printf("\n. Calculation of pairwise distances completed in %ld sec.\n",(long int)(tend-tbegin));
  }
  else
  {
    PhyML_Printf("\n. Calculation of pairwise distances completed in less than 1 sec.\n");
  }
  
  mod->n_catg = n_catg;
  
  free(tmpdata->ambigu);
  free(tmpdata->b_frq);
  free(tmpdata->c_seq);
  free(tmpdata);
  Free_Eigen(eigen_struct);
  free(F);
  
  return mat;
}
/*********************************************************/

void Unconstraint_Lk(t_tree *tree)
{
  int i;

  tree->unconstraint_lk = .0;

  For(i,tree->data->crunch_len)
    {
      tree->unconstraint_lk +=
	tree->data->wght[i]*(phydbl)LOG(tree->data->wght[i]);
    }
  tree->unconstraint_lk -=
    tree->data->init_len*(phydbl)LOG(tree->data->init_len);
}

/*********************************************************/

void Make_Tree_4_Lk(t_tree *tree, calign *cdata, int n_site)
{
  int i;
  
  tree->c_lnL_sorted = (phydbl *)mCalloc(tree->n_pattern, sizeof(phydbl));
  tree->cur_site_lk  = (phydbl *)mCalloc(cdata->crunch_len,sizeof(phydbl));
  tree->old_site_lk  = (phydbl *)mCalloc(cdata->crunch_len,sizeof(phydbl));
    
  tree->site_lk_cat          = (phydbl *)mCalloc(tree->mod->n_catg,sizeof(phydbl));  
  tree->log_site_lk_cat      = (phydbl **)mCalloc(tree->mod->n_catg,sizeof(phydbl *));
  For(i,tree->mod->n_catg) tree->log_site_lk_cat[i] = (phydbl *)mCalloc(cdata->crunch_len,sizeof(phydbl));
      
  tree->log_lks_aLRT = (phydbl **)mCalloc(3,sizeof(phydbl *));
  For(i,3) tree->log_lks_aLRT[i] = (phydbl *)mCalloc(tree->data->init_len,sizeof(phydbl));

  For(i,2*tree->n_otu-3)
  {
    Make_Edge_Lk(tree->t_edges[i],tree);
    Make_Edge_NNI(tree->t_edges[i]);
  }
    
  //Does nothing? Ken 4/1/2017
  For(i,2*tree->n_otu-2) Make_Node_Lk(tree->noeud[i]);

  if(tree->mod->s_opt->greedy) 
  {
    Init_P_Lk_Tips_Double(tree);
  }
  else 
  {
    Init_P_Lk_Tips_Int(tree);
  }
}

/*********************************************************/

void Init_P_Lk_Tips_Double(t_tree *tree)
{
  int curr_site,i,j,k,dim1,dim2;
  
  dim1 = tree->mod->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;

  For(curr_site,tree->data->crunch_len)
    {
      For(i,tree->n_otu)
	{
	  if (tree->io->datatype == NT)
	    Init_Tips_At_One_Site_Nucleotides_Float(tree->data->c_seq[i]->state[curr_site],
						    curr_site*dim1+0*dim2,
						    tree->noeud[i]->b[0]->p_lk_rght);

	  else if(tree->io->datatype == AA)
	    Init_Tips_At_One_Site_AA_Float(tree->data->c_seq[i]->state[curr_site],
					   curr_site*dim1+0*dim2,
					   tree->noeud[i]->b[0]->p_lk_rght);

	  else if(tree->io->datatype == GENERIC)
	    Init_Tips_At_One_Site_Generic_Float(tree->data->c_seq[i]->state+curr_site*tree->mod->state_len,
						tree->mod->ns,
						tree->mod->state_len,
						curr_site*dim1+0*dim2,
						tree->noeud[i]->b[0]->p_lk_rght);
	  else if (tree->io->datatype == CODON)//!< Added by Marcelo.
	    Init_Tips_At_One_Site_Codons_Float(tree->data->c_seq[i]->state[curr_site],
						    curr_site*dim1+0*dim2,//!< Changed by Marcelo.
						    tree->noeud[i]->b[0]->p_lk_rght,
					             tree->data->c_seq[i]->alternativeCodons[curr_site]);

	 for(j=1;j<tree->mod->n_catg;j++)
	    {
	      For(k,tree->mod->ns)
		{
		  tree->noeud[i]->b[0]->p_lk_rght[curr_site*dim1+j*dim2+k] = 
		    tree->noeud[i]->b[0]->p_lk_rght[curr_site*dim1+0*dim2+k];
		}
	    }
	}
    }

  #ifdef M4
  if(tree->mod->m4mod) M4_Init_P_Lk_Tips_Double(tree);
  #endif
}

/*********************************************************/

void Init_P_Lk_Tips_Int(t_tree *tree)
{
  int curr_site,i,dim1;

  dim1 = tree->mod->ns;

  For(curr_site,tree->data->crunch_len)
    {
      For(i,tree->n_otu)
	{
	  if(tree->io->datatype == NT)
	    Init_Tips_At_One_Site_Nucleotides_Int(tree->data->c_seq[i]->state[curr_site],
						  curr_site*dim1,
						  tree->noeud[i]->b[0]->p_lk_tip_r);

	  else if(tree->io->datatype == AA)
	    Init_Tips_At_One_Site_AA_Int(tree->data->c_seq[i]->state[curr_site],
					 curr_site*dim1,					   
					 tree->noeud[i]->b[0]->p_lk_tip_r);

	  else if(tree->io->datatype == GENERIC)
	    {
	      Init_Tips_At_One_Site_Generic_Int(tree->data->c_seq[i]->state+curr_site*tree->mod->state_len,
						tree->mod->ns,
						tree->mod->state_len,
						curr_site*dim1,
						tree->noeud[i]->b[0]->p_lk_tip_r);
	    }
	    else if(tree->io->datatype == CODON)//!< Added by Marcelo.
	    {
	      Init_Tips_At_One_Site_Codons_Int(tree->data->c_seq[i]->state[curr_site],
						  curr_site*dim1,
						  tree->noeud[i]->b[0]->p_lk_tip_r,
					          tree->data->c_seq[i]->alternativeCodons[curr_site]);
	    }
	}
    }
  #ifdef M4
  if(tree->mod->m4mod) M4_Init_P_Lk_Tips_Int(tree);
  #endif
}

/*********************************************************/

void Update_PMat_At_Given_Edge(t_edge *b_fcus, t_tree *tree)
{
  phydbl len;  
  int i,k;
  
  if(b_fcus->l < BL_MIN)      b_fcus->l = BL_MIN;
  else if(b_fcus->l > BL_MAX) b_fcus->l = BL_MAX;
   
  if(tree->io->datatype==CODON) //!< Added by Marcelo.
  {
    if(b_fcus->has_zero_br_len) 
    {
    	printf("zero branch length\nKen should have fixed this by now.\nPlease complain to him.\n");
    	exit(EXIT_FAILURE);
      For(k,tree->mod->n_w_catg) For(i,tree->mod->ns*tree->mod->ns) b_fcus->Pij_rr[i+k*tree->mod->ns*tree->mod->ns]=tree->mod->A0_part[0][i];//Ken 22/8 Need to fix
    } 
    else
    {
      if(tree->mod->n_catg>1 && tree->mod->n_w_catg==1)
      {
	For(i,tree->mod->n_catg)
	{
	  len = b_fcus->l*tree->mod->gamma_rr[i];	  
	  if(len < BL_MIN)      len = BL_MIN;
	  else if(len > BL_MAX) len = BL_MAX;
	  PMat_CODON(len,tree->mod,0,b_fcus->Pij_rr+tree->mod->ns*tree->mod->ns*i);
	}
      }
      else
      {
    	  len = b_fcus->l;
    	  if(tree->io->testcondition){len=BL_MIN;}
    	  int modeli;
    	  for(modeli=0;modeli<tree->mod->nparts;modeli++){
    		  PMat_CODON_part(len,tree->mod,0,tree->mod->qmat_part[modeli],b_fcus->bPmat_part[modeli],modeli);
    	  }
      }
    }
  }
  else
  {
    
    len = -1.0;
    
    For(i,tree->mod->n_catg)
    {
      if(b_fcus->has_zero_br_len) len = -1.0;
      else
      {
	len = b_fcus->l*tree->mod->gamma_rr[i];	  
	if(len < BL_MIN)      len = BL_MIN;
	else if(len > BL_MAX) len = BL_MAX;
      }
      PMat(len,tree->mod,tree->mod->ns*tree->mod->ns*i,b_fcus->Pij_rr);
    }
  }
}

/*********************************************************/

void Update_P_Lk_Along_A_Path(t_node **path, int path_length, t_tree *tree)
{
  int i,j;

  For(i,path_length-1)
    {
      For(j,3)
	if(path[i]->v[j] == path[i+1])
	  {
	    if(path[i] == path[i]->b[j]->left)
	      {
		Update_P_Lk(tree,path[i]->b[j],path[i]->b[j]->left);		    
	      }

	    else if(path[i] == path[i]->b[j]->rght)
	      {
		Update_P_Lk(tree,path[i]->b[j],path[i]->b[j]->rght);
	      }
	    else
	      {
		PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
		Exit("");
	      }
	    break;
	  }      
#ifdef DEBUG
      if(j == 3)
	{
	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  Exit("");
	}
#endif
    }
}

/*********************************************************/

phydbl Lk_Dist(phydbl *F, phydbl dist, model *mod)
{
  int i,j,k;
  phydbl lnL,len;
  int dim1,dim2;

  For(k,mod->n_catg)
    {
      len = dist*mod->gamma_rr[k];
      if(len < BL_MIN)      len = BL_MIN;
      else if(len > BL_MAX) len = BL_MAX;
      PMat(len,mod,mod->ns*mod->ns*k,mod->Pij_rr);
    }
  
  dim1 = mod->ns*mod->ns;
  dim2 = mod->ns;
  lnL = .0;
  For(i,mod->ns)
    {
      For(j,mod->ns)
	{
	  For(k,mod->n_catg)
	    {
 	      lnL += F[dim1*k+dim2*i+j] * LOG(mod->Pij_rr[dim1*k+dim2*i+j]);
	    }
	}
    }

  return lnL;
}

/*********************************************************/
//changed by Ken
phydbl Update_Lk_At_Given_Edge(t_edge *b_fcus, t_tree *tree)
{
  if(!b_fcus->left->tax) Update_P_Lk(tree,b_fcus,b_fcus->left);
  if(!b_fcus->rght->tax) Update_P_Lk(tree,b_fcus,b_fcus->rght);
  tree->c_lnL = Lk_At_Given_Edge(b_fcus,tree);
  return tree->c_lnL;
}

/*********************************************************/

/*       root
           \
           /
          a
	  |
	  |
	  d
	 / \
        /   \
       w     x    	

       d->t has changed and we need to compute
       the likelihood.
       (1) update the three branch lengths l(ad), l(dw) and l(dx)
       (2) update the change proba matrices along these branches
       (3) update the likelihood of subtree (w,x) (WARNING: (a,x) and (a,w) are not updated)
*/
phydbl Lk_Triplet(t_node *a, t_node *d, t_tree *tree)
{
  int i;
  phydbl max_height;
  phydbl up_bound, low_bound;

  if(d->tax)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  up_bound = low_bound = -1.0;
  max_height = -1.0;
  For(i,3)
    {
      if((d->v[i] != a) && (d->b[i] != tree->e_root))
	{
	  if(tree->rates->nd_t[d->v[i]->num] > max_height)
	    {
	      max_height = tree->rates->nd_t[d->v[i]->num];
	    }
	}
      else
	{
	  up_bound = 
	    (a == tree->n_root)?
	    (tree->rates->nd_t[a->num]):
	    (tree->rates->nd_t[d->v[i]->num]);
	}
    }

  low_bound = max_height;

  if(up_bound < low_bound - 1.E-10)
    {
      PhyML_Printf("\n. a->num=%d d->num=%d",a->num,d->num);
      PhyML_Printf("\n. up_bound = %f, low_bound = %f",up_bound,low_bound);
      Warn_And_Exit("\n");
    }

  if(tree->rates->nd_t[d->num] < low_bound) tree->rates->nd_t[d->num] = low_bound;
  else if(tree->rates->nd_t[d->num] > up_bound) tree->rates->nd_t[d->num] = up_bound;

  /* Step (1) */
  For(i,3)
    {
      if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	{
	  d->b[i]->l = 
	    (tree->rates->nd_t[d->num] - tree->rates->nd_t[d->v[i]->num]) * 
	    tree->rates->clock_r * 
	    tree->rates->br_r[d->b[i]->num];
	}
      else
	{
	  if(a == tree->n_root)
	    {
	      d->b[i]->l = 
		(tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->n_root->v[0]->num] + 
		 tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->n_root->v[1]->num]) * tree->rates->clock_r;
	    }
	  else
	    {
	      d->b[i]->l = (tree->rates->nd_t[a->num] - tree->rates->nd_t[d->num]) * tree->rates->clock_r * tree->rates->br_r[d->b[i]->num];	    
	    }
	}
    }
  
  /* Step (2) */
  For(i,3) Update_PMat_At_Given_Edge(d->b[i],tree);
  
  For(i,3) 
    if((d->v[i] == a) || (d->b[i] == tree->e_root))
      {
	Update_P_Lk(tree,d->b[i],d); 
	Lk_At_Given_Edge(d->b[i],tree);
	break;
      }

  return tree->c_lnL;
}

/*********************************************************/

void Print_Lk_Given_Edge_Recurr(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  PhyML_Printf("\n___ Edge %3d (left=%3d rght=%3d) lnL=%f",
	 b->num,
	 b->left->num,
	 b->rght->num,
	 Lk_At_Given_Edge(b,tree));

  if(d->tax) return;
  else
    {
      int i;
      For(i,3)
	if(d->v[i] != a)
	  Print_Lk_Given_Edge_Recurr(d,d->v[i],d->b[i],tree);
    }
}

/*********************************************************/
void Init_Tips_At_One_Site_Codons_Float(char state, int pos, phydbl *p_lk, char * alternatives) //!< Added by Marcelo.
{
  int i,j,numSenseCodons;
  numSenseCodons=Genetic_code_ns();
  
  For(i,numSenseCodons) p_lk[pos+i]=0.0;
  
  if((state>=(char)0) && (state<(char)64)) p_lk[pos+indexSenseCodons[(int)state]]=1.0; 
  else
    if(state==(char)88)
    {
      j=-1;
      while(alternatives[++j]<(char)64) p_lk[pos+indexSenseCodons[(int)alternatives[j]]]=1.0;
    }else 
    {
      PhyML_Printf("\n. Unknown character state : %c\n",state);
      Exit("\n. Init failed (check the data type)\n");
    }
}
/*********************************************************/

void Init_Tips_At_One_Site_Codons_Int(char state, int pos, short int *p_pars, char * alternatives) //!< Added by Marcelo.
{
  int i,j,numSenseCodons;
  numSenseCodons=Genetic_code_ns();
  
  For(i,numSenseCodons) p_pars[pos+i]=0;
  
  if((state>=(char)0) && (state<(char)64)) p_pars[pos+indexSenseCodons[(int)state]]=1; 
  else
    if(state==(char)88)
    {
      j=-1;
      while(alternatives[++j]<(char)64)
      {
	p_pars[pos+indexSenseCodons[(int)alternatives[j]]]=1;
      }    
    }else 
    {
      PhyML_Printf("\n. Unknown character state : %c\n",state);
      Exit("\n. Init failed (check the data type)\n");
    }
}
/*********************************************************/
matrix *ML_CODONDist_Pairwise(calign *data, option *io) //!<Added by Marcelo.
{
  int j,k,tmpkappa;
  model *mod_tmp;
  optimiz *s_opt;
  matrix *mat;
  calign *twodata,*tmpdata;
  phydbl lnL, *expt, *uexpt, brLen, mr;
  time_t tbegin, tend;
 

  mod_tmp = Make_Model_Basic();
  s_opt = Make_Optimiz();
  Set_Defaults_Model(mod_tmp);
  mod_tmp->genetic_code = io->genCode;
  Genetic_code_index_stopCodons(mod_tmp->genetic_code);
  Set_Defaults_Optimiz(s_opt);
  mod_tmp->io = io;
  mod_tmp->s_opt = s_opt;
  mod_tmp->ns = io->mod->ns;
  mod_tmp->whichmodel = GY;
  mod_tmp->whichrealmodel = GY;
  mod_tmp->initqrates = io->init_DistanceTreeCD;
  mod_tmp->freq_model = FMODEL;
  mod_tmp->omegaSiteVar = NOOMEGA;
  tmpkappa = io->kappaECM;
  io->kappaECM = kap1;
  mod_tmp->nomega_part=io->mod->nomega_part;
  mod_tmp->nparts=io->mod->nparts;
  Make_Model_Complete(mod_tmp);
  mod_tmp->s_opt->opt_kappa = NO;

  expt            = (phydbl *)mCalloc(io->mod->ns,sizeof(phydbl));
  uexpt           = (phydbl *)mCalloc(io->mod->ns*io->mod->ns,sizeof(phydbl));
  
  tmpdata         = (calign *)mCalloc(1,sizeof(calign));
  tmpdata->c_seq  = (align **)mCalloc(2,sizeof(align *));
  tmpdata->b_frq  = (phydbl *)mCalloc(io->mod->ns,sizeof(phydbl));
  tmpdata->ambigu = (short int *)mCalloc(data->crunch_len,sizeof(short int));
  
  tmpdata->n_otu  = 2;
  tmpdata->crunch_len = data->crunch_len;
  tmpdata->init_len   = data->init_len;
  
  mat = Make_Mat(data->n_otu);
  Init_Mat(mat,data);
  
  //Added by Ken
  //Make Bmat array equivalent to GY94
  if(io->mod->whichrealmodel == HLP17){
	  int combinations = 3721;
	  phydbl hot[combinations];
	  phydbl *tfr = ecmK07freq;
	  int i1, i2;
	  for(i1=0;i1<61;i1++){
		  for(i2=0;i2<61;i2++){
			  hot[i1*61+i2]=0;
		  }
	  }
	  mod_tmp->Bmat = hot;

	  mod_tmp->nmotifs = 1;
	  mod_tmp->nhotness = 1;
	  mod_tmp->hotspotcmps = malloc(sizeof(int *) * mod_tmp->nmotifs);
	  mod_tmp->hotness = malloc(sizeof(double ) * mod_tmp->nhotness);
	  mod_tmp->motif_hotness = malloc(sizeof(int ) * mod_tmp->nmotifs);
	  mod_tmp->hotness[0] = 0.0;
	  mod_tmp->motif_hotness[0]=0;
	  mod_tmp->hotspotcmps[0] = io->mod->hotspotcmps[0];
  }

  //Add in stuff for partitioned model
  int modeli;
  for(modeli=0;modeli<io->mod->nparts;modeli++){
	  mod_tmp->nparts=1;
	  mod_tmp->nomega_part=mod_tmp->nparts;
	  mod_tmp->qmat_part = (phydbl**)mCalloc(1,sizeof(phydbl*));
	  mod_tmp->partNames = (char**)mCalloc(1,sizeof(char*));
	  mod_tmp->omega_part = (phydbl*)mCalloc(1,sizeof(phydbl));
	  mod_tmp->qmat_part[0]=(phydbl *)mCalloc(io->mod->n_w_catg*io->mod->ns*io->mod->ns,sizeof(phydbl));
	  mod_tmp->partNames[0]="SINGLE";
	  mod_tmp->omega_part[0]=0.4;
  }

  mr = Update_Qmat_Codons(mod_tmp, 0,0);

  printf("\n. FYI: Using single partitioned GY94 for initial tree.\n\n");
  EigenQREV(mod_tmp->qmat_part[0], mod_tmp->pi, io->mod->ns, mod_tmp->eigen->e_val, mod_tmp->eigen->r_e_vect, mod_tmp->eigen->l_e_vect, mod_tmp->eigen->space);

  For(j, io->mod->ns) mod_tmp->eigen->e_val[j] /= mr;
  
  brLen=1.0;
  
  #if defined OMP || defined BLAS_OMP

  tbegin= omp_get_wtime();
  
  #else

  time(&tbegin);
 
  #endif
  int progress = 0; //!Added by Marcelo.
  double loops = ((data->n_otu)*(data->n_otu)-(data->n_otu))/2.00;//!Added by Marcelo.
  For( j, data->n_otu-1 ) {
    tmpdata->c_seq[0]       = data->c_seq[j];
    tmpdata->wght           = data->wght;
    
    for(k=j+1;k<data->n_otu;k++) {
      tmpdata->c_seq[1]       = data->c_seq[k];
      twodata = Compact_Cdata(tmpdata,io);
                                                            
      lnL=Br_Len_Brent_Codon_Pairwise(BL_MIN, brLen, BL_MAX, io->mod->s_opt->min_diff_lk_local, &brLen, mod_tmp->Pij_rr, mod_tmp->pi, mod_tmp->eigen, twodata, io->mod->ns, io->mod->s_opt->brent_it_max, io->mod->s_opt->quickdirty, uexpt, expt);
            
      mat->dist[j][k] = brLen;
      mat->dist[k][j] = mat->dist[j][k];
      Free_Cseq(twodata);
      progress++;
      if(!io->quiet) PhyML_Printf("\r. Computing pairwise distances...%3.2f%c concluded.", progress*100.00/loops, '%');
    }
  }
  if(!io->quiet) PhyML_Printf("\n");
  #if defined OMP || defined BLAS_OMP

  tend= omp_get_wtime();
  
  #else

  time(&tend);
 
  #endif
  
  if((tend-tbegin)>=1) {
    PhyML_Printf("\n. Calculation of pairwise distances completed in %ld sec.\n",(long int)(tend-tbegin));
  } else {
    PhyML_Printf("\n. Calculation of pairwise distances completed in less than 1 sec.\n");
  }
  
  tmpdata->c_seq[0] = NULL;
  tmpdata->c_seq[1] = NULL;
  tmpdata->wght     = NULL;
  
  free(tmpdata->ambigu);
  free(tmpdata->b_frq);
  free(tmpdata->c_seq);
  free(tmpdata);
  free(expt); 
  free(uexpt);
  Free_Optimiz(s_opt);
  Free_Model(mod_tmp);
  
  io->kappaECM = tmpkappa;

  return mat;
}

/*********************************************************/
phydbl LK_Codon_Pairwise(calign *data, phydbl *Pij, phydbl *pi, int ns, phydbl len, eigen *eigenStruct, phydbl *uexpt, phydbl *expt) //!< Added by Marcelo.
{
  phydbl cn_lk;
  // phydbl p_sum0;
  int i;
  // int j,k,m,n;
  char *seq0, *seq1;
  // char state0, state1;
  
  PMat_CODON_Pairwise(len, Pij, eigenStruct->r_e_vect, eigenStruct->l_e_vect, eigenStruct->e_val, ns, uexpt, expt);

  cn_lk=0.0;
  seq0=data->c_seq[0]->state;
  seq1=data->c_seq[1]->state;
  
  #if defined OMP || defined BLAS_OMP
  
  #pragma omp parallel for reduction(+:cn_lk)
  
  #endif
  
  For(i,data->crunch_len)
  {
    if((seq0[i]!=(char)88)&&(seq1[i]!=(char)88))
    {
      cn_lk+=data->wght[i]*LOG(pi[indexSenseCodons[(int)seq0[i]]]*Pij[indexSenseCodons[(int)seq0[i]]*ns+indexSenseCodons[(int)seq1[i]]]);
    }
  } 
  

  return cn_lk;
}


/*********************************************************/
phydbl LK_BFGS_from_CODEML(t_tree* tree, phydbl *x, int n)
{
  int  i, numParams=0;
  
  if (tree->io->datatype==CODON)
  {
    if(tree->mod->s_opt->opt_kappa)
    {
      if(tree->mod->whichmodel!=GYECMK07WK && tree->mod->whichmodel!=GYECMK07WKF && 
	tree->mod->whichmodel!=GYECMS05WK && tree->mod->whichmodel!=GYECMS05WKF
	&& tree->mod->whichmodel!=MGECMK07WK && tree->mod->whichmodel!=MGECMK07WKF && 
	tree->mod->whichmodel!=MGECMS05WK && tree->mod->whichmodel!=MGECMS05WKF
	&& tree->mod->whichmodel!=YAPECMK07WK && tree->mod->whichmodel!=YAPECMK07WKF && 
	tree->mod->whichmodel!=YAPECMS05WK && tree->mod->whichmodel!=YAPECMS05WKF &&
	tree->mod->whichmodel!=GYECMUSRWK && tree->mod->whichmodel!=GYECMUSRWKF &&
	   tree->mod->whichmodel!=MGECMUSRWKF && tree->mod->whichmodel!=YAPECMUSRWKF &&
	   tree->mod->whichmodel!=MGECMUSRWK &&tree->mod->whichmodel!=YAPECMUSRWK )
      {
	tree->mod->kappa = x[numParams++];
      }
      else
      {
	switch(tree->io->kappaECM)
	{
	  case kap6:
	  case kap2:
	  case kap3:          tree->mod->pkappa[0] = x[numParams++];break;
	  case kap4:          tree->mod->pkappa[0] = x[numParams++]; tree->mod->pkappa[1] = x[numParams++]; break;
	  case kap5: For(i,tree->mod->nkappa-1) tree->mod->unspkappa[i] = x[numParams++]; break;
	  default: break;
	}
      }
    }
    
    if(tree->mod->s_opt->opt_omega)
    {
      if(tree->mod->omegaSiteVar==DM0){
    	  int omegai; //added by Ken 18/8/2016
    	  for(omegai=0;omegai<tree->mod->nomega_part;omegai++){
    		  tree->mod->omega_part[omegai] = x[numParams++];
    	  }
      }
      else if(tree->mod->omegaSiteVar==DMODELK) 
      {	
	For(i,tree->mod->n_w_catg) tree->mod->omegas[i] = x[numParams++];
	For(i,tree->mod->n_w_catg-1) tree->mod->prob_omegas_uns[i] = x[numParams++];
      }
      else if(tree->mod->omegaSiteVar==DGAMMAK)
      {
	tree->mod->alpha                = x[numParams++];
	tree->mod->beta                 = x[numParams++];
      }
    }

    if(tree->mod->opthotness){ //added by Kenneth Hoehn 3/6/2016
    		int c;
    		for(c=0;c<tree->mod->nhotness;c++){
    			if(tree->mod->hoptindex[c] == 1){
    				tree->mod->hotness[c] = x[numParams++];
    			}
    		}
     }

    if(tree->mod->s_opt->opt_state_freq)
    {
      switch(tree->mod->freq_model)
      {
	case F1XSENSECODONS: For(i,tree->mod->num_base_freq-1) tree->mod->pi_unscaled[i] = x[numParams++]; break;
	case F1X4: For(i,tree->mod->num_base_freq-1) tree->mod->uns_base_freq[i] = x[numParams++]; break;
	case F3X4:
	case CF3X4: 
	{
	  for(i=0;i<3;i++) tree->mod->uns_base_freq[i] = x[numParams++]; 
	  for(i=4;i<7;i++) tree->mod->uns_base_freq[i] = x[numParams++]; 
	  for(i=8;i<11;i++)tree->mod->uns_base_freq[i] = x[numParams++]; 
	  break;
	}
	default:
	  break;
      }
    }

    
    if(tree->mod->s_opt->opt_pinvar) tree->mod->pinvar = x[numParams++];
    
    if(tree->mod->n_catg>1 && tree->mod->n_w_catg==1 && tree->mod->s_opt->opt_alphaCD) tree->mod->alpha = x[numParams++];
    
    if(tree->mod->pcaModel==1) For(i,tree->mod->npcs) tree->mod->pcsC[i] = x[numParams++];
    
    if(n!=numParams) Warn_And_Exit("Number of parameters optimized does not match.");
    
  } 
  else if(tree->io->datatype==AA)
  {
    
    if(tree->mod->s_opt->opt_state_freq) For(i,tree->mod->ns-1) tree->mod->pi_unscaled[i] = x[numParams++]; 
    
    if(n!=numParams) Warn_And_Exit("Number of parameters optimized does not match.");
  }

  return -Lk(tree);
}

/*********************************************************/
