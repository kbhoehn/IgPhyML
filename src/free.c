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

#include "free.h"


/*********************************************************/
#include "utilities.h"

void Free_All_Nodes_Light(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-1) //hereeeee
    Free_Node(tree->noeud[i]);
  free(tree->noeud);
}

/*********************************************************/

void Free_All_Edges_Light(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-2) 
    if(tree->t_edges[i])
      Free_Edge(tree->t_edges[i]);
  free(tree->t_edges);
}


/*********************************************************/

void Free_Mat(matrix *mat)
{
  int i;

  For(i,mat->n_otu)
    {
      free(mat->P[i]);
      free(mat->Q[i]);
      free(mat->dist[i]);
      free(mat->name[i]);
    }

  free(mat->P);
  free(mat->Q);
  free(mat->dist);
  free(mat->name);
  free(mat->tip_node);
      
  free(mat->on_off);
  free(mat);
}

/*********************************************************/

void Free_Partial_Lk(phydbl *p_lk, int len, int n_catg)
{
  free(p_lk);

/*   int i,j; */
/*   For(i,len) */
/*     { */
/*       For(j,n_catg) free((*p_lk)[i][j]); */
/*       free((*p_lk)[i]); */
/*     } */
/*   free((*p_lk)); */
/*   (*p_lk) = NULL; */
}

/*********************************************************/

void Free_Tree(t_tree *tree)
{
  // int i,j,k;

  free(tree->t_dir);

  Free_Bip(tree);

  free(tree->curr_path);

  Free_All_Edges_Light(tree);
  Free_All_Nodes_Light(tree);

  free(tree);
}

/*********************************************************/

void Free_Bip(t_tree *tree)
{
  int i,j;

  if(tree->has_bip)
    {
      For(i,2*tree->n_otu-2)
	{
	  free(tree->noeud[i]->bip_size);
	  For(j,3) free(tree->noeud[i]->bip_node[j]);
	  free(tree->noeud[i]->bip_node);
	}
    }
  tree->has_bip = NO;
}

/*********************************************************/

void Free_Edge_Labels(t_edge *b)
{
  int i;
  For(i,b->n_labels+b->n_labels%BLOCK_LABELS) free(b->labels[i]);
  free(b->labels);
  b->labels = NULL;
}

/*********************************************************/

void Free_Edge(t_edge *b)
{
  Free_Edge_Labels(b);
  free(b);
}

/*********************************************************/

void Free_Node(t_node *n)
{
  //int i;

  free(n->b);
  free(n->v);
  free(n->l);
  free(n->score);
  
  free(n->ori_name);//!< Added by Marcelo.
  free(n->name);//!< Added by Marcelo.
  
  free(n);
}

/*********************************************************/

void Free_Cseq(calign *data)
{
  //int i,j;
  
  free(data->invar);
  free(data->wght);
  free(data->ambigu);
  free(data->b_frq);
  free(data->sitepatt);
  Free_Seq(data->c_seq, data->n_otu);
  free(data);
}

/*********************************************************/

void Free_Seq(align **d, int n_otu)//, int datatype)
{
  int i,j;
  For(i,n_otu)
  {
    free(d[i]->name);
    free(d[i]->state);
    if(d[i]->is_ambigu) free(d[i]->is_ambigu);
    if(d[i]->ntStates) free(d[i]->ntStates);
    if(d[i]->alternativeCodons)
    {
      For(j,d[i]->len) if(d[i]->alternativeCodons[j]) free(d[i]->alternativeCodons[j]);
      free(d[i]->alternativeCodons);
    }
    free(d[i]);
  }
  free(d);
}

/*********************************************************/

void Free_All(align **d, calign *cdata, t_tree *tree)
{
  Free_Cseq(cdata);
  Free_Seq(d,tree->n_otu);
  Free_Tree(tree);
}      

/*********************************************************/
void Free_SubTree(t_edge *b_fcus, t_node *a, t_node *d, t_tree *tree)
{
  int i;

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    {
	      Free_SubTree(d->b[i],d,d->v[i],tree);
	      Free_Edge(d->b[i]);
	      Free_Node(d->v[i]);
	    }
	}
    }
}

/*********************************************************/
void Free_Tree_Ins_Tar(t_tree *tree)
{
  return;
}

/*********************************************************/

void Free_Tree_Pars(t_tree *tree)
{
  int i;
  
  free(tree->step_mat);
  free(tree->site_pars);
  For(i,2*tree->n_otu-3)
    Free_Edge_Pars(tree->t_edges[i],tree);
}

/*********************************************************/

void Free_Edge_Pars(t_edge *b, t_tree *tree)
{
/*   int i; */

  free(b->pars_l);
  free(b->pars_r);
  
/*   For(i,tree->data->crunch_len)  */
/*     { */
/*       free(b->p_pars_l[i]); */
/*       free(b->p_pars_r[i]); */
/*     } */
  
  free(b->ui_l);
  free(b->ui_r);
  free(b->p_pars_l);
  free(b->p_pars_r);
}

/*********************************************************/

void Free_Tree_Lk(t_tree *tree)
{
  int i;
  t_edge *b;
  t_node *n;

  b = NULL;
  n = NULL;

  For(i,3) free(tree->log_lks_aLRT[i]);
  free(tree->log_lks_aLRT);

  free(tree->c_lnL_sorted);
  free(tree->cur_site_lk);
  free(tree->old_site_lk);
  free(tree->site_lk_cat);

  For(i,tree->mod->n_catg) free(tree->log_site_lk_cat[i]);
  free(tree->log_site_lk_cat);
				
  For(i,2*tree->n_otu-3)
    {
      b = tree->t_edges[i];      
      Free_Edge_Lk(tree,b);
    }
}


/*********************************************************/


void Free_Node_Lk(t_node *n)
{
/*   free(n->n_ex_nodes); */
}

/*********************************************************/

void Free_Edge_Lk(t_tree *tree, t_edge *b)
{
  
  free(b->nni);

  free(b->div_post_pred_left);
  free(b->div_post_pred_rght);

  if(b->p_lk_left)
    {
      free(b->p_lk_left);
      if(b->sum_scale_left) free(b->sum_scale_left);
    }

  if(b->p_lk_tip_l) free(b->p_lk_tip_l);


  if(b->p_lk_rght)
    {
      free(b->p_lk_rght);
      if(b->sum_scale_rght) free(b->sum_scale_rght);
    }
  
  if(b->p_lk_tip_r) free(b->p_lk_tip_r);

  free(b->sum_scale_left_cat);
  free(b->sum_scale_rght_cat);

  free(b->Pij_rr);
}

/*********************************************************/

void Free_Model_Complete(model *mod)
{
  Free_Eigen(mod->eigen); 
  free(mod->gamma_r_proba);
  free(mod->gamma_rr);
  free(mod->qmat);
  free(mod->Pij_rr);
  free(mod->pi_unscaled);                           
  free(mod->pi);                                   
  if(mod->n_rr_branch)
  {
    free(mod->rr_branch);                           
    free(mod->p_rr_branch);                         
  }
  if(mod->io->datatype==CODON) //!< Added by Marcelo. 
  {
    free(mod->prob_omegas_uns);                   
    free(mod->base_freq);                           
    free(mod->uns_base_freq);                        
    
    if(mod->io->heuristicExpm)
    {
      free(mod->A2_part[0]);
    }
    
    if(mod->io->expm==SSPADE)
    {
    int modeli; //Added by Ken 22/8
    for(modeli=0;modeli<mod->nparts;modeli++){
      free(mod->U_part[modeli]);
      free(mod->V_part[modeli]);
      free(mod->ipiv_part[modeli]);
      free(mod->A2_part[modeli]);
      free(mod->A4_part[modeli]);
      free(mod->A6_part[modeli]);
      free(mod->A8_part[modeli]);
      free(mod->matAux_part[modeli]);
      free(mod->Apowers_part[modeli]);
      free(mod->A0_part[modeli]);
      free(mod->qmat_buff_part[modeli]);
    }
    }
    
    free(mod->qmatScaled);
    free(mod->mr_w);
    
//    if(mod->whichmodel<-27 && mod->whichmodel>-38) //!< Added by Marcelo. 
//    {
//      int i;
//      For (i,64) free(mod->userRates[i]);
//      free(mod->userRates);
//      free(mod->userfreq);
//    }
  }
  
 
}

/*********************************************************/

void Free_Model_Basic(model *mod)
{
  free(mod->modelname);
  free(mod->custom_mod_string);
  free(mod->user_b_freq);
  
  free(mod->pkappa); //!< Added by Marcelo.
  free(mod->unspkappa); //!< Added by Marcelo.
  
  if(mod->io->datatype==CODON) //!< Added by Marcelo. 
  {
    free(mod->omegas);
    free(mod->prob_omegas);
  }
  free(mod); //!< Added by Marcelo.
}

/*********************************************************/

void Free_Custom_Model(model *mod)
{
  if(mod->rr)
    {
      free(mod->rr_num);
      free(mod->rr);
      free(mod->rr_val);
      free(mod->n_rr_per_cat);
    }
}

/*********************************************************/
void Free_Model(model *mod)
{

  Free_Model_Complete(mod);
  
  if(mod->whichmodel==GTR) Free_Custom_Model(mod);
  
  #ifdef M4
  M4_Free_M4_Model(mod->m4mod);
  #endif 
  Free_Model_Basic(mod);
}

/*********************************************************/

void Free_Input(option *io)
{
  int i;
  free(io->in_align_file);
  free(io->in_tree_file);
  free(io->out_tree_file);
  free(io->out_trees_file);
  free(io->out_boot_tree_file);
  free(io->out_boot_stats_file);
  free(io->out_stats_file);
  free(io->out_lk_file); 
  free(io->out_ps_file);
  //free(io->out_trace_file);
  free(io->out_trace_stats_file);
  free(io->out_trace_tree_file);
  free(io->nt_or_cd);
  free(io->run_id_string);
  free(io->clade_list_file);
  For(i,T_MAX_ALPHABET) free(io->alphabet[i]);
  free(io->alphabet);
  if(io->short_tax_names)
    {
      For(i,io->size_tax_names) 
	{
	  free(io->short_tax_names[i]);
	  free(io->long_tax_names[i]);
	}
      free(io->long_tax_names);
      free(io->short_tax_names);
    }
  Free_Tree_List(io->treelist);
  
  if(io->lon) free(io->lon);
  if(io->lat) free(io->lat);
  if(io->z_scores) free(io->z_scores);//!< Added by Marcelo.
  free(io->structTs_and_Tv);          //!< Added by Marcelo.
  free(io);
}

/*********************************************************/

void Free_Tree_List(t_treelist *list)
{
  free(list->tree);
  free(list);
}

/*********************************************************/

void Free_St(supert_tree *st)
{
  int i;

  For(i,2*st->tree->n_otu-3) 
    free(st->tree->t_edges[i]->nni);

  For(i,st->n_part) free(st->match_st_node_in_gt[i]);

  free(st->match_st_node_in_gt);

  Free_Tree(st->tree);
  
  free(st);
}

/*********************************************************/

void Free_Eigen(eigen *eigen_struct)
{
  free(eigen_struct->space_int);
  free(eigen_struct->space);
  free(eigen_struct->e_val);
  free(eigen_struct->e_val_im);
  free(eigen_struct->r_e_vect);
  free(eigen_struct->r_e_vect_im);
  free(eigen_struct->l_e_vect);
  free(eigen_struct->q);
  free(eigen_struct);
}

/*********************************************************/

void Free_One_Spr(spr *this_spr)
{
  free(this_spr->path);
  free(this_spr);
}

/*********************************************************/

void Free_Spr_List(t_tree *tree)
{
  int i;

  For(i,tree->size_spr_list+1) Free_One_Spr(tree->spr_list[i]);
  free(tree->spr_list);

}

/*********************************************************/


void Free_Triplet(triplet *t)
{
  int i,j,k;

  free(t->F_bc);
  free(t->F_cd);
  free(t->F_bd);
  free(t->pi_bc);
  free(t->pi_cd);
  free(t->pi_bd);
  
  For(k,t->mod->n_catg) 
    {
      For(i,t->size) 
	{
	  For(j,t->size) free(t->core[k][i][j]);  
	  free(t->core[k][i]);
	}
      free(t->core[k]);	  
    }
  free(t->core);

  For(i,t->size) 
    {
      For(j,t->size) free(t->p_one_site[i][j]);  
      free(t->p_one_site[i]);
    }
  free(t->p_one_site);

  For(i,t->size) 
    {
      For(j,t->size) free(t->sum_p_one_site[i][j]);  
      free(t->sum_p_one_site[i]);
    }
  free(t->sum_p_one_site);

  Free_Eigen(t->eigen_struct);
  
  free(t);
}

/*********************************************************/

void Free_Actual_CSeq(calign *data)
{
  int i;
  For(i,data->n_otu)
    {
      free(data->c_seq[i]->state);
      data->c_seq[i]->state = NULL;
    }
}

/*********************************************************/

void Free_Prefix_Tree(pnode *n, int size)
{
  int i;
  
  For(i,size)
    {
      if(n->next[i])
	{
	  Free_Prefix_Tree(n->next[i],size);
	}
    }
  Free_Pnode(n);
}

/*********************************************************/

void Free_Pnode(pnode *n)
{
  free(n->next);
  free(n);
}

/*********************************************************/

void Free_Optimiz(optimiz *s_opt)
{
  free(s_opt);
}

/*********************************************************/
void Free_Nexus(option *io)
{
  int i,j;
  
  For(i,N_MAX_NEX_COM)
    {
      For(j,io->nex_com_list[i]->nparm) Free_Nexus_Parm(io->nex_com_list[i]->parm[j]);
      free(io->nex_com_list[i]->parm);
      free(io->nex_com_list[i]->name);
      free(io->nex_com_list[i]);      
    }
  free(io->nex_com_list);
}

/*********************************************************/

void Free_Nexus_Com(nexcom **com)
{
  int i;

  For(i,N_MAX_NEX_COM)
    {
      free(com[i]->parm);
      free(com[i]->name);
      free(com[i]);
    }
  free(com);
}

/*********************************************************/

void Free_Nexus_Parm(nexparm *parm)
{
  free(parm->value);
  free(parm->name);
  free(parm);
}

/*********************************************************/

void Free_OutList(token *root) {
    token *curr = root;
    token *tmp;
    while(curr->next != NULL) {
        tmp = curr->next;
        Free_OutToken(curr);
        curr = tmp;
    }
    Free_OutToken(curr);
}

void Free_OutToken(token *t) {
    free(t->full);
    free(t->tag);
    free(t->val);
    free(t);
}

/*********************************************************/

void Free_DarwinList(dtoken *root) {
    dtoken *curr = root;
    dtoken *tmp;
    while(curr->next != NULL) {
        tmp = curr->next;
        Free_DarwinToken(curr);
        curr = tmp;
    }
    Free_DarwinToken(curr);
}

void Free_DarwinToken(dtoken *t) {
    free(t->name);
    free(t->val);
    free(t);
}

/*********************************************************/
/*********************************************************/
