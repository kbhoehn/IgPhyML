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

/*
** spr.c: Routines for performing SPR moves on the tree.
**
** Wim Hordijk   Last modified: 28 August 2006
** Stephane Guindon 2007
*/

#include "spr.h"



/*********************************************************/
/* Below are my functions for SPR search (Stephane Guindon, 2007) */


void Randomize_Spr_List(t_tree *tree)
{
  int i,j;
  spr *buff;

  For(i,tree->size_spr_list)
    {
      j = (int)FLOOR(rand()/(RAND_MAX+1.)*tree->size_spr_list);
      buff              = tree->spr_list[i];
      tree->spr_list[i] = tree->spr_list[j];
      tree->spr_list[j] = buff;
    }
}

/*********************************************************/

int Spr(phydbl init_lnL, t_tree *tree)
{
  int br;
  int pars_diff, max_pars_diff, new_pars, old_pars;
  t_edge *b;

  tree->both_sides = 1;
  pars_diff        = -1;
  max_pars_diff    = -1;
  new_pars         = -1;
  old_pars         = -1;

  Reset_Spr_List(tree);

  For(br,2*tree->n_otu-3){ //loop over all edges in tree
      b = tree->t_edges[br];
      if(tree->mod->whichrealmodel==HLP17){
    	    Update_Ancestors_Edge(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree->noeud[tree->mod->startnode]->b[0],tree); //added by Ken 7/11
    	    Lk(tree);
    		Get_UPP(tree->noeud[tree->mod->startnode], tree->noeud[tree->mod->startnode]->v[0], tree);
    		old_pars = tree->c_pars;
    		Spr_Subtree(b,b->anc_node,tree); //don't consider subtrees that prune out the root
    		new_pars = tree->c_pars;
    		pars_diff =  new_pars - old_pars;
    		if(pars_diff > max_pars_diff) max_pars_diff = pars_diff;
      }else{
    	  	  old_pars = tree->c_pars;
    	      Spr_Subtree(b,b->left,tree);
    	      new_pars = tree->c_pars;

    	      pars_diff =  new_pars - old_pars;
    	      if(pars_diff > max_pars_diff) max_pars_diff = pars_diff;

    	      old_pars = tree->c_pars;
    	      Spr_Subtree(b,b->rght,tree);
    	      new_pars = tree->c_pars;

    	      pars_diff = new_pars - old_pars;
    	      if(pars_diff > max_pars_diff) max_pars_diff = pars_diff;
      }
    }

  return 1;
}

/*********************************************************/
//b = target edge, link = ancestral node
void Spr_Subtree(t_edge *b, t_node *link, t_tree *tree) 
{
  int i;
  int n_moves_pars, n_moves, curr_pars, min_pars, best_move;
  spr *best_pars_move;
  t_edge *target, *residual;
  int ret_val;

  best_move     = -1;
  tree->n_moves = 0;
  curr_pars     = tree->c_pars;
  ret_val       = 0;
  if((link != b->left) && (link != b->rght)){
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");			        
    }else{
      if(!link->tax){
    	  Test_All_Spr_Targets(b,link,tree); //check through all potential SPR targets
      }
      if(tree->n_moves){ //if at least one potential move was found
    	  n_moves_pars = 0;
    	  n_moves      = 0;
    	  For(i,tree->n_moves){
    		  if(curr_pars - tree->spr_list[i]->pars >= -tree->mod->s_opt->pars_thresh){
    			  n_moves_pars++;
    		  }
    	  }
    	  n_moves_pars = MAX(n_moves_pars,1);
	  
    	  //set number of moves to at least 15
    	  if(tree->mod->s_opt->spr_lnL) n_moves = 15;
    	  else                          n_moves = n_moves_pars;
    	  //n_moves = tree->n_moves;
    	  n_moves = MIN(n_moves,2*tree->n_otu-3);

    	  if(tree->mod->s_opt->spr_pars){
    		  printf("Doing parsimony spr - not supported\n");
    		  exit(EXIT_FAILURE);
    	  }else{
    		  best_move = Evaluate_List_Of_Regraft_Pos_Triple(tree->spr_list,n_moves,tree); //evaluate potential moves
    		  if(tree->spr_list[best_move]->lnL > tree->best_lnL + tree->mod->s_opt->min_diff_lk_move){
    			  Try_One_Spr_Move_Triple(tree->spr_list[best_move],tree);
    			  ret_val = 1;
    		  }else{
    			  Pars(tree);
    		  }
    	  }
      }
      Reset_Spr_List(tree);
  }
}


/*********************************************************/
//b_pulled = target edge, n_link = ancestral node

/*
		n_oppo_to_link
		|
		| b_pulled
		|
		n_link (b_pulled->anc)
	   / \
	  /   \
 --dir1  dir2---

After pruning:

 	 	n_oppo_to_link
		|
		| b_pulled
		|
		n_link (b_pulled->anc)
(pruned) \ b_residual
	      \


  --dir1 ----- dir2--
	(b_target - edge that was between n_link and dir1)

	Note - whether this is on dir1 or 2 depends on whether
	the number of dir1/2 was higher
	*/
int Test_All_Spr_Targets(t_edge *b_pulled, t_node *n_link, t_tree *tree)
{
  t_node *n_opp_to_link,*n_v1,*n_v2,*n_up;
  t_edge *b_target,*b_residual;
  int i,dir1,dir2;
  phydbl init_len_v1, init_len_v2, init_len_pulled;
  int best_found,approx;
  phydbl init_lnL;

//Lk here?

  init_lnL = tree->c_lnL;
  n_up = NULL;
  b_target = b_residual = NULL;
  n_opp_to_link  = (n_link == b_pulled->rght)?(b_pulled->left):(b_pulled->rght); //node on the other side of b_pulled
  approx = 1;

  init_len_pulled = b_pulled->l; //initial length of target edge
  dir1 = dir2 = -1; //nodes connected to n_link other than n_oppo_to_link
  For(i,3)
    if(n_link->v[i] != n_opp_to_link){
	if(dir1<0) dir1 = i;
	else       dir2 = i;
    }

  if(n_link->v[dir1]->num < n_link->v[dir2]->num){
      n_v1        = n_link->v[dir1];
      n_v2        = n_link->v[dir2];
      init_len_v1 = n_link->b[dir1]->l;
      init_len_v2 = n_link->b[dir2]->l;
    }else{
      n_v1        = n_link->v[dir2];
      n_v2        = n_link->v[dir1];
      init_len_v1 = n_link->b[dir2]->l;
      init_len_v2 = n_link->b[dir1]->l;
    }
  /* Pruning is meaningless otherwise */
  if(!(n_v1->tax && n_v2->tax) ){
      Prune_Subtree(n_link,n_opp_to_link,&b_target,&b_residual,tree);
      if(tree->mod->whichrealmodel==HLP17){ //update ancestors on pruned tree
    	  Update_Ancestors_Edge(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree->noeud[tree->mod->startnode]->b[0],tree);
    	  tree->both_sides=1;
    	  Lk(tree);
    	  Get_UPP(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree);
      }
      if(tree->mod->s_opt->spr_lnL){
    	  Fast_Br_Len(b_target,tree,0);

      }
      if(tree->mod->whichrealmodel==HLP17){ //update upper likelihoods on pruned tree
    	  tree->both_sides=1;
    	  Lk(tree);
    	  Get_UPP(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree);
      }
      //search through both sides of b_target to find a potential target within a certain depth range
      //importantly, branch length optimization is not performed in this search
      //recurse along a path only until you find a graft target above a certain cutoff
      best_found = 0;
      tree->depth_curr_path = 0; 
      tree->curr_path[0] = b_target->left;
      Test_One_Spr_Target_Recur(b_target->rght,
				b_target->left,
				b_pulled,n_link,b_residual,&best_found,tree);
      
      tree->depth_curr_path = 0; 
      tree->curr_path[0] = b_target->rght;
      Test_One_Spr_Target_Recur(b_target->left,
				b_target->rght,
				b_pulled,n_link,b_residual,&best_found,tree);
      
      Graft_Subtree(b_target,n_link,b_residual,tree); //graft tree back together after searching

      if((n_link->v[dir1] != n_v1) || (n_link->v[dir2] != n_v2)){
    	  PhyML_Printf("\n. Warning : -- SWITCH NEEDED -- ! \n");
      }
      
      n_link->b[dir1]->l = init_len_v1; //reset edge lengths
      n_link->b[dir2]->l = init_len_v2; 
      b_pulled->l = init_len_pulled;

      Update_PMat_At_Given_Edge(n_link->b[dir1],tree);//reset pmats on edges
      Update_PMat_At_Given_Edge(n_link->b[dir2],tree);
      Update_PMat_At_Given_Edge(b_pulled,tree);
      
      if(tree->mod->s_opt->spr_lnL){
    	  Update_P_Lk(tree,b_pulled,  n_link);//reset partial lhoods
    	  Update_P_Lk(tree,b_target,  n_link);
    	  Update_P_Lk(tree,b_residual,n_link);
      }else{
    	  Update_P_Pars(tree,b_pulled,  n_link);
    	  Update_P_Pars(tree,b_target,  n_link);
    	  Update_P_Pars(tree,b_residual,n_link);
      }

      For(i,3){
      	  if(n_link->v[i] != n_opp_to_link){ //reset partial lhoods on subtree, now that we've re-grafted the pruned tree back
      		  if(tree->mod->s_opt->spr_lnL) Pre_Order_Lk(n_link,n_link->v[i],tree);
      		  else                          Pre_Order_Pars(n_link,n_link->v[i],tree);
      	  }
      }
      if(tree->mod->whichrealmodel==HLP17){ //update ancestors and upper likelihoods
    	  Update_Ancestors_Edge(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree->noeud[tree->mod->startnode]->b[0],tree);
    	  Get_UPP(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree);//update upper lhoods across tree
      }
  }
  tree->c_lnL = init_lnL;

  return 0;

}

/*********************************************************/
//a,d: nodes on either side of b_target on larger tree
//pulled: pruned edge, link: node that was pruned, residual: residual edge, best_found:
//only recurse within the boundaries of min_depth_path and max_depth_path
//doesn't do branch length optimization when it tests a site
void Test_One_Spr_Target_Recur(t_node *a, t_node *d, t_edge *pulled, t_node *link, t_edge *residual, int *best_found, t_tree *tree)
{
  int i;

  if(*best_found) return;

  if(d->tax){
	  return; //return if you hit a tip
  }else{
      phydbl move_lnL;

      For(i,3){
       if(d->v[i] != a){ //check targets away from a
	      if(tree->mod->s_opt->spr_lnL) Update_P_Lk(tree,d->b[i],d); //update partial lhoods as you go along
	      else                          Update_P_Pars(tree,d->b[i],d);
	      
	      tree->depth_curr_path++;
	      tree->curr_path[tree->depth_curr_path] = d->v[i];
	      
	      if((tree->depth_curr_path <= tree->mod->s_opt->max_depth_path) && 
	    		  (tree->depth_curr_path >= tree->mod->s_opt->min_depth_path)){
	    	  move_lnL = Test_One_Spr_Target(d->b[i],pulled,link,residual,tree); //test pruning the subtree at this site without br len optimization
	    	  if(move_lnL > tree->best_lnL + tree->mod->s_opt->min_diff_lk_move){
	    		  *best_found = 1; //if you find a move that is above a certain cutoff, pick that one and stop recusring
	    	  }
	      }
	      if(tree->depth_curr_path < tree->mod->s_opt->max_depth_path){//keep recursing
	  		Test_One_Spr_Target_Recur(d,d->v[i],pulled,link,residual,best_found,tree);
	      }
	      tree->depth_curr_path--;
       }
     }
   }
}

/*********************************************************/
//b_target = edge we're grafting on to, b_arrow: edge that was pruned
//tests SPR grafting but doesn't do branch length optimization
//kind of a rough search first, but unclear if the best approach
phydbl Test_One_Spr_Target(t_edge *b_target, t_edge *b_arrow, t_node *n_link, t_edge *b_residual, t_tree *tree)
{
  phydbl init_target_len, init_arrow_len, init_residual_len;
  int i,dir_v0,dir_v1,dir_v2;
  phydbl l0,l1,l2;
  t_node *v1, *v2;
  phydbl init_lnL, move_lnL;
  int init_pars,move_pars;
  int approx;

  tree->n_moves++; //tally up that this is a potential move

  approx    = 1;
  move_lnL  = UNLIKELY;
  init_lnL  = tree->c_lnL;
  init_pars = tree->c_pars;
  Graft_Subtree(b_target,n_link,b_residual,tree);

  Update_PMat_At_Given_Edge(b_arrow,tree);

  init_target_len   = b_target->l;
  init_arrow_len    = b_arrow->l;
  init_residual_len = b_residual->l;

  if(tree->mod->s_opt->spr_lnL){ //upper likelihood of an edge not affected by grafting on a subtree
/*       move_lnL = Triple_Dist(n_link,tree,1); */
      Update_PMat_At_Given_Edge(b_target,tree);
      Update_PMat_At_Given_Edge(b_arrow,tree);
      Update_P_Lk(tree,b_residual,n_link);
      move_lnL = Lk_At_Given_Edge(b_residual,tree);
   //   printf("%d %d %lf %lf\n",b_target->num,b_arrow->num,move_lnL,tree->mod->omega_part[0]);
  }else{
	  printf("doing pars\n");
      Update_P_Pars(tree,b_residual,n_link);
      move_pars = Pars_At_Given_Edge(b_residual,tree);
  }

  //update direections after grafting
  v1 = (b_residual->left == n_link)?(b_residual->rght):(b_residual->left);
  v2 = (b_target->left   == n_link)?(b_target->rght):(b_target->left);
  dir_v1 = dir_v2 = dir_v0 = -1;
  For(i,3){
      if(n_link->v[i]      == v1) dir_v1 = i;
      else if(n_link->v[i] == v2) dir_v2 = i;
      else                        dir_v0 = i;
  }
  l0 = n_link->b[dir_v0]->l;
  if(n_link->v[dir_v1]->num > n_link->v[dir_v2]->num){
      l1 = n_link->b[dir_v2]->l;
      l2 = n_link->b[dir_v1]->l;
  }else{
      l1 = n_link->b[dir_v1]->l;
      l2 = n_link->b[dir_v2]->l;
  }

  //record info in SPR list
  For(i,tree->depth_curr_path+1) tree->spr_list[tree->size_spr_list]->path[i] = tree->curr_path[i];
  tree->spr_list[tree->size_spr_list]->depth_path    = tree->depth_curr_path;
  tree->spr_list[tree->size_spr_list]->pars          = tree->c_pars;
  tree->spr_list[tree->size_spr_list]->lnL           = tree->c_lnL;
  tree->spr_list[tree->size_spr_list]->b_target      = b_target;
  tree->spr_list[tree->size_spr_list]->n_link        = n_link;
  tree->spr_list[tree->size_spr_list]->n_opp_to_link = (n_link==b_arrow->left)?(b_arrow->rght):(b_arrow->left);
  tree->spr_list[tree->size_spr_list]->b_opp_to_link = b_arrow;
  tree->spr_list[tree->size_spr_list]->l0            = l0;
  tree->spr_list[tree->size_spr_list]->l1            = l1;
  tree->spr_list[tree->size_spr_list]->l2            = l2;
  tree->spr_list[tree->size_spr_list]->dist          = b_target->topo_dist_btw_edges;

  Include_One_Spr_To_List_Of_Spr(tree->spr_list[tree->size_spr_list],tree);

  //reset both trees
  b_target->l   = init_target_len;
  b_arrow->l    = init_arrow_len;
  b_residual->l = init_residual_len;

  Prune_Subtree(n_link,
		(n_link==b_arrow->left)?(b_arrow->rght):(b_arrow->left),
		&b_target,
		&b_residual,
		tree);

  if(tree->mod->s_opt->spr_lnL){
	  Update_PMat_At_Given_Edge(b_target,tree); //reset pmat at b_target
  }

  tree->c_lnL   = init_lnL; //reset tree likelihood
  tree->c_pars  = init_pars;
  return move_lnL; //return likelihood from grafting
}

/*********************************************************/
//search through list of potential targets (spr_list)
int Evaluate_List_Of_Regraft_Pos_Triple(spr **spr_list, int list_size, t_tree *tree)
{
  spr *move;
  t_edge *init_target, *b_residual;
  int i,j,best_move;
  int dir_v0, dir_v1, dir_v2;
  phydbl recorded_l;
  phydbl best_lnL,init_lnL,delta_lnL;
  phydbl max_improv;

  best_lnL = UNLIKELY;
  delta_lnL = 0.0;
  max_improv = 0.0;
  init_target = b_residual = NULL;
  best_move = -1;
  init_lnL = tree->c_lnL;
  
  if(!list_size){
      PhyML_Printf("\n. List size is 0 !");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
  }
  
  recorded_l = -1.0;
  For(i,list_size){ //loop through all moves recorded
      move = spr_list[i];
      if(!move){
    	  PhyML_Printf("\n. move is NULL\n");
    	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    	  Warn_And_Exit("");
      }
    if(move->b_target){

      upAllPmats(tree);
	  Record_Br_Len(NULL,tree);
	  /* Prune subtree */
	  if(tree->mod->whichrealmodel==HLP17){//check to make sure subtree directio nis correct
		  if(move->n_link->anc->num == move->n_opp_to_link->num){
			  printf("incorrect subtree direction %d %d\n",move->n_link->num,move->n_opp_to_link->num);
		 	  exit(EXIT_FAILURE);
		 }
      }
	  Prune_Subtree(move->n_link,move->n_opp_to_link,&init_target,&b_residual,tree);//prune subtree from move
	  if(tree->mod->whichrealmodel==HLP17){//update ancestors and upper likelihoods
		  Update_Ancestors_Edge(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree->noeud[tree->mod->startnode]->b[0],tree); //added by Ken 7/11
		  Get_UPP(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree);//update upper lhoods across tree
		  tree->both_sides=1;
		  Lk(tree);
      }
	  if(recorded_l < 0.0){
      /* Rough optimisation of the branch length at prune site
       * We only need to perform this optimisation for the first
       * element of spr_list because the pruned subtree is the
       * same across all the elements of spr_list. It would not
       * be true in the general case */
		  Fast_Br_Len(init_target,tree,0);
		  if(tree->mod->whichrealmodel==HLP17){
		  	  tree->both_sides=1;
		  	  Lk(tree);
		  }
      /* Record branch length at prune site */
		  move->init_target_l = init_target->l;
		  recorded_l          = init_target->l;
	  }else{
		  init_target->l      = recorded_l;
		  move->init_target_l = recorded_l;
	  }
	  //Update the change proba matrix at prune position
	  Update_PMat_At_Given_Edge(init_target,tree);
	  //  Update conditional likelihoods along the path from the prune to the regraft position
	  //Note - does NOT update upper likelihoods
	  Update_P_Lk_Along_A_Path(move->path,move->depth_path+1,tree);

	  /* Regraft subtree at target edge.*/
	  Graft_Subtree(move->b_target,move->n_link,b_residual,tree);
	  dir_v1 = dir_v2 = dir_v0 = -1;
	  For(j,3){
	      if(move->n_link->v[j] == move->n_opp_to_link) dir_v0 = j;
	      else if(dir_v1 < 0)                           dir_v1 = j;
	      else                                          dir_v2 = j;
	    }
	  
	  move->n_link->b[dir_v0]->l = move->l0;
	  if(move->n_link->v[dir_v1]->num > move->n_link->v[dir_v2]->num){
	      move->n_link->b[dir_v2]->l = move->l1;
	      move->n_link->b[dir_v1]->l = move->l2;
	  }else{
	      move->n_link->b[dir_v1]->l = move->l1;
	      move->n_link->b[dir_v2]->l = move->l2;
	  }

	  move->lnL = Triple_Dist(move->n_link,tree,0); //do approximate branch length optimization, and return resulting log likelihood

	  /* Record updated branch lengths for this move */
	  move->l0 = move->n_link->b[dir_v0]->l;
	  
	  if(move->n_link->v[dir_v1]->num > move->n_link->v[dir_v2]->num){
	      move->l1 = move->n_link->b[dir_v2]->l;
	      move->l2 = move->n_link->b[dir_v1]->l;
	    }else{
	      move->l1 = move->n_link->b[dir_v1]->l;
	      move->l2 = move->n_link->b[dir_v2]->l;
	    }
	  if(move->lnL > best_lnL + tree->mod->s_opt->min_diff_lk_move){ //only record a better move if it is better by a very low threshold
	      best_lnL  = move->lnL;
	      best_move = i;
	      //if(delta_lnL > max_improv) max_improv = delta_lnL;
	    }

	  /* Regraft the subtree at its original position */
	  Prune_Subtree(move->n_link,
			move->n_opp_to_link,
			&move->b_target,
			&b_residual,
			tree);
	  Graft_Subtree(init_target,
			move->n_link,
			b_residual,
			tree);
	  /* Restore branch lengths */
	  Restore_Br_Len(NULL,tree);

	  Update_PMat_At_Given_Edge(move->b_target,tree);

	  tree->c_lnL = init_lnL;
	}
      /* Bail out as soon as you've found a true improvement */
      if(move->lnL > tree->best_lnL + tree->mod->s_opt->min_diff_lk_move) break;
    }

  //restore all of the partial lhoods and things that were disturbed
  	For(i,list_size){
      move = spr_list[i];
      if(move->b_target){
	  For(j,3) Update_PMat_At_Given_Edge(move->n_link->b[j],tree);
	  if(tree->mod->whichrealmodel==HLP17){
		  Update_Ancestors_Edge(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree->noeud[tree->mod->startnode]->b[0],tree); //added by Ken 7/11
	  }
	  For(j,3) Update_P_Lk(tree,move->n_link->b[j],move->n_link);
	  /* TO DO : we don't need to update all these partial likelihoods here.
	     Would need to record only those that were along the paths examined
	     above */
	  For(j,3)
	    if(move->n_link->v[j] != move->n_opp_to_link)
	      Pre_Order_Lk(move->n_link,move->n_link->v[j],tree);
	  	  break;
		}
    }
  	//update ancestors and upper likelihoods
    if(tree->mod->whichrealmodel==HLP17){ //update ancestors
      Update_Ancestors_Edge(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree->noeud[tree->mod->startnode]->b[0],tree);
      Get_UPP(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree);//update upper lhoods across tree
    }

#ifdef DEBUG
  if(best_move < 0)
    {
      PhyML_Printf("\n\n. Best_move < 0 !");

      PhyML_Printf("\n. List size = %d",list_size);
      For(i,list_size)
	{
	  move = spr_list[i];
	  PhyML_Printf("\n. %p %p",move,move->b_target);
	}

      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
#endif
  //exit(EXIT_FAILURE);
  return best_move;
}

/*********************************************************/

int Try_One_Spr_Move_Triple(spr *move, t_tree *tree)
{
  t_edge *init_target, *b_residual;
  int j;
  int dir_v0, dir_v1, dir_v2;

  upAllPmats(tree);
  Record_Br_Len(NULL,tree);

  //perform best SPR move identified before
  Prune_Subtree(move->n_link,
		move->n_opp_to_link,
		&init_target,
		&b_residual,
		tree);

  init_target->l = move->init_target_l;

  Graft_Subtree(move->b_target,move->n_link,b_residual,tree);

  dir_v1 = dir_v2 = dir_v0 = -1;
  For(j,3){
      if(move->n_link->v[j] == move->n_opp_to_link) dir_v0 = j;
      else if(dir_v1 < 0)                           dir_v1 = j;
      else                                          dir_v2 = j;
    }

  move->n_link->b[dir_v0]->l = move->l0;

  if(move->n_link->v[dir_v1]->num > move->n_link->v[dir_v2]->num){
      move->n_link->b[dir_v2]->l = move->l1;
      move->n_link->b[dir_v1]->l = move->l2;
 }else{
      move->n_link->b[dir_v1]->l = move->l1;
      move->n_link->b[dir_v2]->l = move->l2;
 }

  if(move->lnL > tree->best_lnL + tree->mod->s_opt->min_diff_lk_move) /* Apply the move */
    {
      #if defined OMP || defined BLAS_OMP

      tree->t_current=omp_get_wtime();

      #else
  
      time(&(tree->t_current));
  
      #endif
      tree->both_sides = 1;
      Lk(tree);
      Pars(tree);
      if(tree->mod->whichrealmodel==HLP17){
        Update_Ancestors_Edge(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree->noeud[tree->mod->startnode]->b[0],tree);
        Get_UPP(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree);//update upper lhoods across tree
      }

      if((tree->mod->s_opt->print) && (!tree->io->quiet)){
    	  Print_Lk(tree,"[Topology           ]");
	  	  PhyML_Printf("[depth=%5d]",move->depth_path); fflush(NULL);
	  	  Print_Trace(tree);
      }

      upAllPmats(tree);
      tree->n_improvements++;
      tree->best_lnL = tree->c_lnL;
      Record_Br_Len(NULL,tree);

      //increase the deepest path searched for
      if(move->depth_path > tree->mod->s_opt->deepest_path){
    	  tree->mod->s_opt->deepest_path = move->depth_path;
      }
      return 1;
    }


  Prune_Subtree(move->n_link,
		move->n_opp_to_link,
		&move->b_target,
		&b_residual,
		tree);

  Graft_Subtree(init_target,
		move->n_link,
		b_residual,
		tree);

  Restore_Br_Len(NULL,tree);
  tree->both_sides = 1;
  Lk(tree);
  Pars(tree);
  if(tree->mod->whichrealmodel==HLP17){
    Update_Ancestors_Edge(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree->noeud[tree->mod->startnode]->b[0],tree);
    Get_UPP(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree);//update upper lhoods across tree
  }
  return 0;
}
/*********************************************************/

void Speed_Spr_Loop(t_tree *tree)
{
  phydbl lk_old;
  int init_thresh;
  init_thresh                      = tree->mod->s_opt->pars_thresh;
  tree->best_pars                  = 1E+8;
  tree->mod->s_opt->spr_lnL        = 1; //changed by Ken 17/2/2017
  tree->mod->s_opt->spr_pars       = 0;
  tree->mod->s_opt->quickdirty     = 0;

  if((tree->mod->s_opt->print) && (!tree->io->quiet)) PhyML_Printf("\n. Maximizing likelihood (using SPR moves)...\n");

  Lk(tree);//!< Added by Marcelo.
  Print_Lk(tree,"[Initial Tree       ]");//!< Added by Marcelo.
  if(tree->io->print_trace){
	  Print_Trace(tree);
  }

  if(tree->io->opt_heuristic_manuel==NO) //!Added by Marcelo
  {
  Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));
  tree->best_lnL = tree->c_lnL;
  if(tree->io->print_trace){
	  Print_Trace(tree);
  }
  //added by Ken 8/3/2017
  /* Optimise branch lengths */
	  int startnode = 0;
	  if(tree->mod->whichrealmodel==HLP17){
		  startnode=tree->mod->startnode;
		  tree->both_sides = 1;
		  Lk(tree);
		  Update_Ancestors_Edge(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree->noeud[tree->mod->startnode]->b[0],tree); //added by Ken 7/11
		  Get_UPP(tree->noeud[startnode], tree->noeud[startnode]->v[0], tree);
	  }

	  Optimize_Br_Len_Serie(tree->noeud[startnode],
			tree->noeud[startnode]->v[0],
			tree->noeud[startnode]->b[0],
			tree,
			tree->data);

	  /* Update partial likelihoods */
	  tree->both_sides = 1;
	  Lk(tree);
	  if(tree->mod->whichrealmodel==HLP17){
		  Get_UPP(tree->noeud[startnode], tree->noeud[startnode]->v[0], tree);
	  }
	  if(tree->io->print_trace){
		  Print_Trace(tree);
	  }
	  Print_Lk(tree,"[Branch lengths     ]");

  /*****************************/
  lk_old = UNLIKELY;
  tree->mod->s_opt->max_depth_path = 2*tree->n_otu-3;
  tree->mod->s_opt->spr_lnL        = 1;
  do
    {
      lk_old = tree->c_lnL;
      Speed_Spr(tree,1);
      if(tree->n_improvements){ //if there were improvements after the loop, optimize parameters
    	  Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));
      }
      if((!tree->n_improvements) || (FABS(lk_old-tree->c_lnL) < 1.)){
    	  break;//if there are no improvements, break
      }
    }
  while(1);


    /*****************************/
    lk_old = UNLIKELY;
    do
    {
      lk_old = tree->c_lnL;
      Simu(tree,10);
      Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));
    }
    while(FABS(lk_old - tree->c_lnL) > tree->mod->s_opt->min_diff_lk_global);
    /*****************************/

    /*****************************/
    do
    {
      if(!Check_NNI_Five_Branches(tree)) break;
    }while(1);
    /*****************************/

    /*   if((tree->mod->s_opt->print) && (!tree->io->quiet)) PhyML_Printf("\n"); */
  }
  else if(tree->io->opt_heuristic_manuel==YES) //!Added by Marcelo
  {

    Round_Optimize(tree,tree->data,tree->io->roundMax_start); //! last parameter: max number of rounds

    tree->best_lnL = tree->c_lnL;
    /*****************************/
    lk_old = UNLIKELY;
    tree->mod->s_opt->max_depth_path = 2*tree->n_otu-3;
    tree->mod->s_opt->spr_lnL        = 0;
    do
    {
      lk_old = tree->c_lnL;
      Speed_Spr(tree,1);
      if((!tree->n_improvements) || (FABS(lk_old-tree->c_lnL) < 1.)) break;
    }
    while(1);
    /*****************************/

    /*****************************/
    lk_old = UNLIKELY;
    do
    {
      lk_old = tree->c_lnL;
      Simu(tree,10);
    }
    while(FABS(lk_old - tree->c_lnL) > tree->mod->s_opt->min_diff_lk_global);
    /*****************************/

    Round_Optimize(tree,tree->data,tree->io->roundMax_end); //! last parameter: max number of rounds

    /*****************************/
    do
    {
      if(!Check_NNI_Five_Branches(tree)) break;
    }while(1);
    /*****************************/


    /*   if((tree->mod->s_opt->print) && (!tree->io->quiet)) PhyML_Printf("\n"); */

  }

}
/*********************************************************/

void Speed_Spr(t_tree *tree, int max_cycles)
{
  int step,old_pars;
  phydbl old_lnL;

  if(tree->lock_topo)
    {
      PhyML_Printf("\n. The tree topology is locked.");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }


  tree->both_sides = 1;
  Pars(tree);
  Lk(tree);
  if(tree->mod->whichrealmodel==HLP17){
	  int startnode=tree->mod->startnode;
	  Update_Ancestors_Edge(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree->noeud[tree->mod->startnode]->b[0],tree); //added by Ken 7/11
	  Get_UPP(tree->noeud[startnode], tree->noeud[startnode]->v[0], tree);
  }
  Record_Br_Len(NULL,tree);

  tree->mod->s_opt->deepest_path  = 0;
  tree->best_pars                 = tree->c_pars;
  tree->best_lnL                  = tree->c_lnL;
  old_lnL                         = tree->c_lnL;
  old_pars                        = tree->c_pars;
  step                            = 0;

  do
    {
      ++step;

      old_lnL  = tree->c_lnL;
      old_pars = tree->c_pars;


      tree->n_improvements         = 0;
      tree->perform_spr_right_away = 1;
      Spr(UNLIKELY,tree);
      if(!tree->mod->s_opt->spr_pars){

	  /* Optimise branch lengths */
    	  int startnode = 0;
    	  if(tree->mod->whichrealmodel==HLP17){
    		  startnode=tree->mod->startnode;
    		  tree->both_sides = 1;
    		  Lk(tree);
    		  Update_Ancestors_Edge(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree->noeud[tree->mod->startnode]->b[0],tree); //added by Ken 7/11
    		  Get_UPP(tree->noeud[startnode], tree->noeud[startnode]->v[0], tree);
    	  }

    	  Optimize_Br_Len_Serie(tree->noeud[startnode],
				tree->noeud[startnode]->v[0],
				tree->noeud[startnode]->b[0],
				tree,
				tree->data);

    	  /* Update partial likelihoods */
    	  tree->both_sides = 1;
    	  Lk(tree);

    	  /* Print log-likelihood and parsimony scores */
    	  if((tree->mod->s_opt->print) && (!tree->io->quiet)) Print_Lk(tree,"[Branch lengths     ]");
      }else{
    	  if((tree->mod->s_opt->print) && (!tree->io->quiet)) Print_Pars(tree);
      }

      Pars(tree);
      if(tree->io->print_trace){Print_Trace(tree);}

      /* Record the current best log-likelihood and parsimony */
      tree->best_lnL  = tree->c_lnL;
      tree->best_pars = tree->c_pars;

      if(!tree->mod->s_opt->spr_pars){
    	  if(tree->io->datatype==CODON){
    		  if(tree->c_lnL < old_lnL-tree->io->mod->s_opt->min_diff_lk_codonModels){
    			  PhyML_Printf("\n. old_lnL = %f c_lnL = %f",old_lnL,tree->c_lnL);
    			  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    			  Warn_And_Exit("");
    		  }
    	  }else{
    		  if(tree->c_lnL < old_lnL-tree->mod->s_opt->min_diff_lk_local){
    			  PhyML_Printf("\n. old_lnL = %f c_lnL = %f",old_lnL,tree->c_lnL);
    			  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    			  Warn_And_Exit("");
    		  }
    	  }
      }else{
    	  if(tree->c_pars > old_pars){
    		  PhyML_Printf("\n. old_pars = %d c_pars = %d",old_pars,tree->c_pars);
    		  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    		  Warn_And_Exit("");
    	  }
      }

      /* Record the current best branch lengths  */
      upAllPmats(tree);
      Record_Br_Len(NULL,tree);

      /* Exit if no improvements after complete optimization */
      if(step+1 > max_cycles) break;
      if((!tree->mod->s_opt->spr_pars) && (FABS(old_lnL-tree->c_lnL) < tree->mod->s_opt->min_diff_lk_global)) break;
/*    if(( tree->mod->s_opt->spr_pars) && (FABS(old_pars-tree->c_pars) < 1.)) break; */
      if(!tree->n_improvements) break;
    }while(1);
}


void Include_One_Spr_To_List_Of_Spr(spr *move, t_tree *tree)
{
  int i;
  spr *buff_spr;

  if((( tree->mod->s_opt->spr_lnL) && (move->lnL  > tree->spr_list[tree->size_spr_list-1]->lnL)) ||
     ((!tree->mod->s_opt->spr_lnL) && (move->pars <= tree->spr_list[tree->size_spr_list-1]->pars)))
    {
      tree->spr_list[tree->size_spr_list-1]->depth_path    = move->depth_path;
      tree->spr_list[tree->size_spr_list-1]->pars          = move->pars;
      tree->spr_list[tree->size_spr_list-1]->lnL           = move->lnL;
      tree->spr_list[tree->size_spr_list-1]->b_target      = move->b_target;
      tree->spr_list[tree->size_spr_list-1]->n_link        = move->n_link;
      tree->spr_list[tree->size_spr_list-1]->n_opp_to_link = move->n_opp_to_link;
      tree->spr_list[tree->size_spr_list-1]->b_opp_to_link = move->b_opp_to_link;
      tree->spr_list[tree->size_spr_list-1]->l0            = move->l0;
      tree->spr_list[tree->size_spr_list-1]->l1            = move->l1;
      tree->spr_list[tree->size_spr_list-1]->l2            = move->l2;
      tree->spr_list[tree->size_spr_list-1]->dist          = move->dist;

      For(i,tree->spr_list[tree->size_spr_list-1]->depth_path+1)
	tree->spr_list[tree->size_spr_list-1]->path[i] = move->path[i];

      for(i=tree->size_spr_list-1;i>0;i--)
	{
	  if((( tree->mod->s_opt->spr_lnL) && (tree->spr_list[i]->lnL > tree->spr_list[i-1]->lnL)) ||
	     ((!tree->mod->s_opt->spr_lnL) && (tree->spr_list[i]->pars <=  tree->spr_list[i-1]->pars)))
	    {
	      buff_spr            = tree->spr_list[i-1];
	      tree->spr_list[i-1] = tree->spr_list[i];
	      tree->spr_list[i]   = buff_spr;
	    }
	  else  break;
	}
    }
}

/*********************************************************/

void Reset_Spr_List(t_tree *tree)
{
  int i;

  For(i,tree->size_spr_list)
    {
      tree->spr_list[i]->depth_path     = 0;
      tree->spr_list[i]->pars           = MAX_PARS;
      tree->spr_list[i]->lnL            = UNLIKELY;
      tree->spr_list[i]->n_link         = NULL;
      tree->spr_list[i]->n_opp_to_link  = NULL;
      tree->spr_list[i]->b_target       = NULL;
    }
}

/*********************************************************/

void Make_Best_Spr(t_tree *tree)
{
  tree->best_spr = Make_One_Spr(tree);
  Init_One_Spr(tree->best_spr);
}

/*********************************************************/

void Make_Spr_List(t_tree *tree)
{
  int i;

  tree->size_spr_list = 2*tree->n_otu-3;
  tree->spr_list = (spr **)mCalloc(2*tree->n_otu-2,sizeof(spr *));

  For(i,2*tree->n_otu-2)
    {
      tree->spr_list[i] = Make_One_Spr(tree);
      Init_One_Spr(tree->spr_list[i]);
    }
  tree->perform_spr_right_away = 0;
}

/*********************************************************/

void Init_One_Spr(spr *a_spr)
{
  a_spr->lnL             = UNLIKELY;
  a_spr->pars            = 1E+5;
  a_spr->depth_path      = 0;
  a_spr->dist            = 0;
  a_spr->init_target_l   = -1.;
  a_spr->l0              = -1.;
  a_spr->l1              = -1.;
  a_spr->l2              = -1.;
  a_spr->n_link          = NULL;
  a_spr->n_opp_to_link   = NULL;
  a_spr->b_opp_to_link   = NULL;
  a_spr->b_target        = NULL;
  a_spr->b_init_target   = NULL;
}

/*********************************************************/

spr *Make_One_Spr(t_tree *tree)
{
  spr *a_spr;
  a_spr       = (spr *)mCalloc(1,sizeof(spr));
  a_spr->path = (t_node **)mCalloc(tree->n_otu,sizeof(t_node *));
  return a_spr;
}

/*********************************************************/

void Spr_Pars(t_tree *tree)
{

  PhyML_Printf("\n. Minimizing parsimony...\n");

  tree->best_pars = 1E+8;
  tree->best_lnL  = UNLIKELY;
  tree->mod->s_opt->spr_lnL = 0;
  
  tree->mod->s_opt->spr_pars = 1;
  do
    {
      Speed_Spr(tree,1);
    }while(tree->n_improvements);  
  tree->mod->s_opt->spr_pars = 0;
  
  PhyML_Printf("\n");
}

/*********************************************************/


void Print_Trace(t_tree *tree){

	PhyML_Fprintf(tree->io->fp_out_tree_trace,"[%d,%f]%s\n",tree->io->tracecount,tree->c_lnL,Write_Tree(tree)); fflush(tree->io->fp_out_tree_trace);
	//if((tree->io->print_site_lnl) && (!tree->mod->s_opt->spr_pars)) Print_Site_Lk(tree,tree->io->fp_out_lk);
	fflush(tree->io->fp_out_lk);

	//print header of stats file
	if(tree->io->tracecount == 0){
		PhyML_Fprintf(tree->io->fp_out_stats_trace,"Index\tTime\tLnL\tKappa");
		if(tree->mod->whichrealmodel == HLP17){
			int omegai;
			for(omegai=0;omegai<tree->mod->nomega_part;omegai++){
				PhyML_Fprintf(tree->io->fp_out_stats_trace,"\tOmega%d",omegai);
			}

			int mot;
			for(mot=0;mot<tree->io->mod->nmotifs;mot++){
				PhyML_Fprintf(tree->io->fp_out_stats_trace,"\t%s:%d;%d",tree->io->mod->motifs[mot],tree->io->mod->motif_hotness[mot],tree->io->mod->hoptindex[tree->io->mod->motif_hotness[mot]]);
          	}

		}
		PhyML_Fprintf(tree->io->fp_out_stats_trace,"\tT1\tC1\tA1\tG1\tT2\tC2\tA2\tG2\tT3\tC3\tA3\tG3\n");

	}

	//print out stats
	PhyML_Fprintf(tree->io->fp_out_stats_trace,"\t%d\t%f\t%f",tree->io->tracecount,(omp_get_wtime()-tree->t_beg),tree->c_lnL);
	PhyML_Fprintf(tree->io->fp_out_stats_trace,"\t%f",tree->mod->kappa);
	if(tree->mod->whichrealmodel ==HLP17){
		int omegai;
		for(omegai=0;omegai<tree->mod->nomega_part;omegai++){
			PhyML_Fprintf(tree->io->fp_out_stats_trace,"\t%lf",tree->mod->omega_part[omegai]);
		}

		int mot;
		for(mot=0;mot<tree->io->mod->nmotifs;mot++){
			PhyML_Fprintf(tree->io->fp_out_stats_trace,"\t%lf",tree->io->mod->hotness[tree->io->mod->motif_hotness[mot]]);
         }

	}
	int c;
	for(c=0;c<12;c++){
		PhyML_Fprintf(tree->io->fp_out_stats_trace,"\t%lf",tree->mod->base_freq[c]);

	}
	PhyML_Fprintf(tree->io->fp_out_stats_trace,"\n");
	fflush(tree->io->fp_out_stats_trace);

	tree->io->tracecount++;
}



/*********************************************************/








/*
** EOF: spr.c
*/
