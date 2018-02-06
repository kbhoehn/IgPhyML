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

#include "utilities.h"

int cf3x4Ind[9]={0,1,2,4,5,6,8,9,10}; //!< Very useful for CF3X4 ... see below.

/*********************************************************/

FILE* GetTreeFile(option *io) {
    FILE *treefile = NULL;
    char *filename;
    char *filenr;
    if(io->logtree == 0) {
        Warn_And_Exit("(in GetTreeFile), should not happen");
    } else if(io->logtree == 1) {
        treefile = openOutputFile(io->out_boot_tree_file, "_igphyml_logtree", ".txt", io);
    } else if(io->logtree == 2) {
        filename = mCalloc(30, sizeof(char));
        strcat(filename, "_igphyml_logtree");
        filenr = mCalloc(8, sizeof(char));
        sprintf(filenr, "_%d", io->treecounter);
        strcat(filename, filenr);
        treefile = openOutputFile(io->out_boot_tree_file, filename, ".txt", io);
        io->treecounter++;
        free(filename);
        free(filenr);
    }
    return(treefile);
}

t_tree *Read_Tree(char *s_tree)
{
  char **subs;
  int i,n_ext,n_int,n_otu;
  t_tree *tree;
  int degree;
  t_node *root_node;
  
  n_int = n_ext = 0;
  
  n_otu=0;
  For(i,(int)strlen(s_tree)) if(s_tree[i] == ',') n_otu++;
  n_otu+=1;
  
  tree = (t_tree *)Make_Tree(n_otu);
  Init_Tree(tree,tree->n_otu);
  Make_All_Tree_Nodes(tree);
  Make_All_Tree_Edges(tree);
  Make_Tree_Path(tree);
  
  subs = Sub_Trees(s_tree,&degree);
  Clean_Multifurcation(subs,degree,3);
  if(degree == 2) 
  {
    /*       Unroot_Tree(subs); */
    /*       degree = 3; */
    /*       root_node = tree->noeud[n_otu]; */
    root_node      = tree->noeud[2*n_otu-2];
    root_node->num = 2*n_otu-2;
    tree->n_root   = root_node;
    n_int         -= 1;
  }
  else
  {
    root_node      = tree->noeud[n_otu];
    root_node->num = n_otu;
    tree->n_root   = NULL;
  }
  
  root_node->tax = 0;
  
  tree->has_branch_lengths = 0;
  tree->num_curr_branch_available = 0;
  For(i,degree) R_rtree(s_tree,subs[i],root_node,tree,&n_int,&n_ext);
  
  if(tree->n_root)
  {
    tree->e_root = tree->t_edges[tree->num_curr_branch_available];
    
    tree->n_root->v[0]->v[0] = tree->n_root->v[1];
    tree->n_root->v[1]->v[0] = tree->n_root->v[0];
    
    Connect_One_Edge_To_Two_Nodes(tree->n_root->v[0],
                                  tree->n_root->v[1],
                                  tree->e_root,
                                  tree);      
    tree->e_root->l = tree->n_root->l[0] + tree->n_root->l[1];
    tree->n_root_pos = tree->n_root->l[0] / tree->e_root->l;
  }
  
  For(i,NODE_DEG_MAX) free(subs[i]);
  free(subs);
  return tree;
}

/*********************************************************/


/*********************************************************/
/* 'a' in t_node a stands for ancestor. 'd' stands for descendant */ 
void R_rtree(char *s_tree_a, char *s_tree_d, t_node *a, t_tree *tree, int *n_int, int *n_ext)
{
  int i;
  t_node *d;
  int n_otu = tree->n_otu;
  
  if(strstr(s_tree_a," ")) 
  {
    PhyML_Printf("\n. [%s]",s_tree_a);
    Warn_And_Exit("\n. Err: the tree must not contain a ' ' character\n");
  }
  
  if(s_tree_d[0] == '(')
  {
    char **subs;
    int degree;
    
    (*n_int)+=1;
    d      = tree->noeud[n_otu+*n_int];
    d->num = n_otu+*n_int;
    d->tax = 0;
    
    Read_Branch_Label(s_tree_d,s_tree_a,tree->t_edges[tree->num_curr_branch_available]);
    Read_Branch_Length(s_tree_d,s_tree_a,tree);
    
    For(i,3)
    {
      if(!a->v[i])
	    {
	      a->v[i]=d;
	      d->l[0]=tree->t_edges[tree->num_curr_branch_available]->l;
	      a->l[i]=tree->t_edges[tree->num_curr_branch_available]->l;
	      break;
	    }
    }
    d->v[0]=a;
    
    if(a != tree->n_root)
    {
      Connect_One_Edge_To_Two_Nodes(a,d,tree->t_edges[tree->num_curr_branch_available],tree);
      tree->num_curr_branch_available++;
    }
    
    subs=Sub_Trees(s_tree_d,&degree);
    Clean_Multifurcation(subs,degree,2);
    R_rtree(s_tree_d,subs[0],d,tree,n_int,n_ext);
    R_rtree(s_tree_d,subs[1],d,tree,n_int,n_ext);
    For(i,NODE_DEG_MAX) free(subs[i]);
    free(subs);
  }
  
  else
  {
    int i;
    
    d      = tree->noeud[*n_ext];
    d->tax = 1;
    
    Read_Node_Name(d,s_tree_d,tree);
    Read_Branch_Label(s_tree_d,s_tree_a,tree->t_edges[tree->num_curr_branch_available]); 
    Read_Branch_Length(s_tree_d,s_tree_a,tree);
    
    For(i,3)
    {
      if(!a->v[i])
      {
        a->v[i]=d;
        d->l[0]=tree->t_edges[tree->num_curr_branch_available]->l;
        a->l[i]=tree->t_edges[tree->num_curr_branch_available]->l;
        break;
      }
    }
    d->v[0]=a;
    
    if(a != tree->n_root)
    {
      Connect_One_Edge_To_Two_Nodes(a,d,tree->t_edges[tree->num_curr_branch_available],tree);
      tree->num_curr_branch_available++;
    }
    
    d->num=*n_ext;
    (*n_ext)+=1;
  }
}

/*********************************************************/

void Read_Branch_Label(char *s_d, char *s_a, t_edge *b)
{
  char *sub_tp;
  char *p;
  int i,pos;
  
  sub_tp = (char *)mCalloc(3+(int)strlen(s_d)+1,sizeof(char));
  
  sub_tp[0] = '(';
  sub_tp[1] = '\0';
  strcat(sub_tp,s_d);
  strcat(sub_tp,"#");
  p = strstr(s_a,sub_tp);
  if(!p)
  {
    sub_tp[0] = ',';
    sub_tp[1] = '\0';
    strcat(sub_tp,s_d);
    strcat(sub_tp,"#");
    p = strstr(s_a,sub_tp);
  }
  
  i = 0;
  b->n_labels = 0;
  if(p)
  {
    if(!(b->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(b);
    b->n_labels++;
    
    pos = 0;
    do 
    {
      b->labels[b->n_labels-1][pos] = p[i+strlen(s_d)+1];
      i++;
      pos++;
      if(p[i+strlen(s_d)+1] == '#') 
	    { 
	      b->labels[b->n_labels-1][pos] = '\0';
	      b->n_labels++;
	      if(!(b->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(b);
	      i++;
	      pos=0;
	    }
    }
    while((p[i+strlen(s_d)+1] != ':') && 
          (p[i+strlen(s_d)+1] != ',') && 
          (p[i+strlen(s_d)+1] != '('));
    
    b->labels[b->n_labels-1][pos] = '\0';
  }
  
  if(p)
  {
    if(b->n_labels == 1)
      PhyML_Printf("\n. Found label '%s' on t_edge %3d.",b->labels[0],b->num);
    else
    {
      PhyML_Printf("\n. Found labels ");
      For(i,b->n_labels) PhyML_Printf("'%s' ",b->labels[i]);
      PhyML_Printf("on t_edge %3d.",b->num);
    }
  }
  
  free(sub_tp);
}

/*********************************************************/

void Read_Branch_Length(char *s_d, char *s_a, t_tree *tree)
{
  char *sub_tp;
  char *p;
  t_edge *b;
  int i;
  
  b = tree->t_edges[tree->num_curr_branch_available];
  
  /*   sub_tp = (char *)mCalloc(T_MAX_LINE,sizeof(char)); */
  sub_tp = (char *)mCalloc(4+(int)strlen(s_d)+1,sizeof(char));
  
  For(i,b->n_labels)
  {
    strcat(s_d,"#");
    strcat(s_d,b->labels[i]);
  }
  
  sub_tp[0] = '(';
  sub_tp[1] = '\0';
  strcat(sub_tp,s_d);
  strcat(sub_tp,":");
  p = strstr(s_a,sub_tp);
  if(!p)
  {
    sub_tp[0] = ',';
    sub_tp[1] = '\0';
    strcat(sub_tp,s_d);
    strcat(sub_tp,":");
    p = strstr(s_a,sub_tp);
  }
  
  if(p)
  {
    b->l = atof((char *)p+(int)strlen(sub_tp));
    tree->has_branch_lengths = 1;
  }      
  free(sub_tp);
}

/*********************************************************/

void Read_Node_Name(t_node *d, char *s_tree_d, t_tree *tree)
{
  int i;
  
  if(!tree->t_edges[tree->num_curr_branch_available]->n_labels)
  {
    /*if(!d->name) d->name = (char *)mCalloc(strlen(s_tree_d)+1,sizeof(char ));*///!< Commented by Marcelo.
    strcpy(d->name,s_tree_d);
  }
  else
  {
    i = 0;
    do
    {
      /*if(!d->name) d->name = (char *)realloc(d->name,(i+1)*sizeof(char ));*///!< Commented by Marcelo.
      d->name[i] = s_tree_d[i];
      i++;
    }
    while(s_tree_d[i] != '#');
    d->name[i] = '\0';
  }
  strcpy(d->ori_name , d->name);//!< Changed by Marcelo.
  
}
/*********************************************************/

void Clean_Multifurcation(char **subtrees, int current_deg, int end_deg)
{
  
  if(current_deg <= end_deg) return;
  else
  {
    char *s_tmp;
    int i;
    
    /*       s_tmp = (char *)mCalloc(T_MAX_LINE,sizeof(char)); */
    s_tmp = (char *)mCalloc(6+
                            (int)strlen(subtrees[0])+1+
                            (int)strlen(subtrees[1])+1,
                            sizeof(char));
    
    strcat(s_tmp,"(\0");
    strcat(s_tmp,subtrees[0]);
    strcat(s_tmp,",\0");
    strcat(s_tmp,subtrees[1]);
    strcat(s_tmp,")\0");
    free(subtrees[0]);
    subtrees[0] = s_tmp;
    
    for(i=1;i<current_deg-1;i++) strcpy(subtrees[i],subtrees[i+1]);
    
    Clean_Multifurcation(subtrees,current_deg-1,end_deg);
  }
}

/*********************************************************/

char **Sub_Trees(char *tree, int *degree)
{
  char **subs;
  int posbeg,posend;
  int i;
  
  if(tree[0] != '(') {*degree = 1; return NULL;}
  
  subs=(char **)mCalloc(NODE_DEG_MAX,sizeof(char *));
  
  For(i,NODE_DEG_MAX) subs[i]=(char *)mCalloc((int)strlen(tree)+1,sizeof(char));
  
  
  posbeg=posend=1;
  (*degree)=0;
  do
  {
    posbeg = posend;
    if(tree[posend] != '(')
    {
      while((tree[posend] != ',' ) &&
            (tree[posend] != ':' ) &&
            (tree[posend] != '#' ) &&
            (tree[posend] != ')' )) 
	    {
	      posend++ ;
	    }
      posend -= 1;
    }
    else posend=Next_Par(tree,posend);
    
    while((tree[posend+1] != ',') &&
          (tree[posend+1] != ':') &&
          (tree[posend+1] != '#') &&
          (tree[posend+1] != ')')) {posend++;}
    
    
    strncpy(subs[(*degree)],tree+posbeg,posend-posbeg+1);
    /*       strcat(subs[(*degree)],"\0"); */
    subs[(*degree)][posend-posbeg+1]='\0'; /* Thanks to Jean-Baka Domelevo-Entfellner */
    
    posend += 1;
    while((tree[posend] != ',') &&
          (tree[posend] != ')')) {posend++;}
    posend+=1;
    
    
    (*degree)++;
    if((*degree) == NODE_DEG_MAX)
    {
      For(i,(*degree))
	    PhyML_Printf("\n. Subtree %d : %s\n",i+1,subs[i]);
      
      PhyML_Printf("\n. The degree of a t_node cannot be greater than %d\n",NODE_DEG_MAX);
      Warn_And_Exit("\n");
    }
  }
  while(tree[posend-1] != ')');
  
  return subs;
}


/*********************************************************/
int Next_Par(char *s, int pos)
{
  int curr;
  
  curr=pos+1;
  
  while(*(s+curr) != ')')
  {
    if(*(s+curr) == '(') curr=Next_Par(s,curr);
    curr++;
  }
  
  return curr;
}

/*********************************************************/

void Print_Tree(FILE *fp, t_tree *tree)
{
  char *s_tree;
  int i;
  
  s_tree = (char *)Write_Tree(tree);
  
  if(OUTPUT_TREE_FORMAT == 0) PhyML_Fprintf(fp,"%s\n",s_tree);
  else if(OUTPUT_TREE_FORMAT == 1)
  {
    PhyML_Fprintf(fp,"#NEXUS\n");
    PhyML_Fprintf(fp,"BEGIN TREES;\n");
    PhyML_Fprintf(fp,"\tTRANSLATE\n");
    For(i,tree->n_otu) PhyML_Fprintf(fp,"\t%3d\t%s,\n",i+1,tree->noeud[i]->name);
    PhyML_Fprintf(fp,"\tUTREE PAUP_1=\n");
    PhyML_Fprintf(fp,"%s\n",s_tree);
    PhyML_Fprintf(fp,"ENDBLOCK;");
  }
  free(s_tree);
}

/*********************************************************/

char *Write_Tree(t_tree *tree)
{
  char *s;
  int i,available;
  int pos;
  
  s=(char *)mCalloc((int)T_MAX_NAME,sizeof(char));
  available = (int)T_MAX_NAME-1;
  
  
  s[0]='(';
  pos = 1;
  
#ifdef PHYML 
  tree->n_root = NULL;
  tree->e_root = NULL;
#endif
  
  if(tree->mod->whichrealmodel == HLP17){ //added by Ken 16/2/2017 to output as a tree ith rooted branch length of zero
	  t_node* r = Make_Node_Light(-1);
	  tree->noeud[tree->mod->startnode]->v[1] = r;
	  t_edge* blank = (t_edge *)mCalloc(1,sizeof(t_edge));
	  blank->l=0.00001;
	  blank->left=tree->noeud[tree->mod->startnode];
	  blank->rght=r;
	  tree->noeud[tree->mod->startnode]->b[1] = blank;
	  r->b[0] = blank;
	  r->tax=1;
	  tree->n_root=r;
	  tree->e_root=blank;
	  tree->n_root_pos=tree->mod->startnode;
	  r->name = (char*)mCalloc(T_MAX_OPTION,sizeof(char*));
	  strcpy(r->name,tree->noeud[tree->mod->startnode]->name);
	  R_wtree(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],&available,&s,&pos,tree);
	  R_wtree(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[1],&available,&s,&pos,tree);

	  free(r);
	  free(r->name);
	  free(blank);
	  tree->noeud[tree->mod->startnode]->b[1]=NULL;
	  tree->noeud[tree->mod->startnode]->v[1]=NULL;
	  tree->n_root=NULL;
	  tree->e_root=NULL;
  }else{
	  if(!tree->n_root){
	  i = 0;
    	while((!tree->noeud[tree->n_otu+i]->v[0]) ||
    		(!tree->noeud[tree->n_otu+i]->v[1]) ||
			(!tree->noeud[tree->n_otu+i]->v[2])) i++;

    	R_wtree(tree->noeud[tree->n_otu+i],tree->noeud[tree->n_otu+i]->v[0],&available,&s,&pos,tree);
    	R_wtree(tree->noeud[tree->n_otu+i],tree->noeud[tree->n_otu+i]->v[1],&available,&s,&pos,tree);
    	R_wtree(tree->noeud[tree->n_otu+i],tree->noeud[tree->n_otu+i]->v[2],&available,&s,&pos,tree);
  	  }else{
  		R_wtree(tree->n_root,tree->n_root->v[0],&available,&s,&pos,tree);
    	R_wtree(tree->n_root,tree->n_root->v[1],&available,&s,&pos,tree);
  	  }
  }
  
  s[pos-1]=')';
  s[pos]=';';
  s[pos+1]='\0';
  
  return s;
}

/*********************************************************/
//modified by ken
void R_wtree(t_node *pere, t_node *fils, int *available, char **s_tree, int *pos, t_tree *tree)
{
  int i,p,ori_len;
  
  p = -1;
  if(fils->tax)
  {
    /*       printf("\n- Writing on %p",*s_tree); */
    ori_len = *pos;
    
    if(OUTPUT_TREE_FORMAT == 0){
      if(tree->io->long_tax_names){
	      strcat(*s_tree,tree->io->long_tax_names[fils->num]);
	      (*pos) += (int)strlen(tree->io->long_tax_names[fils->num]);
	    }else{
	      strcat(*s_tree,fils->name);
	      (*pos) += (int)strlen(fils->name);
	    }		  
    }
    else{
      (*pos) += sprintf(*s_tree+*pos,"%d",fils->num+1);
    }
    if((fils->b) && (fils->b[0]) && (fils->b[0]->l > -1.)){
      if(tree->print_labels){
	      if(fils->b[0]->n_labels < 10)
          For(i,fils->b[0]->n_labels){
          (*pos) += sprintf(*s_tree+*pos,"#%s",fils->b[0]->labels[i]);
	      }else{
	    	  (*pos) += sprintf(*s_tree+*pos,"#%d_labels",fils->b[0]->n_labels);
	      }
	    }
      
      strcat(*s_tree,":");
      (*pos)++;
      
#ifndef MC
      if(!tree->n_root){
	      (*pos) += sprintf(*s_tree+*pos,"%.10f",fils->b[0]->l);
	    }else{
	      if(pere == tree->n_root){
	    	  phydbl root_pos = (fils == tree->n_root->v[0])?(tree->n_root_pos):(1.-tree->n_root_pos);
	    	  (*pos) += sprintf(*s_tree+*pos,"%.10f",tree->e_root->l * root_pos);
	      }else{
          (*pos) += sprintf(*s_tree+*pos,"%.10f",fils->b[0]->l);
	      }
	    }		
#else
      if(!tree->n_root){
	      (*pos) += sprintf(*s_tree+*pos,"%.10f",fils->b[0]->l);
	    }else{
	      (*pos) += sprintf(*s_tree+*pos,"%.10f",tree->rates->cur_l[fils->num]);
	    }
#endif
    }
    strcat(*s_tree,",");
    (*pos)++;
    
    (*available) = (*available) - (*pos - ori_len);
    
    if(*available < 0)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
    
    if(*available < (int)T_MAX_NAME/2)
    {
      (*s_tree) = (char *)mRealloc(*s_tree,*pos+(int)T_MAX_NAME,sizeof(char));
      (*available) = (int)T_MAX_NAME;
    }
    /*       printf(" %s [%d,%d]",*s_tree,(int)strlen(*s_tree),*available); */
  }
  else
  {
    
    (*s_tree)[(*pos)]='(';
    (*s_tree)[(*pos)+1]='\0';
    (*pos)++;
    (*available)--;
    
    if(tree->n_root)
    {
      For(i,3)
	    {
	      if((fils->v[i] != pere) && (fils->b[i] != tree->e_root))
          R_wtree(fils,fils->v[i],available,s_tree,pos,tree);
	      else p=i;
	    }
    }
    else
    {
      For(i,3)
	    {
	      if(fils->v[i] != pere)
          R_wtree(fils,fils->v[i],available,s_tree,pos,tree);
	      else p=i;
	    }
    }
    
    ori_len = *pos;
    
    /*       printf("\n+ Writing on %p",*s_tree); */
    (*s_tree)[(*pos)-1] = ')';
    (*s_tree)[(*pos)]   = '\0';
    
    if((fils->b) && (fils->b[0]->l > -1.))
    {
      if(tree->print_boot_val)
	    {
	      (*pos) += sprintf(*s_tree+*pos,"%d",fils->b[p]->bip_score);
	    }
      else if(tree->print_alrt_val)
	    {
	      (*pos) += sprintf(*s_tree+*pos,"%.10f",fils->b[p]->ratio_test);
	    }
      if(tree->print_labels)
	    {
	      if(fils->b[p]->n_labels < 10)
          For(i,fils->b[p]->n_labels) 
        {
          (*pos) += sprintf(*s_tree+*pos,"#%s",fils->b[p]->labels[i]);
        }
	      else
        {
          (*pos) += sprintf(*s_tree+*pos,"#%d_labels",fils->b[p]->n_labels);
        }
	    }
      
      strcat(*s_tree,":");
      (*pos)++;
      
#ifndef MC
      if(!tree->n_root)
	    {
	      (*pos) += sprintf(*s_tree+*pos,"%.10f",fils->b[p]->l);
	    }
      else
	    {
	      if(pere == tree->n_root)
        {
          phydbl root_pos = (fils == tree->n_root->v[0])?(tree->n_root_pos):(1.-tree->n_root_pos);
          (*pos) += sprintf(*s_tree+*pos,"%.10f",tree->e_root->l * root_pos);
        }
	      else
        {
          (*pos) += sprintf(*s_tree+*pos,"%.10f",fils->b[p]->l);
        }
	    }
#else
      if(!tree->n_root)
	    {
	      (*pos) += sprintf(*s_tree+*pos,"%.10f",fils->b[p]->l);
	    }
      else
	    {
	      (*pos) += sprintf(*s_tree+*pos,"%.10f",tree->rates->cur_l[fils->num]);
	    }
#endif
    }
    strcat(*s_tree,",");
    (*pos)++;
    (*available) = (*available) - (*pos - ori_len);
    
    if(*available < 0)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
    
    if(*available < (int)T_MAX_NAME/2)
    {
      (*s_tree) = (char *)mRealloc(*s_tree,*pos+(int)T_MAX_NAME,sizeof(char));
      (*available) = (int)T_MAX_NAME;
    }
    /*       printf(" %s [%d,%d]",*s_tree,(int)strlen(*s_tree),*available); */
  }
}

/*********************************************************/

void Init_Tree(t_tree *tree, int n_otu)
{
  tree->n_otu                     = n_otu;
  tree->mat                       = NULL;
  tree->n_root                    = NULL;
  tree->e_root                    = NULL;
  tree->ps_tree                   = NULL;
  
  tree->depth_curr_path           = 0;
  tree->has_bip                   = NO;
  tree->n_moves                   = 0;
  tree->n_improvements            = 0;
  tree->bl_from_node_stamps       = 0;
  tree->lock_topo                 = 0;
  tree->ps_page_number            = 0;
  tree->init_lnL                  = UNLIKELY;
  tree->best_lnL                  = UNLIKELY;
  tree->old_lnL                   = UNLIKELY;
  tree->c_lnL                     = UNLIKELY;
  tree->sum_min_sum_scale         = .0;
  tree->n_swap                    = 0;
  tree->best_pars                 = 1E+5;
  
  tree->n_pattern                 = -1;
  tree->n_root_pos                = -1.;
  tree->print_labels              = 1;
  
  tree->print_boot_val            = 0;
  tree->print_alrt_val            = 0;
  tree->num_curr_branch_available = 0;
  
  tree->tip_order_score           = .0;
}

/*********************************************************/
void Make_New_Edge_Label(t_edge *b)
{
  int i;
  
  b->labels = (char **)realloc(b->labels,(b->n_labels+BLOCK_LABELS)*sizeof(char *));
  
  if(!b->labels)
  {
    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    Warn_And_Exit("");
  }
  else
  {
    for(i=b->n_labels;i<b->n_labels+100;i++) b->labels[i] = (char *)mCalloc(T_MAX_LABEL,sizeof(char));
  }
}

/*********************************************************/

t_edge *Make_Edge_Light(t_node *a, t_node *d, int num)
{
  t_edge *b;
  
  b = (t_edge *)mCalloc(1,sizeof(t_edge));
  
  Init_Edge_Light(b,num);
  
  if(a && b)
  {
    b->left = a;  b->rght = d;
    if(a->tax) {b->rght = a; b->left = d;} /* root */
    /* a tip is necessary on the right side of the t_edge */
    
    (b->left == a)?
    (Make_Edge_Dirs(b,a,d)):
    (Make_Edge_Dirs(b,d,a));
    
    b->l                    = a->l[b->l_r];
    if(a->tax) b->l         = a->l[b->r_l];
    if(b->l < BL_MIN)  b->l = BL_MIN;
    else if(b->l > BL_MAX) b->l = BL_MAX;
    b->l_old                = b->l;
  }
  else
  {
    b->left = NULL;
    b->rght = NULL;
  }
  
  return b;
  
}

/*********************************************************/

void Init_Edge_Light(t_edge *b, int num)
{
  b->num                  = num;
  b->bip_score            = 0;
  b->dist_btw_edges       = .0;
  b->topo_dist_btw_edges  = 0;
  b->has_zero_br_len      = 0;
  b->n_jumps              = 0;
  
  b->p_lk_left            = NULL;
  b->p_lk_rght            = NULL;
  b->Pij_rr               = NULL;
}

/*********************************************************/

void Init_Node_Light(t_node *n, int num)
{
  n->num                    = num;
  n->tax                    = -1;
  n->dist_to_root           = .0;
  n->common                 = 1;
  n->ext_node               = NULL;
  n->name                   = NULL;
  n->ori_name               = NULL;
  n->y_rank                 = 0.;
  n->y_rank_ori             = 0.;
  n->y_rank_max             = 0.;
  n->y_rank_min             = 0.;
  n->anc                    = NULL;
}

/*********************************************************/

void Make_Edge_Dirs(t_edge *b, t_node *a, t_node *d)
{
  int i;
  
  if(a == b->rght)
  {
    PhyML_Printf("\n. a->num = %3d ; d->num = %3d",a->num,d->num);
    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    Warn_And_Exit("");
  }
  if(d == b->left)
  {
    PhyML_Printf("\n. a->num = %3d ; d->num = %3d",a->num,d->num);
    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    Warn_And_Exit("");
  }
  
  b->l_r = b->r_l = -1;
  For(i,3)
  {
    if((a->v[i]) && (a->v[i] == d))
    {
      b->l_r  = i; /* we consider here that 'a' is on the left handside of 'b'*/
      a->b[i] = b;
    }
    if((d->v[i]) && (d->v[i] == a))
    {
      b->r_l  = i; /* we consider here that 'd' is on the right handside of 'b'*/
      d->b[i] = b;
    }
  }
  
  if(a->tax) {b->r_l = 0; For(i,3) if(d->v[i]==a) {b->l_r = i; break;}}
  
  b->l_v1 = b->l_v2 = b->r_v1 = b->r_v2 = -1;
  For(i,3)
  {
    if(b->left->v[i] != b->rght)
    {
      if(b->l_v1 < 0) b->l_v1 = i;
      else            b->l_v2 = i;
    }
    
    if(b->rght->v[i] != b->left)
    {
      if(b->r_v1 < 0) b->r_v1 = i;
      else            b->r_v2 = i;
    }
  }
}

/*********************************************************/


void Make_Edge_Pars(t_edge *b, t_tree *tree)
{
  /*   int site; */
  
  b->pars_l = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
  b->pars_r = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
  
  
  b->ui_l = (unsigned int *)mCalloc(tree->data->crunch_len,sizeof(unsigned int));
  b->ui_r = (unsigned int *)mCalloc(tree->data->crunch_len,sizeof(unsigned int));
  
  
  b->p_pars_l = (int *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(int ));
  b->p_pars_r = (int *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(int ));
  
}

/*********************************************************/
void Make_Edge_Lk(t_edge *b, t_tree *tree)
{
  int ns;
  int i;
  ns = -1;
  
  ns = tree->mod->ns;
  
  b->l_old = b->l;
  
  b->div_post_pred_left = (short int *)mCalloc(ns,sizeof(short int));
  b->div_post_pred_rght = (short int *)mCalloc(ns,sizeof(short int));
  
  //Added by Ken 17/8/2016

  b->bPmat_part = (phydbl**)mCalloc(tree->mod->nparts,sizeof(phydbl*));
  b->up_upp = 1;
  for(i=0;i<tree->mod->nparts;i++){
	phydbl* tp = (phydbl *)mCalloc(tree->mod->n_catg*tree->mod->ns*tree->mod->ns,sizeof(phydbl));
	b->bPmat_part[i]=tp;
  }

  b->upp = (phydbl**)malloc(tree->n_pattern*sizeof(phydbl*));
  for(i=0;i<tree->n_pattern;i++){
	  b->upp[i] = (phydbl *)malloc(tree->mod->ns*sizeof(phydbl));
  }
  
  b->sum_scale_left_cat = (int *)mCalloc(tree->mod->n_catg,sizeof(int));  
  b->sum_scale_rght_cat = (int *)mCalloc(tree->mod->n_catg,sizeof(int));
  b->sum_scale_upp_cat = (int *)mCalloc(tree->mod->n_catg,sizeof(int));
  b->sum_scale_upp = (int *)mCalloc(tree->data->crunch_len*tree->mod->n_catg,sizeof(int));

  if(!b->left->tax)
    b->sum_scale_left = (int *)mCalloc(tree->data->crunch_len*tree->mod->n_catg,sizeof(int));
  else
    b->sum_scale_left = NULL;
  
  if(!b->rght->tax)
    b->sum_scale_rght = (int *)mCalloc(tree->data->crunch_len*tree->mod->n_catg,sizeof(int));
  else
    b->sum_scale_rght = NULL;
  

  
  if((!b->left->tax) || (tree->mod->s_opt->greedy))
  {
    b->p_lk_left = (phydbl *)mCalloc(tree->data->crunch_len*tree->mod->n_catg*tree->mod->ns,sizeof(phydbl));
    b->p_lk_tip_l = NULL;
  }
  else if(b->left->tax)
  {
    b->p_lk_left   = NULL;      
    b->p_lk_tip_l  = (short int *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(short int ));
  }  
  
  if((!b->rght->tax) || (tree->mod->s_opt->greedy))
  {
    b->p_lk_rght = (phydbl *)mCalloc(tree->data->crunch_len*tree->mod->n_catg*tree->mod->ns,sizeof(phydbl));
    b->p_lk_tip_r = NULL;
  }
  else if(b->rght->tax)
  {
    b->p_lk_rght = NULL;      
    b->p_lk_tip_r  = (short int *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(short int));
  }
}
/*********************************************************/

void Make_Edge_NNI(t_edge *b)
{
  b->nni    = Make_NNI();
  b->nni->b = b;
  b->nni->left = b->left;
  b->nni->rght = b->rght;
}

/*********************************************************/

nni *Make_NNI()
{
  nni *a_nni;
  a_nni = (nni *)mCalloc(1,sizeof(nni ));
  Init_NNI(a_nni);
  return a_nni;
}

/*********************************************************/

void Init_NNI(nni *a_nni)
{
  a_nni->left         = NULL;
  a_nni->rght         = NULL;
  a_nni->b            = NULL;
  a_nni->init_l       = -1.;
  a_nni->init_lk      = .0;
  a_nni->score        = +1.0;
  a_nni->best_l       = -1.;
  a_nni->swap_node_v1 = NULL;
  a_nni->swap_node_v2 = NULL;
  a_nni->swap_node_v3 = NULL;
  a_nni->swap_node_v4 = NULL;
  a_nni->lk0          = UNLIKELY;
  a_nni->lk1          = UNLIKELY;
  a_nni->lk2          = UNLIKELY;
  a_nni->l0           = -1.0;
  a_nni->l1           = -1.0;
  a_nni->l2           = -1.0;
}

/*********************************************************/

t_node *Make_Node_Light(int num)
{
  t_node *n;
  n        = (t_node *)mCalloc(1,sizeof(t_node));
  n->v     = (t_node **)mCalloc(3,sizeof(t_node *));
  n->b     = (t_edge **)mCalloc(3,sizeof(t_edge *));
  n->l     = (phydbl *)mCalloc(3,sizeof(phydbl));
  n->score = (phydbl *)mCalloc(3,sizeof(phydbl));
  
  Init_Node_Light(n,num);
  
  n->name  = (char *)mCalloc(T_MAX_NAME,sizeof(char)); //!< Added by Marcelo.
  n->ori_name  = (char *)mCalloc(T_MAX_NAME,sizeof(char)); //!< Added by Marcelo.
  return n;
}

/*********************************************************/


void Make_Node_Lk(t_node *n)
{
  /*   n->n_ex_nodes = (int *)mCalloc(2,sizeof(int)); */
  return;
}

/**********************************************************/
// Scan across FASTA file to determine number and length of sequences
void scanFasta(model *mod, FILE* f){
    rewind(f);

  int nseq=0;
  int seql=0;
  char c;
  int maxnl=-1;
  int firstseql=-1;
  int seq=0;
    int id=0;
    int otu=-1;
    int idpos=0;
    do{
      c=(char)fgetc(mod->io->fp_in_align);
      if(c == '>' || c == EOF){
        id=1;
        otu++;
        idpos=0;
        if(otu==1)firstseql=seql;
        if(firstseql != seql && otu > 0){
          printf("\nSome sequences in fasta file are not the same length! %d vs %d\n",firstseql,seql);
          exit(EXIT_FAILURE);
        }
        if(idpos > T_MAX_NAME){
          printf("An ID is too long. Maximum sequence name is %d\n",T_MAX_NAME);
          exit(EXIT_FAILURE);
        }
        seql=0;
      }else if(c == '\n' && id){
        id=0;
      }
      if(c == EOF)break;

      if(id){
        if(c != '>') idpos++;
      }else if(c!='\n')seql++;
    }while(1);

  mod->io->n_otu=otu;
  mod->io->init_len=firstseql;
  rewind(f);
  //printf("taxa and length %d\t%d\n",mod->io->n_otu,mod->io->init_len);
}


/**********************************************************/

align **Read_Seq_Fasta(option *io, model *mod)
{
  int i,pos;
  char *line;
  align **data;
  char c;
  data   = (align **)mCalloc(mod->io->n_otu,sizeof(align *));
  For(i,mod->io->n_otu){
     data[i]        = (align *)mCalloc(1,sizeof(align));
     data[i]->name  = (char *)mCalloc(T_MAX_NAME,sizeof(char));
     data[i]->state = (char *)mCalloc(mod->io->init_len*mod->state_len+1,sizeof(char));
     data[i]->is_ambigu = NULL;
     data[i]->len = 0;
   }

  int seq=0;
  int id=0;
  int otu=-1;
  int idpos=0;
  do{
    c=(char)fgetc(mod->io->fp_in_align);
    if(c == '>'){
      id=1;
      otu++;
      idpos=0;
    }else if(c == '\n' && id){
      id=0;
    }else if(c == EOF)break;

    if(id){
      if(c != '>'){
        data[otu]->name[idpos]=c;
        idpos++;
      }
    }else if(c!='\n'){
      Uppercase(&c);
      data[otu]->state[data[otu]->len++]=c;
    }
  }while(1);

  return data;
}

/*********************************************************/

void Detect_Align_File_Format(option *io)
{
  int c;char *m;
  fpos_t curr_pos;
  model* mod=io->mod;
  fgetpos(mod->io->fp_in_align,&curr_pos);
  
  errno = 0;
  
  while((c=fgetc(mod->io->fp_in_align)) != EOF)
  {
    if(errno) io->data_file_format = PHYLIP;
    else if(c == '>'){
      io->data_file_format = FASTA;
      if(!io->quiet)printf("\nLooks like a fasta file..");
      return;
    }
    else if(c == '#')
    {
      char s[10],t[6]="NEXUS";
      m=fgets(s,6,mod->io->fp_in_align);
      if(!strcmp(t,s)) 
      {
        fsetpos(mod->io->fp_in_align,&curr_pos);
        io->data_file_format = NEXUS;
        return;
      }
    }
  }
  
  fsetpos(mod->io->fp_in_align,&curr_pos);
}


/*********************************************************/

align **Get_Seq(option *io)
{
  io->data = NULL;
  
  Detect_Align_File_Format(io);
  
  switch(io->data_file_format)
  {
    case PHYLIP: 
    {
      io->data = Get_Seq_Phylip(io);
      break;
    }
    case FASTA:
    {
      scanFasta(io->mod,io->fp_in_align);
      io->data=Read_Seq_Fasta(io,io->mod);
      break;
    }
    case NEXUS:
    {
        Warn_And_Exit("\n. Err: NEXUS cannot be used. Please switch to Phylip format.\n");
//      io->nex_com_list = Make_Nexus_Com();
//      Init_Nexus_Format(io->nex_com_list);
//      Get_Nexus_Data(io->fp_in_align,io);
//      Free_Nexus(io);
//      break;
    }
    default:
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
      break;
    }
  }
  if(io->n_otu < 3) {
      Warn_And_Exit("\n. Err: we need at least 3 sequences to analyze.\n");
  }
  
  if(!io->data)
  {
    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    Warn_And_Exit("");
  }
  else 
  {  
    int i,j;
    char **buff;
    int *remove;
    int n_unkn,n_removed,pos;
    
    buff = (char **)mCalloc(io->n_otu,sizeof(char *));
    For(i,io->n_otu) buff[i] = (char *)mCalloc(io->data[0]->len,sizeof(char));
    remove = (int *)mCalloc(io->data[0]->len,sizeof(int));
    
    n_removed = 0;
    
    For(i,io->data[0]->len)
    {
      For(j,io->n_otu)
      {
        if((io->data[j]->state[i] == '?') || (io->data[j]->state[i] == '-')) io->data[j]->state[i] = 'X';
        if((io->datatype == NT) && (io->data[j]->state[i] == 'N')) io->data[j]->state[i] = 'X';
        if(io->data[j]->state[i] == 'U') io->data[j]->state[i] = 'T';
      }
      
      n_unkn = 0;
      For(j,io->n_otu) if(io->data[j]->state[i] == 'X') n_unkn++;
      
      if(n_unkn == io->n_otu)
      {
        remove[i] = 1;
        n_removed++;
      }
      
      For(j,io->n_otu) buff[j][i] = io->data[j]->state[i];
    }
    
    pos = 0;
    For(i,io->data[0]->len)
    {
      /* 	  if(!remove[i]) */
      /* 	    { */
      For(j,io->n_otu) io->data[j]->state[pos] = buff[j][i];
      pos++;
      /* 	    } */
    }
    
    For(i,io->n_otu) io->data[i]->len = pos;
    For(i,io->n_otu) free(buff[i]);
    free(buff);
    free(remove);
  }
  
  
  return io->data;
}

/* /\*********************************************************\/ */
/*********************************************************/
/*********************************************************/

align **Get_Seq_Phylip(option *io)
{
  Read_Ntax_Len_Phylip(io->fp_in_align,&io->n_otu,&io->init_len);
  
  if(io->interleaved) io->data = Read_Seq_Interleaved(io);
  else                io->data = Read_Seq_Sequential(io);
  
  return io->data;
}

/*********************************************************/
void Read_Ntax_Len_Phylip(FILE *fp ,int *n_otu, int *n_tax)
{
  char *line;
  int readok;
  
  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));  
  
  readok = 0;
  do
  {
    if(fscanf(fp,"%s",line) == EOF)
    {
      free(line); 
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
    else
    {
      if(strcmp(line,"\n") && strcmp(line,"\r") && strcmp(line,"\t"))
	    {
	      sscanf(line,"%d",n_otu);
	      if(*n_otu <= 0) Warn_And_Exit("\n. The number of taxa cannot be negative.\n");
        
	      if(!fscanf(fp,"%s",line)) Exit("\n");
	      sscanf(line,"%d",n_tax);
	      if(*n_tax <= 0) Warn_And_Exit("\n. The sequence length cannot be negative.\n");
	      else readok = 1;
	    }
    }
  }while(!readok);
  
  free(line);
}

/*********************************************************/

align **Read_Seq_Sequential(option *io)
{
  int i;
  char *line;
  align **data;
  /*   char c; */
  char *format;
  
  
  format = (char *)mCalloc(T_MAX_NAME,sizeof(char));
  line   = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  data   = (align **)mCalloc(io->n_otu,sizeof(align *));
  
  /*   while((c=fgetc(in))!='\n'); */
  /*  while(((c=fgetc(io->fp_in_align))!='\n') && (c != ' ') && (c != '\r') && (c != '\t')); */
  
  For(i,io->n_otu)
  {
    data[i]        = (align *)mCalloc(1,sizeof(align));
    data[i]->name  = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    data[i]->state = (char *)mCalloc(io->init_len*io->mod->state_len+1,sizeof(char));
    
    data[i]->is_ambigu = NULL;
    data[i]->len = 0;
    
    sprintf(format, "%%%ds", T_MAX_NAME);
    if(!fscanf(io->fp_in_align,format,data[i]->name)) Exit("\n");
    
    while(data[i]->len < io->init_len * io->mod->state_len) Read_One_Line_Seq(&data,i,io->fp_in_align);
    
    if(data[i]->len != io->init_len * io->mod->state_len)
    {
      PhyML_Printf("\n. Err: Problem with species %s's sequence (check the format).\n",data[i]->name);
      PhyML_Printf("\n. Observed sequence length: %d, expected length: %d\n",data[i]->len, io->init_len * io->mod->state_len);
      Warn_And_Exit("");
    }
  }
  
  For(i,io->n_otu) data[i]->state[data[i]->len] = '\0';
  
  free(format);
  free(line);
  
  return data;
}

/*********************************************************/

align **Read_Seq_Interleaved(option *io)
{
  int i,end,num_block;
  char *line;
  align **data;//!< array of align pointers.
  /*   char c; */
  char *format;
  
  line   = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  format = (char *)mCalloc(T_MAX_NAME, sizeof(char));
  data   = (align **)mCalloc(io->n_otu,sizeof(align *));
  
  /*   while(((c=fgetc(io->fp_in_align))!='\n') && (c != ' ') && (c != '\r') && (c != '\t')); */
  
  end = 0;
  For(i,io->n_otu)
  {
    data[i]        = (align *)mCalloc(1,sizeof(align));
    data[i]->name  = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    data[i]->state = (char *)mCalloc(io->init_len*io->mod->state_len+1,sizeof(char));
    
    data[i]->len       = 0;
    data[i]->is_ambigu = NULL;
    
    sprintf(format, "%%%ds", T_MAX_NAME);
    /*       sprintf(format, "%%%ds", 10); */
    
    if(!fscanf(io->fp_in_align,format,data[i]->name)) Exit("\n");
    
    if(!Read_One_Line_Seq(&data,i,io->fp_in_align))
    {
      end = 1;
      if((i != io->n_otu) && (i != io->n_otu-1))
	    {
	      PhyML_Printf("\n. Err: Problem with species %s's sequence.\n",data[i]->name);
	      PhyML_Printf("\n. Observed sequence length: %d, expected length: %d\n",data[i]->len, io->init_len * io->mod->state_len);
	      Warn_And_Exit("");
	    }
      break;
    }
  }
  
  if(data[0]->len == io->init_len * io->mod->state_len) end = 1;
  
  /*   if(end) printf("\n. finished yet '%c'\n",fgetc(io->fp_in_align)); */
  if(!end)
  {
    
    end = 0;
    
    num_block = 1;
    do
    {
      num_block++;
      
      /* interblock */
      if(!fgets(line,T_MAX_LINE,io->fp_in_align)) break;
      
      if(line[0] != 13 && line[0] != 10)
	    {
	      PhyML_Printf("\n. One or more missing sequences in block %d.\n",num_block-1);
	      Warn_And_Exit("");
	    }
      
      For(i,io->n_otu) if(data[i]->len != io->init_len * io->mod->state_len) break;
      
      if(i == io->n_otu) break;
      
      For(i,io->n_otu)
	    {
	      if(data[i]->len > io->init_len * io->mod->state_len)
        {
          PhyML_Printf("\n. Observed length=%d expected length=%d.\n",data[i]->len,io->init_len * io->mod->state_len);
          PhyML_Printf("\n. Err: Problem with species %s's sequence.\n",data[i]->name);
          Warn_And_Exit("");
        }
	      else if(!Read_One_Line_Seq(&data,i,io->fp_in_align))
        {
          end = 1;
          if((i != io->n_otu) && (i != io->n_otu-1))
          {
            PhyML_Printf("\n. Err: Problem with species %s's sequence.\n",data[i]->name);
            PhyML_Printf("\n. Observed sequence length: %d, expected length: %d.\n",data[i]->len, io->init_len * io->mod->state_len);
            Warn_And_Exit("");
          }
          break;
        }
	    }
    }while(!end);
  }
  
  For(i,io->n_otu) data[i]->state[data[i]->len] = '\0';
  
  For(i,io->n_otu)
  {
    if(data[i]->len != io->init_len * io->mod->state_len)
    {
      PhyML_Printf("\n. Check sequence '%s' length...\n",data[i]->name);
      Warn_And_Exit("");
    }
  }
  
  free(format);
  free(line);
  
  return data;
}

/*********************************************************/
int Read_One_Line_Seq(align ***data, int num_otu, FILE *in)
{
  char c = ' ';
  int nchar = 0;
  
  while(1)
  {
    //       if((c == EOF) || (c == '\n') || (c == '\r')) break;
    
    if((c == 13) || (c == 10)) 
    {
      // 	  PhyML_Printf("[%d %d]\n",c,nchar); fflush(NULL);
      if(!nchar)
	    {
	      c=(char)fgetc(in);
	      continue;
	    }
      else 
	    { 
        // 	      PhyML_Printf("break\n");
	      break; 
	    }
    }
    else if(c == EOF)
    {
      // 	  PhyML_Printf("EOL\n");
      break;
    }
    else if((c == ' ') || (c == '\t') || (c == 32)) 
    {
      // 	  PhyML_Printf("[%d]",c);
      c=(char)fgetc(in); 
      continue;
    }
    
    nchar++;
    Uppercase(&c);
    
    if(c == '.')
    {
      c = (*data)[0]->state[(*data)[num_otu]->len];
      if(!num_otu)
        Warn_And_Exit("\n. Err: Symbol \".\" should not appear in the first sequence\n");
    }
    (*data)[num_otu]->state[(*data)[num_otu]->len]=c;
    (*data)[num_otu]->len++;
    //       PhyML_Printf("%c",c);
    c = (char)fgetc(in);
    if(c == ';') break;
  }
  
  if(c == EOF) return 0;
  else return 1;  
}

/*********************************************************/

void Uppercase(char *ch)
{
  /* nvert ch to upper case -- either ASCII or EBCDIC */
  *ch = isupper((int)*ch) ? *ch : toupper((int)*ch);
}

/*********************************************************/

calign *Compact_Data(align **data, option *io)
{
  calign *cdata_tmp,*cdata;
  int i,j,k,site;
  // int numsenseCodons;
  int n_patt,which_patt,n_invar;
  char **sp_names;
  int n_otu, n_sites;
  pnode *proot;
  int compress;
  int n_ambigu,is_ambigu;
  
  n_otu      = io->n_otu;
  n_patt     = 0;
  which_patt = 0;
  
  sp_names = (char **)mCalloc(n_otu,sizeof(char *));
  For(i,n_otu)
  {
    sp_names[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    strcpy(sp_names[i],data[i]->name);
  }
  
  if(io->datatype==CODON) Nucleotides2Codons(data, io);//!< Added by Marcelo. Translates nucleotides into codons.
  else{
    if(io->datatype==AA && io->convert_NT_to_AA){
	//nothing yet
    	printf("\n\nCan't use nucleotides or AAs\n\n"); //Ken 4/1/2017
    	exit(EXIT_FAILURE);
      }
  }
  
  cdata_tmp = Make_Cseq(n_otu,data[0]->len,io->mod->state_len,data[0]->len,sp_names, io);
  proot     = (pnode *)Create_Pnode(T_MAX_ALPHABET);
  
  For(i,n_otu) free(sp_names[i]);
  free(sp_names);
  
  if(data[0]->len%io->mod->state_len)
  {
    PhyML_Printf("\n. Sequence length is not a multiple of %d\n",io->mod->state_len);
    Warn_And_Exit("");
  }
  
  compress = io->colalias;
  n_ambigu = 0;
  is_ambigu = 0;
  
  if(!io->quiet && !compress) { PhyML_Printf("\n. WARNING: sequences are not compressed !\n");}
  
  Fors(site,data[0]->len,io->mod->state_len)
  {
    if(io->rm_ambigu) //remove columns with ambiguous data if requested
    {
      is_ambigu = 0;
      For(j,n_otu) if(Is_Ambigu(data[j]->state+site,io->datatype,io->mod->state_len)) break;
      if(j != n_otu)
      {
        is_ambigu = 1;
        n_ambigu++;
      }
    }
    
    if(!is_ambigu)
    {
      if(compress)
      {
        which_patt = -1;
        Traverse_Prefix_Tree(site,-1,&which_patt,&n_patt,data,io,proot);
        if(which_patt == n_patt-1) /* New pattern found */
        {
          n_patt--;
          k=n_patt;
        }
        else
        {
          k = n_patt-10;
        }
      }
      else
      {
        k = n_patt;
      }
      
      if(k == n_patt) /* add a new site pattern */
      {
        For(j,n_otu)
        { 
          Copy_One_State(data[j]->state+site, cdata_tmp->c_seq[j]->state+n_patt*io->mod->state_len, io->mod->state_len);
          if(io->datatype==CODON) CopyExtraFieldsCodon(data[j], site, cdata_tmp->c_seq[j],  n_patt);//!< Added by Marcelo.The name explains everything.
        }
        
        For(i,n_otu)
        {
          For(j,n_otu)
          {
            if (io->datatype==CODON)
            {
              if(!(Are_Compatible(cdata_tmp->c_seq[i]->state+n_patt*io->mod->state_len,
                                  cdata_tmp->c_seq[j]->state+n_patt*io->mod->state_len,
                                  io->mod->state_len,
                                  io->datatype,cdata_tmp->c_seq[i]->alternativeCodons[n_patt],cdata_tmp->c_seq[j]->alternativeCodons[n_patt]))) break;
            }else if(!(Are_Compatible(cdata_tmp->c_seq[i]->state+n_patt*io->mod->state_len,
                                      cdata_tmp->c_seq[j]->state+n_patt*io->mod->state_len,
                                      io->mod->state_len,
                                      io->datatype,0,0))) break;
          }
          if(j != n_otu) break;
        }
        
        if((j == n_otu) && (i == n_otu)) /* all characters at that site are compatible with one another the site may be invariant */
        {
          if (io->datatype==CODON) cdata_tmp->invar[n_patt] = Intersect_Site_Alternatives(data, n_otu, site); //!< Added by Marcelo.The name explains everything.
          else 
            For(j,n_otu)
          {
            cdata_tmp->invar[n_patt] = Assign_State(cdata_tmp->c_seq[j]->state+n_patt*io->mod->state_len, io->datatype, io->mod->state_len);
            if(cdata_tmp->invar[n_patt] > -1.) break; /*!< It is not actually (at least one state in the column is ambiguous).*/
          }
        }
        else cdata_tmp->invar[n_patt] = -1;
        
        cdata_tmp->sitepatt[site] = n_patt;
        cdata_tmp->wght[n_patt]  += 1;
        n_patt                   += 1;
      }
      else
      {
        cdata_tmp->sitepatt[site]    = which_patt;
        cdata_tmp->wght[which_patt] += 1;
      }
    }
  }
  
  data[0]->len -= n_ambigu;
  
  cdata_tmp->init_len                   = data[0]->len;
  cdata_tmp->crunch_len                 = n_patt;
  For(i,n_otu) cdata_tmp->c_seq[i]->len = n_patt;
  
  if(!io->quiet) PhyML_Printf("\n. %d patterns found (out of a total of %d sites). \n",n_patt,data[0]->len);
  
  if((io->rm_ambigu) && (n_ambigu)) PhyML_Printf("\n. Removed %d columns of the alignment as they contain ambiguous characters (e.g., gaps) \n",n_ambigu);
  
  /*   For(i,n_otu) */
  /*     { */
  /*       For(j,cdata_tmp->crunch_len) */
  /* 	{ */
  /* 	  printf("%c",cdata_tmp->c_seq[i]->state[j*io->mod->state_len+1]); */
  /* 	} */
  /*       printf("\n"); */
  /*     } */
  
  n_invar=0;
  For(i,cdata_tmp->crunch_len) 
  {
    if(cdata_tmp->invar[i] > -1.) n_invar+=(int)cdata_tmp->wght[i];
    
  }
  
  if(!io->quiet) PhyML_Printf("\n. %d sites without polymorphism (%.2f%c).\n",n_invar,100.*(phydbl)n_invar/data[0]->len,'%');
  
  cdata_tmp->obs_pinvar = (phydbl)n_invar/data[0]->len;
  
  n_sites = 0;
  For(i,cdata_tmp->crunch_len) n_sites += cdata_tmp->wght[i];
  if(n_sites != data[0]->len / io->mod->state_len)
  {
    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    Warn_And_Exit("");
  }
  else if(io->datatype == CODON)                                                                           //!< Added by Marcelo.
  {
    if(!(io->mod->s_opt->user_state_freq)) {
      switch(io->mod->freq_model) {
        case F1XSENSECODONS:                                                                                 //!< Added by Marcelo.	
        case F1X4:
        case F3X4: 
        case CF3X4: Get_Base_Freqs_CODONS_FaXb(cdata_tmp, data, io->mod->freq_model, io->mod);                //!< Added by Marcelo.
          break;
        case FMODEL:
          break;
        default:  
          PhyML_Printf("\n. Frequency model not implemented.\n",__FILE__,__LINE__);
          Warn_And_Exit("");
          break;
      }
    }
    Free_Prefix_Tree(proot,T_MAX_ALPHABET);
    return cdata_tmp;
  }
  else {/* Uniform state frequency distribution.*/}
  
  cdata = Copy_Cseq(cdata_tmp,io);                                                                       
  
  Free_Cseq(cdata_tmp);                                                                                   
  Free_Prefix_Tree(proot,T_MAX_ALPHABET);
  return cdata;
}
/*********************************************************/
calign *Compact_Cdata(calign *data, option *io)
{
  calign *cdata;
  int i,j,k,site;
  int n_patt,which_patt;
  int n_otu;
  
  n_otu = data->n_otu;
  
  cdata         = (calign *)mCalloc(1,sizeof(calign));
  cdata->n_otu  = n_otu;
  cdata->c_seq  = (align **)mCalloc(n_otu,sizeof(align *));
  cdata->wght   = (int *)mCalloc(data->crunch_len,sizeof(int));
  cdata->b_frq  = (phydbl *)mCalloc(io->mod->ns,sizeof(phydbl));
  cdata->ambigu = (short int *)mCalloc(data->crunch_len,sizeof(short int));
  cdata->invar  = (short int *)mCalloc(data->crunch_len,sizeof(short int));
  
  cdata->crunch_len = cdata->init_len = -1;
  For(j,n_otu)
  {
    cdata->c_seq[j]            = (align *)mCalloc(1,sizeof(align));
    cdata->c_seq[j]->name      = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    strcpy(cdata->c_seq[j]->name,data->c_seq[j]->name);
    cdata->c_seq[j]->state     = (char *)mCalloc(data->crunch_len,sizeof(char));
    cdata->c_seq[j]->is_ambigu = (short int *)mCalloc(data->crunch_len,sizeof(short int));
    cdata->c_seq[j]->state[0]  = data->c_seq[j]->state[0];
    //cdata->c_seq[j]->state  = data->c_seq[j]->state;
    if(io->datatype==CODON)//! Added by Marcelo.
    {
      cdata->c_seq[j]->alternativeCodons=(char **)mCalloc(data->crunch_len,sizeof(char*));
    }
  }
  
  
  n_patt = which_patt =  0;
  
  For(site,data->crunch_len)
  {
    if(data->wght[site])
    {
      For(k,n_patt)
	    {
	      For(j,n_otu)
        {
          if(io->datatype==CODON) //! Added by Marcelo.
          {
            if(cdata->c_seq[j]->state+k != data->c_seq[j]->state+site) break;
          }
          else
          {
            if(strncmp(cdata->c_seq[j]->state+k*io->mod->state_len,
                       data->c_seq[j]->state+site*io->mod->state_len,
                       io->mod->state_len))
              break;
          }
        }
	      
	      if(j == n_otu)
        {
          which_patt = k;
          break;
        }
	    }
      
      /*       /\* TO DO *\/ */
      /*       k = n_patt; */
      
      if(k == n_patt)
	    {
	      For(j,n_otu)
	      { 
          Copy_One_State(data->c_seq[j]->state+site*io->mod->state_len,
                         cdata->c_seq[j]->state+n_patt*io->mod->state_len,
                         io->mod->state_len);
          if(io->datatype==CODON) CopyExtraFieldsCodon(data->c_seq[j], site, cdata->c_seq[j],  n_patt);//!< Added by Marcelo.The name explains everything.
	      }
	      For(i,n_otu)
	      {
          For(j,n_otu)
          {
            if(io->datatype==CODON)
            {
              if(!(Are_Compatible(cdata->c_seq[i]->state+n_patt*io->mod->state_len,
                                  cdata->c_seq[j]->state+n_patt*io->mod->state_len,
                                  io->mod->state_len,
                                  io->datatype,cdata->c_seq[i]->alternativeCodons[n_patt], cdata->c_seq[j]->alternativeCodons[n_patt]))) break;
            }
            else if(!(Are_Compatible(cdata->c_seq[i]->state+n_patt*io->mod->state_len,
                                     cdata->c_seq[j]->state+n_patt*io->mod->state_len,
                                     io->mod->state_len,
                                     io->datatype,0, 0))) break;
          }
          if(j != n_otu) break;
	      }
	      
	      if((j == n_otu) && (i == n_otu)) 
	      {
          if (io->datatype==CODON) cdata->invar[n_patt] = Intersect_Site_Alternatives(data->c_seq, n_otu, site); //!< Added by Marcelo.The name explains everything.
          else 
            For(j,n_otu)
          {
            cdata->invar[n_patt] = Assign_State(cdata->c_seq[j]->state+n_patt*io->mod->state_len,
                                                io->datatype,
                                                io->mod->state_len);
            if(cdata->invar[n_patt] > -1.) break;
          }
	      }
	      else cdata->invar[n_patt] = -1;
	      
	      cdata->wght[n_patt] += data->wght[site];
	      n_patt+=1;
	    }
      else cdata->wght[which_patt] += data->wght[site];
      
      /*       Print_Site(cdata,k,n_otu,"\n",io->stepsize); */
    }
  }
  
  cdata->init_len   = data->crunch_len;
  cdata->crunch_len = n_patt;
  For(i,n_otu) cdata->c_seq[i]->len = n_patt;
  
  
  return cdata;
}

/*********************************************************/
/*********************************************************/

void Traverse_Prefix_Tree(int site, int seqnum, int *patt_num, int *n_patt, align **data, option *io, pnode *n)
{
  int ret_val;
  
  ret_val = -1;
  
  if(seqnum == io->n_otu-1)
  {
    n->weight++;
    if(n->weight == 1)
    {
      n->num = *n_patt;
      (*n_patt) += 1;
    }
    (*patt_num) = n->num;
    return;
  }
  else
  {
    int next_state;
    
    next_state = -1;
    next_state = Assign_State_With_Ambiguity(data[seqnum+1]->state+site, io->datatype, io->mod->state_len);
    
    if(!n->next[next_state]) n->next[next_state] = Create_Pnode(T_MAX_ALPHABET);
    Traverse_Prefix_Tree(site,seqnum+1,patt_num,n_patt,data,io,n->next[next_state]);
  }
}

/*********************************************************/

pnode *Create_Pnode(int size)
{
  pnode *n;
  int i;
  
  n = (pnode *)mCalloc(1,sizeof(pnode ));
  n->next = (pnode **)mCalloc(size,sizeof(pnode *));
  For(i,size) n->next[i] = NULL;
  n->weight = 0;
  n->num = -1;
  return n;
}
/*********************************************************/



/*********************************************************/
/*********************************************************/

t_tree *Read_Tree_File_Phylip(FILE *fp_input_tree)
{
  char *line;
  t_tree *tree;
  int i;
  char c;
  
  
  do
  {
    c=fgetc(fp_input_tree);
  }
  while((c != '(') && (c != EOF));
  
  if(c==EOF) return NULL;
  
  line = (char *)mCalloc(1,sizeof(char));
  
  i=0;
  for(;;)
  {
    if((c == ' ') || (c == '\n'))
    {
      c=fgetc(fp_input_tree);
      if(c == EOF || c == ';') break;
      else continue;
    }
    
    if(c == '[')
    {
      Skip_Comment(fp_input_tree);
      c = fgetc(fp_input_tree);
      if(c == EOF || c == ';') break;
    }
    
    line = (char *)mRealloc(line,i+2,sizeof(char));
    
    line[i]=c;
    i++;
    c=fgetc(fp_input_tree);
    if(c==EOF || c==';') break;
  }
  line[i] = '\0';
  
  tree = Read_Tree(line);
  free(line);
  return tree;
}

/*********************************************************/

void Connect_Edges_To_Nodes_Recur(t_node *a, t_node *d, t_tree *tree)
{
  int i;
  
  Connect_One_Edge_To_Two_Nodes(a,d,tree->t_edges[tree->num_curr_branch_available],tree);
  tree->num_curr_branch_available += 1;
  
  if(d->tax) return;
  else For(i,3) if(d->v[i] != a) Connect_Edges_To_Nodes_Recur(d,d->v[i],tree);
}

/*********************************************************/

void Connect_One_Edge_To_Two_Nodes(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  int i,dir_a_d;
  
  dir_a_d = -1;
  For(i,3) if(a->v[i] == d) {dir_a_d = i; break;}
  
  if(dir_a_d == -1)
  {
    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    Warn_And_Exit("");
  }
  
  a->b[dir_a_d] = b;
  b->num        = tree->num_curr_branch_available;
  b->left       = a;
  b->rght       = d;
  if(a->tax) {b->rght = a; b->left = d;} /* root */
  /* a tip is necessary on the right hand side of the t_edge */
  
  (b->left == a)?
  (Make_Edge_Dirs(b,a,d)):
  (Make_Edge_Dirs(b,d,a));
  
  b->l                    = a->l[b->l_r];
  if(a->tax) b->l         = a->l[b->r_l];
  if(b->l < BL_MIN)  b->l = BL_MIN;
  else if(b->l > BL_MAX) b->l = BL_MAX;
  b->l_old                = b->l;
}

/*********************************************************/

void Update_Dirs(t_tree *tree)
{
  int i;
  int buff;
  t_edge *b;
  
  b = NULL;
  buff = -1;
  For(i,2*tree->n_otu-3)
  {
    b = tree->t_edges[i];
    
    if((!b->left->tax) && (b->left->v[b->l_v1]->num < b->left->v[b->l_v2]->num))
    {
      buff    = b->l_v1;
      b->l_v1 = b->l_v2;
      b->l_v2 = buff;
    }
    if((!b->rght->tax) && (b->rght->v[b->r_v1]->num < b->rght->v[b->r_v2]->num))
    {
      buff    = b->r_v1;
      b->r_v1 = b->r_v2;
      b->r_v2 = buff;
    }
  }
  
}

/*********************************************************/

void Exit(char *message)
{
  fflush(NULL);
  PhyML_Fprintf(stderr,"%s",message);
  exit(1);
}

/*********************************************************/

void *mCalloc(int nb, size_t size)
{
  void *allocated;
  
  if((allocated = calloc((size_t)nb,size)) != NULL)
  /*   if((allocated = malloc((size_t)nb*(size_t)size)) != NULL) */
  {
    return allocated;
  }
  else
    Warn_And_Exit("\n. Err: low memory\n");
  
  return NULL;
}

/*********************************************************/

void *mRealloc(void *p,int nb, size_t size)
{
  if((p = realloc(p,(size_t)nb*size)) != NULL)
    return p;
  else
    Warn_And_Exit("\n. Err: low memory\n");
  
  return NULL;
}


/*********************************************************/

void Qksort_Matrix(phydbl **A, int col, int ilo, int ihi)
{
  phydbl pivot;	// pivot value for partitioning array
  int ulo, uhi;	// indices at ends of unpartitioned region
  int ieq;		// least index of array entry with value equal to pivot
  phydbl *tempEntry;	// temporary entry used for swapping
  
  tempEntry = NULL;
  
  if (ilo >= ihi) {
    return;
  }
  // Select a pivot value.
  pivot = A[(ilo + ihi)/2][col];
  // Initialize ends of unpartitioned region and least index of entry
  // with value equal to pivot.
  ieq = ulo = ilo;
  uhi = ihi;
  // While the unpartitioned region is not empty, try to reduce its size.
  while (ulo <= uhi) {
    if (A[uhi][col] > pivot) {
	    // Here, we can reduce the size of the unpartitioned region and
	    // try again.
	    uhi--;
    } else {
	    // Here, A[uhi] <= pivot, so swap entries at indices ulo and
	    // uhi.
	    tempEntry = A[ulo];
	    A[ulo] = A[uhi];
	    A[uhi] = tempEntry;
	    // After the swap, A[ulo] <= pivot.
	    if (A[ulo][col] < pivot) {
        // Swap entries at indices ieq and ulo.
        tempEntry = A[ieq];
        A[ieq] = A[ulo];
        A[ulo] = tempEntry;
        // After the swap, A[ieq] < pivot, so we need to change
        // ieq.
        ieq++;
        // We also need to change ulo, but we also need to do
        // that when A[ulo] = pivot, so we do it after this if
        // statement.
	    }
	    // Once again, we can reduce the size of the unpartitioned
	    // region and try again.
	    ulo++;
    }
  }
  // Now, all entries from index ilo to ieq - 1 are less than the pivot
  // and all entries from index uhi to ihi + 1 are greater than the
  // pivot.  So we have two regions of the array that can be sorted
  // recursively to put all of the entries in order.
  Qksort_Matrix(A, col, ilo, ieq - 1);
  Qksort_Matrix(A, col, uhi + 1, ihi);
}

/********************************************************/

void Print_Site_Lk(t_tree *tree, FILE *fp)
{
    int site;
    int catg;
    int el;
    phydbl invar;
    phydbl postmean;
    
    token *rootTkn = Emit_Root_Token();
    token *currTkn = rootTkn;
    
    if(!tree->io->print_site_lnl) {
        PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
        Warn_And_Exit("");
    }
    
    if(!tree->io->print_trace) {
        currTkn = Emit_Out_Token(currTkn, "Header", "Header", COMMENTTKN, TFSTRING, "%s", "Note : P(D|M) is the probability of site D given the model M (i.e., the site likelihood)");
        if(tree->mod->n_catg > 1 || tree->mod->invar) {
            currTkn = Emit_Out_Token(currTkn, "Header", "Header", COMMENTTKN, TFSTRING, "%s", "P(D|M,rr[x]) is the probability of site D given the model M and the relative rate");
            currTkn = Emit_Out_Token(currTkn, "Header", "Header", COMMENTTKN, TFSTRING, "%s", "of evolution rr[x], where x is the class of rate to be considered.");
            currTkn = Emit_Out_Token(currTkn, "Header", "Header", COMMENTTKN, TFSTRING, "%s", "We have P(D|M) = \\sum_x P(x) x P(D|M,rr[x]).");
        }
        
        currTkn = Emit_Out_Token(currTkn, "TblHeader", "TblHeader", COMMENTTKN, TFSTRING, "%s", "Fields are as follows:");
        currTkn = Emit_Out_Token(currTkn, "TblHeader", "TblHeader", COMMENTTKN, TFSTRING, "%s", "Site, P(D|M)");
        
        if(tree->mod->n_catg > 1) {
            For(catg,tree->mod->n_catg) {
                char *s;
                int temp = asprintf(&s,"P(D|M,rr[%d]=%5.4f)",catg+1,tree->mod->gamma_rr[catg]);//changed by Ken 12/1/2017
                currTkn = Emit_Out_Token(currTkn, "TblHeader", "TblHeader", COMMENTTKN, TFSTRING, "%s", s);
            }
            
            currTkn = Emit_Out_Token(currTkn, "TblHeader", "TblHeader", COMMENTTKN, TFSTRING, "%s", "Posterior Mean");
        }
        
        
        if(tree->mod->invar) {
            currTkn = Emit_Out_Token(currTkn, "TblHeader", "TblHeader", COMMENTTKN, TFSTRING, "%s", "P(D|M,rr[0]=0)");
        }
        
        currTkn = Emit_Out_Token(currTkn, "Sites", "Sites", SETSTARTTKN, TFSTRING, "%s", "");
        For(site,tree->data->init_len) {
            char *t;
            int temp2 = asprintf(&t, "%i", site+1); //changed by Ken 12/1/2017
            currTkn = Emit_Out_Token(currTkn, t, t, SETSTARTTKN, TFNUMERIC, "%s", "");
            
            el = 2;
            if(tree->mod->n_catg > 1) {
                el += tree->mod->n_catg + 1;
            }
            if(tree->mod->invar) {
                el += 1;
            }
            currTkn = Emit_Out_Token(currTkn, "", "", TBLSTARTTKN, TFNUMERIC, "%i", el);
            
            currTkn = Emit_Out_Token(currTkn, "P(D|M)", "P(D|M)", SCALARTKN, TFNUMERIC, "%-16g", (phydbl)tree->cur_site_lk[tree->data->sitepatt[site]]);
            
            if(tree->mod->n_catg > 1) {
                for(catg = 0; catg < tree->mod->n_catg; catg++) {
                    char *s;
                    int temp = asprintf(&s,"P(D|M,rr[%d]=%5.4f)",catg+1,tree->mod->gamma_rr[catg]); //changed by Ken 12/1/2017
                    currTkn = Emit_Out_Token(currTkn, s, s, SCALARTKN, TFNUMERIC, "%-22g", (phydbl)tree->log_site_lk_cat[catg][tree->data->sitepatt[site]]);
                }
                postmean = .0;
                For(catg, tree->mod->n_catg)
                    postmean += tree->mod->gamma_rr[catg] *
                                tree->log_site_lk_cat[catg][tree->data->sitepatt[site]] *
                                tree->mod->gamma_r_proba[catg];
                postmean /= tree->cur_site_lk[tree->data->sitepatt[site]];
                
                currTkn = Emit_Out_Token(currTkn, "Posterior Mean", "PostMean", SCALARTKN, TFNUMERIC, "%-22g", postmean);
            }
            if(tree->mod->invar) {
                if((phydbl)tree->data->invar[tree->data->sitepatt[site]] > -0.5) {
                    invar = tree->mod->pi[tree->data->invar[tree->data->sitepatt[site]]];
                } else {
                    invar = 0.0;
                }
                currTkn = Emit_Out_Token(currTkn, "P(D|M,rr[0]=0)", "P(D|M,rr[0]=0)", SCALARTKN, TFNUMERIC, "%-16g", invar);
            }
            
            currTkn = Emit_Out_Token(currTkn, "", "", TBLENDTKN, TFNUMERIC, "%i", el);
            currTkn = Emit_Out_Token(currTkn, t, t, SETENDTKN, TFNUMERIC, "%s", "");
        }
    } else {
        For(site,tree->data->init_len)
        currTkn = Emit_Out_Token(currTkn, "", "", SCALARTKN, TFNUMERIC, "%-16g", tree->cur_site_lk[tree->data->sitepatt[site]]);
    }
    currTkn = Emit_Out_Token(currTkn, "Sites", "Sites", SETENDTKN, TFSTRING, "%s", "");
    
    Output_Tokens(rootTkn, tree->io->out_stats_format, fp);
}


/*********************************************************/

void Order_Tree_CSeq(t_tree *tree, calign *cdata)
{
  int i,j,n_otu_tree,n_otu_cdata;
  align *buff;
  
  n_otu_tree  = tree->n_otu;
  n_otu_cdata = cdata->n_otu;
  
  if(n_otu_tree != n_otu_cdata) 
  {
    PhyML_Printf("\n. Number of taxa in the tree: %d, number of sequences: %d.\n",n_otu_tree,n_otu_cdata);
    Warn_And_Exit("\n. The number of tips in the tree is not the same as the number of sequences\n");
  }
  
  For(i,MAX(n_otu_tree,n_otu_cdata))
  {
    For(j,MIN(n_otu_tree,n_otu_cdata))
	  {
	    if(!strcmp(tree->noeud[i]->name,cdata->c_seq[j]->name))
	      break;
	  }
    
    if(j==MIN(n_otu_tree,n_otu_cdata))
	  {
	    PhyML_Printf("\n. Err: %s is not found in sequence data set\n",tree->noeud[i]->name);
	    Warn_And_Exit("");
	  }
    
    buff            = cdata->c_seq[j];
    cdata->c_seq[j] = cdata->c_seq[i];
    cdata->c_seq[i] = buff;
  }
}

/*********************************************************/

matrix *Make_Mat(int n_otu)
{
  matrix *mat;
  int i;
  
  mat = (matrix *)mCalloc(1,sizeof(matrix));
  
  mat->n_otu = n_otu;
  
  mat->P        = (phydbl **)mCalloc(n_otu,sizeof(phydbl *));
  mat->Q        = (phydbl **)mCalloc(n_otu,sizeof(phydbl *));
  mat->dist     = (phydbl **)mCalloc(n_otu,sizeof(phydbl *));
  mat->on_off   = (int *)mCalloc(n_otu,sizeof(int));
  mat->name     = (char **)mCalloc(n_otu,sizeof(char *));
  mat->tip_node = (t_node **)mCalloc(n_otu,sizeof(t_node *));
  
  
  For(i,n_otu)
  {
    mat->P[i]    = (phydbl *)mCalloc(n_otu,sizeof(phydbl));
    mat->Q[i]    = (phydbl *)mCalloc(n_otu,sizeof(phydbl));
    mat->dist[i] = (phydbl *)mCalloc(n_otu,sizeof(phydbl));
    mat->name[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
  }
  
  return mat;
}

/*********************************************************/

void Init_Mat(matrix *mat, calign *data)
{
  int i;
  
  mat->n_otu = data->n_otu;
  mat->r = mat->n_otu;
  mat->curr_int = mat->n_otu;
  mat->method = 1;
  
  For(i,data->n_otu)
  {
    strcpy(mat->name[i],data->c_seq[i]->name);
    mat->on_off[i] = 1;
  }
}

/*********************************************************/

t_tree *Make_Tree_From_Scratch(int n_otu, calign *data)
{
  t_tree *tree;
  
  tree = Make_Tree(n_otu);
  Make_All_Tree_Nodes(tree);
  Make_All_Tree_Edges(tree);
  Make_Tree_Path(tree);
  if(data)
  {
    Copy_Tax_Names_To_Tip_Labels(tree,data);
    tree->data = data;
  }
  return tree;
}

/*********************************************************/
/*********************************************************/

t_tree *Make_Tree(int n_otu)
{
  t_tree *tree;
  // int i;
  tree = (t_tree *)mCalloc(1,sizeof(t_tree ));
  Init_Tree(tree,n_otu);
  tree->t_dir = (short int *)mCalloc((2*n_otu-2)*(2*n_otu-2),sizeof(int));
  return tree;
}

/*********************************************************/

void Make_Tree_Path(t_tree *tree)
{
  tree->curr_path = (t_node **)mCalloc(tree->n_otu,sizeof(t_node *));
}

/*********************************************************/

void Make_All_Tree_Nodes(t_tree *tree)
{
  int i;
  
  tree->noeud = (t_node **)mCalloc(2*tree->n_otu-1,sizeof(t_node *));
  
  For(i,2*tree->n_otu-1)
  {
    tree->noeud[i] = (t_node *)Make_Node_Light(i);
    if(i < tree->n_otu) tree->noeud[i]->tax = 1;
    else                tree->noeud[i]->tax = 0;
  }
}

/*********************************************************/

void Make_All_Tree_Edges(t_tree *tree)
{
  int i;
  tree->t_edges      = (t_edge **)mCalloc(2*tree->n_otu-2,sizeof(t_edge *));
  For(i,2*tree->n_otu-2) tree->t_edges[i] = (t_edge *)Make_Edge_Light(NULL,NULL,i);
}

/*********************************************************/

void Copy_Tax_Names_To_Tip_Labels(t_tree *tree, calign *data)
{
  int i;
  
  For(i,tree->n_otu)
  {
    /*tree->noeud[i]->name = (char *)mCalloc((int)strlen(data->c_seq[i]->name)+1,sizeof(char));*///!<Commented by Marcelo.
    strcpy(tree->noeud[i]->ori_name, tree->noeud[i]->name);//!<Changed by Marcelo.
    strcpy(tree->noeud[i]->name,data->c_seq[i]->name);
    tree->noeud[i]->tax = 1;
    tree->noeud[i]->num = i;
  }
}

/*********************************************************/

void Share_Lk_Struct(t_tree *t_full, t_tree *t_empt)
{
  int i,j,n_otu;
  t_edge *b_e,*b_f;
  t_node *n_e, *n_f;
  
  n_otu                   = t_full->n_otu;
  t_empt->n_root          = t_full->n_root;
  t_empt->e_root          = t_full->e_root;
  t_empt->c_lnL_sorted    = t_full->c_lnL_sorted;
  t_empt->log_site_lk_cat = t_full->log_site_lk_cat;
  t_empt->cur_site_lk     = t_full->cur_site_lk;
  t_empt->old_site_lk     = t_full->old_site_lk;
  t_empt->triplet_struct  = t_full->triplet_struct;
  t_empt->log_lks_aLRT    = t_full->log_lks_aLRT;
  t_empt->site_lk_cat     = t_full->site_lk_cat;
  
  For(i,2*n_otu-3)
  {
    b_f = t_full->t_edges[i];
    b_e = t_empt->t_edges[i];
    
    b_e->Pij_rr = b_f->Pij_rr;
    
    b_e->nni = b_f->nni;
  }
  
  
  for(i=n_otu;i<2*n_otu-2;i++)
  {
    n_f = t_full->noeud[i];
    n_e = t_empt->noeud[i];
    
    For(j,3)
    {
      if(n_f->b[j]->left == n_f)
	    {
	      if(n_e->b[j]->left == n_e)
        {
          n_e->b[j]->p_lk_left          = n_f->b[j]->p_lk_left;
          n_e->b[j]->sum_scale_left     = n_f->b[j]->sum_scale_left;
          n_e->b[j]->sum_scale_left_cat = n_f->b[j]->sum_scale_left_cat;
          n_e->b[j]->p_lk_tip_l         = n_f->b[j]->p_lk_tip_l;
        }
	      else
        {
          n_e->b[j]->p_lk_rght          = n_f->b[j]->p_lk_left;
          n_e->b[j]->sum_scale_rght     = n_f->b[j]->sum_scale_left;
          n_e->b[j]->sum_scale_rght_cat = n_f->b[j]->sum_scale_left_cat;
          n_e->b[j]->p_lk_tip_r         = n_f->b[j]->p_lk_tip_l;
        }
	    }
      else
	    {
	      if(n_e->b[j]->rght == n_e)
        {
          n_e->b[j]->p_lk_rght          = n_f->b[j]->p_lk_rght;
          n_e->b[j]->sum_scale_rght     = n_f->b[j]->sum_scale_rght;
          n_e->b[j]->sum_scale_rght_cat = n_f->b[j]->sum_scale_rght_cat;
          n_e->b[j]->p_lk_tip_r         = n_f->b[j]->p_lk_tip_r;
        }
	      else
        {
          n_e->b[j]->p_lk_left          = n_f->b[j]->p_lk_rght;
          n_e->b[j]->sum_scale_left     = n_f->b[j]->sum_scale_rght;
          n_e->b[j]->sum_scale_left_cat = n_f->b[j]->sum_scale_rght_cat;
          n_e->b[j]->p_lk_tip_l         = n_f->b[j]->p_lk_tip_r;
        }
	    }
    }
  }
  
  For(i,n_otu)
  {
    n_f = t_full->noeud[i];
    n_e = t_empt->noeud[i];
    
    if(n_f->b[0]->rght == n_f)
    {
      n_e->b[0]->p_lk_rght          = n_f->b[0]->p_lk_rght;
      n_e->b[0]->sum_scale_rght     = n_f->b[0]->sum_scale_rght;
      n_e->b[0]->sum_scale_rght_cat = n_f->b[0]->sum_scale_rght_cat;
      n_e->b[0]->p_lk_tip_r         = n_f->b[0]->p_lk_tip_r;
    }
    else
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  }
}

/*********************************************************/

void Share_Spr_Struct(t_tree *t_full, t_tree *t_empt)
{
  t_empt->size_spr_list = t_full->size_spr_list;
  t_empt->spr_list      = t_full->spr_list;
  t_empt->best_spr      = t_full->best_spr;
}

/*********************************************************/

void Share_Pars_Struct(t_tree *t_full, t_tree *t_empt)
{
  int i;
  
  t_empt->site_pars = t_full->site_pars;
  t_empt->step_mat  = t_full->step_mat;
  
  For(i,2*t_full->n_otu-3)
  {
    t_empt->t_edges[i]->ui_l     = t_full->t_edges[i]->ui_l;
    t_empt->t_edges[i]->ui_r     = t_full->t_edges[i]->ui_r;
    
    t_empt->t_edges[i]->pars_l   = t_full->t_edges[i]->pars_l;
    t_empt->t_edges[i]->pars_r   = t_full->t_edges[i]->pars_r;
    
    t_empt->t_edges[i]->p_pars_l = t_full->t_edges[i]->p_pars_l;
    t_empt->t_edges[i]->p_pars_r = t_full->t_edges[i]->p_pars_r;
  }
}

/*********************************************************/

int Sort_Edges_NNI_Score(t_tree *tree, t_edge **sorted_edges, int n_elem)
{
  int i,j;
  t_edge *buff;
  
  For(i,n_elem-1)
  {
    for(j=i+1;j<n_elem;j++)
    {
      if(sorted_edges[j]->nni->score  < sorted_edges[i]->nni->score)
	    {
	      buff = sorted_edges[j];
	      sorted_edges[j] = sorted_edges[i];
	      sorted_edges[i] = buff;
	    }
    }
  }
  return 1;
}

/*********************************************************/

void NNI(t_tree *tree, t_edge *b_fcus, int do_swap)
{
  int l_r, r_l, l_v1, l_v2, r_v3, r_v4;
  t_node *v1,*v2,*v3,*v4;
  phydbl lk0, lk1, lk2;
  phydbl lk0_init, lk1_init, lk2_init;
  phydbl bl_init;
  phydbl l0,l1,l2;
  phydbl l_infa, l_infb, l_max;
  /*   phydbl lk_infa, lk_infb, lk_max; */
  phydbl lk_init;
  
  bl_init                = b_fcus->l;
  lk_init                = tree->c_lnL;
  
  b_fcus->nni->init_l    = b_fcus->l;
  b_fcus->nni->init_lk   = tree->c_lnL;;
  
  b_fcus->nni->best_conf = 0;
  b_fcus->nni->score     = +1.0;
  
  lk0 = lk1 = lk2        = UNLIKELY;
  v1 = v2 = v3 = v4      = NULL;
  
  l_r = r_l = l_v1 = l_v2 = r_v3 = r_v4 = -1;
  
  l_r                    = b_fcus->l_r;
  r_l                    = b_fcus->r_l;
  
  v1                     = b_fcus->left->v[b_fcus->l_v1];
  v2                     = b_fcus->left->v[b_fcus->l_v2];
  v3                     = b_fcus->rght->v[b_fcus->r_v1];
  v4                     = b_fcus->rght->v[b_fcus->r_v2];
  
  
  if(v1->num < v2->num)
  {
    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    Warn_And_Exit("");
  }
  if(v3->num < v4->num)
  {
    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    Warn_And_Exit("");
  }
  
  l0 = l1 = l2 = -1.;
  
  /***********/
  Swap(v2,b_fcus->left,b_fcus->rght,v3,tree);
  tree->both_sides = 1;
  
 // int br;
  	//      	  int n_edges=2*tree->n_otu-3;
  	//      	  For(br,n_edges) Update_PMat_At_Given_Edge(tree->t_edges[br],tree);
  	//      Get_UPP(tree->noeud[tree->mod->startnode], tree->noeud[tree->mod->startnode]->v[0], tree);
 /* if(tree->mod->whichrealmodel == HLP17){
	  if(b_fcus->anc_node->num != tree->mod->startnode){
 	  	  Fill_UPP_single(tree,b_fcus);
  	  }else{
	  	  Fill_UPP_root(tree,b_fcus);
  	  }
  }*/
  lk1_init = Update_Lk_At_Given_Edge(b_fcus,tree);
  
  l_infa = 10.*b_fcus->l;
  l_max  = b_fcus->l;
  l_infb = BL_MIN;
  
  if(tree->mod->s_opt->fast_nni){
    Fast_Br_Len(b_fcus,tree,1);
    lk1 = Lk_At_Given_Edge(b_fcus,tree);
  }else{
	 //  int br;
	//      	  int n_edges=2*tree->n_otu-3;
	//      	  For(br,n_edges) Update_PMat_At_Given_Edge(tree->t_edges[br],tree);
/*  if(tree->mod->whichrealmodel == HLP17){
	  if(b_fcus->anc_node->num != tree->mod->startnode){
		  Fill_UPP_single(tree,b_fcus);
	  }else{
		  Fill_UPP_root(tree,b_fcus);
	  }
  }*/

   lk1 = Br_Len_Brent(l_infa,l_max,l_infb,
                       tree->mod->s_opt->min_diff_lk_local,
                       b_fcus,tree,
                       tree->mod->s_opt->brent_it_max,
                       tree->mod->s_opt->quickdirty);
    Update_PMat_At_Given_Edge(b_fcus,tree);

 /*  if(tree->mod->whichrealmodel == HLP17){
    if(b_fcus->anc_node->num != tree->mod->startnode){
   	  Fill_UPP_single(tree,b_fcus);
    }else{
   	  Fill_UPP_root(tree,b_fcus);
    }
   }*/
  }
  
  if(lk1 < lk1_init - tree->mod->s_opt->min_diff_lk_local){
    PhyML_Printf("%f %f %f %G\n",l_infa,l_max,l_infb,b_fcus->l);
    PhyML_Printf("%f -- %f \n",lk1_init,lk1);
    PhyML_Printf("\n. Err. in NNI (1)\n");
  }
  
  l1  = b_fcus->l;
  Swap(v3,b_fcus->left,b_fcus->rght,v2,tree);
  /***********/
  
  
  /***********/
  Swap(v2,b_fcus->left,b_fcus->rght,v4,tree);
  b_fcus->l = bl_init;
  tree->both_sides = 1;
  	 // n_edges=2*tree->n_otu-3;
  	 //     	  For(br,n_edges) Update_PMat_At_Given_Edge(tree->t_edges[br],tree);
  	 //     Get_UPP(tree->noeud[tree->mod->startnode], tree->noeud[tree->mod->startnode]->v[0], tree);

  /*if(tree->mod->whichrealmodel == HLP17){
	  if(b_fcus->anc_node->num != tree->mod->startnode){
		  Fill_UPP_single(tree,b_fcus);
	  }else{
		  Fill_UPP_root(tree,b_fcus);
	  }
  }*/

  lk2_init = Update_Lk_At_Given_Edge(b_fcus,tree);
  
  l_infa = 10.*b_fcus->l;
  l_max  = b_fcus->l;
  l_infb = BL_MIN;
  
  if(tree->mod->s_opt->fast_nni){
    Fast_Br_Len(b_fcus,tree,1);
    lk2 = Lk_At_Given_Edge(b_fcus,tree);
  }else{
	int br;
	int n_edges=2*tree->n_otu-3;
	 //   For(br,n_edges) Update_PMat_At_Given_Edge(tree->t_edges[br],tree);
	//      Get_UPP(tree->noeud[tree->mod->startnode], tree->noeud[tree->mod->startnode]->v[0], tree);
    
    lk2 = Br_Len_Brent(l_infa,l_max,l_infb,
                       tree->mod->s_opt->min_diff_lk_local,
                       b_fcus,tree,
                       tree->mod->s_opt->brent_it_max,
                       tree->mod->s_opt->quickdirty);
    Update_PMat_At_Given_Edge(b_fcus,tree);
  /*  if(tree->mod->whichrealmodel == HLP17){
    	if(b_fcus->anc_node->num != tree->mod->startnode){
    		Fill_UPP_single(tree,b_fcus);
    	}else{
    		Fill_UPP_root(tree,b_fcus);
    	}
    }*/

  }
  
  if(lk2 < lk2_init - tree->mod->s_opt->min_diff_lk_local)
  {
    PhyML_Printf("%f %f %f %G\n",l_infa,l_max,l_infb,b_fcus->l);
    PhyML_Printf("%f -- %f \n",lk2_init,lk2);
    PhyML_Printf("\n. Err. in NNI (2)\n");
  }
  
  l2  = b_fcus->l;
  Swap(v4,b_fcus->left,b_fcus->rght,v2,tree);
  /***********/
  
  
  
  /***********/
  b_fcus->l = bl_init;
  if(b_fcus->l < BL_MIN) 
  {
    b_fcus->l = BL_MIN;
  }
  
  tree->both_sides = 1;
 // n_edges=2*tree->n_otu-3;
  	//      	  For(br,n_edges) Update_PMat_At_Given_Edge(tree->t_edges[br],tree);
  	//      Get_UPP(tree->noeud[tree->mod->startnode], tree->noeud[tree->mod->startnode]->v[0], tree);
  /*if(tree->mod->whichrealmodel == HLP17){
	  if(b_fcus->anc_node->num != tree->mod->startnode){
		  Fill_UPP_single(tree,b_fcus);
	  }else{
		  Fill_UPP_root(tree,b_fcus);
	  }
  }*/

  lk0_init = Update_Lk_At_Given_Edge(b_fcus,tree);
  
  if(tree->io->datatype==CODON) //!< Added by Marcelo.This was a sanity check ... we can include it again after the likelihood calculation is revised.
  {
    //     if(lk0_init < lk_init)
    //     if(FABS(lk0_init - lk_init) > tree->io->mod->s_opt->min_diff_lk_codonModels)
    //     { 
    //       PhyML_Printf("\n. lk_init = %.20f; lk_curr = %.20f diff = %.20f l = %G\n",
    //       lk_init,
    //       lk0_init,
    //       lk_init-lk0_init,
    //       b_fcus->l);
    //       PhyML_Printf("\n. Curr_lnL = %.20f\n",Lk(tree));
    //       Warn_And_Exit("\n. Err. in NNI (3)\n");
    //     }
  }
  else
  {
    //     if(FABS(lk0_init - lk_init) > tree->mod->s_opt->min_diff_lk_local)
    //     {
    //       PhyML_Printf("\n. lk_init = %.20f; lk = %.20f diff = %.20f l = %G\n",
    // 		    lk_init,
    // 		    lk0_init,
    // 		    lk_init-lk0_init,
    // 		    b_fcus->l);
    // 		    PhyML_Printf("\n. Curr_lnL = %.20f\n",Lk(tree));
    // 		    Warn_And_Exit("\n. Err. in NNI (3)\n");
    //     }
  }
  
  l_infa = 10.*b_fcus->l;
  l_max  = b_fcus->l;
  l_infb = BL_MIN;
  
  if(tree->mod->s_opt->fast_nni)
  {
    Fast_Br_Len(b_fcus,tree,1);
    lk0 = Lk_At_Given_Edge(b_fcus,tree);
  }
  else
  {
	   //int br;
	    //  	  int n_edges=2*tree->n_otu-3;
	   //   	  For(br,n_edges) Update_PMat_At_Given_Edge(tree->t_edges[br],tree);
	   //   Get_UPP(tree->noeud[tree->mod->startnode], tree->noeud[tree->mod->startnode]->v[0], tree);
	/*  if(tree->mod->whichrealmodel == HLP17){
		  if(b_fcus->anc_node->num != tree->mod->startnode){
	       	  Fill_UPP_single(tree,b_fcus);
	      }else{
	      	  Fill_UPP_root(tree,b_fcus);
	      }
	  }*/

    lk0 = Br_Len_Brent(l_infa,l_max,l_infb,
                       tree->mod->s_opt->min_diff_lk_local,
                       b_fcus,tree,
                       tree->mod->s_opt->brent_it_max,
                       tree->mod->s_opt->quickdirty);
    Update_PMat_At_Given_Edge(b_fcus,tree);

    /*if(tree->mod->whichrealmodel == HLP17){
    	if(b_fcus->anc_node->num != tree->mod->startnode){
    		Fill_UPP_single(tree,b_fcus);
    	}else{
    		Fill_UPP_root(tree,b_fcus);
    	}
    }*/
    
  }
  if(tree->io->datatype==CODON) //!< Changed by Marcelo.//!< Added by Marcelo.This was a sanity check ... we can include it again after the likelihood calculation is revised.
  {
    //     if(lk0 < lk_init)
    //     if(lk0 < lk_init - tree->mod->s_opt->min_diff_lk_codonModels)
    //     {
    //       PhyML_Printf("\n\n%f %f %f %f\n",l_infa,l_max,l_infb,b_fcus->l);
    //       PhyML_Printf("%f -- %f \n",lk0_init,lk0);
    //       PhyML_Printf("\n. Err. in NNI (3)\n");
    //       Warn_And_Exit("\n");
    //     }
  }
  else
  {
    //     if(lk0 < lk_init - tree->mod->s_opt->min_diff_lk_local)
    //     {
    //       PhyML_Printf("\n\n%f %f %f %f\n",l_infa,l_max,l_infb,b_fcus->l);
    //       PhyML_Printf("%f -- %f \n",lk0_init,lk0);
    //       PhyML_Printf("\n. Err. in NNI (3)\n");
    //       Warn_And_Exit("\n");
    //     }
  }
  l0  = b_fcus->l;
  /***********/
  
  b_fcus->nni->lk0 = lk0;
  b_fcus->nni->lk1 = lk1;
  b_fcus->nni->lk2 = lk2;
  
  b_fcus->nni->l0  = l0;
  b_fcus->nni->l1  = l1;
  b_fcus->nni->l2  = l2;  
  
  b_fcus->nni->score = lk0 - MAX(lk1,lk2);
  
  if((b_fcus->nni->score <  tree->mod->s_opt->min_diff_lk_local) &&
     (b_fcus->nni->score > -tree->mod->s_opt->min_diff_lk_local))
  {
    b_fcus->nni->score = .0;
    b_fcus->nni->lk1 = b_fcus->nni->lk0;
    b_fcus->nni->lk2 = b_fcus->nni->lk0;
  }
  
  if(lk0 > MAX(lk1,lk2))
  {
    b_fcus->nni->best_conf    = 0;
    b_fcus->nni->best_l       = l0;
    b_fcus->nni->swap_node_v1 = NULL;
    b_fcus->nni->swap_node_v2 = NULL;
    b_fcus->nni->swap_node_v3 = NULL;
    b_fcus->nni->swap_node_v4 = NULL;
  }
  else if(lk1 > MAX(lk0,lk2))
  {
    b_fcus->nni->best_conf    = 1;
    b_fcus->nni->best_l       = l1;
    b_fcus->nni->swap_node_v1 = v2;
    b_fcus->nni->swap_node_v2 = b_fcus->left;
    b_fcus->nni->swap_node_v3 = b_fcus->rght;
    b_fcus->nni->swap_node_v4 = v3;
  }
  else if(lk2 > MAX(lk0,lk1))
  {
    b_fcus->nni->best_conf    = 2;
    b_fcus->nni->best_l       = l2;
    b_fcus->nni->swap_node_v1 = v2;
    b_fcus->nni->swap_node_v2 = b_fcus->left;
    b_fcus->nni->swap_node_v3 = b_fcus->rght;
    b_fcus->nni->swap_node_v4 = v4;
  }
  else
  {
    b_fcus->nni->score        = +1.0;
    b_fcus->nni->best_conf    = 0;
    b_fcus->nni->best_l       = l0;
    b_fcus->nni->swap_node_v1 = NULL;
    b_fcus->nni->swap_node_v2 = NULL;
    b_fcus->nni->swap_node_v3 = NULL;
    b_fcus->nni->swap_node_v4 = NULL;
  }
  
  if((do_swap) && ((lk1 > lk0) || (lk2 > lk0)))
  {
    tree->n_swap++;
    PhyML_Printf("Swap t_edge %d -> %f\n",b_fcus->num,MAX(lk1,lk2));
    
    if(lk1 > lk2)
    {
      tree->best_lnL = lk1;
      Swap(v2,b_fcus->left,b_fcus->rght,v3,tree);
      b_fcus->l = l1;
      tree->both_sides = 1;
      Lk(tree);
    }
    else
    {
      tree->best_lnL = lk2;
      Swap(v2,b_fcus->left,b_fcus->rght,v4,tree);
      b_fcus->l = l2;
      tree->both_sides = 1;
      Lk(tree);
    }
  }
  else
  {
    b_fcus->l = bl_init;
    Update_PMat_At_Given_Edge(b_fcus,tree);
    tree->c_lnL = lk_init;
  }
}

/*********************************************************/

void Swap(t_node *a, t_node *b, t_node *c, t_node *d, t_tree *tree)
{
  int ab, ba, cd, dc, bc;
  int i;
  
  
  /* \             /d      \             /a
   *  \           /         \           /
   *   \b__...__c/    ->     \b__...__c/
   *   /         \	   		 /		   \
   *  /           \	        /	        \
   * /a            \  	   /d            \
   *
   * nodes b and c are not necessarily on the same branch 
   */
  
  
#ifdef DEBUG
  if(!a || !b || !c || !d)
  {
    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    Warn_And_Exit("");
  }
#endif
  
  
  ab = ba = cd = dc = bc = -1;
  
  For(i,3) if(a->v[i] == b) { ab = i; break; }
  For(i,3) if(b->v[i] == a) { ba = i; break; }
  For(i,3) if(c->v[i] == d) { cd = i; break; }
  For(i,3) if(d->v[i] == c) { dc = i; break; }
  For(i,3) if(b->v[i] == c) { bc = i; break; }
  

  //figure out ancestral branches

#ifdef DEBUG
  if(ab < 0 || ba < 0 || cd < 0 || dc < 0)
  {
    PhyML_Printf("\n. Nodes %d %d %d %d\n",a->num,b->num,c->num,d->num);
    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    Warn_And_Exit("");
  }
#endif
  
  a->v[ab] = c; //change connected nodes
  d->v[dc] = b;
  b->v[ba] = d;
  c->v[cd] = a;
  b->b[ba] = d->b[dc]; //change connected edges from b and c
  c->b[cd] = a->b[ab];
  	  	  	  	  	  	  //still updating edges
  (a->b[ab]->left == b)?  //if b was to the left of former edge ab, change c to that edge's left, otherwise change c to its right
  (a->b[ab]->left = c):
  (a->b[ab]->rght = c);
  
  (d->b[dc]->left == c)? //same thing with d
  (d->b[dc]->left = b):
  (d->b[dc]->rght = b);
  
  For(i,3)
  {
    if(a->b[ab]->left->v[i] == a->b[ab]->rght) a->b[ab]->l_r = i; //update left/right directions
    if(a->b[ab]->rght->v[i] == a->b[ab]->left) a->b[ab]->r_l = i;
    if(d->b[dc]->left->v[i] == d->b[dc]->rght) d->b[dc]->l_r = i;
    if(d->b[dc]->rght->v[i] == d->b[dc]->left) d->b[dc]->r_l = i;
  }
  
  
  a->b[ab]->l_v1 = a->b[ab]->l_v2 =
  a->b[ab]->r_v1 = a->b[ab]->r_v2 =
  d->b[dc]->l_v1 = d->b[dc]->l_v2 =
  d->b[dc]->r_v1 = d->b[dc]->r_v2 = -1;
  
  
  For(i,3)
  {
    if(i != a->b[ab]->l_r)
    {
      if(a->b[ab]->l_v1 < 0) a->b[ab]->l_v1 = i;
      else a->b[ab]->l_v2 = i;
    }
    if(i != a->b[ab]->r_l)
    {
      if(a->b[ab]->r_v1 < 0) a->b[ab]->r_v1 = i;
      else a->b[ab]->r_v2 = i;
    }
    if(i != d->b[dc]->l_r)
    {
      if(d->b[dc]->l_v1 < 0) d->b[dc]->l_v1 = i;
      else d->b[dc]->l_v2 = i;
    }
    if(i != d->b[dc]->r_l)
    {
      if(d->b[dc]->r_v1 < 0) d->b[dc]->r_v1 = i;
      else d->b[dc]->r_v2 = i;
    }
  }
  Update_Dirs(tree);
  if(tree->mod->whichrealmodel == HLP17){
	  Update_Ancestors_Edge(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree->noeud[tree->mod->startnode]->b[0],tree); //added by Ken 7/11
  }
 // printf("\n\n");

 /* int br;
  int n_edges=2*tree->n_otu-3;
  For(br,n_edges) Update_PMat_At_Given_Edge(tree->t_edges[br],tree);
  Get_UPP(tree->noeud[tree->mod->startnode], tree->noeud[tree->mod->startnode]->v[0], tree);*/

}

/*********************************************************/

calign *Make_Cseq(int n_otu, int crunch_len, int state_len, int init_len, char **sp_names, option *io)
{
  calign *cdata;
  int j;
  // int i;
  
  cdata           = (calign *)mCalloc(1,sizeof(calign));
  cdata->n_otu    = n_otu;
  cdata->c_seq    = (align **)mCalloc(n_otu,sizeof(align *));
  cdata->b_frq    = (phydbl *)mCalloc(T_MAX_ALPHABET,sizeof(phydbl));
  cdata->wght     = (int *)mCalloc(crunch_len,sizeof(int));
  cdata->ambigu   = (short int *)mCalloc(crunch_len,sizeof(short int));
  cdata->invar    = (short int *)mCalloc(crunch_len,sizeof(short int));
  cdata->sitepatt = (int *)mCalloc(init_len,sizeof(int ));
  cdata->format   = 0;
  
  cdata->crunch_len = crunch_len;
  cdata->init_len   = init_len;
  cdata->obs_pinvar = .0;
  
  For(j,n_otu)
  {
    cdata->c_seq[j]            = (align *)mCalloc(1,sizeof(align));
    cdata->c_seq[j]->name      = (char *)mCalloc((int)(strlen(sp_names[j])+1),sizeof(char));
    strcpy(cdata->c_seq[j]->name,sp_names[j]);
    cdata->c_seq[j]->state     = (char *)mCalloc(crunch_len*state_len,sizeof(char));
    cdata->c_seq[j]->is_ambigu = (short int *)mCalloc(crunch_len,sizeof(short int));
    
    if(io->datatype==CODON)//!< Added by Marcelo.//!< 
    {
      cdata->c_seq[j]->alternativeCodons=(char **)mCalloc(crunch_len,sizeof(char *));
      //For(i,crunch_len) cdata->c_seq[j]->alternativeCodons[i]=(char *)mCalloc(io->mod->ns+1,sizeof(char));  
    }
  }
  return cdata;
}

/*********************************************************/
/*********************************************************/

calign *Copy_Cseq(calign *ori, option *io)
{
  calign *new;
  int i,j,k,n_otu,c_len;
  char **sp_names;
  
  n_otu = ori->n_otu;
  c_len = ori->crunch_len;
  
  sp_names = (char **)mCalloc(n_otu,sizeof(char *));
  For(i,n_otu)
  {
    sp_names[i] = (char *)mCalloc((int)strlen(ori->c_seq[i]->name)+1,sizeof(char));
    strcpy(sp_names[i],ori->c_seq[i]->name);
  }
  
  new = Make_Cseq(n_otu,c_len+1,io->mod->state_len,ori->init_len,sp_names, io);
  
  new->obs_pinvar = ori->obs_pinvar;
  
  For(i,ori->init_len) new->sitepatt[i] = ori->sitepatt[i];
  
  For(j,ori->crunch_len)
  {
    For(i,ori->n_otu) 
    {
      For(k,io->mod->state_len) new->c_seq[i]->state[j*io->mod->state_len+k] = ori->c_seq[i]->state[j*io->mod->state_len+k];
      new->c_seq[i]->is_ambigu[j] = ori->c_seq[i]->is_ambigu[j];
    }
    
    new->wght[j]   = ori->wght[j];
    new->ambigu[j] = ori->ambigu[j];
    new->invar[j]  = ori->invar[j];
  }
  
  For(i,ori->n_otu)
  {
    new->c_seq[i]->len = ori->c_seq[i]->len;
    strcpy(new->c_seq[i]->name,ori->c_seq[i]->name);
  }
  
  if(io->datatype==CODON)
  {
    For(i,ori->c_seq[0]->len)
    {
      j=-1;
      while(ori->c_seq[i]->alternativeCodons[i][++j]<64);
      For(k,j+1) new->c_seq[i]->alternativeCodons[i][k]=ori->c_seq[i]->alternativeCodons[i][k];
    }
  }
  
  For(i,ori->n_otu) new->c_seq[i]->state[c_len*io->mod->state_len] = '\0'; //!< Why this?
  
  For(i,io->mod->ns) new->b_frq[i] = ori->b_frq[i];
  
  new->init_len           = ori->init_len;
  new->clean_len          = ori->clean_len;
  new->crunch_len         = ori->crunch_len;
  new->n_otu              = ori->n_otu;
  
  For(i,n_otu) free(sp_names[i]);
  free(sp_names);
  
  return new;
}

/*********************************************************/

//! Allocate memory for the optimization options.

/*!
 *    
 */
optimiz *Make_Optimiz()
{
  optimiz *s_opt;                                //!< pointer.
  s_opt = (optimiz *)mCalloc(1,sizeof(optimiz)); //!< Allocate memory.
  return s_opt;                                  //!< Return pointer to the optimization options.
}

/*********************************************************/


int Filexists(char *filename)
{
  FILE *fp;
  fp =fopen(filename,"r");
  if (fp) {
    fclose(fp);
    return 1;
  } else
    return 0;
}

/*********************************************************/

FILE *Openfile(char *filename, int mode)
{
  /* mode = 0 -> read */
  /* mode = 1 -> write */
  /* mode = 2 -> append */
  
  FILE *fp;
  char *s;
  int open_test=0;
  
  /*   s = (char *)mCalloc(T_MAX_FILE,sizeof(char)); */
  
  /*   strcpy(s,filename); */
  
  s = filename;
  
  fp = NULL;
  
  switch(mode)
  {
    case 0 :
    {
      while(!(fp = (FILE *)fopen(s,"r")) && ++open_test<10)
      {
        PhyML_Printf("\n. Can't open file '%s', enter a new name : ",s);
        Getstring_Stdin(s);
      }
      break;
    }
    case 1 :
    {
      fp = (FILE *)fopen(s,"w");
      break;
    }
    case 2 :
    {
      fp = (FILE *)fopen(s,"a");
      break;
    }
      
    default : break;
      
  }
  
  /*   free(s); */
  
  return fp;
}

/*********************************************************/

void Print_Fp_Out(FILE *fp_out, time_t t_beg, time_t t_end, t_tree *tree, option *io, int n_data_set, int num_tree) {
    char *s,*r,*t;
    s=(char *)mCalloc(T_MAX_NAME,sizeof(char));//!< Added by Marcelo.
    r=(char *)mCalloc(T_MAX_NAME,sizeof(char));//!< Added by Marcelo.
    t=(char *)mCalloc(T_MAX_NAME,sizeof(char));//!< Added by Marcelo.
    div_t hour,min;
    int i,j,c;
    char *nucs[] = {"A", "C", "G", "T"};
    token *rootTkn = Emit_Root_Token();
    token *currTkn = rootTkn;
    
    /*   int i; */
    
    /*   For(i,2*tree->n_otu-3) fprintf(fp_out,"\n. * Edge %3d: %f",i,tree->t_edges[i]->l); */
    
    if((!n_data_set) || (!num_tree)) {
        currTkn = Print_Banner_Small(currTkn);
    }
    
    currTkn = Emit_Out_Token(currTkn, "Sequence filename", "inputfile", SCALARTKN, TFSTRING, "%s", Basename(io->in_align_file));
    currTkn = Emit_Out_Token(currTkn, "Data set", "dataset", SCALARTKN, TFNUMERIC, "%g", (double) n_data_set);
    
    if(io->mod->s_opt->random_input_tree) {
        currTkn = Emit_Out_Token(currTkn, "Random init tree #", "RandInitTree", SCALARTKN, TFNUMERIC, "%g", (double)num_tree+1);
    } else if(io->n_trees > 1) {
        currTkn = Emit_Out_Token(currTkn, "Starting tree number", "StartTree", SCALARTKN, TFNUMERIC, "%g", (double)num_tree+1);
    }
    
    if(io->mod->s_opt->opt_topo) {
        if(io->mod->s_opt->topo_search == NNI_MOVE) {
            currTkn = Emit_Out_Token(currTkn, "Tree topology search", "Topologysearch", SCALARTKN, TFSTRING, "%s", "NNIs");
        } else if(io->mod->s_opt->topo_search == SPR_MOVE) {
            currTkn = Emit_Out_Token(currTkn, "Tree topology search", "Topologysearch", SCALARTKN, TFSTRING, "%s", "SSRs");
        } else if(io->mod->s_opt->topo_search == BEST_OF_NNI_AND_SPR) {
            currTkn = Emit_Out_Token(currTkn, "Tree topology search", "Topologysearch", SCALARTKN, TFSTRING, "NNIs and SSRs");
        }
    } else {
        currTkn = Emit_Out_Token(currTkn, "Tree topology search", "Topologysearch", SCALARTKN, TFSTRING, "%s", "fixed");
    }
    
    
    /* was after Sequence file ; moved here FLT */
    
    if(io->in_tree == 2) {
        strcat(strcat(strcat(s,"user tree ("),io->in_tree_file),")");
    } else {
        if(!io->mod->s_opt->random_input_tree) {
            if(io->in_tree == 0)
                strcat(s,"BioNJ");
            if(io->in_tree == 1)
                strcat(s,"parsimony");
        } else {
            strcat(s,"random tree");
        }
    }
    
    if(io->datatype==CODON) {//!< Added by Marcelo
        switch(io->init_DistanceTreeCD) {
            case KOSI07: strcat(s," + GYECMK07\0");break;
            case SCHN05: strcat(s," + GYECMS05\0");break;
            case NUCLEO: strcat(s," + JC69\0");break;
            case ECMUSR: strcat(s," + ECMUSR\0");break;
            default: break;
        }
    }
    currTkn = Emit_Out_Token(currTkn, "Initial tree", "InitTree", SCALARTKN, TFSTRING, "%s", s);
    
    if(tree->io->datatype == NT) {
        currTkn = Emit_Out_Token(currTkn, "Model name", "name", SCALARTKN, TFSTRING, "%s %s", io->mod->modelname, t);
        if(io->mod->whichmodel == CUSTOM) {
            currTkn = Emit_Out_Token(currTkn, "Custom model name", "custname", SCALARTKN, TFSTRING, "%s", io->mod->custom_mod_string, t);
        }
    } else if(tree->io->datatype == AA) {
        currTkn = Emit_Out_Token(currTkn, "Model name", "name", SCALARTKN, TFSTRING, "%s %s", io->mod->modelname, t);
    } else if(tree->io->datatype == CODON) {
        switch(io->mod->freq_model) {
            case F1XSENSECODONS:
                strcpy(t,"F1x");
                sprintf(r,"%d",io->mod->ns);
                strcat(t,r);
                break;
            case F1X4: strcpy(t,"F1x4");
                break;
            case F3X4: strcpy(t,"F3x4");
                break;
            case CF3X4: strcpy(t,"CF3x4");
                break;
            default:
                break;
        }
        
        currTkn = Emit_Out_Token(currTkn, "Model name", "name", SCALARTKN, TFSTRING, "%s %s", io->mod->modelname, t);
        if(io->modeltypeOpt == HLP17){
          currTkn = Emit_Out_Token(currTkn, "Hotspots", "motifs", SCALARTKN, TFSTRING, "%s", io->mod->motifstring); 
          currTkn = Emit_Out_Token(currTkn, "h optimization", "motifs", SCALARTKN, TFSTRING, "%s", io->mod->hotnessstring);
          currTkn = Emit_Out_Token(currTkn, "Partition file", "motifs", SCALARTKN, TFSTRING, "%s", io->mod->partfile);
        }


        
        switch(io->mod->genetic_code)
        {
            case   STANDARD: strcpy(s,"Standard\0");break;
            case       TVMC: strcpy(s,"Vertebrate Mitochondrial\0"); break;
            case       TYMC: strcpy(s,"Yeast Mitochondrial\0"); break;
            case THMPCMCMSC: strcpy(s,"Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma\0"); break;
            case      THIMC: strcpy(s,"Invertebrate Mitochondrial\0"); break;
            case    THCDHNC: strcpy(s,"Ciliate, Dasycladacean and Hexamita Nuclear\0"); break;
            case     THEFMC: strcpy(s,"Echinoderm and Flatworm Mitochondrial\0"); break;
            case      THENC: strcpy(s,"Euplotid Nuclear\0"); break;
            case    THBAPPC: strcpy(s,"Bacterial, Archaeal and Plant Plastid\0"); break;
            case     THAYNC: strcpy(s,"Alternative Yeast Nuclear\0"); break;
            case      THAMC: strcpy(s,"Ascidian Mitochondrial\0"); break;
            case     THAFMC: strcpy(s,"Alternative Flatworm Mitochondrial\0"); break;
            case       BLNC: strcpy(s,"Blepharisma Nuclear\0"); break;
            case       CHMC: strcpy(s,"Chlorophycean Mitochondrial\0"); break;
            case       TRMC: strcpy(s,"Trematode Mitochondrial\0"); break;
            case      SCOMC: strcpy(s,"Scenedesmus obliquus mitochondrial\0"); break;
            case       THMC: strcpy(s,"Thraustochytrium Mitochondrial\0"); break;
            default : Warn_And_Exit("Genetic code not implemented."); break;
        }
        currTkn = Emit_Out_Token(currTkn, "Genetic code", "GenCode", SCALARTKN,TFSTRING,  "%s", s );
    } else {
        currTkn = Emit_Out_Token(currTkn, "Substitution model", "SubstModel", SCALARTKN, TFSTRING, "%s", io->mod->modelname );
    }
    

    currTkn = Emit_Out_Token(currTkn, "Number of taxa", "Taxa", SCALARTKN, TFNUMERIC, "%g", (double)tree->n_otu);
    
    if(io->ratio_test == 1) {
        currTkn = Emit_Out_Token(currTkn, "Branch support", "BranchSupport", SCALARTKN, TFSTRING, "%s", "aLRT statistics");
    } else if(io->ratio_test == 2) {
        currTkn = Emit_Out_Token(currTkn, "Branch support", "BranchSupport", SCALARTKN, TFSTRING, "%s", "aLRT (Chi2-based parametric)");
    } else if(io->ratio_test == 3) {
        currTkn = Emit_Out_Token(currTkn, "Branch support", "BranchSupport", SCALARTKN, TFSTRING, "%s", "aLRT (Minimum of SH-like and Chi2-based)");
    } else if(io->ratio_test == 4) {
        currTkn = Emit_Out_Token(currTkn, "Branch support", "BranchSupport", SCALARTKN, TFSTRING, "%s", "SH-like");
    } else if(io->ratio_test == 5) {
        currTkn = Emit_Out_Token(currTkn, "Branch support", "BranchSupport", SCALARTKN, TFSTRING, "%s", "aBayes");
    }
    
    currTkn = Emit_Out_Token(currTkn, "Log-likelihood", "logL", SCALARTKN, TFNUMERIC, "%.6f", tree->c_lnL);/*was last ; moved here FLT*/
    
    if(tree->io->datatype==AA) {
        if(tree->io->convert_NT_to_AA) {
            currTkn = Emit_Out_Token(currTkn, "Codon Log-likelihood (no gaps)", "CodonLogLNoGaps", SCALARTKN, TFNUMERIC, "%.6f", tree->c_lnL+tree->logLk_correction);
            currTkn = Emit_Out_Token(currTkn, "Codon Log-likelihood (with gaps)", "CodonLogLGaps", SCALARTKN, TFNUMERIC, "%.6f", tree->c_lnL+tree->logLk_correction_gaps+tree->logLk_correction);
            currTkn = Emit_Out_Token(currTkn, "Codon Log-lk. (approx, no gaps)", "CodonLogLApproxNoGaps", SCALARTKN, TFNUMERIC, "%.6f", tree->c_lnL+tree->logLk_correction_approx);
            currTkn = Emit_Out_Token(currTkn, "Codon Log-lk. (approx, with gaps)", "CodonLogLApproxGaps", SCALARTKN, TFNUMERIC, "%.6f", tree->c_lnL+tree->logLk_correction_gaps_approx+tree->logLk_correction_approx);
        } else {
            currTkn = Emit_Out_Token(currTkn, "Codon Log-lk. (approx, no gaps)", "CodonLogLApproxNoGaps", SCALARTKN, TFNUMERIC, "%.6f", tree->c_lnL+tree->logLk_correction_approx);
            currTkn = Emit_Out_Token(currTkn, "Codon Log-lk. (approx, with gaps)", "CodonLogLApproxGaps", SCALARTKN, TFNUMERIC, "%.6f", tree->c_lnL+tree->logLk_correction_gaps_approx+tree->logLk_correction_approx);
        }
    } else if(tree->io->datatype==NT) {
        currTkn = Emit_Out_Token(currTkn, "Codon Log-lk. (equivalent)", "CodonLogLEquiv", SCALARTKN, TFNUMERIC, "%.6f", tree->c_lnL);
    }
    Unconstraint_Lk(tree);
    
    currTkn = Emit_Out_Token(currTkn, "Unconstrained Log-likelihood", "Unconstrained lnL", SCALARTKN, TFNUMERIC, "%.6f", tree->unconstraint_lk);
    
    currTkn = Emit_Out_Token(currTkn, "Parsimony", "Parsimony", SCALARTKN, TFNUMERIC, "%g", (double)tree->c_pars);
    
    currTkn = Emit_Out_Token(currTkn, "Tree size", "Tree size", SCALARTKN, TFNUMERIC, "%.5f", tree->size);
    
    if(tree->mod->io->datatype==CODON && tree->mod->whichrealmodel==PCM) {
        currTkn = Emit_Out_Token(currTkn, "Principal components", "PCs", TBLSTARTTKN, TFNUMERIC, "%g", (double)tree->mod->npcs);
        for(i = 0; i < tree->mod->npcs; i++ ) {
            char * buf;
            int temp2 = asprintf(&buf, "PC%i", i); //changed by Ken 12/1/2017
            currTkn = Emit_Out_Token(currTkn, buf, buf, SCALARTKN, TFNUMERIC, "%10f", tree->mod->pcsC[i]);
        }
        currTkn = Emit_Out_Token(currTkn, "Principal components", "PCs", TBLENDTKN, TFSTRING, "%s", "");
    }
    
    if(tree->mod->io->datatype==CODON) {
        if((tree->mod->n_w_catg == 1)&&(tree->mod->n_catg > 1)) {
            currTkn = Emit_Out_Token(currTkn, "Discrete gamma model", "DiscGammaMod", SETSTARTTKN, TFSTRING, "%s", "Yes");
            currTkn = Emit_Out_Token(currTkn, "Number of categories", "NrCat", SCALARTKN, TFNUMERIC, "%g", (double)tree->mod->n_catg );
            currTkn = Emit_Out_Token(currTkn, "Gamma shape parameter", "Alpha", SCALARTKN, TFNUMERIC, "%.3f", tree->mod->alpha );
            currTkn = Emit_Out_Token(currTkn, "Discrete gamma model", "DiscGammaMod", SETENDTKN, TFSTRING, "%s", "Yes");
        } else if((tree->mod->n_w_catg > 1) && (tree->mod->omegaSiteVar==DGAMMAK)) {
            currTkn = Emit_Out_Token(currTkn, "Discrete gamma model", "DiscGammaMod", SETSTARTTKN, TFSTRING, "%s", "Yes");
            currTkn = Emit_Out_Token(currTkn, "Number of categories", "NrCat", SCALARTKN, TFNUMERIC, "%g", (double)tree->mod->n_w_catg );
            currTkn = Emit_Out_Token(currTkn, "Gamma shape (alpha)", "Alpha", SCALARTKN, TFNUMERIC, "%.3f", tree->mod->alpha );
            currTkn = Emit_Out_Token(currTkn, "Gamma shape (beta)", "Beta", SCALARTKN, TFNUMERIC, "%.3f", tree->mod->beta );
            currTkn = Emit_Out_Token(currTkn, "Discrete gamma model", "DiscGammaMod", SETENDTKN, TFSTRING, "%s", "Yes");
        }
    } else {
        if(tree->mod->n_catg>1) {
            currTkn = Emit_Out_Token(currTkn, "Discrete gamma model", "DiscGammaMod", SETSTARTTKN, TFSTRING, "%s", "Yes");
            currTkn = Emit_Out_Token(currTkn, "Number of categories", "NrCat", SCALARTKN, TFNUMERIC, "%g", (double)tree->mod->n_catg );
            currTkn = Emit_Out_Token(currTkn, "Gamma shape parameter", "Alpha", SCALARTKN, TFNUMERIC, "%.3f", (double)tree->mod->alpha );
            currTkn = Emit_Out_Token(currTkn, "Discrete gamma model", "DiscGammaMod", SETENDTKN, TFSTRING, "%s", "Yes");
        }
    }
    
    currTkn = Emit_Out_Token(currTkn, "Proportion of invariant", "Pinvar", SCALARTKN, TFNUMERIC, "%.3f", tree->mod->pinvar );
    
    /* was before Discrete gamma model ; moved here FLT */
    if((tree->mod->whichmodel == K80) || (tree->mod->whichmodel == HKY85) || (tree->mod->whichmodel == F84)) {
        currTkn = Emit_Out_Token(currTkn, "Transition/transversion ratio", "Kappa", SCALARTKN, TFNUMERIC, "%.6f", tree->mod->kappa);
    } else if(tree->mod->whichmodel == TN93) {
        currTkn = Emit_Out_Token(currTkn, "Transition/transversion ratio for purines", "KappaPurines", SCALARTKN, TFNUMERIC, "%.5f",
                                 tree->mod->kappa*2.*tree->mod->lambda/(1.+tree->mod->lambda) );
        currTkn = Emit_Out_Token(currTkn, "Transition/transversion ratio for pyrimidines", "KappaPyrimidines", SCALARTKN, TFNUMERIC, "%.5f",
                                 tree->mod->kappa*2./(1.+tree->mod->lambda) );
    } else if(tree->io->datatype == CODON) {
        if((tree->mod->whichmodel != GYECMK07) && tree->mod->whichmodel != GYECMK07F  && tree->mod->whichmodel != GYECMK07WK  && tree->mod->whichmodel != GYECMK07WKF  && tree->mod->whichmodel != GYECMS05  && tree->mod->whichmodel != GYECMS05F  && tree->mod->whichmodel != GYECMS05WK  && tree->mod->whichmodel != GYECMS05WKF &&
           tree->mod->whichmodel != MGECMK07  && tree->mod->whichmodel != MGECMK07F  && tree->mod->whichmodel != MGECMK07WK  && tree->mod->whichmodel != MGECMK07WKF  && tree->mod->whichmodel != MGECMS05  && tree->mod->whichmodel != MGECMS05F  && tree->mod->whichmodel != MGECMS05WK  && tree->mod->whichmodel != MGECMS05WKF &&
           tree->mod->whichmodel != YAPECMK07 && tree->mod->whichmodel != YAPECMK07F && tree->mod->whichmodel != YAPECMK07WK && tree->mod->whichmodel != YAPECMK07WKF && tree->mod->whichmodel != YAPECMS05 && tree->mod->whichmodel != YAPECMS05F && tree->mod->whichmodel != YAPECMS05WK && tree->mod->whichmodel != YAPECMS05WKF  &&
           tree->mod->whichmodel != GYECMUSR && tree->mod->whichmodel != GYECMUSRF  && tree->mod->whichmodel != GYECMUSRWK  && tree->mod->whichmodel != MGECMUSRWK  &&tree->mod->whichmodel != GYECMUSRWKF  && tree->mod->whichmodel != MGECMUSRF  && tree->mod->whichmodel != MGECMUSRWKF && tree->mod->whichmodel != YAPECMUSR && tree->mod->whichmodel != YAPECMUSRF && tree->mod->whichmodel != YAPECMUSRWK && tree->mod->whichmodel != YAPECMUSRWKF && !tree->mod->pcaModel)
        {
            currTkn = Emit_Out_Token(currTkn, "Transition/transversion ratio", "Kappa", SCALARTKN, TFNUMERIC, "%.9f", tree->mod->kappa);
            
            if(tree->mod->n_w_catg == 1) {
            		int omegai; //Added by Ken 23/8
            		for(omegai=0;omegai<tree->mod->nomega_part;omegai++){
            			char *buf = mCalloc(T_MAX_OPTION,sizeof(char));
            			int temp = asprintf(&buf, "Omega %d %s", omegai, tree->mod->partNames[omegai]);//changed by Ken 12/1/2017
            			currTkn = Emit_Out_Token(currTkn, buf, buf, SCALARTKN, TFNUMERIC, "%.9f", tree->mod->omega_part[omegai]);
            		}

            } else {
                currTkn = Emit_Out_Token(currTkn, "Nonsynonmous/synonymous ratio", "Omega", SETSTARTTKN, TFSTRING, "%s", "");
                for(i = 0; i < tree->mod->n_w_catg; i++ ) {
                    currTkn = Emit_Out_Token(currTkn, "Omega / Probability", "OmegaProb", TBLSTARTTKN, TFNUMERIC, "%i", 2);
                    
                    char *buf;
                    int temp = asprintf(&buf, "o%i", i);//changed by Ken 12/1/2017
                    currTkn = Emit_Out_Token(currTkn, buf, buf, SCALARTKN, TFNUMERIC, "%.6f", tree->mod->omegas[i]);
                    
                    char *buf2;
                    int temp2 = asprintf(&buf2, "p%i", i);//changed by Ken 12/1/2017
                    currTkn = Emit_Out_Token(currTkn, buf2, buf2, SCALARTKN, TFNUMERIC, "%.6f", tree->mod->prob_omegas[i]);
                    
                    currTkn = Emit_Out_Token(currTkn, "Omega / Probability", "OmegaProb", TBLENDTKN, TFNUMERIC, "%i", 2);
                }
                currTkn = Emit_Out_Token(currTkn, "Nonsynonmous/synonymous ratio", "Omega", SETENDTKN, TFSTRING, "%s", "");
            }

            //Prints out hotspot model information and h values
            //Modified by Ken 22/7/2016
            if(io->modeltypeOpt==HLP17){
              int partsite=0;
              printf("\n");
              for(partsite=0;partsite<io->tree->n_pattern;partsite++){
                printf("%d ",tree->mod->partIndex[partsite]);
              }
              printf("\n");


            	  currTkn = Emit_Out_Token(currTkn, "Hotspot model\t\th_index\toptimized?\th_value", "motifs", SETSTARTTKN, TFSTRING, "%s", "" );
            	  int mot;
            	  for(mot=0;mot<io->mod->nmotifs;mot++){
                       currTkn = Emit_Out_Token(currTkn, "", "", TBLSTARTTKN, TFNUMERIC, "%i", 0);
                       char *info = malloc(100);
                       char motifh[10];
                       sprintf(motifh,"%d",io->mod->motif_hotness[mot]);
                       char motife[10];
                       sprintf(motife,"%d",io->mod->hoptindex[io->mod->motif_hotness[mot]]);
                       strcpy(info, "Motif:\t");
                       strcat(info, io->mod->motifs[mot]);
                       strcat(info, "\t\t");
                       strcat(info, motifh);
                       strcat(info, "\t\t");
                       strcat(info, motife);
                       strcat(info, "\t\t");
                       currTkn = Emit_Out_Token(currTkn, info, "h", SCALARTKN, TFNUMERICNEQ, "%.5f", io->mod->hotness[io->mod->motif_hotness[mot]]);
                      currTkn = Emit_Out_Token(currTkn, "", "", TBLENDTKN, TFNUMERIC, "%i", 4 );
            	  }
            	currTkn = Emit_Out_Token(currTkn, "Motif model", "motifs", SETENDTKN, TFSTRING, "%s", "" );

            }
        }

        else if(!tree->mod->pcaModel)
        {
            switch( io->kappaECM ) {
                case kap1: currTkn = Emit_Out_Token(currTkn, "Emp. Transition/transverion ratio", "EmpKappa", SCALARTKN, TFSTRING, "%s", "Embedded in the Q matrix." );
                    break;
                case kap2: {
                    currTkn = Emit_Out_Token(currTkn, "Emp. Transition/transverion ratio", "EmpKappa", TBLSTARTTKN, TFNUMERIC, "%i", 1 );
                    currTkn = Emit_Out_Token(currTkn, "nts", "nts", SCALARTKN, TFNUMERIC, "%.4f", io->mod->pkappa[0]);
                    currTkn = Emit_Out_Token(currTkn, "Emp. Transition/transverion ratio", "EmpKappa", TBLENDTKN, TFNUMERIC, "%i", 1 );
                    break;
                }
                case kap3: {
                    currTkn = Emit_Out_Token(currTkn, "Emp. Transition/transverion ratio", "EmpKappa", TBLSTARTTKN, TFNUMERIC, "%i", 1 );
                    currTkn = Emit_Out_Token(currTkn, "ntv", "ntv", SCALARTKN, TFNUMERIC, "%.4f", io->mod->pkappa[0]);
                    currTkn = Emit_Out_Token(currTkn, "Emp. Transition/transverion ratio", "EmpKappa", TBLENDTKN, TFNUMERIC, "%i", 1 );
                    break;
                }
                case kap4: {
                    currTkn = Emit_Out_Token(currTkn, "Emp. Transition/transverion ratio", "EmpKappa", TBLSTARTTKN, TFNUMERIC, "%i", 2 );
                    currTkn = Emit_Out_Token(currTkn, "nts", "nts", SCALARTKN, TFNUMERIC, "%.4f", io->mod->pkappa[0]);
                    currTkn = Emit_Out_Token(currTkn, "ntv", "ntv", SCALARTKN, TFNUMERIC, "%.4f", io->mod->pkappa[1]);
                    currTkn = Emit_Out_Token(currTkn, "Emp. Transition/transverion ratio", "EmpKappa", TBLENDTKN, TFNUMERIC, "%i", 2 );
                    break;
                }
                case kap5: {
                    currTkn = Emit_Out_Token(currTkn, "Emp. Transition/transverion ratio", "EmpKappa", TBLSTARTTKN, TFSTRING, "%s", 3 );
                    for(c = 0; c < 9; c++) {
                        char *buf;
                        int temp = asprintf(&buf, "Case%i", c);//changed by Ken 12/1/2017
                        currTkn = Emit_Out_Token(currTkn, buf, buf, SCALARTKN, TFNUMERIC, "%.4f", io->mod->pkappa[c]);
                    }
                    currTkn = Emit_Out_Token(currTkn, "Emp. Transition/transverion ratio", "EmpKappa", TBLENDTKN, TFNUMERIC, "%i", 3 );
                    break;
                }
                case kap6: currTkn = Emit_Out_Token(currTkn, "Emp. Transition/transverion ratio", "EmpKappa", SCALARTKN, TFNUMERIC, "%.4f", tree->mod->pkappa[0] );
                    break;
                default:break;
            }
            
            //      if((tree->mod->s_opt->opt_omega==NO && tree->mod->s_opt->opt_prob_omega==NO ) || tree->mod->whichmodel==YAPECMUSR|| tree->mod->whichmodel==YAPECMUSRF|| tree->mod->whichmodel==MGECMUSRF ||tree->mod->whichmodel==GYECMUSRF || tree->mod->whichmodel==GYECMUSR|| tree->mod->whichmodel==MGECMUSR || tree->mod->whichmodel==GYECMK07 || tree->mod->whichmodel==GYECMK07F || tree->mod->whichmodel==GYECMS05 || tree->mod->whichmodel==GYECMK07F || tree->mod->whichmodel==MGECMK07 || tree->mod->whichmodel==MGECMK07F || tree->mod->whichmodel==MGECMS05 || tree->mod->whichmodel==MGECMK07F || tree->mod->whichmodel==YAPECMK07 || tree->mod->whichmodel==YAPECMK07F || tree->mod->whichmodel==YAPECMS05 || tree->mod->whichmodel==YAPECMK07F)
            if( (tree->mod->s_opt->opt_omega == NO) &&
               (tree->mod->s_opt->opt_prob_omega == NO) &&
               (tree->mod->omega_part[0] == 1.0 ) ) { //modified by Ken 18/8
                currTkn = Emit_Out_Token(currTkn, "Emp. nonsyn./syn.", "EmpOmega", SCALARTKN, TFSTRING, "%s", "Embedded in the Q matrix");
            }
            else
            {
                if(tree->mod->n_w_catg == 1)
                {
                	if(tree->mod->nparts > 1){printf("options not compatible with partitioned model error 3\n");exit(EXIT_FAILURE);}
                    currTkn = Emit_Out_Token(currTkn, "Emp. corrected nonsyn./syn. ratio", "EmpCorrOmega", SCALARTKN, TFNUMERIC, "%.6f", Omega_ECMtoMmechModels(tree->mod->pi, tree->mod->qmat_part[0], tree->mod->qmat_buff_part[0], tree->mod->ns, tree->mod->n_w_catg));
                    if( (tree->mod->omegaSiteVar != NOOMEGA) & (tree->mod->s_opt->opt_omega == NO) ) {
                        currTkn = Emit_Out_Token(currTkn, "User defined nonsyn./syn. ratio", "UserdefOmega", SCALARTKN, TFNUMERIC, "%.6f", tree->mod->omega_part[0] ); //modified by Ken 17/8
                    }
                }
                else
                {
                	if(tree->mod->nparts > 1){printf("options not compatible with partitioned model error 4\n");exit(EXIT_FAILURE);}
                    currTkn = Emit_Out_Token(currTkn, "Estimated dn/ds rate ratio (w) values", "EstimatedOmega", TBLSTARTTKN, TFNUMERIC, "%i", 2);
                    
                    for(i = 0; i < tree->mod->n_w_catg; i++ ) {
                        char *buf;
                        int temp2 = asprintf(&buf, "o%i", i);//changed by Ken 12/1/2017
                        currTkn = Emit_Out_Token(currTkn, buf, buf, SCALARTKN, TFNUMERIC, "%.10f",
                                                 Omega_ECMtoMmechModels( tree->mod->pi, tree->mod->qmat +
                                                                        i*tree->mod->ns * tree->mod->ns,
                                                                        tree->mod->qmat_buff_part[0] +
                                                                        i*tree->mod->ns*tree->mod->ns,
                                                                        tree->mod->ns, tree->mod->n_w_catg )
                                                 );
                        
                        char *buf2;
                        int temp = asprintf(&buf2, "p%i", i);//changed by Ken 12/1/2017
                        currTkn = Emit_Out_Token(currTkn, buf2, buf2, SCALARTKN, TFNUMERIC, "%.10f", tree->mod->prob_omegas[i]);
                    }
                    currTkn = Emit_Out_Token(currTkn, "Estimated dn/ds rate ratio (w) values", "EstimatedOmega", TBLENDTKN, TFNUMERIC, "%i", 2);
                }
            }
        }
    }
    if(tree->io->datatype == NT)
    {
        currTkn = Emit_Out_Token(currTkn, "Nucleotides frequencies", "NucFreqs", TBLSTARTTKN, TFNUMERIC, "%i", 4);
        currTkn = Emit_Out_Token(currTkn, "f(A)", "f(A)", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->pi[0] );
        currTkn = Emit_Out_Token(currTkn, "f(C)", "f(C)", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->pi[1] );
        currTkn = Emit_Out_Token(currTkn, "f(G)", "f(G)", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->pi[2] );
        currTkn = Emit_Out_Token(currTkn, "f(T)", "f(T)", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->pi[3] );
        currTkn = Emit_Out_Token(currTkn, "Nucleotides frequencies", "NucFreqs", TBLENDTKN, TFNUMERIC, "%i", 4);
    }else if(tree->io->datatype == CODON)
    {
        char s1[4];
        phydbl val1;
        
        switch( io->mod->freq_model ) {
            case F1XSENSECODONS:
                break;
            case F1X4:
                currTkn = Emit_Out_Token(currTkn, "Nucleotides frequencies", "NucFreqs", TBLSTARTTKN, TFNUMERIC, "%i", 4);
                currTkn = Emit_Out_Token(currTkn, "f(T)", "f(T)", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->base_freq[0] );
                currTkn = Emit_Out_Token(currTkn, "f(C)", "f(C)", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->base_freq[1] );
                currTkn = Emit_Out_Token(currTkn, "f(A)", "f(A)", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->base_freq[2] );
                currTkn = Emit_Out_Token(currTkn, "f(G)", "f(G)", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->base_freq[3] );
                currTkn = Emit_Out_Token(currTkn, "Nucleotides frequencies", "NucFreqs", TBLENDTKN, TFNUMERIC, "%i", 4);
                break;
            case F3X4:
            case CF3X4:
                currTkn = Emit_Out_Token(currTkn, "Nucleotides frequencies", "NucFreqs", SETSTARTTKN, TFSTRING, "%s", "");
                currTkn = Emit_Out_Token(currTkn, "Position 1", "Pos1", TBLSTARTTKN, TFNUMERIC, "%i", 4);
                currTkn = Emit_Out_Token(currTkn, "f(T1)", "f(T1)", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->base_freq[0] );
                currTkn = Emit_Out_Token(currTkn, "f(C1)", "f(C1)", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->base_freq[1] );
                currTkn = Emit_Out_Token(currTkn, "f(A1)", "f(A1)", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->base_freq[2] );
                currTkn = Emit_Out_Token(currTkn, "f(G1)", "f(G1)", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->base_freq[3] );
                currTkn = Emit_Out_Token(currTkn, "Position 1", "Pos1", TBLENDTKN, TFNUMERIC, "%i", 4);
                currTkn = Emit_Out_Token(currTkn, "Position 2", "Pos2", TBLSTARTTKN, TFNUMERIC, "%i", 4);
                currTkn = Emit_Out_Token(currTkn, "f(T2)", "f(T2)", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->base_freq[4] );
                currTkn = Emit_Out_Token(currTkn, "f(C2)", "f(C2)", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->base_freq[5] );
                currTkn = Emit_Out_Token(currTkn, "f(A2)", "f(A2)", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->base_freq[6] );
                currTkn = Emit_Out_Token(currTkn, "f(G2)", "f(G2)", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->base_freq[7] );
                currTkn = Emit_Out_Token(currTkn, "Position 2", "Pos2", TBLENDTKN, TFNUMERIC, "%i", 4);
                currTkn = Emit_Out_Token(currTkn, "Position 3", "Pos3", TBLSTARTTKN, TFNUMERIC, "%i", 4);
                currTkn = Emit_Out_Token(currTkn, "f(T3)", "f(T3)", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->base_freq[8] );
                currTkn = Emit_Out_Token(currTkn, "f(C3)", "f(C3)", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->base_freq[9] );
                currTkn = Emit_Out_Token(currTkn, "f(A3)", "f(A3)", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->base_freq[10] );
                currTkn = Emit_Out_Token(currTkn, "f(G3)", "f(G3)", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->base_freq[11] );
                currTkn = Emit_Out_Token(currTkn, "Position 3", "Pos3", TBLENDTKN, TFNUMERIC, "%i", 4);
                currTkn = Emit_Out_Token(currTkn, "Nucleotides frequencies", "NucFreqs", SETENDTKN, TFSTRING, "%s", "");
                break;
            default:
                break;
        }
        currTkn = Emit_Out_Token(currTkn, "Codon frequencies", "CodonFreqs", SETSTARTTKN, TFSTRING, "%s", "" );
        for(i = 0; i < 64; i++) {
            if(i == 0 || i%4 == 0) currTkn = Emit_Out_Token(currTkn, "", "", TBLSTARTTKN, TFNUMERIC, "%i", 4 );
            char *buf;
            Sprint_codon(s1, i);
            int temp = asprintf(&buf, "f(%s)", s1);//changed by Ken 12/1/2017
            val1 = (stopCodons[i])?0.0:tree->mod->pi[indexSenseCodons[i]];
            currTkn = Emit_Out_Token(currTkn, buf, buf, SCALARTKN, TFNUMERIC, "%8.8f", val1);
            if((i-3)%4 == 0) currTkn = Emit_Out_Token(currTkn, "", "", TBLENDTKN, TFNUMERIC, "%i", 4 );
        }
        if(currTkn->type != TBLENDTKN) {
            currTkn = Emit_Out_Token(currTkn, "", "", TBLENDTKN, TFNUMERIC, "%i", 4 );
        }
        currTkn = Emit_Out_Token(currTkn, "Codon frequencies", "CodonFreqs", SETENDTKN, TFSTRING, "%s", "");
    }

    /*****************************************/
    if((tree->mod->whichmodel == GTR) ||
       (tree->mod->whichmodel == CUSTOM))
    {
        
        Update_Qmat_GTR(tree->mod->rr,
                        tree->mod->rr_val,
                        tree->mod->rr_num,
                        tree->mod->pi,
                        tree->mod->qmat);
        
        currTkn = Emit_Out_Token(currTkn, "GTR relative rate parameters", "GTRRelRates", TBLSTARTTKN, TFNUMERIC, "%i", 3 );
        currTkn = Emit_Out_Token(currTkn, "A <-> C", "AC", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->rr[0]);
        currTkn = Emit_Out_Token(currTkn, "A <-> G", "AG", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->rr[1]);
        currTkn = Emit_Out_Token(currTkn, "A <-> T", "AT", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->rr[2]);
        currTkn = Emit_Out_Token(currTkn, "C <-> G", "CG", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->rr[3]);
        currTkn = Emit_Out_Token(currTkn, "C <-> T", "CT", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->rr[4]);
        currTkn = Emit_Out_Token(currTkn, "G <-> T", "GT", SCALARTKN, TFNUMERIC, "%8.5f", tree->mod->rr[5]);
        currTkn = Emit_Out_Token(currTkn, "GTR relative rate parameters", "GTRRelRates", TBLENDTKN, TFNUMERIC, "%i", 3 );
        
        currTkn = Emit_Out_Token(currTkn, "Instantaneous rate matrix", "InstRateMatrix", SETSTARTTKN, TFSTRING, "%s", "" );
        currTkn = Emit_Out_Token(currTkn, "[A-C-G-T]", "A-C-G-T", COMMENTTKN, TFSTRING, "%s", "" );
        
        for(i = 0; i < 4; i++) {
            currTkn = Emit_Out_Token(currTkn, nucs[i], nucs[i], TBLSTARTTKN, TFNUMERIC, "%i", 4);
            for(j = 0; j < 4; j++) {
                currTkn = Emit_Out_Token(currTkn, "", "", SCALARTKN, TFNUMERIC, "%8.5f  ", tree->mod->qmat[i*4+j]);
            }
            currTkn = Emit_Out_Token(currTkn, nucs[i], nucs[i], TBLENDTKN, TFNUMERIC, "%i", 4);
        }
        currTkn = Emit_Out_Token(currTkn, "Instantaneous rate matrix", "InstRateMatrix", SETENDTKN, TFSTRING, "%s", "" );
    }
    /*****************************************/
    
    
    //     if(io->ratio_test == 1)
    //     {
    //       PhyML_Fprintf(fp_out,". aLRT statistics to test branches");
    //     }
    //     else if(io->ratio_test == 2)
    //     {
    //       PhyML_Fprintf(fp_out,". aLRT branch supports (cubic approximation, mixture of Chi2s distribution)");
    //     }
    
    
    
    
    hour = div((int)(t_end-t_beg),3600);
    min  = div((int)(t_end-t_beg),60  );
    
    min.quot -= hour.quot*60;
    
    if(tree->io->datatype==CODON)
    {
#ifdef BLAS
        Emit_Out_Token(currTkn, "Code optimization", "CodeOptim", SCALARTKN, TFSTRING, "%s", "BLAS and LAPACK");
#elif BLAS_OMP
        Emit_Out_Token(currTkn, "Code optimization", "CodeOptim", SCALARTKN, TFSTRING, "%s", "BLAS, LAPACK and OpenMP");
#elif OMP
        Emit_Out_Token(currTkn, "Code optimization", "CodeOptim", SCALARTKN, TFSTRING, "%s", "OpenMP");
#else
        Emit_Out_Token(currTkn, "Code optimization", "CodeOptim", SCALARTKN, TFSTRING, "%s", "none");
#endif
    }
    
    currTkn = Emit_Out_Token(currTkn, "Time used", "TimeUsed", SCALARTKN, TFSTRING, "%dh%dm%ds", hour.quot,min.quot,(int)(t_end-t_beg)%60);
    currTkn = Emit_Out_Token(currTkn, "Seconds", "TimeUsedSec", SCALARTKN, TFNUMERIC, "%g", (double)(t_end-t_beg));
    
    if(io->fp_out_compare){
        char m[100];
        char n[100];
        char z[100];
        
        if(io->datatype==CODON)
        {
#ifdef BLAS
            strcpy(n,"BLAS\0");
#elif BLAS_OMP
            strcpy(n,"BLAS+OMP\0");
#elif OMP
            strcpy(n,"OMP\0");
#else
            strcpy(n,"Null\0");
#endif
            
            switch(tree->mod->s_opt->topo_search)
            {
                case SPR_MOVE: strcpy(m,"SPR\0");break;
                case NNI_MOVE: strcpy(m,"NNI\0");break;
                case BEST_OF_NNI_AND_SPR: strcpy(m,"BEST\0");break;
                default: break;
            }
            
            if(tree->io->heuristicExpm) strcpy(z,"TAYLOR\0");
            else if(tree->io->expm==SSPADE) strcpy(z,"PADE\0");
            else strcpy(z,"EIGEN\0");
            
            fprintf(io->fp_out_compare,"%d;",tree->n_otu); //!< Number of Taxa
            fprintf(io->fp_out_compare,"%d;",tree->data->init_len); //!< Sequence length
            fprintf(io->fp_out_compare,"%s;",io->mod->modelname); //!< Codon Model
            fprintf(io->fp_out_compare,"%s;",t); //!< Frequency model
            switch(io->eq_freq_handling)
            {
                case EMPIRICAL: fprintf(io->fp_out_compare,"empirical;");break;
                case OPTIMIZE:fprintf(io->fp_out_compare,"optimize;");break;
                case MODEL:fprintf(io->fp_out_compare,"model;");break;
                case USER:fprintf(io->fp_out_compare,"user;");break;
            }
            //if(io->mod->s_opt->opt_state_freq) fprintf(io->fp_out_compare,"optimized;"); //!< Opt frequencies?
            //else fprintf(io->fp_out_compare,"empiric;"); //!< Opt frequencies?
            
            fprintf(io->fp_out_compare,"%s;",m); //!< SEARCH method
            fprintf(io->fp_out_compare,"%s;",n); //!< OPT method
            fprintf(io->fp_out_compare,"%s;",z); //!< EXPM method
            fprintf(io->fp_out_compare,"%.6f;",tree->c_lnL); //!< Tree Likelihood
            fprintf(io->fp_out_compare,"%.6f;",tree->c_lnL); //!< comparable Tree Likelihood
            fprintf(io->fp_out_compare,"%.6f;",tree->size); //!< Tree length
            
            if(io->mod->n_w_catg>1) fprintf(io->fp_out_compare,"1;"); //!< Substitution rate categories
            else if (io->mod->n_w_catg==1) fprintf(io->fp_out_compare,"%d;",io->mod->n_catg); //!< Substitution rate categories
            
            if(io->mod->n_catg>1 && io->mod->n_w_catg==1) fprintf(io->fp_out_compare,"%.3f;",io->mod->alpha); //!< Substitution rate categories Gamma param.
            else fprintf(io->fp_out_compare,";"); //!< Substitution rate categories Gamma param.
            
            if(
               tree->mod->whichmodel!=GYECMK07  && tree->mod->whichmodel!=GYECMK07WK  && tree->mod->whichmodel!=GYECMK07F  && tree->mod->whichmodel!=GYECMK07WKF  && tree->mod->whichmodel!=GYECMS05  && tree->mod->whichmodel!=GYECMS05WK  && tree->mod->whichmodel!=GYECMS05F  && tree->mod->whichmodel!=GYECMS05WKF &&
               tree->mod->whichmodel!=MGECMK07  && tree->mod->whichmodel!=MGECMK07WK  && tree->mod->whichmodel!=MGECMK07F  && tree->mod->whichmodel!=MGECMK07WKF  && tree->mod->whichmodel!=MGECMS05  && tree->mod->whichmodel!=MGECMS05WK  && tree->mod->whichmodel!=MGECMS05F  && tree->mod->whichmodel!=MGECMS05WKF &&
               tree->mod->whichmodel!=YAPECMK07 && tree->mod->whichmodel!=YAPECMK07WK && tree->mod->whichmodel!=YAPECMK07F && tree->mod->whichmodel!=YAPECMK07WKF && tree->mod->whichmodel!=YAPECMS05 && tree->mod->whichmodel!=YAPECMS05WK && tree->mod->whichmodel!=YAPECMS05F && tree->mod->whichmodel!=YAPECMS05WKF &&
               tree->mod->whichmodel != GYECMUSR && tree->mod->whichmodel != GYECMUSRF  && tree->mod->whichmodel != GYECMUSRWK  && tree->mod->whichmodel != GYECMUSRWKF  && tree->mod->whichmodel != MGECMUSRF  && tree->mod->whichmodel != MGECMUSRWKF && tree->mod->whichmodel != MGECMUSR && tree->mod->whichmodel != MGECMUSRWK  && tree->mod->whichmodel != YAPECMUSR && tree->mod->whichmodel != YAPECMUSRF && tree->mod->whichmodel != YAPECMUSRWK && tree->mod->whichmodel != YAPECMUSRWKF
               )
            {
                fprintf(io->fp_out_compare,"%.4f;",tree->mod->kappa); //!< Kappa
            }
            else
            {
                if(tree->mod->whichmodel==YAPECMUSRF ||tree->mod->whichmodel==YAPECMUSR ||tree->mod->whichmodel==MGECMUSR ||tree->mod->whichmodel==MGECMUSRF ||tree->mod->whichmodel==GYECMUSRF ||tree->mod->whichmodel==GYECMUSR ||tree->mod->whichmodel==GYECMK07 || tree->mod->whichmodel==GYECMK07F || tree->mod->whichmodel==GYECMS05 || tree->mod->whichmodel==GYECMS05F || tree->mod->whichmodel==MGECMK07 || tree->mod->whichmodel==MGECMK07F || tree->mod->whichmodel==MGECMS05 || tree->mod->whichmodel==MGECMS05F || tree->mod->whichmodel==YAPECMK07 || tree->mod->whichmodel==YAPECMK07F || tree->mod->whichmodel==YAPECMS05 || tree->mod->whichmodel==YAPECMS05F )
                {
                    fprintf(io->fp_out_compare,"1;"); //!< Kappa
                }
                else
                {
                    switch(tree->io->kappaECM) //!< Kappa ECM Kosiol 2007.
                    {
                        case kap1: fprintf(io->fp_out_compare,"1;"); break;
                        case kap2: fprintf(io->fp_out_compare,"%.4f;",tree->mod->pkappa[0]); break;
                        case kap3: fprintf(io->fp_out_compare,"%.4f;",tree->mod->pkappa[0]); break;
                        case kap4: fprintf(io->fp_out_compare,"%.4f %.4f;",tree->mod->pkappa[0],tree->mod->pkappa[1]); break;
                        case kap5: For(i,9) fprintf(io->fp_out_compare,"%.4f ",tree->mod->pkappa[i]);  fprintf(io->fp_out_compare,";");  break;
                        case kap6: fprintf(io->fp_out_compare,"%.4f;",tree->mod->pkappa[0]); break;
                        default:break;
                    }
                }
            }
            
            switch(tree->mod->omegaSiteVar)
            {
                case DM0: fprintf(io->fp_out_compare,"M0;");break;
                case DMODELK: fprintf(io->fp_out_compare,"M3;");break;
                case DGAMMAK: fprintf(io->fp_out_compare,"M5;");break;
                default: break;
            }
            
            if(tree->mod->n_w_catg==1)
            {
                if(
                   tree->mod->whichmodel!=GYECMK07  && tree->mod->whichmodel!=GYECMK07F  && tree->mod->whichmodel!=GYECMS05  && tree->mod->whichmodel!=GYECMS05F &&
                   tree->mod->whichmodel!=MGECMK07  && tree->mod->whichmodel!=MGECMK07F  && tree->mod->whichmodel!=MGECMS05  && tree->mod->whichmodel!=MGECMS05F &&
                   tree->mod->whichmodel!=YAPECMK07 && tree->mod->whichmodel!=YAPECMK07F && tree->mod->whichmodel!=YAPECMS05 && tree->mod->whichmodel!=YAPECMS05F &&
                   tree->mod->whichmodel!=YAPECMUSRF &&tree->mod->whichmodel!=YAPECMUSR && tree->mod->whichmodel!=MGECMUSR &&tree->mod->whichmodel!=MGECMUSRF &&tree->mod->whichmodel!=GYECMUSRF &&tree->mod->whichmodel!=GYECMUSR)
                {
                    fprintf(io->fp_out_compare,"%.6f;",tree->mod->omega_part[0]); //!< # omega
                    fprintf(io->fp_out_compare,"1;;"); //!< # prob omega and and empty space corresponding to alpha and beta
                }
                else
                {
                	if(io->mod->nparts > 1){printf("options not compatible with partitioned model error 1\n");exit(EXIT_FAILURE);}//Ken 22/8
                    fprintf(io->fp_out_compare,"%.6f;",Omega_ECMtoMmechModels(tree->mod->pi, tree->mod->qmat_part[0], tree->mod->qmat_buff_part[0], tree->mod->ns, tree->mod->n_w_catg)); //!< # omega//modified by Ken 22/8
                    fprintf(io->fp_out_compare,"1;;"); //!< # prob omega and and empty space corresponding to alpha and beta
                }
            }
            else
            {
                if(
                   tree->mod->whichmodel!=GYECMK07WK  && tree->mod->whichmodel!=GYECMK07WKF  && tree->mod->whichmodel!=GYECMS05WK  && tree->mod->whichmodel!=GYECMS05WKF &&
                   tree->mod->whichmodel!=MGECMK07WK  && tree->mod->whichmodel!=MGECMK07WKF  && tree->mod->whichmodel!=MGECMS05WK  && tree->mod->whichmodel!=MGECMS05WKF &&
                   tree->mod->whichmodel!=YAPECMK07WK && tree->mod->whichmodel!=YAPECMK07WKF && tree->mod->whichmodel!=YAPECMS05WK && tree->mod->whichmodel!=YAPECMS05WKF &&
                   tree->mod->whichmodel!=YAPECMUSRWKF &&tree->mod->whichmodel!=YAPECMUSRWK&&tree->mod->whichmodel!=MGECMUSRWK &&tree->mod->whichmodel!=MGECMUSRWKF &&tree->mod->whichmodel!=GYECMUSRWKF &&tree->mod->whichmodel!=GYECMUSRWK
                   )
                {
                    For(i,tree->mod->n_w_catg) fprintf(io->fp_out_compare,"%.6f ",tree->mod->omegas[i]);
                    fprintf(io->fp_out_compare,";");
                }
                else
                {
                	if(tree->mod->nparts > 1){printf("options not compatible with partitioned model error 2\n");exit(EXIT_FAILURE);}
                    For(i,tree->mod->n_w_catg) fprintf(io->fp_out_compare,"%.6f ",Omega_ECMtoMmechModels(tree->mod->pi, tree->mod->qmat_part[0]+i*tree->mod->ns*tree->mod->ns, tree->mod->qmat_buff_part[0]+i*tree->mod->ns*tree->mod->ns, tree->mod->ns,tree->mod->n_w_catg));
                    fprintf(io->fp_out_compare,";");
                }
                For(i,tree->mod->n_w_catg) fprintf(io->fp_out_compare,"%.6f ",tree->mod->prob_omegas[i]);
                fprintf(io->fp_out_compare,";");
                
                if(tree->mod->omegaSiteVar==DGAMMAK) fprintf(io->fp_out_compare,"%.6f %.6f;",tree->mod->alpha,tree->mod->beta);
                else fprintf(io->fp_out_compare,";");
            }
            
            fprintf(io->fp_out_compare,"**TIME**;");//!< Time elapsed //,(int)(t_end-t_beg));
            
            For(i,tree->mod->ns) fprintf(io->fp_out_compare,"%.6f ",tree->mod->pi[i]); //!< Codon frequencies.
            fprintf(io->fp_out_compare,";");
        }
        else
        {
            switch(tree->mod->s_opt->topo_search)
            {
                case SPR_MOVE: strcpy(m,"SPR\0");break;
                case NNI_MOVE: strcpy(m,"NNI\0");break;
                case BEST_OF_NNI_AND_SPR: strcpy(m,"BEST\0");break;
                default: break;
            }
            
            fprintf(io->fp_out_compare,"%d;",tree->n_otu); //!< Number of Taxa
            fprintf(io->fp_out_compare,"%d;",tree->data->init_len); //!< Sequence length
            fprintf(io->fp_out_compare,"%s;",io->mod->modelname); //!< Codon Model
            fprintf(io->fp_out_compare,";"); //!< Frequency model
            
            switch(io->eq_freq_handling)
            {
                case EMPIRICAL: fprintf(io->fp_out_compare,"empirical;");break;
                case OPTIMIZE:fprintf(io->fp_out_compare,"optimize;");break;
                case MODEL:fprintf(io->fp_out_compare,"model;");break;
                case USER:fprintf(io->fp_out_compare,"user;");break;
            }
            
            // 	if(io->datatype==NT)
            // 	{
            // 	  if(io->mod->s_opt->opt_state_freq) fprintf(io->fp_out_compare,"optimized;"); //!< Opt frequencies?
            // 	  else fprintf(io->fp_out_compare,"estimated;"); 
            // 	}
            // 	else
            // 	{
            // 	  if(io->mod->s_opt->opt_state_freq==NO && io->mod->s_opt->opt_state_freq_AAML==NO) fprintf(io->fp_out_compare,"empiric;"); //!< Opt frequencies?
            // 	  else if(io->mod->s_opt->opt_state_freq==YES && io->mod->s_opt->opt_state_freq_AAML==NO) fprintf(io->fp_out_compare,"estimated;");
            // 	  else if(io->mod->s_opt->opt_state_freq==YES && io->mod->s_opt->opt_state_freq_AAML==YES) fprintf(io->fp_out_compare,"optimized;"); 
            // 	}
            
            fprintf(io->fp_out_compare,"%s;",m); //!< SEARCH method
            fprintf(io->fp_out_compare,";"); //!< OPT method
            fprintf(io->fp_out_compare,";"); //!< EXPM method
            fprintf(io->fp_out_compare,"%.6f;",tree->c_lnL); //!< Tree Likelihood
            if(tree->io->datatype==AA)
            {
                if(tree->io->convert_NT_to_AA) fprintf(io->fp_out_compare,"%.6f %.6f %.6f %.6f;",tree->c_lnL+tree->logLk_correction, tree->c_lnL+tree->logLk_correction_gaps+tree->logLk_correction, tree->c_lnL+tree->logLk_correction_approx, tree->c_lnL+tree->logLk_correction_gaps_approx+tree->logLk_correction_approx); //!< comparable Tree Likelihood
                else fprintf(io->fp_out_compare,"%.6f %.6f;", tree->c_lnL+tree->logLk_correction_approx, tree->c_lnL+tree->logLk_correction_gaps_approx+tree->logLk_correction_approx);;
            }
            else fprintf(io->fp_out_compare,"%.6f;",tree->c_lnL); //!< comparable Tree Likelihood
            fprintf(io->fp_out_compare,"%.6f;",tree->size); //!< Tree length
            
            if(io->mod->n_catg>1) fprintf(io->fp_out_compare,"%d;",io->mod->n_catg); //!< Substitution rate categories
            else fprintf(io->fp_out_compare,"1;"); //!< Substitution rate categories
            
            if(io->mod->n_catg>1) fprintf(io->fp_out_compare,"%.3f;",io->mod->alpha); //!< Substitution rate categories Gamma param.
            else fprintf(io->fp_out_compare,";"); //!< Substitution rate categories Gamma param.
            
            fprintf(io->fp_out_compare,";"); //!< Kappa 
            fprintf(io->fp_out_compare,";"); //!< # Omega Model
            fprintf(io->fp_out_compare,";"); //!< # Omega 
            fprintf(io->fp_out_compare,";;"); //!< # Prob omega and and empty space corresponding to alpha and beta
            
            fprintf(io->fp_out_compare,"**TIME**;"); //!< Time elapsed,(int)(t_end-t_beg)); //!< Time elapsed
            
            For(i,tree->mod->ns) fprintf(io->fp_out_compare,"%.6f ",tree->mod->pi[i]); //!< Codon frequencies.
            fprintf(io->fp_out_compare,";");
        }
    }
    
    if(tree->c_lnL > UNLIKELY) {
	    currTkn = Emit_Out_Token(currTkn, "MLTree", "MLTree", SCALARTKN, TFSTRING, "%s", Write_Tree(tree));
	    currTkn = Emit_Out_Token(currTkn, "MLTreeSize", "MLTreeSize", SCALARTKN, TFNUMERIC, "%f", Get_Tree_Size(tree));
    }

    if(io->modeltypeOpt == HLP17){ //Added by Ken 23/8
    	 char *info = mCalloc(tree->n_pattern*2+1,sizeof(char));
    	 int indexi;
    	for(indexi=0;indexi<tree->n_pattern;indexi++){
    		char tempi[10];
        	sprintf(tempi, "%d ",io->mod->partIndex[indexi]);
        	strcat(info,tempi);
    	}
        currTkn = Emit_Out_Token(currTkn, "Omega partition indexes", "part", SCALARTKN, TFSTRING, "%s", info);
    }

    //Modified by Ken
    currTkn = Emit_Out_Token(currTkn, "Citation", "Citation", COMMENTTKN, TFSTRING,"%s", "");
    currTkn = Emit_Out_Token(currTkn, "Citation", "Citation", COMMENTTKN, TFSTRING, "%s", "If you use IgPhyML, please cite:");
    currTkn = Emit_Out_Token(currTkn, "Citation", "Citation", COMMENTTKN, TFSTRING,"%s", "K.B. Hoehn, G Lunter, O.G. Pybus");
    currTkn = Emit_Out_Token(currTkn, "Citation", "Citation", COMMENTTKN, TFSTRING,"%s", "A phylogenetic codon substitution model for antibody lineages");
    currTkn = Emit_Out_Token(currTkn, "Citation", "Citation", COMMENTTKN, TFSTRING,"%s", "Under review");
    currTkn = Emit_Out_Token(currTkn, "Citation", "Citation", COMMENTTKN, TFSTRING,"%s", "");
    currTkn = Emit_Out_Token(currTkn, "Citation", "Citation", COMMENTTKN, TFSTRING, "%s", "and");
    currTkn = Emit_Out_Token(currTkn, "Citation", "Citation", COMMENTTKN, TFSTRING,"%s", "");
    currTkn = Emit_Out_Token(currTkn, "Citation", "Citation", COMMENTTKN, TFSTRING,"%s", "M. Gil, M.S. Zanetti, S. Zoller and M. Anisimova 2013");
    currTkn = Emit_Out_Token(currTkn, "Citation", "Citation", COMMENTTKN, TFSTRING,"%s", "CodonPhyML: Fast Maximum Likelihood Phylogeny Estimation under Codon Substitution Models");
    currTkn = Emit_Out_Token(currTkn, "Citation", "Citation", COMMENTTKN, TFSTRING,"%s", "Molecular Biology and Evolution, pages 1270-1280, volume 30, number 6");
    currTkn = Emit_Out_Token(currTkn, "Citation", "Citation", COMMENTTKN, TFSTRING,"%s", "");
    currTkn = Emit_Out_Token(currTkn, "Citation", "Citation", COMMENTTKN, TFSTRING,"%s", "If you use aBayes branch supports please cite:");
    currTkn = Emit_Out_Token(currTkn, "Citation", "Citation", COMMENTTKN, TFSTRING,"%s", "M. Anisimova, M. Gil, J.F. Dufayard, C. Dessimoz and O. Gascuel 2011");
    currTkn = Emit_Out_Token(currTkn, "Citation", "Citation", COMMENTTKN, TFSTRING,"%s", "Survey of branch support methods demonstrates accuracy, power, and robustness of fast likelihood-based approximation schemes. Syst Biol 60:685-699");
    
    Output_Tokens(rootTkn, io->out_stats_format, io->fp_out_stats);
    
    free(s);//!< Added by MArcelo.
    free(r);//!< Added by MArcelo.
    free(t);//!< Added by MArcelo.

}

/*********************************************************/
/*FLT wrote this function*/
void Print_Fp_Out_Lines(FILE *fp_out, time_t t_beg, time_t t_end, t_tree *tree, option *io, int n_data_set)
{
  char *s;
  /*div_t hour,min;*/
  
  if (n_data_set==1)
  {
    
    PhyML_Fprintf(fp_out,". Sequence file : [%s]\n\n", Basename(io->in_align_file));
    
    if((tree->io->datatype == NT) || (tree->io->datatype == AA))
	  {
	    (tree->io->datatype == NT)?
      (PhyML_Fprintf(fp_out,". Model of nucleotides substitution : %s\n\n",io->mod->modelname)):
      (PhyML_Fprintf(fp_out,". Model of amino acids substitution : %s\n\n",io->mod->modelname));
	  }
    
    s = (char *)mCalloc(100,sizeof(char));
    
    switch(io->in_tree)
	  {
      case 0: { strcpy(s,"BioNJ");     break; }
      case 1: { strcpy(s,"parsimony"); break; }
      case 2: { strcpy(s,"user tree ("); 
        strcat(s,io->in_tree_file); 
        strcat(s,")");         break; }
	  }
    
    PhyML_Fprintf(fp_out,". Initial tree : [%s]\n\n",s);
    
    free(s);
    
    PhyML_Fprintf(fp_out,"\n");
    
    /*headline 1*/
    PhyML_Fprintf(fp_out, ". Data\t");
    
    PhyML_Fprintf(fp_out,"Nb of \t");
    
    PhyML_Fprintf(fp_out,"Likelihood\t");
    
    PhyML_Fprintf(fp_out, "Discrete   \t");
    
    if(tree->mod->n_catg > 1)
      PhyML_Fprintf(fp_out, "Number of \tGamma shape\t");
    
    PhyML_Fprintf(fp_out,"Proportion of\t");
    
    if(tree->mod->whichmodel <= 6)
      PhyML_Fprintf(fp_out,"Transition/ \t");
    
    PhyML_Fprintf(fp_out,"Nucleotides frequencies               \t");
    
    if((tree->mod->whichmodel == GTR) ||
       (tree->mod->whichmodel == CUSTOM))
      PhyML_Fprintf(fp_out,"Instantaneous rate matrix              \t");
    
    /*    PhyML_Fprintf(fp_out,"Time\t");*/
    
    PhyML_Fprintf(fp_out, "\n");
    
    
    /*headline 2*/
    PhyML_Fprintf(fp_out, "  set\t");
    
    PhyML_Fprintf(fp_out,"taxa\t");
    
    PhyML_Fprintf(fp_out,"loglk     \t");
    
    PhyML_Fprintf(fp_out, "gamma model\t");
    
    if(tree->mod->n_catg > 1)
      PhyML_Fprintf(fp_out, "categories\tparameter  \t");
    
    PhyML_Fprintf(fp_out,"invariant    \t");
    
    if(tree->mod->whichmodel <= 6)
      PhyML_Fprintf(fp_out,"transversion\t");
    
    PhyML_Fprintf(fp_out,"f(A)      f(C)      f(G)      f(T)    \t");
    
    if((tree->mod->whichmodel == GTR) ||
       (tree->mod->whichmodel == CUSTOM))
      PhyML_Fprintf(fp_out,"[A---------C---------G---------T------]\t");
    
    /*    PhyML_PhyML_Fprintf(fp_out,"used\t");*/
    
    PhyML_Fprintf(fp_out, "\n");
    
    
    /*headline 3*/
    if(tree->mod->whichmodel == TN93)
	  {
	    PhyML_Fprintf(fp_out,"    \t      \t          \t           \t");
	    if(tree->mod->n_catg > 1) PhyML_Fprintf(fp_out,"         \t         \t");
	    PhyML_Fprintf(fp_out,"             \t");
	    PhyML_Fprintf(fp_out,"purines pyrimid.\t");
	    PhyML_Fprintf(fp_out, "\n");
    }
    
    PhyML_Fprintf(fp_out, "\n");
  }
  
  
  /*line items*/
  
  PhyML_Fprintf(fp_out,"  #%d\t",n_data_set);
  
  PhyML_Fprintf(fp_out,"%d   \t",tree->n_otu);
  
  PhyML_Fprintf(fp_out,"%.5f\t",tree->c_lnL);
  
  PhyML_Fprintf(fp_out,"%s        \t",
                (tree->mod->n_catg>1)?("Yes"):("No "));
  if(tree->mod->n_catg > 1)
  {
    PhyML_Fprintf(fp_out,"%d        \t",tree->mod->n_catg);
    PhyML_Fprintf(fp_out,"%.3f    \t",tree->mod->alpha);
  }
  
  /*if(tree->mod->invar)*/
  PhyML_Fprintf(fp_out,"%.3f    \t",tree->mod->pinvar);
  
  if(tree->mod->whichmodel <= 5)
  {
    PhyML_Fprintf(fp_out,"%.3f     \t",tree->mod->kappa);
  }
  else if(tree->mod->whichmodel == TN93)
  {
    PhyML_Fprintf(fp_out,"%.3f   ",
                  tree->mod->kappa*2.*tree->mod->lambda/(1.+tree->mod->lambda));
    PhyML_Fprintf(fp_out,"%.3f\t",
                  tree->mod->kappa*2./(1.+tree->mod->lambda));
  }
  
  
  if(tree->io->datatype == NT)
  {
    PhyML_Fprintf(fp_out,"%8.5f  ",tree->mod->pi[0]);
    PhyML_Fprintf(fp_out,"%8.5f  ",tree->mod->pi[1]);
    PhyML_Fprintf(fp_out,"%8.5f  ",tree->mod->pi[2]);
    PhyML_Fprintf(fp_out,"%8.5f\t",tree->mod->pi[3]);
  }
  /*
   hour = div(t_end-t_beg,3600);
   min  = div(t_end-t_beg,60  );
   
   min.quot -= hour.quot*60;
   
   PhyML_Fprintf(fp_out,"%dh%dm%ds\t", hour.quot,min.quot,(int)(t_end-t_beg)%60);
   if(t_end-t_beg > 60)
   PhyML_Fprintf(fp_out,". -> %d seconds\t",(int)(t_end-t_beg));
   */
  
  /*****************************************/
  if((tree->mod->whichmodel == GTR) || (tree->mod->whichmodel == CUSTOM))
  {
    int i,j;
    
    For(i,4)
    {
      if (i!=0) {
        /*format*/
        PhyML_Fprintf(fp_out,"      \t     \t          \t           \t");
        if(tree->mod->n_catg > 1) PhyML_Fprintf(fp_out,"          \t           \t");
        PhyML_Fprintf(fp_out,"             \t                                      \t");
      }
      For(j,4)
	    PhyML_Fprintf(fp_out,"%8.5f  ",tree->mod->qmat[i*4+j]);
      if (i<3) PhyML_Fprintf(fp_out,"\n");
    }
  }
  /*****************************************/
  
  PhyML_Fprintf(fp_out, "\n\n");
}

/*********************************************************/
matrix *JC69_Dist_Codon(calign *data, model *mod)
{
  int site,i,j,k;
  phydbl unc_len;
  matrix *mat;
  phydbl **len;
  int datatype;
  
  
  len = (phydbl **)mCalloc(data->n_otu,sizeof(phydbl *));
  For(i,data->n_otu)
  len[i] = (phydbl *)mCalloc(data->n_otu,sizeof(phydbl));
  
  unc_len = .0;
  
  mat = Make_Mat(data->n_otu);
  Init_Mat(mat,data);
  
  datatype = mod->io->datatype;
  
  For(i,data->n_otu) For(j,data->n_otu) {mat->P[i][j]=0; len[i][j]=0;}
  
  
  For(site,data->c_seq[0]->len)
  {
    For(j,data->n_otu-1)
    {
      for(k=j+1;k<data->n_otu;k++)
	    {
	      if((!Is_Ambigu(data->c_seq[j]->state+site*mod->io->mod->state_len,datatype,mod->io->mod->state_len)) &&
           (!Is_Ambigu(data->c_seq[k]->state+site*mod->io->mod->state_len,datatype,mod->io->mod->state_len)))
        {
          len[j][k]+=data->wght[site];
          len[k][j]=len[j][k];
          
          
          if(data->c_seq[j]->state+site != data->c_seq[k]->state+site)	mat->P[j][k] += data->wght[site];
        }
	    }
    }
  }
  
  
  For(i,data->n_otu-1)
  for(j=i+1;j<data->n_otu;j++)
  {
    if(len[i][j] > .0) mat->P[i][j] /= len[i][j];
    else               mat->P[i][j] = 1.;
    
    mat->P[j][i] = mat->P[i][j];
    
    if((1.-(mod->ns)/(mod->ns-1.)*mat->P[i][j]) < .0) mat->dist[i][j] = -1.;
    else
      mat->dist[i][j] = -(mod->ns-1.)/(mod->ns)*(phydbl)LOG(1.-(mod->ns)/(mod->ns-1.)*mat->P[i][j]);
    
    /* 	PhyML_Printf("\n. Incorrect JC distances"); */
    /* 	mat->dist[i][j] = len[i][j]; */
    
    if(mat->dist[i][j] > DIST_MAX) mat->dist[i][j] = DIST_MAX;
    
    mat->dist[j][i] = mat->dist[i][j];
  }
  
  For(i,data->n_otu) free(len[i]);
  free(len);
  
  return mat;
}

/*********************************************************/
// If state is (char)88, mark as ambiguous
int Is_Ambigu(char *state, int datatype, int stepsize)
{
  int val,i;
  
  val = -1;
  if(datatype == NT)
  {
    For(i,stepsize)
    {
      switch(state[i])
	    {
        case 'A' : case 'C' : case 'G' : case 'T' : case 'U' : { val=0; break; }
        default : { val=1; break; }
	    }
      if(val == 1) break;
    }
  }
  else if(datatype == AA)
  {
    switch(state[0])
    {
      case 'X' : case '?' : case '-' : case '.' : {val=1; break; }
      default : { val=0; break; }
    }
  }
  else if(datatype == GENERIC)
  {
    int i;
    For(i,stepsize) if(!isdigit(state[i])) break;
    if(i == stepsize) val = 0;
    else              val = 1;
  }
  else if(datatype == CODON)
  {
    if( (state[0]>=(char)0) && (state[0]<(char)64) ) val=0;
    else if(state[0]==(char)88) val=1;
    else{printf("Not allowed codon %c.\n",state[0]); Warn_And_Exit("Impossible to continue.\n");}
  }
  return val;
}

/*********************************************************/

void Check_Ambiguities(calign *data, int datatype, int stepsize)
{
  int i,j;
  
  For(j,data->crunch_len) 
  {
    For(i,data->n_otu)
    {
      data->ambigu[j]                  = 0;
      data->c_seq[i]->is_ambigu[j]     = 0;
    }
    
    For(i,data->n_otu)
    {
      if(Is_Ambigu(data->c_seq[i]->state+j*stepsize,
                   datatype,
                   stepsize))
	    {
	      data->ambigu[j]              = 1;
	      data->c_seq[i]->is_ambigu[j] = 1;
	    }
    }
  }
}

/*********************************************************/

int Assign_State(char *c, int datatype, int stepsize)
{
  int state[3];
  int i;
  
  state[0] = state[1] = state[2] = -1;
  if(datatype == NT)
  {
    For(i,stepsize)
    {
      switch(c[i])
      {
        case 'A' : {state[i]=0;  break;}
        case 'C' : {state[i]=1;  break;}
        case 'G' : {state[i]=2;  break;}
        case 'T' : {state[i]=3;  break;}
        case 'U' : {state[i]=3;  break;}
        default  : {state[i]=-1; break;}
      }
    }
    return (stepsize>1)?(state[0]*16+state[1]*4+state[2]):(state[0]);
	}
	else if(datatype == AA)
	{
	  switch(c[0])
	  {
	    case 'A' : {state[0]=0 ; break;}
	    case 'R' : {state[0]=1 ; break;}
	    case 'N' : {state[0]=2 ; break;}
	    case 'D' : {state[0]=3 ; break;}
	    case 'C' : {state[0]=4 ; break;}
	    case 'Q' : {state[0]=5 ; break;}
	    case 'E' : {state[0]=6 ; break;}
	    case 'G' : {state[0]=7 ; break;}
	    case 'H' : {state[0]=8 ; break;}
	    case 'I' : {state[0]=9 ; break;}
	    case 'L' : {state[0]=10; break;}
	    case 'K' : {state[0]=11; break;}
	    case 'M' : {state[0]=12; break;}
	    case 'F' : {state[0]=13; break;}
	    case 'P' : {state[0]=14; break;}
	    case 'S' : {state[0]=15; break;}
	    case 'T' : {state[0]=16; break;}
	    case 'W' : {state[0]=17; break;}
	    case 'Y' : {state[0]=18; break;}
	    case 'V' : {state[0]=19; break;}
        
	    case 'B' : {state[0] = 2; break;}
	    case 'Z' : {state[0] = 5; break;}
	    default  : {state[0]=-1;  break;}
	  }
	  return state[0];
	}
	else if(datatype == GENERIC)
	{
	  char format[6];
	  int ret;
	  
	  sprintf(format,"%%%dd",stepsize);
	  ret = sscanf(c,format,state);
	  if(!ret) state[0] = -1;      
	  return state[0];
	}
	else if(datatype==CODON)
	{
	  if(*c>=(char)0 && *c<(char)64) state[0]=(int)*c; //=*c;
	  else if (*c==(char)88) state[0]=-1;//(char)88;
	  else 
	  {
	    PhyML_Printf("\n. Unknown character state : '%c'\n",*c);
	    Warn_And_Exit("\n. Init failed (check the data type)\n");
	  }
	  return (int)state[0];
	}  
	else
	{
	  PhyML_Printf("\n. Not implemented yet.\n");
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
	
	return -1;
}

/*********************************************************/

int Assign_State_With_Ambiguity(char *c, int datatype, int stepsize)
{
  int state[3];
  int i;
  
  state[0] = state[1] = state[2] = -1;
  if(datatype == NT)
  {
    For(i,stepsize)
    {
      switch(c[i])
	    {
        case 'A' : {state[i]= 0;  break;}
        case 'C' : {state[i]= 1;  break;}
        case 'G' : {state[i]= 2;  break;}
        case 'T' : {state[i]= 3;  break;}
        case 'U' : {state[i]= 3;  break;}
        case 'M' : {state[i]= 4;  break;}
        case 'R' : {state[i]= 5;  break;}
        case 'W' : {state[i]= 6;  break;}
        case 'S' : {state[i]= 7;  break;}
        case 'Y' : {state[i]= 8;  break;}
        case 'K' : {state[i]= 9;  break;}
        case 'B' : {state[i]=10;  break;}
        case 'D' : {state[i]=11;  break;}
        case 'H' : {state[i]=12;  break;}
        case 'V' : {state[i]=13;  break;}
        case 'N' : case 'X' : case '?' : case 'O' : case '-' : {state[i]=T_MAX_ALPHABET-1;  break;}
        default :
	      {
          PhyML_Printf("\n. Unknown character state : '%c'\n",c[i]);
          Warn_And_Exit("\n. Init failed (check the data type)\n");
          break;
	      }
	    }
      return (stepsize>1)?(state[0]*16+state[1]*4+state[2]):(state[0]);
    }
  }
  else if(datatype == AA)
  {
    switch(c[0])
    {
      case 'A' : {state[0]= 0; break;}
      case 'R' : {state[0]= 1; break;}
      case 'N' : {state[0]= 2; break;}
      case 'D' : {state[0]= 3; break;}
      case 'C' : {state[0]= 4; break;}
      case 'Q' : {state[0]= 5; break;}
      case 'E' : {state[0]= 6; break;}
      case 'G' : {state[0]= 7; break;}
      case 'H' : {state[0]= 8; break;}
      case 'I' : {state[0]= 9; break;}
      case 'L' : {state[0]=10; break;}
      case 'K' : {state[0]=11; break;}
      case 'M' : {state[0]=12; break;}
      case 'F' : {state[0]=13; break;}
      case 'P' : {state[0]=14; break;}
      case 'S' : {state[0]=15; break;}
      case 'T' : {state[0]=16; break;}
      case 'W' : {state[0]=17; break;}
      case 'Y' : {state[0]=18; break;}
      case 'V' : {state[0]=19; break;}
      case 'B' : {state[0]= 2; break;}
      case 'Z' : {state[0]= 5; break;}
      case 'X' : case '?' : case '-' : {state[0]=T_MAX_ALPHABET-1; break;}
      default  : 
      {
        PhyML_Printf("\n. Unknown character state : %c\n",state[0]);
        Warn_And_Exit("\n. Init failed (check the data type)\n");
        break;
      }
    }
    return state[0];
  }
  else if(datatype == GENERIC)
  {
    if(Is_Ambigu(c,GENERIC,stepsize)) state[0] = T_MAX_ALPHABET-1;
    else
    {
      char format[6];
      sprintf(format,"%%%dd",stepsize);
      if(!sscanf(c,format,state))
	    {
	      PhyML_Printf("\n. Error reading character. Was expecting an integer, got '%c' instead.\n",c[0]);
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
    }
    return state[0];
  }
  else if(datatype==CODON)
  {
    if(*c>=(char)0 && *c<(char)64) state[0]=*c;
    else if (*c==(char)88) state[0]=(char)88;
    else 
    {
      PhyML_Printf("\n. Unknown character state : '%c'\n",*c);
      Warn_And_Exit("\n. Init failed (check the data type)\n");
    }
    return (int)state[0];
  }  
  
  return -1;
}
/*********************************************************/



void Clean_Tree_Connections(t_tree *tree)
{
  
  int i;
  For(i,2*tree->n_otu-2)
  {
    tree->noeud[i]->v[0] = NULL;
    tree->noeud[i]->v[1] = NULL;
    tree->noeud[i]->v[2] = NULL;
    tree->noeud[i]->b[0] = NULL;
    tree->noeud[i]->b[1] = NULL;
    tree->noeud[i]->b[2] = NULL;
  }
}

/*********************************************************/

void Bootstrap(t_tree *tree)
{
  int *site_num, n_site;
  int replicate,j,k;
  int position,init_len;
  calign *boot_data;
  t_tree *boot_tree;
  model *boot_mod;
  matrix *boot_mat;
  char *s;
  /*   phydbl rf; */
  
  tree->print_boot_val = 1;
  tree->print_alrt_val = 0;
  boot_tree            = NULL;
  
  site_num = (int *)mCalloc(tree->data->init_len,sizeof(int));
  
  Free_Bip(tree);
  Alloc_Bip(tree);
  Get_Bip(tree->noeud[0],tree->noeud[0]->v[0],tree);
  
  n_site = 0;
  For(j,tree->data->crunch_len) For(k,tree->data->wght[j])
  {
    site_num[n_site] = j;
    n_site++;
  }
  
  boot_data = Copy_Cseq(tree->data,tree->io);
  
  PhyML_Printf("\n\n. Non parametric bootstrap analysis \n\n");
  PhyML_Printf("  ["); 
  
  For(replicate,tree->mod->bootstrap)
  {
    For(j,boot_data->crunch_len) boot_data->wght[j] = 0;
    
    init_len = 0;
    For(j,boot_data->init_len)
    {
      position = Rand_Int(0,(int)(tree->data->init_len-1.0));
      boot_data->wght[site_num[position]] += 1;
      init_len++;
    }
    
    if(init_len != tree->data->init_len) Warn_And_Exit("\n. Pb when copying sequences\n");
    
    init_len = 0;
    For(j,boot_data->crunch_len) init_len += boot_data->wght[j];
    
    if(init_len != tree->data->init_len) Warn_And_Exit("\n. Pb when copying sequences\n");
    
    
    if(tree->io->random_boot_seq_order) Randomize_Sequence_Order(boot_data);
    
    boot_mod        = Copy_Model(tree->mod);
    boot_mod->s_opt = tree->mod->s_opt; /* WARNING: re-using the same address here instead of creating a copying
                                         requires to leave the value of s_opt unchanged during the boostrap. */
    boot_mod->io    = tree->io; /* WARNING: re-using the same address here instead of creating a copying
                                 requires to leave the value of io unchanged during the boostrap. */
    Init_Model(boot_data,boot_mod,tree->io);
    
    if(tree->io->in_tree == 2)
    {
      rewind(tree->io->fp_in_tree);
      boot_tree = Read_Tree_File_Phylip(tree->io->fp_in_tree);
    }
    else
    {
      boot_mat = ML_Dist(boot_data,boot_mod);
      boot_mat->tree = Make_Tree_From_Scratch(boot_data->n_otu,boot_data);
      Fill_Missing_Dist(boot_mat);
      Bionj(boot_mat);
      boot_tree = boot_mat->tree;
      boot_tree->mat = boot_mat;
    }
    
    
    boot_tree->mod                = boot_mod;
    boot_tree->io                 = tree->io;
    boot_tree->data               = boot_data;
    boot_tree->both_sides         = 1;
    boot_tree->mod->s_opt->print  = 0;
    boot_tree->n_pattern          = boot_tree->data->crunch_len;
    boot_tree->io->print_site_lnl = 0;
    boot_tree->io->print_trace    = 0;
    
    if((boot_tree->mod->s_opt->random_input_tree) && (boot_tree->mod->s_opt->topo_search == SPR_MOVE)) Random_Tree(boot_tree);
    Order_Tree_CSeq(boot_tree,boot_data);
    Share_Lk_Struct(tree,boot_tree);
    Share_Spr_Struct(tree,boot_tree);
    Share_Pars_Struct(tree,boot_tree);
    Fill_Dir_Table(boot_tree);
    Update_Dirs(boot_tree);
    
    if(tree->mod->s_opt->greedy) Init_P_Lk_Tips_Double(boot_tree);
    else                         Init_P_Lk_Tips_Int(boot_tree);
    Init_Ui_Tips(boot_tree);
    Init_P_Pars_Tips(boot_tree);
    Br_Len_Not_Involving_Invar(boot_tree);
    
    if(boot_tree->mod->s_opt->opt_topo)
    {
      if(boot_tree->mod->s_opt->topo_search == NNI_MOVE) 
	    {
	      Simu_Loop(boot_tree);
	    }
      else if((boot_tree->mod->s_opt->topo_search == SPR_MOVE) ||
              (boot_tree->mod->s_opt->topo_search == BEST_OF_NNI_AND_SPR))
	    {
	      Speed_Spr_Loop(boot_tree);
	    }
    }
    else
    {
      if(boot_tree->mod->s_opt->opt_subst_param || boot_tree->mod->s_opt->opt_bl)
        Round_Optimize(boot_tree,boot_tree->data,ROUND_MAX);
      else
        Lk(boot_tree);
    }
    
    Free_Bip(boot_tree);
    
    Alloc_Bip(boot_tree);
    
    Match_Tip_Numbers(tree,boot_tree);
    
    Get_Bip(boot_tree->noeud[0],
            boot_tree->noeud[0]->v[0],
            boot_tree);
    
    
    Compare_Bip(tree,boot_tree);
    
    Br_Len_Involving_Invar(boot_tree);
    
    
    if(tree->io->print_boot_trees)
    {
      s = Write_Tree(boot_tree);
      PhyML_Fprintf(tree->io->fp_out_boot_tree,"%s\n",s);
      free(s);
      Print_Fp_Out_Lines(tree->io->fp_out_boot_stats,0,0,boot_tree,tree->io,replicate+1);
    }
    
    /*       rf = .0; */
    /*       For(j,2*tree->n_otu-3)  */
    /* 	rf += tree->t_edges[j]->bip_score; */
    
    
    PhyML_Printf("."); 
#ifndef QUIET
    fflush(stdout);
#endif
    if(!((replicate+1)%20))
    {
      PhyML_Printf("] %4d/%4d\n  ",replicate+1,tree->mod->bootstrap);
      if(replicate != tree->mod->bootstrap-1) PhyML_Printf("[");
    }
    
    if(boot_tree->mat) Free_Mat(boot_tree->mat);
    Free_Tree(boot_tree);      
    Free_Model(boot_mod);
  }
  
  if(((replicate)%20)) PhyML_Printf("] %4d/%4d\n ",replicate,tree->mod->bootstrap);
  
  tree->lock_topo = 1; /* Topology should not be modified afterwards */
  
  if(tree->io->print_boot_trees)
  {
    fclose(tree->io->fp_out_boot_tree);
    fclose(tree->io->fp_out_boot_stats);
  }
  
  Free_Cseq(boot_data);
  free(site_num);
}

/*********************************************************/

void Br_Len_Involving_Invar(t_tree *tree)
{
  if(!tree->br_len_invar_set) {
    int i;
    For(i,2*tree->n_otu-3) tree->t_edges[i]->l *= (1.0-tree->mod->pinvar);
    tree->br_len_invar_set = YES;
  }
}

/*********************************************************/

void Br_Len_Not_Involving_Invar(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-3) tree->t_edges[i]->l /= (1.0-tree->mod->pinvar);
}

/*********************************************************/

void Getstring_Stdin(char *file_name)
{
  if(!fgets(file_name,T_MAX_LINE,stdin)) Exit("");
  if (strchr(file_name, '\n') != NULL)
    *strchr(file_name, '\n') = '\0';
}

/*********************************************************/
/*********************************************************/

void Copy_One_State(char *from, char *to, int state_size)
{
  int i;
  For(i,state_size) to[i] = from[i];
}

/*********************************************************/

//! Allocate memory for basic model options.

/*!
 *
 */
model *Make_Model_Basic()
{
  model *mod;
  
  mod                     = (model *)mCalloc(1,sizeof(model));              /*!< Allocate memory for a model structure.*/
  mod->modelname          = (char *)mCalloc(T_MAX_NAME,sizeof(char));       /*!< Allocate memory for string for the name of the model.*/
  mod->custom_mod_string  = (char *)mCalloc(T_MAX_OPTION,sizeof(char));     /*!< Allocate memory for string to store user custom model.*/
  mod->user_b_freq        = (phydbl *)mCalloc(T_MAX_ALPHABET,sizeof(phydbl)); /*!< Allocate memory to store the values of equilibrium frequencies provided by the user.*/
  
  mod->pkappa=(phydbl *)mCalloc(9,sizeof(phydbl));    //!< Added by Marcelo.
  mod->unspkappa=(phydbl *)mCalloc(9,sizeof(phydbl)); //!< Added by Marcelo.
  
  return mod;
}

/*********************************************************/

void Make_Model_Complete(model *mod) 
{
  mod->gamma_r_proba                         = (phydbl *)mCalloc(mod->n_catg,sizeof(phydbl));
  mod->gamma_rr                              = (phydbl *)mCalloc(mod->n_catg,sizeof(phydbl));
  mod->Pij_rr                                = (phydbl *)mCalloc(mod->n_catg*mod->ns*mod->ns,sizeof(phydbl));
  mod->eigen                                 = (eigen  *)Make_Eigen_Struct(mod);       
  mod->pi_unscaled                           = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));
  mod->pi                                    = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));

  //added by Ken 17/8/2016
  mod->qmat_part = (phydbl **)mCalloc(mod->nparts,sizeof(phydbl*));
  mod->Pmat_part = (phydbl **)mCalloc(mod->nparts,sizeof(phydbl*));
  mod->qmat_buff_part = (phydbl **)mCalloc(mod->nparts,sizeof(phydbl*));

	  int i;
	  for(i=0;i<mod->nparts;i++){
		  mod->qmat_part[i]=(phydbl *)mCalloc(mod->n_w_catg*mod->ns*mod->ns,sizeof(phydbl));
		  mod->Pmat_part[i]=(phydbl *)mCalloc(mod->n_catg*mod->ns*mod->ns,sizeof(phydbl));
		  mod->qmat_buff_part[i] = (phydbl *)mCalloc(mod->n_w_catg*mod->ns*mod->ns,sizeof(phydbl));
	  }


  if(mod->n_rr_branch)
  {
    mod->rr_branch                           = (phydbl *)mCalloc(mod->n_rr_branch,sizeof(phydbl));
    mod->p_rr_branch                         = (phydbl *)mCalloc(mod->n_rr_branch,sizeof(phydbl));
  }
  
  if(mod->io->datatype==CODON) //!< Added by Marcelo. 
  {
    int i,j;
    mod->prob_omegas_uns                     = (phydbl *)mCalloc(mod->n_w_catg,sizeof(phydbl));   
    mod->mr_w                                = (phydbl *)mCalloc(mod->n_w_catg,sizeof(phydbl));   
    mod->base_freq                           = (phydbl *)mCalloc(mod->num_base_freq,sizeof(phydbl));   
    mod->uns_base_freq                       = (phydbl *)mCalloc(mod->num_base_freq,sizeof(phydbl));   
    
    mod->A0_part                                  = (phydbl **)mCalloc(mod->nomega_part,sizeof(phydbl*));
    int modeli = 0;
    for(modeli=0;modeli<mod->nomega_part;modeli++){
    	mod->A0_part[modeli]                                  = (phydbl *)mCalloc(mod->ns*mod->ns,sizeof(phydbl)); //modified by Ken 22/8
        Fors(i,mod->ns*mod->ns,mod->ns+1) mod->A0_part[modeli][i]=1.0; //!< create a identity matrix.
    }

    mod->qmatScaled                          = (phydbl *)mCalloc(mod->n_w_catg*mod->ns*mod->ns,sizeof(phydbl));
    
    if(mod->io->heuristicExpm)
    {
    mod->A2_part                                  = (phydbl **)mCalloc(mod->nomega_part,sizeof(phydbl*));
     for(modeli=0;modeli<mod->nomega_part;modeli++){
    	mod->A2_part[modeli]                               = (phydbl *)mCalloc(15*mod->ns*mod->ns*mod->n_w_catg,sizeof(phydbl));
     }
      //Fors(i,mod->ns*mod->ns*mod->n_w_catg*15,mod->ns*mod->ns*15) For(j,mod->ns*mod->ns) mod->A2[i+j]=mod->A0[j];
    }
    
    if(mod->io->expm==SSPADE)
    {

      mod->ipiv_part                                = (int    **)mCalloc(mod->nomega_part,sizeof(int*));
      mod->U_part                                   = (phydbl **)mCalloc(mod->nomega_part,sizeof(phydbl*));
      mod->V_part                                   = (phydbl **)mCalloc(mod->nomega_part,sizeof(phydbl*));
      mod->A2_part                                  = (phydbl **)mCalloc(mod->nomega_part,sizeof(phydbl*));
      mod->A4_part                                  = (phydbl **)mCalloc(mod->nomega_part,sizeof(phydbl*));
      mod->A6_part                                  = (phydbl **)mCalloc(mod->nomega_part,sizeof(phydbl*));
      mod->A8_part                                  = (phydbl **)mCalloc(mod->nomega_part,sizeof(phydbl*));
      mod->Apowers_part                             = (phydbl **)mCalloc(mod->nomega_part,sizeof(phydbl*));
      mod->matAux_part                              = (phydbl **)mCalloc(mod->nomega_part,sizeof(phydbl*));
      for(modeli=0;modeli<mod->nomega_part;modeli++){
    	        mod->ipiv_part[modeli]                                = (int    *)mCalloc(mod->ns*mod->n_w_catg,sizeof(int));
    	        mod->U_part[modeli]                                   = (phydbl *)mCalloc(mod->ns*mod->ns*mod->n_w_catg,sizeof(phydbl));
    	        mod->V_part[modeli]                                   = (phydbl *)mCalloc(mod->ns*mod->ns*mod->n_w_catg,sizeof(phydbl));
    	        mod->A2_part[modeli]                                  = (phydbl *)mCalloc(mod->ns*mod->ns*mod->n_w_catg,sizeof(phydbl));
    	        mod->A4_part[modeli]                                  = (phydbl *)mCalloc(mod->ns*mod->ns*mod->n_w_catg,sizeof(phydbl));
    	        mod->A6_part[modeli]                                  = (phydbl *)mCalloc(mod->ns*mod->ns*mod->n_w_catg,sizeof(phydbl));
    	        mod->A8_part[modeli]                                  = (phydbl *)mCalloc(mod->ns*mod->ns*mod->n_w_catg,sizeof(phydbl));
    	        mod->Apowers_part[modeli]                             = (phydbl *)mCalloc(mod->ns*mod->ns*mod->n_w_catg*5,sizeof(phydbl));
    	        mod->matAux_part[modeli]                              = (phydbl *)mCalloc(mod->ns*mod->ns*mod->n_w_catg,sizeof(phydbl));
    	        Fors(i,mod->ns*mod->ns*mod->n_w_catg*5,mod->ns*mod->ns*5) For(j,mod->ns*mod->ns) mod->Apowers_part[modeli][i+j]=mod->A0_part[modeli][j];
      }


    }
  }
}

/*********************************************************/

void Copy_Dist(phydbl **cpy, phydbl **orig, int n)
{
  int i,j;
  For(i,n) For(j,n) cpy[i][j] = orig[i][j];
}

/*********************************************************/

model *Copy_Model(model *ori)
{
  model *cpy;
  
  cpy                = Make_Model_Basic();
  cpy->ns            = ori->ns;
  cpy->n_catg        = ori->n_catg;
  
  cpy->num_base_freq = ori->num_base_freq;    //!< Added by Marcelo.
  cpy->io            = ori->io;               //!< Added by Marcelo.
  cpy->s_opt         = ori->s_opt;            //!< Added by Marcelo.
  cpy->n_w_catg      = ori->n_w_catg;         //!< Added by Marcelo.
  
  Make_Model_Complete(cpy);
  
  
  if(ori->io->datatype==CODON) //!< Added by Marcelo.
  {
    cpy->prob_omegas = (phydbl *)mCalloc(cpy->n_w_catg,sizeof(phydbl));
    cpy->omegas      = (phydbl *)mCalloc(cpy->n_w_catg,sizeof(phydbl));
  }
  
  Record_Model(ori,cpy);
  
#ifdef M4
  if(ori->m4mod) cpy->m4mod = M4_Copy_M4_Model(ori, ori->m4mod);
#endif
  
  return cpy;
}

/*********************************************************/

void Record_Model(model *ori, model *cpy)
{
  int i;
  if(ori->nparts > 1){
	  printf("Record_Model doesn't work with parititions yet\n");
	  exit(EXIT_FAILURE);
  }
  
  cpy->alpha_old    = ori->alpha_old;
  cpy->kappa_old    = ori->alpha_old;
  cpy->lambda_old   = ori->lambda_old;
  cpy->pinvar_old   = ori->pinvar_old;
  cpy->whichmodel   = ori->whichmodel;
  cpy->update_eigen = ori->update_eigen;
  cpy->kappa        = ori->kappa;
  cpy->alpha        = ori->alpha;
  cpy->lambda       = ori->lambda;
  cpy->bootstrap    = ori->bootstrap;
  cpy->invar        = ori->invar;
  cpy->pinvar       = ori->pinvar;
  cpy->n_diff_rr    = ori->n_diff_rr;
  
  if((ori->whichmodel == CUSTOM)||(ori->whichmodel == GTR))
  {
    For(i,ori->ns*(ori->ns-1)/2)
    {
      cpy->rr_num[i]      = ori->rr_num[i];
      cpy->rr_val[i]      = ori->rr_val[i];
      cpy->rr[i]          = ori->rr[i]; //!< Corrected by Marcelo.
      cpy->n_rr_per_cat[i]= ori->n_rr_per_cat[i];
    }
  }
  
  For(i,cpy->ns)
  {
    cpy->pi[i]          = ori->pi[i];
    cpy->pi_unscaled[i] = ori->pi_unscaled[i];
    cpy->user_b_freq[i] = ori->user_b_freq[i];
  }
  
  For(i,ori->n_w_catg*cpy->ns*cpy->ns)
  { 
    cpy->qmat[i]      = ori->qmat[i]; 
    cpy->qmat_buff_part[0][i] = ori->qmat_buff_part[0][i];
  }
  
  For(i,ori->n_w_catg)
  {
    cpy->eigen->size = ori->eigen->size;
    For(i,2*ori->ns*ori->n_w_catg)       cpy->eigen->space[i]       = ori->eigen->space[i];
    For(i,2*ori->ns*ori->n_w_catg)       cpy->eigen->space_int[i]   = ori->eigen->space_int[i];
    For(i,ori->ns*ori->n_w_catg)         cpy->eigen->e_val[i]       = ori->eigen->e_val[i];
    For(i,ori->ns*ori->n_w_catg)         cpy->eigen->e_val_im[i]    = ori->eigen->e_val_im[i];
    For(i,ori->ns*ori->ns*ori->n_w_catg) cpy->eigen->r_e_vect[i]    = ori->eigen->r_e_vect[i];
    For(i,ori->ns*ori->ns*ori->n_w_catg) cpy->eigen->r_e_vect[i]    = ori->eigen->r_e_vect[i];
    For(i,ori->ns*ori->ns*ori->n_w_catg) cpy->eigen->r_e_vect_im[i] = ori->eigen->r_e_vect_im[i];
    For(i,ori->ns*ori->ns*ori->n_w_catg) cpy->eigen->l_e_vect[i]    = ori->eigen->l_e_vect[i];
    For(i,ori->ns*ori->ns*ori->n_w_catg) cpy->eigen->q[i]           = ori->eigen->q[i];
  }
  
  For(i,cpy->n_catg)
  {
    cpy->gamma_r_proba[i] = ori->gamma_r_proba[i];
    cpy->gamma_rr[i]      = ori->gamma_rr[i];
  }
  
#ifndef PHYML
  cpy->use_m4mod = ori->use_m4mod;
#endif 
  
  if(ori->io->datatype==CODON) //!< Added by Marcelo.
  {
    int j;
    cpy->omegaSiteVar = ori->omegaSiteVar; 
    //cpy->omega        = ori->omega;

    int omegai; //added by Ken 17/8/2016
    for(omegai=0;omegai<ori->nomega_part;omegai++){
    	cpy->omega_part[omegai]    = ori->omega_part[omegai];
    }
    cpy->nomega_part = ori->nomega_part;

    cpy->omega_old    = ori->omega_old; 
    cpy->freq_model   = ori->freq_model;
    cpy->genetic_code = ori->genetic_code;
    strcpy(cpy->modelname,ori->modelname);
    
    For(i,cpy->num_base_freq)
    { 
      cpy->base_freq[i]     = ori->base_freq[i];  
      cpy->uns_base_freq[i] = ori->uns_base_freq[i];
    }
    
    cpy->nkappa             = ori->nkappa;
    For(i,ori->nkappa)
    { 
      cpy->pkappa[i]        = ori->pkappa[i];
      cpy->unspkappa[i]     = ori->unspkappa[i];   
    }
    
    For(j,ori->n_w_catg)
    {
      cpy->mr_w[j]              = ori->mr_w[j];
      cpy->prob_omegas_uns[j]   = ori->prob_omegas_uns[j];
      cpy->omegas[i]            = ori->omegas[i];
      cpy->prob_omegas[i]       = ori->prob_omegas[i];
    }
  }
}
/*********************************************************/
//! Allocate memory for pointers to strings that contain the names of input and output files.
/*!
 *  Input: files of sequence alignement, initial tree file (if any ); Output: tree file, statistic file, bootstrap file, likelihood file, clade list.
 */  
option *Make_Input()
{
  int i;
  option* io               = (option *)mCalloc(1,sizeof(option)); 
  
  io->in_align_file        = (char *)mCalloc(T_MAX_FILE,sizeof(char)); //!< Name of input aligned sequences file.
  io->in_tree_file         = (char *)mCalloc(T_MAX_FILE,sizeof(char)); //!< Name of input tree file.
  io->out_tree_file        = (char *)mCalloc(T_MAX_FILE,sizeof(char)); //!< Name of output tree file.
  io->out_trees_file       = (char *)mCalloc(T_MAX_FILE,sizeof(char)); //!< Name of output tree file.
  io->out_boot_tree_file   = (char *)mCalloc(T_MAX_FILE,sizeof(char)); //!< Name of output boostrapped tree file.
  io->out_boot_stats_file  = (char *)mCalloc(T_MAX_FILE,sizeof(char)); //!< Name of output boostrap statistics file.
  io->out_stats_file       = (char *)mCalloc(T_MAX_FILE,sizeof(char)); //!< Name of output statistics file.
  io->out_lk_file          = (char *)mCalloc(T_MAX_FILE,sizeof(char)); //!< Name of output likelihood file.
  io->out_ps_file          = (char *)mCalloc(T_MAX_FILE,sizeof(char)); //!< Name of output parsimony score file.
 // io->out_trace_file       = (char *)mCalloc(T_MAX_FILE,sizeof(char)); //!< Name of output file containing each phylogeny explored during tree search.
  io->out_trace_tree_file  = (char *)mCalloc(T_MAX_FILE,sizeof(char)); //!< Name of output file containing each phylogeny explored during tree search.  //added by Ken 9/2/2017
  io->out_trace_stats_file = (char *)mCalloc(T_MAX_FILE,sizeof(char)); //!< Name of output file containing each phylogeny explored during tree search. //added by Ken 9/2/2017
  io->tracecount = 0;
  io->nt_or_cd             = (char *)mCalloc(T_MAX_FILE,sizeof(char)); //!< Nucleotide or codon data?.
  io->run_id_string        = (char *)mCalloc(T_MAX_OPTION,sizeof(char)); //!< Name of run id.
  io->clade_list_file      = (char *)mCalloc(T_MAX_FILE,sizeof(char)); //!< Name of output clade list file.
  io->alphabet             = (char **)mCalloc(T_MAX_ALPHABET,sizeof(char *)); //!< Alphabet of symbols.
  For(i,T_MAX_ALPHABET) io->alphabet[i] = (char *)mCalloc(T_MAX_STATE,sizeof(char ));
  io->treelist             = (t_treelist *)mCalloc(1,sizeof(t_treelist)); //!< Stores each phylogeny explored during tree search.
  io->structTs_and_Tv      = (ts_and_tv *)mCalloc(64*64,sizeof(ts_and_tv));//!< Added by Marcelo. 64 possible codons ... will be initialized together with the model parameters in set model default.
  return io;
}

/*********************************************************/
//! Initialize options with default inputs.

/*!
 *
 */
void Set_Defaults_Input(option* io)
{
  io->fp_in_align                = NULL; //!< File pointer.
  io->fp_in_tree                 = NULL; //!< File pointer.
  io->fp_out_tree                = NULL; //!< File pointer.
  io->fp_out_trees               = NULL; //!< File pointer.
  io->fp_out_boot_tree           = NULL; //!< File pointer.
  io->fp_out_boot_stats          = NULL; //!< File pointer.
  io->fp_out_stats               = NULL; //!< File pointer.
  io->out_stats_format           = OUTTXT;
  io->long_tax_names             = NULL; //!< What is the meaning for two taxa name styles?.
  io->short_tax_names            = NULL; //!< What is the meaning for two taxa name styles?.
  io->lon                        = NULL; //!< What is this?.
  io->lat                        = NULL; //!< What is this?.
  io->z_scores                   = NULL; //!< What is this?.
  io->tree                       = NULL;
  io->mod                        = NULL;
  strcpy(io->nt_or_cd,"nucleotides");
  io->n_data_sets                = 1;
  io->interleaved                = 1;
  io->in_tree                    = 0;
  io->out_tree_file_open_mode    = 1;
  io->out_stats_file_open_mode   = 1;
  io->init_len                   = -1;
  io->n_otu                      = -1;
  io->n_data_set_asked           = -1;
  io->print_boot_trees           = 1;
  io->n_part                     = 1;
  io->ratio_test		 = 4;
  io->multigene                  = 0;
  io->config_multigene           = 0;
  io->curr_interface             = 0;
  io->r_seed                     = -1000;
  io->collapse_boot              = 0;
  io->random_boot_seq_order      = 1;
  io->print_trace                = 0;
  io->print_site_lnl             = 0;
  io->m4_model                   = NO;
  io->rm_ambigu                  = 0;
  io->append_run_ID              = 0;
  io->quiet                      = 0;
  io->datatype                   = NT;
  io->colalias                   = YES;
  io->data_file_format           = PHYLIP;
  io->tree_file_format           = PHYLIP;
  
  io->init_DistanceTreeCD        = KOSI07; //!< Added by Marcelo.
  io->kappaECM                   = kap1; //!< Added by Marcelo.
  io->fp_out_compare             = NULL; //!< Added by Marcelo.
  io->expm                       = EIGEN; //!< Added by Marcelo.
  io->lkExperiment               = NO; //!< Added by Marcelo.
  io->heuristicExpm              = NO; //!< Added by Marcelo. 
  io->threshold_exp              = NO; //!< Added by Marcelo. 
  io->testInitTree               = NO; //!< Added by Marcelo. 
  io->testcondition              = NO; //!< Added by Marcelo.
  io->convert_NT_to_AA           = NO; //!< Added by Marcelo.
  io->convert_AAfreq2CDfreq      = NO; //!< Added by Marcelo.
  io->fp_in_usrECM               = NULL; //!< Added by Marcelo.
  io->opt_heuristic_manuel       = NO; //!< Added by Marcelo.
  io->n_termsTaylor              = 3; //!< Added by Marcelo.
  io->init_run                   = 0; //!< Added by Marcelo.
  io->optParam                   = 0; //!< Added by Marcelo.
  io->roundMax_start             = 5; //!< Added by Marcelo.
  io->roundMax_end               = 10; //!< Added by Marcelo.
  io->eq_freq_handling           = NOFEQ; //!< Added by Marcelo.
  io->nobl                       = 0; //!< Added by Marcelo.
  
  io->genCode                    = STANDARD;
  io->omegaOpt                   = NOOMEGA;
  io->userWantsKappa             = NO;
  io->wclasses                   = 1;
  io->userOmega                  = (char *)mCalloc(T_MAX_NAME, sizeof(char)); 
  io->userOmega[0]               = 0;
  io->modeltypeOpt               = HKY85; // default value
  io->freqmodelOpt               = CF3X4;
  io->userFreqs                  = (char *)mCalloc(T_MAX_NAME, sizeof(char));
  io->userFreqs[0]               = 0;
  io->npcs                       = 1;
  io->writemode                  = 1;
  io->open_ps_file               = 0;
  io->confformat                 = -1;
  io->logtree                    = 0; // no logging per default
  io->treecounter                = 0; // reset counter
} 

/*********************************************************/
//! Initialize options with default inputs.

/*!
 *
 */
void Set_Defaults_Model(model *mod)
{
  int i;
  strcpy(mod->modelname,"HKY85");
  strcpy(mod->custom_mod_string,"000000");
  mod->whichmodel              = HKY85;
  mod->n_catg                  = 1;
  mod->kappa                   = 1.0;
  mod->alpha                   = 1.0;
  mod->lambda                  = 1.0;
  mod->bootstrap               = 0;
  mod->invar                   = 0;
  mod->pinvar                  = 0.0;
  mod->ns                      = 4;
  mod->n_diff_rr               = 0;
  mod->use_m4mod               = 0;
  mod->n_rr_branch             = 0;
  mod->rr_branch_alpha         = 0.1;
  mod->gamma_median            = 0;
  mod->state_len               = 1;
  mod->m4mod                   = NULL;
  mod->rr                      = NULL;
  mod->rr_val                  = NULL;
  mod->n_rr_per_cat            = NULL;
  mod->io                      = NULL;
  
  mod->genetic_code            = STANDARD; //!< Added by Marcelo.
  Genetic_code_index_stopCodons(mod->genetic_code); //!< Added by Marcelo.
  mod->freq_model              = CF3X4; //!< Added by Marcelo.
  mod->num_base_freq           = 12; //!< Added by Marcelo.
  //mod->omega_part                   = 1.0; //!< Added by Marcelo. //Ken 18/8
  mod->nkappa                  = 1; //!< Added by Marcelo.
  mod->pkappa[0]               = 1.0; //!< Added by Marcelo.
  mod->n_w_catg                = 1; //!< Added by Marcelo. 
  mod->beta                    = 1.0; //!< Added by Marcelo. 
  mod->omegaSiteVar            = DM0; //!< Added by Marcelo. 
  mod->codon_model_nature      = NOCHOICE; //!< Added by Marcelo. 
  mod->calculate_init_tree     = 0; //!< Added by Marcelo. 
  mod->npcs                    = 1; //!< Added by Marcelo. 
  For(i,100) mod->pcsC[i]      = 1.0; //!< Added by Marcelo. 
  mod->pcaModel                = 0; //!< Added by Marcelo.
  
  mod->initqrates              = NOINITMAT;
  mod->freq_model              = CF3X4;
}
/*********************************************************/

//! Initialize options with default inputs.

/*!
 *
 */
void Set_Defaults_Optimiz(optimiz *s_opt)
{
  s_opt->print                = 1;
  s_opt->last_opt             = 1;
  s_opt->opt_subst_param      = 1;
  s_opt->opt_alpha            = 1;
  s_opt->opt_bl               = 1;
  s_opt->opt_lambda           = 0;
  s_opt->opt_pinvar           = 0;
  s_opt->opt_cov_delta        = 0;
  s_opt->opt_cov_alpha        = 0;
  s_opt->opt_cov_free_rates   = 0;
  s_opt->opt_rr               = 0;
  s_opt->init_lk              = UNLIKELY;
  s_opt->n_it_max             = 1000;
  s_opt->opt_topo             = 1;
  s_opt->topo_search          = NNI_MOVE;
  s_opt->random_input_tree    = 0;
  s_opt->n_rand_starts        = 5;
  s_opt->brent_it_max         = 500;
  s_opt->steph_spr            = 1;
  s_opt->user_state_freq      = 0;
  s_opt->min_diff_lk_local    = 1.E-4;
  s_opt->min_diff_lk_global   = 1.E-3;
  s_opt->min_diff_lk_move     = 1.E-2;
  s_opt->p_moves_to_examine   = 0.15;
  s_opt->fast_nni             = 0;
  s_opt->greedy               = 0;
  s_opt->general_pars         = 0;
  s_opt->tree_size_mult       = 1;
  s_opt->opt_five_branch      = 1;
  s_opt->pars_thresh          = 5;
  s_opt->hybrid_thresh        = 0;
  s_opt->quickdirty           = 0;
  s_opt->spr_pars             = 1;
  s_opt->spr_lnL              = 0;
  s_opt->min_depth_path       = 0;
  s_opt->max_depth_path       = 20;
  s_opt->deepest_path         = 20;
  s_opt->max_delta_lnL_spr    = 50.;
  s_opt->wim_n_rgrft          = -1;
  s_opt->wim_n_globl          = -1;
  s_opt->wim_max_dist         = -1;
  s_opt->wim_n_optim          = -1;
  s_opt->wim_n_best           = -1;
  s_opt->wim_inside_opt       =  0;
  
  s_opt->opt_omega               = YES; //!< Added By Marcelo.
  s_opt->opt_state_freq          = NO; //!< Added By Marcelo.
  s_opt->opt_beta                = YES; //!< Added by Marcelo. 
  s_opt->opt_prob_omega          = YES; //!< Added By Marcelo.
  s_opt->opt_alphaCD             = YES; //!< Added By Marcelo.
  s_opt->min_diff_lk_codonModels = 1e-3; //!< Added By Marcelo.
  s_opt->opt_method              = optPAML; //!< Added By Marcelo.
  s_opt->nBrentCycles            = 2;  //!< Added By Marcelo.
  s_opt->opt_state_freq_AAML     = NO;
  
  s_opt->opt_kappa               = YES; // HKY85 is default
}

/*********************************************************/

void Get_Bip(t_node *a, t_node *d, t_tree *tree)
{
  int i,j;
  t_node *tmp;
  int swapped;
  
  if(d->tax)
  {
    if(d->common)
    {
      d->bip_node[0] = (t_node **)mCalloc(1,sizeof(t_node *));
      d->bip_node[0][0] = d;
      d->bip_size[0]    = 1;
      
      For(i,3)
	    {
	      if(a->v[i] == d)
        {
          a->bip_size[i] = 0;
          For(j,tree->n_otu)
          {
            if(strcmp(tree->noeud[j]->name,d->name))
            {
              a->bip_node[i] = (t_node **)realloc(a->bip_node[i],(a->bip_size[i]+1)*sizeof(t_node *));
              a->bip_node[i][a->bip_size[i]] = tree->noeud[j];
              a->bip_size[i]++;
            }
          }
          
          /* Sort bipartition */
          do
          {
            swapped = NO;
            For(j,a->bip_size[i]-1)
            {
              if(a->bip_node[i][j]->num > a->bip_node[i][j+1]->num)
              {
                swapped = YES;
                tmp                 = a->bip_node[i][j];
                a->bip_node[i][j]   = a->bip_node[i][j+1];
                a->bip_node[i][j+1] = tmp;
              }
            }
          }while(swapped == YES);
          
          break;
          
        }
	    }
    }
    return;
  }
  else
  {
    int k;
    int d_a;
    
    d_a = -1;
    
    For(i,3)
    {
      if(d->v[i] != a) Get_Bip(d,d->v[i],tree);
      else if(d->v[i] == a) d_a = i;
    }
    
    d->bip_size[d_a] = 0;
    For(i,3)
    if(d->v[i] != a)
	  {
	    For(j,3)
      {
        if(d->v[i]->v[j] == d)
        {
          For(k,d->v[i]->bip_size[j])
		      {
            d->bip_node[d_a] = (t_node **)realloc(d->bip_node[d_a],(d->bip_size[d_a]+1)*sizeof(t_node *));
            d->bip_node[d_a][d->bip_size[d_a]] = d->v[i]->bip_node[j][k];
            d->bip_size[d_a]++;
		      }
          break;
        }
      }
	  }
    
    do
    {
      swapped = NO;
      For(j,d->bip_size[d_a]-1)
	    {
	      if(d->bip_node[d_a][j]->num > d->bip_node[d_a][j+1]->num)
        {
          swapped = YES;
          tmp                   = d->bip_node[d_a][j];
          d->bip_node[d_a][j]   = d->bip_node[d_a][j+1];
          d->bip_node[d_a][j+1] = tmp;
		    }
	    }
    }while(swapped == YES);
    
    
    For(i,3)
    if(a->v[i] == d)
	  {
	    a->bip_size[i] = 0;
	    For(j,tree->n_otu)
      {
        For(k,d->bip_size[d_a])
        {
          if(d->bip_node[d_a][k] == tree->noeud[j])
            break;
        }
        
        if((k == d->bip_size[d_a]) && (tree->noeud[j]->common))
        {
          a->bip_node[i] = (t_node **)realloc(a->bip_node[i],(a->bip_size[i]+1)*sizeof(t_node *));
          a->bip_node[i][a->bip_size[i]] = tree->noeud[j];
          a->bip_size[i]++;
        }
      }
      
	    do
      {
        swapped = NO;
        For(j,a->bip_size[i]-1)
        {
          if(a->bip_node[i][j]->num > a->bip_node[i][j+1]->num)
		      {
            swapped = YES;
            tmp                 = a->bip_node[i][j];
            a->bip_node[i][j]   = a->bip_node[i][j+1];
            a->bip_node[i][j+1] = tmp;
		      }
        }
      }while(swapped == YES);
	    
	    if(a->bip_size[i] != tree->n_otu - d->bip_size[d_a])
      {
        PhyML_Printf("%d %d \n",a->bip_size[i],tree->n_otu - d->bip_size[d_a]);
        Warn_And_Exit("\n. Problem in counting bipartitions \n");
      }
	    break;
	  }
  }
}

/*********************************************************/

void Alloc_Bip(t_tree *tree)
{
  int i;
  // int j,k;
  
  if(tree->has_bip) return;
  
  tree->has_bip = YES;
  
  For(i,2*tree->n_otu-2)
  {
    tree->noeud[i]->bip_size = (int *)mCalloc(3,sizeof(int));
    tree->noeud[i]->bip_node = (t_node ***)mCalloc(3,sizeof(t_node **));
    
    /*       For(j,3) */
    /* 	{ */
    /* 	  tree->noeud[i]->bip_node[j] = (t_node **)mCalloc(tree->n_otu,sizeof(t_node *)); */
    /* 	} */
  }
}

/*********************************************************/

void Compare_Bip(t_tree *tree1, t_tree *tree2)
{
  int i,j,k;
  t_edge *b1,*b2;
  /*   char **bip1,**bip2; */
  /*   int *bip1,*bip2; */
  t_node **bip1, **bip2;
  int bip_size1, bip_size2, bip_size;
  
  
  
  For(i,2*tree1->n_otu-3)
  {
    b1 = tree1->t_edges[i];     
    bip_size1 = MIN(b1->left->bip_size[b1->l_r],b1->rght->bip_size[b1->r_l]);
    
    if(bip_size1 > 1)
    {
      For(j,2*tree2->n_otu-3)
	    {
	      b2 = tree2->t_edges[j];	      
	      bip_size2 = MIN(b2->left->bip_size[b2->l_r],b2->rght->bip_size[b2->r_l]);
        
	      if(bip_size2 > 1)
        {
          if(bip_size1 == bip_size2)
          {
            bip_size = bip_size1;
            
            if(b1->left->bip_size[b1->l_r] == b1->rght->bip_size[b1->r_l])
            {
              /* 			  if(b1->left->bip_name[b1->l_r][0][0] < b1->rght->bip_name[b1->r_l][0][0]) */
              if(b1->left->bip_node[b1->l_r][0]->num < b1->rght->bip_node[b1->r_l][0]->num)
              {
                /* 			      bip1 = b1->left->bip_name[b1->l_r]; */
                bip1 = b1->left->bip_node[b1->l_r];
              }
              else
              {
                /* 			      bip1 = b1->rght->bip_name[b1->r_l]; */
                bip1 = b1->rght->bip_node[b1->r_l];
              }
            }
            else if(b1->left->bip_size[b1->l_r] < b1->rght->bip_size[b1->r_l])
            {
              /* 			  bip1 = b1->left->bip_name[b1->l_r]; */
              bip1 = b1->left->bip_node[b1->l_r];
            }
            else
            {
              /* 			  bip1 = b1->rght->bip_name[b1->r_l]; */
              bip1 = b1->rght->bip_node[b1->r_l];
            }
            
            
            if(b2->left->bip_size[b2->l_r] == b2->rght->bip_size[b2->r_l])
            {
              /* 			  if(b2->left->bip_name[b2->l_r][0][0] < b2->rght->bip_name[b2->r_l][0][0]) */
              if(b2->left->bip_node[b2->l_r][0]->num < b2->rght->bip_node[b2->r_l][0]->num)
              {
                /* 			      bip2 = b2->left->bip_name[b2->l_r]; */
                bip2 = b2->left->bip_node[b2->l_r];
              }
              else
              {
                /* 			      bip2 = b2->rght->bip_name[b2->r_l]; */
                bip2 = b2->rght->bip_node[b2->r_l];
              }
            }
            else if(b2->left->bip_size[b2->l_r] < b2->rght->bip_size[b2->r_l])
            {
              /* 			  bip2 = b2->left->bip_name[b2->l_r]; */
              bip2 = b2->left->bip_node[b2->l_r];
            }
            else
            {
              /* 			  bip2 = b2->rght->bip_name[b2->r_l]; */
              bip2 = b2->rght->bip_node[b2->r_l];
            }
            
            if(bip_size == 1) Warn_And_Exit("\n. Problem in Compare_Bip\n");
            
            For(k,bip_size)
            {
              /* 			  if(strcmp(bip1[k],bip2[k])) break; */
              if(bip1[k]->num != bip2[k]->num) break;
            }
            
            if(k == bip_size)
            {
              b1->bip_score++;
              b2->bip_score++;
              break;
            }
          }
        }
	    }
    }
  }
}

/*********************************************************/

void Match_Tip_Numbers(t_tree *tree1, t_tree *tree2)
{
  int i,j;
  
  For(i,tree1->n_otu)
  {
    For(j,tree2->n_otu)
    {
      if(!strcmp(tree1->noeud[i]->name,tree2->noeud[j]->name))
	    {
	      tree2->noeud[j]->num = tree1->noeud[i]->num;
	      break;
	    }
    }
  }
  
}

/*********************************************************/

void Test_Multiple_Data_Set_Format(option *io)
{
  char *line;
  
  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  
  io->n_trees = 0;
  
  while(fgets(line,T_MAX_LINE,io->fp_in_tree)) if(strstr(line,";")) io->n_trees++;
  
  free(line);
  
  if((io->mod->bootstrap > 1) && (io->n_trees > 1))
    Warn_And_Exit("\n. Bootstrap option is not allowed with multiple input trees !\n");
  
  rewind(io->fp_in_tree);
  
  return;
}

/*********************************************************/

int Are_Compatible(char *statea, char *stateb, int stepsize, int datatype, char * alta, char * altb)
{
  int i,j;
  char a,b;
  
  if(datatype == NT)
  {
    For(i,stepsize)
    {
      a = statea[i];
      For(j,stepsize)
	    {
	      b = stateb[j];
        
	      switch(a)
        {
          case 'A':
          {
            switch(b)
            {
              case 'A' :
              case 'M' :
              case 'R' :
              case 'W' :
              case 'D' :
              case 'H' :
              case 'V' :
              case 'X' : {b=b; break;}
              default : return 0;
            }
            break;
          }
          case 'G':
          {
            switch(b)
            {
              case 'G' :
              case 'R' :
              case 'S' :
              case 'K' :
              case 'B' :
              case 'D' :
              case 'V' :
              case 'X' : {b=b; break;}
              default : return 0;
            }
            break;
          }
          case 'C':
          {
            switch(b)
            {
              case 'C' :
              case 'M' :
              case 'S' :
              case 'Y' :
              case 'B' :
              case 'H' :
              case 'V' :
              case 'X' : {b=b; break;}
              default : return 0;
            }
            break;
          }
          case 'T':
          {
            switch(b)
            {
              case 'T' :
              case 'W' :
              case 'Y' :
              case 'K' :
              case 'B' :
              case 'D' :
              case 'H' :
              case 'X' :
              {b=b; break;}
              default : return 0;
            }
            break;
          }
          case 'M' :
          {
            switch(b)
            {
              case 'M' :
              case 'A' :
              case 'C' :
              case 'R' :
              case 'W' :
              case 'S' :
              case 'Y' :
              case 'B' :
              case 'D' :
              case 'H' :
              case 'V' :
              case 'X' :
              {b=b; break;}
              default : return 0;
            }
            break;
          }
          case 'R' :
          {
            switch(b)
            {
              case 'R' :
              case 'A' :
              case 'G' :
              case 'M' :
              case 'W' :
              case 'S' :
              case 'K' :
              case 'B' :
              case 'D' :
              case 'H' :
              case 'V' :
              case 'X' : {b=b; break;}
              default : return 0;
            }
            break;
          }
            
          case 'W' :
          {
            switch(b)
            {
              case 'W' :
              case 'A' :
              case 'T' :
              case 'M' :
              case 'R' :
              case 'Y' :
              case 'K' :
              case 'B' :
              case 'D' :
              case 'H' :
              case 'V' :
              case 'X' : {b=b; break;}
              default : return 0;
            }
            break;
          }
            
          case 'S' :
          {
            switch(b)
            {
              case 'S' :
              case 'C' :
              case 'G' :
              case 'M' :
              case 'R' :
              case 'Y' :
              case 'K' :
              case 'B' :
              case 'D' :
              case 'H' :
              case 'V' :
              case 'X' : {b=b; break;}
              default : return 0;
            }
            break;
          }
            
          case 'Y' :
          {
            switch(b)
            {
              case 'Y' :
              case 'C' :
              case 'T' :
              case 'M' :
              case 'W' :
              case 'S' :
              case 'K' :
              case 'B' :
              case 'D' :
              case 'H' :
              case 'V' :
              case 'X' : {b=b; break;}
              default : return 0;
            }
            break;
          }
            
          case 'K' :
          {
            switch(b)
            {
              case 'K' :
              case 'G' :
              case 'T' :
              case 'R' :
              case 'W' :
              case 'S' :
              case 'Y' :
              case 'B' :
              case 'D' :
              case 'H' :
              case 'V' :
              case 'X' : {b=b; break;}
              default : return 0;
            }
            break;
          }
          case 'B' :
          {
            switch(b)
            {
              case 'B' :
              case 'C' :
              case 'G' :
              case 'T' :
              case 'M' :
              case 'R' :
              case 'W' :
              case 'S' :
              case 'Y' :
              case 'K' :
              case 'D' :
              case 'H' :
              case 'V' :
              case 'X' : {b=b; break;}
              default : return 0;
            }
            break;
          }
          case 'D' :
          {
            switch(b)
            {
              case 'D' :
              case 'A' :
              case 'G' :
              case 'T' :
              case 'M' :
              case 'R' :
              case 'W' :
              case 'S' :
              case 'Y' :
              case 'K' :
              case 'B' :
              case 'H' :
              case 'V' :
              case 'X' : {b=b; break;}
              default : return 0;
            }
            break;
          }
          case 'H' :
          {
            switch(b)
            {
              case 'H' :
              case 'A' :
              case 'C' :
              case 'T' :
              case 'M' :
              case 'R' :
              case 'W' :
              case 'S' :
              case 'Y' :
              case 'K' :
              case 'B' :
              case 'D' :
              case 'V' :
              case 'X' : {b=b; break;}
              default : return 0;
            }
            break;
          }
          case 'V' :
          {
            switch(b)
            {
              case 'V' :
              case 'A' :
              case 'C' :
              case 'G' :
              case 'M' :
              case 'R' :
              case 'W' :
              case 'S' :
              case 'Y' :
              case 'K' :
              case 'B' :
              case 'D' :
              case 'H' :
              case 'X' : {b=b; break;}
              default : return 0;
            }
            break;
          }
          case 'X' :
          {
            switch(b)
            {
              case 'X' :
              case 'A' :
              case 'C' :
              case 'G' :
              case 'T' :
              case 'M' :
              case 'R' :
              case 'W' :
              case 'S' :
              case 'Y' :
              case 'K' :
              case 'B' :
              case 'D' :
              case 'H' :
              case 'V' : {b=b; break;}
              default : return 0;
            }
            break;
          }
          default :
          {
            PhyML_Printf("\n. Err. in Are_Compatible\n");
            PhyML_Printf("\n. Please check that characters `%c` and `%c`\n",a,b);
            PhyML_Printf("  correspond to existing nucleotides.\n");
            Warn_And_Exit("\n");
            return 0;
          }
        }
	    }
    }
  }
  else if(datatype == AA)
  {
    a = statea[0]; b = stateb[0];
    switch(a)
    {
      case 'A' :
      {
        switch(b)
	      {
          case 'A' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'R' :
      {
        switch(b)
	      {
          case 'R' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'N' :
      {
        switch(b)
	      {
          case 'N' :
          case 'B' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'B' :
      {
        switch(b)
	      {
          case 'N' :
          case 'B' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'D' :
      {
        switch(b)
	      {
          case 'D' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'C' :
      {
        switch(b)
	      {
          case 'C' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'Q' :
      {
        switch(b)
	      {
          case 'Q' :
          case 'Z' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'Z' :
      {
        switch(b)
	      {
          case 'Q' :
          case 'Z' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'E' :
      {
        switch(b)
	      {
          case 'E' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'G' :
      {
        switch(b)
	      {
          case 'G' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'H' :
      {
        switch(b)
	      {
          case 'H' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'I' :
      {
        switch(b)
	      {
          case 'I' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'L' :
      {
        switch(b)
	      {
          case 'L' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'K' :
      {
        switch(b)
	      {
          case 'K' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'M' :
      {
        switch(b)
	      {
          case 'M' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'F' :
      {
        switch(b)
	      {
          case 'F' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'P' :
      {
        switch(b)
	      {
          case 'P' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'S' :
      {
        switch(b)
	      {
          case 'S' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'T' :
      {
        switch(b)
	      {
          case 'T' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'W' :
      {
        switch(b)
	      {
          case 'W' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'Y' :
      {
        switch(b)
	      {
          case 'Y' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'V' :
      {
        switch(b)
	      {
          case 'V' :
          case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      case 'X' :
      {
        switch(b)
	      {
          case 'A':case 'R':case 'N' :case 'B' :case 'D' :
          case 'C':case 'Q':case 'Z' :case 'E' :case 'G' :
          case 'H':case 'I':case 'L' :case 'K' :case 'M' :
          case 'F':case 'P':case 'S' :case 'T' :case 'W' :
          case 'Y':case 'V': case 'X' : {b=b; break;}
          default : return 0;
	      }
        break;
      }
      default :
      {
        PhyML_Printf("\n. Err. in Are_Compatible\n");
        PhyML_Printf("\n. Please check that characters `%c` and `%c`\n",a,b);
        PhyML_Printf("  correspond to existing amino-acids.\n");
        Warn_And_Exit("\n");
        return 0;
      }
    }
  }
  else if(datatype == GENERIC)    
  {
    if(Is_Ambigu(statea,GENERIC,stepsize) || Is_Ambigu(stateb,GENERIC,stepsize)) return 1;
    else
    {
      int a1,b1;
      char format[6];      
      
      sprintf(format,"%%%dd",stepsize);      
      
      if(!sscanf(statea,format,&a1))
	    {	    
	      PhyML_Printf("\n. statea = %s",statea);
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
      if(!sscanf(stateb,format,&b1))
	    {	    
	      PhyML_Printf("\n. statea = %s",stateb);
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
      
      /* 	  PhyML_Printf("\n. %s %d a=%d b=%d ",__FILE__,__LINE__,a,b);  */
      
      if(a1 == b1) return 1;
    }
    return 0;
  }
  else if(datatype == CODON)//!< Added by Marcelo.
  {
    int k,m;
    // int n;
    a = *statea;
    b = *stateb;
    
    if (((a>=(char)0)&&(a<(char)64))&&((b>=(char)0)&&(b<(char)64)))
    {
      if(a==b) 
        return 1; 
      else 
        return 0;
    }
    else if (((a>=(char)0)&&(a<(char)64))&&((b==(char)88)))
    {
      k=-1;
      while(altb[++k]<64);
      For(i,k) if(a==altb[i]) return 1; 
      return 0;
    }
    else if (((b>=(char)0)&&(b<(char)64))&&((a==(char)88)))
    {
      k=-1;
      while(alta[++k]<(char)64);
      For(i,k) if(b==alta[i]) return 1; 
      return 0;
    }
    else if (((a==(char)88))&&((b==(char)88)))
    {
      k=-1;
      while(alta[++k]<(char)64);
      m=-1;
      while(altb[++m]<(char)64);
      
      For(i,k) For(j,m) if(alta[i]==altb[j]) return 1; 
      return 0;
    }
    else 
    {
      PhyML_Printf("\n. Err. in Are_Compatible\n");
      PhyML_Printf("\n. Please check that characters `%c` and `%c`\n",a,b);
      PhyML_Printf("  correspond to existing codons (int)[0,63].\n");
      Warn_And_Exit("\n");
    }
  }
  return 1;
}

/*********************************************************/

void Hide_Ambiguities(calign *data)
{
  int i;
  For(i,data->crunch_len) if(data->ambigu[i]) data->wght[i] = 0;
}

/*********************************************************/

void Copy_Tree(t_tree *ori, t_tree *cpy)
{
  int i,j;
  
  
  For(i,2*ori->n_otu-2)
  {
    For(j,3)
    {
      if(ori->noeud[i]->v[j])
	    {
	      cpy->noeud[i]->v[j] = cpy->noeud[ori->noeud[i]->v[j]->num];
	      cpy->noeud[i]->l[j] = ori->noeud[i]->l[j];
	      cpy->noeud[i]->b[j] = cpy->t_edges[ori->noeud[i]->b[j]->num];
	    }
      else
	    {
	      cpy->noeud[i]->v[j] = NULL;
	      cpy->noeud[i]->b[j] = NULL;
	    }
    }
  }
  
  For(i,2*ori->n_otu-3) 
  {
    cpy->t_edges[i]->l    = ori->t_edges[i]->l;
    cpy->t_edges[i]->left = cpy->noeud[ori->t_edges[i]->left->num];
    cpy->t_edges[i]->rght = cpy->noeud[ori->t_edges[i]->rght->num];
    cpy->t_edges[i]->l_v1 = ori->t_edges[i]->l_v1;
    cpy->t_edges[i]->l_v2 = ori->t_edges[i]->l_v2;
    cpy->t_edges[i]->r_v1 = ori->t_edges[i]->r_v1;
    cpy->t_edges[i]->r_v2 = ori->t_edges[i]->r_v2;
    cpy->t_edges[i]->l_r  = ori->t_edges[i]->l_r;
    cpy->t_edges[i]->r_l  = ori->t_edges[i]->r_l;
  }
  
  For(i,ori->n_otu)
  {
    cpy->noeud[i]->tax = 1;
    strcpy(cpy->noeud[i]->name,ori->noeud[i]->name);
  }
  
  cpy->num_curr_branch_available = 0;
  /*   Connect_Edges_To_Nodes_Recur(cpy->noeud[0],cpy->noeud[0]->v[0],cpy); */
  /*   Update_Dirs(cpy); */
  
  cpy->c_lnL=ori->c_lnL;//!< Added by Marcelo.
  
  
}

/*********************************************************/
//a: n_link, d: n_oppo_to_link
void Prune_Subtree(t_node *a, t_node *d, t_edge **target, t_edge **residual, t_tree *tree)
{
  t_node *v1, *v2;
  t_edge *b1, *b2;
  int dir_v1, dir_v2;
  int i;
  /*   phydbl ***buff_p_lk; */
  phydbl *buff_p_lk;
  int *buff_scale;
  int *buff_p_pars, *buff_pars;
  unsigned int *buff_ui;
  short int *buff_p_lk_tip;
  
  if(a->tax)
  {
    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    Warn_And_Exit("");
  }
  
  dir_v1 = dir_v2 = -1;
  For(i,3){ //v1 and v2 are the pendant edges connected via a->d
    if(a->v[i] != d){
      if(dir_v1 < 0) dir_v1 = i;
      else           dir_v2 = i;
    }
  }
  
  if(a->v[dir_v1]->num < a->v[dir_v2]->num) //assign v1 or v2 based on which node number is higher?
  {
    v1 = a->v[dir_v1];
    v2 = a->v[dir_v2];
    b1 = a->b[dir_v1];
    b2 = a->b[dir_v2];
  }else{
    v1 = a->v[dir_v2];
    v2 = a->v[dir_v1];
    b1 = a->b[dir_v2];
    b2 = a->b[dir_v1];
  }
  
  if(v1->tax && v2->tax) PhyML_Printf("\n. Pruning is meaningless here.\n");
  
  a->v[dir_v1] = NULL; //clear pointers in a node (link node)
  a->v[dir_v2] = NULL;
  a->b[dir_v1] = NULL;
  a->b[dir_v2] = NULL;
  
  if(v1 == b1->left)
  {
    b1->rght = v2;
    if(v2 == b2->left)
    {
      buff_p_lk            = b1->p_lk_rght;
      b1->p_lk_rght        = b2->p_lk_left;
      b2->p_lk_left        = buff_p_lk;
      
      buff_p_lk_tip        = b1->p_lk_tip_r;
      b1->p_lk_tip_r       = b2->p_lk_tip_l;
      b2->p_lk_tip_l       = buff_p_lk_tip;
      
      buff_scale           = b1->sum_scale_rght;
      b1->sum_scale_rght   = b2->sum_scale_left;
      b2->sum_scale_left   = buff_scale;
      
      buff_scale             = b1->sum_scale_rght_cat;
      b1->sum_scale_rght_cat = b2->sum_scale_left_cat;
      b2->sum_scale_left_cat = buff_scale;
      
      buff_pars            = b1->pars_r;
      b1->pars_r           = b2->pars_l;
      b2->pars_l           = buff_pars;
      
      buff_ui              = b1->ui_r;
      b1->ui_r             = b2->ui_l;
      b2->ui_l             = buff_ui;
      
      buff_p_pars          = b1->p_pars_r;
      b1->p_pars_r         = b2->p_pars_l;
      b2->p_pars_l         = buff_p_pars;
    }
    else
    {
      buff_p_lk            = b1->p_lk_rght; /* b1->p_lk_rght = NULL if b1->rght->tax */
      b1->p_lk_rght        = b2->p_lk_rght; /* b2->p_lk_rght = NULL if b2->rght->tax */ 
      b2->p_lk_rght        = buff_p_lk;
      
      buff_p_lk_tip        = b1->p_lk_tip_r;
      b1->p_lk_tip_r       = b2->p_lk_tip_r;
      b2->p_lk_tip_r       = buff_p_lk_tip;
      
      buff_scale           = b1->sum_scale_rght;
      b1->sum_scale_rght = b2->sum_scale_rght;
      b2->sum_scale_rght = buff_scale;
      
      buff_pars            = b1->pars_r;
      b1->pars_r           = b2->pars_r;
      b2->pars_r           = buff_pars;
      
      buff_ui              = b1->ui_r;
      b1->ui_r             = b2->ui_r;
      b2->ui_r             = buff_ui;
      
      buff_p_pars          = b1->p_pars_r;
      b1->p_pars_r         = b2->p_pars_r;
      b2->p_pars_r         = buff_p_pars;
    }
  }
  else
  {
    b1->left = v2;
    
    if(v2 == b2->left)
    {
      buff_p_lk            = b1->p_lk_left;
      b1->p_lk_left        = b2->p_lk_left;
      b2->p_lk_left        = buff_p_lk;
      
      buff_p_lk_tip        = b1->p_lk_tip_l;
      b1->p_lk_tip_l       = b2->p_lk_tip_l;
      b2->p_lk_tip_l       = buff_p_lk_tip;
      
      buff_scale           = b1->sum_scale_left;
      b1->sum_scale_left = b2->sum_scale_left;
      b2->sum_scale_left = buff_scale;
      
      buff_scale             = b1->sum_scale_left_cat;
      b1->sum_scale_left_cat = b2->sum_scale_left_cat;
      b2->sum_scale_left_cat = buff_scale;
      
      buff_pars            = b1->pars_l;
      b1->pars_l           = b2->pars_l;
      b2->pars_l           = buff_pars;
      
      buff_ui              = b1->ui_l;
      b1->ui_l             = b2->ui_l;
      b2->ui_l             = buff_ui;
      
      buff_p_pars          = b1->p_pars_l;
      b1->p_pars_l         = b2->p_pars_l;
      b2->p_pars_l         = buff_p_pars;
    }
    else
    {
      buff_p_lk            = b1->p_lk_left;
      b1->p_lk_left        = b2->p_lk_rght; /* b2->p_lk_rght = NULL if b2->rght->tax */
      b2->p_lk_rght        = buff_p_lk;
      
      buff_p_lk_tip        = b1->p_lk_tip_l;
      b1->p_lk_tip_l       = b2->p_lk_tip_r;
      b2->p_lk_tip_r       = buff_p_lk_tip;
      
      buff_scale           = b1->sum_scale_left;
      b1->sum_scale_left = b2->sum_scale_rght;
      b2->sum_scale_rght = buff_scale;
      
      buff_scale             = b1->sum_scale_left_cat;
      b1->sum_scale_left_cat = b2->sum_scale_rght_cat;
      b2->sum_scale_rght_cat = buff_scale;
      
      buff_pars            = b1->pars_l;
      b1->pars_l           = b2->pars_r;
      b2->pars_r           = buff_pars;
      
      buff_ui              = b1->ui_l;
      b1->ui_l             = b2->ui_r;
      b2->ui_r             = buff_ui;
      
      buff_p_pars          = b1->p_pars_l;
      b1->p_pars_l         = b2->p_pars_r;
      b2->p_pars_r         = buff_p_pars;
    }
  }
  
  For(i,3)
  if(v2->v[i] == a)
  {
    v2->v[i] = v1;
    v2->b[i] = b1;
    break;
  }
  
#ifdef DEBUG
  if(i == 3)
  {
    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    Warn_And_Exit("");
  }
#endif
  
  For(i,3)
  if(v1->v[i] == a)
  {
    v1->v[i] = v2;
    break;
  }
  
#ifdef DEBUG
  if(i == 3)
  {
    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    Warn_And_Exit("");
  }
#endif
  
  b1->l += b2->l;
  
  (v1 == b1->left)?
  (Make_Edge_Dirs(b1,v1,v2)):
  (Make_Edge_Dirs(b1,v2,v1));
  
  if(target)   (*target)   = b1; //one of the edges is left behind, the other connects the two
  if(residual) (*residual) = b2;
  
  
#ifdef DEBUG
  if(b1->left->tax)
  {
    PhyML_Printf("\n. b1->left->num = %d",b1->left->num);
    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    Warn_And_Exit("");
  }
#endif
  
  if(tree->mod->whichrealmodel==HLP17){ //added by Ken 17/11
 	  Update_Ancestors_Edge(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree->noeud[tree->mod->startnode]->b[0],tree); //added by Ken 7/11
	  if(tree->mod->slowSPR==1){
	  	  tree->both_sides=1;
	  	  Lk(tree);
	  	  Get_UPP(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree);
	  }
   }
}

/*********************************************************/

void Graft_Subtree(t_edge *target, t_node *link, t_edge *residual, t_tree *tree)
{
  t_node *v1, *v2;
  int i, dir_v1, dir_v2;
  phydbl *buff_p_lk;
  int *buff_scale;
  int *buff_p_pars, *buff_pars; 
  short int *buff_p_lk_tip;
  unsigned int *buff_ui;
  t_edge *b_up;
  
  dir_v1 = dir_v2 = -1;
  b_up = NULL;
  For(i,3)
  {
    if(!link->v[i])
    {
      if(dir_v1 < 0) dir_v1 = i;
      else           dir_v2 = i;
    }
    else b_up = link->b[i];
  }
  
  if(target->left->num < target->rght->num)
  {
    v1                           = target->left;
    v2                           = target->rght;
    
    buff_p_lk                    = residual->p_lk_rght;
    residual->p_lk_rght          = target->p_lk_rght;
    target->p_lk_rght            = buff_p_lk;
    
    buff_p_lk_tip                = residual->p_lk_tip_r;
    residual->p_lk_tip_r         = target->p_lk_tip_r;
    target->p_lk_tip_r           = buff_p_lk_tip;
    
    buff_scale                   = residual->sum_scale_rght;
    residual->sum_scale_rght   = target->sum_scale_rght;
    target->sum_scale_rght     = buff_scale;
    
    buff_pars                    = residual->pars_r;
    residual->pars_r             = target->pars_r;
    target->pars_r               = buff_pars;
    
    buff_ui                      = residual->ui_r;
    residual->ui_r               = target->ui_r;
    target->ui_r                 = buff_ui;
    
    buff_p_pars                  = residual->p_pars_r;
    residual->p_pars_r           = target->p_pars_r;
    target->p_pars_r             = buff_p_pars;
  }
  else
  {
    v1                           = target->rght;
    v2                           = target->left;
    
    buff_p_lk                    = residual->p_lk_rght;
    residual->p_lk_rght          = target->p_lk_left;
    target->p_lk_left            = buff_p_lk;
    
    buff_p_lk_tip                = residual->p_lk_tip_r;
    residual->p_lk_tip_r         = target->p_lk_tip_l;
    target->p_lk_tip_l           = buff_p_lk_tip;
    
    buff_scale                   = residual->sum_scale_rght;
    residual->sum_scale_rght     = target->sum_scale_left;
    target->sum_scale_left       = buff_scale;
    
    buff_scale                   = residual->sum_scale_rght_cat;
    residual->sum_scale_rght_cat = target->sum_scale_left_cat;
    target->sum_scale_left_cat   = buff_scale;
    
    buff_pars                    = residual->pars_r;
    residual->pars_r             = target->pars_l;
    target->pars_l               = buff_pars;
    
    buff_ui                      = residual->ui_r;
    residual->ui_r               = target->ui_l;
    target->ui_l                 = buff_ui;
    
    buff_p_pars                  = residual->p_pars_r;
    residual->p_pars_r           = target->p_pars_l;
    target->p_pars_l             = buff_p_pars;
  }
  
  For(i,3)
  if(v2->b[i] == target)
  {
    v2->v[i] = link;
    v2->b[i] = residual;
    break;
  }
  
  link->v[dir_v2] = v2;
  link->b[dir_v2] = residual;
  
  residual->left  = link;
  residual->rght  = v2;
  
  (v1 == target->left)?(target->rght = link):(target->left = link);
  
  link->v[dir_v1] = v1;
  link->b[dir_v1] = target;
  
  For(i,3)
  if(v1->v[i] == v2)
  {
    v1->v[i] = link;
    break;
  }
  
  target->l /= 2.;
  residual->l = target->l;
  
  Make_Edge_Dirs(target,target->left,target->rght);
  Make_Edge_Dirs(residual,residual->left,residual->rght);
  Make_Edge_Dirs(b_up,b_up->left,b_up->rght);

  if(tree->mod->whichrealmodel==HLP17){ //added by Ken 17/11
	  Update_Ancestors_Edge(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree->noeud[tree->mod->startnode]->b[0],tree); //added by Ken 7/11
	  if(tree->mod->slowSPR==1){
	  	  tree->both_sides=1;
	  	  Lk(tree);
	  	  Get_UPP(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree);
	  }else{
	  	  Update_PMat_At_Given_Edge(target,tree);
	  	  Update_PMat_At_Given_Edge(residual,tree);
	  }
  }else{
	//  tree->both_sides=1;
	//  Lk(tree);
  }
}

/*********************************************************/

void Find_Mutual_Direction(t_node *n1, t_node *n2, short int *dir_n1_to_n2, short int *dir_n2_to_n1)
{
  int scores[3][3];
  // int n_zero_line, n_zero_col;
  int i,j,k,l;
  // int n_otu;
  
  if(n1 == n2) return;
  
  
  For(i,3)
  {
    For(j,3)
    {
      scores[i][j] = 0;
      
      For(k,n1->bip_size[i])
	    {
	      For(l,n2->bip_size[j])
        {
          if(n1->bip_node[i][k] == n2->bip_node[j][l])
          {
            scores[i][j]++;
            break;
          }
        }
	    }
    }
  }
  
  For(i,3)
  {
    For(j,3)
    {
      if(!scores[i][j]) 
	    {
	      *dir_n1_to_n2 = i; 
	      *dir_n2_to_n1 = j; 
	      return;
	    } 
    }
  }
  
  PhyML_Printf("\n. n1=%d n2=%d",n1->num,n2->num);
  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
  Warn_And_Exit("");
  
  
  
  /*   For(i,3) */
  /*     { */
  /*       n_zero_line = 0; */
  /*       For(j,3) */
  /* 	{ */
  /* 	  if(!scores[i][j]) n_zero_line++; */
  /* 	} */
  /*       if(n_zero_line != 2) {*dir_n1_to_n2 = i; break;} */
  /*     } */
  
  
  /*   For(i,3) */
  /*     { */
  /*       n_zero_col = 0; */
  /*       For(j,3) */
  /* 	{ */
  /* 	  if(!scores[j][i]) n_zero_col++; */
  /* 	} */
  /*       if(n_zero_col != 2) {*dir_n2_to_n1 = i; break;} */
  /*     } */
  
}

/*********************************************************/

void Update_Dir_To_Tips(t_node *a, t_node *d, t_tree *tree)
{
  int i,j,k;
  short int *inout;
  int d_a;
  int dim;
  
  dim = 2*tree->n_otu-2;
  
  inout = (short int *)mCalloc(tree->n_otu,sizeof(short int));
  
  For(i,3)
  {
    if(a->v[i] == d)
    {
      For(j,tree->n_otu) inout[j] = 1;
      For(k,a->bip_size[i]) inout[a->bip_node[i][k]->num] = 0;
      For(j,tree->n_otu) if(inout[tree->noeud[j]->num]) tree->t_dir[a->num*dim+tree->noeud[j]->num] = i;
      break;
    }
  }
  
  
  if(!d->tax)
  {
    
    d_a = -1;
    
    For(i,3)
    {
      if(d->v[i] != a && d->b[i] != tree->e_root) Update_Dir_To_Tips(d,d->v[i],tree);
      else if(d->v[i] == a) d_a = i;
    }
    
    For(j,tree->n_otu) inout[j] = 1;
    For(k,d->bip_size[d_a]) inout[d->bip_node[d_a][k]->num] = 0;
    For(j,tree->n_otu) if(inout[tree->noeud[j]->num]) tree->t_dir[d->num*dim+tree->noeud[j]->num] = d_a;
  }
  
  free(inout);
  
}

/*********************************************************/

void Fill_Dir_Table(t_tree *tree)
{
  int i,j;
  int dim;
  dim = 2*tree->n_otu-2;
  For(i,dim*dim) tree->t_dir[i] = 0;
  Free_Bip(tree);
  Alloc_Bip(tree);
  Get_Bip(tree->noeud[0],tree->noeud[0]->v[0],tree);
  Update_Dir_To_Tips(tree->noeud[0],tree->noeud[0]->v[0],tree);
  
  for(i=tree->n_otu;i<2*tree->n_otu-2;i++)
    for(j=i;j<2*tree->n_otu-2;j++)
    {
      Find_Mutual_Direction(tree->noeud[i],tree->noeud[j],
                            &(tree->t_dir[i*dim+j]),
                            &(tree->t_dir[j*dim+i]));
    }
}

/*********************************************************/

void Fast_Br_Len(t_edge *b, t_tree *tree, int approx)
{
  phydbl sum;
  phydbl *prob, *F;
  int i, j, k, site;
  phydbl v_rght;
  int dim1,dim2,dim3;
  phydbl eps_bl,old_l,new_l;
  int n_iter;
  phydbl scale_rght;

  n_iter = 0;
  dim1   = tree->mod->ns * tree->mod->n_catg;
  dim2   = tree->mod->ns ;
  dim3   = tree->mod->ns * tree->mod->ns;
  eps_bl = BL_MIN;
  
  F    = tree->triplet_struct->F_bc;
  prob = tree->triplet_struct->F_cd;

  Update_PMat_At_Given_Edge(b,tree);

  For(i,dim1*dim2) F[i] = .0;
  int *scale_down;
  phydbl *p_lk;

  if(tree->mod->whichrealmodel != HLP17){
   For(site,tree->n_pattern){
    // Joint probabilities of the states at the two ends of the t_edge
    v_rght = -1.;
    For(i,tree->mod->ns){
      For(j,tree->mod->ns){
	      For(k,tree->mod->n_catg){
          v_rght = (b->rght->tax)?((phydbl)(b->p_lk_tip_r[site*dim2+j])):(b->p_lk_rght[site*dim1+k*dim2+j]);
          scale_rght = (b->rght->tax)?(0.0):(b->sum_scale_rght[k*tree->n_pattern+site]);
          prob[dim3*k+dim2*i+j]              =
          tree->mod->pi[i]                 *
		  b->bPmat_part[tree->mod->partIndex[site]][k*dim3+i*dim2+j] *
          b->p_lk_left[site*dim1+k*dim2+i] *
          v_rght *
          POW(2.,-(b->sum_scale_left[k*tree->n_pattern+site] + scale_rght));
        }
	  }
    }
    // Scaling
    sum = .0;
    For(k,tree->mod->n_catg) For(i,tree->mod->ns) For(j,tree->mod->ns) sum += prob[dim3*k+dim2*i+j];
    For(k,tree->mod->n_catg) For(i,tree->mod->ns) For(j,tree->mod->ns) prob[dim3*k+dim2*i+j] /= sum;

    // Expected number of each pair of states
    For(i,tree->mod->ns) For(j,tree->mod->ns) For(k,tree->mod->n_catg)
    F[dim3*k+dim2*i+j] += tree->data->wght[site] * prob[dim3*k+dim2*i+j];
   }
   old_l = b->l;
   Opt_Dist_F(&(b->l),F,tree->mod);
   new_l = b->l;
   n_iter++;
   if(b->l < BL_MIN)      b->l = BL_MIN;
   else if(b->l > BL_MAX) b->l = BL_MAX;
  }
  if(!approx || tree->mod->whichmodel == HLP17){
	 // if(!approx || tree->mod->whichmodel == HLP17){
	  //includes lk_at_given_edge
    Br_Len_Brent(0.02*b->l,b->l,50.*b->l,
                 tree->mod->s_opt->min_diff_lk_local,
                 b,tree,
                 tree->mod->s_opt->brent_it_max,
                 tree->mod->s_opt->quickdirty);
  }else{
    Lk_At_Given_Edge(b,tree);
  }
  Update_PMat_At_Given_Edge(b,tree);
}

/*********************************************************/

eigen *Make_Eigen_Struct(model* mod)
{
  eigen *eig;
  int ns=mod->ns;
  
  eig              = (eigen *)mCalloc(                    1,sizeof(eigen));
  eig->size        = ns;
  eig->space       = (phydbl *)mCalloc( 2*ns*mod->n_w_catg,sizeof(phydbl));
  eig->space_int   = (int    *)mCalloc( 2*ns*mod->n_w_catg,sizeof(int));
  eig->e_val       = (phydbl *)mCalloc(   ns*mod->n_w_catg,sizeof(phydbl));
  eig->e_val_im    = (phydbl *)mCalloc(   ns*mod->n_w_catg,sizeof(phydbl));
  eig->r_e_vect    = (phydbl *)mCalloc(ns*ns*mod->n_w_catg,sizeof(phydbl));
  eig->r_e_vect_im = (phydbl *)mCalloc(ns*ns*mod->n_w_catg,sizeof(phydbl));
  eig->l_e_vect    = (phydbl *)mCalloc(ns*ns*mod->n_w_catg,sizeof(phydbl));
  eig->q           = (phydbl *)mCalloc(ns*ns*mod->n_w_catg,sizeof(phydbl));
  return eig;
}

/*********************************************************/

triplet *Make_Triplet_Struct(model *mod)
{
  int i,j,k;
  // int m;
  triplet *triplet_struct;
  
  triplet_struct                  = (triplet *)mCalloc(1,sizeof(triplet));
  triplet_struct->size            = mod->ns;
  triplet_struct->pi_bc           = (phydbl *)mCalloc(mod->ns,sizeof(phydbl ));
  triplet_struct->pi_cd           = (phydbl *)mCalloc(mod->ns,sizeof(phydbl ));
  triplet_struct->pi_bd           = (phydbl *)mCalloc(mod->ns,sizeof(phydbl ));
  triplet_struct->F_bd            = (phydbl *)mCalloc(mod->ns*mod->ns,sizeof(phydbl));
  triplet_struct->p_one_site      = (phydbl ***)mCalloc(mod->ns,sizeof(phydbl **));
  For(i,mod->ns)
  {
    triplet_struct->p_one_site[i]          = (phydbl **)mCalloc(mod->ns,sizeof(phydbl *));
    For(j,mod->ns)
    triplet_struct->p_one_site[i][j]     = (phydbl  *)mCalloc(mod->ns,sizeof(phydbl ));
  }
  
  triplet_struct->sum_p_one_site  = (phydbl ***)mCalloc(mod->ns,sizeof(phydbl **));
  For(i,mod->ns)
  {
    triplet_struct->sum_p_one_site[i]      = (phydbl **)mCalloc(mod->ns,sizeof(phydbl *));
    For(j,mod->ns)
    triplet_struct->sum_p_one_site[i][j] = (phydbl  *)mCalloc(mod->ns,sizeof(phydbl ));
  }
  
  triplet_struct->F_bc            = (phydbl *)mCalloc(mod->ns*mod->ns*mod->n_catg,sizeof(phydbl));
  triplet_struct->F_cd            = (phydbl *)mCalloc(mod->ns*mod->ns*mod->n_catg,sizeof(phydbl));
  triplet_struct->core            = (phydbl ****)mCalloc(mod->n_catg,sizeof(phydbl ***));
  For(k,mod->n_catg)
  {
    triplet_struct->core[k]            = (phydbl ***)mCalloc(mod->ns,sizeof(phydbl **));
    For(i,mod->ns)
    {
      triplet_struct->core[k][i]       = (phydbl **)mCalloc(mod->ns,sizeof(phydbl *));
      For(j,mod->ns)
      triplet_struct->core[k][i][j]    = (phydbl  *)mCalloc(mod->ns,sizeof(phydbl ));
    }
  }
  
  triplet_struct->eigen_struct    = (eigen *)Make_Eigen_Struct(mod);
  
  triplet_struct->mod             = mod;
  
  return triplet_struct;
}

/*********************************************************/

phydbl Triple_Dist(t_node *a, t_tree *tree, int approx)
{
  if(a->tax) return UNLIKELY;
  else
  {
	Update_PMat_At_Given_Edge(a->b[1],tree);
	Update_PMat_At_Given_Edge(a->b[2],tree);
	if(tree->mod->whichrealmodel==HLP17){
	  	t_edge* ance = a->anc_edge;
    	int dir1 = -1;
		int dir2 = -1;
		int i=0;
		For(i,3){
			if(a->b[i]->num != ance->num){
				if(dir1 == -1) dir1=i;
				else dir2 = i;
			}
		}
		//printf("%d %d %d\n",ance->anc_node->num,a->b[dir1]->anc_node->num,a->b[dir2]->anc_node->num);

    	Update_P_Lk(tree,ance,a);
    	if(ance->anc_node->num != tree->mod->startnode) Fill_UPP_single(tree,ance);
    	else 											Fill_UPP_root(tree,ance);

    	Fast_Br_Len(ance,tree,approx);
    	Fill_UPP_single(tree,a->b[dir1]);
    	Fill_UPP_single(tree,a->b[dir2]);

    	Update_P_Lk(tree,a->b[dir1],a);
    	Fast_Br_Len(a->b[dir1],tree,approx);
    	Fill_UPP_single(tree,a->b[dir2]);
    
    	Update_P_Lk(tree,a->b[dir2],a);
    	Fast_Br_Len(a->b[dir2],tree,approx);
    	Fill_UPP_single(tree,a->b[dir1]);

    	Update_P_Lk(tree,a->b[dir1],a);
    	Update_P_Lk(tree,a->b[dir2],a);
    	Update_P_Lk(tree,ance,a);
  	}else{
    	Update_P_Lk(tree,a->b[0],a);
    	Fast_Br_Len(a->b[0],tree,approx);
    	/*       Br_Len_Brent (BL_MAX, a->b[0]->l,BL_MIN, 1.e-10,a->b[0],tree,50,0); */
    	Update_P_Lk(tree,a->b[1],a);
    	Fast_Br_Len(a->b[1],tree,approx);
    	/*       Br_Len_Brent (BL_MAX, a->b[1]->l,BL_MIN, 1.e-10,a->b[1],tree,50,0); */
    	Update_P_Lk(tree,a->b[2],a);
    	Fast_Br_Len(a->b[2],tree,approx);
    	/*       Br_Len_Brent (BL_MAX, a->b[2]->l,BL_MIN, 1.e-10,a->b[2],tree,50,0); */
    	Update_P_Lk(tree,a->b[1],a);
    	Update_P_Lk(tree,a->b[0],a);
    }

  }
  
  return tree->c_lnL;
  
}

/*********************************************************/

void Fix_All(t_tree *tree)
{
  int i;
  
  tree->mod->pinvar_old = tree->mod->pinvar;
  tree->mod->alpha_old  = tree->mod->alpha;
  tree->mod->kappa_old  = tree->mod->kappa;
  tree->mod->lambda_old = tree->mod->lambda;
  
  for(i=tree->n_otu;i<2*tree->n_otu-2;i++)
  {
    tree->noeud[i]->b[0]->l_old = tree->noeud[i]->b[0]->l;
    tree->noeud[i]->b[1]->l_old = tree->noeud[i]->b[1]->l;
    tree->noeud[i]->b[2]->l_old = tree->noeud[i]->b[2]->l;
  }
}

/*********************************************************/

void Record_Br_Len(phydbl *where, t_tree *tree)
{

	// upAllPmats(tree); //added by Ken 19/2/2017
  int i;
  
  if(!where)
  {
    For(i,2*tree->n_otu-3) tree->t_edges[i]->l_old = tree->t_edges[i]->l;
  }
  else
  {
    For(i,2*tree->n_otu-3) where[i] = tree->t_edges[i]->l;
  }
}

/*********************************************************/

void Restore_Br_Len(phydbl *from, t_tree *tree)
{
  int i;
  
  if(!from)
  {
    For(i,2*tree->n_otu-3) tree->t_edges[i]->l = tree->t_edges[i]->l_old;
  }
  else
  {
    For(i,2*tree->n_otu-3) tree->t_edges[i]->l = from[i];
  }
}


/*********************************************************/

void Random_Tree(t_tree *tree)
{
  int *is_available,*list_of_nodes;
  int i,node_num,step,n_available;
  phydbl min_edge_len;
  
  min_edge_len = 1.E-3;
  
  PhyML_Printf("\n. Randomising the tree...\n");
  
  is_available  = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));
  list_of_nodes = (int *)mCalloc(tree->n_otu,    sizeof(int));
  
  For(i,tree->n_otu) is_available[i]  = 1;
  For(i,tree->n_otu) list_of_nodes[i] = i;
  
  step = 0;
  do
  {
    /*       node_num = (int)RINT(rand()/(phydbl)(RAND_MAX+1.0)*(tree->n_otu-1-step)); */
    node_num = Rand_Int(0,tree->n_otu-1-step);
    node_num = list_of_nodes[node_num];
    is_available[node_num] = 0;
    For(i,tree->n_otu) list_of_nodes[i] = -1;
    n_available = 0;
    For(i,2*tree->n_otu-2) if(is_available[i]) {list_of_nodes[n_available++] = i;}
    
    tree->noeud[node_num]->v[0] = tree->noeud[tree->n_otu+step];
    tree->noeud[tree->n_otu+step]->v[1] = tree->noeud[node_num];
    
    /*       node_num = (int)RINT(rand()/(phydbl)(RAND_MAX+1.0)*(tree->n_otu-2-step)); */
    node_num = Rand_Int(0,tree->n_otu-2-step);
    node_num = list_of_nodes[node_num];
    is_available[node_num] = 0;
    For(i,tree->n_otu) list_of_nodes[i] = -1;
    n_available = 0;
    For(i,2*tree->n_otu-2) if(is_available[i]) {list_of_nodes[n_available++] = i;}
    
    tree->noeud[node_num]->v[0] = tree->noeud[tree->n_otu+step];
    tree->noeud[tree->n_otu+step]->v[2] = tree->noeud[node_num];
    
    is_available[tree->n_otu+step] = 1;
    For(i,tree->n_otu) list_of_nodes[i] = -1;
    n_available = 0;
    For(i,2*tree->n_otu-2) if(is_available[i]) list_of_nodes[n_available++] = i;
    
    step++;
  }while(step < tree->n_otu-2);
  
  tree->noeud[list_of_nodes[0]]->v[0] = tree->noeud[list_of_nodes[1]];
  tree->noeud[list_of_nodes[1]]->v[0] = tree->noeud[list_of_nodes[0]];
  
  tree->num_curr_branch_available = 0;
  Connect_Edges_To_Nodes_Recur(tree->noeud[0],tree->noeud[0]->v[0],tree);
  
  For(i,2*tree->n_otu-3) if(tree->t_edges[i]->l < min_edge_len) tree->t_edges[i]->l = min_edge_len;
  
  Fill_Dir_Table(tree);
  Update_Dirs(tree);
  
  free(is_available);
  free(list_of_nodes);
}

/*********************************************************/



void Print_Settings(option *io)
{
  int answer;
  char *s,*r;
  
  s = (char *)mCalloc(100,sizeof(char));
  r = (char *)mCalloc(100,sizeof(char));
  PhyML_Printf("\n\n\n");
  PhyML_Printf("\n\n"); 
  PhyML_Printf("                               ..................                               \n");
  PhyML_Printf("oooooooooooooooooooooooo       CURRENT   SETTINGS       oooooooooooooooooooooooo\n");
  PhyML_Printf("                               ..................                               \n");
  
  PhyML_Printf("\n. Sequence filename:\t\t\t\t %s", Basename(io->in_align_file));
  
  PhyML_Printf("\n. Number of taxa:\t\t\t\t %d", io->n_otu); //!< Added by Marcelo.
  if(io->datatype == CODON ) PhyML_Printf("\n. Sequence length:\t\t\t\t %d", io->init_len/3); //!< Added by Marcelo.
  else PhyML_Printf("\n. Sequence length:\t\t\t\t %d", io->init_len); //!< Added by Marcelo.
  
  if(io->datatype == CODON)  strcpy(s,"Codon"); //!< Added by Marcelo.                                                       
  else if(io->datatype == NT)       strcpy(s,"DNA");
  else if(io->datatype == AA)        strcpy(s,"Amino Acids");
  else strcpy(s,"generic");
  
  PhyML_Printf("\n. Data type:\t\t\t\t\t %s",s);
  if(io->datatype == CODON) 
  {
    switch(io->mod->genetic_code)
    {
      case   STANDARD: strcpy(s,"Standard\0");break;
      case       TVMC: strcpy(s,"Vertebrate Mitochondrial\0"); break; 
      case       TYMC: strcpy(s,"Yeast Mitochondrial\0"); break;
      case THMPCMCMSC: strcpy(s,"Mold, Protozoan, Coelenterate ...\0"); break; 
      case      THIMC: strcpy(s,"Invertebrate Mitochondrial\0"); break; 
      case    THCDHNC: strcpy(s,"Ciliate, Dasycladacean ...\0"); break;
      case     THEFMC: strcpy(s,"Echinoderm and Flatworm Mit. ...\0"); break;
      case      THENC: strcpy(s,"Euplotid Nuclear\0"); break;
      case    THBAPPC: strcpy(s,"Bacterial, Archaeal and Plant ...\0"); break;
      case     THAYNC: strcpy(s,"Alternative Yeast Nuclear\0"); break;
      case      THAMC: strcpy(s,"Ascidian Mitochondrial\0"); break;
      case     THAFMC: strcpy(s,"Alternative Flatworm Mit. ...\0"); break;
      case       BLNC: strcpy(s,"Blepharisma Nuclear\0"); break;
      case       CHMC: strcpy(s,"Chlorophycean Mitochondrial\0"); break;
      case       TRMC: strcpy(s,"Trematode Mitochondrial\0"); break;
      case      SCOMC: strcpy(s,"Scenedesmus Obliquus Mit. ...\0"); break;
      case       THMC: strcpy(s,"Thraustochytrium Mitochondrial\0"); break;
      default : Warn_And_Exit("Genetic code not implemented."); break; 
        break;
    }
    PhyML_Printf("\n. Alphabet size:\t\t\t\t %d",io->mod->ns); 
    PhyML_Printf("\n. Genetic code:\t\t\t\t\t %s",s);
    
  }
  else 
    PhyML_Printf("\n. Alphabet size:\t\t\t\t %d",io->mod->ns);
  
  PhyML_Printf("\n. Sequence format:\t\t\t\t %s", io->interleaved ? "Interleaved": "Sequential");
  PhyML_Printf("\n. Number of data sets:\t\t\t\t %d", io->n_data_sets);
  
  PhyML_Printf("\n. Number of bootstrapped data sets:\t\t %d", io->mod->bootstrap);
  
  if (io->mod->bootstrap > 0)
    PhyML_Printf("\n. Compute aLRT:\t\t\t\t\t No");
  else
  {
    if(io->ratio_test == 1)
      PhyML_Printf("\n. Compute aLRT:\t\t\t\t\t Yes (aLRT statistics)");
    else if(io->ratio_test == 2)
      PhyML_Printf("\n. Compute aLRT:\t\t\t\t\t Yes (Chi^2-based br. supports)");
    else if(io->ratio_test == 3)
      PhyML_Printf("\n. Compute aLRT:\t\t\t\t\t Yes (Min[SH-like, Chi^2-based]");
    else if(io->ratio_test == 4)
      PhyML_Printf("\n. Compute aLRT:\t\t\t\t\t Yes (SH-like branch supports)");
    else if(io->ratio_test == 5)
      PhyML_Printf("\n. Compute aLRT:\t\t\t\t\t Yes (aBayes branch supports)");
    
  }
  if(io->datatype == CODON)
  {
    switch(io->mod->freq_model)
    {
      case F1XSENSECODONS: 
        strcpy(s,"F1x"); 
        sprintf(r,"%d",io->mod->ns);
        strcat(s,r);
        break;
      case F1X4: strcpy(s,"F1x4");
        break;
      case F3X4: strcpy(s,"F3x4");
        break;
      case CF3X4: strcpy(s,"CF3x4");
        break;
      default: strcpy(s,"\0");break;
    }
    if(io->mod->whichmodel==MGECMUSR||io->mod->whichmodel==YAPECMUSR||io->mod->whichmodel==GYECMUSR ||io->mod->whichmodel==GYECMK07 || io->mod->whichmodel==GYECMS05 || io->mod->whichmodel==YAPECMK07 || io->mod->whichmodel==YAPECMS05 || io->mod->whichmodel==MGECMUSRWK||io->mod->whichmodel==YAPECMUSRWK||io->mod->whichmodel==GYECMUSRWK ||io->mod->whichmodel==GYECMK07WK || io->mod->whichmodel==GYECMS05WK || io->mod->whichmodel==YAPECMK07WK || io->mod->whichmodel==YAPECMS05WK || io->eq_freq_handling == MODEL)
    {
      strcpy(r,"Model-Defined");
    }
    else{
      if(io->mod->s_opt->opt_state_freq) strcpy(r,"Optimized");
      else if(io->mod->s_opt->user_state_freq) strcpy(r,"User-Defined");
      else  strcpy(r,"Empirical");
    }
    PhyML_Printf("\n. Model type:\t\t\t\t\t %s", io->mod->modelname);
    if(io->mod->pcaModel) PhyML_Printf("\n. Number of PCs:\t\t\t\t %d", io->mod->npcs);
    PhyML_Printf("\n. Equilibrium frequencies model:\t\t %s %s", r,s);
    
    if(io->modeltypeOpt == HLP17){
        PhyML_Printf("\n. Motifs:\t\t\t\t\t %s", io->mod->motifstring);
        PhyML_Printf("\n. h parameter(s):\t\t\t\t %s", io->mod->hotnessstring);
        PhyML_Printf("\n. Root ID:\t\t\t\t\t %s", io->mod->rootname);
    }
    
    
    
    
  }else PhyML_Printf("\n. Model name:\t\t\t\t\t %s", io->mod->modelname);
  
  
  if (io->datatype == NT)
  {
    if ((io->mod->whichmodel == K80)  ||
        (io->mod->whichmodel == HKY85)||
        (io->mod->whichmodel == F84)  ||
        (io->mod->whichmodel == TN93))
    {
      if (io->mod->s_opt->opt_kappa)
        PhyML_Printf("\n. Ts/tv ratio:\t\t\t\t\t Estimated");
      else
        PhyML_Printf("\n. Ts/tv ratio:\t\t\t\t\t %f", io->mod->kappa);
    }
  }
  
  if (io->mod->s_opt->opt_pinvar)
    PhyML_Printf("\n. Prop. of invariable sites:\t\t\t Estimated");
  else
    PhyML_Printf("\n. Prop. of invariable sites:\t\t\t %.2f", io->mod->pinvar);
  
  if(io->mod->n_w_catg==1 && io->mod->n_catg>1) //!< Added By Marcelo. 
  {
    PhyML_Printf("\n. Number of subst. rate categs:\t\t\t %d", io->mod->n_catg);
    if(io->mod->s_opt->opt_alpha) PhyML_Printf("\n. Gamma distribution parameter:\t\t\t Estimated");
    else PhyML_Printf("\n. Gamma distribution parameter:\t\t\t %f", io->mod->alpha);
    PhyML_Printf("\n. 'Middle' of each rate class:\t\t\t %s",(io->mod->gamma_median)?("Median"):("Mean"));
  }
  else
  {
    PhyML_Printf("\n. Number of subst. rate categs:\t\t\t 1");
  }
  
  if(io->datatype==CODON) //!< Added By Marcelo.
  {
    if(io->mod->n_w_catg==1) PhyML_Printf("\n. Number of dn/ds rate categs:\t\t\t %d - (M0)", io->mod->n_w_catg);
    else if(io->mod->omegaSiteVar==DGAMMAK) 
    {
      PhyML_Printf("\n. Number of dn/ds rate ratio categories:\t %d - Discrete Gamma model (M5)", io->mod->n_w_catg);
      if(io->mod->s_opt->opt_alpha) PhyML_Printf("\n. Gamma distribution parameter:\t\t\t Estimated");
      else PhyML_Printf("\n. Gamma distribution parameter:\t\t\t %f and %f", io->mod->alpha, io->mod->beta);
      PhyML_Printf("\n. 'Middle' of each rate class:\t\t\t %s",(io->mod->gamma_median)?("Median"):("Mean")); 
    }
    else PhyML_Printf("\n. Number of dn/ds rate ratio categories:\t %d - Discrete model (M3)", io->mod->n_w_catg);
  }
  
  if(io->datatype == AA)
  {
    char ll[30];
    if(io->mod->s_opt->user_state_freq) { strcpy(ll,"User-Defined\0"); }
    else if ((io->mod->s_opt->opt_state_freq==NO)&&(io->mod->s_opt->opt_state_freq_AAML==NO)){ strcpy(ll,"Model-Defined\0"); }
    else if ((io->mod->s_opt->opt_state_freq==YES)&&(io->mod->s_opt->opt_state_freq_AAML==NO)) { strcpy(ll,"Empirical\0");} //herere
    else if ((io->mod->s_opt->opt_state_freq==YES)&&(io->mod->s_opt->opt_state_freq_AAML==YES)){ strcpy(ll,"Optimized\0"); }
    PhyML_Printf("\n. AA equilibrium frequencies:\t\t\t %s", ll);
  }
  else if(io->datatype == NT)
  {
    if((io->mod->whichmodel != JC69) &&
       (io->mod->whichmodel != K80)  &&
       (io->mod->whichmodel != F81))
    {
      if(!io->mod->s_opt->user_state_freq)
	    {
	      PhyML_Printf("\n. NT equilibrium frequencies:\t\t\t %s", (io->mod->s_opt->opt_state_freq) ? ("Optimized"):("Empirical"));
	    }
      else
	    {
	      PhyML_Printf("\n. Nucleotide equilibrium frequencies:\t %s","User-defined");
	    }
    }
  }
  
  PhyML_Printf("\n. Optimise tree topology:\t\t\t %s", (io->mod->s_opt->opt_topo) ? "Yes": "No");
  
  if(io->datatype==CODON)
  {
    PhyML_Printf("\n. Heuristic during tree search:\t\t\t %s", (io->opt_heuristic_manuel) ? "Modified": "Original");
    PhyML_Printf("\n. Parameter optimization strategy:\t\t %s", (io->mod->s_opt->opt_method) ? "Single variate (BRENT)": "Multi-variate (BFGS)");
  }
    
   switch(io->in_tree)
   {
	case 0: { strcpy(s,"BioNJ");     break; }
	case 1: { strcpy(s,"parsimony"); break; }
	case 2: { strcpy(s,"user tree"); break; }
	default: break;
   }
   
   if(io->datatype==CODON)//!< Added by Marcelo
   {
     
      if((io->in_tree==2 && io->nobl) || io->in_tree!=2)
      {
	
	switch(io->init_DistanceTreeCD)
	{
	  case NUCLEO: strcat(s," + JC69\0");break;
	  case KOSI07: strcat(s," + GYECMK07\0");break;
	  case SCHN05: strcat(s," + GYECMS05\0");break;
	  case ECMUSR: strcat(s," + ECMUSR\0");break;
	  default: break;
	}
      }
    }
    
  if(io->mod->s_opt->opt_topo)
  {
    if(io->mod->s_opt->topo_search == NNI_MOVE) PhyML_Printf("\n. Tree topology search:\t\t\t\t NNIs");
    else if(io->mod->s_opt->topo_search == SPR_MOVE) PhyML_Printf("\n. Tree topology search:\t\t\t\t SPRs");
    else if(io->mod->s_opt->topo_search == BEST_OF_NNI_AND_SPR) PhyML_Printf("\n. Tree topology search:\t\t\t\t Best of NNIs and SPRs");
  
    PhyML_Printf("\n. Starting tree:\t\t\t\t %s",s);
    
    PhyML_Printf("\n. Add random input tree:\t\t\t %s", (io->mod->s_opt->random_input_tree) ? "Yes": "No");
    if(io->mod->s_opt->random_input_tree)
      PhyML_Printf("\n. Number of random starting trees:\t %d", io->mod->s_opt->n_rand_starts);	
    
  }
  else
  {
    if(!io->mod->s_opt->random_input_tree)
      PhyML_Printf("\n. Evaluated tree:\t\t\t\t %s",s);
  }
  
  PhyML_Printf("\n. Optimise branch lengths:\t\t\t %s", (io->mod->s_opt->opt_bl) ? "Yes": "No");
  
  answer = 0;
  if(io->mod->s_opt->opt_alpha      ||
     io->mod->s_opt->opt_kappa      ||
     io->mod->s_opt->opt_lambda     ||
     io->mod->s_opt->opt_pinvar     ||
     io->mod->s_opt->opt_rr         ||
     io->mod->s_opt->opt_omega      || //!< Added by Marcelo.
     io->mod->s_opt->opt_beta       || //!< Added by Marcelo.
     io->mod->s_opt->opt_prob_omega || //!< Added by Marcelo.
     io->mod->s_opt->opt_alphaCD) answer = 1;
  
  PhyML_Printf("\n. Optimise subst. model params:\t\t\t %s", (answer) ? "Yes": "No");
  if(io->expm==SSPADE) 
  {
    PhyML_Printf("\n. Likelihood o.f. for subst. mod. param. opt.:\t %s", "Pade approximation\0");
  }
  else 
  {
    if(answer) PhyML_Printf("\n. Matrix exponential:\t\t\t\t %s", io->heuristicExpm ? "Taylor approximation": "Eigenvalue problem");
  }
  PhyML_Printf("\n. Run ID:\t\t\t\t\t %s", (io->append_run_ID) ? (io->run_id_string): ("None"));
  PhyML_Printf("\n. Random seed:\t\t\t\t\t %d", io->r_seed);
  
  if(io->datatype==CODON)
  {
#ifdef BLAS_OMP
    
    PhyML_Printf("\n. Code optimization:\t\t\t\t BLAS, LAPACK and OpenMP");
    
#elif defined BLAS
    
    PhyML_Printf("\n. Code optimization:\t\t\t\t BLAS and LAPACK");
    
#elif defined OMP 
    
    PhyML_Printf("\n. Code optimization:\t\t\t\t OpenMP");
    
#else
    
    PhyML_Printf("\n. Code optimization:\t\t\t\t None");
    
#endif
  }
  
  PhyML_Printf("\n. Version:\t\t\t\t\t %s", VERSION);
  PhyML_Printf("\n\noooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
  PhyML_Printf("\n\n");
  fflush(NULL);
  
  free(s);free(r);
}

/*********************************************************/

void Fill_Missing_Dist(matrix *mat)
{
  int i,j;
  For(i,mat->n_otu)
  {
    for(j=i+1;j<mat->n_otu;j++)
    {
      if(i != j)
	    {
	      if(mat->dist[i][j] < .0) 
        {
          Fill_Missing_Dist_XY(i,j,mat);
          mat->dist[j][i] = mat->dist[i][j];
        }
	    }
    }
  }
}

/*********************************************************/

void Fill_Missing_Dist_XY(int x, int y, matrix *mat)
{
  
  int i,j;
  phydbl *local_mins,**S1S2;
  int cpt;
  int pos_best_estimate;
  phydbl min_crit, curr_crit;
  
  local_mins = (phydbl *)mCalloc(mat->n_otu*mat->n_otu,sizeof(phydbl ));
  S1S2       = (phydbl **)mCalloc(mat->n_otu*mat->n_otu,sizeof(phydbl *));
  For(i,mat->n_otu*mat->n_otu) S1S2[i] = (phydbl *)mCalloc(2,sizeof(phydbl));
  
  cpt = 0;
  For(i,mat->n_otu)
  {
    if((mat->dist[i][x] > .0) && (mat->dist[i][y] > .0))
    {
      For(j,mat->n_otu)
	    {
	      if((mat->dist[j][x] > .0) && (mat->dist[j][y] > .0))
        {
          if((i != j) && (i != x) && (i != y) && (j != x) && (j != y))
          {
            S1S2[cpt][0] = MIN(mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j] , mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]);
            S1S2[cpt][1] = MAX(mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j] , mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]);
            cpt++;
          }
        }
	    }
    }
  }
  
  Qksort_Matrix(S1S2,0,0,cpt-1);
  
  local_mins[0] = S1S2[0][1];
  for(i=1;i<cpt;i++) local_mins[i] = (i*local_mins[i-1] + S1S2[i][1])/(phydbl)(i+1);
  
  pos_best_estimate = 0;
  min_crit = curr_crit = BIG;
	
  For(i,cpt-1)
  {
    if((local_mins[i] < S1S2[i+1][0]) && (local_mins[i] > S1S2[i][0]))
    {
      curr_crit = Least_Square_Missing_Dist_XY(x,y,local_mins[i],mat);
      if(curr_crit < min_crit)
	    {
	      min_crit = curr_crit;
	      pos_best_estimate = i;
	    }
    }
  }
  
  mat->dist[x][y] = local_mins[pos_best_estimate];
  mat->dist[y][x] = mat->dist[x][y];
  
  For(i,mat->n_otu*mat->n_otu) free(S1S2[i]);
  free(S1S2);
  free(local_mins);
}

/*********************************************************/

phydbl Least_Square_Missing_Dist_XY(int x, int y, phydbl dxy, matrix *mat)
{
  int i,j;
  phydbl fit;
  
  fit = .0;
  For(i,mat->n_otu)
  {
    if((mat->dist[i][x] > .0) && (mat->dist[i][y] > .0))
    {
      For(j,mat->n_otu)
	    {
	      if((mat->dist[j][x] > .0) && (mat->dist[j][y] > .0))
        {
          if((i != j) && (i != x) && (i != y) && (j != x) && (j != y))
          {
            if(dxy < MIN(mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j] , mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]))
            {
              fit += POW((mat->dist[i][x] + mat->dist[j][y]) - (mat->dist[i][y] + mat->dist[j][x]),2);
            }
            else if((mat->dist[i][x] + mat->dist[j][y]) < (mat->dist[i][y] + mat->dist[j][x]))
            {
              fit += POW(dxy - (mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]),2);
            }
            else
            {
              fit += POW(dxy - (mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j]),2);
            }
          }
        }
	    }
    }
  }
  return fit;
}

/*********************************************************/

token * Print_Banner_Small(token *t) {                                      
    t = Emit_Out_Token(t, "Header", "Header", COMMENTTKN, TFSTRING, "%s", "oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo");
    t = Emit_Out_Token(t, "Header", "Header", COMMENTTKN, TFSTRING, "%s%s ---", "                      --- IgPhyML ", VERSION);
    t = Emit_Out_Token(t, "Header", "Header", COMMENTTKN, TFSTRING, "%s", "               Kenneth B Hoehn, Gerton Lunter, Oliver G Pybus");
    t = Emit_Out_Token(t, "Header", "Header", COMMENTTKN, TFSTRING, "%s", "                                ");
    t = Emit_Out_Token(t, "Header", "Header", COMMENTTKN, TFSTRING, "%s", "                                ");
    t = Emit_Out_Token(t, "Header", "Header", COMMENTTKN, TFSTRING, "%s", "                         Based off of codonPhyML");
    t = Emit_Out_Token(t, "Header", "Header", COMMENTTKN, TFSTRING, "%s", "Marcelo S Zanetti, Stefan Zoller, Manuel Gil, Louis du Plessis, Maria Anisimova");
    t = Emit_Out_Token(t, "Header", "Header", COMMENTTKN, TFSTRING, "%s", "                 http://sourceforge.net/projects/codonphyml/");
    t = Emit_Out_Token(t, "Header", "Header", COMMENTTKN, TFSTRING, "%s", "oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo");
    return(t);
}

/*********************************************************/



void Check_Memory_Amount(t_tree *tree)
{
  /* Rough estimate of the amount of memory that has to be used */
  
  long long int nbytes;
  int n_otu;
  model *mod;
  
  mod    = tree->mod;
  n_otu  = tree->io->n_otu;
  nbytes = 0;
  
  /* Partial Pars */
  nbytes += (2*n_otu-3) * 2 * tree->data->crunch_len * sizeof(int);
  nbytes += (2*n_otu-3) * 2 * tree->data->crunch_len * sizeof(unsigned int);
  nbytes += (2*n_otu-3) * 2 * tree->data->crunch_len * mod->ns * sizeof(int);
  
  /* Pmat */
  nbytes += (2*n_otu-3) * mod->n_catg * mod->ns * mod->ns * sizeof(phydbl);
  
  /* Partial Lk */
  nbytes += ((2*n_otu-3) * 2 - tree->n_otu) * tree->data->crunch_len * mod->n_catg * mod->ns * sizeof(phydbl);
  
  /* Scaling factors */
  nbytes += ((2*n_otu-3) * 2 - tree->n_otu) * tree->data->crunch_len * sizeof(int);
  
  if(tree->io->datatype==CODON) { //!<Added by Marcelo.
    nbytes += n_otu * tree->data->crunch_len * (mod->ns) * sizeof(int); //!< Codon alternatives.  
    nbytes += sizeof(tree->mod->eigen)+4*mod->ns*mod->ns*sizeof(phydbl)*mod->n_w_catg+64*64*
              sizeof(tree->io->structTs_and_Tv)+mod->ns*mod->ns*sizeof(phydbl)+2*mod->num_base_freq*
              mod->num_base_freq*sizeof(phydbl)+2*mod->n_w_catg*sizeof(phydbl)+mod->n_w_catg*
              mod->ns*mod->ns*sizeof(phydbl);
    nbytes += (sizeof(tree->triplet_struct)+ 3*mod->ns*sizeof(phydbl )+ mod->ns*mod->ns*sizeof(phydbl)+ 
              2*(mod->ns* sizeof(phydbl **)+mod->ns*mod->ns*sizeof(phydbl *)+mod->ns*mod->ns*mod->ns*
              sizeof(phydbl ))+2*mod->ns*mod->ns*mod->n_catg*sizeof(phydbl)+ mod->n_catg*
              sizeof(phydbl ***) + mod->n_catg*mod->ns*sizeof(phydbl **) + mod->n_catg*mod->ns*
              mod->ns*sizeof(phydbl *)+mod->n_catg*mod->ns*mod->ns*mod->ns*sizeof(phydbl ))+
              2*(sizeof(eigen)+6*mod->ns*sizeof(phydbl)*mod->n_w_catg + 4*mod->ns*mod->ns*
              sizeof(phydbl)*mod->n_w_catg);
    if(tree->io->expm==SSPADE) { 
        nbytes +=12*mod->ns*mod->ns*sizeof(phydbl)*mod->n_w_catg;
    }
    if(tree->io->heuristicExpm) { 
        nbytes +=15*mod->ns*mod->ns*sizeof(phydbl)*mod->n_w_catg;
    }
  }
  
  if(((phydbl)nbytes/(1.E+06)) > 256.)
  {
    char answer;
    PhyML_Printf("\n. WARNING: this analysis requires at least %.0f MB of memory space.\n",(phydbl)nbytes/(1.E+06));
#ifndef BATCH
//     if (! tree->io->quiet) {
//       PhyML_Printf("\n. Do you really want to proceed? [Y/n] ");
//       if(scanf("%c", &answer))
//       {
//         if(answer == '\n') answer = 'Y';
//         else if(answer == 'n' || answer == 'N') Warn_And_Exit("\n");
//         else getchar();
//       }
//       else
//       {
//         Warn_And_Exit("\n\n");
//       }
//     }  //COMMENTED OUT BY MARCELO ... 24.07.2014
#endif
  }
  else if(((phydbl)nbytes/(1.E+06)) > 100.)
  {
    if(!tree->io->quiet) PhyML_Printf("\n. WARNING: this analysis will use at least %.0f Mo of memory space...\n",(phydbl)nbytes/(1.E+06));
  }
  else if(((phydbl)nbytes/(1.E+06)) > 1.)
  {
    if(!tree->io->quiet) PhyML_Printf("\n. This analysis requires at least %.0f Mo of memory space.\n",(phydbl)nbytes/(1.E+06));
  }
}

/*********************************************************/

int Get_State_From_P_Pars(short int *p_pars, int pos, t_tree *tree)
{
  int i;
  For(i,tree->mod->ns) if(p_pars[pos+i] > .0) return i;
  return -1;
}

/*********************************************************/

void Print_Lk(t_tree *tree, char *string)
{
#if defined OMP || defined BLAS_OMP
  
  tree->t_current=omp_get_wtime();
  
#else
  
  time(&(tree->t_current));
  
#endif
  PhyML_Printf("\n. (%5d sec) [%15.4f] %s",(int)(tree->t_current-tree->t_beg),tree->c_lnL,string);
#ifndef QUIET 
  fflush(NULL);
#endif
}

/*********************************************************/

void Print_Pars(t_tree *tree)
{
#if defined OMP || defined BLAS_OMP
  
  tree->t_current=omp_get_wtime();
  
#else
  
  time(&(tree->t_current));
  
#endif
  PhyML_Printf("\n. (%5d sec) [%5d]",(int)(tree->t_current-tree->t_beg),tree->c_pars);
#ifndef QUIET
  fflush(NULL);
#endif
}

/*********************************************************/

void Check_Dirs(t_tree *tree)
{
  int i;
  
  For(i,2*tree->n_otu-3)
  {
    if(!tree->t_edges[i]->left->tax)
    {
      if(tree->t_edges[i]->left->v[tree->t_edges[i]->l_v1]->num <
         tree->t_edges[i]->left->v[tree->t_edges[i]->l_v2]->num)
	    {
	      PhyML_Printf("\n. Edge %d ; v1=%d v2=%d",
                     tree->t_edges[i]->num,
                     tree->t_edges[i]->left->v[tree->t_edges[i]->l_v1]->num,
                     tree->t_edges[i]->left->v[tree->t_edges[i]->l_v2]->num);
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
    }
    
    if(!tree->t_edges[i]->rght->tax)
    {
      if(tree->t_edges[i]->rght->v[tree->t_edges[i]->r_v1]->num <
         tree->t_edges[i]->rght->v[tree->t_edges[i]->r_v2]->num)
	    {
	      PhyML_Printf("\n. Edge %d ; v3=%d v4=%d",
                     tree->t_edges[i]->num,
                     tree->t_edges[i]->rght->v[tree->t_edges[i]->r_v1]->num,
                     tree->t_edges[i]->rght->v[tree->t_edges[i]->r_v2]->num);
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
    }
  }
}

/*********************************************************/

void Warn_And_Exit(char *s)
{
  PhyML_Fprintf(stdout,"%s",s);
  fflush(NULL);
#ifndef BATCH
  //  if (! tree->io->quiet) {
  char c;
  PhyML_Fprintf(stdout,"\n. Type enter to exit.\n");
  if(!fscanf(stdin,"%c",&c)) Exit("");
  //  }
#endif
  Exit("\n");
}

/*********************************************************/

void Read_userRatesAndFreqs( phydbl *daa, phydbl *pi, int ns, FILE *fp ) { //! Added by Marcelo.
  int i,j;
  phydbl sum;
  double val;
  
  For( i, ns ) {
    for( j = 0; j < i; j++ ) {
      if( !fscanf( fp, "%lf", &val ) ) Exit("\n");
      daa[i*ns + j] = (phydbl) val;
      daa[j*ns + i] = daa[i*ns + j];
    }
  }
  
  For( i, ns ) { 
    if( !fscanf( fp, "%lf", &val ) ) Exit("\n");
    pi[i] = (phydbl) val;
  }
  sum = .0;
  For( i, ns ) sum += pi[i];
  if( FABS( sum - 1. ) > 1.E-06 ) {
    PhyML_Printf("\n. Scaling codon frequencies...\n");
    For( i, ns ) pi[i] /= sum;
  }
}
/*********************************************************/
void Read_userRatesAndFreqsMG(phydbl *daa, phydbl *pi, phydbl *bfreqs, int nbfreq, int ns, FILE *fp) //! Added by Marcelo.
{
  int i,j;
  // phydbl sum;
  double val;
  
  For(i,ns)
  {
    for(j=0;j<i;j++)
    {
      if(!fscanf(fp,"%lf",&val)) Exit("\n");
      daa[i*ns + j] = (phydbl)val;
      daa[j*ns + i] = daa[i*ns + j];
    }
  }
  
  For(i,ns) 
  { 
    if(!fscanf(fp,"%lf",&val)) Exit("\n");
    pi[i] = (phydbl)val;
  }
  
  Scale_freqs(pi, ns);
  
  For(i,nbfreq) 
  { 
    if(!fscanf(fp,"%lf",&val)) Exit("\n");
    bfreqs[i] = (phydbl)val;
  }
  
  Scale_freqs(bfreqs,   4);
  Scale_freqs(bfreqs+4, 4);
  Scale_freqs(bfreqs+8, 4);
  
}
/*********************************************************/


void Randomize_Sequence_Order(calign *cdata)
{
  int i,exchange_with;
  phydbl buff_dbl;
  char *buff_name,*buff_state;
  short int *buff_ambigu;
  
  exchange_with = -1;
  For(i,cdata->n_otu)
  {
    buff_dbl  = rand();
    buff_dbl /= (RAND_MAX+1.);
    buff_dbl *= cdata->n_otu;
    exchange_with = (int)FLOOR(buff_dbl);
    
    buff_name                         = cdata->c_seq[i]->name;
    cdata->c_seq[i]->name             = cdata->c_seq[exchange_with]->name;
    cdata->c_seq[exchange_with]->name = buff_name;
    
    buff_state                         = cdata->c_seq[i]->state;
    cdata->c_seq[i]->state             = cdata->c_seq[exchange_with]->state;
    cdata->c_seq[exchange_with]->state = buff_state;
    
    buff_ambigu                            = cdata->c_seq[i]->is_ambigu;
    cdata->c_seq[i]->is_ambigu             = cdata->c_seq[exchange_with]->is_ambigu;
    cdata->c_seq[exchange_with]->is_ambigu = buff_ambigu;
  }
}

/*********************************************************/

//Modified by Ken 3/8/2016
void Update_Ancestors(t_node *a, t_node *d, t_tree *tree)
{
  d->anc = a;


  if(d->tax) return;
  else
  {
    int i;
    For(i,3)
    if((d->v[i] != d->anc) && (d->b[i] != tree->e_root))
      Update_Ancestors(d,d->v[i],tree);
  }
}

//Modified by Ken 3/8/2016
void Update_Ancestors_Edge(t_node *a, t_node *d, t_edge *e, t_tree *tree)
{
  d->anc = a;
  d->anc_edge=e;
  a->lupdate=0;//intialize lupdate=0 on all nodes
  d->lupdate=0;
  e->anc_node = a;
  e->des_node = d;

 // printf("a:%d\t%d\td:%d\t%d\t%d\t%lf\n",a->num,a->lupdate,d->num,d->lupdate,d->anc_edge->num,d->anc_edge->l);
 //printf("a:%d\td:%d\t%d\n",a->num,d->num,d->anc_edge->num);

  if(d->tax) return;
  else
  {
    int i;
    /*For(i,3){
    	printf("checking edges %d\n",d->b[i]->num);
        if(d->v[i]->num != d->anc->num){
        	d->b[i]->anc_edge=e;
        }
    }*/

    For(i,3){
	//	printf("%d\n",i);
	//	printf("%d\n",d->v[i]->num);
    	if((d->v[i] != d->anc)){
    		Update_Ancestors_Edge(d,d->v[i],d->b[i],tree);
    	}
    }
  }
}
/*********************************************************/
//Print ambiguous character matrixs
//Added by Ken 3/8/2016
void Print_Ambig_States(t_node *d, t_tree *tree, FILE *ambigfile)
{
	//if(d->tax && d->num != tree->mod->startnode){
	if(d->tax){
		//see if any character states are ambiguous, and then print off their partial lk matrices
		int site,j;
		For(site,tree->n_pattern){
			if(tree->data->c_seq[d->num]->is_ambigu[site]){
				if(d->num == tree->mod->startnode){
					printf("\n. Ambiguous character at root site %d",site);
				}
				For(j,tree->mod->ns){
					  phydbl temp;
					  if(d->num == tree->mod->startnode){
						  temp = (phydbl) d->b[0]->p_lk_tip_r[site*tree->mod->ns+j];
					  }else{
						  temp = (phydbl) d->anc_edge->p_lk_tip_r[site*tree->mod->ns+j];
					  }
					fprintf(ambigfile,"%s %d %d %lf\n",d->name,site,j,temp);
				}
			}
		}
	}
	if(!d->tax || d->num == tree->mod->startnode){
	    int i;
	    For(i,3)
	    if((d->v[i] != d->anc))
	    Print_Ambig_States(d->v[i],tree,ambigfile);
	  }
}

/*********************************************************/

void Best_Of_NNI_And_SPR(t_tree *tree)
{
  if(tree->mod->s_opt->random_input_tree) Speed_Spr_Loop(tree); /* Don't do simultaneous NNIs if starting tree is random */
  else
  {
    t_tree *ori_tree,*best_tree;
    model *ori_mod,*best_mod;
    phydbl *ori_bl,*best_bl;
    phydbl best_lnL,ori_lnL,nni_lnL,spr_lnL;
    
    ori_bl = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));
    best_bl = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));
    
    ori_mod   = Copy_Model(tree->mod);
    best_mod  = Copy_Model(tree->mod);
    
    ori_tree = Make_Tree(tree->n_otu);
    //Init_Tree(ori_tree,ori_tree->n_otu);
    Make_All_Tree_Nodes(ori_tree);
    Make_All_Tree_Edges(ori_tree);
    
    best_tree = Make_Tree(tree->n_otu);
    //Init_Tree(best_tree,best_tree->n_otu);
    Make_All_Tree_Nodes(best_tree);
    Make_All_Tree_Edges(best_tree);
    
    Lk(tree);
    Copy_Tree(tree,ori_tree);
    Record_Br_Len(ori_bl,tree);
    ori_lnL = tree->c_lnL; /* Record likelihood of the starting tree */
    
    Simu_Loop(tree); /* Perform simultaneous NNIs */
    best_lnL = tree->c_lnL; /* Record the likelihood */
    nni_lnL = tree->c_lnL;
    Copy_Tree(tree,best_tree); /* Record the tree topology and branch lengths */
    Record_Br_Len(best_bl,tree);
    Restore_Br_Len(best_bl,best_tree);
    Record_Model(tree->mod,best_mod);
    
    Copy_Tree(ori_tree,tree); /* Back to the original tree topology */
    Restore_Br_Len(ori_bl,tree); /* Back to the original branch lengths */
    Record_Model(ori_mod,tree->mod); /* Back to the original model */
    
    /* Make sure the tree is in its original form */
    Lk(tree);
//     if(FABS(tree->c_lnL - ori_lnL) > tree->mod->s_opt->min_diff_lk_global)
//     {
//       PhyML_Printf("\n. ori_lnL = %f, c_lnL = %f",ori_lnL,tree->c_lnL);
//       PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
//       Warn_And_Exit("");
//     }  //COMMENTED OUT BY MARCELO ... 24.07.2014
    
    Speed_Spr_Loop(tree);
    spr_lnL = tree->c_lnL;
    if(tree->c_lnL > best_lnL)
    {
      best_lnL = tree->c_lnL;
      Copy_Tree(tree,best_tree); /* Record tree topology, branch lengths and model parameters */
      Record_Br_Len(best_bl,tree);
      Restore_Br_Len(best_bl,best_tree);
      Record_Model(tree->mod,best_mod);
    }
    
    Copy_Tree(best_tree,tree);
    Restore_Br_Len(best_bl,tree);
    Record_Model(best_mod,tree->mod);
    
    /* Make sure the current tree has the best topology, branch lengths and model parameters */
    Lk(tree);
//     if(FABS(tree->c_lnL - best_lnL) > tree->mod->s_opt->min_diff_lk_global)
//     {
//       PhyML_Printf("\n. best_lnL = %f, c_lnL = %f",best_lnL,tree->c_lnL);
//       PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
//       Warn_And_Exit("");
//     } //COMMENTED OUT BY MARCELO ... 24.07.2014
    
    if(tree->mod->s_opt->print)
    {
      PhyML_Printf("\n\n. Log likelihood obtained after NNI moves : %f",nni_lnL);
      PhyML_Printf("\n. Log likelihood obtained after SPR moves : %f",spr_lnL);
    }
    
    free(ori_bl);
    free(best_bl);
    
    Free_Model(ori_mod);
    Free_Model(best_mod);
    
    Free_Tree(ori_tree);
    Free_Tree(best_tree);
  }
}

/*********************************************************/

t_tree *Dist_And_BioNJ(calign *cdata, model *mod, option *io)
{
  t_tree *tree;
  matrix *mat;
  // matrix *mat2;
  // int i,j;
  if(!io->quiet) PhyML_Printf("\n");
  if(!io->quiet) PhyML_Printf(". Computing pairwise distances...");
  mod->calculate_init_tree = 1; //!< Added by Marcelo. Used in the calculation of the probability transition matrix.
  if((io->datatype==CODON)&&(io->init_DistanceTreeCD!=NUCLEO)) //!< Added by Marcelo. 
  {
    mat = ML_CODONDist_Pairwise(cdata,mod->io);
    Fill_Missing_Dist(mat);
  }
  else
  {
    
    mat = ML_Dist(cdata,mod);
    Fill_Missing_Dist(mat);
    
  }
  mod->calculate_init_tree = 0; //!< Added by Marcelo. 
  if(!io->quiet) PhyML_Printf("\n. Building BioNJ tree...\n");
  mat->tree = Make_Tree_From_Scratch(cdata->n_otu,cdata);
  Bionj(mat);
  tree      = mat->tree;
  tree->mat = mat;
  return tree;
}

/*********************************************************/

void Add_BioNJ_Branch_Lengths(t_tree *tree, calign *cdata, model *mod)
{
  matrix *mat;
  
  PhyML_Printf("\n. Computing branch length estimates...\n");
  
  Order_Tree_CSeq(tree,cdata);
  
  mod->calculate_init_tree = 1; //!< Added by Marcelo. Used in the calculation of the probability transition matrix.
  if((mod->io->datatype==CODON)&&(mod->io->init_DistanceTreeCD!=NUCLEO)) //!< Added by Marcelo. 
  {
    mat = ML_CODONDist_Pairwise(cdata,mod->io);
    Fill_Missing_Dist(mat);
  }
  else
  {
    mat = ML_Dist(cdata,mod);
    Fill_Missing_Dist(mat);
  }
  mod->calculate_init_tree = 0; //!< Added by Marcelo.
  
  mat->tree = tree;
  mat->method = 0;
  Bionj_Br_Length(mat);
  
  Free_Mat(mat);
}

/*********************************************************/

t_tree *Read_User_Tree(calign *cdata, model *mod, option *io)
{
  t_tree *tree;
  
  
  PhyML_Printf("\n. Reading user tree...\n"); fflush(NULL);
  if(io->n_trees == 1) rewind(io->fp_in_tree);
  tree = Read_Tree_File_Phylip(io->fp_in_tree);
  if(!tree) Exit("\n. Input tree not found...\n");
  /* Add branch lengths if necessary */
  if(!tree->has_branch_lengths || io->nobl){ PhyML_Printf("\n. Tree does not have branch lengths or you chose to ignore them...\n"); fflush(NULL);Add_BioNJ_Branch_Lengths(tree,cdata,mod);}
  else PhyML_Printf("\n. The user tree does have branch lengths...\n");
  return tree;
}

/*********************************************************/
void Print_Time_Info(time_t t_beg, time_t t_end)
{
  div_t hour,min;
  
  hour = div((int)t_end-(int)t_beg,3600);
  min  = div((int)t_end-(int)t_beg,60  );
  min.quot -= hour.quot*60;
  
  PhyML_Printf("\n. Time used %dh%dm%ds\n", hour.quot,min.quot,(int)(t_end-t_beg)%60);
  PhyML_Printf("\noooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
}

/*********************************************************/

char *Bootstrap_From_String(char *s_tree, calign *cdata, model *mod, option *io)
{
  t_tree *tree;
  
  tree = Read_Tree(s_tree);
  
  if(!tree)
  {
    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    Warn_And_Exit("");
  }
  
  tree->mod         = mod;
  tree->io          = io;
  tree->data        = cdata;
  tree->both_sides  = 1;
  tree->n_pattern   = tree->data->crunch_len;
  
  Order_Tree_CSeq(tree,cdata);
  if(tree->mod->s_opt->random_input_tree) Random_Tree(tree);
  Fill_Dir_Table(tree);
  Update_Dirs(tree);
  Make_Tree_4_Pars(tree,cdata,cdata->init_len);
  Make_Tree_4_Lk(tree,cdata,cdata->init_len);
  tree->triplet_struct = Make_Triplet_Struct(mod);
  //Br_Len_Not_Involving_Invar(tree);
  Make_Spr_List(tree);
  Make_Best_Spr(tree);
  
#ifdef MPI
  Bootstrap_MPI(tree);
#else
  Bootstrap(tree);
#endif
  
  free(s_tree);  
  s_tree = Write_Tree(tree);
  
  Free_Spr_List(tree);
  Free_One_Spr(tree->best_spr);
  Free_Triplet(tree->triplet_struct);
  Free_Tree_Pars(tree);
  Free_Tree_Lk(tree);
  Free_Tree(tree);
  
  return s_tree;
}

/*********************************************************/

char *aLRT_From_String(char *s_tree, calign *cdata, model *mod, option *io)
{
  t_tree *tree;
  
  tree = Read_Tree(s_tree);
  
  if(!tree)
  {
    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    Warn_And_Exit("");
  }
  
  tree->mod         = mod;
  tree->io          = io;
  tree->data        = cdata;
  tree->both_sides  = 1;
  tree->n_pattern   = tree->data->crunch_len;
  
  Order_Tree_CSeq(tree,cdata);
  if(tree->mod->s_opt->random_input_tree) Random_Tree(tree);
  Fill_Dir_Table(tree);
  Update_Dirs(tree);
  Make_Tree_4_Pars(tree,cdata,cdata->init_len);
  Make_Tree_4_Lk(tree,cdata,cdata->init_len);
  tree->triplet_struct = Make_Triplet_Struct(mod);
  Make_Spr_List(tree);
  Make_Best_Spr(tree);
  
  tree->both_sides = 1;
  Lk(tree);
  
  aLRT(tree);
  
  free(s_tree);  
  s_tree = Write_Tree(tree);
  
  Free_Spr_List(tree);
  Free_One_Spr(tree->best_spr);
  Free_Triplet(tree->triplet_struct);
  Free_Tree_Pars(tree);
  Free_Tree_Lk(tree);
  Free_Tree(tree);
  
  return s_tree;
}

/*********************************************************/

void Prepare_Tree_For_Lk(t_tree *tree)
{

  Order_Tree_CSeq(tree,tree->data);
  Fill_Dir_Table(tree);
  Update_Dirs(tree);
  Make_Tree_4_Pars(tree,tree->data,tree->data->init_len);
  Make_Tree_4_Lk(tree,tree->data,tree->data->init_len);
  tree->triplet_struct = Make_Triplet_Struct(tree->mod);
  Br_Len_Not_Involving_Invar(tree);
  Make_Spr_List(tree);
  Make_Best_Spr(tree);
}

/*********************************************************/

void PhyML_Printf(char *format, ...)
{
  va_list ptr;
  
#ifdef MPI
  if(Global_myRank == 0)
  {
    va_start (ptr, format);
    vprintf (format, ptr);
    va_end(ptr);
  }
#else
  va_start (ptr, format);
  vprintf (format, ptr);
  va_end(ptr);
#endif
  
  fflush (NULL);
}

/*********************************************************/

void PhyML_Fprintf(FILE *fp, char *format, ...)
{
  va_list ptr;
  
#ifdef MPI
  if(Global_myRank == 0)
  {
    va_start (ptr, format);
    vfprintf (fp,format, ptr);
    va_end(ptr);
  }
#else
  va_start (ptr, format);
  vfprintf (fp,format, ptr);
  va_end(ptr);
#endif
  
  fflush (NULL);
}

/*********************************************************/

phydbl Get_Tree_Size(t_tree *tree)
{
  int i;
  phydbl tree_size;
  
  tree_size = 0.0;
  For(i,2*tree->n_otu-3) tree_size += tree->t_edges[i]->l;
  /*   tree_size = 0.0; */
  /*   For(i,2*tree->n_otu-3) tree_size += tree->rates->u_cur_l[i]; */
  
  /*   For(i,2*tree->n_otu-3)  */
  /*     tree_size +=  */
  /*     FABS(tree->rates->nd_t[tree->t_edges[i]->left->num] -  */
  /* 	 tree->rates->nd_t[tree->t_edges[i]->rght->num]); */
  
  tree->size = tree_size;
  return tree_size;
}

/*********************************************************/
/*********************************************************/
/*********************************************************/
/* 'Borrowed' fromn libgen */
char *Basename(char *path)
{
  char *p;
  
  if( path == NULL || *path == '\0' ) return ".";
  
  p = path + strlen(path) - 1;
  
  while( *p == '/' ) 
  {
    if( p == path ) return path;
    *p-- = '\0';
  }
  
  while( p >= path && *p != '/' ) p--;
  
  return p + 1;
}

/*********************************************************/

/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/

void Copy_Tree_Topology_With_Labels(t_tree *ori, t_tree *cpy)
{
  int i,j;
  
  For(i,2*ori->n_otu-2)
  {
    For(j,3)
    {
      if(ori->noeud[i]->v[j])
      {
        cpy->noeud[i]->v[j] = cpy->noeud[ori->noeud[i]->v[j]->num];
        cpy->noeud[i]->l[j] = ori->noeud[i]->l[j];
      }
      else
        cpy->noeud[i]->v[j] = NULL;
    }
    cpy->noeud[i]->num = ori->noeud[i]->num;
    cpy->noeud[i]->tax = 0;
  }
  
  For(i,2*ori->n_otu-3)
  {
    cpy->t_edges[i]->l = ori->t_edges[i]->l;
  }
  
  For(i,ori->n_otu)
  {
    cpy->noeud[i]->tax = 1;
    strcpy(cpy->noeud[i]->name,ori->noeud[i]->name);
  }
  
}

/*********************************************************/
//! Read the simulation options from interface or command line.

option *Get_Input(int argc, char **argv)
{
  
  option *io; 
  model *mod;
  optimiz *s_opt;
  
  io    = (option *)Make_Input();      //!< Allocate memory for pointers of strings that contain the names of input and output files.
  mod   = (model *)Make_Model_Basic(); //!< Read the options from interface or command line.
  s_opt = (optimiz *)Make_Optimiz();   //!< Allocate memory for the optimization options.
  
  Set_Defaults_Input(io);              //!< Initialize options with default inputs.             
  Set_Defaults_Model(mod); 
  Set_Defaults_Optimiz(s_opt);
  
  io->mod    = mod;
  mod->io    = io;
  mod->s_opt = s_opt;
  
  
#ifdef MPI
  Read_Command_Line(io,argc,argv);
#else
  
  putchar('\n');
  
  switch (argc)
  {
    case 1:
    {
    	printf("\n. Interface isn't implemented in IgPhyML yet. Please use command line."); //Added by Ken 4/1/2017
    	exit(EXIT_FAILURE);
//      Launch_Interface(io);//!< Read user defined options from interface.
      break;
    }
      /*
       case 2:
       Usage();
       break;
       */
    default:
      Read_Command_Line(io,argc,argv); //!< Read user defined options from command line.
  }
#endif
  
  if(io->datatype==CODON) Make_Ts_and_Tv_Matrix(io); //!< Added by Marcelo. Used mainly in EMCK07 model (Kosiol 2007);     
  
  return io;
}

/*********************************************************/

void Set_Model_Name(model *mod)
{
  if(mod->io->datatype == NT)
  {
    switch(mod->io->modeltypeOpt)
    {
      case JC69:
      {
        strcpy(mod->modelname, "JC69");
        break;
      }
      case K80:
      {
        strcpy(mod->modelname, "K80");
        break;
      }
      case F81:
      {
        strcpy(mod->modelname, "F81");
        break;
      }
      case HKY85:
      {
        strcpy(mod->modelname, "HKY85");
        break;
      }
      case F84:
      {
        strcpy(mod->modelname, "F84");
        break;
      }
      case TN93:
      {
        strcpy(mod->modelname, "TN93");
        break;
      }
      case GTR:
      {
        strcpy(mod->modelname, "GTR");
        break;
      }
      case CUSTOM:
      {
        strcpy(mod->modelname, "Custom");
        break;
      }
      default:
      {
        PhyML_Printf("\n. Unknown model name.\n");
        PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
        Warn_And_Exit("");
        break;
      }
    }
  }
  else if(mod->io->datatype == AA)
  {
    switch(mod->io->modeltypeOpt)
    {
      case DAYHOFF:
      {
        strcpy(mod->modelname, "Dayhoff");
        break;
      }
      case JTT:
      {
        strcpy(mod->modelname, "JTT");
        break;
      }
      case MTREV:
      {
        strcpy(mod->modelname, "MtREV");
        break;
      }
      case LG:
      {
        strcpy(mod->modelname, "LG");
        break;
      }
      case WAG:
      {
        strcpy(mod->modelname, "WAG");
        break;
      }
      case DCMUT:
      {
        strcpy(mod->modelname, "DCMut");
        break;
      }
      case RTREV:
      {
        strcpy(mod->modelname, "RtREV");
        break;
      }
      case CPREV:
      {
        strcpy(mod->modelname, "CpREV");
        break;
      }
      case VT:
      {
        strcpy(mod->modelname, "VT");
        break;
      }
      case BLOSUM62:
      {
        strcpy(mod->modelname, "Blosum62");
        break;
      }
      case MTMAM:
      {
        strcpy(mod->modelname, "MtMam");
        break;
      }
      case MTART:
      {
        strcpy(mod->modelname, "MtArt");
        break;
      }
      case HIVW:
      {
        strcpy(mod->modelname, "HIVb");
        break;
      }
      case HIVB:
      {
        strcpy(mod->modelname, "HIVb");
        break;
      }
      case CUSTOMAA:
      {
        strcpy(mod->modelname, "Custom");
        break;
      }
      default:
      {
        PhyML_Printf("\n. Unknown model name.\n");
        PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
        Warn_And_Exit("");
        break;
      }
    }
    if(mod->s_opt->opt_state_freq) strcat(mod->modelname,"+F"); //!Added by Marcelo.
  }
  else if(mod->io->datatype == GENERIC) {
    strcpy(mod->modelname, "JC69");
  }
  else if(mod->io->datatype == CODON) {                                         //!< Added by Marcelo.
    char * tempName = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    
    switch( mod->io->modeltypeOpt ) {
      case GY:
        strcpy( tempName, "GY" );
        break;
      case MG:
        strcpy( tempName, "MG" );
        break;
      case YAP:
        strcpy( tempName, "YAP" );
        break;
      case PCM:
        strcpy( tempName, "PCM" );
        break;
      case HLP17:
       strcpy( tempName, "HLP17" );
       break;
      default:
        break;
    }
    
    if((mod->io->modeltypeOpt != PCM) &&
        (mod->codon_model_nature != PARA)) {
      switch( mod->initqrates ) {
        case KOSI07:
          strcat( tempName, "ECMK07" );
          break;
        
        case SCHN05:
          strcat( tempName, "ECMS05" );
          break;
      
        case ECMUSR:
          strcat( tempName, "ECMUSR" );
          break;
        
        default:
          break;
      }
    
      if(mod->io->omegaOpt != NOOMEGA && 
         (mod->codon_model_nature != PARA &&
          mod->codon_model_nature != EMPI)) {
            strcat( tempName, "+W+K" );
      }
    }
    
    if(mod->io->freqmodelOpt != FUNDEFINED && 
        mod->io->freqmodelOpt != FMODEL &&
        mod->codon_model_nature != PARA &&
        mod->io->eq_freq_handling != MODEL) {
      strcat( tempName, "+F" );
    }
    
    strcpy( mod->modelname, tempName );
    
    free( tempName );
  }
}

/*********************************************************/
void Adjust_Min_Diff_Lk(t_tree *tree)
{
  int exponent;
  
  exponent = (int)FLOOR(log10(FABS(tree->c_lnL)));
  
  if(sizeof(phydbl) == 4)
  {
    tree->mod->s_opt->min_diff_lk_global = POW(10.,exponent - FLT_DIG + 1);
    tree->mod->s_opt->min_diff_lk_local  = tree->mod->s_opt->min_diff_lk_global;
    tree->mod->s_opt->min_diff_lk_move   = tree->mod->s_opt->min_diff_lk_global;
  }
  
  //   PhyML_Printf("\n. Exponent = %d Precision = %E DIG = %d",exponent,tree->mod->s_opt->min_diff_lk_global,FLT_DIG);
}


/*********************************************************/

void Skip_Comment(FILE *fp)
{
  int in_comment;
  char c;
  
  in_comment = 1;
  do
  {
    c = fgetc(fp);
    if(c == EOF) break;
    if(c == '[')      in_comment++;
    else if(c == ']') in_comment--;
  }
  while(in_comment);
}
/*********************************************************/

char     aminoAcidmap[65];                        //!< Added by Marcelo.
int        stopCodons[64];                        //!< Added by Marcelo.
int       senseCodons[64];                        //!< Added by Marcelo.
int  indexSenseCodons[64];                        //!< Added by Marcelo.
int            ThegenCode;                        //!< Added by Marcelo.

int Get_Genetic_Code()
{
  return ThegenCode;
}

void Set_Genetic_Code(int gencode) //! Added by Marcelo
{
  ThegenCode=gencode;
}

void Genetic_code_index_stopCodons(int genCode)   //!< Initialize //!< Added by Marcelo.
{
  int i,j,numsenseCodons;
  Set_Genetic_Code(genCode);
  numsenseCodons=Genetic_code_ns();
  For (i,64) 
  {
    stopCodons[i]=0; 
    senseCodons[i]=-1; 
    indexSenseCodons[i]=-1;
  }
  
  switch(genCode) /*!< standard genCode = 0. */
  {
    case STANDARD: //!Standard
    {
      stopCodons[10]=1; //!< Set stop codons according to the genetic code.
      stopCodons[11]=1; 
      stopCodons[14]=1; 
      strcpy(aminoAcidmap,"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\0"); //!< Amino acids according to the canonical order, below.
      break; 
    }
    case TVMC: //!Vertebrate Mithocondrial Code
    {
      stopCodons[10]=1;
      stopCodons[11]=1;
      stopCodons[46]=1;
      stopCodons[47]=1;
      strcpy(aminoAcidmap,"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG\0");
      break;
    }
    case TYMC: //!Yeast Mitochondrial Code
    {
      stopCodons[10]=1;
      stopCodons[11]=1;
      strcpy(aminoAcidmap,"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG\0");
      break;
    }
    case THMPCMCMSC: //!The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code 
    {
      stopCodons[10]=1;
      stopCodons[11]=1;
      strcpy(aminoAcidmap,"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\0");
      break;
    }
    case THIMC: //!Invertebrate Mitochondrial Code 
    {
      stopCodons[10]=1;
      stopCodons[11]=1;
      strcpy(aminoAcidmap,"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG\0");
      break;
    }
    case THCDHNC: //!The Ciliate, Dasycladacean and Hexamita Nuclear Code 
    {
      stopCodons[14]=1;
      strcpy(aminoAcidmap,"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\0");
      break;
    }
    case THEFMC: //!The Echinoderm and Flatworm Mitochondrial Code 
    {
      stopCodons[10]=1;
      stopCodons[11]=1;
      strcpy(aminoAcidmap,"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG\0");
      break;
    }
    case THENC: //!The Euplotid Nuclear Code
    {
      stopCodons[10]=1;
      stopCodons[11]=1;
      strcpy(aminoAcidmap,"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\0");
      break;
    }
    case THBAPPC: //!The Bacterial, Archaeal and Plant Plastid Code
    {
      stopCodons[10]=1;
      stopCodons[11]=1;
      stopCodons[14]=1;
      strcpy(aminoAcidmap,"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\0");
      break;
    }
    case THAYNC: //!The Alternative Yeast Nuclear Code
    {
      stopCodons[10]=1;
      stopCodons[11]=1;
      stopCodons[14]=1;
      strcpy(aminoAcidmap,"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\0");
      break;
    }    
    case THAMC: //!The Ascidian Mitochondrial Code
    {
      stopCodons[10]=1;
      stopCodons[11]=1;
      strcpy(aminoAcidmap,"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG\0");
      break;
    }     
    case THAFMC: //!The Alternative Flatworm Mitochondrial Code
    {
      stopCodons[11]=1;
      strcpy(aminoAcidmap,"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG\0");
      break;
    }      
    case BLNC: //!Blepharisma Nuclear Code 
    {
      stopCodons[10]=1;
      stopCodons[14]=1;
      strcpy(aminoAcidmap,"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\0");
      break;
    }      
    case CHMC: //!Chlorophycean Mitochondrial Code
    {
      stopCodons[10]=1;
      stopCodons[14]=1;
      strcpy(aminoAcidmap,"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\0");
      break;
    }    
    case TRMC: //!Trematode Mitochondrial Code 
    {
      stopCodons[10]=1;
      stopCodons[11]=1;
      strcpy(aminoAcidmap,"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG\0");
      break;
    }     
    case SCOMC: //!Scenedesmus obliquus mitochondrial Code
    {
      stopCodons[6] =1;
      stopCodons[10]=1;
      stopCodons[14]=1;
      strcpy(aminoAcidmap,"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\0");
      break;
    }    
    case THMC: //!Thraustochytrium Mitochondrial Code  
    {
      stopCodons[2] =1;
      stopCodons[10]=1;
      stopCodons[11]=1;
      stopCodons[14]=1;
      strcpy(aminoAcidmap,"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\0");
      break;
    }    
    default: Warn_And_Exit("Not implemented Genetic Code\n");break;   
  }
  
  j=0; 
  For(i,64) //!< Shift the sense codons down the array.
  {
    if(!stopCodons[i]) 
    {
      senseCodons[j++]=i;
    }
  }
  For(i,numsenseCodons) //!< Correct the index for recovery.
  {
    indexSenseCodons[senseCodons[i]]=i;
  }
  
}

//! codon identifier 
//! <"01234567891011121314151617181920212223242526272829303132343536373839404142434445464748495051525354555657585960616263". 
/*! TTT - 0, TTC - 1, TTA - 3, TTG - 4, and so on. */
/*! this makes easier the reference to stop codons. */
/*! reference ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt. */

int Genetic_code_ns() //!< Added by Marcelo.
{      
  int number = 0;
  int genCode = Get_Genetic_Code();
  switch(genCode) /*!< Give the number of stop codons in the current genetic code.*/
  { 
    case   STANDARD: number= 3; break; //! see above for description.
    case       TVMC: number= 4; break; 
    case       TYMC: number= 2; break;
    case THMPCMCMSC: number= 2; break; 
    case      THIMC: number= 2; break; 
    case    THCDHNC: number= 1; break;
    case     THEFMC: number= 2; break;
    case      THENC: number= 2; break;
    case    THBAPPC: number= 3; break;
    case     THAYNC: number= 3; break;
    case      THAMC: number= 2; break;
    case     THAFMC: number= 1; break;
    case       BLNC: number= 2; break;
    case       CHMC: number= 2; break;
    case       TRMC: number= 2; break;
    case      SCOMC: number= 2; break;
    case       THMC: number= 4; break;
    default: Warn_And_Exit("\nGenetic code not implemented.\n");break;                                                  
  }
  return 64-number;
}

/*********************************************************/
void Get_Base_Freqs_CODONS_FaXb(calign *cdata, align **data, int freqModel, model *mod) //!< Added by Marcelo.
{
  char curr_state;
  int i,j,k,l,m, counter;
  phydbl A,C,G,T;
  phydbl fA,fC,fG,fT;
  phydbl A1,C1,G1,T1, A2,C2,G2,T2, A3,C3,G3,T3;
  phydbl fA1,fC1,fG1,fT1, fA2,fC2,fG2,fT2, fA3,fC3,fG3,fT3;
  phydbl *freq,*codons,sum, sumf;
  
  counter=50;
  switch(freqModel)
  {
    case F1XSENSECODONS://!< F1x61
    {
      freq=(phydbl *)mCalloc(mod->ns,sizeof(phydbl));
      codons=(phydbl *)mCalloc(mod->ns,sizeof(phydbl));
      For(i,mod->ns) freq[i]=1.0/mod->ns; 
      
      For(k,1)
      {
        
        For(i,mod->ns) codons[i]=0.0; 
        
        For(i,cdata->n_otu)
        {
          For(j,data[0]->len)
          {
            curr_state=data[i]->state[j];
            if(curr_state>=(char)0 && curr_state<(char)64 )
            {
              codons[ indexSenseCodons[(int)curr_state] ]++;
            }
            else 
            {
              if(curr_state==(char) 88)
              {
                l=-1;
                while(data[i]->alternativeCodons[j][++l]<64);
                sumf=0;
                For(m,l) sumf+=freq[indexSenseCodons[data[i]->alternativeCodons[j][m]]];
                For(m,l) codons[indexSenseCodons[data[i]->alternativeCodons[j][m]]]+=freq[indexSenseCodons[data[i]->alternativeCodons[j][m]]]/sumf;
              }
            }
          }
        }
        sum=0.0;
        For(i,mod->ns) sum+=codons[i];
        For(i,mod->ns) freq[i]=codons[i]/sum;
      }
      
      For(i,mod->ns) cdata->b_frq[i]=freq[i];
      For(i,mod->num_base_freq) mod->base_freq[i]=cdata->b_frq[i];
      
      free(freq);
      free(codons);
      break;
    }
    case F1X4://!< F1x4
    {
      fT = fC = fA = fG = 0.25;	
      
      For(k,counter)
      {
        A = C = G = T = .0;
        For(i,cdata->n_otu)
        {
          For(j,data[i]->ntLen)
          {
            switch(data[i]->ntStates[j])
            {
              case 'A' : A++;
                break;
              case 'C' : C++;
                break;
              case 'G' : G++;
                break;
              case 'T' : T++;
                break;
              case 'U' : T++;
                break;
              case 'M' : C+=fC/(fC+fA); A+=fA/(fA+fC);
                break;
              case 'R' : G+=fG/(fA+fG); A+=fA/(fA+fG);
                break;
              case 'W' : T+=fT/(fA+fT); A+=fA/(fA+fT);
                break;
              case 'S' : C+=fC/(fC+fG); G+=fG/(fC+fG);
                break;
              case 'Y' : C+=fC/(fC+fT); T+=fT/(fT+fC);
                break;
              case 'K' : G+=fG/(fG+fT); T+=fT/(fT+fG);
                break;
              case 'B' : C+=fC/(fC+fG+fT); G+=fG/(fC+fG+fT); T+=fT/(fC+fG+fT);
                break;
              case 'D' : A+=fA/(fA+fG+fT); G+=fG/(fA+fG+fT); T+=fT/(fA+fG+fT);
                break;
              case 'H' : A+=fA/(fA+fC+fT); C+=fC/(fA+fC+fT); T+=fT/(fA+fC+fT);
                break;
              case 'V' : A+=fA/(fA+fC+fG); C+=fC/(fA+fC+fG); G+=fG/(fA+fC+fG);
                break;
              case 'N' : case 'X' : case '?' : case 'O' : case '-' :
                A+=fA; C+=fC; G+=fG; T+=fT; break;
              default : break;
            }
          }
        }
        fA = A/(A+C+G+T);
        fC = C/(A+C+G+T);
        fG = G/(A+C+G+T);
        fT = T/(A+C+G+T);
      }
      
      cdata->b_frq[0] = fT;
      cdata->b_frq[1] = fC;
      cdata->b_frq[2] = fA;
      cdata->b_frq[3] = fG;
      
      For(i,mod->num_base_freq) mod->base_freq[i]=cdata->b_frq[i];
      
      break;
    }
    case F3X4://!< F3x4 and CF3x4.
    case CF3X4:
    {
      fA1=fC1=fG1=fT1=fA2=fC2=fG2=fT2=fA3=fC3=fG3=fT3=0.25;
      
      For(k,counter)
      {
        A1=C1=G1=T1=A2=C2=G2=T2=A3=C3=G3=T3=0.0;
        For(i,cdata->n_otu)
        {
          Fors(j,data[i]->ntLen,3)
          {
            switch(data[i]->ntStates[j])
            {
              case 'A' : A1++;
                break;
              case 'C' : C1++;
                break;
              case 'G' : G1++;
                break;
              case 'T' : T1++;
                break;
              case 'U' : T1++;
                break;
              case 'M' : C1+=fC1/(fC1+fA1); A1+=fA1/(fA1+fC1);
                break;
              case 'R' : G1+=fG1/(fA1+fG1); A1+=fA1/(fA1+fG1);
                break;
              case 'W' : T1+=fT1/(fA1+fT1); A1+=fA1/(fA1+fT1);
                break;
              case 'S' : C1+=fC1/(fC1+fG1); G1+=fG1/(fC1+fG1);
                break;
              case 'Y' : C1+=fC1/(fC1+fT1); T1+=fT1/(fT1+fC1);
                break;
              case 'K' : G1+=fG1/(fG1+fT1); T1+=fT1/(fT1+fG1);
                break;
              case 'B' : C1+=fC1/(fC1+fG1+fT1); G1+=fG1/(fC1+fG1+fT1); T1+=fT1/(fC1+fG1+fT1);
                break;
              case 'D' : A1+=fA1/(fA1+fG1+fT1); G1+=fG1/(fA1+fG1+fT1); T1+=fT1/(fA1+fG1+fT1);
                break;
              case 'H' : A1+=fA1/(fA1+fC1+fT1); C1+=fC1/(fA1+fC1+fT1); T1+=fT1/(fA1+fC1+fT1);
                break;
              case 'V' : A1+=fA1/(fA1+fC1+fG1); C1+=fC1/(fA1+fC1+fG1); G1+=fG1/(fA1+fC1+fG1);
                break;
              case 'N' : case 'X' : case '?' : case 'O' : case '-' :
                A1+=fA1; C1+=fC1; G1+=fG1; T1+=fT1; break;
              default : break;
            }
            switch(data[i]->ntStates[j+1])
            {
              case 'A' : A2++;
                break;
              case 'C' : C2++;
                break;
              case 'G' : G2++;
                break;
              case 'T' : T2++;
                break;
              case 'U' : T2++;
                break;
              case 'M' : C2+=fC2/(fC2+fA2); A2+=fA2/(fA2+fC2);
                break;
              case 'R' : G2+=fG2/(fA2+fG2); A2+=fA2/(fA2+fG2);
                break;
              case 'W' : T2+=fT2/(fA2+fT2); A2+=fA2/(fA2+fT2);
                break;
              case 'S' : C2+=fC2/(fC2+fG2); G2+=fG2/(fC2+fG2);
                break;
              case 'Y' : C2+=fC2/(fC2+fT2); T2+=fT2/(fT2+fC2);
                break;
              case 'K' : G2+=fG2/(fG2+fT2); T2+=fT2/(fT2+fG2);
                break;
              case 'B' : C2+=fC2/(fC2+fG2+fT2); G2+=fG2/(fC2+fG2+fT2); T2+=fT2/(fC2+fG2+fT2);
                break;
              case 'D' : A2+=fA2/(fA2+fG2+fT2); G2+=fG2/(fA2+fG2+fT2); T2+=fT2/(fA2+fG2+fT2);
                break;
              case 'H' : A2+=fA2/(fA2+fC2+fT2); C2+=fC2/(fA2+fC2+fT2); T2+=fT2/(fA2+fC2+fT2);
                break;
              case 'V' : A2+=fA2/(fA2+fC2+fG2); C2+=fC2/(fA2+fC2+fG2); G2+=fG2/(fA2+fC2+fG2);
                break;
              case 'N' : case 'X' : case '?' : case 'O' : case '-' :
                A2+=fA2; C2+=fC2; G2+=fG2; T2+=fT2; break;
              default : break;
            }
            switch(data[i]->ntStates[j+2])
            {
              case 'A' : A3++;
                break;
              case 'C' : C3++;
                break;
              case 'G' : G3++;
                break;
              case 'T' : T3++;
                break;
              case 'U' : T3++;
                break;
              case 'M' : C3+=fC3/(fC3+fA3); A3+=fA3/(fA3+fC3);
                break;
              case 'R' : G3+=fG3/(fA3+fG3); A3+=fA3/(fA3+fG3);
                break;
              case 'W' : T3+=fT3/(fA3+fT3); A3+=fA3/(fA3+fT3);
                break;
              case 'S' : C3+=fC3/(fC3+fG3); G3+=fG3/(fC3+fG3);
                break;
              case 'Y' : C3+=fC3/(fC3+fT3); T3+=fT3/(fT3+fC3);
                break;
              case 'K' : G3+=fG3/(fG3+fT3); T3+=fT3/(fT3+fG3);
                break;
              case 'B' : C3+=fC3/(fC3+fG3+fT3); G3+=fG3/(fC3+fG3+fT3); T3+=fT3/(fC3+fG3+fT3);
                break;
              case 'D' : A3+=fA3/(fA3+fG3+fT3); G3+=fG3/(fA3+fG3+fT3); T3+=fT3/(fA3+fG3+fT3);
                break;
              case 'H' : A3+=fA3/(fA3+fC3+fT3); C3+=fC3/(fA3+fC3+fT3); T3+=fT3/(fA3+fC3+fT3);
                break;
              case 'V' : A3+=fA3/(fA3+fC3+fG3); C3+=fC3/(fA3+fC3+fG3); G3+=fG3/(fA3+fC3+fG3);
                break;
              case 'N' : case 'X' : case '?' : case 'O' : case '-' :
                A3+=fA3; C3+=fC3; G3+=fG3; T3+=fT3; break;
              default : break;
            }
          }
        }
        fA1 = A1/(A1+C1+G1+T1);
        fC1 = C1/(A1+C1+G1+T1);
        fG1 = G1/(A1+C1+G1+T1);
        fT1 = T1/(A1+C1+G1+T1);
        fA2 = A2/(A2+C2+G2+T2);
        fC2 = C2/(A2+C2+G2+T2);
        fG2 = G2/(A2+C2+G2+T2);
        fT2 = T2/(A2+C2+G2+T2);
        fA3 = A3/(A3+C3+G3+T3);
        fC3 = C3/(A3+C3+G3+T3);
        fG3 = G3/(A3+C3+G3+T3);
        fT3 = T3/(A3+C3+G3+T3);
      }
      
      cdata->b_frq[0 ] = fT1;
      cdata->b_frq[1 ] = fC1;
      cdata->b_frq[2 ] = fA1;
      cdata->b_frq[3 ] = fG1;
      cdata->b_frq[4 ] = fT2;
      cdata->b_frq[5 ] = fC2;
      cdata->b_frq[6] = fA2;
      cdata->b_frq[7] = fG2;
      cdata->b_frq[8] = fT3;
      cdata->b_frq[9] = fC3;
      cdata->b_frq[10] = fA3;
      cdata->b_frq[11] = fG3;
      
      if(freqModel==CF3X4) CF3x4(cdata->b_frq, mod->genetic_code);
      
      For(i,mod->num_base_freq) mod->base_freq[i]=cdata->b_frq[i];
      
      break;
    }
    default: 
    {
      Warn_And_Exit("Model for empirical calculation of base frequencies not implemented. Impossible to continue.\n");
      break;
    }
  }
}
/*********************************************************/
void  CopyExtraFieldsCodon(align *from, int site, align *to,  int num_pattern)//!< Added by Marcelo.
{
  int i,j=-1;  
  to->is_ambigu[num_pattern]=from->is_ambigu[site];
  
  while(from->alternativeCodons[site][++j]<64);
  
  if(!(to->alternativeCodons[num_pattern])) to->alternativeCodons[num_pattern]=(char *)mCalloc(j+1,sizeof(char)); 
  For(i,j+1) to->alternativeCodons[num_pattern][i]=from->alternativeCodons[site][i]; 
}
/*********************************************************/
int Intersect_Site_Alternatives(align **data, int n_otu, int site)    /*!< Intersect ambiguities across sequences ... if empty the site is not invariable, if not empty*/
{                                                                     /*! the site can be considered invariable.*/ //!< Added by Marcelo.
  int i,j,k=-1;
  char set[64];
  
  For(i,64) set[i]=0;
  
  For(i,n_otu)                                                        /*!< Calculate intersection. Sum one to every possibility. At the end the elements that */
  {                                                                   /*!< where marked for all taxa, if any, form the nonempty intersection.*/  
    if(data[i]->is_ambigu[site]) 
    {
      k=-1; 
      while(data[i]->alternativeCodons[site][++k]<64);
      For(j,k) set[ data[i]->alternativeCodons[site][j] ]++;
    }else set[ data[i]->state[site] ]++;
  }
  
  For(i,64) if(set[i]==n_otu) return 1;                               /*!< Not empty intersectiong means site without polymorphism.*/
  
  return -1;                                                          /*!< Empty intersection polymorphism.*/
} 
/*********************************************************/            
void CF3x4(phydbl * freq, int genCode) //!< Added by Marcelo.
{ //!< Corrected F3x4 estimation for the equilibrium frequencies.
  int i, n, info, rhs;
  phydbl freqStopcodons , *b, *x, *fq, *J;
  
  n=9;
  rhs=1;
  
  b=(phydbl *)mCalloc(n, sizeof(phydbl));
  x=(phydbl *)mCalloc(n, sizeof(phydbl));
  fq=(phydbl *)mCalloc(n, sizeof(phydbl));
  J=(phydbl *)mCalloc(n*n, sizeof(phydbl)); //!< Jacobian.
  
  freqStopcodons=0.0;
  For(i,64) if(stopCodons[i]) freqStopcodons+=F3x4(i,freq);
  For(i,n) fq[i]=freq[cf3x4Ind[i]]*(1.0-freqStopcodons);
  For(i,n) x[i]=0.25;
  For(i,n) b[i]=1.0;
  switch(genCode)
  {
    case STANDARD: //!< too complicated to program ... therefore matrices are hardcoded ... add similar in the case of using another than standard genetic code. calculations according to kosakovsky et all (plos).
    {
      while(my2norm(b,n)>1e-10)
      {
        
        J[0*n+0]=-1+x[8]*(1-x[3]-x[4]-x[5])+(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[0*n+1]=0;
        J[0*n+2]=0;
        J[0*n+3]=-x[8]*x[0];
        J[0*n+4]=-x[8]*x[0];
        J[0*n+5]=(1-x[6]-x[7]-x[8])*x[0];
        J[0*n+6]=-x[5]*x[0];
        J[0*n+7]=-x[5]*x[0];
        J[0*n+8]=(1-x[3]-x[4]-x[5])*x[0];
        
        
        J[1*n+0]=0;
        J[1*n+1]=-1;
        J[1*n+2]=0;
        J[1*n+3]=0;
        J[1*n+4]=0;
        J[1*n+5]=0;
        J[1*n+6]=0;
        J[1*n+7]=0;
        J[1*n+8]=0;
        
        
        J[2*n+0]=0;
        J[2*n+1]=0;
        J[2*n+2]=-1;
        J[2*n+3]=0;
        J[2*n+4]=0;
        J[2*n+5]=0;
        J[2*n+6]=0;
        J[2*n+7]=0;
        J[2*n+8]=0;
        
        
        J[3*n+0]=0;
        J[3*n+1]=0;
        J[3*n+2]=0;
        J[3*n+3]=-1;
        J[3*n+4]=0;
        J[3*n+5]=0;
        J[3*n+6]=0;
        J[3*n+7]=0;
        J[3*n+8]=0;
        
        
        J[4*n+0]=0;
        J[4*n+1]=0;
        J[4*n+2]=0;
        J[4*n+3]=0;
        J[4*n+4]=-1;
        J[4*n+5]=0;
        J[4*n+6]=0;
        J[4*n+7]=0;
        J[4*n+8]=0;
        
        
        J[5*n+0]=(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[5*n+1]=0;
        J[5*n+2]=0;
        J[5*n+3]=0;
        J[5*n+4]=0;
        J[5*n+5]=-1+(1-x[6]-x[7]-x[8])*x[0]+x[8]*x[0];
        J[5*n+6]=-x[5]*x[0];
        J[5*n+7]=-x[5]*x[0];
        J[5*n+8]=0;
        
        
        J[6*n+0]=0;
        J[6*n+1]=0;
        J[6*n+2]=0;
        J[6*n+3]=0;
        J[6*n+4]=0;
        J[6*n+5]=0;
        J[6*n+6]=-1;
        J[6*n+7]=0;
        J[6*n+8]=0;
        
        
        J[7*n+0]=0;
        J[7*n+1]=0;
        J[7*n+2]=0;
        J[7*n+3]=0;
        J[7*n+4]=0;
        J[7*n+5]=0;
        J[7*n+6]=0;
        J[7*n+7]=-1;
        J[7*n+8]=0;
        
        
        J[8*n+0]=x[8]*(1-x[3]-x[4]-x[5])+x[8]*x[5];
        J[8*n+1]=0;
        J[8*n+2]=0;
        J[8*n+3]=-x[8]*x[0];
        J[8*n+4]=-x[8]*x[0];
        J[8*n+5]=0;
        J[8*n+6]=0;
        J[8*n+7]=0;
        J[8*n+8]=-1+(1-x[3]-x[4]-x[5])*x[0]+x[5]*x[0];
        
        
        b[0]=fq[0] - x[0]*(1 -x[5]*x[8] - x[5]*(1-(x[6]+x[7]+x[8])) - (1-(x[3]+x[4]+x[5]))*x[8] ); //!< f(x).
        b[1]=fq[1] - x[1];
        b[2]=fq[2] - x[2];
        
        b[3]=fq[3] - x[3];
        b[4]=fq[4] - x[4];
        b[5]=fq[5] - x[5]*(1- x[0]*x[8] - x[0]*(1-(x[6]+x[7]+x[8])));
        
        b[6]=fq[6] - x[6];
        b[7]=fq[7] - x[7];
        b[8]=fq[8] - x[8]*(1- x[0]*x[5] - x[0]*(1-(x[3]+x[4]+x[5])));
        
#if defined BLAS || defined BLAS_OMP
        
        int  *ipiv;
        ipiv=(int *)mCalloc(n, sizeof(int));
        
        dgesv_(&n, &rhs, J, &n, ipiv, b, &n, &info);
        
        free(ipiv);
        
#else
        
        info=MyGaussElimination_gNxRHS(J,b,n,rhs);
        
#endif
        
        if(!info) 
          For(i,n) x[i]-=b[i];
        else 
          Warn_And_Exit("Error in Gauss Elimination Scheme.\n");
        
      }
      
      break;
    }
    case TVMC: 
    {
      while(my2norm(b,n)>1e-10)
      {
        
        J[0*n+0]=-1+(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[0*n+1]=0;                                                                                                                                                                                               
        J[0*n+2]=0;
        J[0*n+3]=0;
        J[0*n+4]=0;
        J[0*n+5]=(1-x[6]-x[7]-x[8])*x[0]+x[8]*x[0];
        J[0*n+6]=-x[5]*x[0];                                                                                                                                                                                      
        J[0*n+7]=-x[5]*x[0];                                                                                                                                                                                      
        J[0*n+8]=0;
        
        
        J[1*n+0]=0;
        J[1*n+1]=-1;
        J[1*n+2]=0;
        J[1*n+3]=0;
        J[1*n+4]=0;
        J[1*n+5]=0;
        J[1*n+6]=0;
        J[1*n+7]=0;
        J[1*n+8]=0;
        
        
        J[2*n+0]=0;
        J[2*n+1]=0;
        J[2*n+2]=-1+(1-x[6]-x[7]-x[8])*(1-x[3]-x[4]-x[5])+x[8]*(1-x[3]-x[4]-x[5]);
        J[2*n+3]=-(1-x[6]-x[7]-x[8])*x[2]-x[8]*x[2];
        J[2*n+4]=-(1-x[6]-x[7]-x[8])*x[2]-x[8]*x[2];
        J[2*n+5]=-(1-x[6]-x[7]-x[8])*x[2]-x[8]*x[2];
        J[2*n+6]=-(1-x[3]-x[4]-x[5])*x[2];
        J[2*n+7]=-(1-x[3]-x[4]-x[5])*x[2];
        J[2*n+8]=0;
        
        
        J[3*n+0]=0;
        J[3*n+1]=0;
        J[3*n+2]=0;
        J[3*n+3]=-1;
        J[3*n+4]=0;
        J[3*n+5]=0;
        J[3*n+6]=0;
        J[3*n+7]=0;
        J[3*n+8]=0;
        
        
        J[4*n+0]=0;
        J[4*n+1]=0;
        J[4*n+2]=0;
        J[4*n+3]=0;
        J[4*n+4]=-1;
        J[4*n+5]=0;
        J[4*n+6]=0;
        J[4*n+7]=0;
        J[4*n+8]=0;
        
        
        J[5*n+0]=(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[5*n+1]=0;
        J[5*n+2]=0;
        J[5*n+3]=0;
        J[5*n+4]=0;
        J[5*n+5]=-1+(1-x[6]-x[7]-x[8])*x[0]+x[8]*x[0];
        J[5*n+6]=-x[5]*x[0];
        J[5*n+7]=-x[5]*x[0];
        J[5*n+8]=0;
        
        
        J[6*n+0]=0;
        J[6*n+1]=0;
        J[6*n+2]=0;
        J[6*n+3]=0;
        J[6*n+4]=0;
        J[6*n+5]=0;
        J[6*n+6]=-1;
        J[6*n+7]=0;
        J[6*n+8]=0;
        
        
        J[7*n+0]=0;
        J[7*n+1]=0;
        J[7*n+2]=0;
        J[7*n+3]=0;
        J[7*n+4]=0;
        J[7*n+5]=0;
        J[7*n+6]=0;
        J[7*n+7]=-1;
        J[7*n+8]=0;
        
        
        J[8*n+0]=x[8]*x[5];
        J[8*n+1]=0;
        J[8*n+2]=x[8]*(1-x[3]-x[4]-x[5]);
        J[8*n+3]=-x[8]*x[2];
        J[8*n+4]=-x[8]*x[2];
        J[8*n+5]=-x[8]*x[2]+x[8]*x[0];
        J[8*n+6]=0;
        J[8*n+7]=0;
        J[8*n+8]=-1+(1-x[3]-x[4]-x[5])*x[2]+x[5]*x[0];
        
        
        b[0] = fq[0] - x[0]+((1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+x[8]*x[5]*x[0]*1+0);
        b[1] = fq[1] - x[1];
        b[2] = fq[2] - x[2]+((1 - (x[6]+x[7]+x[8]))*(1 - (x[3]+x[4]+x[5]))*x[2]*1+x[8]*(1 - (x[3]+x[4]+x[5]))*x[2]*1+0);
        b[3] = fq[3] - x[3];
        b[4] = fq[4] - x[4];
        b[5] = fq[5] - x[5]+((1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+x[8]*x[5]*x[0]*1+0);
        b[6] = fq[6] - x[6];
        b[7] = fq[7] - x[7];
        b[8] = fq[8] - x[8]+(x[8]*(1 - (x[3]+x[4]+x[5]))*x[2]*1+x[8]*x[5]*x[0]*1+0);
        
#if defined BLAS || defined BLAS_OMP
        
        int  *ipiv;
        ipiv=(int *)mCalloc(n, sizeof(int));
        
        dgesv_(&n, &rhs, J, &n, ipiv, b, &n, &info);
        
        free(ipiv);
        
#else
        
        info=MyGaussElimination_gNxRHS(J,b,n,rhs);
        
#endif
        
        if(!info) 
          For(i,n) x[i]-=b[i];
        else 
          Warn_And_Exit("Error in Gauss Elimination Scheme.\n");
        
      }
      
      break;
    }
    case TYMC: 
    {
      while(my2norm(b,n)>1e-10)
      {
        
        J[0*n+0]=-1+(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[0*n+1]=0;
        J[0*n+2]=0;
        J[0*n+3]=0;
        J[0*n+4]=0;
        J[0*n+5]=(1-x[6]-x[7]-x[8])*x[0]+x[8]*x[0];
        J[0*n+6]=-x[5]*x[0];
        J[0*n+7]=-x[5]*x[0];
        J[0*n+8]=0;
        
        
        J[1*n+0]=0;
        J[1*n+1]=-1;
        J[1*n+2]=0;
        J[1*n+3]=0;
        J[1*n+4]=0;
        J[1*n+5]=0;
        J[1*n+6]=0;
        J[1*n+7]=0;
        J[1*n+8]=0;
        
        
        J[2*n+0]=0;
        J[2*n+1]=0;
        J[2*n+2]=-1;
        J[2*n+3]=0;
        J[2*n+4]=0;
        J[2*n+5]=0;
        J[2*n+6]=0;
        J[2*n+7]=0;
        J[2*n+8]=0;
        
        
        J[3*n+0]=0;
        J[3*n+1]=0;
        J[3*n+2]=0;
        J[3*n+3]=-1;
        J[3*n+4]=0;
        J[3*n+5]=0;
        J[3*n+6]=0;
        J[3*n+7]=0;
        J[3*n+8]=0;
        
        
        J[4*n+0]=0;
        J[4*n+1]=0;
        J[4*n+2]=0;
        J[4*n+3]=0;
        J[4*n+4]=-1;
        J[4*n+5]=0;
        J[4*n+6]=0;
        J[4*n+7]=0;
        J[4*n+8]=0;
        
        
        J[5*n+0]=(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[5*n+1]=0;
        J[5*n+2]=0;
        J[5*n+3]=0;
        J[5*n+4]=0;
        J[5*n+5]=-1+(1-x[6]-x[7]-x[8])*x[0]+x[8]*x[0];
        J[5*n+6]=-x[5]*x[0];
        J[5*n+7]=-x[5]*x[0];
        J[5*n+8]=0;
        
        
        J[6*n+0]=0;
        J[6*n+1]=0;
        J[6*n+2]=0;
        J[6*n+3]=0;
        J[6*n+4]=0;
        J[6*n+5]=0;
        J[6*n+6]=-1;
        J[6*n+7]=0;
        J[6*n+8]=0;
        
        
        J[7*n+0]=0;
        J[7*n+1]=0;
        J[7*n+2]=0;
        J[7*n+3]=0;
        J[7*n+4]=0;
        J[7*n+5]=0;
        J[7*n+6]=0;
        J[7*n+7]=-1;
        J[7*n+8]=0;
        
        
        J[8*n+0]=x[8]*x[5];
        J[8*n+1]=0;
        J[8*n+2]=0;
        J[8*n+3]=0;
        J[8*n+4]=0;
        J[8*n+5]=x[8]*x[0];
        J[8*n+6]=0;
        J[8*n+7]=0;
        J[8*n+8]=-1+x[5]*x[0];
        
        
        b[0] = fq[0] - x[0]+((1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+x[8]*x[5]*x[0]*1+0);
        b[1] = fq[1] - x[1];
        b[2] = fq[2] - x[2];
        b[3] = fq[3] - x[3];
        b[4] = fq[4] - x[4];
        b[5] = fq[5] - x[5]+((1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+x[8]*x[5]*x[0]*1+0);
        b[6] = fq[6] - x[6];
        b[7] = fq[7] - x[7];
        b[8] = fq[8] - x[8]+(x[8]*x[5]*x[0]*1+0);
        
#if defined BLAS || defined BLAS_OMP
        
        int  *ipiv;
        ipiv=(int *)mCalloc(n, sizeof(int));
        
        dgesv_(&n, &rhs, J, &n, ipiv, b, &n, &info);
        
        free(ipiv);
        
#else
        
        info=MyGaussElimination_gNxRHS(J,b,n,rhs);
        
#endif
        
        if(!info) 
          For(i,n) x[i]-=b[i];
        else 
          Warn_And_Exit("Error in Gauss Elimination Scheme.\n");
        
      }
      
      break;
    }
    case THMPCMCMSC: 
    {
      while(my2norm(b,n)>1e-10)
      {
        
        J[0*n+0]=-1+(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[0*n+1]=0;
        J[0*n+2]=0;
        J[0*n+3]=0;
        J[0*n+4]=0;
        J[0*n+5]=(1-x[6]-x[7]-x[8])*x[0]+x[8]*x[0];
        J[0*n+6]=-x[5]*x[0];
        J[0*n+7]=-x[5]*x[0];
        J[0*n+8]=0;
        
        
        J[1*n+0]=0;
        J[1*n+1]=-1;
        J[1*n+2]=0;
        J[1*n+3]=0;
        J[1*n+4]=0;
        J[1*n+5]=0;
        J[1*n+6]=0;
        J[1*n+7]=0;
        J[1*n+8]=0;
        
        
        J[2*n+0]=0;
        J[2*n+1]=0;
        J[2*n+2]=-1;
        J[2*n+3]=0;
        J[2*n+4]=0;
        J[2*n+5]=0;
        J[2*n+6]=0;
        J[2*n+7]=0;
        J[2*n+8]=0;
        
        
        J[3*n+0]=0;
        J[3*n+1]=0;
        J[3*n+2]=0;
        J[3*n+3]=-1;
        J[3*n+4]=0;
        J[3*n+5]=0;
        J[3*n+6]=0;
        J[3*n+7]=0;
        J[3*n+8]=0;
        
        
        J[4*n+0]=0;
        J[4*n+1]=0;
        J[4*n+2]=0;
        J[4*n+3]=0;
        J[4*n+4]=-1;
        J[4*n+5]=0;
        J[4*n+6]=0;
        J[4*n+7]=0;
        J[4*n+8]=0;
        
        
        J[5*n+0]=(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[5*n+1]=0;
        J[5*n+2]=0;
        J[5*n+3]=0;
        J[5*n+4]=0;
        J[5*n+5]=-1+(1-x[6]-x[7]-x[8])*x[0]+x[8]*x[0];
        J[5*n+6]=-x[5]*x[0];
        J[5*n+7]=-x[5]*x[0];
        J[5*n+8]=0;
        
        
        J[6*n+0]=0;
        J[6*n+1]=0;
        J[6*n+2]=0;
        J[6*n+3]=0;
        J[6*n+4]=0;
        J[6*n+5]=0;
        J[6*n+6]=-1;
        J[6*n+7]=0;
        J[6*n+8]=0;
        
        
        J[7*n+0]=0;
        J[7*n+1]=0;
        J[7*n+2]=0;
        J[7*n+3]=0;
        J[7*n+4]=0;
        J[7*n+5]=0;
        J[7*n+6]=0;
        J[7*n+7]=-1;
        J[7*n+8]=0;
        
        
        J[8*n+0]=x[8]*x[5];
        J[8*n+1]=0;
        J[8*n+2]=0;
        J[8*n+3]=0;
        J[8*n+4]=0;
        J[8*n+5]=x[8]*x[0];
        J[8*n+6]=0;
        J[8*n+7]=0;
        J[8*n+8]=-1+x[5]*x[0];
        
        
        b[0] = fq[0] - x[0]+((1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+x[8]*x[5]*x[0]*1+0);
        b[1] = fq[1] - x[1];
        b[2] = fq[2] - x[2];
        b[3] = fq[3] - x[3];
        b[4] = fq[4] - x[4];
        b[5] = fq[5] - x[5]+((1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+x[8]*x[5]*x[0]*1+0);
        b[6] = fq[6] - x[6];
        b[7] = fq[7] - x[7];
        b[8] = fq[8] - x[8]+(x[8]*x[5]*x[0]*1+0);
        
#if defined BLAS || defined BLAS_OMP
        
        int  *ipiv;
        ipiv=(int *)mCalloc(n, sizeof(int));
        
        dgesv_(&n, &rhs, J, &n, ipiv, b, &n, &info);
        
        free(ipiv);
        
#else
        
        info=MyGaussElimination_gNxRHS(J,b,n,rhs);
        
#endif
        
        if(!info) 
          For(i,n) x[i]-=b[i];
        else 
          Warn_And_Exit("Error in Gauss Elimination Scheme.\n");
        
      }
      
      break;
    }
    case THIMC: 
    {
      while(my2norm(b,n)>1e-10)
      {
        
        J[0*n+0]=-1+(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[0*n+1]=0;
        J[0*n+2]=0;
        J[0*n+3]=0;
        J[0*n+4]=0;
        J[0*n+5]=(1-x[6]-x[7]-x[8])*x[0]+x[8]*x[0];
        J[0*n+6]=-x[5]*x[0];
        J[0*n+7]=-x[5]*x[0];
        J[0*n+8]=0;
        
        
        J[1*n+0]=0;
        J[1*n+1]=-1;
        J[1*n+2]=0;
        J[1*n+3]=0;
        J[1*n+4]=0;
        J[1*n+5]=0;
        J[1*n+6]=0;
        J[1*n+7]=0;
        J[1*n+8]=0;
        
        
        J[2*n+0]=0;
        J[2*n+1]=0;
        J[2*n+2]=-1;
        J[2*n+3]=0;
        J[2*n+4]=0;
        J[2*n+5]=0;
        J[2*n+6]=0;
        J[2*n+7]=0;
        J[2*n+8]=0;
        
        
        J[3*n+0]=0;
        J[3*n+1]=0;
        J[3*n+2]=0;
        J[3*n+3]=-1;
        J[3*n+4]=0;
        J[3*n+5]=0;
        J[3*n+6]=0;
        J[3*n+7]=0;
        J[3*n+8]=0;
        
        
        J[4*n+0]=0;
        J[4*n+1]=0;
        J[4*n+2]=0;
        J[4*n+3]=0;
        J[4*n+4]=-1;
        J[4*n+5]=0;
        J[4*n+6]=0;
        J[4*n+7]=0;
        J[4*n+8]=0;
        
        
        J[5*n+0]=(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[5*n+1]=0;
        J[5*n+2]=0;
        J[5*n+3]=0;
        J[5*n+4]=0;
        J[5*n+5]=-1+(1-x[6]-x[7]-x[8])*x[0]+x[8]*x[0];
        J[5*n+6]=-x[5]*x[0];
        J[5*n+7]=-x[5]*x[0];
        J[5*n+8]=0;
        
        
        J[6*n+0]=0;
        J[6*n+1]=0;
        J[6*n+2]=0;
        J[6*n+3]=0;
        J[6*n+4]=0;
        J[6*n+5]=0;
        J[6*n+6]=-1;
        J[6*n+7]=0;
        J[6*n+8]=0;
        
        
        J[7*n+0]=0;
        J[7*n+1]=0;
        J[7*n+2]=0;
        J[7*n+3]=0;
        J[7*n+4]=0;
        J[7*n+5]=0;
        J[7*n+6]=0;
        J[7*n+7]=-1;
        J[7*n+8]=0;
        
        
        J[8*n+0]=x[8]*x[5];
        J[8*n+1]=0;
        J[8*n+2]=0;
        J[8*n+3]=0;
        J[8*n+4]=0;
        J[8*n+5]=x[8]*x[0];
        J[8*n+6]=0;
        J[8*n+7]=0;
        J[8*n+8]=-1+x[5]*x[0];
        
        
        b[0] = fq[0] - x[0]+((1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+x[8]*x[5]*x[0]*1+0);
        b[1] = fq[1] - x[1];
        b[2] = fq[2] - x[2];
        b[3] = fq[3] - x[3];
        b[4] = fq[4] - x[4];
        b[5] = fq[5] - x[5]+((1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+x[8]*x[5]*x[0]*1+0);
        b[6] = fq[6] - x[6];
        b[7] = fq[7] - x[7];
        b[8] = fq[8] - x[8]+(x[8]*x[5]*x[0]*1+0);
        
        
#if defined BLAS || defined BLAS_OMP
        
        int  *ipiv;
        ipiv=(int *)mCalloc(n, sizeof(int));
        
        dgesv_(&n, &rhs, J, &n, ipiv, b, &n, &info);
        
        free(ipiv);
        
#else
        
        info=MyGaussElimination_gNxRHS(J,b,n,rhs);
        
#endif
        
        if(!info) 
          For(i,n) x[i]-=b[i];
        else 
          Warn_And_Exit("Error in Gauss Elimination Scheme.\n");
        
      }
      
      break;
    }
    case THCDHNC: 
    {
      while(my2norm(b,n)>1e-10)
      {
        
        J[0*n+0]=-1+x[8]*(1-x[3]-x[4]-x[5]);
        J[0*n+1]=0;
        J[0*n+2]=0;
        J[0*n+3]=-x[8]*x[0];
        J[0*n+4]=-x[8]*x[0];
        J[0*n+5]=-x[8]*x[0];
        J[0*n+6]=0;
        J[0*n+7]=0;
        J[0*n+8]=(1-x[3]-x[4]-x[5])*x[0];
        
        
        J[1*n+0]=0;
        J[1*n+1]=-1;
        J[1*n+2]=0;
        J[1*n+3]=0;
        J[1*n+4]=0;
        J[1*n+5]=0;
        J[1*n+6]=0;
        J[1*n+7]=0;
        J[1*n+8]=0;
        
        
        J[2*n+0]=0;
        J[2*n+1]=0;
        J[2*n+2]=-1;
        J[2*n+3]=0;
        J[2*n+4]=0;
        J[2*n+5]=0;
        J[2*n+6]=0;
        J[2*n+7]=0;
        J[2*n+8]=0;
        
        
        J[3*n+0]=0;
        J[3*n+1]=0;
        J[3*n+2]=0;
        J[3*n+3]=-1;
        J[3*n+4]=0;
        J[3*n+5]=0;
        J[3*n+6]=0;
        J[3*n+7]=0;
        J[3*n+8]=0;
        
        
        J[4*n+0]=0;
        J[4*n+1]=0;
        J[4*n+2]=0;
        J[4*n+3]=0;
        J[4*n+4]=-1;
        J[4*n+5]=0;
        J[4*n+6]=0;
        J[4*n+7]=0;
        J[4*n+8]=0;
        
        
        J[5*n+0]=0;
        J[5*n+1]=0;
        J[5*n+2]=0;
        J[5*n+3]=0;
        J[5*n+4]=0;
        J[5*n+5]=-1;
        J[5*n+6]=0;
        J[5*n+7]=0;
        J[5*n+8]=0;
        
        
        J[6*n+0]=0;
        J[6*n+1]=0;
        J[6*n+2]=0;
        J[6*n+3]=0;
        J[6*n+4]=0;
        J[6*n+5]=0;
        J[6*n+6]=-1;
        J[6*n+7]=0;
        J[6*n+8]=0;
        
        
        J[7*n+0]=0;
        J[7*n+1]=0;
        J[7*n+2]=0;
        J[7*n+3]=0;
        J[7*n+4]=0;
        J[7*n+5]=0;
        J[7*n+6]=0;
        J[7*n+7]=-1;
        J[7*n+8]=0;
        
        
        J[8*n+0]=x[8]*(1-x[3]-x[4]-x[5]);
        J[8*n+1]=0;
        J[8*n+2]=0;
        J[8*n+3]=-x[8]*x[0];
        J[8*n+4]=-x[8]*x[0];
        J[8*n+5]=-x[8]*x[0];
        J[8*n+6]=0;
        J[8*n+7]=0;
        J[8*n+8]=-1+(1-x[3]-x[4]-x[5])*x[0];
        
        
        b[0] = fq[0] - x[0]+(x[8]*(1 - (x[3]+x[4]+x[5]))*x[0]*1+0);
        b[1] = fq[1] - x[1];
        b[2] = fq[2] - x[2];
        b[3] = fq[3] - x[3];
        b[4] = fq[4] - x[4];
        b[5] = fq[5] - x[5];
        b[6] = fq[6] - x[6];
        b[7] = fq[7] - x[7];
        b[8] = fq[8] - x[8]+(x[8]*(1 - (x[3]+x[4]+x[5]))*x[0]*1+0);
        
#if defined BLAS || defined BLAS_OMP
        
        int  *ipiv;
        ipiv=(int *)mCalloc(n, sizeof(int));
        
        dgesv_(&n, &rhs, J, &n, ipiv, b, &n, &info);
        
        free(ipiv);
        
#else
        
        info=MyGaussElimination_gNxRHS(J,b,n,rhs);
        
#endif
        
        if(!info) 
          For(i,n) x[i]-=b[i];
        else 
          Warn_And_Exit("Error in Gauss Elimination Scheme.\n");
        
      }
      
      break;
    }
    case THEFMC: 
    {
      while(my2norm(b,n)>1e-10)
      {
        
        J[0*n+0]=-1+(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[0*n+1]=0;
        J[0*n+2]=0;
        J[0*n+3]=0;
        J[0*n+4]=0;
        J[0*n+5]=(1-x[6]-x[7]-x[8])*x[0]+x[8]*x[0];
        J[0*n+6]=-x[5]*x[0];
        J[0*n+7]=-x[5]*x[0];
        J[0*n+8]=0;
        
        
        J[1*n+0]=0;
        J[1*n+1]=-1;
        J[1*n+2]=0;
        J[1*n+3]=0;
        J[1*n+4]=0;
        J[1*n+5]=0;
        J[1*n+6]=0;
        J[1*n+7]=0;
        J[1*n+8]=0;
        
        
        J[2*n+0]=0;
        J[2*n+1]=0;
        J[2*n+2]=-1;
        J[2*n+3]=0;
        J[2*n+4]=0;
        J[2*n+5]=0;
        J[2*n+6]=0;
        J[2*n+7]=0;
        J[2*n+8]=0;
        
        
        J[3*n+0]=0;
        J[3*n+1]=0;
        J[3*n+2]=0;
        J[3*n+3]=-1;
        J[3*n+4]=0;
        J[3*n+5]=0;
        J[3*n+6]=0;
        J[3*n+7]=0;
        J[3*n+8]=0;
        
        
        J[4*n+0]=0;
        J[4*n+1]=0;
        J[4*n+2]=0;
        J[4*n+3]=0;
        J[4*n+4]=-1;
        J[4*n+5]=0;
        J[4*n+6]=0;
        J[4*n+7]=0;
        J[4*n+8]=0;
        
        
        J[5*n+0]=(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[5*n+1]=0;
        J[5*n+2]=0;
        J[5*n+3]=0;
        J[5*n+4]=0;
        J[5*n+5]=-1+(1-x[6]-x[7]-x[8])*x[0]+x[8]*x[0];
        J[5*n+6]=-x[5]*x[0];
        J[5*n+7]=-x[5]*x[0];
        J[5*n+8]=0;
        
        
        J[6*n+0]=0;
        J[6*n+1]=0;
        J[6*n+2]=0;
        J[6*n+3]=0;
        J[6*n+4]=0;
        J[6*n+5]=0;
        J[6*n+6]=-1;
        J[6*n+7]=0;
        J[6*n+8]=0;
        
        
        J[7*n+0]=0;
        J[7*n+1]=0;
        J[7*n+2]=0;
        J[7*n+3]=0;
        J[7*n+4]=0;
        J[7*n+5]=0;
        J[7*n+6]=0;
        J[7*n+7]=-1;
        J[7*n+8]=0;
        
        
        J[8*n+0]=x[8]*x[5];
        J[8*n+1]=0;
        J[8*n+2]=0;
        J[8*n+3]=0;
        J[8*n+4]=0;
        J[8*n+5]=x[8]*x[0];
        J[8*n+6]=0;
        J[8*n+7]=0;
        J[8*n+8]=-1+x[5]*x[0];
        
        
        b[0] = fq[0] - x[0]+((1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+x[8]*x[5]*x[0]*1+0);
        b[1] = fq[1] - x[1];
        b[2] = fq[2] - x[2];
        b[3] = fq[3] - x[3];
        b[4] = fq[4] - x[4];
        b[5] = fq[5] - x[5]+((1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+x[8]*x[5]*x[0]*1+0);
        b[6] = fq[6] - x[6];
        b[7] = fq[7] - x[7];
        b[8] = fq[8] - x[8]+(x[8]*x[5]*x[0]*1+0);
        
        
#if defined BLAS || defined BLAS_OMP
        
        int  *ipiv;
        ipiv=(int *)mCalloc(n, sizeof(int));
        
        dgesv_(&n, &rhs, J, &n, ipiv, b, &n, &info);
        
        free(ipiv);
        
#else
        
        info=MyGaussElimination_gNxRHS(J,b,n,rhs);
        
#endif
        
        if(!info) 
          For(i,n) x[i]-=b[i];
        else 
          Warn_And_Exit("Error in Gauss Elimination Scheme.\n");
        
      }
      
      break;
    }
    case THENC: 
    {
      while(my2norm(b,n)>1e-10)
      {
        
        J[0*n+0]=-1+(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[0*n+1]=0;
        J[0*n+2]=0;
        J[0*n+3]=0;
        J[0*n+4]=0;
        J[0*n+5]=(1-x[6]-x[7]-x[8])*x[0]+x[8]*x[0];
        J[0*n+6]=-x[5]*x[0];
        J[0*n+7]=-x[5]*x[0];
        J[0*n+8]=0;
        
        
        J[1*n+0]=0;
        J[1*n+1]=-1;
        J[1*n+2]=0;
        J[1*n+3]=0;
        J[1*n+4]=0;
        J[1*n+5]=0;
        J[1*n+6]=0;
        J[1*n+7]=0;
        J[1*n+8]=0;
        
        
        J[2*n+0]=0;
        J[2*n+1]=0;
        J[2*n+2]=-1;
        J[2*n+3]=0;
        J[2*n+4]=0;
        J[2*n+5]=0;
        J[2*n+6]=0;
        J[2*n+7]=0;
        J[2*n+8]=0;
        
        
        J[3*n+0]=0;
        J[3*n+1]=0;
        J[3*n+2]=0;
        J[3*n+3]=-1;
        J[3*n+4]=0;
        J[3*n+5]=0;
        J[3*n+6]=0;
        J[3*n+7]=0;
        J[3*n+8]=0;
        
        
        J[4*n+0]=0;
        J[4*n+1]=0;
        J[4*n+2]=0;
        J[4*n+3]=0;
        J[4*n+4]=-1;
        J[4*n+5]=0;
        J[4*n+6]=0;
        J[4*n+7]=0;
        J[4*n+8]=0;
        
        
        J[5*n+0]=(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[5*n+1]=0;
        J[5*n+2]=0;
        J[5*n+3]=0;
        J[5*n+4]=0;
        J[5*n+5]=-1+(1-x[6]-x[7]-x[8])*x[0]+x[8]*x[0];
        J[5*n+6]=-x[5]*x[0];
        J[5*n+7]=-x[5]*x[0];
        J[5*n+8]=0;
        
        
        J[6*n+0]=0;
        J[6*n+1]=0;
        J[6*n+2]=0;
        J[6*n+3]=0;
        J[6*n+4]=0;
        J[6*n+5]=0;
        J[6*n+6]=-1;
        J[6*n+7]=0;
        J[6*n+8]=0;
        
        
        J[7*n+0]=0;
        J[7*n+1]=0;
        J[7*n+2]=0;
        J[7*n+3]=0;
        J[7*n+4]=0;
        J[7*n+5]=0;
        J[7*n+6]=0;
        J[7*n+7]=-1;
        J[7*n+8]=0;
        
        
        J[8*n+0]=x[8]*x[5];
        J[8*n+1]=0;
        J[8*n+2]=0;
        J[8*n+3]=0;
        J[8*n+4]=0;
        J[8*n+5]=x[8]*x[0];
        J[8*n+6]=0;
        J[8*n+7]=0;
        J[8*n+8]=-1+x[5]*x[0];
        
        
        b[0] = fq[0] - x[0]+((1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+x[8]*x[5]*x[0]*1+0);
        b[1] = fq[1] - x[1];
        b[2] = fq[2] - x[2];
        b[3] = fq[3] - x[3];
        b[4] = fq[4] - x[4];
        b[5] = fq[5] - x[5]+((1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+x[8]*x[5]*x[0]*1+0);
        b[6] = fq[6] - x[6];
        b[7] = fq[7] - x[7];
        b[8] = fq[8] - x[8]+(x[8]*x[5]*x[0]*1+0);
        
#if defined BLAS || defined BLAS_OMP
        
        int  *ipiv;
        ipiv=(int *)mCalloc(n, sizeof(int));
        
        dgesv_(&n, &rhs, J, &n, ipiv, b, &n, &info);
        
        free(ipiv);
        
#else
        
        info=MyGaussElimination_gNxRHS(J,b,n,rhs);
        
#endif
        
        if(!info) 
          For(i,n) x[i]-=b[i];
        else 
          Warn_And_Exit("Error in Gauss Elimination Scheme.\n");
        
      }
      
      break;
    }
    case THBAPPC: 
    {
      while(my2norm(b,n)>1e-10)
      {
        
        J[0*n+0]=-1+x[8]*(1-x[3]-x[4]-x[5])+(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[0*n+1]=0;
        J[0*n+2]=0;
        J[0*n+3]=-x[8]*x[0];
        J[0*n+4]=-x[8]*x[0];
        J[0*n+5]=(1-x[6]-x[7]-x[8])*x[0];
        J[0*n+6]=-x[5]*x[0];
        J[0*n+7]=-x[5]*x[0];
        J[0*n+8]=(1-x[3]-x[4]-x[5])*x[0];
        
        
        J[1*n+0]=0;
        J[1*n+1]=-1;
        J[1*n+2]=0;
        J[1*n+3]=0;
        J[1*n+4]=0;
        J[1*n+5]=0;
        J[1*n+6]=0;
        J[1*n+7]=0;
        J[1*n+8]=0;
        
        
        J[2*n+0]=0;
        J[2*n+1]=0;
        J[2*n+2]=-1;
        J[2*n+3]=0;
        J[2*n+4]=0;
        J[2*n+5]=0;
        J[2*n+6]=0;
        J[2*n+7]=0;
        J[2*n+8]=0;
        
        
        J[3*n+0]=0;
        J[3*n+1]=0;
        J[3*n+2]=0;
        J[3*n+3]=-1;
        J[3*n+4]=0;
        J[3*n+5]=0;
        J[3*n+6]=0;
        J[3*n+7]=0;
        J[3*n+8]=0;
        
        
        J[4*n+0]=0;
        J[4*n+1]=0;
        J[4*n+2]=0;
        J[4*n+3]=0;
        J[4*n+4]=-1;
        J[4*n+5]=0;
        J[4*n+6]=0;
        J[4*n+7]=0;
        J[4*n+8]=0;
        
        
        J[5*n+0]=(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[5*n+1]=0;
        J[5*n+2]=0;
        J[5*n+3]=0;
        J[5*n+4]=0;
        J[5*n+5]=-1+(1-x[6]-x[7]-x[8])*x[0]+x[8]*x[0];
        J[5*n+6]=-x[5]*x[0];
        J[5*n+7]=-x[5]*x[0];
        J[5*n+8]=0;
        
        
        J[6*n+0]=0;
        J[6*n+1]=0;
        J[6*n+2]=0;
        J[6*n+3]=0;
        J[6*n+4]=0;
        J[6*n+5]=0;
        J[6*n+6]=-1;
        J[6*n+7]=0;
        J[6*n+8]=0;
        
        
        J[7*n+0]=0;
        J[7*n+1]=0;
        J[7*n+2]=0;
        J[7*n+3]=0;
        J[7*n+4]=0;
        J[7*n+5]=0;
        J[7*n+6]=0;
        J[7*n+7]=-1;
        J[7*n+8]=0;
        
        
        J[8*n+0]=x[8]*(1-x[3]-x[4]-x[5])+x[8]*x[5];
        J[8*n+1]=0;
        J[8*n+2]=0;
        J[8*n+3]=-x[8]*x[0];
        J[8*n+4]=-x[8]*x[0];
        J[8*n+5]=0;
        J[8*n+6]=0;
        J[8*n+7]=0;
        J[8*n+8]=-1+(1-x[3]-x[4]-x[5])*x[0]+x[5]*x[0];
        
        
        b[0] = fq[0] - x[0]+(x[8]*(1 - (x[3]+x[4]+x[5]))*x[0]*1+(1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+x[8]*x[5]*x[0]*1+0);
        b[1] = fq[1] - x[1];
        b[2] = fq[2] - x[2];
        b[3] = fq[3] - x[3];
        b[4] = fq[4] - x[4];
        b[5] = fq[5] - x[5]+((1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+x[8]*x[5]*x[0]*1+0);
        b[6] = fq[6] - x[6];
        b[7] = fq[7] - x[7];
        b[8] = fq[8] - x[8]+(x[8]*(1 - (x[3]+x[4]+x[5]))*x[0]*1+x[8]*x[5]*x[0]*1+0);
        
        
#if defined BLAS || defined BLAS_OMP
        
        int  *ipiv;
        ipiv=(int *)mCalloc(n, sizeof(int));
        
        dgesv_(&n, &rhs, J, &n, ipiv, b, &n, &info);
        
        free(ipiv);
        
#else
        
        info=MyGaussElimination_gNxRHS(J,b,n,rhs);
        
#endif
        
        if(!info) 
          For(i,n) x[i]-=b[i];
        else 
          Warn_And_Exit("Error in Gauss Elimination Scheme.\n");
        
      }
      
      break;
    }
    case THAYNC: 
    {
      while(my2norm(b,n)>1e-10)
      {
        
        J[0*n+0]=-1+x[8]*(1-x[3]-x[4]-x[5])+(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[0*n+1]=0;
        J[0*n+2]=0;
        J[0*n+3]=-x[8]*x[0];
        J[0*n+4]=-x[8]*x[0];
        J[0*n+5]=(1-x[6]-x[7]-x[8])*x[0];
        J[0*n+6]=-x[5]*x[0];
        J[0*n+7]=-x[5]*x[0];
        J[0*n+8]=(1-x[3]-x[4]-x[5])*x[0];
        
        
        J[1*n+0]=0;
        J[1*n+1]=-1;
        J[1*n+2]=0;
        J[1*n+3]=0;
        J[1*n+4]=0;
        J[1*n+5]=0;
        J[1*n+6]=0;
        J[1*n+7]=0;
        J[1*n+8]=0;
        
        
        J[2*n+0]=0;
        J[2*n+1]=0;
        J[2*n+2]=-1;
        J[2*n+3]=0;
        J[2*n+4]=0;
        J[2*n+5]=0;
        J[2*n+6]=0;
        J[2*n+7]=0;
        J[2*n+8]=0;
        
        
        J[3*n+0]=0;
        J[3*n+1]=0;
        J[3*n+2]=0;
        J[3*n+3]=-1;
        J[3*n+4]=0;
        J[3*n+5]=0;
        J[3*n+6]=0;
        J[3*n+7]=0;
        J[3*n+8]=0;
        
        
        J[4*n+0]=0;
        J[4*n+1]=0;
        J[4*n+2]=0;
        J[4*n+3]=0;
        J[4*n+4]=-1;
        J[4*n+5]=0;
        J[4*n+6]=0;
        J[4*n+7]=0;
        J[4*n+8]=0;
        
        
        J[5*n+0]=(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[5*n+1]=0;
        J[5*n+2]=0;
        J[5*n+3]=0;
        J[5*n+4]=0;
        J[5*n+5]=-1+(1-x[6]-x[7]-x[8])*x[0]+x[8]*x[0];
        J[5*n+6]=-x[5]*x[0];
        J[5*n+7]=-x[5]*x[0];
        J[5*n+8]=0;
        
        
        J[6*n+0]=0;
        J[6*n+1]=0;
        J[6*n+2]=0;
        J[6*n+3]=0;
        J[6*n+4]=0;
        J[6*n+5]=0;
        J[6*n+6]=-1;
        J[6*n+7]=0;
        J[6*n+8]=0;
        
        
        J[7*n+0]=0;
        J[7*n+1]=0;
        J[7*n+2]=0;
        J[7*n+3]=0;
        J[7*n+4]=0;
        J[7*n+5]=0;
        J[7*n+6]=0;
        J[7*n+7]=-1;
        J[7*n+8]=0;
        
        
        J[8*n+0]=x[8]*(1-x[3]-x[4]-x[5])+x[8]*x[5];
        J[8*n+1]=0;
        J[8*n+2]=0;
        J[8*n+3]=-x[8]*x[0];
        J[8*n+4]=-x[8]*x[0];
        J[8*n+5]=0;
        J[8*n+6]=0;
        J[8*n+7]=0;
        J[8*n+8]=-1+(1-x[3]-x[4]-x[5])*x[0]+x[5]*x[0];
        
        
        b[0] = fq[0] - x[0]+(x[8]*(1 - (x[3]+x[4]+x[5]))*x[0]+(1 - (x[6]+x[7]+x[8]))*x[5]*x[0]+x[8]*x[5]*x[0]);
        b[1] = fq[1] - x[1];
        b[2] = fq[2] - x[2];
        b[3] = fq[3] - x[3];
        b[4] = fq[4] - x[4];
        b[5] = fq[5] - x[5]+((1 - (x[6]+x[7]+x[8]))*x[5]*x[0]+x[8]*x[5]*x[0]);
        b[6] = fq[6] - x[6];
        b[7] = fq[7] - x[7];
        b[8] = fq[8] - x[8]+(x[8]*(1 - (x[3]+x[4]+x[5]))*x[0]+x[8]*x[5]*x[0]);
        
#if defined BLAS || defined BLAS_OMP
        
        int  *ipiv;
        ipiv=(int *)mCalloc(n, sizeof(int));
        
        dgesv_(&n, &rhs, J, &n, ipiv, b, &n, &info);
        
        free(ipiv);
        
#else
        
        info=MyGaussElimination_gNxRHS(J,b,n,rhs);
        
#endif
        
        if(!info) 
          For(i,n) x[i]-=b[i];
        else 
          Warn_And_Exit("Error in Gauss Elimination Scheme.\n");
        
      }
      
      break;
    }
    case THAMC: 
    {
      while(my2norm(b,n)>1e-10)
      {
        
        J[0*n+0]=-1+(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[0*n+1]=0;
        J[0*n+2]=0;
        J[0*n+3]=0;
        J[0*n+4]=0;
        J[0*n+5]=(1-x[6]-x[7]-x[8])*x[0]+x[8]*x[0];
        J[0*n+6]=-x[5]*x[0];
        J[0*n+7]=-x[5]*x[0];
        J[0*n+8]=0;
        
        
        J[1*n+0]=0;
        J[1*n+1]=-1;
        J[1*n+2]=0;
        J[1*n+3]=0;
        J[1*n+4]=0;
        J[1*n+5]=0;
        J[1*n+6]=0;
        J[1*n+7]=0;
        J[1*n+8]=0;
        
        
        J[2*n+0]=0;
        J[2*n+1]=0;
        J[2*n+2]=-1;
        J[2*n+3]=0;
        J[2*n+4]=0;
        J[2*n+5]=0;
        J[2*n+6]=0;
        J[2*n+7]=0;
        J[2*n+8]=0;
        
        
        J[3*n+0]=0;
        J[3*n+1]=0;
        J[3*n+2]=0;
        J[3*n+3]=-1;
        J[3*n+4]=0;
        J[3*n+5]=0;
        J[3*n+6]=0;
        J[3*n+7]=0;
        J[3*n+8]=0;
        
        
        J[4*n+0]=0;
        J[4*n+1]=0;
        J[4*n+2]=0;
        J[4*n+3]=0;
        J[4*n+4]=-1;
        J[4*n+5]=0;
        J[4*n+6]=0;
        J[4*n+7]=0;
        J[4*n+8]=0;
        
        
        J[5*n+0]=(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[5*n+1]=0;
        J[5*n+2]=0;
        J[5*n+3]=0;
        J[5*n+4]=0;
        J[5*n+5]=-1+(1-x[6]-x[7]-x[8])*x[0]+x[8]*x[0];
        J[5*n+6]=-x[5]*x[0];
        J[5*n+7]=-x[5]*x[0];
        J[5*n+8]=0;
        
        
        J[6*n+0]=0;
        J[6*n+1]=0;
        J[6*n+2]=0;
        J[6*n+3]=0;
        J[6*n+4]=0;
        J[6*n+5]=0;
        J[6*n+6]=-1;
        J[6*n+7]=0;
        J[6*n+8]=0;
        
        
        J[7*n+0]=0;
        J[7*n+1]=0;
        J[7*n+2]=0;
        J[7*n+3]=0;
        J[7*n+4]=0;
        J[7*n+5]=0;
        J[7*n+6]=0;
        J[7*n+7]=-1;
        J[7*n+8]=0;
        
        
        J[8*n+0]=x[8]*x[5];
        J[8*n+1]=0;
        J[8*n+2]=0;
        J[8*n+3]=0;
        J[8*n+4]=0;
        J[8*n+5]=x[8]*x[0];
        J[8*n+6]=0;
        J[8*n+7]=0;
        J[8*n+8]=-1+x[5]*x[0];
        
        
        b[0] = fq[0] - x[0]+((1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+x[8]*x[5]*x[0]*1+0);
        b[1] = fq[1] - x[1];
        b[2] = fq[2] - x[2];
        b[3] = fq[3] - x[3];
        b[4] = fq[4] - x[4];
        b[5] = fq[5] - x[5]+((1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+x[8]*x[5]*x[0]*1+0);
        b[6] = fq[6] - x[6];
        b[7] = fq[7] - x[7];
        b[8] = fq[8] - x[8]+(x[8]*x[5]*x[0]*1+0);
        
#if defined BLAS || defined BLAS_OMP
        
        int  *ipiv;
        ipiv=(int *)mCalloc(n, sizeof(int));
        
        dgesv_(&n, &rhs, J, &n, ipiv, b, &n, &info);
        
        free(ipiv);
        
#else
        
        info=MyGaussElimination_gNxRHS(J,b,n,rhs);
        
#endif
        
        if(!info) 
          For(i,n) x[i]-=b[i];
        else 
          Warn_And_Exit("Error in Gauss Elimination Scheme.\n");
        
      }
      
      break;
    }
    case THAFMC: 
    {
      while(my2norm(b,n)>1e-10)
      {
        
        J[0*n+0]=-1+(1-x[6]-x[7]-x[8])*x[5];
        J[0*n+1]=0;
        J[0*n+2]=0;
        J[0*n+3]=0;
        J[0*n+4]=0;
        J[0*n+5]=(1-x[6]-x[7]-x[8])*x[0];
        J[0*n+6]=-x[5]*x[0];
        J[0*n+7]=-x[5]*x[0];
        J[0*n+8]=-x[5]*x[0];
        
        
        J[1*n+0]=0;
        J[1*n+1]=-1;
        J[1*n+2]=0;
        J[1*n+3]=0;
        J[1*n+4]=0;
        J[1*n+5]=0;
        J[1*n+6]=0;
        J[1*n+7]=0;
        J[1*n+8]=0;
        
        
        J[2*n+0]=0;
        J[2*n+1]=0;
        J[2*n+2]=-1;
        J[2*n+3]=0;
        J[2*n+4]=0;
        J[2*n+5]=0;
        J[2*n+6]=0;
        J[2*n+7]=0;
        J[2*n+8]=0;
        
        
        J[3*n+0]=0;
        J[3*n+1]=0;
        J[3*n+2]=0;
        J[3*n+3]=-1;
        J[3*n+4]=0;
        J[3*n+5]=0;
        J[3*n+6]=0;
        J[3*n+7]=0;
        J[3*n+8]=0;
        
        
        J[4*n+0]=0;
        J[4*n+1]=0;
        J[4*n+2]=0;
        J[4*n+3]=0;
        J[4*n+4]=-1;
        J[4*n+5]=0;
        J[4*n+6]=0;
        J[4*n+7]=0;
        J[4*n+8]=0;
        
        
        J[5*n+0]=(1-x[6]-x[7]-x[8])*x[5];
        J[5*n+1]=0;
        J[5*n+2]=0;
        J[5*n+3]=0;
        J[5*n+4]=0;
        J[5*n+5]=-1+(1-x[6]-x[7]-x[8])*x[0];
        J[5*n+6]=-x[5]*x[0];
        J[5*n+7]=-x[5]*x[0];
        J[5*n+8]=-x[5]*x[0];
        
        
        J[6*n+0]=0;
        J[6*n+1]=0;
        J[6*n+2]=0;
        J[6*n+3]=0;
        J[6*n+4]=0;
        J[6*n+5]=0;
        J[6*n+6]=-1;
        J[6*n+7]=0;
        J[6*n+8]=0;
        
        
        J[7*n+0]=0;
        J[7*n+1]=0;
        J[7*n+2]=0;
        J[7*n+3]=0;
        J[7*n+4]=0;
        J[7*n+5]=0;
        J[7*n+6]=0;
        J[7*n+7]=-1;
        J[7*n+8]=0;
        
        
        J[8*n+0]=0;
        J[8*n+1]=0;
        J[8*n+2]=0;
        J[8*n+3]=0;
        J[8*n+4]=0;
        J[8*n+5]=0;
        J[8*n+6]=0;
        J[8*n+7]=0;
        J[8*n+8]=-1;
        
        
        b[0] = fq[0] - x[0]+((1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+0);
        b[1] = fq[1] - x[1];
        b[2] = fq[2] - x[2];
        b[3] = fq[3] - x[3];
        b[4] = fq[4] - x[4];
        b[5] = fq[5] - x[5]+((1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+0);
        b[6] = fq[6] - x[6];
        b[7] = fq[7] - x[7];
        b[8] = fq[8] - x[8];
        
        
#if defined BLAS || defined BLAS_OMP
        
        int  *ipiv;
        ipiv=(int *)mCalloc(n, sizeof(int));
        
        dgesv_(&n, &rhs, J, &n, ipiv, b, &n, &info);
        
        free(ipiv);
        
#else
        
        info=MyGaussElimination_gNxRHS(J,b,n,rhs);
        
#endif
        
        if(!info) 
          For(i,n) x[i]-=b[i];
        else 
          Warn_And_Exit("Error in Gauss Elimination Scheme.\n");
        
      }
      
      break;
    }
    case BLNC: 
    {
      J[0*n+0]=-1+x[8]*(1-x[3]-x[4]-x[5])+x[8]*x[5];
      J[0*n+1]=0;
      J[0*n+2]=0;
      J[0*n+3]=-x[8]*x[0];
      J[0*n+4]=-x[8]*x[0];
      J[0*n+5]=0;
      J[0*n+6]=0;
      J[0*n+7]=0;
      J[0*n+8]=(1-x[3]-x[4]-x[5])*x[0]+x[5]*x[0];
      
      
      J[1*n+0]=0;
      J[1*n+1]=-1;
      J[1*n+2]=0;
      J[1*n+3]=0;
      J[1*n+4]=0;
      J[1*n+5]=0;
      J[1*n+6]=0;
      J[1*n+7]=0;
      J[1*n+8]=0;
      
      
      J[2*n+0]=0;
      J[2*n+1]=0;
      J[2*n+2]=-1;
      J[2*n+3]=0;
      J[2*n+4]=0;
      J[2*n+5]=0;
      J[2*n+6]=0;
      J[2*n+7]=0;
      J[2*n+8]=0;
      
      
      J[3*n+0]=0;
      J[3*n+1]=0;
      J[3*n+2]=0;
      J[3*n+3]=-1;
      J[3*n+4]=0;
      J[3*n+5]=0;
      J[3*n+6]=0;
      J[3*n+7]=0;
      J[3*n+8]=0;
      
      
      J[4*n+0]=0;
      J[4*n+1]=0;
      J[4*n+2]=0;
      J[4*n+3]=0;
      J[4*n+4]=-1;
      J[4*n+5]=0;
      J[4*n+6]=0;
      J[4*n+7]=0;
      J[4*n+8]=0;
      
      
      J[5*n+0]=x[8]*x[5];
      J[5*n+1]=0;
      J[5*n+2]=0;
      J[5*n+3]=0;
      J[5*n+4]=0;
      J[5*n+5]=-1+x[8]*x[0];
      J[5*n+6]=0;
      J[5*n+7]=0;
      J[5*n+8]=x[5]*x[0];
      
      
      J[6*n+0]=0;
      J[6*n+1]=0;
      J[6*n+2]=0;
      J[6*n+3]=0;
      J[6*n+4]=0;
      J[6*n+5]=0;
      J[6*n+6]=-1;
      J[6*n+7]=0;
      J[6*n+8]=0;
      
      
      J[7*n+0]=0;
      J[7*n+1]=0;
      J[7*n+2]=0;
      J[7*n+3]=0;
      J[7*n+4]=0;
      J[7*n+5]=0;
      J[7*n+6]=0;
      J[7*n+7]=-1;
      J[7*n+8]=0;
      
      
      J[8*n+0]=x[8]*(1-x[3]-x[4]-x[5])+x[8]*x[5];
      J[8*n+1]=0;
      J[8*n+2]=0;
      J[8*n+3]=-x[8]*x[0];
      J[8*n+4]=-x[8]*x[0];
      J[8*n+5]=0;
      J[8*n+6]=0;
      J[8*n+7]=0;
      J[8*n+8]=-1+(1-x[3]-x[4]-x[5])*x[0]+x[5]*x[0];
      
      b[0] = fq[0] - x[0]+(x[8]*(1 - (x[3]+x[4]+x[5]))*x[0]+x[8]*x[5]*x[0]);
      b[1] = fq[1] - x[1]+(0);
      b[2] = fq[2] - x[2]+(0);
      b[3] = fq[3] - x[3]+(0);
      b[4] = fq[4] - x[4]+(0);
      b[5] = fq[5] - x[5]+(x[8]*x[5]*x[0]);
      b[6] = fq[6] - x[6]+(0);
      b[7] = fq[7] - x[7]+(0);
      b[8] = fq[8] - x[8]+(x[8]*(1 - (x[3]+x[4]+x[5]))*x[0]+x[8]*x[5]*x[0]);
      
#if defined BLAS || defined BLAS_OMP
      
      int  *ipiv;
      ipiv=(int *)mCalloc(n, sizeof(int));
      
      dgesv_(&n, &rhs, J, &n, ipiv, b, &n, &info);
      
      free(ipiv);
      
#else
      
      info=MyGaussElimination_gNxRHS(J,b,n,rhs);
      
#endif
      
      if(!info) 
        For(i,n) x[i]-=b[i];
      else 
        Warn_And_Exit("Error in Gauss Elimination Scheme.\n");
      
      
      
      break;
    }
    case CHMC: 
    {
      while(my2norm(b,n)>1e-10)
      {
        
        J[0*n+0]=-1+x[8]*(1-x[3]-x[4]-x[5])+x[8]*x[5];
        J[0*n+1]=0;
        J[0*n+2]=0;
        J[0*n+3]=-x[8]*x[0];
        J[0*n+4]=-x[8]*x[0];
        J[0*n+5]=0;
        J[0*n+6]=0;
        J[0*n+7]=0;
        J[0*n+8]=(1-x[3]-x[4]-x[5])*x[0]+x[5]*x[0];
        
        
        J[1*n+0]=0;
        J[1*n+1]=-1;
        J[1*n+2]=0;
        J[1*n+3]=0;
        J[1*n+4]=0;
        J[1*n+5]=0;
        J[1*n+6]=0;
        J[1*n+7]=0;
        J[1*n+8]=0;
        
        
        J[2*n+0]=0;
        J[2*n+1]=0;
        J[2*n+2]=-1;
        J[2*n+3]=0;
        J[2*n+4]=0;
        J[2*n+5]=0;
        J[2*n+6]=0;
        J[2*n+7]=0;
        J[2*n+8]=0;
        
        
        J[3*n+0]=0;
        J[3*n+1]=0;
        J[3*n+2]=0;
        J[3*n+3]=-1;
        J[3*n+4]=0;
        J[3*n+5]=0;
        J[3*n+6]=0;
        J[3*n+7]=0;
        J[3*n+8]=0;
        
        
        J[4*n+0]=0;
        J[4*n+1]=0;
        J[4*n+2]=0;
        J[4*n+3]=0;
        J[4*n+4]=-1;
        J[4*n+5]=0;
        J[4*n+6]=0;
        J[4*n+7]=0;
        J[4*n+8]=0;
        
        
        J[5*n+0]=x[8]*x[5];
        J[5*n+1]=0;
        J[5*n+2]=0;
        J[5*n+3]=0;
        J[5*n+4]=0;
        J[5*n+5]=-1+x[8]*x[0];
        J[5*n+6]=0;
        J[5*n+7]=0;
        J[5*n+8]=x[5]*x[0];
        
        
        J[6*n+0]=0;
        J[6*n+1]=0;
        J[6*n+2]=0;
        J[6*n+3]=0;
        J[6*n+4]=0;
        J[6*n+5]=0;
        J[6*n+6]=-1;
        J[6*n+7]=0;
        J[6*n+8]=0;
        
        
        J[7*n+0]=0;
        J[7*n+1]=0;
        J[7*n+2]=0;
        J[7*n+3]=0;
        J[7*n+4]=0;
        J[7*n+5]=0;
        J[7*n+6]=0;
        J[7*n+7]=-1;
        J[7*n+8]=0;
        
        
        J[8*n+0]=x[8]*(1-x[3]-x[4]-x[5])+x[8]*x[5];
        J[8*n+1]=0;
        J[8*n+2]=0;
        J[8*n+3]=-x[8]*x[0];
        J[8*n+4]=-x[8]*x[0];
        J[8*n+5]=0;
        J[8*n+6]=0;
        J[8*n+7]=0;
        J[8*n+8]=-1+(1-x[3]-x[4]-x[5])*x[0]+x[5]*x[0];
        
        
        b[0] = fq[0] - x[0]+(x[8]*(1 - (x[3]+x[4]+x[5]))*x[0]+x[8]*x[5]*x[0]);
        b[1] = fq[1] - x[1];
        b[2] = fq[2] - x[2];
        b[3] = fq[3] - x[3];
        b[4] = fq[4] - x[4];
        b[5] = fq[5] - x[5]+(x[8]*x[5]*x[0]+0);
        b[6] = fq[6] - x[6];
        b[7] = fq[7] - x[7];
        b[8] = fq[8] - x[8]+(x[8]*(1 - (x[3]+x[4]+x[5]))*x[0]+x[8]*x[5]*x[0]);
        
#if defined BLAS || defined BLAS_OMP
        
        int  *ipiv;
        ipiv=(int *)mCalloc(n, sizeof(int));
        
        dgesv_(&n, &rhs, J, &n, ipiv, b, &n, &info);
        
        free(ipiv);
        
#else
        
        info=MyGaussElimination_gNxRHS(J,b,n,rhs);
        
#endif
        
        if(!info) 
          For(i,n) x[i]-=b[i];
        else 
          Warn_And_Exit("Error in Gauss Elimination Scheme.\n");
        
      }
      
      break;
    }
    case TRMC: 
    {
      while(my2norm(b,n)>1e-10)
      {
        
        J[0*n+0]=-1+(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[0*n+1]=0;
        J[0*n+2]=0;
        J[0*n+3]=0;
        J[0*n+4]=0;
        J[0*n+5]=(1-x[6]-x[7]-x[8])*x[0]+x[8]*x[0];
        J[0*n+6]=-x[5]*x[0];
        J[0*n+7]=-x[5]*x[0];
        J[0*n+8]=0;
        
        
        J[1*n+0]=0;
        J[1*n+1]=-1;
        J[1*n+2]=0;
        J[1*n+3]=0;
        J[1*n+4]=0;
        J[1*n+5]=0;
        J[1*n+6]=0;
        J[1*n+7]=0;
        J[1*n+8]=0;
        
        
        J[2*n+0]=0;
        J[2*n+1]=0;
        J[2*n+2]=-1;
        J[2*n+3]=0;
        J[2*n+4]=0;
        J[2*n+5]=0;
        J[2*n+6]=0;
        J[2*n+7]=0;
        J[2*n+8]=0;
        
        
        J[3*n+0]=0;
        J[3*n+1]=0;
        J[3*n+2]=0;
        J[3*n+3]=-1;
        J[3*n+4]=0;
        J[3*n+5]=0;
        J[3*n+6]=0;
        J[3*n+7]=0;
        J[3*n+8]=0;
        
        
        J[4*n+0]=0;
        J[4*n+1]=0;
        J[4*n+2]=0;
        J[4*n+3]=0;
        J[4*n+4]=-1;
        J[4*n+5]=0;
        J[4*n+6]=0;
        J[4*n+7]=0;
        J[4*n+8]=0;
        
        
        J[5*n+0]=(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[5*n+1]=0;
        J[5*n+2]=0;
        J[5*n+3]=0;
        J[5*n+4]=0;
        J[5*n+5]=-1+(1-x[6]-x[7]-x[8])*x[0]+x[8]*x[0];
        J[5*n+6]=-x[5]*x[0];
        J[5*n+7]=-x[5]*x[0];
        J[5*n+8]=0;
        
        
        J[6*n+0]=0;
        J[6*n+1]=0;
        J[6*n+2]=0;
        J[6*n+3]=0;
        J[6*n+4]=0;
        J[6*n+5]=0;
        J[6*n+6]=-1;
        J[6*n+7]=0;
        J[6*n+8]=0;
        
        
        J[7*n+0]=0;
        J[7*n+1]=0;
        J[7*n+2]=0;
        J[7*n+3]=0;
        J[7*n+4]=0;
        J[7*n+5]=0;
        J[7*n+6]=0;
        J[7*n+7]=-1;
        J[7*n+8]=0;
        
        
        J[8*n+0]=x[8]*x[5];
        J[8*n+1]=0;
        J[8*n+2]=0;
        J[8*n+3]=0;
        J[8*n+4]=0;
        J[8*n+5]=x[8]*x[0];
        J[8*n+6]=0;
        J[8*n+7]=0;
        J[8*n+8]=-1+x[5]*x[0];
        
        
        b[0] = fq[0] - x[0]+((1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+x[8]*x[5]*x[0]*1+0);
        b[1] = fq[1] - x[1];
        b[2] = fq[2] - x[2];
        b[3] = fq[3] - x[3];
        b[4] = fq[4] - x[4];
        b[5] = fq[5] - x[5]+((1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+x[8]*x[5]*x[0]*1+0);
        b[6] = fq[6] - x[6];
        b[7] = fq[7] - x[7];
        b[8] = fq[8] - x[8]+(x[8]*x[5]*x[0]*1+0);
        
#if defined BLAS || defined BLAS_OMP
        
        int  *ipiv;
        ipiv=(int *)mCalloc(n, sizeof(int));
        
        dgesv_(&n, &rhs, J, &n, ipiv, b, &n, &info);
        
        free(ipiv);
        
#else
        
        info=MyGaussElimination_gNxRHS(J,b,n,rhs);
        
#endif
        
        if(!info) 
          For(i,n) x[i]-=b[i];
        else 
          Warn_And_Exit("Error in Gauss Elimination Scheme.\n");
        
      }
      
      break;
    }
    case SCOMC: 
    {
      while(my2norm(b,n)>1e-10)
      {
        
        J[0*n+0]=0;J[0*n+0]=-1+x[8]*(1-x[3]-x[4]-x[5])+x[8]*x[5]+x[8]*x[4];
        J[0*n+1]=0;
        J[0*n+2]=0;
        J[0*n+3]=-x[8]*x[0];
        J[0*n+4]=0;
        J[0*n+5]=0;
        J[0*n+6]=0;
        J[0*n+7]=0;
        J[0*n+8]=(1-x[3]-x[4]-x[5])*x[0]+x[5]*x[0]+x[4]*x[0];
        
        
        J[1*n+0]=0;
        J[1*n+1]=-1;
        J[1*n+2]=0;
        J[1*n+3]=0;
        J[1*n+4]=0;
        J[1*n+5]=0;
        J[1*n+6]=0;
        J[1*n+7]=0;
        J[1*n+8]=0;
        
        
        J[2*n+0]=0;
        J[2*n+1]=0;
        J[2*n+2]=-1;
        J[2*n+3]=0;
        J[2*n+4]=0;
        J[2*n+5]=0;
        J[2*n+6]=0;
        J[2*n+7]=0;
        J[2*n+8]=0;
        
        
        J[3*n+0]=0;
        J[3*n+1]=0;
        J[3*n+2]=0;
        J[3*n+3]=-1;
        J[3*n+4]=0;
        J[3*n+5]=0;
        J[3*n+6]=0;
        J[3*n+7]=0;
        J[3*n+8]=0;
        
        
        J[4*n+0]=x[8]*x[4];
        J[4*n+1]=0;
        J[4*n+2]=0;
        J[4*n+3]=0;
        J[4*n+4]=-1+x[8]*x[0];
        J[4*n+5]=0;
        J[4*n+6]=0;
        J[4*n+7]=0;
        J[4*n+8]=x[4]*x[0];
        
        
        J[5*n+0]=x[8]*x[5];
        J[5*n+1]=0;
        J[5*n+2]=0;
        J[5*n+3]=0;
        J[5*n+4]=0;
        J[5*n+5]=-1+x[8]*x[0];
        J[5*n+6]=0;
        J[5*n+7]=0;
        J[5*n+8]=x[5]*x[0];
        
        
        J[6*n+0]=0;
        J[6*n+1]=0;
        J[6*n+2]=0;
        J[6*n+3]=0;
        J[6*n+4]=0;
        J[6*n+5]=0;
        J[6*n+6]=-1;
        J[6*n+7]=0;
        J[6*n+8]=0;
        
        
        J[7*n+0]=0;
        J[7*n+1]=0;
        J[7*n+2]=0;
        J[7*n+3]=0;
        J[7*n+4]=0;
        J[7*n+5]=0;
        J[7*n+6]=0;
        J[7*n+7]=-1;
        J[7*n+8]=0;
        
        
        J[8*n+0]=x[8]*(1-x[3]-x[4]-x[5])+x[8]*x[5]+x[8]*x[4];
        J[8*n+1]=0;
        J[8*n+2]=0;
        J[8*n+3]=-x[8]*x[0];
        J[8*n+4]=0;
        J[8*n+5]=0;
        J[8*n+6]=0;
        J[8*n+7]=0;
        J[8*n+8]=-1+(1-x[3]-x[4]-x[5])*x[0]+x[5]*x[0]+x[4]*x[0];
        
        
        b[0] = fq[0] - x[0]+(x[8]*(1 - (x[3]+x[4]+x[5]))*x[0]*1+x[8]*x[5]*x[0]*1+x[8]*x[4]*x[0]*1+0);
        b[1] = fq[1] - x[1];
        b[2] = fq[2] - x[2];
        b[3] = fq[3] - x[3];
        b[4] = fq[4] - x[4]+(x[8]*x[4]*x[0]*1+0);
        b[5] = fq[5] - x[5]+(x[8]*x[5]*x[0]*1+0);
        b[6] = fq[6] - x[6];
        b[7] = fq[7] - x[7];
        b[8] = fq[8] - x[8]+(x[8]*(1 - (x[3]+x[4]+x[5]))*x[0]*1+x[8]*x[5]*x[0]*1+x[8]*x[4]*x[0]*1+0);
        
#if defined BLAS || defined BLAS_OMP
        
        int  *ipiv;
        ipiv=(int *)mCalloc(n, sizeof(int));
        
        dgesv_(&n, &rhs, J, &n, ipiv, b, &n, &info);
        
        free(ipiv);
        
#else
        
        info=MyGaussElimination_gNxRHS(J,b,n,rhs);
        
#endif
        
        if(!info) 
          For(i,n) x[i]-=b[i];
        else 
          Warn_And_Exit("Error in Gauss Elimination Scheme.\n");
        
      }
      
      break;
    }
    case THMC: 
    {
      while(my2norm(b,n)>1e-10)
      {
        
        J[0*n+0]=-1+x[8]*(1-x[3]-x[4]-x[5])+(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5]+x[8]*x[3];
        J[0*n+1]=0;
        J[0*n+2]=0;
        J[0*n+3]=0;
        J[0*n+4]=-x[8]*x[0];
        J[0*n+5]=(1-x[6]-x[7]-x[8])*x[0];
        J[0*n+6]=-x[5]*x[0];
        J[0*n+7]=-x[5]*x[0];
        J[0*n+8]=(1-x[3]-x[4]-x[5])*x[0]+x[3]*x[0];
        
        
        J[1*n+0]=0;
        J[1*n+1]=-1;
        J[1*n+2]=0;
        J[1*n+3]=0;
        J[1*n+4]=0;
        J[1*n+5]=0;
        J[1*n+6]=0;
        J[1*n+7]=0;
        J[1*n+8]=0;
        
        
        J[2*n+0]=0;
        J[2*n+1]=0;
        J[2*n+2]=-1;
        J[2*n+3]=0;
        J[2*n+4]=0;
        J[2*n+5]=0;
        J[2*n+6]=0;
        J[2*n+7]=0;
        J[2*n+8]=0;
        
        
        J[3*n+0]=x[8]*x[3];
        J[3*n+1]=0;
        J[3*n+2]=0;
        J[3*n+3]=-1+x[8]*x[0];
        J[3*n+4]=0;
        J[3*n+5]=0;
        J[3*n+6]=0;
        J[3*n+7]=0;
        J[3*n+8]=x[3]*x[0];
        
        
        J[4*n+0]=0;
        J[4*n+1]=0;
        J[4*n+2]=0;
        J[4*n+3]=0;
        J[4*n+4]=-1;
        J[4*n+5]=0;
        J[4*n+6]=0;
        J[4*n+7]=0;
        J[4*n+8]=0;
        
        
        J[5*n+0]=(1-x[6]-x[7]-x[8])*x[5]+x[8]*x[5];
        J[5*n+1]=0;
        J[5*n+2]=0;
        J[5*n+3]=0;
        J[5*n+4]=0;
        J[5*n+5]=-1+(1-x[6]-x[7]-x[8])*x[0]+x[8]*x[0];
        J[5*n+6]=-x[5]*x[0];
        J[5*n+7]=-x[5]*x[0];
        J[5*n+8]=0;
        
        
        J[6*n+0]=0;
        J[6*n+1]=0;
        J[6*n+2]=0;
        J[6*n+3]=0;
        J[6*n+4]=0;
        J[6*n+5]=0;
        J[6*n+6]=-1;
        J[6*n+7]=0;
        J[6*n+8]=0;
        
        
        J[7*n+0]=0;
        J[7*n+1]=0;
        J[7*n+2]=0;
        J[7*n+3]=0;
        J[7*n+4]=0;
        J[7*n+5]=0;
        J[7*n+6]=0;
        J[7*n+7]=-1;
        J[7*n+8]=0;
        
        
        J[8*n+0]=x[8]*(1-x[3]-x[4]-x[5])+x[8]*x[5]+x[8]*x[3];
        J[8*n+1]=0;
        J[8*n+2]=0;
        J[8*n+3]=0;
        J[8*n+4]=-x[8]*x[0];
        J[8*n+5]=0;
        J[8*n+6]=0;
        J[8*n+7]=0;
        J[8*n+8]=-1+(1-x[3]-x[4]-x[5])*x[0]+x[5]*x[0]+x[3]*x[0];
        
        
        b[0] = fq[0] - x[0]+(x[8]*(1 - (x[3]+x[4]+x[5]))*x[0]*1+(1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+x[8]*x[5]*x[0]*1+x[8]*x[3]*x[0]*1+0);
        b[1] = fq[1] - x[1];
        b[2] = fq[2] - x[2];
        b[3] = fq[3] - x[3]+(x[8]*x[3]*x[0]*1+0);
        b[4] = fq[4] - x[4];
        b[5] = fq[5] - x[5]+((1 - (x[6]+x[7]+x[8]))*x[5]*x[0]*1+x[8]*x[5]*x[0]*1+0);
        b[6] = fq[6] - x[6];
        b[7] = fq[7] - x[7];
        b[8] = fq[8] - x[8]+(x[8]*(1 - (x[3]+x[4]+x[5]))*x[0]*1+x[8]*x[5]*x[0]*1+x[8]*x[3]*x[0]*1+0);
        
#if defined BLAS || defined BLAS_OMP
        
        int  *ipiv;
        ipiv=(int *)mCalloc(n, sizeof(int));
        
        dgesv_(&n, &rhs, J, &n, ipiv, b, &n, &info);
        
        free(ipiv);
        
#else
        
        info=MyGaussElimination_gNxRHS(J,b,n,rhs);
        
#endif
        
        if(!info) 
          For(i,n) x[i]-=b[i];
        else 
          Warn_And_Exit("Error in Gauss Elimination Scheme.\n");
        
      }
      
      break;
    }
    default: Warn_And_Exit("Genetic code not implemented."); break;
  }
  
  For(i,n) if(x[i]<0.00001) x[i]=0.00001;
  freq[0]=x[0];
  freq[1]=x[1];
  freq[2]=x[2];
  freq[3]=1.0-(x[0]+x[1]+x[2]);
  freq[4]=x[3];
  freq[5]=x[4];
  freq[6]=x[5];
  freq[7]=1.0-(x[3]+x[4]+x[5]);
  freq[8]=x[6];
  freq[9]=x[7];
  freq[10]=x[8];
  freq[11]=1.0-(x[6]+x[7]+x[8]);
  
  free(x);free(fq);free(J);free(b);
}
/*********************************************************/
phydbl my2norm(phydbl *x, int len) //!< Added by Marcelo.
{
  phydbl sum;
  sum=0.0;
  int i;
  
  For(i,len) sum+=x[i]*x[i];
  sum=(phydbl)sqrt((double)sum);
  return sum;
}
/*********************************************************/
int MyGaussElimination_gNxRHS(double *A, double *B, int n, int rhs) //!< Added by Marcelo.
{
  int i,j,k,l,maxrow;
  double tmp;
  
  for(i=0;i<n-1;i++) {
    
    /* Find the row with the largest first value */
    maxrow = i;
    for (k=i+1;k<n;k++) if (fabs(A[k*n + i]) > fabs(A[maxrow*n + i])) maxrow = k;
    
    if(maxrow!=i)
    {
      /* Swap the maxrow and ith row */
      for (j=i;j<n;j++) 
      {
        tmp = A[i*n+j];
        A[i*n+j] = A[maxrow*n+j];
        A[maxrow*n+j] = tmp;
      }
      for (j=0;j<rhs;j++) 
      {
        tmp = B[i*rhs+j];
        B[i*rhs+j] = B[maxrow*rhs+j];
        B[maxrow*rhs+j] = tmp; 
      }
    }
    
    /* Singular matrix? */
    if (ABS(A[i*n+i]) <EPS) return 1;
    
    
    /* Eliminate the jth element of the kth row */
    for (k=i+1;k<n;k++) 
    {
      tmp=A[k*n+i]/A[i*n+i];
      for (j=i;j<n;j++) 
      {
        A[k*n+j] -= tmp*A[i*n+j];
      }
      for (j=0;j<rhs;j++) 
      {
        B[k*rhs+j] -= tmp*B[i*rhs+j];
      }
    }
  }
  
  /* Do the back substitution */
  for (k=n-1;k>=0;k--) 
  {
    for(l=0;l<rhs;l++)
    {
      tmp = 0.0;
      for (j=k+1;j<n;j++) 
	    {
	      tmp += A[k*n+j] * B[j*rhs+l];
	    }
      B[k*rhs+l] = (B[k*rhs+l] - tmp) / A[k*n+k];
    }
  }
  return 0;
}
/*********************************************************/
void Sprint_codon(char *inout,int codon) //!< Added by Marcelo.
{
  int i, remainder=0;
  char nt[4]="TCAG";
  For(i,3){
    remainder=codon-((codon>>2)<<2); //!<  n<<k= n*2^k, n>>k=n/2^k.
    codon=codon>>2;
    inout[2-i]=nt[remainder];
  }
  inout[3]=0;
}
/*********************************************************/

void Make_Ts_and_Tv_Matrix(option *io)//!< Added by Marcelo. Store the number of transitions transversions and if synonymous or nonsynonymous substitution for the given genetic code.
{
  int i,j, numSensecodons, k, icodon[3], jcodon[3], nts, ntv, diff, codoni, codonj,sSub, caseK07;
  ts_and_tv *mat;
  
  mat=io->structTs_and_Tv;
  
  numSensecodons=io->mod->ns;
  
  For(i,numSensecodons)
  {
    For(j,numSensecodons)
    {	
      diff=0;
      nts=0;
      ntv=0;
      sSub=1;
      caseK07=-1;
      if(i==j)//!< It is a diagonal element ... jump to the next. 
      {
        mat[numSensecodons*i+j].ndiff=diff;
        mat[numSensecodons*i+j].nts=nts;
        mat[numSensecodons*i+j].ntv=ntv;
        mat[numSensecodons*i+j].sSub=sSub;
        mat[numSensecodons*i+j].caseK07=caseK07;
      }
      else
      {
        
        codoni=senseCodons[i];
        codonj=senseCodons[j];
        
        if(aminoAcidmap[codoni]!=aminoAcidmap[codonj]) sSub=0;   //!< Synonymous substitution?
        
        For(k,3)
        {
          icodon[k]=codoni-((codoni>>2)<<2);                     //!<  n<<k= n*2^k, n>>k=n/2^k.
          codoni=codoni>>2;
          jcodon[k]=codonj-((codonj>>2)<<2);
          codonj=codonj>>2;
          
          //!< Transcribe the codon into the constituent nucleotides-
          if(icodon[k]!=jcodon[k]) 
          {
            diff++;
            
            switch(abs(icodon[k]*icodon[k]-jcodon[k]*jcodon[k])) /*!< The difference of the square of the values yields */ 
            {                                                    /*!< a nice way to differ between transitions and transversions.*/
              case 1:
              case 5: nts++;  //!< Transition.
                break;
              case 3:
              case 4:
              case 8:
              case 9: ntv++; //!< Transversion.
                break;	
              default: Warn_And_Exit("Cannot assign nucleotide substitution.\n");break;
            }
          }
        }
        if(diff==1)/*!< cases 0-8 according to Kosiol 2007.*/
        {
          if(nts) caseK07=0; else caseK07=1;
        }
        else if(diff==2)
        {
          if(nts==2) caseK07=2; else if(nts==1) caseK07=3; else caseK07=4;
        }
        else if(diff==3)
        {
          if(nts==3) caseK07=5; else if(nts==2) caseK07=6; else if(nts==1) caseK07=7; else caseK07=8;
        }
        mat[numSensecodons*i+j].ndiff=diff;
        mat[numSensecodons*i+j].nts=nts;
        mat[numSensecodons*i+j].ntv=ntv;
        mat[numSensecodons*i+j].sSub=sSub;
        mat[numSensecodons*i+j].caseK07=caseK07;
      }
    }
  }
}
/*********************************************************/
void Scale_freqs(phydbl *f, int n)
{
  phydbl sum;
  int i;
  do{
    sum = 0.0;
    For(i,n)
    {
      if(f[i] < 1e-6) f[i]=1e-6;
      else if(f[i] > 0.999) f[i]=0.999;
      sum += f[i];
    }
    For(i,n) f[i]/=sum;
  }
  while((sum > 1.01) || (sum < 0.99));
}
/*********************************************************/
void Scale_freqs_tol(phydbl *f, int n, phydbl tolmin, phydbl tolmax)
{
  phydbl sum;
  int i;
  do{
    sum = 0.0;
    For(i,n)
    {
      if(f[i] < tolmin) f[i]=tolmin;
      else if(f[i] > tolmax) f[i]=tolmax;
      sum += f[i];
    }
    For(i,n) f[i]/=sum;
  }
  while((sum > 1.0001) || (sum < 0.99));
}
/*********************************************************/
void Freq_to_UnsFreq(phydbl *f, phydbl *uf, int n, int f2U)
{
  phydbl sum;
  int i;
  if(f2U)
  {
    sum=1-Sum(f,n-1);
    sum=1/sum;
    For(i,n-1)  uf[i]=log(f[i]*sum);
    uf[n-1]=0;
  }
  else
  {
    for(i=0,sum=1; i<n-1; i++)  sum+=(f[i]=exp(uf[i]));
    For(i,n-1)  f[i]/=sum;
    f[n-1]=1/sum;
  }
}
/*********************************************************/
void lkExperiment(t_tree * tree, int num_tree)
{
  phydbl minW=tree->io->minParam, minK=tree->io->minParam, maxK=tree->io->maxParam, stepSize=tree->io->lkExpStepSize, kappa, omega;
  // phydbl maxW=tree->io->maxParam;
  int gridSize=maxK/stepSize, i, j;
  lkexperiment *mat;
  char nTree[10], z[20];
  FILE *outW, *outK, *outLK;
  
  mat=(lkexperiment *)calloc(gridSize*gridSize,sizeof(lkexperiment));
  
  omega=minW;
  
  tree->mod->update_eigen = YES;
  
  sprintf(nTree, "_%d", num_tree);
  
  if(tree->io->heuristicExpm)
  {
    tree->io->expm=TAYLOR;
    strcpy(z,"WexpTaylor\0");
    strcat(z,nTree);
    outW  = fopen(z, "w");
    strcpy(z,"KexpTaylor\0");
    strcat(z,nTree);
    outK  = fopen(z, "w");
    strcpy(z,"LkexpTaylor\0");
    strcat(z,nTree);
    outLK = fopen(z,"w");
  }
  else
  {
    strcpy(z,"WexpEigen\0");
    strcat(z,nTree);
    outW  = fopen(z, "w");
    strcpy(z,"KexpEigen\0");
    strcat(z,nTree);
    outK  = fopen(z, "w");
    strcpy(z,"LkexpEigen\0");
    strcat(z,nTree);
    outLK = fopen(z,"w");
  }
  
  For(i,gridSize)
  {
    omega += ((i&&1)*stepSize);
    kappa=minK;
    For(j,gridSize)
    {
      printf("i %d j %d\n",i,j);
      kappa += ((j&&1)*stepSize);
      //mat->omega       = omega;
      mat->kappa       = kappa;
      tree->mod->kappa = mat->kappa;
      //tree->mod->omega = mat->omega;

      int omegai; //added by Ken 17/8/2016
      for(omegai=0;omegai<tree->mod->nomega_part;omegai++){
      	mat->omega_part[omegai]    = tree->mod->omega_part[omegai];
        fprintf(outW,"%f ", mat->omega_part[omegai]);
      }
      mat->nomega_part = tree->mod->nomega_part;


      mat->lk          = Lk(tree);
      fprintf(outK,"%f ", mat->kappa);
      fprintf(outLK,"%.30f ",mat->lk);
    }
    fprintf(outW,"\n");
    fprintf(outK,"\n");
    fprintf(outLK,"\n");
  }  
  fprintf(outW,"\n");
  fprintf(outK,"\n");
  fprintf(outLK,"\n");
  
  if(tree->io->heuristicExpm)
  {
    tree->io->expm=EIGEN;
    Lk(tree);
  }
  
  fclose(outW);
  fclose(outK);
  fclose(outLK);
  free(mat);
}
/*******************************************************/
static unsigned int z_rndu=1237;
static int          w_rndu=1237;

void SetSeed (unsigned int seed)
{
  if(sizeof(int) != 4) 
    Warn_And_Exit("oh-oh, we are in trouble.  int not 32-bit?");
  z_rndu = seed;
  //w_rndu = abs(seed);
  w_rndu=seed; //corrected by Ken 12/1, because abs(unsigned int) doesn't make sense
}
/********************************************************/
double rndu (void)
{
  /* U(0,1): AS 183: Appl. Stat. 31:188-190 
   Wichmann BA & Hill ID.  1982.  An efficient and portable
   pseudo-random number generator.  Appl. Stat. 31:188-190
   
   x, y, z are any numbers in the range 1-30000.  Integer operation up
   to 30323 required.
   */
  static unsigned int x_rndu=11, y_rndu=23;
  double r;
  
  x_rndu = 171*(x_rndu%177) -  2*(x_rndu/177);
  y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176);
  z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178);
  /*
   if (x_rndu<0) x_rndu += 30269;
   if (y_rndu<0) y_rndu += 30307;
   if (z_rndu<0) z_rndu += 30323;
   */
  r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0;
  return (r-(int)r);
}
/********************************************************/
int myFactorial(int n) //! Added by Marcelo.
{
  if(n<=1)
    return 1;
  else return n*myFactorial(n-1);
}
/********************************************************/

FILE * openOutputFile( char * outFile, char * appendString, char * suffix, option *io ) {
    strcpy( outFile, io->in_align_file );
    strcat( outFile, appendString );
    if( io->append_run_ID ) {
        strcat( outFile, "_");
        strcat( outFile, io->run_id_string );
    }
    strcat( outFile, suffix );
    FILE * outPtr = Openfile( outFile, 1 );
    return(outPtr);
}


// added by stefan for output format support

token *Emit_Root_Token() {
    token *tkn = (token *)malloc(sizeof(token));
    tkn->type = ROOTTKN;
    tkn->val = mCalloc(1, sizeof(char));
    tkn->full = mCalloc(1, sizeof(char));
    tkn->tag = mCalloc(1, sizeof(char));
    tkn->next = NULL;
    return tkn;
}

token *Emit_Out_Token(token *currTkn, char *full, char *tag, int type, int f, char *format, ...) {
    va_list ptr;
    char *val;
    phydbl dval;
    dval = 0.0;
    
#ifdef MPI
    if(Global_myRank == 0) {
        va_start(ptr, format);
        //modified by Ken 22/7/2016
        if(f == TFNUMERIC || f == TFNUMERICNEQ) {
            dval = va_arg(ptr, double);
        } else {
            int temp = vasprintf(&val,format, ptr);//changed by Ken 12/1/2017
        }
        va_end(ptr);
    }
#else
    va_start(ptr, format);
    //modified by Ken 22/7/2016
    if(f == TFNUMERIC || f == TFNUMERICNEQ) {
        dval = va_arg(ptr, double);
        val = mCalloc(1, sizeof(char));
        strcpy(val, "");
    } else {
    	int temp = vasprintf(&val,format, ptr);//changed by Ken 12/1/2017
    }
    va_end(ptr);
#endif
    fflush (NULL);
    
    token *tkn = (token *)malloc(sizeof(token));
    tkn->full = mCalloc((int)strlen(full)+1, sizeof(char));
    strcpy(tkn->full, full);
    tkn->tag = mCalloc((int)strlen(tag)+1, sizeof(char));
    strcpy(tkn->tag, tag);
    tkn->val = val;
    tkn->dval = dval;
    tkn->type = type;
    tkn->format = f;
    tkn->next = NULL;
    currTkn->next = tkn;
    currTkn = tkn;
    return(currTkn);
}

token * Del_Last_Token(token *r) {
    while(1) {
        if(r->next && r->next->next) {
            r = r->next;
        } else {
            if(r->next) {
                Free_OutToken(r->next);
            }
            r->next = NULL;
            break;
        }
    }
    return(r);
}

void Output_Tokens(token *root, int format, FILE *fp_out) {
    switch(format) {
#ifdef USEYAML
        case OUTYAML: {
            Print_YAML_Tokenlist(fp_out, root);
            break;
        }
#endif
        case OUTDARWIN: {
            Print_Darwin_Tokenlist(fp_out, root);
            break;
        }
        default: {
            Print_Txt_Tokenlist(fp_out, root);
            break;
        }
    }
    Free_OutList(root);
}

void Print_Txt_Tokenlist(FILE *fp_out, token *t) {
    int level = 0;
    int tabsize = 4;
    int isTable = NO;
    char *tablespace = "   ";
    char * tabspace = malloc(256*sizeof(char));
    const char *padding = "\t\t\t\t\t\t\t\t\t\t\t";
    int spacerlen = 0;
    char *spacer = malloc(128*sizeof(char));
    //sprintf(spacer, "");
    //sprintf(tabspace, "");
    spacer[0]= 0; //fixed by Ken 12/1
    tabspace[0]= 0;
    while(1) {
        switch(t->type) {
            case SCALARTKN: {
                if(isTable) {
                    switch(t->format) {
                        case TFNUMERIC: {
                            PhyML_Fprintf(fp_out, "%s=%.8f%s", t->full, t->dval, tablespace);
                            break;
                        }
                        case TFNUMERICNEQ: {//modified by Ken 22/7/2016
                              PhyML_Fprintf(fp_out, "%s%.8f%s", t->full, t->dval, tablespace);
                              break;
                        }
                        case TFSTRING:
                        default: {
                            PhyML_Fprintf(fp_out, "%s=%s%s", t->full, t->val, tablespace);
                            break;
                        }
                    }
                } else {
                    spacerlen = 5 - (int)(strlen(t->full) / 6);
                    sprintf(spacer, "%s: %*.*s", t->full, spacerlen, spacerlen, padding);
                    switch(t->format) {
                        case TFNUMERIC: {
                            PhyML_Fprintf(fp_out, "\n%s. %s%g", tabspace, spacer, t->dval);
                            break;
                        }
                        case TFNUMERICNEQ: {//modified by Ken 22/7/2016
                            PhyML_Fprintf(fp_out, "%s%.8f%s", t->full, t->dval, tablespace);
                            break;
                        }
                        case TFSTRING:
                        default: {
                            PhyML_Fprintf(fp_out, "\n%s. %s%s", tabspace, spacer, t->val);
                            break;
                        }
                    }
                }
                break;
            }
            case SETSTARTTKN: {
                PhyML_Fprintf(fp_out, "\n%s. %s: %s", tabspace, t->full, t->val);
                level += 1;
                sprintf(tabspace, "%*s", level*tabsize, "");
                break;
            }
            case SETENDTKN: {
                level -= 1;
                sprintf(tabspace, "%*s", level*tabsize, "");
                break;
            }
            case TBLSTARTTKN: {
                if(strcmp(t->full, "")) {
                    PhyML_Fprintf(fp_out, "\n%s%s: ", tabspace, t->full);
                } else {
                    PhyML_Fprintf(fp_out, "\n%s", tabspace, t->full);
                }
                isTable = YES;
                break;
            }
            case TBLENDTKN: {
                isTable = NO;
                break;
            }
            case COMMENTTKN: {
                PhyML_Fprintf(fp_out, "\n# %s", t->val);
                break;
            }
        }
        if(t->next != NULL) {
            t = t->next;
        } else {
            break;
        }
    }
    free(tabspace);
}

#ifdef USEYAML
void Print_YAML_Tokenlist(FILE *fp_out, token *t) {
    int level = 0;
    int tabsize = 4;
    char * tabspace = malloc(256*sizeof(char));
    sprintf(tabspace, "");
    while(1) {
        switch(t->type) {
            case SCALARTKN: {
                switch(t->format) {
                    case TFNUMERIC: {
                        PhyML_Fprintf(fp_out, "\n%s%s: %g", tabspace, t->tag, t->dval);
                        break;
                    }
                    case TFSTRING:
                    default: {
                        PhyML_Fprintf(fp_out, "\n%s%s: %s", tabspace, t->tag, t->val);
                        break;
                    }
                }
                break;
            }
            case SETSTARTTKN: {
                PhyML_Fprintf(fp_out, "\n%s%s%s:", tabspace, t->tag, t->val);
                level += 1;
                sprintf(tabspace, "%*s", level*tabsize, "");
                break;
            }
            case SETENDTKN: {
                level -= 1;
                sprintf(tabspace, "%*s", level*tabsize, "");
                break;
            }
            case COMMENTTKN: {
                PhyML_Fprintf(fp_out, "\n# %s", t->val);
                break;
            }
        }
        if(t->next != NULL) {
            t = t->next;
        } else {
            break;
        }
    }
    free(tabspace);
}
#endif

void Print_Darwin_Tokenlist(FILE *fp_out, token *t) {
    PhyML_Fprintf(fp_out, "if not assigned(cpres) then\n");
    PhyML_Fprintf(fp_out, "cpres := table([],[]);\n");
    PhyML_Fprintf(fp_out, "fi;\n");
    double DARWIN_DBL_MAX = 1.7976931348623157e+308;
    double DARWIN_DBL_EPSILON = 2.2204e-16;
    token *dvar = Emit_Root_Token();
    token *dvarc = Emit_Out_Token(dvar, "cpres", "", SCALARTKN, TFSTRING, "%s", "cpres");
    while(1) {
        switch(t->type) {
            case SCALARTKN: {
                switch(t->format) {
                    case TFNUMERIC: {
                        if(fabs(t->dval) >= DARWIN_DBL_MAX) {
                            PhyML_Fprintf(fp_out, "\n%s['%s'] := %e;", getDarwinTable(dvar), t->tag, DARWIN_DBL_MAX);
                        } else if(fabs(t->dval) <= DARWIN_DBL_EPSILON) {
                            PhyML_Fprintf(fp_out, "\n%s['%s'] := %e;", getDarwinTable(dvar), t->tag, 0.0);
                        } else {
                            if(isnan(t->dval)) {
                                PhyML_Fprintf(fp_out, "\n%s['%s'] := %g;", getDarwinTable(dvar), t->tag, 0.0);
                            } else {
                                PhyML_Fprintf(fp_out, "\n%s['%s'] := %g;", getDarwinTable(dvar), t->tag, t->dval);
                            }
                        }
                        break;
                    }
                    case TFSTRING:
                    default: {
                        PhyML_Fprintf(fp_out, "\n%s['%s'] := '%s';", getDarwinTable(dvar), t->tag, t->val);
                        break;
                    }
                }
                break;
            }
            case SETSTARTTKN: {
                if(t->format == TFSTRING) {
                    PhyML_Fprintf(fp_out, "\n%s['%s'] := table();", getDarwinTable(dvar), t->tag);
                    dvarc = Emit_Out_Token(dvarc, "\ncpres", "", SCALARTKN, TFSTRING, "['%s%s']", t->tag, t->val);
                } else {
                    PhyML_Fprintf(fp_out, "\n%s[%s] := table();", getDarwinTable(dvar), t->tag);
                    dvarc = Emit_Out_Token(dvarc, "\ncpres", "", SCALARTKN, TFSTRING, "[%s]", t->tag);
                }
                break;
            }
            case SETENDTKN: {
                dvarc = Del_Last_Token(dvar);
                break;
            }
                
            case TBLSTARTTKN:
            case TBLENDTKN: {
                break;
            }
                
            case COMMENTTKN: {
                PhyML_Fprintf(fp_out, "\n# %s", t->val);
                break;
            }
        }
        if(t->next != NULL) {
            t = t->next;
        } else {
            break;
        }
    }
}

char * getDarwinTable(token *r) {
    int len = 0;
    token *tmp = r;
    while(1) {
        len += strlen(tmp->val);
        if(tmp->next != NULL) {
            tmp = tmp->next;
        } else {
            break;
        }
    }
    char *res = malloc((len+1)*sizeof(char));
    strcpy(res, "");
    while(1) {
        strcat(res, r->val);
        if(r->next != NULL) {
            r = r->next;
        } else {
            break;
        }
    }
    return(res);
}

