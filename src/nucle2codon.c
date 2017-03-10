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

#include "nucle2codon.h"
#include <utilities.h>

extern int     stopCodons[64];
extern int   senseCodons[64];
/*######start implementation of the quad-tree######*/ 
t_node_nucle2codon * Tree_Init(){
  t_node_nucle2codon *rootTreeCodon, *pointerTreeCodon;
  int depthTree,indexCodon;
  // int i;
  Node_stackInit();
  rootTreeCodon=(t_node_nucle2codon *)malloc(sizeof(t_node_nucle2codon));
  Push_nodeStack(rootTreeCodon);
  depthTree=-1;
  indexCodon=0;
  while(!Node_stackEmpty()){
    depthTree++; 
    pointerTreeCodon=Pop_nodeStack();
    pointerTreeCodon->theCodon=-depthTree;                                                /*!< Property of internal nodes*/
    pointerTreeCodon->stopCodons=0;                                                        /*! which happens to be very usefull*/
    pointerTreeCodon->tP=(t_node_nucle2codon *)malloc(sizeof(t_node_nucle2codon));        /*! during the translation nucle2codon.*/
    pointerTreeCodon->cP=(t_node_nucle2codon *)malloc(sizeof(t_node_nucle2codon));
    pointerTreeCodon->aP=(t_node_nucle2codon *)malloc(sizeof(t_node_nucle2codon));
    pointerTreeCodon->gP=(t_node_nucle2codon *)malloc(sizeof(t_node_nucle2codon));
    if (depthTree<2){
      Push_nodeStack(pointerTreeCodon->gP);
      Push_nodeStack(pointerTreeCodon->aP); 
      Push_nodeStack(pointerTreeCodon->cP);
      Push_nodeStack(pointerTreeCodon->tP);
    }else{ 
      pointerTreeCodon->tP->stopCodons=stopCodons[indexCodon];
      pointerTreeCodon->tP->theCodon=indexCodon++;
      pointerTreeCodon->cP->stopCodons=stopCodons[indexCodon];
      pointerTreeCodon->cP->theCodon=indexCodon++;
      pointerTreeCodon->aP->stopCodons=stopCodons[indexCodon];
      pointerTreeCodon->aP->theCodon=indexCodon++;
      pointerTreeCodon->gP->stopCodons=stopCodons[indexCodon];
      pointerTreeCodon->gP->theCodon=indexCodon++;
      if((indexCodon%16==0)&&(indexCodon!=0)){depthTree-=2;}else{depthTree--;}
    } 
  }
  Free_nodeStack();
  return rootTreeCodon;
}

void Delete_Tree(t_node_nucle2codon *root){
  t_node_nucle2codon *pointerTreeCodon;
  int depthTree,indexCodon;
  // int i;
  Node_stackInit();
  Push_nodeStack(root);
  depthTree=-1;
  indexCodon=0;
  while(!Node_stackEmpty()){
    depthTree++;
    pointerTreeCodon=Pop_nodeStack();
    if (depthTree<2){
      Push_nodeStack(pointerTreeCodon->gP);
      Push_nodeStack(pointerTreeCodon->aP);
      Push_nodeStack(pointerTreeCodon->cP);
      Push_nodeStack(pointerTreeCodon->tP);
      free(pointerTreeCodon);
    }else{
      free(pointerTreeCodon->tP);
      free(pointerTreeCodon->cP);
      free(pointerTreeCodon->aP);
      free(pointerTreeCodon->gP);
      free(pointerTreeCodon);
      indexCodon+=4;
      if((indexCodon%16==0)&&(indexCodon!=0)){depthTree-=2;}else{depthTree--;}
    }
  }
  Free_nodeStack();
}

  void Translate_nucle2codons(t_node_nucle2codon *root,char *sequence,int s_length){
  int i,counter;
  t_node_nucle2codon *pointer;
  Node_stackInit();                                     /*!< Initialize the stack.*/
                                                        /*! which means allocate memory for respective pointers (head,tail)*/
  for(i=0;i<s_length-2;i+=3)                            /*! internal nodes in the tree are <0 , leafs are in [0,63].*/
  {
    pointer=root;
    Branch_intheTree(pointer,sequence[i]);
    counter=0;
    while(!Node_stackEmpty())
    {
      pointer=Pop_nodeStack();
      if(pointer->theCodon<0)
      {
	Branch_intheTree(pointer,sequence[i-(pointer->theCodon)]);
      }
      else
      {
	if(!(pointer->stopCodons))
	{
	  if(++counter==1)
	  {
	    Insert_siteCodon(pointer->theCodon);
	  }else
	  { 
	    if(counter==2) MarkasAmbigousSite();
	    Insert_siteAlternative(pointer->theCodon);
	  }
	}
      }
    }
  }
  Free_nodeStack();
}

void Branch_intheTree(t_node_nucle2codon *ppointer,char nucleo){
  switch(nucleo){
  case 't': 
  case 'T': 
  case 'u': 
  case 'U': 
    Push_nodeStack(ppointer->tP);
    break;
    
  case 'c': 
  case 'C': 
    Push_nodeStack(ppointer->cP);
    break; 
    
  case 'a': 
  case 'A': 
    Push_nodeStack(ppointer->aP);
    break;
    
  case 'g':
  case 'G': 
    Push_nodeStack(ppointer->gP);
    break;
    
  case 'm':
  case 'M': 
    Push_nodeStack(ppointer->cP);
    Push_nodeStack(ppointer->aP); 
    break;
    
  case 'r':
  case 'R':
    Push_nodeStack(ppointer->aP);
    Push_nodeStack(ppointer->gP); 
    break;
    
  case 'w':
  case 'W':
    Push_nodeStack(ppointer->tP);
    Push_nodeStack(ppointer->aP); 
    break;
    
  case 's':
  case 'S':
    Push_nodeStack(ppointer->cP);
    Push_nodeStack(ppointer->gP); 
    break;
    
  case 'y':
  case 'Y':
    Push_nodeStack(ppointer->tP);
    Push_nodeStack(ppointer->cP); 
    break;
    
  case 'k':
  case 'K': 
    Push_nodeStack(ppointer->tP);
    Push_nodeStack(ppointer->gP); 
    break;
    
  case 'b':
  case 'B':
    Push_nodeStack(ppointer->tP);
    Push_nodeStack(ppointer->cP);
    Push_nodeStack(ppointer->gP);
    break;
    
  case 'd':
  case 'D':
    Push_nodeStack(ppointer->tP);
    Push_nodeStack(ppointer->aP);
    Push_nodeStack(ppointer->gP);
    break;
    
  case 'h':
  case 'H':
    Push_nodeStack(ppointer->tP);
    Push_nodeStack(ppointer->cP);
    Push_nodeStack(ppointer->aP);
    break;
    
  case 'v':
  case 'V':
    Push_nodeStack(ppointer->cP);
    Push_nodeStack(ppointer->aP);
    Push_nodeStack(ppointer->gP);
    break;
    
  case '-':
  case 'n':
  case 'N':
  case 'x':
  case 'X':
  case '?':
    Push_nodeStack(ppointer->tP);
    Push_nodeStack(ppointer->cP);
    Push_nodeStack(ppointer->aP);
    Push_nodeStack(ppointer->gP);
    break;
  default: printf("Forbiden character in input sequence.\n"); exit(1);break;
  }
}
/*######end implementation of the quad-tree######*/ 

/*######start implementation of the stack to hold nodes of the quad-tree######*/ 
static t_node_stackUnknown *Head_nodeStack, *Tail_nodeStack;
void Node_stackInit(){
  Head_nodeStack=(t_node_stackUnknown *)malloc(sizeof(t_node_stackUnknown));
  Tail_nodeStack=(t_node_stackUnknown *)malloc(sizeof(t_node_stackUnknown));
  Head_nodeStack->nextCodon=Tail_nodeStack;
  Head_nodeStack->key=NULL;
  Tail_nodeStack->nextCodon=Tail_nodeStack;
}

void Push_nodeStack(t_node_nucle2codon *key){
  t_node_stackUnknown *Pointer_nodeStack;
  Pointer_nodeStack=(t_node_stackUnknown *)malloc(sizeof(t_node_stackUnknown));
  Pointer_nodeStack->key=key;
  Pointer_nodeStack->nextCodon=Head_nodeStack->nextCodon;
  Head_nodeStack->nextCodon=Pointer_nodeStack;
}
 
t_node_nucle2codon * Pop_nodeStack(){
  t_node_nucle2codon *x;
  t_node_stackUnknown *Pointer_nodeStack;
  Pointer_nodeStack=Head_nodeStack->nextCodon;
  Head_nodeStack->nextCodon=Pointer_nodeStack->nextCodon;
  Pointer_nodeStack->nextCodon=Pointer_nodeStack;
  x=Pointer_nodeStack->key;
  free(Pointer_nodeStack);  
  return x;
}  
 
int Node_stackEmpty(){
  return Head_nodeStack->nextCodon==Tail_nodeStack;
}
   
void Free_nodeStack(){
  free(Head_nodeStack);
  free(Tail_nodeStack);
}
/*######end implementation of the stack to hold nodes of the quad-tree######*/ 
 
/*######start implementation of the l-list to hold the resulting codons and alternatives in the case of gaps/unknowns######*/
static t_node_listCodons *Head_listCodons, *Tail_listCodons;
void Node_listInit(){
  Head_listCodons=(t_node_listCodons *)malloc(sizeof(t_node_listCodons));
  Tail_listCodons=(t_node_listCodons *)malloc(sizeof(t_node_listCodons));
  Head_listCodons->nextSite=Tail_listCodons;
  Head_listCodons->previousSite=Head_listCodons;
  Head_listCodons->nextCodon=NULL;
  Head_listCodons->theCodon=-1;
  Tail_listCodons->nextSite=Tail_listCodons;
  Tail_listCodons->previousSite=Head_listCodons;
  Tail_listCodons->nextCodon=NULL;
  Tail_listCodons->theCodon=-1;
} 
 
void Free_listCodons(){
  free(Head_listCodons);
  free(Tail_listCodons); 
}

void Insert_siteCodon(int theCodon){
  t_node_listCodons * pointer_listCodons;
  pointer_listCodons=(t_node_listCodons * )malloc(sizeof(t_node_listCodons));
  pointer_listCodons->theCodon=theCodon;
  pointer_listCodons->moreCombinations=0;
  Tail_listCodons->previousSite->nextSite=pointer_listCodons;
  pointer_listCodons->previousSite=Tail_listCodons->previousSite;
  pointer_listCodons->nextSite=Tail_listCodons;
  Tail_listCodons->previousSite=pointer_listCodons;
  pointer_listCodons->nextCodon=NULL;
} 

void MarkasAmbigousSite()
{
  Tail_listCodons->previousSite->moreCombinations=1;
}

void Insert_siteAlternative(int theCodon){ 
  t_node_listCodons * pointer_local;
  pointer_local=(t_node_listCodons * )malloc(sizeof(t_node_listCodons));
  pointer_local->theCodon=theCodon;
  pointer_local->nextCodon=Tail_listCodons->previousSite->nextCodon;
  Tail_listCodons->previousSite->nextCodon=pointer_local;
}  
  
int Read_nextSite(){
  t_node_listCodons * pointer_local;
  int x,y;
  if(!Site_listEmpty()){
    pointer_local=Head_listCodons->nextSite;
    x=pointer_local->theCodon;
    y=pointer_local->moreCombinations;
    if(y==0){                                                                       /*!< if it is a codon [0,63] just remove it from the list*/
      Head_listCodons->nextSite=pointer_local->nextSite;                            /*! and return its value, and release the pointer.*/ 
      pointer_local->nextSite->previousSite=Head_listCodons;
      free(pointer_local);
    }else{
      if(pointer_local->nextCodon!=NULL){
	x=x+63;                                                                     /*!< If it is a indicator for unknowns or gap then we will*/
	                                                                            /*! We do nothing, just return a value > 63 and the caller */
	                                                                            /*! will process it.*/
      }else{
	Head_listCodons->nextSite=pointer_local->nextSite;                          /*! If there was just one element we delete the pointer.*/
	Head_listCodons->nextSite->previousSite=Head_listCodons;    
	free(pointer_local);  
      }
    }
    return x; 
  }
  else{
    return -1;                                                                      /*!< This means the site is a stop codon. This should not happen in real data!.*/
  }
}
   
t_node_listCodons * Read_listAlternative(){
  t_node_listCodons * pointer_local;
  if(!Site_listEmpty()){
    pointer_local=Head_listCodons->nextSite;
    Head_listCodons->nextSite=pointer_local->nextSite;
    pointer_local->nextSite->previousSite=pointer_local->previousSite;
    pointer_local->nextSite=pointer_local;
    pointer_local->previousSite=pointer_local;
    return pointer_local;
  }
  else{
    return NULL;
  }
}

int Read_nextAlternative(t_node_listCodons *head_list){
  t_node_listCodons * pointer_local;
  int x;
  x=head_list->theCodon;
  if(head_list->nextCodon==NULL){
    free(head_list);
  }else{
    pointer_local=head_list->nextCodon;
    head_list->theCodon=pointer_local->theCodon;
    head_list->nextCodon=pointer_local->nextCodon;
    free(pointer_local); 
  }
  return x;
}  

int Site_listEmpty(){
  return Head_listCodons->nextSite==Tail_listCodons;
}
/*######end implementation of the l-list to hold the resulting codons and alternatives in the case of gaps/unknowns######*/
 
/*######start implementation of the interface of nucleotide to codon translation######*/  
//!< Translate nucleotides into codons.
void Nucleotides2Codons(align **data, option *io)
{
  t_node_nucle2codon * root;
  t_node_listCodons * pointer;
  int i=0,j=0,k=0,condition=0,codon=0;
  // int endStopCodons=0;
  char **alternativeCodons;
  int seqLen;
  
  seqLen=data[0]->len;
  if(seqLen%3)
  {
    Warn_And_Exit("Nucleotide sequence length must be a multiple of 3. Check your input files.\n");
  }
  
  root=Tree_Init();                                                                         /*!< Initialize the tree.*/
  Node_listInit();                                                                          /*!< Initialize the linked-list of codons.*/
  
  alternativeCodons=(char **)mCalloc((io->mod->ns)+1,sizeof(char *));                         //!< Create one matrix that holds the alternative codons for all sites for all taxa.
  if(!alternativeCodons)                                                                    //!< +1 for the marquer at the end of each list of alternatives.
  {
    printf("Impossible to translate codons for taxon %s\n",data[i]->name);
    Warn_And_Exit("Not enough memory");
  }
  For(j,(io->mod->ns)+1)
  {
    alternativeCodons[j]=(char *)mCalloc(seqLen/3,sizeof(char));                          //!< Size is #taxa X #senseCodons.
    if(!alternativeCodons[j])
    {
      printf("Impossible to translate codons for taxon %s\n",data[i]->name);
      Warn_And_Exit("Not enough memory");
    }
  }
  For(i,io->n_otu)
  {
    For(k,(io->mod->ns)+1) For(j,seqLen/3) alternativeCodons[k][j]=(char)99;                //!< Initialize the matrix. 
    
    Translate_nucle2codons(root, data[i]->state, seqLen);                                   //!< Translate a whole sequence of the current taxon.
    j=-1;                                                                                   //!< Initialize the index j.
    while(!Site_listEmpty())
    {
      codon=Read_nextSite();
      j++;
      if(codon>63)
      {
	pointer=Read_listAlternative();
	if(pointer!=NULL)
	{ 
	  condition=1;
	  k=-1;
	  while(condition)
	  {
	    k++;
	    if(pointer->nextCodon==NULL)
	    {
	      condition=0;
	    }
	    codon=Read_nextAlternative(pointer);
	    alternativeCodons[k][j]=(char)codon; 
	  }
	}
      }else
      { 
	alternativeCodons[0][j]=(char)codon; 
      }
    }
    
    if( (j+1) != (seqLen/3) ) 
    { 
      printf("Sequence of taxon %s contains stop codon(s). ", data[i]->name);
      Warn_And_Exit("Impossible to continue. Remove stop codon sites across all sequences and try again.\n");
    } 
    
    data[i]->alternativeCodons=CompactCodonAlternatives(seqLen/3, data[i], alternativeCodons);
    
  }
  For(j,(io->mod->ns)+1) free(alternativeCodons[j]);//!< Delete old matrix;
  free(alternativeCodons);
  
  Free_listCodons();
  Delete_Tree(root);
}


char ** CompactCodonAlternatives(int seqLen, align *data, char ** alternatives)                       //!< Remove memory which is not necessary to hold.
{
  int i,j,k;
  char **altTMP, *state;
  short int *is_ambigu;
    
  altTMP=(char **)mCalloc(seqLen,sizeof(char *)); 
  if(!altTMP) Warn_And_Exit("Impossible to compact alternative codons. Not enough memory");
  
  free(data->is_ambigu);
    
  data->ntStates=data->state;                                                                         //!< Save the nucleotide data for later.
  data->ntLen=data->len;
  
  
  is_ambigu=(short int *)mCalloc(seqLen, sizeof(short int));
  if(!is_ambigu) Warn_And_Exit("Impossible to compact alternative codons. Not enough memory");
  state=(char *)mCalloc(seqLen, sizeof(char));
  if(!state) Warn_And_Exit("Impossible to compact alternative codons. Not enough memory");
  
  For(j,seqLen)
  {
    i=-1;
    while(alternatives[++i][j]<64);
    
    altTMP[j]=(char *)mCalloc(i+1,sizeof(char));                                                      //!< Create one matrix that holds the alternative codons for all sites for all taxa.
    if(!altTMP[j]) Warn_And_Exit("Impossible to compact alternative codons. Not enough memory");
    
    For(k,i+1) altTMP[j][k]=alternatives[k][j];
    
    if(!(i-1))
    {
      state[j]=altTMP[j][0];
      altTMP[j][0]=(char)99;
      is_ambigu[j]=(short int)0;
    }
    else
    {
      state[j]='X';
      is_ambigu[j]=(short int)1;
    }
  }
  data->len=seqLen;
  data->state=state;
  data->is_ambigu=is_ambigu;
  
  return altTMP;
}

/*######end implementation of the interface of nucleotide to codon translation######*/  
