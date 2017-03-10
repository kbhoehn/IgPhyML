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
 
#ifndef NUCLE2CODON_H
#define NUCLE2CODON_H

#include "utilities.h"
 
/*######start definition of the quad-tree######*/ 
typedef struct __Node_nucle2codon{               /*!< Quad-Tree Node. */
  int theCodon;                                  /*!< Codons are enumerated from 0-63. */
  int stopCodons;                                 /*!< 1 in the case of stop codon. */ 
  struct __Node_nucle2codon *tP, *cP, *aP, *gP;  
}t_node_nucle2codon;

t_node_nucle2codon * Tree_Init();                /*!< Create and Initialize the Tree. */
void Delete_Tree(t_node_nucle2codon *);          /*!< Delete a full quad-tree of 4 levels. */
void Translate_nucle2codons(t_node_nucle2codon *,char *,int); /*!< Translate triples of nucleotides {T,C,A,G} into the respective codon [0,63] */
void Branch_intheTree(t_node_nucle2codon *,char);/*!< 4 options T,C,A,G ... follow the pointer accordinly. */
/*######end definition of the quad-tree######*/

/*######start definition of the stack to hold nodes of the quad-tree######*/ 
typedef struct __Node_stackUnknown{              /*!< Build and Depth-First Search for Quad-Tree uses a stack. */     
  t_node_nucle2codon *key;        
  struct __Node_stackUnknown *nextCodon;    
}t_node_stackUnknown; 

void Node_stackInit();                           /*!< Initialize the Stack. */
void Push_nodeStack( t_node_nucle2codon*);       /*!< Push. */
t_node_nucle2codon * Pop_nodeStack();            /*!< Pop. */
int Node_stackEmpty();                           /*!< Is stack empty?. */
void Free_nodeStack();                           /*!< Release memory allocated for head and tail pointers. */
/*######end definition of the stack to hold nodes of the quad-tree######*/ 

/*######start definition of the l-list to hold the resulting codons and alternatives in the case of gaps/unknowns######*/
typedef struct __Node_listCodons{                /*!< The resulting Codons and alternatives are grouped in a Linked-list. */       
  int theCodon;
  int moreCombinations;                          /*!< In a ambigous site there can be many combinations.*/ 
  struct __Node_listCodons *nextSite;
  struct __Node_listCodons *previousSite; 
  struct __Node_listCodons *nextCodon;
}t_node_listCodons;

void Node_listInit();                            /*!< Initialize the Linked-list. */
void Free_listCodons();                          /*!< Release memory allocated for head and tail pointers. */
void Insert_siteCodon(int);                      /*!< Insert a codon [0-63] at the site or indicate 64. */
void Insert_siteAlternative(int);                /*!< and start adding the alternatives that will be included in. */
                                                 /*!< the sum for the total likelihood of site i. */
int Read_nextSite();                             /*!< Return the next codon in the sequence. */
void MarkasAmbigousSite();                       /*!< Mark the first node as ambigous in order to retrieve a list of possibilities later.*/
t_node_listCodons * Read_listAlternative();      /*!< Return the pointer to the list of alternative codon. */
int Read_nextAlternative(t_node_listCodons *);   /*!< Return the next codon in the list of alternatives. */
int Site_listEmpty();
/*######end definition of the l-list to hold the resulting codons and alternatives in the case of gaps/unknowns######*/

/*######start implementation of the interface of nucleotide to codon translation######*/
void Nucleotides2Codons(align **, option *);
char **CompactCodonAlternatives(int lenght, align *data, char **alternativeCodons);
/*######end implementation of the interface of nucleotide to codon translation######*/

#endif
     
  
