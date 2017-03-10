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

#ifndef ALRT_H
#define ALRT_H

#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "models.h"
#include "free.h"
#include "simu.h"


void aLRT(t_tree *tree);
int Check_NNI_Five_Branches(t_tree *tree);
int Compute_Likelihood_Ratio_Test(t_edge *tested_edge, t_tree *tree);
int NNI_Neigh_BL(t_edge *b_fcus, t_tree *tree);
void Make_Target_Swap(t_tree *tree, t_edge *b_fcus, int swaptodo);
phydbl Statistics_To_Probabilities(phydbl in);
phydbl Statistics_To_RELL(t_tree *tree);
phydbl Statistics_To_SH(t_tree *tree);
phydbl Update_Lk_At_Given_Edge_Excluding(t_edge *b_fcus, t_tree *tree, t_node *exclude);

#endif
