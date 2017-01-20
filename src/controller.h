
/*
 
 controller.h
 codonPhyML
 
 Created by Stefan Zoller on 3/30/12.
 Copyright (c) 2012 CBRG. All rights reserved.
 
 
 All parts of  the source except where indicated  are distributed under
 the GNU public licence.  See http://www.opensource.org for details.
 
 */

#ifndef codonPhyML_controller_h
#define codonPhyML_controller_h


#include "utilities.h"


void finishOptions(option * io);

int mainOptionSwitch(int opt, char * optarg, option * io);

void checkForRandStartTree(option * io);
void checkModelCombinations(option * io);
void adaptForM4(option * io);
void createOutFiles(option * io);
void cleanupParameters(option * io);
void setupGeneticCode(option * io);
void setupFreqs(option * io);
void setupInitialRateMats(option * io);
void setupOmegaCats(option * io);
void setupPCM(option * io);
void setupModelIdentifier(option * io);
void setupFreqHandling(option * io);
void setupKappa(option * io);

#endif
