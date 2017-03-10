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

#ifndef CL_H
#define CL_H

#include <unistd.h>
#include <getopt.h>
#include "utilities.h"
#include "help.h"
#include "models.h"
#include "free.h"
/* unused ken 5/1
#include "interface.h"
#include "configfile.h"
*/
#include "controller.h"

void Read_Command_Line(option *input, int argc, char **argv);

#endif
