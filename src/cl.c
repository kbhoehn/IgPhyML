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
#include "cl.h"

extern int stopCodons[64];
extern int senseCodons[64];
extern phydbl pcafreq[61];

extern struct option longopts[];

/*********************************************************/
/**
 * Fill the Option fields, with the argc array
 */ 
void Read_Command_Line( option *io, int argc, char **argv )
{
    int c;
    
    int switchResult;

    io->mod->motifstringopt=0;
    io->mod->hotnessstringopt=0;
    io->mod->partfilespec=0;
    io->mod->rootfound=0;
    io->mod->partfile="NONE";
    io->mod->ambigprint=0;
    io->mod->nomega_part=1;
    io->mod->nparts=1;
    io->mod->ambigprint=0;
    io->mod->startnode=0;
	io->mod->slowSPR=0;
	io->mod->stretch=1.0;


    io->mod->rootname = mCalloc(T_MAX_OPTION,sizeof(char));
    io->mod->hotnessstring = mCalloc(T_MAX_OPTION,sizeof(char));
    io->mod->aamodel = mCalloc(T_MAX_OPTION,sizeof(char));
    io->mod->partfile = mCalloc(T_MAX_FILE,sizeof(char));
    io->mod->motifstring = mCalloc(T_MAX_FILE,sizeof(char));
    io->mod->ambigfile = mCalloc(T_MAX_FILE,sizeof(char));


    while((c = getopt_long_only(argc,argv,"qi:d:g:m:b:n:w:t:f:v:c:a:u:ho:s:p",longopts,NULL)) != -1)            //! Removed zk:x:l:e since they are not implemented  (Louis)
    {
        if(c == 139 || c == 143) {
            // config file
#ifdef USEYAML
            if(c == 139) {
                io->confformat = OUTYAML;
            } else
#endif
            if(c == 143) {
                io->confformat = OUTDARWIN;
            }
            char *tmp;
            int result = -1;
            tmp = (char *) mCalloc (T_MAX_FILE, sizeof(char));
            if( strlen (optarg) > T_MAX_FILE -16 ) {
                strcpy( tmp, "\n. The file name'" );
                strcat( tmp, optarg );
                strcat( tmp, "' is too long.\n" );
                Warn_And_Exit( tmp );
            } else if( !Filexists( optarg ) ) {
                strcpy( tmp, "\n. The file '" );
                strcat( tmp, optarg );
                strcat( tmp, "' does not exist.\n" );
                Warn_And_Exit( tmp );
            } else {
            	printf("\n\nIgPhyML Doesn't take config files right now.\n\n");
            	exit(EXIT_FAILURE);
                //result = handleConfigFile( optarg, io );
            }
            free( tmp );
            if( result != 0 ) {
                Warn_And_Exit( "There was an error when parsing the config file.\n" );
            }
        } else {
            // command line interface
            switchResult = mainOptionSwitch( c, optarg, io );
            
            if( c == '?' ) {
                if( isprint( optopt ) ) {
                    PhyML_Printf("\n. Unknown option `-%c'.\n", optopt );
                } else {
                    PhyML_Printf("\n. Unknown option character `\\x%x'.\n", optopt );
                }
                Warn_And_Exit( "" );
            }
            
            if( switchResult != 0 ) {
                Usage();
            }
        }
    }
    
    /*   if((io->mod->whichmodel == K80) || (io->mod->whichmodel == JC69)) */
    /*     { */
    /*       if(io->mod->s_opt->opt_state_freq) */
    /* 	{ */
    /* 	  char c; */
    /* 	  PhyML_Printf("\n. WARNING: nucleotide frequencies must be set to 1/4 with this model.\n"); */
    /* 	  PhyML_Printf("\n. Type the enter key to resume the analysis.\n"); */
    /* 	  scanf("%c",&c); */
    /* 	} */
    /*       io->mod->s_opt->opt_state_freq = 0; */
    /*     } */
    
    
#ifndef PHYML
    if( (io->open_ps_file) || (io->m4_model == YES) ) {
        strcpy( io->out_ps_file, io->in_align_file );
        strcat( io->out_ps_file, "_mc_tree.ps" );
        io->fp_out_ps = Openfile( io->out_ps_file, 1 );
    }
#endif 
    
    finishOptions( io );
    
	return;
}

/*********************************************************/
