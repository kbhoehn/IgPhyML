
/*
 
 controller.c
 codonPhyML
 
 Created by Stefan Zoller on 3/30/12.
 Copyright (c) 2012 CBRG. All rights reserved.
 
 
 All parts of  the source except where indicated  are distributed under
 the GNU public licence.  See http://www.opensource.org for details.
 
 */

#include "controller.h"
#include <stdio.h>
#include <stdlib.h>

struct option longopts[] =
{
    {"n_rgrft",           required_argument,NULL,0},
    {"n_globl",           required_argument,NULL,1},
    {"max_dist",          required_argument,NULL,2},
    {"n_optim",           required_argument,NULL,3},
    {"n_best",            required_argument,NULL,4},
    {"model",             required_argument,NULL,5},	//! -m
    {"search",            required_argument,NULL,6},	//! -s
    {"datatype",          required_argument,NULL,7},	//! -d
    {"multiple",          required_argument,NULL,8},	//! -n
    {"input",             required_argument,NULL,9},	//! -i
    {"inputfile",         required_argument,NULL,9},	//! -i
    {"bootstrap",         required_argument,NULL,10},	//! -b
    {"ts/tv",             required_argument,NULL,11},	//! -t
    {"kappa",             required_argument,NULL,11},	//! -t
    {"nclasses",          required_argument,NULL,12},	//! -c
    {"pinv",              required_argument,NULL,13},	//! -v
    {"alpha",             required_argument,NULL,14},	//! -a
    {"inputtree",         required_argument,NULL,15},	//! -u
    {"min_diff_lk_local", required_argument,NULL,16},
    {"min_diff_lk_global",required_argument,NULL,17},
    {"steph_spr",         no_argument,NULL,18},
    {"brent_it_max",      required_argument,NULL,19},
    {"rand_start",        no_argument,NULL,20},
    {"n_rand_starts",     required_argument,NULL,21},
    {"sequential",        no_argument,NULL,22},	//! -q
    {"inside_opt",        no_argument,NULL,23},
    {"p_moves",           required_argument,NULL,24},
    {"fast_nni",          no_argument,NULL,25},
    {"g_pars",            no_argument,NULL,26},
    {"r_seed",            required_argument,NULL,27},
    {"collapse_boot",     required_argument,NULL,28},
    {"random_boot",       required_argument,NULL,29},
    {"print_trace",       no_argument,NULL,30},
    {"print_site_lnl",    no_argument,NULL,31},
    {"cov",               no_argument,NULL,32},
    {"cov_delta",         required_argument,NULL,33},
    {"cov_alpha",         required_argument,NULL,34},
    {"cov_ncats",         required_argument,NULL,35},
    {"ps",                no_argument,NULL,36},
    {"cov_free",          no_argument,NULL,37},
    {"no_gap",            no_argument,NULL,38},
    {"n_rr_branch",       required_argument,NULL,39},
    {"append",            no_argument,NULL,40},
    {"no_five_branch",    no_argument,NULL,41},
    {"pars_thresh",       required_argument,NULL,42},
    {"min_diff_lk_move",  required_argument,NULL,43},
    {"hybrid",            no_argument,NULL,44},
    {"use_median",        no_argument,NULL,45},
    {"run_id",            required_argument,NULL,46},
    {"pars",              no_argument,NULL,47},       //! -p
    {"quiet",             no_argument,NULL,48},
    {"version",           no_argument,NULL,49},
    {"calibration",       required_argument,NULL,50},
    
    {"test_init_tree",          no_argument,NULL,51}, //!<Added by Marcelo.
    {"dist_tree_model",   required_argument,NULL,52}, //!<Added by Marcelo.
    {"expm",              required_argument,NULL,53}, //!<Added by Marcelo.
    {"compare_runs",      required_argument,NULL,54}, //!<Added by Marcelo.
    {"lk_experiment",     required_argument,NULL,55}, //!<Added by Marcelo.
    {"threshold_test",    required_argument,NULL,56}, //!<Added by Marcelo.
    {"optBrent",          required_argument,NULL,57}, //!<Added by Marcelo.
    {"testCond",                no_argument,NULL,59}, //!<Added by Marcelo.
    {"NT2AA",                   no_argument,NULL,60}, //!<Added by Marcelo.
    {"pcs",               required_argument,NULL,61}, //!<Added by Marcelo.
    {"optHeuristic",      required_argument,NULL,62}, //!<Added by Marcelo.    
    {"AA2CD",             required_argument,NULL,64}, //!<Added by Marcelo.
    {"NOBL",                    no_argument,NULL,65}, //!<Added by Marcelo.
    
    {"frequencies",       required_argument,NULL,130}, //!<Added by Louis
    {"optimize",          required_argument,NULL,131}, //!<Added by Louis
    {"omega",             required_argument,NULL,132}, //!<Added by Louis
    {"genetic_code",      required_argument,NULL,133}, //!<Added by Louis
    {"help",                    no_argument,NULL,134}, //!<Added by Louis
    
    {"wvals",             required_argument,NULL,135}, //!<Added by Stefan
    {"wclasses",          required_argument,NULL,136}, //!<Added by Stefan
    {"qrates",            required_argument,NULL,137}, //!<Added by Stefan
    {"fmodel",            required_argument,NULL,138}, //!<Added by Stefan
#ifdef USEYAML
    {"yamlconfig",        required_argument,NULL,139}, //!<Added by Stefan
#endif
    
    {"oformat",           required_argument,NULL,142}, //!<Added by Stefan
    {"darwinconfig",      required_argument,NULL,143}, //!<Added by Stefan
   
    {"logtree",           required_argument,NULL,146},  //!<Added by Stefan
    {"hotness",           required_argument,NULL,147},  //!<Added by Ken
    {"root",           required_argument,NULL,148},  //!<Added by Ken
    {"motifs",           required_argument,NULL,149},  //!<Added by Ken
    {"partfile",           required_argument,NULL,150},  //!<Added by Ken
	{"ambigfile",		  required_argument,NULL,151},  //!<Added by Ken
	{"threads",		  required_argument,NULL,152},  //!<Added by Ken
	{"slowSPR",		  required_argument,NULL,153},  //!<Added by Ken
	{"stretch",		  required_argument,NULL,154},  //!<Added by Ken
    {0,0,0,0}
};


void finishOptions(option * io)
{
    checkForRandStartTree(io);
    checkModelCombinations(io);
    adaptForM4(io);
    createOutFiles(io);
    cleanupParameters(io);
    setupGeneticCode(io);
    setupFreqs(io);
    if(io->datatype == CODON) {
        setupKappa(io);
        setupInitialRateMats(io);
        setupOmegaCats(io);
    }
    setupModelIdentifier(io);
    
    Set_Model_Name(io->mod);
	
	io->fp_out_tree  = Openfile(io->out_tree_file, io->writemode);
	io->fp_out_stats = Openfile(io->out_stats_file, io->writemode);
    setupFreqHandling(io);

    //set up HLP17 model
    //Added by Ken 12/7/2016
    if(io->modeltypeOpt == HLP17){
    	io->mod->opthotness=1;
    	char *minfo1,*minfo2, *mtemp1;
    	int c,d;
    	int nh = 0;

    	//Default values
    	if(io->mod->motifstringopt==0){
    		io->mod->motifstring = "WRC_2:0,GYW_0:0";
    	}
    	if(io->mod->hotnessstringopt==0){
    	    io->mod->hotnessstring = "e";
    	}
        if((strcmp(io->mod->motifstring,"FCH")==0)){
            io->mod->motifstring = "WRC_2:0,GYW_0:1,WA_1:2,TW_0:3,SYC_2:4,GRS_0:5";
            io->mod->hotnessstring = "e,e,e,e,e,e";
        }
    	if(io->mod->rootfound==0){
    		printf("\nError: Root sequence ID must be specified using --root\n\n");
    		exit(EXIT_FAILURE);
    	}

    	//parse motif string
    	minfo1 = strdup(io->mod->motifstring);
    	minfo2 = strdup(io->mod->motifstring);

   	    	io->mod->nmotifs=0;    	//determine number of motifs
   	    	while ((mtemp1 = strsep(&minfo1, ",")) != NULL){io->mod->nmotifs++;}
  	    	io->mod->motifs = malloc(sizeof(char *) * io->mod->nmotifs);
  	    	io->mod->motif_hotness = malloc(sizeof(int ) * io->mod->nmotifs);
  	    	for(c=0;c<io->mod->nmotifs;c++){
   	    	   mtemp1 = strsep(&minfo2, ",");
   	    	   char* strtemp = strsep(&mtemp1, ":");
   	    	   io->mod->motifs[c] = malloc(sizeof(char)*10);
   	    	   strcpy(io->mod->motifs[c],strtemp);
  	    	   io->mod->motif_hotness[c] = atoi(strsep(&mtemp1, ":"));
   	    	   if(io->mod->motif_hotness[c] > nh){nh=io->mod->motif_hotness[c];}
  	    	}
  	    	nh = nh + 1;

        //parse hotness string
    	minfo1 = strdup(io->mod->hotnessstring);
    	minfo2 = strdup(io->mod->hotnessstring);
    	io->mod->nhotness=0;//determine number of h parameters
    	while ((mtemp1 = strsep(&minfo1, ",")) != NULL){io->mod->nhotness++;}
    	io->mod->hotness = (phydbl*)mCalloc(io->mod->nhotness,sizeof(phydbl));
    	io->mod->hoptindex = (int*)mCalloc(io->mod->nhotness,sizeof(int ));
    	for(c=0;c<io->mod->nhotness;c++){
    	    mtemp1 = strsep(&minfo2, ",");
    	      if (strcmp(mtemp1,"e")==0){
    	        io->mod->hoptindex[c]=1;
    	        io->mod->hotness[c]=0.0;
    	      }else{
    	        io->mod->hoptindex[c]=0;
    	        io->mod->hotness[c]=atof(mtemp1);
    	      }
    	 }

    	//catch inconsistencies between hotness and motif specification
    	if(nh != io->mod->nhotness){
    		printf("\nError in h parameter specification! Too many or too few h's.\n\n");
    		exit(EXIT_FAILURE);
    	}

    	//catch errors in hotness index specification
    	for(c=0;c<nh;c++){
    		int found = 0;
    		for(d=0;d<io->mod->nmotifs;d++){
    			if(io->mod->motif_hotness[d] == c){
    				found=1;
    			}
    		}
    		if(found==0){
    			printf("\nError in h index specification! Missing %d\n\n",c);
    			exit(EXIT_FAILURE);
    		}
    	}

    	//Read in hotspot tables
    	int mot;
    	printf(". Loading hotspot tables..\n");
    	io->mod->hotspotcmps = (phydbl**)mCalloc(io->mod->nmotifs,sizeof(phydbl *));
    	for(mot = 0; mot < io->mod->nmotifs; mot++){
    	    char *infile = mCalloc(strlen(HTABLE)+strlen(io->mod->motifs[mot])+1,sizeof(char));
    	    strcpy(infile, HTABLE);
    	    strcat(infile, io->mod->motifs[mot]);
    	    printf("%s\n",infile);
    	    FILE *file = fopen(infile, "r");
    	    if(file == NULL){
    	    	printf("\n\nCouldn't open %s\n\n",infile);
    	    	exit(EXIT_FAILURE);
    	    }
    	    phydbl *hot;
    	    int combinations = 13845841;
    	    hot = (phydbl *)mCalloc(combinations,sizeof(phydbl));//checked 15/7/2016
    	    int i=0;
    	    double num;
    	    while(i < combinations) {
    	        int fscn = fscanf(file, "%lf\n",&num);
    	        hot[i] = num;
    	        i++;
    	    }
    	    fclose(file);
    	    io->mod->hotspotcmps[mot] = hot;
    	}

    	//partition model stuff
    	if(io->mod->partfilespec==1){
	    FILE *file = fopen(io->mod->partfile, "r");
	    int npart;
	    int fscn =  fscanf(file, "%d",&io->mod->nparts);
        int nsite;
        fscn = fscanf(file, " %d\n",&nsite);

        io->mod->nomega_part=io->mod->nparts;
        io->mod->partIndex = (int *)mCalloc(nsite,sizeof(double));
        io->mod->partNames = (char**)mCalloc(io->mod->nparts,sizeof(char*));
        io->mod->omega_part = (phydbl*)mCalloc(io->mod->nomega_part,sizeof(phydbl));
        for(c=0;c<io->mod->nomega_part;c++){
        	io->mod->omega_part[c]=0.4;
        }

        int indexi;
        for(indexi=0;indexi<nsite;indexi++){
        	io->mod->partIndex[indexi]=-1;
        }

        ssize_t read;
        size_t len=0;
        char *linepart = NULL;

    	int parti;
    	for(parti=0;parti<io->mod->nparts;parti++){
    		read = getline(&linepart,&len,file);
    		char* l1 = strsep(&linepart,":");
    		printf("%s\t",l1);
    		io->mod->partNames[parti] = (char*)mCalloc(T_MAX_OPTION,sizeof(char));
    		strcpy(io->mod->partNames[parti],l1);
    		//io->mod->partNames[parti] = l1;
    		char* l2 = strsep(&linepart,"\n");
    		printf("%s\n",l2);
    		char *ltemp1;
    		while ((ltemp1 = strsep(&l2, ",")) != NULL){
				 int start = atoi((strsep(&ltemp1,".")));
				 (strsep(&ltemp1,"."));
				 int end = atoi((strsep(&ltemp1,"\n")));
				 int tempcount;
				 for(tempcount=start;tempcount<=end;tempcount++){
					 if(io->mod->partIndex[tempcount]!=-1){
						 printf("\nPosition %d specified more than once in partition file!\n",tempcount);
						 exit(EXIT_FAILURE);
					 }
					 io->mod->partIndex[tempcount]=parti;
				 }
    		}
    	}

    	for(indexi=0;indexi<nsite;indexi++){
    	   	if(io->mod->partIndex[indexi] == -1){
    	   		printf("\nPosition %d not specified in partition file!\n",indexi);
    	   		exit(EXIT_FAILURE);
    	   	}
    		printf("%d ",io->mod->partIndex[indexi]);
    	}
    	printf("\n");

     }else{
    	 io->mod->nparts=1;
    	 io->mod->nomega_part=io->mod->nparts;
    	 io->mod->omega_part[0]=0.4;
     }//partfilespec==1
    }
}

void checkForRandStartTree(option * io)
{
    if((io->mod->s_opt->n_rand_starts)           && 
       (io->mod->s_opt->topo_search == NNI_MOVE) && 
       (io->mod->s_opt->random_input_tree)) {
        Warn_And_Exit("\n. The random starting tree option is only compatible with SPR based search options.\n"); 
    }
}

void checkModelCombinations(option * io)
{
    if((io->datatype == NT) && (io->modeltypeOpt > 10)) {
        Warn_And_Exit("\n. Err: model incompatible with the data type. Please use JC69, K80, F81, HKY, F84, TN93 or GTR\n");
    } else if((io->datatype == AA) && ((io->modeltypeOpt < 11) || io->modeltypeOpt < 0)) {
        Warn_And_Exit("\n. Err: model incompatible with the data type. Please use LG, Dayhoff, JTT, MtREV, WAG, DCMut, RtREV, CpREV, VT, Blosum62, MtMam, MtArt, HIVw or HIVb.\n");
    } else if((io->datatype == CODON) && (io->modeltypeOpt > 0)) { //!< Added by Marcelo.
        Warn_And_Exit("\n. Err: model incompatible with the data type. Please use GY.\n");
    }
}

void adaptForM4(option * io)
{
    if(io->m4_model == YES) {
#ifdef M4
        io->mod->ns *= io->mod->m4mod->n_h;
        io->mod->use_m4mod = 1;
        M4_Make_Complete(io->mod->m4mod->n_h,
                         io->mod->m4mod->n_o,
                         io->mod->m4mod);
#endif
    } else {
        io->mod->s_opt->opt_cov_delta      = 0;
        io->mod->s_opt->opt_cov_alpha      = 0;
        io->mod->s_opt->opt_cov_free_rates = 0;
    }
    
    if((io->mod->s_opt->opt_cov_free_rates) && (io->mod->s_opt->opt_cov_alpha))
    {
        io->mod->s_opt->opt_cov_free_rates = 0;
        io->mod->m4mod->use_cov_alpha      = 0;
        io->mod->m4mod->use_cov_free       = 1;
    }
}

void createOutFiles(option * io)
{
    if(io->print_site_lnl) {
        char *ext;
#ifdef USEYAML
        if(io->out_stats_format == OUTYAML) {
            ext = ".yml";
        } else
#endif
        if(io->out_stats_format == OUTDARWIN) {
            ext = ".drw";
        } else {
            ext = ".txt";
        }
        io->fp_out_lk = openOutputFile(io->out_lk_file, "_igphyml_lk", ext, io);
    }
    if(io->print_trace) {
        io->fp_out_tree_trace = openOutputFile(io->out_trace_tree_file, "_igphyml_tree_trace", ".txt", io);
        io->fp_out_stats_trace = openOutputFile(io->out_trace_stats_file, "_igphyml_stats_trace", ".txt", io);
    }
    if(io->mod->s_opt->random_input_tree) {
        io->fp_out_trees = openOutputFile(io->out_trees_file, "_igphyml_rand_trees", ".txt", io);
    }
    if((io->print_boot_trees) && (io->mod->bootstrap > 0)) {
        io->fp_out_boot_tree = openOutputFile(io->out_boot_tree_file, "_igphyml_boot_trees", ".txt", io);
        io->fp_out_boot_stats = openOutputFile(io->out_boot_stats_file, "_igphyml_boot_stats", ".txt", io);
    }
    
    if(io->append_run_ID) {
        strcat(io->out_tree_file, "_");
        strcat(io->out_stats_file, "_");
        strcat(io->out_tree_file, io->run_id_string);
        strcat(io->out_stats_file, io->run_id_string);
    }
    
#ifdef USEYAML
    if(io->out_stats_format == OUTYAML) {
        strcat(io->out_stats_file, ".yml");
    } else
#endif
    if(io->out_stats_format == OUTDARWIN) {
        strcat(io->out_stats_file, ".drw");
    } else {
        strcat(io->out_stats_file, "");
    }
    strcat(io->out_tree_file, "");
}

void cleanupParameters(option * io)
{
    if(io->datatype == NT) {
        io->mod->ns         = 4;
        io->mod->state_len  = 1;
    } else if(io->datatype == AA) {
        io->mod->state_len        = 1;
        io->mod->s_opt->opt_kappa = 0;
        io->mod->ns               = 20;
        
        if(io->modeltypeOpt != K80 && 
           io->modeltypeOpt != HKY85 && 
           io->modeltypeOpt != F84 &&
           io->modeltypeOpt != TN93) {
            io->mod->s_opt->opt_kappa = 0;
            io->mod->s_opt->opt_omega = 0;
        }
    } else if(io->datatype == CODON) {
        io->mod->state_len        = 1;
    }
    
    if(io->mod->n_catg == 1) {
        io->mod->s_opt->opt_alphaCD = io->mod->s_opt->opt_alpha = 0;
    }
    
    if(!io->mod->s_opt->opt_subst_param) {
        io->mod->s_opt->opt_alpha  = 0;
        io->mod->s_opt->opt_kappa  = 0;
        io->mod->s_opt->opt_lambda = 0;
        io->mod->s_opt->opt_pinvar = 0;
        io->mod->s_opt->opt_rr     = 0;	
        io->mod->s_opt->opt_omega  = 0; //!< Added by Marcelo.
        io->mod->s_opt->opt_beta   = 0; //!< Added by Marcelo.
        io->mod->s_opt->opt_alphaCD= 0; //!< Added by Marcelo.
    }
}

// we are already using codon models at this point
void setupKappa(option *io)
{
    if(io->userWantsKappa == NO) {
        io->mod->s_opt->opt_kappa      = 0;
    }
}

void setupGeneticCode(option * io)
{
    if(io->datatype == CODON) {
        if((io->mod->initqrates == KOSI07)   ||
           (io->mod->initqrates == SCHN05)  ||
           (io->modeltypeOpt == PCM)) { //!< Those models assume the standard genetic code.
            io->mod->genetic_code = STANDARD;
	    io->genCode = STANDARD;
	    Set_Genetic_Code(io->mod->genetic_code);
        } else {
            if(io->genCode == NOGENCODE) { io->genCode = STANDARD; }
            io->mod->genetic_code = io->genCode;
	    Set_Genetic_Code(io->genCode);
	    if((io->genCode!=STANDARD) &&  
	       (io->init_DistanceTreeCD == KOSI07 || 
	        io->init_DistanceTreeCD == SCHN05)) { 
		io->init_DistanceTreeCD = NUCLEO; 
		if(!io->quiet) PhyML_Printf("\n. [Warning] Use --dist_tree_model ECMUSR if you need to use your own model to estimate initial pairwise distances under a non-standard genetic code.\n");  
		}
        } 
        
        Genetic_code_index_stopCodons(io->mod->genetic_code); 
        io->mod->ns = Genetic_code_ns();
    } else if(io->datatype == NT) {
        io->mod->ns = 4;
        io->mod->genetic_code = STANDARD;
        io->genCode = STANDARD;
    } else if(io->datatype == AA) {
        io->mod->ns = 20;
        io->mod->genetic_code = STANDARD;
        io->genCode = STANDARD;
    }
}

void setupFreqs(option * io)
{
    if(io->freqmodelOpt == FUNDEFINED) { 
        io->freqmodelOpt = FMODEL;
    }
    io->mod->freq_model = io->freqmodelOpt;
    switch(io->mod->freq_model) {
        case F1XSENSECODONS: {
            io->mod->ns            = Genetic_code_ns();
            io->mod->num_base_freq = io->mod->ns;
            break;
        }
        case F1X4: {
            io->mod->num_base_freq = 4;
            break;
        }
        case FMODEL: {
            io->mod->num_base_freq = 0;
            break;
        }
        default: { // F3X4, CF3X4
            io->mod->num_base_freq = 12;
            break;
        }
    }
    
    

    
    // check for user defined frequencies
    if(io->userFreqs[0] != 0) {
        phydbl sum;
        double val1, val2, val3, val4;
        double val11, val12, val13, val14, val21, val22, val23, val24, val31, val32, val33, val34;
        // double val41, val42, val43, val44;
        phydbl freqs[64] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; 
        int i, j, k;
        switch(io->mod->freq_model) {
            case F1XSENSECODONS: {
                k = sscanf(io->userFreqs, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                           &freqs[0], &freqs[1], &freqs[2], &freqs[3], &freqs[4], &freqs[5], &freqs[6], &freqs[7], &freqs[8], &freqs[9], &freqs[10],
                           &freqs[11], &freqs[12], &freqs[13], &freqs[14], &freqs[15], &freqs[16], &freqs[17], &freqs[18], &freqs[19], &freqs[20],
                           &freqs[21], &freqs[22], &freqs[23], &freqs[24], &freqs[25], &freqs[26], &freqs[27], &freqs[28], &freqs[29], &freqs[30],
                           &freqs[31], &freqs[32], &freqs[33], &freqs[34], &freqs[35], &freqs[36], &freqs[37], &freqs[38], &freqs[39], &freqs[40],
                           &freqs[41], &freqs[42], &freqs[43], &freqs[44], &freqs[45], &freqs[46], &freqs[47], &freqs[48], &freqs[49], &freqs[50],
                           &freqs[51], &freqs[52], &freqs[53], &freqs[54], &freqs[55], &freqs[56], &freqs[57], &freqs[58], &freqs[59], &freqs[60],
                           &freqs[61], &freqs[62], &freqs[63]
                           );
                if(k != 64) {
                    PhyML_Printf("Too few codon frequencies provided. Expected 64 given %d", k);
                    Warn_And_Exit("");
                }
                sum = 0.0;
                For(i, 64) sum += freqs[i];
                j = -1;
                For(i, 64) { 
                	PhyML_Printf("Read in %i: %f",i,freqs[i]); //added by Ken
                    freqs[i] /= sum;
                    if(freqs[i] < .0 || freqs[i] > 1.) {
                        Warn_And_Exit("\n. Invalid base frequencies.\n");
                    }
                    if(!stopCodons[i]) io->mod->user_b_freq[++j] = freqs[i];
                }
                


                break;
            }
                
            case F1X4: {
                k = sscanf(io->userFreqs,"%lf,%lf,%lf,%lf",&val1,&val2,&val3,&val4); //!< FT FC FA FG.
                if(k!=4) {
                    PhyML_Printf("Too few codon frequencies provided. Expected 4 given %d", k);
                    Warn_And_Exit("");
                }
                
                io->mod->user_b_freq[0] = (phydbl)val1; 
                io->mod->user_b_freq[1] = (phydbl)val2; 
                io->mod->user_b_freq[2] = (phydbl)val3; 
                io->mod->user_b_freq[3] = (phydbl)val4; 
                
                sum = (io->mod->user_b_freq[0] +
                       io->mod->user_b_freq[1] +
                       io->mod->user_b_freq[2] +
                       io->mod->user_b_freq[3]);
                
                io->mod->user_b_freq[0] /= sum;
                io->mod->user_b_freq[1] /= sum;
                io->mod->user_b_freq[2] /= sum;
                io->mod->user_b_freq[3] /= sum;
                
                
                if(io->mod->user_b_freq[0] < .0 ||
                   io->mod->user_b_freq[1] < .0 ||
                   io->mod->user_b_freq[2] < .0 ||
                   io->mod->user_b_freq[3] < .0 ||
                   io->mod->user_b_freq[0] > 1. ||
                   io->mod->user_b_freq[1] > 1. ||
                   io->mod->user_b_freq[2] > 1. ||
                   io->mod->user_b_freq[3] > 1.) {
                    Warn_And_Exit("\n. Invalid base frequencies.\n");
                }
                break;
            }
                
            case F3X4: case CF3X4: {
                k=sscanf(io->userFreqs, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                         &val11, &val12, &val13, &val14, 
                         &val21, &val22, &val23, &val24,
                         &val31, &val32, &val33, &val34); //!< FT1 FC1 FA1 FG1 FT2 FC2 FA2 FG2 FT3 FC3 FA3 FG3.
                if(k!=12) {
                    PhyML_Printf("Too few codon frequencies provided. Expected 12 given %d", k);
                    Warn_And_Exit("");
                }
                io->mod->user_b_freq[0] = (phydbl)val11; 
                io->mod->user_b_freq[1] = (phydbl)val12; 
                io->mod->user_b_freq[2] = (phydbl)val13; 
                io->mod->user_b_freq[3] = (phydbl)val14; 
                io->mod->user_b_freq[4] = (phydbl)val21; 
                io->mod->user_b_freq[5] = (phydbl)val22; 
                io->mod->user_b_freq[6] = (phydbl)val23; 
                io->mod->user_b_freq[7] = (phydbl)val24;
                io->mod->user_b_freq[8] = (phydbl)val31; 
                io->mod->user_b_freq[9] = (phydbl)val32; 
                io->mod->user_b_freq[10] = (phydbl)val33; 
                io->mod->user_b_freq[11] = (phydbl)val34; 		  
                
                sum = 0.0;
                For(i, 4) sum += io->mod->user_b_freq[i];
                For(i, 4) io->mod->user_b_freq[i] /= sum;
                
                sum = 0.0;
                for(i=4; i<8; i++) sum += io->mod->user_b_freq[i];
                for(i=4; i<8; i++) io->mod->user_b_freq[i] /= sum;		  
                
                sum = 0.0;
                for(i=8; i<12; i++) sum += io->mod->user_b_freq[i];
                for(i=8; i<12; i++) io->mod->user_b_freq[i] /= sum;		  
                
                For(i, 12) {
                    if(io->mod->user_b_freq[i] < .0 || io->mod->user_b_freq[i] > 1.) {
                        Warn_And_Exit("\n. Invalid base frequencies.\n");
                    }
                }
                
                break;
            }		
            default:
                break;
        }
    }
    
    switch(io->eq_freq_handling) {
        case MODEL: {
            io->mod->freq_model = FMODEL;
            io->mod->s_opt->opt_state_freq = NO; 
            io->mod->s_opt->user_state_freq = NO;
            break;
        }
            
        case OPTIMIZE: {
            io->mod->s_opt->opt_state_freq = YES;
            io->mod->s_opt->user_state_freq = NO;
            break;
        }
            
        case USER: {
            io->mod->s_opt->opt_state_freq  = NO;
            io->mod->s_opt->user_state_freq = YES;
            break;
        }
            
        case NOFEQ:
        case EMPIRICAL:
        default: {
            io->mod->s_opt->opt_state_freq = NO;
            io->mod->s_opt->user_state_freq = NO;
            break;
        }
    }
    
    free(io->userFreqs);  
}

void setupInitialRateMats(option * io)
{ 
    if((io->mod->initqrates == KOSI07 || io->mod->initqrates == SCHN05) &&
       (io->mod->genetic_code != STANDARD)) {
        PhyML_Printf("This init rate matrix demands for STANDARD genetic code (changed).\n");
        io->mod->genetic_code = STANDARD;
    }
    
    if(io->mod->initqrates == ECMUSR) {
        // read user defined matrix
        int i,j;
        if (!Filexists ("usermatrix.ecm")) {
            char tmp[256];
            char choix;
            strcpy (tmp, "\n. The file '");
            strcat (tmp, "usermatrix.ecm");
            strcat (tmp, "' does not exist.\n");
            PhyML_Printf("%s",tmp);
            PhyML_Printf("\n. Type any key to exit.\n");
            if(!scanf("%c",&choix)) Exit("\n");
            Exit("\n");
        } else {
            io->fp_in_usrECM = Openfile("usermatrix.ecm", 0);
        }
        
        For(i, 64) {
            For(j,64) io->mod->userRates[i][j] = 0.0;
            io->mod->userfreq[i] = 0.0;
            io->mod->userbfreq[i] = 0.0;
        }
        
        if(io->modeltypeOpt == MG && io->mod->freq_model == FMODEL) { //is this correct? Ken 4/1/2017
            Read_userRatesAndFreqs((phydbl*)io->mod->userRates, (phydbl*)io->mod->userfreq, 64, io->fp_in_usrECM);
        } else {
            Read_userRatesAndFreqsMG((phydbl*)io->mod->userRates, (phydbl*)io->mod->userfreq, (phydbl*)io->mod->userbfreq, 12, 64, io->fp_in_usrECM);
        }
        
        For(i,64) {
            For(j,64) io->mod->userRatesT[i][j] = io->mod->userRates[i][j];
            io->mod->userfreqT[i] = io->mod->userfreq[i];
        }
    }
}

void setupOmegaCats(option * io)
{
    io->mod->s_opt->opt_omega = 0;
    
    if (io->omegaOpt == DM0 && io->wclasses > 1) {
        io->wclasses = 1;
    }
    
    io->mod->n_w_catg = io->wclasses;
    if(io->mod->n_w_catg > 1) {
        io->mod->n_catg = io->mod->n_w_catg;
    } else if(io->mod->n_catg > 1) {
        io->mod->n_w_catg = 1;
    }
    
    io->mod->omegaSiteVar = io->omegaOpt;
    
    io->mod->prob_omegas           = (phydbl *)mCalloc(io->mod->n_w_catg,sizeof(phydbl));
    io->mod->omegas                = (phydbl *)mCalloc(io->mod->n_w_catg,sizeof(phydbl));
    
    switch(io->mod->omegaSiteVar) {
        case DMODELK: {
            int i;
            For(i,io->mod->n_w_catg) {
                io->mod->omegas[i]           = io->mod->n_w_catg * (.5+i*2./io->mod->n_w_catg*(0.8+0.4*rndu()));
                io->mod->prob_omegas[i]      = FABS(rndu()); //(phydbl)i+0.1;
	    }
            io->mod->s_opt->opt_omega      = 1;
            io->mod->s_opt->opt_prob_omega = 1;
            io->mod->s_opt->opt_alpha      = 0;
            io->mod->s_opt->opt_beta       = 0;
            break;
        }
        case DGAMMAK: {
            int i;
            For(i,io->mod->n_w_catg) {
                io->mod->omegas[i]           = 1.0;
                io->mod->prob_omegas[i]      = 1.0/(phydbl)io->mod->n_w_catg;
            }
            io->mod->alpha                 = 1.1;
            io->mod->beta                  = 1.1;
            io->mod->s_opt->opt_omega      = 1;
            io->mod->s_opt->opt_prob_omega = 0;
            io->mod->s_opt->opt_alpha      = 1;
            io->mod->s_opt->opt_beta       = 1;
            break;
        }
        case DM0: {
        	io->mod->omegas[0]               = 1.0;
            io->mod->prob_omegas[0]          = 1.0;
            io->mod->omega_part = (phydbl*)mCalloc(1,sizeof(phydbl));
           	io->mod->omega_part[0]=0.4;
            io->mod->s_opt->opt_omega        = 1;
            io->mod->s_opt->opt_prob_omega   = 0;
            io->mod->s_opt->opt_alpha        = 0;
            io->mod->s_opt->opt_beta         = 0;

            break;
        }
        case NOOMEGA: default: {
            io->mod->omegas[0]               = 1.0;
            io->mod->prob_omegas[0]          = 1.0;
           io->mod->omega_part[0]                   = 1.0; //Ken 18/8
            io->mod->s_opt->opt_omega        = 0;
            io->mod->s_opt->opt_prob_omega   = 0;
            io->mod->s_opt->opt_alpha        = 0;
            io->mod->s_opt->opt_beta         = 0;
        }
    }
    
    // check for user defined values for omega
    if(io->userOmega[0] != 0) {
        if(io->mod->omegaSiteVar == DM0) {
            io->mod->omega_part[0] = atof(io->userOmega);
            io->mod->s_opt->opt_omega = 0;
        } else if(io->mod->omegaSiteVar == DMODELK) {
            char *tempStr;
            tempStr=(char *)mCalloc(T_MAX_NAME,sizeof(char));
            int j,k,m,i = 0;
            if((io->userOmega[i]!='w')||(io->userOmega[i]!='p')||(io->userOmega[i]!='b')) {
                phydbl *values;
                int n;
                //if((io->userOmega[i]!='w')||(io->userOmega[i]!='p')) n=1;else n=2;
                n = 2;
                m=i;
                values=(phydbl *)mCalloc(n*io->mod->n_w_catg,sizeof(phydbl));
                
                For(k,n*io->mod->n_w_catg) {
                    j=-1;
                    while((io->userOmega[i]!=',')&&(io->userOmega[i]!=0) ) {
                        tempStr[++j]=io->userOmega[i++];
                    }
                    i++;
                    tempStr[++j]=0;
                    values[k]=(phydbl)atof(tempStr);
                }
                
                phydbl sum=0.0;
                
                if(io->userOmega[m]=='p') {
                    for(i = 0; i < io->mod->n_w_catg; i++) { sum                    += FABS(values[i]); }
                    for(i = 0; i < io->mod->n_w_catg; i++) { io->mod->prob_omegas[i] = FABS(values[i])/sum; }
                    io->mod->s_opt->opt_omega                                        = 1;
                    io->mod->s_opt->opt_prob_omega                                   = 0;
                    for(i=0; i < io->mod->n_w_catg; i++) { io->mod->omegas[i]        = 1.0; }
                } else if(io->userOmega[m]=='w') {
                    for(i=0; i < io->mod->n_w_catg; i++) { io->mod->prob_omegas[i]   = 1.0/(phydbl)io->mod->n_w_catg; }
                    for(i=0; i < io->mod->n_w_catg; i++) { io->mod->omegas[i]        = values[i]; }
                    io->mod->s_opt->opt_omega                                        = 0;
                    io->mod->s_opt->opt_prob_omega                                   = 1;
                } else {
                    for(i=0; i < io->mod->n_w_catg; i++) { io->mod->omegas[i]      = values[i]; }
                    for(j=i; j<2*i; j++) { sum                                     += FABS(values[j]); }
                    for(j=i; j<2*i; j++) { io->mod->prob_omegas[j-i]                = FABS(values[j])/sum; }
                    io->mod->s_opt->opt_omega                                        = 0;
                    io->mod->s_opt->opt_prob_omega                                   = 0;
                }
                
                io->mod->s_opt->opt_alpha                                      = 0;
                io->mod->s_opt->opt_beta                                       = 0;
                io->mod->omegaSiteVar                                          = DMODELK;
                
                free(values);
            }
            free(tempStr);
        } else if(io->mod->omegaSiteVar == DGAMMAK) {
            char *tempStr;
            tempStr=(char *)mCalloc(T_MAX_NAME,sizeof(char));
            int j,k,i=0;
            phydbl * values;
            values=(phydbl *)mCalloc(2,sizeof(phydbl));
            
            For(k,2) {
                j=-1;
                while((io->userOmega[++i]!=',')&&(io->userOmega[i]!=0) ) {
                    tempStr[++j]                                  = io->userOmega[i];
                }
                tempStr[++j]                                    = 0;
                values[k]=(phydbl)atof(tempStr);
            }	
            io->mod->alpha                                    = values[0];
            io->mod->beta                                     = values[1];
            
            io->mod->s_opt->opt_omega                         = 0;
            io->mod->s_opt->opt_prob_omega                    = 0;
            io->mod->s_opt->opt_alpha                         = 0;
            io->mod->s_opt->opt_beta                          = 0;
            io->mod->omegaSiteVar                             = DGAMMAK;
            free(values);
            free(tempStr);
        }
    }
    free(io->userOmega);
}
void setupModelIdentifier(option * io)
{
    // THIS SHOULD NOT BE HERE. THIS SHOULD NOT BE NEEDED - BETTER TO GET RID OF THIS
    // SOONER THAN LATER. BUT IT HAS CONSEQUENCES IN LOTS OF OTHER PLACES, SO...
    //
    if(io->datatype == CODON) {
        io->mod->whichrealmodel = io->modeltypeOpt;
        if(io->modeltypeOpt == GY) {
            io->mod->whichmodel = GYSIMP;
            if(io->mod->initqrates == KOSI07) {
                io->mod->whichmodel -= 3;
            } else if(io->mod->initqrates == SCHN05) {
                io->mod->whichmodel -= 7;
            } else if(io->mod->initqrates == ECMUSR) {
                io->mod->whichmodel -= 31;
            }
        } else if(io->modeltypeOpt == MG) {
            io->mod->whichmodel = MGSIMP;
            if(io->mod->initqrates == KOSI07) {
                io->mod->whichmodel -= 10;
            } else if(io->mod->initqrates == SCHN05) {
                io->mod->whichmodel -= 14;
            } else if(io->mod->initqrates == ECMUSR) {
                io->mod->whichmodel -= 34;
            }
        } else if(io->modeltypeOpt == YAP) {
            io->mod->whichmodel = YAPSIMP;
            if(io->mod->initqrates == KOSI07) {
                io->mod->whichmodel -= 17;
            } else if(io->mod->initqrates == SCHN05) {
                io->mod->whichmodel -= 21;
            } else if(io->mod->initqrates == ECMUSR) {
                io->mod->whichmodel -= 25;
            }
        }
        if(io->mod->initqrates != NOINITMAT) {
            if(io->mod->freq_model != FMODEL) {
                io->mod->whichmodel -= 2;
            }
            if(io->mod->s_opt->opt_omega == 1) {
                io->mod->whichmodel -= 1;
            }
        }
        if(io->modeltypeOpt == PCM) {
            io->mod->whichmodel = GYECMUSRF;
            io->mod->pcaModel = YES;
        } else {
            io->mod->pcaModel = NO;
        }
    } else {
        io->mod->whichmodel = io->modeltypeOpt;
        io->mod->genetic_code = io->genCode;
    }
}

void setupFreqHandling(option * io)
{
    if(io->datatype==NT) { 
        switch(io->eq_freq_handling) {
            case EMPIRICAL: 
                io->mod->s_opt->opt_state_freq = NO; 
                break;
                
            case OPTIMIZE:  
                io->mod->s_opt->opt_state_freq = YES; 
                break;
                
            case MODEL:   
                io->mod->freq_model = FMODEL;
                io->mod->s_opt->opt_state_freq = NO; 
                io->mod->s_opt->user_state_freq = NO;
                break;
                
            case USER: 
                break;
                
            default:
                PhyML_Printf("\nOption for frequency (-f option) not set; assume empirical frequencies.\n");
                io->eq_freq_handling = EMPIRICAL;
                io->mod->s_opt->opt_state_freq = NO;
                break;
        }	  
	} else if(io->datatype==AA) {
        switch(io->eq_freq_handling) {
            case EMPIRICAL: 
                io->mod->s_opt->opt_state_freq = YES; 
                io->mod->s_opt->opt_state_freq_AAML = NO;
                break;
                
            case OPTIMIZE:  
                io->mod->s_opt->opt_state_freq = YES; 
                io->mod->s_opt->opt_state_freq_AAML = YES;
                break;
                
            case MODEL:   
                io->mod->s_opt->opt_state_freq = NO; 
                io->mod->s_opt->user_state_freq = NO;
                break;
                
            case USER: 
                io->mod->s_opt->opt_state_freq = NO; 
                io->mod->s_opt->user_state_freq = YES;
                break;
                
            default:
                PhyML_Printf("\nOption for frequency (-f option) not set; assume empirical frequencies.\n");
                io->eq_freq_handling = EMPIRICAL;
                io->mod->s_opt->opt_state_freq = YES; 
                io->mod->s_opt->opt_state_freq_AAML = NO;
                break;
        }	
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// MAIN OPTION SWITCH
//

int mainOptionSwitch(int opt, char * optarg, option * io)
{
    int result = 0;
    int i = 0;
    char *errorMsg = NULL;
    char *tmp = NULL;
    

    switch(opt) {
            //////////////////////////////////////////////////////////////////////////////////////
            // --optHeuristic
            //
        case 62: {
            int k;
            if((isalpha(optarg[0])) || (atoi(optarg) < 0) || (atoi(optarg) > 10)) {
                errorMsg = "\n. Unknown argument to --optHeuristic option.\n";
            } else { 						 
                io->opt_heuristic_manuel = YES;
                k = sscanf(optarg,"%d,%d",&io->roundMax_start, &io->roundMax_end);
                if(k != 2) errorMsg = "Optimization heuristic requires paremeters n,m.";
            }
            break;
        }
           
                //////////////////////////////////////////////////////////////////////////////////////
            // --NOBL ignore given branch lengths
            //!< Added by Marcelo.
            //
        case 65: {
            io->nobl=1;
            break;
        }
           
            //////////////////////////////////////////////////////////////////////////////////////
            // --pcs
            //!< Added by Marcelo.
            //
        case 61: {
            if((!atoi(optarg)) || (atoi(optarg) < 1) || (atoi(optarg) > 100)) {
                errorMsg = "\n. Unknown argument to --pcs option: the number of principal components must be a positive integer between 1 and 100\n";
            } else { 						 
                io->npcs = atoi(optarg);
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --NT2AA
            //!< Added by Marcelo.
            //
        case 60: {
            if(optarg == NULL || !strcmp(optarg, "true")) {
                io->convert_NT_to_AA = YES;
            }
            break;	  
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --testCond
            //!< Added by Marcelo.
            //
        case 59: {
            if(optarg == NULL || !strcmp(optarg, "true")) {
                io->testcondition = YES;
            }
            break;	  
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --optPHYML
            //!< Added by Marcelo.
            //
        case 57: {
            io->mod->s_opt->opt_method = optPHYML;
            
            if ((!atoi(optarg)) || (atoi(optarg) < 1)) {
                errorMsg = "\n. Unknown argument to --optBrent option: the number of cycles with Brent must be a positive integer between 1 and 100\n";
            } else  { 						 
                io->mod->s_opt->nBrentCycles = atoi(optarg);
            }
            break;	  
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --test_init_tree
            //!< Added by Marcelo.
            //
        case 51: {
            if(optarg == NULL || !strcmp(optarg, "true")) {
                io->testInitTree = YES;
            }
            break;	  
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --threshold_test
            //!< Added by Marcelo.
            //
        case 56: {
            io->threshold_exp = YES;
            io->dataset       = atoi (optarg);
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --lk_experiment
            //!< Added by Marcelo.
            //
        case 55: {
            io->lkExperiment = YES;
            sscanf(optarg,"%lf,%lf,%lf", &io->lkExpStepSize, &io->minParam, &io->maxParam);
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --expm
            //!< Added by Marcelo.
            //
        case 53: {
            if(strcmp(optarg, "EIGEN\0") == 0) {
                io->expm = EIGEN;
            } else if(strcmp(optarg, "SSPADE\0") == 0) {
                io->expm = SSPADE; //!< See MATLAB EXPM. 
            } else if(strcmp(optarg, "TAYLOR\0") == 0) {
                //io->expm = TAYLOR; //!< The simplest way.
                io->heuristicExpm = YES; //!< The simplest way.
            } else errorMsg = "Exponential matrix calculation method not implemented. Impossible to Continue.";
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --compare_runs
            //!< Added by Marcelo.
            //
        case 54: { 
            io->fp_out_compare=fopen(optarg,"a+"); 
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --dist_tree_model
            //!< Added by Marcelo.
            //
        case 52: {
            if(strcmp(optarg, "KOSI07") == 0) {
                io->init_DistanceTreeCD = KOSI07; //!< See Kosiol 2007.
                io->mod->genetic_code = STANDARD;
                Genetic_code_index_stopCodons(io->mod->genetic_code);
                io->mod->ns=Genetic_code_ns();
            } else if(strcmp(optarg, "SCHN05") == 0) {
                io->init_DistanceTreeCD= SCHN05; //!< See Schneider 2005.
                io->mod->genetic_code = STANDARD;
                Genetic_code_index_stopCodons(io->mod->genetic_code);
                io->mod->ns=Genetic_code_ns();
            } else if(strcmp(optarg, "JC69") == 0) {
                io->init_DistanceTreeCD = NUCLEO; //!< JC69 codon adapted
            } else if(strcmp(optarg, "ECMUSR") == 0) {
                io->init_DistanceTreeCD = ECMUSR; //!< USR matrix and frequencies	  
            } else {
                errorMsg = "Pairwise comparison method not implemented. Impossible to Continue.";
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --genetic_code
            //!< Added by Marcelo.
            //
        case 'g': case 133: {
            if(strcmp(optarg, "STANDARD") == 0) {
                io->genCode = STANDARD;
            } else if(strcmp(optarg, "TVMC") == 0) {
                io->genCode = TVMC;
            } else if(strcmp(optarg, "TYMC") == 0) {
                io->genCode = TYMC;
            } else if(strcmp(optarg, "THMPCMCMSC") == 0) {
                io->genCode = THMPCMCMSC;
            } else if(strcmp(optarg, "THIMC") == 0) {
                io->genCode = THIMC;
            } else if(strcmp(optarg, "THCDHNC") == 0) {
                io->genCode = THCDHNC;
            } else if(strcmp(optarg, "THEFMC") == 0) {
                io->genCode = THEFMC;
            } else if(strcmp(optarg, "THENC") == 0) {
                io->genCode = THENC;
            } else if(strcmp(optarg, "THBAPPC") == 0) {
                io->genCode = THBAPPC;
            } else if(strcmp(optarg, "THAYNC") == 0) {
                io->genCode = THAYNC;
            } else if(strcmp(optarg, "THAMC") == 0) {
                io->genCode = THAMC;
            } else if(strcmp(optarg, "THAFMC") == 0) {
                io->genCode = THAFMC;
            } else if(strcmp(optarg, "BLNC") == 0) {
                io->genCode = BLNC;
            } else if(strcmp(optarg, "CHMC") == 0) {
                io->genCode = CHMC;
            } else if(strcmp(optarg, "TRMC") == 0) {
                io->genCode = TRMC;
            } else if(strcmp(optarg, "SCOMC") == 0) {
                io->genCode = SCOMC;
            } else if(strcmp(optarg, "THMC") == 0) {
                io->genCode = THMC;
            } else {
                errorMsg = "Genetic code not implemented. Impossible to Continue.";
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --wclasses
            // number of classes for DGAMMA or DMODEL
            //!< Added by Stefan.
            //
        case 136: {
            io->wclasses = atoi(optarg);
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --omega
            //!< Added by Marcelo.
            //
        case 'w': case 132: {
            if(strcmp(optarg, "d") == 0 || strcmp(optarg, "DMODEL") == 0) { // Discrete unconstrained distribution model
                io->omegaOpt = DMODELK;
                io->wclasses = 4;
            } else if(strcmp(optarg, "g") == 0 || strcmp(optarg, "DGAMMA") == 0) { // Discrete gamma model
                io->omegaOpt = DGAMMAK;
                io->wclasses = 4;
            } else {  // single omega value, default case
            	io->omegaOpt = DM0;
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --wvals
            // user supplied values for omega(s)
            //!< Added by Stefan.
            //
        case 135: {
            for(i=0; optarg[i] != 0; i++) {
            	printf("omegas: %d\n",optarg[i]);
                io->userOmega[i] = optarg[i];
            }
        	printf("omegas s: %s\n",io->userOmega);

            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --calibration
            //
        case 50: {
            strcpy(io->clade_list_file, optarg);
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --version
            //
        case 49: {
            PhyML_Printf("IgPhyML %s\n", VERSION);
            Exit("");
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --quiet
            //
        case 48: {
            if(optarg == NULL || !strcmp(optarg, "true")) {
                io->quiet = 1;
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --pars
            //
        case 'p': case 47: {
            if(optarg == NULL || !strcmp(optarg, "true")) {
                io->in_tree = 1;
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --runid
            //
        case 46: {
            io->append_run_ID = 1;
            strcpy(io->run_id_string, optarg);
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --use_median
            //
        case 45: {
            if(optarg == NULL || !strcmp(optarg, "true")) {
                io->mod->gamma_median = 1;
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --hybrid
            //
        case 44: {
            if(optarg == NULL || !strcmp(optarg, "true")) {
                io->mod->s_opt->hybrid_thresh = 0;
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --min_diff_lk_move
            //
        case 43: {
            io->mod->s_opt->min_diff_lk_move = atof(optarg);
            if(io->mod->s_opt->min_diff_lk_move < 0) {
                errorMsg = "\n. Min_diff_lk_move must be a double greater than 0.\n";
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --pars_thresh
            //
        case 42: {
            io->mod->s_opt->pars_thresh = (int)atoi(optarg);
            /* 	      if(io->mod->s_opt->pars_thresh < 0) */
            /* 		{ */
            /* 		  PhyML_Printf("\n. The parsimony threshold must be an integer greater than 0.\n"); */
            /* 		  PhyML_Printf("\n. Type any key to exit.\n"); */
            /* 		  Exit("\n"); */
            /* 		} */
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --no_five_branch
            //
        case 41: {
            if(optarg == NULL || !strcmp(optarg, "true")) {
                io->mod->s_opt->opt_five_branch = 0;
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --append
            //
        case 40: {
            if(optarg == NULL || !strcmp(optarg, "true")) {
                io->writemode = 2;
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --n_rr_branch
            //
        case 39: {
            io->mod->n_rr_branch = (int)atoi(optarg);
            if(io->mod->n_rr_branch < 1) {
                errorMsg = "\n. The number of classes must be an integer greater than 0.\n";
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --no_gap
            //
        case 38: {
            if(optarg == NULL || !strcmp(optarg, "true")) {
                io->rm_ambigu = 1;
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --cov_free
            //
        case 37: {
            if(optarg == NULL || !strcmp(optarg, "true")) {
                io->mod->s_opt->opt_cov_free_rates = 1;
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --ps
            //
        case 36: {
            if(optarg == NULL || !strcmp(optarg, "true")) {
                io->open_ps_file = 1;
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --cov_ncats
            //
        case 35: {
#ifdef M4
            io->m4_model = YES;
            if(!io->mod->m4mod) {
                int ns;
                if(io->datatype == NT)      ns = 4;
                else if(io->datatype == AA) ns = 20;
                else {
                    errorMsg = "\n. Cov_ncats for codons not implemented yet.\n");
                }
                io->mod->m4mod = (m4 *)M4_Make_Light(ns);
            }
            io->mod->m4mod->n_h = (int)atoi(optarg);
            if(io->mod->m4mod->n_h < 1) {
                errorMsg = "\n. The number of classes must be greater than 0.\n";
            }
#endif
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --cov_alpha
            //
        case 34: {
#ifdef M4
            io->m4_model = YES;
            if(!io->mod->m4mod) {
                int ns;
                if(io->datatype == NT)      ns = 4;
                else if(io->datatype == AA) ns = 20;
                else {
                    errorMsg = "\n. Cov-alpha for codons not implemented yet.\n";
                }
                io->mod->m4mod = (m4 *)M4_Make_Light(ns);
            }
            
            io->mod->m4mod->use_cov_alpha = 1;
            io->mod->m4mod->use_cov_free  = 0;
            
            if(!strcmp(optarg, "e") || !strcmp(optarg, "E") ||
               !strcmp(optarg, "estimated") || !strcmp(optarg, "ESTIMATED")) {
                io->mod->s_opt->opt_cov_alpha = 1;
                io->mod->m4mod->alpha         = 1.0;
            } else {
                io->mod->m4mod->alpha = (phydbl)atof(optarg);
                if(io->mod->m4mod->alpha < 1.E-5) {
                    errorMsg("\n. The value of alpha must be greater than 1.E-5.\n");
                }
            }
#endif	      
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --cov_delta
            //
        case 33 : {
#ifdef M4
            io->m4_model = YES;
            if(!io->mod->m4mod) {
                int ns;
                if(io->datatype == NT)      ns = 4;
                else if(io->datatype == AA) ns = 20;
                else {
                    errorMsg = "\n. Cov-delta for codons not implemented yet.\n";
                }
                io->mod->m4mod = (m4 *)M4_Make_Light(ns);
            }
            
            if(!strcmp(optarg, "e") || !strcmp(optarg, "E") ||
               !strcmp(optarg, "estimated") || !strcmp(optarg, "ESTIMATED")) {
                io->mod->s_opt->opt_cov_delta = 1;
                io->mod->m4mod->delta         = 1.0;
            } else {
                io->mod->m4mod->delta = (phydbl)atof(optarg);
                if(atof(optarg) < 1.E-10) {
                    errorMsg = "\n. The value of delta must be larger than 1.E-10.\n";
                }
            }
#endif
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --cov
            //
        case 32 : {
#ifdef M4
            if(optarg == NULL || !strcmp(optarg, "true")) {
                io->m4_model = YES;
            }
#endif	      
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --print_site_lnl
            //
        case 31 : {
            if(optarg == NULL || !strcmp(optarg, "true") || !strcmp(optarg, "1")) {
                io->print_site_lnl = 1;
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --print_trace
            //
        case 30 : {
            if(optarg == NULL || !strcmp(optarg, "true")) {
                io->print_trace = 1;
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --random_boot
            //
        case 29 : {
            io->random_boot_seq_order = (int)atoi(optarg);
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --collapse_boot
            //
        case 28 : {
            io->collapse_boot = (int)atoi(optarg);
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --r_seed
            //
        case 27 : {
            io->r_seed = (int)atoi(optarg);
            srand(io->r_seed);
	    SetSeed(io->r_seed);
	    break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --g_pars
            //
        case 26 : {
            if(optarg == NULL || !strcmp(optarg, "true")) {
                io->mod->s_opt->general_pars = 1;
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --fast_nni
            //
        case 25 : {
            if(optarg == NULL || !strcmp(optarg, "true")) {
                io->mod->s_opt->fast_nni = 1;
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --p_moves
            //
        case 24 : {
            io->mod->s_opt->p_moves_to_examine = (phydbl)atof(optarg);
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --inside_opt
            //
        case 23 : {
            if(optarg == NULL || !strcmp(optarg, "true")) {
                io->mod->s_opt->wim_inside_opt = 1;
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --n_rgrft
            //
        case 0 : {
            io->mod->s_opt->wim_n_rgrft = atoi(optarg);
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --n_globl
            //
        case 1 : {
            io->mod->s_opt->wim_n_globl = atoi(optarg);
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --max_dist
            //
        case 2 : {
            io->mod->s_opt->wim_max_dist = atoi(optarg);
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --n_optim
            //
        case 3 : {
            io->mod->s_opt->wim_n_optim = atoi(optarg);
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --n_best
            //
        case 4 : {
            io->mod->s_opt->wim_n_best = atoi(optarg);
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --min_diff_lk_local
            //
        case 16 : {
            io->mod->s_opt->min_diff_lk_local = atof(optarg);
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --min_diff_lk_global
            //
        case 17 : {
            io->mod->s_opt->min_diff_lk_global = atof(optarg);
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --steph_spr
            //
        case 18 : {
            io->mod->s_opt->steph_spr = 0;
            io->mod->s_opt->greedy    = 1;
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --brent_it_max
            //
        case 19 : {
            io->mod->s_opt->brent_it_max = atoi(optarg);
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --rand_start
            //
        case 20 : {
            if(optarg == NULL || !strcmp(optarg, "true")) {
                io->mod->s_opt->random_input_tree = 1;
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --n_rand_starts
            //
        case 21 : {
            io->mod->s_opt->random_input_tree = 1;
            io->mod->s_opt->n_rand_starts = atoi(optarg);
            if(io->mod->s_opt->n_rand_starts < 1) errorMsg = "\n. Number of random starting trees must be > 0.\n\n";
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --search
            //
        case 's':case 6: {
            if((!strcmp(optarg, "spr")) || (!strcmp(optarg, "SPR"))) {
                io->mod->s_opt->topo_search         = SPR_MOVE;
                io->mod->s_opt->greedy              = (io->mod->s_opt->steph_spr)?(0):(1);
                io->print_trace=1;
            } else if((!strcmp(optarg, "nni")) || (!strcmp(optarg, "NNI"))) {
                io->mod->s_opt->topo_search         = NNI_MOVE;
                io->mod->s_opt->random_input_tree   = 0;
            } else if((!strcmp(optarg, "best")) || (!strcmp(optarg, "BEST"))) {
                io->mod->s_opt->topo_search         = BEST_OF_NNI_AND_SPR;
                io->mod->s_opt->greedy              = (io->mod->s_opt->steph_spr)?(0):(1);
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --datatype
            //
        case 'd':case 7: {
            if(strcmp(optarg, "nt") == 0 || strcmp(optarg, "NT") == 0) {
                io->datatype        = NT;
                io->mod->ns         = 4;
                io->mod->state_len  = 1;
                
                if(
                   (io->modeltypeOpt == LG)        ||
                   (io->modeltypeOpt == WAG)       ||
                   (io->modeltypeOpt == DAYHOFF)   ||
                   (io->modeltypeOpt == JTT)       ||
                   (io->modeltypeOpt == BLOSUM62)  ||
                   (io->modeltypeOpt == MTREV)     ||
                   (io->modeltypeOpt == RTREV)     ||
                   (io->modeltypeOpt == CPREV)     ||
                   (io->modeltypeOpt == DCMUT)     ||
                   (io->modeltypeOpt == VT)        ||
                   (io->modeltypeOpt == MTMAM)     ||
                   (io->modeltypeOpt == MTART)     ||
                   (io->modeltypeOpt == HIVW)      ||
                   (io->modeltypeOpt == HIVB)      ||
                   (io->modeltypeOpt == CUSTOMAA)  ||
                   (io->modeltypeOpt < 0)                //! Added by Louis (codon models)
                   ) {
                    io->modeltypeOpt = HKY85;
                    strcpy(io->mod->modelname, "HKY85\0");
                }
            } else if(strcmp(optarg,"aa") == 0 || strcmp(optarg, "AA") == 0) {
                io->datatype              = AA;
                io->mod->state_len        = 1;
                io->mod->s_opt->opt_kappa = 0;
                io->mod->ns               = 20;
                if(
                   (io->modeltypeOpt == JC69)   ||
                   (io->modeltypeOpt == K80)    ||
                   (io->modeltypeOpt == F81)    ||
                   (io->modeltypeOpt == HKY85)  ||
                   (io->modeltypeOpt == F84)    ||
                   (io->modeltypeOpt == TN93)   ||
                   (io->modeltypeOpt == GTR)    ||
                   (io->modeltypeOpt == CUSTOM) ||
                   (io->modeltypeOpt < 0)              //! Added by Louis (codon models) 
                   ) {
                    io->modeltypeOpt = LG;
                    strcpy(io->mod->modelname, "LG\0");
                }
            } else if((!strcmp(optarg,"generic")) || (!strcmp(optarg,"gen"))) {
                io->datatype = GENERIC;
            } else if(strcmp(optarg,"codon") == 0 || strcmp(optarg, "CODON") == 0) {          //! Still not correct!
                io->datatype              = CODON;
                io->mod->state_len        = 1;
                if((io->modeltypeOpt == JC69)   ||        //! Added by Louis (Otherwise model is arbitrarily changed)
                   (io->modeltypeOpt == K80)    ||
                   (io->modeltypeOpt == F81)    ||
                   (io->modeltypeOpt == HKY85)  ||
                   (io->modeltypeOpt == F84)    ||
                   (io->modeltypeOpt == TN93)   ||
                   (io->modeltypeOpt == GTR)    ||
                   (io->modeltypeOpt == CUSTOM)    ||
                   (io->modeltypeOpt == LG)        ||
                   (io->modeltypeOpt == WAG)       ||
                   (io->modeltypeOpt == DAYHOFF)   ||
                   (io->modeltypeOpt == JTT)       ||
                   (io->modeltypeOpt == BLOSUM62)  ||
                   (io->modeltypeOpt == MTREV)     ||
                   (io->modeltypeOpt == RTREV)     ||
                   (io->modeltypeOpt == CPREV)     ||
                   (io->modeltypeOpt == DCMUT)     ||
                   (io->modeltypeOpt == VT)        ||
                   (io->modeltypeOpt == MTMAM)     ||
                   (io->modeltypeOpt == MTART)     ||
                   (io->modeltypeOpt == HIVW)      ||
                   (io->modeltypeOpt == HIVB)      ||
                   (io->modeltypeOpt == CUSTOMAA)) {
                    io->mod->freq_model =CF3X4; 
                    io->modeltypeOpt = GY;
                    strcpy(io->mod->modelname, "GY\0");
                }
            } else {
                errorMsg = "\n. Unknown argument to -d option: please use `nt' for DNA or `aa' for Amino-Acids\n";
            }
            break;
        }
            
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --model
            //
        case 'm': case 5 : {
            int i;
            
            For(i,strlen(optarg)) Uppercase(optarg+i);
            
            if(!isalpha(optarg[0])) {
                strcpy(io->mod->custom_mod_string,optarg);
                
                if(strlen(io->mod->custom_mod_string) != 6) {
                    errorMsg = "\n. The string should be of length 6.\n";
                } else {
                }
                
                io->datatype              = NT;
                io->modeltypeOpt       = CUSTOM;
                strcpy(io->mod->modelname, "custom");
                io->mod->s_opt->opt_kappa = 0;
                io->mod->s_opt->opt_rr    = 1;
            }
            
            if(strcmp(optarg, "JC69") == 0) {
                io->datatype              = NT;
                io->modeltypeOpt       = JC69;
            } else if(strcmp(optarg, "K80") == 0) {
                io->datatype              = NT;
                io->modeltypeOpt       = K80;
            } else if(strcmp(optarg, "F81") == 0) {
                io->datatype              = NT;
                io->modeltypeOpt       = F81;
            } else if(strcmp(optarg, "HKY85") == 0) {
                io->datatype              = NT;
                io->modeltypeOpt       = HKY85;
            } else if(strcmp(optarg, "F84") == 0) {
                io->datatype              = NT;
                io->modeltypeOpt       = F84;
            } else if(strcmp (optarg,"TN93") == 0) {
                io->datatype              = NT;
                io->modeltypeOpt       = TN93;
            } else if(strcmp (optarg, "GTR") == 0) {
                io->datatype              = NT;
                io->modeltypeOpt       = GTR;
            } else if(strcmp(optarg, "DAYHOFF") == 0) {
                io->datatype              = AA;
                io->modeltypeOpt       = DAYHOFF;
            } else if(strcmp (optarg, "JTT") == 0) {
                io->datatype              = AA;
                io->modeltypeOpt       = JTT;
            } else if(strcmp(optarg, "MTREV") == 0) {
                io->datatype             = AA;
                io->modeltypeOpt      = MTREV;
            } else if(strcmp (optarg, "LG") == 0) {
                io->datatype              = AA;
                io->modeltypeOpt       = LG;
            } else if(strcmp (optarg, "WAG") == 0) {
                io->datatype              = AA;
                io->modeltypeOpt       = WAG;
            } else if(strcmp(optarg, "DCMUT") == 0) {
                io->datatype              = AA;
                io->modeltypeOpt       = DCMUT;
            } else if(strcmp (optarg, "RTREV") == 0) {
                io->datatype              = AA;
                io->modeltypeOpt       = RTREV;
            } else if(strcmp(optarg, "CPREV") == 0) {
                io->datatype              = AA;
                io->modeltypeOpt       = CPREV;
            } else if(strcmp(optarg, "VT") == 0) {
                io->datatype              = AA;
                io->modeltypeOpt       = VT;
            } else if(strcmp(optarg, "BLOSUM62") == 0) {
                io->datatype              = AA;
                io->modeltypeOpt       = BLOSUM62;
            } else if(strcmp(optarg, "MTMAM") == 0) {
                io->datatype              = AA;
                io->modeltypeOpt       = MTMAM;
            } else if(strcmp(optarg,"MTART") == 0) {
                io->datatype              = AA;
                io->modeltypeOpt       = MTART;
            } else if(strcmp(optarg,"HIVW") == 0) {
                io->datatype              = AA;
                io->modeltypeOpt       = HIVW;
            } else if(strcmp(optarg, "HIVB") == 0) {
                io->datatype              = AA;
                io->modeltypeOpt       = HIVB;
            } else if(strcmp(optarg, "CUSTOM") == 0) {
                io->datatype              = AA;
                io->modeltypeOpt       = CUSTOMAA;
            } else if(strcmp(optarg, "GY") == 0) {
                strcpy(io->nt_or_cd, "codons");
                io->datatype              = CODON;
                io->modeltypeOpt          = GY;
            } else if(strcmp(optarg, "MG") == 0) {
                strcpy(io->nt_or_cd, "codons");
                io->datatype              = CODON;       
                io->modeltypeOpt          = MG;
            } else if(strcmp(optarg, "YAP") == 0) {
                strcpy(io->nt_or_cd, "codons");
                io->datatype              = CODON;       
                io->modeltypeOpt          = YAP;
            } else if(strcmp(optarg, "PCM") == 0) {
                strcpy(io->nt_or_cd, "codons");
                io->datatype              = CODON;       
                io->modeltypeOpt          = PCM;
            }else if(strcmp(optarg, "HLP17") == 0) { //set up HLP17 model with default params
                strcpy(io->nt_or_cd, "codons");      //Added by Ken
                io->datatype              = CODON;
                io->modeltypeOpt          = HLP17;
                io->expm				  = SSPADE;
                io->mod->opthotness       = 1; //optimize h
                io->mod->kappa            = 1.0; //optimize K
                io->mod->s_opt->opt_kappa = 1;
                io->userWantsKappa        = YES;
                io->eq_freq_handling      = OPTIMIZE; // optimize eq freqs
                io->mod->freq_model       = CF3X4;
                io->mod->whichmodel       = HLP17;
                io->freqmodelOpt		  = CF3X4;
                //Initialize Bmat
                io->mod->Bmat = (phydbl *)mCalloc(3721,sizeof(phydbl));

                //set up omega (DM0) model
                io->omegaOpt = DM0;
                io->wclasses = 1;
                io->mod->n_w_catg = io->wclasses;
                io->mod->omegaSiteVar = io->omegaOpt;
                io->mod->prob_omegas           = (phydbl *)mCalloc(io->mod->n_w_catg,sizeof(phydbl));
                io->mod->omegas                = (phydbl *)mCalloc(io->mod->n_w_catg,sizeof(phydbl));
                io->mod->omegas[0]               = 1.0;
                io->mod->prob_omegas[0]          = 1.0;
                //io->mod->omega                   = 0.4;
                io->mod->s_opt->opt_omega        = 1;
                io->mod->s_opt->opt_prob_omega   = 0;
                io->mod->s_opt->opt_alpha        = 0;
                io->mod->s_opt->opt_beta         = 0;
            } else {
                Warn_And_Exit("Please check your model string (-m option).");
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --qrates
            // init rate matrices
            //
        case 137 : {
            if(strcmp(optarg, "KOSI07") == 0) {
                io->mod->initqrates   = KOSI07;
            } else if(strcmp(optarg, "SCHN05") == 0) {
                io->mod->initqrates   = SCHN05;
            } else if(strcmp(optarg, "ECMUSR") == 0) {
                io->mod->initqrates   = ECMUSR;
            } else {
                errorMsg = "Init rate matrix must be either KOSI07, SCHN05 or ECMUSR with user defined matrix.";
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --fmodel
            // frequency model: F1XCODONS, F1X4, F3X4, CF3X4, FMODEL
        case 138 : {
            if(strcmp(optarg, "F1XCODONS") == 0) {
                io->freqmodelOpt  = F1XSENSECODONS;
            } else if(strcmp(optarg, "F1X4") == 0) {
                io->freqmodelOpt  = F1X4;
            } else if(strcmp(optarg, "F3X4") == 0) {
                io->freqmodelOpt  = F3X4;
            } else if(strcmp(optarg, "FMODEL") == 0) {
                io->freqmodelOpt  = FMODEL;
            } else { // default value if option is not present
                io->freqmodelOpt  = CF3X4;
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --alpha
            //
        case 'a':case 14 : {
            if((strcmp (optarg, "e") == 0) ||
               (strcmp (optarg, "E") == 0) ||
               (strcmp (optarg, "estimated") == 0) ||
               (strcmp (optarg, "ESTIMATED") == 0)) {
                io->mod->s_opt->opt_alpha     = 1;
                io->mod->s_opt->opt_alphaCD   = 1; //!< Added by Marcelo.
            } else if(atof(optarg) < 1.E-10) {
                errorMsg = "\n. Alpha must be > 1.E-10.\n";
            } else {
                io->mod->alpha = (phydbl)atof(optarg);
                io->mod->s_opt->opt_alpha   = 0;
                io->mod->s_opt->opt_alphaCD = 0; //!< Added by Marcelo.
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --bootstrap
            //
        case 'b':case 10: {
            if (atoi(optarg) < -5) {
                errorMsg = "\n. Branch test value must be a positive integer for bootstrap, or between -1 and -5 for aLRT branch test\n";
            } else {
                if((int)atoi(optarg) > 0) {
                    io->ratio_test       = 0;
                    io->mod->bootstrap   = (int)atoi(optarg);
                    io->print_boot_trees = 1;
                    
                    if(io->n_data_sets > 1) {
                        errorMsg = "\n. Bootstrap option is not allowed with multiple data sets\n";
                    }
                } else if(atoi(optarg)==0) {
                    io->mod->bootstrap = 0;
                    io->ratio_test     = 0;
                } else {
                    io->mod->bootstrap = 0;
                    io->ratio_test     = -(int)atoi(optarg);
                }
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --nclasses
            //
        case 'c':case 12: {
            if((!atoi(optarg)) || (atoi(optarg) < 0)) {
                errorMsg = "\n. Unknown argument to -c option: the number of rate categories must be a positive integer\n";
            } else { 						 
                io->mod->n_catg = atoi(optarg);
                if(io->mod->n_catg < 1) {
                    errorMsg = "\n. The number of rate categories must be a positive integer\n";
                }
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // -f, --frequencies
            //
        case 'f': case 130: {
            if(strcmp(optarg,"empirical") == 0 || 
               strcmp(optarg,"EMPIRICAL") == 0 || 
               strcmp(optarg,"e") == 0 ||
               strcmp(optarg,"E") == 0) {
                io->eq_freq_handling = EMPIRICAL;
            } else if(strcmp(optarg,"model") == 0 ||
                      strcmp(optarg,"MODEL") == 0 || 
                      strcmp(optarg,"m") == 0 ||
                      strcmp(optarg,"M") == 0) {
                io->eq_freq_handling = MODEL;
            } else if(strcmp(optarg,"optimize") == 0 || 
                      strcmp(optarg,"OPTIMIZE") == 0 ||
                      strcmp(optarg,"o") == 0 ||
                      strcmp(optarg,"O") == 0) {
                io->eq_freq_handling = OPTIMIZE;
            } else if(!isalpha(optarg[0])) {
                phydbl sum;
                double val1, val2, val3, val4;
                phydbl freqs[64] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; 
                int i, j, k;
                
                io->mod->s_opt->opt_state_freq  = NO;
                io->mod->s_opt->user_state_freq = YES;
                io->eq_freq_handling = USER;
                
                if(io->datatype == NT) {
                    k=sscanf(optarg,"%lf,%lf,%lf,%lf",&val4,&val2,&val1,&val3); //fT,fC,fA,fG
                    if(k!=4) errorMsg = "Too few nt frequencies provided.\n";
                    io->mod->user_b_freq[0] = (phydbl)val1; 
                    io->mod->user_b_freq[1] = (phydbl)val2; 
                    io->mod->user_b_freq[2] = (phydbl)val3; 
                    io->mod->user_b_freq[3] = (phydbl)val4; 
                    
                    sum = io->mod->user_b_freq[0] +
                    io->mod->user_b_freq[1] +
                    io->mod->user_b_freq[2] +
                    io->mod->user_b_freq[3];
                    
                    io->mod->user_b_freq[0] /= sum;
                    io->mod->user_b_freq[1] /= sum;
                    io->mod->user_b_freq[2] /= sum;
                    io->mod->user_b_freq[3] /= sum;
                    
                    
                    if(io->mod->user_b_freq[0] < .0 ||
                       io->mod->user_b_freq[1] < .0 ||
                       io->mod->user_b_freq[2] < .0 ||
                       io->mod->user_b_freq[3] < .0 ||
                       io->mod->user_b_freq[0] > 1. ||
                       io->mod->user_b_freq[1] > 1. ||
                       io->mod->user_b_freq[2] > 1. ||
                       io->mod->user_b_freq[3] > 1.) {
                        errorMsg = "\n. Invalid base frequencies.\n";
                    }
                } else if(io->datatype == CODON) {
                    for(i=0; optarg[i] != 0; i++) {
                        io->userFreqs[i] = optarg[i];
                        io->eq_freq_handling = USER;
                    }
                } else if(io->datatype==AA) {
                    k=sscanf(optarg,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",&freqs[0],&freqs[1],&freqs[2],&freqs[3],&freqs[4],&freqs[5],&freqs[6],&freqs[7],&freqs[8],&freqs[9],&freqs[10],&freqs[11],&freqs[12],&freqs[13],&freqs[14],&freqs[15],&freqs[16],&freqs[17],&freqs[18],&freqs[19]);
                    if(k != 20) {
                        errorMsg = "Too few aa frequencies provided.\n";
                    }
                    sum = 0.0;
                    For(i, 20) sum += freqs[i];
                    j = -1;
                    For(i, 20) { 
                        freqs[i] /= sum;
                        if(freqs[i] < .0 || freqs[i] > 1.) {
                            errorMsg = "\n. Invalid base frequencies.\n";
                        }
                        io->mod->user_b_freq[i]=freqs[i];
                    }
                }
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // -h, --help
            //
        case 'h': case 134: {
            Usage();
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // -i, --input
            //
        case 'i':case 9: {
            tmp = (char *) mCalloc (T_MAX_FILE, sizeof(char));
            if(strlen (optarg) > T_MAX_FILE -16) {
                strcpy(tmp, "\n. The file name'");
                strcat(tmp, optarg);
                strcat(tmp, "' is too long.\n");
                errorMsg = tmp;
            } else if(!Filexists(optarg)) {
                strcpy(tmp, "\n. The file '");
                strcat(tmp, optarg);
                strcat(tmp, "' does not exist.\n");
                errorMsg = tmp;
            } else {
                strcpy(io->in_align_file, optarg);
                io->fp_in_align = Openfile(io->in_align_file, 0);
                strcpy(io->out_tree_file, optarg);
#ifdef PHYML
                strcat(io->out_tree_file, "_igphyml_tree.txt");
#else
                strcat(io->out_tree_file, "_mc_tree.txt");
#endif
                strcpy(io->out_stats_file, optarg);
#ifdef PHYML
                strcat(io->out_stats_file, "_igphyml_stats.txt");
#else
                strcat(io->out_stats_file, "_mc_stats.txt");
#endif
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // -t, --ts/tv
            //
        case 't': case  11: {        
            if(io->modeltypeOpt < 0) { // if its a GY, MG or YAP model - PCM gets fixed later, no worries
                if ((strcmp(optarg, "e") == 0)         ||
                    (strcmp(optarg, "E") == 0)         ||
                    (strcmp(optarg, "estimated") == 0) ||
                    (strcmp(optarg, "ESTIMATED") == 0)) {
                    io->mod->kappa                 = 1.0;
                    io->mod->s_opt->opt_kappa      = 1;
                    io->userWantsKappa             = YES;
                } else if(strcmp(optarg, "KAP1") == 0 || strcmp(optarg,"kap1") == 0) { //!< For all KAPn, See Kosiol 2007.
                    io->kappaECM              = kap1;
                    io->mod->s_opt->opt_kappa = 0;
                    io->mod->nkappa           = 1;
                    io->mod->kappa            = 1;
                    io->mod->pkappa[0]        = 1;
                } else if(strcmp(optarg, "KAP2") == 0 || strcmp(optarg,"kap2") == 0) {
                    io->kappaECM              = kap2;
                    io->mod->s_opt->opt_kappa = 1;
                    io->userWantsKappa        = YES;
                    io->mod->nkappa           = 1;
                    io->mod->pkappa[0]        = 1;
                } else if(strcmp(optarg, "KAP3") == 0 || strcmp(optarg,"kap3") == 0) {
                    io->kappaECM              = kap3;
                    io->mod->s_opt->opt_kappa = 1;
                    io->userWantsKappa        = YES;
                    io->mod->nkappa           = 1;
                    io->mod->pkappa[0]        = 1;
                } else if(strcmp(optarg, "KAP4") == 0 || strcmp(optarg,"kap4") == 0) {
                    io->kappaECM              = kap4;
                    io->mod->s_opt->opt_kappa = 1;
                    io->userWantsKappa        = YES;
                    io->mod->nkappa           = 2;
                    io->mod->pkappa[0]        = 1;
                    io->mod->pkappa[1]        = 1;
                } else if(strcmp(optarg, "KAP5") == 0 || strcmp(optarg,"kap5") == 0) {
                    io->kappaECM              = kap5;
                    io->mod->s_opt->opt_kappa = 1;
                    io->userWantsKappa        = YES;
                    io->mod->nkappa           = 9;
                    int i;                                            
                    For(i, 9) io->mod->pkappa[i] = 1.0/9.0;              
                } else if(strcmp(optarg, "KAP6") == 0 || strcmp(optarg,"kap6") == 0) {
                    io->kappaECM              = kap6;
                    io->mod->s_opt->opt_kappa = 1;
                    io->userWantsKappa        = YES;
                    io->mod->nkappa           = 1;
                    io->mod->pkappa[0]        = 1.0;              
                } else if (atof(optarg) < .0) {
                    errorMsg =  "\n. The ts/tv ratio must be a positive number\n";
                } else {
                    io->mod->kappa = (phydbl)atof(optarg);
                    io->mod->s_opt->opt_kappa  = 0;
                    io->mod->s_opt->opt_lambda = 0;
                }
            } else if((io->modeltypeOpt != JC69) &&
                      (io->modeltypeOpt != F81)  &&
                      (io->modeltypeOpt != GTR)) {
                if((strcmp(optarg, "e") == 0)         ||
                   (strcmp(optarg, "E") == 0)         ||
                   (strcmp(optarg, "estimated") == 0) ||
                   (strcmp(optarg, "ESTIMATED") == 0)) {
                    io->mod->kappa                 = 1.0;
                    io->mod->s_opt->opt_kappa      = 1;
                    io->userWantsKappa        = YES;
                    if (io->modeltypeOpt == TN93)
                        io->mod->s_opt->opt_lambda   = 1;
                } else {
                    if (atof(optarg) < .0) {
                        errorMsg =  "\n. The ts/tv ratio must be a positive number\n";
                    } else {
                        io->mod->kappa = (phydbl)atof(optarg);
                        io->mod->s_opt->opt_kappa  = 0;
                        io->mod->s_opt->opt_lambda = 0;
                    }
                }
            }
            break;
        }

        //Added by Kenneth Hoehn 3/6/2016
        //--hotness
        case 147: {
        	strcpy(io->mod->hotnessstring,optarg);
        	io->mod->hotnessstringopt=1;
        	break;
        }
        //--root
        case 148: {
           strcpy(io->mod->rootname,optarg);
           io->mod->rootfound=1;
           break;
        }
        //--motifs
        //1 split first by , and then by : to get information on the
        //type, hotspot position, and number of h parameters
        //strand symmetry e.g. WRC_2:0,GYW_0:0
        //strand asymetry e.g. WRC_2:0,GYW_0:1
        case 149: {
        	strcpy(io->mod->motifstring,optarg);
        	io->mod->motifstringopt=1;
        	break;
        }
        case 150: {
           	strcpy(io->mod->partfile,optarg);
           	io->mod->partfilespec=1;
          	break;
        }
        case 151: {
           	strcpy(io->mod->ambigfile,optarg);
           	io->mod->ambigprint=1;
          	break;
        }
        case 152: {
           	int threads = atoi(strdup(optarg));

           	#if defined OMP || defined BLAS_OMP
           	omp_set_dynamic(0);
           	omp_set_num_threads(threads);
			#else
           	printf("\n. Can't specify number of threads unless compiled with OMP!\n");
           	exit(EXIT_FAILURE);
			#endif
           	break;
        }
        case 153: {
        	io->mod->slowSPR=1;
        	break;
        }
        case 154: {
           	io->mod->stretch=atof(strdup(optarg));
           	printf("stretch by %lf\n",io->mod->stretch);
           	break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --multiple
            //
        case 'n': case 8: {
            if((!atoi(optarg)) || (atoi(optarg) < 0)) {
                errorMsg = "\n. The number of alignments must be a positive integer\n";
            }
            else io->n_data_sets = atoi (optarg);
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --sequential
            //
        case 'q':case 22: {
            if(optarg == NULL || !strcmp(optarg, "true")) {
                io->interleaved = 0;
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --inputtree
            //
        case 'u': case 15: {
            tmp = (char *)mCalloc(T_MAX_FILE, sizeof(char));
            if(strlen(optarg) > (T_MAX_FILE - 11)) {
                strcpy (tmp, "\n. The file name'");
                strcat (tmp, optarg);
                strcat (tmp, "' is too long.\n");
                errorMsg = tmp;
            } else if(!Filexists (optarg)) {
                strcpy (tmp, "\n. The file '");
                strcat (tmp, optarg);
                strcat (tmp, "' doesn't exist.\n");
                errorMsg = tmp;
            } else {
                strcpy(io->in_tree_file, optarg);
                io->in_tree = 2;
                io->fp_in_tree = Openfile(io->in_tree_file,0);
            }
            break;
        }
            
            //////////////////////////////////////////////////////////////////////////////////////
            // oformat
            // output format for stats file
            //
        case 142: {
#ifdef USEYAML
            if ((strcmp (optarg, "yaml") == 0) ||
                (strcmp (optarg, "yml") == 0) ||
                (strcmp (optarg, "YAML") == 0) ||
                (strcmp (optarg, "YML") == 0)) {
                io->out_stats_format = OUTYAML;
            } else
#endif
            if ((strcmp (optarg, "darwin") == 0) ||
                       (strcmp (optarg, "Darwin") == 0)) {
                io->out_stats_format = OUTDARWIN;
            } else {
                io->out_stats_format = OUTTXT;
            }
            break;
        }

        case 146: {
            io->logtree = atoi(optarg);
            if(io->logtree > 2 || io->logtree < 0) {
                Warn_And_Exit("logtree option has to be either 0, 1 or 2");
            }
            break;
        }

            
            //////////////////////////////////////////////////////////////////////////////////////
            // --pinv
            //
        case 'v':case 13: {
            if ((strcmp (optarg, "e") == 0) ||
                (strcmp (optarg, "E") == 0) ||
                (strcmp (optarg, "estimated") == 0) ||
                (strcmp (optarg, "ESTIMATED") == 0)) {
                io->mod->s_opt->opt_pinvar    = 1;
                io->mod->invar                = 1;
            } else if((atof(optarg) < 0.0) || (atof(optarg) > 1.0)) {
                errorMsg = "\n. The proportion of invariable sites must be a number between 0.0 and 1.0\n";
            } else {
                io->mod->pinvar = (phydbl)atof(optarg);
                if(io->mod->pinvar > 0.0+SMALL) {
                    io->mod->invar = 1;
                } else {
                    io->mod->invar = 0;
                }
                io->mod->s_opt->opt_pinvar = 0;
            }
            break;
        }
            
            
            //////////////////////////////////////////////////////////////////////////////////////
            // --optimize
            //
        case 'o': case 131: {
            if(!strcmp(optarg, "tlr")) {
                io->mod->s_opt->opt_topo        = 1;
                io->mod->s_opt->opt_bl          = 1;
                io->mod->s_opt->opt_subst_param = 1;
            } else if(!strcmp(optarg, "tl")) {
                io->mod->s_opt->opt_topo        = 1;
                io->mod->s_opt->opt_bl          = 1;
                io->mod->s_opt->opt_subst_param = 0;
            } else if(!strcmp(optarg, "t")) {
                errorMsg = "\n. You can't optimize the topology without adjusting branch length too...\n";
            } else if(!strcmp(optarg, "lr")) {
                io->mod->s_opt->opt_topo        = 0;
                io->mod->s_opt->opt_bl          = 1;
                io->mod->s_opt->opt_subst_param = 1;
            } else if(!strcmp(optarg, "l")) {
                io->mod->s_opt->opt_topo        = 0;
                io->mod->s_opt->opt_bl          = 1;
                io->mod->s_opt->opt_subst_param = 0;
            } else if(!strcmp(optarg, "r")) {
                io->mod->s_opt->opt_topo        = 0;
                io->mod->s_opt->opt_bl          = 0;
                io->mod->s_opt->opt_subst_param = 1;
            } else if(!strcmp(optarg, "none") || !strcmp(optarg, "n")) {
                io->mod->s_opt->opt_topo        = 0;
                io->mod->s_opt->opt_bl          = 0;
                io->mod->s_opt->opt_subst_param = 0;
            } else {
                errorMsg = "\n. The optimization parameter must be 'tlr' or 'tl' or 'lr' or 'l' or 'r' or ''.\n";
            }
            break;
        }
            
        default:
            result = -1;
            break;
    }
    
    if(errorMsg != NULL) {
        Warn_And_Exit(errorMsg);
    }
    
    if(tmp != NULL) {
        free(tmp);
    }
    
    return(result);
}
