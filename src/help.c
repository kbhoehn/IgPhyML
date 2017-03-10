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


/* Help that needs to be checked: 
 --ts/tv
 --frequencies
 --model
 --nclasses
 
 */


#include "help.h"


/* int  T_MAX_FILE; */
/* phydbl SMALL; */
/* phydbl UNLIKELY; */

/*********************************************************/

void Usage()
{
    
    char *BOLD=(char *)mCalloc(10,sizeof(char));
    char *FLAT=(char *)mCalloc(10,sizeof(char));
    char *LINE=(char *)mCalloc(10,sizeof(char));
    char *cha;
    
    
    cha =getenv("OS");
    
    if(cha!=NULL) 
    {
        strcpy(BOLD, "");
        strcpy(FLAT, "");
        strcpy(LINE, "");
    } 
    else 
    {
        strcpy(BOLD, "\033[00;01m");
        strcpy(FLAT, "\033[00;00m");
        strcpy(LINE, "\033[00;04m");
    }
    
    //! Modified by Louis and Marcelo 25.10.2012
    PhyML_Printf("%sNAME\n"
    			 "%s\t- IgPhyML %s - \n\n"
    		     "%s\tKenneth B. Hoehn, Gerton Lunter, Oliver G. Pybus.\n\n"
                 "%s\t- CodonPhyML %s - \n\n"
                 "%s\tMarcelo S Zanetti, Manuel Gil, Stefan Zoller, Louis Du Plessis and Maria Anisimova.\n\n"
                 "%s\t- PhyML %s - \n\n"
                 "%s\tStephane Guindon and Olivier Gascuel,\n"
                 "%s\tSystematic Biology 52(5):696-704, 2003.\n\n"
                 "%s\tPlease cite this paper if you use this software in your publications.\n",BOLD,FLAT,VERSION,FLAT,FLAT,VERSION,FLAT,FLAT,FLAT,FLAT,FLAT);
    
    PhyML_Printf("%s\nSYNOPSIS:\n\n"
                 "%s\tigphyml %s[command args]\n",BOLD,BOLD,BOLD);
    PhyML_Printf("%s\n\tAll the options below are optional (except '%s-i%s' if you want to use the command-line interface, \n",FLAT,BOLD,FLAT); 
    PhyML_Printf("%s\n\nNOTE: USE OF OPTIONS OTHER THAN THOSE SPECIFIED IN THE IGPHYML MANUAL IS NOT CURRENTLY SUPPORTED.\n\n",FLAT,BOLD,FLAT);
    
    PhyML_Printf("%s\nCOMMAND OPTIONS:\n%s",BOLD,FLAT);
    
    PhyML_Printf("\n\t%sCommon options:%s\n",BOLD,FLAT);
    
    PhyML_Printf("\n\t%s-i (or --input) %sseq_file_name%s (Required)%s\n",BOLD,LINE,BOLD,FLAT); 
    PhyML_Printf("\t\t%sseq_file_name%s is the name of the nucleotide or amino-acid sequence file in PHYLIP format.\n",LINE,FLAT);  //! Modified by Louis
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t-q (or --sequential)\n",BOLD);
    PhyML_Printf("%s\t\tChanges interleaved format (default) to sequential format.\n",FLAT);
    
    PhyML_Printf("\n");
    
   /* PhyML_Printf("%s\n\t-b (or --bootstrap) %sint%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%sint%s >  0: %sint%s is the number of bootstrap replicates.\n",LINE,FLAT,LINE,FLAT);
    PhyML_Printf("\t\t%sint%s =  0: neither approximate likelihood ratio test nor bootstrap values are computed.\n",LINE,FLAT);
    PhyML_Printf("\t\t%sint%s = -1: approximate likelihood ratio test returning aLRT statistics.\n",LINE,FLAT);
    PhyML_Printf("\t\t%sint%s = -2: approximate likelihood ratio test returning Chi2-based parametric branch supports.\n",LINE,FLAT);
    //   PhyML_Printf("\t\t%sint%s = -3 : minimum of Chi2-based parametric and SH-like branch supports.\n",LINE,FLAT);
    PhyML_Printf("\t\t%sint%s = -4: (default) SH-like branch supports alone.\n",LINE,FLAT);
    PhyML_Printf("\t\t%sint%s = -5: approximate Bayes (aBayes) support.\n",LINE,FLAT);
    PhyML_Printf("\n");*/
    
    PhyML_Printf("%s\n\t-m (or --model) %smodel%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%smodel%s : substitution model name.\n",LINE,FLAT);
    PhyML_Printf("\n");
    PhyML_Printf("\t\t%s- %sCodon%s based models:%s HLP17%s | %sGY%s (default)\n",FLAT,LINE,FLAT,LINE,FLAT,
                 LINE,FLAT,
                 LINE,FLAT);
    
  /*  PhyML_Printf("\n");
    PhyML_Printf("%s\n\t--qrates %sratematrix%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\tThis option sets an initial rate matrix to create semiempirical codon models.\n");
    PhyML_Printf("\t\tIt only has an effect if a codon model is used.\n"); 
    PhyML_Printf("\t\t%sratematrix%s = %sKOSI07%s | %sSCHN05%s | %sECMUSR%s\n",LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT);*/

//    PhyML_Printf("\n");
//    PhyML_Printf("\tFor more information on model definition see the manual\n");
//    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t-f (or --frequencies) %sempirical%s, %smodel%s, %soptimized%s, %sfT,fC,fA,fG%s,\n\t\t%sfT1,fC1,fA1,fG1,fT2,fC2,fA2,fG2,fT3,fC3,fA3,fG3%s\n",
                 BOLD,LINE,BOLD,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT);
    PhyML_Printf("\t                       %sor%s %sfC1,fC2, ... ,fC64%s\n",BOLD,FLAT,LINE,FLAT);
    
    PhyML_Printf("\t\t%sempirical%s : (default) the character frequencies are determined as follows: \n",LINE,FLAT);
    PhyML_Printf("%s\t\t- %sCodon%s sequences     : (Empirical) the equilibrium codon frequencies are estimated by counting\n"
                 "\t\t  the occurence of bases or codons in the alignment according to the frequency model that is selected.\n",FLAT,LINE,FLAT);
    PhyML_Printf("\n");
    /*PhyML_Printf("\t\t%smodel%s : the character frequencies are determined as follows: \n",LINE,FLAT);
    PhyML_Printf("%s\t\t- %sAmino-acid%s sequences: (Model) the equilibrium amino-acid frequencies are estimated using\n"
                 "\t\t  the frequencies defined by the substitution model.\n",FLAT,LINE,FLAT);
    PhyML_Printf("%s\t\t- %sCodon%s sequences     : (Model) the equilibrium codon frequencies are estimated using\n"
                 "\t\t  the frequencies defined by the substitution model.\n",FLAT,LINE,FLAT);
    PhyML_Printf("\n");		*/
    PhyML_Printf("\t\t%soptimize%s : the character frequencies are determined as follows: \n",LINE,FLAT);
    PhyML_Printf("%s\t\t- %sNucleotide%s sequences: (ML) the equilibrium base frequencies are estimated using maximum likelihood \n",FLAT,LINE,FLAT);
    PhyML_Printf("%s\t\t- %sAmino-acid%s sequences: (ML) the equilibrium amino-acid frequencies are estimated using maximum likelihood\n",FLAT,LINE,FLAT);
    PhyML_Printf("%s\t\t- %sCodon%s sequences     : (ML) the equilibrium codon frequencies are estimated using maximum likelihood\n",FLAT,LINE,FLAT);
    
    /*PhyML_Printf("\n");
    PhyML_Printf("\t\t%s\"fT,fC,fA,fG\"%s : only valid for codon models using the F1X4 frequency model. \n",LINE,FLAT);
    PhyML_Printf("\t\t  fA, fC, fG and fT are floating point numbers that correspond to the frequencies of A, C, G and T \n\t\t  respectively \n");
    PhyML_Printf("\t\t%s\"fT1,fC1,fA1,fG1,fT2,fC2,fA2,fG2,fT3,fC3,fA3,fG3\"%s : only valid for codon models using the F3X4 or\n",LINE,FLAT);
    PhyML_Printf("\t\t  CF3X4 frequency models.  In this case the numbers correspond to the frequencies of bases at different\n");
    PhyML_Printf("\t\t  positions in the codons.\n");
    PhyML_Printf("\t\t%s\"fC1,fC2, ... ,fC64\"%s : only valid for codon models using the F1XCODONS frequency model. In this case \n",LINE,FLAT);
    PhyML_Printf("\t\t  the numbers correspond to the frequencies of the different codons.\n");
    PhyML_Printf("\t\t(WARNING: do not use any blank space between your values of nucleotide, codon or amino acid frequencies, only commas!)\n");*/
    
    
    PhyML_Printf("\n");
    PhyML_Printf("%s\n\t--fmodel %sfrequency model%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\tWhich frequency model to use.\n");
    PhyML_Printf("\t\t%sfrequency model%s = %sF1XCODONS%s | %sF1X4%s | %sF3X4%s | %sCF3X4%s (default)\n",LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT);
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t-t (or --ts/tv) %sts/tv_ratio%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\tThis option sets the transition/transversion ratio.\n");
    PhyML_Printf("\t\t%sts/tv_ratio%s > 0: Set transition/transversion ratio to the value.\n",LINE,FLAT);
    PhyML_Printf("\t\t%sts/tv_ratio%s = e: Get maximum likelihood estimate of transition/transversion ratio (default under HLP17).\n",LINE,FLAT);
    PhyML_Printf("\t\t%sts/tv_ratio%s = 1: Fix estimate (default under GY94).\n",LINE,FLAT);
    //PhyML_Printf("\n");
    /*PhyML_Printf("\t\tIf the model is a semi-parametric codon model the parameter refers to the kappa(i,j) in Kosiol et al 2007:\n");
    PhyML_Printf("\t\t%sts/tv_ratio%s = KAP1: kappa(i,j) = 1\n",LINE,FLAT);
    PhyML_Printf("\t\t%sts/tv_ratio%s = KAP2: kappa(i,j) = K(nts) where nts is the number of transitions between codon i and j.\n",LINE,FLAT);
    PhyML_Printf("\t\t%sts/tv_ratio%s = KAP3: kappa(i,j) = K(ntv) where ntv is the number of transversions between codon i and j.\n",LINE,FLAT);
    PhyML_Printf("\t\t%sts/tv_ratio%s = KAP4: kappa(i,j) = K1(nts)K2(ntv) that is two different parameters.\n",LINE,FLAT);
    PhyML_Printf("\t\t%sts/tv_ratio%s = KAP5: kappa(i,j) = K(k) where k is one of 9 ts/tv combinations. \n",LINE,FLAT);*/
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t-w (or --omega) %smodel%s (Required if a codon model other than HLP17 is used)%s\n",BOLD,LINE,BOLD,FLAT);
    PhyML_Printf("\t\tThe omega parameter or nonsynonymous/synonymous rate ratio.\n");
    PhyML_Printf("\n");
    PhyML_Printf("\t\t%smodel%s = DM0   : (default) Single omega model.\n",LINE,FLAT);
   /* PhyML_Printf("\t\t%smodel%s = DMODEL: Discrete unconstrained distribution model.\n",LINE,FLAT);
    PhyML_Printf("\t\t%smodel%s = DGAMMA: Discrete gamma model.\n",LINE,FLAT);*/
    
    //PhyML_Printf("\n");
    
    /*PhyML_Printf("%s\n\t--wvals %svalues%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%svalues%s are the values for each omega category, separated by commas if more than one.\n",LINE,FLAT);
    PhyML_Printf("\t\tIf this option is set, omega values are not estimated by ML but fixed to these values.\n",LINE,FLAT);
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t--wclasses %snr of classes%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\tDefines the number of omega rate categories if DMODEL or DGAMMA is chosen as omega model.\n",LINE,FLAT);
    PhyML_Printf("\t\tDefaults to 4 if DMODEL or DGAMMA is chosen.\n",LINE,FLAT);*/
    
    /*PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t-v (or --pinv) %sprop_invar%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%sprop_invar%s : proportion of invariable sites.\n",LINE,FLAT);
    PhyML_Printf("\t\tCan be a fixed value in the [0,1] range or %se%s to get the maximum likelihood estimate (default = 0.0).\n",LINE,FLAT);
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t-c (or --nclasses) %snb_subst_cat%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%snb_subst_cat%s : number of relative substitution rate categories.\n",LINE,FLAT);
    PhyML_Printf("\t\tMust be a positive integer  (default = 4).\n");
    PhyML_Printf("\t\t(WARNING: do not use this option with codon models unless with constant selection pressure [omega])\n");
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t-a (or --alpha) %salpha%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%salpha%s : distribution of the gamma distribution shape parameter.\n",LINE,FLAT);
    PhyML_Printf("\t\tCan be a fixed positive value or %se%s to get the maximum likelihood estimate. (default = e)\n",LINE,FLAT);
    
    PhyML_Printf("\n");*/
    
    PhyML_Printf("%s\n\t-s (or --search) %smove%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\tTree topology search operation option.\n");
    PhyML_Printf("\t\tCan be either %sNNI%s (default, fast) or %sSPR%s (a bit slower than NNI) or %sBEST%s (best of NNI and SPR search).\n",LINE,FLAT,LINE,FLAT,LINE,FLAT);
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t-u (or --inputtree) %suser_tree_file%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%suser_tree_file%s : starting tree filename. The tree must be in Newick format.\n",LINE,FLAT);
    
    PhyML_Printf("\n");
    
    /*PhyML_Printf("%s\n\t--dist_tree_model %smodel%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\tThis option sets the codon model that is used to build the initial distance tree.\n");
    PhyML_Printf("\t\tThis option only has an effect if a codon model is used and no starting tree is provided.\n"); 
    PhyML_Printf("\t\t%smodel%s = %sGYECMS05%s | %sGYECMK07%s (default) | %sJC69%s\n",LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT);*/
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t-o (or --optimize) %sparams%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\tThis option focuses on specific parameter optimisation.\n");
    PhyML_Printf("\t\t%sparams%s = tlr : (default) tree topology (t), branch length (l) and rate parameters (r) are optimised.\n",LINE,FLAT);
    PhyML_Printf("\t\t%sparams%s = tl  : tree topology and branch length are optimised.\n",LINE,FLAT);
    PhyML_Printf("\t\t%sparams%s = lr  : branch length and rate parameters are optimised.\n",LINE,FLAT);
    PhyML_Printf("\t\t%sparams%s = l   : branch length are optimised.\n",LINE,FLAT);
    PhyML_Printf("\t\t%sparams%s = r   : rate parameters are optimised.\n",LINE,FLAT);
    PhyML_Printf("\t\t%sparams%s = n   : no parameter is optimised.\n",LINE,FLAT);
    
    PhyML_Printf("\n");
    
   /* PhyML_Printf("%s\n\t--rand_start%s\n",BOLD,FLAT);
    PhyML_Printf("\t\tThis option sets the initial tree to random.\n");
    PhyML_Printf("\t\tIt is only valid if SPR searches are to be performed.\n");
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t--n_rand_starts %snum%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%snum%s is the number of initial random trees to be used.\n",LINE,FLAT);
    PhyML_Printf("\t\tIt is only valid if SPR searches are to be performed (default = 5).\n");
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t--r_seed %snum%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%snum%s is the seed used to initiate the random number generator.\n",LINE,FLAT);
    PhyML_Printf("\t\tMust be an integer.\n");
    
    PhyML_Printf("\n");*/
    
    PhyML_Printf("%s\n\t--run_id %sID_string%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%sAppend the string %sID_string%s at the end of each PhyML output file.\n",FLAT,LINE,FLAT);
    PhyML_Printf("\t\t%sThis option may be useful when running simulations involving PhyML.\n",FLAT);
    
    PhyML_Printf("\n");
    
    /*PhyML_Printf("%s\n\t--quiet%s\n",BOLD,FLAT);
    PhyML_Printf("\t\t%sNo interactive questions (for running in batch mode).\n",FLAT);
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t--version%s\n",BOLD,FLAT);
    PhyML_Printf("\t\t%sPrint the version and exit.\n",FLAT);
    
    PhyML_Printf("\n");*/
    
    PhyML_Printf("%s\n\t-h (or --help)%s\n",BOLD,FLAT);
    PhyML_Printf("\t\t%sShow this help message and exit.\n",FLAT);

    PhyML_Printf("\n");
/*
    PhyML_Printf("%s\n\t--oformat %soption%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%sDefines the output format. Possible values are:\n", FLAT);   
    PhyML_Printf("\t\t%stxt%s (text format, default)\n", LINE, FLAT);
    PhyML_Printf("\t\t%sdarwin%s (darwin format)\n", LINE, FLAT);
    PhyML_Printf("\t\t%syaml%s (YAML output)\n", LINE, FLAT);
    PhyML_Printf("\t\tYAML output can be chosen independent of compiling codonphyml with libyaml.\n",FLAT);

    PhyML_Printf("\n");

    PhyML_Printf("%s\n\t--darwinconfig %sfilepath%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%sPath to an optional configuration file in darwin format.\n", FLAT);
    PhyML_Printf("\t\t%sIf this is set, no other command line option is parsed.\n",FLAT);
    
    PhyML_Printf("\n");

    PhyML_Printf("%s\n\t--yamlconfig %sfilepath%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%sPath to an optional configuration file in YAML format.\n", FLAT);
    PhyML_Printf("\t\t%sIf this is set, no other command line option is parsed.\n", FLAT);
    PhyML_Printf("\t\t%sOnly available if codonphyml has been compiled with libyaml support.\n",BOLD);
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t--logtree %soption%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%soption%s can be 1 (no tree log), 1 (log current tree to disc) or 2 (log all 'best trees' to disc).\n",LINE,FLAT);
    
    PhyML_Printf("\n");*/
    
    PhyML_Printf("\n");
    
    PhyML_Printf("\n\t%sNot so common options:%s\n",BOLD,FLAT);
    
    /*PhyML_Printf("%s\n\t-p (or --pars)%s\n",BOLD,FLAT);
    PhyML_Printf("%s\t\tUse a minimum parsimony starting tree. This option is taken into account when the '-u' option\n",FLAT);
    PhyML_Printf("%s\t\tis absent and when tree topology modifications are to be done.\n",FLAT);
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t--pars_thresh %sthreshold%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%sthreshold%s is the parsimony threshold.\n",LINE,FLAT);
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t--use_median%s\n",BOLD,FLAT);
    PhyML_Printf("\t\t%sThe middle of each substitution rate class in the discrete gamma distribution\n",FLAT);
    PhyML_Printf("\t\t%sis taken as the median.  The mean is used by default.\n",FLAT);
    
    PhyML_Printf("\n");*/
    
    PhyML_Printf("%s\n\t--expm %smethod%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\tThis option sets the method to use for matrix exponentiation.\n");
    PhyML_Printf("\t\t%smethod%s = EIGEN : (default for GY94) Use eigenvalue decomposition (preferable for reversible models).\n",LINE,FLAT);
    PhyML_Printf("\t\t%smethod%s = SSPADE: (default for HLP17) Pade approximation.\n",LINE,FLAT);
    //PhyML_Printf("\t\t%smethod%s = TAYLOR: Use Taylor expansion (fast, but inaccurate).\n",LINE,FLAT);
    
    /*PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t--testCond%s\n",BOLD,FLAT);
    PhyML_Printf("\t\tCalculate the condition of the matrix\n");
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t-g (or --genetic_code) %scode%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%sThis option sets the genetic code to use.  This only has an effect for codon models and\n",FLAT);
    PhyML_Printf("\t\t%salso during translation of NT data into AA data. Kosiol 2007 and Scheneider 2005 based\n",FLAT); 
    PhyML_Printf("\t\tmodels will require the standard genetic code.\n",FLAT);
    PhyML_Printf("\t\t%scode%s = STANDARD : The standard genetic code is used.\n",LINE,FLAT);
    PhyML_Printf("\t\t%scode%s = TVMC : Vertebrate Mitochondrial.\n",LINE,FLAT);
    PhyML_Printf("\t\t%scode%s = TYMC : Yeast Mitochondrial.\n",LINE,FLAT);
    PhyML_Printf("\t\t%scode%s = THMPCMCMSC : Mold, Protozoan, and Coelenterate Mit. Code and Myco/Spiroplasma.\n",LINE,FLAT);
    PhyML_Printf("\t\t%scode%s = THIMC : Invertebrate Mitochondrial.\n",LINE,FLAT);
    PhyML_Printf("\t\t%scode%s = THCDHNC : Ciliate, Dasycladacean and Hexamita Nuclear.\n",LINE,FLAT);
    PhyML_Printf("\t\t%scode%s = THEFMC : Echinoderm and Flatworm Mitochondrial.\n",LINE,FLAT);
    PhyML_Printf("\t\t%scode%s = THENC : Euplotid Nuclear.\n",LINE,FLAT);
    PhyML_Printf("\t\t%scode%s = THBAPPC : Bacterial, Archaeal and Plant Plastid.\n",LINE,FLAT);
    PhyML_Printf("\t\t%scode%s = THAYNC : Alternative Yeast Nuclear.\n",LINE,FLAT);
    PhyML_Printf("\t\t%scode%s = THAMC : Ascidian Mitochondrial.\n",LINE,FLAT);
    PhyML_Printf("\t\t%scode%s = THAFMC : Alternative Flatworm Mitochondrial.\n",LINE,FLAT);
    PhyML_Printf("\t\t%scode%s = BLNC : Blepharisma Nuclear.\n",LINE,FLAT);
    PhyML_Printf("\t\t%scode%s = CHMC : Chlorophycean Mitochondrial.\n",LINE,FLAT);
    PhyML_Printf("\t\t%scode%s = TRMC : Trematode Mitochondrial.\n",LINE,FLAT);
    PhyML_Printf("\t\t%scode%s = SCOMC : Scenedesmus obliquus mitochondrial.\n",LINE,FLAT);
    PhyML_Printf("\t\t%scode%s = THMC : Thraustochytrium Mitochondrial.\n",LINE,FLAT);*/
    
    PhyML_Printf("\n");
    
    
    
    /*PhyML_Printf("%s\n\t--NT2AA%s\n",BOLD,FLAT);
    PhyML_Printf("\t\tConvert nucleotide input sequences to amino acid sequences.\n");
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t--print_site_lnl%s\n",BOLD,FLAT);
    PhyML_Printf("\t\t%sPrint the likelihood for each site in file %s*_phyml_lk.txt%s.\n\t\tNot implemented for codon models (yet).\n",FLAT,LINE,FLAT);
    
    PhyML_Printf("\n");*/
    
    PhyML_Printf("%s\n\t--print_trace%s\n",BOLD,FLAT);
    PhyML_Printf("\t\t%sPrint each phylogeny explored during the tree search process\n",FLAT);
    PhyML_Printf("\t\t%sin file %s*_phyml_trace.txt%s.\n",FLAT,LINE,FLAT);
    
    /*PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t--append%s\n",BOLD,FLAT);
    PhyML_Printf("\t\tAppend to files that already exist, instead of rewriting.\n");
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t--compare_runs %sfilename%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%sfilename%s specifies a file to which the results of the simulation are appended in a condensed\n",LINE,FLAT);
    PhyML_Printf("\t\tformat.  This option is identical to %s--append%s except that the results are stored in a different\n",BOLD,FLAT);
    PhyML_Printf("\t\tfile and in a condensed format on one line with fields separated by semi-colons.\n");
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t--test_init_tree%s\n",BOLD,FLAT);
    PhyML_Printf("\t\tOnly calculate the likelihood of the initial tree, and exit.  This option is identical to \n");
    PhyML_Printf("\t\t%s-o n%s, except that branch supports are not calculated (which makes it faster).\n",BOLD,FLAT);
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t--n_rgrft %snr_pos%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%snr_pos%s is the number of promising regraft positions to consider when\n",LINE,FLAT);
    PhyML_Printf("\t\tperforming all improving SPR moves (default = 1 + nr_edges / 5).\n",LINE,FLAT);
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t--n_globl %snr_moves%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%snr_moves%s is the number of candidate moves on which to perform local\n",LINE,FLAT);  
    PhyML_Printf("\t\tbranch length optimization (default = 1 + nr_edges / 10).\n",LINE,FLAT);
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t--max_dist %sdistance%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%sdistance%s is the maximum regraft distance to consider (default = 1 + nr_edges / 10).\n",LINE,FLAT);
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t--n_optim %snr_moves%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%snr_moves%s is the number of candidates moves on which to perform global\n",LINE,FLAT);  
    PhyML_Printf("\t\tbranch length optimization (default = 100).\n",LINE,FLAT);
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t--n_best %snr_pos%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%snr_pos%s is the number of promising regraft positions to consider when\n",LINE,FLAT);  
    PhyML_Printf("\t\tperforming only the best SPR move. (default = 1 + nr_edges / 5).\n",LINE,FLAT);
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t--collapse_boot%s\n",BOLD,FLAT);
    PhyML_Printf("\t\tBranch length on bootstrap trees are not collapsed if too small.\n");
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t--random_boot%s\n",BOLD,FLAT);
    PhyML_Printf("\t\tSequence order in bootstrapped data set is random.\n");
    
    PhyML_Printf("\n");
    
    PhyML_Printf("%s\n\t--no_gap%s\n",BOLD,FLAT);
    PhyML_Printf("\t\tColumns with ambiguous characters are discarded.\n");
    
    PhyML_Printf("\n");
    
    
    PhyML_Printf("%s\n\t--optBrent %scycles%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%scycles%s specifies the number of times Brent algorithm will be executed during an optimization turn.\n",LINE,FLAT);
    PhyML_Printf("\t\t(default = 2). Differently from BFGS, the Brent algorithm does single parameter optimization.\n",LINE,FLAT);
    */
    PhyML_Printf("%s\n\n\t\tNote: Use of options other than those here isn't yet supported in IgPhyML.\n\n",FLAT,BOLD,FLAT);

    PhyML_Printf("\n");
    
        Exit("");
}

/*********************************************************/


