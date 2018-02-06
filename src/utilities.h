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
   
#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <errno.h>
#include <float.h>

#include "../config.h"


#define For(i,n)                     for(i=0; i<n; i++)
#define Fors(i,n,s)                  for(i=0; i<n; i+=s)
#define PointGamma(prob,alpha,beta)  PointChi2(prob,2.0*(alpha))/(2.0*(beta))
#define SHFT2(a,b,c)                 (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d)               (a)=(b);(b)=(c);(c)=(d);
#define MAX(a,b)                     ((a)>(b)?(a):(b))
#define MIN(a,b)                     ((a)<(b)?(a):(b))
#define SIGN(a,b)                    ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d)                (a)=(b);(b)=(c);(c)=(d);

#ifndef isnan
# define isnan(x)						 \
  (sizeof (x) == sizeof (long double) ? isnan_ld (x)		 \
   : sizeof (x) == sizeof (double) ? isnan_d (x)		 \
   : isnan_f (x))
static inline int isnan_f  (float       x) { return x != x; }
static inline int isnan_d  (double      x) { return x != x; }
static inline int isnan_ld (long double x) { return x != x; }
#endif

#ifndef isinf
# define isinf(x)						 \
  (sizeof (x) == sizeof (long double) ? isinf_ld (x)		 \
   : sizeof (x) == sizeof (double) ? isinf_d (x)		 \
   : isinf_f (x))
static inline int isinf_f  (float       x) { return isnan (x - x); }
static inline int isinf_d  (double      x) { return isnan (x - x); }
static inline int isinf_ld (long double x) { return isnan (x - x); }
#endif
     
#define N_MAX_NEX_COM   20
#define T_MAX_NEX_COM   100
#define N_MAX_NEX_PARM  50
#define T_MAX_TOKEN     200

#define NEXUS_COM    0
#define NEXUS_PARM   1
#define NEXUS_EQUAL  2
#define NEXUS_VALUE  3
#define NEXUS_SPACE  4

#define  NNI_MOVE            0
#define  SPR_MOVE            1
#define  BEST_OF_NNI_AND_SPR 2

#define  M_1_SQRT_2_PI   0.398942280401432677939946059934
#define  M_SQRT_32       5.656854249492380195206754896838 
#define  PI              3.14159265358979311600
#define  SQRT2PI         2.50662827463100024161
#define  LOG2PI          1.83787706640934533908
#define  LOG2            0.69314718055994528623
#define  LOG_SQRT_2_PI   0.918938533204673

#define  NORMAL 1
#define  EXACT  2

#define  PHYLIP 0
#define  NEXUS  1
#define  FASTA  2


#define  YES 1
#define  NO  0

#define  NT 0 /* nucleotides */
#define  AA 1 /* amino acids */
#define  GENERIC 2 /* custom alphabet */

#define  ACGT 0 /* A,G,G,T encoding */
#define  RY   1 /* R,Y     encoding */

#define INTERFACE_DATA_TYPE      0
#define INTERFACE_MULTIGENE      1
#define INTERFACE_MODEL          2
#define INTERFACE_TOPO_SEARCH    3
#define INTERFACE_BRANCH_SUPPORT 4

#ifndef INFINITY
#define INFINITY HUGE
#endif

#define  N_MAX_OPTIONS        100
 

#define  T_MAX_FILE           500
#define  T_MAX_LINE       2000000
#define  T_MAX_NAME          1000
#define  T_MAX_SEQ        2000000
#define  T_MAX_OPTION         100
#define  T_MAX_LABEL          100
#define  T_MAX_STATE            5
#define  N_MAX_LABEL          100
#define  BLOCK_LABELS         100

#define  NODE_DEG_MAX          10
#define  BRENT_ITMAX        10000
#define  BRENT_CGOLD    0.3819660
#define  BRENT_ZEPS        1.e-10
#define  MNBRAK_GOLD     1.618034
#define  MNBRAK_GLIMIT      100.0
#define  MNBRAK_TINY       1.e-20
#define  ALPHA_MIN           0.04
#define  ALPHA_MAX            100
#define  BL_START          1.e-04
#define  GOLDEN_R      0.61803399
#define  GOLDEN_C  (1.0-GOLDEN_R)
#define  N_MAX_INSERT          20
#define  N_MAX_OTU           4000
#define  UNLIKELY          -1.e10
#define  NJ_SEUIL             0.1
#define  ROUND_MAX            100
#define  DIST_MAX            2.00
#define  AROUND_LK           50.0
#define  PROP_STEP            1.0
#define  T_MAX_ALPHABET       100
#define  MDBL_MIN         FLT_MIN
#define  MDBL_MAX         FLT_MAX
#define  POWELL_ITMAX         200
#define  LINMIN_TOL       2.0E-04
#define  SCALE_POW             10    /* Scaling factor will be 2^SCALE_POW or 2^(-SCALE_POW) [[ WARNING: SCALE_POW < 31 ]]*/ // old 10
#define  DEFAULT_SIZE_SPR_LIST 20
#define  OUTPUT_TREE_FORMAT  0 /* 0-->Newick; 1-->Nexus */
#define  MAX_PARS        1000000000

#define  LIM_SCALE_VAL     1.E-50 /* Scaling limit (deprecated) */

#define  MIN_CLOCK_RATE   1.E-10
#define  MIN_VAR_BL        1.E-7

#define MODELPAREPS   1.E-4

#define TREEBOUNDLOW    1.E-4
#define TREEBOUNDHIGH   999.0

#define JC69       1
#define K80        2
#define F81        3
#define HKY85      4
#define F84        5
#define TN93       6
#define GTR        7
#define CUSTOM     8

#define WAG       11
#define DAYHOFF   12
#define JTT       13
#define BLOSUM62  14
#define MTREV     15
#define RTREV     16
#define CPREV     17
#define DCMUT     18
#define VT        19
#define MTMAM     20
#define MTART     21
#define HIVW      22
#define HIVB      23
#define CUSTOMAA  24
#define LG        25

#define COMPOUND_COR   0
#define COMPOUND_NOCOR 1
#define EXPONENTIAL    2
#define GAMMA          3
#define THORNE         4

#define HDISCRETE	0
#define HSEMIEMP	1

/* #define USE_OLD_LK */

/* Uncomment the lines below to switch to single precision */
/* typedef	float phydbl; */
/* #define LOG logf */
/* #define POW powf */
/* #define EXP expf */
/* #define FABS fabsf */
/* #define SQRT sqrtf */
/* #define CEIL ceilf */
/* #define FLOOR floorf */
/* #define RINT rintf */
/* #define ROUND roundf */
/* #define TRUNC truncf */
/* #define COS cosf */
/* #define SIN sinf */
/* #define TAN tanf */
/* #define SMALL FLT_MIN */
/* #define BIG  FLT_MAX */
/* #define SMALL_PIJ 1.E-10 */
/* #define BL_MIN 1.E-5 */
/* #define  P_LK_LIM_INF   2.168404e-19 /\* 2^-62 *\/ */
/* #define  P_LK_LIM_SUP   4.611686e+18 /\* 2^62 *\/ */

/* Uncomment the line below to switch to double precision */
//typedef	double phydbl;
typedef  double phydbl;
#define LOG log
#define POW pow
#define EXP exp
#define FABS fabs
#define SQRT sqrt
#define CEIL ceil
#define FLOOR floor
#define RINT rint
#define ROUND round
#define TRUNC trunc
#define COS cos
#define SIN sin
#define TAN tan
#define SMALL DBL_MIN
#define BIG  DBL_MAX
#define SMALL_PIJ 1.E-10 //!Adjustable. original value 1e-10 //!< Added by Marcelo; //! Likelihood calculation needs to be revised ... 
#ifndef MC
#define BL_MIN 1.E-8 //!Adjustable. original value 1e-8 //! Likelihood calculation needs to be revised ... 
#define BL_MAX 100.
#else
#define BL_MIN 1.E-5
#define BL_MAX 5.
#endif
#define  P_LK_LIM_INF   3.054936e-150 /* 2^-500 */
#define  P_LK_LIM_SUP   3.273391e+150 /* 2^500 */

/*********************************************************/
typedef struct _Ts_and_Tv{ //!< Added by Marcelo;
int nts;
int ntv;
int ndiff;
int sSub; /*!< 1 if synonymous substitution, 0 otherwise.*/
int caseK07; /*!< 0-8 according to Kosiol 2007.*/
}ts_and_tv;
/*********************************************************/
typedef struct lkexperiment{ //!< Added by Marcelo;
phydbl *omega_part; //modified by ken 18/8
int		nomega_part;
phydbl kappa;
phydbl lk;
}lkexperiment;
/*********************************************************/

typedef struct __Node {
  struct __Node                       **v; /* table of pointers to neighbor nodes. Dimension = 3 */
  struct __Node               ***bip_node; /* three lists of pointer to tip nodes. One list for each direction */
  struct __Edge                       **b; /* table of pointers to neighbor branches */
  struct __Node                      *anc; /* direct ancestor t_node (for rooted tree only) */
  struct __Node                 *ext_node;

  int								lupdate;	//Added by Ken - update likelihood on ancestral edge?
  struct __Edge                     *anc_edge;  //Added by Ken - pointer to ancestral edge


  int                           *bip_size; /* Size of each of the three lists from bip_node */
  int                                 num; /* t_node number */
  int                                 tax; /* tax = 1 -> external node, else -> internal t_node */
  int                        check_branch; /* check_branch=1 is the corresponding branch is labelled with '*' */
  char                              *name; /* taxon name (if exists) */
  char                          *ori_name; /* taxon name (if exists) */

  phydbl                           *score; /* score used in BioNJ to determine the best pair of nodes to agglomerate */
  phydbl                               *l; /* lengths of the (three or one) branche(s) connected this t_node */
  phydbl                     dist_to_root; /* distance to the root t_node */

  short int                        common;
  phydbl                           y_rank;
  phydbl                       y_rank_ori;
  phydbl                       y_rank_min;
  phydbl                       y_rank_max;

}t_node;


/*********************************************************/

typedef struct __Edge {
  /*
    syntax :  (node) [edge]
(left_1) .                   .(right_1)
          \ (left)  (right) /
           \._____________./
           /    [b_fcus]   \
          /                 \
(left_2) .                   .(right_2)

  */

  struct __Node               *left,*rght; /* t_node on the left/right side of the t_edge */
  short int         l_r,r_l,l_v1,l_v2,r_v1,r_v2;
  /* these are directions (i.e., 0, 1 or 2): */
  /* l_r (left to right) -> left[b_fcus->l_r] = right */
  /* r_l (right to left) -> right[b_fcus->r_l] = left */
  /* l_v1 (left t_node to first t_node != from right) -> left[b_fcus->l_v1] = left_1 */
  /* l_v2 (left t_node to secnd t_node != from right) -> left[b_fcus->l_v2] = left_2 */
  /* r_v1 (right t_node to first t_node != from left) -> right[b_fcus->r_v1] = right_1 */
  /* r_v2 (right t_node to secnd t_node != from left) -> right[b_fcus->r_v2] = right_2 */

  struct __NNI                       *nni;

  struct __Edge						*anc_edge;
  struct __Node						*anc_node;
  struct __Node						*des_node;

  int								up_upp;
  int                                 num; /* branch number */
  phydbl                                l; /* branch length */
  phydbl                           best_l; /* best branch length found so far */
  phydbl                            l_old; /* old branch length */

  int                           bip_score; /* score of the bipartition generated by the corresponding edge
					      bip_score = 1 iif the branch is found in both trees to be compared,
					      bip_score = 0 otherwise. */

  phydbl							  **upp; //upp matrix from nhphyml paper

  phydbl            *p_lk_left,*p_lk_rght; /* likelihoods of the subtree on the left and
					      right side (for each site and each relative rate category) */
  short int      *p_lk_tip_r, *p_lk_tip_l; 
  short int           *div_post_pred_left; /* posterior prediction of nucleotide/aa diversity (left-hand subtree) */
  short int           *div_post_pred_rght; /* posterior prediction of nucleotide/aa diversity (rght-hand subtree) */


  phydbl                          *Pij_rr; /* matrix of change probabilities and its first and secnd derivates */

  phydbl						 **bPmat_part; //separate p matrixs Added by Ken 17/8

  int                     *pars_l,*pars_r; /* parsimony of the subtree on the left and right sides (for each site) */
  unsigned int               *ui_l, *ui_r; /* union - intersection vectors used in Fitch's parsimony algorithm */
  int                *p_pars_l, *p_pars_r; /* conditional parsimony vectors */

  int                         num_st_left; /* number of the subtree on the left side */
  int                         num_st_rght; /* number of the subtree on the right side */


  /* Below are the likelihood scaling factors (used in functions
     `Get_All_Partial_Lk_Scale' in lk.c */
  int                 *sum_scale_left_cat;
  int                 *sum_scale_rght_cat;
  int                     *sum_scale_left;
  int                     *sum_scale_rght;
  int                 *sum_scale_upp_cat;
  int                     *sum_scale_upp;

  phydbl                          bootval; /* bootstrap value (if exists) */

  short int                      is_alive; /* is_alive = 1 if this t_edge is used in a tree */

  phydbl                   dist_btw_edges;
  int                 topo_dist_btw_edges;

  int                     has_zero_br_len;

  phydbl                       ratio_test; /* approximate likelihood ratio test */
  phydbl                   alrt_statistic; /* aLRT statistic */

  char                           **labels; /* string of characters that labels the corresponding t_edge */
  int                            n_labels; /* number of labels */
  int                             n_jumps; /* number of jumps of substitution rates */
}t_edge;

/*********************************************************/

typedef struct __Arbre {

  struct __Node                       *n_root; /*!< root t_node.*/
  struct __Edge                       *e_root; /*!< t_edge on which lies the root.*/
  struct __Node                       **noeud; /*!< array of nodes that defines the tree topology.*/
  struct __Edge                     **t_edges; /*!< array of edges.*/
  struct __Model                         *mod; /*!< substitution model. */
  struct __Calign                       *data; /*!< sequences. */
  struct __Option                         *io; /*!< input/output.*/
  struct __Matrix                        *mat; /*!< pairwise distance matrix. */
  struct __Node                   **curr_path; /*!< list of nodes that form a path in the tree. */
  struct __SPR                     **spr_list; /*!< list of the possible spr moves.*/
  struct __SPR                      *best_spr; /*!< best spr move.*/
  struct __Tdraw                     *ps_tree; /*!< structure for drawing trees in postscript format.*/
  struct __Trate                       *rates; /*!< structure for handling rates of evolution. */
  struct __Tmcmc                        *mcmc; /*!<what?. */
  struct __Triplet            *triplet_struct; /*!<what?. */

  int                          ps_page_number; /*!< when multiple trees are printed, this variable give the current page number. */
  int                         depth_curr_path; /*!< depth of the t_node path defined by curr_path. */
  int                                 has_bip; /*!<if has_bip=1, then the structure to compare tree topologies is allocated, has_bip=0 otherwise.*/
  int                                   n_otu; /*!< number of taxa. */
  int                               curr_site; /*!< current site of the alignment to be processed.*/
  int                               curr_catg; /*!< current class of the discrete gamma rate distribution. */
  int                                  n_swap; /*!< number of NNIs performed. */
  int                               n_pattern; /*!< number of distinct site patterns.*/
  int                      has_branch_lengths; /*!< =1 iff input tree displays branch lengths.*/
  int                          print_boot_val; /*!< if print_boot_val=1, the bootstrap values are printed. */
  int                          print_alrt_val; /*!< if print_boot_val=1, the bootstrap values are printed.*/
  int                              both_sides; /*!< both_sides=1 -> a pre-order and a post-order tree traversals are required to compute the likelihood of every subtree in the phylogeny.*/
  int               num_curr_branch_available; /*!<gives the number of the next cell in t_edges that is free to receive a pointer to a branch.*/
  short int                            *t_dir; /*!<what?. */
  int                          n_improvements; /*!<what?. */
  int                                 n_moves; /*!<what?. */

  int                                      dp; /*!< Data partition. */
  int                               s_mod_num; /*!< Substitution model number.*/
  int                               lock_topo; /*!< = 1 any subsequent topological modification will be banished.*/
  int                            print_labels; /*!<what?. */

  phydbl                              init_lnL;/*!<what?. */ 
  phydbl                              best_lnL; /*!< highest value of the loglikelihood found so far. */
  int                                best_pars; /*!< highest value of the parsimony found so far. */
  phydbl                                 c_lnL; /*!< loglikelihood. */
  phydbl                               old_lnL; /*!< old loglikelihood. */
  phydbl                     sum_min_sum_scale; /*!< common factor of scaling factors.*/
  phydbl                         *c_lnL_sorted; /*!< used to compute c_lnL by adding sorted terms to minimize CPU errors.*/
  phydbl                          *cur_site_lk; /*!< vector of likelihoods at individual sites.*/
  phydbl                          *old_site_lk; /*!< vector of likelihoods at individual sites.*/
  phydbl                     **log_site_lk_cat; /*!< loglikelihood at individual sites and for each class of rate.*/
  phydbl                          *site_lk_cat; /*!< loglikelihood at a single site and for each class of rate.*/
  phydbl                       unconstraint_lk; /*!< unconstrained (or multinomial) likelihood  .*/
  phydbl                        **log_lks_aLRT; /*!< used to compute several branch supports .*/
  phydbl                            n_root_pos; /*!< position of the root on its t_edge .*/
  phydbl                                  size; /*!< tree size .*/
  int                              *site_pars; /*!<what?. */
  int                                  c_pars; /*!<what?. */
  int                               *step_mat; /*!<what?. */

  int                           size_spr_list; /*!<what?. */
  int                  perform_spr_right_away; /*!<what?. */

  time_t                                t_beg; /*!<what?. */
  time_t                            t_current; /*!<what?. */
  
  int                     bl_from_node_stamps; /*!< == 1 -> Branch lengths are determined by t_node times.*/
  phydbl                        sum_y_dist_sq; /*!<what?. */
  phydbl                           sum_y_dist; /*!<what?. */
  phydbl                      tip_order_score; /*!<what?. */
  phydbl                      logLk_correction; /*!< correction factor loglikelihood according to Seo Kishino 2009. */
  phydbl                  logLk_correction_gaps; /*!< correction factor (with gaps) loglikelihood according to Seo Kishino 2009. */
  phydbl                      logLk_correction_approx; /*!< correction factor loglikelihood according to Seo Kishino 2009. */
  phydbl                  logLk_correction_gaps_approx; /*!< correction factor (with gaps) loglikelihood according to Seo Kishino 2009. */
  int                       br_len_invar_set; /* is Br_Len_Involving_Invar already done? */
}t_tree;

/*********************************************************/

typedef struct __Super_Arbre {
  struct __Arbre                           *tree;
  struct __List_Arbre                  *treelist; /* list of trees. One tree for each data set to be processed */
  struct __Calign                    *curr_cdata;
  struct __Option                   **optionlist; /* list of pointers to input structures (used in supertrees) */

  struct __Node           ***match_st_node_in_gt;
  /*  match_st_in_gt_node[subdataset number][supertree t_node number]
   *  gives the t_node in tree estimated from 'subdataset number' that corresponds
   *  to 'supertree t_node number' in the supertree
   */

  struct __Node           *****map_st_node_in_gt;
  /*  mat_st_gt_node[gt_num][st_node_num][direction] gives the
   *  t_node in gt gt_num that maps t_node st_node_num in st.
   */

  struct __Edge             ***map_st_edge_in_gt;
  /*  map_st_gt_br[gt_num][st_branch_num] gives the
   *  branch in gt gt_num that maps branch st_branch_num
   *  in st.
   */

  struct __Edge            ****map_gt_edge_in_st;
  /*  mat_gt_st_br[gt_num][gt_branch_num][] is the list of
   *  branches in st that map branch gt_branch_num
   *  in gt gt_num.
   */

  int                   **size_map_gt_edge_in_st;
  /*  size_map_gt_st_br[gt_num][gt_branch_num] gives the
   *  size of the list map_gt_st_br[gt_num][gt_branch_num][]
   */


  struct __Edge             ***match_st_edge_in_gt;
  /* match_st_edge_in_gt[gt_num][st_branch_num] gives the
   * branch in gt gt_num that matches branch st_branch_num
   */

  struct __Edge             ***match_gt_edge_in_st;
  /* match_gt_edge_in_st[gt_num][gt_branch_num] gives the
   * branch in st that matches branch gt_branch_num
   */

  struct __Node                  ****closest_match;
  /* closest_match[gt_num][st_node_num][dir] gives the 
   * closest t_node in st that matches a t_node in gt gt_num
   */

  int                              ***closest_dist;
  /* closest_dist[gt_num][st_node_num][dir] gives the
   * number of edges to traverse to get to node
   * closest_match[gt_num][st_node_num][dir]
   */

  int                                         n_part;
  /* number of trees */

  phydbl                                      **bl;
  /* bl[gt_num][gt_branch_num] gives the length of
   * branch gt_branch_num
   */

  phydbl                                  **bl_cpy;
  /* copy of bl */

  phydbl                                     **bl0;
  /* bl estimated during NNI (original topo) 
   * See Mg_NNI.
   */

  phydbl                                     **bl1;
  /* bl estimated during NNI (topo conf 1) 
   * See Mg_NNI.
   */

  phydbl                                     **bl2;
  /* bl estimated during NNI (topo conf 2) 
   * See Mg_NNI.
   */

  int                                *bl_partition;
  /* partition[gt_num] gives the t_edge partition number 
   * gt_num belongs to.
   */
  int                                   n_bl_part;

  struct __Model                          **s_mod; /* substitution model */

  int                                     n_s_mod;
  int                                 lock_br_len;

}supert_tree;

/*********************************************************/

typedef struct __List_Arbre { /* a list of trees */
  struct __Arbre   **tree;
  int           list_size;                /* number of trees in the list */
}t_treelist;

/*********************************************************/

typedef struct __Align {
  char           *name; /*!< sequence name.  */
  int              len; /*!< sequence length.*/
  char          *state; /*!<sequence itself. */
  short int *is_ambigu; /*!<is_ambigu[site] = 1 if state[site] is an ambiguous character. 0 otherwise. */
  
  char **alternativeCodons;                                                                                         //!< Added by Marcelo.
  int                ntLen;                                                                                         //!< Added by Marcelo. 
  char           *ntStates;                                                                                         //!< Added by Marcelo. 
}align;

/*********************************************************/


typedef struct __Calign{
  struct __Align **c_seq;             /*!< compressed sequences.*/
  phydbl          *b_frq;             /*!< observed state frequencies.*/
  short int       *invar;             /*!< 1 -> states are identical, 0 states vary. */
  int              *wght;             /*!< # of each site in c_align.*/
  short int      *ambigu;             /*!< ambigu[i]=1 if one or more of the sequences at site i display an ambiguous character.*/
  phydbl      obs_pinvar;
  int              n_otu;             /*!< number of taxa.*/
  int          clean_len;             /*!< uncrunched sequences lenghts without gaps.*/
  int         crunch_len;             /*!< crunched sequences lengths.*/
  int           init_len;             /*!< length of the uncompressed sequences.*/
  int          *sitepatt;             /*!< this array maps the position of the patterns in the compressed alignment to the positions in the uncompressed one. */
  int             format;             /*!< 0 (default): PHYLIP. 1: NEXUS. */
  
  
}calign;

/*********************************************************/

typedef struct __Matrix { /* mostly used in BIONJ */
  phydbl    **P,**Q,**dist; /* observed proportions of transition, transverion and  distances
			       between pairs of  sequences */

  t_tree             *tree; /* tree... */
  int              *on_off; /* on_off[i]=1 if column/line i corresponds to a t_node that has not
			       been agglomerated yet */
  int                n_otu; /* number of taxa */
  char              **name; /* sequence names */
  int                    r; /* number of nodes that have not been agglomerated yet */
  struct __Node **tip_node; /* array of pointer to the leaves of the tree */
  int             curr_int; /* used in the NJ/BIONJ algorithms */
  int               method; /* if method=1->NJ method is used, BIONJ otherwise */
}matrix;

/*********************************************************/
/*! \brief Model applied in the phylogeny inference.
 *
 *  Prexisting: Nucleotide and Amino Acid Substitution, Adding Codon models.
*/
typedef struct __Model {
  struct __Optimiz  *s_opt; /*!< pointer to parameters to optimize. */
  struct __Eigen    *eigen;
  struct __M4       *m4mod;
  struct __Option      *io;
  
  int                   ns; /*!< number of states (4 for DNA, 20 for AA, #senseCodons for Codons).*/                             
  char          *modelname;
  char  *custom_mod_string; /*!< string of characters used to define custom models of substitution */
  int              *rr_num; 
  int        *n_rr_per_cat; /*!< number of rate parameters in each category. */
  int            n_diff_rr; /*!< number of different relative substitution rates in the custom model. */
  int         update_eigen; /*!< update_eigen=1-> eigen values/vectors need to be updated. */
  int           whichmodel;
  int               n_catg; /*!< number of categories in the discrete gamma distribution. */
  int                invar; /*!< =1 iff the substitution model takes into account invariable sites. */
  int            bootstrap; /*!< Number of bootstrap replicates (0 : no bootstrap analysis is launched). */
  int            use_m4mod; /*!< Use a Makrkov modulated Markov model ?. */
  int         gamma_median; /*!< 1: use the median of each bin in the discrete gamma distribution. 0: the mean is used. */
  phydbl               *pi; /*!< states frequencies. */
  phydbl      *pi_unscaled; /*!< states frequencies (unscaled). */
  phydbl    *gamma_r_proba; /*!< probabilities of the substitution rates defined by the discrete gamma distribution. */
  phydbl         *gamma_rr; /*!< substitution rates defined by the discrete gamma distribution. */
  phydbl             kappa; /*!< transition/transversion rate. */
  phydbl            lambda; /*!< parameter used to define the ts/tv ratios in the F84 and TN93 models. */
  phydbl             alpha; /*!< gamma shapa parameter.*/
  phydbl            pinvar; /*!< proportion of invariable sites.*/
  phydbl         alpha_old;
  phydbl         kappa_old;
  phydbl        lambda_old;
  phydbl        pinvar_old;
  phydbl               *rr; /*!< relative rate parameters of the GTR or custom model (given by rr_val[rr_num[i]]). */
  phydbl           *rr_val; /*!< relative rate parameters of the GTR or custom model. */
  phydbl           *Pij_rr; /*!< matrix of change probabilities. */
  phydbl                mr; /*!< mean rate = branch length/time interval  mr = -sum(i)(vct_pi[i].mat_Q[ii]).*/
  phydbl      *user_b_freq; /*!< user-defined nucleotide frequencies.*/
  phydbl             *qmat; /*! Q matrix of instantaneous rates.*/
  phydbl        **qmat_buff_part;
  phydbl        *rr_branch; /*!< relative substitution rates on each branch, for the whole set of sites.*/
  phydbl      *p_rr_branch; /*!< corresponding frequencies.*/
  int          n_rr_branch; /*!< number of classes.*/
  phydbl   rr_branch_alpha; /*!< Shape of the gamma distribution that defines the rr_branch and p_rr_branch values. */
  int            state_len; /*!< lenght of the strings that compose the alphabet (e.g. 1 for DNA, CODON and Amino acid).*/
  
  int         genetic_code; /*!< 0: standard. The genetic code for codon models.*///!< Added by Marcelo.
  int           freq_model; /*!Model equilibrium frequencies in codon models.*/   //!< Added by Marcelo.
  phydbl         	*omega_part; /*!< parameter used to define the dn/ds ratio in GY (simplified) model.*///!< Added by Marcelo. //Changed by Ken 18/8/2016
  int				nomega_part;
  phydbl         omega_old;                                     //!< Added by Marcelo.
  char           *nt_or_cd; /*!< nucleotide or codon data ?. */ //!< Added by Marcelo.
  phydbl        *base_freq;                                     //!< Added by Marcelo.
  int        num_base_freq;                                     //!< Added by Marcelo.
  phydbl    *uns_base_freq;                                     //!< Added by Marcelo.
  phydbl           *pkappa;                                     //!< Added by Marcelo.
  phydbl        *unspkappa;                                     //!< Added by Marcelo.
  int               nkappa;                                     //!< Added by Marcelo.
  int             n_w_catg; /*!< number of omega rate ratio categories.*///!< Added by Marcelo.
  struct __Eigen **MW_eigen; /*!< Added by Marcelo.*/
  phydbl              beta; /*!< Added by Marcelo.*/
  phydbl          beta_old; /*!< Added by Marcelo.*/
  phydbl  *prob_omegas_uns; /*!< Added by Marcelo.*/
  int         omegaSiteVar; /*!< Added by Marcelo.*/

  int                **ipiv_part; /*!< Added by Marcelo.*/
  phydbl           **matAux_part; /*!< Added by Marcelo.*/// Modified by Ken 22/8
  phydbl          **Apowers_part; /*!< Added by Marcelo.*/
  phydbl                **U_part; /*!< Added by Marcelo.*/
  phydbl                **V_part; /*!< Added by Marcelo.*/
  phydbl               **A2_part; /*!< Added by Marcelo.*/
  phydbl               **A4_part; /*!< Added by Marcelo.*/
  phydbl               **A6_part; /*!< Added by Marcelo.*/
  phydbl               **A0_part; /*!< Added by Marcelo.*/
  phydbl               **A8_part; /*!< Added by Marcelo.*/

  phydbl               *A3; /*!< Added by Marcelo.*/
  phydbl               *A5; /*!< Added by Marcelo.*/ 
  phydbl       *qmatScaled; /*!< Added by Marcelo.*/
  phydbl             *mr_w; /*!< Added by Marcelo.*/
  phydbl           *omegas; /*!< Added by Marcelo.*/
  phydbl      *prob_omegas; /*!< Added by Marcelo.*/
  phydbl         pcsC[100]; /*!< Added by Marcelo.*/
  int   codon_model_nature; /*!< Added by Marcelo.*/
  int  calculate_init_tree; /*!< Added by Marcelo.*/
  int                 npcs; /*!< Added by Marcelo.*/
  int             pcaModel;
    
  int            initqrates; /*!< Added by Stefan.*/
  int           whichrealmodel;
  
  phydbl  userRates[64][64];
  phydbl       userfreq[64];
  phydbl      userbfreq[64];
  phydbl userRatesT[64][64];
  phydbl      userfreqT[64];
  phydbl     userbfreqT[64];

  //Added by Kenneth Hoehn 3/6/2016
 phydbl		*Bmat;
 int		nmotifs; //number of motifs
 char		**motifs; //motifs to analyze
 int		*motif_hotness; //maps index in motifs array to index in hotness array
 int		nhotness; //number of h parameters to optimize
 double		*hotness;
 int		*hoptindex;


 //partition model stuff
 int		*partIndex; //which q matrix to use for each partition
 phydbl     **qmat_part; /*! array of each Q matrix.*/
 phydbl		**Pmat_part; /*!< array of p matrixes */
 int		nparts; //number of partitions
 char**		partNames; //names of partitions
 char*		partfile; //name of file containing partition model
 int		partfilespec; //=1 if partfile specified


 phydbl		**hotspotcmps; //hotspot comparison tables
 phydbl		*aacmps; //AA exchangeability matrix comparison
 phydbl		aaintercept;
 phydbl 	aaslope;
 char*		aamodel;
 int		opthotness;
 char		*rootname;
 char		*motifstring;
 char		*hotnessstring;
 int		startnode;
 int 		slowSPR;

 phydbl		**rootStates; //61xlength matrix of 1s and 0s describing the character state at the root

 //added by Ken 19/9
 phydbl		hintercept;
 phydbl		hslope;
 int		hmode;
 int		omode;

 int		motifstringopt;
 int		hotnessstringopt;
 int		rootfound;

 int		ambigprint;
 char		*ambigfile;

 phydbl 	stretch; //stretch branch lengths of initial tree

}model;

/*********************************************************/

typedef struct __Eigen{
  int              size;
  phydbl             *q; /* matrix which eigen values and vectors are computed */
  phydbl         *space;
  int        *space_int;
  phydbl         *e_val; /* eigen values (vector), real part. */
  phydbl      *e_val_im; /* eigen values (vector), imaginary part */
  phydbl      *r_e_vect; /* right eigen vector (matrix), real part */
  phydbl   *r_e_vect_im; /* right eigen vector (matrix), imaginary part */
  phydbl      *l_e_vect; /* left eigen vector (matrix), real part */
}eigen;

/*********************************************************/
/*! \brief Run options
 *
 *  mostly used in 'help.c', 'cl.c' and 'interface.c' */
typedef struct __Option { 
  struct __Model                *mod;/*!< pointer to a substitution model. */
  struct __Arbre               *tree;/*!< pointer to the current tree. */
  struct __Align              **data;/*!< pointer to the uncompressed sequences. */
  struct __Calign             *cdata;/*!< pointer to the compressed sequences. */
  struct __Super_Arbre           *st;/*!< pointer to supertree. */
  struct __Tnexcom    **nex_com_list;
  struct __List_Arbre      *treelist;/*!< list of trees. */


  int                    interleaved; /*!< interleaved or sequential sequence file format ?. */
  int                        in_tree; /*!< =1 iff a user input tree is used as input. */

  char                *in_align_file; /*!< alignment file name. */
  FILE                  *fp_in_align; /*!< pointer to the alignment file. */

  char                 *in_tree_file; /*!< input tree file name. */
  FILE                   *fp_in_tree; /*!< pointer to the input tree file. */

  char                *out_tree_file; /*!< name of the tree file. */
  FILE                  *fp_out_tree;

  char               *out_trees_file; /*!< name of the random tree file. */
  FILE                 *fp_out_trees; /*!< pointer to the tree file containing all the trees estimated using random starting trees. */

  char           *out_boot_tree_file; /*!< name of the tree file. */
  FILE             *fp_out_boot_tree; /*!< pointer to the bootstrap tree file. */

  char          *out_boot_stats_file; /*!< name of the tree file. */
  FILE            *fp_out_boot_stats; /*!< pointer to the statistics file. */

  char               *out_stats_file; /*!< name of the statistics file.*/
  FILE                 *fp_out_stats;
  int               out_stats_format;

  //char               *out_trace_file; /*!< name of the file in which the likelihood of the model is written. */
  //FILE                 *fp_out_trace;

  char               *out_trace_tree_file; /*!< name of the file in which the likelihood of the model is written. */
  FILE                 *fp_out_tree_trace;

  char               *out_trace_stats_file; /*!< name of the file in which the likelihood of the model is written. */
  FILE                 *fp_out_stats_trace;
  int					tracecount;

  char                  *out_lk_file; /*!< name of the file in which the likelihood of the model is written. */
  FILE                    *fp_out_lk;

  char                  *out_ps_file; /*!< name of the file in which tree(s) is(are) written. */
  FILE                    *fp_out_ps;

  char              *clade_list_file;
  
  int                       datatype; /*!< 0->DNA, 1->AA, 3->CODON.*/
  int               print_boot_trees; /*!< =1 if the bootstrapped trees are printed in output. */
  int       out_stats_file_open_mode; /*!< opening file mode for statistics file. */
  int        out_tree_file_open_mode; /*!< opening file mode for tree file. */
  int                    n_data_sets; /*!< number of data sets to be analysed. */
  int                        n_trees; /*!< number of trees. */
  int                       init_len; /*!< sequence length. */
  int                          n_otu; /*!< number of taxa. */
  int               n_data_set_asked; /*!< number of bootstrap replicates. */
  char                     *nt_or_cd; /*!< nucleotide or codon data ?. */ 
  int                      multigene; /*!< if=1 -> analyse several partitions. */
  int               config_multigene;
  int                         n_part; /*!< number of data partitions. */
  int                        curr_gt;
  int                     ratio_test; /*!< from 1 to 4 for specific branch supports, 0 of not. */
  int                    ready_to_go;
  int                data_file_format; /*!< Data format: Phylip or Nexus. */
  int                tree_file_format; /*!< Tree format: Phylip or Nexus. */
  int                 curr_interface;
  int                         r_seed; /*!< random seed. */
  int                  collapse_boot; /*!< 0 -> branch length on bootstrap trees are not collapsed if too small. */
  int          random_boot_seq_order; /*!< 0 -> sequence order in bootstrapped data set is random. */
  int                    print_trace;
  int                 print_site_lnl;
  int                       m4_model;
  int                      rm_ambigu; /*!< 0 is the default. 1: columns with ambiguous characters are discarded prior further analysis. */
  int                       colalias;
  int                  append_run_ID;
  char                *run_id_string;
  int                          quiet; /*!< 0 is the default. 1: no interactive question (for batch mode). */
  char                    **alphabet;
  char              **long_tax_names;
  char             **short_tax_names;
  int                 size_tax_names;
  phydbl                   *z_scores;
  phydbl                        *lat;
  phydbl                        *lon;
  
  int             init_DistanceTreeCD; /*!< 0: ML JC69; 1: ML M0; 2: Schneider 2005 CodonPam; 3: Kosiol 2007 Empirical.*/ //!<Added by Marcelo.
  int                        kappaECM; /*!< According to Kosiol 2007.*///!<Added by Marcelo.
  ts_and_tv          *structTs_and_Tv; /*!< Counts for transitions and transversions in codon models.*/ //!<Added by Marcelo.
  FILE                *fp_out_compare; //!< Added by Marcelo. 
  int                            expm; /*!< 0 Eigenvalue; 1 scaling and squaring Pade. approx. */ //!< Added by Marcelo.
  int                    lkExperiment; /*!< plot likelihood for 2 parameters omega and kappa.*/
  phydbl                lkExpStepSize; /*!< Step Size.*/
  phydbl                     minParam; /*!< min Param.*/
  phydbl                     maxParam; /*!< max Param.*/
  int                   heuristicExpm; /*!< TAYLOR in parameter opt?.*/
  int                   threshold_exp; /*!< treshold experiment.*/
  int                         dataset; /*!< treshold experiment.*/
  int                    testInitTree; /*!< Test Initial Tree.*/
  int                   testcondition; /*! Test Condition of the matrix.*/
  int                convert_NT_to_AA; /*! 0 not, 1 yes.*/ 
  phydbl        convert_AAfreq2CDfreq; /*!< 0 approx; 1 solve AX=B.*/ 
  FILE                  *fp_in_usrECM; //!< Added by Marcelo. 
  int            opt_heuristic_manuel; //!< Added by Marcelo. 
  int                   n_termsTaylor; //! Added by Marcelo.
  int                        init_run; //! Added by Marcelo.
  int                        optParam; //! Added by Marcelo.
  int                  roundMax_start; //! Added by Marcelo.
  int                    roundMax_end; //! Added by Marcelo.
  int                eq_freq_handling; //! Added by Marcelo.
  int                            nobl; //! Added by Marcelo.
  int                         genCode;
  int                        omegaOpt;
  int                        wclasses;
  int                  userWantsKappa;
  char *                    userOmega;
  int                    modeltypeOpt;
  int                    freqmodelOpt;
  int                    freqHandling;
  int                   initqratesOpt;
  char *                    userFreqs;
  int                            npcs;
  int                       writemode;
  int                    open_ps_file;
  int                      confformat; // what format does the configuration file have?
  int                         logtree; // should current trees be logged to disc? 0, 1 -> single tree, 2 -> collect all trees
  int                     treecounter; // if logtree = 2, holds a counter for the filename
 }option;

/*********************************************************/
//! Parameters to be optimized by maximum likelihood algorithm.

/*! mostly used in 'optimiz.c' 
*/

typedef struct __Optimiz { 
  int                 print; /*!< =1 -> verbose mode.*/
  int             opt_alpha; /*!< =1 -> the gamma shape parameter is optimised.*/
  int             opt_kappa; /*!< =1 -> the ts/tv ratio parameter is optimised.*/
  int            opt_lambda; /*!< =1 -> the F84|TN93 model specific parameter is optimised.*/
  int            opt_pinvar; /*!< =1 -> the proportion of invariants is optimised.*/
  int        opt_state_freq; /*!< =1 -> the nucleotide frequencies are optimised.*/
  int                opt_rr; /*!< =1 -> the relative rate parameters of the GTR or the custom model are optimised. */
  int       opt_subst_param; /*!< if opt_topo=0 and opt_subst_param=1 -> the numerical parameters of the model are optimised. if opt_topo=0 and opt_subst_param=0 -> no parameter is optimised. */
  int         opt_cov_delta; /*! what is cov_delta?.*/ 
  int         opt_cov_alpha; /*! what is cov_delta?.*/ 
  int    opt_cov_free_rates; /*! what is cov_free_rates?.*/ 
  int                opt_bl; /*!< =1 -> the branch lengths are optimised.*/
  int              opt_topo; /*!< =1 -> the tree topology is optimised. */
  int           topo_search; /*!< generate new topology by 0-NNI, 1-SPR and 2-BEST of NNI and SPR.*/ 
  phydbl            init_lk; /*!< initial loglikelihood value.*/
  int              n_it_max; /*!< maximum number of iteration during an optimisation step.*/
  int              last_opt; /*!< =1 -> the numerical parameters are optimised further while the tree topology remains fixed.*/
  int     random_input_tree; /*!< =true Generate rando input tree.*/
  int         n_rand_starts; /*!< Nunmber of random starting trees.*/
  int          brent_it_max; /*!< Maximum number of iterations for Brent's algorithm.*/ 
  int             steph_spr; /*!< What is steph_spr?.*/
  int       user_state_freq; /*!< =1 Equilibrium frequencies provided by the user.*/
  int       opt_five_branch; /*!< What is opt_five_branch?.*/
  int           pars_thresh; /*!< Parsimony threshold.*/
  int         hybrid_thresh; /*!< What is hybrid threshold?.*/ 
  phydbl     tree_size_mult; /*!< tree size multiplier ... multiplier for what?.*/
  phydbl  min_diff_lk_local; /*!< min_diff_lk_local ... local threshold in which sense ... min improvement until it stops?.*/
  phydbl min_diff_lk_global; /*!< min_diff_lk_global ... global threshold in which sense?.*/
  phydbl   min_diff_lk_move; /*!< min_diff_lk_move ... move threshold in which sense?.*/
  phydbl p_moves_to_examine; /*!< What is p_moves_to_examine?.*/ 
  int              fast_nni; /*!< first improvement is taken?.*/ 
  int                greedy; /*!< best improvement is taken?.*/
  int          general_pars; /*!< What is general_pars?.*/
  int            quickdirty; /*!< What is quickdirty?.*/ 
  int              spr_pars; /*!< What is spr_pars?.*/ 
  int               spr_lnL; /*!< What is spr_lnL?.*/
  int        max_depth_path; /*!< What is max_depth_path?.*/ 
  int        min_depth_path; /*!< What is min_depth_path?.*/ 
  int          deepest_path; /*!< What is deeptest_path?.*/ 
  phydbl  max_delta_lnL_spr; /*!< What is max_delta_lnL_spr?.*/ 
  int           wim_n_rgrft; /*!< Number of promising regraft positions to consider when performing all improving SPR moves.*/
  int           wim_n_globl; /*!< Number of candidate moves on which to perform local branch length optimization.*/
  int          wim_max_dist; /*!< Maximum regraft distance to consider.*/
  int           wim_n_optim; /*!< Number of candidates moves on which to perform global branch length optimization.*/
  int            wim_n_best; /*!< Number of promising regraft positions to consider when performing only the best SPR move.*/
  int        wim_inside_opt; /*!< What is wim_inside_opt?.*/
  
  int                    opt_omega; /*!< =1 -> the dn/ds rate ratio parameter of GY(simplified) codon model is optimised.*///!< Added by Marcelo.
  int                     opt_beta; /*!< 0: no, 1: yes.*///!< Added by Marcelo.
  int               opt_prob_omega; /*!< 0: no, 1: yes.*///!< Added by Marcelo. For discrete model for site variability of omega.
  int                  opt_alphaCD; /*!< 0: no, 1: yes.*///!< Added by Marcelo. For discrete model for site variability of omega.
  phydbl   min_diff_lk_codonModels; /*!< Trial and error threshold.*/
  phydbl min_diff_lk_codonModelsBL; /*!< Trial and error threshold.*/
  int                   opt_method;  /*!< 0- paml, 1-phyml.*/
  int                 nBrentCycles; /*! how many times cycle through parameters using brents?*/
  int          opt_state_freq_AAML; /*! ML AA freqs.*/
}optimiz;
  
/*********************************************************/

typedef struct __NNI{

  struct __Node         *left;
  struct __Node         *rght;
  struct __Edge            *b;

  phydbl                score;
  phydbl               init_l;
  phydbl              init_lk;
  phydbl               best_l;
  phydbl          lk0,lk1,lk2;
  phydbl             l0,l1,l2;

  struct __Node *swap_node_v1;
  struct __Node *swap_node_v2;
  struct __Node *swap_node_v3;
  struct __Node *swap_node_v4;

  int       best_conf;   /* best topological configuration :
			    ((left_1,left_2),right_1,right_2) or
			    ((left_1,right_2),right_1,left_2) or
			    ((left_1,right_1),right_1,left_2)  */
}nni;

/*********************************************************/

typedef struct __SPR{
  struct __Node         *n_link;
  struct __Node  *n_opp_to_link;
  struct __Edge  *b_opp_to_link;
  struct __Edge       *b_target;
  struct __Edge  *b_init_target;
  struct __Node          **path;
  phydbl          init_target_l;
  phydbl               l0,l1,l2;
  phydbl                    lnL;
  int                depth_path;
  int                      pars;
  int                      dist;
}spr;

/*********************************************************/

typedef struct __Triplet{
  int    size;
  phydbl *F_bc;
  phydbl *F_cd;
  phydbl *F_bd;
  phydbl ****core;
  phydbl *****MW_core;//!< Added by Marcelo.
  phydbl ***p_one_site;
  phydbl ***sum_p_one_site;
  phydbl *pi_bc;
  phydbl *pi_cd;
  phydbl *pi_bd;
  struct __Eigen **MW_eigen_struct;//!< Added by Marcelo.
  struct __Eigen *eigen_struct;
  struct __Model *mod;
}triplet;

/*********************************************************/

typedef struct __Pnode{
  struct __Pnode **next;
  int weight;
  int num;
}pnode;

/*********************************************************/

typedef struct __M4 {
  int                  n_h; /* number of hidden states */
  int                  n_o; /* number of observable states  */
  int        use_cov_alpha;
  int         use_cov_free;

  phydbl          **o_mats; /* set of matrices of substitution rates across observable states */
  phydbl          *multipl; /* vector of values that multiply each o_mats matrix */
  phydbl             *o_rr; /* relative rates (symmetric) of substitution between observable states */
  phydbl             *h_rr; /* relative rates (symmetric) of substitution between hidden states */
  phydbl            *h_mat; /* matrix that describes the substitutions between hidden states (aka switches) */
  phydbl             *o_fq; /* equilibrium frequencies for the observable states */
  phydbl             *h_fq; /* equilibrium frequencies for the hidden states */
  phydbl    *h_fq_unscaled; /* unscaled equilibrium frequencies for the hidden states */
  phydbl *multipl_unscaled; /* unscaled  vector of values that multiply each o_mats matrix */

  phydbl             delta; /* switching rate */
  phydbl             alpha; /* gamma shape parameter */
}m4;

/*********************************************************/

typedef struct __Tdraw {
  phydbl          *xcoord; /* t_node coordinates on the x axis */
  phydbl          *ycoord; /* t_node coordinates on the y axis */
  phydbl        *xcoord_s; /* t_node coordinates on the x axis (scaled) */
  phydbl        *ycoord_s; /* t_node coordinates on the y axis (scaled) */
  int          page_width;
  int         page_height;
  int      tree_box_width;

  int         *cdf_mat;
  phydbl       *cdf_mat_x;
  phydbl       *cdf_mat_y;


  phydbl max_dist_to_root;
}tdraw;

/*********************************************************/

typedef struct __Trate {
  phydbl  *br_r; /* Relative substitution rate, i.e., multiplier of mean_r on each branch */
  phydbl  lexp; /* Parameter of the exponential distribution that governs the rate at which substitution between rate classes ocur */
  phydbl alpha;
  phydbl *true_t; /* true t_node times (including root node) */
  phydbl *true_r; /* true t_edge rates (on rooted tree) */
  int *n_jps;
  int *t_jps;
  phydbl *dens; /* Probability densities of mean substitution rates at the nodes */
  phydbl c_lnL; /* Prob(Br len | time stamps, model of rate evolution) */
  phydbl c_lnL_jps; /* Prob(# Jumps | time stamps, rates, model of rate evolution) */
  int adjust_rates; /* if = 1, branch rates are adjusted such that a modification of a given t_node time
		       does not modify any branch lengths */
  int use_rates; /* if = 0, branch lengths are expressed as differences between t_node times */
  phydbl *triplet;
  phydbl less_likely;
  phydbl birth_rate;
  phydbl min_rate;
  phydbl max_rate;

  phydbl clock_r; /* Mean substitution rate, i.e., 'molecular clock' rate */
  phydbl min_clock;
  phydbl max_clock;

  phydbl lbda_nu;
  phydbl min_dt;
  phydbl step_rate;
  phydbl  *nd_r;  /* Current rates at nodes and the corresponding incoming edges */
  phydbl  *old_r; /* Old t_node rates */
  phydbl  *nd_t; /* Current t_node times */
  phydbl  *old_t; /* Old t_node times */
  

  int bl_from_rt; /* if =1, branch lengths are obtained as the product of cur_r and t */
  int approx;
  int model; /* Model number */

  phydbl nu; /* Parameter of the Exponential distribution for the corresponding model */
  phydbl min_nu;
  phydbl max_nu;

  int met_within_gibbs;

  phydbl        *ml_l; /* ML t_edge lengths (rooted) */
  phydbl       *cur_l; /* Current t_edge lengths (rooted) */
  phydbl      *u_ml_l; /* ML t_edge lengths (unrooted) */
  phydbl     *u_cur_l; /* Current t_edge lengths (unrooted) */
  phydbl         *cov;
  phydbl      *invcov;
  phydbl       covdet;
  phydbl       *cov_r;
  phydbl      *mean_r;
  struct __Node **lca; /* 2-way table of common ancestral nodes for each pari of nodes */
  int       lk_approx;

  /* spare vectors */
  phydbl     *_2n_vect1;
  phydbl     *_2n_vect2;
  phydbl     *_2n_vect3;
  phydbl     *_2n_vect4;
  short int  *_2n_vect5;
  phydbl   *_2n2n_vect1;
  phydbl   *_2n2n_vect2;

  phydbl  *cond_var;
  phydbl *reg_coeff;

  phydbl  *trip_cond_cov;
  phydbl *trip_reg_coeff;

  phydbl        *t_prior;
  phydbl    *t_prior_min;
  phydbl    *t_prior_max;
  short int *t_has_prior;
  phydbl         *t_mean;


  phydbl  true_tree_size;

  phydbl           p_max;

}trate;

/*********************************************************/

typedef struct __Tmcmc {
  int run;
  int n_tot_run;
  int sample_interval;
  int acc_lexp;
  int acc_rates;
  int acc_times;
  int acc_nu;
  int n_rate_jumps;
  int randomize;
  int norm_freq;

  phydbl *dt_prop;
  phydbl *p_no_jump;
  phydbl *t_rate_jumps;
  int    *t_rank;
  phydbl *r_path;
  time_t t_beg;
  time_t t_cur;

  char *out_filename;
  FILE *out_fp_stats;
  FILE *out_fp_trees;
  FILE *out_fp_means;
  FILE *out_fp_last;

  phydbl h_times;
  phydbl h_rates;
  phydbl h_nu;
  phydbl h_clock;
  phydbl h_omega;
}tmcmc;

/*********************************************************/

typedef struct __Tpart {
  int *ns;         /* number of states for each partition (e.g., 2, 4, 3) */
  int *cum_ns;     /* cumulative number of states (e.g., 0, 2, 6) */
  int ns_max;      /* maximum number of states */
  int ns_min;      /* minimum number of states */
  int n_partitions; /* number of partitions */
  struct __Eigen    *eigen;
}part;

/*********************************************************/

typedef struct __Tnexcom {
  char *name;
  int nparm;
  int nxt_token_t;
  int cur_token_t;
  struct __Tnexparm **parm;
}nexcom;

/*********************************************************/

typedef struct __Tnexparm {
  char *name;
  char *value;
  int nxt_token_t;
  int cur_token_t;
  int (*fp)(char *, struct __Tnexparm *, struct __Option *);
  struct __Tnexcom *com;
}nexparm;

/*********************************************************/

typedef struct __OutputToken {
    char *full;                 // full name for text output
    char *tag;                  // tag name for yaml/xml/... output
    char *val;                  // value as string
    double dval;                // value as double
    int type;                   // type of token
    int format;                 // format of token
    struct __OutputToken *next; // linked list pointer
}token;

#define ROOTTKN     0
#define SCALARTKN   1
#define SETSTARTTKN 2
#define SETENDTKN   3
#define COMMENTTKN  4
#define TBLSTARTTKN 5
#define TBLENDTKN   6


#define OUTTXT      0
#define OUTYAML     1
#define OUTDARWIN   2

#define TFNUMERIC   0
#define TFSTRING    1
#define TFNUMERICNEQ 3

/*********************************************************/

typedef struct __DarwinToken {
    char *name;                 // name of option
    char *val;                  // value as string
    int type;                   // type of token
    struct __DarwinToken *next; // linked list pointer
}dtoken;

#define DROOT       0
#define DOPTION     1

/*********************************************************/

FILE* GetTreeFile(option *io);
void Plim_Binom(phydbl pH0,int N,phydbl *pinf,phydbl *psup);
t_tree *Read_Tree(char *s_tree);
void Make_All_Edges_Light(t_node *a,t_node *d);
void Make_All_Edges_Lk(t_node *a,t_node *d,t_tree *tree);
void R_rtree(char *s_tree_a, char *s_tree_d, t_node *a, t_tree *tree, int *n_int, int *n_ext);
void Clean_Multifurcation(char **subtrees,int current_deg,int end_deg);
char **Sub_Trees(char *tree,int *degree);
int Next_Par(char *s,int pos);
char *Write_Tree(t_tree *tree);
void R_wtree(t_node *pere,t_node *fils,int *available,char **s_tree, int *pos,t_tree *tree);
void Init_Tree(t_tree *tree, int n_otu);
t_edge *Make_Edge_Light(t_node *a, t_node *d, int num);
void Init_Edge_Light(t_edge *b, int num);
void Make_Edge_Dirs(t_edge *b,t_node *a,t_node *d);
void Make_Edge_Lk(t_edge *b, t_tree *tree);
t_node *Make_Node_Light(int num);
void Make_Node_Lk(t_node *n);
align **Get_Seq(option *input);
align **Get_Seq_Phylip(option *input);
align **Read_Seq_Sequential(option *io);
align **Read_Seq_Interleaved(option *io);
int Read_One_Line_Seq(align ***data,int num_otu,FILE *in);
void Uppercase(char *ch);
void Lowercase(char *ch);
calign *Compact_Data(align **data,option *input);
calign *Compact_Cdata(calign *data, option *io);
void Get_Base_Freqs(calign *data);
void Get_AA_Freqs(calign *data);
t_tree *Read_Tree_File(option *io);
t_tree *Read_Tree_File_Phylip(FILE *fp_input_tree);
void Connect_Edges_To_Nodes(t_node *a,t_node *d,t_tree *tree,int *cur);
void Exit(char *message);
void *mCalloc(int nb,size_t size);
void *mRealloc(void *p,int nb,size_t size);
/* t_tree *Make_Light_Tree_Struct(int n_otu); */
int Sort_Phydbl_Decrease(const void *a, const void *b);
void Qksort(phydbl *A, phydbl *B, int ilo,int ihi);
void Qksort_Int(int *A, int *B, int ilo,int ihi);
void Print_Site(calign *cdata,int num,int n_otu,char *sep,int stepsize);
void Print_Seq(align **data,int n_otu);
void Print_CSeq(FILE *fp,calign *cdata);
void Order_Tree_Seq(t_tree *tree,align **data);
void Order_Tree_CSeq(t_tree *tree,calign *data);
matrix *Make_Mat(int n_otu);
void Init_Mat(matrix *mat,calign *data);
void Print_Dist(matrix *mat);
void Print_Node(t_node *a,t_node *d,t_tree *tree);
void Share_Lk_Struct(t_tree *t_full,t_tree *t_empt);
void Init_Constant();
void Print_Mat(matrix *mat);
int Sort_Edges_NNI_Score(t_tree *tree, t_edge **sorted_edges, int n_elem);
void NNI(t_tree *tree, t_edge *b_fcus, int do_swap);
void Swap(t_node *a,t_node *b,t_node *c,t_node *d,t_tree *tree);
void Update_All_Partial_Lk(t_edge *b_fcus,t_tree *tree);
void Update_SubTree_Partial_Lk(t_edge *b_fcus,t_node *a,t_node *d,t_tree *tree);
calign *Make_Cseq(int n_otu, int crunch_len, int state_len, int init_len, char **sp_names, option *io);
calign *Copy_Cseq(calign *ori, option *io);
optimiz *Make_Optimiz();
int Filexists(char *filename);
FILE *Openfile(char *filename,int mode);
void Print_Fp_Out(FILE *fp_out, time_t t_beg, time_t t_end, t_tree *tree, option *input, int n_data_set, int num_rand_tree);
void Print_Fp_Out_Lines(FILE *fp_out,time_t t_beg,time_t t_end,t_tree *tree,option *input,int n_data_set);
matrix *K80_dist(calign *data,phydbl g_shape);
matrix *JC69_Dist(calign *data, model *mod);
matrix *Hamming_Dist(calign *data,model *mod);
int Is_Ambigu(char *state,int datatype,int stepsize);
void Check_Ambiguities(calign *data,int datatype,int stepsize);
int Assign_State(char *c,int datatype,int stepsize);
void Bootstrap(t_tree *tree);
void Br_Len_Involving_Invar(t_tree *tree);
void Br_Len_Not_Involving_Invar(t_tree *tree);
void Getstring_Stdin(char *file_name);
void Print_Freq(t_tree *tree);
phydbl Num_Derivatives_One_Param(phydbl(*func)(t_tree *tree),t_tree *tree,phydbl f0,phydbl *param,phydbl stepsize,phydbl *err,int precise);
int Num_Derivative_Several_Param(t_tree *tree,phydbl *param,int n_param,phydbl stepsize,phydbl(*func)(t_tree *tree),phydbl *derivatives);
int Compare_Two_States(char *state1,char *state2,int state_size);
void Copy_One_State(char *from,char *to,int state_size);
model *Make_Model_Basic();
void Make_Model_Complete(model *mod); 
model *Copy_Model(model *ori);
void Set_Defaults_Input(option *input);
void Set_Defaults_Model(model *mod);
void Set_Defaults_Optimiz(optimiz *s_opt);
void Get_Bip(t_node *a,t_node *d,t_tree *tree);
void Alloc_Bip(t_tree *tree);
int Sort_Phydbl_Increase(const void *a,const void *b);
int Sort_String(const void *a,const void *b);
void Compare_Bip(t_tree *tree1,t_tree *tree2);
void Test_Multiple_Data_Set_Format(option *input);
int Are_Compatible(char *statea,char *stateb,int stepsize,int datatype,char * alta, char * altb);
void Hide_Ambiguities(calign *data);
void Print_Site_Lk(t_tree *tree, FILE *fp);
option *Make_Input();
t_tree *Make_Tree();
void Make_All_Tree_Nodes(t_tree *tree);
void Make_All_Tree_Edges(t_tree *tree);
void Copy_Tax_Names_To_Tip_Labels(t_tree *tree, calign *data);
t_tree *Make_And_Init_Tree(calign *data);
void Connect_Edges_To_Nodes_Recur(t_node *a, t_node *d, t_tree *tree);
void Connect_One_Edge_To_Two_Nodes(t_node *a, t_node *d, t_edge *b, t_tree *tree);
t_tree *Make_Tree_From_Scratch(int n_otu, calign *data);
t_treelist *Make_Treelist(int list_size);
void Put_Subtree_In_Dead_Objects(t_node *a, t_node *d, t_tree *tree);
void Prune_Subtree(t_node *a, t_node *d, t_edge **target, t_edge **residual, t_tree *tree);
void Reassign_Node_Nums(t_node *a, t_node *d, int *curr_ext_node, int *curr_int_node, t_tree *tree);
void Copy_Tree(t_tree *ori, t_tree *cpy);
void Reassign_Edge_Nums(t_node *a, t_node *d, int *curr_br, t_tree *tree);
void Init_Node_Light(t_node *n, int num);
void Make_All_Edge_Dirs(t_node *a, t_node *d, t_tree *tree);
void Graft_Subtree(t_edge *target, t_node *link, t_edge *add_br, t_tree *tree);
int Get_Subtree_Size(t_node *a, t_node *d);
void Pull_Subtree_From_Dead_Objects(t_node *a, t_node *d, t_tree *tree);
void Make_Edge_NNI(t_edge *b);
nni *Make_NNI();
void Init_NNI(nni *a_nni);
void Insert(t_tree *tree);
void Connect_Two_Nodes(t_node *a, t_node *d);
void Get_List_Of_Target_Edges(t_node *a, t_node *d, t_edge **list, int *list_size, t_tree *tree);
void Fix_All(t_tree *tree);
void Record_Br_Len(phydbl *where, t_tree *tree);
void Restore_Br_Len(phydbl *from, t_tree *tree);
void Get_Dist_Btw_Edges(t_node *a, t_node *d, t_tree *tree);
void Detect_Polytomies(t_edge *b, phydbl l_thresh, t_tree *tree);
void Find_Mutual_Direction(t_node *n1, t_node *n2, short int *dir_n1_to_n2, short int *dir_n2_to_n1);
void Fill_Dir_Table(t_tree *tree);
void Get_List_Of_Nodes_In_Polytomy(t_node *a, t_node *d, t_node ***list, int *size_list);
void NNI_Polytomies(t_tree *tree, t_edge *b_fcus, int do_swap);
void Check_Path(t_node *a, t_node *d, t_node *target, t_tree *tree);
void Get_List_Of_Adjacent_Targets(t_node *a, t_node *d, t_node ***node_list, t_edge ***edge_list, int *list_size);
void Sort_List_Of_Adjacent_Targets(t_edge ***list, int list_size);
t_node *Common_Nodes_Btw_Two_Edges(t_edge *a, t_edge *b);
void Make_Site_Lk_Backup(t_tree *tree);
int KH_Test(phydbl *site_lk_m1, phydbl *site_lk_M2, t_tree *tree);
void Store_P_Lk(phydbl ****ori, phydbl ****cpy, t_tree *tree);
phydbl Triple_Dist(t_node *a, t_tree *tree, int approx);
void Make_Symmetric(phydbl **F, int n);
void Round_Down_Freq_Patt(phydbl **F, t_tree *tree);
phydbl Get_Sum_Of_Cells(phydbl *F, t_tree *tree);
void Divide_Cells(phydbl **F, phydbl div, t_tree *tree);
void Divide_Mat_By_Vect(phydbl **F, phydbl *vect, int size);
void Multiply_Mat_By_Vect(phydbl **F, phydbl *vect, int size);
void Triple_Dist_Recur(t_node *a, t_node *d, t_tree *tree);
triplet *Make_Triplet_Struct(model *mod);
void Fast_Br_Len(t_edge *b, t_tree *tree, int approx);
void Fast_Br_Len_Recur(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void Print_Tree(FILE *fp, t_tree *tree);
void Found_In_Subtree(t_node *a, t_node *d, t_node *target, int *match, t_tree *tree);
void Random_Tree(t_tree *tree);
void Copy_Dist(phydbl **cpy, phydbl **orig, int n);
eigen *Make_Eigen_Struct(model *mod); //!< Changed by Marcelo.
void Random_NNI(int n_moves, t_tree *tree);
void Make_Edge_Pars(t_edge *b, t_tree *tree);
void Make_Tree_Path(t_tree *tree);
void Share_Pars_Struct(t_tree *t_full, t_tree *t_empt);
void Share_Spr_Struct(t_tree *t_full, t_tree *t_empt);
void Clean_Tree_Connections(t_tree *tree);
void Print_Settings(option *input);
void Fill_Missing_Dist(matrix *mat);
void Fill_Missing_Dist_XY(int x, int y, matrix *mat);
phydbl Least_Square_Missing_Dist_XY(int x, int y, phydbl dxy, matrix *mat);
void Update_Dirs(t_tree *tree);
void Print_Banner(FILE *fp);
void Qksort_matrix(phydbl **A, int col, int ilo, int ihi);
void Check_Memory_Amount(t_tree *tree);
int Get_State_From_P_Lk(phydbl *p_lk, int pos, t_tree *tree);
int Get_State_From_P_Pars(short int *p_pars, int pos, t_tree *tree);
void Unroot_Tree(char **subtrees);
void Print_Lk(t_tree *tree, char *string);
void Print_Pars(t_tree *tree);
void Print_Lk_And_Pars(t_tree *tree);
void Check_Dirs(t_tree *tree);
void Warn_And_Exit(char *s);
void Print_Data_Set_Number(option *input, FILE *fp);
phydbl Compare_Bip_On_Existing_Edges(phydbl thresh_len, t_tree *tree1, t_tree *tree2);
void NNI_Pars(t_tree *tree, t_edge *b_fcus, int do_swap);
void Evaluate_One_Regraft_Pos_Triple(spr *move, t_tree *tree);
int Get_State_From_Ui(int ui, int datatype);
void Read_Qmat(phydbl *daa, phydbl *pi, FILE *fp);
void Traverse_Prefix_Tree(int site, int seqnum, int *patt_num, int *n_patt, align **data, option *input, pnode *n);
pnode *Create_Pnode(int size);
int Assign_State_With_Ambiguity(char *c, int datatype, int stepsize);
void Randomize_Sequence_Order(calign *data);
void Dist_To_Node_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void Add_Root(t_edge *target, t_tree *tree);
int Is_Invar(int patt_num, int stepsize, int datatype, calign *data);
void Update_Root_Pos(t_tree *tree);
void Read_Branch_Label(char *sub_part, char *full_part, t_edge *b);
void Read_Branch_Length(char *s_d, char *s_a, t_tree *tree);
void Read_Node_Name(t_node *d, char *s_tree_d, t_tree *tree);
t_tree *Generate_Random_Tree_From_Scratch(int n_otu, int rooted);
void Random_Lineage_Rates(t_node *a, t_node *d, t_edge *b, phydbl stick_prob, phydbl *rates, int curr_rate, int n_rates, t_tree *tree);
t_edge *Find_Edge_With_Label(char *label, t_tree *tree);
void Print_Square_Matrix_Generic(int n, phydbl *mat);
char Reciproc_Assign_State(int i_state, int datatype);
void Evolve(calign *data, model *mod, t_tree *tree);
int Pick_State(int n, phydbl *prob);
void Evolve_Recur(t_node *a, t_node *d, t_edge *b, int a_state, int r_class, int site_num, calign *gen_data, model *mod, t_tree *tree);
void Site_Diversity(t_tree *tree);
void Site_Diversity_Post(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void Site_Diversity_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void Subtree_Union(t_node *n, t_edge *b_fcus, t_tree *tree);
void Binary_Decomposition(int value, int *bit_vect, int size);
void Print_Diversity(FILE *fp, t_tree *tree);
void Print_Diversity_Header(FILE *fp, t_tree *tree);
void Print_Diversity_Pre(t_node *a, t_node *d, t_edge *b, FILE *fp, t_tree *tree);
void Make_New_Edge_Label(t_edge *b);
void Print_Qmat_AA(phydbl *daa, phydbl *pi);
trate *Make_Rate_Struct(t_tree *tree);
void Init_Rate_Struct(trate *rates, t_tree *tree);
phydbl Var(phydbl *x, int n);
phydbl Mean(phydbl *x, int n);
phydbl Multivariate_Kernel_Density_Estimate(phydbl *where, phydbl **x, int sample_size, int vect_size);
void Record_Model(model *ori, model *cpy);
void Best_Of_NNI_And_SPR(t_tree *tree);
int Polint(phydbl *xa, phydbl *ya, int n, phydbl x, phydbl *y, phydbl *dy);
void Reset_Model_Parameters(model *mod);
token *Print_Banner_Small(token *t);
void JF(t_tree *tree);
t_tree *Dist_And_BioNJ(calign *cdata, model *mod, option *io);
void Add_BioNJ_Branch_Lengths(t_tree *tree, calign *cdata, model *mod);
t_tree *Read_User_Tree(calign *cdata, model *mod, option *io);
void Print_Time_Info(time_t t_beg, time_t t_end);
void Prepare_Tree_For_Lk(t_tree *tree);
char *Bootstrap_From_String(char *s_tree, calign *cdata, model *mod, option *io);
char *aLRT_From_String(char *s_tree, calign *cdata, model *mod, option *io);
int Rand_Int(int min, int max);
void PhyML_Printf(char *format, ...);
void PhyML_Fprintf(FILE *fp, char *format, ...);

//Added/Updated by Ken 3/8/2016
void Update_Ancestors(t_node *a, t_node *d, t_tree *tree);
void Update_Ancestors_Edge(t_node *a, t_node *d, t_edge *e, t_tree *tree);
void Update_Ancestors_Edge_Prune(t_node *a, t_node *d, t_edge *e, t_node *stop, int prunedir, t_tree *tree);

void Down_All(t_node *d, t_tree *tree);
void Up_All(t_node *d, t_tree *tree);
void Up_Above(t_node *d, t_tree *tree);
void Down_Above(t_node *d, t_tree *tree);
void Node_Stats(t_node *a,t_node *d, t_tree *tree);
void Print_Ambig_States(t_node *d, t_tree *tree, FILE *ambigfile);


void Find_Common_Tips(t_tree *tree1, t_tree *tree2);
phydbl Get_Tree_Size(t_tree *tree);
int Find_Bipartition(char **target_bip, int bip_size, t_tree *tree);
int Find_Duplication_Node(char **tax_set, int n_tax, t_tree *tree);
void Get_Rid_Of_Prefix(char delim, t_tree *tree);
void Dist_To_Root_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void Dist_To_Root(t_node *n_root, t_tree *tree);
int Is_Duplication_Node(t_node *n, char **tax_set, int n_tax, t_tree *tree);
int Sort_Edges_Depth(t_tree *tree, t_edge **sorted_edges, int n_elem);
char *Basename(char *path);
void Get_List_Of_Ancestors(t_node *ref_nod, t_node **list, int *size, t_tree *tree);
t_node *Find_Lca(t_node *n1, t_node *n2, t_tree *tree);
int Edge_Num_To_Node_Num(int edge_num, t_tree *tree);
void Branch_Lengths_To_Time_Lengths(t_tree *tree);
void Branch_Lengths_To_Time_Lengths_Pre(t_node *a, t_node *d, t_tree *tree);
int Find_Clade(char **tax_name_list, int list_size, t_tree *tree);
void Find_Clade_Pre(t_node *a, t_node *d, int *tax_num_list, int list_size, int *num, t_tree *tree);
void Read_Clade_Priors(char *file_name, t_tree *tree);
t_edge *Find_Root_Edge(FILE *fp_input_tree, t_tree *tree);
void Copy_Tree_Topology_With_Labels(t_tree *ori, t_tree *cpy);
void Make_Custom_Model(model *mod);
option *Get_Input(int argc, char **argv);
nexcom **Make_Nexus_Com();
void Free_Nexus_Com(nexcom **com);
void Init_Nexus_Format(nexcom **com);
void Find_Nexus_Com(char *token, nexcom **found_com, nexparm **default_parm, nexcom **com_list);
void Find_Nexus_Parm(char *token, nexparm **found_parm, nexcom *curr_com);
int Read_Nexus_Dimensions(char *token, nexparm *curr_parm, option *io);
int Read_Nexus_Format(char *token, nexparm *curr_parm, option *io);
int Read_Nexus_Eliminate(char *token, nexparm *curr_parm, option *io);
int Read_Nexus_Taxlabel(char *token, nexparm *curr_parm, option *io);
int Read_Nexus_Charstatelabels(char *token, nexparm *curr_parm, option *io);
int Read_Nexus_Charlabels(char *token, nexparm *curr_parm, option *io);
int Read_Nexus_Statelabels(char *token, nexparm *curr_parm, option *io);
int Read_Nexus_Matrix(char *token, nexparm *curr_parm, option *io);
int Read_Nexus_Begin(char *token, nexparm *curr_parm, option *io);
int Read_Nexus_Matrix(char *token, nexparm *curr_parm, option *io);
int Read_Nexus_Translate(char *token, nexparm *curr_parm, option *io);
int Read_Nexus_Tree(char *token, nexparm *curr_parm, option *io);
int Read_Nexus_Taxa(char *token, nexparm *curr_parm, option *io);

align** Read_Seq_Fasta(option*,model*);
void scanFasta(model*,FILE*);

void Detect_Align_File_Format(option *io);
void Detect_Tree_File_Format(option *io);
int Get_Token(FILE *fp, char *token);
nexparm *Make_Nexus_Parm();
void Free_Nexus_Parm(nexparm *parm);
void Read_Ntax_Len_Phylip(FILE *fp ,int *n_otu, int *n_tax);
void Set_Model_Name(model *mod);
void Adjust_Min_Diff_Lk(t_tree *tree);
void Match_Tip_Numbers(t_tree *tree1, t_tree *tree2);
void Get_Nexus_Data(FILE *fp, option *io);
void Translate_Tax_Names(char **tax_names, t_tree *tree);
void Skip_Comment(FILE *fp);
void Update_Dir_To_Tips(t_node *a, t_node *d, t_tree *tree);
void Test_Node_Table_Consistency(t_tree *tree);

void Genetic_code_index_stopCodons(int gencode); //!< Added by Marcelo.
int  Genetic_code_ns(); //!< Added by Marcelo.
int  Intersect_Site_Alternatives(align **data, int n_otu, int site); //!< Added by Marcelo.
void  CopyExtraFieldsCodon(align *from, int site, align *to,  int num_pattern); //!< Added by Marcelo.
void  Get_Base_Freqs_CODONS_FaXb(calign *cdata, align **data, int freqMod, model *mod); //!< Added by Marcelo.
void  Get_Base_Freqs_CODONS_F1XCODONS(calign *cdata, model *mod); //!< Added by Marcelo.
int CompareCodonStr(char *a,char *b,int len); //!< Added by Marcelo.
void CF3x4(phydbl * freq, int gencode); //!< Added by Marcelo.
phydbl my2norm(phydbl *x,int len); //!< Added by Marcelo.
int MyGaussElimination_gNxRHS(double *A, double *B, int n, int rhs); //!< Added by Marcelo.
void Sprint_codon(char *inout,int codon); //!< Added by Marcelo.
void Make_Ts_and_Tv_Matrix(option *io); //!< Added by Marcelo.
void Change_Tree_Data(int n_otu, calign *data, t_tree * tree); //!< Added by Marcelo.
void Scale_freqs(phydbl *f, int n); //!< Added by Marcelo.
void Freq_to_UnsFreq(phydbl *f, phydbl *uf, int n, int f2U); //!< Added by Marcelo.
void lkExperiment(t_tree * tree, int num_tree); //!< Added by Marcelo.
void Scale_freqs_tol(phydbl *f, int n, phydbl tolmin, phydbl tolmax); //!< Added by Marcelo.
void rel_lk_dimension_reduction(t_tree *tree); //!< Added by Marcelo.
int findAA(char theAA);//!< Added by Marcelo.
void Conv_CD_seq_to_AA_seq(option *io_codon, option *io_aa, option *io);//!< Added by Marcelo.
void Conv_NT_seq_to_AA_seq(option *io);//!< Added by Marcelo.
char Translate_CD_to_AA(char * codon, option * io);//!< Added by Marcelo.
char Translate_CD_to_AA_2(char codon, option * io);//!< Added by Marcelo.
void SetSeed (unsigned int seed);//!< Added by Marcelo.
double rndu (void);//!< Added by Marcelo.
matrix * JC69_Dist_Codon(calign *data, model *mod);//!< Added by Marcelo.
void Read_userRatesAndFreqs( phydbl *daa, phydbl *pi, int ns, FILE *fp ); //! Added by Marcelo.
void Read_userRatesAndFreqsMG( phydbl *daa, phydbl *pi, phydbl *bfreqs, int nbfreq, int ns, FILE *fp ); //! Added by Marcelo
int myFactorial(int n); //! Added by Marcelo.
char findAA_Name(int n); //! Added by Marcelo.

FILE * openOutputFile( char * outFile, char * appendString, char * suffix, option *io );

extern char     aminoAcidmap[65]; //!< Added by Marcelo.
extern int        stopCodons[64]; //!< Added by Marcelo.
extern int       senseCodons[64]; //!< Added by Marcelo.
extern int  indexSenseCodons[64]; //!< Added by Marcelo.

// added by stefan, output format support:
token *Emit_Root_Token();
token *Emit_Out_Token(token *currTkn, char *full, char *tag, int type, int f, char *format, ...);
token *Del_Last_Token(token *r);
void Output_Tokens(token *root, int format, FILE *out);
void Print_Txt_Tokenlist(FILE *fp_out, token *t);
#ifdef USEYAML
void Print_YAML_Tokenlist(FILE *fp_out, token *t);
#endif
void Print_Darwin_Tokenlist(FILE *fp_out, token *t);
char *getDarwinTable(token *r);
// end added stefan

int Get_Genetic_Code();
void Set_Genetic_Code(int gencode);


#include "nucle2codon.h" //!< Added by Marcelo.

#define ABS(x) (x < 0 ? -(x) : (x)) //!< Added by Marcelo.

/*! Codon Models.*/ //!< Added by Marcelo.
#define GY -1
#define MG -2
#define YAP -3
#define PCM -4
#define HLP17 -100 //added by Ken

/*! Opt freq scheme */ //!
#define NOFEQ     -1
#define USER      3
#define OPTIMIZE  2
#define MODEL     1
#define EMPIRICAL 0

/*! Parametric Models (GY, MG, YAP) */
#define GYSIMP -1
#define MGSIMP -2
#define YAPSIMP -3

/*! Kosiol 2007 Empirical and Semi-Empirical with GY */
#define GYECMK07 -4
#define GYECMK07WK -5
#define GYECMK07F -6
#define GYECMK07WKF -7

/*! Schneider 2005 Empircal and Semi-Empirical with GY */ 
#define GYECMS05 -8
#define GYECMS05WK -9
#define GYECMS05F -10
#define GYECMS05WKF -11

/*! Kosiol 2007 Empirical and Semi-Empirical with MG */
#define MGECMK07 -12
#define MGECMK07WK -13
#define MGECMK07F -14
#define MGECMK07WKF -15

/*! Schneider 2005 Empircal and Semi-Empirical with MG */
#define MGECMS05 -16
#define MGECMS05WK -17
#define MGECMS05F -18
#define MGECMS05WKF -19

/*! Kosiol 2007 Empirical and Semi-Empirical with YAP */
#define YAPECMK07 -20
#define YAPECMK07WK -21
#define YAPECMK07F -22
#define YAPECMK07WKF -23

/*! Schneider 2005 Empircal and Semi-Empirical with YAP */
#define YAPECMS05 -24
#define YAPECMS05WK -25
#define YAPECMS05F -26
#define YAPECMS05WKF -27

/*! user empirical matrice, GY, MG and YAP */
#define YAPECMUSR -28
#define YAPECMUSRWK -29
#define YAPECMUSRF -30
#define YAPECMUSRWKF -31
#define GYECMUSR -32
#define GYECMUSRWK -33
#define GYECMUSRF -34
#define GYECMUSRWKF -35
#define MGECMUSRF -36
#define MGECMUSRWKF -37
#define MGECMUSR -38
#define MGECMUSRWK -39
//
/*! Genetic code.*/ //!< Added by Marcelo.
#define  NOGENCODE -1
#define  STANDARD   0   
#define  TVMC       1
#define  TYMC       2
#define  THMPCMCMSC 3
#define  THIMC      4
#define  THCDHNC    5
#define  THEFMC     6
#define  THENC      7
#define  THBAPPC    8
#define  THAYNC     9
#define  THAMC      10
#define  THAFMC     11
#define  BLNC       12
#define  CHMC       13
#define  TRMC       14
#define  SCOMC      15
#define  THMC       16

/*! Equilibrium Frequencies.*/ //!< Added by Marcelo.
#define FUNDEFINED    -1
#define F1XSENSECODONS 0                                                               
#define F1X4           1                                                               
#define F3X4           2                                                             
#define CF3X4          3                                                               
#define FMODEL         4


/*! Codon model nature.*/
#define NOCHOICE -1
#define PARA 0
#define EMPI 1
#define SEMI 2
#define PCA  3

/*! Data type.*/ //!< Added by Marcelo. 
#define CODON 3                                                                    

/*! Initial Distance Tree and/or initial rate matrices.*/ //!< Added by Marcelo.                                                            
#define NOINITMAT -1
#define NUCLEO     0
#define KOSI07     1
#define SCHN05     2
#define ECMUSR     3

/*! Opt Method.*/ //!< Added by Marcelo.                                                            
#define PAML     0
#define PHYML     1


/*! kappa models according to Kosiol 2007.*/ //!< Added by Marcelo.   
#define kap1 0
#define kap2 1
#define kap3 2
#define kap4 3
#define kap5 4
#define kap6 5

/*! Omega models of site variability.*/ //!< Added by Marcelo.
#define NOOMEGA -1
#define DM0      0
#define DMODELK  1
#define DGAMMAK  2

/*! Numerical Libraries.*/ //!< Added by Marcelo.
#if defined BLAS || defined BLAS_OMP
#if defined MACOS
#include <Accelerate/Accelerate.h>
#else
#include <atlas/cblas.h> //edited by Ken 12/1/2017
/*#include "atlas_aux.h"*/ //cut out by Ken 12/1
#include <atlas/clapack.h> //edited by Ken 12/1/2017
#endif
#endif

#if defined OMP || defined BLAS_OMP
#include <omp.h>
#endif

/*!< Exp Matrix Calculation.*/ //!< Added by Marcelo.
#define EIGEN  0
#define SSPADE 1
#define TAYLOR 2

/*!< Opt method .*/
#define optPAML 0
#define optPHYML 1

#include "free.h"
#include "spr.h"
#include "lk.h"
#include "optimiz.h"
#include "models.h"
#include "bionj.h"
#include "simu.h"
#include "eigen.h"
#include "pars.h"
#include "alrt.h"
#include "stats.h"
#include "help.h"

/* unused ken 5/1
#ifdef MPI
#include "mpi_boot.h"
#endif

#ifdef MG
#include "mg.h"
#endif

#ifdef TIME
#include "times.h"
#include "rates.h"
#endif

#ifdef _NOT_NEEDED_A_PRIORI
#include "m4.h"
#endif
*/

#endif
