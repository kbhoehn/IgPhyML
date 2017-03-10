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
#include "models.h"
#include <stdio.h>


//!< Start added by Marcelo.
#include "pcm.h"
#include "ecm.h"

extern int     stopCodons[64]; 
extern int    senseCodons[64]; 
extern char  aminoAcidmap[65]; 


//!< End added by Marcelo.

/*********************************************************/
/* Handle any number of states (>1) */
/*! Jukes Cantor */
void PMat_JC69(phydbl l, int pos, phydbl *Pij, model *mod)
{
    int ns;
    int i,j;
    
    ns = mod->ns;
    
    For(i,ns) Pij[pos+ ns*i+i] = 1. - ((ns - 1.)/ns)*(1. - EXP(-ns*l/(ns - 1.)));
    For(i,ns-1) 
    for(j=i+1;j<ns;j++) 
    {
        Pij[pos+ ns*i+j] = (1./ns)*(1. - EXP(-ns*l/(ns - 1.)));
        if(Pij[pos+ns*i+j] < SMALL_PIJ) Pij[pos+ns*i+j] = SMALL_PIJ;
        Pij[pos+ ns*j+i] = Pij[pos+ ns*i+j];
        
    }
}
/*********************************************************/
/*! Kimura 2-parameter (and JC) */
void PMat_K80(phydbl l, phydbl kappa, int pos, phydbl *Pij)
{
    phydbl Ts,Tv,e1,e2,aux;
    int i,j;
    /*0 => A*/
    /*1 => C*/
    /*2 => G*/
    /*3 => T*/
    
    /* Ts -> transition*/
    /* Tv -> transversion*/
    
    aux = -2*l/(kappa+2);
    e1 = (phydbl)EXP(aux *2);
    
    e2 = (phydbl)EXP(aux *(kappa+1));
    Tv = .25*(1-e1);
    Ts = .25*(1+e1-2*e2);
    
    Pij[pos+ 4*0+0] = Pij[pos+ 4*1+1] = 
    Pij[pos+ 4*2+2] = Pij[pos+ 4*3+3] = 1.-Ts-2.*Tv;
    
    Pij[pos+ 4*0+1] = Pij[pos+ 4*1+0] = Tv;
    Pij[pos+ 4*0+2] = Pij[pos+ 4*2+0] = Ts;
    Pij[pos+ 4*0+3] = Pij[pos+ 4*3+0] = Tv;
    
    Pij[pos+ 4*1+2] = Pij[pos+ 4*2+1] = Tv;
    Pij[pos+ 4*1+3] = Pij[pos+ 4*3+1] = Ts;
    
    Pij[pos+ 4*2+3] = Pij[pos+ 4*3+2] = Tv;
    
    For(i,4) For(j,4)
    if(Pij[pos + 4*i+j] < SMALL_PIJ) Pij[pos + 4*i+j] = SMALL_PIJ;
    
}

/*********************************************************/

/*! Tamura Nei 93 (and Felsenstein 81, 84 and HKY85) */
void PMat_TN93(phydbl l, model *mod, int pos, phydbl *Pij)
{
    int i,j;
    phydbl e1,e2,e3;
    phydbl a1t,a2t,bt;
    phydbl A,C,G,T,R,Y;
    phydbl kappa1,kappa2;
    int kappa_has_changed;
    
    A = mod->pi[0]; C = mod->pi[1]; G = mod->pi[2]; T = mod->pi[3];
    R = A+G;  Y = T+C;
    
    kappa_has_changed = 0;
    if(mod->kappa < .0) mod->kappa = 1.0e-5;
    
    if((mod->whichmodel != F84) && (mod->whichmodel != TN93)) mod->lambda = 1.; 
    else if(mod->whichmodel == F84)
    {
        do
        {
            mod->lambda = (Y+(R-Y)/(2.*mod->kappa))/(R-(R-Y)/(2.*mod->kappa));
            
            if(mod->lambda < .0)
            {
                mod->kappa += mod->kappa/10.;
                kappa_has_changed = 1;
            }
        }while(mod->lambda < .0);
    }
    
    
    if((!mod->s_opt->opt_kappa) && (kappa_has_changed))
    {
        PhyML_Printf("\n. WARNING: This transition/transversion ratio\n");
        PhyML_Printf("  is impossible with these base frequencies!\n");
        PhyML_Printf("  The ratio is now set to %.3f\n",mod->kappa);
    }
    
    kappa2 = mod->kappa*2./(1.+mod->lambda);
    kappa1 = kappa2 * mod->lambda;
    
    
    bt = l/(2.*(A*G*kappa1+C*T*kappa2+R*Y));
    
    a1t = kappa1;
    a2t = kappa2;
    a1t*=bt; a2t*=bt;
    
    e1 = (phydbl)EXP(-a1t*R-bt*Y);
    e2 = (phydbl)EXP(-a2t*Y-bt*R);
    e3 = (phydbl)EXP(-bt);
    
    
    /*A->A*/Pij[pos + 4*0+0] = A+Y*A/R*e3+G/R*e1; 
    /*A->C*/Pij[pos + 4*0+1] = C*(1-e3);
    /*A->G*/Pij[pos + 4*0+2] = G+Y*G/R*e3-G/R*e1;
    /*A->T*/Pij[pos + 4*0+3] = T*(1-e3);
    
    /*C->A*/Pij[pos + 4*1+0] = A*(1-e3);
    /*C->C*/Pij[pos + 4*1+1] = C+R*C/Y*e3+T/Y*e2;
    /*C->G*/Pij[pos + 4*1+2] = G*(1-e3);
    /*C->T*/Pij[pos + 4*1+3] = T+R*T/Y*e3-T/Y*e2;
    
    /*G->A*/Pij[pos + 4*2+0] = A+Y*A/R*e3-A/R*e1;
    /*G->C*/Pij[pos + 4*2+1] = C*(1-e3);
    /*G->G*/Pij[pos + 4*2+2] = G+Y*G/R*e3+A/R*e1;
    /*G->T*/Pij[pos + 4*2+3] = T*(1-e3);
    
    /*T->A*/Pij[pos + 4*3+0] = A*(1-e3);
    /*T->C*/Pij[pos + 4*3+1] = C+R*C/Y*e3-C/Y*e2;
    /*T->G*/Pij[pos + 4*3+2] = G*(1-e3);
    /*T->T*/Pij[pos + 4*3+3] = T+R*T/Y*e3+C/Y*e2;
    
    For(i,4) For(j,4)
    if(Pij[pos + 4*i+j] < SMALL_PIJ) Pij[pos + 4*i+j] = SMALL_PIJ;
    
    /*   /\*A->A*\/(*Pij)[0][0] = A+Y*A/R*e3+G/R*e1;  */
    /*   /\*A->C*\/(*Pij)[0][1] = C*(1-e3); */
    /*   /\*A->G*\/(*Pij)[0][2] = G+Y*G/R*e3-G/R*e1; */
    /*   /\*A->T*\/(*Pij)[0][3] = T*(1-e3); */
    
    /*   /\*C->A*\/(*Pij)[1][0] = A*(1-e3); */
    /*   /\*C->C*\/(*Pij)[1][1] = C+R*C/Y*e3+T/Y*e2; */
    /*   /\*C->G*\/(*Pij)[1][2] = G*(1-e3); */
    /*   /\*C->T*\/(*Pij)[1][3] = T+R*T/Y*e3-T/Y*e2; */
    
    /*   /\*G->A*\/(*Pij)[2][0] = A+Y*A/R*e3-A/R*e1; */
    /*   /\*G->C*\/(*Pij)[2][1] = C*(1-e3); */
    /*   /\*G->G*\/(*Pij)[2][2] = G+Y*G/R*e3+A/R*e1; */
    /*   /\*G->T*\/(*Pij)[2][3] = T*(1-e3); */
    
    /*   /\*T->A*\/(*Pij)[3][0] = A*(1-e3); */
    /*   /\*T->C*\/(*Pij)[3][1] = C+R*C/Y*e3-C/Y*e2; */
    /*   /\*T->G*\/(*Pij)[3][2] = G*(1-e3); */
    /*   /\*T->T*\/(*Pij)[3][3] = T+R*T/Y*e3+C/Y*e2; */
    
    /*   For(i,4) For(j,4) */
    /*     if((*Pij)[i][j] < SMALL) (*Pij)[i][j] = SMALL; */
    
}

/*********************************************************/


/********************************************************************/
/* void PMat_Empirical(phydbl l, model *mod, phydbl ***Pij)         */
/*                                                                  */
/* Computes the substitution probability matrix                     */
/* from the initial substitution rate matrix and frequency vector   */
/* and one specific branch length                                   */
/*                                                                  */
/* input : l , branch length                                        */
/* input : mod , choosen model parameters, qmat and pi             */
/* ouput : Pij , substitution probability matrix                    */
/*                                                                  */
/* matrix P(l) is computed as follows :                             */
/* P(l) = EXP(Q*t) , where :                                        */
/*                                                                  */
/*   Q = substitution rate matrix = Vr*D*inverse(Vr) , where :      */
/*                                                                  */
/*     Vr = right eigenvector matrix for Q                          */
/*     D  = diagonal matrix of eigenvalues for Q                    */
/*                                                                  */
/*   t = time interval = l / mr , where :                           */
/*                                                                  */
/*     mr = mean rate = branch length/time interval                 */
/*        = sum(i)(pi[i]*p(i->j)) , where :                         */
/*                                                                  */
/*       pi = state frequency vector                                */
/*       p(i->j) = subst. probability from i to a different state   */
/*               = -Q[ii] , as sum(j)(Q[ij]) +Q[ii] =0              */
/*                                                                  */
/* the Taylor development of EXP(Q*t) gives :                       */
/* P(l) = Vr*EXP(D*t)        *inverse(Vr)                           */
/*      = Vr*POW(EXP(D/mr),l)*inverse(Vr)                           */
/*                                                                  */
/* for performance we compute only once the following matrixes :    */
/* Vr, inverse(Vr), EXP(D/mr)                                       */
/* thus each time we compute P(l) we only have to :                 */
/* make 20 times the operation POW()                                */
/* make 2 20x20 matrix multiplications , that is :                  */
/*   16000 = 2x20x20x20 times the operation *                       */
/*   16000 = 2x20x20x20 times the operation +                       */
/*   which can be reduced to (the central matrix being diagonal) :  */
/*   8400 = 20x20 + 20x20x20 times the operation *                  */
/*   8000 = 20x20x20 times the operation +                          */
/********************************************************************/

/*! Empirical nucleotide model */
void PMat_Empirical(phydbl l, model *mod, int pos, phydbl *Pij)
{
    int n = mod->ns;
    int i, j, k;
    phydbl *U,*V,*R;
    phydbl *expt; 
    phydbl *uexpt;
    
    expt  = mod->eigen->e_val_im;
    uexpt = mod->eigen->r_e_vect_im;
    U     = mod->eigen->r_e_vect;
    V     = mod->eigen->l_e_vect;
    R     = mod->eigen->e_val; /* exponential of the eigen value matrix */
    
    For(i,n) For(k,n) Pij[pos+mod->ns*i+k] = .0;
    
    /* compute POW(EXP(D/mr),l) into mat_eDmrl */
    For(k,n) expt[k] = (phydbl)POW(R[k],l);
    
    /* multiply Vr*POW(EXP(D/mr),l)*Vi into Pij */
    For (i,n) For (k,n) uexpt[i*n+k] = U[i*n+k] * expt[k];
    
    For (i,n) 
    {
        For (j,n) 
        {
            For(k,n)
            {
                Pij[pos+mod->ns*i+j] += (uexpt[i*n+k] * V[k*n+j]);
            }
            /* 	  if(Pij[pos+mod->ns*i+j] < SMALL) Pij[pos+mod->ns*i+j] = SMALL; */
            if(Pij[pos+mod->ns*i+j] < SMALL_PIJ) Pij[pos+mod->ns*i+j] = SMALL_PIJ;
        }
        
#ifndef PHYML
        phydbl sum;
        sum = .0;
        For (j,n) sum += Pij[pos+mod->ns*i+j];
        if((sum > 1.+.0001) || (sum < 1.-.0001))
        {
            PhyML_Printf("\n");
            PhyML_Printf("\n. Q\n");
            For(i,n) { For(j,n) PhyML_Printf("%7.3f ",mod->eigen->q[i*n+j]); PhyML_Printf("\n"); }
            PhyML_Printf("\n. U\n");
            For(i,n) { For(j,n) PhyML_Printf("%7.3f ",U[i*n+j]); PhyML_Printf("\n"); }
            PhyML_Printf("\n");
            PhyML_Printf("\n. V\n");
            For(i,n) { For(j,n) PhyML_Printf("%7.3f ",V[i*n+j]); PhyML_Printf("\n"); }
            PhyML_Printf("\n");
            PhyML_Printf("\n. Eigen\n");
            For(i,n)  PhyML_Printf("%E ",expt[i]);
            PhyML_Printf("\n");
            PhyML_Printf("\n. Pij\n");
            For(i,n) { For (j,n) PhyML_Printf("%f ",Pij[pos+mod->ns*i+j]); PhyML_Printf("\n"); }
            PhyML_Printf("\n. sum = %f",sum);
            if(mod->m4mod)
            {
                int i;
                PhyML_Printf("\n. mod->m4mod->alpha = %f",mod->m4mod->alpha);
                PhyML_Printf("\n. mod->m4mod->delta = %f",mod->m4mod->delta);
                For(i,mod->m4mod->n_h)
                {
                    PhyML_Printf("\n. mod->m4mod->multipl[%d] = %f",i,mod->m4mod->multipl[i]);
                }
            }
            PhyML_Printf("\n. l=%f",l);
            PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
            Warn_And_Exit("");
        }
#endif
    }
}

/*********************************************************/
void PMat_Gamma(phydbl l, model *mod, int pos, phydbl *Pij)
{
    int n;
    int i, j, k;
    phydbl *U,*V,*R;
    phydbl *expt; 
    phydbl *uexpt;
    phydbl shape;
    
    
    n     = mod->ns;
    expt  = mod->eigen->e_val_im;
    uexpt = mod->eigen->r_e_vect_im;
    U     = mod->eigen->r_e_vect;
    V     = mod->eigen->l_e_vect;
    R     = mod->eigen->e_val; /* exponential of the eigen value matrix */
    
    if(mod->n_catg == 1) shape = 1.E+4;
    else                 shape = mod->alpha;
    
    
    For(i,n) For(k,n) Pij[pos+mod->ns*i+k] = .0;
    
    if(shape < 1.E-10) 
    {
        PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
        Warn_And_Exit("");
    }
    
    /* Formula 13.42, page 220 of Felsenstein's book ``Inferring Phylogenies'' */ 
    For(k,n) expt[k] = POW(shape/(shape-LOG(R[k])*l),shape);
    
    /* multiply Vr*expt*Vi into Pij */
    For(i,n) For(k,n) uexpt[i*n+k] = U[i*n+k] * expt[k];
    
    For (i,n) 
    {
        For (j,n) 
        {
            For(k,n)
            {
                Pij[pos+mod->ns*i+j] += (uexpt[i*n+k] * V[k*n+j]);
            }
            if(Pij[pos+mod->ns*i+j] < SMALL_PIJ) Pij[pos+mod->ns*i+j] = SMALL_PIJ;
        }
        
#ifdef DEBUG
        phydbl sum;
        sum = .0;
        For (j,n) sum += Pij[pos+mod->ns*i+j];
        if((sum > 1.+.0001) || (sum < 1.-.0001))
        {
            PhyML_Printf("\n");
            PhyML_Printf("\n. Q\n");
            For(i,n) { For(j,n) PhyML_Printf("%7.3f ",mod->eigen->q[i*n+j]); PhyML_Printf("\n"); }
            PhyML_Printf("\n. U\n");
            For(i,n) { For(j,n) PhyML_Printf("%7.3f ",U[i*n+j]); PhyML_Printf("\n"); }
            PhyML_Printf("\n");
            PhyML_Printf("\n. V\n");
            For(i,n) { For(j,n) PhyML_Printf("%7.3f ",V[i*n+j]); PhyML_Printf("\n"); }
            PhyML_Printf("\n");
            PhyML_Printf("\n. Eigen\n");
            For(i,n)  PhyML_Printf("%E ",expt[i]);
            PhyML_Printf("\n");
            PhyML_Printf("\n. Pij\n");
            For(i,n) { For (j,n) PhyML_Printf("%f ",Pij[pos+mod->ns*i+j]); PhyML_Printf("\n"); }
            PhyML_Printf("\n. sum = %f",sum);
            if(mod->m4mod)
            {
                int i;
                PhyML_Printf("\n. mod->m4mod->alpha = %f",mod->m4mod->alpha);
                PhyML_Printf("\n. mod->m4mod->delta = %f",mod->m4mod->delta);
                For(i,mod->m4mod->n_h)
                {
                    PhyML_Printf("\n. mod->m4mod->multipl[%d] = %f",i,mod->m4mod->multipl[i]);
                }
            }
            PhyML_Printf("\n. l=%f",l);
            PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
            Warn_And_Exit("");
        }
#endif
    }
}

/*********************************************************/

void PMat_Zero_Br_Len(model  *mod, int pos, phydbl *Pij)
{
    int n = mod->ns;
    int i;
    // int j;
    
    For (i,n*n) Pij[pos+i] = 0.0;  //!< Changed by Marcelo.
    For(i,n) Pij[pos+n*i+i] = 1.0;//!< Changed by Marcelo.
    
}

/*********************************************************/

void PMat(phydbl l, model *mod, int pos, phydbl *Pij)
{
    if(l < BL_MIN-POW(2,-27))
    {
        PMat_Zero_Br_Len(mod,pos,Pij);
    }
    else
    {
        switch(mod->io->datatype)
        {
            case NT :
            {
                if(mod->use_m4mod)
                {
                    PMat_Empirical(l,mod,pos,Pij);
                }
                else
                {
                    if((mod->whichmodel == JC69) ||  
                       (mod->whichmodel == K80))  
                    {
                        /* 		    PMat_JC69(l,pos,Pij,mod); */
                        PMat_K80(l,mod->kappa,pos,Pij);
                    }
                    else
                    {
                        if(
                           (mod->whichmodel == F81)   ||
                           (mod->whichmodel == HKY85) ||
                           (mod->whichmodel == F84)   ||
                           (mod->whichmodel == TN93))
                        {
                            PMat_TN93(l,mod,pos,Pij);
                        }
                        else
                        {
                            PMat_Empirical(l,mod,pos,Pij);
                        }
                    }
                    break;
                }
            case AA : 
                {
                    PMat_Empirical(l,mod,pos,Pij);
                    break;
                }
                
            case CODON ://!< Added by Marcelo.
                {
                    if(mod->calculate_init_tree && mod->io->init_DistanceTreeCD==NUCLEO) PMat_JC69(l,pos,Pij,mod);
                    else PMat_CODON(l,mod,0,Pij);
                    break;
                }
                
            default:
                {
                    PMat_JC69(l,pos,Pij,mod);
                    break;
                    /* 	      PhyML_Printf("\n. Not implemented yet.\n"); */
                    /* 	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
                    /* 	      Warn_And_Exit(""); */
                    /* 	      break; */
                }
            } //incorrect right curly bracket? Whole thing is case NT: CHECK, Ken 4/1/2017
        }
    }
}

/*********************************************************/
/*********************************************************/

void Init_Model(calign *data, model *mod, option *io)
{
    int i,j;
    phydbl sum, aux, mr;
    int result;
    
    if(io->datatype == CODON) { //!<Added by Marcelo.
        if(!mod->invar) {
            For(i, data->crunch_len) {
                data->invar[i] = 0;
            }
        }
        
        if(mod->s_opt->opt_pinvar) {
            mod->pinvar = 0.2;
        }
        
        /*! Omega variation model */
        switch(mod->omegaSiteVar) {
            case DM0: {
                mod->alpha = 1.0;
                mod->omegas[0] = 1.0;
                mod->prob_omegas[0] = 1.0;
                DiscreteGamma(mod->gamma_r_proba, mod->gamma_rr, mod->alpha, mod->alpha, mod->n_catg, mod->gamma_median);
                break;
            }
            case DGAMMAK: {
                DiscreteGamma(mod->prob_omegas, mod->omegas, mod->alpha, mod->beta, mod->n_w_catg, mod->gamma_median);
                For(i, mod->n_w_catg) {
                    if(mod->omegas[i] < MODELPAREPS) mod->omegas[i] = MODELPAREPS;
                }
                break;
            }
            case DMODELK: {
                Scale_freqs_tol(mod->prob_omegas, mod->n_w_catg, MODELPAREPS, 0.99);
                Freq_to_UnsFreq(mod->prob_omegas, mod->prob_omegas_uns, mod->n_w_catg, 1);
                break;
            }
            default:
                break; 
        } 
        
        if((mod->initqrates != NOINITMAT) && (mod->omegaSiteVar != NOOMEGA) && (io->kappaECM == kap5)) {
            Scale_freqs(mod->pkappa, mod->nkappa);
            Freq_to_UnsFreq(mod->pkappa, mod->unspkappa, mod->nkappa, 1);
        }
        if((mod->freq_model != FUNDEFINED) && (mod->freq_model != FMODEL)) {
            if(mod->s_opt->user_state_freq) {
                mod->s_opt->opt_state_freq = NO;
                EqFrequencies(mod->freq_model, mod->pi, mod->user_b_freq, mod->ns);
            } else {
                EqFrequencies(mod->freq_model, mod->pi, mod->base_freq, mod->ns);
                switch(mod->freq_model) {
                    case F1X4: {
                        Scale_freqs(mod->base_freq, mod->num_base_freq);
                        Freq_to_UnsFreq(mod->base_freq, mod->uns_base_freq, mod->num_base_freq, 1);
                        break;
                    }
                    case F3X4:
                    case CF3X4: {
                        Scale_freqs(mod->base_freq, 4);
                        Freq_to_UnsFreq(mod->base_freq, mod->uns_base_freq, 4, 1);
                        
                        Scale_freqs(mod->base_freq+4, 4);
                        Freq_to_UnsFreq(mod->base_freq+4, mod->uns_base_freq+4, 4, 1);
                        
                        Scale_freqs(mod->base_freq+8, 4);
                        Freq_to_UnsFreq(mod->base_freq+8, mod->uns_base_freq+8, 4, 1);
                        
                        break;
                    }
                    default:
                        break;
                } 
            }
            Scale_freqs(mod->pi, mod->ns);
            Freq_to_UnsFreq(mod->pi, mod->pi_unscaled, mod->ns, 1);
        } else { 
            if((mod->initqrates == ECMUSR) && (mod->freq_model == FUNDEFINED || mod->freq_model == FMODEL)) {
                For(i, mod->ns) mod->pi[i] = io->mod->userfreq[senseCodons[i]];
                if(mod->whichrealmodel == MG) {
                    For(i, mod->num_base_freq) mod->base_freq[i] = io->mod->userbfreq[i];
                }
                mod->s_opt->user_state_freq = NO;
            }
            
            if(mod->whichrealmodel == PCM) {
                mod->s_opt->user_state_freq = NO;
            } else if(mod->initqrates == KOSI07) {
                For(i, mod->ns) mod->pi[i] = ecmK07freq[i];
                mod->s_opt->user_state_freq = NO;
            } else if(mod->initqrates == SCHN05) {
                For(i, mod->ns) mod->pi[i] = ecmS05freq[senseCodons[i]];
                mod->s_opt->user_state_freq = NO;
            }
            mod->s_opt->opt_state_freq = NO;
        }
        
        if((mod->s_opt->opt_omega == NO) && 
           (mod->s_opt->opt_kappa == NO) && 
           (mod->s_opt->opt_state_freq == NO)) {
        	if(mod->nparts > 1){printf("options not compatible with partitioned model error 8\n");exit(EXIT_FAILURE);}
            mr = Update_Qmat_Codons(mod, 0, 0); //modified by Ken 19/8
            EigenQREV(mod->qmat, mod->pi, mod->ns, mod->eigen->e_val, mod->eigen->r_e_vect, mod->eigen->l_e_vect, mod->eigen->space);
            For(i, io->mod->ns) mod->eigen->e_val[i] /= mr;
            mod->update_eigen = NO;
        } else {
            mod->update_eigen = YES;
            Set_Model_Parameters(mod);
            mod->update_eigen = NO;
        }
        
        //mod->omega_old  = mod->omega; //Ken 18/8
        mod->beta_old = mod->beta;
    }//!< Finish Added by Marcelo.
    else
    {
        if(io->datatype == GENERIC) mod->whichmodel = JC69;
        if(!mod->invar) For(i,data->crunch_len) data->invar[i] = 0;
        if(mod->s_opt->opt_alpha)   mod->alpha  = 1.0;
        if(mod->s_opt->opt_pinvar)  mod->pinvar = 0.2;
        
        For(i,mod->ns) 
        {
            mod->pi[i] = data->b_frq[i];
            mod->pi_unscaled[i] = mod->pi[i] * 100.;
        }
        
        if(io->datatype == NT)
        {
            /* Set the substitution parameters to their default values
             when they are not fixed by the user */
            if(mod->s_opt->opt_kappa) 
            {
                mod->kappa  = 4.0;
                mod->lambda = 1.0;
            }
            if(mod->s_opt->opt_rr)
            {
                int i;
                For(i,6) 
                {
                    mod->rr[i]     = 1.0;
                    mod->rr_val[i] = 1.0;
                }
            }
            mod->update_eigen = 1;
            mod->lambda       = 1.;
            if(mod->whichmodel == JC69)
            {
                mod->pi[0] = mod->pi[1] = mod->pi[2] = mod->pi[3] = .25;
                mod->kappa = 1.;
                mod->s_opt->opt_state_freq = NO;
                mod->s_opt->opt_kappa      = NO;
                mod->s_opt->opt_lambda     = NO;
                mod->update_eigen          = NO;
            }
            if(mod->whichmodel == K80)
            {
                mod->pi[0] = mod->pi[1] = mod->pi[2] = mod->pi[3] = .25;
                mod->s_opt->opt_state_freq = NO;
                mod->s_opt->opt_lambda     = NO;
                mod->update_eigen          = NO;
            }
            if(mod->whichmodel == F81)
            {
                mod->kappa                 = 1.;
                mod->s_opt->opt_kappa      = NO;
                mod->update_eigen          = NO;
            }
            if(mod->whichmodel == F84)
            {
                aux = ((mod->pi[0]+mod->pi[2])-(mod->pi[1]+mod->pi[3]))/(2.*mod->kappa);
                mod->lambda = ((mod->pi[1]+mod->pi[3]) + aux)/((mod->pi[0]+mod->pi[2]) - aux); 
                mod->update_eigen          = NO;
            }
            if(mod->whichmodel == TN93)
            {
                mod->update_eigen          = NO;
                if(io->mod->s_opt->opt_kappa) io->mod->s_opt->opt_lambda = 1;
            }
            if(mod->whichmodel == GTR)
            {
                mod->kappa = 1.;
                mod->update_eigen          = YES;
                io->mod->s_opt->opt_rr     = YES;
            }
            if(mod->whichmodel == CUSTOM)
            {
                mod->kappa = 1.;
                mod->update_eigen          = YES;
                /* 	  io->mod->s_opt->opt_rr     = YES; */ /* What if the user decided not to optimise the rates? */
            }
            if(mod->whichmodel == GTR)
            {		  
                mod->custom_mod_string[0] = '0';
                mod->custom_mod_string[1] = '1';
                mod->custom_mod_string[2] = '2';
                mod->custom_mod_string[3] = '3';
                mod->custom_mod_string[4] = '4';
                mod->custom_mod_string[5] = '5';
                Translate_Custom_Mod_String(mod);
            }
            if(mod->s_opt->user_state_freq) 
                For(i,4)  
            {
                mod->pi[i] = mod->user_b_freq[i];
            }
            if(!mod->use_m4mod) Set_Model_Parameters(mod);      
            if((mod->whichmodel != GTR)    && 
               (mod->whichmodel != CUSTOM) && 
               (mod->whichmodel != HKY85)) mod->update_eigen = 0;
        }
        else
        {
            PhyML_Printf("\n. Not implemented yet.\n");
            PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
            Warn_And_Exit("");
        }
    }
    mod->alpha_old   = mod->alpha;
    mod->kappa_old   = mod->kappa;
    mod->lambda_old  = mod->lambda;
    mod->pinvar_old  = mod->pinvar;
}

/*********************************************************/

void Update_Qmat_Generic(phydbl *rr, phydbl *pi, int ns, phydbl *qmat)
{
    int i,j;
    phydbl sum,mr;
    
    For(i,ns*ns) qmat[i] = .0;
    
    if(rr[(int)(ns*(ns-1)/2)-1] < 0.00001) {
        PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
        Warn_And_Exit("");
    }
    
    /*   PhyML_Printf("\n"); */
    /*   For(i,(int)(ns*(ns-1)/2))  */
    /*     { */
    /*       PhyML_Printf("\n> rr %d = %f",i,rr[i]); */
    /*     } */
    
    For(i,(int)(ns*(ns-1)/2)) 
    { 
        rr[i] /= rr[(int)(ns*(ns-1)/2)-1];
    }
    
    /* Fill the non-diagonal parts */
    For(i,ns)
    {
        for(j=i+1;j<ns;j++)
        {
            qmat[i*ns+j] = rr[MIN(i,j) * ns + MAX(i,j) -
                              (MIN(i,j)+1+(int)POW(MIN(i,j)+1,2))/2];
            qmat[j*ns+i] = qmat[i*ns+j];
        }
    }
    
    
    /* Multiply by pi */
    For(i,ns)
    {
        For(j,ns)
        {
            qmat[i*ns+j] *= pi[j];
        }
    }
    
    /* Compute diagonal elements */
    mr = .0;
    For(i,ns)
    {
        sum = .0;  
        For(j,ns) {sum += qmat[i*ns+j];}
        qmat[i*ns+i] = -sum;
        mr += sum * pi[i];
    }
    
    /*   For(i,ns) For(j,ns) qmat[i*ns+j] /= mr; */
}

/*********************************************************/

void Update_Qmat_GTR(phydbl *rr, phydbl *rr_val, int *rr_num, phydbl *pi, phydbl *qmat)
{
    int i;
    phydbl mr;
    
    For(i,6) rr[i] = rr_val[rr_num[i]];
    For(i,6) if(rr[i] < 0.001) rr[i] = 0.001;
    For(i,6) rr[i] /= rr[5];
    
    qmat[0*4+1] = (rr[0]*pi[1]);
    qmat[0*4+2] = (rr[1]*pi[2]);
    qmat[0*4+3] = (rr[2]*pi[3]);
    
    qmat[1*4+0] = (rr[0]*pi[0]);
    qmat[1*4+2] = (rr[3]*pi[2]);
    qmat[1*4+3] = (rr[4]*pi[3]);
    
    qmat[2*4+0] = (rr[1]*pi[0]);
    qmat[2*4+1] = (rr[3]*pi[1]);
    qmat[2*4+3] = (rr[5]*pi[3]);
    
    qmat[3*4+0] = (rr[2]*pi[0]);
    qmat[3*4+1] = (rr[4]*pi[1]);
    qmat[3*4+2] = (rr[5]*pi[2]);
    
    qmat[0*4+0] = -(rr[0]*pi[1]+rr[1]*pi[2]+rr[2]*pi[3]);
    qmat[1*4+1] = -(rr[0]*pi[0]+rr[3]*pi[2]+rr[4]*pi[3]);
    qmat[2*4+2] = -(rr[1]*pi[0]+rr[3]*pi[1]+rr[5]*pi[3]);
    qmat[3*4+3] = -(rr[2]*pi[0]+rr[4]*pi[1]+rr[5]*pi[2]);
    
    /* compute diagonal terms of Q and mean rate mr = l/t */
    mr = .0;
    For (i,4) mr += pi[i] * (-qmat[i*4+i]);
    For(i,16) qmat[i] /= mr;
}

/*********************************************************/

void Update_Qmat_HKY(phydbl kappa, phydbl *pi, phydbl *qmat)
{
    int i;
    phydbl mr;
    
    /* A -> C */ qmat[0*4+1] = (phydbl)(pi[1]);
    /* A -> G */ qmat[0*4+2] = (phydbl)(kappa*pi[2]);
    /* A -> T */ qmat[0*4+3] = (phydbl)(pi[3]);
    
    /* C -> A */ qmat[1*4+0] = (phydbl)(pi[0]);
    /* C -> G */ qmat[1*4+2] = (phydbl)(pi[2]);
    /* C -> T */ qmat[1*4+3] = (phydbl)(kappa*pi[3]);
    
    /* G -> A */ qmat[2*4+0] = (phydbl)(kappa*pi[0]);
    /* G -> C */ qmat[2*4+1] = (phydbl)(pi[1]);
    /* G -> T */ qmat[2*4+3] = (phydbl)(pi[3]);
    
    /* T -> A */ qmat[3*4+0] = (phydbl)(pi[0]);
    /* T -> C */ qmat[3*4+1] = (phydbl)(kappa*pi[1]);
    /* T -> G */ qmat[3*4+2] = (phydbl)(pi[2]);
    
    qmat[0*4+0] = (phydbl)(-(qmat[0*4+1]+qmat[0*4+2]+qmat[0*4+3]));
    qmat[1*4+1] = (phydbl)(-(qmat[1*4+0]+qmat[1*4+2]+qmat[1*4+3]));
    qmat[2*4+2] = (phydbl)(-(qmat[2*4+0]+qmat[2*4+1]+qmat[2*4+3]));
    qmat[3*4+3] = (phydbl)(-(qmat[3*4+0]+qmat[3*4+1]+qmat[3*4+2]));
    
    /* compute diagonal terms of Q and mean rate mr = l/t */
    mr = .0;
    For (i,4) mr += pi[i] * (-qmat[i*4+i]);
    For(i,16) qmat[i] /= mr;
}
/*********************************************************/



void Translate_Custom_Mod_String(model *mod)
{
    int i,j;
    
    For(i,6) mod->n_rr_per_cat[i] = 0;
    
    mod->n_diff_rr = 0;
    
    For(i,6)
    {
        For(j,i)
        {
            if(mod->custom_mod_string[i] == mod->custom_mod_string[j])
            {
                break;
            }
        }
        
        if(i == j)
        {
            mod->rr_num[i] = mod->n_diff_rr;
            mod->n_diff_rr++;
        }
        else
        {
            mod->rr_num[i] = mod->rr_num[j];
        }
        
        mod->n_rr_per_cat[mod->rr_num[j]]++;
    }
    
    /*   PhyML_Printf("\n"); */
    /*   For(i,6) PhyML_Printf("%d ",mod->rr_param_num[i]); */
    /*   For(i,mod->n_diff_rr_param) PhyML_Printf("\n. Class %d size %d",i+1,mod->n_rr_param_per_cat[i]); */
}

/*********************************************************/
void Set_Model_Parameters(model *mod) {
    if(mod->io->datatype == CODON) { //!Added by Marcelo.
        phydbl mr;
        int i, j, k, n, nn, n_termsTaylor;
        switch(mod->omegaSiteVar) {
            case DGAMMAK: { 
                DiscreteGamma(mod->prob_omegas, mod->omegas, mod->alpha, mod->beta, mod->n_w_catg, mod->gamma_median);
                For(i, mod->n_w_catg) if(mod->omegas[i] < MODELPAREPS) mod->omegas[i] = MODELPAREPS;
                break;
            }
            case DMODELK: { 
                Freq_to_UnsFreq(mod->prob_omegas, mod->prob_omegas_uns, mod->n_w_catg, 0);
                Scale_freqs_tol(mod->prob_omegas, mod->n_w_catg, MODELPAREPS, 0.99);
                break;
            }
            case DM0: {
                mod->omegas[0] = 1.0;
                mod->prob_omegas[0] = 1.0;
                DiscreteGamma(mod->gamma_r_proba, mod->gamma_rr, mod->alpha, mod->alpha, mod->n_catg, mod->gamma_median);
                break;
            }
            default:
                break;
        }
        
        if((mod->s_opt->opt_state_freq) &&
           (mod->s_opt->opt_omega == NO)) {
            switch(mod->freq_model) {
                case F1XSENSECODONS: {
                    Freq_to_UnsFreq(mod->pi, mod->pi_unscaled, mod->ns, 0);
                    break;
                }
                case F1X4: {
                    Freq_to_UnsFreq(mod->base_freq, mod->uns_base_freq, mod->num_base_freq, 0);
                    EqFrequencies(mod->freq_model, mod->pi, mod->base_freq, mod->ns);
                    break;
                }
                case F3X4:
                case CF3X4: {
                    Freq_to_UnsFreq(mod->base_freq,   mod->uns_base_freq,   4, 0);
                    Freq_to_UnsFreq(mod->base_freq+4, mod->uns_base_freq+4, 4, 0);
                    Freq_to_UnsFreq(mod->base_freq+8, mod->uns_base_freq+8, 4, 0);
                    EqFrequencies(mod->freq_model, mod->pi, mod->base_freq, mod->ns);
                    break;
                }
                default:
                    break;
            }
            
            if(mod->whichrealmodel == MG) {
                switch(mod->freq_model) {
                    case F1X4:{
                        Scale_freqs(mod->base_freq, mod->num_base_freq);
                        break;
                    }
                    case F3X4:
                    case CF3X4: {
                        Scale_freqs(mod->base_freq,   4);
                        Scale_freqs(mod->base_freq+4, 4);
                        Scale_freqs(mod->base_freq+8, 4);
                        break;
                    }   
                    default:
                        break;
                }
            }
            
            Scale_freqs(mod->pi, mod->ns);
        }



        if(mod->update_eigen) {

            if(mod->whichrealmodel == PCM) Update_Rate_Matrix_PCAModel(mod);

            
            if(mod->n_w_catg == 1) {



                if(mod->io->kappaECM==kap5) 
                {
                    Freq_to_UnsFreq(mod->pkappa,   mod->unspkappa,  mod->nkappa, 0);
                    Scale_freqs(mod->pkappa, mod->nkappa);      
                }
                int modeli;
               for(modeli=0;modeli<mod->nparts;modeli++){ //Ken 19/8
            	   //for(modeli=0;modeli<1;modeli++){ //Ken 19/8
                mod->mr_w[0] = Update_Qmat_Codons(mod, 0, modeli); //Ken 19/8


                if(mod->io->expm == EIGEN) {
                	if(mod->nparts > 1){printf("Options not supported with partitioned models error 13\n");exit(EXIT_FAILURE);}
                    EigenQREV(mod->qmat_part[0], mod->pi, mod->ns, mod->eigen->e_val, mod->eigen->r_e_vect, mod->eigen->l_e_vect, mod->eigen->space);
                    For(i,mod->ns) mod->eigen->e_val[i]/=mod->mr_w[0];
                } else if(mod->io->expm == TAYLOR) {
                	if(mod->nparts > 1){printf("Options not supported with partitioned models error 14\n");exit(EXIT_FAILURE);}

                    int nn=mod->ns*mod->ns, n=mod->ns, l;
                    
                    For(i,nn) mod->qmat[i]/=mod->mr_w[0];
                    
                    For(i,nn) mod->A2_part[modeli][0*nn+i]=mod->A0_part[modeli][i];
                    
                    For(i,nn) mod->A2_part[modeli][nn+i]=mod->qmat_part[modeli][i];
                    
                    for(l=2;l<mod->io->n_termsTaylor;l++) {  
#if defined BLAS || defined BLAS_OMP
                        
                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, mod->A2_part[modeli]+(nn*(l-1)), n, mod->A2_part[modeli]+nn,  n, 0.0, mod->A2_part[modeli]+nn*l, n);
                        
#else
                        
                        For(i,nn) mod->A2_part[modeli][nn*l+i]=0.0;
                        
                        For(i,n) For(j,n) For(k,n) mod->A2_part[modeli][nn*l+n*i+j] += mod->A2_part[modeli][nn*(l-1)+i*n+k] * mod->A2_part[modeli][nn+k*n+j];
                        
#endif
                    }
                } else if(mod->io->expm == SSPADE) {

                	//Modified by Ken 17/8/2016
                    n=mod->ns;
                    nn=n*n;
                    For(i,nn) mod->qmat_part[modeli][i]/=mod->mr_w[0];


                    
#if defined BLAS || defined BLAS_OMP


                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, mod->qmat_part[modeli], n, mod->qmat_part[modeli],  n, 0.0, mod->A2_part[modeli], n);
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, mod->A2_part[modeli],   n, mod->A2_part[modeli],    n, 0.0, mod->A4_part[modeli], n);
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, mod->A2_part[modeli],   n, mod->A4_part[modeli],    n, 0.0, mod->A6_part[modeli], n);
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, mod->A2_part[modeli],   n, mod->A6_part[modeli],    n, 0.0, mod->A8_part[modeli], n);
                    
                    cblas_dcopy(nn,mod->A2_part[modeli], 1, mod->Apowers_part[modeli]+  nn,1);
                    cblas_dcopy(nn,mod->A4_part[modeli], 1, mod->Apowers_part[modeli]+2*nn,1);
                    cblas_dcopy(nn,mod->A6_part[modeli], 1, mod->Apowers_part[modeli]+3*nn,1);
                    cblas_dcopy(nn,mod->A8_part[modeli], 1, mod->Apowers_part[modeli]+4*nn,1);
                    

#else
                    
                    For(i,nn) mod->A2_part[modeli][i]=0.0;
                    For(i,nn) mod->A4_part[modeli][i]=0.0;
                    For(i,nn) mod->A6_part[modeli][i]=0.0;
                    For(i,nn) mod->A8_part[modeli][i]=0.0;
                    
                    For(i,n) For(j,n) For(k,n) mod->A2_part[modeli][n*i+j] += mod->qmat_part[modeli][i*n+k] * mod->qmat_part[modeli][k*n+j];
                    For(i,n) For(j,n) For(k,n) mod->A4_part[modeli][n*i+j] += mod->A2_part[modeli][i*n+k]   * mod->A2_part[modeli][k*n+j];
                    For(i,n) For(j,n) For(k,n) mod->A6_part[modeli][n*i+j] += mod->A2_part[modeli][i*n+k]   * mod->A4_part[modeli][k*n+j];
                    For(i,n) For(j,n) For(k,n) mod->A8_part[modeli][n*i+j] += mod->A2_part[modeli][i*n+k]   * mod->A6_part[modeli][k*n+j];
                    
                    For(i,nn) mod->Apowers_part[modeli][nn+i]  =mod->A2_part[modeli][i];
                    For(i,nn) mod->Apowers_part[modeli][2*nn+i]=mod->A4_part[modeli][i];
                    For(i,nn) mod->Apowers_part[modeli][3*nn+i]=mod->A6_part[modeli][i];
                    For(i,nn) mod->Apowers_part[modeli][4*nn+i]=mod->A8_part[modeli][i];
                    
#endif
                }
            }//for(modeli)
            } else if(mod->n_w_catg > 1) { 
            	if(mod->nparts > 1){
            		printf("CAN'T COMBINE PARTITIONS WITH SITE CLASSES\n");
            		exit(EXIT_FAILURE);
            	}
#if defined OMP || defined BLAS_OMP 
                
#pragma omp parallel for
                
#endif
                
                For(i, mod->n_w_catg) {
                    if(mod->io->kappaECM == kap5) {
                        Freq_to_UnsFreq(mod->pkappa,   mod->unspkappa,  mod->nkappa, 0);
                        Scale_freqs(mod->pkappa, mod->nkappa);    
                    }
                    mod->mr_w[i] = Update_Qmat_Codons(mod, i,0);
                }
                
                //!<Scaling Qmat together
                mr = 0.0;
                For(k, mod->n_w_catg) mr +=  mod->mr_w[k]*mod->prob_omegas[k];
                
                
                
                if(mod->io->expm==EIGEN) {
#if defined OMP || defined BLAS_OMP 
                    
#pragma omp parallel for 
                    
#endif
                    For(k, mod->n_w_catg) {
                        EigenQREV(mod->qmat + k*mod->ns*mod->ns,                                     mod->pi, mod->ns, mod->eigen->e_val + k*mod->ns,                                   mod->eigen->r_e_vect + k*mod->ns*mod->ns,                                              mod->eigen->l_e_vect + k*mod->ns*mod->ns, mod->eigen->space + k*2*mod->ns);
                    }
                    For(i, mod->ns * mod->n_w_catg) mod->eigen->e_val[i]/=mr;
                } else if(mod->io->expm == TAYLOR) {
                    int catg, nn=mod->ns*mod->ns, n=mod->ns,l;
                    For(i,nn*mod->n_w_catg) mod->qmat[i]/=mr;
                    For(i,nn) mod->A2_part[0][nn+i]=mod->qmat[i]; //modified by Ken 22/8- doesn't work with partitioned models!
                    
#if defined OMP || defined BLAS_OMP 
                    
#pragma omp parallel for private(i,j,k,l)
                    
#endif
                    
                    For(catg, mod->n_w_catg) {
                        for(l = 2; l < mod->io->n_termsTaylor; l++) {  
#if defined BLAS || defined BLAS_OMP
                        	if(mod->nparts > 1){printf("Options not supported with partitioned models error 15\n");exit(EXIT_FAILURE);}
                            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, mod->A2_part[0]+(nn*(l-1))+catg*nn*15, n, mod->qmat_part[0]+catg*nn,  n, 0.0, mod->A2_part[0]+nn*l+catg*nn*15, n);
                            
#else
                            
                            For(i,nn) mod->A2_part[0][nn*l+i]=0.0; //modified by Ken 22/8
                            For(i,n) For(j,n) For(k,n) mod->A2_part[0][nn*l+n*i+j] += mod->A2_part[0][nn*(l-1)+i*n+k] * mod->qmat_part[0][k*n+j+catg*nn];
                            
#endif
                        }
                    }
                } else if(mod->io->expm == SSPADE) {
                    n=mod->ns;
                    nn=n*n;
                    For(i,nn*mod->n_w_catg) mod->qmat_part[0][i]/=mr; //!< Scaling of the Q matrix.
                    
                    int catg;
                    
#if defined OMP || defined BLAS_OMP 
                    
#pragma omp parallel for private(i,j,k)
                    
#endif
                    
                    For(catg, mod->n_w_catg) {
#if defined BLAS || defined BLAS_OMP
                    	if(mod->nparts > 1){printf("Options not supported with partitioned models error 12\n");exit(EXIT_FAILURE);}
                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, mod->qmat_part[0]+catg*nn,       n, mod->qmat_part[0]+catg*nn,       n, 0.0, mod->A2_part[0]+catg*nn, n);
                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, mod->A2_part[0]+catg*nn,         n, mod->A2_part[0]+catg*nn,         n, 0.0, mod->A4_part[0]+catg*nn, n);
                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, mod->A2_part[0]+catg*nn,         n, mod->A4_part[0]+catg*nn,         n, 0.0, mod->A6_part[0]+catg*nn, n);
                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, mod->A2_part[0]+catg*nn,         n, mod->A4_part[0]+catg*nn,         n, 0.0, mod->A8_part[0]+catg*nn, n);
                        
                        cblas_dcopy(nn,mod->A2_part[0]+catg*nn, 1, mod->Apowers_part[0]+5*nn*catg+  nn,1);
                        cblas_dcopy(nn,mod->A4_part[0]+catg*nn, 1, mod->Apowers_part[0]+5*nn*catg+2*nn,1);
                        cblas_dcopy(nn,mod->A6_part[0]+catg*nn, 1, mod->Apowers_part[0]+5*nn*catg+3*nn,1);
                        cblas_dcopy(nn,mod->A8_part[0]+catg*nn, 1, mod->Apowers_part[0]+5*nn*catg+4*nn,1);
                        
#else
                    	if(mod->nparts > 1){printf("Options not supported with partitioned models error 12\n");exit(EXIT_FAILURE);}
                        For(i,nn) mod->A2_part[0][i+catg*nn]=0.0;//modified by Ken 22/8 won't work with partitioned models
                        For(i,nn) mod->A4_part[0][i+catg*nn]=0.0;
                        For(i,nn) mod->A6_part[0][i+catg*nn]=0.0;
                        For(i,nn) mod->A8_part[0][i+catg*nn]=0.0;
                        
                        For(i,n) For(j,n) For(k,n) mod->A2_part[0][n*i+j+catg*nn] += mod->qmat[i*n+k+catg*nn] * mod->qmat[k*n+j+catg*nn];
                        For(i,n) For(j,n) For(k,n) mod->A4_part[0][n*i+j+catg*nn] += mod->A2_part[0][i*n+k+catg*nn]   * mod->A2_part[0][k*n+j+catg*nn];
                        For(i,n) For(j,n) For(k,n) mod->A6_part[0][n*i+j+catg*nn] += mod->A2_part[0][i*n+k+catg*nn]   * mod->A4_part[0][k*n+j+catg*nn];
                        For(i,n) For(j,n) For(k,n) mod->A8_part[0][n*i+j+catg*nn] += mod->A2_part[0][i*n+k+catg*nn]   * mod->A6_part[0][k*n+j+catg*nn];
                        
                        For(i,nn) mod->Apowers_part[0][5*nn*catg+  nn+i]=mod->A2_part[0][i+catg*nn];
                        For(i,nn) mod->Apowers_part[0][5*nn*catg+2*nn+i]=mod->A4_part[0][i+catg*nn];
                        For(i,nn) mod->Apowers_part[0][5*nn*catg+3*nn+i]=mod->A6_part[0][i+catg*nn];
                        For(i,nn) mod->Apowers_part[0][5*nn*catg+4*nn+i]=mod->A8_part[0][i+catg*nn];
                        
#endif
                    }
                }
            }
        }

    } else {
        phydbl sum;
        int i;
        int result, n_iter;
        phydbl scalar;
        
        DiscreteGamma(mod->gamma_r_proba, mod->gamma_rr, mod->alpha, mod->alpha, mod->n_catg, mod->gamma_median);
        
        if(mod->n_rr_branch > 0)
            DiscreteGamma(mod->p_rr_branch, mod->rr_branch, mod->rr_branch_alpha, mod->rr_branch_alpha, mod->n_rr_branch, mod->gamma_median);
        
        if((mod->io->datatype == NT) && (mod->s_opt->opt_state_freq))
        {
            sum = .0;
            For(i,mod->ns) sum += FABS(mod->pi_unscaled[i]);
            For(i,mod->ns) mod->pi[i] = FABS(mod->pi_unscaled[i])/sum;
            
            do
            {
                sum = .0;
                For(i,mod->ns)
                {
                    if(mod->pi[i] < 0.01) mod->pi[i]=0.01;
                    if(mod->pi[i] > 0.99) mod->pi[i]=0.99;
                    sum += mod->pi[i];
                }
                For(i,mod->ns) mod->pi[i]/=sum;
            }
            while((sum > 1.01) || (sum < 0.99));
        }
        
        if(mod->update_eigen) 
        {
            int j;

            if (mod->io->datatype==AA) //!Added by Marcelo
            {
                For(i,mod->ns*mod->ns) mod->qmat[i] = .0;
                For(i,mod->ns       ) mod->pi[i]   = .0;
                switch(mod->whichmodel)
                {

                    default : break;
                }
                
                Freq_to_UnsFreq(mod->pi, mod->pi_unscaled, mod->ns, 0);
                Scale_freqs(mod->pi, mod->ns);
                
                /*       /\* multiply the nth col of Q by the nth term of pi/100 just as in PAML *\/ */
                For(i,mod->ns) For(j,mod->ns) mod->qmat[i*mod->ns+j] *= mod->pi[j]/100.0;
                /* compute diagonal terms of Q and mean rate mr = l/t */
                mod->mr = .0;
                For (i,mod->ns)
                {
                    sum=.0;
                    For(j, mod->ns) sum += mod->qmat[i*mod->ns+j];
                    mod->qmat[i*mod->ns+i] = -sum;
                    mod->mr += mod->pi[i] * sum;
                }
                /* scale imod->nstantaneous rate matrix so that mu=1 */
                For (i,mod->ns*mod->ns) mod->qmat[i] /= mod->mr;
                /* compute eigenvectors/values */
                result = 0;
            	if(mod->nparts > 1){printf("options not compatible with partitioned model error 10\n");exit(EXIT_FAILURE);}
                For(i,mod->ns*mod->ns) mod->qmat_buff_part[0][i] = mod->qmat_part[0][i];
                if(!Eigen(1,mod->qmat_buff_part[0],mod->eigen->size,mod->eigen->e_val,
                          mod->eigen->e_val_im,mod->eigen->r_e_vect,
                          mod->eigen->r_e_vect_im,mod->eigen->space))
                {
                    /* compute inverse(Vr) into Vi */
                    For (i,mod->ns*mod->ns) mod->eigen->l_e_vect[i] = mod->eigen->r_e_vect[i];
                    if(!Matinv(mod->eigen->l_e_vect,mod->eigen->size,mod->eigen->size))
                    {
                        PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
                        Exit("\n");      
                    }
                    /* compute the diagonal terms of EXP(D) */
                    For(i,mod->ns) mod->eigen->e_val[i] = (phydbl)EXP(mod->eigen->e_val[i]);
                }
                else
                {
                    if (result==-1) PhyML_Printf("\n. Eigenvalues/vectors computation does not converge : computation cancelled");
                    else if (result==1) PhyML_Printf("\n. Complex eigenvalues/vectors : computation cancelled");
                }
            }
            else
            {
                if(!mod->use_m4mod)
                {
                    if(mod->io->datatype == NT)
                    {
                        if(mod->whichmodel == GTR)
                            Update_Qmat_GTR(mod->rr, mod->rr_val, mod->rr_num, mod->pi, mod->qmat);
                        else if(mod->whichmodel == CUSTOM) 
                            Update_Qmat_GTR(mod->rr, mod->rr_val, mod->rr_num, mod->pi, mod->qmat);
                        else if(mod->whichmodel == HKY85)  
                            Update_Qmat_HKY(mod->kappa, mod->pi, mod->qmat);
                        else /* Any other nucleotide-based model */
                            Update_Qmat_HKY(mod->kappa, mod->pi, mod->qmat);
                    }
                    
                }
#ifdef M4
                else 
                {
                    M4_Update_Qmat(mod->m4mod,mod);
                }
#endif
                
                scalar   = 1.0;
                n_iter   = 0;
                result   = 0;
            	if(mod->nparts > 1){printf("options not compatible with partitioned model error 11\n");exit(EXIT_FAILURE);}
                For(i,mod->ns*mod->ns) mod->qmat_buff_part[0][i] = mod->qmat_part[0][i];
                
                /* compute eigenvectors/values */
                /*       if(!EigenRealGeneral(mod->eigen->size,mod->qmat,mod->eigen->e_val, */
                /* 			  mod->eigen->e_val_im, mod->eigen->r_e_vect, */
                /* 			  mod->eigen->space_int,mod->eigen->space)) */
                
                if(!Eigen(1,mod->qmat_buff_part[0],mod->eigen->size,mod->eigen->e_val,mod->eigen->e_val_im,mod->eigen->r_e_vect,mod->eigen->r_e_vect_im,mod->eigen->space))
                {
                    /* compute inverse(Vr) into Vi */
                    For (i,mod->ns*mod->ns) mod->eigen->l_e_vect[i] = mod->eigen->r_e_vect[i];
                    while(!Matinv(mod->eigen->l_e_vect, mod->eigen->size, mod->eigen->size))
                    {
                        PhyML_Printf("\n. Trying Q<-Q*scalar and then Root<-Root/scalar to fix this...\n");
                        scalar += scalar / 3.;
                        For(i,mod->eigen->size*mod->eigen->size) mod->qmat_buff_part[0][i]  = mod->qmat_part[0][i];
                        For(i,mod->eigen->size*mod->eigen->size) mod->qmat_buff_part[0][i] *= scalar;
                        result = Eigen(1,mod->qmat_buff_part[0],mod->eigen->size,mod->eigen->e_val,mod->eigen->e_val_im,mod->eigen->r_e_vect,mod->eigen->r_e_vect_im,mod->eigen->space);
                        if (result == -1)
                            Exit("\n. Eigenvalues/vectors computation did not converge : computation cancelled\n");
                        else if (result == 1)
                            Exit("\n. Complex eigenvalues/vectors : computation cancelled\n");
                        
                        For (i,mod->eigen->size*mod->eigen->size) mod->eigen->l_e_vect[i] = mod->eigen->r_e_vect[i];
                        n_iter++;
                        if(n_iter > 100) Exit("\n. Cannot work out eigen vectors\n");
                    };
                    For(i,mod->eigen->size) mod->eigen->e_val[i] /= scalar;
                    
                    /* compute the diagonal terms of EXP(D) */
                    For(i,mod->ns) mod->eigen->e_val[i] = (phydbl)EXP(mod->eigen->e_val[i]);
                    
                }
                else
                {
                    PhyML_Printf("\n. Eigenvalues/vectors computation does not converge : computation cancelled");
                    Warn_And_Exit("\n");
                }
            }
        }
    }
}

/*********************************************************/

void Switch_From_M4mod_To_Mod(model *mod)
{
    int i;
    
    mod->use_m4mod = 0;
    mod->ns = mod->m4mod->n_o;
    For(i,mod->ns) mod->pi[i] = mod->m4mod->o_fq[i];
    mod->eigen->size = mod->ns;
    mod->update_eigen = 1;
}

/*********************************************************/

void Switch_From_Mod_To_M4mod(model *mod)
{
    int i;
    mod->use_m4mod = 1;
    mod->ns = mod->m4mod->n_o * mod->m4mod->n_h;
    For(i,mod->ns) mod->pi[i] = mod->m4mod->o_fq[i%mod->m4mod->n_o] * mod->m4mod->h_fq[i/mod->m4mod->n_o];
    mod->eigen->size = mod->ns;
    mod->update_eigen = 1;
}


phydbl F1x4(int codon, phydbl *freq) //!<Added by marcelo.
{
    phydbl val=1.0;
    int i, remainder=0;
    For(i,3){
        remainder=codon-((codon>>2)<<2); //!<  n<<k= n*2^k, n>>k=n/2^k.
        codon=codon>>2;
        val*=freq[remainder];
    }
    return val;
}

/*********************************************************/

phydbl F3x4(int codon, phydbl *freq) { //!<Added by marcelo.
    phydbl val=1.0;
    int i, remainder=0;
    for(i=2;i>=0;i--){
        remainder=codon-((codon>>2)<<2); //!<  n<<k= n*2^k, n>>k=n/2^k.
        codon=codon>>2;
        val*=freq[ remainder+i*(2<<1) ];
    }
    return val;
}

/*********************************************************/

void EqFrequencies(int modfreq, phydbl *pi, phydbl *freq, int numSensecodons) //!<Added by marcelo.
{ 
    int i;
    phydbl freqStopcodons;
    
    switch(modfreq) {
        case F1XSENSECODONS: {
            For(i,numSensecodons) pi[i]=freq[i];                     //!< Initialize the values before optimization 
            break;
        }
        case F1X4: {
            freqStopcodons=0.0;                                       //! Calculate the total frequency of stop codons to correct the frequency of the sense codons.
            For(i,64) if(stopCodons[i]) freqStopcodons+=F1x4(i,freq);
            For(i,numSensecodons) pi[i] = F1x4(senseCodons[i], freq)/(1-freqStopcodons); 
            break;
        }
        case F3X4:
        case CF3X4: {
            freqStopcodons=0.0;                                       //! Calculate the total frequency of stop codons to correct the frequency of the sense codons.
            For(i,64) if(stopCodons[i]) freqStopcodons+=F3x4(i,freq);
            For(i,numSensecodons) pi[i] = F3x4(senseCodons[i], freq)/(1-freqStopcodons); 
            break;
        }
        default :{
            Warn_And_Exit("Frequency model not implemented.\n"); break;
        }
    }
}

phydbl Update_Qmat_Codons(model *mod, int cat, int modeli) {
    int numSensecodons, i, j, allocFreqs;
    phydbl sum, mu;
    phydbl *freqs, *qmat, *mat;
    
    //Added by Ken
    //incorporate equilibrium frequencies in case they have been updated
    switch(mod->freq_model) {
        case F1XSENSECODONS: {
          	Freq_to_UnsFreq(mod->pi, mod->pi_unscaled, mod->ns, 0);
          	break;
        }
        case F1X4: {
            Freq_to_UnsFreq(mod->base_freq, mod->uns_base_freq, mod->num_base_freq, 0);
            EqFrequencies(mod->freq_model, mod->pi, mod->base_freq, mod->ns);
            break;
        }
        case F3X4:
        case CF3X4: {
            Freq_to_UnsFreq(mod->base_freq,   mod->uns_base_freq,   4, 0);
            Freq_to_UnsFreq(mod->base_freq+4, mod->uns_base_freq+4, 4, 0);
            Freq_to_UnsFreq(mod->base_freq+8, mod->uns_base_freq+8, 4, 0);
            EqFrequencies(mod->freq_model, mod->pi, mod->base_freq, mod->ns);
            break;
        }
        default:
            break;
     }

    //Added by Ken 17/8/2016
   // int modeli;
   // for(modeli=0;modeli<mod->nomega_parts;modeli++){

    qmat = mod->qmat_part[modeli] + (cat * mod->ns * mod->ns);
    numSensecodons = mod->ns;
    allocFreqs = NO;
    
    // Deal with initial rate matrices
    if(mod->initqrates != NOINITMAT) {
        switch(mod->initqrates) {
            case KOSI07:
                mat = (phydbl *) ecmK07;
                printf("USING KOSI07\n");
                freqs = ecmK07freq;
                break;
                
            case SCHN05:
                mat = (phydbl *) mCalloc(numSensecodons * numSensecodons, sizeof(phydbl));
                printf("USING SHN05n");
                For(i, numSensecodons) {
                    For(j, numSensecodons) {
                        mat[ i*numSensecodons + j ] = ecmS05[senseCodons[i]][senseCodons[j]];
                    }
                }
                freqs = (phydbl *) mCalloc(numSensecodons, sizeof(phydbl));
                For(i, numSensecodons) {
                    freqs[i] = ecmS05freq[senseCodons[i]];
                }
                allocFreqs = YES;
                break;
                
            case ECMUSR:
                if(mod->calculate_init_tree == 1) {
                    mat = (phydbl *) mCalloc(numSensecodons * numSensecodons, sizeof(phydbl));
                    For(i, numSensecodons) {
                        For(j, numSensecodons) mat[i*numSensecodons+j] = mod->userRatesT[senseCodons[i]][senseCodons[j]];
                    }
                    freqs = (phydbl *) mCalloc(numSensecodons, sizeof(phydbl));
                    For(i, numSensecodons) freqs[i] = mod->userfreqT[senseCodons[i]];
                    allocFreqs = YES;
                } else {
                    mat = (phydbl *) mCalloc(numSensecodons * numSensecodons, sizeof(phydbl));
                    For(i, numSensecodons) {
                        For(j, numSensecodons) mat[i*numSensecodons+j] = mod->userRates[senseCodons[i]][senseCodons[j]];
                    }
                    freqs = (phydbl *) mCalloc(numSensecodons, sizeof(phydbl));
                    For(i, numSensecodons) freqs[i] = mod->userfreq[senseCodons[i]];
                    allocFreqs = YES;
                }
                break;
                
            default:
                break;
        }
    } else if(mod->pcaModel) {
        mat = (phydbl *) mCalloc(numSensecodons * numSensecodons, sizeof(phydbl));
        For(i, numSensecodons) {
            For(j, numSensecodons) mat[i*numSensecodons+j] = mod->userRates[senseCodons[i]][senseCodons[j]];
        }
    } else {
        mat = (phydbl *) mCalloc(numSensecodons * numSensecodons, sizeof(phydbl));
        For(i, numSensecodons * numSensecodons) mat[i] = 1.0;
        
        freqs = (phydbl *) mCalloc(numSensecodons, sizeof(phydbl));
        For(i, numSensecodons) freqs[i] = 1.0;
    }
    // deal with frequency models
    if((mod->freq_model != FMODEL && mod->freq_model != FUNDEFINED)) {
        if((mod->initqrates == NOINITMAT && mod->pcaModel == NO) || mod->initqrates == SCHN05) {
            free(freqs);
            allocFreqs = NO;
        }
        freqs = mod->pi;
    } else {
        for(i=0; i<numSensecodons; i++) {
            mod->pi[i] = freqs[i];
        }
    }

    For(i, numSensecodons*numSensecodons) qmat[i] = 0.0;
    
    // calculate the actual Q matrix
    switch(mod->whichrealmodel) {
    	case GY: case PCM:
            Update_Qmat_GY(mat, qmat, freqs, cat, mod);
            break;
        case MG:
            Update_Qmat_MG(mat, qmat, mod->base_freq, cat, mod);
            break;
        case YAP:
            Update_Qmat_YAP(mat, qmat, freqs, cat, mod);
            break;
        case HLP17:
            Update_Qmat_HLP17(mat, qmat, freqs, cat, mod,mod->omega_part[modeli]);
            break;
        default:
            break;
    }


    /*! Calculate the diagonal element which is the negative of the sum of the other elements in a row.*/
    For(i, numSensecodons) {
        sum = 0.0;
        For(j, numSensecodons) sum += qmat[ numSensecodons*i+j ];
        qmat[ numSensecodons*i+i ]= (-1)*sum;    
    }
    
    /*! Normalize Q matrix.*/
    mu = 0;
    For(i, numSensecodons) {
        mu += freqs[i] * (-qmat[numSensecodons*i+i]);
    }

    //Added by Ken - normalizes matrix here.
    if(mod->whichrealmodel == HLP17){
      For(i, numSensecodons) {
        For(j, numSensecodons) {
    		qmat[numSensecodons*i+j] = qmat[numSensecodons*i+j]/mu;
    	}
       }
       mu = 0;
       For(i, numSensecodons) {
           mu += freqs[i] * (-qmat[numSensecodons*i+i]);
       }
    }

    if(mod->initqrates == NOINITMAT) {
       free(mat);
    }
    
    if(mod->initqrates == SCHN05) {
        if(allocFreqs) {
            free(freqs);
        }
        free(mat);
    }
   // }
   // mu=1.0;
    return mu;
}

void Update_Qmat_GY(phydbl *mat, phydbl *qmat, phydbl * freqs, int cat, model *mod) {
    int i, j, numSensecodons;
    phydbl value;
    numSensecodons = mod->ns;
    For(i, numSensecodons) {
        For(j, i) {
            value = mat[i*numSensecodons+j] * Kappa_Omega_Factor(i, j, mod, cat,mod->omega_part[0]);
              qmat[ i*numSensecodons+j ] = value * (freqs[j]);
              qmat[ j*numSensecodons+i ] = value * (freqs[i]);
        }
    }
}

//Modified GY94 by Ken
//Modified to be omega*kappa*pi*(1+b*h) on 12/Jun/2016
void Update_Qmat_HLP17(phydbl *mat, phydbl *qmat, phydbl * freqs, int cat, model *mod,phydbl omega) {
	 int fi,ti,li,ri,hot,c;
	 double htotal[mod->nmotifs];
	 for(fi=0;fi<61;fi++){ //Fill in B matrix
    	for(ti=0;ti<61;ti++){
    		for(c=0;c< mod->nmotifs;c++){ //set htotal array to zero
 	    			htotal[c]=0;
    		}
    		for(li=0;li<61;li++){
    			for(ri=0;ri<61;ri++){
    				for(c=0;c< mod->nmotifs;c++){ //tally up expected number of each type of hotspot mutation
    						htotal[c] += freqs[li]*freqs[ri]*mod->hotspotcmps[c][fi*61*61*61+ti*61*61+li*61+ri];
	 	    		}
    			}
	 	    }
            //additive interaction function
            double hot = 0;
           		for(c=0;c<mod->nmotifs;c++){
               		hot += htotal[c]*mod->hotness[mod->motif_hotness[c]];
           		}
            //constrain total modification to never go below -1
            if(hot < -1){
              	hot=-1;
            }
    		mod->Bmat[fi*61+ti]=hot;
    	}
    }

    int i, j, numSensecodons;
    phydbl value;
    numSensecodons = mod->ns;
    For(i, numSensecodons) {
        For(j, i) {
            	value = mat[i*numSensecodons+j] * Kappa_Omega_Factor(i, j, mod, cat,omega);
              qmat[ i*numSensecodons+j ] = value * freqs[j]*(1+mod->Bmat[ i*numSensecodons+j ]);
              qmat[ j*numSensecodons+i ] = value * freqs[i]*(1+mod->Bmat[ j*numSensecodons+i ]);
        }
    }
}


void Update_Qmat_MG(phydbl *mat, phydbl *qmat, phydbl * freqs, int cat, model *mod) {
    int diff, codoni, codonj, numSensecodons, i, j, k, m;
    int icodon[3], jcodon[3], targetNTi[3], targetNTj[3], posTargetNT[3];
    phydbl value, freqTargetNTi[3], freqTargetNTj[3], freqTNTi, freqTNTj;
    
    numSensecodons = mod->ns;
    For(i, numSensecodons) {
        For(j, i) {
            diff = 0;
            codoni = senseCodons[i];
            codonj = senseCodons[j];
            
            For(k, 3) {
                icodon[k] = codoni - ((codoni >> 2) << 2);                     //!<  n<<k= n*2^k, n>>k=n/2^k.
                codoni = codoni >> 2;
                jcodon[k] = codonj - ((codonj >> 2) << 2);
                codonj = codonj >> 2;
                if(icodon[k] != jcodon[k]) {
                    diff++;
                    targetNTi[diff-1] = icodon[k];                         //!< Store the target nucleotides.
                    targetNTj[diff-1] = jcodon[k];                         //!< Store the target nucleotides.
                    posTargetNT[diff-1] = 2-k;                             //!< Store the codon position of the target nucleotide.  
                }
            }
            
            if(diff > 0) {
                //!< Set the target nucleotide frequency according to the given models.
                switch(mod->freq_model) {
                    case F1X4: {
                        For(m, diff) freqTargetNTi[m] = freqs[targetNTi[m]];
                        For(m, diff) freqTargetNTj[m] = freqs[targetNTj[m]];
                        break;
                    }
                    case F3X4:
                    case CF3X4: {
                        For(m, diff) freqTargetNTi[m] = freqs[ 4*posTargetNT[m] + targetNTi[m] ];
                        For(m, diff) freqTargetNTj[m] = freqs[ 4*posTargetNT[m] + targetNTj[m] ];
                        break;
                    }
                    default: Warn_And_Exit("Frequency model not available in MG framework."); break;
                }
                
                value = mat[i*numSensecodons+j] * Kappa_Omega_Factor(i, j, mod, cat,mod->omega_part[0]);
                
                freqTNTj = 1.0;
                freqTNTi = 1.0;
                
                For(m, diff) freqTNTi *= freqTargetNTi[m];
                For(m, diff) freqTNTj *= freqTargetNTj[m];
                
                if(freqTNTi < MODELPAREPS) freqTNTi = MODELPAREPS;
                if(freqTNTj < MODELPAREPS) freqTNTj = MODELPAREPS;
                
                qmat[numSensecodons*i+j] = value * freqTNTj;
                qmat[numSensecodons*j+i] = value * freqTNTi;
            }
        }
    }
}

void Update_Qmat_YAP(phydbl *mat, phydbl *qmat, phydbl * freqs, int cat, model *mod) {
    int diff, codoni, codonj, numSensecodons, i, j, k, m, n, TC, h, g;
    int icodon[3], jcodon[3], targetCD[3], posTargetNT[3];
    phydbl margFreq, condTCi, condTCj, value;
    
    numSensecodons = mod->ns;
    
    For(i, numSensecodons) {
        For(j, i) {
            diff = 0;
            codoni = senseCodons[i];
            codonj = senseCodons[j];
            For(k, 3) {
                icodon[k] = codoni - ((codoni >> 2) << 2); //!<  n<<k= n*2^k, n>>k=n/2^k.
                codoni = codoni >> 2;
                jcodon[k] = codonj - ((codonj >> 2) << 2);
                codonj = codonj >> 2;
                
                targetCD[2-k] = jcodon[k]; //!< Transcribe the codon into the constituent nucleotides-
                
                if(icodon[k] != jcodon[k]) {
                    diff++;
                    posTargetNT[diff-1] = 2-k; //!< Store the codon position of the target nucleotide.  
                }
            }
            
            if(diff > 0) {
                switch(diff) {
                    case 1: {
                        margFreq = 0.0;
                        switch(posTargetNT[0]) {
                            case 0: {
                                m = 1;
                                n = 2;
                                For(h, 4) {
                                    TC = 16*h + 4*targetCD[m] + targetCD[n];
                                    if(!stopCodons[TC]) margFreq += freqs[indexSenseCodons[TC]]; 
                                }
                                break;
                            }
                            case 1: {
                                m = 0;
                                n = 2;
                                For(h, 4) {
                                    TC = 16*targetCD[m] + 4*h + targetCD[n];
                                    if(!stopCodons[TC]) margFreq += freqs[indexSenseCodons[TC]]; 
                                }
                                break;
                            }
                            case 2: {
                                m = 0;
                                n = 1;
                                For(h, 4) {
                                    TC = 16*targetCD[m] + 4*targetCD[n] + h;
                                    if(!stopCodons[TC]) margFreq += freqs[indexSenseCodons[TC]]; 
                                }
                                break;
                            }
                            default: break;
                        }
                        break;
                    }
                    case 2: {
                        margFreq = 0.0;
                        switch(posTargetNT[0]) {
                            case 0: {
                                if(posTargetNT[1] == 1) {
                                    n = 2;
                                    For(h, 4) {
                                        For(g, 4) {
                                            TC = 16*h + 4*g + targetCD[n];
                                            if(!stopCodons[TC]) margFreq += freqs[indexSenseCodons[TC]];
                                        } 
                                    }
                                } else if(posTargetNT[1]==2) {
                                    n = 1;
                                    For(h, 4) {
                                        For(g, 4) {
                                            TC = 16*h + 4*targetCD[n] + g;
                                            if(!stopCodons[TC]) margFreq += freqs[indexSenseCodons[TC]];
                                        } 
                                    }
                                }
                                break;	    
                            }
                            case 1: {
                                if(posTargetNT[1] == 0) {
                                    n = 2;
                                    For(h, 4) {
                                        For(g, 4) {
                                            TC = 16*g + 4*h + targetCD[n];
                                            if(!stopCodons[TC]) margFreq += freqs[indexSenseCodons[TC]];
                                        } 
                                    }
                                } else if(posTargetNT[1]==2) {
                                    n = 0;
                                    For(h, 4) {
                                        For(g, 4) {
                                            TC = 16*targetCD[n] + 4*h + g;
                                            if(!stopCodons[TC]) margFreq += freqs[indexSenseCodons[TC]];
                                        } 
                                    }
                                }
                                break;	    
                            }
                            case 2: {
                                if(posTargetNT[1]  ==  0) {
                                    n = 1;
                                    For(h, 4) {
                                        For(g, 4) {
                                            TC = 16*g + 4*targetCD[n] + h;
                                            if(!stopCodons[TC]) margFreq += freqs[indexSenseCodons[TC]];
                                        } 
                                    }
                                } else if(posTargetNT[1] == 1) {
                                    n = 0;
                                    For(h, 4) {
                                        For(g, 4) {
                                            TC = 16*targetCD[n] + 4*g + h;
                                            if(!stopCodons[TC]) margFreq += freqs[indexSenseCodons[TC]];
                                        } 
                                    }
                                }
                                break;	    
                            }
                            default:
                                break;
                        }
                        break;
                    }
                    case 3: {
                        margFreq = 1.0;
                        break;
                    }
                    default:
                        break;
                }
                
                condTCi = freqs[i] / margFreq;
                condTCj = freqs[j] / margFreq;
                
                if(condTCi < MODELPAREPS) condTCi = MODELPAREPS;
                if(condTCj < MODELPAREPS) condTCj = MODELPAREPS;
                
                value = mat[i*numSensecodons+j] * Kappa_Omega_Factor(i, j, mod, cat,mod->omega_part[0]);
                
                qmat[numSensecodons*i+j] = value * condTCj;
                qmat[numSensecodons*j+i] = value * condTCi;
            }
        }
    }
}


/*******************************************************/

// Converts ECM omega values to 'normal' ones as given in Kosiol 2007
phydbl Omega_ECMtoMmechModels(phydbl *pi, phydbl *qmat, phydbl *qmatbuff, int numSensecodons, int n_w_catg) //!< Added by Marcelo.
{	
    phydbl val, roa,ros;
    int i,j;
    
    roa = 0.0;  // nonsynonymous
    ros = 0.0;  // synonymous
    
    For(i, numSensecodons*numSensecodons) qmatbuff[i] = qmat[i];
    val = 0.0;
    For(i, numSensecodons) val += pi[i] * (-qmatbuff[i*numSensecodons+i]);
    For(i, numSensecodons * numSensecodons) qmatbuff[i] /= val;
    
    For(i, numSensecodons) {
        For(j, numSensecodons) {
            if(j != i) {
                if(aminoAcidmap[senseCodons[i]] != aminoAcidmap[senseCodons[j]]) {
                    roa += pi[i] * qmatbuff[i*numSensecodons+j];
                } else {
                    ros += pi[i] * qmatbuff[i*numSensecodons+j];
                }
            }
            else continue;
        }
    } if(n_w_catg == 1) {
        val = (roa*0.21)/((1-roa)*0.79);
    } else {
        val = (roa*0.21) / (ros*0.79);
    }
    return val;
}

// Returns a single factor for a given site that includes both values for kappa
// and/or omega. If neither of those is needed, return value will simply be 1.0 .
phydbl Kappa_Omega_Factor(int senseCodoni, int senseCodonj, model* mod, int cat,phydbl omega) {//modified by Ken 18/8
    phydbl val = 1.0;
    phydbl *kappa = mod->pkappa;
   // phydbl omega;
    int numSensecodons = mod->ns;
    ts_and_tv *mat = mod->io->structTs_and_Tv;
   
    switch(mod->omegaSiteVar) {
        case DM0:
            //omega = omega;
            break;
        case NOOMEGA:
            omega = 1.0;
            break;
        default:
            omega = mod->omegas[ cat ];
            break;
    }

    // First deal with kappa...
        if(mod->initqrates != NOINITMAT) {
            switch(mod->io->kappaECM) {
                case kap1: { 
                    val = 1.0;
                    break;
                }
                case kap2: {
                    switch(mat[numSensecodons * senseCodoni + senseCodonj].nts) {
                        case 0: val = 1.0; break;
                        case 1: val = kappa[0]; break;
                        case 2: val = kappa[0]*kappa[0]; break;
                        case 3: val = kappa[0]*kappa[0]*kappa[0]; break;
                        default: break;
                    }
                    break;
                }
                case kap3: {
                    switch(mat[numSensecodons * senseCodoni + senseCodonj].ntv) {
                        case 0: val = 1.0; break;
                        case 1: val = kappa[0]; break;
                        case 2: val = kappa[0] * kappa[0]; break;
                        case 3: val = kappa[0] * kappa[0] * kappa[0]; break;
                        default: break;
                    }
                    break;
                }
                case kap4:  {
                    switch(mat[numSensecodons * senseCodoni + senseCodonj].nts) {
                        case 0: val = 1.0; break;
                        case 1: val = kappa[0]; break;
                        case 2: val = kappa[0] * kappa[0]; break;
                        case 3: val = kappa[0] * kappa[0] * kappa[0]; break;
                        default: break;
                    }
                    switch(mat[numSensecodons * senseCodoni + senseCodonj].ntv) {
                        case 0: val *= 1.0; break;
                        case 1: val *= kappa[1]; break;
                        case 2: val *= kappa[1] * kappa[1]; break;
                        case 3: val *= kappa[1] * kappa[1] * kappa[1]; break;
                        default: break;
                    }
                    break;
                }
                case kap5: {
                    val = kappa[ mat[numSensecodons * senseCodoni + senseCodonj].caseK07 ];
                    break;
                }
                case kap6: {
                    if(mat[numSensecodons * senseCodoni + senseCodonj].ndiff > 1) {
                        val = kappa[0];
                    }
                    break;
                }
                default:
                    break;
            }
        } else {
            if(mat[numSensecodons * senseCodoni + senseCodonj].ndiff == 1) {
                if(mat[numSensecodons * senseCodoni + senseCodonj].nts) {
                    val = mod->kappa;
                } else {
                    val = 1.0;
                }
            } else {
                val = 0.0;
            }
        }
    // ... Then with omega:
    if(aminoAcidmap[ senseCodons[senseCodoni] ] != aminoAcidmap[ senseCodons[senseCodonj] ]) {
        val *= omega;
    }
    return val;
}

/*********************************************************/
void PMat_CODON(phydbl l, model *mod, int pos, phydbl *Pij) //!< Added by Marcelo.
{
    int n = mod->ns, nn = n*n, nw=mod->n_w_catg, i, g, j, k, catg;
    
    if(mod->io->expm==EIGEN)
    {
        phydbl *U, *V, *R, *P;
        phydbl *expt, expt_n, *pP;
        phydbl *uexpt, uexpt_n;
        
        P=Pij;
        
#if defined OMP || defined BLAS_OMP
        
        expt  = (phydbl *)malloc(n*nw*sizeof(phydbl));
        uexpt = (phydbl *)malloc(nn*nw*sizeof(phydbl));
        
#else
        
        expt  = mod->eigen->e_val_im;
        uexpt = mod->eigen->r_e_vect_im;
        
#endif
        
        U     = mod->eigen->r_e_vect;
        V     = mod->eigen->l_e_vect;
        R     = mod->eigen->e_val;
        
        For(catg,nw)
        {
#if defined BLAS || defined BLAS_OMP
            
            /* compute EXP(D/mr * l) into mat_eDmrl */
            For(k,n) expt[k+catg*n] = (phydbl)exp(R[k+catg*n]*l);
            
            /* multiply Vr*EXP(D/mr*l)*Vi into Pij */
            For (i,n) For (k,n) uexpt[i*n+k+catg*nn] = U[i*n+k+catg*nn] * expt[k+catg*n];
            
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, uexpt+catg*nn, n, V+catg*nn, n, 0.0, P+catg*nn, n);
            
#else
            
            for (k=0,zero(P+catg*nn,nn); k<n; k++) for (i=0,pP=P+catg*nn,expt_n=exp(R[k+catg*n]*l); i<n; i++) for (j=0,uexpt_n=U[i*n+k+catg*nn]*expt_n; j<n; j++) *pP++ += uexpt_n*V[k*n+j+catg*nn];
            
#endif 
        }
        
        For(i,nn*nw) if(Pij[i]<SMALL_PIJ) Pij[i]=SMALL_PIJ;
        
#if defined OMP || defined BLAS_OMP
        
        free(expt); 
        free(uexpt);
        
#endif
    }
    else if(mod->io->expm==TAYLOR) 
    {
    	if(mod->nparts > 1){
    		printf("TAYLOR APPROXIMATION DOESN'T WORK WITH PARITTIONED MODELS IN IGPHYML\n");
    		exit(EXIT_FAILURE);
    	}
        phydbl *Q2;
        int m;
        Q2=mod->A2_part[0]; //won't work with partitioned models
        
        For(i,nn*nw) Pij[i]=0.0;
        
        For(catg,nw)
        {
            For(m,mod->io->n_termsTaylor)
            {
                For(i,nn) Pij[i+catg*nn]+=pow(l,m)*Q2[i+catg*nn*15+m*nn]/(phydbl)myFactorial(m);
            }
        } 
        For(i,nn*nw) if(Pij[i]<SMALL_PIJ) Pij[i]=SMALL_PIJ;
    }
    else  if(mod->io->expm==SSPADE)
    {
        phydbl norm_1_Q, *Q, *A, *B, *F, theta[5] = {1.495585217958292e-002, 2.539398330063230e-001, 9.504178996162932e-001, 2.097847961257068e+000, 5.371920351148152e+000}; 
        int m_vals[5]={3, 5, 7, 9, 13}, z=5;
        
#if defined OMP || defined BLAS_OMP
        
        B=(phydbl *)malloc(nn*sizeof(phydbl));
        A=(phydbl *)malloc(nn*sizeof(phydbl));
        
#endif
        
        int modeli;
		for(modeli=0;modeli<mod->nomega_part;modeli++){ //Added by Ken 22/8
        For(catg,nw)
        {

            pos = nn*catg;
            Q   = mod->qmat_part[modeli]+pos;
            F   = Pij+pos;
            
#if defined OMP || defined BLAS_OMP
            
#else
            
            B = mod->matAux_part[modeli]+pos;
            A = mod->qmat_buff_part[modeli]+pos; //? 22/8 Ken
            
#endif
            
            For(i,nn) A[i]=Q[i]*l;
            
#if defined BLAS || defined BLAS_OMP
            
            norm_1_Q=dlange_("1", &n, &n, A, &n, NULL);
            
#else
            
            phydbl *max;
            max=(phydbl*)mCalloc(n,sizeof(phydbl));
            For(i,n) For(j,n) max[j]+=fabs(A[i*n+j]);
            norm_1_Q=DBL_MIN;
            For(i,n) if(max[i]>norm_1_Q) norm_1_Q=max[i];
            free(max);
            
#endif
            
            if(norm_1_Q<=theta[z-1])
            {
                For(i,z)
                {
                    if(norm_1_Q<=theta[i])
                    {
                        PadeApprox(n, nn, A, mod, F, pos, l, m_vals[i],modeli);
                        break;
                    }
                }
            }
            else
            {
                phydbl Mantissa, sFactor, num;
                int Exponent;
                num=norm_1_Q/theta[z-1];
                Mantissa=frexp(num,&Exponent);
                Exponent-=(Mantissa==0.5);
                sFactor=1.0/pow(2,Exponent);
                
                For(i,nn) A[i]*=sFactor;
                
                PadeApprox(n, nn, A, mod, F, pos, l, m_vals[z-1],modeli);
                
#if defined BLAS || defined BLAS_OMP
                
                Fors(i,Exponent,2)
                {
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, F, n, F, n, 0.0, B, n);
                    
                    if(i+1<Exponent)
                    {
                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, B, n, B, n, 0.0, F, n);
                    }
                    else 
                    {
                        cblas_dcopy(nn,B,1,F,1);
                        break;
                    }
                }
                
#else
                
                Fors(g,Exponent,2)
                {
                    For(i,nn) B[i]=0.0;
                    For(i,n) For(j,n) For(k,n) B[n*i+j] += F[i*n+k] * F[k*n+j];
                    
                    if(g+1<Exponent)
                    {
                        For(i,nn) F[i]=0.0;
                        For(i,n) For(j,n) For(k,n) F[n*i+j] += B[i*n+k] * B[k*n+j];
                    }
                    else 
                    {
                        For(i,nn) F[i]=B[i];
                        break;
                    }
                }
                
#endif
            }
        }
        For(i,nn*nw) if(Pij[i]<SMALL_PIJ) Pij[i]=SMALL_PIJ; //!< Correct for too small entries.
        
#if defined OMP || defined BLAS_OMP 
        
        free(A); 
        free(B);
        
#endif
    }// for(modeli=0..
    }
}


/*********************************************************/
void PMat_CODON_part(phydbl l, model *mod, int pos, phydbl *Qmat, phydbl *Pij, int modeli) //!< Added by Marcelo.
{
    int n = mod->ns, nn = n*n, nw=mod->n_w_catg, i, g, j, k, catg;

    if(mod->io->expm==EIGEN)
    {
        phydbl *U, *V, *R, *P;
        phydbl *expt, expt_n, *pP;
        phydbl *uexpt, uexpt_n;

        P=Pij;

#if defined OMP || defined BLAS_OMP

        expt  = (phydbl *)malloc(n*nw*sizeof(phydbl));
        uexpt = (phydbl *)malloc(nn*nw*sizeof(phydbl));

#else

        expt  = mod->eigen->e_val_im;
        uexpt = mod->eigen->r_e_vect_im;

#endif

        U     = mod->eigen->r_e_vect;
        V     = mod->eigen->l_e_vect;
        R     = mod->eigen->e_val;

        For(catg,nw)
        {
#if defined BLAS || defined BLAS_OMP

            /* compute EXP(D/mr * l) into mat_eDmrl */
            For(k,n) expt[k+catg*n] = (phydbl)exp(R[k+catg*n]*l);

            /* multiply Vr*EXP(D/mr*l)*Vi into Pij */
            For (i,n) For (k,n) uexpt[i*n+k+catg*nn] = U[i*n+k+catg*nn] * expt[k+catg*n];

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, uexpt+catg*nn, n, V+catg*nn, n, 0.0, P+catg*nn, n);

#else

            for (k=0,zero(P+catg*nn,nn); k<n; k++) for (i=0,pP=P+catg*nn,expt_n=exp(R[k+catg*n]*l); i<n; i++) for (j=0,uexpt_n=U[i*n+k+catg*nn]*expt_n; j<n; j++) *pP++ += uexpt_n*V[k*n+j+catg*nn];

#endif
        }

        For(i,nn*nw) if(Pij[i]<SMALL_PIJ) Pij[i]=SMALL_PIJ;

#if defined OMP || defined BLAS_OMP

        free(expt);
        free(uexpt);

#endif
    }
    else if(mod->io->expm==TAYLOR)
    {
    	if(mod->nparts > 1){
    		printf("TAYLOR APPROX DOESNT WORK WITH PARTITIONED MODELS IN IGPHYML YET\n");
    		exit(EXIT_FAILURE);
    	}
        phydbl *Q2;
        int m;
        Q2=mod->A2_part[0];

        For(i,nn*nw) Pij[i]=0.0;

        For(catg,nw)
        {
            For(m,mod->io->n_termsTaylor)
            {
                For(i,nn) Pij[i+catg*nn]+=pow(l,m)*Q2[i+catg*nn*15+m*nn]/(phydbl)myFactorial(m);
            }
        }
        For(i,nn*nw) if(Pij[i]<SMALL_PIJ) Pij[i]=SMALL_PIJ;
    }
    else  if(mod->io->expm==SSPADE)
    {
        phydbl norm_1_Q, *Q, *A, *B, *F, theta[5] = {1.495585217958292e-002, 2.539398330063230e-001, 9.504178996162932e-001, 2.097847961257068e+000, 5.371920351148152e+000};
        int m_vals[5]={3, 5, 7, 9, 13}, z=5;

#if defined OMP || defined BLAS_OMP

        B=(phydbl *)malloc(nn*sizeof(phydbl));
        A=(phydbl *)malloc(nn*sizeof(phydbl));

#endif

        For(catg,nw)
        {

            pos = nn*catg;
            Q   = Qmat+pos;
            F   = Pij+pos;

#if defined OMP || defined BLAS_OMP

#else

            B = mod->matAux_part[modeli]+pos;
            A = mod->qmat_buff_part[modeli]+pos;

#endif

            For(i,nn) A[i]=Q[i]*l;

#if defined BLAS || defined BLAS_OMP

            norm_1_Q=dlange_("1", &n, &n, A, &n, NULL);

#else

            phydbl *max;
            max=(phydbl*)mCalloc(n,sizeof(phydbl));
            For(i,n) For(j,n) max[j]+=fabs(A[i*n+j]);
            norm_1_Q=DBL_MIN;
            For(i,n) if(max[i]>norm_1_Q) norm_1_Q=max[i];
            free(max);

#endif

            if(norm_1_Q<=theta[z-1])
            {
                For(i,z)
                {
                    if(norm_1_Q<=theta[i])
                    {
                        PadeApprox(n, nn, A, mod, F, pos, l, m_vals[i],modeli);
                        break;
                    }
                }
            }
            else
            {
                phydbl Mantissa, sFactor, num;
                int Exponent;
                num=norm_1_Q/theta[z-1];
                Mantissa=frexp(num,&Exponent);
                Exponent-=(Mantissa==0.5);
                sFactor=1.0/pow(2,Exponent);

                For(i,nn) A[i]*=sFactor;

                PadeApprox(n, nn, A, mod, F, pos, l, m_vals[z-1],modeli);

#if defined BLAS || defined BLAS_OMP

                Fors(i,Exponent,2)
                {
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, F, n, F, n, 0.0, B, n);

                    if(i+1<Exponent)
                    {
                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, B, n, B, n, 0.0, F, n);
                    }
                    else
                    {
                        cblas_dcopy(nn,B,1,F,1);
                        break;
                    }
                }

#else

                Fors(g,Exponent,2)
                {
                    For(i,nn) B[i]=0.0;
                    For(i,n) For(j,n) For(k,n) B[n*i+j] += F[i*n+k] * F[k*n+j];

                    if(g+1<Exponent)
                    {
                        For(i,nn) F[i]=0.0;
                        For(i,n) For(j,n) For(k,n) F[n*i+j] += B[i*n+k] * B[k*n+j];
                    }
                    else
                    {
                        For(i,nn) F[i]=B[i];
                        break;
                    }
                }

#endif
            }
        }
        For(i,nn*nw) if(Pij[i]<SMALL_PIJ) Pij[i]=SMALL_PIJ; //!< Correct for too small entries.

#if defined OMP || defined BLAS_OMP

        free(A);
        free(B);

#endif
    }
}

/***********************************************************/
void PadeApprox(int n, int nn, phydbl *A, model *mod, phydbl *F, int pos, phydbl len, int m, int modeli) //!<  Added by Marcelo.
{
    phydbl c3[4]={120.0, 60.0, 12.0, 1.0}, c5[6]={30240.0, 15120.0, 3360.0, 420.0, 30.0, 1.0}, c7[8]={17297280.0, 8648640.0, 1995840.0, 277200.0, 25200.0, 1512.0, 56.0, 1.0}, c9[10]={17643225600.0, 8821612800.0, 2075673600.0, 302702400.0, 30270240.0, 2162160.0, 110880.0, 3960.0, 90.0, 1.0}, c13[14]={64764752532480000.0, 32382376266240000.0, 7771770303897600.0, 1187353796428800.0,  129060195264000.0, 10559470521600.0, 670442572800.0, 33522128640.0, 1323241920.0, 40840800.0, 960960.0, 16380.0, 182.0, 1.0};
    phydbl *ptCoeffPade, *U, *V, *A0, *A2, *A4, *A6, *A8, *Apowers, *Aux, len2, len4, len6, len8;
    phydbl lens[5];
    int j, o, p, q, *Ipiv;
    int info=0;
    
    switch(m)
    {	
        case 3:  ptCoeffPade=c3;  break;
        case 5:  ptCoeffPade=c5;  break;
        case 7:  ptCoeffPade=c7;  break;
        case 9:  ptCoeffPade=c9;  break;
        case 13: ptCoeffPade=c13; break;
        default: Warn_And_Exit("Invalid Pade Coefficient.");break;
    }
    
#if defined OMP || defined BLAS_OMP
    
    Aux=(phydbl *)malloc(nn*sizeof(phydbl));
    U=(phydbl *)malloc(nn*sizeof(phydbl));
    V=(phydbl *)malloc(nn*sizeof(phydbl));
    Ipiv=(int *)malloc(n*sizeof(int));
    
#else
    
    Aux =mod->matAux_part[modeli]+pos;
    U   =mod->U_part[modeli]+pos;
    V   =mod->V_part[modeli]+pos;
    Ipiv=mod->ipiv_part[modeli]+(pos/n);
    
#endif
    
    A0      = mod->A0_part[modeli];
    A2      = mod->A2_part[modeli]+pos;
    A4      = mod->A4_part[modeli]+pos;
    A6      = mod->A6_part[modeli]+pos;
    A8      = mod->A8_part[modeli]+pos;
    Apowers = mod->Apowers_part[modeli]+5*pos;
    
    len2=len*len;
    len4=len2*len2;
    len6=len2*len4;
    len8=len2*len6;
    lens[0]=1.0; lens[1]=len2; lens[2]=len4; lens[3]=len6; lens[4]=len8;
    
    switch(m)
    {
        case 3:  
        case 5:  
        case 7:  
        case 9:  
        {
#if defined BLAS || defined BLAS_OMP
            
            //For(j,nn) U[j]=0;
            For(j,nn) U[j]=lens[m/2]*ptCoeffPade[m]*Apowers[j+nn*(m/2)];      
            //for(j=m;j>0;j-=2) cblas_daxpy(nn, lens[j/2]*ptCoeffPade[j], Apowers+nn*(j/2), 1, U, 1);
            for(j=m-2;j>0;j-=2) cblas_daxpy(nn, lens[j/2]*ptCoeffPade[j], Apowers+nn*(j/2), 1, U, 1);
            
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, A, n, U, n, 0.0, Aux, n);
            
            //For(j,nn) V[j]=0;
            For(j,nn) V[j]=lens[(m-1)/2]*ptCoeffPade[m-1]*Apowers[j+nn*((m-1)/2)];
            //for(j=m-1;j>-1;j-=2) cblas_daxpy(nn, lens[j/2]*ptCoeffPade[j], Apowers+nn*(j/2), 1, V, 1);
            for(j=m-3;j>-1;j-=2) cblas_daxpy(nn, lens[j/2]*ptCoeffPade[j], Apowers+nn*(j/2), 1, V, 1);
            
            cblas_dcopy(nn, Aux, 1, F, 1);
            cblas_daxpy(nn,  1, V, 1, F, 1);
            
            cblas_daxpy(nn, -1, Aux, 1, V, 1);
            
            dgesv_(&n, &n, V, &n, Ipiv, F, &n, &info);
            if(info!=0) Warn_And_Exit("Matrix cannot be inverted."); 
            
#else
            
            For(j,nn) U[j]=0;
            
            for(j=m;j>0;j-=2) For(o,nn) U[o]+=lens[j/2]*ptCoeffPade[j]*Apowers[nn*(j/2)+o];
            
            For(o,nn) Aux[o]=0.0;
            For(o,n) For(p,n) For(q,n) Aux[n*o+p] += A[o*n+q] * U[q*n+p];
            
            For(j,nn) V[j]=0;
            
            for(j=m-1;j>-1;j-=2) For(o,nn) V[o]+=lens[j/2]*ptCoeffPade[j]*Apowers[nn*(j/2)+o];
            
            For(o,nn) F[o]=Aux[o];
            For(o,nn) F[o]+=V[o];
            
            For(o,nn) V[o]-=Aux[o];
            
            if(MyGaussElimination_gNxRHS(V,F,n,n)!=0) Warn_And_Exit("Matrix cannot be inverted."); 
            
#endif
            
            break;
        }
        case 13: 
        {
#if defined BLAS || defined BLAS_OMP
            
            //For(j,nn) U[j]=0;
            For(j,nn) U[j]=len6*ptCoeffPade[13]*A6[j];
            //!< MATLAB (index -1 to C std) U = A * (A6*(c(14)*A6 + c(12)*A4 + c(10)*A2) + c(8)*A6 + c(6)*A4 + c(4)*A2 + c(2)*eye(n,classA));
            //cblas_daxpy(nn, len6*ptCoeffPade[13], A6, 1, U, 1);
            cblas_daxpy(nn, len4*ptCoeffPade[11], A4, 1, U, 1);
            cblas_daxpy(nn, len2*ptCoeffPade[ 9], A2, 1, U, 1);
            
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, len6, A6, n, U, n, 0.0, Aux, n);
            
            cblas_daxpy(nn, len6*ptCoeffPade[ 7], A6, 1, Aux, 1);
            cblas_daxpy(nn, len4*ptCoeffPade[ 5], A4, 1, Aux, 1);
            cblas_daxpy(nn, len2*ptCoeffPade[ 3], A2, 1, Aux, 1);
            cblas_daxpy(nn, ptCoeffPade[ 1], A0, 1, Aux, 1);
            cblas_dcopy(nn,Aux,1,U,1);
            
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, A, n, U, n, 0.0, Aux, n);
            
            //For(j,nn) V[j]=0;
            For(j,nn) V[j]=len6*ptCoeffPade[12]*A6[j];
            //!< MATLAB (index -1 to C std) V = A6*(c(13)*A6 + c(11)*A4 + c(9)*A2) + c(7)*A6 + c(5)*A4 + c(3)*A2 + c(1)*eye(n,classA);
            //cblas_daxpy(nn, len6*ptCoeffPade[12], A6, 1, V, 1);
            cblas_daxpy(nn, len4*ptCoeffPade[10], A4, 1, V, 1);
            cblas_daxpy(nn, len2*ptCoeffPade[ 8], A2, 1, V, 1);
            
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, len6, A6, n, V, n, 0.0, U, n);
            
            cblas_daxpy(nn, len6*ptCoeffPade[ 6], A6, 1, U, 1);
            cblas_daxpy(nn, len4*ptCoeffPade[ 4], A4, 1, U, 1);
            cblas_daxpy(nn, len2*ptCoeffPade[ 2], A2, 1, U, 1);
            cblas_daxpy(nn, ptCoeffPade[ 0], A0, 1, U, 1);
            
            cblas_dcopy(nn, Aux, 1, F, 1);
            cblas_daxpy(nn,  1, U, 1, F, 1);
            
            cblas_daxpy(nn, -1, Aux, 1, U, 1);      
            
            dgesv_(&n, &n, U, &n, Ipiv, F, &n, &info);
            if(info!=0) Warn_And_Exit("Matrix cannot be inverted."); 
            
#else
            
            For(j,nn) U[j]=0;
            
            For(j,nn) U[j]+=len6*ptCoeffPade[13]*A6[j];
            For(j,nn) U[j]+=len4*ptCoeffPade[11]*A4[j];
            For(j,nn) U[j]+=len2*ptCoeffPade[ 9]*A2[j];
            
            For(o,nn) Aux[o]=0.0;
            For(o,n) For(p,n) For(q,n) Aux[n*o+p] +=  len6*A6[o*n+q] *  U[q*n+p];
            
            For(j,nn) Aux[j]+=len6*ptCoeffPade[ 7]*A6[j];
            For(j,nn) Aux[j]+=len4*ptCoeffPade[ 5]*A4[j];
            For(j,nn) Aux[j]+=len2*ptCoeffPade[ 4]*A2[j];
            For(j,nn) Aux[j]+=ptCoeffPade[ 1]*A0[j];
            
            For(o,nn) U[o]=0.0;
            For(o,n) For(p,n) For(q,n) U[n*o+p] +=  A[o*n+q] *  Aux[q*n+p];
            
            For(j,nn) V[j]=0;
            
            For(j,nn) V[j]+=len6*ptCoeffPade[12]*A6[j];
            For(j,nn) V[j]+=len4*ptCoeffPade[10]*A4[j];
            For(j,nn) V[j]+=len2*ptCoeffPade[ 8]*A2[j];
            
            For(o,nn) Aux[o]=0.0;
            For(o,n) For(p,n) For(q,n) Aux[n*o+p] +=  len6*A6[o*n+q] *  V[q*n+p];
            
            For(j,nn) Aux[j]+=len6*ptCoeffPade[ 6]*A6[j];
            For(j,nn) Aux[j]+=len4*ptCoeffPade[ 4]*A4[j];
            For(j,nn) Aux[j]+=len2*ptCoeffPade[ 2]*A2[j];
            For(j,nn) Aux[j]+=ptCoeffPade[ 0]*A0[j];
            
            For(o,nn) F[o]=U[o];
            For(o,nn) F[o]+=Aux[o];
            
            For(o,nn) Aux[o]-=U[o]; 
            
            if(MyGaussElimination_gNxRHS(Aux,F,n,n)!=0) Warn_And_Exit("Matrix cannot be inverted.");     
            
#endif
            
            
        }
        default:break;
    }
    
#if defined OMP || defined BLAS_OMP
    
    free(Aux);
    free(U);
    free(V);
    free(Ipiv);
    
#endif
}
/***********************************************************/
phydbl Scale_Q_MatrixManyOmegaModels(phydbl *qmat,  phydbl *pi, phydbl *freq, int n_w_catg, int numSenseCodons, model *mod) //!< Added by Marcelo.
{
    /*! Normalize Q matrix.*/
    int i, k, nn=numSenseCodons*numSenseCodons;
    phydbl mu, sum;
    
    mu=0.0;
    For(k, n_w_catg)
    {
        sum=0.0;
        For(i,numSenseCodons)
        {
            sum+=(-qmat[k*nn + numSenseCodons*i+i])*pi[i];
        }
        mu+=sum*freq[k];
    } 
    
    For(i,n_w_catg*nn) mod->qmatScaled[i]=mod->qmat[i]/mu;
    
    return(mu);
}
/***********************************************************/
phydbl Scale_P_MatrixManyOmegaModels(phydbl *Pij,  phydbl *pi, phydbl *freq, int n_w_catg, int numSenseCodons, model *mod) //!< Added by Marcelo.
{
    /*! Normalize Q matrix.*/
    int i, k, nn=numSenseCodons*numSenseCodons;
    phydbl mu, sum;
    
    mu=0.0;
    For(k, n_w_catg)
    {
        sum=0.0;
        For(i,numSenseCodons)
        {
            sum+=(Pij[k*nn + numSenseCodons*i+i])*pi[i];
        }
        mu+=sum*freq[k];
    } 
    
    For(i,n_w_catg*nn) Pij[i]/=mu;
    
    return(mu);
}
/***********************************************************/
void PMat_CODON_Pairwise(phydbl l, phydbl *Pij,  phydbl *U, phydbl *V, phydbl *R, int n, phydbl *uexpt, phydbl *expt) //!< Added by Marcelo.
{
    int nn=n*n;
    int i, j, k;
    // int g;
    phydbl *pP, expt_n, uexpt_n;
    
#if defined BLAS || defined BLAS_OMP
    
    /* compute POW(EXP(D/mr),l) into mat_eDmrl */
    For(k,n) expt[k] = (phydbl)exp(R[k]*l);
    
    /* multiply Vr*POW(EXP(D/mr),l)*Vi into Pij */
    For (i,n) For (k,n) uexpt[i*n+k] = U[i*n+k] * expt[k];
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, uexpt, n, V, n, 0.0, Pij, n);
    
#else
    
    For(i,nn) Pij[i] = 0.0;
    for (k=0; k<n; k++) for (i=0,pP=Pij,expt_n=exp(R[k]*l); i<n; i++) for (j=0,uexpt_n=U[i*n+k]*expt_n; j<n; j++) *pP++ += uexpt_n*V[k*n+j];
    
#endif
    
    For(i,nn) if(Pij[i]<SMALL_PIJ) Pij[i]=SMALL_PIJ;
    
}

/*********************************************************/

void print_matrix(phydbl *mat, int ns, char *s,int row) //!< Added by Marcelo.
{ //!< Useful for debugging.
    int i,j;
    FILE *out;
    out=fopen(s,"w");
    For(i,ns)
    {
        For(j,ns) fprintf(out,"%.16f ",mat[i*ns+j]);
        if(!row) fprintf(out,"\n ");
    }
    fprintf(out,"\n ");
    fclose(out);
}
/*********************************************************/

void print_array(phydbl *array, int ns, char *s) //!< Added by Marcelo.
{ //!< Useful for debugging.
    int i;
    // int j;
    FILE *out;
    out=fopen(s,"w");
    For(i,ns)
    {
        fprintf(out,"%.11f ",array[i]);
    }
    fprintf(out,"\n ");
    fclose(out);
}

/*******************************************************/

void Update_Rate_Matrix_PCAModel(model *mod) { //!<Added by Marcelo.
    int i, j, k;
    int numSensecodons = mod->ns;
    phydbl value;
    
    For(i, 64) For(j, 64) mod->userRates[i][j] = 0.0;
    
    // printf("%f, %f\n", mod->pcsC[0], mod->pcsC[1]);
    
    For(i, numSensecodons) {
        For(j, i) {
            value = 0;
            For(k, mod->npcs) value += mod->pcsC[k] * pcaP[k*numSensecodons*numSensecodons + i*numSensecodons +j];
            value *= pcaS[i][j];
            value += pcaM[i][j];
            value = (value>0) ? value : 0;
            mod->userRates[senseCodons[i]][senseCodons[j]] = value;
            mod->userRates[senseCodons[j]][senseCodons[i]] = value;
        }
    }
}
/*******************************************************/
