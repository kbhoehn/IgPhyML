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

#include "stats.h"


/*********************************************************/
/* RANDOM VARIATES GENERATORS */
/*********************************************************/

/*********************************************************************/
/* A C-function for TT800 : July 8th 1996 Version */
/* by M. Matsumoto, email: matumoto@math.keio.ac.jp */
/* tt800() generate one pseudorandom number with double precision */
/* which is uniformly distributed on [0,1]-interval */
/* for each call.  One may choose any initial 25 seeds */
/* except all zeros. */

/* See: ACM Transactions on Modelling and Computer Simulation, */
/* Vol. 4, No. 3, 1994, pages 254-266. */

phydbl tt800()
{
  int M=7;
  unsigned long y;
  static int k = 0;
  static unsigned long x[25]={ /* initial 25 seeds, change as you wish */
    0x95f24dab, 0x0b685215, 0xe76ccae7, 0xaf3ec239, 0x715fad23,
    0x24a590ad, 0x69e4b5ef, 0xbf456141, 0x96bc1b7b, 0xa7bdf825,
    0xc1de75b7, 0x8858a9c9, 0x2da87693, 0xb657f9dd, 0xffdc8a9f,
    0x8121da71, 0x8b823ecb, 0x885d05f5, 0x4e20cd47, 0x5a9ad5d9,
    0x512c0c03, 0xea857ccd, 0x4cc1d30f, 0x8891a8a1, 0xa6b7aadb
  };
  static unsigned long mag01[2]={ 
    0x0, 0x8ebfd028 /* this is magic vector `a', don't change */
  };
  if (k==25) { /* generate 25 words at one time */
    int kk;
    for (kk=0;kk<25-M;kk++) {
      x[kk] = x[kk+M] ^ (x[kk] >> 1) ^ mag01[x[kk] % 2];
    }
    for (; kk<25;kk++) {
      x[kk] = x[kk+(M-25)] ^ (x[kk] >> 1) ^ mag01[x[kk] % 2];
    }
    k=0;
  }
  y = x[k];
  y ^= (y << 7) & 0x2b5b2500; /* s and b, magic vectors */
  y ^= (y << 15) & 0xdb8b0000; /* t and c, magic vectors */
  y &= 0xffffffff; /* you may delete this line if word size = 32 */
  /* 
     the following line was added by Makoto Matsumoto in the 1996 version
     to improve lower bit's corellation.
     Delete this line to o use the code published in 1994.
  */
  y ^= (y >> 16); /* added to the 1994 version */
  k++;
  return(phydbl)(y / (unsigned long) 0xffffffff);
}

/*********************************************************************/

phydbl Uni()
{
  phydbl r; 
  r  = (phydbl)rand();
  r /= RAND_MAX;
/*   r = tt800(); */
  return r;
}

/*********************************************************************/

int Rand_Int(int min, int max)
{
/*   phydbl u;   */
/*   u = (phydbl)rand(); */
/*   u /=  (RAND_MAX); */
/*   u *= (max - min + 1); */
/*   u += min; */
/*   return (int)FLOOR(u); */

  int u;
  u = rand();
  return (u%(max+1-min)+min);

}

/*********************************************************/


/********************* random Gamma generator ************************
* Properties:
* (1) X = Gamma(alpha,lambda) = Gamma(alpha,1)/lambda
* (2) X1 = Gamma(alpha1,1), X2 = Gamma(alpha2,1) independent
*     then X = X1+X2 = Gamma(alpha1+alpha2,1)
* (3) alpha = k = integer then
*     X = Gamma(k,1) = Erlang(k,1) = -sum(log(Ui)) = -LOG(prod(Ui))
*     where U1,...Uk iid uniform(0,1)
*
* Decompose alpha = k+delta with k = [alpha], and 0<delta<1
* Apply (3) for Gamma(k,1)
* Apply Ahrens-Dieter algorithm for Gamma(delta,1)
*/
 
phydbl Ahrensdietergamma(phydbl alpha)
{
  phydbl x = 0.;

  if (alpha>0.) 
    {
      phydbl y = 0.;
      phydbl b = (alpha+EXP(1.))/EXP(1.);
      phydbl p = 1./alpha;
      int go = 0;
      while (go==0) 
	{
	  phydbl u = Uni();
	  phydbl w = Uni();
	  phydbl v = b*u;
	  if (v<=1.) 
	    {
	      x = POW(v,p);
	      y = EXP(-x);
	    }
	  else 
	    {
	      x = -LOG(p*(b-v));
	      y = POW(x,alpha-1.);
	    }
	  go = (w<y); // x is accepted when go=1
	}
    }
  return x;
}

/*********************************************************/

phydbl Rgamma(phydbl shape, phydbl scale)
{
  int i;
  phydbl x1 = 0.;
  phydbl delta = shape;
  if (shape>=1.) 
    {
      int k = (int)FLOOR(shape);
      delta = shape - k;
      phydbl u = 1.;
      for (i=0; i<k; i++)
	u *= Uni();
      x1 = -LOG(u);
    }
  phydbl x2 = Ahrensdietergamma(delta);
  return (x1 + x2)*scale;
}

/*********************************************************/

phydbl Rexp(phydbl lambda)
{
  return -LOG(Uni()+1.E-30)/lambda;
}

/*********************************************************/

phydbl Rnorm(phydbl mean, phydbl sd)
{
  /* Box-Muller transformation */
  phydbl u1, u2;

  u1=Uni();
  u2=Uni();
  u1 = SQRT(-2.*LOG(u1))*COS(6.28318530717959f*u2);

  /* Polar */
/*   phydbl d,x,y; */

/*   do */
/*     { */
/*       u1=Uni(); */
/*       u2=Uni(); */
/*       x = 2.*u1-1.; */
/*       y = 2.*u2-1.; */
/*       d = x*x + y*y; */
/*       if(d>.0 && d<1.) break; */
/*     } */
/*   while(1); */
/*   u1 = x*SQRT((-2.*LOG(d))/d); */


  return u1*sd+mean;
}

/*********************************************************/

phydbl *Rnorm_Multid(phydbl *mu, phydbl *cov, int dim)
{
  phydbl *L,*x,*y;
  int i,j;
  
  x = (phydbl *)mCalloc(dim,sizeof(phydbl));
  y = (phydbl *)mCalloc(dim,sizeof(phydbl));

  L = (phydbl *)Cholesky_Decomp(cov,dim);

  For(i,dim) x[i]=Rnorm(0.0,1.0);
  For(i,dim) For(j,dim) y[i] += L[i*dim+j]*x[j];
  For(i,dim) y[i] += mu[i];

  free(L);
  free(x);

  return(y);
}

/*********************************************************/

phydbl Rnorm_Trunc_Inverse(phydbl mean, phydbl sd, phydbl min, phydbl max, int *error)
{

  phydbl u, ret_val,eps;
  phydbl z,rz;
  phydbl z_min,z_max;
  phydbl cdf_min, cdf_max;

  z      = 0.0;
  rz     = 0.0;
  u      = -1.0;
  *error = 0;
  
  if(sd < 1.E-100)
    {
      PhyML_Printf("\n. Small variance detected in Rnorm_Trunc.");
      PhyML_Printf("\n. mean=%f sd=%f min=%f max=%f",mean,sd,min,max);
      *error = 1;
      return -1.0;
    }

  z_min = (min - mean)/sd;
  z_max = (max - mean)/sd;

  eps = (z_max-z_min)/1E+6;


  /*       Simple inversion method. Seems to work well. Needs more thorough testing though... */
  cdf_min = Pnorm(z_min,0.0,1.0);
  cdf_max = Pnorm(z_max,0.0,1.0);
  u = cdf_min + (cdf_max-cdf_min) * Uni();
  z = PointNormal(u);
	
  if((z < z_min-eps) || (z > z_max+eps))
    {
      *error = 1;
      PhyML_Printf("\n. Numerical precision issue detected in Rnorm_Trunc.");
      PhyML_Printf("\n. z = %f",z);
      PhyML_Printf("\n. mean=%f sd=%f z_min=%f z_max=%f min=%f max=%f",mean,sd,z_min,z_max,min,max);
      ret_val = (max - min)/2.;
      Exit("\n");
    }
  
  ret_val = z*sd+mean;
  
  return ret_val;
}

/*********************************************************/

phydbl Rnorm_Trunc(phydbl mean, phydbl sd, phydbl min, phydbl max, int *error)
{

  phydbl u, ret_val,eps;
  int iter;
  phydbl z,rz;
  phydbl z_min,z_max;

  z      = 0.0;
  rz     = 0.0;
  u      = -1.0;
  *error = 0;
  
  if(sd < 1.E-100)
    {
      PhyML_Printf("\n. Small variance detected in Rnorm_Trunc.");
      PhyML_Printf("\n. mean=%f sd=%f min=%f max=%f",mean,sd,min,max);
      *error = 1;
      return -1.0;
    }

  if(max < min)
    {
      PhyML_Printf("\n. Max < Min");
      PhyML_Printf("\n. mean=%f sd=%f min=%f max=%f",mean,sd,min,max);
      *error = 1;
      return -1.0;
    }

  z_min = (min - mean)/sd;
  z_max = (max - mean)/sd;

  eps = (z_max-z_min)/1E+6;

  /* Damien and Walker (2001) method */
  phydbl y,slice_min,slice_max;

/*   if((z_min < -10.) && (z_max > +10.)) /\* cdf < 1.E-6, we should be safe. *\/ */
/*     { */
/*       z = Rnorm(0.0,1.0); */
/*     } */
/*   else */
/*     { */

      iter = 0;
      do
	{
	  y   = Uni()*EXP(-(z*z)/2.);
	  slice_min = MAX(z_min,-SQRT(-2.*LOG(y)));
	  slice_max = MIN(z_max, SQRT(-2.*LOG(y)));
	  z   = Uni()*(slice_max - slice_min) + slice_min;
	  iter++;
	  if(iter > 1000) break;
	}
      while(slice_max < slice_min || iter < 10);

      if(iter > 1000)
	{
	  PhyML_Printf("\n. Too many iterations in Rnorm_Trunc...");
	  *error = 1;
	}

/*     } */

  /* Inverson method */
/*   phydbl cdf_min, cdf_max; */
/*   if((z_min < -10.) && (z_max > +10.)) /\* cdf < 1.E-6, we should be safe. *\/ */
/*     { */
/*       z = Rnorm(0.0,1.0); */
/*     } */
/*   else */
/*     { */
/* /\*       Simple inversion method. Seems to work well. Needs more thorough testing though... *\/ */
/*       cdf_min = Pnorm(z_min,0.0,1.0); */
/*       cdf_max = Pnorm(z_max,0.0,1.0); */
/*       u = cdf_min + (cdf_max-cdf_min) * Uni(); */
/*       z = PointNormal(u); */
/*     } */


   if((z < z_min-eps) || (z > z_max+eps))
    {
      *error = 1;
      PhyML_Printf("\n. Numerical precision issue detected in Rnorm_Trunc.");
      PhyML_Printf("\n. z = %f",z);
      PhyML_Printf("\n. mean=%f sd=%f z_min=%f z_max=%f min=%f max=%f",mean,sd,z_min,z_max,min,max);
      ret_val = (max - min)/2.;
      Exit("\n");
    }

  ret_val = z*sd+mean;

  return ret_val;
}

/*********************************************************/

phydbl *Rnorm_Multid_Trunc(phydbl *mean, phydbl *cov, phydbl *min, phydbl *max, int dim)
{
  int i,j;
  phydbl *L,*x, *u;
  phydbl up, low, rec;
  int err;
  
  u = (phydbl *)mCalloc(dim,sizeof(dim)); 
  x = (phydbl *)mCalloc(dim,sizeof(dim));
 
  L = Cholesky_Decomp(cov,dim);
  
  low = (min[0]-mean[0])/L[0*dim+0];
  up  = (max[0]-mean[0])/L[0*dim+0];
  u[0] = Rnorm_Trunc(0.0,1.0,low,up,&err);

  for(i=1;i<dim;i++)
    {
      rec = .0;
      For(j,i) rec += L[i*dim+j] * u[j];
      low  = (min[i]-mean[i]-rec)/L[i*dim+i];
      up   = (max[i]-mean[i]-rec)/L[i*dim+i];
      u[i] = Rnorm_Trunc(0.0,1.0,low,up,&err);
    }

  x = Matrix_Mult(L,u,dim,dim,dim,1);

/*   PhyML_Printf("\n>>>\n"); */
/*   For(i,dim) */
/*     { */
/*       For(j,dim) */
/* 	{ */
/* 	  PhyML_Printf("%10lf ",L[i*dim+j]); */
/* 	} */
/*       PhyML_Printf("\n"); */
/*     } */
/*   PhyML_Printf("\n"); */

/*   For(i,dim) PhyML_Printf("%f ",u[i]); */
/*   PhyML_Printf("\n"); */

  
/*   PhyML_Printf("\n"); */
/*   For(i,dim) PhyML_Printf("%10lf ",x[i]); */
/*   PhyML_Printf("\n<<<\n"); */

  For(i,dim) x[i] += mean[i];

  free(L);
  free(u);
  
  return x;
}

/*********************************************************/
/* DENSITIES / PROBA */
/*********************************************************/

phydbl Dnorm_Moments(phydbl x, phydbl mean, phydbl var)
{
  phydbl dens,sd,pi;

  pi = 3.141593;
  sd = SQRT(var);

  dens = 1./(SQRT(2*pi)*sd)*EXP(-((x-mean)*(x-mean)/(2.*sd*sd)));

  return dens;
}

/*********************************************************/

phydbl Dnorm(phydbl x, phydbl mean, phydbl sd)
{
  phydbl dens;
/*   dens = -(.5*LOG2PI+LOG(sd))  - .5*POW(x-mean,2)/POW(sd,2); */
/*   return EXP(dens); */
  
  dens = (M_1_SQRT_2_PI * EXP(-0.5*x*x))/sd; 
  return dens;
}

/*********************************************************/

phydbl Log_Dnorm(phydbl x, phydbl mean, phydbl sd, int *err)
{
  phydbl dens;

  *err = NO;

  x = (x-mean)/sd;
  
  dens = -(phydbl)LOG_SQRT_2_PI - x*x*0.5 - LOG(sd);

  if(dens < -BIG)
    {
      PhyML_Printf("\n. dens=%f -- x=%f mean=%f sd=%f\n",dens,x,mean,sd);
      *err = 1;
    }

  return dens;
}

/*********************************************************/

phydbl Log_Dnorm_Trunc(phydbl x, phydbl mean, phydbl sd, phydbl lo, phydbl up, int *err)
{
  phydbl log_dens;
  phydbl cdf_up = 0.0;
  phydbl cdf_lo = 0.0;

  *err = NO;

  log_dens = Log_Dnorm(x,mean,sd,err);

  if(*err == YES)
    {
      PhyML_Printf("\n. mean=%f sd=%f lo=%f up=%f cdf_lo=%G CDF_up=%G",mean,sd,lo,up,cdf_lo,cdf_up);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      *err = YES;
    }

  cdf_up = Pnorm(up,mean,sd);
  cdf_lo = Pnorm(lo,mean,sd);

  log_dens -= LOG(cdf_up - cdf_lo);

  if(isnan(log_dens) || isinf(FABS(log_dens)))
    {
      PhyML_Printf("\n. x=%f mean=%f sd=%f lo=%f up=%f cdf_lo=%G CDF_up=%G",x,mean,sd,lo,up,cdf_lo,cdf_up);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      *err = YES;
   }

  return log_dens;
}

/*********************************************************/

phydbl Dnorm_Trunc(phydbl x, phydbl mean, phydbl sd, phydbl lo, phydbl up)
{
  phydbl dens;
  phydbl cdf_up, cdf_lo;

  dens   = Dnorm(x,mean,sd);
  cdf_up = Pnorm(up,mean,sd);
  cdf_lo = Pnorm(lo,mean,sd);

  dens /= (cdf_up - cdf_lo);

  if(isnan(dens) || isinf(FABS(dens)))
    {
      PhyML_Printf("\n. mean=%f sd=%f lo=%f up=%f cdf_lo=%G CDF_up=%G",mean,sd,lo,up,cdf_lo,cdf_up);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  return dens;
}

/*********************************************************/

phydbl Dnorm_Multi(phydbl *x, phydbl *mu, phydbl *cov, int size, int _log)
{
  phydbl *xmmu,*invcov;
  phydbl *buff1,*buff2;
  int i;
  phydbl det,density;

  xmmu   = (phydbl *)mCalloc(size,sizeof(phydbl));
  invcov = (phydbl *)mCalloc(size*size,sizeof(phydbl));

  For(i,size) xmmu[i] = x[i] - mu[i];
  For(i,size*size) invcov[i] = cov[i];
  
  if(!Matinv(invcov,size,size))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");      
    }

  buff1 = Matrix_Mult(xmmu,invcov,1,size,size,size);
  buff2 = Matrix_Mult(buff1,xmmu,1,size,size,1);
  
  det = Matrix_Det(cov,size,NO);
  /* det_1D(cov,size,&det); */

  density = size * LOG2PI + LOG(det) + buff2[0];
  density /= -2.;

/*   density = (1./(POW(2.*PI,size/2.)*SQRT(FABS(det)))) * EXP(-0.5*buff2[0]); */

  free(xmmu);
  free(invcov);
  free(buff1);
  free(buff2);

  return (_log)?(density):(EXP(density));
}

/*********************************************************/

phydbl Dnorm_Multi_Given_InvCov_Det(phydbl *x, phydbl *mu, phydbl *invcov, phydbl log_det, int size, int _log)
{
  phydbl *xmmu;
  phydbl *buff1,*buff2;
  int i;
  phydbl density;

  xmmu = (phydbl *)mCalloc(size,sizeof(phydbl));

  For(i,size) xmmu[i] = x[i] - mu[i];
  
  buff1 = Matrix_Mult(xmmu,invcov,1,size,size,size);
  buff2 = Matrix_Mult(buff1,xmmu,1,size,size,1);

  density = size * LOG2PI + log_det + buff2[0];
  density /= -2.;

  free(xmmu);
  free(buff1);
  free(buff2);
  
  return (_log)?(density):(EXP(density));
}

/*********************************************************/

phydbl Pbinom(int N, int ni, phydbl p)
{
  return Bico(N,ni)*POW(p,ni)*POW(1-p,N-ni);
}

/*********************************************************/

phydbl Bivariate_Normal_Density(phydbl x, phydbl y, phydbl mux, phydbl muy, phydbl sdx, phydbl sdy, phydbl rho)
{
  phydbl cx, cy;
  phydbl pi;
  phydbl dens;
  phydbl rho2;

  pi = 3.141593;

  cx = x - mux;
  cy = y - muy;

  rho2 = rho*rho;

  dens = 1./(2*pi*sdx*sdy*SQRT(1.-rho2));
  dens *= EXP((-1./(2.*(1.-rho2)))*(cx*cx/(sdx*sdx)+cy*cy/(sdy*sdy)+2*rho*cx*cy/(sdx*sdy)));
	      
  return dens;
}

/*********************************************************/

phydbl Dgamma_Moments(phydbl x, phydbl mean, phydbl var)
{
  phydbl shape, scale;

  if(var < 1.E-20) 
    {
/*       var  = 1.E-20;  */
      PhyML_Printf("\n. var=%f mean=%f",var,mean);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");      
    }

  if(mean < 1.E-20) 
    { 
/*       mean = 1.E-20;  */
      PhyML_Printf("\n. var=%f mean=%f",var,mean);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }


  shape = mean * mean / var;
  scale = var / mean;
  

  return(Dgamma(x,shape,scale));
}

/*********************************************************/

phydbl Dgamma(phydbl x, phydbl shape, phydbl scale)
{
  phydbl v;

  if(x > INFINITY) 
    {
      PhyML_Printf("\n. WARNING: huge value of x -> x = %G",x);
      x = 1.E+10;
    }

  if(x < 1.E-10)
    {
      if(x < 0.0) return 0.0;
      else
	{
	  PhyML_Printf("\n. WARNING: small value of x -> x = %G",x);
	  x = 1.E-10;
	}
    }


  if(scale < 0.0 || shape < 0.0)
    {
      PhyML_Printf("\n. scale=%f shape=%f",scale,shape);
      Exit("\n");
    }


  v = (shape-1.) * LOG(x) - shape * LOG(scale) - x / scale - LnGamma(shape);


  if(v < 500.)
    {
      v = EXP(v);
    }
  else
    {
      PhyML_Printf("\n. WARNING v=%f x=%f shape=%f scale=%f",v,x,shape,scale);
      PhyML_Printf("\n. LOG(x) = %G LnGamma(shape)=%G",LOG(x),LnGamma(shape));
      Exit("\n");
    }

	 
  return v;
}

/*********************************************************/

phydbl Dexp(phydbl x, phydbl param)
{
  return param * EXP(-param * x);
}

/*********************************************************/
phydbl Dpois(phydbl x, phydbl param)
{
  phydbl v;

  if(x < 0) 
    {
      PhyML_Printf("\n. x = %f",x);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  v = x * LOG(param) - param - LnGamma(x+1);

  if(v < 500)
    {
      v = EXP(v);
    }
  else
    {
      PhyML_Printf("\n. WARNING v=%f x=%f param=%f",v,x,param);
      v = EXP(500);
    }
  
/*   PhyML_Printf("\n. Poi %f %f (x=%f param=%f)", */
/* 	 v, */
/* 	 POW(param,x) * EXP(-param) / EXP(LnGamma(x+1)), */
/* 	 x,param); */
/*   return POW(param,x) * EXP(-param) / EXP(LnGamma(x+1)); */
  
  return v;
}

/*********************************************************/


/*********************************************************/
/* CDFs */
/*********************************************************/

phydbl Pnorm(phydbl x, phydbl mean, phydbl sd)
{
/*   const phydbl b1 =  0.319381530; */
/*   const phydbl b2 = -0.356563782; */
/*   const phydbl b3 =  1.781477937; */
/*   const phydbl b4 = -1.821255978; */
/*   const phydbl b5 =  1.330274429; */
/*   const phydbl p  =  0.2316419; */
/*   const phydbl c  =  0.39894228; */
  
  x = (x-mean)/sd;
  
/*   if(x >= 0.0) */
/*     { */
/*       phydbl t = 1.0 / ( 1.0 + p * x ); */
/*       return (1.0 - c * EXP( -x * x / 2.0 ) * t * */
/* 	      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 )); */
/*     } */
/*   else */
/*     { */
/*       phydbl t = 1.0 / ( 1.0 - p * x ); */
/*       return ( c * EXP( -x * x / 2.0 ) * t * */
/* 	       ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 )); */
/*     } */

/* i_tail in {0,1,2} means: "lower", "upper", or "both" :
   if(lower) return  *cum := P[X <= x]
   if(upper) return *ccum := P[X >  x] = 1 - P[X <= x]
*/

/*   return Pnorm_Marsaglia(x); */
  return Pnorm_Ihaka_Derived_From_Cody(x);
}


/* G. Marsaglia. "Evaluating the Normal distribution". Journal of Statistical Software. 2004. Vol. 11. Issue 4. */
phydbl  Pnorm_Marsaglia(phydbl x)
{
  long double s=x,t=0,b=x,q=x*x,i=1;
  while(s!=t) s=(t=s)+(b*=q/(i+=2)); 
  return .5+s*exp(-.5*q-.91893853320467274178L);

}



/* Stolen from R source code */
#define SIXTEN 16

phydbl Pnorm_Ihaka_Derived_From_Cody(phydbl x)
{

    const static double a[5] = {
	2.2352520354606839287,
	161.02823106855587881,
	1067.6894854603709582,
	18154.981253343561249,
	0.065682337918207449113
    };
    const static double b[4] = {
	47.20258190468824187,
	976.09855173777669322,
	10260.932208618978205,
	45507.789335026729956
    };
    const static double c[9] = {
	0.39894151208813466764,
	8.8831497943883759412,
	93.506656132177855979,
	597.27027639480026226,
	2494.5375852903726711,
	6848.1904505362823326,
	11602.651437647350124,
	9842.7148383839780218,
	1.0765576773720192317e-8
    };
    const static double d[8] = {
	22.266688044328115691,
	235.38790178262499861,
	1519.377599407554805,
	6485.558298266760755,
	18615.571640885098091,
	34900.952721145977266,
	38912.003286093271411,
	19685.429676859990727
    };
    const static double p[6] = {
	0.21589853405795699,
	0.1274011611602473639,
	0.022235277870649807,
	0.001421619193227893466,
	2.9112874951168792e-5,
	0.02307344176494017303
    };
    const static double q[5] = {
	1.28426009614491121,
	0.468238212480865118,
	0.0659881378689285515,
	0.00378239633202758244,
	7.29751555083966205e-5
    };

    double xden, xnum, temp, del, eps, xsq, y;
    int i, lower, upper;
    double cum,ccum;
    int i_tail;
    
    i_tail = 0;
    cum = ccum = 0.0;

    if(isnan(x)) { cum = ccum = x; return (phydbl)cum; }

    /* Consider changing these : */
    eps = DBL_EPSILON * 0.5;

    /* i_tail in {0,1,2} =^= {lower, upper, both} */
    lower = i_tail != 1;
    upper = i_tail != 0;

    y = fabs(x);
    if (y <= 0.67448975) { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
	if (y > eps) {
	    xsq = x * x;
	    xnum = a[4] * xsq;
	    xden = xsq;
	    for (i = 0; i < 3; ++i) {
		xnum = (xnum + a[i]) * xsq;
		xden = (xden + b[i]) * xsq;
	    }
	} else xnum = xden = 0.0;

	temp = x * (xnum + a[3]) / (xden + b[3]);
	if(lower)  cum = 0.5 + temp;
	if(upper) ccum = 0.5 - temp;
	}    
    else if (y <= M_SQRT_32) {

	/* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= SQRT(32) ~= 5.657 */

	xnum = c[8] * y;
	xden = y;
	for (i = 0; i < 7; ++i) {
	    xnum = (xnum + c[i]) * y;
	    xden = (xden + d[i]) * y;
	}
	temp = (xnum + c[7]) / (xden + d[7]);

#define do_del(X)							\
	xsq = floor(X * SIXTEN) / SIXTEN;			\
	del = (X - xsq) * (X + xsq);					\
	cum = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp;		\
	ccum = 1.0 - cum;						\
	
#define swap_tail						\
	if (x > 0.) {/* swap  ccum <--> cum */			\
	    temp = cum; if(lower) cum = ccum; ccum = temp;	\
	}

	do_del(y);
	swap_tail;
    }

/* else	  |x| > SQRT(32) = 5.657 :
 * the next two case differentiations were really for lower=T, log=F
 * Particularly	 *not*	for  log_p !

 * Cody had (-37.5193 < x  &&  x < 8.2924) ; R originally had y < 50
 *
 * Note that we do want symmetry(0), lower/upper -> hence use y
 */
    else if((lower && -37.5193 < x  &&  x < 8.2924) || (upper && -8.2924  < x  &&  x < 37.5193)) 
      {
	/* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
	xsq = 1.0 / (x * x);
	xnum = p[5] * xsq;
	xden = xsq;
	for (i = 0; i < 4; ++i) {
	    xnum = (xnum + p[i]) * xsq;
	    xden = (xden + q[i]) * xsq;
	}
	temp = xsq * (xnum + p[4]) / (xden + q[4]);
	temp = (M_1_SQRT_2_PI - temp) / y;

	do_del(x);
	swap_tail;
      }
    else 
      { /* no log_p , large x such that probs are 0 or 1 */
	if(x > 0) {	cum = 1.; ccum = 0.;	}
	else {	        cum = 0.; ccum = 1.;	}
      }

    return (phydbl)cum;


}

/*********************************************************/


phydbl Pgamma(phydbl x, phydbl shape, phydbl scale)
{
  return IncompleteGamma(x/scale,shape,LnGamma(shape));
}

/*********************************************************/

phydbl Ppois(phydbl x, phydbl param)
{
  /* Press et al. (1990) approximation of the CDF for the Poisson distribution */
  if(param < SMALL || x < 0.0) 
    {
      PhyML_Printf("\n. param = %G x=%G",param,x);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  return IncompleteGamma(x,param,LnGamma(param));
}

/*********************************************************/

/*********************************************************/
/* Inverse CDFs */
/*********************************************************/

phydbl PointChi2 (phydbl prob, phydbl v)
{
/* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
   returns -1 if in error.   0.000002<prob<0.999998
   RATNEST FORTRAN by
       Best DJ & Roberts DE (1975) The percentage points of the
       Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
   Converted into C by Ziheng Yang, Oct. 1993.
*/
   double aa=.6931471805, p=prob, g;
   double xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6;
   double e=.5e-6;
   
   if (p<.000002 || p>.999998 || v<=0) return ((phydbl)-1);

   g = (double)LnGamma(v/2);
   xx=v/2;   c=xx-1;
   if (v >= -1.24*log(p)) goto l1;

   ch=pow((p*xx*exp(g+xx*aa)), 1/xx);
   if (ch-e<0) return (ch);
   goto l4;
l1:
   if (v>.32) goto l3;
   ch=0.4;   a=log(1-p);
l2:
   q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch));
   t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
   ch-=(1-exp(a+g+.5*ch+c*aa)*p2/p1)/t;
   if (fabs(q/ch-1)-.01 <= 0) goto l4;
   else                       goto l2;

l3:
   x=(double)PointNormal (p);
   p1=0.222222/v;   ch=v*pow((x*sqrt(p1)+1-p1), 3.0);
   if (ch>2.2*v+6)  ch=-2*(log(1-p)-c*log(.5*ch)+g);
l4:
   q=ch;   p1=.5*ch;
   if ((t=(double)IncompleteGamma (p1, xx, g))<0) {
      PhyML_Printf ("\nerr IncompleteGamma");
      return ((phydbl)-1.);
   }
   p2=p-t;
   t=p2*exp(xx*aa+g+p1-c*log(ch));
   b=t/ch;  a=0.5*t-b*c;

   s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
   s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520;
   s3=(210+a*(462+a*(707+932*a)))/2520;
   s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
   s5=(84+264*a+c*(175+606*a))/2520;
   s6=(120+c*(346+127*c))/5040;
   ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
   if (FABS(q/ch-1) > e) goto l4;

   return (phydbl)(ch);
}

/*********************************************************/

phydbl PointNormal (phydbl p)
{
/* returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
   returns (-9999) if in error
   Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
   Applied Statistics 22: 96-97 (AS70)

   Newer methods:
     Wichura MJ (1988) Algorithm AS 241: the percentage points of the
       normal distribution.  37: 477-484.
     Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage
       points of the normal distribution.  26: 118-121.

*/
/*    phydbl a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245; */
/*    phydbl a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495; */
/*    phydbl b2=.531103462366, b3=.103537752850, b4=.0038560700634; */
/*    phydbl y, z=0, p=prob, p1; */

/*    p1 = (p<0.5 ? p : 1-p); */
/*    if (p1<1e-20) return (-INFINITY); */
/* /\*    if (p1<1e-20) return (-999.); *\/ */

/*    y = sqrt ((phydbl)LOG(1/(p1*p1))); */
/*    z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0); */
/*    return (p<0.5 ? -z : z); */

  static double zero = 0.0, one = 1.0, half = 0.5;
  static double split1 = 0.425, split2 = 5.0;
  static double const1 = 0.180625, const2 = 1.6;
  
  /* coefficients for p close to 0.5 */
  static double a[8] = {
    3.3871328727963666080e0,
    1.3314166789178437745e+2,
    1.9715909503065514427e+3,
    1.3731693765509461125e+4,
    4.5921953931549871457e+4,
    6.7265770927008700853e+4,
    3.3430575583588128105e+4,
    2.5090809287301226727e+3
  };
  static double b[8] = { 
    0.0,
    4.2313330701600911252e+1,
    6.8718700749205790830e+2,
    5.3941960214247511077e+3,
    2.1213794301586595867e+4,
    3.9307895800092710610e+4,
    2.8729085735721942674e+4,
    5.2264952788528545610e+3
  };
  
  /* hash sum ab    55.8831928806149014439 */
  /* coefficients for p not close to 0, 0.5 or 1. */
  static double c[8] = {
    1.42343711074968357734e0,
    4.63033784615654529590e0,
    5.76949722146069140550e0,
    3.64784832476320460504e0,
    1.27045825245236838258e0,
    2.41780725177450611770e-1,
    2.27238449892691845833e-2,
    7.74545014278341407640e-4
  };
  static double d[8] = { 
    0.0,
    2.05319162663775882187e0,
    1.67638483018380384940e0,
    6.89767334985100004550e-1,
    1.48103976427480074590e-1,
    1.51986665636164571966e-2,
    5.47593808499534494600e-4,
    1.05075007164441684324e-9
  };
  
  /* hash sum cd    49.33206503301610289036 */
  /* coefficients for p near 0 or 1. */
  static double e[8] = {
    6.65790464350110377720e0,
    5.46378491116411436990e0,
    1.78482653991729133580e0,
    2.96560571828504891230e-1,
    2.65321895265761230930e-2,
    1.24266094738807843860e-3,
    2.71155556874348757815e-5,
    2.01033439929228813265e-7
  };
  static double f[8] = { 
    0.0,
    5.99832206555887937690e-1,
    1.36929880922735805310e-1,
    1.48753612908506148525e-2,
    7.86869131145613259100e-4,
    1.84631831751005468180e-5,
    1.42151175831644588870e-7,
    2.04426310338993978564e-15
  };
  
  /* hash sum ef    47.52583317549289671629 */
  double q, r, ret;
  
  q = p - half;
  if (fabs(q) <= split1) {
    r = const1 - q * q;
    ret = q * (((((((a[7] * r + a[6]) * r + a[5]) * r + a[4]) * r + a[3])
		 * r + a[2]) * r + a[1]) * r + a[0]) /
      (((((((b[7] * r + b[6]) * r + b[5]) * r + b[4]) * r + b[3])
	 * r + b[2]) * r + b[1]) * r + one);
    
    return (phydbl)ret;
  }
  /* else */
  
  if (q < zero)
    r = p;
  else
    r = one - p;
  
  if (r <= zero)
    return (phydbl)zero;
  
  r = sqrt(-log(r));
  if (r <= split2) {
    r -= const2;
    ret = (((((((c[7] * r + c[6]) * r + c[5]) * r + c[4]) * r + c[3])
	     * r + c[2]) * r + c[1]) * r + c[0]) /
      (((((((d[7] * r + d[6]) * r + d[5]) * r + d[4]) * r + d[3])
	 * r + d[2]) * r + d[1]) * r + one);
  }
  else {
    r -= split2;
    ret = (((((((e[7] * r + e[6]) * r + e[5]) * r + e[4]) * r + e[3])
	     * r + e[2]) * r + e[1]) * r + e[0]) /
      (((((((f[7] * r + f[6]) * r + f[5]) * r + f[4]) * r + f[3])
	 * r + f[2]) * r + f[1]) * r + one);
  }
  
  if (q < zero)
    ret = -ret;
  
  return (phydbl)ret;
}

/*********************************************************/
/* MISCs */
/*********************************************************/

/*********************************************************/

phydbl Bico(int n, int k)
{
  return FLOOR(0.5+EXP(Factln(n)-Factln(k)-Factln(n-k)));
}


/*********************************************************/

phydbl Factln(int n)
{
  static phydbl a[101];
  
  if (n < 0)    { Warn_And_Exit("\n. Err: negative factorial in routine FACTLN"); }
  if (n <= 1)     return 0.0;
  if (n <= 100)   return (a[n]>SMALL) ? a[n] : (a[n]=Gammln(n+1.0));
  else return     Gammln(n+1.0);
}

/*********************************************************/

phydbl Gammln(phydbl xx)
{
  double x,tmp,ser;
  static double cof[6]={76.18009173,-86.50532033,24.01409822,
			-1.231739516,0.120858003e-2,-0.536382e-5};
  int j;
  
  x=xx-1.0;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.0;
  for (j=0;j<=5;j++) 
    {
      x += 1.0;
      ser += cof[j]/x;
    }
  return (phydbl)(-tmp+log(2.50662827465*ser));
}

/*********************************************************/

/* void Plim_Binom(phydbl pH0, int N, phydbl *pinf, phydbl *psup) */
/* { */
/*   *pinf = pH0 - 1.64*SQRT(pH0*(1-pH0)/(phydbl)N); */
/*   if(*pinf < 0) *pinf = .0; */
/*   *psup = pH0 + 1.64*SQRT(pH0*(1-pH0)/(phydbl)N); */
/* } */

/*********************************************************/

phydbl LnGamma (phydbl alpha)
{
/* returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.
   Stirling's formula is used for the central polynomial part of the procedure.
   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
   double x=alpha, f=0, z;
   if (x<7) {
      f=1;  z=x-1;
      while (++z<7)  f*=z;
      x=z;   f=-log(f);
   }
   z = 1/(x*x);
   return (phydbl)(f + (x-0.5)*log(x) - x + .918938533204673
		   + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
		      +.083333333333333)/x);
}

/*********************************************************/

phydbl IncompleteGamma(phydbl x, phydbl alpha, phydbl ln_gamma_alpha)
{
/* returns the incomplete gamma ratio I(x,alpha) where x is the upper
	   limit of the integration and alpha is the shape parameter.
   returns (-1) if in error
   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
   (1) series expansion     if (alpha>x || x<=1)
   (2) continued fraction   otherwise
   RATNEST FORTRAN by
   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
   19: 285-287 (AS32)
*/
   int i;
   double p=alpha, g=ln_gamma_alpha;
   double accurate=1e-8, overflow=1e30;
   double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6];

   if (fabs(x) < SMALL) return ((phydbl).0);
   if (x<0 || p<=0)        return ((phydbl)-1);

   factor=exp(p*log(x)-x-g);
   if (x>1 && x>=p) goto l30;
   /* (1) series expansion */
   gin=1;  term=1;  rn=p;
 l20:
   rn++;
   term*=x/rn;   gin+=term;

   if (term > accurate) goto l20;
   gin*=factor/p;
   goto l50;
 l30:
   /* (2) continued fraction */
   a=1-p;   b=a+x+1;  term=0;
   pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
   gin=pn[2]/pn[3];
 l32:
   a++;  b+=2;  term++;   an=a*term;
   for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
   if (fabs(pn[5]) < .0) goto l35;
   rn=pn[4]/pn[5];   dif=fabs(gin-rn);
   if (dif>accurate) goto l34;
   if (dif<=accurate*rn) goto l42;
 l34:
   gin=rn;
 l35:
   for (i=0; i<4; i++) pn[i]=pn[i+2];
   if (fabs(pn[4]) < overflow) goto l32;
   for (i=0; i<4; i++) pn[i]/=overflow;
   goto l32;
 l42:
   gin=1-factor*gin;

 l50:
   return (phydbl)(gin);
}


/*********************************************************/

int DiscreteGamma (phydbl freqK[], phydbl rK[], phydbl alfa, phydbl beta, int K, int median)
{
  /* discretization of gamma distribution with equal proportions in each
  category
  */
  
  int i;
  phydbl gap05=1.0/(2.0*K), t, factor=alfa/beta*K, lnga1;
  
  if(K==1)
  {
    freqK[0] = 1.0;
    rK[0] = 1.0;
    return 0;
  }
  
  if (median) 
  {
    for (i=0; i<K; i++)     rK[i]=PointGamma((i*2.0+1)*gap05, alfa, beta);
    for (i=0,t=0; i<K; i++) t+=rK[i];
    for (i=0; i<K; i++)     rK[i]*=factor/t;
  }
  else 
  {
    lnga1=LnGamma(alfa+1);
    for (i=0; i<K-1; i++) freqK[i]=PointGamma((i+1.0)/K, alfa, beta);
    for (i=0; i<K-1; i++) freqK[i]=IncompleteGamma(freqK[i]*beta, alfa+1, lnga1);
    rK[0] = freqK[0]*factor;
    rK[K-1] = (1-freqK[K-2])*factor;
    for (i=1; i<K-1; i++) rK[i] = (freqK[i]-freqK[i-1])*factor;
  }
  for (i=0; i<K; i++) freqK[i]=1.0/K;
  return (0);
}

/*********************************************************/

/* Return LOG(n!) */

phydbl LnFact(int n)
{
  int i;
  phydbl res;

  res = 0;
  for(i=2;i<=n;i++) res += LOG(i);
  
  return(res);
}

/*********************************************************/

int Choose(int n, int k)
{
  phydbl accum;
  int i;

  if (k > n) return(0);
  if (k > n/2) k = n-k;
  if(!k) return(1);

  accum = 1.;
  for(i=1;i<k+1;i++) accum = accum * (n-k+i) / i;

  return((int)accum);
}

/*********************************************************/


phydbl *Covariance_Matrix(t_tree *tree)
{
  phydbl *cov, *mean,var_min;
  int *ori_wght,*site_num;
  int dim,i,j,replicate,n_site,position,sample_size;

  sample_size = 100;
  dim = 2*tree->n_otu-3;

  cov      = (phydbl *)mCalloc(dim*dim,sizeof(phydbl));
  mean     = (phydbl *)mCalloc(    dim,sizeof(phydbl));
  ori_wght = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
  site_num = (int *)mCalloc(tree->data->init_len,sizeof(int));
  
  var_min = 1./POW(tree->data->init_len,2);

  For(i,tree->data->crunch_len) ori_wght[i] = tree->data->wght[i];

  n_site = 0;
  For(i,tree->data->crunch_len) For(j,tree->data->wght[i])
    {
      site_num[n_site] = i;
      n_site++;
    }

  		  
  tree->mod->s_opt->print = 0;
  For(replicate,sample_size)
    {
      For(i,2*tree->n_otu-3) tree->t_edges[i]->l = .1;

      For(i,tree->data->crunch_len) tree->data->wght[i] = 0;

      For(i,tree->data->init_len)
	{
	  position = Rand_Int(0,(int)(tree->data->init_len-1.0));
	  tree->data->wght[site_num[position]] += 1;
	}

      Round_Optimize(tree,tree->data,ROUND_MAX);
      
      For(i,2*tree->n_otu-3) For(j,2*tree->n_otu-3) cov[i*dim+j] += LOG(tree->t_edges[i]->l) * LOG(tree->t_edges[j]->l);  
      For(i,2*tree->n_otu-3) mean[i] += LOG(tree->t_edges[i]->l);

      PhyML_Printf("[%3d/%3d]",replicate,sample_size); fflush(NULL);
/*       PhyML_Printf("\n. %3d %12f %12f %12f ", */
/* 	     replicate, */
/*  	     cov[1*dim+1]/(replicate+1)-mean[1]*mean[1]/POW(replicate+1,2), */
/* 	     tree->t_edges[1]->l, */
/* 	     mean[1]/(replicate+1)); */
   }

  For(i,2*tree->n_otu-3) mean[i] /= (phydbl)sample_size;
  For(i,2*tree->n_otu-3) For(j,2*tree->n_otu-3) cov[i*dim+j] /= (phydbl)sample_size;
  For(i,2*tree->n_otu-3) For(j,2*tree->n_otu-3) cov[i*dim+j] -= mean[i]*mean[j];
/*   For(i,2*tree->n_otu-3) if(cov[i*dim+i] < var_min) cov[i*dim+i] = var_min; */
  

/*   PhyML_Printf("\n"); */
/*   For(i,2*tree->n_otu-3) PhyML_Printf("%f %f\n",mean[i],tree->t_edges[i]->l); */
/*   PhyML_Printf("\n"); */
/*   PhyML_Printf("\n"); */
/*   For(i,2*tree->n_otu-3) */
/*     { */
/*       For(j,2*tree->n_otu-3) */
/* 	{ */
/* 	  PhyML_Printf("%G\n",cov[i*dim+j]); */
/* 	} */
/*       PhyML_Printf("\n"); */
/*     } */

  For(i,tree->data->crunch_len) tree->data->wght[i] = ori_wght[i];

  free(mean);
  free(ori_wght);
  free(site_num);

  return cov;
}

/*********************************************************/
/* Work out the Hessian for the likelihood function. Only branch lengths are considered as variable.
   This function is very much inspired from Jeff Thorne's 'hessian' function in his program 'estbranches'. */
phydbl *Hessian(t_tree *tree)
{
  phydbl *hessian;
  phydbl *plus_plus, *minus_minus, *plus_zero, *minus_zero, *plus_minus, zero_zero;
  phydbl *ori_bl,*inc,*buff;
  int *ok_edges,*is_ok;
  int dim;
  int n_ok_edges;
  int i,j;
  phydbl eps;
  phydbl lk;
  phydbl lnL,lnL1,lnL2;

  dim = 2*tree->n_otu-3;
  eps = 0.01;

  hessian     = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  ori_bl      = (phydbl *)mCalloc((int)dim,sizeof(phydbl));
  plus_plus   = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  minus_minus = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  plus_minus  = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  plus_zero   = (phydbl *)mCalloc((int)dim    ,sizeof(phydbl));
  minus_zero  = (phydbl *)mCalloc((int)dim    ,sizeof(phydbl));
  inc         = (phydbl *)mCalloc((int)dim    ,sizeof(phydbl));
  buff        = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  ok_edges    = (int *)mCalloc((int)dim,sizeof(int));
  is_ok       = (int *)mCalloc((int)dim,sizeof(int));

  lnL = lnL1 = lnL2 = UNLIKELY;

  tree->both_sides = 1;
  Lk(tree);

  For(i,dim) ori_bl[i] = tree->t_edges[i]->l;


  n_ok_edges = 0;
  For(i,dim) 
    {
      if(tree->t_edges[i]->l*(1.-eps) > BL_MIN)
	{	  
	  inc[i] = eps * tree->t_edges[i]->l;
	  ok_edges[n_ok_edges] = i;
	  n_ok_edges++;
	  is_ok[i] = 1;
	}
      else
	{
	  inc[i] = -1.0;
	  is_ok[i] = 0;
	}
    }

  /* zero zero */  
  zero_zero = tree->c_lnL;

  /* plus zero */  
  For(i,dim) 
    {
      if(is_ok[i])
	{
	  tree->t_edges[i]->l += inc[i];
	  lk = Lk_At_Given_Edge(tree->t_edges[i],tree);
	  plus_zero[i] = lk;
	  tree->t_edges[i]->l = ori_bl[i];
	}
    }


  /* minus zero */  
  For(i,dim) 
    {
      if(is_ok[i])
	{
	  tree->t_edges[i]->l -= inc[i];
	  lk = Lk_At_Given_Edge(tree->t_edges[i],tree);
	  minus_zero[i] = lk;
	  tree->t_edges[i]->l = ori_bl[i];
	}
    }

  For(i,dim) Update_PMat_At_Given_Edge(tree->t_edges[i],tree);

  /* plus plus  */  
  For(i,dim)
    {
      if(is_ok[i])
	{
	  tree->t_edges[i]->l += inc[i];
	  Update_PMat_At_Given_Edge(tree->t_edges[i],tree);

	  For(j,3)
	    if((!tree->t_edges[i]->left->tax) && (tree->t_edges[i]->left->v[j] != tree->t_edges[i]->rght))
	      Recurr_Hessian(tree->t_edges[i]->left,tree->t_edges[i]->left->v[j],1,inc,plus_plus+i*dim,is_ok,tree);
	  
	  For(j,3)
	    if((!tree->t_edges[i]->rght->tax) && (tree->t_edges[i]->rght->v[j] != tree->t_edges[i]->left))
	      Recurr_Hessian(tree->t_edges[i]->rght,tree->t_edges[i]->rght->v[j],1,inc,plus_plus+i*dim,is_ok,tree);
		      
	  tree->t_edges[i]->l = ori_bl[i];
	  Lk(tree);
	}
    }

  /* plus minus */  
  For(i,dim)
    {
      if(is_ok[i])
	{
	  tree->t_edges[i]->l += inc[i];
	  Update_PMat_At_Given_Edge(tree->t_edges[i],tree);
	  
	  For(j,3)
	    if((!tree->t_edges[i]->left->tax) && (tree->t_edges[i]->left->v[j] != tree->t_edges[i]->rght))
	      Recurr_Hessian(tree->t_edges[i]->left,tree->t_edges[i]->left->v[j],-1,inc,plus_minus+i*dim,is_ok,tree);
	  
	  For(j,3)
	    if((!tree->t_edges[i]->rght->tax) && (tree->t_edges[i]->rght->v[j] != tree->t_edges[i]->left))
	      Recurr_Hessian(tree->t_edges[i]->rght,tree->t_edges[i]->rght->v[j],-1,inc,plus_minus+i*dim,is_ok,tree);
	  
	  tree->t_edges[i]->l = ori_bl[i];
	  Lk(tree);
	}
    }

  /* minus minus */  
  For(i,dim)
    {
      if(is_ok[i])
	{
	  tree->t_edges[i]->l -= inc[i];

	  if(tree->t_edges[i]->l < BL_MIN)
	    {
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Exit("\n");
	    }
	  
	  Update_PMat_At_Given_Edge(tree->t_edges[i],tree);
	  
	  For(j,3)
	    if((!tree->t_edges[i]->left->tax) && (tree->t_edges[i]->left->v[j] != tree->t_edges[i]->rght))
	      Recurr_Hessian(tree->t_edges[i]->left,tree->t_edges[i]->left->v[j],-1,inc,minus_minus+i*dim,is_ok,tree);
	  
	  For(j,3)
	    if((!tree->t_edges[i]->rght->tax) && (tree->t_edges[i]->rght->v[j] != tree->t_edges[i]->left))
	      Recurr_Hessian(tree->t_edges[i]->rght,tree->t_edges[i]->rght->v[j],-1,inc,minus_minus+i*dim,is_ok,tree);
	  
	  tree->t_edges[i]->l = ori_bl[i];
	  Lk(tree);
	}
    }

  
  For(i,dim)
    {
      if(is_ok[i])
	{
	  hessian[i*dim+i] = (plus_zero[i]-2*zero_zero+minus_zero[i])/(POW(inc[i],2));

	  for(j=i+1;j<dim;j++)
	    {
	      if(is_ok[j])
		{
		  hessian[i*dim+j] = 
		    (plus_plus[i*dim+j]-plus_minus[i*dim+j]-plus_minus[j*dim+i]+minus_minus[i*dim+j])/
		    (4*inc[i]*inc[j]);
		  hessian[j*dim+i] = hessian[i*dim+j];
		}
	    }
	}
    }
        
  For(i,n_ok_edges)
    {
      For(j,n_ok_edges)
	{
	  buff[i*n_ok_edges+j] = -1.0*hessian[ok_edges[i]*dim+ok_edges[j]];
	}
    }

  if(!Matinv(buff,n_ok_edges,n_ok_edges))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");      
    }

  For(i,n_ok_edges)
    {
      For(j,n_ok_edges)
	{
	  hessian[ok_edges[i]*dim+ok_edges[j]] = buff[i*n_ok_edges+j];
	}
    }

  /* Approximate variance for very short branches */
  For(i,dim)
    if(inc[i] < 0.0)
      {
	lnL  = tree->c_lnL;
	tree->t_edges[i]->l += eps;
	lnL1 = Lk_At_Given_Edge(tree->t_edges[i],tree);
	tree->t_edges[i]->l += eps;
	lnL2 = Lk_At_Given_Edge(tree->t_edges[i],tree);

	hessian[i*dim+i] = (lnL2 - 2*lnL1 + lnL) / POW(eps,2);
	hessian[i*dim+i] = -1.0 / hessian[i*dim+i];
      }

  For(i,dim)
    if(hessian[i*dim+i] < MIN_VAR_BL)
      {
	lnL  = tree->c_lnL;
	tree->t_edges[i]->l += eps;
	lnL1 = Lk_At_Given_Edge(tree->t_edges[i],tree);
	tree->t_edges[i]->l += eps;
	lnL2 = Lk_At_Given_Edge(tree->t_edges[i],tree);

	hessian[i*dim+i] = (lnL2 - 2*lnL1 + lnL) / POW(eps,2);
	hessian[i*dim+i] = -1.0 / hessian[i*dim+i];
      }

  For(i,dim)
    if(hessian[i*dim+i] < 0.0)
      {
	PhyML_Printf("\n. l=%G var=%G",tree->t_edges[i]->l,hessian[i*dim+i]);
	PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	Exit("\n");
      }

  For(i,dim)
    if(hessian[i*dim+i] < MIN_VAR_BL)
      {
	PhyML_Printf("\n. l=%G var=%G",tree->t_edges[i]->l,hessian[i*dim+i]);
	hessian[i*dim+i] = MIN_VAR_BL;
	PhyML_Printf("\n. Numerical precision issues may alter this analysis...");
      }


  if(!Matinv(hessian,dim,dim))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  For(i,dim*dim) hessian[i] = -1.0*hessian[i];

  For(i,dim)
    {
      For(j,dim)
	{
	  if(FABS(hessian[i*dim+j]-hessian[j*dim+i]) > 1.E-3)
	    {
	      PhyML_Printf("\n. Hessian not symmetrical.");
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Exit("\n");
	    }
	  hessian[i*dim+j] = (hessian[i*dim+j] + hessian[j*dim+i]) / 2.; 
	  hessian[j*dim+i] = hessian[i*dim+j];  
	}
    }
  

/*   For(i,dim) */
/*     { */
/*       PhyML_Printf("[%f] ",tree->t_edges[i]->l); */
/*       For(j,dim) */
/* 	{ */
/* 	  PhyML_Printf("%12lf ",hessian[i*dim+j]); */
/* 	} */
/*       PhyML_Printf("\n"); */
/*     } */

/*   Matinv(hessian,dim,dim); */

/*   PhyML_Printf("\n"); */

/*   For(i,dim) */
/*     { */
/*       PhyML_Printf("[%f] ",tree->t_edges[i]->l); */
/*       For(j,dim) */
/* 	{ */
/* 	  PhyML_Printf("%12lf ",-hessian[i*dim+j]); */
/* 	} */
/*       PhyML_Printf("\n"); */
/*     } */
/*   Exit("\n"); */

  free(ori_bl);
  free(plus_plus);
  free(minus_minus);
  free(plus_zero);
  free(minus_zero);
  free(plus_minus);
  free(inc);
  free(buff);
  free(ok_edges);
  free(is_ok);

  return hessian;

}

/*********************************************************/

void Recurr_Hessian(t_node *a, t_node *d, int plus_minus, phydbl *inc, phydbl *res, int *is_ok, t_tree *tree)
{
  int i;
  phydbl ori_l;

  For(i,3)
    if(a->v[i] == d)
      {
	Update_P_Lk(tree,a->b[i],a);

	ori_l = a->b[i]->l;
	if(is_ok[a->b[i]->num])
	  {
	    if(plus_minus > 0) a->b[i]->l += inc[a->b[i]->num];
	    else               a->b[i]->l -= inc[a->b[i]->num];
	    res[a->b[i]->num] = Lk_At_Given_Edge(a->b[i],tree);
	    a->b[i]->l = ori_l;
	    Update_PMat_At_Given_Edge(a->b[i],tree);
	  }
	break;
      }

  if(d->tax) return;
  else 
    For(i,3) 
      if(d->v[i] != a) 
	Recurr_Hessian(d,d->v[i],plus_minus,inc,res,is_ok,tree);
}

/*********************************************************/

/* Work out the Hessian for the likelihood function. Only LOGARITHM of branch lengths are considered as variable.
   This function is very much inspired from Jeff Thorne's 'hessian' function in his program 'estbranches'. */
phydbl *Hessian_Log(t_tree *tree)
{
  phydbl *hessian;
  phydbl *plus_plus, *minus_minus, *plus_zero, *minus_zero, *plus_minus, *zero_zero;
  phydbl *ori_bl,*inc,*buff;
  int *ok_edges,*is_ok;
  int dim;
  int n_ok_edges;
  int i,j;
  phydbl eps;
  phydbl lk;

  dim = 2*tree->n_otu-3;
  eps = 1.E-4;

  hessian     = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  ori_bl      = (phydbl *)mCalloc((int)dim,    sizeof(phydbl));
  plus_plus   = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  minus_minus = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  plus_minus  = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  plus_zero   = (phydbl *)mCalloc((int)dim    ,sizeof(phydbl));
  minus_zero  = (phydbl *)mCalloc((int)dim    ,sizeof(phydbl));
  zero_zero   = (phydbl *)mCalloc((int)dim    ,sizeof(phydbl));
  inc         = (phydbl *)mCalloc((int)dim    ,sizeof(phydbl));
  buff        = (phydbl *)mCalloc((int)dim*dim,sizeof(phydbl));
  ok_edges    = (int *)mCalloc((int)dim,       sizeof(int));
  is_ok       = (int *)mCalloc((int)dim,       sizeof(int));
  
  tree->both_sides = 1;
  Lk(tree);

  For(i,dim) ori_bl[i] = tree->t_edges[i]->l;

  n_ok_edges = 0;
  For(i,dim) 
    {
      if(tree->t_edges[i]->l > 3.0/(phydbl)tree->data->init_len)
	{
	  inc[i] = eps * tree->t_edges[i]->l;
	  ok_edges[n_ok_edges] = i;
	  n_ok_edges++;
	  is_ok[i] = 1;
	}
      else is_ok[i] = 0;
    }

  /* zero zero */  
  lk = Log_Det(is_ok,tree);
  For(i,dim) if(is_ok[i]) zero_zero[i] = tree->c_lnL+lk;

  /* plus zero */  
  For(i,dim) 
    {
      if(is_ok[i])
	{
	  tree->t_edges[i]->l += inc[i];
	  lk = Lk_At_Given_Edge(tree->t_edges[i],tree);
	  plus_zero[i] = lk+Log_Det(is_ok,tree);
	  tree->t_edges[i]->l = ori_bl[i];
	}
    }


  /* minus zero */
  For(i,dim) 
    {
      if(is_ok[i])
	{
	  tree->t_edges[i]->l -= inc[i];
	  lk = Lk_At_Given_Edge(tree->t_edges[i],tree);
	  minus_zero[i] = lk+Log_Det(is_ok,tree);
	  tree->t_edges[i]->l = ori_bl[i];
	}
    }

  For(i,dim) Update_PMat_At_Given_Edge(tree->t_edges[i],tree);

  /* plus plus  */  
  For(i,dim)
    {
      if(is_ok[i])
	{
	  tree->t_edges[i]->l += inc[i];
	  Update_PMat_At_Given_Edge(tree->t_edges[i],tree);

	  For(j,3)
	    if((!tree->t_edges[i]->left->tax) && (tree->t_edges[i]->left->v[j] != tree->t_edges[i]->rght))
	      Recurr_Hessian_Log(tree->t_edges[i]->left,tree->t_edges[i]->left->v[j],1,inc,plus_plus+i*dim,is_ok,tree);
	  
	  For(j,3)
	    if((!tree->t_edges[i]->rght->tax) && (tree->t_edges[i]->rght->v[j] != tree->t_edges[i]->left))
	      Recurr_Hessian_Log(tree->t_edges[i]->rght,tree->t_edges[i]->rght->v[j],1,inc,plus_plus+i*dim,is_ok,tree);

/* 	  For(j,dim)  */
/* 	    if(j != i) */
/* 	      { */
/* 		if(inc[j] > 0.0) */
/* 		  { */
/* 		    tree->t_edges[j]->l += inc[j]; */
/* 		    Lk(tree); */
/* 		    plus_plus[i*dim+j]=tree->c_lnL; */
/* 		    tree->t_edges[j]->l = ori_bl[j]; */
/* 		  } */
/* 	      } */
		      
	  tree->t_edges[i]->l = ori_bl[i];
	  Lk(tree);
	}
    }

  /* plus minus */  
  For(i,dim)
    {
      if(is_ok[i])
	{
	  tree->t_edges[i]->l += inc[i];
	  Update_PMat_At_Given_Edge(tree->t_edges[i],tree);
	  
	  For(j,3)
	    if((!tree->t_edges[i]->left->tax) && (tree->t_edges[i]->left->v[j] != tree->t_edges[i]->rght))
	      Recurr_Hessian_Log(tree->t_edges[i]->left,tree->t_edges[i]->left->v[j],-1,inc,plus_minus+i*dim,is_ok,tree);
	  
	  For(j,3)
	    if((!tree->t_edges[i]->rght->tax) && (tree->t_edges[i]->rght->v[j] != tree->t_edges[i]->left))
	      Recurr_Hessian_Log(tree->t_edges[i]->rght,tree->t_edges[i]->rght->v[j],-1,inc,plus_minus+i*dim,is_ok,tree);
	  
/* 	  For(j,dim)  */
/* 	    if(j != i) */
/* 	      { */
/* 		if(inc[j] > 0.0) */
/* 		  { */
/* 		    tree->t_edges[j]->l -= inc[j]; */
/* 		    Lk(tree); */
/* 		    plus_minus[i*dim+j] = tree->c_lnL; */
/* 		    tree->t_edges[j]->l = ori_bl[j]; */
/* 		  } */
/* 	      } */

	  tree->t_edges[i]->l = ori_bl[i];
	  Lk(tree);
	}
    }


  /* minus minus */  
  For(i,dim)
    {
      if(is_ok[i])
	{
	  tree->t_edges[i]->l -= inc[i];
	  
	  Update_PMat_At_Given_Edge(tree->t_edges[i],tree);
	  
	  For(j,3)
	    if((!tree->t_edges[i]->left->tax) && (tree->t_edges[i]->left->v[j] != tree->t_edges[i]->rght))
	      Recurr_Hessian_Log(tree->t_edges[i]->left,tree->t_edges[i]->left->v[j],-1,inc,minus_minus+i*dim,is_ok,tree);
	  
	  For(j,3)
	    if((!tree->t_edges[i]->rght->tax) && (tree->t_edges[i]->rght->v[j] != tree->t_edges[i]->left))
	      Recurr_Hessian_Log(tree->t_edges[i]->rght,tree->t_edges[i]->rght->v[j],-1,inc,minus_minus+i*dim,is_ok,tree);
	  
/* 	  For(j,dim)  */
/* 	    if(j != i) */
/* 	      { */
/* 		if(inc[j] > 0.0) */
/* 		  { */
/* 		    tree->t_edges[j]->l -= inc[j]; */
/* 		    Lk(tree); */
/* 		    minus_minus[i*dim+j] = tree->c_lnL; */
/* 		    tree->t_edges[j]->l = ori_bl[j]; */
/* 		  } */
/* 	      } */

	  tree->t_edges[i]->l = ori_bl[i];
	  Lk(tree);
	}
    }

/*   For(i,dim) if(is_ok[i]) inc[i] = POW(tree->t_edges[i]->l+inc[i],2)-POW(tree->t_edges[i]->l,2); */
  For(i,dim) if(is_ok[i]) inc[i] = LOG(tree->t_edges[i]->l+inc[i])-LOG(tree->t_edges[i]->l);
/*   For(i,dim) inc[i] = 2.*inc[i]; */
/*   For(i,dim) if(is_ok[i]) inc[i] = SQRT(tree->t_edges[i]->l+inc[i])-SQRT(tree->t_edges[i]->l); */
  
  For(i,dim)
    {
      if(is_ok[i])
	{
	  hessian[i*dim+i] = (plus_zero[i]-2*zero_zero[i]+minus_zero[i])/(POW(inc[i],2));

	  for(j=i+1;j<dim;j++)
	    {
	      if(is_ok[j])
		{
		  hessian[i*dim+j] = 
		    (plus_plus[i*dim+j]-plus_minus[i*dim+j]-plus_minus[j*dim+i]+minus_minus[i*dim+j])/
		    (4*inc[i]*inc[i]);
		  hessian[j*dim+i] = hessian[i*dim+j];
		}
	    }
	}
    }
        

  For(i,n_ok_edges)
    {
      For(j,n_ok_edges)
	{
	  buff[i*n_ok_edges+j] = -hessian[ok_edges[i]*dim+ok_edges[j]];
	}
    }

  if(!Matinv(buff,n_ok_edges,n_ok_edges))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");      
    }

  For(i,n_ok_edges)
    {
      For(j,n_ok_edges)
	{
	  hessian[ok_edges[i]*dim+ok_edges[j]] = buff[i*n_ok_edges+j];
	}
    }

  /* Approximate variance for very short branches */
  For(i,dim)
    if(!is_ok[i])
      {
	hessian[i*dim+i] = 1./POW(tree->data->init_len,2);
      }

  if(!Matinv(hessian,dim,dim))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");      
    }

  For(i,dim*dim) hessian[i] = -1.0*hessian[i];

/*   For(i,dim) */
/*     { */
/*       PhyML_Printf("[%f] ",tree->t_edges[i]->l); */
/*       For(j,i+1) */
/* 	{ */
/* 	  PhyML_Printf("%12lf ",hessian[i*dim+j]); */
/* 	} */
/*       PhyML_Printf("\n"); */
/*     } */

/*   Matinv(hessian,dim,dim); */

/*   PhyML_Printf("\n"); */

/*   For(i,dim) */
/*     { */
/*       PhyML_Printf("[%f] ",tree->t_edges[i]->l); */
/*       For(j,i+1) */
/* 	{ */
/* 	  PhyML_Printf("%12lf ",hessian[i*dim+j]); */
/* 	} */
/*       PhyML_Printf("\n"); */
/*     } */
/*   Exit("\n"); */


  free(ori_bl);
  free(plus_plus);
  free(minus_minus);
  free(plus_zero);
  free(minus_zero);
  free(plus_minus);
  free(zero_zero);
  free(inc);
  free(buff);
  free(ok_edges);
  free(is_ok);

  return hessian;

}

/*********************************************************/

void Recurr_Hessian_Log(t_node *a, t_node *d, int plus_minus, phydbl *inc, phydbl *res, int *is_ok, t_tree *tree)
{
  int i;
  phydbl ori_l;

  For(i,3)
    if(a->v[i] == d)
      {
	Update_P_Lk(tree,a->b[i],a);

	ori_l = a->b[i]->l;
	if(is_ok[a->b[i]->num])
	  {
	    if(plus_minus > 0) a->b[i]->l += inc[a->b[i]->num];
	    else               a->b[i]->l -= inc[a->b[i]->num];
	    res[a->b[i]->num]  = Lk_At_Given_Edge(a->b[i],tree);
	    res[a->b[i]->num] += Log_Det(is_ok,tree);
	    a->b[i]->l = ori_l;
	    Update_PMat_At_Given_Edge(a->b[i],tree);
	  }
	break;
      }

  if(d->tax) return;
  else 
    For(i,3) 
      if(d->v[i] != a) 
	Recurr_Hessian_Log(d,d->v[i],plus_minus,inc,res,is_ok,tree);
}

/*********************************************************/

phydbl Log_Det(int *is_ok, t_tree *tree)
{
  int i;
  phydbl ldet;

  ldet = 0.0;
/*   For(i,2*tree->n_otu-3) if(is_ok[i]) ldet += LOG(2.*SQRT(tree->t_edges[i]->l)); */
  For(i,2*tree->n_otu-3) if(is_ok[i]) ldet += LOG(tree->t_edges[i]->l);
/*   For(i,2*tree->n_otu-3) if(is_ok[i]) ldet -= LOG(2*tree->t_edges[i]->l); */
  
  return ldet;

}

/*********************************************************/

phydbl Normal_Trunc_Mean(phydbl mu, phydbl sd, phydbl min, phydbl max)
{
  phydbl mean;

  mean = mu + sd * 
    (Dnorm((min-mu)/sd,0.,1.)-Dnorm((max-mu)/sd,0.,1.))/
    (Pnorm((max-mu)/sd,0.,1.)-Pnorm((min-mu)/sd,0.,1.));
  return mean;
}

/*********************************************************/

phydbl Constraint_Normal_Trunc_Mean(phydbl wanted_mu, phydbl sd, phydbl min, phydbl max)
{
  int j;
  phydbl dx,f,fmid,xmid,rtb;
  phydbl x1, x2;

  x1 = min;
  x2 = max;

  f    = Normal_Trunc_Mean(x1,sd,min,max) - wanted_mu;
  fmid = Normal_Trunc_Mean(x2,sd,min,max) - wanted_mu;
  
  if(f*fmid >= 0.0)
    {
      PhyML_Printf("\n. Root must be bracketed for bisection!");
      PhyML_Printf("\n. f=%f fmid=%f",f,fmid);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);

  For(j,100) 
    {
      xmid=rtb+(dx *= 0.5);
      fmid=Normal_Trunc_Mean(xmid,sd,min,max)-wanted_mu;
      if(fmid <= 0.0) rtb=xmid;
      if(fmid > -1.E-10 && fmid < 1.E-10) return rtb;
    }

  Exit("Too many bisections in RTBIS");
  return(-1.);
}

/*********************************************************/

int Matinv(phydbl *x, int n, int m)
{
/* x[n*m]  ... m>=n
*/
   int i,j,k;
   int *irow;
   phydbl ee, t,t1,xmax;
   phydbl det;

   ee = 1.0E-10;
   det = 1.0;
   
   irow = (int *)mCalloc(n,sizeof(int));

   For (i,n)
     {
       xmax = 0.;
       for (j=i; j<n; j++)
         if (xmax < FABS(x[j*m+i]))
	   {
	     xmax = FABS(x[j*m+i]);
	     irow[i]=j;
	   }

      det *= xmax;
      if (xmax < ee)
	{
	  free(irow);
	  PhyML_Printf("\n. Determinant becomes zero at %3d!\t\n", i+1);
	  PhyML_Printf("\n. Failed to invert the matrix.\n");
	  return(0);
	}
      if (irow[i] != i)
	{
	  For (j,m)
	    {
	      t = x[i*m+j];
	      x[i*m+j] = x[irow[i]*m+j];
	      x[irow[i]*m+j] = t;
	    }
	}
      t = 1./x[i*m+i];
      For (j,n)
	{
	  if (j == i) continue;
	  t1 = t*x[j*m+i];
	  For(k,m)  x[j*m+k] -= t1*x[i*m+k];
	  x[j*m+i] = -t1;
	}
      For(j,m)   x[i*m+j] *= t;
      x[i*m+i] = t;
   }                            /* i  */
   for (i=n-1; i>=0; i--)
     {
       if (irow[i] == i) continue;
       For(j,n)
	 {
	   t = x[j*m+i];
	   x[j*m+i] = x[j*m + irow[i]];
	   x[j*m + irow[i]] = t;
	 }
     }

   free(irow);
   return (1);

/*   int i, j, k, lower, upper; */
/*   phydbl temp; */
/*   phydbl *a; */
/*   int nsize; */

/*   nsize = n; */
/*   a = x; */
  
/*   /\*Gauss-Jordan reduction -- invert matrix a in place, */
/*          overwriting previous contents of a.  On exit, matrix a */
/*          contains the inverse.*\/ */
/*   lower = 0; */
/*   upper = nsize-1; */
/*   for(i = lower; i <= upper; i++)  */
/*     { */
/*       temp = 1.0 / a[i*n+i]; */
/*       a[i*n+i] = 1.0; */
/*       for (j = lower; j <= upper; j++)  */
/* 	{ */
/* 	  a[i*n+j] *= temp; */
/* 	} */
/*       for (j = lower; j <= upper; j++)  */
/* 	{ */
/* 	  if (j != i)  */
/* 	    { */
/* 	      temp = a[j*n+i]; */
/* 	      a[j*n+i] = 0.0; */
/* 	      for (k = lower; k <= upper; k++)  */
/* 		{ */
/* 		  a[j*n+k] -= temp * a[i*n+k]; */
/* 		}	       */
/* 	    } */
/* 	} */
/*     } */
  return(1);
}

/*********************************************************/

phydbl *Matrix_Mult(phydbl *A, phydbl *B, int nra, int nca, int nrb, int ncb)
{
  int i,j,k;
  phydbl *C;

  C = (phydbl *)mCalloc(nra*ncb,sizeof(phydbl));

  if(nca != nrb)
    {
      PhyML_Printf("\n. Matrices dimensions don't match.");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");      
    }
  
  For(i,nra)
    For(j,ncb)
       For(k,nca)
         C[i*ncb+j] += A[i*nca+k] * B[k*ncb+j];
  
  return C;
}

/*********************************************************/

phydbl *Matrix_Transpose(phydbl *A, int dim)
{
  phydbl *tA,buff;
  int i,j;

  tA = (phydbl *)mCalloc(dim*dim,sizeof(phydbl));

  For(i,dim*dim) tA[i]=A[i];

  For(i,dim) for(j=i+1;j<dim;j++) 
    {
      buff        = tA[i*dim+j];
      tA[i*dim+j] = tA[j*dim+i];
      tA[j*dim+i]  = buff;
    }

  return tA;
}

/*********************************************************/

phydbl Matrix_Det(phydbl *A, int size, int _log)
{
  phydbl *triA;
  int i;
  phydbl det;

  triA = Cholesky_Decomp(A,size);
  det = 0.0;
  For(i,size) det += LOG(triA[i*size+i]);
  free(triA);
 
  if(_log == NO)
    {
      det = EXP(det);
      return det*det;
    }
  else
    {
      return 2.*det;
    }
}

/*********************************************************/

/* http://en.wikipedia.org/wiki/Multivariate_normal_distribution (Conditional distributions) */
void Normal_Conditional(phydbl *mu, phydbl *cov, phydbl *a, int n, short int *is_1, int n1, phydbl *cond_mu, phydbl *cond_cov)
{
  phydbl *mu1,*mu2;
  phydbl *sig11,*sig12,*sig21,*sig22,*sig12_invsig22,*buff;
  phydbl *ctrd_a;
  phydbl *cond_cov_norder,*cond_mu_norder;
  int    n2;
  int i,j,nr,nc;

  n2 = n-n1;

  mu1             = (phydbl *)mCalloc(n1,   sizeof(phydbl));
  mu2             = (phydbl *)mCalloc(n2,   sizeof(phydbl));
  sig11           = (phydbl *)mCalloc(n1*n1,sizeof(phydbl));
  sig12           = (phydbl *)mCalloc(n1*n2,sizeof(phydbl));
  sig21           = (phydbl *)mCalloc(n2*n1,sizeof(phydbl));
  sig22           = (phydbl *)mCalloc(n2*n2,sizeof(phydbl));
  ctrd_a          = (phydbl *)mCalloc(n2,   sizeof(phydbl)); 
  cond_cov_norder = (phydbl *)mCalloc(n1*n1,sizeof(phydbl));
  cond_mu_norder  = (phydbl *)mCalloc(n1*n1,sizeof(phydbl));

  nr=0;
  For(i,n) { if(!is_1[i]) { ctrd_a[nr] = a[i]-mu[i]; nr++; } }

  nr=0;
  For(i,n) { if( is_1[i]) { mu1[nr] = mu[i]; nr++; } }

  nr=0;
  For(i,n) { if(!is_1[i]) { mu2[nr] = mu[i]; nr++; } }

  nr=0; nc=0;
  For(i,n)
    {
      if(is_1[i])
	{
	  nc = nr;
 	  for(j=i;j<n;j++)
/* 	  nc = 0; */
/* 	  For(j,n) */
	    {
	      if(is_1[j])
		{
		  sig11[nr*n1+nc] = cov[i*n+j];
		  sig11[nc*n1+nr] = cov[i*n+j];
		  nc++;
		}
	    }
	  nr++;
	}
    }


  nr=0; nc=0;
  For(i,n)
    {
      if(is_1[i])
	{
/* 	  nc = nr; */
/*  	  for(j=i;j<n;j++) */
	  nc = 0;
	  For(j,n)
	    {
	      if(!is_1[j])
		{
		  sig12[nr*n2+nc] = cov[i*n+j];
/* 		  sig12[nc*n2+nr] = cov[i*n+j]; */
		  nc++;
		}
	    }
	  nr++;
	}
    }

  nr=0; nc=0;
  For(i,n)
    {
      if(!is_1[i])
	{
/* 	  nc = nr; */
/* 	  for(j=i;j<n;j++) */
	  nc = 0;
	  For(j,n)
	    {
	      if(is_1[j])
		{
		  sig21[nr*n1+nc] = cov[i*n+j];
/* 		  sig21[nc*n1+nr] = cov[i*n+j]; */
		  nc++;
		}
	    }
	  nr++;
	}
    }


  nr=0; nc=0;
  For(i,n)
    {
      if(!is_1[i])
	{
	  nc = nr;
	  for(j=i;j<n;j++)
/* 	  nc = 0; */
/* 	  For(j,n) */
	    {
	      if(!is_1[j])
		{
		  sig22[nr*n2+nc] = cov[i*n+j];
 		  sig22[nc*n2+nr] = cov[i*n+j];
		  nc++;
		}
	    }
	  nr++;
	}
    }

  if(!Matinv(sig22,n2,n2))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");      
    }
  sig12_invsig22 = Matrix_Mult(sig12,sig22,n1,n2,n2,n2);

  buff = Matrix_Mult(sig12_invsig22,ctrd_a,n1,n2,n2,1);
  For(i,n1) cond_mu_norder[i] = mu1[i]+buff[i];
  free(buff);

  buff = Matrix_Mult(sig12_invsig22,sig21,n1,n2,n2,n1);
  For(i,n1) For(j,n1) cond_cov_norder[i*n1+j] = sig11[i*n1+j] - buff[i*n1+j];
  free(buff);

  nr = 0;
  For(i,n) if(is_1[i]) { cond_mu[i] = cond_mu_norder[nr]; nr++; }

  nr = nc = 0;
  For(i,n) 
    {
      if(is_1[i]) 
	{ 
	  nc = 0;
	  For(j,n)
	    {
	      if(is_1[j]) 
		{		  
		  cond_cov[i*n+j] = cond_cov_norder[nr*n1+nc]; 
		  nc++;
		}
	    }
	  nr++;
	}
    }

/*   For(i,n1) */
/*     { */
/*       for(j=i;j<n1;j++) */
/* 	if(FABS(cond_cov_norder[i*n1+j] - cond_cov_norder[j*n1+i]) > 1.E-3) */
/* 	  { */
/* 	    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
/* 	    Warn_And_Exit(""); */
/* 	  } */
/*     } */


  For(i,n)
    {
      for(j=i+1;j<n;j++)
	if(FABS(cond_cov[i*n+j] - cond_cov[j*n+i]) > 1.E-3)
	  {
	    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Warn_And_Exit("");
	  }
    }

  free(mu1);
  free(mu2);
  free(sig11);
  free(sig12);
  free(sig21);
  free(sig22);
  free(ctrd_a);
  free(cond_cov_norder);
  free(cond_mu_norder);
  free(sig12_invsig22);
}


/*********************************************************/

/* http://en.wikipedia.org/wiki/Multivariate_normal_distribution (Conditional distributions) */
void Normal_Conditional_Unsorted(phydbl *mu, phydbl *cov, phydbl *a, int n, short int *is_1, int n1, phydbl *cond_mu, phydbl *cond_cov)
{
  phydbl *mu1,*mu2;
  phydbl *sig11,*sig12,*sig21,*sig22,*sig12_invsig22,*buff;
  phydbl *ctrd_a;
  int    n2;
  int i,j,nr,nc;

  n2 = n-n1;

  mu1             = (phydbl *)mCalloc(n1,   sizeof(phydbl));
  mu2             = (phydbl *)mCalloc(n2,   sizeof(phydbl));
  sig11           = (phydbl *)mCalloc(n1*n1,sizeof(phydbl));
  sig12           = (phydbl *)mCalloc(n1*n2,sizeof(phydbl));
  sig21           = (phydbl *)mCalloc(n2*n1,sizeof(phydbl));
  sig22           = (phydbl *)mCalloc(n2*n2,sizeof(phydbl));
  ctrd_a          = (phydbl *)mCalloc(n2,   sizeof(phydbl)); 

  nr=0;
  For(i,n) { if(!is_1[i]) { ctrd_a[nr] = a[i]-mu[i]; nr++; } }

  nr=0;
  For(i,n) { if( is_1[i]) { mu1[nr] = mu[i]; nr++; } }

  nr=0;
  For(i,n) { if(!is_1[i]) { mu2[nr] = mu[i]; nr++; } }

  nr=0; nc=0;
  For(i,n)
    {
      if(is_1[i])
	{
	  nc = nr;
 	  for(j=i;j<n;j++)
/* 	  nc = 0; */
/* 	  For(j,n) */
	    {
	      if(is_1[j])
		{
		  sig11[nr*n1+nc] = cov[i*n+j];
		  sig11[nc*n1+nr] = cov[i*n+j];
		  nc++;
		}
	    }
	  nr++;
	}
    }


  nr=0; nc=0;
  For(i,n)
    {
      if(is_1[i])
	{
/* 	  nc = nr; */
/*  	  for(j=i;j<n;j++) */
	  nc = 0;
	  For(j,n)
	    {
	      if(!is_1[j])
		{
		  sig12[nr*n2+nc] = cov[i*n+j];
/* 		  sig12[nc*n2+nr] = cov[i*n+j]; */
		  nc++;
		}
	    }
	  nr++;
	}
    }

  nr=0; nc=0;
  For(i,n)
    {
      if(!is_1[i])
	{
/* 	  nc = nr; */
/* 	  for(j=i;j<n;j++) */
	  nc = 0;
	  For(j,n)
	    {
	      if(is_1[j])
		{
		  sig21[nr*n1+nc] = cov[i*n+j];
/* 		  sig21[nc*n1+nr] = cov[i*n+j]; */
		  nc++;
		}
	    }
	  nr++;
	}
    }


  nr=0; nc=0;
  For(i,n)
    {
      if(!is_1[i])
	{
	  nc = nr;
	  for(j=i;j<n;j++)
/* 	  nc = 0; */
/* 	  For(j,n) */
	    {
	      if(!is_1[j])
		{
		  sig22[nr*n2+nc] = cov[i*n+j];
 		  sig22[nc*n2+nr] = cov[i*n+j];
		  nc++;
		}
	    }
	  nr++;
	}
    }

  if(!Matinv(sig22,n2,n2))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");      
    }
  sig12_invsig22 = Matrix_Mult(sig12,sig22,n1,n2,n2,n2);

  buff = Matrix_Mult(sig12_invsig22,ctrd_a,n1,n2,n2,1);
  For(i,n1) cond_mu[i] = mu1[i]+buff[i];
  free(buff);

  buff = Matrix_Mult(sig12_invsig22,sig21,n1,n2,n2,n1);
  For(i,n1) For(j,n1) cond_cov[i*n1+j] = sig11[i*n1+j] - buff[i*n1+j];


  free(mu1);
  free(mu2);
  free(sig11);
  free(sig12);
  free(sig21);
  free(sig22);
  free(ctrd_a);
  free(sig12_invsig22);
}


/*********************************************************/

/* http://en.wikipedia.org/wiki/Multivariate_normal_distribution (Conditional distributions) */
void Get_Reg_Coeff(phydbl *mu, phydbl *cov, phydbl *a, int n, short int *is_1, int n1, phydbl *reg_coeff)
{
  phydbl *sig12,*sig22,*sig12_invsig22;
  int    n2;
  int    i,j,nr,nc;

  n2 = n-n1;

  sig12 = (phydbl *)mCalloc(n1*n2,sizeof(phydbl));
  sig22 = (phydbl *)mCalloc(n2*n2,sizeof(phydbl));

  nr=0; nc=0;
  For(i,n)
    {
      if(is_1[i])
	{
	  nc = 0;
	  For(j,n)
	    {
	      if(!is_1[j])
		{
		  sig12[nr*n2+nc] = cov[i*n+j];
		  nc++;
		}
	    }
	  nr++;
	}
    }


  nr=0; nc=0;
  For(i,n)
    {
      if(!is_1[i])
	{
	  nc = nr;
	  for(j=i;j<n;j++)
	    {
	      if(!is_1[j])
		{
		  sig22[nr*n2+nc] = cov[i*n+j];
 		  sig22[nc*n2+nr] = cov[i*n+j];
		  nc++;
		}
	    }
	  nr++;
	}
    }


  if(!Matinv(sig22,n2,n2))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");      
    }
  sig12_invsig22 = Matrix_Mult(sig12,sig22,n1,n2,n2,n2);


  For(i,n) reg_coeff[i] = 0.0;

/*   nr = 0; */
/*   For(i,n) if(!is_1[i]) { reg_coeff[i] = sig12_invsig22[nr]; nr++; } */

  nc = 0;
  nr = 0;
  For(i,n1) 
    {
      nc = 0;
      For(j,n)
	if(!is_1[j]) 
	  { 
	    reg_coeff[i*n+j] = sig12_invsig22[nr*n2+nc]; 
	    nc++; 
	  }
      nr++;
    }


  if(nc != n2 || nr != n1)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }


  free(sig12);
  free(sig22);
  free(sig12_invsig22);
}


/*********************************************************/

phydbl Norm_Trunc_Sd(phydbl mu, phydbl sd, phydbl a, phydbl b)
{
  phydbl pdfa, pdfb;
  phydbl cdfa, cdfb;
  phydbl ctra, ctrb;
  phydbl cond_var;
  phydbl cdfbmcdfa;

  ctra = (a - mu)/sd;
  ctrb = (b - mu)/sd;

  pdfa = Dnorm(ctra,0.0,1.0);
  pdfb = Dnorm(ctrb,0.0,1.0);

  cdfa = Pnorm(ctra,0.0,1.0);
  cdfb = Pnorm(ctrb,0.0,1.0);

  cdfbmcdfa = cdfb - cdfa;

  if(cdfbmcdfa < SMALL) 
    {
      cdfbmcdfa = SMALL;
      PhyML_Printf("\n. mu=%G sd=%G a=%G b=%G",mu,sd,a,b);
      PhyML_Printf("\n. Numerical precision issue detected.");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    }
	    
  cond_var = sd*sd*(1. + (ctra*pdfa - ctrb*pdfb)/cdfbmcdfa - POW((pdfa - pdfb)/cdfbmcdfa,2));

  return SQRT(cond_var);
}

/*********************************************************/

phydbl Norm_Trunc_Mean(phydbl mu, phydbl sd, phydbl a, phydbl b)
{
  phydbl pdfa, pdfb;
  phydbl cdfa, cdfb;
  phydbl ctra, ctrb;
  phydbl cond_mu;
  phydbl cdfbmcdfa;

  ctra = (a - mu)/sd;
  ctrb = (b - mu)/sd;

  pdfa = Dnorm(ctra,0.0,1.0);
  pdfb = Dnorm(ctrb,0.0,1.0);

  cdfa = Pnorm(ctra,0.0,1.0);
  cdfb = Pnorm(ctrb,0.0,1.0);
  
  cdfbmcdfa = cdfb - cdfa;

  if(cdfbmcdfa < SMALL)
    {
      cdfbmcdfa = SMALL;
      PhyML_Printf("\n. mu=%G sd=%G a=%G b=%G",mu,sd,a,b);
      PhyML_Printf("\n. Numerical precision issue detected.");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    }
  
  cond_mu = mu + sd*(pdfa - pdfb)/cdfbmcdfa;

  return cond_mu;
}

/*********************************************************/

int Norm_Trunc_Mean_Sd(phydbl mu, phydbl sd, phydbl a, phydbl b, phydbl *trunc_mu, phydbl *trunc_sd)
{

  phydbl pdfa, pdfb;
  phydbl cdfa, cdfb;
  phydbl ctra, ctrb;
  phydbl cdfbmcdfa;

  ctra = (a - mu)/sd;
  ctrb = (b - mu)/sd;

  pdfa = Dnorm(ctra,0.0,1.0);
  pdfb = Dnorm(ctrb,0.0,1.0);

  cdfa = Pnorm(ctra,0.0,1.0);
  cdfb = Pnorm(ctrb,0.0,1.0);
  
  cdfbmcdfa = cdfb - cdfa;

  if(cdfbmcdfa < SMALL)
    {
      cdfbmcdfa = SMALL;
      PhyML_Printf("\n. mu=%G sd=%G a=%G b=%G",mu,sd,a,b);
      PhyML_Printf("\n. cdfa=%G cdfb=%G",cdfa,cdfb);
      PhyML_Printf("\n. Numerical precision issue detected.");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      return 0;
    }
  
  *trunc_mu = mu + sd*(pdfa - pdfb)/cdfbmcdfa;
  *trunc_sd = sd*sd*(1. + (ctra*pdfa - ctrb*pdfb)/cdfbmcdfa - POW((pdfa - pdfb)/cdfbmcdfa,2));
  *trunc_sd = SQRT(*trunc_sd);
  return 1;
}
/*********************************************************/
