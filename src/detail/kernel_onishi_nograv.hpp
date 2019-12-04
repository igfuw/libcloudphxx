#pragma once

#if defined(__NVCC__)
#  include <math_constants.h>
#endif

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    { 
      //  adapted from fortran code by Onishi 
      //   returns value of turbulent kernel without gravitational settling
      //   2PI*R^2*<|Wr|>g(R)
      //
      //   g11(R) is formulated in Onishi (2005) and Onishi et al. (2006).
      //   Later the parameters were modified based on the new high-Re data obtained
      //   by Onishi et al. (2013) and Onishi & Vassilicos (2014).
      // Reference
      //   Onishi (2005) phD thesis, Kyoto Univ.
      //   Onishi et al. (2006), Trans. JSME, 72, 2441-2448 (in Japanese)
      //   Onishi et al. (2013), J. Comput. Phys., 242, 809-827
      //   Onishi & Vassilicos (2014), J. Fluid Mech., 745, 279-299
      //
      //   For <|Wr|>: Wang et al. (2000), J. Fluid Mech., 415, 117-153
      //   For g12(R): Zhou et al. (2001), J. Fluid Mech., 433, 77-104
      template<class real_t>
      BOOST_GPU_ENABLED
      real_t kernel_onishi_nograv(const real_t &r1, const real_t &r2, const real_t &Re_l, const real_t &eps, real_t dnu, real_t ratio_den)
      {
        if(eps < 1e-10) return 0.;

#if !defined(__NVCC__)
        using std::max;
        using std::pow;
        using std::sqrt;
        using std::cbrt;
        using std::exp;
#endif

        real_t PI = 
#if !defined(__NVCC__)
        pi<real_t>();
#else
        CUDART_PI;
#endif

        real_t urms;
        real_t leta;
        real_t Wr, WrA2,WrS2;
        real_t taup1,taup2,theta1,theta2,gamma;
        real_t cw;
        real_t CR,tauk,Te,phi;
        real_t A1,A2,A3; //,ReREF,CREF,Cdash;
        real_t CA,CB,alpha;
        real_t StA,StB,St1,St2;
        real_t y11,y21,y31,y12,y22,y32;
        real_t za1,za2,zb1,zb2,xai;
        real_t gR, gR1,gR2,RG12;
    
        urms  = sqrt(Re_l/sqrt(15./dnu/eps)); //RMS of u [m/s];
        CR    = r1+r2;               //collision radius [m];
        taup1 = ratio_den * 4.*r1*r1 /18. /dnu; //particle relaxation time [s];
        taup2 = ratio_den * 4.*r2*r2 /18. /dnu;
        leta  = pow(dnu*dnu*dnu/eps,real_t(1/4));  //Kolmogorov scale [m];
        tauk  = leta*leta/dnu;         //Kolmogorov time  [s];
        Te    = Re_l*tauk/sqrt(15.);
  
  // --- <|Wr|> in Wang et al. (2000)JFM;
        theta1 =2.5*taup1/Te;
        theta2 =2.5*taup2/Te;
        phi = max(theta2/theta1,theta1/theta2);
        cw = 1.+0.6*exp(-pow(phi-1.,1.5));
        gamma = 0.183*urms*urms/(dnu*dnu/leta/leta);
        gamma = phi * gamma;
  
        WrS2 = (dnu*dnu*CR*CR) /(leta*leta*leta*leta) /15.;
        WrA2 = urms*urms * gamma/(gamma-1.)                           
             *( (theta1+theta2)-4.*theta1*theta2/(theta1+theta2)     
             *sqrt( (1.+theta1+theta2)/(1.+theta1)/(1.+theta2) ) ) 
             * ( 1./(1.+theta1)/(1.+theta2)                        
                -1./(1.+gamma*theta1)/(1.+gamma*theta2) );
        WrA2 = cw*WrA2;
        WrA2 = WrA2/3.; //Ayala 2008;
        Wr  =sqrt(2./PI*(WrA2 +WrS2));
  
  // --- g11(R) & g22(R), i.e., monodisperse;
  // --- empirical parameters for g11(R);
  //      A1 = 95.0 //Onishi (2005)phD thesis;
        A1 = 110.0; //2014;
  //      A2 = 0.32.0 //Onishi (2005)phD thesis;
        A2 = 0.38; //2014;
        A3 = 0.16;
  //in Onishi (2005)phD thesis;
  //      ReREF = 54.3;
  //      CREF = 16.2;
  //      Cdash = CREF*ReREF**(1.0/3.0) / sqrt(A1*A2);
  //      alpha = (5.0/6.0 * log10(Re_l) - log10(Cdash)) &;
  //          / log10(2.0)                           //;
        alpha = log10(0.26*sqrt(Re_l))/log10(2.0); //2014;
        alpha = max(alpha, real_t(1.e-20));
  
        CA = 0.06*pow(Re_l, real_t(0.30)); //2014;
        CB = 0.4;
  
        StA = pow(A2/A1 * Re_l, real_t(0.25));
        const real_t hlpr = cbrt(A2/A3);
        StB = hlpr*hlpr * cbrt(Re_l);
  
        St1 = taup1/ tauk; //Stokes number;
        St2 = taup2/ tauk;
  
        if(St2 <= StA)
        {
          y11 = A1 * St1*St1;
          y21 = 0.;
        }
        else
        {
          y11 = 0.;
          y21 = A2 * Re_l / (St1*St1);
        }
        y31 = A3 * sqrt(Re_l/St1);
  
        if(St1 <= StA)
        {
          y12 = A1 * St2*St2;
          y22 = 0.;
        }
        else
        {
          y12 = 0.;
          y22 = A2 * Re_l / (St2*St2);
        }
        y32 = A3 * sqrt(Re_l/St2);
  
        za1 = 0.5 * (1. - tanh((log10(St1) - log10(StA)) / CA) );
        zb1 = 0.5 * (1. + tanh((log10(St1) - log10(StB)) / CB) );
        za2 = 0.5 * (1. - tanh((log10(St2) - log10(StA)) / CA) );
        zb2 = 0.5 * (1. + tanh((log10(St2) - log10(StB)) / CB) );
  
        gR1 = y11*pow(za1,alpha) + y21*pow(real_t(1.)-za1,alpha) + y31*zb1 + 1.;
        gR2 = y12*pow(za2,alpha) + y22*pow(real_t(1.)-za2,alpha) + y32*zb2 + 1.;
  
  //g12=f(g11,g22) in Zhou et al. (2001)JFM      ;
        xai = max(taup2/taup1,taup1/taup2);     //xai:alpha;
        RG12= 2.6*exp(-xai)+0.205*exp(-0.0206*xai) 
            * 0.5*(1.0+tanh(xai-3.0));
        gR = 1. + RG12*sqrt(gR1-1.)*sqrt(gR2-1.);
  
        return  2.*PI*CR*CR *Wr *gR;     //turbulent Kc (no gravity);
      }
    };
  };
};
      

