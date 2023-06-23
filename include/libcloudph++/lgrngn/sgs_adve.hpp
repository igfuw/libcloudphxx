#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
//<listing>
    enum class sgs_adve_t { undefined, GA17, ST_periodic }; 
//</listing>
    // GA17 - SGS velocity modeled with a Langevin equation as in Grabowski & Abade 2017. TKE dissipation rate taken from the diss_rate array; characteristic length scale taken from the SGS_mix_len profile.
    // ST_periodic - synthetic turbulence: SGS velocity is a sum of Fourier modes, periodic within a box. constant and uniform diss_rate; constant and uniform minimum and maximum wavelengths
  };
};
