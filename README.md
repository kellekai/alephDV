alphaDV
======
**alphaDV** evaluates the ALEPH data from the hadronic Tau-Decay with FOPT and CIPT incorporating duality violations (DV).

This code evaluates the ALEPH data from the hadronic Tau-Decay with FOPT and CIPT 
incorporating duality violations (DV).

The code is written in Fortran and uses functions from LAPACK, QUADPACK and MINUIT

The Data can be found at url -> http://aleph.web.lal.in2p3.fr/tau/specfun13.html
  
The evaluation is described in detail in the document Masterthesis2016.pdf in the doc folder

## INSTALLATION
Download or clone repository. Move to the base directory and type:

  * make all

  To remove the libraries and executives type:

  * make remove

## USAGE
To execute the program just type:

 * ./'MODELNAME'.run or ./'MODELNAME'_CIPT.run

## CONFIGURATION
To reset the configuration files, please type:

 * make config
