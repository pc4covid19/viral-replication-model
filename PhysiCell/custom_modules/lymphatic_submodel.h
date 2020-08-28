#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#include "./submodel_data_structures.h" 

// #include "./immune_submodels.h"  
// #include "./receptor_dynamics.h"  
// #include "./internal_viral_dynamics.h"  
// #include "./internal_viral_response.h"   

#ifndef __lymphatic_submodel__
#define __lymphatic_submodel__

extern Submodel_Information lymphatic_submodel_info; 


void lymphatic_submodel_setup( void ); 

// don't put into individual cell models -- needs to be a fast process 
void lymphatic_model( Cell* pCell, Phenotype& phenotype, double dt );

// this needs to be done on faster time scale; 
void lymphatic_main_model( double dt ); 

#endif 