#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#include "./submodel_data_structures.h" 

// #include "./immune_submodels.h"  
// #include "./receptor_dynamics.h"  
// #include "./internal_viral_dynamics.h"  
// #include "./internal_viral_response.h"   

#ifndef __viral_replication_submodel__
#define __viral_replication_submodel__

extern Submodel_Information viral_replication_submodel_info; 


void viral_replication_submodel_setup( void ); 

// don't put into individual cell models -- needs to be a fast process 
void viral_replication_model( Cell* pCell, Phenotype& phenotype, double dt );

// this needs to be done on faster time scale; 
void viral_replication_main_model( double dt ); 

#endif 