#include <iomanip>

#include "./lymphatic_submodel.h" 


using namespace PhysiCell; 

std::string lymphatic_model_version = "0.1.0"; 

Submodel_Information lymphatic_submodel_info; 

void lymphatic_submodel_setup( void )
{
		// set version 
	lymphatic_submodel_info.name = "lymphatic model"; 
	lymphatic_submodel_info.version = lymphatic_model_version; 

		// set functions 
	lymphatic_submodel_info.main_function = lymphatic_main_model;
	lymphatic_submodel_info.phenotype_function = NULL; // pushed into the "main" model  
	lymphatic_submodel_info.mechanics_function = NULL;

		// what microenvironment variables do you need 
	lymphatic_submodel_info.microenvironment_variables.push_back( "virion" );

		// what cell variables and parameters do you need? 
	lymphatic_submodel_info.cell_variables.push_back( "unbound_external_ACE2" ); 
	
	lymphatic_submodel_info.register_model(); 
	
	return; 
}

void lymphatic_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int lung_epithelial_type = get_cell_definition( "lung epithelium" ).type; 
	
	// bookkeeping -- find microenvironment variables we need

	static int nV_external = microenvironment.find_density_index( "virion" ); 
	
	static int nV_internal = pCell->custom_data.find_variable_index( "virion" ); 

	// bookkeeping -- find custom data we need 
	
	static int nR_EU = pCell->custom_data.find_variable_index( "unbound_external_ACE2" ); 
	
	// do nothing if dead 
	if( phenotype.death.dead == true )
	{ return; } 

	// if not lung epithelium, do nothing 
	// if( pCell->type != lung_epithelial_type )
	// { return; } 
	
	// actual model goes here 
	

    //------------------------------------------------------------------
    // dummy code to show how to use libRoadrunner to solve a SBML model

	// static int idx_glucose = 1;
	// static int idx_oxygen = 3;
	rrc::RRVectorPtr vptr;
	rrc::RRCDataPtr result;  // start time, end time, and number of points

	// std::cout << "------ solve SBML model ------" << std::endl;

	// pC->phenotype.molecular.model_rr = rrHandle;  // assign the intracellular model to each cell
	// vptr = rrc::getFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr);
	// std::cout << "--- before updating:" << std::endl;
	// for (int idx=0; idx<vptr->Count; idx++)
	// 	std::cout << idx << ", " << vptr->Data[idx] << std::endl;

	// vptr->Data[idx_oxygen] += 0.1;
	// rrc::setFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr, vptr);

	// vptr = rrc::getFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr);
	// // std::cout << vptr->Count << std::endl;
	// std::cout << "--- after updating oxygen:" << std::endl;
	// for (int idx=0; idx<vptr->Count; idx++)
	// 	std::cout << idx << ", " << vptr->Data[idx] << std::endl;

	// int oxygen_i = microenvironment.find_density_index( "oxygen" ); 
	// int glucose_i = microenvironment.find_density_index( "glucose" ); 

	int vi = microenvironment.nearest_voxel_index(pCell->position);
	// double oxy_val = microenvironment(vi)[oxygen_i];
	// double glucose_val = microenvironment(vi)[glucose_i];
	// std::cout << "oxy_val at voxel of cell = " << oxy_val << std::endl;
	// std::cout << "glucose_val at voxel of cell = " << glucose_val << std::endl;

	// vptr->Data[idx_oxygen] = oxy_val;
	// vptr->Data[idx_glucose] = glucose_val;
	// rrc::setFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr, vptr);

	// result = rrc::simulateEx (pCell->phenotype.molecular.model_rr, 0, 10, 10);  // start time, end time, and number of points
	// int index = 0;
	// // Print out column headers... typically time and species.
	// for (int col = 0; col < result->CSize; col++)
	// {
	// 	std::cout << std::left << std::setw(15) << result->ColumnHeaders[index++];
	// }
	// std::cout << "\n";

	// index = 0;
	// // Print out the data
	// for (int row = 0; row < result->RSize; row++)
	// {
	// 	for (int col = 0; col < result->CSize; col++)
	// 	{
	// 		std::cout << std::left << std::setw(15) << result->Data[index++];
	// 		// if (col < result->CSize -1)
	// 		// {
	// 		// 	// std::cout << "\t";
	// 		// 	std::cout << "  ";
	// 		// }
	// 	}
	// 	std::cout << "\n";
	// }
	// int idx = (result->RSize - 1) * result->CSize + 1;
	// std::cout << "Saving last energy value (cell custom var) = " << result->Data[idx] << std::endl;
	// pCell->custom_data[energy_vi]  = result->Data[idx];

	
	return; 
}

void lymphatic_main_model( double dt )
{
	std::cout << "---------- lymphatic_main_model: " << PhysiCell_globals.current_time << std::endl;
	#pragma omp parallel for 
	for( int n=0; n < (*all_cells).size() ; n++ )
	{
		Cell* pC = (*all_cells)[n]; 
		if( pC->phenotype.death.dead == false )
		{ lymphatic_model( pC, pC->phenotype , dt ); }
	}
	
	return; 
}
	
