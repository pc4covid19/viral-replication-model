# viral-replication-model

Integrate virus replication model (BioNetGen/SBML/ODEs) with PhysiCell.

## Overview

* if you plan to use libRoadrunner as your SBML solver, you should run the script `setup_libroadrunner.py` in the `/beta` directory and follow instructions there. You may then need to edit `Make-sbml` (in the root, `/PhysiCell` directory) to have LIBRR* paths be valid.

* `PhysiCell/config_with_sbml.xml` is a PhysiCell configuration file that contains the
basic "linkage" to the SBML model (search for "sbml" or "SBML" in it).

* to test without the SBML model, compile with `make` and run `COVID19` to see the dummy `viral_replication_main_model` being invoked. Files of interest are `custom_modules/viral_replication_submodel.{h,cpp}` and `main.cpp`.

* to test with the SBML model, compile with `make -f Make-sbml`, have your relevant env var pointing to the Roadrunner libs you installed, and run `COVID19_sbml config_with_sbml.xml`. Note that `Make-sbml` defines the compiler macro `-D LIBROADRUNNER` which allows for conditional compilation of `#ifdef LIBROADRUNNER` code (rf. `intracellular/PhysiCell_intracellular.h`). However, we still need to make meaningful mappings between the SBML species and PhysiCell cells' custom data and, of course, do something meaningful in the custom code.

## Changes from original 3.2 code (and latest immune-response model)
* added `custom_modules/viral_replication_submodel.{h,cpp}`

## additional notes
* added `beta/setup_libroadrunner.py`
* added `intracellular/PhysiCell_intracellular.h`
* added `std::string sbml_filename;` into `class Cell_Definition` (in `core/PhysiCell_cell.h`)
* parse `<molecular>` XML in `core/PhysiCell_cell.cpp` 
* incorporate intracellular (SBML, libRoadrunner) info into `custom.{h,cpp}`
  * see `assign_SBML_model( Cell* pC )` in `custom.cpp`
  * see extra code in `create_cell_types( void )` where we obtain the desired SBML species' indices.
* added `Make-sbml`