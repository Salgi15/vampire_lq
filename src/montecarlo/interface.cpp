//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Adam Laverack and Richard Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>
#include <cstdlib>

// Vampire headers
#include "montecarlo.hpp"
#include "errors.hpp"
#include "vio.hpp"

// montecarlo module headers
#include "internal.hpp"

namespace montecarlo{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for montecarlo module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="montecarlo";
      if(key!=prefix) return false;

      //--------------------------------------------------------------------
      std::string test="algorithm";
      if( word == test ){
         test = "adaptive";
         if(value == test){
            algorithm = adaptive;
            return true;
         }
         test = "spin-flip";
         if( value == test ){
            algorithm = spin_flip;
            return true;
         }
         test = "uniform";
         if( value == test ){
            algorithm = uniform;
            return true;
         }
         test = "angle";
         if( value == test ){
            algorithm = angle;
            return true;
         }
         test = "hinzke-nowak";
         if( value == test ){
            algorithm = hinzke_nowak;
            return true;
         }
         test = "local-quantized";
         if( value == test ){
            algorithm = local_quantized;
            return true;
         }
         test = "local-quantized-heat-bath";
         if( value == test ){
            algorithm = local_quantized_heat_bath;
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error - value for \'montecarlo:" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"adaptive\"" << std::endl;
            std::cerr << "\t\"spin-flip\"" << std::endl;
            std::cerr << "\t\"uniform\"" << std::endl;
            std::cerr << "\t\"angle\"" << std::endl;
            std::cerr << "\t\"hinzke-nowak\"" << std::endl;
            std::cerr << "\t\"local-quantized\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }
      test = "constrain-by-grain";
      if( word == test ){
         // enable cmc with grain level rather than global constraints
         cmc::constrain_by_grain = true;
         return true;
      }
      test = "spin-quantum-number";
      if( word == test ){
         internal::spin_quantum_number = atof(value.c_str());
         return true;
      }

      if(word == "enable-collective-moves"){
         montecarlo::internal::enable_collective_moves = true;
         return true;
      }

      if(word == "collective-move-probability"){
         montecarlo::internal::collective_move_probability = atof(value.c_str());
         return true;
      }

      if(word == "collective-move-amplitude"){
         montecarlo::internal::collective_move_amplitude = atof(value.c_str());
         return true;
      }

      if(word == "collective-move-qmax"){
         montecarlo::internal::collective_move_qmax = atoi(value.c_str());
         return true;
      }

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index){

      // add prefix string
      std::string prefix="material:";

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

   long long get_collective_attempts(){

      return internal::collective_move_attempts;

   }

   long long get_collective_accepts(){

      return internal::collective_move_accepts;

   }

} // end of montecarlo namespace
