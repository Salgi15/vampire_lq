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

#ifndef MONTECARLO_INTERNAL_H_
#define MONTECARLO_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the montecarlo module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <vector>

// Vampire headers
#include "montecarlo.hpp"

// montecarlo module headers
#include "internal.hpp"

namespace montecarlo{

   //-------------------------------------------------------------------------
   // Internal shared variables
   //-------------------------------------------------------------------------


   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------

      //Materials Variables
      extern int num_materials;
      extern std::vector<double> temperature_rescaling_alpha;
      extern std::vector<double> temperature_rescaling_Tc;
      extern std::vector<double> mu_s_SI;

      //MC Variables
      extern double delta_angle;    // Tuned angle for Monte Carlo trial move
      extern double adaptive_sigma; // sigma trial width for adaptive move
      extern double spin_quantum_number;

      extern std::vector<double> Sold;
      extern std::vector<double> Snew;

      //MC-MPI variables
      extern std::vector<std::vector<int> > c_octants; //Core atoms of each octant
      extern std::vector<std::vector<int> > b_octants; //Boundary atoms of each octant

      extern bool enable_collective_moves;

      extern double collective_move_probability;

      extern double collective_move_amplitude;

      extern int collective_move_qmax;
      extern long long collective_move_attempts;
      extern long long collective_move_accepts;

      void mc_collective_rotation_parallel(

         std::vector<double>& x_spin_array,

         std::vector<double>& y_spin_array,

         std::vector<double>& z_spin_array,

         std::vector<int>& type_array

      );
      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------
      void mc_move(const int atom, const std::vector<double>&, std::vector<double>&);

   } // end of internal namespace

} // end of montecarlo namespace

#endif //MONTECARLO_INTERNAL_H_
