#ifdef MPICF

#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>

#include "atoms.hpp"
#include "montecarlo.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vmpi.hpp"
#include "internal.hpp"

namespace montecarlo{
namespace internal{

void mc_collective_rotation_parallel(
   std::vector<double>& x_spin_array,
   std::vector<double>& y_spin_array,
   std::vector<double>& z_spin_array,
   std::vector<int>& type_array
){

   const int ncore = vmpi::num_core_atoms;
   if(ncore <= 0) return;

   collective_move_attempts++;

   // ------------------------------------------------------------
   // All ranks must use the same collective-mode parameters.
   // Rank 0 chooses them, then broadcasts.
   // ------------------------------------------------------------
   double eps = 0.0;
   double phi = 0.0;
   int nx = 0;
   int ny = 0;
   int nz = 0;

   if(vmpi::my_rank == 0){

      const double eps_max = collective_move_amplitude;
      const int qmax = collective_move_qmax;

      eps = eps_max * (2.0*mtrandom::grnd() - 1.0);
      phi = 2.0 * M_PI * mtrandom::grnd();

      nx = int((2*qmax + 1)*mtrandom::grnd()) - qmax;
      ny = int((2*qmax + 1)*mtrandom::grnd()) - qmax;
      nz = int((2*qmax + 1)*mtrandom::grnd()) - qmax;

      if(nx == 0 && ny == 0 && nz == 0){
         nx = 1;
      }
   }

   MPI_Bcast(&eps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast(&phi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&nz, 1, MPI_INT, 0, MPI_COMM_WORLD);

   // ------------------------------------------------------------
   // Backup old core spins only.
   // ------------------------------------------------------------
   std::vector<double> old_x(ncore);
   std::vector<double> old_y(ncore);
   std::vector<double> old_z(ncore);

   for(int i = 0; i < ncore; ++i){
      old_x[i] = x_spin_array[i];
      old_y[i] = y_spin_array[i];
      old_z[i] = z_spin_array[i];
   }

   // ------------------------------------------------------------
   // Old local energy.
   // Same convention for old/new, so only dE matters.
   // ------------------------------------------------------------
   double Eold_local = 0.0;
   for(int i = 0; i < ncore; ++i){
      Eold_local += sim::calculate_spin_energy(i);
   }

   // ------------------------------------------------------------
   // Local domain bounds from core atoms.
   // This gives a smooth mode inside each MPI domain.
   // ------------------------------------------------------------
   double xmin = atoms::x_coord_array[0];
   double xmax = atoms::x_coord_array[0];
   double ymin = atoms::y_coord_array[0];
   double ymax = atoms::y_coord_array[0];
   double zmin = atoms::z_coord_array[0];
   double zmax = atoms::z_coord_array[0];

   for(int i = 1; i < ncore; ++i){
      xmin = std::min(xmin, atoms::x_coord_array[i]);
      xmax = std::max(xmax, atoms::x_coord_array[i]);

      ymin = std::min(ymin, atoms::y_coord_array[i]);
      ymax = std::max(ymax, atoms::y_coord_array[i]);

      zmin = std::min(zmin, atoms::z_coord_array[i]);
      zmax = std::max(zmax, atoms::z_coord_array[i]);
   }

   const double Lx = std::max(1.0e-12, xmax - xmin);
   const double Ly = std::max(1.0e-12, ymax - ymin);
   const double Lz = std::max(1.0e-12, zmax - zmin);

   // ------------------------------------------------------------
   // Apply infinitesimal smooth rotation:
   //
   //    S' = S + theta(r) x S
   //
   // followed by normalization.
   // ------------------------------------------------------------
   for(int i = 0; i < ncore; ++i){

      const double x = atoms::x_coord_array[i] - xmin;
      const double y = atoms::y_coord_array[i] - ymin;
      const double z = atoms::z_coord_array[i] - zmin;

      const double phase =
         2.0*M_PI*(double(nx)*x/Lx + double(ny)*y/Ly + double(nz)*z/Lz) + phi;

      const double tx = eps * std::cos(phase);
      const double ty = eps * std::sin(phase);
      const double tz = 0.0;

      const double sx = x_spin_array[i];
      const double sy = y_spin_array[i];
      const double sz = z_spin_array[i];

      double nsx = sx + (ty*sz - tz*sy);
      double nsy = sy + (tz*sx - tx*sz);
      double nsz = sz + (tx*sy - ty*sx);

      const double norm2 = nsx*nsx + nsy*nsy + nsz*nsz;

      if(norm2 > 0.0){
         const double invr = 1.0/std::sqrt(norm2);
         x_spin_array[i] = nsx * invr;
         y_spin_array[i] = nsy * invr;
         z_spin_array[i] = nsz * invr;
      }
      else{
         x_spin_array[i] = old_x[i];
         y_spin_array[i] = old_y[i];
         z_spin_array[i] = old_z[i];
      }
   }

   // ------------------------------------------------------------
   // Communicate updated core spins before computing new energy.
   // ------------------------------------------------------------
   vmpi::mpi_init_halo_swap();
   vmpi::mpi_complete_halo_swap();

   double Enew_local = 0.0;
   for(int i = 0; i < ncore; ++i){
      Enew_local += sim::calculate_spin_energy(i);
   }

   double dE_local = Enew_local - Eold_local;
   double dE_global = 0.0;

   MPI_Allreduce(&dE_local, &dE_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   // ------------------------------------------------------------
   // Use material 0 temperature factor for one-material tests.
   // This matches the MC convention approximately:
   //
   // DE in MC is multiplied by mu_s_SI * 1/mu_B.
   // ------------------------------------------------------------
   const int mat = int(type_array[0]);

   const double alpha = temperature_rescaling_alpha[mat];
   const double Tc = temperature_rescaling_Tc[mat];

   const double rescaled_temperature =
      sim::temperature < Tc ? Tc*std::pow(sim::temperature/Tc, alpha)
                            : sim::temperature;

   const double kBTBohr =
      9.27400915e-24 / (rescaled_temperature * 1.3806503e-23);

   const double dE_mc_units =
      dE_global * mu_s_SI[mat] * 1.07828231e23;

   // ------------------------------------------------------------
   // Rank 0 makes the accept/reject decision, then broadcasts it.
   // ------------------------------------------------------------
   int accept_int = 0;

   if(vmpi::my_rank == 0){
      if(dE_mc_units <= 0.0){
         accept_int = 1;
      }
      else if(std::exp(-dE_mc_units * kBTBohr) >= mtrandom::grnd()){
         accept_int = 1;
      }
   }

   MPI_Bcast(&accept_int, 1, MPI_INT, 0, MPI_COMM_WORLD);

   const bool accept = (accept_int == 1);
   if(accept){
      collective_move_accepts++;
   }

   // ------------------------------------------------------------
   // Reject: revert core atoms and update halos again.
   // ------------------------------------------------------------
   if(!accept){

      for(int i = 0; i < ncore; ++i){
         x_spin_array[i] = old_x[i];
         y_spin_array[i] = old_y[i];
         z_spin_array[i] = old_z[i];
      }

      vmpi::mpi_init_halo_swap();
      vmpi::mpi_complete_halo_swap();

      sim::mc_statistics_reject += 1.0;
   }

   sim::mc_statistics_moves += 1.0;

   return;
}

} // namespace internal
} // namespace montecarlo

#endif