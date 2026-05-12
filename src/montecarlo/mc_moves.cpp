//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
// standard library header files
#include <vector>

// vampire header files
#include "random.hpp"
#include "sim.hpp"
#include <algorithm>

// Internal header file
#include "internal.hpp"
#include <cmath>
#include <algorithm>

#include "atoms.hpp"
#include "exchange.hpp"

namespace montecarlo{

namespace internal{

// Function declarations
void mc_gaussian(const std::vector<double>&, std::vector<double>&);
void mc_spin_flip(const std::vector<double>&, std::vector<double>&);
void mc_uniform(std::vector<double>&);
void mc_angle(const std::vector<double>&, std::vector<double>&, const double angle);
void mc_hinzke_nowak(const std::vector<double>&, std::vector<double>&);
void mc_adaptive(const std::vector<double>&, std::vector<double>&);
void mc_local_quantized(const int atom, const std::vector<double>&, std::vector<double>&);
void local_quantized_field(const int atom, double& bx, double& by, double& bz);
void mc_local_quantized_heat_bath(const int atom, const std::vector<double>& old_spin, std::vector<double>& new_spin);
///--------------------------------------------------------
///
///  Master function to call desired Monte Carlo move
///
///--------------------------------------------------------
void mc_move(const int atom, const std::vector<double>& old_spin, std::vector<double>& new_spin){

   // Reference enum list for readability
   using namespace montecarlo;

   // Select algorithm using case statement
   switch(algorithm){

      case adaptive:
         mc_adaptive(old_spin, new_spin);
         break;
      case spin_flip:
         mc_spin_flip(old_spin, new_spin);
         break;
      case uniform:
         mc_uniform(new_spin);
         break;
      case angle:
         mc_angle(old_spin, new_spin, montecarlo::internal::delta_angle);
         break;
      case hinzke_nowak:
         mc_hinzke_nowak(old_spin, new_spin);
         break;
      case montecarlo::local_quantized:
         mc_local_quantized(atom, old_spin, new_spin);
         break;
      case montecarlo::local_quantized_heat_bath:
         mc_local_quantized_heat_bath(atom, old_spin, new_spin);
         break;
      default:
         mc_adaptive(old_spin, new_spin);
         break;
   }
   return;
}

/// Angle move
/// Move spin within cone near old position
void mc_angle(const std::vector<double>& old_spin, std::vector<double>& new_spin, const double angle){

   new_spin[0] = old_spin[0] + mtrandom::gaussian() * angle;
   new_spin[1] = old_spin[1] + mtrandom::gaussian() * angle;
   new_spin[2] = old_spin[2] + mtrandom::gaussian() * angle;

   // Calculate new spin length
   const double r = 1.0/sqrt (new_spin[0]*new_spin[0]+new_spin[1]*new_spin[1]+new_spin[2]*new_spin[2]);

   // Apply normalisation
   new_spin[0] *= r;
   new_spin[1] *= r;
   new_spin[2] *= r;

   return;

}

/// Spin flip move
/// Reverse spin direction
void mc_spin_flip(const std::vector<double>& old_spin, std::vector<double>& new_spin){

   new_spin[0]=-old_spin[0];
   new_spin[1]=-old_spin[1];
   new_spin[2]=-old_spin[2];

   return;

}

/// Random move
/// Place spin randomly on unit sphere
void mc_uniform(std::vector<double>& new_spin){

   new_spin[0]=mtrandom::gaussian();
   new_spin[1]=mtrandom::gaussian();
   new_spin[2]=mtrandom::gaussian();

   // Calculate new spin length
   const double r = 1.0/sqrt (new_spin[0]*new_spin[0]+new_spin[1]*new_spin[1]+new_spin[2]*new_spin[2]);

   // Apply normalisation
   for (size_t i=0; i < new_spin.size(); i++) {
      new_spin[i]*=r;
   }

   return;

}

/// Combination move selecting random move from spin_flip, angle and random
///
/// D. Hinzke, U. Nowak, Computer Physics Communications 121–122 (1999) 334–337
/// "Monte Carlo simulation of magnetization switching in a Heisenberg model for small ferromagnetic particles"
///
void mc_hinzke_nowak(const std::vector<double>& old_spin, std::vector<double>& new_spin){

   // Select random move type
   const int pick_move=int(3.0*mtrandom::grnd());

      switch(pick_move){
         case 0:
            mc_spin_flip(old_spin, new_spin);
            break;
         case 1:
            mc_uniform(new_spin);
            break;
         case 2:
            mc_angle(old_spin, new_spin, montecarlo::internal::delta_angle);
            break;
         default:
            mc_angle(old_spin, new_spin, montecarlo::internal::delta_angle);
            break;
      }
      return;
}

//-----------------------------------------------------------------------------------------
/// Generates a new spin from a cone around the old spin, the cone width is
/// derived from the acceptance rate of the previous monte carlo step.
///
/// Adaptive algorithm implemented
/// JD. Alzate-Cardona
/// RFL. Evans
/// D. Sabogal-Suarez
// Implementation by Oscar David Arbeláez E., JD. Alzate-Cardona and R F L Evans 2018
//-----------------------------------------------------------------------------------------
void mc_adaptive(const std::vector<double>& old_spin, std::vector<double>& new_spin){
   // Here we have adaptive_sigma, at least from the monte carlo algorithm
   mc_angle(old_spin, new_spin, montecarlo::internal::adaptive_sigma);
   return;
}

/*void mc_local_quantized(const int atom, const std::vector<double>& old_spin, std::vector<double>& new_spin){

   (void)old_spin; // not needed for this proposal

   const double s = montecarlo::internal::spin_quantum_number;
   const double eps = 1.0e-12;

   // Safety: fall back to uniform move if s is not valid
   if(s <= 0.0){
      mc_uniform(new_spin);
      return;
   }

   // ------------------------------------------------------------
   // Temporary prototype:
   // use sum of all other spins as a crude proxy for local field.
   //
   // This is NOT the final physics you want, but it lets you test
   // the quantized construction before wiring in exchange neighbors.
   // ------------------------------------------------------------
   double bx = 0.0;
   double by = 0.0;
   double bz = 0.0;

   for(int j = 0; j < atoms::num_atoms; ++j){
      if(j == atom) continue;
      bx += atoms::x_spin_array[j];
      by += atoms::y_spin_array[j];
      bz += atoms::z_spin_array[j];
   }

   double bnorm = std::sqrt(bx*bx + by*by + bz*bz);

   // If local field vanishes, choose random orientation
   if(bnorm < eps){
      mc_uniform(new_spin);
      return;
   }

   // Local z-axis along effective field
   double ezx = bx / bnorm;
   double ezy = by / bnorm;
   double ezz = bz / bnorm;

   // Build transverse basis ex, ey
   double ax = 0.0, ay = 0.0, az = 1.0;
   if(std::fabs(ezz) > 0.9){
      ax = 1.0; ay = 0.0; az = 0.0;
   }

   // ex = normalize(a x ez)
   double exx = ay*ezz - az*ezy;
   double exy = az*ezx - ax*ezz;
   double exz = ax*ezy - ay*ezx;
   double enorm = std::sqrt(exx*exx + exy*exy + exz*exz);

   if(enorm < eps){
      mc_uniform(new_spin);
      return;
   }

   exx /= enorm;
   exy /= enorm;
   exz /= enorm;

   // ey = ez x ex
   double eyx = ezy*exz - ezz*exy;
   double eyy = ezz*exx - ezx*exz;
   double eyz = ezx*exy - ezy*exx;

   // Choose m_s from {-s, -s+1, ..., s}
   const int two_s = int(std::round(2.0*s));
   const int nlevels = two_s + 1;
   const int pick = int(nlevels * mtrandom::grnd());

   const double ms = -s + double(pick);

   // Projection and transverse amplitude
   double spar = ms / s;
   if(spar > 1.0) spar = 1.0;
   if(spar < -1.0) spar = -1.0;

   const double sperp = std::sqrt(std::max(0.0, 1.0 - spar*spar));
   const double phi = 2.0 * M_PI * mtrandom::grnd();

   const double c = std::cos(phi);
   const double sn = std::sin(phi);

   new_spin[0] = spar*ezx + sperp*(c*exx + sn*eyx);
   new_spin[1] = spar*ezy + sperp*(c*exy + sn*eyy);
   new_spin[2] = spar*ezz + sperp*(c*exz + sn*eyz);

   // Final normalization for safety
   const double r = 1.0 / std::sqrt(new_spin[0]*new_spin[0] +
                                    new_spin[1]*new_spin[1] +
                                    new_spin[2]*new_spin[2]);

   new_spin[0] *= r;
   new_spin[1] *= r;
   new_spin[2] *= r;

   return;
} */

void local_quantized_field(const int atom, double& bx, double& by, double& bz){

   bx = 0.0;
   by = 0.0;
   bz = 0.0;

   const int start = atoms::neighbour_list_start_index[atom];
   const int end   = atoms::neighbour_list_end_index[atom] + 1;

   const unsigned int ex_type = exchange::get_exchange_type();

   switch(ex_type){

      case exchange::isotropic:
         for(int nn = start; nn < end; ++nn){

            const int natom = atoms::neighbour_list_array[nn];
            const int iid   = atoms::neighbour_interaction_type_array[nn];

            const double Jij = atoms::i_exchange_list[iid].Jij;

            bx += Jij * atoms::x_spin_array[natom];
            by += Jij * atoms::y_spin_array[natom];
            bz += Jij * atoms::z_spin_array[natom];
         }
         break;

      case exchange::vectorial:
         for(int nn = start; nn < end; ++nn){

            const int natom = atoms::neighbour_list_array[nn];
            const int iid   = atoms::neighbour_interaction_type_array[nn];

            bx += atoms::v_exchange_list[iid].Jij[0] * atoms::x_spin_array[natom];
            by += atoms::v_exchange_list[iid].Jij[1] * atoms::y_spin_array[natom];
            bz += atoms::v_exchange_list[iid].Jij[2] * atoms::z_spin_array[natom];
         }
         break;

      case exchange::tensorial:
         for(int nn = start; nn < end; ++nn){

            const int natom = atoms::neighbour_list_array[nn];
            const int iid   = atoms::neighbour_interaction_type_array[nn];

            const double sx = atoms::x_spin_array[natom];
            const double sy = atoms::y_spin_array[natom];
            const double sz = atoms::z_spin_array[natom];

            bx += atoms::t_exchange_list[iid].Jij[0][0] * sx
                + atoms::t_exchange_list[iid].Jij[0][1] * sy
                + atoms::t_exchange_list[iid].Jij[0][2] * sz;

            by += atoms::t_exchange_list[iid].Jij[1][0] * sx
                + atoms::t_exchange_list[iid].Jij[1][1] * sy
                + atoms::t_exchange_list[iid].Jij[1][2] * sz;

            bz += atoms::t_exchange_list[iid].Jij[2][0] * sx
                + atoms::t_exchange_list[iid].Jij[2][1] * sy
                + atoms::t_exchange_list[iid].Jij[2][2] * sz;
         }
         break;

      default:
         // fallback: no field
         bx = 0.0;
         by = 0.0;
         bz = 0.0;
         break;
   }

   return;
}

void mc_local_quantized(const int atom, const std::vector<double>& old_spin, std::vector<double>& new_spin){

   (void)old_spin;

   const double s = montecarlo::internal::spin_quantum_number;
   const double eps = 1.0e-14;

   // guard against invalid input
   if(s <= 0.0){
      mc_uniform(new_spin);
      return;
   }

   // compute real local exchange field
   double bx = 0.0;
   double by = 0.0;
   double bz = 0.0;
   local_quantized_field(atom, bx, by, bz);

   const double bnorm = std::sqrt(bx*bx + by*by + bz*bz);

   // if field vanishes, choose random spin
   if(bnorm < eps){
      mc_uniform(new_spin);
      return;
   }

   // local quantization axis ez = B / |B|
   const double ezx = bx / bnorm;
   const double ezy = by / bnorm;
   const double ezz = bz / bnorm;

   // choose an auxiliary vector not parallel to ez
   double ax = 0.0;
   double ay = 0.0;
   double az = 1.0;

   if(std::fabs(ezz) > 0.9){
      ax = 1.0;
      ay = 0.0;
      az = 0.0;
   }

   // ex = normalize(a x ez)
   double exx = ay*ezz - az*ezy;
   double exy = az*ezx - ax*ezz;
   double exz = ax*ezy - ay*ezx;

   const double exnorm = std::sqrt(exx*exx + exy*exy + exz*exz);

   if(exnorm < eps){
      mc_uniform(new_spin);
      return;
   }

   exx /= exnorm;
   exy /= exnorm;
   exz /= exnorm;

   // ey = ez x ex
   const double eyx = ezy*exz - ezz*exy;
   const double eyy = ezz*exx - ezx*exz;
   const double eyz = ezx*exy - ezy*exx;

   // allowed projections m_s = -s, -s+1, ..., s
   const int two_s = int(std::round(2.0 * s));
   const int nlevels = two_s + 1;

   int pick = int(nlevels * mtrandom::grnd());
   if(pick >= nlevels) pick = nlevels - 1;

   const double ms = -s + double(pick);

   double spar = ms / s;
   if(spar > 1.0) spar = 1.0;
   if(spar < -1.0) spar = -1.0;

   const double sperp = std::sqrt(std::max(0.0, 1.0 - spar*spar));
   const double phi = 2.0 * M_PI * mtrandom::grnd();

   const double c = std::cos(phi);
   const double sn = std::sin(phi);

   new_spin[0] = spar*ezx + sperp*(c*exx + sn*eyx);
   new_spin[1] = spar*ezy + sperp*(c*exy + sn*eyy);
   new_spin[2] = spar*ezz + sperp*(c*exz + sn*eyz);

   // normalize for safety
   const double r = 1.0 / std::sqrt(new_spin[0]*new_spin[0]
                                  + new_spin[1]*new_spin[1]
                                  + new_spin[2]*new_spin[2]);

   new_spin[0] *= r;
   new_spin[1] *= r;
   new_spin[2] *= r;

   return;
}

void mc_local_quantized_heat_bath(
   const int atom,
   const std::vector<double>& old_spin,
   std::vector<double>& new_spin
){

   (void)old_spin;

   const double s = montecarlo::internal::spin_quantum_number;
   const double eps = 1.0e-14;

   if(s <= 0.0){
      mc_uniform(new_spin);
      return;
   }

   // ------------------------------------------------------------
   // Compute local exchange field
   // ------------------------------------------------------------
   double bx = 0.0;
   double by = 0.0;
   double bz = 0.0;

   local_quantized_field(atom, bx, by, bz);

   const double bnorm = std::sqrt(bx*bx + by*by + bz*bz);

   if(bnorm < eps){
      mc_uniform(new_spin);
      return;
   }

   const double ezx = bx / bnorm;
   const double ezy = by / bnorm;
   const double ezz = bz / bnorm;

   // ------------------------------------------------------------
   // Build local transverse basis
   // ------------------------------------------------------------
   double ax = 0.0;
   double ay = 0.0;
   double az = 1.0;

   if(std::fabs(ezz) > 0.9){
      ax = 1.0;
      ay = 0.0;
      az = 0.0;
   }

   double exx = ay*ezz - az*ezy;
   double exy = az*ezx - ax*ezz;
   double exz = ax*ezy - ay*ezx;

   const double exnorm = std::sqrt(exx*exx + exy*exy + exz*exz);

   if(exnorm < eps){
      mc_uniform(new_spin);
      return;
   }

   exx /= exnorm;
   exy /= exnorm;
   exz /= exnorm;

   const double eyx = ezy*exz - ezz*exy;
   const double eyy = ezz*exx - ezx*exz;
   const double eyz = ezx*exy - ezy*exx;

   // ------------------------------------------------------------
   // Effective inverse temperature in VAMPIRE MC convention
   //
   // Existing MC accepts with:
   // exp(-DE * rescaled_material_kBTBohr)
   //
   // Here bnorm has same energy convention as calculate_spin_energy,
   // so use local beta_eff = mu_s/muB * kBTBohr.
   // ------------------------------------------------------------
   const int imaterial = atoms::type_array[atom];

   const double alpha = montecarlo::internal::temperature_rescaling_alpha[imaterial];
   const double Tc    = montecarlo::internal::temperature_rescaling_Tc[imaterial];

   const double rescaled_temperature =
      sim::temperature < Tc
      ? Tc * std::pow(sim::temperature / Tc, alpha)
      : sim::temperature;

   const double kBTBohr =
      9.27400915e-24 / (rescaled_temperature * 1.3806503e-23);

   const double beta_eff =
      montecarlo::internal::mu_s_SI[imaterial] * 1.07828231e23 * kBTBohr;

   // ------------------------------------------------------------
   // Sample m from P(m) ∝ exp(beta_eff * bnorm * m)
   // ------------------------------------------------------------
   const int two_s = int(std::round(2.0 * s));
   const int nlevels = two_s + 1;

   std::vector<double> weights(nlevels, 0.0);

   // stabilize exponent
   double max_arg = -1.0e300;

   for(int p = 0; p < nlevels; ++p){
      const double m = -s + double(p);
      const double arg = beta_eff * bnorm * m;
      if(arg > max_arg) max_arg = arg;
   }

   double Z = 0.0;

   for(int p = 0; p < nlevels; ++p){
      const double m = -s + double(p);
      const double arg = beta_eff * bnorm * m;
      weights[p] = std::exp(arg - max_arg);
      Z += weights[p];
   }

   double rnum = mtrandom::grnd() * Z;

   int pick = nlevels - 1;
   double cumulative = 0.0;

   for(int p = 0; p < nlevels; ++p){
      cumulative += weights[p];
      if(rnum <= cumulative){
         pick = p;
         break;
      }
   }

   const double ms = -s + double(pick);

   double spar = ms / s;
   if(spar > 1.0) spar = 1.0;
   if(spar < -1.0) spar = -1.0;

   const double sperp = std::sqrt(std::max(0.0, 1.0 - spar*spar));

   const double phi = 2.0 * M_PI * mtrandom::grnd();

   const double c = std::cos(phi);
   const double sn = std::sin(phi);

   new_spin[0] = spar*ezx + sperp*(c*exx + sn*eyx);
   new_spin[1] = spar*ezy + sperp*(c*exy + sn*eyy);
   new_spin[2] = spar*ezz + sperp*(c*exz + sn*eyz);

   const double norm =
      std::sqrt(new_spin[0]*new_spin[0]
              + new_spin[1]*new_spin[1]
              + new_spin[2]*new_spin[2]);

   if(norm > 0.0){
      const double inv = 1.0 / norm;
      new_spin[0] *= inv;
      new_spin[1] *= inv;
      new_spin[2] *= inv;
   }
   else{
      mc_uniform(new_spin);
   }

   return;
}

} //end of namespace internal

} //end of namespace montecarlo
