/* ***************************************************************************
 *
 *  Copyright (C) 2013-2016 University of Dundee
 *  All rights reserved. 
 *
 *  This file is part of SAMoS (Soft Active Matter on Surfaces) program.
 *
 *  SAMoS is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  SAMoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ****************************************************************************/


/*!
 * \file pair_actreact_nematic_potential.hpp
 * \author Silke Henkes & Chinmay Pabshettiwar, silkehenkes@gmail.com
 * \date 23-October-2025
 * \brief Computation of the PairActReactNematicPotential class
 */ 

#include "pair_actreact_nematic_potential.hpp"
#include "constraint.hpp"


//! \param dt time step sent by the integrator 
void PairActReactNematicPotential::compute(double dt)
{
  int N = m_system->size();
  double poli, polj; // activities
  double rinti, rintj; // interaction ranges multipliers
  double ai, aj;
  double alpha_i = 1.0;  // phase in factor for particle i
  double alpha_j = 1.0;  // phase in factor for particle j
  double alpha = 1.0; // phase in factor for pair interaction (see below)
  
  for  (int i = 0; i < N; i++)
    {
	  Particle& p = m_system->get_particle(i);
      // Not in fact computing energy for this non-conservative "potential"
	  if (m_system->compute_per_particle_energy())
	  {
		  p.set_pot_energy("actreact_nematic",0.0);
	   }
   }
  
  m_potential_energy = 0.0;
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    // type parameters
    // std::cout << m_has_part_params << "type i " << pi.get_type() << endl;
    if (m_has_part_params)  
    {
      poli = m_particle_params[pi.get_type()-1].p;
      rinti = m_particle_params[pi.get_type()-1].r_int;
    }
    else 
    {
        poli = m_p;
        rinti = m_r_int;
    }
    
    // A note on phasing in: m_val, the value object, has been pre-set with the number of phase in time steps
    // It will give a linear interpolation between 0 and 1 based on the current age of the particle in time steps
    // Second note: I could set the interpolation between 0.5 and 1, however then pair_vertex potential would be dubious
    // Instead this is done manually here
    // The pair ABP inherits the whole phase-in mechanism from the other potentials, so that everything is as compatible as possible
    if (m_phase_in)
      alpha_i = 0.5*(1.0 + m_val->get_val(static_cast<int>(pi.age/dt)));
    ai = pi.get_radius();
    vector<int>& neigh = m_nlist->get_neighbours(i);
    for (unsigned int j = 0; j < neigh.size(); j++)
    {
      Particle& pj = m_system->get_particle(neigh[j]);
      if (m_phase_in)
      {
        alpha_j = 0.5*(1.0 + m_val->get_val(static_cast<int>(pj.age/dt)));
        // Determine global phase in factor: particles start at 0.5 strength (both daugthers of a division replace the mother)
        // Except for the interaction between daugthers which starts at 0
        if (alpha_i < 1.0 && alpha_j < 1.0)
	        alpha = alpha_i + alpha_j - 1.0;
	      else 
	        alpha = alpha_i*alpha_j;
      }
      // data for particle 2
      aj = pj.get_radius();
      // type parameters
      if (m_has_part_params) 
      {
        polj = m_particle_params[pj.get_type()-1].p;
        rintj = m_particle_params[pj.get_type()-1].r_int;
      }
      else 
      {
        polj = m_p;
        rintj = m_r_int;
      }
      

      double dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;
      m_system->apply_periodic(dx,dy,dz);
      double r_sq = dx*dx + dy*dy + dz*dz;
      double r = sqrt(r_sq);

      double ai_p_aj;
      if (!m_use_particle_radii)
        ai_p_aj = 2.0;
      else
        ai_p_aj = ai+aj;

      double rintij = 0.5*(rinti+rintj)*ai_p_aj;
      double b;
      double fax, fay, faz;
      double nidotr;
      double nidotrp;
      double njdotr;
      double njdotrp;

      // local geometry
      double rijx, rijy, rijz; // bond unit vector
      double rijpx, rijpy, rijpz; // perpendicular bond vector
      double Nx, Ny, Nz; // unit normal at middle of bond vector
      double N2;
      double pref;


      // "Parallel" forces from (Q_i + Q_j) . r_ij:
      // F_ij^par = kpar/(\sqrt{3}R) beta(rij) (sigma_i + sigma_j) . hat{rij}
      // where sigma_i = p_i n_i ni_i is the nematic stress tensor
      // "Perpendicular" forces coming from the integral representation of \nabla . \sigma = F:
      // F_{ij}^perp = kperp/(\sqrt{3}R) beta(rij) (sigma_i + sigma_j) . (hat{rij} x N),

      if (r < rintij)
      {
        b = (rintij-r)/rintij;  // force prefactor (positive and between 0 and 1)

        // unit vector along bond
        rijx = dx;
        rijy = dy;
        rijz = dz;
        if (!m_3d) {
          // Unit normal in this location - average it from the two local particles
          Nx = 0.5*(pi.Nx + pj.Nx);
          Ny = 0.5*(pi.Ny + pj.Ny);
          Nz = 0.5*(pi.Nz + pj.Nz);
          // normalise
          N2 = sqrt(Nx*Nx+Ny*Ny+Nz*Nz);
          Nx = Nx/N2;
          Ny = Ny/N2;
          Nz = Nz/N2;
          // Components of rijperp
          rijpx = rijy*Nz-rijz*Ny;
          rijpy = rijz*Nx-rijx*Nz;
          rijpz = rijx*Ny-rijy*Nx;
        }
        pref = b;

        // parallel forces
        nidotr = pi.nx*rijx+pi.ny*rijy+pi.nz*rijz;
        njdotr = pj.nx*rijx+pj.ny*rijy+pj.nz*rijz;

        fax = m_k_par*pref*(poli*nidotr*pi.nx + polj*njdotr*pj.nx);
        fay = m_k_par*pref*(poli*nidotr*pi.ny + polj*njdotr*pj.ny);
        faz = m_k_par*pref*(poli*nidotr*pi.nz + polj*njdotr*pj.nz);

        // perpendicular forces
        if (!m_3d) {
          nidotrp = pi.nx*rijpx+pi.ny*rijpy+pi.nz*rijpz;
          njdotrp = pj.nx*rijpx+pj.ny*rijpy+pj.nz*rijpz;

          fax += m_k_perp*pref*(poli*nidotrp*pi.nx + polj*njdotrp*pj.nx);
          fay += m_k_perp*pref*(poli*nidotrp*pi.ny + polj*njdotrp*pj.ny);
          faz += m_k_perp*pref*(poli*nidotrp*pi.nz + polj*njdotrp*pj.nz);
        }

        // Handle force
        pi.fx += alpha*fax;
        pi.fy += alpha*fay;
        pi.fz += alpha*faz;
        // Use 3d Newton's law: this is an action-reaction and reciprocal active force
        pj.fx -= alpha*fax;
        pj.fy -= alpha*fay;
        pj.fz -= alpha*faz;
        // handle torques
        if (m_torques) {
          pi.tau_x += 0.5*alpha*(dy*faz-dz*fay);
          pi.tau_y += 0.5*alpha*(-dx*faz+dz*fax);
          pi.tau_z += 0.5*alpha*(dx*fay-dy*fax);
          // and the other direction - carefully: this is the *same* torque, not the opposite sign ...
          pj.tau_x += 0.5*alpha*(dy*faz-dz*fay);
          pj.tau_y += 0.5*alpha*(-dx*faz+dz*fax);
          pj.tau_z += 0.5*alpha*(dx*fay-dy*fax);
        }
      }
      
    }
  }
}
