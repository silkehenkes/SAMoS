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
 * \author Silke Henkes + Chinmay & Sander, silkehenkes@gmail.com
 * \date 14-October-2025
 * \brief Declaration of PairActReactNematicPotential class
 */ 

#ifndef __PAIR_ACTREACT_NEMATIC_POTENTIAL_HPP__
#define __PAIR_ACTREACT_NEMATIC_POTENTIAL_HPP__

#include <cmath>

#include "pair_potential.hpp"

using std::make_pair;
using std::sqrt;

//! Structure that handles parameters for the abp action reaction pair "potential"
struct ActReactNematicParameters
{
  double p;
  double r_int;
  double k_par;
  double k_perp;
};


/*! Pair ActReact nematic is an action-reaction preserving implementation of fully nematic cell interactions.
* That makes it by symmetry a nematic active stress (!), and a reciprocal interaction (!).
* Individual directions of particles (n) evolve with whatever the integrator says, e.g nematic alignment or even ABP dynamics
* V1: force density, not force ...
* They lead to the following pair forces:
* "Parallel" forces from (Q_i + Q_j) . r_ij:
* F_ij^par = kpar/(\sqrt{3}R) beta(rij) (sigma_i + sigma_j) . hat{rij}
* where sigma_i = p_i n_i ni_i is the nematic stress tensor
* "Perpendicular" forces coming from the integral representation of \nabla . \sigma = F:
* F_{ij}^perp = kperp/(\sqrt{3}R) beta(rij) (sigma_i + sigma_j) . (hat{rij} x N),
* where N is the local normal to the surface
*
* * They lead to the following pair forces:
* "Parallel" forces from (Q_i + Q_j) . r_ij:
* F_ij^par = kpar beta(rij) (sigma_i + sigma_j) . rij
* where sigma_i = p_i n_i ni_i is the nematic stress tensor
* "Perpendicular" forces coming from the integral representation of \nabla . \sigma = F:
* F_{ij}^perp = kperp beta(rij) (sigma_i + sigma_j) . (rij x N),
* where N is the local normal to the surface
*
* CAREFUL: No clean full 3d version of this interaction. For 3d (keyword 3d), this part is turnef off.
* For results that are as intended, set r_int to the interaction radius of the mechanical potential(s) in the system 
* in particular: r_int = 1+2 eps for soft_attractive, map to re_fact as r_int = 1 + 2*(re_fact-1)  [=1.3 for re_fact=1.15]
* and make use of options use_particle_radii and phase_in. 
 */
class PairActReactNematicPotential : public PairPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param val Value control object (for phasing in)
  //! \param param Contains information about all parameters (k)
  PairActReactNematicPotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, ValuePtr val, pairs_type& param) : PairPotential(sys, msg, nlist, val, param)
  {
    m_known_params.push_back("p");
    m_known_params.push_back("r_int");
    m_known_params.push_back("k_par");
    m_known_params.push_back("k_perp");
    m_known_params.push_back("torques");
    m_known_params.push_back("3d");
    m_known_params.push_back("use_particle_radii");
    m_known_params.push_back("phase_in");
    string param_test = this->params_ok(param);
    if (param_test != "")
    {
      m_msg->msg(Messenger::ERROR,"Parameter \""+param_test+"\" is not a valid parameter for pair actreact nematic potential.");
      throw runtime_error("Unknown parameter \""+param_test+"\" in pair actreact nematic potential.");
    }
    if (param.find("p") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No activity specified for pair actreact nematic potential. Setting it to 0.1.");
      m_p = 0.1;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global activity force specified for pair actreact nematic potential "+param["p"]+".");
      m_p = lexical_cast<double>(param["p"]);
    }
    m_msg->write_config("potential.pair.actreact_nematic.p",lexical_cast<string>(m_p));
    if (param.find("r_int") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No potential range (r_int) specified for pair actreact nematic potential. Setting it to 1.5.");
      m_r_int = 1.5;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global potential range (r_int) for pair actreact nematic is set to "+param["r_int"]+".");
      m_r_int = lexical_cast<double>(param["r_int"]);
    }
    m_msg->write_config("potential.pair.actreact_nematic.r_int",lexical_cast<string>(m_r_int));
    if (param.find("k_par") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No k_par parallel coupling coefficient defined for pair actreact nematic potential. Setting it to 1.0.");
      m_k_par = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global k_par parallel coupling coefficient for pair actreact nematic is set to "+param["k_par"]+".");
      m_k_par = lexical_cast<double>(param["k_par"]);
    }
    m_msg->write_config("potential.pair.actreact_nematic.k_par",lexical_cast<string>(m_k_par));
    if (param.find("k_perp") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No k_perp perpendicular coupling coefficient defined for pair actreact nematic potential. Setting it to 1.0.");
      m_k_perp = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global k_perp perpendicular coupling coefficient for pair actreact nematic is set to "+param["k_perp"]+".");
      m_k_perp = lexical_cast<double>(param["k_perp"]);
    }
    m_msg->write_config("potential.pair.actreact_nematic.k_perp",lexical_cast<string>(m_k_perp));
     if (param.find("torques") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No use of torques specified. Setting torques to False.");
      m_torques = false;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Use of pair interaction: torques set to true");
      m_torques = true;
    }
    if (param.find("3d") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"System is not set to 3d. Assuming 2d surface. Setting 3d to False.");
      m_3d = false;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Computing active nematic forces in 3d");
      m_3d = true;
    }

    if (param.find("use_particle_radii") != param.end())
    {
      m_msg->msg(Messenger::WARNING,"action reaction nematic pair potential is set to use particle radii to control its range. Parameter r_int will be ignored.");
      m_use_particle_radii = true;
      m_msg->write_config("potential.pair.actreact_nematic.use_radii","true");
    }
    if (param.find("phase_in") != param.end())
    {
      m_msg->msg(Messenger::INFO,"action reaction nematic pair potential. Gradually phasing in the potential for new particles.");
      m_phase_in = true;
      m_msg->write_config("potential.pair.actreact_nematic.phase_in","true");
    }    
    
    m_particle_params = new ActReactNematicParameters[m_ntypes];
    //std::cout << "Initial constructor parameter values:" << endl;
    for (int i = 0; i < m_ntypes; i++)
    {
      m_particle_params[i].p = m_p;
      m_particle_params[i].r_int = m_r_int;
      //std::cout << "type " << i << " p " << m_p << " r_int " << m_r_int << endl;
    }  
    
  }

  virtual ~PairActReactNematicPotential()
  {
    delete [] m_particle_params;
  }
  
  //! Set type parameters data for individual particles
  void set_type_parameters(pairs_type& pair_param)
  {
    map<string,double> param;
    int type;
    
    if (pair_param.find("type") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type has not been defined for type specific parameters in pair actreact nematic pair potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    
    type = lexical_cast<int>(pair_param["type"]);
    //std::cout << "Setting pair parametes of type:" << type << endl;
        
    if (pair_param.find("p") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"pair actreact nematic potential. Setting activity "+pair_param["p"]+" for particles of type "+lexical_cast<string>(type)+".");
      param["p"] = lexical_cast<double>(pair_param["p"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"pair actreact nematic potential. Using default activity ("+lexical_cast<string>(m_p)+") for particles of type "+lexical_cast<string>(type)+".");
      param["p"] = m_p;
    }
    m_msg->write_config("potential.pair.actreact_nematic.type_"+pair_param["type"]+".push",lexical_cast<string>(param["p"]));
    if (pair_param.find("r_int") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"pair actreact nematic potential. Setting interaction radius "+pair_param["r_int"]+" for particles of type "+lexical_cast<string>(type)+".");
      param["r_int"] = lexical_cast<double>(pair_param["r_int"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"pair actreact nematic potential. Using default interaction radius ("+lexical_cast<string>(m_r_int)+") for particles of type "+lexical_cast<string>(type)+".");
      param["r_int"] = m_r_int;
    }
    m_msg->write_config("potential.pair.actreact_nematic.type_"+pair_param["type"]+".push",lexical_cast<string>(param["r_int"]));

        
    m_particle_params[type-1].p = param["p"];
    m_particle_params[type-1].r_int = param["r_int"];
    //std::cout << "type " << type << " p " << param["p"] << " r_int " << param["r_int"] << endl;
        
    m_has_part_params = true;
  }

  //! Set pair parameters data for pairwise interactions   .. and do nothing. Otherwise the virtual void deities in the base class are unhappy.
  void set_pair_parameters(pairs_type& pair_param)
  {  
    //std::cout << "Going through empty pair parameter setting" << endl;
  }
                                                                                                                
  
  
  //! Returns true since soft potential needs neighbour list
  bool need_nlist() { return true; }
  
  //! Computes potentials and forces for all particles
  void compute(double);
  
  
private:
       
  double m_p;                       //!< activity
  double m_r_int;                       //!< potential range
  double m_k_par;                       //!< parallel coupling coefficient
  double m_k_perp;                       //!< perpendicular coupling coefficient
  bool m_has_part_params;           //!< true if type specific particle parameters are given
  bool m_torques;               // whether or not to include the torques
  bool m_3d;                   // whether or not to compute 3d forces
  ActReactNematicParameters*  m_particle_params;   //!< type specific particle parameters
     
};

typedef shared_ptr<PairActReactNematicPotential> PairActReactNematicPotentialPtr;

#endif
