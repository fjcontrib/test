//  EnergyCorrelator Package
//  Questions/Comments?  larkoski@mit.edu gavin.salam@cern.ch jthaler@jthaler.net
//
//  Copyright (c) 2013
//  Andrew Larkoski, Gavin Salam, and Jesse Thaler
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#include "EnergyCorrelator.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{


double EnergyCorrelator::result(const PseudoJet& jet) const {

   // get N = 0 case out of the way
   if (_N == 0) return 1.0;

   // find constituents 
   std::vector<fastjet::PseudoJet> particles = jet.constituents();
   double answer = 0.0;

   // take care of N = 1 case.
   if (_N == 1) {
      for (unsigned int i = 0; i < particles.size(); i++) {
         answer += energy(particles[i]);
      }
      return answer;
   }

   // take care of N = 2 case.
   if (_N == 2) {
      for (unsigned int i = 0; i < particles.size(); i++) {
         for (unsigned int j = i + 1; j < particles.size(); j++) {
            answer += energy(particles[i])
                      * energy(particles[j])
                      * pow(angleSquared(particles[i],particles[j]), _beta/2.0);
         }
      }   
      return answer;
   }
   

   // if N > 4, then throw error
   if (_N > 4) {
      std::cerr << "ERROR:  EnergyCorrelator is only hard coded for N = 0,1,2,3,4"  << std::endl;
      assert(_N <= 4);
   }

   
   // Now deal with N = 3 and N = 4.  Different options if storage array is used or not.  
   if (_method == ec_storage_array) {
   
         // For N > 2, fill static storage array to save computation time.

      // Make energy storage
      std::vector<double> energyStore;
      energyStore.resize(particles.size());

      // Make angular storage
      std::vector< std::vector<double> > angleStore;
      angleStore.resize(particles.size());
      for (unsigned int i = 0; i < angleStore.size(); i++) {
         angleStore[i].resize(particles.size());
      }
      
      // Fill storage with energy/angle information
      for (unsigned int i = 0; i < particles.size(); i++) {
         energyStore[i] = energy(particles[i]);
         for (unsigned int j = i+1; j < particles.size(); j++) {
            angleStore[i][j] = pow(angleSquared(particles[i],particles[j]), _beta/2.0);
            angleStore[j][i] = NAN; // no need to store other size.
         }
      }

      // now do recursion
      if (_N == 3) {
         for (unsigned int i = 0; i < particles.size(); i++) {
            for (unsigned int j = i + 1; j < particles.size(); j++) {
               double ans_ij = energyStore[i]
                               * energyStore[j]
                               * angleStore[i][j];
               for (unsigned int k = j + 1; k < particles.size(); k++) {
                  answer += ans_ij
                            * energyStore[k]
                            * angleStore[i][k]
                            * angleStore[j][k];
               }
            }
         }       
      } else if (_N == 4) {
         for (unsigned int i = 0; i < particles.size(); i++) {
            for (unsigned int j = i + 1; j < particles.size(); j++) {
               double ans_ij = energyStore[i]
                               * energyStore[j]
                               * angleStore[i][j];
               for (unsigned int k = j + 1; k < particles.size(); k++) {
                  double ans_ijk = ans_ij
                                 * energyStore[k]
                                 * angleStore[i][k]
                                 * angleStore[j][k];
                  for (unsigned int l = k + 1; l < particles.size(); l++) {
                     answer += ans_ijk
                                       * energyStore[l]
                                       * angleStore[i][l]
                                       * angleStore[j][l]
                                       * angleStore[k][l];
                  }
               }
            }
         } 

      
      } else {
         assert(_N <= 4);
      }

   } else if (_method == ec_simple) {
      if (_N == 3) {
         for (unsigned int i = 0; i < particles.size(); i++) {
            for (unsigned int j = i + 1; j < particles.size(); j++) {
               double ans_ij = energy(particles[i])
                               * energy(particles[j])
                               * pow(angleSquared(particles[i],particles[j]), _beta/2.0);
               for (unsigned int k = j + 1; k < particles.size(); k++) {
                  answer += ans_ij
                            * energy(particles[k])
                            * pow(angleSquared(particles[i],particles[k]), _beta/2.0)
                            * pow(angleSquared(particles[j],particles[k]), _beta/2.0);
               }
            }
         }
      } else if (_N == 4) {
         for (unsigned int i = 0; i < particles.size(); i++) {
            for (unsigned int j = i + 1; j < particles.size(); j++) {
               double ans_ij = energy(particles[i])
                               * energy(particles[j])
                               * pow(angleSquared(particles[i],particles[j]), _beta/2.0);
               for (unsigned int k = j + 1; k < particles.size(); k++) {
                  double ans_ijk = ans_ij
                                   * energy(particles[k])
                                   * pow(angleSquared(particles[i],particles[k]), _beta/2.0)
                                   * pow(angleSquared(particles[j],particles[k]), _beta/2.0);
                  for (unsigned int l = k + 1; l < particles.size(); l++) {
                     answer += ans_ijk
                               * energy(particles[l])
                               * pow(angleSquared(particles[i],particles[l]), _beta/2.0)
                               * pow(angleSquared(particles[j],particles[l]), _beta/2.0)
                               * pow(angleSquared(particles[k],particles[l]), _beta/2.0);
                  }
               }
            }
         }
      } else {
         assert(_N <= 4);
      } 
   } else {
      assert(false);
   }
   
   return answer;
}

} // namespace contrib

FASTJET_END_NAMESPACE
