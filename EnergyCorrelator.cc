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


// initialize storage arrays
double EnergyCorrelator::energyStore[NPARTICLE_STORE];
double EnergyCorrelator::angleStore[NPARTICLE_STORE][NPARTICLE_STORE];


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

   // For N >= 2, fill static storage array to save computation time.

   if (particles.size() >= NPARTICLE_STORE) {
      std::cerr << "ERROR:  EnergyCorrelator is only hard coded to handle " << NPARTICLE_STORE << " particles in an event"  << std::endl;
      assert (particles.size() < NPARTICLE_STORE);
   }
   
   if (_N > 4) {
      std::cerr << "ERROR:  EnergyCorrelator is only hard coded for N = 0,1,2,3,4"  << std::endl;
      assert(_N <= 4);
   }
   
   // create storage array
   for (unsigned int i = 0; i < particles.size(); i++) {
      energyStore[i] = energy(particles[i]);
      for (unsigned int j = i; j < particles.size(); j++) {
         angleStore[i][j] = pow(angle(particles[i],particles[j]), _beta);
         angleStore[j][i] = angleStore[i][j];
      }
   }

   // now do the recursion
   for (unsigned int i = 0; i < particles.size(); i++) {
      for (unsigned int j = i; j < particles.size(); j++) {
         double ans_ij = energyStore[i]
                         * energyStore[j]
                         * angleStore[i][j];
         if (_N == 2) answer += ans_ij;
         else {
            for (unsigned int k = j; k < particles.size(); k++) {
               double ans_ijk = ans_ij
                              * energyStore[k]
                              * angleStore[i][k]
                              * angleStore[j][k];
               if (_N == 3) answer += ans_ijk;
               else {
                  for (unsigned int l = k; l < particles.size(); l++) {
                     double ans_ijkl = ans_ijk
                                       * energyStore[l]
                                       * angleStore[i][l]
                                       * angleStore[j][l]
                                       * angleStore[k][l];
                     assert(_N == 4);
                     answer += ans_ijkl;
                  }
               }
            }
         }
      }
   }   
   return answer;
}

} // namespace contrib

FASTJET_END_NAMESPACE
