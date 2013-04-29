//  EnergyCorrelator Package
//  Questions/Comments?  larkoski@mit.edu gavin.salam@cern.ch jthaler@jthaler.net
//
//  Copyright (c) 2013
//  Andrew Larkoski, Gavin Salam, and Jesse Thaler
//
//  $Id$
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
#include <sstream>
using namespace std;

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
         for (unsigned int j = i + 1; j < particles.size(); j++) { //note offset by one so that angle is never called on identical pairs
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
   if (_strategy == storage_array) {
   
      // For N > 2, fill static storage array to save computation time.

      // we'll need the (signed) integer number of particles
      int n = particles.size();

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
            //angleStore[j][i] = NAN; // no need to store other size.
            angleStore[j][i] = angleStore[i][j]; // no need to store other size.
         }
      }

      // now do recursion
      if (_N == 3) {
         // for (unsigned int i = 0; i < particles.size(); i++) {
         //    for (unsigned int j = i + 1; j < particles.size(); j++) {
         //       double ans_ij = energyStore[i]
         //                       * energyStore[j]
         //                       * angleStore[i][j];
         //       for (unsigned int k = j + 1; k < particles.size(); k++) {
         //          answer += ans_ij
         //                    * energyStore[k]
         //                    * angleStore[i][k]
         //                    * angleStore[j][k];
         //       }
         //    }
         // } 
         // new method
         for (int i = n - 1; i >= 0; i--) {
            for (int j = i - 1; j >= 0 ; j--) {
               double ans_ij = energyStore[i]
                               * energyStore[j]
                               * angleStore[i][j];
               for (int k = j - 1; k >= 0; k--) {
                  answer += ans_ij
                            * energyStore[k]
                            * angleStore[i][k]
                            * angleStore[j][k];
               }
            }
         }       
      } else if (_N == 4) {
         // for (unsigned int i = 0; i < particles.size(); i++) {
         //    for (unsigned int j = i + 1; j < particles.size(); j++) {
         //       double ans_ij = energyStore[i]
         //                       * energyStore[j]
         //                       * angleStore[i][j];
         //       for (unsigned int k = j + 1; k < particles.size(); k++) {
         //          double ans_ijk = ans_ij
         //                         * energyStore[k]
         //                         * angleStore[i][k]
         //                         * angleStore[j][k];
         //          for (unsigned int l = k + 1; l < particles.size(); l++) {
         //             answer += ans_ijk
         //                               * energyStore[l]
         //                               * angleStore[i][l]
         //                               * angleStore[j][l]
         //                               * angleStore[k][l];
         //          }
         //       }
         //    }
         // } 
         // new method
         for (int i = n - 1; i >= 0; i--) {
            for (int j = i - 1; j >= 0 ; j--) {
               double ans_ij = energyStore[i]
                               * energyStore[j]
                               * angleStore[i][j];
               for (int k = j - 1; k >= 0; k--) {
                  double ans_ijk = ans_ij
                                 * energyStore[k]
                                 * angleStore[i][k]
                                 * angleStore[j][k];
                  for (int l = k - 1; l >= 0; l--) {
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

   } else if (_strategy == slow) {
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

double EnergyCorrelator::energy(const PseudoJet& jet) const {
   if (_measure == pt_R) {
      return jet.perp();
   }  else if (_measure == E_theta) {
      return jet.e();
   } else {
      assert(false);
      return NAN;
   }
}

double EnergyCorrelator::angleSquared(const PseudoJet& jet1, const PseudoJet& jet2) const {
   if (_measure == pt_R) {
      return jet1.squared_distance(jet2);
   } else if (_measure == E_theta) {
      // doesn't seem to be a fastjet built in for this
      double dot = jet1.px()*jet2.px() + jet1.py()*jet2.py() + jet1.pz()*jet2.pz();
      double norm1 = sqrt(jet1.px()*jet1.px() + jet1.py()*jet1.py() + jet1.pz()*jet1.pz());
      double norm2 = sqrt(jet2.px()*jet2.px() + jet2.py()*jet2.py() + jet2.pz()*jet2.pz());
      
      double costheta = dot/(norm1 * norm2);
      if (costheta > 1.0) costheta = 1.0; // Need to handle case of numerical overflow
      double theta = acos(costheta);
      return theta*theta;    
        
   } else {
      assert(false);
      return NAN;
   }
}

string EnergyCorrelator::description_parameters() const {
  ostringstream oss;
  oss << "N=" << _N << ", beta=" << _beta;

  if      (_measure == pt_R)    oss << ", pt_R measure";
  else if (_measure == E_theta) oss << ", E_theta measure";
  else throw Error("unrecognized measure");

  if      (_strategy == slow)   oss << " and 'slow' strategy";
  else if (_strategy == storage_array)   oss << " and 'storage_array' strategy";
  else throw Error("unrecognized strategy");

  return oss.str();
}

string EnergyCorrelator::description() const {
  ostringstream oss;
  oss << "Energy Correlator ECF(N,beta) for ";
  oss << description_parameters();
  return oss.str();
}

string EnergyCorrelatorRatio::description() const {
  ostringstream oss;
  oss << "Energy Correlator ratio ECF(N+1,beta)/ECF(N,beta) for ";
  oss << EnergyCorrelator(_N,_beta,_measure,_strategy).description_parameters();
  return oss.str();
}

string EnergyCorrelatorDoubleRatio::description() const {
  ostringstream oss;
  oss << "Energy Correlator double ratio ECF(N-1,beta)ECF(N+1,beta)/ECF(N,beta)^2 for ";
  oss << EnergyCorrelator(_N,_beta,_measure,_strategy).description_parameters();
  return oss.str();
}




} // namespace contrib

FASTJET_END_NAMESPACE
