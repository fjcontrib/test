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

#ifndef __FASTJET_CONTRIB_ENERGYCORRELATOR_HH__
#define __FASTJET_CONTRIB_ENERGYCORRELATOR_HH__

#include <fastjet/internal/base.hh>
#include "fastjet/FunctionOfPseudoJet.hh"

#include <string>
#include <climits>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

//------------------------------------------------------------------------
/// \class EnergyCorrelator
/// Calculates ECF(N,beta).
///
/// EnergyCorrelator(int N, double beta, Measure measure)
/// Called ECF(N,beta) in the publication.
/// N is the multiplicity, beta is the angular exponent, and
/// measure = pt_R (default) or E_theta sets how energies and angles are determined.
class EnergyCorrelator : public FunctionOfPseudoJet<double> {

public:

  enum Measure {
    pt_R,
    E_theta
  };

  enum Strategy {
    slow,
    storage_array
  };

private:

   int _N;
   double _beta;
   Measure _measure;
   Strategy _strategy;

   double energy(const PseudoJet& jet) const {
      if (_measure == pt_R) {
         return jet.perp();
      }  else if (_measure == E_theta) {
         return jet.e();
      } else {
         assert(false);
         return NAN;
      }
   }
   
   double angleSquared(const PseudoJet& jet1, const PseudoJet& jet2) const {
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

public:

   EnergyCorrelator(int N, double beta, Measure measure = pt_R, Strategy strategy = storage_array) : _N(N), _beta(beta), _measure(measure), _strategy(strategy) {};
   ~EnergyCorrelator(){}
   
   double result(const PseudoJet& jet) const;

};

// core EnergyCorrelator::result code in .cc file.



//------------------------------------------------------------------------
/// \class EnergyCorrelatorRatio
/// Calculates ECF(N+1,beta)/ECF(N,beta).
///
/// EnergyCorrelatorRatio(int N, double beta, Measure measure)
/// Called r_N^(beta) in the publication, equal to ECF(N+1,beta)/ECF(N,beta). 
class EnergyCorrelatorRatio : public FunctionOfPseudoJet<double> {

private:

   int _N;
   double _beta;
   EnergyCorrelator::Measure _measure;
   EnergyCorrelator::Strategy _strategy;

public:

   EnergyCorrelatorRatio(int N, double beta, EnergyCorrelator::Measure measure = EnergyCorrelator::pt_R, EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array) : _N(N), _beta(beta), _measure(measure), _strategy(strategy) {};
   ~EnergyCorrelatorRatio() {}
   
   
   double result(const PseudoJet& jet) const;

};


inline double EnergyCorrelatorRatio::result(const PseudoJet& jet) const {

   double numerator = EnergyCorrelator(_N + 1, _beta, _measure, _strategy).result(jet);
   double denominator = EnergyCorrelator(_N, _beta, _measure, _strategy).result(jet);

   return numerator/denominator;

}


//------------------------------------------------------------------------
/// \class EnergyCorrelatorDoubleRatio
/// Calculates ECF(N-1,beta)*ECP(N+1)/ECF(N,beta)^2.
///
/// EnergyCorrelatorDoubleRatio(int N, double beta, Measure measure)
/// Called C_N^(beta) in the publication, equal to r_N/r_{N-1}.
/// This is the recommended function for boosted N-prong object discrimination.
/// (N=1 for quark/gluon, N=2 for boosted W/Z/H, N=3 for boosted top)
class EnergyCorrelatorDoubleRatio : public FunctionOfPseudoJet<double> {

private:

   int _N;
   double _beta;
   EnergyCorrelator::Measure _measure;
   EnergyCorrelator::Strategy _strategy;

public:

   EnergyCorrelatorDoubleRatio(int N, double beta, EnergyCorrelator::Measure measure = EnergyCorrelator::pt_R,  EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array) : _N(N), _beta(beta), _measure(measure), _strategy(strategy) {};
   ~EnergyCorrelatorDoubleRatio() {}
   
   
   double result(const PseudoJet& jet) const;

};


inline double EnergyCorrelatorDoubleRatio::result(const PseudoJet& jet) const {

   double numerator = EnergyCorrelator(_N - 1, _beta, _measure, _strategy).result(jet) * EnergyCorrelator(_N + 1, _beta, _measure, _strategy).result(jet);
   double denominator = pow(EnergyCorrelator(_N, _beta, _measure, _strategy).result(jet), 2.0);

   return numerator/denominator;

}


} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_ENERGYCORRELATOR_HH__
