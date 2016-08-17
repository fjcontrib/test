#ifndef __FASTJET_CONTRIB_ENERGYCORRELATOR_HH__
#define __FASTJET_CONTRIB_ENERGYCORRELATOR_HH__

//  EnergyCorrelator Package
//  Questions/Comments?  larkoski@mit.edu gavin.salam@cern.ch jthaler@jthaler.net lnecib@mit.edu
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

#include <fastjet/internal/base.hh>
#include "fastjet/FunctionOfPseudoJet.hh"

#include <string>
#include <climits>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

/// \mainpage EnergyCorrelator contrib
///
/// The EnergyCorrelator contrib provides an implementation of energy
/// correlators and their ratios as described in arXiv:1305.0007 by
/// Larkoski, Salam and Thaler.  Additionally, the ratio observable
/// D2 described in arXiv:1409.6298 by Larkoski, Moult and Neill
/// is also included in this contrib. Finally, a normalized version of
/// the energy correlation functions is added, as well as new observables
/// N2, M2, a generalized version of D2, the MSeries and the NSeries, all
/// defined in arXiv:16xx.yyyyy by Moult, Necib and Thaler.
///
///
/// <p>There are 11 main classes:
///
/// - EnergyCorrelator
/// - EnergyCorrelatorRatio
/// - EnergyCorrelatorDoubleRatio
/// - EnergyCorrelatorD2
/// - EnergyCorrelatorNormalized
/// - EnergyCorrelatorGeneralizedD2
/// - EnergyCorrelatorNseries
/// - EnergyCorrelatorN2
/// - EnergyCorrelatorN3
/// - EnergyCorrelatorMseries
/// - EnergyCorrelatorM2
///
/// each of which is a FastJet
/// FunctionOfPseudoJet. EnergyCorrelatorDoubleRatio in particular is
/// useful for quark/gluon discrimination and boosted object tagging.
/// EnergyCorrelationD2 has been shown to be the optimal discrimination
/// observable for boosted 2-prong jets.
///
/// See the file example.cc for an illustration of usage.

//------------------------------------------------------------------------
/// \class EnergyCorrelator
/// ECF(N,beta) is the N-point energy correlation function, with an angular exponent beta.
///
/// It is defined as follows
///
///  - ECF(1,\f$ \beta)  = \sum_i E_i \f$
///  - ECF(2,\f$ \beta)  = \sum_{i<j} E_i E_j \theta_{ij}^\beta \f$
///  - ECF(3,\f$ \beta)  = \sum_{i<j<k} E_i E_j E_k (\theta_{ij} \theta_{ik} \theta_{jk})^\beta \f$
///  - ECF(4,\f$ \beta)  = \sum_{i<j<k<l} E_i E_j E_k E_l (\theta_{ij}  \theta_{ik} \theta_{il} \theta_{jk} \theta_{jl} \theta_{kl})^\beta \f$
///  - ...
///
/// The correlation can be determined with energies and angles (as
/// given above) or with transverse momenta and boost invariant angles
/// (the code's default). The choice is controlled by
/// EnergyCorrelator::Measure provided in the constructor.
///
/// The current implementation handles values of N up to and including 5.
/// Run times scale as n^N/N!, where n is the number of particles in a jet.

    class EnergyCorrelator : public FunctionOfPseudoJet<double> {
        friend class EnergyCorrelatorNormalized;
    public:

        enum Measure {
            pt_R,     ///< use transverse momenta and boost-invariant angles,
            ///< eg \f$\mathrm{ECF}(2,\beta) = \sum_{i<j} p_{ti} p_{tj} \Delta R_{ij}^{\beta} \f$
            E_theta,   ///  use energies and angles,
            ///  eg \f$\mathrm{ECF}(2,\beta) = \sum_{i<j} E_{i} E_{j}   \theta_{ij}^{\beta} \f$
            E_inv     ///  use energies and invariant mass,
            ///  eg \f$\mathrm{ECF}(2,\beta) = \sum_{i<j} E_{i} E_{j}   (\frac{2 p_{i} \cdot p_{j}}{E_{i} E_{j}})^{\beta/2} \f$
        };

        enum Strategy {
            slow,          ///< interparticle angles are not cached.
            ///< For N>=3 this leads to many expensive recomputations,
            ///< but has only O(n) memory usage for n particles

            storage_array  /// the interparticle angles are cached. This gives a significant speed
            /// improvement for N>=3, but has a memory requirement of (4n^2) bytes.
        };

    public:

        /// constructs an N-point correlator with angular exponent beta,
        /// using the specified choice of energy and angular measure as well
        /// one of two possible underlying computational Strategy
        EnergyCorrelator(int N,
                         double beta,
                         Measure measure = pt_R,
                         Strategy strategy = storage_array) :
                _N(N), _beta(beta), _measure(measure), _strategy(strategy) {};

        /// destructor
        virtual ~EnergyCorrelator(){}

        /// returns the value of the energy correlator for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

        /// returns the the part of the description related to the parameters
        std::string description_parameters() const;
        std::string description_no_N() const;

    private:

        int _N;
        double _beta;
        Measure _measure;
        Strategy _strategy;

        double energy(const PseudoJet& jet) const;
        double angleSquared(const PseudoJet& jet1, const PseudoJet& jet2) const;

    };

// core EnergyCorrelator::result code in .cc file.



//------------------------------------------------------------------------
/// \class EnergyCorrelatorRatio
/// A class to calculate the ratio of (N+1)-point to N-point energy correlators,
///     ECF(N+1,beta)/ECF(N,beta),
/// called \f$ r_N^{(\beta)} \f$ in the publication.
    class EnergyCorrelatorRatio : public FunctionOfPseudoJet<double> {

    public:

        /// constructs an (N+1)-point to N-point correlator ratio with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorRatio(int N,
                              double  beta,
                              EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                              EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _N(N), _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorRatio() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        int _N;
        double _beta;

        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };


    inline double EnergyCorrelatorRatio::result(const PseudoJet& jet) const {

        double numerator = EnergyCorrelator(_N + 1, _beta, _measure, _strategy).result(jet);
        double denominator = EnergyCorrelator(_N, _beta, _measure, _strategy).result(jet);

        return numerator/denominator;

    }


//------------------------------------------------------------------------
/// \class EnergyCorrelatorDoubleRatio
/// Calculates the double ratio of energy correlators, ECF(N-1,beta)*ECF(N+1)/ECF(N,beta)^2.
///
/// A class to calculate a double ratio of energy correlators,
///     ECF(N-1,beta)*ECF(N+1)/ECF(N,beta)^2,
/// called \f$C_N^{(\beta)}\f$ in the publication, and equal to
/// \f$ r_N^{(\beta)}/r_{N-1}^{(\beta)} \f$.
///
/// Of the different energy correlator classes, this is the one
/// recommended for quark/gluon discrimination (N=1) and for boosted
/// N-prong object discrimination (N=2 for boosted W/Z/H, N=3 for
/// boosted top).
    class EnergyCorrelatorDoubleRatio : public FunctionOfPseudoJet<double> {

    public:

        EnergyCorrelatorDoubleRatio(int N,
                                    double beta,
                                    EnergyCorrelator::Measure measure = EnergyCorrelator::pt_R,
                                    EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _N(N), _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorDoubleRatio() {}


        /// returns the value of the energy correlator double-ratio for a
        /// jet's constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        int _N;
        double _beta;
        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };


    inline double EnergyCorrelatorDoubleRatio::result(const PseudoJet& jet) const {

        double numerator = EnergyCorrelator(_N - 1, _beta, _measure, _strategy).result(jet) * EnergyCorrelator(_N + 1, _beta, _measure, _strategy).result(jet);
        double denominator = pow(EnergyCorrelator(_N, _beta, _measure, _strategy).result(jet), 2.0);

        return numerator/denominator;

    }

//------------------------------------------------------------------------
/// \class EnergyCorrelatorC1
/// A class to calculate the normalized 2-point energy correlators,
///     ECF(2,beta)/ECF(1,beta)^2,
/// called \f$ C_1^{(\beta)} \f$ in the publication.
    class EnergyCorrelatorC1 : public FunctionOfPseudoJet<double> {

    public:

        /// constructs a 2-point correlator ratio with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorC1(double  beta,
                           EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                           EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorC1() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        double _beta;

        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };


    inline double EnergyCorrelatorC1::result(const PseudoJet& jet) const {

        double numerator = EnergyCorrelator(2, _beta, _measure, _strategy).result(jet);
        double denominator = EnergyCorrelator(1, _beta, _measure, _strategy).result(jet);

        return numerator/denominator/denominator;

    }


//------------------------------------------------------------------------
/// \class EnergyCorrelatorC2
/// A class to calculate the double ratio of 3-point to 2-point
/// energy correlators,
///     ECF(3,beta)*ECF(1,beta)/ECF(2,beta)^2,
/// called \f$ C_2^{(\beta)} \f$ in the publication.
    class EnergyCorrelatorC2 : public FunctionOfPseudoJet<double> {

    public:

        /// constructs a 3-point to 2-point correlator double ratio with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorC2(double  beta,
                           EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                           EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorC2() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        double _beta;

        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };


    inline double EnergyCorrelatorC2::result(const PseudoJet& jet) const {

        double numerator3 = EnergyCorrelator(3, _beta, _measure, _strategy).result(jet);
        double numerator1 = EnergyCorrelator(1, _beta, _measure, _strategy).result(jet);
        double denominator = EnergyCorrelator(2, _beta, _measure, _strategy).result(jet);

        return numerator3*numerator1/denominator/denominator;

    }


//------------------------------------------------------------------------
/// \class EnergyCorrelatorD2
/// A class to calculate the observable formed from the ratio of the
/// 3-point and 2-point energy correlators,
///     ECF(3,beta)*ECF(1,beta)^3/ECF(2,beta)^3,
/// called \f$ D_2^{(\beta)} \f$ in the publication.
    class EnergyCorrelatorD2 : public FunctionOfPseudoJet<double> {

    public:

        /// constructs an 3-point to 2-point correlator ratio with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorD2(double  beta,
                           EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                           EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorD2() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        double _beta;

        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };


    inline double EnergyCorrelatorD2::result(const PseudoJet& jet) const {

        double numerator3 = EnergyCorrelator(3, _beta, _measure, _strategy).result(jet);
        double numerator1 = EnergyCorrelator(1, _beta, _measure, _strategy).result(jet);
        double denominator2 = EnergyCorrelator(2, _beta, _measure, _strategy).result(jet);

        return numerator3*numerator1*numerator1*numerator1/denominator2/denominator2/denominator2;

    }


//------------------------------------------------------------------------
/// \class EnergyCorrelatorNormalized
/// ECF(N,beta)/ECF(1,beta) is the normalized N-point energy correlation function, with an angular exponent beta.
/// We then generalize to i number of angles where the definition is ECFN(N, beta, angles). When angles = N choose 2,
/// this reproduces the earlier normalized energy correlation functions. By default, when angles = -1, angles is set to
/// equal N choose 2 and this reproduces the normalized version of the earlier energy correlation functions
/// When angles are not given, or given -1, the energy correlation functions are defined as follows
///
///  - ECFN(1,\f$ \beta)  = 1\f$
///  - ECFN(2,\f$ \beta)  = \sum_{i<j} z_i z_j \theta_{ij}^\beta \f$
///  - ECFN(3,\f$ \beta)  = \sum_{i<j<k} z_i z_j z_k (\theta_{ij} \theta_{ik} \theta_{jk})^\beta \f$
///  - ECFN(4,\f$ \beta)  = \sum_{i<j<k<l} z_i z_j z_k z_l (\theta_{ij}  \theta_{ik} \theta_{il} \theta_{jk} \theta_{jl} \theta_{kl})^\beta \f$
///  - ...
///  where the z_i's are the energy fractions.
///

/// When a new value of angles "a" is given, the ECFN are defined as
///  - ECFN(1,\f$ \beta, a)  = 1\f$
///  - ECFN(2,\f$ \beta, a)  = \sum_{i<j} z_i z_j \theta_{ij}^\beta \f$
///  - ECFN(3,\f$ \beta, a)  = \sum_{i<j<k} z_i z_j z_k min_a( ( \theta_{ij} \theta_{ik}), (\theta_{ik}, \theta_{jk} ),
///  (\theta_{ij}, \theta_{jk}) )^\beta \f$
///  where min_a means the product of a elements of the following list.
///  - ...

/// The correlation can be determined with energies and angles (as
/// given above) or with transverse momenta and boost invariant angles
/// (the code's default). The choice is controlled by
/// EnergyCorrelatorNormalized::Measure provided in the constructor.
///
/// The current implementation handles values of N up to and including 5.
///


    class EnergyCorrelatorNormalized : public FunctionOfPseudoJet<double> {
    public:

        /// constructs an N-point correlator with angular exponent beta,
        /// using the specified choice of energy and angular measure as well
        /// one of two possible underlying computational Strategy
        EnergyCorrelatorNormalized(int N,
                                   double beta,
                                   int angles = -1,
                                   EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                                   EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _N(N), _beta(beta), _angles(angles), _measure(measure), _strategy(strategy) {};

        /// destructor
        virtual ~EnergyCorrelatorNormalized(){}

        /// returns the value of the normalized energy correlator for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;
        std::vector<double> result_all_angles(const PseudoJet& jet) const;

    private:

        int _N;
        double _beta;
        int _angles;
        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;
        //EnergyCorrelator _helper_correlator =  EnergyCorrelator(1,1.0, _measure, _strategy);

        double energy(const PseudoJet& jet) const;
        double angleSquared(const PseudoJet& jet1, const PseudoJet& jet2) const;

    };



//------------------------------------------------------------------------
/// \class EnergyCorrelatorGeneralizedD2
/// A class to calculate the observable formed from the ratio of the
/// 3-point and 2-point energy correlators,
///     ECFN(3,alpha)/ECFN(2,beta)^3 alpha/beta,
/// called \f$ D_2^{(\alpha, \beta)} \f$ in the publication.
    class EnergyCorrelatorGeneralizedD2 : public FunctionOfPseudoJet<double> {

    public:

        /// constructs an 3-point to 2-point correlator ratio with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorGeneralizedD2(
                double alpha,
                double  beta,
                EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _alpha(alpha), _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorGeneralizedD2() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        double _alpha;
        double _beta;

        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };


    inline double EnergyCorrelatorGeneralizedD2::result(const PseudoJet& jet) const {

        double numerator = EnergyCorrelatorNormalized(3, _alpha, -1, _measure, _strategy).result(jet);
        double denominator = EnergyCorrelatorNormalized(2, _beta, -1, _measure, _strategy).result(jet);

        return numerator/pow(denominator, 3.0*_alpha/_beta);

    }


//------------------------------------------------------------------------
/// \class EnergyCorrelatorNseries
/// A class to calculate the observable formed from the ratio of the
/// 3-point and 2-point energy correlators,
///     N_n = ECFN(n+1,beta,2)/ECFN(n,beta,1)^2,
/// called \f$ N_i^{(\alpha, \beta)} \f$ in the publication.
/// By definition, N_1^\beta = ECFN(2, 2*beta, 1)
    class EnergyCorrelatorNseries : public FunctionOfPseudoJet<double> {

    public:

        /// constructs a n 3-point to 2-point correlator ratio with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorNseries(
                int n,
                double  beta,
                EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _n(n), _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorNseries() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        int _n;
        double _beta;
        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;

    };


    inline double EnergyCorrelatorNseries::result(const PseudoJet& jet) const {

        if (_n == 1) return EnergyCorrelatorNormalized(2, 2*_beta, 1, _measure, _strategy).result(jet);
        // By definition, N1 = ECFN(2)^(2 beta)
        double numerator = EnergyCorrelatorNormalized(_n + 1, _beta, 2, _measure, _strategy).result(jet);
        double denominator = EnergyCorrelatorNormalized(_n, _beta, 1, _measure, _strategy).result(jet);

        return numerator/denominator/denominator;

    }



//------------------------------------------------------------------------
/// \class EnergyCorrelatorN2
/// A class to calculate the observable formed from the ratio of the
/// 3-point and 2-point energy correlators,
///     ECFN(3,beta,2)/ECFN(2,beta,1)^2,
/// called \f$ N_2^{(\beta)} \f$ in the publication.
    class EnergyCorrelatorN2 : public FunctionOfPseudoJet<double> {

    public:

        /// constructs an 3-point to 2-point correlator ratio with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorN2(double  beta,
                           EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                           EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorN2() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        double _beta;

        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };


    inline double EnergyCorrelatorN2::result(const PseudoJet& jet) const {

        double numerator = EnergyCorrelatorNormalized(3, _beta, 2, _measure, _strategy).result(jet);
        double denominator = EnergyCorrelatorNormalized(2, _beta, 1, _measure, _strategy).result(jet);

        return numerator/denominator/denominator;

    }


//------------------------------------------------------------------------
/// \class EnergyCorrelatorN3
/// A class to calculate the observable formed from the ratio of the
/// 3-point and 2-point energy correlators,
///     ECFN(4,beta,2)/ECFN(3,beta,1)^2,
/// called \f$ N_3^{(\beta)} \f$ in the publication.
    class EnergyCorrelatorN3 : public FunctionOfPseudoJet<double> {

    public:

        /// constructs an 3-point to 2-point correlator ratio with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorN3(double  beta,
                           EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                           EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorN3() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        double _beta;

        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };


    inline double EnergyCorrelatorN3::result(const PseudoJet& jet) const {

        double numerator = EnergyCorrelatorNormalized(4, _beta, 2, _measure, _strategy).result(jet);
        double denominator = EnergyCorrelatorNormalized(3, _beta, 1, _measure, _strategy).result(jet);

        return numerator/denominator/denominator;

    }


//------------------------------------------------------------------------
/// \class EnergyCorrelatorMseries
/// A class to calculate the observable formed from the ratio of the
/// 3-point and 2-point energy correlators,
///     M_n = ECFN(n+1,beta,1)/ECFN(n,beta,1),
/// called \f$ M_i^{(\alpha, \beta)} \f$ in the publication.
/// By definition, M_1^\beta = ECFN(2, beta, 1)
    class EnergyCorrelatorMseries : public FunctionOfPseudoJet<double> {

    public:

        /// constructs a n 3-point to 2-point correlator ratio with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorMseries(
                int n,
                double  beta,
                EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _n(n), _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorMseries() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        int _n;
        double _beta;
        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;

    };


    inline double EnergyCorrelatorMseries::result(const PseudoJet& jet) const {

        if (_n == 1) return EnergyCorrelatorNormalized(2, _beta, 1, _measure, _strategy).result(jet);

        double numerator = EnergyCorrelatorNormalized(_n + 1, _beta, 1, _measure, _strategy).result(jet);
        double denominator = EnergyCorrelatorNormalized(_n, _beta, 1, _measure, _strategy).result(jet);

        return numerator/denominator;

    }

//------------------------------------------------------------------------
/// \class EnergyCorrelatorM2
/// A class to calculate the observable formed from the ratio of the
/// 3-point and 2-point energy correlators,
///     ECFN(3,beta,1)/ECFN(2,beta,1),
/// called \f$ M_2^{(\beta)} \f$ in the publication.
    class EnergyCorrelatorM2 : public FunctionOfPseudoJet<double> {

    public:

        /// constructs an 3-point to 2-point correlator ratio with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorM2(double  beta,
                           EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                           EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorM2() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        double _beta;

        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };


    inline double EnergyCorrelatorM2::result(const PseudoJet& jet) const {

        double numerator = EnergyCorrelatorNormalized(3, _beta, 1, _measure, _strategy).result(jet);
        double denominator = EnergyCorrelatorNormalized(2, _beta, 1, _measure, _strategy).result(jet);

        return numerator/denominator;

    }


} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_ENERGYCORRELATOR_HH__
