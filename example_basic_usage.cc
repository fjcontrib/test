// Example showing usage of energy correlator classes.
//
// Compile it with "make example" and run it with
//
//   ./example_basic_usage < ../data/single-event.dat
//
// Copyright (c) 2013-2016
// Andrew Larkoski, Lina Necib, Gavin Salam, and Jesse Thaler
//
// $Id: example.cc 958 2016-08-17 00:25:14Z linoush $
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

#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
//#include <time.h>
#include <ctime>
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <string>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"

#include <sstream>
#include "EnergyCorrelator.hh" // In external code, this should be fastjet/contrib/EnergyCorrelator.hh

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);
void analyze(const vector<PseudoJet> & input_particles);

//----------------------------------------------------------------------
int main(){

    //----------------------------------------------------------
    // read in input particles
    vector<PseudoJet> event;
    read_event(event);
    cout << "# read an event with " << event.size() << " particles" << endl;

    //----------------------------------------------------------
    // illustrate how this EnergyCorrelator contrib works

    analyze(event);

    return 0;
}

// read in input particles
void read_event(vector<PseudoJet> &event){
    string line;
    while (getline(cin, line)) {
        istringstream linestream(line);
        // take substrings to avoid problems when there are extra "pollution"
        // characters (e.g. line-feed).
        if (line.substr(0,4) == "#END") {return;}
        if (line.substr(0,1) == "#") {continue;}
        double px,py,pz,E;
        linestream >> px >> py >> pz >> E;
        PseudoJet particle(px,py,pz,E);

        // push event onto back of full_event vector
        event.push_back(particle);
    }
}

////////
//
//  Main Routine for Analysis
//
///////

void analyze(const vector<PseudoJet> & input_particles) {

    /////// EnergyCorrelator /////////////////////////////

    // Initial clustering with anti-kt algorithm
    JetAlgorithm algorithm = antikt_algorithm;
    double jet_rad = 1.00; // jet radius for anti-kt algorithm
    JetDefinition jetDef = JetDefinition(algorithm,jet_rad,E_scheme,Best);
    ClusterSequence clust_seq(input_particles,jetDef);
    vector<PseudoJet> antikt_jets  = sorted_by_pt(clust_seq.inclusive_jets());

    for (int j = 0; j < 1; j++) { // Hardest jet per event
        if (antikt_jets[j].perp() > 200) {

            PseudoJet myJet = antikt_jets[j];


            EnergyCorrelator::Measure measure = EnergyCorrelator::pt_R;
            string measurename = "pt_R";

            // The angularity is set by the value of beta
            double beta;

            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelator:  ECF(N,beta) with " << measurename << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s %14s %14s %14s %15s\n","beta", "N=1 (GeV)", "N=2 (GeV^2)", "N=3 (GeV^3)", "N=4 (GeV^4)", "N=5 (GeV^5)");


            // Define the energy correlator functions
            beta = 1.0;

            EnergyCorrelator ECF0(0,beta,measure);
            EnergyCorrelator ECF1(1,beta,measure);
            EnergyCorrelator ECF2(2,beta,measure);
            EnergyCorrelator ECF3(3,beta,measure);
            EnergyCorrelator ECF4(4,beta,measure);
            EnergyCorrelator ECF5(5,beta,measure);

            // Evaluating the different ECFs.
            printf("%7.3f %14.2f %14.2f %14.2f %14.2f %15.2f \n",beta,ECF1(myJet),ECF2(myJet),ECF3(myJet),ECF4(myJet),ECF5(myJet));

            beta = 2.0;

            EnergyCorrelator ECF0_2(0,beta,measure);
            EnergyCorrelator ECF1_2(1,beta,measure);
            EnergyCorrelator ECF2_2(2,beta,measure);
            EnergyCorrelator ECF3_2(3,beta,measure);
            EnergyCorrelator ECF4_2(4,beta,measure);
            EnergyCorrelator ECF5_2(5,beta,measure);

            printf("%7.3f %14.2f %14.2f %14.2f %14.2f %15.2f \n",beta,ECF1_2(myJet),ECF2_2(myJet),ECF3_2(myJet),ECF4_2(myJet),ECF5_2(myJet));


            cout << "-------------------------------------------------------------------------------------" << endl << endl;

            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorRatio:  r_N^(beta) = ECF(N+1,beta)/ECF(N,beta) with " << measurename << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s %14s %14s %14s %15s \n","beta", "N=0 (GeV)", "N=1 (GeV)", "N=2 (GeV)", "N=3 (GeV)","N=4 (GeV)");


            beta = 1.0;

            EnergyCorrelatorRatio r0(0,beta,measure);
            EnergyCorrelatorRatio r1(1,beta,measure);
            EnergyCorrelatorRatio r2(2,beta,measure);
            EnergyCorrelatorRatio r3(3,beta,measure);
            EnergyCorrelatorRatio r4(4,beta,measure);

            printf("%7.3f %14.4f %14.4f %14.4f %14.4f %15.4f \n",beta,r0(myJet),r1(myJet),r2(myJet),r3(myJet),r4(myJet));

            cout << "-------------------------------------------------------------------------------------" << endl << endl;

            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorDoubleRatio:  C_N^(beta) = r_N^(beta)/r_{N-1}^(beta) with " << measurename << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s %14s %14s %14s \n","beta", "N=1", "N=2", "N=3", "N=4");


            beta = 1.0;

            EnergyCorrelatorDoubleRatio C1(1,beta,measure);
            EnergyCorrelatorDoubleRatio C2(2,beta,measure);
            EnergyCorrelatorDoubleRatio C3(3,beta,measure);
            EnergyCorrelatorDoubleRatio C4(4,beta,measure);

            printf("%7.3f %14.6f %14.6f %14.6f %14.6f \n",beta,C1(myJet),C2(myJet),C3(myJet),C4(myJet));

            cout << "-------------------------------------------------------------------------------------" << endl << endl;

            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorC1:  C_1^(beta) = ECF(2,beta)/ECF(1,beta)^2 with " << measurename << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s \n","beta","C1 obs");


            beta = 1.0;
            EnergyCorrelatorC1 c1(beta,measure);
            printf("%7.3f %14.6f \n",beta,c1(myJet));

            cout << "-------------------------------------------------------------------------------------" << endl << endl;

            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorC2:  C_2^(beta) = ECF(3,beta)*ECF(1,beta)/ECF(2,beta)^2 with " << measurename << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s \n","beta","C2 obs");


            beta = 1.0;
            EnergyCorrelatorC2 c2(beta,measure);
            printf("%7.3f %14.6f \n",beta,c2(myJet));
            cout << "-------------------------------------------------------------------------------------" << endl << endl;


            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorD2:  D_2^(beta) = ECF(3,beta)*ECF(1,beta)^3/ECF(2,beta)^3 with " << measurename << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s \n","beta","D2 obs");


            beta = 1.0;
            EnergyCorrelatorD2 d2(beta,measure);
            printf("%7.3f %14.6f \n",beta,d2(myJet));
            cout << "-------------------------------------------------------------------------------------" << endl << endl;


            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorGeneralizedD2:  D_2^(alpha, beta) = ECFN(3,alpha)/ECFN(2,beta)^(3*alpha/beta) with " << measurename << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %20s\n","beta","alpha = 1.000");

            beta = 2.0;
            printf("%7.3f ", beta);
            double alpha = 1.0;
            EnergyCorrelatorGeneralizedD2 d2_generalized(alpha, beta, measure);
            printf("%20.4f ", d2_generalized(myJet));
            printf("\n");
            cout << "-------------------------------------------------------------------------------------" << endl << endl;



            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorNormalized:  ECFN(N, beta, angles = N Choose 2) with " << measurename << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %7s %18s %18s %18s\n","beta", "N=1", "N=2", "N=3", "N=4");



            beta = 1.0;

            EnergyCorrelatorNormalized ECFN1(1, beta, 1, measure);
            EnergyCorrelatorNormalized ECFN2(2, beta, 1, measure);
            EnergyCorrelatorNormalized ECFN3(3, beta, 3, measure);
            EnergyCorrelatorNormalized ECFN4(4, beta, 6, measure);
            //EnergyCorrelatorNormalized ECFN5(5, 5, beta, measure);

            printf("%7.3f %7.2f %18.13f %18.13f %18.13f \n", beta, ECFN1(myJet), ECFN2(myJet), ECFN3(myJet),
                   ECFN4(myJet));

            cout << "-------------------------------------------------------------------------------------" <<
            endl << endl;


            beta = 2.0;

            EnergyCorrelatorNormalized ECFN1_2(1, beta, 1, measure);
            EnergyCorrelatorNormalized ECFN2_2(2, beta, 1, measure);
            EnergyCorrelatorNormalized ECFN3_2(3, beta, 3, measure);
            EnergyCorrelatorNormalized ECFN4_2(4, beta, 6, measure);
            //EnergyCorrelatorNormalized ECFN5_2(5, 5, beta, measure);

            printf("%7.3f %7.2f %18.13f %18.13f %18.13f \n", beta, ECFN1_2(myJet), ECFN2_2(myJet), ECFN3_2(myJet),
                   ECFN4_2(myJet));

            cout << "-------------------------------------------------------------------------------------" <<
            endl << endl;


            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorNormalized:  ECFN(N, beta=1, angles) with " << measurename << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %7s %18s %18s %18s\n","angles", "N=1", "N=2", "N=3", "N=4");

            beta = 1.0;
            for (unsigned int A = 1; A < 2; A++) {
                double angle = A;

                // Defining different normalized energy correlation functions for different number of angles
                EnergyCorrelatorNormalized ECFN1_A(1, beta, angle, measure);
                EnergyCorrelatorNormalized ECFN2_A(2, beta, angle, measure);
                EnergyCorrelatorNormalized ECFN3_A(3, beta, angle, measure);
                EnergyCorrelatorNormalized ECFN4_A(4, beta, angle, measure);

                printf("%7.0f %7.2f %18.13f %18.13f %18.13f \n", angle, ECFN1_A(myJet), ECFN2_A(myJet), ECFN3_A(myJet),
                       ECFN4_A(myJet));

            }

            for (unsigned int A = 2; A < 4; A++) {
                double angle = A;

                EnergyCorrelatorNormalized ECFN3_B(3, beta, angle, measure);
                EnergyCorrelatorNormalized ECFN4_B(4, beta, angle, measure);


                printf("%7.0f %7s %18s %18.13f %18.13f \n", angle, " " , " " , ECFN3_B(myJet),
                       ECFN4_B(myJet));
            }

            for (unsigned int A = 4; A < 7; A++) {
                double angle = A;

                EnergyCorrelatorNormalized ECFN4_C(4, beta, angle, measure);

                printf("%7.0f %7s %18s %18s %18.13f \n", angle, " ", " ", " ", ECFN4_C(myJet));
            }
            cout << "-------------------------------------------------------------------------------------" <<
            endl << endl;

            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorNseries:  N_i(beta) = ECFN(i+1, beta, 2)/ECFN(i, beta, 1)^2 with " << measurename << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s %14s %14s \n","beta", "N=1", "N=2", "N=3");


            beta = 1.0;
            // Defining the N series for beta = 1.0
            EnergyCorrelatorNseries N1s(1,beta,measure);
            EnergyCorrelatorNseries N2s(2,beta,measure);
            EnergyCorrelatorNseries N3s(3,beta,measure);


            printf("%7.3f %14.6f %14.6f %14.6f \n",beta,N1s(myJet),N2s(myJet),N3s(myJet));

            cout << "-------------------------------------------------------------------------------------" << endl << endl;

            beta = 2.0;
            //Defining the N series for beta = 2.0
            EnergyCorrelatorNseries N1s_2(1,beta,measure);
            EnergyCorrelatorNseries N2s_2(2,beta,measure);
            EnergyCorrelatorNseries N3s_2(3,beta,measure);


            printf("%7.3f %14.6f %14.6f %14.6f \n",beta,N1s_2(myJet),N2s_2(myJet),N3s_2(myJet));

            cout << "-------------------------------------------------------------------------------------" << endl << endl;


            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorN2:  N2(beta) = ECFN(3, beta, 2)/ECFN(2, beta, 1)^2 with " << measurename << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s \n","beta", "N2 obs");


            beta = 1.0;
            // Directly defining the EnergyCorrelator N2
            EnergyCorrelatorN2 N2(beta,measure);
            printf("%7.3f %14.6f \n",beta,N2(myJet));

            beta = 2.0;
            EnergyCorrelatorN2 N2_2(beta,measure);
            printf("%7.3f %14.6f \n",beta,N2_2(myJet));

            cout << "-------------------------------------------------------------------------------------" << endl << endl;


            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorN3:  N3(beta) = ECFN(4, beta, 2)/ECFN(3, beta, 1)^2 with " << measurename << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s \n","beta", "N3 obs");

            beta = 1.0;
            EnergyCorrelatorN3 N3(beta,measure);
            printf("%7.3f %14.6f \n",beta,N3(myJet));
            cout << "-------------------------------------------------------------------------------------" << endl << endl;

            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorMseries:  M_i(beta) = ECFN(i+1, beta, 1)/ECFN(i, beta, 1) with " << measurename << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s %14s %14s \n","beta", "N=1", "N=2", "N=3");


            beta = 1.0;
            //Defining the Mseries for beta= 1.0
            EnergyCorrelatorMseries M1s(1,beta,measure);
            EnergyCorrelatorMseries M2s(2,beta,measure);
            EnergyCorrelatorMseries M3s(3,beta,measure);


            printf("%7.3f %14.6f %14.6f %14.6f \n",beta,M1s(myJet),M2s(myJet),M3s(myJet));

            beta = 2.0;
            //Defining the Mseries for beta= 2.0
            EnergyCorrelatorMseries M1s_2(1,beta,measure);
            EnergyCorrelatorMseries M2s_2(2,beta,measure);
            EnergyCorrelatorMseries M3s_2(3,beta,measure);

            printf("%7.3f %14.6f %14.6f %14.6f \n",beta,M1s_2(myJet),M2s_2(myJet),M3s_2(myJet));

            cout << "-------------------------------------------------------------------------------------" << endl << endl;


            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorM2:  M2(beta) = ECFN(3, beta, 1)/ECFN(3, beta, 1) with " << measurename << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s \n","beta", "M2 obs");


            beta = 1.0;
            //Directly defining M2
            EnergyCorrelatorM2 M2(beta,measure);
            printf("%7.3f %14.6f \n",beta,M2(myJet));

            beta = 2.0;
            EnergyCorrelatorM2 M2_2(beta,measure);
            printf("%7.3f %14.6f \n",beta,M2_2(myJet));
            cout << "-------------------------------------------------------------------------------------" << endl << endl;

            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorUseries:  U_i(beta) = ECFN(i+1, beta, 1) with " << measurename << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s %14s %14s \n","beta", "N=1", "N=2", "N=3");


            beta = 1.0;
            //Defining the Useries for beta= 1.0
            EnergyCorrelatorUseries U1s(1,beta,measure);
            EnergyCorrelatorUseries U2s(2,beta,measure);
            EnergyCorrelatorUseries U3s(3,beta,measure);


            printf("%7.3f %14.6f %14.6f %14.6f \n",beta,U1s(myJet),U2s(myJet),U3s(myJet));

        }
    }
}




