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
            

            // The angularity is set by the value of beta
            double beta;


            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorC1:  C_1^(beta) = ECF(2,beta)/ECF(1,beta)^2 " << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s \n","beta","C1 obs");


            beta = 1.0;
            EnergyCorrelatorC1 c1(beta);
            printf("%7.3f %14.6f \n",beta,c1(myJet));

            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorC2:  C_2^(beta) = ECF(3,beta)*ECF(1,beta)/ECF(2,beta)^2 " << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s \n","beta","C2 obs");


            beta = 1.0;
            EnergyCorrelatorC2 c2(beta);
            printf("%7.3f %14.6f \n",beta,c2(myJet));

            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorD2:  D_2^(beta) = ECF(3,beta)*ECF(1,beta)^3/ECF(2,beta)^3 " << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s \n","beta","D2 obs");


            beta = 1.0;
            EnergyCorrelatorD2 d2(beta);
            printf("%7.3f %14.6f \n",beta,d2(myJet));
            beta = 2.0;
            EnergyCorrelatorD2 d2_2(beta);
            printf("%7.3f %14.6f \n",beta,d2_2(myJet));


            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorGeneralizedD2:  D_2^(alpha, beta) = ECFN(3,alpha)/ECFN(2,beta)^(3*alpha/beta) " << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %20s\n","beta","alpha = 1.000");

            beta = 2.0;
            printf("%7.3f ", beta);
            double alpha = 1.0;
            EnergyCorrelatorGeneralizedD2 d2_generalized(alpha, beta);
            printf("%20.4f ", d2_generalized(myJet));
            printf("\n");

            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorN2:  N2(beta) = ECFN(3, beta, 2)/ECFN(2, beta, 1)^2 " << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s \n","beta", "N2 obs");


            beta = 1.0;
            // Directly defining the EnergyCorrelator N2
            EnergyCorrelatorN2 N2(beta);
            printf("%7.3f %14.6f \n",beta,N2(myJet));

            beta = 2.0;
            EnergyCorrelatorN2 N2_2(beta);
            printf("%7.3f %14.6f \n",beta,N2_2(myJet));

            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorN3:  N3(beta) = ECFN(4, beta, 2)/ECFN(3, beta, 1)^2 " << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s \n","beta", "N3 obs");

            beta = 1.0;
            EnergyCorrelatorN3 N3(beta);
            printf("%7.3f %14.6f \n",beta,N3(myJet));
            beta = 2.0;
            EnergyCorrelatorN3 N3_2(beta);
            printf("%7.3f %14.6f \n",beta,N3_2(myJet));


            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorM2:  M2(beta) = ECFN(3, beta, 1)/ECFN(3, beta, 1) " << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s \n","beta", "M2 obs");


            beta = 1.0;
            //Directly defining M2
            EnergyCorrelatorM2 M2(beta);
            printf("%7.3f %14.6f \n",beta,M2(myJet));

            beta = 2.0;
            EnergyCorrelatorM2 M2_2(beta);
            printf("%7.3f %14.6f \n",beta,M2_2(myJet));

            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorU1, EnergyCorrelatorU2, EnergyCorrelatorU3:  U_i(beta) = ECFN(i+1, beta, 1) " << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s %14s %14s \n","beta", "U1", "U2", "U3");


            beta = 1.0;
            //Defining the Useries for beta= 1.0
            EnergyCorrelatorU1 U1s(beta);
            EnergyCorrelatorU2 U2s(beta);
            EnergyCorrelatorU3 U3s(beta);


            printf("%7.3f %14.6f %14.6f %14.6f \n",beta,U1s(myJet),U2s(myJet),U3s(myJet));

        }
    }
}




