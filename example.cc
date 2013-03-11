// $Id$
//
// Copyright (c) -, 
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

#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
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
   
   for (int j = 0; j < 2; j++) { // Two hardest jets per event
      if (antikt_jets[j].perp() > 200) {
         PseudoJet myJet = antikt_jets[j];
         
         vector<double> betalist;
         betalist.push_back(0.1);
         betalist.push_back(0.2);
         betalist.push_back(0.5);
         betalist.push_back(1.0);
         betalist.push_back(1.5);
         betalist.push_back(2.0);
         
         cout << "-------------------------------------------------------------------------------------" << endl;
         cout << "EnergyCorrelator:  ECF(N,beta)"  << endl;
         cout << "-------------------------------------------------------------------------------------" << endl;
         printf("%7s %15s %15s %15s %15s \n","beta", "N=1 (GeV)", "N=2 (GeV^2)", "N=3 (GeV^3)", "N=4 (GeV^4)");
         
         for (unsigned int B = 0; B < betalist.size(); B++) {
            double beta = betalist[B];
            
            EnergyCorrelator ECF0(0,beta);
            EnergyCorrelator ECF1(1,beta);
            EnergyCorrelator ECF2(2,beta);
            EnergyCorrelator ECF3(3,beta);
            EnergyCorrelator ECF4(4,beta);

            printf("%7.3f %15.2f %15.2f %15.2f %15.2f \n",beta,ECF1(myJet),ECF2(myJet),ECF3(myJet),ECF4(myJet));
         }
         cout << "-------------------------------------------------------------------------------------" << endl << endl;

         cout << "-------------------------------------------------------------------------------------" << endl;
         cout << "EnergyCorrelatorRatio:  r_N^(beta) = ECF(N+1,beta)/ECF(N,beta)"  << endl;
         cout << "-------------------------------------------------------------------------------------" << endl;
         printf("%7s %15s %15s %15s %15s \n","beta", "N=0 (GeV)", "N=1 (GeV)", "N=2 (GeV)", "N=3 (GeV)");
         
         for (unsigned int B = 0; B < betalist.size(); B++) {
            double beta = betalist[B];
            
            EnergyCorrelatorRatio r0(0,beta);
            EnergyCorrelatorRatio r1(1,beta);
            EnergyCorrelatorRatio r2(2,beta);
            EnergyCorrelatorRatio r3(3,beta);

            printf("%7.3f %15.4f %15.4f %15.4f %15.4f \n",beta,r0(myJet),r1(myJet),r2(myJet),r3(myJet));
         }
         cout << "-------------------------------------------------------------------------------------" << endl << endl;

         cout << "-------------------------------------------------------------------------------------" << endl;
         cout << "EnergyCorrelatorDoubleRatio:  C_N^(beta) = r_N^(beta)/r_{N-1}^(beta)"  << endl;
         cout << "-------------------------------------------------------------------------------------" << endl;
         printf("%7s %15s %15s %15s \n","beta", "N=1", "N=2", "N=3");
         
         for (unsigned int B = 0; B < betalist.size(); B++) {
            double beta = betalist[B];
            
            EnergyCorrelatorDoubleRatio C1(1,beta);
            EnergyCorrelatorDoubleRatio C2(2,beta);
            EnergyCorrelatorDoubleRatio C3(3,beta);

            printf("%7.3f %15.6f %15.6f %15.6f \n",beta,C1(myJet),C2(myJet),C3(myJet));
         }
         cout << "-------------------------------------------------------------------------------------" << endl << endl;


      }
   }
}



