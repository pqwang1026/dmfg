//  Created by peiqiw on 10/03/2016.
//  Copyright (c) 2016 peiqiw. All rights reserved.
//


#include "fdsolver.h"
#include "model.h"
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

int main() {
  const double rD = 4.0;
  const double rU = 2.0;
  
  // response rate of defense setup
  const double lambda = 1.0;
  
  // contamination rate of D/U user onto D/U user
  const double bDD = 0.5;
  const double bUD = 0.6;
  const double bDU = 1.5;
  const double bUU = 1.6;
  
  // coef of effectiveness on defended/undefended user
  const double aD = 0.08;
  const double aU = 0.06;
  
  // linear coefficients for transition rate of attacker
  const double cE = 0.25;
  const double dE = 0.01;
  const double cH = 0.05;
  const double dH = 0.005;
  
  const double e0 = 0.1; // income rate for attacker
  const double e = 0.1;// utility rate for user
  const double kI = 0.25; // cost of being infected
  const double kE = 0.25;  // cost of being exposed
  const double kD = 0.1; // cost of security effort
  const double kA = 0.1; // cost of attack effort
  
  const double T = 1.0; // horizon of the game
  
  CyberGameModel* Model = new CyberGameModel(rD, rU, lambda, bDD, bUD, bDU, bUU, aD, aU, \
                                 cE, dE, cH, dH, e0, e, kI, kE, kD, kA, T);                               
  const int nt = 100;
  const int nx = 100;
  
  FiniteDiffSolver* Solver = new FiniteDiffSolver(nx, nt, *Model);
  for (int t = 0; t < 30; t++) {
    Solver->Step();
  }
  //Solver->FlipArrays();
  //Solver->GradientInnerPoint(1, 1, 1);
  //printf("%15.8f %15.8f %15.8f \n", Solver->gE0[0], Solver->gE0[1], Solver->gE0[2]);
  //printf("%15.8f %15.8f \n", Solver->vE0[1][Solver->Index(1,1,1)], Solver->vH0[1][Solver->Index(1,1,1)]);
  //Solver->ComputeOptimalStrategy(1, Solver->Index(1,1,1), 0.01, 0.01, 0.01);
  //double a = Model->AlphaE(Solver->vH0[1][Solver->Index(1,1,1)], Solver->vE0[1][Solver->Index(1,1,1)], \
  //Solver->gE0[0], Solver->gE0[1], Solver->gE0[2], 0.01, 0.01, 0.01);
  //printf("%15.8f \n", a);
  //printf("%15.8f %15.8f \n", Solver->alphaE[1][Solver->Index(1,1,1)], Solver->alphaH[1][Solver->Index(1,1,1)]);
  
  Solver->PrintStrategy();
  //printf("%4.8f \n", Solver->beta34H[0][Solver->Index(50,0,0)]);
  
  delete Solver;
  delete Model;
  
  return 0;
}