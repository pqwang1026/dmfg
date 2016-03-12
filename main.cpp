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
  const double aD = 1.2;
  const double aU = 0.6;
  
  // linear coefficients for transition rate of attacker
  const double cE = 0.5;
  const double dE = 0.5;
  const double cH = 1.0;
  const double dH = 1.0;
  
  const double e0 = 1.0; // income rate for attacker
  const double e = 1.0;// utility rate for user
  const double kI = 0.25; // cost of being infected
  const double kE = 0.25;  // cost of being exposed
  const double kD = 1.0; // cost of security effort
  const double kA = 1.0; // cost of attack effort
  
  const double T = 1.0; // horizon of the game
  
  CyberGameModel* Model = new CyberGameModel(rD, rU, lambda, bDD, bUD, bDU, bUU, aD, aU, \
                                 cE, dE, cH, dH, e0, e, kI, kE, kD, kA, T);                               
  const int nt = 100;
  const int nx = 100;
  
  FiniteDiffSolver* Solver = new FiniteDiffSolver(nx, nt, *Model);
  for (int t = 0; t < 50; t++) {
    Solver->Step();
  }
  printf("%4.8f \n", Solver->v1E[0][Solver->Index(100,0,0)]);
  
  delete Solver;
  delete Model;
  
  return 0;
}