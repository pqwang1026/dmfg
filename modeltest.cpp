//  Created by peiqiw on 09/03/2016.
//  Copyright (c) 2016 peiqiw. All rights reserved.
//

#include "modeltest.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

int ModelTest::testConstructor() {
// Parameters of the model
  // recovery rate for defended/undefended user
  const double rD = 4;
  const double rU = 2;
  
  // response rate of defense setup
  const double lambda = 10;
  
  // contamination rate of D/U user onto D/U user
  const double bDD = 5;
  const double bUD = 6;
  const double bDU = 15;
  const double bUU = 16;
  
  // coef of effectiveness on defended/undefended user
  const double aD = 12;
  const double aU = 6;
  
  // linear coefficients for transition rate of attacker
  const double cE = 5;
  const double dE = 5;
  const double cH = 10;
  const double dH = 10;
  
  const double e0 = 1; // income rate for attacker
  const double e = 1;// utility rate for user
  const double kI = 0.25; // cost of being infected
  const double kE = 0.25;  // cost of being exposed
  const double kD = 1; // cost of security effort
  const double kA = 1; // cost of attack effort
  
  const double T = 5; // horizon of the game
  
  testModel = new CyberGameModel(rD, rU, lambda, bDD, bUD, bDU, bUU, aD, aU, \
                                 cE, dE, cH, dH, e0, e, kI, kE, kD, kA, T);
  return 0;
}

int ModelTest::testLambda() {
  assert(testModel->Lambda(0.5) == 0.5);
  assert(testModel->Lambda(-0.1) == 0);
  assert(testModel->Lambda(1.2) == 1);
}

int ModelTest::testBeta() {
  assert(testModel->Beta(0, 1) == 0);
  assert(testModel->Beta(1, 0) == 1);
  assert(testModel->Beta(0.05, 0) == 0.5);
}
