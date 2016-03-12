//  Created by peiqiw on 09/03/2016.
//  Copyright (c) 2016 peiqiw. All rights reserved.
//

#include "model.h"
#include "fdsolver.h"
#include "solvertest.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

SolverTest::SolverTest() {
  // Parameters of the model
  // recovery rate for defended/undefended user
  const double rD = 4.0;
  const double rU = 2.0;
  
  // response rate of defense setup
  const double lambda = 10.0;
  
  // contamination rate of D/U user onto D/U user
  const double bDD = 5.0;
  const double bUD = 6.0;
  const double bDU = 15.0;
  const double bUU = 16.0;
  
  // coef of effectiveness on defended/undefended user
  const double aD = 12.0;
  const double aU = 6.0;
  
  // linear coefficients for transition rate of attacker
  const double cE = 5.0;
  const double dE = 5.0;
  const double cH = 10.0;
  const double dH = 10.0;
  
  const double e0 = 1.0; // income rate for attacker
  const double e = 1.0;// utility rate for user
  const double kI = 0.25; // cost of being infected
  const double kE = 0.25;  // cost of being exposed
  const double kD = 1.0; // cost of security effort
  const double kA = 1.0; // cost of attack effort
  
  const double T = 5.0; // horizon of the game
  
  testModel = new CyberGameModel(rD, rU, lambda, bDD, bUD, bDU, bUU, aD, aU, \
                                 cE, dE, cH, dH, e0, e, kI, kE, kD, kA, T);
}

int SolverTest::AlmostEqual(double a, double b) {
  assert(abs(a - b) < 0.00001);
  return 0;
}

int SolverTest::testConstructor() {
  const int nt = 3;
  const int nx = 3;
  testSolver = new FiniteDiffSolver(nx, nt, *testModel);
  assert(testSolver->nx_ == 3);
  assert(testSolver->nt_ == 3);
  assert(testSolver->dx_ == 1.0/3.0);
  assert(testSolver->dt_ == 5.0/3.0);
  assert(testSolver->model_.Lambda(0.5) == 0.5);
  assert(0 == testSolver->vH0[3][testSolver->Index(0,0,3)]);
  assert(0 == testSolver->vH0[0][testSolver->Index(0,3,0)]);
  assert(0 == testSolver->vH0[0][testSolver->Index(3,0,0)]);
  assert(0 == testSolver->vH0[1][testSolver->Index(0,2,1)]);
  assert(0 == testSolver->vH0[2][testSolver->Index(0,1,2)]);
  return 0;
}

int SolverTest::PrepareValueArray() {
  for (int k = 0; k <=3; k++) {
    for (int j = 0; j <= 3 - k; j++) {
      for (int i = 0; i <= 3 - k - j; i++) {
        testSolver->vH0_[k][testSolver->Index(i, j, k)] = (i + 1) * (j + 2) * (k + 3) / 3.0;
      }
    }
  }
  AlmostEqual(testSolver->vH0_[0][testSolver->Index(0, 0, 0)], 1 * 2 * 3 / 3.0);
  AlmostEqual(testSolver->vH0_[0][testSolver->Index(1, 0, 0)], 2 * 2 * 3 / 3.0);
  return 0;
}

int SolverTest::testGradient() {
  testSolver->GradientInnerPoint(0, 0, 0);
  // printf(" %4.8f", testSolver->gH0[2]);
  AlmostEqual(testSolver->gH0[0], (2 * 2 * 3 - 1 * 2 * 3));
  AlmostEqual(testSolver->gH0[1], (1 * 3 * 3 - 1 * 2 * 3));
  AlmostEqual(testSolver->gH0[2], (1 * 2 * 4 - 1 * 2 * 3));
  
  testSolver->GradientSurfacePoint(1, 1, 1);
  AlmostEqual(testSolver->gH0[0], (2 * 3 * 4 - 1 * 3 * 4));
  AlmostEqual(testSolver->gH0[1], (2 * 3 * 4 - 2 * 2 * 4));
  AlmostEqual(testSolver->gH0[2], (2 * 3 * 4 - 2 * 3 * 3));
  
  testSolver->GradientEdgePoint1(0, 2, 1);
  AlmostEqual(testSolver->gH0[0], 12.0);
  AlmostEqual(testSolver->gH0[1], 4.0);
  AlmostEqual(testSolver->gH0[2], 4.0);
  
  testSolver->GradientEdgePoint2(2, 0, 1);
  AlmostEqual(testSolver->gH0[0], 8.0);
  AlmostEqual(testSolver->gH0[1], 8.0);
  AlmostEqual(testSolver->gH0[2], 6.0);
  
  testSolver->GradientEdgePoint3(2, 1, 0);
  AlmostEqual(testSolver->gH0[0], 9.0);
  AlmostEqual(testSolver->gH0[1], 9.0);
  AlmostEqual(testSolver->gH0[2], 6.0);
  
  testSolver->GradientVertexPoint1();
  AlmostEqual(testSolver->gH0[0], 6.0);
  AlmostEqual(testSolver->gH0[1], 9.0);
  AlmostEqual(testSolver->gH0[2], 6.0);
  
  testSolver->GradientVertexPoint2();
  AlmostEqual(testSolver->gH0[0], 12.0);
  AlmostEqual(testSolver->gH0[1], 3.0);
  AlmostEqual(testSolver->gH0[2], 4.0);
  
  testSolver->GradientVertexPoint3();
  AlmostEqual(testSolver->gH0[0], 10.0);
  AlmostEqual(testSolver->gH0[1], 5.0);
  AlmostEqual(testSolver->gH0[2], 2.0);
}

int SolverTest::testFlipArrays() {
  PrepareValueArray();
  testSolver->FlipArrays();
  AlmostEqual(testSolver->vH0[0][testSolver->Index(0, 0, 0)], 1 * 2 * 3 / 3.0);
  AlmostEqual(testSolver->vH0_[0][testSolver->Index(0, 0, 0)], 0);
}