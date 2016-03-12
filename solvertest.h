//  Created by peiqiw on 09/03/2016.
//  Copyright (c) 2016 peiqiw. All rights reserved.
//

#ifndef solvertest_h
#define solvertest_h

#include "model.h"
#include "fdsolver.h"
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

class SolverTest {
private:
  FiniteDiffSolver* testSolver;
  CyberGameModel* testModel;
public:
  SolverTest();
  int AlmostEqual(double a, double b);
  int testConstructor();
  int PrepareValueArray();
  int testGradient();
  int testFlipArrays();
};

#endif