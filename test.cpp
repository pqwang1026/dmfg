//  Created by peiqiw on 09/03/2016.
//  Copyright (c) 2016 peiqiw. All rights reserved.
//

#include "modeltest.h"
#include "solvertest.h"
//#include "fdsolver.h"
//#include "model.h"
#include <assert.h>

int main() {
  ModelTest * modelTest = new ModelTest;
  modelTest->testConstructor();
  modelTest->testLambda();
  modelTest->testBeta();
  
  SolverTest * solverTest = new SolverTest;
  solverTest->testConstructor();
  solverTest->PrepareValueArray();
  solverTest->testGradient();
  solverTest->testFlipArrays();
  
  return 0;
}