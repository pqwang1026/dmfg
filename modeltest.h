//  Created by peiqiw on 09/03/2016.
//  Copyright (c) 2016 peiqiw. All rights reserved.
//

#ifndef modeltest_h
#define modeltest_h

#include "model.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

class ModelTest {
private:
  CyberGameModel * testModel;

public:
  int testConstructor();
  int testLambda();
  int testBeta();
  int testPi();
};

#endif