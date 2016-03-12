//  Created by peiqiw on 06/03/2016.
//  Copyright (c) 2016 peiqiw. All rights reserved.
//

#ifndef cybermfg_model_h
#define cybermfg_model_h

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

using namespace std;

class CyberGameModel {
public:
  // Parameters of the model
  // recovery rate for defended/undefended user
  const double rD_;
  const double rU_;
  
  // response rate of defense setup
  const double lambda_;
  
  // contamination rate of D/U user onto D/U user
  const double bDD_;
  const double bUD_;
  const double bDU_;
  const double bUU_;
  
  // coef of effectiveness on defended/undefended user
  const double aD_;
  const double aU_;
  
  // linear coefficients for transition rate of attacker
  const double cE_;
  const double dE_;
  const double cH_;
  const double dH_;
  
  const double e0_; // income rate for attacker
  const double e_; // utility rate for user
  const double kI_; // cost of being infected
  const double kE_; // cost of being exposed
  const double kD_; // cost of security effort
  const double kA_; // cost of attack effort
  
  const double T_; // horizon of the game
  
  // Useful functions
  double Lambda(double x);

  //Constructor
  CyberGameModel();
  
  CyberGameModel(double rD, double rU, double lambda, \
                 double bDD, double bUD, double bDU, double bUU, double aD, double aU, \
                 double cE, double dE, double cH, double dH, double e0, double e, \
                 double kI, double kE, double kD, double kA, double T);
  
  //Destructor
  //~CyberGameModel();
  
  //User's optimal strategy
  double Beta(double v1, double v2);
  
  //Hacker's optimal strategy
  double AlphaH(double vH_0, double vE_0, double p1, double p2, double p3, \
                double x1, double x2, double x3);
  double AlphaE(double vH_0, double vE_0, double p1, double p2, double p3, \
                double x1, double x2, double x3);
                
  //Pi operator
  int Pi(double* result, double x1, double x2, double x3, \
         double beta12, double beta34, double alpha);
        
};

#endif










