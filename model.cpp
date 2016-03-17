//  Created by peiqiw on 06/03/2016.
//  Copyright (c) 2016 peiqiw. All rights reserved.
//

#include "model.h"
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

CyberGameModel::CyberGameModel()
  : rD_(0), rU_(0), lambda_(0), bDD_(0), bUD_(0), bDU_(0), bUU_(0), \
    aD_(0), aU_(0), cE_(0), dE_(0), cH_(0), dH_(0), e0_(0), e_(0), T_(0), \
    kI_(0), kE_(0), kD_(0), kA_(0) {
}

CyberGameModel::CyberGameModel(double rD, double rU, double lambda, \
                 double bDD, double bUD, double bDU, double bUU, double aD, double aU, \
                 double cE, double dE, double cH, double dH, double e0, double e, \
                 double kI, double kE, double kD, double kA, double T)
  : rD_(rD), rU_(rU), lambda_(lambda), bDD_(bDD), bUD_(bUD), bDU_(bDU), bUU_(bUU), \
    aD_(aD), aU_(aU), cE_(cE), dE_(dE), cH_(cH), dH_(dH), e0_(e0), e_(e), T_(T), \
    kI_(kI), kE_(kE), kD_(kD), kA_(kA) {
}
  
double CyberGameModel::Lambda(double x) {
  if (x > 1.0) {
    return 1.0;
  } else if (x  < 0.0) {
    return 0.0;
  } else {
    return x;
  }
}

double CyberGameModel::Beta(double v1, double v2) {
  return Lambda(lambda_ * (v1 - v2) / kD_);
}

double CyberGameModel::AlphaH(double vH_0, double vE_0, double p1, double p2, double p3, \
                              double x1, double x2, double x3) {
  // return Lambda((dH_ * (vE_0 - vH_0) + aD_ * x3 * (p1 - p3) + \
                 aU_ * (1 - x1 - x2 - x3) * p2) / kA_);
  return (dH_ * (vE_0 - vH_0) + aD_ * x3 * (p1 - p3) + aU_ * (1 - x1 - x2 - x3) * p2) / kA_;
}

double CyberGameModel::AlphaE(double vH_0, double vE_0, double p1, double p2, double p3, \
                              double x1, double x2, double x3) {
  // return Lambda((dE_ * (vE_0 - vH_0) + aD_ * x3 * (p1 - p3) + \
                 aU_ * (1 - x1 - x2 - x3) * p2) / kA_);
  return (dE_ * (vE_0 - vH_0) + aD_ * x3 * (p1 - p3) + aU_ * (1 - x1 - x2 - x3) * p2) / kA_;
}

int CyberGameModel::Pi(double* result, double x1, double x2, double x3, \
                       double beta12, double beta34, double alpha) {
  result[0] = lambda_ * beta12 * x2 + (bDD_ * x1 + bUD_ * x2 + aD_ * alpha) * x3 - \
              (rD_ + lambda_ * (1 - beta12)) * x1;
              
  result[1] = lambda_ * (1 - beta12) * x1 - (rU_ + lambda_ * beta12) * x2 + \
              (bDU_ * x1 + bUU_ * x2 + aU_ * alpha) * (1 - x1 - x2 - x3);

  result[2] = lambda_ * beta34 * (1 - x1 - x2 - x3) + rD_ * x1 - \
              (bDD_ * x1 + bUD_ * x2 + aD_ * alpha + lambda_ * (1 - beta34)) * x3;
  
  return 0;
}




















