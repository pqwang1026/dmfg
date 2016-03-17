//  Created by peiqiw on 06/03/2016.
//  Copyright (c) 2016 peiqiw. All rights reserved.
//

#include "model.h"
#include "fdsolver.h"
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

FiniteDiffSolver::FiniteDiffSolver(int nx, int nt, CyberGameModel &model)
  : nx_(nx), nt_(nt), model_(model), dx_(1.0/nx), dt_(model.T_/nt) {
  int n;
  //initialize data structure
  n = nx_ + 1;
  vH0 = new double*[n];  vH0_ = new double*[n];
  vE0 = new double*[n];  vE0_ = new double*[n];
  v1H = new double*[n];  v1H_ = new double*[n];
  v2H = new double*[n];  v2H_ = new double*[n];
  v3H = new double*[n];  v3H_ = new double*[n];
  v4H = new double*[n];  v4H_ = new double*[n];
  v1E = new double*[n];  v1E_ = new double*[n];
  v2E = new double*[n];  v2E_ = new double*[n];
  v3E = new double*[n];  v3E_ = new double*[n];
  v4E = new double*[n];  v4E_ = new double*[n];
  alphaH = new double*[n];  alphaE = new double*[n];
  beta12H = new double*[n];  beta34H = new double*[n];
  beta12E = new double*[n];  beta34E = new double*[n];
  for (int k = 0; k < nx_ + 1; k++){
    n = (nx_ + 1 - k) * (nx_ + 2 - k) / 2;
    vH0[k] = new double[n];  vH0_[k] = new double[n];
    vE0[k] = new double[n];  vE0_[k] = new double[n];
    v1H[k] = new double[n];  v1H_[k] = new double[n];
    v2H[k] = new double[n];  v2H_[k] = new double[n];
    v3H[k] = new double[n];  v3H_[k] = new double[n];
    v4H[k] = new double[n];  v4H_[k] = new double[n];
    v1E[k] = new double[n];  v1E_[k] = new double[n];
    v2E[k] = new double[n];  v2E_[k] = new double[n];
    v3E[k] = new double[n];  v3E_[k] = new double[n];
    v4E[k] = new double[n];  v4E_[k] = new double[n];
    alphaH[k] = new double[n];  alphaE[k] = new double[n];
    beta12H[k] = new double[n];  beta34H[k] = new double[n];
    beta12E[k] = new double[n];  beta34E[k] = new double[n];
    for (int i = 0; i < n; i++) {
      vH0[k][i] = 0.0;  vE0[k][i] = 0.0;
      v1H[k][i] = 0.0;  v2H[k][i] = 0.0;  v3H[k][i] = 0.0;  v4H[k][i] = 0.0;
      v1E[k][i] = 0.0;  v2E[k][i] = 0.0;  v3E[k][i] = 0.0;  v4E[k][i] = 0.0;
      vH0_[k][i] = 0.0;  vE0_[k][i] = 0.0;
      v1H_[k][i] = 0.0;  v2H_[k][i] = 0.0;  v3H_[k][i] = 0.0;  v4H_[k][i] = 0.0;
      v1E_[k][i] = 0.0;  v2E_[k][i] = 0.0;  v3E_[k][i] = 0.0;  v4E_[k][i] = 0.0;
    }
  }
  piE = new double[3];  piH = new double[3];
  gH0 = new double[3];  gE0 = new double[3];
  g1H = new double[3];  g2H = new double[3];  g3H = new double[3];  g4H = new double[3];
  g1E = new double[3];  g2E = new double[3];  g3E = new double[3];  g4E = new double[3];
  
  //Open file for recording the optimal strategy
  strFileBeta12H = fopen("/Users/peiqiw/Documents/DMFG/data/beta12H.txt", "w");
  strFileBeta34H = fopen("/Users/peiqiw/Documents/DMFG/data/beta34H.txt", "w");
  strFileBeta12E = fopen("/Users/peiqiw/Documents/DMFG/data/beta12E.txt", "w");
  strFileBeta34E = fopen("/Users/peiqiw/Documents/DMFG/data/beta34E.txt", "w");
  strFileAlphaH = fopen("/Users/peiqiw/Documents/DMFG/data/alphaH.txt", "w");
  strFileAlphaE = fopen("/Users/peiqiw/Documents/DMFG/data/alphaE.txt", "w");
}

FiniteDiffSolver::~FiniteDiffSolver() {
  for (int k = 0; k < nx_ + 1; k++){
    delete [] vH0[k]; delete [] vE0[k];
    delete [] v1H[k]; delete [] v2H[k]; delete [] v3H[k]; delete [] v4H[k];
    delete [] v1E[k]; delete [] v2E[k]; delete [] v3E[k]; delete [] v4E[k];
  }
  delete [] vH0; delete [] vE0;
  delete [] v1H; delete [] v2H; delete [] v3H; delete [] v4H;
  delete [] v1E; delete [] v2E; delete [] v3E; delete [] v4E;
  delete [] piE; delete [] piH; delete [] gH0; delete [] gE0;
  delete [] g1H; delete [] g2H; delete [] g3H; delete [] g4H;
  delete [] g1E; delete [] g2E; delete [] g3E; delete [] g4E;
  
  fclose(strFileBeta12H); fclose(strFileBeta34H);
  fclose(strFileBeta12E); fclose(strFileBeta34E);
  fclose(strFileAlphaH); fclose(strFileAlphaE);
}

int FiniteDiffSolver::Index(int i, int j, int k) {
  return (i + j * (2 * nx_ + 3 - 2 * k - j) / 2);
}

int FiniteDiffSolver::FlipOneArray(double** a, double** b) {
  double **auxptr = a;
  a = b;
  b = auxptr;
  return 0;
}

int FiniteDiffSolver::FlipArrays() {
  double **auxptr;
  auxptr = vH0; vH0 = vH0_; vH0_ = auxptr;    auxptr = vE0; vE0 = vE0_; vE0_ = auxptr;
  auxptr = v1H; v1H = v1H_; v1H_ = auxptr;    auxptr = v2H; v2H = v2H_; v2H_ = auxptr;
  auxptr = v3H; v3H = v3H_; v3H_ = auxptr;    auxptr = v4H; v4H = v4H_; v4H_ = auxptr;
  auxptr = v1E; v1E = v1E_; v1E_ = auxptr;    auxptr = v2E; v2E = v2E_; v2E_ = auxptr;
  auxptr = v3E; v3E = v3E_; v3E_ = auxptr;    auxptr = v4E; v4E = v4E_; v4E_ = auxptr;
  return 0;
}

// Compute the gradient of an inner point (i,j,k) of function v, write result.
int FiniteDiffSolver::GradientInnerPoint(double* result, double** v, int k, \
          int thisIndex, int voisinIndex1, int voisinIndex2, int voisinIndex3) 
{
  //voisin1 = (i+1, j, k)
  //voisin2 = (i, j+1, k)
  //voisin3 = (i, j, k+1)
  result[0] = (v[k][voisinIndex1] - v[k][thisIndex])/dx_;
  result[1] = (v[k][voisinIndex2] - v[k][thisIndex])/dx_;
  result[2] = (v[k+1][voisinIndex3] - v[k][thisIndex])/dx_;
  return 0;
}

int FiniteDiffSolver::GradientInnerPoint(int i, int j, int k) {
  //voisin1 = (i+1, j, k)
  //voisin2 = (i, j+1, k)
  //voisin3 = (i, j, k+1)
  int a1 = Index(i, j, k);
  int a2 = Index(i + 1, j, k);
  int a3 = Index(i, j + 1, k);
  int a4 = Index(i, j, k + 1);
  GradientInnerPoint(gH0, vH0_, k, a1, a2, a3, a4);  
  GradientInnerPoint(gE0, vH0_, k, a1, a2, a3, a4);
  GradientInnerPoint(g1E, v1H_, k, a1, a2, a3, a4); 
  GradientInnerPoint(g2E, v2H_, k, a1, a2, a3, a4);
  GradientInnerPoint(g1E, v3H_, k, a1, a2, a3, a4); 
  GradientInnerPoint(g2E, v4H_, k, a1, a2, a3, a4);
  GradientInnerPoint(g1E, v1E_, k, a1, a2, a3, a4); 
  GradientInnerPoint(g2E, v2E_, k, a1, a2, a3, a4);
  GradientInnerPoint(g1E, v3E_, k, a1, a2, a3, a4);  
  GradientInnerPoint(g2E, v4E_, k, a1, a2, a3, a4);
  return 0;
}

// Compute the gradient of a surface point (i,j,k) of function v, write result.
int FiniteDiffSolver::GradientSurfacePoint(double* result, double** v, int k, \
                         int thisIndex, int voisinIndex1, int voisinIndex2, int voisinIndex3) 
{
  //voisin1 = (i-1, j, k)
  //voisin2 = (i, j-1, k)
  //voisin3 = (i, j, k-1)
  result[0] = (v[k][thisIndex] - v[k][voisinIndex1])/dx_;
  result[1] = (v[k][thisIndex] - v[k][voisinIndex2])/dx_;
  result[2] = (v[k][thisIndex] - v[k - 1][voisinIndex3])/dx_;
  return 0;
}

int FiniteDiffSolver::GradientSurfacePoint(int i, int j, int k) {
  //voisin1 = (i-1, j, k)
  //voisin2 = (i, j-1, k)
  //voisin3 = (i, j, k-1)
  int a1 = Index(i, j, k);
  int a2 = Index(i - 1, j, k);
  int a3 = Index(i, j - 1, k);
  int a4 = Index(i, j, k - 1);
  GradientSurfacePoint(gH0, vH0_, k, a1, a2, a3, a4);  
  GradientSurfacePoint(gE0, vH0_, k, a1, a2, a3, a4);
  GradientSurfacePoint(g1E, v1H_, k, a1, a2, a3, a4); 
  GradientSurfacePoint(g2E, v2H_, k, a1, a2, a3, a4);
  GradientSurfacePoint(g1E, v3H_, k, a1, a2, a3, a4); 
  GradientSurfacePoint(g2E, v4H_, k, a1, a2, a3, a4);
  GradientSurfacePoint(g1E, v1E_, k, a1, a2, a3, a4); 
  GradientSurfacePoint(g2E, v2E_, k, a1, a2, a3, a4);
  GradientSurfacePoint(g1E, v3E_, k, a1, a2, a3, a4);  
  GradientSurfacePoint(g2E, v4E_, k, a1, a2, a3, a4);
  return 0;
}

// Compute the gradient of a edge point (i,j,k) x1+x2+x3=1, x1=0
int FiniteDiffSolver::GradientEdgePoint1(double* result, double** v, int k, \
                         int thisIndex, int voisinIndex1, int voisinIndex2, int voisinIndex3)
{
  //voisin1 = (i+1, j-1, k)
  //voisin2 = (i, j-1, k)
  //voisin3 = (i, j, k-1)
  result[0] = (v[k][voisinIndex1] - v[k][voisinIndex2])/dx_;
  result[1] = (v[k][thisIndex] - v[k][voisinIndex2])/dx_;
  result[2] = (v[k][thisIndex] - v[k - 1][voisinIndex3])/dx_;
  return 0;
}

int FiniteDiffSolver::GradientEdgePoint1(int i, int j, int k) {
  //voisin1 = (i+1, j-1, k)
  //voisin2 = (i, j-1, k)
  //voisin3 = (i, j, k-1)
  int a1 = Index(i, j, k);
  int a2 = Index(i + 1, j - 1, k);
  int a3 = Index(i, j - 1, k);
  int a4 = Index(i, j, k - 1);
  GradientEdgePoint1(gH0, vH0_, k, a1, a2, a3, a4);  
  GradientEdgePoint1(gE0, vH0_, k, a1, a2, a3, a4);
  GradientEdgePoint1(g1E, v1H_, k, a1, a2, a3, a4); 
  GradientEdgePoint1(g2E, v2H_, k, a1, a2, a3, a4);
  GradientEdgePoint1(g1E, v3H_, k, a1, a2, a3, a4); 
  GradientEdgePoint1(g2E, v4H_, k, a1, a2, a3, a4);
  GradientEdgePoint1(g1E, v1E_, k, a1, a2, a3, a4); 
  GradientEdgePoint1(g2E, v2E_, k, a1, a2, a3, a4);
  GradientEdgePoint1(g1E, v3E_, k, a1, a2, a3, a4);  
  GradientEdgePoint1(g2E, v4E_, k, a1, a2, a3, a4);
  return 0;
}

// Compute the gradient of a edge point (i,j,k) x1+x2+x3=1, x2=0
int FiniteDiffSolver::GradientEdgePoint2(double* result, double** v, int k, \
                         int thisIndex, int voisinIndex1, int voisinIndex2, int voisinIndex3)
{
  //voisin1 = (i-1, j, k)
  //voisin2 = (i-1, j+1, k)
  //voisin3 = (i, j, k-1)
  result[0] = (v[k][thisIndex] - v[k][voisinIndex1])/dx_;
  result[1] = (v[k][voisinIndex2] - v[k][voisinIndex1])/dx_;
  result[2] = (v[k][thisIndex] - v[k - 1][voisinIndex3])/dx_;
  return 0;
}

int FiniteDiffSolver::GradientEdgePoint2(int i, int j, int k) {
  //voisin1 = (i-1, j, k)
  //voisin2 = (i-1, j+1, k)
  //voisin3 = (i, j, k-1)
  int a1 = Index(i, j, k);
  int a2 = Index(i - 1, j, k);
  int a3 = Index(i - 1, j + 1, k);
  int a4 = Index(i, j, k - 1);
  GradientEdgePoint2(gH0, vH0_, k, a1, a2, a3, a4);  
  GradientEdgePoint2(gE0, vH0_, k, a1, a2, a3, a4);
  GradientEdgePoint2(g1E, v1H_, k, a1, a2, a3, a4); 
  GradientEdgePoint2(g2E, v2H_, k, a1, a2, a3, a4);
  GradientEdgePoint2(g1E, v3H_, k, a1, a2, a3, a4); 
  GradientEdgePoint2(g2E, v4H_, k, a1, a2, a3, a4);
  GradientEdgePoint2(g1E, v1E_, k, a1, a2, a3, a4); 
  GradientEdgePoint2(g2E, v2E_, k, a1, a2, a3, a4);
  GradientEdgePoint2(g1E, v3E_, k, a1, a2, a3, a4);  
  GradientEdgePoint2(g2E, v4E_, k, a1, a2, a3, a4);
  return 0;
}

// Compute the gradient of a edge point (i,j,k) x1+x2+x3=1, x3=0
int FiniteDiffSolver::GradientEdgePoint3(double* result, double** v, int k, \
                         int thisIndex, int voisinIndex1, int voisinIndex2, int voisinIndex3)
{
  //voisin1 = (i-1, j, k)
  //voisin2 = (i, j-1, k)
  //voisin3 = (i-1, j, k+1)
  result[0] = (v[k][thisIndex] - v[k][voisinIndex1])/dx_;
  result[1] = (v[k][thisIndex] - v[k][voisinIndex2])/dx_;
  result[2] = (v[k + 1][voisinIndex3] - v[k][voisinIndex1])/dx_;
  return 0;
}

int FiniteDiffSolver::GradientEdgePoint3(int i, int j, int k) {
  //voisin1 = (i-1, j, k)
  //voisin2 = (i, j-1, k)
  //voisin3 = (i-1, j, k+1)
  int a1 = Index(i, j, k);
  int a2 = Index(i - 1, j, k);
  int a3 = Index(i, j - 1, k);
  int a4 = Index(i - 1, j, k + 1);
  GradientEdgePoint3(gH0, vH0_, k, a1, a2, a3, a4);  
  GradientEdgePoint3(gE0, vH0_, k, a1, a2, a3, a4);
  GradientEdgePoint3(g1E, v1H_, k, a1, a2, a3, a4); 
  GradientEdgePoint3(g2E, v2H_, k, a1, a2, a3, a4);
  GradientEdgePoint3(g1E, v3H_, k, a1, a2, a3, a4); 
  GradientEdgePoint3(g2E, v4H_, k, a1, a2, a3, a4);
  GradientEdgePoint3(g1E, v1E_, k, a1, a2, a3, a4); 
  GradientEdgePoint3(g2E, v2E_, k, a1, a2, a3, a4);
  GradientEdgePoint3(g1E, v3E_, k, a1, a2, a3, a4);  
  GradientEdgePoint3(g2E, v4E_, k, a1, a2, a3, a4);
  return 0;
}

// Compute the gradient of V1 (1, 0, 0)
int FiniteDiffSolver::GradientVertexPoint1(double* result, double** v, \
                         int thisIndex, int voisinIndex1, int voisinIndex2, int voisinIndex3)
{
  //voisin1 = (nx_ - 1, 0, 0)
  //voisin2 = (nx_ - 1, 1, 0)
  //voisin3 = (nx_ - 1, 0, 1)
  result[0] = (v[0][thisIndex] - v[0][voisinIndex1]) / dx_;
  result[1] = (v[0][voisinIndex2] - v[0][voisinIndex1]) / dx_;
  result[2] = (v[1][voisinIndex3] - v[0][voisinIndex1]) / dx_;
  return 0;
}

int FiniteDiffSolver::GradientVertexPoint1() {
  //voisin1 = (nx_ - 1, 0, 0)
  //voisin2 = (nx_ - 1, 1, 0)
  //voisin3 = (nx_ - 1, 0, 1)
  int a1 = Index(nx_, 0, 0);
  int a2 = Index(nx_ - 1, 0, 0);
  int a3 = Index(nx_ - 1, 1, 0);
  int a4 = Index(nx_ - 1, 0, 1);
  GradientVertexPoint1(gH0, vH0_, a1, a2, a3, a4);  
  GradientVertexPoint1(gE0, vH0_, a1, a2, a3, a4);
  GradientVertexPoint1(g1E, v1H_, a1, a2, a3, a4); 
  GradientVertexPoint1(g2E, v2H_, a1, a2, a3, a4);
  GradientVertexPoint1(g1E, v3H_, a1, a2, a3, a4); 
  GradientVertexPoint1(g2E, v4H_, a1, a2, a3, a4);
  GradientVertexPoint1(g1E, v1E_, a1, a2, a3, a4); 
  GradientVertexPoint1(g2E, v2E_, a1, a2, a3, a4);
  GradientVertexPoint1(g1E, v3E_, a1, a2, a3, a4);  
  GradientVertexPoint1(g2E, v4E_, a1, a2, a3, a4);
  return 0;
}

// Compute the gradient of V2 (0, 1, 0)
int FiniteDiffSolver::GradientVertexPoint2(double* result, double** v, \
                         int thisIndex, int voisinIndex1, int voisinIndex2, int voisinIndex3)
{
  //voisin1 = (0, nx_ - 1, 0)
  //voisin2 = (1, nx_ - 1, 0)
  //voisin3 = (0, nx_ - 1, 1)
  result[0] = (v[0][voisinIndex2] - v[0][voisinIndex1]) /  dx_;
  result[1] = (v[0][thisIndex] - v[0][voisinIndex1]) / dx_;
  result[2] = (v[1][voisinIndex3] - v[0][voisinIndex1]) / dx_;
  return 0;
}

int FiniteDiffSolver::GradientVertexPoint2() {
  //voisin1 = (0, nx_ - 1, 0)
  //voisin2 = (1, nx_ - 1, 0)
  //voisin3 = (0, nx_ - 1, 1)
  int a1 = Index(0, nx_, 0);
  int a2 = Index(0, nx_ - 1, 0);
  int a3 = Index(1, nx_ - 1, 0);
  int a4 = Index(0, nx_ - 1, 1);
  GradientVertexPoint2(gH0, vH0_, a1, a2, a3, a4);  
  GradientVertexPoint2(gE0, vH0_, a1, a2, a3, a4);
  GradientVertexPoint2(g1E, v1H_, a1, a2, a3, a4); 
  GradientVertexPoint2(g2E, v2H_, a1, a2, a3, a4);
  GradientVertexPoint2(g1E, v3H_, a1, a2, a3, a4); 
  GradientVertexPoint2(g2E, v4H_, a1, a2, a3, a4);
  GradientVertexPoint2(g1E, v1E_, a1, a2, a3, a4); 
  GradientVertexPoint2(g2E, v2E_, a1, a2, a3, a4);
  GradientVertexPoint2(g1E, v3E_, a1, a2, a3, a4);  
  GradientVertexPoint2(g2E, v4E_, a1, a2, a3, a4);
  return 0;
}

// Compute the gradient of V3 (0, 0, 1)
int FiniteDiffSolver::GradientVertexPoint3(double* result, double** v, \
                         int thisIndex, int voisinIndex1, int voisinIndex2, int voisinIndex3)
{
  //voisin1 = (0, 0, nx_ - 1)
  //voisin2 = (1, 0, nx_ - 1)
  //voisin3 = (0, 1, nx_ - 1)
  result[0] = (v[nx_ - 1][voisinIndex2] - v[nx_ - 1][voisinIndex1]) /  dx_;
  result[1] = (v[nx_ - 1][voisinIndex3] - v[nx_ - 1][voisinIndex1]) /  dx_;
  result[2] = (v[nx_][thisIndex] - v[nx_ - 1][voisinIndex1]) / dx_;
  return 0;
}

int FiniteDiffSolver::GradientVertexPoint3() {
  //voisin1 = (0, 0, nx_ - 1)
  //voisin2 = (1, 0, nx_ - 1)
  //voisin3 = (0, 1, nx_ - 1)
  int a1 = Index(0, 0, nx_);
  int a2 = Index(0, 0, nx_ - 1);
  int a3 = Index(1, 0, nx_ - 1);
  int a4 = Index(0, 1, nx_ - 1);
  GradientVertexPoint3(gH0, vH0_, a1, a2, a3, a4);  
  GradientVertexPoint3(gE0, vH0_, a1, a2, a3, a4);
  GradientVertexPoint3(g1E, v1H_, a1, a2, a3, a4); 
  GradientVertexPoint3(g2E, v2H_, a1, a2, a3, a4);
  GradientVertexPoint3(g1E, v3H_, a1, a2, a3, a4); 
  GradientVertexPoint3(g2E, v4H_, a1, a2, a3, a4);
  GradientVertexPoint3(g1E, v1E_, a1, a2, a3, a4); 
  GradientVertexPoint3(g2E, v2E_, a1, a2, a3, a4);
  GradientVertexPoint3(g1E, v3E_, a1, a2, a3, a4);  
  GradientVertexPoint3(g2E, v4E_, a1, a2, a3, a4);
  return 0;
}

int FiniteDiffSolver::ComputeOptimalStrategy(int k, int thisIndex, double x1, double x2, double x3) {
  alphaH[k][thisIndex] = model_.AlphaH(vH0_[k][thisIndex], vE0_[k][thisIndex], \
                                       gH0[0], gH0[1], gH0[2], x1, x2, x3);
  alphaE[k][thisIndex] = model_.AlphaE(vH0_[k][thisIndex], vE0_[k][thisIndex], \
                                       gE0[0], gE0[1], gE0[2], x1, x2, x3);
  beta12H[k][thisIndex] = model_.Beta(v1H_[k][thisIndex], v2H_[k][thisIndex]);
  beta34H[k][thisIndex] = model_.Beta(v3H_[k][thisIndex], v4H_[k][thisIndex]);
  beta12E[k][thisIndex] = model_.Beta(v1E_[k][thisIndex], v2E_[k][thisIndex]);
  beta34E[k][thisIndex] = model_.Beta(v3E_[k][thisIndex], v4E_[k][thisIndex]);
  return 0;
}

int FiniteDiffSolver::FiniteDifferenceScheme(int k, int thisIndex, double x1, double x2, double x12) {
  vH0[k][thisIndex] = vH0_[k][thisIndex] + dt_ * \
                      (model_.e0_ * x12 - model_.kA_ * pow(alphaH[k][thisIndex], 2) / 2 + \
                       (model_.cH_ * x12 + model_.dH_ * alphaH[k][thisIndex]) * (vE0_[k][thisIndex] - vH0_[k][thisIndex]) + \
                       gH0[0] * piH[0] + gH0[1] * piH[1] + gH0[2] * piH[2]);
                       
  vE0[k][thisIndex] = vE0_[k][thisIndex] + dt_ * \
                      (model_.e0_ * x12 - model_.kA_ * pow(alphaE[k][thisIndex], 2) / 2 - model_.kE_ + \
                       (model_.cE_ * (1 - x12) + model_.dE_ * (1 - alphaH[k][thisIndex])) * (vH0_[k][thisIndex] - vE0_[k][thisIndex]) + \
                       gE0[0] * piE[0] + gE0[1] * piE[1] + gE0[2] * piE[2]);
                       
  v1H[k][thisIndex] = v1H_[k][thisIndex] + dt_ * \
                      (model_.e_ * (1 - x12) - model_.kI_ - model_.kD_ * pow(beta12H[k][thisIndex], 2) / 2 + \
                       model_.lambda_ * (1 - beta12H[k][thisIndex]) * (v2H_[k][thisIndex] - v1H_[k][thisIndex]) + \
                       model_.rD_ * (v3H_[k][thisIndex] - v1H_[k][thisIndex]) + \
                       (model_.cH_ * x12 + model_.dH_ * alphaH[k][thisIndex]) * (v1E_[k][thisIndex] - v1H_[k][thisIndex]) + \
                       g1H[0] * piH[0] + g1H[1] * piH[1] + g1H[2] * piH[2]);
  
  v2H[k][thisIndex] = v2H_[k][thisIndex] + dt_ * \
                      (model_.e_ * (1 - x12) - model_.kI_ - model_.kD_ * pow(beta12H[k][thisIndex], 2) / 2 + \
                       model_.lambda_ * beta12H[k][thisIndex] * (v1H_[k][thisIndex] - v2H_[k][thisIndex]) + \
                       model_.rU_ * (v4H_[k][thisIndex] - v2H_[k][thisIndex]) + \
                       (model_.cH_ * x12 + model_.dH_ * alphaH[k][thisIndex]) * (v2E_[k][thisIndex] - v2H_[k][thisIndex]) + \
                       g2H[0] * piH[0] + g2H[1] * piH[1] + g2H[2] * piH[2]);
                       
  v3H[k][thisIndex] = v3H_[k][thisIndex] + dt_ * \
                      (model_.e_ * (1 - x12) - model_.kD_ * pow(beta34H[k][thisIndex], 2) / 2 + \
                       model_.lambda_ * (1 - beta34H[k][thisIndex]) * (v4H_[k][thisIndex] - v3H_[k][thisIndex]) + \
                       (model_.bDD_ * x1 + model_.bUD_ * x2 + model_.aD_ * alphaH[k][thisIndex]) * (v1H_[k][thisIndex] - v3H_[k][thisIndex]) +\
                       (model_.cH_ * x12 + model_.dH_ * alphaH[k][thisIndex]) * (v3E_[k][thisIndex] - v3H_[k][thisIndex]) + \
                       g3H[0] * piH[0] + g3H[1] * piH[1] + g3H[2] * piH[2]);
                       
  v4H[k][thisIndex] = v4H_[k][thisIndex] + dt_ * \
                      (model_.e_ * (1 - x12) - model_.kD_ * pow(beta34H[k][thisIndex], 2) / 2 + \
                       model_.lambda_ * beta34H[k][thisIndex] * (v3H_[k][thisIndex] - v4H_[k][thisIndex]) + \
                       (model_.bDU_ * x1 + model_.bUU_ * x2 + model_.aU_ * alphaH[k][thisIndex]) * (v2H_[k][thisIndex] - v4H_[k][thisIndex]) +\
                       (model_.cH_ * x12 + model_.dH_ * alphaH[k][thisIndex]) * (v4E_[k][thisIndex] - v4H_[k][thisIndex]) + \
                       g4H[0] * piH[0] + g4H[1] * piH[1] + g4H[2] * piH[2]);
                       
  v1E[k][thisIndex] = v1E_[k][thisIndex] + dt_ * \
                      (model_.e_ * (1 - x12) - model_.kI_ - model_.kD_ * pow(beta12E[k][thisIndex], 2) / 2 + \
                       model_.lambda_ * (1 - beta12E[k][thisIndex]) * (v2E_[k][thisIndex] - v1E_[k][thisIndex]) + \
                       model_.rD_ * (v3E_[k][thisIndex] - v1E_[k][thisIndex]) + \
                       (model_.cE_ * (1 - x12) + model_.dE_ * (1 - alphaE[k][thisIndex])) * (v1H_[k][thisIndex] - v1E_[k][thisIndex]) + \
                       g1E[0] * piE[0] + g1E[1] * piE[1] + g1E[2] * piE[2]);
  
  v2E[k][thisIndex] = v2E_[k][thisIndex] + dt_ * \
                      (model_.e_ * (1 - x12) - model_.kI_ - model_.kD_ * pow(beta12E[k][thisIndex], 2) / 2 + \
                       model_.lambda_ * beta12E[k][thisIndex] * (v1E_[k][thisIndex] - v2E_[k][thisIndex]) + \
                       model_.rU_ * (v4E_[k][thisIndex] - v2E_[k][thisIndex]) + \
                       (model_.cE_ * (1 - x12) + model_.dE_ * (1 - alphaE[k][thisIndex])) * (v2H_[k][thisIndex] - v2E_[k][thisIndex]) + \
                       g2E[0] * piE[0] + g2E[1] * piE[1] + g2E[2] * piE[2]);
                       
  v3E[k][thisIndex] = v3E_[k][thisIndex] + dt_ * \
                      (model_.e_ * (1 - x12) - model_.kD_ * pow(beta34E[k][thisIndex], 2) / 2 + \
                       model_.lambda_ * (1 - beta34E[k][thisIndex]) * (v4E_[k][thisIndex] - v3E_[k][thisIndex]) + \
                       (model_.bDD_ * x1 + model_.bUD_ * x2 + model_.aD_ * alphaE[k][thisIndex]) * (v1E_[k][thisIndex] - v3E_[k][thisIndex]) +\
                       (model_.cE_ * (1 - x12) + model_.dE_ * (1 - alphaE[k][thisIndex])) * (v3H_[k][thisIndex] - v3E_[k][thisIndex]) + \
                       g3E[0] * piE[0] + g3E[1] * piE[1] + g3E[2] * piE[2]);
                       
  v4E[k][thisIndex] = v4E_[k][thisIndex] + dt_ * \
                      (model_.e_ * (1 - x12) - model_.kD_ * pow(beta34E[k][thisIndex], 2) / 2 + \
                       model_.lambda_ * beta34E[k][thisIndex] * (v3E_[k][thisIndex] - v4E_[k][thisIndex]) + \
                       (model_.bDU_ * x1 + model_.bUU_ * x2 + model_.aU_ * alphaE[k][thisIndex]) * (v2E_[k][thisIndex] - v4E_[k][thisIndex]) +\
                       (model_.cE_ * (1 - x12) + model_.dE_ * (1 - alphaE[k][thisIndex])) * (v4H_[k][thisIndex] - v4E_[k][thisIndex]) + \
                       g4E[0] * piE[0] + g4E[1] * piE[1] + g4E[2] * piE[2]);
  return 0;
}

//Update the inner points i + j + k < nx
int FiniteDiffSolver::UpdateInnerPoint(int i, int j, int k) {
  thisIndex = Index(i, j, k);
  x1 = i * dx_; x2 = j * dx_; x3 = k * dx_;
  x12 = x1 + x2;
  GradientInnerPoint(i, j, k);
  ComputeOptimalStrategy(k, thisIndex, x1, x2, x3);
  
  //Compute operator piH and piE
  model_.Pi(piH, x1, x2, x3, beta12H[k][thisIndex], beta34H[k][thisIndex], alphaH[k][thisIndex]);
  model_.Pi(piE, x1, x2, x3, beta12E[k][thisIndex], beta34E[k][thisIndex], alphaE[k][thisIndex]);
  
  FiniteDifferenceScheme(k, thisIndex, x1, x2, x12);
  return 0;
}
  
//Update the surface points i + j + k = nx, i,j,k > 0
int FiniteDiffSolver::UpdateSurfacePoint(int i, int j, int k) {
  thisIndex = Index(i, j, k);
  x1 = i * dx_; x2 = j * dx_; x3 = k * dx_;
  x12 = x1 + x2;
  GradientSurfacePoint(i, j, k);
  ComputeOptimalStrategy(k, thisIndex, x1, x2, x3);
  
  //Compute operator piH and piE
  model_.Pi(piH, x1, x2, x3, beta12H[k][thisIndex], beta34H[k][thisIndex], alphaH[k][thisIndex]);
  model_.Pi(piE, x1, x2, x3, beta12E[k][thisIndex], beta34E[k][thisIndex], alphaE[k][thisIndex]);
  
  FiniteDifferenceScheme(k, thisIndex, x1, x2, x12);
  return 0;
}
  
//Update the points on edge 
//i + j + k = nx, i = 0, j > 0, k > 0
int FiniteDiffSolver::UpdateEdgePoint1(int i, int j, int k) {
  thisIndex = Index(i, j, k);
  x1 = i * dx_; x2 = j * dx_; x3 = k * dx_;
  x12 = x1 + x2;
  GradientEdgePoint1(i, j, k);
  ComputeOptimalStrategy(k, thisIndex, x1, x2, x3);
  
  //Compute operator piH and piE
  model_.Pi(piH, x1, x2, x3, beta12H[k][thisIndex], beta34H[k][thisIndex], alphaH[k][thisIndex]);
  model_.Pi(piE, x1, x2, x3, beta12E[k][thisIndex], beta34E[k][thisIndex], alphaE[k][thisIndex]);
  
  FiniteDifferenceScheme(k, thisIndex, x1, x2, x12);
  return 0;
}
  
//i + j + k = nx, j = 0, i > 0, k > 0
int FiniteDiffSolver::UpdateEdgePoint2(int i, int j, int k){
  thisIndex = Index(i, j, k);
  x1 = i * dx_; x2 = j * dx_; x3 = k * dx_;
  x12 = x1 + x2;
  GradientEdgePoint2(i, j, k);
  ComputeOptimalStrategy(k, thisIndex, x1, x2, x3);
  
  //Compute operator piH and piE
  model_.Pi(piH, x1, x2, x3, beta12H[k][thisIndex], beta34H[k][thisIndex], alphaH[k][thisIndex]);
  model_.Pi(piE, x1, x2, x3, beta12E[k][thisIndex], beta34E[k][thisIndex], alphaE[k][thisIndex]);

  FiniteDifferenceScheme(k, thisIndex, x1, x2, x12);
  return 0;
}
  
//i + j + k = nx, k = 0, j > 0, i > 0
int FiniteDiffSolver::UpdateEdgePoint3(int i, int j, int k){
  thisIndex = Index(i, j, k);
  x1 = i * dx_; x2 = j * dx_; x3 = k * dx_;
  x12 = x1 + x2;
  GradientEdgePoint3(i, j, k);
  ComputeOptimalStrategy(k, thisIndex, x1, x2, x3);
  
  //Compute operator piH and piE
  model_.Pi(piH, x1, x2, x3, beta12H[k][thisIndex], beta34H[k][thisIndex], alphaH[k][thisIndex]);
  model_.Pi(piE, x1, x2, x3, beta12E[k][thisIndex], beta34E[k][thisIndex], alphaE[k][thisIndex]);
  
  FiniteDifferenceScheme(k, thisIndex, x1, x2, x12);
  return 0;
}
  
//Update the vertex points
//Update (1,0,0)
int FiniteDiffSolver::UpdateVertexPoint1(){
  thisIndex = Index(nx_, 0, 0);
  x1 = 1; x2 = 0; x3 = 0;
  x12 = 1;
  GradientVertexPoint1();
  ComputeOptimalStrategy(0, thisIndex, x1, x2, x3);
  
  //Compute operator piH and piE
  model_.Pi(piH, 1, 0, 0, beta12H[0][thisIndex], beta34H[0][thisIndex], alphaH[0][thisIndex]);
  model_.Pi(piE, 1, 0, 0, beta12E[0][thisIndex], beta34E[0][thisIndex], alphaE[0][thisIndex]);
  
  FiniteDifferenceScheme(0, thisIndex, x1, x2, x12);
  return 0;
}
  
//Update (0,1,0)
int FiniteDiffSolver::UpdateVertexPoint2(){
  thisIndex = Index(0, nx_, 0);
  x1 = 0; x2 = 1; x3 = 0;
  x12 = 1;
  GradientVertexPoint2();
  ComputeOptimalStrategy(0, thisIndex, x1, x2, x3);
  
  //Compute operator piH and piE
  model_.Pi(piH, 0, 1, 0, beta12H[0][thisIndex], beta34H[0][thisIndex], alphaH[0][thisIndex]);
  model_.Pi(piE, 0, 1, 0, beta12E[0][thisIndex], beta34E[0][thisIndex], alphaE[0][thisIndex]);

  FiniteDifferenceScheme(0, thisIndex, x1, x2, x12);
  return 0;
}
  
//Update (0,0,1)
int FiniteDiffSolver::UpdateVertexPoint3(){
  thisIndex = Index(0, 0, nx_);
  x1 = 0; x2 = 0; x3 = 1;
  x12 = 0;
  
  GradientVertexPoint3();
  ComputeOptimalStrategy(nx_, thisIndex, x1, x2, x3);
  
  //Compute operator piH and piE
  model_.Pi(piH, 0, 0, 1, beta12H[nx_][thisIndex], beta34H[nx_][thisIndex], alphaH[nx_][thisIndex]);
  model_.Pi(piE, 0, 0, 1, beta12E[nx_][thisIndex], beta34E[nx_][thisIndex], alphaE[nx_][thisIndex]);
  
  FiniteDifferenceScheme(nx_, thisIndex, x1, x2, x12);
  return 0;
}

// Update the grid points (i, j, k) where 1 <= k <= (nx - 1)
int FiniteDiffSolver::UpdateMidLayer(int k) {
  for (int j = 0; j < nx_ - k; j++) {
    for (int i = 0; i < nx_ - k - j; i++) {
      UpdateInnerPoint(i, j, k);
    }
  }
  // Deal with points on edge (i = nx_ - k, j = 0) and (i = 0, j = nx_ - k) 
  UpdateEdgePoint1(0, (nx_ - k), k);
  UpdateEdgePoint2((nx_ - k), 0, k);
  
  //Deal with points on surface
  if (k < nx_ - 1) {
    for (int i = 1; i < nx_ - k; i++) {
      UpdateSurfacePoint(i, (nx_ - k - i), k);
    }
  }
  return 0;
}

int FiniteDiffSolver::UpdateBottomLayer() {
  for (int j = 0; j < nx_; j++) {
    for (int i = 0; i < nx_ - j; i++) {
      UpdateInnerPoint(i, j, 0);
    }
  }
  UpdateVertexPoint1();
  UpdateVertexPoint2();
  //Deal with points on edge
  for (int i = 1; i < nx_; i++) {
      UpdateEdgePoint3(i, (nx_ - i), 0);
  }
  return 0;
}

int FiniteDiffSolver::UpdateTop() {
  UpdateVertexPoint3();
  return 0;
}

int FiniteDiffSolver::Step() {
  FlipArrays();
  UpdateBottomLayer();
  for (int k = 1; k < nx_; k++) {
    UpdateMidLayer(k);
  }
  UpdateTop();
  return 0;
}

int FiniteDiffSolver::PrintStrategy() {
  for (int k = 0; k <= nx_; k++) {
    for (int j = 0; j <= nx_ - k; j++) {
      for (int i = 0; i <= nx_ - k - j; i++) {
        fprintf(strFileBeta12H, "%15.8f ", beta12H[k][Index(i,j,k)]);
        fprintf(strFileBeta34H, "%15.8f ", beta34H[k][Index(i,j,k)]);
        fprintf(strFileBeta12E, "%15.8f ", beta12E[k][Index(i,j,k)]);
        fprintf(strFileBeta34E, "%15.8f ", beta34E[k][Index(i,j,k)]);
        fprintf(strFileAlphaH, "%15.8f ", alphaH[k][Index(i,j,k)]);
        fprintf(strFileAlphaE, "%15.8f ", alphaE[k][Index(i,j,k)]);
      }
    }
    fprintf(strFileBeta12H, "\n");  fprintf(strFileBeta34H, "\n");
    fprintf(strFileBeta12E, "\n");  fprintf(strFileBeta34E, "\n");
    fprintf(strFileAlphaH, "\n");  fprintf(strFileAlphaE, "\n");
  }
  return 0;
}