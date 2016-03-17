//  Created by peiqiw on 06/03/2016.
//  Copyright (c) 2016 peiqiw. All rights reserved.
//

#ifndef fd_solver_h
#define fd_solver_h

#include <stdio.h>
#include <stdlib.h>
#include "model.h"

using namespace std;

class FiniteDiffSolver {
    
public:

  // parameters of the solver
  const int nx_; //discretization of space
  const int nt_; //discretization of time
  const double dx_;
  const double dt_;
  
  // model
  CyberGameModel &model_;

  // 2-d array storing the value function at current time t
  double ** vH0;
  double ** vE0;
  double ** v1H;
  double ** v2H;
  double ** v3H;
  double ** v4H;
  double ** v1E;
  double ** v2E;
  double ** v3E;
  double ** v4E;
  
  // 2-d array storing the value function at time t-dt
  double ** vH0_;
  double ** vE0_;
  double ** v1H_;
  double ** v2H_;
  double ** v3H_;
  double ** v4H_;
  double ** v1E_;
  double ** v2E_;
  double ** v3E_;
  double ** v4E_;
  
  // 2-d array storing the optimal strategy at current time t
  double ** alphaH;
  double ** alphaE;
  double ** beta12H;
  double ** beta34H;
  double ** beta12E;
  double ** beta34E;
  
  // 1-d array storing gradient of the value function
  double * gH0, * gE0, * g1H, * g2H, * g3H, * g4H, * g1E, * g2E, * g3E, * g4E;
  
  // 1-d array storing piH and piE operator
  double * piH, * piE;

  // Helper variables
  double x1, x2, x3; //points
  double x12; // sum of x1 and x2;
  int thisIndex; //index of this point
  
  // File for output
  FILE * strFileBeta12H, * strFileBeta34H, * strFileBeta12E, * strFileBeta34E;
  FILE * strFileAlphaH, * strFileAlphaE;
  
  // constructor: initialize the solver
  FiniteDiffSolver(int nx, int nt, CyberGameModel &model);

  // destructor
  ~FiniteDiffSolver();
  
  // Index function
  int Index(int i, int j, int k);
  
  // Flip array
  int FlipOneArray(double** a, double** b);
  
  int FlipArrays();
  
  // Compute the gradient of an inner point (i,j,k) of function v, write result.
  int GradientInnerPoint(double* result, double** v, int k, \
                         int thisIndex, int voisinIndex1, int voisinIndex2, int voisinIndex3);
                         
  int GradientInnerPoint(int i, int j, int k);
  
  // Compute the gradient of a surface point (i,j,k) of function v, write result.
  int GradientSurfacePoint(double* result, double** v, int k, \
                         int thisIndex, int voisinIndex1, int voisinIndex2, int voisinIndex3);
  
  int GradientSurfacePoint(int i, int j, int k);
  
  // Compute the gradient of a edge point (i,j,k) x1+x2+x3=1, x1=0
  int GradientEdgePoint1(double* result, double** v, int k, \
                         int thisIndex, int voisinIndex1, int voisinIndex2, int voisinIndex3);
  
  int GradientEdgePoint1(int i, int j, int k);
  
  // Compute the gradient of a edge point (i,j,k) x1+x2+x3=1, x2=0
  int GradientEdgePoint2(double* result, double** v, int k, \
                         int thisIndex, int voisinIndex1, int voisinIndex2, int voisinIndex3);
  
  int GradientEdgePoint2(int i, int j, int k);
  
  // Compute the gradient of a edge point (i,j,k) x1+x2+x3=1, x3=0
  int GradientEdgePoint3(double* result, double** v, int k, \
                         int thisIndex, int voisinIndex1, int voisinIndex2, int voisinIndex3);
  
  int GradientEdgePoint3(int i, int j, int k);
  
  // Compute the gradient of V1 (1, 0, 0)
  int GradientVertexPoint1(double* result, double** v, \
                           int thisIndex, int voisinIndex1, int voisinIndex2, int voisinIndex3);
  
  int GradientVertexPoint1();
  
  // Compute the gradient of V2 (0, 1, 0)
  int GradientVertexPoint2(double* result, double** v, \
                           int thisIndex, int voisinIndex1, int voisinIndex2, int voisinIndex3);
  
  int GradientVertexPoint2();
  
  // Compute the gradient of V3 (0, 0, 1)
  int GradientVertexPoint3(double* result, double** v, \
                           int thisIndex, int voisinIndex1, int voisinIndex2, int voisinIndex3);
  
  int GradientVertexPoint3();
  
  // Update the whole grid
  int Step();
    
  // Update the grid points (i, j, k) where 1 <= k <= (nx -1)
  int UpdateMidLayer(int k);

  // Update the grid points (i, j, 0)
  int UpdateBottomLayer();

  // Update the point (0, 0, nx)
  int UpdateTop();

  //Update the inner points i + j + k < nx
  int UpdateInnerPoint(int i, int j, int k);
  
  //Update the surface points i + j + k = nx, i,j,k > 0
  int UpdateSurfacePoint(int i, int j, int k);
  
  //Update the points on edge 
  //i + j + k = nx, i = 0, j > 0, k > 0
  int UpdateEdgePoint1(int i, int j, int k);
  
  //i + j + k = nx, j = 0, i > 0, k > 0
  int UpdateEdgePoint2(int i, int j, int k);
  
  //i + j + k = nx, k = 0, j > 0, i > 0
  int UpdateEdgePoint3(int i, int j, int k);
  
  //Update the vertex point
  //Update (1,0,0)
  int UpdateVertexPoint1();
  
  //Update (0,1,0)
  int UpdateVertexPoint2();
  
  //Update (0,0,1)
  int UpdateVertexPoint3();
  
  //Finite difference scheme
  int FiniteDifferenceScheme(int k, int thisIndex, double x1, double x2, double x12);
  
  //Compute optimal strategy
  int ComputeOptimalStrategy(int k, int thisIndex, double x1, double x2, double x3);
  
  //Output optimal strategy
  int PrintStrategy();
};

#endif