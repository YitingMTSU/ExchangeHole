#include<stdio.h>
#include<iostream>
#include<stdlib.h>
#include<vector>
#include<cmath>

#include"sumOfIntegral.h"
#include"coefCalculator.h"

using namespace std;

int main() {
  Int i,j,k;
  Double alpha, beta;
  Double result = 0;
  Double norms;
  vector<Double> p,r;

  //Intial value
  i = 1;
  j = 0;
  k = 0;
  alpha = 10.5;
  beta = 20.5;
 
  p.reserve(DIM);
  r.reserve(DIM);

  norms = Double(1);

  p.push_back(0.0);
  p.push_back(1.0);
  p.push_back(0.0);

  r.push_back(0.0);
  r.push_back(0.0);
  r.push_back(1.0);
  
  cout << "The value of input is following : " << endl;
  printf("The integral input is: (i,j,k) = (%2d,%2d,%2d), The double inp\
ut alpha = %3.4lf, beta = %3.4lf\n", i, j, k, alpha, beta);
  //printf("The vector of s is (%3.4f, %3.4f, %3.4f)\n", s[0], s[1], s[2]);
  printf("The vector of p is (%3.4f, %3.4f, %3.4f)\n", p[0], p[1], p[2]);
  printf("The vector of r is (%3.4f, %3.4f, %3.4f)\n", r[0], r[1], r[2]);

  

  result = gaussian(i, j, k, alpha, beta, norms, p, r);

  cout << "The result of gaussian is " << result << endl;

  Double* PA;
  Double* PB;
  Int l1,l2,m1,m2,n1,n2;
  i = 0;
  j = 0;
  k = 0;
  l1 = 1;
  l2 = 1;
  m1 = 1;
  m2 = 1;
  n1 = 1;
  n2 = 1;
  PA = (Double*)malloc(DIM);
  PB = (Double*)malloc(DIM);
  PA[0] = 1;
  PA[1] = 0;
  PA[2] = 0;
  PB[0] = 1;
  PB[1] = 0;
  PB[2] = 1;
  cout << "The coef is Fijk = "<<coefFijk(i,j,k,l1,l2,m1,m2,n1,n2,PA,PB)<< endl;
  return result;

}

