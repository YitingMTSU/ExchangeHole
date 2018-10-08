#include<stdlib.h>
#include<stdio>
#include<iostream>
#include<vector>


using namespace std;
#include"shellprop.h"
#include"sumOfIntegral.h"
#include"coefCalculator.h"

typedef int             Int;
typedef size_t          UInt;
#ifdef WITH_SINGLE_PRECISION
typedef float           Double;
#else
typedef double          Double;
#endif

void sphereInt(const LInt& Lcode, const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd){
  //loop over grid points
  for(UInt iGrid=0; iGrid<nGrids; iGrid++){

    // get how many base functions from A and B
    Int lAmin,lAmax,lBmin,lBmax;
    decodeL(Lcode,lAmin,lAmax,lBmin,lBmax);
    
    Int totalNumBas = getCartBas(lAmin,lAmax)*getCartBas(lBmin,lBmax);

    for(Int index=0; index<totalNumBas; index++){
      getlmn(lAmin,lAmax,index,lA,mA,nA);
      getlmn(lBmin,lBmax,index,lB,mB,nB);

    for(UInt ip2=0; ip2<inp2; ip2++){
      //Intial number 
      Double ic2 = icoe[ip2];
      Double alpha = ;
      Double beta = ;
      Double norms = 1;
      UInt offsetP = 3*ip2;
      UInt offsetR = 3*iGrid;
      Double PX = P[offsetP ];
      Double PY = P[offsetP+1];
      Double PZ = P[offsetP+2];
      Double RX = R[offsetR ];
      Double RY = R[offsetR+1];
      Double RZ = R[offsetR+2];
      Double PAX = PX - A[0];
      Double PAY = PY - A[1];
      Double PAZ = PZ - A[2];
      Double PBX = PX - B[0];
      Double PBY = PY - B[1];
      Double PBZ = PZ - B[2];
      
      vector<Double> p;
      p.reserve(DIM);
      p.push_back(PX);
      p.push_back(PY);
      p.push_back(PZ);

      vector<Double> r;
      r.reserve(DIM);
      r.push_back(RX);
      r.push_back(RY);
      r.push_back(RZ);

      Double* PA;
      Double* PB;
      PA = (Double*)malloc(DIM);
      PB = (Double*)malloc(DIM);

      PA[0] = PAX;
      PA[1] = PAY;
      PA[2] = PAZ;
      PB[0] = PBX;
      PB[1] = PBY;
      PB[2] = PBZ;

      //


    }
  }
}
