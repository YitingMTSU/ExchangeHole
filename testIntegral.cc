#include<stdio.h>
#include<iostream>
#include<stdlib.h>
#include<vector>
#include<cmath>

#include"integral.h"

using namespace std;
using std::min;
using std::max;

#define DIM 3

Double gaussian(Int i, Int j, Int k, Double alpha, Double beta, Double norms, vector<Double>& p, vector<Double>& r);
Double regressive(Int i, Int j, Int k, Double mu, Double normq, Double norms, vector< vector<Double> >& invA);
void getexpanM(Int i, Int j, Int k, Double normq, vector< vector<Double> >& invA, vector< vector< vector<Double> > >& expanM);
void getInvM(Double theta0, Double phi0, vector< vector<Double> >& invA);
Double expanFun(Int i, Int j, Int k, vector<Int> position, vector< vector<Double> >& invA, Int nVariable, Double normq);

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

  return result;

}

// This is function is to calculate the following function:
// f = x^i * y^j * z^k * exp(-(alpha+beta)*rQ);
// where vector(rQ) = vector(s) + vector(q), vector(q) = vector(r) - vector(p)

Double gaussian(Int i, Int j, Int k, Double alpha, Double beta, Double norms, vector<Double>& p, vector<Double>& r){
  Double normq;
  vector<Double> q;

  q.reserve(DIM);

  for(Int t=0; t<DIM; t++){
    q.push_back(r[t]-p[t]);
  }

  normq = Double(0);

  // Get the mod length for vector s ans q
  for(Int t=0; t<DIM; t++){
    normq += q[t]*q[t];
  }

  normq = sqrt(normq);

  // This part is to get trigonometric functions
  // see formula (1.2)
  Double theta0, phi0;
  theta0 = acos(q[3]/normq);
  phi0 = atan2(q[2],q[1]);

  printf("The vector q = (%.5f,%.5f,%.5f)\n",q[1],q[2],q[3]);

  cout << "the value of theta0 is " << theta0 << " in rad and " << theta0*180/PI << " in degree." << endl;
  
  cout << "the value of phi0 is " << phi0 << " in rad and " << phi0*180/PI << " in degree." << endl;

  vector< vector<Double> > invA;

  //This function is to get the inverse matrix for A
  //see formula 1.4 and 1.5

  getInvM(theta0, phi0, invA);

  //test for expansiom
  //invA[0][0] = 0.1;
  //invA[0][1] = 0.2;
  //invA[0][2] = 0.3;
  //invA[1][0] = 0.2;
  //invA[1][1] = 0.3;
  //invA[1][2] = 0.4;
  //invA[2][0] = 0.3;
  //invA[2][1] = 0.4;
  //invA[2][2] = 0.5;

  for(Int t=0; t<DIM; ++t){
    for(Int s=0; s<DIM; ++s){
      cout << invA[t][s] << " ";
    }
    cout << endl;
  }

  Double result = 0;
  Double gamma;

  gamma = alpha + beta;

  Double factor;

  factor = exp(-gamma*(norms*norms + normq*normq) + 2*gamma*norms*normq);

  cout << "The factor is " << factor << endl;

  Double mu;

  mu = 2*gamma*norms*normq;

  cout << "The value of mu " << mu << endl;

  result = factor * regressive(i, j, k, mu, normq, norms, invA);
      
  return result;
  
}

void getInvM(Double theta0, Double phi0, vector< vector<Double> >& invA){
  invA.reserve(DIM);
  vector<Double> temp;
  temp.reserve(DIM);

  // input the matrix A
  // see formula 1.4
  temp.push_back(cos(theta0) * cos(phi0));
  temp.push_back(-sin(phi0));
  temp.push_back(sin(theta0) * cos(phi0));
  invA.push_back(temp);
  temp.clear();

  temp.push_back(cos(theta0) * sin(phi0));
  temp.push_back(cos(phi0));
  temp.push_back(sin(theta0) * sin(phi0));
  invA.push_back(temp);
  temp.clear();

  temp.push_back(-sin(theta0));
  temp.push_back(Double(0));
  temp.push_back(cos(theta0));
  invA.push_back(temp);
  temp.clear();
}

Double regressive(Int i, Int j, Int k, Double mu, Double normq, Double norms, vector< vector<Double> >& invA){
  vector< vector< vector<Double> > > expanM;
  vector< vector<Double> > tempv;
  vector<Double> temp;
  Double result = 0;

  expanM.reserve(i+j+k+1);
  tempv.reserve(i+j+k+1);
  temp.reserve(i+j+k+1);

  //Initial expansion matrix
  for(Int t=0; t<i+j+k+1; t++){
    temp.push_back(Double(0));
  }
  
  for(Int t=0; t<i+j+k+1; t++){
    tempv.push_back(temp);
  }

  for(Int t=0; t<i+j+k+1; t++){
    expanM.push_back(tempv);
  }
  //cout << "111"<< endl;
  
  getexpanM(i, j, k, normq, invA, expanM);

  //cout << "222" << endl;
  
  Int upper = i+j+k+1;
  for(Int t=0; t<upper; t++){
    for(Int s=0; s<upper; ++s){
      //cout << "(t,s) = "<< t<< s<< " :";
      for(Int r=0; r<upper; ++r){
	//cout<< expanM[t][s][r] << " ";
      }
      //cout<< endl;
    }
  }
  

  for(Int t=0; t<upper; t++){
    for(Int s=0; s<upper-t; ++s){
      for(Int r=0; r<upper-s-t; ++r){
	//printf("The factor expanM[%d,%d,%d] = %lf\n",t,s,r,expanM[t][s][r]);
	//printf("The coef is %lf\n",pow(norms,t+s+r+2));
	//printf("The integrate for integralF(%d,%d,%d,%d,%lf) = %lf\n",t+s+1,r,s,t,mu,integralF(t+s+1,r,s,t,mu));
	result += expanM[t][s][r] * pow(norms,t+s+r+2) * integralF(t+s+1,r,s,t,mu);
      }
    }
  }
  

  return result; 
  
}

void getexpanM(Int i, Int j, Int k, Double normq, vector< vector<Double> >& invA, vector< vector< vector<Double> > >& expanM){

  Int upper = i+j+k+1;
  Int nVariable = DIM;
  for(Int t=0; t<upper; t++){
    for(Int s=0; s<upper-t; ++s){
      for(Int r=0; r<upper-s-t; ++r){
	vector<Int> position;
	position.reserve(DIM);
	position.push_back(t);
	position.push_back(s);
	position.push_back(r);
	//cout << t << " " << s << " " << r << endl;
	expanM[t][s][r] = expanFun(i,j,k,position,invA,nVariable,normq);
      }
    }
  }
}

Double expanFun(Int i, Int j, Int k, vector<Int> position, vector< vector<Double> >& invA, Int nVariable, Double normq){
  Int upperBondl, upperBondm, lowerBondl, lowerBondm;
  Double factor = 0;

  if(nVariable == 0){
    Double a, b, c;
    a = invA[0][2]*normq;
    b = invA[1][2]*normq;
    c = invA[2][2]*normq;
    
    return pow(a,i) * pow(b,j) * pow(c,k);
  }
  
  Double a, b, c;
  Int t;
  a = invA[0][DIM - nVariable];
  b = invA[1][DIM - nVariable];
  c = invA[2][DIM - nVariable];
  t = position[DIM - nVariable];
  upperBondl = min(i,t) + 1;
  lowerBondl = max(0,t-(j+k));
  for(Int l=lowerBondl; l<upperBondl; ++l){
    upperBondm = min(t-l,j) + 1;
    lowerBondm = max(0,t-l-k);
    for(Int m=lowerBondm; m<upperBondm; ++m){
      Double tempfactor;
      tempfactor = pow(a,l) * pow(b,m) * pow(c,t-l-m);
      tempfactor *= factorial(i)/factorial(l)/factorial(i-l);
      tempfactor *= factorial(j)/factorial(m)/factorial(j-m);
      tempfactor *= factorial(k)/factorial(t-l-m)/factorial(k-t+l+m);
      factor += tempfactor*expanFun(i-l, j-m, k+l+m-t, position, invA, nVariable-1, normq);
    }
  }
 
  return factor;
}
