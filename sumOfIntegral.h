#include<vector>
#include"integral.h"

#ifndef DIM
#define DIM 3

using namespace std;

Double gaussian(Int i, Int j, Int k, Double alpha, Double beta, Double norms, vector<Double>& p, vector<Double>& r);
Double regressive(Int i, Int j, Int k, Double mu, Double normq, Double norms, vector< vector<Double> >& invA);
void getexpanM(Int i, Int j, Int k, Double normq, vector< vector<Double> >& invA, vector< vector< vector<Double> > >& expanM);
void getInvM(Double theta0, Double phi0, vector< vector<Double> >& invA);
Double expanFun(Int i, Int j, Int k, vector<Int> position, vector< vector<Double> >& invA, Int nVariable, Double normq);

#endif
