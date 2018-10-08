#include<stdlib.h>
#include<stdio.h>
#include<iostream>
#include<cmath>
#include<vector>
#include<gsl/gsl_sf_bessel.h>

using namespace std;
using std::max;
using std::min;

#define Double double
#define Int int
#define PI 3.1415926

// This function is used to calculate the integral of the following function
// F(theta,phi) = integral[(sin(theta)^i)*(cos(theta)^j)*(sin(phi)^k)*(cos(phi)^l)*exp(-u*cos(theta))]
// There are 5 input variable:
// Int: i,j,k,l
// Double : u
// The output is Double

Double integralF(Int i, Int j, Int k, Int l, Double u);
Double intThetaOdd(Int m, Double mu);
Double intThetaEven(Int m, Double mu);
Double approchBessel0(Double mu);
Double approchBessel1(Double mu);
Double approchBessel2(Double mu);
Double expanFun(Int n, Double mu);
Double erf(Double mu);
Double intPhi(Int n);
Int factorial(Int n);

int main() {
  Int i,j,k,l;
  Double mu;
  Double result;

  //Intial value
  i = 6;
  j = 8;
  k = 4;
  l = 2;
  mu = 200;

  cout << "The value of input is following : " << endl;
  printf("The integral input is: (i,j,k,l) = (%2d,%2d,%2d,%2d), The double input mu = %3.4lf\n",i,j,k,l,mu);

  result = integralF(i,j,k,l,mu);

  cout << "The integral of the function is " << result << endl;

  return 0;

}

Double integralF(Int i, Int j, Int k, Int l, Double mu){
  // Test if i,j,k,l >= 0. If not, return error message.
  if(i < 0 || j < 0 || k < 0 || l < 0){
    cout << "There is at least one negative value in i,j,k,l. Please check again!" << endl;
    return 0;
  }
  else{
    if(k % 2 == 1 || l % 2 == 1){
      return 0;
    }else{
      Double H = 0.0;
      if(i % 2 == 1){
	// This part calculate the H function when i is an odd number
	Int i1 = (i-1)/2;
	vector<Int> fij; // vector to save the coefficient for H
	
	fij.reserve(i1+1);
	fij.push_back(-1);

	// Construct the coef for fij
	for(Int s=0; s<i1; s++){
	  Double tempfij;
	  tempfij = (i1-s)*fij[s]/(s+1);
	  tempfij = (-1)*tempfij;
	  fij.push_back(tempfij);
	}

	// Accumulate the value for H
	for(Int s=0; s<i1+1; s++){
	  //This is a regression fucntion
	  H += fij[s]*intThetaOdd(j+2*s,mu);
	}

	cout << "The value for H is " << H << endl;
	
	
      }else{
	// This part calculate the H function when i is an even number
	Int i1 = i/2;

	vector<Int> fij; // vector to save the coefficient for H

	fij.reserve(i1+1);
	fij.push_back(1);

	// Construct the coef for fij
	for(Int s=0; s<i1; s++){
	  Double tempfij;
	  tempfij = (i1-s)*fij[s]/(s+1);
	  tempfij = (-1)*tempfij;
	  //cout << "s = " << s << " ,factor = " << tempfij << endl;
	  fij.push_back(tempfij);
	}

	// Accumulate the value for H
	for(Int s=0; s<i1+1; s++){
	  //This is a regression fucntion
	  //cout << "j+2s = " << j+2*s << endl;
	  // cout << "s = " << s << " , intThetaEvem = " << intThetaEven(j+2*s,mu)<< endl;;
	  H += fij[s]*intThetaEven(j+2*s,mu);
	}

	printf("The value for H is %20.16e\n", H);
	
      }

      Double G = 0.0;

      Int k1 = k/2;
      Int l1 = l/2;

      vector<Int> gkl; // vector to save the coefficient for G

      gkl.reserve(k1+l1+1);

      //cout << "the size of gkl is " << k1+l1+1 << endl;

      Double temgkl = 1.0/pow(2,k1+l1+1); 
      // Construct the coef for gkl
      for(Int t=0; t<k1+l1+1;t++){
	//cout << "t is "<< t << " and k1+l1+1 is " << k1+l1+1 << endl;
	Double tempgkl = 0;
	//gkl[t] = 0;
	Int uperBond = min(k1,t);
	Int lowerBond = max(0,t-l1);
	//cout << "the uperBond is " << uperBond << endl;
	//cout << "the lowerBond is " << lowerBond << endl;
	for(Int s=lowerBond; s<uperBond+1; ++s){
	  //cout << "the s value is " << s << endl;
	  //printf("the value of (k1,l1,s,k1-s,t-s,l1-t+s) = (%2d,%2d,%2d,%2d,%2d,%2d)\n",k1,l1,s,k1-s,t-s,l1-t+s);
	  tempgkl += pow(-1,s)*factorial(k1)*factorial(l1)/(factorial(s)*factorial(k1-s))/(factorial(t-s)*factorial(l1-t+s));
	  //cout << "the value of tempgkl is " << tempgkl << endl;
	  //cout << "1" << endl;
	  //cout << "the value of pow(-1,s) is " << pow(-1,s) << endl;
	  //cout << "2" << endl;
	  //gkl[t] += factorial(k1)*factorial(l1)/(factorial(s)*factorial(k1-s))/(factorial(t-s)*factorial(l1-t+s));
	  //gkl[t] *= pow(-1,s);
	}
	gkl.push_back(tempgkl);
      }
	
      //cout << "finish coef for G" << endl;

      //for(Int s=0;s<k1+l1+1;++s){
      //cout << "the value of gkl[" << s << "] is " << gkl[s] << endl;
      //}

      // Accumulate the value for G
      for(Int t=0; t<k1+l1+1; t++){
	// This is a regression function
	G += gkl[t]*intPhi(t);
	//cout << "the intPhi[" << t << "] is " << intPhi(t) << endl;
      }

      G *= temgkl;

      //cout << "the value of temgkl is " << temgkl << endl;
	
      cout << "The value for G is " << G << endl;
	
      return H*G;
    }
  }
}

Double intThetaOdd(Int m, Double mu){
  vector<Double> coef;
  
  coef.reserve(m+1);

  Double temCoef;

  temCoef = factorial(m)*1.0/pow(mu,m+1);
  coef.push_back(temCoef);

  for(Int i=1; i<m+1; i++){
    temCoef = coef[i-1]/i*mu;
    coef.push_back(temCoef);    
  }

  Double int1 = 0;
  Double int2 = 0;
  for(Int i=0; i<m+1; i++){
    int1 += coef[i];
    int2 += coef[i]*pow(-1,i);
  }

  int1 *= exp(-2*mu);
  cout << "int1 = " << int1 << endl;
  cout << "int2 = " << int2 << endl;

  return int1 - int2;
}
      
Double intThetaEven(Int m, Double mu){
  if(m > 2){
    return ((m-1)*intThetaEven(m-1,mu) + mu*intThetaEven(m-2,mu) - (m-2)*intThetaEven(m-3,mu)) / mu;
  }else if(m == 2){
    Double ans2;
    if(mu < 1000){
      ans2 = gsl_sf_bessel_In(1,mu)/exp(mu)*PI/mu + gsl_sf_bessel_In(2,mu)/exp(mu)*PI;
    }else{
      ans2 = approchBessel2(mu);
    }
    return ans2;
  }else if(m == 1){
    Double ans1;
    if(mu < 1000){
      ans1 = -gsl_sf_bessel_In(1,mu)/exp(mu)*PI;
    }else{
      ans1 = approchBessel1(mu);
    }
    return ans1;
  }else{
    Double ans0;
    if(mu < 1000){
      ans0 = gsl_sf_bessel_In(0,mu)/exp(mu)*PI;
    }else{
      ans0 = approchBessel0(mu);
    }
    return ans0;
  }    
}

Double approchBessel0(Double mu){
  Int expanNum = 5; // the number of expansion;
  Double ans = 0.0;

  Double factor;
  for(Int n=0; n<expanNum; n++){
    Double tempAns = expanFun(n,mu);
    factor = 1;
    for(Int m=1; m<n+1; ++m){
      factor = factor*(2*m-1);
    }
    factor = factor*1.0/pow(2,n)/factorial(n);
    ans = ans + tempAns*factor;
  }
  ans = ans*2;
  return ans;
}

Double erf(Double mu){
  Double ans;
  ans = sqrt(PI)/2 + 31/200*exp(-2*mu) - 341/8000*exp(-2*2*mu);
  ans *= 2/sqrt(PI)*sqrt(1-exp(-2*mu));
  return ans;
}

Double expanFun(Int n, Double mu){
  Double expan;
  switch(n){
  case 0:
    expan = sqrt(PI/2)*erf(mu)/2/sqrt(mu);
    break;
  case 1:
    expan = (-4*exp(-2*mu)*sqrt(mu) + sqrt(2*PI)*erf(mu)) / (16*pow(mu,1.5));
    break;
  case 2:
    expan = (-4*exp(-2*mu)*sqrt(mu)*(3+4*mu) + 3*sqrt(2*PI)*erf(mu)) / (64*pow(mu,2.5));
    break;
  case 3:
    expan = (-4*exp(-2*mu)*sqrt(mu)*(15+4*mu*(5+4*mu)) + 15*sqrt(2*PI)*erf(mu)) / (256*pow(mu,3.5));
    break;
  case 4:
    expan = (-4*exp(-2*mu)*sqrt(mu)*(105+4*mu*(35+4*mu*(7+4*mu))) + 105*sqrt(2*PI)*erf(mu)) / (1024*pow(mu,4.5));
    break;
  default:
    cout << "Expansion exceed the expansion number." << endl;
    break;
  }
  return expan;
}

Double approchBessel1(Double mu){
  Int expanNum = 5; // the number of expansion;
  Double ans = 0.0;

  Double factor, factor1, factor2;
  for(Int n=0; n<expanNum; n++){
    Double tempAns = expanFun(n,mu);
    if(n == 0){
      factor1 = 0;
    }else{
      factor1 = 1;
      for(Int m=1; m<n; ++m){
	factor1 = factor1*(2*m-1);
      }
      factor1 = 2*factor1*1.0/pow(2,n-1)/factorial(n-1);
    }
    
    factor2 = 1;
    for(Int m=1; m<n+1; ++m){
      factor2 = factor2*(2*m-1);
    }
    factor2 = factor2*1.0/pow(2,n)/factorial(n);
    
    factor = factor1 - factor2;
    ans = ans + tempAns*factor;
  }
  ans = ans*2;
  return ans;
}

Double approchBessel2(Double mu){
  Int expanNum = 5; // the number of expansion;
  Double ans = 0.0;

  Double factor, factor1, factor2, factor3;
  for(Int n=0; n<expanNum; n++){
    Double tempAns = expanFun(n,mu);

    if(n == 0){
      factor1 = 0;
    }else if(n == 1){
      factor1 = 0;
    }else{
      factor1 = 1;
      for(Int m=1; m<n-1; ++m){
	factor1 = factor1*(2*m-1);
      }
      factor1 = 4*factor1*1.0/pow(2,n-2)/factorial(n-2);
    }

    if(n == 0){
      factor2 = 0;
    }else{
      factor2 = 1;
      for(Int m=1; m<n; ++m){
	factor2 = factor2*(2*m-1);
      }
      factor2 = 4*factor2*1.0/pow(2,n-1)/factorial(n-1);
    }

    factor3 = 1;
    for(Int m=1; m<n+1; ++m){
      factor3 = factor3*(2*m-1);
    }
    factor3 = factor3*1.0/pow(2,n)/factorial(n);

    //cout << "factor1 = " << factor1 << endl;
    //cout << "factor2 = " << factor2 << endl;
    //cout << "factor3 = " << factor3 << endl;

    factor = factor1 - factor2 + factor3;
    
    ans = ans + tempAns*factor;
  }
  
  ans = ans*2;

  return ans;
}


Double intPhi(Int n){
  if(n % 2 == 1){
    return 0;
  }else if(n == 0){
    return 4*PI;
  }else{
    return (n-1)*1.0/n*intPhi(n-2);
  }
}

Int factorial(Int n){
  if(n < 0){
    cout << "error, n is less than 0" << endl;
    return 0;
  }
  if(n == 0 || n == 1){
    return 1;
  }else{
    return n*factorial(n-1);
  }
}


