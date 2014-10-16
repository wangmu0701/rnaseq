/*----------------------------------------------------------------------------
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     ADOL-C_NLP.cpp
 Revision: $$
 Contents: class myADOLC_NPL for interfacing with Ipopt
 
 Copyright (c) Andrea Walther
   
 This file is part of ADOL-C. This software is provided as open source.
 Any use, reproduction, or distribution of the software constitutes 
 recipient's acceptance of the terms of the accompanying license file.
 
 This code is based on the file  MyNLP.cpp contained in the Ipopt package
 with the authors:  Carl Laird, Andreas Waechter   
----------------------------------------------------------------------------*/

/** C++ Example NLP for interfacing a problem with IPOPT and ADOL-C.
 *  MyADOL-C_NLP implements a C++ example showing how to interface 
 *  with IPOPT and ADOL-C through the TNLP interface. This class 
 *  implements the Example 5.1 from "Sparse and Parially Separable
 *  Test Problems for Unconstrained and Equality Constrained
 *  Optimization" by L. Luksan and J. Vlcek ignoring sparsity.
 *
 *  no exploitation of sparsity !!
 *
 */
#include <cassert>
#include <stdlib.h>
#include <math.h>
#include "rapsodia_NLP.hpp"
#include "RALib1/RA1include.ipp"
#include "RALib2/RA2include.ipp"
#include "HigherOrderTensor.hpp"
#define Nexon   2
#define Ntran   3
#define Nlen    4000    

#define pois_fun(l,y)   (exp(-(l))*pow((l),(y)))



double ll0;
double l0[Ntran][2];
double p0[Ntran][2];

double y[Nexon][Nlen];
unsigned int elen[Nexon];
unsigned int estart[Nexon];
unsigned int estop[Nexon];
unsigned int trans[Nexon][2];

unsigned int* rind = NULL;
unsigned int* cind = NULL;
double* value = NULL;
int options[2] = {0, 1};
int nnz;

using namespace Ipopt;

void readData(){
    FILE *fp=fopen("./../../y.dat","r");
    char buf[20000];
    unsigned int i,j,k;
    for(i=0;i<Nexon;i++){
        fgets(buf,20000,fp);
        k=0;
        for(j=0;j<Nlen;j++){
            y[i][j]=atoi(&buf[k]);
            while((buf[k]!=' ') && (buf[k]!='\0')){
                k++;
            }
            if (buf[k]!='\0'){
                k++;
            }
        }
    }
    trans[0][0]=0;trans[0][1]=2;
    trans[1][0]=1;trans[1][1]=2;
    elen[0]=2000;estart[0]=0;estop[0]=2000;
    elen[1]=2000;estart[1]=2000;estop[1]=4000;
    
    //Set a good initial value for start
    fclose(fp);
}

/* Constructor. */
MyADOLC_NLP::MyADOLC_NLP()
{}

MyADOLC_NLP::~MyADOLC_NLP(){}

bool MyADOLC_NLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                         Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  n = 7;

  m = 0;

  // in this example the jacobian is dense. Hence, it contains n*m nonzeros
  nnz_jac_g = 0;

  // the hessian is also dense and has n*n total nonzeros, but we
  // only need the lower left corner (since it is symmetric)
  nnz_h_lag = n*(n-1)/2+n;

  // use the C style indexing (0-based)
  index_style = C_STYLE;

  return true;
}

bool MyADOLC_NLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{
  // Set the bounds for the variables
  for (Index i=0; i<4; i++) {
    x_l[i] = 0;
    x_u[i] = 1e20;
  }
  for (Index i=4; i<n; i++) {
    x_l[i] = 0.0;
    x_u[i] = 1.0;
  }
  return true;
}

bool MyADOLC_NLP::get_starting_point(Index n, bool init_x, Number* x,
                               bool init_z, Number* z_L, Number* z_U,
                               Index m, bool init_lambda,
                               Number* lambda)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the others if
  // you wish.
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  // set the starting point
    x[0]=0.05; //ll0
    x[1]=7.0;  //l0[0][0]
    x[2]=8.0;  //l0[0][1]
    x[3]=9.0;  //l0[0][2]
    x[4]=0.4;  //p0[0][0]
    x[5]=0.3;  //p0[1][0]
    x[6]=0.2;  //p0[2][0]

    readData();
  return true;
}

template<class T> bool  MyADOLC_NLP::eval_obj(Index n, const T *x, T& obj_value)
{
    T ll;
    T l[Ntran][2];
    T p[Ntran][2];
    T* eu[Nexon][Nexon][Nexon];
    T sumtmp[Nlen];
    ll=x[0];
    Index i,j,k,e;
    for(i=0;i<Nexon;i++) {
     for(j=0;j<Nexon;j++) {
      for(k=0;k<Nexon;k++) {
        eu[i][j][k] = (T*)malloc(sizeof(T)*Nlen);
      }
     }
    }
    for(i=0;i<3;i++){
        l[i][0]=x[1+i];
    }
    for(i=0;i<3;i++){
        p[i][0]=x[4+i];
        p[i][1]=1.0-p[i][0];
    }
    Index t0;
    Index t1;
    T tmp=0;
    obj_value=0;
    for(i=0;i<Ntran;i++){l[i][1]=ll;}
    for(e=0;e<Nexon;e++){
        t0=trans[e][0];
        t1=trans[e][1];
        //Initial to 0
        for(i=0;i<elen[e];i++){
            sumtmp[i]=0;
        }
        //Compute sum temporary
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                for(k=0;k<elen[e];k++){
                    sumtmp[k]+=pois_fun(l[t0][i]+l[t1][j],y[e][estart[e]+k])*p[t0][i]*p[t1][j];
                    //                printf("sum[%d]=%15.5f\n",k,sumtmp[k].getValue());
                }
            }
        }
        //    exit(-1);
        //Store averages
        tmp=0;
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                for(k=0;k<elen[e];k++){
                    tmp=pois_fun(l[t0][i]+l[t1][j],y[e][estart[e]+k])*p[t0][i]*p[t1][j];
                    eu[e][i][j][k]=tmp/sumtmp[k];
                }
            }
        }
        
    }
    for(e=0;e<Nexon;e++){
        t0=trans[e][0];
        t1=trans[e][1];
        tmp=0;
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                for(k=0;k<elen[e];k++){
                    tmp+=eu[e][i][j][k]*(log(p[t0][i]*p[t1][j])-l[t0][i]-l[t1][j]+y[e][estart[e]+k]*log(l[t0][i]+l[t1][j]));
                }
            }
        }
        obj_value-=tmp;
    }
    for(i=0;i<Nexon;i++) {
     for(j=0;j<Nexon;j++) {
      for(k=0;k<Nexon;k++) {
        free(eu[i][j][k]);
      }
     }
    }
  return true;
}

template<class T> bool  MyADOLC_NLP::eval_constraints(Index n, const T *x, Index m, T* g)
{
  return true;
}

//*************************************************************************
//
//
//         Nothing has to be changed below this point !!
//
//
//*************************************************************************


bool MyADOLC_NLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  eval_obj(n,x,obj_value);

  return true;
}

bool MyADOLC_NLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{

  int o=1;
  RA1floatD *xad = new RA1floatD[n];
  for(Index i=0; i<n; i++) {
    xad[i] = x[i];
  }
  HigherOrderTensor T(n,o);
  Matrix<unsigned int> SeedMatrix = T.getSeedMatrix();
  int dirs = T.getDirectionCount();
  for(Index i=0; i < dirs; i++) {
    for(Index j=0; j < n; j++) {
      xad[j].set(i+1, 1, SeedMatrix[j][i]);
    }
  }
  RA1floatD y;
  eval_obj(n, xad, y);
  Matrix<double> TaylorCoefficients(o, dirs);
  for(Index i=1; i<=o; i++) {
    for(Index j=1; j<=dirs; j++) {
      TaylorCoefficients[i-1][j-1]=y.get(j,i);
    }
  }
  T.setTaylorCoefficients(TaylorCoefficients);
  std::vector<double> compressedTensor = T.getCompressedTensor(o);
  for(Index i=dirs; i>0; i--) {
    grad_f[dirs-i] = compressedTensor[i-1];
  }
/*
  for(Index i=0; i<n; i++) {
    printf("J[%d]=%.10f\n", i, grad_f[i]);
  }
*/
  delete[] xad;
  return true;
}

bool MyADOLC_NLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{

  return true;
}

bool MyADOLC_NLP::eval_jac_g(Index n, const Number* x, bool new_x,
                       Index m, Index nele_jac, Index* iRow, Index *jCol,
                       Number* values)
{
  return true;
}

bool MyADOLC_NLP::eval_h(Index n, const Number* x, bool new_x,
                   Number obj_factor, Index m, const Number* lambda,
                   bool new_lambda, Index nele_hess, Index* iRow,
                   Index* jCol, Number* values)
{
  if (values == NULL) {
    // return the structure. This is a symmetric matrix, fill the lower left
    // triangle only.

    // the hessian for this problem is actually dense
    Index idx=0;
    for (Index row = 0; row < n; row++) {
      for (Index col = row; col <n; col++) {
        iRow[idx] = col;
        jCol[idx] = row;
        idx++;
      }
    }

    assert(idx == nele_hess);
  }
  else {
    // return the values. This is a symmetric matrix, fill the lower left
    // triangle only

  int o=2;
  RA2floatD *xad = new RA2floatD[n];
  for(Index i=0; i<n; i++) {
    xad[i] = x[i];
  }
  HigherOrderTensor T(n,o);
  Matrix<unsigned int> SeedMatrix = T.getSeedMatrix();
  int dirs = T.getDirectionCount();
  for(Index i=0; i < dirs; i++) {
    for(Index j=0; j < n; j++) {
      xad[j].set(i+1, 1, SeedMatrix[j][i]);
    }
  }
  RA2floatD y;
  eval_obj<RA2floatD>(n, xad, y);
  Matrix<double> TaylorCoefficients(o, dirs);
  for(Index i=1; i<=o; i++) {
    for(Index j=1; j<=dirs; j++) {
      TaylorCoefficients[i-1][j-1]=y.get(j,i);
    }
  }
  T.setTaylorCoefficients(TaylorCoefficients);
  // harvest the compressedTensor

  std::vector<double> compressedTensor = T.getCompressedTensor(o);
  for(Index i=dirs; i>0; i--) {
    values[dirs-i] = compressedTensor[i-1]*obj_factor;
  }
/*
  Index l=0;
  for(Index i=0; i<n; i++) {
    for(Index j=i; j<n; j++) {
      printf("H[%d][%d]=%.10f\n", j, i, values[l]);
      l++;
    }
  }
*/
  delete[] xad;


  } // else
  return true;
}

void MyADOLC_NLP::finalize_solution(SolverReturn status,
                              Index n, const Number* x, const Number* z_L, const Number* z_U,
                              Index m, const Number* g, const Number* lambda,
                              Number obj_value,
			      const IpoptData* ip_data,
			      IpoptCalculatedQuantities* ip_cq)
{
  printf("\n\nSolution of the primal variables, x\n");
  for(Index i=0;i<n;i++){
    printf("x[%d]=%e\n",i,x[i]);
  }
  printf("\n\nObjective value\n");
  printf("f(x*) = %e\n", obj_value);

}


