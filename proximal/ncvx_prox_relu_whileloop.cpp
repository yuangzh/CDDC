#include <math.h>
#include "mex.h"
#include <time.h>
#include <stdlib.h>

#define  abs1(a)         ((a) < 0.0 ? -(a) : (a))
#define  sign1(a)        ((a)==0) ? 0 : (((a)>0.0)?1:(-1))
#define  max1(a,b)       ((a) > (b) ? (a) : (b))
#define  min1(a,b)       ((a) < (b) ? (a) : (b))


double hf(double t, double alpha, double beta, double* a, double* b, double* w, mwSize n)
{
    // 0.5*alpha*t^2+beta*t+0.5*norm(max(0,t*a+b))^2-w'*max(0,t*a+b);
    int i;
    double result = 0.5 * alpha * t * t + beta * t;
    for (i = 0; i < n; i++)
    {
        double f0 = max1(0.0, t * a[i] + b[i]);
        double f1 = w[i]*f0;
        double f2 = 0.5 * pow(f0,2);
        result = result - f1 + f2 ;
    }
    return result;
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    double alpha = mxGetScalar(prhs[0]);
    double beta = mxGetScalar(prhs[1]);
    double *a = mxGetPr(prhs[2]);
    double *b = mxGetPr(prhs[3]);
    double *w = mxGetPr(prhs[4]);
    double *z = mxGetPr(prhs[5]);
    double *q1 = mxGetPr(prhs[6]);
    double *q2 = mxGetPr(prhs[7]);
    double *p1 = mxGetPr(prhs[8]);
    double *p2 = mxGetPr(prhs[9]);
    double *a1 = mxGetPr(prhs[10]);
    double *a2 = mxGetPr(prhs[11]);
    double *b1 = mxGetPr(prhs[12]);
    double *b2 = mxGetPr(prhs[13]);
    int len_a = mxGetScalar(prhs[14]);
    int i1 = mxGetScalar(prhs[15]);
    int i2 = mxGetScalar(prhs[16]);
    double r1 = mxGetScalar(prhs[17]);
    double r2 = mxGetScalar(prhs[18]);
    double s1 = mxGetScalar(prhs[19]);
    double s2 = mxGetScalar(prhs[20]);
    
    double best_t = -beta/alpha;
    double best_fobj =  hf(best_t, alpha, beta, a, b, w, len_a);
    int i;
    
    if(len_a>0)
    {
        for (i=1; i<=len_a+1; i++)
        {
            double tt = tt=(z[i-1]+z[i])/2;
            
            while ((i1>0)&&(tt*a1[i1-1]+b1[i1-1]>0))
            {
                r1=r1+q1[i1-1]; s1=s1+p1[i1-1]; i1=i1-1;
            }
            while ((i2>0)&&(tt*a2[i2-1]+b2[i2-1]<=0))
            {
                r2=r2-q2[i2-1]; s2=s2-p2[i2-1]; i2=i2-1;
            }
            tt=(s1+s2-beta)/(r1+r2+alpha);
            if((tt>=z[i-1]-1e-16)&&(tt<=z[i]+1e-16))
            {
                double curr_fobj = hf(tt, alpha, beta, a, b, w, len_a);
                if( curr_fobj < best_fobj)
                {
                    best_fobj = curr_fobj;
                    best_t = tt;
                }
            }
        }
    }
    
    plhs[0] = mxCreateDoubleScalar(best_t);
    
}
