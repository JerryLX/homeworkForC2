/***************************************************************************
 *
 *   File        : tasks.c
 *   Student Id  : <INSERT STUDENT ID HERE>
 *   Name        : <INSERT STUDENT NAME HERE>
 *
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include "tasks.h"
#include <unistd.h>

#define PI 3.14159265

#define MAX_ITERATION 50000
#define MAX_BUF_LEN 2048
#define OUT_LINALSYS "out_linalsys.csv"
#define OUT_SHOCKb "out_shock_b.csv"
#define OUT_SHOCKc "out_shock.csv"
#define OUT_INTERP "out_interp.csv"
#define OUT_LINALSYS "out_linalsys.csv"
#define OUT_heateqn_1 "out_heateqn_explicit_fe.csv"
#define OUT_heateqn_2 "out_heateqn_explicit_ve.csv"
#define OUT_heateqn_3 "out_heateqn_implicit_fe.csv"

void* safe_malloc(size_t num_bytes)
{
    void* ptr = malloc(num_bytes);
    if (ptr == NULL) {
            printf("ERROR: malloc(%lu)\n", num_bytes);
            exit(EXIT_FAILURE);
    }
    memset(ptr,0,num_bytes);
    return ptr;
}

FILE* safe_fopen(const char* path, const char* mode)
{
    FILE* fp = fopen(path, mode);
    if (fp == NULL) {
        printf("ERROR: %s not exist!\n",path);
        exit(EXIT_FAILURE);
    }
    return fp;
}

double eq2(double beta, double theta, double M, double gamma){
    double p1 = M*M*pow((sin(beta*PI/180)),2)-1;
    double p2 = M*M*(gamma+cos(2*beta*PI/180))+2;
    return 2/(tan(beta*PI/180))*p1/p2-tan(theta*PI/180);
}

double deq2(double beta, double M, double gamma){
    double p1 = M*M*pow((sin(beta*PI/180)),2)-1;
    double p2 = M*M*(gamma+cos(2*beta*PI/180))+2;
    double c1 = -2/pow(sin(beta*PI/180),2)*p1/p2;
    double c2 = M*M*sin(2*beta*PI/180)*p2;
    double c3 = p1*(-2*M*M*sin(2*beta*PI/180));
    double c4 = pow(p2,2);
    return c1+(c2-c3)/c4;
}

/*
double newton2(double init, double theta, double M, double gamma){
    double beta = init, beta2 = init+10;
    double e0 = eq2(beta,theta,M,gamma);
    double e1 = eq2(beta2,theta,M,gamma);
    int i;
    for(i=0;i<100;i++){
        double newbeta = (init*e1-beta*e0)/(e1-e0);
        printf("%lf\n",e1-e0);
        e1 = eq2(newbeta,theta,M,gamma);
        beta = newbeta;
    }
    return beta;
}
*/

double newton(double init, double theta, double M, double gamma){
    double beta = init;
    double e0 = eq2(beta,theta,M,gamma);
    double de0 = deq2(beta,M,gamma);
    int i;
    for(i=0;i<MAX_ITERATION;i++){
        double newbeta = beta-e0/de0;
        if(beta == newbeta)  //converge
            break;
        e0 = eq2(newbeta,theta,M,gamma);
        beta = newbeta;
    }
    return beta;
}

void shockwave(const char* q2_file)
{
    FILE *fp, *fpw, *fpw2;
    double M, theta, beta_l, beta_u, gamma;
    double beta_L_init = 5;
    double beta_U_init = 80;
    double theta_max = 90;

    char buf[MAX_BUF_LEN];

    fp = safe_fopen(q2_file, "r");
    fgets(buf, MAX_BUF_LEN,fp); //skip first line

    //part (a) M = 5.0 theta = 20
    fscanf(fp,"%lf,%lf,%lf,%lf,%lf\n",&M,&theta,&beta_l,&beta_u,&gamma);
    double beta_L = newton(beta_L_init,theta,M,gamma);
    double beta_U = newton(beta_U_init,theta,M,gamma);
    printf("beta_L = %.6lf, beta_U = %.6lf\n",beta_L, beta_U);

    //part (b) M = 5.0
    fpw = safe_fopen(OUT_SHOCKb,"w");
    fprintf(fpw, "theta,beta_lower,beta_upper\n");
    fprintf(fpw, "%.6lf,%.6lf,%.6lf\n",(double)0,asin(1/M)*180/PI,(double)90);

    for(int theta_i=1;theta_i<theta_max;theta_i++){
        beta_L = newton(beta_L_init,theta_i,M,gamma);
        beta_U = newton(beta_U_init,theta_i,M,gamma);
        if(beta_U < 0) break;
        fprintf(fpw, "%d,%.6lf,%.6lf\n",theta_i,beta_L,beta_U);
    }

    //part (c)    
    fpw2 = safe_fopen(OUT_SHOCKc,"w");
    fprintf(fpw2,"M,theta,beta_lower,beta_upper\n");
    fgets(buf, MAX_BUF_LEN,fp); //skip first line
    while(!feof(fp)){
        fscanf(fp, "%lf\n", &M);
        fprintf(fpw2, "%.6lf,%d,%.6lf,%.6lf\n",M,0,asin(1/M)*180/PI,(double)90);
        for(int theta_i=1;theta_i<theta_max;theta_i++){
            beta_L = newton(beta_L_init,theta_i,M,gamma);
            beta_U = newton(beta_U_init,theta_i,M,gamma);
            if(beta_U < 0) break;
            fprintf(fpw2, "%.6lf,%d,%.6lf,%.6lf\n",M,theta_i,beta_L,beta_U);
        }
    }
    
    fclose(fp);
    fclose(fpw);
    fclose(fpw2);

    getcwd(buf,sizeof(buf));   
    printf("current working directory: %s\n", buf);
}



//function to calculate tri-diagonal system
void TDMA(double *X, const int n, double *A, double *B, double *C, double *D){
    int i;
    double tmp;

    C[0] = C[0] / B[0];
    D[0] = D[0] / B[0];

    for(i=1;i<n;i++){
        tmp = (B[i]-A[i]*C[i-1]);
        C[i] = C[i]/tmp;
        D[i] = (D[i]-A[i]*D[i-1])/tmp;
    }

    X[n-1] = D[n-1];

    for(i=n-2;i>=0;i--){
        X[i] = D[i]-C[i]*X[i+1];
    }
}

void linalgbsys(const char* q4_file)
{
    FILE *fp, *fpw;
    double *ai, *bi, *ci, *Qi, *xi;
    double a,b,c,q;
    int count = 0;
    char buf[MAX_BUF_LEN];

    fp = safe_fopen(q4_file, "r");
    fgets(buf, MAX_BUF_LEN,fp); //skip first line
    while(fscanf(fp,"%lf,%lf,%lf,%lf\n",&a,&b,&c,&q) == 4) {
        count++;
    }
    if(!count){
        printf("no data in q4_file!\n");
        exit(EXIT_FAILURE);
    }
    rewind(fp); //jump back to first line
    fgets(buf, MAX_BUF_LEN,fp); //skip first line
    ai = safe_malloc(sizeof(double)*count);
    bi = safe_malloc(sizeof(double)*count);
    ci = safe_malloc(sizeof(double)*count);
    Qi = safe_malloc(sizeof(double)*count);
    xi = safe_malloc(sizeof(double)*count);

    
    for(int i=0;i<count;i++){
        fscanf(fp,"%lf,%lf,%lf,%lf\n",&ai[i],&bi[i],&ci[i],&Qi[i]);
    }

    TDMA(xi,count,ci,ai,bi,Qi);

//    bi[0] = bi[0]/ai[0];
//    Qi[0] = Qi[0]/ai[0];
//    for(int i=1;i<count;i++){
//        bi[i] = bi[i]/(ai[i]-ci[i]*bi[i-1]);
//        Qi[i] = (Qi[i]-ci[i]*Qi[i-1])/(ai[i]-ci[i]*bi[i-1]);
//    }
//    xi[count-1] = Qi[count-1];
//    for(int i=count-2;i>=(int)0;i--){
//        xi[i] = (Qi[i]-bi[i]*xi[i+1]);
//    }

    fpw = safe_fopen(OUT_LINALSYS,"w");
    fprintf(fpw, "x\n");
    for(int i=0;i<count;i++){
        fprintf(fpw,"%.6lf\n",xi[i]);
    }

    fclose(fp);
    fclose(fpw);

    free(ai);
    free(bi);
    free(ci);
    free(Qi);
    free(xi);
}

double lagrange(double *x, double *y, int count, double x0)
{
    if(count < 3){
        printf("second order lagrange needs at least 3 points!\n");
        exit(EXIT_FAILURE);
    }
//    x++;y++;
    //use first three value to calculate
    double l0 = (x0-x[1])*(x0-x[2])/(x[0]-x[1])/(x[0]-x[2]);
    double l1 = (x0-x[0])*(x0-x[2])/(x[1]-x[0])/(x[1]-x[2]);
    double l2 = (x0-x[0])*(x0-x[1])/(x[2]-x[0])/(x[2]-x[1]);

    return y[0]*l0+y[1]*l1+y[2]*l2;
}

double cubic_spline(double *x, double *y, int count, double x0)
{
    double *ai = safe_malloc(sizeof(double)* (count-1));
    double *bi = safe_malloc(sizeof(double)* (count-1));
    double *ci = safe_malloc(sizeof(double)* (count-1));
    double *di = safe_malloc(sizeof(double)* (count-1));

    double *h = safe_malloc(sizeof(double) * (count-1));
    
    double *A = safe_malloc(sizeof(double) * (count-2));
    double *B = safe_malloc(sizeof(double) * (count-2));
    double *C = safe_malloc(sizeof(double) * (count-2));
    double *D = safe_malloc(sizeof(double) * (count-2));
    double *E = safe_malloc(sizeof(double) * (count-2));

    double *M = safe_malloc(sizeof(double) * (count));


    for(int i=0;i<count-1;i++){
        h[i] = x[i+1]-x[i];
    }
    
    for(int i=0;i<count-3;i++){
        A[i] = h[i];
        B[i] = 2 *(h[i]+h[i+1]);
        C[i] = h[i+1];
    }

    for(int i=0;i<count-3;i++){
        D[i] = 6 *((y[i+2] - y[1]) / h[i+1] - (y[i+1]-y[i])/h[i]);
    }

    TDMA(E, count-3, A,B,C,D);
    M[0] = 0;
    M[count-1] = 0;
    for(int i=1;i<count-1;i++){
        M[i] = E[i-1];
    }

    for(int i=0;i<count-1;i++){
        ai[i] = y[i];
        bi[i] = (y[i+1]-y[i])/h[i] - h[i]*M[i]/2 - h[i]*(M[i+1]-M[i])/6;
        ci[i] = M[i]/2;
        di[i] = (M[i+1]-M[i])/(6*h[i]);
    }

    int k=0; //the index of the range x0 belongs to
    for(;k<count-1;k++){
        if(x0>=x[k] && x0<x[k+1])
            break;
    }
    if(k == count-1){
        printf("x0=%.6lf cannot be computed use cubic spline\n",x0);
        return -1;
    }

    double res = ai[k]+bi[k]*(x0-x[k])+ci[k]*pow((x0-x[k]),2)+di[k]*pow((x0-x[k]),3);
    free(h);
    free(A);
    free(B);
    free(C);
    free(D);
    free(E);
    free(M);
    free(ai);
    free(bi);
    free(ci);
    free(di);
    return res;
}


void interp(const char* q5_file, const double xo)
{
    double *x, *y;
    double xi, yi;
    FILE *fp, *fpw;
    char buf[MAX_BUF_LEN];
    int count = 0;

    fp = safe_fopen(q5_file,"r");
    fgets(buf, MAX_BUF_LEN,fp); //skip first line
    
    while(fscanf(fp,"%lf,%lf\n",&xi,&yi) == 2) {
        count++;
    }
    rewind(fp); //jump back to first line
    fgets(buf, MAX_BUF_LEN,fp); //skip first line
    x = safe_malloc(sizeof(double)*count);
    y = safe_malloc(sizeof(double)*count);
    
    for(int i=0;i<count;i++){
        fscanf(fp,"%lf,%lf\n",&x[i],&y[i]);
    }
 
    double res1 = lagrange(x,y,count,xo);
    double res2 = cubic_spline(x,y,count,xo);
    fpw = safe_fopen(OUT_INTERP,"w");
    fprintf(fpw,"lagrange\n");
    fprintf(fpw,"%.6lf\n",res1);
    fprintf(fpw,"cubic\n");
    fprintf(fpw,"%.6lf\n",res2);

    fclose(fp);
    fclose(fpw);
    
    free(x);
    free(y);
}

void heateqn(const char* q6_file)
{
   FILE *fp, *fpw1, *fpw2, *fpw3;
   char buf[MAX_BUF_LEN];
   double mu = -1;
   int Nx=-1,Nt=-1;
   fp = safe_fopen(q6_file, "r");
   fgets(buf, MAX_BUF_LEN,fp); //skip first line
   fscanf(fp,"%lf,%d,%d\n",&mu,&Nx,&Nt);
   printf("%lf,%d,%d\n",mu,Nx,Nt);

   double deltax = 1/(double)Nx,deltat = 1/(double)Nt; 
   //printf("double: %lf\n",deltax);   
   
   //explicit-fixed-ends
   double** res1;
   res1 = (double**)malloc(sizeof(double*)*(Nt+1));
   for(int i = 0; i < Nt+1; i++)
     res1[i] = (double*)malloc(sizeof(double)*(Nx+1));
   
   for(int i = 0; i < Nx+1; i++){
	double temp_x = deltax*i;
   	if(temp_x>=0 && temp_x<0.125)
	  res1[0][i] = 0;
	else if(temp_x <=0.375)
	  res1[0][i] = 0.5*(1-cos(8*PI*(temp_x-0.125)));
	else if(temp_x <= 1)
 	  res1[0][i] = 0;
	else
	  printf("error occurs.\n");
   }

   for(int i = 1; i < Nt+1; i++)
     for(int j = 0; j < Nx+1; j++){
	if(j ==0 || j == Nx) 
	  res1[i][j] = res1[0][j];
	else {
          double diff = (res1[i-1][j+1]+res1[i-1][j-1]-2*res1[i-1][j])/(deltax*deltax);
	  res1[i][j] = res1[i-1][j]+deltat*mu*diff;
	}

   }
  
   fpw1 = safe_fopen(OUT_heateqn_1,"w");
   fprintf(fpw1, "x,f(x)\n");
   for(int i = 0; i < Nx+1; i++)
     fprintf(fpw1,"%lf,%lf\n",i*deltax,res1[Nt][i]);
   fclose(fpw1);
   for(int i = 0; i < Nt+1; i++)
     free(res1[i]);
   free(res1);   
   //explicit-variable-ends
   double** res2;
   res2 = (double**)malloc(sizeof(double*)*(Nt+1));
   for(int i = 0; i < Nt+1; i++)
     res2[i] = (double*)malloc(sizeof(double)*(Nx+1));
   
   for(int i = 0; i < Nx+1; i++){
	double temp_x = deltax*i;
   	if(temp_x>=0 && temp_x<0.125)
	  res2[0][i] = 0;
	else if(temp_x <=0.375)
	  res2[0][i] = 0.5*(1-cos(8*PI*(temp_x-0.125)));
	else if(temp_x <= 1)
 	  res2[0][i] = 0;
	else
	  printf("error occurs.\n");
   }

   for(int i = 1; i < Nt+1; i++)
     for(int j = 0; j < Nx+1; j++){
	if(j ==0){
	  double diff = (res2[i-1][j]+res2[i-1][j+2]-2*res2[i-1][j+1])/(deltax*deltax);
	  res2[i][j] = res2[i-1][j] + deltat*mu*diff;
	}
        else if(j == Nx) {
	  double diff = (res2[i-1][j]+res2[i-1][j-2]-2*res2[i-1][j-1])/(deltax*deltax);
	  res2[i][j] = res2[i-1][j] + deltat*mu*diff;
	}
	else {
          double diff = (res2[i-1][j+1]+res2[i-1][j-1]-2*res2[i-1][j])/(deltax*deltax);
	  res2[i][j] = res2[i-1][j]+deltat*mu*diff;
	}

   }
  
   fpw2 = safe_fopen(OUT_heateqn_2,"w");
   fprintf(fpw2, "x,f(x)\n");
   for(int i = 0; i < Nx+1; i++)
     fprintf(fpw2,"%lf,%lf\n",i*deltax,res2[Nt][i]);
   fclose(fpw2);
   for(int i = 0; i < Nt+1; i++)
     free(res2[i]);
   free(res2);   

   //implicit-fixed-ends
   double** res3;
   res3 = (double**)malloc(sizeof(double*)*(Nt+1));
   for(int i = 0; i < Nt+1; i++)
     res3[i] = (double*)malloc(sizeof(double)*(Nx+1));
   
   for(int i = 0; i < Nx+1; i++){
	double temp_x = deltax*i;
   	if(temp_x>=0 && temp_x<0.125)
	  res3[0][i] = 0;
	else if(temp_x <=0.375)
	  res3[0][i] = 0.5*(1-cos(8*PI*(temp_x-0.125)));
	else if(temp_x <= 1)
 	  res3[0][i] = 0;
	else
	  printf("error occurs.\n");
   }
   double a = (deltat*mu)/(deltax*deltax);
   double* A = (double*)malloc(sizeof(double)*(Nx+1));
   double* B = (double*)malloc(sizeof(double)*(Nx+1));
   double* C = (double*)malloc(sizeof(double)*(Nx+1));
   A[0] = 0;
   C[Nx] = 0;
   for(int i = 1; i < Nx+1; i++)
	A[i] = a;
   for(int i = 0; i < Nx+1; i++)
	B[i] = (2*a+1);
   for(int i = 0; i < Nx; i++)
	C[i] = a;
   for(int i = 1; i < Nt+1; i++) {
   	TDMA(res3[i],Nx+1,A,B,C,res3[i-1]);  
   }
   free(A);
   free(B);
   free(C);
  
   fpw3 = safe_fopen(OUT_heateqn_3,"w");
   fprintf(fpw3, "x,f(x)\n");
   for(int i = 0; i < Nx+1; i++)
     fprintf(fpw3,"%lf,%lf\n",i*deltax,res3[Nt][i]);
   fclose(fpw3);
   for(int i = 0; i < Nt+1; i++)
     free(res3[i]);
   free(res3);   

   //close file
   fclose(fp);
  
}
