#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <fftw3.h>
#include <string>



using namespace std;


double sign(double x)
{
  if (x > 0.0)
    return 1;
  else
    return -1;
}

class myComplex {
  fftw_complex val;
public:
  myComplex()
  {
    val[0] = 0;
    val[1] = 0;
  }
  myComplex(double Re, double Im)
  {
    val[0] = Re;
    val[1] = Im;
  }

  double getRe()
  {
    return val[0];
  }

  double getIm()
  {
    return val[1];
  }

  myComplex operator*(myComplex a)
  {
    myComplex tmp(val[0] * a.getRe() - val[1] * a.getIm(), val[0] * a.getIm() + val[1] * a.getRe());
    return tmp;
  }

  myComplex operator+(myComplex a)
  {
    myComplex tmp(val[0] + a.getRe(), val[1] + a.getIm());
    return tmp;
  }

  myComplex operator-(myComplex a)
  {
    myComplex tmp(val[0] - a.getRe(), val[1] - a.getIm());
    return tmp;
  }

  myComplex operator*(double a)
  {
    myComplex tmp(val[0] * a, val[1] * a);
    return tmp;
  }

};


int main() {

  int n, m, h, i, j, k, t;
  double  dt, dx, dy, dz, c = 299792458.0;
  double x, y, z;

  double L = 0.0001;
  double X_Min = -0.03;
  double X_Max = 0.05;
  double Y_Min = -0.01;
  double Y_Max = 0.01;
  double Z_Min = -0.01;
  double Z_Max = 0.01;

  t = 1000;
  n = 128;
  m = 32;
  h = 32;
  dx = (X_Max - X_Min) / n;
  dy = (Y_Max - Y_Min) / m;
  dz = (Z_Max - Z_Min) / h;
  dt = dy / (4.0);

  //ifstream fin("input.txt");
  //  fin >> t;
  //  fin >> dt;
  //  fin >> dx;
  //  fin >> dy;
  //  fin >> dz; 
  //  fin >> n;
  //  fin >> m; 
  //  fin >> h; 

  double *Ex = new double[n*m*h];
  double *Ey = new double[n*m*h];
  double *Ez = new double[n*m*h];
  double *Bx = new double[n*m*h];
  double *By = new double[n*m*h];
  double *Bz = new double[n*m*h];

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < m; j++)
    {
      for (k = 0; k < h; k++)
      {
        x = X_Min + i*dx;
        y = Y_Min + j*dy;
        z = Z_Min + k*dz;
        Ex[i + j*n + k*n*m] = 0;
        Ey[i + j*n + k*n*m] = ((((exp(((((-((y * y) + (z * z))) * 2) * 0.69314718055994529) / 0.000025000000000000001)))) *
          (((exp((((-((((2 * x) / (-0.030000000000000002))) * (((2 * x) / (-0.030000000000000002))))) * 2) * 0.69314718055994529)) *
          (((0.5 * (1 + sign(((-x) - (-0.030000000000000002))))))))))) * cos((2617.9938779914942 * x)));
        Ez[i + j*n + k*n*m] = 0;

        Bx[i + j*n + k*n*m] = 0;
        By[i + j*n + k*n*m] = 0;
        Bz[i + j*n + k*n*m] = ((((exp(((((-((y * y) + (z * z))) * 2) * 0.69314718055994529) / 0.000025000000000000001)))) *
          (((exp((((-((((2 * x) / (-0.030000000000000002))) * (((2 * x) / (-0.030000000000000002))))) * 2) * 0.69314718055994529)) *
          (((0.5 * (1 + sign(((-x) - (-0.030000000000000002))))))))))) * cos((2617.9938779914942 * x)));
      }
    }
  }


  myComplex *Ex_c = new myComplex[n*m*h];
  myComplex *Ey_c = new myComplex[n*m*h];
  myComplex *Ez_c = new myComplex[n*m*h];
  myComplex *Bx_c = new myComplex[n*m*h];
  myComplex *By_c = new myComplex[n*m*h];
  myComplex *Bz_c = new myComplex[n*m*h];



  fftw_plan pEx = fftw_plan_dft_r2c_3d(n, m, h, (double *)Ex, (fftw_complex *)Ex_c, FFTW_MEASURE);
  fftw_plan pEy = fftw_plan_dft_r2c_3d(n, m, h, (double *)Ey, (fftw_complex *)Ey_c, FFTW_MEASURE);
  fftw_plan pEz = fftw_plan_dft_r2c_3d(n, m, h, (double *)Ez, (fftw_complex *)Ez_c, FFTW_MEASURE);
  fftw_plan pBx = fftw_plan_dft_r2c_3d(n, m, h, (double *)Bx, (fftw_complex *)Bx_c, FFTW_MEASURE);
  fftw_plan pBy = fftw_plan_dft_r2c_3d(n, m, h, (double *)By, (fftw_complex *)By_c, FFTW_MEASURE);
  fftw_plan pBz = fftw_plan_dft_r2c_3d(n, m, h, (double *)Bz, (fftw_complex *)Bz_c, FFTW_MEASURE);

  fftw_execute(pEx);
  fftw_execute(pEy);
  fftw_execute(pEz);
  fftw_execute(pBx);
  fftw_execute(pBy);
  fftw_execute(pBz);


  myComplex myI(0, 1);
  myI = myI*dt;


  double *wx = new double[n];
  double *wy = new double[m];
  double *wz = new double[h];

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      for (k = 0; k < h; k++)
      {

        wx[i] = 2 * M_PI*i / ((X_Max - X_Min)*n);
        wy[j] = 2 * M_PI*j / ((Y_Max - Y_Min)*m);
        wz[k] = 2 * M_PI*k / ((Z_Max - Z_Min)*h);
      }


  for (int l = 0; l < t; l++)
  {
    cout << "l: " << l << "/" << t << endl;

    for (i = 0; i < n; i++) //считаем ≈
      for (j = 0; j < m; j++)
        for (k = 0; k < h; k++)
        {
          Ex_c[i + j*n + k*n*m] = Ex_c[i + j*n + k*n*m] + myI*(Bz_c[i + j*n + k*n*m] * wy[j] - By_c[i + j*n + k*n*m] * wz[k]);
          Ey_c[i + j*n + k*n*m] = Ey_c[i + j*n + k*n*m] - myI*(Bz_c[i + j*n + k*n*m] * wx[i] - Bx_c[i + j*n + k*n*m] * wz[k]);
          Ez_c[i + j*n + k*n*m] = Ez_c[i + j*n + k*n*m] + myI*(By_c[i + j*n + k*n*m] * wx[i] - Bx_c[i + j*n + k*n*m] * wy[j]);
        }


    
    for (i = 0; i < n; i++) //считаем B
      for (j = 0; j < m; j++)
        for (k = 0; k < h; k++)
        {
          Bx_c[i + j*n + k*n*m] = Bx_c[i + j*n + k*n*m] - myI*(Ez_c[i + j*n + k*n*m] * wy[j] - Ey_c[i + j*n + k*n*m] * wz[k]);
          By_c[i + j*n + k*n*m] = By_c[i + j*n + k*n*m] + myI*(Ez_c[i + j*n + k*n*m] * wx[i] - Ex_c[i + j*n + k*n*m] * wz[k]);
          Bz_c[i + j*n + k*n*m] = Bz_c[i + j*n + k*n*m] - myI*(Ey_c[i + j*n + k*n*m] * wx[i] - Ex_c[i + j*n + k*n*m] * wy[j]);
        }

    //if (l % 2 == 0)
    //{
    //  string str = "PSTDout"; str += to_string(l); str += ".txt";
    //  ofstream fout(str);
    //  fftw_execute(pEx);
    //  fftw_execute(pEy);
    //  fftw_execute(pEz);
    //  fftw_execute(pBx);
    //  fftw_execute(pBy);
    //  fftw_execute(pBz);
    //  double size = n*m*h;
    //  for (i = 0; i < n; i++) //приводим
    //    for (j = 0; j < m; j++)
    //      for (k = 0; k < h; k++)
    //      {
    //        Ex[i + j*n + k*n*m] = Ex[i + j*n + k*n*m] / size;
    //        Ey[i + j*n + k*n*m] = Ey[i + j*n + k*n*m] / size;
    //        Ey[i + j*n + k*n*m] = Ey[i + j*n + k*n*m] / size;
    //        Bx[i + j*n + k*n*m] = Bx[i + j*n + k*n*m] / size;
    //        By[i + j*n + k*n*m] = By[i + j*n + k*n*m] / size;
    //        By[i + j*n + k*n*m] = By[i + j*n + k*n*m] / size;
    //      }
    //  for (i = 0; i < n; i++)
    //    for (j = 0; j < m; j++)
    //    {
    //      fout << Ey[i + j*n + h*n*m / 2] << "\n" << Bz[i + j*n + h*n*m / 2] << "\n";
    //    }
    //  fout.close();
    //}
  }


  //fftw_execute(pEx);
  //fftw_execute(pEy);
  //fftw_execute(pEz);
  //fftw_execute(pBx);
  //fftw_execute(pBy);
  //fftw_execute(pBz);


  //for (i = 0; i < n; i++) //приводим
  //  for (j = 0; j < m; j++)
  //    for (k = 0; k < h; k++)
  //    {

  //      Ex[i + j*n + k*n*m] = Ex[i + j*n + k*n*m] / n / m / h;
  //      Ey[i + j*n + k*n*m] = Ey[i + j*n + k*n*m] / n / m / h;
  //      Ey[i + j*n + k*n*m] = Ey[i + j*n + k*n*m] / n / m / h;

  //      Bx[i + j*n + k*n*m] = Bx[i + j*n + k*n*m] / n / m / h;
  //      By[i + j*n + k*n*m] = By[i + j*n + k*n*m] / n / m / h;
  //      By[i + j*n + k*n*m] = By[i + j*n + k*n*m] / n / m / h;

  //    }
  //fftw_destroy_plan(pEx);
  //fftw_destroy_plan(pEy);
  //fftw_destroy_plan(pEz);
  //fftw_destroy_plan(pBx);
  //fftw_destroy_plan(pBy);
  //fftw_destroy_plan(pBz);

  //pEx = fftw_plan_dft_c2r_3d(n, m, h, (fftw_complex *)Ex_c, (double *)Ex, FFTW_ESTIMATE);
  //pEy = fftw_plan_dft_c2r_3d(n, m, h, (fftw_complex *)Ey_c, (double *)Ey, FFTW_ESTIMATE);
  //pEz = fftw_plan_dft_c2r_3d(n, m, h, (fftw_complex *)Ez_c, (double *)Ez, FFTW_ESTIMATE);
  //pBx = fftw_plan_dft_c2r_3d(n, m, h, (fftw_complex *)Bx_c, (double *)Bx, FFTW_ESTIMATE);
  //pBy = fftw_plan_dft_c2r_3d(n, m, h, (fftw_complex *)By_c, (double *)By, FFTW_ESTIMATE);
  //pBz = fftw_plan_dft_c2r_3d(n, m, h, (fftw_complex *)Bz_c, (double *)Bz, FFTW_ESTIMATE);

  fftw_destroy_plan(pEx);
  fftw_destroy_plan(pEy);
  fftw_destroy_plan(pEz);
  fftw_destroy_plan(pBx);
  fftw_destroy_plan(pBy);
  fftw_destroy_plan(pBz);

  double *ex = new double[n*m*h];
  double *ey = new double[n*m*h];
  double *ez = new double[n*m*h];
  double *bx = new double[n*m*h];
  double *by = new double[n*m*h];
  double *bz = new double[n*m*h];

  pEx = fftw_plan_dft_c2r_3d(n, m, h, (fftw_complex *)Ex_c, (double *)ex, FFTW_MEASURE);
  pEy = fftw_plan_dft_c2r_3d(n, m, h, (fftw_complex *)Ey_c, (double *)ey, FFTW_MEASURE);
  pEz = fftw_plan_dft_c2r_3d(n, m, h, (fftw_complex *)Ez_c, (double *)ez, FFTW_MEASURE);
  pBx = fftw_plan_dft_c2r_3d(n, m, h, (fftw_complex *)Bx_c, (double *)bx, FFTW_MEASURE);
  pBy = fftw_plan_dft_c2r_3d(n, m, h, (fftw_complex *)By_c, (double *)by, FFTW_MEASURE);
  pBz = fftw_plan_dft_c2r_3d(n, m, h, (fftw_complex *)Bz_c, (double *)bz, FFTW_MEASURE);

  fftw_execute(pEx);
  fftw_execute(pEy);
  fftw_execute(pEz);
  fftw_execute(pBx);
  fftw_execute(pBy);
  fftw_execute(pBz);


  for (i = 0; i < n; i++) //приводим
    for (j = 0; j < m; j++)
      for (k = 0; k < h; k++)
      {

        ex[i + j*n + k*n*m] = ex[i + j*n + k*n*m] / n / m / h;
        ey[i + j*n + k*n*m] = ey[i + j*n + k*n*m] / n / m / h;
        ey[i + j*n + k*n*m] = ey[i + j*n + k*n*m] / n / m / h;

        bx[i + j*n + k*n*m] = bx[i + j*n + k*n*m] / n / m / h;
        by[i + j*n + k*n*m] = by[i + j*n + k*n*m] / n / m / h;
        by[i + j*n + k*n*m] = by[i + j*n + k*n*m] / n / m / h;
      }

  


  ofstream fout("PSTDout.txt");
  fout << n << endl;
  fout << m << endl;
  fout << h << endl;
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
    {
      fout << ey[i + j*n + h*n*m/2] << "\n" << bz[i + j*n + h*n*m / 2] << "\n";
    }
  fout.close();
//myComplex a(2, 3), b(7, 4), i(0, 1);
//cout << a.getRe() << " " << a.getIm() << endl;
//cout << b.getRe() << " " << b.getIm() << endl;
//cout << i.getRe() << " " << i.getIm() << endl;
//a = a + b*5*i;
//cout << a.getRe() << " " << a.getIm() << endl;
//cout << b.getRe() << " " << b.getIm() << endl;
//cout << i.getRe() << " " << i.getIm() << endl;
//b = b + i;
//cout << a.getRe() << " " << a.getIm() << endl;
//cout << b.getRe() << " " << b.getIm() << endl;
//cout << i.getRe() << " " << i.getIm() << endl;


}