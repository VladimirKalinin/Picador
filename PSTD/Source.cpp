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

  void set(double Re, double Im)
  {
    val[0] = Re;
    val[1] = Im;
  }

  void setRe(double Re)
  {
    val[0] = Re;
  }

  void getIm(double Im)
  {
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
  myComplex operator/(double a)
  {
    myComplex tmp(val[0] / a, val[1] / a);
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

  t = 2000;
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

  myComplex *Ey = new myComplex[n*m*h];
  myComplex *Ez = new myComplex[n*m*h];
  myComplex *Bx = new myComplex[n*m*h];
  myComplex *By = new myComplex[n*m*h];
  myComplex *Ex = new myComplex[n*m*h];
  myComplex *Bz = new myComplex[n*m*h];


  for (i = 0; i < n; i++)
  {
    for (j = 0; j < m; j++)
    {
      for (k = 0; k < h; k++)
      {
        x = X_Min + i*dx;
        y = Y_Min + j*dy;
        z = Z_Min + k*dz;
        Ex[k + h*(j + m*i)].set(0, 0);
        Ey[k + h*(j + m*i)].set(sin(4 * M_PI*(x - X_Min) / (X_Max - X_Min)), 0);
        //(((((exp(((((-((y * y) + (z * z))) * 2) * 0.69314718055994529) / 0.000025000000000000001)))) *
        //  (((exp((((-((((2 * x) / (-0.030000000000000002))) * (((2 * x) / (-0.030000000000000002))))) * 2) * 0.69314718055994529)) *
        //  (((0.5 * (1 + sign(((-x) - (-0.030000000000000002))))))))))) * cos((2617.9938779914942 * x))), 0);
        Ez[k + h*(j + m*i)].set(0, 0);

        Bx[k + h*(j + m*i)].set(0, 0);
        By[k + h*(j + m*i)].set(0, 0);
        Bz[k + h*(j + m*i)].set(sin(4 * M_PI*(x - X_Min) / (X_Max - X_Min)), 0);
        //(((((exp(((((-((y * y) + (z * z))) * 2) * 0.69314718055994529) / 0.000025000000000000001)))) *
        //  (((exp((((-((((2 * x) / (-0.030000000000000002))) * (((2 * x) / (-0.030000000000000002))))) * 2) * 0.69314718055994529)) *
        //  (((0.5 * (1 + sign(((-x) - (-0.030000000000000002))))))))))) * cos((2617.9938779914942 * x))), 0);
      }
    }
  }


  myComplex *Ex_c = new myComplex[n*m*h];
  myComplex *Ey_c = new myComplex[n*m*h];
  myComplex *Ez_c = new myComplex[n*m*h];
  myComplex *Bx_c = new myComplex[n*m*h];
  myComplex *By_c = new myComplex[n*m*h];
  myComplex *Bz_c = new myComplex[n*m*h];



  fftw_plan pEx = fftw_plan_dft_3d(n, m, h, (fftw_complex *)Ex, (fftw_complex *)Ex_c, FFTW_FORWARD, FFTW_MEASURE);
  fftw_plan pEy = fftw_plan_dft_3d(n, m, h, (fftw_complex *)Ey, (fftw_complex *)Ey_c, FFTW_FORWARD, FFTW_MEASURE);
  fftw_plan pEz = fftw_plan_dft_3d(n, m, h, (fftw_complex *)Ez, (fftw_complex *)Ez_c, FFTW_FORWARD, FFTW_MEASURE);
  fftw_plan pBx = fftw_plan_dft_3d(n, m, h, (fftw_complex *)Bx, (fftw_complex *)Bx_c, FFTW_FORWARD, FFTW_MEASURE);
  fftw_plan pBy = fftw_plan_dft_3d(n, m, h, (fftw_complex *)By, (fftw_complex *)By_c, FFTW_FORWARD, FFTW_MEASURE);
  fftw_plan pBz = fftw_plan_dft_3d(n, m, h, (fftw_complex *)Bz, (fftw_complex *)Bz_c, FFTW_FORWARD, FFTW_MEASURE);

  fftw_execute(pEx);
  fftw_execute(pEy);
  fftw_execute(pEz);
  fftw_execute(pBx);
  fftw_execute(pBy);
  fftw_execute(pBz);

  fftw_destroy_plan(pEx);
  fftw_destroy_plan(pEy);
  fftw_destroy_plan(pEz);
  fftw_destroy_plan(pBx);
  fftw_destroy_plan(pBy);
  fftw_destroy_plan(pBz);

  //double *ex = new double[n*m*h];
  //double *ey = new double[n*m*h];
  //double *ez = new double[n*m*h];
  //double *bx = new double[n*m*h];
  //double *by = new double[n*m*h];
  //double *bz = new double[n*m*h];

  pEx = fftw_plan_dft_3d(n, m, h, (fftw_complex *)Ex_c, (fftw_complex *)Ex, FFTW_BACKWARD, FFTW_MEASURE);
  pEy = fftw_plan_dft_3d(n, m, h, (fftw_complex *)Ey_c, (fftw_complex *)Ey, FFTW_BACKWARD, FFTW_MEASURE);
  pEz = fftw_plan_dft_3d(n, m, h, (fftw_complex *)Ez_c, (fftw_complex *)Ez, FFTW_BACKWARD, FFTW_MEASURE);
  pBx = fftw_plan_dft_3d(n, m, h, (fftw_complex *)Bx_c, (fftw_complex *)Bx, FFTW_BACKWARD, FFTW_MEASURE);
  pBy = fftw_plan_dft_3d(n, m, h, (fftw_complex *)By_c, (fftw_complex *)By, FFTW_BACKWARD, FFTW_MEASURE);
  pBz = fftw_plan_dft_3d(n, m, h, (fftw_complex *)Bz_c, (fftw_complex *)Bz, FFTW_BACKWARD, FFTW_MEASURE);

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < m; j++)
    {
      for (k = 0; k < h; k++)
      {
        Ex[k + h*(j + m*i)].set(0, 0);
        Ey[k + h*(j + m*i)].set(0, 0);
        Ez[k + h*(j + m*i)].set(0, 0);

        Bx[k + h*(j + m*i)].set(0, 0);
        By[k + h*(j + m*i)].set(0, 0);
        Bz[k + h*(j + m*i)].set(0, 0);
      }
    }
  }


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
          Ex_c[k + h*(j + m*i)] = Ex_c[k + h*(j + m*i)] + myI*(Bz_c[k + h*(j + m*i)] * wy[j] - By_c[k + h*(j + m*i)] * wz[k]);
          Ey_c[k + h*(j + m*i)] = Ey_c[k + h*(j + m*i)] - myI*(Bz_c[k + h*(j + m*i)] * wx[i] - Bx_c[k + h*(j + m*i)] * wz[k]);
          Ez_c[k + h*(j + m*i)] = Ez_c[k + h*(j + m*i)] + myI*(By_c[k + h*(j + m*i)] * wx[i] - Bx_c[k + h*(j + m*i)] * wy[j]);
        }


    
    for (i = 0; i < n; i++) //считаем B
      for (j = 0; j < m; j++)
        for (k = 0; k < h; k++)
        {
          Bx_c[k + h*(j + m*i)] = Bx_c[k + h*(j + m*i)] - myI*(Ez_c[k + h*(j + m*i)] * wy[j] - Ey_c[k + h*(j + m*i)] * wz[k]);
          By_c[k + h*(j + m*i)] = By_c[k + h*(j + m*i)] + myI*(Ez_c[k + h*(j + m*i)] * wx[i] - Ex_c[k + h*(j + m*i)] * wz[k]);
          Bz_c[k + h*(j + m*i)] = Bz_c[k + h*(j + m*i)] - myI*(Ey_c[k + h*(j + m*i)] * wx[i] - Ex_c[k + h*(j + m*i)] * wy[j]);
        }

    if (l % 20 == 0)
    {
      string str = "PSTDout"; str += to_string(l); str += ".txt";
      ofstream fout(str);
      fftw_execute(pEx);
      fftw_execute(pEy);
      fftw_execute(pEz);
      fftw_execute(pBx);
      fftw_execute(pBy);
      fftw_execute(pBz);
      double size = n*m*h;
      for (i = 0; i < n; i++) //приводим
        for (j = 0; j < m; j++)
          for (k = 0; k < h; k++)
          {
            Ex[k + h*(j + m*i)] = Ex[k + h*(j + m*i)] / size;
            Ey[k + h*(j + m*i)] = Ey[k + h*(j + m*i)] / size;
            Ez[k + h*(j + m*i)] = Ez[k + h*(j + m*i)] / size;
            Bx[k + h*(j + m*i)] = Bx[k + h*(j + m*i)] / size;
            By[k + h*(j + m*i)] = By[k + h*(j + m*i)] / size;
            Bz[k + h*(j + m*i)] = Bz[k + h*(j + m*i)] / size;
          }
      fout << n << endl;
      fout << m << endl;
      fout << h << endl;
      for (i = 0; i < n; i++)
        for (j = 0; j < m; j++)
        {
          fout << Ey[h/2 + h*(j + m*i)].getRe() << "\n" << Bz[h/2 + h*(j + m*i)].getRe() << "\n";
        }
      fout.close();
    }
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

  //      Ex[k + h*(j + m*i)] = Ex[k + h*(j + m*i)] / n / m / h;
  //      Ey[k + h*(j + m*i)] = Ey[k + h*(j + m*i)] / n / m / h;
  //      Ey[k + h*(j + m*i)] = Ey[k + h*(j + m*i)] / n / m / h;

  //      Bx[k + h*(j + m*i)] = Bx[k + h*(j + m*i)] / n / m / h;
  //      By[k + h*(j + m*i)] = By[k + h*(j + m*i)] / n / m / h;
  //      By[k + h*(j + m*i)] = By[k + h*(j + m*i)] / n / m / h;
  //    }

  


 /* ofstream fout("PSTDout.txt");
  fout << n << endl;
  fout << m << endl;
  fout << h << endl;
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
    {
      fout << Ey[k + h*(j + m*i)].getRe() << "\n" << Bz[k + h*(j + m*i)].getRe() << "\n";
    }
  fout.close();*/
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