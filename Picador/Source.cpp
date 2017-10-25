#include <stdio.h>
#include <vector>
#include <fstream>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <string>
//#include <fftw3.h>

using namespace std;


double sign(double x)
{
	if (x > 0.0)
		return 1;
	else
		return -1;
}


class vect {
public:
	double x, y, z;
	vect() {
		x = 0;
		y = 0;
		z = 0;
	}
	vect(double a, double b, double c) {
		x = a;
		y = b;
		z = c;
	}
	void set(double a, double b, double c) {
		x = a;
		y = b;
		z = c;
	}
	void operator=(vect a) {
		x = a.x;
		y = a.y;
		z = a.z;
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

  vect ***E;
  E = new vect**[n];
  for (i = 0; i < n; i++) {
    E[i] = new vect*[m];
    for (j = 0; j < m; j++) {
      E[i][j] = new vect[h];
      for (k = 0; k < h; k++)
      {
        x = X_Min + i*dx + 0.5*dx;
        y = Y_Min + j*dy + 0.5*dy;
        z = Z_Min + k*dz + 0.5*dz;
        E[i][j][k].set(0,
          ((((exp(((((-((y * y) + (z * z))) * 2) * 0.69314718055994529) / 0.000025000000000000001)))) *
          (((exp((((-((((2 * x) / (-0.030000000000000002))) * (((2 * x) / (-0.030000000000000002))))) * 2) * 0.69314718055994529)) *
            (((0.5 * (1 + sign(((-x) - (-0.030000000000000002))))))))))) * cos((2617.9938779914942 * x))),
          0);
      }
    }
  }

  vect ***B;
  B = new vect**[n];
  for (i = 0; i < n; i++) {
    B[i] = new vect*[m];
    for (j = 0; j < m; j++) {
      B[i][j] = new vect[h];
      for (k = 0; k < h; k++)
      {
        x = X_Min + i*dx;
        y = Y_Min + j*dy;
        z = Z_Min + k*dz;
        B[i][j][k].set(0, 0,
          ((((exp(((((-((y * y) + (z * z))) * 2) * 0.69314718055994529) / 0.000025000000000000001)))) *
          (((exp((((-((((2 * x) / (-0.030000000000000002))) * (((2 * x) / (-0.030000000000000002))))) * 2) * 0.69314718055994529)) *
            (((0.5 * (1 + sign(((-x) - (-0.030000000000000002))))))))))) * cos((2617.9938779914942 * x))));
      }

    }
  }

  //fin.close(); 

  for (int l = 1; l < t + 1; l++)
  {
    cout << "l: " << l << "/" << t << endl;

    int x, y, z;
    for (i = 0; i < n; i++) //считаем Е
      for (j = 0; j < m; j++)
        for (k = 0; k < h; k++)
        {
          if (i == n - 1) x = 0; else x = i + 1;
          if (j == m - 1) y = 0; else	y = j + 1;
          if (k == h - 1) z = 0; else z = k + 1;

          E[i][j][k].x = E[i][j][k].x + dt*((B[i][y][k].z - B[i][j][k].z) / dy - (B[i][j][z].y - B[i][j][k].y) / dz);
          E[i][j][k].y = E[i][j][k].y - dt*((B[x][j][k].z - B[i][j][k].z) / dx - (B[i][j][z].x - B[i][j][k].x) / dz);
          E[i][j][k].z = E[i][j][k].z + dt*((B[x][j][k].y - B[i][j][k].y) / dx - (B[i][y][k].x - B[i][j][k].x) / dy);
        }



    for (i = 0; i < n; i++) //считаем B
      for (j = 0; j < m; j++)
        for (k = 0; k < h; k++)
        {
          if (i == 0) x = n - 1; else x = i - 1;
          if (j == 0) y = m - 1; else	y = j - 1;
          if (k == 0) z = h - 1; else z = k - 1;

          B[i][j][k].x = B[i][j][k].x - dt*((E[i][j][k].z - E[i][y][k].z) / dy - (E[i][j][k].y - E[i][j][z].y) / dz);
          B[i][j][k].y = B[i][j][k].y + dt*((E[i][j][k].z - E[x][j][k].z) / dx - (E[i][j][k].x - E[i][j][z].x) / dz);
          B[i][j][k].z = B[i][j][k].z - dt*((E[i][j][k].y - E[x][j][k].y) / dx - (E[i][j][k].x - E[i][y][k].x) / dy);
        }

    if (l % 64 == 0)
    {
      string str = "input"; str += to_string(l); str += ".txt";
      ofstream fout(str);
      fout << n << endl;
      fout << m << endl;
      fout << h << endl;
      for (i = 0; i < n; i++)
        for (j = 0; j < m; j++)
        {
          fout << E[i][j][h/2].y << "\n" << B[i][j][h / 2].z << "\n";
        }
      fout.close();
    }

  }

  //ofstream fout("input.txt");
  //for (i = 0; i < n; i++)
  //  for (j = 0; j < m; j++)
  //  {
  //    fout <</* E[i][j][h/2].x << "\n" <<*/ B[i][j][h / 2].z << "\n" /*<< E[i][j][h/2].z << "\n" << B[i][j][h / 2].x << "\n" << B[i][j][h / 2].y << "\n" << B[i][j][h / 2].z << "\n"*/;
  //  }
  //fout.close();
}