#include <stdio.h>
#include <vector>
#include <fstream>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>

using namespace std;

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
  int n, m, h, i, j, k;
  double t, dt, dx, dy, dz, c = 299792458;

  //ifstream fin("input.txt"); // открыли файл для чтения
  //  fin >> t;
  //  fin >> dt;
  //  fin >> dx;
  //  fin >> dy;
  //  fin >> dz; 
  //  fin >> n;
  //  fin >> m; 
  //  fin >> h; 

  t = 200;
  dt = 1/(640*c);
  n = 128;
  m = 4;
  h = 4;
  dx = 1/64.0;
  dy = 1/2.0;
  dz = 1/2.0;

  vect ***E;
  E = new vect** [n];
  for (i = 0; i < n; i++) {
    E[i] = new vect*[m];
    for (j = 0; j < m; j++) {
      E[i][j] = new vect[h];
      for (k = 0; k < h; k++) 
        E[i][j][k].set(0, sin(4.0 * M_PI * (i*dx+dx/2.0)), 0);
    }
  }

  vect ***B;
  B = new vect** [n];
  for (i = 0; i < n; i++) {
    B[i] = new vect*[m];
    for (j = 0; j < m; j++) {
      B[i][j] = new vect[h];
      for (k = 0; k < h; k++) 
        B[i][j][k].set(0, 0, sin(4.0 * M_PI * i*dx));
    }
  }

  //fin.close(); // закрываем файл

  for (int l = 0; l < t; l++) {
    for (i = 0; i < n - 1; i++) //считаем Е
      for (j = 0; j < m - 1; j++)
        for (k = 0; k < h - 1; k++) {
          E[i][j][k].x = E[i][j][k].x + c*dt*((B[i][j + 1][k].z - B[i][j][k].z)/dy - (B[i][j][k + 1].y - B[i][j][k].y)/dz);
          E[i][j][k].y = E[i][j][k].y - c*dt*((B[i + 1][j][k].z - B[i][j][k].z)/dx - (B[i][j][k + 1].x - B[i][j][k].x)/dz);
          E[i][j][k].z = E[i][j][k].z + c*dt*((B[i + 1][j][k].y - B[i][j][k].y)/dx - (B[i][j + 1][k].x - B[i][j][k].x)/dy);
        }
    i = n-1;
    for (j = 0; j < m - 1; j++)
      for (k = 0; k < h - 1; k++) {
        E[i][j][k].x = E[i][j][k].x + c*dt*((B[i][j + 1][k].z - B[i][j][k].z)/dy - (B[i][j][k + 1].y - B[i][j][k].y)/dz);
        E[i][j][k].y = E[i][j][k].y - c*dt*((B[0][j][k].z - B[i][j][k].z)/dx - (B[i][j][k + 1].x - B[i][j][k].x)/dz);
        E[i][j][k].z = E[i][j][k].z + c*dt*((B[0][j][k].y - B[i][j][k].y)/dx - (B[i][j + 1][k].x - B[i][j][k].x)/dy);
      }
    j = m-1;
    for (i = 0; i < n - 1; i++)
      for (k = 0; k < h - 1; k++) {
        E[i][j][k].x = E[i][j][k].x + c*dt*((B[i][0][k].z - B[i][j][k].z)/dy - (B[i][j][k + 1].y - B[i][j][k].y)/dz);
        E[i][j][k].y = E[i][j][k].y - c*dt*((B[i + 1][j][k].z - B[i][j][k].z)/dx - (B[i][j][k + 1].x - B[i][j][k].x)/dz);
        E[i][j][k].z = E[i][j][k].z + c*dt*((B[i + 1][j][k].y - B[i][j][k].y)/dx - (B[i][0][k].x - B[i][j][k].x)/dy);
      }
    k = h-1;
    for (i = 0; i < n - 1; i++)
      for (j = 0; j < m - 1; j++) {
        E[i][j][k].x = E[i][j][k].x + c*dt*((B[i][j + 1][k].z - B[i][j][k].z)/dy - (B[i][j][0].y - B[i][j][k].y)/dz);
        E[i][j][k].y = E[i][j][k].y - c*dt*((B[i + 1][j][k].z - B[i][j][k].z)/dx - (B[i][j][0].x - B[i][j][k].x)/dz);
        E[i][j][k].z = E[i][j][k].z + c*dt*((B[i + 1][j][k].y - B[i][j][k].y)/dx - (B[i][j + 1][k].x - B[i][j][k].x)/dy);
      }
    i = n-1;
    j = m-1;
    for (k = 0; k < h - 1; k++) {
      E[i][j][k].x = E[i][j][k].x + c*dt*((B[i][0][k].z - B[i][j][k].z)/dy - (B[i][j][k + 1].y - B[i][j][k].y)/dz);
      E[i][j][k].y = E[i][j][k].y - c*dt*((B[0][j][k].z - B[i][j][k].z)/dx - (B[i][j][k + 1].x - B[i][j][k].x)/dz);
      E[i][j][k].z = E[i][j][k].z + c*dt*((B[0][j][k].y - B[i][j][k].y)/dx - (B[i][0][k].x - B[i][j][k].x)/dy);
    }
     j = m - 1;
     k = h - 1;
    for (i = 0; i < n - 1; i++) {
      E[i][j][k].x = E[i][j][k].x + c*dt*((B[i][0][k].z - B[i][j][k].z)/dy - (B[i][j][0].y - B[i][j][k].y)/dz);
      E[i][j][k].y = E[i][j][k].y - c*dt*((B[i + 1][j][k].z - B[i][j][k].z)/dx - (B[i][j][0].x - B[i][j][k].x)/dz);
      E[i][j][k].z = E[i][j][k].z + c*dt*((B[i + 1][j][k].y - B[i][j][k].y)/dx - (B[i][0][k].x - B[i][j][k].x)/dy);
    }
    i = n - 1;
    k = h - 1;
    for (j = 0; j < m - 1; j++) {
          E[i][j][k].x = E[i][j][k].x + c*dt*((B[i][j + 1][k].z - B[i][j][k].z)/dy - (B[i][j][0].y - B[i][j][k].y)/dz);
          E[i][j][k].y = E[i][j][k].y - c*dt*((B[0][j][k].z - B[i][j][k].z)/dx - (B[i][j][0].x - B[i][j][k].x)/dz);
          E[i][j][k].z = E[i][j][k].z + c*dt*((B[0][j][k].y - B[i][j][k].y)/dx - (B[i][j + 1][k].x - B[i][j][k].x)/dy);
        }
    i = n - 1;
    j = m - 1;
    k = h - 1;
    E[i][j][k].x = E[i][j][k].x + c*dt*((B[i][0][k].z - B[i][j][k].z)/dy - (B[i][j][0].y - B[i][j][k].y)/dz);
    E[i][j][k].y = E[i][j][k].y - c*dt*((B[0][j][k].z - B[i][j][k].z)/dx - (B[i][j][0].x - B[i][j][k].x)/dz);
    E[i][j][k].z = E[i][j][k].z + c*dt*((B[0][j][k].y - B[i][j][k].y)/dx - (B[i][0][k].x - B[i][j][k].x)/dy);

    for (i = 1; i < n; i++) //считаем B
      for (j = 1; j < m; j++)
        for (k = 1; k < h; k++) {
          B[i][j][k].x = B[i][j][k].x - c*dt*((E[i][j][k].z - E[i][j - 1][k].z)/dy - (E[i][j][k].y - E[i][j][k - 1].y)/dz);
          B[i][j][k].y = B[i][j][k].y + c*dt*((E[i][j][k].z - E[i - 1][j][k].z)/dx - (E[i][j][k].x - E[i][j][k - 1].x)/dz);
          B[i][j][k].z = B[i][j][k].z - c*dt*((E[i][j][k].y - E[i - 1][j][k].y)/dx - (E[i][j][k].x - E[i][j - 1][k].x)/dy);
        }
    i = 0;
    for (j = 1; j < m; j++)
      for (k = 1; k < h; k++) {
        B[i][j][k].x = B[i][j][k].x - c*dt*((E[i][j - 1][k].z - E[i][j][k].z)/dy - (E[i][j][k - 1].y - E[i][j][k].y)/dz);
        B[i][j][k].y = B[i][j][k].y + c*dt*((E[n - 1][j][k].z - E[i][j][k].z)/dx - (E[i][j][k - 1].x - E[i][j][k].x)/dz);
        B[i][j][k].z = B[i][j][k].z - c*dt*((E[n - 1][j][k].y - E[i][j][k].y)/dx - (E[i][j - 1][k].x - E[i][j][k].x)/dy);
      }
    j = 0;
    for (i = 1; i < n; i++)
      for (k = 1; k < h; k++) {
        B[i][j][k].x = B[i][j][k].x - c*dt*((E[i][m - 1][k].z - E[i][j][k].z)/dy - (E[i][j][k - 1].y - E[i][j][k].y)/dz);
        B[i][j][k].y = B[i][j][k].y + c*dt*((E[i - 1][j][k].z - E[i][j][k].z)/dx - (E[i][j][k - 1].x - E[i][j][k].x)/dz);
        B[i][j][k].z = B[i][j][k].z - c*dt*((E[i - 1][j][k].y - E[i][j][k].y)/dx - (E[i][m - 1][k].x - E[i][j][k].x)/dy);
      }
    k = 0;
    for (i = 1; i < n; i++)
      for (j = 1; j < m; j++) {
        B[i][j][k].x = B[i][j][k].x - c*dt*((E[i][j - 1][k].z - E[i][j][k].z)/dy - (E[i][j][h - 1].y - E[i][j][k].y)/dz);
        B[i][j][k].y = B[i][j][k].y + c*dt*((E[i - 1][j][k].z - E[i][j][k].z)/dx - (E[i][j][h - 1].x - E[i][j][k].x)/dz);
        B[i][j][k].z = B[i][j][k].z - c*dt*((E[i - 1][j][k].y - E[i][j][k].y)/dx - (E[i][j - 1][k].x - E[i][j][k].x)/dy);
      }
    i = 0;
    j = 0;
    for (k = 1; k < h; k++) {
      B[i][j][k].x = B[i][j][k].x - c*dt*((E[i][m - 1][k].z - E[i][j][k].z)/dy - (E[i][j][k - 1].y - E[i][j][k].y)/dz);
      B[i][j][k].y = B[i][j][k].y + c*dt*((E[n - 1][j][k].z - E[i][j][k].z)/dx - (E[i][j][k - 1].x - E[i][j][k].x)/dz);
      B[i][j][k].z = B[i][j][k].z - c*dt*((E[n - 1][j][k].y - E[i][j][k].y)/dx - (E[i][m - 1][k].x - E[i][j][k].x)/dy);
    }
     j = 0;
     k = 0;
    for (i = 1; i < n; i++) {
      B[i][j][k].x = B[i][j][k].x - c*dt*((E[i][m - 1][k].z - E[i][j][k].z)/dy - (E[i][j][h - 1].y - E[i][j][k].y)/dz);
      B[i][j][k].y = B[i][j][k].y + c*dt*((E[i - 1][j][k].z - E[i][j][k].z)/dx - (E[i][j][h - 1].x - E[i][j][k].x)/dz);
      B[i][j][k].z = B[i][j][k].z - c*dt*((E[i - 1][j][k].y - E[i][j][k].y)/dx - (E[i][m - 1][k].x - E[i][j][k].x)/dy);
    }
    i = 0;
    k = 0;
    for (j = 1; j < m; j++) {
      B[i][j][k].x = B[i][j][k].x - c*dt*((E[i][j - 1][k].z - E[i][j][k].z)/dy - (E[i][j][h - 1].y - E[i][j][k].y)/dz);
      B[i][j][k].y = B[i][j][k].y + c*dt*((E[n - 1][j][k].z - E[i][j][k].z)/dx - (E[i][j][h - 1].x - E[i][j][k].x)/dz);
      B[i][j][k].z = B[i][j][k].z - c*dt*((E[n - 1][j][k].y - E[i][j][k].y)/dx - (E[i][j - 1][k].x - E[i][j][k].x)/dy);
        }
    i = 0;
    j = 0;
    k = 0;
    E[i][j][k].x = E[i][j][k].x - c*dt*((E[i][m - 1][k].z - E[i][j][k].z)/dy - (E[i][j][h - 1].y - E[i][j][k].y)/dz);
    E[i][j][k].y = E[i][j][k].y + c*dt*((E[n - 1][j][k].z - E[i][j][k].z)/dx - (E[i][j][h - 1].x - E[i][j][k].x)/dz);
    E[i][j][k].z = E[i][j][k].z - c*dt*((E[n - 1][j][k].y - E[i][j][k].y)/dx - (E[i][m - 1][k].x - E[i][j][k].x)/dy);
  }

  ofstream fout("output.txt"); 
  for (i = 0; i < n; i++) 
     {
       fout << E[i][1][1].x << " " << E[i][0][0].y << /*" " << E[i][j][k].z <<*/ "\n";
      }
      //for (i = 0; i < n; i++) 
      //  for (j = 0; j < m; j++)
      //    for (k = 0; k < h; k++)     {
      //      fout << B[i][j][k][1] << " " << B[i][j][k][2] << " " << B[i][j][k][3] << "\n";
      //    }
          fout.close(); // закрываем файл
}