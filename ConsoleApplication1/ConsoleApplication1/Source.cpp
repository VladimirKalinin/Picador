#include <stdio.h>
#include <vector>
#include <fstream>
#include <iostream>
using namespace std;

int main() {
  int n, m, h, i, j, k;
  double t, dt, dx, dy, dz, c = 299792458;

  ifstream fin("input.txt"); // открыли файл для чтения
    fin >> t;
    fin >> dt;
    fin >> dx;
    fin >> dy;
    fin >> dz; 
    fin >> n;
    fin >> m; 
    fin >> h; 

    double ****E;
    E = new double*** [n];
    for (i = 0; i < n; i++) {
      E[i] = new double**[m];
      for (j = 0; j < m; j++) {
        E[i][j] = new double*[h];
        for (k = 0; k < h; k++) {
          E[i][j][k] = new double [3];
          fin >> E[i][j][k][1];
          fin >> E[i][j][k][2];
          fin >> E[i][j][k][3];
        }
      }
    }

    double ****B;
    B = new double*** [n];
    for (i = 0; i < n; i++) {
      B[i] = new double**[m];
      for (j = 0; j < m; j++) {
        B[i][j] = new double*[h];
        for (k = 0; k < h; k++) {
          B[i][j][k] = new double [3];
          fin >> B[i][j][k][1];
          fin >> B[i][j][k][2];
          fin >> B[i][j][k][3];
        }
      }
    }
    fin.close(); // закрываем файл
    for (int l = 0; l < t; l++) {
      for (i = 0; i < n - 1; i++) //считаем Е
        for (j = 0; j < m - 1; j++)
          for (k = 0; k < h - 1; k++) {
            E[i][j][k][1] = E[i][j][k][1] + c*dt*((B[i][j + 1][k][3] - B[i][j][k][3])/dy - (B[i][j][k + 1][2] - B[i][j][k][2])/dz);
            E[i][j][k][2] = E[i][j][k][1] - c*dt*((B[i + 1][j][k][3] - B[i][j][k][3])/dx - (B[i][j][k + 1][1] - B[i][j][k][1])/dz);
            E[i][j][k][3] = E[i][j][k][1] + c*dt*((B[i + 1][j][k][2] - B[i][j][k][2])/dx - (B[i][j + 1][k][1] - B[i][j][k][1])/dy);
          }
      for (i = 0; i < n - 1; i++) //считаем B
        for (j = 0; j < m - 1; j++)
          for (k = 0; k < h - 1; k++) {
            B[i][j][k][1] = B[i][j][k][1] - c*dt*((E[i][j + 1][k][3] - E[i][j][k][3])/dy - (E[i][j][k + 1][2] - E[i][j][k][2])/dz);
            B[i][j][k][2] = B[i][j][k][1] + c*dt*((E[i + 1][j][k][3] - E[i][j][k][3])/dx - (E[i][j][k + 1][1] - E[i][j][k][1])/dz);
            B[i][j][k][3] = B[i][j][k][1] - c*dt*((E[i + 1][j][k][2] - E[i][j][k][2])/dx - (E[i][j + 1][k][1] - E[i][j][k][1])/dy);
          }
    }
  ofstream fout("output.txt"); 
  for (i = 0; i < n; i++) 
      for (j = 0; j < m; j++)
        for (k = 0; k < h; k++) {
          fout << E[i][j][k][1] << " " << E[i][j][k][2] << " " << E[i][j][k][3];
        }
  for (i = 0; i < n; i++) 
      for (j = 0; j < m; j++)
        for (k = 0; k < h; k++) {
          fout << B[i][j][k][1] << " " << B[i][j][k][2] << " " << B[i][j][k][3];
        }
    fout.close(); // закрываем файл
}
