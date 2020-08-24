/*
  3D renderer
  Last modified 23th of August 2020
  Copyright (c) 2020 Arthur Ferreira (arthur.ferreira2@gmail.com)

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.
*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <SDL2/SDL.h>

#define WIDTH 400
#define HEIGHT 400
#define FPS 60               // frames per second

#define X 0                  // arrays indices for better readability
#define Y 1
#define Z 2

#define WIRE 0               // different modes of rendering
#define LINE 1
#define SOLID 2
#define GOURAUD 3
#define PHONG 4

#define AMBIANT_LIGHT 0.1    // Ambiant Light for Phong rendering

typedef struct point{
  int x, y;
} point;

typedef float sommet[4];     // un sommet est un tableau de 3 coord. + nb faces auquelles il appartient
typedef float normale[3];    // pour chaque sommet, chaque face, on pre-calcule la normale
typedef int face[3];         // une face est un tableau de 3 points

struct objet{                // to store all data comming from the file
  int     nb_sommets;        // and pre-calculated normales
  int     nb_faces;
  sommet  *sommets;
  normale *Pnormales;        // normales to points (precalculees)
  face    *faces;
  normale *Fnormales;        // normales to faces (precalculees)
};



//    DECLARATION DE VARIABLES GLOBALES

Uint32* Image;             // our frame buffer
float *Zb = NULL;          // Z-Buffer

float DObj = 0;            // distance objet
float DObs = 200;          // distance observateur, facteur de prespective

int areaL = WIDTH;         // window size
int areaH = HEIGHT;
int MaxX = WIDTH, MaxY = HEIGHT, MinX = 0, MinY = 0;

float luxX = -1000;        // position of the light source
float luxY = -1000;
float luxZ = 1000;

Uint32 Palette[256] = {0};
int nb_colors = 256;



// Cette fonction affiche un point dans une ximage et met a jour le Z-Buffer
void Point(int x, int y, int C, float z){
  if ((x >=  0) && (x < areaL) && (y >=  0) && (y < areaH) && (Zb[x + (y * areaL)] < z)){
     Image[x + (y * areaL)] = (int)C;
     Zb[x + (y * areaL)] = z;
     if (x >= MaxX) MaxX = x + 1;
     if (y >= MaxY) MaxY = y + 1;
     if (x <  MinX) MinX = x;
     if (y <  MinY) MinY = y;
  }
}

// Echange de deux entiers
inline void swapi(int *a, int *b){
  int t = *a;
  *a = *b;
  *b = t;
}

// Echange de deux decimaux
inline void swapf(float *a, float *b){
  float t = *a;
  *a = *b;
  *b = t;
}

//  Trace une droite
void DDA(int x1, int y1, float z1, int x2, int y2, float z2, float C){
  int Color = Palette[(unsigned char)(C * (nb_colors - 1))];
  int x, y;
  float dx, dy, dz, z, m, b;

  if (x2 == x1){             // verticale
    if (y1 > y2){
      swapi(&y1, &y2);
      swapf(&z1, &z2);
    }
    dz = (z2-z1) / (y2-y1);
    z = z1;
    for (y = y1; y<= y2; y++){
      Point(x1, y, Color, z);
      z += dz;
    }
    return;
  }

  if (x1 > x2){              // force le premier point a etre a gauche
    swapi(&x1, &x2);
    swapi(&y1, &y2);
    swapf(&z1, &z2);
  }

  dx = x2 - x1;
  dy = y2 - y1;

  if (fabs(dx) >= fabs(dy)){  // si la pente est entre -1 et 1
    m = dy / dx;
    b = y1 - (m * x1);
    dz = (z2 - z1) / (x2 - x1);
    z = z1;
    for (x = x1;x<= x2;x++){
      y = m * x + b;
      Point(x, y, Color, z);
      z += dz;
    }
  }
  else{                      // dans les autres cas
    if (y1 > y2){            // force le premier point a etre en haut
      swapi(&x1, &x2);
      swapi(&y1, &y2);
      swapf(&z1, &z2);
    }
    m = dx / dy;
    b = x1 - (m * y1);
    dz = (z2 - z1) / (y2 - y1);
    z = z1;
    for (y = y1; y<= y2; y++){
      x = m * y + b;
      Point(x, y, Color, z);
      z += dz;
    }
  }
}

// Trace 3 droites (constituant une face) en "clippant" chacune d'elles
void lines(point *P, float C, float *zed){
  float x1, x2, y1, y2, z1, z2, dx, dy;
  int C1, C2, i;

  for (i = 0;i<= 2;i++){     // Pour chaque droite
    x1 = P[i].x;
    y1 = P[i].y;
    z1 = zed[i];
   if (i == 2){
      x2 = P[0].x;
      y2 = P[0].y;
      z2 = zed[0];
    }
    else{
      x2 = P[i + 1].x;
      y2 = P[i + 1].y;
      z2 = zed[i + 1];
    }
    // Clipping
    C1 = ((x1 < 0) << 3) | ((x1 > areaL) << 2) | ((y1 < 0) << 1) | ((y1 > areaH));
    C2 = ((x2 < 0) << 3) | ((x2 > areaL) << 2) | ((y2 < 0) << 1) | ((y2 > areaH));

    if (C1 & C2) return;

    while(C1 | C2){          // tant qu'il a qqch de visible
      dx = x2 - x1;
      dy = y2 - y1;

      if (C1){
        if (C1 & 8)      { y1 += dy * (0 - x1) / dx;     x1 = 0;     }
        else if (C1 & 4) { y1 += dy * (areaL-x1) / dx;   x1 = areaL; }
        else if (C1 & 2) { x1 += dx * (0 - y1) / dy;     y1 = 0;     }
        else if (C1 & 1) { x1 += dx * (areaH - y1) / dy; y1 = areaH; }
        C1 = ((x1 < 0) << 3) | ((x1 > areaL) << 2) | ((y1 < 0) << 1) | ((y1 > areaH));
      }
      else{
        if (C2 & 8)      { y2 += dy * (0 - x2) / dx;     x2 = 0;     }
        else if (C2 & 4) { y2 += dy * (areaL-x2) / dx;   x2 = areaL; }
        else if (C2 & 2) { x2 += dx * (0 - y2) / dy;     y2 = 0;     }
        else if (C2 & 1) { x2 += dx * (areaH - y2) / dy; y2 = areaH; }
        C2 = ((x2 < 0) << 3) | ((x2 > areaL) << 2) | ((y2 < 0) << 1) | ((y2 > areaH));
      }
    }
    DDA((int)x1, (int)y1, z1, (int)x2, (int)y2, z2, C);     // Tracage
  }
}

//  Echange deux points ayant 3 coordonnees
inline void swap_flat(float *x1, float *y1, float *z1, float *x2, float *y2, float *z2){
  float a = *x1;
  float b = *y1;
  float c = *z1;

  *x1 = *x2;
  *y1 = *y2;
  *z1 = *z2;
  *x2 = a;
  *y2 = b;
  *z2 = c;
}

//  Remplis une face triangulaire d'une couleur C et en mettant a jour le Z-Buffer.
//  Le calcul des profondeurs de chaque point se fait par interpolation
//  des z de chaque sommet grace a zi zf dzi dzf et dz
void flat(point *P, float C, float *zed){
  float  x1, x2, x3, y1, y2, y3, z1, z2, z3,  xi, xf, dx1, dx2,  z, dz, zi, dzi, zf, dzf;
  int  Passe = 0, x, y;
  int Color = Palette[(unsigned char)(C*(nb_colors - 1))];

  x1 = P[0].x; y1 = P[0].y; z1 = zed[0];
  x2 = P[1].x; y2 = P[1].y; z2 = zed[1];
  x3 = P[2].x; y3 = P[2].y; z3 = zed[2];

  if (y1 > y2)
    swap_flat(&x1, &y1, &z1, &x2, &y2, &z2);

  if (y1 > y3)
    swap_flat(&x1, &y1, &z1, &x3, &y3, &z3);         // x1,y1 est le plus ou egal haut des deux autres

  if (x2 > x3)
    swap_flat(&x2, &y2, &z2, &x3, &y3, &z3);         // x2,y2 est a gauche de x3,y3

  if (y1 == y2){
    xi = x1;
    zi = z1;
    y = floor(y1);
    xf = x2;
    zf = z2;
    Passe = 1;
    if (y3 != y1){
      dx1 = (x3 - x1) / (y3 - y1);
      dx2 = (x3 - x2) / (y3 - y2);
      dzi = (z3 - z1) / (y3 - y1);
      dzf = (z3 - z2) / (y3 - y2);
    }
    else
      dx1 = dx2 = dzi = dzf = 0;
  }
  else{
    if (y1 == y3){
      xi = x1;
      zi = z1;
      y = floor(y1);
      xf = x3;
      zf = z3;
      Passe = 1;
      if (y2!= y1){
        dx1 = (x2 - x1) / (y2 - y1);
        dx2 = (x2 - x3) / (y2 - y3);
        dzi = (z2 - z1) / (y2 - y1);
        dzf = (z2 - z3) / (y2 - y3);
      }
      else
       dx1 = dx2 = dzi = dzf = 0;
    }
    else{
      xi = xf = x1;
      zi = zf = z1;
      y = floor(y1);
      dx1 = (x2 - x1) / (y2 - y1);
      dx2 = (x3 - x1) / (y3 - y1);
      dzi = (z2 - z1) / (y2 - y1);
      dzf = (z3 - z1) / (y3 - y1);
    }
  }

  while ((y < y2) || (y < y3)){
    if (xf != xi) dz = (zf - zi) / (xf - xi);
    else dz = 0;
    z = zi;

    if (xi <= xf)
      for (x=(int)xi; x<= (int)xf; x++){
        Point(x, y, Color, z);
        z += dz;
      }
    else
      for (x = (int)xi; x>= (int)xf; x--){
        Point(x, y, Color, z);
        z -= dz;
      }

    y++;

    if(!Passe && (y > y2)){
      if (y3 != y2){
        dx1 = (x3 - x2) / (y3 - y2);
        dzi = (z3 - z2) / (y3 - y2);
      }
      else{
        dx1 = 0;
        dzi = 0;
      }
      Passe = 1;
    }
    if(!Passe && (y > y3)){
      if (y2 != y3){
        dx2 = (x2 - x3) / (y2 - y3);
        dzf = (z2 - z3) / (y2 - y3);
      }
      else{
        dx2 = 0;
        dzf = 0;
      }
      Passe = 1;
    }
    xi += dx1;
    xf += dx2;
    zi += dzi;
    zf += dzf;
  }
}

//  Echange deux points ayant 3 coordonnees ainsi que leur couleur
inline void swap_gouraud(float *x1,float *y1,float *z1,float *c1,float *x2,float *y2,float *z2,float *c2){
  float a = *x1;
  float b = *y1;
  float c = *z1;
  float d = *c1;

  *x1 = *x2;
  *y1 = *y2;
  *z1 = *z2;
  *c1 = *c2;
  *x2 = a;
  *y2 = b;
  *z2 = c;
  *c2 = d;
}

//  Remplis une face triangulaire par interpolation des couleurs
//  et en mettant a jour le Z-Buffer
void gouraud(point *P, float *C, float *zed){
  float x1, x2, x3, y1, y2, y3, z1, z2, z3, xi, xf,  z, dz, zi, dzi, zf, dzf;
  float dx1, dx2, Ci, Cf, dC1, dC2, c1, c2, c3, Coul, dC;
  int x, y, Passe = 0, Color;

  x1 = P[0].x; y1 = P[0].y; c1 = C[0]; z1 = zed[0];
  x2 = P[1].x; y2 = P[1].y; c2 = C[1]; z2 = zed[1];
  x3 = P[2].x; y3 = P[2].y; c3 = C[2]; z3 = zed[2];


  if (y1 > y2)
    swap_gouraud(&x1, &y1, &z1, &c1, &x2, &y2, &z2, &c2);

  if (y1 > y3)
    swap_gouraud(&x1, &y1, &z1, &c1, &x3, &y3, &z3, &c3);

  if (x2 > x3)
    swap_gouraud(&x2, &y2, &z2, &c2, &x3, &y3, &z3, &c3);

  if (y1 == y2){
    Ci = c1;
    xi = x1;
    zi = z1;
    y = floor(y1);
    xf = x2;
    Cf = c2;
    zf = z2;
    Passe = 1;
    if (y3 != y1){
      dx1 = (x3 - x1) / (y3 - y1);
      dx2 = (x3 - x2) / (y3 - y2);
      dC1 = (c3 - c1) / (y3 - y1);
      dC2 = (c3 - c2) / (y3 - y2);
      dzi = (z3 - z1) / (y3 - y1);
      dzf = (z3 - z2) / (y3 - y2);
    }
    else
      dC1 = dC2 = dx1 = dx2 = dzi = dzf = 0;
  }
  else{
    if (y1 == y3){
      Ci = c1;
      xi = x1;
      zi = z1;
      y = floor(y1);
      xf = x3;
      Cf = c3;
      zf = z3;
      Passe = 1;
      if (y2!= y1){
        dx1 = (x2 - x1) / (y2 - y1);
        dx2 = (x2 - x3) / (y2 - y3);
        dC1 = (c2 - c1) / (y2 - y1);
        dC2 = (c2 - c3) / (y2 - y3);
        dzi = (z2 - z1) / (y2 - y1);
        dzf = (z2 - z3) / (y2 - y3);
      }
      else
       dC1 = dC2 = dx1 = dx2 = dzi = dzf = 0;
    }
    else{
      Ci = Cf = c1;
      xi = xf = x1;
      zi = zf = z1;
      y = floor(y1);
      dx1 = (x2 - x1) / (y2 - y1);
      dx2 = (x3 - x1) / (y3 - y1);
      dC1 = (c2 - c1) / (y2 - y1);
      dC2 = (c3 - c1) / (y3 - y1);
      dzi = (z2 - z1) / (y2 - y1);
      dzf = (z3 - z1) / (y3 - y1);
    }
  }

  while ((y < y2) || (y < y3)){
    if (xf != xi)
      dC = (Cf - Ci) / (xf - xi);
    else
      dC = 0;

    Coul = Ci;
    if (xf != xi) dz = (zf - zi) / (xf - xi);
    else dz = 0;
    z = zi;

    if (xi <= xf)
      for (x = (int)xi; x<=(int)xf; x++){
        Color = (int)((Coul + 0.1) * nb_colors);
        if (Color < 0) Color = 0;
        if (Color >= nb_colors) Color = nb_colors - 1;
        Point(x, y, Palette[Color], z);
        Coul += dC;
        z += dz;
      }
    else
      for (x = (int)xi; x>=(int)xf; x--){
        Color = (int)((Coul + 0.1) * nb_colors);
        if (Color < 0) Color = 0;
        if (Color >= nb_colors) Color = nb_colors - 1;
        Point(x, y, Palette[Color], z);
        Coul -= dC;
        z -= dz;
      }

    y++;

    if(!Passe && (y > y2)){
      if (y3 != y2){
        dx1 = (x3 - x2) / (y3 - y2);
        dC1 = (c3 - c2) / (y3 - y2);
        dzi = (z3 - z2) / (y3 - y2);
      }
      else
        dx1 = dzi = dC1 = 0;
      Passe = 1;
    }
    if(!Passe && (y > y3)){
      if (y2!= y3){
        dx2 = (x2 - x3) / (y2 - y3);
        dC2 = (c2 - c3) / (y2 - y3);
        dzf = (z2 - z3) / (y2 - y3);
      }
      else
        dx2 = dzf = dC2 = 0;
      Passe = 1;
    }
    xi += dx1;
    xf += dx2;
    Ci += dC1;
    Cf += dC2;
    zi += dzi;
    zf += dzf;
  }
}

//  Echange deux points ayant 3 coordonnees ainsi que leurs normales
inline void swap_phong(float *x1, float *y1, float *z1, float *N1x, float *N1y, float *N1z, float *x2, float *y2, float *z2, float *N2x, float *N2y, float *N2z){
  float a = *x1;
  float b = *y1;
  float c = *z1;
  float d = *N1x;
  float e = *N1y;
  float f = *N1z;

  *x1 = *x2;
  *y1 = *y2;
  *z1 = *z2;
  *N1x = *N2x;
  *N1y = *N2y;
  *N1z = *N2z;

  *x2 = a;
  *y2 = b;
  *z2 = c;
  *N2x = d;
  *N2y = e;
  *N2z = f;
}

//  Remplis une face triangulaire par interpolation des normales
//  et en mettant a jour le Z-Buffer
void phong(point *P,float *N,float *zed){

  float x1, x2, x3, y1, y2, y3, z1, z2, z3, xi, xf, z, dz, zi, dzi, zf, dzf;
  float dx1, dx2, Nix, Niy, Niz, Nfx, Nfy, Nfz, Nx, Ny, Nz, DY;
  float dN1x, dN1y, dN1z, dN2x, dN2y, dN2z, dNx, dNy, dNz, N1x, N1y, N1z, N2x, N2y, N2z, N3x, N3y, N3z;
  int x, y, Passe = 0, Color;
  float C, cosAlpha;
  float NormeLux = sqrt(luxX * luxX + luxY * luxY + luxZ * luxZ);

  x1 = P[0].x; y1 = P[0].y; N1x = N[0]; N1y = N[1]; N1z = N[2]; z1 = zed[0];
  x2 = P[1].x; y2 = P[1].y; N2x = N[3]; N2y = N[4]; N2z = N[5]; z2 = zed[1];
  x3 = P[2].x; y3 = P[2].y; N3x = N[6]; N3y = N[7]; N3z = N[8]; z3 = zed[2];


  if (y1 > y2)
    swap_phong(&x1, &y1, &z1, &N1x, &N1y, &N1z, &x2, &y2, &z2, &N2x, &N2y, &N2z);

  if (y1 > y3)
    swap_phong(&x1, &y1, &z1, &N1x, &N1y, &N1z, &x3, &y3, &z3, &N3x, &N3y, &N3z);

  if (x2 > x3)
    swap_phong(&x2, &y2, &z2, &N2x, &N2y, &N2z, &x3, &y3, &z3, &N3x, &N3y, &N3z);

  if (y1 == y2){
    Nix = N1x;
    Niy = N1y;
    Niz = N1z;
    xi = x1;
    zi = z1;
    y = floor(y1);
    xf = x2;
    Nfx = N2x;
    Nfy = N2y;
    Nfz = N2z;
    zf = z2;
    Passe = 1;
    if (y3 != y1){
      DY   = y3 - y1;
      dx1  = (x3 - x1) / DY;
      dN1x = (N3x - N1x) / DY;
      dN1y = (N3y - N1y) / DY;
      dN1z = (N3z - N1z) / DY;
      dzi  = (z3 - z1) / DY;

      DY   = y3 - y2;
      dx2  = (x3 - x2) / DY;
      dN2x = (N3x - N2x) / DY;
      dN2y = (N3y - N2y) / DY;
      dN2z = (N3z - N2z) / DY;
      dzf  = (z3 - z2) / DY;
    }
    else
      dN1x = dN1y = dN1z = dN2x = dN2y = dN2z = dx1 = dx2 = dzi = dzf = 0;
  }
  else{
    if (y1 == y3){
      Nix = N1x;
      Niy = N1y;
      Niz = N1z;
      xi = x1;
      zi = z1;
      y = floor(y1);
      xf = x3;
      Nfx = N3x;
      Nfy = N3y;
      Nfz = N3z;
      zf = z3;
      Passe = 1;
      if (y2 != y1){
        DY   = y2 - y1;
        dx1  = (x2 - x1) / DY;
        dx2  = (x2 - x3) / (y2 - y3);
        dN1x = (N2x - N1x) / DY;
        dN1y = (N2y - N1y) / DY;
        dN1z = (N2z - N1z) / DY;
        dzi  = (z2 - z1) / DY;

        DY   = y2 - y3;
        dx2  = (x2 - x3) / DY;
        dN2x = (N2x - N3x) / DY;
        dN2y = (N2y - N3y) / DY;
        dN2z = (N2z - N3z) / DY;
        dzf  = (z2 - z3) / DY;
      }
      else
        dN1x = dN1y = dN1z = dN2x = dN2y = dN2z = dx1 = dx2 = dzi = dzf = 0;
    }
    else{
      Nix = Nfx = N1x;
      Niy = Nfy = N1y;
      Niz = Nfz = N1z;
      xi = xf = x1;
      zi = zf = z1;
      y = floor(y1);

      DY   = y2 - y1;
      dx1  = (x2 - x1) / DY;
      dx2  = (x3 - x1) / (y3 - y1);
      dN1x = (N2x - N1x) / DY;
      dN1y = (N2y - N1y) / DY;
      dN1z = (N2z - N1z) / DY;
      dzi  = (z2 - z1)/DY;

      DY   = y3 - y1;
      dx2  = (x3 - x1) / DY;
      dN2x = (N3x - N1x) / DY;
      dN2y = (N3y - N1y) / DY;
      dN2z = (N3z - N1z) / DY;
      dzf  = (z3 - z1) / DY;
    }
  }

  while ((y < y2) || (y < y3)){
    if (xf != xi){
      DY  = xf - xi;
      dNx = (Nfx - Nix) / DY;
      dNy = (Nfy - Niy) / DY;
      dNz = (Nfz - Niz) / DY;
      dz  = (zf - zi) / DY;
    }
    else{
      dNx = 0;
      dNy = 0;
      dNz = 0;
      dz = 0;
    }

    Nx = Nix;
    Ny = Niy;
    Nz = Niz;
    z = zi;

    if (xi <= xf)
      for (x=(int)xi; x<= (int)xf; x++){
        cosAlpha = (Nx * luxX + Ny * luxY + Nz * luxZ) / (sqrt(Nx * Nx + Ny * Ny + Nz * Nz) * NormeLux);
        C = cosAlpha >= 0 ? cosAlpha : 0;
        C *= C;
        // C *= C;
        // C *= C;
        // C *= C;
        C += AMBIANT_LIGHT;
        Color = (int)(C * nb_colors);
        if (Color < 0) Color = 0;
        if (Color >= nb_colors) Color = nb_colors - 1;
        Point(x, y, Palette[Color], z);
        Nx += dNx;
        Ny += dNy;
        Nz += dNz;
        z  += dz;
      } else
      for (x = (int)xi; x>= (int)xf; x--){
        cosAlpha = (Nx * luxX + Ny * luxY + Nz * luxZ) / (sqrt(Nx * Nx + Ny * Ny + Nz * Nz) * NormeLux);
        C = cosAlpha >= 0 ? cosAlpha : 0;
        C *= C;
        // C *= C;
        // C *= C;
        // C *= C;
        C += AMBIANT_LIGHT;
        Color = (int)(C * nb_colors);
        if (Color < 0) Color = 0;
        if (Color >= nb_colors) Color = nb_colors - 1;
        Point(x, y, Palette[Color], z);
        Nx -= dNx;
        Ny -= dNy;
        Nz -= dNz;
        z  -= dz;
      }

    y++;

    if(!Passe && (y > y2)){
      if (y3 != y2){
        DY   = y3 - y2;
        dx1  = (x3 - x2) / DY;
        dN1x = (N3x - N2x) / DY;
        dN1y = (N3y - N2y) / DY;
        dN1z = (N3z - N2z) / DY;
        dzi  = (z3 - z2) / DY;
      }else
        dx1 = dzi = dN1x = dN1y = dN1z = 0;
      Passe = 1;
    }
    if(!Passe && (y > y3)){
      if (y2 != y3){
        DY   = y2 - y3;
        dx2  = (x2 - x3) / DY;
        dN2x = (N2x - N3x) / DY;
        dN2y = (N2y - N3y) / DY;
        dN2z = (N2z - N3z) / DY;
        dzf  = (z2 - z3) / DY;
      } else
        dx2 = dzf = dN2x = dN2y = dN2z = 0;
      Passe = 1;
    }
    xi  += dx1;
    xf  += dx2;
    Nix += dN1x;
    Niy += dN1y;
    Niz += dN1z;
    Nfx += dN2x;
    Nfy += dN2y;
    Nfz += dN2z;
    zi  += dzi;
    zf  += dzf;
  }
}


void dessine_scene(struct objet O, int mode){
  int i,j;
  long offset;
  point P[3];
  float  C[3], z[3], N[9], Nx, Ny, Nz, cosAlpha;
  float NormeLux = sqrt(luxX * luxX + luxY * luxY + luxZ * luxZ);
  int OX = areaL / 2;          // position of the center
  int OY = areaH / 2;
  int MinXO = areaL/2, MinYO = areaH/2, MaxXO = areaL/2, MaxYO = areaH/2;  // center

  // RAZ

  for (j=MinY; j<MaxY; j++){
    offset = j * areaL + MinX;
    for (i=MinX; i<MaxX; i++){
      Zb[offset] = -100000;
      Image[offset] = 0;
      offset++;
    }
  }

  MinX = MaxX = OX;
  MinY = MaxY = OY;

  //
  for (i=0; i<O.nb_faces; i++){
    Nz = O.Fnormales[i][Z];
    if ((Nz <= 0) && (mode != WIRE)) continue;        // SI LA FACE EST ARRIERE

    z[0] = O.sommets[O.faces[i][0]][Z];
    if (z[0] * DObj > DObs) continue;                 // SI LA FACE EST DERRIERE L'OBSERVATEUR
    z[1] = O.sommets[O.faces[i][1]][Z];
    if (z[1] * DObj > DObs) continue;
    z[2] = O.sommets[O.faces[i][2]][Z];
    if (z[2] * DObj > DObs) continue;


    P[0].x = (O.sommets[O.faces[i][0]][X] * DObs / (DObs - z[0])) * DObj + OX;      // PERSPECTIVE (Thales)
    P[0].y = (O.sommets[O.faces[i][0]][Y] * DObs / (DObs - z[0])) * DObj + OY;      // + centrage et zoom

    P[1].x = (O.sommets[O.faces[i][1]][X] * DObs / (DObs - z[1])) * DObj + OX;
    P[1].y = (O.sommets[O.faces[i][1]][Y] * DObs / (DObs - z[1])) * DObj + OY;

    P[2].x = (O.sommets[O.faces[i][2]][X] * DObs / (DObs - z[2])) * DObj + OX;
    P[2].y = (O.sommets[O.faces[i][2]][Y] * DObs / (DObs - z[2])) * DObj + OY;

    if( ( (P[0].x < 0) && (P[1].x < 0) && (P[2].x < 0) )             // si la face n'est pas dans l'aire de visibilite
     || ( (P[0].y < 0) && (P[1].y < 0) && (P[2].y < 0) )
     || ( (P[0].x > areaL) && (P[1].x > areaL) && (P[2].x > areaL) )
     || ( (P[0].y > areaH) && (P[1].y > areaH) && (P[2].y > areaH) ) ) continue;

    if ((mode == LINE) || (mode == SOLID)){
      Nx = O.Fnormales[i][X];
      Ny = O.Fnormales[i][Y];
      cosAlpha = (Nx * luxX + Ny * luxY + Nz * luxZ)/(sqrt(Nx * Nx + Ny * Ny + Nz * Nz) * NormeLux);
      C[0] = cosAlpha >= 0 ? cosAlpha : 0;
    }

    if (mode == GOURAUD){
      Nx = O.Pnormales[O.faces[i][0]][X];
      Ny = O.Pnormales[O.faces[i][0]][Y];
      Nz = O.Pnormales[O.faces[i][0]][Z];
      cosAlpha = (Nx * luxX + Ny * luxY + Nz * luxZ) / (sqrt(Nx * Nx + Ny * Ny + Nz * Nz) * NormeLux);
      C[0] = cosAlpha >= 0 ? cosAlpha : 0;


      Nx = O.Pnormales[O.faces[i][1]][X];
      Ny = O.Pnormales[O.faces[i][1]][Y];
      Nz = O.Pnormales[O.faces[i][1]][Z];
      cosAlpha = (Nx * luxX + Ny * luxY + Nz * luxZ) / (sqrt(Nx * Nx + Ny * Ny + Nz * Nz)* NormeLux);
      C[1] = cosAlpha>= 0?cosAlpha:0;

      Nx = O.Pnormales[O.faces[i][2]][X];
      Ny = O.Pnormales[O.faces[i][2]][Y];
      Nz = O.Pnormales[O.faces[i][2]][Z];
      cosAlpha = (Nx*luxX+Ny*luxY+Nz*luxZ)/(sqrt(Nx*Nx+Ny*Ny+Nz*Nz)*NormeLux);
      C[2] = cosAlpha>= 0 ? cosAlpha : 0;
    }

    if (mode == PHONG){
      N[0] = O.Pnormales[O.faces[i][0]][X];
      N[1] = O.Pnormales[O.faces[i][0]][Y];
      N[2] = O.Pnormales[O.faces[i][0]][Z];

      N[3] = O.Pnormales[O.faces[i][1]][X];
      N[4] = O.Pnormales[O.faces[i][1]][Y];
      N[5] = O.Pnormales[O.faces[i][1]][Z];

      N[6] = O.Pnormales[O.faces[i][2]][X];
      N[7] = O.Pnormales[O.faces[i][2]][Y];
      N[8] = O.Pnormales[O.faces[i][2]][Z];
    }

    switch (mode){
      case WIRE    : lines(P, 1, z)    ;break;
      case LINE    : lines(P, C[0], z) ;break;
      case SOLID    : flat(P, C[0], z)  ;break;
      case GOURAUD : gouraud(P, C, z)  ;break;
      case PHONG   : phong(P, N, z)    ;break;
    }
  }

  if (MinX < MinXO) MinXO = MinX;
  if (MaxX > MaxXO) MaxXO = MaxX;
  if (MinY < MinYO) MinYO = MinY;
  if (MaxY > MaxYO) MaxYO = MaxY;

  MinXO = MinX;
  MaxXO = MaxX;
  MinYO = MinY;
  MaxYO = MaxY;
}


void rotat(struct objet *O, float AX, float AY, float AZ){
  float rotX, rotY, rotZ, tempX, tempY;
  float cosAX = cos(AX);
  float sinAX = sin(AX);
  float cosAY = cos(AY);
  float sinAY = sin(AY);
  float cosAZ = cos(AZ);
  float sinAZ = sin(AZ);
  int i;

  for (i=0; i<O->nb_faces; i++){
    rotX = O->Fnormales[i][X];
    rotY = O->Fnormales[i][Y];
    rotZ = O->Fnormales[i][Z];

    tempY = rotY;
    rotY = (tempY * cosAX) - (rotZ * sinAX);
    rotZ = (tempY * sinAX) + (rotZ * cosAX);
    tempX = rotX;
    rotX = (tempX * cosAY) - (rotZ * sinAY);
    rotZ = (tempX * sinAY) + (rotZ * cosAY);
    tempX = rotX;
    rotX = (tempX * cosAZ) - (rotY * sinAZ);
    rotY = (tempX * sinAZ) + (rotY * cosAZ);

    O->Fnormales[i][X] = rotX;
    O->Fnormales[i][Y] = rotY;
    O->Fnormales[i][Z] = rotZ;
  }

  for (i = 0;i<O->nb_sommets;i++){
    rotX = O->sommets[i][X];
    rotY = O->sommets[i][Y];
    rotZ = O->sommets[i][Z];

    tempY = rotY;
    rotY = (tempY * cosAX) - (rotZ * sinAX);
    rotZ = (tempY * sinAX) + (rotZ * cosAX);
    tempX = rotX;
    rotX = (tempX * cosAY) - (rotZ * sinAY);
    rotZ = (tempX * sinAY) + (rotZ * cosAY);
    tempX = rotX;
    rotX = (tempX*cosAZ)-(rotY*sinAZ);
    rotY = (tempX*sinAZ)+(rotY*cosAZ);

    O->sommets[i][X] = rotX;
    O->sommets[i][Y] = rotY;
    O->sommets[i][Z] = rotZ;

    rotX = O->Pnormales[i][X];
    rotY = O->Pnormales[i][Y];
    rotZ = O->Pnormales[i][Z];

    tempY = rotY;
    rotY = (tempY * cosAX) - (rotZ * sinAX);
    rotZ = (tempY * sinAX) + (rotZ * cosAX);
    tempX = rotX;
    rotX = (tempX * cosAY) - (rotZ * sinAY);
    rotZ = (tempX * sinAY) + (rotZ * cosAY);
    tempX = rotX;
    rotX = (tempX * cosAZ) - (rotY * sinAZ);
    rotY = (tempX * sinAZ) + (rotY * cosAZ);

    O->Pnormales[i][X] = rotX;
    O->Pnormales[i][Y] = rotY;
    O->Pnormales[i][Z] = rotZ;
  }
}


void translate(struct objet *O, float dX, float dY, float dZ) {
  for (int i=0; i<O->nb_sommets; i++){
    O->sommets[i][X] += dX;
    O->sommets[i][Y] += dY;
    O->sommets[i][Z] += dZ;
  }
}


int lit_fichier(struct objet *O, char *nom_fichier) {
  FILE *F;

  char ligne[128] = {0};
  char type = 'z';
  float x, y, z, OCx = 0, OCy = 0, OCz = 0;
  int a, b, c, i, n;

  float ux, uy, uz, vx, vy, vz;
  float x1, y1, z1, x2, y2, z2, x3, y3, z3;

  if ((F = fopen(nom_fichier, "r")) == NULL) return -1;

  fgets(ligne,128, F);
  sscanf(ligne, "%c%d", &type, &i);
  if (type != 'P') return -2;
  O->nb_sommets = i;

  fgets(ligne, 128, F);
  sscanf(ligne, "%c%d", &type, &i);
  if (type != 'F') return -3;
  O->nb_faces = i;

  free(O->sommets);
  free(O->faces);
  free(O->Pnormales);
  free(O->Fnormales);

  if(( (O->sommets = (sommet *)malloc(O->nb_sommets    * sizeof(sommet)))  ==  NULL)
   ||( (O->Pnormales = (normale *)malloc(O->nb_sommets * sizeof(normale))) ==  NULL)
   ||( (O->faces = (face *)malloc(O->nb_faces          * sizeof(face)))    ==  NULL)
   ||( (O->Fnormales = (normale *)malloc(O->nb_faces   * sizeof(normale))) ==  NULL))
    return -4;


  for (i=0; i<O->nb_sommets; i++){
    fgets(ligne, 128, F);
    if (sscanf(ligne,"%c%f%f%f",&type, &x, &y, &z)!= 4) return -5;
    if (type != 'p') return -6;

    O->sommets[i][X] = x;
    O->sommets[i][Y] = y;
    O->sommets[i][Z] = z;

    OCx += x/O->nb_sommets;
    OCy += y/O->nb_sommets;
    OCz += z/O->nb_sommets;

    O->sommets[i][3] = 0;

    O->Pnormales[i][X] = 0;
    O->Pnormales[i][Y] = 0;
    O->Pnormales[i][Z] = 0;
  }

  for (i=0; i<O->nb_sommets; i++){
    O->sommets[i][X] -= OCx;
    O->sommets[i][Y] -= OCy;
    O->sommets[i][Z] -= OCz;
  }

  for (i=0; i<O->nb_faces; i++){
    fgets(ligne, 128, F);
    if (sscanf(ligne, "%c%d%d%d", &type, &a, &b, &c) != 4) return -7;
    if (type != 'f') return -8;

    O->faces[i][0] = --a;
    O->faces[i][1] = --b;
    O->faces[i][2] = --c;

    O->sommets[a][3]++;     // on incremente le nombre de faces auquelles le point appartient
    O->sommets[b][3]++;
    O->sommets[c][3]++;

    x1 = O->sommets[a][X];
    y1 = O->sommets[a][Y];
    z1 = O->sommets[a][Z];

    x2 = O->sommets[b][X];
    y2 = O->sommets[b][Y];
    z2 = O->sommets[b][Z];

    x3 = O->sommets[c][X];
    y3 = O->sommets[c][Y];
    z3 = O->sommets[c][Z];

    ux = x2 - x1;             // Calcul de la normale a la face
    uy = y2 - y1;
    uz = z2 - z1;
    vx = x3 - x1;
    vy = y3 - y1;
    vz = z3 - z1;

    O->Fnormales[i][X] = (uy * vz) - (vy * uz);
    O->Fnormales[i][Y] = (vx * uz) - (ux * vz);
    O->Fnormales[i][Z] = (ux * vy) - (vx * uy);
  }

  for (i=0; i<O->nb_faces; i++){
    a = O->faces[i][0];
    b = O->faces[i][1];
    c = O->faces[i][2];

    n = O->sommets[a][3];       // nombre de faces
    O->Pnormales[a][X] += (O->Fnormales[i][X] / n);
    O->Pnormales[a][Y] += (O->Fnormales[i][Y] / n);
    O->Pnormales[a][Z] += (O->Fnormales[i][Z] / n);

    n = O->sommets[b][3];
    O->Pnormales[b][X] += (O->Fnormales[i][X] / n);
    O->Pnormales[b][Y] += (O->Fnormales[i][Y] / n);
    O->Pnormales[b][Z] += (O->Fnormales[i][Z] / n);

    n = O->sommets[c][3];
    O->Pnormales[c][X] += (O->Fnormales[i][X] / n);
    O->Pnormales[c][Y] += (O->Fnormales[i][Y] / n);
    O->Pnormales[c][Z] += (O->Fnormales[i][Z] / n);
  }
  return 1;
}


int main (int argc, char** argv){

  SDL_Init(SDL_INIT_VIDEO);
  SDL_Window *window = SDL_CreateWindow("3d renderer", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, areaL, areaH, SDL_WINDOW_OPENGL);
  SDL_Surface *surface = SDL_GetWindowSurface(window);

  Image = (Uint32*) surface->pixels;  // we are ploting directly on the surface of the window

  SDL_Event event;
  SDL_bool running = SDL_TRUE;
  const int frameDelay = 1000 / FPS;
  Uint32 frameStart, frameTime;

  Zb = (float*)malloc(areaL * areaH * sizeof(float));  // Z-buffer
  for (int i=0; i<256; i++) Palette[i] = i * 65793;   // 256 shades of grey
  float ax = 0, ay = 0, az = 0;                      // angles de rotation
  struct objet OBJ = {0};                           // the object we will render
  int mode = PHONG;                                // rendering mode

  while(running){

    frameStart = SDL_GetTicks();
    rotat(&OBJ, ax, ay, az);  // rotate the object

    dessine_scene(OBJ, mode);

    if (SDL_PollEvent(&event)){
      if (event.type == SDL_QUIT) running = SDL_FALSE;

      if (event.type == SDL_DROPFILE){
        char *filename = event.drop.file;
        lit_fichier(&OBJ, filename);
        SDL_free(filename);
        ax = ay = az = 0;  // stop the rotation
        DObj = 1;
      }

      if (event.type == SDL_KEYDOWN){
        switch(event.key.keysym.sym){
          case SDLK_q: translate(&OBJ, 1, 0, 0); break;  // translation
          case SDLK_w: translate(&OBJ, 0, 1, 0); break;
          case SDLK_e: translate(&OBJ, 0, 0, 1); break;
          case SDLK_a: translate(&OBJ,-1, 0, 0); break;
          case SDLK_s: translate(&OBJ, 0,-1, 0); break;
          case SDLK_d: translate(&OBJ, 0, 0,-1); break;

          case SDLK_u: ax += .01; break;  // rotation angles
          case SDLK_i: ay += .01; break;
          case SDLK_o: az += .01; break;
          case SDLK_j: ax -= .01; break;
          case SDLK_k: ay -= .01; break;
          case SDLK_l: az -= .01; break;
          case SDLK_SPACE: ax = ay = az = 0; break;

          case SDLK_r: DObj++; break;
          case SDLK_f: if (DObj > 1) DObj--; break;

          case SDLK_1: mode = WIRE;    break;  // renderer mode
          case SDLK_2: mode = LINE;    break;
          case SDLK_3: mode = SOLID;    break;
          case SDLK_4: mode = GOURAUD; break;
          case SDLK_5: mode = PHONG;   break;

          case SDLK_ESCAPE: running = SDL_FALSE; break;  // exit

          default :
          SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_INFORMATION,
          "About", "3D Renderer\n\
          \nDrop a .obj file in the window and use the following keys :\n\
          \n - q, w, e\ttranslate + on x, y, z axis\
          \n - a, s, d\ttranslate - on x, y, z axis\
          \n - r, f   \tincrease / decrease viewing distance\
          \n - u, i, o\tincrease the rotation angles on x, y, z axis\
          \n - j, k, l\tdecrease the rotation angles on x, y, z axis\
          \n - 1,2,3,4,5 set the rendering mode\
          \n - <space> reset the rotaion angles\
          \n\n - any other key : this help\
          \n\
          \nMore information at github.com/ArthurFerreira2", NULL);
          break;
        }
      }
    }

    frameTime = SDL_GetTicks() - frameStart;  // time since last rendering
    if (frameDelay > frameTime) SDL_Delay(frameDelay - frameTime); // wait if we have spare time

    SDL_UpdateWindowSurface(window);  // render frame
  }

  SDL_FreeSurface(surface);
  SDL_DestroyWindow(window);
  SDL_Quit();
  exit(EXIT_SUCCESS);
}
