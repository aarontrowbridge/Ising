// spin_bodies.h

#ifndef SPIN_BODIES_H
#define SPIN_BODIES_H

#include <array>
#include <math.h>
#include <stdlib.h>
#include <time.h>

namespace SpinBodies
{
  int k_index (int i, int j, int n)
  {
    return n * i + j;
  }

  class SpinBody
  {
  public:
    int s;
    int E;
    int i, j, k;

    SpinBody (int, int, int);
  };

  SpinBody::SpinBody (int i_index, int j_index, int n)
  {
    srand (time(NULL));
    s = rand() % 2;
    E = 0;
    i = i_index;
    j = j_index;
    k = k_index(i_index, j_index, n);
  }

  class SpinLattice
  {
  public:
    static int n;
    static int N;
    SpinBody bs[n][n];
    double T;
    int steps;
    int flips;

    SpinLattice (int, double);
  };

  SpinLattice::SpinLattice (int d, double temp)
  {
    n = d;
    N = d*d;
    for (int i=0, j=0; i < d, j < d; i++, j++)
      {
        bs[i][j] = SpinBody(i, j, d);
      }
    T = temp;
    steps = 0;
    flips = 0;
  }
}


#endif
