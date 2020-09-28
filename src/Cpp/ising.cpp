// ising.cpp

#include "spin_bodies.h"
#include <iostream>

using namespace std;
using namespace SpinBodies;

int main ()
{
  SpinBody b;
  b.init (1, 0, 2);
  cout << "s = " << b.s << "\nk = " << b.k << "\n";
  return 0;
}
