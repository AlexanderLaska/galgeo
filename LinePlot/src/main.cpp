#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <random>
#include <algorithm>
#include <cstdio>
#include <ctime>

#include <GaloisScalar.hpp>
#include <GaloisOneZeroTensor.hpp>
#include <GaloisZeroOneTensor.hpp>
#include <GaloisOneOneTensor.hpp>
#include <GaloisZeroTwoTensor.hpp>
#include <GaloisTwoZeroTensor.hpp>

#include <FiniteProjectiveGeometry.hpp>

using namespace galois;


int main(int argc, char *argv[]) {

  // Use the first command line argument as the prime number
  int primeNumber{std::atoi(argv[1])};
  set_prime(primeNumber);

  int dimension{2};
  set_dimension(dimension+1);

  int squaresCardinality{static_cast<int>((primeNumber - 1)*0.5)};

  GaloisScalar l1{std::atoi(argv[2])};
  GaloisScalar l2{std::atoi(argv[3])};
  GaloisScalar l3{std::atoi(argv[4])};

  GaloisScalar zero{0}; GaloisScalar one{1};

  std::ofstream line, background;
  line.open("line.txt"); background.open("background.txt");

  for(int x{squaresCardinality}; x>-squaresCardinality-1; --x) {
    for(int y{-squaresCardinality}; y<squaresCardinality+1; ++y) {
      if(x*l1 + y*l2 + l3 == zero) {
        line << x << "\t" << y << std::endl;
      }
      background << x << "\t" << y << std::endl;
    }
  }
  
  line.close(); background.close();

  return 0;
}
