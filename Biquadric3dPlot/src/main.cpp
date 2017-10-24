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

  std::ofstream q, q_bar, q_bar2, background, non_squares;
  q.open("q.txt"); q_bar.open("q_bar.txt"); q_bar2.open("q_bar2.txt"); background.open("background.txt"); non_squares.open("non_squares.txt");

  std::vector <GaloisScalar> squares;
  for(GaloisScalar candidate{1}; candidate < squaresCardinality + 1; ++candidate) {
    squares.push_back(candidate*candidate);
  }
  std::vector<GaloisScalar> nonSquares;
  non_squares << "Numbers without a square root for p = " << primeNumber << ":" << std::endl;
  for(int candidate{1}; candidate < primeNumber; ++candidate) {
    if(std::find(squares.begin(), squares.end(), candidate) != squares.end()) {
    } else {
      nonSquares.push_back(candidate);
      non_squares << candidate << " ";
    }
  }

  GaloisScalar m11{std::atoi(argv[2])};
  GaloisScalar m12{std::atoi(argv[3])};
  GaloisScalar m13{std::atoi(argv[4])};
  GaloisScalar m22{std::atoi(argv[5])};
  GaloisScalar m23{std::atoi(argv[6])};
  GaloisScalar m33{std::atoi(argv[7])};

  GaloisScalar zero{0}; GaloisScalar one{1};


  for(int x{squaresCardinality}; x>-squaresCardinality-1; --x) {
    for(int y{-squaresCardinality}; y<squaresCardinality+1; ++y) {
      for(int z{-squaresCardinality}; z<squaresCardinality+1; ++z) {
        if(((x*x)*m11 + 2*((x*y)*m12) + 2*((x*z)*m13) + (y*y)*m22 + 2*((y*z)*m23) + (z*z)*m33 + one) == zero) {
          q << x << "\t" << y << "\t" << z << std::endl;
        }
        if(((x*x)*m11 + 2*((x*y)*m12) + 2*((x*z)*m13) + (y*y)*m22 + 2*((y*z)*m23) + (z*z)*m33 + nonSquares[std::atoi(argv[8])]) == zero) {
          q_bar << x << "\t" << y << "\t" << z << std::endl;
        }
        if(((x*x)*m11 + 2*((x*y)*m12) + 2*((x*z)*m13) + (y*y)*m22 + 2*((y*z)*m23) + (z*z)*m33 + nonSquares[std::atoi(argv[9])]) == zero) {
          q_bar2 << x << "\t" << y << "\t" << z << std::endl;
        }
        background << x << "\t" << y << "\t" << z << std::endl;
        }
      }
    }

  q.close(); q_bar.close(); q_bar2.close(); background.close(); non_squares.close();

  return 0;
}
