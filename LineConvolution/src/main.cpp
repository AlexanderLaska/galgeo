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

using namespace galgeo;


int main(int argc, char *argv[]) {

  // Use the first command line argument as the prime number
  int primeNumber{std::atoi(argv[1])};
  set_prime(primeNumber);

  int dimension{2};
  set_dimension(dimension+1);

  int squaresCardinality{static_cast<int>((primeNumber - 1)*0.5)};

  std::ofstream q, q_bar, /* q_bar2, */  background, non_squares, countMapFinal;
  q.open("points.txt"); background.open("background.txt"); non_squares.open("non_squares.txt"); countMapFinal.open("count_map.txt");

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

  GaloisScalar zero{0}; GaloisScalar one{1};

  std::vector<GaloisOneZeroTensor > line(2*primeNumber+1);
  line[0][0] = 0;
  line[0][1] = 0;
  line[0][2] = 1;
  for(int x{1}; x<primeNumber; ++x) {
    line[x][0] = 1*x+0;
    line[x][1] = (-1)*x+0;
    line[x][2] = 1;

    line[2*x][0] = 1*x+0;
    line[2*x][1] = 1*x+0;
    line[2*x][2] = 1;
  }
  line[primeNumber][0] = 1;
  line[primeNumber][1] = 1;
  line[primeNumber][2] = 0;

  line[2*primeNumber][0] = -1;
  line[2*primeNumber][1] = 1;
  line[2*primeNumber][2] = 0;

  int ceiledSquareRootOfPrimeNumber{static_cast<int>(ceil(sqrt(primeNumber)))};
  std::vector<std::vector<double> > countMap(ceiledSquareRootOfPrimeNumber + primeNumber, std::vector<double>(ceiledSquareRootOfPrimeNumber + primeNumber));

  for(GaloisOneZeroTensor point : line) {
    ++countMap[static_cast<int>(floor(sqrt(static_cast<double>(point[0]*point[0]))))][static_cast<int>(floor(sqrt(static_cast<double>(point[1]*point[1]))))];
  }

  for (auto countColumn : countMap) {
    for (auto count : countColumn) {
      countMapFinal << count << " ";
    }
    countMapFinal << std::endl;
  }
  q.close(); q_bar.close(); /* q_bar2.close(); */ background.close(); non_squares.close(); countMapFinal.close();

  return 0;
}
