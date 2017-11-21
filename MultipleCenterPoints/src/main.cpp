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

  int numberOfPoints{primeNumber*primeNumber + primeNumber + 1};

  int squaresCardinality{static_cast<int>((primeNumber - 1)*0.5)};

  std::vector<GaloisOneZeroTensor> points{set_of_points(primeNumber,dimension+1,numberOfPoints)};
  std::vector<GaloisZeroOneTensor> hyperplanes{set_of_hyperplanes(points)};

  std::ofstream q, q_bar, q_bar2, background, non_squares;
  non_squares.open("non_squares.txt");

  std::vector<GaloisScalar> squares;
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

  GaloisScalar zero{0}; GaloisScalar one{1}; GaloisScalar minusOne{-1};

  std::vector<GaloisZeroTwoTensor> allQuadrics;
  GaloisZeroTwoTensor quadric;
  for (int a{0}; a<primeNumber-1; ++a) {
    for (int b{0}; b<primeNumber-1; ++b) {
      for (int c{0}; c<primeNumber-1; ++c) {
        for (int d{0}; d<primeNumber-1; ++d) {
          for (int e{0}; e<primeNumber-1; ++e) {
            quadric[0][0] = a; quadric[1][0] = b; quadric[2][0] = c;
            quadric[0][1] = b; quadric[1][1] = d; quadric[2][1] = e;
            quadric[0][2] = c; quadric[1][2] = e; quadric[2][2] = 1;
            GaloisScalar determinant{a*d + b*e*c + c*b*e - c*c*d - e*e*a - b*b};
            if (determinant != 0) {
              allQuadrics.push_back(quadric);
            }
          }

          quadric[0][0] = a; quadric[1][0] = b; quadric[2][0] = c;
          quadric[0][1] = b; quadric[1][1] = d; quadric[2][1] = 1;
          quadric[0][2] = c; quadric[1][2] = 1; quadric[2][2] = 0;
          GaloisScalar determinant{b*c + c*b - c*c*d - a};
          if (determinant != 0) {
            allQuadrics.push_back(quadric);
          }

          if (c == 0) {
            quadric[0][0] = a; quadric[1][0] = b; quadric[2][0] = 1;
            quadric[0][1] = b; quadric[1][1] = d; quadric[2][1] = 0;
            quadric[0][2] = 1; quadric[1][2] = 0; quadric[2][2] = 0;
            if (d != 0) {
              allQuadrics.push_back(quadric);
            }
          }
        }
      }
    }
  }

  GaloisZeroTwoTensor emptyQuadric;

  std::vector<GaloisZeroTwoTensor> pair(2, emptyQuadric);

  std::vector<std::vector<GaloisZeroTwoTensor> > biquadrics;
  std::vector<std::vector<GaloisOneZeroTensor> > centerPointSets;
  std::vector<GaloisOneZeroTensor> emptyPointSet;
  size_t n{0};
  for (size_t k{0}; k < allQuadrics.size(); ++k) {
    for (size_t l{k+1}; l < allQuadrics.size(); ++l) {
      pair[0] = allQuadrics[k];
      pair[1] = allQuadrics[l];
      bool areThereNoCenterPoints{true};
      for (GaloisOneZeroTensor centerPointCandidate : points) {
        if (is_center(centerPointCandidate, pair, points, hyperplanes)) {
          if (areThereNoCenterPoints) {
            centerPointSets.push_back(emptyPointSet);
            areThereNoCenterPoints = false;
          }
          centerPointSets[n].push_back(centerPointCandidate);
        }
        if (areThereNoCenterPoints == false) {
          biquadrics.push_back(pair);
          ++n;
          areThereNoCenterPoints = true;
        }
      }
    }
  }

  size_t degeneracyCount{0};
  size_t biquadricCount{0};
  for (std::vector<GaloisZeroTwoTensor> biquadric : biquadrics) {
    std::cout << "Biquadric:\n---\n" << biquadric[0] << "\n" << biquadric[1] << "---" << std::endl;
    for (GaloisOneZeroTensor point : centerPointSets[biquadricCount]) {
      std::cout << point << std::endl;
    }
    if (centerPointSets[biquadricCount].size() > 1) {
      ++degeneracyCount;
    }
    ++biquadricCount;
  }

  std::cout << "Biquadrics with more than one center point: " << degeneracyCount << std::endl;
  non_squares.close();

  return 0;
}
