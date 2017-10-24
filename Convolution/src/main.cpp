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

  std::ofstream q, q_bar, /* q_bar2, */  background, non_squares, countMapFinal, hist;
  q.open("q.txt"); q_bar.open("q_bar.txt"); /* q_bar2.open("q_bar2.txt"); */ background.open("background.txt"); non_squares.open("non_squares.txt"); countMapFinal.open("count_map.txt");
  std::string histFilename("hist_p");
  histFilename.append(std::to_string(primeNumber));
  histFilename.append("_F");
  histFilename.append(argv[5]);
  histFilename.append("_mink.txt");
  hist.open(histFilename);

  std::vector<unsigned int> histogram(2*primeNumber + 2, 0);

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
  GaloisScalar m22{std::atoi(argv[4])};

  GaloisScalar zero{0}; GaloisScalar one{1};

  std::vector<GaloisOneZeroTensor > biquadric(2*primeNumber);

  GaloisScalar FSquare{std::atoi(argv[5])*std::atoi(argv[5])};
  GaloisScalar minusOne{nonSquares[squaresCardinality - 1]};

  int counter{0};
  for(int x{1}; x<primeNumber; ++x) {
    GaloisScalar xGal{x};
    for(int y{1}; y<primeNumber; ++y) {
      GaloisScalar yGal{y};
      if((m11*(xGal*xGal) + 2*(m12*(xGal*yGal)) + m22*(yGal*yGal) + FSquare == 0) || (m11*(xGal*xGal) + 2*(m12*(xGal*yGal)) + m22*(yGal*yGal) + minusOne*FSquare == 0)) {
        biquadric[counter][0] = x;
        biquadric[counter][1] = y;
        biquadric[counter][2] = 1;
        ++counter;
        /*
        biquadric[x][0] = x;
        biquadric[x][1] = (-2)*x;
        biquadric[x][1] = biquadric[x][1].get_inverse();
        biquadric[x][2] = 1;

        biquadric[primeNumber + x][0] = x;
        biquadric[primeNumber + x][1] = 2*x;
        biquadric[primeNumber + x][1] = biquadric[x][1].get_inverse();
        biquadric[primeNumber + x][2] = 1;
        */
      }
    }
    if((m11*(xGal*xGal) + FSquare == 0) || (m11*(xGal*xGal) + minusOne*FSquare == 0)) {
      biquadric[counter][0] = x;
      biquadric[counter][1] = 1;
      biquadric[counter][2] = 0;
      ++counter;
    }
  }

  if(m11 == 0) {
    biquadric[counter][0] = 1;
    biquadric[counter][1] = 0;
    biquadric[counter][2] = 0;
  }

  /*
  biquadric[primeNumber][0] = 1;
  biquadric[primeNumber][1] = 0;
  biquadric[primeNumber][2] = 0;
  */

  int ceiledSquareRootOfPrimeNumber{static_cast<int>(ceil(sqrt(primeNumber)))};
  std::vector<std::vector<double> > countMap(2*ceiledSquareRootOfPrimeNumber + 1, std::vector<double>(2*ceiledSquareRootOfPrimeNumber + 1));

  for(GaloisOneZeroTensor point : biquadric) {
    if (point[0] < squaresCardinality && point[1] < squaresCardinality) {
      ++countMap[static_cast<int>(floor(sqrt(static_cast<double>(point[0]*point[0])))) + ceiledSquareRootOfPrimeNumber][static_cast<int>(floor(sqrt(static_cast<double>(point[1]*point[1])))) + ceiledSquareRootOfPrimeNumber];
    }
    if (point[0] < squaresCardinality && point[1] >= squaresCardinality) {
      ++countMap[static_cast<int>(floor(sqrt(static_cast<double>(point[0]*point[0])))) + ceiledSquareRootOfPrimeNumber][-static_cast<int>(floor(sqrt(static_cast<double>(point[1]*point[1])))) + ceiledSquareRootOfPrimeNumber];
    }
    if (point[0] >= squaresCardinality && point[1] < squaresCardinality) {
      ++countMap[-static_cast<int>(floor(sqrt(static_cast<double>(point[0]*point[0])))) + ceiledSquareRootOfPrimeNumber][static_cast<int>(floor(sqrt(static_cast<double>(point[1]*point[1])))) + ceiledSquareRootOfPrimeNumber];
    }
    if (point[0] >= squaresCardinality && point[1] >= squaresCardinality) {
      ++countMap[-static_cast<int>(floor(sqrt(static_cast<double>(point[0]*point[0])))) + ceiledSquareRootOfPrimeNumber][-static_cast<int>(floor(sqrt(static_cast<double>(point[1]*point[1])))) + ceiledSquareRootOfPrimeNumber];
    }
  }

  unsigned int cellCounter{0};
  for (auto countColumn : countMap) {
    for (auto count : countColumn) {
      countMapFinal << count << " ";
      if ( count != 0) {
        ++cellCounter;
      }
      ++histogram[count];
    }
    countMapFinal << std::endl;
  }
  std::cout << "Number of Cells: " << cellCounter << std::endl;
  int bin_count{0};
  for (int i{1}; i< 2*primeNumber + 2; ++i) {
    if (histogram[i] != 0) {
      hist << i << "\t" << histogram[i] << std::endl;
      ++bin_count;
    }
  }
  std::cout << "Histogram value for '0' = " << histogram[0] << std::endl;

  int check{0};
  double hist_mean{0};
  for (unsigned int i{1}; i < 2*primeNumber + 2; ++i) {
    check += i*histogram[i];
    hist_mean += histogram[i];
  }
  std::cout << "check = " << check << " and mean value of the histogram = " << hist_mean/bin_count << std::endl;

  q.close(); q_bar.close(); /* q_bar2.close(); */ background.close(); non_squares.close(); countMapFinal.close(); hist.close();

  return 0;
}
