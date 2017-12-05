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

  std::ofstream data, non_squares, hist;

  std::string dataFilename("data_p");
  dataFilename.append(std::to_string(primeNumber));
  dataFilename.append(".txt");
  data.open(dataFilename);
  
  std::string nonSquaresFilename("nonSquares_p");
  nonSquaresFilename.append(std::to_string(primeNumber));
  nonSquaresFilename.append(".txt");
  non_squares.open(nonSquaresFilename);

  std::string histFilename("hist_p");
  histFilename.append(std::to_string(primeNumber));
  histFilename.append(".txt");
  hist.open(histFilename);

  std::vector<unsigned int> histogram(primeNumber, 0);

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

  data << "#x\ts(x)\n";
  int intMod;
  for (int n{0}; n < primeNumber; ++n) {
    GaloisScalar nSquared{n*n};
    //if (n < squaresCardinailty) {
      intMod = + static_cast<int>(std::round(std::sqrt(static_cast<float>(nSquared))));
      /*
    }
    else {
      intMod = - static_cast<int>(std::round(std::sqrt(static_cast<float>(nSquared))));
    }
    */
    data << n << "\t" << intMod << "\n";
    ++histogram[intMod];
  }


  for (int k{0}; k < primeNumber; ++k) {
    if (histogram[k] > 0) {
      hist << k << "\t" << histogram[k] << "\n";
    }
  }

  data.close(); non_squares.close(); hist.close();

  return 0;
}
