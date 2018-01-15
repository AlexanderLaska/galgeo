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

  std::ofstream q, q_bar, points, centerPoints, line, background, non_squares;
  q.open("q.txt"); q_bar.open("q_bar.txt"); points.open("points.txt"); centerPoints.open("centerPoints.txt"); line.open("line.txt"); background.open("background.txt"); non_squares.open("non_squares.txt");

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

  int numberOfPoints{primeNumber*primeNumber + primeNumber + 1};
  std::vector<GaloisOneZeroTensor> allPoints{set_of_points(primeNumber,dimension+1,numberOfPoints)};
  std::vector<GaloisZeroOneTensor> allHyperplanes{set_of_hyperplanes(allPoints)};

  GaloisScalar zero{0}; GaloisScalar one{1};

  GaloisZeroTwoTensor M, Mbar;
  M[0][0] = std::atoi(argv[2]); M[0][1] = std::atoi(argv[3]); M[0][2] = std::atoi(argv[4]);
  M[1][0] = std::atoi(argv[3]); M[1][1] = std::atoi(argv[5]); M[1][2] = std::atoi(argv[6]);
  M[2][0] = std::atoi(argv[4]); M[2][1] = std::atoi(argv[6]); M[2][2] = std::atoi(argv[7]);
  Mbar[0][0] = std::atoi(argv[8]); Mbar[0][1] = std::atoi(argv[9]); Mbar[0][2] = std::atoi(argv[10]);
  Mbar[1][0] = std::atoi(argv[9]); Mbar[1][1] = std::atoi(argv[11]); Mbar[1][2] = std::atoi(argv[12]);
  Mbar[2][0] = std::atoi(argv[10]); Mbar[2][1] = std::atoi(argv[12]); Mbar[2][2] = std::atoi(argv[13]);

  std::vector<GaloisOneZeroTensor> pointsOfM{points_in_quadric(M, allPoints)};
  std::vector<GaloisOneZeroTensor> pointsOfMbar{points_in_quadric(Mbar, allPoints)};

  GaloisZeroTwoTensor Mdual{M.get_inverse()}, Mdualbar{Mbar.get_inverse()};

  std::vector<GaloisZeroTwoTensor> quadricPair;
  quadricPair.push_back(M);
  quadricPair.push_back(Mbar);
  std::vector<GaloisZeroTwoTensor> dualQuadricPair;
  dualQuadricPair.push_back(Mdual);
  dualQuadricPair.push_back(Mdualbar);

  points << "The pair of quadrics\n\n" << M << "\n" << Mbar << "\n\n" << "with its dual pair of quadrics\n\n" << Mdual << "\n" << Mdualbar << "\n\n" << "is a biquadric with respect to these center points:\n\n";
  std::vector<GaloisOneZeroTensor> allCenterPoints, alldualCenterPoints;
  std::vector<std::vector<GaloisZeroOneTensor> > allHyperplanesThroughCenters;
  for (GaloisOneZeroTensor point : allPoints) {
    if (is_center(point, quadricPair, allPoints, allHyperplanes)) {
      points << point << "\n";
      allCenterPoints.push_back(point);
      allHyperplanesThroughCenters.push_back(set_of_objects_incident_with(point, allHyperplanes));
    }
  }
  points << "\n and all dual polars:\n\n";
  for (GaloisOneZeroTensor point : allCenterPoints) {
    points << M*point << "\n" << Mbar*point << "\n";
  }

  points << "Points of M:\n";
  for(GaloisOneZeroTensor point : pointsOfM) {
    if (point.norm()[2] != 0) {
      if ((point.norm()[0] <= squaresCardinality) && (point.norm()[1] <= squaresCardinality)) {
        q << point.norm()[0] << "\t" << point.norm()[1] << "\n";
      }
      else if ((point.norm()[0] > squaresCardinality) && (point.norm()[1] <= squaresCardinality)) {
        q << static_cast<int>(point.norm()[0]) - primeNumber << "\t" << point.norm()[1] << "\n";
      }
      else if ((point.norm()[0] <= squaresCardinality) && (point.norm()[1] > squaresCardinality)) {
        q << point.norm()[0] << "\t" << static_cast<int>(point.norm()[1]) - primeNumber << "\n";
      }
      else if ((point.norm()[0] > squaresCardinality) && (point.norm()[1] > squaresCardinality)) {
        q << static_cast<int>(point.norm()[0]) - primeNumber << "\t" << static_cast<int>(point.norm()[1]) - primeNumber << "\n";
      }
    }
    else if ((point.norm()[2] == 0) && (point.norm()[1] != 0)) {
      if (point.norm()[0] <= squaresCardinality) {
        q << squaresCardinality + 2 << "\t" << point.norm()[0] << "\n";
      }
      else {
        q << squaresCardinality + 2 << "\t" << static_cast<int>(point.norm()[0]) - primeNumber << "\n";
      }
    }
    else if ((point.norm()[2] == 0) && (point.norm()[1] == 0) && (point.norm()[0] == 1)) {
      q << squaresCardinality + 2 << "\t" << squaresCardinality + 2 << std::endl;
    }
    points << point << std::endl;
  }
  points << "Points of Mbar:\n";
  for(GaloisOneZeroTensor point : pointsOfMbar) {
    if (point.norm()[2] != 0) {
      if ((point.norm()[0] <= squaresCardinality) && (point.norm()[1] <= squaresCardinality)) {
        q_bar << point.norm()[0] << "\t" << point.norm()[1] << "\n";
      }
      else if ((point.norm()[0] > squaresCardinality) && (point.norm()[1] <= squaresCardinality)) {
        q_bar << static_cast<int>(point.norm()[0]) - primeNumber << "\t" << point.norm()[1] << "\n";
      }
      else if ((point.norm()[0] <= squaresCardinality) && (point.norm()[1] > squaresCardinality)) {
        q_bar << point.norm()[0] << "\t" << static_cast<int>(point.norm()[1]) - primeNumber << "\n";
      }
      else if ((point.norm()[0] > squaresCardinality) && (point.norm()[1] > squaresCardinality)) {
        q_bar << static_cast<int>(point.norm()[0]) - primeNumber << "\t" << static_cast<int>(point.norm()[1]) - primeNumber << "\n";
      }
    }
    else if ((point.norm()[2] == 0) && (point.norm()[1] != 0)) {
      if (point.norm()[0] <= squaresCardinality) {
        q_bar << squaresCardinality + 2 << "\t" << point.norm()[0] << "\n";
      }
      else {
        q_bar << squaresCardinality + 2 << "\t" << static_cast<int>(point.norm()[0]) - primeNumber << "\n";
      }
    }
    else if ((point.norm()[2] == 0) && (point.norm()[1] == 0) && (point.norm()[0] == 1)) {
      q_bar << squaresCardinality + 2 << "\t" << squaresCardinality + 2 << std::endl;
    }
    points << point << std::endl;
  }

  for(GaloisOneZeroTensor point : allCenterPoints) {
    if (point.norm()[2] != 0) {
      if ((point.norm()[0] <= squaresCardinality) && (point.norm()[1] <= squaresCardinality)) {
        centerPoints << point.norm()[0] << "\t" << point.norm()[1] << "\n";
      }
      else if ((point.norm()[0] > squaresCardinality) && (point.norm()[1] <= squaresCardinality)) {
        centerPoints << static_cast<int>(point.norm()[0]) - primeNumber << "\t" << point.norm()[1] << "\n";
      }
      else if ((point.norm()[0] <= squaresCardinality) && (point.norm()[1] > squaresCardinality)) {
        centerPoints << point.norm()[0] << "\t" << static_cast<int>(point.norm()[1]) - primeNumber << "\n";
      }
      else if ((point.norm()[0] > squaresCardinality) && (point.norm()[1] > squaresCardinality)) {
        centerPoints << static_cast<int>(point.norm()[0]) - primeNumber << "\t" << static_cast<int>(point.norm()[1]) - primeNumber << "\n";
      }
    }
    else if ((point.norm()[2] == 0) && (point.norm()[1] != 0)) {
      if (point.norm()[0] <= squaresCardinality) {
        centerPoints << squaresCardinality + 2 << "\t" << point.norm()[0] << "\n";
      }
      else {
        centerPoints << squaresCardinality + 2 << "\t" << static_cast<int>(point.norm()[0]) - primeNumber << "\n";
      }
    }
    else if ((point.norm()[2] == 0) && (point.norm()[1] == 0) && (point.norm()[0] == 1)) {
      centerPoints << squaresCardinality + 2 << "\t" << squaresCardinality + 2 << std::endl;
    }
  }

  std::vector<GaloisOneZeroTensor> pointsOnLine{set_of_objects_incident_with(allHyperplanesThroughCenters[std::atoi(argv[14])][std::atoi(argv[15])], allPoints)};
  for(GaloisOneZeroTensor point : pointsOnLine) {
    if (point.norm()[2] != 0) {
      if ((point.norm()[0] <= squaresCardinality) && (point.norm()[1] <= squaresCardinality)) {
        line << point.norm()[0] << "\t" << point.norm()[1] << "\n";
      }
      else if ((point.norm()[0] > squaresCardinality) && (point.norm()[1] <= squaresCardinality)) {
        line << static_cast<int>(point.norm()[0]) - primeNumber << "\t" << point.norm()[1] << "\n";
      }
      else if ((point.norm()[0] <= squaresCardinality) && (point.norm()[1] > squaresCardinality)) {
        line << point.norm()[0] << "\t" << static_cast<int>(point.norm()[1]) - primeNumber << "\n";
      }
      else if ((point.norm()[0] > squaresCardinality) && (point.norm()[1] > squaresCardinality)) {
        line << static_cast<int>(point.norm()[0]) - primeNumber << "\t" << static_cast<int>(point.norm()[1]) - primeNumber << "\n";
      }
    }
    else if ((point.norm()[2] == 0) && (point.norm()[1] != 0)) {
      if (point.norm()[0] <= squaresCardinality) {
        line << squaresCardinality + 2 << "\t" << point.norm()[0] << "\n";
      }
      else {
        line << squaresCardinality + 2 << "\t" << static_cast<int>(point.norm()[0]) - primeNumber << "\n";
      }
    }
    else if ((point.norm()[2] == 0) && (point.norm()[1] == 0) && (point.norm()[0] == 1)) {
      line << squaresCardinality + 2 << "\t" << squaresCardinality + 2 << std::endl;
    }
  }
  for (int i{-squaresCardinality}; i <= squaresCardinality; ++i) {
    for (int j{-squaresCardinality}; j <= squaresCardinality; ++j) {
      background << i << "\t" << j << "\n";
    }
    background << squaresCardinality + 2 << "\t" << i << "\n";
  }
  background << squaresCardinality + 2 << "\t" << squaresCardinality + 2 << "\n";

  q.close(); q_bar.close(); points.close(); centerPoints.close(); line.close(); background.close(); non_squares.close();

  return 0;
}
