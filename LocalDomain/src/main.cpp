#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <random>
#include <algorithm>
#include <cstdio>
#include <ctime>
#include <tuple>

#include <GaloisScalar.hpp>
#include <GaloisOneZeroTensor.hpp>
#include <GaloisZeroOneTensor.hpp>
#include <GaloisOneOneTensor.hpp>
#include <GaloisZeroTwoTensor.hpp>
#include <GaloisTwoZeroTensor.hpp>

#include <FiniteProjectiveGeometry.hpp>

using namespace galgeo;

int print(std::string message){
  std::cout << message << std::endl;
  return 0;
}
int printint(int i){
  std::cout << i;
  return 0;
}

int main(int argc, char *argv[]) {

/* open output files */
  std::ofstream q, q_bar, qdash, q_bardash, non_squares, background, all_points_with_counter, loc, difference_x, difference_y;
  q.open("q.txt"); q_bar.open("q_bar.txt"); qdash.open("qdash.txt"); q_bardash.open("q_bardash.txt"); non_squares.open("non_squares.txt"); background.open("background.txt"); all_points_with_counter.open("all_points_with_counter.txt"); loc.open("loc.txt"); difference_x.open("difference_x.txt"); difference_y.open("difference_y.txt");
  
  
/* initially define usefull objects */

  // use the first command line argument as the prime number
  int primeNumber = std::atoi(argv[1]); // std::atoi converts strings to ints
  set_prime(primeNumber);

  // choose dimension 2
  int dimension = 2;
  set_dimension(dimension+1);
  
  // secure the zero and one are GaloisScalars and act like ones #safetyFirst
  GaloisScalar zero{0}; GaloisScalar one{1};
  
  // define a vector of all points that we have
  std::vector<GaloisOneZeroTensor> all_points = set_of_points(primeNumber, dimension+1, primeNumber*primeNumber + primeNumber + 1);

  // define size of the set of squares
  int squaresCardinality = static_cast<int>((primeNumber - 1)*0.5); //static_cast<new_type>(arg) safely converts the argument to the new_type
  
  // define a vector of the squares
  std::vector<GaloisScalar> squares;
  for(GaloisScalar candidate = 1; candidate <= squaresCardinality; ++candidate) 
     squares.push_back(candidate*candidate);
  
  // define a vector for the non squares
  std::vector<GaloisScalar> nonSquares;
  
  // fill the vector of non squares and print them to the file
  non_squares << "Numbers without a square root for p = " << primeNumber << ":" << std::endl;
  for(int candidate = 1; candidate < primeNumber; ++candidate) {
    if(std::find(squares.begin(), squares.end(), candidate) != squares.end()) {
    } else {
      nonSquares.push_back(candidate);
      non_squares << candidate << " ";
    }
  }
  
  // define a vector containing all points in the local domain in the e basis, without infinity points
  std::vector<GaloisOneZeroTensor> localDomain_e;
  for(GaloisOneZeroTensor point : all_points) {
    if(point[0] != 0) {
      int p1 = (point[1] <= squaresCardinality) ? static_cast<int>(point[1]) : static_cast<int>(point[1])-primeNumber;
      int p2 = (point[2] <= squaresCardinality) ? static_cast<int>(point[2]) : static_cast<int>(point[2])-primeNumber;
      //if(((p1*p1) + (p2*p2)) < primeNumber) // ! calculate this with the components seen as from Z (as ints) !!!
      if(((p1*p1) < primeNumber) && ((p2*p2) < primeNumber)) // ! calculate this with the components seen as from Z (as ints) !!!
        localDomain_e.push_back(point);
    }
  }

  // define matrix elements of the quadric matrices
  GaloisScalar m11 = std::atoi(argv[2]);
  GaloisScalar m12 = std::atoi(argv[3]);
  GaloisScalar m22 = std::atoi(argv[4]);
  GaloisScalar f = nonSquares[std::atoi(argv[5])];

  
/* define quadric matrices */

  GaloisZeroTwoTensor Qplus; // the [1] and the [2] component are the ones of the affine space
  
  Qplus[0][0] = 1; Qplus[1][0] = 0;   Qplus[2][0] = 0;
  Qplus[0][1] = 0; Qplus[1][1] = m11; Qplus[2][1] = m12;
  Qplus[0][2] = 0; Qplus[1][2] = m12; Qplus[2][2] = m22;

  GaloisZeroTwoTensor Qminus{Qplus};
  Qminus[0][0] = f;
  
/* print determinant of the quadric */

  GaloisScalar determinant = Qplus.determinant();
  std::cout << "the determinant of Qplus ( " << determinant;
  if(std::find(squares.begin(), squares.end(), determinant) != squares.end()) {
    print(" ) is a square.");
  } else {
    print(" ) is a non-square.");
  }
  
/* find points in both quadrics */

  std::vector<GaloisOneZeroTensor> points_in_Qplus = points_in_quadric(Qplus, all_points);
  std::vector<GaloisOneZeroTensor> points_in_Qminus = points_in_quadric(Qminus, all_points);
  
/* print background and quadric points to the file w/out infinity points */

  for(GaloisOneZeroTensor point : all_points) {
    if(point[0] != 0) {
      int p1 = (point[1] <= squaresCardinality) ? static_cast<int>(point[1]) : static_cast<int>(point[1])-primeNumber;
      int p2 = (point[2] <= squaresCardinality) ? static_cast<int>(point[2]) : static_cast<int>(point[2])-primeNumber;
      background << p1 << "\t" << p2 << std::endl;}
    else { }
  }

  for(GaloisOneZeroTensor point : points_in_Qplus) {
    if(point[0] != 0) {
      int p1 = (point[1] <= squaresCardinality) ? static_cast<int>(point[1]) : static_cast<int>(point[1])-primeNumber;
      int p2 = (point[2] <= squaresCardinality) ? static_cast<int>(point[2]) : static_cast<int>(point[2])-primeNumber;
      q << p1 << "\t" << p2 << std::endl;
      qdash << p1 << "\t" << p2 << "\t" << "1" << std::endl;}
    else { }
  }

  for(GaloisOneZeroTensor point : points_in_Qminus) {
    if(point[0] != 0) {
      int p1 = (point[1] <= squaresCardinality) ? static_cast<int>(point[1]) : static_cast<int>(point[1])-primeNumber;
      int p2 = (point[2] <= squaresCardinality) ? static_cast<int>(point[2]) : static_cast<int>(point[2])-primeNumber;
      q_bar << p1 << "\t" << p2 << std::endl;
      q_bardash << p1 << "\t" << p2 << "\t" << "2" << std::endl;}
    else { }
  }
  
  
/* calculate all unit vectors possible from Qplus and Qminus */

  std::vector<GaloisOneZeroTensor> e_from_Qplus(points_in_Qplus);
  std::vector<GaloisOneZeroTensor> e_from_Qminus(points_in_Qminus);
  
  // print("all possible e vectors from Qplus");
  for(GaloisOneZeroTensor point : e_from_Qplus){
    point[1] *= one/f;
    point[2] *= one/f;
    // std::cout << point << std::endl;
  }
  
  // print("all possible e vectors from Qminus");
  for(GaloisOneZeroTensor point : e_from_Qminus){
    point[1] *= one/f;
    point[2] *= one/f;
    // std::cout << point << std::endl;
  }
  
  
/* huge loop over all possible combinations of biquadric point tuples : for one from Qplus and one from Qminus */

  std::vector<int> counter(all_points.size(), 0); // counter has lenght number_of_points and is filled with zeros
  
  for(std::vector<GaloisOneZeroTensor>::size_type i = 0; i != e_from_Qplus.size(); ++i){
    for(std::vector<GaloisOneZeroTensor>::size_type j = 0; j != e_from_Qplus.size(); ++j){
        
      //first basis vector is from Qplus, second one from Qminus
      GaloisOneZeroTensor e0 = e_from_Qplus[i];
      GaloisOneZeroTensor e1 = e_from_Qminus[j];
      
      // coordinate trafo from the local domain in the e basis to the original basis
      std::vector<GaloisOneZeroTensor> localDomain_q;
      for(GaloisOneZeroTensor point : localDomain_e){
        GaloisOneZeroTensor p;
        p[0] = point[0];
        p[1] = point[1] * e0[1] + point[2] * e1[1];
        p[2] = point[1] * e0[2] + point[2] * e1[2];
        localDomain_q.push_back(p);
      }
      
      // loop over all point to see which ones'd be transformed into the local domain, ++ their counter.
      for(std::vector<GaloisOneZeroTensor>::size_type k = 0; k != all_points.size(); ++k){
        if( find(localDomain_q.begin(), localDomain_q.end(), all_points[k]) != localDomain_q.end() ) // check if this condition works as I want it to !
            ++counter[k];
      }
       
       
      // write localDomain_q for some random i, j values to the output to be used as example
      if(i == 3 && j == 2){
        for(GaloisOneZeroTensor point : localDomain_q){
          if(point[0] != 0) {
            int p1 = (point[1] <= squaresCardinality) ? static_cast<int>(point[1]) : static_cast<int>(point[1])-primeNumber;
            int p2 = (point[2] <= squaresCardinality) ? static_cast<int>(point[2]) : static_cast<int>(point[2])-primeNumber;
            loc << p1 << "\t" << p2 << std::endl;
          }
        }
      }
      
    /*
      // first basis vector is from Qminus, second one from Qplus
      e0 = e_from_Qminus[i];
      e1 = e_from_Qplus[j];
      
      // coordinate trafo from the local domain in the e basis to the original basis
      localDomain_q.clear();
      for(GaloisOneZeroTensor point : localDomain_e){
        GaloisOneZeroTensor p;
        p[0] = point[0];
        p[1] = point[1] * e0[1] + point[2] * e1[1];
        p[2] = point[1] * e0[2] + point[2] * e1[2];
        localDomain_q.push_back(p);
      }
      
      // loop over all point to see which ones'd be transformed into the local domain, ++ their counter.
      for(std::vector<GaloisOneZeroTensor>::size_type k = 0; k != all_points.size(); ++k){
        if( find(localDomain_q.begin(), localDomain_q.end(), all_points[k]) != localDomain_q.end() ) // check if this condition works as I want it to !
            ++counter[k];
      }
    */
      
    }
  } 

  
/* huge loop over all possible pais of biquadric points * /

  // define a vector containing all possible basis vectors of the combined biquadric
  std::vector<GaloisOneZeroTensor> e_from_biquadric(e_from_Qplus);
  for(GaloisOneZeroTensor point : e_from_Qminus){
    e_from_biquadric.push_back(point);
  }
  
  // define a counter
  std::vector<int> counter(all_points.size(), 0); // counter has lenght number_of_points and is filled with zeros
  
  // do the loop, j <= i to avoid double counting
  for(std::vector<GaloisOneZeroTensor>::size_type i = 0; i != e_from_Qplus.size(); ++i){
    for(std::vector<GaloisOneZeroTensor>::size_type j = 0; j < i; ++j){
        
      // basis vectors from the biquadric
      GaloisOneZeroTensor e0 = e_from_biquadric[i];
      GaloisOneZeroTensor e1 = e_from_biquadric[j];
    
      if(e0[1] != -e1[1] && e0[2] != -e1[2]){
        // coordinate trafo from the local domain in the e basis to the original basis
        std::vector<GaloisOneZeroTensor> localDomain_q;
        for(GaloisOneZeroTensor point : localDomain_e){
          GaloisOneZeroTensor p;
          p[0] = point[0];
          p[1] = point[1] * e0[1] + point[2] * e1[1];
          p[2] = point[1] * e0[2] + point[2] * e1[2];
          localDomain_q.push_back(p);
        }
      
        // loop over all point to see which ones'd be transformed into the local domain, ++ their counter.
        for(std::vector<GaloisOneZeroTensor>::size_type k = 0; k != all_points.size(); ++k){
          if( find(localDomain_q.begin(), localDomain_q.end(), all_points[k]) != localDomain_q.end() ) // check if this condition works as I want it to !
            ++counter[k];
        }
       
       
        // write localDomain_q for some random i, j values to the output to be used as example
        if(i == 3 && j == 2){
          for(GaloisOneZeroTensor point : localDomain_q){
            if(point[0] != 0) {
              int p1 = (point[1] <= squaresCardinality) ? static_cast<int>(point[1]) : static_cast<int>(point[1])-primeNumber;
              int p2 = (point[2] <= squaresCardinality) ? static_cast<int>(point[2]) : static_cast<int>(point[2])-primeNumber;
              loc << p1 << "\t" << p2 << std::endl;
            }
          }
        } 
      }
      
    }
  } 
*/

  
/* print counter */

  for(std::vector<GaloisOneZeroTensor>::size_type i = 0; i != counter.size(); ++i){
    //printint(i);
    //std::cout << "\t" << counter[i] << std::endl;
    if(all_points[i][0] != 0){
      int p1 = (all_points[i][1] <= squaresCardinality) ? static_cast<int>(all_points[i][1]) : static_cast<int>(all_points[i][1])-primeNumber;
      int p2 = (all_points[i][2] <= squaresCardinality) ? static_cast<int>(all_points[i][2]) : static_cast<int>(all_points[i][2])-primeNumber;
      if(p1==0 && p2==0){
        all_points_with_counter << p1 << "\t" << p2 << "\t"  << 0 << std::endl;}
      else{
        all_points_with_counter << p1 << "\t" << p2 << "\t"  << counter[i] << std::endl;}
    }
  }
  
  
/* see weather the heatmap is symmetric in x or y direction */
  
  for(std::vector<GaloisOneZeroTensor>::size_type i = 0; i != counter.size(); ++i){
    if(all_points[i][0] != 0){
      int p1_1 = (all_points[i][1] <= squaresCardinality) ? static_cast<int>(all_points[i][1]) : static_cast<int>(all_points[i][1])-primeNumber;
      int p2_1 = (all_points[i][2] <= squaresCardinality) ? static_cast<int>(all_points[i][2]) : static_cast<int>(all_points[i][2])-primeNumber;
      //
      for(std::vector<GaloisOneZeroTensor>::size_type j = 0; j != counter.size(); ++j){
        if(all_points[j][0] != 0){
          int p1_2 = (all_points[j][1] <= squaresCardinality) ? static_cast<int>(all_points[j][1]) : static_cast<int>(all_points[j][1])-primeNumber;
          int p2_2 = (all_points[j][2] <= squaresCardinality) ? static_cast<int>(all_points[j][2]) : static_cast<int>(all_points[j][2])-primeNumber;
          //
          if(p1_1 >= 0 && p1_1 == -p1_2 && p2_1 == p2_2){ 
            difference_x << p1_1 << "\t" << p2_1 << "\t"  << counter[i] - counter[j] << std::endl;
            // relative: difference_x << p1_1 << "\t" << p2_1 << "\t"  << static_cast<double>(counter[i]-counter[j])/counter[i] << std::endl;
          }
          if(p2_1 >= 0 && p2_1 == -p2_2 && p1_1 == p1_2){ 
            difference_y << p1_1 << "\t" << p2_1 << "\t"  << counter[i] - counter[j] << std::endl;
            // relative: difference_y << p1_1 << "\t" << p2_1 << "\t"  << static_cast<double>(counter[i]-counter[j])/counter[i] << std::endl;
          }
        }
      }
    }
  }
  
  
/* close the files */

  q.close(); q_bar.close(); qdash.close(); q_bardash.close(); non_squares.close(); background.close(); all_points_with_counter.close(); loc.close(); difference_x.close(); difference_y.close();

  
/* return */

  return 0;
}


















