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
  return 0; }
  
int printint(int i){
  std::cout << i << std::endl;
  return 0; }
  

int main(int argc, char *argv[]) {

/* open output files */

  std::ofstream q, q_bar, qdash, q_bardash, non_squares, background, all_points_with_counter_standard_basis, all_points_with_counter_own_basis, example_set;
  
  q.open("q.txt"); q_bar.open("q_bar.txt"); qdash.open("qdash.txt"); q_bardash.open("q_bardash.txt"); non_squares.open("non_squares.txt"); background.open("background.txt"); all_points_with_counter_standard_basis.open("all_points_with_counter_standard_basis.txt"); all_points_with_counter_own_basis.open("all_points_with_counter_own_basis.txt"); example_set.open("example_set.txt");
  
  
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
  double border = sqrt(squaresCardinality);
  
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
  // this local domain still contains the 'bad' quadric points, they are excluded when transformed into the standard basis
  std::vector<GaloisOneZeroTensor> localDomain_e;
  for(GaloisOneZeroTensor p : all_points) {
    if(p[0] != 0) {
      int p_1 = (p[1] <= squaresCardinality) ? static_cast<int>(p[1]) : static_cast<int>(p[1])-primeNumber;
      int p_2 = (p[2] <= squaresCardinality) ? static_cast<int>(p[2]) : static_cast<int>(p[2])-primeNumber;
      // if(p_\nu <= sqrt((p-1)/2) for all \nu) // ! calculate this with the components seen as from Z (as ints) !!!
      if(abs(p_1) <= border && abs(p_2) <= border) 
        localDomain_e.push_back(p);
    }
  }

  // define matrix elements of the quadric matrices
  GaloisScalar m11 = std::atoi(argv[2]);
  GaloisScalar m12 = std::atoi(argv[3]);
  GaloisScalar m22 = std::atoi(argv[4]);
  GaloisScalar f = nonSquares[std::atoi(argv[5])];

  
/* print infos on console */

  std::cout << "p = " << primeNumber << "\nf = " << f << "\nm11,m12,m22 = " << m11 << "," << m12 << "," << m22 << std::endl;
  
  
/* define quadric matrices */

  GaloisZeroTwoTensor Qplus; // the [1] and the [2] component are the ones of the affine space
  
  Qplus[0][0] = 1; Qplus[1][0] = 0;   Qplus[2][0] = 0;
  Qplus[0][1] = 0; Qplus[1][1] = m11; Qplus[2][1] = m12;
  Qplus[0][2] = 0; Qplus[1][2] = m12; Qplus[2][2] = m22;

  GaloisZeroTwoTensor Qminus(Qplus);
  Qminus[0][0] = f;
  
  
/* find points in both quadrics */

  std::vector<GaloisOneZeroTensor> points_in_Qplus = points_in_quadric(Qplus, all_points);
  std::vector<GaloisOneZeroTensor> points_in_Qminus = points_in_quadric(Qminus, all_points);
  

/* calculate all unit vectors possible from Qplus and Qminus */

  std::vector<GaloisOneZeroTensor> e_from_Qplus(points_in_Qplus);
  std::vector<GaloisOneZeroTensor> e_from_Qminus(points_in_Qminus);
  
  for(GaloisOneZeroTensor point : e_from_Qplus){
    point[1] *= one/f;
    point[2] *= one/f;
  }
  
  for(GaloisOneZeroTensor point : e_from_Qminus){
    point[1] *= one/f;
    point[2] *= one/f;
  }
  
  
/* define the selected local world domain (i=0 and j=0) represented in the standard basis */

  std::vector<GaloisOneZeroTensor> selected_localDomain;
  
  GaloisOneZeroTensor e0_selected_localDomain = e_from_Qplus[0];
  GaloisOneZeroTensor e1_selected_localDomain = e_from_Qminus[0];
      
  // coordinate trafo from the local domain in the e basis to the standard basis
  for(GaloisOneZeroTensor p : localDomain_e){
    // transform into the standard basis
    GaloisOneZeroTensor point;
    point[0] = p[0];
    point[1] = p[1] * e0_selected_localDomain[1] + p[2] * e1_selected_localDomain[1];
    point[2] = p[1] * e0_selected_localDomain[2] + p[2] * e1_selected_localDomain[2];
    // check weather point is a quadric point that does not make the basis, exclude
    bool not_in_Qplus = ( find(points_in_Qplus.begin(), points_in_Qplus.end(), point) == points_in_Qplus.end() );
    bool equal_q_to_e0 = ( point == points_in_Qplus[0] );// here I trust that the order of the points in e_from_Qplus is still the same as in points_in_Qplus
    // if (p notin Qplus  XOR  p == non-normalized vector of e0 )
    if(not_in_Qplus != equal_q_to_e0){
      bool not_in_Qminus = ( find(points_in_Qminus.begin(), points_in_Qminus.end(), point) == points_in_Qminus.end() );
      bool equal_q_to_e1 = ( point == points_in_Qminus[0] );
      // if (p notin Qminus  XOR  p == non-normalized vector of e1 )
      if(not_in_Qminus != equal_q_to_e1){
        selected_localDomain.push_back(point);
      }
    }
  }
  
  /* // write selected basis vectors to file 
  auto e0_selected_localDomain_str = "(" + std::to_string(e0_selected_localDomain[1]) + ", " + std::to_string(e0_selected_localDomain[2]) + ", " + std::to_string(e0_selected_localDomain[0]) + ")";
  auto e1_selected_localDomain_str = "(" + std::to_string(e1_selected_localDomain[1]) + ", " + std::to_string(e1_selected_localDomain[2]) + ", " + std::to_string(e1_selected_localDomain[0]) + ")";
  all_points_with_counter_own_basis << "\"The basis vectors for the selected local world domain are " << e0_selected_localDomain_str << " and " << e1_selected_localDomain_str << ".\"" << std::endl; */
  
  
/* huge loop over all possible combinations of biquadric point tuples : for one from Qplus and one from Qminus */

  std::vector<int> counter(selected_localDomain.size(), 0); // counter has lenght selected_localDomain and is filled with zeros
  
  for(std::vector<GaloisOneZeroTensor>::size_type i = 0; i != e_from_Qplus.size(); ++i){
    for(std::vector<GaloisOneZeroTensor>::size_type j = 0; j != e_from_Qplus.size(); ++j){
      if (i != 0 && j !=0){
          
        // first basis vector is from Qplus, second one from Qminus
        GaloisOneZeroTensor e0 = e_from_Qplus[i];
        GaloisOneZeroTensor e1 = e_from_Qminus[j];
      
        // coordinate trafo from the local domain in the e basis to the standard basis, there exclude the non-special quadric points
        std::vector<GaloisOneZeroTensor> localDomain_q;
        for(GaloisOneZeroTensor p : localDomain_e){
          GaloisOneZeroTensor point;
          point[0] = p[0];
          point[1] = p[1] * e0[1] + p[2] * e1[1];
          point[2] = p[1] * e0[2] + p[2] * e1[2];
          
          bool not_in_Qplus = ( find(points_in_Qplus.begin(), points_in_Qplus.end(), point) == points_in_Qplus.end() );
          bool equal_q_to_e0 = ( point == points_in_Qplus[i] );
          if(not_in_Qplus != equal_q_to_e0){
            bool not_in_Qminus = ( find(points_in_Qminus.begin(), points_in_Qminus.end(), point) == points_in_Qminus.end() );
            bool equal_q_to_e1 = ( point == points_in_Qminus[j] );
            if(not_in_Qminus != equal_q_to_e1){
              localDomain_q.push_back(point);
            }
          }
        }
      
        // loop over the current local domain, and search for each point of it in the selected local domain. If found, increse the counter value of the same position
        for(GaloisOneZeroTensor point : localDomain_q){
          std::vector<GaloisOneZeroTensor>::iterator k = find(selected_localDomain.begin(), selected_localDomain.end(), point);
          if(k != selected_localDomain.end()){
            int index = std::distance(selected_localDomain.begin(), k);
            ++counter[index];
          }
        } 
       
        /* does (should do ^^') the same: / loop over all point of the selected local domain to see which ones'd be transformed into the current local domain, ++ their counter.
        for(std::vector<GaloisOneZeroTensor>::size_type k = 0; k != selected_localDomain.size(); ++k){
          bool in_localDomain_q = ( find(localDomain_q.begin(), localDomain_q.end(), selected_localDomain[k]) != localDomain_q.end() );
          if(in_localDomain_q)
            ++counter[k];
        } */
        
        // write localDomain_q for some random i, j values to the output to be used as example
        if(i == 3 && j == 2){
          auto e0_str = "(" + std::to_string(e0[1]) + ", " + std::to_string(e0[2]) + ", " + std::to_string(e0[0]) + ")";
          auto e1_str = "(" + std::to_string(e1[1]) + ", " + std::to_string(e1[2]) + ", " + std::to_string(e1[0]) + ")";
          example_set << "\"The basis vectors for this example local world domain are " << e0_str << " and " << e1_str << ".\"" << std::endl;
          for(GaloisOneZeroTensor point : localDomain_q){
            if(point[0] != 0) {
              int point_1 = (point[1] <= squaresCardinality) ? static_cast<int>(point[1]) : static_cast<int>(point[1])-primeNumber;
              int point_2 = (point[2] <= squaresCardinality) ? static_cast<int>(point[2]) : static_cast<int>(point[2])-primeNumber;
              example_set << point_1 << "\t" << point_2 << std::endl;
            }
          }
        }
        
      }
    }
  } 

  
/* print counter */

  for(std::vector<GaloisOneZeroTensor>::size_type i = 0; i != counter.size(); ++i){
    GaloisOneZeroTensor point = selected_localDomain[i];
    // in the standard basis
    if(point[0] != 0){
        
      int point_1 = (point[1] <= squaresCardinality) ? static_cast<int>(point[1]) : static_cast<int>(point[1])-primeNumber;
      int point_2 = (point[2] <= squaresCardinality) ? static_cast<int>(point[2]) : static_cast<int>(point[2])-primeNumber;
      
      if(point_1==0 && point_2==0){
        all_points_with_counter_standard_basis << point_1 << "\t" << point_2 << "\t"  << 0 << std::endl;}
      else{
        all_points_with_counter_standard_basis << point_1 << "\t" << point_2 << "\t"  << counter[i] << std::endl;}
        
      // and in the own basis
      // don't trust the order of localDomain_e and selected_localDomain to be the same because of the excluded quadric points -> perform trafo back explicitly
      
      GaloisOneZeroTensor p;
      GaloisScalar det_inv = one / (e0_selected_localDomain[1] * e1_selected_localDomain[2] - e0_selected_localDomain[2] * e1_selected_localDomain[1]);
      
      p[0] = point[0];
      p[1] = det_inv * (point[1] * e1_selected_localDomain[2] - point[2] * e1_selected_localDomain[1]);
      p[2] = det_inv * (point[2] * e0_selected_localDomain[1] - point[1] * e0_selected_localDomain[2]);
      
      int p_1 = (p[1] <= squaresCardinality) ? static_cast<int>(p[1]) : static_cast<int>(p[1])-primeNumber;
      int p_2 = (p[2] <= squaresCardinality) ? static_cast<int>(p[2]) : static_cast<int>(p[2])-primeNumber;
      
      if(p_1==0 && p_2==0){
        all_points_with_counter_own_basis << p_1 << "\t" << p_2 << "\t"  << 0 << std::endl;}
      else{
        all_points_with_counter_own_basis << p_1 << "\t" << p_2 << "\t"  << counter[i] << std::endl;}
    }
  }
  
    
/* print background and quadric points to the file w/out infinity points */

  for(GaloisOneZeroTensor point : all_points) {
    if(point[0] != 0) {
      int point_1 = (point[1] <= squaresCardinality) ? static_cast<int>(point[1]) : static_cast<int>(point[1])-primeNumber;
      int point_2 = (point[2] <= squaresCardinality) ? static_cast<int>(point[2]) : static_cast<int>(point[2])-primeNumber;
      background << point_1 << "\t" << point_2 << std::endl;
    }
  }

  for(GaloisOneZeroTensor point : points_in_Qplus) {
    if(point[0] != 0) {
      int point_1 = (point[1] <= squaresCardinality) ? static_cast<int>(point[1]) : static_cast<int>(point[1])-primeNumber;
      int point_2 = (point[2] <= squaresCardinality) ? static_cast<int>(point[2]) : static_cast<int>(point[2])-primeNumber;
      q << point_1 << "\t" << point_2 << std::endl;
      qdash << point_1 << "\t" << point_2 << "\t" << "1" << std::endl;
    }  
  }

  for(GaloisOneZeroTensor point : points_in_Qminus) {
    if(point[0] != 0) {
      int point_1 = (point[1] <= squaresCardinality) ? static_cast<int>(point[1]) : static_cast<int>(point[1])-primeNumber;
      int point_2 = (point[2] <= squaresCardinality) ? static_cast<int>(point[2]) : static_cast<int>(point[2])-primeNumber;
      q_bar << point_1 << "\t" << point_2 << std::endl;
      q_bardash << point_1 << "\t" << point_2 << "\t" << "2" << std::endl;
    }
  }
  
  
/* close the files */

  q.close(); q_bar.close(); qdash.close(); q_bardash.close(); non_squares.close(); background.close(); all_points_with_counter_standard_basis.close(); all_points_with_counter_own_basis.close(); example_set.close();
  
  
/* return */

  return 0;
}


















