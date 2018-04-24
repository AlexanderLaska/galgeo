#ifndef FINITEPROJECTIVEGEOMETRY_HPP
#define FINITEPROJECTIVEGEOMETRY_HPP

#include <iostream>
#include <math.h>
#include <random>

#include "GaloisScalar.hpp"
#include "GaloisOneZeroTensor.hpp"
#include "GaloisZeroOneTensor.hpp"
#include "GaloisOneOneTensor.hpp"
#include "GaloisZeroTwoTensor.hpp"
#include "GaloisTwoZeroTensor.hpp"


namespace galgeo
{
  // Generate the Galois field
  std::vector<GaloisScalar> set_of_scalars(int& primeNumber);

  // Generate a set of translation matrices
  std::vector<GaloisOneOneTensor> produce_translation_set(const unsigned int& number_of_points, const unsigned int& dimension, const unsigned int& prime_number);

  // Produce a histogram of the number of the number of incidences
  std::vector<unsigned int> get_histogram(std::vector<unsigned int>& number_of_incidences_field);

  // Test whether a point is the center of a biquadric
  bool is_center(GaloisOneZeroTensor& center_point, std::vector<GaloisZeroTwoTensor>& biquadric, std::vector<GaloisOneZeroTensor>& points, std::vector<GaloisZeroOneTensor>& hyperplanes);

  // Set of all objects out of the second argument that are incident with the first object
  template<typename T1, typename T2>
    std::vector<T2> set_of_objects_incident_with(T1& object, std::vector<T2>& dual_objects)
    {
      std::vector<T2> result;

      for(T2 &dual_object : dual_objects)
      {
        if(dual_object*object == 0)
        {
          result.push_back(dual_object);
        }
      }

      return result;
    }

  // Find the set of nonquadratric GaloisScalars
  bool test_nth_root_existence(GaloisScalar& to_be_tested, GaloisScalar& exponent);

  // Count incidences
  std::vector<unsigned int> count_incidences(std::vector<std::vector<GaloisZeroTwoTensor> >& biquadric_distribution, std::vector<GaloisOneZeroTensor>& points);

  // Produce a distribution of point sets
  std::vector<std::vector<GaloisOneZeroTensor> > get_one_point_set_distribution(std::vector<GaloisOneZeroTensor>& initial_point_set, std::vector<GaloisOneOneTensor>& transformation_distribution, std::vector<GaloisOneZeroTensor> points);

  // Produce a distribution of biquadrics
  std::vector<std::vector<GaloisZeroTwoTensor> > get_one_biquadric_distribution(std::vector<GaloisZeroTwoTensor>& initial_biquadric, std::vector<GaloisOneOneTensor>& transformation_distribution, std::vector<GaloisOneZeroTensor> points);

  // Set of random transformations to all points in the projective geometry
  std::vector<GaloisOneOneTensor> get_one_transformation_distribution(std::vector<GaloisOneZeroTensor>& points);

  // Find the points that are incident with a (bi)quadric
  std::vector<GaloisOneZeroTensor> points_in_quadric(GaloisZeroTwoTensor& quadric, std::vector<GaloisOneZeroTensor>& points);
  std::vector<GaloisOneZeroTensor> points_in_biquadric(std::vector<GaloisZeroTwoTensor>& biquadric, std::vector<GaloisOneZeroTensor>& points);

  // Produces all points and hyperplanes in the projective geometry
  std::vector<GaloisOneZeroTensor> set_of_points(const int& prime_number, const unsigned int& dimension, const unsigned int& number_of_points);
  std::vector<GaloisZeroOneTensor> set_of_hyperplanes(std::vector<GaloisOneZeroTensor>& points);

  // Adds a layer of GaloisScalars (cartesian product of a given old_list with the Galois field)
  std::vector<std::vector<int> > add_layer(std::vector<std::vector<int> >& old_list, const int& prime_number);
}

#endif
