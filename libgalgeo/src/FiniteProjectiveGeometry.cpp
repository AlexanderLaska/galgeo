#ifndef FINITEPROJECTIVEGEOMETRY_CPP
#define FINITEPROJECTIVEGEOMETRY_CPP

#include "../include/galgeo/FiniteProjectiveGeometry.hpp"


namespace galgeo
{
  // Generate the Galois field
  std::vector<GaloisScalar> set_of_scalars(int& primeNumber)
  {
    std::vector<GaloisScalar> result;

    for(int i = 0; i < primeNumber; ++i)
    {
      result.push_back(i);
    }

    return result;
  }


  // Dehomogenization map (works only for dim =2)
  /*
  std::vector<std::vector<int> > dehomogenize(GaloisOneZeroTensor& point1, GaloisOneZeroTensor& point2, GaloisZeroTwoTensor& quadric) {
    
  }
  */

  // Generate a set of translation matrices
  std::vector<GaloisOneOneTensor> produce_translation_set(const unsigned int& number_of_points, const unsigned int& dimension, const unsigned int& prime_number)
  {
    std::vector<GaloisOneOneTensor> result(0);

    // This loop produces the affine translations
    for(size_t i=0; i<prime_number; ++i)
    {
      for(size_t j=0; j<prime_number; ++j)
      {
        GaloisOneOneTensor transformation;

        transformation[0][0] = 1; transformation[0][1] = 0; transformation[0][2] = i;
        transformation[1][0] = 0; transformation[1][1] = 1; transformation[1][2] = j;
        transformation[2][0] = 0; transformation[2][1] = 0; transformation[2][2] = 1;

        result.push_back(transformation);
      }
    }

    /* Successive Translations
    for(size_t i=0; i<prime_number; ++i)
    {
      GaloisOneOneTensor transformation;

      transformation[0][0] = 0; transformation[0][1] = 1; transformation[0][2] = i;
      transformation[1][0] = 0; transformation[1][1] = 0; transformation[1][2] = 1;
      transformation[2][0] = 1; transformation[2][1] = 0; transformation[2][2] = 0;

      result.push_back(transformation);
    }
    */

    /* pm t approach
    for(size_t i=0; i<prime_number-1; ++i)
    {
      GaloisOneOneTensor transformation;

      transformation[0][0] = 4*(1+i); transformation[0][1] = -4*(1+i); transformation[0][2] = i-1;
      transformation[1][0] = -4*(i-1); transformation[1][1] = 4*(i-1); transformation[1][2] = i+1;
      transformation[2][0] = 1; transformation[2][1] = 1; transformation[2][2] = 0;

      if(i != 1)
      {
        result.push_back(transformation);
      }
    }
      */

    /* the assymetric aproach
    for(size_t i=1; i<prime_number; ++i)
    {
      GaloisOneOneTensor transformation;

      transformation[0][0] = 1; transformation[0][1] = 0; transformation[0][2] = i;
      transformation[1][0] = 0; transformation[1][1] = 1; transformation[1][2] = 1;
      transformation[2][0] = -i; transformation[2][1] = -1; transformation[2][2] = 0;

      result.push_back(transformation);
    }

    GaloisOneOneTensor transformation;
    */

    /* t inverse involving approach
    for(size_t i=1; i<prime_number; ++i)
    {
      GaloisOneOneTensor transformation;
      GaloisScalar t{i};

      transformation[0][0] = 1; transformation[0][1] = -t.get_inverse(); transformation[0][2] = t;
      transformation[1][0] = -t; transformation[1][1] = 1; transformation[1][2] = 1;
      transformation[2][0] = -t.get_inverse(); transformation[2][1] = -1; transformation[2][2] = 0;

      result.push_back(transformation);
    }
    */

    /* Rotate
    for(size_t i=0; i<prime_number; ++i)
    {
      GaloisOneOneTensor transformation;
      GaloisScalar t{i};

      transformation[0][0] = 0; transformation[0][1] = 1; transformation[0][2] = t;
      transformation[1][0] = 0; transformation[1][1] = -t; transformation[1][2] = 1;
      transformation[2][0] = 1; transformation[2][1] = 1; transformation[2][2] = 0;

      result.push_back(transformation);
    }

    */
    
    GaloisScalar oneHalf{2};
    oneHalf.invert();

    for(size_t i=1; i<prime_number; ++i)
    {
      GaloisOneOneTensor transformation;

      GaloisScalar l{static_cast<unsigned int>(i)};
      /*
      GaloisScalar sign{+1};

      GaloisScalar rootOfTwo{2};
      rootOfTwoInverse.invert();

      if(2*i > prime_number)
      {
        sign = -1;
      }

      GaloisScalar exponent{2};

      if(test_nth_root_existence(rootOfTwo, exponent) != true)
      {
        rootOfTwo *= -1;
        rootOfTwo = rootOfTwo.take_nth_root(exponent);
      }
      else
      {
        rootOfTwo= rootOfTwo.take_nth_root(exponent);
      }

      transformation[0][0] = l; transformation[0][1] = 0; transformation[0][2] = l;
      transformation[1][0] = 0; transformation[1][1] = l.get_inverse(); transformation[1][2] = sign*l.get_inverse();
      transformation[2][0] = sign; transformation[2][1] = 1; transformation[2][2] = 0;

      transformation[0][0] = l*rootOfTwoInverse; transformation[0][1] = -sign*l*rootOfTwoInverse; transformation[0][2] = l;
      transformation[1][0] = -sign*l.get_inverse()*rootOfTwoInverse; transformation[1][1] = sign*l.get_inverse()*rootOfTwoInverse; transformation[1][2] = l.get_inverse();
      transformation[2][0] = sign; transformation[2][1] = 1; transformation[2][2] = 0;
      */

      GaloisScalar bPlus{l*l - l - oneHalf};
      GaloisScalar bMinus{l*l - l + oneHalf};

      transformation[0][0] = l /*l.get_inverse()*/; transformation[0][1] = l /*l.get_inverse()*/; transformation[0][2] = 1;
      transformation[1][0] = bPlus; transformation[1][1] = bMinus; transformation[1][2] = l;
      transformation[2][0] = 1; transformation[2][1] = 1; transformation[2][2] = 0;

      //std::cout << transformation << std::endl;
      result.push_back(transformation);
    }

    GaloisOneOneTensor transformation;

    /* pm t approach
    transformation[0][0] = 4; transformation[0][1] = -4; transformation[0][2] = 1;
    transformation[1][0] = -4; transformation[1][1] = 4; transformation[1][2] = 1;
    transformation[2][0] = 1; transformation[2][1] = 1; transformation[2][2] = 0;

    result.push_back(transformation);

    transformation[0][0] = -1; transformation[0][1] = 0; transformation[0][2] = 1;
    transformation[1][0] = -oneHalf; transformation[1][1] = 1; transformation[1][2] = 0;
    transformation[2][0] = 1; transformation[2][1] = 0; transformation[2][2] = 0;

    result.push_back(transformation);

    transformation[0][0] = 1; transformation[0][1] = 1; transformation[0][2] = 0;
    transformation[1][0] = 0; transformation[1][1] = -oneHalf; transformation[1][2] = 1;
    transformation[2][0] = 0; transformation[2][1] = 1; transformation[2][2] = 0;

    result.push_back(transformation);

    for(GaloisOneOneTensor& matrix : result)
    {
      std::cout << "---\n" << matrix << std::endl;
    }
    */

    return result;
  }

  // Produce a histogram of the number of the number of incidences
  std::vector<unsigned int> get_histogram(std::vector<unsigned int>& number_of_incidences_field)
  {
    std::vector<unsigned int> result;

    unsigned int max_number_of_incidence{0};
    for(unsigned int& number_of_incidences : number_of_incidences_field)
    {
      if(number_of_incidences > max_number_of_incidence)
      {
        max_number_of_incidence = number_of_incidences;
      }
    }

    for(unsigned int i=0; i<max_number_of_incidence + 1; i++)
    {
      unsigned int value_count{0};

      for(unsigned int& number_of_incidences : number_of_incidences_field)
      {
        if(number_of_incidences == i)
        {
          ++value_count;
        }
      }
      result.push_back(value_count);
    }

    return result;
  }

  // Test whether a point is the center of a pair of quadrics
  bool is_center(GaloisOneZeroTensor& center_point, std::vector<GaloisZeroTwoTensor>& biquadric, std::vector<GaloisOneZeroTensor>& points, std::vector<GaloisZeroOneTensor>& hyperplanes)
  {
    std::vector<GaloisZeroOneTensor> hyperplanes_through_center{set_of_objects_incident_with<GaloisOneZeroTensor, GaloisZeroOneTensor>(center_point, hyperplanes)};


    std::vector<GaloisOneZeroTensor> points_of_biquadric{points_in_biquadric(biquadric, points)};

    unsigned int count{0};
    //unsigned int no_alibi{0};

    for(GaloisZeroOneTensor& hyperplane_through_center : hyperplanes_through_center)
    {
      count = 0;

      for(GaloisOneZeroTensor& point_of_biquadric : points_of_biquadric)
      {
        if(hyperplane_through_center(point_of_biquadric) == 0)
        {
          ++count;
        }

      }

      if(count != 2)
      {
        return false;
      }
/*
      if(count < 2)
      {
        ++no_alibi;
      }

      if((count != 2) && (no_alibi >= 2))
      {
        return false;
      }
*/
    }

    return true;
  }

  // Test whether there is a number that is a nth root
  bool test_nth_root_existence(GaloisScalar& to_be_tested, GaloisScalar& exponent)
  {
    for(unsigned int i=0; i<GaloisScalar::prime; i++)
    {
      if(pow(i, exponent) == to_be_tested)
      {
        return true;
      }
    }

    return false;
  }

  // Count incidences
  std::vector<unsigned int> count_incidences(std::vector<std::vector<GaloisZeroTwoTensor> >& biquadric_distribution, std::vector<GaloisOneZeroTensor>& points)
  {
    std::vector<unsigned int> result(biquadric_distribution.size(), 0);

   for(unsigned int i=0; i<points.size(); i++)
   {
     for(unsigned int j=0; j<points.size(); j++)
     {
       if((biquadric_distribution[j][0]*points[i])(points[i]) == 0)
       {
         ++result[i];
       }
       else
       {
         if((biquadric_distribution[j][1]*points[i])(points[i]) == 0)
         {
           ++result[i];
         }
       }
     }
   }

    return result;
  }

  // Produce a distribution of point sets
  std::vector<std::vector<GaloisOneZeroTensor> > get_one_point_set_distribution(std::vector<GaloisOneZeroTensor>& initial_point_set, std::vector<GaloisOneOneTensor>& transformation_distribution, std::vector<GaloisOneZeroTensor> points)
  {
    unsigned int number_of_points{static_cast<unsigned int>(transformation_distribution.size())};

    std::vector<std::vector<GaloisOneZeroTensor> > result;

    result.resize(number_of_points);
    for(unsigned int i=0; i<number_of_points; i++)
    {
      result[i].resize(initial_point_set.size());
      GaloisOneOneTensor temp{transformation_distribution[i].get_inverse()};

      for(unsigned int j=0; j<initial_point_set.size(); j++)
      {
        result[i][j] = transformation_distribution[i]*initial_point_set[j];
      }
    }

    return result;
  }

  // Produce a distribution of biquadrics
  std::vector<std::vector<GaloisZeroTwoTensor> > get_one_biquadric_distribution(std::vector<GaloisZeroTwoTensor>& initial_biquadric, std::vector<GaloisOneOneTensor>& transformation_distribution, std::vector<GaloisOneZeroTensor> points)
  {
    std::vector<std::vector<GaloisZeroTwoTensor> > result;

    result.resize(transformation_distribution.size());
    for(unsigned int i=0; i<transformation_distribution.size(); i++)
    {
      result[i].resize(2);
      GaloisOneOneTensor temp{transformation_distribution[i].get_inverse()};
      result[i][0] = (temp.get_transposed()*initial_biquadric[0])*temp;
      //std::cout << temp.get_transposed() << "*\n" << initial_biquadric[0] << "*\n" << temp << "=\n" << temp.get_transposed()*initial_biquadric[0] << "*\n" << temp << "=\n" << (temp.get_transposed()*initial_biquadric[0])*temp << std::endl;
      result[i][1] = (temp.get_transposed()*initial_biquadric[1])*temp;
      //std::cout << temp.get_transposed() << "*\n" << initial_biquadric[1] << "*\n" << temp << "=\n" << temp.get_transposed()*initial_biquadric[1] << "*\n" << temp << "=\n" << (temp.get_transposed()*initial_biquadric[1])*temp << std::endl;
    }

    return result;
  }

  // Set of random transformations to all points in the projective geometry
  std::vector<GaloisOneOneTensor> get_one_transformation_distribution(std::vector<GaloisOneZeroTensor>& points)
  {
    //unsigned int number_of_points{static_cast<unsigned int>(points.size())};
    //unsigned int dimension{GaloisOneZeroTensor::dimension};

    /*
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, number_of_points-1);
    */

    std::vector<GaloisOneOneTensor> result;

    /*
    for(unsigned int i=0; i<number_of_points; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
        result[i][j][dimension-1] = points[i][j];
      }

      do
      {
        for(unsigned int k=0; k<dimension-1; k++)
        {
          unsigned int random_index{static_cast<unsigned int>(dis(gen))};

          for(unsigned int j=0; j<dimension; j++)
          {
            result[i][j][k] = points[random_index][j];
          }
        }
      }
      while(result[i].determinant() == 0);
    }
    */

    GaloisOneZeroTensor initial_center_point;
    initial_center_point[0] = 0;
    initial_center_point[1] = 0;
    initial_center_point[2] = 1;

    for(GaloisOneZeroTensor& new_center_point : points)
    {
      GaloisOneOneTensor new_transformation;
      
      while((new_transformation*initial_center_point).norm() != new_center_point.norm() || new_transformation.determinant() != 0)
      {
        new_transformation.randomize();
      }

      result.push_back(new_transformation);
    }

    return result;
  }

  // Find the points that are incident with a quadric or a biquadric
  std::vector<GaloisOneZeroTensor> points_in_quadric(GaloisZeroTwoTensor& quadric, std::vector<GaloisOneZeroTensor>& points)
  {
    std::vector<GaloisOneZeroTensor> result;

    for(unsigned int i=0; i<points.size(); i++)
    {
      if((quadric*points[i])(points[i]) == 0)
      {
        result.push_back(points[i]);
      }
    }

    return result;
  }

  std::vector<GaloisOneZeroTensor> points_in_biquadric(std::vector<GaloisZeroTwoTensor>& biquadric, std::vector<GaloisOneZeroTensor>& points)
  {
    std::vector<GaloisOneZeroTensor> result{points_in_quadric(biquadric[0], points)};

    for(unsigned int i=0; i<points.size(); i++)
    {
      if((biquadric[1]*points[i])(points[i]) == 0)
      {
        result.push_back(points[i]);
      }
    }

    return result;
  }

  // Produces all points in the projective geometry
  std::vector<GaloisOneZeroTensor> set_of_points(const int& prime_number, const unsigned int& dimension, const unsigned int& number_of_points)
  {
    std::vector<GaloisOneZeroTensor> points;
    points.resize(number_of_points);
 
    unsigned int i{0};
    for(unsigned int number_of_zeros=0; number_of_zeros<dimension; number_of_zeros++)
    {
      std::vector<std::vector<unsigned int> > same_number_of_zeros;
      same_number_of_zeros.resize(1);
      same_number_of_zeros[0].resize(number_of_zeros + 1);
 
      for(unsigned int j=0; j<number_of_zeros; j++)
      {
        same_number_of_zeros[0][j] = 0;
      }
 
      same_number_of_zeros[0][number_of_zeros] = 1;
 
      for(unsigned int j=number_of_zeros+1; j<dimension; j++)
      {
        same_number_of_zeros = add_layer(same_number_of_zeros, prime_number);
      }
 
      for(unsigned int j=0; j<same_number_of_zeros.size(); j++)
      {
 
        points[i] = same_number_of_zeros[j];
 
        i++;
      }
    }

    return points;
  }

  std::vector<GaloisZeroOneTensor> set_of_hyperplanes(std::vector<GaloisOneZeroTensor>& points)
  {
    std::vector<GaloisZeroOneTensor> hyperplanes;
    hyperplanes.resize(points.size());

    for(unsigned int i=0; i<points.size(); i++)
    {
      for(unsigned int j=0; j<GaloisOneZeroTensor::dimension; j++)
      {
        hyperplanes[i][j] = points[i][j];
      }
    }

    return hyperplanes;
  }

  
  // Adds a layer of prime numbers (cartesian product of a given old_list with the Galois field)
  std::vector<std::vector<unsigned int> > add_layer(std::vector<std::vector<unsigned int> >& old_list, const int& prime_number)
  {
    std::vector<std::vector<unsigned int> > result;
    result.resize(old_list.size()*prime_number);

    unsigned int l{0};
    for(unsigned int i=0; i<old_list.size(); i++)
    {
      for(unsigned int j=0; j<prime_number; j++)
      {
        result[l].resize(old_list[0].size() + 1);

        for(unsigned int k=0; k<result[l].size() - 1; k++)
        {
          result[l][k] = old_list[i][k];
        }
        result[l][result[l].size()-1] = j;

        l++;
      }
    }

    return result;
  }
}

#endif
