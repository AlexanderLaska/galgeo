#ifndef GALOISZEROONETENSOR_HPP
#define GALOISZEROONETENSOR_HPP

#include <vector>
#include "GaloisScalar.hpp"
#include "GaloisOneZeroTensor.hpp"
#include "GaloisOneOneTensor.hpp"


namespace galgeo
{

  class GaloisZeroOneTensor
  {
    private:
      std::vector<GaloisScalar> components;

    public:
      // Constructors and destructor
      GaloisZeroOneTensor();
      GaloisZeroOneTensor(const GaloisZeroOneTensor& rhs);
      virtual ~GaloisZeroOneTensor();

      // Binary operator for the assignment of another GaloisZeroOneTensor or std:: vector (rhs)
      GaloisZeroOneTensor& operator=(const GaloisZeroOneTensor& rhs);
      GaloisZeroOneTensor& operator=(const std::vector<unsigned int>& rhs);

      // Binary operators for checking (in-)equality
      bool operator==(const GaloisZeroOneTensor& rhs);
      bool operator==(const std::vector<int>& rhs);
      bool operator!=(const GaloisZeroOneTensor& rhs);
      bool operator!=(const std::vector<int>& rhs);

      // Binary operations for a GaloisZeroOneTensor on the right hand side (rhs)
      GaloisZeroOneTensor operator+(const GaloisZeroOneTensor& rhs);
      GaloisZeroOneTensor& operator+=(const GaloisZeroOneTensor& rhs);
      GaloisZeroOneTensor operator-(const GaloisZeroOneTensor& rhs);
      GaloisZeroOneTensor& operator-=(const GaloisZeroOneTensor& rhs);

      // Binary operations for a vector of integers on the right hand side (rhs)
      GaloisZeroOneTensor operator+(const std::vector<int>& rhs);
      GaloisZeroOneTensor& operator+=(const std::vector<int>& rhs);
      GaloisZeroOneTensor operator-(const std::vector<int>& rhs);
      GaloisZeroOneTensor& operator-=(const std::vector<int>& rhs);
      
      // Multiply the covector with a GaloisScalar
      GaloisZeroOneTensor operator*(const GaloisScalar& factor);
      GaloisZeroOneTensor& operator*=(const GaloisScalar& factor);

      // Multiply the covector with an integer
      GaloisZeroOneTensor operator*(const int& factor);
      GaloisZeroOneTensor& operator*=(const int& factor);

      // Let the GaloisZeroOneTensor act on GaloisScalar to execute their scalar product
      GaloisScalar operator()(const GaloisOneZeroTensor& argument);
      
      // Making the GaloisZeroOneTensor act on an one_one_tensor
      GaloisZeroOneTensor operator()(const GaloisOneOneTensor& rhs);

      // Norm the GaloisZeroOneTensor (allowed by the homogenity)
      GaloisZeroOneTensor& norm(); 

      // Two methods for the accessing the elements (the second for preventing errors in case a const is expected)
      GaloisScalar& operator[](const unsigned int& component_index);
      const GaloisScalar& operator[](const unsigned int& component_index) const; 
   
      // Access the value of the dimension
      unsigned int get_dimension() const;

      // Binary operators for assigning or reading the values via iostream
      friend std::istream& operator>>(std::istream& is, GaloisZeroOneTensor& covector_input);
      friend std::ostream& operator<<(std::ostream& os, const GaloisZeroOneTensor& covector_output);
  };

  // Vector and Covector Product
  GaloisOneZeroTensor vector_product(std::vector<GaloisZeroOneTensor>& factors);
  GaloisZeroOneTensor covector_product(std::vector<GaloisOneZeroTensor>& factors);

  // Scalar product
  GaloisScalar operator*(GaloisOneZeroTensor lhs, GaloisZeroOneTensor rhs);
  GaloisScalar operator*(GaloisZeroOneTensor lhs, GaloisOneZeroTensor rhs);

}

#endif
