#ifndef GALOISTWOZEROTENSOR_HPP
#define GALOISTWOZEROTENSOR_HPP

#include <vector>
#include "GaloisOneOneTensor.hpp"
#include "GaloisZeroOneTensor.hpp"
#include "GaloisOneZeroTensor.hpp"
#include "GaloisScalar.hpp"


namespace galgeo
{

  class GaloisTwoZeroTensor
  {
    private:
      std::vector<std::vector<GaloisScalar> > components;

    public:
      static int prime;
      static unsigned int dimension;

      // Constructors and destructor
      GaloisTwoZeroTensor();
      GaloisTwoZeroTensor(const GaloisTwoZeroTensor& rhs);
      virtual ~GaloisTwoZeroTensor();

      // Binary operator for the assignment of a one_one_tensor
      GaloisTwoZeroTensor& operator=(const GaloisTwoZeroTensor& rhs);

      // Binary operators for checking (in-)equality
      bool operator==(const GaloisTwoZeroTensor& rhs);
      bool operator==(const std::vector<std::vector<int> >& rhs);
      bool operator!=(const GaloisTwoZeroTensor& rhs);
      bool operator!=(const std::vector<std::vector<int> >& rhs);

      // Binary operations for a one_one_tensor on the right hand side (rhs)
      GaloisTwoZeroTensor operator+(const GaloisTwoZeroTensor& rhs);
      GaloisTwoZeroTensor& operator+=(const GaloisTwoZeroTensor& rhs);
      GaloisTwoZeroTensor operator-(const GaloisTwoZeroTensor& rhs);
      GaloisTwoZeroTensor& operator-=(const GaloisTwoZeroTensor& rhs);

      GaloisTwoZeroTensor operator*(const GaloisOneOneTensor& rhs);
      GaloisTwoZeroTensor& operator*=(const GaloisOneOneTensor& rhs);

      // Making the two_zero_tensor act on a GaloisZeroOneTensor
      GaloisOneZeroTensor operator*(const GaloisZeroOneTensor& rhs);

      // Binary operation for a scalar on the right hand side (rhs)
      GaloisTwoZeroTensor operator*(const GaloisScalar& rhs);
      GaloisTwoZeroTensor& operator*=(const GaloisScalar& rhs);

      // Norm the GaloisTwoZeroTensor (allowed by the homogenity)
      GaloisTwoZeroTensor& norm(); 

      // Four methods for the accessing the elements (the second for preventing errors in case a const is expected)
      GaloisScalar& operator()(const unsigned int& row_index, const unsigned int& column_index);
      const GaloisScalar& operator()(const unsigned int& row_index, const unsigned int& column_index) const;
      std::vector<GaloisScalar>& operator[](const unsigned int& row_index);
      const std::vector<GaloisScalar>& operator[](const unsigned int& row_index) const;

      // Access the value of the dimenions
      unsigned int get_dimension() const;

      // Methods for transposing, inverting and asking for the determinent
      GaloisTwoZeroTensor& transpose();
      GaloisTwoZeroTensor get_transposed();
      GaloisScalar determinant();
      bool is_invertible();
      GaloisTwoZeroTensor& invert();
      GaloisTwoZeroTensor get_inverse();

      // Binary operators for assigning or reading the value via iostream
      friend std::istream& operator>>(std::istream& is, GaloisTwoZeroTensor& one_one_tensor_input);
      friend std::ostream& operator<<(std::ostream& os, const GaloisTwoZeroTensor& one_one_tensor_output);
  };

  //int determinant_of_matrix(const std::vector<std::vector<unsigned int> >& matrix);

  // Multiplying the GaloisTwoZeroTensor form the left with a GaloisOneOneTensor
  GaloisTwoZeroTensor operator*(GaloisOneOneTensor& lhs, GaloisTwoZeroTensor& rhs);
  GaloisTwoZeroTensor operator*(GaloisTwoZeroTensor& lhs, GaloisOneOneTensor& rhs);

}

#endif
