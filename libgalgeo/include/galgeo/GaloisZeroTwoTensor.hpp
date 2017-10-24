#ifndef GALOISZEROTWOTENSOR_HPP
#define GALOISZEROTWOTENSOR_HPP

#include <vector>
#include "GaloisOneOneTensor.hpp"
#include "GaloisZeroOneTensor.hpp"
#include "GaloisOneZeroTensor.hpp"
#include "GaloisScalar.hpp"


namespace galgeo
{

  class GaloisZeroTwoTensor
  {
    private:
      std::vector<std::vector<GaloisScalar> > components;

    public:
      static int prime;
      static unsigned int dimension;

      // Constructors and destructor
      GaloisZeroTwoTensor();
      GaloisZeroTwoTensor(const GaloisZeroTwoTensor& rhs);
      GaloisZeroTwoTensor(const std::vector<std::vector<GaloisScalar> >& rhs);
      GaloisZeroTwoTensor(const std::vector<std::vector<int> >& rhs);
      virtual ~GaloisZeroTwoTensor();

      // Binary operator for the assignment of a one_one_tensor
      GaloisZeroTwoTensor& operator=(const GaloisZeroTwoTensor& rhs);
      GaloisZeroTwoTensor& operator=(const std::vector<std::vector<GaloisScalar> >& rhs);
      GaloisZeroTwoTensor& operator=(const std::vector<std::vector<int > >& rhs);

      // Binary operators for checking (in-)equality
      bool operator==(const GaloisZeroTwoTensor& rhs);
      bool operator==(const std::vector<std::vector<int> >& rhs);
      bool operator!=(const GaloisZeroTwoTensor& rhs);
      bool operator!=(const std::vector<std::vector<int> >& rhs);

      // Binary operations for a one_one_tensor on the right hand side (rhs)
      GaloisZeroTwoTensor operator+(const GaloisZeroTwoTensor& rhs);
      GaloisZeroTwoTensor& operator+=(const GaloisZeroTwoTensor& rhs);
      GaloisZeroTwoTensor operator-(const GaloisZeroTwoTensor& rhs);
      GaloisZeroTwoTensor& operator-=(const GaloisZeroTwoTensor& rhs);

      // Making the two_zero_tensor act on a GaloisOneZeroTensor
      GaloisZeroOneTensor operator*(const GaloisOneZeroTensor& rhs);

      // Binary operation for a scalar on the right hand side (rhs)
      GaloisZeroTwoTensor operator*(const GaloisScalar& rhs);
      GaloisZeroTwoTensor& operator*=(const GaloisScalar& rhs);

      // Norm the GaloisZeroTwoTensor (allowed by the homogenity)
      GaloisZeroTwoTensor& norm(); 

      // Four methods for the accessing the elements (the second for preventing errors in case a const is expected)
      GaloisScalar& operator()(const unsigned int& row_index, const unsigned int& column_index);
      const GaloisScalar& operator()(const unsigned int& row_index, const unsigned int& column_index) const;
      std::vector<GaloisScalar>& operator[](const unsigned int& row_index);
      const std::vector<GaloisScalar>& operator[](const unsigned int& row_index) const;

      // Access the value of the dimenions
      unsigned int get_dimension() const;

      // Methods for transposing, inverting and asking for the determinent
      GaloisZeroTwoTensor& transpose();
      GaloisZeroTwoTensor get_transposed();
      GaloisScalar determinant();
      bool is_invertible();
      GaloisZeroTwoTensor& invert();
      GaloisZeroTwoTensor get_inverse();

      // Binary operators for assigning or reading the value via iostream
      friend std::istream& operator>>(std::istream& is, GaloisZeroTwoTensor& one_one_tensor_input);
      friend std::ostream& operator<<(std::ostream& os, const GaloisZeroTwoTensor& one_one_tensor_output);
  };

  //int determinant_of_matrix(const std::vector<std::vector<unsigned int> >& matrix);

  // Multiplying the GaloisZeroTwoTensor form the left with a GaloisOneOneTensor
  GaloisZeroTwoTensor operator*(const GaloisOneOneTensor lhs, const GaloisZeroTwoTensor rhs);
  GaloisZeroTwoTensor operator*(const GaloisZeroTwoTensor lhs, const GaloisOneOneTensor rhs);

}

#endif
