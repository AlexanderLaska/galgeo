#ifndef GALOISONEONETENSOR_HPP
#define GALOISONEONETENSOR_HPP

#include <vector>
#include "GaloisOneZeroTensor.hpp"
#include "GaloisScalar.hpp"


namespace galgeo
{

  class GaloisOneOneTensor
  {
    private:
      std::vector<std::vector<GaloisScalar> > components;

    public:
      static unsigned int dimension;

      // Constructors and destructor
      GaloisOneOneTensor();
      GaloisOneOneTensor(const GaloisOneOneTensor& rhs);
      virtual ~GaloisOneOneTensor();

      // Binary operator for the assignment of a one_one_tensor
      GaloisOneOneTensor& operator=(const GaloisOneOneTensor& rhs);

      // Binary operators for checking (in-)equality
      bool operator==(const GaloisOneOneTensor& rhs);
      bool operator==(const std::vector<std::vector<int> >& rhs);
      bool operator!=(const GaloisOneOneTensor& rhs);
      bool operator!=(const std::vector<std::vector<int> >& rhs);

      // Binary operations for a one_one_tensor on the right hand side (rhs)
      GaloisOneOneTensor operator+(const GaloisOneOneTensor& rhs);
      GaloisOneOneTensor& operator+=(const GaloisOneOneTensor& rhs);
      GaloisOneOneTensor operator-(const GaloisOneOneTensor& rhs);
      GaloisOneOneTensor& operator-=(const GaloisOneOneTensor& rhs);
      GaloisOneOneTensor operator*(const GaloisOneOneTensor& rhs);
      GaloisOneOneTensor& operator*=(const GaloisOneOneTensor& rhs);

      // Making the one_one_tensor act on a GaloisOneZeroTensor
      GaloisOneZeroTensor operator*(const GaloisOneZeroTensor& rhs);

      // Binary operation for a scalar on the right hand side (rhs)
      GaloisOneOneTensor operator*(const GaloisScalar& rhs);
      GaloisOneOneTensor& operator*=(const GaloisScalar& rhs);

      // Randomize the values uniformely
      GaloisOneOneTensor& randomize(const int lower_limit = 0, const int upper_limit = GaloisScalar::prime - 1);

      // Norm the GaloisOneOneTensor (allowed by the homogenity)
      GaloisOneOneTensor& norm(); 

      // Four methods for the accessing the elements (the second for preventing errors in case a const is expected)
      GaloisScalar& operator()(const unsigned int& row_index, const unsigned int& column_index);
      const GaloisScalar& operator()(const unsigned int& row_index, const unsigned int& column_index) const;
      std::vector<GaloisScalar>& operator[](const unsigned int& row_index);
      const std::vector<GaloisScalar>& operator[](const unsigned int& row_index) const;

      // Access the value of the dimenions
      unsigned int get_dimension() const;

      // Methods for transposing, inverting and asking for the determinent
      GaloisOneOneTensor& transpose();
      GaloisOneOneTensor get_transposed();
      GaloisScalar determinant();
      bool is_invertible();
      GaloisOneOneTensor& invert();
      GaloisOneOneTensor get_inverse();

      // Binary operators for assigning or reading the value via iostream
      friend std::istream& operator>>(std::istream& is, GaloisOneOneTensor& one_one_tensor_input);
      friend std::ostream& operator<<(std::ostream& os, const GaloisOneOneTensor& one_one_tensor_output);
  };

  int determinant_of_matrix(const std::vector<std::vector<unsigned int> >& matrix);

}

#endif
