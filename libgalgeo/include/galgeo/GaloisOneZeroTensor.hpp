#ifndef GALOISONEZEROTENSOR_HPP
#define GALOISONEZEROTENSOR_HPP

#include <vector>
#include "GaloisScalar.hpp"


namespace galgeo
{

  class GaloisOneZeroTensor
  {
    private:
      std::vector<GaloisScalar> components;

    public:
      static unsigned int dimension;

      // Constructors and destructor
      GaloisOneZeroTensor();
      GaloisOneZeroTensor(const GaloisOneZeroTensor& rhs);
      virtual ~GaloisOneZeroTensor();

      // Binary operators for checking (in-)equality
      bool operator==(const GaloisOneZeroTensor& rhs);
      bool operator==(const std::vector<int>& rhs);
      bool operator!=(const GaloisOneZeroTensor& rhs);
      bool operator!=(const std::vector<int>& rhs);

      // Binary operator for the assignment of another GaloisOneZeroTensor or std:: vector (rhs)
      GaloisOneZeroTensor& operator=(const GaloisOneZeroTensor& rhs);
      GaloisOneZeroTensor& operator=(const std::vector<unsigned int>& rhs);

      // Binary operations for a GaloisOneZeroTensor on the right hand side (rhs)
      GaloisOneZeroTensor operator+(const GaloisOneZeroTensor& rhs);
      GaloisOneZeroTensor& operator+=(const GaloisOneZeroTensor& rhs);
      GaloisOneZeroTensor operator-(const GaloisOneZeroTensor& rhs);
      GaloisOneZeroTensor& operator-=(const GaloisOneZeroTensor& rhs);

      // Binary operations for an integer on the right hand side (rhs)
      GaloisOneZeroTensor operator+(const std::vector<int>& rhs);
      GaloisOneZeroTensor& operator+=(const std::vector<int>& rhs);
      GaloisOneZeroTensor operator-(const std::vector<int>& rhs);
      GaloisOneZeroTensor& operator-=(const std::vector<int>& rhs);
      
      // Multiply the vector with a GaloisScalar
      GaloisOneZeroTensor operator*(const GaloisScalar& factor);
      GaloisOneZeroTensor& operator*=(const GaloisScalar& factor);

      // Multiply the vector with an integer
      GaloisOneZeroTensor operator*(const int& factor);
      GaloisOneZeroTensor& operator*=(const int& factor);

      // Norm the GaloisOneZeroTensor (allowed by the homogenity)
      GaloisOneZeroTensor& norm(); 
      GaloisOneZeroTensor get_normed(); 


      // Two methods for the accessing the elements (the second for preventing errors in case a const is expected)
      GaloisScalar& operator[](const unsigned int& component_index);
      const GaloisScalar& operator[](const unsigned int& component_index) const; 
   
      // Access the value of the dimension
      unsigned int get_dimension() const;
      friend void set_dimension(const unsigned int& new_dimension);

      // Binary operators for assigning or reading the values via iostream
      friend std::istream& operator>>(std::istream& is, GaloisOneZeroTensor& vector_input);
      friend std::ostream& operator<<(std::ostream& os, const GaloisOneZeroTensor& vector_output);
  };

  // Set the dimension
  void set_dimension(const unsigned int& new_dimension);

}

#endif
