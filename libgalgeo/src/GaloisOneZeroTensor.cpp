#ifndef GALOISONEZEROTENSOR_CPP
#define GALOISONEZEROTENSOR_CPP

#include "../include/galgeo/GaloisOneZeroTensor.hpp"


namespace galgeo
{

  // Initialize the value as the prime number of the GaloisScalar class
  unsigned int GaloisOneZeroTensor::dimension = 3;

  // Standard Constructor (choosing a dimension of three corresponding to the projective plane)
  GaloisOneZeroTensor::GaloisOneZeroTensor()
  {
    components.resize(dimension);
    for(unsigned int i=0; i<dimension; i++)
    {
      components[i] = 0;
    }
  }

  // Copy Constructor
  GaloisOneZeroTensor::GaloisOneZeroTensor(const GaloisOneZeroTensor& rhs) : components(rhs.components) {}

  // (Virtual) Destructor
  GaloisOneZeroTensor::~GaloisOneZeroTensor() {}

  // Binary operator for the assignment of another GaloisOneZeroTensor
  GaloisOneZeroTensor& GaloisOneZeroTensor::operator=(const GaloisOneZeroTensor& rhs)
  {
    // Test whether the assigned is the current galois_vactor
    if(&rhs == this)
    {
      return *this;
    }

    // ...and fill it with the values of the assigned GaloisOneZeroTensor 
    for(unsigned int i=0; i<dimension; i++)
    {
      components[i] = rhs.components[i];
    }

    return *this;
  }

  GaloisOneZeroTensor& GaloisOneZeroTensor::operator=(const std::vector<int>& rhs)
  {
    // ...and fill it with the values of the assigned vector 
    for(unsigned int i=0; i<dimension; i++)
    {
      components[i] = ((rhs[i] % GaloisScalar::prime) + GaloisScalar::prime) % GaloisScalar::prime;
    }

    return *this;
  }

  // Binary operators for checking (in-)equality
  bool GaloisOneZeroTensor::operator==(const GaloisOneZeroTensor& rhs)
  {
    GaloisOneZeroTensor _rhs{rhs};
    _rhs.norm();

    GaloisOneZeroTensor _lhs{*this};
    _lhs.norm();

    for(unsigned int i=0; i<dimension; i++)
    {
      if(_lhs.components[i] != _rhs.components[i])
      //if(this->components[i] != _rhs.components[i])
      {
        return false;
      }
    }

    return true;
  }

  bool GaloisOneZeroTensor::operator==(const std::vector<int>& rhs)
  {
    GaloisOneZeroTensor _rhs;
    _rhs = rhs;
    _rhs.norm();

    GaloisOneZeroTensor _lhs{*this};
    _lhs.norm();

    for(unsigned int i=0; i<dimension; i++)
    {
      if(_lhs.components[i] != _rhs.components[i])
      //if(this->components[i] != _rhs.components[i])
      {
        return false;
      }
    }

    return true;
  }

  bool GaloisOneZeroTensor::operator!=(const GaloisOneZeroTensor& rhs)
  {
    for(unsigned int i=0; i<dimension; i++)
    {
      if(this->components[i] != rhs.components[i])
      {
        return true;
      }
    }

    return false;
  }

  bool GaloisOneZeroTensor::operator!=(const std::vector<int>& rhs)
  {
    for(unsigned int i=0; i<dimension; i++)
    {
      if(rhs[i] != this->components[i])
      {
        return true;
      }
    }

    return false;
  }

  // Addition of two vectors
  GaloisOneZeroTensor GaloisOneZeroTensor::operator+(const GaloisOneZeroTensor& rhs)
  {
    GaloisOneZeroTensor result;

    for(unsigned int i=0; i<dimension; i++)
    {
      result.components[i] = this->components[i] + rhs.components[i];
    }

    return result;
  }

  // Adding and assigning a second vector to this vector 
  GaloisOneZeroTensor& GaloisOneZeroTensor::operator+=(const GaloisOneZeroTensor& rhs)
  {
    for(unsigned int i=0; i<dimension; i++)
    {
      this->components[i] += rhs.components[i];
    }

    return *this;
  }

  // Subtraction of two matrices
  GaloisOneZeroTensor GaloisOneZeroTensor::operator-(const GaloisOneZeroTensor& rhs)
  {
    GaloisOneZeroTensor result;

    for(unsigned int i=0; i<dimension; i++)
    {
      result.components[i] = this->components[i] - rhs.components[i];
    }

    return result;
  }

  // Subtracting and assigning a second vector to this vector 
  GaloisOneZeroTensor& GaloisOneZeroTensor::operator-=(const GaloisOneZeroTensor& rhs)
  {
    for(unsigned int i=0; i<dimension; i++)
    {
      this->components[i] -= rhs.components[i];
    }

    return *this;
  }

  // Addition of a GaloisOneZeroTensor and a std::vector
  GaloisOneZeroTensor GaloisOneZeroTensor::operator+(const std::vector<int>& rhs)
  {
    GaloisOneZeroTensor result;

    for(unsigned int i=0; i<dimension; i++)
    {
      result.components[i] = this->components[i] + rhs[i];
    }

    return result;
  }

  // Adding and assigning a std::vector to this GaloisOneZeroTensor
  GaloisOneZeroTensor& GaloisOneZeroTensor::operator+=(const std::vector<int>& rhs)
  {
    for(unsigned int i=0; i<dimension; i++)
    {
      this->components[i] += rhs[i];
    }

    return *this;
  }

  // Addition of two matrices
  GaloisOneZeroTensor GaloisOneZeroTensor::operator-(const std::vector<int>& rhs)
  {
    GaloisOneZeroTensor result;

    for(unsigned int i=0; i<dimension; i++)
    {
      result.components[i] = this->components[i] - rhs[i];
    }

    return result;
  }

  // Adding and assigning a second matrix to this matrix
  GaloisOneZeroTensor& GaloisOneZeroTensor::operator-=(const std::vector<int>& rhs)
  {
    for(unsigned int i=0; i<dimension; i++)
    {
      this->components[i] -= rhs[i];
    }

    return *this;
  }

  // Multiply the vector with a GaloisScalar
  GaloisOneZeroTensor GaloisOneZeroTensor::operator*(const GaloisScalar& factor)
  {
    GaloisOneZeroTensor result(*this);

    for(unsigned int i=0; i<this->dimension; i++)
    {
      result.components[i] *= factor;
    }

    return result;
  }

  GaloisOneZeroTensor& GaloisOneZeroTensor::operator*=(const GaloisScalar& factor)
  {
    for(unsigned int i=0; i<this->dimension; i++)
    {
      this->components[i] *= factor;
    }

    return *this;
  }

  // Multiply the vector with an integer
  GaloisOneZeroTensor GaloisOneZeroTensor::operator*(const int& factor)
  {
    GaloisOneZeroTensor result(*this);

    for(unsigned int i=0; i<this->dimension; i++)
    {
      result.components[i] *= factor;
    }

    return result;
  }

  GaloisOneZeroTensor& GaloisOneZeroTensor::operator*=(const int& factor)
  {
    for(unsigned int i=0; i<this->dimension; i++)
    {
      this->components[i] *= factor;
    }

    return *this;
  }

  // Norm the GaloisOneZeroTensor (allowed by the homogenity)
  GaloisOneZeroTensor& GaloisOneZeroTensor::norm()
  {
    for(int i=GaloisOneZeroTensor::dimension - 1; i >= 0; i--)
    {
      if(this->components[i] != 0)
      {
        GaloisScalar factor((this->components[i]).get_inverse());
        for(unsigned int j=0; j<dimension; j++)
        {
          this->components[j] *= factor;
        }

        return *this;
      }
    }

    return *this;
  }

  GaloisOneZeroTensor GaloisOneZeroTensor::get_normed()
  {
    GaloisOneZeroTensor result{this->norm()};

    return result;
  }


  // Two methods for accessing the elements (the second for preventing errors in case a const is expected)
  GaloisScalar& GaloisOneZeroTensor::operator[](const unsigned int& component_index)
  {
    return this->components[component_index];
  }

  const GaloisScalar& GaloisOneZeroTensor::operator[](const unsigned int& component_index) const
  {
    return this->components[component_index];
  }

  // Get the dimension
  unsigned int GaloisOneZeroTensor::get_dimension() const
  {
    return this->dimension;
  }

  // Binary operators for assigning or reading the values via iostream
  std::istream& operator>>(std::istream& is, GaloisOneZeroTensor& vector_input)
  {
    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      is >> vector_input.components[i];
    }

    return is;
  }

  std::ostream& operator<<(std::ostream& os, const GaloisOneZeroTensor& vector_output)
  {
    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      os << vector_output.components[i];
      os << "\n";
    }

    return os;
  }

  // Set the dimension
  void set_dimension(const unsigned int& new_dimension)
  {
    GaloisOneZeroTensor::dimension = new_dimension;
  }

}

#endif
