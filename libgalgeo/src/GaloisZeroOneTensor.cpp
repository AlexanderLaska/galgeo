#ifndef GALOISZEROONETENSOR_CPP
#define GALOISZEROONETENSOR_CPP

#include "../include/galgeo/GaloisZeroOneTensor.hpp"


namespace galgeo
{

  // Standard Constructor (choosing a dimension of three corresponding to the projective plane)
  GaloisZeroOneTensor::GaloisZeroOneTensor()
  {
    components.resize(GaloisOneZeroTensor::dimension);
    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      components[i] = 0;
    }
  }

  // Copy Constructor
  GaloisZeroOneTensor::GaloisZeroOneTensor(const GaloisZeroOneTensor& rhs) : components(rhs.components) {}

  // (Virtual) Destructor
  GaloisZeroOneTensor::~GaloisZeroOneTensor() {}

  // Binary operator for the assignment of another GaloisZeroOneTensor
  GaloisZeroOneTensor& GaloisZeroOneTensor::operator=(const GaloisZeroOneTensor& rhs)
  {
    // Test whether the assigned is the current galois_vactor
    if(&rhs == this)
    {
      return *this;
    }

    // ...and fill it with the values of the assigned GaloisZeroOneTensor 
    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      components[i] = rhs.components[i];
    }

    return *this;
  }

  GaloisZeroOneTensor& GaloisZeroOneTensor::operator=(const std::vector<unsigned int>& rhs)
  {
    // ...and fill it with the values of the assigned vector 
    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      components[i] = ((rhs[i] % GaloisScalar::prime) + GaloisScalar::prime) % GaloisScalar::prime;
    }

    return *this;
  }

  // Binary operators for checking (in-)equality
  bool GaloisZeroOneTensor::operator==(const GaloisZeroOneTensor& rhs)
  {
    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      if(this->components[i] != rhs.components[i])
      {
        return false;
      }
    }

    return true;
  }

  bool GaloisZeroOneTensor::operator==(const std::vector<int>& rhs)
  {
    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      if(rhs[i] != this->components[i])
      {
        return false;
      }
    }

    return true;
  }

  bool GaloisZeroOneTensor::operator!=(const GaloisZeroOneTensor& rhs)
  {
    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      if(this->components[i] != rhs.components[i])
      {
        return true;
      }
    }

    return false;
  }

  bool GaloisZeroOneTensor::operator!=(const std::vector<int>& rhs)
  {
    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      if(rhs[i] != this->components[i])
      {
        return true;
      }
    }

    return false;
  }

  // Addition of two covectors
  GaloisZeroOneTensor GaloisZeroOneTensor::operator+(const GaloisZeroOneTensor& rhs)
  {
    GaloisZeroOneTensor result;

    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      result.components[i] = this->components[i] + rhs.components[i];
    }

    return result;
  }

  // Adding and assigning a second covector to this covector 
  GaloisZeroOneTensor& GaloisZeroOneTensor::operator+=(const GaloisZeroOneTensor& rhs)
  {
    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      this->components[i] += rhs.components[i];
    }

    return *this;
  }

  // Subtraction of two matrices
  GaloisZeroOneTensor GaloisZeroOneTensor::operator-(const GaloisZeroOneTensor& rhs)
  {
    GaloisZeroOneTensor result;

    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      result.components[i] = this->components[i] - rhs.components[i];
    }

    return result;
  }

  // Subtracting and assigning a second covector to this covector 
  GaloisZeroOneTensor& GaloisZeroOneTensor::operator-=(const GaloisZeroOneTensor& rhs)
  {
    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      this->components[i] -= rhs.components[i];
    }

    return *this;
  }

  // Addition of a GaloisZeroOneTensor and a std::vector
  GaloisZeroOneTensor GaloisZeroOneTensor::operator+(const std::vector<int>& rhs)
  {
    GaloisZeroOneTensor result;

    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      result.components[i] = this->components[i] + rhs[i];
    }

    return result;
  }

  // Adding and assigning a std::vector to this GaloisZeroOneTensor
  GaloisZeroOneTensor& GaloisZeroOneTensor::operator+=(const std::vector<int>& rhs)
  {
    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      this->components[i] += rhs[i];
    }

    return *this;
  }

  // Addition of two matrices
  GaloisZeroOneTensor GaloisZeroOneTensor::operator-(const std::vector<int>& rhs)
  {
    GaloisZeroOneTensor result;

    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      result.components[i] = this->components[i] - rhs[i];
    }

    return result;
  }

  // Adding and assigning a second matrix to this matrix
  GaloisZeroOneTensor& GaloisZeroOneTensor::operator-=(const std::vector<int>& rhs)
  {
    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      this->components[i] -= rhs[i];
    }

    return *this;
  }

  // Multiply the covector with a GaloisScalar
  GaloisZeroOneTensor GaloisZeroOneTensor::operator*(const GaloisScalar& factor)
  {
    GaloisZeroOneTensor result(*this);

    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      result.components[i] *= factor;
    }

    return result;
  }

  GaloisZeroOneTensor& GaloisZeroOneTensor::operator*=(const GaloisScalar& factor)
  {
    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      this->components[i] *= factor;
    }

    return *this;
  }

  // Multiply the covector with an integer
  GaloisZeroOneTensor GaloisZeroOneTensor::operator*(const int& factor)
  {
    GaloisZeroOneTensor result(*this);

    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      result.components[i] *= factor;
    }

    return result;
  }

  GaloisZeroOneTensor& GaloisZeroOneTensor::operator*=(const int& factor)
  {
    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      this->components[i] *= factor;
    }

    return *this;
  }

  // Let the GaloisZeroOneTensor act on GaloisScalar to execute their scalar product
  GaloisScalar GaloisZeroOneTensor::operator()(const GaloisOneZeroTensor& argument)
  {
    GaloisScalar result(0);

    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      result += (this->components[i] * argument[i]);
    }

    return result;
  }

  // Making the GaloisZeroOneTensor act on an one_one_tensor
  GaloisZeroOneTensor GaloisZeroOneTensor::operator()(const GaloisOneOneTensor& rhs)
  {
    GaloisZeroOneTensor result;

    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      for(unsigned int j=0; j<GaloisOneZeroTensor::dimension; j++)
          {
            result[i] = this->components[j] * rhs(i,j);
          }
    }

    return result;
  }

  // Norm the GaloisZeroOneTensor (allowed by the homogenity)
  GaloisZeroOneTensor& GaloisZeroOneTensor::norm()
  {
    for(int i=GaloisOneZeroTensor::dimension - 1; i >= 0; i--)
    {
      if(this->components[i] != 0)
      {
        GaloisScalar factor((this->components[i]).get_inverse());
        for(unsigned int j=0; j<GaloisOneZeroTensor::dimension; j++)
        {
          this->components[j] *= factor;
        }

        return *this;
      }
    }

    return *this;
  }

  // Two methods for accessing the elements (the second for preventing errors in case a const is expected)
  GaloisScalar& GaloisZeroOneTensor::operator[](const unsigned int& component_index)
  {
    return this->components[component_index];
  }

  const GaloisScalar& GaloisZeroOneTensor::operator[](const unsigned int& component_index) const
  {
    return this->components[component_index];
  }

  // Get the dimension
  unsigned int GaloisZeroOneTensor::get_dimension() const
  {
    return GaloisOneZeroTensor::dimension;
  }

  // Binary operators for assigning or reading the values via iostream
  std::istream& operator>>(std::istream& is, GaloisZeroOneTensor& covector_input)
  {
    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      is >> covector_input.components[i];
    }

    return is;
  }

  std::ostream& operator<<(std::ostream& os, const GaloisZeroOneTensor& covector_output)
  {
    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      os << covector_output.components[i];
      os << "\n";
    }
      os << "\n";

    return os;
  }

  // Vector and Covector Product
  GaloisOneZeroTensor vector_product(std::vector<GaloisZeroOneTensor>& factors)
  {
    GaloisOneZeroTensor result;

    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      std::vector<std::vector<unsigned int> > sub_matrix;
      sub_matrix.resize(GaloisOneZeroTensor::dimension - 1);
      unsigned int j_sub{0};
      for(unsigned int j=0; j<GaloisOneZeroTensor::dimension - 1; j++)
      {
        sub_matrix[j].resize(factors.size());
        if(j==i)
        {
          j_sub++;
        }
        
        for(unsigned int k=0; k<factors.size(); k++)
        {
          std::cout << sub_matrix[j][k] << std::endl;
          sub_matrix[j][k] = factors[k][j_sub];
        }
        j_sub++;
      }

      if(i%2 == 0)
      {
        result[i] = +determinant_of_matrix(sub_matrix);
      }
      else
      {
        result[i] = -determinant_of_matrix(sub_matrix);
      }
    }

    return result;
  }

  GaloisZeroOneTensor covector_product(std::vector<GaloisOneZeroTensor>& factors)
  {
    GaloisZeroOneTensor result;

    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      std::vector<std::vector<unsigned int> > sub_matrix;
      sub_matrix.resize(GaloisOneZeroTensor::dimension - 1);
      unsigned int j_sub{0};
      for(unsigned int j=0; j<GaloisOneZeroTensor::dimension - 1; j++)
      {
        sub_matrix[j].resize(factors.size());
        if(j==i)
        {
          j_sub++;
        }
        
        for(unsigned int k=0; k<factors.size(); k++)
        {
          sub_matrix[j][k] = factors[k][j_sub];
        }
        j_sub++;
      }

      if(i%2 == 0)
      {
        result[i] = +determinant_of_matrix(sub_matrix);
      }
      else
      {
        result[i] = -determinant_of_matrix(sub_matrix);
      }
    }

    return result;
  }

  // Scalar product
  GaloisScalar operator*(GaloisOneZeroTensor lhs, GaloisZeroOneTensor rhs)
  {
    return rhs(lhs);
  }

  GaloisScalar operator*(GaloisZeroOneTensor lhs, GaloisOneZeroTensor rhs)
  {
    return lhs(rhs);
  }

}

#endif
