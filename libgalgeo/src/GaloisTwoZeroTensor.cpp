#ifndef GALOISTWOZEROTENSOR_CPP
#define GALOISTWOZEROTENSOR_CPP

#include "../include/galgeo/GaloisTwoZeroTensor.hpp"


namespace galgeo
{

  // Initialize the value of the prime number by the same value as the prime number of the GaloisScalar class
  // ... and the dimension is inherited from the GaloisOneZeroTensor class
  int GaloisTwoZeroTensor::prime = GaloisScalar::prime;
  unsigned int GaloisTwoZeroTensor::dimension = GaloisOneZeroTensor::dimension;

  // Standard Constructor
  GaloisTwoZeroTensor::GaloisTwoZeroTensor()
  {
    components.resize(dimension);
    for(unsigned int i=0; i<dimension; i++)
    {
      components[i].resize(dimension);
      for(unsigned int j=0; j<dimension; j++)
      {
        components[i][j] = 0;
      }
    }
  }

  // Copy Constructor
  GaloisTwoZeroTensor::GaloisTwoZeroTensor(const GaloisTwoZeroTensor& rhs) : components(rhs.components) {}

  // (Virtual) Destructor
  GaloisTwoZeroTensor::~GaloisTwoZeroTensor() {}

  // The Operator for assigning another galois one_one_tensor to it
  GaloisTwoZeroTensor& GaloisTwoZeroTensor::operator=(const GaloisTwoZeroTensor& rhs)
  {
    // Test whether the assigned is the current one_one_tensor
    if(&rhs == this)
    {
      return *this;
    }

    for(unsigned int i=0; i<dimension; i++)
    {
      components[i].resize(dimension);
    }

    // ...and fill it with the values of the assigned one_one_tensor
    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
        components[i][j] = rhs(i,j);
      }
    }

    return *this;
  }

  // Binary operators for checking (in-)equality
  bool GaloisTwoZeroTensor::operator==(const GaloisTwoZeroTensor& rhs)
  {
    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
          if(this->components[i][j] == rhs.components[i][j])
          {
            return false;
          }
      }
    }

    return true;
  }

  bool GaloisTwoZeroTensor::operator==(const std::vector<std::vector<int> >& rhs)
  {
    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
          if(rhs[i][j] != this->components[i][j])
          {
            return false;
          }
      }
    }

    return true;
  }

  bool GaloisTwoZeroTensor::operator!=(const GaloisTwoZeroTensor& rhs)
  {
    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
          if(this->components[i][j] != rhs.components[i][j])
          {
            return true;
          }
      }
    }

    return false;
  }

  bool GaloisTwoZeroTensor::operator!=(const std::vector<std::vector<int> >& rhs)
  {
    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
          if(rhs[i][j] != this->components[i][j])
          {
            return true;
          }
      }
    }

    return false;
  }

  // Addition of two GaloisTwoZeroTensors 
  GaloisTwoZeroTensor GaloisTwoZeroTensor::operator+(const GaloisTwoZeroTensor& rhs)
  {
    GaloisTwoZeroTensor result;

    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
        result[i][j] = this->components[i][j] + rhs[i][j];
      }
    }

    return result;
  }

  // Adding and assigning a second one_one_tensor to this one_one_tensor
  GaloisTwoZeroTensor& GaloisTwoZeroTensor::operator+=(const GaloisTwoZeroTensor& rhs)
  {
    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
        this->components[i][j] += rhs[i][j];
      }
    }

    return *this;
  }

  // Subtraction of two matrices
  GaloisTwoZeroTensor GaloisTwoZeroTensor::operator-(const GaloisTwoZeroTensor& rhs)
  {
    GaloisTwoZeroTensor result;

    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
        result[i][j] = this->components[i][j] - rhs[i][j];
      }
    }

    return result;
  }

  // Subtracting and assigning a second one_one_tensor to this one_one_tensor
  GaloisTwoZeroTensor& GaloisTwoZeroTensor::operator-=(const GaloisTwoZeroTensor& rhs)
  {
    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
        this->components[i][j] -= rhs[i][j];
      }
    }

    return *this;
  }

  // Matrix Multplicating this one_one_tensor from the left
  GaloisTwoZeroTensor GaloisTwoZeroTensor::operator*(const GaloisOneOneTensor& rhs)
  {
    GaloisTwoZeroTensor result;

    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
        for(unsigned int k=0; k<dimension; k++)
        {
          result[i][j] += this->components[i][k] * rhs[k][j];
        }
      }
    }

    return result;
  }

  // Executing a one_one_tensor multiplication and assigning it to this one_one_tensor
  GaloisTwoZeroTensor& GaloisTwoZeroTensor::operator*=(const GaloisOneOneTensor& rhs)
  {
    GaloisTwoZeroTensor result = (*this) * rhs;
    (*this) = result;
    return *this;
  }

  // Multiplying the two_zero_one_tensor from the right with a GaloisZeroOneTensor
  GaloisOneZeroTensor GaloisTwoZeroTensor::operator*(const GaloisZeroOneTensor& rhs)
  {
    GaloisOneZeroTensor result;

    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
            result[i] += this->components[i][j] * rhs[j];
      }
    }

    return result;
  }

  // Multiply the one_one_tensor with an integer
  GaloisTwoZeroTensor GaloisTwoZeroTensor::operator*(const GaloisScalar& rhs)
  {
    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
        this->components[i][j] = this->components[i][j] * rhs;
      }
    }

    return *this;
  }

  GaloisTwoZeroTensor& GaloisTwoZeroTensor::operator*=(const GaloisScalar& rhs)
  {
    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
        this->components[i][j] *= rhs;
      }
    }

    return *this;
  }

  // Norm the GaloisTwoZeroTensor (allowed by the homogenity)
  GaloisTwoZeroTensor& GaloisTwoZeroTensor::norm()
  {
    for(int i=dimension - 1; i >= 0; i--)
    {
      if(this->components[dimension - 1][i] != 0)
      {
        GaloisScalar factor((this->components[dimension - 1][i]).get_inverse());
        for(unsigned int j=0; j<dimension; j++)
        {
          for(unsigned int k=0; k<dimension; k++)
          { this->components[j][k] *= factor;
          }
        }

        return *this;
      }
    }

    return *this;
  }

  // Two methods for the accessing the elements (the second for preventing errors in case a const is expected)
  GaloisScalar& GaloisTwoZeroTensor::operator()(const unsigned int& row_index, const unsigned int& column_index)
  {
    return this->components[row_index][column_index];
  }

  const GaloisScalar& GaloisTwoZeroTensor::operator()(const unsigned int& row_index, const unsigned int& column_index) const
  {
    return this->components[row_index][column_index];
  }

  std::vector<GaloisScalar>& GaloisTwoZeroTensor::operator[](const unsigned int& row_index)
  {
    return this->components[row_index];
  }

  const std::vector<GaloisScalar>& GaloisTwoZeroTensor::operator[](const unsigned int& row_index) const
  {
    return this->components[row_index];
  }

  // Access the value of the dimension
  unsigned int GaloisTwoZeroTensor::get_dimension() const
  {
    return this->dimension;
  }

  // Transpose the GaloisTwoZeroTensor
  GaloisTwoZeroTensor& GaloisTwoZeroTensor::transpose()
  {
    GaloisTwoZeroTensor result(*this);
    
    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
        this->components[i][j] = result(j,i);
      }
    }

    return *this;
  }

  GaloisTwoZeroTensor GaloisTwoZeroTensor::get_transposed()
  {
    GaloisTwoZeroTensor _result(*this);
    
    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
        this->components[i][j] = _result(j,i);
      }
    }

    GaloisTwoZeroTensor result(*this);

    return result;
  }

  // Determine the determinant of the GaloisTwoZeroTensor
  GaloisScalar GaloisTwoZeroTensor::determinant()
  {
    std::vector<std::vector<unsigned int> > copy;

    copy.resize(dimension);
    for(unsigned int i=0; i<dimension; i++)
    {
      copy[i].resize(dimension);
      for(unsigned int j=0; j<dimension; j++)
      {
        copy[i][j] = this->components[i][j];
      }
    }

    GaloisScalar result(determinant_of_matrix(copy));

    return result;
  }

  /*
  // Function for the determinant of a matrix
  int determinant_of_matrix(const std::vector<std::vector<unsigned int> >& matrix)
  {
    if(matrix.size() != 1)
    {
      int result{0};

      for(unsigned int i=0; i<matrix.size(); i++)
      {
        std::vector<std::vector<unsigned int> > sub_matrix;
        unsigned int j_sub{0};
        
        sub_matrix.resize(matrix.size() - 1);
        for(unsigned int j=0; j<sub_matrix.size(); j++)
        {
          sub_matrix[j].resize(matrix.size() - 1);
   
          if(j == i)
          {
            j_sub++;
          }
   
          for(unsigned int k=0; k<sub_matrix.size(); k++)
          {
            sub_matrix[j][k] = matrix[j_sub][k + 1] % GaloisTwoZeroTensor::prime;
          }
   
   
          j_sub++;
        }

        if(i%2 == 0)
        {
          result += determinant_of_matrix(sub_matrix)*matrix[i][0];
        }
        else
        {
          result -= determinant_of_matrix(sub_matrix)*matrix[i][0];
        }
      }

      return result;
    }
    else
    {
      return matrix[0][0];
    }
  }

  // Test whether the matrix is singular
  bool GaloisTwoZeroTensor::is_invertible()
  {
    return this->determinant() != 0;
  }

  // Invert the one_one_tensor if possible
  GaloisTwoZeroTensor& GaloisTwoZeroTensor::invert()
  {
    this->components = (this->get_inverse()).components;

    return *this;
  }

  GaloisTwoZeroTensor GaloisTwoZeroTensor::get_inverse()
  {
    GaloisScalar determinant;
    determinant = this->determinant();

    if(determinant != 0)
    {
      GaloisTwoZeroTensor result;
      GaloisScalar inverse_of_determinant;

      inverse_of_determinant = determinant.get_inverse();

      for(unsigned int i=0; i<dimension; i++)
      {
        for(unsigned int j=0; j<dimension; j++)
        {
          std::vector<std::vector<unsigned int> > minor;

          minor.resize(dimension - 1);

          unsigned int k_sub{0};
          for(unsigned int k=0; k<dimension - 1; k++)
          {
            if(k==i)
            {
              k_sub++;
            }

            minor[k].resize(dimension - 1);

            unsigned int l_sub{0};
            for(unsigned int l=0; l<dimension - 1; l++)
            {
              if(l==j)
              {
                l_sub++;
              }

              minor[k][l] = this->components[k_sub][l_sub];

              l_sub++;
            }

            k_sub++;
          }

          if((i+j)%2 == 0)
          {
            result[j][i] = (+1)*inverse_of_determinant*(determinant_of_matrix(minor));
          }
          else
          {
            result[j][i] = (-1)*inverse_of_determinant*(determinant_of_matrix(minor));
          }
        }
      }
        return result;
    }
    else
    {
      return *this;
    }
  }
  */

  // Binary operators for assigning or reading the values via iostream
  std::istream& operator>>(std::istream& is, GaloisTwoZeroTensor& one_one_tensor_input)
  {
    for(unsigned int i=0; i<GaloisTwoZeroTensor::dimension; i++)
    {
      for(unsigned int j=0; j<GaloisTwoZeroTensor::dimension; j++)
      {
        is >> one_one_tensor_input.components[i][j];
        one_one_tensor_input.components[i][j] = ((one_one_tensor_input.components[i][j] % GaloisTwoZeroTensor::prime) + GaloisTwoZeroTensor::prime) % GaloisTwoZeroTensor::prime;
      }
    }

    return is;
  }

  std::ostream& operator<<(std::ostream& os, const GaloisTwoZeroTensor& one_one_tensor_output)
  {
    for(unsigned int i=0; i<GaloisTwoZeroTensor::dimension; i++)
    {
      for(unsigned int j=0; j<GaloisTwoZeroTensor::dimension; j++)
      {
        os << one_one_tensor_output.components[i][j];
        std::cout << " ";
      }
      std::cout << "\n";
    }

    return os;
  }

  // Multiplying the GaloisTwoZeroTensor form the left with a GaloisOneOneTensor
  GaloisTwoZeroTensor operator*(GaloisOneOneTensor& lhs, GaloisTwoZeroTensor& rhs)
  {
    GaloisTwoZeroTensor result;

    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      for(unsigned int j=0; j<GaloisOneZeroTensor::dimension; j++)
      {
        for(unsigned int k=0; k<GaloisOneZeroTensor::dimension; k++)
        {
          result[i][j] += lhs[i][k] * rhs[k][j];
        }
      }
    }

    return result;
  }

  // Multiplying the GaloisTwoZeroTensor form the left with a GaloisOneOneTensor
  GaloisTwoZeroTensor operator*(GaloisTwoZeroTensor& lhs, GaloisOneOneTensor& rhs)
  {
    GaloisTwoZeroTensor result;

    for(unsigned int i=0; i<GaloisOneZeroTensor::dimension; i++)
    {
      for(unsigned int j=0; j<GaloisOneZeroTensor::dimension; j++)
      {
        for(unsigned int k=0; k<GaloisOneZeroTensor::dimension; k++)
        {
          result[i][j] += lhs[i][k] * rhs[k][j];
        }
      }
    }

    return result;
  }

}

#endif
