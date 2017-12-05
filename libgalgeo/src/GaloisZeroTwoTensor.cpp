#ifndef GALOISZEROTWOTENSOR_CPP
#define GALOISZEROTWOTENSOR_CPP

#include "../include/galgeo/GaloisZeroTwoTensor.hpp"


namespace galgeo
{

  // Initialize the value of the prime number by the same value as the prime number of the GaloisScalar class
  // ... and the dimension is inherited from the GaloisOneZeroTensor class
  int GaloisZeroTwoTensor::prime = GaloisScalar::prime;
  unsigned int GaloisZeroTwoTensor::dimension = GaloisOneZeroTensor::dimension;

  // Standard Constructor
  GaloisZeroTwoTensor::GaloisZeroTwoTensor()
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

  // Copy Constructors
  GaloisZeroTwoTensor::GaloisZeroTwoTensor(const GaloisZeroTwoTensor& rhs) : components(rhs.components) {}
  GaloisZeroTwoTensor::GaloisZeroTwoTensor(const std::vector<std::vector<GaloisScalar> >& rhs) : components(rhs) {}
  GaloisZeroTwoTensor::GaloisZeroTwoTensor(const std::vector<std::vector<int> >& rhs)
  {
    components.resize(dimension);
    for(unsigned int i=0; i<dimension; i++)
    {
      components[i].resize(dimension);
      for(unsigned int j=0; j<dimension; j++)
      {
        GaloisScalar component{rhs[i][j]};
        components[i][j] = component;
      }
    }
  }

  // (Virtual) Destructor
  GaloisZeroTwoTensor::~GaloisZeroTwoTensor() {}

  // The Operators for assigning another galois one_one_tensor to it
  GaloisZeroTwoTensor& GaloisZeroTwoTensor::operator=(const GaloisZeroTwoTensor& rhs)
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

  GaloisZeroTwoTensor& GaloisZeroTwoTensor::operator=(const std::vector<std::vector<GaloisScalar> >& rhs)
  {
    for(unsigned int i=0; i<dimension; i++)
    {
      components[i].resize(dimension);
    }

    // ...and fill it with the values of the assigned one_one_tensor
    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
        components[i][j] = rhs[i][j];
      }
    }

    return *this;
  }

  GaloisZeroTwoTensor& GaloisZeroTwoTensor::operator=(const std::vector<std::vector<int > >& rhs)
  {
    for(unsigned int i=0; i<dimension; i++)
    {
      components[i].resize(dimension);
    }

    // ...and fill it with the values of the assigned one_one_tensor
    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
        components[i][j] = rhs[i][j];
      }
    }

    return *this;
  }

  // Binary operators for checking (in-)equality
  bool GaloisZeroTwoTensor::operator==(const GaloisZeroTwoTensor& rhs)
  {
    this->norm();
    GaloisZeroTwoTensor _rhs{rhs};
    _rhs.norm();

    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
        if(this->components[i][j] != _rhs.components[i][j])
        {
          return false;
        }
      }
    }

    return true;
  }

  bool GaloisZeroTwoTensor::operator==(const std::vector<std::vector<int> >& rhs)
  {
    this->norm();
    GaloisZeroTwoTensor _rhs{rhs};
    _rhs.norm();

    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
        if(_rhs[i][j] != this->components[i][j])
        {
          return false;
        }
      }
    }

    return true;
  }

  bool GaloisZeroTwoTensor::operator!=(const GaloisZeroTwoTensor& rhs)
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

  bool GaloisZeroTwoTensor::operator!=(const std::vector<std::vector<int> >& rhs)
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

  // Addition of two GaloisZeroTwoTensors 
  GaloisZeroTwoTensor GaloisZeroTwoTensor::operator+(const GaloisZeroTwoTensor& rhs)
  {
    GaloisZeroTwoTensor result;

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
  GaloisZeroTwoTensor& GaloisZeroTwoTensor::operator+=(const GaloisZeroTwoTensor& rhs)
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
  GaloisZeroTwoTensor GaloisZeroTwoTensor::operator-(const GaloisZeroTwoTensor& rhs)
  {
    GaloisZeroTwoTensor result;

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
  GaloisZeroTwoTensor& GaloisZeroTwoTensor::operator-=(const GaloisZeroTwoTensor& rhs)
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

  // Multiplying the one_one_tensor from the right with a GaloisOneZeroTensor
  GaloisZeroOneTensor GaloisZeroTwoTensor::operator*(const GaloisOneZeroTensor& rhs)
  {
    GaloisZeroOneTensor result;

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
  GaloisZeroTwoTensor GaloisZeroTwoTensor::operator*(const GaloisScalar& rhs)
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

  GaloisZeroTwoTensor& GaloisZeroTwoTensor::operator*=(const GaloisScalar& rhs)
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

  // Norm the GaloisZeroTwoTensor (allowed by the homogenity)
  GaloisZeroTwoTensor& GaloisZeroTwoTensor::norm()
  {
    for(int i=dimension - 1; i >= 0; i--)
    {
      if(this->components[i][dimension - 1] != 0)
      {
        GaloisScalar factor{(this->components[i][dimension - 1]).get_inverse()};
        for(unsigned int j=0; j<dimension; j++)
        {
          for(unsigned int k=0; k<dimension; k++)
          {
            this->components[j][k] *= factor;
          }
        }

        return *this;
      }
    }

    return *this;
  }

  // Two methods for the accessing the elements (the second for preventing errors in case a const is expected)
  GaloisScalar& GaloisZeroTwoTensor::operator()(const unsigned int& row_index, const unsigned int& column_index)
  {
    return this->components[row_index][column_index];
  }

  const GaloisScalar& GaloisZeroTwoTensor::operator()(const unsigned int& row_index, const unsigned int& column_index) const
  {
    return this->components[row_index][column_index];
  }

  std::vector<GaloisScalar>& GaloisZeroTwoTensor::operator[](const unsigned int& row_index)
  {
    return this->components[row_index];
  }

  const std::vector<GaloisScalar>& GaloisZeroTwoTensor::operator[](const unsigned int& row_index) const
  {
    return this->components[row_index];
  }

  // Access the value of the dimension
  unsigned int GaloisZeroTwoTensor::get_dimension() const
  {
    return this->dimension;
  }

  // Transpose the GaloisZeroTwoTensor
  GaloisZeroTwoTensor& GaloisZeroTwoTensor::transpose()
  {
    GaloisZeroTwoTensor result(*this);
    
    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
        this->components[i][j] = result(j,i);
      }
    }

    return *this;
  }

  GaloisZeroTwoTensor GaloisZeroTwoTensor::get_transposed()
  {
    GaloisZeroTwoTensor _result(*this);
    
    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
        this->components[i][j] = _result(j,i);
      }
    }

    GaloisZeroTwoTensor result(*this);

    return result;
  }

  // Determine the determinant of the GaloisZeroTwoTensor
  GaloisScalar GaloisZeroTwoTensor::determinant()
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
            sub_matrix[j][k] = matrix[j_sub][k + 1] % GaloisZeroTwoTensor::prime;
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
  bool GaloisZeroTwoTensor::is_invertible()
  {
    return this->determinant() != 0;
  }

  // Invert the one_one_tensor if possible
  GaloisZeroTwoTensor& GaloisZeroTwoTensor::invert()
  {
    this->components = (this->get_inverse()).components;

    return *this;
  }

  GaloisZeroTwoTensor GaloisZeroTwoTensor::get_inverse()
  {
    GaloisScalar determinant;
    determinant = this->determinant();

    if(determinant != 0)
    {
      GaloisZeroTwoTensor result;
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
  std::istream& operator>>(std::istream& is, GaloisZeroTwoTensor& one_one_tensor_input)
  {
    for(unsigned int i=0; i<GaloisZeroTwoTensor::dimension; i++)
    {
      for(unsigned int j=0; j<GaloisZeroTwoTensor::dimension; j++)
      {
        is >> one_one_tensor_input.components[i][j];
        one_one_tensor_input.components[i][j] = ((one_one_tensor_input.components[i][j] % GaloisZeroTwoTensor::prime) + GaloisZeroTwoTensor::prime) % GaloisZeroTwoTensor::prime;
      }
    }

    return is;
  }

  std::ostream& operator<<(std::ostream& os, const GaloisZeroTwoTensor& zero_two_tensor_output)
  {
    for(unsigned int i=0; i<GaloisZeroTwoTensor::dimension; i++)
    {
      for(unsigned int j=0; j<GaloisZeroTwoTensor::dimension; j++)
      {
        os << zero_two_tensor_output.components[i][j];
        os << " ";
      }
      os << "\n";
    }

    return os;
  }

  // Multiplying the GaloisZeroTwoTensor form the left with a GaloisOneOneTensor
  GaloisZeroTwoTensor operator*(const GaloisOneOneTensor lhs, const GaloisZeroTwoTensor rhs)
  {
    GaloisZeroTwoTensor result;

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

  // Multiplying the GaloisZeroTwoTensor form the right with a GaloisOneOneTensor
  GaloisZeroTwoTensor operator*(const GaloisZeroTwoTensor lhs, const GaloisOneOneTensor rhs)
  {
    GaloisZeroTwoTensor result;

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
