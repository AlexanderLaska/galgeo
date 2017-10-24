#ifndef GALOISONEONETENSOR_CPP
#define GALOISONEONETENSOR_CPP

#include "../include/galgeo/GaloisOneOneTensor.hpp"


namespace galgeo
{

  // Initialize the dimension is inherited from the GaloisOneZeroTensor class
  unsigned int GaloisOneOneTensor::dimension = GaloisOneZeroTensor::dimension;

  // Standard Constructor
  GaloisOneOneTensor::GaloisOneOneTensor()
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
  GaloisOneOneTensor::GaloisOneOneTensor(const GaloisOneOneTensor& rhs) : components(rhs.components) {}

  // (Virtual) Destructor
  GaloisOneOneTensor::~GaloisOneOneTensor() {}

  // The Operator for assigning another galois one_one_tensor to it
  GaloisOneOneTensor& GaloisOneOneTensor::operator=(const GaloisOneOneTensor& rhs)
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
  bool GaloisOneOneTensor::operator==(const GaloisOneOneTensor& rhs)
  {
    GaloisOneOneTensor dummyRhs{rhs};
    dummyRhs.norm();
    this->norm();
    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
          if(this->components[i][j] != dummyRhs.components[i][j])
          {
            return false;
          }
      }
    }

    return true;
  }

  bool GaloisOneOneTensor::operator==(const std::vector<std::vector<int> >& rhs)
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

  bool GaloisOneOneTensor::operator!=(const GaloisOneOneTensor& rhs)
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

  bool GaloisOneOneTensor::operator!=(const std::vector<std::vector<int> >& rhs)
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

  // Addition of two GaloisOneOneTensors 
  GaloisOneOneTensor GaloisOneOneTensor::operator+(const GaloisOneOneTensor& rhs)
  {
    GaloisOneOneTensor result;

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
  GaloisOneOneTensor& GaloisOneOneTensor::operator+=(const GaloisOneOneTensor& rhs)
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
  GaloisOneOneTensor GaloisOneOneTensor::operator-(const GaloisOneOneTensor& rhs)
  {
    GaloisOneOneTensor result;

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
  GaloisOneOneTensor& GaloisOneOneTensor::operator-=(const GaloisOneOneTensor& rhs)
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
  GaloisOneOneTensor GaloisOneOneTensor::operator*(const GaloisOneOneTensor& rhs)
  {
    GaloisOneOneTensor result;

    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
        for(unsigned int k=0; k<dimension; k++)
        {
          result[i][j] += this->components[i][k] * rhs[j][k];
          /*
          if( rhs[0][0] == 1 && rhs[0][1] == 2 && rhs[0][2] == 2 ) {
            std::cout << this->components[i][k] << "*" << rhs[j][k] << " + ";
          }
          */
        }
        /*
        if( rhs[0][0] == 1 && rhs[0][1] == 2 && rhs[0][2] == 2 ) {
          std::cout << " = " << result[i][j] << "\t\t";
        }
        */
      }
      /*
      if( rhs[0][0] == 1 && rhs[0][1] == 2 && rhs[0][2] == 2 ) {
        std::cout << "\n";
      }
      */
    }

    return result;
  }

  // Executing a one_one_tensor multiplication and assigning it to this one_one_tensor
  GaloisOneOneTensor& GaloisOneOneTensor::operator*=(const GaloisOneOneTensor& rhs)
  {
    GaloisOneOneTensor result = (*this) * rhs;
    (*this) = result;
    return *this;
  }

  // Multiplying the one_one_tensor from the left with a GaloisOneZeroTensor
  GaloisOneZeroTensor GaloisOneOneTensor::operator*(const GaloisOneZeroTensor& rhs)
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
  GaloisOneOneTensor GaloisOneOneTensor::operator*(const GaloisScalar& rhs)
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

  GaloisOneOneTensor& GaloisOneOneTensor::operator*=(const GaloisScalar& rhs)
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

  // Randomize the values uniformely
  GaloisOneOneTensor& GaloisOneOneTensor::randomize(const int lower_limit, const int upper_limit)
  {
    for(std::vector<GaloisScalar>& row : this->components)
    {
      for(GaloisScalar& component : row)
      {
        component.randomize(lower_limit, upper_limit);
      }
    }

    return *this;
  }

  // Norm the GaloisOneOneTensor (allowed by the homogenity)
  GaloisOneOneTensor& GaloisOneOneTensor::norm()
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
  GaloisScalar& GaloisOneOneTensor::operator()(const unsigned int& row_index, const unsigned int& column_index)
  {
    return this->components[row_index][column_index];
  }

  const GaloisScalar& GaloisOneOneTensor::operator()(const unsigned int& row_index, const unsigned int& column_index) const
  {
    return this->components[row_index][column_index];
  }

  std::vector<GaloisScalar>& GaloisOneOneTensor::operator[](const unsigned int& row_index)
  {
    return this->components[row_index];
  }

  const std::vector<GaloisScalar>& GaloisOneOneTensor::operator[](const unsigned int& row_index) const
  {
    return this->components[row_index];
  }

  // Access the value of the dimension
  unsigned int GaloisOneOneTensor::get_dimension() const
  {
    return this->dimension;
  }

  // Transpose the GaloisOneOneTensor
  GaloisOneOneTensor& GaloisOneOneTensor::transpose()
  {
    GaloisOneOneTensor result(*this);
    
    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
        this->components[i][j] = result(j,i);
      }
    }

    return *this;
  }

  GaloisOneOneTensor GaloisOneOneTensor::get_transposed()
  {
    GaloisOneOneTensor result{*this};
    
    for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
      {
        result[j][i] = this->components[i][j];
      }
    }

    return result;
  }

  // Determine the determinant of the GaloisOneOneTensor
  GaloisScalar GaloisOneOneTensor::determinant()
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

    GaloisScalar result{determinant_of_matrix(copy)};

    return result;
  }

  // Function for the determinant of a matrix
  int determinant_of_matrix(const std::vector<std::vector<unsigned int> >& matrix)
  {
    if(matrix.size() != 1)
    {
      int result{0};

      /*
      std::cout << "---" << std::endl;
      for(unsigned int j=0; j<matrix.size(); j++)
      {
        for(unsigned int k=0; k<matrix[0].size(); k++)
        {
          
          std::cout << matrix[j][k] << " ";
        }
      std::cout << std::endl;
      }
      */

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
            sub_matrix[j][k] = matrix[j_sub][k + 1]; 
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

      return ((result % GaloisScalar::prime) + GaloisScalar::prime) % GaloisScalar::prime;
    }
    else
    {
      return ((matrix[0][0] % GaloisScalar::prime) + GaloisScalar::prime) % GaloisScalar::prime;
    }
  }

  // Test whether the matrix is singular
  bool GaloisOneOneTensor::is_invertible()
  {
    return this->determinant() != 0;
  }

  // Invert the one_one_tensor if possible
  GaloisOneOneTensor& GaloisOneOneTensor::invert()
  {
    this->components = (this->get_inverse()).components;

    return *this;
  }

  GaloisOneOneTensor GaloisOneOneTensor::get_inverse()
  {
    GaloisScalar determinant;
    determinant = this->determinant();

    if(determinant != 0)
    {
      GaloisOneOneTensor result;

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
            result[j][i] = +determinant_of_matrix(minor);
          }
          else
          {
            result[j][i] = -determinant_of_matrix(minor);
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

  // Binary operators for assigning or reading the values via iostream
  std::istream& operator>>(std::istream& is, GaloisOneOneTensor& one_one_tensor_input)
  {
    for(unsigned int i=0; i<GaloisOneOneTensor::dimension; i++)
    {
      for(unsigned int j=0; j<GaloisOneOneTensor::dimension; j++)
      {
        is >> one_one_tensor_input.components[i][j];
        one_one_tensor_input.components[i][j] = ((one_one_tensor_input.components[i][j] % GaloisScalar::prime) + GaloisScalar::prime) % GaloisScalar::prime;
      }
    }

    return is;
  }

  std::ostream& operator<<(std::ostream& os, const GaloisOneOneTensor& one_one_tensor_output)
  {
    for(unsigned int i=0; i<GaloisOneOneTensor::dimension; i++)
    {
      for(unsigned int j=0; j<GaloisOneOneTensor::dimension; j++)
      {
        os << one_one_tensor_output.components[i][j];
        std::cout << " ";
        //std::cout << i << " " << j << std::endl;
      }
      std::cout << "\n";
    }

    return os;
  }

}

#endif
