#ifndef GALOISSCALAR_CPP
#define GALOISSCALAR_CPP

#include "../include/galgeo/GaloisScalar.hpp"


namespace galgeo
{
  // Set the random device
  std::random_device rd;
  std::mt19937 gen(rd());

  // Sets value of the prime number for all instances of the GaloisScalar class 
  int GaloisScalar::prime = 2;

  // Standard Constructor
  GaloisScalar::GaloisScalar() : value(0) {}

  // Parameter Constructors
  GaloisScalar::GaloisScalar(int _value) : value(static_cast<unsigned int>(((_value % prime) + prime) % prime)) {}
  GaloisScalar::GaloisScalar(unsigned int _value) : value(((_value % prime) + prime) % prime) {}

  // Copy Constructor
  GaloisScalar::GaloisScalar(const GaloisScalar& original_GaloisScalar) : value(original_GaloisScalar.value) {};

  // (Virtual) Destructor
  GaloisScalar::~GaloisScalar() {}


  // Assigning a value via another GaloisScalar or an integer
  GaloisScalar& GaloisScalar::operator=(const GaloisScalar& rhs)
  { 
    // Test whether the assigned is the current GaloisScalar
    if(&rhs == this)
    {
      return *this;
    }

    this->value = rhs.value;

    return *this;
  }

  GaloisScalar& GaloisScalar::operator=(const int& rhs)
  { 
    value = ((rhs % prime) + prime) % prime;

    return *this;
  }


  // Unary prefix operators for adding one or subtracting one
  GaloisScalar& GaloisScalar::operator++()
  {
    *this += 1;

    return *this;
  }

  GaloisScalar GaloisScalar::operator++(int dummy)
  {
    GaloisScalar temp{*this};

    *this += 1;

    return temp;
  }

  GaloisScalar& GaloisScalar::operator--()
  {
    *this -= 1;

    return *this;
  }

  GaloisScalar GaloisScalar::operator--(int dummy)
  {
    GaloisScalar temp{*this};

    *this -= 1;

    return temp;
  }


  // Randomize the value uniformely
  GaloisScalar& GaloisScalar::randomize(const int lower_limit, const int upper_limit)
  {
    std::uniform_int_distribution<> dis(lower_limit, upper_limit);

    this->value = dis(gen);

    return *this;
  }

  // Yield the value of the inverse
  GaloisScalar GaloisScalar::get_inverse()
  {
    if(this->value == 0)
    {
      throw DivideByZeroException();
    }

    GaloisScalar inverse(*this);
    inverse.value = pow(*this, prime - 2);

    return inverse;
  }

  // Invert the GaloisScalar in situ
  GaloisScalar& GaloisScalar::invert()
  {
    if(this->value == 0)
    {
      throw DivideByZeroException();
    }

    this->value = pow(*this, prime - 2);

    return *this;
  }

  // Yielding the value of the nth root as a new GaloisScalar;
  template <typename T>
    GaloisScalar GaloisScalar::get_nth_root(T n)
    {
      GaloisScalar exponent{n};
      GaloisScalar nthRoot{0};
      GaloisScalar result{0};

      size_t i = 0;
      GaloisScalar base{0};

      while(nthRoot != this->value && i != prime + 1)
      {
        base = i;
        nthRoot = pow(base, exponent);

        ++i;
      }
      result = base;
      
      if(i == prime + 1)
      {
        }
      return result;
    }

  template GaloisScalar GaloisScalar::get_nth_root<const unsigned int>(const unsigned int n);
  template GaloisScalar GaloisScalar::get_nth_root<GaloisScalar>(GaloisScalar n);


  // Inverting the GaloisScalar in situ
  template <typename T>
    GaloisScalar& GaloisScalar::take_nth_root(T n)
    {
      this->value = this->get_nth_root(n);

      return *this;
    }

  template GaloisScalar& GaloisScalar::take_nth_root<const unsigned int>(const unsigned int n);
  template GaloisScalar& GaloisScalar::take_nth_root<GaloisScalar>(GaloisScalar n);


  // Implicitely accessing the value
  GaloisScalar::operator int()
  {
    return value;
  }


  // Get and set the prime number
  unsigned int GaloisScalar::get_prime() const
  {
    return this->prime;
  }

  void set_prime(const int& new_prime)
  {
    GaloisScalar::prime = new_prime;
  }


  // Binary operators for assigning or printing the value via iostream
  std::istream& operator>>(std::istream& is, GaloisScalar& scalar_input)
  {
    int read_in;
    is >> read_in;

    scalar_input = read_in;
    return is;
  }

  std::ostream& operator<<(std::ostream& os, const GaloisScalar& scalar_output)
  {
    os << scalar_output.value;

    return os;
  }

  // x= operators

  template <typename T>
    GaloisScalar& GaloisScalar::operator +=(const T& rhs)
    {
      GaloisScalar sum{this->value + rhs};
      *this = sum;

      return *this;
    }

  template GaloisScalar& GaloisScalar::operator+=(const GaloisScalar& rhs);
  template GaloisScalar& GaloisScalar::operator+=(const int& rhs);


  template <typename T>
    GaloisScalar& GaloisScalar::operator -=(const T& rhs)
    {
      GaloisScalar difference{this->value - rhs};
      *this = difference;

      return *this;
    }

  template GaloisScalar& GaloisScalar::operator-=(const GaloisScalar& rhs);
  template GaloisScalar& GaloisScalar::operator-=(const int& rhs);


  template <typename T>
    GaloisScalar& GaloisScalar::operator *=(const T& rhs)
    {
      GaloisScalar product{this->value * rhs};
      *this = product;

      return *this;
    }

  template GaloisScalar& GaloisScalar::operator*=(const GaloisScalar& rhs);
  template GaloisScalar& GaloisScalar::operator*=(const int& rhs);


  template <typename T>
    GaloisScalar& GaloisScalar::operator /=(const T& rhs)
    {
      GaloisScalar quotient{this->value / rhs};
      *this = quotient;

      return *this;
    }

  template GaloisScalar& GaloisScalar::operator/=(const GaloisScalar& rhs);
  template GaloisScalar& GaloisScalar::operator/=(const int& rhs);


  // Pow
  GaloisScalar pow(const GaloisScalar& _base, const GaloisScalar& _exponent)
  {
    unsigned int result_value{1};
    unsigned int remainder;
    unsigned int base{_base.value};
    unsigned int exponent{_exponent.value};

    while (exponent != 0)
    {
      remainder = exponent % 2;
      exponent = exponent / 2;

      if (remainder == 1)
      {
        result_value = (result_value * base) % GaloisScalar::prime;
      }

      base = (base * base) % GaloisScalar::prime;
    }

    return GaloisScalar{result_value};
  }

  GaloisScalar pow(const GaloisScalar& _base, const unsigned int& _exponent)
  {
    unsigned int result_value{1};
    unsigned int remainder;
    unsigned int base{_base.value};
    unsigned int exponent{_exponent};

    while (exponent != 0)
    {
      remainder = exponent % 2;
      exponent = exponent / 2;

      if (remainder == 1)
      {
        result_value = (result_value * base) % GaloisScalar::prime;
      }

      base = (base * base) % GaloisScalar::prime;
    }

    return GaloisScalar{result_value};
  }

  GaloisScalar pow(const unsigned int& _base, const GaloisScalar& _exponent)
  {
    unsigned int result_value{1};
    unsigned int remainder;
    unsigned int base{_base};
    unsigned int exponent{_exponent.value};

    while (exponent != 0)
    {
      remainder = exponent % 2;
      exponent = exponent / 2;

      if (remainder == 1)
      {
        result_value = (result_value * base) % GaloisScalar::prime;
      }

      base = (base * base) % GaloisScalar::prime;
    }

    return GaloisScalar{result_value};
  }

}

#endif
