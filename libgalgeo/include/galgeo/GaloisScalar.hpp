#ifndef GALOISSCALAR_HPP
#define GALOISSCALAR_HPP

#include <iostream>
#include <math.h>
#include <random>
#include "DivideByZeroException.hpp"


namespace galgeo
{

  // The GaloisScalar class models the algebraic structure of a Galois field with respect to a given prime number
  class GaloisScalar 
  {
    private:
      unsigned int value;

    public:
      static int prime; // 'int' (instead of 'unsigned int') in order to make '%' work properly

      // Constructors and destructor
      GaloisScalar();
      GaloisScalar(int _value);
      GaloisScalar(unsigned int _value);
      GaloisScalar(const GaloisScalar& rhs);
      virtual ~GaloisScalar();

      // Assigning a value via another GaloisScalar or an integer
      GaloisScalar& operator=(const GaloisScalar& rhs);
      GaloisScalar& operator=(const int& rhs);

      // Pre- and postfix increment operators
      GaloisScalar& operator++();
      GaloisScalar operator++(int dummy);
      GaloisScalar& operator--();
      GaloisScalar operator--(int dummy);


      // Binary logic operators for checking (in-)equality of GaloisScalars and GaloisScalars or integers
      template <typename T1, typename T2>
        friend bool operator==(const T1& lhs, const T2& rhs);
      template <typename T1, typename T2>
        friend bool operator!=(const T1& lhs, const T2& rhs);

      // Binary field operators
      template <typename T1, typename T2>
        friend GaloisScalar operator+(const T1& lhs, const T2& rhs);
      template <typename T>
        GaloisScalar& operator +=(const T& rhs);

      template <typename T1, typename T2>
        friend GaloisScalar operator-(const T1& lhs, const T2& rhs);
      template <typename T>
        GaloisScalar& operator -=(const T& rhs);
      
      template <typename T1, typename T2>
        friend GaloisScalar operator*(const T1& lhs, const T2& rhs);
      template <typename T>
        GaloisScalar& operator *=(const T& rhs);

      template <typename T1, typename T2>
        friend GaloisScalar operator/(const T1& lhs, const T2& rhs);
      template <typename T>
        GaloisScalar& operator /=(const T& rhs);


      // The power function with GaloisScalars as arguments (or mixed with integers)
      friend GaloisScalar pow(const GaloisScalar& base, const GaloisScalar& _exponent);
      friend GaloisScalar pow(const GaloisScalar& base, const unsigned int& _exponent);
      friend GaloisScalar pow(const unsigned int& base, const GaloisScalar& _exponent);
    

      // Randomize the value uniformely
      GaloisScalar& randomize(const int lower_limit = 0, const int upper_limit = GaloisScalar::prime - 1);

      // Yielding the value of the inverse as a new GaloisScalar;
      GaloisScalar get_inverse();

      // Inverting the GaloisScalar in situ
      GaloisScalar& invert();

      // Yielding the value of the nth root as a new GaloisScalar;
      template <typename T>
        GaloisScalar get_nth_root(T n);

      // Inverting the GaloisScalar in situ
      template <typename T>
        GaloisScalar& take_nth_root(T n);

      // Implicitely accessing the value
      operator int();

      // Accessing the value of the prime number
      unsigned int get_prime() const;
      friend void set_prime(const int& new_prime);

      // Binary operators for assigning or reading the value via iostream
      friend std::istream& operator>>(std::istream& is, GaloisScalar& scalar_input);
      friend std::ostream& operator<<(std::ostream& os, const GaloisScalar& scalar_output);
  };


  // Binary logic operators for checking (in-)equality of GaloisScalars and GaloisScalars or integers (template definition followed by explicit instantiations)
  template <typename T1, typename T2>
    bool operator==(const T1& lhs, const T2& rhs)
    {
      GaloisScalar _rhs{rhs};
      GaloisScalar _lhs{lhs};

      return (_lhs.value == _rhs.value);
    }

  template bool operator==(const GaloisScalar& lhs, const GaloisScalar& rhs);
  template bool operator==(const GaloisScalar& lhs, const int& rhs);
  template bool operator==(const int& lhs, const GaloisScalar& rhs);

  template <typename T1, typename T2>
    bool operator!=(const T1& lhs, const T2& rhs)
    {
      return !(lhs == rhs);
    }

  template bool operator!=(const GaloisScalar& lhs, const GaloisScalar& rhs);
  template bool operator!=(const GaloisScalar& lhs, const int& rhs);
  template bool operator!=(const int& lhs, const GaloisScalar& rhs);


  // Binary field operators (template definition followed by explicit instantiations)
  template <typename T1, typename T2>
    GaloisScalar operator+(const T1& lhs, const T2& rhs)
    {
      GaloisScalar _rhs{rhs};
      GaloisScalar _lhs{lhs};

      return GaloisScalar{static_cast<int>(_lhs.value) + static_cast<int>(_rhs.value)};
    }

  template GaloisScalar operator+(const GaloisScalar& lhs, const GaloisScalar& rhs);
  template GaloisScalar operator+(const GaloisScalar& lhs, const int& rhs);
  template GaloisScalar operator+(const int& lhs, const GaloisScalar& rhs);

  template <typename T1, typename T2>
    GaloisScalar operator-(const T1& lhs, const T2& rhs)
    {
      GaloisScalar _lhs{lhs};
      GaloisScalar _rhs{rhs};

      return GaloisScalar{static_cast<int>(_lhs.value) - static_cast<int>(_rhs.value)};
    }

  template GaloisScalar operator-(const GaloisScalar& lhs, const GaloisScalar& rhs);
  template GaloisScalar operator-(const GaloisScalar& lhs, const int& rhs);
  template GaloisScalar operator-(const int& lhs, const GaloisScalar& rhs);

  template <typename T1, typename T2>
    GaloisScalar operator*(const T1& lhs, const T2& rhs)
    {
      GaloisScalar _rhs{rhs};
      GaloisScalar _lhs{lhs};

      return GaloisScalar{_lhs.value * _rhs.value};
    }

  template GaloisScalar operator*(const GaloisScalar& lhs, const GaloisScalar& rhs);
  template GaloisScalar operator*(const GaloisScalar& lhs, const int& rhs);
  template GaloisScalar operator*(const int& lhs, const GaloisScalar& rhs);

  template <typename T1, typename T2>
    GaloisScalar operator/(const T1& lhs, const T2& rhs)
    {
      GaloisScalar _lhs{lhs};
      GaloisScalar _rhs{rhs};
      _rhs.invert();

      return GaloisScalar{_lhs.value * _rhs.value};
    }

  template GaloisScalar operator/(const GaloisScalar& lhs, const GaloisScalar& rhs);
  template GaloisScalar operator/(const GaloisScalar& lhs, const int& rhs);
  template GaloisScalar operator/(const int& lhs, const GaloisScalar& rhs);


  // The power function with GaloisScalars as arguments (or mixed with integers)
  GaloisScalar pow(const GaloisScalar& base, const GaloisScalar& _exponent);
  GaloisScalar pow(const GaloisScalar& base, const unsigned int& _exponent);
  GaloisScalar pow(const unsigned int& base, const GaloisScalar& _exponent);

  // Set the prime number
  void set_prime(const int& new_prime);

}

#endif
