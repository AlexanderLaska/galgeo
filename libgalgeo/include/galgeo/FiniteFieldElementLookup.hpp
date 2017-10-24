#ifndef GALOIS_FINITE_FIELD_ELEMENT_LOOKUP
#define GALOIS_FINITE_FIELD_ELEMENT_LOOKUP

#include <array>
#include <cmath>

namespace galgeo
{
  template <int prime>
  class FiniteFieldElementLookup
  {
  public:
    typedef long int value_t;

    //! Default constructor
    FiniteFieldElementLookup() : value(0) { }
    //! Constructor with value
    FiniteFieldElementLookup(value_t new_value) : value(new_value) { modulate(); }

    //! Get-accessor for value
    const value_t& get_value() { return value; }
    //! Set-accessor for the value
    void set_value(value_t new_value) { value = new_value; modulate(); }

    //! Equal Comparison operator
    bool operator==(const FiniteFieldElementLookup& rhs) const { return value == rhs.value; }
    //! Unequal Comparison operator
    bool operator!=(const FiniteFieldElementLookup& rhs) const { return !operator==(rhs); }

    FiniteFieldElementLookup& operator++()  { value++; if (value == prime) modulate(); return *this; }
    FiniteFieldElementLookup operator++(int) { FiniteFieldElementLookup tmp(*this); tmp.operator++(); return tmp; }

    FiniteFieldElementLookup& operator+=(const FiniteFieldElementLookup& rhs) { value = lookup_addition[prime*value + rhs.value]; return *this; }
    FiniteFieldElementLookup& operator-=(const FiniteFieldElementLookup& rhs) { value = lookup_addition[prime*value + prime - rhs.value]; return *this; }
    FiniteFieldElementLookup& operator*=(const FiniteFieldElementLookup& rhs) { value = lookup_multiplication[prime*value + rhs.value]; return *this; }
    FiniteFieldElementLookup& operator/=(const FiniteFieldElementLookup& rhs) { value = lookup_multiplication[prime*value + lookup_inverse(rhs.value)]; return *this; }

    FiniteFieldElementLookup& operator+=(const value_t& rhs) { value = lookup_addition[prime*value + rhs]; return *this; }
    FiniteFieldElementLookup& operator-=(const value_t& rhs) { value = lookup_addition[prime*value + prime - rhs]; return *this; }
    FiniteFieldElementLookup& operator*=(const value_t& rhs) { value = lookup_multiplication[prime*value + rhs]; return *this; }
    FiniteFieldElementLookup& operator/=(const value_t& rhs) { value = lookup_multiplication[prime*value + lookup_inverse(rhs)]; return *this; }

    FiniteFieldElementLookup inverse() { return FiniteFieldElementLookup(lookup_inverse[value]); }
    FiniteFieldElementLookup& invert() { value = lookup_inverse[value]; return *this; }

  private:
    //! Member with the actual value
    value_t value;

    //! Static lookup table storing the results of the addition
    static std::array<value_t, prime*prime> lookup_addition;
    //! Static lookup table storing the results of the multiplication
    static std::array<value_t, prime*prime> lookup_multiplication;
    //! Static lookup table storing the inverses of all values
    static std::array<value_t, prime> lookup_inverse;
    //! Static lookup table storing whether a value has a square root
    static std::array<bool, prime> lookup_has_square_root;
    //! Static lookup table for the square roots of the numbers (0 if there is no square root)
    static std::array<value_t, prime> lookup_square_root;

    //! Function to calculate the mod with respect to the prime
    void modulate() { value = value % prime; }

    // Friend functions
    template <class ostream, int fprime> friend ostream& operator<<(ostream& stream, const FiniteFieldElementLookup<fprime>& rhs);
    template <class istream, int fprime> friend istream& operator<<(istream& stream, const FiniteFieldElementLookup<fprime>& rhs);

    template <int fprime> friend FiniteFieldElementLookup<fprime> sqrt(FiniteFieldElementLookup<fprime>& number);
    template <int fprime> friend FiniteFieldElementLookup<fprime> pow(FiniteFieldElementLookup<fprime>& number, int exponent);
  };

  // Binary operators
  template <int prime> 
  FiniteFieldElementLookup<prime> operator+(const FiniteFieldElementLookup<prime>& lhs, const FiniteFieldElementLookup<prime>& rhs) 
  { return FiniteFieldElementLookup<prime>(lhs) += rhs; }
  template <int prime> 
  FiniteFieldElementLookup<prime> operator-(const FiniteFieldElementLookup<prime>& lhs, const FiniteFieldElementLookup<prime>& rhs) 
  { return FiniteFieldElementLookup<prime>(lhs) -= rhs; }
  template <int prime> 
  FiniteFieldElementLookup<prime> operator*(const FiniteFieldElementLookup<prime>& lhs, const FiniteFieldElementLookup<prime>& rhs) 
  { return FiniteFieldElementLookup<prime>(lhs) *= rhs; }
  template <int prime> 
  FiniteFieldElementLookup<prime> operator/(const FiniteFieldElementLookup<prime>& lhs, const FiniteFieldElementLookup<prime>& rhs) 
  { return FiniteFieldElementLookup<prime>(lhs) *= rhs; }

  template <int prime> 
  FiniteFieldElementLookup<prime> operator+(const FiniteFieldElementLookup<prime>& lhs, const typename FiniteFieldElementLookup<prime>::value_t& rhs) 
  { return FiniteFieldElementLookup<prime>(lhs) += rhs; }
  template <int prime> 
  FiniteFieldElementLookup<prime> operator-(const FiniteFieldElementLookup<prime>& lhs, const typename FiniteFieldElementLookup<prime>::value_t& rhs) 
  { return FiniteFieldElementLookup<prime>(lhs) -= rhs; }
  template <int prime> 
  FiniteFieldElementLookup<prime> operator*(const FiniteFieldElementLookup<prime>& lhs, const typename FiniteFieldElementLookup<prime>::value_t& rhs) 
  { return FiniteFieldElementLookup<prime>(lhs) *= rhs; }
  template <int prime> 
  FiniteFieldElementLookup<prime> operator/(const FiniteFieldElementLookup<prime>& lhs, const typename FiniteFieldElementLookup<prime>::value_t& rhs) 
  { return FiniteFieldElementLookup<prime>(lhs) *= rhs; }

  // Other functions
  template <int prime> 
  FiniteFieldElementLookup<prime> sqrt(FiniteFieldElementLookup<prime>& number)
  { return FiniteFieldElementLookup<prime>(FiniteFieldElementLookup<prime>::lookup_square_root(number)); }
  template <int prime> FiniteFieldElementLookup<prime> pow(FiniteFieldElementLookup<prime>& number, int exponent)
  { return FiniteFieldElementLookup<prime>(static_cast<int>(pow(number.value, exponent))); }

  // Initialise the lookup addition table
  template <int prime>
  std::array<typename FiniteFieldElementLookup<prime>::value_t, prime*prime> make_lookup_addition()
  {
    std::array<typename FiniteFieldElementLookup<prime>::value_t, prime*prime> result;
    for (unsigned int i = 0; i < prime; ++i)
      for (unsigned int j = 0; j < prime; ++j)
        result[prime*i + j] = (i + j) % prime;
    return result;
  }
  template <int prime>
  std::array<typename FiniteFieldElementLookup<prime>::value_t, prime*prime> FiniteFieldElementLookup<prime>::lookup_addition = make_lookup_addition<prime>();

  // Initialise the lookup multiplication table
  template <int prime>
  std::array<typename FiniteFieldElementLookup<prime>::value_t, prime*prime> make_lookup_multiplication()
  {
    std::array<typename FiniteFieldElementLookup<prime>::value_t, prime*prime> result;
    for (unsigned int i = 0; i < prime; ++i)
      for (unsigned int j = 0; j < prime; ++j)
        result[prime*i + j] = (i * j) % prime;
    return result;
  }
  template <int prime>
  std::array<typename FiniteFieldElementLookup<prime>::value_t, prime*prime> FiniteFieldElementLookup<prime>::lookup_multiplication = make_lookup_multiplication<prime>();

  // Initialise the lookup inverse table
  template <int prime>
  std::array<typename FiniteFieldElementLookup<prime>::value_t, prime> make_lookup_inverse()
  {
    std::array<typename FiniteFieldElementLookup<prime>::value_t, prime> result;
    for (unsigned int i = 0; i < prime; ++i) result[i] = static_cast<typename FiniteFieldElementLookup<prime>::value_t>(pow(i, prime - 2)) % prime;
    return result;
  }
  template <int prime>
  std::array<typename FiniteFieldElementLookup<prime>::value_t, prime> FiniteFieldElementLookup<prime>::lookup_inverse = make_lookup_inverse<prime>();

  // Initialise the lookup has square root table
  template <int prime>
  std::array<bool, prime> make_lookup_has_square_root()
  {
    std::array<bool, prime> result;
    for (unsigned int i = 0; i < prime; ++i) result[(i*i) % prime] = true;
    return result;
  }
  template <int prime>
  std::array<bool, prime> FiniteFieldElementLookup<prime>::lookup_has_square_root = make_lookup_has_square_root<prime>();

  // Initialise the lookup square root table
  template <int prime>
  std::array<typename FiniteFieldElementLookup<prime>::value_t, prime> make_lookup_square_root()
  {
    std::array<typename FiniteFieldElementLookup<prime>::value_t, prime> result;
    for (unsigned int i = 0; i < prime; ++i) result[(i*i) % prime] = i;
    return result;
  }
  template <int prime>
  std::array<typename FiniteFieldElementLookup<prime>::value_t, prime> FiniteFieldElementLookup<prime>::lookup_square_root = make_lookup_square_root<prime>();
}

#endif
