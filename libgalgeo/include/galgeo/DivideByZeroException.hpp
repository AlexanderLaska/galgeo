#include <stdexcept> 

class DivideByZeroException : public std::runtime_error
{
  public:
    DivideByZeroException() : std::runtime_error("A division of a GaloisScalar by zero had been tried.") {}
};
