#include <iostream>
#include <string>
#include <math.h>
#include <random>

#include "../../../include/galgeo/GaloisScalar.hpp"
#include "../../../include/galgeo/GaloisOneZeroTensor.hpp"
#include "../../../include/galgeo/GaloisZeroOneTensor.hpp"
#include "../../../include/galgeo/GaloisOneOneTensor.hpp"
#include "../../../include/galgeo/GaloisZeroTwoTensor.hpp"
#include "../../../include/galgeo/GaloisTwoZeroTensor.hpp"

#include "../../../include/galgeo/FiniteProjectiveGeometry.hpp"

using namespace galgeo;


int main()
{
  // Save the prime number, the dimension, the number of points and hyperplanes, and the number of runs in easy interpretable identifiers 
  int primeNumber{13};
  set_prime(primeNumber);

  const unsigned int dimension{3};
  set_dimension(dimension);

  std::cout << "TEST!" << std::endl;
  
  return 0;
}
