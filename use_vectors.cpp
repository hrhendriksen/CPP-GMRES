#include <stdlib.h>
#include <iostream>
#include <cassert>
#include "Exception.hpp"
#include "Vector.hpp"
#include "utility"

int main(int argc, char* argv[])
{
  //Test initialisation with array
  double arr[2] = {1,2};
  const Vector a_vector(arr,2);

  // Test reading with parentheses
  std::cout << "first :" << a_vector(1) << " second :" << a_vector(2) << "\n";

  //This would  produce a compiler warning (there is no default constructor)
  // Vector badly_formed;
  Vector new_vector(3);
  new_vector = a_vector;
  //Show that friends and methods can be used to do the same thing
  assert ( a_vector.norm() == norm(a_vector));
  assert ( a_vector.norm(3) == norm(a_vector, 3));

  // Test reshape operator
  std::cout << "a_vector = " << a_vector << "\n";
  std::cout << "reshaped a_vector" << cut(a_vector, 6) << "\n";

  // Test size operator
  std::cout <<  "Size of vector a " << "(" << size(a_vector).first <<", " <<size(a_vector).second<<")\n";

  //Test write to CSV operator
  writetoCSV(a_vector);


  // Test for warnings and exceptions
  Vector bigger_vector(3);
  Vector smaller_vector(1);
  std::cout << "The following produces a warning\n";
  // This produces a warning
  bigger_vector = a_vector;
  std::cout << "The following throws an exception\n";
  // This throws an exception
  try
  {
      smaller_vector = a_vector;
  }
  catch (Exception &ex)
  {
      ex.DebugPrint();
  }
  std::cout << "The following throws an exception\n";
  // This throws an exception
  try
  {
      a_vector/2.0;//Okay  
      a_vector/0.0;
  }
  catch (Exception &ex)
  {
      ex.DebugPrint();
  }
  exit(0);
}
