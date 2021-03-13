#include <iostream>

int main()
{
  char name[256];
  int count = 0;
  
  // Print "Hello, world!" on the standard output...
  std::cout << "Hello, world!"  // <- Notice: no ';' here: 1 statement / 2 lines
				     
  // ...then advance to the next line
	    << std::endl;      // <- 2nd line of the statement

  std::cout << "What is your name?" << std::endl;
					     
  // Read the name from the standard input
  std::cin >> name;
	       
  std::cout << "Give a positive number: ";
  std::cin >> count;

  std::cout << "Hey";
  for (int i = 0; i < count; i++) {   // <-- Here goes the loop!
    std::cout << "-hey";
  }

  std::cout << ", " << name << "!" << std::endl;
  
  return 0;
}

