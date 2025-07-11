//
// Created by oleksii on 11.07.25.
//

#ifndef ANALYSISTREE_MYCLASS1_H
#define ANALYSISTREE_MYCLASS1_H

#include <iostream>
#include <string>

class MyClass1 { // highlights unused class
 public:

  void PrintHello() { std::cout << "Hello"; } // highlights unused function (gray) and proposal to make static (yellow)
  void PrintName() { std::cout << fName; } // highlights unused function (gray) however NO proposal to make const

 private:
  std::string fName;
  int fNumber; // highlights unused class member
};

#endif//ANALYSISTREE_MYCLASS1_H
