#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H
#include <exception>

/* This exception only gets used to prematurely exit the program
 * __in instances when we should__!!! Pretty much exclusively limited
 * to printing Version or Help info from the command line 
 */
struct SilentException : public std::exception {
  const char *what() const throw() {
    return "Have a nice day! :)";
  }
};

#endif
