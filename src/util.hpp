#ifndef UTIL_HPP
#define UTIL_HPP

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <cctype>
#include <cstdbool>

// Utility Macros
// Create Valid C string from tokens (Note: tokens with un-escaped quotes will break this!)
#define STRINGIZE(x) STRINGIZE2(x)
#define STRINGIZE2(x) #x

#endif
