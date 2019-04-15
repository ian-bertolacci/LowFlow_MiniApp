#ifndef UTIL_HPP
#define UTIL_HPP

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <stdbool.h>

// Utility Macros
// Create Valid C string from tokens (Note: tokens with un-escaped quotes will break this!)
#define STRINGIZE(x) STRINGIZE2(x)
#define STRINGIZE2(x) #x

#endif
