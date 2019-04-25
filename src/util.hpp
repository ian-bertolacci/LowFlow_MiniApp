#ifndef UTIL_HPP
#define UTIL_HPP

// Utility Macros
// Create Valid C string from tokens
// Note: tokens with un-escaped quotes *will NOT* break this!
#define STRINGIZE(x) STRINGIZE_NO_PREPROCESS(x)
#define STRINGIZE_NO_PREPROCESS(x) #x

#endif
