// From:
// http://stackoverflow.com/questions/147267/easy-way-to-use-variables-of-enum-types-as-string-in-c
// (http://stackoverflow.com/a/202511/2148414)

// Changes
// Oct 2013 davide.gerbaudo@gmail.com : add the EnumType2str wrapper function

// expansion macro for enum value definition
#define ENUM_VALUE(name,assign) name assign,

// expansion macro for enum to string conversion
#define ENUM_CASE(name,assign) case name: return #name;

// expansion macro for string to enum conversion
#define ENUM_STRCMP(name,assign) if (!strcmp(str,#name)) return name;

/// declare the access function and define enum values
#define DECLARE_ENUM(EnumType,ENUM_DEF) \
  enum EnumType { \
  ENUM_DEF(ENUM_VALUE) \
  }; \
  const char* GetString(EnumType dummy); \
  const char* EnumType ## 2str(EnumType dummy); \
  EnumType GetEnumValue(const char *string); \

/// define the access function names
#define DEFINE_ENUM(EnumType,ENUM_DEF) \
  const char* GetString(EnumType value) \
  { \
  switch(value) \
    { \
  ENUM_DEF(ENUM_CASE) \
 default: return ""; /* handle input error */ \
 } \
  } \
  const char* EnumType ## 2str(EnumType dummy) \
  { \
  return GetString(dummy); \
  } \
  EnumType GetEnumValue(const char *str) \
  { \
  ENUM_DEF(ENUM_STRCMP) \
  return (EnumType)0; /* handle input error */ \
  } \



