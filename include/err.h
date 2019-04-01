#pragma once
#include <stdio.h>

//error handle
#define Fatal_Error(ss)                                                \
    {                                                                  \
        printf("Fatal error '%s' at %s:%d\n", ss, __FILE__, __LINE__); \
        exit(1);                                                       \
    }
    