
#pragma once

#define DIMENSIONS 1


#if DIMENSIONS == 3

#define DLIST(a,b,c) a , b , c 

#endif


#if DIMENSIONS == 2
#define DLIST(a,b,c) a , b  
#endif

#if DIMENSIONS == 1
#define DLIST(a, b, c ) a   
#endif
