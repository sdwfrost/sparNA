#ifndef MEMMAC_H
#define MEMMAC_H

#define ALLOCMEMORY(X,S,T,N) allocmemory(__FILE__,__LINE__, X, S, sizeof(T), N)
#define FREEMEMORY(X,P) freememory(__FILE__, __LINE__, X, P)

#endif
