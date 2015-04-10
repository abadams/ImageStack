#ifndef IMAGESTACK_IMAGESTACK_H
#define IMAGESTACK_IMAGESTACK_H

// We never want SDL when used as a library
#define NO_SDL

// includes that don't survive well in the namespace
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <list>
#include <set>

#ifdef WIN32
#include <winsock2.h>
#else
#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <netdb.h>
#endif

#include <math.h>

$INCLUDES

#undef NO_SDL

#endif
