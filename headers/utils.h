#ifndef __UTILS_H__
#define __UTILS_H__

#include <iostream>
#include <string>
#include "../include/zstr.hpp"

using std::string;
using kmer = uint64_t;

uint64_t asm_log2(const uint64_t x);
char get_data_type(const string& filename);
uint64_t nuc2int(char c);
uint64_t nuc2intrc(char c);
kmer str2numstrand(const string& str);
uint64_t revhash64(uint64_t x);
uint64_t unrevhash64(uint64_t x);
uint64_t hash_family(const uint64_t x, const uint factor);



#endif