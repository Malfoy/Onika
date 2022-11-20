#ifndef __ONIKA_SKETCH_H__
#define __ONIKA_SKETCH_H__

/*
 * #include <codecfactory.h>  To use Compression Lib
 * #include <intersection.h>
 */
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include <algorithm>
#include <mutex>
#include <unordered_set>
#include <sys/types.h>
#include <sys/stat.h>
#include <omp.h>
#include "zstr.hpp"
#include "genome.h"
#include <math.h>



using namespace std;

struct Info {
	uint32_t ID;
	uint32_t position;
};


class Sketch {
	public:
		string filename;
		vector<string> filenames;
		vector<Info> info_array;   // Array containing ID and Position

		/**
		 * \brief Default Sketch constructor.
		 */
		Sketch(char* inString);

		/**
		 * \brief Sketch Destructor.
		 */
		~Sketch();

		void fasta_sketch(const string& filestr);
		void insert_sketch(void);
		void query_sketch(Info isketch);

	private:
		char* iString;

};

#endif
