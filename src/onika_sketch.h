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
using kmer = uint64_t;

struct Info {
	uint32_t ID, position;
};


class Sketch {
	public:
		string filename;
		vector<string> filenames;
		vector<Sketch> info;   // Array containing ID and Position

		/**
		 * \brief Default Sketch constructor.
		 */
		Sketch();

		/**
		 * \brief Sketch Destructor.
		 */
		~Sketch();

		void fasta_sketch(const string& filestr);
		void insert_sketch(void);
		void query_sketch(Info sketch);

};

#endif
