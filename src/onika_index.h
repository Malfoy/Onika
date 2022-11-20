#ifndef __ONIKA_INDEX_H__
#define __ONIKA_INDEX_H__

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
#include "onika_sketch.h"
#include <math.h>



using namespace std;
using kmer = uint64_t;
using gid = uint32_t;
const uint32_t mutex_number=65536;

class Index {
	public:
		//CONSTANTS
		uint32_t K;//kmer size
		uint32_t F;//fingerprint used
		uint32_t W;//fingerprint size
		uint32_t M;//Minhash size
		uint32_t mask_M;//minhash mask
		uint32_t H;//Hll size (M+H=W)
		uint32_t maximal_remainder;//2^H-1
		uint32_t lF;//log2(F)
		int32_t fingerprint_range;//2^w
		uint64_t mask_fingerprint;//2^(64-lf)-1
		uint64_t expected_gemome_size;
		uint64_t offsetUpdatekmer;
		uint32_t min_score;

		//VARIABLE
		string filename;
		uint32_t genome_numbers;//Number of genomes
		vector<gid>* Buckets;
		omp_lock_t lock[mutex_number];
		vector<string> filenames;
		zstr::ofstream* outfile;
		std::vector<Genome> infos;   // Array containing all Genomes informations.

		/**
		 * \brief Default constructor.
		 */
		Index(uint32_t F, uint32_t K, uint32_t W, uint32_t H, string filename);
		Index(const string& filestr, const string ifilename);


		/**
		 * \brief Destructor.
		 */
		~Index();

		/*
		 * void compress_index() const{
		 */

		/**
		 * \brief Returns the informations about the genome number i (starting from 0).
		 *
		 * Be aware that there is no bound verification.
		 *
		 * \return Returns the informations about the genome number i (starting from 0).
		 */
		inline const Genome &operator[](size_t i) const {
			return infos[i];
		}
		/**
		 * \brief Returns the number of indexed genomes.
		 *
		 * \return Returns the number of indexed genomes.
		 */
		inline size_t getNbGenomes() const {
			return genome_numbers;
		}


		inline bool exists_test (const std::string& name)const {
			struct stat buffer;
			return (stat (name.c_str(), &buffer) == 0);
		}		

		void get_filename(const string& filestr);
		void insert_file(const string& filestr,uint32_t identifier, Sketch& iSketch);

		uint64_t asm_log2(const uint64_t x) const;
		void Biogetline(zstr::ifstream* in,string& result,char type)const;
		void Biogetline(zstr::ifstream* in,string& result,char type,string& header)const ;
		char get_data_type(const string& filename)const;
};

#endif
