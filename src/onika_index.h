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
#include <math.h>



using namespace std;
using kmer = uint64_t;

struct Info {
	uint32_t ID, position;
};


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

		void insert_file(const string& filestr);
		void query_file(const string& filestr);

		int32_t get_fingerprint(uint64_t hashed)const;
		void merge_sketch( vector<int32_t>& sketch1,const vector<int32_t>& sketch2)const;
		uint64_t asm_log2(const uint64_t x) const;
		uint64_t revhash64 ( uint64_t x ) const;
		uint64_t unrevhash64 ( uint64_t x ) const;
		void Biogetline(zstr::ifstream* in,string& result,char type)const;
		void Biogetline(zstr::ifstream* in,string& result,char type,string& header)const ;
		char get_data_type(const string& filename)const;
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
