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
		uint32_t E;//Expected genome size (5000000)
		uint32_t maximal_remainder;//2^H-1
		uint32_t lF;//log2(F)
		int32_t fingerprint_range;//2^w
		uint64_t mask_fingerprint;//2^(64-lf)-1
		uint64_t mi; // -1
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
		Index(uint32_t F, uint32_t K, uint32_t W, uint32_t E, string filename);
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
		void insert_file(const string& filestr,uint32_t identifier);

		uint64_t asm_log2(const uint64_t x) const;
		void Biogetline(zstr::ifstream* in,string& result,char type)const;
		void Biogetline(zstr::ifstream* in,string& result,char type,string& header)const ;
		char get_data_type(const string& filename)const;

		kmer rcb(kmer min)const;
		kmer str2numstrand(const string& str)const;
		uint64_t revhash64 ( uint64_t x ) const;
		uint64_t nuc2int(char c) const;
		uint64_t nuc2intrc(char c) const;
		void update_kmer(kmer& min, char nuc)const;

		void update_kmer_RC(kmer& min, char nuc)const;

		uint64_t unrevhash64 ( uint64_t x ) const;
		uint64_t hash_family(const uint64_t x, const uint factor)const;


		void fasta_sketch(const string& filestr);
		void merge_sketch(vector<int32_t>& sketch1,const vector<int32_t>& sketch2)const;
		void sketch_densification(vector<uint64_t>& sketch, uint empty_cell) const;
		void compute_sketch_kmer(const string& reference, vector<uint64_t>& sketch) const;
		void compute_sketch(const string& reference, vector<uint64_t>& sketch) const;
		void insert_sketch(const vector<uint64_t>& sketch,uint32_t genome_id);

		uint64_t get_perfect_fingerprint(uint64_t hashed) const;


};

#endif
