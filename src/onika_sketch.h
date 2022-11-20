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


using kmer = uint64_t;
using namespace std;
using gid = uint32_t;

struct Info {
	uint32_t ID;
	uint32_t position;
};


class Sketch {
	public:
		uint64_t offsetUpdatekmer;
		uint32_t K;//kmer size
		uint32_t F;//fingerprint used
		uint32_t lF;//log2(F) 
			    //string filename;
			    //vector<string> filenames;
			    //vector<Info> info_array;   // Array containing ID and Position
		vector<gid>* Buckets;
		uint64_t expected_gemome_size;
		int32_t fingerprint_range;//2^w
		omp_lock_t lock[65536];
		
		
		/**
		 * \brief Default Sketch constructor.
		 */
		Sketch();
		Sketch(uint32_t F, uint32_t K);
		/**
		 * \brief Sketch Destructor.
		 */
		~Sketch();
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
		void insert_sketch(void);
		void merge_sketch(vector<int32_t>& sketch1,const vector<int32_t>& sketch2)const;
		void sketch_densification(vector<int32_t>& sketch, uint empty_cell) const;
		void compute_sketch(const string& reference, vector<int32_t>& sketch) const;
		//HERE we only select the minimal hashes without computing the HMH fingerprint
		void compute_sketch_kmer(const string& reference, vector<uint64_t>& sketch) const;
		void insert_sketch(const vector<int32_t>& sketch,uint32_t genome_id);

		uint64_t get_perfect_fingerprint(uint64_t hashed) const;
};

#endif
