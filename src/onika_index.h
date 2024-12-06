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
#include <map>
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
#include <iomanip>


using namespace std;
using kmer = uint64_t;
using gid = uint32_t;
const uint32_t mutex_number=65536;



/**
 * @class Index
 * @brief A class representing a genomic index
 *
 * This class is involved in creating and manipulating a genomic index.
 * It includes methods for inserting and querying files and sketches,
 * hash functions, and utility methods for processing genomic sequences.
 */
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
		uint64_t fingerprint_range;//2^w
		uint64_t mask_fingerprint;//2^(64-lf)-1
		uint64_t mi; // -1
		uint64_t offsetUpdatekmer;
		uint32_t min_score;

		//VARIABLE
		string filename;
		uint32_t genome_numbers;//Number of genomes
		vector<gid>* Buckets;
		vector<uint16_t>* Buckets_pos;
		omp_lock_t lock[mutex_number];
		vector<string> filenames;
		zstr::ofstream* outfile;
		std::vector<Genome> infos;   // Array containing all Genomes informations.
		std::map<std::string, std::vector<uint64_t>> file_sketches; // map to store sketches
		/**
		 * @brief Construct a new Index object with several parameters.
		 *
		 * @param F Some parameter F
		 * @param K Kmer size
		 * @param W Some parameter W
		 * @param E Expected genome size
		 * @param filename Name of the file
		 */
		Index(uint32_t F, uint32_t K, uint32_t W, uint32_t E, string filename);

		/**
		 * @brief Construct a new Index object with a file string and an index filename.
		 *
		 * @param filestr Some file string
		 * @param ifilename Index filename
		 */
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
		/// @brief Returns the number of indexed genomes.
		inline size_t getNbGenomes() const {
			return genome_numbers;
		}

		/// @brief Test if a file with given name exists.
		inline bool exists_test (const std::string& name)const {
			struct stat buffer;
			return (stat (name.c_str(), &buffer) == 0);
		}

		/**
		 * @brief Extract filename from the given string.
		 *
		 * @param filestr The string containing the filename.
		 */
		void get_filename(const string& filestr);
		/**
		 * @brief Insert a file into the genomic index.
		 *
		 * @param filestr The file to insert.
		 * @param identifier The identifier for the file in the index.
		 */
		void insert_file(const string& filestr,uint32_t identifier);
		/** 
		 * @brief Compute the base-2 logarithm of x.
		 *
		 * @param x The number to compute the logarithm for.
		 * @return The base-2 logarithm of x.
		 */
		uint64_t asm_log2(const uint64_t x) const;
		void Biogetline(zstr::ifstream* in,string& result,char type)const;
		void Biogetline(zstr::ifstream* in,string& result,char type,string& header)const ;
		char get_data_type(const string& filename)const;
		// Similar Doxygen comments can be written for the rest of the methods.

		/** 
		 * @brief Get the reverse complement of a kmer.
		 *
		 * @param min The kmer to compute the reverse complement for.
		 * @return The reverse complement of the kmer.
		 */
		kmer rcb(kmer min)const;
		kmer str2numstrand(const string& str)const;
		uint64_t revhash64 ( uint64_t x ) const;
		uint64_t nuc2int(char c) const;
		uint64_t nuc2intrc(char c) const;
		void update_kmer(kmer& min, char nuc)const;

		void update_kmer_RC(kmer& min, char nuc)const;

		uint64_t unrevhash64 ( uint64_t x ) const;
		uint64_t hash_family(const uint64_t x, const uint factor)const;
		/**
		 * @brief Query the genomic index with a file.
		 *
		 * @param filestr The file to query the index with.
		 * @return A vector of identifiers of genomes related to the file.
		 */
		void query_file(const string& filestr);
		void query_file_of_file(const string& filestr);
		vector<uint32_t> query_sketch(const vector<uint64_t>& sketch) const;


		/**
		 * @brief Compute a sketch from a FASTA file and store it in the index.
		 *
		 * @param filestr Path to the FASTA file.
		 */
		void fasta_sketch(const string& filestr);

		/**
		 * @brief Merge two sketches into a single one.
		 *
		 * @param sketch1 The first sketch to be merged.
		 * @param sketch2 The second sketch to be merged.
		 */
		void merge_sketch(vector<int32_t>& sketch1, const vector<int32_t>& sketch2)const;

		/**
		 * @brief Densify the sketch by filling the empty cells.
		 *
		 * @param sketch The sketch to densify.
		 * @param empty_cell The value to fill the empty cells with.
		 */
		void sketch_densification(vector<uint64_t>& sketch, uint empty_cell) const;

		/**
		 * @brief Compute a sketch from a reference string using kmers and store it in the index.
		 *
		 * @param reference The reference string.
		 * @param sketch The sketch to store the result.
		 */
		void compute_sketch_kmer(const string& reference, vector<uint64_t>& sketch) const;

		/**
		 * @brief Compute a sketch from a reference string and store it in the index.
		 *
		 * @param reference The reference string.
		 * @param sketch The sketch to store the result.
		 */
		void compute_sketch(const string& reference, vector<uint64_t>& sketch) const;

		/**
		 * @brief Insert a sketch into the genomic index.
		 *
		 * @param sketch The sketch to insert.
		 * @param genome_id The identifier for the genome the sketch is related to.
		 */
		void insert_sketch(const vector<uint64_t>& sketch,uint32_t genome_id);

		/**
		 * @brief Write the current state of the index to a file on disk.
		 *
		 * @param filestr The path to the file to write to.
		 */
		void dump_index_disk(const string& filestr)const ;

		/**
		 * @brief Compute the perfect fingerprint of a hashed value.
		 *
		 * @param hashed The hashed value.
		 * @return The perfect fingerprint of the hashed value.
		 */
		uint64_t get_perfect_fingerprint(uint64_t hashed) const;
		uint64_t get_perfect_fingerprint2(uint64_t hashed) const;

		/**
		 * @brief Print the internal state of the index as a matrix.
		 */
		void print_matrix() const;


};

#endif
