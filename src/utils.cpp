#include "../headers/utils.h"

using namespace std;

/** 
 * @brief Compute the base-2 logarithm of an unsigned 64-bit integer
 * 
 * This function uses the Bit Scan Reverse (bsr) instruction to find the position of the most significant bit set.
 * 
 * @param x The number to compute the logarithm for.
 * @return The base-2 logarithm of x, or undefined if x is zero.
 */
uint64_t asm_log2(const uint64_t x){
	uint64_t y;
	asm ( "\tbsr %1, %0\n"
			: "=r"(y)
			: "r" (x)
	    );
	return y;
}


/** 
 * @brief Determines the data type of a file based on its extension.
 *
 * This function checks the file extension and returns 'Q' for fastq files and 'A' for fasta files.
 * 
 * @param filename The name of the file.
 * @return The type of file ('Q' for fastq or 'A' for fasta).
 */
char get_data_type(const string& filename){
	if(filename.find(".fq")!=string::npos){
		return 'Q';
	}
	if(filename.find(".fastq")!=string::npos){
		return 'Q';
	}
	return 'A';
}


/** 
 * @brief Converts a nucleotide to an integer.
 *
 * This function maps each nucleotide to an integer as follows: A->0, C->1, G->2, T->3. 
 * If the input character does not represent a valid nucleotide, the function will return 0.
 * 
 * @param c The nucleotide to convert.
 * @return The integer representation of the nucleotide.
 */
uint64_t nuc2int(char c){
	switch(c){
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
		default: return 0;
	}
	exit(0);
	return 0;
}


/** 
 * @brief Converts a nucleotide to an integer, considering the reverse complement.
 *
 * This function maps each nucleotide to an integer as follows: A->3, C->2, G->1, T->0. 
 * If the input character does not represent a valid nucleotide, the function will return 0.
 * 
 * @param c The nucleotide to convert.
 * @return The integer representation of the reverse complement of the nucleotide.
 */
uint64_t nuc2intrc(char c){
	switch(c) {
		case 'A': return 3;
		case 'C': return 2;
		case 'G': return 1;
		default: return 0;
	}
	exit(0);
	return 0;
}


/** 
 * @brief Converts a string of nucleotides to a k-mer.
 *
 * This function converts a string of nucleotides to a k-mer by shifting the current k-mer 
 * and adding the integer representation of the next nucleotide in the string.
 * 
 * @param str The string of nucleotides to convert.
 * @return The k-mer representation of the string.
 */
kmer str2numstrand(const string& str){
	uint64_t res(0);
	for(uint i(0);i<str.size();i++) {
		res<<=2;
		switch (str[i]){
			case 'A':res+=0;break;
			case 'C':res+=1;break;
			case 'G':res+=2;break;
			case 'T':res+=3;break;

			case 'a':res+=0;break;
			case 'c':res+=1;break;
			case 'g':res+=2;break;
			case 't':res+=3;break;
			default: return 0 ;break;
		}
	}
	return (res);
}


/** 
 * @brief Reverses the hash of a 64-bit integer.
 *
 * This function reverses the hash of a 64-bit integer by applying a reverse operation.
 * 
 * @param x The 64-bit integer to reverse the hash of.
 * @return The reversed hash.
 */
uint64_t revhash64(uint64_t x){
        // Perform xorshift operations
        x ^= x >> 12;  // Shift right and XOR
        x ^= x << 25;  // Shift left and XOR
        x ^= x >> 27;  // Shift right and XOR
        return x * 0x2545F4914F6CDD1D; // Multiply by a large prime constant
    }


/** 
 * @brief Unreverses the hash of a 64-bit integer.
 *
 * This function unreverses the hash of a 64-bit integer by applying an unreversing operation.
 * 
 * @param x The 64-bit integer to unreverse the hash of.
 * @return The unreversed hash.
 */
uint64_t unrevhash64(uint64_t x){
        x ^= (x >> 33);  // XOR with right shift by 33
        x *= 0xff51afd7ed558ccd; // Multiply by a constant
        x ^= (x >> 33);  // XOR with right shift by 33
        x *= 0xc4ceb9fe1a85ec53; // Multiply by another constant
        x ^= (x >> 33);  // Final XOR with right shift by 33
	return x;
}


/** 
 * @brief Hashing function that depends on a factor.
 *
 * This function calculates a new hash value based on the input integer and a factor.
 * The hashing function is a combination of unrevhash64 and revhash64.
 * 
 * @param x The input integer to hash.
 * @param factor The factor that influences the hash function.
 * @return The hashed value.
 */
uint64_t hash_family(const uint64_t x, const uint factor){
	return unrevhash64(x)+factor*revhash64(x);
}