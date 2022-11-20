#include "strict_fstream.hpp"
#include "zstr.hpp"
#include "onika_index.h"
#include "onika_sketch.h"
#include "common.h"


using namespace std;
const int bufferSize = 10000;



Index::Index(uint32_t ilF=10, uint32_t iK=31,uint32_t iW=8,uint32_t iH=4, const string ifilename="onikaOutput.gz") {
	filename=ifilename;
	lF=ilF;
	K=iK;
	W=iW;
	H=iH;
	F=1<<lF;
	M=W-H;
	fingerprint_range=1<<W;
	mask_M=(1<<M)-1;
	maximal_remainder=(1<<H)-1;
	genome_numbers=0;
	Buckets=new vector<gid>[fingerprint_range*F];
	offsetUpdatekmer=1;
	offsetUpdatekmer<<=2*K;
	for(uint32_t i=0; i<mutex_number;++i) {
		omp_init_lock(&lock[i]);
	}
	outfile=new zstr::ofstream(filename);
}

Index::Index(const string& filestr, const string ifilename) {
	filename=ifilename;
	zstr::ifstream in(filestr,ios::binary);
	in.read(reinterpret_cast< char*>(&lF), sizeof(lF));
	in.read(reinterpret_cast< char*>(&K), sizeof(K));
	in.read(reinterpret_cast< char*>(&H), sizeof(H));
	in.read(reinterpret_cast< char*>(&W), sizeof(W));
	in.read(reinterpret_cast< char*>(&min_score), sizeof(min_score));
	in.read(reinterpret_cast< char*>(&genome_numbers), sizeof(genome_numbers));
	F=1<<lF;
	M=W-H;
	fingerprint_range=1<<W;
	mask_M=(1<<M)-1;
	maximal_remainder=(1<<H)-1;
	Buckets=new vector<gid>[fingerprint_range*F];
	offsetUpdatekmer=1;
	offsetUpdatekmer<<=2*K;
	for(uint32_t i=0; i<mutex_number;++i) {
		omp_init_lock(&lock[i]);
	}
	for(uint i(0);i<fingerprint_range*F;++i){
		uint32_t size;
		in.read(reinterpret_cast< char*>(&size), sizeof(size));
		if(size!=0){
			Buckets[i].resize(size,0);
			in.read(reinterpret_cast< char*>(&(Buckets[i][0])), size*sizeof(gid));
		}
	}
	string genome_name;
	for(uint i(0);i<genome_numbers;++i){
		getline(in,genome_name);
		filenames.push_back(genome_name);
	}

	outfile=new zstr::ofstream(filename,ios::binary);

}


Index::~Index() {
	//delete[] Buckets;
	outfile->close();
}


uint64_t Index::asm_log2(const uint64_t x) const {
	uint64_t y;
	asm ( "\tbsr %1, %0\n"
			: "=r"(y)
			: "r" (x)
	    );
	return y;
}

uint64_t Index::revhash64 ( uint64_t x ) const {
	x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
	x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
	x = ( ( x >> 32 ) ^ x );
	return x;
}



uint64_t Index::unrevhash64(uint64_t x) const{
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x);
	return x;
}




void Index::get_filename(const string& filestr) {
	DEBUG_MSG("Opening file : '"<<filestr<<"'");
	ifstream in(filestr);
	if (!in) {
		cout << "Unable to open the file '" << filestr << "'" << endl;
		exit(0);
	}
#pragma omp parallel
	{
		string ref;

		Sketch iSketch;
		uint32_t id;
		while(not in.eof()) {
			bool go = false;
#pragma omp critical (input)
			{
				getline(in,ref);
			}
			DEBUG_MSG("Getline from file :'"<<filestr<<"' = '"<<ref<<"'");
#pragma omp critical (genome_numbers)
			{
				if(ref.size()>2){
					if(exists_test(ref)) {
						id=genome_numbers;
						DEBUG_MSG("Genome numbers: "<<genome_numbers);
						genome_numbers++;
						filenames.push_back(ref);
						go=true;
					}
				}
			}
			if(go) {
				DEBUG_MSG("Adding file :'"<<ref<<"'");
				insert_file(ref,id,iSketch);
				DEBUG_MSG("File: '"<<ref<<"' added");
			}

			ref.clear();
		}
	}
}



void Index::insert_file(const string& filestr,uint32_t identifier,Sketch iSketch) {
	char type=get_data_type(filestr);
	zstr::ifstream in(filestr);
	string ref;
	vector<uint64_t> kmer_sketch;
	vector<int32_t> sketch(F,-1);
	while(not in.eof()) {
		Biogetline(&in,ref,type);
		if(ref.size()>K) {
			iSketch.compute_sketch(ref,sketch);
		}
		ref.clear();
	}   
	iSketch.insert_sketch(sketch,identifier);
}



/**
 * \note The get_fingerprint take hashed a uint64_t.
 * \brief Get the fingerprint.
 *
 * \param idx hashed to set.
 * \brief Returns the result if any.
 *
 */

int32_t Index::get_fingerprint(uint64_t hashed)const {
	int32_t result;
	result = hashed&mask_M;//we keep the last bits for the minhash part
	uint32_t ll=asm_log2(hashed);//we compute the log of the hash
	uint32_t size_zero_trail(64-ll-1);
	int remaining_nonzero=maximal_remainder-size_zero_trail;
	remaining_nonzero=max(0,remaining_nonzero);
	// if the log is too high can cannont store it on H bit here we sature
	result+=remaining_nonzero<<M;// we concatenant the hyperloglog part with the minhash part
	return result;
}





atomic<uint32_t> genomes_downloaded(0);
atomic<uint64_t> bases_downloaded(0);

void Index::Biogetline(zstr::ifstream* in,string& result,char type,string& header)const {
	string discard;
	result.clear();
	switch(type){
		case 'Q':
			getline(*in,header);
			getline(*in,result);
			getline(*in,discard);
			getline(*in,discard);
			break;
		case 'A':
			getline(*in,header);
			char c=in->peek();
			while(c!='>' and c!=EOF){
				getline(*in,discard);
				result+=discard;
				c=in->peek();
			}
			break;
	}
	if(result.size()< K){
		result.clear();
		header.clear();
	}
}



void Index::Biogetline(zstr::ifstream* in,string& result,char type)const {
	string discard;
	result.clear();
	switch(type){
		case 'Q':
			getline(*in,discard);
			getline(*in,result);
			getline(*in,discard);
			getline(*in,discard);
			break;
		case 'A':
			getline(*in,discard);
			char c=in->peek();
			while(c!='>' and c!=EOF){
				getline(*in,discard);
				result+=discard;
				c=in->peek();
			}
			break;
	}
	if(result.size()< K){
		result.clear();
	}
}


char Index::get_data_type(const string& filename)const{
	if(filename.find(".fq")!=string::npos){
		return 'Q';
	}
	if(filename.find(".fastq")!=string::npos){
		return 'Q';
	}
	return 'A';
}

