#include "strict_fstream.hpp"
#include "zstr.hpp"
#include "onika_sketch.h"
#include "onika_index.h"
#include "common.h"



using namespace std;
const int bufferSize = 10000;

/**
 * \brief Default constructor.
 */
Sketch::Sketch() {
}


Sketch::Sketch(uint32_t ilF=10, uint32_t iK=31) {
	lF=ilF;
	K=iK;
	offsetUpdatekmer=1;
	offsetUpdatekmer<<=2*K;
}

Sketch::~Sketch() {
	//delete[] Buckets;
}


void Sketch::fasta_sketch(const string& filename) {

}

void Sketch::insert_sketch(void) {

}

//void Sketch::query_sketch(Info iSketch) {
//}
//



kmer Sketch::rcb(kmer min)const {
	kmer res(0);
	kmer offset(1);
	offset<<=(2*K-2);
	for(uint i(0); i<K;++i){
		res+=(3-(min%4))*offset;
		min>>=2;
		offset>>=2;
	}
	return res;
}



uint64_t Sketch::nuc2int(char c)const {
	switch(c){
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
		default: return 0;
	}
	exit(0);
	return 0;
}
uint64_t Sketch::nuc2intrc(char c)const {
	switch(c) {
		case 'A': return 3;
		case 'C': return 2;
		case 'G': return 1;
		default: return 0;
	}
	//~ cout<<"Unknow nucleotide: "<<c<<"!"<<endl;
	exit(0);
	return 0;
}

kmer Sketch::str2numstrand(const string& str)const {
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


uint64_t Sketch::revhash64 ( uint64_t x ) const {
	x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
	x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
	x = ( ( x >> 32 ) ^ x );
	return x;
}



uint64_t Sketch::unrevhash64(uint64_t x) const{
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x);
	return x;
}





void Sketch::update_kmer(kmer& min, char nuc)const {
	min<<=2;
	min+=nuc2int(nuc);
	min%=offsetUpdatekmer;
}



void Sketch::update_kmer_RC(kmer& min, char nuc)const {
	min>>=2;
	min+=(nuc2intrc(nuc)<<(2*K-2));
}


uint64_t Sketch::hash_family(const uint64_t x, const uint factor)const{
	return unrevhash64(x)+factor*revhash64(x);
}

void Sketch::sketch_densification(vector<int32_t>& sketch, uint empty_cell) const {
	uint step(0);
	uint size(sketch.size());
	while(empty_cell!=0){
		for(uint i(0);i<size;i++){
			if(sketch[i]!=-1){
				uint64_t hash=hash_family(sketch[i],step)%size;
				if(sketch[hash]==-1){
					sketch[hash]=sketch[i];
					empty_cell--;
					if(empty_cell==0){
						return;
					}
				}
			}
		}
		step++;
	}
}

uint64_t Sketch::get_perfect_fingerprint(uint64_t hashed)const {
	uint64_t result(hashed);
	result<<=1;
	result>>=1;
	uint64_t twopowern(1);
	twopowern<<=63;
	double frac=((double)(twopowern-result))/twopowern;
	frac=pow(frac,expected_gemome_size/F);
	frac=1-frac;
	return fingerprint_range*frac;
}


void Sketch::compute_sketch(const string& reference, vector<int32_t>& sketch) const {
	if(sketch.size()!=F) {
		sketch.resize(F,-1);
	}
	uint empty_cell(F);
	kmer S_kmer(str2numstrand(reference.substr(0,K-1)));//get the first kmer (k-1 bases)
	kmer RC_kmer(rcb(S_kmer));//The reverse complement
	for(uint i(0);i+K<reference.size();++i) {// all kmer in the genome
		Sketch::update_kmer(S_kmer,reference[i+K-1]);
		Sketch::update_kmer_RC(RC_kmer,reference[i+K-1]);
		kmer canon(min(S_kmer,RC_kmer));//Kmer min, the one selected
		uint64_t hashed=revhash64(canon);
		uint64_t bucket_id(unrevhash64(canon)>>(64-lF));//Which Bucket 
		int32_t fp=get_perfect_fingerprint(hashed);
		//MINHASH
		if(sketch[bucket_id]==-1){
			empty_cell--;
			sketch[bucket_id]=fp;
		}else if(sketch[bucket_id] > fp) {
			sketch[bucket_id]=fp;
		}
	}
	sketch_densification(sketch,empty_cell);
}



void Sketch::insert_sketch(const vector<int32_t>& sketch,uint32_t genome_id) {
	for(uint i(0);i<F;++i) {
		if(sketch[i]<fingerprint_range and sketch[i]>=0) {
			omp_set_lock(&lock[(sketch[i]+i*fingerprint_range)%mutex_number]);
			Buckets[sketch[i]+i*fingerprint_range].push_back(genome_id);
			omp_unset_lock(&lock[(sketch[i]+i*fingerprint_range)%mutex_number]);
		}
	}
}


void Sketch::merge_sketch( vector<int32_t>& sketch1,const vector<int32_t>& sketch2) const {
	for(uint i(0);i<sketch1.size();++i) {
		sketch1[i]=min(sketch1[i],sketch2[i]);
	}
}

