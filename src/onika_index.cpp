#include "strict_fstream.hpp"
#include "zstr.hpp"
#include "onika_index.h"
#include "onika_sketch.h"
#include "common.h"


using namespace std;
const int bufferSize = 10000;



Index::Index(uint32_t ilF=10, uint32_t iK=31,uint32_t iW=8,uint32_t iE=5000000, const string ifilename="onikaOutput.gz") {
	filename=ifilename;
	lF=ilF;
	K=iK;
	W=iW;
	E=iE;
	F=1<<lF;
	fingerprint_range=1<<W;
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
	in.read(reinterpret_cast< char*>(&W), sizeof(W));
	in.read(reinterpret_cast< char*>(&min_score), sizeof(min_score));
	in.read(reinterpret_cast< char*>(&genome_numbers), sizeof(genome_numbers));
	F=1<<lF;
	fingerprint_range=1<<W;
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
	delete[] Buckets;
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
				insert_file(ref,id);
				DEBUG_MSG("File: '"<<ref<<"' added");
			}

			ref.clear();
		}
	}
}



void Index::insert_file(const string& filestr,uint32_t identifier) {
	char type=get_data_type(filestr);
	zstr::ifstream in(filestr);
	string ref;
	vector<uint64_t> kmer_sketch;
	vector<int32_t> sketch(F,-1);
	while(not in.eof()) {
		Biogetline(&in,ref,type);
		if(ref.size()>K) {
			compute_sketch(ref,sketch);
		}
		ref.clear();
	}   
	insert_sketch(sketch,identifier);
}



/**
 * \note The get_fingerprint take hashed a uint64_t.
 * \brief Get the fingerprint.
 *
 * \param idx hashed to set.
 * \brief Returns the result if any.
 *
 */


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



kmer Index::rcb(kmer min)const {
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



uint64_t Index::nuc2int(char c)const {
	switch(c){
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
		default: return 0;
	}
	exit(0);
	return 0;
}
uint64_t Index::nuc2intrc(char c)const {
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

kmer Index::str2numstrand(const string& str)const {
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





void Index::update_kmer(kmer& min, char nuc)const {
	min<<=2;
	min+=nuc2int(nuc);
	min%=offsetUpdatekmer;
}



void Index::update_kmer_RC(kmer& min, char nuc)const {
	min>>=2;
	min+=(nuc2intrc(nuc)<<(2*K-2));
}


uint64_t Index::hash_family(const uint64_t x, const uint factor)const{
	return unrevhash64(x)+factor*revhash64(x);
}

void Index::sketch_densification(vector<int32_t>& sketch, uint empty_cell) const {
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

uint64_t Index::get_perfect_fingerprint(uint64_t hashed)const {
	uint64_t result(hashed);
	result<<=1;
	result>>=1;
	uint64_t twopowern(1);
	twopowern<<=63;
	double frac=((double)(twopowern-result))/twopowern;
	DEBUG_MSG("-----------------------------------------------------------");
	DEBUG_MSG("Frac value : '"<<frac<<"'");
	frac=pow(frac,E/F);
	DEBUG_MSG("E/F = '"<<E<<"/"<<F<<" = "<<E/F<<"'");
	DEBUG_MSG("New frac value : '"<<frac<<"'");
	frac=1-frac;
	DEBUG_MSG("Fingerprint range : '"<<fingerprint_range<<"'");
	DEBUG_MSG("Fingerprint range * frac =  '"<<fingerprint_range*frac<<"'");
	DEBUG_MSG("-----------------------------------------------------------");
	return fingerprint_range*frac;
}


void Index::compute_sketch(const string& reference, vector<int32_t>& sketch) const {
	if(sketch.size()!=F) {
		sketch.resize(F,-1);
	}
	uint empty_cell(F);
	kmer S_kmer(str2numstrand(reference.substr(0,K-1)));//get the first kmer (k-1 bases)
	kmer RC_kmer(rcb(S_kmer));//The reverse complement
	for(uint i(0);i+K<reference.size();++i) {// all kmer in the genome
		Index::update_kmer(S_kmer,reference[i+K-1]);
		Index::update_kmer_RC(RC_kmer,reference[i+K-1]);
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



void Index::insert_sketch(const vector<int32_t>& sketch,uint32_t genome_id) {
	for(uint i(0);i<F;++i) {
		if(sketch[i]<fingerprint_range and sketch[i]>=0) {
			omp_set_lock(&lock[(sketch[i]+i*fingerprint_range)%mutex_number]);
			Buckets[sketch[i]+i*fingerprint_range].push_back(genome_id);
			omp_unset_lock(&lock[(sketch[i]+i*fingerprint_range)%mutex_number]);
		}
	}
}


void Index::merge_sketch( vector<int32_t>& sketch1,const vector<int32_t>& sketch2) const {
	for(uint i(0);i<sketch1.size();++i) {
		sketch1[i]=min(sketch1[i],sketch2[i]);
	}
}
