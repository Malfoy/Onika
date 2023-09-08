#include "strict_fstream.hpp"
#include "zstr.hpp"
#include "onika_index.h"
#include "onika_sketch.h"
#include "common.h"



using namespace std;
const int bufferSize = 10000;


/** 
 * @brief Constructor for the Index class
 *
 * @param ilF Default value is 10
 * @param iK Default value is 31
 * @param iW Default value is 8
 * @param iE Default value is 5000000
 * @param ifilename Default value is "onikaOutput.gz"
 */
Index::Index(uint32_t ilF=10, uint32_t iK=31,uint32_t iW=8,uint32_t iE=5000000, const string ifilename="onikaOutput.gz") {
	filename=ifilename;
	lF=ilF;
	W=iW;
	K=iK;
	E=iE;
	F=(uint64_t)1<<lF;
	fingerprint_range=(uint64_t)1<<W;
	mi = -1; //mi =-1
	genome_numbers=0;
	//cout<<fingerprint_range*F<<endl;
	Buckets=new vector<gid>[fingerprint_range*F];
	Buckets_pos=new vector<uint16_t>[fingerprint_range*F];
	offsetUpdatekmer=1;
	offsetUpdatekmer<<=2*K;
	for(uint32_t i=0; i<mutex_number;++i) {
		omp_init_lock(&lock[i]);
	}
	outfile=new zstr::ofstream(filename);
}


/** 
 * @brief Constructor for the Index class that initializes the class with a binary file and filename.
 *
 * This constructor reads data from a binary file to initialize member variables and setup the index object.
 * 
 * @param filestr The name of the binary file from which to read the data.
 * @param ifilename The name of the output file.
 */
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
	Buckets_pos=new vector<uint16_t>[fingerprint_range*F];
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
			in.read(reinterpret_cast< char*>(&(Buckets_pos[i][0])), size*sizeof(uint16_t));
		}
	}
	string genome_name;
	for(uint i(0);i<genome_numbers;++i){
		getline(in,genome_name);
		filenames.push_back(genome_name);
	}

	outfile=new zstr::ofstream(filename,ios::binary);

}


/** 
 * @brief Destructor for the Index class
 *
 * This destructor frees memory allocated to the Buckets array.
 */
Index::~Index() {
	delete[] Buckets;
	delete[] Buckets_pos;
	outfile->close();
}


/** 
 * @brief Compute the base-2 logarithm of an unsigned 64-bit integer
 * 
 * This function uses the Bit Scan Reverse (bsr) instruction to find the position of the most significant bit set.
 * 
 * @param x The number to compute the logarithm for.
 * @return The base-2 logarithm of x, or undefined if x is zero.
 */
uint64_t Index::asm_log2(const uint64_t x) const {
	uint64_t y;
	asm ( "\tbsr %1, %0\n"
			: "=r"(y)
			: "r" (x)
	    );
	return y;
}

/** 
 * @brief Process the given file, and for each line in the file (which should represent a filename), 
 * it checks whether the file exists and if it does, it inserts it into the index.
 * 
 * @param filestr The name of the file that contains the list of filenames to process.
 */
void Index::get_filename(const string& filestr) {
	DEBUG_MSG("Opening file : '" << filestr << "'");
	ifstream in(filestr);
	if (!in) {
		cout << "Unable to open the file '" << filestr << "'" << endl;
		exit(0);
	}

	string ref;
	uint32_t id;

	while (getline(in, ref)) {
		DEBUG_MSG("Getline from file :'" << filestr << "' = '" << ref << "'");

		if (ref.size() > 2 && exists_test(ref)) {
			#pragma omp critical (genome_numbers)
			{
				id = genome_numbers;
				DEBUG_MSG("Genome numbers: " << genome_numbers);
				genome_numbers++;
				filenames.push_back(ref);
			}
			DEBUG_MSG("Adding file :'" << ref << "'");
			insert_file(ref, id);
			DEBUG_MSG("File: '" << ref << "' added");
		}

		ref.clear();
	}

}


/** 
 * @brief Inserts a file into the index.
 * 
 * This function reads the file, generates a sketch of each line, and inserts the sketch into the index.
 * 
 * @param filestr The name of the file to insert into the index.
 * @param identifier A unique identifier for the file. 
 */
void Index::insert_file(const string& filestr,uint32_t identifier) {
	char type=get_data_type(filestr);
	zstr::ifstream in(filestr);
	string ref;
	vector<uint64_t> kmer_sketch;
	vector<uint64_t> sketch(F,-1);
	while(not in.eof()) {
		Biogetline(&in,ref,type);
		if(ref.size()>K) {
			compute_sketch(ref,sketch);
		}
		ref.clear();
	}
	file_sketches[filestr] = sketch;
	insert_sketch(sketch,identifier);
}


/** 
 * @brief Query the index with the contents of a file.
 * 
 * This function reads the file, generates a sketch of each line, and queries the index with each sketch.
 * The results of each query are written to the output file and also returned as a vector of integers.
 * 
 * @param filestr The name of the file to query.
 * @return A vector of integers representing the results of the queries.
 */
void Index::query_file(const string& filestr) {
	cout<<"Query:"<<filestr<<endl;
	char type=get_data_type(filestr);
	zstr::ifstream in(filestr);
	string ref;
	vector<uint64_t> kmer_sketch;
	vector<uint64_t> sketch(F,-1);
	while(not in.eof()) {
		Biogetline(&in,ref,type);
		if(ref.size()>K) {
			compute_sketch(ref,sketch);
		}
		ref.clear();
	}
	auto result(query_sketch(sketch));
	for(uint i(0);i<result.size();++i){
		*outfile<<result[i]<<" ";
		cout<<result[i]<<" ";
	}
	*outfile<<endl;
	cout<<endl;
	return;
}

void Index::query_file_of_file(const string& filestr) {
	cout<<"Query fof:"<<filestr<<endl;
	zstr::ifstream in(filestr);
	string ref;
	while(not in.eof()) {
		getline(in,ref);
		if(ref.size()>3) {
			query_file(ref);
		}
		ref.clear();
	}
	return;
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


/** 
 * @brief Reads a line from the file depending on the data type.
 *
 * This function reads a line or a group of lines from the given file.
 * For fastq files (type 'Q'), it discards three lines and reads the fourth line.
 * For fasta files (type 'A'), it reads lines until it reaches a line starting with '>'.
 * 
 * @param in Pointer to the input file stream.
 * @param result The line or group of lines read from the file.
 * @param type The type of file (either 'Q' for fastq or 'A' for fasta).
 */
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


/** 
 * @brief Determines the data type of a file based on its extension.
 *
 * This function checks the file extension and returns 'Q' for fastq files and 'A' for fasta files.
 * 
 * @param filename The name of the file.
 * @return The type of file ('Q' for fastq or 'A' for fasta).
 */
char Index::get_data_type(const string& filename)const{
	if(filename.find(".fq")!=string::npos){
		return 'Q';
	}
	if(filename.find(".fastq")!=string::npos){
		return 'Q';
	}
	return 'A';
}


/** 
 * @brief Computes the reverse complement of a k-mer.
 *
 * This function calculates the reverse complement of a k-mer by reversing the k-mer and switching each nucleotide with its complement.
 * 
 * @param min The k-mer to reverse complement.
 * @return The reverse complement of the k-mer.
 */
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


/** 
 * @brief Converts a nucleotide to an integer.
 *
 * This function maps each nucleotide to an integer as follows: A->0, C->1, G->2, T->3. 
 * If the input character does not represent a valid nucleotide, the function will return 0.
 * 
 * @param c The nucleotide to convert.
 * @return The integer representation of the nucleotide.
 */
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


/** 
 * @brief Converts a nucleotide to an integer, considering the reverse complement.
 *
 * This function maps each nucleotide to an integer as follows: A->3, C->2, G->1, T->0. 
 * If the input character does not represent a valid nucleotide, the function will return 0.
 * 
 * @param c The nucleotide to convert.
 * @return The integer representation of the reverse complement of the nucleotide.
 */
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


/** 
 * @brief Converts a string of nucleotides to a k-mer.
 *
 * This function converts a string of nucleotides to a k-mer by shifting the current k-mer 
 * and adding the integer representation of the next nucleotide in the string.
 * 
 * @param str The string of nucleotides to convert.
 * @return The k-mer representation of the string.
 */
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


/** 
 * @brief Reverses the hash of a 64-bit integer.
 *
 * This function reverses the hash of a 64-bit integer by applying a reverse operation.
 * 
 * @param x The 64-bit integer to reverse the hash of.
 * @return The reversed hash.
 */
uint64_t Index::revhash64 ( uint64_t x ) const {
	x = ((x >> 32) ^ x) * 0xD6E8FEB86659FD93;
	x = ((x >> 32) ^ x) * 0xD6E8FEB86659FD93;
	x = ((x >> 32) ^ x);
	return x;
}


/** 
 * @brief Unreverses the hash of a 64-bit integer.
 *
 * This function unreverses the hash of a 64-bit integer by applying an unreversing operation.
 * 
 * @param x The 64-bit integer to unreverse the hash of.
 * @return The unreversed hash.
 */
uint64_t Index::unrevhash64(uint64_t x) const{
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x);
	return x;
}


/** 
 * @brief Updates a k-mer by adding a new nucleotide.
 *
 * This function updates a k-mer by shifting the k-mer, adding the integer representation of 
 * the new nucleotide, and applying a modulus operation with 'offsetUpdatekmer'.
 * 
 * @param min The k-mer to update.
 * @param nuc The new nucleotide to add to the k-mer.
 */
void Index::update_kmer(kmer& min, char nuc)const {
	min<<=2;
	min+=nuc2int(nuc);
	min%=offsetUpdatekmer;
}


/** 
 * @brief Updates a k-mer for reverse complement string by removing a nucleotide and adding a new nucleotide at the beginning.
 *
 * This function updates a k-mer by shifting the k-mer, adding the integer representation of 
 * the reverse complement of the new nucleotide at the beginning.
 * 
 * @param min The k-mer to update.
 * @param nuc The new nucleotide to add to the k-mer.
 */
void Index::update_kmer_RC(kmer& min, char nuc)const {
	min>>=2;
	min+=(nuc2intrc(nuc)<<(2*K-2));
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
uint64_t Index::hash_family(const uint64_t x, const uint factor)const{
	return unrevhash64(x)+factor*revhash64(x);
}


/** 
 * @brief Densifies a sparse MinHash sketch.
 *
 * This function fills the empty cells of a MinHash sketch vector by iteratively
 * applying a hash function and replacing empty cells with hashed values.
 * 
 * @param sketch The MinHash sketch to densify.
 * @param empty_cell The number of empty cells in the sketch.
 */
void Index::sketch_densification(vector<uint64_t>& sketch, uint empty_cell) const {
	uint step(0);
	uint size(sketch.size());
	while(empty_cell!=0){
		for(uint i(0);i<size;i++){
			if(sketch[i]!=mi){
				uint64_t hash=hash_family(sketch[i],step)%size;
				if(sketch[hash]==mi){
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


/** 
 * @brief Calculates a perfect fingerprint from a hashed value.
 *
 * This function calculates a perfect fingerprint, which is an integer representation of the density 
 * of the hashed value in the hash space, by applying a transformation to the hashed value.
 * 
 * @param hashed The hashed value.
 * @return The perfect fingerprint of the hashed value.
 */
uint64_t Index::get_perfect_fingerprint(uint64_t hashed)const {
	uint64_t B(hashed);
	DEBUG_MSG("B    = '"<<B<<"'");
	B<<=1;
	B>>=1;
	DEBUG_MSG("B>>1 = '"<<B<<"'");
	uint64_t twopowern(1);
	twopowern<<=63;
	long double frac=((long double)(twopowern-B))/twopowern;
	DEBUG_MSG("-----------------------------------------------------------");
	DEBUG_MSG("Frac value : '"<<frac<<"'");
	// x = taille du sketch/taille du génome
	frac=pow(frac,E/F);
	DEBUG_MSG("E/F = '"<<E<<"/"<<F<<" = "<<E/F<<"'");
	DEBUG_MSG("New frac value : '"<<frac<<"'");
	frac=1-frac;
	DEBUG_MSG("Fingerprint range : '"<<fingerprint_range<<"'");
	DEBUG_MSG("Fingerprint range * frac =  '"<<fingerprint_range*frac<<"'");
	DEBUG_MSG("-----------------------------------------------------------");
	return fingerprint_range*frac;
}


/**
 * \brief Computes a minhash sketch for a given reference string.
 *
 * This method calculates the minhash sketch of the provided reference string and stores it in the provided sketch vector.
 * The sketch is computed by generating all the kmers in the reference, hashing them, and storing the smallest hash in each bucket.
 * The sketch is then densified to remove any empty cells.
 *
 * \param reference The reference string to compute the sketch for.
 * \param sketch A reference to a vector of uint64_t that will hold the computed sketch. 
 */
void Index::compute_sketch(const string& reference, vector<uint64_t>& sketch) const {
	if(sketch.size()!=F) {
		sketch.resize(F,-1);
	}
	DEBUG_MSG("------------------SKETCH--------------------");
	uint empty_cell(F);
	kmer S_kmer(str2numstrand(reference.substr(0,K-1)));//get the first kmer (k-1 bases)
	kmer RC_kmer(rcb(S_kmer));//The reverse complement
	for(uint i(0);i+K<reference.size();++i) {// all kmer in the genome
		Index::update_kmer(S_kmer,reference[i+K-1]);
		Index::update_kmer_RC(RC_kmer,reference[i+K-1]);
		kmer canon(min(S_kmer,RC_kmer));//Kmer min, the one selected
		uint64_t hashed=revhash64(canon);
		uint64_t bucket_id(unrevhash64(canon)>>(64-lF));//Which Bucket
		uint64_t fp=hashed;
		//MINHASH
		if(sketch[bucket_id]==mi){
			empty_cell--;
			sketch[bucket_id]=fp;
		}else if(sketch[bucket_id] > fp) {
			sketch[bucket_id]=fp;
		}
	}
	sketch_densification(sketch,empty_cell);
	for(uint64_t i=0; i <sketch.size();++i){
		//~ cout<<"avant";
		DEBUG_MSG("SKETCH avant : '"<<sketch[i]<<"'");
		//~ cout<<sketch[i]<<" ";
		sketch[i]=get_perfect_fingerprint(sketch[i]);
		//~ cout<<"Apres"<<" ";
		//~ cout<<sketch[i]<<endl;
		DEBUG_MSG("SKETCH apres : '"<<sketch[i]<<"'");
	}
	DEBUG_MSG("------------------SKETCH END--------------------");
}


/**
 * \brief Inserts a minhash sketch into the Buckets structure.
 *
 * This method inserts the given minhash sketch into the Buckets structure.
 * For each non-negative and within range fingerprint in the sketch, 
 * it locks the corresponding bucket, adds the genome_id to the bucket, 
 * and then unlocks the bucket.
 *
 * \param sketch The minhash sketch to insert into the Buckets structure.
 * \param genome_id The id of the genome that the sketch corresponds to.
 */
void Index::insert_sketch(const vector<uint64_t>& sketch,uint32_t genome_id) {
	for(uint i(0);i<F;++i) {
		if(sketch[i]<fingerprint_range and sketch[i]>=0){
			omp_set_lock(&lock[(sketch[i]+i*fingerprint_range)%mutex_number]);
			Buckets[sketch[i]].push_back(genome_id);
			Buckets_pos[sketch[i]].push_back(i);
			omp_unset_lock(&lock[(sketch[i]+i*fingerprint_range)%mutex_number]);
		}
	}
}


/**
 * \brief Queries the Buckets structure with a minhash sketch.
 *
 * This method queries the Buckets structure with the given minhash sketch and returns a result vector.
 * For each non-negative and within range fingerprint in the sketch, 
 * it increments the corresponding index in the result vector for each matching genome in the bucket.
 *
 * \param sketch The minhash sketch to query the Buckets structure with.
 * \return A vector of uint32_t that contains the query results. 
 *         Each index in the vector corresponds to a genome, 
 *         and the value at that index is the number of matching minhashes for that genome.
 */
vector<uint32_t> Index::query_sketch(const vector<uint64_t>& sketch) const{
	vector<uint32_t> result(genome_numbers,0);
	for(uint i(0);i<F;++i) {
		if(sketch[i]<fingerprint_range and sketch[i]>=0) {
			for(uint j(0);j<Buckets[sketch[i]].size();++j){
				if(Buckets_pos[sketch[i]][j]==i){
					result[Buckets[sketch[i]][j]]++;
				}
			}
		}
	}
	return result;
}


/**
 * \brief Merges two sketches into one.
 *
 * This method merges the given sketches into the first one.
 * The merge operation is performed by replacing each value in the first sketch 
 * with the minimum of the value in the first sketch and the corresponding value in the second sketch.
 *
 * \param sketch1 The first sketch that is to be merged. This will also hold the final merged sketch.
 * \param sketch2 The second sketch that is to be merged.
 */
void Index::merge_sketch( vector<int32_t>& sketch1,const vector<int32_t>& sketch2) const {
	for(uint i(0);i<sketch1.size();++i) {
		sketch1[i]=min(sketch1[i],sketch2[i]);
	}
}


/**
 * \brief Writes the current index to a binary file.
 *
 * This method writes the current index to a binary file. 
 * The file includes information about the index (like lF, K, E, W, and genome_numbers), 
 * the buckets in the index, and the filenames of the genomes in the index.
 * This enables the index to be saved and loaded from disk, 
 * so that it does not need to be computed every time the program runs.
 *
 * \param filestr The name of the binary file to write to.
 */
void Index::dump_index_disk(const string& filestr)const{
	zstr::ofstream out(filestr,ios::binary);
	out.write(reinterpret_cast<const char*>(&lF), sizeof(lF));
	out.write(reinterpret_cast<const char*>(&K), sizeof(K));
	out.write(reinterpret_cast<const char*>(&E), sizeof(E));
	out.write(reinterpret_cast<const char*>(&W), sizeof(W));
	out.write(reinterpret_cast<const char*>(&genome_numbers), sizeof(genome_numbers));
	for(uint i(0);i<fingerprint_range*F;++i){
		uint64_t size(Buckets[i].size());
		out.write(reinterpret_cast<const char*>(&size), sizeof(size));
		out.write(reinterpret_cast<const char*>(&(Buckets[i][0])), size*sizeof(gid));
		out.write(reinterpret_cast<const char*>(&(Buckets_pos[i][0])), size*sizeof(uint16_t));
	}
	for(uint i(0);i<filenames.size();++i){
		out.write((filenames[i]+'\n').c_str(),filenames[i].size()+1);
	}
}




/**
 * \brief Calculates and prints a Jaccard distance matrix for the indexed genomes.
 *
 * This function computes a pairwise Jaccard distance for all indexed genomes and
 * prints the resulting distance matrix. The Jaccard distance is calculated by computing
 * the Jaccard index (intersection over union) between the sketches of two genomes,
 * and then subtracting the result from 1. A Jaccard distance of 0 between two genomes
 * implies that they are identical, while a Jaccard distance of 1 implies that they share no sketches in common.
 * The resulting distance matrix is symmetrical, with a diagonal of zeroes.
 *
 * The output is formatted as a tab-delimited matrix with genome filenames as row and column labels.
 * The Jaccard distances are printed with a fixed precision of three decimal places.
 * For a genome compared with itself, a "-" character is printed instead of the distance.
 *
 * Throws a runtime_error if a sketch for a filename cannot be found in the index.
 */
void Index::print_matrix() const {
    DEBUG_MSG("PRINT MATRIX: ");
    int size = genome_numbers;
    vector<vector<double>> matrix(size, vector<double>(size)); // Changer le type de la matrice en double pour stocker la distance de Jaccard.

    for(int i = 0; i < size; ++i){
        DEBUG_MSG("i:'" << i << "' new");
        auto it_i = file_sketches.find(filenames[i]);
        if (it_i == file_sketches.end()) {
            throw std::runtime_error("Key not found");
        }
        const vector<uint64_t>& sketch_i = it_i->second;
        for(int j = i; j < size; ++j){
            DEBUG_MSG("j:'" << j << "' new");
            if(i == j){
                matrix[i][j] = 0.0; // La distance de Jaccard d'un génome à lui-même est 0.
            } else {
                auto it_j = file_sketches.find(filenames[j]);
                if (it_j == file_sketches.end()) {
                    throw std::runtime_error("Key not found");
                }
                const vector<uint64_t>& sketch_j = it_j->second;

                // Calculez l'intersection et l'union des deux sketches.
                unordered_set<uint64_t> intersection, union_set;
                for (auto& sketch : sketch_i) {
                    if (find(sketch_j.begin(), sketch_j.end(), sketch) != sketch_j.end())
                        intersection.insert(sketch);
                    union_set.insert(sketch);
                }
                for (auto& sketch : sketch_j) {
                    union_set.insert(sketch);
                }

                // Calculez la distance de Jaccard et stockez-la dans la matrice.
                double jaccard_distance = 1.0 - static_cast<double>(intersection.size()) / union_set.size();
                matrix[i][j] = jaccard_distance;
                matrix[j][i] = jaccard_distance;
            }
        }
    }

    cout << "##Names ";
    for (int i = 0; i < size; ++i) {
        cout << filenames[i] << "\t";
    }
    cout << endl;

    for(int i = 0; i < size; ++i) {
        cout << filenames[i] << "\t";
        for(int j = 0; j < size; ++j) {
            if(i == j)
                cout << "-\t";
            else
                cout << fixed << setprecision(3) << matrix[i][j] << "\t"; // Affichez la distance de Jaccard avec 3 chiffres après la virgule.
        }
        cout << endl;
    }
}

