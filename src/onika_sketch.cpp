#include "strict_fstream.hpp"
#include "zstr.hpp"
#include "onika_index.h"
#include "common.h"



using namespace std;
const int bufferSize = 10000;



Sketch::Sketch(uint32_t ilF=10, uint32_t iK=31,uint32_t iW=8,uint32_t iH=4, const string ifilename="onikaOutput.gz") {
}




Sketch::~Sketch() {
	//delete[] Buckets;
}


void Sketch::fasta_sketch(const string& filename) {

}

void Sketch::insert_sketch(void) {

}

void Sketch::query_sketch(Info sketch) {
}

