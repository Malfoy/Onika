#include "strict_fstream.hpp"
#include "zstr.hpp"
#include "onika_sketch.h"
#include "common.h"



using namespace std;
const int bufferSize = 10000;


Sketch::Sketch(char* inString) {
	string = inString;
}


Sketch::~Sketch() {
	//delete[] Buckets;
}


void Sketch::fasta_sketch(const string& filename) {

}

void Sketch::insert_sketch(void) {

}

void Sketch::query_sketch(Info iSketch) {
}

