#include "onika_index.h"
#include "onika_sketch.h"
#include "optionparser.h"
#include "common.h"

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include <fcntl.h> //open
#include <libgen.h> // dirname
#include <iomanip> // std::setw
#include <zlib.h>


using namespace std;
using namespace chrono;

using option::Option;
using option::Descriptor;
using option::Parser;
using option::Stats;
using option::ArgStatus;



struct Arg: public option::Arg {
	static void printError(const char* msg1, const option::Option& opt, const char* msg2) {
		fprintf(stderr, "%s", msg1);
		fwrite(opt.name, opt.namelen, 1, stderr);
		fprintf(stderr, "%s", msg2);
	}
	static option::ArgStatus Unknown(const option::Option& option, bool msg) {
		if (msg) printError("Unknown option '", option, "'\n");
		return option::ARG_ILLEGAL;
	}
	static option::ArgStatus NonEmpty(const option::Option& option, bool msg) {
		if (option.arg != 0 && option.arg[0] != 0)
			return option::ARG_OK;
		if (msg) printError("Option '", option, "' requires a non-empty argument\n");
		return option::ARG_ILLEGAL;
	}
	static option::ArgStatus Numeric(const option::Option& option, bool msg) {
		char* endptr = 0;
		if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
		if (endptr != option.arg && *endptr == 0)
			return option::ARG_OK;
		if (msg) printError("Option '", option, "' requires a numeric argument\n");
		return option::ARG_ILLEGAL;
	}
};



enum  optionIndex {
	UNKNOWN,
	LIST,
	DIST,
	QUERY,
	OUTPUT,
	DUMP,
	LOAD,
	KMER,
	FETCH,
	EGS,
	WORD,
	MIN,
	LOGO,
	HELP
};



const option::Descriptor usage[] = {
	{UNKNOWN, 0,"" , "" , Arg::Unknown,"\n***Input***"},
	{LIST, 0, "I" , "index" ,Arg::NonEmpty,
		"  --index, -I <filename> "
			"\tInput file of files to Index.\v"
	},

	{QUERY, 0, "Q", "query"    , Arg::NonEmpty,
		"  --query, -Q <filename> "
			"\tInput file of file to Query.\v"
	},
	{UNKNOWN, 0,"" , "" , Arg::Unknown,"\n***Main parameters***"},
	{OUTPUT, 0, "O", "output", Arg::NonEmpty,
		"  --output, -O <filename> "
			"\tOutput file (onikaOutput.gz)"
	},

	{UNKNOWN, 0,"" , "" , Arg::Unknown,"\n***Main parameters***"},
	{KMER,  0, "K" , "kmer"  ,Arg::Numeric,
		"  --kmer, -K <int> "
			"\tKmer size (31).\v"
	},
	{FETCH,  0, "S" , "sketch"  ,Arg::Numeric,
		"  --sketch, -S <int> "
			"\tSet sketch size to 2^S (15).\v"
	},
	{UNKNOWN, 0,"" , "" , Arg::Unknown,"\n***Advanced parameters*** (You know what you are doing)"},
	{WORD,  0, "W" , "word"  ,Arg::Numeric,
		"  --word, -W <int> "
			"\tFingerprint size (12). Modify with caution, \
			larger fingerprints enable queries with less false positive but increase EXPONENTIALY\
			the overhead as the index count S*2^W cells. \v"
	},
	{EGS,  0, "E" , "EGS"  ,Arg::Numeric,
		"  --EGS, -E <int> "
			"\t Expected gemome size (5 000 000).  Modify with caution.\v"
	},
	{UNKNOWN, 0,"" , "" , Arg::Unknown,"\n***Index files***"},
	{LOAD, 0, "L", "load", Arg::NonEmpty,
		"  --load, -L <filename> "
			"\tLoad an index to the given file."
	},
	{DUMP, 0, "D", "dump", Arg::NonEmpty,
		"  --dump, -D <filename> "
			"\tDump the current index to the given file."
	},
	{DIST, 0, "", "dist", Arg::None,
		"  --dist "
			"\tGenerate a full, asymmetric distance matrix."
	},


	{LOGO, 0, "",  "logo", Arg::None,
		"  --logo "
			"\tPrint ASCII art logo, then exit."
	},

	{HELP,  0, "h" , "help"  ,Arg::None,
		"  --help, -h  "
			"\tPrint usage and exit." },
	{0,0,0,0,0,0}
};



option::Option *options = NULL, *buffer = NULL;
void deleteOptsArrays() {
	if (options) {
		delete [] options;
	}
	if (buffer) {
		delete [] buffer;
	}
}



static int old_wd;
static void changeDirFromFilename(const char* fname) {
	DEBUG_MSG("CWD is " << getcwd(NULL, 0));
	old_wd = open(".", O_CLOEXEC);
	char fname_copy[PATH_MAX];
	strncpy(fname_copy, fname, PATH_MAX);
	fname_copy[PATH_MAX - 1] = '\0';
	char *dname = dirname(fname_copy);
	DEBUG_MSG("Changing to directory " << dname);
	errno = 0;
	if (chdir(dname)) {
		cout << "Error: " << strerror(errno) << endl;
	}
	DEBUG_MSG("Now CWD is " << getcwd(NULL, 0));
}



static void restoreDir() {
	errno = 0;
	if (fchdir(old_wd)) {
		cout << "Error: " << strerror(errno) << endl;
	}
	DEBUG_MSG("Restore working directory to " << getcwd(NULL, 0));
}



int main(int argc, char * argv[]){
	int F=0,K=0,W=0,E=0;
	string list_file = "";
	string query_file = "";
	string out_file = "";
	argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
	option::Stats stats(true, usage, argc, argv);
	options = new option::Option[stats.options_max];
	buffer = new option::Option[stats.buffer_max];
	atexit(deleteOptsArrays);
	option::Parser parse(usage, argc, argv, options, buffer);

	if (parse.error()) {
		cout << "Bad usage!!!" << endl;
		delete[] options;
		delete[] buffer;
		return EXIT_FAILURE;
	}


	/**********************************/
	/* Check Help             options */
	/**********************************/
	if (options[HELP] || argc == 0) {
		option::printUsage(clog, usage,160);
		delete[] options;
		delete[] buffer;
		return EXIT_SUCCESS;
	}

	/***********************************/
	/* Set the k-mer length and Other  */
	/***********************************/
	K = options[KMER] ? atoi(options[KMER].last()->arg) : 31;
	DEBUG_MSG("K = " << K);
	F = options[FETCH] ? atoi(options[FETCH].last()->arg) : 15;
	DEBUG_MSG("F = " << F);
	E = options[EGS] ? atoi(options[EGS].last()->arg) : 5000000;
	DEBUG_MSG("E = " << E);
	W = options[WORD] ? atoi(options[WORD].last()->arg) : 12;
	DEBUG_MSG("W = " << W);

	/************************************/
	/* Complain about unknown arguments */
	/************************************/
	for (int i = 0; i < parse.nonOptionsCount(); ++i) {
		cout << "Non-option argument #" << i << " is " << parse.nonOption(i)<<endl;
		cout << "Ignoring unknown argument '" << parse.nonOption(i) << "'" << endl;
	}
	if (parse.nonOptionsCount()) {
		cout << "Bad usage!!!" << endl;
		return EXIT_FAILURE;
	}
	if (options[OUTPUT]) {
		out_file = options[OUTPUT] ? (options[OUTPUT].last()->arg) : "onikaOutput.gz";
		DEBUG_MSG("Output file name = " << out_file);
	} else {
		out_file="onikaOutput.gz";
	}


	/**********************************************/
	/* Display the ASCII art logo of the program. */
	/**********************************************/
	if (options[LOGO]) {
		string logo_name="../resources/onika.ascii";
		ifstream logo;
		string line;
		logo.open(logo_name);
		if (logo.is_open()){
			while ( getline (logo,line) ){
				cout << line << '\n';
			}
			logo.close();
		}
		else cout << "Unable to open file :'"<<logo_name<<"'"<<endl;
	}



	cout << "+-------------------------------------------------------------------+" << endl;
	cout << "|                            Informations                           |" << endl;
	cout << "+-----------------------------------+-------------------------------+" << endl;
	
	Index* monindex;

	if(options[LOAD]){
		string indexfile = options[LOAD].last()->arg;
		monindex = new Index(indexfile,out_file);
	}else{
		monindex = new Index(F,K,W,E,out_file);

	}
	time_point<system_clock> start, endindex,end;
	start = std::chrono::system_clock::now();

	/*****************************************/
	/* Add the genomes given in config files */
	/*****************************************/
	if (options[LIST]) {
		//option::Option* opt = options[LIST];
		list_file = options[LIST].last()->arg;
		ifstream ifs(list_file);
		if (!ifs) {
			cout << "Unable to open the file '" << list_file << "'" << endl;
		}
		changeDirFromFilename(list_file.c_str());
		DEBUG_MSG("Opening file : '"<<list_file<<"'");
		//monindex->get_filename(list_file.substr(list_file.find_last_of("/\\") + 1));
		monindex->get_filename(list_file);
		DEBUG_MSG("File added");
		restoreDir();

	}
	
	if(options[DUMP]) {
		string indexfile = options[DUMP].last()->arg;
		monindex->dump_index_disk(indexfile);
	}

	if(options[DIST]) {
		monindex->print_matrix();
	}

	endindex = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = endindex - start;
	cout << "| Indexing lasted (s)               |" << setw(30) << setfill(' ') << elapsed_seconds.count() << " |" << endl;


	/*****************************************/
	/* Add the query file and do the request */
	/*****************************************/


	if (options[QUERY]) {
		query_file = options[QUERY].last()->arg;
		ifstream ifs(query_file);
		if (!ifs) {
			cout << "Unable to open the file '" << query_file << "'" << endl;
		}
		DEBUG_MSG("Opening file...");
		monindex->query_file(query_file);
		DEBUG_MSG("Query done.");
	}

	DEBUG_MSG("Output the index to :'"<<out_file<<"'");

	end = std::chrono::system_clock::now();
	elapsed_seconds = end - endindex;
	cout << "| Query lasted (s)                  |" << setw(30) << setfill(' ') << elapsed_seconds.count() << " |" << endl;
	elapsed_seconds = end - start;
	cout << "| Whole run lasted (s)              |" << setw(30) << setfill(' ') << elapsed_seconds.count() << " |" << endl;
	cout << "+-----------------------------------+-------------------------------+" << endl;
	cout << "| k-mer size                        |" << setw(30) << setfill(' ') << K << " |" << endl
		<< "| S                                 |" << setw(30) << setfill(' ') << F << " |" << endl
		<< "| Number of fingerprints            |" << setw(30) << setfill(' ') << monindex->F<< " |" << endl
		<< "| W                                 |" << setw(30) << setfill(' ') << W << " |" << endl
		<< "| E                                 |" << setw(30) << setfill(' ') << E << " |" << endl
		<< "| Number of indexed genomes         |" << setw(30) << setfill(' ') << monindex->getNbGenomes() << " |" << endl;
	
	if(monindex->outfile->is_open()){
		monindex->outfile->close();
	} else {
		std::cerr << "Failed to close index.\n";
	}
	//delete monindex;

	cout << "+-------------------------------------------------------------------+" << endl;
	cout << "|                                  Done                             |" << endl;
	cout << "+-----------------------------------+-------------------------------+" << endl;
	return EXIT_SUCCESS;
}
