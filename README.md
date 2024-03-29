# ONIKA
ONIKA stand for Optimal Naive Inverted index for _K_-mer Analysis.

TODO

# Table of contents
1. [Command line options](#cmo)
2. [Installation](#install)
3. [Sample usages](#su)
4. [Credits](#cr)

## Command line options <a name="cmo"></a>

### 1. Indexing  :

ONIKA can handle FASTA, MultiLine FASTA  and FASTQ files gziped or not.

`--index`, `-I <file>`             Input file of files to Index. 

Such file of files should contain one file name or path per line. Each  file will be sketched and inserted in the index as a single entry.

Example:

```bash
    onika --index file_of_file.txt
```
    
        
### 2. Query  :

`--query,`, `-Q <file>`             Input file of files to Query.

Example:

```bash
    onika --query my_genomes.txt
```

`--querylines,`, `-q <file>`             Input file where each line is a separate entry to Query.

Example:

```bash
    onika --indexlines my_genomes.fa
```


### 3. Dist Option

The `--dist` option allows you to generate a full, asymmetric distance matrix, similar to the output provided by DASH. This can be particularly useful for visualizing the distances between different elements in your dataset. 

Usage:
```bash
onika --index file_of_file.txt --dist --indexinfo
```

### 4. Others

`--logo`                           Print ASCII art logo, then exit.

`--help`, `-h`                    Print usage and exit.


## Installation <a name="install"></a>

### Requirements

ONIKA requires:

* A modern, C++17 ready compiler such as `g++` version 4.9 or higher or `clang` version 3.2 or higher.
* A 64-bit operating system. Either Mac OS X or Linux are currently supported.
* `zlib` to be already installed on your system (on Debian/Ubuntu it corresponds to the packages `zlib1g-dev`).

### Single user installation

To download and install `ONIKA` into some
user local directory (e.g., `${HOME}/local_install`) , use the
following commands:


#### ONIKA

```sh
git clone https://github.com/Malfoy/ONIKA.git ;

cd ONIKA ;

./configure --prefix=$HOME/local_install ;

make ;

make install ;
```


#### Uninstall

To remove `ONIKA`from your system use the following command:

```sh
cd ONIKA && make uninstall
```


## Sample Usages <a name="su"></a>

To test your installation you can go to the resources folder
```sh
cd resources
```


You can generate a distance matrix from the provided file of files.

```sh
onika --index file_of_file.txt --dist
```

By default ONIKA generate a hit list for each query file
 
 ```sh
onika --index file_of_file.txt --query file_of_file.txt --output hits.gz 
zcat hits.gz
```

You can change the sketch size depending of our usage (here 1,024 fingerprint per sketch)
 ```sh
onika --index file_of_file.txt --query file_of_file.txt --output hitsS10.gz --sketch 10
zcat hitsS10.gz
```
## Credits <a name="cr"></a>

TODO

Authors:
----------------

* Clément AGRET     <clement.agret@cyu.fr>
* Bastien CAZAUX    <bastien.cazaux@univ-lille.fr>
* Antoine LIMASSET  <antoine.limasset@univ-lille.fr>

Programmers:
-------------------------

* Clément AGRET     <clement.agret@cyu.fr>
* Antoine LIMASSET  <antoine.limasset@univ-lille.fr>
