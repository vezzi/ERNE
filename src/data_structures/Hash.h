/*
 *  This file is part of ERNE.
 *  Copyright (c) 2011 by Cristian Del Fabbro <delfabbro@appliedgenomics.org>,
 *  Francesco Vezzi <vezzi@appliedgenomics.org>,
 *  Alexandru Tomescu <alexandru.tomescu@uniud.it>
 *  Nicola Prezza <nicolapr@gmail.com>, and
 *  Alberto Policriti <policriti@uniud.it>
 *
 *   ERNE is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   ERNE is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.

 *   You should have received a copy of the GNU General Public License
 *   along with ERNE.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef HASH_H_
#define HASH_H_

#define _HASH_VERSION 0

#include <vector>
using namespace std;

/*
#include "errors/Data_Exception.h"
#include "errors/Data_Not_Found.h"
#include "errors/Data_Overflow.h"
*/
#include "errors/Generic_Exception.h"
#include "errors/Incorrect_Format.h"
using namespace errors;

#include "common.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "io/Fasta.h"
using namespace fasta;

#include <boost/lexical_cast.hpp>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"

#include <iostream>
#include <fstream>

#include "SetZ.h"
#include "GlobalToLocal.h"
#include "ResultItem.h"
#include "ResultItemGapped.h"

#ifdef HAVE_MPI
#include "mpi.h"

#include "data_structures/Mask.h"

#include <boost/thread/thread.hpp>
using namespace boost;
#endif

typedef vector<ResultItem> Items;
typedef vector<ResultItemGapped> ItemsGapped;

class Hash {
public:

	Hash();
	Hash(int k);
	Hash(int k, int blockLength,bool methyl_hash);
	~Hash();
	void free_hash_memory();

	void createHASH(const vector<string> & fastaFiles) throw (Generic_Exception);

	static int get_filesize(const std::string & path);
	static int calculate_k(const string & filename, bool methyl_hash);
	static int calculate_k(const vector<string> &filenames, bool methyl_hash);

	unsigned long int returnQ() {return this->q;}

	void save(const char * filename) throw (File_Not_Found);
	void load(const char * filename) throw (File_Not_Found);

	inline unsigned long int fill_right(unsigned long int Nq, char base);
	inline unsigned long int roll_right(char first_digit, unsigned long int Nq, char base);
	unsigned long int roll_right_XOR(unsigned long int oldVal, char firstCh, char lastCh, int m);

	unsigned long int hOR(int* P, int bl,int pp);
	unsigned long int hXOR(int* P, int bl);

	void search(const string & read, Items & output, unsigned  int max_errors);
	void search(const string & read, Items & output, unsigned int max_errors, unsigned int & delta, vector<pair<int, int> > & DeltaSolutions, bool indels);
	void search_gapped(const string & read, ItemsGapped & output, t_pattern_length seed_size, int seed_errors, int max_errors, t_length max_gap);
	vector<pair<char, unsigned int> > alignSW(const char * read, int patternLength, unsigned int globalPos);

	void set_indels_max_value(int value) { indels = value; }

	int max( int f1, int f2, int f3, int & ptr );

	int k;
	unsigned long int q;
	unsigned int blockLength;
	unsigned long int h;

	char *TEXT; // pointer to an array that stores the text
	unsigned int textLength; // length of the text
	//int Zlength;
	SetZ *Z;
	//int *errorsInZ;

	unsigned int *HASHvalues;
	unsigned int *HASHcounter;
	unsigned int LAST;
	unsigned int MASK;

	bool methyl_hash;//if true use only 1 bit per base, otherwise use 2 bit per base

	GlobalToLocal globalToLocal;

	typedef boost::mt19937 RNGType;
	RNGType rng;

protected:
	size_t limit;

	int match;
	int mismatch;
	int indel;
	int indels;
	int SM[5][5];

	static bool search_gapped_item_sort(const ResultItem &left, const ResultItem &right);
	int stranded_gap_search_subroutine(ItemsGapped & output, Items & left, Items & right, t_length max_gap, const char * s, int left_characters, int right_characters, int max_errors);

	int alignSW_Improved( const char *read, int patternLength, unsigned int globalPos, int errors, int & editOperations, bool & WithIndel, int BestAchivableScore);
};




#endif /* HASH_H_ */
