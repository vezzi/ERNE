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

#include "Hash.h"
#include <math.h>

//BOOST_CLASS_VERSION(Hash, _HASH_VERSION);

//static
int Hash::get_filesize(const std::string & path) {
	FILE *pFile = NULL;
	// get the file stream
	pFile = fopen(path.c_str(), "rb" );
	if (pFile == NULL) {
		std::cerr << "Can't open file " << path << " to acquire its length!" << std::endl;
		exit(3);
	}
	// set the file pointer to end of file
	fseek( pFile, 0, SEEK_END );
	// get the file size
	int size = ftell( pFile );
	// return the file pointer to begin of file if you want to read it
	fclose( pFile );
	return size;
}

// static
int Hash::calculate_k(const string & filename, bool methyl_hash) {
	int size = get_filesize(filename);

	int k;

	if(methyl_hash){
		k = ceil(log(size)/log(2));
		if (k < 20)
			k = 20;
		if (k > 30)
			k = 30;
	}else{
		k = ceil(log(size)/log(4));
		if (k < 10)
			k = 10;
		if (k > 15)
			k = 15;
	}

	return k;
}

//static
int Hash::calculate_k(const vector<string> &filenames, bool methyl_hash) {
	int size = 0;
	for (vector<string>::const_iterator iter = filenames.begin(); iter != filenames.end(); iter++)
		size += get_filesize(*iter);

	int k;

	if(methyl_hash){
		k = ceil(log(size)/log(2));
		if (k < 20)
			k = 20;
		if (k > 30)
			k = 30;
	}else{
		k = ceil(log(size)/log(4));
		if (k < 10)
			k = 10;
		if (k > 15)
			k = 15;
	}

	return k;
}


Hash::Hash() {
	TEXT = NULL;
	Z = NULL;
	HASHvalues = NULL;
	HASHcounter = NULL;
	match = 1;
	mismatch = -1;
	indel = 2;
	indels = 5;
	methyl_hash = false;
	SM[0][0] = match;    SM[0][1] = mismatch; SM[0][2] = mismatch; SM[0][3] = mismatch; SM[0][4] = mismatch;
	SM[1][0] = mismatch; SM[1][1] = match;    SM[1][2] = mismatch; SM[1][3] = mismatch; SM[1][4] = mismatch;
	SM[2][0] = mismatch; SM[2][1] = mismatch; SM[2][2] = match;    SM[2][3] = mismatch; SM[2][4] = mismatch;
	SM[3][0] = mismatch; SM[3][1] = mismatch; SM[3][2] = mismatch; SM[3][3] = match;    SM[3][4] = mismatch;
	SM[4][0] = mismatch; SM[4][1] = mismatch; SM[4][2] = mismatch; SM[4][3] = mismatch; SM[4][4] = match;
}

Hash::Hash(int k) {

	TEXT = NULL;
	Z = NULL;
	HASHvalues = NULL;
	HASHcounter = NULL;
	methyl_hash = false;

	this->k=k;
	this->q =(1 << this->k) - 1;
//	this->maskRollRight= (1 << (2 * this->k));
	this->blockLength=10;
	int tt=blockLength;
	this->h=1;
	while(tt>1) {
		this->h= (4*this->h)%q;
		tt--;
	}
	limit =   512*1024*1024;
	match = 1;
	mismatch = -1;
	indel = 2;
	indels = 5;
	SM[0][0] = match;    SM[0][1] = mismatch; SM[0][2] = mismatch; SM[0][3] = mismatch; SM[0][4] = mismatch;
	SM[1][0] = mismatch; SM[1][1] = match;    SM[1][2] = mismatch; SM[1][3] = mismatch; SM[1][4] = mismatch;
	SM[2][0] = mismatch; SM[2][1] = mismatch; SM[2][2] = match;    SM[2][3] = mismatch; SM[2][4] = mismatch;
	SM[3][0] = mismatch; SM[3][1] = mismatch; SM[3][2] = mismatch; SM[3][3] = match;    SM[3][4] = mismatch;
	SM[4][0] = mismatch; SM[4][1] = mismatch; SM[4][2] = mismatch; SM[4][3] = mismatch; SM[4][4] = match;
}

Hash::Hash(int k, int blockLength,bool methyl_hash) {

	this->methyl_hash = methyl_hash;
	TEXT = NULL;
	Z = NULL;
	HASHvalues = NULL;
	HASHcounter = NULL;

	this->k=k;
	this->q=0;
	int bit;

	if(methyl_hash)
		bit = 1;
	else
		bit=2;

	this->q =(1 << (this->k*bit)) - 1;

	for(int i = 0; i < k*bit; i++) {
		this->q |= 1 << i;
		//DEFAULT_CHANNEL << i <<" "<< this->q << "\n";
	}

	this->blockLength=blockLength;
	int tt=blockLength;
	this->h=1;
	while(tt>1) {
		this->h= (4*this->h)%q;
		tt--;
	}
	limit =   512*1024*1024;

	match = 1;
	mismatch = -1;
	indel = 2;
	indels = 5;
	SM[0][0] = match;    SM[0][1] = mismatch; SM[0][2] = mismatch; SM[0][3] = mismatch; SM[0][4] = mismatch;
	SM[1][0] = mismatch; SM[1][1] = match;    SM[1][2] = mismatch; SM[1][3] = mismatch; SM[1][4] = mismatch;
	SM[2][0] = mismatch; SM[2][1] = mismatch; SM[2][2] = match;    SM[2][3] = mismatch; SM[2][4] = mismatch;
	SM[3][0] = mismatch; SM[3][1] = mismatch; SM[3][2] = mismatch; SM[3][3] = match;    SM[3][4] = mismatch;
	SM[4][0] = mismatch; SM[4][1] = mismatch; SM[4][2] = mismatch; SM[4][3] = mismatch; SM[4][4] = match;
}


Hash::~Hash() {
	if (TEXT != NULL)
		delete [] TEXT;
	if (Z != NULL)
		Z->~SetZ();
	if (HASHvalues != NULL)
		delete [] HASHvalues;
	if (HASHcounter != NULL)
		delete [] HASHcounter;

}

void Hash::free_hash_memory() {

	if (Z != NULL){ Z->~SetZ(); Z = NULL;}
	if (HASHvalues != NULL){ delete [] HASHvalues; HASHvalues = NULL; }
	if (HASHcounter != NULL){ delete [] HASHcounter; HASHcounter = NULL; }

}

/*inline unsigned long int Hash::fill_right(unsigned long int Nq, char base) {
	return ((Nq << 2) + base) % q ;
}

inline unsigned long int Hash::roll_right( char first_digit, unsigned long int Nq, char base) {
	return  ((q + Nq - first_digit * h) * 4 + base) % q;

}*/

unsigned long int Hash::roll_right_XOR(unsigned long int oldVal, char firstCh, char lastCh, int m){

	int bit=2;

	if(methyl_hash)
		bit=1;
	else
		bit=2;

	int r = m % k;//size of the last block

	unsigned long int mask = (((unsigned long int)1<<(k*bit-1))-1)|((unsigned long int)1<<(k*bit-1));// k times '1'

	unsigned long int x,y;//temporary variables

	oldVal &= mask;

	x = ((oldVal<<bit) & mask)^( r==0?0:(((unsigned long int)lastCh)<<((k-r)*bit)) );
	y = (oldVal >> bit*(k-1)) ^ ((unsigned long int)firstCh) ^ (r==0?((unsigned long int)lastCh):0);

	return (x ^ y);

}//roll_right_XOR


void Hash::createHASH(const vector<string> & fasta_files) throw (Generic_Exception) {


	//pattern in the text -> 30 bit value (hXOR) -> HASHcounter
	//read block b -> hXOR(b) XOR (elements in Z[ hOR(b) ]) -> 30 bit value -> HASHcounter

	string header;
	char c;
	unsigned int size = 0;
	vector<unsigned int> start;
	vector<unsigned int> end;
	int numContig = 0;
	string line;

	for (vector<string>::const_iterator iter = fasta_files.begin(); iter != fasta_files.end(); iter++) {
		Auto_Unzip file(*iter);
		istream & fasta = file.filtered();

		Fasta contig;
		while(!fasta.eof()) {
			fasta >> contig;
			numContig++;
			size+=(contig.length()+1);
		}

	}



	TEXT = new char[size+1];

	//int actualBit=0;
	//int actualByte=0;

	if(methyl_hash)
		VERBOSE_CHANNEL << "building hash for methylation search"<<endl;
	else
		VERBOSE_CHANNEL << "building hash for standard search"<<endl;

	VERBOSE_CHANNEL << "k is " << k << "\n";
	VERBOSE_CHANNEL << "size of the TEXT is " << (size+1) << "\n";
	VERBOSE_CHANNEL << "number of contig is " << numContig << "\n";

	globalToLocal = GlobalToLocal(numContig);

	size = 0;
	int numberOfNs = 0;
	numContig = 0;


	for (vector<string>::const_iterator iter = fasta_files.begin(); iter != fasta_files.end(); iter++) {

		Auto_Unzip file(*iter);
		istream & fasta = file.filtered();


		while(not fasta.eof()) {

			Fasta contig;
			fasta >> contig;

			unsigned int contigSize =contig.length();
			globalToLocal.startPositions[numContig] = size;
			globalToLocal.endPositions[numContig] = size+contigSize-1;
			globalToLocal.contigsNames[numContig] = contig.get_id();

			start.push_back(size);
			end.push_back(size+contigSize);

			unsigned int i=0;
			for (string::const_iterator iter = contig.get_sequence().begin(); iter != contig.get_sequence().end(); iter++) {
				TEXT[size+i] = *iter;
				i++;
				if(*iter != 'A' and *iter != 'a' and *iter != 'C' and *iter != 'c' and *iter != 'G' and *iter != 'g' and *iter != 'T' and *iter != 't') {
						numberOfNs++;
				}
			}
			size = size + contigSize;
			TEXT[size] = '$';
			size++;
			numContig++;

		}
	}

	if (not globalToLocal.detectDuplicates()) {
		ERROR_CHANNEL << endl << "Please, check fasta input(s)!" << endl;
		exit(3);
	}

/*	for(int i=0; i< size; i++) {
		DEFAULT_CHANNEL << TEXT[i];
	}
	DEFAULT_CHANNEL << "\n";
*/
	//array declaration
	HASHcounter = new unsigned int[q];
	for(unsigned int i =0; i< q; i++) {
		HASHcounter[i] = 0;
	}
	int *Nnumbers = new int[numberOfNs];
	numberOfNs=0;

	unsigned long int Tq=0;
	int randNum;
	char lastChar[(int)blockLength];
	char character='A';
	int t = blockLength/k;//number of sub-blocks of size k in the block
	int r = blockLength % k;//size of the last sub-block in the block
	int subblocks,j;
	int bit;

	if(methyl_hash)
		bit = 1;
	else
		bit = 2;

	unsigned long int pw;

	for(int i =0; i< numContig; i++) {
		unsigned int s = start.at(i);
		unsigned int e = end.at(i);
		Tq =0;

		for(subblocks=0;subblocks<t;subblocks++){//for each sub-block(of size k) in the block
			pw = 0;
			for(j=0;j<k;j++){//for each char in a sub-block(of size k) in the block

				c = TEXT[s+subblocks*k + j];

				randNum= (int)(rand()%(bit==1?2:4));
				switch (c) {
					case 'a' : case 'A' : character=0; break;
					case 'c' : case 'C' : character=1; break;
					case 'g' : case 'G' : character=(bit==1?0:2);break;
					case 't' : case 'T' : character=(bit==1?1:3);break;
					default : character=randNum; Nnumbers[numberOfNs]=randNum; numberOfNs++; break;
				}

				lastChar[(s+subblocks*k + j)%blockLength]=character;
				pw = (pw << bit) ^ (unsigned long int)character;
			}//for
			Tq = Tq ^ pw;
		}//for all the sub-blocks(of size k) of the block

		pw = 0;
		for(j=0;j<r;j++){//for each char in the last block (only if there is a remainder r)

			c = TEXT[s+t*k + j];

			randNum= (int)(rand()%(bit==1?2:4));
			switch (c) {
				case 'a' : case 'A' : character=0; break;
				case 'c' : case 'C' : character=1; break;
				case 'g' : case 'G' : character=(bit==1?0:2);break;
				case 't' : case 'T' : character=(bit==1?1:3);break;
				default : character=randNum; Nnumbers[numberOfNs]=randNum; numberOfNs++; break;
			}

			lastChar[(s+t*k + j)%blockLength]=character;
			pw = (pw << bit) ^ (unsigned long int)character;
		}//for

		if(r>0)pw = pw << ((k-r)*bit);
		Tq = Tq ^ pw;

		//VERBOSE_CHANNEL<<endl<<"position = 0, hXOR = "<<Tq<<endl;

		HASHcounter[Tq]++;

		for (unsigned int j = s+1; j < e - (unsigned int)blockLength +1; j++) {
			c = TEXT[j + blockLength -1];
			randNum= (int)(rand()%(bit==1?2:4));
			switch (c) {
				case 'a' : case 'A' : character=0; break;
				case 'c' : case 'C' : character=1; break;
				case 'g' : case 'G' : character=(bit==1?0:2);break;
				case 't' : case 'T' : character=(bit==1?1:3);break;
				default : character=randNum; Nnumbers[numberOfNs]=randNum; numberOfNs++;break;
			}

			Tq = roll_right_XOR(Tq,lastChar[(j-1)%blockLength],character,blockLength);
			HASHcounter[Tq]++;
			lastChar[(j-1)%blockLength]=character;

			//if(j<100)VERBOSE_CHANNEL<<"position = "<< j << " hXOR = "<<Tq<<endl;

		}//for

		TEXT[e] = text_delimitator;
	}//forall contig (i)

	TEXT[size] = '\0';
	textLength=size+1;

	unsigned int t1 = HASHcounter[0];
	HASHcounter[0] = 0;
	for(unsigned int i=1; i< q ; i++) {
		unsigned int t2 = HASHcounter[i];
		HASHcounter[i] = HASHcounter[i-1] + t1;
		t1 = t2;
	}
	LAST = HASHcounter[q-1] + t1;

	HASHvalues = new unsigned int[textLength];

//// now build HASHvalues
	numberOfNs=0;
	for(int i =0; i< numContig; i++) {
		unsigned int s = start.at(i);
		unsigned int e = end.at(i);
		Tq =0;

		for(subblocks=0;subblocks<t;subblocks++){//for each sub-block(of size k) in the block
			pw = 0;
			for(j=0;j<k;j++){//for each char in a sub-block(of size k) in the block

				c = TEXT[s+subblocks*k + j];

				switch (c) {
					case 'a' : case 'A' : character=0; break;
					case 'c' : case 'C' : character=1; break;
					case 'g' : case 'G' : character=(bit==1?0:2);break;
					case 't' : case 'T' : character=(bit==1?1:3);break;
					default : character=Nnumbers[numberOfNs]; numberOfNs++; break;
				}

				lastChar[(s+subblocks*k + j)%blockLength]=character;
				pw = (pw << bit) ^ (unsigned long int)character;
			}//for
			Tq = Tq ^ pw;
		}//for all the sub-blocks(of size k) of the block

		pw = 0;
		for(j=0;j<r;j++){//for each char in the last block (only if there is a remainder r)

			c = TEXT[s+t*k + j];

			switch (c) {
				case 'a' : case 'A' : character=0; break;
				case 'c' : case 'C' : character=1; break;
				case 'g' : case 'G' : character=(bit==1?0:2);break;
				case 't' : case 'T' : character=(bit==1?1:3);break;
				default :  character=Nnumbers[numberOfNs]; numberOfNs++; break;
			}

			lastChar[(s+t*k + j)%blockLength]=character;
			pw = (pw << bit) ^ (unsigned long int)character;
		}//for

		if(r>0)pw = pw << ((k-r)*bit);
		Tq = Tq ^ pw;

		HASHvalues[HASHcounter[Tq]] = s;

		/*if(s==246594140+20) {

			VERBOSE_CHANNEL << "inserted, XOR = " << Tq << endl;

			for(int a=246594140+20; a<246594140+40+20;a++) VERBOSE_CHANNEL << TEXT[a];

			int* text = new int[blockLength];

			for(int kk=0;kk<blockLength;kk++){

				switch (TEXT[246594140+20+kk]) {
					case 'a' : case 'A' : text[kk] = 0; break;
					case 'c' : case 'C' : text[kk] = 1; break;
					case 'g' : case 'G' : text[kk] = 2;break;
					case 't' : case 'T' : text[kk] = 3;break;
					default : break;
				}


			}


			VERBOSE_CHANNEL << endl << "real hXOR = " << hXOR(text,(int)blockLength) << endl;

		}*/


		HASHcounter[Tq]++;
		//VERBOSE_CHANNEL<<endl<<"position = 0, hXOR = "<<Tq<<endl;

		for (unsigned int j = s+1; j < e - (unsigned int)blockLength +1; j++) {
			c = TEXT[j + blockLength -1];

			switch (c) {
				case 'a' : case 'A' : character=0; break;
				case 'c' : case 'C' : character=1; break;
				case 'g' : case 'G' : character=(bit==1?0:2);break;
				case 't' : case 'T' : character=(bit==1?1:3);break;
				default : character=Nnumbers[numberOfNs]; numberOfNs++; break;
			}

			Tq = roll_right_XOR(Tq,lastChar[(j-1)%blockLength],character,blockLength);
			HASHvalues[HASHcounter[Tq]] = j;




			/*if(j==246594140+20) {

				VERBOSE_CHANNEL << "inserted, XOR = " << Tq << endl;

				for(int a=246594140+20; a<246594140+40+20;a++) VERBOSE_CHANNEL << TEXT[a];

				int* text = new int[blockLength];

				for(int kk=0;kk<blockLength;kk++){

					switch (TEXT[246594140+20+kk]) {
						case 'a' : case 'A' : text[kk] = 0; break;
						case 'c' : case 'C' : text[kk] = 1; break;
						case 'g' : case 'G' : text[kk] = 2;break;
						case 't' : case 'T' : text[kk] = 3;break;
						default : break;
					}


				}


				VERBOSE_CHANNEL << endl << "real hXOR = " << hXOR(text,(int)blockLength) << endl;

			}*/




			HASHcounter[Tq]++;
			lastChar[(j-1)%blockLength]=character;
			//VERBOSE_CHANNEL<<"position = "<< j << " hXOR = "<<Tq<<endl;
		}//for
	}//forall contig (i)

	for(unsigned int i=q; i > 0 ; i--) {
		HASHcounter[i]= HASHcounter[i-1];
	}
	HASHcounter[0] = 0;

	VERBOSE_CHANNEL << "HASH table complete\n";


//now calcualte statistics fo rthis fata structure
/*		unsigned long int emptyLists = 0;
		unsigned long int nonEmptyLists = 0;
		float loadFactor = 0;
		unsigned long int meanLength = 0;
		for(int i = 0; i < q; i++) {
			if(HASH[i].length == 0) {
				emptyLists++;
			} else if (HASH[i].length > 0) {
				nonEmptyLists++;
				meanLength += HASH[i].length;
			} else {
				DEFAULT_CHANNEL << "what's up\n";
				return;
			}

		}

		float mean = (float)meanLength/(float)nonEmptyLists;
		float stdLength=0;
		for(int i = 0; i < q; i++) {
			if (HASH[i].length > 0) {
				float t = (HASH[i].length - mean)* (HASH[i].length - mean);
				stdLength += t;
			}
		}
		stdLength=stdLength/(float)meanLength;
		stdLength= sqrt(stdLength);



		DEFAULT_CHANNEL << "Elemnts = "<< meanLength <<"\n";
		loadFactor = (float)q/(float)meanLength;
		float loadFactorTest = (float)q/(float)textLength;
		float D = (float)emptyLists/(float)meanLength;

		DEFAULT_CHANNEL << "q = "<< q <<"\n";
		DEFAULT_CHANNEL << "emptyLists = "<< emptyLists <<"\n";
		DEFAULT_CHANNEL << "nonEmptyLists = "<< nonEmptyLists <<"\n";
		DEFAULT_CHANNEL << "meanLength = "<< mean <<"\n";
		DEFAULT_CHANNEL << "stdLength = "<< stdLength <<"\n";
		DEFAULT_CHANNEL << "loadFactor true = "<< loadFactor <<"\n";
		DEFAULT_CHANNEL << "loadFactor prior = "<< loadFactorTest <<"\n";
		DEFAULT_CHANNEL << "load factor wiki = " << D << "\n";
*/

}//createHASH


void Hash::save(const char * filename) throw (File_Not_Found){

	FILE *fp;
	if((fp=fopen(filename, "wb"))==NULL) {
	    printf("Cannot open file.\n");
	    exit(5);
	}

	bool *var0 = new bool[1];
	var0[0] = methyl_hash;
	int *var1 = new int[1];
	var1[0] = k;
	unsigned long int *var2 = new unsigned long int[1];
	var2[0] = q;
//	unsigned long int *var3 = new unsigned long int[1];
//	var3[0] = maskRollRight;
	unsigned int *var4 = new unsigned int[1];
	var4[0] = blockLength;
	unsigned int *var5 = new unsigned int[1];
	var5[0] = textLength;
	unsigned long int *var6 = new unsigned long int[1];
	var6[0] = h;
//	int *var7 = new int[1];
//	var7[0] = Zlength;
	unsigned int *var8 = new unsigned int[1];
	var8[0] = LAST;
	int *var9 = new int[1];
	var9[0] = globalToLocal.contigs;

	fwrite(var0 , sizeof(bool), 1 , fp); //methyl_hash
	fwrite(var1 , sizeof(int), 1 , fp); //k
	fwrite(var2 , sizeof(unsigned long int), 1 , fp); //q
//	fwrite(var3 , sizeof(unsigned long int), 1 , fp); // maskRollRight
	fwrite(var4 , sizeof(unsigned int), 1 , fp); // blockLength
	fwrite(var5 , sizeof(unsigned int), 1 , fp); //textLength
	fwrite(var6 , sizeof(unsigned long int), 1 , fp); // h
	//fwrite(var7 , sizeof(int), 1 , fp); // Zlength
	fwrite(var8 , sizeof(unsigned int), 1 , fp); // LAST
	fwrite(var9 , sizeof(int), 1 , fp); // Contigs

	//save the text
	fwrite(TEXT , sizeof(char), textLength , fp);

	//salvo le lunghezza degli header
	int *Lengths = new int[globalToLocal.contigs];
	for(int i=0; i < globalToLocal.contigs; i++) {
		Lengths[i] = globalToLocal.contigsNames[i].length();
	}
	fwrite(Lengths, sizeof(int), globalToLocal.contigs, fp);

	//salvo gli header
	for(int i=0; i < globalToLocal.contigs; i++) {
		fwrite(globalToLocal.contigsNames[i].c_str(), sizeof(char), globalToLocal.contigsNames[i].length()+1, fp);
	}
	fwrite(globalToLocal.startPositions, sizeof(unsigned int), globalToLocal.contigs, fp);
	fwrite(globalToLocal.endPositions, sizeof(unsigned int), globalToLocal.contigs, fp);

	fwrite(HASHvalues, sizeof(unsigned int), textLength, fp);
	fwrite(HASHcounter, sizeof(unsigned int), q, fp);

	fclose(fp);

}//save

#define check_numBytes() if (numBytes == 0) { VERBOSE_CHANNEL << "Read 0 bytes when reading file " << filename << endl; exit(7); }

/*
void Hash::loadText(const char * filename) throw (File_Not_Found){
	unsigned long int numBytes;
	FILE *fp;
	if((fp=fopen(filename, "rb"))==NULL) {
		VERBOSE_CHANNEL << "Cannot open file " << filename << endl;
	   exit(2);
	}

	//DEFAULT_CHANNEL << "Loading the Hash table\n";

	int *var1 = new int[1];
		numBytes = fread(var1 , sizeof(int), 1 , fp); //k
		check_numBytes();
		k = var1[0];

	//	DEFAULT_CHANNEL << "k " << k << "\n";

		unsigned long int *var2 = new unsigned long int[1];
		numBytes = fread(var2 , sizeof(unsigned long int), 1 , fp); // q
		check_numBytes();
		q = var2[0];
	//	DEFAULT_CHANNEL << "q "<< q << "\n";
	//	unsigned long int *var3 = new unsigned long int[1];
	//	numBytes = fread(var3 , sizeof(unsigned long int), 1 , fp);
	//	maskRollRight = var3[0];
	//	DEFAULT_CHANNEL << "maskrollRight " << maskRollRight << "\n";
		unsigned int *var4 = new unsigned int[1];
		numBytes = fread(var4 , sizeof(unsigned int), 1 , fp);
		check_numBytes();
		blockLength = var4[0];
	//	DEFAULT_CHANNEL << "blockLength " << blockLength << "\n";
		unsigned int *var5 = new unsigned int[1];
		numBytes = fread(var5 , sizeof(unsigned int), 1 , fp);
		check_numBytes();
		textLength = var5[0];
	//	DEFAULT_CHANNEL << "textLength " << textLength << "\n";
		unsigned long int *var6 = new unsigned long int[1];
		numBytes = fread(var6 , sizeof(unsigned long int), 1 , fp);
		check_numBytes();
		h = var6[0];
	//	DEFAULT_CHANNEL << "h " << h << "\n";
	//	int *var7 = new int[1];
	//	numBytes = fread(var7 , sizeof(int), 1 , fp);
	//	check_numBytes();
	//	Zlength = var7[0];
	//	DEFAULT_CHANNEL <<"ZLength "<< Zlength << "\n";
		unsigned int *var8 = new unsigned int[1];
		numBytes = fread(var8 , sizeof(unsigned int), 1 , fp);
		check_numBytes();
		LAST = var8[0];
	//	DEFAULT_CHANNEL << "LAST " << LAST << "\n";
		int *var9 = new int[1];
		numBytes = fread(var9 , sizeof(int), 1 , fp);
		check_numBytes();
		int Contigs =var9[0];
	//	DEFAULT_CHANNEL << "contigs  " << Contigs << "\n";

		TEXT = new char[textLength+1];
		numBytes = fread(TEXT , sizeof(char), textLength , fp);
		check_numBytes();

		globaltolocal = GlobalToLocal(Contigs);
		int *Lengths = new int[globaltolocal.Contigs];
		numBytes = fread(Lengths, sizeof(int), globaltolocal.Contigs, fp);
		check_numBytes();

		for(int i=0; i < globaltolocal.Contigs; i++) {
			char *contigs = new char[Lengths[i]+1];
			contigs[0]='\0';
			numBytes = fread(contigs, sizeof(char), Lengths[i]+1, fp);
			check_numBytes();
			globaltolocal.contigsNames[i] = contigs;
			delete [] contigs;
		}

		numBytes = fread(globaltolocal.startPositions, sizeof(unsigned int), globaltolocal.Contigs, fp);
		check_numBytes();
		numBytes = fread(globaltolocal.endPositions, sizeof(unsigned int), globaltolocal.Contigs, fp);
		check_numBytes();


		fclose(fp);

		delete [] Lengths;
		delete [] var1;
		delete [] var2;
		delete [] var4;
		delete [] var5;
		delete [] var6;
	//	delete [] var7;
		delete [] var8;
		delete [] var9;

}//loadText
*/
void Hash::load(const char * filename) throw (File_Not_Found){
	unsigned long int numBytes;
	FILE *fp;
	if((fp=fopen(filename, "rb"))==NULL) {
		VERBOSE_CHANNEL << "Cannot open file " << filename << endl;
	   exit(2);
	}

	//DEFAULT_CHANNEL << "Loading the Hash table\n";

	bool *var0 = new bool[1];
		numBytes = fread(var0 , sizeof(bool), 1 , fp); //methyl_hash
		check_numBytes();
		methyl_hash = var0[0];

	int *var1 = new int[1];
		numBytes = fread(var1 , sizeof(int), 1 , fp); //k
		check_numBytes();
		k = var1[0];

	//	DEFAULT_CHANNEL << "k " << k << "\n";

		unsigned long int *var2 = new unsigned long int[1];
		numBytes = fread(var2 , sizeof(unsigned long int), 1 , fp); // q
		check_numBytes();
		q = var2[0];
	//	DEFAULT_CHANNEL << "q "<< q << "\n";
	//	unsigned long int *var3 = new unsigned long int[1];
	//	numBytes = fread(var3 , sizeof(unsigned long int), 1 , fp);
	//	maskRollRight = var3[0];
	//	DEFAULT_CHANNEL << "maskrollRight " << maskRollRight << "\n";
		unsigned int *var4 = new unsigned int[1];
		numBytes = fread(var4 , sizeof(unsigned int), 1 , fp);
		check_numBytes();
		blockLength = var4[0];
	//	DEFAULT_CHANNEL << "blockLength " << blockLength << "\n";
		unsigned int *var5 = new unsigned int[1];
		numBytes = fread(var5 , sizeof(unsigned int), 1 , fp);
		check_numBytes();
		textLength = var5[0];
	//	DEFAULT_CHANNEL << "textLength " << textLength << "\n";
		unsigned long int *var6 = new unsigned long int[1];
		numBytes = fread(var6 , sizeof(unsigned long int), 1 , fp);
		check_numBytes();
		h = var6[0];
	//	DEFAULT_CHANNEL << "h " << h << "\n";
	//	int *var7 = new int[1];
	//	numBytes = fread(var7 , sizeof(int), 1 , fp);
	//	check_numBytes();
	//	Zlength = var7[0];
	//	DEFAULT_CHANNEL <<"ZLength "<< Zlength << "\n";
		unsigned int *var8 = new unsigned int[1];
		numBytes = fread(var8 , sizeof(unsigned int), 1 , fp);
		check_numBytes();
		LAST = var8[0];
	//	DEFAULT_CHANNEL << "LAST " << LAST << "\n";
		int *var9 = new int[1];
		numBytes = fread(var9 , sizeof(int), 1 , fp);
		check_numBytes();
		int Contigs =var9[0];
	//	DEFAULT_CHANNEL << "contigs  " << Contigs << "\n";


		VERBOSE_CHANNEL << "text length = " << textLength << endl;

		TEXT = new char[textLength+1];
		numBytes = fread(TEXT , sizeof(char), textLength , fp);
		check_numBytes();

		globalToLocal = GlobalToLocal(Contigs);
		int *lengths = new int[globalToLocal.contigs];
		numBytes = fread(lengths, sizeof(int), globalToLocal.contigs, fp);
		check_numBytes();

		for(int i=0; i < globalToLocal.contigs; i++) {
			char *contigs = new char[lengths[i]+1];
			//contigs[0]='\0';
			numBytes = fread(contigs, sizeof(char), lengths[i]+1, fp);
			check_numBytes();
			globalToLocal.contigsNames[i] = contigs;
			delete [] contigs;
		}

		numBytes = fread(globalToLocal.startPositions, sizeof(unsigned int), globalToLocal.contigs, fp);
		check_numBytes();
		numBytes = fread(globalToLocal.endPositions, sizeof(unsigned int), globalToLocal.contigs, fp);
		check_numBytes();

	/*
	for(int i=0; i<  globaltolocal.Contigs; i++) {
		DEFAULT_CHANNEL << globaltolocal.contigsNames[i]<< " " << globaltolocal.startPositions[i] << " " << globaltolocal.endPositions[i] << "\n";
	}
	*/

		HASHcounter = new unsigned int[q];
		HASHvalues = new unsigned int[textLength];

		numBytes = fread(HASHvalues, sizeof(unsigned  int), textLength, fp);
		check_numBytes();
		numBytes = fread(HASHcounter, sizeof(unsigned int), q, fp);
		check_numBytes();

		fclose(fp);

		MASK = 0;
		int h=0;

		for(int i = sizeof(unsigned int)*8 -1; i >= 2; i--) {
			h++;
			if(h==3) {
				MASK |= (1 << i);
				h=0;
			}
		}

		//create set Z:

		VERBOSE_CHANNEL << "Building Z Set" << endl;

		Z = new SetZ();

		Z->initZ(2,k,methyl_hash);

		VERBOSE_CHANNEL << "Z has " << Z->ZSize << " elements" << endl;

		//--------------------entropy on 30 bits
		/*double ni = 0;
		double H = 0;
		for(unsigned long int i=0;i<q-1;i++){

			ni = HASHcounter[i+1] - HASHcounter[i];

			if(ni>0)
				H += ni * log(textLength/ni)/log(2);

		}

		ni = LAST - HASHcounter[q-1];
		if(ni>0)
			H += ni * log(textLength/ni)/log(2);

		H /= textLength;

		VERBOSE_CHANNEL << "Hash entropy (/30 bits) = " << H << " bits" << endl;

		//----------------------entropy on 15 bits

		ni = 0;
		H = 0;
		unsigned long int sz=32768;
		unsigned long int *val = new unsigned long int[sz];
		unsigned long int maskR = sz - 1;

		for(unsigned int i=0;i<sz;i++)
			val[i]=0;

		for(unsigned long int i=0;i<q-1;i++){

			ni = HASHcounter[i+1] - HASHcounter[i];

			val[i & maskR]+=ni;
			val[i>>15]+=ni;

		}

		ni = LAST - HASHcounter[q-1];
		val[(q-1) & maskR]+=ni;
		val[(q-1)>>15]+=ni;

		for(unsigned long int i=0;i<sz;i++){

			if(val[i]>0)
				H += val[i] * log(2*textLength/val[i])/log(2);

		}

		H /= (textLength*2);

		VERBOSE_CHANNEL << "Hash entropy (/15 bits) = " << H << " bits" << endl;
*/

//--------------------------------------------------------


		delete [] lengths;
		delete [] var0;
		delete [] var1;
		delete [] var2;
		delete [] var4;
		delete [] var5;
		delete [] var6;
	//	delete [] var7;
		delete [] var8;
		delete [] var9;

}//load

int Hash::stranded_gap_search_subroutine(ItemsGapped & output, Items & left, Items & right, t_length max_gap, const char * s, int left_characters, int right_characters, int max_errors) {
	Items::iterator left_end  = left.end();
	Items::iterator right_end = right.end();
	Items::iterator left_min  = left.begin();
	Items::iterator right_min = right.begin();
	Items::iterator right_max = right.begin();
	const int length  = strlen(s);

	while ((right_min != right_end) and ((left_min->globalPosition+left_characters) > right_min->globalPosition)) {
		if (right_min == right_max)
			right_max++;
		right_min++;
	}

	while ((right_max != right_end) and ((right_max+1) != right_end) and ((left_min->globalPosition+left_characters+max_gap) >= (right_max+1)->globalPosition))
		right_max++;

	while ((left_min != left_end) and (right_min != right_end)) {
		if (((left_min->globalPosition+left_characters) <= right_min->globalPosition) and ((left_min->globalPosition+left_characters + max_gap) >= right_min->globalPosition)) {
			bool more = true;
			for (Items::iterator iter = right_min; more and (iter != right_end); iter++) {
				if (iter == right_max)
					more = false;

				// Check if they are from the same sequence
				int contig = globalToLocal.searchContig(left_min->globalPosition);
				if (contig == globalToLocal.searchContig(iter->globalPosition)) {

					int errors1 = left_min->errors;
					int errors2 = iter->errors;

					int length1 = left_characters;
					int length2 = right_characters;

					int local_position_left = length1-1;
					int local_position_right = length - length2;
					int global_position_left = left_min->globalPosition + length1 - 1;
					int global_position_right = iter->globalPosition;

					bool try_error_left = true;
					while ( (errors1 + errors2 <= max_errors) and (length1 + length2 < length) and (global_position_left < global_position_right)) {
						// try without errors
						bool errors_left = false;
						bool errors_right = false;
						bool try_left = true;
						while ( (not errors_left or not errors_right) and (errors1 + errors2 <= max_errors) and (length1 + length2 < length) ) {
							if (try_left and not errors_left) {
								if (TEXT[global_position_left+1] == s[local_position_left+1]) {
									global_position_left++;
									local_position_left++;
									length1++;
								} else
									errors_left = true;
								try_left = false;
							} else if (not try_left and not errors_right) {
								if (TEXT[global_position_right-1] == s[local_position_right-1]) {
									global_position_right--;
									local_position_right--;
									length2++;
								} else
									errors_right = true;
								try_left = true;
							} else {
								errors_left = false;
								errors_right = false;
							}
						}

						if ( (errors1 + errors2 <= max_errors) and (length1 + length2 < length) ) {
							if (try_error_left) {
								errors1++;
								length1++;
								local_position_left++;
								global_position_left++;
								try_error_left = false;
							} else {
								errors2++;
								length2++;
								local_position_right--;
								global_position_right--;
								try_error_left = true;
							}
						}
					}

					if ( (errors1 + errors2 <= max_errors) and (length1 + length2 == length) and (global_position_left < global_position_right)) {
						if (errors1 + errors2 < max_errors) {
							output.clear();
							max_errors = errors1 + errors2;
						}

						if (left_min->strand)
							output.push_back(ResultItemGapped(left_min->globalPosition,global_position_right,left_min->strand,length1,length2,errors1,errors2,contig));
						else
							output.push_back(ResultItemGapped(left_min->globalPosition-1,global_position_right-1,left_min->strand,length1,length2,errors1,errors2,contig));
					}
				}
			}
			left_min++;
		} else {
			while ((left_min != left_end) and ((left_min->globalPosition+left_characters + max_gap) <= right_min->globalPosition))
				left_min++;

			while ((right_min != right_end) and ((left_min->globalPosition+left_characters) > right_min->globalPosition)) {
				if (right_min == right_max)
					right_max++;
				right_min++;
			}

			while ((right_max != right_end) and ((right_max+1) != right_end) and ((left_min->globalPosition+left_characters+max_gap) >= (right_max+1)->globalPosition))
				right_max++;
		}
	}

	return max_errors;
}

unsigned long int Hash::hXOR(int* P, int bl){

	int bit;

	if(methyl_hash)
		bit=1;
	else
		bit=2;
	int t = bl/k;//number of blocks in P having size k
	int r = bl % k;//size of the last block in P

	unsigned long int hXOR = 0;

	unsigned long int pw;
	int i,j;

	for(i=0;i<t;i++){//for each block in P

		pw = 0;
		for(j=0;j<k;j++)//for each char in a block(of size k) in P
			pw = (pw << bit) ^ (unsigned long int)*(P + i*k + j);

		hXOR = hXOR ^ pw;

	}//for

	pw = 0;

	for(j=0;j<r;j++)//for each char in the last block (only if there is a remainder r)
		pw = (pw << bit) ^ (unsigned long int)*(P + t*k + j);

	pw = pw << ((k-r)*bit);
	hXOR = hXOR ^ pw;

	//VERBOSE_CHANNEL << "hXOR = " << hXOR << endl;

	return hXOR;

}//hXOR

void Hash::search(const string & read, Items & output, unsigned int max_errors, unsigned int & delta, vector<pair<int, int> > & DeltaSolutions, bool ind) {

	unsigned int patternLength = read.length();

	if(patternLength <= blockLength) {
		return;
	}
	unsigned int numberOfBlocks = patternLength/blockLength;
	unsigned int blockErrorsIndel = 2;
	unsigned int blockErrors = (max_errors/numberOfBlocks > 2 ) ? 2 : max_errors/numberOfBlocks;
	//(max_errors/numberOfBlocks > 1 ) ? blockErrorsIndel = 1:blockErrorsIndel=max_errors/numberOfBlocks;

	int pattern[2][patternLength];
	char PATTERN[2][patternLength];

	unsigned int p;
	unsigned int l;
	int bit;

	if(methyl_hash)
		bit=1;
	else
		bit=2;

	int randNum;
	unsigned int maxErrorInRead = max_errors;
	// processing the read
	for (unsigned int plTemp = 0; plTemp < patternLength ; plTemp++) {
		randNum= 1; //(int)(rand()%(bit==1?2:4));
		switch (read.at(plTemp)) {
		case 'A' : case 'a' : pattern[0][plTemp] = 0;pattern[1][patternLength- plTemp-1]= (bit==1?1:3); PATTERN[0][plTemp]='A'; PATTERN[1][patternLength- plTemp-1]='T'; break;
		case 'C' : case 'c' : pattern[0][plTemp] = 1;pattern[1][patternLength- plTemp-1] = (bit==1?0:2); PATTERN[0][plTemp]='C'; PATTERN[1][patternLength- plTemp-1]='G'; break;
		case 'G' : case 'g' : pattern[0][plTemp] = (bit==1?0:2);pattern[1][patternLength- plTemp-1] = 1; PATTERN[0][plTemp]='G'; PATTERN[1][patternLength- plTemp-1]='C'; break;
		case 'T' : case 't' : pattern[0][plTemp] = (bit==1?1:3);pattern[1][patternLength- plTemp-1] = 0; PATTERN[0][plTemp]='T'; PATTERN[1][patternLength- plTemp-1]='A'; break;
		default  : pattern[0][plTemp] = randNum;pattern[1][patternLength- plTemp-1] = randNum ; PATTERN[0][plTemp]='N'; PATTERN[1][patternLength- plTemp-1]='N'; break;
		}

	}

	bool foundExact = false;  // IF AN OCCURENCE IS FOUND WITH 0 ERRORS DON'T SEARCH FOR OTHER OCCURENCES WITH MORE THEN ONE ERROR
	unsigned long int Pq=0;

	unsigned int errors = blockErrors;
	unsigned int diffTemp;
	int idxZ = 0;
	unsigned long int idxText = 0;
	unsigned int j;
	unsigned int numN;
	unsigned int dTemp;

	int partition;
	bool get_z_element;

	if(delta == 0) {
		for(int pp=0; pp<=1; pp++) { // reverse complement
			foundExact = false;

			for (unsigned int ii = 0; ii < numberOfBlocks && !foundExact; ii++) {

				numN = 0;

				Pq = hXOR(pattern[pp] + ii * blockLength,blockLength);//hXOR

				for (j = ii * blockLength; j < (ii + 1) * blockLength; j++)
					if(PATTERN[pp][j]=='N') numN++;

				if(numN <= errors ) {

					get_z_element = true;
					idxZ = 0;
					partition = 0;

					while(get_z_element && !foundExact){//list Z elements

						idxText = Z->getElement(partition,idxZ)^Pq;

						p = HASHcounter[idxText];

						(idxText == q-1) ? l= LAST:l = HASHcounter[idxText+1];

						while(p < l) {//list HASHValues elements

							if (HASHvalues[p] >= ii*blockLength) {

								diffTemp = HASHvalues[p] - ii*blockLength;
								if (diffTemp + patternLength < textLength ) {// && Checked.insert(diffTemp).second) { // we check that ideed *iter - i * pattern_length is a valid shift in the text
									dTemp = 0;

									for (int j = 0;  dTemp <= maxErrorInRead && (j < (int)patternLength); j++) {

										if(pp==0){//forward

											if ((PATTERN[pp][j]=='T' || PATTERN[pp][j]=='t') &&
												(TEXT[diffTemp+j]=='C' || TEXT[diffTemp+j]=='c')) dTemp += (methyl_hash?0:1);
											else if (PATTERN[pp][j] != TEXT[diffTemp+j]) dTemp += 1;

										}else{//reverse

											if ((PATTERN[pp][j]=='A' || PATTERN[pp][j]=='a') &&
												(TEXT[diffTemp+j]=='G' || TEXT[diffTemp+j]=='g')) dTemp += (methyl_hash?0:1);
											else if(PATTERN[pp][j] != TEXT[diffTemp+j]) dTemp += 1;

										}//reverse

										if(TEXT[diffTemp+j]== text_delimitator)	dTemp = max_errors + 1;

									}//compute distance

									if (dTemp <= maxErrorInRead) { // a good hit
										ResultItem K;

										K.globalPosition = diffTemp;
										K.errors = dTemp;//(int)
										K.strand = (pp == 0)?true:false;//forward or reverse complement
										K.indels = false;

										if(dTemp == maxErrorInRead) {
											output.push_back(K);
										} else {
											errors=dTemp/numberOfBlocks;
											if(errors > 2) errors = 2;
											maxErrorInRead = dTemp;
											output.clear();
											output.push_back(K);
											if (dTemp == 0) foundExact = true;
										}
									}
								}
							}
							p++;
						}

						Z->updateIndex(&partition,&idxZ,errors,&get_z_element);

					}
				}
			}
		}
	} else {//delta != 0

		unsigned int totalErrors = max_errors + delta;
		(totalErrors/numberOfBlocks > 2 ) ? blockErrors = 2:blockErrors=totalErrors/numberOfBlocks;
		unsigned int  maxErrorInRead = totalErrors; //unsigned int maxErrorInRead = max_errors;
		unsigned int  currentError = 0;
		unsigned int  bestError = maxErrorInRead;

		Items * foundSolutions = new Items[(int)maxErrorInRead+1];

		for(unsigned int pp=0; pp<=1; pp++) { // reverse complement
			for (unsigned int ii = 0; ii < numberOfBlocks; ii++) {
				numN = 0;

				Pq = hXOR(pattern[pp] + ii * blockLength,blockLength);//hXOR

				for (j = ii * blockLength; j < (ii + 1) * (int)blockLength; j++)
					if(PATTERN[pp][j]=='N') numN++;

				if(numN <= errors ) {

					get_z_element = true;
					idxZ = 0;
					partition = 0;

					while(get_z_element && !foundExact){

						idxText = ( (Z->getElement(partition,idxZ))^Pq );

						p = HASHcounter[idxText];
						(idxText == q-1) ? l= LAST:l = HASHcounter[idxText+1];

						while(p < l) {
							if (HASHvalues[p] >= ii*blockLength) {
								diffTemp = HASHvalues[p] - ii*blockLength;
								if (diffTemp + patternLength < textLength ) {// && Checked.insert(diffTemp).second) { // we check that ideed *iter - i * pattern_length is a valid shift in the text
									currentError = 0;
									for (unsigned int j = 0;  currentError <= maxErrorInRead && (j < patternLength); j++) {

										if(pp==0){//forward

											if ((PATTERN[pp][j]=='T' || PATTERN[pp][j]=='t') &&
													(TEXT[diffTemp+j]=='C' || TEXT[diffTemp+j]=='c')) currentError += (methyl_hash?0:1);
											else if(PATTERN[pp][j] != TEXT[diffTemp+j]) currentError += 1;

										}else{//reverse

											if ((PATTERN[pp][j]=='A' || PATTERN[pp][j]=='a') &&
											(TEXT[diffTemp+j]=='G' || TEXT[diffTemp+j]=='g')) currentError += (methyl_hash?0:1);
											else if(PATTERN[pp][j] != TEXT[diffTemp+j]) currentError += 1;

										}//reverse

										if(TEXT[diffTemp+j]== text_delimitator)	currentError = maxErrorInRead + 1;

									}

									if (currentError <= maxErrorInRead) { // a good hit
										ResultItem K;
										K.globalPosition = diffTemp;
										K.errors = currentError;
										K.strand = (pp == 0)?true:false;
										K.indels = false;

										foundSolutions[(int)currentError].push_back(K);

										if (currentError < bestError) {
											(totalErrors<(currentError+delta))?maxErrorInRead=totalErrors:maxErrorInRead=currentError+delta;
											bestError = currentError;
											errors = maxErrorInRead/numberOfBlocks;
											if(errors > 2) errors = 2;
										}
									}
								}
							}
							p++;
						}

						Z->updateIndex(&partition,&idxZ,errors,&get_z_element);

					}
				}
			}
		}

		if(bestError <= max_errors) {
			for(unsigned int i = bestError; i<= bestError + delta; i++ ) {
				sort(foundSolutions[i].begin(), foundSolutions[i].end(), ResultItem::less()); // sort solutions
				foundSolutions[i].erase(unique(foundSolutions[i].begin(), foundSolutions[i].end(), ResultItem::equal()), foundSolutions[i].end());
				output.insert(output.end(), foundSolutions[i].begin(), foundSolutions[i].end());
				DeltaSolutions.push_back(pair<int, int>(i,foundSolutions[i].size()));
			}
		}
		delete [] foundSolutions;
	}

	if(ind && output.size() == 0) {//gap
		int score;

		int bestScore = (patternLength -indels) - max_errors*match - indels*indel -1;
		int EditOperations = 0;
		errors = blockErrorsIndel;
		bool WithIndel;
		for(int pp=0; pp<=1; pp++) { // reverse complement
			foundExact = false;

			Pq = hXOR(pattern[pp], blockLength);//hXOR

			get_z_element = true;
			idxZ = 0;
			partition = 0;

			while(get_z_element && !foundExact){

				idxText = ( (Z->getElement(partition,idxZ))^Pq );

				p = HASHcounter[idxText];
				(idxText == q-1) ? l= LAST:l = HASHcounter[idxText+1];

				while(p < l) {
					diffTemp = HASHvalues[p];
					if (diffTemp + patternLength < textLength ) {// && Checked.insert(diffTemp).second) { // we check that ideed *iter - i * pattern_length is a valid shift in the text
						score = alignSW_Improved(PATTERN[pp], patternLength, diffTemp , max_errors, EditOperations, WithIndel, bestScore);
						if(score >= bestScore) {
							ResultItem K;
							K.indels= WithIndel;
							K.strand = (pp == 0)?true:false;
							K.errors = EditOperations;
							K.globalPosition = diffTemp;
							if(score == bestScore) {
								output.push_back(K);
							} else if (score > bestScore) {
								output.clear();
								output.push_back(K);
								bestScore = score;
								if(EditOperations == 0) { // means a perfect match has been found
									foundExact = true;
									errors = 0;
								}
							}
						}
					}
					p++;
				}

				Z->updateIndex(&partition,&idxZ,errors,&get_z_element);

			}
		}
	}
}//search

int  Hash::max( int f1, int f2, int f3, int & ptr ) {
	int max;

	if( f1 >= f2 && f1 >= f3 ){
		max = f1 ;
		ptr = 0 ;
	} else if( f2 > f3 ){
		max = f2 ;
		ptr = 1 ;
	}else{
		max = f3 ;
		ptr = 2  ;
	}
	return  max ;
}

int Hash::alignSW_Improved( const char *read, int patternLength, unsigned int globalPos, int errors, int & EditOperations, bool & WithIndel, int BestAchivableScore) {
	int totalErrors = errors + indels;
	int SW[2][patternLength + indels +1][3];
	for(int i = 0; i<= patternLength + indels; i++) { // TO DO: this initializationc an be reduced!!!
		SW[0][i][0]=-100; SW[0][i][1]=SW[0][i][2]=0;
		SW[1][i][0]=-100; SW[1][i][1]=SW[1][i][2]=0;
	}
	SW[0][0][0]=0;

	int start, stop, a1, a2, x, y;
	a1=a2=0;
	char nuc;
	int fD,fU, fL, ptr;
	int maxScoreIteration=0;
	int minMisInd;

	for(int i=1; i<=patternLength; i++){
		minMisInd = totalErrors +1 ;
		maxScoreIteration=-i*match;

		a1 = (i-1) % 2; // compute the array that must be used
		a2 = i %  2 ;
		(i-indels>0) ? start = i-indels : start = 1;
		//(i+indels<=patternLength)? stop = i+indels: stop = patternLength; //check the Matrix boundaries
		stop = i+indels;
		SW[a2][start-1][0]= -100; //initialize the first
		for(int j=start; j<=stop; j++) {
			nuc = TEXT[globalPos + j -1];
			switch( nuc ) {
			case 'A':  x = 0 ;  break ;
			case 'C':  x = 1 ;  break ;
			case 'G':  x = 2 ;  break ;
			case 'T':  x = 3 ; break;
			case 'N':  x = 4 ; break;
			case '$': return -1; break;
			default:
				//cout << nuc << " error in ref\n"; return -1;
				x = 4 ; break;
			}
			nuc = read[i - 1];
			switch( nuc ) {
			case 'A':  y = 0 ;  break ;
			case 'C':  y = 1 ;  break ;
			case 'G':  y = 2 ;  break ;
			case 'T':  y = 3 ;  break;
			case 'N':  y = 4; 	break;
			default:
				//cout << nuc << " error in read\n"; return -1;
				y = 4; 	break;
			}

			fD = SW[a1][j-1][0] + SM[x][y] ;
			fU = SW[a2][j-1][0] - indel;
			fL = SW[a1][j][0] - indel;

			SW[a2][j][0] = max(fD, fU, fL, ptr) ;
			if((SW[a2][j][0] + (patternLength - i)*match) >= maxScoreIteration ) {
				maxScoreIteration = (SW[a2][j][0] + (patternLength - i)*match); // maximum achievable score
			}

			switch(ptr) {
				case 0:				// score in (i,j) stems from a match/mismatch
					(x!=y)? SW[a2][j][1] = SW[a1][j-1][1] + 1 : SW[a2][j][1] = SW[a1][j-1][1]; // update number of mismatches
					SW[a2][j][2] =SW[a1][j-1][2]; // copy number of indels
					break;
				case 1:             // score in (i,j) stems from a deletion in sequence A
					SW[a2][j][2] = SW[a2][j-1][2] + 1; // update number of indels
					SW[a2][j][1] = SW[a2][j-1][1]; // copy number of mismatches
					break;
				case 2:             // score in (i,j) stems from a deletion in sequence B
					SW[a2][j][1] = SW[a1][j][1]; // copy number of mismatches
					SW[a2][j][2] = SW[a1][j][2] + 1; // update number of indels
				break;
			}
			if(SW[a2][j][1] + SW[a2][j][2] < minMisInd )
				minMisInd = SW[a2][j][1] + SW[a2][j][2];

	//		cout << "(" << SW[a2][j][0] << "," << SW[a2][j][1] << "," << SW[a2][j][2] << ") ";


		}
	//	cout << "\n";
		if(maxScoreIteration < BestAchivableScore or minMisInd > totalErrors) {
			return -1; // stop computsation no way to reach a result as good as the one already achived
		}

	}

	int SW_max = 0;
	int SW_errors = 0;
	int SW_indels = 0;

	for(int i=patternLength-indels;i<=patternLength + indels ;i++){
		if(SW[a2][i][0] >SW_max){
			SW_max = SW[a2][i][0];
			SW_errors = SW[a2][i][1];
			SW_indels = SW[a2][i][2];
		}
	}



	if(SW_errors <= errors and SW_indels <= indels) {
		EditOperations = SW_errors + SW_indels;
		(SW_indels > 0)? WithIndel=true : WithIndel=false;
		return SW_max;
	}
	return -1;
}

typedef pair<int, int> cell; // score and mismatches

vector<pair<char, unsigned int> >  Hash::alignSW( const char *read, int patternLength, unsigned int globalPos) {

	vector<pair<char, unsigned int> > output;

	cell SW[patternLength + indels +1][patternLength+1]; // create SW matrix
	int I_i[patternLength + indels +1][patternLength+1],I_j[patternLength + indels +1][patternLength+1];  // Index matrices to remember the 'path' for backtracking

	// initialise the diagonal
	int start, stop;
	for(int i=0;i<=patternLength + indels ;i++){ //
		(i-indels > 0) ? start = i-indels: start = 0;
		(i+indels<patternLength)? stop = i+indels: stop = patternLength; //check the Matrix boundaries
		for(int j=start; j<=stop; j++) {
			SW[i][j].first = 0;
			SW[i][j].second = 0;
		}
		if(start > 0) {
			SW[i][start-1].first = -100;
		}
		if(stop < patternLength) {
			SW[i][stop+1].first = -100;
		}
	}
	for(int i=0;i<= 1+ indels ;i++){ //
		SW[i][0].first = -100;
		SW[0][i].first = -100;
	}
	SW[0][0].first = 0;

	int fU, fD, fL ;
	char nuc;
	int ptr;




	int x, y;
	for(int i=1;i<=patternLength + indels ;i++){ //
		(i-indels > 0) ? start = i-indels: start = 1;
		(i+indels<patternLength)? stop = i+indels: stop = patternLength; //check the Matrix boundaries
		for( int j=start; j<=stop; j++) {
			nuc = TEXT[globalPos + i -1];
			switch( nuc ) {
			case 'A':  x = 0 ; break ;
			case 'C':  x = 1 ; break ;
			case 'G':  x = 2 ; break ;
			case 'T':  x = 3 ; break;
			default:   x = 4; break;
			}
			nuc = read[j - 1];
			switch( nuc ) {
			case 'A':  y = 0 ; break ;
			case 'C':  y = 1 ; break ;
			case 'G':  y = 2 ; break ;
			case 'T':  y = 3 ; break;
			default:   y = 4; break;
			}

			fD = SW[i-1][j-1].first + SM[x][y] ;
			fU = SW[i-1][j].first - indel;
			fL = SW[i][j-1].first - indel;

			SW[i][j].first = max(fD, fU, fL, ptr) ;

			switch(ptr) {
			case 0: 								// score in (i,j) stems from a match/mismatch
				I_i[i][j] = i-1;
				I_j[i][j] = j-1;
//				(x!=y)? SW[i][j].second = SW[i-1][j-1].second + 1 : SW[i][j].second = SW[i-1][j-1].second;
				break;
			case 1:                                 // score in (i,j) stems from a deletion in sequence A
				I_i[i][j] = i-1;
				I_j[i][j] = j;
//				SW[i][j].second = SW[i-1][j].second + 1;
				break;
			case 2:                                 // score in (i,j) stems from a deletion in sequence B
				I_i[i][j] = i;
				I_j[i][j] = j-1;
//				SW[i][j].second = SW[i][j-1].second + 1;
				break;
			}
		}
	}

/*
	for(int i=0;i<=patternLength + indels ;i++){
		for(int j=0;j<=patternLength  ;j++){
			cout << SW[i][j].first << " ";
		}
		cout << "\n";
	}
*/

//compute the max that must be in the slice of the last column
	int SW_max = 0;
	int i_max=0,j_max=patternLength;
	for(int i=patternLength-indels;i<=patternLength + indels ;i++){
		if(SW[i][j_max].first>SW_max){
			SW_max = SW[i][j_max].first;
			i_max = i;
		}
	}

	// Backtracking from H_max
	int current_i=i_max,current_j=j_max;
	int next_i=I_i[current_i][current_j];
	int next_j=I_j[current_i][current_j];
	int tick=0;
//	char consensus_a[patternLength+patternLength + indels + 2],consensus_b[patternLength+patternLength + indels + 2];

	unsigned int operations = 0;
	char currentOp='a';
	char Op = 'n';

	//while(((current_i!=next_i) || (current_j!=next_j)) && (next_j!=0) && (next_i!=0)){
	while(current_i!=0 or current_j != 0 ){
		if(next_i == current_i) {
			currentOp = 'I';
//						consensus_a[tick] = '-';                  // deletion in A
		} else {
			currentOp = 'M';
//						consensus_a[tick] = TEXT[globalPos + current_i-1];   // match/mismatch in A
		}
		if(next_j==current_j) {
			currentOp = 'D';
//						consensus_b[tick] = '-';                  // deletion in B
		} else {
//						consensus_b[tick] = read[current_j-1];   // match/mismatch in B
		}
		if(Op == 'n') {
			Op = currentOp;
			operations=1;
		} else if (Op != currentOp) { // means that I performed a new operation
			pair<char,unsigned int> t(Op,operations);
			output.insert(output.begin(), t);
			Op = currentOp;
			operations = 1;
		} else {
			operations++;
		}

		current_i = next_i;
		current_j = next_j;
		next_i = I_i[current_i][current_j];
		next_j = I_j[current_i][current_j];
		tick++;
	}

	pair<char,unsigned int> t(Op,operations);
	output.insert(output.begin(), t);

/*
	// Output of the consensus motif to the console
		for(unsigned int i=0;i<patternLength + indels ;i++){cout<<TEXT[globalPos + i];}; cout<<"  and"<<endl;
		for(unsigned int i=0;i<patternLength;i++){cout<< read[i];}; cout<<endl<<endl;
		for(int i=tick-1;i>=0;i--) cout<<consensus_a[i];
		cout<<endl;
		for(int j=tick-1;j>=0;j--) cout<<consensus_b[j];
		cout<<endl;
*/
		return output;

}


void Hash::search(const string & read, Items & output, unsigned  int max_errors) {
	unsigned int delta = 0;
	vector<pair<int, int> > delta_solutions;
	search(read,output, max_errors, delta, delta_solutions, false);
}


/* static */
bool Hash::search_gapped_item_sort(const ResultItem &left, const ResultItem &right) {
	return left.globalPosition < right.globalPosition;
}

void Hash::search_gapped(const string & read, ItemsGapped & output, t_pattern_length seed_size, int seed_errors, int max_errors, t_length max_gap) {
	if (read.length() < 2*seed_size)
		return;

	string for_left_sequence = read.substr(0,seed_size);
	string for_right_sequence = read.substr(read.length()-1-seed_size,seed_size);

	unsigned int delta = 0;
	vector<pair<int, int> > delta_solutions;

	Items temp_alignments;
	Items for_left_alignments;
	Items for_right_alignments;
	Items rev_left_alignments;
	Items rev_right_alignments;

	search(for_left_sequence,temp_alignments,seed_errors,delta,delta_solutions,false);
	for (Items::iterator iter = temp_alignments.begin(); iter != temp_alignments.end(); iter++)
		if (iter->strand == true)
			for_left_alignments.push_back(*iter);
		else
			rev_right_alignments.push_back(*iter);
	temp_alignments.clear();
	search(for_right_sequence,temp_alignments,seed_errors,delta,delta_solutions,false);
	for (Items::iterator iter = temp_alignments.begin(); iter != temp_alignments.end(); iter++)
		if (iter->strand == true)
			for_right_alignments.push_back(*iter);
		else
			rev_left_alignments.push_back(*iter);

	if (for_left_alignments.size() > 0 and for_right_alignments.size() > 0) {
		sort(for_left_alignments.begin(), for_left_alignments.begin(), search_gapped_item_sort);
		sort(for_right_alignments.begin(), for_right_alignments.begin(), search_gapped_item_sort);

		max_errors = stranded_gap_search_subroutine(output,for_left_alignments,for_right_alignments,max_gap,read.c_str(),seed_size,seed_size,max_errors);

	}

	if ((rev_left_alignments.size() > 0 and rev_right_alignments.size() > 0)) {
		string rev_left_sequence = reverse_complement_standalone_str(for_right_sequence.c_str());
		string rev_right_sequence = reverse_complement_standalone_str(for_left_sequence.c_str());
		string rev_read = reverse_complement_standalone_str(read.c_str());

		sort(rev_left_alignments.begin(), rev_left_alignments.begin(), search_gapped_item_sort);
		sort(rev_right_alignments.begin(), rev_right_alignments.begin(), search_gapped_item_sort);

		max_errors = stranded_gap_search_subroutine(output,rev_left_alignments,rev_right_alignments,max_gap,rev_read.c_str(),seed_size,seed_size,max_errors);
	}

}

/*
#ifdef HAVE_MPI

void Hash::searchMPI(const string & read, Items & output, unsigned int max_errors, int maskSessionNumber, MPI::Win &win, mutex &m) {

	int patternLength = read.length();
	if (patternLength <= (int) blockLength) {
		return;
	}
	int numberOfBlocks = patternLength / blockLength;
	int blockErrors;
	(max_errors / numberOfBlocks > 2) ? blockErrors = 2 : blockErrors = max_errors / numberOfBlocks;

	int pattern[2][patternLength];
	char PATTERN[2][patternLength];

	unsigned int p;
	unsigned int l;
	int randNum;
	// processing the read
	int maxErrorInRead = max_errors;
	for (int plTemp = 0; plTemp < patternLength; plTemp++) {
		randNum = 1; //(int)(rand()%4);
		switch (read.at(plTemp)) {
		case 'A': case 'a': pattern[0][plTemp] = 0;
		pattern[1][patternLength - plTemp - 1] = 3;
		PATTERN[0][plTemp] = 'A';
		PATTERN[1][patternLength - plTemp - 1] = 'T';
		break;
		case 'C': case 'c': pattern[0][plTemp] = 1;
		pattern[1][patternLength - plTemp - 1] = 2;
		PATTERN[0][plTemp] = 'C';
		PATTERN[1][patternLength - plTemp - 1] = 'G';
		break;
		case 'G': case 'g': pattern[0][plTemp] = 2;
		pattern[1][patternLength - plTemp - 1] = 1;
		PATTERN[0][plTemp] = 'G';
		PATTERN[1][patternLength - plTemp - 1] = 'C';
		break;
		case 'T': case 't': pattern[0][plTemp] = 3;
		pattern[1][patternLength - plTemp - 1] = 0;
		PATTERN[0][plTemp] = 'T';
		PATTERN[1][patternLength - plTemp - 1] = 'A';
		break;
		default: pattern[0][plTemp] = randNum;
		pattern[1][patternLength - plTemp - 1] = randNum;
		PATTERN[0][plTemp] = 'N';
		PATTERN[1][patternLength - plTemp - 1] = 'N';
		}

	}

	bool foundExact = false; // IF AN OCCURENCE IS FOUND WITH 0 ERRORS DON'T SEARCH FOR OTHER OCCURENCES WITH MORE THEN ONE ERROR

	unsigned long int Pq;
	int numN;

	int errors = blockErrors;
	unsigned int diffTemp;
	int idxZ = 0;
	int dTemp = 0;
	unsigned long int idxText = 0;
	int j;
	for (int pp = 0; pp <= 1; pp++) { // reverse complement
		foundExact = false;
		for (int ii = 0; ii < numberOfBlocks && !foundExact; ii++) {

			numN = 0;
			for (j = ii * blockLength, Pq = 0; j < (ii + 1) * (int) blockLength; j++) {
				Pq = ((Pq << 2) + pattern[pp][j]);
				if (Pq >= q) {
					Pq = Pq % q;
				}
				if (PATTERN[pp][j] == 'N') numN++;
			}

			if (numN <= errors) {
				for (idxZ = 0; idxZ < errorsInZ[errors] && !foundExact; idxZ++) { // for each z in Z0

					(Pq >= Z[idxZ]) ? idxText = Pq - Z[idxZ] : idxText = q - Z[idxZ] + Pq;
					p = HASHcounter[idxText];
					(idxText == q - 1) ? l = LAST : l = HASHcounter[idxText + 1];

					while (p < l) {
						if (HASHvalues[p] >= ii * blockLength) {
							diffTemp = HASHvalues[p] - ii*blockLength;
							if (diffTemp + patternLength < textLength) {// && Checked.insert(diffTemp).second) { // we check that ideed *iter - i * pattern_length is a valid shift in the text
								dTemp = 0;
								for (int j = 0; dTemp <= maxErrorInRead && (j < (int) patternLength); j++) {
									if (PATTERN[pp][j] != TEXT[diffTemp + j]) {
										dTemp++;
										if (TEXT[diffTemp + j] == text_delimitator) {
											dTemp = max_errors + 1;
										}
									}
								}

								if (dTemp <= maxErrorInRead) { // a good hit
									ResultItem K;
								K.GlobalPosition = diffTemp;
								K.errors = dTemp;
								K.strand = (pp == 0) ? true : false;
								K.indels = false;

								if (dTemp == maxErrorInRead) {
									output.push_back(K);
								} else {
									errors = dTemp / numberOfBlocks;
									if (errors > 2) errors = 2;
									maxErrorInRead = dTemp;
									output.clear();
									output.push_back(K);
									if (dTemp == 0) foundExact = true;
								}
								}
							}
						}
						p++;
					}

					if ((errorsInZ[errors] > 1000000000 && idxZ % 100000000 == 0) ||
							(errorsInZ[errors] > 100000000 && idxZ % 10000000 == 0) ||
							(errorsInZ[errors] > 10000000 && idxZ % 1000000 == 0) ||
							(errorsInZ[errors] > 1000000 && idxZ % 100000 == 0) ||
							(errorsInZ[errors] > 100000 && idxZ % 10000 == 0) ||
							(errorsInZ[errors] > 10000 && idxZ % 1000 == 0) ||
							idxZ == errorsInZ[errors] - 1 ||
							idxZ == 0) {


						//Sequence alignment info
						//0 = number of mismatch
						//1 = number of alignments
						int seqInfo[2];

						{

							//Only one thread by time can read from Master
							mutex::scoped_lock lock(m);

							//Update alignment info
							//Lock memory in shared mode
							win.Lock(MPI::LOCK_SHARED, 0, 0);

							//Get info
							win.Get(&seqInfo, 2 * sizeof (int), MPI::BYTE, 0, (sizeof (Mask::MaskData) * maskSessionNumber) + offsetof(Mask::MaskData, NM), 2 * sizeof (int), MPI::BYTE);

							//Unlock shared memory
							win.Unlock(0);
						}


						//If exact occurrence was found...
						if (seqInfo[0] == 0 && seqInfo[1] > 0) {

							//Exit from alignment
							output.clear();
							break;

						} else if (seqInfo[0] < errors && seqInfo[1] > 0) { //Else if best occurence was found...

							//Update alignment parameters
							errors = seqInfo[0];
							idxZ = 0;
						}
					}
				}
			}
		}
	}
}

#endif
*/
