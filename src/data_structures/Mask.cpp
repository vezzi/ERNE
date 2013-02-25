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

#include "Mask.h"

Mask::~Mask() {
//	this->DELTA.~vector();
//	this->cigarVector.~vector();

}

Mask::Mask() {
	good_region_start = 0;
	good_region_stop = 0;
	globalPosition = 0;
	algn = 0;
	NM = 0;
	NM_gap = 0;
	low_complexity = false;
	low_quality = false;
	discarded = false;
	masked = false;
	gapped = false;
	contaminated = false;
	strand = true; // needed
	DELTA.clear();
	primary = true;
}


Mask::Mask(Fasta const & fasta) {
	good_region_start = 0;
	good_region_stop = 0;
	globalPosition = 0;
	algn = 0;
	NM = 0;
	NM_gap = 0;
	masked = false;
	low_complexity = false;
	low_quality = false;
	discarded = false;
	sequence = fasta.get_sequence();
	quality = fasta.get_quality();
	gapped = false;
	contaminated = false;
	strand = true; // needed
	DELTA.clear();
	primary = true;
}

ostream& operator<<(ostream& channel, const Mask& mask) {
	if (mask.masked) {
		t_pattern_length l = mask.sequence.length();
		for (t_pattern_length i=1; i <= l; i++)
			if ((i < mask.good_region_start) or (i > mask.good_region_stop))
				channel << (char)(mask.sequence[i-1] - 'A' + 'a');
			else
				channel << mask.sequence[i-1];
	} else
		channel << mask.sequence;
	return channel;
}

string Mask::to_string() const {
	stringstream ss;
	ss << *this;
	return ss.str();
}

string Mask::get_good_sequence() const {
	if (discarded)
		return string();
	if (masked)
		return sequence.substr(good_region_start-1,good_region_stop-good_region_start+1);
	else
		return sequence;
}

string Mask::get_good_quality() const {
	if (discarded)
		return string();
	if (masked)
		return quality.substr(good_region_start-1,good_region_stop-good_region_start+1);
	else
		return quality;
}


string Mask::get_masked() const {
	if (discarded) {
		string temp = sequence;
		for (t_pattern_length i = 0; i < sequence.length(); i++)
			temp[i] = tolower(temp[i]);
		return temp;
	} else if (masked) {
		string temp = sequence;
		for (t_pattern_length i = 1; i < good_region_start; i++)
			temp[i-1] = tolower(temp[i-1]);
		for (t_pattern_length i = good_region_stop+1; i <= temp.size(); i++)
			temp[i-1] = tolower(temp[i-1]);
		return temp;
	} else
		return sequence;
}

int Mask::get_five_prime() {
	if (strand) {
		if (masked)
			return position - good_region_start + 1;
		else
			return position;
	} else {
		if (masked)
			return position - good_region_stop - 1;
		else
			return position + sequence.size() - 1;
	}
}

string Mask::get_MD(const char * reference) {
	stringstream MD;
	int e = 0;

	string for_or_rev_string;
	if (strand)
		for_or_rev_string = get_good_sequence();
	else
		for_or_rev_string = reverse_complement_standalone_str(get_good_sequence());

	const char * for_or_rev = for_or_rev_string.c_str();
	int l = get_good_length();
	reference += globalPosition;
	for (int i = 0; i < l ; i++) {
		if (*reference == *for_or_rev)
			e++;
		else {
			if (e > 0)
				MD << e;
			MD << *for_or_rev;
			e = 0;
		}
		reference++;
		for_or_rev++;
	}

	if (e != 0)
		MD << e;
	return MD.str();
}

Mask::t_CIGAR Mask::get_CIGAR() {
	if (gapped)
		return get_CIGAR_GAP();
	else if (indels)
		return get_CIGAR_SW();
	else
		return get_CIGAR_NORMAL();
}

/* INLINE */
void Mask::uint_and_type_to_CIGAR_stringstream(stringstream & ss, unsigned int i, unsigned short int type) {
	i = (i << 4) + type;
	ss << (char)(i%256);
	i = i / 256;
	ss << (char)(i%256);
	i = i / 256;
	ss << (char)(i%256);
	i = i / 256;
	ss << (char)(i%256);
}

Mask::t_CIGAR Mask::get_CIGAR_NORMAL() {
	stringstream CIGAR;
	unsigned int operations = 0;
	if (strand) {
		if (masked and good_region_start > 1) {
			uint_and_type_to_CIGAR_stringstream(CIGAR ,good_region_start-1, BAM_CSOFT_CLIP);
			operations++;
		}
		uint_and_type_to_CIGAR_stringstream(CIGAR, get_good_length(), BAM_CMATCH);
		operations++;
		if (masked and good_region_stop < sequence.size()) {
			uint_and_type_to_CIGAR_stringstream(CIGAR, sequence.size() - good_region_stop, BAM_CSOFT_CLIP);
			operations++;
		}
	} else {
		if (masked and good_region_stop < sequence.size()) {
			uint_and_type_to_CIGAR_stringstream(CIGAR, sequence.size() - good_region_stop, BAM_CSOFT_CLIP);
			operations++;
		}
		uint_and_type_to_CIGAR_stringstream(CIGAR, get_good_length(), BAM_CMATCH);
		operations++;
		if (masked and good_region_start > 1) {
			uint_and_type_to_CIGAR_stringstream(CIGAR, good_region_start-1, BAM_CSOFT_CLIP);
			operations++;
		}
	}
	return t_CIGAR(CIGAR.str(),operations);
}


Mask::t_CIGAR Mask::get_CIGAR_SW() {
	stringstream CIGAR;
	unsigned int operations = 0;
	if (strand) {
		if (masked and good_region_start > 1) {
			uint_and_type_to_CIGAR_stringstream(CIGAR, good_region_start-1, BAM_CSOFT_CLIP);
			operations++;
		}
		//I have to perform SW in order to compute best score and backtrack to create CIGAR
		for(unsigned int k=0; k< cigarVector.size(); k++) {
			switch( cigarVector.at(k).first) {
			case 'M':
				uint_and_type_to_CIGAR_stringstream(CIGAR, cigarVector.at(k).second, BAM_CMATCH);
				break;
			case 'D':
				uint_and_type_to_CIGAR_stringstream(CIGAR, cigarVector.at(k).second, BAM_CDEL);
				break;
			case 'I':
				uint_and_type_to_CIGAR_stringstream(CIGAR, cigarVector.at(k).second, BAM_CINS);
				break;
			}
			operations++;
		}
		if (masked and good_region_stop < sequence.size()) {
			uint_and_type_to_CIGAR_stringstream(CIGAR, sequence.size() - good_region_stop, BAM_CSOFT_CLIP);
			operations++;
		}
	} else {
		if (masked and good_region_stop < sequence.size()) {
			uint_and_type_to_CIGAR_stringstream(CIGAR, sequence.size() - good_region_stop, BAM_CSOFT_CLIP);
			operations++;
		}

		for(unsigned int k=0; k< cigarVector.size(); k++) {
			switch( cigarVector.at(k).first) {
			case 'M':
				uint_and_type_to_CIGAR_stringstream(CIGAR, cigarVector.at(k).second, BAM_CMATCH);
				break;
			case 'D':
				uint_and_type_to_CIGAR_stringstream(CIGAR, cigarVector.at(k).second, BAM_CDEL);
				break;
			case 'I':
				uint_and_type_to_CIGAR_stringstream(CIGAR, cigarVector.at(k).second, BAM_CINS);
				break;
			}
			operations++;
		}
		if (masked and good_region_start > 1) {
			uint_and_type_to_CIGAR_stringstream(CIGAR, good_region_start-1, BAM_CSOFT_CLIP);
			operations++;
		}
	}
	return t_CIGAR(CIGAR.str(),operations);
}

Mask::t_CIGAR Mask::get_CIGAR_GAP() {
	stringstream CIGAR;
	unsigned int operations = 0;
	if (strand) {
		if (masked and good_region_start > 1) {
			uint_and_type_to_CIGAR_stringstream(CIGAR, good_region_start-1, BAM_CSOFT_CLIP);
			operations++;
		}

		uint_and_type_to_CIGAR_stringstream(CIGAR, length1_gap, BAM_CMATCH);
		uint_and_type_to_CIGAR_stringstream(CIGAR, position_gap - position - length1_gap + 1, BAM_CREF_SKIP);
		uint_and_type_to_CIGAR_stringstream(CIGAR, length2_gap, BAM_CMATCH);
		operations+=3;

		if (masked and good_region_stop < sequence.size()) {
			uint_and_type_to_CIGAR_stringstream(CIGAR, sequence.size() - good_region_stop, BAM_CSOFT_CLIP);
			operations++;
		}
	} else {
		if (masked and good_region_stop < sequence.size()) {
			uint_and_type_to_CIGAR_stringstream(CIGAR, sequence.size() - good_region_stop, BAM_CSOFT_CLIP);
			operations++;
		}

		uint_and_type_to_CIGAR_stringstream(CIGAR, length1_gap, BAM_CMATCH);
		uint_and_type_to_CIGAR_stringstream(CIGAR, position_gap - position - length1_gap +1 , BAM_CREF_SKIP);
		uint_and_type_to_CIGAR_stringstream(CIGAR, length2_gap, BAM_CMATCH);
		operations+=3;

		if (masked and good_region_start > 1) {
			uint_and_type_to_CIGAR_stringstream(CIGAR, good_region_start-1, BAM_CSOFT_CLIP);
			operations++;
		}
	}
	return t_CIGAR(CIGAR.str(),operations);
}

void Mask::extract_fasta_masked(Fasta & output) const {
	output.set_sequence(get_good_sequence().c_str());
	output.set_quality(get_good_quality().c_str());
}


void Mask::set_fasta_masked(Fasta & output) const {
	output.set_id(get_id().c_str());
	output.set_sequence(get_good_sequence().c_str());
	output.set_quality(get_good_quality().c_str());
}

void Mask::set_original_fasta(Fasta & output) const {
	output.set_id(get_id().c_str());
	output.set_sequence(sequence.c_str());
	output.set_quality(quality.c_str());
}


string Mask::print_information() const {
	stringstream ss;
	if (masked)
		ss << "=== sequence is masked === " << endl;
	else
		ss << "=== sequence is NOT masked ===" << endl;
	ss << "REAL_SEQUENCE:\t" << sequence << endl;
	ss << "REAL_QUALITY:\t" << quality << endl;
	ss << "GOOD SEQUENCE:\t" << get_good_sequence() << endl;
	ss << "GOOD QUALITY:\t" << get_good_quality() << endl;
	if (masked)
		ss << "Sequence is MASKED" << endl;
	if (low_complexity)
		ss << "Sequence is LOW COMPLEXITY" << endl;
	if (low_quality)
		ss << "Sequence is LOW QUALITY" << endl;
	if (trimmed)
		ss << "Sequence is TRIMMED" << endl;
	if (discarded)
		ss << "Sequence is DISCARDED" << endl;
	if (contaminated)
		ss << "Sequence is CONTAMINATED" << endl;
	ss << "5':\t " << good_region_start << endl;
	ss << "3':\t " << good_region_stop << endl;

	return ss.str();
}

char * Mask::quality_conversion() {
	const char * original = quality.c_str();
	size_t size = quality.size();
	char * quality_converted = new char[size+1];
	quality_converted[size]= '\0';
	if (strand) {
		for (size_t i = 0; i < size; i++)
			quality_converted[i] = original[i] - 33;
	} else {
		for (size_t i = 0; i < size; i++)
			quality_converted[size-i-1] = original[i] - 33;
	}
	return quality_converted;
}

void Mask::quality_trimming_MOTT(t_min_phred_value_MOTT min_phred_value_MOTT, t_min_phred_value_MOTT min_mean_quality, t_min_phred_value_MOTT min_size) {
	const char *q = quality.c_str();
	size_t l = quality.size();

	// finding left start position
	size_t left = 0;
	bool not_found = true;
	while (not_found and (left < l)) {
		if ((q[left] - 33) >= min_phred_value_MOTT)
			not_found = false;
		else
			left++;
	}

	if (not_found) {
		discarded = true;
		low_quality = true;
	} else {
		// finding right end position
		int a = q[left] - 33 - min_phred_value_MOTT;
		int max = a;
		size_t right = left;
		unsigned int mean_sum = q[left] - 33;
		unsigned int mean_sum_temp = q[left] - 33;
		for (size_t i = left + 1; i < l; i++) {
			a += q[i] - 33 - min_phred_value_MOTT;
			mean_sum_temp += q[i] - 33;
			if (a < 0)
				a = 0;
			if (max <= a) {
				max = a;
				right = i;
				mean_sum = mean_sum_temp;
			}
		}

		const char *s = sequence.c_str();
		while ((left > 0) and (left < l) and s[left] == 'N')
			left++;
		while ((right > 0) and (right < l) and s[right] == 'N')
			right--;

		if ((right - left + 1 < min_size) or (mean_sum / (right - left + 1) < min_mean_quality)) {
			discarded = true;
			low_quality = true;
		} else {
			// good region is [left,right]
			masked = true;
			good_region_start = left+1;
			good_region_stop = right+1;
		}
	}

}

void Mask::compact_DNA(string & compact) {
	const char * dna = sequence.c_str();
	size_t l = sequence.size();
	compact.clear();
	compact.reserve((size_t)((l+1)/2)+1);
	char first;
	char second;
	if (strand)
		for (size_t i = 0; i < l; i+=2) {
			switch (dna[i]) {
			case 'A' :
			case 'a' : first = 1; break;
			case 'C' :
			case 'c' : first = 2; break;
			case 'G' :
			case 'g' : first = 4; break;
			case 'T' :
			case 't' : first = 8; break;
			default  : first = 15; break;
			}
			switch ((i+1)<l ? dna[i+1] : 'N') {
			case 'A' :
			case 'a' : second = 1; break;
			case 'C' :
			case 'c' : second = 2; break;
			case 'G' :
			case 'g' : second = 4; break;
			case 'T' :
			case 't' : second = 8; break;
			default  : second = 15; break;
			}
			compact.push_back((char)(first*16+second));
		}
	else
		for (size_t i = 0; i < l; i+=2) {
			switch (dna[l-i-1]) {
			case 'A' :
			case 'a' : first = 8; break;
			case 'C' :
			case 'c' : first = 4; break;
			case 'G' :
			case 'g' : first = 2; break;
			case 'T' :
			case 't' : first = 1; break;
			default  : first = 15; break;
			}
			switch ((i+1)<l ? dna[l-i-2] : 'N') {
			case 'A' :
			case 'a' : second = 8; break;
			case 'C' :
			case 'c' : second = 4; break;
			case 'G' :
			case 'g' : second = 2; break;
			case 'T' :
			case 't' : second = 1; break;
			default  : second = 15; break;
			}
			compact.push_back((char)(first*16+second));
		}
}

void Mask::toData(MaskData &sd) {

	//Copy id
	if (id.length() > MAX_READ_LENGTH) {
		cerr << "ID " << id << " too long. Max allowed read size is " << MAX_READ_LENGTH << endl;
		exit(5);
	} else {
		strcpy(sd.id, id.c_str());
	}

	//Copy sequence and quality
	if (sequence.length() > MAX_READ_LENGTH) {
		cerr << "Read " << id << " too long. Max allowed read size is " << MAX_READ_LENGTH << endl;
		exit(5);
	} else {
		strcpy(sd.sequence, sequence.c_str());
		strcpy(sd.quality, quality.c_str());
	}

	sd.position = position;
	sd.position_gap = position_gap;

	sd.globalPosition = globalPosition;
	sd.globalPosition_gap = globalPosition_gap;

	sd.contig = contig;
	sd.length1_gap = length1_gap;
	sd.length2_gap = length2_gap;

	sd.NM = NM;
	sd.NM_gap = NM_gap;

	sd.algn = algn;

	sd.good_region_start = good_region_start;
	sd.good_region_stop = good_region_stop;

	sd.primary = primary;
	sd.IH = IH;
	sd.HI = HI;

	sd.strand = strand;
	sd.masked = masked;
	sd.low_complexity = low_complexity;
	sd.low_quality = low_quality;
	sd.trimmed = trimmed;
	sd.discarded = discarded;
	sd.contaminated = contaminated;
	sd.gapped = gapped;
}

void Mask::updateFromData(const MaskData &sd) {

	//Already in C legacy form
	id = string(sd.id);
	sequence = string(sd.sequence);
	quality = string(sd.quality);

	position = sd.position;
	position_gap = sd.position_gap;

	globalPosition = sd.globalPosition;
	globalPosition_gap = sd.globalPosition_gap;

	contig = sd.contig;
	length1_gap = sd.length1_gap;
	length2_gap = sd.length2_gap;

	//This two entries MUST be in this order and in sequence! They are used in error communication between Workers
	NM = sd.NM;
	algn = sd.algn;

	NM_gap = sd.NM_gap;

	good_region_start = sd.good_region_start;
	good_region_stop = sd.good_region_stop;

	primary = sd.primary;
	IH = sd.IH;
	HI = sd.HI;


	strand = sd.strand;
	masked = sd.masked;
	low_complexity = sd.low_complexity;
	low_quality = sd.low_quality;
	trimmed = sd.trimmed;
	discarded = sd.discarded;
	contaminated = sd.contaminated;
	gapped = sd.gapped;
}

void Mask::adjust_transcriptome_coordinates(const char * reference, size_t reference_length) {



}

char complement(char c){

	switch(c){

		case 'A' : case 'a' : return 'T';
		case 'T' : case 't' : return 'A';
		case 'C' : case 'c' : return 'G';
		case 'G' : case 'g' : return 'C';
		default : return 'N';

	}

}

void Mask::reverse_complement(){

	string sequence_rc;
	string quality_rev;

	for(int i = sequence.length()-1;i>=0;i--){
		sequence_rc.push_back(complement(sequence.at(i)));
	}

	for(int i = quality.length()-1;i>=0;i--){
		quality_rev.push_back(quality.at(i));
	}

	sequence = sequence_rc;
	quality = quality_rev;

}
