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

#include "BAM1.h"

namespace samobjects {

BAM1 & BAM1::operator=(const BAM1 & src) {
	if (this != &src) { // protect against invalid self-assignment
		copy(src);
	}
	return *this;
}

void BAM1::copy(const BAM1 & src) {
	clear();
	if (&src != NULL and src.bam != NULL)
		bam = bam_dup1(src.bam);
}

void BAM1::copy(const bam1_t * src) {
	clear();
	if (src != NULL)
		bam = bam_dup1(src);
}

SAM_istream& operator>>(SAM_istream& in, BAM1 & bam) throw (SAM_IO_Error) {
	bam.clear();
	bam.bam = in.read();
	return in;
}

SAM_ostream& operator<<(SAM_ostream& out, const BAM1 & bam) throw (SAM_IO_Error) {
	if (bam.pointer() != NULL) {
		if (bam.pointer()->data == NULL) {
			throw SAM_IO_Error(SAM_IO_Error::incorrect_bam1_t,"data is null");
		}
		out.write(bam.bam);
	}
	return out;
}

BAM1::Tag_Result BAM1::find_and_read_tag(const char tag[2]) const {
	Tag_Result tr;
	find_and_read_tag(tag,tr);
	return tr;
}


void BAM1::find_and_read_tag(const char tag[2], Tag_Result & result) const {
	find_and_read_tag(bam,tag,result);
}

/*static*/
void BAM1::find_and_read_tag(bam1_t * bam, const char tag[2], Tag_Result & result) {
	if (bam == NULL) {
		result.t = not_found;
		return;
	}
	uint8_t * p = bam_aux_get(bam, tag);
	if (p == NULL) {
		result.t = not_found;
		return;
	}
	switch (p[0]) {
	case 'Z':
	case 'H':
		result.t = string_t;
		result.d.string_pointer = bam_aux2Z(p);
		return;
	case 'A':
		result.t = char_t;
		result.d.char_value = bam_aux2A(p);
		return;
	case 'c':
	case 'C':
	case 's':
	case 'S':
	case 'i':
	case 'I':
		result.t = integer_t;
		result.d.integer_value = bam_aux2i(p);
		return;
	case 'f':
		result.t = double_t;
		result.d.double_value = bam_aux2f(p);
		return;
	case 'd':
		result.t = double_t;
		result.d.double_value = bam_aux2d(p);
		return;
	default:
		result.t = unknown_t;
		return;
	}
}

bool BAM1::modify_tag_int32_t(const char tag[2], int32_t new_value) {
	return modify_tag_int32_t(bam,tag,new_value);
}


bool BAM1::modify_tag_int32_t(bam1_t *bam, const char tag[2], int32_t new_value) {

	if (bam == NULL)
		return false;
	uint8_t * p = bam_aux_get(bam, tag);
	if (p == NULL)
		return false;
	size_t old_size;
	switch (p[0]) {
	case 'c':
	case 'C':
		old_size = 2;
		break;
	case 's':
	case 'S':
		old_size = 3;
		break;
	case 'i':
	case 'I':
		old_size = 5;
		break;
	default:
		return false;
	}
	string value = bam_integer(new_value);
	size_t new_size = value.size();

	size_t new_m_data_size = bam->data_len - old_size + new_size +32;
	uint8_t * end = bam->data + bam->data_len;
	uint8_t * old = bam->data;
	uint8_t * n = ((uint8_t*)malloc(new_m_data_size));
	uint8_t * s = n;

	// copy before
	while (old < p)
		*(s++) = *(old++);

	// copy new value
	for (size_t i = 0; i < new_size; i++)
		*(s++) = value[i];

	// skip old value
	for (size_t i = 0; i < old_size; i++)
		old++;

	// copy after
	while (old < end)
		*(s++) = *(old++);

	free(bam->data);
	bam->data = n;
	bam->data_len = bam->data_len - old_size + new_size;
	bam->l_aux = bam->l_aux - old_size + new_size;
	bam->m_data = new_m_data_size;

	return true;
}

BAM1::t_alignment BAM1::alignment_type() const {
	Tag_Result r;
	find_and_read_tag("NH",r);
	if (r.type() == integer_t)
		switch (r.asInt()) {
		case 0:
			return alignments_not_found;
		case 1:
			return unique_alignment;
		default:
			return multiple_alignments;
		}
	return unknown_alignment;
}

size_t BAM1::aligned_region_size() const {
	if (bam == NULL)
		return 0;
	size_t size = 0;
	uint32_t n_cigar = bam->core.n_cigar;
	uint32_t * cigar = bam1_cigar(bam);

	for (uint32_t i = 0; i < n_cigar; i++) {
		uint32_t s = cigar[i] >> 4;
		switch (cigar[i] & 0x0000000f) {
		case 0: //M
		case 2: //D
		case 3: //N
			size += s;
			break;
		case 1: //I
		case 4: //S
		case 5: //H
		case 6: //P
		default:
			break;
		}
	}
	return size;
}

string BAM1::get_sequence() const {
	if (bam == NULL)
		return string();
	uint8_t * seq =  bam1_seq(bam);
	int32_t l = bam->core.l_qseq;
	string ts;
	ts.reserve(l+1);
	for (int32_t i = 0; i < l; i++)
		switch (bam1_seqi(seq, i)) {
		case 1:
			ts.push_back('A');
			break;
		case 2:
			ts.push_back('C');
			break;
		case 4:
			ts.push_back('G');
			break;
		case 8:
			ts.push_back('T');
			break;
		default:
			ts.push_back('N');
			break;
		}
	return ts;
}

string BAM1::get_quality() const {
	return bam == NULL ? string() : string((char *)bam1_qual(bam),bam->core.l_qseq);
}

string BAM1::get_aligned_sequence() const {
	if (bam == NULL)
		return string();
	uint8_t * seq =  bam1_seq(bam);
	int32_t l = bam->core.l_qseq - trimmed_right_size();
	string ts;
	ts.reserve(l+1);
	for (int32_t i = trimmed_left_size(); i < l; i++)
		switch (bam1_seqi(seq, i)) {
		case 1:
			ts.push_back('A');
			break;
		case 2:
			ts.push_back('C');
			break;
		case 4:
			ts.push_back('G');
			break;
		case 8:
			ts.push_back('T');
			break;
		default:
			ts.push_back('N');
			break;
		}
	return ts;
}

size_t BAM1::trimmed_left_size() const {
	if (bam == NULL)
		return 0;
	size_t size = 0;
	uint32_t n_cigar = bam->core.n_cigar;
	uint32_t * cigar = bam1_cigar(bam);
	bool end = false;
	for (uint32_t i = 0; not end and i < n_cigar; i++) {
		uint32_t s = cigar[i] >> 4;
		switch (cigar[i] & 0x0000000f) {
		case 0: //M
		case 1: //I
		case 2: //D
			end = true;
			break;
		case 3: //N
		case 4: //S
		case 5: //H
		case 6: //P
		default:
			size+=s;
			break;
		}
	}
	return size;
}

size_t BAM1::trimmed_right_size() const {
	if (bam == NULL)
		return 0;
	size_t size = 0;
	uint32_t n_cigar = bam->core.n_cigar;
	uint32_t * cigar = bam1_cigar(bam);
	bool end = false;
	for (uint32_t i = n_cigar; not end and i > 0; i++) {
		uint32_t s = cigar[i-1] >> 4;
		switch (cigar[i] & 0x0000000f) {
		case 0: //M
		case 1: //I
		case 2: //D
			end = true;
			break;
		case 3: //N
		case 4: //S
		case 5: //H
		case 6: //P
		default:
			size+=s;
			break;
		}
	}
	return size;
}

void BAM1::revcomp() {
	if (bam == NULL)
		return;
	// QUALITY
	string old_qual = get_quality();
	const char * original = old_qual.c_str();
	uint8_t * converted = bam1_qual(bam);
	size_t size = old_qual.size();
	for (size_t i = 0; i < size; i++)
		converted[size-1-i] = original[i];

	// SEQUENCE
	string new_seq = compact_DNA(get_sequence(),reverse_strand);
	size = new_seq.size();
	converted = bam1_seq(bam);
	for (size_t i = 0; i < size; i++)
		converted[i] = new_seq[i];

	// CIGAR
	uint32_t * cigar = bam1_cigar(bam);
	uint32_t n = bam->core.n_cigar;
	uint32_t temp;
	for (uint32_t i = 1; i <= (n/2); i++) {
		temp = cigar[i-1];
		cigar[i-1] = cigar[n-i];
		cigar[n-i] = temp;
	}
}

/* static */ string BAM1::compact_DNA(const string & sequence, t_strand strand) {
	string compact;
	const char * dna = sequence.c_str();
	size_t l = sequence.size();
	compact.clear();
	compact.reserve((size_t)((l+1)/2)+1);
	char first;
	char second;
	if (strand == reverse_strand)
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
	else
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
	return compact;
}

bool BAM1::is_mapped_univocally() const {
	if (is_mapped()) {
		BAM1::Tag_Result result = find_and_read_tag("NH");
		if (result.type() == BAM1::integer_t and result.asInt() != 1)
			return false;
		else
			return true;
	} else
		return false;
}

bool BAM1::is_mapped_multiple() const {
	if (is_mapped()) {
		BAM1::Tag_Result result = find_and_read_tag("NH");
		if (result.type() == BAM1::integer_t and result.asInt() != 1)
			return true;
		else
			return false;
	} else
		return false;
}

int BAM1::NH(bam1_t * bam){
	Tag_Result r;
	find_and_read_tag(bam,"NH",r);
	if (r.type() == integer_t)
		return r.asInt();
	else
		return -1;
}

int BAM1::HI(bam1_t * bam){
	Tag_Result r;
	find_and_read_tag(bam,"HI",r);
	if (r.type() == integer_t)
		return r.asInt();
	else
		return -1;
}

int BAM1::IH(bam1_t * bam){
	Tag_Result r;
	find_and_read_tag(bam,"IH",r);
	if (r.type() == integer_t)
		return r.asInt();
	else
		return -1;
}

}

