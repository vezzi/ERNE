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

#ifndef BAM1_H_
#define BAM1_H_

#include "SAM_istream.h"
#include "SAM_ostream.h"

#include "common.h"

namespace samobjects {

class BAM1 {
public:

	enum t_alignment { unique_alignment, multiple_alignments, alignments_not_found, quality_discarded, low_complexity, contamination, unknown_alignment };

	enum tag_t { not_found, char_t, string_t, integer_t, double_t, unknown_t};
	union data_t {
		int32_t integer_value;
		double double_value;
		const char * string_pointer;
		char char_value;
	};


	static int NH(bam1_t * bam);
	static int HI(bam1_t * bam);
	static int IH(bam1_t * bam);

	class Tag_Result {
	public:

		Tag_Result() { }
		Tag_Result(BAM1 & bam, const char tag[2]) { bam.find_and_read_tag(tag,*this); }
		Tag_Result(BAM1 * bam, const char tag[2]) { bam->find_and_read_tag(tag,*this); }
		tag_t type() { return t; }

		int32_t asInt() const { return d.integer_value; }
		double asDouble() const { return d.double_value; }
		const char * asString() const { return d.string_pointer; }
		char asChar() const { return d.char_value; }

		friend class BAM1;
	protected:
		tag_t t;
		data_t d;
	};

	BAM1() { bam = NULL; }
	BAM1(const bam1_t *src) { bam = NULL; copy(src); }
	BAM1(const BAM1 & src) { bam = NULL; copy(src); }
	virtual ~BAM1() { clear(); }

	friend SAM_istream& operator>>(SAM_istream& in, BAM1 & bam) throw (SAM_IO_Error);
	friend SAM_ostream& operator<<(SAM_ostream& out, const BAM1 & bam) throw (SAM_IO_Error);
	BAM1& operator=(const BAM1 & src);

	void copy(const bam1_t * bam);
	void copy(const BAM1 & src);

	bam1_t * pointer() const { return bam; }
	bam1_core_t * core() const { return bam == NULL ? NULL : &(bam->core); }
	Tag_Result find_and_read_tag(const char tag[2]) const;

	void find_and_read_tag(const char tag[2], Tag_Result & result) const;
	static void find_and_read_tag(bam1_t * bam, const char tag[2], Tag_Result & result);

	bool modify_tag_int32_t(const char tag[2], int32_t new_value);
	static bool modify_tag_int32_t(bam1_t * bam, const char tag[2], int32_t new_value);

	bool empty() const { return bam == NULL ? true : bam->data == NULL; }
	string get_aligned_sequence() const;
	size_t aligned_region_size() const;
	size_t trimmed_left_size() const;
	size_t trimmed_right_size() const;

	t_alignment alignment_type() const;

	bool is_paired() const { return bam == NULL ? false : bam->core.flag & BAM_FPAIRED; }
	bool is_proper_pair() const { return bam == NULL ? false : bam->core.flag & BAM_FPROPER_PAIR; }
	bool is_first_read() const { return bam == NULL ? false : bam->core.flag & BAM_FREAD1; }
	bool is_second_read() const { return bam == NULL ? false : bam->core.flag & BAM_FREAD2; }
	bool is_mapped() const { return bam == NULL ? false : not (bam->core.flag & BAM_FUNMAP); }
	bool is_mate_mapped() const { return bam == NULL ? false : not (bam->core.flag & BAM_FMUNMAP); }
	bool is_quality_discarded() const { return bam == NULL ? false : bam->core.flag & BAM_FQCFAIL; }

	bool is_mapped_univocally() const;
	bool is_mapped_multiple() const;

	const char * get_id() const { return bam == NULL ? NULL : bam1_qname(bam); }
	int32_t get_pos() const { return bam == NULL ? -1 : bam->core.pos; }
	int32_t get_mpos() const { return bam == NULL ? -1 : bam->core.mpos; }
	int32_t get_tid() const { return bam == NULL ? -1 : bam->core.tid; }
	int32_t get_mtid() const { return bam == NULL ? -1 : bam->core.mtid; }
	t_strand get_strand() const { return bam == NULL ? unknown_strand : bam->core.flag & BAM_FREVERSE ? reverse_strand : forward_strand; }
	t_strand get_mstrand() const { return bam == NULL ? unknown_strand : bam->core.flag & BAM_FMREVERSE ? reverse_strand : forward_strand; }
	uint32_t get_end() const { return bam == NULL ? 0 : bam->core.pos + aligned_region_size() -1; }
	int32_t get_insert_size() const { return bam == NULL ? 0 : bam->core.isize; }
	int32_t get_insert_size_abs() const { return bam == NULL ? 0 : abs(bam->core.isize); }

	string get_sequence() const;
	string get_quality() const;

	void revcomp();

	static string compact_DNA(const string & sequence, t_strand strand = unknown_strand);

protected:
	bam1_t * bam;

	void clear() {	if (bam != NULL) { bam_destroy1(bam); bam = NULL;} }

};

}

#endif /* BAM1_H_ */
