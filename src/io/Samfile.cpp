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

#include "Samfile.h"

#include <sys/stat.h>

#define BAM_DATA_DEFAULT_LENGTH 1024

#ifdef HAVE_MPI
void Samfile::open_file_mpi(const Options & options) {
	{
		stringstream ss;
		ss << time(NULL);
		reading_group_id = ss.str();
	}

	gui_output = options.gui_output;

	// prepare header
	stringstream header_name;
	header_name << options.reference_file << "_h.dht";

	if (not file_exists(header_name.str().c_str())) {
		ERROR_CHANNEL << "=== Cannot open file " << header_name.str() << endl;
		ERROR_CHANNEL << "=== Check parameters! " << endl;
		exit(3);
	}

	sam_header = sam_header_read2(header_name.str().c_str());

	// adding more header information
	string description = header_description(options);

	if (sam_header->text != NULL)
		free(sam_header->text);
	sam_header->text = (char *) malloc(description.size()+1);
	sam_header->l_text = description.size();
	memcpy(sam_header->text,description.c_str(),description.size()+1);

	// open file (and write header);
	if (options.bam_format)
		sam_file = samopen(options.output_file.c_str(),"wb",sam_header);
	else
		sam_file = samopen(options.output_file.c_str(),"wh",sam_header);

	if (sam_file == NULL) {
		ERROR_CHANNEL << "=== Cannot open file for output" << endl;
		ERROR_CHANNEL << "=== Check parameters! " << endl;
		exit(3);
	}

	if (gui_output) {
		DEFAULT_CHANNEL << "\n>H\t" << sam_header->n_targets << '\n';
		for (int32_t i = 0; i < sam_header->n_targets; i++)
			DEFAULT_CHANNEL << ">H\t" << sam_header->target_name[i] << '\t' << sam_header->target_len[i] << '\n';
	}
}
#endif

void Samfile::open_file(const Options & options, const Hash & H, const Hash & CR, int index) {
	if (reading_group_id.compare(string("")) == 0) {
		stringstream ss;
		ss << time(NULL);
		reading_group_id = ss.str();
	}

	gui_output = options.gui_output;
	if (options.insert_size_check)
		set_insert_size(options.insert_size_min, options.insert_size_max);

	// prepare file for header
	string tmp_name("/tmp/ERNE_temp_file_");
	char tmp_num[100];
	srand ( time(NULL) );
	sprintf(tmp_num,"%d",rand());
	tmp_name.append(tmp_num);
	tmp_name.append("_");
	sprintf(tmp_num,"%d",rand());
	tmp_name.append(tmp_num);
	ofstream o(tmp_name.c_str());
	for (int i = 0; i < H.globalToLocal.contigs ; i++)
		o << H.globalToLocal.contigsNames[i] << '\t' << (H.globalToLocal.endPositions[i]-H.globalToLocal.startPositions[i]+1) << endl;
	o.close();

	// prepare header
	sam_header = sam_header_read2(tmp_name.c_str());

	// delete temp file
	remove(tmp_name.c_str());

	// adding more header information
	string description = header_description(options);

	if (sam_header->text != NULL)
		free(sam_header->text);
	sam_header->text = (char *) malloc(description.size()+1);
	sam_header->l_text = description.size();
	memcpy(sam_header->text,description.c_str(),description.size()+1);

	// open file (and write header);

	string tmp_output_file;
	if (index >= 0) {
		stringstream ss;
		if(options.bam_format)
			ss << options.output_file << "_" << index << ".bam" ;
		else
			ss << options.output_file << "_" << index << ".sam" ;

		tmp_output_file = ss.str();
	} else
		tmp_output_file = options.output_file;

	if (options.bam_format)
		sam_file = samopen(tmp_output_file.c_str(),"wb",sam_header);
	else
		sam_file = samopen(tmp_output_file.c_str(),"wh",sam_header);

	if (sam_file == NULL) {
		ERROR_CHANNEL << "Cannot open file for output" << endl;
		exit(3);
	}

	this->CR = &CR;

	if (gui_output) {
		DEFAULT_CHANNEL << ">H\t" << sam_header->n_targets << '\n';
		for (int32_t i = 0; i < sam_header->n_targets; i++)
			DEFAULT_CHANNEL << ">H\t" << sam_header->target_name[i] << '\t' << sam_header->target_len[i] << '\n';
	}
}

void Samfile::close_file() {
	if (sam_file != NULL) {
		samclose(sam_file);
		bam_header_destroy(sam_header);
		sam_file = NULL;
	}
}

/* static */
bool Samfile::file_exists(const char * filename) {
	struct stat stFileInfo;
	return stat(filename, &stFileInfo) == 0;
}

/* static */
size_t Samfile::aligned_region_size(bam1_t * bam) {
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

void Samfile::print_bam_debug(bam1_t * b, const Mask & read) { // Only for debug
	DEFAULT_CHANNEL << "=============" << endl;
	DEFAULT_CHANNEL << read.sequence << ' ' << read.quality << endl;
	unsigned short int intero;
	string s = string((char *)b->data,b->data_len);
			for (size_t u = 0; u < s.size(); u++)
				if (s[u] < 32 or s[u] > 126)
					s[u] = ' ';
			//DEFAULT_CHANNEL << "==================" << endl;
			//DEFAULT_CHANNEL << str << endl;
			DEFAULT_CHANNEL << "------------------" << endl;
			DEFAULT_CHANNEL << " tid: " << b->core.tid
			     << " pos: " << b->core.pos
			     << " bin: " << b->core.bin
			     << " qual " << b->core.qual
			     << " l-qname: " << b->core.l_qname
			     << " flag: " << b->core.flag
			     << " n_cigar: " << b->core.n_cigar
			     << " l_qseq: " << b->core.l_qseq
			     << " mtid: " << b->core.mtid
			     << " mpos: " << b->core.mpos
			     << " isize: " << b->core.isize << endl
			     << " l_aux: " << b->l_aux
			     << " data_len: " << b->data_len
			     << " m_data: " << b->m_data << endl
			     //<< strlen((char *)b->data) << ":"
			     << s << endl << endl;
	for (int j = 0; j < b->data_len; j+=40) {
		for (int i = j; i < b->data_len and i < j+40 ; i++) {
			if (i == (b->data_len-b->l_aux))
				DEFAULT_CHANNEL << " |";
			DEFAULT_CHANNEL.width(3);
			if ((b->data[i] >= 32) and (b->data[i] <= 126)) {
				DEFAULT_CHANNEL << right << (char)b->data[i];
			} else {
				DEFAULT_CHANNEL << right << ' ';
			}
		}
		DEFAULT_CHANNEL << hex << endl;
		for (int i = j; i < b->data_len and i < j+40 ; i++) {
			if (i == (b->data_len-b->l_aux))
				DEFAULT_CHANNEL << " |";
			DEFAULT_CHANNEL.width(3);
			intero = (uint8_t)b->data[i];
			DEFAULT_CHANNEL << right << intero;
		}
		DEFAULT_CHANNEL << dec << endl << endl;
	}
}

string Samfile::header_description(const Options & options) {
	string sep;
	sep = string("@CO\t");
	sep.append(reading_group_id);
	sep.append(": ");
	stringstream output;
	output << "@HD\tVN:1.0\tSO:unsorted\n";
	output << "@PG\tID:ERNE\tVN:" << PACKAGE_VERSION << "\tCL:";
	for (int i = 0; i < options.argc; i++)
		output << ' ' << options.argv[i];
	output << endl << sep << package_description() << endl
			<< sep << "--query1           = " <<  options.query1 << endl
			<< sep << "--query2           = " <<  options.query2 << endl
			<< sep << "--output           = " <<  options.output_file << endl
			<< sep << "--sample           = " <<  options.sample << endl
			<< sep << "--reference        = " <<  options.reference_file << endl
			<< sep << "--contamination-reference = " <<  options.contamination_file << endl
			<< sep << "--auto-errors      = " <<  (options.auto_errors == true? "true" : "false") << endl
			<< sep << "--errors-rate      = " << options.errors_rate << endl
			<< sep << "--errors           = " << options.common_errors_allowed << endl
			//			<< sep << "--errors-delta     = " << errors_delta << endl
			<< sep << "--threads          = " << options.threads_number << endl
			//			<< sep << "--use-bases        = " << use_bases_flag << endl
			//			<< sep << "--trim-left        = " << trim_left << endl
			//			<< sep << "--trim-right       = " << trim_right << endl
			//<< sep << "--vectors-file     = " << options.vectors_file << endl
			<< sep << "--min-size         = " << options.min_size << endl
			<< sep << "--min-phred-value-CLC = " << options.min_phred_value_MOTT << endl
			<< sep << "--min-mean-phread-quality = " << options.min_mean_quality << endl
			//			<< sep << "--no-auto-trim     = " << (not auto_trim_trigger == true? "true" : "false") << endl
			<< sep << "--no-quality-check = " << (not options.quality_check == true? "true" : "false") << endl
			//			<< sep << "--fastest          = " << (fastest == true? "true" : "false") << endl
			;
	output << "@RG\tID:" << reading_group_id << "\tSM:" << options.sample << endl;
	/*
	for(int i=0; i< H.globaltolocal.Contigs; i++) {
		output << "@SQ\tSN:"<< H.globaltolocal.contigsNames[i] <<"\tLN:"<< (H.globaltolocal.endPositions[i] - H.globaltolocal.startPositions[i] +1) << '\n';
	}
	 */
	return output.str();
}

void Samfile::print_output_paired(Mask &first, Mask &second, print_method_t print_method) {
	bam1_t bam1;
	bam1_t bam2;
	bam1.data = new uint8_t[BAM_DATA_DEFAULT_LENGTH];
	bam2.data = new uint8_t[BAM_DATA_DEFAULT_LENGTH];

	/*
	bool proper_pair = false;
	if (not first.contaminated and not second.contaminated and (first.algn > 0) and (second.algn > 0)
			and (first.contig == second.contig) and (
			  (first.strand and not second.strand and (first.position + first.get_good_length()) < second.position )
			  or
			  (second.strand and not first.strand and (second.position + second.get_good_length()) < first.position) ) )
		proper_pair = true;
	else
		proper_pair = false;
	*/

	// compute not common suffix
	string first_id;
	string second_id;
	{
		const char * a = first.id.c_str();
		const char * b = second.id.c_str();
		size_t i = 0;
		while (a[i] != '\0' and b[i] != '\0' and a[i] == b[i])
			i++;
		if (i > 0 and a[i-1] == '/') // ILLUMINA standard
			i--;
		first_id = first.id.substr(0,i);
		second_id = second.id.substr(0,i);
	}

	if (print_method == BOTH or print_method == FIRST_ONLY) {
		// compute first read
		if (not first.contaminated and first.algn > 0) {
			bam1.core.flag = BAM_FPAIRED | BAM_FREAD1;
			if (not first.strand)
				bam1.core.flag |= BAM_FREVERSE;
			if (not first.primary)
				bam1.core.flag |= BAM_FSECONDARY;

			char * quality_converted = first.quality_conversion();

			Mask::t_CIGAR cigar = first.get_CIGAR();
			string compact;
			first.compact_DNA(compact);

			bam1.core.pos = first.position-1;
			bam1.core.tid = first.contig;
			bam1.core.bin = bam_reg2bin(bam1.core.pos,bam1.core.pos+first.get_good_length()); // serve?
			if (second.algn == 0 or second.contaminated) {
				bam1.core.flag |= BAM_FMUNMAP;
				bam1.core.mpos = first.position-1;
				bam1.core.mtid = first.contig;
			} else {
				bam1.core.mpos = second.position-1;
				bam1.core.mtid = second.contig;
				if (not second.strand)
					bam1.core.flag |= BAM_FMREVERSE;
			}
			if (first.algn == 1)
				bam1.core.qual = 60;
			else
				bam1.core.qual = 0;
			bam1.core.l_qname = first_id.size()+1;
			bam1.core.n_cigar = cigar.second;
			bam1.core.l_qseq = first.sequence.size();
			if (not second.contaminated and second.algn > 0 and first.contig == second.contig) {
				bam1.core.isize = (second.get_five_prime()-first.get_five_prime());
				if ( (first.strand and not second.strand and (first.position + first.get_good_length()) < second.position )
						or
						(second.strand and not first.strand and (second.position + second.get_good_length()) < first.position) ) {
					if ((not insert_size_check) or (insert_size_check and (abs(bam1.core.isize) >= insert_size_min) and (abs(bam1.core.isize) <= insert_size_max)
							and ( (first.strand and not second.strand and (first.position + first.get_good_length()) < second.position )
									or
									(second.strand and not first.strand and (second.position + second.get_good_length()) < first.position)
							)))
						bam1.core.flag |= BAM_FPROPER_PAIR;
				}
			} else
				bam1.core.isize = 0;
			bam1.m_data = BAM_DATA_DEFAULT_LENGTH;

			// all variable-length data, concatenated;
			stringstream data;
			stringstream aux;
			data << first_id << '\0' << cigar.first << compact;
			data.write(quality_converted,bam1.core.l_qseq);

			aux << "RGZ" << reading_group_id << '\0' << "PGZERNE" << '\0' << "NM" << bam_integer(first.NM+first.NM_gap)
						<< "NH" << bam_integer(first.algn) << "IH" << bam_integer(first.IH) << "HI" << bam_integer(first.HI); //"MDZ"	<< read.get_MD(H.TEXT) <<'\0';
			if(first.DELTA.size()>0) {
				aux << "XDZ";
				for(unsigned int i = 0; i< first.DELTA.size(); i++) {
					aux <<  first.DELTA.at(i).first << ":" << first.DELTA.at(i).second << ",";
				}
				aux << '\0';
			}

			bam1.l_aux = aux.str().size();
			data << aux.str();
			bam1.data_len = data.str().size();

			memcpy(bam1.data, data.str().c_str(), bam1.data_len);
			delete [] quality_converted;
		} else {
			bam1.core.flag = BAM_FPAIRED | BAM_FREAD1 | BAM_FUNMAP;
			/*
		t_alignment type = first.get_type();
		if (type == quality_discarded or type == low_complexity)
			bam1.core.flag |= BAM_FQCFAIL;
			 */
			if (first.low_quality or first.low_complexity)
				bam1.core.flag |= BAM_FQCFAIL;

			char * quality_converted = first.quality_conversion();
			string compact;
			first.compact_DNA(compact);

			bam1.core.bin = 0; // serve?
			if (second.algn == 0 or second.contaminated) {
				bam1.core.flag |= BAM_FMUNMAP;
				bam1.core.pos = -1;
				bam1.core.tid = -1;
				bam1.core.mpos = -1;
				bam1.core.mtid = -1;
			} else {
				bam1.core.pos = second.position-1;
				bam1.core.tid = second.contig;
				bam1.core.mpos = second.position-1;
				bam1.core.mtid = second.contig;
				if (not second.strand)
					bam1.core.flag |= BAM_FMREVERSE;
			}
			bam1.core.qual = 0;
			bam1.core.l_qname = first_id.size()+1;
			bam1.core.n_cigar = 0;
			bam1.core.l_qseq = first.sequence.size();
			bam1.core.isize = 0;
			bam1.m_data = BAM_DATA_DEFAULT_LENGTH;

			// all variable-length data, concatenated;
			stringstream data;
			stringstream aux;
			data << first_id << '\0' << compact;
			data.write(quality_converted,bam1.core.l_qseq);
			aux << "RGZ" << reading_group_id << '\0' << "PGZERNE"<<'\0'<<"NHC"<<'\0'<<"IHC"<<'\0';
			if (first.contaminated)
				aux << "XCZ" << CR->globalToLocal.contigsNames[first.contig] << '\0';
			bam1.l_aux = aux.str().size();
			data << aux.str();
			bam1.data_len = data.str().size();

			memcpy(bam1.data, data.str().c_str(), bam1.data_len);
			delete [] quality_converted;
		}
	}

	if (print_method == BOTH or print_method == SECOND_ONLY) {
		// now compute second read
		if (not second.contaminated and second.algn > 0) {
			bam2.core.flag = BAM_FPAIRED | BAM_FREAD2;
			if (not second.strand)
				bam2.core.flag |= BAM_FREVERSE;
			if (not first.primary)
				bam2.core.flag |= BAM_FSECONDARY;

			char * quality_converted = second.quality_conversion();

			Mask::t_CIGAR cigar = second.get_CIGAR();
			string compact;
			second.compact_DNA(compact);

			bam2.core.pos = second.position-1;
			bam2.core.tid = second.contig;
			bam2.core.bin = bam_reg2bin(bam2.core.pos,bam2.core.pos+second.get_good_length()); // serve?
			if (first.algn == 0 or first.contaminated) {
				bam2.core.flag |= BAM_FMUNMAP;
				bam2.core.mpos = second.position-1;
				bam2.core.mtid = second.contig;
			} else {
				bam2.core.mpos = first.position-1;
				bam2.core.mtid = first.contig;
				if (not first.strand)
					bam2.core.flag |= BAM_FMREVERSE;
			}
			if (second.algn == 1)
				bam2.core.qual = 60;
			else
				bam2.core.qual = 0;
			bam2.core.l_qname = second_id.size()+1;
			bam2.core.n_cigar = cigar.second;
			bam2.core.l_qseq = second.sequence.size();
			if (not first.contaminated and first.algn > 0 and second.contig == first.contig) {
				bam2.core.isize = (first.get_five_prime()-second.get_five_prime());
				if ( (first.strand and not second.strand and (first.position + first.get_good_length()) < second.position )
						or
						(second.strand and not first.strand and (second.position + second.get_good_length()) < first.position) ) {
					if ((not insert_size_check) or (insert_size_check and (abs(bam2.core.isize) >= insert_size_min) and (abs(bam2.core.isize) <= insert_size_max)
							and ( (first.strand and not second.strand and (first.position + first.get_good_length()) < second.position )
									or
									(second.strand and not first.strand and (second.position + second.get_good_length()) < first.position)
							)))
						bam2.core.flag |= BAM_FPROPER_PAIR;
				}

			} else
				bam2.core.isize = 0;

			bam2.m_data = BAM_DATA_DEFAULT_LENGTH;

			// all variable-length data, concatenated;
			stringstream data;
			stringstream aux;
			data << second_id << '\0' << cigar.first << compact;
			data.write(quality_converted,bam2.core.l_qseq);

			aux << "RGZ" << reading_group_id << '\0' << "PGZERNE" << '\0' << "NM" << bam_integer(second.NM+second.NM_gap)
						<< "NH" << bam_integer(second.algn) << "IH" << bam_integer(second.IH) << "HI" << bam_integer(second.HI); //"MDZ"	<< read.get_MD(H.TEXT) <<'\0';
			if(second.DELTA.size()>0) {
				aux << "XDZ";
				for(unsigned int i = 0; i< second.DELTA.size(); i++) {
					aux << second.DELTA.at(i).first << ":" << second.DELTA.at(i).second << ",";
				}
				aux << '\0';
			}

			bam2.l_aux = aux.str().size();
			data << aux.str();
			bam2.data_len = data.str().size();

			memcpy(bam2.data, data.str().c_str(), bam2.data_len);
			delete [] quality_converted;
		} else {
			bam2.core.flag = BAM_FPAIRED | BAM_FREAD2 | BAM_FUNMAP;
			/*
		t_alignment type = second.get_type();
		if (type == quality_discarded or type == low_complexity)
			bam2.core.flag |= BAM_FQCFAIL;
			 */
			if (second.low_quality or second.low_complexity)
				bam2.core.flag |= BAM_FQCFAIL;

			char * quality_converted = second.quality_conversion();
			string compact;
			second.compact_DNA(compact);

			bam2.core.bin = 0; // serve?
			if (first.algn == 0 or first.contaminated) {
				bam2.core.flag |= BAM_FMUNMAP;
				bam2.core.pos = -1;
				bam2.core.tid = -1;
				bam2.core.mpos = -1;
				bam2.core.mtid = -1;
			} else {
				bam2.core.pos = first.position-1;
				bam2.core.tid = first.contig;
				bam2.core.mpos = first.position-1;
				bam2.core.mtid = first.contig;
				if (not first.strand)
					bam2.core.flag |= BAM_FMREVERSE;
			}
			bam2.core.qual = 0;
			bam2.core.l_qname = second_id.size()+1;
			bam2.core.n_cigar = 0;
			bam2.core.l_qseq = second.sequence.size();
			bam2.core.isize = 0;
			bam2.m_data = BAM_DATA_DEFAULT_LENGTH;

			// all variable-length data, concatenated;
			stringstream data;
			stringstream aux;
			data << second_id << '\0' << compact;
			data.write(quality_converted,bam2.core.l_qseq);
			aux << "RGZ" << reading_group_id << '\0' << "PGZERNE"<<'\0'<<"NHC"<<'\0'<<"IHC"<<'\0';
			if (second.contaminated)
				aux << "XCZ" << CR->globalToLocal.contigsNames[second.contig] << '\0';
			bam2.l_aux = aux.str().size();
			data << aux.str();
			bam2.data_len = data.str().size();

			memcpy(bam2.data, data.str().c_str(), bam2.data_len);
			delete [] quality_converted;
		}
	}

	if (print_method != NO_ONE) {
		mutex::scoped_lock lock(write_mutex);
		if (print_method == BOTH or print_method == FIRST_ONLY)
			samwrite(sam_file,&bam1);
		if (print_method == BOTH or print_method == SECOND_ONLY)
			samwrite(sam_file,&bam2);
	}
	if (gui_output) {
		stringstream ss;
		if (first.discarded)
			ss << ">Q\n";
		else if (first.contaminated)
			ss << ">C\n";
		else if (first.algn == 0)
			ss << ">N\n";
		else if (first.algn > 0) {
			ss << ">A\t" << first.algn << '\t' << bam1.core.tid << '\t' <<  bam1.core.pos << '\t' << aligned_region_size(&bam1) << '\t';
			if (first.DELTA.size() > 0) {
				int min = first.DELTA.at(0).first;
				int max = first.DELTA.at(0).first;
				for(unsigned int i = 1; i< first.DELTA.size(); i++) {
					if (first.DELTA.at(i).second > 0 and min > first.DELTA.at(i).first)
						min = first.DELTA.at(i).first;
					if (first.DELTA.at(i).second > 0 and max < first.DELTA.at(i).first)
						max = first.DELTA.at(i).first;
				}
				ss << min << '\t' << max << ((bam1.core.flag & BAM_FPROPER_PAIR) ? '1' : '0') << '\n';
			} else
				ss << (first.NM+first.NM_gap) << '\t' << (first.NM+first.NM_gap)<< ((bam1.core.flag & BAM_FPROPER_PAIR) ? '1' : '0') << '\n';
		} else
			ss << ">E" << first.id << '\n';
		if (second.discarded)
			ss << ">Q\n";
		else if (second.contaminated)
			ss << ">C\n";
		else if (second.algn == 0)
			ss << ">N\n";
		else if (second.algn > 0) {
			ss << ">A\t" << second.algn << '\t' << bam2.core.tid << '\t' <<  bam2.core.pos << '\t' << aligned_region_size(&bam2) << '\t';
			if (second.DELTA.size() > 0) {
				int min = second.DELTA.at(0).first;
				int max = second.DELTA.at(0).first;
				for(unsigned int i = 1; i< second.DELTA.size(); i++) {
					if (second.DELTA.at(i).second > 0 and min > second.DELTA.at(i).first)
						min = second.DELTA.at(i).first;
					if (second.DELTA.at(i).second > 0 and max < second.DELTA.at(i).first)
						max = second.DELTA.at(i).first;
				}
				ss << min << '\t' << max << ((bam1.core.flag & BAM_FPROPER_PAIR) ? '1' : '0') << '\n';
			} else
				ss << (second.NM+second.NM_gap) << '\t' << (second.NM+second.NM_gap) << ((bam1.core.flag & BAM_FPROPER_PAIR) ? '1' : '0') << '\n';

		} else
			ss << ">E" << second.id << '\n';
		{
			mutex::scoped_lock lock(gui_output_lock);
			DEFAULT_CHANNEL << ss.str();
		}
	}
	delete [] bam1.data;
	delete [] bam2.data;
}

void Samfile::print_output(bam1_t* b) {
	samwrite(sam_file,b);
}

void Samfile::print_output(Mask &read) {
	bam1_t bam;
	bam.data = new uint8_t[BAM_DATA_DEFAULT_LENGTH];
	if (not read.contaminated and read.algn > 0) {
		bam.core.flag = 0;
		if (not read.strand)
			bam.core.flag |= BAM_FREVERSE;

		if (not read.primary)
			bam.core.flag |= BAM_FSECONDARY;

		char * quality_converted = read.quality_conversion();


		Mask::t_CIGAR cigar = read.get_CIGAR();
		string compact;
		read.compact_DNA(compact);

		bam.core.pos = read.position-1;
		bam.core.tid = read.contig;
		bam.core.bin = bam_reg2bin(bam.core.pos,bam.core.pos+read.get_good_length()); // serve?
		bam.core.mpos = -1;
		bam.core.mtid = -1;
		if (read.algn == 1)
			bam.core.qual = 60;
		else
			bam.core.qual = 0;
		bam.core.l_qname = read.id.size()+1;
		bam.core.n_cigar = cigar.second;
		bam.core.l_qseq = read.sequence.size();
		bam.core.isize = 0;
		bam.m_data = BAM_DATA_DEFAULT_LENGTH;

		// all variable-length data, concatenated;
		stringstream data;
		stringstream aux;
		data << read.id << '\0' << cigar.first << compact;
		data.write(quality_converted,bam.core.l_qseq);

		aux << "RGZ" << reading_group_id << '\0' << "PGZERNE" << '\0' << "NM" << bam_integer(read.NM+read.NM_gap)
						<< "NH" << bam_integer(read.algn) << "IH" << bam_integer(read.IH) << "HI" << bam_integer(read.HI); //"MDZ"	<< read.get_MD(H.TEXT) <<'\0';

		if(read.DELTA.size()>0) {
			aux << "XDZ";
			string delta;
			for(unsigned int i = 0; i< read.DELTA.size(); i++) {
				aux << read.DELTA.at(i).first << ":" << read.DELTA.at(i).second << ";";
			}
			aux << '\0';
		}

		bam.l_aux = aux.str().size();
		data << aux.str();
		bam.data_len = data.str().size();

		memcpy(bam.data, data.str().c_str(), bam.data_len);
		delete [] quality_converted;
	} else {
		bam.core.flag = BAM_FUNMAP;
		/*
		t_alignment type = read.get_type();
		if (type == quality_discarded or type == low_complexity)
			bam.core.flag |= BAM_FQCFAIL;
		*/
		if (read.low_quality or read.low_complexity)
			bam.core.flag |= BAM_FQCFAIL;

		char * quality_converted = read.quality_conversion();
		string compact;
		read.compact_DNA(compact);

		bam.core.pos = -1;
		bam.core.tid = -1;
		bam.core.bin = 0; // serve?
		bam.core.mpos = -1;
		bam.core.mtid = -1;
		bam.core.qual = 0;
		bam.core.l_qname = read.id.size()+1;
		bam.core.n_cigar = 0;
		bam.core.l_qseq = read.sequence.size();
		bam.core.isize = 0;
		bam.m_data = BAM_DATA_DEFAULT_LENGTH;

		// all variable-length data, concatenated;
		stringstream data;
		stringstream aux;
		data << read.id << '\0' << compact;
		data.write(quality_converted,bam.core.l_qseq);
		aux << "RGZ" << reading_group_id << '\0' << "PGZERNE"<<'\0'<<"NHC"<<'\0'<<"IHC"<<'\0';
		if (read.contaminated)
			aux << "XCZ" << CR->globalToLocal.contigsNames[read.contig] << '\0';
		bam.l_aux = aux.str().size();
		data << aux.str();
		bam.data_len = data.str().size();

		memcpy(bam.data, data.str().c_str(), bam.data_len);
		delete [] quality_converted;
	}

	{
		mutex::scoped_lock lock(write_mutex);
		samwrite(sam_file,&bam);
	}
	if (gui_output) {
		stringstream ss;
		if (read.discarded)
			ss << ">Q\n";
		else if (read.contaminated)
			ss << ">C\n";
		else if (read.algn == 0)
			ss << ">N\n";
		else if (read.algn > 0) {
			ss << ">A\t" << read.algn << '\t' << bam.core.tid << '\t' <<  bam.core.pos << '\t' << aligned_region_size(&bam)
				<< '\t';
			if (read.DELTA.size() > 0) {
				int min = read.DELTA.at(0).first;
				int max = read.DELTA.at(0).first;
				for(unsigned int i = 1; i< read.DELTA.size(); i++) {
					if (read.DELTA.at(i).second > 0 and min > read.DELTA.at(i).first)
						min = read.DELTA.at(i).first;
					if (read.DELTA.at(i).second > 0 and max < read.DELTA.at(i).first)
						max = read.DELTA.at(i).first;
				}
				ss << min << '\t' << max << '\n';
			} else
				ss << (read.NM+read.NM_gap) << '\t' << read.NM << '\n';
		} else
			ss << ">E" << read.id << '\n';
		{
			mutex::scoped_lock lock(gui_output_lock);
			DEFAULT_CHANNEL << ss.str();
		}
	}
	delete [] bam.data;
}
