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

#include "Module_FILTER.h"

#define SEQUENCES_FOR_BLOCK_FILTERING 131072

mutex ffa_write_mutex;
Hash H;
bool auto_errors;

unsigned long int totalReadNum;
unsigned long int LowQualityReads;
unsigned long int FilteredReads;
unsigned long int ContaminatedReads;
unsigned long int totalReadLength;

ofstream outputFileFW;
ofstream outputFileRV;
ofstream outputFileUNPAIR;

Fasta::FASTQ_encoding fastq_input_format;
Fasta::FASTQ_encoding fastq_output_format;

namespace modules {

bool has_reference = false;

void Module_FILTER::execute(const Options & options) {
	has_reference = options.reference_file.length() > 0;

	if (has_reference)
		H.load(options.reference_file.c_str());

	string output = options.output_file;

	auto_errors = true;

	string fw_out = output+"_1.fastq";
	outputFileFW.open(fw_out.c_str());

	string rv_out = output+"_2.fastq";
	string unpaired = output+"_unpaired.fastq";
	if (options.query2.size() > 0) {
		outputFileRV.open(rv_out.c_str());
		outputFileUNPAIR.open(unpaired.c_str());
	}

	// now open input files
	Auto_Unzip * file_fw = new Auto_Unzip(options.query1);
	Auto_Unzip * file_rv = NULL;
	if (options.query2.size() > 0)
		file_rv = new Auto_Unzip(options.query2);

	if (not options.force_fastqformat) {
		fastq_input_format = Fasta::check_FASTQ_type_file(options.query1.c_str());
		if (fastq_input_format == Fasta::unknown_fastq_encoding) {
			ERROR_CHANNEL << "The program can not autodetect input format file, please use --force-illumina or --force-standard options" << endl;
			exit(4);
		}
		VERBOSE_CHANNEL << "FASTQ input format: " << (fastq_input_format == Fasta::illumina_fastq_encoding ? "ILLUMINA" : "STANDARD") << " (autodetected)" << endl;
	} else {
		fastq_input_format = options.fastqformat;
		VERBOSE_CHANNEL << "FASTQ input format: " << (fastq_input_format == Fasta::illumina_fastq_encoding ? "ILLUMINA" : "STANDARD") << " (forced)" << endl;
	}

	fastq_output_format = options.preserve_encoding ? fastq_input_format : Fasta::standard_fastq_encoding;
	VERBOSE_CHANNEL << "FASTQ output format: " << (fastq_output_format == Fasta::illumina_fastq_encoding ? "ILLUMINA" : "STANDARD") << endl;

 	thread_group threads;

 	totalReadNum=0;
 	LowQualityReads=0;
 	FilteredReads=0;
 	ContaminatedReads=0;
 	totalReadLength=0;

	for (int id = 0; id < options.threads_number; id++)
			threads.create_thread(boost::bind(&filter_reads, id, file_fw, file_rv, &options));
	threads.join_all();
	delete file_fw;
	if (file_rv != NULL)
		delete file_rv;

	VERBOSE_CHANNEL << "total number of reads " << totalReadNum << "\n";
	VERBOSE_CHANNEL << "total number of successfully filtered reads " << FilteredReads << "\n";
	VERBOSE_CHANNEL << "total number of discarded reads " << LowQualityReads << "\n";
	VERBOSE_CHANNEL << "total number of contaminated reads " << ContaminatedReads << "\n";
	VERBOSE_CHANNEL << "total length of quality filtered reads " << totalReadLength << "\n";
}

/* static */
int Module_FILTER::read_sequences_conditional(Auto_Unzip * first, Auto_Unzip * second, int num_seq, Mask sequences[], Fasta::FASTQ_encoding format_type, bool gui_output, bool pe) {
	if (pe)
		return read_sequences(*first, *second, num_seq, sequences, format_type, gui_output);
	else
		return read_sequences(*first, num_seq, sequences, format_type, gui_output);
}

/* static */
void Module_FILTER::filter_reads(int id, Auto_Unzip *reads_fw, Auto_Unzip *reads_rv, const Options * options) {


	Mask * sequences = new Mask[SEQUENCES_FOR_BLOCK_FILTERING];
	int read_seq=0;
	while((read_seq = read_sequences_conditional(reads_fw, reads_rv, SEQUENCES_FOR_BLOCK_FILTERING, sequences, fastq_input_format,options->gui_output, (options->query2.size() > 0))) > 0) {
		//DEFAULT_CHANNEL << "thread " << id << " processing " << read_seq << " reads\n";
		t_errors calculated_errors_first,calculated_errors_second;

		if (options->query2.size() > 0 ) {
			for(int i = 0 ; i < read_seq; i+=2) {
				Mask & masked_fw = sequences[i];
				Mask & masked_rv = sequences[i+1];

				masked_fw.quality_trimming_MOTT(options->min_phred_value_MOTT,options->min_mean_quality,options->min_size);
				masked_rv.quality_trimming_MOTT(options->min_phred_value_MOTT,options->min_mean_quality,options->min_size);

				if (auto_errors)
					calculated_errors_first = round((double)masked_fw.get_good_length() / options->errors_rate);
				else
					calculated_errors_first = options->common_errors_allowed;

				if (auto_errors)
					calculated_errors_second = round((double)masked_rv.get_good_length() / options->errors_rate);
				else
					calculated_errors_second = options->common_errors_allowed;

				if (has_reference and not(masked_fw.discarded or masked_fw.low_quality )) {
					//if quality control passed, check contamination
					Items solution_fw;
					H.search(masked_fw.get_good_sequence(), solution_fw, calculated_errors_first);
					masked_fw.algn = solution_fw.size(); // memorize number of solutions, does not matter if more than one occurence has been found
				} else if(! options->trim) { // if --no-auto-trim active force alignment of low quality reads
					Items solution_fw;
					H.search(masked_fw.get_sequence(), solution_fw, 5);
					masked_fw.algn = solution_fw.size(); // memorize number of solutions, does not matter if more than one occurence has been found
				}
				if (has_reference and not(masked_rv.discarded or  masked_rv.low_quality)) {
					//if quality control passed, check contamination
					Items solution_rv;
					H.search(masked_rv.get_good_sequence(), solution_rv, calculated_errors_second);
					masked_rv.algn = solution_rv.size();
				} else if(! options->trim) { // if --no-auto-trim active force alignment of low quality reads
					Items solution_rv;
					H.search(masked_rv.get_sequence(), solution_rv, 5);
					masked_rv.algn = solution_rv.size();
				}
			}
			// now print results
			{
				mutex::scoped_lock lock(ffa_write_mutex);

				if (options->trim  ) {
					for(int i=0; i < read_seq; i+=2){
						Mask & masked_fw = sequences[i];
						Mask & masked_rv = sequences[i+1];
						if( !(masked_fw.discarded || masked_rv.discarded || masked_fw.low_quality || masked_rv.low_quality)) {
							if(masked_fw.algn == 0 && masked_rv.algn == 0 ) {
								Fasta read_fw;
								masked_fw.set_fasta_masked(read_fw);
								read_fw.set_FASTQ_type(fastq_output_format);
								outputFileFW << read_fw;
								Fasta read_rv;
								masked_rv.set_fasta_masked(read_rv);
								read_rv.set_FASTQ_type(fastq_output_format);
								outputFileRV << read_rv;
								FilteredReads+=2;
								totalReadLength += (read_fw.get_sequence().length() + read_rv.get_sequence().length());
							} else {
								ContaminatedReads+=2;
							}
						} else if  (!(masked_fw.discarded ||  masked_fw.low_quality) ) {
							if(masked_fw.algn == 0) {
								Fasta read_fw;
								masked_fw.set_fasta_masked(read_fw);
								read_fw.set_FASTQ_type(fastq_output_format);
								outputFileUNPAIR << read_fw;
								FilteredReads++;
								totalReadLength += (read_fw.get_sequence().length());
							} else {
								ContaminatedReads++;
							}
							LowQualityReads++;
						} else if  (!(masked_rv.discarded ||  masked_rv.low_quality) ) {
							if(masked_rv.algn == 0) {
								Fasta read_rv;
								masked_rv.set_fasta_masked(read_rv);
								read_rv.set_FASTQ_type(fastq_output_format);
								outputFileUNPAIR << read_rv;
								FilteredReads++;
								totalReadLength += (read_rv.get_sequence().length());
							} else {
								ContaminatedReads++;
							}
							LowQualityReads++;
						} else {
							LowQualityReads+=2;
						}
					}

					totalReadNum +=read_seq;
				} else {
					for(int i=0; i < read_seq; i+=2){
						Mask & masked_fw = sequences[i];
						Mask & masked_rv = sequences[i+1];
						if(masked_fw.algn == 0 && masked_rv.algn == 0 ) {
							Fasta read_fw;
							masked_fw.set_original_fasta(read_fw);
							read_fw.set_FASTQ_type(fastq_output_format);
							outputFileFW << read_fw;

							Fasta read_rv;
							masked_rv.set_original_fasta(read_rv);
							read_rv.set_FASTQ_type(fastq_output_format);
							outputFileRV << read_rv;
							FilteredReads+=2;
							totalReadLength += (read_fw.get_sequence().length() + read_rv.get_sequence().length());
						} else {
							ContaminatedReads+=2;
						}
					}
					totalReadNum += read_seq;
				}
			}
		} else {
			for(int i = 0 ; i < read_seq; i++) {
				Mask & masked_fw = sequences[i];

				masked_fw.quality_trimming_MOTT(options->min_phred_value_MOTT,options->min_mean_quality,options->min_size);

				if (auto_errors)
					calculated_errors_first = round((double)masked_fw.get_good_length() / options->errors_rate);
				else
					calculated_errors_first = options->common_errors_allowed;

				if (has_reference and not(masked_fw.discarded or masked_fw.low_quality )) {
					//if quality control passed, check contamination
					Items solution_fw;
					H.search(masked_fw.get_good_sequence(), solution_fw, calculated_errors_first);
					masked_fw.algn = solution_fw.size(); // memorize number of solutions, does not matter if more than one occurence has been found
				}
			}
			// now print results
			{
				mutex::scoped_lock lock(ffa_write_mutex);

				for(int i=0; i < read_seq; i++){
					Mask & masked_fw = sequences[i];
					Mask & masked_rv = sequences[i+1];
					if( !(masked_fw.discarded || masked_fw.low_quality )) {
						if(masked_fw.algn == 0 && masked_rv.algn == 0 ) {
							Fasta read_fw;
							masked_fw.set_fasta_masked(read_fw);
							read_fw.set_FASTQ_type(fastq_output_format);
							outputFileFW << read_fw;
							FilteredReads++;
							totalReadLength += (read_fw.get_sequence().length());
						} else {
							ContaminatedReads++;
						}
					} else {
						LowQualityReads++;
					}
				}

				totalReadNum +=read_seq;
			}
		}
		delete [] sequences;
		sequences = new Mask[SEQUENCES_FOR_BLOCK_FILTERING];
		read_seq = 0;
	}
}

}
