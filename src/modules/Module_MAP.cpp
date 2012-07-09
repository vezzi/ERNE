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

#include "Module_MAP.h"

namespace modules {

RNGType Module_MAP::rng;
boost::uniform_real<> Module_MAP::double_range(0,1);
boost::variate_generator<RNGType&, boost::uniform_real<> > Module_MAP::guessRes(Module_MAP::rng, Module_MAP::double_range);

void Module_MAP::initialize_parameters(const Options & options) {

	rng.seed(time(0));

	force_fastqformat = options.force_fastqformat;
	fastqformat = options.fastqformat;
	threads_number = options.threads_number;
	gui_output = options.gui_output;

	quality_check = options.quality_check;
	trim = options.trim;
	min_phred_value_MOTT = options.min_phred_value_MOTT;
	min_mean_quality = options.min_mean_quality;
	min_size = options.min_size;
	auto_errors = options.auto_errors;
	common_errors_allowed = options.common_errors_allowed;
	errors_rate = options.errors_rate;
	contamination_check = options.contamination_check;
	printAll = options.printAll;
	toBePrinted = options.toBePrinted;
	gap = options.gap;
	transcriptome = options.transcriptome;
	seed_sizes = options.seed_sizes;
	seed_errors = options.seed_errors;
	max_gap = options.max_gap;

	query1 = options.query1;
	query2 = options.query2;

	delta = options.delta;
	indels = options.indels;
	indels_max_value = options.indels_max_value;
	insert_size_check = options.insert_size_check;
	insert_size_min = options.insert_size_min;
	insert_size_max = options.insert_size_max;

	if(options.program_mode == Options::program_bs5) // TODO: rimuovere
		methyl_reconstruction = true;
	else
		methyl_reconstruction = false;

}


void Module_MAP::execute(const Options & options) {
	clock_t started = clock();
	processed = 0;

	initialize_parameters(options);

	if (options.contamination_check) {
		VERBOSE_CHANNEL << "Reading contamination SA in file " << options.contamination_file << endl;
		clock_t start = clock();
		CR.load(options.contamination_file.c_str());
		CR.set_indels_max_value(indels_max_value);
		clock_t end = clock();

		VERBOSE_CHANNEL << "Hash (contamination) loading time = ";
		print_formatted_time(VERBOSE_CHANNEL, (end - start) / (double)CLOCKS_PER_SEC);
		VERBOSE_CHANNEL << endl;
	}

	VERBOSE_CHANNEL << "Loading Hash Table" << endl;

	{
		clock_t start = clock();

		H.load(options.reference_file.c_str());
		H.set_indels_max_value(indels_max_value);
		clock_t end = clock();

		VERBOSE_CHANNEL << "Hash (reference) loading time = ";
		print_formatted_time(VERBOSE_CHANNEL, (end - start) / (double)CLOCKS_PER_SEC);
		VERBOSE_CHANNEL << endl;
	}

	if(options.program_mode == Options::program_bs5 and not H.methyl_hash){

		ERROR_CHANNEL << "error: this hash is for standard search (erne-map), and cannot be used for bisulfite alignment." << endl;
		ERROR_CHANNEL << "Try \"--help\" for help" << endl;
		exit(1);

	} else 	if(options.program_mode == Options::program_map and H.methyl_hash){

		ERROR_CHANNEL << "error: this hash is for bisulfite search (erne-bs5), and cannot be used standard alignment." << endl;
		ERROR_CHANNEL << "Try \"--help\" for help" << endl;
		exit(1);
	}


	output_samfile.open_file(options,H,CR);

	// PROCESS PAIR ENDS FILES
	if (options.paired_ends) {
		VERBOSE_CHANNEL << "Processing files " << options.query1 << " and " << options.query2 << endl;
		process_pair_reads();
	} else {			// PROCESS SINGLE READ FILES
		VERBOSE_CHANNEL << "Processing file " << options.query1 << endl;
		process_single_reads();
	}

	output_samfile.close_file();

	VERBOSE_CHANNEL << "Finished!" << endl;

	clock_t finished = clock();
	double dif = (finished - started) / (double)CLOCKS_PER_SEC;

	VERBOSE_CHANNEL << "Time spent: " << dif << "s (";
	print_formatted_time(VERBOSE_CHANNEL, dif);
	VERBOSE_CHANNEL << ')' << endl;
	VERBOSE_CHANNEL << "Queries per second: " << (processed/dif)	<<endl;


}

// read from the read file then align SEQUENCES_FOR_BLOCK queries and then print them
void Module_MAP::process_single_reads() {
	if (not force_fastqformat) {
		fastqformat = Fasta::check_FASTQ_type_file(query1);
		if (fastqformat == Fasta::unknown_fastq_encoding) {
			ERROR_CHANNEL << "The program can not autodetect input format file, please use --force-illumina or --force-standard options" << endl;
			exit(4);
		}
		VERBOSE_CHANNEL << "FASTQ format: " << (fastqformat == Fasta::illumina_fastq_encoding ? "ILLUMINA" : "STANDARD") << " (autodetected)" << endl;
	} else
		VERBOSE_CHANNEL << "FASTQ format: " << (fastqformat == Fasta::illumina_fastq_encoding ? "ILLUMINA" : "STANDARD") << " (forced)" << endl;

	Auto_Unzip input_file(query1);

	thread_group threads;
	for (int id=0; id < threads_number; id++)
		threads.create_thread(boost::bind(&this->process_single_reads_thr, this, id, &input_file));
	threads.join_all();
}

void Module_MAP::process_single_reads_thr(Module_MAP * search, int id, Auto_Unzip * input) {
	Mask * sequences = new Mask[SEQUENCES_FOR_BLOCK];
	int read_seq=0;
	unsigned int deltaZero=0;
	vector<Mask> printable_solutions;
	vector<pair<int, int > > DeltaOutput;
	while ((read_seq = read_sequences(*input, SEQUENCES_FOR_BLOCK, sequences, search->fastqformat,search->gui_output)) != 0) {
		Items solutions;
		ItemsGapped solutions_gapped;
		t_errors calculated_errors;

		for(int i=0; i< read_seq; i++){
			Mask & s = sequences[i];
			if (search->trim) {
				//	check_routine(sequences[i], 0);
				s.quality_trimming_MOTT(search->min_phred_value_MOTT,search->min_mean_quality,search->min_size);

				if (s.status_discarded()) {
					if (s.status_low_complexity())
						s.low_complexity = true; //s.set_type(low_complexity);
					else
						s.low_quality = true; //s.set_type(quality_discarded);
					printable_solutions.push_back(s);
					continue;
				}
			}
			if (search->auto_errors)
				calculated_errors = round((double)s.get_good_length() / search->errors_rate);
			else
				calculated_errors = search->common_errors_allowed;

			{
				t_errors count = 0;
				for (t_pattern_length i = s.get_good_region_start()-1; (i < s.get_good_region_stop()) and (count <= calculated_errors); i++)
					if (s.sequence[i] == 'N' or s.sequence[i] == 'n')
						count++;
				if (count > calculated_errors) {
					printable_solutions.push_back(s);
					continue;
				}
			}

			/** ALIGNMENT ERNE **/
			solutions.clear();
			DeltaOutput.clear();

			bool contaminated = false;
			if (search->contamination_check)
				search->CR.search(s.get_good_sequence(),solutions,calculated_errors, deltaZero , DeltaOutput, false);

			if (solutions.size() > 0) {
				contaminated = true;
			} else {
				search->H.search(s.get_good_sequence(),solutions,calculated_errors, search->delta, DeltaOutput, search->indels);
			}

			if (solutions.size() == 0) {
				/** Try gapped **/
				solutions_gapped.clear();
				if (search->gap)
					search->H.search_gapped(s.get_good_sequence(),solutions_gapped,search->seed_sizes,search->seed_errors,calculated_errors,search->max_gap);
				if (solutions_gapped.size() == 0) {// size 0 means no alignment found
					//s.set_type(alignments_not_found);
					printable_solutions.push_back(s);
					continue;
				} else {
					if (not search->printAll) {
						int alignments = solutions_gapped.size();
						const ResultItemGapped & HM = solutions_gapped.at((int)(guessRes() * alignments));
						s.gapped = true;
						s.globalPosition = HM.globalPosition1;
						s.algn = alignments;
						s.HI =1;
						s.IH =1;
						s.primary = true;
						s.strand = HM.strand;
						s.NM = HM.errors1;
						s.NM_gap = HM.errors2;
						s.contig = HM.contig;
						unsigned int offset = search->H.globalToLocal.startPositions[s.contig];
						s.position = HM.globalPosition1 - offset + 1 ;
						s.position_gap = HM.globalPosition2 - offset + 1 ;
						s.length1_gap = HM.length1;
						s.length2_gap = HM.length2;
						s.contaminated = contaminated;
						if (not contaminated)
							s.adjust_transcriptome_coordinates(search->H.TEXT,search->H.textLength);
						printable_solutions.push_back(s);

					} else { // printALL
						unsigned int alignments = (search->toBePrinted < solutions_gapped.size()) ? search->toBePrinted : solutions_gapped.size() ;
						for (unsigned int processed = 0; processed < alignments; processed++) {
							const ResultItemGapped & HM = solutions_gapped.at(processed);
							s.gapped = true;
							s.globalPosition = HM.globalPosition1;
							s.algn = alignments;
							s.HI =1;
							s.IH =1;
							(processed == 0 ) ? s.primary = true : s.primary = false;
							s.strand = HM.strand;
							s.NM = HM.errors1;
							s.NM_gap = HM.errors2;
							s.contig = HM.contig;
							unsigned int offset = search->H.globalToLocal.startPositions[s.contig];
							s.position = HM.globalPosition1 - offset + 1 ;
							s.position_gap = HM.globalPosition2 - offset + 1 ;
							s.length1_gap = HM.length1;
							s.length2_gap = HM.length2;
							s.contaminated = contaminated;
							if (not contaminated)
								s.adjust_transcriptome_coordinates(search->H.TEXT,search->H.textLength);
							printable_solutions.push_back(s);
							processed++;
						}
					}
				}
			} else if (not search->printAll) {
				sort(solutions.begin(), solutions.end(), ResultItem::less()); // sort solutions
				solutions.erase(unique(solutions.begin(), solutions.end(), ResultItem::equal()), solutions.end());
				int alignments = solutions.size();
				const ResultItem & HM = solutions.at((int)(guessRes() * alignments));
				s.HI =1;
				s.IH =1;
				s.primary = true;
				s.globalPosition = HM.globalPosition;
				s.algn = alignments;
				s.strand = HM.strand;
				s.NM = HM.errors;
				s.indels = HM.indels;
				if (HM.indels) {
					if(HM.strand) {
						s.cigarVector = search->H.alignSW(s.get_good_sequence().c_str(), s.get_good_length(), HM.globalPosition);
					} else {
						s.cigarVector = search->H.alignSW(reverse_complement_standalone_str(s.get_good_sequence()).c_str(), s.get_good_length(), HM.globalPosition);
					}
				}
				if (contaminated) {
					s.contig = search->CR.globalToLocal.searchContig(HM.globalPosition); // find the contig/scaffold
					s.position = HM.globalPosition - search->CR.globalToLocal.startPositions[s.contig] + 1;
				} else {
					s.contig = search->H.globalToLocal.searchContig(HM.globalPosition); // find the contig/scaffold
					s.position = HM.globalPosition - search->H.globalToLocal.startPositions[s.contig] + 1;
				}
				s.contaminated = contaminated;
				s.DELTA = DeltaOutput; //.insert(s.DELTA.begin(), DeltaOutput.begin(), DeltaOutput.end());//
				printable_solutions.push_back(s);
			} else { // printAll
				// memorize all printable solutions
				sort(solutions.begin(), solutions.end(), ResultItem::less()); // sort solutions
				solutions.erase(unique(solutions.begin(), solutions.end(), ResultItem::equal()), solutions.end());
				unsigned int alignments = (search->toBePrinted < solutions.size()) ? search->toBePrinted : solutions.size();
				for (unsigned int processed = 0; processed < alignments; processed++) {
					// while I print enough solution or there are no more solutions
					const ResultItem & HM = solutions.at(processed);
					s.globalPosition = HM.globalPosition;
					s.algn = solutions.size();
					s.indels = HM.indels;
					s.IH = alignments;
					s.HI = processed+1;
					(processed == 0 ) ? s.primary = true : s.primary = false;
					s.strand = HM.strand;
					s.NM = HM.errors;
					if (HM.indels) {
						if(HM.strand) {
							s.cigarVector = search->H.alignSW(s.get_good_sequence().c_str(), s.get_good_length(), HM.globalPosition);
						} else {
							s.cigarVector = search->H.alignSW(reverse_complement_standalone_str(s.get_good_sequence()).c_str(), s.get_good_length(), HM.globalPosition);
						}
					}

					if (contaminated) {
						s.contig = search->CR.globalToLocal.searchContig(HM.globalPosition); // find the contig/scaffold
						s.position = HM.globalPosition - search->CR.globalToLocal.startPositions[s.contig] + 1;
					} else {
						s.contig = search->H.globalToLocal.searchContig(HM.globalPosition); // find the contig/scaffold
						s.position = HM.globalPosition - search->H.globalToLocal.startPositions[s.contig] + 1;
					}
					s.contaminated = contaminated;
					s.DELTA = DeltaOutput;
					printable_solutions.push_back(s);

				}
			}
		}
		// now print all
		for(unsigned int i=0; i < printable_solutions.size(); i++)
			search->output_samfile.print_output(printable_solutions.at(i));
		search->processed += read_seq; // TODO mettere un lock
		delete [] sequences;
		sequences = new Mask[SEQUENCES_FOR_BLOCK];
		printable_solutions.clear();
		read_seq = 0;
	}
	delete [] sequences;
}


void Module_MAP::process_pair_reads() {
	if (not force_fastqformat) {
		fastqformat = Fasta::check_FASTQ_type_file(query1);
		VERBOSE_CHANNEL << "FASTQ format: " << (fastqformat == Fasta::illumina_fastq_encoding ? "ILLUMINA" : "STANDARD") << " (autodetected)" << endl;
	}

	Auto_Unzip input_file1(query1);
	Auto_Unzip input_file2(query2);

	thread_group threads;

	for (int id=0; id < threads_number; id++)
		threads.create_thread(boost::bind(&process_pair_reads_thr, this, id, &input_file1, &input_file2));
	threads.join_all();
}

void Module_MAP::process_pair_reads_thr(Module_MAP * search, int id, Auto_Unzip * inputFile1, Auto_Unzip * inputFile2) {
	Mask * sequences = new Mask[SEQUENCES_FOR_BLOCK];
	int read_seq=0;
	unsigned int deltaZero=0;

	while((read_seq = read_sequences(*inputFile1, *inputFile2, SEQUENCES_FOR_BLOCK, sequences, search->fastqformat, search->gui_output)) > 0) {
		Items solutionFirst;
		Items solutionSecond;
		ItemsGapped gapped_solutionFirst;
		ItemsGapped gapped_solutionSecond;
		vector<pair<int, int > > deltaOutputFirst;
		vector<pair<int, int > > deltaOutputSecond;
		t_errors calculated_errors_first,calculated_errors_second;

		for(int i = 0 ; i < read_seq; i += 2) {
			solutionFirst.clear();
			solutionSecond.clear();
			Mask & first = sequences[i];
			Mask & second = sequences[i+1];

			if(search->methyl_reconstruction)  // TODO: remove!
				second.reverse_complement();

			bool first_allowed;
			bool second_allowed;
			//if --auto-trim option is set perform auto trimming
			if (search->trim) {
				first.quality_trimming_MOTT(search->min_phred_value_MOTT,search->min_mean_quality,search->min_size);
				second.quality_trimming_MOTT(search->min_phred_value_MOTT,search->min_mean_quality,search->min_size);
			}
			//compute number of allowed errors
			if (search->auto_errors)
				calculated_errors_first = round((double)first.get_good_length() / search->errors_rate);
			else
				calculated_errors_first = search->common_errors_allowed;
			t_errors count = 0;
			for (t_pattern_length i = first.get_good_region_start()-1; (i < first.get_good_region_stop()) and (count <= calculated_errors_first); i++)
				if ((first.sequence[i] == 'N') or (first.sequence[i] == 'n'))
					count++;
			if (count > calculated_errors_first)
				first_allowed = false;
			else
				first_allowed = not second.status_discarded();

			//preprocessing of first read finished
			//now compute allowed errors for second read
			if (search->auto_errors)
				calculated_errors_second = round((double)second.get_good_length() / search->errors_rate);
			else
				calculated_errors_second = search->common_errors_allowed;
			count = 0;
			for (t_pattern_length i = second.get_good_region_start()-1; (i < second.get_good_region_stop()) and (count <= calculated_errors_second); i++)
				if (second.sequence[i] == 'N' or second.sequence[i] == 'n')
					count++;
			if (count > calculated_errors_second)
				second_allowed = false;
			else
				second_allowed = not second.status_discarded();

			//and set it to allowed or not allowed
			//check if reads are contaminated if a contamination sequence is provided
			bool first_contaminated = false;
			bool second_contaminated = false;


			if (not first_allowed) {
				if (first.status_low_complexity())
					first.low_complexity = true; //first.set_type(low_complexity);
				else
					first.low_quality = true; //first.set_type(quality_discarded);
			} else {
				solutionFirst.clear();
				deltaOutputFirst.clear();
				if (search->contamination_check)
					search->CR.search(first.get_good_sequence(), solutionFirst, calculated_errors_first,deltaZero ,deltaOutputFirst, false );


				if (solutionFirst.size() > 0)
					first_contaminated = true;
				else
					search->H.search(first.get_good_sequence(),solutionFirst,calculated_errors_first, search->delta, deltaOutputFirst, search->indels);

				if (solutionFirst.size() == 0) {
					// Try gapped
					gapped_solutionFirst.clear();
					if (search->gap)
						search->H.search_gapped(first.get_good_sequence(),gapped_solutionFirst,search->seed_sizes,search->seed_errors,calculated_errors_first,search->max_gap);
					sort(gapped_solutionFirst.begin(), gapped_solutionFirst.end(), ResultItemGapped::less() ); // sort solutions
					gapped_solutionFirst.erase(unique(gapped_solutionFirst.begin(), gapped_solutionFirst.end(), ResultItemGapped::equal()), gapped_solutionFirst.end());
				} else {
					sort(solutionFirst.begin(), solutionFirst.end(), ResultItem::less() ); // sort solutions
					solutionFirst.erase(unique(solutionFirst.begin(), solutionFirst.end(), ResultItem::equal()), solutionFirst.end());
				}
			} // first read searched

			if (not second_allowed) {
				if (second.status_low_complexity())
					second.low_complexity = true; //second.set_type(low_complexity);
				else
					second.low_quality = true; //second.set_type(quality_discarded);
			} else {
				solutionSecond.clear();
				deltaOutputSecond.clear();

				if (search->contamination_check)
					search->CR.search(second.get_good_sequence(),solutionSecond, calculated_errors_second, deltaZero, deltaOutputFirst, false);

				if (solutionSecond.size() > 0)
					second_contaminated = true;
				else
					search->H.search(second.get_good_sequence(),solutionSecond,calculated_errors_second, search->delta, deltaOutputSecond, search->indels);

				if (solutionSecond.size() == 0) {
					/** Try gapped **/
					solutionSecond.clear();
					if (search->gap)
						search->H.search_gapped(second.get_good_sequence(),gapped_solutionSecond,search->seed_sizes,search->seed_errors,calculated_errors_second,search->max_gap);
					sort(gapped_solutionSecond.begin(), gapped_solutionSecond.end(), ResultItemGapped::less() ); // sort solutions
					gapped_solutionSecond.erase(unique(gapped_solutionSecond.begin(), gapped_solutionSecond.end(), ResultItemGapped::equal()), gapped_solutionSecond.end());

				} else {
					sort(solutionSecond.begin(), solutionSecond.end(), ResultItem::less() ); // sort solutions
					solutionSecond.erase(unique(solutionSecond.begin(), solutionSecond.end(), ResultItem::equal()), solutionSecond.end());
				}
			} //second read searched

			GenericItem * first_solution = NULL;
			GenericItem * second_solution = NULL;

			GenericItems all_solutions_first;
			GenericItems all_solutions_second;

			for (Items::iterator iter = solutionFirst.begin(); iter != solutionFirst.end(); iter++)
				all_solutions_first.push_back(*iter);
			for (ItemsGapped::iterator iter = gapped_solutionFirst.begin(); iter != gapped_solutionFirst.end(); iter++)
				all_solutions_first.push_back(*iter);
			for (Items::iterator iter = solutionSecond.begin(); iter != solutionSecond.end(); iter++)
				all_solutions_second.push_back(*iter);
			for (ItemsGapped::iterator iter = gapped_solutionSecond.begin(); iter != gapped_solutionSecond.end(); iter++)
				all_solutions_second.push_back(*iter);
			if (solutionFirst.size() > 0 and gapped_solutionFirst.size() > 0)
				sort(all_solutions_first.begin(), all_solutions_first.end(), GenericItem::less() );
			if (solutionSecond.size() > 0 and gapped_solutionSecond.size() > 0)
				sort(all_solutions_second.begin(), all_solutions_second.end(), GenericItem::less() );

			Hash & HCR_first = first_contaminated ? search->CR : search->H;
			Hash & HCR_second = second_contaminated ? search->CR : search->H;

			int sizeFirst = all_solutions_first.size();
			int sizeSecond = all_solutions_second.size();

			size_t chosen_first = 0;
			size_t chosen_second = 0;

			// if first allowed second allowed and they align against the same sequence
			if (first_allowed and second_allowed and (first_contaminated == second_contaminated)) {
				int p = 0;
				if (sizeFirst > 1 and sizeSecond == 1) {
					int q = sizeFirst-1;
					second_solution = &all_solutions_second[0];
					if (all_solutions_first.at(q).min_position() <= second_solution->min_position()) {
						chosen_first = q; // tutte le read sono più piccole => scegli l'ultima
					} else if (all_solutions_first.at(0).min_position() >= second_solution->min_position()){
						chosen_first = 0; // tutte le read sono più grandi => scegli la prima
					} else {
						int m;
						while (p != q) {
							m = p+(int)(q-p)/2+1;
							if (all_solutions_first.at(m).min_position() >= second_solution->min_position()){
								q = m-1;
							} else {
								p = m;
							}
						}
						unsigned int d1 = second_solution->min_position() - all_solutions_first.at(p).min_position();
						unsigned int d2 = all_solutions_first.at(p+1).min_position() - second_solution->min_position();
						chosen_first = (d1 <= d2) ? p : p+1;
						first_solution = (d1 <= d2) ? &all_solutions_first.at(p) : &all_solutions_first.at(p+1);

					}
					first_solution = &all_solutions_first.at(chosen_first);
				} else if (sizeFirst == 1 and sizeSecond > 1){
					int q = sizeSecond-1;
					first_solution = &all_solutions_first[0];
					if (all_solutions_second.at(q).min_position() <= first_solution->min_position()) {
						chosen_second = q; // tutte le read sono più piccole => scegli l'ultima
					} else if (all_solutions_second.at(0).min_position() >= first_solution->min_position()){
						chosen_second = 0; // tutte le read sono più grandi => scegli la prima
					} else {
						int m;
						while (p != q) {
							m = p+(int)(q-p)/2+1;
							if (all_solutions_second.at(m).min_position() >= first_solution->min_position()){
								q = m-1;
							} else {
								p = m;
							}
						}
						unsigned int d1 = first_solution->min_position() - all_solutions_second.at(p).min_position();
						unsigned int d2 = all_solutions_second.at(p+1).min_position() - first_solution->min_position();
						chosen_second = (d1 <= d2) ? p : p+1;
					}
					second_solution = &all_solutions_second.at(chosen_second);
				}
			}

			if (sizeFirst != 0) { // size 0 means no alignment found
				//save all possible solutions here
				if (first_solution == NULL)
					first_solution = &all_solutions_first.at((int)(guessRes() * all_solutions_first.size()));

				first.set_globalPosition(first_solution->min_position());
				first.indels = first_solution->indels;
				if (first_solution->indels and not first_solution->gapped) {
					if(first_solution->strand)
						first.cigarVector = search->H.alignSW(first.get_good_sequence().c_str(), first.get_good_length(), first_solution->globalPosition1);
					else
						first.cigarVector = search->H.alignSW(reverse_complement_standalone_str(first.get_good_sequence()).c_str(), first.get_good_length(), first_solution->globalPosition1);
				}

				first.set_algn(all_solutions_first.size());
				first.set_strand(first_solution->strand); //strand
				first.set_NM(first_solution->errors1);
				first.contaminated = first_contaminated;
				first.IH = 1;
				first.HI = 1;
				size_t contig = HCR_first.globalToLocal.searchContig(first_solution->globalPosition1); // find the contig/scaffold
				first.set_contig(contig); // find the name
				first.set_position(first_solution->globalPosition1 - HCR_first.globalToLocal.startPositions[contig] + 1);
				if (first_solution->gapped) {
					first.gapped = true;
					first.set_position_gap(first_solution->globalPosition2 - HCR_first.globalToLocal.startPositions[contig] + 1);
					first.length1_gap = first_solution->length1;
					first.length2_gap = first_solution->length2;
					first.set_NM_gap(first_solution->errors2);
				}
				first.DELTA = deltaOutputFirst;
			}

			if (sizeSecond != 0) { // size 0 means no alignment found
				//save all possible solutions here
				if (second_solution == NULL)
					second_solution = &all_solutions_second.at((int)(guessRes() * all_solutions_second.size()));


				second.set_globalPosition(second_solution->globalPosition1);
				second.indels = second_solution->indels;
				if (second_solution->indels and not second_solution->gapped) {
					if(second_solution->strand)
						second.cigarVector = search->H.alignSW(second.get_good_sequence().c_str(), second.get_good_length(), second_solution->globalPosition1);
					else
						second.cigarVector = search->H.alignSW(reverse_complement_standalone_str(second.get_good_sequence()).c_str(), second.get_good_length(), second_solution->globalPosition1);
				}
				second.set_algn(all_solutions_second.size());
				second.set_strand(second_solution->strand); //strand
				second.set_NM(second_solution->errors1);
				second.contaminated = second_contaminated;
				second.IH = 1;
				second.HI = 1;
				size_t contig = HCR_second.globalToLocal.searchContig(second_solution->globalPosition1); // find the contig/scaffold
				second.set_contig(contig); // find the name
				second.set_position(second_solution->globalPosition1 - HCR_second.globalToLocal.startPositions[contig] + 1);
				if (second_solution->gapped) {
					second.gapped = true;
					second.set_position_gap(second_solution->globalPosition2 - HCR_second.globalToLocal.startPositions[contig] + 1);
					second.length1_gap = second_solution->length1;
					second.length2_gap = second_solution->length2;
					second.set_NM_gap(second_solution->errors2);
				}
				second.DELTA = deltaOutputSecond;
			}

			if (search->printAll) {
				unsigned int alignments = (search->toBePrinted < all_solutions_first.size()) ? search->toBePrinted : all_solutions_first.size();
				for (unsigned int processed = 0; processed < alignments; processed++) {
					// while I print enough solution or there are no more solutions
					const GenericItem & HM = all_solutions_first.at(processed);
					Mask first_temp = first;
					first_temp.globalPosition = HM.globalPosition1;
					first_temp.algn = sizeFirst;
					first_temp.indels = HM.indels;
					first_temp.IH = alignments;
					first_temp.HI = processed+1;
					(processed == chosen_first ) ? first_temp.primary = true : first_temp.primary = false;
					first_temp.strand = HM.strand;
					first_temp.NM = HM.errors1;
					if (HM.indels) {
						if(HM.strand) {
							first_temp.cigarVector = search->H.alignSW(first_temp.get_good_sequence().c_str(), first_temp.get_good_length(), HM.globalPosition1);
						} else {
							first_temp.cigarVector = search->H.alignSW(reverse_complement_standalone_str(first_temp.get_good_sequence()).c_str(), first_temp.get_good_length(), HM.globalPosition1);
						}
					}

					if (first_contaminated) {
						first_temp.contig = search->CR.globalToLocal.searchContig(HM.globalPosition1); // find the contig/scaffold
						first_temp.position = HM.globalPosition1 - search->CR.globalToLocal.startPositions[first_temp.contig] + 1;
					} else {
						first_temp.contig = search->H.globalToLocal.searchContig(HM.globalPosition1); // find the contig/scaffold
						first_temp.position = HM.globalPosition1 - search->H.globalToLocal.startPositions[first_temp.contig] + 1;
					}
					first_temp.contaminated = first_contaminated;
					first_temp.DELTA = deltaOutputFirst;
					search->output_samfile.print_output_paired(first_temp,second,Samfile::FIRST_ONLY);

				}
				if (alignments == 0)
					search->output_samfile.print_output_paired(first,second,Samfile::FIRST_ONLY);


				alignments = (search->toBePrinted < all_solutions_second.size()) ? search->toBePrinted : all_solutions_second.size();
				for (unsigned int processed = 0; processed < alignments; processed++) {
					// while I print enough solution or there are no more solutions
					const GenericItem & HM = all_solutions_second.at(processed);
					Mask second_temp = second;
					second_temp.globalPosition = HM.globalPosition1;
					second_temp.algn = sizeSecond;
					second_temp.indels = HM.indels;
					second_temp.IH = alignments;
					second_temp.HI = processed+1;
					(processed == chosen_second ) ? second_temp.primary = true : second_temp.primary = false;
					second_temp.strand = HM.strand;
					second_temp.NM = HM.errors1;
					if (HM.indels) {
						if(HM.strand) {
							second_temp.cigarVector = search->H.alignSW(second_temp.get_good_sequence().c_str(), second_temp.get_good_length(), HM.globalPosition1);
						} else {
							second_temp.cigarVector = search->H.alignSW(reverse_complement_standalone_str(second_temp.get_good_sequence()).c_str(), second_temp.get_good_length(), HM.globalPosition1);
						}
					}

					if (second_contaminated) {
						second_temp.contig = search->CR.globalToLocal.searchContig(HM.globalPosition1); // find the contig/scaffold
						second_temp.position = HM.globalPosition1 - search->CR.globalToLocal.startPositions[second_temp.contig] + 1;
					} else {
						second_temp.contig = search->H.globalToLocal.searchContig(HM.globalPosition1); // find the contig/scaffold
						second_temp.position = HM.globalPosition1 - search->H.globalToLocal.startPositions[second_temp.contig] + 1;
					}
					second_temp.contaminated = second_contaminated;
					second_temp.DELTA = deltaOutputSecond;
					search->output_samfile.print_output_paired(first,second_temp,Samfile::SECOND_ONLY);

				}
				if (alignments == 0)
					search->output_samfile.print_output_paired(first,second,Samfile::SECOND_ONLY);
			}
		}

		// now print all (only for not print-all)
		if (not search->printAll) {
			for(int i=0; i<read_seq; i += 2)
				search->output_samfile.print_output_paired(sequences[i],sequences[i+1]);
		}
		search->processed += read_seq;

		delete [] sequences;
		sequences = new Mask[SEQUENCES_FOR_BLOCK];
		read_seq = 0;
	}
	delete [] sequences;
}

void Module_MAP::free_hash_memory(){
	H.free_hash_memory();
	CR.free_hash_memory();
}


}
