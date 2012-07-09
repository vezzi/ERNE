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

#include "Module_BS5.h"

#include <stdlib.h>
#include <stdio.h>

#include <samtools/bam.h>
#include <samtools/sam.h>

#include "Module_MAP.h"
#include "options/Options_BS5.h"

namespace modules {

Module_BS5::Module_BS5(){
}

void Module_BS5::execute(const Options & options) {

	clock_t start = clock();

	error_threshold = options.error_threshold;
	C_cov = options.max_C_cov;//min number of C that must be covered to consider the read as covered (Starting value)

	string alignments_path_current_cycle("");

	string methyl_annotation_path(options.output_file.c_str());
	methyl_annotation_path.append("_methyl_annotations.txt");

	//---------------------------------------1. ERNE alignment---------------------------

	VERBOSE_CHANNEL << endl << "*** STEP 1: ERNE ALIGNMENT *** " << endl <<endl;

	Options_BS5 temp_options(options);

	Module_MAP search;

	if(options.skip_alignment)
		search.H.load(temp_options.reference_file.c_str());
	else{
		if (not options.alignment_only) {
			if(temp_options.bam_format)
				temp_options.output_file.append("_0.bam");
			else
				temp_options.output_file.append("_0.sam");
		}

		search.execute(temp_options);
	}
	search.free_hash_memory(); //maintain only H->TEXT

	if(options.alignment_only)
		exit(0);

	H = &search.H;
	CR = &search.CR;

	methylConsensus_T = new unsigned char[H->textLength];
	methylConsensus_C = new unsigned char[H->textLength];
	mult_coverage = new unsigned char [H->textLength];

	for(unsigned int i=0;i<H->textLength;i++){
		methylConsensus_T[i]=0;
		methylConsensus_C[i]=0;
		mult_coverage[i]=0;
	}

	//--------------------------------------2. cycles-----------------------------

	VERBOSE_CHANNEL << endl << "*** STEP 2: METHYLATION RECONSTRUCTION CYCLES ***" << endl<<endl;

	unsigned int alignments_size;
	unsigned int uniquely_mapped_reads = 0;
	unsigned int unique_reads_cycle_0=0;//unique reads in the initial bam
	unsigned int number_of_reads = 0;
	unsigned int multiple_reads;//multiple reads remaining at each cycle
	unsigned int no_alignments=0;//number of reads with no alignments found

	totCoverage=0;

	bam1_t** ali1 = NULL;//alignments for the read 1
	bam1_t** ali2 = NULL;//alignments for the read 2

	if(!options.skip_alignment){//if skip_alignment == true, then all the alignments are used, and this file is not created since is already present.
		samfile_output.set_reading_group_id(search.get_reading_group_id());
		samfile_output.open_file(options,*H,*CR);
	}

	cycle = 0;

	VERBOSE_CHANNEL << "C_cov = " << options.max_C_cov << endl;

	while(C_cov>=options.min_C_cov && !(cycle==1 && options.skip_alignment)  && !(cycle==1 && options.only_unique)){

		alignments_size = 0;
		uniquely_aligned = 0;
		multiple_reads = 0; //number of reads that can be mapped in multiple positions (analyzed in the next cycle)

		SAM_istream SAM_istream_current_cycle;//from here bam1_t objects are read
		Samfile samfile_next_cycle;//here bam1_t objects are saved for the next cycle
		samfile_next_cycle.set_reading_group_id(search.get_reading_group_id()); // pass reading_group_id to output file

		if(options.skip_alignment){
			alignments_path_current_cycle.clear();
			alignments_path_current_cycle.append(options.output_file);
			SAM_istream_current_cycle.open(alignments_path_current_cycle.c_str());//use out.bam / out.sam if skip_alignment==true

		}else{

			std::ostringstream t;
			if(options.bam_format)
				t  << "_" << cycle << ".bam";
			else
				t  << "_" << cycle << ".sam";

			samfile_next_cycle.open_file(options,*H,*CR,cycle+1);
			alignments_path_current_cycle.clear();
			alignments_path_current_cycle.append(options.output_file);
			alignments_path_current_cycle.append(t.str().data());
			SAM_istream_current_cycle.open(alignments_path_current_cycle.c_str());//use out.bam_cycle / out.sam_cycle if skip_alignment==false

		}

		VERBOSE_CHANNEL << "cycle number " <<cycle<< " ... ";

		if(options.paired_ends){

			bam1_t *b1=NULL;
			//b1 = bam_init1();
			bam1_t *b2=NULL;
			//b2 = bam_init1();
			int IH1,IH2,res1,res2;

			while(!SAM_istream_current_cycle.eof()) {

				b1 = SAM_istream_current_cycle.read();

				//-------------read1--------------------
				IH1 = BAM1::IH(b1);

				if(IH1>0){
					ali1 = new bam1_t*[IH1];
					ali1[0] = b1;

					for(int i=1;i<IH1;i++)//save in ali1 all the alignments for read1
						ali1[i] = SAM_istream_current_cycle.read();

				}

				//--------------read2---------------

				b2 = SAM_istream_current_cycle.read();

				IH2 = BAM1::IH(b2);

				if(IH2>0){
					ali2 = new bam1_t*[IH2];
					ali2[0] = b2;

					for(int i=1;i<IH2;i++)//save in ali1 all the alignments for read2
						ali2[i] = SAM_istream_current_cycle.read();

				}

				if(IH1>0)
					alignments_size ++;

				if(IH2>0)
					alignments_size ++;

				//----------now in ali1 and ali2 we have all the alignments for read1 and read2--------------

				res1=0;
				res2=0;

				if(IH1<=0&&cycle==0)
					no_alignments++;
				if(IH2<=0&&cycle==0)
					no_alignments++;

				if(cycle>0){//in the first cycle only uniquely-mapping reads are used

					if(IH1!=1)
						res1 = check_alignments(ali1, IH1);
					if(IH2!=1)
						res2 = check_alignments(ali2, IH2);

				}else{

					if(IH1>1)
						res1 = MULTIPLE;
					if(IH2>1)
						res2 = MULTIPLE;
					if(IH1<=0)
						res1 = NO_ALIGNMENTS;
					if(IH2<=0)
						res2 = NO_ALIGNMENTS;

				}

				if(res1 == MULTIPLE && res2 == MULTIPLE)
					multiple_reads+=2;
				else if((res1 == NO_ALIGNMENTS && res2 == MULTIPLE)||(res2 == NO_ALIGNMENTS && res1 == MULTIPLE))
					multiple_reads++;

				if( ((res1 == MULTIPLE && res2 == MULTIPLE) || (res1 == NO_ALIGNMENTS && res2 == MULTIPLE) ||
					(res2 == NO_ALIGNMENTS && res1 == MULTIPLE)) && !options.skip_alignment){//one or both reads have to be re-processed: save all the alignments in the temp file

					if(IH1<=0)
						samfile_next_cycle.print_output(b1);//if IH1==0 b1 must be saved to maintain synchronization between read1 and read2!
					else
						for(int i=0;i<IH1;i++)
							samfile_next_cycle.print_output(ali1[i]);
					if(IH2<=0)
						samfile_next_cycle.print_output(b2);
					else
						for(int i=0;i<IH2;i++)
							samfile_next_cycle.print_output(ali2[i]);

				}else if(res1>=0 && res2==MULTIPLE){

					align_and_update_profile(ali1+res1, 1, ali2, IH2);

					if(cycle==0)
						unique_reads_cycle_0+=2;

				}else if(res2>=0 && res1==MULTIPLE){

					align_and_update_profile(ali1, IH1, ali2+res2, 1);

					if(cycle==0)
						unique_reads_cycle_0+=2;

				}else if(res1>=0 && res2==NO_ALIGNMENTS ){//one read can be aligned but the other has no good alignments: align only one read

					align_and_update_profile(ali1[res1]);

					if(cycle==0)
						unique_reads_cycle_0++;

				}else if(res2>=0 && res1==NO_ALIGNMENTS ){//one read can be aligned but the other has no good alignments: align only one read

					align_and_update_profile(ali2[res2]);

					if(cycle==0)
						unique_reads_cycle_0++;

				}else if(res1>=0 && res2>=0 && 	((ali1[res1]->core.flag >> 4)&1) == ((ali2[res2]->core.flag >> 4)&1) ){//same strand

					align_and_update_profile(ali1[res1]);
					align_and_update_profile(ali2[res2]);

					if(cycle==0)
						unique_reads_cycle_0+=2;

				}

				bam_destroy1(b1);
				for(int i=1;i<IH1;i++)
					bam_destroy1(ali1[i]);
				if(IH1>0)
					delete [] ali1;

				bam_destroy1(b2);
				for(int i=1;i<IH2;i++)
					bam_destroy1(ali2[i]);
				if(IH2>0)
					delete [] ali2;

				//b1 = bam_init1();
				//b2 = bam_init1();

				if(cycle==0)
					number_of_reads+=2;

			}

			//bam_destroy1(b1);
			//bam_destroy1(b2);

		}else{//not paired ends

			bam1_t *b1;
			//b1 = bam_init1();
			int IH1,res1;

			while(!SAM_istream_current_cycle.eof()) {

				b1 = SAM_istream_current_cycle.read();

				if(!options.skip_alignment){

					IH1 = BAM1::IH(b1);

					if(IH1>0){
						ali1 = new bam1_t*[IH1];
						ali1[0] = b1;

						for(int i=1;i<IH1;i++)//save in ali1 all the alignments for the read
							ali1[i] = SAM_istream_current_cycle.read();

					}
				}else{//align b1 on the reference

					IH1 = 1;
					ali1 = new bam1_t*[IH1];
					ali1[0] = b1;

				}

				if(IH1>0)
					alignments_size ++;
				//---------now in ali1 we have all the alignments for read1--------------------

				res1=0;

				if(IH1<=0 && cycle==0)
					no_alignments++;

				if(IH1==1 && cycle==0)
					unique_reads_cycle_0++;

				if(cycle>0){

					if(IH1!=1)
						res1 = check_alignments(ali1, IH1);

				}else{

					if(IH1 > 1)
						res1 = MULTIPLE;
					if(IH1 <= 0)
						res1 = NO_ALIGNMENTS;

				}

				if(res1 == MULTIPLE)
					multiple_reads++;

				if(res1 >= 0)
					align_and_update_profile(ali1[res1]);
				else if(res1==MULTIPLE && !options.skip_alignment){
					for(int i=0;i<IH1;i++)
						samfile_next_cycle.print_output(ali1[i]);
				}

				bam_destroy1(b1);
				for(int i=1;i<IH1;i++)
					bam_destroy1(ali1[i]);
				if(IH1>0)
					delete [] ali1;

				//b1 = bam_init1();

				if(cycle==0)
					number_of_reads++;

			}//while

			//bam_destroy1(b1);
		}//else

		uniquely_mapped_reads += uniquely_aligned;
		VERBOSE_CHANNEL << uniquely_aligned << " reads aligned uniquely";

		if(uniquely_aligned ==0){
			C_cov--;
			VERBOSE_CHANNEL << " ... setting C_cov = " << C_cov ;
		}
		VERBOSE_CHANNEL << endl;

		if(!options.skip_alignment && C_cov>=options.min_C_cov && !options.only_unique)
			VERBOSE_CHANNEL << multiple_reads << " multiple reads will be analyzed during next cycle" << endl << endl;

		cycle++;

		SAM_istream_current_cycle.close();//close current cycle file

		if(!options.skip_alignment){

			remove(alignments_path_current_cycle.c_str()); //remove cycle-th temporary bam
			samfile_next_cycle.close_file();//close (cycle+1)-th temporary bam

		}

	}//while C_cov>=min_C_cov

	std::ostringstream t;

	if(options.bam_format)
		t << "_" << cycle << ".bam";
	else
		t << "_" << cycle << ".sam";

	if(!options.skip_alignment){

		remove( string(options.output_file).append(t.str().data()).c_str() );//remove last temp file
		samfile_output.close_file();

	}


	//----------------------filter C with coverage < coverage_threshold

	for(unsigned int i=0;i<H->textLength;i++){

		if(covg(i)<options.coverage_threshold){
			methylConsensus_T[i] = 0;
			methylConsensus_C[i] = 0;
		}

	}

	//--------------------------------NOW CALCULATE METHYLATION LEVELS AND COVERAGE----------------------

	VERBOSE_CHANNEL << endl << "*** STEP 3: FILE OUTPUT AND STATISTICS ***" << endl << endl;

	unsigned int covC_CG=0;//number of covered  in CG context
	unsigned int covC_CHG=0;//number of covered  in CHG context
	unsigned int covC_CHH=0;//number of covered  in CHH context

	unsigned int totC=0;//total number of C in the genome

	unsigned int metCG=0;//number of methylated C in CG context
	unsigned int metCHG=0;//number of methylated C in CHG context
	unsigned int metCHH=0;//number of methylated C in CHH context

	for(unsigned int i=0;i<H->textLength;i++){

		if(H->TEXT[i]=='C' || H->TEXT[i]=='c'||H->TEXT[i]=='G' || H->TEXT[i]=='g') totC++;

		if(i < H->textLength-2){

			if(H->TEXT[i]=='C' || H->TEXT[i]=='c'){//FW (+ strand)

				if(H->TEXT[i+1]=='G' || H->TEXT[i+1]=='g'){//CG
					if(covg(i)>0){
						covC_CG++;
						if (methyl_level(i) > 0.5)
							metCG++;
					}
				}
				else if(H->TEXT[i+2]=='G' || H->TEXT[i+2]=='g'){//CHG
					if(covg(i)>0){
						covC_CHG++;
						if (methyl_level(i) > 0.5)
							metCHG++;
					}
				}
				else {
					if(covg(i)>0){//CHH
						covC_CHH++;
						if (methyl_level(i) > 0.5)
							metCHH++;
					}
				}

			}

		}

		if(i > 2){

			if(H->TEXT[i]=='G' || H->TEXT[i]=='g'){//RC (- strand)

				if(H->TEXT[i-1]=='C'||H->TEXT[i-1]=='c'){//CG
					if(covg(i)>0){
						covC_CG++;
						if(methyl_level(i) > 0.5)
							metCG++;
					}
				}
				else if(H->TEXT[i-2]=='C'||H->TEXT[i-2]=='c'){//CHG
					if(covg(i)>0){
						covC_CHG++;
						if(methyl_level(i) > 0.5)
							metCHG++;
					}
				}
				else{
					if(covg(i)>0){//CHH
						covC_CHH++;
						if(methyl_level(i) > 0.5)
							metCHH++;
					}
				}

			}
		}

	}//for all shifts in H->TEXT

	VERBOSE_CHANNEL << "Number of reads analyzed: " << number_of_reads << endl;
	VERBOSE_CHANNEL << "Number of uniquely mapped reads: " << uniquely_mapped_reads << endl;
	VERBOSE_CHANNEL << "Number of discarded reads: \n \t no alignments found : "<< no_alignments <<endl;
	VERBOSE_CHANNEL <<"\t too large methylation error: " << (int)number_of_reads - (int)multiple_reads - (int)uniquely_mapped_reads -(int)no_alignments<< endl;
	VERBOSE_CHANNEL <<"\t total: " << (int)number_of_reads - (int)multiple_reads - (int)uniquely_mapped_reads << endl;

	VERBOSE_CHANNEL << "Number of multiple reads which cannot be mapped uniquely: " << multiple_reads << endl;

	if(!options.skip_alignment && !options.only_unique)
		VERBOSE_CHANNEL << "Number of multiple reads aligned uniquely on the methylation consensus: " << uniquely_mapped_reads - unique_reads_cycle_0 << " (" << ((double)(uniquely_mapped_reads - unique_reads_cycle_0)/(double)uniquely_mapped_reads)*100 << "% of all uniquely mapped reads)" <<endl;

	VERBOSE_CHANNEL << "Mapping efficiency: " << ((double)uniquely_mapped_reads/(double)number_of_reads)*100<< "%" << endl;

	VERBOSE_CHANNEL << endl << "covered C/number of C: " << (covC_CG+covC_CHG+covC_CHH) << "/" << totC << " (" << ((double)(covC_CG+covC_CHG+covC_CHH)/(double)totC)*100 <<"%)" << endl;
	VERBOSE_CHANNEL << "Coverage average per base: "<< (double)totCoverage/(double)H->textLength << endl<<endl;

	VERBOSE_CHANNEL << "Covered C in CG context: " << covC_CG << " (" << ((double)(covC_CG)/(double)(covC_CG+covC_CHG+covC_CHH))*100 << "% of all covered cytosines)"<<endl;
	VERBOSE_CHANNEL << "Covered C in CHG context: " << covC_CHG << " (" << ((double)(covC_CHG)/(double)(covC_CG+covC_CHG+covC_CHH))*100 << "% of all covered cytosines)"<<endl;
	VERBOSE_CHANNEL << "Covered C in CHH context: " << covC_CHH << " (" << ((double)(covC_CHH)/(double)(covC_CG+covC_CHG+covC_CHH))*100 << "% of all covered cytosines)"<<endl;

	VERBOSE_CHANNEL << endl << "METHYLATION STATISTICS:" << endl << endl;

	VERBOSE_CHANNEL << "methylated C in CG context: " << metCG << " (" << (((double)metCG/(double)(covC_CG))*100) << "% CG cytosines are methylated)" << endl;
	VERBOSE_CHANNEL << "methylated C in CHG context: " << metCHG << " (" << (((double)metCHG/(double)(covC_CHG))*100) << "% CHG cytosines are methylated)" << endl;
	VERBOSE_CHANNEL << "methylated C in CHH context: " << metCHH << " (" << (((double)metCHH/(double)(covC_CHH))*100) << "% CHH cytosines are methylated)" << endl<<endl;

	VERBOSE_CHANNEL << "methylated CG fraction (100 * #methylated in CG context/#methylated): " << 100*(double)metCG/(double)(metCG+metCHG+metCHH) << "%" << endl;
	VERBOSE_CHANNEL << "methylated CHG fraction (100 * #methylated in CHG context/#methylated): " << 100*(double)metCHG/(double)(metCG+metCHG+metCHH) << "%" << endl;
	VERBOSE_CHANNEL << "methylated CHH fraction (100 * #methylated in CHH context/#methylated): " << 100*(double)metCHH/(double)(metCG+metCHG+metCHH) << "%" << endl << endl;

	//-------------------------SAVE METHYL ANNOTATIONS----------------

	FILE* methyl_annotations = fopen(methyl_annotation_path.c_str(),"w");

	const char* context;

	unsigned long int l;

	for(int i=0;i<H->globalToLocal.contigs;i++) {

		l=0;
		for(unsigned int j=H->globalToLocal.startPositions[i];j<=H->globalToLocal.endPositions[i];j++){

			l++;

			if(H->TEXT[j]=='C' || H->TEXT[j]=='c'){//+ strand

				if(j<H->textLength-2){

					if(H->TEXT[j+1]=='G' || H->TEXT[j+1]=='g'){
						context = CG;
					}else if(H->TEXT[j+2]=='G' || H->TEXT[j+2]=='g'){
						context = CHG;
					}else{
						context = CHH;
					}

				}else{
					context = CHH;
				}

				std::ostringstream t;
				t << H->globalToLocal.contigsNames[i].c_str() << "\t" << l << "\t" << l << "\t" << context << ":" << covg(j) << "\t" << methyl_level(j) << "\t+\t" << (unsigned int)mult_coverage[j] << endl;
				fputs(t.str().data(),methyl_annotations);

			}

			if(H->TEXT[j]=='G' || H->TEXT[j]=='g'){//- strand

				if(j>1){

					if(H->TEXT[j-1]=='C' || H->TEXT[j-1]=='c'){
						context = CG;
					}else if(H->TEXT[j-2]=='C' || H->TEXT[j-2]=='c'){
						context = CHG;
					}else{
						context = CHH;
					}

				}else{
					context = CHH;
				}

				std::ostringstream t;
				t << H->globalToLocal.contigsNames[i].c_str() << "\t" <<  l << "\t" << l << "\t" << context << ":" << covg(j) << "\t" << methyl_level(j) << "\t-\t" << (unsigned int)mult_coverage[j] << endl;
				fputs(t.str().data(),methyl_annotations);

			}

		}
	}

	fclose(methyl_annotations);

	//---------------------------------VALIDATION:----------------------
	if(!options.validation_files.empty())
		validation(options);

	clock_t end = clock();
	VERBOSE_CHANNEL << "total time (search + methylation reconstruction) = ";
	print_formatted_time(VERBOSE_CHANNEL, (end - start) / (double)CLOCKS_PER_SEC);
	VERBOSE_CHANNEL << endl;

	delete [] methylConsensus_T;
	delete [] methylConsensus_C;
	delete [] mult_coverage;
}//execute


int Module_BS5::context_in_text(unsigned long i){

	if(i < H->textLength-2){

		if(H->TEXT[i]=='C' || H->TEXT[i]=='c'){//FW (+ strand)

			if(H->TEXT[i+1]=='G' || H->TEXT[i+1]=='g'){//CG
				return cg;
			}
			else if(H->TEXT[i+2]=='G' || H->TEXT[i+2]=='g'){//CHG
				return chg;
			}
			else{
				return chh;
			}
		}

	}

	if(i > 2){

		if(H->TEXT[i]=='G' || H->TEXT[i]=='g'){	//C on (-) strand

			if(H->TEXT[i-1]=='C' || H->TEXT[i-1]=='c'){
				return cg;
			}
			else if(H->TEXT[i-2]=='C'||H->TEXT[i-2]=='c'){//CHG
				return chg;
			}
			else{
				return chh;
			}
		}

	}

	return chh;

}//context

Module_BS5::score *Module_BS5::methyl_distance(bam1_t *b){

	int m;//read length
	unsigned int shift=0;
	bool reverse;

	double md = 0; //methyl distance = average of the methyl distances for every covered C
	unsigned int coverage = 0; //=average of the coverage for every covered C
	unsigned int number_of_covered_C = 0;
	unsigned int number_of_C = 0;

	reverse = ((b->core.flag >> 4)&1) == 1;

	int readLength = b->core.l_qseq;
	char *read = new char[readLength];
	char *s   = (char*)bam1_seq(b);

	for(m=0;m<(b->core.l_qseq);m++) {
		int v = bam1_seqi(s,m);
		read[m] = bam_nt16_rev_table[v];
	}

	if(((b->core.flag >> 2)&1) == 0){
		shift = global_coord(b->core.tid,b->core.pos); //shift in H->TEXT
	} else return NULL; //no alignment

	uint32_t *cigar = bam1_cigar(b);

	unsigned int j=0;//position in the read

	for (unsigned int k = 0; k < b->core.n_cigar; ++k) {

		int op = cigar[k] & BAM_CIGAR_MASK;//operation type
		int l = cigar[k] >> BAM_CIGAR_SHIFT;//number of bases

		switch(op){

			case BAM_CMATCH : 	for(unsigned int i=j; i<j+l; i++){

									if(reverse && H->TEXT[shift] == 'G'){
										md += alpha(read[i],shift);
										coverage += covg(shift);
										if(covg(shift)>0) number_of_covered_C++;
										number_of_C++;
									}
									if(!reverse && H->TEXT[shift] == 'C'){
										md += alpha(read[i],shift);
										coverage += covg(shift);
										if(covg(shift)>0) number_of_covered_C++;
										number_of_C++;
									}

									shift++;

								};j+=l;break;

			case BAM_CINS : j+=l;break;//insertion in the reference
			case BAM_CDEL : shift+=l;break;//deletions from the reference
			case BAM_CREF_SKIP : shift+=l;break;//skip on the reference
			case BAM_CSOFT_CLIP: j+=l;break;//clip on the read with clipped sequence present in read
		    case BAM_CHARD_CLIP: break;//clip on the read with clipped sequence trimmed off
		    case BAM_CPAD: shift+=l; break;

		}

	}//for all cigar pairs

	delete [] read;

	if (number_of_C == 0) return new Module_BS5::score(1000,0,0);//in this case the alignment has no C:discard
	return new Module_BS5::score(md,(double)coverage/(double)number_of_covered_C,number_of_covered_C);

}//methyl_distance(bam1_t b)

unsigned int abs(long a){
	return (a<0?-a:a);
}

void Module_BS5::align_and_update_profile(bam1_t ** ali1, int IH1, bam1_t** ali2, int IH2){

	if(IH1 == 1){
		if(IH2==0){//there are no alignments for read 2: save the unique alignment for read1
			align_and_update_profile(ali1[0]);
		}
		else{//IH2 >0
			unsigned int min_idx=0;
			for(int i=1;i<IH2;i++)
				if(	(abs((long)ali1[0]->core.pos - (long)ali2[i]->core.pos)<abs((long)ali1[0]->core.pos - (long)ali2[min_idx]->core.pos) &&
					((ali1[0]->core.flag >> 4)&1) == ((ali2[i]->core.flag >> 4)&1)) || //same strand
					(((ali1[0]->core.flag >> 4)&1) != ((ali2[min_idx]->core.flag >> 4)&1)) ) //min_idx is wrong: different strands
						min_idx=i;

			if( ((ali1[0]->core.flag >> 4)&1) == ((ali2[min_idx]->core.flag >> 4)&1) ){
				align_and_update_profile(ali1[0]);
				align_and_update_profile(ali2[min_idx]);
			}
		}
	}

	else if(IH2 == 1){
		if(IH1==0){//there are no alignments for read 1: save the unique alignment for read2
			align_and_update_profile(ali2[0]);
		}
		else{//IH1 >0
			unsigned int min_idx=0;
			for(int i=1;i<IH1;i++)
				if(	(abs((long)ali2[0]->core.pos - (long)ali1[i]->core.pos)<abs((long)ali2[0]->core.pos - (long)ali1[min_idx]->core.pos) &&
					((ali2[0]->core.flag >> 4)&1) == ((ali1[i]->core.flag >> 4)&1)) ||
					(((ali2[0]->core.flag >> 4)&1) != ((ali1[min_idx]->core.flag >> 4)&1)) )
						min_idx=i;

			if( ((ali2[0]->core.flag >> 4)&1) == ((ali1[min_idx]->core.flag >> 4)&1) ){
				align_and_update_profile(ali1[min_idx]);
				align_and_update_profile(ali2[0]);
			}
		}
	}

}

int Module_BS5::check_alignments(bam1_t **ali1, int IH){

	int result = MULTIPLE;

	if(IH>0){

		Module_BS5::score **scores = new Module_BS5::score*[IH];
		for(int i=0;i<IH;i++)
			scores[i] = methyl_distance(ali1[i]);

		double covered_ali = 0;

		for(int i=0;i<IH;i++)
			if(scores[i]->C_covered >= C_cov)
				covered_ali++;

		if(covered_ali == IH){

			unsigned int min_idx=0;
			for(int i=1;i<IH;i++)
				if(*scores[i] <= *scores[min_idx] && scores[i]->C_covered >= C_cov)
					min_idx = i;

			if(scores[min_idx]->methyl_distance/(double)(scores[min_idx]->C_covered==0?1:scores[min_idx]->C_covered) <= error_threshold)
				result = min_idx;//best alignment
			else
				result = NO_ALIGNMENTS;

		}

		for(int i=0;i<IH;i++)
			delete scores[i];

		delete [] scores;

	}else
		result = NO_ALIGNMENTS;

	return result;

}

void Module_BS5::align_and_update_profile(bam1_t *b){

	if( &samfile_output != NULL ){//save bam object

		BAM1::modify_tag_int32_t(b,"IH",1);
		BAM1::modify_tag_int32_t(b,"HI",1);

		samfile_output.print_output(b);
	}

	uniquely_aligned++;

	int m;//read length
	unsigned int shift=0;
	bool reverse;

	reverse = ((b->core.flag >> 4)&1) == 1;

	int readLength = b->core.l_qseq;
	char *read = new char[readLength];
	char *s   = (char*)bam1_seq(b);

	for(m=0;m<(b->core.l_qseq);m++) {
		int v = bam1_seqi(s,m);
		read[m] = bam_nt16_rev_table[v];
	}

	if(((b->core.flag >> 2)&1) == 0)
		shift = global_coord(b->core.tid,b->core.pos); //shift in H->TEXT

	uint32_t *cigar = bam1_cigar(b);

	unsigned int j=0;//position in the read

	for (unsigned int k = 0; k < b->core.n_cigar; ++k) {

		int op = cigar[k] & BAM_CIGAR_MASK;//operation type
		int l = cigar[k] >> BAM_CIGAR_SHIFT;//number of bases

		switch(op){

			case BAM_CMATCH : 	for(unsigned int i=j; i<j+l; i++){

									if(reverse){

										if((H->TEXT[shift]=='G'||H->TEXT[shift]=='g') && (read[i]=='A'||read[i]=='a')){

											increment_T(shift);//bisulphite conversion on RC

											if(cycle>0)//this C has been covered using multiple-mapping reads
												increment_mult_coverage(shift);

										}
										else if((H->TEXT[shift]=='G'||H->TEXT[shift]=='g') && (read[i]=='G'||read[i]=='g')){

											increment_C(shift);

											if(cycle>0)
												increment_mult_coverage(shift);

										}

									}else{

										if((H->TEXT[shift]=='C'||H->TEXT[shift]=='c') && (read[i]=='T'||read[i]=='t')){

											increment_T(shift);//bisulphite conversion on FW

											if(cycle>0)
												increment_mult_coverage(shift);

										}
										else if((H->TEXT[shift]=='C'||H->TEXT[shift]=='c') && (read[i]=='C'||read[i]=='c')){

											increment_C(shift);

											if(cycle>0)
												increment_mult_coverage(shift);

										}

									}

									shift++;

									totCoverage++;
								};j+=l;break;

			case BAM_CINS : j+=l;break;//insertion in the reference
			case BAM_CDEL : shift+=l;break;//deletions from the reference
			case BAM_CREF_SKIP : shift+=l;break;//skip on the reference
			case BAM_CSOFT_CLIP: j+=l;break;//clip on the read with clipped sequence present in read
		    case BAM_CHARD_CLIP: break;//clip on the read with clipped sequence trimmed off
		    case BAM_CPAD: shift+=l; break;

		}

	}//for all cigar pairs

	delete [] read;

}//align_and_update_profile(bam1_t b)


void Module_BS5::validation(const Options & options){

	const char* methyl_path;

	methyl_path = string().append(options.output_file.c_str()).append("_methyl_error.txt").c_str();

	FILE* methyl_graph = fopen(methyl_path,"w");

	unsigned int false_negatives=0;
	unsigned int false_positives=0;
	unsigned int true_negatives=0;
	unsigned int true_positives=0;

	unsigned int false_negatives_CG=0;
	unsigned int false_positives_CG=0;
	unsigned int true_negatives_CG=0;
	unsigned int true_positives_CG=0;

	unsigned int false_negatives_CHG=0;
	unsigned int false_positives_CHG=0;
	unsigned int true_negatives_CHG=0;
	unsigned int true_positives_CHG=0;

	unsigned int false_negatives_CHH=0;
	unsigned int false_positives_CHH=0;
	unsigned int true_negatives_CHH=0;
	unsigned int true_positives_CHH=0;

	unsigned int false_negatives_FW=0;
	unsigned int false_positives_FW=0;
	unsigned int true_negatives_FW=0;
	unsigned int true_positives_FW=0;

	unsigned int false_negatives_RC=0;
	unsigned int false_positives_RC=0;
	unsigned int true_negatives_RC=0;
	unsigned int true_positives_RC=0;

	unsigned int totC_CG=0;
	unsigned int totC_CHG=0;
	unsigned int totC_CHH=0;

	unsigned int metCG=0;
	unsigned int metCHG=0;
	unsigned int metCHH=0;//validation file levels
	char context[3];

	unsigned int errSize=51;
	unsigned int* error = new unsigned int[errSize];
	for(unsigned int i=0;i<errSize;i++)
		error[i]=0;

	unsigned int shift = 0;

	for (vector<string>::const_iterator input = options.validation_files.begin(); input != options.validation_files.end(); input++) {

		Auto_Unzip file(*input);
		istream & fasta = file.filtered();

		while(not fasta.eof()) {

			Fasta contig;
			fasta >> contig;

			unsigned int contigSize =contig.length();

			unsigned int i=0;
			for (string::const_iterator iter = contig.get_sequence().begin(); iter != contig.get_sequence().end(); iter++) {

				context[i%3]=*iter;

				if(i > 1){

					if(context[(i-2)%3]=='M'){//methylated in FW (+ strand)

						if(context[(i-1)%3]=='G' || context[(i-1)%3]=='L'){//CG
							metCG++;
							totC_CG++;
						}
						else if(context[i%3]=='G' || context[i%3]=='L'){//CHG
							metCHG++;
							totC_CHG++;
						}
						else{
							metCHH++;
							totC_CHH++;
						}
					}
					if(context[i%3]=='L'){//methylated in RC (- strand)

						if(context[(i-1)%3]=='C' || context[(i-1)%3]=='M'){//CG
							metCG++;
							totC_CG++;
						}
						else if(context[(i-2)%3]=='C' || context[(i-2)%3]=='M'){//CHG
							metCHG++;
							totC_CHG++;
						}
						else{
							metCHH++;
							totC_CHH++;
						}
					}
					if(context[(i-2)%3]=='C'){//unmethylated in FW (+ strand)

						if(context[(i-1)%3]=='G' || context[(i-1)%3]=='L')//CG
							totC_CG++;
						else if(context[i%3]=='G' || context[i%3]=='L')//CHG
							totC_CHG++;
						else
							totC_CHH++;
					}
					if(context[i%3]=='G'){//unmethylated in RC (- strand)

						if(context[(i-1)%3]=='C' || context[(i-1)%3]=='M')//CG
							totC_CG++;
						else if(context[(i-2)%3]=='C' || context[(i-2)%3]=='M')//CHG
							totC_CHG++;
						else
							totC_CHH++;
					}
				}

				if(covg(shift+i)>0 || options.coverage_threshold == 0)//if options.coverage_threshold == 0 then all the Cs are validated
					switch(*iter){

						//C non metilata su FW: idealmente tutte le C devono essere diventate T
						case 'C': case 'c':

							error[(unsigned int)(methyl_level(shift+i)*(errSize-1))]++;

							if(context_in_text(shift+i)==cg){
								if(methyl_level(shift+i) > 0.5){
									false_positives_CG++;
								}
								else{
									true_negatives_CG++;
								}
							}

							if(context_in_text(shift+i)==chg){
								if((methyl_level(shift+i)) > 0.5){
									false_positives_CHG++;
								}
								else{
									true_negatives_CHG++;
								}
							}

							if(context_in_text(shift+i)==chh){
								if((methyl_level(shift+i)) > 0.5){
									false_positives_CHH++;
								}
								else{
									true_negatives_CHH++;
								}
							}

							if((methyl_level(shift+i)) > 0.5){
								false_positives++;
								false_positives_FW++;
							}
							else{
								true_negatives++;
								true_negatives_FW++;
							}

						break;//case C

						//C non metilata su RC: idealmente tutte le C devono essere diventate T
						case 'G': case 'g':

							error[(unsigned int)(methyl_level(shift+i)*(errSize-1))]++;

							if(context_in_text(shift+i)==cg){
								if((methyl_level(shift+i)) > 0.5){
									false_positives_CG++;
								}
								else{
									true_negatives_CG++;
								}
							}

							if(context_in_text(shift+i)==chg){
								if((methyl_level(shift+i)) > 0.5){
									false_positives_CHG++;
								}
								else{
									true_negatives_CHG++;
								}
							}

							if(context_in_text(shift+i)==chh){
								if((methyl_level(shift+i)) > 0.5){
									false_positives_CHH++;
								}
								else{
									true_negatives_CHH++;
								}
							}

							if((methyl_level(shift+i)) > 0.5){
								false_positives++;
								false_positives_RC++;
							}
							else{
								true_negatives++;
								true_negatives_RC++;
							}

						break;//case G

						//C metilata su FW: idealmente tutte le C devono essere rimaste C
						case 'M': case 'm':

							error[(unsigned int)((1-methyl_level(shift+i))*(errSize-1))]++;

							if(context_in_text(shift+i)==cg){
								if((methyl_level(shift+i)) > 0.5){
									true_positives_CG++;
								}
								else{
									false_negatives_CG++;
								}
							}

							if(context_in_text(shift+i)==chg){
								if((methyl_level(shift+i)) > 0.5){
									true_positives_CHG++;
								}
								else{
									false_negatives_CHG++;
								}
							}

							if(context_in_text(shift+i)==chh){
								if((methyl_level(shift+i)) > 0.5)
									true_positives_CHH++;
								else{
									false_negatives_CHH++;
								}
							}

							if((methyl_level(shift+i)) > 0.5){
								true_positives++;
								true_positives_FW++;
							}
							else{
								false_negatives++;
								false_negatives_FW++;
							}

						break;//case M

						//C metilata su RC: idealmente tutte le C devono essere rimaste C
						case 'L': case 'l':

							error[(unsigned int)((1-methyl_level(shift+i))*(errSize-1))]++;

							if(context_in_text(shift+i)==cg){
								if((methyl_level(shift+i)) > 0.5){
									true_positives_CG++;
								}
								else{
									false_negatives_CG++;
								}
							}

							if(context_in_text(shift+i)==chg){
								if((methyl_level(shift+i)) > 0.5){
									true_positives_CHG++;
								}
								else{
									false_negatives_CHG++;
								}
							}

							if(context_in_text(shift+i)==chh){
								if((methyl_level(shift+i)) > 0.5){
									true_positives_CHH++;
								}
								else{
									false_negatives_CHH++;
								}
							}

							if((methyl_level(shift+i)) > 0.5){
								true_positives++;
								true_positives_RC++;
							}
							else{
								false_negatives++;
								false_negatives_RC++;
							}

						break;//case L

					}
				i++;
			}
			shift = shift + contigSize + 1;
		}
	}

	VERBOSE_CHANNEL << "VALIDATION RESULTS:" << endl <<endl;

	if(options.verbose){

		VERBOSE_CHANNEL << "methylated C in CG context (in validation file): " << (((double)metCG/(double)(totC_CG))*100) << "%" << endl;
		VERBOSE_CHANNEL << "methylated C in CHG context (in validation file): " << (((double)metCHG/(double)(totC_CHG))*100)<< "%" << endl;
		VERBOSE_CHANNEL << "methylated C in CHH context (in validation file): " << (((double)metCHH/(double)(totC_CHH))*100)<< "%" << endl;

		VERBOSE_CHANNEL << "\nspecificity and sensitivity - CG contexts" <<endl;

		VERBOSE_CHANNEL << endl << "\tnumber of true positives: " << true_positives_CG << endl;
		VERBOSE_CHANNEL << "\tnumber of true negatives: " << true_negatives_CG << endl;
		VERBOSE_CHANNEL << "\tnumber of false positives: " << false_positives_CG << endl;
		VERBOSE_CHANNEL << "\tnumber of false negatives: " << false_negatives_CG << endl<<endl;

		VERBOSE_CHANNEL << "\tSensivity: "<< (double)true_positives_CG/(double)(true_positives_CG+false_negatives_CG) << endl;
		VERBOSE_CHANNEL << "\tSpecificity: "<< (double)true_negatives_CG/(double)(true_negatives_CG+false_positives_CG) << endl<<endl;

		VERBOSE_CHANNEL << "Specificity and sensitivity - CHG contexts" <<endl;

		VERBOSE_CHANNEL << endl << "\tnumber of true positives: " << true_positives_CHG << endl;
		VERBOSE_CHANNEL << "\tnumber of true negatives: " << true_negatives_CHG << endl;
		VERBOSE_CHANNEL << "\tnumber of false positives: " << false_positives_CHG << endl;
		VERBOSE_CHANNEL << "\tnumber of false negatives: " << false_negatives_CHG << endl<<endl;

		VERBOSE_CHANNEL << "\tSensivity: "<< (double)true_positives_CHG/(double)(true_positives_CHG+false_negatives_CHG) << endl;
		VERBOSE_CHANNEL << "\tSpecificity: "<< (double)true_negatives_CHG/(double)(true_negatives_CHG+false_positives_CHG) << endl<<endl;

		VERBOSE_CHANNEL << "Specificity and sensitivity - CHH contexts" <<endl;

		VERBOSE_CHANNEL << endl << "\tnumber of true positives: " << true_positives_CHH << endl;
		VERBOSE_CHANNEL << "\tnumber of true negatives: " << true_negatives_CHH << endl;
		VERBOSE_CHANNEL << "\tnumber of false positives: " << false_positives_CHH << endl;
		VERBOSE_CHANNEL << "\tnumber of false negatives: " << false_negatives_CHH << endl<<endl;

		VERBOSE_CHANNEL << "\tSensivity: "<< (double)true_positives_CHH/(double)(true_positives_CHH+false_negatives_CHH) << endl;
		VERBOSE_CHANNEL << "\tSpecificity: "<< (double)true_negatives_CHH/(double)(true_negatives_CHH+false_positives_CHH) << endl<<endl;

		VERBOSE_CHANNEL << "specificity and sensitivity - positive (+) strand" <<endl;

		VERBOSE_CHANNEL << endl << "\tnumber of true positives: " << true_positives_FW << endl;
		VERBOSE_CHANNEL << "\tnumber of true negatives: " << true_negatives_FW << endl;
		VERBOSE_CHANNEL << "\tnumber of false positives: " << false_positives_FW << endl;
		VERBOSE_CHANNEL << "\tnumber of false negatives: " << false_negatives_FW << endl<<endl;

		VERBOSE_CHANNEL << "\tSensivity: "<< (double)true_positives_FW/(double)(true_positives_FW+false_negatives_FW) << endl;
		VERBOSE_CHANNEL << "\tSpecificity: "<< (double)true_negatives_FW/(double)(true_negatives_FW+false_positives_FW) << endl<<endl;

		VERBOSE_CHANNEL << "specificity and sensitivity - negative (-) strand" <<endl;

		VERBOSE_CHANNEL << endl << "\tnumber of true positives: " << true_positives_RC << endl;
		VERBOSE_CHANNEL << "\tnumber of true negatives: " << true_negatives_RC << endl;
		VERBOSE_CHANNEL << "\tnumber of false positives: " << false_positives_RC << endl;
		VERBOSE_CHANNEL << "\tnumber of false negatives: " << false_negatives_RC << endl<<endl;

		VERBOSE_CHANNEL << "\tSensivity: "<< (double)true_positives_RC/(double)(true_positives_RC+false_negatives_RC) << endl;
		VERBOSE_CHANNEL << "\tSpecificity: "<< (double)true_negatives_RC/(double)(true_negatives_RC+false_positives_RC) << endl<<endl;

	}

	VERBOSE_CHANNEL << "specificity and sensitivity - GLOBAL" <<endl;

	VERBOSE_CHANNEL << endl << "\tnumber of true positives: " << true_positives << endl;
	VERBOSE_CHANNEL << "\tnumber of true negatives: " << true_negatives << endl;
	VERBOSE_CHANNEL << "\tnumber of false positives: " << false_positives << endl;
	VERBOSE_CHANNEL << "\tnumber of false negatives: " << false_negatives << endl<<endl;

	VERBOSE_CHANNEL << "\tSensivity: "<< (double)true_positives/(double)(true_positives+false_negatives) << endl;
	VERBOSE_CHANNEL << "\tSpecificity: "<< (double)true_negatives/(double)(true_negatives+false_positives) << endl<<endl;

	for(unsigned int i=0;i<errSize;i++) {
		std::ostringstream s;
		s << (double)i/50 << "\t" << error[i] << endl;
		fputs(s.str().data(),methyl_graph);
	}
	fclose(methyl_graph);

	delete [] error;

}//validation


unsigned int Module_BS5::global_coord(unsigned int tid, unsigned int pos){
	return H->globalToLocal.startPositions[tid] + pos;
}//global_coord

unsigned int Module_BS5::covg(unsigned int i){
	return C_coverage(i) + T_coverage(i);
}

unsigned int Module_BS5::C_coverage(unsigned int i){
	return (unsigned int)(methylConsensus_C[i]);
}
unsigned int Module_BS5::T_coverage(unsigned int i){
	return (unsigned int)(methylConsensus_T[i]);
}

void Module_BS5::increment_C(unsigned int i){
	if(C_coverage(i)<255)
		methylConsensus_C[i]++;
}
void Module_BS5::increment_T(unsigned int i){
	if(T_coverage(i)<255)
		methylConsensus_T[i]++;
}

void Module_BS5::increment_mult_coverage(unsigned int i){
	if(mult_coverage[i]<255)
		mult_coverage[i]++;
}

double Module_BS5::methyl_level(unsigned int i){
	return (T_coverage(i)+C_coverage(i))==0?0:(double)(C_coverage(i))/(double)(T_coverage(i)+C_coverage(i));
}

double Module_BS5::alpha(char c,unsigned int shift){

	double alpha = 0;

	if(H->TEXT[shift] == 'G' && c == 'A')
		alpha = methyl_level(shift);//RC strand
	else if (H->TEXT[shift] == 'G' && c == 'G')
		alpha = 1-methyl_level(shift);

	if(H->TEXT[shift] == 'C' && c == 'T')
		alpha = methyl_level(shift);//FW strand
	else if (H->TEXT[shift] == 'C' && c == 'C')
		alpha = 1-methyl_level(shift);

	return alpha;
}//alpha

}//namespace modules
