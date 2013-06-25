#include <omp.h>
#include "rna_server.h"

//COLORS
#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"


#define OP_TYPE 0
#define REFERENCE_LEFT 1
#define REFERENCE_MIDLE 2
#define REFERENCE_INTRON 3
#define REFERENCE_RIGHT 4

#define MINIMUN_CAL_LENGTH 10
#define ERRORS_ZONE 8
#define MINIMUM_INTRON_LENGTH 10
#define MAX_CIGAR_OPERATIONS 20

//===============================================

#define A_NT  0
#define C_NT  1
#define G_NT  2
#define T_NT  3

const unsigned char splice_nt[4] = {'A', 'C', 'G', 'T'};

//--------------------------------------------//
//            Not found splice junction       //
//--------------------------------------------//

#define NOT_SPLICE	-1

//--------------------------------------------//
//      No Cannonical Splice junction         //
//--------------------------------------------//

#define UNKNOWN_SPLICE	0
                                            
//--------------------------------------------//
//        Cannonical Splice Junction          //
//--------------------------------------------//

#define GT_AG_SPLICE  	1 //+
#define CT_AC_SPLICE  	2 //-
  
//--------------------------------------------//
//      Semi-Cannonical Splice Junction       //
//--------------------------------------------//

#define AT_AC_SPLICE  	3 //+
#define GT_AT_SPLICE  	4 //-
#define GC_AG_SPLICE  	5 //+
#define CT_GC_SPLICE  	6 //-

//===============================================


#ifndef MAX
   #define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

//======== SPLICE JUNCTION TYPE ============
#define NOT_SPLICE	0
#define GT_AG_SPLICE  	1
#define CT_AC_SPLICE  	2
//==========================================

#define SEARCH_START 0
#define SEARCH_END 1
#define POSSIBLES_MARKS 2
#define STRANDS_NUMBER 2
#define MAX_FUSION 2028
#define MIN_HARD_CLIPPING 10

int ii = -1;

void cal_fusion_data_init(unsigned int id, size_t start, size_t end, 
			  unsigned int strand, unsigned int chromosome, unsigned int fusion_start, 
			  unsigned int fusion_end, cal_fusion_data_t *cal_fusion_data_p) {

  cal_fusion_data_p->genome_strand = strand;
  cal_fusion_data_p->id = id;
  cal_fusion_data_p->genome_start = start;
  cal_fusion_data_p->genome_chromosome = chromosome;
  cal_fusion_data_p->genome_end = end;
  cal_fusion_data_p->fusion_start = fusion_start;
  cal_fusion_data_p->fusion_end = fusion_end;

}

char cigar_automata_status(unsigned char status) {

  switch (status) {
    case CIGAR_MATCH_MISMATCH:
      return 'M';
      break;
    case CIGAR_INSERTION:
      return 'I';
      break;
    case CIGAR_DELETION:
      return 'D';
      break;
    case CIGAR_SKIPPED:
      return 'N';
      break;
    case CIGAR_PADDING:
      return 'P';
      break;
  }

  return ' ';
}

int cigar_generator(cigar_data_t *cigar_p, unsigned int *max_len, 
		    unsigned int *pos, unsigned int status, int *number, int *distance){
  char cigar_info[20];
  //printf("Store %d in %i position\n", *number, *pos);
  /*if (*pos > MAX_CIGAR_OPERATIONS) {
    return 2;
    }*/
  
  if (*number <= 0) {
    return 1;
  }else{
    //printf("**********Allocate %d\n",*number );
    if ((status != CIGAR_DELETION) || ((status == CIGAR_DELETION) && (*number <= MINIMUM_INTRON_LENGTH))){
	  cigar_p[*pos].type = status;
	  cigar_p[*pos].value = *number;	  
	  *pos += 1;
	  if (status != CIGAR_PADDING && status != CIGAR_MATCH_MISMATCH) {
	    *distance += *number;
	  }
	  *number = 1;
	  if (*pos >= *max_len){
	    *max_len = *max_len * 2;
	    cigar_p = (cigar_data_t *)realloc(cigar_p, sizeof(cigar_data_t)*(*max_len));
	  }
      }
  }
  return 0;
}

//list_t* write_list_p, write_batch_t* write_batch_p, unsigned int write_size,
/*
void search_splice_junctions_sw_output(sw_simd_input_t* input_p, sw_simd_output_t* output_p, 
				       unsigned int depth, allocate_fusion_data_t *depth_cal_fusion_p, 
				       allocate_splice_elements_t *chromosome_avls_p,  
				       sw_channel_t *sw_channels_p, mapping_batch_t *mapping_batch_p,  unsigned int sw_id,
				       size_t *sw_no_valids, float min_score, genome_t *genome_p, int min_intron_length, 
				       float gap_open, float gap_extend, float match) {


  //printf("Process Splice Junctions\n");
  const unsigned char EXTRA_SEARCH = 8;
  char header_id[2048];
  alignment_t *alignment_p;
  
  char *header_match, *read_match, *quality_match, *quality_strand;
  
  char extra_nt[EXTRA_SEARCH + 1];
  unsigned int extra_nt_len = 0;
  unsigned char ref_nt_pos;
  unsigned char extra_nt_found;
  int gap_start = -1; int gap_end = -1;
  unsigned char found = 0;
  unsigned char confirmation_type = 0;
  unsigned int l;
  unsigned char count_extra_search;
  unsigned int position;
  unsigned int start_splice, end_splice;
  unsigned char actual_cal = 0;
  unsigned int start_mapped = 0;
  unsigned int n_deletions_seq = 0;
  unsigned int insertion_number;
  unsigned int deletion_number;
  unsigned int deletions_tot;
  unsigned int displacement_start;
  unsigned int mark_dpl;
  unsigned int strand;
  unsigned char strand_splice;
  int cal;
  unsigned int relative_gap_start;
  unsigned int nt_skipped;
  size_t str_pos;
  unsigned int quality_pos;
  unsigned int reference_position_disp;
  unsigned int start_gap_padding;
  unsigned int start_search_j;

  unsigned int extend_start_sp, extend_end_sp;

  unsigned int splice_number;
  size_t maximum_splice_number = 100;
  sp_data_t *allocate_splice = (sp_data_t *)malloc(sizeof(sp_data_t)*maximum_splice_number);//[maximum_splice_number];

  //=========CIGAR INFO. && BAM/SAM ========//
  int error_cigar;
  int hard_clipping;
  short int large_hard_clipping;
  unsigned int cigar_soft;
  unsigned int value;
  int cigar_value;
  unsigned int cigar_max_len = 2048;
  unsigned int cigar_pos;
  unsigned char automata_status;
  char flag;
  short int num_cigar_op;
  unsigned int header_len, read_len, mapped_len, bytes;
  short int primary_alignment;
  unsigned int total_cals;
  char end_exceded = 0;
  unsigned int orig_pos = 0;
  unsigned int tmp;
  unsigned int real_sw_len = 0;
  unsigned int j_start;
  unsigned int change_cal = 0, change_cal_number = 0;
  char cigar_segment[10];
  //unsigned int change_cal_index[20];
  //list_item_t* item_p = NULL;  
  char *cigar_str;
  cigar_data_t *cigar_p = (cigar_data_t *)malloc(sizeof(cigar_data_t)*cigar_max_len);
  char *optional_fields, *p;
  int optional_fields_length;
  int distance;
  int AS;
  float insertions_penalty;

  if (cigar_p == NULL) { exit(-1); }
  //============================//
  
  
  //============ CANONICAL INTRON MARKS ==============//
  unsigned char canonical_GT_AG[4];
  unsigned char canonical_CT_AC[4];
  
  //START INTRON MARKS 'GT' or 'CT'
  canonical_GT_AG[0] = 'G';
  canonical_GT_AG[1] = 'T';
  canonical_CT_AC[0] = 'C';
  canonical_CT_AC[1] = 'T';
  //END INTRON MARKS 'AG' or 'AC'
  canonical_GT_AG[2] = 'A';
  canonical_GT_AG[3] = 'G';
  canonical_CT_AC[2] = 'A';
  canonical_CT_AC[3] = 'C';
  //================================================//

  printf("======================== Process Output SW %d=========================\n", depth);
  sw_simd_input_display(depth, input_p);
  sw_simd_output_display(depth, output_p);
  printf("======================================================================\n");

  for (int i = 0; i < depth; i++) {    
    splice_number = 0;
    total_cals = depth_cal_fusion_p[i].cal_number;

    //printf("Total cals %d\n", total_cals);
    //if (total_cals == 1) {
    if (output_p->norm_score_p[i] < min_score) {
      *sw_no_valids += 1;
      continue;
    }
    AS =(int)output_p->score_p[i];
    large_hard_clipping = 0;
    fastq_read_t *read_p = array_list_get(sw_channels_p[i].read_index, mapping_batch_p->fq_batch);
    reference_position_disp = 0;

    //printf("Actualization start %d\n", output_p->start_p[i]);
    //====================== Actualization CALS ===========================//
    //Move to init CAL and calculate start padding
    cal = 0;
    displacement_start = 0;
    while (depth_cal_fusion_p[i].allocate_data[cal].fusion_end  <= output_p->start_p[i]) {
      displacement_start = depth_cal_fusion_p[i].allocate_data[cal].fusion_end + 1;
      depth_cal_fusion_p[i].allocate_data[cal].fusion_end = 0;
      depth_cal_fusion_p[i].allocate_data[cal].fusion_start = 0;
      cal++;
      //printf("Delete CAL %d\n", cal);
      if (cal >= total_cals) {
	//printf("ERROR: total cals overflow\n");
	continue;
      }
    }

    actual_cal = cal;
    
    displacement_start = output_p->start_p[i] - displacement_start;

    depth_cal_fusion_p[i].allocate_data[cal].genome_start += displacement_start;
        
    start_mapped = depth_cal_fusion_p[i].allocate_data[cal].genome_start;
    start_gap_padding = displacement_start;

    depth_cal_fusion_p[i].allocate_data[cal].fusion_start = 0;
    depth_cal_fusion_p[i].allocate_data[cal].fusion_end = depth_cal_fusion_p[i].allocate_data[cal].fusion_end - output_p->start_p[i];
    cal++;
    while (cal < total_cals) {
      depth_cal_fusion_p[i].allocate_data[cal].fusion_start = depth_cal_fusion_p[i].allocate_data[cal].fusion_start - output_p->start_p[i];
      depth_cal_fusion_p[i].allocate_data[cal].fusion_end = depth_cal_fusion_p[i].allocate_data[cal].fusion_end - output_p->start_p[i];
      cal++;
    }
    
    cal = total_cals - 1;
    while ((depth_cal_fusion_p[i].allocate_data[cal].fusion_start + 1) > output_p->mapped_len_p[i]) {
      depth_cal_fusion_p[i].allocate_data[cal].fusion_end = 0;
      depth_cal_fusion_p[i].allocate_data[cal].fusion_start = 0;
      cal--;
      if (cal < 0) {
	//printf("ERROR: total cals overflow\n");
	continue;
      }
    }

    depth_cal_fusion_p[i].allocate_data[cal].genome_end -= (depth_cal_fusion_p[i].allocate_data[cal].fusion_end - output_p->mapped_len_p[i]);
    depth_cal_fusion_p[i].allocate_data[cal].fusion_end  = output_p->mapped_len_p[i] - 1;
        
    j_start = 0;
    while (((depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end - depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_start) < 
	   MINIMUN_CAL_LENGTH) && (actual_cal < total_cals)) {
      j_start += (depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end - depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_start) + 1;
      depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_start = 0;
      depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end = 0;
      actual_cal++;
      if (actual_cal >= total_cals) {
	//printf("ERROR: total cals overflow\n");
	continue;
      }
      start_mapped = depth_cal_fusion_p[i].allocate_data[actual_cal].genome_start;      
      //printf("Move to next CAL %i, start_j = %i\n", actual_cal, j_start);
    }

    while (output_p->mapped_ref_p[i][j_start] != output_p->mapped_seq_p[i][j_start]) {
      j_start++;
      if (j_start >= output_p->mapped_len_p[i]){ continue; }
    }

    cal = actual_cal;
    change_cal_number = 0;
    real_sw_len = 0;

    //memset( change_cal_index, 0, sizeof(unsigned int)*total_cals);    
    //printf("Actual CAL %d", actual_cal);
    //printf("=========================== AFTER [Depth %d - Total CALs %d]==================================\n", i, depth_cal_fusion_p[i].cal_number);
    //for(int c = 0; c < depth_cal_fusion_p[i].cal_number; c++){
      //printf("CAL %d :\n", c);
      //printf("\t->Genome Start : %d\n", depth_cal_fusion_p[i].allocate_data[c].genome_start);
      //printf("\t->Genome End : %d\n", depth_cal_fusion_p[i].allocate_data[c].genome_end);
      //printf("\t->Fusion Start : %d\n", depth_cal_fusion_p[i].allocate_data[c].fusion_start);
      //printf("\t->Fusion End : %d\n", depth_cal_fusion_p[i].allocate_data[c].fusion_end);
      //printf("\t->Chromosome : %d\n", depth_cal_fusion_p[i].allocate_data[c].genome_chromosome);
      //printf("\t->Strand : %d\n", depth_cal_fusion_p[i].allocate_data[c].genome_strand);
    //}
    //printf("=============================================================================================================\n");
    
    //====================== Actualization CALS End ===========================//

    //Data Info
    position = 0;
    count_extra_search = 0;
    start_splice = 0;
    end_splice = 0;
    insertion_number = 0;
    deletion_number = 0;
    mark_dpl=0;
    relative_gap_start = 0;
    n_deletions_seq = 0;
    orig_pos = 0;
    deletions_tot = 0;    
    str_pos = 0;
    distance = 0;
    optional_fields = (char *)malloc(sizeof(char)*100);

    read_len = output_p->mapped_len_p[i];
    read_match = (char *)malloc(sizeof(char)*(read_len + 1));
    if (read_match == NULL) { printf("Error to allocate memory\n"); exit(-1); }
   
    quality_pos = output_p->start_seq_p[i];
    quality_match = (char *)malloc(sizeof(char)*(read_len + 1));
    if (quality_match == NULL) { printf("Error to allocate memory\n"); exit(-1); }
    
    found = NOT_SPLICE;
    //===========Position Automata in initial State==========
    cigar_value = 0;
    cigar_pos = 0;
    
    if (output_p->mapped_ref_p[i][j_start] != '-') {
      if (output_p->mapped_seq_p[i][j_start] != '-') {
	automata_status = CIGAR_MATCH_MISMATCH;
      } else {
	automata_status = CIGAR_DELETION;
      }
    } else {
      automata_status = CIGAR_INSERTION;
    }
    //========================================================
    
    change_cal = 0;
    //printf("\tStart Search\n");
    for (int j = j_start; j < output_p->mapped_len_p[i]; j++) {
      if (output_p->mapped_ref_p[i][j] != '-') {
	reference_position_disp++;
	start_gap_padding++;
	relative_gap_start++;
	//TODO: IMPLEMENTS WHEN JUMP TO OTHER CALL WITHOUT GAP
	if (output_p->mapped_seq_p[i][j] != '-') {
	  //Match mismatch area	  
	  if (output_p->mapped_seq_p[i][j] != output_p->mapped_ref_p[i][j]) { distance++;  } 
	  if (automata_status == CIGAR_MATCH_MISMATCH) {
	    cigar_value++;
	    //printf("Study %c nt %d\n", output_p->mapped_ref_p[i][j], cigar_value);
	  } else {
	    //printf("%s: %d in [j=%d - jstar=%d]\n", cigar_automata_status(automata_status), cigar_value, j, j_start);
	    error_cigar = cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, automata_status, &cigar_value, &distance);
	    if (error_cigar) { end_exceded = 1; break;}//printf("MATCH-MISMATCH :: %s\n", read_p->id);printf("Error cigar\n"); }
	    automata_status = CIGAR_MATCH_MISMATCH;
	  }
	}
      } else {
	//Insertion area
	orig_pos++;
	insertion_number++;
	if (automata_status == CIGAR_INSERTION) {
	   cigar_value++;
	} else {
	  //printf("%s: %d in [j=%d - jstar=%d]\n", cigar_automata_status(automata_status), cigar_value, j, j_start);
	  //distance += cigar_value;
	  error_cigar = cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, automata_status, &cigar_value, &distance);
	  if (error_cigar) { end_exceded = 1; break;}//if (error_cigar) { printf("INSERTION :: %s\n", read_p->id);printf("Error cigar"); }
	  automata_status = CIGAR_INSERTION;
	}
      }//end mapped_ref_p[i][j] != '-'

      if (output_p->mapped_seq_p[i][j] == '-') {
	if (gap_start == -1) {
	  gap_start = j;
	}

	n_deletions_seq++;
	deletions_tot++;
	
	if (output_p->mapped_ref_p[i][j] == '-') {
	  //Padding area
	  if (automata_status == CIGAR_PADDING) {
	    cigar_value++;
	  } else {
	    error_cigar = cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, automata_status, &cigar_value, &distance);
	     if (error_cigar) { end_exceded = 1; break;}//if (error_cigar) { printf("%s\n", read_p->id);printf("Error cigar"); }
	    automata_status = CIGAR_PADDING;    
	  }
	} else {
	  //Deletion area
	  if (automata_status == CIGAR_DELETION) {
	    cigar_value++;
	  } else {
	    //printf("%s: %d in [j=%d - jstar=%d]\n", cigar_automata_status(automata_status), cigar_value, j, j_start);
	    //distance += cigar_value;
	    error_cigar = cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, automata_status, &cigar_value, &distance);
	    if (error_cigar) { end_exceded = 1; break;}//if (error_cigar) { printf("DELETION :: %s\n", read_p->id);printf("Error cigar"); }
	    automata_status = CIGAR_DELETION;
	  }
	}
      } else { 
	 read_match[str_pos] = output_p->mapped_seq_p[i][j];
	 quality_match[str_pos++] = read_p->quality[quality_pos];//sw_batch_p->allocate_reads_p[sw_channels_p[i].read_index]->quality[quality_pos];
	 quality_pos++;
	 
	 if (n_deletions_seq > MINIMUM_INTRON_LENGTH) {
	    
	   gap_end = j - 1;   
	   //printf("gap_start=%d , gap_end=%d\n", gap_start, gap_end);

	   //Canonical marks strand + (GT-AG)
	   if ((output_p->mapped_ref_p[i][gap_start] == canonical_GT_AG[0]) && 
	       (output_p->mapped_ref_p[i][gap_start + 1] == canonical_GT_AG[1])) {
	     if ((output_p->mapped_ref_p[i][gap_end - 1] == canonical_GT_AG[2]) && 
		 (output_p->mapped_ref_p[i][gap_end] == canonical_GT_AG[3])) {
	       found = GT_AG_SPLICE; 
	     }
	   }
	   //Canonical marks strand - (CT-AC)
	   if ((output_p->mapped_ref_p[i][gap_start] == canonical_CT_AC[0]) && 
	       (output_p->mapped_ref_p[i][gap_start + 1] == canonical_CT_AC[1])) {
	     if ((output_p->mapped_ref_p[i][gap_end - 1] == canonical_CT_AC[2]) && 
		 (output_p->mapped_ref_p[i][gap_end] == canonical_CT_AC[3])) {
	       found = CT_AC_SPLICE;
	     }
	   }
	   
	   if (!found) {
	     //Search Xnt (+)---->	    
	     confirmation_type = NOT_SPLICE;
	     count_extra_search = 1;
	     position = gap_start + 1;
	     extra_nt[count_extra_search - 1] = output_p->mapped_ref_p[i][position - 1];
	     //printf("Extra Search position %d and gap end %d\n", position, gap_end);
	     while (count_extra_search < EXTRA_SEARCH ) {
	       if ((gap_end + count_extra_search) >= output_p->mapped_len_p[i]) {
		 end_exceded = 1; break;
	       }

	       if ((output_p->mapped_ref_p[i][position] == canonical_GT_AG[0]) && 
		   (output_p->mapped_ref_p[i][position + 1] == canonical_GT_AG[1])) {
		 //printf("-Found GT-");
		 if ((output_p->mapped_ref_p[i][gap_end - 1 + count_extra_search] == canonical_GT_AG[2]) && 
		     (output_p->mapped_ref_p[i][gap_end + count_extra_search] == canonical_GT_AG[3])) {
		   confirmation_type = GT_AG_SPLICE; 
		   //printf("Found GT_AG Extra Search\n");
		 }
		 
	       } else if ((output_p->mapped_ref_p[i][position] == canonical_CT_AC[0]) && 
			  (output_p->mapped_ref_p[i][position + 1] == canonical_CT_AC[1])) {
		 if ((output_p->mapped_ref_p[i][gap_end - 1 + count_extra_search] == canonical_CT_AC[2]) && 
		     (output_p->mapped_ref_p[i][gap_end + count_extra_search] == canonical_CT_AC[3])) {
		   confirmation_type = CT_AC_SPLICE; 
		   //printf("Found CT_AC Extra Search\n");
		 }
	       }
	       
	       if (confirmation_type) {
		 extra_nt[count_extra_search] = '\0';
		 extra_nt_len = strlen(extra_nt);
		 ref_nt_pos = gap_end + 1;
		 extra_nt_found = 1;
		 
		 for (int k = 0; k < extra_nt_len; k++) {
		   if (extra_nt[k] != output_p->mapped_ref_p[i][ref_nt_pos]) {
		     extra_nt_found = 0;
		     break;
		   }
		   ref_nt_pos++;
		 }
		
		 if (extra_nt_found) {
		   start_splice = position;
		   end_splice = gap_end + count_extra_search - 1;
		   found= confirmation_type;
		   mark_dpl = count_extra_search;
		   confirmation_type = NOT_SPLICE;
		 }
	       }
	       
	       count_extra_search++;
	       position++;
	       extra_nt[count_extra_search - 1] = output_p->mapped_ref_p[i][position - 1];
	     }//End while extra search
	     //printf("Extra search end\n");
	   }//End if not found
	   
	   if (found) {	    
	     //printf("====================Found Splice Junction %d Splice number========================\n", splice_number);
	     relative_gap_start--;
	     
	     relative_gap_start -= n_deletions_seq;//((gap_end - gap_start) + 1);
	     start_splice = depth_cal_fusion_p[i].allocate_data[actual_cal].genome_start + (relative_gap_start + mark_dpl);
	     extend_start_sp = depth_cal_fusion_p[i].allocate_data[actual_cal].genome_start;
	    
	     allocate_splice[splice_number].start_sp = start_splice;
	
	     //printf("Before.Value %d\n", allocate_splice[0].start_extend_sp);     

	     if (splice_number >= 1) {
	      allocate_splice[splice_number].start_extend_sp = allocate_splice[splice_number - 1].end_sp;
	      allocate_splice[splice_number - 1].end_extend_sp = start_splice;
	     } else {
	       allocate_splice[splice_number].start_extend_sp = extend_start_sp;
	     }

	     //printf("After.Value %d\n", allocate_splice[0].start_extend_sp);     

	    //printf("(+)START SPLICE :: %d + (%d + %d) :: %d\n", depth_cal_fusion_p[i].allocate_data[actual_cal].genome_start, relative_gap_start, mark_dpl, start_splice);
	   
	    if (gap_end >  depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end) {
	      //printf("\tChange to Next CAL\n");	     
	      tmp = ((depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end - j_start) - 
		     (depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_start - j_start) + 1) - relative_gap_start ;
	      
	      relative_gap_start = n_deletions_seq - tmp;
	     
	      actual_cal++;
	      nt_skipped = 0;

	      if (gap_end >  depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end) {
		while (gap_end >  depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end) {
		  nt_skipped += (depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end - 
				 depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_start);
		  actual_cal++;
		}
		nt_skipped++;
	      }
	      
	      if (actual_cal >= total_cals){ end_exceded = 1; break; }
	      
	      relative_gap_start -= nt_skipped;
	      
	      start_gap_padding = 0;
	      
	      end_splice = depth_cal_fusion_p[i].allocate_data[actual_cal].genome_start + relative_gap_start + mark_dpl;
	      end_splice--;
	
	      relative_gap_start++;
	      insertion_number = 0;
	    } else {
	      //printf("\tStay in same CAL\n");
	      relative_gap_start += n_deletions_seq;
	      
	      end_splice = depth_cal_fusion_p[i].allocate_data[actual_cal].genome_start + (relative_gap_start - 1) + mark_dpl;
	    }

	    extend_end_sp = depth_cal_fusion_p[i].allocate_data[actual_cal].genome_end;
	    allocate_splice[splice_number].end_sp = end_splice;

	    if (found == CT_AC_SPLICE) {
	      strand_splice = 1;
	    } else {
	      strand_splice = 0;
	    }

	    allocate_splice[splice_number].strand_sp = strand_splice;
	    allocate_splice[splice_number].end_extend_sp = extend_end_sp;
	    cigar_value = (int)end_splice - (int)start_splice + 1;

	    if (cigar_value > min_intron_length) {
	      splice_number++;
	      if (splice_number >= maximum_splice_number) {
		maximum_splice_number *= 1.25;
		allocate_splice = (sp_data_t *)realloc(allocate_splice, sizeof(sp_data_t)*maximum_splice_number);
	      }
	      cigar_p[cigar_pos - 1].value += mark_dpl;
	      //cigar_p[cigar_pos + 1].value = -mark_dpl;
	      //printf("Report splice displacement %i, %i, %i\n", mark_dpl, cigar_pos, cigar_p[cigar_pos + 1].value);
	      //printf("%i - %i\n", cigar_p[cigar_pos - 1].value, mark_dpl);
	      error_cigar = cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, CIGAR_SKIPPED, &cigar_value, &distance);
	      //cigar_value += mark_dpl;
	      //printf("j=%i\n", j);
	      cigar_value -= mark_dpl;
	      //printf("start cigar value %d\n", cigar_value);
	      //j += mark_dpl;
	      if (error_cigar) { end_exceded = 1; break; }
	      //printf("::: %d || %d :::\n", cigar_value, mark_dpl);	      

	      //Add to the mapping the score of extend gap in splice junction search
	      //printf("::: %d || %d :::\n", cigar_value, mark_dpl);
	    } else {
	      //distance += cigar_value;
	      error_cigar = cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, CIGAR_BIG_DELETION, &cigar_value, &distance);
	       if (error_cigar) { end_exceded = 1; break; }
	    }
	    
	    found = NOT_SPLICE;
	    deletion_number = 0;
	    //printf("==================== Exact->( %d:%d-%d ) Extend->( %d:%d-%d )  ========================\n", depth_cal_fusion_p[i].allocate_data->genome_chromosome, start_splice, end_splice, depth_cal_fusion_p[i].allocate_data->genome_chromosome, extend_start_sp, extend_end_sp);
	  } else { //IF not found
	    if (reference_position_disp > depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end) {
	      actual_cal++;
	      //printf("SPLICE NOT FOUND report deletion\n");
	      if (actual_cal >= total_cals) {
		//Cals overflow, read no valid		
		end_exceded = 1;
		break;
	      } else {
		
		relative_gap_start = gap_end   - (depth_cal_fusion_p[i].allocate_data[actual_cal - 1].fusion_end + insertion_number);
		
		cigar_value = (depth_cal_fusion_p[i].allocate_data[actual_cal].genome_start + relative_gap_start) -  
		  (depth_cal_fusion_p[i].allocate_data[actual_cal - 1].genome_end - ((gap_end - gap_start) - relative_gap_start));
	      }
	      //error_cigar = cigar_generator(cigar_p, &cigar_max_len, &cigar_pos,  CIGAR_BIG_DELETION, &cigar_value, &distance);
	      //printf("No splice found %i\n", cigar_value);
	      error_cigar = cigar_generator(cigar_p, &cigar_max_len, &cigar_pos,  CIGAR_SKIPPED, &cigar_value, &distance);
	      if (error_cigar) { end_exceded = 1; break;}
	    } else {
	     cigar_value = gap_end - gap_start + 1;
	     //error_cigar = cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, CIGAR_BIG_DELETION, &cigar_value, &distance);
	     //printf("No splice found %i\n", cigar_value);
	     error_cigar = cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, CIGAR_SKIPPED, &cigar_value, &distance);
	     if (error_cigar) { end_exceded = 1; break;}
	    }
	    deletion_number = 0;
	    insertion_number = 0;
	    n_deletions_seq = 0;
	   //CIGAR_BIG_DELETION
	  }

	 } else {
	   //Deletion
	   deletion_number += n_deletions_seq;
	   n_deletions_seq = 0;
	 }//if minimum deletions mark 
	 gap_start = -1;
	 n_deletions_seq = 0;
      }//Not found
    }//End loop process smith-waterman

    //printf("\tEnd process sw %i\n", i);
    //printf("\tStart Search\n");
    if (end_exceded) { free(optional_fields); free(read_match); free(quality_match); end_exceded = 0; continue; }
    //printf("\tEnd Search Splice\n");

    //printf("last cigar value %i\n", cigar_value);
    //=============================== MAKE CIGAR & SAM/BAM =====================================
    cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, automata_status, &cigar_value, &distance);

    if (cigar_pos < 1) { free(optional_fields); free(read_match); free(quality_match); continue; }

    //Generate cigar string
    num_cigar_op = 0;
    cigar_str = (char *)malloc(sizeof(char)*cigar_max_len);
    if (cigar_str == NULL) { exit(-1); }
    cigar_str[0] = '\0';
    
    //Hard and Soft clipped start
    if (output_p->start_seq_p[i] > 0) {
      hard_clipping = output_p->start_seq_p[i] + j_start;
      sprintf(cigar_segment, "%iH", hard_clipping);
      cigar_str = strcat(cigar_str, cigar_segment);
      num_cigar_op++;
      if (MIN_HARD_CLIPPING <= hard_clipping) {
	large_hard_clipping = 1;
      }
    }
    
    if (cigar_p[0].type == CIGAR_MATCH_MISMATCH) {
      cigar_soft = 0;
      value = 0;
      while ((output_p->mapped_ref_p[i][cigar_soft] != '-') && (output_p->mapped_seq_p[i][cigar_soft] != '-') && 
	      (output_p->mapped_seq_p[i][cigar_soft] != output_p->mapped_ref_p[i][cigar_soft])) {
		cigar_soft++;
		value++;
		//printf("Soft (%c!=%c)", output_p->mapped_ref_p[i][cigar_soft - 1], output_p->mapped_seq_p[i][cigar_soft - 1]);
      }
      
      if (value > 0) {
	cigar_p[0].value -= cigar_soft;
	sprintf(cigar_segment, "%iS", value);
	cigar_str = strcat(cigar_str, cigar_segment);
	num_cigar_op++;
      }
    }
    
    //Hard and Soft clipped end
    //printf("Cigar value %d\n", cigar_pos);
    value = 0;
    if (cigar_p[cigar_pos - 1].type == CIGAR_MATCH_MISMATCH) {
      cigar_soft = output_p->mapped_len_p[i] - 1;
      while ((output_p->mapped_ref_p[i][cigar_soft] != '-') && (output_p->mapped_seq_p[i][cigar_soft] != '-') && 
	      (output_p->mapped_seq_p[i][cigar_soft] != output_p->mapped_ref_p[i][cigar_soft])) {
		cigar_soft--;
		value++;
		//printf("(Soft %c!=%c)", output_p->mapped_ref_p[i][cigar_soft], output_p->mapped_seq_p[i][cigar_soft]);
      }
      cigar_p[cigar_pos - 1].value -= value;
    }
    
    for (int cig = 0; cig < cigar_pos; cig++) {
      if (cigar_p[cig].type == CIGAR_MATCH_MISMATCH) {
	flag = 'M';
      } else if (cigar_p[cig].type == CIGAR_DELETION || cigar_p[cig].type == CIGAR_BIG_DELETION) {
	flag = 'D';
      } else if (cigar_p[cig].type == CIGAR_PADDING) {
	flag = 'P';
      } else if (cigar_p[cig].type == CIGAR_INSERTION) {
	flag = 'I';
      } else if (cigar_p[cig].type == CIGAR_SKIPPED) {
	flag = 'N';
      }
      //printf("%u\n", cigar_p[cig].value);
      sprintf(cigar_segment, "%u%c", cigar_p[cig].value, flag);
      //printf("%s\n", cigar_segment);
      cigar_str = strcat(cigar_str, cigar_segment);
      num_cigar_op++;
      //printf("%d%c", cigar_p[cig].value, flag);
    }
    
       
    if (value > 0) {
      sprintf(cigar_segment, "%iS", value);
      cigar_str = strcat(cigar_str, cigar_segment);
      num_cigar_op++;
     }
    
    if (((output_p->mapped_len_p[i] - deletions_tot) + output_p->start_seq_p[i]) < input_p->seq_len_p[i]) {
      hard_clipping = input_p->seq_len_p[i] - ((output_p->mapped_len_p[i] - deletions_tot) + output_p->start_seq_p[i]);
      sprintf(cigar_segment, "%iH", hard_clipping);
      cigar_str = strcat(cigar_str, cigar_segment);
      num_cigar_op++;
      if (MIN_HARD_CLIPPING <= hard_clipping) {
	large_hard_clipping = 2;
      }
    }
    
    //printf("Cigar(%d):%s\n", cigar_pos, cigar_str);
    //end cigar string
//    if (array_list_size(mapping_batch_p->mapping_lists[sw_channels_p[i].read_index])) {
    if (array_list_size(mapping_batch_p->mapping_lists[sw_channels_p[i].read_index]) > 1) {
      primary_alignment = 1;
    } else {
      primary_alignment = 0;
    }

    alignment_p = alignment_new();
    header_len = sw_channels_p[i].header_len;
    //printf("Reallocate %d\n", header_len); 
    get_to_first_blank(read_p->id, header_len, header_id);
    header_match = (char *)malloc(sizeof(char)*header_len);
    if (header_match == NULL) { exit(-1); }
    memcpy(header_match, header_id, header_len);
    
    //printf("Header len = %d\n", header_len);

    read_match[str_pos] = '\0';
    //quality_match[str_pos] = '\0';
    //memcpy(quality_match, &(sw_batch_p->quality_p[sw_batch_p->read_indices_p[sw_channels_p[i].read_index]]), read_len );
    if (depth_cal_fusion_p[i].allocate_data->genome_strand) {
      quality_strand = (char *)malloc(sizeof(char)*(str_pos + 1));
      reverse_str(quality_match, quality_strand, str_pos);
      //memcpy(quality_strand, quality_match, str_pos);
      free(quality_match);
      quality_strand[str_pos] = '\0';
      quality_match = quality_strand;
    }

    //printf("CIGAR: %s\n", cigar_str);
    p = optional_fields;
    sprintf(p, "ASi");
    p += 3;
    memcpy(p, &AS, sizeof(int));
    p += sizeof(int);

    sprintf(p, "NMi");
    p += 3;
    memcpy(p, &distance, sizeof(int));
    p += sizeof(int);

    optional_fields_length = p - optional_fields;
    //printf("Start mapped %i\n", start_mapped);
    alignment_init_single_end(header_match, read_match, quality_match, 
			      depth_cal_fusion_p[i].allocate_data->genome_strand, 
			      depth_cal_fusion_p[i].allocate_data->genome_chromosome - 1, start_mapped - 1, 
			      cigar_str, num_cigar_op, output_p->norm_score_p[i] * 254, 1, 
			      primary_alignment, optional_fields_length, optional_fields, large_hard_clipping, 
			      alignment_p);

    //printf("Score %i\n", alignment_p->map_quality);
    //printf("seq: %s\n", alignment_p->sequence);
    //printf("Insert mapping\n");
    

    //--------->
    //printf("-->Allocate mapping in %i position\n", sw_channels_p[i].read_index);
    array_list_insert(alignment_p, mapping_batch_p->mapping_lists[sw_channels_p[i].read_index]);
    //if (ii == sw_channels_p[i].read_index){ printf("SW READ list size %i - read position %i\n", array_list_size(allocate_mappings), sw_channels_p[i].read_index); 
    //    alignment_print(alignment_p);
    

    //============================================================================================
    
    //printf(" (score = %.2f, norm. score = %.2f)\n", output_p->score_p[i], output_p->norm_score_p[i]);
   
    //printf("Free done 1\n");
    
    //mapping_reads_p[sw_channels_p[i].read_index] = 1;
    
    //Report Splice Junctions
    for (unsigned int s = 0; s < splice_number; s++) {
      pthread_mutex_lock(&(chromosome_avls_p[depth_cal_fusion_p[i].allocate_data[0].genome_chromosome].mutex));
      allocate_new_splice(depth_cal_fusion_p[i].allocate_data->genome_chromosome - 1, allocate_splice[s].strand_sp, 
			  allocate_splice[s].end_sp, allocate_splice[s].start_sp, allocate_splice[s].start_extend_sp, 
			  allocate_splice[s].end_extend_sp, FROM_READ, chromosome_avls_p);
      pthread_mutex_unlock(&(chromosome_avls_p[depth_cal_fusion_p[i].allocate_data[0].genome_chromosome].mutex));
    }
  }//End loop depth
  
  
  for (int i = 0; i < depth; i++) {
    free(depth_cal_fusion_p[i].allocate_data);
    free(input_p->seq_p[i]);
  }

  for (unsigned int i = 0; i < 4; i++) {
    free(output_p->mapped_seq_p[i]);
    free(output_p->mapped_ref_p[i]);
  }

  free(cigar_p);
  free(allocate_splice);

  //printf("Process Splice Junctions END\n");
  //return write_batch_p;
  
}
*/

// Select the splice junction type
inline int splice_junction_type(char nt_start_1, char nt_start_2, char nt_end_1, char nt_end_2) {
  int splice_type = NOT_SPLICE;

  if (nt_start_1 == splice_nt[G_NT]) {
    if (nt_start_2 == splice_nt[T_NT]) {
      if (nt_end_1 == splice_nt[A_NT]) {
	if (nt_end_2 == splice_nt[G_NT]) {
	  //Report GT-AG Splice
	  splice_type = GT_AG_SPLICE;
	} else if (nt_end_2 == splice_nt[T_NT]) {
	  //Report GT-AT Splice
	  splice_type = GT_AT_SPLICE;
	}
      }
    } else if (nt_start_2 == splice_nt[C_NT] && 
	       nt_end_1 == splice_nt[A_NT] && 
	       nt_end_2 == splice_nt[G_NT]) {
      //Report GC -AG
      splice_type = GC_AG_SPLICE;
    }
  } else if (nt_start_1 == splice_nt[C_NT] && 
	     nt_start_2 == splice_nt[T_NT]) {
    if (nt_end_1 == splice_nt[A_NT] && 
	nt_end_2 == splice_nt[C_NT]) {
      //Report CT-AC
      splice_type = CT_AC_SPLICE;
    }else if (nt_end_1 == splice_nt[G_NT] && 
	      nt_end_2 == splice_nt[C_NT]) {
      //Report CT-GC
      splice_type = CT_GC_SPLICE;
    } 
  } else if (nt_start_1 == splice_nt[A_NT] && 
	     nt_start_2 == splice_nt[T_NT] && 
	     nt_end_1 == splice_nt[A_NT] && 
	     nt_end_2 == splice_nt[C_NT] ) {
    //Report AT-AC
    splice_type = AT_AC_SPLICE;
  }

  return splice_type;

}


/*int genoerate_cigar_sw_output(char *seq_sw, 
			      char *ref_sw,
			      int seq_start, 
			      int ref_start,
			      size_t genome_left,
			      size_t read_left,
			      size_t genome_right,
			      size_t read_right,
			      int flank_length,
			      int len_genome_left,
			      int len_genome_right,
			      array_list_t *cigar_ops) {
  int MIN_INTRON_SIZE = 40;  
  int const MIN_GAP_SEARCH = 10;
  int map_sw_len = strlen(seq_sw);
  unsigned char automata_status = CIGAR_MATCH_MISMATCH;
  unsigned char found = NOT_SPLICE;
  int cigar_value = 0;
  int j = 0;
  int start_gap, end_gap;
  cigar_op_t *cigar_op;
  
  int tot_insertion = 0;
  int padding_left = ref_start;

  printf("genome_left=%i, read_left=%i, genome_right=%i, read_right=%i and flank_length=%i\n", genome_left, read_left, genome_right, read_right, flank_length);
  printf("Refer (%i): %s\n", seq_start, seq_sw);
  printf("Query (%i): %s\n", ref_start, ref_sw);

  while (j < map_sw_len) {
    if (ref_sw[j] != '-'  && seq_sw[j] != '-') {
      //Match/Mismatch Area	  
      if (automata_status == CIGAR_MATCH_MISMATCH) {
	cigar_value++;
	padding_left++; 
      } else {
	cigar_op = cigar_op_new(cigar_value, NULL, cigar_automata_status(automata_status));
	array_list_insert(cigar_op, cigar_ops);
	automata_status = CIGAR_MATCH_MISMATCH;
	cigar_value = 1;
      } 
    } else if (ref_sw[j] == '-' && seq_sw[j] != '-') {
      //Insertion Area
      if (automata_status == CIGAR_INSERTION) {
	cigar_value++;
	tot_insertions++;
      } else {
	cigar_op = cigar_op_new(cigar_value, NULL, cigar_automata_status(automata_status));
	array_list_insert(cigar_op, cigar_ops);
	automata_status = CIGAR_INSERTION;
	cigar_value = 1;
      }
    } else if (ref_sw[j] != '-' && seq_sw[j] == '-') {
      cigar_op = cigar_op_new(cigar_value, NULL, cigar_automata_status(automata_status));
      array_list_insert(cigar_op, cigar_ops);
      //Deletion Area. Travel in the deletions gap to found some splice junction
      start_gap = j;
      while (seq_sw[j] == '-' && j < map_sw_len) {	  
	j++;
      }
      end_gap = j - 1;
      cigar_value = end_gad - start_gap;
      //Search gap start and gap end
      if (gap_end - tot_insertion >= len_genome_left) {//¿¿¿¿ cigar_value >= MIN_GAP_SEARCH) { ????
	//Search splice junction
	found = splice_junction_type(ref_sw[start_gap], ref_sw[start_gap + 1], ref_sw[end_gap - 1], ref_sw[end_gap]);
	cigar_value = genome_left + padding_left;
	cigar_op = cigar_op_new(cigar_value, NULL, 'N');
	array_list_insert(cigar_op, cigar_ops);
      } else {
	cigar_op = cigar_op_new(cigar_value, NULL, 'D');
	array_list_insert(cigar_op, cigar_ops);
	padding_left += cigar_value;
      }
      
      if (j == map_sw_len) { return; }
      if (ref_sw[j] != '-') { automata_status = CIGAR_MATCH_MISMATCH; }
      else { automata_status = CIGAR_INSERTION; }
      cigar_value = 1;
    } else {
      //Padding Area
      if (automata_status == CIGAR_PADDING) {
	cigar_value++;
	tot_insertions++;
      } else {
	cigar_op = cigar_op_new(cigar_value, NULL, cigar_automata_status(automata_status));
	array_list_insert(cigar_op, cigar_ops);
	automata_status = CIGAR_PADDING;
	cigar_value = 1;
      }
    }
    j++;
  }

  cigar_op = cigar_op_new(cigar_value, NULL, cigar_automata_status(automata_status));
  array_list_insert(cigar_op, cigar_ops);

  //If we don't have 3 cigar operations, [OK op]--ERROR op--[OK op], we found a SW incoherence
  int number_cigar_ops = array_list_size(cigar_ops);
  if ( number_cigar_ops >= 3) {
    //Extract The first operation and the last to apply flank_length desviation
    cigar_op = array_list_get(0, cigar_ops);
    cigar_op->number -= (flank_length + 1);

    cigar_op = array_list_get(number_cigar_ops - 1, cigar_ops);
    cigar_op->number -= flank_length;
  }

}

*/

typedef struct reference_info {
  seed_region_t *first_region;
  seed_region_t *last_region;
} reference_info_t;

reference_info_t *reference_into_new(seed_region_t *first_region,
				     seed_region_t *last_region) {
  reference_info_t *reference_info = (reference_info_t *)malloc(sizeof(reference_info_t));
  reference_info->first_region = first_region;
  reference_info->last_region = last_region;

  return reference_info;
}

typedef struct array_item {
  void *data;
  int type;
} array_item_t;

array_item_t *array_item_new(void *data, int type) {
  array_item_t *new_item = (array_item_t *)malloc(sizeof(array_item_t));
  new_item->data = data;
  new_item->type = type;

  return new_item;
}

void array_item_free(array_item_t *item) {
  free(item);
}



typedef struct ref_seq_item {
  char *q;
  char *r;
  int pos;
}ref_seq_item_t;


ref_seq_item_t *ref_seq_item_new(char *q, char *r, int pos) {
  ref_seq_item_t *rf_item = (ref_seq_item_t *)malloc(sizeof(ref_seq_item_t));
  rf_item->q = q;
  rf_item->r = r;
  rf_item->pos = pos;

  return rf_item;
}


typedef struct sw_item {
  array_list_t *ref_seq_list;
  array_list_t *cigar_ops_list;
  size_t read_id;
  size_t cal_id;  
  float score;
} sw_item_t;


sw_item_t *sw_item_new(size_t read_id, size_t cal_id) {
  sw_item_t *sw_item = (sw_item_t *)malloc(sizeof(sw_item_t));
  sw_item->ref_seq_list = array_list_new(32,
					 1.5f,
					 COLLECTION_MODE_ASYNCHRONIZED);

  sw_item->cigar_ops_list = array_list_new(32,
					   1.5f,
					   COLLECTION_MODE_ASYNCHRONIZED);
  sw_item->read_id =read_id;
  sw_item->cal_id = cal_id;
  sw_item->score = 0;

  return sw_item;
}


void sw_item_free(sw_item_t *sw_item) {
  array_list_free(sw_item->ref_seq_list, NULL);
  array_list_free(sw_item->cigar_ops_list, NULL);

  free(sw_item);
}

char *generate_cigar_str_2(array_list_t *ops_list) {
  cigar_op_t *cigar_op;
  char *cigar_str = (char *)calloc(512, sizeof(char));
  char cigar_aux[128];

  for (int op = 0; op < array_list_size(ops_list); op++) {
    cigar_op = array_list_get(op, ops_list);
    sprintf(cigar_aux, "%i%c\0", cigar_op->number, cigar_op->name);
    strcat(cigar_str, cigar_aux);
  }
  
  return cigar_str;
}


int generate_cals_score(array_list_t *fusion_cals) {
  int score = 0;
  size_t num_cals = array_list_size(fusion_cals);
  cal_t *cal;
  linked_list_iterator_t itr;  
  seed_region_t *s;
  cigar_code_t *cigar_code;

  for (int j = 0; j < num_cals; j++) {
    cal = array_list_get(j, fusion_cals);
    cigar_code = (cigar_code_t *)cal->info;
    if (cigar_code) {
      score += cigar_read_coverage(cigar_code);
    } else {
      return 0;
    }
  }

  return score;

}


int apply_sw_rna(sw_server_input_t* input_p, batch_t *batch) {  
  size_t max_intron_size = 1000000;
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  sw_optarg_t *sw_optarg = &input_p->sw_optarg;
  genome_t *genome = input_p->genome_p;
  size_t num_targets = mapping_batch->num_targets;
  array_list_t *cals_list, *fusion_cals, *fusion_cals_aux;
  cal_t *cal, *cal_prev, *cal_next;
  fastq_read_t *fq_read;
  linked_list_iterator_t itr;  
  seed_region_t *s;
  cigar_code_t *cigar_code;
  cigar_op_t *cigar_op_start, *cigar_op_end;
  array_list_t *cals_targets[mapping_batch->num_allocated_targets];
  array_list_t *merge_cals;
  int *cals_score = (int *)calloc(100, sizeof(int));
  char reference[2048];
  char query[2048];

  size_t num_new_targets = 0;
  int coverage;
  int read_length;
  int cal_pos = 0;
  int number_of_best;
  int target;
  size_t genome_start, genome_end;

  register size_t num_cals;
  register size_t t;
  register int i;
  

  fill_gaps(mapping_batch, sw_optarg, genome, 20);
  merge_seed_regions(batch->mapping_batch);
  
  LOG_DEBUG_F("%s********************RNA PHASE*******************%s\n", KGRN, KWHT);

  for (t = 0; t < num_targets; t++) {
    target = mapping_batch->targets[t];
    cals_list = mapping_batch->mapping_lists[target];
    fq_read = array_list_get(target, mapping_batch->fq_batch);
    read_length = fq_read->length;
    num_cals = array_list_size(cals_list);
    
    //===== Step-1: Concatenate CALs =====//
    cals_targets[target] = array_list_new(num_cals,
				  1.25f,
				  COLLECTION_MODE_ASYNCHRONIZED);
    merge_cals = array_list_new(10,
				1.25f,
				COLLECTION_MODE_ASYNCHRONIZED);
    cal_pos = 0;
    cal_prev = (cal_t *)array_list_get(cal_pos++, cals_list);
    array_list_insert(cal_prev, merge_cals);    
    while (cal_pos < num_cals) {                                                                                                                 
      cal_next = (cal_t *)array_list_get(cal_pos, cals_list);      
      if (cal_prev->chromosome_id == cal_next->chromosome_id &&                                                                                                       
	  cal_prev->strand == cal_next->strand && 
	  (cal_next->start <= (cal_prev->end + max_intron_size))) {                                                                                                   
	array_list_insert(cal_next, merge_cals);
      } else { 
	array_list_insert(merge_cals, cals_targets[target]);
	merge_cals = array_list_new(10,
				    1.25f,
				    COLLECTION_MODE_ASYNCHRONIZED);
	array_list_insert(cal_next, merge_cals);
      }                                                                                         
      cal_prev = cal_next;
      cal_pos++;
    }
    array_list_insert(merge_cals, cals_targets[target]);         
    array_list_clear(mapping_batch->mapping_lists[target], NULL);
    //===== Step-1: END =====//


    //===== Step-2: Extend to first mismatch and Generate CALs score  =====//
    for (int i = 0; i < array_list_size(cals_targets[target]); i++) {
      fusion_cals = array_list_get(i, cals_targets[target]);
      LOG_DEBUG_F("CAL merge %i\n", i);
      for (int j = 0; j < array_list_size(fusion_cals); j++) {
	LOG_DEBUG_F("\tCAL id %i\n", j);
	cal = array_list_get(j, fusion_cals);	
	s = linked_list_get_first(cal->sr_list);
	if (s) {
	  cigar_code = (cigar_code_t *)cal->info;
	  LOG_DEBUG_F("\t\tFUSION %i:[%lu|%i - %i|%lu]: Distance(%i) %s\n", j, s->genome_start, s->read_start, s->read_end, 
		      s->genome_end, cigar_code->distance, new_cigar_code_string(cigar_code));
		  
	  if (cigar_code) {
	    cigar_op_start = cigar_code_get_first_op(cigar_code);
	    cigar_op_end = cigar_code_get_last_op(cigar_code);
	  }
	  
	  if (cigar_op_start->name == 'H') {	 
	    genome_start = s->genome_start;
	    genome_end = s->genome_start + cigar_op_start->number - 1;
	    genome_read_sequence_by_chr_index(reference, 0, 
					      cal->chromosome_id - 1, &genome_start, &genome_end, genome);
	    memcpy(query, fq_read->sequence, cigar_op_start->number);
	    query[cigar_op_start->number] = '\0';

	    LOG_DEBUG_F("\t\tExtract reference for Start hard clipping [%i:%lu-%lu]: \n", cal->chromosome_id, genome_start, genome_end);
	    LOG_DEBUG_F("REFERENCE: %s\n", reference);
	    LOG_DEBUG_F("QUERY    : %s\n", query);
	  }
		
	  if (cigar_op_end->name == 'H') {
	    genome_end = s->genome_end;
	    genome_start = s->genome_end - cigar_op_end->number + 1;
	    genome_read_sequence_by_chr_index(reference, 0, 
	    				      cal->chromosome_id - 1, &genome_start, &genome_end, genome);
	    memcpy(query, (fq_read->sequence + fq_read->length) - cigar_op_end->number , cigar_op_end->number);
	    query[cigar_op_end->number] = '\0';

	    LOG_DEBUG_F("\t\tExtract reference for End hard clipping [%i:%lu-%lu]: \n", cal->chromosome_id, genome_start, genome_end);
	    LOG_DEBUG_F("REFERENCE: %s\n", reference);
	    LOG_DEBUG_F("QUERY    : %s\n", query);
	  }
	}
      }

      cals_score[i] = generate_cals_score(fusion_cals);

    }

    //===== Step-2: END =====//


    //===== Step-3: Ranking CALs by score=====//
    /*number_of_best = 1;
    for (int i = 0; i < number_of_best; i++) {
      for (int j = i + 1; j < array_list_size(cals_targets[target]); j++) {
	//printf("%i : Compare %i > %i\n", j, cals_score[j], cals_score[i]);
	if (cals_score[j] > cals_score[i]) {
	  array_list_swap(i, j, cals_targets);
	  coverage = cals_score[j];
	  cals_score[j] = cals_score[i];
	  cals_score[i] = coverage;
	}
      }
      }*/
    //===== Step-3: END =====//


    //===== Step-4: Process CALs, fill extern gaps and mount Smith-Watermans =====//
    /*for (int i = 0; i < number_of_best; i++) {
      fusion_cals = array_list_get(i, cals_targets[target]);
      LOG_DEBUG_F("\tPROCESS CALS FUSION %i:\n", array_list_size(fusion_cals));  
      if (array_list_size(fusion_cals) == 1) {
	cal = array_list_get(j, fusion_cals);
	
      }

      for (int j = 0; j < array_list_size(fusion_cals); j++) {
	cal = array_list_get(j, fusion_cals);
	linked_list_iterator_init(cal->sr_list, &itr);
	s = linked_list_iterator_curr(&itr);
	if (s) {
	  cigar_code = (cigar_code_t *)s->info;
	  LOG_DEBUG_F("\t\tFUSION %i:[%lu|%i - %i|%lu]: Distance(%i) %s\n", j, s->genome_start, s->read_start, s->read_end, 
	  	      s->genome_end, cigar_code->distance, new_cigar_code_string(cigar_code));
	}	
      }

      }*/
    //===== Step-4: Process CALs =====//

  } //End loop targets 
  

  exit(-1);
  
  return RNA_POST_PAIR_STAGE;
}
  /*  if (s_prev && s_prev->info) {
	    if (s_prev->read_end == s->read_start - 1) {
	      cigar_op_prev = array_list_get(array_list_size(((cigar_code_t *)s_prev->info)->ops) - 1, ((cigar_code_t *)s_prev->info)->ops);
	      char name_op_prev = cigar_op_prev->name;

	      num_ops = array_list_size(cigar_code->ops);
	      op = 0;
	      cigar_op = array_list_get(0, cigar_code->ops);
	      while (cigar_op->name == name_op_prev) {
		cigar_op_prev->number += cigar_op->number;
		op++;
		if (op >= num_ops) { break; }
		cigar_op = array_list_get(op, cigar_code->ops);
	      }

	      while (op < num_ops) {
		cigar_op = array_list_get(op, cigar_code->ops);
		array_list_insert(cigar_op, ((cigar_code_t *)s_prev->info)->ops);
		op++;
	      }

	      s_prev->read_end = s->read_end;
	      s_prev->genome_end = s->genome_end;
	      LOG_DEBUG_F("\t\t\tFUSION PREV RESULT: Item [%lu|%i - %i|%lu]: \n", s_prev->genome_start, s_prev->read_start, s_prev->read_end, s_prev->genome_end);
	      for (op = 0, cigar_op = array_list_get(op, ((cigar_code_t *)s_prev->info)->ops); 
		   op < array_list_size(((cigar_code_t *)s_prev->info)->ops);
		   op++, cigar_op = array_list_get(op, ((cigar_code_t *)s_prev->info)->ops)) { 
		LOG_DEBUG_F("\t\t\t\t %i%c\n", cigar_op->number, cigar_op->name);
	      }


	  
	      continue;
	    }
	  }
	}

	s_prev = s;

      }

      LOG_DEBUG_F("Fusion Result %i\n", linked_list_size(cal->sr_list));

      linked_list_iterator_init(cal->sr_list, &itr);
      s = linked_list_iterator_curr(&itr);
      while (s) {
	LOG_DEBUG_F("\t\tItem [%lu|%i - %i|%lu]: \n", s->genome_start, s->read_start, s->read_end, s->genome_end);
	cigar_code = (cigar_code_t *)s->info;
	
	if (cigar_code) {
	  for (op = 0, cigar_op = array_list_get(op, cigar_code->ops); op < num_ops; op++, cigar_op = array_list_get(op, cigar_code->ops)) { 
	    if (cigar_op->name == 'D' && cigar_op->number > 40) { printf("Search Intron\n"); }
	    LOG_DEBUG_F("\t\t\t %i%c\n", cigar_op->number, cigar_op->name);
	  }
	} else {
	  LOG_DEBUG("\t\t\t NULL\n");
	}
	s = linked_list_iterator_next(&itr);
      }

    }
  }
  */


  ///////******************************************************************\\\\\\\\\\\\


  //printf("end\n");
  //Initialization of array_list for each CAL of the reads
  //array_list_t ***storage_ops;
  //storage_ops = (array_list_t ***)calloc(num_targets, sizeof(array_list_t **));
  /*
  array_list_t *sw_targets = array_list_new(num_targets*MAX_CALS_PROCESS*32*2,
					    1.5f,
					    COLLECTION_MODE_ASYNCHRONIZED);

  array_list_t *sw_cigar_ops = array_list_new(64,
					      1.5f,
					      COLLECTION_MODE_ASYNCHRONIZED);

  array_list_t **cals_lists = (array_list_t **)calloc(mapping_batch->num_allocated_targets, sizeof(array_list_t*));

  for (size_t i = 0; i < mapping_batch->num_allocated_targets; i++) {
    cals_lists[i] = array_list_new(50,
				   1.25f,
				   COLLECTION_MODE_ASYNCHRONIZED);
  }

  array_list_t *cals_list;
  //sw_simd_context_t *context_p = sw_simd_context_new(input_p->match, input_p->mismatch, input_p->gap_open, input_p->gap_extend);
  sw_item_t *sw_item;
  int cigar_pos;
  ref_seq_item_t *rs_item;

  //printf("I am in RNA Server Score Match %f!!\n", score_match);  
  //printf("Num Targets %i, Num Reads %i\n", num_targets, array_list_size(mapping_batch->fq_batch));

  for (size_t t = 0; t < num_targets; t++) {
    cals_list = mapping_batch->mapping_lists[mapping_batch->targets[t]];
    fq_read = array_list_get(mapping_batch->targets[t], mapping_batch->fq_batch);
    //printf("%s\n", fq_read->id);

    //TODO: Set flag negative strand 
    strcpy(read_rev, fq_read->sequence);
    seq_reverse_complementary(read_rev, fq_read->length);

    num_cals = array_list_size(cals_list);
    score_mapping = 0.0;
    tot_matches = 0;
    //printf("Total CALs %i:\n", num_cals);
    for (int i = 0; i < 1; i++) {
      sw_item = sw_item_new(mapping_batch->targets[t], i);
      cigar_pos = 0;
      cal_t *cal = array_list_get(i, cals_list);
      //printf("\tCAL%i:= Num Seeds: %i, chr %i:(%i)[%lu-%lu], Total Seeds Regions %lu: \n",i, cal->num_seeds,
      //	     cal->chromosome_id, cal->strand, cal->start, cal->end, linked_list_size(cal->sr_list));

      
      //printf("Seeds Travel ...!!\n");
      if (cal->sr_list->first) {
	seed_region_t *s = cal->sr_list->first->item;
	if (s->read_start != 0) {
	  //Seed don't starts in 0 position
	  end_read = s->read_start;
	  start_ref = s->genome_start - end_read;
	  end_ref = s->genome_start - 1;
	  
	  //printf("\tExtract First Reference from [%lu|%i] to [%i|%lu]\n", start_ref, start_read, end_read, end_ref);

	  reference = (char *)malloc(sizeof(char)*(end_read + 5));
	  genome_read_sequence_by_chr_index(reference, 0,
					    cal->chromosome_id - 1, &start_ref, &end_ref, genome);

	  query = (char *)calloc(end_read + 5, sizeof(char));

	  if (cal->strand) {
	    strncpy(query, read_rev, end_read);
	  } else {
	    strncpy(query, fq_read->sequence, end_read);
	  }

	  query[end_read] = '\0';

	  if (!strlen(query)) {
	    printf("ERROR: Query null 0\n");exit(-1);
	  }

	  num_sw++;
	  rs_item = ref_seq_item_new(query, reference, cigar_pos);
	  array_list_insert(rs_item, sw_item->ref_seq_list);

	  //array_list_insert(reference, storage_sw_data);
	  //array_list_insert(query, storage_sw_data);
	  //array_item = array_item_new(NULL, REFERENCE_LEFT);
	  //array_list_insert(array_item, cigar_ops_list);	  
	}
      }

      for (linked_list_item_t *list_item = cal->sr_list->first; list_item != NULL; list_item = list_item->next) {
	seed_region_t *s = list_item->item;

	number_op = (s->genome_end - s->genome_start);	
	sw_item->score += (number_op + 1)*score_match;
	//cigar_op = cigar_op_new(number_op + 1, NULL, 'M');
	//array_list_insert(cigar_op, cigar_ops_list);
	tot_matches += number_op + 1;
	cigar_pos++;
	//printf("\tExact +%i = %iM\n ", number_op + 1, tot_matches);	

	  
	//printf("[%i|%i - %i|%i] -", s->genome_start, s->read_start, s->read_end, s->genome_end);
	//number_op = (s->genome_end - s->genome_start);
	//if  ((s->read_end - s->read_start) == number_op) {
	  //cigar_op = cigar_op_new(number_op + 1, 'M');
	  //array_list_insert(cigar_op, cigar_ops_list);
	  //printf(" Exact %i%c - \n", cigar_op->number, cigar_op->name);
	//} else {
	  //printf(" Inexact %i/%i SW - \n", (s->read_end - s->read_start), (s->genome_end - s->genome_start));
	//}
	  

	if (list_item->next) {
	  seed_region_t *s_next = list_item->next->item;
	  if (s_next->read_start < s->read_end) { 
	    sw_item->score = 0; 
	    array_list_clear(NULL, sw_item->ref_seq_list); 
	    sw_item->ref_seq_list->size = 0;
	    break;
	  }
	  number_op = (s_next->genome_start - s->genome_end);
	  if  ((s_next->read_start - s->read_end) == number_op) {

	    score_mapping += (number_op - 1)*score_match;
	    //cigar_op = cigar_op_new(number_op - 1, NULL, 'M');
	    //array_list_insert(cigar_op , cigar_ops_list);
	    tot_matches += number_op - 1;
	    //printf("\tExact +%i = %iM\n ", number_op - 1, tot_matches);

	      //cigar_op = cigar_op_new(number_op - 1, 'M');
	    //array_list_insert(cigar_op , cigar_ops_list);
	    //printf(" Exact %i%c - \n", cigar_op->number, cigar_op->name);

	  } else {
	    cigar_op = cigar_op_new(tot_matches, NULL, 'M');
	    array_list_insert(cigar_op, sw_item->cigar_ops_list);

	    //array_item = array_item_new(cigar_op, OP_TYPE);
	    //array_list_insert(array_item, cigar_ops_list);
	    //array_list_insert(cigar_op , cigar_ops_list);

	    tot_matches = 0;
	    if (min_intron_length <= s_next->genome_start - s->genome_end) {
	      //printf("\tExtrat Intron! Seed: [%lu|%i] to [%i|%lu] - Exon Right [%lu|%i] to [%i|%lu]\n", s->genome_start, s->read_start, 
	      //     s->read_end, s->genome_end, s_next->genome_start, s_next->read_start, 
	      //     s_next->read_end, s_next->genome_end);
	      
	      start_ref = s->genome_end - flank_length;
	      end_ref = s->genome_end + flank_length;
	      reference = (char *)malloc(sizeof(char)*((end_ref - start_ref)*2 + 5));
	      genome_read_sequence_by_chr_index(reference, 0,
						cal->chromosome_id - 1, &start_ref, &end_ref, genome);
	
	      //printf("\t\tExtract [%lu-%lu]: %s\n", start_ref, end_ref, reference);

	      start_ref = s_next->genome_start - flank_length;
	      end_ref = s_next->genome_start + flank_length;  
	      genome_read_sequence_by_chr_index(ref_aux, 0,
						cal->chromosome_id - 1, &start_ref, &end_ref, genome);
	      //printf("\t\tExtract [%lu-%lu]: %s\n", start_ref, end_ref, ref_aux);

	      strcat(reference, ref_aux);
	      
	      if (s->read_end < flank_length) {
		start_read = 0;
	      }	else {
		start_read = s->read_end - flank_length;
	      }
	      end_read = s_next->read_start + flank_length;
	      //printf("Read length %i - %i = %i\n", end_read, start_read, end_read - start_read);
	      
	      query = (char *)calloc(end_read - start_read + 5, sizeof(char));
	      if (cal->strand) {
		strncpy(query, read_rev + start_read, end_read - start_read);
	      } else {
		strncpy(query, fq_read->sequence + start_read, end_read - start_read);
	      }

	      query[end_read - start_read] = '\0';

	      if (!strlen(query)) {
		printf("ERROR: Query null 1\n");exit(-1);
	      }

	      num_sw++;	  
	      rs_item = ref_seq_item_new(query, reference, cigar_pos);
	      array_list_insert(rs_item, sw_item->ref_seq_list);

	      //array_list_insert(reference, storage_sw_data);
	      //array_list_insert(query, storage_sw_data);
	      //array_item = array_item_new(NULL, REFERENCE_INTRON);
	      //array_list_insert(array_item, cigar_ops_list);	  
	    } else {
	      //printf("\tExtract Middle Reference from Inexact [%lu|%i] to [%i|%lu]\n", s->genome_end, s->read_end, 
	      //     s_next->read_start, s_next->genome_start);

	      start_ref = s->genome_end - 5;
	      end_ref = s_next->genome_start + 5;	  
	      //printf("\tExtract Reference from [%lu|%i] to [%i|%lu]\n", start_ref, start_read, end_read, end_ref);
	      reference = (char *)malloc(sizeof(char)*((end_ref - start_ref) + 15));
	      genome_read_sequence_by_chr_index(reference, 0,
						cal->chromosome_id - 1, &start_ref, &end_ref, genome);
	      start_read = s->read_end - 5;
	      end_read = s_next->read_start + 5;
	      query = (char *)calloc(((end_read - start_read) + 15), sizeof(char));
	      
	      //printf("\t\tExtract [%lu-%lu]\n", start_read, end_read);
	      if (cal->strand) {
		strncpy(query, read_rev + start_read, end_read - start_read);
	      } else {
		strncpy(query, fq_read->sequence + start_read, end_read - start_read);
	      }
	      query[end_read - start_read] = '\0';
	      //printf("\t\tExtract [%lu-%lu]: %s\n", start_read, end_read, query);
	      if (!strlen(query)) {
		printf("ERROR: Query null 2\n");exit(-1);
	      }

	      num_sw++;
	      rs_item = ref_seq_item_new(query, reference, cigar_pos);
	      array_list_insert(rs_item, sw_item->ref_seq_list);

	      //array_list_insert(reference, storage_sw_data);
	      //array_list_insert(query, storage_sw_data);
	      //array_item = array_item_new(NULL, REFERENCE_MIDLE);
	      //array_list_insert(array_item, cigar_ops_list);	  
	    }
	  }
	}
      }// End seed regions loop

      cigar_op = cigar_op_new(tot_matches, NULL, 'M');
      array_list_insert(cigar_op, sw_item->cigar_ops_list);

      //array_list_insert(cigar_op , cigar_ops_list);
      //array_item = array_item_new(cigar_op, OP_TYPE);
      //array_list_insert(array_item, cigar_ops_list);	  


      if (cal->sr_list->last) {
	seed_region_t *s = cal->sr_list->last->item;
	int number_nt = (fq_read->length - 1) - s->read_end;
	
	if (number_nt > 0) { 
	  end_read = s->read_end;
	  start_read = s->read_start;
	  start_ref = s->genome_end;
	  end_ref = s->genome_end + number_nt;
	  
	  //printf("\tExtract Final Reference from [%lu|%i] to [%i|%lu] \n", start_ref, start_read, end_read, end_ref);
	  reference = (char *)malloc(sizeof(char)*(end_ref - start_ref + 5));
	  genome_read_sequence_by_chr_index(reference, 0,
					    cal->chromosome_id - 1, &start_ref, &end_ref, genome);

	  query = (char *)calloc(number_nt + 5, sizeof(char));
	  
	  if (cal->strand) {
	    strncpy(query, read_rev + end_read, number_nt);
	  } else {
	    strncpy(query, fq_read->sequence + end_read, number_nt);
	  }
	  query[number_nt] = '\0';
	  //printf("[%i-%i]", end_read, number_nt);
	  if (!strlen(query)) {
	    printf("ERROR: Query null 3\n");exit(-1);
	  }

	  num_sw++;
	  rs_item = ref_seq_item_new(query, reference, cigar_pos);
	  array_list_insert(rs_item, sw_item->ref_seq_list);

	}
      }
    
      //printf("End of Process\n");
      array_list_insert(sw_item, sw_targets);
      linked_list_free(cal->sr_list, free);
      linked_list_free(cal->sr_duplicate_list, free);
      array_list_insert(cal, cals_lists[mapping_batch->targets[t]]);          
    } //End Loop CALs
    mapping_batch->mapping_lists[mapping_batch->targets[t]]->size = 0;
  }

  //printf("End Extract reference. SW init... %i\n", num_sw);
  char *q[num_sw];
  char *r[num_sw];  
  sw_multi_output_t *output = sw_multi_output_new(num_sw);  
  size_t sw_pos = 0;
  
  for (size_t it = 0; it < array_list_size(sw_targets); it++) {
    sw_item = array_list_get(it, sw_targets);
    fq_read = array_list_get(sw_item->read_id, mapping_batch->fq_batch);

    //printf("%s\n ", fq_red->id);

    for (size_t it2 = 0; it2 < array_list_size(sw_item->ref_seq_list); it2++) {
      rs_item = array_list_get(it2, sw_item->ref_seq_list);
      //printf("Item %i\n", sw_pos);
      q[sw_pos] = rs_item->q;
      //printf("Query: (%i)%s\n", strlen(q[sw_pos]), q[sw_pos]);
      r[sw_pos++] = rs_item->r;

      mapping_batch->histogram_sw[strlen(r[sw_pos - 1])]++;

      //printf("Refer: (%i)%s\n", strlen(r[sw_pos - 1]), r[sw_pos - 1]);
    }
  }
  
  //printf("End Sw innit %i\n", sw_pos);
  pthread_mutex_lock(&sw_mutex);
  total_sw += sw_pos;
  pthread_mutex_unlock(&sw_mutex);

  smith_waterman_mqmr(q, r, sw_pos, sw_optarg, 1, output);

  alignment_t *alignment;
  char header_id[2048];
  char *cigar_str;
  int number_ref_seq;
  int op_target;
  array_list_t *new_cigar;
  int rs_pos;
  float norm_score;

  for (size_t it = 0; it < array_list_size(sw_targets); it++) {
    sw_item = array_list_get(it, sw_targets);
    cals_list = cals_lists[sw_item->read_id];
    fq_read = array_list_get(sw_item->read_id, mapping_batch->fq_batch);
    get_to_first_blank(fq_read->id, strlen(fq_read->id), header_id);
    cal = array_list_get(sw_item->cal_id, cals_list);
    
    //Regenerate Cigar
    sw_pos = 0;
    number_ref_seq = array_list_size(sw_item->ref_seq_list);
    rs_pos = 0;
    
    if (number_ref_seq) {
      rs_item = array_list_get(rs_pos, sw_item->ref_seq_list);
      op_target = rs_item->pos;        
    
      for (int i = 0; i < array_list_size(sw_item->cigar_ops_list); i++) {
	cigar_op = array_list_get(i, sw_item->cigar_ops_list);
	if (i == op_target) {
	  if (op_target == 0) {
	    flank_length = 0;
	  } else {
	    flank_length = 5;
	  }
	  sw_item->score += output->score_p[sw_pos];
	  sw_cigar_ops->size = 0;

	  //genoerate_cigar_sw_output(output->query_map_p[sw_pos], 
	//			    output->ref_map_p[sw_pos], 
	//			    output->query_start_p[sw_pos],
	//			    output->ref_start_p[sw_pos],
	//			    0, 0, 0,
	  //			    0, flank_length, sw_cigar_ops);
	  sw_pos++;
	  rs_pos++;
	  if (rs_pos < number_ref_seq) {
	    rs_item = array_list_get(rs_pos, sw_item->ref_seq_list);
	    op_target = rs_item->pos;
	  }
	}
      }
    }

    norm_score = NORM_SCORE(sw_item->score, fq_read->length, score_match);
    if ( norm_score > 0.7) {
      //cigar_str = generate_cigar_str(sw_item->cigar_ops_list);	
	alignment = alignment_new();
	//printf("%i\n", sw_item->read_id);
	alignment_init_single_end(strdup(header_id), strdup(fq_read->sequence), strdup(fq_read->quality),
				  cal->strand,
				  cal->chromosome_id - 1, cal->start - 1,
				  strdup("100M"), 1, norm_score*254, 1, 1,
				  0, NULL, 0, alignment);
	//alignment_print(alignment);
	array_list_insert(alignment, mapping_batch->mapping_lists[sw_item->read_id]);
    }
    sw_item_free(sw_item);
  }


  for (size_t t = 0; t < mapping_batch->num_allocated_targets; t++) {
    array_list_free(cals_lists[t], cal_free);
  }

  free(cals_lists);

  array_list_free(sw_targets, NULL);
  array_list_free(sw_cigar_ops, NULL);

  //printf("End RNA Server\n");
  sw_multi_output_free(output);

  return POST_PAIR_STAGE;
  
}
  */
/*
int apply_sw_rna_1(sw_server_input_t* input_p, batch_t *batch) {
  /**----------------------------------------------------------------------------------**
   **                                    RNA-SEQ                                       **
   **----------------------------------------------------------------------------------**
  //printf("RNA Server\n");
  struct timeval start, end;
  double time;

  if (time_on) { start_timer(start); }

  mapping_batch_t *mapping_batch_p = batch->mapping_batch;  
  unsigned int seed_max_distance = input_p->seed_max_distance;
  allocate_splice_elements_t *chromosome_avls_p = input_p->chromosome_avls_p;
  unsigned int flank_length = input_p->flank_length;
  unsigned int max_intron_size = input_p->max_intron_size;

  size_t j, z, k, num_reads, read, bytes, num_cals, cals_positive, cals_negative;
  size_t i;

  // genome
  genome_t* genome_p = input_p->genome_p;
  int min_intron_length = input_p->min_intron_size;
  
  // SIMD support for Smith-Waterman
  unsigned int depth = 4, curr_depth = 0;
  sw_simd_input_t* sw_input_p = sw_simd_input_new(depth);
  sw_simd_output_t* sw_output_p = sw_simd_output_new(depth);
  sw_simd_context_t *context_p = sw_simd_context_new(input_p->match, input_p->mismatch, input_p->gap_open, input_p->gap_extend);
  float min_score = input_p->min_score;
	
  array_list_t *cals_list =  array_list_new(MAX_RNA_CALS + 32,
					    1.25f,
					    COLLECTION_MODE_ASYNCHRONIZED);

  sw_channel_t *channel_p, *sw_channels_p = (sw_channel_t*) calloc(depth, sizeof(sw_channel_t));
  size_t header_len, read_len;
  unsigned int total = 0, total_valids = 0, total_reads = 0;
  
  size_t len;
  size_t reference_fusion_len;
  size_t max_reference_fusion = 2048;
  size_t max_reference = 2048;

  char *reference = (char *)malloc(sizeof(char)*max_reference);
  if (reference == NULL) { exit(-1); }
  char *reference_fusion_p = (char *)malloc(sizeof(char)*max_reference_fusion);
  if (reference_fusion_p == NULL) { exit(-1); }

  allocate_fusion_data_t *cals_fusion_p, *allocate_cals_fusion_p = (allocate_fusion_data_t *)calloc(depth, sizeof(allocate_fusion_data_t));
  if (allocate_cals_fusion_p == NULL) { exit(-1); }

  alignment_t *alignment_p;
  cal_t *cal_prev, *cal_next, *cal, *cal_p, *cal_aux, *cal_orig, *cal_target;
  size_t  num_sw_process = 0;
  size_t sw_no_valids = 0;
  unsigned int num_cals_fusion;
  int num_mappings;
  char *read_reverse;
  char *seq_p;
  size_t num_targets;
  char *optional_fields;
  unsigned char *associate_list;
  preprocess_data_t *preprocess_data = (preprocess_data_t *)batch->optional_data;
  unsigned char *negative_strand = preprocess_data->negative_strand;
  int num_cal_targets;
  int *cal_targets;
  size_t c_orig, c_target;
  unsigned char no_store;
  size_t flank;

  curr_depth = 0;
  num_targets = mapping_batch_p->num_targets;//sw_batch_p->num_reads;
  num_reads = array_list_size(mapping_batch_p->fq_batch);
  total_reads += num_targets;
    
  /** for each read **
  for (i = 0; i < num_targets; i++) {      
    num_cal_targets = preprocess_data->num_cal_targets[i];
    cal_targets = preprocess_data->cal_targets[i];
    associate_list = preprocess_data->associate_cals[i];
    //cals_list = mapping_batch_p->mapping_lists[mapping_batch_p->targets[i]];

    fastq_read_t *read_p = array_list_get(mapping_batch_p->targets[i], mapping_batch_p->fq_batch); 
    //rintf("%s\n", read_p->sequence);
    //read_len = read_p->length;
    read_len = strlen(read_p->sequence);

    header_len = strlen(read_p->id);
    num_cals = array_list_size(mapping_batch_p->mapping_lists[mapping_batch_p->targets[i]]);
    if (!num_cals) { continue; }
    //printf("Process %lu CALs\n", num_cals);
    //printf("RNA SERVER :: %d - %d\n", read_len, header_len);
    total += num_cals;
    //printf("Process read(%i): %s\n", i, read_p->sequence);
    //printf("Need Reverse? (%i)\n", negative_strand[i]);
    // Order CALs and count number of CALs per chromosome and initialize cal type array
    if (negative_strand[i]) {
      read_reverse = (char *)calloc((read_len + 1), sizeof(char));
      memcpy(read_reverse, read_p->sequence, read_len);
      seq_reverse_complementary(read_reverse, read_len);
    }
    
    //printf("NUM CALS %i\n", num_cals);
    for (j = 0; j < num_cals; j++) {
      cal = (cal_t *)array_list_get(j, mapping_batch_p->mapping_lists[mapping_batch_p->targets[i]]);
      //printf("%i - %i, chr %i:(%i)[%lu-%lu]\n", j, i, cal->chromosome_id, cal->strand, cal->start, cal->end);
      array_list_insert(cal, cals_list);
    }
    mapping_batch_p->mapping_lists[mapping_batch_p->targets[i]]->size = 0;    

   
    for (int t = 0; t < num_cal_targets; t++) {
      j = cal_targets[t];
      
      cals_fusion_p = &allocate_cals_fusion_p[curr_depth];
      cals_fusion_p->allocate_data = (cal_fusion_data_t *)calloc(num_cals, sizeof(cal_fusion_data_t));      
      num_cals_fusion = 0;
      
      c_orig = j;
      c_target = associate_list[j];
      //printf("\tORGI: %i to TARGET: %i\n", c_orig, c_target);
      cal_p = (cal_t *)array_list_get(j, cals_list);
      
      cal_target = (cal_t *)array_list_get(c_target, cals_list);
      cal_orig = (cal_t *)array_list_get(c_orig, cals_list);

      if ( (cal_p->flank_start + 1) >= cal_p->start) {
	cal_p->start = 0;
      } else {
	cal_p->start -= (cal_p->flank_start + 1); 
      }

      if (c_orig == c_target) {
	flank = cal_p->flank_end + 1; 
      }	else {
	flank = flank_length; 
      }
      
      cal_p->end += flank; 
      if (cal_p->end >= genome_p->chr_size[cal_p->chromosome_id - 1]) {
	cal_p->end = genome_p->chr_size[cal_p->chromosome_id - 1] - 1;
      }
      
      
      //printf("\tNew CAL Pos[%i-%i]\n", cal_p->start, cal_p->end);
      reference_fusion_len = cal_p->end - cal_p->start + 1;
      
      if(reference_fusion_len >= max_reference_fusion) {
	max_reference_fusion += 1.25*reference_fusion_len;
	free(reference_fusion_p);
	reference_fusion_p = (char *)malloc(max_reference_fusion*sizeof(char));
	if (reference_fusion_p == NULL) { exit(-1); }
      }
      
      genome_read_sequence_by_chr_index(reference_fusion_p, 0,
					cal_p->chromosome_id - 1, &cal_p->start, &cal_p->end, genome_p);
      
      cal_fusion_data_init(j, cal_p->start, cal_p->end, cal_p->strand, cal_p->chromosome_id, 
			   0, reference_fusion_len - 1, 
			   &cals_fusion_p->allocate_data[num_cals_fusion]);

      no_store = 0;
      while (j < c_target) {
	j++;
	cal_p = (cal_t *)array_list_get(j, cals_list);
	num_cals_fusion++;
	//printf("\t<--- Extend left %i\n", flank_length);
	if ( flank_length >= cal_p->start) {
	  cal_p->start = 0;
	} else {
	  cal_p->start -= flank_length; 
	}
      
	if (j == c_target) {
	  flank = cal_p->flank_end + 1; 
	} else {
	  flank = flank_length; 
	}

	cal_p->end += flank; 
	if (cal_p->end >= genome_p->chr_size[cal_p->chromosome_id - 1]) {
	  cal_p->end = genome_p->chr_size[cal_p->chromosome_id - 1] - 1;
	}
	
	len = cal_p->end - cal_p->start + 1;
	if(reference_fusion_len + len >= MAX_FUSION) {
	  no_store = 1;
	  break;
	}
	
	//printf("%lu >= %lu\n", len, max_reference);
	if(len >= max_reference) {
	  max_reference += 1.25*len;
	  free(reference);
	  reference = (char *)malloc(max_reference*sizeof(char));
	  if (reference == NULL) { exit(-1); }
	}
	
	//printf("\tCAL fusion data fusion chr:%i-strand%i-[%i-%i] %i...\n", cal_p->chromosome_id, cal_p->strand, cal_p->start, cal_p->end, j);
	
	genome_read_sequence_by_chr_index(reference, 0, 
					  cal_p->chromosome_id - 1, &cal_p->start, &cal_p->end, genome_p);

	cal_fusion_data_init(j, cal_p->start, cal_p->end, cal_p->strand, cal_p->chromosome_id, 
			     reference_fusion_len, reference_fusion_len + len - 1, 
			     &cals_fusion_p->allocate_data[num_cals_fusion]);
	
	reference_fusion_len += len;
	strcat(reference_fusion_p, reference);
      } //End While concatenate
      
      if (no_store){
	no_store = 0;
	while (j <= c_target) { j++; }
	//exit(-1);
	continue;
      }
      
      //printf("\tStoring reference...\n");
      cals_fusion_p->cal_number = num_cals_fusion + 1;   	  
      channel_p = &sw_channels_p[curr_depth];
      
      sw_channel_allocate_ref(reference_fusion_len + 1, channel_p);
      
      memcpy(channel_p->ref_p, reference_fusion_p, reference_fusion_len);
      channel_p->ref_p[reference_fusion_len]= '\0';
      
      sw_channel_update(mapping_batch_p->targets[i], cal_p->strand, read_len, header_len, reference_fusion_len, channel_p);
      seq_p = (char *)calloc((read_len + 1), sizeof(char));
      if(cal_p->strand){
	memcpy(seq_p, read_reverse, read_len);
      }else {
	memcpy(seq_p, read_p->sequence, read_len);
      }

      sw_simd_input_add(seq_p, read_len,
			channel_p->ref_p, channel_p->ref_len, 
			curr_depth, sw_input_p);

      if ((++curr_depth) == depth) {
	smith_waterman_simd(sw_input_p, sw_output_p, context_p);
	//printf("\tProcess sw done\n");
	search_splice_junctions_sw_output(sw_input_p, sw_output_p, curr_depth, 
					  allocate_cals_fusion_p, chromosome_avls_p, sw_channels_p, 
					  mapping_batch_p, 0, &sw_no_valids,
					  min_score, genome_p, min_intron_length, input_p->gap_open, 
					  input_p->gap_extend, input_p->match);
	num_sw_process += curr_depth;
	curr_depth = 0;
      }
      //printf("\tStored reference done!\n");
      //printf("%i\n", associate_list[i]);
    } //End Loop CALs
    //printf("End Loop\n");
    
    if (negative_strand[i]) { free(read_reverse); } 
    //printf("Clear list...\n");
    array_list_clear(cals_list, cal_free);
    //printf("Clear done...\n");
  }// end of for 0..num_reads
  //printf("End process reads\n");

  //printf("%d.Seq(%d): %s\n",cal_p->strand, read_len, sw_batch_p->allocate_reads_p[i]->sequence);
  if (curr_depth > 0) {
    //printf("Current depth %d/%d\n", curr_depth, depth);
    num_sw_process += curr_depth;
    for(k = curr_depth ; k < depth ; k++) {
      //printf("STORE : %s\n", sw_input_p->seq_p[0]);
      sw_simd_input_add(sw_input_p->seq_p[0], sw_input_p->seq_len_p[0],
			channel_p->ref_p, channel_p->ref_len, 
			k, sw_input_p);
      //printf("k=%d : %s\n", k , sw_input_p->seq_p[k]);
    }
    //printf("Run smith waterman... %i\n", i);
    smith_waterman_simd(sw_input_p, sw_output_p, context_p);
      
    search_splice_junctions_sw_output(sw_input_p, sw_output_p, curr_depth, 
				      allocate_cals_fusion_p, chromosome_avls_p,  
				      sw_channels_p, mapping_batch_p,
				      0, &sw_no_valids, min_score, genome_p, 
				      min_intron_length, input_p->gap_open, 
				      input_p->gap_extend, input_p->match);
    
    //found_write_p = process_sw_output(sw_output_p, sw_input_p, min_score, curr_depth, sw_channels_p, sw_batch_p, write_list_p, found_write_p, write_size, sw_id, &total_valids, mapping_reads_p, genome_p);
    curr_depth = 0;
  }
  
  //**************************************************
    //PROCESS OUTPUT ALIGNMENTS
    for (i = 0; i < num_targets; i++) {
      //num_mappings = array_list_size(mapping_batch_p->mapping_lists[mapping_batch_p->targets[i]]);
      /*num_mappings = alignments_filter(input_p->bwt_optarg_p->report_all, 
				       input_p->bwt_optarg_p->report_best, 
				       input_p->bwt_optarg_p->report_n_hits, 
				       mapping_batch_p->mapping_lists[mapping_batch_p->targets[i]]); 
      *
      num_mappings = array_list_size(mapping_batch_p->mapping_lists[mapping_batch_p->targets[i]]); 

      if ( num_mappings > 0 ) {
	array_list_set_flag(1, mapping_batch_p->mapping_lists[mapping_batch_p->targets[i]]);
	for (int m = 0; m < num_mappings; m++) {
	  alignment_p = array_list_get(m, mapping_batch_p->mapping_lists[mapping_batch_p->targets[i]]);
	  //alignment_print(alignment_p);
	  optional_fields = (char *)alignment_p->optional_fields;
	  optional_fields += alignment_p->optional_fields_length;
	  sprintf(optional_fields, "NHi");
	  alignment_p->optional_fields_length += 3;
	  optional_fields += 3;
	  memcpy(optional_fields, &num_mappings, sizeof(int));
	  alignment_p->optional_fields_length += sizeof(int);
	  //alignment_p->optional_fields = (uint8_t *)optional_fields;
	  //printf("::%s\n", optional_fields);
	}
      } else {
	array_list_set_flag(0, mapping_batch_p->mapping_lists[mapping_batch_p->targets[i]]);
      }
    }
    

    //printf("RNA Server process batch finish!\n");

  
    array_list_free(cals_list, NULL);

  // insert or free memory
  /*  if (write_batch_p != NULL) {
    if (write_batch_p->size > 0) {
      item_p = list_item_new(0, WRITE_ITEM, write_batch_p);
      list_insert_item(item_p, write_list_p);
    } else {
      write_batch_free(write_batch_p);
    }
    }*

  for(k=0 ; k<depth ; k++) {
    free(sw_channels_p[k].ref_p);
  }

  
  free(sw_channels_p);
  sw_simd_input_free(sw_input_p);
  sw_simd_output_free(sw_output_p);
  sw_simd_context_free(context_p);

  free(allocate_cals_fusion_p);
  //free(mapping_reads_p);
  
  free(reference_fusion_p);
  free(reference);

  preprocess_data_free(preprocess_data);

  if (time_on) { stop_timer(start, end, time); timing_add(time, RNA_SERVER, timing); }

  //printf("RNA End\n");  
  return POST_PAIR_STAGE;

}

*/
