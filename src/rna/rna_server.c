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

#define FILL_GAP_LEFT  0
#define FILL_GAP_RIGHT 1

#define CALING_DOUBLE_ANCHORS  0
#define CALING_LEFT_ANCHORS    1
#define CALING_RIGHT_ANCHORS   2

#define META_OPEN  0
#define META_CLOSE 1

#define META_ALIGNMENT_NONE   0
#define META_ALIGNMENT_MIDDLE 1
#define META_ALIGNMENT_LEFT   2
#define META_ALIGNMENT_RIGHT  3

#define CIGAR_SW_MIDDLE      0
#define CIGAR_SIMPLE_MIDDLE  1
#define CIGAR_ANCHOR_LEFT    2
#define CIGAR_ANCHOR_RIGHT   3

#define MAX_DEPTH 4

#define SIMPLE_SW       1
#define SP_SW           2
#define EXTREM_SW_LEFT  4
#define EXTREM_SW_RIGHT 5

#define SP_METAEXON     3

#define SW_NORMAL 0
#define SW_FINAL 1

#define OP_TYPE 0
#define REFERENCE_LEFT 1
#define REFERENCE_MIDLE 2
#define REFERENCE_INTRON 3
#define REFERENCE_RIGHT 4

#define MINIMUN_CAL_LENGTH 10
#define ERRORS_ZONE 8
#define MINIMUM_INTRON_LENGTH 10
#define MAX_CIGAR_OPERATIONS 20

#define LEFT_EXTREM 0
#define RIGHT_EXTREM 1

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
//#define NOT_SPLICE	        0
//#define GT_AG_SPLICE  	1
//#define CT_AC_SPLICE  	2
//==========================================

#define SEARCH_START 0
#define SEARCH_END 1
#define POSSIBLES_MARKS 2
#define STRANDS_NUMBER 2
#define MAX_FUSION 2028
#define MIN_HARD_CLIPPING 10


#define SW_LEFT   0
#define SW_RIGHT  1
#define SW_MIDDLE 2

int ii = -1;


cigar_code_t *search_left_single_anchor(int gap_close, 
					cal_t *cal,
					int filter_pos, 
					array_list_t *right_breaks,
					char *query_map,
					metaexons_t *metaexons,
					genome_t *genome);

cigar_code_t *search_right_single_anchor(int gap_close, 
					 cal_t *cal,
					 int filter_pos, 
					 array_list_t *left_breaks,
					 char *query_map, metaexons_t *metaexons,
					 genome_t *genome);

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

// Select the splice junction type
int splice_junction_type(char nt_start_1, char nt_start_2, char nt_end_1, char nt_end_2) {
  LOG_DEBUG_F("SEARCH SPLICE JUNCTION TYPE FOR %c%c - %c%c\n", nt_start_1, nt_start_2, nt_end_1, nt_end_2);

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


cigar_code_t *generate_cigar_sw_output(char *seq_sw, 
				       char *ref_sw,
				       size_t l_exon_start,
				       size_t l_exon_end,
				       size_t r_exon_start,
				       size_t r_exon_end,
				       int chromosome,
				       int strand,
				       int seq_start, 
				       int ref_start,
				       int len_orig_seq,
				       int len_orig_ref,
				       avls_list_t *avls_list,
				       avl_node_t **node_avl_start,
				       avl_node_t **node_avl_end) {

  *node_avl_start = NULL;
  *node_avl_end   = NULL;

  int MIN_INTRON_SIZE = 40;  
  int const MIN_GAP_SEARCH = 15;
  int const EXTRA_SEARCH = 8;
  int map_sw_len = strlen(seq_sw);
  unsigned char automata_status = CIGAR_MATCH_MISMATCH;
  int found = NOT_SPLICE;
  int cigar_value = 0;
  int j = 0;
  int start_gap, end_gap;
  cigar_op_t *cigar_op;  
  int insertions_tot = 0;
  int deletions_tot = 0;
  int tot_matches = 0;
  int padding_left = ref_start;
  //int mode = MODE_EXON;
  int len_ref, len_r_gap;
  int gap_len;
  int cnt_ext = 0;
  int pos;
  int sw_seq_len = 0;
  
  size_t start_splice, end_splice;
  cigar_code_t *cigar_code = cigar_code_new();
  cigar_op_t *op;
  
  int n_splice = 0;
  
  if (seq_start > 0) {
    //Middle or last ref
    if (ref_start == 0) {
      cigar_code_append_op(cigar_op_new(seq_start, 'I'), cigar_code);
    } else {
      if (ref_start == seq_start) {
	cigar_code_append_op(cigar_op_new(seq_start, 'M'), cigar_code);
      } else {
	if (ref_start > seq_start) {
	  cigar_code_append_op(cigar_op_new(ref_start - seq_start, 'D'), cigar_code);
	  cigar_code_append_op(cigar_op_new(seq_start, 'M'), cigar_code);
	  cigar_code->distance += seq_start;
	} else {
	  cigar_code_append_op(cigar_op_new(seq_start - ref_start, 'I'), cigar_code);
	  cigar_code_append_op(cigar_op_new(ref_start, 'M'), cigar_code);
	  cigar_code->distance += ref_start;
	} 
      }
    }
  } else if (ref_start > 0) {
    cigar_code_append_op(cigar_op_new(ref_start, 'D'), cigar_code);
  }

  while (j < map_sw_len) {
    //printf("NT-SEQ %c (%i)\n", seq_sw[j], cigar_value);
    if (ref_sw[j] != '-'  && seq_sw[j] != '-') {
      tot_matches++;
      padding_left++;
      if (ref_sw[j] != seq_sw[j]) {
	cigar_code->distance++;
      }
      //Match/Mismatch Area	  
      if (automata_status == CIGAR_MATCH_MISMATCH) {
	cigar_value++;
      } else {
	//printf("Report %i\n", cigar_value);
	cigar_code_append_new_op(cigar_value, cigar_automata_status(automata_status), cigar_code);
	automata_status = CIGAR_MATCH_MISMATCH;
	cigar_value = 1;
      } 
    } else if (ref_sw[j] == '-' && seq_sw[j] != '-') {
      insertions_tot++;
      //Insertion Area
      if (automata_status == CIGAR_INSERTION) {
	cigar_value++;
      } else {
	//printf("Report %i\n", cigar_value);
	cigar_code_append_new_op(cigar_value, cigar_automata_status(automata_status), cigar_code);
	automata_status = CIGAR_INSERTION;
	cigar_value = 1;
      }
    } else if (ref_sw[j] != '-' && seq_sw[j] == '-') {
      //printf("Report %i\n", cigar_value);
      op = cigar_op_new(cigar_value, cigar_automata_status(automata_status));
      cigar_code_append_op(op, cigar_code);
      //Deletion Area. Travel in the deletions gap to found some splice junction
      start_gap = j;
      while (seq_sw[j] == '-' && j < map_sw_len) {	  
	j++;
	deletions_tot++;
      }
      end_gap = j - 1;
      gap_len = end_gap - start_gap + 1;
      //Search gap start and gap end
      if (gap_len > MIN_GAP_SEARCH) {
	//Search splice junction
	found = splice_junction_type(ref_sw[start_gap], ref_sw[start_gap + 1], ref_sw[end_gap - 1], ref_sw[end_gap]);
	//printf("Found %i == %i\n", found, NOT_SPLICE);

	if (found == NOT_SPLICE) {
	  //Search Xnt (+)---->	
	  cnt_ext = 1;
	  while (cnt_ext < EXTRA_SEARCH ) {	  
	    found = splice_junction_type(ref_sw[start_gap + cnt_ext], ref_sw[start_gap + cnt_ext + 1], 
					 ref_sw[end_gap + cnt_ext - 1], ref_sw[end_gap + cnt_ext]);	       
	    if (found != NOT_SPLICE) {
	      //printf("Found %i\n", found);
	      break;
	    } else {
	      //printf("continue search...\n");
	    }
	    cnt_ext++;
	  }
	}
	
	//if (!found) { printf("Seq: %s\n", seq_sw); printf("Ref: %s\n", ref_sw); }

	if (found == NOT_SPLICE) {
	  cnt_ext = 0;
	}

	len_ref = l_exon_end - l_exon_start + 1;
	len_r_gap = gap_len - (len_ref - (tot_matches + ref_start));
	if (len_r_gap < 0) { assert(len_r_gap); }
	LOG_DEBUG_F("Calculating SP: l_exon_end = %lu, l_exon_start = %lu, gap_len = %i, len_ref = %i, tot_matches = %i, len_r_gap=%i, padding_left = %i, cnt_ext = %i, seq_start = %i, ref_start = %i, r_exon_start = %i\n",
		    l_exon_end, l_exon_start, gap_len, len_ref, tot_matches, len_r_gap, padding_left, cnt_ext, seq_start, ref_start, r_exon_start);
	
	start_splice = l_exon_start + padding_left + cnt_ext;
	end_splice = r_exon_start + len_r_gap - 1 + cnt_ext;
	cigar_value = end_splice - start_splice;
	if (cigar_value < 0) { assert(cigar_value); }
	
	cigar_code_append_new_op(cigar_value + 1, 'N', cigar_code);	
	if (found == CT_AC_SPLICE || found == GT_AT_SPLICE || found == CT_GC_SPLICE ) {
	  strand = 1;
	} else if (found != NOT_SPLICE) { 
	  strand = 0;
	}

	//printf("SP COORDS(%i)=> [%i:%lu-%lu] = %d\n", strand, chromosome, start_splice, end_splice, cigar_value);

	//if (chromosome == 1 && start_splice == 17743 && end_splice == 17914) { printf("%s ::-->\n", id); exit(-1); }
	n_splice++;
	if (n_splice > 1) { 
	  array_list_clear(cigar_code->ops, cigar_op_free);
	  cigar_code_free(cigar_code); 
	  return NULL; 
	}

	allocate_start_node(chromosome - 1,
			    strand,
			    start_splice,
			    end_splice,
			    start_splice,
			    end_splice,
			    FROM_READ, found,
			    NULL, 
			    node_avl_start,
			    node_avl_end, 
			    avls_list);

	//printf("End allocate %i -> %i\n", op->number, op->number + cnt_ext);
	op->number += cnt_ext;
	  //} else {
	  //printf(":( NOT FOUND! [%lu-%lu]---[%lu-%lu]\n", l_exon_start, l_exon_end, r_exon_start, r_exon_end);
	  //cigar_value = gap_len;	
	  //cigar_code_append_new_op(cigar_value, 'D', cigar_code);	
	  //padding_left += cigar_value;
	  //}
      } else {
	cigar_value = gap_len;	
	//printf("Report %i\n", cigar_value);
	cigar_code_append_new_op(cigar_value, 'D', cigar_code);	
	padding_left += cigar_value;
      }
      
      if (j == map_sw_len) { 
	array_list_clear(cigar_code->ops, cigar_op_free);
	cigar_code_free(cigar_code);
	return NULL; 
      }
      if (ref_sw[j] != '-') { automata_status = CIGAR_MATCH_MISMATCH; tot_matches++; }
      else { automata_status = CIGAR_INSERTION; insertions_tot++; }
      cigar_value = 1 - cnt_ext;
      cnt_ext = 0;
    } else {
      //Padding Area
      insertions_tot++;
      deletions_tot++;
      if (automata_status == CIGAR_PADDING) {
	cigar_value++;
      } else {
	//printf("Report %i\n", cigar_value);
	cigar_code_append_new_op(cigar_value, cigar_automata_status(automata_status), cigar_code);
	automata_status = CIGAR_PADDING;
	cigar_value = 1;
      }
    }
    j++;
  }

  cigar_code_append_new_op(cigar_value, cigar_automata_status(automata_status), cigar_code);

  //sw_seq_len = tot_matches + tot_insertions + seq_start;

  int map_seq_len  = ((map_sw_len - deletions_tot) + seq_start);
  int map_ref_len  = ((map_sw_len - insertions_tot) + ref_start);
  int last_h, last_h_aux;

  //printf("query_start = %i, ref_start = %i, map_seq_len = %i, map_ref_len = %i, query_len = %i, ref_len = %i, map_len = %i\n", 
  //	 query_start, ref_start, map_seq_len, map_ref_len, query_len, ref_len, map_len);
  //printf("***Map_seq_len = %i, map_Ref_len = %i, len_orig_seq = %i, len_orig_ref = %i\n", 
  //	 map_seq_len, map_ref_len, len_orig_seq, len_orig_ref);

  if (map_seq_len < len_orig_seq) {
    last_h = len_orig_seq - map_seq_len;
    //printf("last_h = %i\n", last_h);
    //Middle or first ref
    if (map_ref_len == len_orig_ref) {
      cigar_code_append_op(cigar_op_new(last_h, 'I'), cigar_code);
    } else {
      last_h_aux = len_orig_ref - map_ref_len;
      //printf("last_h_aux = %i\n", last_h_aux);
      if (last_h_aux == last_h) {
	cigar_code_append_op(cigar_op_new(last_h, 'M'), cigar_code);
      } else {	  
	if (last_h_aux > last_h) {
	  cigar_code_append_op(cigar_op_new(last_h_aux - last_h, 'D'), cigar_code);
	  cigar_code_append_op(cigar_op_new(last_h, 'M'), cigar_code);
	  cigar_code->distance += last_h;
	} else {
	  cigar_code_append_op(cigar_op_new(last_h - last_h_aux, 'I'), cigar_code);
	  cigar_code_append_op(cigar_op_new(last_h_aux, 'M'), cigar_code);
	  cigar_code->distance += last_h_aux;
	} 
      }
    }  
  } else if (map_ref_len < len_orig_ref) {
    cigar_code_append_op(cigar_op_new(len_orig_ref - map_ref_len, 'D'), cigar_code);
  }

  if (n_splice == 0) {
    array_list_clear(cigar_code->ops, cigar_op_free);
    cigar_code_free(cigar_code);
    return NULL;
  }

  //Refresh distance cigar
  /* for (int i = 0; i < cigar_code_get_num_ops(cigar_code); i++) {
    cigar_op_t *cigar_op = array_list_get(i, cigar_code->ops);
    if (cigar_op->name == 'D' || cigar_op->name == 'I') {
      cigar_code->distance += cigar_op->number;
    }
  }
  */

  return cigar_code;

}

typedef struct fusion_coords {
  size_t l_exon_start;
  size_t l_exon_end;
  size_t r_exon_start;
  size_t r_exon_end;
  size_t read_start;
  size_t read_end;
  int chromosome;
  int strand;
  int type_sw;
  size_t read_id;
  char *id;
  cal_t *cal_ref;
  int cal_group_id;
} fusion_coords_t;

fusion_coords_t *fusion_coords_new(size_t l_exon_start,
				   size_t l_exon_end,
				   size_t r_exon_start,
				   size_t r_exon_end,
				   size_t read_start,
				   size_t read_end, 
				   int chromosome,
				   int strand,
				   int type_sw,
				   char *id, 
				   cal_t *cal_ref,
				   size_t read_id,
				   int cal_group_id) {

  fusion_coords_t *fusion_coords = (fusion_coords_t *)malloc(sizeof(fusion_coords_t));

  fusion_coords->l_exon_start = l_exon_start;
  fusion_coords->l_exon_end   = l_exon_end;
  fusion_coords->r_exon_start = r_exon_start;
  fusion_coords->r_exon_end   = r_exon_end;
  fusion_coords->read_start   = read_start;
  fusion_coords->read_end     = read_end;
  fusion_coords->chromosome   = chromosome;
  fusion_coords->strand       = strand;
  fusion_coords->type_sw      = type_sw;
  fusion_coords->id           = id;
  fusion_coords->cal_ref      = cal_ref;
  fusion_coords->read_id      = read_id;
  fusion_coords->cal_group_id = cal_group_id;

  return fusion_coords;
}

void fusion_coords_free(fusion_coords_t *fusion_coords) {
  free(fusion_coords);
}

    //array_list_size(cals_targets[target])
    /*for (i = 0; i < number_of_best; i++) {
      fusion_cals = array_list_get(i, cals_targets[target]);
      printf("\t FUSION PACKAGE %i(%f):\n", i, cals_score[i]);
      for (j = 0; j < array_list_size(fusion_cals); j++) {
	cal = array_list_get(j, fusion_cals);
	printf("\t\t CAL(%i) [%i:%lu-%lu]:\n", j, cal->chromosome_id, cal->start,
	       cal->end);
	
	linked_list_iterator_init(cal->sr_list, &itr);
	s = (seed_region_t *) linked_list_iterator_curr(&itr);
	while (s != NULL) {
	  printf("\t\t\t seed: [%i|%i - %i|%i]\n",
		 s->genome_start, s->read_start, s->read_end, s->genome_end);
	  s = (seed_region_t *) linked_list_iterator_next(&itr);
	}
      } 
      }*/

void extend_by_mismatches(char *ref_e1, char *ref_e2, char *query, 
			  int r_e1, int r_e2, 
			  int q_e1, int q_e2, int lim_err,
			  int *dsp_e1, int *dsp_e2) {
  int len_ref_e1 = strlen(ref_e1);
  int len_ref_e2 = strlen(ref_e2);
  int len_query = strlen(query);
  int num_err = 0;
  int tot_err = 0;
  *dsp_e1 = 0;
  *dsp_e2 = 0;
  //printf("e1=%i, %i, %i, %c == %c, %i, %i, %i, %i, %i\n", e1, strlen(query), 
  //strlen(reference_prev), reference_prev[e1], query[e1], genome_start, genome_end, read_start, read_end, seeds_nt);
  //printf("Ref Ex1: %s\n", &ref_e1[0]);
  //printf("Ref Ex2: %s\n", &ref_e2[0]);
  //printf("Query  : %s\n", &query[0]);
 
  while (q_e1 < len_query && r_e1 < len_ref_e1) {
    //printf("START: %c != %c (%i)\n", ref_e1[r_e1], query[q_e1], num_err);
    if (ref_e1[r_e1] != query[q_e1]) {
      num_err++;
      if (num_err >= lim_err) { break; }
    }
    r_e1++;
    q_e1++;
    (*dsp_e1)++;
  }

  tot_err += num_err;
  num_err = 0;

  while (r_e2 >= 0 && q_e2 >= 0) { 
    //printf("END: %c != %c(%i)\n", ref_e2[r_e2], query[q_e2], num_err);
    if (ref_e2[r_e2] != query[q_e2]) {
      num_err++;
      if (num_err >= lim_err) { break; }
    }
    r_e2--;
    q_e2--;
    (*dsp_e2)++;
  }

  tot_err += num_err;

}


cigar_code_t *fill_extrem_gap(char *query, 
			      cal_t *cal,
			      int type,
			      genome_t *genome,
			      metaexon_t *metaexon,
			      metaexons_t *metaexons) {
  
  //printf("FILL EXTREM GAPS...%s\n", type == FILL_GAP_LEFT? "LEFT" : "RIGHT");
  int chromosome_id = cal->chromosome_id;
  size_t genome_start, genome_end;
  int read_start, read_end;
  int read_gap;
  char reference[2048];
  int type_search;
  cigar_code_t *cigar_code = NULL;
  int length = strlen(query);

  //FILL_GAP_LEFT  0
  //FILL_GAP_RIGHT 1
  if (type == FILL_GAP_LEFT) {
    //printf("FILL_GAP_LEFT\n");
    seed_region_t *s_prev = linked_list_get_first(cal->sr_list);
    read_start = 0;
    read_end = s_prev->read_start - 1;
    read_gap = s_prev->read_start;
    genome_start = s_prev->genome_start - read_gap + 1;
    genome_end = s_prev->genome_start;
  } else {
    seed_region_t *s_prev = linked_list_get_last(cal->sr_list);
    read_start = s_prev->read_end + 1;
    read_end = length - 1;
    read_gap = length - s_prev->read_end - 1;
    genome_start = s_prev->genome_end + 1;
    genome_end = s_prev->genome_end + read_gap + 1;
    //printf("FILL_GAP_RIGHT %i\n", read_gap);
  }

  if (metaexon != NULL) {
    //printf("METAEXON NOT NULL!\n");
    if (type == FILL_GAP_LEFT) {
      //printf("genome_start(%i) >= metaexon->start(%i)\n", genome_start, metaexon->start);
      if (genome_start >= metaexon->start) {
	type_search = 1;
      } else { type_search = 2; }
    } else {
      //printf("genome_end(%i) >= metaexon->end(%i)\n", genome_end, metaexon->end);
      if (genome_end <= metaexon->end) {
	type_search = 1;
      } else { type_search = 3; }
    }    
    if (type_search == 2 &&
	!metaexon->left_closed) {
      type_search = 0;
    } else if (type_search == 3 &&
	       !metaexon->right_closed) {
      type_search = 0;
    }
  } else { type_search = 0; }
  
  
  if (type_search == 2) {
    //Search other
    //printf("Search WITH RIGHT ANCHOR\n");
    cigar_code = search_right_single_anchor(read_gap, 
					    cal,
					    0, 
					    metaexon->left_breaks,
					    query,
					    metaexons,
					    genome);    
  } else if (type_search == 3) {
    //printf("Search WITH LEFT ANCHOR\n");
    cigar_code = search_left_single_anchor(read_gap, 
					   cal,
					   0, 
					   metaexon->right_breaks,
					   query,
					   metaexons,
					   genome);   
  } else {
    genome_read_sequence_by_chr_index(reference, 0, 
				      chromosome_id - 1,
				      &genome_start, &genome_end,
				      genome);    
    int distance = 0;
    int err_dsp = 0;
    const int num_err = 2;
    
    //printf("REF: %s\n", reference);
    //printf("FILL EXTREM GAP %i\n", read_gap );
    for (int i = 0; i < read_gap; i++) {
      //printf("%c vs %c [%i/%i] (%i)\n", query[read_start + i], reference[i], 
      // i, read_gap, distance);      
      if (query[read_start + i] != reference[i]) {
	distance++;
      }
    }
    
    //printf("distance = %i < max_err = %i\n", distance, ((read_gap / 3) + 1));

    const int MAX_GAP = 10;
    const int MAX_ERR = 3;
   
    if ((read_gap <= MAX_GAP && distance <= MAX_ERR) || 
	(read_gap > MAX_GAP && distance <= ((read_gap / 4) + 1))) {
      cigar_code = cigar_code_new();
      cigar_code->distance = distance;
      cigar_code_append_new_op(read_gap, 'M', cigar_code);
    }
    
    if (type_search == 0 && cigar_code != NULL) {
      metaexon_insert(cal->strand, cal->chromosome_id - 1,
		      genome_start, genome_end, 40,
		      METAEXON_NORMAL, NULL,
		      metaexons);
    }
  }
  
  return cigar_code;

}

int extend_extrem_nt(char *reference, char *query, int extrem_type, int *distance) {
  int pos, action, limit, num_errors = 0;
  int limit_errors = strlen(query)/4, num_matches = 0;
  int consecutive_err = 0;
  int max_consecutive_err = 3;
  int number_nt = 0;

  //printf("Reference: %s\n", reference);
  //printf("Query    : %s\n", query);

  if (extrem_type == LEFT_EXTREM) {
    pos = 0;
    action = 1;
    limit = strlen(query);
  } else {
    pos = strlen(query) - 1;
    action = -1;
    limit = -1;
  }

  while (pos != limit) {
    if (reference[pos] != query[pos]) { 
      num_errors++;
      consecutive_err++;
      if (num_errors >= limit_errors || 
	  consecutive_err > max_consecutive_err) { number_nt = 0; break; }
    } else {
      num_matches++;
      consecutive_err = 0;
    }
    pos += action;
    number_nt++;
  }
  
  //printf("I come to pos %i with %i errors\n", number_nt, num_errors);
  *distance = num_errors;

  return number_nt;

}


float generate_cals_score(array_list_t *cals_list, int read_length) {
  int len_cal = 0;
  size_t num_cals = array_list_size(cals_list);
  cal_t *cal;
  linked_list_iterator_t itr;  
  seed_region_t *s, *s_prev;
  int i;

  //printf("==== CALS SCORE ====\n");
  for (i = 0; i < num_cals; i++) {    
    cal = array_list_get(i, cals_list);
    //printf("\tCAL (%i)[%i:%lu - %lu]:\n", cal->strand, cal->chromosome_id, cal->start, cal->end);
    // s_prev = NULL;

    linked_list_iterator_init(cal->sr_list, &itr);
    s = (seed_region_t *) linked_list_iterator_curr(&itr);
    while (s != NULL) {
      if (s->read_start > s->read_end) { assert(s->read_start); }
      len_cal += s->read_end - s->read_start;
      //if (s_prev) {
      //len_cal += s->read_start - s_prev->read_end;
      //}
      //s_prev = s;
      //printf("\t\tSEED %i:[%lu|%i-%i|%lu]Len CAL %i\n", i, s->genome_start, s->read_start, s->read_end, s->genome_end, len_cal);
      s = (seed_region_t *) linked_list_iterator_next(&itr);
    }

    if (array_list_size(cal->candidates_seeds_start)) {
      len_cal += 16;
      //printf("\tSeeds Candidates %i:\n", array_list_size(cal->candidates_seeds_start));
      for (int j = 0; j < array_list_size(cal->candidates_seeds_start); j++) {
	seed_region_t *seed_region = array_list_get(j, cal->candidates_seeds_start);
	//printf("\t\t\t Candidate Seed S:[%lu|%i-%i|%lu]\n", seed_region->genome_start, seed_region->read_start, seed_region->read_end, seed_region->genome_end);
      }
    }

    if (array_list_size(cal->candidates_seeds_end)) {
      len_cal += 16;
      //printf("\tSeeds Candidates %i:\n", array_list_size(cal->candidates_seeds_end));
      for (int j = 0; j < array_list_size(cal->candidates_seeds_end); j++) {
	seed_region_t *seed_region = array_list_get(j, cal->candidates_seeds_end);
	//printf("\t\t\t Candidate Seed E:[%lu|%i-%i|%lu]\n", seed_region->genome_start, seed_region->read_start, seed_region->read_end, seed_region->genome_end);
      }
    }
    //printf("\t<------- NEW CAL --------\n");
  }
  //printf("(SCORE %i,%i, %f)\n", len_cal, read_length, (float)(len_cal*100)/(float)read_length);
  //printf("<##### END FUNCTION #####>\n");

  return (float)(len_cal*100)/(float)read_length;

}

void order_cals(array_list_t *cals_list) {

  if (array_list_size(cals_list) <= 1) { return; }
 
  cal_t *cal_prev, *cal_next;
  for (int i = 0; i < array_list_size(cals_list) - 1; i++) {
    cal_prev = (cal_t *)array_list_get(i, cals_list);
    for (int j = i + 1; j < array_list_size(cals_list); j++) {
      cal_next = (cal_t *)array_list_get(j, cals_list);
      if (cal_next->strand < cal_prev->strand) {
	//printf("(%i-%i) %i[%i:%lu-%lu] <(strand)> %i[%i:%lu-%lu]\n", i, j, cal_prev->strand, cal_prev->chromosome_id, cal_prev->start, cal_prev->end, 
	//       cal_next->strand, cal_next->chromosome_id, cal_next->start, cal_next->end);
	array_list_swap(i, j, cals_list);
	cal_prev = cal_next;
      } else {
	if (cal_next->chromosome_id < cal_prev->chromosome_id) {
	  //printf("(%i-%i) %i[%i:%lu-%lu] <(chr)> %i[%i:%lu-%lu]\n", i, j, cal_prev->strand, cal_prev->chromosome_id, cal_prev->start, cal_prev->end, 
	  //	 cal_next->strand, cal_next->chromosome_id, cal_next->start, cal_next->end);
	  array_list_swap(i, j, cals_list);
	  cal_prev = cal_next;
	} else {
	  if (cal_next->start < cal_prev->start) {
	    //printf("(%i-%i) %i[%i:%lu-%lu] <(start)> %i[%i:%lu-%lu]\n", i, j, cal_prev->strand, cal_prev->chromosome_id, cal_prev->start, cal_prev->end, 
	    //	   cal_next->strand, cal_next->chromosome_id, cal_next->start, cal_next->end);
	    array_list_swap(i, j, cals_list);
	    cal_prev = cal_next;
	  }
	}
      }
    }
  }
}


int merge_and_filter_cals(array_list_t *cals_targets, 
			  array_list_t *cals_list, 
			  fastq_read_t *fq_read, 
			  bwt_optarg_t *bwt_optarg,
			  bwt_index_t *bwt_index, 
			  genome_t *genome, 
			  float *score_ranking) {

  int num_cals = array_list_size(cals_list);
  int cal_pos;
  cal_t *cal_prev, *cal_next, *cal;
  array_list_t *merge_cals, *fusion_cals;
  seed_region_t *s, *s_prev, *s_next;
  float *cals_score = score_ranking;  
  int number_of_best;
  size_t max_intron_size = 500000;
  float score;
  cigar_code_t *cigar_code = NULL;
  linked_list_iterator_t itr;  
  linked_list_t *linked_list;
  int seed_err_size = 20;
  char query[fq_read->length];
  char reference[fq_read->length];
  char *rev_comp = NULL;
  char *sequence;
  int number_nt;
  size_t genome_start, genome_end;
  int best_cals = num_cals;
  int distance = 0;

  const float MIN_CAL_SCORE = 60.0;
  register int i, j;

  //===== Step-1: Concatenate CALs =====//
  LOG_DEBUG("STEP-1: CONCATENATE CALS");
  merge_cals = array_list_new(10,
			      1.25f,
			      COLLECTION_MODE_ASYNCHRONIZED);
  cal_pos = 0;
  cal_prev = (cal_t *)array_list_get(cal_pos++, cals_list);
  
  array_list_insert(cal_prev, merge_cals);    
  while (cal_pos < num_cals) {
    cal_next = (cal_t *)array_list_get(cal_pos, cals_list);      
    s_prev = linked_list_get_last(cal_prev->sr_list);
    s = linked_list_get_first(cal_next->sr_list);

    if (cal_prev->chromosome_id == cal_next->chromosome_id && 
	cal_prev->strand == cal_next->strand && 
	s_prev->read_end <= s->read_start &&
	(cal_next->start <= (cal_prev->end + max_intron_size))) {
      //printf("Merge!! cal_prev->end = %lu, cal_next->start = %lu\n", cal_prev->end, cal_next->start);
      array_list_insert(cal_next, merge_cals);
    } else { 
      array_list_insert(merge_cals, cals_targets);
      merge_cals = array_list_new(10,
				  1.25f,
				  COLLECTION_MODE_ASYNCHRONIZED);
      array_list_insert(cal_next, merge_cals);
    }                                                                                         
    cal_prev = cal_next;
    cal_pos++;
  }
  array_list_insert(merge_cals, cals_targets);         
  array_list_clear(cals_list, NULL);
  //===== Step-1: END =====//


  //===== Step-2: Generate CALs Score =====//
  LOG_DEBUG("STEP-2: GENERATE CALS SCORE");
  for (i = 0; i < array_list_size(cals_targets); i++) {
    fusion_cals = array_list_get(i, cals_targets);
    cals_score[i] = generate_cals_score(fusion_cals, fq_read->length);
    //printf("SCORE %f\n", cals_score[i]);
  }
  //===== Step-2: END =====//

  /*num_cals = array_list_size(cals_targets);
  if (best_cals > num_cals) {
    number_of_best = num_cals;
  } else {
    number_of_best = best_cals;
    }*/

  number_of_best = array_list_size(cals_targets);
  
  //===== Step-3: Ranking CALs by score (the n best)=====//
  LOG_DEBUG("STEP-3: ORDER CALS");
  for (i = 0; i < number_of_best; i++) {
    for (j = i + 1; j < array_list_size(cals_targets); j++) {
      if (cals_score[j] > cals_score[i]) {
	array_list_swap(i, j, cals_targets);
	score = cals_score[j];
	cals_score[j] = cals_score[i];
	cals_score[i] = score;
      }
    }
  }
  //===== Step-3: END =====//

  //===== Step-4: Search Near Seeds of CALs =====//
  for (i = 0; i < number_of_best; i++) {
    fusion_cals = array_list_get(i, cals_targets);
    // Has this CAL one or more seeds near?
    cal_t *cal_prev = array_list_get(0, fusion_cals);
    cal_t *cal_next = array_list_get(array_list_size(fusion_cals) - 1, fusion_cals);
    
    //TODO: FIRST SEED? NO... SEARCH THE BEST SEEDS
    if (array_list_size(cal_prev->candidates_seeds_start) > 0) {
      seed_region_t *seed_region = array_list_remove_at(0, cal_prev->candidates_seeds_start);
      linked_list_t *linked_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
      linked_list_insert_first(seed_region, linked_list);

      cal = cal_new(cal_prev->chromosome_id, cal_prev->strand, 
		    seed_region->genome_start, seed_region->genome_end, 
		    1, linked_list, 
		    linked_list_new(COLLECTION_MODE_ASYNCHRONIZED));
      array_list_insert_at(0, cal, fusion_cals);
    }

    if (array_list_size(cal_next->candidates_seeds_end) > 0) {
      seed_region_t *seed_region = array_list_remove_at(0, cal_next->candidates_seeds_end);
      linked_list_t *linked_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
      linked_list_insert_first(seed_region, linked_list);

      cal = cal_new(cal_next->chromosome_id, cal_next->strand, 
		    seed_region->genome_start, seed_region->genome_end, 
		    1, linked_list,
		    linked_list_new(COLLECTION_MODE_ASYNCHRONIZED));
      array_list_insert(cal, fusion_cals);
    }

    /*for (j = 0; j < array_list_size(fusion_cals); j++) {
      cal = array_list_get(j, fusion_cals);            
      array_list_insert(cal, cals_list);
    }
    */
  }
  
  //Extend extrem gaps to final    
  /*for (i = 0; i < number_of_best; i++) {
    fusion_cals = array_list_get(i, cals_targets);
    cal_prev = array_list_get(0, fusion_cals);
    cal_next = array_list_get(array_list_size(fusion_cals) - 1, fusion_cals);

    s_prev = linked_list_get_first(cal_prev->sr_list);
    s_next = linked_list_get_last(cal_next->sr_list);

    int nt_tot = 0;
    //printf("Extend [%i:%lu-%lu]??\n", cal_prev->chromosome_id, cal_prev->start, cal_next->end);
    if (s_prev->read_start > 20) {
      //printf("\tExtend left extrem\n");
      if (cal_prev->strand == 1) {
	if (!rev_comp ) {
	  rev_comp = (char *) calloc(fq_read->length + 1, sizeof(char));
	  strcpy(rev_comp, fq_read->sequence);
	  seq_reverse_complementary(rev_comp, fq_read->length);
	}
	sequence = rev_comp;
      } else {
	sequence = fq_read->sequence;
      }

      genome_start = cal_prev->start - s_prev->read_start + 1;
      genome_end = cal_prev->start;
      genome_read_sequence_by_chr_index(reference, 0, 
					cal_prev->chromosome_id - 1, &genome_start, &genome_end, genome);      
      memcpy(query, sequence, s_prev->read_start);
      query[s_prev->read_start] = '\0';
      number_nt = extend_extrem_nt(reference, query, LEFT_EXTREM, &distance);
      if (number_nt == 0) {
	q[*num_sw] = strdup(query);
	r[*num_sw] = strdup(reference);
	fusion_coords[(*num_sw)++] = fusion_coords_new(0, 0, 0, 0, 0, 0, 0, 
						       0, FIRST_SW, fq_read->id, cal_prev, read_id, i);
	//printf("Insert Reference CAL %i, Read %i Start\n", i, read_id);
      } else {
	s_prev->read_start -= number_nt;
	s_prev->genome_start -= number_nt;
	cal_prev->start -= number_nt;      
	nt_tot += number_nt;
	cals_score[i] += (((nt_tot - distance)*100)/fq_read->length);
	//printf("Actualization score\n");
      }
    }
    
    if (s_next->read_end < fq_read->length - 20) {
      //printf("\tExtend right extrem\n");
      if (cal_next->strand == 1) {
	if (!rev_comp ) {
	  rev_comp = (char *) calloc(fq_read->length + 1, sizeof(char));
	  strcpy(rev_comp, fq_read->sequence);
	  seq_reverse_complementary(rev_comp, fq_read->length);
	}
	sequence = rev_comp;
      } else {
	sequence = fq_read->sequence;
      }

      genome_start = cal_next->end;
      genome_end = cal_next->end + (fq_read->length - s_next->read_end) - 1;
      genome_read_sequence_by_chr_index(reference, 0, 
					cal_next->chromosome_id - 1, &genome_start, &genome_end, genome);

      //printf("#############From %i:%lu-%lu: %s\n", cal_next->chromosome_id - 1, genome_start, genome_end, reference);

      memcpy(query, sequence + s_next->read_end, fq_read->length - s_next->read_end);
      query[fq_read->length - s_next->read_end] = '\0';
      number_nt = extend_extrem_nt(reference, query, RIGHT_EXTREM, &distance);
      if (number_nt == 0) {
	q[*num_sw] = strdup(query);
	r[*num_sw] = strdup(reference);
	fusion_coords[(*num_sw)++] = fusion_coords_new(0, 0, 0, 0, 0, 0, 0, 0, 
						       LAST_SW, fq_read->id, cal_next, read_id, i);
	//printf("Insert Reference CAL %i, Read %i End\n", i, read_id);
      } else {
	s_next->read_end += number_nt;
	s_next->genome_end += number_nt;
	cal_next->end += number_nt;	
	nt_tot += number_nt;      
	//printf("Actualization score %i nt, %i errors\n", number_nt, distance);
	cals_score[i] += (((nt_tot - distance)*100)/fq_read->length);
      }
    }    
  }*/
  
  //The best CAL has a start/end big gap?
  //printf("BEST SCORE %f\n", cals_score[0]);
  /*
    if (cals_score[0] <= MIN_CAL_SCORE) {
    fusion_cals = array_list_get(0, cals_targets);
    cal = array_list_get(0, fusion_cals);
    s_prev = linked_list_get_first(cal->sr_list);
    s_next = linked_list_get_last(cal->sr_list);
  
    //printf("$$$Make Seeds with one Error!! FIRST SEED:[%i-%i] LAST_SEED:[%i-%i] %s\n", s_prev->read_start, s_prev->read_end,
    //	   s_next->read_start, s_next->read_end, fq_read->id);
    
    if (s_prev->read_start >= seed_err_size) {
      printf("\t @@@@SEEDs in First positions [%i-%i]\n", 0, s_prev->read_start - 1);
      //Seeds in Start position
      array_list_t *mapping_list_prev = array_list_new(1000,
						       1.25f,
						       COLLECTION_MODE_ASYNCHRONIZED);
      
      bwt_map_inexact_seeds_by_region(0, s_prev->read_start - 1,
				      cal->chromosome_id, cal->start - 500000,
				      cal->start,
				      fq_read->sequence, seed_err_size,
				      seed_err_size,
				      bwt_optarg,
				      bwt_index,
				      mapping_list_prev);
      
      
      for (int r = 0; r < array_list_size(mapping_list_prev); r++) {
	region_t *region = array_list_get(r, mapping_list_prev);
	printf("\t Region [%i:%lu-%lu]\n", region->chromosome_id, region->start, region->end);
      }

      array_list_free(mapping_list_prev, region_bwt_free);

    }

    if (s_next->read_end <= fq_read->length - seed_err_size) {
      printf("\t @@@@SEEDs in Last positions [%i-%i]\n", s_next->read_end, fq_read->length - 1);
      //Seeds in End position
      array_list_t *mapping_list_next = array_list_new(1000,
						       1.25f,
						       COLLECTION_MODE_ASYNCHRONIZED);
     
      bwt_map_inexact_seeds_by_region(s_next->read_end, fq_read->length - 1,
				      cal->chromosome_id, cal->end,
				      cal->end + 500000,
				      fq_read->sequence, seed_err_size,
				      seed_err_size/2,
				      bwt_optarg,
				      bwt_index,
				      mapping_list_next);
      
      for (int r = 0; r < array_list_size(mapping_list_next); r++) {
	region_t *region = array_list_get(r, mapping_list_next);
	printf("\t Region [%i:%lu-%lu]\n", region->chromosome_id, region->start, region->end);
      }

      array_list_free(mapping_list_next, region_bwt_free);
    }
    }*/
  //End

  //Delete other CALs
  /*for (i = array_list_size(cals_targets) - 1; i >= number_of_best; i--) {
    fusion_cals = array_list_remove_at(i, cals_targets);
    for (j = 0; j < array_list_size(fusion_cals); j++) {
      cal = array_list_get(j, fusion_cals);
      cigar_code = (cigar_code_t *)cal->info;
      if (cigar_code != NULL) {
	array_list_clear(cigar_code->ops, cigar_op_free);
	cigar_code_free(cigar_code);
      }
      cal_free(cal);
    }
    array_list_free(fusion_cals, NULL);
  }
  */
  //free(rev_comp);

  return number_of_best;

}




void generate_reference_splice_juntion(array_list_t *cals_targets, char *query_revcomp, 
				       fastq_read_t *fq_read, char **r, char **q, 
				       fusion_coords_t **fusion_coords, int *num_sw, 
				       genome_t *genome, size_t id_read) {
  int n_fusion_cals = array_list_size(cals_targets);

  cal_t *cal_prev, *cal_next, *cal;
  array_list_t *fusion_cals;
  char *query_map;
  seed_region_t *s, *s_prev;
  char reference[2048];
  char reference_prev[2048];
  char reference_next[2048];
  char query[2048];

  size_t genome_start, genome_end;
  int read_start, read_end;
  size_t genome_start2, genome_end2;
  int seeds_nt;
  int flank = 30;

  register int i, j;

  //Delete for debuging, detect splice junctions
  for (i = 0; i < n_fusion_cals; i++) {
    fusion_cals = array_list_get(i, cals_targets);
    j = 0;
    cal_prev = array_list_get(j, fusion_cals);
    s_prev = linked_list_get_last(cal_prev->sr_list);
    if (cal_prev->strand == 1) {
      query_map = query_revcomp;
    } else {
      query_map = fq_read->sequence;
    }

    //printf("SP_CAL_PREV(%i) %i:%lu-%lu %i:%i=%i\n", cal->strand, cal->chromosome_id, cal_prev->end, 
    //     cal->start, read_start, read_end, seeds_nt);
    
    for (j = 1; j < array_list_size(fusion_cals); j++) {
      cal = array_list_get(j, fusion_cals);
      s = linked_list_get_first(cal->sr_list);

      read_start = s_prev->read_end;
      read_end = s->read_start;
      seeds_nt = read_end - read_start;
      
      //printf("SP_CAL(%i) %i:%lu-%lu %i:%i=%i\n", cal->strand, cal->chromosome_id, cal_prev->end, 
      //   cal->start, read_start, read_end, seeds_nt);
      
      if (seeds_nt >= 10) {
	//Extend to Right --> <-- Extend to Left
	genome_start = s_prev->genome_end;
	genome_end = s_prev->genome_end + seeds_nt - 1;
	genome_read_sequence_by_chr_index(reference_prev, 0, 
					  cal->chromosome_id - 1, &genome_start, &genome_end, genome);

	genome_start2 = s->genome_start - seeds_nt;
	genome_end2 = s->genome_start - 1;
	genome_read_sequence_by_chr_index(reference_next, 0, 
					  cal->chromosome_id - 1, &genome_start2, &genome_end2, genome);

	memcpy(query, query_map + read_start, read_end - read_start);
	query[read_end - read_start]  = '\0';

	//CAll new function
	int dsp_e1, dsp_e2;
	int lim_err = 3;

	extend_by_mismatches(reference_prev, reference_next, query, 
			     0, strlen(reference_next) - 1, 
			     0, read_end - read_start - 1, lim_err,
			     &dsp_e1, &dsp_e2);
	//printf("dsp_e1=%i, dsp_e2=%i\n", dsp_e1, dsp_e2);

	cal_prev->end = cal_prev->end + dsp_e1;	
	cal->start = cal->start - dsp_e2;

	s_prev->read_end = s_prev->read_end + dsp_e1;
	s->read_start = s->read_start - dsp_e2;

	if (s_prev->read_end > s->read_start) { 
	  seeds_nt = s_prev->read_end - s->read_start;
	  cal_prev->end -= seeds_nt;
	  cal->start += seeds_nt;
	  s_prev->read_end-= seeds_nt;
	  s->read_start += seeds_nt;
	}
      }
      
      read_start = s_prev->read_end - flank;
      if (read_start < 0) { read_start = 0; }
      read_end = s->read_start + flank;
      if (read_end >= fq_read->length) { read_end = fq_read->length - 1; }

      //Extract and fusion Reference SW
      genome_start = cal_prev->end - flank;
      genome_end = cal_prev->end + flank - 1;
      genome_read_sequence_by_chr_index(reference_prev, 0, 
					cal_prev->chromosome_id - 1, &genome_start, &genome_end, genome);

      //printf("Ref Prev %i\n", strlen(reference_prev));
      genome_start2 = cal->start - flank;
      genome_end2 = cal->start + flank - 1;
      genome_read_sequence_by_chr_index(reference_next, 0, 
					cal->chromosome_id - 1, &genome_start2, &genome_end2, genome);
      
      //printf("Ref Next %i [%i:%i-%i]=%s\n", strlen(reference_next), cal->chromosome_id, genome_start2, genome_end2, reference_next);
      //if (genome_start2 < genome_start) { printf("%s\n", fq_read->id); assert(genome_start); }
      //printf("\tSP_CAL(%i) %i:%lu-%lu %i:%i\n", cal->strand, cal->chromosome_id, cal_prev->end, 
      //   cal->start, read_start, read_end);
      
      strcat(reference_prev, reference_next);

      if (read_start > read_end) { printf("ALERT ERROR:: %s\n", fq_read->id);exit(-1); }
      //printf("Read %i-%i: %s\n", read_start, read_end, query_map);

      memcpy(query, query_map + read_start, read_end - read_start);
      query[read_end -read_start] = '\0';

      r[*num_sw] = strdup(reference_prev);
      q[*num_sw] = strdup(query);

      //printf("Que: %s\n", q[*num_sw]);
      //printf("Ref: %s\n", r[*num_sw]);
      //printf("Ref len %i, Query len %i\n", strlen(r[num_sw]), strlen(q[num_sw]));
      fusion_coords[(*num_sw)++] = fusion_coords_new(genome_start, genome_end,
						     genome_start2, genome_end2, 
						     read_start, read_end,
						     cal->chromosome_id, cal->strand,
						     MIDDLE_SW, fq_read->id, cal_prev, id_read, 0);	
      
      cal_prev = cal;
      s_prev = linked_list_get_last(cal_prev->sr_list);
      
      }
  }
}


/*
array_list_t *fusion_regions (array_list_t *regions_list, int max_distance) {
  int r = 0;
  region_t *region, *region_next;
  array_list_t *merge_regions_list = array_list_new(array_list_size(regions_list),
						   1.25f,
						   COLLECTION_MODE_ASYNCHRONIZED);
  int last_id;

  //printf("Regions Found %i:\n", array_list_size(regions_list));
  while (array_list_size(regions_list) > 0) {
    region_t *region = array_list_get(r, regions_list);
    //printf("\t::: %i-Region[%lu|%i-%i|%lu] (%i/%i)\n", 
    //	   region->id, region->start, region->seq_start, region->seq_end, region->end, 
    //	   r, array_list_size(regions_list));
    if (r + 1 < array_list_size(regions_list)) {
      region_t *region_next = array_list_get(r + 1, regions_list);
      if (region->id == region_next->id) {
	//TODO: Equal seeds id select not do this
	r++;
	continue;
      } else if ((region->end + max_distance) >= region_next->start) {
	//Fusion seeds
	region_next->start = region->start;
	region_next->seq_start = region->seq_start;
      } else {
	//Not fusion
	//printf("\t\tInsert region\n");
	array_list_insert(region_bwt_new(region->chromosome_id,
					 region->strand,
					 region->start,
					 region->end,
					 region->seq_start,
					 region->seq_end,
					 region->seq_len,
					 region->id), 
			  merge_regions_list);
      }
    } else {
      //printf("\t\tInsert region\n");
      array_list_insert(region_bwt_new(region->chromosome_id,
				       region->strand,
				       region->start,
				       region->end,
				       region->seq_start,
				       region->seq_end,
				       region->seq_len,
				       region->id),
			merge_regions_list);	    
      break;
    }
    r++;
  }

  array_list_clear(regions_list, region_bwt_free);

  for (int i = 0; i < array_list_size(merge_regions_list); i++) {
    region_t *region = array_list_get(i, merge_regions_list);
    array_list_insert(region, regions_list);
  }
  
  array_list_free(merge_regions_list, NULL);

  return regions_list;
  
}
*/

//============================ STRUCTURES AND TYPES SECTION =============================//
int search_simple_splice_junction(seed_region_t *s_prev, seed_region_t *s_next,
				  int chromosome_id, int strand, 
				  char *sequence, genome_t *genome, 
				  avls_list_t *avls_list, 
				  avl_node_t **node_avl_start,
				  avl_node_t **node_avl_end,
				  int  *distance);


typedef struct info_sp {
  size_t l_genome_start;
  size_t l_genome_end;
  size_t r_genome_start;
  size_t r_genome_end;
} info_sp_t;


typedef struct sw_item {
  int type_sw;
  int read_id;
  int fusion_id;
  int cal_id;
  seed_region_t *seed_prev;
  seed_region_t *seed_next;
  cal_t *cal_prev;
  cal_t *cal_next;
  meta_alignment_t *meta_alignment;
  void *info;
} sw_item_t;

sw_item_t *sw_item_new(int type_sw, int read_id, 
		       int fusion_id, int cal_id,
		       cal_t *cal_prev, cal_t *cal_next, 
		       meta_alignment_t *meta_alignment, 
		       seed_region_t *seed_prev,
		       seed_region_t *seed_next,
		       void *info);

typedef struct sw_depth {
  char *q[MAX_DEPTH];
  char *r[MAX_DEPTH];
  sw_item_t *items[MAX_DEPTH];
  int depth;
} sw_depth_t;

void sw_depth_insert(char *query, char *reference, sw_item_t *sw_item, 
		     sw_optarg_t *sw_optarg, sw_multi_output_t *output, 
		     avls_list_t *avls_list, metaexons_t *metaexons, 
		     sw_depth_t *sw_depth);

//=======================================================================================//


//=============================== META ALIGNMENT SECTION ================================//
meta_alignment_t *meta_alignment_new() {
  meta_alignment_t *meta_alignment = (meta_alignment_t *)malloc(sizeof(meta_alignment_t));

  meta_alignment->cals_list = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  meta_alignment->status = META_OPEN;
  meta_alignment->sp_sw = 0;
  meta_alignment->num_cigars = 0;
  meta_alignment->cigar_code = cigar_code_new();
  meta_alignment->type = META_ALIGNMENT_NONE;
  meta_alignment->cigar_left = NULL;
  meta_alignment->cigar_right = NULL;
  meta_alignment->score = 0;
  meta_alignment->flag = 0;
  return meta_alignment;
}

void meta_alignment_free(meta_alignment_t *meta_alignment) {
  //if (meta_alignment->cals_list != NULL) { array_list_free(meta_alignment->cals_list, NULL); }
  cigar_code_free(meta_alignment->cigar_code);
  free(meta_alignment);
}

int meta_alignment_num_cals(meta_alignment_t *meta_alignment) {
  return array_list_size(meta_alignment->cals_list);
}

cal_t *meta_alignment_get_first_cal(meta_alignment_t *meta_alignment) {
  return array_list_get(0, meta_alignment->cals_list);
}

cal_t *meta_alignment_get_last_cal(meta_alignment_t *meta_alignment) {
  return array_list_get(meta_alignment_num_cals(meta_alignment) - 1, meta_alignment->cals_list);
}


int meta_alignment_num_cigars(meta_alignment_t *meta_alignment) {
  if (meta_alignment != NULL) {
    return meta_alignment->num_cigars;
  } else {
    return 0;
  }
}

int meta_alignment_set_status(int status, meta_alignment_t *meta_alignment) {
  meta_alignment->status = status;
  return status;
}

int meta_alignment_get_status(meta_alignment_t *meta_alignment) {
  return meta_alignment->status;
}

int meta_alignment_insert_cal(cal_t *cal, meta_alignment_t *meta_alignment) {
  return array_list_insert(cal, meta_alignment->cals_list);
}


void meta_alignments_order_by_score(array_list_t *meta_alignments) {
  meta_alignment_t *meta_prev, *meta_next;
  for (int i = 0; i < array_list_size(meta_alignments) - 1; i++) {
    meta_prev = array_list_get(i, meta_alignments);
    for (int j = i + 1; j < array_list_size(meta_alignments); j++) {
      meta_next = array_list_get(j, meta_alignments);
      if (meta_next->score > meta_prev->score) {
	array_list_swap(i, j, meta_alignments);
      }
    }
  }
  
}

void meta_alignment_insert_cigar(cigar_code_t *cigar, int type, int pos, meta_alignment_t *meta_alignment) {  
  //if (cigar != NULL && (type == CIGAR_ANCHOR_LEFT || type == CIGAR_ANCHOR_RIGHT)) {
    //printf("-----------------------------> INSERT %s. (%i)\n", new_cigar_code_string(cigar), array_list_size(cigar->ops));
  //}

  if (meta_alignment == NULL) { return; }
  else {
    if (type == CIGAR_ANCHOR_LEFT) {
      meta_alignment->cigar_right = cigar;
    } else if (type == CIGAR_ANCHOR_RIGHT) {
      meta_alignment->cigar_left = cigar;
    } else {
      //printf("Insert middle cigar %i => %s\n", pos, new_cigar_code_string(cigar));
      if (pos > 20 || pos < 0) { exit(-1); }
      meta_alignment->middle_cigars[pos] = cigar;
      meta_alignment->type_cigars[pos] = type;
      meta_alignment->num_cigars++;
    }
  }
}

void meta_alignment_calculate_score(meta_alignment_t *meta_alignment) {
  cigar_code_t *cigar_code = meta_alignment->cigar_code;
  int num_di = 0;
  meta_alignment->score = 0;
  if (cigar_code == NULL) { printf("NULL CIGAR\n"); return; }
  //printf("NUM OPS: %i\n", array_list_size(cigar_code->ops));
  for (int i = 0; i < array_list_size(cigar_code->ops); i++) {
    cigar_op_t *op = array_list_get(i, cigar_code->ops);
    if (op->name == 'M' || op->name == 'I') { 
      meta_alignment->score += op->number; 
    }
    //if (op->name == 'I' || op->name == 'D') {
    //num_di += op->number;
    //}
    //printf("SCORE (%i) OP+ : %i%c\n", meta_alignment->score, 
    //	   op->number, op->name);
  }
}


int meta_alignment_get_cals_score(meta_alignment_t *meta_alignment) {
  int num_di = 0, score = 0;
  
  for (int i = 0; i < array_list_size(meta_alignment->cals_list); i++) {
    cal_t *cal = array_list_get(i, meta_alignment->cals_list);
    linked_list_item_t *item = cal->sr_list->first;
    while (item != NULL) {
      seed_region_t *seed_aux = item->item;
      score += seed_aux->read_end - seed_aux->read_start;
      item= item->next;
    }
  }

  return score;

}

void meta_alignment_close(meta_alignment_t *meta_alignment) {
  cal_t *first_cal, *cal;
  cigar_code_t *cigar_code, *cigar_code_aux;
  linked_list_item_t *list_item;
  seed_region_t *s, *s_prev;
  int type;
  int cr_pos = 0;
  cigar_op_t *op;
  int bad_cigar = 0;

  //printf("\n==================== CLOSE META ALIGNMENT ==========================\n");
  assert(meta_alignment_num_cals(meta_alignment) > 0);
  cigar_code = meta_alignment->cigar_code;

  if (cigar_code_get_num_ops(cigar_code) > 0) {
    array_list_clear(cigar_code->ops, cigar_op_free);
    cigar_code->distance = 0;
  }

  //cal_t *cal_prev = NULL;
  if (meta_alignment->type == META_ALIGNMENT_MIDDLE) {
    //printf("META_ALIGNMENT_MIDDLE %i\n", meta_alignment->type);
    if (meta_alignment->cigar_left != NULL) {
      cigar_code_aux = meta_alignment->cigar_left;
      for (int i = 0; i < array_list_size(cigar_code_aux->ops); i++) {
	cigar_op_t *op = array_list_get(i, cigar_code_aux->ops);
	//printf("\t OP-SP-LEFT: %i%c\n", op->number, op->name);
	cigar_code_append_new_op(op->number, op->name, cigar_code);
      } 
      cigar_code->distance += cigar_code_aux->distance;
    }

    for (int i = 0; i < meta_alignment_num_cals(meta_alignment); i++) {
      cal_t *cal = array_list_get(i, meta_alignment->cals_list);    
      //printf(" ***** CLOSE GAP *****\n");
      //cal_print(cal);
      //printf(" *********************\n");
      cigar_code_aux = cal->info;

      if (cr_pos > 0 &&
	  meta_alignment->type_cigars[cr_pos - 1] == CIGAR_SW_MIDDLE) {
	cigar_code_t *c_c = cigar_code_new();
	for (int t = 0; t < array_list_size(cigar_code_aux->ops); t++) {
	  op = array_list_get(t, cigar_code_aux->ops);
	  cigar_code_append_new_op(op->number, op->name, c_c);
	}
	cigar_code_delete_nt(cal->r_flank, 0, c_c);	  
	for (int j = 0; j < cigar_code_get_num_ops(c_c); j++) {
	  op = array_list_get(j, c_c->ops);
	  //printf("\t OP-->1: %i%c\n", op->number, op->name);
	  cigar_code_append_op(op, cigar_code);
	}
      } else {      
	for (int j = 0; j < cigar_code_get_num_ops(cigar_code_aux); j++) {
	  op = array_list_get(j, cigar_code_aux->ops);
	  //printf("\t OP-->2: %i%c\n", op->number, op->name);
	  cigar_code_append_new_op(op->number, op->name, cigar_code);
	  //cigar_code_append_op(op, cigar_code);
	}
      }

      cigar_code->distance += cigar_code_aux->distance;

      if (cr_pos < meta_alignment_num_cigars(meta_alignment)) {
	if (meta_alignment->type_cigars[cr_pos] == CIGAR_SW_MIDDLE) {
	  cigar_code_delete_nt(cal->l_flank, 1, cigar_code);	  
	}

	//printf("CLOSE-META: MIDDLE SPLICE\n");
	cigar_code_aux = meta_alignment->middle_cigars[cr_pos++];
	if (cigar_code_aux == NULL) { bad_cigar = 1; break; }

	for (int i = 0; i < array_list_size(cigar_code_aux->ops); i++) {
	  cigar_op_t *op = array_list_get(i, cigar_code_aux->ops);
	  //cigar_code_append_op(op, cigar_code);
	  cigar_code_append_new_op(op->number, op->name, cigar_code);
	  //printf("\t OP-SP-M : %i%c\n", op->number, op->name);
	}
      }
    }

    if (meta_alignment->cigar_right != NULL) {
      cigar_code_aux = meta_alignment->cigar_right;
      for (int i = 0; i < array_list_size(cigar_code_aux->ops); i++) {
	cigar_op_t *op = array_list_get(i, cigar_code_aux->ops);
	//printf("\t OP-SP-RIGHT: %i%c\n", op->number, op->name);
	cigar_code_append_new_op(op->number, op->name, cigar_code);
      }
      cigar_code->distance += cigar_code_aux->distance;
    }
  } else {    
    //printf("META_ALIGNMENT_OTHER %i\n", meta_alignment->type);
    if (meta_alignment->cigar_left != NULL) {
      cigar_code_aux = meta_alignment->cigar_left;
      for (int i = 0; i < array_list_size(cigar_code_aux->ops); i++) {
	cigar_op_t *op = array_list_get(i, cigar_code_aux->ops);
	//printf("\t OP-SP-RIGHT: %i%c\n", op->number, op->name);
	//cigar_code_append_op(op, cigar_code);
	cigar_code_append_new_op(op->number, op->name, cigar_code);
      } 
      cigar_code->distance += cigar_code_aux->distance;
    }

    cal_t *cal = array_list_get(0, meta_alignment->cals_list);    
    //printf(" ***** CLOSE GAP *****\n");
    //cal_print(cal);
    //printf(" *********************\n");
    cigar_code_aux = cal->info;  
    //printf("SEEDS OPS: %i\n", cigar_code_get_num_ops(cigar_code_aux));
    for (int j = 0; j < cigar_code_get_num_ops(cigar_code_aux); j++) {
      op = array_list_get(j, cigar_code_aux->ops);
      //printf("\t OP: %i%c\n", op->number, op->name);
      //cigar_code_append_op(op, cigar_code);
      cigar_code_append_new_op(op->number, op->name, cigar_code);
    } 

    cigar_code->distance += cigar_code_aux->distance;

    if (meta_alignment->cigar_right != NULL) {
      cigar_code_aux = meta_alignment->cigar_right;
      for (int i = 0; i < array_list_size(cigar_code_aux->ops); i++) {
	cigar_op_t *op = array_list_get(i, cigar_code_aux->ops);
	//printf("\t OP-SP-LEFT: %i%c\n", op->number, op->name);
	//cigar_code_append_op(op, cigar_code);
	cigar_code_append_new_op(op->number, op->name, cigar_code);
      } 
      cigar_code->distance += cigar_code_aux->distance;
    } 
  }

  if (meta_alignment->cigar_code == NULL) { bad_cigar = 1; }
  
  if (!bad_cigar) {
    meta_alignment_set_status(META_CLOSE, meta_alignment);
    meta_alignment_calculate_score(meta_alignment);
  }
  //printf("----- META CLOSE INSERT %s.\n", new_cigar_code_string(meta_alignment->cigar_code));
  //printf("------------------------------------------------------------------\n");
}

void merge_seeds_cal(cal_t *cal) {
  seed_region_t *seed_prev = NULL, *seed_next;
  linked_list_item_t *list_item = cal->sr_list->first, *list_item_prev;
  //if (list_item == NULL) { 
  //printf("ALARMM NULL, SIZE LIST\n");
  //}
  cigar_code_t *cigar_code = cigar_code_new();
  cigar_op_t *op;
  
  while (list_item != NULL) {
    seed_next = (seed_region_t *)list_item->item;
    //printf(" Seed ---> %i-%i\n", seed_next->read_start, seed_next->read_end);
    if (seed_prev != NULL) {
      if (seed_prev->fusion_right == 1 &&
	  seed_next->fusion_left == 1) {
	if (seed_prev->info == NULL) {
	  op = cigar_op_new(seed_prev->read_end - seed_prev->read_start + 1, 'M');
	  cigar_code_append_op(op, cigar_code);
	} else {
	  cigar_code_t *cigar_code_aux = seed_prev->info;
	  for (int i = 0; i < array_list_size(cigar_code_aux->ops); i++) {
	    op = array_list_get(i, cigar_code_aux->ops);
	    cigar_code_append_op(op, cigar_code);
	  }
	  cigar_code->distance += cigar_code_aux->distance;
	  cigar_code_free(cigar_code_aux);
	  seed_prev->info = NULL;
	}
      } else {
	//Detect splice
	op = cigar_op_new(seed_prev->read_end - seed_prev->read_start + 1, 'M');
	cigar_code_append_op(op, cigar_code);

	//printf("MINI SPLICE %lu - %lu = %lu\n", seed_prev->genome_end, 
	//     seed_next->genome_start, 
	//     seed_next->genome_start - seed_prev->genome_end + 1);
	op = cigar_op_new(seed_next->genome_start - seed_prev->genome_end - 1, 'N');
	cigar_code_append_op(op, cigar_code);
      } 
    }
    seed_prev = seed_next;
    list_item = list_item->next;
  }

  if (seed_prev == NULL) { exit(-1); }
  op = cigar_op_new(seed_prev->read_end - seed_prev->read_start + 1, 'M');
  cigar_code_append_op(op, cigar_code);

  //printf("FILL GAPS CLOSE CAL [%i:%lu-%lu]: %s\n", cal->chromosome_id, cal->start, cal->end, new_cigar_code_string(cigar_code));

  cal->info = cigar_code;
  
}

void meta_alignment_fill_gaps(int meta_type,
			      meta_alignment_t *meta_alignment, 
			      char *query_map, genome_t *genome,
			      sw_optarg_t *sw_optarg, sw_multi_output_t *output,			     
			      metaexons_t *metaexons, 
			      sw_depth_t *sw_depth,
			      avls_list_t *avls_list) {
  cigar_code_t *cigar_code;
  char reference[2048];
  int distance;
  int first, last;
  int min_distance;
  int closed;  
  avl_node_t *node_avl_start;
  avl_node_t *node_avl_end;
  int add_seed, sw_add, sp_found;
  sw_item_t *sw_item;
  const int FLANK = 2;	    

  meta_alignment->type = meta_type;
  cigar_code = meta_alignment->cigar_code;

  //printf("FILL GAPS TOTAL CALS %i\n", array_list_size(meta_alignment->cals_list));
  for (int i = 0; i < array_list_size(meta_alignment->cals_list); i++) {
    //printf("Process CAL %i\n", i);
    cal_t *cal = array_list_get(i, meta_alignment->cals_list);
    //printf("===== FILL GAPS CAL: =====\n");
    //cal_print(cal);
    //printf("=========================\n");
    //fill gaps
    cal->num_targets = 0;
    seed_region_t *seed_prev = NULL, *seed_next;
    linked_list_item_t *list_item = cal->sr_list->first, *list_item_prev;    
    while (list_item != NULL) {
      seed_next = (seed_region_t *)list_item->item;
      int num_match = seed_next->read_end - seed_next->read_start + 1;
      if (seed_prev != NULL) {
	//Close nt
	assert(seed_next->read_start > seed_prev->read_end);
	assert(seed_next->genome_start > seed_prev->genome_end);
	size_t gap_read = seed_next->read_start - seed_prev->read_end  - 1;
	//assert(gap_read >= 0);
	size_t gap_genome = seed_next->genome_start - seed_prev->genome_end  - 1; 
	//assert(gap_genome > 0);

	//printf("FILL GAPS ::: (gap_read = %lu), (gap_genome = %lu)\n", gap_read, gap_genome);
	distance = 0;
	closed = 0;
	add_seed = 0;
	sw_add = 0;
	sp_found = 0;
	if (gap_read == gap_genome &&
	    gap_read == 0) {
	  closed = 1;
	} else if (gap_read == gap_genome) {
	  size_t genome_end   = seed_next->genome_start - 1;
          size_t genome_start = seed_prev->genome_end + 1;
	  int read_end        = seed_next->read_start - 1;
	  int read_start      = seed_prev->read_end + 1;
	  char *query         = &query_map[read_start];
	  first = -1;
	  last = -1;
          genome_read_sequence_by_chr_index(reference, 0, cal->chromosome_id - 1,
                                            &genome_start, &genome_end, genome);
	  //printf("[%lu|%i]-GAP-[%i|%lu]: %s\n", genome_start, read_start, read_end, genome_end, reference);
	  for (int k = 0; k < gap_read; k++) {
	    //printf("[q:%c vs r:%c]\n", query[k], reference[k]);
	    if (query[k] != reference[k]) {
	      distance++;
	      if (first == -1) first = k;
	      last = k;
	    }
	  }
	  min_distance = (gap_genome / 3) + 2;
	  //printf("Distance %i <= %i\n", distance, min_distance);
	  if (distance <= min_distance) {
	    cigar_code->distance += distance;
	    closed = 1;
	    add_seed = 1;
	  }
	}

	if (!closed) {
	  int gap = gap_genome - gap_read; 
	  if (gap > 40) {
	    //printf("Search splice");
	    int distance_aux;
	    int nt = search_simple_splice_junction(seed_prev, seed_next,
						   cal->chromosome_id, cal->strand, 
						   query_map, genome, 
						   avls_list, &node_avl_start,
						   &node_avl_end, &distance_aux);
	    if (nt) {
	      cigar_code->distance = distance_aux;
	      metaexon_insert(cal->strand, cal->chromosome_id - 1,
			      seed_prev->genome_start, node_avl_start->position, 40,
			      METAEXON_RIGHT_END, node_avl_start,
			      metaexons);
	      
	      metaexon_insert(cal->strand, cal->chromosome_id - 1,
			      node_avl_end->position, seed_next->genome_end, 40,
			      METAEXON_LEFT_END, node_avl_end,
			      metaexons);
	      closed = 1;
	      sp_found = 1;
	      //TODO: ADD DISTANCE
	    }
	  }
	} 

	if (!closed) {
	  //SMITH-WATERMAN
	  if (gap_read <= 0 || gap_genome <= 0) {
	    seed_prev->genome_end   = seed_prev->genome_end - FLANK;
	    seed_next->genome_start = seed_next->genome_start + FLANK;
	    seed_prev->read_end   = seed_prev->read_end - FLANK;
	    seed_next->read_start = seed_next->read_start + FLANK;	    
	  }
	  size_t genome_start = seed_prev->genome_end + 1;
	  size_t genome_end   = seed_next->genome_start - 1;
	  int read_start      = seed_prev->read_end + 1;
	  int read_end        = seed_next->read_start - 1;

	  genome_read_sequence_by_chr_index(reference, 0, cal->chromosome_id - 1,
					    &genome_start, &genome_end, genome);
	  char query[2048];
	  memcpy(query, &query_map[read_start],  read_end - read_start + 1);
	  query[read_end - read_start + 1] = '\0';

	  //printf("query : %s [%i-%i]\n", query, read_start, read_end);
	  //printf("ref   : %s [%lu-%lu]\n", reference, genome_start, genome_end);
	  seed_region_t *new_seed = seed_region_new(seed_prev->read_end + 1, seed_next->read_start - 1, 
						    seed_prev->genome_end + 1, seed_next->genome_start - 1, 
						    seed_prev->id + 1);
	  new_seed->fusion_left  = 1;
	  new_seed->fusion_right = 1;
	  linked_list_item_t *new_item = linked_list_item_new(new_seed);
	  list_item_prev->next = new_item;
	  new_item->prev = list_item_prev;
	  list_item->prev = new_item;
	  new_item->next = list_item;	  
	  cal->sr_list->size++;

	  sw_item = sw_item_new(SIMPLE_SW, i, 0, 0,
				cal, cal, NULL, 
				new_seed, new_seed,
				NULL);
	
	  //Insert item... and process if depth is full
	  sw_depth_insert(query, reference, sw_item,
			  sw_optarg, output,
			  avls_list, metaexons, sw_depth);	    

	  cal->num_targets++;
	}

	if (add_seed) {
	  seed_region_t *new_seed = seed_region_new(seed_prev->read_end + 1, seed_next->read_start - 1, 
						    seed_prev->genome_end + 1, seed_next->genome_start - 1, 
						    seed_prev->id + 1);
	  new_seed->fusion_left  = 1;
	  new_seed->fusion_right = 1;
	  linked_list_item_t *new_item = linked_list_item_new(new_seed);
	  list_item_prev->next = new_item;
	  new_item->prev = list_item_prev;
	  list_item->prev = new_item;
	  new_item->next = list_item;	  
	  cal->sr_list->size++;
	}

	if (!sp_found) {
	  seed_prev->fusion_right = 1;
	  seed_next->fusion_left  = 1;
	}

      }
      seed_prev      = seed_next;
      list_item_prev = list_item;
      list_item      = list_item->next;
    }

    cal->fill_gaps = 1;
    //printf("#### (%i)[%i:%lu-%lu] : %i ####\n", cal->strand, 
    //	   cal->chromosome_id, cal->start, cal->end, cal->num_targets);

    if (cal->num_targets == 0 && cal->fill_gaps) {
      //Close CAL
      merge_seeds_cal(cal);
    }
    //printf("END CAL %i\n", i); 
  }
  //printf("EXIT LOOP\n");
}


meta_alignment_t *meta_alignment_cal_new(cal_t *cal) {
  meta_alignment_t *meta_alignment = meta_alignment_new();
  meta_alignment_insert_cal(cal, meta_alignment);

  return meta_alignment;
}

meta_alignment_t *meta_alignment_cals_new(array_list_t *cals_list) {
  meta_alignment_t *meta_alignment = meta_alignment_new();
  
  for (int i = 0; i < array_list_size(cals_list); i++) {
    cal_t *cal = array_list_get(i, cals_list);
    meta_alignment_insert_cal(cal, meta_alignment);
  }

  return meta_alignment;
}



//===============================================================================//


//============================= SMITH-WATERMAN SECTION ==========================//
info_sp_t *info_sp_new(size_t l_genome_start, size_t l_genome_end,
		       size_t r_genome_start, size_t r_genome_end) {

  info_sp_t *info_sp = (info_sp_t *)malloc(sizeof(info_sp_t));

  info_sp->l_genome_start = l_genome_start;
  info_sp->l_genome_end = l_genome_end;
  info_sp->r_genome_start = r_genome_start;
  info_sp->r_genome_end = r_genome_end;

  return info_sp;
}

void info_sp_free(info_sp_t *info_sp) {
  free(info_sp);
}

sw_item_t *sw_item_new(int type_sw, int read_id, 
		       int fusion_id, int cal_id,
		       cal_t *cal_prev, cal_t *cal_next, 
		       meta_alignment_t *meta_alignment, 
		       seed_region_t *seed_prev,
		       seed_region_t *seed_next,
		       void *info) {

  sw_item_t *sw_item = (sw_item_t *)malloc(sizeof(sw_item_t));
  
  sw_item->type_sw = type_sw;
  sw_item->read_id = read_id;
  sw_item->fusion_id = fusion_id;
  sw_item->cal_id = cal_id;
  sw_item->cal_prev = cal_prev;
  sw_item->cal_next = cal_next;
  sw_item->info = info;
  sw_item->meta_alignment = meta_alignment;
  sw_item->seed_prev = seed_prev;
  sw_item->seed_next = seed_next;

  return sw_item;

}

void sw_item_free(sw_item_t *sw_item) {
  free(sw_item);
}

void sw_depth_process(sw_optarg_t *sw_optarg, sw_multi_output_t *output, 
		      sw_depth_t *sw_depth, avls_list_t *avls_list,
		      metaexons_t *metaexons, int step) {
  int distance, len;
  float norm_score;
  float match = sw_optarg->subst_matrix['A']['A'];
  if (sw_depth->depth == MAX_DEPTH || 
      (step == SW_FINAL && sw_depth->depth > 0)) {

    smith_waterman_mqmr(sw_depth->q, sw_depth->r, sw_depth->depth, sw_optarg, 1, output);

    for (int i = 0; i < sw_depth->depth; i++) {
      sw_item_t *sw_item = sw_depth->items[i];     
      cal_t *cal_prev = sw_item->cal_prev;
      cal_t *cal_next = sw_item->cal_next;
      //printf("QUE: %s(%i)\n", output->query_map_p[i], output->query_start_p[i]);
      //printf("REF: %s(%i)\n", output->ref_map_p[i], output->ref_start_p[i]);
      
      if (sw_item->type_sw == EXTREM_SW_LEFT) {	
	norm_score = NORM_SCORE(output->score_p[i], strlen(sw_depth->q[i]), match);
	//printf("EXTREM SW LEFT SCORE %f\n", norm_score);
	cigar_code_t *cigar_code = generate_cigar_code(output->query_map_p[i], output->ref_map_p[i],
						       strlen(output->ref_map_p[i]),
						       output->query_start_p[i], output->ref_start_p[i],
						       strlen(sw_depth->q[i]), strlen(sw_depth->r[i]),
						       &distance, FIRST_SW);
	if (norm_score >= 0.4) {
	  //printf("....>%s\n", new_cigar_code_string(cigar_code));
	  meta_alignment_insert_cigar(cigar_code, CIGAR_ANCHOR_RIGHT, sw_item->cal_id, sw_item->meta_alignment); 
	} else {
	  meta_alignment_insert_cigar(NULL, CIGAR_ANCHOR_RIGHT, sw_item->cal_id, sw_item->meta_alignment); 
	}
      } else if (sw_item->type_sw == EXTREM_SW_RIGHT) {
	norm_score = NORM_SCORE(output->score_p[i], strlen(sw_depth->q[i]), match);
	//printf("EXTREM SW LEFT SCORE %f\n", norm_score);
	cigar_code_t *cigar_code = generate_cigar_code(output->query_map_p[i], output->ref_map_p[i],
						       strlen(output->ref_map_p[i]),
						       output->query_start_p[i], output->ref_start_p[i],
						       strlen(sw_depth->q[i]), strlen(sw_depth->r[i]),
						       &distance, LAST_SW);	
	if (norm_score >= 0.4) {
	  meta_alignment_insert_cigar(cigar_code, CIGAR_ANCHOR_LEFT, sw_item->cal_id, sw_item->meta_alignment); 
	} else {
	  meta_alignment_insert_cigar(NULL, CIGAR_ANCHOR_LEFT, sw_item->cal_id, sw_item->meta_alignment); 
	}

      } else if (sw_item->type_sw == SIMPLE_SW) {
	cigar_code_t *cigar_code = generate_cigar_code(output->query_map_p[i], output->ref_map_p[i],
						       strlen(output->ref_map_p[i]),
						       output->query_start_p[i], output->ref_start_p[i],
						       strlen(sw_depth->q[i]), strlen(sw_depth->r[i]),
						       &distance, MIDDLE_SW);
	cal_prev->num_targets--;
	seed_region_t *seed_prev = sw_item->seed_prev;
	seed_prev->info = cigar_code;
	
	//printf("SW CIGAR %s\n", new_cigar_code_string(cigar_code));
	//printf("(%i)[%i:%lu-%lu].....MERGE SEEDS %i??\n", cal_prev->strand, cal_prev->chromosome_id, 
	//     cal_prev->start, cal_prev->end, cal_prev->num_targets);
	if (cal_prev->num_targets == 0 && cal_prev->fill_gaps) {
	  //Close CAL
	  //printf("****MERGE SEEDS\n");
	  merge_seeds_cal(cal_prev);
	}
      } else {
	//fastq_read_t *read = array_list_get(sw_item->read_id, fq_batch);
	//if (read == NULL) { printf("@@@@@(%i)@ %s\n", sw_item->read_id, read->id); exit(-1); }
	info_sp_t *info_sp = sw_item->info;
	avl_node_t *node_avl_start, *node_avl_end;
	norm_score = NORM_SCORE(output->score_p[i], strlen(sw_depth->q[i]), match);
	//printf("QUE: %s(%i)\n", output->query_map_p[i], output->query_start_p[i]);
	//printf("REF: %s(%i)\n", output->ref_map_p[i], output->ref_start_p[i]);
	//printf("%f\n", norm_score);
	cigar_code_t *cigar_code;
	if (norm_score >= 0.6) {
	  cigar_code = generate_cigar_sw_output(output->query_map_p[i], 
						output->ref_map_p[i],
						info_sp->l_genome_start,
						info_sp->l_genome_end,
						info_sp->r_genome_start,
						info_sp->r_genome_end,
						sw_item->cal_prev->chromosome_id,
						sw_item->cal_prev->strand,
						output->query_start_p[i],
						output->ref_start_p[i],
						strlen(sw_depth->q[i]),
						strlen(sw_depth->r[i]),
						avls_list,
						&node_avl_start,
						&node_avl_end);
	  //printf("CIGAR SW : %s\n", new_cigar_code_string(cigar_code));
	  if (cigar_code != NULL) {
	    if (node_avl_start  != NULL) {
	      metaexon_insert(cal_prev->strand, cal_prev->chromosome_id - 1,
			      info_sp->l_genome_start, node_avl_start->position, 40,
			      METAEXON_RIGHT_END, node_avl_start,
			      metaexons);
	      
	      metaexon_insert(cal_next->strand, cal_next->chromosome_id - 1,
			    node_avl_end->position, info_sp->r_genome_end, 40,
			      METAEXON_LEFT_END, node_avl_end,
			      metaexons);
	    }
	    //printf("Actualize splice junction... %i\n", cigar_code->distance);
	    //cigar_code_print(cigar_code);
	  }
	} else {
	  //array_list_clear(cigar_code->ops, cigar_op_free);
	  //cigar_code_free(cigar_code);
	  cigar_code = NULL;
	}
	//TODO: REPORT READS WITH SMALL-MIDDLE-EXON
	/*
	if (cigar_code == NULL) {
	  seed_region_t *s_prev = linked_list_get_last(sw_item->cal_prev->sr_list);
	  seed_region_t *s_next = linked_list_get_first(sw_item->cal_next->sr_list);
	  size_t read_nt = s_next->read_start - s_prev->read_end - 1;
	  size_t genome_nt  = (s_next->genome_start) - (s_prev->genome_end + read_nt) - 1;
	  if (genome_nt > 40) {	    
	    cigar_code = cigar_code_new();
	    cigar_code_append_new_op(read_nt, 'I', cigar_code);
	    cigar_code_append_new_op(genome_nt, 'N', cigar_code);
	  }
	}*/
	info_sp_free(info_sp);
	meta_alignment_insert_cigar(cigar_code, CIGAR_SW_MIDDLE, sw_item->cal_id, sw_item->meta_alignment);	  		
      }
      free(output->query_map_p[i]);
      free(output->ref_map_p[i]);
      output->query_map_p[i] = NULL;
      output->ref_map_p[i] = NULL;
      free(sw_depth->q[i]);
      free(sw_depth->r[i]);
      sw_item_free(sw_item);
    }
    sw_depth->depth = 0;
  }
}

void sw_depth_insert(char *query, char *reference, sw_item_t *sw_item, 
		     sw_optarg_t *sw_optarg, sw_multi_output_t *output, 
		     avls_list_t *avls_list, metaexons_t *metaexons, 
		     sw_depth_t *sw_depth) {

  sw_depth->q[sw_depth->depth]       = strdup(query);
  sw_depth->r[sw_depth->depth]       = strdup(reference);
  sw_depth->items[sw_depth->depth++] = sw_item;

  //printf("\tInsert SW depth %i\n", sw_depth->depth);
  //printf("QUE: %s\n", query);
  //printf("REF: %s\n", reference);

  sw_depth_process(sw_optarg, output, 
		   sw_depth, avls_list, metaexons, SW_NORMAL);

}
//===============================================================================//


info_sp_t* sw_reference_splice_junction(cal_t *cal_prev, cal_t *cal_next,
					char *query_map, genome_t *genome,
					char *q, char *r) {

  //printf("============= M O U N T    S M I T H - W A T E R M A N =============\n");  
  //printf(" ==== CALS INFO ==== \n");
  //cal_print(cal_prev);
  //cal_print(cal_next);
  //printf(" ==== CALS INFO END ==== \n");

  seed_region_t *s_next, *s_prev;
  char reference[2048];
  char reference_prev[2048];
  char reference_next[2048];
  char query[2048];

  size_t genome_start, genome_end;
  int read_start, read_end, read_gap;
  size_t genome_start2, genome_end2;
  int seeds_nt;
  int flank = 30;
  int flank_left, flank_right;

  //Delete for debuging, detect splice junctions
  s_prev = linked_list_get_last(cal_prev->sr_list);    
  s_next = linked_list_get_first(cal_next->sr_list);

  read_start = s_prev->read_end;
  read_end = s_next->read_start;
  seeds_nt = read_end - read_start;

  //printf("SP coords %i - %i = %i\n", read_end, read_start, seeds_nt);

  if (seeds_nt >= 10) {
    //Extend to Right --> <-- Extend to Left
    genome_start = s_prev->genome_end;
    genome_end = s_prev->genome_end + seeds_nt - 1;
    genome_read_sequence_by_chr_index(reference_prev, 0, 
				      cal_prev->chromosome_id - 1, &genome_start, &genome_end, genome);

    genome_start2 = s_next->genome_start - seeds_nt;
    genome_end2 = s_next->genome_start - 1;
    genome_read_sequence_by_chr_index(reference_next, 0, 
				      cal_next->chromosome_id - 1, &genome_start2, &genome_end2, genome);

    memcpy(query, query_map + read_start, read_end - read_start);
    query[read_end - read_start]  = '\0';

    //CAll new function
    int dsp_e1, dsp_e2;
    int lim_err = 2;

    extend_by_mismatches(reference_prev, reference_next, query, 
			 0, strlen(reference_next) - 1, 
			 0, read_end - read_start - 1, lim_err,
			 &dsp_e1, &dsp_e2);
    //printf("dsp_e1=%i, dsp_e2=%i\n", dsp_e1, dsp_e2);

    cal_prev->end = cal_prev->end + dsp_e1;
    cal_next->start = cal_next->start - dsp_e2;

    s_prev->read_end = s_prev->read_end + dsp_e1;
    s_prev->genome_end = s_prev->genome_end + dsp_e1;

    s_next->read_start = s_next->read_start - dsp_e2;
    s_next->genome_start = s_next->genome_start - dsp_e2;

    if (s_prev->read_end > s_next->read_start) { 
      seeds_nt = s_prev->read_end - s_next->read_start;
      cal_prev->end -= seeds_nt;
      cal_next->start += seeds_nt;
      s_prev->read_end-= seeds_nt;
      s_next->read_start += seeds_nt;
    }
  }
    
  //printf("Read_start =: s_prev->read_end = %i - %i + 1 = %i\n", s_prev->read_end, flank, s_prev->read_end - flank + 1);
  
  int seed_size = s_prev->read_end - s_prev->read_start + 1;  
  flank_left = flank;
  if (flank > seed_size) {
    flank_left = seed_size;    
  } 
  read_start = s_prev->read_end - flank_left + 1;
  
  //printf("Read_end =: s_next->read_start = %i + %i - 1 = %i\n", s_next->read_start, flank, s_next->read_start + flank - 1);
  seed_size = s_next->read_end - s_next->read_start + 1;  
  flank_right = flank;
  if (flank > seed_size) {
    flank_right = seed_size;    
  } 
  read_end = s_next->read_start + flank_right - 1;
  

  //Extract and fusion Reference SW
  cal_prev->l_flank = flank_left;
  //cal_prev->r_flank = flank;
  genome_start = cal_prev->end - flank_left + 1;
  //printf("GENOME END %lu - %i + 1 = %lu\n", cal_prev->end, flank_left, genome_start);
  genome_end = cal_prev->end + flank - 1;
  genome_read_sequence_by_chr_index(reference_prev, 0, 
				    cal_prev->chromosome_id - 1, &genome_start, &genome_end, genome);
  //printf("1g[From %lu to %lu](%i): %s(%i)\n", genome_start, genome_end, genome_end - genome_start + 1, 
  //	 reference_prev, strlen(reference_prev));

  //cal_next->l_flank = flank;
  cal_next->r_flank = flank_right;
  genome_start2 = cal_next->start - flank;
  genome_end2 = cal_next->start + flank_right - 1;
  genome_read_sequence_by_chr_index(reference_next, 0, 
				    cal_next->chromosome_id - 1, &genome_start2, &genome_end2, genome);
  //printf("2g[From %lu to %lu](%i): %s(%i)\n", genome_start2, genome_end2, genome_end2 - genome_start2 + 1, 
  //	 reference_next, strlen(reference_next));

  strcat(reference_prev, reference_next);

  if (read_start > read_end) { LOG_FATAL("READ COORDS ERROR\n"); }
  //printf("Read %i-%i: %s\n", read_start, read_end, query_map);

  read_gap = read_end - read_start + 1;
  memcpy(query, &query_map[read_start], read_gap);
  query[read_gap] = '\0';

  //CALs Actualization flank
  //cal_prev->end -= flank;
  //cal_next->start += flank;
  //s_prev->read_end-= flank;
  //s_next->read_start += flank;
  //printf("\tread_start = %i, read_end = %i, genome_start = %lu, genome_end = %lu, flank_left=%i, flank_right=%i, query=%i =? read_gap=%i\n", 
  //	 read_start, read_end, genome_start, genome_end,  
  //	 flank_left, flank_right, strlen(query), read_gap);

  strcpy(q, query);
  strcpy(r, reference_prev);

  //printf("============= M O U N T    S M I T H - W A T E R M A N    E N D =============\n");

  return info_sp_new(genome_start, genome_end,
		     genome_start2, genome_end2);

}

//====================================================================================//


//FOR SEARCH CANNONICAL SPLICE JUNCTION WITHOUT SMITH-WATERMAN ALGORITHM
int search_simple_splice_junction(seed_region_t *s_prev, seed_region_t *s_next,
				  int chromosome_id, int strand, 
				  char *sequence, genome_t *genome, 
				  avls_list_t *avls_list, 
				  avl_node_t **node_avl_start,
				  avl_node_t **node_avl_end,
				  int *distance) {

  assert(s_prev != NULL);
  assert(s_next != NULL);

  int read_start   = s_prev->read_end;
  int read_end     = s_next->read_start;

  //printf("START READ START %i/%lu READ END %i/%lu\n", read_start, s_prev->genome_end, 
  //	 read_end, s_next->genome_start);

  size_t genome_start;
  size_t genome_end;
  
  const int FLANK = 20;
  const int SEQURITY_FLANK = 5;

  int gap_read = read_end - read_start - 1;
  if (gap_read == 0) {
    gap_read = -1;
  }
  
  int read_end_aux = s_prev->read_end;
  int read_start_aux = s_next->read_start;
  size_t genome_start_aux = s_next->genome_start;
  size_t genome_end_aux = s_prev->genome_end;

  LOG_DEBUG_F("SEARCH | %i-%i | %lu-%lu | \n", read_end_aux, read_start_aux, genome_end_aux, genome_start_aux);
  if (gap_read < 0) {
    gap_read = abs(gap_read) + 5;
    read_end_aux     -= gap_read;
    read_start_aux   += gap_read;
    genome_end_aux   -= gap_read;
    genome_start_aux += gap_read;    
  } else {
    read_end_aux     -= SEQURITY_FLANK;
    read_start_aux   += SEQURITY_FLANK;
    genome_end_aux   -= SEQURITY_FLANK;
    genome_start_aux += SEQURITY_FLANK;
  }

  read_start = read_end_aux;
  read_end = read_start_aux;
  gap_read = read_end - read_start - 1;

  char left_exon[2048];
  char right_exon[2048];

  genome_start = genome_end_aux + 1;
  genome_end   = genome_end_aux + gap_read + FLANK;

  LOG_DEBUG_F("GAP READ %i - %i = %i\n", read_end, read_start, gap_read);
  LOG_DEBUG_F("SEQUENCE   : %s\n", sequence);

  genome_read_sequence_by_chr_index(left_exon, 0, 
				    chromosome_id - 1, 
				    &genome_start, &genome_end, genome);		  

  LOG_DEBUG_F("LEFT EXON  (%lu-%lu): %s\n", genome_start, genome_end, left_exon);

  genome_start = genome_start_aux - gap_read - FLANK;
  genome_end   = genome_start_aux - 1;
  
  genome_read_sequence_by_chr_index(right_exon, 0, 
				    chromosome_id - 1, 
				    &genome_start, &genome_end, genome);		  

  LOG_DEBUG_F("RIGHT EXON (%lu-%lu): %s\n", genome_start, genome_end, right_exon);

  int dsp_l, dsp_r, type;

  int breaks_starts[gap_read];
  int found_starts = 0;

  int breaks_ends[gap_read];
  int found_ends = 0;

  int type_starts[gap_read];
  int type_ends[gap_read];
  int c_s, c_e;
  int end_search = gap_read + SEQURITY_FLANK;

  //printf("search start!\n");
  // Search step by step (GT)/(AG) 
  for (c_s = 0, c_e = strlen(right_exon) - 1;
       c_s < end_search; c_s++, c_e--) {
    if (left_exon[c_s] == 'G' && left_exon[c_s + 1] == 'T') {
      type_starts[found_starts] = GT_AG_SPLICE;
      breaks_starts[found_starts++] = c_s;
      LOG_DEBUG_F("FOUND GT (%i)\n", c_s);
    } else if (left_exon[c_s] == 'C' && left_exon[c_s + 1] == 'T') {
      type_starts[found_starts] = CT_AC_SPLICE;
      breaks_starts[found_starts++] = c_s;
      LOG_DEBUG_F("FOUND CT (%i)\n", c_s);
    }

    if (right_exon[c_e] == 'G' && right_exon[c_e - 1] == 'A') {
      type_ends[found_ends] = GT_AG_SPLICE;
      breaks_ends[found_ends++] = strlen(right_exon) - c_e - 1;
      LOG_DEBUG_F("FOUND AG (%i)\n", strlen(right_exon) - c_e - 1);
    } else if (right_exon[c_e] == 'C' && right_exon[c_e - 1] == 'A') {
      type_ends[found_ends] = CT_AC_SPLICE;
      breaks_ends[found_ends++] = strlen(right_exon) - c_e - 1;
      LOG_DEBUG_F("FOUND AC (%i)\n", strlen(right_exon) - c_e - 1);
    }
  }

  //Not found any splice junction
  if (found_starts == 0 || found_ends == 0) {
    return 0;
  } 

  array_list_t *splice_junction = array_list_new(found_starts + found_ends, 
						 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  array_list_t *splice_junction_1 = array_list_new(found_starts + found_ends, 
						   1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  //printf("FOUND STARTS %i, FOUND ENDS %i\n", found_starts, found_ends);

  //If found more than one break...
  for (int i = 0; i < found_starts; i++) {
    for (int j = 0; j < found_ends; j++) {
      //printf("%i(%i) + %i(%i)=%i\n", breaks_starts[i], type_starts[i],
      //     breaks_ends[j], type_ends[j], breaks_starts[i] + breaks_ends[j]);
      if (type_starts[i] == type_ends[j]) {
	int gap_break = breaks_starts[i] + breaks_ends[j];
	if (gap_break == gap_read) {
	  array_list_insert(breaks_starts[i], splice_junction);
	  array_list_insert(breaks_ends[j], splice_junction);
	  array_list_insert(type_ends[j], splice_junction);
	}      
      }
    }
  }

  if (array_list_size(splice_junction) <= 0) {
    array_list_free(splice_junction, NULL);
    return 0;
  }

  //Calculating scores...
  int matches[array_list_size(splice_junction)];
  int mismatches[array_list_size(splice_junction)];
  float max_score = 0.0f;
  float score;
  int max_sp_pos;
  int sp_pos;
  int read_pos;
  int genome_pos;

  //printf("NUM SP (%i)\n", array_list_size(splice_junction)/3);
  for (int i = 0; i < array_list_size(splice_junction); i += 3) {
    int limit_left = array_list_get(i, splice_junction);
    int limit_right = array_list_get(i + 1, splice_junction);      

    sp_pos = 0;
      
    matches[sp_pos] = 0;
    mismatches[sp_pos] = 0;

    read_pos = read_start + 1;
    for (int c_l = 0; c_l < limit_left; c_l++) {
      //printf("l: %c == %c\n", left_exon[c_l], sequence[c_l]);
      if (left_exon[c_l] == sequence[read_pos++]) {
	matches[sp_pos]++; 
      } else {
	mismatches[sp_pos]++; 
      }
    }

    read_pos = read_end - 1;
    genome_pos = strlen(right_exon) - 1;

    for (int c_r = 0; c_r < limit_right; c_r++) {
      //printf("r: %c == %c\n", right_exon[genome_pos], sequence[read_pos]);
      if (right_exon[genome_pos--] == sequence[read_pos--]) {
	matches[sp_pos]++;
      } else {
	mismatches[sp_pos]++;
      }
    }

    score = matches[sp_pos]*0.5 - mismatches[sp_pos]*0.4;
    if (score > max_score) {
      max_sp_pos = sp_pos;
      max_score = score;
    }
   
    //printf("MATCHES %i / MISMATCHES %i\n", matches[sp_pos], mismatches[sp_pos]);
    sp_pos++;
  }

  dsp_l = array_list_get((max_sp_pos * 3), splice_junction);
  dsp_r = array_list_get((max_sp_pos * 3) + 1, splice_junction);
  type  = array_list_get((max_sp_pos * 3) + 2, splice_junction);

  *distance = mismatches[max_sp_pos];
  
  if (type == CT_AC_SPLICE) {
    strand = 1;
  } else {
    strand = 0;
  }

  //assert(array_list_size(splice_junction) != 0);

  //TODO: CALCULATE DISTANCE 
  //printf("dsp_l = %i, dsp_r = %i, read_end_aux = %i, read_start_aux = %i\n", dsp_l, dsp_r, read_end_aux, read_start_aux);

  read_end_aux     += dsp_l;
  read_start_aux   -= dsp_r;
  genome_end_aux   += dsp_l;
  genome_start_aux -= dsp_r;
  
  s_prev->read_end     = read_end_aux;
  s_next->read_start   = read_start_aux;
  s_next->genome_start = genome_start_aux;
  s_prev->genome_end   = genome_end_aux;
  
  size_t start_splice = s_prev->genome_end + 1;
  size_t end_splice   = s_next->genome_start - 1;


  allocate_start_node(chromosome_id - 1,
		      strand,
		      start_splice,
		      end_splice,
		      start_splice,
		      end_splice,
		      FROM_READ,
		      type,
		      NULL, 
		      node_avl_start,
		      node_avl_end, 
		      avls_list);

  //printf("SP :=> [%i:%lu-%lu]\n", chromosome_id, start_splice, end_splice);

  array_list_free(splice_junction, NULL);

  return end_splice - start_splice + 1;

}

cigar_code_t *search_left_single_anchor(int gap_close, 
					cal_t *cal,
					int filter_pos, 
					array_list_t *right_breaks,
					char *query_map,
					metaexons_t *metaexons,
					genome_t *genome) {
  
  int max_nt;
  avl_node_t *node_start_prev, *node_start_next;
  int dist, final_dist = 0;
  int abort = 0, map = 0;
 
  array_list_t *final_positions = array_list_new(50, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *starts_targets  = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *final_starts    = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
 
  metaexon_t *final_metaexon;
  size_t final_pos;
  char reference[2048];
  size_t genome_start, genome_end;
  avl_node_t *node_start;
  int max_dist;
  int s;
  seed_region_t *s_prev = linked_list_get_last(cal->sr_list);

  int read_pos = s_prev->read_end + 1;
  size_t first_cal_end = cal->end;
  int cal_strand = cal->strand;
  int cal_chromosome_id = cal->chromosome_id;

  genome_start = first_cal_end + 1;
  //printf("IN PARAMETERS: read_pos->%i, first_cal_end->%lu, gap_close->%i\n", read_pos, first_cal_end, gap_close, genome_start );
  //array_list_insert(genome_start, final_positions);
  //==== 1st-Select the correct start ====//
  for (int s_0 = 0; s_0 < array_list_size(right_breaks); s_0++) {
    node_start = array_list_get(s_0, right_breaks);
    if (first_cal_end - 10 <= node_start->position) {
      LOG_DEBUG_F("\tNode start %lu\n", node_start->position);
      array_list_insert(node_start, starts_targets);
    }
  }
  
  while (gap_close > 0) {
    //1. Order starts
    int num_targets = array_list_size(starts_targets); 
    if (!num_targets) { break; }
    for (int s0 = 0; s0 < num_targets - 1; s0++) {
      node_start_prev = array_list_get(s0, starts_targets);
      for (int s1 = s0 + 1; s1 <  num_targets; s1++) {
	node_start_next = array_list_get(s1, starts_targets);
	if (node_start_next->position < node_start_prev->position) {
	  array_list_swap(s0, s1, starts_targets);
	}
      } 
    }

    //2. Map the gap to the genome and Select the correct start node
    node_start = array_list_get(num_targets - 1, starts_targets);
    genome_end = node_start->position - 1;
    int lim_ref;
    lim_ref = genome_end - genome_start + 1; 
    //printf("genome_end = %lu, genome_start = %lu, lim_ref = %i\n", genome_end, genome_start, lim_ref);
    if (lim_ref < 0) {
      //For first step, genome_start = first_cal->end
      int dsp = abs(lim_ref);
      gap_close += dsp;
      read_pos  -= dsp;
      array_list_insert(genome_end, final_positions);
      //printf("1):::::::::::: INSERT FINAL POSITIONS: %i, (%i)\n", genome_end, array_list_size(final_positions));
      //printf(":::--::::GENOME GAP %i\n", genome_end);
      //printf("RECALCULATING GAP AND READ POS (gap_close, read_pos)(%i, %i)\n", gap_close, read_pos);
    } else {      
      genome_read_sequence_by_chr_index(reference, 0, 
					cal_chromosome_id - 1, 
					&genome_start, &genome_end, genome);
      //printf("Reference START_SP [[[START_SP]]]----[END_SP]: %s\n", reference);
      max_dist = 3;
      int t, c;
      for (t = 0; t < num_targets; t++) {
	node_start = array_list_get(t, starts_targets);
	lim_ref = node_start->position - genome_start;
	if (lim_ref > gap_close) { lim_ref = gap_close; }
	dist = 0;
	//printf("\tTravel to %lu\n", node_start->position);
	for (c  = 0; c < lim_ref; c++) { 
	  //printf("\t\t[%c vs %c]\n", query_map[read_pos + c], reference[c]);
	  if (query_map[read_pos + c] != reference[c]) { 
	    dist++;
	  }
	}
	if (dist < max_dist) { 
	  //printf("\t\tMAX DISTANCE GOOD %i!\n", dist);
	  max_dist += dist;
	  final_dist += dist;
	} else { 
	  //printf("\t\tExit LOOP %i\n", t - 1);
	  break;
	}
      }
      t--;
      if (t < 0) { break; }
      else {
	node_start = array_list_get(t, starts_targets);
	lim_ref = node_start->position - genome_start;
	gap_close -= lim_ref;
	read_pos += lim_ref;	
	//printf("GAP CLOSE (%i) READ POS (%i)\n", gap_close, read_pos);
	if (gap_close <= 0) {
	  //Final Map
	  map = 1;
	  array_list_insert(genome_start + c - 1, final_positions);
	  //printf("3):::::::::::: INSERT FINAL POSITIONS: %i (%i)\n", genome_start + c, array_list_size(final_positions));
	  break;
	} else {
	  array_list_insert(node_start->position - 1, final_positions);
	  //printf("2):::::::::::: INSERT FINAL POSITIONS: %i (%i)\n", node_start->position - 1, array_list_size(final_positions));
	}
      }
    }
	    
    //3.Select the correct start
    int pos;
    array_list_t *ends_list = ((start_data_t *)node_start->data)->list_ends;

    array_list_clear(final_starts, NULL);
    //Filter ends before cal_next->start
    for (int e = 0; e < array_list_size(ends_list); e++) {
      splice_end_t *splice_end = array_list_get(e, ends_list);
      size_t start = splice_end->end;
      if (filter_pos) {
	if (start < filter_pos) {
	  array_list_insert(start, final_starts);
	}
      } else {
	array_list_insert(start, final_starts);
      }
    }
    if (array_list_size(final_starts) >= 1) {
      //Select best start
      char s_reference[array_list_size(final_starts)][2048];
      //Making reference...
      for (int s = 0; s < array_list_size(final_starts); s++) {
	genome_start = array_list_get(s, final_starts) + 1;
	genome_end = genome_start + gap_close + 1;
	genome_read_sequence_by_chr_index(s_reference[s], 0, 
					  cal_chromosome_id - 1, 
					  &genome_start, &genome_end, genome); 	      
	//printf("Reference END_SP [START_SP]----[[[END_SP]]]  [%i:%lu-%lu](%i): %s\n", cal_chromosome_id, genome_start, genome_end, 
	//   s, s_reference[s]);
      } 
      //Select the best start

      int anchor_nt = 20; 
      if (anchor_nt > gap_close) { 
	anchor_nt = gap_close;
      }
      int *ref_matches = (int *)calloc(array_list_size(final_starts), sizeof(int));
      int *ref_mismatches = (int *)calloc(array_list_size(final_starts), sizeof(int));
      //printf("anchor nt = %i, array_list_size(final_starts)=%i\n", anchor_nt, array_list_size(final_starts));
      for (int c = 0; c < anchor_nt; c++) {
	for (int s = 0; s < array_list_size(final_starts); s++) {
	  if (read_pos + c > strlen(query_map)) { exit(-1); }
	  if (c > strlen(s_reference[s])) { exit(-1); }
	  //printf("\t(%i )[%c vs %c]\n", s, query_map[read_pos + c], s_reference[s][c]);
	  if (query_map[read_pos + c] != s_reference[s][c]) {
	    ref_mismatches[s]++;
	  } else {
	    ref_matches[s]++;
	  }
	}
      }
      pos = 0;
      float score, max_score = 0.0;
      for (int s = 0; s < array_list_size(final_starts); s++) {
	score = ref_matches[s]*0.5 - ref_mismatches[s]*0.4;
	if (score > max_score) { 
	  max_score = score;
	  pos = s;
	}	
      }		
      max_dist = 3;
      if (ref_mismatches[pos] > max_dist) {
	break;
      }

      final_dist += ref_mismatches[pos]; 
      final_pos = array_list_get(pos, final_starts);
      array_list_insert(final_pos, final_positions);
      //printf("4):::::::::::: INSERT FINAL POSITIONS: %i (%i)\n", final_pos, array_list_size(final_positions));
      genome_start = final_pos + anchor_nt + 1;
      gap_close -= anchor_nt;
      read_pos += anchor_nt;
      //printf("FINAL MATCHES (%i) VS FINAL MISMATCHES (%i)\n", ref_matches[s], ref_mismatches[s]);
      //printf("GAP CLOSE (%i) READ POS (%i)\n", gap_close, read_pos);

      free(ref_mismatches);
      free(ref_matches);

      if (anchor_nt < 20) {
	array_list_insert(final_pos + anchor_nt, final_positions);
	//printf("5):::::::::::: INSERT FINAL POSITIONS: %i (%i)\n", final_pos + anchor_nt, array_list_size(final_positions));
	map = 1;
	break; 
      } else {
	if (metaexon_search(cal_strand, cal_chromosome_id - 1, 
			    final_pos, final_pos + anchor_nt, 
			    &final_metaexon, metaexons)) {
	  if (final_metaexon) {
	    if (final_metaexon->right_closed) {
	      array_list_clear(starts_targets, NULL);
	      //printf("FOUND!!\n");
	      for (int s_0 = 0; s_0 < array_list_size(final_metaexon->right_breaks); s_0++) {
		node_start = array_list_get(s_0, final_metaexon->right_breaks);
		if (final_pos <= node_start->position) {
		  array_list_insert(node_start, starts_targets);		    
		}
	      }
	    } else {
	      //We can close the gap
	      int exon_length = final_metaexon->end - final_metaexon->start;
	      dist = 0;
	      if (exon_length >= gap_close) {
		genome_end = genome_start + gap_close + 1;
		genome_read_sequence_by_chr_index(reference, 0, 
						  cal_chromosome_id - 1, 
						  &genome_start, &genome_end, genome);		
		//printf("CLOSE GAP: %s\n", reference);
		for (int c  = 0; c < gap_close; c++) { 
		  //printf("\t[%c vs %c]\n", query_map[read_pos + c], reference[c]);
		  if (query_map[read_pos + c] != reference[c]) { 
		    dist++;
		  }
		}
		if (dist < 5) {
		  array_list_insert(final_pos + gap_close + anchor_nt, final_positions);
		  //printf("6):::::::::::: INSERT FINAL POSITIONS: %i (%i)\n", final_pos + gap_close + anchor_nt,
		  //	 array_list_size(final_positions));
		  final_dist += dist; 
		  map = 1;
		}
	      }
	      break;
	    }
	  }
	} else {
	  //printf("NOT  FOUND\n");
	  break;
	}
      }
    } else {
      break;
    }
  }

  cigar_code_t *cigar_code = NULL;

  if (map) {
    size_t pos_prev = cal->start, pos_next;
    cigar_op_t *op;
    cigar_code = cigar_code_new();
    cigar_code->distance = final_dist;

    pos_prev = cal->end + 1;
    //printf("FINAL POSITION %i\n", array_list_size(final_positions));
    for (int sp = 0; sp < array_list_size(final_positions); sp++) {
      pos_next = array_list_get(sp, final_positions);
      if (sp % 2 == 0) {		
	//printf("%lu - %lu(%i/%i)\n", pos_prev, pos_next, sp, array_list_size(final_positions));
	op = cigar_op_new(pos_next - pos_prev + 1, 'M');
      } else {
	int aux = pos_next - pos_prev + 1;
	op = cigar_op_new(aux, 'N');
      }
      //printf("\tADD %i%c\n", op->number, op->name);
      cigar_code_append_op(op, cigar_code);
      pos_prev = pos_next + 1;
    }
  }

  array_list_free(final_positions, NULL);
  array_list_free(starts_targets, NULL);
  array_list_free(final_starts, NULL);

  return cigar_code;

}

cigar_code_t *search_right_single_anchor(int gap_close, 
					 cal_t *cal,
					 int filter_pos, 
					 array_list_t *left_breaks,
					 char *query_map, metaexons_t *metaexons,
					 genome_t *genome) {
  int max_nt;
  avl_node_t *node_end_prev, *node_end_next;
  int dist, final_dist = 0;
  int abort = 0, map = 0;
 
  array_list_t *final_positions = array_list_new(50, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *ends_targets = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *final_ends   = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
 
  metaexon_t *final_metaexon;
  size_t final_pos;
  char reference[2048];
  size_t genome_start, genome_end;
  avl_node_t *node_end;
  int max_dist;
  int s;
  size_t reference_len;

  seed_region_t *s_prev = linked_list_get_first(cal->sr_list);
  int read_pos = s_prev->read_start - 1; 
  size_t last_cal_start =  cal->start;
  int cal_strand = cal->strand;
  int cal_chromosome_id = cal->chromosome_id;

  //printf("IN PARAMETERS: read_pos->%i, first_cal_end->%lu, gap_close->%i\n", read_pos, last_cal_start - 1, gap_close );

  genome_end = last_cal_start -  1;
  //array_list_insert(genome_start, final_positions);
  //==== 1st-Select the correct start ====//
  for (int s_0 = 0; s_0 < array_list_size(left_breaks); s_0++) {
    node_end = array_list_get(s_0, left_breaks);
    if (last_cal_start + 10 >= node_end->position) {
      //printf("\tNode end %lu\n", node_end->position);
      array_list_insert(node_end, ends_targets);
    }
  }
  
  while (gap_close > 0) {
    //1. Order starts
    int num_targets = array_list_size(ends_targets); 
    if (!num_targets) { break; }
    for (int s0 = 0; s0 < num_targets - 1; s0++) {
      node_end_prev = array_list_get(s0, ends_targets);
      for (int s1 = s0 + 1; s1 <  num_targets; s1++) {
	node_end_next = array_list_get(s1, ends_targets);
	if (node_end_next->position > node_end_prev->position) {
	  array_list_swap(s0, s1, ends_targets);
	}
      } 
    }

    //2. Map the gap to the genome and Select the correct start node
    node_end = array_list_get(num_targets - 1, ends_targets);
    genome_start = node_end->position;
    int lim_ref;
    lim_ref = genome_end - genome_start;	    
    //printf("lim_ref(%i) = genome_end(%lu) - genome_start(%lu) | gap_close = %i\n", lim_ref, genome_end, genome_start, gap_close);
    if (lim_ref < 0) {
      //For first step, genome_start = first_cal->end
      int dsp = abs(lim_ref);
      gap_close += dsp;
      read_pos  += dsp;
      array_list_insert(genome_start + 1, final_positions);
      //printf(" Insert 1) %lu\n", genome_start + 1);
      //printf(":::--::::GENOME GAP %i\n", genome_end);
      //printf("RECALCULATING GAP AND READ POS (%i) (gap_close, read_pos)(%i, %i)\n", dsp, gap_close, read_pos);
    } else {
      genome_read_sequence_by_chr_index(reference, 0, 
					cal_chromosome_id - 1, 
					&genome_start, &genome_end, genome);
      //printf("(%lu-%lu)Reference START_SP [[[START_SP]]]----[END_SP]: %s\n",  genome_start, genome_end, reference);
      max_dist = 3;
      int t, c;
      for (t = 0; t < num_targets; t++) {
	node_end = array_list_get(t, ends_targets);
	lim_ref = genome_end - node_end->position;
	//printf("???What?? lim_ref = %i, gap_close = %i\n", lim_ref, gap_close);
	if (lim_ref > gap_close) { lim_ref = gap_close; }
	dist = 0;
	//printf("\tTravel to %lu\n", node_end->position);
	reference_len = strlen(reference) - 1;
	for (c  = 0; c < lim_ref; c++) { 
	  //printf("\t\t[%c vs %c]\n", query_map[read_pos - c], reference[reference_len - c]);

	  if ((read_pos - c) < 0 || (reference_len - c) < 0) { exit(-1); }

	  if (query_map[read_pos - c] != reference[reference_len - c]) { 
	    dist++;
	  }
	}
	if (dist < max_dist) { 
	  //printf("\t\tMAX DISTANCE GOOD %i!\n", dist);
	  max_dist += dist;
	  final_dist += dist;
	} else { 
	  //printf("\t\tExit LOOP %i\n", t - 1);
	  break;
	}
      }
      t--;
      if (t < 0) { break; }
      else {
	node_end = array_list_get(t, ends_targets);
	lim_ref = genome_end - node_end->position;
	gap_close -= lim_ref;
	read_pos -= lim_ref;
	//printf("GAP CLOSE (%i) READ POS (%i)\n", gap_close, read_pos);
	if (gap_close <= 0) {
	  //Final Map
	  map = 1;
	  array_list_insert(genome_end - c + 1, final_positions);
	  //printf(" Insert 2) %lu\n", genome_end - c + 1);
	  //array_list_insert(node_end->position - lim_ref, final_positions);
	  break;
	} else {
	  array_list_insert(node_end->position + 1, final_positions);
	  //printf(" Insert 3) %lu\n", node_end->position + 1);
	}
      }
    }
	    
    //3.Select the correct start
    int pos;
    array_list_t *starts_list = ((end_data_t *)node_end->data)->list_starts;

    array_list_clear(final_ends, NULL);
    //Filter ends before cal_next->start
    for (int e = 0; e < array_list_size(starts_list); e++) {
      size_t start = array_list_get(e, starts_list);
      if (filter_pos) {
	if (start > filter_pos) {
	  array_list_insert(start, final_ends);
	}
      } else {
	array_list_insert(start, final_ends);
      }
    }
    if (array_list_size(final_ends) >= 1) {
      //Select best start
      char s_reference[array_list_size(final_ends)][2048];
      int references_len[array_list_size(final_ends)];
      //Making reference...
      for (int s = 0; s < array_list_size(final_ends); s++) {
	genome_end = array_list_get(s, final_ends) - 1;
	genome_start = genome_end - gap_close - 1;
	genome_read_sequence_by_chr_index(s_reference[s], 0, 
					  cal_chromosome_id - 1, 
					  &genome_start, &genome_end, genome); 	      
	references_len[s] = strlen(s_reference[s]) - 1;
	//printf("(%lu-%lu)Reference END_SP [START_SP]----[[[END_SP]]](%i): %s\n", genome_start, genome_end, s, s_reference[s]);
      } 
      //Select the best start
      //int lim = (array_list_get(1, final_starts)) - (array_list_get(0, final_starts));
      int anchor_nt = 20; 
      if (anchor_nt > gap_close) { 
	anchor_nt = gap_close;
      }
      pos = 0;
      int *ref_matches = (int *)calloc(array_list_size(final_ends), sizeof(int));
      int *ref_mismatches = (int *)calloc(array_list_size(final_ends), sizeof(int));
      //printf("anchor nt = %i, array_list_size(final_ends) = %i\n", anchor_nt, array_list_size(final_ends));
      for (int c = 0; c < anchor_nt; c++) {
	for (int s = 0; s < array_list_size(final_ends); s++) {	  
	  if (read_pos - c < 0 || references_len[s] - c < 0) { exit(-1); }
	  //printf("\t[%c vs %c]\n", query_map[read_pos - c], s_reference[s][references_len[s] - c]);
	  if (query_map[read_pos - c] != s_reference[s][references_len[s] - c]) {
	    ref_mismatches[s]++;
	  } else {
	    ref_matches[s]++;
	  }
	}
      }
      float score, max_score = 0.0;
      for (int s = 0; s < array_list_size(final_ends); s++) {
	score = ref_matches[s]*0.5 - ref_mismatches[s]*0.4;
	if (score > max_score) { 
	  max_score = score;
	  pos = s;
	}
      }		
      max_dist = 3;
      if (ref_mismatches[pos] > max_dist) {
	break;
      }

      final_dist += ref_mismatches[pos]; 
      final_pos = array_list_get(pos, final_ends);
      array_list_insert(final_pos, final_positions);
      //printf(" Insert 4) %lu\n", final_pos);
      genome_end = final_pos - anchor_nt - 1;
      gap_close -= anchor_nt;
      read_pos -= anchor_nt;
      //printf("FINAL MATCHES (%i) VS FINAL MISMATCHES (%i)\n", ref_matches[s], ref_mismatches[s]);
      //printf("GAP CLOSE (%i) READ POS (%i)\n", gap_close, read_pos);

      free(ref_mismatches);
      free(ref_matches);

      if (anchor_nt < 20) {
	array_list_insert(final_pos - anchor_nt, final_positions);
	//printf(" Insert 5) %lu\n", final_pos - anchor_nt);
	map = 1;
	break;
      } else {
	if (metaexon_search(cal_strand, cal_chromosome_id - 1, 
			    final_pos, final_pos - anchor_nt, 
			    &final_metaexon, metaexons)) {
	  if (final_metaexon) {
	    if (final_metaexon->left_closed) {
	      array_list_clear(ends_targets, NULL);
	      //printf("FOUND!!\n");
	      for (int s_0 = 0; s_0 < array_list_size(final_metaexon->left_breaks); s_0++) {
		node_end = array_list_get(s_0, final_metaexon->left_breaks);
		if (final_pos + 10 >= node_end->position) {
		  array_list_insert(node_end, ends_targets);		    
		}
	      }
	    } else {
	      //We can close the gap
	      int exon_length = final_metaexon->end - final_metaexon->start;
	      dist = 0;
	      if (exon_length >= gap_close) {
		genome_start = genome_end - gap_close - 1;
		genome_read_sequence_by_chr_index(reference, 0, 
						  cal_chromosome_id - 1, 
						  &genome_start, &genome_end, genome);		
		//printf("CLOSE GAP (%i): %s\n", gap_close, reference);
		reference_len = strlen(reference) - 1;
		for (int c  = 0; c < gap_close; c++) { 
		  //printf("\t[%c vs %c]\n", query_map[read_pos - c], reference[reference_len - c]);
		  if (query_map[read_pos - c] != reference[reference_len - c]) { 
		    dist++;
		  }
		}
		if (dist < 5) {
		  array_list_insert(final_pos - gap_close - anchor_nt, final_positions);
		  //printf(" Insert 6) %lu\n", final_pos - gap_close - anchor_nt);
		  final_dist += dist; 
		  map = 1;
		}
	      }
	      break;
	    }
	  }
	} else {
	  //printf("NOT  FOUND\n");
	  break;
	}
      }
    } else {
      break;
    }
  }

  cigar_code_t *cigar_code = NULL;

  if (map) {
    size_t pos_prev = cal->end, pos_next;
    cigar_op_t *op;
    cigar_code = cigar_code_new();
    cigar_code->distance = final_dist;

    pos_prev = cal->start - 1;
    for (int sp = 0; sp < array_list_size(final_positions); sp++) {
      pos_next = array_list_get(sp, final_positions);
      //printf("%lu - %lu(%i/%i)\n", pos_prev, pos_next, sp, array_list_size(final_positions));
      if (sp % 2 == 0) {		
	op = cigar_op_new(pos_prev - pos_next + 1, 'M');
      } else {
	int aux = pos_prev - pos_next + 1;
	op = cigar_op_new(aux, 'N');
      }
      //printf("\tADD %i%c\n", op->number, op->name);
      cigar_code_insert_first_op(op, cigar_code);
      pos_prev = pos_next - 1;
    }

    //printf(" :::::::::::::::::::::::::::: %s\n", new_cigar_code_string(cigar_code));
    //cal->info = cigar_code;
    //s_prev->info = cigar_code;
    //meta_alignment = meta_alignment_new();
    //meta_alignment_insert_cal(cal, meta_alignment);
    //meta_alignment_set_status(META_CLOSE, meta_alignment);
	      
  }

  array_list_free(final_positions, NULL);
  array_list_free(ends_targets, NULL);
  array_list_free(final_ends, NULL);

  return cigar_code;

}

cigar_code_t *search_double_anchors_cal(char *query_map,
					cal_t *first_cal, cal_t *last_cal,
					metaexons_t *metaexons, genome_t *genome,
					fastq_read_t *fq_read, int *type) {
  int l_found = 0;
  int r_found = 0;
  int read_nt, read_orig_nt;
  metaexon_t *first_metaexon, *last_metaexon;
  seed_region_t *s_prev = linked_list_get_last(first_cal->sr_list);    
  seed_region_t *s_next = linked_list_get_first(last_cal->sr_list);	
  //int found_sp = 0;
  size_t genome_start, genome_end;
  char reference[2048];
  //meta_alignment_t *meta_alignment = NULL;
  cigar_code_t *cigar_code = NULL;

  *type = META_ALIGNMENT_MIDDLE;

  if (metaexon_search(first_cal->strand, first_cal->chromosome_id - 1,
		      first_cal->start, first_cal->end, &first_metaexon,
		      metaexons)) {
    //printf("First Meta %i\n", first_metaexon->right_closed);
    l_found = first_metaexon->right_closed;
  }	 

  if (metaexon_search(last_cal->strand, last_cal->chromosome_id - 1,
		      last_cal->start, last_cal->end, &last_metaexon,
		      metaexons)) {
    //printf("Last Meta %i\n", last_metaexon->left_closed);
    r_found = last_metaexon->left_closed;
  }
	  
  //printf("l_found = %i, r_found = %i\n", l_found, r_found);
  if (l_found && r_found) {
    //Is the correct splice?? Ask to length exons...
    //printf("Found Metaexon.... Is the correct???\n");	    
    array_list_t *introns_targets = array_list_new(10, 1.5f, COLLECTION_MODE_ASYNCHRONIZED);
    array_list_t *starts_targets  = array_list_new(10, 1.5f, COLLECTION_MODE_ASYNCHRONIZED);
    array_list_t *ends_targets  = array_list_new(10, 1.5f, COLLECTION_MODE_ASYNCHRONIZED);
    intron_t *intron;
    int size_ex_l, size_ex_r, dsp_l, dsp_r;
    avl_node_t *node_start, *node_end;

    //1st Search posible starts splice junctions
    for (int s_0 = 0; s_0 < array_list_size(first_metaexon->right_breaks); s_0++) {
      node_start = array_list_get(s_0, first_metaexon->right_breaks);
      //printf("\t\t Node start %lu\n", node_start->position);
      if (first_cal->end - 10 <= node_start->position) {
	array_list_insert(node_start, starts_targets);
      }
    }
    
    for (int e_0 = 0; e_0 < array_list_size(last_metaexon->left_breaks); e_0++) {
      node_end = array_list_get(e_0, last_metaexon->left_breaks);
      //printf("\t\t Node end %lu\n", node_end->position);
      if (last_cal->start + 10 >= node_end->position) {
	array_list_insert(node_end, ends_targets);
      }
    } 

    for (int s_0 = 0; s_0 < array_list_size(starts_targets); s_0++) {
      node_start = array_list_get(s_0, starts_targets);
      array_list_t *ends_list = ((start_data_t *)node_start->data)->list_ends;
      for (int e_0 = 0; e_0 < array_list_size(ends_list); e_0++) {
	splice_end_t *splice_end = array_list_get(e_0, ends_list);
	for (int e_1 = 0; e_1 < array_list_size(ends_targets); e_1++) {
	  avl_node_t *node_aux = array_list_get(e_1, ends_targets);
	  //printf("\t\t Node end %lu == %lu\n", splice_end->end, node_aux->position);
	  if (splice_end->end == node_aux->position) {
	    //*******************************************//
	    //     [  L_EXON.1  ]-------[  R_EXON.1  ]   //
	    //     [  L_EXON.2 ]----------[R_EXON.2  ]   //
	    //     [  L_EXON.3]-------------[R_EXON.3]   //
	    //*******************************************//
	    intron = intron_new(last_cal->strand, last_cal->chromosome_id, 
				node_start->position, splice_end->end);
	    array_list_insert(intron, introns_targets);
	  }
	}
      }
    }
	  
    read_nt = s_next->read_start - s_prev->read_end - 1; //ok
    if (read_nt < 0) { read_nt = 0; }	    

    read_orig_nt = read_nt;
    //We not found introns... We have one exon between CALs??  //min-exon 30nt
    if (array_list_size(introns_targets) <= 0) { 
      //printf("INTRON TARGETS, SEARCH LEFT\n");
      read_nt = fq_read->length - (s_prev->read_end + 1);
      cigar_code = search_left_single_anchor(read_nt, 
					     first_cal,
					     0,
					     first_metaexon->right_breaks,
					     query_map, metaexons, genome);
      *type = META_ALIGNMENT_LEFT;
    }
    
    //metaexons_show(metaexons);
    //assert(array_list_size(introns_targets) != 0);
    //printf("read_gap = (%i - %i) %i, introns_targets = %i\n", 
    //	   s_next->read_start, s_prev->read_end, read_nt, array_list_size(introns_targets));

    size_t s_intron, e_intron;
    int distance = 0;
    size_t s_prev_read_end, s_prev_genome_end, first_cal_end;
    size_t s_next_read_start, s_next_genome_start, last_cal_start;

    //Check introns... good luck! :)
    for (int in = 0; in < array_list_size(introns_targets); in++) {
      intron = array_list_get(in, introns_targets);
      //printf("Intron target %i:[%lu-%lu] vs cal[%lu-%lu]\n", in, intron->start, intron->end,
      //     first_cal->end, last_cal->start);
      
      read_nt             = read_orig_nt;
      s_prev_read_end     = s_prev->read_end;
      s_prev_genome_end   = s_prev->genome_end;
      first_cal_end       = first_cal->end;
      s_next_read_start   = s_next->read_start; 
      s_next_genome_start = s_next->genome_start; 
      last_cal_start      = last_cal->start; 

      dsp_l = 0;
      dsp_r = 0;

      size_ex_l = (int)(intron->start - (first_cal->end + 1));
      //printf("\t LEFT = (i)%lu - %lu = %i\n", intron->start, first_cal->end + 1, size_ex_l);
      if (size_ex_l < 0) {
	dsp_l = abs(size_ex_l);
	s_prev_read_end -= dsp_l;
	s_prev_genome_end -= dsp_l;
	first_cal_end -= dsp_l;
	read_nt += dsp_l;
	size_ex_l = 0;
      }
	    
      size_ex_r = (int)(last_cal->start - (intron->end + 1));
      //printf("\t RIGHT = %lu - (i)%lu = %i\n", last_cal->start, intron->end + 1, size_ex_r);
      if (size_ex_r < 0) { 
	dsp_r = abs(size_ex_r); 
	s_next_read_start += dsp_r;
	s_next_genome_start += dsp_r;
	last_cal_start += dsp_r;
	read_nt += dsp_r;
	size_ex_r = 0;
      }
      
      assert(read_nt >= 0);
      
      //printf("size_ex_l = %i, size_ex_r = %i, read_nt = %i\n", size_ex_l, size_ex_r, read_nt);
      if (read_nt == (size_ex_l + size_ex_r)) {
	//Close gap		
	if (size_ex_l > 1) {
	  genome_start = first_cal_end + 1;
	  genome_end   = first_cal_end + size_ex_l;
	  genome_read_sequence_by_chr_index(reference, 0, 
					    first_cal->chromosome_id - 1, 
					    &genome_start, &genome_end, genome);
	  //printf("Ref-l-%s\n", reference);
	  for (int c = 0; c < size_ex_l; c++) {
	    //printf("%c/%c | ", query_map[s_prev->read_end + 1 + c], reference[c]);
	    if (query_map[s_prev_read_end + 1 + c] != reference[c]) {
	      distance++;
	    }
	  }
	  //printf("\n");
	}else if (size_ex_l == 1) {
	  distance++;
	}
	//printf("1st-Distance-l %i\n", distance);
	
	s_prev_read_end += size_ex_l;
	s_prev_genome_end += size_ex_l;
	first_cal_end += size_ex_l;
	
	if (size_ex_r > 1) {
	  genome_start = last_cal_start - size_ex_r;
	  genome_end   = last_cal_start;
	  genome_read_sequence_by_chr_index(reference, 0, 
					    first_cal->chromosome_id - 1, 
					    &genome_start, &genome_end, genome);		  
	  //printf("Ref-r-%s\n", reference);
	  for (int c = 0; c < size_ex_r; c++) {
	    //printf("%c/%c | ", query_map[s_next->read_start - size_ex_r + c], reference[c]);
	    if (query_map[s_next_read_start - size_ex_r + c] != reference[c]) {
	      distance++;
	    }
	  }
	  //printf("\n");
	}else if (size_ex_r == 1) {
	  distance++;
	}

	s_next_read_start -= size_ex_r;
	s_next_genome_start -= size_ex_r;
	last_cal_start -= size_ex_r;
	
	//seed_region = seed_region_new(s_next->read_start, s_next->read_end, 
	//			      s_next->genome_start, s_next->genome_end, 1);
	//linked_list_insert_last(seed_region, first_cal->sr_list);		
	//first_cal->end = last_cal->end;
	//printf("2nd-Distance-r %i\n", distance);
	if (distance > (size_ex_l + size_ex_r)/2) {
	  continue;
	}

	cigar_code = cigar_code_new();

	size_ex_l -= dsp_l;
	size_ex_r -= dsp_r;

	cigar_code_append_new_op(size_ex_l, 'M', cigar_code);
	cigar_code_append_new_op(last_cal_start - first_cal_end - 1, 'N', cigar_code);
	cigar_code_append_new_op(size_ex_r, 'M', cigar_code);		
	cigar_code->distance = distance;
	
	//seed_region_t *s = linked_list_get_first(first_cal->sr_list);
	//s->info = cigar_code;
	//first_cal->info = cigar_code;
	//cal_free(last_cal);	
	//meta_alignment = meta_alignment_supra_new(META_CLOSE, first_cal);
	//array_list_insert(meta_alignment, meta_alignments_list);
	//found_sp = 1;
	break; //If found break loop	
      } else {
	LOG_DEBUG_F("@@@PROBLEM CALCULATING SPLICE JUNCTION id = %s\n", fq_read->id);
      }		
      //printf("Distance to close start splice %i\n", size_ex_l);
      //printf("Distance to close end   splice %i\n", size_ex_r);	      
    } //End loop introns	    
    array_list_free(starts_targets, NULL);
    array_list_free(ends_targets, NULL);
    array_list_free(introns_targets, intron_free);
  }//End if (l_found && r_found)
  
  return cigar_code;

}

//IF MODE == 0, DOUBLE ANCHORS
//IF MODE == 1, LEFT ANCHOR
//IF MODE == 2, RIGHT ANCHOR
int generate_cals_between_anchors (int mode,
				   int num_chromosomes,
				   cal_t *first_cal, 
				   cal_t *last_cal, 
				   fastq_read_t *fq_read, 
				   array_list_t *seeds_list, 
				   array_list_t *cals_list, 
				   cal_optarg_t *cal_optarg) {
  //printf("GENERATE CALS\n");
  seed_region_t *s_prev, *s_next;
  int min_seeds = 0, max_seeds = 1000;
  region_t *region_prev = NULL, *region_next = NULL;
  array_list_t *mapping_list = array_list_new(500, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  int strand, chromosome_id;

  array_list_clear(cals_list, NULL);

  if (mode == CALING_DOUBLE_ANCHORS || 
      mode == CALING_LEFT_ANCHORS) {
    s_prev = linked_list_get_last(first_cal->sr_list);    
    region_prev = region_bwt_new(first_cal->chromosome_id,
				 first_cal->strand,
				 first_cal->start, 
				 first_cal->end,
				 s_prev->read_start,
				 s_prev->read_end,
				 fq_read->length,
				 0);
    array_list_insert(region_prev, mapping_list);
    strand = first_cal->strand;
    chromosome_id = first_cal->chromosome_id;
  }

  //printf("Region FORW(%i)[%i:%lu|%i-%i|%lu]\n ", region_prev->id, region_prev->chromosome_id, 
  //	   region_prev->start, region_prev->seq_start, region_prev->seq_end, region_prev->end);	  
  if (mode == CALING_DOUBLE_ANCHORS || 
      mode == CALING_RIGHT_ANCHORS) {    
    s_next = linked_list_get_first(last_cal->sr_list);
    int id, gap;
    if (mode == CALING_DOUBLE_ANCHORS) {
      gap = s_next->read_start - s_prev->read_end;
      id = (gap / 16) + ((gap % 16) > 0) + 1;
    } else {
      id = 1;
    }
    region_next = region_bwt_new(last_cal->chromosome_id,
				 last_cal->strand,
				 last_cal->start, 
				 last_cal->end,
				 s_next->read_start,
				 s_next->read_end,
				 fq_read->length,
				 id);
    array_list_insert(region_next, mapping_list);	      
    strand = last_cal->strand;
    chromosome_id = last_cal->chromosome_id;
  }
  //printf("Region BACK(%i)[%i:%lu|%i-%i|%lu]\n ", region_next->id, region_next->chromosome_id, 
  //	   region_next->start, region_next->seq_start, region_next->seq_end, region_next->end);
  
  int min_cal_size = cal_optarg->min_cal_size;
  if (mode == CALING_DOUBLE_ANCHORS ) {
    if (last_cal->start - first_cal->end < 50000) {
      cal_optarg->min_cal_size = 15;
    }
  }

  //printf(":::::::::::::: seeds %i\n", array_list_size(seeds_list));
  for (int j = 0; j < array_list_size(seeds_list); j++) {
    region_t *region = array_list_get(j, seeds_list);
    if (region == NULL) { exit(-1); }

    //printf("@@@@@@@(%i)Region [%i:%lu|%i-%i|%lu]: ", region->id, region->chromosome_id, 
    //	     region->start, region->seq_start, region->seq_end, region->end);
    if (region->strand == strand && 
	region->chromosome_id == chromosome_id) {
      if (mode == CALING_DOUBLE_ANCHORS) {
	if (region_prev->end - 16 <= region->start &&
	    region_next->start + 16 >= region->end) {
	  //printf("Insert\n");
	    array_list_insert(region, mapping_list);
	} else { /*printf("Not Insert\n");*/ }
      } else if (mode == CALING_LEFT_ANCHORS) {
	if (region_prev->end - 16 <= region->start) {
	  //printf("Insert\n");
	  array_list_insert(region, mapping_list);
	} else { /*printf("Not Insert\n");*/ } 
      } else {
	if (region_next->start + 16 >= region->end) {
	  //printf("Insert\n");
	  array_list_insert(region, mapping_list);
	} else { /*printf("Not Insert\n");*/ }
      }
    } else {
      //printf("Not Insert\n");
    }
  }

  bwt_generate_cal_list_linked_list(mapping_list,
				    cal_optarg,
				    &min_seeds, &max_seeds,
				    num_chromosomes,
				    cals_list,
				    fq_read->length);

  cal_optarg->min_cal_size = min_cal_size;
  int num_cals = array_list_size(cals_list);
  int founds[num_cals];
  int found = 0;
  cal_t *cal;
  for (size_t j = 0; j < num_cals; j++) {
    cal = array_list_get(j, cals_list);
    LOG_DEBUG_F("\tcal %i of %i: sr_list size = %i (cal->num_seeds = %i) %i:%lu-%lu\n", 
		j, num_cals, cal->sr_list->size, cal->num_seeds,
		cal->chromosome_id, cal->start, cal->end);
    founds[j] = 0;
    if (cal->sr_list->size > 0) {
      int start = 0;
      for (linked_list_item_t *list_item = cal->sr_list->first; list_item != NULL; list_item = list_item->next) {
	seed_region_t *s = list_item->item;		  
	LOG_DEBUG_F("\t\t:: star %lu > %lu s->read_start\n", start, s->read_start);
	if (start > s->read_start) {
	  LOG_DEBUG("\t\t\t:: remove\n");
	  found++;
	  founds[j] = 1;
	}
	start = s->read_end + 1;
      }
    } else {
      found++;
      founds[j] = 1;
    }
  }

  //Check CALs TODO:Delete?Duplicate??
  int start = 0;
  for (int c = 0; c < num_cals; c++) { 
    cal = array_list_get(c, cals_list);
    seed_region_t *s_prev_aux = linked_list_get_first(cal->sr_list);    
    seed_region_t *s_next_aux = linked_list_get_last(cal->sr_list);
    if (s_prev_aux == NULL || s_next_aux == NULL) { 
      founds[c] = 1;
      found++;
    } else {
      //printf("CAL %i: %i < %i\n", c, s_prev_aux->read_start, start);
      if (s_prev_aux->read_start < start) {
	founds[c] = 1;
	found++;	
      } else {
	start = s_next_aux->read_end;
      }
    }
  }
	      
  if (found) {
    for (int j = num_cals - 1; j >= 0; j--) {
      if (founds[j]) {
	cal = array_list_remove_at(j, cals_list);
	array_list_free(cal->candidates_seeds_start, NULL);
	array_list_free(cal->candidates_seeds_end, NULL);
	if (cal->sr_list != NULL) { 
	  linked_list_free(cal->sr_list, NULL); 
	}
	if (cal->sr_duplicate_list != NULL) { 
	  linked_list_free(cal->sr_duplicate_list, NULL); 
	}
	free(cal);
	//array_list_set(j, NULL, cals_list);
	//cal_free(cal);
      }
    }
  }
  //seeds 
  
  cal_free(first_cal);
  cal_free(last_cal);

  //printf(" ==== FINAL CALS (%i) ==== \n", array_list_size(cals_list));
  //for (size_t j = 0; j < array_list_size(cals_list); j++) {
  //cal = array_list_get(j, cals_list);
  //cal_print(cal);
  //}
  //printf(" ==== END FINAL CALS ==== \n");
  array_list_free(mapping_list, NULL);

  if (region_prev) { region_bwt_free(region_prev); }
  if (region_next) { region_bwt_free(region_next); }

  return array_list_size(cals_list);
}

void *meta_alignment_regenerate_cigar(meta_alignment_t *meta_alignment, char *query_map, genome_t *genome) {
  array_list_t *cals_list = meta_alignment->cals_list;

  /*
  printf("==== CALS LIST %i ====\n", array_list_size(cals_list));
  for (int j = 0; j < array_list_size(cals_list); j++) {
    cal_t *cal = array_list_get(j, cals_list);
    cal_print(cal);
  }
  printf("==== CALS LIST ====\n");      
  */

  int num_cals = array_list_size(cals_list);
  if (!num_cals) { return; }

  cigar_code_t *cigar_code_n = cigar_code_new();
  cal_t *cal_prev = array_list_get(0, cals_list);
  cigar_code_t *cigar_code_prev = cal_prev->info;

  if (meta_alignment->cigar_left != NULL) {
    cigar_code_t *cigar_code_aux = meta_alignment->cigar_left;
    for (int i = 0; i < array_list_size(cigar_code_aux->ops); i++) {
      cigar_op_t *op = array_list_get(i, cigar_code_aux->ops);
      cigar_code_append_new_op(op->number, op->name, cigar_code_n);
    } 
    cigar_code_n->distance += cigar_code_aux->distance;
  }

  for (int i = 0; i < array_list_size(cigar_code_prev->ops); i++) {
    cigar_op_t *op = array_list_get(i, cigar_code_prev->ops);
    cigar_code_append_new_op(op->number, op->name, cigar_code_n);
  }

  cal_t *cal_next;
  int read_gap;
  seed_region_t *seed_prev, *seed_next;
  int read_start, read_end;
  size_t genome_start, genome_end;
  char reference_prev[2048];
  char reference_next[2048];
  int max_distance = 3;
  int distance;
  int prev_M, next_M;

  if (num_cals > 1) {
    for (int i = 1; i < num_cals; i++) {
      cal_next = array_list_get(i, cals_list);
      seed_prev = linked_list_get_last(cal_prev->sr_list);
      seed_next = linked_list_get_first(cal_next->sr_list);

      read_gap = seed_next->read_start - seed_prev->read_end - 1;

      read_start = seed_prev->read_end + 1;
      read_end = seed_next->read_start - 1;

      genome_start = seed_prev->genome_end + 1;
      genome_end = seed_prev->genome_end + read_gap + 5;      
      genome_read_sequence_by_chr_index(reference_prev, 0, 
					cal_next->chromosome_id - 1,
					&genome_start, &genome_end,
					genome);

      genome_start = seed_next->genome_start - read_gap - 5;
      genome_end = seed_next->genome_start - 1;
      genome_read_sequence_by_chr_index(reference_next, 0, 
					cal_next->chromosome_id - 1,
					&genome_start, &genome_end,
					genome);
      
      cigar_code_t *cigar_code_next = cal_next->info;

      //printf("REF-PREV: %s\n", reference_prev);
      //printf("REF-NEXT: %s\n", reference_next);
      //printf("read_gap = %i\n", read_gap);

      int c = 0;
      int r_pos = strlen(reference_next) - 1;

      distance = 0;
      while (distance < max_distance && 
	     c < read_gap) { 
	if (query_map[read_start + c] != reference_prev[c]) {	  
	  distance++;
	}
	c++;
      }
      prev_M = c;
      if (prev_M > 0) {
	cigar_code_append_new_op(prev_M, 'M', cigar_code_n);
      }
      //printf("Distance %i and Matches PREV %i\n", distance, prev_M);

      distance = 0;
      while (distance < max_distance && 
	     c < read_gap) {
	if (query_map[read_end - c] != reference_next[r_pos - c]) {	  
	  distance++;
	}
	c++;
      }
      next_M = c - prev_M;
      
      int num_I = read_end - (read_start + prev_M + next_M) + 1;

      genome_start = seed_prev->genome_end + prev_M;
      genome_end = seed_next->genome_start - next_M;

      if (num_I > 0) {
	cigar_code_append_new_op(num_I, 'I', cigar_code_n);
      }

      cigar_code_append_new_op(genome_end - genome_start - 1 - num_I, 'N', cigar_code_n);
            
      if (next_M > 0) {
	cigar_code_append_new_op(next_M, 'M', cigar_code_n);
      }
      //printf("Distance %i and Matches NEXT %i\n", distance, next_M);

      for (int i = 0; i < array_list_size(cigar_code_next->ops); i++) {
	cigar_op_t *op = array_list_get(i, cigar_code_next->ops);
	//printf("NEXT:: %i%c\n", op->number, op->name);
	cigar_code_append_new_op(op->number, op->name, cigar_code_n);
      }
            
      cal_prev = cal_next;
    }

  } else {
    exit(-1);
  }

  if (meta_alignment->cigar_right != NULL) {
    cigar_code_t *cigar_code_aux = meta_alignment->cigar_right;
    for (int i = 0; i < array_list_size(cigar_code_aux->ops); i++) {
      cigar_op_t *op = array_list_get(i, cigar_code_aux->ops);
      //printf("R:: %i%c\n", op->number, op->name);
      cigar_code_append_new_op(op->number, op->name, cigar_code_n);
    } 
    cigar_code_n->distance += cigar_code_aux->distance;
  } 
  
  meta_alignment->cigar_code = cigar_code_n;

  //printf("--->: %s\n", new_cigar_code_string(meta_alignment->cigar_code));

}

//NEW FUNCTION!
int apply_sw_rna(sw_server_input_t* input_p, batch_t *batch) {
  LOG_DEBUG("========= SPLICE JUNCTION SEARCH =========\n");
  size_t max_intron_size = 500000;  
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  sw_optarg_t *sw_optarg = &input_p->sw_optarg;
  genome_t *genome = input_p->genome_p;
  size_t num_reads = array_list_size(mapping_batch->fq_batch);
  size_t num_targets = mapping_batch->num_targets;
  metaexons_t *metaexons = input_p->metaexons;
  cal_optarg_t *cal_optarg = input_p->cal_optarg_p;
  avls_list_t *avls_list = input_p->avls_list;
  bwt_optarg_t *bwt_optarg = input_p->bwt_optarg_p;
  bwt_index_t *bwt_index = input_p->bwt_index_p;
  linked_list_t *buffer = input_p->buffer;
  linked_list_t *buffer_hc = input_p->buffer_hc;

  array_list_t *cals_list, *fusion_cals, *fusion_cals_aux;
  cal_t *cal, *cal_prev, *cal_next, *first_cal, *last_cal;
  fastq_read_t *fq_read;
  linked_list_iterator_t itr;  
  seed_region_t *s, *s_prev, *s_next;
  cigar_code_t *cigar_code, *cigar_code_prev, *cigar_code_aux;
  cigar_code_t *alig_cigar_code;
  cigar_op_t *cigar_op_start, *cigar_op_end, *cigar_op, *cigar_op_prev, *cigar_op_aux;
  int *new_targets = (int *)calloc(mapping_batch->num_allocated_targets, sizeof(int));
  array_list_t *merge_cals;
  linked_list_t *linked_list;
  seed_region_t *seed_region;

  //float *cals_score = (float *)calloc(100, sizeof(float));
  float score;
  char reference[2048];
  char reference_prev[2048];
  char reference_next[2048];
  char reference_aux[2048];
  char query[2048];
  char query_revcomp[2048];
  alignment_t *alignment;
  char q[2048];
  char r[2048];

  char **rev_comp = (char **)calloc(num_reads, sizeof(char *));

  fusion_coords_t *extrem_coords[2*40*mapping_batch->num_allocated_targets];
  fusion_coords_t *sp_coords[2*40*mapping_batch->num_allocated_targets];
  cigar_code_t *extrem_cigars[2*40*mapping_batch->num_allocated_targets];
  cigar_code_t *sp_cigars[2*40*mapping_batch->num_allocated_targets];

  char *sequence;
  char *query_ref;
  char *quality_map, *query_map;
  //float scores_ranking[mapping_batch->num_allocated_targets][50];
  float *cals_score;
  char cigar_str[1024];

  cigar_op_t *first_op;
  char *match_seq, *match_qual, *optional_fields, *p;
  int match_start, match_len, optional_fields_length, AS;
  float norm_score;

  int delete_not_cigar;
  int padding_left, padding_right, len_query;  
  int num_sw = 0;
  int num_sw_sp = 0;
  int num_sw_ext = 0;
  size_t num_new_targets = 0;
  int coverage;
  int read_length;
  int cal_pos = 0;
  int number_of_best;
  //int num_extrem_ok = 0;
  //int num_sp_ok = 0;

  //Change to report most alignments
  int seed_err_size = 20;
  int n_alignments = 1;
  int target;
  int c, exact_nt;
  size_t genome_start, genome_end, genome_start2, genome_end2, read_start, read_end;
  int flank = 20;
  int n_delete;
  int sw_pos;
  int seeds_nt; 
  int e1;
  int e2;
  int ref_pos;
  int nt_2;
  
  int start_seeding, end_seeding;
  int lim_start, lim_end;

  int min_intron_size = 40;

  metaexon_t *first_metaexon, *last_metaexon;
  int l_found, r_found; 
  size_t num_cals;
  register size_t t;
  register int i, j;

  int meta_type;
  meta_alignment_t *meta_alignment;
  sw_item_t *sw_item;  
  sw_depth_t sw_depth;
  sw_depth.depth = 0;
  sw_multi_output_t *output = sw_multi_output_new(MAX_DEPTH);
  //array_list_t **sw_items_list = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
  //array_list_t **alignments_list = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
  array_list_t **meta_alignments_list = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
  array_list_t *seeds_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *final_positions = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED); 
  int make_seeds = 0;

  array_list_t *cals_targets = array_list_new(50, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  int *post_process_reads = (int *)calloc(num_reads, sizeof(int));
  float scores_ranking[num_reads][50];
  int read_nt;

  /*for (int i = 0; i < num_reads; i++) {
    //sw_items_list[i]   = array_list_new(20, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    //alignments_list[i] = array_list_new(20, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    meta_alignments_list[i] = array_list_new(20, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    }*/

  /*
    number_of_best = merge_and_filter_cals(cals_targets[i], 
    cals_list,
    fq_read, num_cals, bwt_optarg,
    bwt_index, genome, 
    scores_ranking[i], NULL, NULL, &num_sw, 
    extrem_coords, target);
  */
  
  // array_list flag: 0 -> Not  BWT Anchors found (NOT_ANCHORS)        *
  //                  1 -> One  BWT Anchors found (SINGLE_ANCHORS)     *
  //                  2 -> Pair BWT Anchors found (DOUBLE_ANCHORS)     *
  //                  3 -> Alignments found       (ALIGNMENTS_FOUND)
  //                  4 -> Alignments exceded     (ALIGNMENTS_EXCEEDED)
  int flag;

  for (i = 0; i < num_reads; i++) {
    cals_list = mapping_batch->mapping_lists[i];
    flag = array_list_get_flag(cals_list);
    fq_read = array_list_get(i, mapping_batch->fq_batch);

    num_cals = array_list_size(cals_list);
    meta_alignments_list[i] = array_list_new(num_cals,
				     1.25f,
				     COLLECTION_MODE_ASYNCHRONIZED);

    //printf("WK_1ph-Process: == (%i CALs)Read %s ==\n", array_list_size(cals_list), fq_read->id);
    if (flag == DOUBLE_ANCHORS) {
      //printf("\tWK_1ph: -- DOUBLE ANCHOR PROOCESS --\n");
      //printf("<<<<@ %s\n", fq_read->id);
      //1st- Merge Double anchors CALs. TODO: make only merge, not generate score	
      //array_list_clear(mapping_batch->mapping_lists[i], cal_free);continue;
      make_seeds = 0;
      int found[array_list_size(cals_list)];
      int founds = 0;
      for (int p = 0; p < array_list_size(cals_list); p += 2) {	
	first_cal = array_list_get(p, cals_list);
	last_cal  = array_list_get(p + 1, cals_list);
	if (first_cal->end + fq_read->length >= last_cal->start) {
	  found[p] = 1;
	  found[p + 1] = 1;
	  founds++;
	} else {
	  found[p] = 0;
	  found[p + 1] = 0;	  
	}
      }

      if (founds) {
	//array_list_t *cals_aux = array_list_new(array_list_size(cals_list), 1.5f, COLLECTION_MODE_ASYNCHRONIZED);
	//printf("REPORT FILL GAPS\n");
	for (int p = 0; p < array_list_size(cals_list); p += 2) {
	  first_cal = array_list_get(p, cals_list);
	  if (first_cal->strand == 1) {
	    if (!rev_comp[i]) {
	      rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
	      strcpy(rev_comp[i], fq_read->sequence);
	      seq_reverse_complementary(rev_comp[i], fq_read->length);
	    }
	    query_map = rev_comp[i];
	  } else {
	    query_map = fq_read->sequence;
	  }
	  last_cal = array_list_get(p + 1, cals_list);
	  if (found[p] == 1) {
	    s = linked_list_get_first(last_cal->sr_list);
	    seed_region = seed_region_new(s->read_start, s->read_end, 
					  s->genome_start, s->genome_end, 1);
	    linked_list_insert_last(seed_region, first_cal->sr_list);
	    first_cal->end = last_cal->end;
	    cal_free(last_cal);
	    meta_alignment = meta_alignment_cal_new(first_cal);
	    meta_alignment_fill_gaps(META_ALIGNMENT_NONE,
				     meta_alignment, query_map, genome,
				     sw_optarg, output, metaexons, 
				     &sw_depth, avls_list);	    
	    array_list_insert(meta_alignment, meta_alignments_list[i]);
	  } else {
	    cal_free(first_cal);
	    cal_free(last_cal);
	  }
	}
	continue;
      } else {
	first_cal = array_list_get(0, cals_list);
	last_cal  = array_list_get(1, cals_list);
	s_prev = linked_list_get_last(first_cal->sr_list);    
	s_next = linked_list_get_first(last_cal->sr_list);	
	int gap = s_next->read_start - s_prev->read_end;
	if (gap > 20) {
	  make_seeds = 1;
	}
      }

      int found_sp;
      int first_cal_len, last_cal_len;
      for (int p = 0; p < array_list_size(cals_list); p += 2) {
	first_cal = array_list_get(p, cals_list);
	first_cal_len = first_cal->end - first_cal->start;
	last_cal  = array_list_get(p + 1, cals_list);
	last_cal_len = last_cal->end - last_cal->start;
	if (last_cal->strand == 1) {
	  if (!rev_comp[i]) {
	    rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
	    strcpy(rev_comp[i], fq_read->sequence);
	    seq_reverse_complementary(rev_comp[i], fq_read->length);
	  }
	  query_map = rev_comp[i];
	} else {
	  query_map = fq_read->sequence;
	}
	
	//printf(" \t(%i) [%i:%lu-%lu] - [%i:%lu-%lu]\n", first_cal->strand, first_cal->chromosome_id, first_cal->start, first_cal->end, 
	//     last_cal->chromosome_id, last_cal->start, last_cal->end);
	//cal_print(first_cal);
	//cal_print(last_cal);
	s_prev = linked_list_get_last(first_cal->sr_list);    
	s_next = linked_list_get_first(last_cal->sr_list);

	
	cigar_code = search_double_anchors_cal(query_map,
					       first_cal, last_cal,
					       metaexons, genome,
					       fq_read, &meta_type);
	if (cigar_code != NULL) { 
	  meta_alignment = meta_alignment_new();
	  meta_alignment_insert_cal(first_cal, meta_alignment);
	  if (meta_type == META_ALIGNMENT_LEFT) {
	    meta_alignment_insert_cigar(cigar_code, CIGAR_ANCHOR_LEFT, 0, meta_alignment);
	  } else {
	    meta_alignment_insert_cal(last_cal, meta_alignment);
	    meta_alignment_insert_cigar(cigar_code, CIGAR_SIMPLE_MIDDLE, 0, meta_alignment);
	  }
	  meta_alignment_fill_gaps(meta_type,
				   meta_alignment, query_map, genome,
				   sw_optarg, output, metaexons, 
				   &sw_depth, avls_list);	  
	  array_list_insert(meta_alignment, meta_alignments_list[i]);
	} else {
	  //====(((((S M I T H - W A T E R M A N    P H A S E)))))====//	
	  //TODO: Extend before caling
	  cal_t *father_cal;
	  int gap = s_next->read_start - s_prev->read_end;
	  array_list_t *init_list = array_list_new(40, 1.25f, 
						   COLLECTION_MODE_ASYNCHRONIZED);
	  array_list_t *aux_list = array_list_new(40, 1.25f, 
						  COLLECTION_MODE_ASYNCHRONIZED);
	  region_t *region_prev, *region_next;
	  if (make_seeds) {
	    int read_start, read_end;
	    if (first_cal->strand == 1) {
	      read_start = fq_read->length - s_next->read_start;
	      read_end = (fq_read->length - (s_prev->read_end + 1)) - 1;
	    } else {
	      read_start = s_prev->read_end + 1;
	      read_end = s_next->read_start - 1;
	    }	  
	    //printf("SEEDING BETWEEN +(%i - %i), -(%i - %i)\n", s_prev->read_end + 1, s_next->read_start - 1, read_start, read_end);
	    bwt_map_exact_seeds_by_region(read_start, read_end,
					  fq_read->sequence, 16, 16,
					  bwt_optarg, bwt_index,
					  seeds_list);
	    make_seeds = 0;
	  }	  


	  //array_list_t *mapping_list = array_list_new(array_list_size(seeds_list) + 2, 1.25f, 
	  //COLLECTION_MODE_ASYNCHRONIZED);
	  generate_cals_between_anchors(CALING_DOUBLE_ANCHORS,
					genome->num_chromosomes,
					first_cal, 
					last_cal, 
					fq_read, 
					seeds_list, 
					aux_list, 
					cal_optarg);
	  
	  if (array_list_size(aux_list) > 0) {
	    first_cal = array_list_get(0, aux_list);
	    //first_cal->fill = 1;
	    if (array_list_size(aux_list) == 1) {
	      meta_alignment = meta_alignment_cal_new(first_cal);
	      //meta_alignment_set_status(META_CLOSE, meta_alignment);
	      meta_type = META_ALIGNMENT_NONE;
	    } else {
	      meta_alignment = meta_alignment_cals_new(aux_list);		
	      meta_type = META_ALIGNMENT_MIDDLE;
	      for (int c = 1; c < array_list_size(aux_list); c++) { 
		last_cal = array_list_get(c, aux_list);
		//last_cal->fill = 1;
		//Generate Fusion Reference...
		//printf("\t->%i/%i (%i)[%i:%lu-%lu] vs (%i)[%i:%lu-%lu], %p, %p, %p, %p\n", c, array_list_size(aux_list), 
		//first_cal->strand, first_cal->chromosome_id, 
		//first_cal->start, first_cal->end, 
		//last_cal->strand, last_cal->chromosome_id,
		//last_cal->start, last_cal->end, first_cal, last_cal, first_cal->sr_list, last_cal->sr_list); 

		avl_node_t *node_avl_start, *node_avl_end;
		seed_region_t *s_prev_aux = linked_list_get_last(first_cal->sr_list);    
		seed_region_t *s_next_aux = linked_list_get_first(last_cal->sr_list);
		int distance_aux;

		int nt = search_simple_splice_junction(s_prev_aux, s_next_aux,
						       first_cal->chromosome_id, 
						       first_cal->strand,
						       query_map, genome, 
						       avls_list, &node_avl_start,
						       &node_avl_end, &distance_aux);
		if (nt) {
		  cigar_code = cigar_code_new();
		  cigar_code->distance = distance_aux;
		  cigar_code_append_new_op(nt, 'N', cigar_code);

		  metaexon_insert(first_cal->strand, first_cal->chromosome_id - 1,
				  s_prev_aux->genome_start, node_avl_start->position, 40,
				  METAEXON_RIGHT_END, node_avl_start,
				  metaexons);
		  
		  metaexon_insert(last_cal->strand, last_cal->chromosome_id - 1,
				  node_avl_end->position, s_next_aux->genome_end, 40,
				  METAEXON_LEFT_END, node_avl_end,
				  metaexons);

		  meta_alignment_insert_cigar(cigar_code, CIGAR_SIMPLE_MIDDLE, c - 1, meta_alignment);

		} else { 	
		  //printf("2.PROCESS.2 $ %s\n", fq_read->id);
		  info_sp_t *info_sp = sw_reference_splice_junction(first_cal, last_cal,
								    query_map, genome,
								    q, r);
		  //New sw item. Storing data...
		  sw_item = sw_item_new(SP_SW, i, 0, c - 1,
					first_cal, last_cal, 
					meta_alignment, NULL, 
					NULL, info_sp);	      
		  //Insert item... and process if depth is full
		  sw_depth_insert(q, r, sw_item,
				  sw_optarg, output,
				  avls_list, metaexons, &sw_depth);	        
		}

		first_cal = last_cal;
	      } //End loop aux list 
	    }
	    meta_alignment_fill_gaps(meta_type,
				     meta_alignment, query_map, genome,
				     sw_optarg, output, metaexons, 
				     &sw_depth, avls_list);
	    array_list_insert(meta_alignment, meta_alignments_list[i]);
	  }
	  array_list_free(aux_list, NULL);
	  array_list_free(init_list, NULL);
	} //Not found splice in metaexon      
	//====(((((S M I T H - W A T E R M A N    P H A S E    E N D)))))====//

      } //End loop double anchors (cals_list)      

      if (array_list_size(meta_alignments_list[i]) == 0) {
	array_list_clear(mapping_batch->mapping_lists[i], NULL);
      }

    } else if (flag == SINGLE_ANCHORS) {
      //printf("\tWK_1ph: -- SINGLE ANCHOR PROOCESS --\n");
      //metaexons_show(metaexons);
      int map = 0;
      
      int read_nt;
      int seeds_process = 0;
      for (int p = 0; p < array_list_size(cals_list); p++) {
	cal = array_list_get(p, cals_list);
	//cal_print(cal);
	make_seeds = 0;
	if (cal->strand == 1) {
	  if (!rev_comp[i]) {
	    rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
	    strcpy(rev_comp[i], fq_read->sequence);
	    seq_reverse_complementary(rev_comp[i], fq_read->length);
	  }
	  query_map = rev_comp[i];
	} else {
	  query_map = fq_read->sequence;
	}
	
	s_prev = linked_list_get_last(cal->sr_list);
	if (s_prev->read_start == 0) {
	  read_nt = fq_read->length - (s_prev->read_end + 1);
	} else {
	  read_nt = s_prev->read_start;
	}
	
	//printf("CAL %i-%i\n", cal->start, cal->end);
	//metaexons_show(metaexons);
	if (metaexon_search(cal->strand, cal->chromosome_id - 1,
			    cal->start, cal->end, &first_metaexon,
			    metaexons) && first_metaexon) {	  
	  if (read_nt <= 20) {
	    meta_alignment = meta_alignment_new();
	    if (s_prev->read_start == 0) {
	      cigar_code = fill_extrem_gap(query_map, 
					   cal,
					   FILL_GAP_RIGHT,
					   genome,
					   first_metaexon,
					   metaexons);
	      if (cigar_code) {
		meta_alignment->cigar_right = cigar_code;
	      }
	    } else {
	      cigar_code = fill_extrem_gap(query_map, 
					   cal,
					   FILL_GAP_LEFT,
					   genome,
					   first_metaexon,
					   metaexons);
	      if (cigar_code) {
		meta_alignment->cigar_left = cigar_code;
	      }
	    }

	    meta_alignment_insert_cal(cal, meta_alignment);
	    meta_alignment_fill_gaps(META_ALIGNMENT_NONE,
				     meta_alignment, query_map, genome,
				     sw_optarg, output, metaexons, 
				     &sw_depth, avls_list);
	    array_list_insert(meta_alignment, meta_alignments_list[i]); 	  
	  } else {
	    int distance;
	    meta_alignment = NULL;
	    //printf("START WITH SEARCH IN METAEXONS!!\n");
	    array_list_clear(final_positions, NULL);
	    //cigar_code = NULL;
	    if (s_prev->read_start == 0) {
	      if (first_metaexon->right_closed) {
		//printf(" LEFT SEARCH\n");
		cigar_code = search_left_single_anchor(read_nt, 
						       cal,
						       0,
						       first_metaexon->right_breaks,
						       query_map, metaexons, genome);
	      } else {
		cigar_code = fill_extrem_gap(query_map, 
					     cal,
					     FILL_GAP_RIGHT,
					     genome,
					     first_metaexon,
					     metaexons);	   
	      }

	      if (cigar_code != NULL) {
		//printf(" :::: LEFT CIGAR %s", new_cigar_code_string(cigar_code));
		meta_alignment = meta_alignment_new();
		meta_alignment_insert_cal(cal, meta_alignment);
		meta_alignment_insert_cigar(cigar_code, CIGAR_ANCHOR_LEFT, 0, meta_alignment);
		meta_alignment_fill_gaps(META_ALIGNMENT_LEFT,
					 meta_alignment, query_map, genome,
					 sw_optarg, output, metaexons, 
					 &sw_depth, avls_list);	
	      }
	    } else {
	      if (first_metaexon->left_closed) { 
		//printf(" RIGHT SEARCH \n");
		cigar_code = search_right_single_anchor(read_nt,
							cal, 
							0,
							first_metaexon->left_breaks,
							query_map, metaexons, genome);
	      } else {
		cigar_code = fill_extrem_gap(query_map, 
					     cal,
					     FILL_GAP_LEFT,
					     genome,
					     first_metaexon,
					     metaexons);	   
	      }

	      if (cigar_code != NULL) {
		meta_alignment = meta_alignment_new();
		meta_alignment_insert_cal(cal, meta_alignment);
		meta_alignment_insert_cigar(cigar_code, CIGAR_ANCHOR_RIGHT, 0, meta_alignment);
		meta_alignment_fill_gaps(META_ALIGNMENT_RIGHT,
					 meta_alignment, query_map, genome,
					 sw_optarg, output, metaexons, 
					 &sw_depth, avls_list);	
	      }
	    }

	    if (meta_alignment != NULL) {
	      //Report mapping
	      array_list_insert(meta_alignment, meta_alignments_list[i]); 	  	    
	    } 
	  }
	}
      }
      
      //---------------------------------//

      if (array_list_size(meta_alignments_list[i]) <= 0) {
	post_process_reads[i] = 1;
	buffer_item_insert_new_item(fq_read, cals_list, 
				    NULL, BITEM_SINGLE_ANCHORS, 
				    buffer, buffer_hc, 0);
	//printf("READ UNMAPPED %i\n", linked_list_size(buffer));
	array_list_clear(cals_targets, NULL);	
	array_list_clear(mapping_batch->mapping_lists[i], NULL);
      }
      
    } else if (flag == NOT_ANCHORS) {
      //printf("\tWK_1ph: -- NOT ANCHORS (CALS PROCESS) --\n");
      if (array_list_size(cals_list) <= 0) {
	post_process_reads[i] = 1;
	buffer_item_insert_new_item(fq_read, cals_list, 
				    NULL, BITEM_NO_CALS,
				    buffer, buffer_hc, 0);
	continue;
      }

      
      /*
      printf("==== ORIG CALS LIST %i ====\n", array_list_size(cals_list));
      for (int j = 0; j < array_list_size(cals_list); j++) {
      cal = array_list_get(j, cals_list);
      cal_print(cal);
      }
      printf("==== CALS LIST ====\n");      
      */


      //printf("NOT_ANCHORS FOUND %i CALs\n", array_list_size(cals_list)); 
      merge_and_filter_cals(cals_targets, cals_list, fq_read, 
			    bwt_optarg, bwt_index, genome, 
			    scores_ranking[i]);
      
      //printf("::===:: BEST SCORE %f WITH %i CALs ::===:: %s\n", scores_ranking[i][0], 
      ///array_list_size(cals_targets), fq_read->id);      
      const int MAX_PROCESS = 5;
      int meta_type;
      if (array_list_size(cals_targets) <= MAX_PROCESS) {
	order_cals(cals_list);
	int limit = MAX_PROCESS > array_list_size(cals_targets) ?
	  array_list_size(cals_targets) : MAX_PROCESS;	
	//printf("::----> Process LIM %i\n", limit);
	for (int j = 0; j < limit; j++) {
	  meta_alignment = NULL;
	  meta_type = -1;
	  array_list_t *fusion_list = array_list_get(j, cals_targets);

	  first_cal = array_list_get(0, fusion_list); 
	  if (first_cal->strand == 1) {
	    if (!rev_comp[i]) {
	      rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
	      strcpy(rev_comp[i], fq_read->sequence);
	      seq_reverse_complementary(rev_comp[i], fq_read->length);
	    }
	    query_map = rev_comp[i];
	  } else {
	    query_map = fq_read->sequence;
	  }
	  
	    //printf("==== FUSION CALS PROCESS ====\n");
	    //for (int t = 0; t < array_list_size(fusion_list); t++) {
	      //cal_t *cal_aux = array_list_get(t, fusion_list);
	    //cal_print(cal_aux);
	    //}
	    //printf("==== ------------------- ====\n");
	  if (array_list_size(fusion_list) > 1) { 
	    //printf("FUSION CALS REPORT\n");
	    cal_t *first_cal = array_list_get(0, fusion_list);
	    cal_t *last_cal  = array_list_get(array_list_size(fusion_list) - 1, fusion_list);
	    
	    //TODO: IF WE HAVE MORE THAN TWO CALS SEARCH SINGLE ANCHOR
	    cigar_code = search_double_anchors_cal(query_map,
						   first_cal, last_cal,
						   metaexons, genome,
						   fq_read, &meta_type);
	    if (cigar_code != NULL) {
	      //printf("FOUND! %s\n", new_cigar_code_string(cigar_code));
	      meta_alignment = meta_alignment_new();
	      meta_alignment_insert_cal(first_cal, meta_alignment);
	      if (meta_type == META_ALIGNMENT_LEFT) {
		meta_alignment_insert_cigar(cigar_code, CIGAR_ANCHOR_LEFT, 0, meta_alignment);
	      } else {
		//if (array_list_size(fusion_list) > 2) {
		//LOG_FATAL("MORE THAN 2 CALS MAPP IN DOUBLE ANCHORS\n");
		//}
		for (int c = 1; c < array_list_size(fusion_list); c++) {
		  cal = array_list_get(c, fusion_list);
		  meta_alignment_insert_cal(cal, meta_alignment);
		}
		meta_alignment_insert_cigar(cigar_code, CIGAR_SIMPLE_MIDDLE, 0, meta_alignment);
	      }
	      //array_list_free(fusion_list, NULL);
	      meta_alignment_fill_gaps(meta_type,
				       meta_alignment, query_map, genome,
				       sw_optarg, output, metaexons, 
				       &sw_depth, avls_list);	
	    } else {
	      //printf("NOT FOUND\n");
	      meta_alignment = meta_alignment_cals_new(fusion_list);  
	      for (int c = 1; c < array_list_size(fusion_list); c++) { 
		last_cal = array_list_get(c, fusion_list);
		avl_node_t *node_avl_start, *node_avl_end;
		seed_region_t *s_prev = linked_list_get_last(first_cal->sr_list);    
		seed_region_t *s_next = linked_list_get_first(last_cal->sr_list);	
		int distance_aux;
		int nt = search_simple_splice_junction(s_prev, s_next, 
						       last_cal->chromosome_id, 
						       last_cal->strand,
						       query_map, genome, 
						       avls_list, &node_avl_start,
						       &node_avl_end, &distance_aux);

		  if (nt) {
		    cigar_code = cigar_code_new();
		    cigar_code->distance = distance_aux;
		    cigar_code_append_new_op(nt, 'N', cigar_code);
		    
		    seed_region_t *s_prev_aux = linked_list_get_last(first_cal->sr_list);    
		    seed_region_t *s_next_aux = linked_list_get_first(last_cal->sr_list);
		    
		    metaexon_insert(first_cal->strand, first_cal->chromosome_id - 1,
				    s_prev_aux->genome_start, node_avl_start->position, 40,
				    METAEXON_RIGHT_END, node_avl_start,
				    metaexons);
		    
		    metaexon_insert(last_cal->strand, last_cal->chromosome_id - 1,
				    node_avl_end->position, s_next_aux->genome_end, 40,
				    METAEXON_LEFT_END, node_avl_end,
				    metaexons);
		    /*if (strcmp(fq_read->id,
			       "@ENST00000415752@ENSG00000241465@protein_coding@X@49216659@49223943@1@KNOWN_359_74_1_0_0_0_3:0:0_1:0:0_3#1789101-1") == 0) {
		      cal_print(first_cal);
		      cal_print(last_cal);
		      printf("  middle cigar : %s\n", 
			     new_cigar_code_string(cigar_code));
			     }*/
		    meta_alignment_insert_cigar(cigar_code, CIGAR_SIMPLE_MIDDLE, c - 1, meta_alignment);
		    
		  } else { 
		    info_sp_t *info_sp = sw_reference_splice_junction(first_cal, last_cal,
								      query_map, genome,
								      q, r);
		    //New sw item. Storing data...
		    sw_item = sw_item_new(SP_SW, i, 0, c - 1,
					  first_cal, last_cal, 
					  meta_alignment, NULL, 
					  NULL, info_sp);	      
		    //Insert item... and process if depth is full
		    sw_depth_insert(q, r, sw_item,
				    sw_optarg, output,
				    avls_list, metaexons, &sw_depth); 
		  }
		  first_cal = last_cal;
		}
		//array_list_free(fusion_cals, NULL);
		meta_alignment_fill_gaps(META_ALIGNMENT_MIDDLE,
					 meta_alignment, query_map, genome,
					 sw_optarg, output, metaexons, 
					 &sw_depth, avls_list);	 
	      }
	      array_list_insert(meta_alignment, meta_alignments_list[i]);
	    } else {
	      //printf(":::: SINGLE MAP (%i)::::\n", array_list_size(cals_list));
	      cal = array_list_get(0, fusion_list);
	      //cal_print(cal);

	      seed_region_t *s_prev = linked_list_get_first(cal->sr_list);
	      seed_region_t *s_next = linked_list_get_last(cal->sr_list);
	      int gap_start = s_prev->read_start;
	      int gap_end = s_next->read_end - fq_read->length - 1;	      
	      int process_left = 0, process_right = 0;

	      metaexon_insert(cal->strand, cal->chromosome_id - 1,
			      s_prev->genome_start, s_next->genome_end, 40,
			      METAEXON_NORMAL, NULL,
			      metaexons);
	      
	      meta_alignment = meta_alignment_new();
	      meta_alignment_insert_cal(cal, meta_alignment);
	      meta_alignment_fill_gaps(META_ALIGNMENT_MIDDLE,
				       meta_alignment, query_map, genome,
				       sw_optarg, output, metaexons, 
				       &sw_depth, avls_list);	 	      
	      array_list_insert(meta_alignment, meta_alignments_list[i]);		
	    }

	    cal_t *first_cal = array_list_get(0, fusion_list);
	    cal_t *last_cal  = array_list_get(array_list_size(fusion_list) - 1, fusion_list);	    
	    seed_region_t *s_prev = linked_list_get_first(first_cal->sr_list);
	    seed_region_t *s_next = linked_list_get_last(last_cal->sr_list);
	    cigar_code_t *cc_left, *cc_right;

	    if (s_prev->read_start != 0) {
	      //printf("LEFT SEARCH %i\n", first_cal->chromosome_id);
	      metaexon_search(first_cal->strand, first_cal->chromosome_id - 1,
			      first_cal->start, first_cal->end, &first_metaexon,
			      metaexons);
	      //printf("FINISH SEARCH\n");
	      cc_left = fill_extrem_gap(query_map, 
					first_cal,
					FILL_GAP_LEFT,
					genome,
					first_metaexon,
					metaexons); 
	      assert(meta_alignment != NULL);

	      /*if (!strcmp(fq_read->id, "@ENST00000372651@ENSG00000066136@protein_coding@1@41175215@41236612@1@KNOWN_353_701_0_1_0_0_3:0:0_1:0:0_0#115271-1")) {
		printf(" :: %s\n", new_cigar_code_string(cc_left));
		exit(-1);
		}*/

	      /*if (strcmp(fq_read->id,
			 "@ENST00000415752@ENSG00000241465@protein_coding@X@49216659@49223943@1@KNOWN_359_74_1_0_0_0_3:0:0_1:0:0_3#1789101-1") == 0) {  
		printf("%s :::::::.. LEFT > %s\n", new_cigar_code_string(meta_alignment->cigar_code), 
		       new_cigar_code_string(cc_left));
		       }*/
	      meta_alignment_insert_cigar(cc_left, CIGAR_ANCHOR_RIGHT, 0, meta_alignment);
	    }

	    if (meta_type != META_ALIGNMENT_LEFT &&
		s_next->read_end != fq_read->length - 1) {
	      //printf("RIGHT SEARCH\n");
	      metaexon_search(last_cal->strand, last_cal->chromosome_id - 1,
			      last_cal->start, last_cal->end, &first_metaexon,
			      metaexons);	    
	      cc_right = fill_extrem_gap(query_map, 
					 last_cal,
					 FILL_GAP_RIGHT,
					 genome,
					 first_metaexon,
					 metaexons); 
	      /*if (strcmp(fq_read->id,
			 "@ENST00000415752@ENSG00000241465@protein_coding@X@49216659@49223943@1@KNOWN_359_74_1_0_0_0_3:0:0_1:0:0_3#1789101-1") == 0) {  
	      printf("%s :::::::.. RIGHT > %s\n", new_cigar_code_string(meta_alignment->cigar_code), 
		     new_cigar_code_string(cc_right));
		     }*/
	      meta_alignment_insert_cigar(cc_right, CIGAR_ANCHOR_LEFT, 0, meta_alignment);
	    }
	}
      }


      if (array_list_size(meta_alignments_list[i]) <= 0) {
	post_process_reads[i] = 1;

	array_list_t *cals_aux = array_list_new(50, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
	for (int j = 0; j < array_list_size(cals_targets); j++) {
	  array_list_t *fusion_list = array_list_get(j, cals_targets);
	  for (int z = 0; z < array_list_size(fusion_list); z++) {
	    cal = array_list_get(z, fusion_list);
	    array_list_insert(cal, cals_aux);
	  }
	  array_list_free(fusion_list, NULL);
	}
	buffer_item_insert_new_item(fq_read, cals_aux, 
				    NULL, BITEM_CALS,
				    buffer, buffer_hc, 0);
	//buffer_item_t *buffer_item = buffer_item_complete_new(fq_read, cals_aux, SINGLE_ANCHORS);
	//linked_list_insert(buffer_item, buffer);
	//array_list_set_flag(BITEM_CALS, buffer_item->items_list);
	array_list_free(cals_aux, NULL);
	//printf("::==:: NO MAP ::==::\n");
      } else {
	for (int j = 0; j < array_list_size(cals_targets); j++) {
	  array_list_t *fusion_list = array_list_get(j, cals_targets);
	  array_list_free(fusion_list, NULL);
	}
      }

      array_list_clear(mapping_batch->mapping_lists[i], NULL);
      array_list_clear(cals_targets, NULL);

    } //else {
      //printf("\tWK_1ph: -- NOTHING --\n");
    //}

    //printf("WK_1ph: READ END PROCESS\n");
    array_list_clear(seeds_list, region_bwt_free);

  } //End loop reads
  
  sw_depth_process(sw_optarg, output, 
		   &sw_depth, avls_list, metaexons, SW_FINAL);


  //Merge splice junctions & order by distance
  //printf("CLOSE META ALIGNMENTS\n");
  for (int i = 0; i < num_reads; i++) {
    fq_read = array_list_get(i, mapping_batch->fq_batch);
    //printf("<<<<CLOSE META (%i) %s:\n", array_list_size(meta_alignments_list[i]), fq_read->id);

    if (array_list_size(meta_alignments_list[i]) == 0) { continue; }

    for (int m = 0; m < array_list_size(meta_alignments_list[i]); m++) {
      meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list[i]);
      if (meta_alignment_get_status(meta_alignment) == META_OPEN) {
	meta_alignment_close(meta_alignment);     
      }
    }

    //TODO: Alert! Reads with middle small exons in middle no map!Add 'I' operation in the cigar
    int no_map = 0;
    for (int m = 0; m < array_list_size(meta_alignments_list[i]); m++) {
      meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list[i]);
      /*if (strcmp(fq_read->id,
		 "@ENST00000415752@ENSG00000241465@protein_coding@X@49216659@49223943@1@KNOWN_359_74_1_0_0_0_3:0:0_1:0:0_3#1789101-1") == 0) {  
	cal_t *cal = array_list_get(0,  meta_alignment->cals_list);
	printf(" META CIGAR(%i): %s | middle %s\n",  cal->chromosome_id, 
	       new_cigar_code_string(meta_alignment->cigar_code), 
	       new_cigar_code_string(meta_alignment->middle_cigars[0]));
	//exit(-1);
	}*/
      if (meta_alignment->score != fq_read->length) {
	no_map = 1;
	break;
      }
    }

    if (no_map) {
      meta_alignment_t *meta_alignment = array_list_get(0, meta_alignments_list[i]);
      meta_alignment->flag = 1;
    }
  }
  
  //printf("WK_1ph: =============== REPORT ALIGNMENTS =====================\n");
  for (int i = 0; i < num_reads; i++) {
    fq_read = array_list_get(i, mapping_batch->fq_batch);
    //printf("(%i) WK_1ph-Meta_Report:  %s >>>>\n", array_list_size(meta_alignments_list[i]), fq_read->id);
    
    if (array_list_size(meta_alignments_list[i]) > 0) {
      array_list_clear(mapping_batch->mapping_lists[i], NULL);
    } else { continue; }

    meta_alignment_t *meta_alignment = array_list_get(0, meta_alignments_list[i]);
    if (meta_alignment->flag == 0) {
      char query[2048];
      char quality[2048];
      int map = 0;
      for (int m = 0; m < array_list_size(meta_alignments_list[i]); m++) {
	meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list[i]);
	//if (meta_alignment_get_status(meta_alignment) == META_CLOSE) {	  
	optional_fields_length = 0;
	optional_fields = NULL;
	
	first_cal = meta_alignment_get_first_cal(meta_alignment);
	s_prev = linked_list_get_first(first_cal->sr_list);
	int h_left = s_prev->read_start;
	
	last_cal = meta_alignment_get_last_cal(meta_alignment);
	s_next = linked_list_get_last(last_cal->sr_list);
	int h_right = (fq_read->length - 1) - s_next->read_end;

	if (first_cal->strand == 1) {
	  query_map = rev_comp[i];
	} else {
	  query_map = fq_read->sequence;
	}

	cigar_code = meta_alignment->cigar_code;	
	assert(cigar_code != NULL);

	//printf("Read length %i - \n", fq_read->length);
	if (h_left > 0 && 
	    meta_alignment->cigar_left == NULL) {
	  //printf("H_LEFT ->  %i -> seed %i\n", h_left, s_prev->read_start);
	  array_list_insert_at(0, cigar_op_new(h_left, 'H'), cigar_code->ops);
	} else {
	  h_left = 0;
	}
	
	if (h_right > 0 &&
	    meta_alignment->cigar_right  == NULL) {
	  //printf("H_RIGHT ->  %i -> seed %i\n", h_right,  s_next->read_end);
	  array_list_insert(cigar_op_new(h_right, 'H'), cigar_code->ops);
	} else {
	  h_right = 0;
	}

	if (h_left > fq_read->length || h_left < 0) { exit(-1); }

	//printf("FINAL H_LEFT = %i, H_RIGHT = %i\n", h_left, h_right);
	int len_read = fq_read->length - h_left - h_right;
	memcpy(query, &query_map[h_left], len_read);
	query[len_read] = '\0';

	memcpy(quality, &fq_read->quality[h_left], len_read);
	quality[len_read] = '\0';

	//printf("* * * %s * * *\n", fq_read->id);
	alignment = alignment_new();

	int header_len = strlen(fq_read->id); 
	char *header_id[header_len + 1];
	get_to_first_blank(fq_read->id, header_len, header_id);
	char *header_match = (char *)malloc(sizeof(char)*header_len);
	if (header_match == NULL) { exit(-1); }
	memcpy(header_match, header_id, header_len);

	if (!cigar_code_validate_(fq_read, cigar_code)) {
	  char cigar_fake[512];
	  sprintf(cigar_fake, "%iM", fq_read->length);
	  //printf("WK_1ph: * * * M E T A    A L I G N M E N T    R E P O R T    F A K E %s* * *\n", cigar_fake);
	  printf("@FAKE :%s\n", fq_read->id);
	  meta_alignment_regenerate_cigar(meta_alignment, query_map, genome);
	  cigar_code = meta_alignment->cigar_code;

	  printf("##NEW_CIGAR: %s\n", new_cigar_code_string(cigar_code));

	  alignment_init_single_end(header_match, 
				    strdup(fq_read->sequence),
				    strdup(fq_read->quality),
				    first_cal->strand, first_cal->chromosome_id - 1, first_cal->start - 1,
				    strdup(new_cigar_code_string(cigar_code))/*strdup(cigar_fake)*/,
				    cigar_code_get_num_ops(cigar_code)/*1*/,
				    norm_score * 254, 1, (array_list_size(meta_alignments_list[i]) >= 1),
				    optional_fields_length, optional_fields, 0, alignment);	  
	} else {
	  //printf("META ALIGNMENT REPORT %i: %s\n", m, new_cigar_code_string(cigar_code));
	  //printf("WK_1ph: * * * M E T A    A L I G N M E N T    R E P O R T    %s* * *\n", new_cigar_code_string(cigar_code));
	  alignment_init_single_end(header_match, 
				    strdup(query),//match_seq
				    strdup(quality),//match_qual
				    first_cal->strand, first_cal->chromosome_id - 1, first_cal->start - 1,
				    strdup(new_cigar_code_string(cigar_code)),//strdup(cigar_fake)
				    cigar_code_get_num_ops(cigar_code),//1
				    norm_score * 254, 1, (array_list_size(meta_alignments_list[i]) >= 1),
				    optional_fields_length, optional_fields, 0, alignment);
	  //alignment_print(alignment);
	}
	array_list_insert(alignment, mapping_batch->mapping_lists[i]);
	map = 1;
	
	for (int o = 0; o < array_list_size(cigar_code->ops); o++) {
	  cigar_op_t *op = array_list_get(o, cigar_code->ops);
	  cigar_op_free(op);
	}

	for (int c = 0; c < array_list_size(meta_alignment->cals_list); c++) {
	  cal = array_list_get(c, meta_alignment->cals_list);	  
	  //linked_list_item_t *item_list = cal->sr_list->first, *item_prev = NULL;
	  while (s = linked_list_remove_first(cal->sr_list)) {
	    if (s->info != NULL) {
	      cigar_code = s->info;
	      array_list_clear(cigar_code->ops, cigar_op_free);
	      cigar_code_free(cigar_code);
	      s->info = NULL;
	    }
	    seed_region_free(s);
	  }
	  linked_list_free(cal->sr_list, NULL);
	  cal->sr_list = NULL;
	  if (cal->info != NULL) { 
	    cigar_code = cal->info;
	    array_list_clear(cigar_code->ops, cigar_op_free);
	    cigar_code_free(cigar_code); 
	  }
	  cal_free(cal);
	}
	array_list_free(meta_alignment->cals_list, NULL);
	
	if (meta_alignment->cigar_left != NULL) {
	  cigar_code = meta_alignment->cigar_left;
	  array_list_clear(cigar_code->ops, cigar_op_free);
	  cigar_code_free(cigar_code);
	}

	if (meta_alignment->cigar_right != NULL) {
	  cigar_code = meta_alignment->cigar_right;
	  array_list_clear(cigar_code->ops, cigar_op_free);
	  cigar_code_free(cigar_code);
	}

	for (int c = 0; c < meta_alignment->num_cigars; c++) {
	  cigar_code = meta_alignment->middle_cigars[c];
	  if (cigar_code != NULL) {
	    cigar_code = meta_alignment->middle_cigars[c];
	    array_list_clear(cigar_code->ops, cigar_op_free);
	    cigar_code_free(cigar_code); 
	  }
	}

	meta_alignment_free(meta_alignment);
	
      }
      
    } else {
      post_process_reads[i] = 1;
 
      buffer_item_insert_new_item(fq_read, meta_alignments_list[i], 
				  NULL, BITEM_META_ALIGNMENTS,
				  buffer, buffer_hc, 0);      
      //buffer_item_t *buffer_item = buffer_item_complete_new(fq_read, meta_alignments_list[i], SINGLE_ANCHORS);
      //linked_list_insert(buffer_item, buffer);
      //array_list_set_flag(BITEM_META_ALIGNMENTS, buffer_item->items_list);
    }    
  }

  //printf("WK_1ph: =============== REPORT ALIGNMENTS END =====================\n");

  array_list_t *new_fq_batch = array_list_new(num_reads, 
					      1.25f, 
					      COLLECTION_MODE_ASYNCHRONIZED);
  int num_new_reads = 0;
  for (i = 0; i < num_reads; i++) {
    fq_read = array_list_get(i, mapping_batch->fq_batch);
    array_list_t *alignments_list = mapping_batch->mapping_lists[i];
    if (post_process_reads[i] == 0) {
      mapping_batch->mapping_lists[num_new_reads++] = alignments_list;
      array_list_insert(fq_read, new_fq_batch);
    } else {
      array_list_free(alignments_list, NULL);
    }
  }

  free(post_process_reads);
  array_list_free(mapping_batch->fq_batch, NULL);
  mapping_batch->fq_batch = new_fq_batch;

  array_list_free(seeds_list, NULL);
  array_list_free(final_positions, NULL);
  free(new_targets);
  sw_multi_output_free(output);
  array_list_free(cals_targets, NULL);

  for (i = 0; i < num_reads; i++) {
    array_list_free(meta_alignments_list[i], NULL);
    if (rev_comp) { free(rev_comp[i]); }
  }

  free(meta_alignments_list);
  free(rev_comp);

  LOG_DEBUG("========= SPLICE JUNCTION SEARCH END =========\n");

  //metaexons_show(metaexons);

  return RNA_POST_PAIR_STAGE;

}

int apply_rna_last(sw_server_input_t* input_p, batch_t *batch) {
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  sw_optarg_t *sw_optarg = &input_p->sw_optarg;
  genome_t *genome = input_p->genome_p;
  size_t num_reads = array_list_size(mapping_batch->fq_batch);
  size_t num_targets = mapping_batch->num_targets;
  metaexons_t *metaexons = input_p->metaexons;
  cal_optarg_t *cal_optarg = input_p->cal_optarg_p;
  avls_list_t *avls_list = input_p->avls_list;
  bwt_optarg_t *bwt_optarg = input_p->bwt_optarg_p;
  bwt_index_t *bwt_index = input_p->bwt_index_p;
  linked_list_t *buffer = input_p->buffer;
  linked_list_t *buffer_hc = input_p->buffer_hc;

  //fprintf(stderr, "APPLY RNA LAST START... %i\n", num_reads);

  array_list_t *cals_list, *fusion_cals, *fusion_cals_aux;
  cal_t *cal, *cal_prev, *cal_next, *first_cal, *last_cal;
  fastq_read_t *fq_read;
  linked_list_iterator_t itr;  
  seed_region_t *s, *s_prev, *s_next;
  cigar_code_t *cigar_code, *cigar_code_prev, *cigar_code_aux;
  cigar_code_t *alig_cigar_code;
  cigar_op_t *cigar_op_start, *cigar_op_end, *cigar_op, *cigar_op_prev, *cigar_op_aux;
  int *new_targets = (int *)calloc(mapping_batch->num_allocated_targets, sizeof(int));
  array_list_t *merge_cals;
  linked_list_t *linked_list;
  seed_region_t *seed_region;

  //float *cals_score = (float *)calloc(100, sizeof(float));
  float score;
  char reference[2048];
  char reference_prev[2048];
  char reference_next[2048];
  char reference_aux[2048];
  char query[2048];
  char query_revcomp[2048];
  alignment_t *alignment;
  char q[2048];
  char r[2048];

  char **rev_comp = (char **)calloc(num_reads, sizeof(char *));

  fusion_coords_t *extrem_coords[2*40*mapping_batch->num_allocated_targets];
  fusion_coords_t *sp_coords[2*40*mapping_batch->num_allocated_targets];
  cigar_code_t *extrem_cigars[2*40*mapping_batch->num_allocated_targets];
  cigar_code_t *sp_cigars[2*40*mapping_batch->num_allocated_targets];

  char *sequence;
  char *query_ref;
  char *quality_map, *query_map;
  //float scores_ranking[mapping_batch->num_allocated_targets][50];
  float *cals_score;
  char cigar_str[1024];

  cigar_op_t *first_op;
  char *match_seq, *match_qual, *optional_fields, *p;
  int match_start, match_len, optional_fields_length, AS;
  float norm_score;

  int delete_not_cigar;
  int padding_left, padding_right, len_query;  
  int num_sw = 0;
  int num_sw_sp = 0;
  int num_sw_ext = 0;
  size_t num_new_targets = 0;
  int coverage;
  int read_length;
  int cal_pos = 0;
  int number_of_best;
  //int num_extrem_ok = 0;
  //int num_sp_ok = 0;

  //Change to report most alignments
  int seed_err_size = 20;
  int n_alignments = 1;
  int target;
  int c, exact_nt;
  size_t genome_start, genome_end, genome_start2, genome_end2, read_start, read_end;
  int flank = 20;
  int n_delete;
  int sw_pos;
  int seeds_nt; 
  int e1;
  int e2;
  int ref_pos;
  int nt_2;
  
  int start_seeding, end_seeding;
  int lim_start, lim_end;

  int min_intron_size = 40;

  metaexon_t *first_metaexon, *last_metaexon;
  int l_found, r_found; 
  size_t num_cals;
  register size_t t;
  register int i, j;

  int meta_type;
  meta_alignment_t *meta_alignment;
  sw_item_t *sw_item;  
  sw_depth_t sw_depth;
  sw_depth.depth = 0;
  sw_multi_output_t *output = sw_multi_output_new(MAX_DEPTH);
  //array_list_t **sw_items_list = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
  //array_list_t **alignments_list = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
  array_list_t **meta_alignments_list = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
  array_list_t *seeds_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *final_positions = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED); 
  int make_seeds = 0;

  array_list_t *cals_targets = array_list_new(50, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  int *post_process_reads = (int *)calloc(num_reads, sizeof(int));
  float scores_ranking[num_reads][50];
  int read_nt;

  for (int i = 0; i < num_reads; i++) {
    meta_alignments_list[i] = array_list_new(20, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  }


  //1st - MAP THE EASY READS
  for (i = 0; i < num_reads; i++) {
    fq_read = array_list_get(i, mapping_batch->fq_batch);
    cals_list = mapping_batch->mapping_lists[i];
    //fprintf(stderr, "@@@@:::: %s\n", fq_read->id); 
    //fprintf(stderr, "@@@@:::: %s\n", fq_read->sequence); 
    //fprintf(stderr, "@@@@:::: +\n"); 
    //fprintf(stderr, "@@@@:::: %s\n", fq_read->quality); 
    //continue;
    
    /*
    if (strcmp(fq_read->id, "@ENST00000415752@ENSG00000241465@protein_coding@X@49216659@49223943@1@KNOWN_359_74_1_0_0_0_3:0:0_1:0:0_3#1789101-1") == 0) {
      printf("<<<Last Read(%i): %s \n", array_list_size(cals_list), fq_read->id);      
    }*/

    //printf("WK_2ph-Process: == (%i CALs)Read %s ==\n", array_list_size(cals_list), fq_read->id);
    if (array_list_size(cals_list) > 0) {
      if (array_list_get_flag(cals_list) == BITEM_SINGLE_ANCHORS) {
	//SINGLE ANCHORS FOUND
	/*if (strcmp(fq_read->id,
		   "@ENST00000415752@ENSG00000241465@protein_coding@X@49216659@49223943@1@KNOWN_359_74_1_0_0_0_3:0:0_1:0:0_3#1789101-1") == 0) {	    
	  printf("\tWK_2ph: -- SINGLE ANCHORS PROCESS --\n");
	  }*/
	//fprintf(stderr, "\tSINGLE ANCHOR\n");
	int map_read = 0, map = 0;      
	int read_nt;
	int seeds_process = 0;
	for (int p = 0; p < array_list_size(cals_list); p++) {
	  cal = array_list_get(p, cals_list);
	  //cal_print(cal);
	  make_seeds = 0;
	  if (cal->strand == 1) {
	    if (!rev_comp[i]) {
	      rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
	      strcpy(rev_comp[i], fq_read->sequence);
	      seq_reverse_complementary(rev_comp[i], fq_read->length);
	    }
	    query_map = rev_comp[i];
	  } else {
	    query_map = fq_read->sequence;
	  }

	  s_prev = linked_list_get_last(cal->sr_list);
	  if (s_prev->read_start == 0) {
	    read_nt = fq_read->length - (s_prev->read_end + 1);
	  } else {
	    read_nt = s_prev->read_start;
	  }
	
	  //printf("CAL %i-%i\n", cal->start, cal->end);
	  //metaexons_show(metaexons);
	  if (metaexon_search(cal->strand, cal->chromosome_id - 1,
			      cal->start, cal->end, &first_metaexon,
			      metaexons) && first_metaexon) {	  
	    if (read_nt <= 20) {
	      meta_alignment = meta_alignment_new();
	      if (s_prev->read_start == 0) {
		cigar_code = fill_extrem_gap(query_map, 
					     cal,
					     FILL_GAP_RIGHT,
					     genome,
					     first_metaexon,
					     metaexons);
		if (cigar_code) {
		  meta_alignment->cigar_right = cigar_code;
		}
	      } else {
		cigar_code = fill_extrem_gap(query_map, 
					     cal,
					     FILL_GAP_LEFT,
					     genome,
					     first_metaexon,
					     metaexons);
		if (cigar_code) {
		  meta_alignment->cigar_left = cigar_code;
		}
	      }

	      meta_alignment_insert_cal(cal, meta_alignment);
	      meta_alignment_fill_gaps(META_ALIGNMENT_NONE,
				       meta_alignment, query_map, genome,
				       sw_optarg, output, metaexons, 
				       &sw_depth, avls_list);
	      array_list_insert(meta_alignment, meta_alignments_list[i]); 	  
	    } else {
	      int distance;
	      meta_alignment = NULL;
	      //printf("START WITH SEARCH IN METAEXONS!!\n");
	      array_list_clear(final_positions, NULL);
	      //cigar_code = NULL;
	      if (s_prev->read_start == 0) {
		if (first_metaexon->right_closed) {
		  //printf(" LEFT SEARCH\n");
		  cigar_code = search_left_single_anchor(read_nt, 
							 cal,
							 0,
							 first_metaexon->right_breaks,
							 query_map, metaexons, genome);
		} else {
		  cigar_code = fill_extrem_gap(query_map, 
					       cal,
					       FILL_GAP_RIGHT,
					       genome,
					       first_metaexon,
					       metaexons);	   
		}

		if (cigar_code != NULL) {
		  //printf(" :::: LEFT CIGAR %s", new_cigar_code_string(cigar_code));
		  meta_alignment = meta_alignment_new();
		  meta_alignment_insert_cal(cal, meta_alignment);
		  meta_alignment_insert_cigar(cigar_code, CIGAR_ANCHOR_LEFT, 0, meta_alignment);
		  meta_alignment_fill_gaps(META_ALIGNMENT_LEFT,
					   meta_alignment, query_map, genome,
					   sw_optarg, output, metaexons, 
					   &sw_depth, avls_list);	
		}
	      } else {
		if (first_metaexon->left_closed) { 
		  //printf(" RIGHT SEARCH \n");
		  cigar_code = search_right_single_anchor(read_nt,
							  cal, 
							  0,
							  first_metaexon->left_breaks,
							  query_map, metaexons, genome);
		} else {
		  cigar_code = fill_extrem_gap(query_map, 
					       cal,
					       FILL_GAP_LEFT,
					       genome,
					       first_metaexon,
					       metaexons);	   
		}
		if (cigar_code != NULL) {
		  meta_alignment = meta_alignment_new();
		  meta_alignment_insert_cal(cal, meta_alignment);
		  meta_alignment_insert_cigar(cigar_code, CIGAR_ANCHOR_RIGHT, 0, meta_alignment);
		  meta_alignment_fill_gaps(META_ALIGNMENT_RIGHT,
					   meta_alignment, query_map, genome,
					   sw_optarg, output, metaexons, 
					   &sw_depth, avls_list);	
		}
	      }

	      if (meta_alignment != NULL) {
		//Report mapping
		array_list_insert(meta_alignment, meta_alignments_list[i]); 	  	    
	      }
	    }
	  }
	}

	if (array_list_size(meta_alignments_list[i]) <= 0) {	
	  //Make Seeding and Caling
	  //printf("NEW SEEDS\n");
	  //array_list_clear(cals_list, NULL);
	  array_list_t *new_cals_list = array_list_new(200, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
	  int seed_size = 16;
	  num_cals = bwt_generate_cals(fq_read->sequence, seed_size, bwt_optarg,
				       bwt_index, new_cals_list);
	  //filter-incoherent CALs
	  int founds[num_cals], found = 0;
	  for (size_t j = 0; j < num_cals; j++) {
	    founds[j] = 0;
	    cal = array_list_get(j, new_cals_list);
	    LOG_DEBUG_F("\tcal %i of %i: sr_list size = %i (cal->num_seeds = %i) %i:%lu-%lu\n", 
			j, num_cals, cal->sr_list->size, cal->num_seeds,
			cal->chromosome_id, cal->start, cal->end);
	    if (cal->sr_list->size > 0) {
	      int start = 0;
	      size_t genome_start = 0;
	      int  first = 1;
	      for (linked_list_item_t *list_item = cal->sr_list->first; list_item != NULL; list_item = list_item->next) {
		seed_region_t *s = list_item->item;
	      
		LOG_DEBUG_F("\t\t:: star %lu > %lu s->read_start\n", start, s->read_start);
		if (start > s->read_start || 
		    s->read_start >= s->read_end) {
		  LOG_DEBUG("\t\t\t:: remove\n");
		  found++;
		  founds[j] = 1;
		}
		if (!first && 
		    ((s->genome_start < genome_start) || 
		     (s->genome_start - genome_start) > 2*fq_read->length)) {
		  //printf("Remove (genome_start = %i s->genome_start = %i)\n", genome_start, s->genome_start);
		  //cal_print(cal);
		  found++;
		  founds[j] = 1;
		}

		first = 0;
		start = s->read_end + 1;
		genome_start = s->genome_end + 1;
	      }
	    } else {
	      found++;
	      founds[j] = 1;
	    }
	  }

	  array_list_t *cal_list_aux;
	  if (found) {
	    int min_seeds = 100000;
	    int max_seeds = 0;
	    cal_list_aux = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
	    for (size_t j = 0; j < num_cals; j++) {
	      if (!founds[j]) {
		cal = array_list_get(j, new_cals_list);
		cal->num_seeds = cal->sr_list->size;
		if (cal->num_seeds > max_seeds) max_seeds = cal->num_seeds;
		if (cal->num_seeds < min_seeds) min_seeds = cal->num_seeds;
		array_list_insert(cal, cal_list_aux);
		array_list_set(j, NULL, new_cals_list);
	      }
	    }
	    array_list_free(new_cals_list, (void *) cal_free);
	    num_cals = array_list_size(cal_list_aux);
	    new_cals_list = cal_list_aux;
	    if (num_cals) {
	      mapping_batch->mapping_lists[i] = new_cals_list;
	    }
	  }	

	  if (!num_cals) {
	    //array_list_set_flag(BITEM_NO_CALS, new_cals_list);
	    array_list_free(new_cals_list, cal_free);
	  } else {
	    array_list_set_flag(BITEM_CALS, new_cals_list);	  
	    cals_list = new_cals_list;
	  }
	}
      }
  
      if (array_list_get_flag(cals_list) == BITEM_CALS) {    
	//CALS FOUND
	//printf("\tWK_2ph: -- CALS PROCESS --\n");
	//fprintf(stderr, "\tCALS PROCESS\n");
	array_list_clear(cals_targets, NULL);

	order_cals(cals_list);

	/*if (!strcmp(fq_read->id, "@ENST00000372651@ENSG00000066136@protein_coding@1@41175215@41236612@1@KNOWN_353_701_0_1_0_0_3:0:0_1:0:0_0#115271-1")) {
	  printf("CALS IN PIPELINE-2\n");
	  }*/

	/*	printf("==== BEFORE ORDER. FUSION CALS PROCESS ====\n");
	for (int t = 0; t < array_list_size(cals_list); t++) {
	  cal_t *cal_aux = array_list_get(t, cals_list);
	  cal_print(cal_aux);
	}	    
	
	
	printf("==== AFTER ORDER. FUSION CALS PROCESS ====\n");
	for (int t = 0; t < array_list_size(cals_list); t++) {
	  cal_t *cal_aux = array_list_get(t, cals_list);
	  cal_print(cal_aux);
	}	    
	*/

	merge_and_filter_cals(cals_targets, cals_list, fq_read, 
			      bwt_optarg, bwt_index, genome, 
			      scores_ranking[i]);

	//int MAX_CALS_PROCESS = 10;
	size_t seed_size = 16;
	float best_score = scores_ranking[i][0];
	if (best_score < 40.0) {
	  //Penalty CALs with overlap seeds 
	  for (int j = 0; j < array_list_size(cals_targets); j++) {
	    array_list_t *fusion_list = array_list_get(j, cals_targets);
	    cal = array_list_get(0, fusion_list);
	    seed_region_t *s_prev = linked_list_get_first(cal->sr_list);    
	    seed_region_t *s_next = linked_list_get_last(cal->sr_list);	
	   int last_nt = fq_read->length % seed_size;
	    int cal_size = s_next->read_end - s_prev->read_start;
	    if (s_prev->read_start == (seed_size / 2) || 
		s_prev->read_end == ((seed_size / 2) + seed_size) || 
		s_prev->read_start == last_nt) {
	      scores_ranking[i][j] = 0;
	    } else {
	      if (last_nt != 0 &&
		  s_next->read_end == (fq_read->length - last_nt - 1) || 
		  s_next->read_start == (read_length - last_nt - seed_size)) {
		scores_ranking[i][j] = 0;
	      }
	    }
	    /*if (cal_size < seed_size*2 &&
		cal_size % seed_size != 0) {
	      scores_ranking[i][j] = 0;
	      }*/
	  }
	}

	int limit;
	//MAX_CALS_PROCESS > array_list_size(cals_targets) ? 
	//array_list_size(cals_targets) : MAX_CALS_PROCESS;
	
	//for (int j = 1; j < array_list_size(cals_targets); j++) {
	//if (scores_ranking[i][j] >= best_score - 30.0) {
	//  limit++;
	//}
	//}
	
	//if (limit > 40) { limit = 40; }	
	
	limit = array_list_size(cals_targets);

	//fprintf(stderr, "%i vs %i : %s\n", limit, array_list_size(cals_targets), 
	//	fq_read->id);
	for (int j = 0; j < limit; j++) {	
	  meta_type = -1;
	  array_list_t *fusion_list = array_list_get(j, cals_targets);
	  //printf("==== FUSION CALS PROCESS ====\n");
	  //for (int t = 0; t < array_list_size(fusion_list); t++) {
	  //  cal_t *cal_aux = array_list_get(t, fusion_list);
	  //  cal_print(cal_aux);
	  //}	    
	
	  
	  first_cal = array_list_get(0, fusion_list);
	  last_cal = array_list_get(array_list_size(fusion_list) - 1, fusion_list);
	  if (first_cal->strand == 1) {
	    if (!rev_comp[i]) {
	      rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
	      strcpy(rev_comp[i], fq_read->sequence);
	      seq_reverse_complementary(rev_comp[i], fq_read->length);
	    }
	    query_map = rev_comp[i];
	  } else {
	    query_map = fq_read->sequence;
	  }

	  s_prev = linked_list_get_first(first_cal->sr_list);    
	  s_next = linked_list_get_last(last_cal->sr_list);		    	  

	  /*if (strcmp(fq_read->id,
		     "@ENST00000415752@ENSG00000241465@protein_coding@X@49216659@49223943@1@KNOWN_359_74_1_0_0_0_3:0:0_1:0:0_3#1789101-1") == 0) {	    
	    printf("==== FUSION CALS PROCESS ====\n");
	    for (int t = 0; t < array_list_size(fusion_list); t++) {
	      cal_t *cal_aux = array_list_get(t, fusion_list);
	      cal_print(cal_aux);
	    }
	    }*/

	  /*printf("(%i)(%i:%lu-%lu) Process CAL %i s_prev_[%lu-%lu] s_next_[%lu-%lu]: %f\n", 
		 first_cal->strand, first_cal->chromosome_id, first_cal->start, first_cal->end,
		 j, scores_ranking[i][j],
	  	 s_prev->read_start, s_prev->read_end,
	  	 s_next->read_start, s_next->read_end);
	  */
	  if (scores_ranking[i][j] == 0) { continue; }

	  //printf("==== ------------------- ====\n");
	  if (array_list_size(fusion_list) > 1) { 
	    //printf("FUSION CALS REPORT\n");
	    cal_t *first_cal = array_list_get(0, fusion_list);
	    cal_t *last_cal  = array_list_get(array_list_size(fusion_list) - 1, fusion_list);

	    //TODO: IF WE HAVE MORE THAN TWO CALS SEARCH SINGLE ANCHOR
	    cigar_code = search_double_anchors_cal(query_map,
						   first_cal, last_cal,
						   metaexons, genome,
						   fq_read, &meta_type);
	    if (cigar_code != NULL) {
	      //printf("FOUND! %s\n", new_cigar_code_string(cigar_code));
	      meta_alignment = meta_alignment_new();
	      meta_alignment_insert_cal(first_cal, meta_alignment);
	      if (meta_type == META_ALIGNMENT_LEFT) {
		meta_alignment_insert_cigar(cigar_code, CIGAR_ANCHOR_LEFT, 0, meta_alignment);
	      } else {
		meta_alignment_insert_cal(last_cal, meta_alignment);
		meta_alignment_insert_cigar(cigar_code, CIGAR_SIMPLE_MIDDLE, 0, meta_alignment);
	      }

	      meta_alignment_fill_gaps(meta_type,
				       meta_alignment, query_map, genome,
				       sw_optarg, output, metaexons, 
				       &sw_depth, avls_list);	
	    } else {
	      //printf("NOT FOUND\n");
	      meta_alignment = meta_alignment_cals_new(fusion_list);  
	      for (int c = 1; c < array_list_size(fusion_list); c++) { 
		last_cal = array_list_get(c, fusion_list);
		avl_node_t *node_avl_start, *node_avl_end;
		seed_region_t *s_prev = linked_list_get_last(first_cal->sr_list);    
		seed_region_t *s_next = linked_list_get_first(last_cal->sr_list);	
		int distance_aux;
		if (s_prev->genome_start >= s_next->genome_start) {
		  for (int t = 0; t < array_list_size(fusion_list); t++) {
		    cal_t *cal_aux = array_list_get(t, fusion_list);
		    //cal_print(cal_aux);
		    fprintf(stderr, "[%i:%lu-%lu]\n", cal_aux->chromosome_id, cal_aux->start, cal_aux->end);
		  }

		  fprintf(stderr, "%s\n", fq_read->id);
		  exit(-1);
		}

		int nt = search_simple_splice_junction(s_prev, s_next, 
						       last_cal->chromosome_id, 
						       last_cal->strand,
						       query_map, genome, 
						       avls_list, &node_avl_start,
						       &node_avl_end, &distance_aux);
		if (nt) {
		  cigar_code = cigar_code_new();
		  cigar_code->distance = distance_aux;
		  cigar_code_append_new_op(nt, 'N', cigar_code);
		
		  seed_region_t *s_prev_aux = linked_list_get_last(first_cal->sr_list);    
		  seed_region_t *s_next_aux = linked_list_get_first(last_cal->sr_list);
		
		  metaexon_insert(first_cal->strand, first_cal->chromosome_id - 1,
				  s_prev_aux->genome_start, node_avl_start->position, 40,
				  METAEXON_RIGHT_END, node_avl_start,
				  metaexons);
		
		  metaexon_insert(last_cal->strand, last_cal->chromosome_id - 1,
				  node_avl_end->position, s_next_aux->genome_end, 40,
				  METAEXON_LEFT_END, node_avl_end,
				  metaexons);
		
		  meta_alignment_insert_cigar(cigar_code, CIGAR_SIMPLE_MIDDLE, c - 1, meta_alignment);
		
		} else {/*
		  info_sp_t *info_sp = sw_reference_splice_junction(first_cal, last_cal,
								    query_map, genome,
								    q, r);
		  //New sw item. Storing data...
		  sw_item = sw_item_new(SP_SW, i, 0, c - 1,
					first_cal, last_cal, 
					meta_alignment, NULL, 
					NULL, info_sp);	      
		  //Insert item... and process if depth is full
		  sw_depth_insert(q, r, sw_item,
				  sw_optarg, output,
				  avls_list, metaexons, &sw_depth); 
			*/
		  meta_alignment_insert_cigar(NULL, CIGAR_SIMPLE_MIDDLE, c - 1, meta_alignment);
		}

		first_cal = last_cal;

	      }
	      //array_list_free(fusion_cals, NULL);
	      meta_alignment_fill_gaps(META_ALIGNMENT_MIDDLE,
				       meta_alignment, query_map, genome,
				       sw_optarg, output, metaexons, 
				       &sw_depth, avls_list);	 
	    }
	  } else {
	    //printf(":::: SINGLE MAP (%i)::::\n", array_list_size(cals_list));
	    cal = array_list_get(0, fusion_list);
	    //cal_print(cal);
	  
	    meta_alignment = meta_alignment_new();
	    meta_alignment_insert_cal(cal, meta_alignment);
	    meta_alignment_fill_gaps(META_ALIGNMENT_MIDDLE,
				     meta_alignment, query_map, genome,
				     sw_optarg, output, metaexons, 
				     &sw_depth, avls_list); 
	  }
	
	  array_list_insert(meta_alignment, meta_alignments_list[i]); 
	
	  cal_t *first_cal = array_list_get(0, fusion_list);
	  cal_t *last_cal  = array_list_get(array_list_size(fusion_list) - 1, fusion_list);	    
	  seed_region_t *s_prev = linked_list_get_first(first_cal->sr_list);
	  seed_region_t *s_next = linked_list_get_last(last_cal->sr_list);
	  cigar_code_t *cc_left, *cc_right;
	  if (s_prev->read_start != 0) {
	    //printf("LEFT SEARCH %i\n", first_cal->chromosome_id);
	    metaexon_search(first_cal->strand, first_cal->chromosome_id - 1,
			    first_cal->start, first_cal->end, &first_metaexon,
			    metaexons);
	    //printf("FINISH SEARCH\n");
	    cc_left = fill_extrem_gap(query_map, 
				      first_cal,
				      FILL_GAP_LEFT,
				      genome,
				      first_metaexon,
				      metaexons); 
	    /*if (strcmp(fq_read->id,
		       "@ENST00000415752@ENSG00000241465@protein_coding@X@49216659@49223943@1@KNOWN_359_74_1_0_0_0_3:0:0_1:0:0_3#1789101-1") == 0) {	    
	      printf(" L:: %s\n", new_cigar_code_string(cc_left));
	      }*/

	    assert(meta_alignment != NULL);
	    meta_alignment_insert_cigar(cc_left, CIGAR_ANCHOR_RIGHT, 0, meta_alignment);
	  }
	
	  if (meta_type != META_ALIGNMENT_LEFT &&
	      s_next->read_end != fq_read->length - 1) {
	    //printf("RIGHT SEARCH\n");
	    metaexon_search(last_cal->strand, last_cal->chromosome_id - 1,
			    last_cal->start, last_cal->end, &first_metaexon,
			    metaexons);	    
	    cc_right = fill_extrem_gap(query_map, 
				       last_cal,
				       FILL_GAP_RIGHT,
				       genome,
				       first_metaexon,
				       metaexons); 
	    /*if (strcmp(fq_read->id,
		       "@ENST00000415752@ENSG00000241465@protein_coding@X@49216659@49223943@1@KNOWN_359_74_1_0_0_0_3:0:0_1:0:0_3#1789101-1") == 0) {	    
	      printf(" R:: %s\n", cc_right);
	      }*/

	    meta_alignment_insert_cigar(cc_right, CIGAR_ANCHOR_LEFT, 0, meta_alignment);
	  }
	}
      } else if (array_list_get_flag(cals_list) == BITEM_META_ALIGNMENTS) {
	/*if (!strcmp(fq_read->id, "@ENST00000372651@ENSG00000066136@protein_coding@1@41175215@41236612@1@KNOWN_353_701_0_1_0_0_3:0:0_1:0:0_0#115271-1")) {
	  printf("META-ALIG IN PIPELINE-2\n");
	  }*/

	//META ALIGNMENTS
	//printf("\tWK_2ph: -- META-ALIGNMENTS PROCESS -- \n");
	//fprintf(stderr, "\tMETA ALIGNMENTS\n");
	for (int j = 0; j < array_list_size(cals_list); j++) {
	  meta_alignment = array_list_get(j, cals_list);

	  /*if (strcmp(fq_read->id,
		     "@ENST00000415752@ENSG00000241465@protein_coding@X@49216659@49223943@1@KNOWN_359_74_1_0_0_0_3:0:0_1:0:0_3#1789101-1") == 0) {	    
	    array_list_t *fusion_list = meta_alignment->cals_list;
	    printf("META ALIGNMENT FOUND \n");
	    printf("==== FUSION CALS PROCESS ====\n");
	    for (int t = 0; t < array_list_size(fusion_list); t++) {
	      cal_t *cal_aux = array_list_get(t, fusion_list);
	      cal_print(cal_aux);
	    }	    
	    }*/

	  //printf("SCORE IN : %i\n", meta_alignment->score);
	  array_list_t *fusion_list = meta_alignment->cals_list;
	  cal_t *first_cal = array_list_get(0, fusion_list);
	  cal_t *last_cal  = array_list_get(array_list_size(fusion_list) - 1, fusion_list);	    
	  seed_region_t *s_prev = linked_list_get_first(first_cal->sr_list);
	  seed_region_t *s_next = linked_list_get_last(last_cal->sr_list);
	  cigar_code_t *cc_left, *cc_right;

	  if (first_cal->strand == 1) {
	    if (!rev_comp[i]) {
	      rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
	      strcpy(rev_comp[i], fq_read->sequence);
	      seq_reverse_complementary(rev_comp[i], fq_read->length);
	    }
	    query_map = rev_comp[i];
	  } else {
	    query_map = fq_read->sequence;
	  }

	  if (s_prev->read_start != 0 && 
	      meta_alignment->cigar_left == NULL) {
	    metaexon_search(first_cal->strand, first_cal->chromosome_id - 1,
			    first_cal->start, first_cal->end, &first_metaexon,
			    metaexons);
	    cc_left = fill_extrem_gap(query_map, 
				      first_cal,
				      FILL_GAP_LEFT,
				      genome,
				      first_metaexon,
				      metaexons); 
	    assert(meta_alignment != NULL);
	    /*if (strcmp(fq_read->id,
		       "@ENST00000415752@ENSG00000241465@protein_coding@X@49216659@49223943@1@KNOWN_359_74_1_0_0_0_3:0:0_1:0:0_3#1789101-1") == 0) {	    
	      printf("FINISH LEFT SEARCH %s\n", new_cigar_code_string(cc_left));
	      }*/
	    meta_alignment_insert_cigar(cc_left, CIGAR_ANCHOR_RIGHT, 0, meta_alignment);
	  }
	
	  if (s_next->read_end != fq_read->length - 1 &&
	      meta_alignment->cigar_right == NULL) {
	    metaexon_search(last_cal->strand, last_cal->chromosome_id - 1,
			    last_cal->start, last_cal->end, &first_metaexon,
			    metaexons);	    
	    cc_right = fill_extrem_gap(query_map, 
				       last_cal,
				       FILL_GAP_RIGHT,
				       genome,
				       first_metaexon,
				       metaexons); 
	    /*if (strcmp(fq_read->id,
		       "@ENST00000415752@ENSG00000241465@protein_coding@X@49216659@49223943@1@KNOWN_359_74_1_0_0_0_3:0:0_1:0:0_3#1789101-1") == 0) {	    
	      printf("FINISH RIGHT SEARCH %s\n", new_cigar_code_string(cc_left));
	      }*/

	    //printf("RESULT CIGAR %s\n", new_cigar_code_string(cc_right));
	    meta_alignment_insert_cigar(cc_right, CIGAR_ANCHOR_LEFT, 0, meta_alignment);
	  }	
	}
      } else {
	//NO CALS FOUND
	//printf("\tWK_2ph: -- NOT CALS FOUND --\n");
	//fprintf(stderr, "\tNO CALS\n");
      }
    }
  }

  sw_depth_process(sw_optarg, output, 
		   &sw_depth, avls_list, metaexons, SW_FINAL);

  
  //fprintf(stderr, "CLOSE META ALIGNMENTS\n");
  for (int i = 0; i < num_reads; i++) {
    fq_read = array_list_get(i, mapping_batch->fq_batch);
    if (array_list_get_flag(mapping_batch->mapping_lists[i]) != BITEM_META_ALIGNMENTS) {
      array_list_clear(mapping_batch->mapping_lists[i], NULL);    
      //printf("<<<<CLOSE META (%i) %s\n", array_list_size(meta_alignments_list[i]), fq_read->id);
      for (int m = 0; m < array_list_size(meta_alignments_list[i]); m++) {
	meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list[i]);
	//printf("Status %i == %i\n", meta_alignment_get_status(meta_alignment), META_OPEN);
	if (meta_alignment_get_status(meta_alignment) == META_OPEN ) {
	  meta_alignment_close(meta_alignment);
	  //printf("CIGAR CLOSE %i: %s\n", m, new_cigar_code_string(meta_alignment->cigar_code));
	}
      }
    } else {      
      for (int j  = array_list_size(mapping_batch->mapping_lists[i]) - 1; j >= 0; j--) {
	meta_alignment_t *meta_alignment = array_list_remove_at(j, mapping_batch->mapping_lists[i]);
	meta_alignment_close(meta_alignment);
	//printf("CLOSE META %i : %s\n", j, new_cigar_code_string(meta_alignment->cigar_code));
	meta_alignment_calculate_score(meta_alignment);
	//printf("\n2: SCORE-META %i\n", meta_alignment->score);
	array_list_insert(meta_alignment, meta_alignments_list[i]);
      }
    }
    
    int no_map = 0;
    for (int m = 0; m < array_list_size(meta_alignments_list[i]); m++) {
      meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list[i]);
      if (meta_alignment->score != fq_read->length) {
	no_map = 1;
	break;
      }
    }

    if (no_map) {
      meta_alignment_t *meta_alignment = array_list_get(0, meta_alignments_list[i]);
      meta_alignment->flag = 1;
    }
  }
  //Score < 60??? New Seeding...
  /*
  for (int i = 0; i < num_reads; i++) {
    int seeding = 1;
    for (int j = 0; j < array_list_size(meta_alignments_list[i]); j++) {
      meta_alignment_t *meta_alignment = array_list_get(j, meta_alignments_list[i]);      
      cal_t *first_cal = array_list_get(0, meta_alignment->cals_list);
      cal_t *last_cal = array_list_get(array_list_size(meta_alignment->cals_list) - 1, meta_alignment->cals_list);

      if (meta_alignment->score >= 60) {
	seeding = 0;
	break;
      }
    }

    if (seeding) {
      bwt_map_inexact_seed_by_region(char *seq, size_t seq_len,
				     size_t seq_start, size_t seq_end,
				     bwt_optarg_t *bwt_optarg,
				     bwt_index_t *index,
				     array_list_t *mapping_list,
				     int seed_id, int chromosome_id,
				     size_t start_pos,
				     size_t end_pos, int strand_target);
	bwt_map_inexact_seeds_by_region(s_next->read_end, fq_read->length - 1,
					cal->chromosome_id, cal->end,
					cal->end + 500000,
					fq_read->sequence, seed_err_size,
					seed_err_size/2,
					bwt_optarg,
					bwt_index,
					mapping_list_next);
    }
  }*/

  //printf("WK_2ph: =============== REPORT ALIGNMENTS =====================\n");
  //fprintf(stderr, "REPORT META ALIGNMENTS\n");
  for (int i = 0; i < num_reads; i++) {
    if (array_list_size(meta_alignments_list[i]) <= 0) { continue; }

    fq_read = array_list_get(i, mapping_batch->fq_batch);
    assert(fq_read->id != NULL);

    //if (strcmp(fq_read->id,
    //	       "@ENST00000415752@ENSG00000241465@protein_coding@X@49216659@49223943@1@KNOWN_359_74_1_0_0_0_3:0:0_1:0:0_3#1789101-1") == 0) {  
    //exit(-1);
    //}

    //fprintf(stderr, ".( %i )Read %i: %s .\n", array_list_size(meta_alignments_list[i]), i, fq_read->id);
    //printf(".( %i )Read %i: %s .\n", array_list_size(meta_alignments_list[i]), i, fq_read->id);

    char query[2048];
    char quality[2048];
    int map = 0;

    //meta_alignment_t *meta_alignment = array_list_get(0, meta_alignments_list[i]);
    //if (meta_alignment->score == fq_read->length) {
    meta_alignment_t *meta_alignment = array_list_get(0, meta_alignments_list[i]);      
    //printf("WK_2ph-Meta_Report:  %s >>>>\n", fq_read->id);
    if (meta_alignment->flag == 0) {
      for (int m = 0; m < array_list_size(meta_alignments_list[i]); m++) {
	//printf("Report meta %i\n", m);
	meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list[i]);      
	if (meta_alignment_get_status(meta_alignment) == META_CLOSE) { 
	  //if (fq_read->id[1] == 'r') { printf("::--> %s\n", fq_read->id); exit(-1); }
	  optional_fields_length = 0;
	  optional_fields = NULL;
	  
	  first_cal = meta_alignment_get_first_cal(meta_alignment);
	  //printf("= SHOWING CALS  =\n");
	  //cal_print(first_cal);
	  //printf("= SHOWING ENDS  =\n");
	  
	  s_prev = linked_list_get_first(first_cal->sr_list);
	  int h_left = s_prev->read_start;
	    
	  last_cal = meta_alignment_get_last_cal(meta_alignment);
	  s_next = linked_list_get_last(last_cal->sr_list);
	  int h_right = (fq_read->length - 1) - s_next->read_end;
	    
	  if (first_cal->strand == 1) {
	    query_map = rev_comp[i];
	  } else {
	    query_map = fq_read->sequence;
	  }
	    
	  cigar_code = meta_alignment->cigar_code;	
	  assert(cigar_code != NULL);
	    
	  //fprintf(stderr, "Read length %i - \n", fq_read->length);
	  if (h_left > 0 && 
	      meta_alignment->cigar_left == NULL &&
	      meta_alignment->type != META_ALIGNMENT_LEFT &&
	      meta_alignment->type != META_ALIGNMENT_RIGHT) {
	    //fprintf(stderr, "H_LEFT ->  %i -> seed %i\n", h_left, s_prev->read_start);
	    array_list_insert_at(0, cigar_op_new(h_left, 'H'), cigar_code->ops);
	  } else {
	    h_left = 0;
	  }
	    
	  if (h_right > 0 &&
	      meta_alignment->cigar_right == NULL &&
	      meta_alignment->type != META_ALIGNMENT_LEFT &&
	      meta_alignment->type != META_ALIGNMENT_RIGHT) {
	    //fprintf(stderr, "H_RIGHT ->  %i -> seed %i\n", h_right,  s_next->read_end);
	    array_list_insert(cigar_op_new(h_right, 'H'), cigar_code->ops);
	  } else {
	    h_right = 0;
	  }
	    
	  //printf("FINAL H_LEFT = %i, H_RIGHT = %i\n", h_left, h_right);
	    
	  if (h_left > fq_read->length || h_left < 0) { exit(-1); }
	  if (h_left + h_right >= fq_read->length) { continue; }
	    
	  int len_read = fq_read->length - h_left - h_right;
	  //printf("(REAL %i): %s (%i)\n", len_read, query, s);
	  memcpy(query, &query_map[h_left], len_read);
	  query[len_read] = '\0';
	    
	  memcpy(quality, &fq_read->quality[h_left], len_read);
	  quality[len_read] = '\0';
	    
	  //printf("* * * %s * * *\n", fq_read->id);
	  //fprintf(stderr, "* * * (%i) M E T A    A L I G N M E N T    R E P O R T    %s* * *\n", strlen(query), new_cigar_code_string(cigar_code));
	  alignment = alignment_new();
	    
	  int header_len = strlen(fq_read->id); 
	  char *header_id[header_len + 1];
	  get_to_first_blank(fq_read->id, header_len, header_id);
	  char *header_match = (char *)malloc(sizeof(char)*header_len);
	  if (header_match == NULL) { exit(-1); }
	  memcpy(header_match, header_id, header_len);
	    
	  if (!cigar_code_validate(len_read, cigar_code)) {
	    char cigar_fake[512];
	    sprintf(cigar_fake, "%iM", fq_read->length);
	    printf("WK_2ph: * * * M E T A    A L I G N M E N T    R E P O R T    F A K E  [%i:%lu]  %s* * *\n", 
	    	   first_cal->chromosome_id, first_cal->start, new_cigar_code_string(cigar_code));
	    //fprintf(stderr, "@@@@@@@@@ :%s\n", fq_read->id);
	    printf("@FAKE :%s\n", fq_read->id);
	    alignment_init_single_end(header_match, 
				      strdup(fq_read->sequence),
				      strdup(fq_read->quality),
				      first_cal->strand, first_cal->chromosome_id - 1, first_cal->start - 1,
				      strdup(cigar_fake),
				      1,
				      norm_score * 254, 1, (array_list_size(meta_alignments_list[i]) >= 1),
				      optional_fields_length, optional_fields, 0, alignment);	  
	    //fprintf(stderr, "OK INSERT\n");
	  } else {
	    //fprintf(stderr, "META ALIGNMENT REPORT %i: %s\n", m, new_cigar_code_string(cigar_code));
	    //printf("WK_2ph: * * * M E T A    A L I G N M E N T    R E P O R T  [%i:%lu]  %s* * *\n", 
	    //	   first_cal->chromosome_id, first_cal->start, new_cigar_code_string(cigar_code));
	    alignment_init_single_end(header_match, 
				      strdup(query)/*match_seq*/,
				      strdup(quality)/*match_qual*/,
				      first_cal->strand, first_cal->chromosome_id - 1, first_cal->start - 1,
				      strdup(new_cigar_code_string(cigar_code))/*strdup(cigar_fake)*/,
				      cigar_code_get_num_ops(cigar_code)/*1*/,
				      norm_score * 254, 1, (array_list_size(meta_alignments_list[i]) >= 1),
				      optional_fields_length, optional_fields, 0, alignment);
	    //fprintf(stderr, "OK INSERT\n");
	    //alignment_print(alignment);
	  }
	  array_list_insert(alignment, mapping_batch->mapping_lists[i]);
	  //printf("Insert ok!\n");
	  map = 1;

	  for (int o = 0; o < array_list_size(cigar_code->ops); o++) {
	    cigar_op_t *op = array_list_get(o, cigar_code->ops);
	    cigar_op_free(op);
	  }
	  
	  for (int c = 0; c < array_list_size(meta_alignment->cals_list); c++) {
	    cal = array_list_get(c, meta_alignment->cals_list);	  
	    //linked_list_item_t *item_list = cal->sr_list->first, *item_prev = NULL;
	    while (s = linked_list_remove_first(cal->sr_list)) {
	      if (s->info != NULL) {
		cigar_code = s->info;
		array_list_clear(cigar_code->ops, cigar_op_free);
		cigar_code_free(cigar_code);
		s->info = NULL;
	      }
	      seed_region_free(s);
	    }
	    linked_list_free(cal->sr_list, NULL);
	    cal->sr_list = NULL;
	    if (cal->info != NULL) { 
	      cigar_code = cal->info;
	      array_list_clear(cigar_code->ops, cigar_op_free);
	      cigar_code_free(cigar_code); 
	    }
	    cal_free(cal);
	  }
	  array_list_free(meta_alignment->cals_list, NULL);
	  
	  if (meta_alignment->cigar_left != NULL) {
	    cigar_code = meta_alignment->cigar_left;
	    array_list_clear(cigar_code->ops, cigar_op_free);
	    cigar_code_free(cigar_code);
	  }
	  
	  if (meta_alignment->cigar_right != NULL) {
	    cigar_code = meta_alignment->cigar_right;
	    array_list_clear(cigar_code->ops, cigar_op_free);
	    cigar_code_free(cigar_code);
	  }
	  
	  for (int c = 0; c < meta_alignment->num_cigars; c++) {
	    cigar_code = meta_alignment->middle_cigars[c];
	    if (cigar_code != NULL) {
	      cigar_code = meta_alignment->middle_cigars[c];
	      array_list_clear(cigar_code->ops, cigar_op_free);
	      cigar_code_free(cigar_code); 
	    }
	  }
	  
	  meta_alignment_free(meta_alignment);
	}
      }
    } else {
      post_process_reads[i] = 1;
      buffer_item_insert_new_item(fq_read, meta_alignments_list[i], 
				  NULL, BITEM_META_ALIGNMENTS,
				  buffer, buffer_hc, 1);       
    }
  }


  //printf("WK_2ph: =============== REPORT ALIGNMENTS END =====================\n");

  array_list_t *new_fq_batch = array_list_new(num_reads, 
					      1.25f, 
					      COLLECTION_MODE_ASYNCHRONIZED);
  int num_new_reads = 0;
  for (i = 0; i < num_reads; i++) {
    fq_read = array_list_get(i, mapping_batch->fq_batch);
    array_list_t *alignments_list = mapping_batch->mapping_lists[i];
    if (post_process_reads[i] == 0) {
      mapping_batch->mapping_lists[num_new_reads++] = alignments_list;
      array_list_insert(fq_read, new_fq_batch);
    } else {
      array_list_free(alignments_list, NULL);
    }
  }

  free(post_process_reads);
  array_list_free(mapping_batch->fq_batch, NULL);
  mapping_batch->fq_batch = new_fq_batch;

  array_list_free(seeds_list, NULL);
  array_list_free(final_positions, NULL);
  free(new_targets);
  sw_multi_output_free(output);
  array_list_free(cals_targets, NULL);

  for (i = 0; i < num_reads; i++) {
    array_list_free(meta_alignments_list[i], NULL);
    if (rev_comp) { free(rev_comp[i]); }
  }

  free(meta_alignments_list);
  free(rev_comp);

  //fprintf(stderr, "APPLY RNA LAST END!\n");

  return LAST_RNA_POST_PAIR_STAGE;

}

int apply_rna_last_hc(sw_server_input_t* input_p, batch_t *batch) {
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  sw_optarg_t *sw_optarg = &input_p->sw_optarg;
  genome_t *genome = input_p->genome_p;
  size_t num_reads = array_list_size(mapping_batch->fq_batch);
  size_t num_targets = mapping_batch->num_targets;
  metaexons_t *metaexons = input_p->metaexons;
  cal_optarg_t *cal_optarg = input_p->cal_optarg_p;
  avls_list_t *avls_list = input_p->avls_list;
  bwt_optarg_t *bwt_optarg = input_p->bwt_optarg_p;
  bwt_index_t *bwt_index = input_p->bwt_index_p;
  linked_list_t *buffer = input_p->buffer;
  linked_list_t *buffer_hc = input_p->buffer_hc;

  //fprintf(stderr, "APPLY RNA LAST START... %i\n", num_reads);

  array_list_t *cals_list, *fusion_cals, *fusion_cals_aux;
  cal_t *cal, *cal_prev, *cal_next, *first_cal, *last_cal;
  fastq_read_t *fq_read;
  linked_list_iterator_t itr;  
  seed_region_t *s, *s_prev, *s_next;
  cigar_code_t *cigar_code, *cigar_code_prev, *cigar_code_aux;
  cigar_code_t *alig_cigar_code;
  cigar_op_t *cigar_op_start, *cigar_op_end, *cigar_op, *cigar_op_prev, *cigar_op_aux;
  int *new_targets = (int *)calloc(mapping_batch->num_allocated_targets, sizeof(int));
  array_list_t *merge_cals;
  linked_list_t *linked_list;
  seed_region_t *seed_region;

  //float *cals_score = (float *)calloc(100, sizeof(float));
  float score;
  char reference[2048];
  char reference_prev[2048];
  char reference_next[2048];
  char reference_aux[2048];
  char query[2048];
  char query_revcomp[2048];
  alignment_t *alignment;
  char q[2048];
  char r[2048];

  char **rev_comp = (char **)calloc(num_reads, sizeof(char *));

  fusion_coords_t *extrem_coords[2*40*mapping_batch->num_allocated_targets];
  fusion_coords_t *sp_coords[2*40*mapping_batch->num_allocated_targets];
  cigar_code_t *extrem_cigars[2*40*mapping_batch->num_allocated_targets];
  cigar_code_t *sp_cigars[2*40*mapping_batch->num_allocated_targets];

  char *sequence;
  char *query_ref;
  char *quality_map, *query_map;
  //float scores_ranking[mapping_batch->num_allocated_targets][50];
  float *cals_score;
  char cigar_str[1024];

  cigar_op_t *first_op;
  char *match_seq, *match_qual, *optional_fields, *p;
  int match_start, match_len, optional_fields_length, AS;
  float norm_score;

  int delete_not_cigar;
  int padding_left, padding_right, len_query;  
  int num_sw = 0;
  int num_sw_sp = 0;
  int num_sw_ext = 0;
  size_t num_new_targets = 0;
  int coverage;
  int read_length;
  int cal_pos = 0;
  int number_of_best;
  //int num_extrem_ok = 0;
  //int num_sp_ok = 0;

  //Change to report most alignments
  int seed_err_size = 20;
  int n_alignments = 1;
  int target;
  int c, exact_nt;
  size_t genome_start, genome_end, genome_start2, genome_end2, read_start, read_end;
  int flank = 20;
  int n_delete;
  int sw_pos;
  int seeds_nt; 
  int e1;
  int e2;
  int ref_pos;
  int nt_2;
  
  int start_seeding, end_seeding;
  int lim_start, lim_end;

  int min_intron_size = 40;

  metaexon_t *first_metaexon, *last_metaexon;
  int l_found, r_found; 
  size_t num_cals;
  register size_t t;
  register int i, j;

  int meta_type;
  meta_alignment_t *meta_alignment;
  sw_item_t *sw_item;  
  sw_depth_t sw_depth;
  sw_depth.depth = 0;
  sw_multi_output_t *output = sw_multi_output_new(MAX_DEPTH);
  //array_list_t **sw_items_list = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
  //array_list_t **alignments_list = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
  array_list_t *meta_alignments_list;
  array_list_t *seeds_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *final_positions = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED); 
  int make_seeds = 0;

  array_list_t *cals_targets = array_list_new(50, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *alignments_list;
  int *post_process_reads = (int *)calloc(num_reads, sizeof(int));
  float scores_ranking[num_reads][50];
  int read_nt;

  //Convert CALs in META_ALIGNMENTS 
  for (int i = 0; i < num_reads; i++) {
    fq_read = array_list_get(i, mapping_batch->fq_batch);
    //fprintf(stderr, "WK_3ph: %s\n", fq_read->id);

    cals_list = mapping_batch->mapping_lists[i];

    if (array_list_size(cals_list) <= 0) { continue; }

    if (array_list_get_flag(cals_list) != BITEM_META_ALIGNMENTS) {
      meta_alignments_list = array_list_new(array_list_size(cals_list),
					    1.25f, COLLECTION_MODE_ASYNCHRONIZED);

      array_list_clear(cals_targets, NULL);

      order_cals(cals_list);

      merge_and_filter_cals(cals_targets, cals_list, fq_read, 
			    bwt_optarg, bwt_index, genome, 
			    scores_ranking[i]);

      int limit = array_list_size(cals_targets);
      const int MAX_META_ALIGNMENTS = 10;

      if (limit > MAX_META_ALIGNMENTS) {
	/*for (int i = limit - 1; i >= MAX_META_ALIGNMENTS; i--) {
	  fusion_cals = array_list_get(j, cals_targets);
	  array_list_free(fusion_cals, cal_free);
	  }*/
	limit = MAX_META_ALIGNMENTS;
      }
      
      for (int j = 0; j < limit; j++) {
	fusion_cals = array_list_get(j, cals_targets);
	cal_t *first_cal = array_list_get(0, fusion_cals);
	cal_t *last_cal = array_list_get(array_list_size(fusion_cals) - 1, fusion_cals);

	if (first_cal->strand == 1) {
	  if (!rev_comp[i]) {
	    rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
	    strcpy(rev_comp[i], fq_read->sequence);
	    seq_reverse_complementary(rev_comp[i], fq_read->length);
	  }
	  query_map = rev_comp[i];
	} else {
	  query_map = fq_read->sequence;
	}

	meta_alignment = meta_alignment_cals_new(fusion_cals);
	meta_alignment_fill_gaps(META_ALIGNMENT_NONE,
				 meta_alignment, query_map, genome,
				 sw_optarg, output, metaexons, 
				 &sw_depth, avls_list);		
	//meta_alignment_close(meta_alignment);

	array_list_insert(meta_alignment, meta_alignments_list);
	s_prev = linked_list_get_first(first_cal->sr_list);
	s_next = linked_list_get_last(last_cal->sr_list);

	if (s_prev->read_start != 0) {
	  metaexon_search(first_cal->strand, first_cal->chromosome_id - 1,
			  first_cal->start, first_cal->end, &first_metaexon,
			  metaexons);
	  //printf("FINISH SEARCH\n");
	  cigar_code_t *cc_left = fill_extrem_gap(query_map, 
						  first_cal,
						  FILL_GAP_LEFT,
						  genome,
						  first_metaexon,
						  metaexons); 
	  assert(meta_alignment != NULL);
	  meta_alignment_insert_cigar(cc_left, CIGAR_ANCHOR_RIGHT, 0, meta_alignment);
	}
	
	if (s_next->read_end != fq_read->length - 1
	    /*meta_alignment->cigar_right == NULL*/) {
	  metaexon_search(last_cal->strand, last_cal->chromosome_id - 1,
			  last_cal->start, last_cal->end, &first_metaexon,
			    metaexons);	    
	    cigar_code_t *cc_right = fill_extrem_gap(query_map, 
						     last_cal,
						     FILL_GAP_RIGHT,
						     genome,
						     first_metaexon,
						     metaexons); 
	    
	    //printf("RESULT CIGAR %s\n", new_cigar_code_string(cc_right));
	    meta_alignment_insert_cigar(cc_right, CIGAR_ANCHOR_LEFT, 0, meta_alignment);
	} 	
      }

      mapping_batch->mapping_lists[i] = meta_alignments_list;

    } else {
      for (int m = 0; m < array_list_size(cals_list); m++) {
	meta_alignment = array_list_get(m, cals_list);
	fusion_cals = meta_alignment->cals_list;

	cal_t *first_cal = array_list_get(0, fusion_cals);
	cal_t *last_cal = array_list_get(array_list_size(fusion_cals) - 1, fusion_cals);

	if (first_cal->strand == 1) {
	  if (!rev_comp[i]) {
	    rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
	    strcpy(rev_comp[i], fq_read->sequence);
	    seq_reverse_complementary(rev_comp[i], fq_read->length);
	  }
	  query_map = rev_comp[i];
	} else {
	  query_map = fq_read->sequence;
	}

	s_prev = linked_list_get_first(first_cal->sr_list);
	s_next = linked_list_get_last(last_cal->sr_list);

	if (s_prev->read_start != 0 && 
	    meta_alignment->cigar_left == NULL) {
	  metaexon_search(first_cal->strand, first_cal->chromosome_id - 1,
			  first_cal->start, first_cal->end, &first_metaexon,
			  metaexons);

	  cigar_code_t *cc_left = fill_extrem_gap(query_map, 
						  first_cal,
						  FILL_GAP_LEFT,
						  genome,
						  first_metaexon,
						  metaexons); 
	  /*if (!strcmp(fq_read->id, "@ENST00000372651@ENSG00000066136@protein_coding@1@41175215@41236612@1@KNOWN_353_701_0_1_0_0_3:0:0_1:0:0_0#115271-1")) {
	    printf(" %s\n", new_cigar_code_string(cc_left));
	    exit(-1);
	    }*/
	  assert(meta_alignment != NULL);
	  meta_alignment_insert_cigar(cc_left, CIGAR_ANCHOR_RIGHT, 0, meta_alignment);
	}
	
	if (s_next->read_end != fq_read->length - 1 &&
	    meta_alignment->cigar_right == NULL) {
	  metaexon_search(last_cal->strand, last_cal->chromosome_id - 1,
			  last_cal->start, last_cal->end, &first_metaexon,
			  metaexons);	    
	  cigar_code_t *cc_right = fill_extrem_gap(query_map, 
						   last_cal,
						   FILL_GAP_RIGHT,
						   genome,
						   first_metaexon,
						   metaexons); 
	  
	  //printf("RESULT CIGAR %s\n", new_cigar_code_string(cc_right));
	  meta_alignment_insert_cigar(cc_right, CIGAR_ANCHOR_LEFT, 0, meta_alignment);
	}	
      }
    }

    meta_alignments_list = mapping_batch->mapping_lists[i];
    //SW for complete meta-alignments extrems
    for (int j = 0; j < array_list_size(meta_alignments_list); j++) {
      meta_alignment = array_list_get(j, meta_alignments_list);
      cal = array_list_get(0, meta_alignment->cals_list);

      /*array_list_t *fusion_cals = meta_alignment->cals_list;
      printf("==== FUSION CALS PROCESS ====\n");
      for (int t = 0; t < array_list_size(fusion_cals); t++) {
	cal_t *cal_aux = array_list_get(t, fusion_cals);
	cal_print(cal_aux);
	}	*/    

      //printf("Meta - %i:%lu-%lu\n", cal->chromosome_id, cal->start, cal->end);
      if (cal->strand == 1) {
	if (!rev_comp[i]) {
	  rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
	  strcpy(rev_comp[i], fq_read->sequence);
	  seq_reverse_complementary(rev_comp[i], fq_read->length);
	}
	query_map = rev_comp[i];
      } else {
	query_map = fq_read->sequence;
      }
      
      if (meta_alignment->cigar_left == NULL) {
	seed_region_t *seed_reg = linked_list_get_first(cal->sr_list);
	if (seed_reg->read_start >= 5) {
	  //SW
	  genome_start = seed_reg->genome_start - seed_reg->read_start ;
	  genome_end = seed_reg->genome_start - 1;
	  genome_read_sequence_by_chr_index(r, 0, 
					    cal->chromosome_id - 1,
					    &genome_start, &genome_end,
					    genome);    
	  memcpy(q, query_map, seed_reg->read_start);
	  q[seed_reg->read_start] = '\0';
	  
	  //printf("query     L ::: %s\n", q);
	  //printf("reference L ::: %s\n", r);
	  //New sw item. Storing data...
	  sw_item = sw_item_new(EXTREM_SW_LEFT, i, j, j,
				cal, cal, 
				meta_alignment, NULL, 
				NULL, NULL);	      
	  //Insert item... and process if depth is full
	  sw_depth_insert(q, r, sw_item,
			  sw_optarg, output,
			  avls_list, metaexons, &sw_depth); 
	} else {
	  cigar_code = cigar_code_new();
	  cigar_code_append_new_op(seed_reg->read_start, 'M', cigar_code);
	  meta_alignment->cigar_left = cigar_code;
	}
      }

      if (meta_alignment->cigar_right == NULL) {
	cal = array_list_get(array_list_size(meta_alignment->cals_list) - 1, meta_alignment->cals_list);
	seed_region_t *seed_reg = linked_list_get_last(cal->sr_list);
	if (seed_reg->read_end <= fq_read->length - 5) {
	  //SW
	  int r_gap = fq_read->length - seed_reg->read_end - 1;
	  genome_start = seed_reg->genome_end + 1;
	  genome_end = genome_start + r_gap;	  
	  genome_read_sequence_by_chr_index(r, 0, 
					    cal->chromosome_id - 1,
					    &genome_start, &genome_end,
					    genome);    
	  memcpy(q, query_map + seed_reg->read_end, r_gap);
	  q[r_gap] = '\0';
	  //printf("query     R ::: %s\n", q);
	  //printf("reference R ::: %s\n", r);

	  //New sw item. Storing data...
	  sw_item = sw_item_new(EXTREM_SW_RIGHT, i, j, j,
				cal, cal, 
				meta_alignment, NULL, 
				NULL, NULL);	      
	  //Insert item... and process if depth is full
	  sw_depth_insert(q, r, sw_item,
			  sw_optarg, output,
			  avls_list, metaexons, &sw_depth); 
	} else {
	  cigar_code = cigar_code_new();
	  cigar_code_append_new_op(fq_read->length - seed_reg->read_end - 1, 'M', cigar_code);
	  meta_alignment->cigar_right = cigar_code;
	}
      }
    }
  }

  sw_depth_process(sw_optarg, output, 
		   &sw_depth, avls_list, metaexons, SW_FINAL);

  const int MIN_SCORE_CAL = 50;
  //printf("=============== REPORT ALIGNMENTS =====================\n");
  //fprintf(stderr, "REPORT META ALIGNMENTS\n");
  //fprintf(stderr, "WK_3ph: =============== REPORT ALIGNMENTS (%i) =====================\n", num_reads);
  for (int i = 0; i < num_reads; i++) {
    fq_read = array_list_get(i, mapping_batch->fq_batch);
    printf("WK_3ph-Meta_Report (%i):  %s >>>>\n", array_list_size(mapping_batch->mapping_lists[i]),
    	   fq_read->id);
    meta_alignments_list = mapping_batch->mapping_lists[i];
    
    alignments_list = array_list_new(array_list_size(mapping_batch->mapping_lists[i]), 
				     1.25f, COLLECTION_MODE_ASYNCHRONIZED);

    if (array_list_size(meta_alignments_list) <= 0) { continue; }
    
    //printf("filter meta-alignments-1\n");

    assert(fq_read->id != NULL);
    //fprintf(stderr, ".( %i )Read %i: %s .\n", array_list_size(meta_alignments_list[i]), i, fq_read->id);
    //printf(".( %i )Read %i: %s .\n", array_list_size(meta_alignments_list[i]), i, fq_read->id);
    
    char query[2048];
    char quality[2048];
    int map = 0;

    for (int m = array_list_size(meta_alignments_list) - 1; m >= 0; m--) {
      meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list);
      meta_alignment_close(meta_alignment);
      if (meta_alignment_get_status(meta_alignment) != META_CLOSE) { 
	array_list_remove_at(m, meta_alignments_list);
	continue;
      }
      meta_alignment_calculate_score(meta_alignment);
    }

    if (array_list_size(meta_alignments_list) <= 0) { continue; }
    
    //printf("Order meta-alignments-1\n");

    //Order by score
    for (int m = 0; m < array_list_size(meta_alignments_list); m++) {
      meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list);
      for (int m1 = m + 1; m1 < array_list_size(meta_alignments_list); m1++) {
	meta_alignment_t *meta_alignment_next = array_list_get(m1, meta_alignments_list);
	if (meta_alignment_next->score > meta_alignment->score) {
	  array_list_swap(m, m1, meta_alignments_list);
	}
      }
    }

    int n_report = array_list_size(meta_alignments_list);
    meta_alignment_t *meta_alignment = array_list_get(0, meta_alignments_list);
    int best_score = array_list_size(meta_alignments_list);//meta_alignment->score;
    
    /*for (int m = 1; m < array_list_size(meta_alignments_list); m++) {
      meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list);
      if (meta_alignment->score >= best_score -  10) {
	n_report++;
      } 
      }*/

    /*if (!strcmp(fq_read->id, "@ENST00000372651@ENSG00000066136@protein_coding@1@41175215@41236612@1@KNOWN_353_701_0_1_0_0_3:0:0_1:0:0_0#115271-1")) {
      for (int m = 0; m < array_list_size(meta_alignments_list); m++) {
	meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list);
	array_list_t *fusion_cals = meta_alignment->cals_list;
	printf("==== META FUSION CALS PROCESS ====\n");
	for (int t = 0; t < array_list_size(fusion_cals); t++) {
	  cal_t *cal_aux = array_list_get(t, fusion_cals);
	  cal_print(cal_aux);
	}
	printf("==== . %s . ====\n", new_cigar_code_string(meta_alignment->cigar_code));
      }
      exit(-1);
      }*/

    //printf("N Report %i\n", n_report);

    for (int m = 0; m < n_report;/*array_list_size(meta_alignments_list);*/ m++) {
      //printf("Report meta %i\n", m);
      meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list);

      //printf("Score meta %i\n", meta_alignment->score);
      int cals_score = meta_alignment_get_cals_score(meta_alignment);
      int report = 0;
      
      if (cals_score < MIN_SCORE_CAL) {	      
	if (meta_alignment->score >= 40) {
	  report = 1;
	}
      } else if (meta_alignment->score >= 50) {
	report = 1;
      }
      
      //report = 1;
      //first_cal = meta_alignment_get_first_cal(meta_alignment);
      //printf("SCORE (%i)[%i:%lu]: %i\n", first_cal->strand, first_cal->chromosome_id, first_cal->start, meta_alignment->score);
      if (report) {	
	//printf("Report...!\n");
	map = 1;
	optional_fields_length = 0;
	optional_fields = NULL;
	
	first_cal = meta_alignment_get_first_cal(meta_alignment);
	if (first_cal->strand == 1) {
	  if (!rev_comp[i]) {
	    rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
	    strcpy(rev_comp[i], fq_read->sequence);
	    seq_reverse_complementary(rev_comp[i], fq_read->length);
	  }
	  query_map = rev_comp[i];
	} else {
	  query_map = fq_read->sequence;
	}

	cigar_code = meta_alignment->cigar_code;	
	assert(cigar_code != NULL);       	    

	s_prev = linked_list_get_first(first_cal->sr_list);
	if (meta_alignment->cigar_left == NULL && 
	    s_prev->read_start > 0) {
	  cigar_op_t *op = cigar_op_new(s_prev->read_start, 'H');
	  array_list_insert_at(0, op, cigar_code->ops);
	} 
	
	last_cal = array_list_get(array_list_size(meta_alignment->cals_list) - 1, meta_alignment->cals_list);
	s_next = linked_list_get_last(last_cal->sr_list);
	//printf("LAST SEED %i, report %i, %p, %p\n", s_next->read_end, fq_read->length - s_next->read_end - 1, s_prev, s_next);

	if (meta_alignment->cigar_right == NULL && 
	    s_next->read_end < fq_read->length - 1) {
	  cigar_op_t *op = cigar_op_new(fq_read->length - s_next->read_end - 1, 'H');
	  array_list_insert(op, cigar_code->ops);
	}
	
	if (first_cal->strand == 1) {
	  query_map = rev_comp[i];
	} else {
	  query_map = fq_read->sequence;
	}
		
	int h_left, h_right;

	cigar_op_t *first_op = array_list_get(0, cigar_code->ops);
	if (first_op->name == 'H') {	  
	  h_left = first_op->number;
	} else  {
	  h_left = 0;
	}
	
	cigar_op_t *last_op = array_list_get(array_list_size(cigar_code->ops) - 1, cigar_code->ops);
	if (last_op->name == 'H') { 
	  h_right = last_op->number;
	} else  {
	  h_right = 0;
	}

	//printf("FINAL H_LEFT = %i, H_RIGHT = %i\n", h_left, h_right);
	
	//if (h_left > fq_read->length || h_left < 0) { exit(-1); }
	//if (h_left + h_right >= fq_read->length) { continue; }
	
	int len_read = fq_read->length - h_left - h_right;
	//printf("(REAL %i): %s (%i)\n", len_read, query, s);
	
	//printf("* * * %s * * *\n", fq_read->id);
	//fprintf(stderr, "* * * (%i) M E T A    A L I G N M E N T    R E P O R T    %s* * *\n", strlen(query), new_cigar_code_string(cigar_code));
	alignment = alignment_new();
	
	int header_len = strlen(fq_read->id); 
	char *header_id[header_len + 1];
	get_to_first_blank(fq_read->id, header_len, header_id);
	char *header_match = (char *)malloc(sizeof(char)*header_len);
	if (header_match == NULL) { exit(-1); }
	memcpy(header_match, header_id, header_len);
	
	if (!cigar_code_validate(len_read, cigar_code)) {
	  char cigar_fake[512];
	  sprintf(cigar_fake, "%iM", fq_read->length);
	  fprintf(stderr, "WK_3ph: * * * M E T A    A L I G N M E N T    R E P O R T    F A K E    %s => [%i:%lu]* * *\n", 
	  	  fq_read->id, first_cal->chromosome_id, first_cal->start);
	  //fprintf(stderr, "@@@@@@@@@ :%s\n", fq_read->id);
	  alignment_init_single_end(header_match, 
				    strdup(fq_read->sequence),
				    strdup(fq_read->quality),
				    first_cal->strand, first_cal->chromosome_id - 1, first_cal->start - 1,
				    strdup(cigar_fake),
				    1,
				    norm_score * 254, 1, (array_list_size(meta_alignments_list) >= 1),
				    optional_fields_length, optional_fields, 0, alignment);	  
	  //fprintf(stderr, "OK INIT");
	} else {
	  //fprintf(stderr, "META ALIGNMENT REPORT %i: %s\n", m, new_cigar_code_string(cigar_code));
	  //fprintf(stderr, "WK_3ph: * * * M E T A    A L I G N M E N T    R E P O R T  [%i:%lu]  %s* * *\n", 
	  //	 first_cal->chromosome_id, first_cal->start, new_cigar_code_string(cigar_code));
	  //printf("h_left = %i, len_Read = %i\n", );
	  memcpy(query, &query_map[h_left], len_read);
	  query[len_read] = '\0';
	  
	  memcpy(quality, &fq_read->quality[h_left], len_read);
	  quality[len_read] = '\0';
	  
	  alignment_init_single_end(header_match, 
				    strdup(query)/*match_seq*/,
				    strdup(quality)/*match_qual*/,
				    first_cal->strand, first_cal->chromosome_id - 1, first_cal->start - 1,
				    strdup(new_cigar_code_string(cigar_code))/*strdup(cigar_fake)*/,
				    cigar_code_get_num_ops(cigar_code)/*1*/,
				    norm_score * 254, 1, (array_list_size(meta_alignments_list) >= 1),
				    optional_fields_length, optional_fields, 0, alignment);
	    //fprintf(stderr, "OK INSERT\n");
	    //alignment_print(alignment);
	}

	//printf("Report CIGAR OK!\n");
	
	array_list_insert(alignment, alignments_list);
	for (int o = 0; o < array_list_size(cigar_code->ops); o++) {
	  cigar_op_t *op = array_list_get(o, cigar_code->ops);
	  cigar_op_free(op);
	}
	
	for (int c = 0; c < array_list_size(meta_alignment->cals_list); c++) {
	  cal = array_list_get(c, meta_alignment->cals_list);	  
	  //linked_list_item_t *item_list = cal->sr_list->first, *item_prev = NULL;
	  while (s = linked_list_remove_first(cal->sr_list)) {
	    if (s->info != NULL) {
	      cigar_code = s->info;
	      array_list_clear(cigar_code->ops, cigar_op_free);
	      cigar_code_free(cigar_code);
	      s->info = NULL;
	    }
	    seed_region_free(s);
	  }
	  linked_list_free(cal->sr_list, NULL);
	  cal->sr_list = NULL;
	  if (cal->info != NULL) { 
	    cigar_code = cal->info;
	    array_list_clear(cigar_code->ops, cigar_op_free);
	    cigar_code_free(cigar_code); 
	  }
	  cal_free(cal);
	}
	array_list_free(meta_alignment->cals_list, NULL);
	
	if (meta_alignment->cigar_left != NULL) {
	  cigar_code = meta_alignment->cigar_left;
	  array_list_clear(cigar_code->ops, cigar_op_free);
	  cigar_code_free(cigar_code);
	}
	
	if (meta_alignment->cigar_right != NULL) {
	  cigar_code = meta_alignment->cigar_right;
	  array_list_clear(cigar_code->ops, cigar_op_free);
	  cigar_code_free(cigar_code);
	}
	
	for (int c = 0; c < meta_alignment->num_cigars; c++) {
	  cigar_code = meta_alignment->middle_cigars[c];
	  if (cigar_code != NULL) {
	    cigar_code = meta_alignment->middle_cigars[c];
	    array_list_clear(cigar_code->ops, cigar_op_free);
	    cigar_code_free(cigar_code); 
	  }
	}	  
	meta_alignment_free(meta_alignment); 
        
      }
    }
    //printf("Final Read\n");
    //if (!map) {
    //printf("##@@:: %s\n", fq_read->id);
    //printf("##@@:: %s\n", fq_read->sequence);
    //printf("##@@:: +\n");
    //printf("##@@:: %s\n", fq_read->quality);
    //}
    mapping_batch->mapping_lists[i] = alignments_list;
  }

  //fprintf(stderr, "WK_3ph: =============== REPORT ALIGNMENTS END =====================\n");

  return LAST_RNA_POST_PAIR_STAGE;

}


  /*
    else { 
    info_sp_t *info_sp = sw_reference_splice_junction(first_cal, last_cal,
    query_map, genome,
    q, r);
    //New sw item. Storing data...
    sw_item = sw_item_new(SP_SW, i, 0, c - 1,
    first_cal, last_cal, 
    meta_alignment, NULL, 
    NULL, info_sp);	      
    //Insert item... and process if depth is full
    sw_depth_insert(q, r, sw_item,
    sw_optarg, output,
    avls_list, metaexons, &sw_depth); 
    }*/
  /*
    printf("Merge filter ok!\n");
    for (int j = 0; j < array_list_size(cals_targets); j++) {
    array_list_t *fusion_list = array_list_get(j, cals_targets);
    //for (int p = 0; p < array_list_size(fusion_list); p++) {
    cal = array_list_get(0, fusion_list);
    printf("DIR REV_COMP = %p\n", rev_comp[i]);
    if (cal->strand == 1) {
    if (!rev_comp[i]) {
    rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
    strcpy(rev_comp[i], fq_read->sequence);
    seq_reverse_complementary(rev_comp[i], fq_read->length);
    }
    query_map = rev_comp[i];
    } else {
    query_map = fq_read->sequence;
    }
    
    
    cal_print(cal);
    meta_alignment = meta_alignment_new();
    meta_alignment_insert_cal(cal, meta_alignment);
    meta_alignment_fill_gaps(META_ALIGNMENT_NONE,
    meta_alignment, query_map, genome,
    sw_optarg, output, metaexons, 
    &sw_depth, avls_list);	
    array_list_insert(meta_alignment, meta_alignments_list[i]);  	
    //}
    }
    array_list_clear(cals_targets, NULL);
    fprintf(stderr, "\t>>>MAP NEW READ END, %s\n", fq_read->id);*/


      /*s_prev = linked_list_get_last(cal->sr_list);
	if (s_prev->read_start == 0) {
	  read_nt = fq_read->length - (s_prev->read_end + 1);
	} else {
	  read_nt = s_prev->read_start;
	}

	if (metaexon_search(cal->strand, cal->chromosome_id - 1,
			    cal->start, cal->end, &first_metaexon,
			    metaexons) && first_metaexon) {
	  int distance;
	  printf("START WITH SEARCH IN METAEXONS!!\n");
	  array_list_clear(final_positions, NULL);
	  //cigar_code = NULL;
	  meta_alignment = NULL;
	  if (s_prev->read_start == 0) {
	    if (first_metaexon->right_closed) {
	      //printf(" LEFT SEARCH\n");
	      cigar_code = search_left_single_anchor(read_nt, 
						     cal,
						     0,
						     first_metaexon->right_breaks,
						     query_map, metaexons, genome);
	      if (cigar_code != NULL) {
		//printf(" :::: LEFT CIGAR %s", new_cigar_code_string(cigar_code));
		meta_alignment = meta_alignment_new();
		meta_alignment_insert_cal(cal, meta_alignment);
		meta_alignment_insert_cigar(cigar_code, CIGAR_ANCHOR_LEFT, 0, meta_alignment);
		meta_alignment_fill_gaps(META_ALIGNMENT_LEFT,
					 meta_alignment, query_map, genome,
					 sw_optarg, output, metaexons, 
					 &sw_depth, avls_list);	
	      } 
	    } 
	  } else {
	    if (first_metaexon->left_closed) { 
	      //printf(" RIGHT SEARCH \n");
	      cigar_code = search_right_single_anchor(read_nt,
						    cal, 
						      0,
						      first_metaexon->left_breaks,
						      query_map, metaexons, genome);
	      if (cigar_code != NULL) {
		meta_alignment = meta_alignment_new();
		meta_alignment_insert_cal(cal, meta_alignment);
		meta_alignment_insert_cigar(cigar_code, CIGAR_ANCHOR_RIGHT, 0, meta_alignment);
		meta_alignment_fill_gaps(META_ALIGNMENT_RIGHT,
					 meta_alignment, query_map, genome,
					 sw_optarg, output, metaexons, 
					 &sw_depth, avls_list);	      
		//printf(" :::: RIGHT CIGAR %s\n", new_cigar_code_string(cigar_code));
	      }
	    }
	  }

	  if (meta_alignment != NULL) {
	    //Report mapping
	    map_read = 1;
	    array_list_insert(meta_alignment, meta_alignments_list[i]); 	  	    
	  } 
	}
      }
      */


      /*if (!map_read) {
	//------ First close gaps ---------//
	for (int c = 0; c < array_list_size(cals_list); c++) {
	  cal = array_list_get(c, cals_list);
	  s_prev = linked_list_get_first(cal->sr_list);
	  make_seeds = 0;
	  if (cal->strand == 1) {
	    if (!rev_comp[i]) {
	      rev_comp[i] = (char *) calloc(fq_read->length + 1, sizeof(char));
	      strcpy(rev_comp[i], fq_read->sequence);
	      seq_reverse_complementary(rev_comp[i], fq_read->length);
	    }
	    query_map = rev_comp[i];
	  } else {
	    query_map = fq_read->sequence;
	  }
	  int pos_err, distance;
	  int read_start, read_end;
	  size_t genome_start, genome_end;
	  int read_gap;


	  distance = fill_extrem_gap(query_map, 
				     cal->chromosome_id,
				     genome_start, genome_end,
				     read_start, read_end,
				     genome, read_gap,
				     &pos_err);
	  if (distance < 5) {
	    printf("(%i) COMPLETE DISTANCE %s\n", distance, fq_read->id);
	    meta_alignment = meta_alignment_new();
	    meta_alignment_insert_cal(cal, meta_alignment);
	    meta_alignment_fill_gaps(META_ALIGNMENT_NONE,
				     meta_alignment, query_map, genome,
				     sw_optarg, output, metaexons, 
				     &sw_depth, avls_list);	
	    array_list_insert(meta_alignment, meta_alignments_list[i]);  
	    map_read = 1;
	    }
	    }*/
	
	/*	printf("MAKE SEEDS %i\n", make_seeds);
	if (make_seeds) {
	  if (!seeds_process) {
	    seeds_process = 1;
	    int read_start, read_end;
	    if (cal->strand == 1) {
	      if (s_prev->read_start == 0) {
		read_start = 0;//fq_read->length - s_next->read_start;
		read_end = (fq_read->length - (s_prev->read_end + 1)) - 1;
	      } else {
		read_start = fq_read->length - s_prev->read_start;
		read_end = fq_read->length - 1;//(fq_read->length - (s_prev->read_end + 1)) - 1;
	      }
	    } else {
	      if (s_prev->read_start == 0) {
		read_start = s_prev->read_end + 1;
		read_end = fq_read->length - 1;
	      } else {
		read_start = 0;
		read_end = s_prev->read_start - 1;
	      }
	    }
	    array_list_clear(seeds_list, NULL);
	    printf("SEEDING BETWEEN (%i - %i)\n", read_start, read_end);
	    bwt_map_exact_seeds_by_region(read_start, read_end,
					  fq_read->sequence, 16, 16,
					  bwt_optarg, bwt_index,
					  seeds_list);
	  }
	  int caling_mode;
	  if (s_prev->read_start == 0) {
	    caling_mode = CALING_LEFT_ANCHORS;
	    last_cal = NULL;
	    first_cal = cal;
	  } else {
	    caling_mode = CALING_RIGHT_ANCHORS;
	    last_cal = cal;
	    first_cal = NULL;
	  }
	  array_list_t *aux_list = array_list_new(40, 1.25f, 
						  COLLECTION_MODE_ASYNCHRONIZED);	  
	  generate_cals_between_anchors(caling_mode,
					genome->num_chromosomes,
					first_cal, 
					last_cal, 
					fq_read, 
					seeds_list, 
					aux_list, 
					cal_optarg); 
	  if (array_list_size(aux_list) > 0) {
	    first_cal = array_list_get(0, aux_list);
	    //first_cal->fill = 1;
	    if (array_list_size(aux_list) == 1) {
	      meta_alignment = meta_alignment_cal_new(first_cal);
	      //meta_alignment_set_status(META_CLOSE, meta_alignment);
	      meta_type = META_ALIGNMENT_NONE;
	    } else {
	      meta_alignment = meta_alignment_cals_new(aux_list);		
	      meta_type = META_ALIGNMENT_MIDDLE;
	      for (int c = 1; c < array_list_size(aux_list); c++) { 
		last_cal = array_list_get(c, aux_list);
		if (last_cal == NULL) { exit(-1); }
		//last_cal->fill = 1;
		//Generate Fusion Reference...
		printf("\t->%i/%i (%i)[%i:%lu-%lu] vs (%i)[%i:%lu-%lu], %p, %p, %p, %p\n", c, array_list_size(aux_list), 
		       first_cal->strand, first_cal->chromosome_id, 
		       first_cal->start, first_cal->end, 
		       last_cal->strand, last_cal->chromosome_id,
		       last_cal->start, last_cal->end, first_cal, last_cal, first_cal->sr_list, last_cal->sr_list); 

		avl_node_t *node_avl_start, *node_avl_end;
		seed_region_t *s_prev_aux = linked_list_get_last(first_cal->sr_list);    
		seed_region_t *s_next_aux = linked_list_get_first(last_cal->sr_list);

		int nt = search_simple_splice_junction(s_prev_aux, s_next_aux,
						       first_cal->chromosome_id, 
						       first_cal->strand,
						       query_map, genome, 
						       avls_list, &node_avl_start,
						       &node_avl_end);
		if (nt) {
		  cigar_code = cigar_code_new();
		  cigar_code_append_new_op(nt, 'N', cigar_code);

		  metaexon_insert(first_cal->strand, first_cal->chromosome_id - 1,
				  s_prev_aux->genome_start, node_avl_start->position, 40,
				  METAEXON_RIGHT_END, node_avl_start,
				  metaexons);
		  
		  metaexon_insert(last_cal->strand, last_cal->chromosome_id - 1,
				  node_avl_end->position, s_next_aux->genome_end, 40,
				  METAEXON_LEFT_END, node_avl_end,
				  metaexons);

		  meta_alignment_insert_cigar(cigar_code, CIGAR_SIMPLE_MIDDLE, c - 1, meta_alignment);

		} else {		
		  //printf("2.PROCESS.2 $ %s\n", fq_read->id);
		  info_sp_t *info_sp = sw_reference_splice_junction(first_cal, last_cal,
								    query_map, genome,
								    q, r);
		  //New sw item. Storing data...
		  sw_item = sw_item_new(SP_SW, i, 0, c - 1,
					first_cal, last_cal, 
					meta_alignment, NULL, 
					NULL, info_sp);	      
		  //Insert item... and process if depth is full
		  sw_depth_insert(q, r, sw_item,
				  sw_optarg, output,
				  avls_list, metaexons, &sw_depth);	        
		}
		first_cal = last_cal;
	      } //End loop aux list 
	    }
	    meta_alignment_fill_gaps(meta_type,
				     meta_alignment, query_map, genome,
				     sw_optarg, output, metaexons, 
				     &sw_depth, avls_list);
	    array_list_insert(meta_alignment, meta_alignments_list[i]);
	  }
	  array_list_free(aux_list, NULL);
	  }*/


/*



	//norm_score = cigar_code_get_score(fq_read->length, cigar_code);
	//score = norm_score * 100;
	/*cal = array_list_get(0, meta_alignment->cals_list);
	  s = linked_list_get_first(cal->sr_list);
	  cigar_code = s->info;
	  cigar_code_print(cigar_code);
*/	  


/*

	  match_start = 0;
	  match_len = cigar_code_nt_length(cigar_code);
	  first_op = cigar_code_get_first_op(cigar_code);
	  match_start = (first_op && first_op->name == 'H' ? first_op->number : 0);
	  
	  match_seq = (char *) malloc((match_len + 1)* sizeof(char));
	  memcpy(match_seq, &fq_read->sequence[match_start], match_len);
	  match_seq[match_len] = 0;
	  
	  match_qual = (char *) malloc((match_len + 1)* sizeof(char));
	  memcpy(match_qual, &fq_read->quality[match_start], match_len);
	  match_qual[match_len] = 0;
	  */
	// set optional fields                                                             
	/*
	  optional_fields_length = 100;
	  optional_fields = (char *) calloc(optional_fields_length, sizeof(char));
	  
	  p = optional_fields;
	  AS = (int) norm_score * 100;
	  
	  sprintf(p, "ASi");
	  p += 3;
	  memcpy(p, &AS, sizeof(int));
	  p += sizeof(int);
	  
	  sprintf(p, "NHi");
	  p += 3;
	  memcpy(p, &num_cals, sizeof(int));
	  p += sizeof(int);
	  
	  sprintf(p, "NMi");
	  p += 3;
	  memcpy(p, &cigar_code->distance, sizeof(int));
	  p += sizeof(int);
	*/	

    /*
    int dst_prev, dst_next;
    if (array_list_size(meta_alignments_list[i]) >= 2) {
      for (int m0 = 0; m0 < array_list_size(meta_alignments_list[i]) - 1; m0++) {
	meta_alignment_t *ma_prev = array_list_get(m0, meta_alignments_list[i]);
	cal = meta_alignment_get_first_cal(ma_prev);
	s = linked_list_get_first(cal->sr_list);
	cigar_code_t *cigar_code = s->info;
	dst_prev = cigar_code->distance;
	for (int m1 = m0 + 1; m1 < array_list_size(meta_alignments_list[i]); m1++) {
	  meta_alignment_t *ma_next = array_list_get(m1, meta_alignments_list[i]);
	  cal = meta_alignment_get_first_cal(ma_next);
	  s = linked_list_get_first(cal->sr_list);
	  cigar_code_t *cigar_code = s->info;
	  dst_next = cigar_code->distance;
	  //printf("Meta %i/%i (%i) vs (%i)\n", m0, array_list_size(meta_alignments_list[i]) - 1, dst_prev, dst_next);
	  if (dst_next < dst_prev) {
	    array_list_swap(m0, m1, meta_alignments_list[i]);
	  }
	}	
      }
      }*/

  //Refresh batch, delete reads to post process
  /*int found;
  num_new_targets = 0;
  for (i = 0; i < num_reads; i++) {
    found = 0;
    if (array_list_size(meta_alignments_list[i]) > 0) {
      array_list_clear(mapping_batch->mapping_lists[i], NULL);
      for (int m = 0; m < array_list_size(meta_alignments_list[i]); m++) { 
	meta_alignment_t *meta_alignment = array_list_get(m, meta_alignments_list[i]);
	for (int c = 0; c < array_list_size(meta_alignment->cals_list); c++) {
	  cal = array_list_get(c, meta_alignment->cals_list);
	  if (cal->fill) {
	    found = 1;
	    array_list_insert(cal, mapping_batch->mapping_lists[i]);
	  }
	}
      }
    }
    if (found) {
      mapping_batch->targets[num_new_targets++] = i;
    }
    }*/

  //mapping_batch->num_targets = num_new_targets;
  //fill_gaps(mapping_batch, sw_optarg, genome, 15, 5);
  //merge_seed_regions(mapping_batch); 

  //Merge Cigars
  /*printf("MERGE CIGARS (%i)\n", num_new_targets);
  for (t = 0; t < num_new_targets; t++) {
    target = mapping_batch->targets[t];
    cals_list = mapping_batch->mapping_lists[target];
    fq_read = array_list_get(target, mapping_batch->fq_batch);    
    printf("<<<<READ UPDATE CIGARS %s\n", fq_read->id);
    num_cals = array_list_size(cals_list);
    for (int c = 0; c < num_cals; c++) {
      cal = array_list_get(c, cals_list);
      s = (seed_region_t *) linked_list_get_first(cal->sr_list);
      cigar_code = (cigar_code_t *) s->info;
      cigar_code_update(cigar_code);
      printf("\tCAL [%i:%lu-%lu]: \n", cal->chromosome_id, cal->start, cal->end);
      cigar_code_print(cigar_code);
    }
  }      
  */


	    /*if (gap >= 20) { 
	      printf("\tGAP SEARCHING... %i\n", gap);
	      //Seeding and caling between anchors cals...
	      int read_start, read_end;
	      if (first_cal->strand == 1) {
		read_start = fq_read->length - s_next->read_start + 1;
		read_end = fq_read->length - s_prev->read_end;
	      } else {
		read_start = s_prev->read_end + 1;
		read_end = s_next->read_start;
	      }
	      printf("\tNOT FOUND SEEDING... CALING... SW... anchor (%i-%i), (%i-%i)\n", s_prev->read_end, s_next->read_start, 
		     read_start, read_end);

	      region = region_bwt_new(first_cal->chromosome_id,
				      first_cal->strand,
				      first_cal->start, 
				      first_cal->end,
				      s_prev->read_start,
				      s_prev->read_end,
				      fq_read->length,
				      0);
	      array_list_insert(region, init_list);
	      region = region_bwt_new(last_cal->chromosome_id,
				      last_cal->strand,
				      last_cal->start, 
				      last_cal->end,
				      s_next->read_start,
				      s_next->read_end,
				      fq_read->length,
				      (gap / 16) + (gap % 16) + 1);
	      array_list_insert(region, init_list);

	      bwt_generate_cals_between_coords(first_cal->strand, first_cal->chromosome_id,
					       first_cal->end - 10, last_cal->start + 10,
					       read_start, read_end,
					       fq_read->sequence, 16, 16,
					       bwt_optarg, bwt_index,
					       init_list, aux_list);

	      if (array_list_size(aux_list) > 4) {
		printf("¡ALERT MAX CALS (%i)! %s\n", array_list_size(aux_list), fq_read->id);
		array_list_clear(aux_list, cal_free); 
	      }
	    } else {
	      array_list_insert(first_cal, aux_list);
	      array_list_insert(last_cal, aux_list);
	    }//End if gap < 20
	    */	    

/*
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
*/


/*
int generate_cals_score(array_list_t *fusion_cals, int read_length) {
  int score = 0;
  size_t num_cals = array_list_size(fusion_cals);
  cal_t *cal;
  linked_list_iterator_t itr;  
  seed_region_t *s;
  cigar_code_t *cigar_code;
  cigar_op_t *cigar_op_start, *cigar_op_end;

  if (num_cals == 1) {
    cal = array_list_get(0, fusion_cals);
    cigar_code = (cigar_code_t *)cal->info;
    if (!cigar_code) { return 0; }
    score = cigar_read_coverage(cigar_code);
    cigar_op_start = cigar_code_get_first_op(cigar_code);
    cigar_op_end = cigar_code_get_last_op(cigar_code);

    if (cigar_op_start->name != 'H' && 
	cigar_op_end->name != 'H') {
      score += score * read_length;
    }
 
    if (cigar_op_start->name == 'H' && 
	cigar_op_start->number < 20 && 
	cigar_op_end->name == 'H' && 
	cigar_op_end->number < 20) {
      score += score * read_length / 2;
    }

    return score;

  }

  for (int j = 0; j < num_cals; j++) {
    cal = array_list_get(j, fusion_cals);
    cigar_code = (cigar_code_t *)cal->info;
    if (cigar_code) {
      score += cigar_read_coverage(cigar_code);
      if (score > read_length ) { score = read_length; }
    } else {
      return 0;
    }
  }

  return score;

}
*/

/*
int apply_sw_rna_old(sw_server_input_t* input_p, batch_t *batch) {  
  //printf("RNA SW\n");
  size_t max_intron_size = 500000;  
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  sw_optarg_t *sw_optarg = &input_p->sw_optarg;
  genome_t *genome = input_p->genome_p;
  size_t num_targets = mapping_batch->num_targets;
  
  //allocate_splice_elements_t *chromosome_avls = input_p->chromosome_avls_p;
  avls_list_t *avls_list = input_p->avls_list;
  bwt_optarg_t *bwt_optarg = input_p->bwt_optarg_p;
  bwt_index_t *bwt_index = input_p->bwt_index_p;
  
  array_list_t *cals_list, *fusion_cals, *fusion_cals_aux;
  cal_t *cal, *cal_prev, *cal_next, *first_cal, *last_cal;
  fastq_read_t *fq_read;
  linked_list_iterator_t itr;  
  seed_region_t *s, *s_prev, *s_next;
  cigar_code_t *cigar_code, *cigar_code_prev, *cigar_code_aux;
  cigar_code_t *alig_cigar_code;
  cigar_op_t *cigar_op_start, *cigar_op_end, *cigar_op, *cigar_op_prev, *cigar_op_aux;
  int *new_targets = (int *)calloc(mapping_batch->num_allocated_targets, sizeof(int));
  array_list_t *merge_cals;
  linked_list_t *linked_list;
  seed_region_t *seed_region;

  //float *cals_score = (float *)calloc(100, sizeof(float));
  float score;
  char reference[2048];
  char reference_prev[2048];
  char reference_next[2048];
  char reference_aux[2048];
  char query[2048];
  char query_revcomp[2048];
  alignment_t *alignment;
  array_list_t **cals_targets = (array_list_t **)malloc(sizeof(array_list_t *)*mapping_batch->num_allocated_targets);
  char *q[2*40*mapping_batch->num_allocated_targets];
  char *r[2*40*mapping_batch->num_allocated_targets];

  fusion_coords_t *extrem_coords[2*40*mapping_batch->num_allocated_targets];
  fusion_coords_t *sp_coords[2*40*mapping_batch->num_allocated_targets];
  cigar_code_t *extrem_cigars[2*40*mapping_batch->num_allocated_targets];
  cigar_code_t *sp_cigars[2*40*mapping_batch->num_allocated_targets];

  char *sequence;
  char *query_ref;
  char *quality_map, *query_map;
  float scores_ranking[mapping_batch->num_allocated_targets][50];
  float *cals_score;
  char cigar_str[1024];

  int delete_not_cigar;
  int padding_left, padding_right, len_query;  
  int num_sw = 0;
  int num_sw_sp = 0;
  int num_sw_ext = 0;
  size_t num_new_targets = 0;
  int coverage;
  int read_length;
  int cal_pos = 0;
  int number_of_best;
  //int num_extrem_ok = 0;
  //int num_sp_ok = 0;

  //Change to report most alignments
  int seed_err_size = 20;
  int n_alignments = 1;
  int target;
  int c, exact_nt;
  size_t genome_start, genome_end, genome_start2, genome_end2, read_start, read_end;
  int flank = 20;
  int n_delete;
  int sw_pos;
  int seeds_nt; 
  int e1;
  int e2;
  int ref_pos;
  int nt_2;
  
  int start_seeding, end_seeding;
  int lim_start, lim_end;

  int min_intron_size = 40;

  register size_t num_cals;
  register size_t t;
  register int i, j;

  //flank = 20;
  flank = 0;
  for (t = 0; t < num_targets; t++) {
    target = mapping_batch->targets[t];
    cals_list = mapping_batch->mapping_lists[target];
    fq_read = array_list_get(target, mapping_batch->fq_batch);
    //printf("FROM RNA SERVER: %s\n", fq_read->id);
    //printf("List of CALs\n");
    //for (int i = 0; i < array_list_size(cals_list); i++)  {
    //cal = array_list_get(i, cals_list);
    //printf("[(%i)%i:%lu-%lu]\n", cal->strand, cal->chromosome_id, cal->start, cal->end);
    //}
    num_cals = array_list_size(cals_list);

    //strcpy(query_revcomp, fq_read->sequence);
    //seq_reverse_complementary(query_revcomp, fq_read->length);
    cals_targets[target] = array_list_new(num_cals,
					  1.25f,
					  COLLECTION_MODE_ASYNCHRONIZED);

    //array_list_clear(mapping_batch->mapping_lists[target], NULL);    
    number_of_best = merge_and_filter_cals(cals_targets[target], 
					   mapping_batch->mapping_lists[target], 
					   fq_read, n_alignments, bwt_optarg,
					   bwt_index, genome, scores_ranking[target], q, r, &num_sw, 
					   extrem_coords, target);     
  }

  float norm_score;
  float min_score = 0.6;
  int distance = 0;

  num_sw_ext = num_sw;
  sw_multi_output_t *output = sw_multi_output_new(num_sw);
  smith_waterman_mqmr(q, r, num_sw, sw_optarg, 1, output);
   
  //printf("\n-------------- SMITH-WATERMAN OUTPUTS EXTREMS --------------\n");
  for (i = 0;  i < num_sw_ext; i++) {
    //printf("Que:%s (Start:%i, Len:%i)\n", output->query_map_p[i], output->query_start_p[i], strlen(output->query_map_p[i]));
    //printf("Ref:%s (Start:%i, Len:%i)\n", output->ref_map_p[i],   output->ref_start_p[i],   strlen(output->ref_map_p[i]));

    score = NORM_SCORE(output->score_p[i], strlen(q[i]), sw_optarg->subst_matrix['A']['A']);

    fq_read = array_list_get(extrem_coords[i]->read_id, mapping_batch->fq_batch);
    cal = extrem_coords[i]->cal_ref;
    //printf("SW(%i/%i) Score %f [%i:%lu-%lu]\n", i, num_sw_ext, score, cal->chromosome_id, 
    //	   cal->start, cal->end);
    if (score >= 0.4) {
      //printf("#####LIMIT[%i:%lu-%lu] %s\n", extrem_coords[i]->cal_ref->chromosome_id, 
      //     extrem_coords[i]->cal_ref->start, extrem_coords[i]->cal_ref->end, extrem_coords[i]->id);
      distance = 0;
      cigar_code_aux = generate_cigar_code(output->query_map_p[i], output->ref_map_p[i], strlen(output->ref_map_p[i]),
					   output->query_start_p[i], output->ref_start_p[i],
					   strlen(q[i]), strlen(r[i]),
					   &distance, extrem_coords[i]->type_sw);

      //printf("\tSW SCORE AFTER %f (%s)\n", scores_ranking[extrem_coords[i]->read_id][extrem_coords[i]->cal_group_id], extrem_coords[i]->id);

      scores_ranking[extrem_coords[i]->read_id][extrem_coords[i]->cal_group_id] += 
	(((strlen(output->ref_map_p[i]) - distance)*100)/fq_read->length);

      //printf("\tSW SCORE NEW %f (%s)\n", scores_ranking[extrem_coords[i]->read_id][extrem_coords[i]->cal_group_id], extrem_coords[i]->id);

      extrem_cigars[i] = cigar_code_aux;
    } else {
      //printf("\tSW BAD %f (%s)\n", scores_ranking[extrem_coords[i]->read_id][extrem_coords[i]->cal_group_id], extrem_coords[i]->id);
      extrem_cigars[i] = NULL;
    }
    //printf("OLD Distance = %i, SW Distance = %i, SWItem = %i [%i:%lu-%lu]\n", cigar_code->distance, distance, 
    //     i, cal->chromosome_id, cal->start, cal->end);        
    free(q[i]);
    free(r[i]);
    //free(fusion_coords[i]);    
  }
  //printf("-------------- SMITH-WATERMAN OUTPUTS EXTREMS END --------------\n\n");

  sw_multi_output_free(output);


  num_sw = 0;
  for (t = 0; t < num_targets; t++) {
    target = mapping_batch->targets[t];
    cals_list = mapping_batch->mapping_lists[target];
    num_cals = array_list_size(cals_list);
    fq_read = array_list_get(target, mapping_batch->fq_batch);
    strcpy(query_revcomp, fq_read->sequence);
    seq_reverse_complementary(query_revcomp, fq_read->length);

    cals_score = scores_ranking[target];

    //Order the news scores
    for (i = 0; i < array_list_size(cals_targets[target]); i++) {
      for (j = i + 1; j < array_list_size(cals_targets[target]); j++) {
	if (cals_score[j] > cals_score[i]) {
	  array_list_swap(i, j, cals_targets[target]);
	  score = cals_score[j];
	  cals_score[j] = cals_score[i];
	  cals_score[i] = score;
	}
      }
    }

    //printf("FROM RNA SERVER 2: %s (Best Score %f)\n", fq_read->id, scores_ranking[target][0]);

    //for (i = 0; i < array_list_size(cals_targets[target]); i++) {
    //printf("\t SCORE %i:%f\n", cals_score[i]);
    //}

    generate_reference_splice_juntion(cals_targets[target], query_revcomp, 
				      fq_read, r, q, sp_coords, &num_sw,
				      genome, target);    
  }
  
  num_sw_sp = num_sw;
  output = sw_multi_output_new(num_sw);
  smith_waterman_mqmr(q, r, num_sw, sw_optarg, 1, output);
  
  LOG_DEBUG("O U T P U T \n");
  for (i = 0;  i < num_sw_sp; i++) {
    //printf("Que:%s (Start:%i, Len:%i)\n", output->query_map_p[i], output->query_start_p[i], strlen(output->query_map_p[i]));
    //printf("Ref:%s (Start:%i, Len:%i)\n", output->ref_map_p[i],   output->ref_start_p[i],   strlen(output->ref_map_p[i]));
    
    output->score_p[i] += 20.0;//Aprox 20nt*0.5(gap extend) + 10(open gap penalty)
    
    
    score = NORM_SCORE(output->score_p[i], strlen(q[i]), sw_optarg->subst_matrix['A']['A']);
    
    fq_read = array_list_get(sp_coords[i]->read_id, mapping_batch->fq_batch);
    if (score > min_score) {
      int tot_matches;
      flank = 30;
      cigar_code_aux = generate_cigar_sw_output(output->query_map_p[i],
						output->ref_map_p[i],
						sp_coords[i]->l_exon_start,
						sp_coords[i]->l_exon_end,
						sp_coords[i]->r_exon_start,
						sp_coords[i]->r_exon_end,
						sp_coords[i]->read_start,
						sp_coords[i]->read_end,
						sp_coords[i]->chromosome,
						sp_coords[i]->strand,
						output->query_start_p[i],
						output->ref_start_p[i],
						sp_coords[i]->type_sw, 
						strlen(q[i]),
						avls_list,
						sp_coords[i]->id, &tot_matches);
      
      scores_ranking[sp_coords[i]->read_id][sp_coords[i]->cal_group_id] -= (flank*100)/fq_read->length; 
      scores_ranking[sp_coords[i]->read_id][sp_coords[i]->cal_group_id] += (tot_matches*100)/fq_read->length; 
      sp_cigars[i] = cigar_code_aux;
    } else {
      sp_cigars[i] = NULL;
    }
    free(q[i]);
    free(r[i]);
  }
  
  sw_multi_output_free(output);

  fill_gaps(mapping_batch, sw_optarg, genome, 20, 5);
  merge_seed_regions(batch->mapping_batch); 

  for (i = 0; i < num_sw_ext; i++) {
    cal = extrem_coords[i]->cal_ref;
    //printf("Merge Cigars Extrems: %i/%i [%i:%lu-%lu] %s\n", i, num_sw_ext, cal->chromosome_id, cal->start, cal->end, extrem_coords[i]->id);
    s = (seed_region_t *) linked_list_get_first(cal->sr_list);
    cigar_code = (cigar_code_t *)s->info;
    cigar_code_aux = extrem_cigars[i];

    if (cigar_code_aux != NULL) {
      //printf("\tMerge\n");
      //Fusion Cigars
      if (extrem_coords[i]->type_sw == FIRST_SW) { 
	//printf("After loopp ops %i\n", array_list_size(cigar_code->ops));
	for (int c = 0; c < array_list_size(cigar_code->ops); c++) {
	  cigar_op = array_list_get(c, cigar_code->ops);
	  //printf("From Loop op: %i%c\n", cigar_op->number, cigar_op->name);
	  array_list_insert(cigar_op, cigar_code_aux->ops);
	}
	cigar_code_aux->distance += cigar_code->distance;
	cigar_code_free(cigar_code);
	cal->info = cigar_code_aux;
      } else {
	//array_list_remove_at(array_list_size(cigar_code->ops) - 1, cigar_code->ops);
	for (int c = 0; c < array_list_size(cigar_code_aux->ops); c++) {
	  cigar_op = array_list_get(c, cigar_code_aux->ops);
	  array_list_insert(cigar_op, cigar_code->ops);
	}
	cigar_code->distance += cigar_code_aux->distance;
	cigar_code_free(cigar_code_aux);
	cal->info = cigar_code;
	//printf("\tLast Final Distance %i\n", cigar_code->distance);
      }
      s->info = cal->info;
    }// else { printf("\t NULL\n"); }

    fusion_coords_free(extrem_coords[i]);

  }

  for (i = 0; i < num_sw_sp; i++) {
    cal = sp_coords[i]->cal_ref;
    cigar_code_aux = sp_cigars[i];
    //printf("Merge Cigars SP: %i/%i [%i:%lu-%lu] %s\n", i, num_sw_sp, cal->chromosome_id, cal->start, cal->end, sp_coords[i]->id);
    if (cigar_code_aux) {
      //printf("Pass Score with %f: [%i:%lu-%lu]", score, cal->chromosome_id, cal->start, cal->end);
      s = (seed_region_t *) linked_list_get_first(cal->sr_list);
      cigar_code = (cigar_code_t *)s->info;
      
      for (int c = 0; c < array_list_size(cigar_code_aux->ops); c++) {
	cigar_op = array_list_get(c, cigar_code_aux->ops);
	array_list_insert(cigar_op, cigar_code->ops);
      }
      cigar_code_free(cigar_code_aux);
    } //else { printf("NULL\n"); }

    fusion_coords_free(sp_coords[i]);
  }


  float reads_score[100];
  int n_process;
  int n_best = 0;

  min_score = 60.0;

  //printf("Report results\n");

  //printf("Last Section\n");
  for (t = 0; t < num_targets; t++) {
    target = mapping_batch->targets[t];
    fq_read = array_list_get(target, mapping_batch->fq_batch);
    strcpy(query_revcomp, fq_read->sequence);
    seq_reverse_complementary(query_revcomp, fq_read->length);
    mapping_batch->mapping_lists[target]->size = 0;

    //printf("RNA READ FINAL SCORE: %s\n", fq_read->id);
    n_best = 0;
    //==== New CAL order by new Score ====//
    for (j = 0; j < array_list_size(cals_targets[target]); j++) {
      fusion_cals = array_list_get(j, cals_targets[target]);
      norm_score = 0.0;
      
      for (int z = 0; z < array_list_size(fusion_cals); z++) {
	// get cal and read index
	cal = array_list_get(z, fusion_cals);
	//printf("CAL (%i)[%i:%lu-%lu] %i:\n", cal->strand, cal->chromosome_id, cal->start, cal->end, j);
	if (cal->sr_list->size == 0) continue;	
	s = (seed_region_t *) linked_list_get_first(cal->sr_list);
	cigar_code = (cigar_code_t *) s->info;
	norm_score += cigar_code_get_score(fq_read->length, cigar_code);

	for (int o = 0; o < array_list_size(cigar_code->ops); o++) {
	  cigar_op = array_list_get(o, cigar_code->ops);
	  if (cigar_op->name == 'D' && cigar_op->number >= min_intron_size) {
	    //printf("EXTRA DELETION[ %i:%lu-%lu ] %s, %s\n", cal->chromosome_id, cal->start, cal->end, 
	    //	   new_cigar_code_string(cal->info), fq_read->id);
	    norm_score += ((cigar_op->number*100)/fq_read->length);
	  }
	}
	//printf("\t Distance %i: %s\n", cigar_code->distance, new_cigar_code_string(cigar_code));
      }
      //if (j == 0) { printf("Max score %f...", norm_score);norm_score += 20.0; }
      reads_score[j] = norm_score;
      if (norm_score > 80.0) { n_best++; }
      //printf("\tCALs Score %f\n", reads_score[j]);
    }

    //====================================//


    num_cals = array_list_size(cals_targets[target]);
    //if (n_alignments > num_cals) {
      n_process = num_cals;
      //} else {
      //if (n_best > n_alignments) { n_process = n_best; }
      //else { n_process = n_alignments; }
      //}

    //n_process = array_list_size(cals_targets[target]);
    //===== Order by score =====//
    for (i = 0; i < n_process; i++) {
      for (j = i + 1; j < array_list_size(cals_targets[target]); j++) {
	if (reads_score[j] > reads_score[i]) {
	  array_list_swap(i, j, cals_targets[target]);
	  norm_score = reads_score[j];
	  reads_score[j] = reads_score[i];
	  reads_score[i] = norm_score;
	}
      }
    }
    //===== Step-3: END =====//
      
    //printf("Norm Score = %f > min_score = %f\n", norm_score, min_score);
    //==== Report the n best alignments ====//
    for (i = 0; i < n_process; i++) {
      if (reads_score[i] >= min_score) {
	fusion_cals = array_list_get(i, cals_targets[target]);
	cal = array_list_get(0, fusion_cals);      
	alignment = alignment_new();
	sprintf(cigar_str, "%iM", strlen(fq_read->sequence));
	alignment_init_single_end(strdup(&fq_read->id[1]), strdup(fq_read->sequence), strdup(fq_read->quality),
				  cal->strand,
				  cal->chromosome_id - 1, cal->start - 1,
				  strdup(cigar_str), 1, 254, 1, 1,
				  0, NULL, 0, alignment);
	array_list_insert(alignment, mapping_batch->mapping_lists[target]);     
      }
    }
    
    if (array_list_size(mapping_batch->mapping_lists[target]) == 0) {
      fusion_cals = array_list_get(0, cals_targets[target]);
      cal = array_list_get(0, fusion_cals);      
      
      if (fq_read->id[1] != 'r' && fq_read->id[2] != 'a') {
	//Extra
	//exit(-1);
      }
    } else { printf("Yes!\n"); }
    //Is alignment in a good region? Delete
    
    //==== ==== ==== ==== ==== ==== ==== ====//
    for (i = 0; i < array_list_size(cals_targets[target]); i++) {
      fusion_cals = array_list_get(i, cals_targets[target]);
      for (int z = 0; z < array_list_size(fusion_cals); z++) {
	cal = array_list_get(z, fusion_cals);
	s = (seed_region_t *) linked_list_get_first(cal->sr_list);
	cigar_code = (cigar_code_t *) s->info;
	array_list_clear(cigar_code->ops, cigar_op_free);
	cigar_code_free(cigar_code);
      }
      array_list_free(fusion_cals, cal_free);
    }
    
    //printf("Last Section End\n");

    array_list_free(cals_targets[target], NULL);
  }
  
  free(new_targets);
  //free(cals_score);
  free(cals_targets);


  //apply_sw(input_p, batch);

  return RNA_POST_PAIR_STAGE;
  }*/


      /*} else {
      //fusion_cals = array_list_get(0, cals_targets[target]);
      //cal = array_list_get(0, fusion_cals);      
      
      //Is alignment in a good region?
      if (fq_read->id[1] != 'r' && fq_read->id[2] != 'a') {
	int found;
	int pos = 0, n_arroba = 0;
	char *seq_id = fq_read->id;
	char num[124];
	int p_num = 0;
	int chromosome_ok;
	size_t start_ok, end_ok;
	while (n_arroba < 4) { if (seq_id[pos++] == '@') { n_arroba++; } }
	p_num = 0;
	
	while (seq_id[pos] != '@') { num[p_num++] = seq_id[pos++]; }
	num[p_num] = '\0';
	if (num[0] == 'X') { chromosome_ok = 23; }
	else if (num[0] == 'Y') { chromosome_ok = 24; }
	else { chromosome_ok = atoi(num); }
	
	p_num = 0;pos++;
	while (seq_id[pos] != '@') { num[p_num++] = seq_id[pos++]; }
	num[p_num] = '\0';
	start_ok = atof(num);
	
	p_num = 0;pos++;
	while (seq_id[pos] != '@') { num[p_num++] = seq_id[pos++]; }
	num[p_num] = '\0';
	end_ok = atof(num);
	
	found = 0;
	for (int a = 0; a < array_list_size(mapping_batch->mapping_lists[target]); a++) {
	  alignment = array_list_get(a, mapping_batch->mapping_lists[target]);
	  if (alignment->chromosome  + 1 == chromosome_ok && 
	      alignment->position >= start_ok &&
	      alignment->position <= end_ok) {
	    found = 1;
	  }
	}
	
	if (!found) {
	  printf(" :::::?%s\n", fq_read->id);
	  printf(" :::::?%s\n", fq_read->sequence);
	  printf(" :::::?+\n", fq_read->id);
	  printf(" :::::?%s\n", fq_read->quality);
	  for (int a = 0; a < array_list_size(mapping_batch->mapping_lists[target]); a++) {
	    alignment = array_list_get(a, mapping_batch->mapping_lists[target]);
	    printf(" ::::: [%i:%lu]\n", alignment->chromosome, alignment->position);
	  }
	}      
      //Extra
      }
    }
    */


    /* 
    for (int i = 0; i < array_list_size(cals_targets[target]); i++) {
      //Only save the seeds in range for the first CAL, other CALs need seeds??? or is the first CAL the correct
      fusion_cals = array_list_get(i, cals_targets[target]);    
      if (scores_ranking[target][i] <= 70.0) {
	first_cal = array_list_get(0, fusion_cals);
	last_cal = array_list_get(array_list_size(fusion_cals) - 1, fusion_cals);

	s_prev = linked_list_get_first(first_cal->sr_list);
	s_next = linked_list_get_last(last_cal->sr_list);
	
	if (s_prev->read_start >= seed_err_size) {
	  printf("\t @@@@SEEDs in First positions [%i-%i](%i)\n", 0, s_prev->read_start - 1, first_cal->strand);
	  if (first_cal->strand == 0) {
	    start_seeding = 0;
	    end_seeding = s_prev->read_start - 1;
	  } else {
	    start_seeding = fq_read->length - s_prev->read_start;
	    end_seeding = fq_read->length - 1;
	  }
	  
	  lim_start = first_cal->start - 500000;
	  lim_end = first_cal->start;
	  
	  sequence = fq_read->sequence;
	  //Seeds in Start position
	  array_list_t *mapping_list_prev = array_list_new(1000,
							   1.25f,
							   COLLECTION_MODE_ASYNCHRONIZED);
	  
	  bwt_map_inexact_seeds_by_region(start_seeding, end_seeding,
					  first_cal->strand,
					  first_cal->chromosome_id,
					  lim_start,
					  lim_end,
					  sequence, seed_err_size,
					  seed_err_size,
					  bwt_optarg,
					  bwt_index,
					  mapping_list_prev);	
	  //First Merge Seeds
	  fusion_regions(mapping_list_prev, fq_read->length);	
	  //if (array_list_size(mapping_list_prev) > 1) { array_list_clear(mapping_list_prev, region_bwt_free); }
	  region_t *region, *region_prev = NULL;
	  for (int r = 0; r < array_list_size(mapping_list_prev); r++) {
	    region = array_list_get(r, mapping_list_prev);
	    if (region_prev == NULL || region_prev->id != region->id) { 
	      printf("\t HELP CAL [%i:%lu|%i-%i|%lu]\n", region->chromosome_id, region->start, 
	      	   region->seq_start, region->seq_end, 
	      	   region->end);
	      linked_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	      seed_region = seed_region_new(region->seq_start, region->seq_end,
					    region->start, region->end, region->id);
	      linked_list_insert_last(seed_region, linked_list);
	      cal = cal_new(first_cal->chromosome_id, first_cal->strand, 
			    region->start, region->end, 
			    1, linked_list, 
			    linked_list_new(COLLECTION_MODE_ASYNCHRONIZED));
	      //if (cal->strand == 0) {
	      array_list_insert_at(r, cal, fusion_cals);
	      array_list_insert_at(r, cal, cals_list);
	      //} else {
	      //array_list_insert(cal, fusion_cals);
	      //array_list_insert(cal, cals_list);
	      //}
	    } else {
	      printf("ERROR SEEDS 1 ERROR: %s\n", fq_read->id);
	      assert(region_prev->id);
	    }
	    region_prev = region;
	  }
	  array_list_free(mapping_list_prev, region_bwt_free);	
	}      
	
	if (s_next->read_end <= fq_read->length - seed_err_size) {
	  if (first_cal->strand == 0) {
	    start_seeding = s_next->read_end;
	    end_seeding = fq_read->length - 1;
	  } else {
	    start_seeding = 0;
	    end_seeding = fq_read->length - s_next->read_end - 1;
	  }
	  
	  lim_start = last_cal->end;
	  lim_end = last_cal->end + 500000;
	  
	  sequence = fq_read->sequence;
	  
	  printf("\t @@@@SEEDs in Last positions [%i-%i](%i)\n", s_next->read_end, fq_read->length - 1, last_cal->strand);
	  //Seeds in End position
	  array_list_t *mapping_list_next = array_list_new(1000,
							   1.25f,
							   COLLECTION_MODE_ASYNCHRONIZED);
	  
	  bwt_map_inexact_seeds_by_region(start_seeding, end_seeding, 
					  last_cal->strand,
					  last_cal->chromosome_id,
					  lim_start,
					  lim_end,
					  sequence, seed_err_size,
					  seed_err_size/2,
					  bwt_optarg,
					  bwt_index,
					  mapping_list_next);
	  
	  //First Merge Seeds
	  fusion_regions(mapping_list_next, fq_read->length);	
	  //if (array_list_size(mapping_list_next) > 1) { array_list_clear(mapping_list_next, region_bwt_free); }
	  
	  region_t *region, *region_prev = NULL;
	  for (int r = 0; r < array_list_size(mapping_list_next); r++) {
	    region = array_list_get(r, mapping_list_next);
	    if (region_prev == NULL || region_prev->id != region->id) { 
	      printf("\t HELP CAL [%i:%lu|%i-%i|%lu]\n", region->chromosome_id, region->start, 
	      	   region->seq_start, region->seq_end, region->end);
	      linked_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	      seed_region = seed_region_new(region->seq_start, region->seq_end,
					    region->start, region->end, region->id);
	      linked_list_insert_last(seed_region, linked_list);
	      cal = cal_new(first_cal->chromosome_id, first_cal->strand, 
			    region->start, region->end, 
			    1, linked_list, 
			    linked_list_new(COLLECTION_MODE_ASYNCHRONIZED));
	      //if (cal->strand == 0) {
	      array_list_insert(cal, fusion_cals);
	      array_list_insert(cal, cals_list);
	      //} else {
	      //array_list_insert_at(r, cal, fusion_cals);
	      //array_list_insert_at(r, cal, cals_list);	      
	      //}
	    } else {
	      printf("ERROR SEEDS 1 ERROR: %s\n", fq_read->id);
	      assert(region_prev->id);
	    }
	    region_prev = region;
	  }
	  array_list_free(mapping_list_next, region_bwt_free);	  
	}

	printf("================ FUSION CALS MERGE     ===================\n");
	for (int p = 0; p < array_list_size(fusion_cals); p++) {
	  cal = array_list_get(p, fusion_cals);
	  s_prev = linked_list_get_first(cal->sr_list);
	  s_next = linked_list_get_last(cal->sr_list);
	  printf("\t(%i)CAL FUSION(%i):[%i:%lu|%i-%i|%lu]\n", cal->strand, p, cal->chromosome_id, cal->start, 
		 s_prev->read_start, s_next->read_end, cal->end);
	}
	printf("================ FUSION CALS MERGE END ===================\n");
      }
    }
    */


      /*

  /*
    //error = 0;//delete
    if (error) {
      printf("%s\n", fq_read->id);
      printf("Merge CALs %i\n", i);
      for (j = 0; j < array_list_size(fusion_cals); j++) {
	cal = array_list_get(j, fusion_cals);
	printf("\t\t CAL(%i) [%i:%lu-%lu]:\n", j, cal->chromosome_id, cal->start,
	       cal->end);
	
	linked_list_iterator_init(cal->sr_list, &itr);
	s = (seed_region_t *) linked_list_iterator_curr(&itr);
	  while (s != NULL) {
	    printf("\t\t\t seed: [%i|%i - %i|%i]\n",
		   s->genome_start, s->read_start, s->read_end, s->genome_end);
	    s = (seed_region_t *) linked_list_iterator_next(&itr);
	  }
      }
      exit(-1);
    }*/
    //=====================================================================================
      /*for (int c = 0; c < array_list_size(fusion_cals); c++) {
	cal = array_list_get(c, fusion_cals);
	cigar_code = (cigar_code_t *)cal->info;
	printf("CIGAR: ");
	for (int c_o = 0; c_o < array_list_size(cigar_code->ops); c_o++) {
	  cigar_op = array_list_get(c_o, cigar_code->ops);
	  printf("%i%c", cigar_op->number, cigar_op->name);
	  if (cigar_op->name == 'D' && cigar_op->number > 40) { 
	    //printf("%s: %i:%lu-%lu\n", fq_read->id, cal->chromosome_id, cal->start, cal->end);
	    //exit(-1);
	  }
	}
	printf("\n");
	}*/


  /*  
  LOG_DEBUG_F("%s********************RNA PHASE*******************%s\n", KGRN, KWHT);

  for (t = 0; t < num_targets; t++) {
    target = mapping_batch->targets[t];
    cals_list = mapping_batch->mapping_lists[target];
    fq_read = array_list_get(target, mapping_batch->fq_batch);    
    read_length = fq_read->length;
    strcpy(query_revcomp, fq_read->sequence);
    seq_reverse_complementary(query_revcomp, fq_read->length);

    LOG_DEBUG_F("Process Read %s\n", fq_read->id);

    //===== Step-0: Concatenate CALs =====//
    /*for (i = array_list_size(cals_list) - 1; i >= 0; i--) {
      cal = array_list_get(i, cals_list);
      cigar_code = cal->info;
      if (!cigar_code) { array_list_remove_at(i, cals_list); }
      }
    //===== End Step-0 =====//

    num_cals = array_list_size(cals_list);

    //===== Step-1: Concatenate CALs =====//
    LOG_DEBUG("STEP-1");
    cals_targets[target] = array_list_new(num_cals,
					  1.25f,
					  COLLECTION_MODE_ASYNCHRONIZED);
    merge_cals = array_list_new(10,
				1.25f,
				COLLECTION_MODE_ASYNCHRONIZED);
    cal_pos = 0;
    cal_prev = (cal_t *)array_list_get(cal_pos++, cals_list);
    LOG_DEBUG_F("\t CAL (%i)[%i:%lu-%lu]: %s\n", cal_prev->strand, cal_prev->chromosome_id, cal_prev->start,
		cal_prev->end, new_cigar_code_string(cal_prev->info));
    array_list_insert(cal_prev, merge_cals);    
    while (cal_pos < num_cals) {
      cal_next = (cal_t *)array_list_get(cal_pos, cals_list);      
      LOG_DEBUG_F("\t CAL (%i)[%i:%lu-%lu]: %s\n", cal_next->strand, cal_next->chromosome_id, cal_next->start,
		  cal_next->end, new_cigar_code_string(cal_next->info));

      cigar_code_prev = cal_prev->info;
      cigar_op_prev = cigar_code_get_last_op(cigar_code_prev);

      cigar_code = cal->info;
      cigar_op = cigar_code_get_first_op(cigar_code);

      if (cal_prev->chromosome_id == cal_next->chromosome_id && 
	  cal_prev->strand == cal_next->strand && 
	  (cal_next->start <= (cal_prev->end + max_intron_size)) &&
	  cigar_op_prev->name == 'H' && cigar_op->name == 'H') {
	LOG_DEBUG("\t Fusion!\n");
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
    
    
    //===== Step-2: Generate CALs Score and filter reads that were process =====//
    LOG_DEBUG("STEP-2");
    for (i = 0; i < array_list_size(cals_targets[target]); i++) {
      fusion_cals = array_list_get(i, cals_targets[target]);
      for (j = 0; j < array_list_size(fusion_cals); j++) {	
	cals_score[i] = generate_cals_score(fusion_cals, fq_read->length);
      }
    }
    //===== Step-2: END =====//

    number_of_best = array_list_size(cals_targets[target]);
        
    //===== Step-3: Ranking CALs by score (the n best)=====//
    LOG_DEBUG("STEP-3");
    for (i = 0; i < number_of_best; i++) {
      for (j = i + 1; j < array_list_size(cals_targets[target]); j++) {
	if (cals_score[j] > cals_score[i]) {
	  array_list_swap(i, j, cals_targets[target]);
	  coverage = cals_score[j];
	  cals_score[j] = cals_score[i];
	  cals_score[i] = coverage;
	}
      }
    }
    //===== Step-3: END =====//


    //===== Step-4: Remove the complete best scores from the list  =====//
    LOG_DEBUG("STEP-4");
    n_delete = 0;
    delete_not_cigar = 0;
    for (i = 0; i < number_of_best; i++) {
      if (cals_score[i] == 0) { delete_not_cigar = 1; LOG_DEBUG("Score 0! Alert!\n"); break; }
      if (cals_score[i] > fq_read->length && 
	  array_list_size(cals_targets[target]) == 1) {
	//CAL Mapped
	printf("%i - %i\n", array_list_size(cals_targets[target]), target);
	fusion_cals = array_list_remove_at(i, cals_targets[target]);
	cal = array_list_get(0, fusion_cals);
	cigar_code = (cigar_code_t *)cal->info;
	LOG_DEBUG_F("+++++ Read Mapped %s +++++\n", new_cigar_code_string(cigar_code));
	n_delete++;

	cigar_op = cigar_code_get_first_op(cigar_code);
	if (cigar_op->name == 'H') { 
	  padding_left =  cigar_op->number;
	} else {
	  padding_left = 0;
	}

	cigar_op = cigar_code_get_last_op(cigar_code);
	if (cigar_op->name == 'H') { 
	  padding_right = cigar_op->number;
	} else {
	  padding_right = 0;
	}

	len_query = strlen(fq_read->sequence) - padding_left - padding_right;
	query_map = (char *)calloc(len_query + 1, sizeof(char));
	strncpy(query_map, fq_read->sequence + padding_left, len_query);
	quality_map = (char *)calloc(len_query + 1, sizeof(char));

	if (cal->strand == 1) {
	  seq_reverse_complementary(query_map, len_query);
	  reverse_str(fq_read->quality + padding_left, quality_map, len_query);
	} else {
	  strncpy(quality_map, fq_read->quality + padding_left, len_query);
	}

	sprintf(cigar_str, "%iM", strlen(fq_read->sequence));

	alignment = alignment_new();
	alignment_init_single_end(strdup(fq_read->id), strdup(fq_read->sequence), strdup(fq_read->quality),
				  cal->strand,
				  cal->chromosome_id - 1, cal->start - 1,
				  strdup(cigar_str), 1, 254, 1, 1,
				  0, NULL, 0, alignment);
	//alignment_print(alignment);
	
	//TODO:For debugug
	//alignment_print(alignment);
	if (cigar_read_coverage(cigar_code) != strlen(alignment->sequence)) { assert(alignment->cigar); }

	array_list_insert(alignment, mapping_batch->mapping_lists[target]);
	cal_free(cal);
      }
    }

    if (delete_not_cigar) { continue; }

    if (n_delete != number_of_best && array_list_size(cals_targets[target]) ) {
      new_targets[num_new_targets++] = target;
    } else {
      continue;
    }
    //===== Step-4: END =====//
    

    //===== Step-5: Extend to first mismatch and Generate CALs score  =====//
    LOG_DEBUG("STEP-5");
    for (i = 0; i < number_of_best; i++) {
      fusion_cals = array_list_get(i, cals_targets[target]);
      //LOG_DEBUG_F("CAL merge %i\n", i);
      for (j = 0; j < array_list_size(fusion_cals); j++) {
	//LOG_DEBUG_F("\tCAL id %i\n", j);
	cal = array_list_get(j, fusion_cals);

	if (cal->strand == 1) {
	  query_map = query_revcomp;
	} else {
	  query_map = fq_read->sequence;
	}

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
	    memcpy(query, query_map, cigar_op_start->number);
	    query[cigar_op_start->number] = '\0';

	    //LOG_DEBUG_F("\t\tExtract reference for Start hard clipping [%i:%lu-%lu]: \n", cal->chromosome_id, genome_start, genome_end);
	    //LOG_DEBUG_F("REFERENCE: %s\n", reference);
	    //LOG_DEBUG_F("QUERY    : %s\n", query);
	    
	    c = cigar_op_start->number - 1;
	    exact_nt = 0;
	    while (reference[c] == query[c]) { exact_nt++; c--; }
	    cigar_op = cigar_code_get_op(1, cigar_code);
	    cigar_op->number += exact_nt;
	    cigar_op_start->number -= exact_nt;
	  }
		
	  if (cigar_op_end->name == 'H') {
	    genome_end = s->genome_end;
	    genome_start = s->genome_end - cigar_op_end->number + 1;
	    genome_read_sequence_by_chr_index(reference, 0, 
	    				      cal->chromosome_id - 1, &genome_start, &genome_end, genome);
	    memcpy(query, (query_map + fq_read->length) - cigar_op_end->number , cigar_op_end->number);
	    query[cigar_op_end->number] = '\0';

	    LOG_DEBUG_F("\t\tExtract reference for End hard clipping [%i:%lu-%lu]: \n", cal->chromosome_id, genome_start, genome_end);
	    LOG_DEBUG_F("REFERENCE: %s\n", reference);
	    LOG_DEBUG_F("QUERY    : %s\n", query);

	    c = 0;
	    exact_nt = 0;
	    while (reference[c] == query[c]) { exact_nt++; c++; }
	    cigar_op = cigar_code_get_op(cigar_code_get_num_ops(cigar_code) - 2, cigar_code);
	    cigar_op->number += exact_nt;
	    cigar_op_end->number -= exact_nt;
	  }
	  //LOG_DEBUG_F("\t\tNew Cigar %s\n", new_cigar_code_string(cigar_code));
	}
      }

      cals_score[i] = generate_cals_score(fusion_cals, fq_read->length);

    }
    //===== Step-5: END =====//


    //===== Step-6: Process CALs and mount Smith-Watermans =====//
    LOG_DEBUG("STEP-6");
    for (i = 0; i < number_of_best; i++) {
      fusion_cals = array_list_get(i, cals_targets[target]);
      
      j = 0;
      cal = array_list_get(j, fusion_cals);
      if (cal->strand == 1) {
	query_map = query_revcomp;
      } else {
	query_map = fq_read->sequence;
      }

      cigar_code = (cigar_code_t *)cal->info;
      if (!cigar_code) { assert(cigar_code); }

      s = linked_list_get_first(cal->sr_list);
      cigar_op = cigar_code_get_first_op(cigar_code);
      LOG_DEBUG_F("\tCAL CIGAR %i: %s\n", j, new_cigar_code_string(cigar_code));
      
      if (cigar_op->name == 'H' && cigar_op->number > 20) {
	flank = 10;
	printf("%i\n", cigar_op->number);
	genome_start = s->genome_start;
	genome_end = s->genome_start + cigar_op->number + flank - 1;
	flank = 10;
	genome_read_sequence_by_chr_index(reference, 0, 
					  cal->chromosome_id - 1, &genome_start, &genome_end, genome);
	memcpy(query, query_map, cigar_op->number + flank);
	query[cigar_op->number + flank] = '\0';

	r[num_sw] = strdup(reference);
	q[num_sw] = strdup(query);
	fusion_coords[num_sw++] = fusion_coords_new(0, 0, 0, 0, 0, 0, 0, 0, SW_LEFT, NULL, NULL);

	LOG_DEBUG_F("START-REFERENCE: %s\n", reference);
	LOG_DEBUG_F("START-QUERY    : %s\n", query);	
	  
	cal->l_flank = flank;
      }

      j++;
      cal_prev = cal;
      flank = 20;
      for (; j < array_list_size(fusion_cals); j++) {
	read_start = 0;
	s = linked_list_get_first(cal_prev->sr_list);
	cigar_code = (cigar_code_t *)cal_prev->info;
	cigar_op = cigar_code_get_last_op(cigar_code);
	if (cigar_op->name != 'H') { assert(cigar_op); }

	cigar_op = cigar_code_get_first_op(cigar_code);
	genome_start = s->genome_start;
	if (cigar_op->name == 'H') { 
	  genome_start += cigar_op->number; 
	  read_start += cigar_op->number;	  
	}
	
	read_start += cigar_read_coverage(cigar_code);

	if (cigar_read_coverage(cigar_code) < flank ) {
	  assert(flank);
	}

	read_start -= flank;
	cal_prev->r_flank = flank;
	genome_start += cigar_genome_coverage(cigar_code);
	genome_end = genome_start + flank;
	genome_start -= flank;
	genome_read_sequence_by_chr_index(reference, 0,
					  cal->chromosome_id - 1, &genome_start, &genome_end, genome);

	LOG_DEBUG_F("MIDDLE-REFERENCE [%i:%lu-%lu]: %s\n", cal->chromosome_id, genome_start, genome_end, reference);
	
	cal = array_list_get(j, fusion_cals);
	if (cal->strand == 1) {
	  query_map = query_revcomp;
	} else {
	  query_map = fq_read->sequence;
	}
	
	s = linked_list_get_first(cal->sr_list);
	//flank = 20;
	cigar_code = (cigar_code_t *)cal->info;
	cigar_op = cigar_code_get_first_op(cigar_code);
	LOG_DEBUG_F("\tCAL CIGAR %i: %s\n", j, new_cigar_code_string(cigar_code));

	genome_start2 = s->genome_start;
	if (cigar_op->name == 'H') { genome_start2 += cigar_op->number; }
	else { assert(cigar_op->name); }

	if (cigar_read_coverage(cigar_code) < flank ) {
	  assert(flank);
	}
	read_end = cigar_op->number;
	read_end += flank;
	cal->l_flank = flank;
	genome_end2 = genome_start2 + flank;
	genome_start2 -= flank;
	genome_read_sequence_by_chr_index(reference_aux, 0,
					  cal->chromosome_id - 1, &genome_start2, &genome_end2, genome);
	LOG_DEBUG_F("MIDDLE-REFERENCE [%i:%lu-%lu]: %s\n", cal->chromosome_id, genome_start2, genome_end2, reference_aux);

	strcat(reference, reference_aux);	
	memcpy(query, query_map + read_start, read_end - read_start); 
	query[read_end - read_start] = '\0';
	LOG_DEBUG_F("MIDDLE-QUERY     [%lu-%lu]: %s\n", read_start, read_end, query);	

	r[num_sw] = strdup(reference);
	q[num_sw] = strdup(query);
	fusion_coords[num_sw++] = fusion_coords_new(genome_start, genome_end, 
						    genome_start2, genome_end2, 
						    read_start, read_end,
						    cal->chromosome_id, cal->strand,
						    SW_MIDDLE, NULL, NULL);
	cal_prev = cal;
      }

      cigar_op = cigar_code_get_last_op(cigar_code);

      if (cigar_op->name == 'H' && cigar_op->number > 20) {
	genome_end = s->genome_end + flank;
	flank = 10;
	genome_start = s->genome_end - cigar_op->number - flank + 1;
	genome_read_sequence_by_chr_index(reference, 0, 
					  cal->chromosome_id - 1, &genome_start, &genome_end, genome);
	memcpy(query, (query_map + fq_read->length) - (cigar_op->number + flank), cigar_op->number + flank);
	query[cigar_op->number + flank] = '\0';
	cal->r_flank = flank;
	r[num_sw] = strdup(reference);
	q[num_sw] = strdup(query);
	fusion_coords[num_sw++] = fusion_coords_new(0, 0, 0, 0, 0, 0, 0, 0, SW_RIGHT, NULL, NULL);
	LOG_DEBUG_F("END-REFERENCE: %s\n", reference);
	LOG_DEBUG_F("END-QUERY    : %s\n", query);
      }
    }
    //===== Step-6: Process CALs =====//
  } //End loop targets 

  LOG_DEBUG("========== S M I T H - W A T E R M A N ==========\n");

  LOG_DEBUG("I N P U T \n");
  for (i = 0;  i < num_sw; i++) {
    LOG_DEBUG_F("Query (%i): %s\n", strlen(q[i]), q[i]);
    LOG_DEBUG_F("Refer (%i): %s\n", strlen(r[i]), r[i]);
  }
  
  output = sw_multi_output_new(num_sw);
  smith_waterman_mqmr(q, r, num_sw, sw_optarg, 1, output);
  //cigar_code_t *cigar_codes_list[num_sw];

  LOG_DEBUG("O U T P U T \n");
  for (i = 0;  i < num_sw; i++) {
    LOG_DEBUG_F("Q:%s (Start:%i, Len:%i)\n", output->query_map_p[i], output->query_start_p[i], strlen(output->query_map_p[i]));
    LOG_DEBUG_F("R:%s (Start:%i, Len:%i)\n", output->ref_map_p[i],   output->ref_start_p[i],   strlen(output->ref_map_p[i]));
      
    if (fusion_coords[i]->type_sw == SW_LEFT ||
	fusion_coords[i]->type_sw == SW_RIGHT) {
      score = NORM_SCORE(output->score_p[i], strlen(q[i]), sw_optarg->subst_matrix['A']['A']);
    } else {      
      score = 1.0;
    }

    if (score > 0.7) {
      cigar_codes_list[i] = genoerate_cigar_sw_output(output->query_map_p[i],
						      output->ref_map_p[i],
						      fusion_coords[i]->l_exon_start,
						      fusion_coords[i]->l_exon_end,
						      fusion_coords[i]->r_exon_start,
						      fusion_coords[i]->r_exon_end,
						      fusion_coords[i]->read_start,
						      fusion_coords[i]->read_end,
						      fusion_coords[i]->chromosome,
						      fusion_coords[i]->strand,
						      output->query_start_p[i],
						      output->ref_start_p[i],
						      fusion_coords[i]->type_sw, 
						      strlen(q[i]),
						      chromosome_avls, NULL);
      LOG_DEBUG_F("Generate CIGAR from SW output str %f: %s\n", score, new_cigar_code_string(cigar_codes_list[i]));
    } else {
      cigar_codes_list[i] = NULL;
      LOG_DEBUG_F("Generate CIGAR from SW output str: %s, score=%f\n", "Not CIGAR", score);    
    }
  }

  LOG_DEBUG("========== = = = = = = = = = = = = = = ==========\n");


  sw_pos = 0;
  //===== Last Phase: Process all results and report alignments =====//
  printf("Num new targets %i\n", num_new_targets);
  for (t = 0; t < num_new_targets; t++) {
    target = new_targets[t];
    cals_list = cals_targets[target];
    fq_read = array_list_get(target, mapping_batch->fq_batch);
    read_length = fq_read->length;
    num_cals = array_list_size(cals_list);
    LOG_DEBUG_F("CAL TARGET %i: , %s\n", target, fq_read->id);

    //number_of_best = num_cals;
    number_of_best = array_list_size(cals_targets[target]);
    
    for (i = 0; i < number_of_best; i++) {
      fusion_cals = array_list_get(i, cals_targets[target]);
      j = 0;
      alig_cigar_code = cigar_code_new();
      cal = array_list_get(j, fusion_cals);
      first_cal = cal;
      cigar_code_prev = (cigar_code_t *)cal->info;
      cigar_op = cigar_code_get_first_op(cigar_code_prev);

      int z = 0;

      //===== Step-1: Merge starts cigars and Select SW-Cigar or CAL-Cigar =====//
      if (cigar_op->name == 'H' && cigar_op->number > 20) {
	cigar_code = cigar_codes_list[sw_pos++];
	if (cigar_code) {
	  cigar_code_merge(alig_cigar_code, cigar_code);
	  z = 1;	
	  //===== Refresh value if we use cal->flank =====//
	  if (cal->l_flank) {
	    cigar_op = array_list_get(z, cigar_code_prev->ops);
	    cigar_op->number -= cal->l_flank;
	    cigar_op_aux = array_list_remove_at(array_list_size(cigar_code_prev->ops) - 1, cigar_code_prev);
	    cigar_op->number += cigar_op->number;
	  }
	}
      }

      for (; z < array_list_size(cigar_code_prev->ops); z++) {
	cigar_op = array_list_get(z, cigar_code_prev->ops);
	array_list_insert(cigar_op, alig_cigar_code->ops);
      }
      //===== Step-1: End =====//

      //===== Step-2: Merge SpliceJunctions-Cigar =====//
      j++;
      cal_prev = cal;
      for (; j < array_list_size(fusion_cals); j++) {      
	cal = array_list_get(j, fusion_cals);
	cigar_code = (cigar_code_t *)cal->info;
	cigar_code_aux = cigar_codes_list[sw_pos++];	

	//===== Extract the last cigar opration ('H') =====//
	cigar_op = array_list_remove_at(array_list_size(alig_cigar_code->ops) - 1, alig_cigar_code->ops);
	if (cigar_op->name != 'H') { assert(cigar_op); }

	//===== Get the new last cigar operation ('M') to refresh value =====//
	cigar_op = array_list_get(array_list_size(alig_cigar_code->ops) - 1, alig_cigar_code->ops);
	if (cigar_op->name != 'M') { assert(cigar_op); }
	cigar_op->number -= cal_prev->r_flank;
	
	//===== Get first operation SW-Cigar And Merge operations with final CIGAR =====//
	cigar_op_aux = cigar_code_get_first_op(cigar_code_aux);
	if (cigar_op_aux->name != 'M') { assert(cigar_op_aux->name ); }
	cigar_op->number += cigar_op_aux->number;
	//-->Merge Splice<--
	for (int z = 1; z < array_list_size(cigar_code_aux->ops) - 1; z++) {
	  cigar_op_aux = array_list_get(z, cigar_code_aux->ops);
	  array_list_insert(cigar_op_aux, alig_cigar_code->ops);
	}
	
	//===== Merge last Exon. First, get the last operation SW-Cigar =====//
	cigar_op_aux = cigar_code_get_last_op(cigar_code_aux);
	if (cigar_op_aux->name != 'M' || array_list_size(cigar_code->ops) < 2) { assert(cigar_op_aux->name); }
        
	//===== Second, get  second get second CAL-Cigar operation ('M'), First ('H') and refresh value =====//
	cigar_op = cigar_code_get_op(1, cigar_code);
	cigar_op->number -= cal->l_flank;
	cigar_op->number += cigar_op_aux->number;

	//===== Merge Cigars =====//
	for (int z = 1; z < array_list_size(cigar_code->ops); z++) {
	  cigar_op = array_list_get(z, cigar_code->ops);
	  array_list_insert(cigar_op, alig_cigar_code->ops);
	}

	cal_prev = cal;
	//cigar_code_prev = (cigar_code_t *)cal->info;
      }
      //===== Step-2: End =====//

      cigar_op = cigar_code_get_last_op(alig_cigar_code);

      //===== Step-3: Merge ends cigars and Select SW-Cigar or CAL-Cigar =====//	
      if (cigar_op->name == 'H' && cigar_op->number > 20) {
	cigar_code = cigar_codes_list[sw_pos++];
	if (cigar_code) {
	  //===== Extract the last cigar opration ('H') =====//
	  cigar_op = array_list_remove_at(array_list_size(alig_cigar_code->ops) - 1, alig_cigar_code->ops);
	  if (cigar_op->name != 'H') { assert(cigar_op); }

	   //===== Refresh value if we use cal->flank =====//
	  if (cal->r_flank) {
	    cigar_op = cigar_code_get_first_op(cigar_code);
	    cigar_op->number -= cal->r_flank;
	    if (cigar_op->number != 0) {	   
	      cigar_op_aux = array_list_get(array_list_size(alig_cigar_code->ops) - 1, alig_cigar_code->ops);
	      cigar_op_aux->number += cigar_op->number;
	    }
	    array_list_remove_at(0, cigar_code->ops); 
	  }
	  //===== Merge SW-Cigar with Final-Cigar =====//
	  cigar_code_merge(alig_cigar_code, cigar_code);
	  
	}
      }
      //===== Step-3: End =====//
      LOG_DEBUG_F("Generate FINAL CIGAR: %s\n", new_cigar_code_string(alig_cigar_code));

      cigar_op = cigar_code_get_first_op(alig_cigar_code);
      if (cigar_op->name == 'H') { 
	padding_left =  cigar_op->number;
      } else {
	padding_left = 0;
      }

      cigar_op = cigar_code_get_last_op(alig_cigar_code);
      if (cigar_op->name == 'H') { 
	padding_right = cigar_op->number;
      } else {
	padding_right = 0;
      }

      if (cal->strand == 1) {
	query_map = query_revcomp;
      } else {
	query_map = fq_read->sequence;
      }


      //len_query = strlen(fq_read->sequence) - padding_left - padding_right;
      //query_map = (char *)calloc(len_query + 1, sizeof(char));
      //strncpy(query_map, fq_read->sequence + padding_left, len_query);
      //quality_map = (char *)calloc(len_query + 1, sizeof(char));
      //strncpy(quality_map, fq_read->quality + padding_left, len_query);
      sprintf(cigar_str, "%iM", strlen(fq_read->sequence));

      alignment = alignment_new();
      alignment_init_single_end(strdup(fq_read->id), strdup(fq_read->sequence), strdup(fq_read->quality),
				first_cal->strand,
				first_cal->chromosome_id - 1, first_cal->start - 1,
				strdup(cigar_str), 
				1, 254, 1, 1,
				0, NULL, 0, alignment);

      //TODO:For debugug
      //alignment_print(alignment);
      if (cigar_read_coverage(alig_cigar_code) != strlen(alignment->sequence)) { assert(alignment->cigar); }

      array_list_insert(alignment, mapping_batch->mapping_lists[target]);
      //alignment_print(alignment);
      cigar_code_free(alig_cigar_code);

    }
  }

  LOG_DEBUG_F("%s********************RNA PHASE END *******************%s\n", KGRN, KWHT);
  
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
