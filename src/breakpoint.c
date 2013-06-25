#include "breakpoint.h"

//--------------------------------------------------------------------------------------

array_list_t *breakpoint_list = NULL;

//--------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------

cigar_op_t *cigar_op_new(int number, char name) {
  cigar_op_t *p = (cigar_op_t *) malloc(sizeof(cigar_op_t));
  p->number = number;
  p->name = name;
  return p;
}

void cigar_op_free(cigar_op_t *p) {
  if (p) free(p);
}

//--------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------

cigar_code_t *cigar_code_new() {
  cigar_code_t *p = (cigar_code_t *) calloc(1, sizeof(cigar_code_t));

  p->distance = 0;
  p->cigar_str = NULL;
  p->ops = array_list_new(20, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  return p;
}

cigar_code_t *cigar_code_new_by_string(char *cigar_str) {
  cigar_code_t *p = cigar_code_new();

  /*
  int cigar_len = strlen(cigar_str);

  p->cigar_str = strdup(cigar_str);
  p->num_allocated_ops = (cigar_len > 10 ? cigar_len / 2 : cigar_len);
  p->ops = (cigar_op_t*) calloc(p->num_allocated_ops, sizeof(cigar_op_t));

  int cigar_op_counter = 0;
  int c = 0;
  char op;
  for (int j = 0; j < cigar_len; j++) {
    op = cigar_str[j];
    if (op < 58) {
      p->ops[cigar_op_counter].number_str[c++] = op;
    } else {
      p->ops[cigar_op_counter].number_str[c++] = '\0';
      p->ops[cigar_op_counter].number = atoi(p->ops[cigar_op_counter].number_str);
      c = 0;
      
      p->ops[cigar_op_counter].name = op;
      cigar_op_counter++;
    }
  }
  
  p->num_ops = cigar_op_counter;
  */

  return p;
}

//--------------------------------------------------------------------------------------

void cigar_code_free(cigar_code_t* p) {
  if (p) {
    if (p->ops) array_list_free(p->ops, (void *) cigar_op_free);
    if (p->cigar_str) free(p->cigar_str);
    free(p);
  }
}

//--------------------------------------------------------------------------------------

int cigar_code_get_num_ops(cigar_code_t *p) {
  int num = 0;
  if (p && p->ops) { 
    return array_list_size(p->ops);
  }
  return num;
}

//--------------------------------------------------------------------------------------

cigar_op_t *cigar_code_get_first_op(cigar_code_t *p) {
  return array_list_get(0, p->ops);
}

cigar_op_t *cigar_code_get_op(int index, cigar_code_t *p) {
  int num_ops = cigar_code_get_num_ops(p);
  if (num_ops > 0 && index < num_ops) {
    return array_list_get(index, p->ops);
  }
  return NULL;
}

cigar_op_t *cigar_code_get_last_op(cigar_code_t *p) {
  int num_ops = cigar_code_get_num_ops(p);
  //if (num_ops > 0) {
  return array_list_get(num_ops - 1, p->ops);
    //}
  //return NULL;
}

//--------------------------------------------------------------------------------------

void cigar_code_append_op(cigar_op_t *op, cigar_code_t *p) {
  if (p && p->ops) {
    cigar_op_t *last = cigar_code_get_last_op(p);
    if (last && last->name == op->name) {
      last->number += op->number;
    } else {
      array_list_insert(op, p->ops);
    }
  }
}

//--------------------------------------------------------------------------------------

void cigar_code_inc_distance(int distance, cigar_code_t *p) {
  if (p && p->ops) {
    p->distance += distance;
  }
}

//--------------------------------------------------------------------------------------

char *new_cigar_code_string(cigar_code_t *p) {
  if (p->cigar_str) {
    free(p->cigar_str);
  }

  int num_ops = cigar_code_get_num_ops(p);
  if (!num_ops) { return NULL; }

  char *str = malloc(num_ops * 5 * sizeof(char));
  *str = 0;

  cigar_op_t *op;
  for (int i = 0; i < num_ops; i++) {
    op = array_list_get(i, p->ops);
    sprintf(str, "%s%i%c", str, op->number, op->name);
  }

  p->cigar_str = str;
  
  return p->cigar_str;
}

//--------------------------------------------------------------------------------------

int cigar_read_coverage(cigar_code_t *p) {
  int coverage = 0;

  if (p) {    
    size_t num_ops = array_list_size(p->ops);
    for (size_t i = 0; i < num_ops; i++) {
      cigar_op_t *cigar_op = array_list_get(i, p->ops);
      if (cigar_op->name == 'M' || cigar_op->name == 'I') {
	coverage += cigar_op->number;
      }
    }
  }

  return coverage;
}

//--------------------------------------------------------------------------------------

void init_cigar_string(cigar_code_t *p) {
  return;

  /*if (p->cigar_str) {
    free(p->cigar_str);
  }

  int num_ops = cigar_code_get_num_ops(p);
  char *str = malloc(num_ops * 5 * sizeof(char));
  *str = 0;

  cigar_op_t *op;
  for (int i = 0; i < num_ops; i++) {
    op = array_list_get(i, p->ops);
    sprintf(str, "%s%i%c", str, op->number, op->name);
  }

  p->cigar_str = str;*/
}

//--------------------------------------------------------------------------------------

cigar_code_t *generate_cigar_code(char *query_map, char *ref_map, unsigned int map_len,
				  unsigned int query_start, unsigned int query_len, 
				  int *distance) {
  
  cigar_code_t *p = cigar_code_new();

  char operation_number[map_len * 2];

  unsigned char status;
  unsigned char transition;
  short int cigar_soft;
  short int value = 0;
  unsigned int number_op = 0;
  char operation;
  unsigned int perfect = 0;  
  unsigned int deletions_tot = 0;
  
  int dist = 0;

  //  printf("seq(%d) start::%d : %s\n", length, start_seq, str_seq_p );
  //  printf("ref(%d): %s\n", length, str_ref_p);
  
  // hard clipping start
  if (query_start > 0){
    cigar_code_append_op(cigar_op_new(query_start, 'H'), p);
  }
  
  // first Status
  if (query_map[0] != '-' && ref_map[0] != '-') {
    status = CIGAR_MATCH_MISMATCH;
    // soft clipping
    cigar_soft = 0;
    while ((ref_map[cigar_soft] != '-') && (query_map[cigar_soft] != '-') && 
	   (ref_map[cigar_soft] != query_map[cigar_soft])) {
      cigar_soft++;
      value++;
    }
    if (value > 0) {
      cigar_code_append_op(cigar_op_new(value, 'S'), p);
    } 
  } else if (query_map[0] == '-') {
    if (ref_map[0] == '-') {
      status = CIGAR_PADDING;
    } else {
      status = CIGAR_DELETION;
    }
  } else if(ref_map[0] == '-') {
    status = CIGAR_INSERTION;
  }
  
  for (int i = value; i < map_len; i++) {
    // transition
    if (query_map[i] != '-' && ref_map[i] != '-') {
      transition = CIGAR_MATCH_MISMATCH;
      if (query_map[i] == ref_map[i]) {
        perfect++;
      } else {
	dist++;
      }
    } else if(query_map[i] == '-') {
      if (ref_map[i] == '-') {
        transition = CIGAR_PADDING;
      } else {
        transition = CIGAR_DELETION;
        deletions_tot++;
	dist++;
      }
    } else if (ref_map[i] == '-') {
      transition = CIGAR_INSERTION;
      dist++;
    }
    
    if (transition != status) {
      // insert operation in cigar string
      operation = select_op(status);
      cigar_code_append_op(cigar_op_new(number_op, operation), p);
      number_op = 1;
      status = transition;
    } else {
      number_op++;
    }
  }
  
  if ((map_len == perfect) && (perfect == query_len)) {
    status = CIGAR_PERFECT_MATCH;
  }
  
  operation = select_op(status);
  
  // hard and Soft clipped end
  if (status == CIGAR_MATCH_MISMATCH) {
    cigar_soft = map_len - 1;
    value = 0;
    while ((ref_map[cigar_soft] != '-') && (query_map[cigar_soft] != '-') && 
	   (query_map[cigar_soft] != ref_map[cigar_soft])){
      cigar_soft--;
      value++;
      //printf("(Soft %c!=%c)", output_p->mapped_ref_p[i][cigar_soft], output_p->mapped_seq_p[i][cigar_soft]);
    }
    
    cigar_code_append_op(cigar_op_new(number_op - value, operation), p);
    
    if (value > 0) {
      number_op -= value;
      cigar_code_append_op(cigar_op_new(value, 'S'), p);
    }
  } else {
    cigar_code_append_op(cigar_op_new(number_op, operation), p);
  }
  //printf("%d+%d < %d\n", length - deletions_tot, start_seq, seq_orig_len);
  if (((map_len - deletions_tot) + query_start) < query_len) {
    cigar_code_append_op(cigar_op_new(query_len - ((map_len - deletions_tot) + query_start), 'H'), p);
  }
  
  //printf("%d-%d\n", length, *number_op_tot);
  *distance = dist;

  init_cigar_string(p);

  return p;
}

/*
//--------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------

breakpoint_info_t *breakpoint_info_new(int side, int num_M, int index_M,
				       cigar_code_t *cigar) {
  breakpoint_info_t *p = calloc(1, sizeof(breakpoint_info_t));
  
  p->side = side;
  p->num_M = num_M;
  p->index_M = index_M;
  if (cigar) p->cigar = cigar;
  //  if (seq) p->seq = strdup(seq);
  
  return p;
}

//--------------------------------------------------------------------------------------

void breakpoint_info_free(breakpoint_info_t *p) {
  if (p) {
    if (p->cigar) cigar_code_free(p->cigar);
    //    if (p->seq) free(p->seq);

    free(p);
  }
}

//--------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------


breakpoint_t *breakpoint_new(int chr_index, int position) {
  breakpoint_t *p = calloc(1, sizeof(breakpoint_t));

  p->chr_index = chr_index;
  p->position = position;
  p->coverage = 1;
  p->info_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  return p;
}

//--------------------------------------------------------------------------------------

void breakpoint_free(breakpoint_t *p) {
  if (p) {
    if (p->info_list) array_list_free(p->info_list, (void *) breakpoint_info_free);

    free(p);
  }
}

//--------------------------------------------------------------------------------------

static inline void add_breakpoint_info(breakpoint_info_t *info, breakpoint_t *b) {
  array_list_insert((void *) info, b->info_list);
}

//--------------------------------------------------------------------------------------

breakpoint_t *compute_breakpoint(int chr_index, int alig_position, cigar_code_t *cigar,
				 array_list_t *list) {

  int m_left, m_left_index, m_right, m_right_index;
  get_matching_sides(&m_left, &m_left_index, &m_right, &m_right_index, cigar);

  int pos;
  int side = LEFT_SIDE;
  int num_M, index_M, index = cigar->num_ops - m_right_index -1;

  if ( (index < m_left_index) || 
       (index == m_left_index && m_right > m_left) ) {
    side = RIGHT_SIDE;
    pos = alig_position; //ends[i] - m_right;
    for (int i = 0; i < m_right_index; i++) {
      if (cigar->ops[i].name == 'M' ||
	  cigar->ops[i].name == 'S' ||
	  cigar->ops[i].name == 'D'   ) {
	pos += cigar->ops[i].number;
      }
    }
    num_M = m_right;
    index_M = m_right_index;
    
    
  } else {
    side == LEFT_SIDE;
    pos = alig_position + m_left;	  
    num_M = m_left;
    index_M = m_left_index;
  }

  if (!list) {
    list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  }

  breakpoint_t *breakpoint = NULL, *b = NULL;
  if (array_list_size(list)) {
    size_t num_items = array_list_size(list);
    for (size_t i = 0; i < num_items; i++) {
      b = (breakpoint_t *) array_list_get(i, list);
      if (b->chr_index == chr_index && b->position == pos) {
	b->coverage++;
	breakpoint = b;
	break;
      }
    }
  }
  if (!breakpoint) {
    breakpoint = breakpoint_new(chr_index, pos);
    array_list_insert((void *) breakpoint, list);
  }

  add_breakpoint_info(breakpoint_info_new(side, num_M, index_M, cigar), breakpoint);

  return breakpoint;
}

//--------------------------------------------------------------------------------------

breakpoint_t *get_breakpoint(int chr_index, size_t pos, array_list_t *list) {
  breakpoint_t *b, *breakpoint = NULL;

  if (list) {
    size_t num_items = array_list_size(list);
    for (size_t i = 0; i < num_items; i++) {
      b = (breakpoint_t *) array_list_get(i, list);
      if (b->chr_index == chr_index && b->position == pos) {
	b->coverage++;
	breakpoint = b;
	break;
      }
    }
  } 

  if (!breakpoint) {
    breakpoint = breakpoint_new(chr_index, pos);
    array_list_insert((void *) breakpoint, list);
  }
  return breakpoint;
}

//--------------------------------------------------------------------------------------

void display_breakpoints(array_list_t *list)  {

  if (list) {
    breakpoint_t *b;
    breakpoint_info_t *info;
    
    size_t num_infos, num_breaks = array_list_size(list);
    
    for (size_t i = 0; i < num_breaks; i++) {
      b = (breakpoint_t *) array_list_get(i, list);
      LOG_DEBUG_F("breakpoint (chr. index, pos.) = (%i, %lu)\tcoverage = %i\n",
		  b->chr_index, b->position, b->coverage);
      
      num_infos = array_list_size(b->info_list);
      
      for (size_t j = 0; j < num_infos; j++) {
	info = (breakpoint_info_t *) array_list_get(j, b->info_list);
	LOG_DEBUG_F("\t%iM at %i (%s)\t%s\t%s\n", 
		    info->num_M, info->index_M, (info->side == LEFT_SIDE ? "left" : "right"), 
		    info->cigar->cigar_str, ""); //info->seq);

      }
    }
  } else {
    LOG_DEBUG("Breakpont list is NULL !\n");
  }
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
/*
void get_matching_sides(int *m_left, int *m_left_index, 
			int *m_right, int *m_right_index, 
			cigar_code_t *cigar_code) {
  *m_left = 0;
  *m_right = 0;

  *m_left_index;
  *m_right_index;
  for (int j = 0; j < cigar_code->num_ops; j++) {
    if (cigar_code->ops[j].name == 'M') {
      *m_left = cigar_code->ops[j].number;
      *m_left_index = j;
      break;
    }
  }
  
  for (int j = cigar_code->num_ops - 1; j >= 0; j--) {
    if (cigar_code->ops[j].name == 'M') {
      *m_right = cigar_code->ops[j].number;
      *m_right_index = j;
      break;
    }
  } 
}
*/
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------


/*

    ///////////////////////////////////////////////////////////////////////////
    //
    // search for breakpoints
    {
      int v = -1;
      size_t alig_pos = starts[i] + output->ref_start_p[i] - 1;

      if (alig_pos >= 30000 && alig_pos <= 30100 && chromosomes[i] == 4 && strands[i] == 0) {
	printf("alig_pos = %lu\n", alig_pos);

      //      if (strands[i] == 0) {


	if (breakpoint_list == NULL) {
	  breakpoint_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
	}

	//      if (norm_score < min_score) {
	//      if (strstr(fq_read->id, "10000") != NULL) {
	
	int distance, num_cigar_ops;
	char *cigar, *p1, *p2;
	cigar =  generate_cigar_str(output->query_map_p[i], 
				    output->ref_map_p[i], 
				    output->query_start_p[i], 
				    read_len, 
				    strlen(output->query_map_p[i]), 
				    &distance,
				    (int *) &num_cigar_ops);

	cigar_code_t *cigar_code = cigar_code_new(cigar);


//	for (int j = 0; j < cigar_op_counter; j++) {
//	  printf("%i -> number: %i, operation: %c\n", j, cigar_ops[j].number, cigar_ops[j].name);
//	}

	int m_left = 0, m_right = 0;
	int m_left_index, m_right_index;
	for (int j = 0; j < cigar_code->num_ops; j++) {
	  if (cigar_code->ops[j].name == 'M') {
	    m_left = cigar_code->ops[j].number;
	    m_left_index = j;
	    break;
	  }
	}

	for (int j = cigar_code->num_ops - 1; j >= 0; j--) {
	  if (cigar_code->ops[j].name == 'M') {
	    m_right = cigar_code->ops[j].number;
	    m_right_index = cigar_code->num_ops - 1 - j;
	    break;
	  }
	}
	
	size_t breakpoint_pos;
	int side = LEFT_SIDE;
	int num_M = 0, index_M;

	if ( (m_right_index < m_left_index) || (m_right_index == m_left_index && m_right > m_left) ) {
	  side = RIGHT_SIDE;
	  breakpoint_pos = ends[i] - m_right;
	  num_M = m_right;
	  index_M = m_right_index;
	} else {
	  side == LEFT_SIDE;
	  breakpoint_pos = alig_pos + m_left;	  
	  num_M = m_left;
	  index_M = m_left_index;
	}

	breakpoint_t *breakpoint = get_breakpoint(chromosomes[i], breakpoint_pos, breakpoint_list);
	add_breakpoint_info(breakpoint_info_new(side, num_M, index_M, cigar_code, fq_read->sequence), breakpoint);

	if (1) {
	  //	if (breakpoint_pos == 20072) {
	  //	if (breakpoint_pos >= 20072 && breakpoint_pos >= 20112) {

	  if (v != read_index) {
	    LOG_DEBUG("\n");
	    LOG_DEBUG_F("%s\n", fq_read->id);
	    LOG_DEBUG_F("%s\n", fq_read->sequence);
	  }
	  LOG_DEBUG("\n");

	  if (side == LEFT_SIDE) {
	    LOG_DEBUG_F("\t%s\n", &(fq_read->sequence[output->query_start_p[i] + m_left]));
	    LOG_DEBUG_F("\tchr. index %i\tbreakpoint pos. %lu (left)\n", chromosomes[i], breakpoint_pos);	  
	  } else {
	    LOG_DEBUG_F("\tchr. index %i\tbreakpoint pos. %lu (right)\n", chromosomes[i], breakpoint_pos);
	  }

	  LOG_DEBUG("\n");
	  
	  LOG_DEBUG_F("\t\tcigar = %s (num. ops. %i), score = %0.2f, norm. score = %0.2f\n",
		      cigar, num_cigar_ops, output->score_p[i], norm_score);
	  LOG_DEBUG_F("\t\tCAL (chr, strand, start, end) = (%i, %i, %lu, %lu) -> alignment pos. = %lu\n", 
		      chromosomes[i], strands[i], starts[i], ends[i], alig_pos);
	  LOG_DEBUG_F("\t\tM-left (%i, index %i), M-rigth (%i, index %i)\n",
		      m_left, m_left_index, m_right, m_right_index);
	  LOG_DEBUG_F("\t\t\tref. map : %s (start: %i)\n", 
		      output->ref_map_p[i], output->ref_start_p[i]);
	  LOG_DEBUG_F("\t\t\tquery map: %s (start: %i)\n\n", 
		      output->query_map_p[i], output->query_start_p[i]);
	}

	free(cigar);

	v = read_index;
      }
    }
    ///////////////////////////////////////////////////////////////////////////


 */
