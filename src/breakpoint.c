#include "breakpoint.h"

//--------------------------------------------------------------------------------------

array_list_t *breakpoint_list = NULL;

//--------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------

cigar_op_t *cigar_op_new(int number, char number_str[], char name) {
  cigar_op_t *cigar_op = (cigar_op_t *)malloc(sizeof(cigar_op_t));
  cigar_op->number = number;
  if (number_str) {
    strcpy(cigar_op->number_str, number_str);
  }
  cigar_op->name = name;

  return cigar_op;
}

void cigar_op_free(cigar_op_t *cigar_op) {
  free(cigar_op);
}

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
		    info->cigar->cigar_str, ""/*info->seq*/);
      }
    }
  } else {
    LOG_DEBUG("Breakpont list is NULL !\n");
  }
}

//--------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------

cigar_code_t *cigar_code_new(char *cigar_str) {

  cigar_code_t *p = (cigar_code_t *) calloc(1, sizeof(cigar_code_t));

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

  return p;
}

//--------------------------------------------------------------------------------------

void cigar_code_free(cigar_code_t* p) {
  if (p) {
    if (p->ops) free(p->ops);
    if (p->cigar_str) free(p->cigar_str);
    free(p);
  }
}

//--------------------------------------------------------------------------------------

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
