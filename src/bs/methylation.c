#include "methylation.h"

const size_t margin_bwt = 10;
const size_t margin_sw  = 25;

//====================================================================================

void replace(char * refs, int len, int type){
  char c1, c2;
  switch (type) {
  case ACGT:
    return;
  case AGT:
    c1 = 'C';
    c2 = 'T';
    break;
  case ACT:
    c1 = 'G';
    c2 = 'A';
    break;
  case AT:
    replace(refs, len, 1);
    replace(refs, len, 2);
    return;
  default:
    return;
  }

  for (int j = 0; j < len; j++) {
    if (refs[j] == c1) {
      refs[j] = c2;
    }
  }

  return;
}

//====================================================================================

char complement (char c) {
  switch (c) {
  case 'A':
    return 'T';
  case 'T':
    return 'A';
  case 'C':
    return 'G';
  case 'G':
    return 'C';
  default:
    return c;
  }
}

//====================================================================================

void rev_comp(char *query, char *seq, int len) {
  for (int i = 0; i < len; i++) {
    seq[len - i - 1] = complement(query[i]);
  }
  seq[len] = '\0';
}

void comp(char *seq, int len) {
  for (int i = 0; i < len - 1; i++) {
    seq[i] = complement(seq[i]);
  }
}

//====================================================================================

void replace_array(array_list_t *reads, int type) {
  size_t num_reads = array_list_size(reads);
  fastq_read_t* fq_read;
  
  for (size_t i = 0; i < num_reads; i++) {
    fq_read = (fastq_read_t *) array_list_get(i, reads);

    replace(fq_read->sequence, fq_read->length, type);
  }
}

//====================================================================================

void rev_comp_array(array_list_t *dest, array_list_t *src) {
  size_t num_reads = array_list_size(src);
  fastq_read_t* fq_read_src;
  fastq_read_t* fq_read_dest;
  
  for (size_t i = 0; i < num_reads; i++) {
    fq_read_src  = (fastq_read_t *) array_list_get(i, src);
    fq_read_dest = (fastq_read_t *) array_list_get(i, dest);

    rev_comp(fq_read_src->sequence, fq_read_dest->sequence, fq_read_src->length);
  }
}

//====================================================================================

void copy_array(array_list_t *dest, array_list_t *src) {
  
}
