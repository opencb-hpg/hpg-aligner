#include "methylation.h"

const size_t margin_bwt = 10;
const size_t margin_sw  = 25;

//====================================================================================

void replace(char * refs, int len, int type) {
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

  //printf("replace\nleng = %lu\ntype = %lu\n", len, type);
  //printf("c1 = %c\nc2 = %c\n\n", c1, c2);

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

void rev_comp(char *orig, char *dest, int len) {
  for (int i = 0; i < len; i++) {
    dest[len - i - 1] = complement(orig[i]);
  }
  dest[len] = '\0';
}

//====================================================================================

void comp(char *seq, int len) {
  for (int i = 0; i < len - 1; i++) {
    seq[i] = complement(seq[i]);
  }
}

//====================================================================================

void replace_array(array_list_t *reads, int type) {
  size_t num_reads = array_list_size(reads);
  fastq_read_t* fq_read;

  //printf("reads = %lu\n", num_reads);
  
  for (size_t i = 0; i < num_reads; i++) {
    fq_read = (fastq_read_t *) array_list_get(i, reads);

    //printf("read = %lu\n", i);
    //printf("src = %s\n", fq_read->sequence);
    replace(fq_read->sequence, fq_read->length, type);
    //printf("src = %s\ntam = %lu\ntype = %lu\n\n", fq_read->sequence, fq_read->length, type);
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

    //printf("read = %lu\n", i);
    //printf("src = %s\ndst = %s\ntam = %lu\n\n", fq_read_src->sequence, fq_read_dest->sequence, fq_read_src->length);
    rev_comp(fq_read_src->sequence, fq_read_dest->sequence, fq_read_src->length);
    //printf("src = %s\ndst = %s\ntam = %lu\n\n", fq_read_src->sequence, fq_read_dest->sequence, fq_read_src->length);
  }
}

//====================================================================================

void copy_array(array_list_t *dest, array_list_t *src) {
  size_t num_reads = array_list_size(src);
  fastq_read_t* fq_read_src;
  fastq_read_t* fq_read_dest;
  
  for (size_t i = 0; i < num_reads; i++) {
    fq_read_src  = (fastq_read_t *) array_list_get(i, src);
    //fq_read_dest = fastq_read_new(fq_read_src);
    //insert_in_array(dest, fq_read_dest);
  }
}

//====================================================================================

void cpy_array_bs(array_list_t *src, array_list_t *dest1, array_list_t *dest2, array_list_t *dest3, array_list_t *dest4) {
  size_t num_reads = array_list_size(src);
  fastq_read_t* fq_read_src;
  fastq_read_t* fq_read_dest;
  
  for (size_t i = 0; i < num_reads; i++) {
    fq_read_src  = (fastq_read_t *) array_list_get(i, src);
    fq_read_dest = fastq_read_dup(fq_read_src);
    array_list_insert(fq_read_dest, dest1);
    fq_read_dest = fastq_read_dup(fq_read_src);
    array_list_insert(fq_read_dest, dest2);
    fq_read_dest = fastq_read_dup(fq_read_src);
    array_list_insert(fq_read_dest, dest3);
    fq_read_dest = fastq_read_dup(fq_read_src);
    array_list_insert(fq_read_dest, dest4);
  }
}

//====================================================================================

void cpy_transform_array_bs(array_list_t *src, array_list_t *dest_ct, array_list_t *dest_ct_rev, array_list_t *dest_ga, array_list_t *dest_ga_rev) {
  size_t num_reads = array_list_size(src);
  fastq_read_t *fq_read_src;
  fastq_read_t *fq_read_dest;
  fastq_read_t *fq_read_tmp;

  for (size_t i = 0; i < num_reads; i++) {
    fq_read_src  = (fastq_read_t *) array_list_get(i, src);

    fq_read_dest = fastq_read_dup(fq_read_src);
    replace(fq_read_dest->sequence, fq_read_dest->length, AGT);
    array_list_insert(fq_read_dest, dest_ct);

    fq_read_tmp = fastq_read_dup(fq_read_src);
    rev_comp(fq_read_dest->sequence, fq_read_tmp->sequence, fq_read_dest->length);
    array_list_insert(fq_read_tmp, dest_ct_rev);

    fq_read_dest = fastq_read_dup(fq_read_src);
    replace(fq_read_dest->sequence, fq_read_dest->length, ACT);
    array_list_insert(fq_read_dest, dest_ga);

    fq_read_tmp = fastq_read_dup(fq_read_src);
    rev_comp(fq_read_dest->sequence, fq_read_tmp->sequence, fq_read_dest->length);
    array_list_insert(fq_read_tmp, dest_ga_rev);
  }
}

//====================================================================================

void insert_mappings_array(array_list_t **dest, array_list_t **src) {
  size_t num_reads = array_list_size(src);
  size_t num_mappings;
  alignment_t *align_tmp;

  for (size_t i = 0; i < num_reads; i++) {
    //isnert_mappings(dest[i], src[i]);
    num_mappings = array_list_size(src[i]);
    for (size_t j = 0; j < num_mappings; j++) {
      align_tmp  = (alignment_t *) array_list_get(j, src[i]);
      // insert the alignments from the 'src' into 'dest'
      array_list_insert(align_tmp, dest[i]);
    }
  }
}

//====================================================================================

void insert_mappings(array_list_t *dest, array_list_t *src) {
  size_t num_mappings;
  alignment_t *align_tmp;

  num_mappings = array_list_size(src);
  for (size_t j = 0; j < num_mappings; j++) {
    align_tmp  = (alignment_t *) array_list_get(j, src);
    // insert the alignments from the 'src' into the 'dest'
    array_list_insert(align_tmp, dest);
  }
}

//====================================================================================

void insert_regions(array_list_t *dest, array_list_t *src) {
  size_t num_mappings;
  region_t *region_tmp;

  num_mappings = array_list_size(src);
  for (size_t j = 0; j < num_mappings; j++) {
    region_tmp  = (region_t *) array_list_get(j, src);
    // insert the alignments from the 'src' into the 'dest'
    array_list_insert(region_tmp, dest);
  }
}

//====================================================================================

void transform_mappings_array(array_list_t **src){
  size_t num_reads = array_list_size(src);
  size_t num_mappings;
  alignment_t *align_tmp;

  for (size_t i = 0; i < num_reads; i++) {
    num_mappings = array_list_size(src[i]);
    for (size_t j = 0; j < num_mappings; j++) {
      align_tmp  = (alignment_t *) array_list_get(j, src[i]);
      // modify the strand, and make it reverse
      align_tmp->seq_strand = 1;
    }
  }
}

//====================================================================================

void transform_mappings(array_list_t *src){
  size_t num_mappings;
  alignment_t *align_tmp;

  num_mappings = array_list_size(src);
  for (size_t j = 0; j < num_mappings; j++) {
    align_tmp  = (alignment_t *) array_list_get(j, src);
    // modify the strand, and make it reverse
    align_tmp->seq_strand = 1;
  }
}

//====================================================================================
void transform_regions(array_list_t *src){
  size_t num_mappings;
  region_t *region_tmp;

  num_mappings = array_list_size(src);
  for (size_t j = 0; j < num_mappings; j++) {
    region_tmp  = (region_t *) array_list_get(j, src);
    // modify the strand, and make it reverse
    region_tmp->strand = 1;
  }
}

//====================================================================================

void select_targets_bs(size_t *num_unmapped, size_t *unmapped_indices,
		       size_t *indices_1, size_t num_reads,
		       array_list_t **lists, array_list_t **lists2) {

  *num_unmapped = 0;

  for (size_t i = 0; i < num_reads; i++) {
    if (indices_1[i] == 4) {
      unmapped_indices[(*num_unmapped)++] = indices_1[i];
      array_list_set_flag(0, lists[indices_1[i]]);
      array_list_set_flag(0, lists2[indices_1[i]]);
    }
  }
}

//====================================================================================

void update_targets_bs(size_t num_reads, size_t *unmapped_indices,
		       size_t unmapped, size_t *indices) {
  for (size_t i = 0; i < unmapped; i++) {
    if (indices[i] < num_reads) {
      unmapped_indices[indices[i]]++;
    }
  }
}

//====================================================================================

void revert_mappings_seqs(array_list_t **src1, array_list_t **src2, array_list_t *orig) {
  //printf("+++++++++++++++ is orig NULL ? %i, size %i\n", (orig == NULL), array_list_size(orig));

  size_t num_mappings;
  alignment_t *align_tmp;
  fastq_read_t *fastq_orig;
  size_t num_reads = array_list_size(orig);

  //printf("num reads = %lu\t", num_reads);

  //printf("num reads = %lu\t", array_list_size(orig));
  //printf("num reads src1 = %lu\t", array_list_size(src1));
  //printf("num reads src2 = %lu\t", array_list_size(src2));

  for (size_t i = 0; i < num_reads; i++) {
    fastq_orig  = (fastq_read_t *) array_list_get(i, orig);
 
    // mappings for reads transformations 1
    num_mappings = array_list_size(src1[i]);
    //printf("mapps to convert 1 (read %lu) = %lu\n", i, num_mappings);
    for (size_t j = 0; j < num_mappings; j++) {
      align_tmp  = (alignment_t *) array_list_get(j, src1[i]);
      // modify the strand, and make it reverse
      if (align_tmp->sequence != NULL) free(align_tmp->sequence);

      align_tmp->sequence = strdup(fastq_orig->sequence);
    }

    // mappings for reads transformations 2
    num_mappings = array_list_size(src2[i]);
    //printf("mapps to convert 1 (read %lu) = %lu\n", i, num_mappings);
    for (size_t j = 0; j < num_mappings; j++) {
      align_tmp  = (alignment_t *) array_list_get(j, src2[i]);
      // modify the strand, and make it reverse
      if (align_tmp->sequence) {
	free(align_tmp->sequence);
      }
      align_tmp->sequence = strdup(fastq_orig->sequence);
    }
  }
}

//====================================================================================

