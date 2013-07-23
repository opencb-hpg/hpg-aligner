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

    /*
    if (strcmp(fq_read_src->id, fq_read_dest->id)) {
      printf("read %lu\n");
    }
    */

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

char *obtain_seq(alignment_t *alig) {
  char *read = alig->sequence;
  char *cigar = strdup(alig->cigar);
  //char *cigar = strdup("100M1D");
  int num;
  char car;
  int cont, pos, pos_read;
  int operations;
  int len = strlen(cigar) - 1;
  char *seq = (char *)calloc(1024, sizeof(char));
  //printf("cigar %s\n", cigar);
  //printf("read  %s\n", read);
  
  pos = 0;
  pos_read = 0;
  
  for (operations = 0; operations < alig->num_cigar_operations; operations++) {
    sscanf(cigar, "%i%c%s", &num, &car, cigar);
    //printf("%3i %c %s\n",  num, car, cigar);
    if (car == 'M' || car == '=') {
      for (cont = 0; cont < num; cont++, pos++, pos_read++) {
	seq[pos] = read[pos_read];
      }
    }
    else {
      if (car == 'D') {
	pos_read += num - 1;
      }
      else {
	if (car == 'I') {
	  for (cont = 0; cont < num; cont++, pos++) {
	    seq[pos] = 'N';
	  }
	}
      }
    }
  }
  seq[pos] = '\0';

  //printf("seq   %s\n", seq);

  free(cigar);
  return seq;
}

//====================================================================================

void write_metilation_status(array_list_t *array_list, metil_file_t *metil_file) {

  //printf("Init metilation status\n");
  
  size_t num_items = array_list_size(array_list);
  alignment_t *alig;
  genome_t *genome = metil_file->genome;
  hash_table_t *table_bs = metil_file->table_isles;
  char *seq, *gen, *read;
  size_t len, end, start;
  char *key = (char *)malloc(3 * sizeof(char));
  int data;

  char *cigar;
  int num;
  char car;
  int cont, pos, pos_read;

  for (size_t j = 0; j < num_items; j++) {
    alig = (alignment_t *) array_list_get(j, array_list);
    if (alig != NULL && alig->is_seq_mapped) {
      /*
      printf("alignment %lu\n", j);
      printf("query_name %s\n", alig->query_name);
      printf("sequence   %s\n", alig->sequence);
      printf("cigar      %s\n", alig->cigar);
      printf("position   %i\n", alig->position);
      printf("mate pos   %i\n", alig->mate_position);
      printf("temp len   %i\n", alig->template_length);
      printf("chromo     %i\n", alig->chromosome);
      printf("strand     %i\n", alig->seq_strand);
      printf("mate  chro %i\n", alig->mate_chromosome);
      printf("map qual   %i\n", alig->map_quality);
      printf("n cigar    %i\n", alig->num_cigar_operations);
      */      
      
      //seq = obtain_seq(alig);
      //len = strlen(seq);
      len = strlen(alig->sequence);
      gen = (char *)calloc(len + len + 2, sizeof(char));

      start = alig->position + 1;
      end = start + len - 1;
      if (end >= genome->chr_size[alig->chromosome]) {
        end = genome->chr_size[alig->chromosome] - 1;
      }

      genome_read_sequence_by_chr_index(gen, alig->seq_strand, alig->chromosome, &start, &end, genome);

      //printf("\nseq %s\ngen %s\n", seq, gen);

      read = alig->sequence;
      cigar = strdup(alig->cigar);
      for (int operations = 0; operations < alig->num_cigar_operations; operations++) {
	sscanf(cigar, "%i%c%s", &num, &car, cigar);
	//printf("%3i %c %s\n",  num, car, cigar);
	if (car == 'M' || car == '=') {
	  for (cont = 0; cont < num; cont++, pos++, pos_read++) {
	    //seq[pos] = read[pos_read];
	    if (read[pos_read] == 'C') {
	      if (gen[pos] == 'C') {
		// methylated cytosine
		if (gen[pos + 1] == 'G') {
		  // CpG zone
		  ;
		  //printf("%s\t+\t%i %i\t%lu\tZ\n", alig->query_name, alig->chromosome, alig->seq_strand, start + pos);
		  //fprintf(metil_file->CpG, "%s\t+\t%i %i\t%lu\tZ\n",
		  //	  alig->query_name, alig->chromosome, alig->seq_strand, start + pos);
		}
		else {
		  if (gen[pos + 2] == 'G') {
		    // CHG zone
		    ;
		    //printf("%s\t+\t%i %i\t%lu\tX\n", alig->query_name, alig->chromosome, alig->seq_strand, start + pos);
		    //fprintf(metil_file->CHG, "%s\t+\t%i %i\t%lu\tX\n",
		    //	    alig->query_name, alig->chromosome, alig->seq_strand, start + pos);
		  }
		  else {
		    // CHH zone
		    ;
		    //printf("%s\t+\t%i %i\t%lu\tH\n", alig->query_name, alig->chromosome, alig->seq_strand, start + pos);
		    //fprintf(metil_file->CHH, "%s\t+\t%i %i\t%lu\tH\n",
		    //	    alig->query_name, alig->chromosome, alig->seq_strand, start + pos);
		  }
		}
	      }
	      else {
		// mutated cytosine
		;
		//printf("%s\t+\t%i %i\t%lu\tM\n", alig->query_name, alig->chromosome, alig->seq_strand, start + pos);
		//fprintf(metil_file->MUT, "%s\t+\t%i %i\t%lu\tM\n",
		//	alig->query_name, alig->chromosome, alig->seq_strand, start + pos);
	      }
	    }
	  }
	}
	else {
	  if (car == 'D') {
	    pos_read += num - 1;
	  }
	  else {
	    if (car == 'I') {
	      for (cont = 0; cont < num; cont++, pos++) {
		//seq[pos] = 'N';
	      }
	    }
	  }
	}
      }

      if (seq) free(seq);
      if (gen) free(gen);
    }
  }
  if (key) free(key);
  //printf("end status\n");
}

//====================================================================================

void metil_file_init(metil_file_t *metil_file, char *CpG, char *CHG, char *CHH, char *MUT, genome_t *genome) {
  metil_file->filenameCpG = strdup(CpG);
  metil_file->filenameCHG = strdup(CHG);
  metil_file->filenameCHH = strdup(CHH);
  metil_file->filenameMUT = strdup(MUT);

  metil_file->genome = genome;

  metil_file->CpG = fopen(metil_file->filenameCpG, "w");
  metil_file->CHG = fopen(metil_file->filenameCHG, "w");
  metil_file->CHH = fopen(metil_file->filenameCHH, "w");
  metil_file->MUT = fopen(metil_file->filenameMUT, "w");

  metil_file->table_isles = hash_create(jenkins_one_at_a_time_hash,
					scmp,
					destr_key,
					nulldes,
					200);
  hash_init(metil_file->table_isles, "ACGTN-");
}

//====================================================================================

void metil_file_free(metil_file_t *metil_file) {
  free(metil_file->filenameCpG);
  free(metil_file->filenameCHG);
  free(metil_file->filenameCHH);
  free(metil_file->filenameMUT);

  fclose(metil_file->CpG);
  fclose(metil_file->CHG);
  fclose(metil_file->CHH);
  fclose(metil_file->MUT);

  hash_destroy(metil_file->table_isles);

  free(metil_file);
}

//====================================================================================

void write_metilation_status_old(array_list_t *array_list, metil_file_t *metil_file) {

  printf("Init metilation status\n");
  
  size_t num_items = array_list_size(array_list);
  alignment_t *alig;
  genome_t *genome = metil_file->genome;
  hash_table_t *table_bs = metil_file->table_isles;
  char *seq, *gen;
  size_t len, end, start;
  char *key = (char *)malloc(3 * sizeof(char));
  int data;

  for (size_t j = 0; j < num_items; j++) {
    alig = (alignment_t *) array_list_get(j, array_list);
    if (alig != NULL && alig->is_seq_mapped) {
      /*
      printf("alignment %lu\n", j);
      printf("query_name %s\n", alig->query_name);
      printf("sequence   %s\n", alig->sequence);
      printf("cigar      %s\n", alig->cigar);
      printf("position   %i\n", alig->position);
      printf("mate pos   %i\n", alig->mate_position);
      printf("temp len   %i\n", alig->template_length);
      printf("chromo     %i\n", alig->chromosome);
      printf("strand     %i\n", alig->seq_strand);
      printf("mate  chro %i\n", alig->mate_chromosome);
      printf("map qual   %i\n", alig->map_quality);
      printf("n cigar    %i\n", alig->num_cigar_operations);
      */      
      
      seq = obtain_seq(alig);
      //printf("seq %s\n", seq);

      len = strlen(seq);
      gen = (char *)calloc(len + 2, sizeof(char));

      start = alig->position;
      end = start + len - 1;
      if (end >= genome->chr_size[alig->chromosome]) {
	//printf("    end %lu\n", end);
	//printf("    gen %lu\n", genome->chr_size[alig->chromosome]);
        end = genome->chr_size[alig->chromosome] - 1;
	//printf("new end %lu\n", end);
      }


      //printf("seq %s\n", seq);
      //printf("chromo %i, strand %i, begin %lu, end %lu\n",
      //     alig->chromosome, alig->seq_strand, alig->position, end);
      
      genome_read_sequence_by_chr_index(gen, alig->seq_strand, alig->chromosome, &start, &end, genome);

      printf("\nseq %s\ngen %s\n", seq, gen);

      for (size_t i = 0; i < len - 2; i++) {
	if (gen[i] == 'C') {
	  if (seq[i] == 'C' || seq[i] == 'T') {
	    key[0] = gen[i];
	    key[1] = gen[i + 1];
	    key[2] = gen[i + 2];
	    data = (int)hash_find_data(table_bs, key);
	    
	    //printf("Candidata (%i) en %lu\n", data, start + i);

	    switch (data) {
	    case ZONE_CpG:
	      if (seq[i] == 'C') {
		printf("%s\t+\t%i %i\t%lu\tZ\n", alig->query_name, alig->chromosome, alig->seq_strand, start + i);
	      } else {
		printf("%s\t-\t%i %i\t%lu\tz\n", alig->query_name, alig->chromosome, alig->seq_strand, start + i);
	      }
	      break;
	    case ZONE_CHG:
	      if (seq[i] == 'C') {
		printf("%s\t+\t%i %i\t%lu\tX\n", alig->query_name, alig->chromosome, alig->seq_strand, start + i);
	      } else {
		printf("%s\t-\t%i %i\t%lu\tx\n", alig->query_name, alig->chromosome, alig->seq_strand, start + i);
	      }
	      break;
	    case ZONE_CHH:
	      if (seq[i] == 'C') {
		printf("%s\t+\t%i %i\t%lu\tH\n", alig->query_name, alig->chromosome, alig->seq_strand, start + i);
	      } else {
		printf("%s\t-\t%i %i\t%lu\th\n", alig->query_name, alig->chromosome, alig->seq_strand, start + i);
	      }
	      break;
	    case ZONE_OTHER:
	      if (seq[i] == 'C') {
		printf("%s\t+\t%i %i\t%lu\tM\n", alig->query_name, alig->chromosome, alig->seq_strand, start + i);
	      } else {
		printf("%s\t-\t%i %i\t%lu\tm\n", alig->query_name, alig->chromosome, alig->seq_strand, start + i);
	      }
	      break;
	    default:
	      printf("%s\t-\t%i %i\t%lu\t???\n", alig->query_name, alig->chromosome, alig->seq_strand, start + i);
	    }
	  }
	  else {
	    printf("%s\t-\t%i %i\t%lu\tm\n", alig->query_name, alig->chromosome, alig->seq_strand, start + i);
	    //printf("Mutacion en %lu\n", start + i);
	  }
	}
	;
      }

      if (seq) free(seq);
      if (gen) free(gen);
    }
  }
  free(key);
}

//====================================================================================

char *obtain_seq_old(alignment_t *alig) {
  char *read = alig->sequence;
  char *cigar = strdup(alig->cigar);
  //char *cigar = strdup("100M1D");
  int num;
  char car;
  int cont, pos, pos_read;
  int len = strlen(cigar) - 1;
  char *seq = (char *)calloc(1024, sizeof(char));
  //printf("cigar %s\n", cigar);
  //printf("read  %s\n", read);
  
  pos = 0;
  pos_read = 0;
  while(strlen(cigar) > 0 && strlen(cigar) != len) {
    len = strlen(cigar);
    sscanf(cigar, "%i%c%s", &num, &car, cigar);
    //printf("%3i %c %s\n",  num, car, cigar);
    if (car == 'M' || car == '=') {
      for (cont = 0; cont < num; cont++, pos++, pos_read++) {
	seq[pos] = read[pos_read];
      }
    }
    else {
      if (car == 'D') {
	pos_read += num - 1;
      }
      else {
	if (car == 'I') {
	  for (cont = 0; cont < num; cont++, pos++) {
	    seq[pos] = 'N';
	  }
	}
      }
    }
  }
  seq[pos] = '\0';

  //printf("seq   %s\n", seq);

  free(cigar);
  return seq;
}

//====================================================================================

