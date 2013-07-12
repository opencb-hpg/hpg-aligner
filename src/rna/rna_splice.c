#include "rna_splice.h"

size_t junction_id = 0;
size_t total_splice = 0;

int search_end_splice(node_element_splice_t *node, size_t end, unsigned char strand) {
  for (int i = 0; i < node->number_allocate_ends; i++) {
    //printf("%i == %i, %i == %i\n", node->allocate_ends[i]->end, end, node->allocate_ends[i]->strand, strand);
    if ((node->allocate_ends[i]->end == end) && 
	(node->allocate_ends[i]->strand == strand)) { 
      return 1;
    }
  }
  return 0;
}


void load_intron_file(genome_t *genome, char* intron_filename, allocate_splice_elements_t *avls) {
  FILE *fd_intron = fopen(intron_filename, "r");
  size_t max_size = 2048;
  size_t line_size;
  char line[max_size];
  char value[max_size];
  int pos, value_pos;
  size_t chr, start, end;
  unsigned char strand;
  int found;
  node_element_splice_t *node;

  if (fd_intron == NULL) {
    return;
    //LOG_FATAL_F("Error opening file: %s, mode (%s)\n", intron_filename, "r");
  }
  
  while (!feof(fd_intron)) {
    fgets(line, max_size, fd_intron);
    line_size = strlen(line);
    pos = 0;
    value_pos = 0;
    
    while (line[pos] != '\t' && pos < line_size) {
      value[value_pos++] = line[pos++];
    }
    value[value_pos] = '\0';

    found = 0;
    for (chr = 0; chr < genome->num_chromosomes; chr++) {
      //printf("%s == %s\n", value, genome->chr_name[chr]);
      if (strcmp(value, genome->chr_name[chr]) == 0) { found = 1; break; }
    }

    if (!found) { continue; }
    
    pos++;
    value_pos = 0;
    while (line[pos] != '\t' && pos < line_size) {
      value[value_pos++] = line[pos++];
    }
    value[value_pos] = '\0';
    start = atoi(value);

    pos++;
    value_pos = 0;
    while (line[pos] != '\t' && pos < line_size) {
      value[value_pos++] = line[pos++];
    }
    value[value_pos] = '\0';
    end = atoi(value);

    pos++;
    value_pos = 0;
    while (line[pos] != '\t' && pos < line_size) {
      value[value_pos++] = line[pos++];
    }
    value[value_pos] = '\0';
    if (strcmp(value, "-1") == 0 || strcmp(value, "-")) {
      strand = 1;
    } else { strand = 0; }

    node = cp_avltree_get(avls[chr].avl_splice, (void *)start);
    found = 0;
    if (node) {
      found = search_end_splice(node, end, strand);
    } 

    if (!found) {    
      allocate_new_splice(chr, strand, end, start, start, end, FROM_FILE, avls);
    }
  }

}

int node_compare(node_element_splice_t* a, size_t b) {
  if(a->splice_start == b){ 
    return 0;
  }else if(a->splice_start < b){
   return -1;
  }else{
    return 1;
  }
}

node_element_splice_t* node_copy(size_t b) {
 node_element_splice_t *node = (node_element_splice_t *)malloc(sizeof(node_element_splice_t));

 //assert(node != NULL);
 
 node->maximum_allocate_ends = 10;
 node->number_allocate_ends = 0;
 node->allocate_ends = (splice_end_t **)malloc(node->maximum_allocate_ends * sizeof(splice_end_t *));
 if (node->allocate_ends == NULL){ exit(-1); }
 
 node->splice_start = b;

 return node;
}

void node_free(node_element_splice_t* a) {  
  for (unsigned int i = 0; i < a->number_allocate_ends; i++) {
    free_splice_end(a->allocate_ends[i]);
  }

  free(a->allocate_ends);
  free(a);
}


node_element_splice_t* insert_end_splice(splice_end_t *splice_end_p, node_element_splice_t *element_p){
  element_p->allocate_ends[element_p->number_allocate_ends] = splice_end_p;
  element_p->number_allocate_ends++;

  if(element_p->number_allocate_ends >= element_p->maximum_allocate_ends){
    element_p->maximum_allocate_ends = element_p->maximum_allocate_ends * 2; 
    element_p->allocate_ends = (splice_end_t **) realloc (element_p->allocate_ends, element_p->maximum_allocate_ends * sizeof (splice_end_t *));
    if(element_p->allocate_ends == NULL){exit(-1);}
  }

  return element_p;
}


node_element_splice_t* search_and_insert_end_splice(unsigned int chromosome, unsigned char strand, 
						    size_t end, size_t splice_start, 
						    size_t splice_end, int type_orig, 
						    node_element_splice_t *element_p){
  unsigned int i;

  if(element_p->splice_start_extend > splice_start){
    element_p->splice_start_extend = splice_start;
  }	  
  
  for(i = 0; i < element_p->number_allocate_ends; i++){
    if( (element_p->allocate_ends[i]->end == end) && (element_p->allocate_ends[i]->strand == strand) ) { 

      element_p->allocate_ends[i]->reads_number++;
      
      if(element_p->allocate_ends[i]->splice_end_extend < splice_end){
	element_p->allocate_ends[i]->splice_end_extend = splice_end;
      }  
	  
      return element_p;
    }
  }
  
  splice_end_t *splice_end_p = new_splice_end(strand, end, type_orig, splice_end);
  
  return insert_end_splice(splice_end_p, element_p);
}


allocate_splice_elements_t* init_allocate_splice_elements(allocate_splice_elements_t* chromosomes_avls_p, size_t nchromosomes){
   int i;
   for(i = 0; i < nchromosomes; i++){
     chromosomes_avls_p[i].avl_splice = cp_avltree_create_by_option(COLLECTION_MODE_NOSYNC | 
								    COLLECTION_MODE_COPY   |
								    COLLECTION_MODE_DEEP, 
								    (cp_compare_fn) node_compare, 
								    (cp_copy_fn) node_copy, 
								    (cp_destructor_fn)node_free, 
								    (cp_copy_fn) node_copy, 
								    (cp_destructor_fn)node_free);
     
     if(chromosomes_avls_p[i].avl_splice == NULL) {exit(-1);}
     pthread_mutex_init(&(chromosomes_avls_p[i].mutex), NULL);
   }
   
   return chromosomes_avls_p;
}

allocate_splice_elements_t* allocate_new_splice(unsigned int chromosome, unsigned char strand, 
						size_t end, size_t start, 
						size_t splice_start, size_t splice_end, int type_orig,
						allocate_splice_elements_t* chromosome_avls_p){
  node_element_splice_t *node;

  node = (node_element_splice_t *)cp_avltree_get(chromosome_avls_p[chromosome].avl_splice, (void *)start);

  if(node == NULL) {
    node = cp_avltree_insert(chromosome_avls_p[chromosome].avl_splice, (void *)start, (void *)start);
    node->splice_start_extend = splice_start;
  }
  
  node = search_and_insert_end_splice(chromosome, strand, end, splice_start, splice_end, type_orig, node);

  return chromosome_avls_p;
}

splice_end_t* new_splice_end(unsigned char strand, size_t end, 
			     int type_orig, size_t splice_end) {

  splice_end_t* splice_end_p = (splice_end_t *)malloc(sizeof(splice_end_t));

  if(splice_end_p == NULL) {exit(-1);}

  splice_end_p->strand = strand;
  splice_end_p->end = end;
  splice_end_p->splice_end_extend = splice_end;
  
  if (type_orig == FROM_READ) {
    splice_end_p->reads_number = 1;
  } else {
    splice_end_p->reads_number = 0;
  }

  return splice_end_p;
}


void free_splice_end(splice_end_t *splice_end_p){
  free(splice_end_p);
}


allocate_buffers_t * process_avlnode_in_order(cp_avlnode *node, unsigned int chromosome, 
					      list_t* write_list_p, unsigned int write_size,   allocate_buffers_t *allocate_batches){  
  if (node->left) { 
    allocate_batches = process_avlnode_in_order(node->left, chromosome, write_list_p, write_size, allocate_batches);
  }
  
  allocate_batches = process_avlnode_ends_in_order((node_element_splice_t *)node->value, chromosome, write_list_p, write_size,
						   allocate_batches);
  
  if (node->right) { 
    allocate_batches = process_avlnode_in_order(node->right, chromosome,  write_list_p, write_size, allocate_batches);
  }
  
  return allocate_batches;
  // return exact_splice_write_p;

}

allocate_buffers_t* process_avlnode_ends_in_order(node_element_splice_t *node, unsigned int chromosome,
						  list_t* write_list_p, unsigned int write_size, 
						  allocate_buffers_t *allocate_batches) {
  int i;
  char strand[2] = {'+', '-'};
  list_item_t* item_p = NULL;
  unsigned int bytes_exact, bytes_extend;
  allocate_batches->write_exact_sp;
  //  write_batch_t* extend_splice_write_p = write_batch_new(write_size, SPLICE_EXTEND_FLAG);

  //printf("----------------->%i\n", chromosome);
  for(i = 0; i < node->number_allocate_ends; i++){
    if(( allocate_batches->write_exact_sp->size + 100) > write_size) {
      //item_p = list_item_new(0, WRITE_ITEM,  allocate_batches->write_exact_sp);
      //list_insert_item(item_p, write_list_p);
      //allocate_batches->write_exact_sp = write_batch_new(write_size, SPLICE_EXACT_FLAG);
      fwrite((char *)allocate_batches->write_exact_sp->buffer_p, allocate_batches->write_exact_sp->size, 1, allocate_batches->fd_exact);
      fwrite((char *)allocate_batches->write_extend_sp->buffer_p, allocate_batches->write_extend_sp->size, 1, allocate_batches->fd_extend);
      allocate_batches->write_exact_sp->size = 0;
      allocate_batches->write_extend_sp->size = 0;
    } 
    /*
    if(( allocate_batches->write_extend_sp->size + 100) > write_size) {
      //item_p = list_item_new(0, WRITE_ITEM,  allocate_batches->write_extend_sp);
      //list_insert_item(item_p, write_list_p);
      //allocate_batches->write_extend_sp = write_batch_new(write_size, SPLICE_EXTEND_FLAG);
      } */
    
    if (node->allocate_ends[i]->reads_number) {
      bytes_exact = pack_junction(chromosome, node->allocate_ends[i]->strand, 
				  node->splice_start, node->allocate_ends[i]->end, 
				  junction_id, node->allocate_ends[i]->reads_number, 
				  &(((char *)allocate_batches->write_exact_sp->buffer_p)[allocate_batches->write_exact_sp->size]));
      
      bytes_extend = pack_junction(chromosome, node->allocate_ends[i]->strand, node->splice_start_extend, 
				   node->allocate_ends[i]->splice_end_extend, junction_id, node->allocate_ends[i]->reads_number, 
				   &(((char *)allocate_batches->write_extend_sp->buffer_p)[allocate_batches->write_extend_sp->size])); 
      
      allocate_batches->write_exact_sp->size += bytes_exact;
      allocate_batches->write_extend_sp->size += bytes_extend;
      
      total_splice += node->allocate_ends[i]->reads_number;
      junction_id++;
    }
  }
  return allocate_batches;
  //return exact_splice_write_p;
}


void write_chromosome_avls(allocate_splice_elements_t *chromosome_avls, 
			   list_t* write_list_p, char *extend_sp, char *exact_sp, 
			   unsigned int write_size, size_t nchromosomes) {
  int c, chr;
  allocate_buffers_t *allocate_batches = (allocate_buffers_t *)malloc(sizeof(allocate_buffers_t));
  //write_batch_t *exact_splice_write_p;
  //write_batch_t *extend_splice_write_p;
  //printf("Open splice\n");
  FILE *fd_exact = fopen(exact_sp, "w");
  FILE *fd_extend = fopen(extend_sp, "w");
  //  printf("Open splice\n");

  allocate_batches->fd_exact = fd_exact;
  allocate_batches->fd_extend = fd_extend;

  allocate_batches->write_exact_sp  = write_batch_new(write_size, SPLICE_EXACT_FLAG);
  allocate_batches->write_extend_sp  = write_batch_new(write_size, SPLICE_EXTEND_FLAG);
      
  for(c = 0; c < nchromosomes; c++){
    //printf("Chromosome %i:\n", c);
    if(chromosome_avls[c].avl_splice->root != NULL) {
      chr = c + 1;
      //printf("\tYes\n");
      //allocate_batches->write_extend_sp  = write_batch_new(1000, SPLICE_EXTEND_FLAG);

      allocate_batches = process_avlnode_in_order(chromosome_avls[c].avl_splice->root, chr, write_list_p, write_size, allocate_batches);
      
      //exact_splice_write_p = allocate_batches->write_exact_sp;
      //extend_splice_write_p = allocate_batches->write_extend_sp;
      
      if(allocate_batches->write_extend_sp != NULL) {
	if(allocate_batches->write_extend_sp->size > 0) {
	  //item_p = list_item_new(0, WRITE_ITEM, exact_splice_write_p);
	  //list_insert_item(item_p, write_list_p);
	  //} else {
	  //write_batch_free(exact_splice_write_p);
	  //}
	  fwrite((char *)allocate_batches->write_exact_sp->buffer_p, allocate_batches->write_exact_sp->size, 1, allocate_batches->fd_exact);
	  fwrite((char *)allocate_batches->write_extend_sp->buffer_p, allocate_batches->write_extend_sp->size, 1, allocate_batches->fd_extend);
	  allocate_batches->write_exact_sp->size = 0;
	  allocate_batches->write_extend_sp->size = 0;
	}
      }
	/*
      if(extend_splice_write_p != NULL) {
	list_item_t* item_p = NULL;
	if(extend_splice_write_p->size > 0) {
	  item_p = list_item_new(0, WRITE_ITEM, extend_splice_write_p);
	  list_insert_item(item_p, write_list_p);
	  } else {
	  write_batch_free(extend_splice_write_p);
	}
	}*/
      
    } //end IF chromosome splice not NULL
    cp_avltree_destroy(chromosome_avls[c].avl_splice);
  }

  write_batch_free(allocate_batches->write_exact_sp);
  write_batch_free(allocate_batches->write_extend_sp);
  fclose(allocate_batches->fd_extend);
  fclose(allocate_batches->fd_exact);

  free(allocate_batches);
  basic_statistics_sp_init(total_splice, junction_id, basic_st);
  /*
  if (statistics_on) { 
    statistics_set(TOTAL_ST, 3, total_splice, statistics_p);
  }
  */
  
}

//====================================================================================================

void cp_avlnode_print_new(cp_avlnode *node, int level){
  int i;
  if (node->right) cp_avlnode_print_new(node->right, level + 1);
  for (i = 0; i < level; i++) printf("  . ");
  printf("(%d) [%lu => %lu]", node->balance, ((node_element_splice_t *)node->key)->splice_start, ((node_element_splice_t *)node->value)->splice_start);
  node_list_print((node_element_splice_t *)node->value);
  if (node->left) cp_avlnode_print_new(node->left, level + 1);
}

void cp_avlnode_print_in_order(cp_avlnode *node){
    
  if (node->left) cp_avlnode_print_in_order(node->left);
  printf("(%d) [%lu => %lu]", node->balance, ((node_element_splice_t *)node->key)->splice_start, ((node_element_splice_t *)node->value)->splice_start);
  node_list_print((node_element_splice_t *)node->value);
  if (node->right) cp_avlnode_print_in_order(node->right);
}

void node_list_print(node_element_splice_t *node){
  int i;
  printf("::(%lu)Ends{", node->number_allocate_ends);
  for(i = 0; i < node->number_allocate_ends; i++)
    printf("|%lu-%lu|#", node->allocate_ends[i]->end, node->allocate_ends[i]->reads_number);
    printf("}\n");
}
