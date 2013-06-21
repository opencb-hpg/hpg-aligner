#ifndef BREAKPOINT_H
#define BREAKPOINT_H

#include "containers/array_list.h"

#include "bioformats/bam/alignment.h"

//--------------------------------------------------------------------------------------

//====================================================================================
//  Input structure for CIGAR format
//====================================================================================

typedef struct cigar_op {
  int number;
  //  char number_str[16];
  char name;
} cigar_op_t;
//cigar_op_t *cigar_op_new(int number, char number_str[], char name);
cigar_op_t *cigar_op_new(int number, char name);
void cigar_op_free(cigar_op_t *cigar_op);

typedef struct cigar_code {
  //  int num_ops;
  //  int num_allocated_ops;
  //  cigar_op_t *ops;
  char *cigar_str;
  array_list_t *ops;
} cigar_code_t;

cigar_code_t *cigar_code_new();
cigar_code_t *cigar_code_new_by_string(char *cigar_str);
void cigar_code_free(cigar_code_t* p);

char *cigar_code_get_string(cigar_code_t *p);
int cigar_code_get_num_ops(cigar_code_t *p);
cigar_op_t *cigar_code_get_last_op(cigar_code_t *p);
void cigar_code_append_op(cigar_op_t *op, cigar_code_t *p);

cigar_code_t *generate_cigar_code(char *query_map, char *ref_map, unsigned int map_len,
				  unsigned int query_start, unsigned int query_len, 
				  int *distance);

//====================================================================================
//  Input structure for breakpoint definition
//====================================================================================

#define LEFT_SIDE  0
#define RIGHT_SIDE 1

typedef struct breakpoint_info {
  int side;
  int num_M;
  int index_M;
  cigar_code_t *cigar;
  //  char *seq;
} breakpoint_info_t;

breakpoint_info_t *breakpoint_info_new(int side, int num_M, int index_M,
				       cigar_code_t *cigar);
void breakpoint_info_free(breakpoint_info_t *p);

//--------------------------------------------------------------------------------------

typedef struct breakpoint {
  int chr_index;
  int position;
  int coverage;
  array_list_t *info_list;
} breakpoint_t;

breakpoint_t *breakpoint_new(int chr_index, int position);
void breakpoint_free(breakpoint_t *p);

breakpoint_t *compute_breakpoint(int chr_index, int alig_position, 
				 cigar_code_t *cigar, array_list_t *list);

//breakpoint_t *get_breakpoint(int chr_index, size_t pos, array_list_t *list);
void display_breakpoints(array_list_t *list);

//====================================================================================
//  Input structure for CIGAR format
//====================================================================================

extern array_list_t *breakpoint_list;

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

#endif  // SW_SERVER_H
