#ifndef OPTIONS_H
#define OPTIONS_H

/*
 * hpg-fastq.h
 *
 *  Created on: Aug 29, 2012
 *      Author: imedina
 */

#include <stdlib.h>
#include <string.h>

#include "argtable2.h"
#include "libconfig.h"
#include "commons/log.h"
#include "commons/system_utils.h"
#include "commons/file_utils.h"

//============================ DEFAULT VALUES ============================
#define DEFAULT_GPU_THREADS		32
#define DEFAULT_CPU_THREADS		1
#define DEFAULT_CAL_SEEKER_ERRORS	0
#define DEFAULT_MIN_SEED_PADDING_LEFT	5
#define DEFAULT_MIN_SEED_PADDING_RIGHT	5
#define DEFAULT_WRITE_BATCH_SIZE	500000
#define DEFAULT_NUM_CAL_SEEKERS		1
#define DEFAULT_REGION_THREADS		1
#define DEFAULT_NUM_SW_THREADS		1
#define DEFAULT_MIN_NUM_SEEDS		10
#define DEFAULT_MAX_NUM_SEEDS		20
#define DEFAULT_MAX_INTRON_LENGTH	800000
#define DEFAULT_MIN_INTRON_LENGTH	40
#define DEFAULT_SW_MIN_SCORE		0.6
#define DEFAULT_SW_MATCH		5
#define DEFAULT_SW_MISMATCH		-4
#define DEFAULT_SW_GAP_OPEN		10
#define DEFAULT_SW_GAP_EXTEND		0.5
#define DEFAULT_PAIR_MODE	        0
#define DEFAULT_PAIR_MIN_DISTANCE	0
#define DEFAULT_PAIR_MAX_DISTANCE	800
#define MINIMUM_BATCH_SIZE              10000
//=====================================================================

#define NUM_OPTIONS			42

typedef struct options {
  char mode[64];
  unsigned char bwt_set;
  unsigned char reg_set;
  unsigned char cal_set;
  unsigned char sw_set;
  int min_intron_length;
  int num_gpu_threads;
  int num_cpu_threads;
  int min_cal_size; 
  int min_num_seeds; 
  int max_num_seeds; 
  int seeds_max_distance;
  int min_seed_padding_right;
  int min_seed_padding_left;
  int batch_size;
  int write_size;
  int num_cal_seekers;
  int min_seed_size;
  int seed_size;
  int max_intron_length;
  int flank_length;
  int timming;
  int statistics;
  int rna_seq; 
  int help;
  int cal_seeker_errors;
  int pair_mode;
  int pair_min_distance;
  int pair_max_distance;
  int report_all;
  int report_best;
  int report_n_hits;
  int report_unpaired;
  int gpu_process;
  int log_level;
  int index_ratio;
  double min_score;
  double match;
  double mismatch;
  double gap_open;
  double gap_extend;
  char* prefix_name;
  char* in_filename;
  char* in_filename2;
  char* bwt_dirname;
  char* genome_filename;
  char* output_name;
  char* header_filename;
  char* intron_filename;
} options_t;


options_t *options_new(void);

void options_free(options_t *options);

void options_display(options_t *options);

void usage_cli();

void validate_options(options_t *options, char *mode);

/**
 * @brief Initializes an global_options_t structure mandatory members.
 * @return A new global_options_t structure.
 *
 * Initializes the only mandatory member of a global_options_t, which is the output directory.
 */
void** argtable_options_new(void);


/**
 * @brief Free memory associated to a global_options_data_t structure.
 * @param options_data the structure to be freed
 *
 * Free memory associated to a global_options_data_t structure, including its text buffers.
 */
void argtable_options_free(void **argtable_options);


/**
 * @brief Initializes an global_options_data_t structure mandatory members.
 * @return A new global_options_data_t structure.
 *
 * Initializes the only mandatory member of a global_options_data_t, which is the output directory.
 */
options_t *read_CLI_options(void **argtable_options, options_t *options);



/**
 * @brief Reads the configuration parameters of the application.
 * @param filename file the options data are read from
 * @param options_data options values (vcf filename, output directory...)
 * @return Zero if the configuration has been successfully read, non-zero otherwise
 *
 * Reads the basic configuration parameters of the application. If the configuration file can't be
 * read, these parameters should be provided via the command-line interface.
 */
int read_config_file(const char *filename, options_t *options);



options_t *parse_options(int argc, char **argv);


void usage(void **argtable);

#endif

//typedef struct argtable_options {
//	/*	IO options	*/
//	struct arg_file *fastq_file; /**< VCF file used as input. */
//	struct arg_file *fastq1_file; /**< PED file used as input. */
//	struct arg_file *fastq2_file; /**< PED file used as input. */
//	struct arg_file *config_file; /**< Path to the configuration file */
//	struct arg_file *genomic_signature_input; /**< Filename template for the main output file. */
//	struct arg_str *output_directory; /**< Directory where the output files will be stored. */
//
//	//    struct arg_int qc_flag;
//	//    struct arg_int *filter_flag;
//	//    struct arg_int *prepro_flag;
//	struct arg_int *min_read_length;	// 50
//	struct arg_int *max_read_length;	// 200
//	struct arg_int *min_quality;  	// 20
//	struct arg_int *max_quality;  	// 60
//	struct arg_int *max_nts_out_quality;
//	struct arg_int *max_n_per_read;
//	struct arg_int *start_quality_nt;
//	struct arg_int *end_quality_nt;
//	struct arg_int *phred_quality;
//	struct arg_int *rtrim_nts;
//	struct arg_int *ltrim_nts;
//	struct arg_int *rfilter_nts;
//	struct arg_int *lfilter_nts;
//	struct arg_lit *kmers_flag;
//
//	struct arg_lit *cg_flag; //Chaos Game flag
//	struct arg_int *k_cg;    //7
//
//	// variables to store the value options
//	//    struct arg_int *gpu_num_blocks =  DEFAULT_GPU_NUM_BLOCKS;   // 16
//	//    struct arg_int *gpu_num_threads =  DEFAULT_GPU_NUM_THREADS;  // 512
//	//    struct arg_int *gpu_num_devices =  DEFAULT_GPU_NUM_DEVICES;  // -1
//	struct arg_int *cpu_num_threads;  			// 2
//	//    struct arg_int *cpu_qc_calc_num_threads; 	//0
//	struct arg_int *batch_size; 				// 64MB
//	struct arg_int *batch_list_size;  			// 4
//
//	struct arg_int *log_level; /**< desc */
//	struct arg_file *log_file; /**< desc */
//	struct arg_lit *verbose; /**< desc */
//	struct arg_lit *help; /**< desc */
//	struct arg_lit *time; /**< desc */
//
//	int num_options;
//} argtable_options_t;
