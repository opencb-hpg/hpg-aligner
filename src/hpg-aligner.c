#include "dna/dna_aligner.h"
#include "rna/rna_aligner.h"
#include "build-index/index_builder.h"


double emboss_matrix_t = 0.0f, emboss_tracking_t = 0.0f;
double sse_matrix_t = 0.0f, sse_tracking_t = 0.0f;
double sse1_matrix_t = 0.0f, sse1_tracking_t = 0.0f;
double avx_matrix_t = 0.0f, avx_tracking_t = 0.0f;
double avx1_matrix_t = 0.0f, avx1_tracking_t = 0.0f;

//--------------------------------------------------------------------
// constants
//--------------------------------------------------------------------

#define OPTIONS 			30
#define MIN_ARGC  			5
#define NUM_SECTIONS_STATISTICS 	5
#define NUM_SECTIONS_STATISTICS_SB	21
//#define REQUIRED 1
//#define NO_REQUIRED 0

//--------------------------------------------------------------------
// global variables for log functions
//--------------------------------------------------------------------

//int log_level = DEFAULT_LOG_LEVEL;
//bool verbose = true;

//--------------------------------------------------------------------
// global variables for timing and capacity meausures
//--------------------------------------------------------------------

char time_on = 0;
char statistics_on = 0;
timing_t* timing = NULL;
basic_statistics_t *basic_st;
cal_st_t cal_st;
double kl_time;

pthread_cond_t cond_sp;
pthread_mutex_t mutex_sp;

size_t bwt_correct = 0;
size_t bwt_error = 0;
size_t seeding_reads = 0;
pthread_mutex_t bwt_mutex;

// timing
double main_time;
size_t TOTAL_SW, 
  TOTAL_READS_PROCESS,
  TOTAL_READS_SEEDING,
  TOTAL_READS_SEEDING2,
  TOTAL_READS_SA;

struct timeval time_start_alig, time_end_alig;
double time_alig;

size_t w2_r = 0;
size_t w3_r = 0;
size_t w2_3_r = 0;

size_t tot_reads_in = 0;
size_t tot_reads_out = 0;

FILE *fd_log;
size_t junction_id;

float min_score = 300.0, tot_score;
size_t tot_rep;

/*
pthread_mutex_t sw_mutex;
size_t *histogram_sw;
size_t num_reads_map;
size_t num_reads;
double seeding_time_2;
size_t reads_cals, reads_single, reads_h;
size_t tot_reads;
*/

/*
void run_dna_aligner(genome_t *genome, bwt_index_t *bwt_index, 
		     bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
		     pair_mng_t *pair_mng, report_optarg_t *report_optarg,
		     options_t *options);


void run_rna_aligner(genome_t *genome, bwt_index_t *bwt_index, pair_mng_t *pair_mng,
		     bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg,
		     report_optarg_t *report_optarg, options_t *options);
*/

//--------------------------------------------------------------------
// main parameters support
//--------------------------------------------------------------------
int main(int argc, char* argv[]) {

  pthread_mutex_init(&mutex_sp, NULL);
  TOTAL_SW = 0;
  TOTAL_READS_PROCESS = 0;
  TOTAL_READS_SEEDING = 0;
  TOTAL_READS_SEEDING2 = 0;
  TOTAL_READS_SA = 0;
  
  const char HEADER_FILE[1024] = "Human_NCBI37.hbam\0";
  basic_st = basic_statistics_new();

  // init logs, after parsing the command-line
  // logs will be re-set according to the command-line
  log_level = LOG_FATAL_LEVEL;
  log_verbose = 1;
  log_file = NULL;

  if (argc <= 1) {
    LOG_FATAL("Missing command.\nValid commands are:\n\tdna: to map DNA sequences\n\trna: to map RNA sequences\n\tbuild-index: to create the genome index.\nUse -h or --help to display hpg-aligner options.");
  }

  if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
    usage_cli();
  }

  char *command = argv[1];  
  // We need to consume command: {dna | rna | bs | build-index}
  argc -= 1;
  argv += 1;

  if(strcmp(command, "dna") != 0 && 
     strcmp(command, "rna") != 0 &&
     //     strcmp(command, "bs" ) != 0 && 
     strcmp(command, "build-index") != 0) {
    LOG_FATAL("Command unknown.\nValid commands are:\n\tdna: to map DNA sequences\n\trna: to map RNA sequences\n\tbuild-index: to create the genome index.\nUse -h or --help to display hpg-aligner options.");
  }

  // parsing options
  options_t *options = parse_options(argc, argv);

  // now, we can set logs according to the command-line
  init_log_custom(options->log_level, 1, "hpg-aligner.log", "w");

  validate_options(options, command);
  LOG_DEBUG_F("Command Mode: %s\n", command);

  if (!strcmp(command, "build-index")) {
       run_index_builder(options->genome_filename, options->bwt_dirname, options->index_ratio);
       LOG_DEBUG("Done !!\n");
       exit(0);
  }

  time_on =  (unsigned int) options->timming;
  statistics_on =  (unsigned int) options->statistics;

  struct timeval time_genome_s, time_genome_e;
  double time_genome;

  start_timer(time_genome_s);
  // Genome parameters
  LOG_DEBUG("Reading genome...");
  genome_t* genome = genome_new("dna_compression.bin", options->bwt_dirname);
  LOG_DEBUG("Done !!");
  
  // Metaexons structure
  metaexons_t *metaexons = metaexons_new(genome);

  // BWT index
  LOG_DEBUG("Reading bwt index...");
  
  //TODO: CLI parameters, bwt_index_new(bwt_index, inverse_sa (more than one error), duplicate_strand (search strand + & - ))
  bwt_index_t *bwt_index = bwt_index_new(options->bwt_dirname, false, false);
  
  stop_timer(time_genome_s, time_genome_e, time_genome);
  LOG_DEBUG("Reading bwt index done\n");
  
  //BWT parameters
  bwt_optarg_t *bwt_optarg = bwt_optarg_new(1, 0,
					    options->filter_read_mappings, 
					    options->filter_seed_mappings);
  
  // CAL parameters
  cal_optarg_t *cal_optarg = cal_optarg_new(options->min_cal_size, options->seeds_max_distance, 
					    options->num_seeds, options->min_num_seeds_in_cal,
					    options->seed_size, options->min_seed_size, 
					    options->cal_seeker_errors);
  
  // paired mode parameters
  pair_mng_t *pair_mng = pair_mng_new(options->pair_mode, options->pair_min_distance, 
				      options->pair_max_distance, options->report_only_paired);
  
  // report parameters
  report_optarg_t *report_optarg = report_optarg_new(options->report_all,
						     options->report_n_best,
						     options->report_n_hits, 
						     options->report_only_paired,
						     options->report_best);  
  LOG_DEBUG("init table...");
  initTable();
  LOG_DEBUG("init table done !!");

  if (!strcmp(command, "rna")) {
    run_rna_aligner(genome, bwt_index, pair_mng, bwt_optarg, cal_optarg, report_optarg, metaexons, options);
  } else {
    // DNA version
    run_dna_aligner(genome, bwt_index, bwt_optarg, cal_optarg, pair_mng, report_optarg, metaexons, options);
  }

  LOG_DEBUG("main done !!");

  //show_metaexons(metaexons);
  // Free memory
  
  bwt_index_free(bwt_index);
  genome_free(genome);
  bwt_optarg_free(bwt_optarg);
  cal_optarg_free(cal_optarg);
  pair_mng_free(pair_mng);
  report_optarg_free(report_optarg);

  if (time_on) { 
    //timing_display(timing);
  }


  basic_statistics_display(basic_st, !strcmp(command, "rna"), time_alig / 1000000, time_genome / 1000000);


  if (!strcmp(command, "rna") && fd_log != NULL) {
    char str_log[2048];
    if (options->in_filename) {
      sprintf(str_log, "FILE INPUT 0: %s\n\0", options->in_filename);
      fwrite(str_log, sizeof(char), strlen(str_log), fd_log);
    }

    if (options->in_filename2) {
      sprintf(str_log, "FILE INPUT 1: %s\n\0", options->in_filename2);
      fwrite(str_log, sizeof(char), strlen(str_log), fd_log);
    }

    sprintf(str_log, "TOTAL READS PROCESSED: %i\n\0", basic_st->total_reads);
    fwrite(str_log, sizeof(char), strlen(str_log), fd_log);

    sprintf(str_log, "    TOTAL READS MAPPED: %0.2f%%\n\0", basic_st->num_mapped_reads * 100.0 / basic_st->total_reads);
    fwrite(str_log, sizeof(char), strlen(str_log), fd_log);

    sprintf(str_log, "    TOTAL READS UNMAPPED: %0.2f%%\n\0", (basic_st->total_reads - basic_st->num_mapped_reads) * 100.0 / basic_st->total_reads);
    fwrite(str_log, sizeof(char), strlen(str_log), fd_log);

    sprintf(str_log, "    TOTAL SPLICE JUNCTIONS FOUND: %i\n\0", junction_id);
    fwrite(str_log, sizeof(char), strlen(str_log), fd_log);

    sprintf(str_log, "TIME ALIGNER: %0.2f (s)\n\0", time_alig / 1000000);
    fwrite(str_log, sizeof(char), strlen(str_log), fd_log);

    sprintf(str_log, "TOTAL TIME: %0.2f (s)\n\0", (time_alig + time_genome) / 1000000);
    fwrite(str_log, sizeof(char), strlen(str_log), fd_log);

    fclose(fd_log);

  }

  if (time_on){ timing_free(timing); }

  options_free(options);
  
  metaexons_free(metaexons);
  //printf("Reads mapped in BWT Section %lu (%lu correct(%f) and %lu misses(%f))\n", bwt_correct + bwt_error, bwt_correct, (float)(bwt_correct * 100)/(float)(bwt_correct + bwt_error), bwt_error, (float)(bwt_error * 100)/(float)(bwt_correct + bwt_error));


  //printf("Total SW Processed: %i\n", total_sw);
  //printf("Length Ref; Number\n");
  //for (int i = 0; i < 1024; i++) {
  //printf(" %i; %i\n", i, histogram_sw[i]);
  //}

  //printf("Total SW Processed: %i\n", total_sw);
  /*
  printf("Length Ref; Number\n");
  for (int i = 0; i < 1024; i++) {
    printf(" %i; %i\n", i, histogram_sw[i]);
  }
  */
  //printf("TOTAL READS SEEDING %i\n", seeding_reads);
  //printf("Time seeding second phase %f(s)\n", seeding_time_2/1000000);

  return 0;
}

//--------------------------------------------------------------------
