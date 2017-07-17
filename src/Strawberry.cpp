#include <iostream>
#include <algorithm>
#include <getopt.h>
#include <chrono>
#include <fstream>

#include "fasta.h"
#include "gff.h"
#include "alignments.h"
#include "logger.hpp"
#include "StrawberryConfig.hpp"
#include "interval.hpp"
#include "isoform.h"
#include "kmer.h"

//#include "kmer.h"
//#include "qp.h"


//#define OPT_MIN_DEPTH_4_ASSEMBLY 260
#define OPT_MIN_DEPTH_4_TRANSCRIPT     261
#define OPT_MIN_SUPPORT_4_INTRON   262
#define OPT_ALLOW_MULTIPLE_HITS   263
#define OPT_MIN_EXON_COV   264
//#define OPT_NO_ASSEMBLY 260
using namespace std;


static struct option long_options[] = {
      {"output-dir",                      required_argument,      0,       'o'},
      {"verbose",                         no_argument,            0,       'v'},
      {"max-insert-size",                 required_argument,      0,       'I'},
      {"max-junction-splice-size",        required_argument,      0,       'J'},
      {"min-junction-splice-size",        required_argument,      0,       'j'},
      {"min-mapping-qual",                required_argument,      0,       'q'},
      {"num-reads-4-prerun",              required_argument,      0,       'n'},
      {"allow-multiple-his",              no_argument,            0,       OPT_ALLOW_MULTIPLE_HITS},
#if ENABLE_THREADS
      {"num-threads",                      required_argument,      0,       'p'},
#endif
//assembly
      {"GTF",                             required_argument,      0,       'g'},
      {"no-assembly",                     no_argument,            0,       'r'},
      {"min-transcript-size",             required_argument,      0,       't'},
      {"max-overlap-distance",            required_argument,      0,       'd'},
      {"small-anchor-size",               required_argument,      0,       's'},
      {"small-anchor-alpha",              required_argument,      0,       'a'},
      {"min-support-4-intron",            required_argument,      0,       OPT_MIN_SUPPORT_4_INTRON},
      //{"min-depth-4-assembly",            required_argument,      0,       OPT_MIN_DEPTH_4_ASSEMBLY},
      {"min-depth-4-transcript",              required_argument,      0,       OPT_MIN_DEPTH_4_TRANSCRIPT},
      {"combine-short-transfrag",          no_argument,            0,       'c'},
//quantification
      {"insert-size-mean-and-sd",         required_argument,      0,       'i'},
      {"bias-correction",                 required_argument,      0,       'b'},
      {"min-isoform-frac",                required_argument,      0,       'm'},
      {"fragment-context",                required_argument,      0,       'f'},
      {"filter-low-expression",           required_argument,      0,       'e'},
      {"min-exon-cov",                    required_argument,      0,       OPT_MIN_EXON_COV},
      {0, 0, 0, 0} // terminator
};

const char *short_options = "m:q:p:o:i:j:J:n:g:t:d:s:a:b:f:e:crvGc";

void print_help()
{
   fprintf(stderr, "\nstrawberry v%s\n", strawberry::version);
   fprintf(stderr, "--------------------------------------\n");
   fprintf(stderr, "Usage: strawberry [options] <input.bam> \n");
   fprintf(stderr, "General Options:\n");
   fprintf(stderr, "   -o/--output-dir                       Output files directory.                                                                              [default:     ./strawberry_out ]\n");
   fprintf(stderr, "   -g/--GTF                              Reference transcripts annotation file. Current support gff3 and gtf format.                          [default:     NULL]\n");
   fprintf(stderr, "   -r/--no-assembly                      Skip assembly and use reference annotation to quantify transcript abundance (only use with -g)       [default:     false]\n");
   fprintf(stderr, "   -p/--num-threads                      number of threads used for Strawberry                                                                [default:     1]\n");
   fprintf(stderr, "   -v/--verbose                          Strawberry starts to gives more information.                                                         [default:     false]\n");
   fprintf(stderr, "   -q/--min-mapping-qual                 Minimum mapping quality to be included in the analyses.                                              [default:     0]\n");
   fprintf(stderr, "   -J/--max-junction-splice-size         Maximum spliced junction.                                                                            [default:     200000]\n");
   fprintf(stderr, "   -j/--min-junction-splice-size         Minimum spliced junction size.                                                                       [default:     50]\n");
   fprintf(stderr, "   -m/--min-isoform-frac                 Minimum isoform fraction.                                                                            [default:     0.05]\n");
   //fprintf(stderr, "   -n/--num-read-4-prerun                Use this number of reads to calculate empirical insert size distribution.                            [default:     500000]\n");
   fprintf(stderr, "   --allow-multiple-his                  By default, Strawberry only use reads which map to unique position in the genome.                    [default:     false]\n");
   fprintf(stderr, "\n Assembly Options:\n");
   fprintf(stderr, "   -t/--min-transcript-size              Minimun transcript size to be assembled.                                                             [default:     200]\n");
   fprintf(stderr, "   -d/--max-overlap-distance             Maximum distance between read clusters to be merged.                                                 [default:     30]\n");
   fprintf(stderr, "   -s/--small-anchor-size                Read overhang less than this value is subject to Binomial test.                                      [default:     4]\n");
   fprintf(stderr, "   -a/--small-anchor-alpha               Threshold alpha for junction binomial test filter.                                                   [default:     0]\n");
   fprintf(stderr, "   --min-support-4-intron                Minimum number of spliced aligned read required to support a intron.                                 [default:     1.0] \n");
   fprintf(stderr, "   --min-exon-cov                        Minimum exon coverage.                                                                               [default:     1.0] \n");
   fprintf(stderr, "   -c/-combine-short-transfrag           merging non-overlap short transfrags.                                                                [default:     false]\n");
//   fprintf(stderr, "   --min-depth-4-assembly                Minimum read depth for a locus to be assembled.                                                      [default:     1]\n");
   fprintf(stderr, "   --min-depth-4-transcript              Minimum average read depth for transcript.                                                           [default:     1.0]\n");

   fprintf(stderr, "\n Quantification Options:\n");
   fprintf(stderr, "   -f/--fragment-context                 Print fragment context for differential expression to this file.                                     [default:     Disabled]\n");
   fprintf(stderr, "   -i/--insert-size-mean-and-sd          User specified insert size mean and standard deviation, format: mean/sd, e.g., 300/25.               [default:     Disabled]\n");
   fprintf(stderr, "                                         This will disable empirical insert distribution learning.                                            [default:     NULL]\n");
   fprintf(stderr, "   -b/--bias-correction                  Specify reference genome for bias correction.                                                        [default:     NULL]\n");
   //fprintf(stderr, "  --infer-missing-end                Disable infering the missing end for a pair of reads.                                                [default:     true]\n" );
   fprintf(stderr, "   -e/--filter-low-expression            Skip isoforms whose relative expression (within locus) are less than this number.                    [default:     0.]\n" );
}

int parse_options(int argc, char** argv)
{
   int option_index = 0;
   int next_option;
   do {
      next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
      switch(next_option)
            {
               case  -1:
                        break;
               case 'o':
                        output_dir = optarg;
                        break;
               case OPT_ALLOW_MULTIPLE_HITS:
                        use_only_unique_hits = false;
                        break;
               case 'p':
                        num_threads = parseInt(optarg, 1, "-p/--num-threads must be at least 1", print_help);
                        use_threads = true;
                        break;
               case 'f':
                        print_frag_context = true;
                        frag_context_out = optarg;
               case 'v':
                        verbose = true;
                        break;
               case 'q':
                        kMinMapQual = parseInt(optarg, 0, "-q/--min-mapping-qual must be at least 0", print_help);
                        break;
               case 'J':
                        kMaxIntronLength = parseInt(optarg, 1, "-J/--max-intron-size must be at least 1", print_help);
                        break;
               case 'j':
                        kMinIntronLength = parseInt(optarg, 1, "-j/--min-intron-size must be at least 1", print_help);
                        break;
               case 'n':
                        kMaxReadNum4FD = parseInt(optarg, 50000, "-n/--num-read-4-prerun is suggested to be at least 50000", print_help);
                        break;
               case 'g':
                        ref_gtf_filename = optarg;
                        utilize_ref_models = true;
                        break;
               case 'r':
                        no_assembly = true;
                        enforce_ref_models = true;
                        break;
               case 'e':
                        kMinIsoformFrac = parseFloat(optarg, 0, 1.0, "-e/--filter-low-expression must be between 0-1.0", print_help);
                        filter_by_expression = true;
                        break;
               case 't':
                        kMinTransLen = parseInt(optarg, 1, "-t/--min-trancript-size must be at least 1", print_help);
                        break;
               case 'd':
                        kMaxOlapDist = parseInt(optarg, 1, "-d/--max-overlap-distance must be at least 1", print_help);
                        break;
               case 's':
                        kMaxSmallAnchor = parseInt(optarg, 1, "-s/--small-anchor-size must be at least 1", print_help);
                        break;
               case 'a':
                        kBinomialOverHangAlpha = parseFloat(optarg, 0, 1.0, "-a/--small-anchor-alpha must be between 0-1.0", print_help);
                        break;
               case OPT_MIN_SUPPORT_4_INTRON:
                        kMinJuncSupport = parseInt(optarg, 1, "--min-support-4-intron must be at least 1", print_help);
                        break;
               case OPT_MIN_EXON_COV:
                        kMinExonDoc= parseFloat(optarg, 0, 999999.0, "--min-exon-cov must be at least 0", print_help);
                        break;
               case OPT_MIN_DEPTH_4_TRANSCRIPT:
                        kMinDepth4Contig = parseFloat(optarg, 0.1, 999999.0, "--min-depth-4-quant must be at least 0.1", print_help);
                        break;
               case 'c':
                        kCombineShrotTransfrag = true;
                        break;
               case 'm':
                        kMinIsoformFrac = parseFloat(optarg, 0.001, 999999.0, "--min-isoform-frac must be at least 0.001", print_help);
                        break;
               case 'b':
                        ref_fasta_file = optarg;
                        BIAS_CORRECTION = true;
                        break;
               case 'i':
                       {
                          vector<string> mean_and_sd;
                          split(optarg, "/", mean_and_sd);
                          if(mean_and_sd.size() != 2 ) {
                             cerr<<"Wrong format for specifying insert size mean and sd!"<<endl;
                             print_help();
                             exit(1);
                          }
                          cerr<<"Insert size mean: "<<mean_and_sd[0]<<endl;
                          cerr<<"Insert size sd: "<<mean_and_sd[1]<<endl;
                          kInsertSizeMean = (double) parseInt(mean_and_sd[0].c_str(), 1, "minimun insert size mean must be at least 1", print_help);
                          kInsertSizeSD = (double) parseInt(mean_and_sd[1].c_str(), 1, "minimun insert size sd must be at least 1", print_help);
                          break;
                       }
               default:
                        print_help();
                        return 1;
            }
   }while(next_option != -1);

   return 0;
}



int driver(const char* const bam_file, FILE* pFile, FILE* plogfile, FILE* pfragfile){
   auto start = chrono::steady_clock::now();
   ReadTable read_table;
   RefSeqTable ref_seq_table(true);
   shared_ptr<HitFactory> hf(new BAMHitFactory(bam_file, read_table, ref_seq_table));
   hf->inspect_header();
   Sample read_sample(move(hf));

   GffReader* greader= NULL;
   if(ref_gtf_filename != ""){
      FILE* gff = fopen(ref_gtf_filename.c_str(), "r");
      if(gff == NULL){
         fprintf(stderr, "Error: cannot open refernce gtf file %s for reading\n", ref_gtf_filename.c_str());
         exit(1);
      }
      greader = new GffReader(ref_gtf_filename.c_str(), gff);
      greader->readAll(); fclose(gff);
      // count 2-isoform genes
//      int num_iso2gene = 0;
//      for (auto const& each : greader->_g_seqs) {
//         for(auto const& gene: each->_genes) {
//            cout<<gene->_mrnas.size()<<endl;
//            if (gene->_mrnas.size() == 2) {
//               num_iso2gene++;
//            }
//         }
//      }
//      std::cerr<<"number of two-isoform genes: "<<num_iso2gene<<std::endl;
//      exit(0);
      greader->sortExonOrderInMinusStrand();
      read_sample.loadRefmRNAs(greader->_g_seqs, ref_seq_table);
   }


   if(ref_fasta_file != "") {
      //if(ref_gtf_filename != ""){
         //read_sample.loadRefmRNAs(greader->_g_seqs, ref_seq_table);
      //}
      //else{
         shared_ptr<FaInterface> fa_api(new FaInterface());
         fa_api->initiate(ref_fasta_file.c_str());
         shared_ptr<FaSeqGetter> fsg(new FaSeqGetter());
         read_sample._fasta_interface = move(fa_api);
         read_sample._fasta_getter = move(fsg);
      //}
   }
   if(verbose){
      cerr<<"Inspecting sample......"<<endl;
   }
   if (no_assembly) read_sample.preProcess(plogfile);
   else read_sample.assembleSample(plogfile);

   if(verbose){
      cerr<<"Total number of mapped reads is: "<<read_sample.total_mapped_reads()<<endl;
   }


   if(SINGLE_END_EXP){
      kInsertSizeMean = 200;
      kInsertSizeSD = 80;
      infer_the_other_end = false;
   }

   if(kInsertSizeMean !=0 && kInsertSizeSD != 0){
      if(verbose){
         cerr<<"Using user specified insert size mean: "<<kInsertSizeMean<<" and standard deviation: "<<kInsertSizeSD<<endl;
      }
      unique_ptr<InsertSize> insert_size(new InsertSize(kInsertSizeMean, kInsertSizeSD));
      read_sample._insert_size_dist = move(insert_size);
   }
   else{
      const vector<int> & fd = read_sample._hit_factory->_reads_table._frag_dist;
      unique_ptr<InsertSize> insert_size(new InsertSize(fd));
      if(verbose){
         cerr<<"Using empirical insert size distribution "<<endl;
         //cerr<<"Mean insert size mean: "<<insert_size->_mean<<endl;
         //cerr<<"Mean insert size standard deviation: "<<insert_size->_sd<<endl;
      }
      read_sample._insert_size_dist = move(insert_size);
   }


   read_sample.procSample(pFile, plogfile, pfragfile);

   fclose(pFile);
   fclose(plogfile);
   if (pfragfile != NULL) fclose(pfragfile);
   pFile = NULL;
   plogfile = NULL;
   pfragfile = NULL;
   delete greader;
   greader = NULL;
   auto end = chrono::steady_clock::now();
   auto diff = end - start;
   cerr << "Finished in " << chrono::duration <double, milli> (diff).count() << " ms" << endl;
   return 0;
}


INITIALIZE_EASYLOGGINGPP
int main(int argc, char** argv){
   string cmdline;
   for(int i=0; i<argc; i++){
      cmdline += argv[i];
      cmdline+=" ";
   }
   int parse_ret = parse_options(argc, argv);
   if(parse_ret) return 1;
   if(optind >= argc){
      print_help();
      return 1;
   }
   int ret = mkpath(output_dir.c_str(), 0777);
   if(ret == -1){
      if(errno != EEXIST){
         fprintf(stderr, "ERROR: cannot create directory %s\n", output_dir.c_str());
         exit(1);
      }
   }
   string exec_log_file = output_dir + "/" + "execution.log";
   CreateLogger(exec_log_file);
   string assembled_file = output_dir + "/";// assembled_transcripts.gtf 25 characters
   assembled_file += string("assembled_transcripts.gtf");
   string tracker = output_dir + "/" + tracking_log;
   string fragfile = output_dir + "/" + frag_context_out;

   if(verbose){
      fprintf(stderr, "OUTPUT gtf file: \n%s\n", assembled_file.c_str());
      fprintf(stderr, "see %s for the progress of the program \n", tracker.c_str());
   }
   FILE *pFile = fopen(assembled_file.c_str(), "w");
   fprintf(pFile, "#%s\n", cmdline.c_str());
   fprintf(pFile, "#########################################\n");
   FILE *plogfile = fopen(tracker.c_str(), "w");

   FILE* pfragfile = NULL;
   if (print_frag_context) pfragfile = fopen(fragfile.c_str(), "w");

   char* bam_file = argv[optind++];
   driver(bam_file, pFile, plogfile, pfragfile);
   fprintf(stdout, "Program finished\n");
   return 0;
}
