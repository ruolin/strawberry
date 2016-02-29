/*
>HEADER
    Copyright (c) 2015 Ruolin Liu rliu0606@gmail.com
    This file is part of Strawberry.
    Strawberry is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Strawberry is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Strawberry.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
*/
#ifdef DEBUG
   #include <iostream>
#endif

#include <algorithm>
#include <getopt.h>
#include <chrono>

#include "fasta.h"
#include "gff.h"
#include "alignments.h"
#include "logger.hpp"
#include "StrawberryConfig.hpp"
//#include "qp.h"


#define OPT_MIN_DEPTH_4_ASSEMBLY 260
#define OPT_MIN_DEPTH_4_CONTIG     261
#define OPT_MIN_SUPPORT_4_INTRON   262
#define OPT_ALLOW_MULTIPLE_HITS   263
using namespace std;


static struct option long_options[] = {
      {"output-dir",                      required_argument,      0,       'o'},
      {"verbose",                         no_argument,            0,       'v'},
      {"max-insert-size",                 required_argument,      0,       'I'},
      {"max-junction-splice-size",        required_argument,      0,       'J'},
      {"min-junction-splice-size",        required_argument,      0,       'j'},
      {"num-reads-4-prerun",              required_argument,      0,       'n'},
      {"allow-multiple-his",              no_argument,            0,       OPT_ALLOW_MULTIPLE_HITS},
#if ENABLE_THREADS
      {"num-threads",                      required_argument,      0,       'p'},
#endif
//assembly
      {"GTF",                             required_argument,      0,       'g'},
//      {"enforce-ref-model",               no_argument,            0,       'G'},
      {"min-transcript-size",             required_argument,      0,       't'},
      {"max-overlap-distance",            required_argument,      0,       'd'},
      {"small-anchor-size",               required_argument,      0,       's'},
      {"small-anchor-alpha",              required_argument,      0,       'a'},
      {"min-support-4-intron",            required_argument,      0,       OPT_MIN_SUPPORT_4_INTRON},
      {"min-depth-4-assembly",            required_argument,      0,       OPT_MIN_DEPTH_4_ASSEMBLY},
      {"combine-short-transfrag",          no_argument,            0,       'c'},
//quantification
      {"insert-size-mean-and-sd",         required_argument,      0,       'i'},
      {"bias-correction",                 required_argument,      0,       'b'},
      {"min-depth-4-contig",              required_argument,      0,       OPT_MIN_DEPTH_4_CONTIG},
      {"infer-missing-end",               no_argument,            0,       'm'},
};

#if ENABLE_THREADS
const char *short_options = "p:o:i:j:J:n:g:t:d:s:a:c:b:vGcm";
#else
const char *short_options = "o:i:j:J:n:g:t:d:s:a:c:b:vGcm";
#endif

void print_help()
{
   fprintf(stderr, "\n\nstrawberry v%s\n", strawberry::version);
   fprintf(stderr, "--------------------------------------\n");
   fprintf(stderr, "Usage: strawberry [options] <input.bam> \n");
   fprintf(stderr, "General Options:\n");
   fprintf(stderr, "   -o/--output-dir                       Output files directory.                                                                              [default:     ./out ]\n");
#if ENABLE_THREADS
   fprintf(stderr, "   -p/--num-threads                      number of threads used for Strawberry                                                                [default:     1]\n");
#endif
   fprintf(stderr, "   -v/--verbose                          Strawberry starts to gives more information.                                                         [default:     false]\n");
   fprintf(stderr, "   -J/--max-junction-splice-size         Maximum spliced junction.                                                                            [default:     200000]\n");
   fprintf(stderr, "   -j/--min-junction-splice-size         Minimum spliced junction size.                                                                       [default:     50]\n");
   fprintf(stderr, "   -n/--num-read-4-prerun                Use this number of reads to calculate empirical insert size distribution.                            [default:     500000]\n");
   fprintf(stderr, "   --allow-multiple-his                  By default, Strawberry only use reads which map to unique position in the genome.                    [default:     false]\n");
   fprintf(stderr, "\n Assembly Options:\n");
   fprintf(stderr, "   -g/--GTF                              Use reference transcripts annotation to guide assembly.                                              [default:     NULL]\n");
//   fprintf(stderr, "   -G/--enforce-ref-model                Omit assembled transcripts that are not in the reference.                    [default:     false]\n");
   fprintf(stderr, "   -t/--min-transcript-size              Minimun transcript size to be assembled.                                                             [default:     200]\n");
   fprintf(stderr, "   -d/--max-overlap-distance             Maximum distance between read clusters to be merged.                                                 [default:     30]\n");
   fprintf(stderr, "   -s/--small-anchor-size                Read overhang less than this value is subject to Binomial test.                                      [default:     4]\n");
   fprintf(stderr, "   -a/--small-anchor-alpha               Threshold alpha for junction binomial test filter.                                                   [default:     0]\n");
   fprintf(stderr, "   --min-support-4-intron                Minimum number of spliced aligned read required to support a intron.                                 [default:     1] \n");
   fprintf(stderr, "   -c/--combine-short-transfrag          Disable merging non-overlap short transfrags .                                                       [default:     true]\n");
   fprintf(stderr, "   --min-depth-4-assembly                Minimum read depth for a locus to be assembled.                                                      [default:     1]\n");

   fprintf(stderr, "\n Quantification Options:\n");
   fprintf(stderr, "   -i/--insert-size-mean-and-sd          User specified insert size mean and standard deviation, format: mean/sd, e.g., 300/25.               [default:     Disabled]\n");
   fprintf(stderr, "                                         This will disable empirical insert distribution learning.                                            [default:     NULL]\n");
   fprintf(stderr, "   -b/--bias-correction                  Use bias correction.                                                                                 [default:     false]\n");
   fprintf(stderr, "   --min-depth-4-transcript              Minimum average read delpth for transcript.                                                          [default:     1]\n");
   fprintf(stderr, "   -m/--infer-missing-end                Disable infering the missing end for a pair of reads.                                                [default:     true]\n" );
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
                        use_unique_hits = false;
                        break;
               case '-p':
                        num_threads = parseInt(optarg, 1, "-p/--num-threads must be at least 1", print_help);
                        break;
               case 'v':
                        verbose = true;
                        break;
               case 'J':
                        kMaxIntronLength = parseInt(optarg, 1, "-J/--max-intron-size must be at least 1", print_help);
                        break;
               case 'j':
                        kMinIntronLength = parseInt(optarg, 1, "-j/--min-intron-size must be at least 1", print_help);
                        break;
               case 'n':
                        kMaxReadNum4FD = parseInt(optarg, 5000, "-n/--num-read-4-prerun is suggested to be at least 5000", print_help);
                        break;
               case 'g':
                        ref_gtf_filename = optarg;
                        utilize_ref_models = true;
                        break;
//               case 'G':
//                        enforce_ref_models = true;
//                        break;
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
               case OPT_MIN_DEPTH_4_ASSEMBLY:
                        kMinDepth4Locus = parseFloat(optarg, 0, 999999.0, "--min-depth-4-assembly must be at least 0", print_help);
                        break;
               case OPT_MIN_DEPTH_4_CONTIG:
                        kMinDepth4Contig = parseFloat(optarg, 0.1, 999999.0, "--min-depth-4-quant must be at least 0.1", print_help);
                        break;
               case 'c':
                        kCombineShrotTransfrag = false;
                        break;
               case 'm':
                        infer_the_other_end = false;
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
                        //cerr<<"Invalid option"<<endl;
                        print_help();
                        return 1;
            }
   }while(next_option != -1);

   return 0;
}


int main(int argc, char** argv){
   auto start = chrono::steady_clock::now();
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

   if(output_dir != "./"){
      int ret = mkpath(output_dir.c_str(), 0777);
      if(ret == -1){
         if(errno != EEXIST){
            fprintf(stderr, "ERROR: cannot create directory %s\n", output_dir.c_str());
            exit(1);
         }
      }
   }
   string assembled_file = output_dir;// assembled_transcripts.gtf 25 characters
   assembled_file += string("/assembled_transcripts.gtf");
   string tracker = output_dir + tracking_log;
   if(verbose){
      fprintf(stderr, "OUTPUT gtf file: \n%s\n", assembled_file.c_str());
      fprintf(stderr, "see %s for the progress of the program \n", tracker.c_str());
   }
   FILE *pFile = fopen(assembled_file.c_str(), "w");
   FILE *plogfile = fopen(tracker.c_str(), "w");
   char* bam_file = argv[optind++];
   ReadTable read_table;
   RefSeqTable ref_seq_table(true);
   unique_ptr<HitFactory> hf(new BAMHitFactory(bam_file, read_table, ref_seq_table));
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
      greader->readAll();
      greader->reverseExonOrderInMinusStrand();
   }


   if(ref_fasta_file != "") {
      if(ref_gtf_filename != ""){
         read_sample.loadRefmRNAs(greader->_g_seqs, ref_seq_table, ref_fasta_file.c_str());
      }
      else{
         shared_ptr<FaInterface> fa_api(new FaInterface());
         fa_api->initiate(ref_fasta_file.c_str());
         shared_ptr<FaSeqGetter> fsg(new FaSeqGetter());
         read_sample._fasta_interface = fa_api;
         read_sample._fasta_getter = fsg;
      }
   }
   if(verbose){
      cerr<<"Inspecting sample......"<<endl;
   }
   read_sample.inspectSample(plogfile);
   if(verbose){
      cerr<<"Total number of mapped reads is: "<<read_sample.total_mapped_reads()<<endl;
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

   read_sample.procSample(pFile, plogfile );

   //const char *path = "/home/ruolin/Dropbox/Strawberry/Arabidopsis";
   //const char *ara_gtf = "/home/ruolin/Dropbox/Strawberry/TAIR10_GFF3_genes-1.gff";
   //const char *human_gtf = "/home/ruolin/Downloads/gencode.v21.annotation.gff3";
   //const char *bam_file = "/home/ruolin/Dropbox/Strawberry/accepted_hits.bam";
   //char *bam_file = "/home/ruolin/git/CompareTransAbun/assembly_comp/RD100/RD100.control_r1.concordant_uniq.sort.bam";
   //char *bam_file = "/home/ruolin/git/CompareTransAbun/assembly_comp/RD100/accepted_hits.bam";
   //const char* bam_dir = stripFileName(bam_file);
   //size_t len1 = strlen(bam_dir);

   //const char *bam_file = "/home/ruolin/Dropbox/Strawberry/WetFT1.sm.bam";

   //FaInterface fa_api(path);
   //FaSeqGetter fsg;
   //fa_api.load2FaSeqGetter(fsg, "mitochondria");
   //cout<<"success\t"<<fsg.loadSeq()<<endl;
   //cout<<fsg.fetchSeq(80,4)<<endl;
   //GffReader greader(ara_gtf);
   //greader.readAll();
   //greader.reverseExonOrderInMinusStrand();
   //ReadTable read_table;
   //RefSeqTable ref_seq_table(true);
   //unique_ptr<HitFactory> hf(new BAMHitFactory(bam_file, read_table, ref_seq_table));
   //hf->inspect_header();
   //Sample read_sample(move(hf));
   //read_sample.loadRefmRNAs(greader._g_seqs, ref_seq_table, path);

   //QpSolver qps;

   //read_sample.inspectSample();
   //double mean, sd;

   fclose(pFile);
   fclose(plogfile);
   pFile = NULL;
   plogfile = NULL;
   delete greader;
   greader = NULL;
   auto end = chrono::steady_clock::now();
   auto diff = end - start;
   cout << "Finished in " << chrono::duration <double, milli> (diff).count() << " ms" << endl;
}


