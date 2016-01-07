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

#include <iostream>
#include <algorithm>
#include "fasta.h"
#include "gff.h"
#include "alignments.h"
#include "logger.hpp"
//#include "qp.h"
#include <chrono>
using namespace std;


int main(){
   const char *path = "/home/ruolin/Dropbox/Strawberry/Arabidopsis";
   const char *ara_gtf = "/home/ruolin/Dropbox/Strawberry/TAIR10_GFF3_genes-1.gff";
   //const char *human_gtf = "/home/ruolin/Downloads/gencode.v21.annotation.gff3";
   //const char *bam_file = "/home/ruolin/Dropbox/Strawberry/accepted_hits.bam";
   //char *bam_file = "/home/ruolin/git/Strawberry/RD100.control_r1.concordant_uniq.sort.bam";
   char *bam_file = "/home/ruolin/git/CompareTransAbun/assembly_comp/RD100/accepted_hits.bam";
   const char* bam_dir = stripFileName(bam_file);
   size_t len1 = strlen(bam_dir);
   char *output = (char*) malloc(len1+27);// assembled_transcripts.gtf 25 characters
   char *suffix = "/assembled_transcripts.gtf";
   strcpy(output, bam_dir);
   strcat(output, suffix);

   //const char *bam_file = "/home/ruolin/Dropbox/Strawberry/WetFT1.sm.bam";

   //FaInterface fa_api(path);
   //FaSeqGetter fsg;
   //fa_api.load2FaSeqGetter(fsg, "mitochondria");
   //cout<<"success\t"<<fsg.loadSeq()<<endl;
   //cout<<fsg.fetchSeq(80,4)<<endl;
   auto start = chrono::steady_clock::now();
   GffReader greader(ara_gtf);
   greader.readAll();
   greader.reverseExonOrderInMinusStrand();
   ReadTable read_table;
   RefSeqTable ref_seq_table(true);
   unique_ptr<HitFactory> hf(new BAMHitFactory(bam_file, read_table, ref_seq_table));
   hf->inspect_header();
   ClusterFactory read_clusters(move(hf));
   read_clusters.loadRefmRNAs(greader._g_seqs, ref_seq_table, path);
   FILE *pFile;
   pFile = fopen(output, "w");
   //QpSolver qps;

   read_clusters.inspectCluster();
   double mean, sd;
   const vector<int> & fd = read_clusters._hit_factory->_reads_table._frag_dist;
   unique_ptr<InsertSize> insert_size(new InsertSize(fd));
   //unique_ptr<InsertSize> insert_size(new InsertSize(300.0, 70.0));
   read_clusters._insert_size_dist = move(insert_size);
   cout<<"Total number of mapped reads is: "<<read_clusters.total_mapped_reads()<<endl;
   read_clusters.parseClusters(pFile);


   fclose(pFile);

   auto end = chrono::steady_clock::now();
   auto diff = end - start;
   cout << "Finished in " << chrono::duration <double, milli> (diff).count() << " ms" << endl;
}


