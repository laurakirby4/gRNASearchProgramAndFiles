/*---------------------------------------------------------------
 *
 *   main.c++ for nw program.
 *
 *   Implements the Smith-Waterman algorithm
 *   for local nucleotide sequence alignment.
 *
 *   adpated from Rolf Muertter,  rolf@dslextreme.com
 *   12/10/2011
 *
 *   sort all alignments in the order of alignment starting positions
 *   in mRNA; output all gRNAs in the same order

 *   2/6: corrected an error in file openning. Will this affect previous output?
 *   9/2/2012: note that the alignment is conducted between the "prefix" and the 
     edited mRNA. The prefix is substring ending at TTTT. Sometimes, TTTT occurs towards 5'
     end and can lead to shorter alignments than the actual ones. 
     Alignment using only the prefix is used throughout other programs (such as LCS).

 *   3/14/2014: prefixU = read.substr(0, pos1+4). This will include the prefix of the read
     up to TTTT. Currently, prefixU is the prefix up to "TTT". 
 ---------------------------------------------------------------*/

#include <fstream>
#include <vector>

#include "nw.h"
#include "fasta.h"

using namespace std;

struct record{
  int mRNAs;
  int mRNAe;
  int gRNAs;
  int gRNAe;
  string gRNA;
  int score;
  int L;
  string gname;
  string mname;
  int prefixL;
  string seq_1_al;
  string seq_2_al;
  string mid;
};

bool comparison(struct record one, struct record two)
{
  if(one.mRNAs != two.mRNAs)
    return (one.mRNAs < two.mRNAs);
  else if(one.mRNAs == two.mRNAs)
    return (one.mRNAe < two.mRNAe);
}

int  main( int argc, char ** argv )
{
  char *  program = *argv ;
  bool    prm = false;

  if( argc < 6 )
    {
      cerr << "\n   Usage:   " << program << " gRNA_fasta_file edited_transcript_file scoreT ali_L_cutoff gRNA-file gRNA-fasta-file [-p]\n";
      cerr << "\n   -p:       Print matrices\n";
      cerr << "\n   Output:   alignment\n\n";
      exit( 1 ) ;
    }                       
 
  if( argc == 7 )
    {
      string  pr_arg  =  argv[ 6 ] ;
      if( pr_arg == "-p" )  prm = true;   // Print matrices                                                      
    }
       
  FASTAFILE *ffp;
  char *seq;
  char *name;
  int L;
  vector<string> editseq;
  vector<string> editname;

  int Tcutoff=atoi(argv[3]);
  float Lcutoff=atof(argv[4]);
  //int Lcutoff=atoi(argv[4]);

  ffp = OpenFASTA(argv[2]);
  while (ReadFASTA(ffp, &seq, &name, &L))
    {

#ifdef VERBOSE
      printf("%s %d %s\n", name, L, seq);
#endif
      string edit(seq);
      editseq.push_back(edit);
      editname.push_back(name);

#ifdef VERBOSE
      cout<<"the size of the set is "<<editseq.size()<<endl;
#endif
      free(seq);
      free(name);
    }
  CloseFASTA(ffp);

  //read the putative gRNA fasta file; align each sequence with edit sequence
  vector<struct record> aliarray;
  FASTAFILE *ffp1;
  ffp1 = OpenFASTA(argv[1]);
  while (ReadFASTA(ffp1, &seq, &name, &L))
    {

#ifdef VERBOSE
      printf("%s %d %s\n", name, L, seq);
#endif
      string read(seq);
 
      // Aligned sequences                                                                             
      string  seq_1_al="";
      string  seq_2_al="";

      // size of gRNA without polyU tail                                                               
      size_t pos1;
      string prefix1=read;
      pos1=read.rfind("TTTT");
      if (pos1==string::npos)
	{
	    cerr<<"Error in finding TTTT\n";                                                             
	    exit(1);                                                                                     
	}
      else
	{
	  int startT1=pos1;
	  while(startT1>=0 and read[startT1]=='T')
	    {
	      startT1--;
	    }
	  if(startT1>=0)
	    prefix1=read.substr(0, startT1+1);
	}

      //          3'-----------5' gRNA
      // 5'------------------------------------3' edited sequence
      // process the putative gRNA: reverse the prefix up to the polyT tail
      //reverse(read.begin(), read.end());
      string prefixU=read.substr(0, pos1+3);
      reverse(prefixU.begin(), prefixU.end());

      // Get alignment
      for(int i=0; i<editseq.size(); i++)
	{
	  string seq_2=editseq[i];
	  string seq_2_name = editname[i];
	  int gRNAs, gRNAe;
	  int mRNAs, mRNAe;
	 
#ifdef VERBOSE
	  cout<<"!conduct alignment between "<<prefixU<<" and "<<seq_2<<endl;
#endif
	  int score = nw(prefixU, seq_2, seq_1_al, seq_2_al, prm, gRNAs, gRNAe, mRNAs, mRNAe ) ;
	  
	  // Print alignment when the alignment score/length are above given cutoffs
	  if(score >= Tcutoff and (float)seq_1_al.size()/(float)prefix1.size()>=Lcutoff)
	    //if(score >= Tcutoff and seq_1_al.size()>=Lcutoff)  
	  {
	    struct record onealign;
	    onealign.mRNAs = mRNAs;
	    onealign.mRNAe = mRNAe;
	    onealign.gRNAs = gRNAs;
	    onealign.gRNAe = gRNAe;
	    onealign.gRNA = read;
	    onealign.score = score;
	    onealign.L = seq_1_al.size();
	    string namestr(name);
	    onealign.gname=namestr;
	    onealign.mname=seq_2_name;
	    onealign.prefixL=prefix1.size();
	    onealign.seq_1_al=seq_1_al;
	    onealign.seq_2_al=seq_2_al;

	    //cout<<"> "<<seq_2_name<<" "<<name<<" "<<prefix1.size()<<endl;
	    onealign.mid=print_al( seq_1_al, seq_2_al, score, gRNAs, gRNAe, mRNAs, mRNAe ) ;
	    
	    aliarray.push_back(onealign);
	    
	    }
	  //#endif
	}

      free(seq);
      free(name);
    }
  CloseFASTA(ffp1);

  FILE* gRNAffp;
  gRNAffp=fopen(argv[5], "w");

  //print using sorted information
  sort(aliarray.begin(), aliarray.end(), comparison);
  for(int i = 0; i<aliarray.size(); i++)
    {
      struct record onealign=aliarray[i];
      fprintf(gRNAffp, ">%s\n", onealign.gname.c_str());
      fprintf(gRNAffp,"%s\n", onealign.gRNA.c_str());
      cout<<"> "<<onealign.mname<<" "<<onealign.gname<<" "<<onealign.prefixL<<endl;
      cout<<"mRNA: "<<onealign.mRNAs<<" "<<onealign.mRNAe<<" "<<" gRNA: "<<onealign.gRNAs<<" "
	  <<onealign.gRNAe<<" S: "<<onealign.score<<" L: "<<onealign.L<<endl;
      cout<<onealign.seq_2_al<<endl;
      cout<<onealign.mid<<endl;
      cout<<onealign.seq_1_al<<endl;
    }
  fclose(gRNAffp);
     
  return  0 ;
}
