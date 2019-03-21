/*
 * nw.h for program nw.
 *
 */

#include <iostream>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <vector>

#define DEBUG 0

using namespace std;

/*
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
*/
/*
extern void nw_LCS(
		string, string, int, int, int,
		string, char*, string, int, bool,
               vector< struct record>&
		);
*/
extern int  nw( 
	       string, string, 
	       string&, string&,
	       bool, int&, int&, int&, int&
		);

extern int  nw_align( 
		     int **, char **,
		     string, string, 
		     string&, string&,
		     int, int&, int&, int&, int& 
		      );

extern void  dpm_init        ( int **, char **, int, int, int );
//extern void dpm_init_LCS(int ***,  char **, int, int);

//extern void  print_al        ( string&, string&, int, int, int, int, int );
extern string  print_al        ( string&, string&, int, int, int, int, int );
extern void  print_matrix    ( int ** const, string, string );
extern void  print_traceback ( char ** const, string, string );
extern int   max             ( int, int, int, char * );
//extern void  print_3Dmatrix  (int *** const, string, string );


