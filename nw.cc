/*-------------------------------------------
 * 
 *        nw.c++ for program nw
 * 
 * 12/8/2011: ys nw-->smith waterman between guideRNA and edited transcripts
 -------------------------------------------*/

#include "nw.h"

using namespace std;


int nw(                                                          
       string       seq_1,          /*  Smith-Waterman   */
       string       seq_2,          /*  algorithm for      */
       string&      seq_1_al,       /*  local alignment   */
       string&      seq_2_al,       /*  of nt sequence.    */
       bool         prm,
       int& gRNAs,
       int& gRNAe,
       int& mRNAs,
       int& mRNAe
								 )
{
  int  d = 100000 ;                 /* gap penalty */

  int  L1 = seq_1.length();
  int  L2 = seq_2.length();

  // Dynamic programming matrix
  int ** F = new int * [ L2+1 ];
  for( int i = 0; i <= L2; i++ )  F[ i ] = new int [ L1 +1 ];

  // Traceback matrix
  char ** traceback = new char * [ L2+1 ];
  for( int i = 0; i <= L2; i++ )  traceback[ i ] = new char [ L1 +1 ];

  // Initialize traceback and F matrix (fill in first row and column)
  dpm_init( F, traceback, L1, L2, d );

  // Create alignment
  int align_score = nw_align( F, traceback, seq_1, seq_2, seq_1_al, seq_2_al, d, gRNAs, gRNAe, mRNAs, mRNAe );

#if DEBUG
  int  L_al = seq_1_al.length();
  cout << "Length after alignment: " << L_al << endl;
#endif

  if( prm )
    {
      cout << "\nDynamic programming matrix: " << "\n\n";
      print_matrix( F, seq_1, seq_2 );

      cout << "\nTraceback matrix: " << "\n\n";
      print_traceback( traceback, seq_1, seq_2 );

      cout << endl;
    }

  for( int i = 0; i <= L2; i++ )  delete F[ i ];
  delete[] F;
  for( int i = 0; i <= L2; i++ )  delete traceback[ i ];
  delete[] traceback;

  return  align_score ;
}

void  dpm_init( int ** F, char ** traceback, int L1, int L2, int d )
{
  F[ 0 ][ 0 ] =  0 ;
  traceback[ 0 ][ 0 ] = 'n' ;

  int i=0, j=0;

  for( j = 1; j <= L1; j++ )
    {
      F[ 0 ][ j ] = 0; // -j * d ;
      traceback[ 0 ][ j ] =  '-' ;
    }
  for( i = 1; i <= L2; i++ )
    {
      F[ i ][ 0 ] =  0; //-i * d ;
      traceback[ i ][ 0 ] =  '|' ;
    }
}

int nw_align(                  // Needleman-Wunsch algorithm
	     int **     F,
	     char **    traceback,
	     string     seq_1,
	     string     seq_2,
	     string&    seq_1_al,
	     string&    seq_2_al,
	     int        d,         // Gap penalty
	     int& gRNAs,
	     int& gRNAe,
	     int& mRNAs,
	     int& mRNAe
			       )
{
  int        k = 0, x = 0, y = 0;
  int        fU, fD, fL ;
  char       ptr, nuc ;
  int        i = 0, j = 0;

  const int  a =  -100000;   // Match
  const int  c = 1;
  const int  b = 2;   // Mismatch

  //A-U: 2 G-C: 2 G-U: 1 others -1
  const int  s[ 4 ][ 4 ] = { { a, a, a, b },    /* substitution matrix */
			     { a, a, b, a },
			     { a, b, a, c },
			     { b, a, c, a } } ;

  int  L1 = seq_1.length();
  int  L2 = seq_2.length();

  int max_score=-1;
  int maxi=-1;
  int maxj=-1;

  for( i = 1; i <= L2; i++ )
    {
      for( j = 1; j <= L1; j++ )
	{
	  nuc = seq_1[ j-1 ] ;

	  switch( nuc )
	    {
	    case 'A':  x = 0 ;  break ;
	    case 'a':  x = 0 ;  break ;
	    case 'C':  x = 1 ;  break ;
	    case 'c': x=1; break;
	    case 'G':  x = 2 ;  break ;		
	    case 'g': x=2; break;
	    case 'U':  x = 3 ;  break ;
	    case 'u': x=3; break;
	    case 'T':  x = 3 ;  break ;
	    case 't': x = 3; break;
	    case '*':  x = 4 ;
	    }

	  nuc = seq_2[ i-1 ] ;
	  
	  switch( nuc )
	    {
	    case 'A':  y = 0 ;  break ;
	    case 'a': y=0; break;
	    case 'C':  y = 1 ;  break ;
	    case 'c': y=1; break;
	    case 'G':  y = 2 ;  break ;
	    case 'g':y=2; break;
	    case 'T':  y = 3 ;  break ;
	    case 't':y=3; break;
	    case 'U': y=3; break;
	    case 'u': y=3; break;
	    case '*':  y = 4;
	    }

	  // * is aligned with gaps without penalty
	  if(x == 4)
	    {
	      F[i][j]=F[i][j-1];
	      ptr = '-';
	    }
	  else if(y == 4)
	    {
	      F[i][j] = F[i-1][j];
	      ptr = '|';
	    }
	  else
	    {
	      fU = F[ i-1 ][ j ] - d ;
	      fD = F[ i-1 ][ j-1 ] + s[ x ][ y ] ;
	      fL = F[ i ][ j-1 ] - d ;
	      F[ i ][ j ] = max( fU, fD, fL, &ptr ) ;
	    }
	  
	  if(F[i][j] < 0)
	    F[i][j]=0;

	  if(F[i][j]>max_score)
	    {
	      max_score=F[i][j];
	      maxi=i;
	      maxj=j;
	      //cerr<<"max_score = "<<max_score<<endl;
	    }
	  traceback[ i ][ j ] =  ptr ;
	}
    }
  //i-- ; j-- ;

  gRNAe=maxj-1;
  mRNAe=maxi-1;
  i=maxi; j=maxj;
  seq_1_al="";
  seq_2_al="";

#ifdef VERBOSE
  cout<<"maxscore "<<max_score<<" is at ("<<maxi<<","<<maxj<<")\n";
#endif

  //while( i > 0 || j > 0 )
  while(F[i][j]>0)
    {
      switch( traceback[ i ][ j ] )
	{
	case '|' :      seq_1_al += '-' ; 
	  seq_2_al += seq_2[ i-1 ] ; 
	  i-- ;
	  break ;

	case '\\':      seq_1_al += seq_1[ j-1 ] ; 
	  seq_2_al += seq_2[ i-1 ] ; 
	  i-- ;  j-- ;
	  break ;

	case '-' :      seq_1_al += seq_1[ j-1 ] ; 
	  seq_2_al += '-' ; 
	  j-- ;
	}
#ifdef VERBOSE
      cout<<"seq_1_al is "<<seq_1_al<<endl;
      cout<<"seq_2_al is "<<seq_2_al<<endl;
#endif

      k++ ;
    }

  if(F[i][j]==0)
    {
	gRNAs=j; //0-based counting
	mRNAs=i;
    }
        
  reverse( seq_1_al.begin(), seq_1_al.end() );
  reverse( seq_2_al.begin(), seq_2_al.end() );

  //cerr<<"!!max_score = "<<max_score<<endl;
  return  max_score ;
}


int  max( int f1, int f2, int f3, char * ptr )
{
  int  max = 0 ;

  if( f1 >= f2 && f1 >= f3 )  
    {
      max = f1 ;
      *ptr = '|' ;
    }
  else if( f2 > f3 )              
    {
      max = f2 ;
      *ptr = '\\' ;
    }
  else
    {
      max = f3 ;
      *ptr = '-' ;
    }
        
  return  max ;   
}


void  print_matrix( int ** F, string seq_1, string seq_2 )
{
  int  L1 = seq_1.length();
  int  L2 = seq_2.length();

  cout << "        ";
  for( int j = 0; j < L1; j++ )
    {
      cout << seq_1[ j ] << "   ";
    }
  cout << "\n  ";

  for( int i = 0; i <= L2; i++ )
    {
      if( i > 0 )
	{
	  cout << seq_2[ i-1 ] << " ";
	}
      for( int j = 0; j <= L1; j++ )
	{
	  cout.width( 3 );
	  cout << F[ i ][ j ] << " ";
	}
      cout << endl;
    }
}


void  print_traceback( char ** traceback, string seq_1, string seq_2 )
{
  int  L1 = seq_1.length();
  int  L2 = seq_2.length();

  cout << "    ";
  for( int j = 0; j < L1; j++ )
    {
      cout << seq_1[ j ] << " ";
    }
  cout << "\n  ";

  for( int i = 0; i <= L2; i++ )
    {
      if( i > 0 )
	{
	  cout << seq_2[ i-1 ] << " ";
	}
      for( int j = 0; j <= L1; j++ )
	{
	  cout << traceback[ i ][ j ] << " ";
	}
      cout << endl;
    }
}

string  print_al( string& seq_1_al, string& seq_2_al, int score, int gRNAs, int gRNAe, 
		int mRNAs, int mRNAe)
{
  //cout<<"mRNA: "<<mRNAs<<" "<<mRNAe<<" "<<"gRNA: "<<gRNAs<<" "<<gRNAe<<" S: "<<score<<" L: "<<seq_1_al.size()<<endl;
  //cout << seq_2_al << endl;
  
  string retstr=""; 
  for(int i=0; i<seq_1_al.size(); i++)
    {
      if(seq_1_al[i]=='A' and seq_2_al[i]=='U' or
	 seq_1_al[i]=='A' and seq_2_al[i]=='u' or
	 seq_1_al[i]=='U' and seq_2_al[i]=='A' or
	 seq_1_al[i]=='T' and seq_2_al[i]=='A' or
	 seq_1_al[i]=='G' and seq_2_al[i]=='C' or
	 seq_1_al[i]=='C' and seq_2_al[i]=='G')
	//cout<<'|';
	retstr+="|";
      else if(seq_1_al[i]=='G' and seq_2_al[i]=='U' or
	      seq_1_al[i]=='G' and seq_2_al[i]=='u' or
	      seq_1_al[i]=='T' and seq_2_al[i]=='G' or
	      seq_1_al[i]=='U' and seq_2_al[i]=='G')
	//cout<<':';
	retstr+=":";
      else if(seq_1_al[i]=='*' or seq_2_al[i]=='*')
	//cout<<" ";
	retstr+=" ";
      else
	//cout<<'.';
	retstr+=".";
    }
  //cout<<endl;
  //cout << seq_1_al << endl;
  return retstr;
}
