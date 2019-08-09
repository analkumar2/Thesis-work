/******************************************************************************
 * 
 * file to read in neuron morphology file
 * prints out the radius, length, and distance from soma for each dendrite
 *
 *******************************************************************************/

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <iomanip>
#include <string.h>

using namespace std;
const double PI = 3.14159;

int main(int argc, char * const argv[]){

  FILE * cfptr;
  cfptr = fopen("n123.txt", "r");
  
  FILE * cfptr2;
  cfptr2 = fopen("n123_trim.txt", "w");


  /* Axial resistivity (?) for calculating axial resistance of each compartment */
  double Ra = 50.0e4;  // ohm*um; baseline somatic resistivity from poirazi et al 
  double resist;

  /* variables for reading in and manipulating string */
  char line[50];
  char character_start = '(';
  const char * character_end = ")}";
  char * tokenPtr;
  char * start;
  char * end;

  /* variables for output */
  int segment=0;
  double x, y, z;
  double xnew, ynew, znew, radiusnew, radiusold;
  double dx, distance, radius;
  double count;
  int count_seg = 0;

  fscanf(cfptr, "%s", line);

  while( !feof( cfptr ) ) {
    
    //cout << "line: " << line << endl; 

    /* check if line is a number, indictating the segment number: */
    if ( isdigit( line[0] ) ) {

      /* write previous segment and stats to file */
      if( count_seg > 0) fprintf( cfptr2, "%d\n%.2f\n%.2f\n%.2f\n\n", segment, dx, radius, resist );

      segment = atoi( line );
      //cout << "segment: " << segment << endl;
      count = 0.;
      dx = 0.0; 
      radius = 0.0;
      resist = 0.0;
      count_seg++;
    }
    fscanf(cfptr, "%s", line);

    while( ! isdigit( line[0]) && (!feof( cfptr )) ){
      start = strchr( line, character_start );
      
      if ( start != NULL ){
	start++;
	int size = strcspn( start, character_end );
	
	if (size > 0 ){
	  
	  //cout << "start isn't null, start: " << start << endl;
	  tokenPtr = strtok( start, "," );
	  xnew = atof( tokenPtr );
	  //cout << "x: " << xnew << "\t";
	  tokenPtr = strtok( NULL, "," );
	  ynew = atof( tokenPtr );
	  //cout << "y: " << ynew << "\t";
	  tokenPtr = strtok( NULL, ",");
	  znew = atof( tokenPtr );
	  //cout << "z: " << znew << "\t";
	  tokenPtr = strtok(NULL, ")");
	  radiusnew = atof( tokenPtr );
	  //cout << "radius: " << radiusnew << endl;
	  
	  if ( count < 1. ){
	    x = xnew;
	    y = ynew;
	    z = znew;
	    radius = radiusnew;
	    radiusold = radiusnew;
	  }
	  if ( count > 0.) {
	    double deltax = sqrt( pow((xnew-x),2) + pow((ynew-y),2) + pow((znew-z),2) );
	    dx = dx + deltax;
	    x = xnew;
	    y = ynew;
	    z = znew;
	    radius = radius*(count/(count+1.)) + radiusnew*(1./(count+1.));
	    resist = resist + Ra*deltax/(PI*radiusnew*radiusold);  // resistance = resistivity*length/area;  ohm = ohm*um*um/um2
	    radiusold = radiusnew;
	    //cout << "radius: " << radius << " radius new: " << radiusnew<< endl;
	  }
	  count++;
	 
	}

      }
      fscanf( cfptr, "%s", line );
    } // ends check for whether first character in the line is a digit

  } // matches while(!feof( cfptr))  
  
  /* write the last segment info to file */
  fprintf( cfptr2, "%d\n%.2f\n%.2f\n%.2f\n\n", segment, dx, radius, resist);
  
  fclose(cfptr);
  fclose(cfptr2);

  return 0;

}
