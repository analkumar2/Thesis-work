/*
 *  get_morphology.cpp
 *  
 *
 *  Created by Naomi Lewin on 11/14/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include "get_morphology.h"
#include "neuron_structures.h"

using namespace std;

void calc_distance(Compartment * neuron){

  int seg;
  int prev_neighbor = 0;
  double dist;

	FILE * dist_list = fopen("dist_list", "w"); 

  for(seg = 0; seg < 5; seg++){
    neuron[seg].distance = 0.0; // these are the somatic compartments therefore the distance from soma = 0;
	fprintf(dist_list, "%d\t%f\t%f\t%f\n", seg, neuron[seg].distance, neuron[seg].radius, neuron[seg].dx);
  }

  for(seg = 5; seg < 183; seg++){
    if(neuron[seg].num_neighbors == 1) prev_neighbor = neuron[seg].neighbor_list[0];
    if(neuron[seg].num_neighbors == 2) prev_neighbor = neuron[seg].neighbor_list[1];
    if(neuron[seg].num_neighbors == 3) prev_neighbor = neuron[seg].neighbor_list[2];
    if(prev_neighbor > seg) cout << "ERROR: previous neighbor has higher value than segment number at segment " << seg << endl;
    neuron[seg].distance = .5*(neuron[seg].dx + neuron[prev_neighbor].dx) + neuron[prev_neighbor].distance;
	fprintf(dist_list, "%d\t%f\t%f\t%f\n", seg, neuron[seg].distance, neuron[seg].radius, neuron[seg].dx);
    //cout << "seg: " << seg << " previous neighbor: " << prev_neighbor << endl << " distance: " << neuron[seg].distance << endl;;
  }

}

void reset_grid(Compartment * neuron, Compartment * new_neuron, int * create_new){

    // function to reset the maximum dx of the compartments
    // uses previously calculated distance, dx, neighbor list
    // creates new compartments
    // resets dx
    // resets neighbors
    // resets total number of compartments
    // count_new is the number of new compartments created in order to reach desired maximum dx
	// set max size of compartments:


    int counter = 0;
	int final_index;
    for(int seg=0; seg<183; seg++){
        new_neuron[seg].dx = neuron[seg].dx/(create_new[seg]+1);  // new length based on max grid
        new_neuron[seg].radius = neuron[seg].radius;
        new_neuron[seg].distance = neuron[seg].distance;
        new_neuron[seg].num_neighbors = neuron[seg].num_neighbors;
        //if(create_new[seg] == 0){
		// copy neighbors:
		if(neuron[seg].num_neighbors == 1) new_neuron[seg].neighbor_list[0] = neuron[seg].neighbor_list[0];
		if(neuron[seg].num_neighbors == 2) {
                new_neuron[seg].neighbor_list[0] = neuron[seg].neighbor_list[0];
                new_neuron[seg].neighbor_list[1] = neuron[seg].neighbor_list[1];                
		}
		if(neuron[seg].num_neighbors == 3) {
                new_neuron[seg].neighbor_list[0] = neuron[seg].neighbor_list[0];
                new_neuron[seg].neighbor_list[1] = neuron[seg].neighbor_list[1]; 
                new_neuron[seg].neighbor_list[2] = neuron[seg].neighbor_list[2];
		}
        //}  // copied variables for unchanged compartments into new structure
	}
	for(int seg = 0; seg<183; seg++){
		// fill in the variables for the new compartments corresponding to each of the original compartments
        if(create_new[seg]!=0){
			int number_ofnews = create_new[seg];
            int new_seg[number_ofnews];
            for(int i=0; i<create_new[seg]; i++){
                int index = 183 + counter;
                new_seg[i] = index;  // list of compartment numbers of new compartments corresponding to the old segment number we're currently dealing with
                counter++;
                new_neuron[index].dx = new_neuron[seg].dx;
                new_neuron[index].radius = new_neuron[seg].radius;
                new_neuron[index].distance = new_neuron[seg].distance;
            }
                
            // for each condition (1, 2, or 3 neighbors) of the original compartment, 
            // set neighbors for original compartment:
            // set neighbors for inital new compartment:
            // set neighbors for middle new compartments:
            // set neighbors for final new compartment:
            if(neuron[seg].num_neighbors == 1) {
                new_neuron[seg].num_neighbors = 1;
                new_neuron[seg].neighbor_list[0] = new_seg[0];  // new neighbor is first new compartment
                int neighbor_index = neuron[seg].neighbor_list[0];
                for(int i=0; i<create_new[seg]; i++){
                    int index = new_seg[i];
                    if(i == 0){
                        new_neuron[index].neighbor_list[1]=seg; //[0] = seg;
                        if(create_new[seg]>1) new_neuron[index].neighbor_list[0]=new_seg[i+1]; //[1] = new_seg[i+1];
                        else new_neuron[index].neighbor_list[0] = neuron[seg].neighbor_list[0]; //[1] = neuron[seg].neighbor_list[0];
                        new_neuron[index].num_neighbors = 2;
                    }
                    else if(i == create_new[seg]-1){
                        new_neuron[index].neighbor_list[0] = neuron[seg].neighbor_list[0]; //new_seg[i-1];
                        new_neuron[index].neighbor_list[1] = new_seg[i-1]; //neuron[seg].neighbor_list[0];
                        new_neuron[index].num_neighbors = 2;
                    }
                    else {
                        new_neuron[index].neighbor_list[0] = new_seg[i+1]; //new_seg[i-1];
                        new_neuron[index].neighbor_list[1] = new_seg[i-1]; //new_seg[i+1];
                        new_neuron[index].num_neighbors = 2;
                    }
                }
                // replace the index in the original neighbor's list:
                for(int i=0; i<neuron[neighbor_index].num_neighbors; i++){
                    if(neuron[neighbor_index].neighbor_list[i] == seg) new_neuron[neighbor_index].neighbor_list[i] = new_seg[number_ofnews-1];
                }
            
            }

            else if(neuron[seg].num_neighbors == 2) {
                new_neuron[seg].num_neighbors = 2;
                new_neuron[seg].neighbor_list[0] = new_seg[0]; //neuron[seg].neighbor_list[0];
                //new_neuron[seg].neighbor_list[1] = neuron[seg].neighbor_list[1]; //new_seg[0];  
                int neighbor_index = neuron[seg].neighbor_list[0]; //neuron[seg].neighbor_list[1];
                for(int i=0; i<create_new[seg]; i++){
                    int index = new_seg[i];
                    if(i == 0){
                        new_neuron[index].neighbor_list[1] = seg; //new_neuron[index].neighbor_list[0] = seg;
                        if(create_new[seg]>1) new_neuron[index].neighbor_list[0] = new_seg[i+1]; //neighbor_list[1] = new_seg[i+1];
                        else new_neuron[index].neighbor_list[0] = neuron[seg].neighbor_list[0]; //neighbor_list[1] = neuron[seg].neighbor_list[1];
                        new_neuron[index].num_neighbors = 2;
                    }
                    else if(i == create_new[seg]-1){
                        new_neuron[index].neighbor_list[1] = new_seg[i-1]; //[0] = new_seg[i-1];
                        new_neuron[index].neighbor_list[0] = neuron[seg].neighbor_list[0]; //[1] = neuron[seg].neighbor_list[1];
                        new_neuron[index].num_neighbors = 2;
                    }
                    else {
                        new_neuron[index].neighbor_list[0] = new_seg[i+1]; //new_seg[i-1];
                        new_neuron[index].neighbor_list[1] = new_seg[i-1]; //new_seg[i+1];
                        new_neuron[index].num_neighbors = 2;
                    }
                }
                // replace the index in the original neighbor's list:
                for(int i=0; i<neuron[neighbor_index].num_neighbors; i++){
                    if(neuron[neighbor_index].neighbor_list[i] == seg) new_neuron[neighbor_index].neighbor_list[i] = new_seg[number_ofnews-1];
                }
            }
            
            else if(neuron[seg].num_neighbors == 3) {
                new_neuron[seg].num_neighbors = 2;
                new_neuron[seg].neighbor_list[0] = new_seg[0]; //neuron[seg].neighbor_list[0];
                new_neuron[seg].neighbor_list[1] = new_neuron[seg].neighbor_list[2]; //new_seg[0];  
                int neighbor_index1 = neuron[seg].neighbor_list[0]; //[1];
                int neighbor_index2 = neuron[seg].neighbor_list[1]; //[2];
                for(int i=0; i<create_new[seg]; i++){
                    int index = new_seg[i];
                    if(i == 0){
                        if(create_new[seg]>1) {
                            new_neuron[index].num_neighbors = 2;
							new_neuron[index].neighbor_list[1]=seg; //[0] = seg;
                            new_neuron[index].neighbor_list[0]=new_seg[i+1]; //[1] = new_seg[i+1];
                        }
                        else {
							new_neuron[index].neighbor_list[0] = neuron[seg].neighbor_list[0];
                            new_neuron[index].neighbor_list[1] = neuron[seg].neighbor_list[1];
                            new_neuron[index].neighbor_list[2] = seg; //neuron[seg].neighbor_list[2];
                            new_neuron[index].num_neighbors = 3;
							final_index = index;
                        }
                    }
                    else if(i == create_new[seg]-1){
                        new_neuron[index].neighbor_list[0] = neuron[seg].neighbor_list[0]; //new_seg[i-1];
                        new_neuron[index].neighbor_list[1] = neuron[seg].neighbor_list[1];
                        new_neuron[index].neighbor_list[2] = new_seg[i-1]; //neuron[seg].neighbor_list[2];
                        new_neuron[index].num_neighbors = 3;
						final_index = index;
                    }
                    else {
                        new_neuron[index].neighbor_list[0] = new_seg[i+1]; //new_seg[i-1];
                        new_neuron[index].neighbor_list[1] = new_seg[i-1]; //new_seg[i+1];
                        new_neuron[index].num_neighbors = 2;
                    }
                }
                // replace the index in the original neighbor's list:
                for(int i=0; i<neuron[neighbor_index1].num_neighbors; i++){
					cout << "segment: " << seg << endl;
                    if(neuron[neighbor_index1].neighbor_list[i] == seg) {
						new_neuron[neighbor_index1].neighbor_list[i] = new_seg[number_ofnews-1];
						cout << "for segment " << neighbor_index1 << " new neighbor " << final_index << " at neighbor #: " << i  << endl;
					}
				}
                for(int i=0; i<neuron[neighbor_index2].num_neighbors; i++){
                    if(neuron[neighbor_index2].neighbor_list[i] == seg){
						new_neuron[neighbor_index2].neighbor_list[i] = new_seg[number_ofnews-1];
						cout << "for segment " << neighbor_index2 << " new neighbor " << final_index << " at neighbor #: " << i  << endl;
					}
				}
			}
		
		}  // filled in variables for changed and new compartments into new structure
    }
    
}


void read_dimensions(Compartment * neuron){

  FILE * fptr;
  fptr = fopen("n123_trim.txt", "r");
  char string[30];
  int segment;
  double length, rad, resist;

  fscanf(fptr, "%s", string);
  while( !feof( fptr ) ){
    segment = atoi( string );
    fscanf(fptr, "%s", string);
    length = atof( string );
    neuron[segment].dx = length;
    fscanf(fptr, "%s", string);
    rad = atof( string );
	  neuron[segment].radius = rad/2.; //rad/2.;
    fscanf(fptr, "%s", string);
    resist = atof( string );
    neuron[segment].ra = resist; // this isn't used; axial resistance set later based on distance from soma
    fscanf(fptr, "%s", string);
  }

  fclose( fptr );
 
}

void read_neighborlist(Compartment * neuron){
  
  FILE * fptr1;
  fptr1 = fopen("n123_neighborslist.txt", "r");

  char string[30];
  char * tokenPtr;
  int i, segment, neighbor[3];

  fscanf( fptr1, "%s", string);
  while( !feof( fptr1 ) ){

    for(i = 0; i<3; i++) neighbor[i] = -1;
    tokenPtr = strtok( string, ",");
    segment = atoi( tokenPtr );
    i = 0;
    while (tokenPtr != NULL ) {
      tokenPtr = strtok( NULL, ",");
      if (tokenPtr != NULL ) {
	neighbor[i] = atoi( tokenPtr );
	i++;
      }
    }
  
    neuron[segment].num_neighbors = i;
    for( int i = 0; i<3; i++) neuron[segment].neighbor_list[i] = neighbor[i];
    fscanf( fptr1, "%s", string);
  }


  /* do last line: */
  for(i = 0; i<3; i++) neighbor[i] = -1;
  tokenPtr = strtok( string, ",");
  segment = atoi( tokenPtr );
  i = 0;
  while (tokenPtr != NULL ) {
    tokenPtr = strtok( NULL, ",");
    if (tokenPtr != NULL ) {
      neighbor[i] = atoi( tokenPtr );
      i++;
    }
  }

  neuron[segment].num_neighbors = i;
  for( int i = 0; i<3; i++) neuron[segment].neighbor_list[i] = neighbor[i];

  fclose( fptr1 );

}
