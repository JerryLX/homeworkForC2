/***************************************************************************
 *
 *   File        : main.c
 *   Student Id  : <INSERT STUDENT ID HERE>
 *   Name        : <INSERT STUDENT NAME HERE>
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include "tasks.h"

int main(int argc, char *argv[]) {
	
	/* TODO: Parse Command Line Arguments
	DONOT explicitly set arguments to filenames */
	char* q2_file = NULL;
	char* q4_file = NULL;
	char* q5_file = NULL;
	double xo;
	char* q6_file = NULL;
    
    struct timeval timemark1,timemark2,timemark3,timemark4,timemark5;
    float elapsed1,elapsed2,elapsed3,elapsed4;
    
    if(argc < 6){
        printf("USAGE: %s <q2_file><q4_file><q5_file><xo><q6_file>\n", argv[0]);
        return (EXIT_FAILURE);
    }
    q2_file = argv[1];
    q4_file = argv[2];
    q5_file = argv[3];
    xo = atoi(argv[4]);
    q6_file = argv[5];
    
	/* Question 2 */
    gettimeofday(&timemark1,NULL);
	shockwave(q2_file);
    gettimeofday(&timemark2,NULL);
    elapsed1 = (timemark2.tv_usec - timemark1.tv_usec) / 1000.0 
        + (timemark2.tv_sec - timemark1.tv_sec) * 1000.0;
    printf("Question 2: %.2f milliseconds\n", elapsed1);
	
	/* Question 4 */
	linalgbsys(q4_file);
    gettimeofday(&timemark3,NULL);
    elapsed2 = (timemark3.tv_usec - timemark2.tv_usec) / 1000.0 
        + (timemark3.tv_sec - timemark2.tv_sec) * 1000.0;
    printf("Question 4: %.2f milliseconds\n", elapsed2);
	
	/* Question 5 */
	interp(q5_file,xo);
    gettimeofday(&timemark4,NULL);
    elapsed3 = (timemark4.tv_usec - timemark3.tv_usec) / 1000.0 
        + (timemark4.tv_sec - timemark3.tv_sec) * 1000.0;
    printf("Question 5: %.2f milliseconds\n", elapsed3);
	
	/* Question 6 */
	heateqn(q6_file);
    gettimeofday(&timemark5,NULL);
    elapsed4 = (timemark5.tv_usec - timemark4.tv_usec) / 1000.0 
        + (timemark5.tv_sec - timemark4.tv_sec) * 1000.0;
    printf("Question 6: %.2f milliseconds\n", elapsed4);
    
	return (EXIT_SUCCESS);
}
