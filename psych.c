/*
 ============================================================================
 Name        : psych.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include "psych.h"

int main(void) {
	printf("75DB 65WB %fRH", psych(14.7,75,65,1,3,0));
	return EXIT_SUCCESS;
}
