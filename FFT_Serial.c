/*
Braden Harrelson
Kevin Ellis
Midwestern State University
High Performance Computing
FFT SERIAL IMPLEMETATION
*/
/* This program sequentialy computes the FFT
	Runs the program X amount of times and computes the time
	for each iteration. Then averages the time taken.
	
	COMPILE: mpicc -o seqverexe sequentialVersion.c -lm *NOTE: -lm * allows use of library that uses cos and sin
	RUN: qsub mpiJob1N1P
*/
#include <stdio.h>
#include <complex.h> //to use complex numbers
#include <math.h>	//for cos() and sin()
#include "Timer.h" //to use timer

#define PI 3.14159265
#define bigN 16384 //Problem Size
#define howmanytimesavg 3 //How many times do I wanna run for the AVG?
int main()
{	
	double avgtime = 0;
	int h;
	FILE *outfile;
	outfile = fopen("sequentialVersionOutput.txt", "w"); //oepn from current directory

	for(h = 0;h < howmanytimesavg; h++ )  //loop how many times you want to avg over
	{
		double start,finish; //For time
		GET_TIME(start); //start the timer
											
		double table[bigN][3] = 
								{
								 0,3.6,2.6, //n, Real,Imaginary CREATES TABLE
								 1,2.9,6.3,
								 2,5.6,4.0,
								 3,4.8,9.1,
								 4,3.3,0.4,
								 5,5.9,4.8,
								 6,5.0,2.6,
								 7,4.3,4.1,
								 };
							  
		double complex evenpart[bigN/2]; //array to save the data for EVENHALF
		double complex oddpart[bigN/2]; //array to save the data for ODDHALF
		double storeKsumreal[bigN]; //store the K real variable so we can abuse symmerty
		double storeKsumimag[bigN]; //store the K imaginary variable so we can abuse symmerty
		int k, i ,j; 
		if(bigN > 8)  //Everything after row 8 is all 0's
		{
			for(i = 8; i < bigN; i++)
			{
				table[i][0] = i;
				for(j = 1; j < 3;j++)
				{
					table[i][j] = 0.0; //set to 0.0
				}
			}
		}

		for (k = 0; k < bigN / 2; k++ ) //K loop
		{	
			/* Variables used for the computation */
			double sumrealeven = 0.0; //sum of real numbers for even
			double sumimageven = 0.0; //sum of imaginary numbers for even
			double sumrealodd = 0.0; //sum of real numbers for odd
			double sumimagodd = 0.0; //sum of imaginary numbers for odd
			
			for (i = 0; i <= (bigN/2 - 1); i++) //loop for series 0->N/2 -1
			{
				/* -------- EVEN PART -------- */
				double realeven = table[2*i][1]; //Access table for real number at spot 2i
				double complex imaginaryeven = table[2*i][2]; //Access table for imaginary number at spot 2i
				double complex componeeven = (realeven + imaginaryeven * I); //Create the first component from table

				double factoreven = ((2*PI)*((2*i)*k))/bigN; //Calculates the even factor for Cos() and Sin()
				
				double complex comptwoeven = (cos(factoreven) - (sin(factoreven)*I)); //Create the second component

				evenpart[i] = (componeeven * comptwoeven); //store in the evenpart array
				
				/* -------- ODD PART -------- */
				double realodd = table[2*i + 1][1]; //Access table for real number at spot 2i+1
				double complex imaginaryodd = table[2*i + 1][2]; //Access table for imaginary number at spot 2i+1
				double complex componeodd = (realodd + imaginaryodd * I); //Create the first component from table
				
				double factorodd = ((2*PI)*((2*i+1)*k))/bigN;//Calculates the odd factor for Cos() and Sin()
															
				double complex comptwoodd = (cos(factorodd) - (sin(factorodd)*I));//Create the second component

				oddpart[i] = (componeodd * comptwoodd); //store in the oddpart array
				
			}
			
			for(i = 0; i < bigN/2; i++) //loop to sum the EVEN and ODD parts
			{
				sumrealeven += creal(evenpart[i]); //sums the realpart of the even half
				sumimageven += cimag(evenpart[i]); //sums the imaginarypart of the even half
				
				sumrealodd += creal(oddpart[i]); //sums the realpart of the odd half
				sumimagodd += cimag(oddpart[i]); //sums the imaginary part of the odd half
			}
			
			storeKsumreal[k] = sumrealeven + sumrealodd; //add the calculated reals from even and odd
			storeKsumimag[k] = sumimageven + sumimagodd; //add the calculated imaginary from even and odd
			
			storeKsumreal[k + bigN/2] = sumrealeven - sumrealodd; //ABUSE symmetry Xkreal + N/2 = Evenk - OddK
			storeKsumimag[k + bigN/2] = sumimageven - sumimagodd; //ABUSE symmetry Xkimag + N/2 = Evenk - OddK
			if(k <= 10) //Do the first 10 K's
			{
				if(k == 0)
				{
					fprintf(outfile," \n\nTOTAL PROCESSED SAMPLES : %d\n",bigN);
				}
				fprintf(outfile,"================================\n");
				fprintf(outfile,"XR[%d]: %.4f XI[%d]: %.4f \n",k,storeKsumreal[k],k,storeKsumimag[k]);
				fprintf(outfile,"================================\n");
			}
		}
		GET_TIME(finish); //stop timer
		double timeElapsed = finish-start; //Time for that iteration
		avgtime = avgtime + timeElapsed; //AVG the time 
		fprintf(outfile,"Time Elaspsed on Iteration %d: %f Seconds\n", (h+1),timeElapsed);
	}
	avgtime = avgtime / howmanytimesavg;
	fprintf(outfile,"\nAverage Time Elaspsed: %f Seconds", avgtime);
	fclose(outfile); //close file
	return 0;
}
