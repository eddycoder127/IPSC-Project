# Fast Fourier Transform
## Overall Description
Contains Serial and Parallel Version (openMPI) of the Fast Fourier Transform algorithm. Used to converts a signal from its original domain (often time or space) to a representation in the frequency domain and vice versa. A performance analysis was performed to show display execution time, speed-up, and efficiency.

## Serial Code
Sequentially performs the FFT algorithm and records time. Execution time was recorded using the TimeSource.c and timer.h files.
    
    *LOOK AT FFT-Serial Version.c
## Parallel Code
Parallel implementation of the FFT written in C using openMPI. Execution time was recorded using MPI's time function MPI_Wtime(). 
    
    *LOOK AT FFT-ParallelVersion.c
### Master Process
The master process creates a table of real and imaginary numbers. The master process then distributes a portion of the table to each individual process. The master process then gathers all the calculated information from the slave processes. Then calculates the last part of the FFT which requires all the imaginary and real numbers to be summed. The master process is in charge of recording the execution time.
### Slave Processes
The slaves processes recieves a portion of the table that hold imaginary and real values. The slave processes then perform the FFT algorithm on thier portion. This portion is later gathered by the master process.

## Performance Analysis Breakdown
### Execution Time
Execution time goes down with the more processors that are being utilized. The exeception is when number of processors is > 16. This is due to each node in the cluster that I used has 16 possible cores. When there is more than 16 processors the overhead of communication between nodes slows then execution down.
![](https://i.imgur.com/dFeWrAB.png)
### Speed-Up
Speedup is greatest when using the max amount of processes in a single node. In this case speedup is greatest using 16 proceses. 
![](https://i.imgur.com/fDoZf3y.png)
### Efficiency
Efficiency decreases in these examples due to the relativly small size of the problem. In this case problem size never goes beyond 16384. Efficiency decreases due to each process spending more time waiting than computing. If the problem size was bigger then efficiency would increase due to each process spending more time computing than waiting. 
![](https://i.imgur.com/YJjKvIp.png)
