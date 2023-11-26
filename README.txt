Steps to run:
1. Run "make build" in terminal
2. An executable called "pagerank" will be created"
3. Run this executable in the following format:
    "./pagerank {input file} {output file} {type of synchronization} {number of threads}"

Types of synchronization:
The type of synchronization can either be:
0: Serial - No multi-threading
1: Mutex - Synchronization using mutexes
2: Spinlock - Synchronization using spinlocks
3: Atomic instructions - Synchronization using atomic instructions
4: Atomic instruction + load balancing - Same as atomic instructions with load 
   balancing by edges

Use the number corresponding to the desired type of synchronization in the place
of {type of synchronization} when running the executable

Input files:
The input files are in the "inputs" subdirectory. Output files should be placed
in the "outputs" subdirectory, however no restrictions are in place to enforce
that. 

Examples:
"./pagerank inputs/rmat15.dimacs outputs/rmat15/mutex/4_thread 1 4
 - This command will take the rmat15.dimacs file as input, and compute the page rank
   using mutex synchronization (1) with 4 threads and output the result to
   outputs/rmat15/mutex/4_thread

"./pagerank inputs/roadNY.dimacs outputs/roadNY/atomic/8_thread 3 8
 - This command will take the roadNY.dimacs file as input, and compute the page rank
   using atomic instructions synchronization (3) with 8 threads and output the result to
   outputs/roadNY/atomic/8_thread