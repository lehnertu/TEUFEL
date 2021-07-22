/*
 * 
 * compile with :
 *     mpic++ -fopenmp hybrid.c -o hybrid
 * 
 * run with :
 *     export OMP_NUM_THREADS=4
 *     mpirun -np 2 -x OMP_NUM_THREADS ./hybrid
 */

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <mpi.h>
#include <omp.h>

int main(int argc, char *argv[]) {
    int numprocs, rank, namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int iam = 0, np = 1;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processor_name, &namelen);

    if (rank==0)
    {
        std::string line;
        std::string delimiter = ":";
        std::ifstream myfile;
        myfile.open("/proc/meminfo");
        if(!myfile.is_open()) {
            std::cout << "cannot access /proc/meminfo - memory size unknown." << std::endl;
        } else {
            while(getline(myfile, line)) {
                // std::cout << line << std::endl;
                size_t pos = line.find(delimiter);
                std::string name = line.substr(0, pos);
                std::string rest = line.erase(0,pos+1);
                if (name=="MemTotal") {
                    int size_kB = std::stoi(rest);
                    std::cout << " available memory ";
                    std::cout << size_kB/1024 << " MB" << std::endl;
                }
            }
        }        
    }
    
    #pragma omp parallel default(shared) private(iam, np)
    {
        np = omp_get_num_threads();
        iam = omp_get_thread_num();
        printf("Hello from thread %d out of %d from process %d out of %d on %s\n",
            iam, np, rank, numprocs, processor_name);
    }

    MPI_Finalize();
}
