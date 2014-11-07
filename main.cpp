#include <iostream>
#include <ctime>
#include <mpi/mpi.h>

#define RESULT_TAG 100501
#define START_TAG 100502
#define END_TAG 100503

#define INT_STEP 0.000001

using namespace std;

double func (double x) {
    return x*x;
}

//формула левых прямоугольников
double calcInt(double start, double end){
    double res = 0;
    for(start += INT_STEP; start <= end; start += INT_STEP)
        res += INT_STEP*func(start);
    return res;
}

void mainProc(const int procCount){
    double start, end;
    cout << "Enter start and end of the interval" << endl;
    cin >> start;
    cin >> end;

    clock_t startTime = clock();
    MPI_Status status;

    double cur_a, cur_b;
    double step = (end - start) / procCount;
    cur_a = start + step;

    for (int i = 1; i < procCount; ++i) {
            cur_b = cur_a + step;
            MPI_Send(&cur_a, 1, MPI_DOUBLE, i, START_TAG, MPI_COMM_WORLD);
            MPI_Send(&cur_b, 1, MPI_DOUBLE, i, END_TAG, MPI_COMM_WORLD);
            cur_a += step;
    }

    double res = calcInt(start, start + step);

    for (int i = 1; i < procCount; ++i){
        double tmp;
        MPI_Recv(&tmp, 1, MPI_DOUBLE, i, RESULT_TAG, MPI_COMM_WORLD, &status);
        res += tmp;
    }

    time_t endTime = clock();
    cout.flush();
    cout << "Time used: " << (endTime - startTime) << " ms." << endl;
    cout << "Result: " << res << endl;
}

void slaveProc(){

    MPI_Status status;
    double a, b;
    MPI_Recv(&a, 1, MPI_DOUBLE, 0, START_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&b, 1, MPI_DOUBLE, 0, END_TAG, MPI_COMM_WORLD, &status);

    double res = calcInt(a, b);

    MPI_Send(&res, 1, MPI_DOUBLE, 0, RESULT_TAG, MPI_COMM_WORLD);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int procRank, procCount;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procCount);

    if (procRank){
        slaveProc();
    } else {
        mainProc(procCount);
    }

    MPI_Finalize();
    return 0;
}