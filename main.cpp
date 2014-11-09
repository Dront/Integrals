#include <iostream>
#include <ctime>
#include <mpi/mpi.h>


#define START_TAG_Y 100500
#define END_TAG_Y 100501
#define START_TAG_X 100502
#define END_TAG_X 100503

#define INT_STEP 0.00002

using namespace std;

static double result;

//static double func(const double x) {
//    return x;
//}

static double func(const double x, const double y){
    return x * x * y * y;
}


//формула симпсона
//double simpson(const double a, const double b) {
//    long stepCount = (long)((b - a) / INT_STEP);
//    if (stepCount % 2){
//        stepCount++;
//    }
//
//    double I2 = 0, I4 = 0;
//    for(long k = 1; k < stepCount; k ++ )
//    {
//        if (k % 2){
//            I2 += func(a + k * INT_STEP);
//        } else {
//            I4 += func(a + (k + 1) * INT_STEP);
//        }
//    }
//    return INT_STEP / 3 * (func(a) + func(b) + 4 * I4 + 2 * I2);
//}

double simpson(const double a, const double b, const double x) {
    long stepCount = (long)((b - a) / INT_STEP);
    if (stepCount % 2){
        stepCount++;
    }

    double I2 = 0, I4 = 0;
    for(long k = 1; k < stepCount; k ++ )
    {
        if (k % 2){
            I2 += func(x, a + k * INT_STEP);
        } else {
            I4 += func(x, a + (k + 1) * INT_STEP);
        }
    }
    return INT_STEP / 3 * (func(x, a) + func(x, b) + 4 * I4 + 2 * I2);
}

double doubleSimpson(const double aX, const double bX, const double aY, const double bY){
    long stepCountX = (long)((bX - aX) / INT_STEP);
    if (stepCountX % 2){
        stepCountX++;
    }

    //проходим по x
    double I2 = 0, I4 = 0;
    for(long k = 1; k < stepCountX; k ++ )
    {

        //суммируем функции, в каждой проходим обычным симпсоном по y
        if (k % 2){
            I2 += simpson(aY, bY, aX + k * INT_STEP);
        } else {
            I4 += simpson(aY, bY, aX + (k + 1) * INT_STEP);
        }
    }
    double res = simpson(aY, bY, aX) + simpson(aY, bY, bX) + 4 * I4 + 2 * I2;
    res *= INT_STEP / 3;

    return res;
}


void mainProc(const int procCount){
    double startX, endX;
    cout << "Enter interval for X: ";
    cin >> startX;
    cin >> endX;

    double startY, endY;
    cout << "Enter interval for Y: ";
    cin >> startY;
    cin >> endY;

    double cur_a, cur_b;
    double step = (endX - startX) / procCount;
    cur_a = startX + step;

    time_t startTime = time(NULL);

    for (int i = 1; i < procCount; ++i) {
        cur_b = cur_a + step;
        MPI_Send(&cur_a, 1, MPI_DOUBLE, i, START_TAG_X, MPI_COMM_WORLD);
        MPI_Send(&cur_b, 1, MPI_DOUBLE, i, END_TAG_X, MPI_COMM_WORLD);
        MPI_Send(&startY, 1, MPI_DOUBLE, i, START_TAG_Y, MPI_COMM_WORLD);
        MPI_Send(&endY, 1, MPI_DOUBLE, i, END_TAG_Y, MPI_COMM_WORLD);
        cur_a += step;
    }

    //double res = Rectangles(startX, startX + step);
    double res = doubleSimpson(startX, startX + step, startY, endY);

    MPI_Reduce(&res, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    time_t endTime = time(NULL);

    cout.flush();
    cout << "Time used: " << (endTime - startTime) << " sec." << endl;
    cout << "Result: " << result << endl;
}

void slaveProc(){

    MPI_Status status;
    double aX, bX, aY, bY;
    MPI_Recv(&aX, 1, MPI_DOUBLE, 0, START_TAG_X, MPI_COMM_WORLD, &status);
    MPI_Recv(&bX, 1, MPI_DOUBLE, 0, END_TAG_X, MPI_COMM_WORLD, &status);
    MPI_Recv(&aY, 1, MPI_DOUBLE, 0, START_TAG_Y, MPI_COMM_WORLD, &status);
    MPI_Recv(&bY, 1, MPI_DOUBLE, 0, END_TAG_Y, MPI_COMM_WORLD, &status);

    double res = doubleSimpson(aX, bX, aY, bY);

    MPI_Reduce(&res, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
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