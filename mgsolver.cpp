#include <iostream>
#include "Data.h"
#include <fstream>
#include <sys/time.h>

using namespace std;

int main(int argc, char** argv){

    if (argc != 3)
	{
		cout << "Input: (number of levels) (number of cycles)" << endl;
		exit(0);
	}

	int nlvls = atoi(argv[1]);
	int ncycles = atoi(argv[2]);

	struct timeval tim;

    //define array to store u and f on different levels
	Data<double> *u = new Data<double>[nlvls];
    Data<double> *f = new Data<double>[nlvls];

    //memory allocation on each level
	for (int i = 0; i < nlvls; i++)
	{
		u[i].grid_resize(pow(2,i+1)+1,pow(2,i+1)+1);
		f[i].grid_resize(pow(2,i+1)+1,pow(2,i+1)+1);
	}

    int i = 1;
    double res = 0, tempres;

	gettimeofday(&tim, NULL);
	double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);

	
    while (i <= ncycles)
    {
        tempres = res;
        //SetBoundary(nlvls-1,u);
        //For bonus task, comment out the 'SetBoundary' function above and use the following:
        SetBonusBoundary(nlvls-1,u);
        f[nlvls-1].grid_fill(2);
        V_cycle(nlvls,u,f);
        res = residual(nlvls-1,u,f);
        cout << "After V-cycle " << i << ":" << endl;
        cout << "L2-norm of residual is: " << res << endl;
        if (i != 1)
        {
            cout << "Convergence rate 'q': " << double(res/tempres) << endl;
        }		
        i++;
    }
	
	
	gettimeofday(&tim, NULL);
	double t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
		
	//task 4: print out the time and write the solution to a txt file
	cout << "Time elapsed (seconds): " << t2 - t1 << endl;
	
	ofstream SaveFile("solution.txt");
	for (int k = 0; k < u[nlvls-1].ncols(); k++)
    {
		for (int j = 0; j < u[nlvls-1].nrows(); j++)
		{
			SaveFile << double(k)/(u[nlvls-1].ncols()-1) << "\t" << double(j)/(u[nlvls-1].nrows()-1) << "\t" << u[nlvls-1](k,j) << endl;
		}
		SaveFile << endl;
	}
	SaveFile.close();

	//Task 5: print out the error on different level of grid spacing
	if (nlvls >= 3)
	{
		for (int q = 3; q <= nlvls; q++)
		{
			cout << "L2 norm of error e(h) for grid size 1/" << pow(2,q) << " is: " << grid_error(q-1,u) << endl;
		}
	}

	//free the memory
	delete[] u;
    delete[] f;
    
	return 0;

}
