#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>

using namespace std;

template <class T>
class Data {
private:
	vector<T> grid;     // one-direction grid point data including the boundary

public:
	inline T& operator()(int i, int j)
	{
		return grid[i*this->ncols() + j];
	}
	int nrows();
	int ncols();
	void grid_resize(int i, int j);
	void grid_fill(double d);
};


template <class T>
int Data<T>::nrows(void)
{
	return sqrt(grid.size());
}

template <class T>
void Data<T>::grid_fill(double d)
{
    fill(grid.begin(),grid.end(),d);
}

template <class T>
int Data<T>::ncols(void)
{
	return sqrt(grid.size());
}

template <class T>
void Data<T>::grid_resize(int i, int j)
{
	grid.resize(i*j);
}

void SetBoundary(int lvl, Data<double> *u)
{
	//left and right boundary (x,0) and (x,1)
	for (int i = 0; i < u[lvl].nrows(); i++)
	{
		u[lvl](i, 0) = 0;
		u[lvl](i, u[lvl].ncols() - 1) = sin(M_PI*i/pow(double(2.0),lvl+1))*sinh(M_PI);
	}
	//up and down boundary (0,y) and (1,y)
	for (int j = 0; j < u[lvl].ncols(); j++)
	{
		u[lvl](0, j) = 0;
		u[lvl](u[lvl].nrows() - 1, j) = 0;
	}
}

void SetBonusBoundary(int lvl, Data<double> *u)
{
	//horizontal boundary (x,0) and (x,1), Dirichlet BC
	for (int i = 0; i < u[lvl].nrows(); i++)
	{
		u[lvl](i, 0) = i/pow(double(2.0),lvl+1)*(1-i/pow(double(2.0),lvl+1));
		u[lvl](i, u[lvl].ncols() - 1) = i/pow(double(2.0),lvl+1)*(1-i/pow(double(2.0),lvl+1));
	}
	//vertical boundary (0,y) and (1,y), Neumann BC
	for (int j = 0; j < u[lvl].ncols(); j++)
	{
		u[lvl](0, j) = u[lvl](1, j) - 1/pow(double(2.0),lvl+1);
		u[lvl](u[lvl].nrows() - 1, j) = u[lvl](u[lvl].nrows() - 2, j) - 1/pow(double(2.0),lvl+1);
	}
}

void RBGauss(int lvl, Data<double> *u, Data<double> *f)
{
	//This is the additional Red-Black Gauss-Seidel Algorithm. Feel free to switch this with 'GaussSeidel' in the V_cycle function.
	for (int i = 1; i < u[lvl].nrows() - 1; i++)
	{
		for (int j = 1; j < u[lvl].ncols() - 1; j++)
		{
			if ((i+j)%2 == 0)
			{
				u[lvl](i, j) = 0.25*(pow(double(2.0), -2 * (lvl+1))*f[lvl](i, j) + u[lvl](i+1, j) + u[lvl](i-1, j) + u[lvl](i, j+1) + u[lvl](i, j-1));
			}
		}
	}
	for (int i = 1; i < u[lvl].nrows() - 1; i++)
	{
		for (int j = 1; j < u[lvl].ncols() - 1; j++)
		{
			if ((i+j)%2 == 1)
			{
				u[lvl](i, j) = 0.25*(pow(double(2.0), -2 * (lvl+1))*f[lvl](i, j) + u[lvl](i+1, j) + u[lvl](i-1, j) + u[lvl](i, j+1) + u[lvl](i, j-1));
			}
		}
	}
}

void GaussSeidel(int lvl, Data<double> *u, Data<double> *f)
{
	//the default smoother
	for (int i = 1; i < u[lvl].nrows() - 1; i++)
	{
		for (int j = 1; j < u[lvl].ncols() - 1; j++)
		{
			u[lvl](i, j) = 0.25*(pow(double(2.0), -2 * (lvl+1))*f[lvl](i, j) + u[lvl](i+1, j) + u[lvl](i-1, j) + u[lvl](i, j+1) + u[lvl](i, j-1));
		}
	}
}

//residual definition: r = f - Lu, here u is the twice-iterated approximation.
#define r(i,j) f[lvl](i,j)+(u[lvl](i+1,j)+u[lvl](i-1,j)+u[lvl](i,j+1)+u[lvl](i,j-1)-4*u[lvl](i,j))/pow(double(2.0),-2*(lvl+1))

void Restrict_residual_from(int lvl, Data<double> *u, Data<double> *f)
{
	//loop over the coarser grid
    for (int i = 1; i < f[lvl-1].nrows() - 1; i++)
	{
        for (int j = 1; j < f[lvl-1].ncols() - 1; j++)
		{
            f[lvl-1](i, j) = double(0.25)*r(2*i, 2*j) + 
			double(0.125)*(r(2*i+1, 2*j) + r(2*i-1, 2*j) + r(2*i, 2*j+1) + r(2*i, 2*j-1)) + 
            double(0.0625)*(r(2*i+1, 2*j+1) + r(2*i-1, 2*j+1) + r(2*i+1, 2*j-1) + r(2*i-1, 2*j-1));
		}
	}
}

void Interpolate_correct_from(int lvl, Data<double> *u)
{
	//loop over the coarser grid
	for (int i = 1; i < u[lvl].nrows(); i++)
	{
		for (int j = 1; j < u[lvl].ncols(); j++)
		{
			u[lvl + 1](2 * i, 2 * j) += u[lvl](i, j);
			u[lvl + 1](2 * i - 1, 2 * j) += double(0.5)*(u[lvl](i, j) + u[lvl](i-1, j));
			u[lvl + 1](2 * i - 1, 2 * j - 1) += double(0.25)*(u[lvl](i, j) + u[lvl](i-1, j) + u[lvl](i, j-1) + u[lvl](i-1, j-1));
			u[lvl + 1](2 * i, 2 * j - 1) += double(0.5)*(u[lvl](i, j) + u[lvl](i, j-1));
		}
	}
}

double residual(int lvl, Data<double> *u, Data<double> *f)
{
    double res = double (0.0);
    double rf;

    for (int i = 1; i<u[lvl].nrows()-1; i++)
    {
        for (int j = 1; j<u[lvl].ncols()-1; j++)
        {
            rf = r(i,j);
            res += rf*rf;
        }
    }
    return sqrt(res/((u[lvl].nrows()-2)*(u[lvl].ncols()-2)));
}

double grid_error(int lvl, Data<double> *u)
{
    double eh = double (0.0);
    double diff;

    for (int i = 1; i<u[lvl].nrows()-1; i++)
    {
        for (int j = 1; j<u[lvl].ncols()-1; j++)
        {
            diff = u[lvl](i,j)-sin(M_PI*i/pow(double(2.0),lvl+1))*sinh(M_PI*j/pow(double(2.0),lvl+1));
            eh += diff*diff;
        }
    }
    return sqrt(eh/((u[lvl].nrows()-2)*(u[lvl].ncols()-2)));
}

void V_cycle(int lvl, Data<double> *u, Data<double> *f)
{
    //solve problem on coarsest grid
    if (lvl == 1)
        GaussSeidel(lvl-1,u,f);
    else
        //or do the cycle recursively
    {
        //pre-smoothing twice
        GaussSeidel(lvl-1,u,f);
        GaussSeidel(lvl-1,u,f);
		
		//restriction
		Restrict_residual_from(lvl-1,u,f);
		
        //initialize the solution of next level to zero
		u[lvl - 2].grid_fill(0);

        V_cycle(lvl-1,u,f);

        //interpolate error and update the solution
        Interpolate_correct_from(lvl-2,u);

        //post-smoothing once
        GaussSeidel(lvl-1,u,f);
    }
}
