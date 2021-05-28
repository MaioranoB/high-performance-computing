/* A pure C/C++ version of a Gauss-Siedel Laplacian solver to test the
   speed of a C program versus that of doing it with
   Python/Numeric/Weave. */

#include <iostream>
#include <cmath>
#include <time.h>
#include <omp.h>

typedef double Real;

inline Real SQR(const Real &x)
{
    return (x*x);
}

inline Real BC(Real x, Real y)
{
    return (x*x - y*y);
}

struct Grid {
    Real dx, dy;
    int nx, ny;
    Real **u;
    
    Grid(const int n_x=10, const int n_y=10);
    ~Grid();
    
    void setBCFunc(Real (*f)(const Real, const Real));
    void resetGrid(const int n_x, const int n_y);
    /*Real computeError();*/
};

Grid :: Grid(const int n_x, const int n_y) : nx(n_x), ny(n_y)
{
    dx = 1.0/Real(nx - 1);
    dy = 1.0/Real(ny - 1);

    u = new Real* [nx];
    for (int i=0; i<nx; ++i) {
        u[i] = new double [ny];
    }

    for (int i=0; i<nx; ++i) {
        for (int j=0; j<ny; ++j) {
            u[i][j] = 0.0;
        }
    }
    
}

Grid :: ~Grid()
{
    for (int i=0; i<nx; ++i) {
        delete [] u[i];
    }
    delete [] u;
}

void Grid :: resetGrid(const int n_x, const int n_y){
    for (int i=0; i<n_x; ++i) {
        for (int j=0; j<n_y; ++j) {
            u[i][j] = 0.0;
        }
    }
}

void Grid :: setBCFunc(Real (*f)(const Real, const Real))
{
    Real xmin, ymin, xmax, ymax, x, y;
    xmin = 0.0;
    ymin = 0.0;
    xmax = 1.0;
    ymax = 1.0;
    /* Left and right sides. */
    for (int j=0; j<ny; ++j) {
        y = j*dy;
        u[0][j] = f(xmin, y);
        u[nx-1][j] = f(xmax, y);
    }
    /* Top and bottom sides. */
    for (int i=0; i<nx; ++i) {
        x = i*dx;
        u[i][0] = f(x, ymin);
        u[i][ny-1] = f(x, ymax);
    }
}

struct LaplaceSolver{
    Grid *g;

    LaplaceSolver(Grid *g);
    ~LaplaceSolver();
    void initialize();
    Real timeStep(const Real dt=0.0);
    Real solve(const int n_iter=0, const Real eps=1e-16);
};

LaplaceSolver :: LaplaceSolver(Grid *grid)
{
    g = grid;
    initialize();
}
LaplaceSolver:: ~LaplaceSolver()
{
}

void LaplaceSolver :: initialize()
{
}

Real LaplaceSolver :: timeStep(const Real dt)
{
    Real dx2 = g->dx*g->dx;
    Real dy2 = g->dy*g->dy;
    Real tmp;
    Real err = 0.0;
    int nx = g->nx;
    int ny = g->ny;
    int i,j;
    Real **u = g->u;
    Real c = 0.5/(dx2 + dy2);

    #pragma omp parallel for shared(u) private(i,j,tmp) reduction(+:err)
    for (i=1; i<nx-1; ++i) {
        for (j=1; j<ny-1; ++j) {
            tmp = u[i][j];
            u[i][j] = ((u[i-1][j] + u[i+1][j])*dy2 + (u[i][j-1] + u[i][j+1])*dx2)*c;
            // err += SQR(u[i][j] - tmp);
            err += (u[i][j] - tmp) * (u[i][j] - tmp);            
        }
    }
    return sqrt(err);
}

Real LaplaceSolver :: solve(const int n_iter, const Real eps)
{
    Real err = timeStep();
    int count = 1;
    while (err > eps) {
        if (n_iter && (count >= n_iter)) {
            return err;
        }
        err = timeStep();
        ++count;        
    }
    return Real(count);
}

int main(int argc, char * argv[])
{
    int nx, n_iter,max_threads;
    Real eps, result;
    Real t_start, t_end;

    nx = atoi(argv[1]);
    n_iter = 1000;
    eps = 10e-16;

    Grid *g = new Grid(nx, nx);
    g->setBCFunc(BC);
    LaplaceSolver s = LaplaceSolver(g);
    
    max_threads = omp_get_max_threads();
    for (int n_threads=1; n_threads <= max_threads; n_threads++){
        omp_set_num_threads(n_threads);
        

        t_start = omp_get_wtime();
        result = s.solve(n_iter, eps);
        t_end = omp_get_wtime();

        g->resetGrid(nx,nx);
        g->setBCFunc(BC);

        std::cout
            <<argv[0]<<","<<nx<<","<<n_threads<<","<<t_end - t_start<<","<<result<<std::endl;
    }
    return 0;
}
