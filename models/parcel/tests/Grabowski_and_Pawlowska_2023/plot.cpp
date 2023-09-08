#include <gnuplot-iostream.h>

int main()
{
    Gnuplot gp;
    gp << "set term svg size " << 500 << "," << 500 << "\n";
    gp << "set output 'test.svg'"<<"\n";
    gp << "plot sin(x)" << "\n";
}
