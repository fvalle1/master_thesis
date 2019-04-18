//
// Created by Filippo Valle on 2019-04-11.
//

#include <boost/python.hpp>
#include <boost/timer/timer.hpp>
#include <Python.h>
#include "GraphGenerator.h"


void makegraph(){
    boost::timer::auto_cpu_timer stopwatch;
    auto G = new GraphGenerator(3000, 0.9 ,false, false);
    G->MakeGraph();
    delete G;
}

BOOST_PYTHON_MODULE(graphgenerator_ext)
{
    using namespace boost::python;

    def("makegraph", makegraph);
}
