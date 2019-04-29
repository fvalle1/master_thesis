//
// Created by Filippo Valle on 2019-04-11.
//

#include <boost/python.hpp>
#include <boost/timer/timer.hpp>
#include <Python.h>
#include "GraphGenerator.h"
#include "MainTable.h"
#include "SamplingModel.h"

void makegraph(){
    boost::timer::auto_cpu_timer stopwatch;
    auto G = new GraphGenerator(3000, 0.9 ,true, true);
    G->MakeGraph();
    delete G;
}

void sampling(int statisticsRepeat=1){
    boost::timer::auto_cpu_timer stopwatch;
    auto model = new SamplingModel();
    model->GenerateNullData(statisticsRepeat);
    model->~SamplingModel();
}

BOOST_PYTHON_FUNCTION_OVERLOADS(sampling_overloads, sampling, 0, 1)

void statistics(bool saveOccurrences=true){
    auto model = new MainTable();
    model->read("mainTable.csv", saveOccurrences);
    model->~MainTable();
}

BOOST_PYTHON_FUNCTION_OVERLOADS(statistics_overloads, statistics, 0, 1)


BOOST_PYTHON_MODULE(graphgenerator_ext)
{
    using namespace boost::python;

    def("statistics", statistics, statistics_overloads(
            boost::python::arg("saveAbundancesOccurrences")=true)
    );
    def("sampling", sampling, sampling_overloads(
            boost::python::arg("averageOver")=1)
    );
    def("makegraph", makegraph);
}