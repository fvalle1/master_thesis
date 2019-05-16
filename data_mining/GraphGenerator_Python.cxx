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
    auto G = new GraphGenerator(5000, 1 ,false, true);
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

void statistics(bool saveOccurrences=true, bool considerZeros=true){
    auto model = new MainTable();
    model->read("mainTable.csv", saveOccurrences);
    model->SaveMeansVariances("mainTable.csv", considerZeros);
    model->~MainTable();
}

BOOST_PYTHON_FUNCTION_OVERLOADS(statistics_overloads, statistics, 0, 2)


BOOST_PYTHON_MODULE(tacos)
{
    using namespace boost::python;

    def("statistics", statistics, statistics_overloads(
            (
                    boost::python::arg("saveAbundancesOccurrences")=true,
                            boost::python::arg("considerZeros")=true)
        )
    );
    def("sampling", sampling, sampling_overloads(
            boost::python::arg("averageOver")=1)
    );
    def("makegraph", makegraph);
}
