//
// Created by Filippo Valle on 2019-02-28.
//

#ifndef THESIS_DATA_MINING_GRAPHGENERATOR_H
#define THESIS_DATA_MINING_GRAPHGENERATOR_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

using namespace std;
using boost::property_tree::ptree;
typedef uint64_t checkable;

class GraphGenerator {
public:
    void MakeGraph(uint64_t maxStorableDocs = 100);

private:
    void addKeyAttrs(ptree &graphml) const;
    void addWordNode(ptree &graph, string name, uint64_t id) const;
    void addDocumentNode(ptree &graph, string title, uint64_t id) const;
    void addEdge(ptree &graph, uint64_t idSource, uint64_t idTarget, uint64_t weight) const;
    std::vector<std::string> tokenize(const std::string &);

    bool testWord(string word, checkable parameter);
};


#endif //THESIS_DATA_MINING_GRAPHGENERATOR_H
