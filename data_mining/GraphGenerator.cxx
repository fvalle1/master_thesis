//
// Created by Filippo Valle on 2019-02-28.
//

#include "GraphGenerator.h"


GraphGenerator::GraphGenerator(uint64_t maxStorableDocs, bool weighted, bool countPropertyNode):
        fCounts(countPropertyNode),
        fWeighted(weighted),
        fMaxStorableDocs(maxStorableDocs){
    if(!fWeighted) fCounts=true;
}

void GraphGenerator::MakeGraph() {
    cout << "Making graph.xml" << endl;
    fstream graphxmlfile("graph.xml", std::ios::out);
    fstream tableFile("mainTable.csv", std::ios::in);
    fstream geneInfoFile("meanVariances.csv", std::ios::in);
    if (!tableFile.is_open() || !geneInfoFile.is_open()) {
        cerr << "some files not found" << endl;
        return;
    }

    //structure for nodes
    typedef pair<string, int64_t> TitleId;
    typedef pair<int64_t, checkable> IdCheckable;
    typedef pair<uint64_t,pair<string,string>> Weight_WordTitle;
    map<string, int64_t> titles; //file_id, node
    map<string, IdCheckable> words;
    vector<Weight_WordTitle> edges;
    string line;
    line.reserve(1000000);

    cout << "loading words and documents.." << endl;

    //read header line and generate titles.txt
    //NB begins()+1 because first row is just "gene" word
    std::string header;
    vector<string> FullFilesList;
    if (getline(tableFile, header).good()) {
        FullFilesList = tokenize(header);
        if (FullFilesList.size() < fMaxStorableDocs + 1) fMaxStorableDocs = FullFilesList.size() - 1;

        cout << "acceptable documents: " << fMaxStorableDocs << endl;
        std::for_each(FullFilesList.begin() + 1, FullFilesList.begin() + 1 + fMaxStorableDocs, [&](string name) {
            titles.insert(TitleId(name.substr(0, 36), -2)); //TCGA file-id has 36 chars
        });
    } else {
        cerr << "error parsing header" << endl;
        return;
    }

    /*
    * Load gene (aka words) infos to put condition in the next passage
    */
    //Reads a file with ENSG,mean,type_of_gene,variance,occurrence

    //trash out header
    getline(geneInfoFile, line);
    while (getline(geneInfoFile, line).good()) {
        auto tokenizedLine = tokenize(line);
        auto genename = *(tokenizedLine.begin());
        auto geneocc = std::stoull(*(tokenizedLine.begin() + 4));
        words.insert(pair<string, IdCheckable>(genename.substr(0, 15), IdCheckable(-2, geneocc))); //ENSG has 15 characters
    }
    geneInfoFile.close();

    cout << "computing graph.." << endl;
    uint64_t currentDoc = 1;
    while (getline(tableFile, line).good()) {
        auto tokenizedline = tokenize(line);
        auto genename = (tokenizedline.begin())->substr(0, 15);
        printf("\r%s", genename.c_str());
        auto geneIt = words.find(genename);
        if (geneIt != words.end()) {//if have info about word..can look for edges
            std::for_each(tokenizedline.begin() + 1, tokenizedline.begin() + 1 + fMaxStorableDocs,
                          [&](string fpkm_string) {
                              auto currentRead = std::stoull(fpkm_string); //weight
                              //condition to work on
                              if (geneIt->second.second < 0.25 * fMaxStorableDocs) { //check checkable condition
                                  if (currentRead >= 1) {
                                      geneIt->second.first = -1; //add to nodes
                                      auto currentDocTitle = FullFilesList[currentDoc].substr(0, 36); //TCGA file_name has 36 chars
                                      auto currentDocIterator = titles.find(currentDocTitle);
                                      if (currentDocIterator != titles.end()) {
                                          currentDocIterator->second = -1; //add to nodes
                                          edges.emplace_back(Weight_WordTitle(currentRead,
                                                                              pair<string, string>(genename,
                                                                                                   currentDocIterator->first)));
                                      }
                                  }
                              }

                              currentDoc++;
                          }
            );
        }
        currentDoc = 1;
    }

    tableFile.close();
    line.clear();

    printf("\rcomputing completed..\n");
    cout<<"adding nodes to graph.."<<endl;

    ptree xmlstructure;
    ptree graphml;
    addKeyAttrs(graphml);

    ptree graph;
    graph.put("<xmlattr>.id", "G");
    graph.put("<xmlattr>.edgedefault", "undirected");
    graph.put(boost::property_tree::ptree::path_type("<xmlattr>|parse.nodeids",'|'),"canonical");
    graph.put(boost::property_tree::ptree::path_type("<xmlattr>|parse.edgeids",'|'),"canonical");
    graph.put(boost::property_tree::ptree::path_type("<xmlattr>|parse.order",'|'),"nodesfirst");

    long long nodescount = 0;

    /*
    * Add documents nodes
    */
    cout<<"adding documents.."<<endl;
    BOOST_FOREACH(string title, titles | boost::adaptors::map_keys){
                    auto cDoc = titles.find(title);
                    if (cDoc->second==-1){
                        cDoc->second=nodescount;
                        addDocumentNode(graph, title, nodescount++);
                    }
                }

    /*
    * Add word nodes
    */
    cout<<"adding words.."<<endl;
    BOOST_FOREACH(string gene, words | boost::adaptors::map_keys){
                    auto geneInfo = words.find(gene);
                    if(geneInfo->second.first==-1) {//have to add as node
                        geneInfo->second.first=nodescount;
                        addWordNode(graph, gene, nodescount++);
                    }
                }
    /*
    * Add edges
    */

    cout<<"adding edges.."<<endl;
    BOOST_FOREACH(Weight_WordTitle e,edges){
                    auto wordIndex = words.find(e.second.first)->second.first;
                    auto docIndex = titles.find(e.second.second)->second;
                    if(fCounts){
                        addEdge(graph, docIndex, wordIndex, e.first);
                    }else{
                        for(uint64_t w = 0; w < e.first; w++) addEdge(graph, docIndex, wordIndex, e.first);
                    }
                }

    //clean vectors
    titles.clear();
    words.clear();
    edges.clear();

    graphml.add_child("graph", graph);
    xmlstructure.add_child("graphml", graphml);

    cout<<"writing xml file.."<<endl;
    write_xml(graphxmlfile, xmlstructure, boost::property_tree::xml_writer_make_settings<std::string>(' ', 4));
    graphxmlfile.close();

    //compress jast generate file
    system("rm -f graph.xml.gz");
    system("gzip graph.xml");
    cout<<"graph.xml.gz ready..."<<endl;

    cout<<endl;
}

void GraphGenerator::addEdge(ptree &graph, uint64_t idSource, uint64_t idTarget, uint64_t weight) const {
    static uint64_t intEdgeId = 0;
    std::ostringstream edgename, startname, targetname, weightname;
    startname<<"n"<<idSource;
    targetname<<"n"<<idTarget;
    edgename<<"e"<<intEdgeId++;
    weightname<<weight;

    ptree edgedata;
    if(fCounts) {
        edgedata.put("<xmlattr>.key", "key0");
        edgedata.put("", weightname.str());
    }

    ptree edge;
    edge.put("<xmlattr>.id", edgename.str());
    edge.put("<xmlattr>.source", startname.str());
    edge.put("<xmlattr>.target", targetname.str());
    edge.add_child("data", edgedata);
    graph.add_child("edge", edge);
}

void GraphGenerator::addDocumentNode(ptree &graph, string title, uint64_t id) const {
    ptree docdata;
    docdata.put("<xmlattr>.key", "key1");
    docdata.put("", "0");

    std::ostringstream nodeid;
    ptree node;
    nodeid << "n" << id;

    node.put("<xmlattr>.id", nodeid.str());
    node.add_child("data", docdata);

    ptree nodedata;
    nodedata.put("", title);
    nodedata.put("<xmlattr>.key", "key2");
    node.add_child("data", nodedata);

    graph.add_child("node", node);
}

void GraphGenerator::addWordNode(ptree &graph, string name, uint64_t id) const {
    ptree worddata;
    worddata.put("<xmlattr>.key", "key1");
    worddata.put("", "1");

    ptree node;
    std::ostringstream nodeid;
    nodeid << "n" << id;

    node.put("<xmlattr>.id", nodeid.str());
    node.add_child("data", worddata);

    ptree nodedata;
    nodedata.put("", name);
    nodedata.put("<xmlattr>.key", "key2");
    node.add_child("data", nodedata);

    graph.add_child("node", node);
}

void GraphGenerator::addKeyAttrs(ptree &graphml) const {
    graphml.put("<xmlattr>.xmlns", "http://graphml.graphdrawing.org/xmlns");
    graphml.put("<xmlattr>.xmlns:xsi","http://www.w3.org/2001/XMLSchema-instance");
    graphml.put("<xmlattr>.xsi:schemaLocation","http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd");

    ptree key0, key1, key2;
    key0.put("<xmlattr>.id", "key0");
    key0.put("<xmlattr>.for", "edge");
    key0.put(boost::property_tree::ptree::path_type("<xmlattr>|attr.name", '|'), "count");
    key0.put(boost::property_tree::ptree::path_type("<xmlattr>|attr.type", '|'), "int");
    graphml.add_child("key", key0);

    key1.put("<xmlattr>.id", "key1");
    key1.put("<xmlattr>.for", "node");
    key1.put(boost::property_tree::ptree::path_type("<xmlattr>|attr.name", '|'), "kind");
    key1.put(boost::property_tree::ptree::path_type("<xmlattr>|attr.type", '|'), "int");
    graphml.add_child("key", key1);

    key2.put("<xmlattr>.id", "key2");
    key2.put("<xmlattr>.for", "node");
    key2.put(boost::property_tree::ptree::path_type("<xmlattr>|attr.name", '|'), "name");
    key2.put(boost::property_tree::ptree::path_type("<xmlattr>|attr.type", '|'), "string");
    graphml.add_child("key", key2);
}


std::vector<std::string> GraphGenerator::tokenize( const std::string& line )
{
    // escape char is \ , fields are seperated by , , some fields may be quoted with
    boost::escaped_list_separator<char> sep( '\\', ',', '\0' ) ;
    boost::tokenizer< boost::escaped_list_separator<char> > tokenizer( line, sep ) ;
    return std::vector<std::string>(tokenizer.begin(), tokenizer.end()) ;
}

