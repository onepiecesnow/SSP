#include<vector>
#include<algorithm>
#include"Util.h"

using namespace std;



class graph{

public:
   int num_vertices;
   int num_edges;


   vector<vertex> vertice;
   vector<edge *> *adj;


   graph(signature_t *q, signature_t  *p);
   ~graph();

   vertex  newVertex(float weight, int index);
   //edge* newEdge(int dindex, float capacity, float resflow, float cost);


};


