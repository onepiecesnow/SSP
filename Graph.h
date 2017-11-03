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
graph::graph(signature_t *q, signature_t  *p){

    num_vertices = q->num + p->num+2;
    num_edges = 2*(q->num+p->num+q->num * p-> num);
    this->adj = new vector<edge *> [num_vertices];

    this->adj->reserve(num_vertices);

    vertex s,t;
    s.index=0;
    t.index=num_vertices-1;
    s.potential = 0;
    t.potential = 0;
    s.excess = calTotalFlow(q);
    t.excess = -1*calTotalFlow(p);
    s.weight =0;
    t.weight =0;

    this->vertice.push_back(s);

    //init vertex of q
    for(int i =0; i<q->num; i++){

       vertex  tempVertex = newVertex(q->weights[i],i+1);
       this->vertice.push_back(tempVertex);

    }
    //init vertex of p
    for(int i = 0; i<p->num; i++){
       vertex  tempVertex = newVertex(p->weights[i],i+q->num+1);
       this->vertice.push_back(tempVertex);
    }

    this->vertice.push_back(t);
   /* print vertex
    for(int i =0; i<num_vertices; i++){
        cout<<i<<this->vertice[i].index<<"weight"<<vertice[i].weight<<"  excess"<<vertice[i].excess<<endl;
    }*/

    //init edges
    for(int i=0; i<q->num; i++){

        if(vertice[i+1].weight == 0){
             continue;
        }

        //init s-q edges

       edge * tempEdge = new edge(vertice[i+1].index,vertice[i+1].weight,vertice[i+1].weight,0);
       edge *counterEdge = new edge(0,vertice[i+1].weight,0,0);

       tempEdge->counterEdge = counterEdge;
       counterEdge->counterEdge = tempEdge;

       this->adj[0].push_back(tempEdge);
       this->adj[i+1].push_back(counterEdge);



       //init q-p edges
       for(int j =0; j<p->num; j++){

      if(vertice[q->num+j+1].weight==0){
          continue;
       }
       float dist = Dist(q->features[i],p->features[j]);
       edge * tempEdge = new edge(vertice[q->num+j+1].index,min(vertice[i+1].weight,vertice[q->num+j+1].weight),min(vertice[i+1].weight,vertice[q->num+j+1].weight),dist);
       edge *counterEdge = new edge(vertice[i+1].index,min(vertice[i+1].weight,vertice[q->num+j+1].weight),0,-1*dist);

       tempEdge->counterEdge = counterEdge;
       counterEdge->counterEdge = tempEdge;

       this->adj[i+1].push_back(tempEdge);
       this->adj[q->num+j+1].push_back(counterEdge);



       }
    }
        //init p-t edges
    for(int i = 0; i<p->num;i++){

      if(vertice[q->num+i+1].weight == 0){
            continue;
      }


       edge * tempEdge = new edge(t.index,vertice[q->num+i+1].weight,vertice[q->num+i+1].weight,0);
       edge * counterEdge = new edge(vertice[q->num+i+1].index,vertice[q->num+i+1].weight,0,0);

       tempEdge->counterEdge = counterEdge;
       counterEdge->counterEdge = tempEdge;

       this->adj[q->num+i+1].push_back(tempEdge);
       this ->adj[num_vertices-1].push_back(counterEdge);


}
  /*print initial graph g
      for(int i=0; i<num_vertices; i++){
        for(int j=0; j<adj[i].size(); j++){
            cout << i << " " << adj[i][j]->endIndex << " " << adj[i][j]->resflow << " " << adj[i][j]->cost << endl;
        }
    }*/


}


vertex graph:: newVertex(float weight, int index){
    vertex  v ;
    v.index =index ;
    v.weight =weight ;
    v.excess =0 ;
    v.potential = 0;
    return v;
}
graph::~graph(){
    int Size0 = this->adj[0].size();

    for(int j =0; j<Size0;j++){
            delete this->adj[0][j];
            this->adj[0][j] = NULL;


        }

    for(int i =  1; i < num_vertices-1; i++){

        if(this->vertice[i].weight == 0){
            continue;
        }

        int Size = this->adj[i].size();
        for(int j =0; j<Size;j++){
            delete this->adj[i][j];
            this->adj[i][j] = NULL;


        }
    }
     int Sizet = this->adj[num_vertices-1].size();
    for(int j =0; j<Sizet;j++){
            delete this->adj[num_vertices-1][j];
            this->adj[num_vertices-1][j] = NULL;

        }

    // delete this->adj;
    // this->adj = NULL;
     //cout<<"destruction called"<<endl;
}
