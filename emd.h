#include<queue>
#include<math.h>
#include"Graph.h"

#include<x86intrin.h>
#include<immintrin.h>
#include<xmmintrin.h>


 graph  initGraph(signature_t * s1, signature_t * s2);//Transform histogram to graph


 bool DijSSP(graph * g, int parentVertex[], edge * parentEdge[]);// Find the shortest path in residual Graph

 void updateReducedCost(graph * g);//After finding the SSP, reduce the path's reduced cost
 void updateNodePot(graph * g, node d[]);//Update node potential

 float augFlow();//Augment the flow of the shortest path

float calEMD(signature_t * Signature1, signature_t * Signature2);
float calSSEEMD(signature_t * Signature1, signature_t * Signature2);
float calAVXEMD(signature_t * Signature1, signature_t * Signature2);
//
void updateModifiedPot(graph * g,  int parentVertex[], node d[],int label []);

//Using SSE;
void SSE_updateNodePot(graph* g, node d[] );
void SSE_updateReducedCost(graph * g);

//Using AVX;
void AVX_updateNodePot(graph* g, node d[] );
void AVX_updateReducedCost(graph * g);

/************************Init Graph***************************/
 graph  initGraph(signature_t * s1, signature_t * s2){
    graph   g (s1 , s2);
    return g;
}

/*******************Main Four Steps****************************/
 bool DijSSP(graph * g, int parentVertex[], edge * parentEdge[], node d[]){
    int num_vertices = g->num_vertices; //total number of vertex in graph g
    int source = 0; //index of s
    int sink = num_vertices-1; //index of t
    vector<edge *> * adj = g->adj;
    bool visited[num_vertices];
    //initialize visited, parent vertex and distance
    for(int i=0; i<num_vertices; i++){
        parentVertex[i]=-1;
        d[i].w=1000000;
        d[i].id =i;
        visited[i]=false;
    }
    priority_queue<node> q;
    d[source].w=0;
    q.push(d[source]);
    while(!q.empty()){
        node temp = q.top();

        q.pop();
        int u = temp.id;

        if(visited[u]){

            continue;
        }

        visited[u]=true;
        int edgeSize = adj[u].size();
        for(int e =0; e< edgeSize;e++){

            if(adj[u][e]->resflow >1e-6){
                 int v = adj[u][e]->endIndex;
                 float w = adj[u][e]->reducedcost;
                 //cout<<"d v w"<<d[v].w<<endl;
                 if(d[v].w>d[u].w+w&&!visited[v]){
                     d[v].w=d[u].w+w;
                     parentVertex[v]=u;
                     parentEdge[v]=adj[u][e];
                     q.push(d[v]);
                 }
            }

        }
    }


    if(parentVertex[sink]==-1)
        return false;
    else
        return true;

}
/*************************************LabelDijkstra****************************/
 bool LabelDijSSP(graph * g, int parentVertex[], edge * parentEdge[], node d[],int label[]){
    int num_vertices = g->num_vertices; //total number of vertex in graph g
    int source = 0; //index of s
    int sink = num_vertices-1; //index of t
    vector<edge *> * adj = g->adj;
    bool visited[num_vertices];
    //initialize visited, parent vertex and distance
    for(int i=0; i<num_vertices; i++){
        parentVertex[i]=-1;
        d[i].w=1000000;
        d[i].id =i;
        visited[i]=false;
    }
    priority_queue<node> q;
    d[source].w=0;
    q.push(d[source]);
    while(!q.empty()){
        node temp = q.top();

        q.pop();
        int u = temp.id;

		label[u] = 1;

        if(visited[u]){

            continue;
        }

        visited[u]=true;
        int edgeSize = adj[u].size();
        for(int e =0; e< edgeSize;e++){

            if(adj[u][e]->resflow >1e-6){
                 int v = adj[u][e]->endIndex;
                 float w = adj[u][e]->reducedcost;
                 //cout<<"d v w"<<d[v].w<<endl;
                 if(d[v].w>d[u].w+w&&!visited[v]){
                     d[v].w=d[u].w+w;
                     parentVertex[v]=u;
                     parentEdge[v]=adj[u][e];
                     q.push(d[v]);
                 }
            }

        }
    }


    if(parentVertex[sink]==-1)
        return false;
    else
        return true;

}
/*************************************Potential&Cost**************************/
void updateReducedCost(graph * g){
    int num_vertices = g->num_vertices;
    //Update the cost of every edge
    for(int i=0; i<num_vertices; i++){
            int edgeSize = g->adj[i].size();
            for(int j =0; j< edgeSize;j++){

                g->adj[i][j]->reducedcost = g->adj[i][j]->cost - g->vertice[i].potential + g->vertice[g->adj[i][j]->endIndex].potential;

            }

    }
}
void updateNodePot(graph* g, node d[] ){
    int num_vertices = g->num_vertices;
    for(int i = 0; i< num_vertices;i++){
        g->vertice[i].potential = g->vertice[i].potential -d[i].w ;
    }
}


//A more efficient way to update the node potential
void updateModifiedPot(graph * g,  int parentVertex[], node d[],int label[]){


     int num = g->num_vertices;

     for(int i = 0; i< num;i++){
	 if(label[i]==1)
        g->vertice[i].potential = g->vertice[i].potential - d[i].w + d[num-1].w;

     }

}

/**********************************Augflow****************************************/
float augFlow(graph * g, int parentVertex[], edge * parentEdge[]){
    float path_flow = 1000000;
    float runningCost = 0;
    int s =0;
    int t = g->num_vertices - 1;
    int v;
    edge *tel1,*counterEdge;
    //find the min edge flow of one path
    for(v=t; v!=s; v= parentVertex[v]){

        tel1 = parentEdge[v];

        path_flow = min(path_flow, tel1->resflow);

    }


    path_flow = min(path_flow, g->vertice[s].excess);
   //update the excess of source and sink
    g->vertice[s].excess -=path_flow;
    g->vertice[t].excess +=path_flow;

  //update the residual graph
    for(v=t;v!=s;v = parentVertex[v]){

        tel1 = parentEdge[v];
        counterEdge =  tel1->counterEdge;

        tel1->resflow -= path_flow;

        counterEdge->resflow += path_flow ;

        runningCost += path_flow * (tel1->cost);


    }

    return runningCost;
}
/************************************SSE*********************************/
//Using SSE to accelerate computing;
void SSE_updateReducedCost(graph * g){
    for(int i=0; i<g->num_vertices; i++){
            int num = g->adj[i].size();
            int num4Mod = g->adj[i].size()%4;

            for(int j =0; j<num-num4Mod;j=j+4){

               edge * edgeOne = g->adj[i][j];
               edge * edgeTwo = g->adj[i][j+1];
               edge * edgeThree = g->adj[i][j+2];
               edge * edgeFour = g->adj[i][j+3];



              // __attribute__((aligned(16))) float m_fRcost[4] = {edgeOne->reducedcost,edgeTwo->reducedcost,edgeThree->reducedcost,edgeFour->reducedcost};
               __attribute__((aligned(16))) float m_fCost[4] ={edgeOne->cost,edgeTwo->cost,edgeThree->cost,edgeFour->cost};
               __attribute__((aligned(16))) float m_fEPot[4] ={  g->vertice[edgeOne->endIndex].potential, g->vertice[edgeTwo->endIndex].potential, g->vertice[edgeThree->endIndex].potential, g->vertice[edgeFour->endIndex].potential};
               __attribute__((aligned(16))) float m_Potental[4] = {g->vertice[i].potential,g->vertice[i].potential,g->vertice[i].potential,g->vertice[i].potential};

             /* __declspec(align(16)) float m_fRcost[4] = {g->adj[i][j]->reducedcost,g->adj[i][j+1]->reducedcost,g->adj[i][j+2]->reducedcost,g->adj[i][j+3]->reducedcost};
              __declspec(align(16)) float m_fCost[4] ={g->adj[i][j]->cost,g->adj[i][j+1]->cost,g->adj[i][j+2]->cost,g->adj[i][j+3]->cost};
              __declspec(align(16)) float m_fEPot[4] ={  g->vertice[g->adj[i][j]->endIndex].potential, g->vertice[g->adj[i][j+1]->endIndex].potential, g->vertice[g->adj[i][j+2]->endIndex].potential, g->vertice[g->adj[i][j+3]->endIndex].potential};
              __declspec(align(16)) float m_Potental[4] = {g->vertice[i].potential,g->vertice[i].potential,g->vertice[i].potential,g->vertice[i].potential};*/



                __m128 c_vec, e_vec, p_vec;
                c_vec = _mm_load_ps(m_fCost);
                //r_vec = _mm_load_ps(m_fRcost);
                e_vec = _mm_load_ps(m_fEPot);
                p_vec = _mm_load_ps(m_Potental);

                __m128 temp = _mm_sub_ps(c_vec,  p_vec);
                __m128 result = _mm_add_ps(temp,e_vec);

               // _mm_store_ss(m_fRcost,result);

                g->adj[i][j]->reducedcost = result[0];
                g->adj[i][j+1]->reducedcost = result[1];
                g->adj[i][j+2]->reducedcost = result[2];
                g->adj[i][j+3]->reducedcost = result[3];

            }
            for(int j = num-1; j>num-1-num4Mod ;  j--){
                 g->adj[i][j]->reducedcost = g->adj[i][j]->cost - g->vertice[i].potential + g->vertice[g->adj[i][j]->endIndex].potential;

            }

    }
}

void SSE_updateNodePot(graph* g, node d[] ){
    int numSize = g->num_vertices;
    int num4Mod = numSize%4;

    for(int i = 0; i<numSize-num4Mod;i= i+4){
            //align node d[] and  node potential
           __attribute__((aligned(16))) float m_Dist[4] = {d[i].w,d[i+1].w,d[i+2].w,d[i+3].w} ;
           __attribute__((aligned(16))) float m_Pot[4] = {g->vertice[i].potential,g->vertice[i+1].potential,g->vertice[i+2].potential,g->vertice[i+3].potential};

         __m128 p_vec,d_vec;
         p_vec = _mm_load_ps(m_Pot);
         d_vec = _mm_load_ps(m_Dist);

         __m128 result = _mm_sub_ps(p_vec,d_vec);

        // _mm_store_ps(m_Pot,result);

        g->vertice[i].potential = result[0];
        g->vertice[i+1].potential = result[1];
        g->vertice[i+2].potential = result[2];
        g->vertice[i+3].potential = result[3];
    }

    for(int j = numSize-1; j> numSize-1-num4Mod; j--){
         g->vertice[j].potential = g->vertice[j].potential -d[j].w ;

    }

}
/********************************AVX************************************/

//Using AVX to accelerate computing;
/*void AVX_updateReducedCost(graph * g){
    for(int i=0; i<g->num_vertices; i++){
            int num = g->adj[i].size();
            int num4Mod = g->adj[i].size()%8;

            for(int j =0; j<num-num4Mod;j=j+8){

               edge * edgeOne = g->adj[i][j];
               edge * edgeTwo = g->adj[i][j+1];
               edge * edgeThree = g->adj[i][j+2];
               edge * edgeFour = g->adj[i][j+3];

               edge * edgeFive = g->adj[i][j+4];
               edge * edgeSix = g->adj[i][j+5];
               edge * edgeSeven = g->adj[i][j+6];
               edge * edgeEight = g->adj[i][j+7];


               __attribute__((aligned(32))) float m_fRcost[8] = {edgeOne->reducedcost,edgeTwo->reducedcost,edgeThree->reducedcost,edgeFour->reducedcost,
                                                                 edgeFive->reducedcost,edgeSix->reducedcost,edgeSeven->reducedcost,edgeEight->reducedcost   };
               __attribute__((aligned(32))) float m_fCost[8] ={edgeOne->cost,edgeTwo->cost,edgeThree->cost,edgeFour->cost,
                                                                edgeFive->cost,edgeSix->cost,edgeSeven->cost,edgeEight->cost };
               __attribute__((aligned(32))) float m_fEPot[8] ={  g->vertice[edgeOne->endIndex].potential, g->vertice[edgeTwo->endIndex].potential, g->vertice[edgeThree->endIndex].potential, g->vertice[edgeFour->endIndex].potential,
                                                                 g->vertice[edgeFive->endIndex].potential, g->vertice[edgeSix->endIndex].potential, g->vertice[edgeSeven->endIndex].potential, g->vertice[edgeEight->endIndex].potential};
              // __attribute__((aligned(32))) float m_Potental[8] = {g->vertice[i].potential,g->vertice[i].potential,g->vertice[i].potential,g->vertice[i].potential,};
              __attribute__((aligned(32))) float m_Potental[8] = {g->vertice[i].potential};



                __m256 c_vec, e_vec, p_vec;
                c_vec = _mm256_load_ps(m_fCost);
                //r_vec = _mm_load_ps(m_fRcost);
                e_vec = _mm256_load_ps(m_fEPot);
                p_vec = _mm256_load_ps(m_Potental);
                __m256 temp = _mm256_sub_ps(c_vec,  p_vec);
                __m256 result = _mm256_add_ps(temp,e_vec);

               // _mm_store_ss(m_fRcost,result);

                g->adj[i][j]->reducedcost = result[0];
                g->adj[i][j+1]->reducedcost = result[1];
                g->adj[i][j+2]->reducedcost = result[2];
                g->adj[i][j+3]->reducedcost = result[3];

                g->adj[i][j+4]->reducedcost = result[4];
                g->adj[i][j+5]->reducedcost = result[5];
                g->adj[i][j+6]->reducedcost = result[6];
                g->adj[i][j+7]->reducedcost = result[7];
            }
            for(int j = num-1; j>num-1-num4Mod ;  j--){
                 g->adj[i][j]->reducedcost = g->adj[i][j]->cost - g->vertice[i].potential + g->vertice[g->adj[i][j]->endIndex].potential;

            }

    }
}

void AVX_updateNodePot(graph* g, node d[] ){
    int numSize = g->num_vertices;
    int num4Mod = numSize%8;

    for(int i = 0; i<numSize-num4Mod;i= i+8){
           __attribute__((aligned(32))) float m_Dist[8] = {d[i].w,d[i+1].w,d[i+2].w,d[i+3].w,d[i+4].w,d[i+5].w,d[i+6].w,d[i+7].w} ;
           __attribute__((aligned(32))) float m_Pot[8] = {g->vertice[i].potential,g->vertice[i+1].potential,g->vertice[i+2].potential,g->vertice[i+3].potential,
                                                          g->vertice[i+4].potential,g->vertice[i+5].potential,g->vertice[i+6].potential,g->vertice[i+7].potential};

         __m256 p_vec,d_vec;
         p_vec = _mm256_load_ps(m_Pot);
         d_vec = _mm256_load_ps(m_Dist);

         __m256 result = _mm256_sub_ps(p_vec,d_vec);

        // _mm_store_ps(m_Pot,result);

        g->vertice[i].potential = result[0];
        g->vertice[i+1].potential = result[1];
        g->vertice[i+2].potential = result[2];
        g->vertice[i+3].potential = result[3];

        g->vertice[i+4].potential = result[4];
        g->vertice[i+5].potential = result[5];
        g->vertice[i+6].potential = result[6];
        g->vertice[i+7].potential = result[7];
    }

    for(int j = numSize-1; j> numSize-1-num4Mod; j--){
         g->vertice[j].potential = g->vertice[j].potential -d[j].w ;

    }

}*/


/************************************EMD**********************************/
float calEMD(signature_t * Signature1, signature_t * Signature2)
{
       graph   g =  initGraph( Signature1, Signature2);

       int num = g.num_vertices;
       int parentVertex[num];

       edge * parentEdge[num];
       node d[num];
       int label[num];

       float totalCost =0;
       int source = 0;

       float norm = min(g.vertice[source].excess,-1*g.vertice[num-1].excess);
       g.vertice[source].excess = norm;

       while(g.vertice[source].excess > 1e-6){

          DijSSP(&g,parentVertex,parentEdge,d);

       // LabelDijSSP(&g,parentVertex,parentEdge,d,label);
        /*print shortest path
           for(int i =0;i<g.num_vertices;i++ ){
            cout<<"index :   "<<i<<"  "<<parentVertex[i]<<endl;
           }*/

          updateNodePot(&g,d);

        //  updateModifiedPot(&g,parentVertex, d, label);


           float runningCost = augFlow(&g, parentVertex, parentEdge);


         updateReducedCost(&g);

           totalCost += runningCost;


       }

       return totalCost/norm;
}
/***************************SSEEMD*******************************/

float calSSEEMD(signature_t * Signature1, signature_t * Signature2)
{
       graph   g =  initGraph( Signature1, Signature2);

       int num = g.num_vertices;
       int parentVertex[num];

       edge * parentEdge[num];
       node d[num];


       float totalCost =0;
       int source = 0;

       float norm = min(g.vertice[source].excess,-1*g.vertice[num-1].excess);
       g.vertice[source].excess = norm;

       while(g.vertice[source].excess > 1e-6){

        DijSSP(&g,parentVertex,parentEdge,d);

        /*print shortest path
           for(int i =0;i<g.num_vertices;i++ ){
            cout<<"index :   "<<i<<"  "<<parentVertex[i]<<endl;
           }*/

        SSE_updateNodePot(&g,d);



        float runningCost = augFlow(&g, parentVertex, parentEdge);


        SSE_updateReducedCost(&g);

        totalCost += runningCost;


       }

       return totalCost/norm;
}

/*float calAVXEMD(signature_t * Signature1, signature_t * Signature2){

       graph   g =  initGraph( Signature1, Signature2);

       int num = g.num_vertices;
       int parentVertex[num];
       //float norm = g->vertice[0].excess;
       edge * parentEdge[num];
       node d[num];
       int label[num]={0};

       float totalCost =0;
       int source = 0;

       float norm = min(g.vertice[source].excess,-1*g.vertice[num-1].excess);
       g.vertice[source].excess = norm;

       while(g.vertice[source].excess > 1e-6){

        DijSSP(&g,parentVertex,parentEdge,d);
        //  LabelDijSSP(&g,parentVertex,parentEdge,d,label);


        AVX_updateNodePot(&g,d);


         // updateModifiedPot(&g,parentVertex, d);


        float runningCost = augFlow(&g, parentVertex, parentEdge);


        AVX_updateReducedCost(&g);

        totalCost += runningCost;


       }

       return totalCost/norm;
}*/



