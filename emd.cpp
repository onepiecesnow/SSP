#include "emd.h"


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
/*************************************SIMDDijkstra***************************/
void SimdDijSSP(graph * g, int parentVertex[], edge * parentEdge[], node d[]){
    int num_vertices = g->num_vertices;
    int num_subvertices = (num_vertices-2)/2;
    int source = 0;
    int sink = num_vertices-1;
    vector<edge *> * adj = g->adj;
   // bool visited[num_vertices];
    for(int i=0; i<num_vertices; i++){
        parentVertex[i]=-1;
        d[i].w=1000000;
        d[i].id =i;
        }

    d[source].w = 0;


    int S_size = adj[source].size();
    int S_num8Mod = adj[source].size()%8;

    //Find shortest distance from s to u
    for(int i = 0; i < S_size-S_num8Mod; i=i+8){

        __m256 d_vec, w_vec, f_vec, res_vec, pv_vec, p_vec, s_vec;
        w_vec = _mm256_setr_ps(adj[source][i]->reducedcost,adj[source][i+1]->reducedcost,adj[source][i+2]->reducedcost,adj[source][i+3]->reducedcost,
                               adj[source][i+4]->reducedcost,adj[source][i+5]->reducedcost,adj[source][i+6]->reducedcost,adj[source][i+7]->reducedcost);

        pv_vec =_mm256_set1_ps(source);

        d_vec = _mm256_setr_ps (d[i+1].w,d[i+2].w,d[i+3].w,d[i+4].w,d[i+5].w,d[i+6].w,d[i+7].w,d[i+8].w);

        f_vec = _mm256_setr_ps(adj[source][i]->resflow,adj[source][i+1]->resflow,adj[source][i+2]->resflow,adj[source][i+3]->resflow,
                               adj[source][i+4]->resflow,adj[source][i+5]->resflow,adj[source][i+6]->resflow,adj[source][i+7]->resflow);
        s_vec = _mm256_set1_ps(d[source].w);


        res_vec = _mm256_set1_ps(1e-6);

        p_vec = _mm256_setr_ps(  parentVertex[i+1], parentVertex[i+2], parentVertex[i+3], parentVertex[i+4],
                                 parentVertex[i+5], parentVertex[i+6], parentVertex[i+7], parentVertex[i+8]);

        __m256 v_vec, min_vec,add_vec;
        add_vec = _mm256_add_ps(w_vec,s_vec);
        min_vec = _mm256_min_ps(add_vec, d_vec);

        __m256 mask1;
        mask1 = _mm256_cmp_ps(f_vec,res_vec,_CMP_GT_OQ);


        v_vec = _mm256_blendv_ps(d_vec,min_vec,mask1);



       d[i+1].w= v_vec[0]; d[i+2].w= v_vec[1]; d[i+3].w= v_vec[2]; d[i+4].w= v_vec[3];
       d[i+5].w= v_vec[4]; d[i+6].w= v_vec[5]; d[i+7].w= v_vec[6]; d[i+8].w= v_vec[7];


        __m256 mask2;
        mask2 = _mm256_cmp_ps(d_vec,add_vec, _CMP_GT_OQ);

        __m256 mask3;
        mask3 = _mm256_and_ps(mask1, mask2);




        __m256 upv_vec;
        upv_vec = _mm256_blendv_ps(p_vec, pv_vec, mask3);



        parentVertex[i+1] = upv_vec[0]; parentVertex[i+2] = upv_vec[1]; parentVertex[i+3] = upv_vec[2]; parentVertex[i+4] = upv_vec[3];
        parentVertex[i+5] = upv_vec[4]; parentVertex[i+6] = upv_vec[5]; parentVertex[i+7] = upv_vec[6]; parentVertex[i+8] = upv_vec[7];


     for(int k = 0;k < 8; k++){

       if(mask3[k]!=0){
            int v = adj[source][i+k]->endIndex;

           parentEdge[v] = adj[source][i+k];
          }
       }

    }

    for ( int j = S_size -1; j>S_size-1-S_num8Mod;j--){

        if(adj[source][j]->resflow > 1e-6){
                int v = adj[source][j]->endIndex;
                 float w = adj[source][j]->reducedcost;

                 if(d[v].w>d[source].w+w){
                     d[v].w=d[source].w+w;
                     parentVertex[v]= source;
                     parentEdge[v]= adj[source][j];
                 }
         }
    }

 while(parentVertex[sink]== -1){


   for(int i = 1; i<num_subvertices+1;i++){
                      int uv_size = adj[i].size();
                      int U_num8Mod = uv_size%8;
                     //Find shortest distance from s to u

                     for(int j = 0; j < uv_size - U_num8Mod; j=j+8){
                      __m256 d_vec, w_vec, f_vec, u_vec,  res_vec, p_vec;
                      w_vec = _mm256_setr_ps (adj[i][j]->reducedcost,adj[i][j+1]->reducedcost,adj[i][j+2]->reducedcost,adj[i][j+3]->reducedcost,
                                              adj[i][j+4]->reducedcost,adj[i][j+5]->reducedcost,adj[i][j+6]->reducedcost,adj[i][j+7]->reducedcost);

                      int ld1 = adj[i][j]->endIndex; int ld2 = adj[i][j+1]->endIndex; int ld3 = adj[i][j+2]->endIndex; int ld4 =adj[i][j+3]->endIndex ;
                      int ld5 = adj[i][j+4]->endIndex; int ld6 = adj[i][j+5]->endIndex; int ld7 = adj[i][j+6]->endIndex; int ld8 = adj[i][j+7]->endIndex;

                      d_vec = _mm256_setr_ps (d[ld1].w,d[ld2].w,d[ld3].w,d[ld4].w,
                                              d[ld5].w,d[ld6].w,d[ld7].w,d[ld8].w);

                      u_vec = _mm256_set1_ps(d[i].w);

                      f_vec = _mm256_setr_ps(adj[i][j]->resflow,adj[i][j+1]->resflow,adj[i][j+2]->resflow,adj[i][j+3]->resflow,
                               adj[i][j+4]->resflow,adj[i][j+5]->resflow,adj[i][j+6]->resflow,adj[i][j+7]->resflow);


                     // p_vec = _mm256_set1_ps(-1);
                     p_vec = _mm256_setr_ps( parentVertex[ld1], parentVertex[ld2], parentVertex[ld3], parentVertex[ld4],
                                             parentVertex[ld5], parentVertex[ld6], parentVertex[ld7], parentVertex[ld8]);


                      res_vec = _mm256_set1_ps(1e-6);

                      __m256 v_vec,add_vec,min_vec;
                      add_vec = _mm256_add_ps(w_vec,u_vec);
                      min_vec = _mm256_min_ps(add_vec,d_vec);

                       __m256 mask1;
                       mask1 = _mm256_cmp_ps(f_vec,res_vec,_CMP_GT_OQ);
                       v_vec = _mm256_blendv_ps(d_vec,min_vec,mask1);

                      d[ld1].w= v_vec[0]; d[ld2].w= v_vec[1]; d[ld3].w= v_vec[2]; d[ld4].w= v_vec[3];
                      d[ld5].w= v_vec[4]; d[ld6].w= v_vec[5]; d[ld7].w= v_vec[6]; d[ld8].w= v_vec[7];




                       __m256 mask2;
                       mask2 = _mm256_cmp_ps(d_vec,add_vec, _CMP_GT_OQ);
                       __m256 mask3;
                       mask3 = _mm256_and_ps(mask1, mask2);

                       __m256 pv_vec;
                       pv_vec =_mm256_set1_ps(i);

                       __m256 upv_vec;
                       upv_vec = _mm256_blendv_ps(p_vec, pv_vec, mask3);

                      parentVertex[ld1] = upv_vec[0]; parentVertex[ld2] = upv_vec[1]; parentVertex[ld3] = upv_vec[2]; parentVertex[ld4] = upv_vec[3];
                      parentVertex[ld5] = upv_vec[4]; parentVertex[ld6] = upv_vec[5]; parentVertex[ld7] = upv_vec[6]; parentVertex[ld8] = upv_vec[7];

                      for(int k = 0; k < 8; k++){
                        if(mask3[k]!=0){
                            int v = adj[i][j+k]->endIndex;
                            parentEdge[v] = adj[i][j+k];
                        }
                      }

                }

        for ( int k = uv_size -1; k>uv_size-1-U_num8Mod;k--){

             if(adj[i][k]->resflow > 1e-6){
                      int v = adj[i][k]->endIndex;
                      float w = adj[i][k]->reducedcost;
                     if(d[v].w>d[i].w+w){
                          d[v].w=d[i].w+w;
                          parentVertex[v]=i;
                          parentEdge[v]=adj[i][k];
                       }
                 }
              }
          }

     for(int i = num_subvertices+1; i<num_vertices-1;i++){

                      int v_size = adj[i].size();
                      int  V_num8Mod = v_size%8;

                     //Find shortest distance from s to u
                     for(int j = 0; j < v_size - V_num8Mod; j=j+8){
                      __m256 d_vec, w_vec, f_vec, u_vec, res_vec, p_vec;
                      w_vec = _mm256_setr_ps (adj[i][j]->reducedcost,adj[i][j+1]->reducedcost,adj[i][j+2]->reducedcost,adj[i][j+3]->reducedcost,
                                              adj[i][j+4]->reducedcost,adj[i][j+5]->reducedcost,adj[i][j+6]->reducedcost,adj[i][j+7]->reducedcost);

                      int ld1 = adj[i][j]->endIndex; int ld2 = adj[i][j+1]->endIndex; int ld3 = adj[i][j+2]->endIndex; int ld4 =adj[i][j+3]->endIndex ;
                      int ld5 = adj[i][j+4]->endIndex; int ld6 = adj[i][j+5]->endIndex; int ld7 = adj[i][j+6]->endIndex; int ld8 = adj[i][j+7]->endIndex;

                      d_vec = _mm256_setr_ps (d[ld1].w,d[ld2].w,d[ld3].w,d[ld4].w,
                                              d[ld5].w,d[ld6].w,d[ld7].w,d[ld8].w);
                      u_vec = _mm256_set1_ps(d[i].w);

                      f_vec = _mm256_setr_ps(adj[i][j]->resflow,adj[i][j+1]->resflow,adj[i][j+2]->resflow,adj[i][j+3]->resflow,
                               adj[i][j+4]->resflow,adj[i][j+5]->resflow,adj[i][j+6]->resflow,adj[i][j+7]->resflow);


                     res_vec = _mm256_set1_ps(1e-6);

                     p_vec = _mm256_setr_ps( parentVertex[ld1], parentVertex[ld2], parentVertex[ld3], parentVertex[ld4],
                                             parentVertex[ld5], parentVertex[ld6], parentVertex[ld7], parentVertex[ld8]);

                      __m256 v_vec,add_vec, min_vec;
                      add_vec = _mm256_add_ps(w_vec,u_vec);
                      min_vec = _mm256_min_ps(add_vec,d_vec);

                       __m256 mask1;
                       mask1 = _mm256_cmp_ps(f_vec,res_vec,_CMP_GT_OQ);
                       v_vec = _mm256_blendv_ps(d_vec,min_vec,mask1);

                      d[ld1].w= v_vec[0]; d[ld2].w= v_vec[1]; d[ld3].w= v_vec[2]; d[ld4].w= v_vec[3];
                      d[ld5].w= v_vec[4]; d[ld6].w= v_vec[5]; d[ld7].w= v_vec[6]; d[ld8].w= v_vec[7];

                       __m256 mask2;
                       mask2 = _mm256_cmp_ps(d_vec,add_vec, _CMP_GT_OQ);
                       __m256 mask3;
                       mask3 = _mm256_and_ps(mask1, mask2);

                       __m256 pv_vec;
                       pv_vec =_mm256_set1_ps(i);


                       __m256 upv_vec;
                       upv_vec = _mm256_blendv_ps(p_vec, pv_vec, mask3);

                      parentVertex[ld1] = upv_vec[0]; parentVertex[ld2] = upv_vec[1]; parentVertex[ld3] = upv_vec[2]; parentVertex[ld4] = upv_vec[3];
                      parentVertex[ld5] = upv_vec[4]; parentVertex[ld6] = upv_vec[5]; parentVertex[ld7] = upv_vec[6]; parentVertex[ld8] = upv_vec[7];

                     for(int k = 0; k < 8; k++){
                        if(mask3[k]!=0){
                            int v = adj[i][j+k]->endIndex;
                            parentEdge[v] = adj[i][j+k];
                         }
                      }

                }

        for ( int k = v_size -1; k>v_size-1-V_num8Mod;k--){

             if(adj[i][k]->resflow > 1e-6){
                      int v = adj[i][k]->endIndex;
                      float w = adj[i][k]->reducedcost;

                     if(d[v].w>d[i].w+w){
                          d[v].w=d[i].w+w;

                          parentVertex[v]=i;
                          parentEdge[v]=adj[i][k];
                       }
                 }
              }


        }

     }
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


                __m128 c_vec, e_vec, p_vec;

                c_vec = _mm_setr_ps(g->adj[i][j]->cost,g->adj[i][j+1]->cost,g->adj[i][j+2]->cost,g->adj[i][j+3]->cost);
                e_vec = _mm_setr_ps(g->vertice[g->adj[i][j]->endIndex].potential, g->vertice[g->adj[i][j+1]->endIndex].potential, g->vertice[g->adj[i][j+2]->endIndex].potential, g->vertice[g->adj[i][j+3]->endIndex].potential);
                p_vec = _mm_set1_ps(g->vertice[i].potential);

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


         __m128 p_vec,d_vec;

        d_vec = _mm_setr_ps(d[i].w,d[i+1].w,d[i+2].w,d[i+3].w);
        p_vec = _mm_setr_ps(g->vertice[i].potential,g->vertice[i+1].potential,g->vertice[i+2].potential,g->vertice[i+3].potential);

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
void AVX_updateReducedCost(graph * g){
    for(int i=0; i<g->num_vertices; i++){
            int num = g->adj[i].size();
            int num4Mod = g->adj[i].size()%8;

            for(int j =0; j<num-num4Mod;j=j+8){



                __m256 c_vec, e_vec, p_vec;

                c_vec = _mm256_setr_ps(g->adj[i][j]->cost,g->adj[i][j+1]->cost,g->adj[i][j+2]->cost,g->adj[i][j+3]->cost,
                                       g->adj[i][j+4]->cost,g->adj[i][j+5]->cost,g->adj[i][j+6]->cost,g->adj[i][j+7]->cost);
                p_vec = _mm256_set1_ps(g->vertice[i].potential);
                e_vec = _mm256_setr_ps(g->vertice[g->adj[i][j]->endIndex].potential, g->vertice[g->adj[i][j+1]->endIndex].potential, g->vertice[g->adj[i][j+2]->endIndex].potential, g->vertice[g->adj[i][j+3]->endIndex].potential,
                                       g->vertice[g->adj[i][j+4]->endIndex].potential, g->vertice[g->adj[i][j+5]->endIndex].potential, g->vertice[g->adj[i][j+6]->endIndex].potential, g->vertice[g->adj[i][j+7]->endIndex].potential);

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


         __m256 p_vec,d_vec;

         d_vec =  _mm256_setr_ps(d[i].w,d[i+1].w,d[i+2].w,d[i+3].w,d[i+4].w,d[i+5].w,d[i+6].w,d[i+7].w);
         p_vec = _mm256_setr_ps(g->vertice[i].potential,g->vertice[i+1].potential,g->vertice[i+2].potential,g->vertice[i+3].potential,g->vertice[i+4].potential,g->vertice[i+5].potential,g->vertice[i+6].potential,g->vertice[i+7].potential);
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

}


/************************************EMD**********************************/
float calEMD(signature_t * Signature1, signature_t * Signature2)
{
       graph   g =  initGraph( Signature1, Signature2);

       int num = g.num_vertices;
       int parentVertex[num];

       edge * parentEdge[num];
       node d[num];
      // int label[num];

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

float calAVXEMD(signature_t * Signature1, signature_t * Signature2){

       graph   g =  initGraph( Signature1, Signature2);

       int num = g.num_vertices;
       int parentVertex[num];
       //float norm = g->vertice[0].excess;
       edge * parentEdge[num];
       node d[num];
      // int label[num]={0};

       float totalCost =0;
       int source = 0;

       float norm = min(g.vertice[source].excess,-1*g.vertice[num-1].excess);
       g.vertice[source].excess = norm;

       while(g.vertice[source].excess > 1e-6){

       SimdDijSSP(&g,parentVertex,parentEdge,d);
        //  LabelDijSSP(&g,parentVertex,parentEdge,d,label);

        AVX_updateNodePot(&g,d);




        float runningCost = augFlow(&g, parentVertex, parentEdge);


        AVX_updateReducedCost(&g);

        totalCost += runningCost;


       }

       return totalCost/norm;
}

