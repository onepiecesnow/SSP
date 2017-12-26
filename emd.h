#include<queue>
#include<math.h>
#include"Graph.h"
#include<iostream>
#include<x86intrin.h>
#include<immintrin.h>
#include<xmmintrin.h>

using namespace std;

 graph  initGraph(signature_t * s1, signature_t * s2);//Transform histogram to graph


bool DijSSP(graph * g, int parentVertex[], edge * parentEdge[], node d[]);// Find the shortest path in residual Graph
void SimdDijSSP(graph * g, int parentVertex[], edge * parentEdge[], node d[]);
bool LabelDijSSP(graph * g, int parentVertex[], edge * parentEdge[], node d[],int label[]);

void updateReducedCost(graph * g);//After finding the SSP, reduce the path's reduced cost

void updateNodePot(graph * g, node d[]);//Update node potential
void updateModifiedPot(graph * g,  int parentVertex[], node d[],int label []);

 float augFlow();//Augment the flow of the shortest path

float calEMD(signature_t * Signature1, signature_t * Signature2);
float calSSEEMD(signature_t * Signature1, signature_t * Signature2);
float calAVXEMD(signature_t * Signature1, signature_t * Signature2);



//Using SSE;
void SSE_updateNodePot(graph* g, node d[] );
void SSE_updateReducedCost(graph * g);

//Using AVX;
void AVX_updateNodePot(graph* g, node d[] );
void AVX_updateReducedCost(graph * g);

