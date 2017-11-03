#include<math.h>
#include<queue>

using namespace std;

typedef struct {
    int  x,y,z;
}feature_t;

typedef struct {
  int num; //The number of features
  feature_t * features;
  float * weights;
}signature_t;



//To store the index and distance from source to every other node
typedef struct node {
   int id;
   float w;
   friend bool operator <(node a, node b){
   return a.w > b.w;
   }
}node;

typedef struct {
   int index;
   float weight;
   float excess;
   float potential;
}vertex;

typedef struct edge {
   //int sIndex;
   int endIndex;
   float capacity;
   float resflow;
   float cost;
   float reducedcost;
   edge * counterEdge;

   edge(){
    this->endIndex = -1;
    this->capacity = 0;
    this->resflow = 0;
    this->cost = 0;
    this->reducedcost =0;
   }

   edge(int dindex, float capacity, float resflow, float cost){
    this->endIndex = dindex;
    this->capacity = capacity;
    this->resflow = resflow;
    this->cost = cost;
    this->reducedcost =cost;
   }
   ~edge(){};
}edge;

 float calTotalFlow(signature_t * s);
 float Dist(feature_t  f1, feature_t  f2);


 float Dist(feature_t  f1, feature_t  f2){
    float dX = f1.x - f2.x;
    float dY = f1.y - f2.y;
    float dZ = f1.z - f2.z;
    return sqrt(dX*dX+dY*dY+dZ*dZ);

}
 float calTotalFlow(signature_t * s){
    float sum =0;
    for(int i =0; i<s->num; i++)
        sum = sum + s->weights[i];
    return sum;
}
