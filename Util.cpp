#include "Util.h"


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
