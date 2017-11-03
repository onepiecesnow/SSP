#include<iostream>
#include<fstream>
#include<string>
#include<stdio.h>
#include<iomanip>
#include"emd.h"

using namespace std;


int main(){

    int num = 5001;
    int bin = 256;

    //output file
    ofstream outfile("Emd.txt");
    if(!outfile){
        cerr<<"Open EMDoutfile Failed!"<<endl;
        exit(0);
    }
    double avrTime = 0;

   for (int i = 2; i< num; i++){
    int totalFlow1 =0;
    int totalFlow2 =0;

    int j = 0;
    int k = 0;

    feature_t     f1[bin],
                  f2[bin];
    float         w1[bin],
                  w2[bin];

     //input image1 data
        char infilename1[50] ;
		sprintf(infilename1, "C:\\Users\\dell\\Desktop\\EMD\\Corel-LAB\\image_data_%d.txt", i );
        ifstream infile1(infilename1);
        if(!infile1){
           cerr<<"Open infile1 Failed!"<<endl;
           exit(0);
       }

     //input image2 data
        char infilename2[50] ;
		sprintf(infilename2, "C:\\Users\\dell\\Desktop\\EMD\\VOC-LAB\\image_data_%d.txt", i );
        ifstream infile2(infilename2);
        if(!infile2){
           cerr<<"Open infile2 Failed!"<<endl;
           exit(0);
       }


     infile1 >>totalFlow1;
     while(j!=bin){
        infile1>>f1[j].x;
        infile1>>f1[j].y;
        infile1>>f1[j].z;
        //cout<<f1[j].z<<'\t';
        infile1 >>w1[j];
        //w1[j]/=totalFlow1;
        //w1[j] = (int)(w1[j]/totalFlow1 * 1000);
        //cout<<w1[j]<<endl;
         j++;
     }
     infile1.close();


     infile2 >> totalFlow2;
     while(k!=bin){

        infile2>>f2[k].x;
        infile2>>f2[k].y;
        infile2>>f2[k].z;
       //cout<<f2[k].z<<'\t';
        infile2>>w2[k];
       // cout<<w2[k]<<endl;
       // w2[k] =(int)(w2[k]/totalFlow2 * 1000);
        //w2[j]/=totalFlow2;
        k++;
     }
     infile2.close();

       signature_t s1 = {bin,f1,w1},
                   s2 = {bin,f2,w2};

        float e;

        clock_t start, end;
        start = clock();

        e = calEMD(&s1,&s2);
       // e = calSSEEMD(&s1,&s2);
       // e = calAVXEMD(&s1,&s2);

        end = clock();

        avrTime += (double)(end-start)/CLOCKS_PER_SEC;

        cout<<setiosflags(ios::fixed)<<setprecision(4)<<"emd=" <<e<<endl;

       // outfile<<setiosflags(ios::fixed)<<setprecision(4)<<e<<'\n';
        outfile<<e<<'\n';


   }
   cout<<"average Time:"<<avrTime/num<<endl;
   //cout<<"done"<<endl;
    outfile.close();

}

