#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <conio.h>
#include <fstream>
#include "flowrate.h"
#include <iomanip>
#include <algorithm>
#include "apvector.h"

using namespace std;

int main (){


    int N;
    float L,H,VIS,DEN,SG,Qi,CF;
    //double C,Q,D,D1,Dh,Pi,Ah,Dr,Cd,As,Vs,Area,nDia;
    double C,Q,D,D1,Dh,Ah,Dr,As,Vs,Area,nDia;
    const double Pi=(22.0/7.0);
    std::cin.clear();

    begin1:
    std::cout<<"Enter a value for Length of Sparger in(m) and less than 30m"<<endl;
    std::cin>>L;
    if(L>30){
        std::cout<<"Enter a value Less than 30"<<endl;
        goto begin1;
    }

    begin2:
    std::cout<<"Enter a value for Diameter of Sparger in(mm) and less than 1500 mm"<<endl;
    std::cin>>D1;
    D=D1/1000;
    if(D1>1500){
        std::cout<<"Enter a value Less than 1500"<<endl;
        goto begin2;
    }

    begin3:
    std::cout<<"Enter the Number of holes in Sparger"<<endl;
    std::cin>>N;
    if(N>1000){
        std::cout<<"Enter a value Less than 1000"<<endl;
        goto begin3;
    }

    begin4:
    std::cout<<"Enter the Head of liquid above Sparger in (m)  and less than 50m"<<endl;
    std::cin>>H;
    if(H>50){
        std::cout<<"Enter a value Less than 50"<<endl;
        goto begin4;
    }

    begin5:
    std::cout<<"Enter the Flow Rate of Air in cms"<<endl;
    std::cin>>Q;
    if(Q>5000){
        std::cout<<"Check the value again. Enter a value Less than 5000"<<endl;
        goto begin5;
    }
       /*Defined Constants*/
    VIS = 0.0000173; //viscocity
    DEN = 1.225; //Density
    SG = 1.3; //Specific Gravity
    C  = 0.593; //Constant

    /*Write to File*/
    string filename = "D:/desktop-files-amitray/Desktop_folders/Dissolver-FRFCF/NEW-Dissolver/recent/Dissolver_Latest_Drawings&DOCS/doc-revisions-a/Dissolver-HT-Documents/sparger-design/Sparger-C++/Sparger/data.txt";
    ofstream outFile(filename.c_str(), std::ios_base::app | std::ios_base::out);

    outFile<<"The value for Length of Sparger is : "<<L<<"m"<<endl;
    outFile<<"Value for Diameter of Sparger  is : "<<D<<"m"<<endl;
    outFile<<"Number of holes in Sparger is : "<<N<<endl;
    outFile<<"Head of liquid above Sparger is : "<<H<<endl;
    outFile<<"The Flow rate of air in Sparger is : "<<Q<<"cms"<<endl;
    outFile<<"The Viscocity of air used in Sparger is : "<<VIS<<"Pa.s"<<endl;
    outFile<<"The Density of air used in Sparger is : "<<DEN<<"kg/m3"<<endl;
    outFile<<"The Specific Gravity of Fluid above sparger is : "<<SG<<endl;
    outFile<<"The Discharge Coefficient of sparger hole is : "<<C<<endl;


    Qi = Q/N;
    std::cout<<"Flow per Orifice is : "<<Qi<<"cms"<<endl;
    Ah=Qi/(C*(sqrt(2*9.81*H)));
    outFile<<"The Area of each hole is : "<<Ah<<"sqm"<<endl;

    Dh = sqrt(Ah*4/Pi);
    outFile<<"The Diameter of each hole is : "<<Dh*1000<<"mm"<<endl;
    Dr = (Dh/D);
    std::cout<<"The Diameter ratio is : "<<Dr<<endl;
    /*if(Dr<0.5){
        Cd = 0.593;
    } else {
        Cd = 0.623;
    }
    outFile<<"The Calculated Discharge Coefficient is : "<<Cd<<endl;
    if(C==Cd){
        std::cout<<"The Assumed Discharge Coefficient is correct."<<endl;
    } else {
        std::cout<<"Assumed Discharge Coefficient is wrong. Change the C value and recalculate"<<endl;
    }*/
    /*Sparger Area Calculation*/
    As = Pi*pow(D,2)/4.0;

     /*Vs Initial Velocity*/
    Vs = Q/As;

    outFile<<"Area of Sparger is : "<<As<<endl;
    double ReynoldsNo;
    flowrate flow;
    outFile<<"Value of N is : "<<N<<endl;
    double S[N];
    S[0] = 0.5*L/N;
    for(int i=1;i<N;i++){
        S[i] = L/N;
    }
    /*std::cout<<"S[1] ="<<S[1]<<"S[2]="<<S[2]<<"S[3]="<<S[3]<<"S[4]="<<S[4]<<"S[5]="<<S[5]<<endl;
    std::cout<<"S[6] ="<<S[6]<<"S[7]="<<S[7]<<"S[8]="<<S[8]<<"S[9]="<<S[9]<<"S[10]="<<S[10]<<endl;
    std::cout<<"S[11] ="<<S[11]<<"S[12]="<<S[12]<<"S[13]="<<S[13]<<"S[14]="<<S[14]<<"S[0]="<<S[0]<<endl;*/

    std::cout<<"Qi is : "<<Qi<<endl;
    double Qn[N];
    Qn[0]=Q;
    for(int j=1;j<N;j++){
        Qn[j] = Q-(j*Qi);
    }
    double Vi[N];
    for(int k=0;k<=N;k++){
        Vi[k] = Qn[k]/As;
    }
    /*double Vj[N];
    Vj[0]=0;
    for(int k=1;k<N;k++){
        Vj[k] = Qn[k]/As;
    }*/
    double IntFactor[N];
    double fri;
    double HL[N];
    double Rootterm[N]; //Energy Balance Term
    double Qin[N];
    double QSum;
    double SumPrLoss=0;
    double rough=0.0018; //roughness in inch converted from 0.015 mm
    double LTemp=0;
    for(int x=0;x<N;x++){
        /*Reynolds Number*/
        ReynoldsNo = flow.Reynolds(D,Vi[x],DEN,VIS); //Get Reynolds number from Velocity
        fri = flow.Friction(ReynoldsNo,rough,(D*33.79)); // get friction factor from reynolds number, roughness factor and diameter
        //double fri2;
        //fri2 = flow.Frictionold(ReynoldsNo);
        cout<<"Reynolds Number"<<ReynoldsNo<<"  roughness: "<<rough<<"  Diameter: "<<D<<endl;
        //cout<<"OLD Friction factor "<<fri2<<endl;
        //cout<<"NEW Friction factor "<<fri<<endl;
        //Get length of travel

        for(int y=0;y<=x;y++){
            LTemp = S[y]+LTemp;
        }

        double Q = Qn[x]*2118.88 ; //m3/s to ft3/min
        double density = DEN*0.062428 ; //kg/m3 to lb/ft3
        double LEN = LTemp*3.28084; //m to feet
        double dia = D*39.3701; //m to inch
        //get the friction loss in line from the obtained friction factor
        IntFactor[x] = flow.FrictionLoss(fri,density,LEN,Q,dia);//Line pressure drop
        //IntFactor[x] = S[x]*fri[x]*pow(Vj[x],2);
        HL[x]=flow.HLength(IntFactor[x]); //Get the head loss in m; convert psi to m
        //cout<<"The Head Loss for "<<x<<" is "<<HL[x]<<endl;

        Rootterm[x]=H+((pow(Vs,2))/(2*9.81))-(pow(Vi[x],2)/(2*9.81))-HL[x]; //head
        SumPrLoss = SumPrLoss + HL[x];
        //cout<<"The Pressure Loss for "<<x<<" is "<<SumPrLoss<<endl;
        Area = Ah;
        //Qin[x]=Ah*C*sqrt(2*9.81*Rootterm[x]);
        Qin[x]=flow.FlowInd(Area,C,Rootterm[x]);

    }
    outFile<<"Total Pressure Drop in 'm' is "<<SumPrLoss<<endl;
    /*std::cout<<"First Value"<<endl;
    std::cout<<"Qin[1] ="<<Qin[1]<<"Qin[2]="<<Qin[2]<<"Qin[3]="<<Qin[3]<<"Qin[4]="<<Qin[4]<<"Qin[5]="<<Qin[5]<<endl;
    std::cout<<"First Value Ends"<<endl;*/
    QSum=Qin[0];
    for(int y=1;y<N;y++){
        QSum=QSum+Qin[y];
    }
    //std::cout<<"Rootterm[1] ="<<Rootterm[1]<<"Rootterm[2]="<<Rootterm[2]<<"Rootterm[3]="<<Rootterm[3]<<"Rootterm[4]="<<Rootterm[4]<<"Rootterm[5]="<<Rootterm[5]<<endl;
    std::cout<<"Total Calculated Flow rate is : "<<QSum<<"cms"<<endl;
    CF = Q/QSum;
    std::cout<<"Hole Area correction Factor is : "<<CF<<endl;
    int iter = 0;
    /*Calculation to get difference*/

    while(CF < 0.98 || CF > 1){
        cout<<"Initial Area is: "<<Ah<<endl;
        Area = CF*Ah;
        double *ptr1,Qint[N];
        ptr1 = &Ah;
        *ptr1 = Area;
        cout<<"Change in Area is: "<<Ah<<endl;
        cout<<"CF Value : "<<CF<<endl;
        for(int x=0;x<N;x++){
            Qint[x]=flow.FlowInd(Area,C,Rootterm[x]);
            Qin[x]=Qint[x];
        }
        QSum=Qint[0];
        for(int y=1;y<N;y++){
            QSum=QSum+Qint[y];
        }
        CF = Q/QSum;
        iter++;
    }
    cout<<"CF final value is "<<CF<<endl;
    double maxval = Qin[0], minval = Qin[0];
    for(int i = 0; i<N; i++)
    {
          if(Qin[i] > maxval){
                maxval = Qin[i];
          }
    }
    for(int i = 0; i<N; i++)
    {
          if(Qin[i] < minval){
                minval = Qin[i];
          }
          //outFile<<"minimum value after iteration "<<i<<" is "<<minval<<"cms"<<endl;
    }
    outFile<<"maximum is "<<maxval<<"cms"<<endl;
    outFile<<"Minimum is "<<minval<<"cms"<<endl;
    double diffperc = flow.getDiff(maxval, minval);//Calculate percentage difference between flows
    outFile<<"Difference in max and min orifice flow rate "<<diffperc<<"%"<<endl;
    nDia = sqrt(Area*4/Pi);
    outFile<<"No of iterations is : "<<iter<<endl;
    outFile<<"Final CF Value : "<<CF<<endl;
    //outFile<<"Qin[1] ="<<Qin[1]<<"Qin[2]="<<Qin[2]<<"Qin[3]="<<Qin[3]<<"Qin[4]="<<Qin[4]<<"Qin[5]="<<Qin[5]<<endl;
    outFile<<"New Hole Area is : "<<setprecision(3)<<Area<<"sqm"<<endl;
    outFile<<"New Diameter is : "<<setprecision(3)<<nDia*1000<<"mm"<<endl;
    outFile<<"Total Calculated Final Flow rate is : "<<QSum<<"cms"<<endl;

    double NL[N];
    for(int y=0;y<N;y++){
        NL[y]=(Qin[y]/Q)*L;
    }
    //std::cout<<"NL[0] ="<<NL[0]<<"NL[1] ="<<NL[1]<<"NL[2]="<<NL[2]<<"NL[3]="<<NL[3]<<"NL[4]="<<NL[4]<<"NL[5]="<<NL[5]<<endl;
    double SL[N];
    SL[0]=0.5*NL[0];
    for(int y=1;y<N;y++){
        SL[y]=0.5*(NL[y-1]+NL[y]);
    }
    /*outFile<<"The head is is -"<<endl;
    for(int i=0;i<N;i++){
        outFile<<"No-"<<i+1<<" = "<<setprecision(3)<<Rootterm[i]<<"m"<<endl;
    }*/

    outFile<<"Compiled Data is"<<endl;
        outFile<<"Hole Sl No "<<"\t"<<"Hole Distance"<<"\t"<<"Flow rate cms"<<"\t"<<endl;
    for(int i=0;i<N;i++){
        outFile<<"No-"<<i+1<<"\t"<<setprecision(3)<<SL[i]<<"\t"<<setprecision(3)<<Qin[i]<<"\t"<<endl;
    }
    //std::cout<<"SL[1] ="<<setprecision(3)<<SL[1]<<"SL[2] ="<<setprecision(3)<<SL[2]<<"SL[3] ="<<setprecision(3)<<SL[3]<<"SL[4] ="<<setprecision(3)<<SL[4]<<"SL[5] ="<<setprecision(3)<<SL[5]<<"SL[6] ="<<setprecision(3)<<SL[6]<<endl;
    double QS=SL[1];
    for(int y=2;y<=N;y++){
        QS=QS+SL[y];
    }
    std::cout<<"The final Length for confirmation is : "<<QS<<"m"<<endl;

    outFile.close();

    if (!outFile)
    {
            cout << "Error opening file" << endl;
            return -1;
    }

    /*while(!inFile.eof())
    {
        cout << "Error opening file" << endl;
    }*/
    outFile<<"\n";
    outFile<<"\n";
    outFile<<"\n";
}
