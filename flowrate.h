#ifndef FLOWRATE_H
#define FLOWRATE_H
#include "apvector.h"

class flowrate
{
    //double Ps,Ga,Vs,Vi,Hl,Sj,Fr,As;
    double f,Re,d,v,r,m;
    public:
    /*Function to Calculate Reynolds number*/
    double Reynolds(double d,double v,double r,double m){
        return (d*v*r)/m;
    }

    /*Function to Calculate Friction Factor from Reynolds NumberReynolds number*/
    double Friction(double R, double e,double d){
        //For Gas flow
        if(R>2000){
            f = 0.25/(pow((log10((e/3.7*d)+(5.74/(pow(R,0.9))))),2));//swamee and Jain
        } else {
            f= 64/R;
        }
        //f=(1.0/(-1.8*log10(6.9/Re)));
        return f;
    }
    double Frictionold(double Re){
        //For Gas flow
        //f = 0.25/(pow((log10(e/3.7*d+5.74/(pow(R,0.9)))),2));
        if(Re>2000){
            f=(1.0/(-1.8*log10(6.9/Re)));
        } else {
            f = 64/Re;
        }
        return f;
    }
    double FrictionLoss(double f,double den,double l,double q, double d){
        double fl = (f*den*l*(pow(q,2)))/(82.76*(pow(d,5)));
        return fl;
    }
    double HLength(double FL){
        return FL*0.70307; //convert psi to meters of water column
    }

    double FlowInd(double a,double c,double x){
        return a*c*sqrt(2*9.81*x);
    }

    /*double MaxValue(apvector<double> &array)
    {
         int length = array.length( );  // establish size of array
         double maxval = array[0];       // start with max = first element

         for(int i = 1; i<length; i++)
         {
              if(array[i] > maxval){
                    maxval = array[i];
              }
         }
         return maxval;                // return highest value in array
    }

    double MinValue(apvector<double> &array)
    {
         int length = array.length( );  // establish size of array
         double minval = array[0];       // start with max = first element

         for(int i = 1; i<length; i++)
         {
              if(array[i] < minval){
                    minval = array[i];
              }
         }
         return minval;                // return highest value in array
    }*/

    double getDiff(double &maxval, double &minval){
        double diff = (maxval - minval);
        double diffperc = (diff/maxval)*100;
        return diffperc;
    }
    /*array WriteFile(char filename, data)
    {

    	//ofstream outFile  ("H:\\workspace\\Windfarms2.csv");
    	ofstream outFile  (filename);
    	outFile << "Wind Farm Data, Wind Speed, Power Generated" << endl;
        for(int j=0; j<count(data); j++){
            for(int i=0; i<count(headers); i++){
                outFile << data[j][i] <<","<< endl;
            }
    	}
    	outFile.close();

    	if (!outFile)
    	{
                cout << "Error opening file" << endl;
    			return -1;
    	}

    	while(!inFile.eof())
    	{
    		cout << "Error opening file" << endl;
    	}
    }*/
    protected:
    private:
};

#endif // FLOWRATE_H
