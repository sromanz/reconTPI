#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <omp.h> 

int calcADCSampling(double **pADCSampleMatrix, double *pRampSamplingTrajArray);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{


 // Arrays for transporting the trajectory
  double ** m_dk;
  double *  m_dDiffTraj;
  long lI;
  long m_iNColMeas = 256;
  long m_iNEco=1;

    // ADC samples
    m_dk          = new double * [m_iNColMeas + 1];
    for (lI=0; lI < m_iNColMeas + 1; lI++)  m_dk[lI]    = new double [m_iNEco];    

    // Ramp Sampling Trajectory for phase correction
    m_dDiffTraj   = new double [m_iNColMeas + 1];
  
  if (nrhs != 1 ) mexErrMsgTxt("usage: [traces dirs] = mexLoadSiemensTraces(filename)\n");
    
    mxArray *tmp;
    double *pre;
    int N;
    int i;


    pre = mxGetPr(prhs[0]);
    N = mxGetNumberOfElements(prhs[0]);
    
  

    // m_dGradDelayTime = prhs[0];//  -0.76;
    // double m_dRORampTime = prhs[1]60.0;
    // long m_iNEco = 1;
    // long m_iNColMeas = 256;
    // bool m_bIsUTE = true;
    // bool m_bMonopolarReadOut = true;
    // int m_viDwellTime[3] = {10000,10000,10000};

  /* Construct outdata */
  mwSize outdims[1];
  outdims[0] = m_iNColMeas;

  plhs[0] = mxCreateNumericArray(1, outdims, mxDOUBLE_CLASS, mxREAL);
  double* outputBuffer = (double*)mxGetData(plhs[0]);

  plhs[1] = mxCreateNumericArray(1, outdims, mxDOUBLE_CLASS, mxREAL);
  double* outputBuffer2 = (double*)mxGetData(plhs[1]);
      
      calcADCSampling(m_dk,m_dDiffTraj);


      //float *outdata = new float[m_iNColMeas];

      for (int k=0; k<m_iNColMeas;k++){

        outputBuffer[k] = m_dk[k][0];
        outputBuffer2[k] =m_dDiffTraj[k];

    }


      mexPrintf("program finished\n");
}

int calcADCSampling(double **pADCSampleMatrix, double *pRampSamplingTrajArray)
{
   
    const bool bDebug = false;
    
    double  dTimeSample;    
    double  dTime;
    double  dkFT;
    double  dIndex;
    long    lIndex;
    long    lI = 0;
    long    lEco;
    
    //hard coded functor internal variables
    double m_dGradDelayTime = -0.76;
    double m_dRORampTime = 60.0;
    long m_iNEco = 1;
    long m_iNColMeas = 256;
    bool m_bIsUTE = true;
    bool m_bMonopolarReadOut = true;
    int m_viDwellTime[3] = {7000,7000,7000};

    double  dGradDelayTime = - m_dGradDelayTime; // in microseconds
    double  dRampTime      = m_dRORampTime;      // ramp time in microseconds
    double  dGradientInt   = 10;                 // gradient interval in microseconds (the gradient has 10 microseconds intervals)
    double  dSamplingInt   = 0;                  // sampling interval in microseconds
    double  * adSamplingInt;                     // array of sampling interval in microseconds
    adSamplingInt = new double [m_iNEco];
    for (lEco = 0; lEco < m_iNEco; lEco++)  
    {
        adSamplingInt[lEco] = 1.0E-3 * (double) m_viDwellTime[lEco]; 
        dSamplingInt = (dSamplingInt < adSamplingInt[lEco]) ? adSamplingInt[lEco] : dSamplingInt;
    }   

    long lGradSamples = ceil( double(m_iNColMeas + 10) * dSamplingInt / dGradientInt ); 
    long lRampRasterPoints = long (floor (dRampTime / dGradientInt));
    if (!( lGradSamples > lRampRasterPoints ))  lRampRasterPoints = lGradSamples - 1; 
    //std::cout<<"ramprasterpoints: "<<lRampRasterPoints<<"\n";
           
    // Gradient Samples
    double ** dr;
    dr            = new double * [lGradSamples];
    for (lI=0; lI < lGradSamples; lI++)    dr[lI]      = new double [m_iNEco];    
    
    // Linear Samples for Phase Correction
    double * dr_lin;
    dr_lin        = new double [lGradSamples];    
    
    for (lEco = 0; lEco < m_iNEco; lEco++)
    {
        dSamplingInt = adSamplingInt[lEco];


        if (lEco == 0 && m_bIsUTE)// Ramp Sampling
        {   
            // std::cout<<"ramp sampling:\n";

            for (lI = 0; lI < lRampRasterPoints + 1; lI++)
            {
                
                dTime = double (lI) * dGradientInt;// time evolution of the gradient                
                dr[lI][lEco] = (dTime * dTime) / (2.0 * dRampTime * dSamplingInt);// ramp
                // std::cout<< dr[lI][lEco]<<" ";

            }
                // std::cout<<"\n";
            
            for (lI = lRampRasterPoints + 1; lI < lGradSamples; lI++)
            {
                dTime = double (lI) * dGradientInt;// time evolution of the gradient
                dr[lI][lEco] = ( dTime / dSamplingInt ) - dr[lRampRasterPoints][lEco];// constant gradient
                
                // std::cout<< dr[lI][lEco]<<" ";
            }
            for (lI = 0; lI < lGradSamples; lI++)
            {
                dTime = double (lI) * dGradientInt;// time evolution of the gradient
                dr_lin[lI] = dTime / dSamplingInt; // constant gradient - asymmetric echo
            }

        }
        else// No Ramp Sampling
        {
            for (lI = 0; lI < lGradSamples; lI++)
            {
                dTime = double (lI) * dGradientInt;// time evolution of the gradient
                dr[lI][lEco] = dTime / dSamplingInt - (m_iNColMeas / 2.0);// constant gradient - symmetric echo
            }

            dGradDelayTime = 0.0;// For symmetric echoes the correction takes place in the sequence
        }                


        for (lI = 0; lI < m_iNColMeas + 1; lI++)
        {
            dTimeSample = lI * dSamplingInt;// sampling time evolution of the ADC
            dIndex = ( dTimeSample - dGradDelayTime ) / dGradientInt;
            if (dIndex < 0 ) dIndex = 0.0;
            if (dIndex > lGradSamples - 2 ) dIndex = lGradSamples - 2;
            lIndex = (long) floor (dIndex);// interpolate radius
            pADCSampleMatrix[lI][lEco] = dr[lIndex][lEco] + (dr[lIndex+1][lEco]-dr[lIndex][lEco]) * (dIndex - (double) lIndex);
            if (m_bIsUTE && (lEco == 0))// Ramp Sampling
            {
                //dkFT = dr_lin[lIndex] + (dr_lin[lIndex+1]-dr_lin[lIndex]) * (dIndex - (double) lIndex);
                dkFT = dr_lin[lIndex];
                pRampSamplingTrajArray[lI] = dkFT - pADCSampleMatrix[lI][lEco];
            }
        }
        
        // Polarity of the second echo
        if ((lEco > 0) && !m_bMonopolarReadOut)
        {

            for (lI = 0; lI < m_iNColMeas + 1; lI++)
            {
                pADCSampleMatrix[lI][lEco] = pADCSampleMatrix[lI][lEco] * pow(-1.0,lEco);
            }
        }


           
    }// endfor lEco

    if( adSamplingInt != NULL )
    {
        delete adSamplingInt;
    }
    
    if ( dr != NULL )
    {  
        for (lI=0 ; lI < lGradSamples; lI++)    delete [] dr[lI] ;
        delete [] dr;
        dr = NULL;
    }

    if ( dr_lin != NULL )
    {
        delete dr_lin;
    }    

    return 0;   
}
