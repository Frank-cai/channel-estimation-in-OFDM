#include "OFDM.h"

/* Generates pilot sequence (random QPSK symbols)
 * Parameters
 *  pltI    - pilot sequence I-component samples
 *  pltQ    - pilot sequence Q-component samples
 *  len     - pilot sequence length in symbols
 */

void generate_pilot_qpsk(double pltI[], double pltQ[], int len) {
    // Randomly generate QPSK symbols
	double tmp1, tmp2;
	//***Student_code_start*** 
	for(int i = 0; i < len; i++){
		tmp1 = rand_ra();
		tmp2 = rand_ra();
		pltI[i] = (tmp1 > 0.5) ? SCALE_QPSK : -SCALE_QPSK;
		pltQ[i] = (tmp2 > 0.5) ? SCALE_QPSK : -SCALE_QPSK;
	}
	//***Student_code_end*****
}


/* Inserts TDM pilot sequence into the frame (OFDM symbol)
 * Parameters
 *  frmI    - frame sequence (OFDM symbol) I-component
 *  frmQ    - frame sequence (OFDM symbol) Q-component
 *  pltI    - pilot sequence I-component samples
 *  pltQ    - pilot sequence Q-component samples
 *  len     - pilot sequence length in symbols
 */

void insert_pilot_tdm(double frmI[], double frmQ[], double pltI[], double pltQ[], int len) {
	//***Student_code_start*** 
	for(int i = 0; i < len; i++){
		frmI[i] = pltI[i];
		frmQ[i] = pltQ[i];
	}
	//***Student_code_end*****
}


/* Get data and pilot subcarrier indexes.
 *
 * Parameters
 *  dataIndx    - array to be filled with indexes of data subcarries
 *  pltIndx     - array to be filled with indexes of pilot subcarries
 *  numCarriers - total number of carriers (data + pilots)
 *  numPilots   - number of pilot subcarriers
 *  lastPilot   - indicator if the last carruer should be pilot (1) or not (0)
 *  
 * >>
 * NOTE: Edge pilots are better for polynomial interpolation as there is no
 *       extrapolation error at the edges. However, the high-resolution 
 *       interpolation based on FFT fails in that case.
 * >>
 * 
 */


void get_data_pilot_indexes(int dataIndx[], int pltIndx[], int numCarriers, int numPilots, int lastPilot){
	//***Student_code_start*** 
	int m = 0, n = 0;
	if(lastPilot == 1){
		for(int i = 0; i < numCarriers; i++){
			if(i == (floor(m*(numCarriers-1.0)/(numPilots-1)))){pltIndx[m] = i;m++;}
			else{dataIndx[n] = i;n++;}
		}
	}
	else{
		for(int i = 0; i < numCarriers; i++){
			if(i % (numCarriers/numPilots) == 0){pltIndx[m] = i;m++;}
			else{dataIndx[n] = i;n++;}
		}
	}

	//***Student_code_end***

}


/* Insert symbols onto subcarriers with given indexes.
 *
 * Parameters
 *  frmI    - frame sequence (OFDM symbol) I-component
 *  frmQ    - frame sequence (OFDM symbol) Q-component
 *  symIndx - indexes of subcarries to insert symbols
 *  symI    - symbol sequence I-component samples
 *  symQ    - symbol sequence Q-component samples
 *  symLen  - symbol sequence length
 * 
 * >>
 * NOTE: Used for both pilot and data symbols.
 * >>
 */

void insert_symbols_fdm(double frmI[], double frmQ[], double symI[], double symQ[], int symIndx[], int symLen) {
	//***Student_code_start*** 
	for(int i = 0; i < symLen; i++){
		frmI[symIndx[i]] = symI[i];
		frmQ[symIndx[i]] = symQ[i];
	}
	//***Student_code_end*****
}

////////////////////
///Receiver/////////
////////////////////
//////////////////////////
//  CHANNEL ESTIMATION  //
//////////////////////////

/* Estimate channel (frequency) transfer function (CTF), for
 * TDM pilot, i.e. when all carriers in a frame are pilots.
 * 
 * Parameters
 *  Hi      - channel transfer function's I-component
 *  Hq      - channel transfer function's Q-component
 *  pltI    - Rx symbol sequence I-component
 *  pltQ    - Rx symbol sequence Q-component
 *  pltQ    - pilot sequence Q-component samples
 *  pltI    - pilot sequence I-component samples
 *  pltQ    - pilot sequence Q-component samples
 *  len     - pilot sequence length in symbols
 */


void ch_estimation(double Hi[], double Hq[], double rxSigI[], double rxSigQ[], double pltI[], double pltQ[], int len) {
	//***Student_code_start*** 
	for(int i = 0; i < len; i++){
		Hi[i] = (rxSigI[i]*pltI[i]+rxSigQ[i]*pltQ[i])/(pltI[i]*pltI[i]+pltQ[i]*pltQ[i]);
		Hq[i] = (rxSigQ[i]*pltI[i]-rxSigI[i]*pltQ[i])/(pltI[i]*pltI[i]+pltQ[i]*pltQ[i]);
	}
	//***Student_code_end*****
}


/* Filter noise from the estimated CTF.
 *
 * Parameters
 *  Hi  - CTF I-component
 *  Hq  - CTF Q-component
 *  len - common sequence length
 * NOTE: Memory allocation on every run is not efficient.
 *       If implemented as a structure, preallocated array
 *       could be kept internally to avoid repeated allocations.
 */

void filter_noise(double Hi[], double Hq[], int len){
	//***Student_code_start*** 
	double* hi = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));
    double* hq = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));
	for(int i = 0; i < len; i++){
		hi[i] = 0;
		hq[i] = 0;
		for(int j = 0; j < len; j++){
			hi[i] += (Hi[j]*cos(2 * PI * i * j / len)-Hq[j]*sin(2 * PI * i * j / len));
			hq[i] += (Hi[j]*sin(2 * PI * i * j / len)+Hq[j]*cos(2 * PI * i * j / len));
		}
	}

	for(int i = 0; i < len; i++){
		Hi[i] = 0;
		Hq[i] = 0;
		for(int j = 0; j < MAX_TAPS; j++){
			Hi[i] += (hi[j]*cos(- 2 * PI * i * j / len)-hq[j]*sin(- 2 * PI * i * j / len));
			Hq[i] += (hi[j]*sin(- 2 * PI * i * j / len)+hq[j]*cos(- 2 * PI * i * j / len));
		}
	}
	free(hi);
	free(hq);

	//***Student_code_end*****
}


// >> FDM PILOT <<

/* Estimate channel (frequency) transfer function (CTF)
 * at frequency points determined by pilot subcarrier indxes.
 *
 * Parameters
 *  Hi      - CTF I-component
 *  Hq      - CTF Q-component
 *  rxSigI  - Rx symbol sequence I-component
 *  rxSigQ  - Rx symbol sequence Q-component
 *  pltI    - pilot sequence I-component samples
 *  pltQ    - pilot sequence Q-component samples
 *  pltIndx - indexes of pilot subcarriers
 *  pltLen  - pilot sequence length in symbols
 * 
 * >>
 * NOTE: This function can be also used for TDM pilot, but an array
 *       with all indexes (0, 1, ..., numCarriers-1) has tu be provided
 *       in place of pltIndx.
 * >>
 */


void ch_estimation_fdm(double Hi[], double Hq[], double rxSigI[], double rxSigQ[], double pltI[], double pltQ[], int pltIndx[], int pltLen) {
	//***Student_code_start*** 
	for(int i = 0; i < pltLen; i++){
		Hi[pltIndx[i]] = (rxSigI[pltIndx[i]]*pltI[i]+rxSigQ[pltIndx[i]]*pltQ[i])/(pltI[i]*pltI[i]+pltQ[i]*pltQ[i]);
		Hq[pltIndx[i]] = (rxSigQ[pltIndx[i]]*pltI[i]-rxSigI[pltIndx[i]]*pltQ[i])/(pltI[i]*pltI[i]+pltQ[i]*pltQ[i]);
	}
	//***Student_code_end*****
}

/* Linear interpolation of the CTF based on available samples
 * at frequency points determined by pilot subcarrier indxes.
 *
 * Parameters
 *  Hi          - CTF I-component
 *  Hq          - CTF Q-component
 *  pltIndx     - indexes of pilot subcarriers
 *  numCarriers - total number of subcarriers (data + pilots)
 */


void interp_linear(double Hi[], double Hq[], int pltIndx[], int numCarriers, int numPilots) {
	//***Student_code_start*** 
	int n = 0;
	for(int i = 0; i < numCarriers; i++){
		if(pltIndx[n] == i){n++;continue;}
		else{
			Hi[i] = (Hi[pltIndx[n]]-Hi[pltIndx[n-1]])/(pltIndx[n]-pltIndx[n-1])*(i-pltIndx[n-1])+Hi[pltIndx[n-1]];
			Hq[i] = (Hq[pltIndx[n]]-Hq[pltIndx[n-1]])/(pltIndx[n]-pltIndx[n-1])*(i-pltIndx[n-1])+Hq[pltIndx[n-1]];
		}
	}
	//***Student_code_end*****
}


/* Second-order polynomial (quadratic) interpolation of the CTF based on
 * available samples at frequency points determined by pilot subcarrier indxes.
 *
 * Parameters
 *  Hi          - CTF I-component
 *  Hq          - CTF Q-component
 *  pltIndx     - indexes of pilot subcarriers
 *  numCarriers - total number of subcarriers (data + pilots)
 */


void interp_quadratic(double Hi[], double Hq[], int pltIndx[], int numCarriers, int numPilots) {
  	//***Student_code_start*** 
	int n = 0;
	double c0 = 0, c1 = 0, c2 = 0;
	int k;
	k = numCarriers / numPilots;
	for(int i = 0; i < numCarriers; i++){
		if(i == pltIndx[n]){n++;continue;}
		else if (n == numPilots){
			c0 = 1.0*(pltIndx[n]-i-2*k)*(pltIndx[n]-i-k)/(2*k*k);
			c1 = 1.0*(pltIndx[n]-i)*(2*k-pltIndx[n]+i)/(k*k);
			c2 = 1.0*(pltIndx[n]-i)*(pltIndx[n]-i-k)/(2*k*k);
			Hi[i] = c0*Hi[pltIndx[n]]+c1*Hi[pltIndx[n-1]]+c2*Hi[pltIndx[n-2]];
			Hq[i] = c0*Hq[pltIndx[n]]+c1*Hq[pltIndx[n-1]]+c2*Hq[pltIndx[n-2]];
		}
		else{
			c0 = 1.0*(i-pltIndx[n-1]-2*k)*(i-pltIndx[n-1]-k)/(2*k*k);
			c1 = 1.0*(i-pltIndx[n-1])*(2*k-i+pltIndx[n-1])/(k*k);
			c2 = 1.0*(i-pltIndx[n-1])*(i-pltIndx[n-1]-k)/(2*k*k);
			Hi[i] = c0*Hi[pltIndx[n-1]]+c1*Hi[pltIndx[n]]+c2*Hi[pltIndx[n+1]];
			Hq[i] = c0*Hq[pltIndx[n-1]]+c1*Hq[pltIndx[n]]+c2*Hq[pltIndx[n+1]];
		}
	}
	//***Student_code_end*****
}


/* High-resolution interpolation based on FFT.
 *
 * Parameters
 *  Hi           - CTF sequence I-component samples
 *  Hq           - CTF sequence Q-component samples
 *  pltIndx      - indexes of pilot subcarriers
 *  numCarriers  - total number of subcarriers (data + pilots)
 *  numPilots    - total number of pilots
 */
// Interpolates the CTF at missing points by high-resolution DFT method.


void interp_fft(double Hi[], double Hq[], int pltIndx[], int numCarriers, int numPilots) {
  	//***Student_code_start*** 
	int n = 0;
	int k;
	k = numCarriers / numPilots;
	double phiI,phiQ;
	for(int i = 0; i < numCarriers; i++){
		if(i == pltIndx[n]){n++;continue;}
		else{
			Hi[i] = 0;
			Hq[i] = 0;
			for(int q = 0; q < numPilots; q++){
				phiI = 1.0/numPilots*sin(PI*(i-k*q)/k)/sin(PI*(i-k*q)/numCarriers)*cos(-PI*(i-k*q)*(numCarriers-k)/numCarriers/k);
				phiQ = 1.0/numPilots*sin(PI*(i-k*q)/k)/sin(PI*(i-k*q)/numCarriers)*sin(-PI*(i-k*q)*(numCarriers-k)/numCarriers/k);
				Hi[i] += phiI*Hi[pltIndx[q]]-phiQ*Hq[pltIndx[q]];
				Hq[i] += phiQ*Hi[pltIndx[q]]+phiI*Hq[pltIndx[q]];
			}
		}	
	}
	filter_noise(Hi, Hq, numCarriers);
	//***Student_code_end*****
}





