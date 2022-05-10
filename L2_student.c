#include "OFDM.h"

/* Generates pilot sequence (random QPSK symbols)
 * Parameters
 *  pltI    - pilot sequence I-component samples
 *  pltQ    - pilot sequence Q-component samples
 *  len     - pilot sequence length in symbols
 */

void generate_pilot_qpsk(double pltI[], double pltQ[], int len) {
    // Randomly generate QPSK symbols

	//***Student_code_start*** 
	//    Your code...
	//    ...
	//    ...
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
	//    Your code...
	//    ...
	//    ...
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
	//    Your code...
	//    ...
	//    ...
	//***Student_code_end*****lse dataIndx[j++] = i;

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
	//    Your code...
	//    ...
	//    ...
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
	//    Your code...
	//    ...
	//    ...
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
	//    Your code...
	//    ...
	//    ...
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
	//    Your code...
	//    ...
	//    ...
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
	//    Your code...
	//    ...
	//    ...
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
	//    Your code...
	//    ...
	//    ...
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
	//    Your code...
	//    ...
	//    ...
	//***Student_code_end*****
}





