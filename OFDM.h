#ifndef OFDM_H
#define OFDM_H

#include <stdio.h>
// #include <stdbool.h> // locical data type (boolean)
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>

#define LIGHT_SPEED 3e8 // m/s
#define PI 3.1415926535
#define SQRT_OF_2 1.41421356237
    
#define SCALE_BPSK 1.0    
#define SCALE_QPSK 0.707106781186547
#define SCALE_16QAM 0.316227766016838
#define SCALE_64QAM 0.154303349962092 
#define SCALE_256QAM 0.076696498884737
    
#define BPSK 1
#define QPSK 2
#define QAM16 4
#define QAM64 6
#define QAM256 8

#define CLIPPING_OFF -1

#define fDT 1e-4        // 30 m/s, 1 GHz, 1 Mb/s
#define SYMBOLS_PER_BLOCK 64
#define TAPS 8          // selected number of taps
#define MAX_TAPS 16     // maximum allowed number of taps

#define FILE_NAME_SIZE 50 //50
#define FILE_BER "./results/BER" // basis filename for BER
#define FILE_SER "./results/SER" // for SER

// Number of incoming waves (multipaths) for each channel tap
#define N 64

// Cyclic prefix length
#define CP_LEN 16

#define FADING_FILE "MpathFading.txt"

// enumerated list for equalizer type
typedef enum {ZF, MRC, MMSE} equalizerType;


// Simulation termination criteria
#define MAX_ERRS 200   // max. errors
#define MAX_BITS 1.0e8 // max. bits


int set_equalizer_type(int bitsPerSymbol, char* userInput, equalizerType* equalizer);
double rand_ra();
void box_muller(double* real, double* imag);
void fft(double* inI, double* inQ, int len, double* outI, double* outQ);
void ifft(double* inI, double* inQ, int len, double* outI, double* outQ);

void generate_info_bits(int len, unsigned char* bits);
void print_info_bits(unsigned char* bits, int len);
void insert_cp(int signalLen, int giLen, double* txI, double* txQ);
void clipping(double clipFactor, double* txI, double* txQ, int len);
void generate_symbol_bpsk(unsigned char* bits, double* symbol_I, double* symbol_Q);
void generate_symbol_qpsk(unsigned char* bits, double* symbol_I, double* symbol_Q);
void generate_symbol_16qam(unsigned char* bits, double* symbol_I, double* symbol_Q);
void generate_symbol_64qam(unsigned char* bits, double* symbol_I, double* symbol_Q);
void generate_symbol_256qam(unsigned char* bits, double* symbol_I, double* symbol_Q);
int generate_asymbol(unsigned char* bits, int len, double* symbol_I, double* symbol_Q);
int generate_symbols(unsigned char* txBits, int bitsPerSymbol, int howManySymbols, 
        double* txSymI, double* txSymQ);

// CHANNEL
void path_loss(double* symI, double* symQ, int howManySymbols, double sqrtSNR);
void add_noise(double* sI, double* sQ, int len);
void rayleigh_fading_init(int waves, double fadingFactor, 
        double* alpha, double* phi);
void rayleigh_fading(double* alpha, double* phi, int waves, int symbolIndex,
        double* hi, double* hq);
void multipath_fading_init(int waves, int howManyTaps, double fadingFactor,
        double* alpha, double* phi);
void multipath_fading(double* alpha, double* phi, int waves, int howManyTaps,
        int realization, double* hi, double* hq);
void conv_with_gi(double* hI, double* hQ, int lenH, 
        double* inI, double* inQ, int lenIN, int lenGI,
        double* outI, double* outQ);

// CHANNEL ESTIMATION
// Tx side
void generate_pilot_qpsk(double pltI[], double pltQ[], int len);
void insert_pilot_tdm(double frmI[], double frmQ[], double pltI[], double pltQ[], int len);
void get_data_pilot_indexes(int dataIndx[], int pltIndx[], int numCarriers, int numPilots, int lastPilot);
void insert_symbols_fdm(double frmI[], double frmQ[], double symI[], double symQ[], int symIndx[], int symLen);

// Rx side
void ch_estimation(double Hi[], double Hq[], double rxSigI[], double rxSigQ[], double pltI[], double pltQ[], int len);
void filter_noise(double Hi[], double Hq[], int len);
void ch_estimation_fdm(double Hi[], double Hq[], double rxSigI[], double rxSigQ[], double pltI[],
                       double pltQ[], int pltIndx[], int pltLen);
void interp_linear(double Hi[], double Hq[], int pltIndx[], int numCarriers, int numPilots);
void interp_quadratic(double Hi[], double Hq[], int pltIndx[], int numCarriers, int numPilots);
void interp_fft(double Hi[], double Hq[], int pltIndx[], int numCarriers, int numPilots);

// EQUALIZATION (at the Rx)
int fde(equalizerType equalizer, double snr, double* Hi, double* Hq,
        double* RXSigI, double* RXSigQ, int len);

// DEMODULATION/DETECTION
void decode_bpsk(double* recSymI, double* recSymQ, unsigned char* recBit);
void decode_qpsk(double* recSymI, double* recSymQ, unsigned char* recBit);
void decode_16qam(double* recSymI, double* recSymQ, unsigned char* recBit);
void decode_64qam(double* recSymI, double* recSymQ, unsigned char* recBit);
void decode_256qam(double* recSymI, double* recSymQ, unsigned char* recBit);
int decode_asymbol(double* recSymI, double* recSymQ,
        unsigned char* recBits, int bitsPerSymbol); 
void decode_symbols(double* rxSymI, double* rxSymQ, int howManySymbols,
        int bitsPerSymbol, unsigned char* rxBits);

// UTILITY FUNCTIONS
int save_asignal(double* t, double* signal, int len, char* fileName);
// int save_complex_signal(int* k, double* Si, double* Sq, int len, char* fileName);

#endif /* OFDM_H */

