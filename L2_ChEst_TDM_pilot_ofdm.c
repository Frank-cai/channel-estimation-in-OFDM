/*
 COMPILE with:
    gcc -o L2_ChEst_TDM_pilot_ofdm.exe L2_ChEst_TDM_pilot_ofdm.c Transmitter.c Receiver.c Channel.c Auxiliary.c -Wall

 Usage example:
    L2_ChEst_TDM_pilot_ofdm.exe 8 5000 2 0
*/

#include "OFDM.h"



void usage(char* progName) {
    printf("\nUsage: %s <taps> <symbol blocks> <pilot rate>\n\n", progName);
    printf("taps: length of the channel delay profile\n");
    printf("symbol blocks: number of symbol blocks to simulate\n");
    printf("Pilot rate (>= 2): Distance between TDM pilot frames\n");
    printf("Perfect channel information flag:\n");
    printf("   1 - perfect channel knowledge is assumed;\n");
    printf("   0 - channel estimation procedure is employed.\n\n");
}

void output_file_name(char const* filePath, char* fileName, 
        char* howManyTaps, char* pilotRate, int PCSI) {

    strcpy(fileName, filePath);
    strcat(fileName, "BER_QPSK_MMSE");

    strcat(fileName, "_L");
    strcat(fileName, howManyTaps);


    if (PCSI) {
        strcat(fileName, "_PCSI.txt");
    }
    else {
        strcat(fileName, "_PltRate");
        strcat(fileName, pilotRate);
        strcat(fileName, ".txt");
    }
}



int main(int argc, char** argv) {
    
    // Initialize the random generator. Must be done once per run of main
    srand((unsigned) time(NULL));
    
    // set if SNR-per-bit values for BER/SER evalation
    double EbN0[] = {0, 5, 10, 15, 20, 25, 30};

    // check required if all required parameters are provided
    if (argc != 5) {
        usage(argv[0]);
        return(EXIT_FAILURE);
    }
    
    // set number of bits per symbol (corresponding to QPSK)
    int bitsPerSymbol = QPSK;

    // set equalizer type (MMSE)
    equalizerType equalizer = MMSE;
    
    // read the number of channel taps from parameter list
    int taps = atoi(argv[1]);
    
    // set number of blocks to simulate
    int howManyBlocks = atoi(argv[2]);
    int howManyBits = SYMBOLS_PER_BLOCK * bitsPerSymbol; // bits per block (OFDM symbol)

    // set pilot insertion rate
    int pltInsRate = atoi(argv[3]);

    // check if appropriate value is specified
    if (pltInsRate < 2) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    // set perfect channel state information flag
    int PCSI = atoi(argv[4]);


    // Throughput loss due to GI insertion
    double cpLoss = (double) SYMBOLS_PER_BLOCK / (SYMBOLS_PER_BLOCK + CP_LEN);

    // Throughput loss due to pilot insertion
    double pltLoss;
    if (PCSI) pltLoss = 1.0; // no loss for perfect CSI
    else pltLoss = (double) (pltInsRate-1)/pltInsRate;
 
    // Signal magnitude adjustment to match specified EbNo
    // NOTE: Noise is assumed to be of unit power per I/Q component
    double* snr = (double*) malloc(sizeof(EbN0));
    double* sqrtSNR = (double*) malloc(sizeof(EbN0));
    for (int i = 0; i < sizeof(EbN0) / sizeof(double); i++) {
        snr[i] = pow( 10, EbN0[i] / 10); // SNR per bit
        
        /* 
         * To normalize energy per bit. Models "modulation loss".
         * Makes comparison of different modulations fair.
         */ 
        snr[i] *= bitsPerSymbol;    // SNR per symbol
        snr[i] *= cpLoss;           // GI (cyclic prefix) insertion
        snr[i] *= pltLoss;       // pilot insertion
        sqrtSNR[i] = sqrt(snr[i]);  // amplitude scaling
    }

    // arrays to keep Tx and Rx bits in each iteration (i.e. one OFDM symbol)
    unsigned char txBits[howManyBits];
    unsigned char rxBits[howManyBits];

    // Tx symbols
    double* txSymI = (double*) malloc( SYMBOLS_PER_BLOCK * sizeof(double));
    double* txSymQ = (double*) malloc( SYMBOLS_PER_BLOCK * sizeof(double));
    
    // Tx OFDM symbol (modulated)
    double* txModI = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    double* txModQ = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    
    // Rx OFDM symbol (after channel)
    double* rxSymI = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    double* rxSymQ = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    
    // Rx symbols
    double* rxEstI = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));
    double* rxEstQ = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));
    
    // multipath channel
    double* hi = (double*) malloc(taps * sizeof(double));
    double* hq = (double*) malloc(taps * sizeof(double));
    
    double* alpha = (double*) malloc(N * taps * sizeof(double));
    double* phi = (double*) malloc(N * taps * sizeof(double));

    // Pilot sequence
    double* pltI = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));
    double* pltQ = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));

    // Channel transfer function estimate
    double* HestI = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));
    double* HestQ = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));

    
    // Bit and symbol error calculation
    double* ber = (double*) malloc(sizeof(EbN0));


    // >> SIMULATION <<

    // terminal output header
    printf("\nEbNo\tBER\n");

    // generate TDM pilot sequence (whole OFDM block)
    generate_pilot_qpsk(pltI, pltQ, SYMBOLS_PER_BLOCK);

    // pilot frame indicator
    unsigned int pltFlag = 0;

    // Repeat simulation for each given SNR point
    for (int i = 0; i < sizeof(EbN0) / sizeof(double); i++) {
        
        // initialize bit/block error counters
        ber[i] = 0;

        // Initialize/reset multipath fading simulator
        multipath_fading_init(N, taps, fDT, alpha, phi);

        // OFDM simulation loop (block by block)
        for (int j = 0; j < howManyBlocks; j++) {

            ///////////////////
            //  TRANSMITTER  //
            ///////////////////

            // set pilot flag according to pilot insertion rate
            if (j % pltInsRate == 0) { pltFlag = 1; }
            else { pltFlag = 0; }

            // if perfect CSI is available, no need for pilot frames
            if (PCSI) pltFlag = 0;

            
            if (pltFlag) {
                // insert pilot
                insert_pilot_tdm(txSymI, txSymQ, pltI, pltQ, SYMBOLS_PER_BLOCK);
            }
            else {
                // generate data bits
                generate_info_bits(howManyBits, txBits);
                
                // modulate bit sequence
                generate_symbols(txBits, bitsPerSymbol, SYMBOLS_PER_BLOCK, 
                        txSymI, txSymQ);
            }

            // OFDM modulation
            ifft(txSymI, txSymQ, SYMBOLS_PER_BLOCK, 
                    txModI + CP_LEN, txModQ + CP_LEN);
            
            // insert cyclic prefix
            insert_cp(SYMBOLS_PER_BLOCK, CP_LEN, txModI, txModQ);
            

            // scale Tx symbols to match SNR                    
            path_loss(txModI, txModQ, SYMBOLS_PER_BLOCK + CP_LEN, sqrtSNR[i]);
            
            ///////////////
            //  CHANNEL  //
            ///////////////

            // MULTIPATH FADING
            // update channel impulse response at the beginning of each block
            multipath_fading(alpha, phi, N, taps, j*(SYMBOLS_PER_BLOCK + CP_LEN), hi, hq);
            
            // Convolution with channel impulse response
            conv_with_gi(hi, hq, taps, txModI, txModQ, SYMBOLS_PER_BLOCK + CP_LEN, 
                    CP_LEN, rxSymI, rxSymQ);
            
            // AWGN
            add_noise(rxSymI, rxSymQ, SYMBOLS_PER_BLOCK + CP_LEN);
            
            ////////////////
            //  RECEIVER  //
            ////////////////

            // Automatic Gain Control (AGC)
            // SNR estimation is typically done based on preamble (perfect estimate is assumed here)
            // NOTE: Reverts the pilot scaling introduced in path_loss function
            for (int k = CP_LEN; k < SYMBOLS_PER_BLOCK+CP_LEN; k++) {
                rxSymI[k] /= (sqrtSNR[i]*SQRT_OF_2);
                rxSymQ[k] /= (sqrtSNR[i]*SQRT_OF_2);
            }

            // OFDM demodulation (CI is dropped)
            fft(rxSymI + CP_LEN, rxSymQ + CP_LEN, SYMBOLS_PER_BLOCK, rxEstI, rxEstQ);
            

            // estimate channel or process data frame
            if (pltFlag) {
                // Channel transfer function estimate
                ch_estimation(HestI, HestQ, rxEstI, rxEstQ, pltI, pltQ, SYMBOLS_PER_BLOCK);

                // filter noise from CTF
                filter_noise(HestI, HestQ, SYMBOLS_PER_BLOCK);
            }
            else {
                if (PCSI) {
                    // Perfect channel knowledge
                    for (int k = 0; k < SYMBOLS_PER_BLOCK; k++) {
                        HestI[k] = 0;
                        HestQ[k] = 0;
                        for (int m = 0; m < taps; m++) {
                            HestI[k] += hi[m] * cos(-2 * PI * k * m / SYMBOLS_PER_BLOCK) 
                                    - hq[m] * sin(-2 * PI * k * m / SYMBOLS_PER_BLOCK);
                            HestQ[k] += hq[m] * cos(-2 * PI * k * m / SYMBOLS_PER_BLOCK) 
                                    + hi[m] * sin(-2 * PI * k * m / SYMBOLS_PER_BLOCK);
                        }
                    }
                }

                // Equalization (takes channel transfer function as input)
                fde(equalizer, snr[i], HestI, HestQ, rxEstI, rxEstQ, SYMBOLS_PER_BLOCK);
                
                // Demodulate (detect symbols and output data bits)
                decode_symbols(rxEstI, rxEstQ, SYMBOLS_PER_BLOCK, bitsPerSymbol, rxBits);
                
                // Count bit and symbol errors in the block (only for data symbols)
                for (int m = 0; m < SYMBOLS_PER_BLOCK; m++) {
                    for (int k = 0; k < bitsPerSymbol; k++) {
                        if (txBits[m * bitsPerSymbol + k] != rxBits[m * bitsPerSymbol + k]) {
                            ber[i]++;
                        }
                    }
                }
            }
        }
        
        // Calculate BER
        ber[i] /= (howManyBits * howManyBlocks);
        printf("%f\t%f\n", EbN0[i], ber[i]);
    }
    
    // generate output file name
    char fileNameBER[FILE_NAME_SIZE];
    output_file_name("./results_tdm/", fileNameBER, argv[1], argv[3], PCSI);
    
    // save BER to file
    save_asignal(EbN0, ber, sizeof(EbN0)/sizeof(double), fileNameBER);
    printf("\nBER saved to %s\n", fileNameBER);
    
    // Release all allocated memory
    free(txSymI);
    free(txSymQ);
    free(txModI);
    free(txModQ);
    free(rxSymI);
    free(rxSymQ);
    free(rxEstI);
    free(rxEstQ);
    free(snr);
    free(sqrtSNR);
    free(alpha);
    free(phi);
    free(hi);
    free(hq);
    free(ber);
    //
    free(pltI);
    free(pltQ);
    free(HestI);
    free(HestQ);
    //
    return EXIT_SUCCESS;
}

