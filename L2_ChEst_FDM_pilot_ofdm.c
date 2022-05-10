/*
 COMPILE with:
    gcc -o L2_ChEst_FDM_pilot_ofdm.exe L2_ChEst_FDM_pilot_ofdm.c Transmitter.c Receiver.c Channel.c Auxiliary.c -Wall

 Usage example:
    L2_ChEst_FDM_pilot_ofdm.exe 8 5000 16 lin 0
*/

#include "OFDM.h"


void usage(char* progName) {
    printf("\nUsage: %s <modulation> <taps> <equalizer> <symbol blocks> <clipping> <pilots> <interp>\n\n", progName);
    printf("taps: length of the channel delay profile\n");
    printf("symbol blocks: number of symbol blocks to simulate\n");
    printf("pilots: Number of subcarriers dedicated for channel estimation\n");
    printf("        between 3 and < 32\n");
    printf("        power of two (for interp='FFT')\n");
    printf("interp: [ LIN | lin | QUAD | quad | FFT | fft]\n\n");
    printf("Perfect channel information flag:\n");
    printf("   1 - perfect channel knowledge is assumed;\n");
    printf("   0 - channel estimation procedure is employed.\n\n");
}

void output_file_name(char const* filePath, char* fileName, char* howManyTaps, 
                      char* numPilots, char* interpMethod, int PCSI) {

    // char tmp[FILE_NAME_SIZE];
    strcpy(fileName, filePath);
    strcat(fileName, "BER_QPSK_MMSE");

    strcat(fileName, "_L");
    strcat(fileName, howManyTaps);

    if (PCSI) {
        strcat(fileName, "_PCSI.txt");
    }
    else {
        strcat(fileName, "_plt");
        strcat(fileName, numPilots);

        strcat(fileName, "_");
        strcat(fileName, interpMethod);
        strcat(fileName, ".txt");
    }
}



int main(int argc, char** argv) {
    
    /*
     * Initialize the random generator. Must be done once per run of main
     */
    srand((unsigned) time(NULL));
    
    // set if SNR-per-bit values for BER/SER evalation
    double EbN0[] = {0, 5, 10, 15, 20, 25, 30};

    // check required if all required parameters are provided
    if (argc != 6) {
        usage(argv[0]);
        return(EXIT_FAILURE);
    }
    
    // set number of bits per symbol (corresponding to QPSK)
    int bitsPerSymbol = QPSK;

    // set equalizer type (MMSE)
    equalizerType equalizer = MMSE;
    
    // read the number of channel taps from parameter list
    int taps = atoi(argv[1]);
    
    // number of OFDM blocks to simulate
    int howManyBlocks = atoi(argv[2]);

    // set number of pilot carriers
    int num_pilots = atoi(argv[3]);

    // set interpolation method
    char interp[10];
    strcpy(interp, argv[4]);    

    // accociate appropriate interpolation function
    void (*interpolate)(double[], double[], int[], int, int); // function pointer
    if (!strcmp(interp, "LIN") | !strcmp(interp, "lin")) {
        interpolate = interp_linear;
    } else if (!strcmp(interp, "QUAD") | !strcmp(interp, "quad")) {
        interpolate = interp_quadratic;
    } else if (!strcmp(interp, "FFT") | !strcmp(interp, "fft")) {
        interpolate = interp_fft;
    } else {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    // check if constraints on number of pilots are satisfied
    if (num_pilots < 3 || num_pilots > SYMBOLS_PER_BLOCK/2) {
        usage(argv[0]);
        return(EXIT_FAILURE);
    } else if ( !strcmp(interp, "FFT") | !strcmp(interp, "fft") ) {
        // for FFT interpolation, number should be a power of 2
        int num = num_pilots;
        while( num != 1) {
            if(num % 2 != 0) {
                usage(argv[0]);
                return(EXIT_FAILURE);
            }
            num /= 2;
        }
    }

    // number of data subcarriers
    int num_data = (SYMBOLS_PER_BLOCK-num_pilots);

    // set perfect channel state information flag
    int PCSI = atoi(argv[5]);

    // Throughput loss due to GI insertion
    double cpLoss = (double) SYMBOLS_PER_BLOCK / (SYMBOLS_PER_BLOCK + CP_LEN);

    // Throughput loss due to pilot insertion
    double pltLoss;
    if (PCSI) pltLoss = 1.0; // no loss for perfect CSI
    else pltLoss = (double) num_data/SYMBOLS_PER_BLOCK;
 
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
        snr[i] *= pltLoss;          // pilot insertion
        sqrtSNR[i] = sqrt(snr[i]);  // amplitude scaling
    }

    // arrays to keep Tx and Rx bits in each iteration (i.e. one OFDM symbol)
    int howManyBits = num_data*bitsPerSymbol;
    unsigned char *txBits = (unsigned char*) malloc(howManyBits * sizeof(unsigned char));
    unsigned char *rxBits = (unsigned char*) malloc(howManyBits * sizeof(unsigned char));

    // Tx data symbols
    double* dataI = (double*) malloc(num_data * sizeof(double));
    double* dataQ = (double*) malloc(num_data * sizeof(double));

    // Tx symbols (whole OFDM frame)
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


    // data and pilot subcarriers' indexes
    int* dataIndx = (int*) malloc(num_data * sizeof(int));
    int* pltIndx = (int*) malloc(num_pilots * sizeof(int));

    // pilot sequence
    double* pltI = (double*) malloc((num_pilots) * sizeof(double));
    double* pltQ = (double*) malloc((num_pilots) * sizeof(double));

    // channel transfer function estimate
    double* HestI = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));
    double* HestQ = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));

    
    // Bit and symbol error calculation
    double* ber = (double*) malloc(sizeof(EbN0));


    // >> SIMULATION <<

    // terminal output header
    printf("\nEbNo\tBER\n");

    // get data and pilot indexes (for FFT, last carrier is not pilot)
    if (!strcmp(interp, "FFT") | !strcmp(interp, "fft")) {
        get_data_pilot_indexes(dataIndx, pltIndx, SYMBOLS_PER_BLOCK, num_pilots, 0);
    }
    else {
        get_data_pilot_indexes(dataIndx, pltIndx, SYMBOLS_PER_BLOCK, num_pilots, 1);
    }

    // generate TDM pilot sequence (whole OFDM block)
    generate_pilot_qpsk(pltI, pltQ, num_pilots);

    // insert pilot (only needs to be done once)
    insert_symbols_fdm(txSymI, txSymQ, pltI, pltQ, pltIndx, num_pilots);

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

            // generate data bits
            generate_info_bits(howManyBits, txBits);
            
            // modulate bit sequence
            generate_symbols(txBits, bitsPerSymbol, num_data, dataI, dataQ);

            // insert data into data carriers
            insert_symbols_fdm(txSymI, txSymQ, dataI, dataQ, dataIndx, num_data);

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
            // NOTE: Reverts the pilot scaling in path_loss function
            for (int k = CP_LEN; k < SYMBOLS_PER_BLOCK+CP_LEN; k++) {
                rxSymI[k] /= (sqrtSNR[i]*SQRT_OF_2);
                rxSymQ[k] /= (sqrtSNR[i]*SQRT_OF_2);
            }

            // OFDM demodulation (CI is dropped)
            fft(rxSymI + CP_LEN, rxSymQ + CP_LEN, SYMBOLS_PER_BLOCK, rxEstI, rxEstQ);


            // check if perfect CSI flag is on
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
            else {
                // Channel Transfer Function (CTF) estimate
                ch_estimation_fdm(HestI, HestQ, rxEstI, rxEstQ, pltI, pltQ, pltIndx, num_pilots);

                // interpolate
                interpolate(HestI, HestQ, pltIndx, SYMBOLS_PER_BLOCK, num_pilots);
            }

            // Equalization (takes channel transfer function as input)
            fde(equalizer, snr[i], HestI, HestQ, rxEstI, rxEstQ, SYMBOLS_PER_BLOCK);

            
            // Demodulate (detect symbols and output data bits)
            for (int k = 0; k < num_data; k++) {
                decode_asymbol(rxEstI + dataIndx[k], rxEstQ + dataIndx[k],
                               rxBits+(k*bitsPerSymbol), bitsPerSymbol);
            }
            
            // Count bit and symbol errors in the block (only for data symbols)
            for (int m = 0; m < num_data; m++) {
                for (int k = 0; k < bitsPerSymbol; k++) {
                    if (txBits[m*bitsPerSymbol + k] != rxBits[m*bitsPerSymbol + k]) {
                        ber[i]++;
                    }
                }
            }
        }
        
        // Calculate BER
        ber[i] /= (howManyBits * howManyBlocks);
        printf("%.1f\t%f\n", EbN0[i], ber[i]);
    }
    
    // generate output file name
    char fileNameBER[FILE_NAME_SIZE];
    output_file_name("./results_fdm/", fileNameBER, argv[1], argv[3], argv[4], PCSI);
    
    // save BER to file
    save_asignal(EbN0, ber, sizeof(EbN0)/sizeof(double), fileNameBER);
    printf("\nBER saved to %s\n", fileNameBER);
    
    // Release all allocated memory
    free(txBits);
    free(rxBits);
    //
    free(dataI);
    free(dataQ);
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
    free(dataIndx);
    free(pltIndx);
    free(pltI);
    free(pltQ);
    free(HestI);
    free(HestQ);
    //
    return EXIT_SUCCESS;
}

