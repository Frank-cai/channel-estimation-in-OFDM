tdm:
	gcc -std=c99 -o L2_ChEst_TDM_pilot_ofdm L2_ChEst_TDM_pilot_ofdm.c L2_student.c -I include lib/libOFDM.so -lm -Wall
fdm:
	gcc -std=c99 -o L2_ChEst_FDM_pilot_ofdm L2_ChEst_FDM_pilot_ofdm.c L2_student.c -I include lib/libOFDM.so -lm -Wall
