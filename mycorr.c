#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "sacio.h"

int main(int argc, char *argv[]) {
    int i, npts, tap_npts;
    float *data1, *data2, *data_cor, scale, mean = 0.;
    fftw_complex *in1, *in2, *out1, *out2, *cor_in, *cor_out;
    fftw_plan p1, p2, p3;
    SACHEAD hd1, hd2;

    if ( argc != 4 ) {
        fprintf(stderr,"Usage: my_corr <sacfile1> <sacfile2> <cor-file>\n");
        fprintf(stderr,"       return cross-correlation file of sacfile1 and sacfile2\n");
        fprintf(stderr,"       <sacfile1> The first inputing SAC format file;\n");
        fprintf(stderr,"       <sacfile2> The second inputing SAC foramt file;\n");
        fprintf(stderr,"       <cor-file> The cross-correlation SAC foramt file.\n");
        exit(1);
    }

    data1 = read_sac(argv[1],&hd1);
    data2 = read_sac(argv[2],&hd2);
    npts = 2*hd1.npts-1;
    tap_npts = (int)(npts*5./220.);
    scale = npts;

    in1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);
    in2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);
    out1 =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);
    out2 =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);
    cor_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);
    cor_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);
    data_cor = (float*) malloc(sizeof(float) * npts);

    for (i = 0; i < npts; i ++) {
        if ( i < hd1.npts ) {
            in1[i][0] = data1[hd1.npts-i-1];
            in2[i][0] = data2[hd2.npts-i-1];
        }
        else {
            in1[i][0] = data1[npts-i+1];
            in2[i][0] = data2[npts-i+1];
        }
        in1[i][1] = 0.;
        in2[i][1] = 0.;
    }

    p1 = fftw_plan_dft_1d(npts, in1, out1, FFTW_FORWARD, FFTW_ESTIMATE);
    p2 = fftw_plan_dft_1d(npts, in2, out2, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(p1);
    fftw_execute(p2);

    for ( i = 0; i < npts; i ++ ) {
        cor_in[i][0] = out1[i][0]*out2[i][0] + out1[i][1]*out2[i][1];
        cor_in[i][1] = out1[i][1]*out2[i][0] - out1[i][0]*out2[i][1];
    }

    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_free(in1); fftw_free(in2); fftw_free(out1); fftw_free(out2);

    p3 = fftw_plan_dft_1d(npts, cor_in, cor_out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p3);

    for ( i = 0; i < npts; i ++ ) data_cor[i] = cor_out[i][0]/scale;
    hd1.npts = npts; hd1.b = 0. - hd1.e;

    for ( i = 0; i < tap_npts; i ++ ) mean += data_cor[i]/tap_npts;

    for ( i = 0; i < tap_npts; i ++ ) {
        data_cor[i] = sin(i*hd1.delta/100.)*mean;
        data_cor[npts-tap_npts+i] = cos((npts-tap_npts+i)*hd1.delta/100.)*mean;
    }

    write_sac(argv[3],hd1,data_cor);
    fftw_destroy_plan(p3);
    fftw_free(cor_in); fftw_free(cor_out); free(data1); free(data2); free(data_cor);


    return 0;



}
