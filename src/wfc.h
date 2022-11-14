#ifndef WAVEFUNCOLLAPSE_H
#define WAVEFUNCOLLAPSE_H

#include <stdbool.h>
#include <utopia.h>
#include <imgtool.h>

#define WFC_ENTROPY 0
#define WFC_MRV 1
#define WFC_SCANLINE 2

typedef struct int2 {
    int data[2];
} int2;

typedef struct int4 {
    int data[4];
} int4;

typedef struct wfc_model_t {
    unsigned int N;
    unsigned int width;
    unsigned int height;
    array_t colors;
    array_t patterns;
    array_t weights;
    array_t stack;
    bool* wave;
    int4* compatible;
    int* sumofones;
    double* weightLog;
    double* distribution;
    double* sumsofweights;
    double* sumsofweightsLog;
    double* entropies;
    double sumweights;
    double sumweightlogs;
    double start_entropy;
    array_t* propagator;
} wfc_model_t;

void wfc_model_save(const wfc_model_t* model, uint8_t* pixbuf, bool periodic);
bool wfc_model_step(wfc_model_t* model, int heuristic, bool periodic);
bool wfc_model_run(wfc_model_t* model, int heuristic, bool periodic, const int steps);
bool wfc_model_collapse(wfc_model_t* model, const int node, const int r);
bool wfc_model_observe(wfc_model_t* model, const int node);
int  wfc_model_scan(const wfc_model_t* model, int heuristic, bool periodic);
bool wfc_model_clear(wfc_model_t* model, bool periodic, bool ground);
bool wfc_model_propagate(wfc_model_t* model, bool periodic);
bool wfc_model_ban(wfc_model_t* model, const int i, const int t);

wfc_model_t wfc_model_create(const bmp_t* bmp, unsigned width, unsigned height, unsigned N, bool periodic, int symmetry);
void wfc_model_destroy(wfc_model_t* model);

#endif /* WAVEFUNCOLLAPSE_H */