#include "wfc.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* random floating point operations */

#define WFC_RNDNUM 10000000
#define WFC_FRNDNUM 10000000.0

static double wfc_random_double(void)
{
    int i = rand() % WFC_RNDNUM;
    return (double)i / WFC_FRNDNUM;
}

static int wfc_random_distribution(double* weights, const double r, const int len)
{
    double sum = 0.0;
    for (int i = 0; i < len; i++) {
        sum += weights[i];
    }
    double threshold = r * sum;

    double partialSum = 0.0;
    for (int i = 0; i < len; i++) {
        partialSum += weights[i];
        if (partialSum >= threshold) {
            return i;
        }
    }
    return 0;
}

/* pattern operations */

static uint8_t* wfc_protate(const uint8_t* p, const int N)
{
    uint8_t* result = malloc(N * N);
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            result[y * N + x] = p[N - 1 - y + x * N];
        }
    }
    return result;
}

static uint8_t* wfc_preflect(const uint8_t* p, const int N)
{
    uint8_t* result = malloc(N * N);
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            result[y * N + x] = p[N - 1 - x + y * N];
        }
    }
    return result;
}

static uint8_t* wfc_psampler(
    const uint8_t* p,   const int width,    const int height, 
    const int N,        const int xpos,     const int ypos)
{
    uint8_t* result = malloc(N * N);
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            result[y * N + x] = p[(x + xpos) % width + ((y + ypos) % height) * width];
        }
    }
    return result;
}

static bool wfc_pagrees(
    const uint8_t* p1,  const uint8_t* p2, 
    const int dx,       const int dy,       const int N)
{
    const int xmin = dx < 0 ? 0 : dx, xmax = dx < 0 ? dx + N : N, ymin = dy < 0 ? 0 : dy, ymax = dy < 0 ? dy + N : N;
    for (int y = ymin; y < ymax; y++) {
        for (int x = xmin; x < xmax; x++) {
            if (p1[x + N * y] != p2[x - dx + N * (y - dy)]) {
                return false;
            }
        }
    }
    return true;
}

static long wfc_phash(const uint8_t* p, const int len, const long C)
{
    long result = 0, power = 1;
    for (int i = 0; i < len; i++) {
        result += p[len - 1 - i] * power;
        power *= C;
    }
    return result;
}

/* static model initialization  */

static int dirx[4] = { -1, 0, 1, 0 };
static int diry[4] = { 0, 1, 0, -1 };
static int opposite[4] = { 2, 3, 0, 1 };

static void wfc_model_init(wfc_model_t* model)
{
    const int T = model->weights.size;
    const int wavelen = model->width * model->height;
    const int len = wavelen * T;
    
    model->wave = malloc(len * sizeof(bool));
    model->compatible = malloc(len * sizeof(int4));

    model->distribution = malloc(sizeof(double) * T);
    model->weightLog = malloc(sizeof(double) * T);
    
    model->sumweights = 0.0;
    model->sumweightlogs = 0.0;

    const double* w = model->weights.data;
    for (int t = 0; t < T; t++) {
        model->weightLog[t] = w[t] * log(w[t]);
        model->sumweights += w[t];
        model->sumweightlogs += model->weightLog[t];
    }

    assert(model->sumweights != 0.0);
    model->start_entropy = log(model->sumweights) - model->sumweightlogs / model->sumweights;

    model->sumofones = malloc(sizeof(int) * wavelen);
    model->sumsofweights = malloc(sizeof(double) * wavelen);
    model->sumsofweightsLog = malloc(sizeof(double) * wavelen);
    model->entropies = malloc(sizeof(double) * wavelen);
    model->stack = array_create(sizeof(int2));
}

static void wfc_model_parse(wfc_model_t* model, const bmp_t* bmp, const int symmetry, bool periodic_input)
{
    const uint32_t* px = (uint32_t*)bmp->pixels;
    const uint32_t* cols = (uint32_t*)model->colors.data;
    const int len = bmp->width * bmp->height;
    const int N = model->N;
    
    uint8_t* sample = malloc(bmp->width * bmp->height);

    for (int i = 0, j; i < len; ++i) {

        for (j = 0; j < (int)model->colors.size; ++j) {
            if (px[i] == cols[j]) {
                break;
            }
        }
        
        if (j == (int)model->colors.size) {
            array_push(&model->colors, px + i);
            cols = model->colors.data;
        }

        sample[i] = (uint8_t)j;
    }

    map_t pindices = map_create(sizeof(long), sizeof(size_t));
    map_overload(&pindices, &hash_uint);

    const int nsqr = N * N;
    const int xmax = periodic_input ? bmp->width : bmp->width - N + 1;
    const int ymax = periodic_input ? bmp->height : bmp->height - N + 1;

    double* w = model->weights.data, one = 1.0;
    for (int y = 0; y < ymax; ++y) {
        for (int x = 0; x < xmax; ++x) {
            
            uint8_t* ps[8];
            ps[0] = wfc_psampler(sample, bmp->width, bmp->height, N, x, y);
            ps[1] = wfc_preflect(ps[0], N);
            ps[2] = wfc_protate(ps[0], N);
            ps[3] = wfc_preflect(ps[2], N);
            ps[4] = wfc_protate(ps[2], N);
            ps[5] = wfc_preflect(ps[4], N);
            ps[6] = wfc_protate(ps[4], N);
            ps[7] = wfc_preflect(ps[6], N);

            for (int i = 0; i < 8; i++) {
                if (i < symmetry) {
                    uint8_t* p = ps[i];
                    const long h = wfc_phash(p, nsqr, model->colors.size);
                    size_t index = map_search(&pindices, &h);
                    if (index--) {
                        w[index] += 1.0;
                        free(ps[i]);
                    } else {
                        map_push(&pindices, &h, &model->weights.size);
                        array_push(&model->weights, &one);
                        array_push(&model->patterns, &p);
                        w = model->weights.data;
                    }
                } else {
                    free(ps[i]);
                }
            }
        }
    }

    map_free(&pindices);
    free(sample);

    const int T = model->weights.size;
    model->propagator = malloc(4 * T * sizeof(array_t));

    const uint8_t** ps = model->patterns.data;
    for (int d = 0; d < 4; d++) {
        for (int t = 0; t < T; t++) {
            array_t arr = array_create(sizeof(int));
            for (int t2 = 0; t2 < T; t2++) {
                if (wfc_pagrees(ps[t], ps[t2], dirx[d], diry[d], N)) {
                    array_push(&arr, &t2);
                }
            }
            model->propagator[d * T + t] = arr;
        }
    }
}

/* low-level wfc operations */

bool wfc_model_ban(wfc_model_t* model, const int i, const int t)
{
    const int T = model->weights.size;
    double* weights = model->weights.data;

    model->wave[i * T + t] = false;
    memset(model->compatible + i * T + t, 0, sizeof(int4));
    
    int2 rec = {i, t};
    array_push(&model->stack, &rec);

    model->sumofones[i] -= 1;
    model->sumsofweights[i] -= weights[t];
    model->sumsofweightsLog[i] -= model->weightLog[t];
    
    if (model->sumofones[i] < 1) {
        return false;
    }

    const double sum = model->sumsofweights[i];
    model->entropies[i] = log(sum) - model->sumsofweightsLog[i] / sum;
    return true;
}

bool wfc_model_propagate(wfc_model_t* model, bool periodic)
{
    const int T = model->weights.size;
    const int N = model->N;
    while (model->stack.size > 0) {
        
        const int2 n = *(int2*)array_pop(&model->stack);
        const int i1 = n.data[0];
        const int t1 = n.data[1];

        const int x1 = i1 % model->width;
        const int y1 = i1 / model->width;

        for (int d = 0; d < 4; d++) {
            int x2 = x1 + dirx[d];
            int y2 = y1 + diry[d];
            if (!periodic && (x2 < 0 || y2 < 0 || 
                x2 + N > (int)model->width || y2 + N > (int)model->height)) {
                continue;
            }

            if (x2 < 0) {
                x2 += model->width;
            } else if (x2 >= (int)model->width) {
                x2 -= model->width;
            }

            if (y2 < 0) {
                y2 += model->height;
            } else if (y2 >= (int)model->height) {
                y2 -= model->height;
            }

            int i2 = x2 + y2 * model->width;
            const array_t pro = model->propagator[d * T + t1];
            const int* p = pro.data;
            for (size_t l = 0; l < pro.size; l++) {
                int t2 = p[l];
                int4* comp = model->compatible + i2 * T + t2;
                --comp->data[d];
                if (!comp->data[d] && !wfc_model_ban(model, i2, t2)) {
                    return false;
                }
            }
        }
    }

    return true;
}

int wfc_model_scan(const wfc_model_t* model, int heuristic, bool periodic)
{
    const int wavelen = model->width * model->height;
    if (heuristic == WFC_SCANLINE) {
        static int observedsofar = 0;
        for (int i = observedsofar; i < wavelen; i++) {
            if (!periodic && 
                (i % model->width + model->N > model->width ||
                 i / model->width + model->N > model->height)) {
                continue;
            }
            
            if (model->sumofones[i] > 1) {
                observedsofar = i + 1;
                return i;
            }
        }
        return -1;
    }

    double min = 1E+4;
    int argmin = -1;
    for (int i = 0; i < wavelen; i++) {
        if (!periodic && 
            (i % model->width + model->N > model->width || 
             i / model->width + model->N > model->height)) {
            continue;
        }
        
        int remainingValues = model->sumofones[i];
        double entropy = heuristic == WFC_ENTROPY ? model->entropies[i] : (double)remainingValues;
        
        if (remainingValues > 1 && entropy <= min) {
            double noise = 1E-6 * wfc_random_double();
            if (entropy + noise < min) {
                min = entropy + noise;
                argmin = i;
            }
        }
    }
    return argmin;
}

bool wfc_model_collapse(wfc_model_t* model, const int node, const int r)
{
    const int T = model->weights.size;
    const int K = node * T;
    for (int t = 0; t < T; t++) {
        if (model->wave[K + t] != (t == r) && !wfc_model_ban(model, node, t)) {
            return false;
        }
    }
    return true;
}

bool wfc_model_observe(wfc_model_t* model, const int node)
{
    const int T = model->weights.size;
    const int K = node * T;
    
    const double* weights = model->weights.data;
    for (int t = 0; t < T; t++) {
        model->distribution[t] = model->wave[K + t] ? weights[t] : 0.0;
    }
    
    const int r = wfc_random_distribution(model->distribution, wfc_random_double(), T);
    return wfc_model_collapse(model, node, r);
}

/* wfc model handling functions */

bool wfc_model_clear(wfc_model_t* model, bool periodic, bool ground)
{
    const int T = model->weights.size;
    const int wavelen = model->width * model->height;
    for (int i = 0; i < wavelen; i++) {
        for (int t = 0; t < T; t++) {
            model->wave[i * T + t] = true;
            for (int d = 0; d < 4; d++) {
                model->compatible[i * T + t].data[d] = (int)model->propagator[opposite[d] * T + t].size;
            }
        }

        model->sumofones[i] = T;
        model->sumsofweights[i] = model->sumweights;
        model->sumsofweightsLog[i] = model->sumweightlogs;
        model->entropies[i] = model->start_entropy;
    }

    if (ground) {
        for (int x = 0; x < (int)model->width; x++) {
            for (int t = 0; t < T - 1; t++) {
                if (!wfc_model_ban(model, x + (model->height - 1) * model->width, t)) {
                    return false;
                }
            }
            
            for (int y = 0; y < (int)model->height - 1; y++) {
                if (!wfc_model_ban(model, x + y * model->width, T - 1)) {
                    return false;
                }
            }
        }
        if (!wfc_model_propagate(model, periodic)) {
            return false;
        }
    }
    return true;
}

wfc_model_t wfc_model_create(   
    const bmp_t* bmp,   unsigned width, unsigned height, 
    unsigned N,         bool periodic,  int symmetry)
{
    wfc_model_t model;
    model.width = width;
    model.height = height;
    model.N = N;
    model.colors = array_create(sizeof(uint32_t));
    model.patterns = array_create(sizeof(uint8_t*));
    model.weights = array_create(sizeof(double));
    wfc_model_parse(&model, bmp, symmetry, periodic);
    wfc_model_init(&model);
    return model;
}

void wfc_model_destroy(wfc_model_t* model)
{
    uint8_t** ps = model->patterns.data;
    for (size_t i = 0; i < model->patterns.size; ++i) {
        free(ps[i]);
    }

    const size_t propcount = model->weights.size * 4;
    for (size_t i = 0; i < propcount; ++i) {
        array_free(model->propagator + i);
    }

    array_free(&model->patterns);
    array_free(&model->colors);
    array_free(&model->weights);
    array_free(&model->stack);

    free(model->weightLog);
    free(model->entropies);
    free(model->distribution);
    free(model->sumsofweights);
    free(model->sumsofweightsLog);
    free(model->sumofones);
    free(model->propagator);
    free(model->compatible);
    free(model->wave);
}

bool wfc_model_step(wfc_model_t* model, int heuristic, bool periodic)
{
    const int node = wfc_model_scan(model, heuristic, periodic);
    return node >= 0 && wfc_model_observe(model, node) && wfc_model_propagate(model, periodic);
}

bool wfc_model_run(wfc_model_t* model, int heuristic, bool periodic, const int steps)
{
    for (int i = 0; i < steps || steps < 0; ++i) {
        if (!wfc_model_step(model, heuristic, periodic)) {
            return false;
        }
    }
    return true;
}

void wfc_model_save(const wfc_model_t* model, uint8_t* pixbuf, bool periodic)
{
    const int MX = model->width, MY = model->height;
    const int T = model->weights.size, N = model->N, wavelen = MX * MY;
    const uint32_t* cols = model->colors.data;
    const uint8_t** ps = model->patterns.data;

    for (int i = 0; i < wavelen; i++) {
        
        int x = i % MX, y = i / MX;
        int contributors = 0, r = 0, g = 0, b = 0;
        
        for (int dy = 0; dy < N; dy++) {
            for (int dx = 0; dx < N; dx++) {
                
                int sx = x - dx;
                if (sx < 0) {
                    sx += MX;
                }

                int sy = y - dy;
                if (sy < 0) {
                    sy += MY;
                }

                int s = sx + sy * MX;
                if (!periodic && (sx + N > MX || sy + N > MY || sx < 0 || sy < 0)) {
                    continue;
                }

                for (int t = 0; t < T; t++) {
                    if (model->wave[s * T + t]) {
                        contributors++;
                        const uint8_t* p = ps[t];
                        const uint8_t* px = (uint8_t*)(cols + p[dx + dy * N]);
                        r += px[0];
                        g += px[1];
                        b += px[2];
                    }
                }
            }
        }

        if (contributors) {
            uint8_t px[4] = {r / contributors, g / contributors, b / contributors, 255};
            memcpy(pixbuf + i * 4, px, 4);
        }
    }
}