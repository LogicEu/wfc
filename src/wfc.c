#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include <imgtool.h>
#include <utopia.h>

#define SPXE_APPLICATION
#include <spxe.h>

#define bug printf("%s, %s, %d\n", __FILE__, __func__, __LINE__)

typedef struct int2 {
    int data[2];
} int2;

typedef struct int4 {
    int data[4];
} int4;

typedef enum heuristic_t {
    Entropy, 
    MRV, 
    Scanline
} heuristic_t;

typedef struct model_t {
    bmp_t bmp;
    uint8_t* sample;
    array_t colors;
    array_t patterns;
    array_t weights;
    array_t stack;
    map_t pindices;
    bool* wave;
    int4* compatible;
    double* weightLog;
    double* distribution;
    double* sumsofweights;
    double* sumsofweightsLog;
    double* entropies;
    double sumweights;
    double sumweightlogs;
    double start_entropy;
    int* sumofones;
    int* observed;
    int observedsofar;
    array_t* propagator;
    unsigned int width;
    unsigned int height;
    unsigned int N;
    int symmetry;
    bool periodic;
    bool ground;
    bool init;
    heuristic_t heuristic;
} model_t;

static Px* fb;
static bmp_t outbmp;

static int dirx[4] = { -1, 0, 1, 0 };
static int diry[4] = { 0, 1, 0, -1 };
static int opposite[4] = { 2, 3, 0, 1 };

/* random floating point operations */

#define RNDNUM 10000000
#define FRNDNUM 10000000.0

static double random_double(void)
{
    int i = rand() % RNDNUM;
    return (double)i / FRNDNUM;
}

static int random_distribution(double* weights, double r, const int len)
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

static uint8_t* rotate(const uint8_t* p, const int N)
{
    uint8_t* result = malloc(N * N);
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            result[y * N + x] = p[N - 1 - y + x * N];
        }
    }
    return result;
}

static uint8_t* reflect(const uint8_t* p, const int N)
{
    uint8_t* result = malloc(N * N);
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            result[y * N + x] = p[N - 1 - x + y * N];
        }
    }
    return result;
}

static uint8_t* sampler(const uint8_t* p,   const int width,    const int height, 
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

static bool agrees(const uint8_t* p1,   const uint8_t* p2, 
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

static long phash(const uint8_t* p, const int len, const long C)
{
    long result = 0, power = 1;
    for (int i = 0; i < len; i++) {
        result += p[len - 1 - i] * power;
        power *= C;
    }
    return result;
}

/* create & destroy model */

static void model_parse(model_t* model)
{
    const uint32_t* px = (uint32_t*)model->bmp.pixels;
    const uint32_t* cols = (uint32_t*)model->colors.data;
    const int len = model->bmp.width * model->bmp.height;
    const int N = model->N;

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

        model->sample[i] = (uint8_t)j;
    }

    const int nsqr = N * N;
    const int xmax = model->periodic ? model->bmp.width : model->bmp.width - N + 1;
    const int ymax = model->periodic ? model->bmp.height : model->bmp.height - N + 1;

    double* w = model->weights.data, one = 1.0;
    for (int y = 0; y < ymax; ++y) {
        for (int x = 0; x < xmax; ++x) {
            uint8_t* ps[8];
            ps[0] = sampler(model->sample, model->bmp.width, model->bmp.height, N, x, y);
            ps[1] = reflect(ps[0], N);
            ps[2] = rotate(ps[0], N);
            ps[3] = reflect(ps[2], N);
            ps[4] = rotate(ps[2], N);
            ps[5] = reflect(ps[4], N);
            ps[6] = rotate(ps[4], N);
            ps[7] = reflect(ps[6], N);

            for (int i = 0; i < 8; i++) {
                if (i < model->symmetry) {
                    uint8_t* p = ps[i];
                    const long h = phash(p, nsqr, model->colors.size);
                    size_t index = map_search(&model->pindices, &h);
                    if (index--) {
                        w[index] += 1.0;
                        free(ps[i]);
                    } else {
                        map_push(&model->pindices, &h, &model->weights.size);
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

    const int T = model->weights.size;
    model->propagator = malloc(4 * T * sizeof(array_t));

    const uint8_t** ps = model->patterns.data;
    for (int d = 0; d < 4; d++) {
        for (int t = 0; t < T; t++) {
            array_t arr = array_create(sizeof(int));
            for (int t2 = 0; t2 < T; t2++) {
                if (agrees(ps[t], ps[t2], dirx[d], diry[d], N)) {
                    array_push(&arr, &t2);
                    //printf("%d, %d, %d\n", d, t, t2);
                }
            }
            model->propagator[d * T + t] = arr;
        }
    }
}

static void 
model_init(model_t* model)
{
    const int T = model->weights.size;
    const int wavelen = model->width * model->height;
    const int len = wavelen * T;
    
    model->wave = malloc(len * sizeof(bool));
    model->compatible = malloc(len * sizeof(int4));

    model->observed = malloc(sizeof(int) * wavelen);

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
    model->init = true;
}

model_t 
model_create(bmp_t bmp, unsigned width, unsigned height, unsigned N, 
            bool periodic, bool ground, heuristic_t heuristic, int symmetry)
{
    model_t model;
    model.width = width;
    model.height = height;
    model.N = N;
    model.bmp = bmp;
    model.heuristic = heuristic;
    model.periodic = periodic;
    model.symmetry = symmetry;
    model.ground = ground;
    model.init = false;
    model.sample = malloc(model.bmp.width * model.bmp.height);
    model.colors = array_create(sizeof(uint32_t));
    model.patterns = array_create(sizeof(uint8_t*));
    model.weights = array_create(sizeof(double));
    model.pindices = map_create(sizeof(long), sizeof(size_t));
    map_overload(&model.pindices, &hash_uint);
    model_parse(&model);
    return model;
}

void 
model_destroy(model_t* model)
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

    free(model->sample);
    free(model->weightLog);
    free(model->entropies);
    free(model->distribution);
    free(model->sumsofweights);
    free(model->sumsofweightsLog);
    free(model->sumofones);
    free(model->observed);
    free(model->propagator);
    free(model->compatible);
    free(model->wave);

    map_free(&model->pindices);
    bmp_free(&model->bmp);
}

/* wave function collapse functions */

static int model_scan(model_t* model)
{
    const int wavelen = model->width * model->height;
    if (model->heuristic == Scanline) {
        for (int i = model->observedsofar; i < wavelen; i++) {
            if (!model->periodic && 
                (i % model->width + model->N > model->width ||
                 i / model->width + model->N > model->height)) {
                continue;
            }
            
            if (model->sumofones[i] > 1) {
                model->observedsofar = i + 1;
                return i;
            }
        }
        return -1;
    }

    double min = 1E+4;
    int argmin = -1;
    for (int i = 0; i < wavelen; i++) {
        if (!model->periodic && 
            (i % model->width + model->N > model->width || 
             i / model->width + model->N > model->height)) {
            continue;
        }
        
        int remainingValues = model->sumofones[i];
        double entropy = model->heuristic == Entropy ? model->entropies[i] : (double)remainingValues;
        
        if (remainingValues > 1 && entropy <= min) {
            double noise = 1E-6 * random_double();
            if (entropy + noise < min) {
                min = entropy + noise;
                argmin = i;
            }
        }
    }
    return argmin;
}

static void model_clear(model_t* model);
bool model_run(model_t* model, const int seed, const int limit);

static int model_waveat(const bool* wave, const int i, const int T)
{
    int count = 0;
    for (int n = 0; n < T; ++n) {
        if (wave[i * T + n]) {
            ++count;
        }
    }
    return count;
}

static bool model_propagate(model_t* model);

static void model_ban(model_t* model, const int i, const int t)
{
    const int T = model->weights.size;
    double* weights = model->weights.data;

    printf("BAN: %d, %d, %d\n", i, t, model->sumofones[i]);
    
    int count = 0, m = -1;
    for (int n = 0; n < T; ++n) {
        if (model->wave[i * T + n]) {
            ++count;
            m = n;
        }
    }
    
    assert(count);
    assert(count == model->sumofones[i]);

    if (model->sumofones[i] == 1 && m == t) {
        //model_clear(model);
        //model_run(model, 0, -1);
        printf("Paradox [%d]\n", t);
        assert(0);
        return;
    }

    model->wave[i * T + t] = false;
    memset(model->compatible + i * T + t, 0, sizeof(int4));
    
    int2 rec = {i, t};
    if (!array_push_if(&model->stack, &rec)) {

        model->sumofones[i] -= 1;
        model->sumsofweights[i] -= weights[t];
        model->sumsofweightsLog[i] -= model->weightLog[t];

        double sum = model->sumsofweights[i];
        assert(sum != 0.0);
        model->entropies[i] = log(sum) - model->sumsofweightsLog[i] / sum;
    }
}

static void model_observe(model_t* model, const int node)
{
    const int T = model->weights.size;
    const int K = node * T;
    
    bool* wave = model->wave;
    double* weights = model->weights.data;
    int count = 0;
    for (int t = 0; t < T; t++) {
        if (wave[K + t]) ++count;
        model->distribution[t] = wave[K + t] ? weights[t] : 0.0;
        //printf("%lf, %lf\n", weights[t], model->distribution[t]);
    }
    printf("Observe: %d, %d\n", count, node);
    //assert(count);
    //assert(count != 1);
    
    int r = random_distribution(model->distribution, random_double(), T);
    int precount = 0, postcount = 0;
    for (int t = 0; t < T; t++) {
        precount += wave[K + t];
        if (wave[K + t] != (t == r)) {
            //printf("1[%d]\n", model_waveat(wave, node, T));
            model_ban(model, node, t);
        } 
        //if (wave[K + t]) printf("%d, %lf, %d\n", t, weights[t], model->sumofones[node]);
        postcount += wave[K + t];
    }
    //printf("\n");
    assert(precount != 0);
    assert(postcount == 1);
}

void model_save(const model_t* model, bmp_t* bmp);
static void fbprnt(Px* fb, const bmp_t* bmp, const int x, const int y);

static bool model_propagate(model_t* model)
{
    const int T = model->weights.size;
    const int N = model->N;
    while (model->stack.size > 0) {
        
        const int2 n = *(int2*)array_pop(&model->stack);
        const int i1 = n.data[0];
        const int t1 = n.data[1];

        printf("PRO: %d, %d, %d\n", i1, t1, model->sumofones[i1]);

        const int x1 = i1 % model->width;
        const int y1 = i1 / model->width;

        for (int d = 0; d < 4; d++) {
            int x2 = x1 + dirx[d];
            int y2 = y1 + diry[d];
            if (!model->periodic && (x2 < 0 || y2 < 0 || 
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
                if (comp->data[d] == 0) {
                    //printf("2[%d], %d\n", model_waveat(model->wave, i2, T), model->sumofones[i2]);
                    model_ban(model, i2, t2);
                }
            }
        }
    }

    return model->sumofones[0] > 0;
}

static void model_clear(model_t* model)
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
        model->observed[i] = -1;
    }

    model->observedsofar = 0;

    if (model->ground) {
        for (int x = 0; x < (int)model->width; x++) {
            for (int t = 0; t < T - 1; t++) {
                model_ban(model, x + (model->height - 1) * model->width, t);
            }
            
            for (int y = 0; y < (int)model->height - 1; y++) {
                model_ban(model, x + y * model->width, T - 1);
            }
        }
        model_propagate(model);
    }
}

bool model_step(model_t* model)
{
    const int T = model->weights.size;
    const int wavelen = model->width* model->height;
    int node = model_scan(model);
    if (node >= 0) {
        model_observe(model, node);
        if (!model_propagate(model)) {
            return false;
        }
    } else {
        for (int i = 0; i < wavelen; i++) {
            int count = 0;
            for (int t = 0; t < T; t++) {
                if (model->wave[i * T + t]) {
                    ++count;
                    model->observed[i] = t;
                    break;
                }
            }
            assert(count == 1);
        }
        printf("Finished?\n");
    }
#if 1
    model_save(model, &outbmp);
    fbprnt(fb, &outbmp, 0, 0);
    if (spxeKeyPressed(ESCAPE) || !spxeRun(fb)) {
        model_destroy(model);
        bmp_free(&outbmp);
        exit(spxeEnd(fb));
    }
#endif
    return true;
}

bool model_run(model_t* model, const int seed, const int limit)
{
    if (!model->init) {
        srand(seed);
        model_init(model);
    }

    model_clear(model);

    for (int l = 0; l < limit || limit < 0; l++) {
        if (!model_step(model)) {
            return false;
        }
    }

    return true;
}

void model_save(const model_t* model, bmp_t* bmp)
{
    const int MX = bmp->width, MY = bmp->height;
    const int T = model->weights.size, N = model->N, wavelen = MX * MY;
    const uint32_t* cols = model->colors.data;
    const uint8_t** ps = model->patterns.data;

#if 0
    if (model->observed[0] >= 0) {
        for (int y = 0; y < MY; y++) {
            int dy = y < MY - N + 1 ? 0 : N - 1;
            for (int x = 0; x < MX; x++) {
                int dx = x < MX - N + 1 ? 0 : N - 1;
                assert(model->observed[x - dx + (y - dy) * MX] >= 0);
                const uint8_t* p = ps[model->observed[x - dx + (y - dy) * MX]];
                const uint32_t* c = cols + p[dx + dy * N];
                memcpy(px_at(bmp, x, y), c, sizeof(uint32_t));
            }
        }
    } else 
#endif
    {
        for (int i = 0; i < wavelen; i++) {
            int x = i % MX, y = i / MX;
#if 0
            if (model->sumofones[i] == 1) {
                int dx = x < MX - N + 1 ? 0 : N - 1;
                int dy = y < MY - N + 1 ? 0 : N - 1;
                const uint8_t** ps = model->patterns.data;
                for (int t = 0; t < T; ++t) {
                    if (model->wave[i * T + t]) {
                        const uint8_t* p = ps[t];
                        const uint32_t* c = cols + p[dx + dy * N];
                        memcpy(px_at(bmp, x, y), c, sizeof(uint32_t));
                        break;
                    }
                }
                continue;
            }
#endif
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
                    if (!model->periodic && (sx + N > MX || sy + N > MY || sx < 0 || sy < 0)) {
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

            assert(contributors);
            uint8_t px[4] = {r / contributors, g / contributors, b / contributors, 255};
            memcpy(bmp->pixels + i * 4, px, 4);
        }
    }
}

static void fbprnt(Px* fb, const bmp_t* bmp, const int x, const int y)
{
    const int xpos = spxe.scrres.width / 2 - bmp->width / 2 + x;
    const int ypos = spxe.scrres.height / 2 - bmp->height / 2 + y;
    for (unsigned int y = 0; y < bmp->height; ++y) {
        memcpy(
            fb + (bmp->height - y - 1 + ypos) * spxe.scrres.width + xpos, 
            px_at(bmp, 0, y), 
            bmp->width * sizeof(Px)
        );
    }
}

int main(int argc, char** argv)
{
    const char* path = NULL;
    if (argc > 1) {
        path = argv[1];
    }

    if (!path) {
        printf("Missing input image file.\n");
        return EXIT_FAILURE;
    }

    const unsigned int seed = (unsigned int)time(NULL);

    bmp_t bmp = bmp_load(path);
    if (!bmp.pixels) {
        printf("Could not load image file '%s'.\n", path);
        return EXIT_FAILURE;
    }

    outbmp = bmp_new(bmp.width, bmp.height, bmp.channels);
    fb = spxeStart("Wave Function Collapse", 800, 600, 100, 75);

    model_t model = model_create(bmp, 32, 32, 3, false, false, Entropy, 1);
    model_init(&model);
    model_clear(&model);
    if (!model_run(&model, seed, -1)) {
        printf("quantum contradiction found!!\n");
        model_destroy(&model);
        return EXIT_FAILURE;
    }
    
    model_save(&model, &outbmp);
    fbprnt(fb, &outbmp, 0, 0);

    while (spxeRun(fb)) {
        
        if (spxeKeyPressed(ESCAPE)) {
            break;
        }
        
        if (spxeKeyPressed(X)) {
            model_step(&model);
        }

        if (spxeKeyPressed(R)) {
            /*if (!model_run(&model, seed, -1)) {
                printf("quantum contradiction found!!\n");
                continue;
            }*/
            model_clear(&model);
            model_save(&model, &outbmp);
            fbprnt(fb, &outbmp, 0, 0);
        }

        if (spxeKeyPressed(SPACE)) {
            static int index = 0;
            uint8_t px[model.N * model.N * 4];
            uint8_t** ps = model.patterns.data;
            for (int y = 0; y < model.N; ++y) {
                for (int x = 0; x < model.N; ++x) {
                    memcpy(&px[(y * model.N + x) * 4], array_index(&model.colors, ps[index][y * model.N + x]), 4);
                }
            }
            bmp_t b = {model.N, model.N, 4, px};
            fbprnt(fb, &b, -32, -32);
            ++index;
            index %= model.patterns.size;
        }
    }
    
    model_destroy(&model);
    bmp_free(&outbmp);
    return spxeEnd(fb);
}
