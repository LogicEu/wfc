#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <limits.h>

#include <utopia.h>
#include <imgtool.h>
#define SPXE_APPLICATION
#include <spxe.h>

#define RNDNUM 1000000
#define bmp_outside(bmp, x, y)\
    ((x) < 0 || (x) >= (int)bmp->width || (y) < 0 || (y) >= (int)bmp->height)
#define swap(a, b, type) do { type tmp = a; a = b; b = tmp; } while (0)
//#define bug printf("%s, %s, %d\n", __FILE__, __func__, __LINE__)

static Px* fb;
static array_t freearr;
static bmp_t* freebmp;

static int sidetable[8][2] = {
    {1, 0}, {1, -1}, {0, -1}, {-1, -1}, {-1, 0}, {-1, 1}, {0, 1}, {1, 1}
};

typedef enum side_t {
    RIGHT, RIGHT_DOWN, DOWN, LEFT_DOWN, LEFT, LEFT_UP, UP, RIGHT_UP
} side_t;

typedef struct rule_t {
    size_t self;
    size_t next;
    size_t side;
} rule_t;

static void wfc_free(void)
{
    array_t* a = freearr.data;
    for (size_t i = 0; i < freearr.size; ++i) {
        array_free(a + i);
    }
    array_free(&freearr);
    bmp_free(freebmp);
    spxeEnd(fb);
}

static array_t array_set_weights(array_t* arr)
{
    array_t cpy = array_copy(arr);
    array_set(arr);
    
    array_t weights = array_reserve(sizeof(float), arr->size);
    const float n = 1.0f / (float)cpy.size;
    
    const size_t arrcount = arr->size, cpycount = cpy.size, bytes = arr->bytes;
    char* d1 = arr->data, *d2;
    for (size_t i = 0; i < arrcount; ++i) {
        d2 = cpy.data;
        float f = 0.0f;
        for (size_t j = 0; j < cpycount; ++j) {
            f += n * !memcmp(d1, d2, bytes);
            d2 += bytes;
        }
        array_push(&weights, &f);
        d1 += bytes;
    }

    array_free(&cpy);
    return weights;
}

static array_t tiles_get(const bmp_t* bmp, array_t* weights)
{
    array_t tiles = array_create(bmp->channels);
    for (unsigned y = 0; y < bmp->height; ++y) {
        for (unsigned x = 0; x < bmp->width; ++x) {
            array_push(&tiles, px_at(bmp, x, y));
        }
    }
    *weights = array_set_weights(&tiles);
    return tiles;
}

static array_t rules_get(const bmp_t* bmp, const array_t* tiles, array_t* weights)
{
    array_t rules = array_create(sizeof(rule_t));
    for (unsigned y = 0; y < bmp->height; ++y) {
        for (unsigned x = 0; x < bmp->width; ++x) {
            size_t self = array_search(tiles, px_at(bmp, x, y)) - 1;
            for (int i = 0; i < 8; ++i) {
                int xpos = (int)x + sidetable[i][0];
                int ypos = (int)y + sidetable[i][1];
                if (bmp_outside(bmp, xpos, ypos)) {
                    continue;
                }

                size_t next = array_search(tiles, px_at(bmp, xpos, ypos)) - 1;

                rule_t rule = {self, next, i};
                array_push(&rules, &rule);
                rule = (rule_t){next, self, (i + 4) % 8};
                array_push(&rules, &rule);
            }
        }
    }
    *weights = array_set_weights(&rules);
    return rules;
}

static ptrdiff_t qrsort_partition(rule_t* arr, const ptrdiff_t low, const ptrdiff_t high)
{
    ptrdiff_t pivot = arr[high].self;
    ptrdiff_t i = (low - 1);
    for (ptrdiff_t j = low; j < high; ++j) {
        if ((ptrdiff_t)arr[j].self < pivot) {
            ++i;
            swap(arr[i], arr[j], rule_t);
        }
    }
    swap(arr[i + 1], arr[high], rule_t);
    return i + 1;
}

static void qrsort_dnq(rule_t* arr, const ptrdiff_t low, const ptrdiff_t high)
{
    if (low < high) {
        ptrdiff_t p = qrsort_partition(arr, low, high);
        qrsort_dnq(arr, low, p - 1);
        qrsort_dnq(arr, p + 1, high);
    }
}

static void rules_sort(const array_t* rulearr)
{
    return qrsort_dnq(rulearr->data, 0, (ptrdiff_t)rulearr->size - 1);
}

static void fbprnt(const bmp_t* bmp, const int x, const int y)
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


static size_t 
rules_search_self(const rule_t* rules, const size_t self)
{
    size_t i;
    for (i = 0; rules[i].self != self; ++i);
    return i;
}

static array_t
bmp_entropy_scan_at(const bmp_t* bmp,       const array_t* tiles,
                    const array_t* rules,   const int x,        const int y)
{
    array_t windx = array_create(sizeof(size_t));
    array_t indx = array_create(sizeof(size_t));
    const rule_t* r = rules->data;
    for (size_t i = 0; i < 8; ++i) {
        const int xpos = x + sidetable[i][0], ypos = y + sidetable[i][1];
        if (bmp_outside(bmp, xpos, ypos)) {
            continue;
        }
        
        const uint8_t* p = px_at(bmp, xpos, ypos);
        if (!p[3]) {
            continue;
        }
        
        size_t self = array_search(tiles, p);
        if (!self--) {
            continue;
        }
        
        const size_t rindex = rules_search_self(r, self), opside = (i + 4) % 8;
        for (size_t j = rindex; r[j].self == self; ++j) {
            if (r[j].side == opside) {
                array_push(&indx, &r[j].next);
            }
        }
        
        if (!indx.size) {
            continue;
        }

        if (!windx.size) {
            windx = indx;
            indx = array_create(sizeof(size_t));
        } else {
            const size_t* indices = windx.data;
            for (size_t j = 0; j < windx.size; ++j) {
                if (!array_search(&indx, indices + j)) {
                    array_remove(&windx, j--);
                }
            }
            array_clear(&indx);
        }
    }
    
    array_free(&indx);
    return windx;
}

static void
bmp_entropy_scan(   const bmp_t* bmp,       const array_t* tiles, 
                    const array_t* rules,   int* xptr,          int* yptr)
{
    size_t min = SIZE_MAX;
    for (size_t y = 0; y < bmp->height; ++y) {
        for (size_t x = 0; x < bmp->width; ++x) {
            const uint8_t* P = px_at(bmp, x, y);
            if (P[3]) {
                continue;
            }

            array_t w = bmp_entropy_scan_at(bmp, tiles, rules, x, y);
            if (!w.size) {
                continue;
            }

            if (w.size < min) {
                min = w.size;
                *xptr = x;
                *yptr = y;
            }

            array_free(&w);
        }
    }

    if (min == SIZE_MAX) {
        *xptr = -1;
        *yptr = -1;
    }
}

static void 
bmp_collapse(   const bmp_t* bmp,       const array_t* tiles, 
                const array_t* rules,   const array_t* weights,
                const int x,            const int y)
{
    array_t windx = array_create(sizeof(size_t));
    array_t indx = array_create(sizeof(size_t));
    array_t wtiles = array_reserve(sizeof(float), weights->size);
    memset(wtiles.data, 0, sizeof(float) * weights->size);
    wtiles.size = weights->size;

    const rule_t* r = rules->data;
    for (size_t i = 0; i < 8; ++i) {
        const int xpos = x + sidetable[i][0], ypos = y + sidetable[i][1];
        if (bmp_outside(bmp, xpos, ypos)) {
            continue;
        }

        const uint8_t* p = px_at(bmp, xpos, ypos);
        if (!p[3]) {
            continue;
        }

        size_t self = array_search(tiles, p);
        if (!self--) {
            continue;
        }
 
        const size_t rindex = rules_search_self(r, self), opside = (i + 4) % 8;
        for (size_t j = rindex; r[j].self == self; ++j) {
            if (r[j].side == opside) {
                array_push(&indx, &r[j].next);
            }
        }

        if (!windx.size) {
            windx = indx;
            indx = array_create(sizeof(size_t));
        } else {
            const size_t* indices = windx.data;
            for (size_t j = 0; j < windx.size; ++j) {
                if (!array_search(&indx, indices + j)) {
                    array_remove(&windx, j--);
                }
            }
            array_clear(&indx);
        }
    }

    array_free(&indx);
    assert(windx.size);

    size_t *wn = windx.data;
    const float *w = weights->data;
    float* fw = wtiles.data, sum = 0.0f; 
    for (size_t i = 0; i < windx.size; ++i) {
        const size_t index = wn[i];
        const float f = w[index];
        fw[index] += f;
        sum += f;
    }

    array_free(&windx);
    assert(sum != 0.0f);
    const float s = 1.0f / sum;
    
    const size_t rnd = rand() % RNDNUM;
    size_t n = 0, m = 0, found = 0;
    for (size_t i = 0; i < wtiles.size; ++i) {
        m += (size_t)((fw[i] * s) * RNDNUM);
        if (rnd >= n && rnd < m) {
            found = i + 1;
            break;
        }
        n = m;
    }

    array_free(&wtiles);
    assert(found);
    --found;

    memcpy(px_at(bmp, x, y), array_index(tiles, found), sizeof(Px));
    if (spxeKeyPressed(ESCAPE)) {
        exit(EXIT_SUCCESS);
    }
    
    fbprnt(bmp, 0, 0);
    if (!spxeRun(fb)) {
        exit(EXIT_FAILURE);
    }
}

static bmp_t
bmp_wfc(    const bmp_t* img,       const array_t* tiles, 
            const array_t* rules,   const array_t* weights) 
{
    bmp_t bmp = bmp_new(img->width, img->height, img->channels);
    memset(bmp.pixels, 0, bmp.width * bmp.height * bmp.channels);

#if 1
    memcpy(
        px_at(&bmp, 0, bmp.height - 1), 
        px_at(img, 0, img->height - 1),
        sizeof(Px) * bmp.width
    );
    memcpy(bmp.pixels, img->pixels, sizeof(Px) * bmp.width);
    for (size_t y = 0; y < bmp.height; ++y) {
        memcpy(px_at(&bmp, 0, y), px_at(img, 0, y), sizeof(Px));
        memcpy(px_at(&bmp, bmp.width - 1, y), px_at(img, img->width - 1, y), sizeof(Px));
    }
#else
    memcpy(bmp.pixels, img->pixels, sizeof(Px));
#endif

    int x, y;
    for (size_t i = 0; i < bmp.width * bmp.height; ++i) {
        bmp_entropy_scan(&bmp, tiles, rules, &x, &y);
        if (x == -1 || y == -1) {
            continue;
        }
        bmp_collapse(&bmp, tiles, rules, weights, x, y);
    }
    
    return bmp;
}

int main(int argc, char** argv)
{
    const char* path = NULL;
    int w = 200, h = 150;
    if (argc > 1) {
        path = argv[1];
    }
    if (argc > 2) {
        w = h = atoi(argv[1]);
    }
    if (argc > 3) {
        h = atoi(argv[2]);
    }

    if (!path) {
        printf("Missing input image file.\n");
        return EXIT_FAILURE;
    }

    srand(time(NULL));

    bmp_t bmp = bmp_load(path);
    if (!bmp.pixels) {
        printf("Could not load image file '%s'.\n", path);
        return EXIT_FAILURE;
    }
    freebmp = &bmp;
    
    fb = spxeStart("Wave Function", 800, 600, w, h);
    if (!fb) {
        printf("Could not init window.\n");
        bmp_free(&bmp);
        return EXIT_FAILURE;
    }

    freearr = array_create(sizeof(array_t));
    atexit(&wfc_free);

    fbprnt(&bmp, 0, bmp.height + 16);

    array_t weights, rweights;
    array_t tiles = tiles_get(&bmp, &weights);

    array_t rules = rules_get(&bmp, &tiles, &rweights);
    rules_sort(&rules);
    //rules_print(&rules);
    //tiles_print(&tiles);
    //weights_print(&rweights);
    //weights_print(&weights);

    array_push(&freearr, &tiles);
    array_push(&freearr, &weights);
    array_push(&freearr, &rules);
    array_push(&freearr, &rweights);

    bmp_t img = bmp_wfc(&bmp, &tiles, &rules, &weights);
    //freebmp = &img;

    while(spxeRun(fb)) {
        if (spxeKeyPressed(ESCAPE)) {
            break;
        }

        if (spxeKeyPressed(R)) {
            bmp_free(&img);
            img = bmp_wfc(&bmp, &tiles, &rules, &weights);
        }
    }

    bmp_free(&img);
    return EXIT_SUCCESS;
}
