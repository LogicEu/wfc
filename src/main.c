#include "wfc.h"
#include <stdio.h>
#include <time.h>
#define SPXE_APPLICATION
#include <spxe.h>

static void wfc_plot(Px* fb, const bmp_t* bmp, const int x, const int y)
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

static void wfc_model_flush(const wfc_model_t* model, bmp_t* bmp, Px* fb, bool periodic)
{
    wfc_model_save(model, bmp->pixels, periodic);
    wfc_plot(fb, bmp, 0, 0);
}

int main(const int argc, const char** argv)
{
    const char* path = NULL;
    if (argc > 1) {
        path = argv[1];
    }

    if (!path) {
        printf("Missing input image file.\n");
        return EXIT_FAILURE;
    }

    bmp_t bmp = bmp_load(path);
    if (!bmp.pixels) {
        printf("Could not load image file '%s'.\n", path);
        return EXIT_FAILURE;
    }

    const unsigned int seed = (unsigned int)time(NULL), w = 48, h = 48;
    srand(seed);

    bmp_t outbmp = bmp_new(w, h, 4);
    Px* fb = spxeStart("Wave Function Collapse", 800, 600, 100, 75);
    wfc_plot(fb, &bmp, -33, 0);

    bool periodic = true;
    bool ground = false;
    int heuristic = WFC_ENTROPY;

    wfc_model_t model = wfc_model_create(&bmp, w, h, 5, true, 2);
    wfc_model_clear(&model, periodic, ground);
    wfc_model_flush(&model, &outbmp, fb, periodic);

    const int hx = spxe.scrres.width / 2, hy = spxe.scrres.height / 2;
    const int mhx = model.width / 2, mhy = model.height / 2;

    while (spxeRun(fb)) {

        if (spxeKeyPressed(ESCAPE)) {
            break;
        }
        else if (spxeKeyPressed(R)) {
            wfc_model_clear(&model, periodic, ground);
        }
        else if (spxeKeyPressed(X)) {
            wfc_model_run(&model, heuristic, periodic, -1);
        }
        else if (spxeKeyPressed(N)) {
            wfc_model_clear(&model, periodic, ground);
            wfc_model_run(&model, heuristic, periodic, -1);
        }
        else if (spxeKeyDown(SPACE) || spxeKeyPressed(Z)) {
            wfc_model_step(&model, heuristic, periodic);
        }
        
        if (spxeMouseDown(1)) {
            int x, y, dx, dy;
            spxeMousePos(&x, &y);
            dx = x - hx, dy = y - hy;
            if (dx >= -mhx && dy >= -mhx && dx < mhx && dy < mhy) {
                const int i = (model.height - (dy + mhy)) * model.width + dx + mhx;
                wfc_model_observe(&model, i);
                wfc_model_propagate(&model, periodic);
            }
        }

        wfc_model_flush(&model, &outbmp, fb, periodic);
    }
    
    wfc_model_destroy(&model);
    bmp_free(&bmp);
    bmp_free(&outbmp);
    return spxeEnd(fb);
}