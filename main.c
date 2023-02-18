/*  Copyright (c) 2023 Eugenio Arteaga A.

Permission is hereby granted, free of charge, to any 
person obtaining a copy of this software and associated 
documentation files (the "Software"), to deal in the 
Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to 
permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice 
shall be included in all copies or substantial portions
of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS
OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  */

/**************************
 * WAVE FUNCTION COLLAPSE *
 *************************/

#define WFC_APPLICATION
#include "wfc.h"
#define SPXE_APPLICATION
#include <spxe.h>
#include <imgtool.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#define WIDTH 100
#define HEIGHT 75
#define STRSIZE 1024

static char wfc_input_file[STRSIZE] = {0};

static void wfc_print_help(const char* name)
{
    printf("usage: %s\n", name);
    printf("<arg>\t\t: choose a file as input to analize its patterns\n");
    printf("-p\t\t: set periodic flag to true, output image warps edges\n");
    printf("-g\t\t: set ground flag, set output borders equal to input image\n");
    printf("-D\t\t: set output image width and height equal to window screen\n");
    printf("-w <arg>\t: set width of output image equal to <arg>\n");
    printf("-h <arg>\t: set height of output image equal to <arg>\n");
    printf("-s <arg>\t: set seed of wave function collapse randomess\n");
    printf("-N <arg>\t: set size of the analize patterns to NxN\n");
    printf("-S <arg>\t: 8 is the highest level of symmetry, 0 is no symmetry\n");
    printf("-H <arg>\t: heuristic to observe next pixel (0 = entropy, 1 = MRV, 2 = scanline)\n");
}

static void wfc_print_controls(void)
{
    printf("controls:\n");
    printf("ESCAPE\t\t: exit\n");
    printf("SPACE\t\t: collapse the model progresively\n");
    printf("LEFT CLICK\t: collapse pixel under mouse position\n");
    printf("R\t\t: reset model to unobserved state\n");
    printf("X\t\t: run current model until contradiction is met or image finished\n");
    printf("N\t\t: go to next seed and run model until contradiction is met or image finished\n");
    printf("P\t\t: save current output image as 'output.png' in current directory\n");
}

static void bmp_plot(Px* fb, const bmp_t* bmp, const int x, const int y)
{
    const int xpos = spxe.scrres.width / 2 - bmp->width / 2 + x;
    const int ypos = spxe.scrres.height / 2 - bmp->height / 2 + y;
    if (xpos >= 0 && xpos < spxe.scrres.width && 
        ypos >= 0 && ypos < spxe.scrres.height) {
        for (unsigned int y = 0; y < bmp->height; ++y) {
            const int width =   (int)bmp->width < spxe.scrres.width - xpos ? 
                                (int)bmp->width : 
                                spxe.scrres.width - xpos;
            memcpy(
                fb + (bmp->height - y - 1 + ypos) * spxe.scrres.width + xpos, 
                px_at(bmp, 0, y), 
                width * sizeof(Px)
            );
        }
    }
}

static void wfc_file_drop_input(GLFWwindow* window, int count, const char** paths)
{
    (void)window;
    if (count) {
        strncpy(wfc_input_file, paths[0], STRSIZE - 1);
    }
}

static void wfc_model_flush(
        const struct wfc_model* model, bmp_t* bmp, Px* fb, bool periodic)
{
    wfc_model_save(model, bmp->pixels, periodic);
    bmp_plot(fb, bmp, 0, 0);
}

int main(const int argc, const char** argv)
{
    const char* path = NULL;
    const int scrwidth = 100, scrheight = 75;
    int heuristic = WFC_ENTROPY, symmetry = 1;
    int width = 48, height = 48, N = 3, seed = time(NULL);
    bool periodic = false, ground = false, dimension = false;
    
    for (int i = 1; i < argc; ++i) {
        if (!strcmp("-help", argv[i]) || !strcmp("--help", argv[i])) {
            wfc_print_help(argv[0]);
            return EXIT_SUCCESS;
        }
        if (argv[i][0] == '-' && !argv[i][2]) {
            int* n = NULL;
            switch (argv[i][1]) {
                case 'p': { periodic = true; break; }
                case 'g': { ground = true; break; }
                case 'D': { dimension = true; break; }
                case 'w': { n = &width; break; }
                case 'h': { n = &height; break; }
                case 's': { n = &seed; break; }
                case 'N': { n = &N; break; }
                case 'S': { n = &symmetry; break; }
                case 'H': { n = &heuristic; break; }
                default: {
                    fprintf(
                        stderr, 
                        "Unrecognized option %s. See -help for more info\n", 
                        argv[i]
                    );
                    return EXIT_FAILURE;
                }
            }
            
            if (n) {
                if (i + 1 == argc) {
                    fprintf(
                        stderr, 
                        "Missing argument for option %s. See -help for more info\n",
                        argv[i]
                    );
                    return EXIT_FAILURE;
                }
                *n = atoi(argv[++i]);
            }
        }
        else path = argv[i];
    }

    if (!path) {
        fprintf(stderr, "Missing input image file. See -help for more information\n");
        return EXIT_FAILURE;
    }

    if (dimension) {
        width = scrwidth;
        height = scrheight;
    }

    bmp_t bmp = bmp_load(path);
    if (!bmp.pixels) {
        fprintf(stderr, "Could not load image file '%s'.\n", path);
        return EXIT_FAILURE;
    }

    srand(seed);

    bmp_t outbmp = bmp_new(width, height, 4);
    Px* fb = spxeStart("Wave Function Collapse", 800, 600, WIDTH, HEIGHT);
    glfwSetDropCallback(spxe.window, &wfc_file_drop_input);
    memset(fb, 125, WIDTH * HEIGHT * 4);
    bmp_plot(fb, &bmp, -33, 0);

    struct wfc_model model = wfc_model_create(
        bmp.pixels, 
        bmp.width, 
        bmp.height, 
        width, 
        height, 
        N, 
        periodic, 
        symmetry
    );

    wfc_model_clear(&model, periodic, ground);
    wfc_model_flush(&model, &outbmp, fb, periodic);
    wfc_print_controls();

    const int hx = spxe.scrres.width / 2, hy = spxe.scrres.height / 2;
    const int mhx = model.width / 2, mhy = model.height / 2;

    while (spxeRun(fb)) {

        if (spxeKeyPressed(ESCAPE)) {
            break;
        }
        else if (spxeKeyPressed(R)) {
            wfc_model_clear(&model, periodic, ground);
        }
        else if (spxeKeyPressed(P)) {
            bmp_write("output.png", &outbmp);
            printf("saved output image as 'output.png'\n");
        }
        else if (spxeKeyPressed(X)) {
            wfc_model_run(&model, heuristic, periodic, -1);
        }
        else if (spxeKeyPressed(N)) {
            wfc_model_clear(&model, periodic, ground);
            wfc_model_run(&model, heuristic, periodic, -1);
        }
        else if (spxeKeyDown(SPACE)) {
            wfc_model_step(&model, heuristic, periodic);
        }
        
        if (spxeMouseDown(LEFT)) {
            int x, y, dx, dy;
            spxeMousePos(&x, &y);
            dx = x - hx, dy = y - hy;
            if (dx >= -mhx && dy >= -mhx && dx < mhx && dy < mhy) {
                const int i = (model.height - (dy + mhy)) * model.width + dx + mhx;
                wfc_model_observe(&model, i);
                wfc_model_propagate(&model, periodic);
            }
        }

        if (wfc_input_file[0]) {
            bmp_t new_bmp = bmp_load(wfc_input_file);
            if (new_bmp.pixels) {
                bmp_free(&bmp);
                bmp = new_bmp;
                
                memset(fb, 125, WIDTH * HEIGHT * 4);
                bmp_plot(fb, &bmp, -33, 0);
                
                wfc_model_destroy(&model);
                model = wfc_model_create(
                    bmp.pixels, 
                    bmp.width, 
                    bmp.height, 
                    width, 
                    height,
                    N, 
                    periodic, 
                    symmetry
                );
                
                wfc_model_clear(&model, periodic, ground);
                wfc_model_flush(&model, &outbmp, fb, periodic);
            }
            wfc_input_file[0] = 0;
        }

        wfc_model_flush(&model, &outbmp, fb, periodic);
    }
    
    wfc_model_destroy(&model);
    bmp_free(&bmp);
    bmp_free(&outbmp);
    return spxeEnd(fb);
}
