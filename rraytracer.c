#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "rtga/rtga.h"

#define HEIGHT 400
#define WIDTH 600

int main(void) {
    TgaImage image;
    int success;

    // Set image struct to all 0s
    memset(&image, 0, sizeof(image));

    // Allocate image
    success = tga_alloc(UNCOMPRESSED_TRUE_COLOR_IMAGE, WIDTH, HEIGHT, 24, &image);
    if (success != TGA_SUCCESS) {
        printf("Error %i occured while allocating image.\n", success);
        return 1;
    }

    // Trace rays for each pixel in image
    for (int y = 0; y < HEIGHT; ++y) {
        printf("\r\e[2K");
        printf("\rScanlines remaining: %i", HEIGHT - y);
        fflush(stdout);
        for (int x = 0; x < WIDTH; ++x) {
            TgaColor color = COLOR24((uint8_t)((float)x / WIDTH * 255.0), (uint8_t)((float)y / HEIGHT * 255.0), 128);
            tga_set_pixel(&image, x, y, color);
        }
    }
    printf("\r\e[2K");
    printf("\rRaytracer finished!\n");


    // Write image to file
    success = tga_write_file(&image, "raytracer.tga");
    if (success != TGA_SUCCESS) {
        printf("Error occured while writing image to a file.\n");
        return 2;
    }

    // Free image
    tga_free(&image);

    return 0;
}
