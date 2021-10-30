#include <stdio.h>
#include <string.h>

#include "rtga/rtga.h"

int main(void) {
    TgaImage image;
    int success;

    memset(&image, 0, sizeof(image));

    success = tga_alloc(UNCOMPRESSED_TRUE_COLOR_IMAGE, 600, 400, 24, &image);
    if (success != TGA_SUCCESS) {
        printf("Error %i occured while allocating image.\n", success);
        return 1;
    }

    success = tga_write_file(&image, "raytracer.tga");
    if (success != TGA_SUCCESS) {
        printf("Error occured while writing image to a file.\n");
        return 2;
    }

    tga_free(&image);
    return 0;
}
