#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define image(x,y) pixels[y*width+x]

typedef struct {
    char red, green, blue, alpha;
} RGBA;

int main(void) {
    srand(time(NULL));
    unsigned long long int size;
    unsigned long long int x, y;

    FILE *out;

    out = fopen("image.in", "wb");
    if (out == NULL) {
        perror("image.in");
        exit(EXIT_FAILURE);
    }

    short width, height;
    width = 7680;
    height = 4320;

    fwrite(&width, sizeof (width), 1, out);
    fwrite(&height, sizeof (height), 1, out);

    RGBA *pixels = (RGBA*) malloc(width * height * sizeof (RGBA));
    if (pixels == NULL) {
        perror("malloc pixels");
        exit(EXIT_FAILURE);
    }

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            image(x, y).red = '1';
            image(x, y).green = '1';
            image(x, y).blue = '1';
            image(x, y).alpha = '1';
        }
    }
    fwrite(pixels, height * sizeof (RGBA) * width, 1, out);

    fflush(out);
    fclose(out);
    free(pixels);

    return EXIT_SUCCESS;
}
