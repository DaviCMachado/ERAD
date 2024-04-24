#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define image(x,y) pixels[y*width+x]
#define smooth(x,y) filtered[y*width+x]

typedef struct {
    char red, green, blue, alpha;
} RGBA;

int main(int argc, char *argv[]) {

    FILE *in;
    FILE *out;

    in = fopen("image.in", "rb");
    if (in == NULL) {
        perror("image.in");
        exit(EXIT_FAILURE);
    }

    out = fopen("image.out", "wb");
    if (out == NULL) {
        perror("image.out");
        exit(EXIT_FAILURE);
    }

    short width, height;

    fread(&width, sizeof (width), 1, in);
    fread(&height, sizeof (height), 1, in);

    fwrite(&width, sizeof (width), 1, out);
    fwrite(&height, sizeof (height), 1, out);

    RGBA *pixels = (RGBA *) malloc(height * width * sizeof (RGBA));
    RGBA *filtered = (RGBA *) malloc(height * width * sizeof (RGBA));

    if (pixels == NULL || filtered == NULL) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    int DY[] = {-2, -2, -2, -2, -2, -1, -1, -1, -1, -1, +0, +0, +0, +0, +0, +1, +1, +1, +1, +1, +2, +2, +2, +2, +2};
    int DX[] = {-2, -1, +0, +1, +2, -2, -1, +0, +1, +2, -2, -1, +0, +1, +2, -2, -1, +0, +1, +2, -2, -1, +0, +1, +2};
    int x, y, d, dx, dy, i;

    do {
        if (!fread(pixels, height * width * sizeof (RGBA), 1, in))
            break;

        for (y = 0; y < height; y++) {
            for (x = 0; x < width; x++) {
                for (i = 0; i < 4; i++) {
                    long long int sum = 0;
                    for (d = 0; d < 25; d++) {
                        dx = x + DX[d];
                        dy = y + DY[d];
                        if (dx >= 0 && dx < width && dy >= 0 && dy < height) {
                            sum += *(((char*) (&image(dx, dy))) + i);
                        }
                    }
                    (*(((char*) (&smooth(x, y)) + i))) = sum / 25;
                }
            }
        }

        fwrite(filtered, height * width * sizeof (RGBA), 1, out);

    } while (!feof(in));

    free(pixels);
    free(filtered);

    fclose(out);
    fclose(in);

    return EXIT_SUCCESS;
}
