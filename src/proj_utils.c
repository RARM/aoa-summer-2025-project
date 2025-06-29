#include "proj_utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>

Point* generate_random_points(size_t num_points, const char* filename) {
  Point* arr = malloc(sizeof(Point) * num_points);
  if (!arr) fprintf(stderr, "[Project Utils] Error: Can't allocate memory for points.\n");

  FILE* fp;
  if (filename) {
    fp = fopen(filename, "w");
    if (!fp) fprintf(stderr, "[Project Utils] Error: Couldn't open the file: \"%s\"\n", filename);
  }
  
  srand(time(NULL));
  for (size_t i = 0; i < num_points; i++) {
    (arr + i)->x = rand() | (rand() % 2) * INT_MIN;
    (arr + i)->y = rand() | (rand() % 2) * INT_MIN;

    if (filename && fp) fprintf(fp, "%d %d\n", (arr + i)->x, (arr + i)->y);
  }

  if (filename && fp) fclose(fp);

  return arr;
}

void free_points(Point* points) {
  free(points);
}