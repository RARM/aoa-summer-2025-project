#include "proj_utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>

Point* get_points(size_t num_points, const char* filename) {
  if (!filename) { // If no filename, generate random points without saving.
    return generate_random_points(num_points, NULL);
  }

  FILE* fp = fopen(filename, "r");
  if (fp) {
    Point* points = malloc(sizeof(Point) * num_points);
    if (!points) {
      fprintf(stderr, "[Project Utils] Error: Can't allocate memory for points.\n");
      fclose(fp);
      return NULL;
    }
    
    for (size_t i = 0; i < num_points; i++) {
      if (fscanf(fp, "%d %d", &points[i].x, &points[i].y) != 2) {
        fprintf(stderr, "[Project Utils] Error: Failed to read point from file.\n");
        free(points);
        fclose(fp);
        return NULL;
      }
    }
    
    fclose(fp);
    return points;
  } else {
    return generate_random_points(num_points, filename);
  }
}

Point* generate_random_points(size_t num_points, const char* filename) {
  Point* arr = malloc(sizeof(Point) * num_points);
  if (!arr) fprintf(stderr, "[Project Utils] Error: Can't allocate memory for points.\n");

  FILE* fp;
  if (filename) {
    // FIXME: The program fails if the directory does not exist.
    fp = fopen(filename, "w");
    if (!fp) fprintf(stderr, "[Project Utils] Error: Couldn't open the file: \"%s\"\n", filename);
  }

  srand(time(NULL));
  
  // Generate unique random points.
  size_t generated = 0;
  while (generated < num_points) {
    int x = rand() | (rand() % 2) * INT_MIN;
    int y = rand() | (rand() % 2) * INT_MIN;

    // Check for uniqueness.
    int unique = 1;
    for (size_t j = 0; j < generated; j++) {
      if (arr[j].x == x && arr[j].y == y) {
        unique = 0;
        break;
      }
    }
    if (!unique) continue;

    arr[generated].x = x;
    arr[generated].y = y;
    if (filename && fp) fprintf(fp, "%d %d\n", x, y);
    generated++;
  }

  if (filename && fp) fclose(fp);

  return arr;
}

void free_points(Point* points) {
  free(points);
}