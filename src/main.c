#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "proj_utils.h"

// === Util Function Prototypes ===

/**
 * @brief Calculates the Euclidean distance between two points.
 *
 * Time complexity: O(1)
 * Space complexity: O(1)
 *
 * @param p1 Pointer to the first point.
 * @param p2 Pointer to the second point.
 * @return The Euclidean distance between p1 and p2.
 */
double points_dist(Point* p1, Point* p2);

/**
 * @brief Finds the closest pair of points using the brute-force approach.
 *
 * The brute-force approach consists in checking the distance of every point
 * with everything other point. It uses two nested loops.
 *
 * Time complexity: O(n²)
 * Space complexity: O(1)
 *
 * @param points_arr Array of points.
 * @param points_num Number of points in the array.
 * @return A pointer to an array of integers representing the 0-indexed indices
 *         of the closets pair of points. E.g., { 1, 5 }.
 */
int* bf_closest_pair(Point* points_arr, int points_num);

/**
 * @brief Find closest pair of points using the divide-and-conquer approach.
 *
 * The divide-and-conquer approach consists in dividing the points into two
 * halves, finding the closest pair in each half, and then checking if there is
 * a closer pair that crosses the division line.
 *
 * Time complexity: O(n log n)
 * Space complexity: O(n)
 *
 * @param points_arr Array of points.
 * @param points_num Number of points in the array.
 * @return A pointer to an array of integers representing the 0-indexed indices
 *         of the closest pair of points. E.g., { 1, 5 }.
 */
int* dnc_closest_pair(Point* points_arr, int points_num);

// Heap Sort

/**
 * Return the index of the parent of the given index.
 * 
 * It performs the operation by shifting the number to the right by one bit.
 * 
 * @param i Index of the element in the heap.
 * @return Index of the parent of the element in i.
 */
inline int heap_parent(int i) { return (i - 1) >> 1; }

/**
 * Get the index of the left child of a node in a heap given the element index.
 * 
 * It performs the operation by shifting the index one bit to the left.
 * 
 * @param i Index of the element in the heap.
 * @return Index of the left child in the heap.
 */
inline int heap_left(int i) {return (i << 1) + 1; }

/**
 * Get the index of the right child of a node in a heap given the index.
 * 
 * The operation is performed shifting the index by the left one bit and adding
 * one.
 * 
 * @param i Index of the element in the heap.
 * @return Index of the right child.
 */
inline int heap_right(int i) { return (i << 1) + 2; }

/**
 * Maintain the max-heap property after adding an item to the array.
 * 
 * This function assume that the left and right subtrees are max-heaph (they
 * maintain the property of a max-heap). It performs the operation in-place.
 * 
 * Note that the function orders the points by the x-coordinates.
 * 
 * @param heap Pointer to the root of the max-heap.
 * @param s Size of the max-heap.
 * @param i Index of the element to heapify.
 */
void max_heapify_points_by_x(Point* heap, int s, int i);

/**
 * Given an array, this function generates a max heap in place.
 * 
 * Note that this functions creates a max heap using the array of points and
 * their x value as reference.
 * 
 * @param heap The contiguous array to convert to a max heap.
 * @param s Size of the input array.
 */
void build_max_heap_points_by_x(Point* arr, int s);

/**
 * Sort an array of points by their x coordinate using the heap sort algorithm.
 * 
 * Time compexity: O(n log n)
 * Space complexity: O(1)
 * 
 * @param arr Array of points to sort.
 * @param s Size of the array.
 */
void heap_sort_points_by_x(Point* arr, int s);

// === Implementation ===

int main(void) {
  // Get the points (generate or retrieve from file).
  for (int i = 0; i < 10; i++) {
    int points_num = 1000 * (i + 1);
    printf("\n\n=== [Main] Generating %d points.\n", points_num);

    for (int j = 0; j < 10; j++) { // 10 variations of points_num.
      char filename[50];
      snprintf(filename, sizeof(filename), "build/bin/data/points-%d-var-%d.txt", points_num, j + 1);
      Point* points_arr = get_points(points_num, filename);
      if (!points_arr) {
        fprintf(stderr, "[Main] Error: Failed to generate points.\n");
        return EXIT_FAILURE;
      }

      heap_sort_points_by_x(points_arr, points_num);

      printf("\nPrinting 10 points sorted by their x-coordinate:\n");
      for (int k = 0; k < 10; k++) printf("%2d. (%d, %d)\n",
        k + 1,
        (points_arr + k)->x,
        (points_arr + k)->y
      );

      free_points(points_arr); // Free the allocated memory for points.
      // printf("[Main] Successfully generated %d points.\n", points_num);
    }
  }
  return 0;
}

double points_dist(Point* p1, Point* p2) {
  double p1x = (double) p1->x;
  double p1y = (double) p1->y;
  double p2x = (double) p2->x;
  double p2y = (double) p2->y;

  return sqrt( // √[(x₂ - x₁)² + (y₂ - y₁)²]
    pow((double) (p2x - p1x), 2.0) +
    pow((double) (p2y - p1y), 2.0)
  );
}

int* bf_closest_pair(Point* points_arr, int points_num) {
  if (points_num < 2) return NULL;

  int* pairs = calloc(2, sizeof(int));
  pairs[0] = 0, pairs[1] = 1;

  // d_min is not inf, the min distance must be the dist of one pair.
  double d_min = points_dist(points_arr, points_arr + 1);
  double d_temp = d_min;

  for (int i = 0; i < points_num - 1; i++) {
    for (int j = i + 1; j < points_num; j++){
      d_temp = points_dist(points_arr + i, points_arr + j);
      if (d_temp < d_min) {
        d_min = d_temp;
        pairs[0] = i;
        pairs[1] = j;
      }
    }
  }

  return pairs;
}

void max_heapify_points_by_x(Point* heap, int s, int i) {
  int l = heap_left(i);
  int r = heap_right(i);
  int largest;

  if (l < s && (heap + l)->x > (heap + i)->x) {
    largest = l;
  } else {
    largest = i;
  }

  if (r < s && (heap + r)->x > (heap + largest)->x) {
    largest = r;
  }

  if (largest != i) {
    Point temp = *(heap + i);
    *(heap + i) = *(heap + largest);
    *(heap + largest) = temp;

    max_heapify_points_by_x(heap, s, largest);
  }
}

void build_max_heap_points_by_x(Point* arr, int s) {
  for (int i = (s >> 2); i > 0; i--) {
    max_heapify_points_by_x(arr, s, i);
  }
}

void heap_sort_points_by_x(Point* arr, int s) {
  build_max_heap_points_by_x(arr, s);
  int heap_s = s;

  for (int i = s - 1; i > 0; i--) {
    Point temp = *arr;
    *arr = *(arr + i);
    *(arr + i) = temp;

    heap_s--;
    max_heapify_points_by_x(arr, heap_s, 0);
  }
}