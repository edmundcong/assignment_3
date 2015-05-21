#ifndef MATRIX_H
#define MATRIX_H
#define index(x,y) ((x)+g_width*(y))

#include <stdint.h>

struct scalar_struct {
    const uint32_t *matrix;
    uint32_t *result;
    ssize_t scalar;
    ssize_t array_start;
    ssize_t array_limit;
};

struct mul_struct {
  const uint32_t *mat_a;
  const uint32_t *mat_b;
  uint32_t *result;
  size_t index;
  size_t chunk;
  size_t threads;
  size_t width;
};

struct cmp_struct {
  const uint32_t *matrix;
  uint32_t *result;
  size_t index;
  // size_t value;
  ssize_t array_start;
  ssize_t array_limit;
};

/* utility functions */

uint32_t fast_rand(void);

void set_seed(uint32_t value);
void set_nthreads(ssize_t count);
void set_dimensions(ssize_t width);

void display(const uint32_t *matrix);
void display_row(const uint32_t *matrix, ssize_t row);
void display_column(const uint32_t *matrix, ssize_t column);
void display_element(const uint32_t *matrix, ssize_t row, ssize_t column);

/* matrix operations */

uint32_t *new_matrix(void);
uint32_t *identity_matrix(void);
uint32_t *random_matrix(uint32_t seed);
uint32_t *uniform_matrix(uint32_t value);
uint32_t *sequence_matrix(uint32_t start, uint32_t step);

uint32_t *cloned(const uint32_t *matrix);
uint32_t *reversed(const uint32_t *matrix);
uint32_t *transposed(const uint32_t *matrix);

uint32_t *scalar_add(const uint32_t *matrix, uint32_t scalar);
uint32_t *scalar_mul(const uint32_t *matrix, uint32_t scalar);
uint32_t *matrix_pow(const uint32_t *matrix, uint32_t exponent);
uint32_t *matrix_add(const uint32_t *matrix_a, const uint32_t *matrix_b);
uint32_t *matrix_mul(const uint32_t *matrix_a, const uint32_t *matrix_b);

/* worker functions */

// static void *worker_matrix_mul(void *matrix);

/* compute operations */

uint32_t get_sum(const uint32_t *matrix);
uint32_t get_trace(const uint32_t *matrix);
uint32_t get_minimum(const uint32_t *matrix);
uint32_t get_maximum(const uint32_t *matrix);
uint32_t get_frequency(const uint32_t *matrix, uint32_t value);


#endif
