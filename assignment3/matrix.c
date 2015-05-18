#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>
#include <pthread.h>
#include <inttypes.h>
#include "matrix.h"

static uint32_t g_seed = 0;

static ssize_t g_width    = 0;
static ssize_t g_height   = 0;
static ssize_t g_elements = 0;

static ssize_t g_nthreads = 1;

int temp = 0;


////////////////////////////////
///     UTILITY FUNCTIONS    ///
////////////////////////////////

/**
 * Returns pseudorandom number determined by the seed
 */
uint32_t fast_rand(void) {

    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & 0x7FFF;
}

/**
 * Sets the seed used when generating pseudorandom numbers
 */
void set_seed(uint32_t seed) {

    g_seed = seed;
}

/**
 * Sets the number of threads available
 */
void set_nthreads(ssize_t count) {

    g_nthreads = count;
}

/**
 * Sets the dimensions of the matrix
 */
void set_dimensions(ssize_t order) {

    g_width  = order;
    g_height = order;

    g_elements = g_width * g_height;
}

/**
 * Displays given matrix
 */
void display(const uint32_t *matrix) {

    for (ssize_t y = 0; y < g_height; y++) {
        for (ssize_t x = 0; x < g_width; x++) {
            if (x > 0) printf(" ");
            printf("%" PRIu32, matrix[y * g_width + x]);
        }

        printf("\n");
    }
}

/**
 * Displays given matrix row
 */
void display_row(const uint32_t *matrix, ssize_t row) {

    for (ssize_t x = 0; x < g_width; x++) {
        if (x > 0) printf(" ");
        printf("%" PRIu32, matrix[row * g_width + x]);
    }

    printf("\n");
}

/**
 * Displays given matrix column
 */
void display_column(const uint32_t *matrix, ssize_t column) {

    for (ssize_t y = 0; y < g_height; y++) {
        printf("%" PRIu32 "\n", matrix[y * g_width + column]);
    }
}

/**
 * Displays the value stored at the given element index
 */
void display_element(const uint32_t *matrix, ssize_t row, ssize_t column) {

    printf("%" PRIu32 "\n", matrix[row * g_width + column]);
}

////////////////////////////////
///   MATRIX INITALISATIONS  ///
////////////////////////////////

/**
 * Returns new matrix with all elements set to zero
 */
uint32_t *new_matrix(void) {

    return calloc(g_elements, sizeof(uint32_t));
}

/**
 * Returns new identity matrix
 */
uint32_t *identity_matrix(void) {

    uint32_t *matrix = new_matrix();

    for (ssize_t i = 0; i < g_width; i++) {
        matrix[i * g_width + i] = 1;
    }

    return matrix;
}

/**
 * Returns new matrix with elements generated at random using given seed
 */
uint32_t *random_matrix(uint32_t seed) {

    uint32_t *matrix = new_matrix();

    set_seed(seed);

    for (ssize_t i = 0; i < g_elements; i++) {
        matrix[i] = fast_rand();
    }

    return matrix;
}

/**
 * Returns new matrix with all elements set to given value
 */
uint32_t *uniform_matrix(uint32_t value) {

    uint32_t *matrix = new_matrix();

    for (ssize_t i = 0; i < g_elements; i++) {
        matrix[i] = value;
    }

    return matrix;
}

/**
 * Returns new matrix with elements in sequence from given start and step
 */
uint32_t *sequence_matrix(uint32_t start, uint32_t step) {

    uint32_t *matrix = new_matrix();
    uint32_t current = start;

    for (ssize_t i = 0; i < g_elements; i++) {
        matrix[i] = current;
        current  += step;
    }

    return matrix;
}

////////////////////////////////
///     MATRIX OPERATIONS    ///
////////////////////////////////

/**
 * Returns new matrix with elements cloned from given matrix
 */
uint32_t *cloned(const uint32_t *matrix) {

    uint32_t *result = new_matrix();

    for (ssize_t i = 0; i < g_elements; i++) {
        result[i] = matrix[i];
    }

    return result;
}

/**
 * Returns new matrix with elements ordered in reverse
 */
uint32_t *reversed(const uint32_t *matrix) {

    uint32_t *result = new_matrix();

    for (ssize_t i = 0; i < g_elements; i++) {
        result[i] = matrix[g_elements - 1 - i];
    }

    return result;
}

/**
 * Returns new transposed matrix
 */
uint32_t *transposed(const uint32_t *matrix) {

    uint32_t *result = new_matrix();

    for (ssize_t y = 0; y < g_height; y++) {
        for (ssize_t x = 0; x < g_width; x++) {
            result[x * g_width + y] = matrix[y * g_width + x];
        }
    }

    return result;
}

void *scalar_add_worker(void *scalar_add_struct){

    struct scalar_add_struct m_row = *(struct scalar_add_struct *)scalar_add_struct;//casting

    const int start = m_row.index*m_row.array_limit;
    //e.g. 0 < 16, 8 < 16,
    for (int i = start; i < g_elements; i++) {
        m_row.result[i] = m_row.matrix[i]+m_row.scalar;
    }

    return NULL;

}

/**
 * Returns new matrix with scalar added to each element
 */
uint32_t *scalar_add(const uint32_t *matrix, uint32_t scalar) {

    uint32_t *result = new_matrix();

    struct scalar_add_struct *m_row = malloc(sizeof(struct scalar_add_struct) * g_nthreads);


    int array_limit = g_elements/g_nthreads;

    pthread_t t_id[g_nthreads]; // store thread IDs

    for (ssize_t i = 0; i < g_nthreads; i++) {
        m_row[i].matrix      = matrix;
        m_row[i].result      = result;
        m_row[i].index       = i;
        m_row[i].scalar      = scalar;
        m_row[i].array_limit = array_limit;
    }

    for (ssize_t i = 0; i < g_nthreads; i++) {
        pthread_create(t_id+i, NULL, scalar_add_worker, &m_row[i]);
    }


    for (ssize_t i = 0; i < g_nthreads; ++i) {
        pthread_join(t_id[i], NULL);
    }

    free(m_row);

    return result;
}

/**
 * Returns new matrix with scalar multiplied to each element
 */
uint32_t *scalar_mul(const uint32_t *matrix, uint32_t scalar) {

    uint32_t *result = new_matrix();

    for (ssize_t i = 0; i < g_elements; i++) {
        result[i] = matrix[i] * scalar;
    }

    return result;
}

/**
 * Returns new matrix with elements added at the same index
 */
uint32_t *matrix_add(const uint32_t *matrix_a, const uint32_t *matrix_b) {

    uint32_t *result = new_matrix();

    for (ssize_t i = 0; i < g_elements; i++) {
        result[i] = matrix_a[i] + matrix_b[i];
    }

    return result;
}

static void *worker_matrix_mul(void *matrix){
    struct matrix_mul_struct res_mat = *(struct matrix_mul_struct *)matrix;

    const size_t start = res_mat.index * res_mat.chunk;
    const size_t end   = res_mat.index == g_nthreads - 1 ? g_width : (res_mat.index + 1)*res_mat.chunk;


    for (size_t i = start; i < end; ++i) {
        for (size_t j = 0; j < g_width; j++) {
            for (size_t k = 0; k < g_width; k++) {
                res_mat.result[IDX(j, i)] += res_mat.mat_a[IDX(k, i)]*res_mat.mat_b[IDX(j, k)];
            }
        }
    }
    return NULL;
}

/**
 * Returns new matrix, multiplying the two matrices together
 */
uint32_t *matrix_mul(const uint32_t *matrix_a, const uint32_t *matrix_b) {

    struct matrix_mul_struct *res_mat = malloc(sizeof(struct matrix_mul_struct)*g_nthreads);

    uint32_t *result = new_matrix();

    for (size_t i = 0; i < g_nthreads; ++i) {
        res_mat[i] = (struct matrix_mul_struct){
            .mat_a  = matrix_a,
            .mat_b  = matrix_b,
            .result = result,
            .index  = i,
            .chunk  = g_width/g_nthreads,
        };
    }

    pthread_t t_id[g_nthreads]; //creating array to store thread IDs

    for (size_t i = 0; i < g_nthreads; ++i) {
        //passing thread id + i, no attributes, function pointer, matrix struct
        pthread_create(t_id + i, NULL, worker_matrix_mul, res_mat + i);
    }

    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_join(t_id[i], NULL);
    }

    free(res_mat);
    return result;

}

/**
 * Returns new matrix, powering the matrix to the exponent
 */
uint32_t *matrix_pow(const uint32_t *matrix, uint32_t exponent) {

    uint32_t *result = new_matrix();

    if (exponent == 0) {
        return result = identity_matrix();
    }

    for (int i = 0; i < g_elements; i++) {
        result[i] = matrix[i];
    }


    if (exponent == 1) {
        return result;
    }

    if (exponent == 2) {
        return result = matrix_mul(matrix, matrix);
    }

    for (int i = 0; i < exponent-1; i++) { //-1 because ^1 is itself
        result = matrix_mul(result, matrix);
    }


    return result;


}

////////////////////////////////
///       COMPUTATIONS       ///
////////////////////////////////

/**
 * Returns the sum of all elements
 */
uint32_t get_sum(const uint32_t *matrix) {

    uint32_t sum = 0;

    for (int i = 0; i < g_elements; i++) {
        sum += matrix[i];
    }

    return sum;
}

/**
 * Returns the trace of the matrix
 */
uint32_t get_trace(const uint32_t *matrix) {

    uint32_t trace = 0;

    for (int i = 0; i < g_height; i++) {
        trace += matrix[i*g_width+i];  //row we're at + i places in
    }

    return trace;
}

/**
 * Returns the smallest value in the matrix
 */
uint32_t get_minimum(const uint32_t *matrix) {

    temp = matrix[0];

    for (int i = 1; i < g_elements; i++) {
        if (matrix[i] < temp) {
            temp = matrix[i];
        }

    }

    return temp;


}

/**
 * Returns the largest value in the matrix
 */
uint32_t get_maximum(const uint32_t *matrix) {

    temp = 0;

    for (int i = 0; i < g_elements; i++) {
        if (matrix[i] > temp) {
            temp = matrix[i];
        }

    }

    return temp;
}

/**
 * Returns the frequency of the value in the matrix
 */
uint32_t get_frequency(const uint32_t *matrix, uint32_t value) {

    uint32_t count = 0;

    for (int i = 0; i < g_elements; i++) {
        if (matrix[i] == value) {
            count++;
        }
    }

    return count;
}
