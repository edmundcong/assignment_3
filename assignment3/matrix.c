#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>
#include <pthread.h>
#include <inttypes.h>
#include <string.h>
#include "matrix.h"

static uint32_t g_seed = 0;

static ssize_t g_width    = 0;
static ssize_t g_height   = 0;
static ssize_t g_elements = 0;

static ssize_t g_nthreads = 1;

uint32_t temp = 0;


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
inline uint32_t *cloned(const uint32_t *matrix) {

    uint32_t *result = new_matrix();

    memcpy(result, matrix, g_elements * sizeof(*result));

    return result;
}

/**
 * Returns new matrix with elements ordered in reverse
 */
inline uint32_t *reversed(const uint32_t *matrix) {

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
            result[index(y, x)] = matrix[index(x, y)];
        }
    }

    return result;
}

/**
 * Worker function for scalar add.
 * Performs scalar addition on each entry using pthreads.
 */
void *scalar_add_worker(void *add_struct){

    struct scalar_struct m_row = *(struct scalar_struct *)add_struct;//casting

    for (int i = m_row.array_start; i < m_row.array_limit; i++) {
        m_row.result[i] = m_row.matrix[i]+m_row.scalar;
    }

    return NULL;

}

/**
 * Returns new matrix with scalar added to each element
 */
uint32_t *scalar_add(const uint32_t *matrix, uint32_t scalar) {

    uint32_t *result = new_matrix();

    //structs allow us to pass more than 1 variable to the worker function
    struct scalar_struct *m_row = malloc(sizeof(struct scalar_struct) * g_nthreads);

    pthread_t t_id[g_nthreads]; // store thread IDs

    int chunk = g_elements/g_nthreads;

    for (ssize_t i = 0; i < g_nthreads; i++) {
        m_row[i].matrix      = matrix;
        m_row[i].result      = result;
        m_row[i].scalar      = scalar;
        m_row[i].array_start = i*chunk;
        m_row[i].array_limit = (i*chunk) + chunk;
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
 * Worker function for scalar multiply.
 * Performs scalar multiplication on each entry using pthreads.
 */
void *scalar_mul_worker(void *mult_struct){

    //cast the void struct as a scalar_struct data type
    struct scalar_struct m_row = *(struct scalar_struct *)mult_struct;

    for (int i = m_row.array_start; i < m_row.array_limit; i++) {
        m_row.result[i] = m_row.matrix[i]*m_row.scalar;
    }

    return NULL;

}

/**
 * Returns new matrix with scalar multiplied to each element
 */
uint32_t *scalar_mul(const uint32_t *matrix, uint32_t scalar) {

    uint32_t *result = new_matrix();

    struct scalar_struct *m_row = malloc(sizeof(struct scalar_struct) * g_nthreads);

    pthread_t t_id[g_nthreads]; // store thread IDs

    int chunk = g_elements/g_nthreads;

    for (ssize_t i = 0; i < g_nthreads; i++) {
        m_row[i].matrix      = matrix;
        m_row[i].result      = result;
        m_row[i].scalar      = scalar;
        m_row[i].array_start = i*chunk;
        m_row[i].array_limit = (i*chunk) + chunk;
    }

    for (ssize_t i = 0; i < g_nthreads; i++) {
        pthread_create(t_id+i, NULL, scalar_mul_worker, &m_row[i]);
    }

    for (ssize_t i = 0; i < g_nthreads; ++i) {
        pthread_join(t_id[i], NULL);
    }

    free(m_row);

    return result;

}

/**
 * Returns new matrix with elements added at the same index
 */
inline uint32_t *matrix_add(const uint32_t *matrix_a, const uint32_t *matrix_b) {

    uint32_t *result = new_matrix();

    int res = 0;

    for (ssize_t i = 0; i < g_elements; i++) {

        res       = matrix_a[i] + matrix_b[i];
        result[i] = res;
    }

    return result;
}

/**
 * Worker function for scalar multiply.
 * Performs scalar multiplication on each entry using pthreads.
 */
inline static void *worker_matrix_mul(void *mult_matrix){
    struct mul_struct res_mat = *(struct mul_struct *)mult_matrix;

    const size_t start = res_mat.index * res_mat.chunk;
    //if index == threads-1 then end is equal to width, otherwise equal to (index+1)*chunk
    const size_t end = res_mat.index == res_mat.threads - 1 ? res_mat.width : (res_mat.index + 1)*res_mat.chunk;

    int res = 0;

    //performs matrix multiplicaton on two matrices and assigns to result
    for (size_t i = start; i < end; ++i) {
        for (size_t j = 0; j < res_mat.width; j++) {
            res = 0;
            for (size_t k = 0; k < res_mat.width; k++) {
                res                        += res_mat.mat_a[index(k, i)]*res_mat.mat_b[index(j, k)];
                res_mat.result[index(j, i)] = res;
            }
        }
    }
    return NULL;
}

/**
 * Returns new matrix, multiplying the two matrices together
 */
uint32_t *matrix_mul(const uint32_t *matrix_a, const uint32_t *matrix_b) {
    int width   = g_width; //local
    int threads = g_nthreads;

    struct mul_struct *res_mat = malloc(sizeof(struct mul_struct)*threads);

    uint32_t *result = new_matrix();

    for (size_t i = 0; i < threads; ++i) {
        res_mat[i] = (struct mul_struct){
            .mat_a   = matrix_a,
            .mat_b   = matrix_b,
            .result  = result, //operates on what result points to
            .index   = i,
            .chunk   = width/threads, //do this at the start of this function
            .width   = width,
            .threads = threads,
        };
    }

    pthread_t t_id[threads]; //creating array to store thread IDs

    for (size_t i = 0; i < threads; ++i) {
        // create p_threads
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

    if (exponent == 0) {
        uint32_t *result = identity_matrix();
        return result;
    }

    uint32_t *result = new_matrix();

    memcpy(result, matrix, g_elements * sizeof(*result));

    // return if A^1
    if (exponent == 1) {
        return result;
    }

    int floor_exp = exponent;

    if (floor_exp % 2 != 0) {
        floor_exp = (floor_exp-1) >> 1;
    } else {
        floor_exp = floor_exp >> 1;
    }

    for (int i = 0; i < floor_exp-1; i++) { //-1 because A^1 is itself
        result = matrix_mul(result, matrix);
    }

    if (exponent % 2 == 0) {
        return result = matrix_mul(result, result);
    } else {
        result = matrix_mul(result, result);
        return matrix_mul(result, matrix);
    }


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

    temp = 0;

    for (int i = 0; i < g_height; i++) {
        temp += matrix[i*g_width+i];  //row we're at + i places in
    }

    return temp;
}

/**
 * Worker function for finding min.
 */
void *min_worker(void *cmp_struct){

    struct cmp_struct cmp_mat = *(struct cmp_struct *)cmp_struct; //casting
    //initial min
    cmp_mat.result[cmp_mat.index] = cmp_mat.matrix[0];

    for (int i = cmp_mat.array_start; i < cmp_mat.array_limit; i++) {
        if (cmp_mat.matrix[i] < cmp_mat.result[cmp_mat.index]) {
            cmp_mat.result[cmp_mat.index] = cmp_mat.matrix[i];
        }
    }
    return NULL;

}

/**
 * Returns the largest value in the matrix
 */
uint32_t get_minimum(const uint32_t *matrix) {
    struct cmp_struct *cmp_mat = malloc(sizeof(struct cmp_struct)*g_nthreads);
    uint32_t          *result  = new_matrix();

    int chunk = g_elements/g_nthreads;

    for (size_t i = 0; i < g_nthreads; ++i) {
        cmp_mat[i] = (struct cmp_struct){
            .matrix      = matrix,
            .result      = result,
            .index       = i,
            .array_start = i * chunk,
            .array_limit = (i * chunk)+chunk,
        };
    }

    pthread_t t_id[g_nthreads]; //creating array to store thread IDs

    for (size_t i = 0; i < g_nthreads; ++i) {

        pthread_create(t_id + i, NULL, min_worker, cmp_mat + i);
    }

    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_join(t_id[i], NULL);
    }

    int min = cmp_mat->result[0];

    for (size_t i = 0; i < g_nthreads; ++i) {
        if (cmp_mat->result[i] < min) {
            min = cmp_mat->result[i];
        }
    }

    free(cmp_mat);
    return min;

}

/**
 * Worker function for finding max.
 */
void *max_worker(void *cmp_struct){

    struct cmp_struct cmp_mat = *(struct cmp_struct *)cmp_struct;//casting

    for (int i = cmp_mat.array_start; i < cmp_mat.array_limit; i++) {
        if (cmp_mat.matrix[i] > cmp_mat.result[cmp_mat.index]) {
            cmp_mat.result[cmp_mat.index] = cmp_mat.matrix[i];
        }
    }
    return NULL;

}

/**
 * Returns the largest value in the matrix
 */
uint32_t get_maximum(const uint32_t *matrix) {
    struct cmp_struct *cmp_mat = malloc(sizeof(struct cmp_struct)*g_nthreads);
    uint32_t          *result  = new_matrix();

    int max   = 0;
    int chunk = g_elements/g_nthreads;

    for (size_t i = 0; i < g_nthreads; ++i) {
        cmp_mat[i] = (struct cmp_struct){
            .matrix      = matrix,
            .result      = result,
            .index       = i,
            .array_start = i*chunk,
            .array_limit = (i*chunk)+chunk,
        };
    }

    pthread_t t_id[g_nthreads]; //creating array to store thread IDs

    for (size_t i = 0; i < g_nthreads; ++i) {

        pthread_create(t_id + i, NULL, max_worker, cmp_mat + i);
    }

    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_join(t_id[i], NULL);
    }

    for (size_t i = 0; i < g_nthreads; ++i) {
        if (cmp_mat->result[i] > max) {
            max = cmp_mat->result[i];
        }
    }

    free(cmp_mat);
    return max;

}

/**
 * Returns the frequency of the value in the matrix
 */
inline uint32_t get_frequency(const uint32_t *matrix, uint32_t value) {

    temp = 0;

    for (int i = 0; i < g_elements; i++) {
        if (matrix[i] == value) {
            temp++;
        }
    }

    return temp;
}
