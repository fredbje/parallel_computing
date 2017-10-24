/*
 * GENANN - Minimal C Artificial Neural Network
 *
 * Copyright (c) 2015, 2016 Lewis Van Winkle
 *
 * http://CodePlea.com
 *
 * This software is provided 'as-is', without any express or implied
 * warranty. In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgement in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 *
 */

#include "genann_blas.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <sys/time.h>
#include "cblas.h"


#define LOOKUP_SIZE 4096

double walltime() {
  static struct timeval t;
  gettimeofday(&t, NULL);
  return (t.tv_sec + 1e-6 * t.tv_usec);
}

double genann_act_sigmoid(double a) {
    if (a < -45.0) return 0;
    if (a > 45.0) return 1;
    return 1.0 / (1 + exp(-a));
}


double genann_act_sigmoid_cached(double a) {
    /* If you're optimizing for memory usage, just
     * delete this entire function and replace references
     * of genann_act_sigmoid_cached to genann_act_sigmoid
     */
    const double min = -15.0;
    const double max = 15.0;
    static double interval;
    static int initialized = 0;
    static double lookup[LOOKUP_SIZE];

    /* Calculate entire lookup table on first run. */
    if (!initialized) {
        interval = (max - min) / LOOKUP_SIZE;
        int i;
        for (i = 0; i < LOOKUP_SIZE; ++i) {
            lookup[i] = genann_act_sigmoid(min + interval * i);
        }
        /* This is down here to make this thread safe. */
        initialized = 1;
    }

    int i;
    i = (int)((a-min)/interval+0.5);
    if (i <= 0) return lookup[0];
    if (i >= LOOKUP_SIZE) return lookup[LOOKUP_SIZE-1];
    return lookup[i];
}


double genann_act_threshold(double a) {
    return a > 0;
}


double genann_act_linear(double a) {
    return a;
}


genann *genann_init(int inputs, int hidden_layers, int hidden, int outputs) {
    if (hidden_layers < 0) return 0;
    if (inputs < 1) return 0;
    if (outputs < 1) return 0;
    if (hidden_layers > 0 && hidden < 1) return 0;


    const int hidden_weights = hidden_layers ? (inputs+1) * hidden + (hidden_layers-1) * (hidden+1) * hidden : 0;
    const int output_weights = (hidden_layers ? (hidden+1) : (inputs+1)) * outputs;
    const int total_weights = (hidden_weights + output_weights);

    const int total_neurons = (inputs + hidden * hidden_layers + outputs);

    /* Allocate extra size for weights, outputs, and deltas. */
    const int size = sizeof(genann) + sizeof(double) * (total_weights + total_neurons + (total_neurons - inputs));
    
    genann *ret = malloc(size);
    if (!ret) return 0;

    ret->inputs = inputs;
    ret->hidden_layers = hidden_layers;
    ret->hidden = hidden;
    ret->outputs = outputs;

    ret->total_weights = total_weights;
    ret->total_neurons = total_neurons;

    /* Set pointers. */
    ret->weight = (double*)((char*)ret + sizeof(genann));
    ret->output = ret->weight + ret->total_weights;
    ret->delta = ret->output + ret->total_neurons;

    genann_randomize(ret);

    ret->activation_hidden = genann_act_sigmoid_cached;
    ret->activation_output = genann_act_sigmoid_cached;

    printf("sizeof(genann): %d \n", (int)sizeof(genann));
    printf("inputs: %d \n", inputs);
    printf("hidden_layers: %d\n", hidden_layers);
    printf("hidden: %d\n", hidden);
    printf("outputs: %d\n", outputs);
    printf("total_weights: %d\n", total_weights);
    printf("total_neurons: %d\n", total_neurons);
    printf("total_deltas: %d\n", total_neurons - inputs);

    return ret;
}


genann *genann_read(FILE *in) {
    int inputs, hidden_layers, hidden, outputs;
    fscanf(in, "%d %d %d %d", &inputs, &hidden_layers, &hidden, &outputs);

    genann *ann = genann_init(inputs, hidden_layers, hidden, outputs);

    int i;
    for (i = 0; i < ann->total_weights; ++i) {
        fscanf(in, " %le", ann->weight + i);
    }

    return ann;
}


genann *genann_copy(genann const *ann) {
    const int size = sizeof(genann) + sizeof(double) * (ann->total_weights + ann->total_neurons + (ann->total_neurons - ann->inputs));
    genann *ret = malloc(size);
    if (!ret) return 0;

    memcpy(ret, ann, size);

    /* Set pointers. */
    ret->weight = (double*)((char*)ret + sizeof(genann));
    ret->output = ret->weight + ret->total_weights;
    ret->delta = ret->output + ret->total_neurons;

    return ret;
}


void genann_randomize(genann *ann) {
    int i;
    for (i = 0; i < ann->total_weights; ++i) {
        double r = GENANN_RANDOM();
        /* Sets weights from -0.5 to 0.5. */
        ann->weight[i] = r - 0.5;
    }
}


void genann_free(genann *ann) {
    /* The weight, output, and delta pointers go to the same buffer. */
    free(ann);
}


double const *genann_run(genann const *ann, double const *inputs) {
    double const *w = ann->weight;
    double *o = ann->output + ann->inputs;
    double const *i = ann->output;

    /* Copy the inputs to the scratch area, where we also store each neuron's
     * output, for consistency. This way the first layer isn't a special case. */
    memcpy(ann->output, inputs, sizeof(double) * ann->inputs);

    int h, j;

    const genann_actfun act = ann->activation_hidden;
    const genann_actfun acto = ann->activation_output;

    //##########################################################################
    //##########################################################################
    //##########################################################################
    //TODO 1.a The inner double for loops in the following section are doing
    //a matrix vector multiplication(dgemv or sgemv).
    //Optimize it by calling a BLAS library.
    //##########################################################################
    //##########################################################################
    //##########################################################################

    ///////////////////// matrix vector multiplication implementation - not working /////////////////////////
    // double *x = (double*)malloc(sizeof(double) * ((ann->inputs >= ann->hidden ? ann->inputs : ann->hidden) + 1));
    // x[0] = -1.0;
    // for (h = 0; h < ann->hidden_layers; ++h) {
    //     memcpy(x + 1, i, (h == 0 ? ann->inputs : ann->hidden));
    //     cblas_dgemv(CblasRowMajor, CblasNoTrans, ann->hidden, (h == 0 ? ann->inputs : ann->hidden) + 1, 
    //                                 1.0, w, (h == 0 ? ann->inputs : ann->hidden) + 1, x, 1, 0.0, o, 1);
    //     for (j = 0; j < ann->hidden; ++j) {
    //         *o = act(*o);
    //         //printf("h = %d, j = %d, o = %f \n", h, j, *o);
    //         o++;
    //     }
    //     w += ((h == 0 ? ann->inputs : ann->hidden) + 1) * ann->hidden;
    //     i += (h == 0 ? ann->inputs : ann->hidden);
    // }
    // free(x);

    ///////////////////// dot product implementation - working /////////////////////////
    for (h = 0; h < ann->hidden_layers; ++h) {
        for (j = 0; j < ann->hidden; ++j) {
            double sum = *w++ * -1.0;
            sum += cblas_ddot((h == 0 ? ann->inputs : ann->hidden), w, 1, i, 1); 
            w += (h == 0 ? ann->inputs : ann->hidden);
            *o++ = act(sum);
        }
        i += (h == 0 ? ann->inputs : ann->hidden);
    }

    ///////////////////// original implementation /////////////////////////
    // for (h = 0; h < ann->hidden_layers; ++h) {
    //     for (j = 0; j < ann->hidden; ++j) {
    //         double sum = *w++ * -1.0;
    //         for (int k = 0; k < (h == 0 ? ann->inputs : ann->hidden); ++k) {
    //             sum += *w++ * i[k];
    //         }
    //         *o++ = act(sum);
    //     }
    //     i += (h == 0 ? ann->inputs : ann->hidden);
    // }

    

    double const *ret = o;

    //TODO 1.b Double for loop -> blas calls

    ///////////////////// matrix vector product implementation - not working /////////////////////////
    // double *x = (double*)malloc(sizeof(double) * ((ann->hidden_layers ? ann->hidden : ann->inputs) + 1));
    // x[0] = -1.0;
    // memcpy(x + 1, i, (ann->hidden_layers ? ann->hidden : ann->inputs));
    // cblas_dgemv(CblasRowMajor, CblasNoTrans, ann->outputs, (ann->hidden_layers ? ann->hidden : ann->inputs) + 1,
    //                              1.0, w, (ann->hidden_layers ? ann->hidden : ann->inputs) + 1, x, 1, 0.0, o, 1);
    // for (j = 0; j < ann->outputs; ++j) {
    //     *o = acto(*o);
    //     o++;
    // }
    // w += ((ann->hidden_layers ? ann->hidden : ann->inputs) + 1) * ann->outputs;
    // free(x);


    ///////////////////// dot product implementation - working /////////////////////////
    for (j = 0; j < ann->outputs; ++j) {
        double sum = *w++ * -1.0;
        sum += cblas_ddot((ann->hidden_layers ? ann->hidden : ann->inputs), w, 1, i, 1); 
        w += (ann->hidden_layers ? ann->hidden : ann->inputs);
        *o++ = acto(sum);
    }

    ///////////////////// original implementation /////////////////////////
    // for (j = 0; j < ann->outputs; ++j) {
    //     double sum = *w++ * -1.0;
    //     for (int k = 0; k < (ann->hidden_layers ? ann->hidden : ann->inputs); ++k) {
    //         sum += *w++ * i[k];
    //     }
    //     *o++ = acto(sum);
    // }

    /* Sanity check that we used all weights and wrote all outputs. */

    //##########################################################################
    // Commented out because it is not needed when doing the matrix vector
    //multiplication with BLAS
    //##########################################################################
    //assert(w - ann->weight == ann->total_weights);
    //assert(o - ann->output == ann->total_neurons);

    return ret;
}


void genann_train(genann const *ann, double const *inputs, double const *desired_outputs, double learning_rate) {
    /* To begin with, we must run the network forward. */
    genann_run(ann, inputs);

    int h, j, k;

    /* First set the output layer deltas. */
    {
        double const *o = ann->output + ann->inputs + ann->hidden * ann->hidden_layers; /* First output. */
        double *d = ann->delta + ann->hidden * ann->hidden_layers; /* First delta. */
        double const *t = desired_outputs; /* First desired output. */


        /* Set output layer deltas. */
        if (ann->activation_output == genann_act_linear) {
            for (j = 0; j < ann->outputs; ++j) {
                *d++ = *t++ - *o++;
            }
        } else {
            for (j = 0; j < ann->outputs; ++j) {
                *d++ = (*t - *o) * *o * (1.0 - *o);
                ++o; ++t;
            }
        }
    }


    /* Set hidden layer deltas, start on last layer and work backwards. */
    /* Note that loop is skipped in the case of hidden_layers == 0. */
    for (h = ann->hidden_layers - 1; h >= 0; --h) {

        /* Find first output and delta in this layer. */
        double const *o = ann->output + ann->inputs + (h * ann->hidden);
        double *d = ann->delta + (h * ann->hidden);

        /* Find first delta in following layer (which may be hidden or output). */
        double const * const dd = ann->delta + ((h+1) * ann->hidden);

        /* Find first weight in following layer (which may be hidden or output). */
        double const * const ww = ann->weight + ((ann->inputs+1) * ann->hidden) + ((ann->hidden+1) * ann->hidden * (h));

    //##########################################################################
    //##########################################################################
    //##########################################################################
    //TODO 2 The double for loops in the following section are doing
    //a matrix vector multiplication(see TODO 1).
    //Optimize it by calling a BLAS library.
    //##########################################################################
    //##########################################################################
    //##########################################################################

        ///////////////////// matrix vector product implementation - not working /////////////////////////



        ///////////////////// dot product implementation - working /////////////////////////
        for (j = 0; j < ann->hidden; ++j) {
            double delta = cblas_ddot((h == ann->hidden_layers-1 ? ann->outputs : ann->hidden), 
                ww + (j + 1), (ann->hidden + 1), dd, 1);
            *d = *o * (1.0-*o) * delta;
            ++d; ++o;
        }

        ///////////////////// original implementation /////////////////////////
        // for (j = 0; j < ann->hidden; ++j) {
        //     double delta = 0;

        //     for (k = 0; k < (h == ann->hidden_layers-1 ? ann->outputs : ann->hidden); ++k) {
        //         const double forward_delta = dd[k];
        //         const int windex = k * (ann->hidden + 1) + (j + 1);
        //         const double forward_weight = ww[windex];
        //         delta += forward_delta * forward_weight;
        //         // delta += dd[k] * ww[k * (ann->hidden + 1) + (j + 1)];
        //     }

        //     *d = *o * (1.0-*o) * delta;
        //     ++d; ++o;
        // }
    }


    /* Train the outputs. */
    {
        /* Find first output delta. */
        double const *d = ann->delta + ann->hidden * ann->hidden_layers; /* First output delta. */

        /* Find first weight to first output delta. */
        double *w = ann->weight + (ann->hidden_layers
                ? ((ann->inputs+1) * ann->hidden + (ann->hidden+1) * ann->hidden * (ann->hidden_layers-1))
                : (0));

        /* Find first output in previous layer. */
        double const * const i = ann->output + (ann->hidden_layers
                ? (ann->inputs + (ann->hidden) * (ann->hidden_layers-1))
                : 0);

        //##########################################################################
        //##########################################################################
        //##########################################################################
        //TODO 3.a You can optimize the following using the BLAS level 1 SAXPY or DAXPY
        //##########################################################################
        //##########################################################################
        //##########################################################################

        ///////////////////// daxpy implementation - working /////////////////////////
        for (j = 0; j < ann->outputs; ++j) {            
            *w++ += *d * learning_rate * -1.0;
            cblas_daxpy((ann->hidden_layers ? ann->hidden : ann->inputs), *d * learning_rate, i, 1, w, 1);
            w += (ann->hidden_layers ? ann->hidden : ann->inputs);
            ++d;
        }   


        ///////////////////// original implementation /////////////////////////
        /* Set output layer weights. */
        // for (j = 0; j < ann->outputs; ++j) {            
        //     for (k = 0; k < (ann->hidden_layers ? ann->hidden : ann->inputs) + 1; ++k) {
        //         if (k == 0) {
        //             *w++ += *d * learning_rate * -1.0;
        //         } else {
        //             *w++ += *d * learning_rate * i[k-1];
        //         }
        //     }
        //     ++d;
        // }

        // assert(w - ann->weight == ann->total_weights);
    }


    /* Train the hidden layers. */
    for (h = ann->hidden_layers - 1; h >= 0; --h) {

        /* Find first delta in this layer. */
        double const *d = ann->delta + (h * ann->hidden);

        /* Find first input to this layer. */
        double const *i = ann->output + (h
                ? (ann->inputs + ann->hidden * (h-1))
                : 0);

        /* Find first weight to this layer. */
        double *w = ann->weight + (h
                ? ((ann->inputs+1) * ann->hidden + (ann->hidden+1) * (ann->hidden) * (h-1))
                : 0);


        //##########################################################################
        //##########################################################################
        //##########################################################################
        //TODO 3.b You can optimize the following using the BLAS level 1 SAXPY or DAXPY
        //##########################################################################
        //##########################################################################
        //##########################################################################

        ///////////////////// daxpy implementation - not working /////////////////////////
        for (j = 0; j < ann->hidden; ++j) {            
            *w++ += *d * learning_rate * -1.0;
            cblas_daxpy((h == 0 ? ann->inputs : ann->hidden), *d * learning_rate, i, 1, w, 1);
            w += (h == 0 ? ann->inputs : ann->hidden);
            ++d;
        } 

        ///////////////////// original implementation /////////////////////////
        // for (j = 0; j < ann->hidden; ++j) {
        //     for (k = 0; k < (h == 0 ? ann->inputs : ann->hidden) + 1; ++k) {
        //         if (k == 0) {
        //             *w++ += *d * learning_rate * -1.0;
        //         } else {
        //             *w++ += *d * learning_rate * i[k-1];
        //         }
        //     }
        //     ++d;
        // }

    }

}


void genann_write(genann const *ann, FILE *out) {
    fprintf(out, "%d %d %d %d", ann->inputs, ann->hidden_layers, ann->hidden, ann->outputs);

    int i;
    for (i = 0; i < ann->total_weights; ++i) {
        fprintf(out, " %.20e", ann->weight[i]);
    }
}
