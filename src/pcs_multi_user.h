#include <gmp.h>
#include <omp.h>
#include <inttypes.h>
#include <time.h>
#include <mpi.h>
#include "pcs.h"

//#define __NB_USERS__ ((1<<16)-4)
#define __NB_USERS__ 16


//#define SEED_ (time(NULL))
#define SEED_ (time(NULL)^0xcafebabe)
//#define SEED_ 0x9469bfc9
//#define SEED_ 0x94685a71


long int SEED;

void set_seed();

void pcs_mu_init(point_t P_init,
                   point_t Q_init[__NB_USERS__],
                   elliptic_curve_t E_init,
                   mpz_t n_init,
		               mpz_t *A_init,
                   uint8_t nb_bits_init,
                   uint8_t trailling_bits_init,
                   int type_struct,
                   int nb_threads,
                   uint8_t level);

long long int pcs_mu_run(mpz_t x_res,
                           int nb_threads,
                           int nb_collisions);

long long int pcs_mu_run_order(mpz_t x_res[__NB_USERS__],
			       int nb_threads,
			       unsigned long long int times[__NB_USERS__],
			       unsigned long int pts_per_users[__NB_USERS__]);



void pcs_mu_clear();

int share_mpz_array(mpz_t *a, int size, int world_rank, int init);
int share_mpz_var(mpz_t a, int world_rank, int init);
int send_mpz_var(mpz_t a, int dest, int tag);
int recv_mpz_var(mpz_t a, int src, int tag, int init);



long long int pcs_mu_run_order_server(mpz_t x_res[__NB_USERS__], int nb_threads, unsigned long long int times[__NB_USERS__],unsigned long int pts_per_users[__NB_USERS__]);
long long int pcs_mu_run_order_client(int nb_threads, int world_rank);



void pcs_mu_init_server(
									 point_t  P_init,
									 point_t Q_init[__NB_USERS__],
									 elliptic_curve_t E_init,
									 mpz_t n_init,
									 mpz_t *A_init,
                   uint8_t nb_bits_init,
                   uint8_t trailling_bits_init,
                   int type_struct,
                   int nb_threads,
                   uint8_t level,
                   int world_size);



void pcs_mu_init_client(
                    point_t P_init,
                    point_t Q_init[__NB_USERS__],
                    elliptic_curve_t E_init,
                    mpz_t n_init,
 		                mpz_t *A_init,
                    uint8_t nb_bits_init,
                    uint8_t trailling_bits_init,
                    int nb_threads,
                    uint8_t level);

void pcs_mu_clear_server();
void pcs_mu_clear_client();
void dbg_init_xtrue(mpz_t xtrue_init[__NB_USERS__]);
