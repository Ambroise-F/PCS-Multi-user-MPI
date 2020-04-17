#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <gmp.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>
#include "pcs_elliptic_curve_operations.h"
#include "pcs_pollard_rho.h"
#include "pcs_storage.h"
#include "pcs.h"
#include "pcs_multi_user.h"

#define FF fflush(stdout)
#define qp(a) printf("%d",a);fflush(stdout)
#define verbose 0
#define ATT 0
//#define SEED (time(null))
//#define SEED 0xE5CA1ADE


#define TAG_START 1

#define TAG_END 2

#define TAG_THREAD_OFFSET 3 // start of the tag range for threads


elliptic_curve_t E;
point_t P;
point_t Q[__NB_USERS__]; // One Q per user
mpz_t n;
mpz_t *A;
mpz_t *B;
mpz_t X_res[__NB_USERS__]; // One x per user
point_t M[__NB_ENSEMBLES__]; // precomputed a*P values
uint8_t trailling_bits;
uint8_t nb_bits;
int world_size;

mpz_t xtrue[__NB_USERS__];




void set_seed()
{
	SEED = SEED_;
	printf("seed : %lx\n",SEED);
}

/** Determines whether a point is a distinguished one.
*
*  @param[in]	R				A point on an elliptic curve.
*  @param[in]	trailling_bits	Number of trailling zero bits in a ditinguished point.
*  @param[out]	q				The x-coordinate, without the trailling zeros.
*  @return 	1 if the point is distinguished, 0 otherwise.
*/
int is_distinguished_mu(point_t R, int trailling_bits, mpz_t *q)
{
	int res;
	mpz_t r;
	mpz_inits(r, NULL);
	mpz_tdiv_qr_ui(*q, r, R.x, (unsigned long int)pow(2, trailling_bits));
	res=(mpz_sgn(r) == 0);
	mpz_clears(r, NULL);
	return (res);
}


/** Checks if the linear combination aP+bQ is equal to R or its inverse.
*
*  @param[in]	R	A point on an elliptic curve.
*  @param[in]	a	a coefficient.
*  @param[in]	b	b coefficient.
*  @param[in]  user    user id
*  @return 	1 if aP+bQ[user] = R, 0 if aP+bQ[user] = -R.
*/
int same_point_mu(point_t R, mpz_t a, mpz_t b, uint16_t user)
{
	int res;
	point_t S1, S2, S;
	mpz_inits(S1.x, S1.y, S1.z, S2.x, S2.y, S2.z, S.x, S.y, S.z, NULL);
	double_and_add(&S1, P, a, E); // S1 = aP
	double_and_add(&S2, Q[user], b, E); // S2 = bQ[user]
	add(&S, S1, S2, E); // S2 = aP + bQ
	res=(mpz_cmp(R.y, S.y) == 0); //
	mpz_clears(S1.x, S1.y, S1.z, S2.x, S2.y, S2.z, S.x, S.y, S.z, NULL);
	return res;
}

/** Computes R = a * P  on E.
*
*  @param[out]	R	Resulting point.
*  @param[in]	a	a coefficient.
*/
void lin_comb_mu(point_t * R, mpz_t a)
{
	double_and_add(R, P, a, E);
}

/** Checks if there is a collision.
*
*/
//int is_collision(mpz_t x, mpz_t a1, mpz_t a2, int trailling_bits)
int is_collision_mu(mpz_t x, mpz_t b1, uint16_t userid1, mpz_t b2, uint16_t userid2, int trailling_bits)
{
	uint8_t r;
	mpz_t xDist_;
	int retval = 0;
	mpz_t a1, a2;
	point_t R;
	// mpz_t xR, xI; // debug


	point_init(&R);
	mpz_inits(a1, a2, xDist_, NULL);

	mpz_set_ui(a2, 0); // a1 = a2 = 0
	mpz_set_ui(a1, 0);


	//recompute first a,b pair

	double_and_add(&R, Q[userid1], b1, E); // R = b1 * Qi


	while(!is_distinguished_mu(R, trailling_bits, &xDist_))
	{
		r = hash(R.y);
		compute_a(a1, A[r], n); // a1 = a1 + A[r] % n   <=> ajout de A[r] au coeff a1
		f(R, M[r], &R, E);      // R  = R  + M[r]       <=> calcul du point R = R + M[r] = R + A[r]*P
	}

	//recompute second a,b pair

	double_and_add(&R, Q[userid2], b2, E);

	while(!is_distinguished_mu(R, trailling_bits, &xDist_))
	{
		r = hash(R.y);
		compute_a(a2, A[r], n);
		f(R, M[r], &R, E);
	}

	//gmp_printf("(a1,b1,a2,b2,n) = (%Zd,%Zd,%Zd,%Zd,%Zd)\n",a1,b1,a2,b2,n);
	/*
	if (userid1==2771)
	  {
	    gmp_printf("a1 = %-10Zd, b1 = %-10Zd, a2 = %-10Zd, b2 = %-10Zd, x1 = %-10Zd, x2 = %-10Zd\n",a1,b1,a2,b2,x_true1,x_true2);
	    gmp_printf("n = %Zd\n",n);

	    mpz_inits(xR, xI, NULL);

	    mpz_mul(xR,b2,x_true2);  // b2*x2
	    mpz_mmod(xR,xR,n);
	    gmp_printf("xR = %Zd\n",xR);
	    mpz_add(xR,xR,a2);  // +a2
	    mpz_mmod(xR,xR,n);
	    gmp_printf("xR = %Zd\n",xR);
	    mpz_sub(xR,xR,a1);  // -a1
	    mpz_mmod(xR,xR,n);
	    gmp_printf("xR = %Zd\n",xR);
	    mpz_invert(xI,b1,n);// b1^-1 mod n
	    gmp_printf("xI = %Zd\n",xI);
	    mpz_mul(xR,xR,xI);  // (b2*x2+a2-a1) * (b1^-1)
	    mpz_mmod(xR,xR,n);
	    gmp_printf("x = %Zd\n",xR);


	    mpz_clears(xR,xI, NULL);
	  }

	// end debug
	*/
	/*
	printf("%d-%d\n",userid1,userid2);
	if (userid1==userid2)
	{
		gmp_printf("(a1,b1,a2,b2,n,x) = (%Zd,%Zd,%Zd,%Zd,%Zd,",a1,b1,a2,b2,n);
	}
	else
	{
		gmp_printf("(a1,b1,a2,b2,n,x2,x) = (%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,",a1,b1,a2,b2,n,X_res[userid2]);
	}
	*/

	if(userid1==userid2 && (mpz_cmp(b1, b2) != 0)) //two different pairs with the same Q, so collision
	{
		if(!same_point_mu(R, a1, b1,userid1)) //it's the inverse point // to be modified
		{
			mpz_neg(a2, a2);
			mpz_mmod(a2, a2, n);
			mpz_neg(b2, b2);
			mpz_mmod(b2, b2, n);
		}
		compute_x(x, a1, a2, b1, b2, n);
		retval = 1;
	}
	else if(userid1!=userid2) //two different Qs, so collision - userid2 has to be known
	{
		if(!same_point_mu(R, a1, b1,userid1)) //it's the inverse point // to be modified
		{
			mpz_neg(a2, a2);
			mpz_mmod(a2, a2, n);
			mpz_neg(b2, b2);
			mpz_mmod(b2, b2, n);
		}
		compute_x_2users(x, a1, a2, b1, b2, X_res[userid2], n);
		retval = 1;
	}

	if (mpz_cmp(xtrue[userid1],x))
	{
		printf("erreur : \n");
		gmp_printf("(a1,b1,a2,b2,n,x1,x2,xT) = (%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd)\n",a1,b1,a2,b2,n,x,X_res[userid2],xtrue[userid1]);
	}




	point_clear(&R);
	mpz_clears(a1, a2, xDist_, NULL);


	return retval;
}


char * pack_(size_t * size_v,int thread_num, uint16_t userid1, mpz_t b, mpz_t xDist, int nb_bits, int trailing_bits)
{
	size_t size_vect,size_b,size_x;
	char *vect;
	int i;

	*size_v = sizeof(int);
	*size_v+= sizeof(uint16_t);

	size_b = (int)((nb_bits-1)/8)+1;
	size_x = (int)((nb_bits-1-trailing_bits)/8)+1;

	*size_v+= size_b + size_x;

	vect = malloc(*size_v);

	mpz_export(vect+sizeof(int)+sizeof(uint16_t),NULL,1,1,-1,0,b);
	mpz_export(vect+sizeof(int)+sizeof(uint16_t)+size_b,NULL,1,1,-1,0,xDist);

	for (i=0;i<sizeof(int);i++)
	{
		vect[i] = (thread_num>>((sizeof(int)-1-i)*8))&0xff; // 0-3 1-4
	}
	for (i=0;i<sizeof(uint16_t);i++)
	{
		vect[i+sizeof(int)] = (userid1>>((sizeof(uint16_t)-1-i)*8))&0xff;
	}
	return vect;
}

char * pack(size_t * size_v,int thread_num, uint16_t userid1, mpz_t b, mpz_t xDist, int nb_bits, int trailing_bits)
{
	/*
	 [----]        [--]     [-][-]     [.?.]    [.?.]
	  int         uint16    2*char     mpz_t    mpz_t
	thread_num    userid1  size_b/x      b        x
	*/
	size_t size_vect,size_b,size_x;
	size_t bsize_b,bsize_x;
	char *vect;
	int i;

	bsize_b = mpz_sizeinbase(b,2);
	bsize_x = mpz_sizeinbase(xDist,2);

	size_b = (unsigned char)((bsize_b-1)/8)+1;
	size_x = (unsigned char)((bsize_x-1)/8)+1;

	*size_v = 2*sizeof(unsigned char);
	*size_v+= sizeof(int);
	*size_v+= sizeof(uint16_t);
	*size_v+= size_b + size_x;

	vect = malloc(*size_v);

	mpz_export(vect+2*sizeof(unsigned char)+sizeof(int)+sizeof(uint16_t),NULL,1,1,-1,0,b);
	mpz_export(vect+2*sizeof(unsigned char)+sizeof(int)+sizeof(uint16_t)+size_b,NULL,1,1,-1,0,xDist);

	for (i=0;i<sizeof(int);i++)
	{
		vect[i] = (thread_num>>((sizeof(int)-1-i)*8))&0xff; // 0-3 1-4
	}
	for (i=0;i<sizeof(uint16_t);i++)
	{
		vect[i+sizeof(int)] = (userid1>>((sizeof(uint16_t)-1-i)*8))&0xff;
	}
	for (i=0;i<sizeof(unsigned char);i++)
	{
		vect[i+sizeof(int)+sizeof(uint16_t)] = 0xff&(unsigned char)size_b;
	}
	for (i=0;i<sizeof(unsigned char);i++)
	{
		vect[i+sizeof(int)+sizeof(uint16_t)+sizeof(unsigned char)] = 0xff&(unsigned char)size_x;
	}
	return vect;
}




int unpack(char * vect ,int *thread_num, uint16_t *userid1, mpz_t b, mpz_t xDist, int nb_bits, int trailing_bits)
{
	/*
	 [----]        [--]     [-][-]     [.?.]    [.?.]
		int         uint16    2*char     mpz_t    mpz_t
	thread_num    userid1  size_b/x      b        x
	*/
	size_t size_b,size_x;
	int i;
	i=0;
//	size_b = (int)((nb_bits-1)/8)+1;
//	size_x = (int)((nb_bits-1-trailing_bits)/8)+1;
	*thread_num = 0;
	for (i;i<sizeof(int);i++)
	{
		*thread_num <<= 8;
		*thread_num += (vect[i]&0xff);
	}
	*userid1 = 0;
	for (i;i<sizeof(int)+sizeof(uint16_t);i++)
	{
		*userid1 <<= 8;
		*userid1 += vect[i]&0xff;
	}
	for (i;i<sizeof(int)+sizeof(uint16_t)+sizeof(unsigned char);i++)
	{
		size_b = (size_t)vect[i]&0xff;
	}
	for (i;i<sizeof(int)+sizeof(uint16_t)+2*sizeof(unsigned char);i++)
	{
		size_x = (size_t)vect[i]&0xff;
	}
	mpz_import(b,size_b,1,1,-1,0,vect+i);
	i+=size_b;
	mpz_import(xDist,size_x,1,1,-1,0,vect+i);
	return i+size_x;
}

int generate_random_b(mpz_t b, int nb_bits, gmp_randstate_t r_state)
{
	mpz_t p2,nbits,two;
	mpz_inits(p2,nbits,two,NULL);
	mpz_set_ui(nbits, nb_bits-2);
	mpz_set_ui(two,2);
	mpz_powm(p2,two,nbits,n);

	mpz_urandomb(b, r_state, nb_bits-2);
	mpz_add(b, b, p2); // b = b + 2**nb_bits-2
	while(mpz_cmp(b,n)>0)
	{
		mpz_urandomb(b, r_state, nb_bits-2);
		mpz_add(b, b, p2); // b = 2 * b
	} // fonction pour random b ??? -> pile nb_bits-1 bits mais < n
	mpz_clears(p2,nbits,two,NULL);
}

/** Initialize all variables needed to do a PCS algorithm.
*
*/
void pcs_mu_init(point_t  P_init,
	point_t Q_init[__NB_USERS__],
	elliptic_curve_t E_init,
	mpz_t n_init,
	mpz_t *A_init,
	uint8_t nb_bits_init,
	uint8_t trailling_bits_init,
	int type_struct,
	int nb_threads,
	uint8_t level)
{
		uint8_t i; //  __NB_ENSEMBLES__
		int j; // __NB_USERS__
		uint16_t user;


		point_init(&P);
		//point_init(&Q);
		curve_init(&E);
		mpz_init(n);


		mpz_set(P.x, P_init.x);
		mpz_set(P.y, P_init.y);
		mpz_set(P.z, P_init.z);

		//mpz_set(Q.x, Q_init.x);
		//mpz_set(Q.y, Q_init.y);
		//mpz_set(Q.z, Q_init.z);

		mpz_set(E.A, E_init.A);
		mpz_set(E.B, E_init.B);
		mpz_set(E.p, E_init.p);

		mpz_set(n, n_init);

		A = A_init;



		for(j=0; j<__NB_USERS__; j++) // Q init
		{
			mpz_set(Q[j].x,Q_init[j].x);
			mpz_set(Q[j].y,Q_init[j].y);
			mpz_set(Q[j].z,Q_init[j].z);
			user = (uint16_t) j;

		}
		for(i=0; i<__NB_ENSEMBLES__; i++) // has to be after Q inits at index j since it is needed for lin_comb_mu
		{
			mpz_inits(M[i].x,M[i].y,M[i].z,NULL);
			lin_comb_mu(&M[i],A[i]);
		}

		trailling_bits = trailling_bits_init;
		nb_bits = nb_bits_init;
		struct_init_mu(type_struct, n, trailling_bits, nb_bits, nb_threads, level); // TODO : mu adaptation : done?
}

/** Initialize all variables needed to do a PCS algorithm.
*
*/
void pcs_mu_init_client(point_t  P_init,
	point_t Q_init[__NB_USERS__],
	elliptic_curve_t E_init,
	mpz_t n_init,
	mpz_t *A_init,
	uint8_t nb_bits_init,
	uint8_t trailling_bits_init,
	int nb_threads,
	uint8_t level)
{
		uint8_t i; //  __NB_ENSEMBLES__
		int j; // __NB_USERS__
		uint16_t user;


		point_init(&P);
		//point_init(&Q);
		curve_init(&E);
		mpz_init(n);


		mpz_set(P.x, P_init.x);
		mpz_set(P.y, P_init.y);
		mpz_set(P.z, P_init.z);

		//mpz_set(Q.x, Q_init.x);
		//mpz_set(Q.y, Q_init.y);
		//mpz_set(Q.z, Q_init.z);

		mpz_set(E.A, E_init.A);
		mpz_set(E.B, E_init.B);
		mpz_set(E.p, E_init.p);

		mpz_set(n, n_init);

		A = A_init;



		for(j=0; j<__NB_USERS__; j++) // Q init
		{
			mpz_set(Q[j].x,Q_init[j].x);
			mpz_set(Q[j].y,Q_init[j].y);
			mpz_set(Q[j].z,Q_init[j].z);
			user = (uint16_t) j;

		}
		for(i=0; i<__NB_ENSEMBLES__; i++) // has to be after Q inits at index j since it is needed for lin_comb_mu
		{
			mpz_inits(M[i].x,M[i].y,M[i].z,NULL);
			lin_comb_mu(&M[i],A[i]);
		}

		trailling_bits = trailling_bits_init;
		nb_bits = nb_bits_init;
		printf("Init client\n" );
}

/** Initialize all variables needed to do a PCS algorithm.
*
*/
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
                   int world_size_init)
{
		uint8_t i; //  __NB_ENSEMBLES__
		int j; // __NB_USERS__
		uint16_t user;
		world_size = world_size_init;


		point_init(&P);
		//point_init(&Q);
		curve_init(&E);
		mpz_init(n);


		mpz_set(P.x, P_init.x);
		mpz_set(P.y, P_init.y);
		mpz_set(P.z, P_init.z);

		//mpz_set(Q.x, Q_init.x);
		//mpz_set(Q.y, Q_init.y);
		//mpz_set(Q.z, Q_init.z);

		mpz_set(E.A, E_init.A);
		mpz_set(E.B, E_init.B);
		mpz_set(E.p, E_init.p);

		mpz_set(n, n_init);

		A = A_init;

		for(j=0; j<__NB_USERS__; j++) // Q init
		{
			mpz_set(Q[j].x,Q_init[j].x);
			mpz_set(Q[j].y,Q_init[j].y);
			mpz_set(Q[j].z,Q_init[j].z);
			user = (uint16_t) j;

		}
		for(i=0; i<__NB_ENSEMBLES__; i++) // has to be after Q inits at index j since it is needed for lin_comb_mu
		{
			mpz_inits(M[i].x,M[i].y,M[i].z,NULL);
			lin_comb_mu(&M[i],A[i]);
		}

		trailling_bits = trailling_bits_init;
		nb_bits = nb_bits_init;


		printf("nb_threads = %d\n",nb_threads);
		printf("world_size = %d\n",world_size);
		struct_init_mu(type_struct, n, trailling_bits, nb_bits, nb_threads*(world_size-1), level); // TODO : mu adaptation : done?

		printf("Init server\n" );
}



/** Run the PCS algorithm.
*
*/
long long int pcs_mu_run_order(mpz_t x_res[__NB_USERS__], int nb_threads, unsigned long long int times[__NB_USERS__],unsigned long int pts_per_users[__NB_USERS__])
{
	point_t R;
	mpz_t b, b2;
	mpz_t x, xDist;
	uint8_t r;
	int trail_length;
	int col;
	int trail_length_max = pow(2, trailling_bits) * 20; // 20 * 1<<trailling_bits
	int collision_count = 0;
	unsigned long long int time1,time2;
	// int userid1;
	// int* userid2;
	uint16_t userid1,userid2;
	char xDist_str[50];
	struct timeval tv1,tv2;
	unsigned long int nb_pts;
	//unsigned long long int times[__NB_USERS__];
	//nb_threads = 1;
	//nb_threads = omp_get_max_threads();
	//nb_threads = __NB_USERS__; // for testing purposes : userid1 = thread number

	for (userid1=0; userid1<__NB_USERS__; userid1++)
	{
		pts_per_users[userid1] = 0;
		gettimeofday(&tv1,NULL);
		#pragma omp parallel private(nb_pts,userid2, R, b, b2, x, r, xDist, xDist_str, trail_length,col) shared(collision_count, X_res, trail_length_max,pts_per_users) num_threads(nb_threads)
		{
			nb_pts = 0;
			col = 0;
			point_init(&R);
			mpz_inits(x, b, b2, xDist, NULL);


			//Initialize a starting point
			gmp_randstate_t r_state;
			gmp_randinit_default(r_state);
			gmp_randseed_ui(r_state, SEED * (userid1+1) * (omp_get_thread_num() + 1));
			mpz_urandomb(b, r_state, nb_bits); // random b
			double_and_add(&R, Q[userid1], b, E); // R = bQi
			trail_length = 0;
			collision_count = 0;
			while(collision_count < 1)
			{
				if(is_distinguished_mu(R, trailling_bits, &xDist)) // xDist = R.x >> trailling_bits
				{
					//printf("dist point!\n");fflush(stdout);

					userid2 = __NB_USERS__; // debug / useless
					nb_pts++;
					//printf(".");FF;

					if(struct_add_mu(b2, &userid2, b, userid1, xDist, xDist_str)) // ajout de b dans la mémoire, b2 = b d'un autre point avec collision
					{
						//printf("added\n)");fflush(stdout);
						if(is_collision_mu(x, b, userid1, b2, userid2, trailling_bits)) // si b et b2 forment une vraie collision
						{
							//printf("\nThread num %d :\n",omp_get_thread_num());
							//printf("True collision %2hu - %2hu",userid1,userid2);
							if(verbose)
							{
								printf("coll(%hu-%hu);",userid1,userid2);FF;
							}
							//if(userid1!=userid2) printf(" ---- different origin");
							col = 1;
							//printf("\n");

							#pragma omp critical
							{
								collision_count++;
								mpz_init_set(X_res[userid1],x);
							}
						}
					}
					//              else // pt distingué ajouté, pas de collision
					if(col==0)
					{
						mpz_urandomb(b, r_state, nb_bits);
						double_and_add(&R, Q[userid1], b, E); // new start, R = bQi
						trail_length = 0;
					}
				}
				else // R n'est pas un pt dist.
				{


					r=hash(R.y); // y%20
					//printf("f...\n");FF;
					//printf("%d,%d\n",userid1,r);FF;
					//gmp_printf("%Zd\n",M[userid1][r].x);FF;
					f(R, M[r], &R, E); // 1 step (among 20) of the path
					//gmp_printf("%Zd.",M[userid1][r].x);
					//gmp_printf("%Zd-",M[userid1][r].y);FF;
					//printf("f - ok\n");FF;

					trail_length++;
					if(trail_length > trail_length_max)
					{
						mpz_urandomb(b, r_state, nb_bits); // new random start
						double_and_add(&R, Q[userid1], b, E);
						trail_length = 0;
					}

				}

			} // end while
			//printf("thread %d : %lu\n",omp_get_thread_num(),nb_pts);
			#pragma omp critical
			{
				pts_per_users[userid1]+=nb_pts;
			}

			point_clear(&R);
			mpz_clears(b, b2, x, xDist, NULL);
			gmp_randclear(r_state);
		} // end omp parallel

		//printf("user %d : %lu pts\n",userid1,pts_per_users[userid1]);
		gettimeofday(&tv2,NULL);
		// times.append(tv2-tv1)
		time1 = (tv1.tv_sec) * 1000000 + tv1.tv_usec;
		time2 = (tv2.tv_sec) * 1000000 + tv2.tv_usec;
		times[userid1] = time2 - time1;
		if (verbose && !userid1%((int)(__NB_USERS__/100)))
		{
			printf("#");FF;
		}


	} // end for
	if (verbose)
	printf("\n");
	for (userid1=0; userid1<__NB_USERS__; userid1++)
	{
		mpz_init_set(x_res[userid1],X_res[userid1]);
	}
	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

long long int pcs_mu_run_order_server(mpz_t x_res[__NB_USERS__], int nb_threads, unsigned long long int times[__NB_USERS__],unsigned long int pts_per_users[__NB_USERS__])
{
	char * payload;
	MPI_Status status;
	int tag,thread_num,coll,i;
	uint16_t userid1,userid2,current_user;
	int resp;
	mpz_t b, b2;
	mpz_t x, xDist;
  char xDist_str[50];
	struct timeval tv1,tv2;
	unsigned long int nb_pts;
	unsigned long long int time1,time2;
	unsigned char end;
	MPI_Request req;
	int reqflag;

  mpz_inits(b,b2,xDist,NULL);
	current_user = 0;
	nb_pts = 0;
	end = 0;
	int payload_size;
	payload_size = 2*sizeof(unsigned char)+sizeof(int)+sizeof(uint16_t)+(int)((nb_bits-1)/8)+1+(int)((nb_bits-1-trailling_bits)/8)+1;
	payload = malloc(payload_size); // size of thread_num (int), userid (uint_16), b and xDist (mpz_t)

	gettimeofday(&tv1,NULL);

	while(!end)
	{
		//MPI_Recv(payload,payload_size,MPI_CHAR, MPI_ANY_SOURCE, TAG_START, MPI_COMM_WORLD,&status);
		MPI_Irecv(payload,payload_size,MPI_CHAR, MPI_ANY_SOURCE, TAG_START, MPI_COMM_WORLD,&req);
		reqflag = 0;
		while(!reqflag)
		{
			MPI_Test(&req,&reqflag,&status);
			if (ATT) printf("ATTENTE IRECV TO ANY THREAD                               \r");FF;
		}


		unpack(payload,&thread_num,&userid1,b,xDist,nb_bits,trailling_bits);
	//printf("%d",thread_num);
		tag = thread_num + TAG_THREAD_OFFSET;
			//gmp_printf("payload<: %-10x,%-10x,%-10Zx,%-10Zd\n",thread_num,userid1,b,xDist);

		resp = current_user;
		if (userid1 == current_user)
		{
			coll = struct_add_mu(b2,&userid2,b,userid1,xDist,xDist_str);
		  //printf("coll : %d\n",coll );
			nb_pts++;
			if (coll)
			{
			//printf("coll %d-%d\n",userid1,userid2);
			//gmp_printf("(b,b2,n,x2) = (%Zd,%Zd,%Zd,%Zd)\n",b,b2,n,x_res[userid2]);
				if(is_collision_mu(x, b, userid1, b2, userid2, trailling_bits))
				{
					//gmp_printf("b = %Zd, b2 = %Zd, lp = %Zd\n",b,b2,n);
					mpz_init_set(x_res[userid1],x);
					gettimeofday(&tv2,NULL);
					mpz_init_set(X_res[userid1],x);
					pts_per_users[userid1] = nb_pts;
					time1 = (tv1.tv_sec) * 1000000 + tv1.tv_usec;
					time2 = (tv2.tv_sec) * 1000000 + tv2.tv_usec;
					times[userid1] = time2 - time1;
					gettimeofday(&tv1,NULL);

					nb_pts = 0;
					printf("found user %d\r",current_user);FF;
					current_user++;
					resp = current_user;
				}
			}
		}
		if (current_user==__NB_USERS__)
		{
			end = 1;
			resp = -1;

		}
		//MPI_Send(&resp,1,MPI_INT,status.MPI_SOURCE, tag, MPI_COMM_WORLD);

		MPI_Isend(&resp,1,MPI_INT,status.MPI_SOURCE,tag, MPI_COMM_WORLD, &req);
		reqflag = 0;
		while(!reqflag)
		{
			MPI_Test(&req,&reqflag,NULL);
			if (ATT) printf("ATTENTE ISEND TO THREAD %d                \r",status.MPI_SOURCE);FF;
		}

	}
	int end_int;
	end_int = 1;
	for (i=1;i<world_size;i++)
	{
		MPI_Send(&end_int,1,MPI_INT,i,TAG_END,MPI_COMM_WORLD);
	}
	//printf("\n");
}

/** Run the PCS algorithm.
*
*/
long long int pcs_mu_run_order_client(int nb_threads, int world_rank)
{
	point_t R;
	mpz_t b, xDist;
	uint8_t r;
	int trail_length;
	int trail_length_max = pow(2, trailling_bits) * 20; // 20 * 1<<trailling_bits
	// int userid1;
	// int* userid2;
	uint16_t userid1;
	char xDist_str[50];
	uint8_t end;
	int thread_num;
	int tag;
	MPI_Request req;
	int flagreq;
	//unsigned long long int times[__NB_USERS__];
	//nb_threads = 1;
	//nb_threads = omp_get_max_threads();
	//nb_threads = __NB_USERS__; // for testing purposes : userid1 = thread number

	//for (userid1=0; userid1<__NB_USERS__; userid1++)

	end = 0;


	#pragma omp parallel private(userid1,R, b, r, xDist, xDist_str, trail_length,thread_num,tag,req,flagreq) shared(end, trail_length_max) num_threads(nb_threads+1)
	{
		thread_num = omp_get_thread_num();
		if (!thread_num)
		{
			int end_resp,end_flag;
			MPI_Request end_req;
			MPI_Status end_status;
			MPI_Irecv(&end_resp, 1, MPI_INT, 0, TAG_END, MPI_COMM_WORLD,&end_req);
			end_flag = 0;
			while(!end_flag && !end)
			{
				MPI_Test(&end_req,&end_flag,&end_status);
				//printf("ATTENTE THREAD STOP                                   \r");FF;
			}
			end=1;
		}
		else
		{
			userid1=0;
			int resp;
			point_init(&R);
			mpz_inits(b, xDist, NULL);

			//Initialize a starting point
			gmp_randstate_t r_state;
			gmp_randinit_default(r_state);
			//printf("seed = %d\n",(SEED<<16)+(omp_get_thread_num()<<8)+(world_rank));
			gmp_randseed_ui(r_state, (SEED<<16)+(omp_get_thread_num()<<8)+(world_rank));
			//printf("stseed = %d\n", (SEED<<16)+(omp_get_thread_num()<<8)+(world_rank));
			//mpz_urandomb(b, r_state, nb_bits); // random b
			//mpz_mod(b,b,n);
			generate_random_b(b,nb_bits,r_state);
			double_and_add(&R, Q[userid1], b, E); // R = bQi
			trail_length = 0;
			while(!end)
			{
				if(is_distinguished_mu(R, trailling_bits, &xDist)) // xDist = R.x >> trailling_bits
				{
					/*
					// build package with thread_num, userid1, b, xDist and send it at once with MPI_Isend
					// receive the reply with MPI_Recv
					*/
					char * payload;
					size_t size_vect;


						//gmp_printf("payload : %-10x,%-10x,%-10Zx,%-10Zd\n",thread_num,userid1,b,xDist);

					payload = pack(&size_vect,thread_num, userid1, b, xDist, nb_bits, trailling_bits);
					/*
					//printf("size_vect = %ld (should be 4+2+4+3 = 13)\n",size_vect);
					//mpz_inits(bb,xx,NULL);
					//unpack(payload,&thr,&user,bb,xx,nb_bits,trailling_bits);
					//gmp_printf("payload : %-10x,%-10x,%-10Zx,%-10Zx\n",thr,user,bb,xx);
					*/
					MPI_Isend(payload, size_vect,   MPI_CHAR,      0,   TAG_START,   MPI_COMM_WORLD,&req);
					flagreq = 0;
					while(!flagreq && !end) // while not sent and not finished
					{
						MPI_Test(&req,&flagreq,NULL);
						//printf("ATTENTE ISEND THREAD %d                      \r",thread_num);FF;
					}
					if(end)
					{
						break;
					}
					free(payload);
					tag = thread_num + TAG_THREAD_OFFSET;
					MPI_Irecv(&resp,1,MPI_INT,0,tag,MPI_COMM_WORLD,&req);
					flagreq = 0;
					while(!flagreq && !end) // while not received and not finished
					{
						MPI_Test(&req,&flagreq,NULL);
						//printf("ATTENTE IRECV THREAD %d                    \r",thread_num);FF;
					}
					if(end)
					{
						break;
					}

					if (resp==-1) // si resp est à -1 c'est terminé, sinon userid1 = resp
					{
						#pragma omp critical
						{
							end = 1;
						}
					}
					else
					{
						userid1 = (uint16_t)resp;
						generate_random_b(b,nb_bits,r_state);
						double_and_add(&R, Q[userid1], b, E); // new start, R = bQi
						trail_length = 0;
					}
				}
				else // R n'est pas un pt dist.
				{
					r=hash(R.y); // y%20
					f(R, M[r], &R, E); // 1 step (among 20) of the path

					trail_length++;
					if(trail_length > trail_length_max)
					{
						//mpz_urandomb(b, r_state, nb_bits); // new random start
						//mpz_mod(b,b,n);
						generate_random_b(b,nb_bits,r_state);

						double_and_add(&R, Q[userid1], b, E);
						trail_length = 0;
					}
				}
	    }

			point_clear(&R);
			mpz_clears(b, xDist, NULL);
			gmp_randclear(r_state);
		}
		//printf("fin thread n°%d\n",thread_num);
	}
// end parallel search

	//printf("user %d : %lu pts\n",userid1,pts_per_users[userid1]);
	return 0;
}

/** Free all variables used in the previous PCS run.
*
*/
void pcs_mu_clear()
{

	uint32_t i;
	point_clear(&P);
	/*
	for (i=0;i<__NB_USERS__;i++)
	{
	point_clear(&Q[i]);
	//mpz_clears(Q[i].x, Q[i].y, Q[i].z, NULL);
}
free(Q);
printf("cleared\n");
printf("clearing M... ");fflush(stdout);
*/
for(i = 0; i < __NB_ENSEMBLES__; i++)
{
	point_clear(&M[i]);
	//mpz_clears(M[i].x, M[i].y, M[i].z, NULL);
}

curve_clear(&E);
mpz_clear(n);
struct_free_mu();
}

/** Free all variables used in the previous PCS run.
*
*/
void pcs_mu_clear_server()
{
	mpz_clear(n);
	struct_free_mu();
}

/** Free all variables used in the previous PCS run.
*
*/
void pcs_mu_clear_client()
{

	uint32_t i;
	point_clear(&P);

	for(i = 0; i < __NB_ENSEMBLES__; i++)
	{
		point_clear(&M[i]);
	}

	curve_clear(&E);
	mpz_clear(n);
}

int share_mpz_var(mpz_t a, int world_rank, int init)
{
	char * a_char;
	int count_int;

	if (!world_rank)
	{
		size_t count;
		a_char = mpz_export(NULL,&count,1,1,1,0,a);
		count_int = (int)count;
	}

	MPI_Bcast(&count_int, 1 , MPI_INT, 0, MPI_COMM_WORLD);
	if(world_rank)
	{
		a_char = (char*)malloc(count_int);
	}

	MPI_Bcast(a_char, count_int, MPI_CHAR, 0, MPI_COMM_WORLD);

	if(world_rank)
	{
		if(!init)
		mpz_init(a);
		mpz_import(a, count_int, 1, 1, 1, 0, a_char);
	}
	return 0;
}

int share_mpz_array(mpz_t *a, int size, int world_rank, int init) // a already malloc'ed
{
	int i,j;
	char * a_char;
	char * a_;
	int count_int;

	if (!world_rank)
	{
		size_t count;
		for (i=0;i<size;i++)
		{
			a_ = mpz_export(NULL,&count,1,1,1,0,a[i]);
			count_int = (int)count;
			if(i==0)
			{
				a_char = malloc(size*count_int);
			}
			for (j=0;j<count_int;j++)
			{
				a_char[i*count_int+j] = a_[j];
			}
			free(a_);
		}
	}

	MPI_Bcast(&count_int, 1 , MPI_INT, 0, MPI_COMM_WORLD);

	if(world_rank)
	{
		a_char = (char*)malloc(count_int*size);
	}

	MPI_Bcast(a_char, count_int*size, MPI_CHAR, 0, MPI_COMM_WORLD);

	if(world_rank)
	{
		a_ = (char*)malloc(count_int);
		for (i=0; i<size; i++)
		{
			for (j=0;j<count_int;j++)
			{
				a_[j] = a_char[i*count_int+j];
			}
			if(!init)
			mpz_init(a[i]);
			mpz_import(a[i],count_int,1,1,1,0,a_);
		}
		free(a_);
	}
	return 0;
}

int send_mpz_var(mpz_t a, int dest, int tag)
{
	char * a_char;
	size_t count;
	int count_int;
	a_char = mpz_export(NULL,&count,1,1,1,0,a);
	count_int = (int) count;
	MPI_Send(&count_int,         1,  MPI_INT, dest, tag, MPI_COMM_WORLD);
	MPI_Send(    a_char, count_int, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
	return 0;
}

int recv_mpz_var(mpz_t a, int src, int tag, int init)
{
	size_t count;
	char* a_char;
	MPI_Recv(&count,          1,  MPI_INT, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	a_char = (char*)malloc((int)count);
	MPI_Recv(a_char, (int)count, MPI_CHAR, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	if (!init)
	mpz_init(a);
	mpz_import(a,(int)count,1,1,1,0,a_char);
	free(a_char);
	return 0;
}

void dbg_init_xtrue(mpz_t xtrue_init[__NB_USERS__])
{
	int i;
	for (i=0;i<__NB_USERS__;i++)
	{
		mpz_init(xtrue[i]);
		mpz_set(xtrue[i],xtrue_init[i]);
	}
}
