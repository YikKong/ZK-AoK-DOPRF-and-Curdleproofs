/****************************************************************************
this hpp implements twisted ElGamal PKE scheme
*****************************************************************************
* @author     This file is part of PGC, developed by Yu Chen
* @paper      https://eprint.iacr.org/2019/319
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/
#include "global.hpp"
#include "hash.hpp"
#include "print.hpp"
#include "routines.hpp"

#include "calculate_dlog.hpp"
//#include "fast_mul.hpp"

const string hashmap_file = "h_point2index.table"; // name of hashmap file

// define the structure of PP
struct ElGamal_PP
{
    size_t IO_THREAD_NUM; // optimized number of threads for faster building hash map 
    size_t DEC_THREAD_NUM; // optimized number of threads for faster decryption: CPU dependent
    EC_POINT* g;
    EC_POINT* h; // two random generators 
};

// define the structure of keypair
struct ElGamal_KP
{
    EC_POINT* pk;  // define pk
    BIGNUM* sk;    // define sk
};

// define the structure of ciphertext
struct ElGamal_CT
{
    EC_POINT* X; // X = pk，r 
    EC_POINT* Y; // Y = g，r+m
};

/* allocate memory for PP */
void ElGamal_PP_new(ElGamal_PP& pp)
{
    pp.g = EC_POINT_new(group);
    pp.h = EC_POINT_new(group);
}

/* free memory of PP */
void ElGamal_PP_free(ElGamal_PP& pp)
{
    EC_POINT_free(pp.g);
    EC_POINT_free(pp.h);
}

void ElGamal_KP_new(ElGamal_KP& keypair)
{
    keypair.pk = EC_POINT_new(group);
    keypair.sk = BN_new();
}

void ElGamal_KP_free(ElGamal_KP& keypair)
{
    EC_POINT_free(keypair.pk);
    BN_free(keypair.sk);
}

void ElGamal_CT_new(ElGamal_CT& CT)
{
    CT.X = EC_POINT_new(group);
    CT.Y = EC_POINT_new(group);
}

void ElGamal_CT_free(ElGamal_CT& CT)
{
    EC_POINT_free(CT.X);
    EC_POINT_free(CT.Y);
}


void ElGamal_PP_print(ElGamal_PP& pp)
{
    ECP_print(pp.g, "pp.g");
    ECP_print(pp.h, "pp.h");
}

void ElGamal_KP_print(ElGamal_KP& keypair)
{
    ECP_print(keypair.pk, "pk");
    BN_print(keypair.sk, "sk");
}

void ElGamal_CT_print(ElGamal_CT& CT)
{
    ECP_print(CT.X, "CT.X");
    ECP_print(CT.Y, "CT.Y");
}

void ElGamal_CT_serialize(ElGamal_CT& CT, ofstream& fout)
{
    ECP_serialize(CT.X, fout);
    ECP_serialize(CT.Y, fout);
}

void ElGamal_CT_deserialize(ElGamal_CT& CT, ifstream& fin)
{
    ECP_deserialize(CT.X, fin);
    ECP_deserialize(CT.Y, fin);
}

/* Setup algorithm */
void ElGamal_Setup(ElGamal_PP& pp,size_t IO_THREAD_NUM, size_t DEC_THREAD_NUM)
{
    pp.IO_THREAD_NUM = IO_THREAD_NUM;
    pp.DEC_THREAD_NUM = DEC_THREAD_NUM;

    EC_POINT_copy(pp.g, generator);
    /* generate pp.h via deterministic manner */
    Hash_ECP_to_ECP(pp.g, pp.h);
#ifdef DEBUG
    cout << "generate the public parameters for ElGamal >>>" << endl;
    ElGamal_PP_print(pp);
#endif
}

/* KeyGen algorithm */
void ElGamal_KeyGen(ElGamal_PP& pp, ElGamal_KP& keypair)
{
    BN_random(keypair.sk); // sk \sample Z_p
    EC_POINT_mul(group, keypair.pk, keypair.sk, NULL, NULL, bn_ctx); // pk = g，sk  

#ifdef DEBUG
    cout << "key generation finished >>>" << endl;
    ElGamal_KP_print(keypair);
#endif
}

/* Encryption algorithm: compute CT = Enc(pk, m; r) */
void ElGamal_Enc(ElGamal_PP& pp,
    EC_POINT*& pk,
    EC_POINT*& m,
    ElGamal_CT& CT)
{
    // generate the random coins 
    BIGNUM* r = BN_new();
    BN_random(r);

    // begin encryption
    // x=r，pk
    EC_POINT_mul(group, CT.X, NULL, pk, r, bn_ctx); 
    // y=g，r+m
    EC_POINT_mul(group, CT.Y, r, NULL, NULL, bn_ctx); 
    EC_POINT_add(group, CT.Y, CT.Y, m, bn_ctx);

    BN_free(r);

#ifdef DEBUG
    cout << " ElGamal encryption finishes >>>" << endl;
    ElGamal_CT_print(CT);
#endif
}

/* Decryption algorithm: compute m = Dec(sk, CT) */
void ElGamal_Dec(ElGamal_PP& pp,
    BIGNUM*& sk,
    ElGamal_CT& CT,
    EC_POINT*& m)
{
    //begin decryption  
    BIGNUM* sk_inverse = BN_new();
    // sk_inverse = sk議剃圷
    BN_mod_inverse(sk_inverse, sk, order, bn_ctx);  // compute the inverse of sk in Z_q^* 

    EC_POINT* M = EC_POINT_new(group);
    // M=pk，r，sk^(-1)=g，r
    EC_POINT_mul(group, M, NULL, CT.X, sk_inverse, bn_ctx); 
    // M=-g，r
    EC_POINT_invert(group, M, bn_ctx);          
    // m+g，r-g，r
    EC_POINT_add(group, M, CT.Y, M, bn_ctx);   
    m = M;
    BN_free(sk_inverse);
    //EC_POINT_free(M);

}



/* homomorphic add */
void ElGamal_HomoAdd(ElGamal_CT& CT_result, ElGamal_CT& CT1, ElGamal_CT& CT2)
{
    EC_POINT_add(group, CT_result.X, CT1.X, CT2.X, bn_ctx);
    EC_POINT_add(group, CT_result.Y, CT1.Y, CT2.Y, bn_ctx);
}

/* homomorphic sub */
void ElGamal_HomoSub(ElGamal_CT& CT_result, ElGamal_CT& CT1, ElGamal_CT& CT2)
{
    EC_POINT_sub(CT_result.X, CT1.X, CT2.X);
    EC_POINT_sub(CT_result.Y, CT1.Y, CT2.Y);
}

/* scalar operation */
void ElGamal_ScalarMul(ElGamal_CT& CT_result, ElGamal_CT& CT, BIGNUM*& k)
{
    EC_POINT_mul(group, CT_result.X, NULL, CT.X, k, bn_ctx);
    EC_POINT_mul(group, CT_result.Y, NULL, CT.Y, k, bn_ctx);
}



















