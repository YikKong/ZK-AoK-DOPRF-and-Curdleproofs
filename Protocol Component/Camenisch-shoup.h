#pragma once
/* ****************************************************************************************************************************************************************** *\
 * @file    csve.h                                                                                                                                                    *
 * @author  Myungsun Kim                                                                                                                                              *
 * @date    Feb. 10, 2021                                                                                                                                             *
 * @brief   Camenisch-Shoup verifiable encryption                                                                                                                     *
 *                                                                                                                                                                    *
\* ****************************************************************************************************************************************************************** */

#ifndef __CAMENISCH_SHOUP_VE_H
#define __CAMENISCH_SHOUP_VE_H

#include <iostream>
#include <string>
#include <chrono>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include <cstddef>
#include <cassert>

#include <openssl/sha.h>

#include <NTL/ZZ.h>

using namespace NTL;

namespace _CS {
    const uint32_t  OPENSSL_HASH_BYTES = 32;
    const uint32_t  primeBits = 768;
    const uint32_t  zeta = 2; //zeta = 2
    const uint32_t  s = 4; //zeta = 4
    const uint32_t  BBits = 768;//BBits=768
    const uint32_t  t = 1;
    const ZZ        one(1);
    const ZZ        zero(0);
    const ZZ        two(2);
    const ZZ        BBit(BBits);

    /**
     * @brief key structures
     *
     */
    struct _public_key {
        ZZ  _N,
            _N_prime,
            _N_zeta_1,
            _N_zeta;
        ZZ  _T;      //< T = (1 + N) mod N_zeta_1    
        ZZ  _g;
        ZZ  _y;
        ZZ  _2_B;
        
    };

    struct _secret_key {
        ZZ  x;
            //x_2_;
            //x_3_;
    };

    struct _commit_key {
        std::vector<ZZ> g;
        ZZ h;
    };

    /**
     * @brief Camenisch-Shoup ciphertext
     *
     */
    struct _ciphertext {
        ZZ  _u,
            _e;
            //_e2;
            //_v;
    };

    /**
     * @brief   encryption proof
     *
     */
    struct _proof {
        ZZ      _v_tilde,                   //< commitments
            _u_prime,
            _e1_prime,
            _e2_prime,
            _v_tilde_prime;
        uint8_t _c[OPENSSL_HASH_BYTES];    //< challenge
        ZZ      _r_bar,
            _m11_bar,
            _m12_bar,
            _m13_bar,
            _m21_bar,
            _m22_bar,
            _m23_bar,
            _s_bar;                     //< responses
    };

    void _random_oracle(uint8_t challenge[OPENSSL_HASH_BYTES], const std::string& input);
    void _setup(_public_key& pk);
    void _key_generation(_public_key& pk, _secret_key& sk,_commit_key &ck);
    //void //_prove_encryption(const _public_key pk, std::vector<ZZ>& message,ZZ &r,_ciphertext& c);
        //_encrypt(const _public_key pk, std::vector<ZZ>& message, ZZ& r, _ciphertext& c)
    void _encrypt(const _public_key pk, std::vector<ZZ>& message, ZZ& r, _ciphertext& c);
    void _L(_public_key pk, const ZZ m, ZZ& t1);
    void _DLog(_public_key pk, const ZZ m, ZZ& p);
    void _decrypt(const _public_key pk, const _secret_key sk, const _ciphertext& c, std::vector<ZZ>& message);
    void CS_HomoAdd(const _public_key pk, _ciphertext& ciphertext_result, _ciphertext ciphertext_1, _ciphertext ciphertext_2);
    void CS_Homosub(const _public_key pk, _ciphertext& ciphertext_result, _ciphertext ciphertext_1, _ciphertext ciphertext_2);
    void CS_ScalarMul(const _public_key pk, _ciphertext& ciphertext_result, _ciphertext ciphertext_1, ZZ k);
    //bool //_verify_encryption(const _public_key pk, const _ciphertext cpx, const _proof pf,int i);
}

#endif //__CAMENISCH_SHOUP_VE_H