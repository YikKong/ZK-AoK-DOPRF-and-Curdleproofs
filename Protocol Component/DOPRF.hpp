#pragma once
#include "global.hpp"
#include "hash.hpp"
#include "print.hpp"
#include "routines.hpp"
#include "calculate_dlog.hpp"
#include "Camenisch-shoup.h"
#include <iostream>
#include <string>
#include <chrono>
#include <fstream>
#include <sstream>
#include <cstddef>
#include <cassert>
#include <openssl/sha.h>
#include <NTL/ZZ.h>
using namespace NTL;

namespace _CS {
    /**
     * @brief random oracle
     *
     */
    void
        _random_oracle(uint8_t challenge[OPENSSL_HASH_BYTES], const std::string& input)
    {
        SHA256_CTX ctx;
        SHA256_Init(&ctx);
        SHA256_Update(&ctx, input.data(), input.size());
        SHA256_Final(challenge, &ctx);
        return;
    }

    /**
     * @brief setup
     *
     */
    void
        _setup(_public_key& pk)
    {
        printf("setup\n");
        NTL::ZZ p, q;
        //< safe primes
        GenGermainPrime(p, primeBits);
        GenGermainPrime(q, primeBits);
        //< the RSA modulus
        pk._N = p * q;
        pk._N_prime = ((p - 1) / 2) * ((q - 1) / 2);
        pk._N_zeta = 1;
        pk._N_zeta_1 = pk._N;
        // n_sqr=n^s  n_sAdd1=n^(s+1)
        for (int i = 1; i <= zeta; ++i) {
            pk._N_zeta = pk._N_zeta * pk._N;
            pk._N_zeta_1 = pk._N_zeta_1 * pk._N;
        }
        //< a default generator
        //< zeta = (1 + n) mod n^(s+1)
        pk._T = AddMod(one, pk._N, pk._N_zeta_1);
        pk._2_B = PowerMod(two, BBit - 5, pk._N_zeta_1);
        return;
    }

    /**
     * @brief key generation
     *
     */
    void
        _key_generation(_public_key& pk, _secret_key& sk, _commit_key& ck)
    {
        NTL::ZZ g_prime;
#ifdef _DEBUG

#endif
        //< choose randoms secret key
        RandomBits(sk.x, 2 * primeBits - 2);
#ifdef _DEBUG
#endif
        while (true) {
            RandomLen(g_prime, 2 * primeBits - 2);
            NTL::ZZ d = GCD(g_prime, pk._N_zeta_1);
            if (1 == d) {
                break;
            }
        }
#ifdef _DEBUG
#endif    
        pk._g = PowerMod(g_prime, 2 * pk._N_zeta, pk._N_zeta_1);
        pk._y = PowerMod(pk._g, sk.x, pk._N_zeta_1);

        // Pedersen commitment key: gi/hΪRSAģN��ѭ����Ⱥ������Ԫ
        for (int i = 0; i < s; ++i) {
            NTL::ZZ gtemp;
            while (true) {
                RandomBits(gtemp, 2 * primeBits - 2);
                NTL::ZZ a = PowerMod(gtemp, pk._N_prime, pk._N);
                if (1 == a) {
                    break;
                }
            }
            ck.g.emplace_back(gtemp);
        }
        while (true) {
            RandomBits(ck.h, 2 * primeBits - 2);
            NTL::ZZ a = PowerMod(ck.h, pk._N_prime, pk._N);
            if (1 == a) {
                break;
            }
        }
        return;
    }

    /**
     * @brief encryption
     *
     */

    void _encrypt(const _public_key pk, std::vector<NTL::ZZ>& message, NTL::ZZ& r, _ciphertext& c)
    {
        if (message.size() > 4) {
            std::cout << "encrypt message size error " << std::endl;
            return;
        }
        NTL::ZZ          m, etemp;

        m = zero;
        RandomBits(r, 2 * primeBits - 2);
        for (int i = 0; i < message.size(); ++i) {
            NTL::ZZ m_add_weight;
            mul(m_add_weight, message[i], power(pk._2_B, i));
            add(m, m, m_add_weight);
        }
        c._u = PowerMod(pk._g, r, pk._N_zeta_1);
        c._e = PowerMod(pk._y, r, pk._N_zeta_1);
        etemp = PowerMod(pk._T, m, pk._N_zeta_1);
        c._e = MulMod(c._e, etemp, pk._N_zeta_1);

        return;
    }

    void _encrypt(const _public_key pk, std::vector<NTL::ZZ>& message, NTL::ZZ& r, _ciphertext& c,bool flag)
    {
        if (message.size() > 4) {
            std::cout << "encrypt message size error " << std::endl;
            return;
        }
        NTL::ZZ          m, etemp;

        m = zero;
        for (int i = 0; i < message.size(); ++i) {
            NTL::ZZ m_add_weight;
            mul(m_add_weight, message[i], power(pk._2_B, i));
            add(m, m, m_add_weight);
        }
        c._u = PowerMod(pk._g, r, pk._N_zeta_1);
        c._e = PowerMod(pk._y, r, pk._N_zeta_1);
        etemp = PowerMod(pk._T, m, pk._N_zeta_1);
        c._e = MulMod(c._e, etemp, pk._N_zeta_1);

        return;
    }

    void _L(_public_key pk, const NTL::ZZ m, NTL::ZZ& t1) {
        NTL::ZZ f;
        f = m - 1;
        div(t1, f, pk._N);
        return;
    }

    void _DLog(_public_key pk, const NTL::ZZ m, NTL::ZZ& p) {
        NTL::ZZ i, j, k, t1, t2, nj, nk1, k1, t3, t4;
        int a;
        i = zero;
        for (j = 1; j <= zeta; ++j) {
            nj = one;
            nk1 = one;
            k1 = one;
            _L(pk, m, t1);
            t2 = i;
            for (a = 1; a <= j; ++a)
            {
                nj = nj * pk._N;
            }
            for (k = 2; k <= j; ++k) {
                for (a = 1; a < k; ++a)
                {
                    nk1 = nk1 * pk._N;
                    k1 = a * k1;
                }
                k1 = k1 * k;
                i = i - 1;
                t2 = MulMod(t2, i, nj);
                t3 = MulMod(t2, nk1, nj);
                t4 = InvMod(k1, nj);
                t3 = MulMod(t3, t4, nj);
                t1 = SubMod(t1, t3, nj);
            }
            i = t1;
        }
        p = i;
        return;
    }


    /**
     * @brief decryption
     *
     */
    void
        _decrypt(const _public_key pk, const _secret_key sk, const _ciphertext& c, std::vector<NTL::ZZ>& message)
    {
        NTL::ZZ          z, m;
        //< recover the message
        // (1+N)^m = e/(u^x) mod n^(s+1)
        z = InvMod(c._u, pk._N_zeta_1);
        z = PowerMod(z, sk.x, pk._N_zeta_1);
        z = MulMod(c._e, z, pk._N_zeta_1);

        // m = DLog((1+N)^m)
        _DLog(pk, z, m);
        for (int i = 0; i < s; ++i) {
            NTL::ZZ mtemp;
            div(mtemp, m, power(pk._2_B, s - i - 1));
            message.emplace_back(mtemp);
            mul(mtemp, mtemp, power(pk._2_B, s - i - 1));
            sub(m, m, mtemp);
        }
        std::reverse(message.begin(), message.end());
        return;
    }

    void CS_HomoAdd(const _public_key pk, _ciphertext& ciphertext_result, _ciphertext ciphertext_1, _ciphertext ciphertext_2) {
        ciphertext_result._u = MulMod(ciphertext_1._u, ciphertext_2._u, pk._N_zeta_1);
        ciphertext_result._e = MulMod(ciphertext_1._e, ciphertext_2._e, pk._N_zeta_1);
        return;
    }

    void CS_Homosub(const _public_key pk, _ciphertext& ciphertext_result, _ciphertext ciphertext_1, _ciphertext ciphertext_2) {
        ciphertext_2._u = PowerMod(ciphertext_2._u, -1, pk._N_zeta_1);
        ciphertext_2._e = PowerMod(ciphertext_2._e, -1, pk._N_zeta_1);
        ciphertext_result._u = MulMod(ciphertext_1._u, ciphertext_2._u, pk._N_zeta_1);
        ciphertext_result._e = MulMod(ciphertext_1._e, ciphertext_2._e, pk._N_zeta_1);
        return;
    }

    void CS_ScalarMul(const _public_key pk, _ciphertext& ciphertext_result, _ciphertext ciphertext_1, NTL::ZZ k) {
        ZZ temp;
        ciphertext_result._u = PowerMod(ciphertext_1._u, k, pk._N_zeta_1);
        ciphertext_result._e = PowerMod(ciphertext_1._e, k, pk._N_zeta_1);
        return;
    }

}
//ʮ�����ַ���ת��Ϊʮ���Ƶ�ZZ���ݣ�BN_bn2dec(bnm)��ʮ�����Ƶ�BIGNUM����ת��Ϊʮ���Ƶ��ַ���
NTL::ZZ stringToZZ(std::string str)
{
    NTL::ZZ number = conv<NTL::ZZ>(str[0]) - 48;
    long len = str.length();
    for (long i = 1; i < len; i++)
    {
        number *= 10;
        number += conv<NTL::ZZ>(str[i]) - 48;
    }

    return number;
}

//ʮ���Ƶ�ZZ���ݷ��ض�Ӧʮ��������ɵ��ַ��������ͨ��BN_dec2bn(bn(�洢BIGNUM* ���ݵĵ�ַ), &str[0](�ַ����׵�ַ))
//���ת��Ϊʮ�����Ƶ�BIGNUM����
std::string ZZToString(NTL::ZZ num)
{

    std::vector<char> str;
    long t;

    while (num != 0)
    {
        conv(t, num % 10);
        char temp = int(t) + '0';
        str.push_back(temp);
        num /= 10;
    }
    reverse(str.begin(), str.end());
    std::string str1(str.begin(), str.end());
    return str1;
}

BIGNUM* ZZtoBIGNUM(NTL::ZZ m) {
    std::string str = ZZToString(m);
    BIGNUM* bnm = BN_new();
    BIGNUM** bn = &bnm;
    BN_dec2bn(bn, &str[0]);
    return bnm;
}

NTL::ZZ BIGNUMtoZZ(BIGNUM* m) {
    std::string str = BN_bn2dec(m);
    return stringToZZ(str);
}

NTL::ZZ BIGNUMtoZZ(const BIGNUM* m) {
    std::string str = BN_bn2dec(m);
    return stringToZZ(str);
}

struct EG_Public
{
    EC_POINT* PRF_g;
    EC_POINT* EG_g;
    EC_POINT* EXP_EG_g;
    vector<EC_POINT*> com_g;
    EC_POINT* com_h;
};

struct EG_Keypair
{
    EC_POINT* pk;  // define pk
    BIGNUM* sk;    // define sk
};

struct EG_Cipher
{
    EC_POINT* X; // X = pk��r 
    EC_POINT* Y; // Y = g��r+m
};

void EG_Public_new(EG_Public& eg_pp)
{
    eg_pp.PRF_g = EC_POINT_new(group);
    eg_pp.EG_g = EC_POINT_new(group);
    eg_pp.EXP_EG_g = EC_POINT_new(group);
    eg_pp.com_g.resize(_CS::s, NULL);
    ECP_vec_new(eg_pp.com_g);
    eg_pp.com_h = EC_POINT_new(group);
}

/* free memory of PP */
void EG_Public_free(EG_Public& eg_pp)
{
    EC_POINT_free(eg_pp.PRF_g);
    EC_POINT_free(eg_pp.EG_g);
    EC_POINT_free(eg_pp.EXP_EG_g);
    ECP_vec_free(eg_pp.com_g);
    EC_POINT_free(eg_pp.com_h);
}

void EG_Keypair_new(EG_Keypair& eg_kp)
{
    eg_kp.pk = EC_POINT_new(group);
    eg_kp.sk = BN_new();
}

void EG_Keypair_free(EG_Keypair& eg_kp)
{
    EC_POINT_free(eg_kp.pk);
    BN_free(eg_kp.sk);
}

void EG_Cipher_new(EG_Cipher& eg_cipher)
{
    eg_cipher.X = EC_POINT_new(group);
    eg_cipher.Y = EC_POINT_new(group);
}

void EG_Cipher_free(EG_Cipher& eg_cipher)
{
    EC_POINT_free(eg_cipher.X);
    EC_POINT_free(eg_cipher.Y);
}

void EG_Public_setup(EG_Public& eg_pp)
{
    ECP_random(eg_pp.PRF_g);
    ECP_random(eg_pp.EG_g);
    ECP_random(eg_pp.EXP_EG_g);
    return;
}

void EG_Keypair_keygen(EG_Public eg_pp, EG_Keypair& eg_kp)
{
    BN_random(eg_kp.sk); // sk \sample Z_p
    EC_POINT_mul(group, eg_kp.pk, NULL, eg_pp.EG_g, eg_kp.sk, bn_ctx); // pk = g��sk  
    return;
}

void ElGamal_Enc(EG_Public eg_pp,
    EC_POINT*& m,
    BIGNUM*& r,
    EC_POINT*& pk,
    EG_Cipher& eg_cipher)
{
    EC_POINT_mul(group, eg_cipher.X, NULL, eg_pp.EG_g, r, bn_ctx);
    EC_POINT_mul(group, eg_cipher.Y, NULL, pk, r, bn_ctx);
    EC_POINT_add(group, eg_cipher.Y, eg_cipher.Y, m, bn_ctx);
    BN_free(r);
}

void ElGamal_Dec(EG_Public eg_pp,
    BIGNUM*& sk,
    EG_Cipher eg_cipher,
    EC_POINT*& m)
{
    EC_POINT_mul(group, m, NULL, eg_cipher.X, sk, bn_ctx);
    EC_POINT_invert(group, m, bn_ctx);
    EC_POINT_add(group, m, eg_cipher.Y, m, bn_ctx);
}

//ȺԪ�صļ��ܣ�����̬ͬ����Ϊm_1+m_2
void ElGamal_HomoAdd(EG_Cipher& eg_cipher_result, EG_Cipher eg_cipher1, EG_Cipher eg_cipher2)
{
    EC_POINT_add(group, eg_cipher_result.X, eg_cipher1.X, eg_cipher2.X, bn_ctx);
    EC_POINT_add(group, eg_cipher_result.Y, eg_cipher1.Y, eg_cipher2.Y, bn_ctx);
}

/* homomorphic sub */
void ElGamal_HomoSub(EG_Cipher& eg_cipher_result, EG_Cipher eg_cipher1, EG_Cipher eg_cipher2)
{
    EC_POINT_sub(eg_cipher_result.X, eg_cipher1.X, eg_cipher2.X);
    EC_POINT_sub(eg_cipher_result.Y, eg_cipher1.Y, eg_cipher2.Y);
}

// �����˷����� ���Լ���ȺԪ�أ�����Ԫ����Կ������һ��ֵ
void ElGamal_ScalarMul(EG_Cipher& eg_cipher_result, EG_Cipher eg_cipher, BIGNUM*& k)
{
    EC_POINT_mul(group, eg_cipher_result.X, NULL, eg_cipher.X, k, bn_ctx);
    EC_POINT_mul(group, eg_cipher_result.Y, NULL, eg_cipher.Y, k, bn_ctx);
}

struct Party {
    vector<vector<long>> input_set;
    vector<NTL::ZZ> input_set_com;
    vector<NTL::ZZ> input_set_com_random;
    vector<long> input_related_set;
    NTL::ZZ prf_key;
    _CS::_public_key cs_pk;
    _CS::_secret_key cs_sk;
    _CS::_commit_key cs_ck;
    EG_Keypair eg_kp;
};

void cs_com_gen(const Party& party, vector<NTL::ZZ> message, NTL::ZZ r, vector<NTL::ZZ>& cs_com) {
    NTL::ZZ com;
    com = PowerMod(party.cs_ck.h, r, party.cs_pk._N * party.cs_pk._N);
    NTL::ZZ comtemp;
    for (int i = 0; i < message.size(); ++i) {
        comtemp = PowerMod(party.cs_ck.g[i], message[i], party.cs_pk._N*party.cs_pk._N);
        com = MulMod(com, comtemp, party.cs_pk._N * party.cs_pk._N);
    }
    cs_com.emplace_back(com);
    return;
}

void cs_com_gen(const Party& party, vector<NTL::ZZ> message, NTL::ZZ r, NTL::ZZ& cs_com) {
    NTL::ZZ com;
    com = PowerMod(party.cs_ck.h, r, party.cs_pk._N * party.cs_pk._N);
    NTL::ZZ comtemp;

    for (int i = 0; i < message.size(); ++i) {
        comtemp = PowerMod(party.cs_ck.g[i], message[i], party.cs_pk._N * party.cs_pk._N);
        com = MulMod(com, comtemp, party.cs_pk._N * party.cs_pk._N);
    }
    cs_com = com;
    return;
}

void cs_com_gen(const Party& party, vector<long> message, NTL::ZZ r, vector<NTL::ZZ>& cs_com) {
    
    NTL::ZZ com;
    com = PowerMod(party.cs_ck.h, r, party.cs_pk._N * party.cs_pk._N);
    NTL::ZZ comtemp;
    for (int i = 0; i < message.size(); ++i) {
        comtemp = PowerMod(party.cs_ck.g[i], message[i], party.cs_pk._N * party.cs_pk._N);
        com = MulMod(com, comtemp, party.cs_pk._N * party.cs_pk._N);
    }
    cs_com.emplace_back(com);
    return;
}

void eg_com_gen(const EG_Public eg_pp,const Party party,vector<NTL::ZZ>message,BIGNUM* r,EC_POINT* eg_com) {
    EC_POINT_mul(group, eg_com, NULL, eg_pp.com_h, r, bn_ctx);
    EC_POINT* eg_com_temp=EC_POINT_new(group);
    BIGNUM* message_temp = BN_new();
    for (int i = 0; i < message.size(); ++i) {
        message_temp = ZZtoBIGNUM(message[i]);
        EC_POINT_mul(group, eg_com_temp, NULL, eg_pp.com_g[i], message_temp, bn_ctx);
        EC_POINT_add(group, eg_com, eg_com, eg_com_temp, bn_ctx);
    }
    EC_POINT_free(eg_com_temp);
    BN_free(message_temp);
    return;
}

void Party_setup(EG_Public& eg_pp, Party& p) {
    RandomBits(p.prf_key, BN_num_bits(order) - 10);
    _CS::_setup(p.cs_pk);
    _CS::_key_generation(p.cs_pk,p.cs_sk,p.cs_ck);

    p.eg_kp.pk = EC_POINT_new(group);
    p.eg_kp.sk = BN_new();
    EG_Keypair_keygen(eg_pp, p.eg_kp);

    int begin = 0;
    while (begin < p.input_set.size()) {
        NTL::ZZ comtemp;
        NTL::ZZ r;
        RandomBits(r, NumBits(p.cs_pk._N) - 10);
        p.input_set_com_random.emplace_back(r);
        cs_com_gen(p, p.input_set[begin], r, p.input_set_com);
        begin++;
    }
    return;
}

void Party_free(Party& p) {
    EG_Keypair_free(p.eg_kp);
    return;
}

struct Round_1_witness {
    NTL::ZZ prf_key; //round_1
    NTL::ZZ com_random; //round_1
    NTL::ZZ ciphertext_random; //round_1
};

struct Round_1_proof_random {
    NTL::ZZ prf_key_random; //round_1_proof
    NTL::ZZ com_random_random; //round_1_proof
    NTL::ZZ ciphertext_random_random; //round_1_proof
};

struct Round_1_message {
    _CS::_ciphertext k_ciphertext; //round_1
    _CS::_ciphertext k_ciphertext_random; //round_1_proof
    NTL::ZZ k_com; //round_1
    NTL::ZZ k_com_random; //round_1_proof
    NTL::ZZ prf_key_random; //round_1_proof
    NTL::ZZ com_random_random; //round_1_proof
    NTL::ZZ ciphertext_random_random; //round_1_proof
    
    NTL::ZZ c;
};

void round_1_proof(const Party party, 
    Round_1_witness& round_1_witness,
    Round_1_proof_random& round_1_proof_random, 
    Round_1_message& round_1_message) {
    std::cout << "round_1_proof " << std::endl;
    RandomBits(round_1_proof_random.prf_key_random, BN_num_bits(order) - 10);
    RandomBits(round_1_proof_random.com_random_random, NumBits(party.cs_pk._N) - 10);
    RandomBits(round_1_proof_random.ciphertext_random_random, NumBits(party.cs_pk._N) - 10);
    vector<NTL::ZZ> prf_key_temp;
    prf_key_temp.emplace_back(round_1_proof_random.prf_key_random);
    // k_ciphertext  ciphertext_Random
    _encrypt(party.cs_pk, prf_key_temp, round_1_proof_random.ciphertext_random_random, round_1_message.k_ciphertext_random,true);
    for (int i = 1; i < _CS::s; ++i) {
        prf_key_temp.emplace_back(round_1_proof_random.prf_key_random);
    }
    cs_com_gen(party, prf_key_temp, round_1_proof_random.com_random_random, round_1_message.k_com_random);

    size_t      vlen = 0;
    uint8_t* vstr = 0;
    std::string input;
    uint8_t c[_CS::OPENSSL_HASH_BYTES];
    vlen = NumBytes(round_1_message.k_com);
    vstr = new uint8_t[vlen];
    std::memset(vstr, 0, sizeof(uint8_t) * vlen);
    NTL::BytesFromZZ(vstr, round_1_message.k_com, vlen);
    input.append(reinterpret_cast<const char*>(vstr), vlen);
    _CS::_random_oracle(c, input);
    round_1_message.c = NTL::ZZFromBytes(c, _CS::OPENSSL_HASH_BYTES);
    round_1_message.prf_key_random = round_1_proof_random.prf_key_random + round_1_message.c * round_1_witness.prf_key;
    round_1_message.com_random_random = round_1_proof_random.com_random_random + round_1_message.c * round_1_witness.com_random;
    round_1_message.ciphertext_random_random = round_1_proof_random.ciphertext_random_random + round_1_message.c * round_1_witness.ciphertext_random;
}

bool round_1_verify(
    const Party party,
    Round_1_witness& round_1_witness,
    Round_1_proof_random& round_1_proof_random,
    Round_1_message& round_1_message)
{
    std::cout << "round_1_verify " << std::endl;
    vector<NTL::ZZ> result;
    _CS::_ciphertext verify_ciphertext_left, verify_ciphertext_right;
    vector<NTL::ZZ> verify_message_left;
    verify_message_left.emplace_back(round_1_message.prf_key_random);
    _CS::_encrypt(party.cs_pk, verify_message_left, round_1_message.ciphertext_random_random, verify_ciphertext_left,true);
    _CS::CS_ScalarMul(party.cs_pk, verify_ciphertext_right, round_1_message.k_ciphertext, round_1_message.c);
    _CS::CS_HomoAdd(party.cs_pk, verify_ciphertext_right, verify_ciphertext_right, round_1_message.k_ciphertext_random);
    result.emplace_back(SubMod(verify_ciphertext_left._u, verify_ciphertext_right._u, party.cs_pk._N_zeta_1));
    result.emplace_back(SubMod(verify_ciphertext_left._e, verify_ciphertext_right._e, party.cs_pk._N_zeta_1));

    for (int i = 1; i < _CS::s; ++i) {
        verify_message_left.emplace_back(round_1_message.prf_key_random);
    }
    
    NTL::ZZ verify_com_left, verify_com_right;
    cs_com_gen(party, verify_message_left, round_1_message.com_random_random, verify_com_left);
    verify_com_right = PowerMod(round_1_message.k_com, round_1_message.c, party.cs_pk._N * party.cs_pk._N);
    verify_com_right = MulMod(verify_com_right, round_1_message.k_com_random, party.cs_pk._N * party.cs_pk._N);
    result.emplace_back(SubMod(verify_com_left, verify_com_right, party.cs_pk._N * party.cs_pk._N));


    for (int i = 0; i < result.size(); ++i) {
        if (result[i] != 0) {
            return false;
        }
    }
    return true;
}

void round_1(const Party party,
    Round_1_witness& round_1_witness,
    Round_1_proof_random& round_1_proof_random,
    Round_1_message& round_1_message) {
    std::cout << "round_1 " << std::endl;
    round_1_witness.prf_key = party.prf_key;
    // prf_key_random
    
    
    vector<NTL::ZZ> prf_key_temp;
    prf_key_temp.emplace_back(round_1_witness.prf_key);
    // k_ciphertext  ciphertext_Random
    _encrypt(party.cs_pk, prf_key_temp, round_1_witness.ciphertext_random, round_1_message.k_ciphertext);
    for (int i = 1; i < _CS::s; ++i) {
        prf_key_temp.emplace_back(round_1_witness.prf_key);
    }


    RandomBits(round_1_witness.com_random, NumBits(party.cs_pk._N) - 10);
    cs_com_gen(party, prf_key_temp, round_1_witness.com_random, round_1_message.k_com);
    
    round_1_proof(party, round_1_witness, round_1_proof_random, round_1_message);
    return;
}

struct Round_2_witness {
    vector<vector<NTL::ZZ>> a, b, alpha; //round_2
    vector<NTL::ZZ> a_com_random, b_com_random, alpha_com_random; //round_2
    vector<NTL::ZZ> beta_ciphertext_random; //round_2
};

struct Round_2_proof_random {
    vector<NTL::ZZ> random_a, random_b, random_alpha; //round_2_proof
    NTL::ZZ random_a_com_random, random_b_com_random, random_alpha_com_random; //round_2_proof
    NTL::ZZ random_beta_ciphertext_random; //round_2_proof
};

struct Round_2_message {
    vector<_CS::_ciphertext> beta_ciphertext; //round_2
    _CS::_ciphertext random_beta_ciphertext; //round_2_proof

    vector<NTL::ZZ> a_com, b_com, alpha_com; //round_2
    NTL::ZZ random_a_com, random_b_com, random_alpha_com; //round_2_proof

    vector<NTL::ZZ> random_a, random_b, random_alpha; //round_2_proof
    NTL::ZZ random_a_com_random, random_b_com_random, random_alpha_com_random; //round_2_proof
    NTL::ZZ random_beta_ciphertext_random; //round_2_proof

    vector<NTL::ZZ> c; //round_2_proof

    vector<vector<EC_POINT*>> g_a; // round_2
    vector<EC_POINT*> g_random_a; //round_2_proof
};

void round_2_proof(const Party& P_input,
    const Party& P_compute,
    EG_Public& eg_public,
    Round_1_message& P_compute_round_1_message,
    Round_1_witness& P_input_round_1_witness,
    Round_2_witness& round_2_witness,
    Round_2_proof_random& round_2_proof_random,
    Round_2_message& round_2_message) {
    std::cout << "round_2_proof " << std::endl;
    
    
    for (int i = 0; i < _CS::s; ++i) {
        ZZ temp;
        EC_POINT* g_a = EC_POINT_new(group);
        RandomBits(temp, BN_num_bits(order) - 5);
        EC_POINT_mul(group, g_a, NULL, eg_public.PRF_g, ZZtoBIGNUM(temp), bn_ctx);
        round_2_message.g_random_a.emplace_back(g_a);
        round_2_proof_random.random_a.emplace_back(temp);
        round_2_message.random_a.emplace_back(temp);

        RandomBits(temp, BN_num_bits(order) - 5);
        round_2_proof_random.random_b.emplace_back(temp);
        round_2_message.random_b.emplace_back(temp);

        RandomBits(temp, BN_num_bits(order) - 5);
        round_2_proof_random.random_alpha.emplace_back(temp);
        round_2_message.random_alpha.emplace_back(temp);
    }
    RandomBits(round_2_proof_random.random_a_com_random, NumBits(P_compute.cs_pk._N) - 10);
    
    RandomBits(round_2_proof_random.random_b_com_random, NumBits(P_compute.cs_pk._N) - 10);
    RandomBits(round_2_proof_random.random_alpha_com_random, NumBits(P_compute.cs_pk._N) - 10);
    cs_com_gen(P_compute, round_2_proof_random.random_a, round_2_proof_random.random_a_com_random, round_2_message.random_a_com);
    cs_com_gen(P_compute, round_2_proof_random.random_b, round_2_proof_random.random_b_com_random, round_2_message.random_b_com);
    cs_com_gen(P_compute, round_2_proof_random.random_alpha, round_2_proof_random.random_alpha_com_random, round_2_message.random_alpha_com);

    NTL::ZZ k_asum = _CS::zero;
    for (int i = 0; i < _CS::s; ++i) {
        NTL::ZZ m_add_weight;
        mul(m_add_weight, round_2_proof_random.random_a[i], power(P_compute.cs_pk._2_B, i));
        add(k_asum, k_asum, m_add_weight);
    }
    _CS::_ciphertext ciphertext_k_a_sum;
    CS_ScalarMul(P_compute.cs_pk, ciphertext_k_a_sum, P_compute_round_1_message.k_ciphertext, k_asum);
    RandomBits(round_2_proof_random.random_beta_ciphertext_random, NumBits(P_compute.cs_pk._N) - 10);

    vector<NTL::ZZ> random_message_temp;
    for (int i = 0; i < _CS::s; ++i) {
        ZZ add_temp;
        NTL::mul(add_temp, round_2_proof_random.random_b[i], BIGNUMtoZZ(order));
        NTL::add(add_temp, round_2_proof_random.random_alpha[i], add_temp);
        random_message_temp.emplace_back(add_temp);
    }
    _CS::_encrypt(P_compute.cs_pk, random_message_temp, round_2_proof_random.random_beta_ciphertext_random, round_2_message.random_beta_ciphertext,true);
    CS_HomoAdd(P_compute.cs_pk, round_2_message.random_beta_ciphertext, ciphertext_k_a_sum, round_2_message.random_beta_ciphertext);

    round_2_message.random_a_com_random = round_2_proof_random.random_a_com_random;
    round_2_message.random_b_com_random = round_2_proof_random.random_b_com_random;
    round_2_message.random_alpha_com_random = round_2_proof_random.random_alpha_com_random;
    round_2_message.random_beta_ciphertext_random = round_2_proof_random.random_beta_ciphertext_random;

    ZZ temp;
    for (int i = 0; i < P_input.input_set.size(); i++)
    {
        size_t vlen = 0;
        uint8_t* vstr = 0;
        std::string input;
        uint8_t c[_CS::OPENSSL_HASH_BYTES];
        vlen = NumBytes(round_2_message.a_com[i]);
        vstr = new uint8_t[vlen];
        std::memset(vstr, 0, sizeof(uint8_t) * vlen);
        NTL::BytesFromZZ(vstr, round_2_message.a_com[i], vlen);
        input.append(reinterpret_cast<const char*>(vstr), vlen);
        _CS::_random_oracle(c, input);
        round_2_message.c.emplace_back(NTL::ZZFromBytes(c, _CS::OPENSSL_HASH_BYTES));   
        round_2_message.random_a_com_random = round_2_message.random_a_com_random + round_2_message.c[i] * round_2_witness.a_com_random[i];
        round_2_message.random_b_com_random = round_2_message.random_b_com_random + round_2_message.c[i] * round_2_witness.b_com_random[i];
        round_2_message.random_alpha_com_random = round_2_message.random_alpha_com_random + round_2_message.c[i] * round_2_witness.alpha_com_random[i];
        round_2_message.random_beta_ciphertext_random = round_2_message.random_beta_ciphertext_random + round_2_message.c[i] * round_2_witness.beta_ciphertext_random[i];
        temp = temp + round_2_message.c[i] * round_2_witness.beta_ciphertext_random[i];
        
    }

    for (int i = 0; i < _CS::s; ++i) {
        NTL::ZZ temp1,temp2,temp3;
        for (int j = 0; j < P_input.input_set.size(); ++j) {
            temp1 = temp1 + round_2_message.c[j] * round_2_witness.a[j][i];
            temp2 = temp2 + round_2_message.c[j] * round_2_witness.b[j][i];
            temp3 = temp3 + round_2_message.c[j] * round_2_witness.alpha[j][i];
            round_2_message.random_a[i] = round_2_message.random_a[i] + round_2_message.c[j] * round_2_witness.a[j][i];
            round_2_message.random_b[i] = round_2_message.random_b[i] + round_2_message.c[j] * round_2_witness.b[j][i];
            round_2_message.random_alpha[i] = round_2_message.random_alpha[i] + round_2_message.c[j] * round_2_witness.alpha[j][i];
        }
    }

}

bool round_2_verify(const Party& P_input,
    const Party& P_compute,
    EG_Public& eg_public,
    Round_1_message& P_compute_round_1_message,
    Round_1_witness& P_input_round_1_witness,
    Round_2_witness& round_2_witness,
    Round_2_proof_random& round_2_proof_random,
    Round_2_message& round_2_message) {
    std::cout << "round_2_verify " << std::endl;
    vector<NTL::ZZ> result;
    NTL::ZZ verify_a_com_left, verify_b_com_left, verify_alpha_com_left;
    NTL::ZZ verify_a_com_right, verify_b_com_right, verify_alpha_com_right;
    cs_com_gen(P_compute, round_2_message.random_a, round_2_message.random_a_com_random, verify_a_com_left);
    cs_com_gen(P_compute, round_2_message.random_b, round_2_message.random_b_com_random, verify_b_com_left);
    cs_com_gen(P_compute, round_2_message.random_alpha, round_2_message.random_alpha_com_random, verify_alpha_com_left);

    verify_a_com_right = round_2_message.random_a_com;
    verify_b_com_right = round_2_message.random_b_com;
    verify_alpha_com_right = round_2_message.random_alpha_com;
    for (int i = 0; i < P_input.input_set.size(); ++i) {
        verify_a_com_right = MulMod(verify_a_com_right, PowerMod(round_2_message.a_com[i], round_2_message.c[i], P_compute.cs_pk._N * P_compute.cs_pk._N), P_compute.cs_pk._N * P_compute.cs_pk._N);
        verify_b_com_right = MulMod(verify_b_com_right, PowerMod(round_2_message.b_com[i], round_2_message.c[i], P_compute.cs_pk._N * P_compute.cs_pk._N), P_compute.cs_pk._N * P_compute.cs_pk._N);
        verify_alpha_com_right = MulMod(verify_alpha_com_right, PowerMod(round_2_message.alpha_com[i], round_2_message.c[i], P_compute.cs_pk._N * P_compute.cs_pk._N), P_compute.cs_pk._N * P_compute.cs_pk._N);
    }
    result.emplace_back(SubMod(verify_a_com_left, verify_a_com_right, P_compute.cs_pk._N * P_compute.cs_pk._N));
    result.emplace_back(SubMod(verify_b_com_left, verify_b_com_right, P_compute.cs_pk._N * P_compute.cs_pk._N));
    result.emplace_back(SubMod(verify_alpha_com_left, verify_alpha_com_right, P_compute.cs_pk._N * P_compute.cs_pk._N));
    _CS::_ciphertext verify_ciphertext_left, verify_ciphertext_right;
    NTL::ZZ k_asum = _CS::zero;
    for (int i = 0; i < _CS::s; ++i) {
        NTL::ZZ m_add_weight;
        mul(m_add_weight, round_2_message.random_a[i], power(P_compute.cs_pk._2_B, i));
        add(k_asum, k_asum, m_add_weight);
    }
    _CS::_ciphertext ciphertext_k_a_sum;
    CS_ScalarMul(P_compute.cs_pk, ciphertext_k_a_sum, P_compute_round_1_message.k_ciphertext, k_asum);

    vector<NTL::ZZ> random_message_temp;
    for (int i = 0; i < _CS::s; ++i) {
        ZZ add_temp;
        NTL::mul(add_temp, round_2_message.random_b[i], BIGNUMtoZZ(order));
        NTL::add(add_temp, round_2_message.random_alpha[i], add_temp);
        random_message_temp.emplace_back(add_temp);
    }

    _CS::_encrypt(P_compute.cs_pk, random_message_temp, round_2_message.random_beta_ciphertext_random, verify_ciphertext_left,true);
    CS_HomoAdd(P_compute.cs_pk, verify_ciphertext_left, ciphertext_k_a_sum, verify_ciphertext_left);

    verify_ciphertext_right = round_2_message.random_beta_ciphertext;
    for (int i = 0; i < round_2_message.beta_ciphertext.size(); ++i) {
        _CS::_ciphertext ciphertext_temp;
        _CS::CS_ScalarMul(P_compute.cs_pk, ciphertext_temp, round_2_message.beta_ciphertext[i], round_2_message.c[i]);
        _CS::CS_HomoAdd(P_compute.cs_pk, verify_ciphertext_right, ciphertext_temp, verify_ciphertext_right);
    }
    result.emplace_back(SubMod(verify_ciphertext_left._u, verify_ciphertext_right._u, P_compute.cs_pk._N_zeta_1));
    result.emplace_back(SubMod(verify_ciphertext_left._e, verify_ciphertext_right._e, P_compute.cs_pk._N_zeta_1));

    for (int i = 0; i < _CS::s; ++i) {
        EC_POINT* g_a_left = EC_POINT_new(group);
        EC_POINT* g_a_right = EC_POINT_new(group);
        EC_POINT_mul(group, g_a_left, NULL, eg_public.PRF_g, ZZtoBIGNUM(round_2_message.random_a[i]), bn_ctx);
        for (int j = 0; j < round_2_message.g_a.size(); ++j) {
            EC_POINT* g_a_temp = EC_POINT_new(group);
            EC_POINT_mul(group, g_a_temp, NULL, round_2_message.g_a[j][i], ZZtoBIGNUM(round_2_message.c[j]), bn_ctx);
            EC_POINT_add(group, g_a_right, g_a_right, g_a_temp, bn_ctx);
            EC_POINT_free(g_a_temp);
            
        }
        EC_POINT_add(group, g_a_right, g_a_right, round_2_message.g_random_a[i], bn_ctx);
        result.emplace_back(EC_POINT_cmp(group, g_a_left, g_a_right, bn_ctx));
        
        EC_POINT_free(g_a_left);
        EC_POINT_free(g_a_right);
    }
    
    for (int i = 0; i < result.size(); ++i) {
        if (result[i] != _CS::zero) {
            return false;
        }
    }
    return true;
}

void round_2(const Party& P_input,
    const Party& P_compute, 
    EG_Public& eg_public,
    Round_1_message& P_compute_round_1_message,
    Round_1_witness& P_input_round_1_witness, 
    Round_2_witness& round_2_witness, 
    Round_2_proof_random& round_2_proof_random, 
    Round_2_message& round_2_message) {
    std::cout << "round_2 " << std::endl;
    long inputsize = P_input.input_set.size();
    long begin = 0;
    EG_Public_setup(eg_public);
    while (begin < inputsize) {
        vector<NTL::ZZ> atemp,btemp,alphatemp;
        NTL::ZZ temp;
        vector<EC_POINT*> g_a_temp;
        for (int i = 0; i < P_input.input_set[begin].size(); ++i) {
            EC_POINT* g_a = EC_POINT_new(group);
            RandomBits(temp, BN_num_bits(order) - 10);
            EC_POINT_mul(group, g_a, NULL, eg_public.PRF_g, ZZtoBIGNUM(temp), bn_ctx);
            g_a_temp.emplace_back(g_a);
            atemp.emplace_back(temp);
            RandomBits(temp, BN_num_bits(order) - 10);
            btemp.emplace_back(temp);

            NTL::add(temp, P_input.input_set[begin][i], P_input.prf_key);
            NTL:mul(temp, temp, atemp[i]);
            alphatemp.emplace_back(temp); 
        }
        round_2_witness.a.emplace_back(atemp);
        round_2_witness.b.emplace_back(btemp);
        round_2_witness.alpha.emplace_back(alphatemp);
        round_2_message.g_a.emplace_back(g_a_temp);
        NTL::ZZ com_random_temp;
        RandomBits(com_random_temp, NumBits(P_compute.cs_pk._N) - 10);
        round_2_witness.a_com_random.emplace_back(com_random_temp);
        cs_com_gen(P_compute, atemp, com_random_temp, round_2_message.a_com);
        RandomBits(com_random_temp, NumBits(P_compute.cs_pk._N) - 10);
        round_2_witness.b_com_random.emplace_back(com_random_temp);
        cs_com_gen(P_compute, btemp, com_random_temp, round_2_message.b_com);

        NTL::add(com_random_temp, P_input.input_set_com_random[begin], P_input_round_1_witness.com_random);
        NTL::mul(com_random_temp, com_random_temp, round_2_witness.a_com_random[begin]);
        round_2_witness.alpha_com_random.emplace_back(com_random_temp);
        cs_com_gen(P_compute, alphatemp, round_2_witness.alpha_com_random[begin], round_2_message.alpha_com);

        // ct_beta
        NTL::ZZ k_asum = _CS::zero;
        for (int i = 0; i < round_2_witness.a[begin].size(); ++i) {
            NTL::ZZ m_add_weight;
            mul(m_add_weight, round_2_witness.a[begin][i], power(P_compute.cs_pk._2_B, i));
            add(k_asum, k_asum, m_add_weight);
        }
        // r_k2��a_sum
        _CS::_ciphertext ciphertext_k_a_sum;
        CS_ScalarMul(P_compute.cs_pk, ciphertext_k_a_sum, P_compute_round_1_message.k_ciphertext, k_asum);

        vector<NTL::ZZ> message_temp;
        for (int i = 0; i < round_2_witness.a[begin].size(); ++i) {
            ZZ add_temp;
            NTL::mul(add_temp, round_2_witness.b[begin][i], BIGNUMtoZZ(order));
            NTL::add(add_temp, round_2_witness.alpha[begin][i], add_temp);
            message_temp.emplace_back(add_temp);
        }
        _CS::_ciphertext ciphertext_alpha_bq;

        NTL::ZZ ciphertext_random;
        _CS::_encrypt(P_compute.cs_pk, message_temp, ciphertext_random, ciphertext_alpha_bq);
        round_2_witness.beta_ciphertext_random.emplace_back(ciphertext_random);

        _CS::_ciphertext ciphertext_beta;
        _CS::CS_HomoAdd(P_compute.cs_pk, ciphertext_beta, ciphertext_k_a_sum, ciphertext_alpha_bq);
        round_2_message.beta_ciphertext.emplace_back(ciphertext_beta);
       

        begin = begin + 1;
    }

    round_2_proof(P_input, P_compute, eg_public, P_compute_round_1_message, P_input_round_1_witness, round_2_witness, round_2_proof_random, round_2_message);

}

struct Round_3_witness {
    vector<vector<NTL::ZZ>> beta; //round_3
    vector<NTL::ZZ> beta_cs_com_random; //round_3
    vector<BIGNUM*> beta_eg_com_random; //round_3
};

struct Round_3_proof_random {
    vector<NTL::ZZ> random_beta; //round_3_proof
    NTL::ZZ random_beta_cs_com_random; //round_3_proof
    BIGNUM* random_beta_eg_com_random; //round_3_proof
};

struct Round_3_message {
    _CS::_ciphertext random_beta_ciphertext; 
    _CS::_ciphertext random_beta_beta_ciphertext; 
    vector<NTL::ZZ> beta_cs_com; //round_3
    vector<EC_POINT* > beta_eg_com;//round_3
    NTL::ZZ random_beta_cs_com; //round_3_proof
    EC_POINT* random_beta_eg_com; //round_3_proof

    vector<NTL::ZZ> random_beta; //round_3_proof
    NTL::ZZ random_beta_cs_com_random; //round_3_proof
    BIGNUM* random_beta_eg_com_random; //round_3_proof

    vector<NTL::ZZ> c; //round_3_proof
};

void round_3_proof(
    const EG_Public eg_pp,
    const Party& P_input,
    const Party& P_compute,
    Round_2_message& round_2_message,
    Round_3_witness& round_3_witness,
    Round_3_proof_random& round_3_proof_random,
    Round_3_message& round_3_message
) {
    for (int i = 0; i < _CS::s; ++i) {
        ZZ temp;
        RandomBits(temp, BN_num_bits(order) - 5);

        round_3_proof_random.random_beta.emplace_back(temp);
        round_3_message.random_beta.emplace_back(temp);
    }
    RandomBits(round_3_proof_random.random_beta_cs_com_random, NumBits(P_compute.cs_pk._N) - 10);

    round_3_proof_random.random_beta_eg_com_random = BN_new();
    BN_random(round_3_proof_random.random_beta_eg_com_random);
    round_3_message.random_beta_eg_com_random = BN_new();
    BN_copy(round_3_message.random_beta_eg_com_random, round_3_proof_random.random_beta_eg_com_random);

    cs_com_gen(P_compute, round_3_proof_random.random_beta, round_3_proof_random.random_beta_cs_com_random, round_3_message.random_beta_cs_com);

    round_3_message.random_beta_eg_com = EC_POINT_new(group);
    eg_com_gen(eg_pp, P_input, round_3_proof_random.random_beta, round_3_proof_random.random_beta_eg_com_random, round_3_message.random_beta_eg_com);

    round_3_message.random_beta_cs_com_random = round_3_proof_random.random_beta_cs_com_random;

    for (int i = 0; i < P_input.input_set.size(); i++)
    {
        size_t vlen = 0;
        uint8_t* vstr = 0;
        std::string input;
        uint8_t c[_CS::OPENSSL_HASH_BYTES];
        vlen = NumBytes(round_3_message.beta_cs_com[i]);
        vstr = new uint8_t[vlen];
        std::memset(vstr, 0, sizeof(uint8_t) * vlen);
        NTL::BytesFromZZ(vstr, round_3_message.beta_cs_com[i], vlen);
        input.append(reinterpret_cast<const char*>(vstr), vlen);
        _CS::_random_oracle(c, input);
        round_3_message.c.emplace_back(NTL::ZZFromBytes(c, _CS::OPENSSL_HASH_BYTES));
        round_3_message.c[i] = i;
        
        round_3_message.random_beta_cs_com_random = round_3_message.random_beta_cs_com_random + round_3_message.c[i] * round_3_witness.beta_cs_com_random[i];
        
        BIGNUM* temp_bn = BN_new();
        BN_mul(temp_bn, ZZtoBIGNUM(round_3_message.c[i]), round_3_witness.beta_eg_com_random[i], bn_ctx);
        BN_mod_add(round_3_message.random_beta_eg_com_random, 
                  round_3_message.random_beta_eg_com_random, 
                  temp_bn, 
                  order, 
                  bn_ctx);
        BN_free(temp_bn);

    }

    for (int i = 0; i < _CS::s; ++i) {
        for (int j = 0; j < P_input.input_set.size(); ++j) {
            round_3_message.random_beta[i] = round_3_message.random_beta[i] + round_3_message.c[j] * round_3_witness.beta[j][i];
        }
    }
    
}

void round_3(const EG_Public eg_pp,
    const Party& P_input,
    const Party& P_compute,
    Round_2_message& round_2_message,
    Round_3_witness& round_3_witness,
    Round_3_proof_random& round_3_proof_random,
    Round_3_message& round_3_message
) {
    std::cout << "round_3 " << std::endl;
    for (int i = 0; i < round_2_message.beta_ciphertext.size(); ++i) {
        vector<NTL::ZZ> message_temp;
        _CS::_decrypt(P_compute.cs_pk, P_compute.cs_sk, round_2_message.beta_ciphertext[i], message_temp);
        round_3_witness.beta.emplace_back(message_temp);

        NTL::ZZ cs_com_random_temp;
        RandomBits(cs_com_random_temp, NumBits(P_compute.cs_pk._N) - 10);
        round_3_witness.beta_cs_com_random.emplace_back(cs_com_random_temp);
        cs_com_gen(P_input, message_temp, cs_com_random_temp, round_3_message.beta_cs_com);
        
        BIGNUM* beta_eg_random = BN_new();
        EC_POINT* beta_eg_com=EC_POINT_new(group);
        BN_random(beta_eg_random);
        round_3_witness.beta_eg_com_random.emplace_back(beta_eg_random);
        eg_com_gen(eg_pp,P_input,message_temp,beta_eg_random,beta_eg_com);
        round_3_message.beta_eg_com.emplace_back(beta_eg_com);
    }
    round_3_proof(eg_pp, P_input, P_compute, round_2_message, round_3_witness, round_3_proof_random, round_3_message);
}

bool round_3_verify(
    const EG_Public eg_pp,
    const Party& P_input,
    const Party& P_compute,
    Round_2_message& round_2_message,
    Round_3_witness& round_3_witness,
    Round_3_proof_random& round_3_proof_random,
    Round_3_message& round_3_message
) {
    std::cout << "round_3_verify" << std::endl;
    vector<NTL::ZZ> result;

    NTL::ZZ verify_cs_com_left, verify_cs_com_right;
    cs_com_gen(P_compute, round_3_message.random_beta, round_3_message.random_beta_cs_com_random, verify_cs_com_left);

    verify_cs_com_right = round_3_message.random_beta_cs_com;
    for (int i = 0; i < P_input.input_set.size(); ++i) {
        verify_cs_com_right = MulMod(verify_cs_com_right, PowerMod(round_3_message.beta_cs_com[i], round_3_message.c[i], P_compute.cs_pk._N * P_compute.cs_pk._N), P_compute.cs_pk._N * P_compute.cs_pk._N);
    }
    result.emplace_back(SubMod(verify_cs_com_left, verify_cs_com_right, P_compute.cs_pk._N * P_compute.cs_pk._N));


    EC_POINT* verify_eg_com_left = EC_POINT_new(group);
    EC_POINT* verify_eg_com_right = EC_POINT_new(group);
    
    eg_com_gen(eg_pp, P_input, round_3_message.random_beta, round_3_message.random_beta_eg_com_random, verify_eg_com_left);

    EC_POINT_copy(verify_eg_com_right, round_3_message.random_beta_eg_com);
    for (int i = 0; i < round_3_message.beta_eg_com.size(); ++i) {
        EC_POINT* temp_point = EC_POINT_new(group);
        EC_POINT_mul(group, temp_point, NULL, round_3_message.beta_eg_com[i], ZZtoBIGNUM(round_3_message.c[i]), bn_ctx);
        EC_POINT_add(group, verify_eg_com_right, verify_eg_com_right, temp_point, bn_ctx);
        EC_POINT_free(temp_point);
    }
    
    if (EC_POINT_cmp(group, verify_eg_com_left, verify_eg_com_right, bn_ctx) != 0) {
        return false;
    }

    for (int i = 0; i < result.size(); ++i) {
        if (result[i] != _CS::zero) {
            return false;
        }
    }

    EC_POINT_free(verify_eg_com_left);
    EC_POINT_free(verify_eg_com_right);
    return true;
}