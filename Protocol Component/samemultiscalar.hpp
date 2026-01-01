#ifndef __IP__
#define __IP__

#include "global.hpp"
#include "hash.hpp"
#include "print.hpp"
#include "routines.hpp"
#include <algorithm>
#include <random>

struct PedersenCM_KP
{
    size_t VECTOR_LEN;      // denotes the size of witness (witness is upto l = 2^VECTOR_LEN)
    size_t LOG_VECTOR_LEN;  // LOG_VECTOR_LEN = log(VECTOR_LEN) 

    // size of the vector = VECTOR_LEN
    vector<EC_POINT*> vec_g;
    EC_POINT* h;
};

// define the structure of InnerProduct Proof
struct InnerProduct_KP
{
    size_t VECTOR_LEN;      // denotes the size of witness (witness is upto l = 2^VECTOR_LEN)
    size_t LOG_VECTOR_LEN;  // LOG_VECTOR_LEN = log(VECTOR_LEN) 

    // size of the vector = VECTOR_LEN
    vector<EC_POINT*> vec_G;
    vector<EC_POINT*> vec_U;
    vector<EC_POINT*> vec_E;

};


struct InnerProduct_Witness
{
    // size of the vector = VECTOR_LEN
    vector<BIGNUM*> vec_m;
};

struct InnerProduct_Proof
{
    // size of the vector = LOG_VECTOR_LEN

    EC_POINT* BA;
    EC_POINT* BU;
    EC_POINT* BE;
    vector<EC_POINT*> vec_LA;
    vector<EC_POINT*> vec_LU;
    vector<EC_POINT*> vec_LE;
    vector<EC_POINT*> vec_RA;
    vector<EC_POINT*> vec_RU;
    vector<EC_POINT*> vec_RE;
    BIGNUM* m;
};

void PedersenCM_KP_print(PedersenCM_KP& kp)
{

    // size of the vector = VECTOR_LEN
    ECP_vec_print(kp.vec_g, "g");
    ECP_print(kp.h, "h");

};

void InnerProduct_KP_print(InnerProduct_KP& kp)
{

    // size of the vector = VECTOR_LEN
    ECP_vec_print(kp.vec_G, "G");
    ECP_vec_print(kp.vec_G, "U");
    ECP_vec_print(kp.vec_G, "E");

};

void InnerProduct_Witness_print(InnerProduct_Witness& witness)
{
    BN_vec_print(witness.vec_m, "m");
};


void InnerProduct_Proof_print(InnerProduct_Proof& proof)
{
    ECP_vec_print(proof.vec_LA, "LA");
    ECP_vec_print(proof.vec_RA, "RA");
    ECP_vec_print(proof.vec_LU, "LU");
    ECP_vec_print(proof.vec_RU, "RU");
    ECP_vec_print(proof.vec_LE, "LE");
    ECP_vec_print(proof.vec_RE, "RE");
    BN_print(proof.m, "proof.m");
    ECP_print(proof.BA, "BA");
    ECP_print(proof.BU, "BU");
    ECP_print(proof.BE, "BE");
};

void PedersenCM_KP_new(PedersenCM_KP& kp, size_t VECTOR_LEN)
{
    kp.vec_g.resize(VECTOR_LEN-3); ECP_vec_new(kp.vec_g);
    kp.h = EC_POINT_new(group);
}

void PedersenCM_KP_free(PedersenCM_KP& kp)
{
    ECP_vec_free(kp.vec_g);
    EC_POINT_free(kp.h);
}

void InnerProduct_KP_new(InnerProduct_KP& kp, size_t VECTOR_LEN)
{
    kp.vec_G.resize(VECTOR_LEN); ECP_vec_new(kp.vec_G);
    kp.vec_U.resize(VECTOR_LEN); ECP_vec_new(kp.vec_U);
    kp.vec_E.resize(VECTOR_LEN); ECP_vec_new(kp.vec_E);
    
}

void InnerProduct_KP_free(InnerProduct_KP& kp)
{
    ECP_vec_free(kp.vec_G);
    ECP_vec_free(kp.vec_U);
    ECP_vec_free(kp.vec_E);
}

void InnerProduct_Witness_new(InnerProduct_Witness& witness, uint64_t VECTOR_LEN)
{
    witness.vec_m.resize(VECTOR_LEN);
    BN_vec_new(witness.vec_m);
}

void InnerProduct_Witness_free(InnerProduct_Witness& witness)
{
    BN_vec_free(witness.vec_m);
}

void InnerProduct_Proof_new(InnerProduct_Proof& proof)
{
    proof.BA = EC_POINT_new(group);
    proof.BU = EC_POINT_new(group);
    proof.BE = EC_POINT_new(group);
    proof.m = BN_new();
}

void InnerProduct_Proof_free(InnerProduct_Proof& proof)
{
    BN_free(proof.m);
    EC_POINT_free(proof.BA);
    EC_POINT_free(proof.BU);
    EC_POINT_free(proof.BE);

    ECP_vec_free(proof.vec_LA);
    ECP_vec_free(proof.vec_RA);
    ECP_vec_free(proof.vec_LU);
    ECP_vec_free(proof.vec_RU);
    ECP_vec_free(proof.vec_LE);
    ECP_vec_free(proof.vec_RE);

    proof.vec_LA.resize(0);
    proof.vec_RA.resize(0);
    proof.vec_LU.resize(0);
    proof.vec_RU.resize(0);
    proof.vec_LE.resize(0);
    proof.vec_RE.resize(0);
}


/* compute the jth bit of a big integer i (count from little endian to big endian) */
inline uint64_t BN_parse_binary(BIGNUM* BN_i, uint64_t j)
{
    BIGNUM* BN_bit = BN_new();
    BN_copy(BN_bit, BN_i);

    BN_rshift(BN_bit, BN_bit, j);
    BN_mod(BN_bit, BN_bit, BN_2, bn_ctx);

    uint64_t bit;
    if (BN_is_one(BN_bit)) bit = 1;
    else bit = 0;
    BN_free(BN_bit);
    return bit;
}


/* compute the jth bit of a small integer num \in [0, 2^{m-1}] (count from big endian to little endian) */
inline uint64_t int_parse_binary(size_t num, size_t j, size_t m)
{
    size_t cursor = 1 << (m - 1); // set cursor = 2^{m-1} = 1||0...0---(m-1)

    for (auto i = 0; i < j; i++)
    {
        cursor = cursor >> 1;
    }
    if ((num & cursor) != 0) return 1;
    else return 0;
}

/* generate a^n = (a^0, a^1, a^2, ..., a^{n-1}) */
inline void BN_vec_gen_power(vector<BIGNUM*>& result, BIGNUM*& a)
{
    BN_one(result[0]); // set result[0] = 1
    for (auto i = 1; i < result.size(); i++)
    {
        BN_mod_mul(result[i], a, result[i - 1], order, bn_ctx); // result[i] = result[i-1]*a % order
    }
}

/* assign left or right part of a Zn vector */
inline void BN_vec_assign(vector<BIGNUM*>& result, vector<BIGNUM*>& vec_a, string selector)
{
    size_t start_index;
    if (selector == "left") start_index = 0;
    if (selector == "right") start_index = vec_a.size() / 2;

    for (auto i = 0; i < result.size(); i++) {
        BN_copy(result[i], vec_a[start_index + i]);
    }
}

// assign left or right part of an ECn vector
inline void ECP_vec_assign(vector<EC_POINT*>& result, vector<EC_POINT*>& vec_g, string selector)
{
    size_t start_index;
    if (selector == "left") start_index = 0;
    if (selector == "right") start_index = vec_g.size() / 2;

    for (auto i = 0; i < result.size(); i++) {
        EC_POINT_copy(result[i], vec_g[start_index + i]);
    }
}

/* sum_i^n a[i]*b[i] */
inline void BN_vec_inner_product(BIGNUM*& result, vector<BIGNUM*>& vec_a, vector<BIGNUM*>& vec_b)
{
    BN_zero(result); // set result = 0

    BIGNUM* product = BN_new();

    if (vec_a.size() != vec_b.size())
    {
        cout << "vector size does not match!" << endl;
        exit(EXIT_FAILURE);
    }
    for (auto i = 0; i < vec_a.size(); i++)
    {
        BN_mul(product, vec_a[i], vec_b[i], bn_ctx); // product = (vec_a[i]*vec_b[i]) mod order
        BN_add(result, result, product);     // result = (result+product) mod order
    }
    BN_mod(result, result, order, bn_ctx);
}

/* g[i] = g[i]+h[i] */
inline void ECP_vec_add(vector<EC_POINT*>& result, vector<EC_POINT*>& vec_A, vector<EC_POINT*>& vec_B)
{
    if (vec_A.size() != vec_B.size())
    {
        cout << "vector size does not match!" << endl;
        exit(EXIT_FAILURE);
    }
    for (auto i = 0; i < vec_A.size(); i++)
    {
        EC_POINT_add(group, result[i], vec_A[i], vec_B[i], bn_ctx);
    }
}

/* a[i] = (a[i]+b[i]) mod order */
inline void BN_vec_add(vector<BIGNUM*>& result, vector<BIGNUM*>& vec_a, vector<BIGNUM*>& vec_b)
{
    if (vec_a.size() != vec_b.size())
    {
        cout << "vector size does not match!" << endl;
        exit(EXIT_FAILURE);
    }
    for (auto i = 0; i < vec_a.size(); i++)
    {
        BN_mod_add(result[i], vec_a[i], vec_b[i], order, bn_ctx);
    }
}

/* a[i] = (a[i]-b[i]) mod order */
inline void BN_vec_sub(vector<BIGNUM*>& result, vector<BIGNUM*>& vec_a, vector<BIGNUM*>& vec_b)
{
    if (vec_a.size() != vec_b.size())
    {
        cout << "vector size does not match!" << endl;
        exit(EXIT_FAILURE);
    }
    for (auto i = 0; i < vec_a.size(); i++)
    {
        BN_mod_sub(result[i], vec_a[i], vec_b[i], order, bn_ctx);
    }
}

/* c[i] = a[i]*b[i] mod order */
inline void BN_vec_product(vector<BIGNUM*>& result, vector<BIGNUM*>& vec_a, vector<BIGNUM*>& vec_b)
{
    if (vec_a.size() != vec_b.size())
    {
        cout << "vector size does not match!" << endl;
        exit(EXIT_FAILURE);
    }

    for (auto i = 0; i < vec_a.size(); i++)
    {
        BN_mod_mul(result[i], vec_a[i], vec_b[i], order, bn_ctx); // product = (vec_a[i]*vec_b[i]) mod order
    }
}

/* compute the inverse of a[i] */
inline void BN_vec_inverse(vector<BIGNUM*>& vec_a_inverse, vector<BIGNUM*>& vec_a)
{
    for (auto i = 0; i < vec_a.size(); i++)
    {
        BN_mod_inverse(vec_a_inverse[i], vec_a[i], order, bn_ctx);
    }
}

/* vec_g = c * vec_g */
inline void ECP_vec_scalar(vector<EC_POINT*>& result, vector<EC_POINT*>& vec_A, BIGNUM*& c)
{
    for (auto i = 0; i < vec_A.size(); i++)
    {
        EC_POINT_mul(group, result[i], NULL, vec_A[i], c, bn_ctx); // result[i] = vec_g[i]^c
    }
}

/* result[i] = c * a[i] */
inline void BN_vec_scalar(vector<BIGNUM*>& result, vector<BIGNUM*>& vec_a, BIGNUM*& c)
{
    for (auto i = 0; i < vec_a.size(); i++)
    {
        BN_mod_mul(result[i], vec_a[i], c, order, bn_ctx);
    }
}

/* result[i] = -result[i] */
inline void BN_vec_negative(vector<BIGNUM*>& result)
{
    for (auto i = 0; i < result.size(); i++)
    {
        BN_mod_negative(result[i]);
    }
}

/* result[i] = A[i]*a[i] */
inline void ECP_vec_product(vector<EC_POINT*>& result, vector<EC_POINT*>& vec_A, vector<BIGNUM*>& vec_a)
{
    if (vec_A.size() != vec_a.size())
    {
        cout << "vector size does not match!" << endl;
        exit(EXIT_FAILURE);
    }
    for (auto i = 0; i < vec_A.size(); i++)
    {
        EC_POINT_mul(group, result[i], NULL, vec_A[i], vec_a[i], bn_ctx);
    }
}

/* result = sum_{i=1^n} a[i]*A[i] */
inline void ECP_vec_mul(EC_POINT*& result, vector<EC_POINT*>& vec_A, vector<BIGNUM*>& vec_a)
{
    if (vec_A.size() != vec_a.size()) {
        cout << "vector size does not match!" << endl;
        exit(EXIT_FAILURE);
    }
    EC_POINTs_mul(group, result, NULL, vec_A.size(),
        (const EC_POINT**)vec_A.data(), (const BIGNUM**)vec_a.data(), bn_ctx);
}

/* this module is used to enable fast verification (cf pp.15) */
void compute_vec_ss(vector<BIGNUM*>& vec_s, vector<BIGNUM*>& vec_x, vector<BIGNUM*>& vec_x_inverse)
{
    size_t m = vec_x.size();
    size_t n = vec_s.size(); //int n = pow(2, m); 

    // compute s[0], ..., s[i-1]
    // vector<BIGNUM *> vec_s(n); 
    uint64_t flag;
    for (auto i = 0; i < n; i++)
    {
        BN_one(vec_s[i]); // set s[i] = 1
        for (auto j = 0; j < m; j++)
        {
            flag = int_parse_binary(i, j, m);
            if (flag == 1) {
                BN_mod_mul(vec_s[i], vec_s[i], vec_x[j], order, bn_ctx);
            }
            else {
                BN_mod_mul(vec_s[i], vec_s[i], vec_x_inverse[j], order, bn_ctx);
            }
        }
    }
}


/* (Protocol 2 on pp.15) */
void InnerProduct_Setup(InnerProduct_KP& kp, size_t VECTOR_LEN, bool INITIAL_FLAG)
{
    kp.VECTOR_LEN = VECTOR_LEN;
    kp.LOG_VECTOR_LEN = log2(VECTOR_LEN);

    if (INITIAL_FLAG == true) {
        ECP_vec_random(kp.vec_G);
        ECP_vec_random(kp.vec_U);
        ECP_vec_random(kp.vec_E);
    }
}

void InnerProduct_Prove(InnerProduct_KP pp,
    InnerProduct_Witness witness,
    InnerProduct_Proof& proof)
{
    
    size_t n = pp.VECTOR_LEN;

    if (n == 1)
    {
        BN_copy(proof.m, witness.vec_m[0]);

#ifdef DEBUG
        cout << "Inner Product Proof Generation Finishes >>>" << endl;
#endif 

        return;
    }
    else
    {
        n = n / 2;

        vector<BIGNUM*> vec_xL(n), vec_xR(n);
        BN_vec_new(vec_xL);
        BN_vec_new(vec_xR);

        vector<EC_POINT*> vec_GL(n), vec_GR(n), vec_UL(n), vec_UR(n), vec_EL(n), vec_ER(n);
        ECP_vec_new(vec_GL);
        ECP_vec_new(vec_GR);
        ECP_vec_new(vec_UL);
        ECP_vec_new(vec_UR);
        ECP_vec_new(vec_EL);
        ECP_vec_new(vec_ER);

        // prepare aL, aR, bL, bR
        BN_vec_assign(vec_xL, witness.vec_m, "left");
        BN_vec_assign(vec_xR, witness.vec_m, "right");

        ECP_vec_assign(vec_GL, pp.vec_G, "left");
        ECP_vec_assign(vec_GR, pp.vec_G, "right");
        ECP_vec_assign(vec_UL, pp.vec_U, "left");
        ECP_vec_assign(vec_UR, pp.vec_U, "right");
        ECP_vec_assign(vec_EL, pp.vec_E, "left");
        ECP_vec_assign(vec_ER, pp.vec_E, "right");


        EC_POINT* LG = EC_POINT_new(group);
        EC_POINT* LU = EC_POINT_new(group);
        EC_POINT* LE = EC_POINT_new(group);
        ECP_vec_mul(LG, vec_GR, vec_xL);
        ECP_vec_mul(LU, vec_UR, vec_xL);
        ECP_vec_mul(LE, vec_ER, vec_xL);

        EC_POINT* RG = EC_POINT_new(group);
        EC_POINT* RU = EC_POINT_new(group);
        EC_POINT* RE = EC_POINT_new(group);
        ECP_vec_mul(RG, vec_GL, vec_xR);
        ECP_vec_mul(RU, vec_UL, vec_xR);
        ECP_vec_mul(RE, vec_EL, vec_xR);

        proof.vec_LA.push_back(LG);
        proof.vec_LU.push_back(LU);
        proof.vec_LE.push_back(LE);
        proof.vec_RA.push_back(RG);
        proof.vec_RU.push_back(RU);
        proof.vec_RE.push_back(RE);

        string transcript_str;
        transcript_str = ECP_ep2string(LG) + ECP_ep2string(RG);
        BIGNUM* c = BN_new();
        Hash_String_to_BN(transcript_str, c); // compute the n-th round challenge Eq (26,27)
        BIGNUM* c_inverse = BN_new();
        BN_mod_inverse(c_inverse, c, order, bn_ctx);

        // generate new pp
        InnerProduct_KP kp_sub;
        // pp_sub.VECTOR_LEN = pp.VECTOR_LEN/2; 
        // pp_sub.LOG_VECTOR_LEN = pp.LOG_VECTOR_LEN - 1; 
        InnerProduct_KP_new(kp_sub, pp.VECTOR_LEN / 2);
        InnerProduct_Setup(kp_sub, pp.VECTOR_LEN / 2, false);

        // compute vec_g
        ECP_vec_scalar(vec_GR, vec_GR, c);
        ECP_vec_add(kp_sub.vec_G, vec_GL, vec_GR); // Eq (29)
        ECP_vec_scalar(vec_UR, vec_UR, c);
        ECP_vec_add(kp_sub.vec_U, vec_UL, vec_UR); // Eq (29)
        ECP_vec_scalar(vec_ER, vec_ER, c);
        ECP_vec_add(kp_sub.vec_E, vec_EL, vec_ER); // Eq (29)
        
        InnerProduct_Witness witness_sub;
        InnerProduct_Witness_new(witness_sub, kp_sub.VECTOR_LEN);
        BN_vec_scalar(vec_xR, vec_xR, c_inverse);
        BN_vec_add(witness_sub.vec_m, vec_xL, vec_xR); // Eq (33)

        InnerProduct_Prove(kp_sub, witness_sub, proof);

        InnerProduct_KP_free(kp_sub);
        InnerProduct_Witness_free(witness_sub);

        BN_free(c), BN_free(c_inverse);

        BN_vec_free(vec_xL);
        BN_vec_free(vec_xR);

        ECP_vec_free(vec_GL);
        ECP_vec_free(vec_GR);
        ECP_vec_free(vec_UL);
        ECP_vec_free(vec_UR);
        ECP_vec_free(vec_EL);
        ECP_vec_free(vec_ER);
    }
}

void Samemultiscalar_Prove(InnerProduct_KP& kp,
    InnerProduct_Witness& witness,
    InnerProduct_Proof& proof,
    EC_POINT*& A,
    EC_POINT*& ZU,
    EC_POINT*& ZE,
    int i,
    long long &total
    ) {

    size_t n = kp.VECTOR_LEN;

    vector<BIGNUM*> vec_m(n);
    BN_vec_new(vec_m);
    BN_vec_random(vec_m);
    ECP_vec_mul(A, kp.vec_G, vec_m);
    ECP_vec_mul(ZU, kp.vec_U, vec_m);
    ECP_vec_mul(ZE, kp.vec_E, vec_m);

    

    vector<BIGNUM*> vec_mr(n);
    BN_vec_new(vec_mr);
    BN_vec_random(vec_mr);

    auto begin_prover = std::chrono::high_resolution_clock::now();
    ECP_vec_mul(proof.BA, kp.vec_G, vec_mr);
    ECP_vec_mul(proof.BU, kp.vec_U, vec_mr);
    ECP_vec_mul(proof.BE, kp.vec_E, vec_mr);

    string transcript_str;
    transcript_str = ECP_ep2string(A) + ECP_ep2string(ZU) + ECP_ep2string(ZE);
    BIGNUM* c = BN_new();
    Hash_String_to_BN(transcript_str, c);

    BN_vec_scalar(witness.vec_m, vec_m, c);
    BN_vec_add(witness.vec_m, vec_mr, witness.vec_m); // Eq (33)
    InnerProduct_Prove(kp, witness, proof);

}

//void Samemultiscalar_Prove(InnerProduct_KP& kp,
//    InnerProduct_Witness& witness,
//    InnerProduct_Proof& proof,
//    EC_POINT*& A,
//    EC_POINT*& ZU,
//    EC_POINT*& ZE, 
//    EC_POINT*& SA,
//    EC_POINT*& SZU,
//    EC_POINT*& SZE) {
//
//    size_t n = kp.VECTOR_LEN;
//
//    vector<BIGNUM*> vec_m(n);
//    BN_vec_new(vec_m);
//    BN_vec_random(vec_m);
//    ECP_vec_mul(A, kp.vec_G, vec_m);
//    ECP_vec_mul(ZU, kp.vec_U, vec_m);
//    ECP_vec_mul(ZE, kp.vec_E, vec_m);
//
//    vector<BIGNUM*> vec_ms(n);
//    BN_vec_new(vec_ms);
//    BN_vec_copy(vec_ms, vec_m);
//    vector<EC_POINT*> vec_GS;
//    vector<EC_POINT*> vec_US;
//    vector<EC_POINT*> vec_ES;
//    vec_GS.resize(kp.VECTOR_LEN); ECP_vec_new(vec_GS);
//    vec_US.resize(kp.VECTOR_LEN); ECP_vec_new(vec_US);
//    vec_ES.resize(kp.VECTOR_LEN); ECP_vec_new(vec_ES);
//    ECP_vec_copy(vec_GS, kp.vec_G);
//    ECP_vec_copy(vec_US, kp.vec_U);
//    ECP_vec_copy(vec_ES, kp.vec_E);
//    std::srand(1);
//    std::random_shuffle(vec_ms.begin(), vec_ms.end());
//    std::srand(1);
//    std::random_shuffle(vec_GS.begin(), vec_GS.end());
//    std::srand(1);
//    std::random_shuffle(vec_US.begin(), vec_US.end());
//    std::srand(1);
//    std::random_shuffle(vec_ES.begin(), vec_ES.end());
//    ECP_vec_mul(SA, vec_GS, vec_ms);
//    ECP_vec_mul(SZU, vec_US, vec_ms);
//    ECP_vec_mul(SZE, vec_ES, vec_ms);
//    
//    std::ofstream       prover_time;
//    std::stringstream   ptime;
//
//    auto begin_prover = std::chrono::high_resolution_clock::now();
//    vector<BIGNUM*> vec_mr(n);
//    BN_vec_new(vec_mr);
//    BN_vec_random(vec_mr);
//    ECP_vec_mul(proof.BA, kp.vec_G, vec_mr);
//    ECP_vec_mul(proof.BU, kp.vec_U, vec_mr);
//    ECP_vec_mul(proof.BE, kp.vec_E, vec_mr);
//
//    string transcript_str;
//    transcript_str = ECP_ep2string(A) + ECP_ep2string(ZU)+ECP_ep2string(ZE);
//    BIGNUM* c = BN_new();
//    Hash_String_to_BN(transcript_str, c);
//
//    BN_vec_scalar(witness.vec_m, vec_m, c);
//    BN_vec_add(witness.vec_m, vec_mr, witness.vec_m); // Eq (33)
//    InnerProduct_Prove(kp, witness, proof);
//    auto end_prover = std::chrono::high_resolution_clock::now();
//    auto elapsed_prover = std::chrono::duration_cast<std::chrono::nanoseconds>(end_prover - begin_prover);
//    prover_time.open("samescalar_prove_time.txt", std::ofstream::out | std::ofstream::app);
//    ptime << "samescalar_prove = " << elapsed_prover.count() * 1e-6 << " msec" << std::endl;
//    prover_time << ptime.str();
//    prover_time.close();
//}

/* Check if PI is a valid proof for inner product statement (G1^w = H1 and G2^w = H2) */
bool Samemultiscalar_Verify(InnerProduct_KP& kp,
    InnerProduct_Proof& proof,
    EC_POINT*& A,
    EC_POINT*& ZU,
    EC_POINT*& ZE)
{
    bool Validity;

    // recover the challenge
    vector<BIGNUM*> vec_x(kp.LOG_VECTOR_LEN); // the vector of challenge 
    vector<BIGNUM*> vec_x_inverse(kp.LOG_VECTOR_LEN); // the vector of challenge inverse
    size_t n = kp.VECTOR_LEN;

    BN_vec_new(vec_x);
    BN_vec_new(vec_x_inverse);

    string transcript_str;
    transcript_str = ECP_ep2string(A) + ECP_ep2string(ZU) + ECP_ep2string(ZE);
    BIGNUM* c = BN_new();
    Hash_String_to_BN(transcript_str, c);

    EC_POINT* VA = EC_POINT_new(group);
    EC_POINT* VZU = EC_POINT_new(group);
    EC_POINT* VZE = EC_POINT_new(group);
    EC_POINT_mul(group, VA, NULL, A, c, bn_ctx);
    EC_POINT_add(group, VA, VA, proof.BA,bn_ctx);
    EC_POINT_mul(group, VZU, NULL, ZU, c, bn_ctx);
    EC_POINT_add(group, VZU, VZU, proof.BU, bn_ctx);
    EC_POINT_mul(group, VZE, NULL, ZE, c, bn_ctx);
    EC_POINT_add(group, VZE, VZE, proof.BE, bn_ctx);

    EC_POINT* L = EC_POINT_new(group);
    EC_POINT* R = EC_POINT_new(group);

    vector<EC_POINT*> vec_G_cp(n), vec_U_cp(n), vec_E_cp(n);
    ECP_vec_new(vec_G_cp);
    ECP_vec_new(vec_U_cp);
    ECP_vec_new(vec_E_cp);
    ECP_vec_copy(vec_G_cp, kp.vec_G);
    ECP_vec_copy(vec_U_cp, kp.vec_U);
    ECP_vec_copy(vec_E_cp, kp.vec_E);


    for (auto i = 0; i < kp.LOG_VECTOR_LEN; i++)
    {
        transcript_str = ECP_ep2string(proof.vec_LA[i]) + ECP_ep2string(proof.vec_RA[i]);
        Hash_String_to_BN(transcript_str, vec_x[i]); // reconstruct the challenge
        BN_mod_inverse(vec_x_inverse[i], vec_x[i], order, bn_ctx);

        EC_POINT_mul(group, L, NULL, proof.vec_LA[i], vec_x[i], bn_ctx);
        EC_POINT_mul(group, R, NULL, proof.vec_RA[i], vec_x_inverse[i], bn_ctx);
        EC_POINT_add(group, VA, L, VA, bn_ctx);
        EC_POINT_add(group, VA, R, VA, bn_ctx);

        EC_POINT_mul(group, L, NULL, proof.vec_LU[i], vec_x[i], bn_ctx);
        EC_POINT_mul(group, R, NULL, proof.vec_RU[i], vec_x_inverse[i], bn_ctx);
        EC_POINT_add(group, VZU, L, VZU, bn_ctx);
        EC_POINT_add(group, VZU, R, VZU, bn_ctx);

        EC_POINT_mul(group, L, NULL, proof.vec_LE[i], vec_x[i], bn_ctx);
        EC_POINT_mul(group, R, NULL, proof.vec_RE[i], vec_x_inverse[i], bn_ctx);
        EC_POINT_add(group, VZE, L, VZE, bn_ctx);
        EC_POINT_add(group, VZE, R, VZE, bn_ctx);

        n = n / 2;
        vector<EC_POINT*> vec_GL(n), vec_GR(n), vec_UL(n), vec_UR(n), vec_EL(n), vec_ER(n);
        ECP_vec_new(vec_GL);
        ECP_vec_new(vec_GR);
        ECP_vec_new(vec_UL);
        ECP_vec_new(vec_UR);
        ECP_vec_new(vec_EL);
        ECP_vec_new(vec_ER);

        ECP_vec_assign(vec_GL, vec_G_cp, "left");
        ECP_vec_assign(vec_GR, vec_G_cp, "right");
        ECP_vec_assign(vec_UL, vec_U_cp, "left");
        ECP_vec_assign(vec_UR, vec_U_cp, "right");
        ECP_vec_assign(vec_EL, vec_E_cp, "left");
        ECP_vec_assign(vec_ER, vec_E_cp, "right");
        vec_G_cp.resize(n);
        vec_U_cp.resize(n);
        vec_E_cp.resize(n);
        ECP_vec_scalar(vec_GR, vec_GR, vec_x[i]);
        ECP_vec_add(vec_G_cp, vec_GL, vec_GR); // Eq (29)
        ECP_vec_scalar(vec_UR, vec_UR, vec_x[i]);
        ECP_vec_add(vec_U_cp, vec_UL, vec_UR); // Eq (29)
        ECP_vec_scalar(vec_ER, vec_ER, vec_x[i]);
        ECP_vec_add(vec_E_cp, vec_EL, vec_ER); // Eq (29)
        ECP_vec_free(vec_GL);
        ECP_vec_free(vec_GR);
        ECP_vec_free(vec_UL);
        ECP_vec_free(vec_UR);
        ECP_vec_free(vec_EL);
        ECP_vec_free(vec_ER);
        
        
    }
    bool V1, V2, V3;
    
    EC_POINT* G1;
    EC_POINT* U1;
    EC_POINT* E1;
    G1 = vec_G_cp[0];
    U1 = vec_U_cp[0];
    E1 = vec_E_cp[0];
    EC_POINT_mul(group, G1, NULL, G1, proof.m, bn_ctx);
    EC_POINT_mul(group, U1, NULL, U1, proof.m, bn_ctx);
    EC_POINT_mul(group, E1, NULL, E1, proof.m, bn_ctx);
    V1 = (EC_POINT_cmp(group, VA, G1, bn_ctx) == 0);
    V2 = (EC_POINT_cmp(group, VZU, U1, bn_ctx) == 0);
    V3 = (EC_POINT_cmp(group, VZE, E1, bn_ctx) == 0);
    //V3 = (EC_POINT_cmp(group, VZE, U1, bn_ctx) == 0);
    Validity = V1 && V2 && V3;
    return Validity;
}

void Sigma_InnerProduct_Prove(InnerProduct_KP pp,
    InnerProduct_Witness witness,
    InnerProduct_Proof& proof,
    EC_POINT*& hA,
    EC_POINT*& hU,
    EC_POINT*& hE,
    BIGNUM*& proofr)
{

    size_t n = pp.VECTOR_LEN;

    if (n == 1)
    {
        BIGNUM* xr = BN_new();
        BN_random(xr);
        BIGNUM* rr = BN_new();
        BN_random(rr);
        EC_POINT* S = EC_POINT_new(group);

        string transcript_str;
        transcript_str = ECP_ep2string(proof.BA) + ECP_ep2string(proof.BU) + ECP_ep2string(proof.BE);
        BIGNUM* c = BN_new();
        Hash_String_to_BN(transcript_str, c);

        EC_POINT_mul(group, proof.BA, NULL, pp.vec_G[0], xr, bn_ctx);
        EC_POINT_mul(group, S, NULL, hA, rr, bn_ctx);
        EC_POINT_add(group, proof.BA, S, proof.BA, bn_ctx);
        EC_POINT_mul(group, proof.BU, NULL, pp.vec_U[0], xr, bn_ctx);
        EC_POINT_mul(group, S, NULL, hU, rr, bn_ctx);
        EC_POINT_add(group, proof.BU, S, proof.BU, bn_ctx);
        EC_POINT_mul(group, proof.BE, NULL, pp.vec_E[0], xr, bn_ctx);
        EC_POINT_mul(group, S, NULL, hE, rr, bn_ctx);
        EC_POINT_add(group, proof.BE, S, proof.BE, bn_ctx);


        BIGNUM* x = BN_new();
        BN_mod_mul(x, c, witness.vec_m[0], order, bn_ctx);
        BN_mod_add(xr, xr, x, order, bn_ctx);
        BN_copy(proof.m, xr);

        BN_mod_mul(proofr, c, proofr, order, bn_ctx);
        BN_mod_add(proofr, rr, proofr, order, bn_ctx);


#ifdef DEBUG
        cout << "Inner Product Proof Generation Finishes >>>" << endl;
#endif 

        return;
    }
    else
    {

        n = n / 2;

        vector<BIGNUM*> vec_xL(n), vec_xR(n);
        BN_vec_new(vec_xL);
        BN_vec_new(vec_xR);
        BIGNUM* rL = BN_new();
        BIGNUM* rR = BN_new();
        BN_random(rL);
        BN_random(rR);

        vector<EC_POINT*> vec_GL(n), vec_GR(n), vec_UL(n), vec_UR(n), vec_EL(n), vec_ER(n);
        ECP_vec_new(vec_GL);
        ECP_vec_new(vec_GR);
        ECP_vec_new(vec_UL);
        ECP_vec_new(vec_UR);
        ECP_vec_new(vec_EL);
        ECP_vec_new(vec_ER);

        // prepare aL, aR, bL, bR
        BN_vec_assign(vec_xL, witness.vec_m, "left");
        BN_vec_assign(vec_xR, witness.vec_m, "right");

        ECP_vec_assign(vec_GL, pp.vec_G, "left");
        ECP_vec_assign(vec_GR, pp.vec_G, "right");
        ECP_vec_assign(vec_UL, pp.vec_U, "left");
        ECP_vec_assign(vec_UR, pp.vec_U, "right");
        ECP_vec_assign(vec_EL, pp.vec_E, "left");
        ECP_vec_assign(vec_ER, pp.vec_E, "right");


        EC_POINT* LG = EC_POINT_new(group);
        EC_POINT* LU = EC_POINT_new(group);
        EC_POINT* LE = EC_POINT_new(group);
        ECP_vec_mul(LG, vec_GR, vec_xL);
        ECP_vec_mul(LU, vec_UR, vec_xL);
        ECP_vec_mul(LE, vec_ER, vec_xL);
        EC_POINT* SLR = EC_POINT_new(group);
        EC_POINT_mul(group, SLR, NULL, hA, rL, bn_ctx);
        EC_POINT_add(group, LG, SLR, LG, bn_ctx);
        EC_POINT_mul(group, SLR, NULL, hU, rL, bn_ctx);
        EC_POINT_add(group, LU, SLR, LU, bn_ctx);
        EC_POINT_mul(group, SLR, NULL, hE, rL, bn_ctx);
        EC_POINT_add(group, LE, SLR, LE, bn_ctx);

        EC_POINT* RG = EC_POINT_new(group);
        EC_POINT* RU = EC_POINT_new(group);
        EC_POINT* RE = EC_POINT_new(group);
        ECP_vec_mul(RG, vec_GL, vec_xR);
        ECP_vec_mul(RU, vec_UL, vec_xR);
        ECP_vec_mul(RE, vec_EL, vec_xR);
        EC_POINT_mul(group, SLR, NULL, hA, rR, bn_ctx);
        EC_POINT_add(group, RG, SLR, RG, bn_ctx);
        EC_POINT_mul(group, SLR, NULL, hU, rR, bn_ctx);
        EC_POINT_add(group, RU, SLR, RU, bn_ctx);
        EC_POINT_mul(group, SLR, NULL, hE, rR, bn_ctx);
        EC_POINT_add(group, RE, SLR, RE, bn_ctx);

        proof.vec_LA.push_back(LG);
        proof.vec_LU.push_back(LU);
        proof.vec_LE.push_back(LE);
        proof.vec_RA.push_back(RG);
        proof.vec_RU.push_back(RU);
        proof.vec_RE.push_back(RE);

        string transcript_str;
        transcript_str = ECP_ep2string(LG) + ECP_ep2string(RG);
        BIGNUM* c = BN_new();
        Hash_String_to_BN(transcript_str, c); // compute the n-th round challenge Eq (26,27)
        BIGNUM* c_inverse = BN_new();
        BN_mod_inverse(c_inverse, c, order, bn_ctx);

        // generate new pp
        InnerProduct_KP kp_sub;
        // pp_sub.VECTOR_LEN = pp.VECTOR_LEN/2; 
        // pp_sub.LOG_VECTOR_LEN = pp.LOG_VECTOR_LEN - 1; 
        InnerProduct_KP_new(kp_sub, pp.VECTOR_LEN / 2);
        InnerProduct_Setup(kp_sub, pp.VECTOR_LEN / 2, false);

        // compute vec_g
        ECP_vec_scalar(vec_GR, vec_GR, c);
        ECP_vec_add(kp_sub.vec_G, vec_GL, vec_GR); // Eq (29)
        ECP_vec_scalar(vec_UR, vec_UR, c);
        ECP_vec_add(kp_sub.vec_U, vec_UL, vec_UR); // Eq (29)
        ECP_vec_scalar(vec_ER, vec_ER, c);
        ECP_vec_add(kp_sub.vec_E, vec_EL, vec_ER); // Eq (29)

        InnerProduct_Witness witness_sub;
        InnerProduct_Witness_new(witness_sub, kp_sub.VECTOR_LEN);
        BN_vec_scalar(vec_xR, vec_xR, c_inverse);
        BN_vec_add(witness_sub.vec_m, vec_xL, vec_xR); // Eq (33)

        BIGNUM* temp = BN_new();
        BN_mod_mul(temp, rL, c, order, bn_ctx);
        BN_mod_add(proofr, temp, proofr, order, bn_ctx);
        BN_mod_mul(temp, rR, c_inverse, order, bn_ctx);
        BN_mod_add(proofr, temp, proofr, order, bn_ctx);

        Sigma_InnerProduct_Prove(kp_sub, witness_sub, proof, hA, hU, hE, proofr);

        InnerProduct_KP_free(kp_sub);
        InnerProduct_Witness_free(witness_sub);

        BN_free(c), BN_free(c_inverse);

        BN_vec_free(vec_xL);
        BN_vec_free(vec_xR);

        ECP_vec_free(vec_GL);
        ECP_vec_free(vec_GR);
        ECP_vec_free(vec_UL);
        ECP_vec_free(vec_UR);
        ECP_vec_free(vec_EL);
        ECP_vec_free(vec_ER);
        EC_POINT_free(SLR);
    }
}


void Sigma_Samemultiscalar_Prove(InnerProduct_KP& kp,
    InnerProduct_Witness& witness,
    InnerProduct_Proof& proof,
    EC_POINT*& A,
    EC_POINT*& ZU,
    EC_POINT*& ZE,
    EC_POINT*& hA,
    EC_POINT*& hU,
    EC_POINT*& hE,
    BIGNUM*& proofr,
    int i,
    long long& total
) {

    size_t n = kp.VECTOR_LEN;

    vector<BIGNUM*> vec_m(n);
    BN_vec_new(vec_m);
    BN_vec_random(vec_m);
    ECP_vec_mul(A, kp.vec_G, vec_m);
    ECP_vec_mul(ZU, kp.vec_U, vec_m);
    ECP_vec_mul(ZE, kp.vec_E, vec_m);

    EC_POINT* SA = EC_POINT_new(group);
    EC_POINT* SZU = EC_POINT_new(group);
    EC_POINT* SZE = EC_POINT_new(group);
    EC_POINT_mul(group, SA, NULL, hA, proofr, bn_ctx);
    EC_POINT_mul(group, SZU, NULL, hU, proofr, bn_ctx);
    EC_POINT_mul(group, SZE, NULL, hE, proofr, bn_ctx);
    EC_POINT_add(group, proof.BA, SA, A, bn_ctx);
    EC_POINT_add(group, proof.BU, SZU, ZU, bn_ctx);
    EC_POINT_add(group, proof.BE, SZE, ZE, bn_ctx);
    EC_POINT_add(group, A, SA, A, bn_ctx);
    EC_POINT_add(group, ZU, SZU, ZU, bn_ctx);
    EC_POINT_add(group, ZE, SZE, ZE, bn_ctx);
 
    BN_vec_copy(witness.vec_m, vec_m);

    Sigma_InnerProduct_Prove(kp, witness, proof,hA,hU,hE,proofr);
}


bool Sigma_Samemultiscalar_Verify(InnerProduct_KP& kp,
    InnerProduct_Proof& proof,
    EC_POINT*& A,
    EC_POINT*& ZU,
    EC_POINT*& ZE,
    EC_POINT*& hA,
    EC_POINT*& hU,
    EC_POINT*& hE,
    BIGNUM*& proofr
    )
{
    bool Validity;

    // recover the challenge
    vector<BIGNUM*> vec_x(kp.LOG_VECTOR_LEN); // the vector of challenge 
    vector<BIGNUM*> vec_x_inverse(kp.LOG_VECTOR_LEN); // the vector of challenge inverse
    size_t n = kp.VECTOR_LEN;

    BN_vec_new(vec_x);
    BN_vec_new(vec_x_inverse);

    string transcript_str;
    transcript_str = ECP_ep2string(A) + ECP_ep2string(ZU) + ECP_ep2string(ZE);
    BIGNUM* c = BN_new();
    Hash_String_to_BN(transcript_str, c);

    EC_POINT* VA = EC_POINT_new(group);
    EC_POINT* VZU = EC_POINT_new(group);
    EC_POINT* VZE = EC_POINT_new(group);
    

    EC_POINT* L = EC_POINT_new(group);
    EC_POINT* R = EC_POINT_new(group);

    vector<EC_POINT*> vec_G_cp(n), vec_U_cp(n), vec_E_cp(n);
    ECP_vec_new(vec_G_cp);
    ECP_vec_new(vec_U_cp);
    ECP_vec_new(vec_E_cp);
    ECP_vec_copy(vec_G_cp, kp.vec_G);
    ECP_vec_copy(vec_U_cp, kp.vec_U);
    ECP_vec_copy(vec_E_cp, kp.vec_E);


    for (auto i = 0; i < kp.LOG_VECTOR_LEN; i++)
    {
        transcript_str = ECP_ep2string(proof.vec_LA[i]) + ECP_ep2string(proof.vec_RA[i]);
        Hash_String_to_BN(transcript_str, vec_x[i]); // reconstruct the challenge
        BN_mod_inverse(vec_x_inverse[i], vec_x[i], order, bn_ctx);

        EC_POINT_mul(group, L, NULL, proof.vec_LA[i], vec_x[i], bn_ctx);
        EC_POINT_mul(group, R, NULL, proof.vec_RA[i], vec_x_inverse[i], bn_ctx);
        EC_POINT_add(group, A, L, A, bn_ctx);
        EC_POINT_add(group, A, R, A, bn_ctx);

        EC_POINT_mul(group, L, NULL, proof.vec_LU[i], vec_x[i], bn_ctx);
        EC_POINT_mul(group, R, NULL, proof.vec_RU[i], vec_x_inverse[i], bn_ctx);
        EC_POINT_add(group, ZU, L, ZU, bn_ctx);
        EC_POINT_add(group, ZU, R, ZU, bn_ctx);

        EC_POINT_mul(group, L, NULL, proof.vec_LE[i], vec_x[i], bn_ctx);
        EC_POINT_mul(group, R, NULL, proof.vec_RE[i], vec_x_inverse[i], bn_ctx);
        EC_POINT_add(group, ZE, L, ZE, bn_ctx);
        EC_POINT_add(group, ZE, R, ZE, bn_ctx);

        n = n / 2;
        vector<EC_POINT*> vec_GL(n), vec_GR(n), vec_UL(n), vec_UR(n), vec_EL(n), vec_ER(n);
        ECP_vec_new(vec_GL);
        ECP_vec_new(vec_GR);
        ECP_vec_new(vec_UL);
        ECP_vec_new(vec_UR);
        ECP_vec_new(vec_EL);
        ECP_vec_new(vec_ER);

        ECP_vec_assign(vec_GL, vec_G_cp, "left");
        ECP_vec_assign(vec_GR, vec_G_cp, "right");
        ECP_vec_assign(vec_UL, vec_U_cp, "left");
        ECP_vec_assign(vec_UR, vec_U_cp, "right");
        ECP_vec_assign(vec_EL, vec_E_cp, "left");
        ECP_vec_assign(vec_ER, vec_E_cp, "right");
        vec_G_cp.resize(n);
        vec_U_cp.resize(n);
        vec_E_cp.resize(n);
        ECP_vec_scalar(vec_GR, vec_GR, vec_x[i]);
        ECP_vec_add(vec_G_cp, vec_GL, vec_GR); // Eq (29)
        ECP_vec_scalar(vec_UR, vec_UR, vec_x[i]);
        ECP_vec_add(vec_U_cp, vec_UL, vec_UR); // Eq (29)
        ECP_vec_scalar(vec_ER, vec_ER, vec_x[i]);
        ECP_vec_add(vec_E_cp, vec_EL, vec_ER); // Eq (29)
        ECP_vec_free(vec_GL);
        ECP_vec_free(vec_GR);
        ECP_vec_free(vec_UL);
        ECP_vec_free(vec_UR);
        ECP_vec_free(vec_EL);
        ECP_vec_free(vec_ER);


    }
    bool V1, V2, V3;

    EC_POINT* temp = EC_POINT_new(group);
   
    
    EC_POINT_mul(group, VA, NULL, vec_G_cp[0], proof.m, bn_ctx);
    EC_POINT_mul(group, temp, NULL, hA, proofr, bn_ctx);
    EC_POINT_add(group, VA, temp, VA, bn_ctx);

    EC_POINT_mul(group, VZU, NULL, vec_U_cp[0], proof.m, bn_ctx);
    EC_POINT_mul(group, temp, NULL, hU, proofr, bn_ctx);
    EC_POINT_add(group,VZU, temp, VZU, bn_ctx);

    EC_POINT_mul(group, VZE, NULL, vec_E_cp[0], proof.m, bn_ctx);
    EC_POINT_mul(group, temp, NULL, hE, proofr, bn_ctx);
    EC_POINT_add(group, VZE, temp, VZE, bn_ctx);
    
    EC_POINT* VAR = EC_POINT_new(group);
    EC_POINT* VZUR = EC_POINT_new(group);
    EC_POINT* VZER = EC_POINT_new(group);

    EC_POINT_mul(group, VAR, NULL, A, c, bn_ctx);
    EC_POINT_add(group, VAR, proof.BA, VAR, bn_ctx);

    EC_POINT_mul(group, VZUR, NULL, ZU, c, bn_ctx);
    EC_POINT_add(group, VZUR, proof.BU, VZUR, bn_ctx);

    EC_POINT_mul(group, VZER, NULL, ZE, c, bn_ctx);
    EC_POINT_add(group, VZER, proof.BE, VZER, bn_ctx);

    V1 = (EC_POINT_cmp(group, VA, VAR, bn_ctx) == 0);
    V2 = (EC_POINT_cmp(group, VZU, VZUR, bn_ctx) == 0);
    V3 = (EC_POINT_cmp(group, VZE, VZER, bn_ctx) == 0);
    //V3 = (EC_POINT_cmp(group, VZE, VZUR, bn_ctx) == 0);
    Validity = V1 && V2 && V3;
    return Validity;
}


#endif



