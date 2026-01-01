//#define DEBUG

#include "elgamal_pke.hpp"
#include "group_commit.hpp"
#include "samescalar.hpp"
#include "samemultiscalar.hpp"
#include <algorithm>
#include <random>

#define TEST_NUM 30000
//
//
void test_elgamal(
    size_t IO_THREAD_NUM, size_t DEC_THREAD_NUM,size_t VECTOR_LEN)
{
    cout << "begin the basic correctness test >>>" << endl;

    ElGamal_PP pp;
    ElGamal_PP_new(pp);

    ElGamal_Setup(pp,  IO_THREAD_NUM, DEC_THREAD_NUM);

    ElGamal_KP keypair;
    ElGamal_KP_new(keypair);

    ElGamal_KeyGen(pp, keypair);

    ElGamal_CT CT;
    ElGamal_CT_new(CT);

    /* random test */
    SplitLine_print('-');
    cout << "begin the random test >>>" << endl;
    EC_POINT* m = EC_POINT_new(group);
    ECP_random(m);
    ECP_print(m,"m");
    ElGamal_Enc(pp, keypair.pk, m, CT);
    EC_POINT* m_dec = EC_POINT_new(group);
    ElGamal_Dec(pp, keypair.sk, CT, m_dec);
    ECP_print(m_dec,"m_dec");



    /* shuffle test */
    /*vector<BIGNUM*> vec_m(VECTOR_LEN);
    BN_vec_new(vec_m);
    vector<BIGNUM*> vec_mr(VECTOR_LEN);
    BN_vec_new(vec_mr);
    BN_vec_random(vec_mr);
    BN_vec_copy(vec_m, vec_mr);
    BN_vec_print(vec_m,"vec_m");
    BN_vec_print(vec_mr,"vec_mr");
    cout << "begin shuffle" << endl;
    std::srand(1);
    std::random_shuffle(vec_m.begin(), vec_m.end());
    std::srand(1);
    std::random_shuffle(vec_mr.begin(), vec_mr.end());
    BN_vec_print(vec_m, "vec_m");
    BN_vec_print(vec_mr, "vec_mr");*/

    /*samescalar test*/
    // 使用Elgamal加密密文模拟生成群元素，零知识证明samescalar argument
    SplitLine_print('-');
    Group_Commit_GCM cmu;
    Group_Commit_GCM_new(cmu);
    Group_Commit_GCM cme;
    Group_Commit_GCM_new(cme);
    Group_Commit_GCM cmur;
    Group_Commit_GCM_new(cmur);
    Group_Commit_GCM cmer;
    Group_Commit_GCM_new(cmer);
    Group_Commit_KP kp;
    Group_Commit_KP_new(kp);
    Group_Commit_KP_KeyGen(kp);
    BIGNUM* k= BN_new();
    BN_random(k);
    BIGNUM* ru= BN_new();
    BN_random(ru);
    BIGNUM* re= BN_new();
    BN_random(re);
    BN_print(re);
    EC_POINT_mul(group, m, NULL, CT.X, k, bn_ctx);
    ECP_print(CT.X);
    Group_Commit(kp.gu, kp.h, m, ru, cmu);
    EC_POINT_mul(group, m, NULL, CT.Y, k, bn_ctx);
    Group_Commit(kp.ge, kp.h, m, re, cme);
    BN_print(re, "field element");
    EC_POINT* g = EC_POINT_new(group);
    ECP_random(g);
    ECP_print(g, "group element");
    ECP_print(CT.X, "Group element");
    BIGNUM* zk = BN_new();
    BIGNUM* zru = BN_new();
    BIGNUM* zre = BN_new();
    Prove_SameScalar(kp, CT.X, CT.Y, k, ru, re, zk, zru, zre, cmur, cmer);
    bool success=false;
    cout << "success= >>>" << success << endl;
    success = Verify_SameScalar(kp, CT.X, CT.Y, cmu, cme, cmur, cmer, zk, zru, zre);
    cout << "SameScalar verify is success? >>>" <<success<< endl;

    ElGamal_PP_free(pp);
    ElGamal_KP_free(keypair);
    ElGamal_CT_free(CT);
    Group_Commit_GCM_free(cmu);
    Group_Commit_GCM_free(cme);
    Group_Commit_GCM_free(cmur);
    Group_Commit_GCM_free(cmer);
    Group_Commit_KP_free(kp);
    EC_POINT_free(m);
    EC_POINT_free(m_dec);
    BN_free(k);
    BN_free(ru);
    BN_free(re);
    BN_free(zk);
    BN_free(zru);
    BN_free(zre);

    /*samemultiscalar test*/
    // 生成随机群元素向量、随机数向量、随机数掩码向量构造samemultiscalar argument
    // A=xG ZU=xU ZE=xE SA=shuffle(x)shuffle(G) SZU=shuffle(x)shuffle(U) SZE=shuffle(x)shuffle(E)
    // 证明向量x与向量G、U、E一一对应且进行内积运算
    SplitLine_print('-');
    
    long long total;

    for (int i = 0; i < 1; i++)
    {
        InnerProduct_KP ipkp;
        InnerProduct_KP_new(ipkp, VECTOR_LEN);
        InnerProduct_Witness witness;
        InnerProduct_Witness_new(witness, VECTOR_LEN);
        InnerProduct_Proof proof;
        InnerProduct_Proof_new(proof);

        InnerProduct_Setup(ipkp, VECTOR_LEN, true);
        EC_POINT* A = EC_POINT_new(group);
        EC_POINT* ZU = EC_POINT_new(group);
        EC_POINT* ZE = EC_POINT_new(group);

        Samemultiscalar_Prove(ipkp, witness, proof, A, ZU, ZE, i, total);
        ECP_print(A, "A element");
        success = false;
        //cout << "success= >>>" << success << endl;
        success = Samemultiscalar_Verify(ipkp, proof, A, ZU, ZE);
        if (success == false)
        {
            cout << "SameMultiScalar verify is failed >>>" << endl;
        }
        else
        {
            cout << "SameMultiScalar verify is success >>>" << endl;
        }
        InnerProduct_KP_free(ipkp);
        InnerProduct_Witness_free(witness);
        InnerProduct_Proof_free(proof);
        EC_POINT_free(A);
        EC_POINT_free(ZU);
        EC_POINT_free(ZE);
    }
    
   

    for (int i = 0; i < 1; i++)
    {
        InnerProduct_KP ipkp1;
        InnerProduct_KP_new(ipkp1, VECTOR_LEN);
        InnerProduct_Witness witness1;
        InnerProduct_Witness_new(witness1, VECTOR_LEN);
        InnerProduct_Proof proof1;
        InnerProduct_Proof_new(proof1);

        InnerProduct_Setup(ipkp1, VECTOR_LEN, true);
        EC_POINT* A1 = EC_POINT_new(group);
        EC_POINT* ZU1 = EC_POINT_new(group);
        EC_POINT* ZE1 = EC_POINT_new(group);
        EC_POINT* hA1 = EC_POINT_new(group);
        EC_POINT* hU1 = EC_POINT_new(group);
        EC_POINT* hE1 = EC_POINT_new(group);
        ECP_random(hA1);
        ECP_random(hU1);
        ECP_random(hE1);
        BIGNUM* proofr1 = BN_new();
        BN_random(proofr1);

        Sigma_Samemultiscalar_Prove(ipkp1, witness1, proof1, A1, ZU1, ZE1, hA1, hU1, hE1, proofr1, i, total);
        success = false;
        //cout << "success= >>>" << success << endl;
        success = Sigma_Samemultiscalar_Verify(ipkp1, proof1, A1, ZU1, ZE1, hA1, hU1, hE1, proofr1);
        if (success == false)
        {
            cout << "Sigma_SameMultiScalar verify is failed >>>" << endl;
        }
        else
        {
            cout << "Sigma_SameMultiScalar verify is success >>>" << endl;
        }
        InnerProduct_KP_free(ipkp1);
        InnerProduct_Witness_free(witness1);
        InnerProduct_Proof_free(proof1);
        EC_POINT_free(A1);
        EC_POINT_free(ZU1);
        EC_POINT_free(ZE1);
    }

}


int main()
{
    global_initialize(NID_X9_62_prime256v1);




    SplitLine_print('-');
    cout << "ElGamal PKE test begins >>>>>>" << endl;
    SplitLine_print('-');

     size_t IO_THREAD_NUM = 4; 
     size_t DEC_THREAD_NUM = 4;  
     size_t VECTOR_LEN = 2048;

    test_elgamal(IO_THREAD_NUM, DEC_THREAD_NUM,VECTOR_LEN);

     SplitLine_print('-'); 
     cout << "ElGamal PKE test finishes <<<<<<" << endl; 
     SplitLine_print('-'); 

    global_finalize();

    return 0;
}
