/* ****************************************************************************************************************************************************************** *\
 * @file    test.cc                                                                                                                                                   *
 * @author  Myungsun Kim                                                                                                                                              *
 * @date    Feb. 10, 2021                                                                                                                                             *
 * @brief   main                                                                                                                                                      *
 *                                                                                                                                                                    *
\* ****************************************************************************************************************************************************************** */
#include <iostream>

//#include "Camenisch-shoup.h"
#include "DOPRF.hpp"

using namespace NTL;

/**
 * @brief main
 *
 */
int
main(int argc, char** argv)
{

    global_initialize(NID_X9_62_prime256v1);
    EG_Public eg_pp;
    EG_Public_new(eg_pp);
    EG_Public_setup(eg_pp);

    //encryption_test();
    //return 1;

    long inputsize = 4;
    //inputsize = 16;
    // zeta = 2, s = 4
    Party P1, P2;
    long begin = 0;
    // zeta = 2, s = 4
    while (begin < inputsize) {
        vector<long> inputtemp;
        inputtemp.emplace_back(begin);
        inputtemp.emplace_back(begin+1);
        inputtemp.emplace_back(begin+2);
        inputtemp.emplace_back(begin+3);
        P1.input_set.emplace_back(inputtemp);
        P2.input_set.emplace_back(inputtemp);
        begin = begin + 4;
    }
    // zeta = 1, s = 1
    /*while (begin < inputsize) {
        vector<long> inputtemp;
        inputtemp.emplace_back(begin);
        P1.input_set.emplace_back(inputtemp);
        P2.input_set.emplace_back(inputtemp);
        begin = begin + 1;
    }*/
    Party_setup(eg_pp, P1);
    Party_setup(eg_pp, P2);
    std::cout << "party setup finish " << std::endl;
    Round_1_witness round_1_witness;
    Round_1_proof_random round_1_proof_random;
    Round_1_message round_1_message;
    round_1(P2, round_1_witness, round_1_proof_random, round_1_message);
    std::cout << "round_1_verify " << round_1_verify(P2,round_1_witness,round_1_proof_random,round_1_message) << std::endl;
    Round_2_witness round_2_witness;
    Round_2_proof_random round_2_proof_random;
    Round_2_message round_2_message;
    round_2(P1, P2, eg_pp, round_1_message, round_1_witness, round_2_witness, round_2_proof_random, round_2_message);
    std::cout << "round_2_verify " << round_2_verify(P1,P2,eg_pp,round_1_message,round_1_witness,round_2_witness,round_2_proof_random,round_2_message) << std::endl;
    
    Round_3_witness round_3_witness;
    Round_3_proof_random round_3_proof_random;
    Round_3_message round_3_message;
    round_3(eg_pp, P1, P2, round_2_message, round_3_witness, round_3_proof_random, round_3_message);
    std::cout << "round_3_verify " << round_3_verify(eg_pp, P1, P2, round_2_message, round_3_witness, round_3_proof_random, round_3_message) << std::endl;
    
    EC_POINT_free(round_3_message.random_beta_eg_com);
    BN_free(round_3_message.random_beta_eg_com_random);
    for (size_t i = 0; i < round_3_message.beta_eg_com.size(); ++i) {
        EC_POINT_free(round_3_message.beta_eg_com[i]);
        BN_free(round_3_witness.beta_eg_com_random[i]);
    }
    
    EG_Public_free(eg_pp);
    EG_Keypair_free(P1.eg_kp);
    EG_Keypair_free(P2.eg_kp);
    for (int i = 0; i < round_2_message.g_a.size(); ++i) {
        ECP_vec_free(round_2_message.g_a[i]);
    }
    ECP_vec_free(round_2_message.g_random_a);
    global_finalize();
    return 1;
}