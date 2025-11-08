#pragma once
#include "global.hpp"
#include "hash.hpp"
#include "print.hpp"
#include "routines.hpp"

struct Group_Commit_KP
{
    EC_POINT* gu;
    EC_POINT* ge;
    EC_POINT* h; 
};

struct Group_Commit_GCM
{
    EC_POINT* X; 
    EC_POINT* Y; 
};

void Group_Commit_KP_new(Group_Commit_KP& kp)
{
    kp.gu = EC_POINT_new(group);
    kp.ge = EC_POINT_new(group);
    kp.h = EC_POINT_new(group);
}

void Group_Commit_KP_free(Group_Commit_KP& kp)
{
    EC_POINT_free(kp.gu);
    EC_POINT_free(kp.ge);
    EC_POINT_free(kp.h);
}

void Group_Commit_KP_print(Group_Commit_KP& kp)
{
    ECP_print(kp.gu, "kp.gu");
    ECP_print(kp.ge, "kp.ge");
    ECP_print(kp.h, "kp.h");
}

void Group_Commit_GCM_new(Group_Commit_GCM& gcm)
{
    gcm.X = EC_POINT_new(group);
    gcm.Y = EC_POINT_new(group);
}

void Group_Commit_GCM_free(Group_Commit_GCM& gcm)
{
    EC_POINT_free(gcm.X);
    EC_POINT_free(gcm.Y);
}

void Group_Commit_GCM_print(Group_Commit_GCM& gcm)
{
    ECP_print(gcm.X, "gcm.X");
    ECP_print(gcm.Y, "gcm.Y");
}

void Group_Commit_KP_KeyGen(Group_Commit_KP& kp)
{
    ECP_random(kp.gu);
    ECP_random(kp.ge);
    ECP_random(kp.h);
}

void Group_Commit(EC_POINT*&g, EC_POINT*& h,EC_POINT*& u,BIGNUM*& r,Group_Commit_GCM&gcm)
{
    EC_POINT_mul(group, gcm.X, NULL, g, r, bn_ctx);
    EC_POINT_mul(group, gcm.Y, NULL, h, r, bn_ctx);
    EC_POINT_add(group, gcm.Y, gcm.Y, u, bn_ctx);
#ifdef DEBUG
    cout << " Group Commit finishes >>>" << endl;
    //ElGamal_CT_print(CT);
#endif
}