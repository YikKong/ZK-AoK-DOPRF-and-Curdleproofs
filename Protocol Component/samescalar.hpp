#pragma once
#include "global.hpp"
#include "hash.hpp"
#include "print.hpp"
#include "routines.hpp"
#include "group_commit.hpp"

struct CRS_SameScalar
{
	EC_POINT* G_u;
	EC_POINT* G_e;
	EC_POINT* H;
};

struct CM_SameScalar
{
	EC_POINT* u;
	EC_POINT* e;
	EC_POINT* cm_u;
	EC_POINT* cm_e;
};

struct WT_Samescalar {
	BIGNUM* k;
	BIGNUM* r_u;
	BIGNUM* r_e;
};

void Prove_SameScalar(Group_Commit_KP& kp,
	EC_POINT*& U,EC_POINT*& E,
	BIGNUM*& k,
	BIGNUM*& ru,
	BIGNUM*& re,
	BIGNUM*& zk,
	BIGNUM*& zru,
	BIGNUM*& zre,
	Group_Commit_GCM& cmur,
	Group_Commit_GCM& cmer
	)
{
	BIGNUM* zr_temp = BN_new();
	BN_random(zr_temp);
	BIGNUM* zu_temp = BN_new();
	BN_random(zu_temp);
	BIGNUM* ze_temp = BN_new();
	BN_random(ze_temp);
	BIGNUM* temp = BN_new();
	BIGNUM* c_temp = BN_new();
	string c;

	// cmur;cmer
	EC_POINT* m = EC_POINT_new(group);
	EC_POINT_mul(group, m, NULL, U, zr_temp, bn_ctx);
	Group_Commit(kp.gu, kp.h, m, zu_temp, cmur);

	EC_POINT_mul(group, m, NULL, E, zr_temp, bn_ctx);
	Group_Commit(kp.ge, kp.h, m, ze_temp, cmer);


	// challenge
	c = ECP_ep2string(cmur.X) + ECP_ep2string(cmer.X);
	Hash_String_to_BN(c, c_temp);

	// zk;zru;zre
	BN_mod_mul(temp, c_temp, k, order, bn_ctx);
	BN_mod_add(zk, zr_temp, temp, order, bn_ctx);
	BN_mod_mul(temp, c_temp, ru, order, bn_ctx);
	BN_mod_add(zru, zu_temp, temp, order, bn_ctx);
	BN_mod_mul(temp, c_temp, re, order, bn_ctx);
	BN_mod_add(zre, ze_temp, temp, order, bn_ctx);

	BN_free(zr_temp);
	BN_free(zu_temp);
	BN_free(ze_temp);
	BN_free(temp);
	BN_free(c_temp);
	EC_POINT_free(m);
}

bool Verify_SameScalar(Group_Commit_KP& kp,
	EC_POINT*& U,
	EC_POINT*& E,
	Group_Commit_GCM& cmu,
	Group_Commit_GCM& cme,
	Group_Commit_GCM& cmur,
	Group_Commit_GCM& cmer,
	BIGNUM*& zk,
	BIGNUM*& zru,
	BIGNUM*& zre
)
{
	bool V1, V2, V3, V4;
	string c;
	BIGNUM* c_temp = BN_new();
	// challenge c_temp
	c = ECP_ep2string(cmur.X) + ECP_ep2string(cmer.X);
	Hash_String_to_BN(c, c_temp);

	Group_Commit_GCM cm_temp;
	Group_Commit_GCM_new(cm_temp);
	EC_POINT* m = EC_POINT_new(group);
	EC_POINT* M = EC_POINT_new(group);

	// cm_temp = Groupcommit(zk¡¤U£¬zru)
	EC_POINT_mul(group, m, NULL, U, zk, bn_ctx);
	Group_Commit(kp.gu, kp.h, m, zru, cm_temp);

	// c_temp¡¤cmu.X+cmur.X == cm_temp.X
	EC_POINT_mul(group, M, NULL, cmu.X, c_temp, bn_ctx);
	EC_POINT_add(group, M, M, cmur.X, bn_ctx);
	V1 = (EC_POINT_cmp(group, M, cm_temp.X, bn_ctx) == 0);

	// c_temp¡¤cmu.Y + cmur.Y == cm_temp.Y
	EC_POINT_mul(group, M, NULL, cmu.Y, c_temp, bn_ctx);
	EC_POINT_add(group, M, M, cmur.Y, bn_ctx);
	V2 = (EC_POINT_cmp(group, M, cm_temp.Y, bn_ctx) == 0);

	// cm_temp = Groupcommit(zk¡¤E£¬zre)
	EC_POINT_mul(group, m, NULL, E, zk, bn_ctx);
	Group_Commit(kp.ge, kp.h, m, zre, cm_temp);

	// c_temp¡¤cme.X+cmer.X == cm_temp.X
	EC_POINT_mul(group, M, NULL, cme.X, c_temp, bn_ctx);
	EC_POINT_add(group, M, M, cmer.X, bn_ctx);
	V3 = (EC_POINT_cmp(group, M, cm_temp.X, bn_ctx) == 0);

	// c_temp¡¤cme.Y+cmer.Y == cm_temp.Y
	EC_POINT_mul(group, M, NULL, cme.Y, c_temp, bn_ctx);
	EC_POINT_add(group, M, M, cmer.Y, bn_ctx);
	V4 = (EC_POINT_cmp(group, M, cm_temp.Y, bn_ctx) == 0);
	//V4 = (EC_POINT_cmp(group, M, cm_temp.X, bn_ctx) == 0);

	bool Validity = V1 && V2 && V3 && V4;

	BN_free(c_temp);
	EC_POINT_free(m);
	EC_POINT_free(M);
	Group_Commit_GCM_free(cm_temp);
	return Validity;
}