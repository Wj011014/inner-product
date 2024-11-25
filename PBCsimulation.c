#include <pbc/pbc.h>
#include <pbc/pbc_test.h>

#define VECTOR_SIZE 10  // 定义向量大小



// 假设初始化向量的函数
void init_vector(element_t *vec, pairing_t pairing) {
    for (int i = 0; i < VECTOR_SIZE; i++) {
        element_init_Zr(vec[i], pairing);  // 初始化为零元素
        element_random(vec[i]);  // 随机赋值
    }
}

int main(int argc, char **argv) {
    pairing_t pairing;
    element_t rho[VECTOR_SIZE], tk[VECTOR_SIZE], t[VECTOR_SIZE], x[VECTOR_SIZE], y[VECTOR_SIZE];
    element_t alpha, beta, gamma, u, h, g, g1, g2, g3, r, ID, d, d1, d2, d3, uh, uh1, dg3, C1, C2, C3, C4,eggx,egg,C51,mu,H;

    // 初始化配对
    pbc_demo_pairing_init(pairing, argc, argv);
	
	
	double t0 = pbc_get_time();
    // 调用 init_elements 初始化其他元素
    element_init_G1(u, pairing);
	element_init_G1(g1, pairing);
	element_init_G1(g2, pairing);
	element_init_G1(g3, pairing);
    element_init_G1(g, pairing);
    element_init_G1(h, pairing);
    element_init_Zr(alpha, pairing);
	
    element_init_Zr(beta, pairing);
    element_init_Zr(gamma, pairing);
	init_vector(tk, pairing);
	init_vector(rho, pairing);
	element_random(g);
	element_random(alpha);
    element_random(beta);
    element_random(gamma);
	element_set1(g3);
	for (int i = 0; i < VECTOR_SIZE; i++) {
        element_t temp;
        element_init_Zr(temp, pairing);
        element_pow_zn(temp, g, rho[i]);
        element_mul_zn(g3, g3, temp);
    }
	
	for (int i = 0; i < VECTOR_SIZE; i++) {
        element_t temp;
        element_init_Zr(temp, pairing);
        element_pow_zn(temp, g, rho[i]);
        element_mul_zn(tk[i], gamma, rho[i]);  // 修改为 tk[i]
    }
    // 记录 t0 时间
    double t1 = pbc_get_time();
    printf("Setup Phase-----= %fs\n", t1-t0);  // 修正为 t0
	
	// KeyGen Phase
    double t2 = pbc_get_time();
    init_vector(t, pairing);
    init_vector(y, pairing);
	
    element_init_Zr(ID, pairing);
	element_init_Zr(d, pairing);
    element_init_Zr(r, pairing);
    element_init_G1(d1, pairing);
    element_init_G1(d2, pairing);
    element_init_G1(d3, pairing);
    element_init_G1(uh, pairing);
	element_init_G1(uh1, pairing);
    element_init_G1(dg3, pairing);
    element_random(ID);
    element_random(r);
	element_init_GT(eggx, pairing);
	
    element_from_hash(d, "ID", 1);
    
    // 使用 element_neg 来取负值
    element_t neg_r, neg_t;
    element_init_Zr(neg_r, pairing);
    element_init_Zr(neg_t, pairing);
    element_neg(neg_r, r);
    element_neg(neg_t, t[0]);  // 示例，仅取 t[0] 为负数

    element_pow_zn(d3, g, neg_r);  // 使用 neg_r 代替 -r
    element_pow_zn(dg3, g3, neg_t);  // 使用 neg_t 代替 -t[i]
    
    element_pow_zn(uh, u, ID);
    element_mul(uh, uh, h);
    element_pow_zn(uh1, uh, r);

    for (int i = 0; i < VECTOR_SIZE; i++) {
        element_t temp;
        element_init_G1(temp, pairing);
        element_pow_zn(temp, uh1, y[i]);
        element_mul_zn(uh1, uh1, temp);
    }

    element_mul_zn(dg3, dg3, d);
    element_pow_zn(d1, dg3, alpha);
    element_mul_zn(d1, d1, uh1);
    element_pow_zn(d2, dg3, beta);
    element_mul_zn(d2, d2, uh1);

    double t3 = pbc_get_time();
    printf("KeyGen Phase-----= %fs\n", t3 - t2);
	
	// Enc Phase
    double t4 = pbc_get_time();
    element_t s, dg1, dg2, eta, Cx, C5,dg21;
    element_init_Zr(s, pairing);
    element_init_GT(dg1, pairing);
    element_init_GT(dg2, pairing);
    element_init_Zr(eta, pairing);
    element_init_GT(Cx, pairing);
    element_init_GT(C5, pairing);
	element_init_G1(C1, pairing);
	element_init_G1(C2, pairing);
	element_init_GT(C3, pairing);
	element_init_GT(C4, pairing);
	element_init_G1(H, pairing);
	
    element_random(s);
    element_random(eta);
    element_pow_zn(C1, g, s);
    element_pow_zn(C2, uh, s);
    pairing_apply(C3, g1, g3, pairing);
    element_pow_zn(C3, C3, s);
    pairing_apply(C4, g2, g2, pairing);
    element_pow_zn(C4, C4, s);
    pairing_apply(dg1, d, g1, pairing);
    element_pow_zn(dg1, dg1, s);
    pairing_apply(dg2, d, g2, pairing);
    element_pow_zn(dg2, dg2, s);
    element_init_GT(dg21, pairing);
    element_pow_zn(dg21, dg2, eta);
	element_init_GT(egg, pairing);
    pairing_apply(egg, g, g, pairing);
    init_vector(x, pairing);
    element_pow_zn(eggx, egg, x[0]);

    for (int i = 1; i < VECTOR_SIZE; i++) {
        element_t temp;
        element_init_GT(temp, pairing);
        element_pow_zn(temp, egg, x[i]);
        element_mul_zn(eggx, eggx, temp);
    }

    element_mul_zn(Cx, dg1, dg21);
    element_mul_zn(Cx, Cx, eggx);
    element_add(H, C1, C2);
    element_add(H, H, C3);
    element_add(H, H, C4);
    element_add(H, H, Cx);
    element_add(H, H, eta);
	element_init_Zr(mu, pairing);
    element_from_hash(mu, "H", 6);
	element_init_GT(C51, pairing);
    element_pow_zn(C51, dg1, mu);
    element_mul_zn(C5, C51, dg2);
    double t5 = pbc_get_time();
    printf("Enc Phase-----= %fs\n", t5 - t4);

	// Update Phase
    double t6 = pbc_get_time();
    element_t d11, d21, d31, d41, r1, uhr1;
    element_init_Zr(r1, pairing);
    element_init_G1(d11, pairing);
    element_init_G1(d21, pairing);
    element_init_G1(d31, pairing);
    element_init_G1(d41, pairing);
    element_init_G1(uhr1, pairing);
    element_random(r1);

    element_pow_zn(d11, dg1, r1);
    element_pow_zn(d21, dg2, r1);
    element_pow_zn(d31, dg3, r1);
    element_pow_zn(uhr1, u, r1);
    element_mul_zn(uhr1, uhr1, h);
    element_pow_zn(d41, uhr1, r1);

    for (int i = 0; i < VECTOR_SIZE; i++) {
        element_t temp;
        element_init_G1(temp, pairing);
        element_pow_zn(temp, uhr1, t[i]);
        element_mul_zn(d41, d41, temp);
    }

    element_mul_zn(d41, d41, dg3);
    double t7 = pbc_get_time();
    printf("Update Phase-----= %fs\n", t7 - t6);
    //Dec
    double t8 = pbc_get_time();
    element_t Cxy,w1,w2,sumdg1,sumdg2,C1d1,C2d3,C1d2,sumC3,sumC4,sum,xy,ty;
    element_init_GT(Cxy, pairing);
	element_init_Zr(xy, pairing);
	element_init_Zr(ty, pairing);
	element_set0(xy);
	element_set0(ty);
    element_init_GT(w1, pairing);
    element_init_GT(C1d1, pairing);
    element_init_GT(sum, pairing);
    element_init_GT(C1d2, pairing);
    element_init_GT(sumC3, pairing);
    element_init_GT(sumC4, pairing);
    element_init_GT(w2, pairing);
    element_init_GT(C2d3, pairing);
    element_init_GT(sumdg1, pairing);
    element_init_GT(sumdg2, pairing);
	for (int i = 0; i < VECTOR_SIZE; i++) {
        element_t temp;
		element_init_Zr(temp, pairing);
		element_mul(temp, x[i], y[i]);
		element_add(xy, xy, temp);  // 累加结果
        element_clear(temp);  // 清理临时变量
   }

    element_printf("Dot product result: %B\n", xy);
    
    element_pow_zn(sumdg1, dg1, y[0]);
    element_pow_zn(sumdg2, dg21, y[0]);
    for (int i = 1; i < VECTOR_SIZE; i++) {
        element_t temp,temp1;
        element_init_GT(temp, pairing);
        element_init_GT(temp1, pairing);
        element_pow_zn(temp, dg1, y[i]);
        element_pow_zn(temp1, dg21, y[i]);
        element_mul_zn(sumdg1, sumdg1,temp);
        element_mul_zn(sumdg2, sumdg2,temp1);
    }

    element_pow_zn(Cxy, egg, xy);
    element_mul_zn(Cxy, Cxy,sumdg1);
    element_mul_zn(Cxy, Cxy,sumdg2);
    pairing_apply(C1d1, C1, d1, pairing);
    pairing_apply(C1d2, C1, d2, pairing);
    pairing_apply(C2d3, C2, d3, pairing);
	for (int i = 0; i < VECTOR_SIZE; i++) {
        element_t temp;
		element_init_Zr(temp, pairing);
		element_mul(temp, t[i], y[i]);
		element_add(ty, ty, temp);  // 累加结果
        element_clear(temp);  // 清理临时变量
   }
element_printf("Dot product result: %B\n", ty);
	
    
   
    element_mul_zn(w1, C1d1,C2d3);
    element_pow_zn(sumC3, C3, ty);
    element_pow_zn(sumC4, C4, ty);
    element_mul_zn(w1, w1,sumC3);
    element_mul_zn(w2, C1d2,C2d3);
    element_mul_zn(w2, w2,sumC4);
    element_pow_zn(w2, w2,eta);
    element_mul_zn(sum, w2,w1);
    element_pow_zn(egg, egg,xy);
    element_mul_zn(sum, egg,sum);
    if(Cxy==sum)
        printf("OK\n");
    double t9 = pbc_get_time();
    printf("Dec Phase-----= %fs\n", t9 - t8);
    
    pairing_clear(pairing);
}
