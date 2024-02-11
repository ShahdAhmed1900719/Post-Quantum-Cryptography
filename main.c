#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "toy.h"


int main() {
    short A[TK_K * TK_K * TK_N], t[TK_K * TK_N], s[TK_K * TK_N];
    short u[TK_K * TK_N], v[TK_N] ;
    short mes=4,plain;

    toy_gen(A, t, s);
    toy_enc(A, t, mes, u, v);
    plain = toy_dec(s, u, v);
    printf("message before encryption is %2d\nmessage after decryption is %2d", mes,plain);
  return 0;
    }


//test
/*int main() {
    int data[] = {1, 2, 3, 4};

    printf("Original data: ");
    for (int i = 0; i < TK_N; i++) {
        printf("%d ", data[i]);
    }
    printf("\n");

    reference_ntt(data, 1);  // Forward NTT

    printf("Data after NTT (forward): ");
    for (int i = 0; i < TK_N; i++) {
        printf("%d ", data[i]);
    }
    printf("\n");

    reference_ntt(data, 0);  // Inverse NTT

    printf("Data after NTT (inverse): ");
    for (int i = 0; i < TK_N; i++) {
        printf("%d ", data[i]);
    }
    printf("\n");

    return 0;
}
*/
//test for ntt
/*int main()
{

     short poly1[] = {1, 2, 3, 4};
     short poly2[] = {5, 6, 7, 8};
    // Transform to NTT domain
    ntt(poly1, 1,1);
    ntt(poly2, 1,1);
    // Element-wise multiplication in NTT domain
    for (int i = 0; i < TK_N; i++) {
        poly1[i] = (poly1[i] * poly2[i]) % TK_Q;
    }
    // Inverse NTT to get product polynomial
    ntt(poly1, 0,1);
    // Print the product polynomial coefficients
    printf("Product polynomial coefficients: ");
    for (int i = 0; i < TK_N; i++) {
        printf("%d ", poly1[i]);
    }
    printf("\n");
    return 0;

}*/






