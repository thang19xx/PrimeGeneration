/*thangdn - 17/5/2015
Big prime number - Miller-Rabin tests*/
#include<iostream>
#include <bitset>
#include <gmpxx.h>
#include <time.h>
#include <math.h>
#include <fstream>
using namespace std;
#define N 3072
#define LENGTH (int)(0.693*N*2)
#define _TIME int starts,finishs;
#define STARTS_TIME starts=clock();
#define FINISHS_TIME finishs=clock(); cout<<"Time ~ "<<(double)(finishs-starts)/CLOCKS_PER_SEC<<endl;


typedef mpz_class ZZ;
ZZ prime[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197,
199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331,337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461,
463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751,
757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049,
1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307,
1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601,
1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889,
1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203,
2207, 2213, 2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 2441, 2447, 2459, 2467, 2473,
2477, 2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741 };
ZZ modulo(const ZZ &a,const ZZ&b,const ZZ&c)//a^b%c
{
	ZZ x; x = 1;
	ZZ y = a;
	ZZ k = b,m;
	while (k > 0){
		m = k & 1;
		if (m == 1)
			x = (x*y) % c;
		y = (y*y) % c;
		k >>=1;
	}
	return x%c;
}

bool WillerRabin(const ZZ& n, const ZZ& x, long &k, const ZZ&m)
{
	ZZ z,y;
	//z = modulo(x, m, n);
	mpz_powm(z.get_mpz_t(),x.get_mpz_t(),m.get_mpz_t(),n.get_mpz_t());
	if (z == 1 || z == n - 1) return 0;

	long j=0;
	do {
		y = z;
		//z=y*y%n;
		mpz_powm_ui(z.get_mpz_t(),y.get_mpz_t(),2,n.get_mpz_t());
		j++;
	} while (j < k && z != 1);
	return (z != 1 || y != n - 1)?true:false;
	//return true;
}
bitset<LENGTH> BC(const ZZ& n)
{
	bitset<LENGTH> S;
	int j = 0;
	while (j<400) {
		ZZ p = prime[j++];
		ZZ r = n%p;
		for (ZZ i = p - r; i < LENGTH; i += p)
		{
			S[i.get_ui()] = true;							//S |= (1 << i);
		}
	}
	return S;
}

bool PrimeTest(const ZZ& n)
{
    ZZ y, m = n^1;                    //sub(m, n, 1); because n is odd
	long k;					// MakeOdd(m);
	k = 0;
	do{
		m >>= 1;				//= > m = m / 2;
		y = m & 1;				//check m is even number
		k++;
	} while (y == 0);

    gmp_randclass r(gmp_randinit_default);
    ZZ x;
	for (int i = 0; i < 20; i++) {
		x = prime[i];
		if (WillerRabin(n, x, k, m))return false;
	}
	for (int i = 0; i < 44; i++)
	{
		do {
        x= r.get_z_range(n); //RandomBnd(x, n);
		} while (x <= 73);
		if (WillerRabin(n, x, k, m))return false;
	}

	return true;
}


int main()
{
	ZZ number;
    gmp_randclass rr(gmp_randinit_default);

    _TIME
	while (true)
	{
		STARTS_TIME
		rr.seed(time(NULL));
        number =rr.get_z_bits(N);
		ZZ num = number & 1;
		if (num == 0)number = number | 1;
		bitset<LENGTH> bitGet = BC(number);

		int i = 0;
		while (i < LENGTH)
		{
			if (bitGet[i] == 0)
			{
			    ZZ k=number+i;
				if (PrimeTest(k)){
					cout<<k.get_str(10)<<endl;
					FINISHS_TIME
					STARTS_TIME
				}
			}
			i+=2;
		}
	}

}
