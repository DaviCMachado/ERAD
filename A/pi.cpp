#include <iostream>
#include <sstream>
#include <cstring>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <iomanip>


using namespace std;

struct bignumber {
	char *digits;
	long unsigned int length;
	bool negative;
	long int point; // 10^point
};

char _n1[1] = { 1 };
bignumber n1 = { _n1, 1, false, 1 };

char _n16[2] = { 1, 6 };
bignumber n16 = { _n16, 2, false, 2 };

#ifdef MODE_DEBUG
	#ifndef DEBUG
		#define DEBUG
	#endif
#endif

#ifdef DEBUG
void print(const char *str, bignumber &number, const long unsigned int d) {
	cout << str << ": [";
	cout << (number.negative ? "-" : "+");

	if(!number.point)
		cout << "0";
	if(number.point >= 1) {
		for (long int i=0; i<number.point; ++i)
			cout << int(number.digits[i]);
		cout << ".";
		for (long int i=number.point; i<=number.length; ++i)
			cout << int(number.digits[i]);
	} else {
		cout << "0.";
		for (long int i=0; i<-number.point; ++i)
			cout << "0";
		for (long int i=0; i<=number.length+number.point; ++i)
			cout << int(number.digits[i]);
	}
	cout << "]:e" << number.point << ":l" << number.length << endl;
}
#endif

void mult_trim(bignumber &result) {
	long int i = 0;
	while(result.digits[i] == 0 && i < result.length) ++i;
	result.point -= i;
	for(long int j=0; i && j<result.length-i; ++i, ++j)
		 result.digits[j] = result.digits[i];
	i = result.length-1;
	while(result.digits[i] == 0 && i >= result.point) --i;
	result.length = i + 1;
}

void mult(bignumber &result, bignumber &number1, bignumber &number2) {
	long int i, j, carry, k, n;
	memset(result.digits, 0, result.length);
	for (i = number1.length - 1; i >= 0; i--) {
		for (j = number2.length - 1, k = i + j + 1, carry = 0; j >= 0; j--, k--) {
			n = number1.digits[i] * number2.digits[j] + result.digits[k] + carry;
			carry = n / 10;
			result.digits[k] = n % 10;
		}
		result.digits[k] += carry;
	}
	result.point = number1.point + number2.point;
	mult_trim(result);
#ifdef DEBUG
	print("mult", result, result.length);
#endif
}

void mult(bignumber &result, bignumber &number1, long int number2) {
	long int i, carry=0, n;
	for (i = number1.length; i >= 0; i--) {
		n = number1.digits[i] * number2 + carry;
		carry = n / 10;
		result.digits[i] = n % 10;
	}
	result.point = number1.point;
#ifdef DEBUG
	print("mult", result, result.length);
#endif
}

void pow16(bignumber &result, long unsigned int k) {
	if(k == 0) {
		result.digits[0] = n1.digits[0];
		result.length = n1.length;
		result.point = n1.point;
	} else {
		result.digits[0] = n16.digits[0];
		result.digits[1] = n16.digits[1];
		result.length = n16.length;
		result.point = n16.point;
		for(long unsigned int i=2; i<=k; ++i) {
			char digits[k * 2 + 10];
			bignumber tmp = { digits, k * 2 + 10, false, 0 };
			mult(tmp, result, n16);
			memcpy(result.digits, tmp.digits, k*2+10);
			result.point = tmp.point;
			result.length = tmp.length;
		}
	}
#ifdef DEBUG
	print("pow16", result, result.length);
#endif
}

void int2bignumber(bignumber &result, long unsigned int num) {
	long int c = 0;
	memset(result.digits, 0, result.length);
	do {
		result.digits[c] = num % 10;
		num = num / 10;
		c++;
	} while(num);
	char tmp;
	for(long int i=0; i<c/2; ++i) {
		tmp = result.digits[i];
		result.digits[i] = result.digits[c-i-1];
		result.digits[c-i-1] = tmp;
	}
	result.length=c;
	result.point=c;
	result.negative=false;
#ifdef DEBUG
	print("int2bignumber", result, result.length);
#endif
}

void add(bignumber &result, bignumber &number1, bignumber &number2, const long unsigned int d) {
	long int carry=0;
	for(long int i=d; i>=0; i--) {
		long int sum = number1.digits[i] + number2.digits[i] + carry;
		result.digits[i] = sum % 10;
		carry = sum / 10;
	}
#ifdef DEBUG
	print("add", result, d);
#endif
}

void sub(bignumber &result, bignumber &minuend, bignumber &subtrahend, const long unsigned int d) {
	long int borrow, diff;
	if(minuend.negative != subtrahend.negative) {
		subtrahend.negative = !subtrahend.negative;
		add(result, minuend, subtrahend, d);
	} else {
		borrow = 0;
		for(long unsigned int i = minuend.length; i>0; --i) {
			diff = minuend.digits[i] - (subtrahend.digits[i] + borrow);
			if(diff < 0) {
				diff = minuend.digits[i] + 10 - (subtrahend.digits[i] + borrow);
				borrow = 1;
			} else {
				borrow = 0;
			}
			result.digits[i] = diff;
		}
		diff = minuend.digits[0] - (subtrahend.digits[0] + borrow);
		if(diff < 0) {
			result.negative = true;
			result.digits[0] = -diff;
		} else {
			result.negative = false;
			result.digits[0] = diff;
		}
#ifdef MORE_DEBUG
		print("sub", result, d);
#endif
	}
}

void one_over_n(bignumber &result, bignumber &dividend, bignumber &divisor, const long unsigned int d) {
	char _one[divisor.length+d+11];
	memset(_one, 0, divisor.length+d+11);
	bignumber one = { _one, divisor.length+1, false, -divisor.length};
	one.digits[0] = dividend.digits[0];

	if(divisor.length == 1 && divisor.digits[0] == 1) {
		memcpy(result.digits, one.digits, one.length);
		result.length = d;
		result.point = 0;
		result.negative = false;
	} else {
		for(long int i=divisor.length;i>0;--i) {
			divisor.digits[i] = divisor.digits[i-1];
		}
		divisor.digits[0]=0;
		divisor.length++;

		for(long int i=0; i<=d; i++) {
			char _part_sub[divisor.length+d+11];
			bignumber part_sub = { _part_sub, divisor.length, false, 0};
			long int c = 0;
			do {
				sub(part_sub, one, divisor, d);
				if(!part_sub.negative) {
					memcpy(one.digits, part_sub.digits, one.length);
					c++;
				} else {
					for(long int i=0;i<one.length;i++)
						one.digits[i] = one.digits[i+1];
					one.digits[one.length-1]=0;
				}
			} while(!part_sub.negative);
			result.digits[i] = c;
		}
		result.length = d+1;
		result.point = 2-divisor.length;
		for(long int i=result.length-result.point+1; i>=-result.point; --i) {
			result.digits[i] = result.digits[i+result.point-1];
		}
		memset(result.digits, 0, -result.point);
		result.point=1;
	}
#ifdef DEBUG
	print("one_over_n", result, d);
#endif
}

void sum(bignumber &output, long unsigned int(*denominator)(const long unsigned int), const long unsigned int d, const long unsigned int n) {
	char digits[d + 11];
	bignumber result = { digits, d+1, false, 1 };
	memset(output.digits, 0, d+11);
	output.point=1;
	output.length=d+1;
	output.negative = false;
	memset(result.digits, 0, d+11);
	for (long unsigned int i = 0; i <= n; ++i) {
		long int den = (*denominator)(i);

#ifdef DEBUG
		cout << i << " - " << den << endl;
#endif

		char _int_den[int(log10(den))+10];
		bignumber big_int_den = { _int_den, int(log10(den))+10, false, 0 };
		int2bignumber(big_int_den, den);

		char _pow[i*2+10];
		bignumber pow = { _pow, i*2+10, false, 0};
		memset(pow.digits, 0, i*2+10);
		pow16(pow, i);

		char _den[d+11];
		bignumber big_den = { _den, d+11, false, 0 };
		memset(big_den.digits, 0, d+11);
		mult(big_den, pow, big_int_den);

		char _term[pow.length+d+11];
		bignumber term = {_term, pow.length+d+11, false, 0};
		memset(term.digits, 0, pow.length+d+11);
		one_over_n(term, n1, big_den, d);

		memcpy(result.digits, output.digits, d+11);
		add(output, result, term, d);
	}
#ifdef DEBUG
	print("sum", output, d);
#endif
}

void convert_to_str(bignumber &output, const long unsigned int d) {
	for(long int i=d+2;i>1;--i)
		output.digits[i] = output.digits[i-1] + '0';
	output.digits[1] = '.';
	output.digits[0] += '0';
	output.digits[d+2] = 0;
}

inline long unsigned int p1(const long unsigned int k) {
  return 8*k+1;
} 

inline long unsigned int p2(const long unsigned int k) {
  return 8*k+4;
} 

inline long unsigned int p3(const long unsigned int k) {
  return 8*k+5;
} 

inline long unsigned int p4(const long unsigned int k) {
  return 8*k+6;
} 

int main() {

	long unsigned int d, n;

	cin >> d >> n;


    char digits[d + 10];
	bignumber output = { digits, d+1, false, 1 };

	char digits_p1[d + 10];
	bignumber output_p1 = { digits_p1 };
	char digits_p2[d + 10];
	bignumber output_p2 = { digits_p2 };
	char digits_p3[d + 10];
	bignumber output_p3 = { digits_p3 };
	char digits_p4[d + 10];
	bignumber output_p4 = { digits_p4 };

	sum(output_p1, &p1, d, n);
	mult(output_p1, output_p1, 4);
	sum(output_p2, &p2, d, n);
	mult(output_p2, output_p2, 2);

	sum(output_p3, &p3, d, n);
	sum(output_p4, &p4, d, n);

	sub(output, output_p1, output_p2, d);
	sub(output, output, output_p3, d);
	sub(output, output, output_p4, d);

	convert_to_str(output, d);

	cout << output.digits << endl;
    end = clock();




    return 0;
}
