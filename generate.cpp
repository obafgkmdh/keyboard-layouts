#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <vector>

#define pack(a, b) (((paired)(a) << 32) | (paired)(b))
#define unpack1(a) ((a) >> 32)
#define unpack2(a) ((a) & 0xffffffff)

typedef unsigned int bitset;
typedef unsigned int count;
typedef unsigned long long paired;

constexpr unsigned int ALPHA_LEN = 30;
constexpr bitset ALPHA_SHIFT = (1 << ALPHA_LEN);
constexpr bitset ALPHA_MASK = ALPHA_SHIFT - 1;

const char alpha[ALPHA_LEN + 1] = "etoainshrlduymwcgf.pb,vk'jx;zq";

bitset tobits(char c1, char c2) {
	bitset bits = 0;
	for (unsigned int i = 0; i < ALPHA_LEN; i++) {
		if (alpha[i] == c1 || alpha[i] == c2) {
			bits |= (1 << i);
		}
	}
	return bits;
}

inline unsigned int ctz(bitset x) {
	return __builtin_ctz(x);
}

inline unsigned int clz(bitset x) {
	return __builtin_clz(x);
}

unsigned int combs[ALPHA_LEN][10] = {0};
void computeCombs() {
	for (int i = 0; i < ALPHA_LEN; i++) {
		combs[i][0] = 1;
		for (int j = 1; j < 10 && j <= i; j++) {
			combs[i][j] = combs[i][j - 1] * (i - j + 1) / j;
		}
	}
}

unsigned int hash9(bitset x) {
	unsigned int h = 4686824;
	unsigned int n = ALPHA_LEN;
	for (int i = 0; i < 9; i++) {
		unsigned int k = ctz(x) + 1;
		n -= k;
		h -= combs[n][9 - i];
		x >>= k;
	}
	return h;
}

inline unsigned int hash2(bitset x) {
	return clz(x) | (ctz(x) << 5);
}

std::vector<count> bigrams(928);
std::vector<count> skipgrams(928);
count totalBigrams, totalSkipgrams;
count cutoffSfb, cutoffSfs;

void loadCorpus(std::string fn) {
	std::ifstream corpus(fn);
	if (!corpus.is_open()) {
		std::cerr << "Could not open file: " << fn << std::endl;
		std::exit(EXIT_FAILURE);
	}
	std::string line;
	while (std::getline(corpus, line)) {
		char prev = 0, prev2 = 0;
		for (char c : line) {
			switch (c) {
				// unshifted keys
				case 'a' ... 'z':
				case ',':
				case '.':
				case ';':
				case '\'':
					break;
				// uppercase letters
				case 'A' ... 'Z':
					c |= 0x20;
					break;
				// shifted symbols
				case '"': c = '\''; break;
				case ':': c = ';'; break;
				default:
					// invalid char
					prev2 = prev;
					prev = 0;
					continue;
			}
			if (prev != 0) {
				totalBigrams++;
				bigrams[hash2(tobits(c, prev))]++;
			}
			if (prev2 != 0) {
				totalSkipgrams++;
				skipgrams[hash2(tobits(c, prev2))]++;
			}
			prev2 = prev;
			prev = c;
		}
	}
	corpus.close();
}

std::vector<std::vector<paired>> groups3(32);
std::vector<count> groups3_sfs(28672);
inline unsigned int hash3(bitset x) {
	return clz(x) | (ctz(x) << 5) | (ctz(x >> (ctz(x) + 1)) << 10);
}

void findGroups3() {
	count sfb2, sfs2, sfb3, sfs3;
	for (bitset i = 1; i < (ALPHA_SHIFT >> 2); i <<= 1) {
		for (bitset j = i << 1; j < (ALPHA_SHIFT >> 1); j <<= 1) {
			if ((sfb2 = bigrams[hash2(i | j)]) > cutoffSfb) { continue; }
			if ((sfs2 = skipgrams[hash2(i | j)]) > cutoffSfs) { continue; }
			for (bitset k = j << 1; k < ALPHA_SHIFT; k <<= 1) {
				if ((sfb3 = sfb2 + bigrams[hash2(i | k)] + bigrams[hash2(j | k)]) > cutoffSfb) { continue; }
				if ((sfs3 = sfs2 + skipgrams[hash2(i | k)] + skipgrams[hash2(j | k)]) > cutoffSfs) { continue; }
				
				bitset ijk = i | j | k;
				groups3[ctz(i)].push_back(pack(sfb3, ijk));
				groups3_sfs[hash3(ijk)] = sfs3;
			}
		}
	}
}

std::vector<std::vector<paired>> groups6;
void groups6Helper(unsigned int n, bitset start, bitset bits, count sfb, count sfs) {
	if (n == 0) {
		groups6[sfb].push_back(pack(sfs, bits));
		return;
	}
	
	for (bitset i = start << 1; i < ALPHA_SHIFT; i <<= 1) {
		count new_sfb = sfb;
		count new_sfs = sfs;
		
		bitset mask = 1;
		while (mask != i) {
			if (bits & mask) {
				new_sfb += bigrams[hash2(mask | i)];
				new_sfs += skipgrams[hash2(mask | i)];
				if (new_sfb > cutoffSfb || new_sfs > cutoffSfs) {
					break;
				}
			}
			mask <<= 1;
		}
		
		if (mask == i) {
			groups6Helper(n-1, i, bits | i, new_sfb, new_sfs);
		}
	}
}

void findGroups6() {
	for (bitset i = 1; i <= (ALPHA_SHIFT >> 6); i <<= 1) {
		groups6Helper(5, i, i, 0, 0);
	}
}

std::vector<std::vector<std::pair<paired, paired>>> groups9(4686825);
count g9count = 0;
void findGroups3N(unsigned int n, bitset start, bitset bits, count sfb, count sfs, paired acc) {
	count sfb2, sfs2, sfb3, sfs3;
	paired p;
	for (bitset i = start << 1; i <= (ALPHA_SHIFT >> n); i <<= 1) {
		if (bits & i) { continue; }
		for (bitset j = i << 1; j <= (ALPHA_SHIFT >> (n-1)); j <<= 1) {
			if (bits & j) { continue; }
			if ((sfb2 = sfb + bigrams[hash2(i | j)]) > cutoffSfb) { continue; }
			if ((sfs2 = sfs + skipgrams[hash2(i | j)]) > cutoffSfs) { continue; }
			for (bitset k = j << 1; k <= (ALPHA_SHIFT >> (n-2)); k <<= 1) {
				if (bits & k) { continue; }
				if ((sfb3 = sfb2 + bigrams[hash2(i | k)] + bigrams[hash2(j | k)]) > cutoffSfb) { continue; }
				if ((sfs3 = sfs2 + skipgrams[hash2(i | k)] + skipgrams[hash2(j | k)]) > cutoffSfs) { continue; }
				
				bitset ijk = i | j | k;
				if (n == 3) {
					unsigned int hijk = hash9(ijk | bits);
					groups9[hijk].push_back(std::make_pair(pack(sfb3, sfs3), acc));
					g9count += 1;
				} else {
					findGroups3N(n-3, i, bits | ijk, sfb3, sfs3, acc | ((paired)(ijk) << (n == 9 ? 32 : 0)));
				}
			}
		}
	}
}


unsigned long long total = 0;
bitset bits_arr[5];
std::ofstream outfile;

void write_bits(bitset p0, bitset p1) {
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < ALPHA_LEN; j++) {
			if (bits_arr[i] & (1 << j)) {
				outfile << alpha[j];
			}
		}
	}
	for (int i = 2; i < 5; i++) {
		for (int j = 0, c = 0; j < ALPHA_LEN; j++) {
			if (bits_arr[i] & (1 << j)) {
				if (c > 0) {
					outfile << alpha[j];
				}
				c++;
			}
		}
	}
	for (int j = 0, c = 0; j < ALPHA_LEN; j++) {
		if (p0 & (1 << j)) {
			if (c > 0) {
				outfile << alpha[j];
			}
			c++;
		}
	}
	for (int j = 0, c = 0; j < ALPHA_LEN; j++) {
		if (p1 & (1 << j)) {
			if (c > 0) {
				outfile << alpha[j];
			}
			c++;
		}
	}
}

void findLayouts(bitset bits, count sfb, count sfs, int n) {
	count new_sfb, new_sfs, sfb9, sfs9;
	bitset b;
	for (const auto& p : groups3[ctz(~bits)]) {
		if (sfb + unpack1(p) > cutoffSfb) { break; }
		if ((p & bits) ||
			(new_sfs = sfs + groups3_sfs[
				hash3(bits_arr[n] = unpack2(p))
			]) > cutoffSfs) { continue; }
		b = bits | bits_arr[n];
		if (n == 4) {
			for (const auto& result : groups9[hash9(b ^ ALPHA_MASK)]) {
				sfb9 = unpack1(result.first);
				sfs9 = unpack2(result.first);
				if (sfb + unpack1(p) + sfb9 <= cutoffSfb && new_sfs + sfs9 <= cutoffSfs) {
					write_bits(unpack1(result.second), unpack2(result.second));
					total += 1;
				}
			}
		} else {
			findLayouts(b, sfb + unpack1(p), new_sfs, n+1);
		}
	}
}

int main(int argc, char* argv[]) {
	if (argc < 2) {
		std::cerr << "usage: " << argv[0] << " corpus_file out_file" << std::endl;
		return 1;
	}
	loadCorpus(argv[1]);
	std::cout << "loaded " << totalBigrams << " bigrams, " << totalSkipgrams << " skipgrams" << std::endl;
	
	cutoffSfb = totalBigrams * 1 / 100;
	cutoffSfs = totalSkipgrams * 7 / 100;
	
	findGroups3();
	for (auto& kv : groups3) { std::sort(kv.begin(), kv.end()); }
	std::cout << "found 3-groups" << std::endl;
	
	groups6.resize(cutoffSfb + 1);
	findGroups6();
	for (auto& v : groups6) { std::sort(v.begin(), v.end()); }
	std::cout << "found 6-groups" << std::endl;
	
	computeCombs();
	findGroups3N(9, 1 << 2, 0, 0, 0, 0);
	std::cout << "found 3-3-3-groups " << g9count << std::endl;
	
	unsigned long long found = 0, done = 0;
	count sfs1, sfs2;
	
	std::cout.precision(4);
	std::cout << std::endl << std::fixed;
	
	outfile.open(argv[2], std::ios::out | std::ios::trunc);
	if (!outfile.is_open()) {
		std::cerr << "Could not open out file" << std::endl;
		std::exit(EXIT_FAILURE);
	}
	
	unsigned long long n_iterations = 0;
	for (count targetSfb = 0; targetSfb <= cutoffSfb; targetSfb++) {
		for (count sfb1 = 0; sfb1 <= targetSfb / 2; sfb1++) {
			for (const auto& entry : groups6[sfb1]) {
				sfs1 = unpack1(entry);
				bits_arr[0] = unpack2(entry);
				for (const auto& entry2 : groups6[targetSfb - sfb1]) {
					sfs2 = sfs1 + unpack1(entry2);
					if (sfs2 > cutoffSfs) { break; }
					bits_arr[1] = unpack2(entry2);
					if (bits_arr[0] & bits_arr[1]) { continue; }
					n_iterations++;
				}
			}
		}
	}
	
	std::cout << "total work: " << n_iterations << std::endl;
	
	for (count targetSfb = 0; targetSfb <= cutoffSfb; targetSfb++) {
		for (count sfb1 = 0; sfb1 <= targetSfb / 2; sfb1++) {
			for (const auto& entry : groups6[sfb1]) {
				sfs1 = unpack1(entry);
				bits_arr[0] = unpack2(entry);
				for (const auto& entry2 : groups6[targetSfb - sfb1]) {
					sfs2 = sfs1 + unpack1(entry2);
					if (sfs2 > cutoffSfs) { break; }
					bits_arr[1] = unpack2(entry2);
					if (bits_arr[0] & bits_arr[1]) { continue; }
					
					findLayouts(bits_arr[0] | bits_arr[1], targetSfb, sfs2, 2);
					done++;
				}
			}
		}
		std::cout << "\r" << (100 * done / (double) n_iterations) << "%, found " << total;
		//if (done > 100000) { break; }
	}
	outfile.close();
	std::cout << "\ntotal: " << total << std::endl;
}