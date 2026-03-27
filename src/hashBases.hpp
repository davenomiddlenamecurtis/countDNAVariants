#ifndef HASHBASESHPP
#define HASHBASESHPP

class baseHasher {
	unsigned char table['T'+1]; // long enough for all letters
	char bpTable['T'+1];
public:
	baseHasher();
	unsigned int hashBases(const char* in, int len);
	int getComplement(char* out, const char* in, int len);
	int getVariantComplement(char* out, const char* in);
	int convertToSig(char* out, const char* in);
};

extern const char* baseNames,*basePairs[2];
#endif