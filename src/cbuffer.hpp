/*
 * Copyright [2018] [XiaoFei Zhao]
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#ifndef CBUFFER_HPP
#define CBUFFER_HPP

#ifndef IS_INPUT_GZIPPED
#define IS_INPUT_GZIPPED 1
#endif

#if IS_INPUT_GZIPPED
#include "zlib.h"
#define XFILE  gzFile
#define XOPEN  gzopen
#define XCLOSE gzclose
//#define XEOF   gzeof
//#define XGETC  gzgetc
#else
#define XFILE  FILE*
#define XOPEN  fopen
#define XCLOSE fclose
//#define XEOF   feof
//#define XGETC  fgetc
#endif

#include <string>

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFSIZE (128*1024)

const unsigned char RCMAP[256] = {
  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64,
'T', 66,'G', 68, 69, 70,'C',
 72, 73, 74, 75, 76, 77, 78,
 79, 80, 81,     82, 83,'A',
 85, 86, 87,     88, 89, 90,
 91, 92, 93,     94, 95, 96,
't', 98,'g',100,101,102,'c',
104,105,106,107,108,109,110,
111,112,113,    114,115,'a',
117,118,119,    120,121,122,
123,124,125,126,127,
128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,
144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,
160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,
176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,
192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,
208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,
224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,
240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255};

class CBuf {
private:
	char buffer[BUFSIZE+1];
	int bufidx = BUFSIZE;
	int readsize = BUFSIZE;
	XFILE file;
public:
	size_t idx = 0;
	const bool iscasepreserved;
	unsigned char out = '\0';
	const unsigned int size;
	uint64_t slen = 0;
	std::string buf;
	uint64_t chfreqs[256] = {0};
	
	CBuf(std::string fname, int s, bool icp) : iscasepreserved(icp), size(s), buf(std::string(size, 0)) {
		file = XOPEN(fname.c_str(), "r");
		if (NULL == file) {
			fprintf(stderr, "Unable to open the file %s\n", fname.c_str());
			exit(-1);
		}
	}
	~CBuf(void) {
		XCLOSE(file);
	}
	
	unsigned char fgetc_buf();
	bool ceof();
	unsigned char fgetc_visible();
	uint64_t eatnext();
	unsigned char getnewest();
	unsigned char getout();
	std::string tostring();
	unsigned char getith(size_t i);
};

unsigned char CBuf::fgetc_buf() {
	if (bufidx + 1 < readsize) {
		bufidx++;			
		return buffer[bufidx];
	} else {
#if IS_INPUT_GZIPPED
		readsize = gzread(file, buffer, BUFSIZE);
#else
		readsize = fread(buffer, sizeof(char), BUFSIZE, file);
#endif	
		// buffer[readsize] = EOF;
		memset(&buffer[readsize], EOF, BUFSIZE-readsize+1);
		bufidx = 0;
		return buffer[bufidx];
	}
}

bool CBuf::ceof() {
	return bufidx + 1 >= readsize && readsize < BUFSIZE;
}

unsigned char CBuf::fgetc_visible() {
	unsigned char ch = fgetc_buf();
	while (!ceof() && ch < 33) {
		ch = fgetc_buf();
	}
	if (!iscasepreserved) {
		ch = toupper(ch);
	}
	assert(0 <= ch || ceof());
	out = buf[idx];
	buf[idx] = ch;
	idx = (idx + 1) % size;
	return ch;
};

uint64_t CBuf::eatnext() {
	unsigned char ch2 = fgetc_visible();
	if ('>' == ch2 || '@' == ch2) {
		while (!ceof() && (ch2 = fgetc_buf()) != '\n'); // { printf("%c,", ch2); }
		slen = 0;
	} else if ('+' == ch2) {
		for (unsigned int i = 0; i < 3; i++) {
			while (!ceof() && (ch2 = fgetc_buf()) != '\n'); // { printf("(%u,%c), ", i, ch2); }
		}
		slen = 0;
	} else {
		slen++;
		chfreqs[(unsigned int)ch2]++;
	}
	return slen;
}


unsigned char CBuf::getith(size_t i) {
	return buf[(idx + size + i) % size];
};

unsigned char CBuf::getnewest() {
	return buf[(idx + size - 1) % size];
};

unsigned char CBuf::getout() {
	return out;
};

std::string CBuf::tostring() {
	std::string ret;
	for (size_t i = idx; i < idx + size; i++) {
		ret += buf[i%size];
	}
	return ret;
}

#endif

