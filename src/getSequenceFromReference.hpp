#if 0
Copyright 2018 David Curtis

This file is part of the geneVarAssoc package.

geneVarAssoc is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

geneVarAssoc is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with geneVarAssoc.If not, see <http://www.gnu.org/licenses/>.
#endif

#ifndef GETSEQUENCEFROMREFERENCEHPP
#define GETSEQUENCEFROMREFERENCEHPP
#include <stdio.h>

class faSequenceFile {
	FILE *fp;
	long startPos;
	int lineLength,gapLength;
public:
	faSequenceFile();
	~faSequenceFile() { fp && fclose(fp); }
	void close() { fp && fclose(fp); fp = 0; }
	int inited() { return fp!=0; }
	int init(const char *fn);
	int getSequence(char *s,int ePos,int len);
	FILE* getFp() { return fp; }
	int getPos(); // return coordinate of current position in file
};

#endif