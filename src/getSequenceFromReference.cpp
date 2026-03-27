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

#include <ctype.h>
#include <string.h>
#include "getSequenceFromReference.hpp"
#include "dcerror.hpp"
char compBases[5][2]=
{
	{'c','g'},
	{'a','t'},
	{'C','G'},
	{'A','T'},
	{'*','*'}
};

char getCompBase(char b)
{
	int i,j;
	for (i=0;i<5;++i)
		for (j=0;j<2;++j)
			if (b==compBases[i][j])
				return compBases[i][(j+1)%2];
	return '?';
}

int faSequenceFile::init(const char *fn)
{
	int c;
	// idea is to use ftell to get position so it does not matter whether fgetc counts \n as 1 or 2 characters
	long lineEnd,lineStart;
	if (fp!=0)
		fclose (fp);
	fp=fopen(fn,"rb");
	if (fp==0)
	{
		dcerror(1,"Could not open file %s\n",fn);
		return 0;
	}
	while (!isspace(c=fgetc(fp))) ; // skip first line
	do {
		startPos=ftell(fp);
	} while (isspace(c=fgetc(fp)));
	lineStart=startPos;
	do {
		lineEnd=ftell(fp);
	} while (!isspace(c=fgetc(fp)));
	lineLength=lineEnd-lineStart;
	do {
		lineStart=ftell(fp);
	} while (isspace(c=fgetc(fp)));
	gapLength=lineStart-lineEnd;
	// first base has position 1
	return 1;
}

int faSequenceFile::getPos()
// return coordinate of current character
{
	long p;
	long inLine, lineNo;
	int ePos;
	p = ftell(fp);
	p -= startPos;
	lineNo = p / (lineLength + gapLength);
	inLine= p % (lineLength + gapLength);
	ePos = lineNo * lineLength + inLine + 1;
	return ePos;
}

int faSequenceFile::getSequence(char *s,int ePos,int len)
{
	long inLine,lineNo;
	inLine=(ePos-1)%lineLength;
	lineNo=(ePos-1)/lineLength;
	fseek(fp,startPos+lineNo*(lineLength+gapLength)+inLine,SEEK_SET);
	for (;len>0;--len)
	{
		*s++=fgetc(fp);
		if (++inLine==lineLength)
		{
			inLine=0;
			++lineNo;
			fseek(fp,startPos+lineNo*(lineLength+gapLength)+inLine,SEEK_SET);
		}
	}
	*s='\0';
	return 1;
}

faSequenceFile::faSequenceFile()
{
	fp=0;
}

