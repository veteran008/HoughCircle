#include "stdafx.h"
#include "Vis_assist.h"


VisBuf::VisBuf()
{
}

VisBuf::~VisBuf()
{
}

int VisBuf::set_IMG_UBBUF(IMG_UBBUF &srcBuf, IMG_UBYTE *ptr, IMG_SIZE size, IMG_UWORD linestep)
{
	srcBuf.ptr = ptr;
	srcBuf.size = size;
	srcBuf.linestep = linestep;

	return 0;
}

int VisBuf::set_IMG_WBUF(IMG_WBUF & srcBuf, IMG_WORD * ptr, IMG_SIZE size, IMG_UWORD linestep)
{
	srcBuf.ptr = ptr;
	srcBuf.size = size;
	srcBuf.linestep = linestep;
	return 0;
}

int VisBuf::set_IMG_RBUF(IMG_RBUF & srcBuf, IMG_REAL * ptr, IMG_SIZE size, IMG_UWORD linestep)
{
	srcBuf.ptr = ptr;
	srcBuf.size = size;
	srcBuf.linestep = linestep;
	return 0;
}



